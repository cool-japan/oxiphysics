// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Operational-space control: task-space inertia matrix Λ = (J · M⁻¹ · Jᵀ)⁻¹.
//!
//! Reference: Khatib 1987 "A unified approach for motion and force control of robot
//! manipulators: The operational space formulation."

use crate::{crba::compute_mass_matrix_crba, error::ArticulatedError, model::ArticulatedModel};

/// Relative pivot tolerance used to declare `J·M⁻¹·Jᵀ` singular.
const SINGULAR_REL_TOL: f64 = 1e-12;

/// Compute the operational-space (task-space) inertia matrix.
///
/// Given the joint-space mass matrix `M(q)` (computed via CRBA) and a task-space
/// Jacobian `J ∈ ℝ^{6×n}`, returns Λ = (J · M⁻¹ · Jᵀ)⁻¹ ∈ ℝ^{6×6}.
///
/// `jacobian` is a `6×n` matrix represented as a slice of 6 rows, each row being
/// a `Vec<f64>` of length `model.total_dof()`. `q` has length
/// `model.total_nq()`.
///
/// # Errors
///
/// - [`ArticulatedError::DimensionMismatch`] if `q`, the number of Jacobian
///   rows (must be 6), or any row length (must be `model.total_dof()`) is wrong.
/// - [`ArticulatedError::NotPositiveDefinite`] if `M(q)` is not SPD.
/// - [`ArticulatedError::SingularTaskInertia`] if `J·M⁻¹·Jᵀ` is singular —
///   always the case for robots with fewer than 6 DoF or at a kinematic
///   singularity. Use [`compute_lambda_damped`] to obtain a regularised Λ.
pub fn compute_lambda(
    model: &ArticulatedModel,
    q: &[f64],
    jacobian: &[Vec<f64>],
) -> Result<[[f64; 6]; 6], ArticulatedError> {
    let a = task_inertia_inverse(model, q, jacobian)?;
    invert_6x6(a).ok_or(ArticulatedError::SingularTaskInertia)
}

/// Damped operational-space inertia `Λ_λ = (J · M⁻¹ · Jᵀ + λ²·I)⁻¹`.
///
/// The Tikhonov term `λ²·I` keeps the inverse well defined for rank-deficient
/// tasks (fewer than 6 DoF, kinematic singularities); `λ = 0` reduces to
/// [`compute_lambda`]. Arguments are as for [`compute_lambda`].
///
/// # Errors
///
/// - [`ArticulatedError::InvalidParameter`] if `damping` is negative or not finite.
/// - Otherwise the same errors as [`compute_lambda`]
///   ([`ArticulatedError::SingularTaskInertia`] can only occur for `λ = 0` or
///   a `λ` negligible against the scale of `J·M⁻¹·Jᵀ`).
pub fn compute_lambda_damped(
    model: &ArticulatedModel,
    q: &[f64],
    jacobian: &[Vec<f64>],
    damping: f64,
) -> Result<[[f64; 6]; 6], ArticulatedError> {
    if !damping.is_finite() || damping < 0.0 {
        return Err(ArticulatedError::InvalidParameter {
            what: "damping",
            value: damping,
        });
    }
    let mut a = task_inertia_inverse(model, q, jacobian)?;
    let lambda2 = damping * damping;
    for (i, a_row) in a.iter_mut().enumerate() {
        a_row[i] += lambda2;
    }
    invert_6x6(a).ok_or(ArticulatedError::SingularTaskInertia)
}

/// Validate inputs and form `A = J · M⁻¹ · Jᵀ` (6×6).
fn task_inertia_inverse(
    model: &ArticulatedModel,
    q: &[f64],
    jacobian: &[Vec<f64>],
) -> Result<[[f64; 6]; 6], ArticulatedError> {
    let n = model.total_dof();
    if q.len() != model.total_nq() {
        return Err(ArticulatedError::DimensionMismatch {
            what: "q",
            expected: model.total_nq(),
            got: q.len(),
        });
    }
    if jacobian.len() != 6 {
        return Err(ArticulatedError::DimensionMismatch {
            what: "jacobian rows",
            expected: 6,
            got: jacobian.len(),
        });
    }
    for row in jacobian {
        if row.len() != n {
            return Err(ArticulatedError::DimensionMismatch {
                what: "jacobian row length",
                expected: n,
                got: row.len(),
            });
        }
    }

    // Step 1: mass matrix M (n×n)
    let m = compute_mass_matrix_crba(model, q);

    // Step 2: Cholesky factorize M = L · L^T
    let l = cholesky_lower(&m)?;

    // Step 3: M^{-1} · J^T — solve L · L^T · x = j_col for each column of J^T
    // j_col is the k-th column of J^T = the k-th row of J
    let mut m_inv_jt = vec![vec![0.0f64; n]; 6]; // m_inv_jt[k][i] = (M^-1 J^T)[i,k]
    for (k, jk) in jacobian.iter().enumerate() {
        let x = cholesky_solve(&l, jk);
        m_inv_jt[k].copy_from_slice(&x);
    }

    // Step 4: J · (M^{-1} · J^T) → 6×6 matrix A
    // A[r][c] = sum_i J[r][i] * (M^-1 J^T)[c][i] = J[r] · m_inv_jt[c]
    let mut a = [[0.0f64; 6]; 6];
    for r in 0..6 {
        for c in 0..6 {
            a[r][c] = jacobian[r]
                .iter()
                .zip(m_inv_jt[c].iter())
                .map(|(j_ri, m_ci)| j_ri * m_ci)
                .sum();
        }
    }
    Ok(a)
}

/// Dense lower-triangular Cholesky factorization: A = L · L^T.
///
/// Returns [`ArticulatedError::NotPositiveDefinite`] if `A` is not
/// positive-definite (diagonal entry ≤ 0 after subtraction).
fn cholesky_lower(a: &[Vec<f64>]) -> Result<Vec<Vec<f64>>, ArticulatedError> {
    let n = a.len();
    let mut l = vec![vec![0.0f64; n]; n];
    // Double-indexed loops (i, j, k) access two or three different Vecs simultaneously;
    // iterator form would require unsafe split_at_mut or significantly more complex code.
    for i in 0..n {
        for j in 0..=i {
            let mut sum = a[i][j];
            for (l_ik, l_jk) in l[i][..j].iter().zip(l[j][..j].iter()) {
                sum -= l_ik * l_jk;
            }
            if i == j {
                if sum <= 0.0 || !sum.is_finite() {
                    return Err(ArticulatedError::NotPositiveDefinite { index: i });
                }
                l[i][j] = sum.sqrt();
            } else {
                l[i][j] = sum / l[j][j];
            }
        }
    }
    Ok(l)
}

/// Solve L · L^T · x = b using forward/backward substitution.
fn cholesky_solve(l: &[Vec<f64>], b: &[f64]) -> Vec<f64> {
    let n = l.len();
    // Forward substitution: L · y = b
    // Index i used to access l[i], l[i][j], y[j] simultaneously.
    let mut y = vec![0.0f64; n];
    for i in 0..n {
        let mut s = b[i];
        for j in 0..i {
            s -= l[i][j] * y[j];
        }
        y[i] = s / l[i][i];
    }
    // Backward substitution: L^T · x = y
    // Index i reversed, accesses l[j][i] (column i of L^T).
    let mut x = vec![0.0f64; n];
    for i in (0..n).rev() {
        let mut s = y[i];
        for j in i + 1..n {
            s -= l[j][i] * x[j];
        }
        x[i] = s / l[i][i];
    }
    x
}

/// Invert a 6×6 matrix via Gauss-Jordan elimination with partial pivoting.
///
/// Returns `None` when a pivot falls below `SINGULAR_REL_TOL · max|a_ij|`
/// (numerically singular relative to the matrix scale), or for a zero /
/// non-finite matrix.
fn invert_6x6(a: [[f64; 6]; 6]) -> Option<[[f64; 6]; 6]> {
    let scale = a
        .iter()
        .flat_map(|row| row.iter())
        .fold(0.0f64, |acc, v| acc.max(v.abs()));
    if !scale.is_finite() || scale == 0.0 {
        return None;
    }
    let pivot_tol = SINGULAR_REL_TOL * scale;
    let mut aug = [[0.0f64; 12]; 6];
    // Build augmented matrix [A | I]
    // Inner loops access fixed-size arrays by two indices simultaneously.
    for i in 0..6 {
        for j in 0..6 {
            aug[i][j] = a[i][j];
        }
        aug[i][6 + i] = 1.0;
    }

    for col in 0..6 {
        // Partial pivot: find row with largest absolute value in this column
        let mut max_row = col;
        let mut max_val = aug[col][col].abs();
        for (row, aug_row) in aug[col + 1..]
            .iter()
            .enumerate()
            .map(|(k, r)| (col + 1 + k, r))
        {
            if aug_row[col].abs() > max_val {
                max_val = aug_row[col].abs();
                max_row = row;
            }
        }
        if max_val <= pivot_tol {
            return None; // singular relative to the matrix scale
        }
        aug.swap(col, max_row);

        // Normalise pivot row
        let pivot = aug[col][col];
        for entry in aug[col].iter_mut() {
            *entry /= pivot;
        }
        // Eliminate column from all other rows
        for row in 0..6 {
            if row != col {
                let factor = aug[row][col];
                let pivot_row = aug[col]; // copy (fixed-size array)
                for j in 0..12 {
                    aug[row][j] -= factor * pivot_row[j];
                }
            }
        }
    }

    let mut result = [[0.0f64; 6]; 6];
    for i in 0..6 {
        for j in 0..6 {
            result[i][j] = aug[i][6 + j];
        }
    }
    Some(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cholesky_identity() {
        let a = vec![
            vec![1.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
            vec![0.0, 0.0, 1.0],
        ];
        let l = cholesky_lower(&a).expect("identity should factor");
        // L should be identity
        for (i, l_row) in l.iter().enumerate() {
            for (j, &l_ij) in l_row.iter().enumerate() {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((l_ij - expected).abs() < 1e-14, "L[{i}][{j}] = {}", l_ij);
            }
        }
    }

    #[test]
    fn test_cholesky_solve_simple() {
        // A = [[4,2],[2,3]], L = [[2,0],[1,√2]]
        let a = vec![vec![4.0, 2.0], vec![2.0, 3.0]];
        let l = cholesky_lower(&a).expect("should factor");
        let b = vec![1.0, 0.0];
        let x = cholesky_solve(&l, &b);
        // Verify A · x ≈ b
        let ax0 = 4.0 * x[0] + 2.0 * x[1];
        let ax1 = 2.0 * x[0] + 3.0 * x[1];
        assert!((ax0 - b[0]).abs() < 1e-12, "A·x[0]={ax0}");
        assert!((ax1 - b[1]).abs() < 1e-12, "A·x[1]={ax1}");
    }

    #[test]
    fn test_invert_6x6_identity() {
        let mut a = [[0.0f64; 6]; 6];
        for (i, a_row) in a.iter_mut().enumerate() {
            a_row[i] = 1.0;
        }
        let inv = invert_6x6(a).expect("identity should invert");
        for (i, inv_row) in inv.iter().enumerate() {
            for (j, &inv_ij) in inv_row.iter().enumerate() {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (inv_ij - expected).abs() < 1e-14,
                    "inv[{i}][{j}]={}",
                    inv_ij
                );
            }
        }
    }
}
