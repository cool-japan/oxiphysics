// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Operational-space control: task-space inertia matrix Λ = (J · M⁻¹ · Jᵀ)⁻¹.
//!
//! Reference: Khatib 1987 "A unified approach for motion and force control of robot
//! manipulators: The operational space formulation."

use crate::{crba::compute_mass_matrix_crba, model::ArticulatedModel};

/// Compute the operational-space (task-space) inertia matrix.
///
/// Given the joint-space mass matrix `M(q)` (computed via CRBA) and a task-space
/// Jacobian `J ∈ ℝ^{6×n}`, returns Λ = (J · M⁻¹ · Jᵀ)⁻¹ ∈ ℝ^{6×6}.
///
/// `jacobian` is a `6×n` matrix represented as a slice of 6 rows, each row being
/// a `Vec<f64>` of length `model.total_dof()`.
///
/// # Panics
///
/// Panics if `jacobian` does not have exactly 6 rows, or if any row does not have
/// length `model.total_dof()`.
pub fn compute_lambda(
    model: &mut ArticulatedModel,
    q: &[f64],
    jacobian: &[Vec<f64>],
) -> [[f64; 6]; 6] {
    let n = model.total_dof();
    assert_eq!(jacobian.len(), 6, "Jacobian must have 6 rows");
    for row in jacobian {
        assert_eq!(
            row.len(),
            n,
            "Each Jacobian row must have length total_dof()={n}"
        );
    }

    // Step 1: mass matrix M (n×n)
    let m = compute_mass_matrix_crba(model, q);

    // Step 2: Cholesky factorize M = L · L^T
    let l = cholesky_lower(&m).expect("mass matrix must be positive-definite");

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

    // Step 5: invert A (6×6) via Gauss-Jordan with partial pivoting
    invert_6x6(a).expect("operational-space inertia matrix must be invertible")
}

/// Dense lower-triangular Cholesky factorization: A = L · L^T.
///
/// Returns `None` if `A` is not positive-definite (diagonal entry ≤ 0 after
/// subtraction). `A` must be symmetric positive-definite.
#[allow(clippy::needless_range_loop)]
fn cholesky_lower(a: &[Vec<f64>]) -> Option<Vec<Vec<f64>>> {
    let n = a.len();
    let mut l = vec![vec![0.0f64; n]; n];
    // Double-indexed loops (i, j, k) access two or three different Vecs simultaneously;
    // iterator form would require unsafe split_at_mut or significantly more complex code.
    for i in 0..n {
        for j in 0..=i {
            let mut sum = a[i][j];
            for k in 0..j {
                sum -= l[i][k] * l[j][k];
            }
            if i == j {
                if sum <= 0.0 {
                    return None;
                }
                l[i][j] = sum.sqrt();
            } else {
                l[i][j] = sum / l[j][j];
            }
        }
    }
    Some(l)
}

/// Solve L · L^T · x = b using forward/backward substitution.
#[allow(clippy::needless_range_loop)]
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
#[allow(clippy::needless_range_loop)]
fn invert_6x6(a: [[f64; 6]; 6]) -> Option<[[f64; 6]; 6]> {
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
        for row in col + 1..6 {
            if aug[row][col].abs() > max_val {
                max_val = aug[row][col].abs();
                max_row = row;
            }
        }
        if max_val < 1e-15 {
            return None; // singular
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
#[allow(clippy::needless_range_loop)]
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
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (l[i][j] - expected).abs() < 1e-14,
                    "L[{i}][{j}] = {}",
                    l[i][j]
                );
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
        for i in 0..6 {
            a[i][i] = 1.0;
        }
        let inv = invert_6x6(a).expect("identity should invert");
        for i in 0..6 {
            for j in 0..6 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (inv[i][j] - expected).abs() < 1e-14,
                    "inv[{i}][{j}]={}",
                    inv[i][j]
                );
            }
        }
    }
}
