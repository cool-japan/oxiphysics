// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Sparse symmetric-positive-definite linear algebra used by the Biot
//! [`ConsolidationSolver`](super::ConsolidationSolver): CSR validation and
//! Dirichlet masking, IC(0) / Jacobi preconditioners and preconditioned CG.

use crate::error::{Error, Result};
use crate::sparse::CsrMatrix;

/// Extract the entries of a CSR matrix as triplets after validating its
/// structure and values.
pub(super) fn csr_triplets_checked(matrix: &CsrMatrix) -> Result<Vec<(usize, usize, f64)>> {
    let n_rows = matrix.nrows;
    if matrix.row_ptr.len() != n_rows + 1 {
        return Err(Error::InvalidInput(format!(
            "row_ptr has length {}, expected {}",
            matrix.row_ptr.len(),
            n_rows + 1
        )));
    }
    let nnz = matrix.row_ptr.get(n_rows).copied().unwrap_or(0);
    if matrix.row_ptr.first().copied() != Some(0)
        || matrix.col_indices.len() != nnz
        || matrix.values.len() != nnz
        || matrix.row_ptr.windows(2).any(|w| w[0] > w[1])
    {
        return Err(Error::InvalidInput(
            "malformed CSR arrays (row_ptr / col_indices / values)".to_string(),
        ));
    }
    let mut triplets = Vec::with_capacity(nnz);
    for row in 0..n_rows {
        for idx in matrix.row_ptr[row]..matrix.row_ptr[row + 1] {
            let col = matrix.col_indices[idx];
            let value = matrix.values[idx];
            if col >= matrix.ncols {
                return Err(Error::DofOutOfRange {
                    dof: col,
                    n: matrix.ncols,
                });
            }
            if !value.is_finite() {
                return Err(Error::InvalidInput(format!(
                    "non-finite value {value} at ({row}, {col})"
                )));
            }
            triplets.push((row, col, value));
        }
    }
    Ok(triplets)
}

/// Copy of `matrix` restricted to free rows and columns; optionally with a unit
/// diagonal on constrained rows (so that the operator stays non-singular).
pub(super) fn mask_matrix(
    matrix: &CsrMatrix,
    row_free: &[bool],
    col_free: &[bool],
    unit_constrained_diagonal: bool,
) -> CsrMatrix {
    let mut triplets = Vec::with_capacity(matrix.nnz() + row_free.len());
    for (row, &free_row) in row_free.iter().enumerate().take(matrix.nrows) {
        if free_row {
            for idx in matrix.row_ptr[row]..matrix.row_ptr[row + 1] {
                let col = matrix.col_indices[idx];
                if col_free.get(col).copied().unwrap_or(false) {
                    triplets.push((row, col, matrix.values[idx]));
                }
            }
        } else if unit_constrained_diagonal && row < matrix.ncols {
            triplets.push((row, row, 1.0));
        }
    }
    CsrMatrix::from_triplets(matrix.nrows, matrix.ncols, &triplets)
}

/// Symmetric positive definite preconditioner for the stiffness block.
#[derive(Debug, Clone)]
pub(super) enum SpdPreconditioner {
    /// Incomplete Cholesky factor `L` (IC(0), possibly diagonally shifted).
    IncompleteCholesky(IncompleteCholeskyFactor),
    /// Inverse diagonal (fallback when IC(0) breaks down for every shift).
    Jacobi(Vec<f64>),
}

impl SpdPreconditioner {
    /// IC(0) of `a`, retried with growing diagonal shifts on breakdown
    /// (Manteuffel); falls back to Jacobi.  `diag` must be positive.
    pub(super) fn new(a: &CsrMatrix, diag: &[f64]) -> Self {
        let mut shift = 0.0_f64;
        for _ in 0..12 {
            if let Some(factor) = IncompleteCholeskyFactor::new(a, shift) {
                return Self::IncompleteCholesky(factor);
            }
            shift = if shift == 0.0 { 1e-3 } else { shift * 4.0 };
        }
        Self::Jacobi(diag.iter().map(|d| 1.0 / d).collect())
    }

    pub(super) fn apply(&self, r: &[f64], z: &mut [f64]) {
        match self {
            Self::IncompleteCholesky(factor) => factor.solve(r, z),
            Self::Jacobi(inv_diag) => {
                for ((zi, ri), di) in z.iter_mut().zip(r).zip(inv_diag) {
                    *zi = ri * di;
                }
            }
        }
    }
}

/// Lower-triangular IC(0) factor in CSR form; the diagonal is the last entry
/// of every row.
#[derive(Debug, Clone)]
pub(super) struct IncompleteCholeskyFactor {
    row_ptr: Vec<usize>,
    col_indices: Vec<usize>,
    values: Vec<f64>,
}

impl IncompleteCholeskyFactor {
    /// IC(0) of `A + shift·diag(A)` on the lower-triangular pattern of `A`
    /// (rows must have sorted column indices).  Returns `None` on breakdown.
    fn new(a: &CsrMatrix, shift: f64) -> Option<Self> {
        let n = a.nrows;
        let mut row_ptr = Vec::with_capacity(n + 1);
        row_ptr.push(0);
        let mut col_indices = Vec::new();
        let mut values = Vec::new();
        for i in 0..n {
            for idx in a.row_ptr[i]..a.row_ptr[i + 1] {
                let c = a.col_indices[idx];
                if c <= i {
                    col_indices.push(c);
                    values.push(if c == i {
                        a.values[idx] * (1.0 + shift)
                    } else {
                        a.values[idx]
                    });
                }
            }
            if col_indices.last().copied() != Some(i) {
                return None;
            }
            row_ptr.push(col_indices.len());
        }
        for i in 0..n {
            let start = row_ptr[i];
            let diag_pos = row_ptr[i + 1] - 1;
            let a_ii = values[diag_pos];
            for p in start..diag_pos {
                let k = col_indices[p];
                let k_start = row_ptr[k];
                let k_diag = row_ptr[k + 1] - 1;
                // dot(L_i[0..k), L_k[0..k)) over the common pattern.
                let mut dot = 0.0;
                let mut pi = start;
                let mut pk = k_start;
                while pi < p && pk < k_diag {
                    let ci = col_indices[pi];
                    let ck = col_indices[pk];
                    if ci == ck {
                        dot += values[pi] * values[pk];
                        pi += 1;
                        pk += 1;
                    } else if ci < ck {
                        pi += 1;
                    } else {
                        pk += 1;
                    }
                }
                values[p] = (values[p] - dot) / values[k_diag];
            }
            let sum_sq: f64 = values[start..diag_pos].iter().map(|v| v * v).sum();
            let pivot = a_ii - sum_sq;
            if !pivot.is_finite() || pivot <= 1e-12 * a_ii.abs() {
                return None;
            }
            values[diag_pos] = pivot.sqrt();
        }
        Some(Self {
            row_ptr,
            col_indices,
            values,
        })
    }

    /// Solve `L Lᵀ z = r`.
    fn solve(&self, r: &[f64], z: &mut [f64]) {
        let n = self.row_ptr.len().saturating_sub(1);
        for i in 0..n {
            let diag_pos = self.row_ptr[i + 1] - 1;
            let mut s = r[i];
            for p in self.row_ptr[i]..diag_pos {
                s -= self.values[p] * z[self.col_indices[p]];
            }
            z[i] = s / self.values[diag_pos];
        }
        for i in (0..n).rev() {
            let diag_pos = self.row_ptr[i + 1] - 1;
            z[i] /= self.values[diag_pos];
            let zi = z[i];
            for p in self.row_ptr[i]..diag_pos {
                z[self.col_indices[p]] -= self.values[p] * zi;
            }
        }
    }
}

/// Outcome of a converged CG solve.
#[derive(Debug, Clone, Copy)]
pub(super) struct CgOutcome {
    pub(super) iterations: usize,
    pub(super) relative_residual: f64,
}

pub(super) fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

/// Preconditioned conjugate gradients for an SPD operator.
///
/// Iterates from the initial guess in `x` until `‖b − A x‖ ≤ tol·‖b‖`.
pub(super) fn preconditioned_cg<A, P>(
    mut apply_a: A,
    mut apply_m: P,
    b: &[f64],
    x: &mut [f64],
    tolerance: f64,
    max_iterations: usize,
    label: &str,
) -> Result<CgOutcome>
where
    A: FnMut(&[f64], &mut [f64]) -> Result<()>,
    P: FnMut(&[f64], &mut [f64]),
{
    let n = b.len();
    let b_norm = dot(b, b).sqrt();
    if b_norm == 0.0 {
        x.iter_mut().for_each(|v| *v = 0.0);
        return Ok(CgOutcome {
            iterations: 0,
            relative_residual: 0.0,
        });
    }
    let mut ax = vec![0.0; n];
    apply_a(x, &mut ax)?;
    let mut r: Vec<f64> = b.iter().zip(&ax).map(|(bi, ai)| bi - ai).collect();
    let mut r_norm = dot(&r, &r).sqrt();
    if r_norm <= tolerance * b_norm {
        return Ok(CgOutcome {
            iterations: 0,
            relative_residual: r_norm / b_norm,
        });
    }
    let mut z = vec![0.0; n];
    apply_m(&r, &mut z);
    let mut direction = z.clone();
    let mut rz = dot(&r, &z);
    let mut a_dir = vec![0.0; n];
    for iteration in 1..=max_iterations {
        apply_a(&direction, &mut a_dir)?;
        let curvature = dot(&direction, &a_dir);
        if !(curvature.is_finite() && curvature > 0.0 && rz.is_finite() && rz > 0.0) {
            return Err(Error::SingularSystem);
        }
        let alpha = rz / curvature;
        for i in 0..n {
            x[i] += alpha * direction[i];
            r[i] -= alpha * a_dir[i];
        }
        r_norm = dot(&r, &r).sqrt();
        if r_norm <= tolerance * b_norm {
            return Ok(CgOutcome {
                iterations: iteration,
                relative_residual: r_norm / b_norm,
            });
        }
        apply_m(&r, &mut z);
        let rz_new = dot(&r, &z);
        let beta = rz_new / rz;
        rz = rz_new;
        for i in 0..n {
            direction[i] = z[i] + beta * direction[i];
        }
    }
    Err(Error::NotConverged {
        solver: label.to_string(),
        iterations: max_iterations,
        residual: r_norm / b_norm,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn laplacian_1d(n: usize) -> CsrMatrix {
        let mut t = Vec::new();
        for i in 0..n {
            t.push((i, i, 2.0));
            if i + 1 < n {
                t.push((i, i + 1, -1.0));
                t.push((i + 1, i, -1.0));
            }
        }
        CsrMatrix::from_triplets(n, n, &t)
    }

    #[test]
    fn test_ic0_is_exact_for_tridiagonal_matrix() {
        // A tridiagonal SPD matrix has no fill-in, so IC(0) is the exact Cholesky factor.
        let a = laplacian_1d(8);
        let factor = IncompleteCholeskyFactor::new(&a, 0.0).expect("SPD factor");
        let x_true: Vec<f64> = (0..8).map(|i| (i as f64).sin() + 1.0).collect();
        let b = a.mul_vec(&x_true);
        let mut x = vec![0.0; 8];
        factor.solve(&b, &mut x);
        for (xi, ti) in x.iter().zip(&x_true) {
            assert!((xi - ti).abs() < 1e-12, "{xi} vs {ti}");
        }
    }

    #[test]
    fn test_ic0_reports_breakdown_for_indefinite_matrix() {
        let a =
            CsrMatrix::from_triplets(2, 2, &[(0, 0, 1.0), (0, 1, 2.0), (1, 0, 2.0), (1, 1, 1.0)]);
        assert!(IncompleteCholeskyFactor::new(&a, 0.0).is_none());
        // The SpdPreconditioner still yields a usable (shifted IC or Jacobi) operator.
        let pc = SpdPreconditioner::new(&a, &[1.0, 1.0]);
        let mut z = [0.0; 2];
        pc.apply(&[1.0, 1.0], &mut z);
        assert!(z.iter().all(|v| v.is_finite()));
    }

    #[test]
    fn test_preconditioned_cg_solves_spd_system() {
        let a = laplacian_1d(30);
        let x_true: Vec<f64> = (0..30).map(|i| 0.1 * i as f64).collect();
        let b = a.mul_vec(&x_true);
        let mut x = vec![0.0; 30];
        let outcome = preconditioned_cg(
            |v, y| {
                y.copy_from_slice(&a.mul_vec(v));
                Ok(())
            },
            |r, z| z.copy_from_slice(r),
            &b,
            &mut x,
            1e-12,
            100,
            "test CG",
        )
        .expect("converges");
        // Unpreconditioned CG terminates in at most n steps in exact arithmetic.
        assert!(outcome.iterations <= 31);
        assert!(outcome.relative_residual <= 1e-12);
        for (xi, ti) in x.iter().zip(&x_true) {
            assert!((xi - ti).abs() < 1e-9);
        }
    }

    #[test]
    fn test_preconditioned_cg_reports_non_convergence() {
        let a = laplacian_1d(50);
        let b = vec![1.0; 50];
        let mut x = vec![0.0; 50];
        let result = preconditioned_cg(
            |v, y| {
                y.copy_from_slice(&a.mul_vec(v));
                Ok(())
            },
            |r, z| z.copy_from_slice(r),
            &b,
            &mut x,
            1e-14,
            3,
            "test CG",
        );
        assert!(matches!(
            result,
            Err(Error::NotConverged { iterations: 3, .. })
        ));
    }

    #[test]
    fn test_mask_matrix_eliminates_constrained_rows_and_columns() {
        let a = laplacian_1d(3);
        let free = [true, false, true];
        let masked = mask_matrix(&a, &free, &free, true);
        assert_eq!(masked.get(1, 1), 1.0);
        assert_eq!(masked.get(0, 1), 0.0);
        assert_eq!(masked.get(1, 0), 0.0);
        assert_eq!(masked.get(0, 0), 2.0);
        assert_eq!(masked.get(2, 2), 2.0);
        let rect = mask_matrix(&a, &free, &free, false);
        assert_eq!(rect.get(1, 1), 0.0);
    }

    #[test]
    fn test_csr_triplets_checked_rejects_malformed_arrays() {
        let mut a = laplacian_1d(3);
        assert_eq!(csr_triplets_checked(&a).expect("valid").len(), a.nnz());
        a.col_indices[0] = 9;
        assert!(matches!(
            csr_triplets_checked(&a),
            Err(Error::DofOutOfRange { dof: 9, n: 3 })
        ));
        let mut b = laplacian_1d(3);
        b.row_ptr.pop();
        assert!(matches!(
            csr_triplets_checked(&b),
            Err(Error::InvalidInput(_))
        ));
    }
}
