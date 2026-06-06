// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Gauss-Seidel smoothers for AMG cycles.

use crate::parallel_solver::CsrMatrix;

/// Forward Gauss-Seidel: sweep rows 0..nrows in ascending order.
///
/// For each row `i`: `x[i] = (b[i] - Σ_{j≠i} A[i,j]*x[j]) / A[i,i]`
///
/// Updates are applied in-place so later rows benefit from earlier updates.
pub fn forward_gs(a: &CsrMatrix, b: &[f64], x: &mut [f64], n_sweeps: usize) {
    let n = a.nrows;
    for _ in 0..n_sweeps {
        for i in 0..n {
            let rs = a.row_offsets[i];
            let re = a.row_offsets[i + 1];
            let mut sigma = 0.0f64;
            let mut diag = 1.0f64;
            for k in rs..re {
                let j = a.col_indices[k];
                if j == i {
                    diag = a.values[k];
                } else {
                    sigma += a.values[k] * x[j];
                }
            }
            if diag.abs() > 1e-300 {
                x[i] = (b[i] - sigma) / diag;
            }
        }
    }
}

/// Backward Gauss-Seidel: sweep rows nrows-1..0 in descending order.
pub fn backward_gs(a: &CsrMatrix, b: &[f64], x: &mut [f64], n_sweeps: usize) {
    let n = a.nrows;
    for _ in 0..n_sweeps {
        for i in (0..n).rev() {
            let rs = a.row_offsets[i];
            let re = a.row_offsets[i + 1];
            let mut sigma = 0.0f64;
            let mut diag = 1.0f64;
            for k in rs..re {
                let j = a.col_indices[k];
                if j == i {
                    diag = a.values[k];
                } else {
                    sigma += a.values[k] * x[j];
                }
            }
            if diag.abs() > 1e-300 {
                x[i] = (b[i] - sigma) / diag;
            }
        }
    }
}

/// Symmetric Gauss-Seidel: one forward sweep followed by one backward sweep per call.
///
/// Each call to `symmetric_gs` with `n_sweeps` performs `n_sweeps` complete SGS
/// sweeps (each consisting of one forward + one backward pass).
pub fn symmetric_gs(a: &CsrMatrix, b: &[f64], x: &mut [f64], n_sweeps: usize) {
    for _ in 0..n_sweeps {
        forward_gs(a, b, x, 1);
        backward_gs(a, b, x, 1);
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// Build 1D Poisson tridiagonal of size n.
    fn make_1d_poisson(n: usize) -> CsrMatrix {
        let mut row_offsets = vec![0usize; n + 1];
        let mut col_indices = Vec::new();
        let mut values = Vec::new();

        for i in 0..n {
            if i > 0 {
                col_indices.push(i - 1);
                values.push(-1.0);
            }
            col_indices.push(i);
            values.push(2.0);
            if i + 1 < n {
                col_indices.push(i + 1);
                values.push(-1.0);
            }
            row_offsets[i + 1] = col_indices.len();
        }

        CsrMatrix {
            nrows: n,
            ncols: n,
            row_offsets,
            col_indices,
            values,
        }
    }

    fn residual_norm(a: &CsrMatrix, b: &[f64], x: &[f64]) -> f64 {
        let n = a.nrows;
        let mut r = vec![0.0f64; n];
        a.spmv(x, &mut r);
        let mut norm = 0.0f64;
        for i in 0..n {
            let ri = b[i] - r[i];
            norm += ri * ri;
        }
        norm.sqrt()
    }

    #[test]
    fn test_gs_reduces_residual_1d_poisson() {
        let n = 5;
        let a = make_1d_poisson(n);
        let b = vec![1.0f64; n];
        let mut x = vec![0.0f64; n];

        let r0 = residual_norm(&a, &b, &x);
        symmetric_gs(&a, &b, &mut x, 1);
        let r1 = residual_norm(&a, &b, &x);

        assert!(
            r1 < r0,
            "SGS did not reduce residual: before={r0}, after={r1}"
        );
    }

    #[test]
    fn test_gs_converges_tiny() {
        // 3×3 strictly diagonally dominant system
        // A = [[4,-1,0],[-1,4,-1],[0,-1,4]], b = [1,1,1]
        // Row offsets: row 0 has 2 entries, row 1 has 3 entries, row 2 has 2 entries
        let a = CsrMatrix {
            nrows: 3,
            ncols: 3,
            row_offsets: vec![0, 2, 5, 7],
            col_indices: vec![0, 1, 0, 1, 2, 1, 2],
            values: vec![4.0, -1.0, -1.0, 4.0, -1.0, -1.0, 4.0],
        };
        let b = vec![1.0f64, 1.0, 1.0];
        let mut x = vec![0.0f64; 3];

        // Run 20 GS sweeps
        symmetric_gs(&a, &b, &mut x, 20);
        let r = residual_norm(&a, &b, &x);
        assert!(
            r < 1e-10,
            "GS did not converge in 20 sweeps: residual = {r}"
        );
    }

    #[test]
    fn test_forward_and_backward_reduce_residual() {
        let n = 8;
        let a = make_1d_poisson(n);
        let b: Vec<f64> = (0..n).map(|i| (i + 1) as f64 * 0.1).collect();
        let mut x_fwd = vec![0.0f64; n];
        let mut x_bwd = vec![0.0f64; n];

        let r0 = residual_norm(&a, &b, &x_fwd);
        forward_gs(&a, &b, &mut x_fwd, 5);
        backward_gs(&a, &b, &mut x_bwd, 5);
        let r_fwd = residual_norm(&a, &b, &x_fwd);
        let r_bwd = residual_norm(&a, &b, &x_bwd);

        assert!(r_fwd < r0, "Forward GS did not reduce: {r0} -> {r_fwd}");
        assert!(r_bwd < r0, "Backward GS did not reduce: {r0} -> {r_bwd}");
    }
}
