// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! CSR × CSR sparse matrix multiplication (SpMM) and Galerkin triple product.

use crate::parallel_solver::CsrMatrix;
use std::collections::HashSet;

/// Sparse matrix-matrix multiplication: C = A × B.
///
/// Uses a two-pass algorithm:
/// - Pass 1 (symbolic): compute the column-index sets for each row of C.
/// - Pass 2 (numeric): accumulate values using a dense temporary array.
///
/// Column indices in each output row are sorted.
pub fn spmm(a: &CsrMatrix, b: &CsrMatrix) -> CsrMatrix {
    assert_eq!(a.ncols, b.nrows, "spmm: inner dimensions must match");

    let n_rows = a.nrows;
    let n_cols = b.ncols;

    // Pass 1 (symbolic): collect column index sets per row
    let mut row_col_sets: Vec<HashSet<usize>> = (0..n_rows).map(|_| HashSet::new()).collect();
    for (i, row_set) in row_col_sets.iter_mut().enumerate() {
        for k in a.row_offsets[i]..a.row_offsets[i + 1] {
            let ak_col = a.col_indices[k]; // row in B
            for bk in b.row_offsets[ak_col]..b.row_offsets[ak_col + 1] {
                row_set.insert(b.col_indices[bk]);
            }
        }
    }

    // Build row_offsets and sorted col_indices
    let mut row_offsets = vec![0usize; n_rows + 1];
    for i in 0..n_rows {
        row_offsets[i + 1] = row_offsets[i] + row_col_sets[i].len();
    }
    let nnz = row_offsets[n_rows];

    let mut col_indices = vec![0usize; nnz];
    // For each row, sort the column indices
    for i in 0..n_rows {
        let base = row_offsets[i];
        let mut cols: Vec<usize> = row_col_sets[i].iter().copied().collect();
        cols.sort_unstable();
        for (j, &c) in cols.iter().enumerate() {
            col_indices[base + j] = c;
        }
    }

    // Pass 2 (numeric): compute values using a dense accumulator
    let mut values = vec![0.0f64; nnz];
    // Dense accumulator and "touched" marker to avoid clearing the full array
    let mut dense = vec![0.0f64; n_cols];
    let mut touched: Vec<usize> = Vec::new();

    for i in 0..n_rows {
        // Accumulate contributions
        for k in a.row_offsets[i]..a.row_offsets[i + 1] {
            let a_val = a.values[k];
            let ak_col = a.col_indices[k];
            for bk in b.row_offsets[ak_col]..b.row_offsets[ak_col + 1] {
                let b_col = b.col_indices[bk];
                if dense[b_col] == 0.0 {
                    touched.push(b_col);
                }
                dense[b_col] += a_val * b.values[bk];
            }
        }

        // Scatter from dense into the CSR values array
        for k in row_offsets[i]..row_offsets[i + 1] {
            values[k] = dense[col_indices[k]];
        }

        // Reset touched entries
        for &c in &touched {
            dense[c] = 0.0;
        }
        touched.clear();
    }

    CsrMatrix {
        nrows: n_rows,
        ncols: n_cols,
        row_offsets,
        col_indices,
        values,
    }
}

/// Transpose a CSR matrix.
///
/// Returns `A^T` in CSR format with sorted column indices per row.
pub fn csr_transpose(a: &CsrMatrix) -> CsrMatrix {
    let n_rows = a.nrows;
    let n_cols = a.ncols;
    let nnz = a.nnz();

    // Count entries per column of A (= per row of A^T)
    let mut row_counts = vec![0usize; n_cols];
    for &c in &a.col_indices {
        row_counts[c] += 1;
    }

    // Build row_offsets for A^T
    let mut row_offsets = vec![0usize; n_cols + 1];
    for i in 0..n_cols {
        row_offsets[i + 1] = row_offsets[i] + row_counts[i];
    }

    // Fill col_indices and values of A^T
    let mut col_indices = vec![0usize; nnz];
    let mut values = vec![0.0f64; nnz];
    let mut write_pos = row_offsets.clone();

    for i in 0..n_rows {
        for k in a.row_offsets[i]..a.row_offsets[i + 1] {
            let j = a.col_indices[k];
            let pos = write_pos[j];
            col_indices[pos] = i;
            values[pos] = a.values[k];
            write_pos[j] += 1;
        }
    }

    // Sort each row of A^T by column index (ensures canonical ordering)
    for r in 0..n_cols {
        let rs = row_offsets[r];
        let re = row_offsets[r + 1];
        if re > rs + 1 {
            // Sort the slice by col_indices using a permutation
            let mut perm: Vec<usize> = (rs..re).collect();
            perm.sort_unstable_by_key(|&k| col_indices[k]);
            // Apply permutation in-place
            let cols_copy: Vec<usize> = col_indices[rs..re].to_vec();
            let vals_copy: Vec<f64> = values[rs..re].to_vec();
            for (out_k, &in_k) in perm.iter().enumerate() {
                col_indices[rs + out_k] = cols_copy[in_k - rs];
                values[rs + out_k] = vals_copy[in_k - rs];
            }
        }
    }

    CsrMatrix {
        nrows: n_cols,
        ncols: n_rows,
        row_offsets,
        col_indices,
        values,
    }
}

/// Compute the Galerkin coarse-level operator: A_c = P^T × A × P.
///
/// Steps:
/// 1. Compute `AP = A * P`
/// 2. Compute `A_c = P^T * AP`
pub fn galerkin_coarse(a: &CsrMatrix, p: &CsrMatrix) -> CsrMatrix {
    let ap = spmm(a, p);
    let pt = csr_transpose(p);
    spmm(&pt, &ap)
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// Extract value at (row, col) from a CsrMatrix.
    fn get_val(m: &CsrMatrix, row: usize, col: usize) -> f64 {
        for k in m.row_offsets[row]..m.row_offsets[row + 1] {
            if m.col_indices[k] == col {
                return m.values[k];
            }
        }
        0.0
    }

    /// Build a small tridiagonal CsrMatrix for testing.
    fn make_tridiag(n: usize, diag: f64, off: f64) -> CsrMatrix {
        let mut row_offsets = vec![0usize; n + 1];
        let mut col_indices = Vec::new();
        let mut values = Vec::new();

        for i in 0..n {
            if i > 0 {
                col_indices.push(i - 1);
                values.push(off);
            }
            col_indices.push(i);
            values.push(diag);
            if i + 1 < n {
                col_indices.push(i + 1);
                values.push(off);
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

    #[test]
    fn test_spmm_identity() {
        let n = 6;
        let a = make_tridiag(n, 2.0, -1.0);
        let id = CsrMatrix::identity(n);

        // I * A should equal A
        let ia = spmm(&id, &a);
        assert_eq!(ia.nrows, a.nrows);
        assert_eq!(ia.ncols, a.ncols);
        for i in 0..n {
            for k in a.row_offsets[i]..a.row_offsets[i + 1] {
                let j = a.col_indices[k];
                let expected = a.values[k];
                let got = get_val(&ia, i, j);
                assert!(
                    (got - expected).abs() < 1e-13,
                    "I*A[{i},{j}] = {got}, expected {expected}"
                );
            }
        }

        // A * I should equal A
        let ai = spmm(&a, &id);
        for i in 0..n {
            for k in a.row_offsets[i]..a.row_offsets[i + 1] {
                let j = a.col_indices[k];
                let expected = a.values[k];
                let got = get_val(&ai, i, j);
                assert!(
                    (got - expected).abs() < 1e-13,
                    "A*I[{i},{j}] = {got}, expected {expected}"
                );
            }
        }
    }

    #[test]
    fn test_galerkin_preserves_symmetry() {
        // A = 4x4 tridiagonal SPD
        let n = 4;
        let a = make_tridiag(n, 2.0, -1.0);

        // P = 4x2 prolongation from coarse (C-points 0,2) to fine
        // C-points: 0 and 2 (direct injection); F-points: 1 and 3 (averaging)
        let p = CsrMatrix {
            nrows: 4,
            ncols: 2,
            row_offsets: vec![0, 1, 3, 4, 5],
            col_indices: vec![0, 0, 1, 1, 1],
            values: vec![1.0, 0.5, 0.5, 1.0, 0.5],
        };

        let ac = galerkin_coarse(&a, &p);
        assert_eq!(ac.nrows, 2);
        assert_eq!(ac.ncols, 2);

        // Check symmetry: ac[0,1] == ac[1,0]
        let v01 = get_val(&ac, 0, 1);
        let v10 = get_val(&ac, 1, 0);
        assert!(
            (v01 - v10).abs() < 1e-12,
            "Galerkin product not symmetric: ac[0,1]={v01}, ac[1,0]={v10}"
        );

        // Diagonal should be positive (SPD preserved for SPD input with valid P)
        let v00 = get_val(&ac, 0, 0);
        let v11 = get_val(&ac, 1, 1);
        assert!(v00 > 0.0, "ac[0,0] should be positive, got {v00}");
        assert!(v11 > 0.0, "ac[1,1] should be positive, got {v11}");
    }

    #[test]
    fn test_csr_transpose_involution() {
        let n = 5;
        let a = make_tridiag(n, 3.0, -1.0);
        let at = csr_transpose(&a);
        let att = csr_transpose(&at);

        // (A^T)^T should equal A
        assert_eq!(att.nrows, a.nrows);
        assert_eq!(att.ncols, a.ncols);
        for i in 0..n {
            for k in a.row_offsets[i]..a.row_offsets[i + 1] {
                let j = a.col_indices[k];
                let expected = a.values[k];
                let got = get_val(&att, i, j);
                assert!(
                    (got - expected).abs() < 1e-13,
                    "(A^T)^T[{i},{j}] = {got}, expected {expected}"
                );
            }
        }
    }

    #[test]
    fn test_spmm_small_explicit() {
        // A = [[1,2],[3,4]], B = [[5,6],[7,8]]
        // AB = [[19,22],[43,50]]
        let a = CsrMatrix {
            nrows: 2,
            ncols: 2,
            row_offsets: vec![0, 2, 4],
            col_indices: vec![0, 1, 0, 1],
            values: vec![1.0, 2.0, 3.0, 4.0],
        };
        let b = CsrMatrix {
            nrows: 2,
            ncols: 2,
            row_offsets: vec![0, 2, 4],
            col_indices: vec![0, 1, 0, 1],
            values: vec![5.0, 6.0, 7.0, 8.0],
        };
        let ab = spmm(&a, &b);
        assert!((get_val(&ab, 0, 0) - 19.0).abs() < 1e-13);
        assert!((get_val(&ab, 0, 1) - 22.0).abs() < 1e-13);
        assert!((get_val(&ab, 1, 0) - 43.0).abs() < 1e-13);
        assert!((get_val(&ab, 1, 1) - 50.0).abs() < 1e-13);
    }
}
