//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
use super::functions::invert_small_dense;
#[allow(unused_imports)]
use super::functions::*;
#[allow(unused_imports)]
use super::functions_2::*;

/// Sparse LU factorization result.
///
/// Stores L and U factors as CSR matrices.
pub struct SparseLu {
    /// Lower triangular factor (unit diagonal).
    pub l: CsrMatrix,
    /// Upper triangular factor.
    pub u: CsrMatrix,
    /// Pivot permutation vector (row i of A is in perm\[i\]).
    pub perm: Vec<usize>,
}
impl SparseLu {
    /// Compute an ILU(0) factorization (incomplete LU with no fill-in).
    ///
    /// This is a preconditioner-quality factorization.
    /// The sparsity pattern of L+U matches that of A.
    pub fn ilu0(a: &CsrMatrix) -> Self {
        let n = a.nrows;
        assert_eq!(a.nrows, a.ncols, "ILU(0) requires square matrix");
        let mut lval = vec![vec![0.0f64; n]; n];
        for r in 0..n {
            for k in a.row_ptr[r]..a.row_ptr[r + 1] {
                lval[r][a.col_idx[k]] = a.values[k];
            }
        }
        for k in 0..n {
            let akk = lval[k][k];
            if akk.abs() < 1e-30 {
                continue;
            }
            for i in (k + 1)..n {
                if lval[i][k].abs() > 0.0 {
                    lval[i][k] /= akk;
                    for j in (k + 1)..n {
                        if lval[i][j].abs() > 0.0 || lval[k][j].abs() > 0.0 {
                            lval[i][j] -= lval[i][k] * lval[k][j];
                        }
                    }
                }
            }
        }
        let mut l_rows = Vec::new();
        let mut l_cols = Vec::new();
        let mut l_vals = Vec::new();
        let mut u_rows = Vec::new();
        let mut u_cols = Vec::new();
        let mut u_vals = Vec::new();
        for i in 0..n {
            l_rows.push(i);
            l_cols.push(i);
            l_vals.push(1.0);
            for j in 0..i {
                if lval[i][j].abs() > 1e-30 {
                    l_rows.push(i);
                    l_cols.push(j);
                    l_vals.push(lval[i][j]);
                }
            }
            for j in i..n {
                if lval[i][j].abs() > 1e-30 {
                    u_rows.push(i);
                    u_cols.push(j);
                    u_vals.push(lval[i][j]);
                }
            }
        }
        let l = CsrMatrix::from_triplets(n, n, &l_rows, &l_cols, &l_vals);
        let u = CsrMatrix::from_triplets(n, n, &u_rows, &u_cols, &u_vals);
        let perm: Vec<usize> = (0..n).collect();
        SparseLu { l, u, perm }
    }
    /// Solve A x = b using the LU factors via forward/backward substitution.
    pub fn solve(&self, b: &[f64]) -> Vec<f64> {
        let n = b.len();
        let mut y = vec![0.0f64; n];
        for i in 0..n {
            let mut s = b[i];
            for k in self.l.row_ptr[i]..self.l.row_ptr[i + 1] {
                let j = self.l.col_idx[k];
                if j < i {
                    s -= self.l.values[k] * y[j];
                }
            }
            y[i] = s;
        }
        let mut x = vec![0.0f64; n];
        for i in (0..n).rev() {
            let mut s = y[i];
            let mut diag = 1.0;
            for k in self.u.row_ptr[i]..self.u.row_ptr[i + 1] {
                let j = self.u.col_idx[k];
                if j > i {
                    s -= self.u.values[k] * x[j];
                } else if j == i {
                    diag = self.u.values[k];
                }
            }
            x[i] = if diag.abs() < 1e-30 { 0.0 } else { s / diag };
        }
        x
    }
}
/// Sparse matrix in Compressed Sparse Column (CSC) format.
///
/// Dual of CSR: `col_ptr[j]..col_ptr[j+1]` index the entries for column `j`.
pub struct CscMatrix {
    /// Number of rows.
    pub nrows: usize,
    /// Number of columns.
    pub ncols: usize,
    /// Column pointer array (length `ncols + 1`).
    pub col_ptr: Vec<usize>,
    /// Row indices of stored entries.
    pub row_idx: Vec<usize>,
    /// Values of stored entries.
    pub values: Vec<f64>,
}
impl CscMatrix {
    /// Build a `CscMatrix` from a `CsrMatrix` by transposing the index structure.
    pub fn from_csr(a: &CsrMatrix) -> Self {
        let nrows = a.nrows;
        let ncols = a.ncols;
        let nnz = a.nnz();
        let mut col_counts = vec![0usize; ncols];
        for &c in &a.col_idx {
            col_counts[c] += 1;
        }
        let mut col_ptr = vec![0usize; ncols + 1];
        for j in 0..ncols {
            col_ptr[j + 1] = col_ptr[j] + col_counts[j];
        }
        let mut row_idx = vec![0usize; nnz];
        let mut values = vec![0.0f64; nnz];
        let mut pos = col_ptr.clone();
        for i in 0..nrows {
            for k in a.row_ptr[i]..a.row_ptr[i + 1] {
                let j = a.col_idx[k];
                let p = pos[j];
                row_idx[p] = i;
                values[p] = a.values[k];
                pos[j] += 1;
            }
        }
        CscMatrix {
            nrows,
            ncols,
            col_ptr,
            row_idx,
            values,
        }
    }
    /// Return all (row, value) pairs stored in column `j`.
    pub fn get_column(&self, j: usize) -> Vec<(usize, f64)> {
        assert!(j < self.ncols, "column index out of bounds");
        let start = self.col_ptr[j];
        let end = self.col_ptr[j + 1];
        (start..end)
            .map(|k| (self.row_idx[k], self.values[k]))
            .collect()
    }
    /// Sparse matrix-vector multiply: y = A * x (CSC format).
    pub fn matvec(&self, x: &[f64]) -> Vec<f64> {
        assert_eq!(x.len(), self.ncols, "matvec: dimension mismatch");
        let mut y = vec![0.0f64; self.nrows];
        for j in 0..self.ncols {
            for k in self.col_ptr[j]..self.col_ptr[j + 1] {
                y[self.row_idx[k]] += self.values[k] * x[j];
            }
        }
        y
    }
    /// Number of stored non-zeros.
    pub fn nnz(&self) -> usize {
        self.values.len()
    }
}
/// Block Jacobi preconditioner: splits the matrix into diagonal blocks and
/// inverts each block independently.
///
/// `block_size` must divide `n`.
pub struct BlockJacobi {
    /// Inverse of each diagonal block, stored as flat row-major matrices.
    #[allow(dead_code)]
    pub block_invs: Vec<Vec<f64>>,
    /// Dimension of each block.
    #[allow(dead_code)]
    pub block_size: usize,
    /// Number of rows / columns.
    #[allow(dead_code)]
    pub n: usize,
}
#[allow(dead_code)]
impl BlockJacobi {
    /// Build a block Jacobi preconditioner from a CSR matrix.
    pub fn new(a: &CsrMatrix, block_size: usize) -> Self {
        let n = a.nrows;
        assert_eq!(n % block_size, 0, "block_size must divide matrix dimension");
        let n_blocks = n / block_size;
        let mut block_invs = Vec::with_capacity(n_blocks);
        for b in 0..n_blocks {
            let start = b * block_size;
            let mut block = vec![0.0f64; block_size * block_size];
            for i in 0..block_size {
                for j in 0..block_size {
                    block[i * block_size + j] = a.get(start + i, start + j);
                }
            }
            let inv = invert_small_dense(&block, block_size).unwrap_or(block);
            block_invs.push(inv);
        }
        BlockJacobi {
            block_invs,
            block_size,
            n,
        }
    }
    /// Apply the block Jacobi preconditioner: x = M^{-1} b.
    pub fn apply(&self, b: &[f64]) -> Vec<f64> {
        let mut x = vec![0.0f64; self.n];
        let bs = self.block_size;
        for (blk, inv) in self.block_invs.iter().enumerate() {
            let start = blk * bs;
            for i in 0..bs {
                let mut s = 0.0f64;
                for j in 0..bs {
                    s += inv[i * bs + j] * b[start + j];
                }
                x[start + i] = s;
            }
        }
        x
    }
}
/// Sparse matrix in Compressed Sparse Row (CSR) format.
///
/// # Layout
/// - `row_ptr[i]..row_ptr[i+1]` are the indices into `col_idx` / `values` for row `i`.
/// - `col_idx[k]` and `values[k]` give the column index and value of the k-th stored entry.
pub struct CsrMatrix {
    /// Number of rows.
    pub nrows: usize,
    /// Number of columns.
    pub ncols: usize,
    /// Row pointer array (length `nrows + 1`).
    pub row_ptr: Vec<usize>,
    /// Column indices of stored entries.
    pub col_idx: Vec<usize>,
    /// Values of stored entries.
    pub values: Vec<f64>,
}
impl CsrMatrix {
    /// Create an empty sparse matrix with no stored entries.
    pub fn new(nrows: usize, ncols: usize) -> Self {
        CsrMatrix {
            nrows,
            ncols,
            row_ptr: vec![0usize; nrows + 1],
            col_idx: Vec::new(),
            values: Vec::new(),
        }
    }
    /// Build a CSR matrix from COO triplets (row, col, val).
    ///
    /// Duplicate (row, col) entries are summed.
    pub fn from_triplets(
        nrows: usize,
        ncols: usize,
        rows: &[usize],
        cols: &[usize],
        vals: &[f64],
    ) -> Self {
        assert_eq!(rows.len(), cols.len());
        assert_eq!(rows.len(), vals.len());
        let mut row_ptr = vec![0usize; nrows + 1];
        for &r in rows {
            assert!(r < nrows, "row index out of bounds");
            row_ptr[r + 1] += 1;
        }
        for i in 0..nrows {
            row_ptr[i + 1] += row_ptr[i];
        }
        let nnz = rows.len();
        let mut col_idx = vec![0usize; nnz];
        let mut values = vec![0.0f64; nnz];
        let mut pos = row_ptr[..nrows].to_vec();
        for k in 0..nnz {
            let r = rows[k];
            let p = pos[r];
            col_idx[p] = cols[k];
            values[p] = vals[k];
            pos[r] += 1;
        }
        let mut mat = CsrMatrix {
            nrows,
            ncols,
            row_ptr,
            col_idx,
            values,
        };
        mat.sort_and_sum_duplicates();
        mat
    }
    /// Sort column indices within each row and collapse duplicate (row, col) entries by summing.
    fn sort_and_sum_duplicates(&mut self) {
        let mut new_col_idx = Vec::with_capacity(self.col_idx.len());
        let mut new_values = Vec::with_capacity(self.values.len());
        let mut new_row_ptr = vec![0usize; self.nrows + 1];
        for r in 0..self.nrows {
            let start = self.row_ptr[r];
            let end = self.row_ptr[r + 1];
            let mut pairs: Vec<(usize, f64)> = self.col_idx[start..end]
                .iter()
                .copied()
                .zip(self.values[start..end].iter().copied())
                .collect();
            pairs.sort_unstable_by_key(|&(c, _)| c);
            let row_start = new_col_idx.len();
            for (c, v) in pairs {
                if let Some(last_c) = new_col_idx.last().copied()
                    && last_c == c
                    && new_col_idx.len() > row_start
                {
                    if let Some(last) = new_values.last_mut() {
                        *last += v;
                    }
                    continue;
                }
                new_col_idx.push(c);
                new_values.push(v);
            }
            new_row_ptr[r + 1] = new_col_idx.len();
        }
        self.row_ptr = new_row_ptr;
        self.col_idx = new_col_idx;
        self.values = new_values;
    }
    /// Return the value at position `(row, col)`, or 0.0 if not stored.
    pub fn get(&self, row: usize, col: usize) -> f64 {
        let start = self.row_ptr[row];
        let end = self.row_ptr[row + 1];
        for k in start..end {
            if self.col_idx[k] == col {
                return self.values[k];
            }
        }
        0.0
    }
    /// Number of stored (non-zero) entries.
    pub fn nnz(&self) -> usize {
        self.values.len()
    }
    /// Matrix-vector product: returns `y = A * x`.
    pub fn matvec(&self, x: &[f64]) -> Vec<f64> {
        assert_eq!(x.len(), self.ncols, "matvec: x length must equal ncols");
        let mut y = vec![0.0f64; self.nrows];
        for r in 0..self.nrows {
            let mut acc = 0.0;
            for k in self.row_ptr[r]..self.row_ptr[r + 1] {
                acc += self.values[k] * x[self.col_idx[k]];
            }
            y[r] = acc;
        }
        y
    }
    /// Sparse matrix addition: returns `self + other`.
    pub fn add(&self, other: &CsrMatrix) -> CsrMatrix {
        assert_eq!(self.nrows, other.nrows);
        assert_eq!(self.ncols, other.ncols);
        let mut rows = Vec::new();
        let mut cols = Vec::new();
        let mut vals = Vec::new();
        for r in 0..self.nrows {
            for k in self.row_ptr[r]..self.row_ptr[r + 1] {
                rows.push(r);
                cols.push(self.col_idx[k]);
                vals.push(self.values[k]);
            }
            for k in other.row_ptr[r]..other.row_ptr[r + 1] {
                rows.push(r);
                cols.push(other.col_idx[k]);
                vals.push(other.values[k]);
            }
        }
        CsrMatrix::from_triplets(self.nrows, self.ncols, &rows, &cols, &vals)
    }
    /// Scalar multiplication: returns `alpha * self`.
    pub fn scale(&self, alpha: f64) -> CsrMatrix {
        CsrMatrix {
            nrows: self.nrows,
            ncols: self.ncols,
            row_ptr: self.row_ptr.clone(),
            col_idx: self.col_idx.clone(),
            values: self.values.iter().map(|&v| alpha * v).collect(),
        }
    }
    /// Returns the transpose of `self`.
    pub fn transpose(&self) -> CsrMatrix {
        let rows: Vec<usize> = (0..self.nrows)
            .flat_map(|r| {
                let start = self.row_ptr[r];
                let end = self.row_ptr[r + 1];
                vec![r; end - start]
            })
            .collect();
        let cols = self.col_idx.clone();
        let vals = self.values.clone();
        CsrMatrix::from_triplets(self.ncols, self.nrows, &cols, &rows, &vals)
    }
    /// Extract the main diagonal. Returns a vector of length `min(nrows, ncols)`.
    pub fn diagonal(&self) -> Vec<f64> {
        let n = self.nrows.min(self.ncols);
        let mut d = vec![0.0f64; n];
        for i in 0..n {
            d[i] = self.get(i, i);
        }
        d
    }
    /// Convert to dense row-major representation.
    pub fn to_dense(&self) -> Vec<Vec<f64>> {
        let mut dense = vec![vec![0.0f64; self.ncols]; self.nrows];
        for r in 0..self.nrows {
            for k in self.row_ptr[r]..self.row_ptr[r + 1] {
                dense[r][self.col_idx[k]] = self.values[k];
            }
        }
        dense
    }
    /// Build a CSR matrix from a dense 2-D slice (rows of equal length).
    pub fn from_dense(a: &[Vec<f64>]) -> Self {
        let nrows = a.len();
        let ncols = if nrows == 0 { 0 } else { a[0].len() };
        let mut rows = Vec::new();
        let mut cols = Vec::new();
        let mut vals = Vec::new();
        for (r, row) in a.iter().enumerate() {
            for (c, &v) in row.iter().enumerate() {
                if v != 0.0 {
                    rows.push(r);
                    cols.push(c);
                    vals.push(v);
                }
            }
        }
        CsrMatrix::from_triplets(nrows, ncols, &rows, &cols, &vals)
    }
    /// Check if the matrix is symmetric within `tol`.
    pub fn is_symmetric(&self, tol: f64) -> bool {
        if self.nrows != self.ncols {
            return false;
        }
        for r in 0..self.nrows {
            for k in self.row_ptr[r]..self.row_ptr[r + 1] {
                let c = self.col_idx[k];
                let v = self.values[k];
                if (v - self.get(c, r)).abs() > tol {
                    return false;
                }
            }
        }
        true
    }
    /// Compute the Kronecker product `A ⊗ B` of two CSR matrices.
    ///
    /// The result is an `(m*p) × (n*q)` matrix where `A` is `m×n` and `B` is `p×q`.
    #[allow(dead_code)]
    pub fn kronecker_product(&self, b: &CsrMatrix) -> CsrMatrix {
        let (m, n) = (self.nrows, self.ncols);
        let (p, q) = (b.nrows, b.ncols);
        let nrows = m * p;
        let ncols = n * q;
        let mut rows: Vec<usize> = Vec::new();
        let mut cols: Vec<usize> = Vec::new();
        let mut vals: Vec<f64> = Vec::new();
        for ra in 0..m {
            for ka in self.row_ptr[ra]..self.row_ptr[ra + 1] {
                let ca = self.col_idx[ka];
                let va = self.values[ka];
                for rb in 0..p {
                    for kb in b.row_ptr[rb]..b.row_ptr[rb + 1] {
                        let cb = b.col_idx[kb];
                        let vb = b.values[kb];
                        rows.push(ra * p + rb);
                        cols.push(ca * q + cb);
                        vals.push(va * vb);
                    }
                }
            }
        }
        CsrMatrix::from_triplets(nrows, ncols, &rows, &cols, &vals)
    }
    /// Compute the IC(0) incomplete Cholesky factorization of a symmetric
    /// positive-definite CSR matrix.
    ///
    /// Returns the lower-triangular factor `L` (as a CSR matrix with the same
    /// sparsity pattern as the lower triangle of `self`) such that `L * L^T ≈ A`.
    ///
    /// # Panics
    /// Panics if the matrix is not square, or if a zero pivot is encountered.
    #[allow(dead_code)]
    pub fn incomplete_cholesky(&self) -> CsrMatrix {
        assert_eq!(self.nrows, self.ncols, "IC(0): matrix must be square");
        let n = self.nrows;
        let mut diag = vec![0.0f64; n];
        let mut lower: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
        for r in 0..n {
            for k in self.row_ptr[r]..self.row_ptr[r + 1] {
                let c = self.col_idx[k];
                let v = self.values[k];
                if c == r {
                    diag[r] += v;
                } else if c < r {
                    lower[r].push((c, v));
                }
            }
        }
        let mut l_diag = vec![0.0f64; n];
        let mut l_lower: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
        for i in 0..n {
            let mut sum_sq = 0.0f64;
            for &(j, _) in &lower[i] {
                let l_ij = l_lower[i]
                    .iter()
                    .find(|&&(c, _)| c == j)
                    .map_or(0.0, |&(_, v)| v);
                sum_sq += l_ij * l_ij;
            }
            let pivot = diag[i] - sum_sq;
            assert!(pivot > 0.0, "IC(0): non-positive pivot at row {}", i);
            l_diag[i] = pivot.sqrt();
            for k in (i + 1)..n {
                let a_ki = self.get(k, i);
                if a_ki.abs() < 1e-300 {
                    continue;
                }
                let mut dot = 0.0f64;
                for &(j, l_ij) in &l_lower[i] {
                    let l_kj = l_lower[k]
                        .iter()
                        .find(|&&(c, _)| c == j)
                        .map_or(0.0, |&(_, v)| v);
                    dot += l_kj * l_ij;
                }
                let l_ki = (a_ki - dot) / l_diag[i];
                l_lower[k].push((i, l_ki));
            }
        }
        let mut rows: Vec<usize> = Vec::new();
        let mut cols: Vec<usize> = Vec::new();
        let mut vals: Vec<f64> = Vec::new();
        for i in 0..n {
            rows.push(i);
            cols.push(i);
            vals.push(l_diag[i]);
            for &(j, v) in &l_lower[i] {
                rows.push(i);
                cols.push(j);
                vals.push(v);
            }
        }
        CsrMatrix::from_triplets(n, n, &rows, &cols, &vals)
    }
    /// Estimate the dominant eigenvalue using power iteration.
    ///
    /// Runs for at most `max_iter` iterations and stops when the relative change
    /// in the Rayleigh quotient drops below `tol`.
    ///
    /// Returns `(eigenvalue, eigenvector)`.  The eigenvector has unit L2 norm.
    ///
    /// # Panics
    /// Panics if the matrix is not square or if `n == 0`.
    #[allow(dead_code)]
    pub fn power_iteration(&self, max_iter: usize, tol: f64) -> (f64, Vec<f64>) {
        assert_eq!(
            self.nrows, self.ncols,
            "power_iteration: matrix must be square"
        );
        let n = self.nrows;
        assert!(n > 0, "power_iteration: matrix must be non-empty");
        let mut v: Vec<f64> = vec![1.0 / (n as f64).sqrt(); n];
        let mut lambda = 0.0f64;
        for _ in 0..max_iter {
            let av = self.matvec(&v);
            let new_lambda: f64 = v.iter().zip(av.iter()).map(|(&vi, &avi)| vi * avi).sum();
            let norm: f64 = av.iter().map(|&x| x * x).sum::<f64>().sqrt();
            if norm < 1e-300 {
                break;
            }
            v = av.into_iter().map(|x| x / norm).collect();
            if (new_lambda - lambda).abs() < tol * (lambda.abs().max(1.0)) {
                lambda = new_lambda;
                break;
            }
            lambda = new_lambda;
        }
        (lambda, v)
    }
}
/// Sparse QR factorization result using Householder reflections.
///
/// Stores Q and R as dense matrices (suitable for small-to-medium problems).
pub struct SparseQr {
    /// Orthogonal matrix Q (nrows × nrows).
    pub q: Vec<Vec<f64>>,
    /// Upper triangular matrix R (nrows × ncols).
    pub r: Vec<Vec<f64>>,
    /// Number of rows.
    pub nrows: usize,
    /// Number of columns.
    pub ncols: usize,
}
impl SparseQr {
    /// Compute a QR factorization of the sparse matrix `a`.
    ///
    /// Uses Householder reflections on the dense representation.
    pub fn factor(a: &CsrMatrix) -> Self {
        let m = a.nrows;
        let n = a.ncols;
        let mut q: Vec<Vec<f64>> = (0..m)
            .map(|i| (0..m).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
            .collect();
        let mut r: Vec<Vec<f64>> = a.to_dense();
        while r.len() < m {
            r.push(vec![0.0; n]);
        }
        for k in 0..m.min(n) {
            let col_len = m - k;
            let mut v: Vec<f64> = (k..m).map(|i| r[i][k]).collect();
            let norm_v: f64 = v.iter().map(|x| x * x).sum::<f64>().sqrt();
            if norm_v < 1e-14 {
                continue;
            }
            let sign = if v[0] >= 0.0 { 1.0 } else { -1.0 };
            v[0] += sign * norm_v;
            let norm_v2: f64 = v.iter().map(|x| x * x).sum::<f64>();
            if norm_v2 < 1e-28 {
                continue;
            }
            for j in k..n {
                let dot: f64 = (0..col_len).map(|i| v[i] * r[i + k][j]).sum();
                let factor = 2.0 * dot / norm_v2;
                for i in 0..col_len {
                    r[i + k][j] -= factor * v[i];
                }
            }
            for i in 0..m {
                let dot: f64 = (0..col_len).map(|l| q[i][l + k] * v[l]).sum();
                let factor = 2.0 * dot / norm_v2;
                for l in 0..col_len {
                    q[i][l + k] -= factor * v[l];
                }
            }
        }
        SparseQr {
            q,
            r,
            nrows: m,
            ncols: n,
        }
    }
    /// Solve the least-squares problem A x = b using QR.
    ///
    /// Returns the least-squares solution x of length `ncols`.
    pub fn solve_least_squares(&self, b: &[f64]) -> Vec<f64> {
        let m = self.nrows;
        let n = self.ncols;
        let mut c = vec![0.0f64; m];
        for i in 0..m {
            for j in 0..m.min(b.len()) {
                c[i] += self.q[j][i] * b[j];
            }
        }
        let mut x = vec![0.0f64; n];
        for i in (0..n).rev() {
            let mut s = c[i];
            for j in (i + 1)..n {
                s -= self.r[i][j] * x[j];
            }
            let diag = self.r[i][i];
            x[i] = if diag.abs() < 1e-14 { 0.0 } else { s / diag };
        }
        x
    }
}
/// Incomplete Cholesky factorization result.
pub struct IncompleteCholesky {
    /// Lower triangular factor L such that A ≈ L * L^T.
    pub l: CsrMatrix,
}
impl IncompleteCholesky {
    /// Compute the IC(0) factorization of a symmetric positive-definite matrix.
    ///
    /// Only the lower-triangular part of `a` is used.
    pub fn factor(a: &CsrMatrix) -> Option<Self> {
        let n = a.nrows;
        assert_eq!(a.nrows, a.ncols, "IC(0) requires square matrix");
        let mut lval = vec![vec![0.0f64; n]; n];
        for r in 0..n {
            for k in a.row_ptr[r]..a.row_ptr[r + 1] {
                let c = a.col_idx[k];
                if c <= r {
                    lval[r][c] = a.values[k];
                }
            }
        }
        for j in 0..n {
            let mut sum_sq = 0.0f64;
            for k in 0..j {
                sum_sq += lval[j][k] * lval[j][k];
            }
            let d = lval[j][j] - sum_sq;
            if d <= 0.0 {
                return None;
            }
            lval[j][j] = d.sqrt();
            let ljj = lval[j][j];
            for i in (j + 1)..n {
                if lval[i][j].abs() < 1e-30 {
                    continue;
                }
                let mut dot = 0.0;
                for k in 0..j {
                    dot += lval[i][k] * lval[j][k];
                }
                lval[i][j] = (lval[i][j] - dot) / ljj;
            }
        }
        let mut rows = Vec::new();
        let mut cols = Vec::new();
        let mut vals = Vec::new();
        for i in 0..n {
            for j in 0..=i {
                if lval[i][j].abs() > 1e-30 {
                    rows.push(i);
                    cols.push(j);
                    vals.push(lval[i][j]);
                }
            }
        }
        let l = CsrMatrix::from_triplets(n, n, &rows, &cols, &vals);
        Some(IncompleteCholesky { l })
    }
    /// Apply L^{-1} (forward substitution).
    pub fn solve_lower(&self, b: &[f64]) -> Vec<f64> {
        let n = b.len();
        let mut x = vec![0.0f64; n];
        for i in 0..n {
            let mut s = b[i];
            for k in self.l.row_ptr[i]..self.l.row_ptr[i + 1] {
                let j = self.l.col_idx[k];
                if j < i {
                    s -= self.l.values[k] * x[j];
                }
            }
            let diag = self.l.get(i, i);
            x[i] = if diag.abs() < 1e-30 { 0.0 } else { s / diag };
        }
        x
    }
    /// Apply L^T^{-1} (backward substitution).
    pub fn solve_upper(&self, b: &[f64]) -> Vec<f64> {
        let n = b.len();
        let mut x = vec![0.0f64; n];
        for i in (0..n).rev() {
            let mut s = b[i];
            for k in self.l.row_ptr[i]..self.l.row_ptr[i + 1] {
                let j = self.l.col_idx[k];
                if j < i {
                    s -= self.l.values[k] * x[j];
                }
            }
            let diag = self.l.get(i, i);
            x[i] = if diag.abs() < 1e-30 { 0.0 } else { s / diag };
        }
        x
    }
    /// Apply the IC preconditioner: solve L L^T x = b.
    pub fn apply(&self, b: &[f64]) -> Vec<f64> {
        let y = self.solve_lower(b);
        let n = b.len();
        let mut x = vec![0.0f64; n];
        for i in (0..n).rev() {
            let mut s = y[i];
            for row in (i + 1)..n {
                let val = self.l.get(row, i);
                if val.abs() > 1e-30 {
                    s -= val * x[row];
                }
            }
            let diag = self.l.get(i, i);
            x[i] = if diag.abs() < 1e-30 { 0.0 } else { s / diag };
        }
        x
    }
}
/// Symmetric SOR (SSOR) preconditioner.
///
/// Applies one forward and one backward SOR sweep.
pub struct SsorPreconditioner {
    /// The matrix.
    pub a: CsrMatrix,
    /// Relaxation parameter ω (typically in (0, 2)).
    pub omega: f64,
}
impl SsorPreconditioner {
    /// Create a new SSOR preconditioner.
    pub fn new(a: CsrMatrix, omega: f64) -> Self {
        Self { a, omega }
    }
    /// Apply the SSOR preconditioner to a vector `b`.
    ///
    /// Computes M^{-1} b where M ≈ (D/ω + L) D^{-1} (D/ω + U) / ω.
    pub fn apply(&self, b: &[f64]) -> Vec<f64> {
        let n = b.len();
        let a = &self.a;
        let w = self.omega;
        let mut x = vec![0.0f64; n];
        for i in 0..n {
            let aii = a.get(i, i);
            if aii.abs() < f64::EPSILON {
                continue;
            }
            let mut sigma = 0.0;
            for k in a.row_ptr[i]..a.row_ptr[i + 1] {
                let j = a.col_idx[k];
                if j < i {
                    sigma += a.values[k] * x[j];
                }
            }
            x[i] = (b[i] - w * sigma) * w / aii;
        }
        for i in 0..n {
            let aii = a.get(i, i);
            x[i] *= aii / w;
        }
        for i in (0..n).rev() {
            let aii = a.get(i, i);
            if aii.abs() < f64::EPSILON {
                continue;
            }
            let mut sigma = 0.0;
            for k in a.row_ptr[i]..a.row_ptr[i + 1] {
                let j = a.col_idx[k];
                if j > i {
                    sigma += a.values[k] * x[j];
                }
            }
            x[i] = (x[i] - w * sigma) * w / aii;
        }
        x
    }
}
