#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! GPU-ready sparse matrix formats and iterative solvers with data-oriented layouts.
//!
//! Provides CSR, ELLPACK, Hybrid (ELL+COO), and Block-CSR formats along with
//! iterative solvers (CG, BiCGSTAB, preconditioned CG) suitable for GPU offload.
#![allow(missing_docs)]
#![allow(dead_code)]

// ---------------------------------------------------------------------------
// Sparse vector operations
// ---------------------------------------------------------------------------

/// Dot product of two vectors.
pub fn dot(x: &[f64], y: &[f64]) -> f64 {
    x.iter().zip(y.iter()).map(|(a, b)| a * b).sum()
}

/// `y + alpha * x` (AXPY).
pub fn axpy(alpha: f64, x: &[f64], y: &[f64]) -> Vec<f64> {
    y.iter()
        .zip(x.iter())
        .map(|(yi, xi)| yi + alpha * xi)
        .collect()
}

/// Euclidean norm of a vector.
pub fn norm2(x: &[f64]) -> f64 {
    dot(x, x).sqrt()
}

/// Scale every element of `x` by `s`.
pub fn scale_vec(x: &[f64], s: f64) -> Vec<f64> {
    x.iter().map(|v| v * s).collect()
}

// ---------------------------------------------------------------------------
// SparseTriplet – coordinate (COO) format for assembly
// ---------------------------------------------------------------------------

/// Coordinate (COO) format sparse matrix for incremental assembly.
pub struct SparseTriplet {
    pub rows: Vec<usize>,
    pub cols: Vec<usize>,
    pub vals: Vec<f64>,
}

impl SparseTriplet {
    /// Create an empty triplet store.
    pub fn new() -> Self {
        Self {
            rows: Vec::new(),
            cols: Vec::new(),
            vals: Vec::new(),
        }
    }

    /// Push a single entry `(row, col, val)`.
    pub fn add(&mut self, row: usize, col: usize, val: f64) {
        self.rows.push(row);
        self.cols.push(col);
        self.vals.push(val);
    }

    /// Convert to [`CsrMatrix`], sorting by (row, col) and summing duplicates.
    pub fn to_csr(&self, n_rows: usize, n_cols: usize) -> CsrMatrix {
        // Sort indices by (row, col)
        let mut order: Vec<usize> = (0..self.rows.len()).collect();
        order.sort_by_key(|&i| (self.rows[i], self.cols[i]));

        // Accumulate into (row, col, val) triples, summing duplicates
        let mut entries: Vec<(usize, usize, f64)> = Vec::new();
        for &i in &order {
            let r = self.rows[i];
            let c = self.cols[i];
            let v = self.vals[i];
            if let Some(last) = entries.last_mut()
                && last.0 == r
                && last.1 == c
            {
                last.2 += v;
                continue;
            }
            entries.push((r, c, v));
        }

        // Build CSR
        let nnz = entries.len();
        let mut row_ptr = vec![0usize; n_rows + 1];
        let mut col_idx = Vec::with_capacity(nnz);
        let mut values = Vec::with_capacity(nnz);

        for &(r, c, v) in &entries {
            row_ptr[r + 1] += 1;
            col_idx.push(c);
            values.push(v);
        }
        for i in 0..n_rows {
            row_ptr[i + 1] += row_ptr[i];
        }

        CsrMatrix {
            n_rows,
            n_cols,
            row_ptr,
            col_idx,
            values,
        }
    }
}

impl Default for SparseTriplet {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// CsrMatrix – Compressed Sparse Row
// ---------------------------------------------------------------------------

/// Compressed Sparse Row (CSR) matrix stored as plain `f64` arrays.
pub struct CsrMatrix {
    pub n_rows: usize,
    pub n_cols: usize,
    /// Row start indices, length `n_rows + 1`.
    pub row_ptr: Vec<usize>,
    /// Column indices of non-zeros.
    pub col_idx: Vec<usize>,
    /// Non-zero values.
    pub values: Vec<f64>,
}

impl CsrMatrix {
    /// Create an empty (all-zero) CSR matrix.
    pub fn new(n_rows: usize, n_cols: usize) -> Self {
        Self {
            n_rows,
            n_cols,
            row_ptr: vec![0; n_rows + 1],
            col_idx: Vec::new(),
            values: Vec::new(),
        }
    }

    /// Build from a dense row-major matrix (rows of length `n_cols`).
    pub fn from_dense(m: &[Vec<f64>]) -> Self {
        let n_rows = m.len();
        let n_cols = if n_rows > 0 { m[0].len() } else { 0 };
        let mut row_ptr = vec![0usize; n_rows + 1];
        let mut col_idx = Vec::new();
        let mut values = Vec::new();
        for (r, row) in m.iter().enumerate() {
            for (c, &v) in row.iter().enumerate() {
                if v != 0.0 {
                    col_idx.push(c);
                    values.push(v);
                }
            }
            row_ptr[r + 1] = col_idx.len();
        }
        Self {
            n_rows,
            n_cols,
            row_ptr,
            col_idx,
            values,
        }
    }

    /// Number of stored non-zeros.
    pub fn nnz(&self) -> usize {
        self.values.len()
    }

    /// Sparse matrix–vector product `y = A * x`.
    pub fn spmv(&self, x: &[f64]) -> Vec<f64> {
        let mut y = vec![0.0f64; self.n_rows];
        for r in 0..self.n_rows {
            let start = self.row_ptr[r];
            let end = self.row_ptr[r + 1];
            let mut sum = 0.0;
            for k in start..end {
                sum += self.values[k] * x[self.col_idx[k]];
            }
            y[r] = sum;
        }
        y
    }

    /// Get the value at `(row, col)`, returning `0.0` if not stored.
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

    /// Return `A^T` as a new `CsrMatrix`.
    pub fn transpose(&self) -> CsrMatrix {
        // Count non-zeros per column (which become rows of A^T)
        let mut row_ptr = vec![0usize; self.n_cols + 1];
        for &c in &self.col_idx {
            row_ptr[c + 1] += 1;
        }
        for i in 0..self.n_cols {
            row_ptr[i + 1] += row_ptr[i];
        }

        let nnz = self.values.len();
        let mut col_idx = vec![0usize; nnz];
        let mut values = vec![0.0f64; nnz];
        let mut pos = row_ptr[..self.n_cols].to_vec();

        for r in 0..self.n_rows {
            let start = self.row_ptr[r];
            let end = self.row_ptr[r + 1];
            for k in start..end {
                let c = self.col_idx[k];
                let dest = pos[c];
                col_idx[dest] = r;
                values[dest] = self.values[k];
                pos[c] += 1;
            }
        }

        CsrMatrix {
            n_rows: self.n_cols,
            n_cols: self.n_rows,
            row_ptr,
            col_idx,
            values,
        }
    }

    /// Add two CSR matrices of identical shape.
    pub fn add(&self, other: &CsrMatrix) -> CsrMatrix {
        assert_eq!(self.n_rows, other.n_rows);
        assert_eq!(self.n_cols, other.n_cols);
        // Assemble via triplet
        let mut trip = SparseTriplet::new();
        for r in 0..self.n_rows {
            for k in self.row_ptr[r]..self.row_ptr[r + 1] {
                trip.add(r, self.col_idx[k], self.values[k]);
            }
            for k in other.row_ptr[r]..other.row_ptr[r + 1] {
                trip.add(r, other.col_idx[k], other.values[k]);
            }
        }
        trip.to_csr(self.n_rows, self.n_cols)
    }

    /// Return a new matrix with every value multiplied by `s`.
    pub fn scale(&self, s: f64) -> CsrMatrix {
        CsrMatrix {
            n_rows: self.n_rows,
            n_cols: self.n_cols,
            row_ptr: self.row_ptr.clone(),
            col_idx: self.col_idx.clone(),
            values: self.values.iter().map(|v| v * s).collect(),
        }
    }
}

// ---------------------------------------------------------------------------
// EllMatrix – ELLPACK / ITPACK format (GPU-friendly padded layout)
// ---------------------------------------------------------------------------

/// ELLPACK-format sparse matrix: rows padded to `max_nnz_per_row`.
pub struct EllMatrix {
    pub n_rows: usize,
    pub n_cols: usize,
    pub max_nnz_per_row: usize,
    /// Column indices, row-major `[n_rows × max_nnz_per_row]`.
    pub col_idx: Vec<usize>,
    /// Values, row-major `[n_rows × max_nnz_per_row]`.
    pub values: Vec<f64>,
}

impl EllMatrix {
    /// Convert a [`CsrMatrix`] to ELLPACK.
    pub fn from_csr(csr: &CsrMatrix) -> Self {
        let n_rows = csr.n_rows;
        let n_cols = csr.n_cols;
        let max_nnz_per_row = (0..n_rows)
            .map(|r| csr.row_ptr[r + 1] - csr.row_ptr[r])
            .max()
            .unwrap_or(0);

        let size = n_rows * max_nnz_per_row;
        let mut col_idx = vec![0usize; size];
        let mut values = vec![0.0f64; size];

        for r in 0..n_rows {
            let start = csr.row_ptr[r];
            let end = csr.row_ptr[r + 1];
            for (j, k) in (start..end).enumerate() {
                col_idx[r * max_nnz_per_row + j] = csr.col_idx[k];
                values[r * max_nnz_per_row + j] = csr.values[k];
            }
        }

        Self {
            n_rows,
            n_cols,
            max_nnz_per_row,
            col_idx,
            values,
        }
    }

    /// Sparse matrix–vector product `y = A * x`.
    pub fn spmv(&self, x: &[f64]) -> Vec<f64> {
        let mut y = vec![0.0f64; self.n_rows];
        for r in 0..self.n_rows {
            let mut sum = 0.0;
            for j in 0..self.max_nnz_per_row {
                let v = self.values[r * self.max_nnz_per_row + j];
                if v != 0.0 {
                    let c = self.col_idx[r * self.max_nnz_per_row + j];
                    sum += v * x[c];
                }
            }
            y[r] = sum;
        }
        y
    }
}

// ---------------------------------------------------------------------------
// HybridMatrix – ELL + COO for irregular sparsity
// ---------------------------------------------------------------------------

/// Hybrid ELL+COO matrix: regular rows stored in ELL, overflow in COO.
pub struct HybridMatrix {
    pub ell: EllMatrix,
    pub coo_row: Vec<usize>,
    pub coo_col: Vec<usize>,
    pub coo_val: Vec<f64>,
}

impl HybridMatrix {
    /// Sparse matrix–vector product `y = (ELL + COO) * x`.
    pub fn spmv(&self, x: &[f64]) -> Vec<f64> {
        let mut y = self.ell.spmv(x);
        for k in 0..self.coo_val.len() {
            y[self.coo_row[k]] += self.coo_val[k] * x[self.coo_col[k]];
        }
        y
    }
}

// ---------------------------------------------------------------------------
// BlockCsrMatrix – block sparse row for FEM
// ---------------------------------------------------------------------------

/// Block-CSR matrix where every stored entry is a `block_size × block_size` dense tile.
pub struct BlockCsrMatrix {
    pub block_size: usize,
    pub n_block_rows: usize,
    pub n_block_cols: usize,
    /// Block row start indices, length `n_block_rows + 1`.
    pub row_ptr: Vec<usize>,
    /// Column (block) indices.
    pub col_idx: Vec<usize>,
    /// Dense tiles, each of length `block_size * block_size`.
    pub blocks: Vec<Vec<f64>>,
}

impl BlockCsrMatrix {
    /// Sparse matrix–vector product treating the matrix as `(n_block_rows * block_size)` × `(n_block_cols * block_size)`.
    pub fn spmv_block(&self, x: &[f64]) -> Vec<f64> {
        let bs = self.block_size;
        let n = self.n_block_rows * bs;
        let mut y = vec![0.0f64; n];
        for br in 0..self.n_block_rows {
            let row_start = self.row_ptr[br];
            let row_end = self.row_ptr[br + 1];
            for k in row_start..row_end {
                let bc = self.col_idx[k];
                let blk = &self.blocks[k];
                // Multiply dense block into y
                for i in 0..bs {
                    let mut s = 0.0;
                    for j in 0..bs {
                        s += blk[i * bs + j] * x[bc * bs + j];
                    }
                    y[br * bs + i] += s;
                }
            }
        }
        y
    }

    /// Convert to a flat [`CsrMatrix`].
    pub fn to_csr(&self) -> CsrMatrix {
        let bs = self.block_size;
        let n_rows = self.n_block_rows * bs;
        let n_cols = self.n_block_cols * bs;
        let mut trip = SparseTriplet::new();
        for br in 0..self.n_block_rows {
            for k in self.row_ptr[br]..self.row_ptr[br + 1] {
                let bc = self.col_idx[k];
                let blk = &self.blocks[k];
                for i in 0..bs {
                    for j in 0..bs {
                        let v = blk[i * bs + j];
                        if v != 0.0 {
                            trip.add(br * bs + i, bc * bs + j, v);
                        }
                    }
                }
            }
        }
        trip.to_csr(n_rows, n_cols)
    }
}

// ---------------------------------------------------------------------------
// Iterative solvers
// ---------------------------------------------------------------------------

/// Conjugate Gradient solver for symmetric positive-definite systems `A x = b`.
///
/// Returns `(solution, iterations_used)`.
pub fn cg_solve(
    a: &CsrMatrix,
    b: &[f64],
    x0: &[f64],
    max_iter: usize,
    tol: f64,
) -> (Vec<f64>, usize) {
    let n = b.len();
    let mut x = x0.to_vec();
    // r = b - A*x
    let ax = a.spmv(&x);
    let mut r: Vec<f64> = (0..n).map(|i| b[i] - ax[i]).collect();
    let mut p = r.clone();
    let mut rs_old = dot(&r, &r);

    for iter in 0..max_iter {
        if rs_old.sqrt() < tol {
            return (x, iter);
        }
        let ap = a.spmv(&p);
        let alpha = rs_old / dot(&p, &ap);
        x = axpy(alpha, &p, &x);
        r = axpy(-alpha, &ap, &r);
        let rs_new = dot(&r, &r);
        let beta = rs_new / rs_old;
        p = axpy(beta, &p, &r);
        rs_old = rs_new;
    }
    (x, max_iter)
}

/// BiCGSTAB solver for general (possibly non-symmetric) systems `A x = b`.
///
/// Returns `(solution, iterations_used)`.
pub fn bicgstab_solve(
    a: &CsrMatrix,
    b: &[f64],
    x0: &[f64],
    max_iter: usize,
    tol: f64,
) -> (Vec<f64>, usize) {
    let n = b.len();
    let mut x = x0.to_vec();
    let ax = a.spmv(&x);
    let mut r: Vec<f64> = (0..n).map(|i| b[i] - ax[i]).collect();
    let r_hat = r.clone();
    let mut rho = 1.0_f64;
    let mut alpha_s = 1.0_f64;
    let mut omega = 1.0_f64;
    let mut v = vec![0.0f64; n];
    let mut p = vec![0.0f64; n];

    for iter in 0..max_iter {
        if norm2(&r) < tol {
            return (x, iter);
        }
        let rho_new = dot(&r_hat, &r);
        let beta = (rho_new / rho) * (alpha_s / omega);
        p = axpy(beta, &p, &axpy(-beta * omega, &v, &r));
        v = a.spmv(&p);
        let denom = dot(&r_hat, &v);
        if denom.abs() < 1e-300 {
            return (x, iter);
        }
        alpha_s = rho_new / denom;
        let s: Vec<f64> = axpy(-alpha_s, &v, &r);
        if norm2(&s) < tol {
            x = axpy(alpha_s, &p, &x);
            return (x, iter + 1);
        }
        let t = a.spmv(&s);
        let tt = dot(&t, &t);
        omega = if tt.abs() < 1e-300 {
            0.0
        } else {
            dot(&t, &s) / tt
        };
        x = axpy(omega, &s, &axpy(alpha_s, &p, &x));
        r = axpy(-omega, &t, &s);
        rho = rho_new;
    }
    (x, max_iter)
}

/// Conjugate Gradient with Jacobi (diagonal) preconditioner for `A x = b` (SPD).
///
/// Returns `(solution, iterations_used)`.
pub fn jacobi_preconditioned_cg(
    a: &CsrMatrix,
    b: &[f64],
    max_iter: usize,
    tol: f64,
) -> (Vec<f64>, usize) {
    let n = b.len();
    // Build inverse diagonal preconditioner M^{-1}
    let mut m_inv = vec![1.0f64; n];
    for r in 0..n {
        let d = a.get(r, r);
        if d.abs() > 1e-300 {
            m_inv[r] = 1.0 / d;
        }
    }

    let mut x = vec![0.0f64; n];
    let ax = a.spmv(&x);
    let mut r: Vec<f64> = (0..n).map(|i| b[i] - ax[i]).collect();
    // z = M^{-1} r
    let z: Vec<f64> = (0..n).map(|i| m_inv[i] * r[i]).collect();
    let mut p = z.clone();
    let mut rz_old = dot(&r, &z);

    for iter in 0..max_iter {
        if norm2(&r) < tol {
            return (x, iter);
        }
        let ap = a.spmv(&p);
        let alpha = rz_old / dot(&p, &ap);
        x = axpy(alpha, &p, &x);
        r = axpy(-alpha, &ap, &r);
        let z_new: Vec<f64> = (0..n).map(|i| m_inv[i] * r[i]).collect();
        let rz_new = dot(&r, &z_new);
        let beta = rz_new / rz_old;
        p = axpy(beta, &p, &z_new);
        rz_old = rz_new;
    }
    (x, max_iter)
}

// ---------------------------------------------------------------------------
// GPU simulation utilities
// ---------------------------------------------------------------------------

/// Estimate SpMV throughput in GFLOPS given matrix dimensions and nnz count.
///
/// Uses a simple roofline model: bandwidth-bound at 100 GB/s with 12 bytes/flop.
pub fn simulate_spmv_throughput(n: usize, nnz: usize) -> f64 {
    // Memory traffic estimate: col_idx (8 bytes) + values (8 bytes) + x access (~random, 8 bytes)
    // = 24 bytes per nnz, plus row_ptr read (n * 8 bytes)
    let _ = n; // size used for context but dominated by nnz traffic
    let bytes_transferred = (nnz * 24) as f64;
    let bandwidth_gb_s = 100.0_f64; // typical GPU HBM bandwidth
    let time_s = bytes_transferred / (bandwidth_gb_s * 1e9);
    let flops = 2.0 * nnz as f64; // one multiply + one add per nnz
    flops / time_s / 1e9 // GFLOPS
}

/// Choose the ELLPACK row width (max non-zeros per row) to minimize padding waste.
///
/// Uses the 75th-percentile of the per-row nnz distribution.
pub fn optimal_ell_row_width(nnz_distribution: &[usize]) -> usize {
    if nnz_distribution.is_empty() {
        return 0;
    }
    let mut sorted = nnz_distribution.to_vec();
    sorted.sort_unstable();
    let idx = (sorted.len() * 3) / 4; // 75th percentile index
    sorted[idx]
}

// ---------------------------------------------------------------------------
// SpMV – segmented (row-parallel) variant
// ---------------------------------------------------------------------------

/// Segmented SpMV: processes each row independently to prepare for
/// GPU-style parallel execution. Functionally identical to `CsrMatrix::spmv`
/// but structured for row-parallel dispatch.
pub fn spmv_segmented(a: &CsrMatrix, x: &[f64]) -> Vec<f64> {
    let mut y = vec![0.0_f64; a.n_rows];
    for r in 0..a.n_rows {
        let start = a.row_ptr[r];
        let end = a.row_ptr[r + 1];
        let mut acc = 0.0_f64;
        for k in start..end {
            acc += a.values[k] * x[a.col_idx[k]];
        }
        y[r] = acc;
    }
    y
}

// ---------------------------------------------------------------------------
// Sparse matrix assembly helpers
// ---------------------------------------------------------------------------

/// Assemble a 1D Laplacian matrix of size `n × n` (tridiagonal: 2 on diag, -1 off-diag).
pub fn assemble_1d_laplacian(n: usize) -> CsrMatrix {
    let mut trip = SparseTriplet::new();
    for i in 0..n {
        trip.add(i, i, 2.0);
        if i > 0 {
            trip.add(i, i - 1, -1.0);
        }
        if i + 1 < n {
            trip.add(i, i + 1, -1.0);
        }
    }
    trip.to_csr(n, n)
}

// ---------------------------------------------------------------------------
// CSR-to-ELL conversion (alternative entry point)
// ---------------------------------------------------------------------------

/// Convert a CSR matrix to ELLPACK format (convenience wrapper).
pub fn csr_to_ell(csr: &CsrMatrix) -> EllMatrix {
    EllMatrix::from_csr(csr)
}

// ---------------------------------------------------------------------------
// Extract diagonal
// ---------------------------------------------------------------------------

/// Extract the main diagonal of a CSR matrix.
pub fn extract_diagonal(a: &CsrMatrix) -> Vec<f64> {
    let n = a.n_rows.min(a.n_cols);
    let mut diag = vec![0.0_f64; n];
    for r in 0..n {
        diag[r] = a.get(r, r);
    }
    diag
}

// ---------------------------------------------------------------------------
// Per-row nnz distribution
// ---------------------------------------------------------------------------

/// Compute the number of non-zeros per row for a CSR matrix.
pub fn compute_nnz_per_row(a: &CsrMatrix) -> Vec<usize> {
    (0..a.n_rows)
        .map(|r| a.row_ptr[r + 1] - a.row_ptr[r])
        .collect()
}

// ---------------------------------------------------------------------------
// Frobenius norm
// ---------------------------------------------------------------------------

/// Compute the Frobenius norm of a sparse matrix: `sqrt(sum(a_ij^2))`.
pub fn frobenius_norm(a: &CsrMatrix) -> f64 {
    let sum_sq: f64 = a.values.iter().map(|v| v * v).sum();
    sum_sq.sqrt()
}

// ---------------------------------------------------------------------------
// Sparse triangular solves
// ---------------------------------------------------------------------------

/// Forward-substitution solve `L x = b` where `L` is lower-triangular (CSR).
///
/// Assumes `L` has non-zero diagonal entries. The diagonal entry for row `i`
/// is taken as `L.get(i, i)`.
pub fn sparse_lower_triangular_solve(l: &CsrMatrix, b: &[f64]) -> Vec<f64> {
    let n = b.len();
    let mut x = vec![0.0_f64; n];
    for i in 0..n {
        let mut sum = b[i];
        let start = l.row_ptr[i];
        let end = l.row_ptr[i + 1];
        let mut diag = 1.0_f64;
        for k in start..end {
            let c = l.col_idx[k];
            if c < i {
                sum -= l.values[k] * x[c];
            } else if c == i {
                diag = l.values[k];
            }
        }
        x[i] = sum / diag;
    }
    x
}

/// Back-substitution solve `U x = b` where `U` is upper-triangular (CSR).
///
/// Assumes `U` has non-zero diagonal entries.
pub fn sparse_upper_triangular_solve(u: &CsrMatrix, b: &[f64]) -> Vec<f64> {
    let n = b.len();
    let mut x = vec![0.0_f64; n];
    for i in (0..n).rev() {
        let mut sum = b[i];
        let start = u.row_ptr[i];
        let end = u.row_ptr[i + 1];
        let mut diag = 1.0_f64;
        for k in start..end {
            let c = u.col_idx[k];
            if c > i {
                sum -= u.values[k] * x[c];
            } else if c == i {
                diag = u.values[k];
            }
        }
        x[i] = sum / diag;
    }
    x
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_csr_from_dense_nnz() {
        let m = vec![
            vec![1.0, 0.0, 2.0],
            vec![0.0, 3.0, 0.0],
            vec![4.0, 5.0, 6.0],
        ];
        let csr = CsrMatrix::from_dense(&m);
        assert_eq!(csr.n_rows, 3);
        assert_eq!(csr.n_cols, 3);
        assert_eq!(csr.nnz(), 6);
    }

    #[test]
    fn test_csr_spmv_identity() {
        // 3x3 identity
        let m = vec![
            vec![1.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
            vec![0.0, 0.0, 1.0],
        ];
        let csr = CsrMatrix::from_dense(&m);
        let x = vec![3.0, 7.0, -2.0];
        let y = csr.spmv(&x);
        assert_eq!(y, x);
    }

    #[test]
    fn test_csr_spmv_known_3x3() {
        // A = [[2,1,0],[1,3,1],[0,1,2]]
        let m = vec![
            vec![2.0, 1.0, 0.0],
            vec![1.0, 3.0, 1.0],
            vec![0.0, 1.0, 2.0],
        ];
        let csr = CsrMatrix::from_dense(&m);
        let x = vec![1.0, 2.0, 3.0];
        let y = csr.spmv(&x);
        // y[0] = 2+2 = 4, y[1] = 1+6+3 = 10, y[2] = 2+6 = 8
        assert!((y[0] - 4.0).abs() < 1e-12);
        assert!((y[1] - 10.0).abs() < 1e-12);
        assert!((y[2] - 8.0).abs() < 1e-12);
    }

    #[test]
    fn test_cg_solve_diagonal_spd() {
        // A = diag(1,2,3,4), b = [1,2,3,4], exact solution x = [1,1,1,1]
        let m = vec![
            vec![1.0, 0.0, 0.0, 0.0],
            vec![0.0, 2.0, 0.0, 0.0],
            vec![0.0, 0.0, 3.0, 0.0],
            vec![0.0, 0.0, 0.0, 4.0],
        ];
        let a = CsrMatrix::from_dense(&m);
        let b = vec![1.0, 2.0, 3.0, 4.0];
        let x0 = vec![0.0; 4];
        let (x, _iters) = cg_solve(&a, &b, &x0, 100, 1e-12);
        for v in &x {
            assert!((v - 1.0).abs() < 1e-10, "x value {v} not close to 1.0");
        }
    }

    #[test]
    fn test_sparse_triplet_to_csr_duplicate_sum() {
        let mut trip = SparseTriplet::new();
        trip.add(0, 0, 1.0);
        trip.add(0, 0, 2.0); // duplicate → should sum to 3.0
        trip.add(1, 1, 5.0);
        let csr = trip.to_csr(2, 2);
        assert!((csr.get(0, 0) - 3.0).abs() < 1e-12);
        assert!((csr.get(1, 1) - 5.0).abs() < 1e-12);
        assert_eq!(csr.nnz(), 2);
    }

    #[test]
    fn test_ell_spmv_matches_csr() {
        let m = vec![
            vec![2.0, 1.0, 0.0],
            vec![1.0, 3.0, 1.0],
            vec![0.0, 1.0, 2.0],
        ];
        let csr = CsrMatrix::from_dense(&m);
        let ell = EllMatrix::from_csr(&csr);
        let x = vec![1.0, -1.0, 2.0];
        let y_csr = csr.spmv(&x);
        let y_ell = ell.spmv(&x);
        for (a, b) in y_csr.iter().zip(y_ell.iter()) {
            assert!((a - b).abs() < 1e-12, "ELL mismatch: {a} vs {b}");
        }
    }

    #[test]
    fn test_csr_transpose() {
        // A = [[1,2,3],[4,5,6]]  →  A^T = [[1,4],[2,5],[3,6]]
        let m = vec![vec![1.0, 2.0, 3.0], vec![4.0, 5.0, 6.0]];
        let csr = CsrMatrix::from_dense(&m);
        let at = csr.transpose();
        assert_eq!(at.n_rows, 3);
        assert_eq!(at.n_cols, 2);
        assert!((at.get(0, 0) - 1.0).abs() < 1e-12);
        assert!((at.get(0, 1) - 4.0).abs() < 1e-12);
        assert!((at.get(1, 0) - 2.0).abs() < 1e-12);
        assert!((at.get(1, 1) - 5.0).abs() < 1e-12);
        assert!((at.get(2, 0) - 3.0).abs() < 1e-12);
        assert!((at.get(2, 1) - 6.0).abs() < 1e-12);
    }

    // ── spmv_segmented ─────────────────────────────────────────────────────

    #[test]
    fn test_spmv_segmented_identity() {
        let m = vec![
            vec![1.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
            vec![0.0, 0.0, 1.0],
        ];
        let csr = CsrMatrix::from_dense(&m);
        let x = vec![3.0, 7.0, -2.0];
        let y = spmv_segmented(&csr, &x);
        for (a, b) in y.iter().zip(x.iter()) {
            assert!((a - b).abs() < 1e-12);
        }
    }

    #[test]
    fn test_spmv_segmented_matches_csr() {
        let m = vec![
            vec![2.0, 1.0, 0.0],
            vec![1.0, 3.0, 1.0],
            vec![0.0, 1.0, 2.0],
        ];
        let csr = CsrMatrix::from_dense(&m);
        let x = vec![1.0, -1.0, 2.0];
        let y_std = csr.spmv(&x);
        let y_seg = spmv_segmented(&csr, &x);
        for (a, b) in y_std.iter().zip(y_seg.iter()) {
            assert!((a - b).abs() < 1e-12);
        }
    }

    // ── assemble_1d_laplacian ──────────────────────────────────────────────

    #[test]
    fn test_assemble_1d_laplacian_3x3() {
        let l = assemble_1d_laplacian(3);
        // Expected: [[2,-1,0],[-1,2,-1],[0,-1,2]]
        assert!((l.get(0, 0) - 2.0).abs() < 1e-12);
        assert!((l.get(0, 1) - (-1.0)).abs() < 1e-12);
        assert!((l.get(0, 2)).abs() < 1e-12);
        assert!((l.get(1, 0) - (-1.0)).abs() < 1e-12);
        assert!((l.get(1, 1) - 2.0).abs() < 1e-12);
        assert!((l.get(1, 2) - (-1.0)).abs() < 1e-12);
        assert!((l.get(2, 2) - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_assemble_1d_laplacian_spd() {
        // 1D Laplacian is symmetric positive definite
        let n = 5;
        let l = assemble_1d_laplacian(n);
        // Check symmetry
        for i in 0..n {
            for j in 0..n {
                assert!((l.get(i, j) - l.get(j, i)).abs() < 1e-12);
            }
        }
        // CG should converge for SPD systems
        let b = vec![1.0; n];
        let x0 = vec![0.0; n];
        let (x, iters) = cg_solve(&l, &b, &x0, 200, 1e-10);
        assert!(iters < 200);
        // Verify Ax ≈ b
        let ax = l.spmv(&x);
        for i in 0..n {
            assert!((ax[i] - b[i]).abs() < 1e-8);
        }
    }

    // ── csr_to_ell ─────────────────────────────────────────────────────────

    #[test]
    fn test_csr_to_ell_spmv() {
        let m = vec![
            vec![5.0, 0.0, 1.0, 0.0],
            vec![0.0, 3.0, 0.0, 2.0],
            vec![1.0, 0.0, 4.0, 0.0],
            vec![0.0, 0.0, 0.0, 6.0],
        ];
        let csr = CsrMatrix::from_dense(&m);
        let ell = csr_to_ell(&csr);
        let x = vec![1.0, 2.0, 3.0, 4.0];
        let y_csr = csr.spmv(&x);
        let y_ell = ell.spmv(&x);
        for (a, b) in y_csr.iter().zip(y_ell.iter()) {
            assert!((a - b).abs() < 1e-12, "mismatch: {a} vs {b}");
        }
    }

    #[test]
    fn test_csr_to_ell_max_nnz() {
        let m = vec![
            vec![1.0, 0.0, 0.0],
            vec![1.0, 2.0, 3.0], // 3 nnz → max
            vec![0.0, 0.0, 1.0],
        ];
        let csr = CsrMatrix::from_dense(&m);
        let ell = csr_to_ell(&csr);
        assert_eq!(ell.max_nnz_per_row, 3);
    }

    // ── BlockCsrMatrix additional tests ────────────────────────────────────

    #[test]
    fn test_block_csr_spmv_2x2() {
        // Single 2x2 block [[1,2],[3,4]] at block position (0,0)
        let bcsr = BlockCsrMatrix {
            block_size: 2,
            n_block_rows: 1,
            n_block_cols: 1,
            row_ptr: vec![0, 1],
            col_idx: vec![0],
            blocks: vec![vec![1.0, 2.0, 3.0, 4.0]],
        };
        let x = vec![1.0, 1.0];
        let y = bcsr.spmv_block(&x);
        assert!((y[0] - 3.0).abs() < 1e-12); // 1+2
        assert!((y[1] - 7.0).abs() < 1e-12); // 3+4
    }

    #[test]
    fn test_block_csr_to_csr_roundtrip() {
        let bcsr = BlockCsrMatrix {
            block_size: 2,
            n_block_rows: 2,
            n_block_cols: 2,
            row_ptr: vec![0, 1, 2],
            col_idx: vec![0, 1],
            blocks: vec![vec![1.0, 2.0, 3.0, 4.0], vec![5.0, 6.0, 7.0, 8.0]],
        };
        let csr = bcsr.to_csr();
        let x = vec![1.0, 1.0, 1.0, 1.0];
        let y_block = bcsr.spmv_block(&x);
        let y_csr = csr.spmv(&x);
        for (a, b) in y_block.iter().zip(y_csr.iter()) {
            assert!((a - b).abs() < 1e-12);
        }
    }

    // ── sparse_lower_triangular_solve ──────────────────────────────────────

    #[test]
    fn test_lower_tri_solve_identity() {
        let m = vec![
            vec![1.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
            vec![0.0, 0.0, 1.0],
        ];
        let l = CsrMatrix::from_dense(&m);
        let b = vec![3.0, 7.0, -2.0];
        let x = sparse_lower_triangular_solve(&l, &b);
        for (a, bv) in x.iter().zip(b.iter()) {
            assert!((a - bv).abs() < 1e-12);
        }
    }

    #[test]
    fn test_lower_tri_solve_3x3() {
        // L = [[2,0,0],[1,3,0],[4,2,5]]
        let m = vec![
            vec![2.0, 0.0, 0.0],
            vec![1.0, 3.0, 0.0],
            vec![4.0, 2.0, 5.0],
        ];
        let l = CsrMatrix::from_dense(&m);
        let b = vec![4.0, 7.0, 26.0];
        let x = sparse_lower_triangular_solve(&l, &b);
        // x[0] = 4/2 = 2
        // x[1] = (7 - 1*2)/3 = 5/3
        // x[2] = (26 - 4*2 - 2*(5/3))/5
        assert!((x[0] - 2.0).abs() < 1e-10);
        assert!((x[1] - 5.0 / 3.0).abs() < 1e-10);
        let expected_x2 = (26.0 - 8.0 - 10.0 / 3.0) / 5.0;
        assert!((x[2] - expected_x2).abs() < 1e-10);
    }

    #[test]
    fn test_lower_tri_solve_verify_lx_eq_b() {
        let m = vec![
            vec![3.0, 0.0, 0.0, 0.0],
            vec![1.0, 2.0, 0.0, 0.0],
            vec![0.0, 4.0, 5.0, 0.0],
            vec![2.0, 0.0, 1.0, 6.0],
        ];
        let l = CsrMatrix::from_dense(&m);
        let b = vec![9.0, 8.0, 22.0, 29.0];
        let x = sparse_lower_triangular_solve(&l, &b);
        // Verify L*x = b
        let lx = l.spmv(&x);
        for i in 0..4 {
            assert!(
                (lx[i] - b[i]).abs() < 1e-10,
                "row {i}: {} vs {}",
                lx[i],
                b[i]
            );
        }
    }

    // ── sparse_upper_triangular_solve ──────────────────────────────────────

    #[test]
    fn test_upper_tri_solve_identity() {
        let m = vec![
            vec![1.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
            vec![0.0, 0.0, 1.0],
        ];
        let u = CsrMatrix::from_dense(&m);
        let b = vec![3.0, 7.0, -2.0];
        let x = sparse_upper_triangular_solve(&u, &b);
        for (a, bv) in x.iter().zip(b.iter()) {
            assert!((a - bv).abs() < 1e-12);
        }
    }

    #[test]
    fn test_upper_tri_solve_verify_ux_eq_b() {
        // U = [[2,1,3],[0,4,2],[0,0,5]]
        let m = vec![
            vec![2.0, 1.0, 3.0],
            vec![0.0, 4.0, 2.0],
            vec![0.0, 0.0, 5.0],
        ];
        let u = CsrMatrix::from_dense(&m);
        let b = vec![13.0, 14.0, 10.0];
        let x = sparse_upper_triangular_solve(&u, &b);
        let ux = u.spmv(&x);
        for i in 0..3 {
            assert!(
                (ux[i] - b[i]).abs() < 1e-10,
                "row {i}: {} vs {}",
                ux[i],
                b[i]
            );
        }
    }

    // ── CsrMatrix::add & scale ─────────────────────────────────────────────

    #[test]
    fn test_csr_add() {
        let m1 = vec![vec![1.0, 0.0], vec![0.0, 2.0]];
        let m2 = vec![vec![0.0, 3.0], vec![4.0, 0.0]];
        let a = CsrMatrix::from_dense(&m1);
        let b = CsrMatrix::from_dense(&m2);
        let c = a.add(&b);
        assert!((c.get(0, 0) - 1.0).abs() < 1e-12);
        assert!((c.get(0, 1) - 3.0).abs() < 1e-12);
        assert!((c.get(1, 0) - 4.0).abs() < 1e-12);
        assert!((c.get(1, 1) - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_csr_scale() {
        let m = vec![vec![1.0, 2.0], vec![3.0, 4.0]];
        let a = CsrMatrix::from_dense(&m);
        let b = a.scale(2.0);
        assert!((b.get(0, 0) - 2.0).abs() < 1e-12);
        assert!((b.get(0, 1) - 4.0).abs() < 1e-12);
        assert!((b.get(1, 0) - 6.0).abs() < 1e-12);
        assert!((b.get(1, 1) - 8.0).abs() < 1e-12);
    }

    // ── BiCGSTAB ───────────────────────────────────────────────────────────

    #[test]
    fn test_bicgstab_diagonal() {
        let m = vec![
            vec![2.0, 0.0, 0.0],
            vec![0.0, 3.0, 0.0],
            vec![0.0, 0.0, 4.0],
        ];
        let a = CsrMatrix::from_dense(&m);
        let b = vec![4.0, 9.0, 16.0];
        let x0 = vec![0.0; 3];
        let (x, _iters) = bicgstab_solve(&a, &b, &x0, 100, 1e-10);
        assert!((x[0] - 2.0).abs() < 1e-8);
        assert!((x[1] - 3.0).abs() < 1e-8);
        assert!((x[2] - 4.0).abs() < 1e-8);
    }

    #[test]
    fn test_bicgstab_nonsymmetric() {
        // Non-symmetric but diagonally dominant
        let m = vec![
            vec![4.0, 1.0, 0.0],
            vec![0.0, 3.0, 1.0],
            vec![0.0, 0.0, 5.0],
        ];
        let a = CsrMatrix::from_dense(&m);
        let b = vec![5.0, 4.0, 5.0];
        let x0 = vec![0.0; 3];
        let (x, _iters) = bicgstab_solve(&a, &b, &x0, 200, 1e-10);
        // Verify Ax ≈ b
        let ax = a.spmv(&x);
        for i in 0..3 {
            assert!(
                (ax[i] - b[i]).abs() < 1e-6,
                "row {i}: {} vs {}",
                ax[i],
                b[i]
            );
        }
    }

    // ── Jacobi preconditioned CG ───────────────────────────────────────────

    #[test]
    fn test_jacobi_pcg_laplacian() {
        let n = 8;
        let a = assemble_1d_laplacian(n);
        let b: Vec<f64> = (0..n).map(|i| (i as f64 + 1.0).sin()).collect();
        let (x, iters) = jacobi_preconditioned_cg(&a, &b, 500, 1e-10);
        assert!(iters < 500, "PCG should converge, used {iters} iterations");
        let ax = a.spmv(&x);
        for i in 0..n {
            assert!(
                (ax[i] - b[i]).abs() < 1e-6,
                "row {i}: {} vs {}",
                ax[i],
                b[i]
            );
        }
    }

    // ── GPU simulation utilities ───────────────────────────────────────────

    #[test]
    fn test_simulate_spmv_throughput() {
        let gflops = simulate_spmv_throughput(100, 1000);
        assert!(gflops > 0.0);
        // With 100 GB/s bandwidth and 24 bytes/nnz, bigger nnz should give higher GFLOPS
        let gflops2 = simulate_spmv_throughput(1000, 10000);
        assert!((gflops2 - gflops).abs() < 1.0); // roofline is flat w.r.t. nnz ratio
    }

    #[test]
    fn test_optimal_ell_row_width_empty() {
        assert_eq!(optimal_ell_row_width(&[]), 0);
    }

    #[test]
    fn test_optimal_ell_row_width_uniform() {
        // All rows have 5 nnz
        let dist = vec![5; 10];
        assert_eq!(optimal_ell_row_width(&dist), 5);
    }

    #[test]
    fn test_optimal_ell_row_width_varied() {
        let dist = vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
        let w = optimal_ell_row_width(&dist);
        // 75th percentile index = (10 * 3) / 4 = 7 → sorted[7] = 8
        assert_eq!(w, 8);
    }

    // ── SparseTriplet default ──────────────────────────────────────────────

    #[test]
    fn test_sparse_triplet_default() {
        let t = SparseTriplet::default();
        assert!(t.rows.is_empty());
        assert!(t.cols.is_empty());
        assert!(t.vals.is_empty());
    }

    // ── vector operations ──────────────────────────────────────────────────

    #[test]
    fn test_dot_product() {
        let x = vec![1.0, 2.0, 3.0];
        let y = vec![4.0, 5.0, 6.0];
        assert!((dot(&x, &y) - 32.0).abs() < 1e-12);
    }

    #[test]
    fn test_axpy() {
        let x = vec![1.0, 2.0, 3.0];
        let y = vec![10.0, 20.0, 30.0];
        let z = axpy(2.0, &x, &y);
        assert_eq!(z, vec![12.0, 24.0, 36.0]);
    }

    #[test]
    fn test_norm2() {
        let x = vec![3.0, 4.0];
        assert!((norm2(&x) - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_scale_vec() {
        let x = vec![1.0, 2.0, 3.0];
        let s = scale_vec(&x, 3.0);
        assert_eq!(s, vec![3.0, 6.0, 9.0]);
    }

    // ── HybridMatrix ──────────────────────────────────────────────────────

    #[test]
    fn test_hybrid_spmv() {
        let m = vec![
            vec![1.0, 2.0, 0.0],
            vec![0.0, 3.0, 0.0],
            vec![0.0, 0.0, 4.0],
        ];
        let csr = CsrMatrix::from_dense(&m);
        let ell = EllMatrix::from_csr(&csr);
        let hybrid = HybridMatrix {
            ell,
            coo_row: vec![],
            coo_col: vec![],
            coo_val: vec![],
        };
        let x = vec![1.0, 1.0, 1.0];
        let y = hybrid.spmv(&x);
        assert!((y[0] - 3.0).abs() < 1e-12);
        assert!((y[1] - 3.0).abs() < 1e-12);
        assert!((y[2] - 4.0).abs() < 1e-12);
    }

    #[test]
    fn test_hybrid_spmv_with_coo() {
        // ELL part: identity, COO part: adds off-diagonal
        let m = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let csr = CsrMatrix::from_dense(&m);
        let ell = EllMatrix::from_csr(&csr);
        let hybrid = HybridMatrix {
            ell,
            coo_row: vec![0],
            coo_col: vec![1],
            coo_val: vec![5.0],
        };
        let x = vec![1.0, 2.0];
        let y = hybrid.spmv(&x);
        assert!((y[0] - 11.0).abs() < 1e-12); // 1*1 + 5*2
        assert!((y[1] - 2.0).abs() < 1e-12);
    }

    // ── CsrMatrix from empty ──────────────────────────────────────────────

    #[test]
    fn test_csr_empty() {
        let csr = CsrMatrix::new(3, 3);
        assert_eq!(csr.nnz(), 0);
        let y = csr.spmv(&[1.0, 2.0, 3.0]);
        assert_eq!(y, vec![0.0, 0.0, 0.0]);
    }

    // ── extract_diagonal ──────────────────────────────────────────────────

    #[test]
    fn test_extract_diagonal() {
        let m = vec![
            vec![5.0, 1.0, 0.0],
            vec![0.0, 3.0, 2.0],
            vec![0.0, 0.0, 7.0],
        ];
        let csr = CsrMatrix::from_dense(&m);
        let diag = extract_diagonal(&csr);
        assert!((diag[0] - 5.0).abs() < 1e-12);
        assert!((diag[1] - 3.0).abs() < 1e-12);
        assert!((diag[2] - 7.0).abs() < 1e-12);
    }

    #[test]
    fn test_extract_diagonal_missing() {
        // Matrix with zero diagonal
        let m = vec![vec![0.0, 1.0], vec![1.0, 0.0]];
        let csr = CsrMatrix::from_dense(&m);
        let diag = extract_diagonal(&csr);
        assert!((diag[0]).abs() < 1e-12);
        assert!((diag[1]).abs() < 1e-12);
    }

    // ── compute_nnz_per_row ────────────────────────────────────────────────

    #[test]
    fn test_compute_nnz_per_row() {
        let m = vec![
            vec![1.0, 2.0, 0.0],
            vec![0.0, 3.0, 0.0],
            vec![4.0, 5.0, 6.0],
        ];
        let csr = CsrMatrix::from_dense(&m);
        let nnz = compute_nnz_per_row(&csr);
        assert_eq!(nnz, vec![2, 1, 3]);
    }

    // ── frobenius_norm ─────────────────────────────────────────────────────

    #[test]
    fn test_frobenius_norm() {
        // Identity 3x3 → frobenius = sqrt(3)
        let m = vec![
            vec![1.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
            vec![0.0, 0.0, 1.0],
        ];
        let csr = CsrMatrix::from_dense(&m);
        let f = frobenius_norm(&csr);
        assert!((f - 3.0_f64.sqrt()).abs() < 1e-12);
    }

    #[test]
    fn test_frobenius_norm_known() {
        let m = vec![vec![3.0, 4.0], vec![0.0, 0.0]];
        let csr = CsrMatrix::from_dense(&m);
        let f = frobenius_norm(&csr);
        // sqrt(9 + 16) = 5
        assert!((f - 5.0).abs() < 1e-12);
    }
}
