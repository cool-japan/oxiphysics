#![allow(clippy::needless_range_loop, clippy::ptr_arg)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! GPU sparse linear system solver (CPU mock implementation).
//!
//! Provides CSR sparse matrix storage, SpMV, dot product, AXPY, conjugate
//! gradient (CG), and preconditioned CG solvers that mirror GPU dispatch
//! semantics while executing on the CPU.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

// ── CSR Sparse Matrix ────────────────────────────────────────────────────────

/// Compressed Sparse Row (CSR) matrix stored in GPU-friendly layout.
///
/// Rows have entries in `col_idx[row_ptr[i]..row_ptr[i+1]]` with
/// corresponding values in `values[row_ptr[i]..row_ptr[i+1]]`.
#[derive(Debug, Clone)]
pub struct SparseMatrixGpu {
    /// Number of rows (and columns — assumed square).
    pub n: usize,
    /// Number of non-zero entries.
    pub nnz: usize,
    /// Row pointer array (length n+1).
    pub row_ptr: Vec<usize>,
    /// Column index array (length nnz).
    pub col_idx: Vec<usize>,
    /// Value array (length nnz).
    pub values: Vec<f64>,
}

impl SparseMatrixGpu {
    /// Create a new empty n×n sparse matrix.
    pub fn new(n: usize) -> Self {
        Self {
            n,
            nnz: 0,
            row_ptr: vec![0usize; n + 1],
            col_idx: Vec::new(),
            values: Vec::new(),
        }
    }

    /// Construct from pre-built CSR arrays.
    ///
    /// # Panics
    /// Panics in debug mode if `row_ptr.len() != n + 1`.
    pub fn from_csr(n: usize, row_ptr: Vec<usize>, col_idx: Vec<usize>, values: Vec<f64>) -> Self {
        debug_assert_eq!(row_ptr.len(), n + 1);
        let nnz = values.len();
        Self {
            n,
            nnz,
            row_ptr,
            col_idx,
            values,
        }
    }

    /// Append a single entry `(row, col, val)` in COO style.
    ///
    /// This invalidates the CSR structure until [`finalize`](Self::finalize) is called.
    pub fn add_entry(&mut self, row: usize, col: usize, val: f64) {
        // Store as a pending triplet encoded in the arrays.
        self.col_idx.push(col);
        self.values.push(val);
        // Use row_ptr[row+1] as a count (finalize will prefix-sum them).
        if row + 1 < self.row_ptr.len() {
            self.row_ptr[row + 1] += 1;
        }
        self.nnz += 1;
    }

    /// Convert count-per-row in `row_ptr` to proper CSR offsets.
    ///
    /// Must be called after a series of [`add_entry`](Self::add_entry) calls
    /// and before any SpMV or solver operations.
    pub fn finalize(&mut self) {
        // Prefix-sum row_ptr in-place.
        for i in 1..=self.n {
            self.row_ptr[i] += self.row_ptr[i - 1];
        }
        self.nnz = *self.row_ptr.last().unwrap_or(&0);
    }

    /// Return the number of non-zero entries.
    pub fn nnz(&self) -> usize {
        self.nnz
    }

    /// Return the number of rows.
    pub fn rows(&self) -> usize {
        self.n
    }

    /// Extract the main diagonal as a dense vector.
    ///
    /// Missing diagonal entries are filled with zero.
    pub fn diagonal(&self) -> Vec<f64> {
        let mut diag = vec![0.0f64; self.n];
        for row in 0..self.n {
            let start = self.row_ptr[row];
            let end = self.row_ptr[row + 1];
            for k in start..end {
                if self.col_idx[k] == row {
                    diag[row] = self.values[k];
                }
            }
        }
        diag
    }

    /// Check whether the matrix is structurally symmetric (A\[i,j\] ↔ A\[j,i\]).
    pub fn is_symmetric(&self) -> bool {
        for row in 0..self.n {
            let start = self.row_ptr[row];
            let end = self.row_ptr[row + 1];
            for k in start..end {
                let col = self.col_idx[k];
                let val = self.values[k];
                // Look for the transpose entry.
                let mut found = false;
                let cs = self.row_ptr[col];
                let ce = self.row_ptr[col + 1];
                for m in cs..ce {
                    if self.col_idx[m] == row {
                        if (self.values[m] - val).abs() > 1e-12 {
                            return false;
                        }
                        found = true;
                        break;
                    }
                }
                if !found {
                    return false;
                }
            }
        }
        true
    }

    /// Frobenius norm: √(∑ aᵢⱼ²).
    pub fn frobenius_norm(&self) -> f64 {
        self.values.iter().map(|v| v * v).sum::<f64>().sqrt()
    }
}

// ── Utility constructors ─────────────────────────────────────────────────────

/// Build an n×n identity matrix in CSR format.
pub fn sparse_identity(n: usize) -> SparseMatrixGpu {
    let row_ptr: Vec<usize> = (0..=n).collect();
    let col_idx: Vec<usize> = (0..n).collect();
    let values = vec![1.0f64; n];
    SparseMatrixGpu::from_csr(n, row_ptr, col_idx, values)
}

/// Build a diagonal matrix from a dense vector.
pub fn sparse_diagonal_matrix(diag: &[f64]) -> SparseMatrixGpu {
    let n = diag.len();
    let row_ptr: Vec<usize> = (0..=n).collect();
    let col_idx: Vec<usize> = (0..n).collect();
    let values = diag.to_vec();
    SparseMatrixGpu::from_csr(n, row_ptr, col_idx, values)
}

// ── BLAS-like primitives ─────────────────────────────────────────────────────

/// Sparse matrix–vector product: y = A · x.
///
/// # Panics
/// Panics if `x.len() != mat.n`.
pub fn gpu_spmv(mat: &SparseMatrixGpu, x: &[f64]) -> Vec<f64> {
    assert_eq!(x.len(), mat.n, "gpu_spmv: x length mismatch");
    let mut y = vec![0.0f64; mat.n];
    for row in 0..mat.n {
        let start = mat.row_ptr[row];
        let end = mat.row_ptr[row + 1];
        let mut sum = 0.0f64;
        for k in start..end {
            sum += mat.values[k] * x[mat.col_idx[k]];
        }
        y[row] = sum;
    }
    y
}

/// Dense dot product: ∑ aᵢ bᵢ.
///
/// Returns 0 for empty slices.
pub fn gpu_dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(ai, bi)| ai * bi).sum()
}

/// In-place AXPY: y ← a·x + y.
///
/// # Panics
/// Panics if `x.len() != y.len()`.
pub fn gpu_axpy(a: f64, x: &[f64], y: &mut Vec<f64>) {
    assert_eq!(x.len(), y.len(), "gpu_axpy: length mismatch");
    for (yi, &xi) in y.iter_mut().zip(x.iter()) {
        *yi += a * xi;
    }
}

// ── Solvers ──────────────────────────────────────────────────────────────────

/// Conjugate gradient solver for symmetric positive-definite systems A·x = b.
///
/// Returns `(x, iterations, final_residual_norm)`.
/// Terminates when the relative residual ‖r‖/‖b‖ < `tol` or after `max_iter`
/// steps.
pub fn gpu_cg_solver(
    mat: &SparseMatrixGpu,
    b: &[f64],
    max_iter: usize,
    tol: f64,
) -> (Vec<f64>, usize, f64) {
    let n = mat.n;
    let mut x = vec![0.0f64; n];
    let mut r = b.to_vec(); // r = b - A*x; x=0 so r=b
    let mut p = r.clone();
    let mut rr = gpu_dot(&r, &r);
    let b_norm = rr.sqrt().max(1e-100);

    for iter in 0..max_iter {
        if rr.sqrt() / b_norm < tol {
            return (x, iter, rr.sqrt());
        }
        let ap = gpu_spmv(mat, &p);
        let pap = gpu_dot(&p, &ap);
        if pap.abs() < 1e-300 {
            break;
        }
        let alpha = rr / pap;
        gpu_axpy(alpha, &p, &mut x);
        gpu_axpy(-alpha, &ap, &mut r);
        let rr_new = gpu_dot(&r, &r);
        let beta = rr_new / rr.max(1e-300);
        // p = r + beta*p
        for i in 0..n {
            p[i] = r[i] + beta * p[i];
        }
        rr = rr_new;
    }
    (x, max_iter, rr.sqrt())
}

/// Build the Jacobi (diagonal) preconditioner for matrix `mat`.
///
/// Returns a vector `M⁻¹` where `M⁻¹[i] = 1/A[i,i]` (or 1 if the diagonal
/// entry is zero).
pub fn gpu_jacobi_preconditioner(mat: &SparseMatrixGpu) -> Vec<f64> {
    mat.diagonal()
        .iter()
        .map(|&d| if d.abs() > 1e-15 { 1.0 / d } else { 1.0 })
        .collect()
}

/// Preconditioned conjugate gradient solver.
///
/// `precond` should be a vector of reciprocal diagonal entries (from
/// [`gpu_jacobi_preconditioner`]).  Returns `(x, iterations, residual_norm)`.
pub fn gpu_pcg_solver(
    mat: &SparseMatrixGpu,
    b: &[f64],
    precond: &[f64],
    max_iter: usize,
    tol: f64,
) -> (Vec<f64>, usize, f64) {
    let n = mat.n;
    let mut x = vec![0.0f64; n];
    let mut r = b.to_vec();
    // z = M⁻¹ r
    let mut z: Vec<f64> = r
        .iter()
        .zip(precond.iter())
        .map(|(ri, mi)| ri * mi)
        .collect();
    let mut p = z.clone();
    let mut rz = gpu_dot(&r, &z);
    let b_norm = gpu_dot(b, b).sqrt().max(1e-100);

    for iter in 0..max_iter {
        if gpu_dot(&r, &r).sqrt() / b_norm < tol {
            return (x, iter, gpu_dot(&r, &r).sqrt());
        }
        let ap = gpu_spmv(mat, &p);
        let pap = gpu_dot(&p, &ap);
        if pap.abs() < 1e-300 {
            break;
        }
        let alpha = rz / pap;
        gpu_axpy(alpha, &p, &mut x);
        gpu_axpy(-alpha, &ap, &mut r);
        z = r
            .iter()
            .zip(precond.iter())
            .map(|(ri, mi)| ri * mi)
            .collect();
        let rz_new = gpu_dot(&r, &z);
        let beta = rz_new / rz.max(1e-300);
        for i in 0..n {
            p[i] = z[i] + beta * p[i];
        }
        rz = rz_new;
    }
    (x, max_iter, gpu_dot(&r, &r).sqrt())
}

// ── Stats ────────────────────────────────────────────────────────────────────

/// Performance statistics for a GPU sparse solve.
#[derive(Debug, Clone)]
pub struct GpuSparseSolverStats {
    /// Number of solver iterations performed.
    pub iterations: usize,
    /// ‖r‖ after the final iteration.
    pub final_residual: f64,
    /// Whether the solver converged within the requested tolerance.
    pub converged: bool,
    /// Wall-clock time in milliseconds (mock — always 0.0 in CPU simulation).
    pub time_ms: f64,
}

impl GpuSparseSolverStats {
    /// Create a new stats record.
    pub fn new(iterations: usize, final_residual: f64, converged: bool, time_ms: f64) -> Self {
        Self {
            iterations,
            final_residual,
            converged,
            time_ms,
        }
    }
}

// ============================================================
// Tests
// ============================================================
#[cfg(test)]
mod tests {
    use super::*;

    // ── helpers ─────────────────────────────────────────────────────────────

    /// Build a simple 3×3 SPD matrix:
    /// \[ 4 -1  0 \]
    /// \[-1  4 -1 \]
    /// \[ 0 -1  4 \]
    fn build_3x3_spd() -> SparseMatrixGpu {
        let row_ptr = vec![0, 2, 5, 7];
        let col_idx = vec![0, 1, 0, 1, 2, 1, 2];
        let values = vec![4.0, -1.0, -1.0, 4.0, -1.0, -1.0, 4.0];
        SparseMatrixGpu::from_csr(3, row_ptr, col_idx, values)
    }

    // ── identity ─────────────────────────────────────────────────────────────

    #[test]
    fn test_identity_spmv_returns_input() {
        let id = sparse_identity(4);
        let x = vec![1.0, 2.0, 3.0, 4.0];
        let y = gpu_spmv(&id, &x);
        for (yi, xi) in y.iter().zip(x.iter()) {
            assert!((yi - xi).abs() < 1e-12);
        }
    }

    #[test]
    fn test_identity_nnz() {
        let id = sparse_identity(5);
        assert_eq!(id.nnz(), 5);
    }

    #[test]
    fn test_identity_rows() {
        let id = sparse_identity(7);
        assert_eq!(id.rows(), 7);
    }

    #[test]
    fn test_identity_diagonal() {
        let id = sparse_identity(4);
        let diag = id.diagonal();
        assert_eq!(diag, vec![1.0; 4]);
    }

    #[test]
    fn test_identity_frobenius() {
        let n = 9usize;
        let id = sparse_identity(n);
        let expected = (n as f64).sqrt();
        assert!((id.frobenius_norm() - expected).abs() < 1e-10);
    }

    #[test]
    fn test_identity_is_symmetric() {
        assert!(sparse_identity(4).is_symmetric());
    }

    // ── diagonal matrix ──────────────────────────────────────────────────────

    #[test]
    fn test_diagonal_matrix_spmv() {
        let d = vec![2.0, 3.0, 5.0];
        let mat = sparse_diagonal_matrix(&d);
        let x = vec![1.0, 1.0, 1.0];
        let y = gpu_spmv(&mat, &x);
        assert!((y[0] - 2.0).abs() < 1e-12);
        assert!((y[1] - 3.0).abs() < 1e-12);
        assert!((y[2] - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_diagonal_matrix_frobenius() {
        let d = vec![3.0, 4.0];
        let mat = sparse_diagonal_matrix(&d);
        assert!((mat.frobenius_norm() - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_diagonal_matrix_is_symmetric() {
        let mat = sparse_diagonal_matrix(&[1.0, 2.0, 3.0]);
        assert!(mat.is_symmetric());
    }

    #[test]
    fn test_sparse_identity_zero_size() {
        let id = sparse_identity(0);
        assert_eq!(id.nnz(), 0);
        assert_eq!(id.rows(), 0);
    }

    // ── dot product ──────────────────────────────────────────────────────────

    #[test]
    fn test_gpu_dot_basic() {
        assert!((gpu_dot(&[1.0, 2.0, 3.0], &[4.0, 5.0, 6.0]) - 32.0).abs() < 1e-12);
    }

    #[test]
    fn test_gpu_dot_empty() {
        assert!((gpu_dot(&[], &[])).abs() < 1e-15);
    }

    #[test]
    fn test_gpu_dot_orthogonal() {
        assert!((gpu_dot(&[1.0, 0.0], &[0.0, 1.0])).abs() < 1e-15);
    }

    // ── axpy ─────────────────────────────────────────────────────────────────

    #[test]
    fn test_gpu_axpy_basic() {
        let x = vec![1.0, 2.0, 3.0];
        let mut y = vec![4.0, 5.0, 6.0];
        gpu_axpy(2.0, &x, &mut y);
        assert!((y[0] - 6.0).abs() < 1e-12);
        assert!((y[1] - 9.0).abs() < 1e-12);
        assert!((y[2] - 12.0).abs() < 1e-12);
    }

    #[test]
    fn test_gpu_axpy_zero_alpha() {
        let x = vec![1.0, 2.0];
        let mut y = vec![3.0, 4.0];
        gpu_axpy(0.0, &x, &mut y);
        assert!((y[0] - 3.0).abs() < 1e-12);
        assert!((y[1] - 4.0).abs() < 1e-12);
    }

    // ── spmv ─────────────────────────────────────────────────────────────────

    #[test]
    fn test_spmv_3x3_spd() {
        let mat = build_3x3_spd();
        let x = vec![1.0, 0.0, 0.0];
        let y = gpu_spmv(&mat, &x);
        assert!((y[0] - 4.0).abs() < 1e-12);
        assert!((y[1] + 1.0).abs() < 1e-12);
        assert!((y[2]).abs() < 1e-12);
    }

    #[test]
    fn test_spmv_zeros_input() {
        let mat = build_3x3_spd();
        let y = gpu_spmv(&mat, &[0.0, 0.0, 0.0]);
        for yi in y {
            assert!(yi.abs() < 1e-15);
        }
    }

    // ── from_csr / add_entry / finalize ──────────────────────────────────────

    #[test]
    fn test_from_csr_nnz() {
        let mat = build_3x3_spd();
        assert_eq!(mat.nnz(), 7);
    }

    #[test]
    fn test_add_entry_finalize() {
        let mut mat = SparseMatrixGpu::new(2);
        mat.add_entry(0, 0, 2.0);
        mat.add_entry(0, 1, -1.0);
        mat.add_entry(1, 0, -1.0);
        mat.add_entry(1, 1, 2.0);
        mat.finalize();
        assert_eq!(mat.nnz(), 4);
        let y = gpu_spmv(&mat, &[1.0, 0.0]);
        assert!((y[0] - 2.0).abs() < 1e-12);
        assert!((y[1] + 1.0).abs() < 1e-12);
    }

    // ── CG solver ────────────────────────────────────────────────────────────

    #[test]
    fn test_cg_converges_3x3() {
        let mat = build_3x3_spd();
        let b = vec![1.0, 0.0, 0.0];
        let (x, iters, res) = gpu_cg_solver(&mat, &b, 100, 1e-10);
        assert!(iters < 100, "CG did not converge: iters={iters}");
        // Verify A*x ≈ b
        let ax = gpu_spmv(&mat, &x);
        for (ai, bi) in ax.iter().zip(b.iter()) {
            assert!((ai - bi).abs() < 1e-8, "residual too large");
        }
        let _ = res; // residual checked via A*x
    }

    #[test]
    fn test_cg_identity_system() {
        let id = sparse_identity(3);
        let b = vec![1.0, 2.0, 3.0];
        let (x, _iters, _res) = gpu_cg_solver(&id, &b, 50, 1e-10);
        for (xi, bi) in x.iter().zip(b.iter()) {
            assert!((xi - bi).abs() < 1e-8);
        }
    }

    // ── Jacobi preconditioner ────────────────────────────────────────────────

    #[test]
    fn test_jacobi_preconditioner_diagonal() {
        let d = vec![2.0, 4.0, 5.0];
        let mat = sparse_diagonal_matrix(&d);
        let prec = gpu_jacobi_preconditioner(&mat);
        assert!((prec[0] - 0.5).abs() < 1e-12);
        assert!((prec[1] - 0.25).abs() < 1e-12);
        assert!((prec[2] - 0.2).abs() < 1e-12);
    }

    #[test]
    fn test_jacobi_preconditioner_identity() {
        let id = sparse_identity(3);
        let prec = gpu_jacobi_preconditioner(&id);
        for p in prec {
            assert!((p - 1.0).abs() < 1e-12);
        }
    }

    // ── PCG solver ───────────────────────────────────────────────────────────

    #[test]
    fn test_pcg_converges_3x3() {
        let mat = build_3x3_spd();
        let b = vec![1.0, 0.0, 1.0];
        let prec = gpu_jacobi_preconditioner(&mat);
        let (x, iters, _res) = gpu_pcg_solver(&mat, &b, &prec, 100, 1e-10);
        assert!(iters < 100);
        let ax = gpu_spmv(&mat, &x);
        for (ai, bi) in ax.iter().zip(b.iter()) {
            assert!((ai - bi).abs() < 1e-7);
        }
    }

    #[test]
    fn test_pcg_identity_trivial() {
        let id = sparse_identity(4);
        let b = vec![1.0, 2.0, 3.0, 4.0];
        let prec = gpu_jacobi_preconditioner(&id);
        let (x, _iters, _res) = gpu_pcg_solver(&id, &b, &prec, 10, 1e-12);
        for (xi, bi) in x.iter().zip(b.iter()) {
            assert!((xi - bi).abs() < 1e-8);
        }
    }

    // ── stats ────────────────────────────────────────────────────────────────

    #[test]
    fn test_solver_stats_fields() {
        let s = GpuSparseSolverStats::new(5, 1e-8, true, 0.0);
        assert_eq!(s.iterations, 5);
        assert!(s.converged);
        assert!((s.final_residual - 1e-8).abs() < 1e-20);
    }

    // ── is_symmetric ────────────────────────────────────────────────────────

    #[test]
    fn test_asymmetric_matrix() {
        // Upper triangular only — not symmetric
        let row_ptr = vec![0, 2, 3, 3];
        let col_idx = vec![0, 1, 1, 2];
        let values = vec![1.0, 2.0, 3.0, 4.0];
        let mat = SparseMatrixGpu::from_csr(3, row_ptr, col_idx, values);
        assert!(!mat.is_symmetric());
    }

    #[test]
    fn test_frobenius_zero_matrix() {
        let mat = SparseMatrixGpu::new(4);
        assert!((mat.frobenius_norm()).abs() < 1e-15);
    }

    // ── new() ────────────────────────────────────────────────────────────────

    #[test]
    fn test_new_empty_matrix() {
        let mat = SparseMatrixGpu::new(3);
        assert_eq!(mat.n, 3);
        assert_eq!(mat.nnz, 0);
        assert_eq!(mat.row_ptr, vec![0; 4]);
    }

    #[test]
    fn test_diagonal_of_empty_matrix() {
        let mat = SparseMatrixGpu::new(3);
        let diag = mat.diagonal();
        assert_eq!(diag, vec![0.0; 3]);
    }

    #[test]
    fn test_cg_zero_rhs() {
        let id = sparse_identity(3);
        let b = vec![0.0; 3];
        let (x, _iters, res) = gpu_cg_solver(&id, &b, 20, 1e-12);
        for xi in &x {
            assert!(xi.abs() < 1e-12);
        }
        assert!(res < 1e-10);
    }
}
