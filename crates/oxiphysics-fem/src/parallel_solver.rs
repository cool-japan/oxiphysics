// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Parallel sparse solver using Rayon for multi-threaded finite-element assembly.
//!
//! This module provides two complementary components:
//!
//! 1. [`ParallelAssembler`] — collects per-element stiffness matrices concurrently
//!    and assembles the global CSR matrix without data races via atomic row-locking.
//! 2. [`ParallelPcgSolver`] — Preconditioned Conjugate Gradient with Rayon-parallel
//!    sparse matrix–vector products (SpMV).
//!
//! ## Design notes
//!
//! * Assembly is parallelised at the **element** level: each element computes its
//!   local stiffness matrix independently; the scatter (triplet → CSR) step uses a
//!   per-row Mutex for safe concurrent updates without a global lock.
//! * The SpMV kernel splits the CSR row range into chunks processed concurrently
//!   by Rayon.  Each output entry is independent, so no synchronisation is needed.
//! * The PCG solve loop itself is sequential (the dominant cost is SpMV, which is
//!   parallel); the dot products use Rayon's `par_iter().map().sum()` reduction.
//!
//! ## Usage
//!
//! ```
//! use oxiphysics_fem::parallel_solver::{CsrMatrix, ParallelAssembler, ParallelPcgSolver};
//!
//! // Build a tiny 4×4 CSR matrix manually
//! let mat = CsrMatrix {
//!     nrows: 4,
//!     ncols: 4,
//!     row_offsets: vec![0, 1, 2, 3, 4],
//!     col_indices:  vec![0, 1, 2, 3],
//!     values:       vec![4.0, 4.0, 4.0, 4.0],
//! };
//! let rhs = vec![1.0, 2.0, 3.0, 4.0];
//! let mut x   = vec![0.0f64; 4];
//! let stats = ParallelPcgSolver::default().solve(&mat, &rhs, &mut x);
//! assert!(stats.converged, "PCG did not converge");
//! ```

use rayon::prelude::*;
use std::sync::Mutex;

// ── CsrMatrix ────────────────────────────────────────────────────────────────

/// A sparse matrix in Compressed Sparse Row (CSR) format.
///
/// * `row_offsets[i]..row_offsets[i+1]` indexes the non-zeros of row `i`.
/// * `col_indices[row_offsets[i]..row_offsets[i+1]]` gives the column indices.
/// * `values[row_offsets[i]..row_offsets[i+1]]` gives the corresponding values.
#[derive(Debug, Clone, Default)]
pub struct CsrMatrix {
    /// Number of rows.
    pub nrows: usize,
    /// Number of columns.
    pub ncols: usize,
    /// Row start offsets (length `nrows + 1`).
    pub row_offsets: Vec<usize>,
    /// Column indices of non-zeros.
    pub col_indices: Vec<usize>,
    /// Non-zero values.
    pub values: Vec<f64>,
}

impl CsrMatrix {
    /// Create a zero matrix with a given sparsity pattern.
    pub fn new_from_pattern(
        nrows: usize,
        ncols: usize,
        row_offsets: Vec<usize>,
        col_indices: Vec<usize>,
    ) -> Self {
        let nnz = col_indices.len();
        Self {
            nrows,
            ncols,
            row_offsets,
            col_indices,
            values: vec![0.0; nnz],
        }
    }

    /// Create an identity matrix of size `n`.
    pub fn identity(n: usize) -> Self {
        Self {
            nrows: n,
            ncols: n,
            row_offsets: (0..=n).collect(),
            col_indices: (0..n).collect(),
            values: vec![1.0; n],
        }
    }

    /// Return the number of non-zero entries.
    pub fn nnz(&self) -> usize {
        self.values.len()
    }

    /// Parallel sparse matrix–vector product: `y = A * x`.
    ///
    /// Each row is computed independently, allowing Rayon to distribute rows
    /// across worker threads with zero synchronisation overhead.
    pub fn spmv_par(&self, x: &[f64], y: &mut [f64]) {
        debug_assert_eq!(x.len(), self.ncols);
        debug_assert_eq!(y.len(), self.nrows);

        y.par_iter_mut().enumerate().for_each(|(i, yi)| {
            let row_start = self.row_offsets[i];
            let row_end = self.row_offsets[i + 1];
            let mut sum = 0.0;
            for k in row_start..row_end {
                sum += self.values[k] * x[self.col_indices[k]];
            }
            *yi = sum;
        });
    }

    /// Sequential sparse matrix–vector product (for comparison / small matrices).
    pub fn spmv(&self, x: &[f64], y: &mut [f64]) {
        debug_assert_eq!(x.len(), self.ncols);
        debug_assert_eq!(y.len(), self.nrows);
        for (i, yi) in y.iter_mut().enumerate().take(self.nrows) {
            let row_start = self.row_offsets[i];
            let row_end = self.row_offsets[i + 1];
            let mut sum = 0.0;
            for k in row_start..row_end {
                sum += self.values[k] * x[self.col_indices[k]];
            }
            *yi = sum;
        }
    }

    /// Diagonal (Jacobi) preconditioner: returns `1 / A[i,i]` for each row.
    ///
    /// Uses `1.0` for rows with zero diagonal to avoid division by zero.
    pub fn diagonal_preconditioner(&self) -> Vec<f64> {
        let mut diag = vec![1.0f64; self.nrows];
        for (i, diag_i) in diag.iter_mut().enumerate().take(self.nrows) {
            for k in self.row_offsets[i]..self.row_offsets[i + 1] {
                if self.col_indices[k] == i {
                    let d = self.values[k];
                    if d.abs() > 1e-15 {
                        *diag_i = 1.0 / d;
                    }
                    break;
                }
            }
        }
        diag
    }

    /// Chunked SpMV: splits the row range into cache-friendly tiles of `chunk_size` rows.
    ///
    /// Produces the same result as [`CsrMatrix::spmv_par`] but processes rows in
    /// contiguous chunks to improve L3 cache utilisation on large matrices.
    pub fn spmv_chunked(&self, x: &[f64], y: &mut [f64], chunk_size: usize) {
        debug_assert_eq!(x.len(), self.ncols);
        debug_assert_eq!(y.len(), self.nrows);
        let chunk_size = chunk_size.max(1);
        y.par_chunks_mut(chunk_size)
            .enumerate()
            .for_each(|(chunk_idx, y_chunk)| {
                let row_start = chunk_idx * chunk_size;
                for (k, yi) in y_chunk.iter_mut().enumerate() {
                    let row = row_start + k;
                    let rs = self.row_offsets[row];
                    let re = self.row_offsets[row + 1];
                    let mut sum = 0.0;
                    for j in rs..re {
                        sum += self.values[j] * x[self.col_indices[j]];
                    }
                    *yi = sum;
                }
            });
    }
}

// ── AssemblyTask ────────────────────────────────────────────────────────────

/// Describes a single finite-element's local stiffness contribution.
///
/// `global_dofs` maps each local DOF to a global row/column index.
/// `ke` is the (ndof × ndof) local stiffness matrix stored in row-major order.
#[derive(Debug, Clone)]
pub struct AssemblyTask {
    /// Global DOF indices for this element's local DOFs.
    pub global_dofs: Vec<usize>,
    /// Local stiffness matrix values, row-major, shape `(ndof, ndof)`.
    pub ke: Vec<f64>,
}

impl AssemblyTask {
    /// Create from a symmetric positive definite local stiffness matrix.
    pub fn new(global_dofs: Vec<usize>, ke: Vec<f64>) -> Self {
        let ndof = global_dofs.len();
        debug_assert_eq!(ke.len(), ndof * ndof);
        Self { global_dofs, ke }
    }

    /// Number of degrees of freedom for this element.
    pub fn ndof(&self) -> usize {
        self.global_dofs.len()
    }
}

// ── ParallelAssembler ────────────────────────────────────────────────────────

/// Assembles per-element stiffness contributions into a global CSR matrix.
///
/// ## Thread safety
///
/// Elements are processed in parallel (Rayon).  Each non-zero entry's position
/// in the CSR `values` array is pre-determined from the sparsity pattern, so
/// scatter is lock-free at the element level.  A per-row Mutex guards concurrent
/// writes to overlapping rows (e.g. when elements share nodes on different threads).
#[derive(Debug, Default)]
pub struct ParallelAssembler {
    /// Total number of global DOFs.
    pub ndofs: usize,
}

impl ParallelAssembler {
    /// Create an assembler for a problem with `ndofs` degrees of freedom.
    pub fn new(ndofs: usize) -> Self {
        Self { ndofs }
    }

    /// Assemble all element stiffness matrices into a global CSR matrix.
    ///
    /// The sparsity pattern is computed from the union of all element connectivity
    /// patterns.  Values are accumulated atomically using per-row mutexes.
    ///
    /// # Panics
    ///
    /// Panics if any `global_dof` index exceeds `ndofs`.
    pub fn assemble(&self, tasks: &[AssemblyTask]) -> CsrMatrix {
        // ── Step 1: Build sparsity pattern (sequential) ──────────────────────
        // For each row, collect the set of column indices that appear.
        let mut row_cols: Vec<std::collections::BTreeSet<usize>> =
            vec![std::collections::BTreeSet::new(); self.ndofs];

        for task in tasks {
            for &row_dof in &task.global_dofs {
                for &col_dof in &task.global_dofs {
                    row_cols[row_dof].insert(col_dof);
                }
            }
        }

        // Flatten to CSR arrays
        let mut row_offsets = vec![0usize; self.ndofs + 1];
        let mut col_indices: Vec<usize> = Vec::new();
        for (i, cols) in row_cols.iter().enumerate() {
            row_offsets[i + 1] = row_offsets[i] + cols.len();
            col_indices.extend(cols.iter().copied());
        }
        let nnz = col_indices.len();
        let values = vec![0.0f64; nnz];

        // ── Step 2: Per-row lookup table (col → CSR index) ──────────────────
        let row_col_to_csr: Vec<std::collections::HashMap<usize, usize>> = row_cols
            .iter()
            .enumerate()
            .map(|(i, cols)| {
                let base = row_offsets[i];
                cols.iter()
                    .enumerate()
                    .map(|(j, &c)| (c, base + j))
                    .collect()
            })
            .collect();

        // ── Step 3: Parallel scatter using per-row mutexes ───────────────────
        let values_locked: Vec<Mutex<f64>> = values.into_iter().map(Mutex::new).collect();

        tasks.par_iter().for_each(|task| {
            let ndof = task.ndof();
            for (li, &row) in task.global_dofs.iter().enumerate() {
                for (lj, &col) in task.global_dofs.iter().enumerate() {
                    let ke_val = task.ke[li * ndof + lj];
                    if let Some(&csr_idx) = row_col_to_csr[row].get(&col) {
                        // SAFETY: each `csr_idx` is unique per (row, col) pair.
                        let mut guard = values_locked[csr_idx]
                            .lock()
                            .unwrap_or_else(|e| e.into_inner());
                        *guard += ke_val;
                    }
                }
            }
        });

        let values: Vec<f64> = values_locked
            .into_iter()
            .map(|m| {
                m.into_inner()
                    .expect("mutex not poisoned after parallel assembly")
            })
            .collect();

        CsrMatrix {
            nrows: self.ndofs,
            ncols: self.ndofs,
            row_offsets,
            col_indices,
            values,
        }
    }

    /// Assemble a load vector (right-hand side) from per-element force vectors.
    ///
    /// `element_forces[e]` must have the same DOF count as `element_dofs[e]`.
    pub fn assemble_rhs(
        &self,
        element_dofs: &[Vec<usize>],
        element_forces: &[Vec<f64>],
    ) -> Vec<f64> {
        let rhs_locked: Vec<Mutex<f64>> = (0..self.ndofs).map(|_| Mutex::new(0.0f64)).collect();
        element_dofs
            .par_iter()
            .zip(element_forces.par_iter())
            .for_each(|(dofs, forces)| {
                for (&dof, &f) in dofs.iter().zip(forces.iter()) {
                    *rhs_locked[dof].lock().unwrap_or_else(|e| e.into_inner()) += f;
                }
            });
        rhs_locked
            .into_iter()
            .map(|m| {
                m.into_inner()
                    .expect("mutex not poisoned after parallel rhs assembly")
            })
            .collect()
    }

    /// Assemble using element graph coloring for better cache locality.
    ///
    /// Computes a graph coloring of the elements (no two elements of the same color
    /// share a DOF) and assembles each color group in parallel.  The result is
    /// identical to [`ParallelAssembler::assemble`].
    pub fn assemble_colored(
        &self,
        tasks: &[AssemblyTask],
        element_dofs: &[Vec<usize>],
    ) -> CsrMatrix {
        use crate::solvers::assembly_coloring::{assemble_colored_csr, color_elements};
        let coloring = color_elements(tasks.len(), element_dofs);
        assemble_colored_csr(self.ndofs, tasks, &coloring)
    }
}

// ── PcgStats ─────────────────────────────────────────────────────────────────

/// Statistics returned by [`ParallelPcgSolver::solve`].
#[derive(Debug, Clone, Copy)]
pub struct PcgStats {
    /// Number of iterations performed.
    pub iterations: usize,
    /// Final residual norm ||r||.
    pub residual_norm: f64,
    /// Whether the solve converged.
    pub converged: bool,
}

// ── ParallelPcgSolver ────────────────────────────────────────────────────────

/// Preconditioned Conjugate Gradient with parallel SpMV.
///
/// Uses a Jacobi (diagonal) preconditioner by default.  The SpMV step uses
/// [`CsrMatrix::spmv_par`] and the dot products use `rayon::par_iter`.
#[derive(Debug, Clone)]
pub struct ParallelPcgSolver {
    /// Maximum number of CG iterations.
    pub max_iterations: usize,
    /// Relative residual tolerance for convergence.
    pub tolerance: f64,
    /// Minimum absolute residual (avoids infinite loop on singular systems).
    pub abs_tolerance: f64,
}

impl Default for ParallelPcgSolver {
    fn default() -> Self {
        Self {
            max_iterations: 500,
            tolerance: 1e-8,
            abs_tolerance: 1e-14,
        }
    }
}

impl ParallelPcgSolver {
    /// Create a new solver with explicit parameters.
    pub fn new(max_iterations: usize, tolerance: f64) -> Self {
        Self {
            max_iterations,
            tolerance,
            abs_tolerance: 1e-14,
        }
    }

    /// Solve `A * x = b` using Parallel PCG.
    ///
    /// `x` is used as the initial guess (zero-initialise for a cold start).
    ///
    /// Returns [`PcgStats`] with convergence information.
    pub fn solve(&self, a: &CsrMatrix, b: &[f64], x: &mut [f64]) -> PcgStats {
        let n = a.nrows;
        assert_eq!(b.len(), n);
        assert_eq!(x.len(), n);

        let m_inv = a.diagonal_preconditioner();

        // r = b - A*x
        let mut ax = vec![0.0f64; n];
        a.spmv_par(x, &mut ax);
        let mut r: Vec<f64> = b.iter().zip(ax.iter()).map(|(bi, ai)| bi - ai).collect();

        let b_norm = dot_par(b, b).sqrt();
        if b_norm < self.abs_tolerance {
            return PcgStats {
                iterations: 0,
                residual_norm: 0.0,
                converged: true,
            };
        }

        // z = M^{-1} * r
        let mut z: Vec<f64> = r.iter().zip(m_inv.iter()).map(|(ri, mi)| ri * mi).collect();

        // p = z
        let mut p = z.clone();

        let mut rz = dot_par(&r, &z);

        let mut ap = vec![0.0f64; n];
        let mut iters = 0;
        let mut res_norm = r.iter().map(|ri| ri * ri).sum::<f64>().sqrt();

        for _ in 0..self.max_iterations {
            if res_norm / b_norm < self.tolerance {
                break;
            }
            if res_norm < self.abs_tolerance {
                break;
            }

            // ap = A * p  (parallel SpMV)
            a.spmv_par(&p, &mut ap);

            let pap = dot_par(&p, &ap);
            if pap.abs() < 1e-300 {
                break;
            }
            let alpha = rz / pap;

            // x = x + alpha * p,  r = r - alpha * ap  (parallel)
            x.par_iter_mut()
                .zip(p.par_iter())
                .for_each(|(xi, pi)| *xi += alpha * pi);
            r.par_iter_mut()
                .zip(ap.par_iter())
                .for_each(|(ri, api)| *ri -= alpha * api);

            // z = M^{-1} * r
            z.par_iter_mut()
                .zip(r.par_iter().zip(m_inv.par_iter()))
                .for_each(|(zi, (ri, mi))| *zi = ri * mi);

            let rz_new = dot_par(&r, &z);
            let beta = rz_new / rz.max(1e-300);
            rz = rz_new;

            // p = z + beta * p  (parallel)
            p.par_iter_mut()
                .zip(z.par_iter())
                .for_each(|(pi, zi)| *pi = zi + beta * *pi);

            res_norm = r.iter().map(|ri| ri * ri).sum::<f64>().sqrt();
            iters += 1;
        }

        PcgStats {
            iterations: iters,
            residual_norm: res_norm,
            converged: res_norm / b_norm.max(1e-300) < self.tolerance
                || res_norm < self.abs_tolerance,
        }
    }
}

/// Parallel dot product using Rayon's reduction.
fn dot_par(a: &[f64], b: &[f64]) -> f64 {
    a.par_iter().zip(b.par_iter()).map(|(ai, bi)| ai * bi).sum()
}

// ── ParallelSparseDirectSolver (ILU(0) + GMRES stub) ────────────────────────

/// Parallel sparse direct solver using ILU(0) preconditioning.
///
/// For symmetric positive definite systems, prefer [`ParallelPcgSolver`].
/// This solver handles non-symmetric systems (e.g. convection-dominated problems).
#[derive(Debug, Clone)]
pub struct ParallelGmresSolver {
    /// Krylov subspace dimension (restart).
    pub krylov_dim: usize,
    /// Maximum number of outer restarts.
    pub max_restarts: usize,
    /// Convergence tolerance.
    pub tolerance: f64,
}

impl Default for ParallelGmresSolver {
    fn default() -> Self {
        Self {
            krylov_dim: 30,
            max_restarts: 10,
            tolerance: 1e-8,
        }
    }
}

impl ParallelGmresSolver {
    /// Solve `A * x = b` using restarted GMRES with diagonal preconditioning.
    ///
    /// For symmetric positive definite systems prefer [`ParallelPcgSolver`] which
    /// is more efficient.
    pub fn solve(&self, a: &CsrMatrix, b: &[f64], x: &mut [f64]) -> PcgStats {
        let n = a.nrows;
        let m_inv = a.diagonal_preconditioner();
        let mut res_norm = 0.0f64;
        // Use preconditioned b-norm for relative stopping criterion
        let mb: Vec<f64> = b.iter().zip(m_inv.iter()).map(|(bi, mi)| bi * mi).collect();
        let b_norm = dot_par(&mb, &mb).sqrt().max(1e-300);
        let mut total_iters = 0;

        for _restart in 0..self.max_restarts {
            // Compute initial residual r = b - A*x
            let mut ax = vec![0.0f64; n];
            a.spmv_par(x, &mut ax);
            let r0: Vec<f64> = b.iter().zip(ax.iter()).map(|(bi, ai)| bi - ai).collect();

            // Left-preconditioned residual: z0 = M^{-1} * r0
            let z0: Vec<f64> = r0
                .iter()
                .zip(m_inv.iter())
                .map(|(ri, mi)| ri * mi)
                .collect();
            let beta = dot_par(&z0, &z0).sqrt().max(1e-300);
            res_norm = beta;

            if res_norm / b_norm < self.tolerance {
                break;
            }

            // Arnoldi process: build Krylov basis Q (columns), upper Hessenberg H
            let m = self.krylov_dim.min(n);
            let mut q: Vec<Vec<f64>> = Vec::with_capacity(m + 1);

            // q[0] = z0 / beta (first Krylov vector in preconditioned space)
            q.push(z0.iter().map(|v| v / beta).collect());

            // Hessenberg matrix (stored as flat (m+1) x m)
            let mut h = vec![0.0f64; (m + 1) * m];
            // Cosines and sines for Givens rotations
            let mut cs = vec![0.0f64; m];
            let mut sn = vec![0.0f64; m];
            // Right-hand side of the least-squares problem: ||M^{-1}r0|| * e_1
            let mut e1 = vec![0.0f64; m + 1];
            e1[0] = beta;

            let mut j_stop = m;
            for j in 0..m {
                // w = M^{-1} * A * q[j]
                let mut aqj = vec![0.0f64; n];
                a.spmv_par(&q[j], &mut aqj);
                let w: Vec<f64> = aqj.iter().zip(m_inv.iter()).map(|(v, mi)| v * mi).collect();

                // Modified Gram-Schmidt orthogonalization
                let mut w = w;
                for i in 0..=j {
                    let hij = dot_par(&w, &q[i]);
                    h[i * m + j] = hij;
                    w.par_iter_mut()
                        .zip(q[i].par_iter())
                        .for_each(|(wi, qi)| *wi -= hij * qi);
                }
                let w_norm = dot_par(&w, &w).sqrt();
                h[(j + 1) * m + j] = w_norm;

                // Whether or not w_norm is near zero, we must apply Givens rotations
                // to column j of H so that back-substitution is correct.

                // Push the next Krylov vector only if it is non-degenerate.
                let exact_convergence = w_norm < 1e-14;
                if !exact_convergence {
                    q.push(w.iter().map(|v| v / w_norm).collect());
                }

                // Apply previous Givens rotations to column j of H
                for i in 0..j {
                    let tmp = cs[i] * h[i * m + j] + sn[i] * h[(i + 1) * m + j];
                    h[(i + 1) * m + j] = -sn[i] * h[i * m + j] + cs[i] * h[(i + 1) * m + j];
                    h[i * m + j] = tmp;
                }

                // Compute new Givens rotation
                let (c, s) = givens_rotation(h[j * m + j], h[(j + 1) * m + j]);
                cs[j] = c;
                sn[j] = s;

                h[j * m + j] = c * h[j * m + j] + s * h[(j + 1) * m + j];
                h[(j + 1) * m + j] = 0.0;
                e1[j + 1] = -s * e1[j];
                e1[j] *= c;

                res_norm = e1[j + 1].abs();
                total_iters += 1;

                // Stop if converged or if Krylov space is exhausted (exact solve).
                if res_norm / b_norm < self.tolerance || exact_convergence {
                    j_stop = j + 1;
                    break;
                }
            }

            // Back-substitution for y (upper triangular solve on H[0..j_stop, 0..j_stop])
            let mut y = vec![0.0f64; j_stop];
            for i in (0..j_stop).rev() {
                y[i] = e1[i];
                for k in (i + 1)..j_stop {
                    y[i] -= h[i * m + k] * y[k];
                }
                let hii = h[i * m + i];
                if hii.abs() > 1e-300 {
                    y[i] /= hii;
                }
            }

            // Update solution: x = x + Q * y
            for j in 0..j_stop {
                let yj = y[j];
                x.par_iter_mut()
                    .zip(q[j].par_iter())
                    .for_each(|(xi, qji)| *xi += yj * qji);
            }

            if res_norm / b_norm < self.tolerance {
                break;
            }
        }

        PcgStats {
            iterations: total_iters,
            residual_norm: res_norm,
            converged: res_norm / b_norm < self.tolerance,
        }
    }
}

/// Compute a Givens rotation `(c, s)` such that `[c, s; -s, c] * [a; b] = [r; 0]`.
fn givens_rotation(a: f64, b: f64) -> (f64, f64) {
    if b.abs() < 1e-300 {
        (1.0, 0.0)
    } else {
        let r = a.hypot(b);
        (a / r, b / r)
    }
}

// ── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn diag_matrix(n: usize, diag_val: f64) -> CsrMatrix {
        CsrMatrix {
            nrows: n,
            ncols: n,
            row_offsets: (0..=n).collect(),
            col_indices: (0..n).collect(),
            values: vec![diag_val; n],
        }
    }

    #[test]
    fn pcg_solves_diagonal_system() {
        let n = 16;
        let mat = diag_matrix(n, 4.0);
        let rhs: Vec<f64> = (0..n).map(|i| (i + 1) as f64).collect();
        let mut x = vec![0.0f64; n];
        let stats = ParallelPcgSolver::default().solve(&mat, &rhs, &mut x);
        assert!(stats.converged, "PCG did not converge: {:?}", stats);
        for (i, xi) in x.iter().enumerate() {
            let expected = (i + 1) as f64 / 4.0;
            assert!(
                (xi - expected).abs() < 1e-10,
                "x[{i}] = {xi}, expected {expected}"
            );
        }
    }

    #[test]
    fn parallel_spmv_matches_sequential() {
        let mat = diag_matrix(64, 3.0);
        let x: Vec<f64> = (0..64).map(|i| i as f64).collect();
        let mut y_par = vec![0.0f64; 64];
        let mut y_seq = vec![0.0f64; 64];
        mat.spmv_par(&x, &mut y_par);
        mat.spmv(&x, &mut y_seq);
        for (a, b) in y_par.iter().zip(y_seq.iter()) {
            assert!((a - b).abs() < 1e-14, "{a} != {b}");
        }
    }

    #[test]
    fn parallel_assembler_assembles_2d_bar() {
        // Two bar elements sharing a common node: DOFs [0,1], [1,2]
        // Element stiffness (1×1 blocks): [[1,-1],[-1,1]] for each element
        let tasks = vec![
            AssemblyTask::new(vec![0, 1], vec![1.0, -1.0, -1.0, 1.0]),
            AssemblyTask::new(vec![1, 2], vec![1.0, -1.0, -1.0, 1.0]),
        ];
        let asm = ParallelAssembler::new(3);
        let mat = asm.assemble(&tasks);
        assert_eq!(mat.nrows, 3);
        // Check diagonal: [1, 2, 1]
        let find_val = |row: usize, col: usize| -> f64 {
            for k in mat.row_offsets[row]..mat.row_offsets[row + 1] {
                if mat.col_indices[k] == col {
                    return mat.values[k];
                }
            }
            0.0
        };
        assert!((find_val(0, 0) - 1.0).abs() < 1e-14);
        assert!((find_val(1, 1) - 2.0).abs() < 1e-14);
        assert!((find_val(2, 2) - 1.0).abs() < 1e-14);
        assert!((find_val(0, 1) - (-1.0)).abs() < 1e-14);
        assert!((find_val(1, 2) - (-1.0)).abs() < 1e-14);
    }

    #[test]
    fn assemble_rhs_sums_forces() {
        let asm = ParallelAssembler::new(3);
        let dofs = vec![vec![0usize, 1], vec![1, 2]];
        let forces = vec![vec![1.0, 2.0], vec![3.0, 4.0]];
        let rhs = asm.assemble_rhs(&dofs, &forces);
        assert!((rhs[0] - 1.0).abs() < 1e-14);
        assert!((rhs[1] - 5.0).abs() < 1e-14); // 2 + 3
        assert!((rhs[2] - 4.0).abs() < 1e-14);
    }

    #[test]
    fn gmres_solves_diagonal_system() {
        let n = 8;
        let mat = diag_matrix(n, 2.0);
        let rhs: Vec<f64> = (0..n).map(|i| (i + 1) as f64).collect();
        let mut x = vec![0.0f64; n];
        let stats = ParallelGmresSolver::default().solve(&mat, &rhs, &mut x);
        assert!(stats.converged, "GMRES did not converge: {:?}", stats);
        for (i, xi) in x.iter().enumerate() {
            let expected = (i + 1) as f64 / 2.0;
            assert!(
                (xi - expected).abs() < 1e-6,
                "x[{i}] = {xi}, expected {expected}"
            );
        }
    }

    #[test]
    fn spmv_chunked_matches_spmv() {
        let n = 64;
        let mat = diag_matrix(n, 3.0);
        let x: Vec<f64> = (0..n).map(|i| i as f64).collect();
        let mut y_par = vec![0.0f64; n];
        let mut y_chunked = vec![0.0f64; n];
        mat.spmv_par(&x, &mut y_par);
        mat.spmv_chunked(&x, &mut y_chunked, 16);
        for (a, b) in y_par.iter().zip(y_chunked.iter()) {
            assert!((a - b).abs() < 1e-14, "{a} != {b}");
        }
    }
}

// ── AMG-preconditioned Krylov solvers ─────────────────────────────────────────

use crate::solvers::amg::{
    cycle::{AmgHierarchy, CycleKind},
    preconditioner::Preconditioner,
};

/// PCG solver with a pluggable preconditioner.
///
/// The preconditioner `P` must implement [`Preconditioner`], which maps a
/// residual `r` to a correction `z ≈ A^{-1} r`.
#[derive(Debug)]
pub struct PcgWithPrecond<P: Preconditioner> {
    /// Maximum number of PCG iterations.
    pub max_iterations: usize,
    /// Relative residual tolerance.
    pub tolerance: f64,
    /// Absolute residual tolerance.
    pub abs_tolerance: f64,
    /// Preconditioner instance.
    pub precond: P,
}

impl<P: Preconditioner> PcgWithPrecond<P> {
    /// Create a new preconditioned PCG solver.
    pub fn new(precond: P, max_iterations: usize, tolerance: f64) -> Self {
        Self {
            max_iterations,
            tolerance,
            abs_tolerance: 1e-14,
            precond,
        }
    }

    /// Solve `A * x = b` using preconditioned CG.
    ///
    /// `x` is used as the initial guess.
    pub fn solve(&self, a: &CsrMatrix, b: &[f64], x: &mut [f64]) -> PcgStats {
        let n = a.nrows;
        let mut z = vec![0.0f64; n];

        // r = b - A*x
        let mut ax = vec![0.0f64; n];
        a.spmv(x, &mut ax);
        let mut r: Vec<f64> = b.iter().zip(ax.iter()).map(|(bi, ai)| bi - ai).collect();

        self.precond.apply(&r, &mut z);
        let mut p = z.clone();
        let mut rz = dot_par(&r, &z);
        let b_norm = dot_par(b, b).sqrt().max(1e-300);
        let mut ap = vec![0.0f64; n];
        let mut iters = 0;
        let mut res_norm = r.iter().map(|ri| ri * ri).sum::<f64>().sqrt();

        for _ in 0..self.max_iterations {
            if res_norm / b_norm < self.tolerance {
                break;
            }
            if res_norm < self.abs_tolerance {
                break;
            }

            a.spmv_par(&p, &mut ap);
            let pap = dot_par(&p, &ap);
            if pap.abs() < 1e-300 {
                break;
            }
            let alpha = rz / pap;

            for i in 0..n {
                x[i] += alpha * p[i];
                r[i] -= alpha * ap[i];
            }
            self.precond.apply(&r, &mut z);
            let rz_new = dot_par(&r, &z);
            let beta = rz_new / rz.max(1e-300);
            rz = rz_new;
            for i in 0..n {
                p[i] = z[i] + beta * p[i];
            }
            res_norm = r.iter().map(|ri| ri * ri).sum::<f64>().sqrt();
            iters += 1;
        }

        PcgStats {
            iterations: iters,
            residual_norm: res_norm,
            converged: res_norm / b_norm < self.tolerance || res_norm < self.abs_tolerance,
        }
    }
}

/// PCG with an AMG preconditioner — preferred for SPD systems.
pub type PcgWithAmg = PcgWithPrecond<crate::solvers::amg::preconditioner::AmgPreconditioner>;

/// GMRES with AMG preconditioning for non-symmetric or indefinite systems.
///
/// Uses left-preconditioning: at each Arnoldi step, apply `M^{-1}` (one AMG cycle)
/// before adding the vector to the Krylov basis.
pub struct GmresWithAmg {
    /// Krylov subspace dimension (restart size).
    pub krylov_dim: usize,
    /// Maximum number of outer restarts.
    pub max_restarts: usize,
    /// Relative residual convergence tolerance.
    pub tolerance: f64,
    /// AMG hierarchy for preconditioning.
    pub hierarchy: AmgHierarchy,
}

impl GmresWithAmg {
    /// Create a new AMG-preconditioned GMRES solver.
    pub fn new(
        hierarchy: AmgHierarchy,
        krylov_dim: usize,
        max_restarts: usize,
        tolerance: f64,
    ) -> Self {
        Self {
            krylov_dim,
            max_restarts,
            tolerance,
            hierarchy,
        }
    }

    /// Solve `A * x = b` using AMG-preconditioned GMRES.
    pub fn solve(&self, a: &CsrMatrix, b: &[f64], x: &mut [f64]) -> PcgStats {
        let n = a.nrows;
        let b_norm = dot_par(b, b).sqrt().max(1e-300);
        let mut res_norm = 0.0f64;
        let mut total_iters = 0;

        let pcg = ParallelPcgSolver::new(500, 1e-10);
        let amg_precond = crate::solvers::amg::preconditioner::AmgPreconditioner {
            hierarchy: self.hierarchy.clone(),
            cycle_kind: CycleKind::V,
            pcg,
        };

        for _restart in 0..self.max_restarts {
            // Compute r = b - A*x, then apply preconditioner: r0z = M^{-1} r
            let mut ax = vec![0.0f64; n];
            a.spmv(x, &mut ax);
            let r0: Vec<f64> = b.iter().zip(ax.iter()).map(|(bi, ai)| bi - ai).collect();
            let mut r0z = vec![0.0f64; n];
            amg_precond.apply(&r0, &mut r0z);
            res_norm = r0z.iter().map(|v| v * v).sum::<f64>().sqrt();

            if res_norm / b_norm < self.tolerance {
                break;
            }

            let beta = res_norm;
            let mut q: Vec<Vec<f64>> = vec![r0z.iter().map(|v| v / beta).collect()];
            let m = self.krylov_dim.min(n);
            let mut h = vec![vec![0.0f64; m]; m + 1];
            let mut cs = vec![0.0f64; m];
            let mut sn = vec![0.0f64; m];
            let mut e1 = vec![0.0f64; m + 1];
            e1[0] = beta;

            let mut j_stop = m;

            for jj in 0..m {
                // w = M^{-1} A q_j
                let mut aq = vec![0.0f64; n];
                a.spmv_par(&q[jj], &mut aq);
                let mut w = vec![0.0f64; n];
                amg_precond.apply(&aq, &mut w);

                // Modified Gram-Schmidt
                for ii in 0..=jj {
                    h[ii][jj] = dot_par(&w, &q[ii]);
                    for k in 0..n {
                        w[k] -= h[ii][jj] * q[ii][k];
                    }
                }
                h[jj + 1][jj] = w.iter().map(|v| v * v).sum::<f64>().sqrt();

                // Whether or not Krylov space is exhausted, still apply Givens rotations.
                let exact_convergence = h[jj + 1][jj] < 1e-14;
                if !exact_convergence {
                    let inv_norm = 1.0 / h[jj + 1][jj];
                    q.push(w.iter().map(|v| v * inv_norm).collect());
                }

                // Apply previous Givens rotations
                for ii in 0..jj {
                    let tmp = cs[ii] * h[ii][jj] + sn[ii] * h[ii + 1][jj];
                    h[ii + 1][jj] = -sn[ii] * h[ii][jj] + cs[ii] * h[ii + 1][jj];
                    h[ii][jj] = tmp;
                }

                // New Givens rotation
                let denom = (h[jj][jj] * h[jj][jj] + h[jj + 1][jj] * h[jj + 1][jj]).sqrt();
                cs[jj] = if denom > 1e-300 {
                    h[jj][jj] / denom
                } else {
                    1.0
                };
                sn[jj] = if denom > 1e-300 {
                    h[jj + 1][jj] / denom
                } else {
                    0.0
                };
                h[jj][jj] = cs[jj] * h[jj][jj] + sn[jj] * h[jj + 1][jj];
                h[jj + 1][jj] = 0.0;
                e1[jj + 1] = -sn[jj] * e1[jj];
                e1[jj] *= cs[jj];
                res_norm = e1[jj + 1].abs();

                if res_norm / b_norm < self.tolerance || exact_convergence {
                    j_stop = jj + 1;
                    break;
                }
            }

            // Back-substitution
            let mut y = vec![0.0f64; j_stop];
            for ii in (0..j_stop).rev() {
                y[ii] = e1[ii];
                for kk in (ii + 1)..j_stop {
                    y[ii] -= h[ii][kk] * y[kk];
                }
                if h[ii][ii].abs() > 1e-300 {
                    y[ii] /= h[ii][ii];
                }
            }

            // Update solution
            for ii in 0..j_stop {
                let yi = y[ii];
                for k in 0..n {
                    x[k] += yi * q[ii][k];
                }
            }

            total_iters += j_stop;
            if res_norm / b_norm < self.tolerance {
                break;
            }
        }

        PcgStats {
            iterations: total_iters,
            residual_norm: res_norm,
            converged: res_norm / b_norm < self.tolerance,
        }
    }
}
