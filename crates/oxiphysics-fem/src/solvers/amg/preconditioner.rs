// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Preconditioner trait and AMG-based preconditioner.

#![allow(dead_code)]

use crate::parallel_solver::ParallelPcgSolver;
use crate::solvers::amg::cycle::{AmgHierarchy, CycleKind, v_cycle, w_cycle};

/// Generic preconditioner trait.
///
/// A preconditioner approximates the inverse of a system matrix.
/// It maps a residual `r` to a correction `z ≈ A^{-1} r`.
pub trait Preconditioner: Send + Sync {
    /// Apply the preconditioner: compute `z ≈ A^{-1} r`.
    fn apply(&self, r: &[f64], z: &mut [f64]);

    /// Return the problem dimension.
    fn n(&self) -> usize;
}

/// Jacobi (diagonal) preconditioner.
///
/// Computes `z[i] = r[i] / A[i,i]` (stored as `diag_inv[i] = 1/A[i,i]`).
#[derive(Debug, Clone)]
pub struct JacobiPreconditioner {
    /// Reciprocal diagonal values: `1/A[i,i]`.
    pub diag_inv: Vec<f64>,
}

impl Preconditioner for JacobiPreconditioner {
    fn apply(&self, r: &[f64], z: &mut [f64]) {
        for i in 0..r.len() {
            z[i] = r[i] * self.diag_inv[i];
        }
    }

    fn n(&self) -> usize {
        self.diag_inv.len()
    }
}

/// AMG preconditioner: applies one V-cycle or W-cycle as a preconditioner.
///
/// Used in PCG or GMRES to provide a near-optimal preconditioner for
/// elliptic problems.
pub struct AmgPreconditioner {
    /// The AMG level hierarchy.
    pub hierarchy: AmgHierarchy,
    /// Which cycle type to apply (V or W).
    pub cycle_kind: CycleKind,
    /// PCG solver used for coarse-level solves.
    pub pcg: ParallelPcgSolver,
}

impl Preconditioner for AmgPreconditioner {
    fn apply(&self, r: &[f64], z: &mut [f64]) {
        // Start from zero initial guess
        for zi in z.iter_mut() {
            *zi = 0.0;
        }
        // Run one cycle with b=r
        match self.cycle_kind {
            CycleKind::V => v_cycle(&self.hierarchy, 0, r, z, &self.pcg),
            CycleKind::W => w_cycle(&self.hierarchy, 0, r, z, &self.pcg),
        }
    }

    fn n(&self) -> usize {
        self.hierarchy.levels[0].a.nrows
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parallel_solver::CsrMatrix;
    use crate::solvers::amg::classical::AmgClassical;

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

    #[test]
    fn test_jacobi_preconditioner_apply() {
        // 4×4 diagonal matrix
        let a = CsrMatrix {
            nrows: 4,
            ncols: 4,
            row_offsets: vec![0, 1, 2, 3, 4],
            col_indices: vec![0, 1, 2, 3],
            values: vec![2.0, 3.0, 4.0, 5.0],
        };
        let diag_inv = a.diagonal_preconditioner();
        let jac = JacobiPreconditioner { diag_inv };

        let r = vec![2.0, 6.0, 8.0, 15.0];
        let mut z = vec![0.0f64; 4];
        jac.apply(&r, &mut z);

        // z[i] = r[i] / A[i,i]
        assert!((z[0] - 1.0).abs() < 1e-14, "z[0] = {}", z[0]);
        assert!((z[1] - 2.0).abs() < 1e-14, "z[1] = {}", z[1]);
        assert!((z[2] - 2.0).abs() < 1e-14, "z[2] = {}", z[2]);
        assert!((z[3] - 3.0).abs() < 1e-14, "z[3] = {}", z[3]);
    }

    #[test]
    fn test_amg_preconditioner_reduces() {
        // Build AMG on 64-node 1D Poisson with a real multilevel hierarchy.
        // coarse_cutoff=8 produces ~3 levels: 64→32→16→8.
        let n = 64;
        let a = make_1d_poisson(n);

        let mut amg_builder = AmgClassical::new();
        amg_builder.coarse_cutoff = 8; // real multilevel hierarchy
        let hier = amg_builder.build(&a);

        let pcg = ParallelPcgSolver::new(500, 1e-10);
        let precond = AmgPreconditioner {
            hierarchy: hier,
            cycle_kind: CycleKind::V,
            pcg,
        };

        // Use b = A*ones as the RHS so exact solution is x_true = ones.
        // A * ones for 1D Poisson (diag=2, off=-1) = [1, 0, 0, ..., 0, 1].
        let b: Vec<f64> = (0..n)
            .map(|i| if i == 0 || i == n - 1 { 1.0 } else { 0.0 })
            .collect();

        let mut z = vec![0.0f64; n];
        precond.apply(&b, &mut z);

        // After one V-cycle: z is a preconditioner correction, not the exact solution.
        // Verify: (1) z is non-trivial, (2) residual ||b - A*z|| < ||b|| (reduction).
        let z_norm: f64 = z.iter().map(|v| v * v).sum::<f64>().sqrt();
        assert!(
            z_norm > 1e-10,
            "AMG preconditioner produced trivially zero result: ||z|| = {z_norm}"
        );

        // Compute residual ||b - A*z||
        let mut az = vec![0.0f64; n];
        a.spmv(&z, &mut az);
        let res_norm: f64 = b
            .iter()
            .zip(az.iter())
            .map(|(bi, ai)| (bi - ai).powi(2))
            .sum::<f64>()
            .sqrt();
        let b_norm: f64 = b.iter().map(|v| v * v).sum::<f64>().sqrt();

        // One V-cycle should reduce the residual significantly (at least 10×)
        assert!(
            res_norm < b_norm * 0.1,
            "AMG preconditioner did not reduce residual: ||r||/||b|| = {:.3e}",
            res_norm / b_norm
        );
    }
}
