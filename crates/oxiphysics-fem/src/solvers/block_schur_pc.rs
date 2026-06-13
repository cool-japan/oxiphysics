// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Block-diagonal (Silvester-Wathen) preconditioner for Stokes saddle systems.
//!
//! For the saddle-point system arising from a Taylor-Hood (or any inf-sup
//! stable) Stokes discretisation,
//!
//! ```text
//! K = | A   Bᵀ |
//!     | B   0  |
//! ```
//!
//! the ideal block-diagonal preconditioner is `P = diag(A, S)` with the Schur
//! complement `S = B·A⁻¹·Bᵀ`.  By the Murphy-Golub-Wathen theorem the
//! preconditioned operator `P⁻¹K` then has only the three distinct eigenvalues
//! `{1, (1 ± √5)/2}`, so a Krylov method converges in a mesh-independent number
//! of iterations.
//!
//! This module provides a practical approximation:
//!
//! * the velocity block `A` is solved approximately by one AMG V-cycle
//!   ([`AmgPreconditioner`]);
//! * the Schur complement is replaced by the (lumped) pressure mass matrix,
//!   `S ≈ (1/ν)·Mp`, which is spectrally equivalent to `B·A⁻¹·Bᵀ` for inf-sup
//!   stable elements (Silvester & Wathen 1994).
//!
//! The accompanying [`gmres_left_preconditioned`] driver runs restarted GMRES
//! with any [`Preconditioner`], including [`BlockSchurPreconditioner`] (the
//! built-in [`ParallelGmresSolver`](crate::parallel_solver::ParallelGmresSolver)
//! is hard-wired to a Jacobi preconditioner and cannot host a custom one).

use crate::parallel_solver::{CsrMatrix, ParallelPcgSolver, PcgStats};
use crate::solvers::amg::classical::AmgClassical;
use crate::solvers::amg::cycle::CycleKind;
use crate::solvers::amg::preconditioner::{AmgPreconditioner, Preconditioner};

/// Threshold below which a value is treated as numerically zero.
const TINY: f64 = 1e-300;

/// Block-diagonal preconditioner for the Stokes saddle system `[[A, Bᵀ], [B, 0]]`.
///
/// The unknown vector is laid out velocity-first: indices `0..n_vel` are the
/// velocity DOFs and `n_vel..n_vel + n_pres` are the pressure DOFs.
///
/// * Velocity block: one AMG V-cycle approximating `A⁻¹`.
/// * Pressure block: `Ŝ⁻¹ ≈ (1/ν)·Mp⁻¹` using the lumped pressure-mass diagonal.
pub struct BlockSchurPreconditioner {
    /// Number of velocity degrees of freedom.
    n_vel: usize,
    /// Number of pressure degrees of freedom.
    n_pres: usize,
    /// Preconditioner for the velocity block (AMG V-cycle by default).
    vel_pc: Box<dyn Preconditioner>,
    /// Lumped pressure mass matrix diagonal (`Mp`), length `n_pres`.
    mp_diag: Vec<f64>,
    /// Kinematic viscosity `ν`.
    nu: f64,
}

impl BlockSchurPreconditioner {
    /// Build a preconditioner with an AMG V-cycle on the velocity block.
    ///
    /// `velocity_block` is the dense velocity-velocity block `A` (SPD,
    /// Laplacian-type); `mp_diag` is the lumped pressure mass diagonal of length
    /// `n_pres`; `nu` is the viscosity `ν`.
    pub fn new(
        velocity_block: &[Vec<f64>],
        mp_diag: Vec<f64>,
        nu: f64,
        n_vel: usize,
        n_pres: usize,
    ) -> Self {
        let a_csr = dense_to_csr(velocity_block);
        let hierarchy = AmgClassical::new().build(&a_csr);
        let pcg = ParallelPcgSolver::new(500, 1e-12);
        let vel_pc = AmgPreconditioner {
            hierarchy,
            cycle_kind: CycleKind::V,
            pcg,
        };
        Self {
            n_vel,
            n_pres,
            vel_pc: Box::new(vel_pc),
            mp_diag,
            nu,
        }
    }

    /// Build a preconditioner with a caller-supplied velocity-block preconditioner.
    ///
    /// Useful when a Jacobi (diagonal) velocity solve suffices, or to inject a
    /// pre-built AMG hierarchy.
    pub fn with_velocity_preconditioner(
        vel_pc: Box<dyn Preconditioner>,
        mp_diag: Vec<f64>,
        nu: f64,
        n_vel: usize,
        n_pres: usize,
    ) -> Self {
        Self {
            n_vel,
            n_pres,
            vel_pc,
            mp_diag,
            nu,
        }
    }
}

impl Preconditioner for BlockSchurPreconditioner {
    fn apply(&self, r: &[f64], z: &mut [f64]) {
        let (z_vel, z_pres) = z.split_at_mut(self.n_vel);
        let (r_vel, r_pres) = r.split_at(self.n_vel);

        // Velocity block: one AMG V-cycle ≈ A⁻¹.
        self.vel_pc.apply(r_vel, z_vel);

        // Pressure block: Ŝ⁻¹ ≈ (1/ν)·Mp⁻¹ via the lumped mass diagonal.
        for (zp, (rp, &mp)) in z_pres
            .iter_mut()
            .zip(r_pres.iter().zip(self.mp_diag.iter()))
        {
            let denom = self.nu * mp;
            *zp = if denom.abs() > TINY { *rp / denom } else { *rp };
        }
    }

    fn n(&self) -> usize {
        self.n_vel + self.n_pres
    }
}

/// Extract the dense velocity-velocity block `A` (top-left `n_vel × n_vel`) from
/// a saddle-point matrix stored as a flat `n_total × n_total` dense matrix.
pub fn extract_velocity_block(saddle_mat: &[Vec<f64>], n_vel: usize) -> Vec<Vec<f64>> {
    saddle_mat
        .iter()
        .take(n_vel)
        .map(|row| row.iter().take(n_vel).copied().collect())
        .collect()
}

/// Lump a (consistent) pressure mass matrix into its row-sum diagonal.
///
/// Row-sum lumping `Mp_lump[i] = Σⱼ Mp[i,j]` preserves the total mass and is
/// spectrally equivalent to the consistent mass matrix for low-order elements.
pub fn lump_pressure_mass(mp_block: &[Vec<f64>]) -> Vec<f64> {
    mp_block.iter().map(|row| row.iter().sum()).collect()
}

/// Convert a dense row-major matrix to CSR, dropping exact zeros.
fn dense_to_csr(a: &[Vec<f64>]) -> CsrMatrix {
    let n = a.len();
    let ncols = a.first().map_or(0, Vec::len);
    let mut row_offsets = vec![0usize; n + 1];
    let mut col_indices = Vec::new();
    let mut values = Vec::new();
    for (i, row) in a.iter().enumerate() {
        for (j, &v) in row.iter().enumerate() {
            if v != 0.0 {
                col_indices.push(j);
                values.push(v);
            }
        }
        row_offsets[i + 1] = col_indices.len();
    }
    CsrMatrix {
        nrows: n,
        ncols,
        row_offsets,
        col_indices,
        values,
    }
}

/// Sequential dot product.
fn dot(u: &[f64], v: &[f64]) -> f64 {
    u.iter().zip(v.iter()).map(|(a, b)| a * b).sum()
}

/// Restarted GMRES with an arbitrary left preconditioner.
///
/// Solves `A·x = b` (with `x` as the initial guess) using modified Gram-Schmidt
/// Arnoldi and Givens rotations, applying `pc` (e.g. a
/// [`BlockSchurPreconditioner`]) as a left preconditioner.  The relative
/// stopping test uses the preconditioned residual norm `‖M⁻¹r‖ / ‖M⁻¹b‖`.
///
/// Returns the iteration count, final preconditioned residual norm, and whether
/// the relative tolerance `tol` was met.
pub fn gmres_left_preconditioned(
    a: &CsrMatrix,
    b: &[f64],
    x: &mut [f64],
    pc: &dyn Preconditioner,
    krylov_dim: usize,
    max_restarts: usize,
    tol: f64,
) -> PcgStats {
    let n = a.nrows;

    // Preconditioned RHS norm for the relative stopping criterion.
    let mut mb = vec![0.0f64; n];
    pc.apply(b, &mut mb);
    let b_norm = dot(&mb, &mb).sqrt().max(TINY);

    let mut res_norm = 0.0f64;
    let mut total_iters = 0usize;

    for _restart in 0..max_restarts {
        // r0 = b − A·x ; z0 = M⁻¹·r0.
        let mut ax = vec![0.0f64; n];
        a.spmv(x, &mut ax);
        let r0: Vec<f64> = b.iter().zip(ax.iter()).map(|(bi, ai)| bi - ai).collect();
        let mut z0 = vec![0.0f64; n];
        pc.apply(&r0, &mut z0);

        let beta = dot(&z0, &z0).sqrt();
        res_norm = beta;
        if res_norm / b_norm < tol {
            break;
        }

        let m = krylov_dim.min(n).max(1);
        let inv_beta = 1.0 / beta.max(TINY);
        let mut q: Vec<Vec<f64>> = Vec::with_capacity(m + 1);
        q.push(z0.iter().map(|v| v * inv_beta).collect());

        let mut h = vec![vec![0.0f64; m]; m + 1];
        let mut cs = vec![0.0f64; m];
        let mut sn = vec![0.0f64; m];
        let mut e1 = vec![0.0f64; m + 1];
        e1[0] = beta;

        let mut j_stop = m;
        for j in 0..m {
            // w = M⁻¹·A·q_j.
            let mut aq = vec![0.0f64; n];
            a.spmv(&q[j], &mut aq);
            let mut w = vec![0.0f64; n];
            pc.apply(&aq, &mut w);

            // Modified Gram-Schmidt orthogonalisation.
            for i in 0..=j {
                let hij = dot(&w, &q[i]);
                h[i][j] = hij;
                for k in 0..n {
                    w[k] -= hij * q[i][k];
                }
            }
            let w_norm = dot(&w, &w).sqrt();
            h[j + 1][j] = w_norm;

            let exact = w_norm < 1e-14;
            if !exact {
                let inv = 1.0 / w_norm;
                q.push(w.iter().map(|v| v * inv).collect());
            }

            // Apply previous Givens rotations to column j.
            for i in 0..j {
                let tmp = cs[i] * h[i][j] + sn[i] * h[i + 1][j];
                h[i + 1][j] = -sn[i] * h[i][j] + cs[i] * h[i + 1][j];
                h[i][j] = tmp;
            }

            // New Givens rotation eliminating the sub-diagonal.
            let denom = (h[j][j] * h[j][j] + h[j + 1][j] * h[j + 1][j]).sqrt();
            if denom > TINY {
                cs[j] = h[j][j] / denom;
                sn[j] = h[j + 1][j] / denom;
            } else {
                cs[j] = 1.0;
                sn[j] = 0.0;
            }
            h[j][j] = cs[j] * h[j][j] + sn[j] * h[j + 1][j];
            h[j + 1][j] = 0.0;
            e1[j + 1] = -sn[j] * e1[j];
            e1[j] *= cs[j];

            res_norm = e1[j + 1].abs();
            if res_norm / b_norm < tol || exact {
                j_stop = j + 1;
                break;
            }
        }

        // Back-substitution for the least-squares solution y.
        let mut y = vec![0.0f64; j_stop];
        for i in (0..j_stop).rev() {
            y[i] = e1[i];
            for k in (i + 1)..j_stop {
                y[i] -= h[i][k] * y[k];
            }
            if h[i][i].abs() > TINY {
                y[i] /= h[i][i];
            }
        }

        // Update the solution: x += Q·y.
        for i in 0..j_stop {
            let yi = y[i];
            for k in 0..n {
                x[k] += yi * q[i][k];
            }
        }

        total_iters += j_stop;
        if res_norm / b_norm < tol {
            break;
        }
    }

    PcgStats {
        iterations: total_iters,
        residual_norm: res_norm,
        converged: res_norm / b_norm < tol,
    }
}

/// Solve a dense Stokes saddle system with a [`BlockSchurPreconditioner`].
///
/// Convenience wrapper that converts the dense `saddle` matrix to CSR and runs
/// [`gmres_left_preconditioned`].
pub fn block_schur_gmres(
    saddle: &[Vec<f64>],
    b: &[f64],
    x: &mut [f64],
    pc: &BlockSchurPreconditioner,
    krylov_dim: usize,
    max_restarts: usize,
    tol: f64,
) -> PcgStats {
    let csr = dense_to_csr(saddle);
    gmres_left_preconditioned(&csr, b, x, pc, krylov_dim, max_restarts, tol)
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn extract_velocity_block_top_left() {
        let saddle = vec![
            vec![4.0, 1.0, 9.0],
            vec![1.0, 4.0, 8.0],
            vec![9.0, 8.0, 0.0],
        ];
        let a = extract_velocity_block(&saddle, 2);
        assert_eq!(a, vec![vec![4.0, 1.0], vec![1.0, 4.0]]);
    }

    #[test]
    fn lump_pressure_mass_row_sums() {
        let mp = vec![vec![2.0, 0.5], vec![0.5, 2.0]];
        assert_eq!(lump_pressure_mass(&mp), vec![2.5, 2.5]);
    }

    #[test]
    fn preconditioner_dimension() {
        let a = vec![vec![2.0, 0.0], vec![0.0, 2.0]];
        let pc = BlockSchurPreconditioner::new(&a, vec![1.0], 1.0, 2, 1);
        assert_eq!(pc.n(), 3);
    }
}
