// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Symmetric positive definite algebraic-multigrid (AMG) V-cycle.
//!
//! A conjugate-gradient iteration needs a preconditioner `M⁻¹` that is itself
//! a *symmetric positive definite* operator, and it is only worth using if the
//! set-up cost is paid once rather than in every iteration.  [`SymmetricAmg`]
//! therefore builds its grid hierarchy a single time and then applies a fixed
//! V-cycle:
//!
//! * **coarsening** – greedy maximal independent set on the strength-of-connection
//!   graph ([`amg_coarsen`]) with direct interpolation ([`amg_prolongation`]);
//! * **coarse operators** – Galerkin products `A_{l+1} = Pᵀ A_l P`
//!   ([`amg_galerkin_coarse`]);
//! * **smoothing** – the same number of damped-Jacobi sweeps before and after
//!   the coarse-grid correction.  Jacobi is self-adjoint in the energy inner
//!   product, so the error propagation `S^ν (I − P A_c⁻¹ Pᵀ A) S^ν` is
//!   self-adjoint too and `M⁻¹` is symmetric.  The damping `ω = 4 / (3 ρ̂)`
//!   uses the Gershgorin bound `ρ̂ ≥ ρ(D⁻¹A)`, which makes every sweep a
//!   contraction in the energy norm and `M⁻¹` positive definite;
//! * **coarsest level** – an exact dense Cholesky solve once a level has at most
//!   [`DIRECT_SOLVE_LIMIT`] unknowns.  If coarsening stalls on a larger level
//!   (for instance a diagonal matrix, which has nothing to coarsen) the solve
//!   is a fixed number of symmetric Gauss-Seidel sweeps from a zero guess,
//!   which is again a symmetric positive definite operator.
//!
//! The matrix handed to [`SymmetricAmg::new`] must be symmetric positive
//! definite.

use super::amg::preconditioner::Preconditioner;
use super::types::{CholeskyDense, SolverError};
use crate::sparse::{CsrMatrix, amg_coarsen, amg_galerkin_coarse, amg_prolongation};

/// Largest level dimension that is solved exactly with a dense Cholesky factor.
pub const DIRECT_SOLVE_LIMIT: usize = 256;
/// Hard cap on the number of levels of the hierarchy.
const MAX_LEVELS: usize = 32;
/// Coarsening is regarded as stalled if it keeps more than this fraction of the unknowns.
const STALL_RATIO: f64 = 0.9;
/// Symmetric Gauss-Seidel sweeps used when the coarsest level is too large for a direct solve.
const COARSEST_SGS_SWEEPS: usize = 2;
/// Default number of damped-Jacobi sweeps before and after the coarse correction.
const DEFAULT_SWEEPS: usize = 2;

/// Solver used on the coarsest level of the hierarchy.
enum CoarsestSolver {
    /// Lower-triangular Cholesky factor (row-major) of the coarsest matrix.
    Direct(Vec<f64>),
    /// Fixed number of symmetric Gauss-Seidel sweeps from a zero guess.
    SymmetricGaussSeidel,
}

/// One level of the hierarchy.
struct Level {
    /// Level operator `A_l`.
    a: CsrMatrix,
    /// Reciprocal diagonal of `A_l`.
    diag_inv: Vec<f64>,
    /// Jacobi damping factor.
    omega: f64,
    /// Prolongation from the next coarser level (`None` on the coarsest level).
    prolongation: Option<CsrMatrix>,
}

/// Symmetric positive definite AMG V-cycle, usable as a CG preconditioner.
///
/// The hierarchy is built once by [`SymmetricAmg::new`]; every call of
/// [`SymmetricAmg::apply`] then runs one fixed V-cycle starting from a zero
/// guess, so `z = M⁻¹ r` is a linear, symmetric and positive definite function
/// of `r`.
pub struct SymmetricAmg {
    levels: Vec<Level>,
    coarsest: CoarsestSolver,
    sweeps: usize,
}

impl SymmetricAmg {
    /// Build the hierarchy for the SPD matrix `a` with the default smoothing
    /// (two damped-Jacobi sweeps before and after the coarse correction).
    ///
    /// # Errors
    /// * [`SolverError::DimensionMismatch`] if `a` is empty or not square.
    /// * [`SolverError::NotPositiveDefinite`] if a diagonal entry is not
    ///   positive or the coarsest-level Cholesky factorisation breaks down.
    pub fn new(a: &CsrMatrix) -> Result<Self, SolverError> {
        Self::with_sweeps(a, DEFAULT_SWEEPS)
    }

    /// Build the hierarchy with `sweeps` damped-Jacobi sweeps before and after
    /// the coarse correction (at least one sweep is always used: without a
    /// smoother the cycle would only be positive *semi*-definite).
    ///
    /// # Errors
    /// Same as [`SymmetricAmg::new`].
    pub fn with_sweeps(a: &CsrMatrix, sweeps: usize) -> Result<Self, SolverError> {
        if a.nrows == 0 || a.nrows != a.ncols {
            return Err(SolverError::DimensionMismatch);
        }
        let mut levels: Vec<Level> = Vec::new();
        let mut current = a.clone();
        let coarsest = loop {
            let n = current.nrows;
            let (diag_inv, omega) = jacobi_data(&current)?;
            if n <= DIRECT_SOLVE_LIMIT {
                let factor = dense_cholesky(&current)?;
                levels.push(Level {
                    a: current,
                    diag_inv,
                    omega,
                    prolongation: None,
                });
                break CoarsestSolver::Direct(factor);
            }
            let coarse_points = if levels.len() + 1 < MAX_LEVELS {
                amg_coarsen(&current)
            } else {
                Vec::new()
            };
            if coarse_points.is_empty() || coarse_points.len() as f64 >= STALL_RATIO * n as f64 {
                levels.push(Level {
                    a: current,
                    diag_inv,
                    omega,
                    prolongation: None,
                });
                break CoarsestSolver::SymmetricGaussSeidel;
            }
            let p = amg_prolongation(&current, &coarse_points);
            let next = amg_galerkin_coarse(&current, &p);
            levels.push(Level {
                a: current,
                diag_inv,
                omega,
                prolongation: Some(p),
            });
            current = next;
        };
        Ok(Self {
            levels,
            coarsest,
            sweeps: sweeps.max(1),
        })
    }

    /// Dimension of the fine-level system.
    pub fn n(&self) -> usize {
        self.levels.first().map_or(0, |lv| lv.a.nrows)
    }

    /// Number of levels in the hierarchy (a single level means a direct solve).
    pub fn n_levels(&self) -> usize {
        self.levels.len()
    }

    /// Number of unknowns on each level, finest first.
    pub fn level_sizes(&self) -> Vec<usize> {
        self.levels.iter().map(|lv| lv.a.nrows).collect()
    }

    /// Apply the preconditioner: write `z = M⁻¹ r` using one V-cycle from a
    /// zero initial guess.
    ///
    /// # Panics
    /// Panics if `r.len()` or `z.len()` differs from [`SymmetricAmg::n`].
    pub fn apply(&self, r: &[f64], z: &mut [f64]) {
        assert_eq!(
            r.len(),
            self.n(),
            "residual length must equal the system size"
        );
        assert_eq!(
            z.len(),
            self.n(),
            "output length must equal the system size"
        );
        self.cycle(0, r, z);
    }

    /// Recursive V-cycle on `level`, solving `A_level z ≈ r` from `z = 0`.
    fn cycle(&self, level: usize, r: &[f64], z: &mut [f64]) {
        let Some(lv) = self.levels.get(level) else {
            return;
        };
        let Some(p) = lv.prolongation.as_ref() else {
            self.solve_coarsest(lv, r, z);
            return;
        };
        z.fill(0.0);
        let mut scratch = vec![0.0_f64; r.len()];
        jacobi_sweeps(lv, r, z, self.sweeps, &mut scratch);
        // Coarse-grid correction: restrict the residual, solve, interpolate back.
        scratch.copy_from_slice(r);
        lv.a.mul_vec_axpby(z, &mut scratch, -1.0, 1.0);
        let mut r_coarse = vec![0.0_f64; p.ncols];
        for (i, &res_i) in scratch.iter().enumerate() {
            for idx in p.row_ptr[i]..p.row_ptr[i + 1] {
                r_coarse[p.col_indices[idx]] += p.values[idx] * res_i;
            }
        }
        let mut e_coarse = vec![0.0_f64; p.ncols];
        self.cycle(level + 1, &r_coarse, &mut e_coarse);
        for (i, z_i) in z.iter_mut().enumerate() {
            let mut correction = 0.0;
            for idx in p.row_ptr[i]..p.row_ptr[i + 1] {
                correction += p.values[idx] * e_coarse[p.col_indices[idx]];
            }
            *z_i += correction;
        }
        jacobi_sweeps(lv, r, z, self.sweeps, &mut scratch);
    }

    /// Solve the coarsest-level system `A z = r`.
    fn solve_coarsest(&self, lv: &Level, r: &[f64], z: &mut [f64]) {
        match &self.coarsest {
            CoarsestSolver::Direct(factor) => {
                let x = CholeskyDense::new(lv.a.nrows).solve(factor, r);
                z.copy_from_slice(&x);
            }
            CoarsestSolver::SymmetricGaussSeidel => {
                z.fill(0.0);
                for _ in 0..COARSEST_SGS_SWEEPS {
                    for i in 0..lv.a.nrows {
                        relax_row(&lv.a, r, z, i);
                    }
                    for i in (0..lv.a.nrows).rev() {
                        relax_row(&lv.a, r, z, i);
                    }
                }
            }
        }
    }
}

impl Preconditioner for SymmetricAmg {
    fn apply(&self, r: &[f64], z: &mut [f64]) {
        SymmetricAmg::apply(self, r, z);
    }

    fn n(&self) -> usize {
        SymmetricAmg::n(self)
    }
}

/// Reciprocal diagonal and Jacobi damping `ω = 4 / (3 ρ̂)` of `a`, where
/// `ρ̂ = max_i Σ_j |a_ij| / a_ii ≥ ρ(D⁻¹A)` is the Gershgorin bound.
fn jacobi_data(a: &CsrMatrix) -> Result<(Vec<f64>, f64), SolverError> {
    let mut diag_inv = Vec::with_capacity(a.nrows);
    let mut rho = 0.0_f64;
    for i in 0..a.nrows {
        let mut diag = 0.0_f64;
        let mut abs_sum = 0.0_f64;
        for idx in a.row_ptr[i]..a.row_ptr[i + 1] {
            let value = a.values[idx];
            abs_sum += value.abs();
            if a.col_indices[idx] == i {
                diag += value;
            }
        }
        if diag <= 0.0 || !diag.is_finite() {
            return Err(SolverError::NotPositiveDefinite);
        }
        diag_inv.push(1.0 / diag);
        rho = rho.max(abs_sum / diag);
    }
    Ok((diag_inv, 4.0 / (3.0 * rho)))
}

/// Lower Cholesky factor (row-major) of the CSR matrix `a`.
fn dense_cholesky(a: &CsrMatrix) -> Result<Vec<f64>, SolverError> {
    let n = a.nrows;
    let mut dense = vec![0.0_f64; n * n];
    for i in 0..n {
        for idx in a.row_ptr[i]..a.row_ptr[i + 1] {
            dense[i * n + a.col_indices[idx]] += a.values[idx];
        }
    }
    CholeskyDense::new(n).factorize(&dense)
}

/// `sweeps` damped-Jacobi sweeps `z += ω D⁻¹ (r − A z)`; `scratch` is a work vector.
fn jacobi_sweeps(lv: &Level, r: &[f64], z: &mut [f64], sweeps: usize, scratch: &mut [f64]) {
    for _ in 0..sweeps {
        scratch.copy_from_slice(r);
        lv.a.mul_vec_axpby(z, scratch, -1.0, 1.0);
        for ((z_i, &res_i), &d_inv) in z.iter_mut().zip(scratch.iter()).zip(lv.diag_inv.iter()) {
            *z_i += lv.omega * d_inv * res_i;
        }
    }
}

/// Gauss-Seidel relaxation of unknown `i` of `A x = rhs`, in place.
fn relax_row(a: &CsrMatrix, rhs: &[f64], x: &mut [f64], i: usize) {
    let mut diag = 0.0_f64;
    let mut off_diag = 0.0_f64;
    for idx in a.row_ptr[i]..a.row_ptr[i + 1] {
        let j = a.col_indices[idx];
        if j == i {
            diag += a.values[idx];
        } else {
            off_diag += a.values[idx] * x[j];
        }
    }
    if diag > 0.0 {
        x[i] = (rhs[i] - off_diag) / diag;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 2-D five-point Laplacian on an `m × m` grid with Dirichlet boundary.
    fn poisson_2d(m: usize) -> CsrMatrix {
        let mut triplets = Vec::new();
        for iy in 0..m {
            for ix in 0..m {
                let row = iy * m + ix;
                triplets.push((row, row, 4.0));
                if ix > 0 {
                    triplets.push((row, row - 1, -1.0));
                }
                if ix + 1 < m {
                    triplets.push((row, row + 1, -1.0));
                }
                if iy > 0 {
                    triplets.push((row, row - m, -1.0));
                }
                if iy + 1 < m {
                    triplets.push((row, row + m, -1.0));
                }
            }
        }
        CsrMatrix::from_triplets(m * m, m * m, &triplets)
    }

    /// Deterministic pseudo-random vector with entries in (-1, 1).
    fn pseudo_random(n: usize, seed: u64) -> Vec<f64> {
        let mut state = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (0..n)
            .map(|_| {
                state = state
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                ((state >> 11) as f64 / (1u64 << 53) as f64) * 2.0 - 1.0
            })
            .collect()
    }

    fn dot(a: &[f64], b: &[f64]) -> f64 {
        a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
    }

    fn apply_to(amg: &SymmetricAmg, r: &[f64]) -> Vec<f64> {
        let mut z = vec![0.0; r.len()];
        amg.apply(r, &mut z);
        z
    }

    #[test]
    fn hierarchy_coarsens_down_to_a_direct_solve() {
        let amg = SymmetricAmg::new(&poisson_2d(32)).expect("SPD Poisson matrix");
        let sizes = amg.level_sizes();
        assert_eq!(sizes[0], 1024);
        assert!(amg.n_levels() >= 2, "levels = {sizes:?}");
        for pair in sizes.windows(2) {
            assert!(pair[1] < pair[0], "levels must shrink: {sizes:?}");
        }
        assert!(
            sizes[sizes.len() - 1] <= DIRECT_SOLVE_LIMIT,
            "coarsest level too large: {sizes:?}"
        );
        assert_eq!(amg.n(), 1024);
    }

    #[test]
    fn small_system_is_solved_exactly() {
        let a = poisson_2d(8);
        let amg = SymmetricAmg::new(&a).expect("SPD Poisson matrix");
        assert_eq!(amg.n_levels(), 1);
        let r = pseudo_random(64, 7);
        let z = apply_to(&amg, &r);
        let az = a.mul_vec(&z);
        for (azi, ri) in az.iter().zip(r.iter()) {
            assert!((azi - ri).abs() < 1e-12, "A z = {azi}, r = {ri}");
        }
    }

    #[test]
    fn cycle_is_a_symmetric_operator() {
        let amg = SymmetricAmg::new(&poisson_2d(32)).expect("SPD Poisson matrix");
        for seed in 1..6 {
            let u = pseudo_random(1024, seed);
            let v = pseudo_random(1024, seed + 100);
            let u_mv = dot(&u, &apply_to(&amg, &v));
            let v_mu = dot(&v, &apply_to(&amg, &u));
            assert!(
                (u_mv - v_mu).abs() < 1e-12 * (u_mv.abs() + v_mu.abs() + 1.0),
                "seed {seed}: uᵀMv = {u_mv}, vᵀMu = {v_mu}"
            );
        }
    }

    #[test]
    fn cycle_is_positive_definite() {
        let amg = SymmetricAmg::new(&poisson_2d(32)).expect("SPD Poisson matrix");
        // Smooth, oscillatory and random residuals must all give uᵀ M u > 0.
        let smooth = vec![1.0; 1024];
        let checker: Vec<f64> = (0..1024)
            .map(|i| {
                if (i / 32 + i % 32) % 2 == 0 {
                    1.0
                } else {
                    -1.0
                }
            })
            .collect();
        let mut samples = vec![smooth, checker];
        for seed in 1..8 {
            samples.push(pseudo_random(1024, seed));
        }
        for u in &samples {
            let quad = dot(u, &apply_to(&amg, u));
            assert!(quad > 0.0, "uᵀMu = {quad}");
        }
    }

    #[test]
    fn cycle_is_linear_and_stateless() {
        let amg = SymmetricAmg::new(&poisson_2d(24)).expect("SPD Poisson matrix");
        let u = pseudo_random(576, 11);
        let v = pseudo_random(576, 12);
        let combo: Vec<f64> = u
            .iter()
            .zip(v.iter())
            .map(|(a, b)| 2.0 * a - 3.0 * b)
            .collect();
        let mu = apply_to(&amg, &u);
        let mv = apply_to(&amg, &v);
        let mc = apply_to(&amg, &combo);
        for i in 0..576 {
            let expected = 2.0 * mu[i] - 3.0 * mv[i];
            assert!((mc[i] - expected).abs() < 1e-12, "i = {i}");
        }
        assert_eq!(
            apply_to(&amg, &u),
            mu,
            "a second application must repeat exactly"
        );
    }

    #[test]
    fn stationary_iteration_with_the_cycle_converges_geometrically() {
        // x ← x + M⁻¹ (b − A x) contracts when M⁻¹ is a proper multigrid cycle.
        let a = poisson_2d(32);
        let amg = SymmetricAmg::new(&a).expect("SPD Poisson matrix");
        let b = vec![1.0; 1024];
        let mut x = vec![0.0; 1024];
        let mut norms = Vec::new();
        for _ in 0..9 {
            let ax = a.mul_vec(&x);
            let r: Vec<f64> = b.iter().zip(ax.iter()).map(|(bi, ai)| bi - ai).collect();
            norms.push(dot(&r, &r).sqrt());
            let z = apply_to(&amg, &r);
            for (xi, zi) in x.iter_mut().zip(z.iter()) {
                *xi += zi;
            }
        }
        assert!(norms[4] < 0.1 * norms[0], "norms = {norms:?}");
        assert!(norms[8] < 0.1 * norms[4], "norms = {norms:?}");
        assert!(norms[8] < 1e-2 * norms[0], "norms = {norms:?}");
    }

    #[test]
    fn diagonal_matrix_cannot_be_coarsened_but_is_still_handled() {
        // 600 > DIRECT_SOLVE_LIMIT unknowns with no couplings: coarsening
        // stalls and the symmetric Gauss-Seidel coarsest solve is used.
        let n = 600;
        let triplets: Vec<(usize, usize, f64)> = (0..n).map(|i| (i, i, 1.0 + i as f64)).collect();
        let a = CsrMatrix::from_triplets(n, n, &triplets);
        let amg = SymmetricAmg::new(&a).expect("SPD diagonal matrix");
        assert_eq!(amg.n_levels(), 1);
        let r = pseudo_random(n, 3);
        let z = apply_to(&amg, &r);
        for i in 0..n {
            assert!((z[i] - r[i] / (1.0 + i as f64)).abs() < 1e-14, "i = {i}");
        }
        assert!(dot(&r, &z) > 0.0);
    }

    #[test]
    fn rejects_matrices_that_are_not_spd() {
        let negative = CsrMatrix::from_triplets(2, 2, &[(0, 0, 1.0), (1, 1, -1.0)]);
        assert!(matches!(
            SymmetricAmg::new(&negative),
            Err(SolverError::NotPositiveDefinite)
        ));
        let indefinite =
            CsrMatrix::from_triplets(2, 2, &[(0, 0, 1.0), (0, 1, 3.0), (1, 0, 3.0), (1, 1, 1.0)]);
        assert!(matches!(
            SymmetricAmg::new(&indefinite),
            Err(SolverError::NotPositiveDefinite)
        ));
        let rectangular = CsrMatrix::from_triplets(2, 3, &[(0, 0, 1.0), (1, 1, 1.0)]);
        assert!(matches!(
            SymmetricAmg::new(&rectangular),
            Err(SolverError::DimensionMismatch)
        ));
        let empty = CsrMatrix::new(0, 0);
        assert!(matches!(
            SymmetricAmg::new(&empty),
            Err(SolverError::DimensionMismatch)
        ));
    }

    #[test]
    fn zero_sweeps_are_clamped_to_one() {
        let a = poisson_2d(32);
        let amg = SymmetricAmg::with_sweeps(&a, 0).expect("SPD Poisson matrix");
        let u = pseudo_random(1024, 21);
        assert!(dot(&u, &apply_to(&amg, &u)) > 0.0);
    }

    #[test]
    fn works_through_the_preconditioner_trait() {
        let a = poisson_2d(8);
        let amg = SymmetricAmg::new(&a).expect("SPD Poisson matrix");
        let pc: &dyn Preconditioner = &amg;
        assert_eq!(pc.n(), 64);
        let r = vec![1.0; 64];
        let mut z = vec![0.0; 64];
        pc.apply(&r, &mut z);
        assert_eq!(z, apply_to(&amg, &r));
    }
}
