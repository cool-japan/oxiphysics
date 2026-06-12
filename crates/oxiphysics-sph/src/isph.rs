// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Incompressible SPH (ISPH) via the Cummins–Rudman pressure-projection method.
//!
//! Unlike WCSPH, which enforces (weak) incompressibility through a stiff
//! explicit Tait equation of state, ISPH solves a *pressure-Poisson equation*
//! (PPE) implicitly each time step.  The projection splits the update into a
//! predictor that ignores pressure and a corrector that removes the divergent
//! part of the predicted velocity field.  This removes the acoustic CFL limit
//! of WCSPH and yields a smooth, noise-free pressure field.
//!
//! Reference: Cummins, S. J. & Rudman, M. (1999), *An SPH Projection Method*,
//! Journal of Computational Physics 152(2), 584–607.
//!
//! # Algorithm (one step)
//!
//! 1. **Predict** velocity `v* = vⁿ + Δt (g + ν∇²vⁿ)` (no pressure force).
//! 2. **Divergence** `∇·v*ᵢ = Σⱼ (mⱼ/ρⱼ)(v*ⱼ − v*ᵢ)·∇Wᵢⱼ`.
//! 3. **Assemble** the symmetric PPE operator `A p = b` as a CSR matrix:
//!    - `A ≈ −∇·((1/ρ)∇p)`, a symmetric positive-definite M-matrix with one
//!      Dirichlet node per free surface (`ρᵢ < 0.95 ρ₀ ⇒ pᵢ = 0`),
//!    - `bᵢ = −∇·v*ᵢ / Δt` (the projection-consistent Cummins–Rudman source).
//! 4. **Solve** with a self-contained PCG + Jacobi preconditioner.
//! 5. **Pressure accel** `aₚᵢ = −Σⱼ mⱼ (pᵢ/ρᵢ² + pⱼ/ρⱼ²) ∇Wᵢⱼ`.
//! 6. **Correct** velocity `v^(n+1) = v* + Δt aₚ`.
//! 7. **Integrate** positions `x^(n+1) = xⁿ + Δt v^(n+1)`.
//!
//! The Laplacian off-diagonal uses the Brookshaw / Morris finite-difference
//! form `(rᵢⱼ·∇Wᵢⱼ)/(|rᵢⱼ|² + η²)`.  The diagonal is built as the negative row
//! sum so the discrete operator annihilates a constant pressure field, exactly
//! as the continuous Laplacian does.  Walls use the natural (homogeneous
//! Neumann) boundary condition by participating in the PPE with `v* = 0`.
//!
//! The linear algebra (`CsrMatrix`, [`pcg_solve`]) is intentionally
//! self-contained — no external sparse-solver dependency — so the crate stays
//! pure Rust and free of cross-crate coupling.

use crate::kernel::{CubicSplineKernel, SphKernel, grad};

// ---------------------------------------------------------------------------
// Public configuration / result / error types
// ---------------------------------------------------------------------------

/// Configuration for a single ISPH projection step.
#[derive(Debug, Clone)]
pub struct IsphConfig {
    /// Time-step Δt (s).
    pub dt: f64,
    /// Smoothing length `h` (kernel radius parameter; support is `2h`).
    pub smoothing_length: f64,
    /// Reference (rest) density ρ₀ (kg/m³).
    pub rho0: f64,
    /// Kinematic viscosity ν (m²/s) for the Morris viscous predictor term.
    pub viscosity: f64,
    /// Body acceleration (e.g. gravity) applied to every fluid particle (m/s²).
    pub gravity: [f64; 3],
    /// Relative residual tolerance for the PCG solve.
    pub pcg_tol: f64,
    /// Maximum number of PCG iterations.
    pub pcg_max_iters: usize,
}

impl Default for IsphConfig {
    fn default() -> Self {
        Self {
            dt: 1.0e-3,
            smoothing_length: 1.0,
            rho0: 1000.0,
            viscosity: 0.0,
            gravity: [0.0, 0.0, -9.81],
            pcg_tol: 1.0e-4,
            pcg_max_iters: 50,
        }
    }
}

/// Diagnostics returned from a single ISPH step.
#[derive(Debug, Clone)]
pub struct IsphResult {
    /// Per-particle pressure solved from the PPE (Pa).
    pub pressures: Vec<f64>,
    /// Number of PCG iterations performed.
    pub pcg_iters: usize,
    /// Final relative PCG residual `‖b − Ap‖ / ‖b‖`.
    pub pcg_residual: f64,
    /// RMS of `∇·v^(n+1)` over interior fluid particles (post-projection).
    pub divergence_rms: f64,
    /// RMS of `∇·v*` over interior fluid particles (pre-projection).
    pub predicted_divergence_rms: f64,
}

/// Errors that can arise during an ISPH step.
#[derive(Debug, thiserror::Error)]
pub enum IsphError {
    /// The PPE solve produced a non-finite residual or pressure.
    #[error("ISPH PPE PCG did not converge: {0} iters, residual {1:.3e}")]
    PcgNonConvergence(usize, f64),
    /// The particle list was empty.
    #[error("empty particle list")]
    NoParticles,
    /// `smoothing_length`, `dt`, or `rho0` was not strictly positive.
    #[error("smoothing_length must be > 0")]
    InvalidConfig,
}

/// Density fraction below which a particle is treated as a free surface and
/// pinned to `p = 0` (Dirichlet boundary condition).
const FREE_SURFACE_FRACTION: f64 = 0.95;

// ---------------------------------------------------------------------------
// Compressed-sparse-row matrix (self-contained)
// ---------------------------------------------------------------------------

/// Minimal compressed-sparse-row (CSR) matrix.
///
/// Self-contained sparse-matrix building block for the pressure-Poisson solve
/// (no external sparse-solver dependency), consumed by [`pcg_solve`].
#[derive(Debug, Clone)]
pub(crate) struct CsrMatrix {
    /// Number of rows (= columns; the matrix is square).
    nrows: usize,
    /// Row pointers of length `nrows + 1`.
    row_ptr: Vec<usize>,
    /// Column indices, one per stored value.
    col_idx: Vec<usize>,
    /// Stored matrix values, aligned with `col_idx`.
    values: Vec<f64>,
}

impl CsrMatrix {
    /// Assemble a CSR matrix from per-row off-diagonal entries plus an explicit
    /// diagonal.  Each row's columns are emitted in ascending order with the
    /// diagonal inserted at its sorted position.
    fn from_rows(nrows: usize, off_diagonals: &[Vec<(usize, f64)>], diagonal: &[f64]) -> Self {
        let mut row_ptr = Vec::with_capacity(nrows + 1);
        let mut col_idx = Vec::new();
        let mut values = Vec::new();
        row_ptr.push(0);
        for i in 0..nrows {
            // Merge off-diagonals with the diagonal entry, then sort by column.
            let mut entries: Vec<(usize, f64)> = Vec::with_capacity(off_diagonals[i].len() + 1);
            entries.extend_from_slice(&off_diagonals[i]);
            entries.push((i, diagonal[i]));
            entries.sort_by_key(|&(c, _)| c);
            for (c, v) in entries {
                col_idx.push(c);
                values.push(v);
            }
            row_ptr.push(col_idx.len());
        }
        Self {
            nrows,
            row_ptr,
            col_idx,
            values,
        }
    }

    /// Sparse matrix–vector product `y = A x`.
    fn matvec(&self, x: &[f64], y: &mut [f64]) {
        for (i, yi) in y.iter_mut().enumerate() {
            let mut sum = 0.0_f64;
            for k in self.row_ptr[i]..self.row_ptr[i + 1] {
                sum += self.values[k] * x[self.col_idx[k]];
            }
            *yi = sum;
        }
    }

    /// Extract the main diagonal of the matrix.
    fn diag(&self) -> Vec<f64> {
        let mut d = vec![0.0_f64; self.nrows];
        for (i, di) in d.iter_mut().enumerate() {
            for k in self.row_ptr[i]..self.row_ptr[i + 1] {
                if self.col_idx[k] == i {
                    *di = self.values[k];
                    break;
                }
            }
        }
        d
    }

    /// Number of rows (equivalently, columns) in the square matrix.
    #[cfg(test)]
    fn nrows(&self) -> usize {
        self.nrows
    }
}

// ---------------------------------------------------------------------------
// Preconditioned conjugate gradient (Jacobi preconditioner)
// ---------------------------------------------------------------------------

/// Dot product of two equal-length slices.
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(&x, &y)| x * y).sum()
}

/// Solve the symmetric positive-definite system `A x = b` (with `A` stored as a
/// [`CsrMatrix`]) using a Jacobi-preconditioned conjugate gradient method.
///
/// Returns `(iterations, final_relative_residual)`.  `x` is used as the
/// initial guess and overwritten with the solution.
pub(crate) fn pcg_solve(
    a: &CsrMatrix,
    b: &[f64],
    x: &mut [f64],
    tol: f64,
    max_iters: usize,
) -> (usize, f64) {
    let n = b.len();
    let bnorm = dot(b, b).sqrt();
    if bnorm <= f64::MIN_POSITIVE {
        // Zero right-hand side ⇒ zero solution.
        for xi in x.iter_mut() {
            *xi = 0.0;
        }
        return (0, 0.0);
    }

    let inv_diag: Vec<f64> = a
        .diag()
        .iter()
        .map(|&d| if d.abs() > 1e-30 { 1.0 / d } else { 0.0 })
        .collect();

    let mut r = vec![0.0_f64; n];
    a.matvec(x, &mut r);
    for i in 0..n {
        r[i] = b[i] - r[i];
    }

    let mut res = dot(&r, &r).sqrt() / bnorm;
    if res <= tol {
        return (0, res);
    }

    let mut z: Vec<f64> = (0..n).map(|i| r[i] * inv_diag[i]).collect();
    let mut p = z.clone();
    let mut rz_old = dot(&r, &z);
    let mut ap = vec![0.0_f64; n];
    let mut iters = 0;

    for k in 0..max_iters {
        iters = k + 1;
        a.matvec(&p, &mut ap);
        let pap = dot(&p, &ap);
        if pap.abs() < 1e-300 {
            break;
        }
        let alpha = rz_old / pap;
        for i in 0..n {
            x[i] += alpha * p[i];
            r[i] -= alpha * ap[i];
        }
        res = dot(&r, &r).sqrt() / bnorm;
        if res <= tol {
            break;
        }
        for i in 0..n {
            z[i] = r[i] * inv_diag[i];
        }
        let rz_new = dot(&r, &z);
        if rz_old.abs() < 1e-300 {
            break;
        }
        let beta = rz_new / rz_old;
        for i in 0..n {
            p[i] = z[i] + beta * p[i];
        }
        rz_old = rz_new;
    }

    (iters, res)
}

// ---------------------------------------------------------------------------
// Small vector helpers
// ---------------------------------------------------------------------------

#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// RMS of `values[i]` taken only over indices where `mask[i]` is `true`.
fn masked_rms(values: &[f64], mask: &[bool]) -> f64 {
    let mut sum = 0.0_f64;
    let mut count = 0_usize;
    for (i, &v) in values.iter().enumerate() {
        if mask[i] {
            sum += v * v;
            count += 1;
        }
    }
    if count == 0 {
        return 0.0;
    }
    (sum / count as f64).sqrt()
}


/// Shared per-step SPH field context (frozen positions, masses, densities,
/// neighbor lists, kernel, and smoothing length).  Bundling these keeps the
/// per-particle operators to a small argument count.
struct SphField<'a> {
    positions: &'a [[f64; 3]],
    masses: &'a [f64],
    densities: &'a [f64],
    neighbors: &'a [Vec<usize>],
    kernel: CubicSplineKernel,
    h: f64,
}

impl SphField<'_> {
    /// SPH velocity divergence at particle `i`:
    /// `∇·vᵢ = Σⱼ (mⱼ/ρⱼ)(vⱼ − vᵢ)·∇Wᵢⱼ`.
    fn divergence(&self, i: usize, velocities: &[[f64; 3]]) -> f64 {
        let mut div = 0.0_f64;
        for &j in &self.neighbors[i] {
            if j == i {
                continue;
            }
            let rij = sub3(self.positions[i], self.positions[j]);
            if dot3(rij, rij) < 1e-24 {
                continue;
            }
            let gw = grad(&self.kernel, rij, self.h);
            let dv = sub3(velocities[j], velocities[i]);
            div += (self.masses[j] / self.densities[j]) * dot3(dv, gw);
        }
        div
    }

    /// Momentum-conserving symmetric SPH pressure acceleration at particle `i`:
    /// `aₚᵢ = −Σⱼ mⱼ (pᵢ/ρᵢ² + pⱼ/ρⱼ²) ∇Wᵢⱼ ≈ −(1/ρ)∇p`.
    fn pressure_accel(&self, i: usize, pressures: &[f64]) -> [f64; 3] {
        let mut accel = [0.0_f64; 3];
        let pi_term = pressures[i] / (self.densities[i] * self.densities[i]);
        for &j in &self.neighbors[i] {
            if j == i {
                continue;
            }
            let rij = sub3(self.positions[i], self.positions[j]);
            if dot3(rij, rij) < 1e-24 {
                continue;
            }
            let gw = grad(&self.kernel, rij, self.h);
            let pj_term = pressures[j] / (self.densities[j] * self.densities[j]);
            let coeff = self.masses[j] * (pi_term + pj_term);
            accel[0] -= coeff * gw[0];
            accel[1] -= coeff * gw[1];
            accel[2] -= coeff * gw[2];
        }
        accel
    }

    /// Assemble the symmetric Brookshaw/Morris discrete Laplacian
    /// `A ≈ −∇·((1/ρ)∇)` as a CSR matrix.
    ///
    /// Off-diagonal `Aᵢⱼ = −cᵢⱼ` and diagonal `Aᵢᵢ = Σⱼ cᵢⱼ` with the symmetric
    /// coupling `cᵢⱼ = cⱼᵢ = (mᵢ+mⱼ)/(ρᵢρⱼ) · −(rᵢⱼ·∇Wᵢⱼ)/(|rᵢⱼ|²+η²) > 0`.  This
    /// is a positive-definite M-matrix (after pinning at least one Dirichlet
    /// node), so Jacobi-preconditioned CG converges quickly — unlike the exact
    /// divergence-of-gradient operator, whose collocated checkerboard null space
    /// is badly conditioned.  The residual inconsistency between this Laplacian
    /// and the pressure-gradient correction is removed by the outer projection
    /// iteration in [`isph_step_with_walls`].  Dirichlet rows (`p = 0` free
    /// surfaces) are reduced to identity with their couplings eliminated
    /// symmetrically.
    fn brookshaw(&self, dirichlet: &[bool], eta2: f64) -> CsrMatrix {
        let n = self.positions.len();
        let mut off_diagonals: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
        let mut diagonal = vec![0.0_f64; n];
        for i in 0..n {
            if dirichlet[i] {
                diagonal[i] = 1.0;
                continue;
            }
            for &j in &self.neighbors[i] {
                if j == i {
                    continue;
                }
                let rij = sub3(self.positions[i], self.positions[j]);
                let r2 = dot3(rij, rij);
                if r2 < 1e-24 {
                    continue;
                }
                let gw = grad(&self.kernel, rij, self.h);
                let kappa = dot3(rij, gw) / (r2 + eta2); // ≤ 0
                let c =
                    (self.masses[i] + self.masses[j]) / (self.densities[i] * self.densities[j])
                        * (-kappa);
                diagonal[i] += c;
                if !dirichlet[j] {
                    off_diagonals[i].push((j, -c));
                }
            }
            if diagonal[i].abs() < 1e-30 {
                // Isolated particle: trivial p = 0 row.
                diagonal[i] = 1.0;
                off_diagonals[i].clear();
            }
        }
        CsrMatrix::from_rows(n, &off_diagonals, &diagonal)
    }
}

// ---------------------------------------------------------------------------
// Public ISPH step
// ---------------------------------------------------------------------------

/// Advance the fluid one ISPH (Cummins–Rudman) time step.
///
/// All particles are treated as fluid (gravity applied to every particle).
/// Positions and velocities are mutated in place; densities are recomputed by
/// SPH summation and written back.  `neighbor_fn(i, h)` must return the indices
/// of the particles within the kernel support of particle `i`.
///
/// For wall/boundary handling, see [`isph_step_with_walls`].
pub fn isph_step(
    positions: &mut [[f64; 3]],
    velocities: &mut [[f64; 3]],
    masses: &[f64],
    densities: &mut [f64],
    config: &IsphConfig,
    neighbor_fn: &dyn Fn(usize, f64) -> Vec<usize>,
) -> Result<IsphResult, IsphError> {
    let no_walls = vec![false; positions.len()];
    isph_step_with_walls(
        positions,
        velocities,
        masses,
        densities,
        &no_walls,
        config,
        neighbor_fn,
    )
}

/// Advance the fluid one ISPH step with explicit wall particles.
///
/// Particles flagged in `is_wall` are held fixed: they keep their velocity
/// (typically zero), are not advected, and their density is clamped to `ρ₀`.
/// They still participate in the PPE so they exert a pressure reaction on the
/// fluid (the natural homogeneous-Neumann wall treatment).  Wall particles are
/// never classified as free surfaces.
pub fn isph_step_with_walls(
    positions: &mut [[f64; 3]],
    velocities: &mut [[f64; 3]],
    masses: &[f64],
    densities: &mut [f64],
    is_wall: &[bool],
    config: &IsphConfig,
    neighbor_fn: &dyn Fn(usize, f64) -> Vec<usize>,
) -> Result<IsphResult, IsphError> {
    let n = positions.len();
    if n == 0 {
        return Err(IsphError::NoParticles);
    }
    if config.smoothing_length <= 0.0 || config.dt <= 0.0 || config.rho0 <= 0.0 {
        return Err(IsphError::InvalidConfig);
    }

    let h = config.smoothing_length;
    let dt = config.dt;
    let rho0 = config.rho0;
    let eta = 0.01 * h;
    let eta2 = eta * eta;
    let kernel = CubicSplineKernel;

    let wall = |i: usize| is_wall.get(i).copied().unwrap_or(false);

    // Collect neighbor lists at the start-of-step configuration.
    let neighbors: Vec<Vec<usize>> = (0..n).map(|i| neighbor_fn(i, h)).collect();

    // -- Step 0: density by SPH summation (walls clamped to ρ₀). -------------
    let w_self = kernel.w(0.0, h);
    for i in 0..n {
        if wall(i) {
            densities[i] = rho0;
            continue;
        }
        let mut rho = masses[i] * w_self;
        for &j in &neighbors[i] {
            if j == i {
                continue;
            }
            let rij = sub3(positions[i], positions[j]);
            let r = dot3(rij, rij).sqrt();
            rho += masses[j] * kernel.w(r, h);
        }
        densities[i] = rho.max(1e-6);
    }

    // -- Step 1: predicted velocity v* (gravity + Morris viscosity). --------
    let mut vstar = velocities.to_vec();
    for i in 0..n {
        if wall(i) {
            continue;
        }
        let mut visc = [0.0_f64; 3];
        if config.viscosity != 0.0 {
            for &j in &neighbors[i] {
                if j == i {
                    continue;
                }
                let rij = sub3(positions[i], positions[j]);
                let r2 = dot3(rij, rij);
                if r2 < 1e-24 {
                    continue;
                }
                let gw = grad(&kernel, rij, h);
                let kappa = dot3(rij, gw) / (r2 + eta2);
                let vij = sub3(velocities[i], velocities[j]);
                let f = (masses[j] / densities[j]) * 2.0 * config.viscosity * kappa;
                visc[0] += f * vij[0];
                visc[1] += f * vij[1];
                visc[2] += f * vij[2];
            }
        }
        for d in 0..3 {
            vstar[i][d] = velocities[i][d] + dt * (config.gravity[d] + visc[d]);
        }
    }

    // Shared SPH field context over the frozen start-of-step configuration.
    let field = SphField {
        positions: &*positions,
        masses,
        densities: &*densities,
        neighbors: &neighbors,
        kernel,
        h,
    };

    // -- Step 2: divergence of v*. ------------------------------------------
    let div_star: Vec<f64> = (0..n).map(|i| field.divergence(i, &vstar)).collect();

    // -- Step 3: free-surface (Dirichlet) classification. -------------------
    let mut dirichlet = vec![false; n];
    let mut any_dirichlet = false;
    for i in 0..n {
        if !wall(i) && densities[i] < FREE_SURFACE_FRACTION * rho0 {
            dirichlet[i] = true;
            any_dirichlet = true;
        }
    }
    if !any_dirichlet {
        // Pure-Neumann problem is singular: pin the lowest-density fluid
        // particle (or particle 0) to anchor the pressure level.
        let mut pin = 0usize;
        let mut min_rho = f64::MAX;
        for (i, &rho_i) in densities.iter().enumerate() {
            if !wall(i) && rho_i < min_rho {
                min_rho = rho_i;
                pin = i;
            }
        }
        dirichlet[pin] = true;
    }

    // Interior fluid mask = fluid particles that carry a real PPE unknown.
    let interior: Vec<bool> = (0..n).map(|i| !wall(i) && !dirichlet[i]).collect();
    let predicted_divergence_rms = masked_rms(&div_star, &interior);

    // -- Step 4: assemble the well-conditioned Brookshaw Laplacian once. -----
    // Positions are frozen during the step, so the matrix is reused across the
    // outer projection iterations.
    let a = field.brookshaw(&dirichlet, eta2);

    // The running corrected velocity starts at the predicted velocity v*.
    velocities.copy_from_slice(&vstar);

    // -- Steps 5–6: iterated pressure projection. ---------------------------
    // Each outer sweep computes the velocity-divergence residual, solves
    // `A Δp = b`, and applies the symmetric pressure-gradient correction.
    // Because the Brookshaw Laplacian only approximates the divergence of that
    // gradient, one solve leaves a residual; sweeping drives it toward the
    // projection floor (a preconditioned Richardson iteration with `A`
    // preconditioning the true divergence-of-gradient operator).
    let mut pressures = vec![0.0_f64; n];
    let mut iters = 0usize;
    let mut residual = 0.0_f64;
    let max_outer = 25usize;
    let mut prev_metric = f64::MAX;
    let mut target = 0.0_f64;
    for outer in 0..max_outer {
        // Velocity-divergence residual (to be driven to zero) and PPE RHS.
        let mut resid = vec![0.0_f64; n];
        let mut rhs = vec![0.0_f64; n];
        for i in 0..n {
            let d = field.divergence(i, velocities);
            resid[i] = d;
            if !dirichlet[i] {
                rhs[i] = -d / dt;
            }
        }

        let metric = masked_rms(&resid, &interior);
        if outer == 0 {
            target = (1.0e-2 * metric).max(1.0e-9);
        }
        if metric <= target {
            break;
        }
        if outer > 0 && metric >= prev_metric {
            // Projection floor reached (no further improvement).
            break;
        }
        prev_metric = metric;

        let mut dp = vec![0.0_f64; n];
        let (it, res) = pcg_solve(&a, &rhs, &mut dp, config.pcg_tol, config.pcg_max_iters);
        iters = it;
        residual = res;
        if !res.is_finite() || dp.iter().any(|p| !p.is_finite()) {
            return Err(IsphError::PcgNonConvergence(it, res));
        }

        for i in 0..n {
            if !wall(i) {
                let accel = field.pressure_accel(i, &dp);
                for d in 0..3 {
                    velocities[i][d] += dt * accel[d];
                }
            }
            pressures[i] += dp[i];
        }
    }

    // -- Post-projection divergence (frozen positions, corrected velocities). -
    let div_new: Vec<f64> = (0..n).map(|i| field.divergence(i, velocities)).collect();
    let divergence_rms = masked_rms(&div_new, &interior);

    // -- Step 7: integrate positions (fluid only). --------------------------
    for i in 0..n {
        if wall(i) {
            continue;
        }
        for d in 0..3 {
            positions[i][d] += dt * velocities[i][d];
        }
    }

    Ok(IsphResult {
        pressures,
        pcg_iters: iters,
        pcg_residual: residual,
        divergence_rms,
        predicted_divergence_rms,
    })
}

// ===========================================================================
// Tests
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::neighbor::SpatialHash3D;

    /// Build a per-step neighbor closure from a position snapshot (support 2h).
    fn neighbor_closure(snapshot: &[[f64; 3]], h: f64) -> Vec<Vec<usize>> {
        SpatialHash3D::all_neighbors(snapshot, 2.0 * h)
    }

    // ---- CSR / PCG unit tests ---------------------------------------------

    #[test]
    fn csr_matvec_identity() {
        // 3x3 identity: off-diagonals empty, diagonal all 1.
        let off: Vec<Vec<(usize, f64)>> = vec![Vec::new(), Vec::new(), Vec::new()];
        let diag = vec![1.0, 1.0, 1.0];
        let a = CsrMatrix::from_rows(3, &off, &diag);
        assert_eq!(a.nrows(), 3);
        let x = [2.0, -3.0, 5.0];
        let mut y = [0.0; 3];
        a.matvec(&x, &mut y);
        assert!((y[0] - 2.0).abs() < 1e-12);
        assert!((y[1] + 3.0).abs() < 1e-12);
        assert!((y[2] - 5.0).abs() < 1e-12);
    }

    #[test]
    fn csr_diag_extraction() {
        let off = vec![vec![(1usize, -1.0)], vec![(0usize, -1.0)]];
        let diag = vec![2.0, 3.0];
        let a = CsrMatrix::from_rows(2, &off, &diag);
        let d = a.diag();
        assert!((d[0] - 2.0).abs() < 1e-12);
        assert!((d[1] - 3.0).abs() < 1e-12);
    }

    #[test]
    fn pcg_solves_spd_tridiagonal() {
        // 1D Poisson: tridiag(-1, 2, -1), n = 5, b = ones.
        let n = 5;
        let mut off: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
        let diag = vec![2.0; n];
        for (i, row) in off.iter_mut().enumerate() {
            if i > 0 {
                row.push((i - 1, -1.0));
            }
            if i + 1 < n {
                row.push((i + 1, -1.0));
            }
        }
        let a = CsrMatrix::from_rows(n, &off, &diag);
        let b = vec![1.0; n];
        let mut x = vec![0.0; n];
        let (iters, res) = pcg_solve(&a, &b, &mut x, 1e-10, 100);
        assert!(res < 1e-10, "residual {res} too large after {iters} iters");
        // Verify A x ≈ b.
        let mut ax = vec![0.0; n];
        a.matvec(&x, &mut ax);
        for i in 0..n {
            assert!((ax[i] - b[i]).abs() < 1e-8, "row {i}: Ax={}, b=1", ax[i]);
        }
    }

    #[test]
    fn pcg_zero_rhs_returns_zero() {
        let off = vec![vec![(1usize, -1.0)], vec![(0usize, -1.0)]];
        let diag = vec![2.0, 2.0];
        let a = CsrMatrix::from_rows(2, &off, &diag);
        let b = vec![0.0, 0.0];
        let mut x = vec![7.0, -3.0];
        let (iters, res) = pcg_solve(&a, &b, &mut x, 1e-8, 50);
        assert_eq!(iters, 0);
        assert_eq!(res, 0.0);
        assert_eq!(x, vec![0.0, 0.0]);
    }

    // ---- Error / config handling ------------------------------------------

    #[test]
    fn empty_particle_list_errors() {
        let mut pos: Vec<[f64; 3]> = Vec::new();
        let mut vel: Vec<[f64; 3]> = Vec::new();
        let masses: Vec<f64> = Vec::new();
        let mut dens: Vec<f64> = Vec::new();
        let cfg = IsphConfig::default();
        let nfn = |_i: usize, _h: f64| Vec::new();
        let r = isph_step(&mut pos, &mut vel, &masses, &mut dens, &cfg, &nfn);
        assert!(matches!(r, Err(IsphError::NoParticles)));
    }

    #[test]
    fn invalid_config_errors() {
        let mut pos = vec![[0.0, 0.0, 0.0]];
        let mut vel = vec![[0.0, 0.0, 0.0]];
        let masses = vec![1.0];
        let mut dens = vec![1000.0];
        let cfg = IsphConfig {
            smoothing_length: 0.0,
            ..IsphConfig::default()
        };
        let nfn = |_i: usize, _h: f64| Vec::new();
        let r = isph_step(&mut pos, &mut vel, &masses, &mut dens, &cfg, &nfn);
        assert!(matches!(r, Err(IsphError::InvalidConfig)));
    }

    // ---- Test 1: hydrostatic column ---------------------------------------

    /// Build a fluid block of `nx × ny × nz` particles fully enclosed by
    /// `wall`-thick solid layers on all six faces.  The fluid occupies index
    /// box `[0,nx) × [0,ny) × [0,nz)`; everything outside (down to `-wall` and
    /// up to `n+wall`) is a fixed wall particle.  Enclosure prevents the
    /// divergence-free rigid drift that an open free-surface column would
    /// permit, so the projection develops the full hydrostatic pressure.
    fn build_enclosed_box(
        nx: usize,
        ny: usize,
        nz: usize,
        wall: usize,
        dx: f64,
    ) -> (Vec<[f64; 3]>, Vec<bool>) {
        let mut pos = Vec::new();
        let mut is_wall = Vec::new();
        let w = wall as isize;
        let (nxi, nyi, nzi) = (nx as isize, ny as isize, nz as isize);
        for i in -w..nxi + w {
            for j in -w..nyi + w {
                for k in -w..nzi + w {
                    let solid = i < 0 || i >= nxi || j < 0 || j >= nyi || k < 0 || k >= nzi;
                    pos.push([i as f64 * dx, j as f64 * dx, k as f64 * dx]);
                    is_wall.push(solid);
                }
            }
        }
        (pos, is_wall)
    }

    #[test]
    fn hydrostatic_column_pressure_profile() {
        // Fluid 2x2x11 (height H = 10·dx = 1.0) enclosed by 2 wall layers on
        // every face.  mass = rho0·dx³ = 1.0 (consistent with the spec's
        // "masses all 1.0, rho0 = 1000").  The hydrostatic pressure *difference*
        // between the centre-bottom and centre-top fluid particles must match
        // ρ₀·g·H within 20% (a pin-independent check of the hydrostatic gradient
        // p(z) ≈ ρ₀ g (H − z)).
        let dx = 0.1;
        let (nx, ny, nz) = (2usize, 2usize, 11usize);
        let wall = 2usize;
        let (mut pos, is_wall) = build_enclosed_box(nx, ny, nz, wall, dx);
        let n = pos.len();
        let rho0 = 1000.0;
        let g = 9.81;
        let mass = rho0 * dx * dx * dx; // = 1.0
        let masses = vec![mass; n];
        let mut vel = vec![[0.0, 0.0, 0.0]; n];
        let mut dens = vec![rho0; n];

        let cfg = IsphConfig {
            dt: 1.0e-3,
            smoothing_length: 1.3 * dx,
            rho0,
            viscosity: 0.0,
            gravity: [0.0, 0.0, -g],
            pcg_tol: 1.0e-6,
            pcg_max_iters: 200,
        };

        let mut last = None;
        for _ in 0..10 {
            let snapshot = pos.clone();
            let nls = neighbor_closure(&snapshot, cfg.smoothing_length);
            let nfn = |i: usize, _h: f64| nls[i].clone();
            last = Some(
                isph_step_with_walls(
                    &mut pos, &mut vel, &masses, &mut dens, &is_wall, &cfg, &nfn,
                )
                .expect("hydrostatic step"),
            );
        }
        let result = last.expect("ran 10 steps");

        // Locate the centre fluid column (x = y = dx) at the bottom (z = 0) and
        // top (z = (nz-1)·dx).
        let cx = (nx / 2) as f64 * dx;
        let cy = (ny / 2) as f64 * dx;
        let z_top = (nz - 1) as f64 * dx;
        let nearest_fluid = |target: [f64; 3]| -> usize {
            let mut best = f64::MAX;
            let mut idx = 0usize;
            for i in 0..n {
                if is_wall[i] {
                    continue;
                }
                let d = (pos[i][0] - target[0]).powi(2)
                    + (pos[i][1] - target[1]).powi(2)
                    + (pos[i][2] - target[2]).powi(2);
                if d < best {
                    best = d;
                    idx = i;
                }
            }
            idx
        };
        let bottom = nearest_fluid([cx, cy, 0.0]);
        let top = nearest_fluid([cx, cy, z_top]);
        let p_bottom = result.pressures[bottom];
        let p_top = result.pressures[top];
        let dp = p_bottom - p_top;

        let height = z_top; // 1.0
        let expected = rho0 * g * height; // ρ₀ g H
        let rel_err = (dp - expected).abs() / expected;
        assert!(
            dp.is_finite() && p_bottom.is_finite() && p_top.is_finite(),
            "pressures must be finite"
        );
        assert!(
            dp > 0.0,
            "pressure must increase with depth: p_bottom={p_bottom:.1}, p_top={p_top:.1}"
        );
        assert!(
            rel_err < 0.20,
            "hydrostatic Δp {dp:.1} Pa vs expected {expected:.1} Pa (rel err {:.1}%)",
            rel_err * 100.0
        );
    }

    // ---- Test 2: divergence reduction -------------------------------------

    #[test]
    fn projection_reduces_divergence() {
        // Regular 6x6x6 lattice with a smooth compressive velocity field.
        let dx = 0.1;
        let nside = 6;
        let mut pos = Vec::new();
        for i in 0..nside {
            for j in 0..nside {
                for k in 0..nside {
                    pos.push([i as f64 * dx, j as f64 * dx, k as f64 * dx]);
                }
            }
        }
        let n = pos.len();
        let rho0 = 1000.0;
        let mass = rho0 * dx * dx * dx;
        let masses = vec![mass; n];
        let mut dens = vec![rho0; n];

        // Uniform radial compression toward the centre: ∇·v = const < 0.
        let center = [
            (nside - 1) as f64 * dx * 0.5,
            (nside - 1) as f64 * dx * 0.5,
            (nside - 1) as f64 * dx * 0.5,
        ];
        let beta = 0.5;
        let mut vel = vec![[0.0; 3]; n];
        for i in 0..n {
            for d in 0..3 {
                vel[i][d] = -beta * (pos[i][d] - center[d]);
            }
        }

        let cfg = IsphConfig {
            dt: 1.0e-3,
            smoothing_length: 1.3 * dx,
            rho0,
            viscosity: 0.0,
            gravity: [0.0, 0.0, 0.0], // isolate the projection
            pcg_tol: 1.0e-6,
            pcg_max_iters: 200,
        };

        let snapshot = pos.clone();
        let nls = neighbor_closure(&snapshot, cfg.smoothing_length);
        let nfn = |i: usize, _h: f64| nls[i].clone();
        let result = isph_step(&mut pos, &mut vel, &masses, &mut dens, &cfg, &nfn)
            .expect("divergence-reduction step");

        assert!(
            result.predicted_divergence_rms > 1e-3,
            "predicted divergence should be non-trivial: {}",
            result.predicted_divergence_rms
        );
        assert!(
            result.divergence_rms < 0.1 * result.predicted_divergence_rms,
            "projection must cut divergence by >=90%: before={:.4e}, after={:.4e}",
            result.predicted_divergence_rms,
            result.divergence_rms
        );
    }

    // ---- Test 3: PCG convergence ------------------------------------------

    #[test]
    fn pcg_converges_within_iteration_budget() {
        // ~50-particle blob with a divergent velocity field so the PPE is
        // non-trivial; PCG must converge to < 1e-4 in <= 50 iterations.
        let dx = 0.1;
        let mut pos = Vec::new();
        // 4x4x3 = 48 plus 2 extra → 50.
        for i in 0..4 {
            for j in 0..4 {
                for k in 0..3 {
                    pos.push([i as f64 * dx, j as f64 * dx, k as f64 * dx]);
                }
            }
        }
        pos.push([0.5 * dx, 0.5 * dx, 3.0 * dx]);
        pos.push([1.5 * dx, 1.5 * dx, 3.0 * dx]);
        let n = pos.len();
        assert_eq!(n, 50);

        let rho0 = 1000.0;
        let mass = rho0 * dx * dx * dx;
        let masses = vec![mass; n];
        let mut dens = vec![rho0; n];
        let mut vel = vec![[0.0; 3]; n];
        // Shear-like divergent field.
        for i in 0..n {
            vel[i][0] = 0.3 * pos[i][0];
            vel[i][2] = -0.2 * pos[i][2];
        }

        let cfg = IsphConfig {
            dt: 1.0e-3,
            smoothing_length: 1.3 * dx,
            rho0,
            viscosity: 0.0,
            gravity: [0.0, 0.0, 0.0],
            pcg_tol: 1.0e-4,
            pcg_max_iters: 50,
        };

        let snapshot = pos.clone();
        let nls = neighbor_closure(&snapshot, cfg.smoothing_length);
        let nfn = |i: usize, _h: f64| nls[i].clone();
        let result =
            isph_step(&mut pos, &mut vel, &masses, &mut dens, &cfg, &nfn).expect("pcg step");

        assert!(
            result.pcg_iters <= 50,
            "PCG used {} iterations (budget 50)",
            result.pcg_iters
        );
        assert!(
            result.pcg_residual < 1e-4,
            "PCG residual {:.3e} did not reach tolerance",
            result.pcg_residual
        );
    }

    // ---- Test 4: energy monotonicity / stability --------------------------

    #[test]
    fn dam_break_remains_stable() {
        // 200 fluid particles in a box (floor + side walls) under gravity.
        // Loose bound: kinetic energy must stay finite and bounded.
        let dx = 0.1;
        let rho0 = 1000.0;
        let g = 9.81;
        let mass = rho0 * dx * dx * dx;

        let mut pos = Vec::new();
        let mut is_wall = Vec::new();

        // Floor: 10x2 wall layer at z=-0.1 and z=-0.2.
        for layer in 1..=2 {
            let z = -(layer as f64) * dx;
            for i in 0..10 {
                for j in 0..2 {
                    pos.push([i as f64 * dx, j as f64 * dx, z]);
                    is_wall.push(true);
                }
            }
        }
        // Fluid block: 5 (x) x 2 (y) x 20 (z) = 200 particles.
        for k in 0..20 {
            for i in 0..5 {
                for j in 0..2 {
                    pos.push([i as f64 * dx, j as f64 * dx, k as f64 * dx]);
                    is_wall.push(false);
                }
            }
        }
        let n = pos.len();
        let n_fluid = is_wall.iter().filter(|&&w| !w).count();
        assert_eq!(n_fluid, 200);

        let masses = vec![mass; n];
        let mut vel = vec![[0.0; 3]; n];
        let mut dens = vec![rho0; n];

        // Initial potential energy relative to the floor (z = 0 datum).
        let mut initial_pe = 0.0;
        for i in 0..n {
            if !is_wall[i] {
                initial_pe += mass * g * (pos[i][2] + 1.0);
            }
        }

        let cfg = IsphConfig {
            dt: 5.0e-4,
            smoothing_length: 1.3 * dx,
            rho0,
            viscosity: 1.0e-4,
            gravity: [0.0, 0.0, -g],
            pcg_tol: 1.0e-5,
            pcg_max_iters: 100,
        };

        let mut max_ke = 0.0_f64;
        for _ in 0..20 {
            let snapshot = pos.clone();
            let nls = neighbor_closure(&snapshot, cfg.smoothing_length);
            let nfn = |i: usize, _h: f64| nls[i].clone();
            isph_step_with_walls(&mut pos, &mut vel, &masses, &mut dens, &is_wall, &cfg, &nfn)
                .expect("dam-break step");

            let mut ke = 0.0;
            for i in 0..n {
                if !is_wall[i] {
                    let v2 = vel[i][0].powi(2) + vel[i][1].powi(2) + vel[i][2].powi(2);
                    ke += 0.5 * mass * v2;
                }
            }
            assert!(ke.is_finite(), "kinetic energy became non-finite");
            for p in &pos {
                assert!(p[0].is_finite() && p[1].is_finite() && p[2].is_finite());
            }
            if ke > max_ke {
                max_ke = ke;
            }
        }

        assert!(
            max_ke < 10.0 * initial_pe,
            "kinetic energy {max_ke:.3e} exceeded 10x initial PE {initial_pe:.3e} (instability)"
        );
    }

    // ---- Test 5: Taylor–Green-style divergence decay (thin slab) ----------

    #[test]
    fn taylor_green_divergence_decays() {
        // A thin 3D slab (16x16x3) seeded with a 2D Taylor–Green velocity
        // field in the x-y plane.  The incompressible projection should keep
        // the post-projection divergence well below the predicted divergence.
        let dx = 1.0 / 12.0;
        let nxy = 12;
        let nz = 3;
        let mut pos = Vec::new();
        for i in 0..nxy {
            for j in 0..nxy {
                for k in 0..nz {
                    pos.push([i as f64 * dx, j as f64 * dx, k as f64 * dx]);
                }
            }
        }
        let n = pos.len();
        let rho0 = 1000.0;
        let mass = rho0 * dx * dx * dx;
        let masses = vec![mass; n];
        let mut dens = vec![rho0; n];

        // Taylor–Green: u = sin(2πx)cos(2πy), v = −cos(2πx)sin(2πy), w = 0.
        // (Analytically divergence-free, but the discrete SPH operator and the
        // finite slab introduce a non-zero predicted divergence.)
        use std::f64::consts::PI;
        let mut vel = vec![[0.0; 3]; n];
        for i in 0..n {
            let x = pos[i][0];
            let y = pos[i][1];
            vel[i][0] = (2.0 * PI * x).sin() * (2.0 * PI * y).cos();
            vel[i][1] = -(2.0 * PI * x).cos() * (2.0 * PI * y).sin();
        }

        let cfg = IsphConfig {
            dt: 1.0e-3,
            smoothing_length: 1.3 * dx,
            rho0,
            viscosity: 0.0,
            gravity: [0.0, 0.0, 0.0],
            pcg_tol: 1.0e-6,
            pcg_max_iters: 200,
        };

        let mut result = None;
        for _ in 0..5 {
            let snapshot = pos.clone();
            let nls = neighbor_closure(&snapshot, cfg.smoothing_length);
            let nfn = |i: usize, _h: f64| nls[i].clone();
            result = Some(
                isph_step(&mut pos, &mut vel, &masses, &mut dens, &cfg, &nfn).expect("tg step"),
            );
        }
        let result = result.expect("ran TG steps");
        assert!(result.divergence_rms.is_finite());
        assert!(
            result.divergence_rms < result.predicted_divergence_rms,
            "projection should reduce divergence: before={:.4e}, after={:.4e}",
            result.predicted_divergence_rms,
            result.divergence_rms
        );
    }
}
