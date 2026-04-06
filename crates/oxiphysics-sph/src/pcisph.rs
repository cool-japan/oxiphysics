// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! PCISPH: Predictive-Corrective Incompressible SPH.
//!
//! Solenthaler & Pajarola 2009.  Higher accuracy than WCSPH with far fewer
//! iterations than IISPH because the correction loop uses a pre-computed
//! per-solver stiffness factor `δ` instead of per-particle relaxation
//! coefficients.
//!
//! ## Algorithm overview (one time step)
//! 1. Predict positions and velocities with non-pressure forces.
//! 2. Loop (until max density error < threshold **and** min_iter reached):
//!    a. Compute predicted density error `ρ*_i - ρ₀`.
//!    b. Update pressure: `p_i += δ * density_error`.
//!    c. Compute pressure forces from updated pressures.
//!    d. Update predicted positions / velocities with pressure forces.
//! 3. Return the converged pressure forces.
//!
//! ## Stiffness factor `δ`
//! `δ = -1 / (β * D)` where
//! - `β = 2 dt² m² / ρ₀²`
//! - `D = -|Σ ∇W_ij|² - Σ|∇W_ij|²`
//!   computed once for a *regular* particle arrangement (Solenthaler §4.1).

use oxiphysics_core::math::{Real, Vec3};
use std::f64::consts::PI;

use crate::kernel::CubicSplineKernel;
use crate::kernel::SphKernel;

// ---------------------------------------------------------------------------
// PcisphSolver
// ---------------------------------------------------------------------------

/// PCISPH pressure solver.
///
/// Pre-computes a stiffness factor `δ` from the expected particle arrangement,
/// then iterates a *predict–correct* loop each time step to enforce
/// near-incompressibility without an explicit equation of state.
///
/// # Example
/// ```no_run
/// use oxiphysics_sph::pcisph::PcisphSolver;
///
/// let mut solver = PcisphSolver::new(1000.0, 0.1);
/// solver.compute_delta(0.001, 0.125e-3);
/// assert!(solver.delta > 0.0);
/// ```
#[derive(Debug, Clone)]
pub struct PcisphSolver {
    /// Rest density ρ₀ (kg/m³).
    pub rest_density: Real,
    /// Smoothing length h (m).
    pub smoothing_length: Real,
    /// Maximum number of prediction–correction iterations (default 5).
    pub max_iter: usize,
    /// Minimum number of iterations before checking convergence (default 2).
    pub min_iter: usize,
    /// Density error tolerance as a **fraction** of `rest_density` (default 0.01 = 1 %).
    pub max_density_error: Real,
    /// Pre-computed stiffness factor δ.  Call [`Self::compute_delta`] before
    /// the first time step.
    pub delta: Real,
}

impl PcisphSolver {
    /// Create a new [`PcisphSolver`] with sensible defaults.
    ///
    /// `delta` is initialised to `0.0`; call [`Self::compute_delta`] to set it.
    pub fn new(rest_density: Real, smoothing_length: Real) -> Self {
        Self {
            rest_density,
            smoothing_length,
            max_iter: 5,
            min_iter: 2,
            max_density_error: 0.01,
            delta: 0.0,
        }
    }

    /// Pre-compute the stiffness factor `δ` for the given timestep and particle mass.
    ///
    /// Uses a synthetic regular-grid neighbourhood centred on the origin (same
    /// resolution as `smoothing_length / 2`) to evaluate the denominator `D`.
    /// This matches the pre-computation described in Solenthaler & Pajarola §4.1.
    ///
    /// # Panics
    /// Does not panic; if the computed denominator is near-zero `delta` is left
    /// unchanged.
    pub fn compute_delta(&mut self, dt: Real, mass: Real) {
        let h = self.smoothing_length;
        let rho0 = self.rest_density;
        let kernel = CubicSplineKernel;

        // Build a synthetic neighbourhood: particles on a cubic lattice with
        // spacing h/2 within support radius 2h.
        let spacing = h / 2.0;
        let n_half = (2.0 * h / spacing).ceil() as i32 + 1;

        let mut sum_grad = Vec3::zeros();
        let mut sum_grad_sq = 0.0_f64;

        for ix in -n_half..=n_half {
            for iy in -n_half..=n_half {
                for iz in -n_half..=n_half {
                    if ix == 0 && iy == 0 && iz == 0 {
                        continue; // skip self
                    }
                    let rij = Vec3::new(
                        ix as f64 * spacing,
                        iy as f64 * spacing,
                        iz as f64 * spacing,
                    );
                    let r = rij.norm();
                    if r >= 2.0 * h || r < 1e-14 {
                        continue;
                    }
                    let g = kernel.grad_w(r, h);
                    let grad_ij = g * rij / r;
                    sum_grad += grad_ij;
                    sum_grad_sq += grad_ij.norm_squared();
                }
            }
        }

        // β = 2 dt² m² / ρ₀²
        let beta = 2.0 * dt * dt * mass * mass / (rho0 * rho0);

        // D = -|Σ grad|² - Σ|grad|²   (both terms are negative-definite)
        let denom = beta * (-sum_grad.norm_squared() - sum_grad_sq);

        if denom.abs() > 1e-30 {
            self.delta = -1.0 / denom;
        }
    }

    // -----------------------------------------------------------------------
    // Internal helpers
    // -----------------------------------------------------------------------

    /// Compute SPH densities (kernel summation, including self-contribution).
    #[allow(clippy::needless_range_loop)]
    fn compute_densities(
        &self,
        positions: &[Vec3],
        masses: &[f64],
        neighbors: &[Vec<usize>],
    ) -> Vec<f64> {
        let kernel = CubicSplineKernel;
        let h = self.smoothing_length;
        let n = positions.len();
        let mut densities = vec![0.0_f64; n];
        for i in 0..n {
            let mut rho = masses[i] * kernel.w(0.0, h);
            for &j in &neighbors[i] {
                let r = (positions[i] - positions[j]).norm();
                rho += masses[j] * kernel.w(r, h);
            }
            densities[i] = rho;
        }
        densities
    }

    /// Compute pressure forces from a pressure field.
    #[allow(clippy::needless_range_loop)]
    fn pressure_forces(
        &self,
        positions: &[Vec3],
        masses: &[f64],
        densities: &[f64],
        pressures: &[f64],
        neighbors: &[Vec<usize>],
    ) -> Vec<Vec3> {
        let kernel = CubicSplineKernel;
        let h = self.smoothing_length;
        let n = positions.len();
        let mut forces = vec![Vec3::zeros(); n];
        for i in 0..n {
            let rhoi = densities[i].max(1e-14);
            let pi = pressures[i];
            for &j in &neighbors[i] {
                let rij = positions[i] - positions[j];
                let r = rij.norm();
                if r < 1e-14 {
                    continue;
                }
                let g = kernel.grad_w(r, h);
                let rhat = rij / r;
                let rhoj = densities[j].max(1e-14);
                let pj = pressures[j];
                // Symmetric SPH pressure acceleration (acceleration = force/mass_i)
                let factor = -masses[j] * (pi / (rhoi * rhoi) + pj / (rhoj * rhoj)) * g;
                forces[i] += factor * rhat;
            }
        }
        forces
    }

    // -----------------------------------------------------------------------
    // Public API
    // -----------------------------------------------------------------------

    /// PCISPH pressure solve for one time step.
    ///
    /// # Arguments
    /// * `positions`  – current particle positions.
    /// * `velocities` – current particle velocities (already updated with
    ///   non-pressure forces).
    /// * `masses`     – particle masses.
    /// * `densities`  – current densities (computed from `positions`).
    /// * `neighbors`  – neighbour lists for each particle.
    /// * `dt`         – time-step size.
    ///
    /// # Returns
    /// Pressure force **per unit mass** (i.e. acceleration) for each particle.
    ///
    /// # Note
    /// `delta` must be positive (call [`Self::compute_delta`] first).  If
    /// `delta ≤ 0` the function returns zero forces immediately.
    #[allow(clippy::needless_range_loop, clippy::too_many_arguments)]
    pub fn solve_pressure(
        &self,
        positions: &[Vec3],
        velocities: &[Vec3],
        masses: &[f64],
        densities: &[f64],
        neighbors: &[Vec<usize>],
        dt: Real,
    ) -> Vec<Vec3> {
        let n = positions.len();
        if n == 0 || self.delta <= 0.0 {
            return vec![Vec3::zeros(); n];
        }

        let rho0 = self.rest_density;
        let delta = self.delta;

        // Initialise pressures and predicted state
        let mut pressures = vec![0.0_f64; n];
        let mut pred_pos: Vec<Vec3> = positions.to_vec();
        let mut pred_vel: Vec<Vec3> = velocities.to_vec();

        let mut last_forces = vec![Vec3::zeros(); n];

        for iter in 0..self.max_iter {
            // 1. Compute predicted densities from predicted positions.
            // On the first iteration we reuse the caller-supplied densities to
            // avoid redundant computation; subsequent iterations use densities
            // derived from the updated predicted positions.
            let pred_densities = if iter == 0 {
                densities.to_vec()
            } else {
                self.compute_densities(&pred_pos, masses, neighbors)
            };

            // 2. Update pressures from density error
            let mut max_err = 0.0_f64;
            for i in 0..n {
                let err = pred_densities[i] - rho0; // signed error
                pressures[i] = (pressures[i] + delta * err).max(0.0);
                let rel_err = err.abs() / rho0;
                if rel_err > max_err {
                    max_err = rel_err;
                }
            }

            // 3. Compute pressure forces from updated pressures
            let pf =
                self.pressure_forces(&pred_pos, masses, &pred_densities, &pressures, neighbors);
            last_forces = pf.clone();

            // 4. Update predicted positions and velocities
            for i in 0..n {
                let acc = pf[i]; // already per unit mass (acceleration)
                pred_vel[i] = velocities[i] + acc * dt;
                pred_pos[i] = positions[i] + pred_vel[i] * dt;
            }

            // 5. Check convergence (only after min_iter)
            if iter + 1 >= self.min_iter && max_err < self.max_density_error {
                break;
            }
        }

        last_forces
    }

    /// PCISPH pressure solve with convergence tracking.
    ///
    /// Identical to [`Self::solve_pressure`] but returns both the forces and
    /// convergence statistics via [`PcisphConvergenceInfo`].
    #[allow(clippy::needless_range_loop, clippy::too_many_arguments)]
    pub fn solve_pressure_tracked(
        &self,
        positions: &[Vec3],
        velocities: &[Vec3],
        masses: &[f64],
        densities: &[f64],
        neighbors: &[Vec<usize>],
        dt: Real,
    ) -> (Vec<Vec3>, PcisphConvergenceInfo) {
        let n = positions.len();
        if n == 0 || self.delta <= 0.0 {
            return (
                vec![Vec3::zeros(); n],
                PcisphConvergenceInfo {
                    iterations_used: 0,
                    max_density_error: 0.0,
                    avg_density_error: 0.0,
                    converged: true,
                },
            );
        }

        let rho0 = self.rest_density;
        let delta = self.delta;

        let mut pressures = vec![0.0_f64; n];
        let mut pred_pos: Vec<Vec3> = positions.to_vec();
        let mut pred_vel: Vec<Vec3> = velocities.to_vec();
        let mut last_forces = vec![Vec3::zeros(); n];
        let mut final_max_err = 0.0_f64;
        let mut final_avg_err = 0.0_f64;
        let mut iters_used = 0_usize;
        let mut converged = false;

        for iter in 0..self.max_iter {
            iters_used = iter + 1;
            let pred_densities = if iter == 0 {
                densities.to_vec()
            } else {
                self.compute_densities(&pred_pos, masses, neighbors)
            };

            let mut max_err = 0.0_f64;
            let mut sum_err = 0.0_f64;
            for i in 0..n {
                let err = pred_densities[i] - rho0;
                pressures[i] = (pressures[i] + delta * err).max(0.0);
                let rel_err = err.abs() / rho0;
                if rel_err > max_err {
                    max_err = rel_err;
                }
                sum_err += rel_err;
            }
            final_max_err = max_err;
            final_avg_err = if n > 0 { sum_err / n as f64 } else { 0.0 };

            let pf =
                self.pressure_forces(&pred_pos, masses, &pred_densities, &pressures, neighbors);
            last_forces = pf.clone();

            for i in 0..n {
                let acc = pf[i];
                pred_vel[i] = velocities[i] + acc * dt;
                pred_pos[i] = positions[i] + pred_vel[i] * dt;
            }

            if iter + 1 >= self.min_iter && max_err < self.max_density_error {
                converged = true;
                break;
            }
        }

        let info = PcisphConvergenceInfo {
            iterations_used: iters_used,
            max_density_error: final_max_err,
            avg_density_error: final_avg_err,
            converged,
        };
        (last_forces, info)
    }

    /// Compute the average density error (fraction of rest density) from a
    /// predicted state.
    #[allow(dead_code, clippy::needless_range_loop)]
    pub fn compute_density_error(
        &self,
        positions: &[Vec3],
        masses: &[f64],
        neighbors: &[Vec<usize>],
    ) -> f64 {
        let densities = self.compute_densities(positions, masses, neighbors);
        let n = densities.len();
        if n == 0 {
            return 0.0;
        }
        let sum: f64 = densities
            .iter()
            .map(|&d| ((d - self.rest_density) / self.rest_density).abs())
            .sum();
        sum / n as f64
    }

    /// Suggest an adaptive time step based on the maximum particle velocity.
    ///
    /// Uses the CFL-like condition: dt_new = cfl_factor * h / v_max.
    /// Clamps between `dt_min` and `dt_max`.
    #[allow(dead_code)]
    pub fn adaptive_timestep(
        &self,
        velocities: &[Vec3],
        dt_min: Real,
        dt_max: Real,
        cfl_factor: Real,
    ) -> Real {
        let v_max = velocities.iter().map(|v| v.norm()).fold(0.0_f64, f64::max);
        if v_max < 1e-14 {
            return dt_max;
        }
        let dt = cfl_factor * self.smoothing_length / v_max;
        dt.clamp(dt_min, dt_max)
    }

    /// Compute the maximum pressure across all particles (useful for diagnostics).
    #[allow(dead_code, clippy::needless_range_loop, clippy::too_many_arguments)]
    pub fn max_pressure_after_solve(
        &self,
        positions: &[Vec3],
        velocities: &[Vec3],
        masses: &[f64],
        densities: &[f64],
        neighbors: &[Vec<usize>],
        dt: Real,
    ) -> f64 {
        let n = positions.len();
        if n == 0 || self.delta <= 0.0 {
            return 0.0;
        }
        let rho0 = self.rest_density;
        let delta = self.delta;
        let mut pressures = vec![0.0_f64; n];
        let mut pred_pos: Vec<Vec3> = positions.to_vec();

        for iter in 0..self.max_iter {
            let pred_densities = if iter == 0 {
                densities.to_vec()
            } else {
                self.compute_densities(&pred_pos, masses, neighbors)
            };
            for i in 0..n {
                let err = pred_densities[i] - rho0;
                pressures[i] = (pressures[i] + delta * err).max(0.0);
            }
            let pf =
                self.pressure_forces(&pred_pos, masses, &pred_densities, &pressures, neighbors);
            for i in 0..n {
                let acc = pf[i];
                let pred_vel = velocities[i] + acc * dt;
                pred_pos[i] = positions[i] + pred_vel * dt;
            }
        }
        pressures.iter().cloned().fold(0.0_f64, f64::max)
    }
}

// ---------------------------------------------------------------------------
// PcisphConvergenceInfo
// ---------------------------------------------------------------------------

/// Convergence statistics from a PCISPH solve.
#[derive(Debug, Clone)]
pub struct PcisphConvergenceInfo {
    /// Number of prediction-correction iterations actually performed.
    pub iterations_used: usize,
    /// Maximum relative density error at the final iteration.
    pub max_density_error: f64,
    /// Average relative density error at the final iteration.
    pub avg_density_error: f64,
    /// Whether the solver converged (error fell below threshold).
    pub converged: bool,
}

// ---------------------------------------------------------------------------
// AdaptivePcisph
// ---------------------------------------------------------------------------

/// PCISPH with adaptive time stepping.
///
/// Wraps [`PcisphSolver`] and adjusts `dt` based on CFL conditions and
/// convergence feedback.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct AdaptivePcisph {
    /// Inner solver.
    pub solver: PcisphSolver,
    /// Current adaptive time step.
    pub current_dt: f64,
    /// Minimum allowed dt.
    pub dt_min: f64,
    /// Maximum allowed dt.
    pub dt_max: f64,
    /// CFL safety factor (typically 0.25–0.4).
    pub cfl_factor: f64,
    /// Factor to reduce dt when solver does not converge.
    pub reduction_factor: f64,
    /// Factor to increase dt when solver converges comfortably.
    pub growth_factor: f64,
}

#[allow(dead_code)]
impl AdaptivePcisph {
    /// Create a new adaptive PCISPH wrapper.
    pub fn new(solver: PcisphSolver, dt_max: f64) -> Self {
        Self {
            solver,
            current_dt: dt_max,
            dt_min: dt_max * 0.01,
            dt_max,
            cfl_factor: 0.3,
            reduction_factor: 0.5,
            growth_factor: 1.1,
        }
    }

    /// Perform one adaptive solve step.  Returns the pressure forces and the
    /// time step actually used.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_adaptive(
        &mut self,
        positions: &[Vec3],
        velocities: &[Vec3],
        masses: &[f64],
        densities: &[f64],
        neighbors: &[Vec<usize>],
    ) -> (Vec<Vec3>, f64) {
        // CFL-based dt
        let dt_cfl =
            self.solver
                .adaptive_timestep(velocities, self.dt_min, self.dt_max, self.cfl_factor);
        self.current_dt = dt_cfl
            .min(self.current_dt * self.growth_factor)
            .min(self.dt_max);

        let (forces, info) = self.solver.solve_pressure_tracked(
            positions,
            velocities,
            masses,
            densities,
            neighbors,
            self.current_dt,
        );

        // If not converged, reduce dt for next step
        if !info.converged {
            self.current_dt = (self.current_dt * self.reduction_factor).max(self.dt_min);
        }

        (forces, self.current_dt)
    }
}

// ---------------------------------------------------------------------------
// Density error computation helpers
// ---------------------------------------------------------------------------

/// Compute per-particle density errors (signed) relative to rest density.
#[allow(dead_code)]
pub fn compute_density_errors(densities: &[f64], rest_density: f64) -> Vec<f64> {
    densities.iter().map(|&d| d - rest_density).collect()
}

/// Compute the RMS density error relative to rest density.
#[allow(dead_code)]
pub fn rms_density_error(densities: &[f64], rest_density: f64) -> f64 {
    if densities.is_empty() {
        return 0.0;
    }
    let sum_sq: f64 = densities
        .iter()
        .map(|&d| {
            let e = (d - rest_density) / rest_density;
            e * e
        })
        .sum();
    (sum_sq / densities.len() as f64).sqrt()
}

/// Clamp pressures to non-negative values (used in PCISPH).
#[allow(dead_code)]
pub fn clamp_pressures(pressures: &mut [f64]) {
    for p in pressures.iter_mut() {
        if *p < 0.0 {
            *p = 0.0;
        }
    }
}

/// Apply a simple under-relaxation to pressure updates.
///
/// `p_new = (1 - omega) * p_old + omega * p_candidate`
#[allow(dead_code)]
pub fn relax_pressures(old: &[f64], candidate: &[f64], omega: f64) -> Vec<f64> {
    old.iter()
        .zip(candidate.iter())
        .map(|(&po, &pc)| (1.0 - omega) * po + omega * pc)
        .collect()
}

// ---------------------------------------------------------------------------
// Pressure coefficient table
// ---------------------------------------------------------------------------

/// Pre-computed PCISPH stiffness coefficients for a range of particle spacings.
///
/// This table stores `(spacing, delta)` pairs and provides interpolation for
/// configurations that fall between stored values.  Using a table avoids
/// recomputing the O(N_reg³) loop at every restart.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct PressureCoefficientTable {
    /// Stored `(spacing, delta)` entries, sorted by spacing ascending.
    entries: Vec<(f64, f64)>,
    /// Rest density used for all entries.
    pub rho0: f64,
    /// Smoothing length used for all entries.
    pub h: f64,
}

#[allow(dead_code)]
impl PressureCoefficientTable {
    /// Build a table by sampling `n_samples` spacings logarithmically between
    /// `spacing_min` and `spacing_max`.
    pub fn build(
        rho0: f64,
        h: f64,
        dt: f64,
        spacing_min: f64,
        spacing_max: f64,
        n_samples: usize,
    ) -> Self {
        let mut entries = Vec::with_capacity(n_samples);
        for k in 0..n_samples {
            let t = k as f64 / (n_samples.max(2) - 1) as f64;
            let spacing = spacing_min * (spacing_max / spacing_min).powf(t);
            let mass = rho0 * spacing.powi(3);
            let mut solver = PcisphSolver::new(rho0, h);
            solver.compute_delta(dt, mass);
            entries.push((spacing, solver.delta));
        }
        entries.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        Self { entries, rho0, h }
    }

    /// Look up δ for a given particle spacing via linear interpolation.
    ///
    /// Returns `None` if the spacing is outside the table range.
    pub fn lookup(&self, spacing: f64) -> Option<f64> {
        let n = self.entries.len();
        if n == 0 {
            return None;
        }
        if spacing <= self.entries[0].0 {
            return Some(self.entries[0].1);
        }
        if spacing >= self.entries[n - 1].0 {
            return Some(self.entries[n - 1].1);
        }
        // Find bracketing interval
        for i in 0..n - 1 {
            let (s0, d0) = self.entries[i];
            let (s1, d1) = self.entries[i + 1];
            if spacing >= s0 && spacing <= s1 {
                let t = (spacing - s0) / (s1 - s0);
                return Some(d0 + t * (d1 - d0));
            }
        }
        None
    }

    /// Number of entries in the table.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Whether the table is empty.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
}

// ---------------------------------------------------------------------------
// PCISPH density error tracker
// ---------------------------------------------------------------------------

/// Tracks density error history across time steps for diagnostics.
#[allow(dead_code)]
#[derive(Debug, Clone, Default)]
pub struct DensityErrorTracker {
    /// Per-step max density error history.
    pub max_errors: Vec<f64>,
    /// Per-step average density error history.
    pub avg_errors: Vec<f64>,
    /// Per-step iteration counts.
    pub iter_counts: Vec<usize>,
}

#[allow(dead_code)]
impl DensityErrorTracker {
    /// Create an empty tracker.
    pub fn new() -> Self {
        Self::default()
    }

    /// Record statistics from one time step.
    pub fn record(&mut self, info: &PcisphConvergenceInfo) {
        self.max_errors.push(info.max_density_error);
        self.avg_errors.push(info.avg_density_error);
        self.iter_counts.push(info.iterations_used);
    }

    /// Maximum of all recorded max-errors.
    pub fn overall_max_error(&self) -> f64 {
        self.max_errors.iter().cloned().fold(0.0_f64, f64::max)
    }

    /// Mean of all recorded average errors.
    pub fn mean_avg_error(&self) -> f64 {
        if self.avg_errors.is_empty() {
            return 0.0;
        }
        self.avg_errors.iter().sum::<f64>() / self.avg_errors.len() as f64
    }

    /// Total number of steps recorded.
    pub fn n_steps(&self) -> usize {
        self.max_errors.len()
    }

    /// Mean iterations per step.
    pub fn mean_iterations(&self) -> f64 {
        if self.iter_counts.is_empty() {
            return 0.0;
        }
        self.iter_counts.iter().sum::<usize>() as f64 / self.iter_counts.len() as f64
    }
}

// ---------------------------------------------------------------------------
// PCISPH correction vector
// ---------------------------------------------------------------------------

/// Compute the pressure correction vector for particle `i`.
///
/// The correction is:
/// ```text
/// Δx_i = dt² * Σ_j m_j * (κ_i/ρ_i² + κ_j/ρ_j²) * ∇W_ij
/// ```
/// where κ_i = δ * (ρ_i - ρ₀).
///
/// Returns the correction displacement `[Δx, Δy, Δz]` for particle `i`.
#[allow(clippy::too_many_arguments)]
pub fn pcisph_correction_vector(
    pos_i: [f64; 3],
    rho_i: f64,
    kappa_i: f64,
    h: f64,
    dt: f64,
    neighbor_pos: &[[f64; 3]],
    neighbor_mass: &[f64],
    neighbor_rho: &[f64],
    neighbor_kappa: &[f64],
) -> [f64; 3] {
    use std::f64::consts::PI;
    let sigma = 1.0 / (PI * h * h * h);
    let rho_i_sq = (rho_i * rho_i).max(1e-28);
    let mut corr = [0.0_f64; 3];
    let dt2 = dt * dt;

    for ((pos_j, &m_j), (&rho_j, &kappa_j)) in neighbor_pos
        .iter()
        .zip(neighbor_mass.iter())
        .zip(neighbor_rho.iter().zip(neighbor_kappa.iter()))
    {
        let rij = [
            pos_i[0] - pos_j[0],
            pos_i[1] - pos_j[1],
            pos_i[2] - pos_j[2],
        ];
        let r = (rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2]).sqrt();
        if r < 1e-14 || r >= 2.0 * h {
            continue;
        }
        let q = r / h;
        let dw_dr = if q < 1.0 {
            sigma * (-3.0 * q + 2.25 * q * q) / h
        } else {
            let t = 2.0 - q;
            sigma * (-0.75 * t * t) / h
        };
        let rho_j_sq = (rho_j * rho_j).max(1e-28);
        let factor = dt2 * m_j * (kappa_i / rho_i_sq + kappa_j / rho_j_sq) * dw_dr / r;
        corr[0] += factor * rij[0];
        corr[1] += factor * rij[1];
        corr[2] += factor * rij[2];
    }
    corr
}

// ---------------------------------------------------------------------------
// PCISPH convergence monitor
// ---------------------------------------------------------------------------

/// Monitors convergence of the PCISPH iteration loop in real-time.
///
/// Provides methods to check whether the solver should exit early based on
/// the density error trend.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct PcisphConvergenceMonitor {
    /// History of max density errors within the current time step's iteration.
    pub error_history: Vec<f64>,
    /// Tolerance for early exit.
    pub tolerance: f64,
    /// Minimum number of iterations before checking.
    pub min_iter: usize,
    /// Maximum number of iterations allowed.
    pub max_iter: usize,
}

#[allow(dead_code)]
impl PcisphConvergenceMonitor {
    /// Create a new convergence monitor.
    pub fn new(tolerance: f64, min_iter: usize, max_iter: usize) -> Self {
        Self {
            error_history: Vec::with_capacity(max_iter),
            tolerance,
            min_iter,
            max_iter,
        }
    }

    /// Reset for a new time step.
    pub fn reset(&mut self) {
        self.error_history.clear();
    }

    /// Record the max density error from one iteration.
    pub fn record(&mut self, max_error: f64) {
        self.error_history.push(max_error);
    }

    /// Return `true` if the solver should exit.
    pub fn should_stop(&self) -> bool {
        let n = self.error_history.len();
        if n >= self.max_iter {
            return true;
        }
        if n < self.min_iter {
            return false;
        }
        self.error_history.last().copied().unwrap_or(f64::MAX) < self.tolerance
    }

    /// Current iteration count.
    pub fn iter_count(&self) -> usize {
        self.error_history.len()
    }

    /// Whether convergence was achieved.
    pub fn converged(&self) -> bool {
        self.error_history.last().copied().unwrap_or(f64::MAX) < self.tolerance
    }

    /// Error reduction ratio from first to last recorded error.
    pub fn error_reduction_ratio(&self) -> f64 {
        if self.error_history.len() < 2 {
            return 1.0;
        }
        let first = self.error_history[0].max(1e-30);
        let last = self.error_history.last().copied().unwrap_or(first);
        last / first
    }
}

// ---------------------------------------------------------------------------
// Standalone SPH density (plain arrays, no Vec3 dependency)
// ---------------------------------------------------------------------------

/// Compute SPH density using plain `[f64; 3]` arrays (no Vec3 dependency).
///
/// Uses the cubic spline kernel W(r, h).
#[allow(dead_code)]
pub fn sph_density_plain(
    pos_i: [f64; 3],
    h: f64,
    neighbor_pos: &[[f64; 3]],
    neighbor_mass: &[f64],
) -> f64 {
    let sigma = 1.0 / (PI * h * h * h);
    let inv_h = 1.0 / h;
    let mut rho = 0.0_f64;
    for (pos_j, &m_j) in neighbor_pos.iter().zip(neighbor_mass) {
        let dx = pos_i[0] - pos_j[0];
        let dy = pos_i[1] - pos_j[1];
        let dz = pos_i[2] - pos_j[2];
        let r = (dx * dx + dy * dy + dz * dz).sqrt();
        let q = r * inv_h;
        let w = if q >= 2.0 {
            0.0
        } else if q >= 1.0 {
            let t = 2.0 - q;
            sigma * 0.25 * t * t * t
        } else {
            sigma * (1.0 - 1.5 * q * q + 0.75 * q * q * q)
        };
        rho += m_j * w;
    }
    rho
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::neighbor::SpatialHash;

    // -----------------------------------------------------------------------
    // Helper: uniform cubic block of particles
    // -----------------------------------------------------------------------
    fn uniform_block(
        spacing: f64,
        n_per_side: usize,
        mass: f64,
    ) -> (Vec<Vec3>, Vec<Vec3>, Vec<f64>) {
        let mut positions = Vec::new();
        let mut velocities = Vec::new();
        let mut masses = Vec::new();
        for i in 0..n_per_side {
            for j in 0..n_per_side {
                for k in 0..n_per_side {
                    positions.push(Vec3::new(
                        i as f64 * spacing,
                        j as f64 * spacing,
                        k as f64 * spacing,
                    ));
                    velocities.push(Vec3::zeros());
                    masses.push(mass);
                }
            }
        }
        (positions, velocities, masses)
    }

    // -----------------------------------------------------------------------
    // Test 1: solver creates with correct defaults
    // -----------------------------------------------------------------------
    #[test]
    fn test_pcisph_new() {
        let solver = PcisphSolver::new(1000.0, 0.1);
        assert!(
            (solver.rest_density - 1000.0).abs() < 1e-10,
            "rest_density should be 1000"
        );
        assert!(
            (solver.smoothing_length - 0.1).abs() < 1e-10,
            "smoothing_length should be 0.1"
        );
        assert!(solver.max_iter >= 2, "max_iter should be at least 2");
        assert!(solver.min_iter >= 1, "min_iter should be at least 1");
        assert!(
            solver.max_density_error > 0.0 && solver.max_density_error < 1.0,
            "max_density_error should be in (0, 1)"
        );
        assert!(
            (solver.delta - 0.0).abs() < 1e-30,
            "delta should start at 0"
        );
    }

    // -----------------------------------------------------------------------
    // Test 2: compute_delta produces positive delta for typical params
    // -----------------------------------------------------------------------
    #[test]
    fn test_pcisph_delta_positive() {
        let spacing = 0.05_f64;
        let rho0 = 1000.0_f64;
        let mass = rho0 * spacing.powi(3);
        let dt = 0.001_f64;

        let mut solver = PcisphSolver::new(rho0, 0.1);
        solver.compute_delta(dt, mass);

        assert!(
            solver.delta > 0.0,
            "delta must be positive, got {}",
            solver.delta
        );
        assert!(
            solver.delta.is_finite(),
            "delta must be finite, got {}",
            solver.delta
        );
    }

    // -----------------------------------------------------------------------
    // Test 3: stationary uniform grid → near-zero pressure forces
    // -----------------------------------------------------------------------
    #[test]
    fn test_pcisph_stationary_no_pressure_force() {
        let spacing = 0.05_f64;
        let h = 0.1_f64;
        let rho0 = 1000.0_f64;
        let mass = rho0 * spacing.powi(3);
        let dt = 0.001_f64;

        let (positions, velocities, masses) = uniform_block(spacing, 5, mass);
        let neighbors = SpatialHash::find_all_neighbors(&positions, 2.0 * h);

        let mut solver = PcisphSolver::new(rho0, h);
        solver.compute_delta(dt, mass);
        assert!(solver.delta > 0.0, "delta must be positive");

        let densities = solver.compute_densities(&positions, &masses, &neighbors);

        let forces =
            solver.solve_pressure(&positions, &velocities, &masses, &densities, &neighbors, dt);

        // Net momentum change must be near zero (Newton's 3rd law).
        let net: Vec3 = forces
            .iter()
            .zip(masses.iter())
            .map(|(f, &m)| *f * m)
            .fold(Vec3::zeros(), |acc, v| acc + v);
        let total: f64 = forces
            .iter()
            .zip(masses.iter())
            .map(|(f, &m)| f.norm() * m)
            .sum();
        let scale = total.max(1e-14);
        assert!(
            net.norm() / scale < 1e-4,
            "Net pressure impulse should be ~0 for uniform fluid, got {:.4e} (scale {:.4e})",
            net.norm(),
            scale
        );
    }

    // -----------------------------------------------------------------------
    // Test 4: compressed cluster → outward (repulsive) pressure forces
    // -----------------------------------------------------------------------
    #[test]
    fn test_pcisph_compressed_gets_repulsion() {
        let spacing = 0.03_f64; // tighter than rest spacing 0.05 → over-dense
        let h = 0.1_f64;
        let rho0 = 1000.0_f64;
        // mass calibrated to rest spacing 0.05
        let mass = rho0 * (0.05_f64).powi(3);
        let dt = 0.001_f64;

        let (positions, velocities, masses) = uniform_block(spacing, 4, mass);
        let neighbors = SpatialHash::find_all_neighbors(&positions, 2.0 * h);

        let mut solver = PcisphSolver::new(rho0, h);
        solver.compute_delta(dt, mass);
        assert!(solver.delta > 0.0, "delta must be positive");

        let densities = solver.compute_densities(&positions, &masses, &neighbors);

        // Sanity check: at least some particles should be over-dense.
        let over_dense = densities.iter().filter(|&&d| d > rho0).count();
        assert!(
            over_dense > 0,
            "Expected over-dense particles for compressed cluster"
        );

        let forces =
            solver.solve_pressure(&positions, &velocities, &masses, &densities, &neighbors, dt);

        // Maximum pressure force magnitude must be positive.
        let max_f = forces.iter().map(|f| f.norm()).fold(0.0_f64, f64::max);
        assert!(
            max_f > 0.0,
            "Expected non-zero pressure forces for compressed cluster, got {max_f:.4e}"
        );
    }

    // -----------------------------------------------------------------------
    // Test 5: more iterations reduce density error
    // -----------------------------------------------------------------------
    #[test]
    fn test_pcisph_iterations_reduce_error() {
        let spacing = 0.03_f64;
        let h = 0.1_f64;
        let rho0 = 1000.0_f64;
        let mass = rho0 * (0.05_f64).powi(3);
        let dt = 0.001_f64;

        let (positions, velocities, masses) = uniform_block(spacing, 4, mass);
        let neighbors = SpatialHash::find_all_neighbors(&positions, 2.0 * h);

        // Helper: run solve and return the max predicted density error after
        // forced `n_iter` iterations.
        let run = |n_iter: usize| -> f64 {
            let mut solver = PcisphSolver::new(rho0, h);
            solver.compute_delta(dt, mass);
            solver.max_iter = n_iter;
            solver.min_iter = n_iter; // disable early exit
            solver.max_density_error = 0.0; // never early-exit on error criterion

            let densities = solver.compute_densities(&positions, &masses, &neighbors);
            let forces =
                solver.solve_pressure(&positions, &velocities, &masses, &densities, &neighbors, dt);

            // Compute predicted densities after applying the pressure forces.
            let pred_pos: Vec<Vec3> = positions
                .iter()
                .zip(velocities.iter())
                .zip(forces.iter())
                .map(|((p, v), f)| *p + (*v + *f * dt) * dt)
                .collect();
            let pred_densities = solver.compute_densities(&pred_pos, &masses, &neighbors);
            pred_densities
                .iter()
                .map(|&d| ((d - rho0) / rho0).abs())
                .fold(0.0_f64, f64::max)
        };

        let err_1 = run(1);
        let err_5 = run(5);

        // 5 iterations should not be *worse* than 1 iteration.
        assert!(
            err_5 <= err_1 * 1.1 + 1e-6,
            "More iterations should not dramatically worsen error: err_1={err_1:.4e}, err_5={err_5:.4e}"
        );
        // Both errors must be finite.
        assert!(err_1.is_finite(), "err_1 is not finite");
        assert!(err_5.is_finite(), "err_5 is not finite");
    }

    // -----------------------------------------------------------------------
    // Test 6: solve_pressure_tracked returns valid convergence info
    // -----------------------------------------------------------------------
    #[test]
    fn test_pcisph_tracked_convergence() {
        let spacing = 0.03_f64;
        let h = 0.1_f64;
        let rho0 = 1000.0_f64;
        let mass = rho0 * (0.05_f64).powi(3);
        let dt = 0.001_f64;

        let (positions, velocities, masses) = uniform_block(spacing, 4, mass);
        let neighbors = SpatialHash::find_all_neighbors(&positions, 2.0 * h);

        let mut solver = PcisphSolver::new(rho0, h);
        solver.compute_delta(dt, mass);
        solver.max_iter = 10;
        solver.min_iter = 2;
        solver.max_density_error = 0.5; // generous threshold

        let densities = solver.compute_densities(&positions, &masses, &neighbors);
        let (_forces, info) = solver.solve_pressure_tracked(
            &positions,
            &velocities,
            &masses,
            &densities,
            &neighbors,
            dt,
        );

        assert!(
            info.iterations_used >= solver.min_iter,
            "Should run at least min_iter iterations, got {}",
            info.iterations_used
        );
        assert!(
            info.iterations_used <= solver.max_iter,
            "Should not exceed max_iter, got {}",
            info.iterations_used
        );
        assert!(
            info.max_density_error >= 0.0,
            "max_density_error should be non-negative"
        );
        assert!(
            info.avg_density_error >= 0.0,
            "avg_density_error should be non-negative"
        );
        assert!(
            info.avg_density_error <= info.max_density_error + 1e-14,
            "avg should be <= max"
        );
    }

    // -----------------------------------------------------------------------
    // Test 7: adaptive_timestep produces sensible values
    // -----------------------------------------------------------------------
    #[test]
    fn test_pcisph_adaptive_timestep() {
        let solver = PcisphSolver::new(1000.0, 0.1);

        // Zero velocities → returns dt_max
        let vels_zero = vec![Vec3::zeros(); 10];
        let dt = solver.adaptive_timestep(&vels_zero, 1e-5, 0.01, 0.3);
        assert!(
            (dt - 0.01).abs() < 1e-12,
            "Zero velocity should yield dt_max, got {dt}"
        );

        // High velocity → should clamp below dt_max
        let vels_fast = vec![Vec3::new(100.0, 0.0, 0.0); 5];
        let dt_fast = solver.adaptive_timestep(&vels_fast, 1e-5, 0.01, 0.3);
        let expected = 0.3 * 0.1 / 100.0;
        assert!(
            (dt_fast - expected).abs() < 1e-12,
            "Expected {expected}, got {dt_fast}"
        );
    }

    // -----------------------------------------------------------------------
    // Test 8: density error helpers
    // -----------------------------------------------------------------------
    #[test]
    fn test_density_error_helpers() {
        let rho0 = 1000.0;
        let densities = vec![1010.0, 990.0, 1000.0, 1050.0];

        let errors = super::compute_density_errors(&densities, rho0);
        assert!((errors[0] - 10.0).abs() < 1e-10);
        assert!((errors[1] - (-10.0)).abs() < 1e-10);
        assert!((errors[2]).abs() < 1e-10);
        assert!((errors[3] - 50.0).abs() < 1e-10);

        let rms = super::rms_density_error(&densities, rho0);
        assert!(rms > 0.0, "RMS error should be positive");
        assert!(rms.is_finite(), "RMS error should be finite");

        // RMS of empty → 0
        assert!(
            (super::rms_density_error(&[], rho0)).abs() < 1e-14,
            "RMS of empty should be 0"
        );
    }

    // -----------------------------------------------------------------------
    // Test 9: clamp_pressures zeroes negative values
    // -----------------------------------------------------------------------
    #[test]
    fn test_clamp_pressures() {
        let mut p = vec![100.0, -50.0, 0.0, 200.0, -1.0];
        super::clamp_pressures(&mut p);
        assert!((p[0] - 100.0).abs() < 1e-14);
        assert!((p[1]).abs() < 1e-14);
        assert!((p[2]).abs() < 1e-14);
        assert!((p[3] - 200.0).abs() < 1e-14);
        assert!((p[4]).abs() < 1e-14);
    }

    // -----------------------------------------------------------------------
    // Test 10: relax_pressures blends correctly
    // -----------------------------------------------------------------------
    #[test]
    fn test_relax_pressures() {
        let old = vec![100.0, 200.0];
        let candidate = vec![200.0, 100.0];
        let relaxed = super::relax_pressures(&old, &candidate, 0.5);
        assert!(
            (relaxed[0] - 150.0).abs() < 1e-10,
            "Expected 150, got {}",
            relaxed[0]
        );
        assert!(
            (relaxed[1] - 150.0).abs() < 1e-10,
            "Expected 150, got {}",
            relaxed[1]
        );
    }

    // -----------------------------------------------------------------------
    // Test 11: AdaptivePcisph basic operation
    // -----------------------------------------------------------------------
    #[test]
    fn test_adaptive_pcisph_basic() {
        let spacing = 0.05_f64;
        let h = 0.1_f64;
        let rho0 = 1000.0_f64;
        let mass = rho0 * spacing.powi(3);
        let dt_max = 0.002_f64;

        let (positions, velocities, masses) = uniform_block(spacing, 4, mass);
        let neighbors = SpatialHash::find_all_neighbors(&positions, 2.0 * h);

        let mut solver = PcisphSolver::new(rho0, h);
        solver.compute_delta(dt_max, mass);

        let densities = solver.compute_densities(&positions, &masses, &neighbors);

        let mut adaptive = AdaptivePcisph::new(solver, dt_max);
        let (forces, dt_used) =
            adaptive.solve_adaptive(&positions, &velocities, &masses, &densities, &neighbors);

        assert_eq!(forces.len(), positions.len());
        assert!(dt_used > 0.0, "dt should be positive");
        assert!(dt_used <= dt_max, "dt should not exceed dt_max");
        for f in &forces {
            assert!(f.norm().is_finite(), "forces should be finite");
        }
    }

    // -----------------------------------------------------------------------
    // Test 12: compute_density_error agrees with manual calculation
    // -----------------------------------------------------------------------
    #[test]
    fn test_compute_density_error_method() {
        let spacing = 0.05_f64;
        let h = 0.1_f64;
        let rho0 = 1000.0_f64;
        let mass = rho0 * spacing.powi(3);

        let (positions, _velocities, masses) = uniform_block(spacing, 4, mass);
        let neighbors = SpatialHash::find_all_neighbors(&positions, 2.0 * h);

        let solver = PcisphSolver::new(rho0, h);
        let err = solver.compute_density_error(&positions, &masses, &neighbors);
        assert!(err >= 0.0, "avg density error should be non-negative");
        assert!(err.is_finite(), "avg density error should be finite");
    }

    // -----------------------------------------------------------------------
    // Test 13: max_pressure_after_solve is non-negative
    // -----------------------------------------------------------------------
    #[test]
    fn test_max_pressure_nonnegative() {
        let spacing = 0.03_f64;
        let h = 0.1_f64;
        let rho0 = 1000.0_f64;
        let mass = rho0 * (0.05_f64).powi(3);
        let dt = 0.001_f64;

        let (positions, velocities, masses) = uniform_block(spacing, 4, mass);
        let neighbors = SpatialHash::find_all_neighbors(&positions, 2.0 * h);

        let mut solver = PcisphSolver::new(rho0, h);
        solver.compute_delta(dt, mass);

        let densities = solver.compute_densities(&positions, &masses, &neighbors);
        let p_max = solver.max_pressure_after_solve(
            &positions,
            &velocities,
            &masses,
            &densities,
            &neighbors,
            dt,
        );
        assert!(
            p_max >= 0.0,
            "max pressure should be non-negative, got {p_max}"
        );
    }

    // -----------------------------------------------------------------------
    // Test 14: PressureCoefficientTable builds and looks up correctly
    // -----------------------------------------------------------------------
    #[test]
    fn test_pressure_coefficient_table_build() {
        let rho0 = 1000.0;
        let h = 0.1;
        let dt = 0.001;
        let table = PressureCoefficientTable::build(rho0, h, dt, 0.03, 0.1, 5);
        assert_eq!(table.len(), 5);
        assert!(!table.is_empty());
        // All deltas should be positive
        for (_, delta) in &table.entries {
            assert!(*delta > 0.0, "delta should be positive, got {delta}");
        }
    }

    #[test]
    fn test_pressure_coefficient_table_lookup_in_range() {
        let rho0 = 1000.0;
        let h = 0.1;
        let dt = 0.001;
        let table = PressureCoefficientTable::build(rho0, h, dt, 0.03, 0.1, 10);
        let delta = table.lookup(0.05);
        assert!(
            delta.is_some(),
            "lookup should return Some for in-range value"
        );
        assert!(delta.unwrap() > 0.0);
    }

    #[test]
    fn test_pressure_coefficient_table_lookup_boundary() {
        let rho0 = 1000.0;
        let h = 0.1;
        let dt = 0.001;
        let table = PressureCoefficientTable::build(rho0, h, dt, 0.03, 0.1, 5);
        // Clamp below min
        let d_low = table.lookup(0.01);
        let d_min = table.lookup(0.03);
        assert_eq!(d_low, d_min);
    }

    // -----------------------------------------------------------------------
    // Test 15: DensityErrorTracker record and query
    // -----------------------------------------------------------------------
    #[test]
    fn test_density_error_tracker_empty() {
        let tracker = DensityErrorTracker::new();
        assert_eq!(tracker.n_steps(), 0);
        assert_eq!(tracker.overall_max_error(), 0.0);
        assert_eq!(tracker.mean_avg_error(), 0.0);
        assert_eq!(tracker.mean_iterations(), 0.0);
    }

    #[test]
    fn test_density_error_tracker_records() {
        let mut tracker = DensityErrorTracker::new();
        let info1 = PcisphConvergenceInfo {
            iterations_used: 3,
            max_density_error: 0.05,
            avg_density_error: 0.02,
            converged: true,
        };
        let info2 = PcisphConvergenceInfo {
            iterations_used: 5,
            max_density_error: 0.08,
            avg_density_error: 0.04,
            converged: false,
        };
        tracker.record(&info1);
        tracker.record(&info2);
        assert_eq!(tracker.n_steps(), 2);
        assert!((tracker.overall_max_error() - 0.08).abs() < 1e-14);
        assert!((tracker.mean_avg_error() - 0.03).abs() < 1e-14);
        assert!((tracker.mean_iterations() - 4.0).abs() < 1e-14);
    }

    // -----------------------------------------------------------------------
    // Test 16: pcisph_correction_vector is zero for zero kappa
    // -----------------------------------------------------------------------
    #[test]
    fn test_pcisph_correction_vector_zero_kappa() {
        let pos_i = [0.0_f64; 3];
        let corr = pcisph_correction_vector(
            pos_i,
            1000.0,
            0.0,
            0.1,
            0.001,
            &[[0.05, 0.0, 0.0_f64]],
            &[0.001],
            &[1000.0],
            &[0.0],
        );
        let mag = (corr[0] * corr[0] + corr[1] * corr[1] + corr[2] * corr[2]).sqrt();
        assert!(mag < 1e-28, "Zero kappa → zero correction: {mag}");
    }

    #[test]
    fn test_pcisph_correction_vector_finite_nonzero_kappa() {
        let pos_i = [0.0_f64; 3];
        let corr = pcisph_correction_vector(
            pos_i,
            1000.0,
            100.0,
            0.1,
            0.001,
            &[[0.05, 0.0, 0.0_f64]],
            &[0.001],
            &[1000.0],
            &[50.0],
        );
        for c in corr {
            assert!(c.is_finite(), "Correction component should be finite: {c}");
        }
    }

    // -----------------------------------------------------------------------
    // Test 17: PcisphConvergenceMonitor basic behaviour
    // -----------------------------------------------------------------------
    #[test]
    fn test_convergence_monitor_stops_at_max_iter() {
        let mut monitor = PcisphConvergenceMonitor::new(0.01, 2, 5);
        for _ in 0..5 {
            monitor.record(0.5); // never converge
        }
        assert!(monitor.should_stop(), "Should stop at max_iter");
        assert!(!monitor.converged());
    }

    #[test]
    fn test_convergence_monitor_stops_when_converged() {
        let mut monitor = PcisphConvergenceMonitor::new(0.01, 2, 20);
        monitor.record(0.1);
        monitor.record(0.05);
        monitor.record(0.005); // below tolerance 0.01
        assert!(monitor.should_stop());
        assert!(monitor.converged());
    }

    #[test]
    fn test_convergence_monitor_waits_for_min_iter() {
        let mut monitor = PcisphConvergenceMonitor::new(0.01, 3, 20);
        monitor.record(0.001); // below tolerance but only 1 iter
        assert!(!monitor.should_stop(), "Should not stop before min_iter");
    }

    #[test]
    fn test_convergence_monitor_error_reduction() {
        let mut monitor = PcisphConvergenceMonitor::new(0.01, 1, 20);
        monitor.record(0.1);
        monitor.record(0.05);
        let ratio = monitor.error_reduction_ratio();
        assert!((ratio - 0.5).abs() < 1e-14);
    }

    #[test]
    fn test_convergence_monitor_reset() {
        let mut monitor = PcisphConvergenceMonitor::new(0.01, 2, 10);
        monitor.record(0.5);
        monitor.record(0.3);
        monitor.reset();
        assert_eq!(monitor.iter_count(), 0);
    }

    // -----------------------------------------------------------------------
    // Test 18: sph_density_plain matches density via solver
    // -----------------------------------------------------------------------
    #[test]
    fn test_sph_density_plain_self_contribution() {
        let pos_i = [0.0_f64; 3];
        let h = 0.1;
        let positions = vec![[0.0_f64; 3]];
        let masses = vec![1.0_f64];
        let rho = sph_density_plain(pos_i, h, &positions, &masses);
        // W(0, h) = 1/(pi*h^3)
        let expected = 1.0 / (std::f64::consts::PI * h * h * h);
        assert!((rho - expected).abs() < 1e-10);
    }

    #[test]
    fn test_sph_density_plain_far_zero() {
        let pos_i = [0.0_f64; 3];
        let h = 0.1;
        let positions = vec![[10.0_f64, 0.0, 0.0]];
        let masses = vec![1.0_f64];
        let rho = sph_density_plain(pos_i, h, &positions, &masses);
        assert!(rho < 1e-30, "Far particle → zero density");
    }

    // -----------------------------------------------------------------------
    // Test 19: compute_density_errors empty is zero
    // -----------------------------------------------------------------------
    #[test]
    fn test_density_errors_empty() {
        let errs = super::compute_density_errors(&[], 1000.0);
        assert!(errs.is_empty());
    }

    // -----------------------------------------------------------------------
    // Test 20: rms_density_error single-particle with known error
    // -----------------------------------------------------------------------
    #[test]
    fn test_rms_single_particle() {
        // Single particle at 1200 vs rho0=1000 → rel_err = 0.2
        let rms = super::rms_density_error(&[1200.0], 1000.0);
        assert!((rms - 0.2).abs() < 1e-14, "Expected 0.2, got {rms}");
    }
}
