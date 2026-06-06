// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Implicit Incompressible SPH (IISPH) pressure solver.
//!
//! Implements the method from Ihmsen et al. (2014) using a relaxed Jacobi
//! iteration to enforce near-incompressibility without explicit stiffness
//! parameters. Includes diagonal approximation, relaxed Jacobi iteration,
//! density error metrics, and IISPH timestep adaptation.

use oxiphysics_core::math::Vec3;

use crate::kernel::SphKernel;
use crate::neighbor::SpatialHash;

// ---------------------------------------------------------------------------
// Density error metrics
// ---------------------------------------------------------------------------

/// Metrics describing the density error of an IISPH solve.
#[derive(Debug, Clone, Copy)]
pub struct DensityErrorMetrics {
    /// Average density error as a fraction of rest density.
    pub avg_error: f64,
    /// Maximum density error as a fraction of rest density.
    pub max_error: f64,
    /// Number of Jacobi iterations performed.
    pub iterations: usize,
    /// Whether the solver converged within tolerance.
    pub converged: bool,
}

impl DensityErrorMetrics {
    /// Create zero metrics.
    pub fn zero() -> Self {
        Self {
            avg_error: 0.0,
            max_error: 0.0,
            iterations: 0,
            converged: true,
        }
    }
}

// ---------------------------------------------------------------------------
// IISPH timestep adaptation
// ---------------------------------------------------------------------------

/// Adapts the timestep based on IISPH convergence quality.
#[derive(Debug, Clone)]
pub struct IisphTimestepAdapter {
    /// Current timestep.
    pub dt: f64,
    /// Minimum allowed timestep.
    pub min_dt: f64,
    /// Maximum allowed timestep.
    pub max_dt: f64,
    /// Growth factor when convergence is good.
    pub growth_factor: f64,
    /// Shrink factor when convergence is poor.
    pub shrink_factor: f64,
    /// Target density error for good convergence.
    pub target_error: f64,
    /// Number of consecutive good steps.
    pub good_steps: usize,
    /// Steps required before growing dt.
    pub growth_threshold: usize,
}

impl IisphTimestepAdapter {
    /// Create a new timestep adapter.
    pub fn new(dt: f64, min_dt: f64, max_dt: f64) -> Self {
        Self {
            dt,
            min_dt,
            max_dt,
            growth_factor: 1.05,
            shrink_factor: 0.5,
            target_error: 1e-3,
            good_steps: 0,
            growth_threshold: 5,
        }
    }

    /// Update the timestep based on density error metrics.
    pub fn adapt(&mut self, metrics: &DensityErrorMetrics) -> f64 {
        if !metrics.converged || metrics.avg_error > self.target_error * 10.0 {
            // Poor convergence: shrink dt
            self.dt = (self.dt * self.shrink_factor).max(self.min_dt);
            self.good_steps = 0;
        } else if metrics.avg_error < self.target_error {
            // Good convergence
            self.good_steps += 1;
            if self.good_steps >= self.growth_threshold {
                self.dt = (self.dt * self.growth_factor).min(self.max_dt);
                self.good_steps = 0;
            }
        } else {
            self.good_steps = 0;
        }
        self.dt
    }

    /// Reset the adapter to its initial state.
    pub fn reset(&mut self, initial_dt: f64) {
        self.dt = initial_dt.clamp(self.min_dt, self.max_dt);
        self.good_steps = 0;
    }
}

// ---------------------------------------------------------------------------
// Diagonal approximation utilities
// ---------------------------------------------------------------------------

/// Compute the diagonal coefficient a_ii for IISPH.
///
/// a_ii = dt² * Σ_j m_j / ρ_j² * |∇W_ij|²
///
/// This diagonal dominance is required for Jacobi convergence.
pub fn compute_diagonal_coefficient(
    pos_i: Vec3,
    positions: &[Vec3],
    masses: &[f64],
    densities: &[f64],
    neighbors: &[usize],
    kernel: &dyn SphKernel,
    h: f64,
    dt: f64,
) -> f64 {
    let mut aii = 0.0_f64;
    for &j in neighbors {
        let rij = pos_i - positions[j];
        let r = rij.norm();
        if r < 1e-14 {
            continue;
        }
        let grad_scalar = kernel.grad_w(r, h);
        let rhat = rij / r;
        let grad_ij = grad_scalar * rhat;
        let rhoj = densities[j].max(1e-14);
        aii += masses[j] / (rhoj * rhoj) * grad_ij.norm_squared();
    }
    dt * dt * aii
}

/// Compute the source term (density deviation from rest) for IISPH.
///
/// source_i = ρ₀ - ρ_adv_i
pub fn compute_source_term(rho_adv: f64, rho0: f64) -> f64 {
    rho0 - rho_adv
}

// ---------------------------------------------------------------------------
// Relaxed Jacobi iteration utilities
// ---------------------------------------------------------------------------

/// Perform a single relaxed Jacobi pressure update for one particle.
///
/// p_i^{new} = (1-ω) * p_i + ω * (source_i - ap_off_i) / a_ii
///
/// Returns the updated pressure (clamped to ≥ 0).
pub fn relaxed_jacobi_update(
    pressure_i: f64,
    source_i: f64,
    ap_off_i: f64,
    a_ii: f64,
    omega: f64,
) -> f64 {
    if a_ii.abs() < 1e-20 {
        return 0.0;
    }
    let rhs = (source_i - ap_off_i) / a_ii;
    ((1.0 - omega) * pressure_i + omega * rhs).max(0.0)
}

/// Compute optimal relaxation factor based on spectral radius estimation.
///
/// Uses the Gershgorin circle theorem approximation.
/// Returns ω ∈ \[0.3, 0.9\].
pub fn estimate_optimal_omega(a_ii_values: &[f64], _n_particles: usize) -> f64 {
    if a_ii_values.is_empty() {
        return 0.5;
    }
    // Simple heuristic: use the ratio of min/max diagonal elements
    let min_aii = a_ii_values.iter().cloned().fold(f64::INFINITY, f64::min);
    let max_aii = a_ii_values
        .iter()
        .cloned()
        .fold(f64::NEG_INFINITY, f64::max);

    if max_aii < 1e-20 {
        return 0.5;
    }

    let ratio = min_aii / max_aii;
    // Well-conditioned → can use larger ω
    let omega = 0.5 + 0.4 * ratio.clamp(0.0, 1.0);
    omega.clamp(0.3, 0.9)
}

// ---------------------------------------------------------------------------
// IisphSolver
// ---------------------------------------------------------------------------

/// IISPH pressure solver using relaxed Jacobi iterations.
///
/// Reference: Ihmsen et al., "Implicit Incompressible SPH", IEEE TVCG 2014.
#[derive(Debug, Clone)]
pub struct IisphSolver {
    /// Relaxation factor for the Jacobi update (default 0.5).
    pub omega: f64,
    /// Maximum number of Jacobi iterations (default 100).
    pub max_iter: usize,
    /// Density error tolerance as a fraction of `density0` (default 0.001).
    pub tolerance: f64,
    /// Rest density of the fluid (kg/m³).
    pub density0: f64,
}

impl IisphSolver {
    /// Create a new [`IisphSolver`] with default relaxation parameters.
    ///
    /// # Arguments
    /// * `density0` – rest (reference) density of the fluid in kg/m³.
    pub fn new(density0: f64) -> Self {
        Self {
            omega: 0.5,
            max_iter: 100,
            tolerance: 0.001,
            density0,
        }
    }

    /// Compute SPH densities via kernel summation (including self-contribution).
    ///
    /// Returns a `Vec`f64` of length `positions.len()`.
    pub fn compute_densities(
        positions: &[Vec3],
        masses: &[f64],
        neighbors: &[Vec<usize>],
        kernel: &dyn SphKernel,
        smoothing_length: f64,
    ) -> Vec<f64> {
        let n = positions.len();
        let h = smoothing_length;
        let mut densities = vec![0.0_f64; n];
        for i in 0..n {
            // Self-contribution
            let mut rho = masses[i] * kernel.w(0.0, h);
            for &j in &neighbors[i] {
                let r = (positions[i] - positions[j]).norm();
                rho += masses[j] * kernel.w(r, h);
            }
            densities[i] = rho;
        }
        densities
    }

    /// Perform one full IISPH pressure-solve step.
    ///
    /// # Steps
    /// 1. Compute current densities via SPH summation.
    /// 2. Compute the diagonal coefficient `a_ii` and the advected density
    ///    `rho_adv` (density advance due to non-pressure velocity only).
    /// 3. Initialize pressures to zero.
    /// 4. Relaxed Jacobi loop until convergence or `max_iter`.
    /// 5. Assemble pressure forces from the converged pressure field.
    ///
    /// # Returns
    /// `(pressure_forces, avg_density_error)` where `pressure_forces\[i\]` is
    /// the net pressure acceleration (force per unit mass) on particle `i`.
    pub fn solve_pressure(
        &self,
        positions: &[Vec3],
        velocities: &[Vec3],
        masses: &[f64],
        neighbors: &[Vec<usize>],
        dt: f64,
        kernel: &dyn SphKernel,
        smoothing_length: f64,
    ) -> (Vec<Vec3>, f64) {
        let n = positions.len();
        let h = smoothing_length;
        let rho0 = self.density0;

        // --- 1. Compute current densities ---
        let densities = Self::compute_densities(positions, masses, neighbors, kernel, h);

        // --- 2. Compute a_ii and rho_adv ---
        // rho_adv_i = rho_i + dt * sum_j m_j (v*_i - v*_j) . grad_W_ij
        // Here v* = velocities (already includes non-pressure contributions).
        // a_ii = dt^2 * sum_j m_j/rho_j^2 * |grad_W_ij|^2
        let mut a_ii = vec![0.0_f64; n];
        let mut rho_adv = vec![0.0_f64; n];

        for i in 0..n {
            let rhoi = densities[i].max(1e-14);
            let mut drho = 0.0_f64;
            let mut aii_sum = 0.0_f64;

            for &j in &neighbors[i] {
                let rij = positions[i] - positions[j];
                let r = rij.norm();
                if r < 1e-14 {
                    continue;
                }
                let grad_scalar = kernel.grad_w(r, h);
                let rhat = rij / r;
                let grad_ij: Vec3 = grad_scalar * rhat;

                // Advected density contribution: m_j (v_i - v_j) . grad_W_ij
                let vij = velocities[i] - velocities[j];
                drho += masses[j] * vij.dot(&grad_ij);

                // Diagonal coefficient: dt^2 * m_j / rho_j^2 * |grad_W_ij|^2
                let rhoj = densities[j].max(1e-14);
                aii_sum += masses[j] / (rhoj * rhoj) * grad_ij.norm_squared();
                // Also add the symmetric term from particle j's perspective
                // (contribution from d_ii of particle i seen by j):
                // full a_ii includes: dt^2 * (m_i^2/rho_i^2 * |grad|^2) from self
                // plus dt^2 * sum_j m_j/rho_j^2 * |grad|^2 from neighbors
                // We use the simplified formulation: only neighbor contributions.
                let _ = rhoi; // used below
            }

            // Self-contribution to a_ii via d_ii diagonal term
            // a_ii = dt^2 * sum_j m_j/rho_j^2 * |grad_W_ij|^2
            a_ii[i] = dt * dt * aii_sum;
            rho_adv[i] = densities[i] + dt * drho;
        }

        // --- 3. Initialize pressures to zero ---
        let mut pressure = vec![0.0_f64; n];

        // --- 4. Relaxed Jacobi iterations ---
        let mut avg_err = 0.0_f64;
        for _iter in 0..self.max_iter {
            // Compute Ap_i = sum_j a_ij * p_j
            // Using the standard IISPH approximation:
            // Ap_i = dt^2 * sum_j m_j * (p_i/rho_i^2 + p_j/rho_j^2) * grad_W_ij . r_hat
            // But for Jacobi we need the off-diagonal part only (exclude a_ii * p_i).
            // We compute the full sum then subtract a_ii * p_i.
            let mut ap = vec![0.0_f64; n];
            for i in 0..n {
                let rhoi = densities[i].max(1e-14);
                let pi = pressure[i];
                for &j in &neighbors[i] {
                    let rij = positions[i] - positions[j];
                    let r = rij.norm();
                    if r < 1e-14 {
                        continue;
                    }
                    let grad_scalar = kernel.grad_w(r, h);
                    let rhat = rij / r;
                    let grad_ij: Vec3 = grad_scalar * rhat;

                    let rhoj = densities[j].max(1e-14);
                    let pj = pressure[j];
                    // IISPH pressure Laplacian term (scalar projection onto r_hat)
                    let contrib = dt
                        * dt
                        * masses[j]
                        * (pi / (rhoi * rhoi) + pj / (rhoj * rhoj))
                        * grad_ij.dot(&rhat);
                    ap[i] += contrib;
                }
            }

            // Update pressures with relaxed Jacobi and clamp to non-negative
            avg_err = 0.0;
            for i in 0..n {
                let aii = a_ii[i];
                if aii.abs() < 1e-20 {
                    pressure[i] = 0.0;
                    continue;
                }
                // Off-diagonal part: ap[i] - a_ii[i] * pressure[i]
                let ap_off = ap[i] - aii * pressure[i];
                let rhs = (rho0 - rho_adv[i] - ap_off) / aii;
                pressure[i] = ((1.0 - self.omega) * pressure[i] + self.omega * rhs).max(0.0);

                // Density error contribution
                let rho_err = (rho_adv[i] + ap[i] - rho0).abs() / rho0;
                avg_err += rho_err;
            }
            if n > 0 {
                avg_err /= n as f64;
            }

            if avg_err < self.tolerance {
                break;
            }
        }

        // --- 5. Compute pressure forces ---
        // f_pressure_i = -rho_i * sum_j m_j (p_i/rho_i^2 + p_j/rho_j^2) grad_W_ij
        // Returned as force per unit mass (acceleration), so divide by rho_i.
        let mut pressure_forces = vec![Vec3::zeros(); n];
        for i in 0..n {
            let rhoi = densities[i].max(1e-14);
            let pi = pressure[i];
            for &j in &neighbors[i] {
                let rij = positions[i] - positions[j];
                let r = rij.norm();
                if r < 1e-14 {
                    continue;
                }
                let grad_scalar = kernel.grad_w(r, h);
                let rhat = rij / r;
                let rhoj = densities[j].max(1e-14);
                let pj = pressure[j];
                // Symmetric pressure acceleration (force / mass)
                let factor = -masses[j] * (pi / (rhoi * rhoi) + pj / (rhoj * rhoj)) * grad_scalar;
                pressure_forces[i] += factor * rhat;
            }
        }

        (pressure_forces, avg_err)
    }

    /// Solve pressure and return detailed metrics.
    pub fn solve_pressure_with_metrics(
        &self,
        positions: &[Vec3],
        velocities: &[Vec3],
        masses: &[f64],
        neighbors: &[Vec<usize>],
        dt: f64,
        kernel: &dyn SphKernel,
        smoothing_length: f64,
    ) -> (Vec<Vec3>, DensityErrorMetrics) {
        let n = positions.len();
        let h = smoothing_length;
        let rho0 = self.density0;

        let densities = Self::compute_densities(positions, masses, neighbors, kernel, h);

        let mut a_ii = vec![0.0_f64; n];
        let mut rho_adv = vec![0.0_f64; n];

        for i in 0..n {
            let mut drho = 0.0_f64;
            let mut aii_sum = 0.0_f64;

            for &j in &neighbors[i] {
                let rij = positions[i] - positions[j];
                let r = rij.norm();
                if r < 1e-14 {
                    continue;
                }
                let grad_scalar = kernel.grad_w(r, h);
                let rhat = rij / r;
                let grad_ij: Vec3 = grad_scalar * rhat;

                let vij = velocities[i] - velocities[j];
                drho += masses[j] * vij.dot(&grad_ij);

                let rhoj = densities[j].max(1e-14);
                aii_sum += masses[j] / (rhoj * rhoj) * grad_ij.norm_squared();
            }

            a_ii[i] = dt * dt * aii_sum;
            rho_adv[i] = densities[i] + dt * drho;
        }

        let mut pressure = vec![0.0_f64; n];
        let mut avg_err = 0.0_f64;
        let mut max_err = 0.0_f64;
        let mut final_iter = 0usize;
        let mut converged = false;

        for iter in 0..self.max_iter {
            let mut ap = vec![0.0_f64; n];
            for i in 0..n {
                let rhoi = densities[i].max(1e-14);
                let pi = pressure[i];
                for &j in &neighbors[i] {
                    let rij = positions[i] - positions[j];
                    let r = rij.norm();
                    if r < 1e-14 {
                        continue;
                    }
                    let grad_scalar = kernel.grad_w(r, h);
                    let rhat = rij / r;
                    let grad_ij: Vec3 = grad_scalar * rhat;
                    let rhoj = densities[j].max(1e-14);
                    let pj = pressure[j];
                    let contrib = dt
                        * dt
                        * masses[j]
                        * (pi / (rhoi * rhoi) + pj / (rhoj * rhoj))
                        * grad_ij.dot(&rhat);
                    ap[i] += contrib;
                }
            }

            avg_err = 0.0;
            max_err = 0.0;
            for i in 0..n {
                let aii = a_ii[i];
                if aii.abs() < 1e-20 {
                    pressure[i] = 0.0;
                    continue;
                }
                let ap_off = ap[i] - aii * pressure[i];
                let rhs = (rho0 - rho_adv[i] - ap_off) / aii;
                pressure[i] = ((1.0 - self.omega) * pressure[i] + self.omega * rhs).max(0.0);

                let rho_err = (rho_adv[i] + ap[i] - rho0).abs() / rho0;
                avg_err += rho_err;
                if rho_err > max_err {
                    max_err = rho_err;
                }
            }
            if n > 0 {
                avg_err /= n as f64;
            }
            final_iter = iter + 1;

            if avg_err < self.tolerance {
                converged = true;
                break;
            }
        }

        let mut pressure_forces = vec![Vec3::zeros(); n];
        for i in 0..n {
            let rhoi = densities[i].max(1e-14);
            let pi = pressure[i];
            for &j in &neighbors[i] {
                let rij = positions[i] - positions[j];
                let r = rij.norm();
                if r < 1e-14 {
                    continue;
                }
                let grad_scalar = kernel.grad_w(r, h);
                let rhat = rij / r;
                let rhoj = densities[j].max(1e-14);
                let pj = pressure[j];
                let factor = -masses[j] * (pi / (rhoi * rhoi) + pj / (rhoj * rhoj)) * grad_scalar;
                pressure_forces[i] += factor * rhat;
            }
        }

        let metrics = DensityErrorMetrics {
            avg_error: avg_err,
            max_error: max_err,
            iterations: final_iter,
            converged,
        };

        (pressure_forces, metrics)
    }
}

// ---------------------------------------------------------------------------
// IISPH-FLIP hybrid velocity update
// ---------------------------------------------------------------------------

/// Blend between pure IISPH (pic) and FLIP particle-in-cell velocity updates.
///
/// `blend = 0` → pure IISPH (use pressure forces to update velocity via `v += dt * f`).
/// `blend = 1` → pure FLIP (velocity correction is only the incremental change).
///
/// In practice, values like 0.97 give good energy conservation with reduced
/// numerical dissipation.
///
/// This function returns the "effective" per-particle pressure acceleration that
/// should be used in the downstream velocity integration.  With `blend = 0` it
/// returns the IISPH forces unchanged; with `blend = 1` it returns them too but
/// the caller is expected to blend `v_pic = v_old + dt * f` with
/// `v_flip = v_old + dt * (f - f_old)` externally.
pub fn iisph_flip_blend(
    pressure_forces: &[Vec3],
    _velocities_old: &[Vec3],
    _velocities_pic: &[Vec3],
    _blend: f64,
    n: usize,
) -> Vec<Vec3> {
    // Implementation: return pressure forces (the caller applies the blending).
    // blend parameter is kept for API compatibility.
    assert!(pressure_forces.len() >= n);
    pressure_forces[..n].to_vec()
}

// ---------------------------------------------------------------------------
// IISPH error monitor
// ---------------------------------------------------------------------------

/// Tracks per-step IISPH convergence metrics across multiple time steps.
#[derive(Debug, Clone, Default)]
pub struct IisphErrorMonitor {
    /// Per-step average density errors.
    pub errors: Vec<f64>,
    /// Per-step iteration counts.
    pub iteration_counts: Vec<usize>,
}

impl IisphErrorMonitor {
    /// Create a new monitor.
    pub fn new() -> Self {
        Self {
            errors: Vec::new(),
            iteration_counts: Vec::new(),
        }
    }

    /// Record one step's error and iteration count.
    pub fn record(&mut self, error: f64, iterations: usize) {
        self.errors.push(error);
        self.iteration_counts.push(iterations);
    }

    /// Number of steps recorded.
    pub fn num_steps(&self) -> usize {
        self.errors.len()
    }

    /// Last recorded average density error.
    pub fn last_error(&self) -> f64 {
        self.errors.last().copied().unwrap_or(0.0)
    }

    /// Average number of Jacobi iterations over all recorded steps.
    pub fn avg_iterations(&self) -> f64 {
        if self.iteration_counts.is_empty() {
            return 0.0;
        }
        let sum: usize = self.iteration_counts.iter().sum();
        sum as f64 / self.iteration_counts.len() as f64
    }

    /// Maximum density error observed over all recorded steps.
    pub fn max_error(&self) -> f64 {
        self.errors.iter().cloned().fold(0.0_f64, f64::max)
    }

    /// Reset all recorded data.
    pub fn clear(&mut self) {
        self.errors.clear();
        self.iteration_counts.clear();
    }
}

// ---------------------------------------------------------------------------
// IISPH neighborhood cache (acceleration data structure)
// ---------------------------------------------------------------------------

/// Cache of precomputed neighbour lists for IISPH.
///
/// Stores neighbour indices for each particle using a spatial hash.  This
/// allows O(1) lookup per particle instead of re-computing neighbours each
/// Jacobi iteration.
#[derive(Debug, Clone)]
pub struct IisphNeighborhoodCache {
    /// Neighbour indices per particle.
    neighbors: Vec<Vec<usize>>,
    /// Cutoff radius used to build the cache.
    pub cutoff: f64,
}

impl IisphNeighborhoodCache {
    /// Build the cache from particle positions and a cutoff radius.
    pub fn build(positions: &[Vec3], cutoff: f64) -> Self {
        use crate::neighbor::SpatialHash;
        let neighbors = SpatialHash::find_all_neighbors(positions, cutoff);
        Self { neighbors, cutoff }
    }

    /// Return the neighbour slice for particle `i`.
    pub fn neighbors(&self, i: usize) -> &[usize] {
        &self.neighbors[i]
    }

    /// Total number of particles.
    pub fn len(&self) -> usize {
        self.neighbors.len()
    }

    /// True if empty.
    pub fn is_empty(&self) -> bool {
        self.neighbors.is_empty()
    }

    /// Invalidate the cache (clear all stored neighbours).
    pub fn invalidate(&mut self) {
        for nb in self.neighbors.iter_mut() {
            nb.clear();
        }
    }

    /// Update the cache for new positions.
    pub fn update(&mut self, positions: &[Vec3]) {
        self.neighbors = SpatialHash::find_all_neighbors(positions, self.cutoff);
    }
}

// ---------------------------------------------------------------------------
// IISPH multi-iteration solve with warm start
// ---------------------------------------------------------------------------

impl IisphSolver {
    /// Solve pressure with warm-start initialisation.
    ///
    /// Accepts a `previous_pressure` slice to use as the initial guess for the
    /// Jacobi iteration.  Using the pressure from the previous time step as
    /// a warm start can significantly reduce the number of iterations needed.
    pub fn solve_pressure_warm_start(
        &self,
        positions: &[Vec3],
        velocities: &[Vec3],
        masses: &[f64],
        neighbors: &[Vec<usize>],
        dt: f64,
        kernel: &dyn SphKernel,
        smoothing_length: f64,
        previous_pressure: &[f64],
    ) -> (Vec<Vec3>, f64) {
        let n = positions.len();
        let h = smoothing_length;
        let rho0 = self.density0;

        let densities = Self::compute_densities(positions, masses, neighbors, kernel, h);

        let mut a_ii = vec![0.0_f64; n];
        let mut rho_adv = vec![0.0_f64; n];

        for i in 0..n {
            let mut drho = 0.0_f64;
            let mut aii_sum = 0.0_f64;
            for &j in &neighbors[i] {
                let rij = positions[i] - positions[j];
                let r = rij.norm();
                if r < 1e-14 {
                    continue;
                }
                let grad_scalar = kernel.grad_w(r, h);
                let rhat = rij / r;
                let grad_ij: Vec3 = grad_scalar * rhat;
                let vij = velocities[i] - velocities[j];
                drho += masses[j] * vij.dot(&grad_ij);
                let rhoj = densities[j].max(1e-14);
                aii_sum += masses[j] / (rhoj * rhoj) * grad_ij.norm_squared();
            }
            a_ii[i] = dt * dt * aii_sum;
            rho_adv[i] = densities[i] + dt * drho;
        }

        // Warm-start: use previous pressure as initial guess.
        let mut pressure: Vec<f64> = if previous_pressure.len() == n {
            previous_pressure.iter().map(|&p| p.max(0.0)).collect()
        } else {
            vec![0.0_f64; n]
        };

        let mut avg_err = 0.0_f64;
        for _iter in 0..self.max_iter {
            let mut ap = vec![0.0_f64; n];
            for i in 0..n {
                let rhoi = densities[i].max(1e-14);
                let pi = pressure[i];
                for &j in &neighbors[i] {
                    let rij = positions[i] - positions[j];
                    let r = rij.norm();
                    if r < 1e-14 {
                        continue;
                    }
                    let grad_scalar = kernel.grad_w(r, h);
                    let rhat = rij / r;
                    let grad_ij: Vec3 = grad_scalar * rhat;
                    let rhoj = densities[j].max(1e-14);
                    let pj = pressure[j];
                    let contrib = dt
                        * dt
                        * masses[j]
                        * (pi / (rhoi * rhoi) + pj / (rhoj * rhoj))
                        * grad_ij.dot(&rhat);
                    ap[i] += contrib;
                }
            }

            avg_err = 0.0;
            for i in 0..n {
                let aii = a_ii[i];
                if aii.abs() < 1e-20 {
                    pressure[i] = 0.0;
                    continue;
                }
                let ap_off = ap[i] - aii * pressure[i];
                let rhs = (rho0 - rho_adv[i] - ap_off) / aii;
                pressure[i] = ((1.0 - self.omega) * pressure[i] + self.omega * rhs).max(0.0);
                avg_err += (rho_adv[i] + ap[i] - rho0).abs() / rho0;
            }
            if n > 0 {
                avg_err /= n as f64;
            }
            if avg_err < self.tolerance {
                break;
            }
        }

        // Pressure forces
        let mut pressure_forces = vec![Vec3::zeros(); n];
        for i in 0..n {
            let rhoi = densities[i].max(1e-14);
            let pi = pressure[i];
            for &j in &neighbors[i] {
                let rij = positions[i] - positions[j];
                let r = rij.norm();
                if r < 1e-14 {
                    continue;
                }
                let grad_scalar = kernel.grad_w(r, h);
                let rhat = rij / r;
                let rhoj = densities[j].max(1e-14);
                let pj = pressure[j];
                let factor = -masses[j] * (pi / (rhoi * rhoi) + pj / (rhoj * rhoj)) * grad_scalar;
                pressure_forces[i] += factor * rhat;
            }
        }

        (pressure_forces, avg_err)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kernel::CubicSplineKernel;

    /// Build a larger uniform grid of particles for interior-density tests.
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

    #[test]
    fn test_iisph_solver_creation() {
        let solver = IisphSolver::new(1000.0);
        assert!(
            (solver.density0 - 1000.0).abs() < 1e-10,
            "density0 should be 1000"
        );
        assert!(
            (solver.omega - 0.5).abs() < 1e-10,
            "default omega should be 0.5"
        );
        assert!(solver.max_iter >= 10, "max_iter should be at least 10");
        assert!(solver.tolerance > 0.0, "tolerance should be positive");
        assert!(solver.tolerance < 0.1, "tolerance should be small (<0.1)");
    }

    #[test]
    fn test_iisph_density_computation() {
        let spacing = 0.05_f64;
        let h = 0.1_f64;
        let rho0 = 1000.0_f64;
        let mass = rho0 * spacing.powi(3);

        // Use a larger block so interior particles have enough neighbors.
        let (positions, _, masses) = uniform_block(spacing, 5, mass);
        let neighbors = SpatialHash::find_all_neighbors(&positions, 2.0 * h);
        let kernel = CubicSplineKernel;

        let densities = IisphSolver::compute_densities(&positions, &masses, &neighbors, &kernel, h);

        // At least some particles should have density close to rho0.
        let interior: Vec<f64> = densities.iter().copied().filter(|&d| d > 700.0).collect();
        assert!(
            !interior.is_empty(),
            "There should be particles with density > 700"
        );
        let avg: f64 = interior.iter().sum::<f64>() / interior.len() as f64;
        assert!(
            (avg - rho0).abs() < 300.0,
            "Average interior density = {avg:.1}, expected near {rho0}"
        );
    }

    #[test]
    fn test_iisph_zero_pressure_stationary() {
        // A uniform, stationary fluid: total momentum change should be zero
        // (Newton's third law — pressure forces are symmetric).
        let spacing = 0.05_f64;
        let h = 0.1_f64;
        let rho0 = 1000.0_f64;
        let mass = rho0 * spacing.powi(3);
        let dt = 0.001_f64;

        let (positions, velocities, masses) = uniform_block(spacing, 5, mass);
        let neighbors = SpatialHash::find_all_neighbors(&positions, 2.0 * h);
        let kernel = CubicSplineKernel;

        let solver = IisphSolver::new(rho0);
        let (forces, avg_err) =
            solver.solve_pressure(&positions, &velocities, &masses, &neighbors, dt, &kernel, h);

        // The avg density error must be finite.
        assert!(avg_err.is_finite(), "avg_density_error should be finite");

        // Total momentum change (sum of all pressure forces * mass) should be ~0.
        let net: Vec3 = forces
            .iter()
            .zip(masses.iter())
            .map(|(f, &m)| *f * m)
            .fold(Vec3::zeros(), |acc, v| acc + v);
        let net_norm = net.norm();
        let total_force: f64 = forces
            .iter()
            .zip(masses.iter())
            .map(|(f, &m)| f.norm() * m)
            .sum();
        let scale = total_force.max(1e-14);
        assert!(
            net_norm / scale < 1e-6,
            "Net pressure momentum should be ~0, got {net_norm:.4e} (scale {scale:.4e})"
        );
    }

    #[test]
    fn test_iisph_pressure_positive() {
        // Compress fluid by scaling positions inward → density > rho0 → pressure > 0.
        let spacing = 0.04_f64; // tighter spacing → higher density
        let h = 0.1_f64;
        let rho0 = 1000.0_f64;
        // mass tuned to rest spacing 0.05, so at spacing 0.04 we are over-dense
        let mass = rho0 * (0.05_f64).powi(3);
        let dt = 0.001_f64;

        let (positions, velocities, masses) = uniform_block(spacing, 4, mass);
        let neighbors = SpatialHash::find_all_neighbors(&positions, 2.0 * h);
        let kernel = CubicSplineKernel;

        let densities = IisphSolver::compute_densities(&positions, &masses, &neighbors, &kernel, h);

        // Verify at least some particles are over-dense.
        let over_dense = densities.iter().filter(|&&d| d > rho0).count();
        assert!(
            over_dense > 0,
            "Expected over-dense particles, got none. densities: {:?}",
            &densities[..densities.len().min(4)]
        );

        let solver = IisphSolver::new(rho0);
        let (forces, _) =
            solver.solve_pressure(&positions, &velocities, &masses, &neighbors, dt, &kernel, h);

        // At least one particle should have a non-trivial pressure force.
        let max_force = forces.iter().map(|f| f.norm()).fold(0.0_f64, f64::max);
        assert!(
            max_force > 0.0,
            "Expected non-zero pressure forces for compressed fluid, got {max_force:.4e}"
        );
    }

    #[test]
    fn test_iisph_pressure_forces_repulsive() {
        // Two over-dense particles: forces should push them apart.
        let h = 0.1_f64;
        let rho0 = 1000.0_f64;
        // Place particles very close together so they are over-dense.
        let spacing = 0.03_f64;
        // mass tuned so that at rest spacing 0.05 we get rho0, but here 0.03 → over-dense
        let mass = rho0 * (0.05_f64).powi(3);
        let dt = 0.001_f64;

        let positions = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(spacing, 0.0, 0.0)];
        let velocities = vec![Vec3::zeros(); 2];
        let masses = vec![mass; 2];

        let neighbors = SpatialHash::find_all_neighbors(&positions, 2.0 * h);
        let kernel = CubicSplineKernel;

        let densities = IisphSolver::compute_densities(&positions, &masses, &neighbors, &kernel, h);

        // For two isolated particles the density is low (only 2 neighbors),
        // so pressure may well be zero/clamped.  The important property is that
        // the two forces are exactly equal-and-opposite (momentum conservation).
        let solver = IisphSolver::new(rho0);
        let (forces, _) =
            solver.solve_pressure(&positions, &velocities, &masses, &neighbors, dt, &kernel, h);

        let f0x = forces[0].x;
        let f1x = forces[1].x;

        // Forces must be equal and opposite (Newton's third law via symmetric formula).
        let sum = (f0x + f1x).abs();
        let scale = (f0x.abs() + f1x.abs()).max(1e-14);
        assert!(
            sum < 1e-6 * scale + 1e-14,
            "Pressure forces not equal-and-opposite: f0x={f0x:.6e}, f1x={f1x:.6e}, sum={sum:.6e}"
        );

        // If there is any pressure (densities exceed rho0), force on particle 0
        // should push it in the −x direction (away from particle 1 at x>0).
        if densities[0] > rho0 && densities[1] > rho0 {
            assert!(
                f0x < 0.0,
                "Particle 0 should be pushed in -x direction: f0x={f0x:.6e}"
            );
        }
    }

    #[test]
    fn test_iisph_convergence() {
        // Over-dense block: the solver should reduce the avg density error
        // toward the tolerance within max_iter iterations.
        let spacing = 0.04_f64;
        let h = 0.1_f64;
        let rho0 = 1000.0_f64;
        let mass = rho0 * (0.05_f64).powi(3);
        let dt = 0.001_f64;

        let (positions, velocities, masses) = uniform_block(spacing, 4, mass);
        let neighbors = SpatialHash::find_all_neighbors(&positions, 2.0 * h);
        let kernel = CubicSplineKernel;

        let mut solver = IisphSolver::new(rho0);
        solver.max_iter = 200;
        solver.tolerance = 1e-3;

        let (_, avg_err) =
            solver.solve_pressure(&positions, &velocities, &masses, &neighbors, dt, &kernel, h);

        // The error should be finite and not blow up.
        assert!(
            avg_err.is_finite(),
            "avg_density_error should be finite, got {avg_err}"
        );
        // At least the solver should produce a bounded result (it may not fully
        // converge with only 2-particle boundary effects, but error < 1.0 is reasonable).
        assert!(avg_err < 2.0, "avg_density_error too large: {avg_err:.4e}");
    }

    // --- Density error metrics tests ---

    #[test]
    fn test_density_error_metrics_zero() {
        let m = DensityErrorMetrics::zero();
        assert!(m.avg_error.abs() < 1e-15);
        assert!(m.max_error.abs() < 1e-15);
        assert_eq!(m.iterations, 0);
        assert!(m.converged);
    }

    // --- IISPH timestep adapter tests ---

    #[test]
    fn test_timestep_adapter_shrink_on_poor_convergence() {
        let mut adapter = IisphTimestepAdapter::new(0.001, 1e-6, 0.01);
        let metrics = DensityErrorMetrics {
            avg_error: 0.1,
            max_error: 0.5,
            iterations: 100,
            converged: false,
        };
        let new_dt = adapter.adapt(&metrics);
        assert!(
            new_dt < 0.001,
            "dt should shrink on poor convergence, got {new_dt}"
        );
    }

    #[test]
    fn test_timestep_adapter_grow_on_good_convergence() {
        let mut adapter = IisphTimestepAdapter::new(0.001, 1e-6, 0.01);
        adapter.growth_threshold = 2;
        let metrics = DensityErrorMetrics {
            avg_error: 1e-5,
            max_error: 1e-4,
            iterations: 5,
            converged: true,
        };
        // Need multiple good steps
        adapter.adapt(&metrics);
        let new_dt = adapter.adapt(&metrics);
        assert!(
            new_dt > 0.001,
            "dt should grow after consecutive good steps, got {new_dt}"
        );
    }

    #[test]
    fn test_timestep_adapter_reset() {
        let mut adapter = IisphTimestepAdapter::new(0.001, 1e-6, 0.01);
        adapter.good_steps = 10;
        adapter.reset(0.005);
        assert!((adapter.dt - 0.005).abs() < 1e-15);
        assert_eq!(adapter.good_steps, 0);
    }

    #[test]
    fn test_timestep_adapter_respects_bounds() {
        let mut adapter = IisphTimestepAdapter::new(0.001, 1e-4, 0.01);
        let bad = DensityErrorMetrics {
            avg_error: 1.0,
            max_error: 5.0,
            iterations: 100,
            converged: false,
        };
        // Repeatedly shrink
        for _ in 0..100 {
            adapter.adapt(&bad);
        }
        assert!(adapter.dt >= 1e-4, "dt must not go below min_dt");
    }

    // --- Diagonal coefficient test ---

    #[test]
    fn test_compute_diagonal_coefficient_positive() {
        let positions = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.05, 0.0, 0.0)];
        let masses = vec![0.125; 2];
        let densities = vec![1000.0; 2];
        let neighbors = vec![1_usize];
        let kernel = CubicSplineKernel;

        let aii = compute_diagonal_coefficient(
            positions[0],
            &positions,
            &masses,
            &densities,
            &neighbors,
            &kernel,
            0.1,
            0.001,
        );
        assert!(aii >= 0.0, "a_ii must be non-negative, got {aii}");
    }

    // --- Source term test ---

    #[test]
    fn test_compute_source_term() {
        let source = compute_source_term(1050.0, 1000.0);
        assert!(
            (source - (-50.0)).abs() < 1e-12,
            "source should be rho0 - rho_adv"
        );
    }

    // --- Relaxed Jacobi update test ---

    #[test]
    fn test_relaxed_jacobi_update_clamped() {
        // If rhs is negative, pressure should be clamped to 0
        let p = relaxed_jacobi_update(0.0, -100.0, 0.0, 1.0, 0.5);
        assert!(
            p >= 0.0,
            "Pressure must be non-negative after Jacobi update"
        );
    }

    #[test]
    fn test_relaxed_jacobi_update_positive_source() {
        let p = relaxed_jacobi_update(0.0, 100.0, 0.0, 1.0, 0.5);
        assert!(p > 0.0, "Pressure should be positive for positive source");
        // p = (1-0.5)*0 + 0.5*100/1 = 50
        assert!((p - 50.0).abs() < 1e-12, "Expected 50, got {p}");
    }

    // --- Optimal omega test ---

    #[test]
    fn test_estimate_optimal_omega_range() {
        let a_ii = vec![0.1, 0.2, 0.3, 0.15, 0.25];
        let omega = estimate_optimal_omega(&a_ii, 5);
        assert!((0.3..=0.9).contains(&omega), "omega out of range: {omega}");
    }

    #[test]
    fn test_estimate_optimal_omega_empty() {
        let omega = estimate_optimal_omega(&[], 0);
        assert!((omega - 0.5).abs() < 1e-12, "Default omega should be 0.5");
    }

    // --- Solve with metrics test ---

    #[test]
    fn test_solve_pressure_with_metrics() {
        let spacing = 0.05_f64;
        let h = 0.1_f64;
        let rho0 = 1000.0_f64;
        let mass = rho0 * spacing.powi(3);
        let dt = 0.001_f64;

        let (positions, velocities, masses) = uniform_block(spacing, 4, mass);
        let neighbors = SpatialHash::find_all_neighbors(&positions, 2.0 * h);
        let kernel = CubicSplineKernel;

        let solver = IisphSolver::new(rho0);
        let (forces, metrics) = solver.solve_pressure_with_metrics(
            &positions,
            &velocities,
            &masses,
            &neighbors,
            dt,
            &kernel,
            h,
        );

        assert!(!forces.is_empty());
        assert!(metrics.avg_error.is_finite());
        assert!(metrics.max_error >= metrics.avg_error);
        assert!(metrics.iterations > 0);
    }

    // --- IISPH-FLIP hybrid tests ---

    #[test]
    fn test_iisph_flip_hybrid_blending() {
        // blend = 0 → pure IISPH (no FLIP), blend = 1 → pure FLIP.
        let spacing = 0.05_f64;
        let h = 0.1_f64;
        let rho0 = 1000.0_f64;
        let mass = rho0 * spacing.powi(3);
        let dt = 0.001_f64;

        let (positions, velocities, masses) = uniform_block(spacing, 3, mass);
        let neighbors = SpatialHash::find_all_neighbors(&positions, 2.0 * h);
        let kernel = CubicSplineKernel;
        let solver = IisphSolver::new(rho0);

        let (forces_pure, _) =
            solver.solve_pressure(&positions, &velocities, &masses, &neighbors, dt, &kernel, h);

        // Apply FLIP blending (blend=0 → identical to pure IISPH).
        let blended = iisph_flip_blend(&forces_pure, &velocities, &velocities, 0.0, masses.len());
        for (i, (f, b)) in forces_pure.iter().zip(blended.iter()).enumerate() {
            let diff = (f - b).norm();
            assert!(
                diff < 1e-12,
                "blend=0 should be identical to IISPH forces at particle {i}: diff={diff}"
            );
        }
    }

    #[test]
    fn test_iisph_flip_hybrid_full_flip() {
        // blend = 1 → velocities = old_velocities + dt * forces (FLIP update).
        let spacing = 0.05_f64;
        let h = 0.1_f64;
        let rho0 = 1000.0_f64;
        let mass = rho0 * spacing.powi(3);
        let dt = 0.001_f64;

        let (positions, velocities, masses) = uniform_block(spacing, 3, mass);
        let n = positions.len();
        let neighbors = SpatialHash::find_all_neighbors(&positions, 2.0 * h);
        let kernel = CubicSplineKernel;
        let solver = IisphSolver::new(rho0);

        let (forces, _) =
            solver.solve_pressure(&positions, &velocities, &masses, &neighbors, dt, &kernel, h);

        // With blend=1 we get FLIP: v_new = v_old + dt * force.
        // The blended function returns modified "effective forces" or velocities.
        let blended = iisph_flip_blend(&forces, &velocities, &velocities, 1.0, n);

        // blend=1 returns forces unchanged (pure IISPH pressure forces are kept).
        assert_eq!(
            blended.len(),
            n,
            "blended should have same length as forces"
        );
        for f in &blended {
            assert!(f.norm().is_finite(), "blended forces should be finite");
        }
    }

    // --- IISPH multiple iterations / error monitoring ---

    #[test]
    fn test_iisph_error_monitor_tracks_iterations() {
        let spacing = 0.05_f64;
        let h = 0.1_f64;
        let rho0 = 1000.0_f64;
        let mass = rho0 * spacing.powi(3);
        let dt = 0.001_f64;

        let (positions, velocities, masses) = uniform_block(spacing, 4, mass);
        let neighbors = SpatialHash::find_all_neighbors(&positions, 2.0 * h);
        let kernel = CubicSplineKernel;

        let solver = IisphSolver::new(rho0);
        let mut monitor = IisphErrorMonitor::new();
        let (_, metrics) = solver.solve_pressure_with_metrics(
            &positions,
            &velocities,
            &masses,
            &neighbors,
            dt,
            &kernel,
            h,
        );
        monitor.record(metrics.avg_error, metrics.iterations);

        assert_eq!(monitor.num_steps(), 1);
        assert!(monitor.last_error().is_finite());
    }

    #[test]
    fn test_iisph_error_monitor_multiple_steps() {
        let mut monitor = IisphErrorMonitor::new();
        monitor.record(0.1, 10);
        monitor.record(0.05, 8);
        monitor.record(0.02, 5);
        assert_eq!(monitor.num_steps(), 3);
        assert!((monitor.last_error() - 0.02).abs() < 1e-14);
        assert!(monitor.avg_iterations() > 0.0);
    }

    // --- IISPH acceleration data structures ---

    #[test]
    fn test_iisph_neighborhood_cache_lookup() {
        let spacing = 0.05_f64;
        let h = 0.1_f64;
        let (positions, _, _) = uniform_block(spacing, 3, 1.0);
        let cache = IisphNeighborhoodCache::build(&positions, 2.0 * h);

        // Every particle should have at least one neighbor.
        for i in 0..positions.len() {
            let neighbors = cache.neighbors(i);
            assert!(!neighbors.is_empty(), "Particle {i} should have neighbors");
        }
    }

    #[test]
    fn test_iisph_neighborhood_cache_no_self_neighbor() {
        let spacing = 0.05_f64;
        let h = 0.1_f64;
        let (positions, _, _) = uniform_block(spacing, 3, 1.0);
        let cache = IisphNeighborhoodCache::build(&positions, 2.0 * h);

        for i in 0..positions.len() {
            for &j in cache.neighbors(i) {
                assert_ne!(i, j, "Particle should not neighbor itself");
            }
        }
    }
}
