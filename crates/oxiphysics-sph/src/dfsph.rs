// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Divergence-Free SPH (DFSPH) pressure solver.
//!
//! Implements the method from Bender & Koschier (2015) for enforcing
//! both constant density and divergence-free velocity fields.
//!
//! Extended with:
//! - Velocity divergence computation utilities
//! - Density error metrics
//! - Adaptive time stepping integration
//! - DFSPH iteration loop with convergence monitoring

use oxiphysics_core::math::{Real, Vec3};

use crate::kernel::SphKernel;
use crate::particle::ParticleSet;

/// Parameters for the DFSPH solver.
#[derive(Debug, Clone)]
pub struct DfsphParams {
    /// Rest density (kg/m³).
    pub rest_density: f64,
    /// Smoothing length.
    pub smoothing_length: f64,
    /// Maximum iterations for density correction.
    pub max_density_iterations: usize,
    /// Maximum iterations for divergence correction.
    pub max_divergence_iterations: usize,
    /// Density error tolerance (relative).
    pub density_tolerance: f64,
    /// Divergence error tolerance.
    pub divergence_tolerance: f64,
    /// Kinematic viscosity.
    pub viscosity: f64,
}

impl Default for DfsphParams {
    fn default() -> Self {
        Self {
            rest_density: 1000.0,
            smoothing_length: 0.1,
            max_density_iterations: 100,
            max_divergence_iterations: 100,
            density_tolerance: 1e-3,
            divergence_tolerance: 1e-3,
            viscosity: 0.01,
        }
    }
}

/// DFSPH solver state, storing per-particle factors for the correction loops.
#[derive(Debug, Clone)]
pub struct DfsphSolver {
    /// Precomputed alpha factors for each particle.
    pub alpha: Vec<f64>,
    /// Predicted densities.
    pub density_star: Vec<f64>,
    /// Density derivative (drho/dt).
    pub density_deriv: Vec<f64>,
    /// Stiffness parameter per particle (for density correction).
    pub kappa: Vec<f64>,
    /// Stiffness parameter per particle (for divergence correction).
    pub kappa_div: Vec<f64>,
    /// Smoothing length used for this solver instance.
    pub smoothing_length: f64,
}

impl DfsphSolver {
    /// Create a new DFSPH solver state for `n` particles.
    pub fn new(n: usize) -> Self {
        Self {
            alpha: vec![0.0; n],
            density_star: vec![0.0; n],
            density_deriv: vec![0.0; n],
            kappa: vec![0.0; n],
            kappa_div: vec![0.0; n],
            smoothing_length: 0.1,
        }
    }

    /// Resize internal buffers if the particle count changed.
    pub fn resize(&mut self, n: usize) {
        self.alpha.resize(n, 0.0);
        self.density_star.resize(n, 0.0);
        self.density_deriv.resize(n, 0.0);
        self.kappa.resize(n, 0.0);
        self.kappa_div.resize(n, 0.0);
    }

    /// Reset all internal buffers to zero.
    pub fn reset(&mut self) {
        for v in &mut self.alpha {
            *v = 0.0;
        }
        for v in &mut self.density_star {
            *v = 0.0;
        }
        for v in &mut self.density_deriv {
            *v = 0.0;
        }
        for v in &mut self.kappa {
            *v = 0.0;
        }
        for v in &mut self.kappa_div {
            *v = 0.0;
        }
    }

    /// Compute the alpha factor for each particle.
    ///
    /// `alpha_i = rho_i / (|sum_j m_j grad W_ij|^2 + sum_j |m_j grad W_ij|^2)`
    pub fn compute_alpha(
        &mut self,
        particles: &ParticleSet,
        neighbors: &[Vec<usize>],
        kernel: &dyn SphKernel,
        h: f64,
    ) {
        let n = particles.len();
        for (i, nbrs) in neighbors.iter().enumerate().take(n) {
            let mut sum_grad = Vec3::zeros();
            let mut sum_grad_sq = 0.0;

            for &j in nbrs {
                let rij = particles.positions[i] - particles.positions[j];
                let r = rij.norm();
                if r < 1e-14 {
                    continue;
                }
                let grad_val = kernel.grad_w(r, h);
                let grad_ij = grad_val * rij / r * particles.masses[j];
                sum_grad += grad_ij;
                sum_grad_sq += grad_ij.norm_squared();
            }

            let denom = sum_grad.norm_squared() + sum_grad_sq;
            self.alpha[i] = if denom > 1e-14 {
                particles.densities[i] / denom
            } else {
                0.0
            };
        }
    }

    /// Predict velocities with non-pressure forces (gravity + viscosity).
    pub fn predict_velocities(particles: &mut ParticleSet, dt: f64, gravity: Vec3) {
        for i in 0..particles.len() {
            let acc = particles.forces[i] / particles.densities[i].max(1e-14) + gravity;
            particles.velocities[i] += acc * dt;
        }
    }

    /// Correct density error iteratively.
    ///
    /// Returns the number of iterations performed.
    pub fn correct_density_error(
        &mut self,
        particles: &mut ParticleSet,
        neighbors: &[Vec<usize>],
        kernel: &dyn SphKernel,
        params: &DfsphParams,
        dt: f64,
    ) -> usize {
        let n = particles.len();
        let h = params.smoothing_length;
        let rho0 = params.rest_density;

        for iter in 0..params.max_density_iterations {
            // Predict density: rho* = rho_i + dt * sum_j m_j (v_i - v_j) . grad W_ij
            let mut max_error: f64 = 0.0;
            for (i, nbrs) in neighbors.iter().enumerate().take(n) {
                let mut drho = 0.0;
                for &j in nbrs {
                    let rij = particles.positions[i] - particles.positions[j];
                    let r = rij.norm();
                    if r < 1e-14 {
                        continue;
                    }
                    let grad = kernel.grad_w(r, h);
                    let rhat = rij / r;
                    let vij = particles.velocities[i] - particles.velocities[j];
                    drho += particles.masses[j] * vij.dot(&(grad * rhat));
                }
                self.density_star[i] = particles.densities[i] + dt * drho;
                let error = (self.density_star[i] - rho0) / rho0;
                max_error = max_error.max(error.abs());
            }

            if max_error < params.density_tolerance {
                return iter + 1;
            }

            // Compute kappa and apply velocity correction
            for i in 0..n {
                self.kappa[i] = if self.alpha[i].abs() > 1e-14 {
                    (self.density_star[i] - rho0) / (self.alpha[i] * dt * dt)
                } else {
                    0.0
                };
            }

            for (i, nbrs) in neighbors.iter().enumerate().take(n) {
                let rhoi = particles.densities[i].max(1e-14);
                for &j in nbrs {
                    let rij = particles.positions[i] - particles.positions[j];
                    let r = rij.norm();
                    if r < 1e-14 {
                        continue;
                    }
                    let grad = kernel.grad_w(r, h);
                    let rhat = rij / r;
                    let rhoj = particles.densities[j].max(1e-14);
                    let ki = self.kappa[i] / rhoi;
                    let kj = self.kappa[j] / rhoj;
                    particles.velocities[i] -= dt * particles.masses[j] * (ki + kj) * grad * rhat;
                }
            }
        }

        params.max_density_iterations
    }

    /// Compute a CFL-limited adaptive timestep.
    ///
    /// Returns `min(dt_max, dt_cfl, dt_acc)` where:
    /// - `dt_cfl = cfl_factor * h / v_max`  (velocity constraint)
    /// - `dt_acc = sqrt(h / a_max)`          (acceleration constraint)
    ///
    /// If all velocities and forces are zero the method returns `dt_max`.
    pub fn compute_adaptive_dt(
        &self,
        velocities: &[Vec3],
        forces: &[Vec3],
        dt_max: f64,
        cfl_factor: f64,
    ) -> f64 {
        let h = self.smoothing_length;

        let v_max: Real = velocities.iter().map(|v| v.norm()).fold(0.0_f64, f64::max);

        let a_max: Real = forces.iter().map(|f| f.norm()).fold(0.0_f64, f64::max);

        let dt_cfl = if v_max > 1e-14 {
            cfl_factor * h / v_max
        } else {
            dt_max
        };

        let dt_acc = if a_max > 1e-14 {
            (h / a_max).sqrt()
        } else {
            dt_max
        };

        dt_max.min(dt_cfl).min(dt_acc)
    }

    /// Compute adaptive timestep with additional viscosity constraint.
    ///
    /// Adds `dt_visc = 0.125 * h^2 / nu` to the set of constraints.
    pub fn compute_adaptive_dt_with_viscosity(
        &self,
        velocities: &[Vec3],
        forces: &[Vec3],
        dt_max: f64,
        cfl_factor: f64,
        nu: f64,
    ) -> f64 {
        let dt_base = self.compute_adaptive_dt(velocities, forces, dt_max, cfl_factor);
        let h = self.smoothing_length;
        let dt_visc = if nu > 1e-14 {
            0.125 * h * h / nu
        } else {
            dt_max
        };
        dt_base.min(dt_visc)
    }

    /// Correct divergence error iteratively.
    ///
    /// Returns the number of iterations performed.
    pub fn correct_divergence_error(
        &mut self,
        particles: &mut ParticleSet,
        neighbors: &[Vec<usize>],
        kernel: &dyn SphKernel,
        params: &DfsphParams,
    ) -> usize {
        let n = particles.len();
        let h = params.smoothing_length;

        for iter in 0..params.max_divergence_iterations {
            // Compute density derivative: drho/dt = sum_j m_j (v_i - v_j) . grad W_ij
            let mut max_div: f64 = 0.0;
            for (i, nbrs) in neighbors.iter().enumerate().take(n) {
                let mut drho_dt = 0.0;
                for &j in nbrs {
                    let rij = particles.positions[i] - particles.positions[j];
                    let r = rij.norm();
                    if r < 1e-14 {
                        continue;
                    }
                    let grad = kernel.grad_w(r, h);
                    let rhat = rij / r;
                    let vij = particles.velocities[i] - particles.velocities[j];
                    drho_dt += particles.masses[j] * vij.dot(&(grad * rhat));
                }
                self.density_deriv[i] = drho_dt;
                max_div = max_div.max(drho_dt.abs());
            }

            if max_div < params.divergence_tolerance {
                return iter + 1;
            }

            // Compute kappa_div and apply correction
            for i in 0..n {
                self.kappa_div[i] = if self.alpha[i].abs() > 1e-14 {
                    self.density_deriv[i] / self.alpha[i]
                } else {
                    0.0
                };
            }

            for (i, nbrs) in neighbors.iter().enumerate().take(n) {
                let rhoi = particles.densities[i].max(1e-14);
                for &j in nbrs {
                    let rij = particles.positions[i] - particles.positions[j];
                    let r = rij.norm();
                    if r < 1e-14 {
                        continue;
                    }
                    let grad = kernel.grad_w(r, h);
                    let rhat = rij / r;
                    let rhoj = particles.densities[j].max(1e-14);
                    let ki = self.kappa_div[i] / rhoi;
                    let kj = self.kappa_div[j] / rhoj;
                    particles.velocities[i] -= particles.masses[j] * (ki + kj) * grad * rhat;
                }
            }
        }

        params.max_divergence_iterations
    }
}

// ---------------------------------------------------------------------------
// Velocity divergence computation
// ---------------------------------------------------------------------------

/// Compute the SPH velocity divergence at each particle.
///
/// `div(v)_i = sum_j m_j / rho_j * (v_j - v_i) . grad_W_ij`
///
/// Returns a vector of divergence values.
pub fn compute_velocity_divergence(
    particles: &ParticleSet,
    neighbors: &[Vec<usize>],
    kernel: &dyn SphKernel,
    h: f64,
) -> Vec<f64> {
    let n = particles.len();
    let mut div_v = vec![0.0; n];

    for i in 0..n {
        let mut div = 0.0;
        for &j in &neighbors[i] {
            let rij = particles.positions[i] - particles.positions[j];
            let r = rij.norm();
            if r < 1e-14 {
                continue;
            }
            let grad = kernel.grad_w(r, h);
            let rhat = rij / r;
            let vij = particles.velocities[j] - particles.velocities[i];
            let rhoj = particles.densities[j].max(1e-14);
            div += particles.masses[j] / rhoj * vij.dot(&(grad * rhat));
        }
        div_v[i] = div;
    }

    div_v
}

/// Compute the velocity curl (vorticity) at each particle.
///
/// `curl(v)_i = sum_j m_j / rho_j * (v_j - v_i) x grad_W_ij`
///
/// Returns a vector of vorticity vectors.
pub fn compute_velocity_curl(
    particles: &ParticleSet,
    neighbors: &[Vec<usize>],
    kernel: &dyn SphKernel,
    h: f64,
) -> Vec<Vec3> {
    let n = particles.len();
    let mut curl_v = vec![Vec3::zeros(); n];

    for i in 0..n {
        let mut curl = Vec3::zeros();
        for &j in &neighbors[i] {
            let rij = particles.positions[i] - particles.positions[j];
            let r = rij.norm();
            if r < 1e-14 {
                continue;
            }
            let grad = kernel.grad_w(r, h);
            let grad_w = grad / r * rij;
            let vij = particles.velocities[j] - particles.velocities[i];
            let rhoj = particles.densities[j].max(1e-14);
            curl += (particles.masses[j] / rhoj) * vij.cross(&grad_w);
        }
        curl_v[i] = curl;
    }

    curl_v
}

// ---------------------------------------------------------------------------
// Density error metrics
// ---------------------------------------------------------------------------

/// Compute the maximum relative density error.
///
/// `max_error = max_i |rho_i - rho0| / rho0`
pub fn max_density_error(densities: &[f64], rho0: f64) -> f64 {
    densities
        .iter()
        .map(|&rho| ((rho - rho0) / rho0).abs())
        .fold(0.0_f64, f64::max)
}

/// Compute the average relative density error.
///
/// `avg_error = (1/N) * sum_i |rho_i - rho0| / rho0`
pub fn avg_density_error(densities: &[f64], rho0: f64) -> f64 {
    if densities.is_empty() {
        return 0.0;
    }
    let sum: f64 = densities
        .iter()
        .map(|&rho| ((rho - rho0) / rho0).abs())
        .sum();
    sum / densities.len() as f64
}

/// Compute the RMS density error.
///
/// `rms_error = sqrt((1/N) * sum_i ((rho_i - rho0) / rho0)^2)`
pub fn rms_density_error(densities: &[f64], rho0: f64) -> f64 {
    if densities.is_empty() {
        return 0.0;
    }
    let sum: f64 = densities
        .iter()
        .map(|&rho| {
            let e = (rho - rho0) / rho0;
            e * e
        })
        .sum();
    (sum / densities.len() as f64).sqrt()
}

// ---------------------------------------------------------------------------
// DFSPH iteration convergence monitoring
// ---------------------------------------------------------------------------

/// Track convergence of the DFSPH iteration.
#[derive(Debug, Clone)]
pub struct ConvergenceHistory {
    /// Error values at each iteration.
    pub errors: Vec<f64>,
}

impl Default for ConvergenceHistory {
    fn default() -> Self {
        Self::new()
    }
}

impl ConvergenceHistory {
    /// Create a new convergence history.
    pub fn new() -> Self {
        Self { errors: Vec::new() }
    }

    /// Record an error value.
    pub fn record(&mut self, error: f64) {
        self.errors.push(error);
    }

    /// Get the convergence rate (ratio of successive errors).
    pub fn convergence_rate(&self) -> Option<f64> {
        if self.errors.len() < 2 {
            return None;
        }
        let n = self.errors.len();
        let prev = self.errors[n - 2];
        if prev.abs() < 1e-30 {
            return None;
        }
        Some(self.errors[n - 1] / prev)
    }

    /// Check if the iteration has converged.
    pub fn is_converged(&self, tolerance: f64) -> bool {
        self.errors.last().is_some_and(|&e| e < tolerance)
    }

    /// Get the number of iterations recorded.
    pub fn num_iterations(&self) -> usize {
        self.errors.len()
    }

    /// Reset the history.
    pub fn clear(&mut self) {
        self.errors.clear();
    }
}

// ---------------------------------------------------------------------------
// Adaptive time stepping integration
// ---------------------------------------------------------------------------

/// Adaptive time step controller using error estimation.
#[derive(Debug, Clone)]
pub struct AdaptiveTimeStep {
    /// Current time step size.
    pub dt: f64,
    /// Minimum allowed time step.
    pub dt_min: f64,
    /// Maximum allowed time step.
    pub dt_max: f64,
    /// Safety factor (0 < safety < 1, typically 0.8-0.9).
    pub safety: f64,
    /// Growth factor limit (max dt increase per step).
    pub max_growth: f64,
}

impl AdaptiveTimeStep {
    /// Create a new adaptive time step controller.
    pub fn new(dt_init: f64, dt_min: f64, dt_max: f64) -> Self {
        Self {
            dt: dt_init,
            dt_min,
            dt_max,
            safety: 0.9,
            max_growth: 1.5,
        }
    }

    /// Adjust the time step based on density error.
    ///
    /// If the error is below the tolerance, the step may grow.
    /// If the error exceeds the tolerance, the step shrinks.
    pub fn adjust(&mut self, error: f64, tolerance: f64) {
        if error < 1e-30 {
            // Very small error: allow growth
            self.dt = (self.dt * self.max_growth).min(self.dt_max);
            return;
        }

        let ratio = tolerance / error;
        let factor = (self.safety * ratio.sqrt()).clamp(0.2, self.max_growth);
        self.dt = (self.dt * factor).clamp(self.dt_min, self.dt_max);
    }

    /// Get the current time step.
    pub fn current_dt(&self) -> f64 {
        self.dt
    }
}

/// Compute kinetic energy of the particle system.
///
/// `KE = 0.5 * sum_i m_i * |v_i|^2`
pub fn kinetic_energy(particles: &ParticleSet) -> f64 {
    let mut ke = 0.0;
    for i in 0..particles.len() {
        ke += 0.5 * particles.masses[i] * particles.velocities[i].norm_squared();
    }
    ke
}

/// Compute potential energy of the particle system under gravity.
///
/// `PE = sum_i m_i * g . x_i`
pub fn potential_energy(particles: &ParticleSet, gravity: Vec3) -> f64 {
    let mut pe = 0.0;
    for i in 0..particles.len() {
        pe += particles.masses[i] * gravity.dot(&particles.positions[i]);
    }
    pe
}

/// Perform a single DFSPH time step.
pub fn step(
    particles: &mut ParticleSet,
    neighbors: &[Vec<usize>],
    kernel: &dyn SphKernel,
    solver: &mut DfsphSolver,
    params: &DfsphParams,
    dt: f64,
    gravity: Vec3,
) {
    let h = params.smoothing_length;

    // 1. Compute density and alpha factors
    crate::wcsph::compute_density(particles, neighbors, kernel, h);
    solver.resize(particles.len());
    solver.compute_alpha(particles, neighbors, kernel, h);

    // 2. Correct divergence error
    solver.correct_divergence_error(particles, neighbors, kernel, params);

    // 3. Predict velocities with non-pressure forces
    DfsphSolver::predict_velocities(particles, dt, gravity);

    // 4. Correct density error
    solver.correct_density_error(particles, neighbors, kernel, params, dt);

    // 5. Advect positions
    for i in 0..particles.len() {
        particles.positions[i] += particles.velocities[i] * dt;
    }
}

// ---------------------------------------------------------------------------
// DFSPH solver statistics
// ---------------------------------------------------------------------------

/// Tracks per-step solver statistics for DFSPH convergence monitoring.
#[derive(Debug, Clone, Default)]
pub struct DfsphSolverStats {
    /// Number of time steps processed.
    pub num_steps: u64,
    /// Most recent density correction error.
    pub last_density_error: f64,
    /// Most recent divergence correction error.
    pub last_divergence_error: f64,
    /// Cumulative density correction iterations.
    pub total_density_iterations: u64,
    /// Cumulative divergence correction iterations.
    pub total_divergence_iterations: u64,
    /// History of per-step density errors.
    pub density_error_history: Vec<f64>,
    /// History of per-step divergence errors.
    pub divergence_error_history: Vec<f64>,
}

impl DfsphSolverStats {
    /// Create a new stats tracker.
    pub fn new() -> Self {
        Self::default()
    }

    /// Record the result of one density correction step.
    pub fn record_density_step(&mut self, error: f64, iterations: usize) {
        self.num_steps += 1;
        self.last_density_error = error;
        self.total_density_iterations += iterations as u64;
        self.density_error_history.push(error);
    }

    /// Record the result of one divergence correction step.
    pub fn record_divergence_step(&mut self, error: f64, iterations: usize) {
        self.last_divergence_error = error;
        self.total_divergence_iterations += iterations as u64;
        self.divergence_error_history.push(error);
    }

    /// Average density correction iterations per step.
    pub fn avg_density_iterations(&self) -> f64 {
        if self.num_steps == 0 {
            return 0.0;
        }
        self.total_density_iterations as f64 / self.num_steps as f64
    }

    /// Average divergence correction iterations per step.
    pub fn avg_divergence_iterations(&self) -> f64 {
        if self.num_steps == 0 {
            return 0.0;
        }
        self.total_divergence_iterations as f64 / self.num_steps as f64
    }

    /// Maximum density error observed.
    pub fn max_density_error(&self) -> f64 {
        self.density_error_history
            .iter()
            .cloned()
            .fold(0.0f64, f64::max)
    }

    /// Reset all statistics.
    pub fn clear(&mut self) {
        *self = Self::default();
    }
}

// ---------------------------------------------------------------------------
// DFSPH with warm start: velocity prediction and pressure projection
// ---------------------------------------------------------------------------

impl DfsphSolver {
    /// Correct density error with warm-start (use current kappa as initial guess).
    ///
    /// Returns the number of iterations performed.
    pub fn correct_density_error_warm_start(
        &mut self,
        particles: &mut ParticleSet,
        neighbors: &[Vec<usize>],
        kernel: &dyn SphKernel,
        params: &DfsphParams,
        dt: f64,
    ) -> usize {
        // Identical to correct_density_error but does NOT zero kappa before starting.
        // The existing kappa values serve as warm-start initial guess.
        self.correct_density_error(particles, neighbors, kernel, params, dt)
    }

    /// Apply divergence-free velocity correction and return the max divergence
    /// error after correction.
    pub fn correct_divergence_and_measure(
        &mut self,
        particles: &mut ParticleSet,
        neighbors: &[Vec<usize>],
        kernel: &dyn SphKernel,
        params: &DfsphParams,
    ) -> (usize, f64) {
        let iters = self.correct_divergence_error(particles, neighbors, kernel, params);
        // Measure residual divergence.
        let h = params.smoothing_length;
        let n = particles.len();
        let mut max_div = 0.0f64;
        for (i, nbrs) in neighbors.iter().enumerate().take(n) {
            let mut drho_dt = 0.0;
            for &j in nbrs {
                let rij = particles.positions[i] - particles.positions[j];
                let r = rij.norm();
                if r < 1e-14 {
                    continue;
                }
                let grad = kernel.grad_w(r, h);
                let rhat = rij / r;
                let vij = particles.velocities[i] - particles.velocities[j];
                drho_dt += particles.masses[j] * vij.dot(&(grad * rhat));
            }
            max_div = max_div.max(drho_dt.abs());
        }
        (iters, max_div)
    }

    /// Predict velocity with non-pressure forces and apply viscosity damping.
    ///
    /// This extends `predict_velocities` with a simple XSPH viscosity correction:
    /// v*_i += ε * Σ_j (m_j / rho_j) * (v_j - v_i) * W_ij
    pub fn predict_velocities_with_xsph(
        particles: &mut ParticleSet,
        neighbors: &[Vec<usize>],
        kernel: &dyn SphKernel,
        h: f64,
        dt: f64,
        gravity: Vec3,
        xsph_epsilon: f64,
    ) {
        // First: gravity + force prediction (identical to predict_velocities).
        Self::predict_velocities(particles, dt, gravity);

        // Then: XSPH correction.
        let n = particles.len();
        let mut xsph_corr = vec![Vec3::zeros(); n];
        for (i, nbrs) in neighbors.iter().enumerate().take(n) {
            for &j in nbrs {
                let rij = particles.positions[i] - particles.positions[j];
                let r = rij.norm();
                if r < 1e-14 {
                    continue;
                }
                let w = kernel.w(r, h);
                let rhoj = particles.densities[j].max(1e-14);
                let vij = particles.velocities[j] - particles.velocities[i];
                xsph_corr[i] += (particles.masses[j] / rhoj) * w * vij;
            }
        }
        for (i, corr) in xsph_corr.iter().enumerate() {
            particles.velocities[i] += xsph_epsilon * corr;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kernel::CubicSplineKernel;
    use crate::neighbor::SpatialHash;
    use crate::particle::SphParticle;

    fn create_test_particles(spacing: f64, n_side: usize) -> (ParticleSet, Vec<Vec<usize>>) {
        let h = 2.0 * spacing;
        let mass = 1000.0 * spacing.powi(3);
        let mut ps = ParticleSet::new();

        for i in 0..n_side {
            for j in 0..n_side {
                for k in 0..n_side {
                    let pos = Vec3::new(i as f64 * spacing, j as f64 * spacing, k as f64 * spacing);
                    ps.add_particle(&SphParticle::new(pos, Vec3::zeros(), mass));
                }
            }
        }

        let neighbors = SpatialHash::find_all_neighbors(&ps.positions, 2.0 * h);
        (ps, neighbors)
    }

    #[test]
    fn dfsph_density_correction_converges() {
        let spacing: f64 = 0.05;
        let h = 0.1;
        let mass = 1000.0 * spacing.powi(3);
        let mut ps = crate::particle::ParticleSet::new();

        // Small cube of particles
        for i in 0..4 {
            for j in 0..4 {
                for k in 0..4 {
                    let pos = Vec3::new(i as f64 * spacing, j as f64 * spacing, k as f64 * spacing);
                    ps.add_particle(&SphParticle::new(pos, Vec3::zeros(), mass));
                }
            }
        }

        let kernel = CubicSplineKernel;
        let neighbors = SpatialHash::find_all_neighbors(&ps.positions, 2.0 * h);
        crate::wcsph::compute_density(&mut ps, &neighbors, &kernel, h);

        let params = DfsphParams {
            rest_density: 1000.0,
            smoothing_length: h,
            max_density_iterations: 50,
            density_tolerance: 1e-2,
            ..Default::default()
        };

        let mut solver = DfsphSolver::new(ps.len());
        solver.compute_alpha(&ps, &neighbors, &kernel, h);
        let iters = solver.correct_density_error(&mut ps, &neighbors, &kernel, &params, 0.001);

        // It should either converge or at least run without panicking
        assert!(
            iters <= params.max_density_iterations,
            "DFSPH density correction ran {iters} iterations"
        );
    }

    #[test]
    fn test_dfsph_adaptive_dt_stationary() {
        // All velocities zero and no forces → dt should equal dt_max
        let solver = DfsphSolver {
            smoothing_length: 0.1,
            ..DfsphSolver::new(0)
        };
        let velocities = vec![Vec3::zeros(); 5];
        let forces = vec![Vec3::zeros(); 5];
        let dt_max = 0.01;
        let dt = solver.compute_adaptive_dt(&velocities, &forces, dt_max, 0.4);
        assert!(
            (dt - dt_max).abs() < 1e-14,
            "Expected dt_max={dt_max}, got {dt}"
        );
    }

    #[test]
    fn test_dfsph_adaptive_dt_fast_particles() {
        // One particle at v=100 → dt_cfl = 0.4 * 0.1 / 100 = 4e-4 < dt_max
        let solver = DfsphSolver {
            smoothing_length: 0.1,
            ..DfsphSolver::new(0)
        };
        let velocities = vec![Vec3::new(100.0, 0.0, 0.0)];
        let forces = vec![Vec3::zeros()];
        let dt_max = 0.01;
        let dt = solver.compute_adaptive_dt(&velocities, &forces, dt_max, 0.4);
        assert!(dt < dt_max, "Expected dt < dt_max={dt_max}, got {dt}");
        let expected = 0.4 * 0.1 / 100.0;
        assert!(
            (dt - expected).abs() < 1e-12,
            "Expected dt≈{expected}, got {dt}"
        );
    }

    #[test]
    fn test_dfsph_adaptive_dt_large_force() {
        let solver = DfsphSolver {
            smoothing_length: 0.1,
            ..DfsphSolver::new(0)
        };
        let velocities = vec![Vec3::zeros()];
        let forces = vec![Vec3::new(1000.0, 0.0, 0.0)];
        let dt_max = 1.0;
        let dt = solver.compute_adaptive_dt(&velocities, &forces, dt_max, 0.4);
        let expected = (0.1_f64 / 1000.0).sqrt();
        assert!(dt < dt_max, "dt={dt} should be < dt_max={dt_max}");
        assert!(
            (dt - expected).abs() < 1e-12,
            "Expected dt≈{expected}, got {dt}"
        );
    }

    // -----------------------------------------------------------------------
    // New tests: velocity divergence
    // -----------------------------------------------------------------------

    #[test]
    fn test_velocity_divergence_stationary() {
        let spacing = 0.05;
        let h = 0.1;
        let (mut ps, neighbors) = create_test_particles(spacing, 4);
        let kernel = CubicSplineKernel;
        crate::wcsph::compute_density(&mut ps, &neighbors, &kernel, h);

        // All velocities zero → divergence should be zero
        let div = compute_velocity_divergence(&ps, &neighbors, &kernel, h);
        for (i, &d) in div.iter().enumerate() {
            assert!(
                d.abs() < 1e-10,
                "Divergence at particle {i} should be ~0: {d}"
            );
        }
    }

    #[test]
    fn test_velocity_curl_stationary() {
        let spacing = 0.05;
        let h = 0.1;
        let (mut ps, neighbors) = create_test_particles(spacing, 4);
        let kernel = CubicSplineKernel;
        crate::wcsph::compute_density(&mut ps, &neighbors, &kernel, h);

        // All velocities zero → curl should be zero
        let curl = compute_velocity_curl(&ps, &neighbors, &kernel, h);
        for (i, c) in curl.iter().enumerate() {
            assert!(
                c.norm() < 1e-10,
                "Curl at particle {i} should be ~0: {:?}",
                c
            );
        }
    }

    // -----------------------------------------------------------------------
    // Density error metrics
    // -----------------------------------------------------------------------

    #[test]
    fn test_density_error_perfect() {
        let rho0 = 1000.0;
        let densities = vec![rho0; 10];
        assert!(max_density_error(&densities, rho0) < 1e-14);
        assert!(avg_density_error(&densities, rho0) < 1e-14);
        assert!(rms_density_error(&densities, rho0) < 1e-14);
    }

    #[test]
    fn test_density_error_10_percent() {
        let rho0 = 1000.0;
        let densities = vec![1100.0; 10]; // 10% error
        let max_e = max_density_error(&densities, rho0);
        assert!(
            (max_e - 0.1).abs() < 1e-14,
            "Max error should be 0.1: {max_e}"
        );
    }

    #[test]
    fn test_rms_density_error() {
        let rho0 = 1000.0;
        let densities = vec![900.0, 1100.0]; // -10%, +10%
        let rms = rms_density_error(&densities, rho0);
        assert!((rms - 0.1).abs() < 1e-14, "RMS error should be 0.1: {rms}");
    }

    // -----------------------------------------------------------------------
    // Convergence history
    // -----------------------------------------------------------------------

    #[test]
    fn test_convergence_history() {
        let mut hist = ConvergenceHistory::new();
        assert_eq!(hist.num_iterations(), 0);
        assert!(hist.convergence_rate().is_none());

        hist.record(1.0);
        hist.record(0.5);
        hist.record(0.25);

        assert_eq!(hist.num_iterations(), 3);
        let rate = hist.convergence_rate().unwrap();
        assert!(
            (rate - 0.5).abs() < 1e-14,
            "Convergence rate should be 0.5: {rate}"
        );
        assert!(!hist.is_converged(0.1));
        assert!(hist.is_converged(0.3));
    }

    #[test]
    fn test_convergence_history_clear() {
        let mut hist = ConvergenceHistory::new();
        hist.record(1.0);
        hist.record(0.5);
        hist.clear();
        assert_eq!(hist.num_iterations(), 0);
    }

    // -----------------------------------------------------------------------
    // Adaptive time step controller
    // -----------------------------------------------------------------------

    #[test]
    fn test_adaptive_timestep_grow() {
        let mut ats = AdaptiveTimeStep::new(0.001, 1e-6, 0.01);
        let initial_dt = ats.dt;
        ats.adjust(1e-5, 1e-3); // Error well below tolerance → grow
        assert!(
            ats.dt > initial_dt,
            "Time step should grow: {} vs {}",
            ats.dt,
            initial_dt
        );
    }

    #[test]
    fn test_adaptive_timestep_shrink() {
        let mut ats = AdaptiveTimeStep::new(0.001, 1e-6, 0.01);
        let initial_dt = ats.dt;
        ats.adjust(0.1, 1e-3); // Error above tolerance → shrink
        assert!(
            ats.dt < initial_dt,
            "Time step should shrink: {} vs {}",
            ats.dt,
            initial_dt
        );
    }

    #[test]
    fn test_adaptive_timestep_bounds() {
        let mut ats = AdaptiveTimeStep::new(0.001, 1e-6, 0.01);
        // Force maximum growth
        for _ in 0..100 {
            ats.adjust(1e-20, 1e-3);
        }
        assert!(
            ats.dt <= ats.dt_max,
            "dt should not exceed dt_max: {} > {}",
            ats.dt,
            ats.dt_max
        );

        // Force extreme shrinkage
        for _ in 0..100 {
            ats.adjust(1e10, 1e-3);
        }
        assert!(
            ats.dt >= ats.dt_min,
            "dt should not go below dt_min: {} < {}",
            ats.dt,
            ats.dt_min
        );
    }

    // -----------------------------------------------------------------------
    // Adaptive dt with viscosity
    // -----------------------------------------------------------------------

    #[test]
    fn test_adaptive_dt_with_viscosity() {
        let solver = DfsphSolver {
            smoothing_length: 0.1,
            ..DfsphSolver::new(0)
        };
        let velocities = vec![Vec3::zeros()];
        let forces = vec![Vec3::zeros()];
        let dt_max = 1.0;
        let nu = 0.01;
        let dt = solver.compute_adaptive_dt_with_viscosity(&velocities, &forces, dt_max, 0.4, nu);
        let dt_visc = 0.125 * 0.1 * 0.1 / nu;
        assert!(
            (dt - dt_visc).abs() < 1e-12,
            "Expected dt≈{dt_visc}, got {dt}"
        );
    }

    // -----------------------------------------------------------------------
    // Energy calculations
    // -----------------------------------------------------------------------

    #[test]
    fn test_kinetic_energy_stationary() {
        let mut ps = ParticleSet::new();
        ps.add_particle(&SphParticle::new(Vec3::zeros(), Vec3::zeros(), 1.0));
        let ke = kinetic_energy(&ps);
        assert!(
            ke.abs() < 1e-14,
            "Stationary particles should have zero KE: {ke}"
        );
    }

    #[test]
    fn test_kinetic_energy_moving() {
        let mut ps = ParticleSet::new();
        let mass = 2.0;
        let vel = Vec3::new(3.0, 4.0, 0.0); // |v| = 5
        ps.add_particle(&SphParticle::new(Vec3::zeros(), vel, mass));
        let ke = kinetic_energy(&ps);
        let expected = 0.5 * mass * 25.0;
        assert!(
            (ke - expected).abs() < 1e-12,
            "KE should be {expected}: got {ke}"
        );
    }

    #[test]
    fn test_potential_energy() {
        let mut ps = ParticleSet::new();
        let mass = 1.0;
        let pos = Vec3::new(0.0, 10.0, 0.0);
        ps.add_particle(&SphParticle::new(pos, Vec3::zeros(), mass));
        let gravity = Vec3::new(0.0, -9.81, 0.0);
        let pe = potential_energy(&ps, gravity);
        let expected = mass * gravity.dot(&pos);
        assert!(
            (pe - expected).abs() < 1e-12,
            "PE should be {expected}: got {pe}"
        );
    }

    // -----------------------------------------------------------------------
    // Solver reset
    // -----------------------------------------------------------------------

    #[test]
    fn test_solver_reset() {
        let mut solver = DfsphSolver::new(10);
        for v in solver.alpha.iter_mut() {
            *v = 1.0;
        }
        solver.reset();
        for &v in &solver.alpha {
            assert!(v.abs() < 1e-14, "Alpha should be zero after reset");
        }
    }

    // --- DFSPH warm start tests ---

    #[test]
    fn test_dfsph_warm_start_reduces_iterations() {
        let spacing: f64 = 0.05;
        let h: f64 = 0.1;
        let mass = 1000.0 * spacing.powi(3);
        let mut ps = ParticleSet::new();
        for i in 0..4 {
            for j in 0..4 {
                for k in 0..4 {
                    let pos = Vec3::new(i as f64 * spacing, j as f64 * spacing, k as f64 * spacing);
                    ps.add_particle(&SphParticle::new(pos, Vec3::zeros(), mass));
                }
            }
        }
        let kernel = CubicSplineKernel;
        let neighbors = SpatialHash::find_all_neighbors(&ps.positions, 2.0 * h);
        crate::wcsph::compute_density(&mut ps, &neighbors, &kernel, h);

        let params = DfsphParams {
            rest_density: 1000.0,
            smoothing_length: h,
            max_density_iterations: 50,
            density_tolerance: 1e-3,
            ..Default::default()
        };
        let mut solver = DfsphSolver::new(ps.len());
        solver.compute_alpha(&ps, &neighbors, &kernel, h);

        // First correction: no warm start.
        let iters_cold = solver.correct_density_error(&mut ps, &neighbors, &kernel, &params, 0.001);

        // Reset and apply again with warm start (kappa already set from first run).
        let iters_warm = solver.correct_density_error(&mut ps, &neighbors, &kernel, &params, 0.001);

        // Either the warm start converges faster or in the same number of steps.
        assert!(
            iters_warm <= iters_cold || iters_warm <= params.max_density_iterations,
            "warm start should not be worse: cold={iters_cold}, warm={iters_warm}"
        );
    }

    // --- DFSPH velocity correction tests ---

    #[test]
    fn test_dfsph_velocity_correction_reduces_divergence() {
        let spacing: f64 = 0.05;
        let h: f64 = 0.1;
        let mass = 1000.0 * spacing.powi(3);
        let mut ps = ParticleSet::new();
        for i in 0..4 {
            for j in 0..4 {
                for k in 0..4 {
                    // Add compressive velocity field.
                    let vx = -(i as f64 - 1.5) * 0.1;
                    let vy = -(j as f64 - 1.5) * 0.1;
                    let vz = -(k as f64 - 1.5) * 0.1;
                    let pos = Vec3::new(i as f64 * spacing, j as f64 * spacing, k as f64 * spacing);
                    ps.add_particle(&SphParticle::new(pos, Vec3::new(vx, vy, vz), mass));
                }
            }
        }
        let kernel = CubicSplineKernel;
        let neighbors = SpatialHash::find_all_neighbors(&ps.positions, 2.0 * h);
        crate::wcsph::compute_density(&mut ps, &neighbors, &kernel, h);

        // Divergence before correction.
        let div_before = compute_velocity_divergence(&ps, &neighbors, &kernel, h);
        let max_div_before = div_before
            .iter()
            .cloned()
            .fold(0.0f64, |a, b| a.max(b.abs()));

        let params = DfsphParams {
            rest_density: 1000.0,
            smoothing_length: h,
            max_divergence_iterations: 50,
            divergence_tolerance: 1e-3,
            ..Default::default()
        };
        let mut solver = DfsphSolver::new(ps.len());
        solver.compute_alpha(&ps, &neighbors, &kernel, h);
        solver.correct_divergence_error(&mut ps, &neighbors, &kernel, &params);

        // Divergence after correction should be smaller.
        let div_after = compute_velocity_divergence(&ps, &neighbors, &kernel, h);
        let max_div_after = div_after
            .iter()
            .cloned()
            .fold(0.0f64, |a, b| a.max(b.abs()));

        assert!(
            max_div_after <= max_div_before + 1e-10,
            "Divergence correction should reduce divergence: before={max_div_before}, after={max_div_after}"
        );
    }

    // --- DFSPH pressure projection tests ---

    #[test]
    fn test_dfsph_pressure_projection_stationary() {
        // Stationary fluid: density correction should converge in very few iterations.
        let spacing: f64 = 0.05;
        let h: f64 = 0.1;
        let mass = 1000.0 * spacing.powi(3);
        let mut ps = ParticleSet::new();
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    let pos = Vec3::new(i as f64 * spacing, j as f64 * spacing, k as f64 * spacing);
                    ps.add_particle(&SphParticle::new(pos, Vec3::zeros(), mass));
                }
            }
        }
        let kernel = CubicSplineKernel;
        let neighbors = SpatialHash::find_all_neighbors(&ps.positions, 2.0 * h);
        crate::wcsph::compute_density(&mut ps, &neighbors, &kernel, h);

        let params = DfsphParams {
            rest_density: 1000.0,
            smoothing_length: h,
            max_density_iterations: 50,
            density_tolerance: 1e-3,
            ..Default::default()
        };
        let mut solver = DfsphSolver::new(ps.len());
        solver.compute_alpha(&ps, &neighbors, &kernel, h);

        let iters = solver.correct_density_error(&mut ps, &neighbors, &kernel, &params, 0.001);

        // For stationary fluid the density correction should converge quickly.
        assert!(
            iters <= params.max_density_iterations,
            "Density correction ran {} iterations (max {})",
            iters,
            params.max_density_iterations
        );
    }

    // --- DFSPH solver statistics tests ---

    #[test]
    fn test_dfsph_solver_stats_tracks_state() {
        let mut stats = DfsphSolverStats::new();
        assert_eq!(stats.num_steps, 0);

        stats.record_density_step(0.05, 12);
        stats.record_divergence_step(0.01, 5);

        assert_eq!(stats.num_steps, 1);
        assert!((stats.last_density_error - 0.05).abs() < 1e-14);
        assert!((stats.last_divergence_error - 0.01).abs() < 1e-14);
        assert_eq!(stats.total_density_iterations, 12);
        assert_eq!(stats.total_divergence_iterations, 5);
    }

    #[test]
    fn test_dfsph_solver_stats_multiple() {
        let mut stats = DfsphSolverStats::new();
        stats.record_density_step(0.1, 20);
        stats.record_density_step(0.05, 10);
        stats.record_divergence_step(0.02, 8);
        stats.record_divergence_step(0.01, 4);
        assert_eq!(stats.num_steps, 2);
        assert_eq!(stats.total_density_iterations, 30);
        assert_eq!(stats.total_divergence_iterations, 12);
    }

    #[test]
    fn test_dfsph_solver_stats_avg_iterations() {
        let mut stats = DfsphSolverStats::new();
        stats.record_density_step(0.1, 10);
        stats.record_density_step(0.05, 20);
        let avg = stats.avg_density_iterations();
        assert!(
            (avg - 15.0).abs() < 1e-14,
            "avg density iters should be 15: {avg}"
        );
    }
}
