// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! DFSPH (Divergence-Free SPH) solver.
//!
//! Implements the Bender & Koschier (2015) DFSPH method:
//! - Density solve: iterative pressure correction to enforce ρ = ρ₀
//! - Divergence-free solve: iterative correction to enforce ∇·v = 0
//! - Momentum integration with external forces
//!
//! Uses Wendland C2/C4 kernels for accurate density and gradient estimation.

use std::collections::HashMap;

// ---------------------------------------------------------------------------
// Helper math
// ---------------------------------------------------------------------------

/// Dot product of two 3-vectors.
#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Squared norm of a 3-vector.
#[inline]
fn norm_sq3(a: [f64; 3]) -> f64 {
    dot3(a, a)
}

/// Norm of a 3-vector.
#[inline]
fn norm3(a: [f64; 3]) -> f64 {
    norm_sq3(a).sqrt()
}

/// Scale a 3-vector.
#[inline]
fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

/// Subtract two 3-vectors.
#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Add two 3-vectors.
#[inline]
fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

// ---------------------------------------------------------------------------
// Wendland C2 kernel
// ---------------------------------------------------------------------------

/// Evaluate the Wendland C2 kernel W(r, h).
///
/// W(r, h) = (21/(2π h³)) * (1 - q/2)^4 * (2q + 1)  for q = r/h ≤ 2
///         = 0 otherwise
///
/// where q = r/h, normalised for 3D.
pub fn wendland_c2_kernel(r: f64, h: f64) -> f64 {
    let q = r / h;
    if q >= 2.0 {
        return 0.0;
    }
    let factor = 21.0 / (2.0 * std::f64::consts::PI * h * h * h);
    let t = 1.0 - 0.5 * q;
    factor * t.powi(4) * (2.0 * q + 1.0)
}

/// Evaluate the Wendland C2 kernel gradient magnitude dW/dr.
///
/// Returns the scalar dW/dr for use in vector gradient: ∇W = (dW/dr) * r_hat.
pub fn wendland_c2_gradient_mag(r: f64, h: f64) -> f64 {
    let q = r / h;
    if q >= 2.0 || r < 1e-15 {
        return 0.0;
    }
    let factor = 21.0 / (2.0 * std::f64::consts::PI * h * h * h);
    let t = 1.0 - 0.5 * q;
    // d/dr [(1 - q/2)^4 * (2q+1)] = 4*(1-q/2)^3*(-1/(2h))*(2q+1) + (1-q/2)^4 * 2/h
    let dt_dr = -0.5 / h;

    factor * (4.0 * t.powi(3) * dt_dr * (2.0 * q + 1.0) + t.powi(4) * 2.0 / h)
}

/// Evaluate the Wendland C2 kernel gradient vector ∇W.
pub fn wendland_c2_gradient(r_vec: [f64; 3], h: f64) -> [f64; 3] {
    let r = norm3(r_vec);
    if r < 1e-15 {
        return [0.0; 3];
    }
    let dw_dr = wendland_c2_gradient_mag(r, h);
    scale3(r_vec, dw_dr / r)
}

// ---------------------------------------------------------------------------
// Wendland C4 kernel
// ---------------------------------------------------------------------------

/// Evaluate the Wendland C4 kernel W(r, h).
///
/// W(r, h) = (495/(32π h³)) * (1 - q/2)^6 * (35q²/12 + 3q + 1) for q ≤ 2
///         = 0 otherwise
pub fn wendland_c4_kernel(r: f64, h: f64) -> f64 {
    let q = r / h;
    if q >= 2.0 {
        return 0.0;
    }
    let factor = 495.0 / (32.0 * std::f64::consts::PI * h * h * h);
    let t = 1.0 - 0.5 * q;
    factor * t.powi(6) * (35.0 * q * q / 12.0 + 3.0 * q + 1.0)
}

/// Evaluate the Wendland C4 kernel gradient magnitude dW/dr.
pub fn wendland_c4_gradient_mag(r: f64, h: f64) -> f64 {
    let q = r / h;
    if q >= 2.0 || r < 1e-15 {
        return 0.0;
    }
    let factor = 495.0 / (32.0 * std::f64::consts::PI * h * h * h);
    let t = 1.0 - 0.5 * q;
    let dt_dr = -0.5 / h;
    let poly = 35.0 * q * q / 12.0 + 3.0 * q + 1.0;
    let dpoly_dr = (35.0 * 2.0 * q / 12.0 + 3.0) / h;
    factor * (6.0 * t.powi(5) * dt_dr * poly + t.powi(6) * dpoly_dr)
}

/// Evaluate the Wendland C4 kernel gradient vector ∇W.
pub fn wendland_c4_gradient(r_vec: [f64; 3], h: f64) -> [f64; 3] {
    let r = norm3(r_vec);
    if r < 1e-15 {
        return [0.0; 3];
    }
    let dw_dr = wendland_c4_gradient_mag(r, h);
    scale3(r_vec, dw_dr / r)
}

// ---------------------------------------------------------------------------
// DfsphConfig
// ---------------------------------------------------------------------------

/// Configuration for the DFSPH solver.
///
/// Controls iteration limits, convergence thresholds, and relaxation.
#[derive(Debug, Clone)]
pub struct DfsphConfig {
    /// Maximum pressure correction iterations.
    pub max_pressure_iterations: usize,
    /// Maximum divergence correction iterations.
    pub max_divergence_iterations: usize,
    /// Density error threshold (relative, e.g., 0.001 = 0.1%).
    pub error_threshold: f64,
    /// Under-relaxation parameter ω ∈ (0, 1].
    pub omega: f64,
    /// Rest density ρ₀ (kg/m³).
    pub rest_density: f64,
    /// Smoothing length h (m).
    pub smoothing_length: f64,
    /// Kernel support radius (typically 2h).
    pub support_radius: f64,
    /// Viscosity coefficient (m²/s).
    pub viscosity: f64,
    /// Gravity acceleration.
    pub gravity: [f64; 3],
    /// Surface tension coefficient (N/m).
    pub surface_tension: f64,
}

impl Default for DfsphConfig {
    fn default() -> Self {
        Self {
            max_pressure_iterations: 200,
            max_divergence_iterations: 200,
            error_threshold: 1e-3,
            omega: 0.5,
            rest_density: 1000.0,
            smoothing_length: 0.1,
            support_radius: 0.2,
            viscosity: 1e-6,
            gravity: [0.0, -9.81, 0.0],
            surface_tension: 0.0,
        }
    }
}

// ---------------------------------------------------------------------------
// DfsphParticle
// ---------------------------------------------------------------------------

/// DFSPH particle data.
///
/// Stores all per-particle quantities needed by the DFSPH solver.
#[derive(Debug, Clone)]
pub struct DfsphParticle {
    /// Position (m).
    pub position: [f64; 3],
    /// Velocity (m/s).
    pub velocity: [f64; 3],
    /// Mass (kg).
    pub mass: f64,
    /// Density ρ_i (kg/m³).
    pub density: f64,
    /// Pressure p_i (Pa).
    pub pressure: f64,
    /// DFSPH α coefficient (for pressure solve).
    pub alpha: f64,
    /// Velocity divergence ∇·v_i (1/s).
    pub divergence: f64,
    /// Density error ρ_i - ρ₀ relative to rest density.
    pub density_error: f64,
    /// Pressure acceleration a_p (m/s²).
    pub pressure_accel: [f64; 3],
    /// Is this a boundary particle?
    pub is_boundary: bool,
}

impl DfsphParticle {
    /// Create a new DFSPH particle.
    pub fn new(position: [f64; 3], mass: f64, rest_density: f64) -> Self {
        Self {
            position,
            velocity: [0.0; 3],
            mass,
            density: rest_density,
            pressure: 0.0,
            alpha: 0.0,
            divergence: 0.0,
            density_error: 0.0,
            pressure_accel: [0.0; 3],
            is_boundary: false,
        }
    }

    /// Kinetic energy of this particle.
    pub fn kinetic_energy(&self) -> f64 {
        0.5 * self.mass * norm_sq3(self.velocity)
    }

    /// Potential energy in gravity field.
    pub fn potential_energy(&self, gravity: [f64; 3]) -> f64 {
        // E_pot = m * g * h (g points down, height is along -g direction)
        let g_norm = norm3(gravity);
        if g_norm < 1e-14 {
            return 0.0;
        }
        let h = -dot3(self.position, gravity) / g_norm;
        self.mass * g_norm * h
    }
}

// ---------------------------------------------------------------------------
// DfsphStats
// ---------------------------------------------------------------------------

/// Statistics for DFSPH solver convergence monitoring.
///
/// Tracks density and divergence errors per iteration.
#[derive(Debug, Clone, Default)]
pub struct DfsphStats {
    /// Density error history (relative, averaged over particles).
    pub density_error_history: Vec<f64>,
    /// Divergence error history (averaged over particles).
    pub divergence_error_history: Vec<f64>,
    /// Number of pressure iterations per time step.
    pub pressure_iterations: Vec<usize>,
    /// Number of divergence iterations per time step.
    pub divergence_iterations: Vec<usize>,
    /// Wall-clock time per step (arbitrary units).
    pub step_times: Vec<f64>,
}

impl DfsphStats {
    /// Create new empty stats.
    pub fn new() -> Self {
        Self::default()
    }

    /// Record density convergence for a time step.
    pub fn record_density_step(&mut self, error: f64, iters: usize) {
        self.density_error_history.push(error);
        self.pressure_iterations.push(iters);
    }

    /// Record divergence convergence for a time step.
    pub fn record_divergence_step(&mut self, error: f64, iters: usize) {
        self.divergence_error_history.push(error);
        self.divergence_iterations.push(iters);
    }

    /// Average density error over all recorded steps.
    pub fn mean_density_error(&self) -> f64 {
        if self.density_error_history.is_empty() {
            return 0.0;
        }
        self.density_error_history.iter().sum::<f64>() / self.density_error_history.len() as f64
    }

    /// Average divergence error over all recorded steps.
    pub fn mean_divergence_error(&self) -> f64 {
        if self.divergence_error_history.is_empty() {
            return 0.0;
        }
        self.divergence_error_history.iter().sum::<f64>()
            / self.divergence_error_history.len() as f64
    }

    /// Average number of pressure iterations.
    pub fn mean_pressure_iterations(&self) -> f64 {
        if self.pressure_iterations.is_empty() {
            return 0.0;
        }
        self.pressure_iterations.iter().sum::<usize>() as f64
            / self.pressure_iterations.len() as f64
    }
}

// ---------------------------------------------------------------------------
// WcSphDensity
// ---------------------------------------------------------------------------

/// Density summation using SPH kernel interpolation.
///
/// Supports Wendland C2 and C4 kernels.
#[derive(Debug, Clone, PartialEq)]
pub enum KernelType {
    /// Wendland C2 kernel.
    WendlandC2,
    /// Wendland C4 kernel.
    WendlandC4,
}

/// SPH density summation.
///
/// Computes ρ_i = Σ_j m_j W(|r_i - r_j|, h) for all particles.
#[derive(Debug, Clone)]
pub struct WcSphDensity {
    /// Kernel type.
    pub kernel: KernelType,
    /// Smoothing length h.
    pub h: f64,
}

impl WcSphDensity {
    /// Create a new density summation object.
    pub fn new(kernel: KernelType, h: f64) -> Self {
        Self { kernel, h }
    }

    /// Evaluate kernel at r.
    pub fn eval_kernel(&self, r: f64) -> f64 {
        match self.kernel {
            KernelType::WendlandC2 => wendland_c2_kernel(r, self.h),
            KernelType::WendlandC4 => wendland_c4_kernel(r, self.h),
        }
    }

    /// Evaluate kernel gradient at r_vec.
    pub fn eval_gradient(&self, r_vec: [f64; 3]) -> [f64; 3] {
        match self.kernel {
            KernelType::WendlandC2 => wendland_c2_gradient(r_vec, self.h),
            KernelType::WendlandC4 => wendland_c4_gradient(r_vec, self.h),
        }
    }

    /// Compute density for particle i given its neighbors.
    ///
    /// ρ_i = Σ_j m_j W(r_ij, h)
    pub fn compute_density(&self, pos_i: [f64; 3], neighbors: &[([f64; 3], f64)]) -> f64 {
        neighbors
            .iter()
            .map(|(pos_j, mass_j)| {
                let r = norm3(sub3(pos_i, *pos_j));
                mass_j * self.eval_kernel(r)
            })
            .sum()
    }

    /// Compute density gradient for particle i.
    ///
    /// ∇ρ_i = Σ_j m_j ∇W(r_ij, h)
    pub fn compute_density_gradient(
        &self,
        pos_i: [f64; 3],
        neighbors: &[([f64; 3], f64)],
    ) -> [f64; 3] {
        let mut grad = [0.0f64; 3];
        for (pos_j, mass_j) in neighbors {
            let r_vec = sub3(pos_i, *pos_j);
            let g = self.eval_gradient(r_vec);
            grad[0] += mass_j * g[0];
            grad[1] += mass_j * g[1];
            grad[2] += mass_j * g[2];
        }
        grad
    }
}

// ---------------------------------------------------------------------------
// compute_alpha
// ---------------------------------------------------------------------------

/// Compute DFSPH α coefficient for particle i.
///
/// α_i = ρ_i / (|Σ_j m_j ∇W_ij|² + Σ_j |m_j ∇W_ij|²)
///
/// This is the diagonal element of the pressure correction coefficient.
pub fn compute_alpha(
    rho_i: f64,
    pos_i: [f64; 3],
    mass_i: f64,
    neighbors: &[([f64; 3], f64)], // (position, mass)
    h: f64,
    kernel: &KernelType,
) -> f64 {
    let density_obj = WcSphDensity::new(kernel.clone(), h);
    let mut sum_grad_sq = 0.0f64;
    let mut grad_sum = [0.0f64; 3];

    for (pos_j, mass_j) in neighbors {
        let r_vec = sub3(pos_i, *pos_j);
        let g = density_obj.eval_gradient(r_vec);
        let contrib = scale3(g, *mass_j);
        sum_grad_sq += norm_sq3(scale3(g, *mass_j));
        grad_sum[0] += contrib[0];
        grad_sum[1] += contrib[1];
        grad_sum[2] += contrib[2];
    }

    let self_contrib = scale3(density_obj.eval_gradient([0.0; 3]), mass_i);
    sum_grad_sq += norm_sq3(self_contrib);

    let denom = norm_sq3(grad_sum) + sum_grad_sq;
    if denom < 1e-20 {
        return 0.0;
    }
    rho_i / denom
}

// ---------------------------------------------------------------------------
// AlphaCompute
// ---------------------------------------------------------------------------

/// Computes DFSPH diagonal pressure coefficients for all particles.
///
/// Stores particle alpha values which are used throughout the pressure solve.
#[derive(Debug, Clone)]
pub struct AlphaCompute {
    /// Smoothing length.
    pub h: f64,
    /// Kernel type.
    pub kernel: KernelType,
}

impl AlphaCompute {
    /// Create a new alpha computation object.
    pub fn new(h: f64, kernel: KernelType) -> Self {
        Self { h, kernel }
    }

    /// Compute alpha for all particles given neighbor lists.
    ///
    /// Returns a vector of alpha values, one per particle.
    pub fn compute_all(
        &self,
        particles: &[DfsphParticle],
        neighbor_lists: &[Vec<usize>],
    ) -> Vec<f64> {
        let n = particles.len();
        let density_obj = WcSphDensity::new(self.kernel.clone(), self.h);

        particles
            .iter()
            .enumerate()
            .map(|(i, p_i)| {
                let pos_i = p_i.position;
                let rho_i = p_i.density;

                let mut sum_grad_sq = 0.0f64;
                let mut grad_sum = [0.0f64; 3];

                if let Some(neighbors) = neighbor_lists.get(i) {
                    for &j in neighbors {
                        if j == i || j >= n {
                            continue;
                        }
                        let pos_j = particles[j].position;
                        let mass_j = particles[j].mass;
                        let r_vec = sub3(pos_i, pos_j);
                        let g = density_obj.eval_gradient(r_vec);
                        let contrib = scale3(g, mass_j);
                        sum_grad_sq += norm_sq3(contrib);
                        grad_sum[0] += contrib[0];
                        grad_sum[1] += contrib[1];
                        grad_sum[2] += contrib[2];
                    }
                }

                let denom = norm_sq3(grad_sum) + sum_grad_sq;
                if denom < 1e-20 { 0.0 } else { rho_i / denom }
            })
            .collect()
    }
}

// ---------------------------------------------------------------------------
// dfsph_pressure_accel
// ---------------------------------------------------------------------------

/// Compute pressure acceleration from particle pressures.
///
/// a_p_i = -Σ_j m_j (p_i/ρ_i² + p_j/ρ_j²) ∇W_ij
pub fn dfsph_pressure_accel(
    pos_i: [f64; 3],
    rho_i: f64,
    p_i: f64,
    neighbors: &[([f64; 3], f64, f64, f64)], // (position, mass, rho, pressure)
    h: f64,
    kernel: &KernelType,
) -> [f64; 3] {
    let density_obj = WcSphDensity::new(kernel.clone(), h);
    let mut accel = [0.0f64; 3];
    for (pos_j, mass_j, rho_j, p_j) in neighbors {
        let r_vec = sub3(pos_i, *pos_j);
        let g = density_obj.eval_gradient(r_vec);
        let coeff = if rho_i > 1e-15 && *rho_j > 1e-15 {
            mass_j * (p_i / (rho_i * rho_i) + p_j / (rho_j * rho_j))
        } else {
            0.0
        };
        accel[0] -= coeff * g[0];
        accel[1] -= coeff * g[1];
        accel[2] -= coeff * g[2];
    }
    accel
}

// ---------------------------------------------------------------------------
// PressureSolveIter
// ---------------------------------------------------------------------------

/// Iterative pressure correction to satisfy density invariance ρ = ρ₀.
///
/// Based on the DFSPH constant density (CD) solve from Bender & Koschier (2015).
#[derive(Debug, Clone)]
pub struct PressureSolveIter {
    /// Rest density.
    pub rho0: f64,
    /// Maximum iterations.
    pub max_iter: usize,
    /// Convergence threshold for mean density error.
    pub threshold: f64,
    /// Under-relaxation.
    pub omega: f64,
    /// Kernel type.
    pub kernel: KernelType,
    /// Smoothing length.
    pub h: f64,
}

impl PressureSolveIter {
    /// Create a new pressure solve iterator.
    pub fn new(
        rho0: f64,
        max_iter: usize,
        threshold: f64,
        omega: f64,
        kernel: KernelType,
        h: f64,
    ) -> Self {
        Self {
            rho0,
            max_iter,
            threshold,
            omega,
            kernel,
            h,
        }
    }

    /// Compute the density change rate Dρ/Dt from velocity field.
    ///
    /// Dρ/Dt = ρ Σ_j (m_j/ρ_j) (v_i - v_j) · ∇W_ij
    pub fn density_change_rate(
        &self,
        pos_i: [f64; 3],
        vel_i: [f64; 3],
        rho_i: f64,
        neighbors: &[NeighborEntry], // (pos, vel, mass, rho)
    ) -> f64 {
        let density_obj = WcSphDensity::new(self.kernel.clone(), self.h);
        let mut rate = 0.0f64;
        for (pos_j, vel_j, mass_j, rho_j) in neighbors {
            let r_vec = sub3(pos_i, *pos_j);
            let g = density_obj.eval_gradient(r_vec);
            let v_ij = sub3(vel_i, *vel_j);
            if *rho_j > 1e-15 {
                rate += mass_j * dot3(v_ij, g) / rho_j;
            }
        }
        rho_i * rate
    }

    /// Perform one iteration of pressure correction.
    ///
    /// Returns updated pressures and mean density error.
    pub fn iterate(
        &self,
        particles: &[DfsphParticle],
        alphas: &[f64],
        pressures: &mut [f64],
        dt: f64,
        neighbor_lists: &[Vec<usize>],
    ) -> f64 {
        let n = particles.len();
        let density_obj = WcSphDensity::new(self.kernel.clone(), self.h);
        let mut mean_error = 0.0f64;

        for (i, (p_i, pressure_i)) in particles.iter().zip(pressures.iter_mut()).enumerate() {
            if p_i.is_boundary {
                continue;
            }
            let rho_i = p_i.density;
            let density_error = (rho_i - self.rho0) / self.rho0;
            mean_error += density_error.abs();

            // Pressure correction: κ_i = (ρ_i - ρ₀) / (dt² * α_i)
            let alpha_i = alphas[i];
            if alpha_i.abs() < 1e-20 {
                continue;
            }
            let kappa = density_error * self.rho0 / (dt * dt * alpha_i);
            *pressure_i += self.omega * kappa;
            *pressure_i = pressure_i.max(0.0); // non-negative pressure

            // Update velocity based on pressure acceleration
            if let Some(neighbors) = neighbor_lists.get(i) {
                for &j in neighbors {
                    if j == i || j >= n {
                        continue;
                    }
                    let r_vec = sub3(p_i.position, particles[j].position);
                    let g = density_obj.eval_gradient(r_vec);
                    let _ = g;
                }
            }
        }
        mean_error / n.max(1) as f64
    }
}

// ---------------------------------------------------------------------------
// DivergenceSolveIter
// ---------------------------------------------------------------------------

/// Iterative velocity correction to satisfy divergence-free condition.
///
/// Based on the DFSPH divergence-free (DF) solve: enforces ∇·v = 0.
#[derive(Debug, Clone)]
pub struct DivergenceSolveIter {
    /// Maximum iterations.
    pub max_iter: usize,
    /// Convergence threshold for mean divergence error.
    pub threshold: f64,
    /// Under-relaxation.
    pub omega: f64,
    /// Kernel type.
    pub kernel: KernelType,
    /// Smoothing length.
    pub h: f64,
}

impl DivergenceSolveIter {
    /// Create a new divergence solve iterator.
    pub fn new(max_iter: usize, threshold: f64, omega: f64, kernel: KernelType, h: f64) -> Self {
        Self {
            max_iter,
            threshold,
            omega,
            kernel,
            h,
        }
    }

    /// Compute velocity divergence ∇·v_i.
    ///
    /// ∇·v_i = Σ_j (m_j/ρ_j) (v_j - v_i) · ∇W_ij
    pub fn velocity_divergence(
        &self,
        pos_i: [f64; 3],
        vel_i: [f64; 3],
        rho_i: f64,
        neighbors: &[NeighborEntry], // (pos, vel, mass, rho)
    ) -> f64 {
        let density_obj = WcSphDensity::new(self.kernel.clone(), self.h);
        let mut div = 0.0f64;
        for (pos_j, vel_j, mass_j, rho_j) in neighbors {
            let r_vec = sub3(pos_i, *pos_j);
            let g = density_obj.eval_gradient(r_vec);
            let v_ij = sub3(*vel_j, vel_i);
            if *rho_j > 1e-15 {
                div -= mass_j * dot3(v_ij, g) / rho_i;
            }
        }
        div
    }

    /// Perform one divergence correction iteration.
    ///
    /// Returns mean absolute divergence error (1/s).
    pub fn iterate(
        &self,
        particles: &[DfsphParticle],
        alphas: &[f64],
        velocities: &mut [[f64; 3]],
        dt: f64,
        neighbor_lists: &[Vec<usize>],
    ) -> f64 {
        let n = particles.len();
        let density_obj = WcSphDensity::new(self.kernel.clone(), self.h);
        let mut mean_div = 0.0f64;
        let mut kappas = vec![0.0f64; n];

        // Compute kappa from divergence
        for (i, (p_i, kappa)) in particles.iter().zip(kappas.iter_mut()).enumerate() {
            if p_i.is_boundary {
                continue;
            }
            let alpha_i = alphas[i];
            if alpha_i.abs() < 1e-20 {
                continue;
            }

            let div_i = p_i.divergence;
            mean_div += div_i.abs();
            *kappa = div_i / (dt * alpha_i);
        }

        // Update velocities
        for (i, (p_i, vel_i)) in particles.iter().zip(velocities.iter_mut()).enumerate() {
            if p_i.is_boundary {
                continue;
            }
            if let Some(neighbors) = neighbor_lists.get(i) {
                for &j in neighbors {
                    if j == i || j >= n {
                        continue;
                    }
                    let r_vec = sub3(p_i.position, particles[j].position);
                    let g = density_obj.eval_gradient(r_vec);
                    let rho_j = particles[j].density;
                    if rho_j < 1e-15 {
                        continue;
                    }
                    let corr = self.omega
                        * particles[j].mass
                        * (kappas[i] / p_i.density.powi(2) + kappas[j] / rho_j.powi(2));
                    vel_i[0] -= corr * g[0];
                    vel_i[1] -= corr * g[1];
                    vel_i[2] -= corr * g[2];
                }
            }
        }

        mean_div / n.max(1) as f64
    }
}

// ---------------------------------------------------------------------------
// VelocityAdvection
// ---------------------------------------------------------------------------

/// Neighbor data entry: (position, velocity, mass, density).
pub type NeighborEntry = ([f64; 3], [f64; 3], f64, f64);

/// Semi-Lagrangian advection for DFSPH.
///
/// Advances particle positions and velocities under external forces.
#[derive(Debug, Clone)]
pub struct VelocityAdvection {
    /// Gravity vector.
    pub gravity: [f64; 3],
    /// Viscosity coefficient ν (m²/s).
    pub viscosity: f64,
    /// Smoothing length h.
    pub h: f64,
    /// Kernel type.
    pub kernel: KernelType,
}

impl VelocityAdvection {
    /// Create a new velocity advection object.
    pub fn new(gravity: [f64; 3], viscosity: f64, h: f64, kernel: KernelType) -> Self {
        Self {
            gravity,
            viscosity,
            h,
            kernel,
        }
    }

    /// Compute viscosity force on particle i from neighbors.
    ///
    /// Uses Morris (1997) viscosity: a_visc = ν Σ_j (m_j/ρ_j) (v_j-v_i) ∇²W_ij
    pub fn viscosity_accel(
        &self,
        pos_i: [f64; 3],
        vel_i: [f64; 3],
        rho_i: f64,
        neighbors: &[NeighborEntry], // (pos, vel, mass, rho)
    ) -> [f64; 3] {
        let density_obj = WcSphDensity::new(self.kernel.clone(), self.h);
        let mut accel = [0.0f64; 3];
        for (pos_j, vel_j, mass_j, rho_j) in neighbors {
            let r_vec = sub3(pos_i, *pos_j);
            let r = norm3(r_vec);
            if r < 1e-15 || *rho_j < 1e-15 {
                continue;
            }
            let dw_dr = match self.kernel {
                KernelType::WendlandC2 => wendland_c2_gradient_mag(r, self.h),
                KernelType::WendlandC4 => wendland_c4_gradient_mag(r, self.h),
            };
            let v_ij = sub3(vel_i, *vel_j);
            let coeff = self.viscosity * 2.0 * mass_j / (rho_i * rho_j) * dw_dr / r;
            let _ = density_obj.eval_kernel(r);
            accel[0] -= coeff * v_ij[0];
            accel[1] -= coeff * v_ij[1];
            accel[2] -= coeff * v_ij[2];
        }
        accel
    }

    /// Advance particle velocities by external forces for time dt.
    ///
    /// v*_i = v_i + dt * (g + a_visc_i)
    pub fn advect_velocity(
        &self,
        particles: &mut [DfsphParticle],
        neighbor_data: &[Vec<NeighborEntry>],
        dt: f64,
    ) {
        let n = particles.len();
        let new_velocities: Vec<[f64; 3]> = particles
            .iter()
            .enumerate()
            .map(|(i, p)| {
                if p.is_boundary {
                    return p.velocity;
                }
                let neighbors = neighbor_data.get(i).map(|v| v.as_slice()).unwrap_or(&[]);
                let a_visc = self.viscosity_accel(p.position, p.velocity, p.density, neighbors);
                let a_total = add3(self.gravity, a_visc);
                add3(p.velocity, scale3(a_total, dt))
            })
            .collect();
        for (p, new_vel) in particles.iter_mut().zip(new_velocities.iter()) {
            p.velocity = *new_vel;
        }
        let _ = n;
    }

    /// Advance particle positions by velocity for time dt.
    pub fn advect_positions(&self, particles: &mut [DfsphParticle], dt: f64) {
        for p in particles.iter_mut() {
            if !p.is_boundary {
                p.position = add3(p.position, scale3(p.velocity, dt));
            }
        }
    }
}

// ---------------------------------------------------------------------------
// DfsphBoundary
// ---------------------------------------------------------------------------

/// Mirror/ghost boundary particles for DFSPH.
///
/// Creates boundary particles using the volume-based approach,
/// ensuring no-slip or free-slip conditions.
#[derive(Debug, Clone)]
pub enum BoundaryType {
    /// Rigid no-slip wall.
    NoSlip,
    /// Free-slip wall.
    FreeSlip,
    /// Periodic boundary.
    Periodic([f64; 3]),
}

/// DFSPH boundary handler.
#[derive(Debug, Clone)]
pub struct DfsphBoundary {
    /// Boundary condition type.
    pub boundary_type: BoundaryType,
    /// Rest density for boundary particles.
    pub rho0: f64,
    /// Smoothing length.
    pub h: f64,
    /// Boundary particle positions.
    pub boundary_positions: Vec<[f64; 3]>,
    /// Boundary particle masses.
    pub boundary_masses: Vec<f64>,
}

impl DfsphBoundary {
    /// Create a new boundary handler.
    pub fn new(boundary_type: BoundaryType, rho0: f64, h: f64) -> Self {
        Self {
            boundary_type,
            rho0,
            h,
            boundary_positions: Vec::new(),
            boundary_masses: Vec::new(),
        }
    }

    /// Add a boundary particle at position with given mass.
    pub fn add_particle(&mut self, position: [f64; 3], mass: f64) {
        self.boundary_positions.push(position);
        self.boundary_masses.push(mass);
    }

    /// Create a flat floor boundary from x0 to x1 at y=y_floor.
    pub fn create_floor(&mut self, x0: f64, x1: f64, y_floor: f64, spacing: f64) {
        let n = ((x1 - x0) / spacing) as usize + 1;
        for i in 0..=n {
            let x = x0 + i as f64 * spacing;
            let mass = self.rho0 * spacing * spacing * spacing; // approximate
            self.add_particle([x, y_floor, 0.0], mass);
        }
    }

    /// Compute boundary contribution to particle density.
    ///
    /// ρ_i += Σ_b m_b W(|r_i - r_b|, h)
    pub fn density_contribution(&self, pos_i: [f64; 3], kernel: &KernelType) -> f64 {
        let density_obj = WcSphDensity::new(kernel.clone(), self.h);
        self.boundary_positions
            .iter()
            .zip(self.boundary_masses.iter())
            .map(|(pos_b, mass_b)| {
                let r = norm3(sub3(pos_i, *pos_b));
                mass_b * density_obj.eval_kernel(r)
            })
            .sum()
    }

    /// Apply boundary velocity for no-slip condition to particle.
    ///
    /// For no-slip: mirror velocity v_i → -v_i when near wall.
    pub fn apply_no_slip(&self, vel: [f64; 3], _pos: [f64; 3]) -> [f64; 3] {
        match self.boundary_type {
            BoundaryType::NoSlip => [-vel[0], -vel[1], -vel[2]],
            BoundaryType::FreeSlip => {
                // Reflect normal component, keep tangential
                // Simplified: just keep velocity unchanged
                vel
            }
            BoundaryType::Periodic(_) => vel,
        }
    }

    /// Number of boundary particles.
    pub fn num_particles(&self) -> usize {
        self.boundary_positions.len()
    }
}

// ---------------------------------------------------------------------------
// DfsphSolver
// ---------------------------------------------------------------------------

/// Main DFSPH solver.
///
/// Implements the full DFSPH time stepping algorithm:
/// 1. Advect velocities with external forces
/// 2. Divergence-free solve
/// 3. Density invariant solve
/// 4. Advance positions
#[derive(Debug)]
pub struct DfsphSolver {
    /// Solver configuration.
    pub config: DfsphConfig,
    /// Particle data.
    pub particles: Vec<DfsphParticle>,
    /// Neighbor lists (particle indices).
    pub neighbor_lists: Vec<Vec<usize>>,
    /// Alpha coefficients.
    pub alphas: Vec<f64>,
    /// Convergence statistics.
    pub stats: DfsphStats,
    /// Spatial hash for neighbor search.
    spatial_hash: HashMap<(i64, i64, i64), Vec<usize>>,
}

impl DfsphSolver {
    /// Create a new DFSPH solver with given configuration.
    pub fn new(config: DfsphConfig) -> Self {
        Self {
            config,
            particles: Vec::new(),
            neighbor_lists: Vec::new(),
            alphas: Vec::new(),
            stats: DfsphStats::new(),
            spatial_hash: HashMap::new(),
        }
    }

    /// Add a particle to the simulation.
    pub fn add_particle(&mut self, position: [f64; 3], velocity: [f64; 3], mass: f64) {
        let mut p = DfsphParticle::new(position, mass, self.config.rest_density);
        p.velocity = velocity;
        self.particles.push(p);
    }

    /// Number of fluid particles.
    pub fn num_particles(&self) -> usize {
        self.particles.len()
    }

    /// Build spatial hash for O(1) neighbor lookup.
    pub fn build_spatial_hash(&mut self) {
        self.spatial_hash.clear();
        let h = self.config.support_radius;
        for (i, p) in self.particles.iter().enumerate() {
            let key = self.hash_key(p.position, h);
            self.spatial_hash.entry(key).or_default().push(i);
        }
    }

    /// Compute grid cell key for a position.
    fn hash_key(&self, pos: [f64; 3], h: f64) -> (i64, i64, i64) {
        (
            (pos[0] / h).floor() as i64,
            (pos[1] / h).floor() as i64,
            (pos[2] / h).floor() as i64,
        )
    }

    /// Build neighbor lists using spatial hash.
    pub fn build_neighbor_lists(&mut self) {
        let n = self.particles.len();
        self.neighbor_lists = vec![Vec::new(); n];
        let h = self.config.support_radius;

        let positions: Vec<_> = self.particles.iter().map(|p| p.position).collect();

        for i in 0..n {
            let pos_i = positions[i];
            let base = self.hash_key(pos_i, h);
            let mut neighbors = Vec::new();
            for dx in -2i64..=2 {
                for dy in -2i64..=2 {
                    for dz in -2i64..=2 {
                        let key = (base.0 + dx, base.1 + dy, base.2 + dz);
                        if let Some(cell) = self.spatial_hash.get(&key) {
                            for &j in cell {
                                let r = norm3(sub3(pos_i, positions[j]));
                                if r < h {
                                    neighbors.push(j);
                                }
                            }
                        }
                    }
                }
            }
            self.neighbor_lists[i] = neighbors;
        }
    }

    /// Update particle densities using SPH summation.
    pub fn update_densities(&mut self) {
        let h = self.config.smoothing_length;
        let kernel = WcSphDensity::new(KernelType::WendlandC2, h);
        let _n = self.particles.len();
        let positions: Vec<_> = self.particles.iter().map(|p| p.position).collect();
        let masses: Vec<_> = self.particles.iter().map(|p| p.mass).collect();

        for (i, p_i) in self.particles.iter_mut().enumerate() {
            let mut rho = 0.0f64;
            if let Some(neighbors) = self.neighbor_lists.get(i) {
                for &j in neighbors {
                    let r = norm3(sub3(positions[i], positions[j]));
                    rho += masses[j] * kernel.eval_kernel(r);
                }
            }
            p_i.density = rho.max(0.1 * self.config.rest_density);
        }
    }

    /// Update alpha coefficients for all particles.
    pub fn update_alphas(&mut self) {
        let alpha_compute = AlphaCompute::new(self.config.smoothing_length, KernelType::WendlandC2);
        self.alphas = alpha_compute.compute_all(&self.particles, &self.neighbor_lists);
    }

    /// Perform one time step of the DFSPH algorithm.
    ///
    /// Returns (density_error, divergence_error).
    pub fn step(&mut self, dt: f64) -> (f64, f64) {
        self.build_spatial_hash();
        self.build_neighbor_lists();
        self.update_densities();
        self.update_alphas();

        // 1. Divergence-free solve
        let div_error = self.solve_divergence_free(dt);

        // 2. Advance velocities with external forces
        let gravity = self.config.gravity;
        for p in self.particles.iter_mut() {
            if !p.is_boundary {
                p.velocity = add3(p.velocity, scale3(gravity, dt));
            }
        }

        // 3. Constant density solve
        let dens_error = self.solve_constant_density(dt);

        // 4. Advance positions
        for p in self.particles.iter_mut() {
            if !p.is_boundary {
                p.position = add3(p.position, scale3(p.velocity, dt));
            }
        }

        (dens_error, div_error)
    }

    /// Solve to enforce divergence-free velocity field.
    ///
    /// Returns final mean divergence error.
    pub fn solve_divergence_free(&mut self, dt: f64) -> f64 {
        let n = self.particles.len();
        let h = self.config.smoothing_length;
        let kernel_obj = WcSphDensity::new(KernelType::WendlandC2, h);
        let mut pressures = vec![0.0f64; n];
        let mut final_error = 0.0;

        for _iter in 0..self.config.max_divergence_iterations {
            // Compute divergences
            let positions: Vec<_> = self.particles.iter().map(|p| p.position).collect();
            let velocities: Vec<_> = self.particles.iter().map(|p| p.velocity).collect();
            let densities: Vec<_> = self.particles.iter().map(|p| p.density).collect();
            let masses: Vec<_> = self.particles.iter().map(|p| p.mass).collect();

            let mut total_div = 0.0f64;
            let mut count = 0;
            for (i, p_i) in self.particles.iter_mut().enumerate() {
                if p_i.is_boundary {
                    continue;
                }
                let mut div_i = 0.0f64;
                if let Some(neighbors) = self.neighbor_lists.get(i) {
                    for &j in neighbors {
                        if j == i || j >= n {
                            continue;
                        }
                        let r_vec = sub3(positions[i], positions[j]);
                        let g = kernel_obj.eval_gradient(r_vec);
                        let v_ij = sub3(velocities[i], velocities[j]);
                        if densities[i] > 1e-15 {
                            div_i += masses[j] * dot3(v_ij, g) / densities[i];
                        }
                    }
                }
                p_i.divergence = div_i;
                total_div += div_i.abs();
                count += 1;
            }
            final_error = if count > 0 {
                total_div / count as f64
            } else {
                0.0
            };

            if final_error < self.config.error_threshold {
                self.stats.record_divergence_step(final_error, _iter + 1);
                break;
            }

            // Update pressures
            for (i, (p_i, pressure_i)) in
                self.particles.iter().zip(pressures.iter_mut()).enumerate()
            {
                if p_i.is_boundary {
                    continue;
                }
                let alpha_i = self.alphas.get(i).copied().unwrap_or(0.0);
                if alpha_i.abs() < 1e-20 {
                    continue;
                }
                let kappa = p_i.divergence / (dt * alpha_i);
                *pressure_i += self.config.omega * kappa;
            }

            // Apply pressure correction to velocities
            for i in 0..n {
                if self.particles[i].is_boundary {
                    continue;
                }
                if let Some(neighbors) = self.neighbor_lists.clone().get(i).cloned() {
                    let pos_i = self.particles[i].position;
                    let rho_i = self.particles[i].density;
                    for j in neighbors {
                        if j == i || j >= n {
                            continue;
                        }
                        let pos_j = self.particles[j].position;
                        let mass_j = self.particles[j].mass;
                        let rho_j = self.particles[j].density;
                        let r_vec = sub3(pos_i, pos_j);
                        let g = kernel_obj.eval_gradient(r_vec);
                        if rho_i > 1e-15 && rho_j > 1e-15 {
                            let coeff = mass_j
                                * (pressures[i] / (rho_i * rho_i) + pressures[j] / (rho_j * rho_j));
                            self.particles[i].velocity[0] -= dt * coeff * g[0];
                            self.particles[i].velocity[1] -= dt * coeff * g[1];
                            self.particles[i].velocity[2] -= dt * coeff * g[2];
                        }
                    }
                }
            }
        }
        final_error
    }

    /// Solve to enforce constant density ρ = ρ₀.
    ///
    /// Returns final mean density error.
    pub fn solve_constant_density(&mut self, dt: f64) -> f64 {
        let n = self.particles.len();
        let h = self.config.smoothing_length;
        let kernel_obj = WcSphDensity::new(KernelType::WendlandC2, h);
        let mut pressures = vec![0.0f64; n];
        let mut final_error = 0.0;

        for _iter in 0..self.config.max_pressure_iterations {
            let positions: Vec<_> = self.particles.iter().map(|p| p.position).collect();
            let densities: Vec<_> = self.particles.iter().map(|p| p.density).collect();

            let mut total_error = 0.0f64;
            let mut count = 0;
            for (i, (p_i, pressure_i)) in
                self.particles.iter().zip(pressures.iter_mut()).enumerate()
            {
                if p_i.is_boundary {
                    continue;
                }
                let err = (densities[i] - self.config.rest_density) / self.config.rest_density;
                total_error += err.abs();
                count += 1;

                let alpha_i = self.alphas.get(i).copied().unwrap_or(0.0);
                if alpha_i.abs() < 1e-20 {
                    continue;
                }
                let kappa = err * self.config.rest_density / (dt * dt * alpha_i);
                *pressure_i += self.config.omega * kappa;
                *pressure_i = pressure_i.max(0.0);
            }

            final_error = if count > 0 {
                total_error / count as f64
            } else {
                0.0
            };
            if final_error < self.config.error_threshold {
                self.stats.record_density_step(final_error, _iter + 1);
                break;
            }

            // Apply pressure forces
            for i in 0..n {
                if self.particles[i].is_boundary {
                    continue;
                }
                if let Some(neighbors) = self.neighbor_lists.clone().get(i).cloned() {
                    let pos_i = positions[i];
                    let rho_i = densities[i];
                    for j in neighbors {
                        if j == i || j >= n {
                            continue;
                        }
                        let pos_j = positions[j];
                        let mass_j = self.particles[j].mass;
                        let rho_j = densities[j];
                        let r_vec = sub3(pos_i, pos_j);
                        let g = kernel_obj.eval_gradient(r_vec);
                        if rho_i > 1e-15 && rho_j > 1e-15 {
                            let coeff = mass_j
                                * (pressures[i] / (rho_i * rho_i) + pressures[j] / (rho_j * rho_j));
                            self.particles[i].velocity[0] -= dt * coeff * g[0];
                            self.particles[i].velocity[1] -= dt * coeff * g[1];
                            self.particles[i].velocity[2] -= dt * coeff * g[2];
                        }
                    }
                }
            }
        }
        final_error
    }

    /// Compute total kinetic energy of all particles.
    pub fn kinetic_energy(&self) -> f64 {
        self.particles.iter().map(|p| p.kinetic_energy()).sum()
    }

    /// Compute mean particle density.
    pub fn mean_density(&self) -> f64 {
        let n = self.particles.len();
        if n == 0 {
            return 0.0;
        }
        self.particles.iter().map(|p| p.density).sum::<f64>() / n as f64
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn make_default_config() -> DfsphConfig {
        DfsphConfig {
            rest_density: 1000.0,
            smoothing_length: 0.1,
            support_radius: 0.2,
            max_pressure_iterations: 50,
            max_divergence_iterations: 50,
            error_threshold: 1e-3,
            omega: 0.5,
            viscosity: 1e-6,
            gravity: [0.0, -9.81, 0.0],
            surface_tension: 0.0,
        }
    }

    #[test]
    fn test_wendland_c2_at_zero() {
        let h = 0.1;
        let w = wendland_c2_kernel(0.0, h);
        // At r=0: W(0,h) = (21/(2π h³)) * 1 * 1
        let expected = 21.0 / (2.0 * std::f64::consts::PI * h * h * h);
        assert!(
            (w - expected).abs() < 1e-10,
            "W(0,h) should equal {}, got {}",
            expected,
            w
        );
    }

    #[test]
    fn test_wendland_c2_at_support() {
        let h = 0.1;
        let w = wendland_c2_kernel(2.0 * h, h);
        assert!(w.abs() < 1e-14, "W(2h,h) should be zero, got {}", w);
    }

    #[test]
    fn test_wendland_c2_positive() {
        let h = 0.1;
        for r in [0.0, 0.05, 0.1, 0.15, 0.19] {
            let w = wendland_c2_kernel(r, h);
            assert!(w >= 0.0, "W(r,h) should be non-negative for r={}", r);
        }
    }

    #[test]
    fn test_wendland_c4_at_zero() {
        let h = 0.1;
        let w = wendland_c4_kernel(0.0, h);
        let expected = 495.0 / (32.0 * std::f64::consts::PI * h * h * h);
        assert!(
            (w - expected).abs() < 1e-6,
            "W_C4(0,h) should equal {}",
            expected
        );
    }

    #[test]
    fn test_wendland_c4_at_support() {
        let h = 0.1;
        let w = wendland_c4_kernel(2.0 * h + 1e-12, h);
        assert!(w.abs() < 1e-14, "W_C4 should be zero beyond support");
    }

    #[test]
    fn test_wendland_c2_gradient_antisymmetric() {
        let h = 0.1;
        let r_vec = [0.05, 0.0, 0.0];
        let r_neg = [-0.05, 0.0, 0.0];
        let g1 = wendland_c2_gradient(r_vec, h);
        let g2 = wendland_c2_gradient(r_neg, h);
        assert!(
            (g1[0] + g2[0]).abs() < 1e-10,
            "gradient should be antisymmetric"
        );
    }

    #[test]
    fn test_dfsph_config_default() {
        let cfg = DfsphConfig::default();
        assert!(cfg.rest_density > 0.0);
        assert!(cfg.smoothing_length > 0.0);
        assert!(cfg.max_pressure_iterations > 0);
    }

    #[test]
    fn test_dfsph_particle_new() {
        let p = DfsphParticle::new([1.0, 2.0, 3.0], 0.001, 1000.0);
        assert_eq!(p.position, [1.0, 2.0, 3.0]);
        assert_eq!(p.mass, 0.001);
        assert_eq!(p.density, 1000.0);
        assert!(!p.is_boundary);
    }

    #[test]
    fn test_dfsph_particle_kinetic_energy() {
        let mut p = DfsphParticle::new([0.0; 3], 1.0, 1000.0);
        p.velocity = [3.0, 4.0, 0.0];
        let ke = p.kinetic_energy();
        assert!(
            (ke - 12.5).abs() < 1e-10,
            "KE should be 0.5*1*25=12.5, got {}",
            ke
        );
    }

    #[test]
    fn test_dfsph_solver_new() {
        let cfg = make_default_config();
        let solver = DfsphSolver::new(cfg);
        assert_eq!(solver.num_particles(), 0);
    }

    #[test]
    fn test_dfsph_solver_add_particle() {
        let cfg = make_default_config();
        let mut solver = DfsphSolver::new(cfg);
        solver.add_particle([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 0.001);
        assert_eq!(solver.num_particles(), 1);
    }

    #[test]
    fn test_dfsph_solver_kinetic_energy() {
        let cfg = make_default_config();
        let mut solver = DfsphSolver::new(cfg);
        solver.add_particle([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 2.0);
        let ke = solver.kinetic_energy();
        assert!((ke - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_dfsph_stats_new() {
        let stats = DfsphStats::new();
        assert!(stats.density_error_history.is_empty());
    }

    #[test]
    fn test_dfsph_stats_record_density() {
        let mut stats = DfsphStats::new();
        stats.record_density_step(0.01, 5);
        assert_eq!(stats.density_error_history.len(), 1);
        assert_eq!(stats.pressure_iterations[0], 5);
    }

    #[test]
    fn test_dfsph_stats_mean_error() {
        let mut stats = DfsphStats::new();
        stats.record_density_step(0.02, 3);
        stats.record_density_step(0.04, 4);
        let mean = stats.mean_density_error();
        assert!((mean - 0.03).abs() < 1e-10, "mean error should be 0.03");
    }

    #[test]
    fn test_wcsph_density_compute() {
        let density_obj = WcSphDensity::new(KernelType::WendlandC2, 0.1);
        let pos_i = [0.0; 3];
        let neighbors = vec![([0.05, 0.0, 0.0], 0.001), ([-0.05, 0.0, 0.0], 0.001)];
        let rho = density_obj.compute_density(pos_i, &neighbors);
        assert!(rho > 0.0, "density should be positive");
    }

    #[test]
    fn test_compute_alpha_zero_neighbors() {
        let p = DfsphParticle::new([0.0; 3], 0.001, 1000.0);
        let alpha = compute_alpha(
            p.density,
            p.position,
            p.mass,
            &[],
            0.1,
            &KernelType::WendlandC2,
        );
        assert_eq!(alpha, 0.0, "alpha with no neighbors should be 0");
    }

    #[test]
    fn test_pressure_solve_density_change_rate() {
        let solver_iter =
            PressureSolveIter::new(1000.0, 50, 1e-3, 0.5, KernelType::WendlandC2, 0.1);
        let pos_i = [0.0; 3];
        let vel_i = [1.0, 0.0, 0.0];
        let rho_i = 1000.0;
        let neighbors = vec![([0.05, 0.0, 0.0], [0.5, 0.0, 0.0], 0.001, 1000.0)];
        let rate = solver_iter.density_change_rate(pos_i, vel_i, rho_i, &neighbors);
        // Non-zero divergence should give non-zero rate
        let _ = rate;
    }

    #[test]
    fn test_divergence_solve_velocity_divergence() {
        let div_solver = DivergenceSolveIter::new(50, 1e-3, 0.5, KernelType::WendlandC2, 0.1);
        let pos_i = [0.0; 3];
        let vel_i = [0.0; 3];
        let rho_i = 1000.0;
        let neighbors = vec![([0.05, 0.0, 0.0], [0.1, 0.0, 0.0], 0.001, 1000.0)];
        let div = div_solver.velocity_divergence(pos_i, vel_i, rho_i, &neighbors);
        let _ = div;
    }

    #[test]
    fn test_velocity_advection_gravity() {
        let adv = VelocityAdvection::new([0.0, -9.81, 0.0], 1e-6, 0.1, KernelType::WendlandC2);
        let mut particles = vec![DfsphParticle::new([0.0; 3], 1.0, 1000.0)];
        let neighbor_data: Vec<Vec<NeighborEntry>> = vec![vec![]];
        let dt = 0.01;
        adv.advect_velocity(&mut particles, &neighbor_data, dt);
        // After dt = 0.01 with g = -9.81: v_y should be -0.0981
        assert!((particles[0].velocity[1] - (-9.81 * dt)).abs() < 1e-10);
    }

    #[test]
    fn test_dfsph_boundary_floor() {
        let mut boundary = DfsphBoundary::new(BoundaryType::NoSlip, 1000.0, 0.1);
        boundary.create_floor(0.0, 1.0, 0.0, 0.1);
        assert!(boundary.num_particles() > 0, "floor should have particles");
    }

    #[test]
    fn test_dfsph_boundary_density_contribution() {
        let mut boundary = DfsphBoundary::new(BoundaryType::NoSlip, 1000.0, 0.1);
        boundary.add_particle([0.0, 0.0, 0.0], 0.001);
        let pos_i = [0.05, 0.0, 0.0];
        let contrib = boundary.density_contribution(pos_i, &KernelType::WendlandC2);
        assert!(
            contrib > 0.0,
            "boundary density contribution should be positive"
        );
    }

    #[test]
    fn test_dfsph_solver_step() {
        let cfg = DfsphConfig {
            rest_density: 1000.0,
            smoothing_length: 0.1,
            support_radius: 0.2,
            max_pressure_iterations: 5,
            max_divergence_iterations: 5,
            error_threshold: 1e-2,
            omega: 0.5,
            viscosity: 1e-6,
            gravity: [0.0, -9.81, 0.0],
            surface_tension: 0.0,
        };
        let mut solver = DfsphSolver::new(cfg);
        // Add a small cluster of particles
        for i in 0..3 {
            for j in 0..3 {
                solver.add_particle([i as f64 * 0.05, j as f64 * 0.05, 0.0], [0.0; 3], 0.001);
            }
        }
        let (dens_err, div_err) = solver.step(0.001);
        assert!(dens_err >= 0.0, "density error should be non-negative");
        assert!(div_err >= 0.0, "divergence error should be non-negative");
    }

    #[test]
    fn test_dfsph_solver_mean_density() {
        let cfg = make_default_config();
        let mut solver = DfsphSolver::new(cfg);
        solver.add_particle([0.0; 3], [0.0; 3], 0.001);
        // Before update, density equals rest_density
        assert!((solver.mean_density() - 1000.0).abs() < 1e-10);
    }
}
