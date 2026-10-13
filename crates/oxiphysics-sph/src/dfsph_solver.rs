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

use crate::kernel::{SphKernel, WendlandC2Kernel, WendlandC4Kernel};

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
/// W(r, h) = (21/(16π h³)) * (1 - q/2)^4 * (2q + 1)  for q = r/h < 2
///         = 0 otherwise
///
/// The support radius is `2h` and the kernel is normalised in 3-D, i.e.
/// `∫ W dV = 1` over the support sphere. This is the canonical
/// `21/(2π H³)` Wendland constant written for the support radius `H = 2h`;
/// it delegates to [`WendlandC2Kernel`] so there is a single definition of the
/// kernel in the crate.
pub fn wendland_c2_kernel(r: f64, h: f64) -> f64 {
    WendlandC2Kernel.w(r, h)
}

/// Evaluate the Wendland C2 kernel gradient magnitude dW/dr.
///
/// Returns the scalar dW/dr for use in vector gradient: ∇W = (dW/dr) * r_hat.
pub fn wendland_c2_gradient_mag(r: f64, h: f64) -> f64 {
    WendlandC2Kernel.grad_w(r, h)
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
/// W(r, h) = (495/(256π h³)) * (1 - q/2)^6 * (35q²/12 + 3q + 1) for q = r/h < 2
///         = 0 otherwise
///
/// The support radius is `2h` and the kernel is normalised in 3-D
/// (`∫ W dV = 1`); it delegates to [`WendlandC4Kernel`].
pub fn wendland_c4_kernel(r: f64, h: f64) -> f64 {
    WendlandC4Kernel.w(r, h)
}

/// Evaluate the Wendland C4 kernel gradient magnitude dW/dr.
pub fn wendland_c4_gradient_mag(r: f64, h: f64) -> f64 {
    WendlandC4Kernel.grad_w(r, h)
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
    /// Density-change rate over density `Dρ_i/Dt / ρ_i` (1/s): the negative of
    /// the velocity divergence `∇·v_i`, positive where the particle is being
    /// compressed.
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

    /// Perform one iteration of the DFSPH constant-density solve.
    ///
    /// This is the canonical, self-correcting DFSPH constant-density (CD) step
    /// of Bender & Koschier (2015), matching the reference
    /// `Dfsph::correct_density_error`:
    ///
    /// 1. Predict the density from the **current** velocity field,
    ///    `ρ*_i = ρ_i + dt · Σ_j m_j (v_i − v_j) · ∇W_ij`.
    /// 2. Form the per-particle stiffness (a pressure, in Pa)
    ///    `κ_i = ω · max(ρ*_i − ρ₀, 0) · ρ_i · α_i / dt²` (written into
    ///    `pressures` so the caller can inspect it; non-negative). With
    ///    `α_i = ρ_i / (|Σ_j m_j ∇W_ij|² + Σ_j |m_j ∇W_ij|²)` the factor
    ///    `ρ_i α_i / dt²` is the inverse of the diagonal response of `ρ*_i` to its
    ///    own pressure, so a full step (`ω = 1`) removes the particle's own
    ///    density error to first order.
    /// 3. Project the velocity field with the symmetric pressure acceleration
    ///    `v_i -= dt · Σ_j m_j (κ_i/ρ_i² + κ_j/ρ_j²) ∇W_ij`.
    ///
    /// Because the predicted density is re-evaluated from the corrected
    /// velocities each call, repeated application drives the density error to
    /// `ρ₀` and the iteration is stable (it does not accumulate unbounded
    /// pressure). This mirrors the divergence-free sibling
    /// [`DivergenceSolveIter::iterate`].
    ///
    /// Returns the mean compressive relative predicted-density error
    /// `mean_i max(ρ*_i − ρ₀, 0) / ρ₀` over the fluid particles, evaluated
    /// before the projection, so a caller can use the return value as a
    /// convergence criterion. Only compression is counted because a positive
    /// pressure cannot repair an under-dense particle (for example one on a
    /// free surface, whose kernel sum is truncated); counting it would keep the
    /// criterion above any tolerance forever.
    pub fn iterate(
        &self,
        particles: &[DfsphParticle],
        alphas: &[f64],
        pressures: &mut [f64],
        velocities: &mut [[f64; 3]],
        dt: f64,
        neighbor_lists: &[Vec<usize>],
    ) -> f64 {
        let n = particles.len();
        let density_obj = WcSphDensity::new(self.kernel.clone(), self.h);
        let mut mean_error = 0.0f64;
        let mut fluid_count = 0usize;

        // 1./2. Predict the density from the live velocity field and turn the
        //       residual density error into the per-particle stiffness κ.
        for (i, (p_i, pressure_i)) in particles.iter().zip(pressures.iter_mut()).enumerate() {
            if p_i.is_boundary {
                *pressure_i = 0.0;
                continue;
            }
            let rho_i = p_i.density;

            // Dρ/Dt from the corrected velocities: ρ Σ_j m_j (v_i − v_j)·∇W_ij.
            let mut drho_dt = 0.0f64;
            if let Some(neighbors) = neighbor_lists.get(i) {
                for &j in neighbors {
                    if j == i || j >= n {
                        continue;
                    }
                    let p_j = &particles[j];
                    let rho_j = p_j.density;
                    if rho_j <= 1e-15 {
                        continue;
                    }
                    let r_vec = sub3(p_i.position, p_j.position);
                    let g = density_obj.eval_gradient(r_vec);
                    let v_ij = sub3(velocities[i], velocities[j]);
                    drho_dt += p_j.mass * dot3(v_ij, g);
                }
            }

            // Predicted density and its compressive error against rest density.
            let rho_star = rho_i + dt * drho_dt;
            let compression = (rho_star - self.rho0).max(0.0);
            mean_error += compression / self.rho0;
            fluid_count += 1;

            // Stiffness: κ_i = ω max(ρ*_i − ρ₀, 0) ρ_i α_i / dt²; only compresses
            // (non-negative), so expansion is left to the divergence solve.
            let alpha_i = alphas.get(i).copied().unwrap_or(0.0);
            if alpha_i.abs() < 1e-20 {
                *pressure_i = 0.0;
                continue;
            }
            let kappa = self.omega * compression * rho_i * alpha_i / (dt * dt);
            *pressure_i = kappa.max(0.0);
        }

        // 3. Apply the pressure acceleration to the velocity field so that the
        //    incompressibility constraint is actually enforced. Without this
        //    projection the computed stiffnesses would have no dynamical effect.
        let kappas: &[f64] = pressures;
        for (i, (p_i, vel_i)) in particles.iter().zip(velocities.iter_mut()).enumerate() {
            if p_i.is_boundary {
                continue;
            }
            let rho_i = p_i.density;
            if rho_i <= 1e-15 {
                continue;
            }
            let kappa_i = kappas.get(i).copied().unwrap_or(0.0);
            if let Some(neighbors) = neighbor_lists.get(i) {
                for &j in neighbors {
                    if j == i || j >= n {
                        continue;
                    }
                    let p_j = &particles[j];
                    let rho_j = p_j.density;
                    if rho_j <= 1e-15 {
                        continue;
                    }
                    let r_vec = sub3(p_i.position, p_j.position);
                    let g = density_obj.eval_gradient(r_vec);
                    let kappa_j = kappas.get(j).copied().unwrap_or(0.0);
                    let corr = p_j.mass * (kappa_i / (rho_i * rho_i) + kappa_j / (rho_j * rho_j));
                    vel_i[0] -= dt * corr * g[0];
                    vel_i[1] -= dt * corr * g[1];
                    vel_i[2] -= dt * corr * g[2];
                }
            }
        }

        mean_error / fluid_count.max(1) as f64
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

    /// Compute the density-change rate over density of particle `i`,
    /// `Dρ_i/Dt / ρ_i = (1/ρ_i) Σ_j m_j (v_i − v_j) · ∇W_ij`.
    ///
    /// This is the **negative** of the velocity divergence `∇·v_i`: it is
    /// positive where the particle is being compressed and negative where it is
    /// expanding, the sign convention [`iterate`](Self::iterate) works with.
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
    /// This is the DFSPH divergence-free (DF) step of Bender & Koschier
    /// (2015, 2017):
    ///
    /// 1. Evaluate the density-change rate of the **current** `velocities`,
    ///    `div_i = (1/ρ_i) Σ_j m_j (v_i − v_j) · ∇W_ij` (`= Dρ_i/Dt / ρ_i`,
    ///    positive in compression, i.e. the negative of `∇·v`). The
    ///    `divergence` field of the particles is not read, so repeated calls
    ///    re-measure what the previous call left behind instead of re-applying a
    ///    stale correction.
    /// 2. Form the per-particle stiffness (a pressure, in Pa)
    ///    `κ_i = max(div_i, 0) · ρ_i² · α_i / dt`, with
    ///    `α_i = ρ_i / (|Σ_j m_j ∇W_ij|² + Σ_j |m_j ∇W_ij|²)`. The factor
    ///    `ρ_i² α_i / dt` inverts the diagonal response of `Dρ_i/Dt` to its own
    ///    pressure. Only compression is corrected, as in the reference DFSPH
    ///    solver: a positive pressure cannot repair an expanding particle, and
    ///    a negative one would pull particles together (tensile clumping); the
    ///    two-sided diagonally preconditioned Jacobi iteration is also unstable
    ///    on compact clusters (it diverges already at `ω = 0.3` on a compressed
    ///    5×5×5 block).
    /// 3. Project the velocities with the ω-relaxed symmetric pressure
    ///    acceleration `v_i -= ω · dt · Σ_j m_j (κ_i/ρ_i² + κ_j/ρ_j²) ∇W_ij`.
    ///
    /// Returns the mean compressive divergence `mean_i max(div_i, 0)` (1/s) of
    /// the fluid particles, evaluated before the projection.
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
        let mut fluid_count = 0usize;
        let mut kappas = vec![0.0f64; n];

        // 1./2. Divergence of the live velocity field and the stiffness κ it implies.
        for (i, (p_i, kappa)) in particles.iter().zip(kappas.iter_mut()).enumerate() {
            if p_i.is_boundary {
                continue;
            }
            let rho_i = p_i.density;
            if rho_i <= 1e-15 {
                continue;
            }
            let mut drho_dt = 0.0f64;
            if let Some(neighbors) = neighbor_lists.get(i) {
                for &j in neighbors {
                    if j == i || j >= n {
                        continue;
                    }
                    let p_j = &particles[j];
                    if p_j.density <= 1e-15 {
                        continue;
                    }
                    let g = density_obj.eval_gradient(sub3(p_i.position, p_j.position));
                    let v_ij = sub3(velocities[i], velocities[j]);
                    drho_dt += p_j.mass * dot3(v_ij, g);
                }
            }
            let compression = (drho_dt / rho_i).max(0.0);
            mean_div += compression;
            fluid_count += 1;

            let alpha_i = alphas.get(i).copied().unwrap_or(0.0);
            if alpha_i.abs() < 1e-20 {
                continue;
            }
            *kappa = compression * rho_i * rho_i * alpha_i / dt;
        }

        // 3. Project the velocities with the symmetric pressure acceleration.
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
                        * dt
                        * particles[j].mass
                        * (kappas[i] / p_i.density.powi(2) + kappas[j] / rho_j.powi(2));
                    vel_i[0] -= corr * g[0];
                    vel_i[1] -= corr * g[1];
                    vel_i[2] -= corr * g[2];
                }
            }
        }

        mean_div / fluid_count.max(1) as f64
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
    ///
    /// [`build_spatial_hash`](Self::build_spatial_hash) must have been called
    /// for the current positions. The hash cell size equals the support radius,
    /// so every particle closer than the support radius lies in the 3×3×3 block
    /// of cells around the query cell. Each list contains the particle itself
    /// (distance zero), which [`update_densities`](Self::update_densities) uses
    /// for the self contribution; the solver loops skip `j == i`. A particle
    /// with non-finite coordinates has an empty list (no distance to it is
    /// below the support radius) and does not disturb the others; cell indices
    /// that would overflow `i64` for such a particle are skipped.
    pub fn build_neighbor_lists(&mut self) {
        let n = self.particles.len();
        self.neighbor_lists = vec![Vec::new(); n];
        let h = self.config.support_radius;
        let support_sq = h * h;

        let positions: Vec<_> = self.particles.iter().map(|p| p.position).collect();

        for i in 0..n {
            let pos_i = positions[i];
            let base = self.hash_key(pos_i, h);
            let mut neighbors = Vec::new();
            for dx in -1i64..=1 {
                for dy in -1i64..=1 {
                    for dz in -1i64..=1 {
                        let (Some(kx), Some(ky), Some(kz)) = (
                            base.0.checked_add(dx),
                            base.1.checked_add(dy),
                            base.2.checked_add(dz),
                        ) else {
                            continue;
                        };
                        if let Some(cell) = self.spatial_hash.get(&(kx, ky, kz)) {
                            for &j in cell {
                                if norm_sq3(sub3(pos_i, positions[j])) < support_sq {
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

    /// Rate of density change `Dρ_i/Dt = Σ_j m_j (v_i − v_j) · ∇W_ij` of
    /// particle `i`, evaluated from the current velocity field.
    fn density_rate(&self, i: usize, kernel_obj: &WcSphDensity) -> f64 {
        let n = self.particles.len();
        let p_i = &self.particles[i];
        let mut rate = 0.0f64;
        if let Some(neighbors) = self.neighbor_lists.get(i) {
            for &j in neighbors {
                if j == i || j >= n {
                    continue;
                }
                let p_j = &self.particles[j];
                let g = kernel_obj.eval_gradient(sub3(p_i.position, p_j.position));
                rate += p_j.mass * dot3(sub3(p_i.velocity, p_j.velocity), g);
            }
        }
        rate
    }

    /// Evaluate the divergence-free constraint on the current velocities.
    ///
    /// Stores `Dρ_i/Dt / ρ_i` (positive in compression) in the `divergence`
    /// field of every fluid particle and returns the per-particle stiffnesses
    /// (pressures, Pa) `κ_i = ω · max(div_i, 0) · ρ_i² · α_i / dt` together
    /// with the mean compressive divergence `mean_i max(div_i, 0)` (1/s) of the
    /// fluid particles. Only compression is corrected (reference DFSPH): a
    /// positive pressure cannot repair an expanding particle, and a negative
    /// one pulls particles together and makes the Jacobi iteration unstable on
    /// compact clusters.
    fn divergence_stiffnesses(&mut self, dt: f64, kernel_obj: &WcSphDensity) -> (Vec<f64>, f64) {
        let n = self.particles.len();
        let omega = self.config.omega;
        let mut kappas = vec![0.0f64; n];
        let mut total = 0.0f64;
        let mut count = 0usize;
        for (i, kappa) in kappas.iter_mut().enumerate() {
            if self.particles[i].is_boundary {
                continue;
            }
            let rho_i = self.particles[i].density;
            let div_i = if rho_i > 1e-15 {
                self.density_rate(i, kernel_obj) / rho_i
            } else {
                0.0
            };
            self.particles[i].divergence = div_i;
            let compression = div_i.max(0.0);
            total += compression;
            count += 1;

            let alpha_i = self.alphas.get(i).copied().unwrap_or(0.0);
            if alpha_i.abs() < 1e-20 {
                continue;
            }
            *kappa = omega * compression * rho_i * rho_i * alpha_i / dt;
        }
        let mean = if count > 0 { total / count as f64 } else { 0.0 };
        (kappas, mean)
    }

    /// Evaluate the constant-density constraint on the current velocities.
    ///
    /// Predicts `ρ*_i = ρ_i + dt · Dρ_i/Dt` and returns the per-particle
    /// stiffnesses (pressures, Pa) `κ_i = ω · max(ρ*_i − ρ₀, 0) · ρ_i · α_i / dt²`
    /// together with the mean compressive relative density error
    /// `mean_i max(ρ*_i − ρ₀, 0) / ρ₀` of the fluid particles. Only compression
    /// is penalised: a positive pressure cannot repair an under-dense (for
    /// example free-surface) particle, so counting it would keep the error above
    /// any tolerance forever.
    fn density_stiffnesses(&self, dt: f64, kernel_obj: &WcSphDensity) -> (Vec<f64>, f64) {
        let n = self.particles.len();
        let omega = self.config.omega;
        let rho0 = self.config.rest_density;
        let mut kappas = vec![0.0f64; n];
        let mut total = 0.0f64;
        let mut count = 0usize;
        for (i, kappa) in kappas.iter_mut().enumerate() {
            if self.particles[i].is_boundary {
                continue;
            }
            let rho_i = self.particles[i].density;
            let rho_star = rho_i + dt * self.density_rate(i, kernel_obj);
            let compression = (rho_star - rho0).max(0.0);
            total += compression / rho0;
            count += 1;

            let alpha_i = self.alphas.get(i).copied().unwrap_or(0.0);
            if alpha_i.abs() < 1e-20 {
                continue;
            }
            *kappa = omega * compression * rho_i * alpha_i / (dt * dt);
        }
        let mean = if count > 0 { total / count as f64 } else { 0.0 };
        (kappas, mean)
    }

    /// Apply the symmetric pressure acceleration of the stiffnesses `kappas`
    /// (Pa) to the fluid velocities:
    /// `v_i -= dt · Σ_j m_j (κ_i/ρ_i² + κ_j/ρ_j²) ∇W_ij`.
    ///
    /// Every correction is computed from the pre-projection stiffnesses into a
    /// scratch buffer and only then applied (Jacobi update), so the result does
    /// not depend on the particle order and no per-particle copy of the
    /// neighbour lists is needed.
    fn apply_stiffness_projection(&mut self, kappas: &[f64], dt: f64, kernel_obj: &WcSphDensity) {
        let n = self.particles.len();
        let mut dv = vec![[0.0f64; 3]; n];
        for (i, dv_i) in dv.iter_mut().enumerate() {
            let p_i = &self.particles[i];
            if p_i.is_boundary || p_i.density <= 1e-15 {
                continue;
            }
            let k_i = kappas[i] / (p_i.density * p_i.density);
            if let Some(neighbors) = self.neighbor_lists.get(i) {
                for &j in neighbors {
                    if j == i || j >= n {
                        continue;
                    }
                    let p_j = &self.particles[j];
                    if p_j.density <= 1e-15 {
                        continue;
                    }
                    let g = kernel_obj.eval_gradient(sub3(p_i.position, p_j.position));
                    let coeff = p_j.mass * (k_i + kappas[j] / (p_j.density * p_j.density));
                    dv_i[0] -= dt * coeff * g[0];
                    dv_i[1] -= dt * coeff * g[1];
                    dv_i[2] -= dt * coeff * g[2];
                }
            }
        }
        for (p, d) in self.particles.iter_mut().zip(dv.iter()) {
            if !p.is_boundary {
                p.velocity = add3(p.velocity, *d);
            }
        }
    }

    /// Solve to enforce divergence-free velocity field.
    ///
    /// Each iteration measures the density-change rate `div_i = Dρ_i/Dt / ρ_i`
    /// (positive in compression) of the **current** velocities, turns its
    /// compressive part into per-particle stiffnesses
    /// `κ_i = ω · max(div_i, 0) · ρ_i² · α_i / dt` and projects the velocities
    /// with `v_i -= dt · Σ_j m_j (κ_i/ρ_i² + κ_j/ρ_j²) ∇W_ij`. The stiffnesses
    /// are re-derived from the corrected velocities every iteration (they are
    /// assigned, not accumulated), so the iteration converges instead of
    /// re-applying the whole accumulated pressure each time.
    ///
    /// The number of iterations used (the converging pass counts, so
    /// `max_divergence_iterations` when the solve did not converge) and the
    /// final error are recorded in [`DfsphStats`] either way. Returns the final
    /// mean compressive divergence (1/s) of the fluid particles.
    pub fn solve_divergence_free(&mut self, dt: f64) -> f64 {
        let kernel_obj = WcSphDensity::new(KernelType::WendlandC2, self.config.smoothing_length);
        let max_iter = self.config.max_divergence_iterations;
        let threshold = self.config.error_threshold;

        for iter in 0..max_iter {
            let (kappas, error) = self.divergence_stiffnesses(dt, &kernel_obj);
            if error < threshold {
                self.stats.record_divergence_step(error, iter + 1);
                return error;
            }
            self.apply_stiffness_projection(&kappas, dt, &kernel_obj);
        }

        // Not converged: report the residual that is actually left after the
        // last projection rather than the one measured before it.
        let (_, error) = self.divergence_stiffnesses(dt, &kernel_obj);
        self.stats.record_divergence_step(error, max_iter);
        error
    }

    /// Solve to enforce constant density ρ = ρ₀.
    ///
    /// Each iteration predicts the density `ρ*_i = ρ_i + dt · Dρ_i/Dt` from the
    /// **current** velocities, turns the compressive part of the error into
    /// per-particle stiffnesses `κ_i = ω · max(ρ*_i − ρ₀, 0) · ρ_i · α_i / dt²`
    /// and projects the velocities with
    /// `v_i -= dt · Σ_j m_j (κ_i/ρ_i² + κ_j/ρ_j²) ∇W_ij`. The stiffnesses are
    /// assigned, not accumulated, and are re-derived from the corrected
    /// velocities every iteration, so the predicted density really approaches
    /// ρ₀ and the convergence test can fire.
    ///
    /// The number of iterations used (the converging pass counts, so
    /// `max_pressure_iterations` when the solve did not converge) and the final
    /// error are recorded in [`DfsphStats`] either way. Returns the final mean
    /// compressive relative density error of the fluid particles.
    pub fn solve_constant_density(&mut self, dt: f64) -> f64 {
        let kernel_obj = WcSphDensity::new(KernelType::WendlandC2, self.config.smoothing_length);
        let max_iter = self.config.max_pressure_iterations;
        let threshold = self.config.error_threshold;

        for iter in 0..max_iter {
            let (kappas, error) = self.density_stiffnesses(dt, &kernel_obj);
            if error < threshold {
                self.stats.record_density_step(error, iter + 1);
                return error;
            }
            self.apply_stiffness_projection(&kappas, dt, &kernel_obj);
        }

        // Not converged: report the residual that is actually left after the
        // last projection rather than the one measured before it.
        let (_, error) = self.density_stiffnesses(dt, &kernel_obj);
        self.stats.record_density_step(error, max_iter);
        error
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
mod tests;
