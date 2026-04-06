#![allow(clippy::needless_range_loop, clippy::ptr_arg)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Incompressible SPH (ISPH) pressure solver and simulation step.
//!
//! Implements the projection-based ISPH method:
//! 1. Predict velocity using explicit integration (ignoring pressure).
//! 2. Solve the pressure Poisson equation ∇²P = ρ/Δt · ∇·v*.
//! 3. Correct velocity with the pressure gradient.
//!
//! # Included types
//! - [`IsphParticle`] — per-particle state: position, velocity, pressure, intermediate velocity
//! - [`IsphPressureSolver`] — Poisson solver (Jacobi/Gauss-Seidel iteration)
//! - [`IsphStep`] — full ISPH timestep
//! - [`DivergenceFreeCondition`] — check ∇·v = 0 after pressure correction
//! - [`VelocityDivergence`] — compute ∇·v at each particle via SPH summation
//! - [`PressureGradient`] — compute ∇P (symmetric formulation)
//! - [`IsphBoundary`] — mirror particles for solid wall enforcement
//! - [`IsphMultiphase`] — two-fluid ISPH: density ratio and surface tension
//! - [`IsphTurbulence`] — Smagorinsky SGS model for ISPH-LES
//! - [`IsphStats`] — runtime divergence, pressure convergence, CFL tracking

use std::f64::consts::PI;

// ─────────────────────────────────────────────────────────────────────────────
// Simple 3-D vector helpers (no nalgebra)
// ─────────────────────────────────────────────────────────────────────────────

/// A 3-D position or velocity stored as `[f64; 3]`.
pub type Vec3 = [f64; 3];

#[allow(dead_code)]
fn vec3_add(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[allow(dead_code)]
fn vec3_sub(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[allow(dead_code)]
fn vec3_scale(a: Vec3, s: f64) -> Vec3 {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[allow(dead_code)]
fn vec3_dot(a: Vec3, b: Vec3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[allow(dead_code)]
fn vec3_norm(a: Vec3) -> f64 {
    vec3_dot(a, a).sqrt()
}

#[allow(dead_code)]
fn vec3_norm2(a: Vec3) -> f64 {
    vec3_dot(a, a)
}

// ─────────────────────────────────────────────────────────────────────────────
// Kernel (cubic spline, 3-D)
// ─────────────────────────────────────────────────────────────────────────────

/// Cubic spline kernel value W(r, h).
#[allow(dead_code)]
fn kernel_w(r: f64, h: f64) -> f64 {
    let sigma = 1.0 / (PI * h * h * h);
    let q = r / h;
    if q >= 2.0 {
        0.0
    } else if q >= 1.0 {
        let t = 2.0 - q;
        sigma * 0.25 * t * t * t
    } else {
        sigma * (1.0 - 1.5 * q * q + 0.75 * q * q * q)
    }
}

/// Gradient of cubic spline kernel dW/dr (scalar; multiply by r̂ to get vector).
#[allow(dead_code)]
fn kernel_dw_dr(r: f64, h: f64) -> f64 {
    let sigma = 1.0 / (PI * h * h * h);
    let q = r / h;
    if !(1e-14..2.0).contains(&q) {
        0.0
    } else if q >= 1.0 {
        let t = 2.0 - q;
        sigma * (-0.75 * t * t) / h
    } else {
        sigma * (-3.0 * q + 2.25 * q * q) / h
    }
}

/// Gradient vector ∇W(rᵢⱼ, h) where rᵢⱼ = xᵢ − xⱼ.
#[allow(dead_code)]
fn kernel_grad(rij: Vec3, h: f64) -> Vec3 {
    let r = vec3_norm(rij);
    if r < 1e-14 {
        return [0.0; 3];
    }
    let dw = kernel_dw_dr(r, h);
    vec3_scale(rij, dw / r)
}

// ─────────────────────────────────────────────────────────────────────────────
// IsphParticle
// ─────────────────────────────────────────────────────────────────────────────

/// State of a single ISPH fluid particle.
#[derive(Debug, Clone)]
pub struct IsphParticle {
    /// Current position (m).
    pub position: Vec3,
    /// Current velocity (m/s).
    pub velocity: Vec3,
    /// Pressure (Pa).
    pub pressure: f64,
    /// Intermediate (predicted) velocity before pressure correction.
    pub velocity_star: Vec3,
    /// Particle mass (kg).
    pub mass: f64,
    /// Rest density (kg/m³).
    pub rest_density: f64,
    /// Computed density (kg/m³).
    pub density: f64,
    /// Velocity divergence ∇·v (1/s).
    pub div_velocity: f64,
    /// Fluid phase index (0 or 1) for multiphase simulations.
    pub phase: usize,
}

impl IsphParticle {
    /// Create a new ISPH particle at rest.
    pub fn new(position: Vec3, mass: f64, rest_density: f64) -> Self {
        Self {
            position,
            velocity: [0.0; 3],
            pressure: 0.0,
            velocity_star: [0.0; 3],
            mass,
            rest_density,
            density: rest_density,
            div_velocity: 0.0,
            phase: 0,
        }
    }

    /// Kinetic energy ½ m |v|².
    pub fn kinetic_energy(&self) -> f64 {
        0.5 * self.mass * vec3_dot(self.velocity, self.velocity)
    }

    /// Particle volume = mass / density.
    pub fn volume(&self) -> f64 {
        if self.density > 1e-300 {
            self.mass / self.density
        } else {
            0.0
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// VelocityDivergence
// ─────────────────────────────────────────────────────────────────────────────

/// Computes ∇·v at each particle via SPH summation.
///
/// ∇·vᵢ = Σⱼ (mⱼ/ρⱼ) (vⱼ − vᵢ) · ∇Wᵢⱼ
#[allow(dead_code)]
pub struct VelocityDivergence {
    /// Smoothing length h.
    pub h: f64,
}

impl VelocityDivergence {
    /// Create a velocity-divergence computer with smoothing length `h`.
    pub fn new(h: f64) -> Self {
        Self { h }
    }

    /// Compute ∇·vᵢ for all particles in-place.
    pub fn compute(&self, particles: &mut Vec<IsphParticle>) {
        let n = particles.len();
        let mut divs = vec![0.0_f64; n];

        for i in 0..n {
            let xi = particles[i].position;
            let vi = particles[i].velocity;
            let mut div = 0.0_f64;
            for j in 0..n {
                if i == j {
                    continue;
                }
                let xj = particles[j].position;
                let vj = particles[j].velocity;
                let rij = vec3_sub(xi, xj);
                let r = vec3_norm(rij);
                if r >= 2.0 * self.h {
                    continue;
                }
                let grad_w = kernel_grad(rij, self.h);
                let dv = vec3_sub(vj, vi);
                let vdotg = vec3_dot(dv, grad_w);
                let vol_j = particles[j].mass / particles[j].density.max(1e-300);
                div += vol_j * vdotg;
            }
            divs[i] = div;
        }
        for (p, d) in particles.iter_mut().zip(divs.iter()) {
            p.div_velocity = *d;
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PressureGradient
// ─────────────────────────────────────────────────────────────────────────────

/// Computes ∇P at each particle using the symmetric SPH formulation.
///
/// (∇P/ρ)ᵢ = Σⱼ mⱼ (Pᵢ/ρᵢ² + Pⱼ/ρⱼ²) ∇Wᵢⱼ
#[allow(dead_code)]
pub struct PressureGradient {
    /// Smoothing length h.
    pub h: f64,
}

impl PressureGradient {
    /// Create a pressure-gradient computer.
    pub fn new(h: f64) -> Self {
        Self { h }
    }

    /// Compute the pressure-acceleration −∇P/ρ for particle `i`.
    pub fn acceleration(&self, i: usize, particles: &[IsphParticle]) -> Vec3 {
        let xi = particles[i].position;
        let pi = particles[i].pressure;
        let rhoi = particles[i].density.max(1e-300);
        let mut acc = [0.0_f64; 3];

        for (j, pj_part) in particles.iter().enumerate() {
            if i == j {
                continue;
            }
            let xj = pj_part.position;
            let rij = vec3_sub(xi, xj);
            let r = vec3_norm(rij);
            if r >= 2.0 * self.h {
                continue;
            }
            let grad_w = kernel_grad(rij, self.h);
            let pj = pj_part.pressure;
            let rhoj = pj_part.density.max(1e-300);
            let coeff = pj_part.mass * (pi / (rhoi * rhoi) + pj / (rhoj * rhoj));
            acc = vec3_sub(acc, vec3_scale(grad_w, coeff));
        }
        acc
    }

    /// Compute pressure-acceleration for all particles.
    pub fn all_accelerations(&self, particles: &[IsphParticle]) -> Vec<Vec3> {
        (0..particles.len())
            .map(|i| self.acceleration(i, particles))
            .collect()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// IsphPressureSolver
// ─────────────────────────────────────────────────────────────────────────────

/// Parameters for the ISPH pressure Poisson solver.
#[derive(Debug, Clone)]
pub struct IsphPressureSolverParams {
    /// Maximum Jacobi/GS iterations.
    pub max_iter: usize,
    /// Convergence tolerance on relative pressure change.
    pub tolerance: f64,
    /// Relaxation parameter ω ∈ (0, 2).
    pub relaxation: f64,
    /// Rest density (kg/m³).
    pub rest_density: f64,
}

impl Default for IsphPressureSolverParams {
    fn default() -> Self {
        Self {
            max_iter: 100,
            tolerance: 1e-4,
            relaxation: 0.5,
            rest_density: 1000.0,
        }
    }
}

/// Iterative pressure Poisson solver for ISPH.
///
/// Solves ∇²P = ρ/Δt · ∇·v* using a relaxed Jacobi iteration with
/// diagonal approximation of the Laplacian operator.
#[allow(dead_code)]
pub struct IsphPressureSolver {
    /// Solver parameters.
    pub params: IsphPressureSolverParams,
    /// Smoothing length.
    pub h: f64,
}

impl IsphPressureSolver {
    /// Create a new pressure solver.
    pub fn new(h: f64, params: IsphPressureSolverParams) -> Self {
        Self { params, h }
    }

    /// Run one pressure-Poisson solve; updates `p.pressure` for each particle.
    ///
    /// Returns `(iterations, max_pressure_change)`.
    pub fn solve(&self, particles: &mut Vec<IsphParticle>, dt: f64) -> (usize, f64) {
        let n = particles.len();
        let h = self.h;
        let rho0 = self.params.rest_density;
        let omega = self.params.relaxation;

        let mut iters = 0usize;
        let mut max_change = f64::INFINITY;

        while iters < self.params.max_iter && max_change > self.params.tolerance {
            max_change = 0.0_f64;
            let old_pressures: Vec<f64> = particles.iter().map(|p| p.pressure).collect();

            for i in 0..n {
                let xi = particles[i].position;
                let div_vi = particles[i].div_velocity;

                // Diagonal coefficient aᵢᵢ ≈ Σⱼ (mⱼ/ρⱼ) |∇Wᵢⱼ|² × 2/ρᵢ
                let mut a_diag = 0.0_f64;
                let mut rhs = rho0 / dt * div_vi;

                // Off-diagonal contribution to residual
                let mut off_diag_sum = 0.0_f64;
                for j in 0..n {
                    if i == j {
                        continue;
                    }
                    let xj = particles[j].position;
                    let rij = vec3_sub(xi, xj);
                    let r = vec3_norm(rij);
                    if r >= 2.0 * h {
                        continue;
                    }
                    let grad_w = kernel_grad(rij, h);
                    let vol_j = particles[j].mass / particles[j].density.max(1e-300);
                    let g2 = vec3_dot(grad_w, grad_w);
                    let rhoi = particles[i].density.max(1e-300);
                    let rhoj = particles[j].density.max(1e-300);
                    a_diag += vol_j * g2 * (1.0 / rhoi + 1.0 / rhoj);
                    off_diag_sum += vol_j * g2 * old_pressures[j] * (1.0 / rhoi + 1.0 / rhoj);
                }

                rhs -= off_diag_sum;
                let new_p = if a_diag.abs() > 1e-300 {
                    (1.0 - omega) * old_pressures[i] + omega * rhs / a_diag
                } else {
                    old_pressures[i]
                };
                let new_p = new_p.max(0.0); // Non-negative pressure enforcement

                let dp = (new_p - old_pressures[i]).abs();
                if dp > max_change {
                    max_change = dp;
                }
                particles[i].pressure = new_p;
            }
            iters += 1;
        }
        (iters, max_change)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// IsphStep
// ─────────────────────────────────────────────────────────────────────────────

/// Full ISPH timestep: predict → solve pressure → correct.
///
/// 1. Predict: v* = v + Δt (g + f_visc)
/// 2. Divergence: compute ∇·v*
/// 3. Pressure: solve Poisson equation
/// 4. Correct: v = v* − Δt/ρ ∇P,  x = x + Δt v
#[allow(dead_code)]
pub struct IsphStep {
    /// Smoothing length.
    pub h: f64,
    /// Pressure solver.
    pub pressure_solver: IsphPressureSolver,
    /// Kinematic viscosity ν (m²/s).
    pub viscosity: f64,
    /// Gravity vector (m/s²).
    pub gravity: Vec3,
}

impl IsphStep {
    /// Create a new ISPH step integrator.
    pub fn new(h: f64, pressure_solver: IsphPressureSolver, viscosity: f64, gravity: Vec3) -> Self {
        Self {
            h,
            pressure_solver,
            viscosity,
            gravity,
        }
    }

    /// Predict velocity using explicit forces (gravity + viscosity).
    fn predict(&self, particles: &mut Vec<IsphParticle>, dt: f64) {
        let n = particles.len();
        let h = self.h;
        let nu = self.viscosity;
        let g = self.gravity;

        let pos: Vec<Vec3> = particles.iter().map(|p| p.position).collect();
        let vel: Vec<Vec3> = particles.iter().map(|p| p.velocity).collect();
        let mas: Vec<f64> = particles.iter().map(|p| p.mass).collect();
        let rho: Vec<f64> = particles.iter().map(|p| p.density.max(1e-300)).collect();

        for i in 0..n {
            let xi = pos[i];
            let vi = vel[i];
            let rhoi = rho[i];

            // Viscous acceleration via SPH Laplacian approximation
            // a_visc = ν Σⱼ mⱼ/ρⱼ (vⱼ - vᵢ) · ∇²W
            let mut a_visc = [0.0_f64; 3];
            for j in 0..n {
                if i == j {
                    continue;
                }
                let xj = pos[j];
                let rij = vec3_sub(xi, xj);
                let r = vec3_norm(rij);
                if r >= 2.0 * h {
                    continue;
                }
                let dv = vec3_sub(vel[j], vi);
                let dw = kernel_dw_dr(r, h);
                let laplacian_coeff = 2.0 * dw / (r.max(1e-14) * rhoi);
                let vol_j = mas[j] / rho[j];
                a_visc = vec3_add(a_visc, vec3_scale(dv, nu * vol_j * laplacian_coeff));
            }

            let v_star = vec3_add(vec3_add(vi, vec3_scale(g, dt)), vec3_scale(a_visc, dt));
            particles[i].velocity_star = v_star;
        }
    }

    /// Correct velocity using pressure gradient.
    fn correct(&self, particles: &mut Vec<IsphParticle>, dt: f64) {
        let pg = PressureGradient::new(self.h);
        let accs = pg.all_accelerations(particles);
        for (i, p) in particles.iter_mut().enumerate() {
            p.velocity = vec3_add(p.velocity_star, vec3_scale(accs[i], dt));
            p.position = vec3_add(p.position, vec3_scale(p.velocity, dt));
        }
    }

    /// Advance the particle system by one timestep `dt`.
    ///
    /// Returns solver statistics from the pressure solve.
    pub fn advance(&self, particles: &mut Vec<IsphParticle>, dt: f64) -> (usize, f64) {
        self.predict(particles, dt);
        // Copy velocity_star → velocity for divergence computation
        for p in particles.iter_mut() {
            p.velocity = p.velocity_star;
        }
        let div_computer = VelocityDivergence::new(self.h);
        div_computer.compute(particles);
        let stats = self.pressure_solver.solve(particles, dt);
        self.correct(particles, dt);
        stats
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// DivergenceFreeCondition
// ─────────────────────────────────────────────────────────────────────────────

/// Checks and reports the divergence-free condition ∇·v ≈ 0.
#[allow(dead_code)]
pub struct DivergenceFreeCondition {
    /// Tolerance threshold for acceptable divergence.
    pub tolerance: f64,
}

impl DivergenceFreeCondition {
    /// Create a divergence checker with given tolerance.
    pub fn new(tolerance: f64) -> Self {
        Self { tolerance }
    }

    /// Compute the maximum absolute divergence over all particles.
    pub fn max_divergence(&self, particles: &[IsphParticle]) -> f64 {
        particles
            .iter()
            .map(|p| p.div_velocity.abs())
            .fold(0.0_f64, f64::max)
    }

    /// Mean-square divergence error.
    pub fn mean_squared_divergence(&self, particles: &[IsphParticle]) -> f64 {
        if particles.is_empty() {
            return 0.0;
        }
        let sum: f64 = particles
            .iter()
            .map(|p| p.div_velocity * p.div_velocity)
            .sum();
        sum / particles.len() as f64
    }

    /// Returns `true` if the divergence-free condition is satisfied within tolerance.
    pub fn is_satisfied(&self, particles: &[IsphParticle]) -> bool {
        self.max_divergence(particles) <= self.tolerance
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// IsphBoundary
// ─────────────────────────────────────────────────────────────────────────────

/// Solid wall boundary treatment via mirror (ghost) particles.
///
/// Each fluid particle near a wall is mirrored across the wall plane to
/// enforce no-penetration and no-slip conditions.
#[allow(dead_code)]
pub struct IsphBoundary {
    /// Wall plane origin.
    pub plane_point: Vec3,
    /// Outward unit normal of the wall (pointing into the fluid).
    pub plane_normal: Vec3,
    /// Smoothing length used to determine cutoff for mirroring.
    pub h: f64,
}

impl IsphBoundary {
    /// Create a planar boundary.
    pub fn new(plane_point: Vec3, plane_normal: Vec3, h: f64) -> Self {
        // Normalise the normal
        let n = vec3_norm(plane_normal);
        let normal = if n > 1e-14 {
            vec3_scale(plane_normal, 1.0 / n)
        } else {
            [0.0, 1.0, 0.0]
        };
        Self {
            plane_point,
            plane_normal: normal,
            h,
        }
    }

    /// Signed distance of a point from the wall (positive = inside fluid).
    pub fn signed_distance(&self, x: Vec3) -> f64 {
        let d = vec3_sub(x, self.plane_point);
        vec3_dot(d, self.plane_normal)
    }

    /// Mirror position of a fluid particle across the wall.
    pub fn mirror_position(&self, x: Vec3) -> Vec3 {
        let dist = self.signed_distance(x);
        vec3_sub(x, vec3_scale(self.plane_normal, 2.0 * dist))
    }

    /// Mirror velocity (flip normal component for no-slip wall).
    pub fn mirror_velocity(&self, v: Vec3) -> Vec3 {
        let vn = vec3_dot(v, self.plane_normal);
        // Reflect normal component, negate tangential for no-slip
        let v_normal = vec3_scale(self.plane_normal, vn);
        let v_tangent = vec3_sub(v, v_normal);
        // No-slip: tangential also flipped; no-penetration: normal flipped

        vec3_sub(vec3_scale(v_normal, -1.0), v_tangent)
    }

    /// Generate ghost particles for all fluid particles within 2h of the wall.
    ///
    /// Returns a list of ghost [`IsphParticle`] objects.
    pub fn generate_ghosts(&self, particles: &[IsphParticle]) -> Vec<IsphParticle> {
        particles
            .iter()
            .filter(|p| {
                let d = self.signed_distance(p.position);
                d >= 0.0 && d < 2.0 * self.h
            })
            .map(|p| {
                let mut ghost = p.clone();
                ghost.position = self.mirror_position(p.position);
                ghost.velocity = self.mirror_velocity(p.velocity);
                ghost.velocity_star = self.mirror_velocity(p.velocity_star);
                ghost
            })
            .collect()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// IsphMultiphase
// ─────────────────────────────────────────────────────────────────────────────

/// Two-fluid ISPH parameters for immiscible fluids.
#[derive(Debug, Clone)]
pub struct IsphMultiphaseParams {
    /// Density of phase 0 (kg/m³).
    pub density_phase0: f64,
    /// Density of phase 1 (kg/m³).
    pub density_phase1: f64,
    /// Surface tension coefficient σ (N/m).
    pub surface_tension: f64,
    /// Interface thickness parameter (in units of h).
    pub interface_thickness: f64,
}

impl Default for IsphMultiphaseParams {
    fn default() -> Self {
        Self {
            density_phase0: 1000.0,
            density_phase1: 800.0,
            surface_tension: 0.07,
            interface_thickness: 1.5,
        }
    }
}

/// Two-fluid ISPH: handles density ratio and surface tension.
#[allow(dead_code)]
pub struct IsphMultiphase {
    /// Multiphase parameters.
    pub params: IsphMultiphaseParams,
    /// Smoothing length.
    pub h: f64,
}

impl IsphMultiphase {
    /// Create a multiphase ISPH object.
    pub fn new(h: f64, params: IsphMultiphaseParams) -> Self {
        Self { h, params }
    }

    /// Assign rest density to each particle based on its phase.
    pub fn assign_densities(&self, particles: &mut Vec<IsphParticle>) {
        for p in particles.iter_mut() {
            p.rest_density = match p.phase {
                0 => self.params.density_phase0,
                _ => self.params.density_phase1,
            };
            p.density = p.rest_density;
        }
    }

    /// Density ratio ρ₁/ρ₀.
    pub fn density_ratio(&self) -> f64 {
        self.params.density_phase1 / self.params.density_phase0.max(1e-300)
    }

    /// Compute the interface color function cᵢ ∈ \[0, 1\] for each particle
    /// (0 = pure phase 0, 1 = pure phase 1) via SPH interpolation.
    pub fn color_function(&self, particles: &[IsphParticle]) -> Vec<f64> {
        let h = self.h;
        let n = particles.len();
        let mut color = vec![0.0_f64; n];
        for i in 0..n {
            let xi = particles[i].position;
            let ci = particles[i].phase as f64;
            let mut wsum = 0.0_f64;
            let mut csum = 0.0_f64;
            for j in 0..n {
                let xj = particles[j].position;
                let r = vec3_norm(vec3_sub(xi, xj));
                if r >= 2.0 * h {
                    continue;
                }
                let w = kernel_w(r, h);
                csum += w * particles[j].phase as f64;
                wsum += w;
            }
            color[i] = if wsum > 1e-300 {
                (csum / wsum).clamp(0.0, 1.0)
            } else {
                ci
            };
        }
        color
    }

    /// Continuum surface force (CSF) acceleration on particle i.
    ///
    /// Returns the surface-tension-induced acceleration vector.
    pub fn surface_tension_acceleration(&self, i: usize, particles: &[IsphParticle]) -> Vec3 {
        let h = self.h;
        let sigma = self.params.surface_tension;
        let xi = particles[i].position;
        let mut grad_c = [0.0_f64; 3];

        for (j, pj) in particles.iter().enumerate() {
            if i == j {
                continue;
            }
            let rij = vec3_sub(xi, pj.position);
            let r = vec3_norm(rij);
            if r >= 2.0 * h {
                continue;
            }
            let gw = kernel_grad(rij, h);
            let dc = pj.phase as f64 - particles[i].phase as f64;
            grad_c = vec3_add(
                grad_c,
                vec3_scale(gw, dc * pj.mass / pj.density.max(1e-300)),
            );
        }

        let kappa = vec3_norm(grad_c); // Approximate curvature via |∇c|
        let rho_i = particles[i].density.max(1e-300);
        vec3_scale(grad_c, sigma * kappa / rho_i)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// IsphTurbulence
// ─────────────────────────────────────────────────────────────────────────────

/// Smagorinsky SGS turbulence model for ISPH-LES.
///
/// Computes the sub-grid scale (SGS) eddy viscosity:
/// νₜ = (Cₛ Δ)² |S̃|
/// where |S̃| is the magnitude of the resolved strain-rate tensor and
/// Δ is the filter width (≈ h).
#[allow(dead_code)]
pub struct IsphTurbulence {
    /// Smagorinsky constant Cₛ (typically 0.1–0.2).
    pub cs: f64,
    /// Filter width Δ (usually the smoothing length h).
    pub delta: f64,
    /// Smoothing length for SPH operators.
    pub h: f64,
}

impl IsphTurbulence {
    /// Create a Smagorinsky model.
    pub fn new(cs: f64, delta: f64, h: f64) -> Self {
        Self { cs, delta, h }
    }

    /// Compute the strain-rate magnitude |S̃| for particle i via SPH.
    ///
    /// |S̃|² = 2 Sₐβ Sₐβ, where Sₐβ = ½(∂vα/∂xβ + ∂vβ/∂xα).
    pub fn strain_rate_magnitude(&self, i: usize, particles: &[IsphParticle]) -> f64 {
        let h = self.h;
        let xi = particles[i].position;
        let vi = particles[i].velocity;
        // Build velocity gradient tensor G_αβ = ∂vα/∂xβ via SPH
        let mut g = [[0.0_f64; 3]; 3];
        for (j, pj) in particles.iter().enumerate() {
            if i == j {
                continue;
            }
            let rij = vec3_sub(xi, pj.position);
            let r = vec3_norm(rij);
            if r >= 2.0 * h {
                continue;
            }
            let gw = kernel_grad(rij, h);
            let dv = vec3_sub(pj.velocity, vi);
            let vol_j = pj.mass / pj.density.max(1e-300);
            for alpha in 0..3 {
                for beta in 0..3 {
                    g[alpha][beta] += vol_j * dv[alpha] * gw[beta];
                }
            }
        }
        // Symmetric part S = ½(G + Gᵀ)
        let mut s2 = 0.0_f64;
        for alpha in 0..3 {
            for beta in 0..3 {
                let s_ab = 0.5 * (g[alpha][beta] + g[beta][alpha]);
                s2 += 2.0 * s_ab * s_ab;
            }
        }
        s2.sqrt()
    }

    /// SGS eddy viscosity νₜ = (Cₛ Δ)² |S̃| for particle i.
    pub fn eddy_viscosity(&self, i: usize, particles: &[IsphParticle]) -> f64 {
        let s_mag = self.strain_rate_magnitude(i, particles);
        (self.cs * self.delta).powi(2) * s_mag
    }

    /// SGS viscous acceleration for particle i.
    ///
    /// a_sgs = ∇ · (2 νₜ S̃) ≈ Σⱼ mⱼ/ρⱼ νₜ (vⱼ−vᵢ) laplacian_W
    pub fn sgs_acceleration(&self, i: usize, particles: &[IsphParticle]) -> Vec3 {
        let nu_t = self.eddy_viscosity(i, particles);
        let h = self.h;
        let xi = particles[i].position;
        let vi = particles[i].velocity;
        let rhoi = particles[i].density.max(1e-300);
        let mut acc = [0.0_f64; 3];
        for (j, pj) in particles.iter().enumerate() {
            if i == j {
                continue;
            }
            let rij = vec3_sub(xi, pj.position);
            let r = vec3_norm(rij);
            if r >= 2.0 * h {
                continue;
            }
            let dw = kernel_dw_dr(r, h);
            let dv = vec3_sub(pj.velocity, vi);
            let coeff = 2.0 * nu_t * pj.mass / pj.density.max(1e-300) * dw / (r.max(1e-14) * rhoi);
            acc = vec3_add(acc, vec3_scale(dv, coeff));
        }
        acc
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// IsphStats
// ─────────────────────────────────────────────────────────────────────────────

/// Runtime statistics for an ISPH simulation.
#[derive(Debug, Clone, Default)]
pub struct IsphStats {
    /// Maximum divergence error over all particles.
    pub max_divergence: f64,
    /// Mean divergence error.
    pub mean_divergence: f64,
    /// Pressure solver iterations at last timestep.
    pub pressure_iterations: usize,
    /// Final pressure residual at last timestep.
    pub pressure_residual: f64,
    /// CFL number (u_max Δt / h).
    pub cfl: f64,
    /// Simulation time.
    pub time: f64,
    /// Timestep index.
    pub step: usize,
}

impl IsphStats {
    /// Create zero-initialised stats.
    pub fn new() -> Self {
        Self::default()
    }

    /// Update stats from current particle state and solver output.
    pub fn update(
        &mut self,
        particles: &[IsphParticle],
        pressure_iters: usize,
        pressure_residual: f64,
        dt: f64,
        h: f64,
    ) {
        // Divergence
        let n = particles.len();
        if n > 0 {
            let divs: Vec<f64> = particles.iter().map(|p| p.div_velocity.abs()).collect();
            self.max_divergence = divs.iter().cloned().fold(0.0_f64, f64::max);
            self.mean_divergence = divs.iter().sum::<f64>() / n as f64;
        }
        // CFL
        let u_max = particles
            .iter()
            .map(|p| vec3_norm(p.velocity))
            .fold(0.0_f64, f64::max);
        self.cfl = u_max * dt / h.max(1e-300);
        // Pressure solver
        self.pressure_iterations = pressure_iters;
        self.pressure_residual = pressure_residual;
        self.time += dt;
        self.step += 1;
    }

    /// Returns `true` if the CFL condition is satisfied (CFL < 1).
    pub fn cfl_ok(&self) -> bool {
        self.cfl < 1.0
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── Helpers ───────────────────────────────────────────────────────────────

    fn make_particle(pos: Vec3) -> IsphParticle {
        IsphParticle::new(pos, 1.0, 1000.0)
    }

    fn two_particles_separated(sep: f64) -> Vec<IsphParticle> {
        vec![
            make_particle([0.0, 0.0, 0.0]),
            make_particle([sep, 0.0, 0.0]),
        ]
    }

    // ── Vec3 utilities ────────────────────────────────────────────────────────

    #[test]
    fn test_vec3_add() {
        let a = [1.0, 2.0, 3.0];
        let b = [4.0, 5.0, 6.0];
        let c = vec3_add(a, b);
        assert_eq!(c, [5.0, 7.0, 9.0]);
    }

    #[test]
    fn test_vec3_sub() {
        let r = vec3_sub([3.0, 5.0, 7.0], [1.0, 2.0, 3.0]);
        assert_eq!(r, [2.0, 3.0, 4.0]);
    }

    #[test]
    fn test_vec3_dot() {
        assert!((vec3_dot([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]) - 32.0).abs() < 1e-12);
    }

    #[test]
    fn test_vec3_norm() {
        assert!((vec3_norm([3.0, 4.0, 0.0]) - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_vec3_scale() {
        let v = vec3_scale([1.0, 2.0, 3.0], 2.0);
        assert_eq!(v, [2.0, 4.0, 6.0]);
    }

    // ── Kernel ────────────────────────────────────────────────────────────────

    #[test]
    fn test_kernel_w_zero_at_cutoff() {
        assert_eq!(kernel_w(2.1, 1.0), 0.0);
    }

    #[test]
    fn test_kernel_w_positive_inside() {
        assert!(kernel_w(0.5, 1.0) > 0.0);
    }

    #[test]
    fn test_kernel_dw_dr_zero_at_cutoff() {
        assert_eq!(kernel_dw_dr(2.0, 1.0), 0.0);
    }

    #[test]
    fn test_kernel_grad_zero_at_origin() {
        let g = kernel_grad([0.0, 0.0, 0.0], 1.0);
        assert_eq!(g, [0.0, 0.0, 0.0]);
    }

    // ── IsphParticle ──────────────────────────────────────────────────────────

    #[test]
    fn test_particle_new() {
        let p = make_particle([1.0, 2.0, 3.0]);
        assert_eq!(p.position, [1.0, 2.0, 3.0]);
        assert_eq!(p.velocity, [0.0, 0.0, 0.0]);
        assert_eq!(p.pressure, 0.0);
        assert_eq!(p.phase, 0);
    }

    #[test]
    fn test_particle_kinetic_energy_at_rest() {
        let p = make_particle([0.0; 3]);
        assert_eq!(p.kinetic_energy(), 0.0);
    }

    #[test]
    fn test_particle_kinetic_energy_moving() {
        let mut p = make_particle([0.0; 3]);
        p.velocity = [1.0, 0.0, 0.0];
        assert!((p.kinetic_energy() - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_particle_volume() {
        let p = make_particle([0.0; 3]);
        let vol = p.volume();
        assert!((vol - 1.0 / 1000.0).abs() < 1e-12);
    }

    // ── VelocityDivergence ────────────────────────────────────────────────────

    #[test]
    fn test_velocity_divergence_zero_velocity() {
        let h = 0.5;
        let mut particles = two_particles_separated(0.3);
        let div = VelocityDivergence::new(h);
        div.compute(&mut particles);
        for p in &particles {
            assert!(p.div_velocity.is_finite());
        }
    }

    #[test]
    fn test_velocity_divergence_far_apart() {
        let h = 0.1;
        let mut particles = two_particles_separated(10.0);
        for p in particles.iter_mut() {
            p.velocity = [1.0, 0.0, 0.0];
        }
        let div = VelocityDivergence::new(h);
        div.compute(&mut particles);
        // Far apart → no interaction → divergence = 0
        for p in &particles {
            assert!(p.div_velocity.abs() < 1e-10);
        }
    }

    #[test]
    fn test_velocity_divergence_single_particle() {
        let h = 1.0;
        let mut particles = vec![make_particle([0.0; 3])];
        VelocityDivergence::new(h).compute(&mut particles);
        assert_eq!(particles[0].div_velocity, 0.0);
    }

    // ── PressureGradient ──────────────────────────────────────────────────────

    #[test]
    fn test_pressure_gradient_zero_pressure() {
        let particles = two_particles_separated(0.3);
        let pg = PressureGradient::new(0.5);
        let acc = pg.acceleration(0, &particles);
        assert_eq!(acc, [0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_pressure_gradient_far_apart() {
        let particles = two_particles_separated(10.0);
        let pg = PressureGradient::new(0.1);
        let acc = pg.acceleration(0, &particles);
        assert_eq!(acc, [0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_pressure_gradient_all_accelerations_length() {
        let particles = two_particles_separated(0.3);
        let pg = PressureGradient::new(0.5);
        let accs = pg.all_accelerations(&particles);
        assert_eq!(accs.len(), 2);
    }

    // ── IsphPressureSolver ────────────────────────────────────────────────────

    #[test]
    fn test_pressure_solver_default_params() {
        let p = IsphPressureSolverParams::default();
        assert_eq!(p.max_iter, 100);
        assert!(p.tolerance > 0.0);
        assert!(p.relaxation > 0.0 && p.relaxation < 2.0);
    }

    #[test]
    fn test_pressure_solver_zero_divergence() {
        let mut particles = two_particles_separated(0.3);
        let params = IsphPressureSolverParams::default();
        let solver = IsphPressureSolver::new(0.5, params);
        // All divergence = 0 → pressures should stay near 0
        let (iters, _res) = solver.solve(&mut particles, 0.01);
        assert!(iters <= 100);
        for p in &particles {
            assert!(p.pressure >= 0.0);
        }
    }

    #[test]
    fn test_pressure_solver_returns_nonneg_pressure() {
        let mut particles: Vec<IsphParticle> = (0..5)
            .map(|i| make_particle([i as f64 * 0.1, 0.0, 0.0]))
            .collect();
        for p in particles.iter_mut() {
            p.div_velocity = 0.01;
        }
        let params = IsphPressureSolverParams {
            rest_density: 1000.0,
            ..Default::default()
        };
        let solver = IsphPressureSolver::new(0.3, params);
        solver.solve(&mut particles, 0.01);
        for p in &particles {
            assert!(p.pressure >= 0.0, "pressure must be non-negative");
        }
    }

    // ── DivergenceFreeCondition ───────────────────────────────────────────────

    #[test]
    fn test_divergence_free_zero() {
        let particles = vec![make_particle([0.0; 3])];
        let cond = DivergenceFreeCondition::new(1e-4);
        assert!(cond.is_satisfied(&particles));
    }

    #[test]
    fn test_divergence_free_violated() {
        let mut p = make_particle([0.0; 3]);
        p.div_velocity = 1.0;
        let cond = DivergenceFreeCondition::new(1e-4);
        assert!(!cond.is_satisfied(&[p]));
    }

    #[test]
    fn test_divergence_max_and_mean() {
        let mut particles: Vec<IsphParticle> = (0..4).map(|_| make_particle([0.0; 3])).collect();
        particles[0].div_velocity = 0.5;
        particles[1].div_velocity = 0.3;
        let cond = DivergenceFreeCondition::new(1.0);
        assert!((cond.max_divergence(&particles) - 0.5).abs() < 1e-12);
        let msd = cond.mean_squared_divergence(&particles);
        assert!(msd >= 0.0);
    }

    #[test]
    fn test_divergence_free_empty() {
        let cond = DivergenceFreeCondition::new(1e-4);
        assert_eq!(cond.mean_squared_divergence(&[]), 0.0);
    }

    // ── IsphBoundary ──────────────────────────────────────────────────────────

    #[test]
    fn test_boundary_signed_distance() {
        let b = IsphBoundary::new([0.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.1);
        let d = b.signed_distance([0.0, 0.5, 0.0]);
        assert!((d - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_boundary_mirror_position() {
        let b = IsphBoundary::new([0.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.1);
        let m = b.mirror_position([0.0, 0.05, 0.0]);
        assert!((m[1] + 0.05).abs() < 1e-10);
    }

    #[test]
    fn test_boundary_ghost_generation() {
        let b = IsphBoundary::new([0.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.5);
        let mut p = make_particle([0.0, 0.1, 0.0]); // within 2h=1.0 of wall
        p.density = 1000.0;
        let ghosts = b.generate_ghosts(&[p]);
        assert_eq!(ghosts.len(), 1);
        assert!(ghosts[0].position[1] < 0.0);
    }

    #[test]
    fn test_boundary_no_ghosts_far_away() {
        let b = IsphBoundary::new([0.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.1);
        let p = make_particle([0.0, 5.0, 0.0]); // far from wall
        let ghosts = b.generate_ghosts(&[p]);
        assert!(ghosts.is_empty());
    }

    // ── IsphMultiphase ────────────────────────────────────────────────────────

    #[test]
    fn test_multiphase_density_ratio() {
        let params = IsphMultiphaseParams {
            density_phase0: 1000.0,
            density_phase1: 800.0,
            ..Default::default()
        };
        let mp = IsphMultiphase::new(0.1, params);
        assert!((mp.density_ratio() - 0.8).abs() < 1e-10);
    }

    #[test]
    fn test_multiphase_assign_densities() {
        let params = IsphMultiphaseParams::default();
        let mp = IsphMultiphase::new(0.1, params.clone());
        let mut particles = vec![make_particle([0.0; 3]), {
            let mut p = make_particle([0.2, 0.0, 0.0]);
            p.phase = 1;
            p
        }];
        mp.assign_densities(&mut particles);
        assert_eq!(particles[0].density, params.density_phase0);
        assert_eq!(particles[1].density, params.density_phase1);
    }

    #[test]
    fn test_multiphase_color_function_uniform() {
        let mp = IsphMultiphase::new(0.5, IsphMultiphaseParams::default());
        let particles = vec![make_particle([0.0; 3]), make_particle([0.1, 0.0, 0.0])];
        let color = mp.color_function(&particles);
        for &c in &color {
            assert!((0.0..=1.0).contains(&c), "color out of [0,1]: {c}");
        }
    }

    // ── IsphTurbulence ────────────────────────────────────────────────────────

    #[test]
    fn test_turbulence_eddy_viscosity_zero_velocity() {
        let turb = IsphTurbulence::new(0.1, 0.1, 0.1);
        let particles = two_particles_separated(0.05);
        let nu_t = turb.eddy_viscosity(0, &particles);
        assert!(nu_t >= 0.0);
    }

    #[test]
    fn test_turbulence_strain_rate_finite() {
        let turb = IsphTurbulence::new(0.1, 0.1, 0.5);
        let mut particles = two_particles_separated(0.3);
        particles[1].velocity = [1.0, 0.0, 0.0];
        let sr = turb.strain_rate_magnitude(0, &particles);
        assert!(sr.is_finite() && sr >= 0.0);
    }

    #[test]
    fn test_turbulence_sgs_acceleration_finite() {
        let turb = IsphTurbulence::new(0.1, 0.1, 0.5);
        let mut particles: Vec<IsphParticle> = (0..3)
            .map(|i| make_particle([i as f64 * 0.2, 0.0, 0.0]))
            .collect();
        particles[1].velocity = [0.5, 0.0, 0.0];
        let acc = turb.sgs_acceleration(0, &particles);
        for &a in &acc {
            assert!(a.is_finite());
        }
    }

    // ── IsphStats ─────────────────────────────────────────────────────────────

    #[test]
    fn test_stats_default() {
        let s = IsphStats::new();
        assert_eq!(s.step, 0);
        assert_eq!(s.time, 0.0);
        assert!(s.cfl_ok()); // cfl=0 < 1
    }

    #[test]
    fn test_stats_update_step_count() {
        let mut s = IsphStats::new();
        let particles = vec![make_particle([0.0; 3])];
        s.update(&particles, 5, 1e-5, 0.01, 0.1);
        assert_eq!(s.step, 1);
        assert!((s.time - 0.01).abs() < 1e-12);
    }

    #[test]
    fn test_stats_cfl_computed() {
        let mut s = IsphStats::new();
        let mut p = make_particle([0.0; 3]);
        p.velocity = [2.0, 0.0, 0.0]; // |v| = 2
        s.update(&[p], 10, 1e-6, 0.05, 1.0); // CFL = 2 * 0.05 / 1.0 = 0.1
        assert!((s.cfl - 0.1).abs() < 1e-10);
        assert!(s.cfl_ok());
    }

    #[test]
    fn test_stats_cfl_violation() {
        let mut s = IsphStats::new();
        let mut p = make_particle([0.0; 3]);
        p.velocity = [10.0, 0.0, 0.0]; // |v| = 10
        s.update(&[p], 10, 1e-6, 0.2, 0.1); // CFL = 10 * 0.2 / 0.1 = 20
        assert!(!s.cfl_ok());
    }

    #[test]
    fn test_stats_divergence_tracking() {
        let mut s = IsphStats::new();
        let mut p = make_particle([0.0; 3]);
        p.div_velocity = 0.1;
        s.update(&[p], 3, 1e-4, 0.01, 0.1);
        assert!((s.max_divergence - 0.1).abs() < 1e-12);
        assert!((s.mean_divergence - 0.1).abs() < 1e-12);
    }

    // ── IsphStep (smoke test) ─────────────────────────────────────────────────

    #[test]
    fn test_isph_step_smoke() {
        let h = 0.5;
        let params = IsphPressureSolverParams {
            max_iter: 5,
            ..Default::default()
        };
        let solver = IsphPressureSolver::new(h, params);
        let step = IsphStep::new(h, solver, 1e-3, [0.0, -9.81, 0.0]);
        let mut particles: Vec<IsphParticle> = (0..3)
            .map(|i| {
                let mut p = make_particle([i as f64 * 0.2, 0.5, 0.0]);
                p.density = 1000.0;
                p
            })
            .collect();
        let (iters, _) = step.advance(&mut particles, 0.001);
        assert!(iters <= 5);
        for p in &particles {
            assert!(p.position[0].is_finite());
            assert!(p.position[1].is_finite());
        }
    }

    #[test]
    fn test_isph_step_gravity_moves_particles_down() {
        let h = 0.5;
        let params = IsphPressureSolverParams {
            max_iter: 3,
            ..Default::default()
        };
        let solver = IsphPressureSolver::new(h, params);
        let step = IsphStep::new(h, solver, 0.0, [0.0, -9.81, 0.0]);
        let mut particles = vec![{
            let mut p = make_particle([0.0, 1.0, 0.0]);
            p.density = 1000.0;
            p
        }];
        let y0 = particles[0].position[1];
        step.advance(&mut particles, 0.01);
        let y1 = particles[0].position[1];
        assert!(y1 < y0, "particle should move downward under gravity");
    }
}
