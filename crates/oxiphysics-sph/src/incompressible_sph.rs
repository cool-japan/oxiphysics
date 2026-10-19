// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Incompressible SPH (ISPH) pressure solver and simulation step.
//!
//! Implements the projection-based ISPH method:
//! 1. Predict velocity using explicit integration (ignoring pressure).
//! 2. Solve the pressure Poisson equation ∇²P = ρ/Δt · ∇·v* (Brookshaw SPH
//!    Laplacian, free-surface particles pinned to P = 0).
//! 3. Correct velocity with the pressure gradient.
//!
//! # Included types
//! - [`IsphParticle`] — per-particle state: position, velocity, pressure, intermediate velocity
//! - [`IsphPressureSolver`] — Poisson solver (relaxed Gauss–Seidel/SOR iteration)
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

fn vec3_add(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn vec3_sub(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn vec3_scale(a: Vec3, s: f64) -> Vec3 {
    [a[0] * s, a[1] * s, a[2] * s]
}

fn vec3_dot(a: Vec3, b: Vec3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn vec3_norm(a: Vec3) -> f64 {
    vec3_dot(a, a).sqrt()
}

// ─────────────────────────────────────────────────────────────────────────────
// Kernel (cubic spline, 3-D)
// ─────────────────────────────────────────────────────────────────────────────

/// Cubic spline kernel value W(r, h).
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
    pub fn compute(&self, particles: &mut [IsphParticle]) {
        let n = particles.len();
        let mut divs = vec![0.0_f64; n];

        for (i, particle_i) in particles.iter().enumerate() {
            let xi = particle_i.position;
            let vi = particle_i.velocity;
            let mut div = 0.0_f64;
            for (j, pj) in particles.iter().enumerate() {
                if i == j {
                    continue;
                }
                let xj = pj.position;
                let vj = pj.velocity;
                let rij = vec3_sub(xi, xj);
                let r = vec3_norm(rij);
                if r >= 2.0 * self.h {
                    continue;
                }
                let grad_w = kernel_grad(rij, self.h);
                let dv = vec3_sub(vj, vi);
                let vdotg = vec3_dot(dv, grad_w);
                let vol_j = pj.mass / pj.density.max(1e-300);
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
    /// Maximum relaxed Gauss–Seidel sweeps.
    pub max_iter: usize,
    /// Convergence tolerance on the pressure change of one sweep, relative to
    /// the largest pressure magnitude of the field.
    pub tolerance: f64,
    /// Relaxation parameter ω ∈ (0, 2) of the successive over-relaxation
    /// (`1` is plain Gauss–Seidel, `< 1` under-relaxes).
    pub relaxation: f64,
    /// Rest density (kg/m³).
    pub rest_density: f64,
    /// Clamp negative pressures to zero after every sweep.
    ///
    /// `false` (default) solves the pressure Poisson equation exactly, so the
    /// projection removes expansion as well as compression.  `true` suppresses
    /// tension, i.e. the fluid may separate where it would be stretched, which
    /// avoids the tensile (particle-clumping) instability at the price of a
    /// divergence that is left uncorrected in those regions.
    pub clamp_negative_pressure: bool,
}

impl Default for IsphPressureSolverParams {
    fn default() -> Self {
        Self {
            max_iter: 100,
            tolerance: 1e-4,
            relaxation: 0.5,
            rest_density: 1000.0,
            clamp_negative_pressure: false,
        }
    }
}

/// A particle whose SPH position divergence
/// `(∇·r)ᵢ = Σⱼ Vⱼ (xⱼ − xᵢ)·∇ᵢWᵢⱼ = −Σⱼ Vⱼ r dW/dr` falls below this value has a
/// truncated kernel support and is treated as a free-surface particle (Dirichlet
/// condition p = 0).
///
/// A particle with a complete support has `∇·r ≈ 3` (2.8–2.95 on a lattice with
/// `h = 1.2 Δx`, 2.7–3.0 when the particles are displaced randomly by 10 % of the
/// spacing); the outermost layer of such a lattice has at most 2.1 (2.2 with the
/// disorder) and an isolated particle 0.  `2.4 = 0.8 · 3` lies between the two
/// classes, whereas a kernel-sum criterion (a complete support sums to
/// 1 ± 10 % under the same disorder) misclassifies interior particles.
const FREE_SURFACE_POSITION_DIVERGENCE: f64 = 2.4;

/// Brookshaw/Morris SPH Laplacian weight of the neighbour `j` of a particle at
/// distance `r`,
///
/// `cᵢⱼ = 2 Vⱼ (−dW/dr) r / (r² + η²) ≥ 0`,  `η² = 0.01 h²`,
///
/// such that `∇²Aᵢ ≈ Σⱼ cᵢⱼ (Aⱼ − Aᵢ)` for any scalar or vector field `A`.  The
/// weight is positive because `dW/dr < 0`: the Laplacian of a field diffuses the
/// difference `Aⱼ − Aᵢ`.
fn laplacian_weight(vol_j: f64, r: f64, h: f64) -> f64 {
    let eta_sq = 0.01 * h * h;
    2.0 * vol_j * (-kernel_dw_dr(r, h)) * r / (r * r + eta_sq)
}

/// Iterative pressure Poisson solver for ISPH.
///
/// Solves `∇²P = ρ₀/Δt · ∇·v*` with the Brookshaw SPH Laplacian
/// `Σⱼ cᵢⱼ (Pⱼ − Pᵢ) = bᵢ`, `bᵢ = ρ₀/Δt · (∇·v*)ᵢ` (see `laplacian_weight`),
/// using relaxed Gauss–Seidel (SOR) sweeps
///
/// `Pᵢ ← (1 − ω) Pᵢ + ω (Σⱼ cᵢⱼ Pⱼ − bᵢ) / Σⱼ cᵢⱼ`.
///
/// The Laplacian is singular for pure Neumann data, so particles with a
/// truncated kernel support (the free surface, see
/// [`detect_free_surface`](Self::detect_free_surface)) are pinned to `P = 0`.
/// Compression (`∇·v* < 0`) then produces positive pressure inside the fluid.
/// The neighbour couplings are assembled once per solve, not once per sweep.
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

    /// Flag the free-surface particles of `particles`.
    ///
    /// A particle is a free-surface particle when the SPH divergence of the
    /// position field, `(∇·r)ᵢ = Σⱼ (mⱼ/ρⱼ) (xⱼ − xᵢ)·∇ᵢW(rᵢⱼ, h)` over the
    /// particles within `2h`, is below `2.4`: it equals the space dimension 3 for a
    /// complete kernel support and drops towards 0 when part of the support is
    /// empty.  The indicator keeps a wide margin against the number-density
    /// disorder of a real particle distribution.  The result is index-aligned
    /// with `particles`; an invalid smoothing length flags every particle.
    pub fn detect_free_surface(&self, particles: &[IsphParticle]) -> Vec<bool> {
        let h = self.h;
        if !(h > 0.0 && h.is_finite()) {
            return vec![true; particles.len()];
        }
        let support = 2.0 * h;
        particles
            .iter()
            .enumerate()
            .map(|(i, particle_i)| {
                let mut position_divergence = 0.0_f64;
                for (j, particle_j) in particles.iter().enumerate() {
                    if i == j {
                        continue;
                    }
                    let r = vec3_norm(vec3_sub(particle_i.position, particle_j.position));
                    if r < support {
                        // (xⱼ − xᵢ)·∇ᵢWᵢⱼ = −r dW/dr, positive for a decreasing kernel.
                        position_divergence -= particle_j.volume() * r * kernel_dw_dr(r, h);
                    }
                }
                position_divergence < FREE_SURFACE_POSITION_DIVERGENCE
            })
            .collect()
    }

    /// Run one pressure-Poisson solve; updates `p.pressure` for each particle.
    ///
    /// The right-hand side is read from `p.div_velocity` (see
    /// [`VelocityDivergence`]); free-surface particles
    /// ([`detect_free_surface`](Self::detect_free_surface)) get `p = 0`.
    ///
    /// Returns `(iterations, max_pressure_change)`, where the second value is the
    /// largest absolute pressure change of the last sweep.
    pub fn solve(&self, particles: &mut [IsphParticle], dt: f64) -> (usize, f64) {
        let free_surface = self.detect_free_surface(particles);
        self.solve_with_dirichlet(particles, dt, &free_surface)
    }

    /// [`solve`](Self::solve) with caller-supplied Dirichlet particles.
    ///
    /// `dirichlet[i] == true` pins particle `i` to `p = 0` (missing entries count
    /// as `false`).  Particles without any neighbour inside `2h` are pinned as
    /// well.  When no particle is pinned (a closed domain) the pressure is only
    /// determined up to a constant: the right-hand side is then made compatible
    /// (zero volume-weighted mean) and the pressure normalised to zero
    /// volume-weighted mean.  Non-positive or non-finite `dt` or `h` leave the
    /// particles untouched and return `(0, 0.0)`.
    pub fn solve_with_dirichlet(
        &self,
        particles: &mut [IsphParticle],
        dt: f64,
        dirichlet: &[bool],
    ) -> (usize, f64) {
        let n = particles.len();
        let h = self.h;
        if n == 0 || !(h > 0.0 && h.is_finite()) || !(dt > 0.0 && dt.is_finite()) {
            return (0, 0.0);
        }
        let support = 2.0 * h;
        let is_pinned = |i: usize| dirichlet.get(i).copied().unwrap_or(false);
        let any_pinned = (0..n).any(is_pinned);

        // Assemble the sparse Laplacian once: couplings cᵢⱼ to the unpinned
        // neighbours (pinned neighbours contribute cᵢⱼ·0) and the diagonal Σⱼ cᵢⱼ.
        let volumes: Vec<f64> = particles.iter().map(IsphParticle::volume).collect();
        let mut couplings: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
        let mut diagonal = vec![0.0_f64; n];
        for i in 0..n {
            if is_pinned(i) {
                continue;
            }
            for j in 0..n {
                if i == j {
                    continue;
                }
                let r = vec3_norm(vec3_sub(particles[i].position, particles[j].position));
                if r >= support {
                    continue;
                }
                let c_ij = laplacian_weight(volumes[j], r, h);
                if c_ij <= 0.0 {
                    continue;
                }
                diagonal[i] += c_ij;
                if !is_pinned(j) {
                    couplings[i].push((j, c_ij));
                }
            }
        }
        let active: Vec<bool> = (0..n).map(|i| !is_pinned(i) && diagonal[i] > 0.0).collect();

        let rho0 = self.params.rest_density;
        let mut rhs: Vec<f64> = particles
            .iter()
            .map(|p| rho0 / dt * p.div_velocity)
            .collect();
        let mut pressure: Vec<f64> = particles
            .iter()
            .zip(active.iter())
            .map(|(p, &is_active)| if is_active { p.pressure } else { 0.0 })
            .collect();

        // Closed domain: pure Neumann problem, solvable only for a compatible
        // right-hand side and unique only up to an additive constant.
        let singular = !any_pinned && active.iter().any(|&a| a);
        let volume_total: f64 = (0..n).filter(|&i| active[i]).map(|i| volumes[i]).sum();
        let weighted_mean = |field: &[f64]| -> f64 {
            if volume_total > 0.0 {
                (0..n)
                    .filter(|&i| active[i])
                    .map(|i| volumes[i] * field[i])
                    .sum::<f64>()
                    / volume_total
            } else {
                0.0
            }
        };
        if singular {
            let mean_rhs = weighted_mean(&rhs);
            for (i, b) in rhs.iter_mut().enumerate() {
                if active[i] {
                    *b -= mean_rhs;
                }
            }
        }

        let omega = self.params.relaxation;
        let clamp = self.params.clamp_negative_pressure;
        let mut iters = 0usize;
        let mut max_change = f64::INFINITY;
        let mut relative_change = f64::INFINITY;

        while iters < self.params.max_iter && relative_change > self.params.tolerance {
            max_change = 0.0_f64;
            for i in 0..n {
                if !active[i] {
                    continue;
                }
                let neighbour_sum: f64 = couplings[i].iter().map(|&(j, c)| c * pressure[j]).sum();
                let target = (neighbour_sum - rhs[i]) / diagonal[i];
                let mut updated = (1.0 - omega) * pressure[i] + omega * target;
                if clamp {
                    updated = updated.max(0.0);
                }
                max_change = max_change.max((updated - pressure[i]).abs());
                pressure[i] = updated;
            }
            if singular {
                let mean_pressure = weighted_mean(&pressure);
                for (i, p) in pressure.iter_mut().enumerate() {
                    if active[i] {
                        *p -= mean_pressure;
                    }
                }
            }
            let scale = pressure.iter().fold(0.0_f64, |m, p| m.max(p.abs()));
            relative_change = if scale > 0.0 {
                max_change / scale
            } else {
                max_change
            };
            iters += 1;
        }

        for (particle, &value) in particles.iter_mut().zip(pressure.iter()) {
            particle.pressure = value;
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
    fn predict(&self, particles: &mut [IsphParticle], dt: f64) {
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

            // Viscous acceleration a_visc = ν ∇²v with the Brookshaw/Morris SPH
            // Laplacian ∇²vᵢ ≈ Σⱼ cᵢⱼ (vⱼ − vᵢ), cᵢⱼ = 2 Vⱼ (−dW/dr) r / (r² + η²) ≥ 0,
            // which diffuses velocity differences.  ν is the kinematic viscosity,
            // so no division by the density.
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
                let c_ij = laplacian_weight(mas[j] / rho[j], r, h);
                a_visc = vec3_add(a_visc, vec3_scale(dv, nu * c_ij));
            }

            let v_star = vec3_add(vec3_add(vi, vec3_scale(g, dt)), vec3_scale(a_visc, dt));
            particles[i].velocity_star = v_star;
        }
    }

    /// Correct velocity using pressure gradient.
    fn correct(&self, particles: &mut [IsphParticle], dt: f64) {
        let pg = PressureGradient::new(self.h);
        let accs = pg.all_accelerations(particles);
        for (i, p) in particles.iter_mut().enumerate() {
            p.velocity = vec3_add(p.velocity_star, vec3_scale(accs[i], dt));
            p.position = vec3_add(p.position, vec3_scale(p.velocity, dt));
        }
    }

    /// Advance the particle system by one timestep `dt`.
    ///
    /// Free-surface particles are pinned to `p = 0` by the pressure solve (see
    /// [`IsphPressureSolver::detect_free_surface`]).  The SPH divergence and
    /// pressure-gradient operators have a truncated support there, so the
    /// corrected velocity is divergence-free in the interior but keeps a divergence
    /// error within about two smoothing lengths of the free surface.
    ///
    /// Returns solver statistics from the pressure solve.
    pub fn advance(&self, particles: &mut [IsphParticle], dt: f64) -> (usize, f64) {
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
    pub fn assign_densities(&self, particles: &mut [IsphParticle]) {
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
        for (i, particle_i) in particles.iter().enumerate() {
            let xi = particle_i.position;
            let ci = particle_i.phase as f64;
            let mut wsum = 0.0_f64;
            let mut csum = 0.0_f64;
            for pj in particles.iter() {
                let xj = pj.position;
                let r = vec3_norm(vec3_sub(xi, xj));
                if r >= 2.0 * h {
                    continue;
                }
                let w = kernel_w(r, h);
                csum += w * pj.phase as f64;
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
            for (g_row, dv_a) in g.iter_mut().zip(dv.iter()) {
                for (g_ab, gw_b) in g_row.iter_mut().zip(gw.iter()) {
                    *g_ab += vol_j * dv_a * gw_b;
                }
            }
        }
        // Symmetric part S = ½(G + Gᵀ)
        let mut s2 = 0.0_f64;
        for (alpha, g_row) in g.iter().enumerate() {
            for (beta, g_ab) in g_row.iter().enumerate() {
                let s_ab = 0.5 * (g_ab + g[beta][alpha]);
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

    /// SGS eddy viscosity νₜ of every particle, `result[i] == eddy_viscosity(i, …)`.
    ///
    /// Each νₜ is evaluated exactly once, which is what
    /// [`sgs_accelerations`](Self::sgs_accelerations) needs for the pair-averaged
    /// viscosity.
    pub fn eddy_viscosities(&self, particles: &[IsphParticle]) -> Vec<f64> {
        (0..particles.len())
            .map(|i| self.eddy_viscosity(i, particles))
            .collect()
    }

    /// SGS viscous acceleration for particle i.
    ///
    /// `a_sgs = ∇·(νₜ ∇v) ≈ Σⱼ cᵢⱼ ν̄ᵢⱼ (vⱼ − vᵢ)` with the Brookshaw/Morris weight
    /// `cᵢⱼ = 2 Vⱼ (−dW/dr) r / (r² + η²) ≥ 0` and the pair-averaged eddy
    /// viscosity `ν̄ᵢⱼ = ½ (νₜᵢ + νₜⱼ)`.  νₜ is kinematic, so there is no division
    /// by the density, and the term diffuses velocity differences.
    ///
    /// Only νₜ of particle `i` and of its neighbours within `2h` is evaluated;
    /// use [`sgs_accelerations`](Self::sgs_accelerations) for all particles.
    pub fn sgs_acceleration(&self, i: usize, particles: &[IsphParticle]) -> Vec3 {
        self.sgs_acceleration_with(i, particles, |j| self.eddy_viscosity(j, particles))
    }

    /// SGS viscous acceleration of every particle (see
    /// [`sgs_acceleration`](Self::sgs_acceleration)); the eddy viscosities are
    /// computed once for all particles.
    pub fn sgs_accelerations(&self, particles: &[IsphParticle]) -> Vec<Vec3> {
        let nu_t = self.eddy_viscosities(particles);
        (0..particles.len())
            .map(|i| self.sgs_acceleration_with(i, particles, |j| nu_t[j]))
            .collect()
    }

    /// SGS acceleration of particle `i` given the eddy viscosity `nu_t(j)` of any
    /// particle `j`.
    fn sgs_acceleration_with(
        &self,
        i: usize,
        particles: &[IsphParticle],
        nu_t: impl Fn(usize) -> f64,
    ) -> Vec3 {
        let h = self.h;
        let xi = particles[i].position;
        let vi = particles[i].velocity;
        let nu_i = nu_t(i);
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
            let nu_ij = 0.5 * (nu_i + nu_t(j));
            let c_ij = laplacian_weight(pj.mass / pj.density.max(1e-300), r, h);
            let dv = vec3_sub(pj.velocity, vi);
            acc = vec3_add(acc, vec3_scale(dv, nu_ij * c_ij));
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

    // ── Lattice fixtures for the projection tests ─────────────────────────────

    /// `n`³ lattice of fluid particles (spacing `s`, centred on the origin) with
    /// density `rho0` and mass `rho0 s³`, so every particle owns the volume `s³`.
    fn fluid_lattice(n: usize, s: f64, rho0: f64) -> Vec<IsphParticle> {
        let half = (n as f64 - 1.0) / 2.0;
        let mut particles = Vec::with_capacity(n * n * n);
        for ix in 0..n {
            for iy in 0..n {
                for iz in 0..n {
                    let position = [
                        (ix as f64 - half) * s,
                        (iy as f64 - half) * s,
                        (iz as f64 - half) * s,
                    ];
                    particles.push(IsphParticle::new(position, rho0 * s * s * s, rho0));
                }
            }
        }
        particles
    }

    fn lattice_centre_index(particles: &[IsphParticle]) -> usize {
        particles
            .iter()
            .position(|p| p.position.iter().all(|c| c.abs() < 1e-12))
            .expect("lattice has a centre particle")
    }

    /// Particle whose complete kernel support (and that of its neighbours) lies
    /// inside the lattice, far from the pinned free-surface shell.
    fn is_deep_interior(p: &IsphParticle, s: f64) -> bool {
        p.position.iter().all(|c| c.abs() <= s + 1e-12)
    }

    fn converged_params(rho0: f64) -> IsphPressureSolverParams {
        IsphPressureSolverParams {
            max_iter: 3000,
            tolerance: 1e-12,
            relaxation: 1.0,
            rest_density: rho0,
            clamp_negative_pressure: false,
        }
    }

    /// Residual `Σⱼ cᵢⱼ (pⱼ − pᵢ) − bᵢ` of the pressure Poisson equation at
    /// particle `i` (all particles inside `2h` contribute).
    fn poisson_residual(particles: &[IsphParticle], i: usize, h: f64, b_i: f64) -> f64 {
        let mut sum = 0.0;
        for (j, pj) in particles.iter().enumerate() {
            if i == j {
                continue;
            }
            let r = vec3_norm(vec3_sub(particles[i].position, pj.position));
            if r < 2.0 * h {
                sum += laplacian_weight(pj.volume(), r, h) * (pj.pressure - particles[i].pressure);
            }
        }
        sum - b_i
    }

    /// Lattice with the strained velocity `v* = (sign·x, 0, 0)` (∇·v* = sign, so
    /// −1 compresses and +1 expands the block) and its SPH divergence evaluated.
    fn strained_lattice(s: f64, h: f64, rho0: f64, sign: f64) -> Vec<IsphParticle> {
        let mut particles = fluid_lattice(9, s, rho0);
        for p in particles.iter_mut() {
            p.velocity = [sign * p.position[0], 0.0, 0.0];
        }
        VelocityDivergence::new(h).compute(&mut particles);
        particles
    }

    #[test]
    fn test_detect_free_surface_flags_only_the_outer_shell() {
        let s = 0.1;
        let h = 1.2 * s;
        let particles = fluid_lattice(9, s, 1000.0);
        let solver = IsphPressureSolver::new(h, IsphPressureSolverParams::default());
        let flags = solver.detect_free_surface(&particles);
        assert_eq!(flags.len(), particles.len());
        assert_eq!(
            flags.iter().filter(|&&f| f).count(),
            9 * 9 * 9 - 7 * 7 * 7,
            "exactly the outermost lattice layer is free surface"
        );
        for (p, &flag) in particles.iter().zip(flags.iter()) {
            let on_shell = p.position.iter().any(|c| c.abs() > 3.5 * s);
            assert_eq!(flag, on_shell, "particle at {:?}", p.position);
        }
    }

    #[test]
    fn test_detect_free_surface_flags_everything_for_invalid_h() {
        let particles = fluid_lattice(3, 0.1, 1000.0);
        let solver = IsphPressureSolver::new(0.0, IsphPressureSolverParams::default());
        assert!(solver.detect_free_surface(&particles).iter().all(|&f| f));
    }

    #[test]
    fn test_pressure_solver_compression_matches_cube_poisson_solution() {
        // v* = (−x, 0, 0) compresses the block (∇·v* = −1).  The PPE
        // ∇²p = ρ₀/Δt ∇·v* with p = 0 on the free-surface shell is the Poisson
        // problem of a cube of side L = 8 s, whose centre value is
        // 0.0562 L² ρ₀/Δt (positive: compression raises the pressure).  The former
        // operator solved Σ c (pᵢ + pⱼ) = b, whose solution is negative here and
        // was clamped to zero.
        let (s, rho0, dt) = (0.1, 1000.0, 1e-3);
        let h = 1.2 * s;
        let mut particles = strained_lattice(s, h, rho0, -1.0);
        let solver = IsphPressureSolver::new(h, converged_params(rho0));
        let (iters, change) = solver.solve(&mut particles, dt);
        assert!(iters < 3000, "solver did not converge: {iters} sweeps");
        assert!(change.is_finite());

        let centre = lattice_centre_index(&particles);
        let expected = 0.0562 * (8.0 * s).powi(2) * rho0 / dt;
        let p_centre = particles[centre].pressure;
        assert!(
            p_centre > 0.0,
            "compression must give p > 0, got {p_centre}"
        );
        assert!(
            (p_centre - expected).abs() < 0.12 * expected,
            "p(centre) = {p_centre}, cube solution {expected}"
        );
        let p_max = particles.iter().map(|p| p.pressure).fold(0.0, f64::max);
        assert!(
            (p_max - p_centre).abs() < 1e-9 * p_max,
            "maximum at the centre"
        );
        for p in &particles {
            assert!(
                p.pressure >= -1e-9 * p_max,
                "negative pressure {}",
                p.pressure
            );
        }

        // Free-surface shell pinned to zero.
        let flags = solver.detect_free_surface(&particles);
        for (p, &flag) in particles.iter().zip(flags.iter()) {
            if flag {
                assert_eq!(p.pressure, 0.0);
            }
        }

        // The converged field solves the Poisson equation at every free particle.
        let b_scale = particles
            .iter()
            .map(|p| (rho0 / dt * p.div_velocity).abs())
            .fold(0.0, f64::max);
        for (i, p) in particles.iter().enumerate() {
            if flags[i] {
                continue;
            }
            let residual = poisson_residual(&particles, i, h, rho0 / dt * p.div_velocity);
            assert!(
                residual.abs() < 1e-6 * b_scale,
                "PPE residual {residual} at {:?} (|b| ≤ {b_scale})",
                p.position
            );
        }
    }

    #[test]
    fn test_pressure_solver_default_params_reach_the_pressure_level() {
        // With the stock parameters (100 under-relaxed sweeps) the pressure is
        // already within a percent of the converged cube solution.
        let (s, rho0, dt) = (0.1, 1000.0, 1e-3);
        let h = 1.2 * s;
        let mut particles = strained_lattice(s, h, rho0, -1.0);
        let solver = IsphPressureSolver::new(h, IsphPressureSolverParams::default());
        let (iters, _) = solver.solve(&mut particles, dt);
        assert!(iters <= 100);
        let centre = lattice_centre_index(&particles);
        let expected = 0.0562 * (8.0 * s).powi(2) * rho0 / dt;
        assert!(particles[centre].pressure > 0.85 * expected);
        assert!(particles[centre].pressure < 1.12 * expected);
    }

    #[test]
    fn test_isph_advance_removes_interior_divergence() {
        // After one projection step the deep-interior divergence of the corrected
        // field is below a tenth of the predicted one.  (Particles next to the
        // pinned free-surface shell keep a divergence error: their kernel support
        // is truncated, so they are not part of the measurement.)
        let (s, rho0, dt) = (0.1, 1000.0, 1e-3);
        let h = 1.2 * s;
        let mut particles = strained_lattice(s, h, rho0, -1.0);
        let max_deep_div = |ps: &[IsphParticle]| {
            ps.iter()
                .filter(|p| is_deep_interior(p, s))
                .map(|p| p.div_velocity.abs())
                .fold(0.0, f64::max)
        };
        let before = max_deep_div(&particles);
        assert!(before > 0.9, "predicted divergence {before}");

        let solver = IsphPressureSolver::new(h, converged_params(rho0));
        let step = IsphStep::new(h, solver, 0.0, [0.0; 3]);
        let (iters, _) = step.advance(&mut particles, dt);
        assert!(iters < 3000);
        VelocityDivergence::new(h).compute(&mut particles);
        let after = max_deep_div(&particles);
        assert!(
            after < 0.1 * before,
            "interior divergence {before} -> {after}, expected a > 90 % drop"
        );
        let centre = lattice_centre_index(&particles);
        assert!(particles[centre].pressure > 0.0);
    }

    #[test]
    fn test_pressure_solver_expansion_gives_tension_unless_clamped() {
        // v* = (+x, 0, 0) (∇·v* = +1): the exact projection has negative
        // pressure inside; the clamped solver keeps p = 0 (fluid may separate).
        let (s, rho0, dt) = (0.1, 1000.0, 1e-3);
        let h = 1.2 * s;
        let centre_pressure = |clamp: bool| {
            let mut particles = strained_lattice(s, h, rho0, 1.0);
            let params = IsphPressureSolverParams {
                clamp_negative_pressure: clamp,
                ..converged_params(rho0)
            };
            IsphPressureSolver::new(h, params).solve(&mut particles, dt);
            let min_p = particles.iter().map(|p| p.pressure).fold(0.0, f64::min);
            (particles[lattice_centre_index(&particles)].pressure, min_p)
        };
        let expected = 0.0562 * (8.0 * s).powi(2) * rho0 / dt;
        let (exact_centre, _) = centre_pressure(false);
        assert!(
            (exact_centre + expected).abs() < 0.12 * expected,
            "p(centre) = {exact_centre}, expected {}",
            -expected
        );
        let (clamped_centre, clamped_min) = centre_pressure(true);
        assert_eq!(clamped_centre, 0.0);
        assert!(clamped_min >= 0.0, "clamped pressure {clamped_min} < 0");
    }

    #[test]
    fn test_pressure_solver_closed_domain_fixes_the_pressure_level() {
        // Without any pinned particle the Neumann problem is only defined up to a
        // constant and for a zero-mean right-hand side: the solver removes the
        // mean of b and normalises the pressure to zero mean.
        let (s, rho0, dt) = (0.1, 1000.0, 1e-3);
        let h = 1.2 * s;
        let mut particles = strained_lattice(s, h, rho0, -1.0);
        let n = particles.len();
        let solver = IsphPressureSolver::new(h, converged_params(rho0));
        let (iters, _) = solver.solve_with_dirichlet(&mut particles, dt, &vec![false; n]);
        assert!(iters < 3000, "closed-domain solve did not converge");
        assert!(particles.iter().all(|p| p.pressure.is_finite()));
        let p_max = particles
            .iter()
            .map(|p| p.pressure.abs())
            .fold(0.0, f64::max);
        assert!(p_max > 0.0);
        let mean = particles.iter().map(|p| p.pressure).sum::<f64>() / n as f64;
        assert!(mean.abs() < 1e-9 * p_max, "mean pressure {mean}");

        // Poisson equation with the mean-compensated right-hand side.
        let mean_b = particles
            .iter()
            .map(|p| rho0 / dt * p.div_velocity)
            .sum::<f64>()
            / n as f64;
        let b_scale = particles
            .iter()
            .map(|p| (rho0 / dt * p.div_velocity).abs())
            .fold(0.0, f64::max);
        for (i, p) in particles.iter().enumerate() {
            let residual = poisson_residual(&particles, i, h, rho0 / dt * p.div_velocity - mean_b);
            assert!(
                residual.abs() < 1e-6 * b_scale,
                "closed-domain residual {residual} at particle {i}"
            );
        }
    }

    #[test]
    fn test_pressure_solver_rejects_non_positive_dt() {
        let s = 0.1;
        let mut particles = strained_lattice(s, 1.2 * s, 1000.0, -1.0);
        for p in particles.iter_mut() {
            p.pressure = 5.0;
        }
        let solver = IsphPressureSolver::new(1.2 * s, IsphPressureSolverParams::default());
        for dt in [0.0, -1e-3, f64::NAN, f64::INFINITY] {
            assert_eq!(solver.solve(&mut particles, dt), (0, 0.0), "dt = {dt}");
            assert!(particles.iter().all(|p| p.pressure == 5.0));
        }
    }

    #[test]
    fn test_pressure_solver_isolated_particles_stay_at_zero() {
        let mut particles = two_particles_separated(10.0);
        for p in particles.iter_mut() {
            p.div_velocity = -1.0;
            p.pressure = 3.0;
        }
        let solver = IsphPressureSolver::new(0.1, IsphPressureSolverParams::default());
        solver.solve(&mut particles, 1e-3);
        for p in &particles {
            assert_eq!(p.pressure, 0.0);
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

    // ── Viscosity and SGS regression tests ────────────────────────────────────

    /// Lattice carrying the shear-wave velocity `v = (shear·y + amplitude·cos(k y), 0, 0)`.
    fn shear_wave_lattice(
        s: f64,
        rho0: f64,
        shear: f64,
        amplitude: f64,
        wavenumber: f64,
    ) -> Vec<IsphParticle> {
        let mut particles = fluid_lattice(9, s, rho0);
        for p in particles.iter_mut() {
            let y = p.position[1];
            p.velocity = [shear * y + amplitude * (wavenumber * y).cos(), 0.0, 0.0];
        }
        particles
    }

    fn viscous_step(h: f64, nu: f64) -> IsphStep {
        let solver = IsphPressureSolver::new(h, IsphPressureSolverParams::default());
        IsphStep::new(h, solver, nu, [0.0; 3])
    }

    #[test]
    fn test_predict_viscosity_diffuses_a_velocity_wave() {
        // vₓ = cos(k y), λ = 16 s: at the crest ∇²vₓ = −k², so one viscous step
        // changes vₓ = 1 by −ν k² Δt.  The former anti-diffusive sign gave +ν k² Δt
        // (scaled by 1/ρ).  The Brookshaw Laplacian is accurate to a few percent.
        let (s, nu, dt) = (0.1, 0.01, 1e-3);
        let h = 1.2 * s;
        let k = 2.0 * PI / (16.0 * s);
        let mut particles = shear_wave_lattice(s, 1000.0, 0.0, 1.0, k);
        viscous_step(h, nu).predict(&mut particles, dt);
        let centre = lattice_centre_index(&particles);
        let dv = particles[centre].velocity_star[0] - particles[centre].velocity[0];
        let expected = -nu * k * k * dt;
        assert!(dv < 0.0, "viscosity must damp the crest, Δvₓ = {dv}");
        assert!(
            (dv - expected).abs() < 0.1 * expected.abs(),
            "Δvₓ = {dv}, expected ≈ {expected}"
        );
        assert!(particles[centre].velocity_star[1].abs() < 1e-15);
        assert!(particles[centre].velocity_star[2].abs() < 1e-15);
    }

    #[test]
    fn test_predict_viscosity_is_kinematic_and_independent_of_density() {
        // ν is the kinematic viscosity: the same wave in water and in air (volumes
        // unchanged) relaxes identically.
        let (s, nu, dt) = (0.1, 0.01, 1e-3);
        let h = 1.2 * s;
        let k = 2.0 * PI / (16.0 * s);
        let change = |rho0: f64| {
            let mut particles = shear_wave_lattice(s, rho0, 0.0, 1.0, k);
            viscous_step(h, nu).predict(&mut particles, dt);
            let centre = lattice_centre_index(&particles);
            particles[centre].velocity_star[0] - particles[centre].velocity[0]
        };
        let water = change(1000.0);
        let air = change(1.2);
        assert!(
            (water - air).abs() < 1e-12 * water.abs(),
            "Δvₓ(water) = {water}, Δvₓ(air) = {air}"
        );
    }

    #[test]
    fn test_predict_viscosity_leaves_uniform_flow_unchanged() {
        let s = 0.1;
        let mut particles = shear_wave_lattice(s, 1000.0, 0.0, 0.0, 1.0);
        for p in particles.iter_mut() {
            p.velocity = [0.3, -0.2, 0.1];
        }
        viscous_step(1.2 * s, 0.05).predict(&mut particles, 1e-3);
        for p in &particles {
            for (star, v) in p.velocity_star.iter().zip(p.velocity.iter()) {
                assert!((star - v).abs() < 1e-14);
            }
        }
    }

    #[test]
    fn test_sgs_acceleration_is_diffusive_and_kinematic() {
        // vₓ = S₀ y + a cos(k y): |S̃| = |∂ᵧvₓ|, so the SGS flux is (CₛΔ)² (∂ᵧvₓ)²
        // and at y = 0 the acceleration is ∂ᵧ((CₛΔ)² (∂ᵧvₓ)²) = −2 (CₛΔ)² S₀ a k²
        // (negative: it damps the cosine perturbation).  The former implementation
        // had the opposite sign and divided by the density.
        let (s, rho0) = (0.1, 1000.0);
        let h = 1.2 * s;
        let k = 2.0 * PI / (16.0 * s);
        let (shear, amplitude) = (1.0, 0.05);
        let (cs, delta) = (0.15, h);
        let turbulence = IsphTurbulence::new(cs, delta, h);
        let expected = -2.0 * (cs * delta).powi(2) * shear * amplitude * k * k;

        let acceleration = |rho: f64| {
            let particles = shear_wave_lattice(s, rho, shear, amplitude, k);
            let centre = lattice_centre_index(&particles);
            turbulence.sgs_acceleration(centre, &particles)
        };
        let water = acceleration(rho0);
        assert!(
            water[0] < 0.0,
            "SGS term must be diffusive, a_x = {}",
            water[0]
        );
        assert!(
            (water[0] - expected).abs() < 0.2 * expected.abs(),
            "a_x = {}, expected ≈ {expected}",
            water[0]
        );
        assert!(water[1].abs() < 1e-12 && water[2].abs() < 1e-12);
        let air = acceleration(1.2);
        assert!(
            (water[0] - air[0]).abs() < 1e-9 * water[0].abs(),
            "kinematic eddy viscosity must not depend on ρ: {} vs {}",
            water[0],
            air[0]
        );
    }

    #[test]
    fn test_sgs_batch_methods_match_single_particle_calls() {
        let (s, rho0) = (0.1, 1000.0);
        let h = 1.2 * s;
        let k = 2.0 * PI / (16.0 * s);
        let particles = shear_wave_lattice(s, rho0, 1.0, 0.05, k);
        let turbulence = IsphTurbulence::new(0.15, h, h);
        let nu_t = turbulence.eddy_viscosities(&particles);
        let batch = turbulence.sgs_accelerations(&particles);
        assert_eq!(nu_t.len(), particles.len());
        assert_eq!(batch.len(), particles.len());
        let centre = lattice_centre_index(&particles);
        assert!((nu_t[centre] - turbulence.eddy_viscosity(centre, &particles)).abs() < 1e-18);
        let single = turbulence.sgs_acceleration(centre, &particles);
        for (b, s1) in batch[centre].iter().zip(single.iter()) {
            assert!((b - s1).abs() < 1e-15, "{b} vs {s1}");
        }
        assert!(nu_t.iter().all(|&v| v >= 0.0 && v.is_finite()));
    }

    #[test]
    fn test_sgs_acceleration_vanishes_for_uniform_flow() {
        let s = 0.1;
        let mut particles = shear_wave_lattice(s, 1000.0, 0.0, 0.0, 1.0);
        for p in particles.iter_mut() {
            p.velocity = [0.5, 0.5, -0.25];
        }
        let turbulence = IsphTurbulence::new(0.15, 1.2 * s, 1.2 * s);
        for acc in turbulence.sgs_accelerations(&particles) {
            for a in acc {
                assert!(a.abs() < 1e-14, "uniform flow has no SGS force, got {a}");
            }
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

#[cfg(test)]
mod disorder_tests;
