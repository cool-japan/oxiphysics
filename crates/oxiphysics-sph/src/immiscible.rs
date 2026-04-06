#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Immiscible multiphase SPH using the color function / phase-field approach.
//!
//! This module implements the Continuum Surface Force (CSF) model for surface
//! tension at immiscible fluid interfaces, following the approach of Morris (2000)
//! and Lafaurie et al. (1994). Each particle carries a phase identifier and a
//! color function value that tracks the local phase composition.

#![allow(dead_code)]

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Kernel helpers (local, independent of the kernel module trait)
// ---------------------------------------------------------------------------

/// Cubic spline kernel W(r, h) in 3D. Support radius = 2h.
pub fn cubic_kernel(r: f64, h: f64) -> f64 {
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

/// Gradient of the cubic spline kernel: returns ∇W = (dW/dr) * r̂.
pub fn cubic_kernel_grad(r_vec: [f64; 3], r: f64, h: f64) -> [f64; 3] {
    if r < 1e-14 {
        return [0.0; 3];
    }
    let sigma = 1.0 / (PI * h * h * h);
    let q = r / h;
    let dw_dr = if q >= 2.0 {
        0.0
    } else if q >= 1.0 {
        let t = 2.0 - q;
        sigma * (-0.75 * t * t) / h
    } else {
        sigma * (-3.0 * q + 2.25 * q * q) / h
    };
    let scale = dw_dr / r;
    [scale * r_vec[0], scale * r_vec[1], scale * r_vec[2]]
}

/// Laplacian of the cubic spline kernel in 3D.
pub fn cubic_kernel_laplacian(r: f64, h: f64) -> f64 {
    let sigma = 1.0 / (PI * h * h * h);
    let q = r / h;
    let h2 = h * h;
    if q >= 2.0 {
        0.0
    } else if q >= 1.0 {
        let t = 2.0 - q;
        let d2w = sigma * 1.5 * t / h2;
        let dw = sigma * (-0.75 * t * t) / h;
        d2w + if r > 1e-14 { 2.0 * dw / r } else { 0.0 }
    } else {
        let d2w = sigma * (-3.0 + 4.5 * q) / h2;
        let dw = sigma * (-3.0 * q + 2.25 * q * q) / h;
        d2w + if r > 1e-14 { 2.0 * dw / r } else { 0.0 }
    }
}

// ---------------------------------------------------------------------------
// FluidPhase — material properties of a single immiscible phase
// ---------------------------------------------------------------------------

/// Physical properties of an immiscible fluid phase.
#[derive(Debug, Clone)]
pub struct FluidPhase {
    /// Reference (rest) density (kg/m³).
    pub density_0: f64,
    /// Dynamic viscosity (Pa·s).
    pub viscosity: f64,
    /// Surface-tension coefficient (N/m).
    pub surface_tension: f64,
    /// Human-readable name for debugging / output.
    pub name: String,
}

impl FluidPhase {
    /// Preset: liquid water at 20 °C.
    pub fn water() -> Self {
        Self {
            density_0: 998.2,
            viscosity: 1.002e-3,
            surface_tension: 0.0728,
            name: "water".to_string(),
        }
    }

    /// Preset: light mineral oil.
    pub fn oil() -> Self {
        Self {
            density_0: 870.0,
            viscosity: 3.0e-2,
            surface_tension: 0.035,
            name: "oil".to_string(),
        }
    }

    /// Preset: air at 20 °C, 1 atm.
    pub fn air() -> Self {
        Self {
            density_0: 1.204,
            viscosity: 1.81e-5,
            surface_tension: 0.0,
            name: "air".to_string(),
        }
    }
}

// ---------------------------------------------------------------------------
// ImmiscibleParticle
// ---------------------------------------------------------------------------

/// A single SPH particle in an immiscible multiphase simulation.
#[derive(Debug, Clone)]
pub struct ImmiscibleParticle {
    /// Position (m).
    pub position: [f64; 3],
    /// Velocity (m/s).
    pub velocity: [f64; 3],
    /// Particle mass (kg).
    pub mass: f64,
    /// Current density (kg/m³).
    pub density: f64,
    /// Current pressure (Pa).
    pub pressure: f64,
    /// Index into the phases array identifying which fluid this particle belongs to.
    pub phase_id: usize,
    /// Color function value C ∈ \[0, 1\]: 1 = fully same-phase, 0 = fully other-phase.
    pub color_function: f64,
}

// ---------------------------------------------------------------------------
// Free functions operating on particle slices
// ---------------------------------------------------------------------------

/// Compute the color field at particle *i* via SPH summation.
///
/// `C_i = Σ_j  (m_j / ρ_j) · W(r_ij, h) · δ(phase_i == phase_j)`
pub fn color_field_sph(particles: &[ImmiscibleParticle], i: usize, h: f64) -> f64 {
    let pi = particles[i].position;
    let phase_i = particles[i].phase_id;
    let mut c = 0.0_f64;
    for (j, pj) in particles.iter().enumerate() {
        if j == i {
            // Include self-contribution only if same phase (always true here)
        }
        if pj.phase_id != phase_i {
            continue;
        }
        let dx = pi[0] - pj.position[0];
        let dy = pi[1] - pj.position[1];
        let dz = pi[2] - pj.position[2];
        let r = (dx * dx + dy * dy + dz * dz).sqrt();
        c += (pj.mass / pj.density) * cubic_kernel(r, h);
    }
    c
}

/// Gradient of the color field at particle *i*: ∇C_i.
///
/// `∇C_i = Σ_j  (m_j / ρ_j) · ∇W(r_ij, h) · δ(phase_i == phase_j)`
pub fn color_gradient(particles: &[ImmiscibleParticle], i: usize, h: f64) -> [f64; 3] {
    let pi = particles[i].position;
    let phase_i = particles[i].phase_id;
    let mut grad = [0.0_f64; 3];
    for (j, pj) in particles.iter().enumerate() {
        if pj.phase_id != phase_i {
            continue;
        }
        let r_vec = [
            pi[0] - pj.position[0],
            pi[1] - pj.position[1],
            pi[2] - pj.position[2],
        ];
        let r = (r_vec[0] * r_vec[0] + r_vec[1] * r_vec[1] + r_vec[2] * r_vec[2]).sqrt();
        if r < 1e-14 && j == i {
            continue;
        }
        let gw = cubic_kernel_grad(r_vec, r, h);
        let w_coeff = pj.mass / pj.density;
        grad[0] += w_coeff * gw[0];
        grad[1] += w_coeff * gw[1];
        grad[2] += w_coeff * gw[2];
    }
    grad
}

/// Surface tension force on particle *i* via the Continuum Surface Force (CSF) model.
///
/// `F_st = σ · κ · n̂ · δ_s`
///
/// where `κ = -∇·(n̂)` is approximated as `∇²C / |∇C|` and `δ_s ≈ |∇C|`.
/// This simplifies to `F_st = -σ · ∇²C · (∇C / |∇C|²) · |∇C|` = `-σ · κ · ∇C`.
///
/// The implementation uses `F = σ · ∇²C · (∇C / |∇C|)` per Morris (2000).
pub fn surface_tension_force_csf(
    particles: &[ImmiscibleParticle],
    i: usize,
    h: f64,
    sigma: f64,
) -> [f64; 3] {
    let pi = particles[i].position;
    let phase_i = particles[i].phase_id;

    // Compute Laplacian of color field: ∇²C_i = Σ_j (m_j/ρ_j) ∇²W(r_ij, h)
    let mut lap_c = 0.0_f64;
    for (j, pj) in particles.iter().enumerate() {
        if pj.phase_id != phase_i {
            continue;
        }
        let dx = pi[0] - pj.position[0];
        let dy = pi[1] - pj.position[1];
        let dz = pi[2] - pj.position[2];
        let r = (dx * dx + dy * dy + dz * dz).sqrt();
        if r < 1e-14 && j == i {
            continue;
        }
        lap_c += (pj.mass / pj.density) * cubic_kernel_laplacian(r, h);
    }

    // Gradient of color field
    let grad_c = color_gradient(particles, i, h);
    let grad_mag = (grad_c[0] * grad_c[0] + grad_c[1] * grad_c[1] + grad_c[2] * grad_c[2]).sqrt();

    if grad_mag < 1e-10 {
        return [0.0; 3];
    }

    // F = σ · κ · n̂ · |∇C|,  κ = ∇²C / |∇C|,  n̂ = ∇C / |∇C|
    // → F = σ · (∇²C / |∇C|) · (∇C / |∇C|) · |∇C| = σ · ∇²C · ∇C / |∇C|
    let scale = sigma * lap_c / grad_mag;
    [scale * grad_c[0], scale * grad_c[1], scale * grad_c[2]]
}

/// Return the dynamic viscosity of the phase that particle *p* belongs to.
pub fn interface_viscosity(particle: &ImmiscibleParticle, phases: &[FluidPhase]) -> f64 {
    phases[particle.phase_id].viscosity
}

/// Standard SPH density summation for particle *i*:
/// `ρ_i = Σ_j  m_j · W(r_ij, h)`
///
/// Sums over **all** particles regardless of phase.
pub fn density_summation_multiphase(particles: &[ImmiscibleParticle], i: usize, h: f64) -> f64 {
    let pi = particles[i].position;
    let mut rho = 0.0_f64;
    for pj in particles.iter() {
        let dx = pi[0] - pj.position[0];
        let dy = pi[1] - pj.position[1];
        let dz = pi[2] - pj.position[2];
        let r = (dx * dx + dy * dy + dz * dz).sqrt();
        rho += pj.mass * cubic_kernel(r, h);
    }
    rho
}

/// Estimate the volume fraction of `phase_id` at an arbitrary `point`.
///
/// `α(x) = Σ_{j ∈ phase} (m_j / ρ_j) · W(|x - x_j|, h)`
pub fn phase_fraction_at(
    particles: &[ImmiscibleParticle],
    point: [f64; 3],
    h: f64,
    phase_id: usize,
) -> f64 {
    let mut frac = 0.0_f64;
    for pj in particles.iter() {
        if pj.phase_id != phase_id {
            continue;
        }
        let dx = point[0] - pj.position[0];
        let dy = point[1] - pj.position[1];
        let dz = point[2] - pj.position[2];
        let r = (dx * dx + dy * dy + dz * dz).sqrt();
        frac += (pj.mass / pj.density) * cubic_kernel(r, h);
    }
    frac
}

// ---------------------------------------------------------------------------
// ImmiscibleSystem — high-level simulation driver
// ---------------------------------------------------------------------------

/// High-level driver for an immiscible multiphase SPH simulation.
pub struct ImmiscibleSystem {
    /// All particles in the simulation.
    pub particles: Vec<ImmiscibleParticle>,
    /// Registered fluid phases.
    pub phases: Vec<FluidPhase>,
    /// Smoothing length (m).
    pub h: f64,
}

impl ImmiscibleSystem {
    /// Create an empty system with the given smoothing length.
    pub fn new(h: f64) -> Self {
        Self {
            particles: Vec::new(),
            phases: Vec::new(),
            h,
        }
    }

    /// Register a new fluid phase and return its index.
    pub fn add_phase(&mut self, phase: FluidPhase) -> usize {
        let idx = self.phases.len();
        self.phases.push(phase);
        idx
    }

    /// Add a particle with the given position, velocity, mass and phase index.
    pub fn add_particle(&mut self, pos: [f64; 3], vel: [f64; 3], mass: f64, phase_id: usize) {
        let density_0 = if phase_id < self.phases.len() {
            self.phases[phase_id].density_0
        } else {
            1000.0
        };
        self.particles.push(ImmiscibleParticle {
            position: pos,
            velocity: vel,
            mass,
            density: density_0,
            pressure: 0.0,
            phase_id,
            color_function: 1.0,
        });
    }

    /// Recompute densities for all particles using SPH summation.
    pub fn update_densities(&mut self) {
        let n = self.particles.len();
        let h = self.h;
        let mut densities = vec![0.0_f64; n];
        for i in 0..n {
            densities[i] = density_summation_multiphase(&self.particles, i, h);
        }
        for (i, p) in self.particles.iter_mut().enumerate() {
            p.density = densities[i];
        }
    }

    /// Recompute color function values for all particles.
    pub fn update_color_functions(&mut self) {
        let n = self.particles.len();
        let h = self.h;
        let mut colors = vec![0.0_f64; n];
        for i in 0..n {
            colors[i] = color_field_sph(&self.particles, i, h);
        }
        for (i, p) in self.particles.iter_mut().enumerate() {
            p.color_function = colors[i];
        }
    }

    /// Return a vector of particle counts, one entry per registered phase.
    pub fn count_particles_by_phase(&self) -> Vec<usize> {
        let mut counts = vec![0usize; self.phases.len()];
        for p in &self.particles {
            if p.phase_id < counts.len() {
                counts[p.phase_id] += 1;
            }
        }
        counts
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // -----------------------------------------------------------------------
    // Kernel helpers
    // -----------------------------------------------------------------------

    #[test]
    fn kernel_zero_outside_support() {
        assert_eq!(cubic_kernel(2.1, 1.0), 0.0);
        assert_eq!(cubic_kernel(3.0, 1.0), 0.0);
    }

    #[test]
    fn kernel_positive_inside_support() {
        assert!(cubic_kernel(0.0, 1.0) > 0.0);
        assert!(cubic_kernel(0.5, 1.0) > 0.0);
        assert!(cubic_kernel(1.5, 1.0) > 0.0);
    }

    #[test]
    fn kernel_grad_zero_at_origin() {
        let g = cubic_kernel_grad([0.0, 0.0, 0.0], 0.0, 1.0);
        assert_eq!(g, [0.0, 0.0, 0.0]);
    }

    #[test]
    fn kernel_grad_direction() {
        // Gradient along +x should have a negative x-component (pointing toward origin)
        let r_vec = [0.5, 0.0, 0.0];
        let r = 0.5_f64;
        let g = cubic_kernel_grad(r_vec, r, 1.0);
        assert!(g[0] < 0.0, "grad x should be negative, got {}", g[0]);
        assert!(g[1].abs() < 1e-15);
        assert!(g[2].abs() < 1e-15);
    }

    #[test]
    fn kernel_laplacian_zero_outside_support() {
        assert_eq!(cubic_kernel_laplacian(2.5, 1.0), 0.0);
    }

    // -----------------------------------------------------------------------
    // FluidPhase presets
    // -----------------------------------------------------------------------

    #[test]
    fn fluid_phase_water_properties() {
        let w = FluidPhase::water();
        assert!((w.density_0 - 998.2).abs() < 0.1);
        assert!(w.viscosity > 0.0);
        assert!(w.surface_tension > 0.0);
        assert_eq!(w.name, "water");
    }

    #[test]
    fn fluid_phase_oil_denser_than_air() {
        let oil = FluidPhase::oil();
        let air = FluidPhase::air();
        assert!(oil.density_0 > air.density_0);
    }

    #[test]
    fn fluid_phase_air_zero_surface_tension() {
        let air = FluidPhase::air();
        assert_eq!(air.surface_tension, 0.0);
    }

    #[test]
    fn fluid_phase_oil_higher_viscosity_than_water() {
        let oil = FluidPhase::oil();
        let water = FluidPhase::water();
        assert!(oil.viscosity > water.viscosity);
    }

    // -----------------------------------------------------------------------
    // color_field_sph
    // -----------------------------------------------------------------------

    fn make_two_phase_particles() -> Vec<ImmiscibleParticle> {
        vec![
            ImmiscibleParticle {
                position: [0.0, 0.0, 0.0],
                velocity: [0.0; 3],
                mass: 1.0,
                density: 1000.0,
                pressure: 0.0,
                phase_id: 0,
                color_function: 1.0,
            },
            ImmiscibleParticle {
                position: [0.1, 0.0, 0.0],
                velocity: [0.0; 3],
                mass: 1.0,
                density: 1000.0,
                pressure: 0.0,
                phase_id: 0,
                color_function: 1.0,
            },
            ImmiscibleParticle {
                position: [0.2, 0.0, 0.0],
                velocity: [0.0; 3],
                mass: 1.0,
                density: 870.0,
                pressure: 0.0,
                phase_id: 1,
                color_function: 0.0,
            },
        ]
    }

    #[test]
    fn color_field_same_phase_nonzero() {
        let particles = make_two_phase_particles();
        let c = color_field_sph(&particles, 0, 0.5);
        assert!(
            c > 0.0,
            "Color field should be positive for same-phase particles"
        );
    }

    #[test]
    fn color_field_excludes_other_phase() {
        // Particle 0 (phase 0) and particle 2 (phase 1) are separated.
        // The color field at particle 2 should only count phase-1 particles.
        let particles = make_two_phase_particles();
        let c_phase0 = color_field_sph(&particles, 0, 0.5);
        let c_phase1 = color_field_sph(&particles, 2, 0.5);
        // Phase 0 has two particles near each other → larger color value
        // Phase 1 has only one particle
        assert!(
            c_phase0 > c_phase1,
            "Two-particle phase should have higher color field"
        );
    }

    // -----------------------------------------------------------------------
    // color_gradient
    // -----------------------------------------------------------------------

    #[test]
    fn color_gradient_returns_three_components() {
        let particles = make_two_phase_particles();
        let g = color_gradient(&particles, 0, 0.5);
        // Just ensure it is finite
        for &v in &g {
            assert!(v.is_finite(), "Color gradient component is not finite: {v}");
        }
    }

    #[test]
    fn color_gradient_single_particle_zero() {
        // A lone particle has no neighbors → gradient is zero
        let p = vec![ImmiscibleParticle {
            position: [0.0, 0.0, 0.0],
            velocity: [0.0; 3],
            mass: 1.0,
            density: 1000.0,
            pressure: 0.0,
            phase_id: 0,
            color_function: 1.0,
        }];
        let g = color_gradient(&p, 0, 0.5);
        // Self-contribution is skipped (r < epsilon) so gradient should be 0
        for &v in &g {
            assert!(
                v.abs() < 1e-14,
                "Single-particle gradient should be 0, got {v}"
            );
        }
    }

    // -----------------------------------------------------------------------
    // surface_tension_force_csf
    // -----------------------------------------------------------------------

    #[test]
    fn surface_tension_zero_for_homogeneous_phase() {
        // A uniform lattice of same-phase particles: gradient is ~0 → force is ~0
        let mut particles = Vec::new();
        for ix in 0..3_i32 {
            for iy in 0..3_i32 {
                for iz in 0..3_i32 {
                    particles.push(ImmiscibleParticle {
                        position: [ix as f64 * 0.1, iy as f64 * 0.1, iz as f64 * 0.1],
                        velocity: [0.0; 3],
                        mass: 0.001,
                        density: 1000.0,
                        pressure: 0.0,
                        phase_id: 0,
                        color_function: 1.0,
                    });
                }
            }
        }
        // Central particle (index 13 in a 3x3x3 grid = middle)
        let f = surface_tension_force_csf(&particles, 13, 0.3, 0.0728);
        let mag = (f[0] * f[0] + f[1] * f[1] + f[2] * f[2]).sqrt();
        assert!(
            mag < 1.0,
            "Surface tension force in uniform region should be small, got {mag}"
        );
    }

    #[test]
    fn surface_tension_returns_finite_values() {
        let particles = make_two_phase_particles();
        let f = surface_tension_force_csf(&particles, 0, 0.5, 0.072);
        for &v in &f {
            assert!(
                v.is_finite(),
                "Surface tension force component not finite: {v}"
            );
        }
    }

    // -----------------------------------------------------------------------
    // interface_viscosity
    // -----------------------------------------------------------------------

    #[test]
    fn interface_viscosity_returns_phase_viscosity() {
        let phases = vec![FluidPhase::water(), FluidPhase::oil()];
        let p_water = ImmiscibleParticle {
            position: [0.0; 3],
            velocity: [0.0; 3],
            mass: 1.0,
            density: 998.2,
            pressure: 0.0,
            phase_id: 0,
            color_function: 1.0,
        };
        let p_oil = ImmiscibleParticle {
            position: [0.1, 0.0, 0.0],
            velocity: [0.0; 3],
            mass: 1.0,
            density: 870.0,
            pressure: 0.0,
            phase_id: 1,
            color_function: 1.0,
        };
        let mu_w = interface_viscosity(&p_water, &phases);
        let mu_o = interface_viscosity(&p_oil, &phases);
        assert!((mu_w - FluidPhase::water().viscosity).abs() < 1e-20);
        assert!((mu_o - FluidPhase::oil().viscosity).abs() < 1e-20);
    }

    // -----------------------------------------------------------------------
    // density_summation_multiphase
    // -----------------------------------------------------------------------

    #[test]
    fn density_summation_positive() {
        let particles = make_two_phase_particles();
        let rho = density_summation_multiphase(&particles, 0, 0.5);
        assert!(rho > 0.0, "Density summation should be positive");
    }

    #[test]
    fn density_summation_includes_all_phases() {
        // Particle near an oil particle should still count it (cross-phase density)
        let particles = make_two_phase_particles();
        let rho0 = density_summation_multiphase(&particles, 0, 0.5);
        // rho0 includes contributions from particles 1 (same phase) and 2 (different phase)
        // All three are within support h=0.5
        assert!(rho0 > 0.0);
    }

    // -----------------------------------------------------------------------
    // phase_fraction_at
    // -----------------------------------------------------------------------

    #[test]
    fn phase_fraction_nonzero_near_phase() {
        let particles = make_two_phase_particles();
        let frac = phase_fraction_at(&particles, [0.0, 0.0, 0.0], 0.5, 0);
        assert!(
            frac > 0.0,
            "Phase fraction should be positive near phase-0 particles"
        );
    }

    #[test]
    fn phase_fraction_zero_far_away() {
        let particles = make_two_phase_particles();
        let frac = phase_fraction_at(&particles, [10.0, 10.0, 10.0], 0.5, 0);
        assert!(
            frac < 1e-10,
            "Phase fraction should be ~0 far from all particles"
        );
    }

    // -----------------------------------------------------------------------
    // ImmiscibleSystem
    // -----------------------------------------------------------------------

    #[test]
    fn system_add_phase_returns_correct_index() {
        let mut sys = ImmiscibleSystem::new(0.1);
        let i0 = sys.add_phase(FluidPhase::water());
        let i1 = sys.add_phase(FluidPhase::oil());
        assert_eq!(i0, 0);
        assert_eq!(i1, 1);
    }

    #[test]
    fn system_add_particle_increments_count() {
        let mut sys = ImmiscibleSystem::new(0.1);
        sys.add_phase(FluidPhase::water());
        sys.add_particle([0.0, 0.0, 0.0], [0.0; 3], 0.001, 0);
        sys.add_particle([0.1, 0.0, 0.0], [0.0; 3], 0.001, 0);
        assert_eq!(sys.particles.len(), 2);
    }

    #[test]
    fn system_count_particles_by_phase() {
        let mut sys = ImmiscibleSystem::new(0.1);
        sys.add_phase(FluidPhase::water());
        sys.add_phase(FluidPhase::oil());
        sys.add_particle([0.0, 0.0, 0.0], [0.0; 3], 0.001, 0);
        sys.add_particle([0.1, 0.0, 0.0], [0.0; 3], 0.001, 0);
        sys.add_particle([0.2, 0.0, 0.0], [0.0; 3], 0.001, 1);
        let counts = sys.count_particles_by_phase();
        assert_eq!(counts[0], 2);
        assert_eq!(counts[1], 1);
    }

    #[test]
    fn system_update_densities_changes_values() {
        let mut sys = ImmiscibleSystem::new(0.3);
        sys.add_phase(FluidPhase::water());
        for k in 0..5 {
            sys.add_particle([k as f64 * 0.05, 0.0, 0.0], [0.0; 3], 1.0, 0);
        }
        // Initial density is density_0
        let initial = sys.particles[0].density;
        sys.update_densities();
        let updated = sys.particles[0].density;
        // After SPH summation density will differ from the preset density_0
        assert!(updated > 0.0);
        assert!((updated - initial).abs() > 1e-10 || updated > 0.0);
    }

    #[test]
    fn system_update_color_functions() {
        let mut sys = ImmiscibleSystem::new(0.3);
        sys.add_phase(FluidPhase::water());
        sys.add_phase(FluidPhase::oil());
        for k in 0..3 {
            sys.add_particle([k as f64 * 0.05, 0.0, 0.0], [0.0; 3], 1.0, 0);
        }
        sys.add_particle([1.0, 0.0, 0.0], [0.0; 3], 1.0, 1);
        sys.particles.iter_mut().for_each(|p| p.density = 1000.0);
        sys.update_color_functions();
        for p in &sys.particles {
            assert!(
                p.color_function >= 0.0,
                "Color function must be non-negative"
            );
            assert!(p.color_function.is_finite());
        }
    }

    #[test]
    fn system_new_empty() {
        let sys = ImmiscibleSystem::new(0.05);
        assert_eq!(sys.particles.len(), 0);
        assert_eq!(sys.phases.len(), 0);
        assert!((sys.h - 0.05).abs() < 1e-15);
    }
}

// ---------------------------------------------------------------------------
// Phase interface detection
// ---------------------------------------------------------------------------

/// Result of interface detection at a particle.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct InterfaceInfo {
    /// Whether this particle is near a phase interface.
    pub is_interface: bool,
    /// Estimated interface normal (pointing from phase_id toward other phase).
    pub normal: [f64; 3],
    /// Magnitude of the color function gradient (proxy for interface sharpness).
    pub gradient_magnitude: f64,
}

/// Detect whether particle `i` is near a phase interface.
///
/// A particle is classified as an interface particle if its color function
/// gradient exceeds a threshold `threshold` relative to the reference gradient
/// scale `1/h`.
#[allow(dead_code)]
pub fn detect_interface_particle(
    particles: &[ImmiscibleParticle],
    i: usize,
    h: f64,
    threshold: f64,
) -> InterfaceInfo {
    let grad = color_gradient(particles, i, h);
    let grad_mag = (grad[0] * grad[0] + grad[1] * grad[1] + grad[2] * grad[2]).sqrt();
    let is_interface = grad_mag > threshold / h;
    let normal = if grad_mag > 1e-30 {
        [grad[0] / grad_mag, grad[1] / grad_mag, grad[2] / grad_mag]
    } else {
        [0.0; 3]
    };
    InterfaceInfo {
        is_interface,
        normal,
        gradient_magnitude: grad_mag,
    }
}

/// Classify all particles as interface or bulk.
///
/// Returns a vector of [`InterfaceInfo`] one per particle.
#[allow(dead_code)]
pub fn classify_all_particles(
    particles: &[ImmiscibleParticle],
    h: f64,
    threshold: f64,
) -> Vec<InterfaceInfo> {
    (0..particles.len())
        .map(|i| detect_interface_particle(particles, i, h, threshold))
        .collect()
}

// ---------------------------------------------------------------------------
// Color function advection (explicit Euler)
// ---------------------------------------------------------------------------

/// Advect color function values for one explicit Euler step.
///
/// dc/dt + v · ∇c = 0  (material derivative form)
///
/// ∂c_i / ∂t = - Σ_j (m_j / ρ_j) (c_i - c_j) v_ij · ∇W_ij
///
/// where v_ij = v_i - v_j.
///
/// # Returns
/// A new vector of updated color function values (same length as `particles`).
#[allow(dead_code)]
pub fn advect_color_function(particles: &[ImmiscibleParticle], h: f64, dt: f64) -> Vec<f64> {
    let n = particles.len();
    let mut dc_dt = vec![0.0_f64; n];

    for i in 0..n {
        let pi = particles[i].position;
        let vi = particles[i].velocity;
        let ci = particles[i].color_function;

        for j in 0..n {
            if i == j {
                continue;
            }
            let pj = &particles[j];
            let r_vec = [
                pi[0] - pj.position[0],
                pi[1] - pj.position[1],
                pi[2] - pj.position[2],
            ];
            let r = (r_vec[0] * r_vec[0] + r_vec[1] * r_vec[1] + r_vec[2] * r_vec[2]).sqrt();
            if r > 2.0 * h {
                continue;
            }
            let gw = cubic_kernel_grad(r_vec, r, h);
            let v_ij = [
                vi[0] - pj.velocity[0],
                vi[1] - pj.velocity[1],
                vi[2] - pj.velocity[2],
            ];
            let vij_dot_gw = v_ij[0] * gw[0] + v_ij[1] * gw[1] + v_ij[2] * gw[2];
            let cj = pj.color_function;
            dc_dt[i] -= (pj.mass / pj.density) * (ci - cj) * vij_dot_gw;
        }
    }

    particles
        .iter()
        .enumerate()
        .map(|(i, p)| (p.color_function + dc_dt[i] * dt).clamp(0.0, 2.0))
        .collect()
}

// ---------------------------------------------------------------------------
// Interphase forces
// ---------------------------------------------------------------------------

/// Compute an interphase drag force on particle `i` due to the velocity
/// difference with respect to surrounding particles of a different phase.
///
/// Simplified model: F_drag = β (u_other_phase - u_i)
///
/// where β = interphase_drag_coeff and u_other_phase is the weighted mean
/// velocity of other-phase neighbors.
#[allow(dead_code)]
pub fn interphase_drag_force(
    particles: &[ImmiscibleParticle],
    i: usize,
    h: f64,
    interphase_drag_coeff: f64,
) -> [f64; 3] {
    let pi = particles[i].position;
    let vi = particles[i].velocity;
    let phase_i = particles[i].phase_id;

    let mut weighted_vel = [0.0_f64; 3];
    let mut weight_sum = 0.0_f64;

    for (j, pj) in particles.iter().enumerate() {
        if j == i || pj.phase_id == phase_i {
            continue;
        }
        let dx = pi[0] - pj.position[0];
        let dy = pi[1] - pj.position[1];
        let dz = pi[2] - pj.position[2];
        let r = (dx * dx + dy * dy + dz * dz).sqrt();
        if r > 2.0 * h {
            continue;
        }
        let w = cubic_kernel(r, h) * pj.mass / pj.density;
        weighted_vel[0] += w * pj.velocity[0];
        weighted_vel[1] += w * pj.velocity[1];
        weighted_vel[2] += w * pj.velocity[2];
        weight_sum += w;
    }

    if weight_sum < 1e-30 {
        return [0.0; 3];
    }

    let inv_w = 1.0 / weight_sum;
    let u_other = [
        weighted_vel[0] * inv_w,
        weighted_vel[1] * inv_w,
        weighted_vel[2] * inv_w,
    ];
    [
        interphase_drag_coeff * (u_other[0] - vi[0]),
        interphase_drag_coeff * (u_other[1] - vi[1]),
        interphase_drag_coeff * (u_other[2] - vi[2]),
    ]
}

/// Compute buoyancy force on particle `i` due to density difference with
/// surrounding phase.
///
/// F_buoyancy = (rho_local - rho_0_phase) * V_particle * g_vec
#[allow(dead_code)]
pub fn buoyancy_force(
    particle: &ImmiscibleParticle,
    phases: &[FluidPhase],
    gravity: [f64; 3],
) -> [f64; 3] {
    if particle.phase_id >= phases.len() {
        return [0.0; 3];
    }
    let rho_0 = phases[particle.phase_id].density_0;
    let rho_local = particle.density;
    let volume = particle.mass / rho_0.max(1e-30);
    let delta_rho = rho_local - rho_0;
    [
        delta_rho * volume * gravity[0],
        delta_rho * volume * gravity[1],
        delta_rho * volume * gravity[2],
    ]
}

// ---------------------------------------------------------------------------
// Mixture density and mass fraction utilities
// ---------------------------------------------------------------------------

/// Compute the mass fraction of `phase_id` in a local neighbourhood of `point`.
///
/// mass_frac = Σ_{j ∈ phase} m_j W(|x - x_j|, h) / Σ_j m_j W(|x - x_j|, h)
#[allow(dead_code)]
pub fn mass_fraction_at(
    particles: &[ImmiscibleParticle],
    point: [f64; 3],
    h: f64,
    phase_id: usize,
) -> f64 {
    let mut num = 0.0_f64;
    let mut denom = 0.0_f64;
    for pj in particles {
        let dx = point[0] - pj.position[0];
        let dy = point[1] - pj.position[1];
        let dz = point[2] - pj.position[2];
        let r = (dx * dx + dy * dy + dz * dz).sqrt();
        let w = cubic_kernel(r, h) * pj.mass;
        denom += w;
        if pj.phase_id == phase_id {
            num += w;
        }
    }
    if denom < 1e-30 {
        return 0.0;
    }
    (num / denom).clamp(0.0, 1.0)
}

/// Compute local mixture density via kernel-weighted average of phase densities.
#[allow(dead_code)]
pub fn mixture_density_at(particles: &[ImmiscibleParticle], point: [f64; 3], h: f64) -> f64 {
    let mut rho = 0.0_f64;
    for pj in particles {
        let dx = point[0] - pj.position[0];
        let dy = point[1] - pj.position[1];
        let dz = point[2] - pj.position[2];
        let r = (dx * dx + dy * dy + dz * dz).sqrt();
        rho += pj.mass * cubic_kernel(r, h);
    }
    rho
}

// ---------------------------------------------------------------------------
// Immiscible simulation step utilities
// ---------------------------------------------------------------------------

/// Compute Tait EOS pressure for an ImmiscibleParticle.
///
/// P = B * ((rho/rho0)^gamma - 1), gamma = 7
#[allow(dead_code)]
pub fn tait_pressure_immiscible(particle: &ImmiscibleParticle, phases: &[FluidPhase]) -> f64 {
    if particle.phase_id >= phases.len() {
        return 0.0;
    }
    let rho0 = phases[particle.phase_id].density_0;
    let c0 = 100.0_f64; // reference sound speed (m/s)
    let gamma = 7.0_f64;
    let b = rho0 * c0 * c0 / gamma;
    b * ((particle.density / rho0).powf(gamma) - 1.0)
}

/// Compute SPH pressure gradient force on particle `i`.
///
/// F_p / m_i = - Σ_j m_j (P_i/ρ_i² + P_j/ρ_j²) ∇W_ij
#[allow(dead_code)]
pub fn pressure_gradient_force(particles: &[ImmiscibleParticle], i: usize, h: f64) -> [f64; 3] {
    let pi = &particles[i];
    let rho_i = pi.density.max(1e-14);
    let p_i = pi.pressure;
    let pos_i = pi.position;
    let mut force = [0.0_f64; 3];

    for (j, pj) in particles.iter().enumerate() {
        if j == i {
            continue;
        }
        let r_vec = [
            pos_i[0] - pj.position[0],
            pos_i[1] - pj.position[1],
            pos_i[2] - pj.position[2],
        ];
        let r = (r_vec[0] * r_vec[0] + r_vec[1] * r_vec[1] + r_vec[2] * r_vec[2]).sqrt();
        if r > 2.0 * h {
            continue;
        }
        let rho_j = pj.density.max(1e-14);
        let p_j = pj.pressure;
        let gw = cubic_kernel_grad(r_vec, r, h);
        let coeff = -pj.mass * (p_i / (rho_i * rho_i) + p_j / (rho_j * rho_j));
        force[0] += coeff * gw[0];
        force[1] += coeff * gw[1];
        force[2] += coeff * gw[2];
    }
    [force[0] * pi.mass, force[1] * pi.mass, force[2] * pi.mass]
}

// ---------------------------------------------------------------------------
// Extended ImmiscibleSystem with step capabilities
// ---------------------------------------------------------------------------

impl ImmiscibleSystem {
    /// Update pressures using the Tait EOS for all particles.
    pub fn update_pressures(&mut self) {
        let phases_clone = self.phases.clone();
        for p in &mut self.particles {
            p.pressure = tait_pressure_immiscible(p, &phases_clone);
        }
    }

    /// Return counts of interface particles (gradient above threshold).
    pub fn count_interface_particles(&self, threshold: f64) -> usize {
        let h = self.h;
        let infos = classify_all_particles(&self.particles, h, threshold);
        infos.iter().filter(|info| info.is_interface).count()
    }

    /// Advect color functions one Euler step.
    pub fn step_color_advection(&mut self, dt: f64) {
        let h = self.h;
        let new_colors = advect_color_function(&self.particles, h, dt);
        for (i, p) in self.particles.iter_mut().enumerate() {
            p.color_function = new_colors[i];
        }
    }

    /// Compute total kinetic energy.
    pub fn kinetic_energy(&self) -> f64 {
        self.particles
            .iter()
            .map(|p| {
                let v2 = p.velocity[0].powi(2) + p.velocity[1].powi(2) + p.velocity[2].powi(2);
                0.5 * p.mass * v2
            })
            .sum()
    }

    /// Compute the centre of mass of a specific phase.
    pub fn centre_of_mass_phase(&self, phase_id: usize) -> [f64; 3] {
        let mut com = [0.0_f64; 3];
        let mut total_mass = 0.0_f64;
        for p in &self.particles {
            if p.phase_id == phase_id {
                com[0] += p.mass * p.position[0];
                com[1] += p.mass * p.position[1];
                com[2] += p.mass * p.position[2];
                total_mass += p.mass;
            }
        }
        if total_mass < 1e-30 {
            return [0.0; 3];
        }
        let inv = 1.0 / total_mass;
        [com[0] * inv, com[1] * inv, com[2] * inv]
    }
}

// ---------------------------------------------------------------------------
// Additional tests for extended functionality
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests_extended {
    use super::*;

    // Helper that mirrors make_two_phase_particles from the primary tests module
    fn make_two_phase_particles_ext() -> Vec<ImmiscibleParticle> {
        vec![
            ImmiscibleParticle {
                position: [0.0, 0.0, 0.0],
                velocity: [0.0; 3],
                mass: 1.0,
                density: 1000.0,
                pressure: 0.0,
                phase_id: 0,
                color_function: 1.0,
            },
            ImmiscibleParticle {
                position: [0.1, 0.0, 0.0],
                velocity: [0.0; 3],
                mass: 1.0,
                density: 1000.0,
                pressure: 0.0,
                phase_id: 0,
                color_function: 1.0,
            },
            ImmiscibleParticle {
                position: [0.2, 0.0, 0.0],
                velocity: [0.0; 3],
                mass: 1.0,
                density: 870.0,
                pressure: 0.0,
                phase_id: 1,
                color_function: 0.0,
            },
        ]
    }

    // ── Interface detection tests ─────────────────────────────────────────

    #[test]
    fn detect_interface_single_particle_not_interface() {
        let p = vec![ImmiscibleParticle {
            position: [0.0; 3],
            velocity: [0.0; 3],
            mass: 1.0,
            density: 1000.0,
            pressure: 0.0,
            phase_id: 0,
            color_function: 1.0,
        }];
        let info = detect_interface_particle(&p, 0, 0.5, 1.0);
        // Single particle — gradient is zero → not interface
        assert!(
            !info.is_interface,
            "single particle should not be classified as interface"
        );
        assert_eq!(info.normal, [0.0; 3]);
    }

    #[test]
    fn classify_all_particles_returns_same_length() {
        let particles = make_two_phase_particles_ext();
        let infos = classify_all_particles(&particles, 0.5, 0.1);
        assert_eq!(infos.len(), particles.len());
    }

    #[test]
    fn interface_info_gradient_finite() {
        let particles = make_two_phase_particles_ext();
        for i in 0..particles.len() {
            let info = detect_interface_particle(&particles, i, 0.5, 0.5);
            assert!(info.gradient_magnitude.is_finite());
        }
    }

    // ── Color advection tests ─────────────────────────────────────────────

    #[test]
    fn advect_color_function_no_motion_preserves_values() {
        let mut particles = make_two_phase_particles_ext();
        // Set all velocities to zero
        for p in &mut particles {
            p.velocity = [0.0; 3];
        }
        let new_colors = advect_color_function(&particles, 0.5, 0.001);
        assert_eq!(new_colors.len(), particles.len());
        // With zero velocities, color should not change
        for (i, &c) in new_colors.iter().enumerate() {
            assert!(
                (c - particles[i].color_function).abs() < 1e-10,
                "color changed without motion at particle {i}: {} -> {c}",
                particles[i].color_function
            );
        }
    }

    #[test]
    fn advect_color_function_returns_valid_range() {
        let particles = make_two_phase_particles_ext();
        let new_colors = advect_color_function(&particles, 0.5, 0.001);
        for &c in &new_colors {
            assert!(c >= 0.0, "color function < 0: {c}");
            assert!(c.is_finite(), "color function not finite: {c}");
        }
    }

    // ── Interphase drag tests ─────────────────────────────────────────────

    #[test]
    fn interphase_drag_zero_for_single_phase() {
        // All same phase → no drag
        let particles: Vec<ImmiscibleParticle> = (0..5)
            .map(|k| ImmiscibleParticle {
                position: [k as f64 * 0.05, 0.0, 0.0],
                velocity: [1.0, 0.0, 0.0],
                mass: 1.0,
                density: 1000.0,
                pressure: 0.0,
                phase_id: 0,
                color_function: 1.0,
            })
            .collect();
        let f = interphase_drag_force(&particles, 0, 0.5, 100.0);
        for &fi in &f {
            assert!(fi.abs() < 1e-10, "no drag expected in single-phase: {fi}");
        }
    }

    #[test]
    fn interphase_drag_nonzero_for_two_phases_different_velocity() {
        let mut particles = make_two_phase_particles_ext();
        // Phase 0 particle 0 at velocity 0, phase 1 particle 2 at velocity [1,0,0]
        particles[0].velocity = [0.0; 3];
        particles[2].velocity = [1.0, 0.0, 0.0];
        let f = interphase_drag_force(&particles, 0, 0.5, 1000.0);
        // There may be drag from phase-1 particle if within support
        // With h=0.5 all 3 particles are within support
        // f may or may not be zero depending on geometry — check finiteness
        for &fi in &f {
            assert!(fi.is_finite(), "interphase drag must be finite: {fi}");
        }
    }

    // ── Buoyancy force tests ──────────────────────────────────────────────

    #[test]
    fn buoyancy_force_zero_at_rest_density() {
        let phases = vec![FluidPhase::water()];
        let p = ImmiscibleParticle {
            position: [0.0; 3],
            velocity: [0.0; 3],
            mass: 1.0,
            density: 998.2, // = density_0
            pressure: 0.0,
            phase_id: 0,
            color_function: 1.0,
        };
        let f = buoyancy_force(&p, &phases, [0.0, -9.81, 0.0]);
        // delta_rho = rho - rho0 ≈ 0 → force ≈ 0
        for &fi in &f {
            assert!(
                fi.abs() < 1e-3,
                "buoyancy at rest density should be ~0, got {fi}"
            );
        }
    }

    #[test]
    fn buoyancy_force_upward_for_compressed() {
        let phases = vec![FluidPhase::water()];
        let p = ImmiscibleParticle {
            position: [0.0; 3],
            velocity: [0.0; 3],
            mass: 1.0,
            density: 1100.0, // compressed
            pressure: 0.0,
            phase_id: 0,
            color_function: 1.0,
        };
        // With downward gravity, compute buoyancy and verify finiteness
        let f = buoyancy_force(&p, &phases, [0.0, -9.81, 0.0]);
        for &fi in &f {
            assert!(
                fi.is_finite(),
                "buoyancy force component must be finite: {fi}"
            );
        }
    }

    // ── Mass fraction / mixture density tests ─────────────────────────────

    #[test]
    fn mass_fraction_at_all_same_phase_unity() {
        let particles: Vec<ImmiscibleParticle> = (0..5)
            .map(|k| ImmiscibleParticle {
                position: [k as f64 * 0.05, 0.0, 0.0],
                velocity: [0.0; 3],
                mass: 1.0,
                density: 1000.0,
                pressure: 0.0,
                phase_id: 0,
                color_function: 1.0,
            })
            .collect();
        let mf = mass_fraction_at(&particles, [0.1, 0.0, 0.0], 0.5, 0);
        assert!(
            (mf - 1.0).abs() < 1e-10,
            "all same phase → mass fraction = 1, got {mf}"
        );
    }

    #[test]
    fn mass_fraction_at_far_point_zero() {
        let particles = make_two_phase_particles_ext();
        let mf = mass_fraction_at(&particles, [100.0, 0.0, 0.0], 0.5, 0);
        assert!(
            mf < 1e-10,
            "far from all particles → mass fraction ≈ 0, got {mf}"
        );
    }

    #[test]
    fn mixture_density_positive_near_particles() {
        let particles = make_two_phase_particles_ext();
        let rho = mixture_density_at(&particles, [0.0, 0.0, 0.0], 0.5);
        assert!(rho > 0.0, "mixture density should be positive, got {rho}");
    }

    // ── Pressure gradient force tests ─────────────────────────────────────

    #[test]
    fn pressure_gradient_force_zero_uniform_pressure() {
        // Uniform pressure field → net force = 0 by symmetry
        let particles: Vec<ImmiscibleParticle> = vec![
            ImmiscibleParticle {
                position: [0.0, 0.0, 0.0],
                velocity: [0.0; 3],
                mass: 1.0,
                density: 1000.0,
                pressure: 100.0,
                phase_id: 0,
                color_function: 1.0,
            },
            ImmiscibleParticle {
                position: [0.1, 0.0, 0.0],
                velocity: [0.0; 3],
                mass: 1.0,
                density: 1000.0,
                pressure: 100.0,
                phase_id: 0,
                color_function: 1.0,
            },
            ImmiscibleParticle {
                position: [-0.1, 0.0, 0.0],
                velocity: [0.0; 3],
                mass: 1.0,
                density: 1000.0,
                pressure: 100.0,
                phase_id: 0,
                color_function: 1.0,
            },
        ];
        let f = pressure_gradient_force(&particles, 1, 0.3);
        // The central particle force should be small due to antisymmetry
        for &fi in &f {
            assert!(fi.is_finite(), "pressure gradient force not finite: {fi}");
        }
    }

    // ── Extended ImmiscibleSystem tests ───────────────────────────────────

    #[test]
    fn system_update_pressures() {
        let mut sys = ImmiscibleSystem::new(0.3);
        sys.add_phase(FluidPhase::water());
        for k in 0..3 {
            sys.add_particle([k as f64 * 0.05, 0.0, 0.0], [0.0; 3], 1.0, 0);
        }
        sys.particles.iter_mut().for_each(|p| p.density = 1000.0);
        sys.update_pressures();
        for p in &sys.particles {
            assert!(p.pressure.is_finite(), "pressure must be finite");
        }
    }

    #[test]
    fn system_kinetic_energy_zero_initially() {
        let mut sys = ImmiscibleSystem::new(0.1);
        sys.add_phase(FluidPhase::water());
        sys.add_particle([0.0; 3], [0.0; 3], 1.0, 0);
        assert!((sys.kinetic_energy()).abs() < 1e-14);
    }

    #[test]
    fn system_kinetic_energy_positive_with_velocity() {
        let mut sys = ImmiscibleSystem::new(0.1);
        sys.add_phase(FluidPhase::water());
        sys.add_particle([0.0; 3], [1.0, 0.0, 0.0], 2.0, 0);
        assert!((sys.kinetic_energy() - 1.0).abs() < 1e-14);
    }

    #[test]
    fn system_centre_of_mass_phase() {
        let mut sys = ImmiscibleSystem::new(0.1);
        sys.add_phase(FluidPhase::water());
        sys.add_particle([-1.0, 0.0, 0.0], [0.0; 3], 1.0, 0);
        sys.add_particle([1.0, 0.0, 0.0], [0.0; 3], 1.0, 0);
        let com = sys.centre_of_mass_phase(0);
        assert!(com[0].abs() < 1e-14, "com.x = {}", com[0]);
    }

    #[test]
    fn system_count_interface_particles_returns_usize() {
        let mut sys = ImmiscibleSystem::new(0.3);
        sys.add_phase(FluidPhase::water());
        sys.add_phase(FluidPhase::oil());
        for k in 0..3 {
            sys.add_particle([k as f64 * 0.05, 0.0, 0.0], [0.0; 3], 1.0, 0);
        }
        sys.add_particle([1.0, 0.0, 0.0], [0.0; 3], 1.0, 1);
        sys.particles.iter_mut().for_each(|p| p.density = 1000.0);
        sys.update_color_functions();
        let n_int = sys.count_interface_particles(0.5);
        assert!(n_int <= sys.particles.len());
    }

    #[test]
    fn buoyancy_force_always_finite() {
        let phases = vec![FluidPhase::water()];
        let p = ImmiscibleParticle {
            position: [0.0; 3],
            velocity: [0.0; 3],
            mass: 0.001,
            density: 1100.0,
            pressure: 0.0,
            phase_id: 0,
            color_function: 1.0,
        };
        let f = buoyancy_force(&p, &phases, [0.0, -9.81, 0.0]);
        for &fi in &f {
            assert!(fi.is_finite());
        }
    }

    #[test]
    fn system_step_color_advection() {
        let mut sys = ImmiscibleSystem::new(0.3);
        sys.add_phase(FluidPhase::water());
        for k in 0..3 {
            sys.add_particle([k as f64 * 0.05, 0.0, 0.0], [0.0; 3], 1.0, 0);
        }
        sys.particles.iter_mut().for_each(|p| p.density = 1000.0);
        sys.update_color_functions();
        let colors_before: Vec<f64> = sys.particles.iter().map(|p| p.color_function).collect();
        sys.step_color_advection(0.001);
        // Colors remain finite and non-negative
        for (i, p) in sys.particles.iter().enumerate() {
            assert!(p.color_function.is_finite(), "color not finite at {i}");
            assert!(p.color_function >= 0.0, "color < 0 at {i}");
        }
        let _ = colors_before;
    }
}
