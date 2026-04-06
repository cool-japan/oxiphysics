// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Particle transport in LBM flow fields.
//!
//! Implements one-way coupling of passive particles (tracer/heavy) in an LBM
//! flow field using:
//! - Stokes drag: `F = 3π μ d Δv`
//! - Buoyancy/gravity: `F = (ρ_fluid − ρ_particle) * V * g`
//! - Explicit Euler time integration
//! - Stokes terminal (sedimentation) velocity

#![allow(dead_code)]

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// LbmParticle
// ---------------------------------------------------------------------------

/// A passive particle transported by an LBM flow field.
#[derive(Clone, Debug, PartialEq)]
pub struct LbmParticle {
    /// Position \[x, y\] in lattice units.
    pub position: [f64; 2],
    /// Particle velocity \[ux, uy\].
    pub velocity: [f64; 2],
    /// Particle diameter (lattice units).
    pub diameter: f64,
    /// Ratio ρ_particle / ρ_fluid.
    pub density_ratio: f64,
}

impl LbmParticle {
    /// Create a new `LbmParticle`.
    ///
    /// # Arguments
    /// * `x`, `y` — initial position
    /// * `diameter` — particle diameter in lattice units
    /// * `density_ratio` — ρ_p / ρ_f
    pub fn new(x: f64, y: f64, diameter: f64, density_ratio: f64) -> Self {
        Self {
            position: [x, y],
            velocity: [0.0, 0.0],
            diameter,
            density_ratio,
        }
    }

    /// Particle radius.
    #[inline]
    pub fn radius(&self) -> f64 {
        self.diameter * 0.5
    }

    /// Particle volume (2D: area of a circle with the given diameter).
    #[inline]
    pub fn volume_2d(&self) -> f64 {
        PI * self.radius() * self.radius()
    }
}

// ---------------------------------------------------------------------------
// Stokes drag
// ---------------------------------------------------------------------------

/// Stokes drag force on a spherical particle.
///
/// `F = 3π μ d (u_fluid − u_particle)`
///
/// # Arguments
/// * `mu` — dynamic viscosity of the fluid
/// * `d` — particle diameter
/// * `u_fluid` — fluid velocity at particle position \[ux, uy\]
/// * `u_particle` — particle velocity \[ux, uy\]
///
/// Returns force \[fx, fy\].
pub fn drag_force_stokes(mu: f64, d: f64, u_fluid: [f64; 2], u_particle: [f64; 2]) -> [f64; 2] {
    let coeff = 3.0 * PI * mu * d;
    [
        coeff * (u_fluid[0] - u_particle[0]),
        coeff * (u_fluid[1] - u_particle[1]),
    ]
}

// ---------------------------------------------------------------------------
// Buoyancy force
// ---------------------------------------------------------------------------

/// Buoyancy (gravity) force on a particle in a 2D LBM simulation.
///
/// `F_y = (ρ_fluid − ρ_particle) * V * g`
///
/// # Arguments
/// * `rho_fluid` — fluid reference density
/// * `rho_particle` — particle density
/// * `volume` — particle volume (use `LbmParticle::volume_2d()`)
/// * `gravity` — gravitational acceleration (positive downward)
///
/// Returns \[0.0, F_y\] — force acts only in the y (vertical) direction.
pub fn buoyancy_force(rho_fluid: f64, rho_particle: f64, volume: f64, gravity: f64) -> [f64; 2] {
    let fy = (rho_fluid - rho_particle) * volume * gravity;
    [0.0, fy]
}

// ---------------------------------------------------------------------------
// Terminal sedimentation velocity (Stokes)
// ---------------------------------------------------------------------------

/// Terminal settling velocity for a sphere under Stokes drag and gravity.
///
/// `v_s = (ρ_particle − ρ_fluid) * d² * g / (18 μ)`
///
/// # Arguments
/// * `rho_particle` — particle density
/// * `rho_fluid` — fluid density
/// * `d` — particle diameter
/// * `mu` — dynamic viscosity
/// * `gravity` — gravitational acceleration
///
/// Returns the magnitude of the terminal velocity (positive = downward settling
/// when ρ_particle > ρ_fluid).
pub fn particle_sedimentation_velocity(
    rho_particle: f64,
    rho_fluid: f64,
    d: f64,
    mu: f64,
    gravity: f64,
) -> f64 {
    (rho_particle - rho_fluid) * d * d * gravity / (18.0 * mu)
}

// ---------------------------------------------------------------------------
// Bilinear interpolation helper
// ---------------------------------------------------------------------------

/// Bilinear interpolation of a 2D grid field at a sub-grid position.
///
/// # Arguments
/// * `field` — row-major grid of shape `ny × nx` (index `iy * nx + ix`)
/// * `nx`, `ny` — grid dimensions
/// * `x`, `y` — query position in grid coordinates
///
/// Returns the interpolated value clamped to the grid interior.
pub fn bilinear_interp(field: &[f64], nx: usize, ny: usize, x: f64, y: f64) -> f64 {
    let x = x.clamp(0.0, (nx - 1) as f64);
    let y = y.clamp(0.0, (ny - 1) as f64);
    let ix = x.floor() as usize;
    let iy = y.floor() as usize;
    let ix1 = (ix + 1).min(nx - 1);
    let iy1 = (iy + 1).min(ny - 1);
    let tx = x - ix as f64;
    let ty = y - iy as f64;
    let v00 = field[iy * nx + ix];
    let v10 = field[iy * nx + ix1];
    let v01 = field[iy1 * nx + ix];
    let v11 = field[iy1 * nx + ix1];
    v00 * (1.0 - tx) * (1.0 - ty) + v10 * tx * (1.0 - ty) + v01 * (1.0 - tx) * ty + v11 * tx * ty
}

// ---------------------------------------------------------------------------
// Single-particle advection step
// ---------------------------------------------------------------------------

/// Advance a single particle by one explicit Euler step.
///
/// The particle is subject to Stokes drag and buoyancy forces.
/// Particle mass is approximated as `ρ_particle * volume`.
///
/// # Arguments
/// * `particle` — mutable particle to update
/// * `ux_grid`, `uy_grid` — LBM velocity grids (row-major, `ny × nx`)
/// * `nx`, `ny` — grid dimensions
/// * `rho_fluid` — fluid reference density
/// * `mu` — dynamic viscosity
/// * `gravity` — gravitational acceleration (y-direction)
/// * `dt` — time step
#[allow(clippy::too_many_arguments)]
pub fn advect_particle(
    particle: &mut LbmParticle,
    ux_grid: &[f64],
    uy_grid: &[f64],
    nx: usize,
    ny: usize,
    rho_fluid: f64,
    mu: f64,
    gravity: f64,
    dt: f64,
) {
    let x = particle.position[0];
    let y = particle.position[1];
    let uf = [
        bilinear_interp(ux_grid, nx, ny, x, y),
        bilinear_interp(uy_grid, nx, ny, x, y),
    ];
    let drag = drag_force_stokes(mu, particle.diameter, uf, particle.velocity);
    let buoy = buoyancy_force(
        rho_fluid,
        rho_fluid * particle.density_ratio,
        particle.volume_2d(),
        gravity,
    );
    let rho_p = rho_fluid * particle.density_ratio;
    let mass = rho_p * particle.volume_2d();
    let inv_mass = if mass > 0.0 { 1.0 / mass } else { 0.0 };
    let ax = (drag[0] + buoy[0]) * inv_mass;
    let ay = (drag[1] + buoy[1]) * inv_mass;
    particle.velocity[0] += dt * ax;
    particle.velocity[1] += dt * ay;
    particle.position[0] += dt * particle.velocity[0];
    particle.position[1] += dt * particle.velocity[1];
}

// ---------------------------------------------------------------------------
// Multi-particle tracking
// ---------------------------------------------------------------------------

/// Advance all particles by one time step.
///
/// # Arguments
/// * `particles` — mutable slice of particles
/// * `ux_grid`, `uy_grid` — LBM velocity grids
/// * `nx`, `ny` — grid dimensions
/// * `rho_fluid` — fluid reference density
/// * `mu` — dynamic viscosity
/// * `gravity` — gravitational acceleration
/// * `dt` — time step
#[allow(clippy::too_many_arguments)]
pub fn track_particles(
    particles: &mut [LbmParticle],
    ux_grid: &[f64],
    uy_grid: &[f64],
    nx: usize,
    ny: usize,
    rho_fluid: f64,
    mu: f64,
    gravity: f64,
    dt: f64,
) {
    for p in particles.iter_mut() {
        advect_particle(p, ux_grid, uy_grid, nx, ny, rho_fluid, mu, gravity, dt);
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- LbmParticle ---

    #[test]
    fn test_particle_new() {
        let p = LbmParticle::new(1.0, 2.0, 0.5, 2.0);
        assert_eq!(p.position, [1.0, 2.0]);
        assert_eq!(p.velocity, [0.0, 0.0]);
        assert_eq!(p.diameter, 0.5);
        assert_eq!(p.density_ratio, 2.0);
    }

    #[test]
    fn test_particle_radius() {
        let p = LbmParticle::new(0.0, 0.0, 4.0, 1.0);
        assert!((p.radius() - 2.0).abs() < 1e-14);
    }

    #[test]
    fn test_particle_volume_2d() {
        let p = LbmParticle::new(0.0, 0.0, 2.0, 1.0);
        // V = π r² = π
        assert!((p.volume_2d() - PI).abs() < 1e-12);
    }

    #[test]
    fn test_particle_clone() {
        let p = LbmParticle::new(3.0, 4.0, 1.0, 1.5);
        let p2 = p.clone();
        assert_eq!(p, p2);
    }

    // --- drag_force_stokes ---

    #[test]
    fn test_stokes_drag_zero_relative_velocity() {
        let f = drag_force_stokes(1.0, 1.0, [1.0, 2.0], [1.0, 2.0]);
        assert_eq!(f, [0.0, 0.0]);
    }

    #[test]
    fn test_stokes_drag_positive_x() {
        // F = 3π μ d Δu; with μ=1, d=1, Δu=1: F = 3π
        let f = drag_force_stokes(1.0, 1.0, [2.0, 0.0], [1.0, 0.0]);
        assert!((f[0] - 3.0 * PI).abs() < 1e-12, "fx = {}", f[0]);
        assert!(f[1].abs() < 1e-14);
    }

    #[test]
    fn test_stokes_drag_scales_with_viscosity() {
        let f1 = drag_force_stokes(1.0, 1.0, [1.0, 0.0], [0.0, 0.0]);
        let f2 = drag_force_stokes(2.0, 1.0, [1.0, 0.0], [0.0, 0.0]);
        assert!((f2[0] / f1[0] - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_stokes_drag_scales_with_diameter() {
        let f1 = drag_force_stokes(1.0, 1.0, [1.0, 0.0], [0.0, 0.0]);
        let f2 = drag_force_stokes(1.0, 3.0, [1.0, 0.0], [0.0, 0.0]);
        assert!((f2[0] / f1[0] - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_stokes_drag_negative_relative() {
        let f = drag_force_stokes(1.0, 1.0, [0.0, 0.0], [1.0, 0.0]);
        assert!(f[0] < 0.0);
    }

    // --- buoyancy_force ---

    #[test]
    fn test_buoyancy_no_x_component() {
        let f = buoyancy_force(1.0, 0.5, 1.0, 9.8);
        assert_eq!(f[0], 0.0);
    }

    #[test]
    fn test_buoyancy_lighter_particle_floats() {
        // ρ_fluid > ρ_particle → upward force → F_y > 0 for positive gravity
        let f = buoyancy_force(1.0, 0.5, 1.0, 1.0);
        assert!(f[1] > 0.0);
    }

    #[test]
    fn test_buoyancy_heavier_particle_sinks() {
        // ρ_particle > ρ_fluid → F_y < 0
        let f = buoyancy_force(1.0, 2.0, 1.0, 1.0);
        assert!(f[1] < 0.0);
    }

    #[test]
    fn test_buoyancy_neutral() {
        let f = buoyancy_force(1.0, 1.0, 1.0, 9.8);
        assert_eq!(f[1], 0.0);
    }

    // --- terminal sedimentation velocity ---

    #[test]
    fn test_sedimentation_velocity_positive_for_heavy_particle() {
        let v = particle_sedimentation_velocity(2.0, 1.0, 1.0, 1.0, 1.0);
        assert!(v > 0.0);
    }

    #[test]
    fn test_sedimentation_velocity_zero_for_neutral() {
        let v = particle_sedimentation_velocity(1.0, 1.0, 1.0, 1.0, 9.8);
        assert_eq!(v, 0.0);
    }

    #[test]
    fn test_sedimentation_velocity_stokes_formula() {
        // v_s = (ρ_p - ρ_f) d² g / (18 μ)
        // With ρ_p=2, ρ_f=1, d=0.1, g=9.8, μ=0.001 → v_s = 1*0.01*9.8/18*0.001 = 5.44
        let v = particle_sedimentation_velocity(2.0, 1.0, 0.1, 0.001, 9.8);
        let expected = 1.0 * 0.01 * 9.8 / (18.0 * 0.001);
        assert!((v - expected).abs() < 1e-12, "v={v}, expected={expected}");
    }

    #[test]
    fn test_sedimentation_velocity_scales_d_squared() {
        let v1 = particle_sedimentation_velocity(2.0, 1.0, 1.0, 1.0, 1.0);
        let v2 = particle_sedimentation_velocity(2.0, 1.0, 2.0, 1.0, 1.0);
        assert!((v2 / v1 - 4.0).abs() < 1e-12);
    }

    // --- bilinear interpolation ---

    #[test]
    fn test_bilinear_interp_at_node() {
        let nx = 3;
        let ny = 3;
        let field: Vec<f64> = (0..9).map(|i| i as f64).collect();
        let v = bilinear_interp(&field, nx, ny, 1.0, 1.0);
        assert!((v - field[nx + 1]).abs() < 1e-12);
    }

    #[test]
    fn test_bilinear_interp_midpoint() {
        let nx = 2;
        let ny = 2;
        let field = vec![0.0, 1.0, 1.0, 2.0];
        let v = bilinear_interp(&field, nx, ny, 0.5, 0.5);
        assert!((v - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_bilinear_interp_clamped() {
        let nx = 3;
        let ny = 3;
        let field = vec![1.0_f64; 9];
        // Out-of-bounds query should clamp.
        let v = bilinear_interp(&field, nx, ny, 10.0, 10.0);
        assert!((v - 1.0).abs() < 1e-12);
    }

    // --- advect_particle ---

    #[test]
    fn test_advect_particle_moves_with_flow() {
        let nx = 10;
        let ny = 10;
        let ux = vec![0.1_f64; nx * ny];
        let uy = vec![0.0_f64; nx * ny];
        let mut p = LbmParticle::new(5.0, 5.0, 0.5, 1.001);
        let x0 = p.position[0];
        advect_particle(&mut p, &ux, &uy, nx, ny, 1.0, 0.01, 0.0, 1.0);
        // Particle should move in x.
        assert!(p.position[0] != x0 || p.position[1] != 5.0);
    }

    #[test]
    fn test_advect_particle_no_flow_no_gravity() {
        let nx = 5;
        let ny = 5;
        let ux = vec![0.0_f64; nx * ny];
        let uy = vec![0.0_f64; nx * ny];
        let mut p = LbmParticle::new(2.0, 2.0, 0.5, 1.0);
        advect_particle(&mut p, &ux, &uy, nx, ny, 1.0, 1.0, 0.0, 0.1);
        // Neutral density + no flow → no motion.
        assert_eq!(p.position, [2.0, 2.0]);
    }

    // --- track_particles ---

    #[test]
    fn test_track_particles_does_not_panic() {
        let nx = 8;
        let ny = 8;
        let ux = vec![0.05_f64; nx * ny];
        let uy = vec![0.0_f64; nx * ny];
        let mut particles = vec![
            LbmParticle::new(3.0, 4.0, 0.5, 1.5),
            LbmParticle::new(5.0, 5.0, 0.3, 0.8),
        ];
        track_particles(&mut particles, &ux, &uy, nx, ny, 1.0, 0.01, 0.0, 0.1);
    }

    #[test]
    fn test_track_particles_empty() {
        let nx = 4;
        let ny = 4;
        let ux = vec![0.0_f64; nx * ny];
        let uy = vec![0.0_f64; nx * ny];
        let mut particles: Vec<LbmParticle> = vec![];
        // Should not panic on empty slice.
        track_particles(&mut particles, &ux, &uy, nx, ny, 1.0, 1.0, 0.0, 0.1);
    }
}
