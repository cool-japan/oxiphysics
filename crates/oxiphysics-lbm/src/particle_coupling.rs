#![allow(clippy::needless_range_loop)]
#![allow(clippy::manual_range_contains)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Particle-fluid coupling for the Lattice Boltzmann Method.
//!
//! Provides two-way coupling between discrete particles and the LBM fluid,
//! including drag forces, void fraction corrections, wall reflections, and
//! hard-sphere collisions.

use std::f64::consts::PI;

// ─────────────────────────────────────────────────────────────────────────────
// Core data structures
// ─────────────────────────────────────────────────────────────────────────────

/// A coupled particle tracked in the fluid domain.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct CoupledParticle {
    /// Position \[x, y\] in lattice units.
    pub position: [f64; 2],
    /// Velocity \[ux, uy\] in lattice units / step.
    pub velocity: [f64; 2],
    /// Net force \[fx, fy\] acting on the particle.
    pub force: [f64; 2],
    /// Particle mass.
    pub mass: f64,
    /// Particle radius (lattice units).
    pub radius: f64,
    /// Unique identifier.
    pub id: usize,
}

impl CoupledParticle {
    /// Create a new coupled particle.
    #[allow(dead_code)]
    pub fn new(position: [f64; 2], velocity: [f64; 2], mass: f64, radius: f64, id: usize) -> Self {
        Self {
            position,
            velocity,
            force: [0.0, 0.0],
            mass,
            radius,
            id,
        }
    }
}

/// Macroscopic fluid state at a single cell.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct FluidCell {
    /// Fluid density.
    pub rho: f64,
    /// x-component of fluid velocity.
    pub ux: f64,
    /// y-component of fluid velocity.
    pub uy: f64,
}

impl FluidCell {
    /// Create a new fluid cell state.
    #[allow(dead_code)]
    pub fn new(rho: f64, ux: f64, uy: f64) -> Self {
        Self { rho, ux, uy }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Drag force models
// ─────────────────────────────────────────────────────────────────────────────

/// Collection of drag force models.
#[allow(dead_code)]
pub struct DragForce;

impl DragForce {
    /// Stokes drag: **F = 6 π μ r (u_fluid − u_particle)**.
    ///
    /// Valid for low Reynolds number (Re ≪ 1).
    #[allow(dead_code)]
    pub fn stokes_drag(particle: &CoupledParticle, fluid: &FluidCell, viscosity: f64) -> [f64; 2] {
        let factor = 6.0 * PI * viscosity * particle.radius;
        let fx = factor * (fluid.ux - particle.velocity[0]);
        let fy = factor * (fluid.uy - particle.velocity[1]);
        [fx, fy]
    }

    /// Particle Reynolds number: **Re = 2 r |u_slip| / ν**.
    #[allow(dead_code)]
    pub fn reynolds_number(particle: &CoupledParticle, fluid: &FluidCell, viscosity: f64) -> f64 {
        let dux = fluid.ux - particle.velocity[0];
        let duy = fluid.uy - particle.velocity[1];
        let slip_speed = (dux * dux + duy * duy).sqrt();
        let diameter = 2.0 * particle.radius;
        if viscosity.abs() < f64::EPSILON {
            0.0
        } else {
            diameter * slip_speed / viscosity
        }
    }

    /// Schiller–Naumann drag force.
    ///
    /// Drag coefficient:
    /// - Re < 1000: `Cd = 24/Re * (1 + 0.15 Re^0.687)`
    /// - Re ≥ 1000: `Cd = 0.44`
    ///
    /// Force: `F = 0.5 * Cd * rho_f * A * |slip|² * slip_hat`
    #[allow(dead_code)]
    pub fn schiller_naumann_drag(
        particle: &CoupledParticle,
        fluid: &FluidCell,
        viscosity: f64,
        rho_f: f64,
    ) -> [f64; 2] {
        let dux = fluid.ux - particle.velocity[0];
        let duy = fluid.uy - particle.velocity[1];
        let slip_speed = (dux * dux + duy * duy).sqrt();

        if slip_speed < f64::EPSILON {
            return [0.0, 0.0];
        }

        let re = Self::reynolds_number(particle, fluid, viscosity);

        let cd = if re < 1000.0 {
            if re < f64::EPSILON {
                // Stokes limit: Cd → 24/Re but avoid divide-by-zero
                return Self::stokes_drag(particle, fluid, viscosity);
            }
            24.0 / re * (1.0 + 0.15 * re.powf(0.687))
        } else {
            0.44
        };

        // Cross-sectional area of a sphere (circle in 2-D projection).
        let area = PI * particle.radius * particle.radius;
        let magnitude = 0.5 * cd * rho_f * area * slip_speed * slip_speed;

        let fx = magnitude * dux / slip_speed;
        let fy = magnitude * duy / slip_speed;
        [fx, fy]
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Two-way coupling
// ─────────────────────────────────────────────────────────────────────────────

/// Two-way particle–fluid coupling on a 2-D Cartesian grid.
#[allow(dead_code)]
pub struct TwoWayCoupling {
    /// Number of grid cells in x.
    pub nx: usize,
    /// Number of grid cells in y.
    pub ny: usize,
    /// Grid spacing (lattice units).
    pub dx: f64,
}

impl TwoWayCoupling {
    /// Create a new coupling manager.
    #[allow(dead_code)]
    pub fn new(nx: usize, ny: usize, dx: f64) -> Self {
        Self { nx, ny, dx }
    }

    /// Bilinear interpolation of the fluid velocity at a continuous position.
    ///
    /// `fluid_u` is a flat `nx × ny` array of `(ux, uy)` pairs indexed as
    /// `[iy * nx + ix]`.  Positions outside the domain are clamped.
    #[allow(dead_code)]
    pub fn interpolate_fluid_velocity(&self, fluid_u: &[(f64, f64)], pos: [f64; 2]) -> [f64; 2] {
        let nx = self.nx as f64;
        let ny = self.ny as f64;

        // Convert position to cell index space.
        let xf = (pos[0] / self.dx).clamp(0.0, nx - 1.0 - f64::EPSILON);
        let yf = (pos[1] / self.dx).clamp(0.0, ny - 1.0 - f64::EPSILON);

        let ix0 = xf.floor() as usize;
        let iy0 = yf.floor() as usize;
        let ix1 = (ix0 + 1).min(self.nx - 1);
        let iy1 = (iy0 + 1).min(self.ny - 1);

        let tx = xf - ix0 as f64;
        let ty = yf - iy0 as f64;

        let idx = |ix: usize, iy: usize| iy * self.nx + ix;

        let (u00x, u00y) = fluid_u[idx(ix0, iy0)];
        let (u10x, u10y) = fluid_u[idx(ix1, iy0)];
        let (u01x, u01y) = fluid_u[idx(ix0, iy1)];
        let (u11x, u11y) = fluid_u[idx(ix1, iy1)];

        let ux = (1.0 - tx) * (1.0 - ty) * u00x
            + tx * (1.0 - ty) * u10x
            + (1.0 - tx) * ty * u01x
            + tx * ty * u11x;

        let uy = (1.0 - tx) * (1.0 - ty) * u00y
            + tx * (1.0 - ty) * u10y
            + (1.0 - tx) * ty * u01y
            + tx * ty * u11y;

        [ux, uy]
    }

    /// Distribute particle reaction forces to the nearest grid cell.
    ///
    /// Returns a flat `nx × ny` array of `[fx, fy]` force densities.
    #[allow(dead_code)]
    pub fn spread_particle_force(
        &self,
        particles: &[CoupledParticle],
        nx: usize,
        ny: usize,
    ) -> Vec<[f64; 2]> {
        let mut grid_force = vec![[0.0_f64; 2]; nx * ny];

        for p in particles {
            let ix = ((p.position[0] / self.dx).round() as usize).min(nx - 1);
            let iy = ((p.position[1] / self.dx).round() as usize).min(ny - 1);
            let k = iy * nx + ix;
            // Reaction: fluid receives −force (Newton's 3rd law).
            grid_force[k][0] -= p.force[0];
            grid_force[k][1] -= p.force[1];
        }

        grid_force
    }

    /// Advance a particle using the semi-implicit Euler scheme.
    ///
    /// Velocity update: `v^{n+1} = v^n + (drag + gravity) / m * dt`
    /// Position update: `x^{n+1} = x^n + v^{n+1} * dt`
    #[allow(dead_code)]
    pub fn update_particle(
        &self,
        particle: &mut CoupledParticle,
        drag: [f64; 2],
        gravity: [f64; 2],
        dt: f64,
    ) {
        let ax = (drag[0] + gravity[0]) / particle.mass;
        let ay = (drag[1] + gravity[1]) / particle.mass;

        particle.velocity[0] += ax * dt;
        particle.velocity[1] += ay * dt;

        particle.position[0] += particle.velocity[0] * dt;
        particle.position[1] += particle.velocity[1] * dt;

        // Store the total force on the particle for two-way feedback.
        particle.force[0] = drag[0] + gravity[0];
        particle.force[1] = drag[1] + gravity[1];
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Void fraction
// ─────────────────────────────────────────────────────────────────────────────

/// Void-fraction (fluid volume fraction) utilities.
#[allow(dead_code)]
pub struct VoidFraction;

impl VoidFraction {
    /// Compute the void fraction field (fraction of fluid per cell).
    ///
    /// Each particle contributes its projected area `π r²` to the cell it
    /// occupies.  The void fraction is `ε = 1 − solid_fraction`, clamped to
    /// `[0, 1]`.
    #[allow(dead_code)]
    pub fn compute_void_fraction(
        particles: &[CoupledParticle],
        nx: usize,
        ny: usize,
        dx: f64,
    ) -> Vec<f64> {
        let cell_area = dx * dx;
        let mut solid_frac = vec![0.0_f64; nx * ny];

        for p in particles {
            let ix = ((p.position[0] / dx).round() as usize).min(nx - 1);
            let iy = ((p.position[1] / dx).round() as usize).min(ny - 1);
            let k = iy * nx + ix;
            let particle_area = PI * p.radius * p.radius;
            solid_frac[k] += particle_area / cell_area;
        }

        solid_frac
            .into_iter()
            .map(|sf| (1.0 - sf).clamp(0.0, 1.0))
            .collect()
    }

    /// Apply the Wen–Yu void-fraction correction to drag: **F_corr = F · ε^{−2.65}**.
    #[allow(dead_code)]
    pub fn corrected_drag(drag: [f64; 2], void_fraction: f64) -> [f64; 2] {
        let eps = void_fraction.clamp(f64::EPSILON, 1.0);
        let factor = eps.powf(-2.65);
        [drag[0] * factor, drag[1] * factor]
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Particle collisions
// ─────────────────────────────────────────────────────────────────────────────

/// Particle–wall and particle–particle collision handling.
#[allow(dead_code)]
pub struct ParticleCollision;

impl ParticleCollision {
    /// Return `true` if the particle overlaps with any domain wall.
    ///
    /// Domain walls are at `x = 0`, `x = (nx-1)*dx`, `y = 0`, and
    /// `y = (ny-1)*dx`.
    #[allow(dead_code)]
    pub fn detect_particle_wall(particle: &CoupledParticle, nx: usize, ny: usize, dx: f64) -> bool {
        let r = particle.radius;
        let x = particle.position[0];
        let y = particle.position[1];
        let x_max = (nx - 1) as f64 * dx;
        let y_max = (ny - 1) as f64 * dx;

        x < r || x > x_max - r || y < r || y > y_max - r
    }

    /// Reflect the particle from domain walls with a coefficient of restitution.
    ///
    /// The particle velocity is reversed in the penetrated direction, and the
    /// position is corrected so the particle sits tangent to the wall.
    #[allow(dead_code)]
    pub fn reflect_from_wall(
        particle: &mut CoupledParticle,
        nx: usize,
        ny: usize,
        dx: f64,
        restitution: f64,
    ) {
        let r = particle.radius;
        let x_max = (nx - 1) as f64 * dx;
        let y_max = (ny - 1) as f64 * dx;

        // Left wall
        if particle.position[0] < r {
            particle.position[0] = r;
            particle.velocity[0] = particle.velocity[0].abs() * restitution;
        }
        // Right wall
        if particle.position[0] > x_max - r {
            particle.position[0] = x_max - r;
            particle.velocity[0] = -particle.velocity[0].abs() * restitution;
        }
        // Bottom wall
        if particle.position[1] < r {
            particle.position[1] = r;
            particle.velocity[1] = particle.velocity[1].abs() * restitution;
        }
        // Top wall
        if particle.position[1] > y_max - r {
            particle.position[1] = y_max - r;
            particle.velocity[1] = -particle.velocity[1].abs() * restitution;
        }
    }

    /// Hard-sphere elastic (or inelastic) collision between two particles.
    ///
    /// Only applies if the particles overlap (`|r12| < r1 + r2`).
    /// Uses the 1-D impulse formula along the line of centres, with
    /// coefficient of restitution `e`.
    #[allow(dead_code)]
    pub fn particle_particle_collision(
        p1: &mut CoupledParticle,
        p2: &mut CoupledParticle,
        restitution: f64,
    ) {
        let dx = p2.position[0] - p1.position[0];
        let dy = p2.position[1] - p1.position[1];
        let dist = (dx * dx + dy * dy).sqrt();
        let min_dist = p1.radius + p2.radius;

        if dist >= min_dist || dist < f64::EPSILON {
            return;
        }

        // Unit normal from p1 → p2.
        let nx = dx / dist;
        let ny = dy / dist;

        // Relative velocity along the normal.
        let dvx = p1.velocity[0] - p2.velocity[0];
        let dvy = p1.velocity[1] - p2.velocity[1];
        let vrel_n = dvx * nx + dvy * ny;

        // Only resolve if approaching.
        if vrel_n <= 0.0 {
            return;
        }

        let m1 = p1.mass;
        let m2 = p2.mass;
        let impulse = (1.0 + restitution) * vrel_n / (1.0 / m1 + 1.0 / m2);

        p1.velocity[0] -= impulse / m1 * nx;
        p1.velocity[1] -= impulse / m1 * ny;
        p2.velocity[0] += impulse / m2 * nx;
        p2.velocity[1] += impulse / m2 * ny;

        // Positional correction: separate overlapping particles.
        let overlap = min_dist - dist;
        let correction = overlap / 2.0;
        p1.position[0] -= correction * nx;
        p1.position[1] -= correction * ny;
        p2.position[0] += correction * nx;
        p2.position[1] += correction * ny;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Analytical utilities
// ─────────────────────────────────────────────────────────────────────────────

/// Stokes settling velocity in a quiescent fluid.
///
/// **v_s = 2 r² (ρ_p − ρ_f) g / (9 μ)**
///
/// Positive value means downward (gravity direction).
#[allow(dead_code)]
pub fn settling_velocity_stokes(
    radius: f64,
    rho_particle: f64,
    rho_fluid: f64,
    viscosity: f64,
) -> f64 {
    let g = 9.81_f64;
    2.0 * radius * radius * (rho_particle - rho_fluid) * g / (9.0 * viscosity)
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn make_particle(pos: [f64; 2], vel: [f64; 2], r: f64) -> CoupledParticle {
        CoupledParticle::new(pos, vel, 1.0, r, 0)
    }

    // ── DragForce tests ──────────────────────────────────────────────────────

    #[test]
    fn test_stokes_drag_zero_slip() {
        let p = make_particle([5.0, 5.0], [0.1, 0.2], 0.5);
        let f = FluidCell::new(1.0, 0.1, 0.2);
        let drag = DragForce::stokes_drag(&p, &f, 1.0 / 6.0);
        assert!(
            drag[0].abs() < 1e-14 && drag[1].abs() < 1e-14,
            "Zero slip → zero drag: {:?}",
            drag
        );
    }

    #[test]
    fn test_stokes_drag_direction() {
        // Fluid faster than particle in x → force should be in +x.
        let p = make_particle([5.0, 5.0], [0.0, 0.0], 0.5);
        let f = FluidCell::new(1.0, 0.1, 0.0);
        let drag = DragForce::stokes_drag(&p, &f, 1.0 / 6.0);
        assert!(
            drag[0] > 0.0,
            "Drag should accelerate particle: {}",
            drag[0]
        );
        assert!(drag[1].abs() < 1e-14, "No y-drag expected: {}", drag[1]);
    }

    #[test]
    fn test_stokes_drag_magnitude() {
        // F = 6 π μ r Δu
        let mu = 1.0 / 6.0;
        let r = 1.0;
        let delta_u = 1.0;
        let p = make_particle([0.0, 0.0], [0.0, 0.0], r);
        let f = FluidCell::new(1.0, delta_u, 0.0);
        let drag = DragForce::stokes_drag(&p, &f, mu);
        let expected = 6.0 * PI * mu * r * delta_u;
        assert!(
            (drag[0] - expected).abs() < 1e-12,
            "Stokes drag magnitude: got {}, expected {}",
            drag[0],
            expected
        );
    }

    #[test]
    fn test_reynolds_number_zero_slip() {
        let p = make_particle([0.0, 0.0], [0.5, 0.3], 0.5);
        let f = FluidCell::new(1.0, 0.5, 0.3);
        let re = DragForce::reynolds_number(&p, &f, 1.0 / 6.0);
        assert!(re < 1e-14, "Re should be 0 for zero slip: {re}");
    }

    #[test]
    fn test_reynolds_number_positive() {
        let p = make_particle([0.0, 0.0], [0.0, 0.0], 0.5);
        let f = FluidCell::new(1.0, 1.0, 0.0);
        let re = DragForce::reynolds_number(&p, &f, 1.0 / 6.0);
        // Re = 2*0.5*1.0 / (1/6) = 1.0 * 6 = 6
        assert!((re - 6.0).abs() < 1e-12, "Re = {re}, expected 6.0");
    }

    #[test]
    fn test_schiller_naumann_low_re_approaches_stokes() {
        // At very low Re, SN should closely match Stokes drag.
        let mu = 1.0 / 6.0;
        let r = 0.01; // tiny radius → tiny Re
        let p = make_particle([0.0, 0.0], [0.0, 0.0], r);
        let f = FluidCell::new(1.0, 0.001, 0.0);
        let stokes = DragForce::stokes_drag(&p, &f, mu);
        let sn = DragForce::schiller_naumann_drag(&p, &f, mu, 1.0);
        // They won't be identical (SN adds 0.15*Re^0.687 correction) but
        // both should have the same sign and similar magnitude.
        assert!(
            stokes[0] > 0.0 && sn[0] > 0.0,
            "Both drags must be positive"
        );
    }

    #[test]
    fn test_schiller_naumann_high_re() {
        // Re > 1000: Cd = 0.44.
        let mu = 1e-6; // very low viscosity → very high Re
        let r = 0.5;
        let p = make_particle([0.0, 0.0], [0.0, 0.0], r);
        let f = FluidCell::new(1.0, 10.0, 0.0);
        let sn = DragForce::schiller_naumann_drag(&p, &f, mu, 1.0);
        // Should be non-zero and in the drag direction.
        assert!(sn[0] > 0.0, "High-Re SN drag should be positive: {}", sn[0]);
    }

    #[test]
    fn test_schiller_naumann_zero_slip() {
        let p = make_particle([0.0, 0.0], [0.5, 0.0], 0.5);
        let f = FluidCell::new(1.0, 0.5, 0.0);
        let sn = DragForce::schiller_naumann_drag(&p, &f, 1.0 / 6.0, 1.0);
        assert!(
            sn[0].abs() < 1e-14 && sn[1].abs() < 1e-14,
            "Zero slip → zero SN drag: {:?}",
            sn
        );
    }

    // ── TwoWayCoupling tests ─────────────────────────────────────────────────

    #[test]
    fn test_interpolate_exact_cell_center() {
        let coupling = TwoWayCoupling::new(4, 4, 1.0);
        let mut fluid = vec![(0.0_f64, 0.0_f64); 16];
        // Set cell (2, 1) to a known velocity.
        fluid[6] = (0.7, -0.3);
        let vel = coupling.interpolate_fluid_velocity(&fluid, [2.0, 1.0]);
        assert!((vel[0] - 0.7).abs() < 1e-12, "ux mismatch: {}", vel[0]);
        assert!((vel[1] - (-0.3)).abs() < 1e-12, "uy mismatch: {}", vel[1]);
    }

    #[test]
    fn test_interpolate_uniform_field() {
        let coupling = TwoWayCoupling::new(5, 5, 1.0);
        let fluid = vec![(0.5_f64, -0.2_f64); 25];
        // Any position in a uniform field should return the same velocity.
        let vel = coupling.interpolate_fluid_velocity(&fluid, [2.3, 1.7]);
        assert!((vel[0] - 0.5).abs() < 1e-12, "ux = {}", vel[0]);
        assert!((vel[1] - (-0.2)).abs() < 1e-12, "uy = {}", vel[1]);
    }

    #[test]
    fn test_spread_particle_force_single() {
        let coupling = TwoWayCoupling::new(4, 4, 1.0);
        let mut p = make_particle([2.0, 1.0], [0.0, 0.0], 0.5);
        p.force = [3.0, -1.0];
        let forces = coupling.spread_particle_force(&[p], 4, 4);
        let k = 6;
        // The reaction is −force on the fluid cell.
        assert!(
            (forces[k][0] - (-3.0)).abs() < 1e-12,
            "Grid fx: {}",
            forces[k][0]
        );
        assert!(
            (forces[k][1] - 1.0).abs() < 1e-12,
            "Grid fy: {}",
            forces[k][1]
        );
    }

    #[test]
    fn test_update_particle_velocity() {
        let mut p = make_particle([0.0, 0.0], [0.0, 0.0], 0.5);
        let coupling = TwoWayCoupling::new(10, 10, 1.0);
        let drag = [2.0, 0.0];
        let gravity = [0.0, -1.0];
        coupling.update_particle(&mut p, drag, gravity, 1.0);
        // v = (2, -1) * 1 / mass=1.0
        assert!(
            (p.velocity[0] - 2.0).abs() < 1e-12,
            "vx = {}",
            p.velocity[0]
        );
        assert!(
            (p.velocity[1] - (-1.0)).abs() < 1e-12,
            "vy = {}",
            p.velocity[1]
        );
    }

    #[test]
    fn test_update_particle_position() {
        let mut p = make_particle([1.0, 2.0], [0.0, 0.0], 0.5);
        let coupling = TwoWayCoupling::new(10, 10, 1.0);
        // No acceleration, but give it an initial velocity via update.
        coupling.update_particle(&mut p, [0.0, 0.0], [0.0, 0.0], 1.0);
        // Position should not change (v was zero).
        assert!((p.position[0] - 1.0).abs() < 1e-12);
        assert!((p.position[1] - 2.0).abs() < 1e-12);
    }

    // ── VoidFraction tests ───────────────────────────────────────────────────

    #[test]
    fn test_void_fraction_no_particles() {
        let vf = VoidFraction::compute_void_fraction(&[], 4, 4, 1.0);
        for &v in &vf {
            assert!((v - 1.0).abs() < 1e-14, "Empty domain: void fraction = {v}");
        }
    }

    #[test]
    fn test_void_fraction_reduces_with_particle() {
        let p = make_particle([2.0, 2.0], [0.0, 0.0], 0.3);
        let vf = VoidFraction::compute_void_fraction(&[p], 5, 5, 1.0);
        // Cell (2,2): void fraction < 1.
        let k = 2 * 5 + 2;
        assert!(
            vf[k] < 1.0,
            "Void fraction should be < 1 where particle sits: {}",
            vf[k]
        );
    }

    #[test]
    fn test_void_fraction_clamped() {
        // Very large particle → solid fraction > 1 without clamping.
        let p = make_particle([2.0, 2.0], [0.0, 0.0], 5.0);
        let vf = VoidFraction::compute_void_fraction(&[p], 5, 5, 1.0);
        for &v in &vf {
            assert!(v >= 0.0 && v <= 1.0, "Void fraction out of range: {v}");
        }
    }

    #[test]
    fn test_corrected_drag_unit_fraction() {
        // ε = 1 → correction factor = 1.
        let drag = [1.0, -2.0];
        let corrected = VoidFraction::corrected_drag(drag, 1.0);
        assert!((corrected[0] - 1.0).abs() < 1e-12);
        assert!((corrected[1] - (-2.0)).abs() < 1e-12);
    }

    #[test]
    fn test_corrected_drag_increases_for_lower_eps() {
        // Lower void fraction → larger correction.
        let drag = [1.0, 0.0];
        let c_high = VoidFraction::corrected_drag(drag, 0.9);
        let c_low = VoidFraction::corrected_drag(drag, 0.5);
        assert!(
            c_low[0] > c_high[0],
            "Lower ε should give larger corrected drag: {} vs {}",
            c_low[0],
            c_high[0]
        );
    }

    // ── ParticleCollision tests ──────────────────────────────────────────────

    #[test]
    fn test_detect_wall_inside() {
        // Particle well inside the domain: no collision.
        let p = make_particle([5.0, 5.0], [0.0, 0.0], 0.5);
        assert!(!ParticleCollision::detect_particle_wall(&p, 10, 10, 1.0));
    }

    #[test]
    fn test_detect_wall_outside() {
        // Particle touching left wall.
        let p = make_particle([0.3, 5.0], [0.0, 0.0], 0.5);
        assert!(ParticleCollision::detect_particle_wall(&p, 10, 10, 1.0));
    }

    #[test]
    fn test_reflect_from_left_wall() {
        let mut p = make_particle([0.1, 5.0], [-1.0, 0.0], 0.5);
        ParticleCollision::reflect_from_wall(&mut p, 10, 10, 1.0, 1.0);
        assert!(
            p.velocity[0] > 0.0,
            "Velocity should reverse: {}",
            p.velocity[0]
        );
        assert!(
            p.position[0] >= 0.5,
            "Position should be corrected: {}",
            p.position[0]
        );
    }

    #[test]
    fn test_reflect_from_bottom_wall() {
        let mut p = make_particle([5.0, 0.2], [0.0, -1.0], 0.5);
        ParticleCollision::reflect_from_wall(&mut p, 10, 10, 1.0, 0.8);
        assert!(p.velocity[1] > 0.0, "vy should reverse: {}", p.velocity[1]);
    }

    #[test]
    fn test_reflect_restitution() {
        let mut p = make_particle([0.1, 5.0], [-2.0, 0.0], 0.5);
        ParticleCollision::reflect_from_wall(&mut p, 10, 10, 1.0, 0.5);
        // Speed after = 0.5 * |original speed|
        assert!(
            (p.velocity[0] - 1.0).abs() < 1e-12,
            "Restitution not applied: {}",
            p.velocity[0]
        );
    }

    #[test]
    fn test_particle_particle_no_overlap() {
        let mut p1 = make_particle([0.0, 0.0], [1.0, 0.0], 0.5);
        let mut p2 = make_particle([5.0, 0.0], [-1.0, 0.0], 0.5);
        let v1_before = p1.velocity[0];
        let v2_before = p2.velocity[0];
        ParticleCollision::particle_particle_collision(&mut p1, &mut p2, 1.0);
        // No overlap → no velocity change.
        assert!((p1.velocity[0] - v1_before).abs() < 1e-14);
        assert!((p2.velocity[0] - v2_before).abs() < 1e-14);
    }

    #[test]
    fn test_particle_particle_elastic_collision() {
        // Equal masses head-on → velocities exchange.
        let mut p1 = make_particle([0.0, 0.0], [1.0, 0.0], 0.4);
        let mut p2 = make_particle([0.5, 0.0], [-1.0, 0.0], 0.4);
        ParticleCollision::particle_particle_collision(&mut p1, &mut p2, 1.0);
        // After elastic collision with equal masses: velocities should have swapped signs.
        assert!(
            p1.velocity[0] < 0.0,
            "p1 should reverse: {}",
            p1.velocity[0]
        );
        assert!(
            p2.velocity[0] > 0.0,
            "p2 should reverse: {}",
            p2.velocity[0]
        );
    }

    // ── settling_velocity_stokes ─────────────────────────────────────────────

    #[test]
    fn test_settling_velocity_stokes_positive_for_dense_particle() {
        // Denser particle than fluid → positive settling velocity (downward).
        let vs = settling_velocity_stokes(0.5, 2.0, 1.0, 1.0 / 6.0);
        assert!(vs > 0.0, "Settling velocity should be positive: {vs}");
    }

    #[test]
    fn test_settling_velocity_stokes_negative_for_buoyant() {
        // Lighter particle than fluid → rises (negative velocity).
        let vs = settling_velocity_stokes(0.5, 0.5, 1.0, 1.0 / 6.0);
        assert!(vs < 0.0, "Buoyant particle: {vs}");
    }

    #[test]
    fn test_settling_velocity_stokes_formula() {
        let r = 1.0;
        let rho_p = 2.0;
        let rho_f = 1.0;
        let mu = 1.0;
        let vs = settling_velocity_stokes(r, rho_p, rho_f, mu);
        let expected = 2.0 * r * r * (rho_p - rho_f) * 9.81 / (9.0 * mu);
        assert!(
            (vs - expected).abs() < 1e-12,
            "Formula mismatch: {vs} vs {expected}"
        );
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 3-D coupled particle with rotation
// ─────────────────────────────────────────────────────────────────────────────

/// A 3-D coupled particle that tracks position, velocity, angular velocity,
/// and all associated forces/torques.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct CoupledParticle3D {
    /// Position \[x, y, z\] in lattice units.
    pub position: [f64; 3],
    /// Translational velocity \[ux, uy, uz\].
    pub velocity: [f64; 3],
    /// Angular velocity \[wx, wy, wz\].
    pub angular_velocity: [f64; 3],
    /// Net force \[fx, fy, fz\].
    pub force: [f64; 3],
    /// Net torque \[tx, ty, tz\].
    pub torque: [f64; 3],
    /// Particle mass.
    pub mass: f64,
    /// Particle radius.
    pub radius: f64,
    /// Unique identifier.
    pub id: usize,
}

#[allow(dead_code)]
impl CoupledParticle3D {
    /// Create a new 3-D coupled particle.
    pub fn new(position: [f64; 3], velocity: [f64; 3], mass: f64, radius: f64, id: usize) -> Self {
        Self {
            position,
            velocity,
            angular_velocity: [0.0; 3],
            force: [0.0; 3],
            torque: [0.0; 3],
            mass,
            radius,
            id,
        }
    }

    /// Moment of inertia for a solid sphere: I = 2/5 m r².
    pub fn moment_of_inertia(&self) -> f64 {
        0.4 * self.mass * self.radius * self.radius
    }

    /// Translational kinetic energy.
    pub fn translational_ke(&self) -> f64 {
        let v2 = self.velocity[0].powi(2) + self.velocity[1].powi(2) + self.velocity[2].powi(2);
        0.5 * self.mass * v2
    }

    /// Rotational kinetic energy: ½ I ω².
    pub fn rotational_ke(&self) -> f64 {
        let w2 = self.angular_velocity[0].powi(2)
            + self.angular_velocity[1].powi(2)
            + self.angular_velocity[2].powi(2);
        0.5 * self.moment_of_inertia() * w2
    }

    /// Total kinetic energy.
    pub fn total_ke(&self) -> f64 {
        self.translational_ke() + self.rotational_ke()
    }

    /// Update translational velocity and position using semi-implicit Euler.
    ///
    /// v^{n+1} = v^n + (F/m) * dt
    /// x^{n+1} = x^n + v^{n+1} * dt
    pub fn integrate_translation(&mut self, dt: f64) {
        let inv_m = 1.0 / self.mass.max(f64::EPSILON);
        for k in 0..3 {
            self.velocity[k] += self.force[k] * inv_m * dt;
            self.position[k] += self.velocity[k] * dt;
        }
    }

    /// Update angular velocity from torque: ω^{n+1} = ω^n + (T/I) * dt.
    pub fn integrate_rotation(&mut self, dt: f64) {
        let inv_i = 1.0 / self.moment_of_inertia().max(f64::EPSILON);
        for k in 0..3 {
            self.angular_velocity[k] += self.torque[k] * inv_i * dt;
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Particle lift force models
// ─────────────────────────────────────────────────────────────────────────────

/// Particle lift force utilities.
#[allow(dead_code)]
pub struct LiftForce;

#[allow(dead_code)]
impl LiftForce {
    /// Saffman lift force on a particle in a shear flow.
    ///
    /// F_L = C_L * rho_f * (2r)² * |shear_rate| * slip_speed
    ///
    /// The direction is perpendicular to the slip velocity, in the
    /// direction of increasing velocity gradient.
    ///
    /// # Arguments
    /// * `particle` – the coupled particle (2-D)
    /// * `fluid_ux_gradient` – du_x / dy (shear rate)
    /// * `rho_f` – fluid density
    /// * `c_l` – lift coefficient (≈ 1.615 for Saffman)
    pub fn saffman_lift(
        particle: &CoupledParticle,
        fluid: &FluidCell,
        fluid_ux_gradient: f64,
        rho_f: f64,
        c_l: f64,
    ) -> [f64; 2] {
        let dux = fluid.ux - particle.velocity[0];
        let duy = fluid.uy - particle.velocity[1];
        let slip_speed = (dux * dux + duy * duy).sqrt();
        let diameter = 2.0 * particle.radius;
        // Lift magnitude: F_L = c_l * rho_f * d² * |G| * |u_slip|
        let lift_mag = c_l * rho_f * diameter * diameter * fluid_ux_gradient.abs() * slip_speed;
        // Direction: perpendicular to slip, sign from sign of shear × slip_x
        let lift_dir = if fluid_ux_gradient >= 0.0 { 1.0 } else { -1.0 };
        // Perpendicular to slip = (-duy, dux) / slip_speed
        if slip_speed < f64::EPSILON {
            return [0.0, 0.0];
        }
        [
            lift_mag * lift_dir * (-duy / slip_speed),
            lift_mag * lift_dir * (dux / slip_speed),
        ]
    }

    /// Magnus lift force due to particle rotation in a fluid.
    ///
    /// F_Magnus = (1/2) * rho_f * V_p * (Ω × u_slip)
    ///
    /// In 2-D: Ω is the z-component of angular velocity.
    /// F = (1/2) * rho_f * V_p * Ω_z * \[-duy, dux\]
    ///
    /// # Arguments
    /// * `omega_z` – angular velocity z-component (rad/s)
    /// * `rho_f` – fluid density
    pub fn magnus_lift_2d(
        particle: &CoupledParticle,
        fluid: &FluidCell,
        omega_z: f64,
        rho_f: f64,
    ) -> [f64; 2] {
        let volume = (4.0 / 3.0) * PI * particle.radius.powi(3);
        let dux = fluid.ux - particle.velocity[0];
        let duy = fluid.uy - particle.velocity[1];
        let factor = 0.5 * rho_f * volume * omega_z;
        // Ω × u_slip in 2-D: (Ω_z e_z) × (dux, duy, 0) = (-Ω_z duy, Ω_z dux, 0)
        [factor * (-duy), factor * dux]
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Particle torque from fluid
// ─────────────────────────────────────────────────────────────────────────────

/// Torque on a particle due to fluid velocity gradients.
#[allow(dead_code)]
pub struct FluidTorque;

#[allow(dead_code)]
impl FluidTorque {
    /// Hydrodynamic torque on a sphere in a linear shear flow.
    ///
    /// T = 8π μ r³ (Ω_fluid - Ω_particle)
    ///
    /// where Ω_fluid = (1/2) * (duy/dx - dux/dy) is the fluid rotation rate.
    ///
    /// Returns the z-component of torque in 2-D.
    pub fn shear_torque_2d(
        particle: &CoupledParticle,
        fluid_rotation_rate: f64,
        particle_omega_z: f64,
        viscosity: f64,
    ) -> f64 {
        8.0 * PI * viscosity * particle.radius.powi(3) * (fluid_rotation_rate - particle_omega_z)
    }

    /// Pitching torque estimate from asymmetric drag due to rotation.
    ///
    /// T ≈ -mu_r * r * F_n * sign(omega)
    pub fn rolling_drag_torque(omega: f64, fn_mag: f64, mu_r: f64, radius: f64) -> f64 {
        if omega.abs() < f64::EPSILON {
            return 0.0;
        }
        -mu_r * radius * fn_mag * omega.signum()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Particle rotation update
// ─────────────────────────────────────────────────────────────────────────────

/// Integrate the angular velocity of a set of 3-D particles one timestep.
///
/// Uses the semi-implicit Euler: ω^{n+1} = ω^n + (T/I) * dt
#[allow(dead_code)]
pub fn integrate_rotations(particles: &mut [CoupledParticle3D], dt: f64) {
    for p in particles.iter_mut() {
        p.integrate_rotation(dt);
    }
}

/// Integrate the translations of a set of 3-D particles one timestep.
#[allow(dead_code)]
pub fn integrate_translations(particles: &mut [CoupledParticle3D], dt: f64) {
    for p in particles.iter_mut() {
        p.integrate_translation(dt);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Particle settling: analytical and simulated
// ─────────────────────────────────────────────────────────────────────────────

/// Terminal settling velocity including buoyancy.
///
/// v_terminal = 2 r² g (rho_p - rho_f) / (9 mu)
///
/// Identical to [`settling_velocity_stokes`] but explicitly returns the
/// terminal velocity signed positive downward.
#[allow(dead_code)]
pub fn terminal_settling_velocity(
    radius: f64,
    rho_particle: f64,
    rho_fluid: f64,
    viscosity: f64,
    g: f64,
) -> f64 {
    2.0 * radius * radius * (rho_particle - rho_fluid) * g / (9.0 * viscosity)
}

/// Simulate Stokes settling of a single particle in a quiescent fluid.
///
/// Returns a vector of (time, y_position, y_velocity) tuples.
///
/// The equation of motion is:
/// m * dv/dt = F_drag + F_buoyancy + F_gravity
///            = -6π μ r v + (rho_p - rho_f) V g_mag downward
#[allow(dead_code)]
#[allow(clippy::too_many_arguments)]
pub fn simulate_settling(
    radius: f64,
    rho_particle: f64,
    rho_fluid: f64,
    viscosity: f64,
    g: f64,
    total_time: f64,
    dt: f64,
) -> Vec<(f64, f64, f64)> {
    let mass = (4.0 / 3.0) * PI * radius.powi(3) * rho_particle;
    let drag_coeff = 6.0 * PI * viscosity * radius;
    let volume = (4.0 / 3.0) * PI * radius.powi(3);
    let buoyancy_net = (rho_particle - rho_fluid) * volume * g;

    let mut t = 0.0;
    let mut y = 0.0;
    let mut vy = 0.0;
    let mut history = Vec::new();

    while t <= total_time {
        history.push((t, y, vy));
        let drag = -drag_coeff * vy;
        let net_force = drag + buoyancy_net;
        let ay = net_force / mass;
        vy += ay * dt;
        y += vy * dt;
        t += dt;
    }
    history
}

// ─────────────────────────────────────────────────────────────────────────────
// Two-way coupling: two-phase flow corrections
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the fluid momentum source density from particle reaction forces.
///
/// For each cell (ix, iy), the momentum source is:
///   S_ij = - Σ_{particles in cell} F_p / (cell_volume)
///
/// Returns a flat `nx × ny` array of `[sx, sy]`.
#[allow(dead_code)]
pub fn fluid_momentum_source(
    particles: &[CoupledParticle],
    nx: usize,
    ny: usize,
    dx: f64,
) -> Vec<[f64; 2]> {
    let cell_area = dx * dx;
    let mut source = vec![[0.0_f64; 2]; nx * ny];
    for p in particles {
        let ix = ((p.position[0] / dx).round() as usize).min(nx - 1);
        let iy = ((p.position[1] / dx).round() as usize).min(ny - 1);
        let k = iy * nx + ix;
        source[k][0] -= p.force[0] / cell_area;
        source[k][1] -= p.force[1] / cell_area;
    }
    source
}

/// Compute local particle volume fraction (solid fraction) on the grid.
///
/// Returns the particle volume fraction ε_p = V_particle / V_cell for each cell.
#[allow(dead_code)]
pub fn particle_volume_fraction(
    particles: &[CoupledParticle],
    nx: usize,
    ny: usize,
    dx: f64,
) -> Vec<f64> {
    let cell_volume = dx * dx; // 2D cell area treated as volume
    let mut pf = vec![0.0_f64; nx * ny];
    for p in particles {
        let ix = ((p.position[0] / dx).round() as usize).min(nx - 1);
        let iy = ((p.position[1] / dx).round() as usize).min(ny - 1);
        let k = iy * nx + ix;
        let vp = PI * p.radius * p.radius; // 2-D projected area
        pf[k] += vp / cell_volume;
    }
    pf.iter().map(|&v| v.min(1.0)).collect()
}

// ─────────────────────────────────────────────────────────────────────────────
// Multi-particle system manager
// ─────────────────────────────────────────────────────────────────────────────

/// Manager for a collection of 3-D coupled particles with fluid interaction.
#[allow(dead_code)]
pub struct ParticleSystem3D {
    /// All tracked particles.
    pub particles: Vec<CoupledParticle3D>,
    /// Gravity vector.
    pub gravity: [f64; 3],
    /// Fluid dynamic viscosity.
    pub viscosity: f64,
    /// Fluid density.
    pub rho_fluid: f64,
}

#[allow(dead_code)]
impl ParticleSystem3D {
    /// Create an empty system.
    pub fn new(gravity: [f64; 3], viscosity: f64, rho_fluid: f64) -> Self {
        Self {
            particles: Vec::new(),
            gravity,
            viscosity,
            rho_fluid,
        }
    }

    /// Add a particle.
    pub fn add_particle(&mut self, p: CoupledParticle3D) {
        self.particles.push(p);
    }

    /// Number of particles.
    pub fn len(&self) -> usize {
        self.particles.len()
    }

    /// True if no particles.
    pub fn is_empty(&self) -> bool {
        self.particles.is_empty()
    }

    /// Apply gravity body force to all particles.
    pub fn apply_gravity(&mut self) {
        for p in &mut self.particles {
            for k in 0..3 {
                p.force[k] = p.mass * self.gravity[k];
            }
        }
    }

    /// Apply Stokes drag to each particle given a uniform background fluid velocity.
    pub fn apply_stokes_drag(&mut self, fluid_vel: [f64; 3]) {
        let mu = self.viscosity;
        for p in &mut self.particles {
            let factor = 6.0 * PI * mu * p.radius;
            for k in 0..3 {
                p.force[k] += factor * (fluid_vel[k] - p.velocity[k]);
            }
        }
    }

    /// Advance all particles one step.
    pub fn step(&mut self, dt: f64) {
        for p in &mut self.particles {
            p.integrate_translation(dt);
            p.integrate_rotation(dt);
        }
    }

    /// Total translational kinetic energy.
    pub fn total_translational_ke(&self) -> f64 {
        self.particles.iter().map(|p| p.translational_ke()).sum()
    }

    /// Total rotational kinetic energy.
    pub fn total_rotational_ke(&self) -> f64 {
        self.particles.iter().map(|p| p.rotational_ke()).sum()
    }

    /// Centre of mass position.
    pub fn centre_of_mass(&self) -> [f64; 3] {
        let total_mass: f64 = self.particles.iter().map(|p| p.mass).sum();
        if total_mass < f64::EPSILON {
            return [0.0; 3];
        }
        let mut com = [0.0_f64; 3];
        for p in &self.particles {
            for k in 0..3 {
                com[k] += p.mass * p.position[k];
            }
        }
        let inv = 1.0 / total_mass;
        [com[0] * inv, com[1] * inv, com[2] * inv]
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests for new functionality
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests_extended {
    use super::*;

    fn make_3d_particle() -> CoupledParticle3D {
        CoupledParticle3D::new([0.0; 3], [0.0; 3], 1.0, 0.5, 0)
    }

    // ── CoupledParticle3D tests ───────────────────────────────────────────

    #[test]
    fn particle3d_moment_of_inertia() {
        let p = make_3d_particle();
        let expected = 0.4 * 1.0 * 0.25;
        assert!((p.moment_of_inertia() - expected).abs() < 1e-14);
    }

    #[test]
    fn particle3d_ke_zero_at_rest() {
        let p = make_3d_particle();
        assert!((p.total_ke()).abs() < 1e-14);
    }

    #[test]
    fn particle3d_translational_ke_moving() {
        let mut p = make_3d_particle();
        p.velocity = [1.0, 0.0, 0.0];
        assert!((p.translational_ke() - 0.5).abs() < 1e-14);
    }

    #[test]
    fn particle3d_rotational_ke_spinning() {
        let mut p = make_3d_particle();
        p.angular_velocity = [1.0, 0.0, 0.0];
        let i = p.moment_of_inertia();
        assert!((p.rotational_ke() - 0.5 * i).abs() < 1e-14);
    }

    #[test]
    fn particle3d_integrate_translation_gravity() {
        let mut p = make_3d_particle();
        p.force = [0.0, -9.81 * p.mass, 0.0];
        p.integrate_translation(0.1);
        // vy = -9.81 * 0.1 = -0.981
        assert!(
            (p.velocity[1] - (-0.981)).abs() < 1e-12,
            "vy = {}",
            p.velocity[1]
        );
        // y = vy * dt = -0.981 * 0.1 = -0.0981
        assert!(
            (p.position[1] - (-0.0981)).abs() < 1e-12,
            "y = {}",
            p.position[1]
        );
    }

    #[test]
    fn particle3d_integrate_rotation_torque() {
        let mut p = make_3d_particle();
        let i = p.moment_of_inertia();
        p.torque = [i, 0.0, 0.0]; // angular acceleration = 1 rad/s²
        p.integrate_rotation(0.1);
        assert!(
            (p.angular_velocity[0] - 0.1).abs() < 1e-12,
            "omega_x = {}",
            p.angular_velocity[0]
        );
    }

    // ── LiftForce tests ───────────────────────────────────────────────────

    #[test]
    fn saffman_lift_zero_slip() {
        let p = CoupledParticle::new([0.0, 0.0], [0.5, 0.3], 1.0, 0.1, 0);
        let f = FluidCell::new(1.0, 0.5, 0.3);
        let lift = LiftForce::saffman_lift(&p, &f, 0.1, 1.0, 1.615);
        assert!(
            lift[0].abs() < 1e-14 && lift[1].abs() < 1e-14,
            "Zero slip → zero lift: {:?}",
            lift
        );
    }

    #[test]
    fn saffman_lift_finite_for_nonzero_shear() {
        let p = CoupledParticle::new([0.0, 0.0], [0.0, 0.0], 1.0, 0.1, 0);
        let f = FluidCell::new(1.0, 0.5, 0.0);
        let lift = LiftForce::saffman_lift(&p, &f, 1.0, 1.0, 1.615);
        for &li in &lift {
            assert!(li.is_finite(), "Saffman lift must be finite: {li}");
        }
    }

    #[test]
    fn magnus_lift_zero_angular_velocity() {
        let p = CoupledParticle::new([0.0, 0.0], [0.0, 0.0], 1.0, 0.1, 0);
        let f = FluidCell::new(1.0, 1.0, 0.0);
        let lift = LiftForce::magnus_lift_2d(&p, &f, 0.0, 1.0);
        assert!(
            lift[0].abs() < 1e-14 && lift[1].abs() < 1e-14,
            "Zero omega → zero Magnus lift: {:?}",
            lift
        );
    }

    #[test]
    fn magnus_lift_nonzero_for_spinning_particle() {
        let p = CoupledParticle::new([0.0, 0.0], [0.0, 0.0], 1.0, 0.1, 0);
        let f = FluidCell::new(1.0, 1.0, 0.0);
        let lift = LiftForce::magnus_lift_2d(&p, &f, 10.0, 1.0);
        // With ux_slip=1.0 and omega_z=10.0, there should be a nonzero lift
        assert!(
            lift[0].abs() + lift[1].abs() > 0.0,
            "Magnus lift must be nonzero: {:?}",
            lift
        );
    }

    // ── FluidTorque tests ─────────────────────────────────────────────────

    #[test]
    fn shear_torque_no_relative_rotation() {
        let p = CoupledParticle::new([0.0, 0.0], [0.0, 0.0], 1.0, 0.5, 0);
        let t = FluidTorque::shear_torque_2d(&p, 1.0, 1.0, 0.1);
        assert!(t.abs() < 1e-14, "Equal rotation rates → zero torque: {t}");
    }

    #[test]
    fn shear_torque_positive_for_faster_fluid() {
        let p = CoupledParticle::new([0.0, 0.0], [0.0, 0.0], 1.0, 0.5, 0);
        // Fluid rotating faster than particle → positive torque
        let t = FluidTorque::shear_torque_2d(&p, 5.0, 0.0, 0.1);
        assert!(t > 0.0, "Faster fluid rotation → positive torque: {t}");
    }

    #[test]
    fn rolling_drag_torque_zero_for_no_rotation() {
        let t = FluidTorque::rolling_drag_torque(0.0, 100.0, 0.1, 0.5);
        assert_eq!(t, 0.0);
    }

    #[test]
    fn rolling_drag_torque_opposes_rotation() {
        let t = FluidTorque::rolling_drag_torque(1.0, 100.0, 0.1, 0.5);
        assert!(t < 0.0, "Rolling drag should oppose positive rotation: {t}");
        let t_neg = FluidTorque::rolling_drag_torque(-1.0, 100.0, 0.1, 0.5);
        assert!(
            t_neg > 0.0,
            "Rolling drag should oppose negative rotation: {t_neg}"
        );
    }

    // ── Terminal settling velocity ─────────────────────────────────────────

    #[test]
    fn terminal_settling_matches_stokes() {
        let r = 0.5;
        let rho_p = 2.0;
        let rho_f = 1.0;
        let mu = 1.0 / 6.0;
        let g = 9.81;
        let v_t = terminal_settling_velocity(r, rho_p, rho_f, mu, g);
        let v_s = settling_velocity_stokes(r, rho_p, rho_f, mu);
        // Both should give the same result
        assert!((v_t - v_s).abs() < 1e-10, "terminal: {v_t}, stokes: {v_s}");
    }

    // ── simulate_settling tests ───────────────────────────────────────────

    #[test]
    fn simulate_settling_dense_particle_falls() {
        let traj = simulate_settling(0.1, 2000.0, 1000.0, 1e-3, 9.81, 0.1, 0.001);
        assert!(!traj.is_empty());
        let last = traj.last().unwrap();
        // Dense particle should have moved downward (positive y since g>0)
        assert!(last.1 > 0.0, "dense particle should fall: y = {}", last.1);
    }

    #[test]
    fn simulate_settling_buoyant_particle_rises() {
        let traj = simulate_settling(0.1, 100.0, 1000.0, 1e-3, 9.81, 0.1, 0.001);
        let last = traj.last().unwrap();
        // Buoyant particle: rho_p < rho_f → net upward → negative y
        assert!(last.1 < 0.0, "buoyant particle should rise: y = {}", last.1);
    }

    #[test]
    fn simulate_settling_approaches_terminal_velocity() {
        // After sufficient time, velocity should approach terminal value
        let r = 0.001;
        let rho_p = 2000.0;
        let rho_f = 1000.0;
        let mu = 1e-3;
        let g = 9.81;
        let v_terminal = terminal_settling_velocity(r, rho_p, rho_f, mu, g);
        let traj = simulate_settling(r, rho_p, rho_f, mu, g, 10.0, 0.01);
        let last_vy = traj.last().unwrap().2;
        let rel_err = (last_vy - v_terminal).abs() / v_terminal.abs().max(1e-30);
        assert!(
            rel_err < 0.01,
            "velocity should reach ~terminal: {last_vy} vs {v_terminal}"
        );
    }

    // ── fluid_momentum_source tests ──────────────────────────────────────

    #[test]
    fn fluid_momentum_source_empty_particles() {
        let src = fluid_momentum_source(&[], 4, 4, 1.0);
        for &s in &src {
            assert_eq!(s, [0.0, 0.0]);
        }
    }

    #[test]
    fn fluid_momentum_source_single_particle() {
        let mut p = CoupledParticle::new([2.0, 2.0], [0.0, 0.0], 1.0, 0.3, 0);
        p.force = [1.0, -2.0];
        let src = fluid_momentum_source(&[p], 5, 5, 1.0);
        let k = 2 * 5 + 2;
        // Reaction: fluid receives -force / cell_area
        assert!((src[k][0] - (-1.0)).abs() < 1e-12, "src_x = {}", src[k][0]);
        assert!((src[k][1] - 2.0).abs() < 1e-12, "src_y = {}", src[k][1]);
    }

    // ── particle_volume_fraction tests ───────────────────────────────────

    #[test]
    fn particle_volume_fraction_empty() {
        let pf = particle_volume_fraction(&[], 4, 4, 1.0);
        for &v in &pf {
            assert_eq!(v, 0.0);
        }
    }

    #[test]
    fn particle_volume_fraction_nonzero_where_particle_sits() {
        let p = CoupledParticle::new([2.0, 2.0], [0.0, 0.0], 1.0, 0.3, 0);
        let pf = particle_volume_fraction(&[p], 5, 5, 1.0);
        let k = 2 * 5 + 2;
        assert!(
            pf[k] > 0.0,
            "volume fraction should be > 0 at particle location"
        );
        assert!(pf[k] <= 1.0, "volume fraction should be <= 1");
    }

    // ── ParticleSystem3D tests ────────────────────────────────────────────

    #[test]
    fn particle_system_add_and_count() {
        let mut sys = ParticleSystem3D::new([0.0, -9.81, 0.0], 1e-3, 1000.0);
        sys.add_particle(make_3d_particle());
        sys.add_particle(make_3d_particle());
        assert_eq!(sys.len(), 2);
        assert!(!sys.is_empty());
    }

    #[test]
    fn particle_system_gravity_increases_ke() {
        let mut sys = ParticleSystem3D::new([0.0, -9.81, 0.0], 0.0, 1000.0);
        sys.add_particle(make_3d_particle());
        let ke_before = sys.total_translational_ke();
        sys.apply_gravity();
        sys.step(0.1);
        let ke_after = sys.total_translational_ke();
        assert!(ke_after > ke_before, "gravity should increase KE");
    }

    #[test]
    fn particle_system_stokes_drag_decelerates() {
        let mut sys = ParticleSystem3D::new([0.0; 3], 1.0, 1000.0);
        let mut p = CoupledParticle3D::new([0.0; 3], [1.0, 0.0, 0.0], 1.0, 0.1, 0);
        p.force = [0.0; 3];
        sys.add_particle(p);
        // Fluid at rest → drag opposes particle velocity
        sys.apply_stokes_drag([0.0; 3]);
        assert!(
            sys.particles[0].force[0] < 0.0,
            "drag should be in -x for +x velocity"
        );
    }

    #[test]
    fn particle_system_centre_of_mass_symmetric() {
        let mut sys = ParticleSystem3D::new([0.0; 3], 1e-3, 1000.0);
        let mut p1 = make_3d_particle();
        p1.position = [-1.0, 0.0, 0.0];
        let mut p2 = make_3d_particle();
        p2.position = [1.0, 0.0, 0.0];
        sys.add_particle(p1);
        sys.add_particle(p2);
        let com = sys.centre_of_mass();
        assert!(com[0].abs() < 1e-14, "com.x = {}", com[0]);
    }

    #[test]
    fn integrate_rotations_function() {
        let mut particles = vec![make_3d_particle()];
        let i = particles[0].moment_of_inertia();
        particles[0].torque = [i, 0.0, 0.0];
        integrate_rotations(&mut particles, 0.1);
        assert!((particles[0].angular_velocity[0] - 0.1).abs() < 1e-12);
    }

    #[test]
    fn integrate_translations_function() {
        let mut particles = vec![make_3d_particle()];
        particles[0].force = [1.0, 0.0, 0.0]; // F=1, m=1 → a=1
        integrate_translations(&mut particles, 0.5);
        assert!((particles[0].velocity[0] - 0.5).abs() < 1e-12);
        assert!((particles[0].position[0] - 0.25).abs() < 1e-12); // x = v*dt = 0.5*0.5
    }
}
