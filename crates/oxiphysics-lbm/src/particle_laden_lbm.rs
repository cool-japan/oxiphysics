// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Particle-laden flow simulation using the Lattice Boltzmann Method.
//!
//! # Overview
//!
//! This module provides Lagrangian particle tracking coupled with an LBM
//! velocity field, including:
//!
//! - **Lagrangian particle tracking** in LBM velocity fields
//! - **Stokes drag** on particles
//! - **Saffman lift force** for particles in shear flow
//! - **Magnus force** for rotating particles
//! - **Particle-particle collision** via soft-sphere DEM
//! - **Two-way coupling** (particle back-reaction on fluid)
//! - **Particle size distributions** (log-normal, Rosin-Rammler)
//! - **Settling velocity** computation
//! - **Particle clustering and agglomeration**
//! - **Erosion and deposition** models
//! - **Residence time distribution** analysis

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Particle data structures
// ---------------------------------------------------------------------------

/// A Lagrangian particle tracked in the LBM flow field.
#[derive(Debug, Clone)]
pub struct LagrangianParticle {
    /// Position \[x, y, z\].
    pub pos: [f64; 3],
    /// Velocity \[vx, vy, vz\].
    pub vel: [f64; 3],
    /// Angular velocity \[wx, wy, wz\].
    pub angular_vel: [f64; 3],
    /// Particle diameter.
    pub diameter: f64,
    /// Particle density.
    pub density: f64,
    /// Particle mass (computed from density and diameter).
    pub mass: f64,
    /// Moment of inertia (for a sphere).
    pub inertia: f64,
    /// Force accumulator \[fx, fy, fz\].
    pub force: [f64; 3],
    /// Torque accumulator \[tx, ty, tz\].
    pub torque: [f64; 3],
    /// Temperature (optional, for heat transfer).
    pub temperature: f64,
    /// Unique particle ID.
    pub id: usize,
    /// Whether the particle is active.
    pub active: bool,
    /// Residence time (time since injection).
    pub residence_time: f64,
    /// Agglomerate ID (0 = none).
    pub agglomerate_id: usize,
}

impl LagrangianParticle {
    /// Create a new particle with given position, diameter, and density.
    pub fn new(pos: [f64; 3], diameter: f64, density: f64, id: usize) -> Self {
        let mass = density * PI * diameter.powi(3) / 6.0;
        let inertia = 0.1 * mass * diameter * diameter; // I = (1/10)*m*d^2 for sphere
        Self {
            pos,
            vel: [0.0; 3],
            angular_vel: [0.0; 3],
            diameter,
            density,
            mass,
            inertia,
            force: [0.0; 3],
            torque: [0.0; 3],
            temperature: 300.0,
            id,
            active: true,
            residence_time: 0.0,
            agglomerate_id: 0,
        }
    }

    /// Return the particle radius.
    pub fn radius(&self) -> f64 {
        self.diameter / 2.0
    }

    /// Return the particle volume.
    pub fn volume(&self) -> f64 {
        PI * self.diameter.powi(3) / 6.0
    }

    /// Return the particle cross-sectional area.
    pub fn cross_section_area(&self) -> f64 {
        PI * self.diameter * self.diameter / 4.0
    }

    /// Reset the force and torque accumulators.
    pub fn reset_forces(&mut self) {
        self.force = [0.0; 3];
        self.torque = [0.0; 3];
    }

    /// Return the speed (magnitude of velocity).
    pub fn speed(&self) -> f64 {
        (self.vel[0].powi(2) + self.vel[1].powi(2) + self.vel[2].powi(2)).sqrt()
    }
}

/// Parameters for the fluid environment.
#[derive(Debug, Clone)]
pub struct FluidParams {
    /// Fluid density.
    pub rho_f: f64,
    /// Fluid kinematic viscosity.
    pub nu: f64,
    /// Fluid dynamic viscosity (rho_f * nu).
    pub mu: f64,
    /// Gravitational acceleration \[gx, gy, gz\].
    pub gravity: [f64; 3],
}

impl FluidParams {
    /// Create fluid parameters.
    pub fn new(rho_f: f64, nu: f64, gravity: [f64; 3]) -> Self {
        Self {
            rho_f,
            nu,
            mu: rho_f * nu,
            gravity,
        }
    }

    /// Return the dynamic viscosity.
    pub fn dynamic_viscosity(&self) -> f64 {
        self.mu
    }
}

/// A simple 2D LBM velocity field for particle interpolation.
#[derive(Debug, Clone)]
pub struct VelocityField2D {
    /// Velocity x-component at grid points, row-major \[ny\]\[nx\].
    pub ux: Vec<f64>,
    /// Velocity y-component at grid points, row-major \[ny\]\[nx\].
    pub uy: Vec<f64>,
    /// Number of grid points in x.
    pub nx: usize,
    /// Number of grid points in y.
    pub ny: usize,
    /// Grid spacing.
    pub dx: f64,
}

impl VelocityField2D {
    /// Create a velocity field from given data.
    pub fn new(ux: Vec<f64>, uy: Vec<f64>, nx: usize, ny: usize, dx: f64) -> Self {
        Self { ux, uy, nx, ny, dx }
    }

    /// Create a uniform flow field.
    pub fn uniform(ux_val: f64, uy_val: f64, nx: usize, ny: usize, dx: f64) -> Self {
        let n = nx * ny;
        Self {
            ux: vec![ux_val; n],
            uy: vec![uy_val; n],
            nx,
            ny,
            dx,
        }
    }

    /// Bilinear interpolation of velocity at a point (x, y).
    pub fn interpolate(&self, x: f64, y: f64) -> [f64; 2] {
        let fi = x / self.dx;
        let fj = y / self.dx;
        let i0 = (fi.floor() as usize).min(self.nx.saturating_sub(2));
        let j0 = (fj.floor() as usize).min(self.ny.saturating_sub(2));
        let i1 = (i0 + 1).min(self.nx - 1);
        let j1 = (j0 + 1).min(self.ny - 1);
        let fx = fi - i0 as f64;
        let fy = fj - j0 as f64;

        let idx00 = j0 * self.nx + i0;
        let idx10 = j0 * self.nx + i1;
        let idx01 = j1 * self.nx + i0;
        let idx11 = j1 * self.nx + i1;

        let ux = self.ux[idx00] * (1.0 - fx) * (1.0 - fy)
            + self.ux[idx10] * fx * (1.0 - fy)
            + self.ux[idx01] * (1.0 - fx) * fy
            + self.ux[idx11] * fx * fy;

        let uy = self.uy[idx00] * (1.0 - fx) * (1.0 - fy)
            + self.uy[idx10] * fx * (1.0 - fy)
            + self.uy[idx01] * (1.0 - fx) * fy
            + self.uy[idx11] * fx * fy;

        [ux, uy]
    }

    /// Compute the velocity gradient tensor at a point using finite differences.
    ///
    /// Returns `[[du/dx, du/dy\], [dv/dx, dv/dy]]` as `[f64; 4]`.
    pub fn velocity_gradient(&self, x: f64, y: f64) -> [f64; 4] {
        let eps = self.dx * 0.5;
        let [ux_plus, uy_plus] = self.interpolate(x + eps, y);
        let [ux_minus, uy_minus] = self.interpolate(x - eps, y);
        let [ux_up, uy_up] = self.interpolate(x, y + eps);
        let [ux_down, uy_down] = self.interpolate(x, y - eps);

        let dudx = (ux_plus - ux_minus) / (2.0 * eps);
        let dudy = (ux_up - ux_down) / (2.0 * eps);
        let dvdx = (uy_plus - uy_minus) / (2.0 * eps);
        let dvdy = (uy_up - uy_down) / (2.0 * eps);

        [dudx, dudy, dvdx, dvdy]
    }
}

// ---------------------------------------------------------------------------
// Drag force models
// ---------------------------------------------------------------------------

/// Compute the particle Reynolds number.
pub fn particle_reynolds_number(rel_vel_mag: f64, diameter: f64, nu: f64) -> f64 {
    if nu < 1e-300 {
        return 0.0;
    }
    rel_vel_mag * diameter / nu
}

/// Compute the Stokes drag coefficient: C_D = 24/Re for Re << 1.
pub fn stokes_drag_coefficient(re_p: f64) -> f64 {
    if re_p < 1e-10 {
        return 24.0 / 1e-10;
    }
    24.0 / re_p
}

/// Compute the Schiller-Naumann drag coefficient (valid for Re < 1000).
///
/// C_D = (24/Re) * (1 + 0.15 * Re^0.687)
pub fn schiller_naumann_drag_coefficient(re_p: f64) -> f64 {
    if re_p < 1e-10 {
        return 24.0 / 1e-10;
    }
    (24.0 / re_p) * (1.0 + 0.15 * re_p.powf(0.687))
}

/// Compute the Stokes drag force on a spherical particle.
///
/// F_drag = 3 * pi * mu * d * (u_f - u_p)
///
/// Returns the drag force vector `[fx, fy, fz]`.
pub fn stokes_drag_force(
    fluid_vel: &[f64; 3],
    particle_vel: &[f64; 3],
    diameter: f64,
    mu: f64,
) -> [f64; 3] {
    let factor = 3.0 * PI * mu * diameter;
    [
        factor * (fluid_vel[0] - particle_vel[0]),
        factor * (fluid_vel[1] - particle_vel[1]),
        factor * (fluid_vel[2] - particle_vel[2]),
    ]
}

/// Compute the drag force using the Schiller-Naumann correlation.
///
/// Returns the drag force vector `[fx, fy, fz]`.
pub fn drag_force_schiller_naumann(
    fluid_vel: &[f64; 3],
    particle_vel: &[f64; 3],
    diameter: f64,
    rho_f: f64,
    nu: f64,
) -> [f64; 3] {
    let rel = [
        fluid_vel[0] - particle_vel[0],
        fluid_vel[1] - particle_vel[1],
        fluid_vel[2] - particle_vel[2],
    ];
    let rel_mag = (rel[0].powi(2) + rel[1].powi(2) + rel[2].powi(2)).sqrt();
    if rel_mag < 1e-30 {
        return [0.0; 3];
    }
    let re_p = particle_reynolds_number(rel_mag, diameter, nu);
    let cd = schiller_naumann_drag_coefficient(re_p);
    let area = PI * diameter * diameter / 4.0;
    let f_mag = 0.5 * cd * rho_f * area * rel_mag;
    [
        f_mag * rel[0] / rel_mag,
        f_mag * rel[1] / rel_mag,
        f_mag * rel[2] / rel_mag,
    ]
}

// ---------------------------------------------------------------------------
// Saffman lift force
// ---------------------------------------------------------------------------

/// Compute the Saffman lift force on a particle in shear flow.
///
/// F_saffman = 1.615 * mu * d * (u_f - u_p) * sqrt(rho_f * |omega| / mu)
///
/// where omega is the fluid vorticity (approximately du/dy for 2D shear).
///
/// Returns the lift force vector `[fx, fy, fz]`.
pub fn saffman_lift_force(
    fluid_vel: &[f64; 3],
    particle_vel: &[f64; 3],
    diameter: f64,
    rho_f: f64,
    mu: f64,
    shear_rate: f64,
) -> [f64; 3] {
    let rel = [
        fluid_vel[0] - particle_vel[0],
        fluid_vel[1] - particle_vel[1],
        fluid_vel[2] - particle_vel[2],
    ];
    if mu < 1e-300 || shear_rate.abs() < 1e-300 {
        return [0.0; 3];
    }
    let factor = 1.615 * mu * diameter * (rho_f * shear_rate.abs() / mu).sqrt();
    // Lift acts perpendicular to relative velocity (simplified: in y for 2D)
    // Cross product of shear direction (z-axis) and relative velocity
    [factor * (-rel[1]), factor * rel[0], 0.0]
}

// ---------------------------------------------------------------------------
// Magnus force (rotating particles)
// ---------------------------------------------------------------------------

/// Compute the Magnus force on a rotating spherical particle.
///
/// F_magnus = (pi/8) * rho_f * d^3 * (omega_rel x (u_f - u_p))
///
/// where omega_rel = omega_f - omega_p.
///
/// Returns the Magnus force vector `[fx, fy, fz]`.
pub fn magnus_force(
    fluid_vel: &[f64; 3],
    particle_vel: &[f64; 3],
    particle_angular_vel: &[f64; 3],
    fluid_angular_vel: &[f64; 3],
    diameter: f64,
    rho_f: f64,
) -> [f64; 3] {
    let rel = [
        fluid_vel[0] - particle_vel[0],
        fluid_vel[1] - particle_vel[1],
        fluid_vel[2] - particle_vel[2],
    ];
    let omega_rel = [
        fluid_angular_vel[0] - particle_angular_vel[0],
        fluid_angular_vel[1] - particle_angular_vel[1],
        fluid_angular_vel[2] - particle_angular_vel[2],
    ];
    // Cross product: omega_rel x rel
    let cross = [
        omega_rel[1] * rel[2] - omega_rel[2] * rel[1],
        omega_rel[2] * rel[0] - omega_rel[0] * rel[2],
        omega_rel[0] * rel[1] - omega_rel[1] * rel[0],
    ];
    let factor = PI / 8.0 * rho_f * diameter.powi(3);
    [factor * cross[0], factor * cross[1], factor * cross[2]]
}

// ---------------------------------------------------------------------------
// Buoyancy force
// ---------------------------------------------------------------------------

/// Compute the buoyancy force on a submerged particle.
///
/// F_buoy = (rho_p - rho_f) * V_p * g
pub fn buoyancy_force(particle: &LagrangianParticle, fluid: &FluidParams) -> [f64; 3] {
    let vol = particle.volume();
    let delta_rho = particle.density - fluid.rho_f;
    [
        -delta_rho * vol * fluid.gravity[0],
        -delta_rho * vol * fluid.gravity[1],
        -delta_rho * vol * fluid.gravity[2],
    ]
}

// ---------------------------------------------------------------------------
// Particle-particle collision (soft sphere DEM)
// ---------------------------------------------------------------------------

/// Soft-sphere DEM collision parameters.
#[derive(Debug, Clone)]
pub struct DemParams {
    /// Normal spring stiffness.
    pub kn: f64,
    /// Tangential spring stiffness.
    pub kt: f64,
    /// Normal damping coefficient.
    pub cn: f64,
    /// Tangential damping coefficient.
    pub ct: f64,
    /// Friction coefficient.
    pub mu_friction: f64,
}

impl DemParams {
    /// Create default DEM parameters.
    pub fn new(kn: f64, kt: f64, cn: f64, ct: f64, mu_friction: f64) -> Self {
        Self {
            kn,
            kt,
            cn,
            ct,
            mu_friction,
        }
    }
}

/// Compute the contact force between two particles using soft-sphere DEM.
///
/// Returns `(force_on_p1, force_on_p2)` as `([f64; 3], [f64; 3])`.
pub fn soft_sphere_collision(
    p1: &LagrangianParticle,
    p2: &LagrangianParticle,
    params: &DemParams,
) -> ([f64; 3], [f64; 3]) {
    let dx = [
        p2.pos[0] - p1.pos[0],
        p2.pos[1] - p1.pos[1],
        p2.pos[2] - p1.pos[2],
    ];
    let dist = (dx[0].powi(2) + dx[1].powi(2) + dx[2].powi(2)).sqrt();
    let contact_dist = p1.radius() + p2.radius();

    if dist >= contact_dist || dist < 1e-30 {
        return ([0.0; 3], [0.0; 3]);
    }

    let overlap = contact_dist - dist;
    let normal = [dx[0] / dist, dx[1] / dist, dx[2] / dist];

    // Relative velocity
    let rel_vel = [
        p1.vel[0] - p2.vel[0],
        p1.vel[1] - p2.vel[1],
        p1.vel[2] - p2.vel[2],
    ];
    let vn = rel_vel[0] * normal[0] + rel_vel[1] * normal[1] + rel_vel[2] * normal[2];

    // Normal force: spring + damping
    let fn_mag = params.kn * overlap - params.cn * vn;
    let fn_mag = fn_mag.max(0.0); // Only repulsive

    // Tangential velocity
    let vt = [
        rel_vel[0] - vn * normal[0],
        rel_vel[1] - vn * normal[1],
        rel_vel[2] - vn * normal[2],
    ];
    let vt_mag = (vt[0].powi(2) + vt[1].powi(2) + vt[2].powi(2)).sqrt();

    // Tangential force (limited by friction)
    let ft_mag = (params.kt * overlap + params.ct * vt_mag).min(params.mu_friction * fn_mag);

    let mut force1 = [0.0; 3];
    let mut force2 = [0.0; 3];
    for i in 0..3 {
        let fn_i = fn_mag * normal[i];
        let ft_i = if vt_mag > 1e-30 {
            -ft_mag * vt[i] / vt_mag
        } else {
            0.0
        };
        force1[i] = -(fn_i + ft_i);
        force2[i] = fn_i + ft_i;
    }

    (force1, force2)
}

/// Check all particle pairs for collisions and accumulate forces.
pub fn compute_dem_collisions(particles: &mut [LagrangianParticle], params: &DemParams) {
    let n = particles.len();
    // Collect forces first to avoid borrow issues
    let mut forces = vec![[0.0_f64; 3]; n];
    for i in 0..n {
        if !particles[i].active {
            continue;
        }
        for j in (i + 1)..n {
            if !particles[j].active {
                continue;
            }
            let (f1, f2) = soft_sphere_collision(&particles[i], &particles[j], params);
            for k in 0..3 {
                forces[i][k] += f1[k];
                forces[j][k] += f2[k];
            }
        }
    }
    for (p, pf) in particles.iter_mut().zip(forces.iter()) {
        for (fk, &pfk) in p.force.iter_mut().zip(pf.iter()) {
            *fk += pfk;
        }
    }
}

// ---------------------------------------------------------------------------
// Two-way coupling
// ---------------------------------------------------------------------------

/// Compute the reaction force from a particle on the fluid.
///
/// By Newton's third law, the force on the fluid is the negative of the
/// force on the particle.
pub fn particle_reaction_force(particle: &LagrangianParticle) -> [f64; 3] {
    [-particle.force[0], -particle.force[1], -particle.force[2]]
}

/// Distribute a point force to the nearby grid nodes using bilinear interpolation.
///
/// Returns a list of (grid_index, force_x, force_y) tuples.
pub fn distribute_force_to_grid(
    pos: [f64; 2],
    force: [f64; 2],
    nx: usize,
    ny: usize,
    dx: f64,
) -> Vec<(usize, f64, f64)> {
    let fi = pos[0] / dx;
    let fj = pos[1] / dx;
    let i0 = (fi.floor() as usize).min(nx.saturating_sub(2));
    let j0 = (fj.floor() as usize).min(ny.saturating_sub(2));
    let i1 = (i0 + 1).min(nx - 1);
    let j1 = (j0 + 1).min(ny - 1);
    let fx = fi - i0 as f64;
    let fy = fj - j0 as f64;

    let w00 = (1.0 - fx) * (1.0 - fy);
    let w10 = fx * (1.0 - fy);
    let w01 = (1.0 - fx) * fy;
    let w11 = fx * fy;

    vec![
        (j0 * nx + i0, force[0] * w00, force[1] * w00),
        (j0 * nx + i1, force[0] * w10, force[1] * w10),
        (j1 * nx + i0, force[0] * w01, force[1] * w01),
        (j1 * nx + i1, force[0] * w11, force[1] * w11),
    ]
}

// ---------------------------------------------------------------------------
// Particle size distributions
// ---------------------------------------------------------------------------

/// Log-normal probability density function.
///
/// f(d) = (1/(d*sigma*sqrt(2*pi))) * exp(-(ln(d)-mu)^2 / (2*sigma^2))
pub fn lognormal_pdf(d: f64, mu: f64, sigma: f64) -> f64 {
    if d <= 0.0 || sigma <= 0.0 {
        return 0.0;
    }
    let ln_d = d.ln();
    let exponent = -((ln_d - mu).powi(2)) / (2.0 * sigma * sigma);
    exponent.exp() / (d * sigma * (2.0 * PI).sqrt())
}

/// Log-normal cumulative distribution function (numerical approximation).
pub fn lognormal_cdf(d: f64, mu: f64, sigma: f64) -> f64 {
    if d <= 0.0 {
        return 0.0;
    }
    if sigma <= 0.0 {
        return if d.ln() >= mu { 1.0 } else { 0.0 };
    }
    let z = (d.ln() - mu) / sigma;
    0.5 * (1.0 + erf_approx(z / std::f64::consts::SQRT_2))
}

/// Approximate error function using Abramowitz & Stegun formula.
fn erf_approx(x: f64) -> f64 {
    let sign = x.signum();
    let x = x.abs();
    let t = 1.0 / (1.0 + 0.3275911 * x);
    let poly = t
        * (0.254829592
            + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
    sign * (1.0 - poly * (-x * x).exp())
}

/// Rosin-Rammler cumulative distribution function.
///
/// R(d) = 1 - exp(-(d/d_star)^n)
///
/// where `d_star` is the characteristic diameter and `n` is the spread parameter.
pub fn rosin_rammler_cdf(d: f64, d_star: f64, n: f64) -> f64 {
    if d <= 0.0 || d_star <= 0.0 {
        return 0.0;
    }
    1.0 - (-(d / d_star).powf(n)).exp()
}

/// Rosin-Rammler probability density function.
pub fn rosin_rammler_pdf(d: f64, d_star: f64, n: f64) -> f64 {
    if d <= 0.0 || d_star <= 0.0 {
        return 0.0;
    }
    let ratio = d / d_star;
    (n / d_star) * ratio.powf(n - 1.0) * (-ratio.powf(n)).exp()
}

/// Generate diameters from a log-normal distribution using the inverse CDF.
///
/// `n_particles` samples between quantiles 0.01 and 0.99.
pub fn lognormal_diameters(mu: f64, sigma: f64, n_particles: usize) -> Vec<f64> {
    let mut diameters = Vec::with_capacity(n_particles);
    for i in 0..n_particles {
        let p = 0.01 + 0.98 * (i as f64 + 0.5) / n_particles as f64;
        // Inverse CDF: d = exp(mu + sigma * sqrt(2) * erfinv(2p-1))
        let z = erfinv_approx(2.0 * p - 1.0);
        let d = (mu + sigma * std::f64::consts::SQRT_2 * z).exp();
        diameters.push(d);
    }
    diameters
}

/// Approximate inverse error function using a rational approximation.
fn erfinv_approx(x: f64) -> f64 {
    if x.abs() >= 1.0 {
        return x.signum() * 6.0;
    }
    let a = 0.147;
    let ln_term = ((1.0 - x * x).max(1e-300)).ln();
    let term1 = 2.0 / (PI * a) + ln_term / 2.0;
    let result = (term1 * term1 - ln_term / a).sqrt() - term1;
    x.signum() * result.max(0.0).sqrt()
}

/// Compute the mean diameter from a set of particle diameters.
pub fn mean_diameter(diameters: &[f64]) -> f64 {
    if diameters.is_empty() {
        return 0.0;
    }
    diameters.iter().sum::<f64>() / diameters.len() as f64
}

/// Compute the Sauter mean diameter (d32 = sum(d^3) / sum(d^2)).
pub fn sauter_mean_diameter(diameters: &[f64]) -> f64 {
    let sum_d3: f64 = diameters.iter().map(|&d| d.powi(3)).sum();
    let sum_d2: f64 = diameters.iter().map(|&d| d.powi(2)).sum();
    if sum_d2 < 1e-300 {
        return 0.0;
    }
    sum_d3 / sum_d2
}

// ---------------------------------------------------------------------------
// Settling velocity
// ---------------------------------------------------------------------------

/// Compute the Stokes settling velocity for a sphere.
///
/// v_s = (rho_p - rho_f) * g * d^2 / (18 * mu)
pub fn stokes_settling_velocity(diameter: f64, rho_p: f64, rho_f: f64, mu: f64, g: f64) -> f64 {
    if mu < 1e-300 {
        return 0.0;
    }
    (rho_p - rho_f) * g * diameter * diameter / (18.0 * mu)
}

/// Compute the settling velocity with drag correction (Schiller-Naumann).
///
/// Iterative solution for particles with Re_p > 1.
pub fn settling_velocity_corrected(
    diameter: f64,
    rho_p: f64,
    rho_f: f64,
    nu: f64,
    g: f64,
    max_iter: usize,
) -> f64 {
    let mu = rho_f * nu;
    let mut vs = stokes_settling_velocity(diameter, rho_p, rho_f, mu, g);

    for _iter in 0..max_iter {
        let re_p = particle_reynolds_number(vs.abs(), diameter, nu);
        let cd = schiller_naumann_drag_coefficient(re_p);
        // Force balance: (4/3)*d*g*(rho_p-rho_f) = cd*rho_f*vs^2/2
        let vs_new_sq = (4.0 * diameter * g * (rho_p - rho_f).abs()) / (3.0 * cd * rho_f);
        let vs_new = vs_new_sq.max(0.0).sqrt();
        if (vs_new - vs.abs()).abs() < 1e-10 * vs.abs().max(1e-30) {
            return vs_new * (rho_p - rho_f).signum();
        }
        vs = vs_new * (rho_p - rho_f).signum();
    }
    vs
}

// ---------------------------------------------------------------------------
// Particle clustering / agglomeration
// ---------------------------------------------------------------------------

/// Identify particle clusters using a distance threshold.
///
/// Returns a vector of agglomerate IDs (one per particle).
pub fn identify_clusters(particles: &[LagrangianParticle], distance_threshold: f64) -> Vec<usize> {
    let n = particles.len();
    let mut cluster_ids = vec![0_usize; n];
    let mut next_cluster = 1_usize;

    for i in 0..n {
        if !particles[i].active || cluster_ids[i] != 0 {
            continue;
        }
        cluster_ids[i] = next_cluster;
        // BFS to find all connected particles
        let mut queue = vec![i];
        let mut head = 0;
        while head < queue.len() {
            let current = queue[head];
            head += 1;
            for j in 0..n {
                if j == current || !particles[j].active || cluster_ids[j] != 0 {
                    continue;
                }
                let dx = particles[j].pos[0] - particles[current].pos[0];
                let dy = particles[j].pos[1] - particles[current].pos[1];
                let dz = particles[j].pos[2] - particles[current].pos[2];
                let dist = (dx * dx + dy * dy + dz * dz).sqrt();
                if dist < distance_threshold {
                    cluster_ids[j] = next_cluster;
                    queue.push(j);
                }
            }
        }
        next_cluster += 1;
    }
    cluster_ids
}

/// Compute the number of clusters and their sizes.
pub fn cluster_statistics(cluster_ids: &[usize]) -> Vec<usize> {
    if cluster_ids.is_empty() {
        return vec![];
    }
    let max_id = *cluster_ids.iter().max().unwrap_or(&0);
    let mut sizes = vec![0_usize; max_id + 1];
    for &id in cluster_ids {
        if id > 0 {
            sizes[id] += 1;
        }
    }
    // Return only non-zero cluster sizes
    sizes.into_iter().filter(|&s| s > 0).collect()
}

/// Simple agglomeration: merge two particles into one.
///
/// Returns a new particle at the center of mass with combined mass.
pub fn agglomerate_particles(
    p1: &LagrangianParticle,
    p2: &LagrangianParticle,
    new_id: usize,
) -> LagrangianParticle {
    let total_mass = p1.mass + p2.mass;
    let inv_mass = if total_mass > 0.0 {
        1.0 / total_mass
    } else {
        0.0
    };

    let pos = [
        (p1.pos[0] * p1.mass + p2.pos[0] * p2.mass) * inv_mass,
        (p1.pos[1] * p1.mass + p2.pos[1] * p2.mass) * inv_mass,
        (p1.pos[2] * p1.mass + p2.pos[2] * p2.mass) * inv_mass,
    ];
    let vel = [
        (p1.vel[0] * p1.mass + p2.vel[0] * p2.mass) * inv_mass,
        (p1.vel[1] * p1.mass + p2.vel[1] * p2.mass) * inv_mass,
        (p1.vel[2] * p1.mass + p2.vel[2] * p2.mass) * inv_mass,
    ];

    // Combined volume -> equivalent diameter
    let combined_vol = p1.volume() + p2.volume();
    let eq_diameter = (6.0 * combined_vol / PI).cbrt();
    let density = total_mass / combined_vol;

    let mut result = LagrangianParticle::new(pos, eq_diameter, density, new_id);
    result.vel = vel;
    result
}

// ---------------------------------------------------------------------------
// Erosion / Deposition models
// ---------------------------------------------------------------------------

/// Compute the Shields parameter for sediment transport.
///
/// theta = tau_b / ((rho_s - rho_f) * g * d)
pub fn shields_parameter(
    wall_shear_stress: f64,
    rho_s: f64,
    rho_f: f64,
    g: f64,
    diameter: f64,
) -> f64 {
    let denom = (rho_s - rho_f) * g * diameter;
    if denom.abs() < 1e-300 {
        return 0.0;
    }
    wall_shear_stress / denom
}

/// Check if erosion occurs based on the Shields criterion.
///
/// Critical Shields parameter is approximately 0.047 for turbulent flow.
pub fn shields_erosion_check(
    wall_shear_stress: f64,
    rho_s: f64,
    rho_f: f64,
    g: f64,
    diameter: f64,
    theta_critical: f64,
) -> bool {
    shields_parameter(wall_shear_stress, rho_s, rho_f, g, diameter) > theta_critical
}

/// Compute the erosion rate using the Partheniades formula.
///
/// E = M * (tau/tau_c - 1) for tau > tau_c, else 0.
pub fn partheniades_erosion_rate(
    wall_shear_stress: f64,
    critical_shear_stress: f64,
    erosion_constant: f64,
) -> f64 {
    if wall_shear_stress <= critical_shear_stress {
        return 0.0;
    }
    erosion_constant * (wall_shear_stress / critical_shear_stress - 1.0)
}

/// Compute the deposition rate using the Krone formula.
///
/// D = c * v_s * (1 - tau/tau_d) for tau < tau_d, else 0.
pub fn krone_deposition_rate(
    concentration: f64,
    settling_vel: f64,
    wall_shear_stress: f64,
    critical_deposition_stress: f64,
) -> f64 {
    if wall_shear_stress >= critical_deposition_stress {
        return 0.0;
    }
    concentration * settling_vel * (1.0 - wall_shear_stress / critical_deposition_stress)
}

// ---------------------------------------------------------------------------
// Residence time distribution
// ---------------------------------------------------------------------------

/// Compute the residence time distribution from exit times.
///
/// Returns `(time_bins, pdf_values)`.
pub fn residence_time_distribution(exit_times: &[f64], n_bins: usize) -> (Vec<f64>, Vec<f64>) {
    if exit_times.is_empty() || n_bins == 0 {
        return (vec![], vec![]);
    }
    let t_min = exit_times.iter().cloned().fold(f64::INFINITY, f64::min);
    let t_max = exit_times.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let range = t_max - t_min;
    if range < 1e-300 {
        return (vec![t_min], vec![1.0]);
    }
    let bin_width = range / n_bins as f64;

    let mut counts = vec![0_usize; n_bins];
    for &t in exit_times {
        let bin = ((t - t_min) / bin_width).floor() as usize;
        let bin = bin.min(n_bins - 1);
        counts[bin] += 1;
    }

    let n_total = exit_times.len() as f64;
    let times: Vec<f64> = (0..n_bins)
        .map(|i| t_min + (i as f64 + 0.5) * bin_width)
        .collect();
    let pdf: Vec<f64> = counts
        .iter()
        .map(|&c| c as f64 / (n_total * bin_width))
        .collect();

    (times, pdf)
}

/// Compute the mean residence time.
pub fn mean_residence_time(exit_times: &[f64]) -> f64 {
    if exit_times.is_empty() {
        return 0.0;
    }
    exit_times.iter().sum::<f64>() / exit_times.len() as f64
}

/// Compute the variance of residence times.
pub fn residence_time_variance(exit_times: &[f64]) -> f64 {
    if exit_times.len() < 2 {
        return 0.0;
    }
    let mean = mean_residence_time(exit_times);
    let n = exit_times.len() as f64;
    exit_times.iter().map(|&t| (t - mean).powi(2)).sum::<f64>() / (n - 1.0)
}

// ---------------------------------------------------------------------------
// Particle time integration
// ---------------------------------------------------------------------------

/// Advance a particle by one Euler step.
pub fn euler_step(p: &mut LagrangianParticle, dt: f64) {
    if !p.active {
        return;
    }
    // Velocity update: v += (F/m) * dt
    if p.mass > 1e-300 {
        for i in 0..3 {
            p.vel[i] += (p.force[i] / p.mass) * dt;
        }
    }
    // Position update: x += v * dt
    for i in 0..3 {
        p.pos[i] += p.vel[i] * dt;
    }
    // Angular velocity update
    if p.inertia > 1e-300 {
        for i in 0..3 {
            p.angular_vel[i] += (p.torque[i] / p.inertia) * dt;
        }
    }
    p.residence_time += dt;
}

/// Advance a particle using the velocity Verlet method.
pub fn verlet_step(p: &mut LagrangianParticle, old_force: &[f64; 3], dt: f64) {
    if !p.active || p.mass < 1e-300 {
        return;
    }
    let inv_mass = 1.0 / p.mass;
    // Position: x += v*dt + 0.5*(F/m)*dt^2
    for ((pos_k, &vel_k), &of_k) in p.pos.iter_mut().zip(p.vel.iter()).zip(old_force.iter()) {
        *pos_k += vel_k * dt + 0.5 * of_k * inv_mass * dt * dt;
    }
    // Velocity: v += 0.5*(F_old + F_new)/m * dt
    for ((vel_k, &of_k), &nf_k) in p.vel.iter_mut().zip(old_force.iter()).zip(p.force.iter()) {
        *vel_k += 0.5 * (of_k + nf_k) * inv_mass * dt;
    }
    p.residence_time += dt;
}

/// Track all particles for one time step in a 2D velocity field.
///
/// Applies Stokes drag and gravity, then advances with Euler integration.
pub fn track_particles_2d(
    particles: &mut [LagrangianParticle],
    field: &VelocityField2D,
    fluid: &FluidParams,
    dt: f64,
) {
    for p in particles.iter_mut() {
        if !p.active {
            continue;
        }
        p.reset_forces();

        // Interpolate fluid velocity at particle position
        let uv = field.interpolate(p.pos[0], p.pos[1]);
        let fluid_vel = [uv[0], uv[1], 0.0];

        // Stokes drag
        let drag = stokes_drag_force(&fluid_vel, &p.vel, p.diameter, fluid.mu);
        for (fk, &dk) in p.force.iter_mut().zip(drag.iter()) {
            *fk += dk;
        }

        // Buoyancy + gravity
        let buoy = buoyancy_force(p, fluid);
        for (fk, &bk) in p.force.iter_mut().zip(buoy.iter()) {
            *fk += bk;
        }

        // Integrate
        euler_step(p, dt);
    }
}

// ---------------------------------------------------------------------------
// Particle Stokes number
// ---------------------------------------------------------------------------

/// Compute the particle Stokes number.
///
/// St = tau_p / tau_f, where tau_p = rho_p * d^2 / (18 * mu).
pub fn stokes_number(diameter: f64, rho_p: f64, mu: f64, char_time: f64) -> f64 {
    if mu < 1e-300 || char_time < 1e-300 {
        return 0.0;
    }
    let tau_p = rho_p * diameter * diameter / (18.0 * mu);
    tau_p / char_time
}

/// Compute the particle relaxation time.
///
/// tau_p = rho_p * d^2 / (18 * mu).
pub fn particle_relaxation_time(diameter: f64, rho_p: f64, mu: f64) -> f64 {
    if mu < 1e-300 {
        return 0.0;
    }
    rho_p * diameter * diameter / (18.0 * mu)
}

// ---------------------------------------------------------------------------
// Concentration and volume fraction
// ---------------------------------------------------------------------------

/// Compute the local particle volume fraction on a 2D grid.
///
/// Returns a grid-sized vector of volume fractions.
pub fn compute_volume_fraction(
    particles: &[LagrangianParticle],
    nx: usize,
    ny: usize,
    dx: f64,
) -> Vec<f64> {
    let n = nx * ny;
    let mut vf = vec![0.0; n];
    let cell_vol = dx * dx * dx; // 3D cell volume

    for p in particles {
        if !p.active {
            continue;
        }
        let i = (p.pos[0] / dx).floor() as usize;
        let j = (p.pos[1] / dx).floor() as usize;
        if i < nx && j < ny {
            let idx = j * nx + i;
            vf[idx] += p.volume() / cell_vol;
        }
    }
    vf
}

/// Compute the total mass loading ratio (particle mass / fluid mass in domain).
pub fn mass_loading_ratio(particles: &[LagrangianParticle], rho_f: f64, domain_volume: f64) -> f64 {
    let total_particle_mass: f64 = particles.iter().filter(|p| p.active).map(|p| p.mass).sum();
    let fluid_mass = rho_f * domain_volume;
    if fluid_mass < 1e-300 {
        return 0.0;
    }
    total_particle_mass / fluid_mass
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- Particle tests ---

    #[test]
    fn test_particle_creation() {
        let p = LagrangianParticle::new([1.0, 2.0, 3.0], 0.01, 2500.0, 0);
        assert!(p.mass > 0.0);
        assert!(p.active);
        assert_eq!(p.id, 0);
    }

    #[test]
    fn test_particle_radius() {
        let p = LagrangianParticle::new([0.0; 3], 0.1, 1000.0, 0);
        assert!((p.radius() - 0.05).abs() < 1e-10);
    }

    #[test]
    fn test_particle_volume() {
        let p = LagrangianParticle::new([0.0; 3], 1.0, 1000.0, 0);
        let expected = PI / 6.0;
        assert!((p.volume() - expected).abs() < 1e-10);
    }

    #[test]
    fn test_particle_cross_section() {
        let p = LagrangianParticle::new([0.0; 3], 1.0, 1000.0, 0);
        let expected = PI / 4.0;
        assert!((p.cross_section_area() - expected).abs() < 1e-10);
    }

    #[test]
    fn test_particle_speed() {
        let mut p = LagrangianParticle::new([0.0; 3], 0.01, 1000.0, 0);
        p.vel = [3.0, 4.0, 0.0];
        assert!((p.speed() - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_particle_reset_forces() {
        let mut p = LagrangianParticle::new([0.0; 3], 0.01, 1000.0, 0);
        p.force = [1.0, 2.0, 3.0];
        p.reset_forces();
        assert_eq!(p.force, [0.0, 0.0, 0.0]);
    }

    // --- Fluid params tests ---

    #[test]
    fn test_fluid_params() {
        let fp = FluidParams::new(1000.0, 1e-6, [0.0, -9.81, 0.0]);
        assert!((fp.mu - 1e-3).abs() < 1e-10);
        assert!((fp.dynamic_viscosity() - 1e-3).abs() < 1e-10);
    }

    // --- Velocity field tests ---

    #[test]
    fn test_uniform_velocity_field() {
        let field = VelocityField2D::uniform(1.0, 0.0, 10, 10, 0.1);
        let v = field.interpolate(0.5, 0.5);
        assert!((v[0] - 1.0).abs() < 1e-10);
        assert!((v[1]).abs() < 1e-10);
    }

    #[test]
    fn test_velocity_gradient_uniform() {
        let field = VelocityField2D::uniform(1.0, 0.0, 20, 20, 0.1);
        let grad = field.velocity_gradient(0.5, 0.5);
        // Uniform field -> zero gradient
        assert!(grad[0].abs() < 1e-6);
        assert!(grad[1].abs() < 1e-6);
    }

    // --- Drag force tests ---

    #[test]
    fn test_stokes_drag_coefficient() {
        let cd = stokes_drag_coefficient(1.0);
        assert!((cd - 24.0).abs() < 1e-10);
    }

    #[test]
    fn test_schiller_naumann_reduces_to_stokes() {
        // At very low Re, S-N should approach Stokes
        let cd_sn = schiller_naumann_drag_coefficient(0.001);
        let cd_stokes = stokes_drag_coefficient(0.001);
        assert!((cd_sn - cd_stokes).abs() / cd_stokes < 0.01);
    }

    #[test]
    fn test_stokes_drag_force_direction() {
        let drag = stokes_drag_force(&[1.0, 0.0, 0.0], &[0.0, 0.0, 0.0], 0.01, 1e-3);
        assert!(drag[0] > 0.0, "drag should be in fluid direction");
        assert!(drag[1].abs() < 1e-15);
    }

    #[test]
    fn test_stokes_drag_zero_relative_velocity() {
        let drag = stokes_drag_force(&[1.0, 2.0, 3.0], &[1.0, 2.0, 3.0], 0.01, 1e-3);
        assert!(drag[0].abs() < 1e-15);
        assert!(drag[1].abs() < 1e-15);
        assert!(drag[2].abs() < 1e-15);
    }

    #[test]
    fn test_drag_force_schiller_naumann() {
        let drag =
            drag_force_schiller_naumann(&[1.0, 0.0, 0.0], &[0.0, 0.0, 0.0], 0.01, 1000.0, 1e-6);
        assert!(drag[0] > 0.0);
    }

    // --- Saffman lift tests ---

    #[test]
    fn test_saffman_lift_nonzero() {
        let lift = saffman_lift_force(&[1.0, 0.0, 0.0], &[0.0, 0.0, 0.0], 0.01, 1000.0, 1e-3, 10.0);
        // Lift should be perpendicular to flow
        assert!(lift[1].abs() > 0.0 || lift[0].abs() > 0.0);
    }

    #[test]
    fn test_saffman_lift_zero_shear() {
        let lift = saffman_lift_force(&[1.0, 0.0, 0.0], &[0.0, 0.0, 0.0], 0.01, 1000.0, 1e-3, 0.0);
        assert_eq!(lift, [0.0, 0.0, 0.0]);
    }

    // --- Magnus force tests ---

    #[test]
    fn test_magnus_force_no_rotation() {
        let f = magnus_force(
            &[1.0, 0.0, 0.0],
            &[0.0, 0.0, 0.0],
            &[0.0, 0.0, 0.0],
            &[0.0, 0.0, 0.0],
            0.01,
            1000.0,
        );
        assert_eq!(f, [0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_magnus_force_with_spin() {
        let f = magnus_force(
            &[1.0, 0.0, 0.0],
            &[0.0, 0.0, 0.0],
            &[0.0, 0.0, 10.0], // spinning around z
            &[0.0, 0.0, 0.0],
            0.01,
            1000.0,
        );
        // omega_rel = [0,0,-10], rel = [1,0,0]
        // cross = [-10*0 - 0*0, 0*1 - (-10)*0, (-10)*0 - 0*1] = [0, 0, 0]
        // Actually: omega_rel x rel = [oy*rz - oz*ry, oz*rx - ox*rz, ox*ry - oy*rx]
        // = [0*0 - (-10)*0, (-10)*1 - 0*0, 0*0 - 0*1] = [0, -10, 0]
        assert!(f[1].abs() > 0.0, "should have y-component");
    }

    // --- Buoyancy test ---

    #[test]
    fn test_buoyancy_heavier_particle() {
        let p = LagrangianParticle::new([0.0; 3], 0.01, 2500.0, 0);
        let fluid = FluidParams::new(1000.0, 1e-6, [0.0, -9.81, 0.0]);
        let f = buoyancy_force(&p, &fluid);
        // Heavier particle -> buoyancy should push upward (opposing gravity)
        // delta_rho = 2500-1000 > 0, F = -delta_rho*V*g => Fy = +delta_rho*V*9.81
        assert!(
            f[1] > 0.0,
            "buoyancy should oppose gravity for heavy particle"
        );
    }

    // --- DEM collision tests ---

    #[test]
    fn test_no_collision_distant_particles() {
        let p1 = LagrangianParticle::new([0.0, 0.0, 0.0], 0.01, 1000.0, 0);
        let p2 = LagrangianParticle::new([1.0, 0.0, 0.0], 0.01, 1000.0, 1);
        let params = DemParams::new(1000.0, 500.0, 10.0, 5.0, 0.3);
        let (f1, f2) = soft_sphere_collision(&p1, &p2, &params);
        assert_eq!(f1, [0.0, 0.0, 0.0]);
        assert_eq!(f2, [0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_collision_overlapping_particles() {
        let p1 = LagrangianParticle::new([0.0, 0.0, 0.0], 0.1, 1000.0, 0);
        let p2 = LagrangianParticle::new([0.08, 0.0, 0.0], 0.1, 1000.0, 1);
        let params = DemParams::new(1000.0, 500.0, 10.0, 5.0, 0.3);
        let (f1, f2) = soft_sphere_collision(&p1, &p2, &params);
        // Particles overlap, should repel
        assert!(f1[0] < 0.0, "p1 should be pushed left");
        assert!(f2[0] > 0.0, "p2 should be pushed right");
    }

    #[test]
    fn test_compute_dem_collisions_runs() {
        let mut particles = vec![
            LagrangianParticle::new([0.0, 0.0, 0.0], 0.1, 1000.0, 0),
            LagrangianParticle::new([0.08, 0.0, 0.0], 0.1, 1000.0, 1),
        ];
        let params = DemParams::new(1000.0, 500.0, 10.0, 5.0, 0.3);
        compute_dem_collisions(&mut particles, &params);
        // Forces should be nonzero due to overlap
        assert!(particles[0].force[0].abs() > 0.0);
    }

    // --- Two-way coupling tests ---

    #[test]
    fn test_reaction_force_newton_third_law() {
        let mut p = LagrangianParticle::new([0.0; 3], 0.01, 1000.0, 0);
        p.force = [1.0, -2.0, 3.0];
        let reaction = particle_reaction_force(&p);
        assert_eq!(reaction, [-1.0, 2.0, -3.0]);
    }

    #[test]
    fn test_distribute_force_weights_sum_to_one() {
        let nodes = distribute_force_to_grid([0.15, 0.15], [1.0, 0.0], 10, 10, 0.1);
        let total_fx: f64 = nodes.iter().map(|n| n.1).sum();
        assert!((total_fx - 1.0).abs() < 1e-10, "weights should sum to 1");
    }

    // --- Size distribution tests ---

    #[test]
    fn test_lognormal_pdf_positive() {
        let pdf = lognormal_pdf(0.01, -4.0_f64, 0.5);
        assert!(pdf > 0.0);
    }

    #[test]
    fn test_lognormal_pdf_zero_for_negative() {
        assert_eq!(lognormal_pdf(-1.0, 0.0, 1.0), 0.0);
    }

    #[test]
    fn test_lognormal_cdf_monotone() {
        let c1 = lognormal_cdf(0.01, 0.0, 1.0);
        let c2 = lognormal_cdf(0.1, 0.0, 1.0);
        let c3 = lognormal_cdf(1.0, 0.0, 1.0);
        assert!(c1 < c2);
        assert!(c2 < c3);
    }

    #[test]
    fn test_rosin_rammler_cdf_limits() {
        assert!(rosin_rammler_cdf(0.0, 1.0, 2.0) < 1e-10);
        assert!((rosin_rammler_cdf(100.0, 1.0, 2.0) - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_rosin_rammler_pdf_positive() {
        let pdf = rosin_rammler_pdf(0.5, 1.0, 2.0);
        assert!(pdf > 0.0);
    }

    #[test]
    fn test_lognormal_diameters_count() {
        let diams = lognormal_diameters(0.0, 0.5, 20);
        assert_eq!(diams.len(), 20);
        for &d in &diams {
            assert!(d > 0.0);
        }
    }

    #[test]
    fn test_sauter_mean_diameter() {
        let diameters = vec![1.0, 1.0, 1.0];
        let d32 = sauter_mean_diameter(&diameters);
        assert!((d32 - 1.0).abs() < 1e-10);
    }

    // --- Settling velocity tests ---

    #[test]
    fn test_stokes_settling_positive_for_heavy_particle() {
        let vs = stokes_settling_velocity(0.001, 2500.0, 1000.0, 1e-3, 9.81);
        assert!(vs > 0.0);
    }

    #[test]
    fn test_stokes_settling_negative_for_light_particle() {
        let vs = stokes_settling_velocity(0.001, 500.0, 1000.0, 1e-3, 9.81);
        assert!(vs < 0.0);
    }

    #[test]
    fn test_settling_velocity_corrected() {
        let vs = settling_velocity_corrected(0.001, 2500.0, 1000.0, 1e-6, 9.81, 50);
        assert!(vs > 0.0);
    }

    // --- Clustering tests ---

    #[test]
    fn test_identify_clusters_separate() {
        let particles = vec![
            LagrangianParticle::new([0.0, 0.0, 0.0], 0.01, 1000.0, 0),
            LagrangianParticle::new([10.0, 0.0, 0.0], 0.01, 1000.0, 1),
        ];
        let ids = identify_clusters(&particles, 0.1);
        assert_ne!(ids[0], ids[1]);
    }

    #[test]
    fn test_identify_clusters_close() {
        let particles = vec![
            LagrangianParticle::new([0.0, 0.0, 0.0], 0.01, 1000.0, 0),
            LagrangianParticle::new([0.05, 0.0, 0.0], 0.01, 1000.0, 1),
        ];
        let ids = identify_clusters(&particles, 0.1);
        assert_eq!(ids[0], ids[1]);
    }

    #[test]
    fn test_cluster_statistics() {
        let ids = vec![1, 1, 1, 2, 2, 3];
        let sizes = cluster_statistics(&ids);
        assert!(sizes.contains(&3));
        assert!(sizes.contains(&2));
        assert!(sizes.contains(&1));
    }

    #[test]
    fn test_agglomerate_particles() {
        let p1 = LagrangianParticle::new([0.0, 0.0, 0.0], 0.1, 1000.0, 0);
        let p2 = LagrangianParticle::new([0.1, 0.0, 0.0], 0.1, 1000.0, 1);
        let merged = agglomerate_particles(&p1, &p2, 10);
        assert!(merged.mass > p1.mass);
        assert!(merged.diameter > 0.1);
        assert_eq!(merged.id, 10);
    }

    // --- Erosion/deposition tests ---

    #[test]
    fn test_shields_parameter() {
        let theta = shields_parameter(0.5, 2650.0, 1000.0, 9.81, 0.001);
        assert!(theta > 0.0);
    }

    #[test]
    fn test_shields_erosion_check() {
        assert!(shields_erosion_check(
            1.0, 2650.0, 1000.0, 9.81, 0.001, 0.047
        ));
        assert!(!shields_erosion_check(
            0.001, 2650.0, 1000.0, 9.81, 0.001, 0.047
        ));
    }

    #[test]
    fn test_partheniades_no_erosion_below_threshold() {
        let rate = partheniades_erosion_rate(0.5, 1.0, 0.001);
        assert_eq!(rate, 0.0);
    }

    #[test]
    fn test_partheniades_erosion_above_threshold() {
        let rate = partheniades_erosion_rate(2.0, 1.0, 0.001);
        assert!(rate > 0.0);
    }

    #[test]
    fn test_krone_deposition() {
        let rate = krone_deposition_rate(0.1, 0.01, 0.3, 1.0);
        assert!(rate > 0.0);
    }

    #[test]
    fn test_krone_no_deposition_above_threshold() {
        let rate = krone_deposition_rate(0.1, 0.01, 1.5, 1.0);
        assert_eq!(rate, 0.0);
    }

    // --- Residence time tests ---

    #[test]
    fn test_residence_time_distribution() {
        let times = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let (bins, pdf) = residence_time_distribution(&times, 4);
        assert_eq!(bins.len(), 4);
        assert_eq!(pdf.len(), 4);
        let total: f64 = pdf.iter().sum::<f64>() * (4.0 / 4.0);
        assert!((total - 1.0).abs() < 0.5, "PDF should integrate to ~1");
    }

    #[test]
    fn test_mean_residence_time() {
        let times = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let mean = mean_residence_time(&times);
        assert!((mean - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_residence_time_variance() {
        let times = vec![2.0, 2.0, 2.0];
        let var = residence_time_variance(&times);
        assert!(var.abs() < 1e-10, "identical times -> zero variance");
    }

    // --- Integration tests ---

    #[test]
    fn test_euler_step_free_fall() {
        let mut p = LagrangianParticle::new([0.0, 10.0, 0.0], 0.01, 1000.0, 0);
        p.force = [0.0, -p.mass * 9.81, 0.0];
        euler_step(&mut p, 0.01);
        assert!(p.vel[1] < 0.0, "should accelerate downward");
        assert!(p.pos[1] < 10.0, "should move down");
    }

    #[test]
    fn test_verlet_step() {
        let mut p = LagrangianParticle::new([0.0; 3], 0.01, 1000.0, 0);
        let old_f = [0.0, -p.mass * 9.81, 0.0];
        p.force = [0.0, -p.mass * 9.81, 0.0];
        verlet_step(&mut p, &old_f, 0.01);
        assert!(p.vel[1] < 0.0);
    }

    #[test]
    fn test_track_particles_2d() {
        let field = VelocityField2D::uniform(0.1, 0.0, 20, 20, 0.1);
        let fluid = FluidParams::new(1000.0, 1e-6, [0.0, -9.81, 0.0]);
        let mut particles = vec![LagrangianParticle::new([0.5, 0.5, 0.0], 0.001, 2500.0, 0)];
        track_particles_2d(&mut particles, &field, &fluid, 0.001);
        // Particle should have moved
        assert!(particles[0].residence_time > 0.0);
    }

    // --- Stokes number tests ---

    #[test]
    fn test_stokes_number() {
        let st = stokes_number(0.001, 2500.0, 1e-3, 0.1);
        assert!(st > 0.0);
    }

    #[test]
    fn test_particle_relaxation_time() {
        let tau = particle_relaxation_time(0.001, 2500.0, 1e-3);
        // tau = 2500 * 1e-6 / 18e-3 = 2.5e-3 / 18e-3 = 0.1389
        assert!((tau - 2500.0 * 1e-6 / (18.0 * 1e-3)).abs() < 1e-10);
    }

    // --- Volume fraction tests ---

    #[test]
    fn test_volume_fraction_single_particle() {
        let particles = vec![LagrangianParticle::new([0.05, 0.05, 0.0], 0.01, 1000.0, 0)];
        let vf = compute_volume_fraction(&particles, 10, 10, 0.1);
        let total: f64 = vf.iter().sum();
        assert!(total > 0.0);
    }

    #[test]
    fn test_mass_loading_ratio() {
        let particles = vec![LagrangianParticle::new([0.0; 3], 0.01, 2500.0, 0)];
        let ml = mass_loading_ratio(&particles, 1000.0, 1.0);
        assert!(ml > 0.0);
        assert!(ml < 1.0, "single small particle should have small loading");
    }

    // --- Particle Reynolds number test ---

    #[test]
    fn test_particle_reynolds_number() {
        let re = particle_reynolds_number(0.1, 0.001, 1e-6);
        assert!((re - 100.0).abs() < 1e-6);
    }
}
