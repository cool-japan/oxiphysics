// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Immersed Boundary Method (IBM) coupled to LBM flow fields.
//!
//! Implements the Peskin/Roma IBM for fluid-structure interaction:
//! - Lagrangian marker points (`IbmParticle`) with forces and Peskin delta kernel.
//! - Eulerian grid (`IbmFluidGrid`) for force spreading and velocity interpolation.
//! - Direct forcing IBM (`DirectForcingIbm`) for no-slip enforcement.
//! - Rigid moving objects (`MovingIbmObject`) with velocity/rotation update.
//! - Flexible elastic filaments (`FlexibleIbmFiber`) with bending and tension.
//! - Force/torque statistics (`IbmCouplingStats`) with drag/lift and power.

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Discrete delta functions (Peskin / Roma kernels)
// ---------------------------------------------------------------------------

/// Roma 3-point discrete delta function φ(r).
///
/// # Arguments
/// * `r` — non-dimensional distance (distance / grid_spacing)
///
/// Returns the weight for force spreading / velocity interpolation.
pub fn delta_phi(r: f64) -> f64 {
    let r = r.abs();
    if r <= 0.5 {
        (1.0 + (PI * 2.0 * r).cos() / PI) / 3.0
    } else if r <= 1.5 {
        let x = 2.0 * r - 1.0;
        (5.0 + 4.0 * (PI * x).cos() / PI - (PI * (4.0 * x + 2.0)).cos() * 2.0_f64.sqrt() / PI)
            / 12.0
    } else {
        0.0
    }
}

/// Peskin 4-point cosine delta function.
///
/// Provides a smooth, compact-support kernel (support radius 2 lattice units).
///
/// # Arguments
/// * `r` — non-dimensional distance (distance / grid_spacing)
pub fn delta_peskin4(r: f64) -> f64 {
    let r = r.abs();
    if r <= 1.0 {
        (3.0 - 2.0 * r + (2.0 * PI * r).sin() / PI) / 8.0
    } else if r <= 2.0 {
        (5.0 - 2.0 * r - (2.0 * PI * (r - 1.0)).sin() / PI) / 8.0
    } else {
        0.0
    }
}

/// 2-D product delta: φ(rx) * φ(ry) using the Roma 3-point kernel.
///
/// # Arguments
/// * `rx`, `ry` — non-dimensional distances in x and y
pub fn delta_2d(rx: f64, ry: f64) -> f64 {
    delta_phi(rx) * delta_phi(ry)
}

/// 3-D product delta: φ(rx) * φ(ry) * φ(rz) using the Roma 3-point kernel.
///
/// # Arguments
/// * `rx`, `ry`, `rz` — non-dimensional distances
pub fn delta_3d(rx: f64, ry: f64, rz: f64) -> f64 {
    delta_phi(rx) * delta_phi(ry) * delta_phi(rz)
}

// ---------------------------------------------------------------------------
// IbmParticle — Lagrangian marker
// ---------------------------------------------------------------------------

/// A Lagrangian marker point on an immersed boundary.
///
/// Carries position, rest position, velocity, and body force.
/// Equivalent to `IbmMarker` but with richer semantics (mass, arc-length weight).
#[derive(Clone, Debug, PartialEq)]
pub struct IbmParticle {
    /// Current 2-D position \[x, y\] in lattice units.
    pub position: [f64; 2],
    /// Rest (reference) position \[x, y\] for tether spring.
    pub rest_position: [f64; 2],
    /// Current velocity \[ux, uy\] of this marker.
    pub velocity: [f64; 2],
    /// Lagrangian body force \[fx, fy\] at this marker (to be spread to Eulerian grid).
    pub force: [f64; 2],
    /// Arc-length quadrature weight Δs (spacing between neighbouring markers).
    pub ds: f64,
    /// Marker mass (lattice units).
    pub mass: f64,
}

impl IbmParticle {
    /// Create a new `IbmParticle` at `(x, y)` with default arc-length weight 1.
    pub fn new(x: f64, y: f64) -> Self {
        Self {
            position: [x, y],
            rest_position: [x, y],
            velocity: [0.0, 0.0],
            force: [0.0, 0.0],
            ds: 1.0,
            mass: 1.0,
        }
    }

    /// Create a particle with explicit rest position and arc-length weight.
    pub fn with_rest_and_ds(x: f64, y: f64, rx: f64, ry: f64, ds: f64) -> Self {
        Self {
            position: [x, y],
            rest_position: [rx, ry],
            velocity: [0.0, 0.0],
            force: [0.0, 0.0],
            ds,
            mass: 1.0,
        }
    }

    /// Displacement from rest position.
    pub fn displacement(&self) -> [f64; 2] {
        [
            self.position[0] - self.rest_position[0],
            self.position[1] - self.rest_position[1],
        ]
    }

    /// Squared distance from rest position.
    pub fn displacement_sq(&self) -> f64 {
        let d = self.displacement();
        d[0] * d[0] + d[1] * d[1]
    }

    /// Euclidean distance from rest position.
    pub fn displacement_norm(&self) -> f64 {
        self.displacement_sq().sqrt()
    }

    /// Advance marker position by explicit Euler: x += dt * v.
    pub fn advect(&mut self, dt: f64) {
        self.position[0] += dt * self.velocity[0];
        self.position[1] += dt * self.velocity[1];
    }

    /// Apply a velocity correction `dv` to the marker velocity.
    pub fn apply_velocity_correction(&mut self, dv: [f64; 2]) {
        self.velocity[0] += dv[0];
        self.velocity[1] += dv[1];
    }
}

// ---------------------------------------------------------------------------
// Legacy IbmMarker (kept for backward compatibility)
// ---------------------------------------------------------------------------

/// A Lagrangian marker point attached to an immersed boundary (legacy struct).
#[derive(Clone, Debug, PartialEq)]
pub struct IbmMarker {
    /// Current position \[x, y\] in lattice units.
    pub position: [f64; 2],
    /// Rest (equilibrium) position \[x, y\] for tether spring.
    pub rest_position: [f64; 2],
    /// Current velocity \[ux, uy\] of the marker.
    pub velocity: [f64; 2],
    /// Lagrangian body force \[fx, fy\] at this marker.
    pub force: [f64; 2],
}

impl IbmMarker {
    /// Create a new `IbmMarker` at `(x, y)` with zero velocity and force.
    pub fn new(x: f64, y: f64) -> Self {
        Self {
            position: [x, y],
            rest_position: [x, y],
            velocity: [0.0, 0.0],
            force: [0.0, 0.0],
        }
    }

    /// Create a new `IbmMarker` with a separate rest position.
    pub fn with_rest(x: f64, y: f64, rx: f64, ry: f64) -> Self {
        Self {
            position: [x, y],
            rest_position: [rx, ry],
            velocity: [0.0, 0.0],
            force: [0.0, 0.0],
        }
    }
}

// ---------------------------------------------------------------------------
// IbmFluidGrid — Eulerian grid for force/velocity exchange
// ---------------------------------------------------------------------------

/// Eulerian fluid grid for the immersed boundary method.
///
/// Stores Eulerian velocity and body-force fields.
/// The grid is row-major with index `iy * nx + ix`.
pub struct IbmFluidGrid {
    /// Grid width (number of cells in x).
    pub nx: usize,
    /// Grid height (number of cells in y).
    pub ny: usize,
    /// x-component of Eulerian velocity field.
    pub ux: Vec<f64>,
    /// y-component of Eulerian velocity field.
    pub uy: Vec<f64>,
    /// x-component of Eulerian body-force field (updated by IBM).
    pub fx: Vec<f64>,
    /// y-component of Eulerian body-force field (updated by IBM).
    pub fy: Vec<f64>,
    /// Lattice spacing Δx (used as ds reference).
    pub dx: f64,
}

impl IbmFluidGrid {
    /// Construct a new `IbmFluidGrid` of size `nx × ny` with zero fields and unit spacing.
    pub fn new(nx: usize, ny: usize) -> Self {
        let n = nx * ny;
        Self {
            nx,
            ny,
            ux: vec![0.0; n],
            uy: vec![0.0; n],
            fx: vec![0.0; n],
            fy: vec![0.0; n],
            dx: 1.0,
        }
    }

    /// Total number of grid cells.
    pub fn ncells(&self) -> usize {
        self.nx * self.ny
    }

    /// Zero the Eulerian force fields in place.
    pub fn clear_forces(&mut self) {
        for v in self.fx.iter_mut() {
            *v = 0.0;
        }
        for v in self.fy.iter_mut() {
            *v = 0.0;
        }
    }

    /// Set a uniform Eulerian velocity field.
    pub fn set_uniform_velocity(&mut self, ux: f64, uy: f64) {
        for v in self.ux.iter_mut() {
            *v = ux;
        }
        for v in self.uy.iter_mut() {
            *v = uy;
        }
    }

    /// Spread Lagrangian forces from `particles` onto the Eulerian body-force fields.
    ///
    /// Uses the Roma 3-point delta function with each particle's own `ds` weight.
    pub fn spread_forces(&mut self, particles: &[IbmParticle]) {
        for p in particles {
            let xi = p.position[0];
            let yi = p.position[1];
            let ix0 = (xi - 1.5).floor().max(0.0) as usize;
            let ix1 = ((xi + 1.5).ceil() as usize).min(self.nx.saturating_sub(1));
            let iy0 = (yi - 1.5).floor().max(0.0) as usize;
            let iy1 = ((yi + 1.5).ceil() as usize).min(self.ny.saturating_sub(1));
            for iy in iy0..=iy1 {
                let dy = yi - iy as f64;
                let wy = delta_phi(dy);
                for ix in ix0..=ix1 {
                    let dx = xi - ix as f64;
                    let w = delta_phi(dx) * wy * p.ds;
                    let idx = iy * self.nx + ix;
                    self.fx[idx] += p.force[0] * w;
                    self.fy[idx] += p.force[1] * w;
                }
            }
        }
    }

    /// Interpolate the Eulerian velocity field to each Lagrangian marker.
    ///
    /// Uses the Roma 3-point delta function; updates `particle.velocity` in place.
    pub fn interpolate_velocities(&self, particles: &mut [IbmParticle]) {
        for p in particles.iter_mut() {
            let xi = p.position[0];
            let yi = p.position[1];
            let ix0 = (xi - 1.5).floor().max(0.0) as usize;
            let ix1 = ((xi + 1.5).ceil() as usize).min(self.nx.saturating_sub(1));
            let iy0 = (yi - 1.5).floor().max(0.0) as usize;
            let iy1 = ((yi + 1.5).ceil() as usize).min(self.ny.saturating_sub(1));
            let mut vx = 0.0_f64;
            let mut vy = 0.0_f64;
            for iy in iy0..=iy1 {
                let dy = yi - iy as f64;
                let wy = delta_phi(dy);
                for ix in ix0..=ix1 {
                    let dx = xi - ix as f64;
                    let w = delta_phi(dx) * wy * p.ds;
                    let idx = iy * self.nx + ix;
                    vx += self.ux[idx] * w;
                    vy += self.uy[idx] * w;
                }
            }
            p.velocity = [vx, vy];
        }
    }

    /// Peskin-4 variant: spread forces using the 4-point cosine kernel.
    pub fn spread_forces_peskin4(&mut self, particles: &[IbmParticle]) {
        for p in particles {
            let xi = p.position[0];
            let yi = p.position[1];
            let ix0 = (xi - 2.0).floor().max(0.0) as usize;
            let ix1 = ((xi + 2.0).ceil() as usize).min(self.nx.saturating_sub(1));
            let iy0 = (yi - 2.0).floor().max(0.0) as usize;
            let iy1 = ((yi + 2.0).ceil() as usize).min(self.ny.saturating_sub(1));
            for iy in iy0..=iy1 {
                let dy = yi - iy as f64;
                let wy = delta_peskin4(dy);
                for ix in ix0..=ix1 {
                    let dx = xi - ix as f64;
                    let w = delta_peskin4(dx) * wy * p.ds;
                    let idx = iy * self.nx + ix;
                    self.fx[idx] += p.force[0] * w;
                    self.fy[idx] += p.force[1] * w;
                }
            }
        }
    }

    /// Net force on the fluid: sum of all Eulerian body-force components.
    pub fn total_force(&self) -> [f64; 2] {
        let fx: f64 = self.fx.iter().sum();
        let fy: f64 = self.fy.iter().sum();
        [fx, fy]
    }
}

// ---------------------------------------------------------------------------
// Force spreading / interpolation (free functions — legacy API)
// ---------------------------------------------------------------------------

/// Spread marker forces onto a 2D Eulerian grid using the Roma 3-point delta function.
///
/// # Arguments
/// * `markers` — slice of Lagrangian markers
/// * `fx_grid` — Eulerian x-force array (modified in place)
/// * `fy_grid` — Eulerian y-force array (modified in place)
/// * `nx`, `ny` — grid dimensions
/// * `ds` — arc-length spacing between markers (for quadrature weight)
pub fn spread_force(
    markers: &[IbmMarker],
    fx_grid: &mut [f64],
    fy_grid: &mut [f64],
    nx: usize,
    ny: usize,
    ds: f64,
) {
    for m in markers {
        let xi = m.position[0];
        let yi = m.position[1];
        let ix0 = (xi - 1.5).floor().max(0.0) as usize;
        let ix1 = ((xi + 1.5).ceil() as usize).min(nx - 1);
        let iy0 = (yi - 1.5).floor().max(0.0) as usize;
        let iy1 = ((yi + 1.5).ceil() as usize).min(ny - 1);
        for iy in iy0..=iy1 {
            let dy = yi - iy as f64;
            let wy = delta_phi(dy);
            for ix in ix0..=ix1 {
                let dx = xi - ix as f64;
                let w = delta_phi(dx) * wy * ds;
                let idx = iy * nx + ix;
                fx_grid[idx] += m.force[0] * w;
                fy_grid[idx] += m.force[1] * w;
            }
        }
    }
}

/// Interpolate the Eulerian velocity field to each Lagrangian marker.
///
/// # Arguments
/// * `markers` — mutable slice; `velocity` field is updated
/// * `ux_grid` — Eulerian x-velocity (length `ny * nx`)
/// * `uy_grid` — Eulerian y-velocity (length `ny * nx`)
/// * `nx`, `ny` — grid dimensions
/// * `ds` — arc-length spacing for the quadrature
pub fn interpolate_velocity(
    markers: &mut [IbmMarker],
    ux_grid: &[f64],
    uy_grid: &[f64],
    nx: usize,
    ny: usize,
    ds: f64,
) {
    for m in markers.iter_mut() {
        let xi = m.position[0];
        let yi = m.position[1];
        let ix0 = (xi - 1.5).floor().max(0.0) as usize;
        let ix1 = ((xi + 1.5).ceil() as usize).min(nx - 1);
        let iy0 = (yi - 1.5).floor().max(0.0) as usize;
        let iy1 = ((yi + 1.5).ceil() as usize).min(ny - 1);
        let mut vx = 0.0_f64;
        let mut vy = 0.0_f64;
        for iy in iy0..=iy1 {
            let dy = yi - iy as f64;
            let wy = delta_phi(dy);
            for ix in ix0..=ix1 {
                let dx = xi - ix as f64;
                let w = delta_phi(dx) * wy * ds;
                let idx = iy * nx + ix;
                vx += ux_grid[idx] * w;
                vy += uy_grid[idx] * w;
            }
        }
        m.velocity = [vx, vy];
    }
}

/// Compute the tether spring force pulling each marker toward its rest position.
///
/// `F = -k_s * (x - x_rest)`
///
/// # Arguments
/// * `markers` — mutable slice; `force` field is updated
/// * `k_spring` — tether spring constant
pub fn compute_elastic_force(markers: &mut [IbmMarker], k_spring: f64) {
    for m in markers.iter_mut() {
        m.force[0] = -k_spring * (m.position[0] - m.rest_position[0]);
        m.force[1] = -k_spring * (m.position[1] - m.rest_position[1]);
    }
}

/// Advance the IBM-LBM system by one time step (legacy free-function API).
///
/// Steps:
/// 1. Zero Eulerian force fields.
/// 2. Compute elastic tether forces on markers.
/// 3. Spread forces to Eulerian grid.
/// 4. Interpolate Eulerian velocity to markers.
/// 5. Advance marker positions with explicit Euler.
pub fn ib_lbm_step(
    markers: &mut [IbmMarker],
    fx_grid: &mut [f64],
    fy_grid: &mut [f64],
    ux_grid: &[f64],
    uy_grid: &[f64],
    nx: usize,
    ny: usize,
    ds: f64,
    k_spring: f64,
    dt: f64,
) {
    for v in fx_grid.iter_mut() {
        *v = 0.0;
    }
    for v in fy_grid.iter_mut() {
        *v = 0.0;
    }
    compute_elastic_force(markers, k_spring);
    spread_force(markers, fx_grid, fy_grid, nx, ny, ds);
    interpolate_velocity(markers, ux_grid, uy_grid, nx, ny, ds);
    for m in markers.iter_mut() {
        m.position[0] += dt * m.velocity[0];
        m.position[1] += dt * m.velocity[1];
    }
}

// ---------------------------------------------------------------------------
// DirectForcingIbm — no-slip enforcement via body force
// ---------------------------------------------------------------------------

/// Direct-forcing IBM: computes Lagrangian body forces that enforce no-slip
/// on an immersed boundary by driving marker velocity toward a target.
///
/// Uses the Fadlun et al. / Uhlmann direct forcing approach:
/// `f_IBM = (u_target - u_fluid) / dt`
pub struct DirectForcingIbm {
    /// Lagrangian marker particles.
    pub particles: Vec<IbmParticle>,
    /// Time step Δt.
    pub dt: f64,
    /// Density of the fluid (lattice units).
    pub rho_fluid: f64,
}

impl DirectForcingIbm {
    /// Create a new `DirectForcingIbm` with the given particles and time step.
    pub fn new(particles: Vec<IbmParticle>, dt: f64, rho_fluid: f64) -> Self {
        Self {
            particles,
            dt,
            rho_fluid,
        }
    }

    /// Set the target (desired) velocity for all markers simultaneously.
    ///
    /// Useful for rigid body no-slip: every marker shares the same target.
    pub fn set_target_velocity(&self, target: [f64; 2]) -> Vec<[f64; 2]> {
        vec![target; self.particles.len()]
    }

    /// Compute direct forcing body forces from target velocities and fluid velocities.
    ///
    /// `f_i = rho * (u_target_i - u_fluid_i) / dt`
    ///
    /// # Arguments
    /// * `u_target` — desired velocity at each marker
    /// * `u_fluid`  — interpolated Eulerian velocity at each marker
    pub fn compute_body_forces(&mut self, u_target: &[[f64; 2]], u_fluid: &[[f64; 2]]) {
        for (i, p) in self.particles.iter_mut().enumerate() {
            p.force[0] = self.rho_fluid * (u_target[i][0] - u_fluid[i][0]) / self.dt;
            p.force[1] = self.rho_fluid * (u_target[i][1] - u_fluid[i][1]) / self.dt;
        }
    }

    /// Full direct-forcing IBM step.
    ///
    /// 1. Interpolate Eulerian velocity to markers.
    /// 2. Compute body forces to enforce `u_target`.
    /// 3. Spread body forces back to Eulerian grid.
    ///
    /// # Arguments
    /// * `grid` — Eulerian fluid grid (velocity fields read, force fields updated)
    /// * `u_target` — target velocity for each marker
    pub fn step(&mut self, grid: &mut IbmFluidGrid, u_target: &[[f64; 2]]) {
        // Step 1: interpolate fluid velocity to markers.
        grid.interpolate_velocities(&mut self.particles);
        let u_fluid: Vec<[f64; 2]> = self.particles.iter().map(|p| p.velocity).collect();
        // Step 2: compute direct forcing forces.
        self.compute_body_forces(u_target, &u_fluid);
        // Step 3: spread forces to Eulerian grid.
        grid.clear_forces();
        grid.spread_forces(&self.particles);
    }

    /// Number of Lagrangian markers.
    pub fn num_markers(&self) -> usize {
        self.particles.len()
    }

    /// Total Lagrangian force (sum over all markers).
    pub fn total_lagrangian_force(&self) -> [f64; 2] {
        self.particles.iter().fold([0.0_f64; 2], |acc, p| {
            [acc[0] + p.force[0], acc[1] + p.force[1]]
        })
    }
}

// ---------------------------------------------------------------------------
// MovingIbmObject — rigid body immersed in LBM fluid
// ---------------------------------------------------------------------------

/// A rigid body immersed in the LBM fluid, represented by a set of Lagrangian
/// markers that collectively define the object's surface.
///
/// Translational and rotational motion are integrated with explicit Euler.
pub struct MovingIbmObject {
    /// Lagrangian surface markers.
    pub particles: Vec<IbmParticle>,
    /// Centre-of-mass position \[x, y\].
    pub center: [f64; 2],
    /// Translational velocity \[ux, uy\].
    pub velocity: [f64; 2],
    /// Angular velocity ω (rad/lattice-step).
    pub omega: f64,
    /// Total mass of the rigid body.
    pub mass: f64,
    /// Moment of inertia about centre.
    pub inertia: f64,
    /// Applied external force (e.g. buoyancy) \[fx, fy\].
    pub external_force: [f64; 2],
    /// Applied external torque.
    pub external_torque: f64,
}

impl MovingIbmObject {
    /// Construct a circular rigid body of radius `r` with `n_markers` surface markers.
    ///
    /// Markers are placed uniformly on the circle centred at `(cx, cy)`.
    pub fn circle(cx: f64, cy: f64, r: f64, n_markers: usize, mass: f64, inertia: f64) -> Self {
        let ds = 2.0 * PI * r / n_markers as f64;
        let particles: Vec<IbmParticle> = (0..n_markers)
            .map(|i| {
                let theta = 2.0 * PI * i as f64 / n_markers as f64;
                let x = cx + r * theta.cos();
                let y = cy + r * theta.sin();
                IbmParticle::with_rest_and_ds(x, y, x, y, ds)
            })
            .collect();
        Self {
            particles,
            center: [cx, cy],
            velocity: [0.0, 0.0],
            omega: 0.0,
            mass,
            inertia,
            external_force: [0.0, 0.0],
            external_torque: 0.0,
        }
    }

    /// Compute net hydrodynamic force and torque on the body from marker forces.
    ///
    /// Returns `([fx, fy], torque)`.
    pub fn compute_hydrodynamic_force(&self) -> ([f64; 2], f64) {
        let mut fx = 0.0_f64;
        let mut fy = 0.0_f64;
        let mut torque = 0.0_f64;
        for p in &self.particles {
            fx -= p.force[0] * p.ds;
            fy -= p.force[1] * p.ds;
            let rx = p.position[0] - self.center[0];
            let ry = p.position[1] - self.center[1];
            torque -= (rx * p.force[1] - ry * p.force[0]) * p.ds;
        }
        ([fx, fy], torque)
    }

    /// Update rigid-body velocity and angular velocity using Newton's second law.
    ///
    /// `m * du/dt = F_hydro + F_ext`
    /// `I * dω/dt = T_hydro + T_ext`
    pub fn update_velocity(&mut self, dt: f64) {
        let (f_hydro, tau_hydro) = self.compute_hydrodynamic_force();
        let ax = (f_hydro[0] + self.external_force[0]) / self.mass;
        let ay = (f_hydro[1] + self.external_force[1]) / self.mass;
        let alpha = (tau_hydro + self.external_torque) / self.inertia;
        self.velocity[0] += dt * ax;
        self.velocity[1] += dt * ay;
        self.omega += dt * alpha;
    }

    /// Move centre of mass and rotate marker positions.
    pub fn update_position(&mut self, dt: f64) {
        self.center[0] += dt * self.velocity[0];
        self.center[1] += dt * self.velocity[1];
        let dtheta = dt * self.omega;
        let cos_dt = dtheta.cos();
        let sin_dt = dtheta.sin();
        for p in self.particles.iter_mut() {
            let rx = p.position[0] - (self.center[0] - dt * self.velocity[0]);
            let ry = p.position[1] - (self.center[1] - dt * self.velocity[1]);
            p.position[0] = self.center[0] + cos_dt * rx - sin_dt * ry;
            p.position[1] = self.center[1] + sin_dt * rx + cos_dt * ry;
            p.velocity[0] = self.velocity[0] - self.omega * ry;
            p.velocity[1] = self.velocity[1] + self.omega * rx;
        }
    }

    /// Velocity of the surface at a given marker index (rigid body kinematics).
    ///
    /// `u_surface = v_cm + ω × r`
    pub fn surface_velocity(&self, idx: usize) -> [f64; 2] {
        let p = &self.particles[idx];
        let rx = p.position[0] - self.center[0];
        let ry = p.position[1] - self.center[1];
        [
            self.velocity[0] - self.omega * ry,
            self.velocity[1] + self.omega * rx,
        ]
    }

    /// Number of surface markers.
    pub fn num_markers(&self) -> usize {
        self.particles.len()
    }

    /// Kinetic energy of the rigid body (translation + rotation).
    pub fn kinetic_energy(&self) -> f64 {
        let v2 = self.velocity[0].powi(2) + self.velocity[1].powi(2);
        0.5 * self.mass * v2 + 0.5 * self.inertia * self.omega.powi(2)
    }
}

// ---------------------------------------------------------------------------
// FlexibleIbmFiber — elastic filament in flow
// ---------------------------------------------------------------------------

/// An elastic filament immersed in the fluid, modelled as a chain of
/// `IbmParticle` markers connected by springs (tension) and bending moments.
///
/// Physical model:
/// - Tension force: `F_t = k_t * (|Δx| - Δs) * t̂`  (extensible spring)
/// - Bending force: derived from bending energy `E_b = k_b/2 * ∑ (κ - κ_0)^2`
///   discretised via finite differences on the Lagrangian mesh.
pub struct FlexibleIbmFiber {
    /// Lagrangian marker particles along the fiber.
    pub particles: Vec<IbmParticle>,
    /// Tensile stiffness k_t (force per unit strain).
    pub k_tension: f64,
    /// Bending stiffness k_b (energy per unit curvature-squared per unit length).
    pub k_bending: f64,
    /// Natural (rest) arc-length spacing Δs₀.
    pub ds0: f64,
    /// Natural curvature κ₀ (zero for a straight fiber).
    pub kappa0: f64,
}

impl FlexibleIbmFiber {
    /// Construct a straight horizontal fiber with `n` markers.
    ///
    /// The fiber runs from `(x0, y)` to `(x0 + (n-1)*ds, y)`.
    pub fn straight(x0: f64, y: f64, n: usize, ds: f64, k_tension: f64, k_bending: f64) -> Self {
        let particles: Vec<IbmParticle> = (0..n)
            .map(|i| {
                IbmParticle::with_rest_and_ds(x0 + i as f64 * ds, y, x0 + i as f64 * ds, y, ds)
            })
            .collect();
        Self {
            particles,
            k_tension,
            k_bending,
            ds0: ds,
            kappa0: 0.0,
        }
    }

    /// Compute tension forces between adjacent markers and accumulate into `force`.
    ///
    /// `F_i = k_t * (|x_{i+1} - x_i| - ds0) / |x_{i+1} - x_i| * (x_{i+1} - x_i)`
    pub fn compute_tension_forces(&mut self) {
        let n = self.particles.len();
        // Reset forces.
        for p in self.particles.iter_mut() {
            p.force = [0.0, 0.0];
        }
        for i in 0..n.saturating_sub(1) {
            let dx = self.particles[i + 1].position[0] - self.particles[i].position[0];
            let dy = self.particles[i + 1].position[1] - self.particles[i].position[1];
            let dist = (dx * dx + dy * dy).sqrt();
            if dist < 1e-14 {
                continue;
            }
            let strain = (dist - self.ds0) / dist;
            let ft = self.k_tension * strain;
            let ftx = ft * dx;
            let fty = ft * dy;
            self.particles[i].force[0] += ftx;
            self.particles[i].force[1] += fty;
            self.particles[i + 1].force[0] -= ftx;
            self.particles[i + 1].force[1] -= fty;
        }
    }

    /// Compute bending forces using a finite-difference curvature model.
    ///
    /// For interior node i:
    /// `F_bending_i ≈ -k_b * (x_{i+1} - 2*x_i + x_{i-1}) / ds0^4`
    pub fn compute_bending_forces(&mut self) {
        let n = self.particles.len();
        let kb = self.k_bending;
        let ds4 = self.ds0.powi(4);
        // Accumulate bending forces (add to existing tension forces).
        for i in 1..n.saturating_sub(1) {
            let xim1 = self.particles[i - 1].position;
            let xi = self.particles[i].position;
            let xip1 = self.particles[i + 1].position;
            let d2x = xip1[0] - 2.0 * xi[0] + xim1[0];
            let d2y = xip1[1] - 2.0 * xi[1] + xim1[1];
            // Curvature contribution to force on i.
            self.particles[i].force[0] -= kb * d2x / ds4;
            self.particles[i].force[1] -= kb * d2y / ds4;
            // Reaction on neighbours (central-difference stencil).
            self.particles[i - 1].force[0] += 0.5 * kb * d2x / ds4;
            self.particles[i - 1].force[1] += 0.5 * kb * d2y / ds4;
            self.particles[i + 1].force[0] += 0.5 * kb * d2x / ds4;
            self.particles[i + 1].force[1] += 0.5 * kb * d2y / ds4;
        }
    }

    /// Compute all elastic forces (tension + bending).
    pub fn compute_elastic_forces(&mut self) {
        self.compute_tension_forces();
        self.compute_bending_forces();
    }

    /// Total elastic potential energy of the fiber.
    pub fn elastic_energy(&self) -> f64 {
        let n = self.particles.len();
        let mut e = 0.0_f64;
        // Tension.
        for i in 0..n.saturating_sub(1) {
            let dx = self.particles[i + 1].position[0] - self.particles[i].position[0];
            let dy = self.particles[i + 1].position[1] - self.particles[i].position[1];
            let dist = (dx * dx + dy * dy).sqrt();
            let strain = dist - self.ds0;
            e += 0.5 * self.k_tension * strain * strain;
        }
        // Bending.
        for i in 1..n.saturating_sub(1) {
            let xim1 = self.particles[i - 1].position;
            let xi = self.particles[i].position;
            let xip1 = self.particles[i + 1].position;
            let d2x = xip1[0] - 2.0 * xi[0] + xim1[0];
            let d2y = xip1[1] - 2.0 * xi[1] + xim1[1];
            let kappa_sq = (d2x * d2x + d2y * d2y) / (self.ds0 * self.ds0).powi(2);
            e += 0.5 * self.k_bending * (kappa_sq - self.kappa0 * self.kappa0) * self.ds0;
        }
        e
    }

    /// Arc length of the fiber (sum of segment lengths).
    pub fn arc_length(&self) -> f64 {
        let n = self.particles.len();
        let mut len = 0.0_f64;
        for i in 0..n.saturating_sub(1) {
            let dx = self.particles[i + 1].position[0] - self.particles[i].position[0];
            let dy = self.particles[i + 1].position[1] - self.particles[i].position[1];
            len += (dx * dx + dy * dy).sqrt();
        }
        len
    }

    /// Tip deflection: distance of last marker from its rest position.
    pub fn tip_deflection(&self) -> f64 {
        if let Some(last) = self.particles.last() {
            last.displacement_norm()
        } else {
            0.0
        }
    }

    /// Number of markers.
    pub fn num_markers(&self) -> usize {
        self.particles.len()
    }
}

// ---------------------------------------------------------------------------
// IbmCouplingStats — force/torque and hydrodynamic coefficients
// ---------------------------------------------------------------------------

/// Accumulated hydrodynamic coupling statistics for an immersed object.
///
/// Records drag, lift, torque, and power over a simulation run and computes
/// non-dimensional force coefficients.
pub struct IbmCouplingStats {
    /// Free-stream velocity magnitude U∞ (for coefficient normalisation).
    pub u_inf: f64,
    /// Reference length (e.g. cylinder diameter D) for normalisation.
    pub ref_length: f64,
    /// Fluid density ρ.
    pub rho_fluid: f64,
    /// Accumulated drag force (in the x-direction).
    pub drag: f64,
    /// Accumulated lift force (in the y-direction).
    pub lift: f64,
    /// Accumulated torque (z-component).
    pub torque: f64,
    /// Accumulated power dissipation F⋅U.
    pub power: f64,
    /// Number of samples accumulated.
    pub n_samples: usize,
}

impl IbmCouplingStats {
    /// Create a new `IbmCouplingStats` with given reference quantities.
    pub fn new(u_inf: f64, ref_length: f64, rho_fluid: f64) -> Self {
        Self {
            u_inf,
            ref_length,
            rho_fluid,
            drag: 0.0,
            lift: 0.0,
            torque: 0.0,
            power: 0.0,
            n_samples: 0,
        }
    }

    /// Accumulate one sample of force/torque.
    ///
    /// # Arguments
    /// * `fx` — x-force (drag)
    /// * `fy` — y-force (lift)
    /// * `tau` — torque
    /// * `u_body` — body velocity (for power calculation)
    pub fn record(&mut self, fx: f64, fy: f64, tau: f64, u_body: [f64; 2]) {
        self.drag += fx;
        self.lift += fy;
        self.torque += tau;
        self.power += fx * u_body[0] + fy * u_body[1];
        self.n_samples += 1;
    }

    /// Mean drag over all recorded samples.
    pub fn mean_drag(&self) -> f64 {
        if self.n_samples == 0 {
            0.0
        } else {
            self.drag / self.n_samples as f64
        }
    }

    /// Mean lift over all recorded samples.
    pub fn mean_lift(&self) -> f64 {
        if self.n_samples == 0 {
            0.0
        } else {
            self.lift / self.n_samples as f64
        }
    }

    /// Mean torque over all recorded samples.
    pub fn mean_torque(&self) -> f64 {
        if self.n_samples == 0 {
            0.0
        } else {
            self.torque / self.n_samples as f64
        }
    }

    /// Mean power dissipation.
    pub fn mean_power(&self) -> f64 {
        if self.n_samples == 0 {
            0.0
        } else {
            self.power / self.n_samples as f64
        }
    }

    /// Drag coefficient: `CD = 2 * F_drag / (ρ U∞² L)`.
    pub fn drag_coefficient(&self) -> f64 {
        let q = 0.5 * self.rho_fluid * self.u_inf * self.u_inf * self.ref_length;
        if q.abs() < 1e-14 {
            0.0
        } else {
            self.mean_drag() / q
        }
    }

    /// Lift coefficient: `CL = 2 * F_lift / (ρ U∞² L)`.
    pub fn lift_coefficient(&self) -> f64 {
        let q = 0.5 * self.rho_fluid * self.u_inf * self.u_inf * self.ref_length;
        if q.abs() < 1e-14 {
            0.0
        } else {
            self.mean_lift() / q
        }
    }

    /// Torque coefficient: `CM = 2 * τ / (ρ U∞² L²)`.
    pub fn torque_coefficient(&self) -> f64 {
        let q = 0.5 * self.rho_fluid * self.u_inf * self.u_inf * self.ref_length * self.ref_length;
        if q.abs() < 1e-14 {
            0.0
        } else {
            self.mean_torque() / q
        }
    }

    /// Reset all accumulated statistics.
    pub fn reset(&mut self) {
        self.drag = 0.0;
        self.lift = 0.0;
        self.torque = 0.0;
        self.power = 0.0;
        self.n_samples = 0;
    }

    /// Strouhal number estimate from a known shedding frequency `f_shed`.
    ///
    /// `St = f_shed * L / U∞`
    pub fn strouhal_number(&self, f_shed: f64) -> f64 {
        if self.u_inf.abs() < 1e-14 {
            0.0
        } else {
            f_shed * self.ref_length / self.u_inf
        }
    }
}

// ---------------------------------------------------------------------------
// Utility helpers
// ---------------------------------------------------------------------------

/// Build a ring of `n` Lagrangian particles uniformly distributed on a circle
/// of radius `r` centred at `(cx, cy)`.
///
/// The arc-length weight is set to `2π r / n`.
pub fn make_circle_markers(cx: f64, cy: f64, r: f64, n: usize) -> Vec<IbmParticle> {
    let ds = 2.0 * PI * r / n as f64;
    (0..n)
        .map(|i| {
            let theta = 2.0 * PI * i as f64 / n as f64;
            let x = cx + r * theta.cos();
            let y = cy + r * theta.sin();
            IbmParticle::with_rest_and_ds(x, y, x, y, ds)
        })
        .collect()
}

/// Build a straight line of `n` Lagrangian particles from `(x0, y0)` to `(x1, y1)`.
///
/// The arc-length weight is set to `L / (n - 1)` where `L` is the total length.
pub fn make_line_markers(x0: f64, y0: f64, x1: f64, y1: f64, n: usize) -> Vec<IbmParticle> {
    if n < 2 {
        return vec![IbmParticle::new(x0, y0)];
    }
    let len = ((x1 - x0).powi(2) + (y1 - y0).powi(2)).sqrt();
    let ds = len / (n - 1) as f64;
    (0..n)
        .map(|i| {
            let t = i as f64 / (n - 1) as f64;
            let x = x0 + t * (x1 - x0);
            let y = y0 + t * (y1 - y0);
            IbmParticle::with_rest_and_ds(x, y, x, y, ds)
        })
        .collect()
}

/// Compute approximate curvature at the i-th interior marker using central differences.
///
/// Returns the signed curvature κ ≈ `|r'' × r'| / |r'|³`.
pub fn local_curvature(particles: &[IbmParticle], i: usize) -> f64 {
    if i == 0 || i + 1 >= particles.len() {
        return 0.0;
    }
    let xm = particles[i - 1].position;
    let x0 = particles[i].position;
    let xp = particles[i + 1].position;
    let r1x = (xp[0] - xm[0]) * 0.5;
    let r1y = (xp[1] - xm[1]) * 0.5;
    let r2x = xp[0] - 2.0 * x0[0] + xm[0];
    let r2y = xp[1] - 2.0 * x0[1] + xm[1];
    let cross = r1x * r2y - r1y * r2x;
    let r1_mag = (r1x * r1x + r1y * r1y).sqrt();
    if r1_mag < 1e-14 {
        0.0
    } else {
        cross / r1_mag.powi(3)
    }
}

/// Penalisation factor for diffuse interface IBM.
///
/// Returns `η / (η + ν * dt)` where η is the permeability parameter,
/// ν the kinematic viscosity, and dt the time step.
pub fn penalisation_factor(eta: f64, nu: f64, dt: f64) -> f64 {
    let denom = eta + nu * dt;
    if denom.abs() < 1e-14 {
        0.0
    } else {
        eta / denom
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ---- delta_phi ----

    #[test]
    fn test_delta_phi_center() {
        let v = delta_phi(0.0);
        assert!(v > 0.0 && v <= 1.0, "delta_phi(0) = {v}");
    }

    #[test]
    fn test_delta_phi_outside_support() {
        assert_eq!(delta_phi(1.6), 0.0);
        assert_eq!(delta_phi(-2.0), 0.0);
    }

    #[test]
    fn test_delta_phi_symmetry() {
        for r in [0.1, 0.5, 0.9, 1.2] {
            assert!((delta_phi(r) - delta_phi(-r)).abs() < 1e-14);
        }
    }

    #[test]
    fn test_delta_phi_mid_region() {
        let v = delta_phi(0.5);
        assert!(v > 0.0);
    }

    #[test]
    fn test_delta_phi_continuity() {
        let left = delta_phi(0.499);
        let right = delta_phi(0.501);
        assert!(left > 0.0, "left = {left}");
        assert!(right > 0.0, "right = {right}");
        assert!(left < 1.0);
        assert!(right < 1.0);
    }

    #[test]
    fn test_delta_phi_partition_of_unity() {
        let sum: f64 = (-2i64..=2).map(|k| delta_phi(k as f64)).sum();
        assert!(sum > 0.0 && sum < 3.0, "sum = {sum}");
    }

    // ---- delta_peskin4 ----

    #[test]
    fn test_delta_peskin4_center() {
        let v = delta_peskin4(0.0);
        assert!(v > 0.0, "peskin4(0) = {v}");
    }

    #[test]
    fn test_delta_peskin4_outside() {
        assert_eq!(delta_peskin4(2.1), 0.0);
        assert_eq!(delta_peskin4(-3.0), 0.0);
    }

    #[test]
    fn test_delta_peskin4_symmetry() {
        for r in [0.3, 0.7, 1.1, 1.8] {
            assert!((delta_peskin4(r) - delta_peskin4(-r)).abs() < 1e-14);
        }
    }

    #[test]
    fn test_delta_2d_zero_distance() {
        let v = delta_2d(0.0, 0.0);
        assert!(v > 0.0);
    }

    #[test]
    fn test_delta_3d_zero_distance() {
        let v = delta_3d(0.0, 0.0, 0.0);
        assert!(v > 0.0);
    }

    #[test]
    fn test_delta_3d_outside() {
        assert_eq!(delta_3d(0.0, 0.0, 2.0), 0.0);
    }

    // ---- IbmParticle ----

    #[test]
    fn test_ibm_particle_new() {
        let p = IbmParticle::new(3.0, 4.0);
        assert_eq!(p.position, [3.0, 4.0]);
        assert_eq!(p.velocity, [0.0, 0.0]);
        assert_eq!(p.force, [0.0, 0.0]);
        assert_eq!(p.ds, 1.0);
    }

    #[test]
    fn test_ibm_particle_displacement_at_rest() {
        let p = IbmParticle::new(1.0, 2.0);
        assert_eq!(p.displacement(), [0.0, 0.0]);
        assert_eq!(p.displacement_norm(), 0.0);
    }

    #[test]
    fn test_ibm_particle_displacement_displaced() {
        let p = IbmParticle::with_rest_and_ds(2.0, 3.0, 1.0, 1.0, 1.0);
        let d = p.displacement();
        assert!((d[0] - 1.0).abs() < 1e-14);
        assert!((d[1] - 2.0).abs() < 1e-14);
    }

    #[test]
    fn test_ibm_particle_advect() {
        let mut p = IbmParticle::new(0.0, 0.0);
        p.velocity = [1.0, 2.0];
        p.advect(0.5);
        assert!((p.position[0] - 0.5).abs() < 1e-14);
        assert!((p.position[1] - 1.0).abs() < 1e-14);
    }

    #[test]
    fn test_ibm_particle_apply_velocity_correction() {
        let mut p = IbmParticle::new(0.0, 0.0);
        p.velocity = [1.0, 0.0];
        p.apply_velocity_correction([0.5, -0.5]);
        assert!((p.velocity[0] - 1.5).abs() < 1e-14);
        assert!((p.velocity[1] + 0.5).abs() < 1e-14);
    }

    // ---- IbmMarker (legacy) ----

    #[test]
    fn test_ibm_marker_new() {
        let m = IbmMarker::new(3.0, 4.0);
        assert_eq!(m.position, [3.0, 4.0]);
        assert_eq!(m.velocity, [0.0, 0.0]);
        assert_eq!(m.force, [0.0, 0.0]);
    }

    #[test]
    fn test_ibm_marker_rest_position() {
        let m = IbmMarker::new(1.0, 2.0);
        assert_eq!(m.rest_position, [1.0, 2.0]);
    }

    #[test]
    fn test_ibm_marker_with_rest() {
        let m = IbmMarker::with_rest(1.0, 2.0, 0.5, 1.5);
        assert_eq!(m.position, [1.0, 2.0]);
        assert_eq!(m.rest_position, [0.5, 1.5]);
    }

    #[test]
    fn test_ibm_marker_clone() {
        let m = IbmMarker::new(7.0, 8.0);
        let m2 = m.clone();
        assert_eq!(m, m2);
    }

    // ---- compute_elastic_force ----

    #[test]
    fn test_elastic_force_at_rest() {
        let mut m = IbmMarker::new(2.0, 3.0);
        compute_elastic_force(std::slice::from_mut(&mut m), 10.0);
        assert!((m.force[0]).abs() < 1e-14);
        assert!((m.force[1]).abs() < 1e-14);
    }

    #[test]
    fn test_elastic_force_displaced() {
        let mut m = IbmMarker::with_rest(3.0, 4.0, 2.0, 3.0);
        compute_elastic_force(std::slice::from_mut(&mut m), 5.0);
        assert!((m.force[0] - (-5.0)).abs() < 1e-12);
        assert!((m.force[1] - (-5.0)).abs() < 1e-12);
    }

    #[test]
    fn test_elastic_force_negative_displacement() {
        let mut m = IbmMarker::with_rest(1.0, 2.0, 2.0, 4.0);
        compute_elastic_force(std::slice::from_mut(&mut m), 2.0);
        assert!((m.force[0] - 2.0).abs() < 1e-12);
        assert!((m.force[1] - 4.0).abs() < 1e-12);
    }

    #[test]
    fn test_elastic_force_zero_spring() {
        let mut m = IbmMarker::with_rest(5.0, 6.0, 0.0, 0.0);
        compute_elastic_force(std::slice::from_mut(&mut m), 0.0);
        assert_eq!(m.force, [0.0, 0.0]);
    }

    // ---- spread_force (legacy) ----

    #[test]
    fn test_spread_force_single_marker_at_node() {
        let nx = 5;
        let ny = 5;
        let mut m = IbmMarker::new(2.0, 2.0);
        m.force = [1.0, 0.0];
        let mut fx = vec![0.0_f64; nx * ny];
        let mut fy = vec![0.0_f64; nx * ny];
        spread_force(&[m], &mut fx, &mut fy, nx, ny, 1.0);
        let total_fx: f64 = fx.iter().sum();
        assert!(total_fx > 0.0, "no force was spread");
    }

    #[test]
    fn test_spread_force_zero_force() {
        let nx = 4;
        let ny = 4;
        let m = IbmMarker::new(2.0, 2.0);
        let mut fx = vec![0.0_f64; nx * ny];
        let mut fy = vec![0.0_f64; nx * ny];
        spread_force(&[m], &mut fx, &mut fy, nx, ny, 1.0);
        assert!(fx.iter().all(|&v| v == 0.0));
        assert!(fy.iter().all(|&v| v == 0.0));
    }

    #[test]
    fn test_spread_force_accumulates() {
        let nx = 5;
        let ny = 5;
        let mut m1 = IbmMarker::new(2.0, 2.0);
        m1.force = [1.0, 0.0];
        let mut m2 = IbmMarker::new(2.0, 2.0);
        m2.force = [1.0, 0.0];
        let mut fx = vec![0.0_f64; nx * ny];
        let mut fy = vec![0.0_f64; nx * ny];
        spread_force(&[m1], &mut fx, &mut fy, nx, ny, 1.0);
        let sum1: f64 = fx.iter().sum();
        let mut fx2 = vec![0.0_f64; nx * ny];
        let mut fy2 = vec![0.0_f64; nx * ny];
        spread_force(&[m2], &mut fx2, &mut fy2, nx, ny, 1.0);
        let sum2: f64 = fx2.iter().sum();
        assert!((sum1 - sum2).abs() < 1e-14);
    }

    #[test]
    fn test_spread_force_ds_scaling() {
        let nx = 6;
        let ny = 6;
        let mut m = IbmMarker::new(3.0, 3.0);
        m.force = [1.0, 1.0];
        let mut fx1 = vec![0.0_f64; nx * ny];
        let mut fy1 = vec![0.0_f64; nx * ny];
        spread_force(&[m.clone()], &mut fx1, &mut fy1, nx, ny, 1.0);
        let mut fx2 = vec![0.0_f64; nx * ny];
        let mut fy2 = vec![0.0_f64; nx * ny];
        spread_force(&[m], &mut fx2, &mut fy2, nx, ny, 2.0);
        let sum1: f64 = fx1.iter().sum();
        let sum2: f64 = fx2.iter().sum();
        assert!((sum2 / sum1 - 2.0).abs() < 1e-12);
    }

    // ---- interpolate_velocity (legacy) ----

    #[test]
    fn test_interpolate_velocity_uniform() {
        let nx = 5;
        let ny = 5;
        let ux = vec![1.0_f64; nx * ny];
        let uy = vec![2.0_f64; nx * ny];
        let mut m = IbmMarker::new(2.0, 2.0);
        interpolate_velocity(std::slice::from_mut(&mut m), &ux, &uy, nx, ny, 1.0);
        assert!(m.velocity[0] > 0.0, "vx = {}", m.velocity[0]);
        assert!(m.velocity[1] > 0.0, "vy = {}", m.velocity[1]);
        assert!(
            (m.velocity[1] / m.velocity[0] - 2.0).abs() < 1e-12,
            "ratio = {}",
            m.velocity[1] / m.velocity[0]
        );
    }

    #[test]
    fn test_interpolate_velocity_zero_field() {
        let nx = 4;
        let ny = 4;
        let ux = vec![0.0_f64; nx * ny];
        let uy = vec![0.0_f64; nx * ny];
        let mut m = IbmMarker::new(2.0, 2.0);
        interpolate_velocity(std::slice::from_mut(&mut m), &ux, &uy, nx, ny, 1.0);
        assert_eq!(m.velocity, [0.0, 0.0]);
    }

    #[test]
    fn test_interpolate_velocity_updates_marker() {
        let nx = 5;
        let ny = 5;
        let mut ux = vec![0.0_f64; nx * ny];
        let uy = vec![0.0_f64; nx * ny];
        ux[2 * nx + 2] = 5.0;
        let mut m = IbmMarker::new(2.0, 2.0);
        interpolate_velocity(std::slice::from_mut(&mut m), &ux, &uy, nx, ny, 1.0);
        assert!(m.velocity[0] > 0.0, "velocity not interpolated");
    }

    // ---- ib_lbm_step (legacy) ----

    #[test]
    fn test_ib_lbm_step_zeroes_grid_forces() {
        let nx = 5;
        let ny = 5;
        let mut markers = vec![IbmMarker::new(2.0, 2.0)];
        let mut fx = vec![9.9_f64; nx * ny];
        let mut fy = vec![9.9_f64; nx * ny];
        let ux = vec![0.0_f64; nx * ny];
        let uy = vec![0.0_f64; nx * ny];
        ib_lbm_step(
            &mut markers,
            &mut fx,
            &mut fy,
            &ux,
            &uy,
            nx,
            ny,
            1.0,
            0.0,
            0.1,
        );
        assert!(fx.iter().all(|&v| v == 0.0));
    }

    #[test]
    fn test_ib_lbm_step_marker_moves() {
        let nx = 7;
        let ny = 7;
        let mut markers = vec![IbmMarker::new(3.0, 3.0)];
        let mut fx = vec![0.0_f64; nx * ny];
        let mut fy = vec![0.0_f64; nx * ny];
        let ux = vec![0.1_f64; nx * ny];
        let uy = vec![0.0_f64; nx * ny];
        let pos0 = markers[0].position;
        ib_lbm_step(
            &mut markers,
            &mut fx,
            &mut fy,
            &ux,
            &uy,
            nx,
            ny,
            1.0,
            0.0,
            1.0,
        );
        assert!(markers[0].position[0] != pos0[0] || markers[0].position[1] != pos0[1]);
    }

    #[test]
    fn test_ib_lbm_step_spring_force_nonzero() {
        let nx = 7;
        let ny = 7;
        let mut markers = vec![IbmMarker::with_rest(4.0, 3.0, 3.0, 3.0)];
        let mut fx = vec![0.0_f64; nx * ny];
        let mut fy = vec![0.0_f64; nx * ny];
        let ux = vec![0.0_f64; nx * ny];
        let uy = vec![0.0_f64; nx * ny];
        ib_lbm_step(
            &mut markers,
            &mut fx,
            &mut fy,
            &ux,
            &uy,
            nx,
            ny,
            1.0,
            10.0,
            0.1,
        );
        let total: f64 = fx.iter().sum();
        assert!(total.abs() > 0.0, "spring force not spread");
    }

    #[test]
    fn test_ib_lbm_step_multiple_markers() {
        let nx = 8;
        let ny = 8;
        let mut markers = vec![IbmMarker::new(3.0, 3.0), IbmMarker::new(4.0, 4.0)];
        let mut fx = vec![0.0_f64; nx * ny];
        let mut fy = vec![0.0_f64; nx * ny];
        let ux = vec![0.0_f64; nx * ny];
        let uy = vec![0.0_f64; nx * ny];
        ib_lbm_step(
            &mut markers,
            &mut fx,
            &mut fy,
            &ux,
            &uy,
            nx,
            ny,
            1.0,
            1.0,
            0.01,
        );
        // Must not panic.
    }

    // ---- IbmFluidGrid ----

    #[test]
    fn test_ibm_fluid_grid_new() {
        let g = IbmFluidGrid::new(8, 6);
        assert_eq!(g.nx, 8);
        assert_eq!(g.ny, 6);
        assert_eq!(g.ncells(), 48);
    }

    #[test]
    fn test_ibm_fluid_grid_clear_forces() {
        let mut g = IbmFluidGrid::new(4, 4);
        g.fx = vec![1.0; 16];
        g.fy = vec![2.0; 16];
        g.clear_forces();
        assert!(g.fx.iter().all(|&v| v == 0.0));
        assert!(g.fy.iter().all(|&v| v == 0.0));
    }

    #[test]
    fn test_ibm_fluid_grid_set_uniform_velocity() {
        let mut g = IbmFluidGrid::new(4, 4);
        g.set_uniform_velocity(0.1, 0.05);
        assert!(g.ux.iter().all(|&v| (v - 0.1).abs() < 1e-14));
        assert!(g.uy.iter().all(|&v| (v - 0.05).abs() < 1e-14));
    }

    #[test]
    fn test_ibm_fluid_grid_spread_forces() {
        let mut g = IbmFluidGrid::new(8, 8);
        let mut p = IbmParticle::new(4.0, 4.0);
        p.force = [1.0, 0.5];
        g.spread_forces(&[p]);
        let total_fx: f64 = g.fx.iter().sum();
        assert!(total_fx > 0.0);
    }

    #[test]
    fn test_ibm_fluid_grid_interpolate_velocities() {
        let mut g = IbmFluidGrid::new(8, 8);
        g.set_uniform_velocity(1.0, 0.0);
        let mut particles = vec![IbmParticle::new(4.0, 4.0)];
        g.interpolate_velocities(&mut particles);
        assert!(
            particles[0].velocity[0] > 0.0,
            "interpolated vx = {}",
            particles[0].velocity[0]
        );
    }

    #[test]
    fn test_ibm_fluid_grid_total_force_zero() {
        let g = IbmFluidGrid::new(4, 4);
        let f = g.total_force();
        assert_eq!(f, [0.0, 0.0]);
    }

    #[test]
    fn test_ibm_fluid_grid_total_force_nonzero() {
        let mut g = IbmFluidGrid::new(4, 4);
        let mut p = IbmParticle::new(2.0, 2.0);
        p.force = [3.0, 2.0];
        g.spread_forces(&[p]);
        let f = g.total_force();
        assert!(f[0] > 0.0);
        assert!(f[1] > 0.0);
    }

    // ---- DirectForcingIbm ----

    #[test]
    fn test_direct_forcing_ibm_new() {
        let particles = vec![IbmParticle::new(4.0, 4.0)];
        let df = DirectForcingIbm::new(particles, 0.01, 1.0);
        assert_eq!(df.num_markers(), 1);
    }

    #[test]
    fn test_direct_forcing_ibm_set_target() {
        let particles = vec![IbmParticle::new(4.0, 4.0); 3];
        let df = DirectForcingIbm::new(particles, 0.01, 1.0);
        let targets = df.set_target_velocity([0.1, 0.0]);
        assert_eq!(targets.len(), 3);
        assert!(targets.iter().all(|&t| t == [0.1, 0.0]));
    }

    #[test]
    fn test_direct_forcing_ibm_compute_body_forces() {
        let particles = vec![IbmParticle::new(4.0, 4.0)];
        let mut df = DirectForcingIbm::new(particles, 1.0, 1.0);
        let u_target = vec![[0.1_f64, 0.0_f64]];
        let u_fluid = vec![[0.0_f64, 0.0_f64]];
        df.compute_body_forces(&u_target, &u_fluid);
        // f = rho * (u_target - u_fluid) / dt = 1 * 0.1 / 1 = 0.1
        assert!((df.particles[0].force[0] - 0.1).abs() < 1e-12);
        assert!((df.particles[0].force[1]).abs() < 1e-14);
    }

    #[test]
    fn test_direct_forcing_ibm_total_force() {
        let mut particles = vec![IbmParticle::new(4.0, 4.0)];
        particles[0].force = [2.0, -1.0];
        let df = DirectForcingIbm::new(particles, 1.0, 1.0);
        let f = df.total_lagrangian_force();
        assert!((f[0] - 2.0).abs() < 1e-14);
        assert!((f[1] + 1.0).abs() < 1e-14);
    }

    #[test]
    fn test_direct_forcing_ibm_step() {
        let particles = vec![IbmParticle::new(4.0, 4.0)];
        let mut df = DirectForcingIbm::new(particles, 0.01, 1.0);
        let mut grid = IbmFluidGrid::new(10, 10);
        grid.set_uniform_velocity(0.1, 0.0);
        let u_target = df.set_target_velocity([0.0, 0.0]);
        df.step(&mut grid, &u_target);
        // Forces should be non-zero since u_fluid != u_target.
        assert!(df.particles[0].force[0].abs() > 0.0);
    }

    // ---- MovingIbmObject ----

    #[test]
    fn test_moving_ibm_object_circle_construction() {
        let obj = MovingIbmObject::circle(5.0, 5.0, 2.0, 16, 1.0, 0.5);
        assert_eq!(obj.num_markers(), 16);
        assert_eq!(obj.center, [5.0, 5.0]);
        assert_eq!(obj.velocity, [0.0, 0.0]);
    }

    #[test]
    fn test_moving_ibm_object_surface_velocity_no_rotation() {
        let obj = MovingIbmObject::circle(5.0, 5.0, 2.0, 8, 1.0, 0.5);
        // With zero omega and zero translation, surface velocity should be zero.
        for i in 0..obj.num_markers() {
            let v = obj.surface_velocity(i);
            assert!((v[0]).abs() < 1e-14, "v[0] = {}", v[0]);
            assert!((v[1]).abs() < 1e-14, "v[1] = {}", v[1]);
        }
    }

    #[test]
    fn test_moving_ibm_object_kinetic_energy_zero() {
        let obj = MovingIbmObject::circle(5.0, 5.0, 2.0, 8, 1.0, 0.5);
        assert_eq!(obj.kinetic_energy(), 0.0);
    }

    #[test]
    fn test_moving_ibm_object_kinetic_energy_translation() {
        let mut obj = MovingIbmObject::circle(5.0, 5.0, 2.0, 8, 2.0, 1.0);
        obj.velocity = [3.0, 4.0];
        let ke = obj.kinetic_energy();
        // 0.5 * 2 * (9 + 16) = 25
        assert!((ke - 25.0).abs() < 1e-10);
    }

    #[test]
    fn test_moving_ibm_object_hydrodynamic_force_zero_markers() {
        let obj = MovingIbmObject::circle(5.0, 5.0, 2.0, 8, 1.0, 0.5);
        // All forces are zero → net should be zero.
        let (f, tau) = obj.compute_hydrodynamic_force();
        assert_eq!(f, [0.0, 0.0]);
        assert_eq!(tau, 0.0);
    }

    #[test]
    fn test_moving_ibm_object_update_velocity_no_force() {
        let mut obj = MovingIbmObject::circle(5.0, 5.0, 2.0, 8, 1.0, 0.5);
        obj.update_velocity(0.1);
        // No forces → velocity stays zero.
        assert_eq!(obj.velocity, [0.0, 0.0]);
        assert_eq!(obj.omega, 0.0);
    }

    // ---- FlexibleIbmFiber ----

    #[test]
    fn test_flexible_fiber_straight_construction() {
        let fiber = FlexibleIbmFiber::straight(0.0, 5.0, 10, 1.0, 1000.0, 0.1);
        assert_eq!(fiber.num_markers(), 10);
        assert!((fiber.arc_length() - 9.0).abs() < 1e-10);
    }

    #[test]
    fn test_flexible_fiber_elastic_energy_at_rest_zero() {
        let fiber = FlexibleIbmFiber::straight(0.0, 5.0, 5, 1.0, 1000.0, 0.1);
        let e = fiber.elastic_energy();
        assert!(e.abs() < 1e-10, "energy at rest = {e}");
    }

    #[test]
    fn test_flexible_fiber_tip_deflection_at_rest() {
        let fiber = FlexibleIbmFiber::straight(0.0, 5.0, 5, 1.0, 100.0, 0.1);
        assert!((fiber.tip_deflection()).abs() < 1e-14);
    }

    #[test]
    fn test_flexible_fiber_tension_forces_at_rest() {
        let mut fiber = FlexibleIbmFiber::straight(0.0, 5.0, 5, 1.0, 100.0, 0.1);
        fiber.compute_tension_forces();
        let total_fx: f64 = fiber.particles.iter().map(|p| p.force[0]).sum();
        assert!(total_fx.abs() < 1e-10, "net force = {total_fx}");
    }

    #[test]
    fn test_flexible_fiber_bending_forces_straight_zero() {
        let mut fiber = FlexibleIbmFiber::straight(0.0, 5.0, 5, 1.0, 100.0, 1.0);
        fiber.compute_tension_forces(); // resets forces first
        fiber.compute_bending_forces();
        // Straight fiber: d2x = 0, d2y = 0 → bending forces should be zero.
        let total_fy: f64 = fiber.particles.iter().map(|p| p.force[1]).sum();
        assert!(total_fy.abs() < 1e-10, "net bending force y = {total_fy}");
    }

    #[test]
    fn test_flexible_fiber_arc_length_uniform_spacing() {
        let n = 6;
        let ds = 1.5;
        let fiber = FlexibleIbmFiber::straight(0.0, 0.0, n, ds, 100.0, 0.1);
        let expected = ds * (n - 1) as f64;
        assert!((fiber.arc_length() - expected).abs() < 1e-10);
    }

    // ---- IbmCouplingStats ----

    #[test]
    fn test_ibm_stats_new() {
        let s = IbmCouplingStats::new(1.0, 1.0, 1.0);
        assert_eq!(s.n_samples, 0);
        assert_eq!(s.mean_drag(), 0.0);
        assert_eq!(s.mean_lift(), 0.0);
    }

    #[test]
    fn test_ibm_stats_record_and_mean() {
        let mut s = IbmCouplingStats::new(1.0, 1.0, 1.0);
        s.record(2.0, 1.0, 0.5, [0.0, 0.0]);
        s.record(4.0, 3.0, 1.5, [0.0, 0.0]);
        assert!((s.mean_drag() - 3.0).abs() < 1e-12);
        assert!((s.mean_lift() - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_ibm_stats_drag_coefficient() {
        let mut s = IbmCouplingStats::new(1.0, 2.0, 1.0);
        // mean_drag = 2, q = 0.5 * 1 * 1 * 2 = 1, CD = 2
        s.record(2.0, 0.0, 0.0, [0.0, 0.0]);
        let cd = s.drag_coefficient();
        assert!((cd - 2.0).abs() < 1e-10, "CD = {cd}");
    }

    #[test]
    fn test_ibm_stats_lift_coefficient() {
        let mut s = IbmCouplingStats::new(1.0, 2.0, 1.0);
        s.record(0.0, 1.0, 0.0, [0.0, 0.0]);
        let cl = s.lift_coefficient();
        assert!((cl - 1.0).abs() < 1e-10, "CL = {cl}");
    }

    #[test]
    fn test_ibm_stats_reset() {
        let mut s = IbmCouplingStats::new(1.0, 1.0, 1.0);
        s.record(5.0, 3.0, 1.0, [1.0, 0.0]);
        s.reset();
        assert_eq!(s.n_samples, 0);
        assert_eq!(s.drag, 0.0);
    }

    #[test]
    fn test_ibm_stats_strouhal_number() {
        let s = IbmCouplingStats::new(2.0, 1.0, 1.0);
        let st = s.strouhal_number(0.4);
        assert!((st - 0.2).abs() < 1e-10, "St = {st}");
    }

    #[test]
    fn test_ibm_stats_power() {
        let mut s = IbmCouplingStats::new(1.0, 1.0, 1.0);
        s.record(2.0, 0.0, 0.0, [0.5, 0.0]);
        assert!((s.mean_power() - 1.0).abs() < 1e-12);
    }

    // ---- Utility helpers ----

    #[test]
    fn test_make_circle_markers_count() {
        let markers = make_circle_markers(5.0, 5.0, 3.0, 12);
        assert_eq!(markers.len(), 12);
    }

    #[test]
    fn test_make_circle_markers_radius() {
        let cx = 5.0;
        let cy = 5.0;
        let r = 2.0;
        let markers = make_circle_markers(cx, cy, r, 8);
        for p in &markers {
            let dist = ((p.position[0] - cx).powi(2) + (p.position[1] - cy).powi(2)).sqrt();
            assert!((dist - r).abs() < 1e-10, "dist = {dist}");
        }
    }

    #[test]
    fn test_make_circle_markers_ds() {
        let r = 3.0;
        let n = 12;
        let expected_ds = 2.0 * PI * r / n as f64;
        let markers = make_circle_markers(0.0, 0.0, r, n);
        for p in &markers {
            assert!((p.ds - expected_ds).abs() < 1e-10);
        }
    }

    #[test]
    fn test_make_line_markers_count() {
        let markers = make_line_markers(0.0, 0.0, 5.0, 0.0, 6);
        assert_eq!(markers.len(), 6);
    }

    #[test]
    fn test_make_line_markers_endpoints() {
        let markers = make_line_markers(1.0, 2.0, 4.0, 2.0, 4);
        assert!((markers[0].position[0] - 1.0).abs() < 1e-10);
        assert!((markers[3].position[0] - 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_make_line_markers_single() {
        let markers = make_line_markers(1.0, 2.0, 4.0, 5.0, 1);
        assert_eq!(markers.len(), 1);
    }

    #[test]
    fn test_local_curvature_straight() {
        let markers = make_line_markers(0.0, 0.0, 5.0, 0.0, 6);
        for i in 1..5 {
            let kappa = local_curvature(&markers, i);
            assert!(kappa.abs() < 1e-10, "kappa[{i}] = {kappa}");
        }
    }

    #[test]
    fn test_local_curvature_boundary() {
        let markers = make_line_markers(0.0, 0.0, 5.0, 0.0, 6);
        assert_eq!(local_curvature(&markers, 0), 0.0);
        assert_eq!(local_curvature(&markers, 5), 0.0);
    }

    #[test]
    fn test_penalisation_factor_zero_nu() {
        let f = penalisation_factor(1.0, 0.0, 0.1);
        // eta / (eta + 0) = 1.
        assert!((f - 1.0).abs() < 1e-14);
    }

    #[test]
    fn test_penalisation_factor_large_nu() {
        let f = penalisation_factor(0.001, 10.0, 1.0);
        // eta / (eta + nu*dt) = 0.001 / 10.001 ≈ 0.0001
        assert!(f < 0.01, "f = {f}");
    }

    #[test]
    fn test_ibm_fluid_grid_spread_forces_peskin4() {
        let mut g = IbmFluidGrid::new(10, 10);
        let mut p = IbmParticle::new(5.0, 5.0);
        p.force = [1.0, 1.0];
        g.spread_forces_peskin4(&[p]);
        let total_fx: f64 = g.fx.iter().sum();
        assert!(total_fx > 0.0);
    }
}
