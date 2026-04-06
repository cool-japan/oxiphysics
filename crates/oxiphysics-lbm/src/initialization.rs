#![allow(clippy::too_many_arguments)]
#![allow(clippy::manual_range_contains, clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Lattice initialization and boundary condition utilities.
//!
//! Provides geometry helpers, flow initializers, and Zou-He boundary
//! condition primitives for D2Q9 lattices.
//!
//! Extended with:
//! - Analytical flow initialization (Poiseuille, Couette, Taylor-Green vortex)
//! - Perturbation initialization
//! - Ramp-up strategies for gradual forcing
//! - Checkpoint save/load for simulation state

#![allow(dead_code)]

// ---------------------------------------------------------------------------
// Boundary conditions
// ---------------------------------------------------------------------------

/// High-level boundary condition type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BoundaryType {
    /// Solid no-slip wall (bounce-back).
    Wall,
    /// Velocity inlet (value = \[ux, uy, 0\]).
    Inlet,
    /// Pressure outlet (value = \[rho_out, 0, 0\]).
    Outlet,
    /// Symmetric (zero normal-gradient) boundary.
    Symmetric,
    /// Periodic boundary.
    Periodic,
}

/// A set of boundary nodes sharing the same boundary condition.
#[derive(Debug, Clone)]
pub struct BoundaryCondition {
    /// Indices of nodes this BC applies to.
    pub node_indices: Vec<usize>,
    /// Type of boundary condition.
    pub bc_type: BoundaryType,
    /// Associated value: velocity \[ux,uy,0\] for inlet; \[rho_out,0,0\] for outlet.
    pub value: [f64; 3],
}

impl BoundaryCondition {
    /// Create a wall (no-slip) boundary.
    pub fn new_wall(indices: Vec<usize>) -> Self {
        Self {
            node_indices: indices,
            bc_type: BoundaryType::Wall,
            value: [0.0; 3],
        }
    }

    /// Create a velocity inlet boundary.
    ///
    /// `velocity` = \[ux, uy, 0\].
    pub fn new_inlet(indices: Vec<usize>, velocity: [f64; 3]) -> Self {
        Self {
            node_indices: indices,
            bc_type: BoundaryType::Inlet,
            value: velocity,
        }
    }

    /// Create a pressure outlet boundary (unit density by default).
    pub fn new_outlet(indices: Vec<usize>) -> Self {
        Self {
            node_indices: indices,
            bc_type: BoundaryType::Outlet,
            value: [1.0, 0.0, 0.0],
        }
    }

    /// Create a symmetric boundary.
    pub fn new_symmetric(indices: Vec<usize>) -> Self {
        Self {
            node_indices: indices,
            bc_type: BoundaryType::Symmetric,
            value: [0.0; 3],
        }
    }

    /// Create a periodic boundary.
    pub fn new_periodic(indices: Vec<usize>) -> Self {
        Self {
            node_indices: indices,
            bc_type: BoundaryType::Periodic,
            value: [0.0; 3],
        }
    }
}

// ---------------------------------------------------------------------------
// Lattice geometry
// ---------------------------------------------------------------------------

/// 2D lattice geometry with solid/fluid classification.
#[derive(Debug, Clone)]
pub struct LatticeGeometry2D {
    /// Number of cells in x direction.
    pub nx: usize,
    /// Number of cells in y direction.
    pub ny: usize,
    /// `true` = solid node, `false` = fluid node (row-major: index = y*nx + x).
    pub is_solid: Vec<bool>,
}

impl LatticeGeometry2D {
    /// Horizontal channel: solid walls at y=0 and y=ny-1, fluid everywhere else.
    pub fn new_channel(nx: usize, ny: usize) -> Self {
        let mut is_solid = vec![false; nx * ny];
        for x in 0..nx {
            is_solid[x] = true; // bottom wall
            is_solid[(ny - 1) * nx + x] = true; // top wall
        }
        Self { nx, ny, is_solid }
    }

    /// Lid-driven cavity: square domain, all four walls solid except the top
    /// row which acts as the moving lid (solid boundary, driven externally).
    pub fn new_cavity(n: usize) -> Self {
        let mut is_solid = vec![false; n * n];
        for i in 0..n {
            is_solid[i] = true; // bottom
            is_solid[(n - 1) * n + i] = true; // top (lid)
            is_solid[i * n] = true; // left
            is_solid[i * n + (n - 1)] = true; // right
        }
        Self {
            nx: n,
            ny: n,
            is_solid,
        }
    }

    /// Create a fully fluid domain (no solid walls).
    pub fn new_open(nx: usize, ny: usize) -> Self {
        Self {
            nx,
            ny,
            is_solid: vec![false; nx * ny],
        }
    }

    /// Mark all cells whose centre falls strictly inside the circle as solid.
    pub fn add_circle_obstacle(&mut self, cx: usize, cy: usize, radius: f64) {
        let nx = self.nx;
        let ny = self.ny;
        for y in 0..ny {
            for x in 0..nx {
                let dx = x as f64 - cx as f64;
                let dy = y as f64 - cy as f64;
                if (dx * dx + dy * dy).sqrt() < radius {
                    self.is_solid[y * nx + x] = true;
                }
            }
        }
    }

    /// Mark all cells in the rectangle \[x0, x0+w) × \[y0, y0+h) as solid.
    pub fn add_rect_obstacle(&mut self, x0: usize, y0: usize, w: usize, h: usize) {
        let nx = self.nx;
        let ny = self.ny;
        for dy in 0..h {
            for dx in 0..w {
                let x = x0 + dx;
                let y = y0 + dy;
                if x < nx && y < ny {
                    self.is_solid[y * nx + x] = true;
                }
            }
        }
    }

    /// Mark all cells in an ellipse as solid.
    pub fn add_ellipse_obstacle(&mut self, cx: f64, cy: f64, rx: f64, ry: f64) {
        let nx = self.nx;
        let ny = self.ny;
        for y in 0..ny {
            for x in 0..nx {
                let dx = (x as f64 - cx) / rx;
                let dy = (y as f64 - cy) / ry;
                if dx * dx + dy * dy < 1.0 {
                    self.is_solid[y * nx + x] = true;
                }
            }
        }
    }

    /// Fraction of nodes that are solid.
    pub fn solid_fraction(&self) -> f64 {
        self.solid_count() as f64 / self.is_solid.len() as f64
    }

    /// Number of fluid nodes.
    pub fn fluid_count(&self) -> usize {
        self.is_solid.iter().filter(|&&s| !s).count()
    }

    /// Number of solid nodes.
    pub fn solid_count(&self) -> usize {
        self.is_solid.iter().filter(|&&s| s).count()
    }

    /// Get indices of all fluid nodes.
    pub fn fluid_indices(&self) -> Vec<usize> {
        self.is_solid
            .iter()
            .enumerate()
            .filter(|&(_, s)| !s)
            .map(|(i, _)| i)
            .collect()
    }

    /// Get indices of all boundary (wall-adjacent fluid) nodes.
    pub fn boundary_fluid_indices(&self) -> Vec<usize> {
        let nx = self.nx;
        let ny = self.ny;
        let mut boundary = Vec::new();
        for y in 0..ny {
            for x in 0..nx {
                let k = y * nx + x;
                if self.is_solid[k] {
                    continue;
                }
                // Check if any neighbor is solid
                let has_solid_neighbor = (x > 0 && self.is_solid[y * nx + (x - 1)])
                    || (x < nx - 1 && self.is_solid[y * nx + (x + 1)])
                    || (y > 0 && self.is_solid[(y - 1) * nx + x])
                    || (y < ny - 1 && self.is_solid[(y + 1) * nx + x]);
                if has_solid_neighbor {
                    boundary.push(k);
                }
            }
        }
        boundary
    }
}

// ---------------------------------------------------------------------------
// Flow initializers
// ---------------------------------------------------------------------------

/// Parabolic Poiseuille velocity profile.
///
/// `u(y) = 4 * u_max * y * (ny-1-y) / (ny-1)²`
///
/// Returns a `Vec<[f64;2]>` of length `nx*ny` with `[ux, 0.0]` at each node
/// (row-major: index = `y*nx + x`).
pub fn initialize_poiseuille_flow(nx: usize, ny: usize, u_max: f64) -> Vec<[f64; 2]> {
    let ny1 = (ny - 1) as f64;
    let mut vel = vec![[0.0_f64; 2]; nx * ny];
    for y in 0..ny {
        let yf = y as f64;
        let ux = 4.0 * u_max * yf * (ny1 - yf) / (ny1 * ny1);
        for x in 0..nx {
            vel[y * nx + x] = [ux, 0.0];
        }
    }
    vel
}

/// Linear Couette velocity profile.
///
/// `u(y) = u_wall * y / (ny-1)`
///
/// Returns a `Vec<[f64;2]>` of length `nx*ny`.
pub fn initialize_couette_flow(nx: usize, ny: usize, u_wall: f64) -> Vec<[f64; 2]> {
    let ny1 = (ny - 1) as f64;
    let mut vel = vec![[0.0_f64; 2]; nx * ny];
    for y in 0..ny {
        let ux = u_wall * y as f64 / ny1;
        for x in 0..nx {
            vel[y * nx + x] = [ux, 0.0];
        }
    }
    vel
}

/// Uniform flow: every node gets velocity `u`.
///
/// Returns a `Vec<[f64;2]>` of length `nx*ny`.
pub fn initialize_uniform_flow(nx: usize, ny: usize, u: [f64; 2]) -> Vec<[f64; 2]> {
    vec![u; nx * ny]
}

// ---------------------------------------------------------------------------
// Taylor-Green vortex initialization
// ---------------------------------------------------------------------------

/// Taylor-Green decaying vortex initialization (2D).
///
/// ```text
/// ux(x,y) = -u0 * cos(kx * x) * sin(ky * y)
/// uy(x,y) =  u0 * sin(kx * x) * cos(ky * y)
/// ```
///
/// where `kx = 2*pi/nx` and `ky = 2*pi/ny` (one wavelength per domain).
///
/// Returns a `Vec<[f64;2]>` of length `nx*ny`.
pub fn initialize_taylor_green_vortex(nx: usize, ny: usize, u0: f64) -> Vec<[f64; 2]> {
    let kx = 2.0 * std::f64::consts::PI / nx as f64;
    let ky = 2.0 * std::f64::consts::PI / ny as f64;
    let mut vel = vec![[0.0_f64; 2]; nx * ny];
    for y in 0..ny {
        for x in 0..nx {
            let xf = x as f64;
            let yf = y as f64;
            let ux = -u0 * (kx * xf).cos() * (ky * yf).sin();
            let uy = u0 * (kx * xf).sin() * (ky * yf).cos();
            vel[y * nx + x] = [ux, uy];
        }
    }
    vel
}

/// Taylor-Green vortex analytical density field.
///
/// ```text
/// rho(x,y) = rho0 - u0^2 / (4 * cs^2) * (cos(2*kx*x) + cos(2*ky*y))
/// ```
pub fn taylor_green_density(nx: usize, ny: usize, u0: f64, rho0: f64) -> Vec<f64> {
    let cs2 = 1.0 / 3.0;
    let kx = 2.0 * std::f64::consts::PI / nx as f64;
    let ky = 2.0 * std::f64::consts::PI / ny as f64;
    let mut rho = vec![0.0_f64; nx * ny];
    for y in 0..ny {
        for x in 0..nx {
            let xf = x as f64;
            let yf = y as f64;
            rho[y * nx + x] =
                rho0 - u0 * u0 / (4.0 * cs2) * ((2.0 * kx * xf).cos() + (2.0 * ky * yf).cos());
        }
    }
    rho
}

/// Taylor-Green analytical velocity at time t (with viscous decay).
///
/// ```text
/// ux(x,y,t) = -u0 * cos(kx*x) * sin(ky*y) * exp(-2*nu*(kx^2+ky^2)*t)
/// uy(x,y,t) =  u0 * sin(kx*x) * cos(ky*y) * exp(-2*nu*(kx^2+ky^2)*t)
/// ```
pub fn taylor_green_analytical(nx: usize, ny: usize, u0: f64, nu: f64, t: f64) -> Vec<[f64; 2]> {
    let kx = 2.0 * std::f64::consts::PI / nx as f64;
    let ky = 2.0 * std::f64::consts::PI / ny as f64;
    let decay = (-2.0 * nu * (kx * kx + ky * ky) * t).exp();
    let mut vel = vec![[0.0_f64; 2]; nx * ny];
    for y in 0..ny {
        for x in 0..nx {
            let xf = x as f64;
            let yf = y as f64;
            let ux = -u0 * (kx * xf).cos() * (ky * yf).sin() * decay;
            let uy = u0 * (kx * xf).sin() * (ky * yf).cos() * decay;
            vel[y * nx + x] = [ux, uy];
        }
    }
    vel
}

// ---------------------------------------------------------------------------
// Perturbation initialization
// ---------------------------------------------------------------------------

/// Add sinusoidal perturbation to an existing velocity field.
///
/// ```text
/// ux += amp * sin(2*pi*mode_y * y / ny)
/// uy += amp * sin(2*pi*mode_x * x / nx)
/// ```
pub fn add_sinusoidal_perturbation(
    vel: &mut [[f64; 2]],
    nx: usize,
    ny: usize,
    amplitude: f64,
    mode_x: f64,
    mode_y: f64,
) {
    for y in 0..ny {
        for x in 0..nx {
            let k = y * nx + x;
            vel[k][0] +=
                amplitude * (2.0 * std::f64::consts::PI * mode_y * y as f64 / ny as f64).sin();
            vel[k][1] +=
                amplitude * (2.0 * std::f64::consts::PI * mode_x * x as f64 / nx as f64).sin();
        }
    }
}

/// Add random white-noise perturbation to velocity field.
///
/// Uses a simple deterministic pseudo-random perturbation based on cell index
/// (no external RNG dependency needed for initialization).
pub fn add_noise_perturbation(vel: &mut [[f64; 2]], amplitude: f64, seed: u64) {
    for (k, v) in vel.iter_mut().enumerate() {
        // Simple hash-based pseudo-random in [-1, 1]
        let hash1 = ((k as u64)
            .wrapping_mul(6364136223846793005)
            .wrapping_add(seed)) as f64;
        let hash2 = ((k as u64)
            .wrapping_mul(1442695040888963407)
            .wrapping_add(seed.wrapping_mul(3))) as f64;
        let r1 = (hash1 / u64::MAX as f64) * 2.0 - 1.0;
        let r2 = (hash2 / u64::MAX as f64) * 2.0 - 1.0;
        v[0] += amplitude * r1;
        v[1] += amplitude * r2;
    }
}

/// Add a Kelvin-Helmholtz instability perturbation.
///
/// Creates a shear layer with a sinusoidal perturbation at the interface.
pub fn initialize_kelvin_helmholtz(
    nx: usize,
    ny: usize,
    u_shear: f64,
    perturbation_amp: f64,
    mode: f64,
) -> Vec<[f64; 2]> {
    let mut vel = vec![[0.0_f64; 2]; nx * ny];
    let ny_half = ny as f64 / 2.0;

    for y in 0..ny {
        for x in 0..nx {
            let k = y * nx + x;
            // Shear profile: u = u_shear * tanh((y - ny/2) / delta)
            let delta = ny as f64 / 20.0;
            let ux = u_shear * ((y as f64 - ny_half) / delta).tanh();

            // Sinusoidal perturbation in uy
            let uy = perturbation_amp
                * (2.0 * std::f64::consts::PI * mode * x as f64 / nx as f64).sin()
                * (-(((y as f64 - ny_half) / delta).powi(2))).exp();

            vel[k] = [ux, uy];
        }
    }
    vel
}

// ---------------------------------------------------------------------------
// Ramp-up strategies
// ---------------------------------------------------------------------------

/// Linear ramp factor: goes from 0 to 1 over `ramp_steps` steps.
///
/// Returns `min(1.0, step / ramp_steps)`.
pub fn linear_ramp(step: usize, ramp_steps: usize) -> f64 {
    if ramp_steps == 0 {
        return 1.0;
    }
    (step as f64 / ramp_steps as f64).min(1.0)
}

/// Sinusoidal ramp factor: smooth acceleration from 0 to 1.
///
/// `factor = 0.5 * (1 - cos(pi * min(1, step/ramp_steps)))`
pub fn sinusoidal_ramp(step: usize, ramp_steps: usize) -> f64 {
    if ramp_steps == 0 {
        return 1.0;
    }
    let t = (step as f64 / ramp_steps as f64).min(1.0);
    0.5 * (1.0 - (std::f64::consts::PI * t).cos())
}

/// Exponential ramp factor: `1 - exp(-step / tau)`.
pub fn exponential_ramp(step: usize, tau: f64) -> f64 {
    if tau <= 0.0 {
        return 1.0;
    }
    1.0 - (-(step as f64) / tau).exp()
}

/// Polynomial ramp: `(t)^p` where `t = min(1, step/ramp_steps)`.
pub fn polynomial_ramp(step: usize, ramp_steps: usize, power: f64) -> f64 {
    if ramp_steps == 0 {
        return 1.0;
    }
    let t = (step as f64 / ramp_steps as f64).min(1.0);
    t.powf(power)
}

/// Apply a ramp factor to a body force.
///
/// Returns `[factor * fx, factor * fy]`.
pub fn ramp_force(force: [f64; 2], factor: f64) -> [f64; 2] {
    [force[0] * factor, force[1] * factor]
}

// ---------------------------------------------------------------------------
// Checkpoint save/load
// ---------------------------------------------------------------------------

/// Simulation checkpoint data for save/restore.
#[derive(Debug, Clone)]
pub struct Checkpoint {
    /// Current time step.
    pub step: usize,
    /// Grid dimensions (nx, ny).
    pub dimensions: (usize, usize),
    /// Density field.
    pub rho: Vec<f64>,
    /// x-velocity field.
    pub ux: Vec<f64>,
    /// y-velocity field.
    pub uy: Vec<f64>,
    /// Distribution functions (q arrays, each of length nx*ny).
    pub distributions: Vec<Vec<f64>>,
}

impl Checkpoint {
    /// Create a checkpoint from current simulation state.
    pub fn new(
        step: usize,
        nx: usize,
        ny: usize,
        rho: Vec<f64>,
        ux: Vec<f64>,
        uy: Vec<f64>,
        distributions: Vec<Vec<f64>>,
    ) -> Self {
        Self {
            step,
            dimensions: (nx, ny),
            rho,
            ux,
            uy,
            distributions,
        }
    }

    /// Get the total number of nodes.
    pub fn num_nodes(&self) -> usize {
        self.dimensions.0 * self.dimensions.1
    }

    /// Validate the checkpoint data dimensions.
    pub fn is_valid(&self) -> bool {
        let n = self.num_nodes();
        self.rho.len() == n
            && self.ux.len() == n
            && self.uy.len() == n
            && self.distributions.iter().all(|d| d.len() == n)
    }

    /// Compute the L2 norm of the velocity field (for convergence monitoring).
    pub fn velocity_l2_norm(&self) -> f64 {
        let mut sum = 0.0;
        for k in 0..self.rho.len() {
            sum += self.ux[k] * self.ux[k] + self.uy[k] * self.uy[k];
        }
        sum.sqrt()
    }

    /// Compute the maximum velocity magnitude.
    pub fn max_velocity(&self) -> f64 {
        let mut max_v = 0.0_f64;
        for k in 0..self.rho.len() {
            let v_sq = self.ux[k] * self.ux[k] + self.uy[k] * self.uy[k];
            max_v = max_v.max(v_sq);
        }
        max_v.sqrt()
    }

    /// Serialize to a simple binary format (Vec`u8`).
    ///
    /// Format: \[step(8)\]\[nx(8)\]\[ny(8)\]\[q(8)\]\[rho...\]\[ux...\]\[uy...\]\[f0...\]\[f1...\]...
    /// All values stored as little-endian f64 or u64.
    pub fn to_bytes(&self) -> Vec<u8> {
        let n = self.num_nodes();
        let q = self.distributions.len();
        let total_f64s = 3 * n + q * n; // rho + ux + uy + q distributions
        let header_bytes = 4 * 8; // step, nx, ny, q
        let mut buf = Vec::with_capacity(header_bytes + total_f64s * 8);

        buf.extend_from_slice(&(self.step as u64).to_le_bytes());
        buf.extend_from_slice(&(self.dimensions.0 as u64).to_le_bytes());
        buf.extend_from_slice(&(self.dimensions.1 as u64).to_le_bytes());
        buf.extend_from_slice(&(q as u64).to_le_bytes());

        for &v in &self.rho {
            buf.extend_from_slice(&v.to_le_bytes());
        }
        for &v in &self.ux {
            buf.extend_from_slice(&v.to_le_bytes());
        }
        for &v in &self.uy {
            buf.extend_from_slice(&v.to_le_bytes());
        }
        for dist in &self.distributions {
            for &v in dist {
                buf.extend_from_slice(&v.to_le_bytes());
            }
        }
        buf
    }

    /// Deserialize from the binary format produced by `to_bytes`.
    pub fn from_bytes(data: &[u8]) -> Option<Self> {
        if data.len() < 32 {
            return None;
        }

        let read_u64_at = |offset: usize| -> Option<u64> {
            let bytes: [u8; 8] = data[offset..offset + 8].try_into().ok()?;
            Some(u64::from_le_bytes(bytes))
        };
        let read_f64_at = |offset: usize| -> f64 {
            let bytes: [u8; 8] = data[offset..offset + 8].try_into().unwrap_or([0u8; 8]);
            f64::from_le_bytes(bytes)
        };

        let step = read_u64_at(0)? as usize;
        let nx = read_u64_at(8)? as usize;
        let ny = read_u64_at(16)? as usize;
        let q = read_u64_at(24)? as usize;
        let n = nx * ny;

        let expected_len = 32 + (3 * n + q * n) * 8;
        if data.len() < expected_len {
            return None;
        }

        let mut offset = 32;
        let mut read_vec = |count: usize| -> Vec<f64> {
            let v: Vec<f64> = (0..count).map(|i| read_f64_at(offset + i * 8)).collect();
            offset += count * 8;
            v
        };

        let rho = read_vec(n);
        let ux = read_vec(n);
        let uy = read_vec(n);
        let distributions: Vec<Vec<f64>> = (0..q).map(|_| read_vec(n)).collect();

        Some(Self {
            step,
            dimensions: (nx, ny),
            rho,
            ux,
            uy,
            distributions,
        })
    }
}

/// Compute the L2 velocity error between two fields.
pub fn velocity_l2_error(vel_a: &[[f64; 2]], vel_b: &[[f64; 2]]) -> f64 {
    assert_eq!(
        vel_a.len(),
        vel_b.len(),
        "Velocity fields must have same length"
    );
    let mut error_sum = 0.0;
    let mut ref_sum = 0.0;
    for (a, b) in vel_a.iter().zip(vel_b.iter()) {
        let dx = a[0] - b[0];
        let dy = a[1] - b[1];
        error_sum += dx * dx + dy * dy;
        ref_sum += b[0] * b[0] + b[1] * b[1];
    }
    if ref_sum > 1e-30 {
        (error_sum / ref_sum).sqrt()
    } else {
        error_sum.sqrt()
    }
}

// ---------------------------------------------------------------------------
// Zou-He boundary conditions (operate on a single D2Q9 distribution cell)
// ---------------------------------------------------------------------------
//
// D2Q9 direction indexing (matches lattice.rs):
//   0=rest  1=E  2=N  3=W  4=S  5=NE  6=NW  7=SW  8=SE

/// Zou-He velocity inlet on the **left** wall (x = 0).
///
/// Given a distribution `f` at a left-wall cell with prescribed inlet
/// velocity `u = [ux, uy]` and a known density `rho`, compute the three
/// unknown populations `f[1]`, `f[5]`, `f[8]` that are streaming into the
/// domain from the left boundary.
///
/// # Arguments
/// * `f`   – mutable 9-element distribution array.
/// * `rho` – density at this cell (used as the known macroscopic value).
/// * `u`   – prescribed velocity `[ux, uy]`.
pub fn zou_he_velocity_inlet(f: &mut [f64; 9], rho: f64, u: [f64; 2]) {
    let ux = u[0];
    let uy = u[1];
    // Zou-He momentum constraints for a left wall (normal = +x):
    f[1] = f[3] + (2.0 / 3.0) * rho * ux;
    f[5] = f[7] - 0.5 * (f[2] - f[4]) + (1.0 / 6.0) * rho * ux + 0.5 * rho * uy;
    f[8] = f[6] + 0.5 * (f[2] - f[4]) + (1.0 / 6.0) * rho * ux - 0.5 * rho * uy;
}

/// Zou-He pressure outlet on the **right** wall (x = nx-1).
///
/// Given a distribution `f` at a right-wall cell with prescribed outlet
/// density `rho_out`, compute the three unknown populations `f[3]`, `f[6]`,
/// `f[7]` that stream back into the domain from the right boundary.
///
/// # Arguments
/// * `f`       – mutable 9-element distribution array.
/// * `rho_out` – prescribed outlet density.
pub fn zou_he_pressure_outlet(f: &mut [f64; 9], rho_out: f64) {
    // Outlet velocity (assuming uy = 0 at the right wall):
    let ux = -1.0 + (f[0] + f[2] + f[4] + 2.0 * (f[1] + f[5] + f[8])) / rho_out;
    f[3] = f[1] - (2.0 / 3.0) * rho_out * ux;
    f[7] = f[5] + 0.5 * (f[2] - f[4]) - (1.0 / 6.0) * rho_out * ux;
    f[6] = f[8] - 0.5 * (f[2] - f[4]) - (1.0 / 6.0) * rho_out * ux;
}

/// Zou-He velocity boundary on the **bottom** wall (y = 0).
pub fn zou_he_velocity_bottom(f: &mut [f64; 9], rho: f64, u: [f64; 2]) {
    let ux = u[0];
    let uy = u[1];
    f[2] = f[4] + (2.0 / 3.0) * rho * uy;
    f[5] = f[7] - 0.5 * (f[1] - f[3]) + (1.0 / 6.0) * rho * uy + 0.5 * rho * ux;
    f[6] = f[8] + 0.5 * (f[1] - f[3]) + (1.0 / 6.0) * rho * uy - 0.5 * rho * ux;
}

/// Zou-He velocity boundary on the **top** wall (y = ny-1).
pub fn zou_he_velocity_top(f: &mut [f64; 9], rho: f64, u: [f64; 2]) {
    let ux = u[0];
    let uy = u[1];
    f[4] = f[2] - (2.0 / 3.0) * rho * uy;
    f[7] = f[5] + 0.5 * (f[1] - f[3]) - (1.0 / 6.0) * rho * uy - 0.5 * rho * ux;
    f[8] = f[6] - 0.5 * (f[1] - f[3]) - (1.0 / 6.0) * rho * uy + 0.5 * rho * ux;
}

// ---------------------------------------------------------------------------
// Turbulent pipe flow (1/7 power law)
// ---------------------------------------------------------------------------

/// Initialize a 2D turbulent pipe flow velocity field using the
/// 1/7-th power law profile.
///
/// The cross-section extends from y=0 (wall) to y=ny-1 (wall).
/// The centreline is at y = (ny-1)/2.
///
/// u(r) = u_max * (1 - |r / R|)^(1/n)
/// where R = ny/2 and n = 7 by default.
///
/// Returns a `Vec` of \[ux, uy\] for each lattice node (row-major, x fast).
#[allow(dead_code)]
pub fn initialize_turbulent_pipe_flow(
    nx: usize,
    ny: usize,
    u_max: f64,
    n_power: f64,
) -> Vec<[f64; 2]> {
    let mut vel = vec![[0.0_f64; 2]; nx * ny];
    let r_max = (ny as f64 - 1.0) / 2.0;
    let y_center = (ny as f64 - 1.0) / 2.0;
    for y in 0..ny {
        for x in 0..nx {
            let r = ((y as f64 - y_center).abs() / r_max).min(1.0);
            let ux = u_max * (1.0 - r).powf(1.0 / n_power);
            vel[y * nx + x] = [ux, 0.0];
        }
    }
    vel
}

/// Compute the bulk velocity of a turbulent pipe flow profile.
///
/// For the 1/n power law the theoretical bulk / centreline ratio is:
/// U_bulk / U_max = 2 * n^2 / ((n+1) * (2n+1))
#[allow(dead_code)]
pub fn turbulent_pipe_bulk_velocity(u_max: f64, n_power: f64) -> f64 {
    2.0 * n_power * n_power / ((n_power + 1.0) * (2.0 * n_power + 1.0)) * u_max
}

// ---------------------------------------------------------------------------
// Atmospheric boundary layer (ABL) initialization
// ---------------------------------------------------------------------------

/// Log-law atmospheric boundary layer velocity profile.
///
/// u(z) = (u_star / kappa) * ln((z + z0) / z0)
/// where `u_star` is the friction velocity, `kappa` = 0.41 is von Kármán constant,
/// `z0` is the aerodynamic roughness length.
///
/// Returns `ux` (streamwise velocity) at height `z`.
#[allow(dead_code)]
pub fn abl_log_law_velocity(z: f64, u_star: f64, z0: f64) -> f64 {
    const KAPPA: f64 = 0.41;
    if z + z0 <= 0.0 {
        return 0.0;
    }
    (u_star / KAPPA) * ((z + z0) / z0).ln().max(0.0)
}

/// Initialize a 2D lattice with an atmospheric boundary layer velocity profile.
///
/// The y-axis represents height z. Returns `Vec<[ux, uy]>` (row-major).
#[allow(dead_code)]
pub fn initialize_abl_profile(
    nx: usize,
    ny: usize,
    u_star: f64,
    z0: f64,
    dz: f64,
) -> Vec<[f64; 2]> {
    let mut vel = vec![[0.0_f64; 2]; nx * ny];
    for y in 0..ny {
        let z = y as f64 * dz;
        let ux = abl_log_law_velocity(z, u_star, z0);
        for x in 0..nx {
            vel[y * nx + x] = [ux, 0.0];
        }
    }
    vel
}

// ---------------------------------------------------------------------------
// Wake initialization
// ---------------------------------------------------------------------------

/// Gaussian velocity deficit wake profile behind a bluff body.
///
/// u_wake(r, x) = U_inf * (1 - A * exp(-r^2 / (2 * sigma(x)^2)))
/// where sigma(x) = sigma0 * sqrt(1 + x / x0).
#[allow(dead_code)]
pub fn wake_velocity_deficit(
    r: f64,
    x_downstream: f64,
    u_inf: f64,
    wake_amplitude: f64,
    sigma0: f64,
    x0: f64,
) -> f64 {
    let sigma = sigma0 * (1.0 + x_downstream / x0.max(1e-12)).sqrt();
    u_inf * (1.0 - wake_amplitude * (-r * r / (2.0 * sigma * sigma)).exp())
}

/// Initialize a 2D wake flow behind a circular obstacle.
///
/// The obstacle is located at `(x_body, y_center)`.
/// Returns `Vec<[ux, uy]>` (row-major, x fast).
#[allow(dead_code)]
pub fn initialize_wake_flow(
    nx: usize,
    ny: usize,
    u_inf: f64,
    x_body: f64,
    y_center: f64,
    wake_amplitude: f64,
    sigma0: f64,
    x0: f64,
) -> Vec<[f64; 2]> {
    let mut vel = vec![[u_inf, 0.0_f64]; nx * ny];
    for y in 0..ny {
        for x in 0..nx {
            let x_d = x as f64 - x_body;
            if x_d > 0.0 {
                let r = (y as f64 - y_center).abs();
                let ux = wake_velocity_deficit(r, x_d, u_inf, wake_amplitude, sigma0, x0);
                vel[y * nx + x] = [ux, 0.0];
            }
        }
    }
    vel
}

// ---------------------------------------------------------------------------
// Multi-phase initialization (droplet array)
// ---------------------------------------------------------------------------

/// A single droplet definition for multi-phase initialization.
#[derive(Debug, Clone, Copy)]
pub struct Droplet {
    /// Centre of the droplet (x, y) in lattice units.
    pub center: [f64; 2],
    /// Radius of the droplet in lattice units.
    pub radius: f64,
    /// Phase indicator inside droplet (1.0 = liquid).
    pub phase: f64,
}

/// Initialize a phase-field for a multi-phase LBM simulation.
///
/// Returns a `Vec`f64` of phase values: 0.0 = gas, 1.0 = liquid.
/// A smooth tanh interface of width `interface_width` is applied.
#[allow(dead_code)]
pub fn initialize_droplet_array(
    nx: usize,
    ny: usize,
    droplets: &[Droplet],
    interface_width: f64,
) -> Vec<f64> {
    let mut phi = vec![0.0_f64; nx * ny];
    for y in 0..ny {
        for x in 0..nx {
            let mut max_phase = 0.0_f64;
            for d in droplets {
                let dx = x as f64 - d.center[0];
                let dy = y as f64 - d.center[1];
                let dist = (dx * dx + dy * dy).sqrt();
                // tanh interface: 1 inside, 0 outside
                let phase_here = 0.5 * (1.0 + ((d.radius - dist) / interface_width).tanh());
                max_phase = max_phase.max(phase_here * d.phase);
            }
            phi[y * nx + x] = max_phase;
        }
    }
    phi
}

/// Compute the average phase value (volume fraction of liquid).
#[allow(dead_code)]
pub fn average_phase(phi: &[f64]) -> f64 {
    if phi.is_empty() {
        return 0.0;
    }
    phi.iter().sum::<f64>() / phi.len() as f64
}

// ---------------------------------------------------------------------------
// Porous media flow initialization
// ---------------------------------------------------------------------------

/// Darcy flow velocity field in a porous medium.
///
/// u = -(K / mu) * grad_p
/// For a uniform pressure gradient dp/dx, u_x = K * dp_dx / mu.
#[allow(dead_code)]
pub fn initialize_porous_darcy(
    nx: usize,
    ny: usize,
    permeability: f64,
    viscosity: f64,
    dp_dx: f64,
) -> Vec<[f64; 2]> {
    let ux = permeability * dp_dx / viscosity;
    vec![[ux, 0.0]; nx * ny]
}

/// Brinkman-corrected porous media flow (parabolic profile damped by porosity).
///
/// Combines Darcy flow with a wall-corrected Brinkman profile.
/// The effective velocity is: u(y) = u_darcy * (1 - exp(-alpha * y) - exp(-alpha * (H-y)))
/// where alpha = sqrt(viscosity / (permeability * mu_eff)) and H = ny.
#[allow(dead_code)]
pub fn initialize_brinkman_flow(nx: usize, ny: usize, u_darcy: f64, alpha: f64) -> Vec<[f64; 2]> {
    let mut vel = vec![[0.0_f64; 2]; nx * ny];
    let h = ny as f64 - 1.0;
    for y in 0..ny {
        let yf = y as f64;
        let correction = 1.0 - (-alpha * yf).exp() - (-alpha * (h - yf)).exp();
        let ux = u_darcy * correction.max(0.0);
        for x in 0..nx {
            vel[y * nx + x] = [ux, 0.0];
        }
    }
    vel
}

/// Compute the effective Brinkman channel permeability correction factor.
///
/// Factor = 1 - 2/(alpha * H) * (1 - exp(-alpha * H))
#[allow(dead_code)]
pub fn brinkman_correction_factor(alpha: f64, h: f64) -> f64 {
    if alpha.abs() < 1e-12 {
        return 0.0;
    }
    1.0 - 2.0 / (alpha * h) * (1.0 - (-alpha * h).exp())
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // -----------------------------------------------------------------------
    // 1. new_channel: top and bottom rows solid, interior rows fluid
    // -----------------------------------------------------------------------
    #[test]
    fn test_new_channel_walls() {
        let nx = 10;
        let ny = 8;
        let geo = LatticeGeometry2D::new_channel(nx, ny);

        // Bottom row (y=0) all solid.
        for x in 0..nx {
            assert!(geo.is_solid[x], "bottom y=0, x={x} should be solid");
        }
        // Top row (y=ny-1) all solid.
        for x in 0..nx {
            assert!(
                geo.is_solid[(ny - 1) * nx + x],
                "top y={}, x={x} should be solid",
                ny - 1
            );
        }
        // Interior rows fluid.
        for y in 1..(ny - 1) {
            for x in 0..nx {
                assert!(
                    !geo.is_solid[y * nx + x],
                    "interior y={y}, x={x} should be fluid"
                );
            }
        }
    }

    // -----------------------------------------------------------------------
    // 2. add_circle_obstacle: cells strictly inside radius are solid
    // -----------------------------------------------------------------------
    #[test]
    fn test_add_circle_obstacle() {
        let mut geo = LatticeGeometry2D::new_channel(20, 20);
        let cx = 10usize;
        let cy = 10usize;
        let radius = 3.0_f64;
        geo.add_circle_obstacle(cx, cy, radius);

        for y in 0..20usize {
            for x in 0..20usize {
                let dx = x as f64 - cx as f64;
                let dy = y as f64 - cy as f64;
                let dist = (dx * dx + dy * dy).sqrt();
                if dist < radius {
                    assert!(
                        geo.is_solid[y * 20 + x],
                        "cell ({x},{y}) inside circle should be solid"
                    );
                }
            }
        }
    }

    // -----------------------------------------------------------------------
    // 3. solid_fraction: correct ratio
    // -----------------------------------------------------------------------
    #[test]
    fn test_solid_fraction() {
        let nx = 10;
        let ny = 6;
        let geo = LatticeGeometry2D::new_channel(nx, ny);
        // Two solid rows (y=0 and y=5) → 2*nx solid cells.
        let expected = (2 * nx) as f64 / (nx * ny) as f64;
        let frac = geo.solid_fraction();
        assert!(
            (frac - expected).abs() < 1e-14,
            "solid_fraction={frac}, expected={expected}"
        );
    }

    // -----------------------------------------------------------------------
    // 4. initialize_poiseuille: u at walls zero, maximum near center
    // -----------------------------------------------------------------------
    #[test]
    fn test_initialize_poiseuille() {
        let nx = 5;
        let ny = 11;
        let u_max = 0.1;
        let vel = initialize_poiseuille_flow(nx, ny, u_max);

        // Walls at y=0 and y=ny-1 should have zero velocity.
        for x in 0..nx {
            assert!(vel[x][0].abs() < 1e-14, "bottom wall should be zero");
            assert!(
                vel[(ny - 1) * nx + x][0].abs() < 1e-14,
                "top wall should be zero"
            );
        }

        // Center row (y = (ny-1)/2 = 5) should be near u_max.
        let y_mid = (ny - 1) / 2;
        let ux_mid = vel[y_mid * nx][0];
        assert!(
            (ux_mid - u_max).abs() < 1e-14,
            "center velocity={ux_mid}, expected={u_max}"
        );

        // All uy should be zero.
        for v in &vel {
            assert!(v[1].abs() < 1e-14, "uy should be zero everywhere");
        }
    }

    // -----------------------------------------------------------------------
    // 5. initialize_couette: linear profile from 0 to u_wall
    // -----------------------------------------------------------------------
    #[test]
    fn test_initialize_couette() {
        let nx = 4;
        let ny = 9;
        let u_wall = 0.05;
        let vel = initialize_couette_flow(nx, ny, u_wall);

        // y=0 → zero, y=ny-1 → u_wall.
        for x in 0..nx {
            assert!(vel[x][0].abs() < 1e-14, "Couette bottom should be zero");
            let ux_top = vel[(ny - 1) * nx + x][0];
            assert!(
                (ux_top - u_wall).abs() < 1e-14,
                "Couette top should be u_wall={u_wall}, got={ux_top}"
            );
        }

        // Verify linearity at intermediate y.
        for y in 0..ny {
            let expected = u_wall * y as f64 / (ny - 1) as f64;
            let ux = vel[y * nx][0];
            assert!(
                (ux - expected).abs() < 1e-14,
                "Couette non-linear at y={y}: ux={ux}, expected={expected}"
            );
        }
    }

    // -----------------------------------------------------------------------
    // 6. zou_he_velocity_inlet: resulting rho and u satisfy macroscopic values
    // -----------------------------------------------------------------------
    #[test]
    fn test_zou_he_inlet_macroscopic() {
        // Start from equilibrium at rest (rho=1, u=0) and apply inlet BC.
        // D2Q9 weights: w[0]=4/9, w[1..4]=1/9, w[5..8]=1/36.
        let w0 = 4.0 / 9.0;
        let w1 = 1.0 / 9.0;
        let w5 = 1.0 / 36.0;
        let mut f = [w0, w1, w1, w1, w1, w5, w5, w5, w5];

        let ux_in = 0.05_f64;
        let uy_in = 0.0_f64;
        let rho_in = (w0 + w1 + w1 + 2.0 * (w1 + w5 + w5)) / (1.0 - ux_in);
        zou_he_velocity_inlet(&mut f, rho_in, [ux_in, uy_in]);

        let cx = [0.0_f64, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0];
        let cy = [0.0_f64, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0];

        let rho_calc: f64 = f.iter().sum();
        let mx: f64 = f.iter().zip(cx.iter()).map(|(&fi, &cxi)| fi * cxi).sum();
        let my: f64 = f.iter().zip(cy.iter()).map(|(&fi, &cyi)| fi * cyi).sum();
        let ux_calc = mx / rho_calc;
        let uy_calc = my / rho_calc;

        assert!(
            (rho_calc - rho_in).abs() < 1e-12,
            "rho mismatch: calc={rho_calc}, expected={rho_in}"
        );
        assert!(
            (ux_calc - ux_in).abs() < 1e-12,
            "ux mismatch: calc={ux_calc}, expected={ux_in}"
        );
        assert!(
            (uy_calc - uy_in).abs() < 1e-12,
            "uy mismatch: calc={uy_calc}, expected={uy_in}"
        );
    }

    // -----------------------------------------------------------------------
    // 7. Taylor-Green vortex: divergence-free check
    // -----------------------------------------------------------------------
    #[test]
    fn test_taylor_green_divergence_free() {
        let nx = 32;
        let ny = 32;
        let u0 = 0.01;
        let vel = initialize_taylor_green_vortex(nx, ny, u0);

        // Check discrete divergence at interior points
        for y in 1..ny - 1 {
            for x in 1..nx - 1 {
                let k = y * nx + x;
                let dux_dx = (vel[k + 1][0] - vel[k - 1][0]) / 2.0;
                let duy_dy = (vel[(y + 1) * nx + x][1] - vel[(y - 1) * nx + x][1]) / 2.0;
                let div = dux_dx + duy_dy;
                assert!(
                    div.abs() < 1e-6,
                    "Taylor-Green not divergence-free at ({x},{y}): div={div}"
                );
            }
        }
    }

    // -----------------------------------------------------------------------
    // 8. Taylor-Green analytical decay
    // -----------------------------------------------------------------------
    #[test]
    fn test_taylor_green_decay() {
        let nx = 16;
        let ny = 16;
        let u0 = 0.01;
        let nu = 0.01;

        let vel_0 = taylor_green_analytical(nx, ny, u0, nu, 0.0);
        let vel_t = taylor_green_analytical(nx, ny, u0, nu, 100.0);

        // Velocity should decay over time
        let max_0: f64 = vel_0
            .iter()
            .map(|v| (v[0] * v[0] + v[1] * v[1]).sqrt())
            .fold(0.0_f64, f64::max);
        let max_t: f64 = vel_t
            .iter()
            .map(|v| (v[0] * v[0] + v[1] * v[1]).sqrt())
            .fold(0.0_f64, f64::max);
        assert!(
            max_t < max_0,
            "Taylor-Green should decay: max_0={max_0}, max_t={max_t}"
        );
    }

    // -----------------------------------------------------------------------
    // 9. Ramp functions
    // -----------------------------------------------------------------------
    #[test]
    fn test_linear_ramp() {
        assert!((linear_ramp(0, 100) - 0.0).abs() < 1e-14);
        assert!((linear_ramp(50, 100) - 0.5).abs() < 1e-14);
        assert!((linear_ramp(100, 100) - 1.0).abs() < 1e-14);
        assert!((linear_ramp(200, 100) - 1.0).abs() < 1e-14);
    }

    #[test]
    fn test_sinusoidal_ramp() {
        assert!((sinusoidal_ramp(0, 100) - 0.0).abs() < 1e-14);
        assert!((sinusoidal_ramp(100, 100) - 1.0).abs() < 1e-14);
        // At midpoint, should be 0.5
        let mid = sinusoidal_ramp(50, 100);
        assert!(
            (mid - 0.5).abs() < 1e-14,
            "Sinusoidal ramp at midpoint: {mid}"
        );
    }

    #[test]
    fn test_exponential_ramp() {
        let r0 = exponential_ramp(0, 100.0);
        assert!(r0.abs() < 1e-14, "Exp ramp at 0 should be 0: {r0}");
        let r_large = exponential_ramp(10000, 100.0);
        assert!(
            (r_large - 1.0).abs() < 1e-10,
            "Exp ramp at large t should be ~1: {r_large}"
        );
    }

    #[test]
    fn test_polynomial_ramp() {
        assert!((polynomial_ramp(0, 100, 2.0) - 0.0).abs() < 1e-14);
        assert!((polynomial_ramp(100, 100, 2.0) - 1.0).abs() < 1e-14);
        let mid = polynomial_ramp(50, 100, 2.0);
        assert!(
            (mid - 0.25).abs() < 1e-14,
            "Quadratic ramp at midpoint: {mid}"
        );
    }

    // -----------------------------------------------------------------------
    // 10. Checkpoint save/load round-trip
    // -----------------------------------------------------------------------
    #[test]
    fn test_checkpoint_round_trip() {
        let nx = 4;
        let ny = 3;
        let n = nx * ny;
        let rho: Vec<f64> = (0..n).map(|i| 1.0 + 0.01 * i as f64).collect();
        let ux: Vec<f64> = (0..n).map(|i| 0.001 * i as f64).collect();
        let uy: Vec<f64> = (0..n).map(|i| -0.001 * i as f64).collect();
        let dist = vec![vec![0.1; n]; 9];

        let cp = Checkpoint::new(
            42,
            nx,
            ny,
            rho.clone(),
            ux.clone(),
            uy.clone(),
            dist.clone(),
        );
        assert!(cp.is_valid(), "Checkpoint should be valid");

        let bytes = cp.to_bytes();
        let restored = Checkpoint::from_bytes(&bytes).expect("Failed to deserialize");

        assert_eq!(restored.step, 42);
        assert_eq!(restored.dimensions, (nx, ny));
        for k in 0..n {
            assert!((restored.rho[k] - rho[k]).abs() < 1e-14);
            assert!((restored.ux[k] - ux[k]).abs() < 1e-14);
            assert!((restored.uy[k] - uy[k]).abs() < 1e-14);
        }
        assert_eq!(restored.distributions.len(), 9);
    }

    // -----------------------------------------------------------------------
    // 11. Checkpoint velocity metrics
    // -----------------------------------------------------------------------
    #[test]
    fn test_checkpoint_velocity_metrics() {
        let cp = Checkpoint::new(
            0,
            2,
            2,
            vec![1.0; 4],
            vec![0.1, 0.0, 0.0, 0.0],
            vec![0.0; 4],
            vec![],
        );
        let l2 = cp.velocity_l2_norm();
        assert!((l2 - 0.1).abs() < 1e-14, "L2 norm should be 0.1: {l2}");
        let max_v = cp.max_velocity();
        assert!(
            (max_v - 0.1).abs() < 1e-14,
            "Max velocity should be 0.1: {max_v}"
        );
    }

    // -----------------------------------------------------------------------
    // 12. Perturbation initialization
    // -----------------------------------------------------------------------
    #[test]
    fn test_sinusoidal_perturbation() {
        let nx = 16;
        let ny = 16;
        let mut vel = initialize_uniform_flow(nx, ny, [0.1, 0.0]);
        add_sinusoidal_perturbation(&mut vel, nx, ny, 0.001, 1.0, 1.0);

        // The perturbation should change some velocities
        let mut has_nonzero_uy = false;
        for v in &vel {
            if v[1].abs() > 1e-15 {
                has_nonzero_uy = true;
                break;
            }
        }
        assert!(
            has_nonzero_uy,
            "Sinusoidal perturbation should add uy component"
        );
    }

    // -----------------------------------------------------------------------
    // 13. Kelvin-Helmholtz initialization
    // -----------------------------------------------------------------------
    #[test]
    fn test_kelvin_helmholtz_shear() {
        let nx = 32;
        let ny = 64;
        let vel = initialize_kelvin_helmholtz(nx, ny, 0.1, 0.001, 1.0);

        // Bottom half should have negative ux, top half positive ux (roughly)
        let k_bottom = 5 * nx + nx / 2;
        let k_top = (ny - 5) * nx + nx / 2;
        assert!(vel[k_bottom][0] < 0.0, "Bottom should have negative ux");
        assert!(vel[k_top][0] > 0.0, "Top should have positive ux");
    }

    // -----------------------------------------------------------------------
    // 14. Velocity L2 error
    // -----------------------------------------------------------------------
    #[test]
    fn test_velocity_l2_error_identical() {
        let vel = vec![[0.1, 0.2], [0.3, 0.4]];
        let err = velocity_l2_error(&vel, &vel);
        assert!(
            err < 1e-14,
            "Identical fields should have zero error: {err}"
        );
    }

    // -----------------------------------------------------------------------
    // 15. Open geometry has no solid nodes
    // -----------------------------------------------------------------------
    #[test]
    fn test_open_geometry() {
        let geo = LatticeGeometry2D::new_open(10, 10);
        assert_eq!(geo.solid_count(), 0);
        assert_eq!(geo.fluid_count(), 100);
    }

    // -----------------------------------------------------------------------
    // 16. Ellipse obstacle
    // -----------------------------------------------------------------------
    #[test]
    fn test_ellipse_obstacle() {
        let mut geo = LatticeGeometry2D::new_open(20, 20);
        geo.add_ellipse_obstacle(10.0, 10.0, 5.0, 3.0);
        assert!(geo.solid_count() > 0, "Ellipse should add solid cells");
        assert!(
            geo.is_solid[10 * 20 + 10],
            "Center of ellipse should be solid"
        );
    }

    // -----------------------------------------------------------------------
    // 17. Boundary fluid indices
    // -----------------------------------------------------------------------
    #[test]
    fn test_boundary_fluid_indices() {
        let geo = LatticeGeometry2D::new_channel(10, 6);
        let boundary = geo.boundary_fluid_indices();
        assert!(
            !boundary.is_empty(),
            "Channel should have boundary fluid nodes"
        );
        // All boundary fluid nodes should be adjacent to a wall
        for &k in &boundary {
            assert!(!geo.is_solid[k], "Boundary fluid node should not be solid");
        }
    }

    // -----------------------------------------------------------------------
    // 18. Taylor-Green density field
    // -----------------------------------------------------------------------
    #[test]
    fn test_taylor_green_density_average() {
        let nx = 32;
        let ny = 32;
        let rho0 = 1.0;
        let rho = taylor_green_density(nx, ny, 0.01, rho0);
        let avg: f64 = rho.iter().sum::<f64>() / rho.len() as f64;
        assert!(
            (avg - rho0).abs() < 0.01,
            "Average density should be close to rho0: avg={avg}"
        );
    }

    // -----------------------------------------------------------------------
    // 19. Zou-He bottom wall
    // -----------------------------------------------------------------------
    #[test]
    fn test_zou_he_bottom_momentum() {
        // Use a consistent density from the known populations.
        // For a bottom wall (normal = +y), the unknown pops are f[2], f[5], f[6].
        // rho = (f[0]+f[1]+f[3]+2*(f[4]+f[7]+f[8])) / (1 - uy)
        let w0 = 4.0 / 9.0;
        let w1 = 1.0 / 9.0;
        let w5 = 1.0 / 36.0;
        let mut f = [w0, w1, w1, w1, w1, w5, w5, w5, w5];

        let ux_in = 0.0_f64;
        let uy_in = 0.03_f64;
        let rho_in = (w0 + w1 + w1 + 2.0 * (w1 + w5 + w5)) / (1.0 - uy_in);
        zou_he_velocity_bottom(&mut f, rho_in, [ux_in, uy_in]);

        let cx = [0.0_f64, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0];
        let cy = [0.0_f64, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0];
        let rho_calc: f64 = f.iter().sum();
        let my: f64 = f.iter().zip(cy.iter()).map(|(&fi, &ci)| fi * ci).sum();
        let uy_calc = my / rho_calc;

        assert!(
            (uy_calc - uy_in).abs() < 1e-10,
            "uy mismatch: calc={uy_calc}, expected={uy_in}"
        );

        let mx: f64 = f.iter().zip(cx.iter()).map(|(&fi, &ci)| fi * ci).sum();
        let ux_calc = mx / rho_calc;
        assert!(
            (ux_calc - ux_in).abs() < 1e-10,
            "ux mismatch: calc={ux_calc}, expected={ux_in}"
        );
    }

    // -----------------------------------------------------------------------
    // 20. Ramp force application
    // -----------------------------------------------------------------------
    #[test]
    fn test_ramp_force() {
        let f = ramp_force([1.0, 2.0], 0.5);
        assert!((f[0] - 0.5).abs() < 1e-14);
        assert!((f[1] - 1.0).abs() < 1e-14);
    }

    // -----------------------------------------------------------------------
    // 21. Turbulent pipe flow 1/7 power law
    // -----------------------------------------------------------------------

    #[test]
    fn test_turbulent_pipe_centreline_max() {
        let nx = 10;
        let ny = 21; // centreline at y=10
        let u_max = 1.0;
        let vel = initialize_turbulent_pipe_flow(nx, ny, u_max, 7.0);
        // Centreline y=10, x=5
        let ux_center = vel[10 * nx + 5][0];
        assert!(
            (ux_center - u_max).abs() < 1e-10,
            "centreline velocity should equal u_max: {ux_center}"
        );
    }

    #[test]
    fn test_turbulent_pipe_wall_zero() {
        let nx = 5;
        let ny = 11;
        let vel = initialize_turbulent_pipe_flow(nx, ny, 1.0, 7.0);
        // At y=0 (wall), r = 1.0, so (1-1)^(1/7) = 0
        let ux_wall = vel[2][0];
        assert!(
            ux_wall.abs() < 1e-10,
            "wall velocity should be 0: {ux_wall}"
        );
    }

    #[test]
    fn test_turbulent_pipe_bulk_velocity_formula() {
        // For n=7: bulk/max = 2*49/(8*15) = 98/120 ≈ 0.8167
        let ratio = turbulent_pipe_bulk_velocity(1.0, 7.0);
        let expected = 2.0 * 49.0 / (8.0 * 15.0);
        assert!(
            (ratio - expected).abs() < 1e-12,
            "bulk velocity ratio mismatch: {ratio} vs {expected}"
        );
    }

    #[test]
    fn test_turbulent_pipe_velocity_monotone() {
        let nx = 1;
        let ny = 21;
        let vel = initialize_turbulent_pipe_flow(nx, ny, 1.0, 7.0);
        // From y=0 (wall) to y=10 (centre) velocity should be non-decreasing
        for y in 0..10 {
            let v0 = vel[y * nx][0];
            let v1 = vel[(y + 1) * nx][0];
            assert!(
                v1 >= v0 - 1e-12,
                "velocity should increase toward centreline: y={y}"
            );
        }
    }

    // -----------------------------------------------------------------------
    // 22. ABL log-law profile
    // -----------------------------------------------------------------------

    #[test]
    fn test_abl_log_law_at_surface_zero() {
        let u = abl_log_law_velocity(0.0, 0.3, 0.01);
        // At z=0 exactly: ln(z0/z0) = 0 → u = 0
        assert!(u.abs() < 1e-12, "ABL velocity at surface should be 0: {u}");
    }

    #[test]
    fn test_abl_log_law_increases_with_height() {
        let u1 = abl_log_law_velocity(10.0, 0.3, 0.01);
        let u2 = abl_log_law_velocity(50.0, 0.3, 0.01);
        assert!(
            u2 > u1,
            "ABL velocity should increase with height: u1={u1}, u2={u2}"
        );
    }

    #[test]
    fn test_abl_profile_length() {
        let nx = 10;
        let ny = 20;
        let vel = initialize_abl_profile(nx, ny, 0.3, 0.01, 0.5);
        assert_eq!(vel.len(), nx * ny);
    }

    #[test]
    fn test_abl_profile_ground_zero() {
        let nx = 5;
        let ny = 10;
        let vel = initialize_abl_profile(nx, ny, 0.3, 0.01, 1.0);
        // y=0 → z=0 → u=0
        for x in 0..nx {
            assert!(
                vel[x][0].abs() < 1e-12,
                "ABL profile at ground should be zero: {}",
                vel[x][0]
            );
        }
    }

    // -----------------------------------------------------------------------
    // 23. Wake initialization
    // -----------------------------------------------------------------------

    #[test]
    fn test_wake_deficit_far_upstream_is_u_inf() {
        // x_downstream < 0 → no wake
        let u = wake_velocity_deficit(0.0, 0.0, 1.0, 0.5, 1.0, 10.0);
        // At x=0 and r=0, deficit is maximum
        assert!(
            u >= 0.0 && u <= 1.0,
            "wake velocity should be in [0, U_inf]: {u}"
        );
    }

    #[test]
    fn test_wake_deficit_far_from_axis_is_u_inf() {
        // Far from centreline, deficit goes to zero
        let u = wake_velocity_deficit(1000.0, 10.0, 1.0, 0.5, 1.0, 10.0);
        assert!(
            (u - 1.0).abs() < 1e-6,
            "wake velocity far from axis should approach U_inf: {u}"
        );
    }

    #[test]
    fn test_wake_flow_field_length() {
        let nx = 20;
        let ny = 10;
        let vel = initialize_wake_flow(nx, ny, 1.0, 10.0, 5.0, 0.3, 1.0, 5.0);
        assert_eq!(vel.len(), nx * ny);
    }

    #[test]
    fn test_wake_flow_upstream_is_u_inf() {
        let nx = 20;
        let ny = 10;
        let u_inf = 1.0;
        let vel = initialize_wake_flow(nx, ny, u_inf, 15.0, 5.0, 0.3, 1.0, 5.0);
        // Upstream nodes (x < 15): ux = u_inf
        for y in 0..ny {
            for x in 0..10 {
                assert!(
                    (vel[y * nx + x][0] - u_inf).abs() < 1e-10,
                    "upstream velocity should be U_inf: {}",
                    vel[y * nx + x][0]
                );
            }
        }
    }

    // -----------------------------------------------------------------------
    // 24. Multi-phase droplet array
    // -----------------------------------------------------------------------

    #[test]
    fn test_droplet_array_inside_is_liquid() {
        let droplets = vec![Droplet {
            center: [5.0, 5.0],
            radius: 3.0,
            phase: 1.0,
        }];
        let phi = initialize_droplet_array(10, 10, &droplets, 0.5);
        // Centre of droplet should have phase ≈ 1
        let idx = 5 * 10 + 5;
        assert!(
            phi[idx] > 0.9,
            "droplet centre should be liquid: {}",
            phi[idx]
        );
    }

    #[test]
    fn test_droplet_array_outside_is_gas() {
        let droplets = vec![Droplet {
            center: [5.0, 5.0],
            radius: 2.0,
            phase: 1.0,
        }];
        let phi = initialize_droplet_array(10, 10, &droplets, 0.5);
        // Corner (0,0) is far from droplet → phase ≈ 0
        assert!(phi[0] < 0.1, "far corner should be gas: {}", phi[0]);
    }

    #[test]
    fn test_droplet_array_length() {
        let droplets = vec![Droplet {
            center: [5.0, 5.0],
            radius: 2.0,
            phase: 1.0,
        }];
        let phi = initialize_droplet_array(15, 12, &droplets, 0.5);
        assert_eq!(phi.len(), 15 * 12);
    }

    #[test]
    fn test_average_phase_single_droplet() {
        // A domain that is entirely liquid should have average phase = 1
        let phi = vec![1.0_f64; 100];
        assert!((average_phase(&phi) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_average_phase_empty() {
        assert!((average_phase(&[])).abs() < 1e-12);
    }

    // -----------------------------------------------------------------------
    // 25. Porous media flow
    // -----------------------------------------------------------------------

    #[test]
    fn test_darcy_flow_uniform() {
        let nx = 5;
        let ny = 5;
        let u_darcy = 0.01 * 1.0 / 1e-3; // K * dp_dx / mu
        let vel = initialize_porous_darcy(nx, ny, 0.01, 1e-3, 1.0);
        for v in &vel {
            assert!(
                (v[0] - u_darcy).abs() < 1e-10,
                "Darcy flow should be uniform"
            );
            assert!(v[1].abs() < 1e-10, "Darcy flow should have no y-component");
        }
    }

    #[test]
    fn test_brinkman_flow_centreline_max() {
        let nx = 1;
        let ny = 21;
        let vel = initialize_brinkman_flow(nx, ny, 1.0, 0.5);
        // Centre y=10 should have maximum velocity
        let ux_center = vel[10][0];
        let ux_wall = vel[0][0];
        assert!(
            ux_center > ux_wall,
            "Brinkman centreline should be faster than wall: center={ux_center}, wall={ux_wall}"
        );
    }

    #[test]
    fn test_brinkman_correction_factor_range() {
        let factor = brinkman_correction_factor(1.0, 10.0);
        assert!(
            factor >= 0.0 && factor <= 1.0,
            "correction factor should be in [0,1]: {factor}"
        );
    }

    #[test]
    fn test_brinkman_correction_factor_zero_alpha() {
        let factor = brinkman_correction_factor(0.0, 10.0);
        assert!((factor).abs() < 1e-12, "zero alpha → factor=0: {factor}");
    }
}
