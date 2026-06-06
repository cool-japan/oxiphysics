//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
use super::functions::{D2Q9_CX, D2Q9_CY, D2Q9_OPP, D2Q9_Q, D2Q9_W};
use std::f64::consts::PI;

/// A full FLIP/PIC fluid simulation (CPU mock of GPU dispatch).
#[derive(Debug, Clone)]
pub struct FluidSimulation {
    /// Simulation configuration.
    pub config: FluidSimConfig,
    /// MAC grid.
    pub grid: MacGrid,
    /// Previous grid (for FLIP delta).
    pub grid_prev: MacGrid,
    /// FLIP particles.
    pub particles: Vec<FlipParticle>,
    /// Level-set function (negative = fluid).
    pub phi: Vec<f64>,
    /// Simulation time.
    pub time: f64,
    /// Step count.
    pub step_count: u64,
}
impl FluidSimulation {
    /// Create a new fluid simulation.
    pub fn new(config: FluidSimConfig) -> Self {
        let [nx, ny, nz] = config.grid_size;
        let dx = config.dx;
        let grid = MacGrid::new(nx, ny, nz, dx);
        let grid_prev = MacGrid::new(nx, ny, nz, dx);
        let phi = vec![1.0f64; nx * ny * nz];
        Self {
            config,
            grid,
            grid_prev,
            particles: Vec::new(),
            phi,
            time: 0.0,
            step_count: 0,
        }
    }
    /// Add particles in a rectangular region.
    pub fn add_fluid_block(&mut self, min: [f64; 3], max: [f64; 3], count: usize) {
        let dx = (max[0] - min[0]) / (count as f64).cbrt().ceil();
        let dy = (max[1] - min[1]) / (count as f64).cbrt().ceil();
        let dz = (max[2] - min[2]) / (count as f64).cbrt().ceil();
        let steps = (count as f64).cbrt().ceil() as usize;
        let mut added = 0;
        'outer: for k in 0..steps {
            for j in 0..steps {
                for i in 0..steps {
                    if added >= count {
                        break 'outer;
                    }
                    let pos = [
                        min[0] + (i as f64 + 0.5) * dx,
                        min[1] + (j as f64 + 0.5) * dy,
                        min[2] + (k as f64 + 0.5) * dz,
                    ];
                    self.particles.push(FlipParticle::new(pos));
                    added += 1;
                }
            }
        }
    }
    /// Advance the simulation by one time step.
    pub fn step(&mut self) {
        let dt = self.config.dt;
        let [nx, ny, nz] = self.config.grid_size;
        p2g_transfer(&self.particles, &mut self.grid);
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..=nx {
                    let u = self.grid.get_u(i.min(nx), j, k);
                    self.grid
                        .set_u(i.min(nx), j, k, u + dt * self.config.gravity[0]);
                }
                for i in 0..nx {
                    let v = self.grid.get_v(i, j.min(ny), k);
                    self.grid
                        .set_v(i, j.min(ny), k, v + dt * self.config.gravity[1]);
                }
            }
        }
        self.grid
            .jacobi_pressure_solve(self.config.density, dt, self.config.pressure_iters);
        self.grid.pressure_project(self.config.density, dt);
        let vorticity = compute_vorticity(&self.grid);
        vorticity_confinement(&mut self.grid, &vorticity, self.config.vorticity_eps, dt);
        let grid_copy = self.grid.clone();
        g2p_transfer(
            &mut self.particles,
            &self.grid,
            &grid_copy,
            self.config.flip_ratio,
        );
        self.grid_prev = self.grid.clone();
        self.time += dt;
        self.step_count += 1;
    }
    /// Return total kinetic energy of the particle system.
    pub fn kinetic_energy(&self) -> f64 {
        let m = 1.0;
        self.particles
            .iter()
            .map(|p| {
                let v2 = p.velocity[0] * p.velocity[0]
                    + p.velocity[1] * p.velocity[1]
                    + p.velocity[2] * p.velocity[2];
                0.5 * m * v2
            })
            .sum()
    }
}
/// Multi-GPU domain decomposition along the X axis.
///
/// Splits the particle set evenly across `n_devices` by X coordinate.
/// In a production system each sub-domain is transferred to a separate device
/// memory and processed by its own command queue.
pub struct MultiGpuDomain;
impl MultiGpuDomain {
    /// Decompose particle positions across `n_devices` by X coordinate.
    pub fn decompose_x(positions: &[[f64; 3]], n_devices: usize) -> Vec<GpuSubDomain> {
        if n_devices == 0 || positions.is_empty() {
            return Vec::new();
        }
        let x_min = positions.iter().map(|p| p[0]).fold(f64::INFINITY, f64::min);
        let x_max = positions
            .iter()
            .map(|p| p[0])
            .fold(f64::NEG_INFINITY, f64::max);
        let range = (x_max - x_min).max(1e-30);
        let slice = range / n_devices as f64;
        let mut domains: Vec<GpuSubDomain> = (0..n_devices)
            .map(|dev| GpuSubDomain {
                device_id: dev,
                x_range: [x_min + dev as f64 * slice, x_min + (dev + 1) as f64 * slice],
                particle_indices: Vec::new(),
            })
            .collect();
        for (idx, pos) in positions.iter().enumerate() {
            let dev = (((pos[0] - x_min) / slice) as usize).min(n_devices - 1);
            domains[dev].particle_indices.push(idx);
        }
        domains
    }
}
/// An SPH fluid particle.
#[derive(Debug, Clone)]
pub struct SphParticle {
    /// Position in world space.
    pub position: [f64; 3],
    /// Velocity.
    pub velocity: [f64; 3],
    /// Mass.
    pub mass: f64,
    /// Computed density.
    pub density: f64,
    /// Computed pressure.
    pub pressure: f64,
    /// Accumulated force.
    pub force: [f64; 3],
}
impl SphParticle {
    /// Create a new SPH particle.
    pub fn new(position: [f64; 3], mass: f64) -> Self {
        Self {
            position,
            velocity: [0.0; 3],
            mass,
            density: 0.0,
            pressure: 0.0,
            force: [0.0; 3],
        }
    }
}
/// Configuration for the GPU-accelerated fluid simulation.
#[derive(Debug, Clone)]
pub struct FluidSimConfig {
    /// Grid cell count on each axis.
    pub grid_size: [usize; 3],
    /// Grid cell size (m).
    pub dx: f64,
    /// Simulation time step (s).
    pub dt: f64,
    /// Fluid density (kg/m³).
    pub density: f64,
    /// Gravity vector.
    pub gravity: [f64; 3],
    /// Number of pressure solver iterations.
    pub pressure_iters: usize,
    /// Vorticity confinement strength.
    pub vorticity_eps: f64,
    /// Surface tension coefficient.
    pub surface_tension: f64,
    /// FLIP/PIC blending ratio.
    pub flip_ratio: f64,
}
/// A spatial hash / uniform grid neighbour list built on a simulated GPU.
///
/// In a real GPU implementation the binning is done in a single prefix-sum
/// pass; here we use a `HashMap` to keep the mock dependency-free.
#[derive(Debug, Clone)]
pub struct GpuNeighborList {
    /// Smoothing radius used when building the list.
    pub h: f64,
    /// Domain size in each axis.
    pub domain: [f64; 3],
    /// Number of grid cells in each axis.
    pub grid_dims: [usize; 3],
    /// Cell edge length.
    pub cell_size: f64,
    /// Neighbour indices for every particle (indexed by particle id).
    pub neighbor_data: Vec<Vec<usize>>,
}
impl GpuNeighborList {
    /// Build the neighbour list from a slice of particle positions.
    pub fn build(positions: &[[f64; 3]], h: f64, domain: [f64; 3]) -> Self {
        let cell_size = h;
        let grid_dims = [
            ((domain[0] / cell_size).ceil() as usize).max(1),
            ((domain[1] / cell_size).ceil() as usize).max(1),
            ((domain[2] / cell_size).ceil() as usize).max(1),
        ];
        let n = positions.len();
        let mut cells: std::collections::HashMap<(usize, usize, usize), Vec<usize>> =
            std::collections::HashMap::new();
        for (idx, &pos) in positions.iter().enumerate() {
            let ci = ((pos[0] / cell_size) as usize).min(grid_dims[0].saturating_sub(1));
            let cj = ((pos[1] / cell_size) as usize).min(grid_dims[1].saturating_sub(1));
            let ck = ((pos[2] / cell_size) as usize).min(grid_dims[2].saturating_sub(1));
            cells.entry((ci, cj, ck)).or_default().push(idx);
        }
        let mut neighbor_data: Vec<Vec<usize>> = vec![Vec::new(); n];
        for (idx, &pos_i) in positions.iter().enumerate() {
            let ci = ((pos_i[0] / cell_size) as isize).max(0) as usize;
            let cj = ((pos_i[1] / cell_size) as isize).max(0) as usize;
            let ck = ((pos_i[2] / cell_size) as isize).max(0) as usize;
            for dz in 0usize..3 {
                for dy in 0usize..3 {
                    for dx in 0usize..3 {
                        let ni = (ci + dx).wrapping_sub(1);
                        let nj = (cj + dy).wrapping_sub(1);
                        let nk = (ck + dz).wrapping_sub(1);
                        if ni >= grid_dims[0] || nj >= grid_dims[1] || nk >= grid_dims[2] {
                            continue;
                        }
                        if let Some(bucket) = cells.get(&(ni, nj, nk)) {
                            for &jdx in bucket {
                                if jdx == idx {
                                    continue;
                                }
                                let r = length3(sub3(pos_i, positions[jdx]));
                                if r <= h {
                                    neighbor_data[idx].push(jdx);
                                }
                            }
                        }
                    }
                }
            }
        }
        Self {
            h,
            domain,
            grid_dims,
            cell_size,
            neighbor_data,
        }
    }
    /// Return the neighbour indices of particle `idx`.
    pub fn neighbors_of(&self, idx: usize) -> &[usize] {
        &self.neighbor_data[idx]
    }
}
/// LBM cell type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LbmCellType {
    /// Fluid cell.
    Fluid,
    /// Solid wall (no-slip bounce-back).
    Solid,
    /// Inlet velocity boundary.
    Inlet,
    /// Outlet pressure boundary.
    Outlet,
}
/// LBM D2Q9 simulation on a 2D grid.
#[derive(Debug, Clone)]
pub struct LbmD2Q9 {
    /// Grid width.
    pub nx: usize,
    /// Grid height.
    pub ny: usize,
    /// Distribution functions f\[y * nx * Q + x * Q + q\].
    pub f: Vec<f64>,
    /// Temporary buffer for streaming.
    pub f_tmp: Vec<f64>,
    /// Cell types.
    pub cell_type: Vec<LbmCellType>,
    /// Relaxation parameter τ (0 < τ < 2, τ=1 → ν=1/6).
    pub tau: f64,
    /// Inverse τ (precomputed).
    pub inv_tau: f64,
}
impl LbmD2Q9 {
    /// Create a new D2Q9 LBM simulation.
    pub fn new(nx: usize, ny: usize, tau: f64) -> Self {
        let n = nx * ny * D2Q9_Q;
        let mut f = vec![0.0f64; n];
        for y in 0..ny {
            for x in 0..nx {
                for q in 0..D2Q9_Q {
                    f[(y * nx + x) * D2Q9_Q + q] = D2Q9_W[q];
                }
            }
        }
        let f_tmp = f.clone();
        let cell_type = vec![LbmCellType::Fluid; nx * ny];
        Self {
            nx,
            ny,
            f,
            f_tmp,
            cell_type,
            tau,
            inv_tau: 1.0 / tau,
        }
    }
    /// Get distribution function f at cell (x, y), direction q.
    pub fn get_f(&self, x: usize, y: usize, q: usize) -> f64 {
        self.f[(y * self.nx + x) * D2Q9_Q + q]
    }
    /// Set distribution function.
    pub fn set_f(&mut self, x: usize, y: usize, q: usize, val: f64) {
        self.f[(y * self.nx + x) * D2Q9_Q + q] = val;
    }
    /// Compute macroscopic density ρ at cell (x, y).
    pub fn density(&self, x: usize, y: usize) -> f64 {
        (0..D2Q9_Q).map(|q| self.get_f(x, y, q)).sum()
    }
    /// Compute macroscopic velocity (ux, uy) at cell (x, y).
    pub fn velocity(&self, x: usize, y: usize) -> [f64; 2] {
        let rho = self.density(x, y);
        if rho < 1e-15 {
            return [0.0, 0.0];
        }
        let ux: f64 = (0..D2Q9_Q)
            .map(|q| self.get_f(x, y, q) * D2Q9_CX[q] as f64)
            .sum::<f64>()
            / rho;
        let uy: f64 = (0..D2Q9_Q)
            .map(|q| self.get_f(x, y, q) * D2Q9_CY[q] as f64)
            .sum::<f64>()
            / rho;
        [ux, uy]
    }
    /// Compute equilibrium distribution f_eq.
    pub fn f_equilibrium(rho: f64, ux: f64, uy: f64, q: usize) -> f64 {
        let cx = D2Q9_CX[q] as f64;
        let cy = D2Q9_CY[q] as f64;
        let cu = cx * ux + cy * uy;
        let u2 = ux * ux + uy * uy;
        D2Q9_W[q] * rho * (1.0 + 3.0 * cu + 4.5 * cu * cu - 1.5 * u2)
    }
    /// BGK collision step: f_i → f_i - (1/τ)(f_i - f_i^eq).
    pub fn collide(&mut self) {
        for y in 0..self.ny {
            for x in 0..self.nx {
                if self.cell_type[y * self.nx + x] == LbmCellType::Solid {
                    continue;
                }
                let rho = self.density(x, y);
                let [ux, uy] = self.velocity(x, y);
                for q in 0..D2Q9_Q {
                    let f_eq = Self::f_equilibrium(rho, ux, uy, q);
                    let f_old = self.get_f(x, y, q);
                    let f_new = f_old - self.inv_tau * (f_old - f_eq);
                    self.set_f(x, y, q, f_new);
                }
            }
        }
    }
    /// Streaming step: propagate f to neighboring cells.
    pub fn stream(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        self.f_tmp.copy_from_slice(&self.f);
        for y in 0..ny {
            for x in 0..nx {
                for q in 0..D2Q9_Q {
                    let nx_cell = (x as i32 + D2Q9_CX[q]).rem_euclid(nx as i32) as usize;
                    let ny_cell = (y as i32 + D2Q9_CY[q]).rem_euclid(ny as i32) as usize;
                    let src = self.f_tmp[(y * nx + x) * D2Q9_Q + q];
                    self.f[(ny_cell * nx + nx_cell) * D2Q9_Q + q] = src;
                }
            }
        }
    }
    /// Bounce-back boundary conditions for solid cells.
    pub fn apply_bounce_back(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let mut updates = Vec::new();
        for y in 0..ny {
            for x in 0..nx {
                if self.cell_type[y * nx + x] != LbmCellType::Solid {
                    continue;
                }
                for (q, &opp) in D2Q9_OPP.iter().enumerate().take(D2Q9_Q) {
                    let val = self.get_f(x, y, q);
                    updates.push((x, y, opp, val));
                }
            }
        }
        for (x, y, q, val) in updates {
            self.set_f(x, y, q, val);
        }
    }
    /// Full LBM step: collide → stream → bounce-back.
    pub fn step(&mut self) {
        self.collide();
        self.stream();
        self.apply_bounce_back();
    }
    /// Set inlet velocity boundary condition on left wall.
    pub fn set_inlet_velocity(&mut self, ux_in: f64, uy_in: f64) {
        for y in 0..self.ny {
            self.cell_type[y * self.nx] = LbmCellType::Inlet;
            let rho = 1.0;
            for q in 0..D2Q9_Q {
                let f_eq = Self::f_equilibrium(rho, ux_in, uy_in, q);
                self.set_f(0, y, q, f_eq);
            }
        }
    }
    /// Set outlet pressure boundary (right wall, ρ = ρ_out).
    pub fn set_outlet_pressure(&mut self, rho_out: f64) {
        let nx = self.nx;
        for y in 0..self.ny {
            self.cell_type[y * nx + nx - 1] = LbmCellType::Outlet;
            let [ux, uy] = self.velocity(nx - 1, y);
            for q in 0..D2Q9_Q {
                let f_eq = Self::f_equilibrium(rho_out, ux, uy, q);
                self.set_f(nx - 1, y, q, f_eq);
            }
        }
    }
}
/// SPH smoothing kernels (Müller et al. 2003).
pub struct SphKernels;
impl SphKernels {
    /// Poly6 kernel for density estimation.
    ///
    /// W_poly6(r, h) = (315 / (64πh^9)) * (h² - r²)³  for r ≤ h
    pub fn poly6(r: f64, h: f64) -> f64 {
        if r > h || r < 0.0 {
            return 0.0;
        }
        let h2 = h * h;
        let r2 = r * r;
        let diff = h2 - r2;
        (315.0 / (64.0 * PI * h.powi(9))) * diff * diff * diff
    }
    /// Gradient of Poly6 kernel.
    pub fn poly6_grad(r_vec: [f64; 3], r: f64, h: f64) -> [f64; 3] {
        if r > h || r < 1e-15 {
            return [0.0; 3];
        }
        let h2 = h * h;
        let r2 = r * r;
        let coeff = -6.0 * (315.0 / (64.0 * PI * h.powi(9))) * (h2 - r2) * (h2 - r2);
        scale3(r_vec, coeff)
    }
    /// Spiky kernel gradient for pressure forces.
    ///
    /// ∇W_spiky(r, h) = -(45 / (πh^6)) * (h - r)² * r̂  for r ≤ h
    pub fn spiky_grad(r_vec: [f64; 3], r: f64, h: f64) -> [f64; 3] {
        if r > h || r < 1e-15 {
            return [0.0; 3];
        }
        let diff = h - r;
        let coeff = -(45.0 / (PI * h.powi(6))) * diff * diff / r;
        scale3(r_vec, coeff)
    }
    /// Viscosity kernel Laplacian.
    ///
    /// ∇²W_visc(r, h) = (45 / (πh^6)) * (h - r)  for r ≤ h
    pub fn viscosity_laplacian(r: f64, h: f64) -> f64 {
        if r > h {
            return 0.0;
        }
        (45.0 / (PI * h.powi(6))) * (h - r)
    }
}
/// A FLIP/PIC particle.
#[derive(Debug, Clone)]
pub struct FlipParticle {
    /// Particle position.
    pub position: [f64; 3],
    /// Particle velocity.
    pub velocity: [f64; 3],
    /// Accumulated grid velocity (for FLIP update).
    pub velocity_grid: [f64; 3],
}
impl FlipParticle {
    /// Create a new FLIP particle at rest.
    pub fn new(position: [f64; 3]) -> Self {
        Self {
            position,
            velocity: [0.0; 3],
            velocity_grid: [0.0; 3],
        }
    }
}
/// Axis-aligned bounding box for GPU boundary enforcement.
#[derive(Debug, Clone)]
pub struct GpuBoundaryBox {
    /// Minimum corner of the box.
    pub min: [f64; 3],
    /// Maximum corner of the box.
    pub max: [f64; 3],
    /// Coefficient of restitution (0 = inelastic, 1 = elastic).
    pub restitution: f64,
}
/// A single GPU sub-domain in a multi-GPU decomposition.
#[derive(Debug, Clone)]
pub struct GpuSubDomain {
    /// GPU device index (logical).
    pub device_id: usize,
    /// X-axis boundary \[x_min, x_max\].
    pub x_range: [f64; 2],
    /// Particle indices assigned to this device.
    pub particle_indices: Vec<usize>,
}
/// Marker-and-Cell (MAC) staggered grid for incompressible fluid.
///
/// u-velocity is stored at cell face centers on the X axis (i+½, j, k).
/// v-velocity is at (i, j+½, k).
/// w-velocity is at (i, j, k+½).
/// Pressure is at cell centers (i, j, k).
#[derive(Debug, Clone)]
pub struct MacGrid {
    /// Grid dimensions (cells).
    pub nx: usize,
    /// Grid dimensions (cells).
    pub ny: usize,
    /// Grid dimensions (cells).
    pub nz: usize,
    /// Cell size.
    pub dx: f64,
    /// u-velocity: (nx+1) × ny × nz.
    pub u: Vec<f64>,
    /// v-velocity: nx × (ny+1) × nz.
    pub v: Vec<f64>,
    /// w-velocity: nx × ny × (nz+1).
    pub w: Vec<f64>,
    /// Pressure: nx × ny × nz.
    pub p: Vec<f64>,
    /// Divergence of velocity field.
    pub div: Vec<f64>,
    /// Cell flags (0=air, 1=fluid, 2=solid).
    pub flags: Vec<u8>,
}
impl MacGrid {
    /// Create a new MAC grid.
    pub fn new(nx: usize, ny: usize, nz: usize, dx: f64) -> Self {
        Self {
            nx,
            ny,
            nz,
            dx,
            u: vec![0.0; (nx + 1) * ny * nz],
            v: vec![0.0; nx * (ny + 1) * nz],
            w: vec![0.0; nx * ny * (nz + 1)],
            p: vec![0.0; nx * ny * nz],
            div: vec![0.0; nx * ny * nz],
            flags: vec![0u8; nx * ny * nz],
        }
    }
    fn cell_idx(&self, i: usize, j: usize, k: usize) -> usize {
        k * self.nx * self.ny + j * self.nx + i
    }
    /// Get u-velocity at face (i+½, j, k).
    pub fn get_u(&self, i: usize, j: usize, k: usize) -> f64 {
        self.u[k * (self.nx + 1) * self.ny + j * (self.nx + 1) + i]
    }
    /// Set u-velocity.
    pub fn set_u(&mut self, i: usize, j: usize, k: usize, val: f64) {
        self.u[k * (self.nx + 1) * self.ny + j * (self.nx + 1) + i] = val;
    }
    /// Get v-velocity at face (i, j+½, k).
    pub fn get_v(&self, i: usize, j: usize, k: usize) -> f64 {
        self.v[k * self.nx * (self.ny + 1) + j * self.nx + i]
    }
    /// Set v-velocity.
    pub fn set_v(&mut self, i: usize, j: usize, k: usize, val: f64) {
        self.v[k * self.nx * (self.ny + 1) + j * self.nx + i] = val;
    }
    /// Get w-velocity at face (i, j, k+½).
    pub fn get_w(&self, i: usize, j: usize, k: usize) -> f64 {
        self.w[(k) * self.nx * self.ny + j * self.nx + i]
    }
    /// Set w-velocity.
    pub fn set_w(&mut self, i: usize, j: usize, k: usize, val: f64) {
        self.w[(k) * self.nx * self.ny + j * self.nx + i] = val;
    }
    /// Compute velocity divergence at each cell.
    pub fn compute_divergence(&mut self) {
        let inv_dx = 1.0 / self.dx;
        for k in 0..self.nz {
            for j in 0..self.ny {
                for i in 0..self.nx {
                    let du = self.get_u(i + 1, j, k) - self.get_u(i, j, k);
                    let dv = self.get_v(i, j + 1, k) - self.get_v(i, j, k);
                    let dw = self.get_w(i, j, k + 1) - self.get_w(i, j, k);
                    let idx = self.cell_idx(i, j, k);
                    self.div[idx] = (du + dv + dw) * inv_dx;
                }
            }
        }
    }
    /// Jacobi pressure solve: solve ∇²p = ρ/dt * ∇·u.
    ///
    /// Runs `iterations` Jacobi sweeps.
    pub fn jacobi_pressure_solve(&mut self, rho: f64, dt: f64, iterations: usize) {
        self.compute_divergence();
        let scale = rho * self.dx * self.dx / dt;
        let nx = self.nx;
        let ny = self.ny;
        let nz = self.nz;
        let mut p_new = self.p.clone();
        for _ in 0..iterations {
            for k in 0..nz {
                for j in 0..ny {
                    for i in 0..nx {
                        let idx = self.cell_idx(i, j, k);
                        if self.flags[idx] != 1 {
                            continue;
                        }
                        let mut neighbor_sum = 0.0;
                        let mut n_count = 0u32;
                        if i + 1 < nx {
                            neighbor_sum += self.p[self.cell_idx(i + 1, j, k)];
                            n_count += 1;
                        }
                        if i > 0 {
                            neighbor_sum += self.p[self.cell_idx(i - 1, j, k)];
                            n_count += 1;
                        }
                        if j + 1 < ny {
                            neighbor_sum += self.p[self.cell_idx(i, j + 1, k)];
                            n_count += 1;
                        }
                        if j > 0 {
                            neighbor_sum += self.p[self.cell_idx(i, j - 1, k)];
                            n_count += 1;
                        }
                        if k + 1 < nz {
                            neighbor_sum += self.p[self.cell_idx(i, j, k + 1)];
                            n_count += 1;
                        }
                        if k > 0 {
                            neighbor_sum += self.p[self.cell_idx(i, j, k - 1)];
                            n_count += 1;
                        }
                        if n_count > 0 {
                            p_new[idx] = (neighbor_sum - scale * self.div[idx]) / n_count as f64;
                        }
                    }
                }
            }
            self.p.copy_from_slice(&p_new);
        }
    }
    /// Apply pressure gradient to velocities (project step).
    pub fn pressure_project(&mut self, rho: f64, dt: f64) {
        let scale = dt / (rho * self.dx);
        let nx = self.nx;
        let ny = self.ny;
        let nz = self.nz;
        for k in 0..nz {
            for j in 0..ny {
                for i in 1..nx {
                    let p_r = self.p[self.cell_idx(i, j, k)];
                    let p_l = self.p[self.cell_idx(i - 1, j, k)];
                    let u_old = self.get_u(i, j, k);
                    self.set_u(i, j, k, u_old - scale * (p_r - p_l));
                }
            }
        }
        for k in 0..nz {
            for j in 1..ny {
                for i in 0..nx {
                    let p_t = self.p[self.cell_idx(i, j, k)];
                    let p_b = self.p[self.cell_idx(i, j - 1, k)];
                    let v_old = self.get_v(i, j, k);
                    self.set_v(i, j, k, v_old - scale * (p_t - p_b));
                }
            }
        }
        for k in 1..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let p_f = self.p[self.cell_idx(i, j, k)];
                    let p_bk = self.p[self.cell_idx(i, j, k - 1)];
                    let w_old = self.get_w(i, j, k);
                    self.set_w(i, j, k, w_old - scale * (p_f - p_bk));
                }
            }
        }
    }
    /// Trilinear interpolation of u-velocity at world position (x, y, z).
    pub fn interp_u(&self, x: f64, y: f64, z: f64) -> f64 {
        let inv_dx = 1.0 / self.dx;
        let ix = (x * inv_dx).floor() as isize;
        let jy = (y * inv_dx - 0.5).floor() as isize;
        let kz = (z * inv_dx - 0.5).floor() as isize;
        let fx = x * inv_dx - ix as f64;
        let fy = y * inv_dx - 0.5 - jy as f64;
        let fz = z * inv_dx - 0.5 - kz as f64;
        let nx = self.nx as isize;
        let ny = self.ny as isize;
        let nz = self.nz as isize;
        let clamp_i = |v: isize| v.clamp(0, nx) as usize;
        let clamp_j = |v: isize| v.clamp(0, ny - 1) as usize;
        let clamp_k = |v: isize| v.clamp(0, nz - 1) as usize;
        let u000 = self.get_u(clamp_i(ix), clamp_j(jy), clamp_k(kz));
        let u100 = self.get_u(clamp_i(ix + 1), clamp_j(jy), clamp_k(kz));
        let u010 = self.get_u(clamp_i(ix), clamp_j(jy + 1), clamp_k(kz));
        let u110 = self.get_u(clamp_i(ix + 1), clamp_j(jy + 1), clamp_k(kz));
        let u001 = self.get_u(clamp_i(ix), clamp_j(jy), clamp_k(kz + 1));
        let u101 = self.get_u(clamp_i(ix + 1), clamp_j(jy), clamp_k(kz + 1));
        let u011 = self.get_u(clamp_i(ix), clamp_j(jy + 1), clamp_k(kz + 1));
        let u111 = self.get_u(clamp_i(ix + 1), clamp_j(jy + 1), clamp_k(kz + 1));
        let lerp = |a: f64, b: f64, t: f64| a + t * (b - a);
        let u_x0 = lerp(u000, u100, fx);
        let u_x1 = lerp(u010, u110, fx);
        let u_x2 = lerp(u001, u101, fx);
        let u_x3 = lerp(u011, u111, fx);
        let u_y0 = lerp(u_x0, u_x1, fy);
        let u_y1 = lerp(u_x2, u_x3, fy);
        lerp(u_y0, u_y1, fz)
    }
}
/// SPH simulation configuration.
#[derive(Debug, Clone)]
pub struct SphConfig {
    /// Smoothing radius.
    pub h: f64,
    /// Rest density (kg/m³).
    pub rest_density: f64,
    /// Pressure stiffness.
    pub pressure_k: f64,
    /// Viscosity coefficient.
    pub viscosity: f64,
    /// Surface tension coefficient.
    pub surface_tension: f64,
    /// Gravity.
    pub gravity: [f64; 3],
    /// Time step.
    pub dt: f64,
}
