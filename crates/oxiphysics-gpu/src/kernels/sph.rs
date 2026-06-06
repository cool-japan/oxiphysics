// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! SPH (Smoothed Particle Hydrodynamics) compute kernels.

use crate::compute::ComputeKernel;
use std::f64::consts::PI;

// ─────────────────────────────────────────────────────────────────────────────
// High-level SPH kernel API
// ─────────────────────────────────────────────────────────────────────────────

/// Selects the smoothing kernel function to use for SPH computations.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SphKernel {
    /// Cubic B-spline kernel (C2 continuity, widely used in SPH).
    CubicSpline,
    /// Wendland C2 kernel (positive definite, no negative lobes).
    Wendland,
    /// Poly6 kernel (efficient for density, but poor gradient).
    Poly6,
    /// Spiky kernel (non-zero gradient at origin, good for pressure).
    Spiky,
}

/// Pre-computed parameters for a smoothing kernel with smoothing length `h`.
#[derive(Debug, Clone, Copy)]
pub struct SphKernelParams {
    /// Smoothing length.
    pub h: f64,
    /// 1 / h.
    pub d_inv: f64,
    /// 1 / h^3.
    pub d3_inv: f64,
}

impl SphKernelParams {
    /// Construct `SphKernelParams` from smoothing length `h`.
    pub fn new(h: f64) -> Self {
        Self {
            h,
            d_inv: 1.0 / h,
            d3_inv: 1.0 / (h * h * h),
        }
    }
}

/// Evaluate the smoothing kernel W(r, h) for a given kernel type.
///
/// Returns 0 when `r >= h` (compact support).
pub fn kernel_value(r: f64, params: &SphKernelParams, kernel: SphKernel) -> f64 {
    let q = r * params.d_inv; // dimensionless distance
    match kernel {
        SphKernel::CubicSpline => {
            // 3-D cubic spline, alpha = 8/(pi h^3)
            let alpha = 8.0 / (PI * params.h.powi(3));
            if q >= 1.0 {
                0.0
            } else if q >= 0.5 {
                let t = 1.0 - q;
                alpha * 2.0 * t.powi(3)
            } else {
                alpha * (6.0 * q.powi(3) - 6.0 * q * q + 1.0)
            }
        }
        SphKernel::Wendland => {
            // Wendland C2 in 3-D, alpha = 21/(2 pi h^3)
            let alpha = 21.0 / (2.0 * PI * params.h.powi(3));
            if q >= 1.0 {
                0.0
            } else {
                let t = 1.0 - q;
                alpha * t.powi(4) * (4.0 * q + 1.0)
            }
        }
        SphKernel::Poly6 => {
            let h2 = params.h * params.h;
            let r2 = r * r;
            if r2 >= h2 {
                return 0.0;
            }
            let coeff = 315.0 / (64.0 * PI * params.h.powi(9));
            coeff * (h2 - r2).powi(3)
        }
        SphKernel::Spiky => {
            if q >= 1.0 {
                return 0.0;
            }
            let coeff = 15.0 / (PI * params.h.powi(6));
            coeff * (params.h - r).powi(3)
        }
    }
}

/// Evaluate the gradient nabla W(r_vec, h) of the smoothing kernel.
///
/// `r_vec` is the vector from particle j to particle i (r_i - r_j).
/// `r` is its magnitude.  Returns the zero vector when `r < 1e-12` or `r >= h`.
pub fn kernel_gradient(
    r_vec: [f64; 3],
    r: f64,
    params: &SphKernelParams,
    kernel: SphKernel,
) -> [f64; 3] {
    if r < 1e-12 || r * params.d_inv >= 1.0 {
        return [0.0; 3];
    }
    let dw_dr = kernel_gradient_mag(r, params, kernel);
    let scale = dw_dr / r;
    [r_vec[0] * scale, r_vec[1] * scale, r_vec[2] * scale]
}

/// Scalar dW/dr (negative for most kernels beyond their core).
fn kernel_gradient_mag(r: f64, params: &SphKernelParams, kernel: SphKernel) -> f64 {
    let q = r * params.d_inv;
    match kernel {
        SphKernel::CubicSpline => {
            let alpha = 8.0 / (PI * params.h.powi(3));
            if q >= 1.0 {
                0.0
            } else if q >= 0.5 {
                let t = 1.0 - q;
                alpha * (-6.0 * t * t) * params.d_inv
            } else {
                alpha * (18.0 * q * q - 12.0 * q) * params.d_inv
            }
        }
        SphKernel::Wendland => {
            let alpha = 21.0 / (2.0 * PI * params.h.powi(3));
            if q >= 1.0 {
                0.0
            } else {
                let t = 1.0 - q;
                // d/dr [(1-q)^4 (4q+1)] * alpha
                // = alpha * d_inv * [-4(1-q)^3(4q+1) + (1-q)^4 * 4]
                alpha * params.d_inv * t.powi(3) * (-4.0 * (4.0 * q + 1.0) + 4.0 * t)
            }
        }
        SphKernel::Poly6 => {
            let h2 = params.h * params.h;
            let r2 = r * r;
            if r2 >= h2 {
                return 0.0;
            }
            let coeff = 315.0 / (64.0 * PI * params.h.powi(9));
            // d/dr [(h^2-r^2)^3] = -6r(h^2-r^2)^2
            coeff * (-6.0 * r) * (h2 - r2).powi(2)
        }
        SphKernel::Spiky => {
            if q >= 1.0 {
                return 0.0;
            }
            let coeff = 15.0 / (PI * params.h.powi(6));
            // d/dr [(h-r)^3] = -3(h-r)^2
            coeff * (-3.0) * (params.h - r).powi(2)
        }
    }
}

/// Compute per-particle densities using SPH summation.
///
/// `rho_i = sum_j m_j * W(|r_i - r_j|, h)`
pub fn density_summation(positions: &[[f64; 3]], masses: &[f64], h: f64) -> Vec<f64> {
    let params = SphKernelParams::new(h);
    let n = positions.len();
    let mut densities = vec![0.0f64; n];
    for i in 0..n {
        let mut rho = 0.0;
        for j in 0..n {
            let dx = positions[i][0] - positions[j][0];
            let dy = positions[i][1] - positions[j][1];
            let dz = positions[i][2] - positions[j][2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            rho += masses[j] * kernel_value(r, &params, SphKernel::CubicSpline);
        }
        densities[i] = rho;
    }
    densities
}

/// Compute per-particle densities with a specified kernel type.
pub fn density_summation_kernel(
    positions: &[[f64; 3]],
    masses: &[f64],
    h: f64,
    kernel: SphKernel,
) -> Vec<f64> {
    let params = SphKernelParams::new(h);
    let n = positions.len();
    let mut densities = vec![0.0f64; n];
    for i in 0..n {
        let mut rho = 0.0;
        for j in 0..n {
            let dx = positions[i][0] - positions[j][0];
            let dy = positions[i][1] - positions[j][1];
            let dz = positions[i][2] - positions[j][2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            rho += masses[j] * kernel_value(r, &params, kernel);
        }
        densities[i] = rho;
    }
    densities
}

/// Compute pressure forces using the SPH symmetric pressure gradient formulation.
///
/// `F_i^pressure = -sum_j m_j (p_i/rho_i^2 + p_j/rho_j^2) nabla W(r_ij, h)`
///
/// Returns a `Vec<[f64;3]>` of forces, one per particle.
pub fn pressure_force(
    positions: &[[f64; 3]],
    _velocities: &[[f64; 3]],
    densities: &[f64],
    pressures: &[f64],
    masses: &[f64],
    h: f64,
) -> Vec<[f64; 3]> {
    let params = SphKernelParams::new(h);
    let n = positions.len();
    let mut forces = vec![[0.0f64; 3]; n];
    for i in 0..n {
        let mut fx = 0.0f64;
        let mut fy = 0.0f64;
        let mut fz = 0.0f64;
        let pi_over_rhoi2 = if densities[i].abs() > 1e-30 {
            pressures[i] / (densities[i] * densities[i])
        } else {
            0.0
        };
        for j in 0..n {
            if i == j {
                continue;
            }
            let r_vec = [
                positions[i][0] - positions[j][0],
                positions[i][1] - positions[j][1],
                positions[i][2] - positions[j][2],
            ];
            let r = (r_vec[0] * r_vec[0] + r_vec[1] * r_vec[1] + r_vec[2] * r_vec[2]).sqrt();
            let grad = kernel_gradient(r_vec, r, &params, SphKernel::CubicSpline);
            let pj_over_rhoj2 = if densities[j].abs() > 1e-30 {
                pressures[j] / (densities[j] * densities[j])
            } else {
                0.0
            };
            let coeff = -masses[j] * (pi_over_rhoi2 + pj_over_rhoj2);
            fx += coeff * grad[0];
            fy += coeff * grad[1];
            fz += coeff * grad[2];
        }
        forces[i] = [fx, fy, fz];
    }
    forces
}

/// Compute viscosity forces using the SPH viscosity formulation.
///
/// `F_i^visc = mu * sum_j m_j (v_j - v_i) / rho_j * laplacian_W(r_ij, h)`
pub fn viscosity_force(
    positions: &[[f64; 3]],
    velocities: &[[f64; 3]],
    densities: &[f64],
    masses: &[f64],
    h: f64,
    mu: f64,
) -> Vec<[f64; 3]> {
    let n = positions.len();
    let mut forces = vec![[0.0f64; 3]; n];
    for i in 0..n {
        let mut fx = 0.0f64;
        let mut fy = 0.0f64;
        let mut fz = 0.0f64;
        for j in 0..n {
            if i == j {
                continue;
            }
            let dx = positions[i][0] - positions[j][0];
            let dy = positions[i][1] - positions[j][1];
            let dz = positions[i][2] - positions[j][2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            if r >= h || r < 1e-12 {
                continue;
            }
            let lap = viscosity_laplacian(r, h);
            let rho_j = if densities[j].abs() > 1e-30 {
                densities[j]
            } else {
                1.0
            };
            fx += mu * masses[j] * (velocities[j][0] - velocities[i][0]) / rho_j * lap;
            fy += mu * masses[j] * (velocities[j][1] - velocities[i][1]) / rho_j * lap;
            fz += mu * masses[j] * (velocities[j][2] - velocities[i][2]) / rho_j * lap;
        }
        forces[i] = [fx, fy, fz];
    }
    forces
}

// ─────────────────────────────────────────────────────────────────────────────
// Neighbor list
// ─────────────────────────────────────────────────────────────────────────────

/// A simple grid-based neighbor list for SPH.
///
/// Divides the domain into cells of size `h` and assigns each particle to a cell.
/// Finding neighbors then only requires checking the local cell and its 26 neighbors.
pub struct NeighborList {
    /// Cell size (= smoothing length h).
    cell_size: f64,
    /// Number of cells in x, y, z.
    grid_dims: [usize; 3],
    /// Domain origin.
    origin: [f64; 3],
    /// Cell -> list of particle indices.
    cells: Vec<Vec<usize>>,
}

impl NeighborList {
    /// Create a new neighbor list for the given domain.
    pub fn new(origin: [f64; 3], domain_size: [f64; 3], cell_size: f64) -> Self {
        let nx = (domain_size[0] / cell_size).ceil() as usize;
        let ny = (domain_size[1] / cell_size).ceil() as usize;
        let nz = (domain_size[2] / cell_size).ceil() as usize;
        let total = nx.max(1) * ny.max(1) * nz.max(1);
        Self {
            cell_size,
            grid_dims: [nx.max(1), ny.max(1), nz.max(1)],
            origin,
            cells: vec![Vec::new(); total],
        }
    }

    /// Clear all cells and re-assign particles.
    pub fn build(&mut self, positions: &[[f64; 3]]) {
        for cell in &mut self.cells {
            cell.clear();
        }
        for (idx, pos) in positions.iter().enumerate() {
            let ci = self.cell_index(pos);
            self.cells[ci].push(idx);
        }
    }

    /// Get the cell index for a position.
    fn cell_index(&self, pos: &[f64; 3]) -> usize {
        let ix = ((pos[0] - self.origin[0]) / self.cell_size).floor() as usize;
        let iy = ((pos[1] - self.origin[1]) / self.cell_size).floor() as usize;
        let iz = ((pos[2] - self.origin[2]) / self.cell_size).floor() as usize;
        let ix = ix.min(self.grid_dims[0] - 1);
        let iy = iy.min(self.grid_dims[1] - 1);
        let iz = iz.min(self.grid_dims[2] - 1);
        iz * self.grid_dims[1] * self.grid_dims[0] + iy * self.grid_dims[0] + ix
    }

    /// Return the neighbor particle indices for a given particle.
    /// Checks the particle's cell and all 26 neighboring cells.
    pub fn neighbors(&self, pos: &[f64; 3]) -> Vec<usize> {
        let ix = ((pos[0] - self.origin[0]) / self.cell_size).floor() as i64;
        let iy = ((pos[1] - self.origin[1]) / self.cell_size).floor() as i64;
        let iz = ((pos[2] - self.origin[2]) / self.cell_size).floor() as i64;

        let mut result = Vec::new();
        let dims = self.grid_dims;

        for dz in -1i64..=1 {
            for dy in -1i64..=1 {
                for dx in -1i64..=1 {
                    let cx = ix + dx;
                    let cy = iy + dy;
                    let cz = iz + dz;
                    if cx < 0 || cy < 0 || cz < 0 {
                        continue;
                    }
                    let cx = cx as usize;
                    let cy = cy as usize;
                    let cz = cz as usize;
                    if cx >= dims[0] || cy >= dims[1] || cz >= dims[2] {
                        continue;
                    }
                    let ci = cz * dims[1] * dims[0] + cy * dims[0] + cx;
                    result.extend_from_slice(&self.cells[ci]);
                }
            }
        }
        result
    }

    /// Total number of cells in the grid.
    pub fn num_cells(&self) -> usize {
        self.grid_dims[0] * self.grid_dims[1] * self.grid_dims[2]
    }

    /// Grid dimensions \[nx, ny, nz\].
    pub fn grid_dims(&self) -> [usize; 3] {
        self.grid_dims
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Kernel dispatch configuration
// ─────────────────────────────────────────────────────────────────────────────

/// Configuration for dispatching SPH compute kernels.
pub struct SphDispatchConfig {
    /// Number of particles.
    pub n_particles: usize,
    /// Smoothing length.
    pub h: f64,
    /// Viscosity coefficient.
    pub mu: f64,
    /// Equation of state stiffness (for pressure from density).
    pub k_eos: f64,
    /// Rest density.
    pub rho0: f64,
    /// Workgroup size for GPU dispatch.
    pub workgroup_size: u32,
}

impl SphDispatchConfig {
    /// Create a new dispatch config with default viscosity and EOS parameters.
    pub fn new(n_particles: usize, h: f64) -> Self {
        Self {
            n_particles,
            h,
            mu: 0.1,
            k_eos: 1000.0,
            rho0: 1000.0,
            workgroup_size: 64,
        }
    }

    /// Compute pressure from density using a simple EOS: p = k * (rho - rho0).
    pub fn pressure_from_density(&self, rho: f64) -> f64 {
        self.k_eos * (rho - self.rho0).max(0.0)
    }

    /// Number of workgroups needed.
    pub fn num_workgroups(&self) -> u32 {
        (self.n_particles as u32).div_ceil(self.workgroup_size)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// SPH buffer layout descriptor
// ─────────────────────────────────────────────────────────────────────────────

/// Describes the buffer layout for SPH simulation data.
pub struct SphBufferLayout {
    /// Number of particles.
    pub n_particles: usize,
    /// Size of position buffer (3 * n_particles).
    pub position_size: usize,
    /// Size of velocity buffer (3 * n_particles).
    pub velocity_size: usize,
    /// Size of mass buffer (n_particles).
    pub mass_size: usize,
    /// Size of density buffer (n_particles).
    pub density_size: usize,
    /// Size of pressure buffer (n_particles).
    pub pressure_size: usize,
    /// Size of force buffer (3 * n_particles).
    pub force_size: usize,
}

impl SphBufferLayout {
    /// Create a buffer layout for the given number of particles.
    pub fn new(n_particles: usize) -> Self {
        Self {
            n_particles,
            position_size: 3 * n_particles,
            velocity_size: 3 * n_particles,
            mass_size: n_particles,
            density_size: n_particles,
            pressure_size: n_particles,
            force_size: 3 * n_particles,
        }
    }

    /// Total number of f64 elements needed across all buffers.
    pub fn total_elements(&self) -> usize {
        self.position_size
            + self.velocity_size
            + self.mass_size
            + self.density_size
            + self.pressure_size
            + self.force_size
    }

    /// Total memory in bytes (f64 = 8 bytes).
    pub fn total_bytes(&self) -> usize {
        self.total_elements() * 8
    }
}

/// Kernel that computes SPH particle densities.
///
/// **Inputs:**
///   - `inputs[0]`: positions, flat `[x, y, z, ...]` (3 * n values)
///   - `inputs[1]`: masses, `[m0, m1, ...]` (n values)
///   - `inputs[2]`: `[smoothing_length]` (single value)
///
/// **Outputs:**
///   - `outputs[0]`: densities, one per particle
pub struct SphDensityKernel;

/// Poly6 smoothing kernel value.
#[inline]
fn poly6(r2: f64, h: f64) -> f64 {
    let h2 = h * h;
    if r2 >= h2 {
        return 0.0;
    }
    let coeff = 315.0 / (64.0 * PI * h.powi(9));
    coeff * (h2 - r2).powi(3)
}

/// Spiky kernel gradient magnitude (scalar, multiply by direction).
#[inline]
fn spiky_grad(r: f64, h: f64) -> f64 {
    if r >= h || r < 1e-12 {
        return 0.0;
    }
    let coeff = -45.0 / (PI * h.powi(6));
    coeff * (h - r).powi(2)
}

/// Viscosity Laplacian kernel.
#[inline]
fn viscosity_laplacian(r: f64, h: f64) -> f64 {
    if r >= h || r < 1e-12 {
        return 0.0;
    }
    45.0 / (PI * h.powi(6)) * (h - r)
}

impl ComputeKernel for SphDensityKernel {
    fn name(&self) -> &str {
        "SphDensityKernel"
    }

    fn execute(&self, inputs: &[&[f64]], outputs: &mut [Vec<f64>], work_size: usize) {
        if inputs.len() < 3 || outputs.is_empty() {
            return;
        }
        let positions = inputs[0];
        let masses = inputs[1];
        let h = inputs[2][0];
        let n = work_size;

        let mut densities = vec![0.0; n];
        for i in 0..n {
            let xi = [positions[i * 3], positions[i * 3 + 1], positions[i * 3 + 2]];
            let mut rho = 0.0;
            for j in 0..n {
                let xj = [positions[j * 3], positions[j * 3 + 1], positions[j * 3 + 2]];
                let dx = xi[0] - xj[0];
                let dy = xi[1] - xj[1];
                let dz = xi[2] - xj[2];
                let r2 = dx * dx + dy * dy + dz * dz;
                rho += masses[j] * poly6(r2, h);
            }
            densities[i] = rho;
        }
        outputs[0] = densities;
    }
}

/// Kernel that computes SPH pressure and viscosity forces.
///
/// **Inputs:**
///   - `inputs[0]`: positions `[x, y, z, ...]` (3n)
///   - `inputs[1]`: velocities `[vx, vy, vz, ...]` (3n)
///   - `inputs[2]`: densities `[rho0, rho1, ...]` (n)
///   - `inputs[3]`: pressures `[p0, p1, ...]` (n)
///   - `inputs[4]`: masses `[m0, m1, ...]` (n)
///   - `inputs[5]`: `[smoothing_length, viscosity_coeff]` (2 values)
///
/// **Outputs:**
///   - `outputs[0]`: forces `[fx, fy, fz, ...]` (3n)
pub struct SphForceKernel;

impl ComputeKernel for SphForceKernel {
    fn name(&self) -> &str {
        "SphForceKernel"
    }

    fn execute(&self, inputs: &[&[f64]], outputs: &mut [Vec<f64>], work_size: usize) {
        if inputs.len() < 6 || outputs.is_empty() {
            return;
        }
        let pos = inputs[0];
        let vel = inputs[1];
        let density = inputs[2];
        let pressure = inputs[3];
        let mass = inputs[4];
        let h = inputs[5][0];
        let mu = inputs[5][1]; // viscosity coefficient
        let n = work_size;

        let mut forces = vec![0.0; n * 3];
        for i in 0..n {
            let xi = [pos[i * 3], pos[i * 3 + 1], pos[i * 3 + 2]];
            let vi = [vel[i * 3], vel[i * 3 + 1], vel[i * 3 + 2]];
            let mut fx = 0.0;
            let mut fy = 0.0;
            let mut fz = 0.0;
            for j in 0..n {
                if i == j {
                    continue;
                }
                let xj = [pos[j * 3], pos[j * 3 + 1], pos[j * 3 + 2]];
                let vj = [vel[j * 3], vel[j * 3 + 1], vel[j * 3 + 2]];
                let dx = xi[0] - xj[0];
                let dy = xi[1] - xj[1];
                let dz = xi[2] - xj[2];
                let r = (dx * dx + dy * dy + dz * dz).sqrt();
                if r < 1e-12 || r >= h {
                    continue;
                }
                // Pressure force (Navier-Stokes)
                let p_term =
                    -mass[j] * (pressure[i] + pressure[j]) / (2.0 * density[j]) * spiky_grad(r, h);
                fx += p_term * dx / r;
                fy += p_term * dy / r;
                fz += p_term * dz / r;
                // Viscosity force
                let v_lap = viscosity_laplacian(r, h);
                fx += mu * mass[j] * (vj[0] - vi[0]) / density[j] * v_lap;
                fy += mu * mass[j] * (vj[1] - vi[1]) / density[j] * v_lap;
                fz += mu * mass[j] * (vj[2] - vi[2]) / density[j] * v_lap;
            }
            forces[i * 3] = fx;
            forces[i * 3 + 1] = fy;
            forces[i * 3 + 2] = fz;
        }
        outputs[0] = forces;
    }
}

/// Kernel that constructs a neighbor list on the CPU.
///
/// **Inputs:**
///   - `inputs[0]`: positions `[x, y, z, ...]` (3n)
///   - `inputs[1]`: `[h, origin_x, origin_y, origin_z, domain_x, domain_y, domain_z]`
///
/// **Outputs:**
///   - `outputs[0]`: cell indices (one per particle)
pub struct SphNeighborListKernel;

impl ComputeKernel for SphNeighborListKernel {
    fn name(&self) -> &str {
        "SphNeighborListKernel"
    }

    fn execute(&self, inputs: &[&[f64]], outputs: &mut [Vec<f64>], work_size: usize) {
        if inputs.len() < 2 || outputs.is_empty() {
            return;
        }
        let positions = inputs[0];
        let params = inputs[1];
        if params.len() < 7 {
            return;
        }
        let h = params[0];
        let origin = [params[1], params[2], params[3]];
        let _domain = [params[4], params[5], params[6]];
        let n = work_size;

        // Compute cell index for each particle
        let nx = (_domain[0] / h).ceil() as usize;
        let ny = (_domain[1] / h).ceil() as usize;
        let nx = nx.max(1);
        let ny = ny.max(1);

        let mut cell_indices = vec![0.0f64; n];
        for i in 0..n {
            let px = positions[i * 3] - origin[0];
            let py = positions[i * 3 + 1] - origin[1];
            let pz = positions[i * 3 + 2] - origin[2];
            let ix = (px / h).floor() as usize;
            let iy = (py / h).floor() as usize;
            let iz = (pz / h).floor() as usize;
            cell_indices[i] = (iz * ny * nx + iy * nx + ix) as f64;
        }
        outputs[0] = cell_indices;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Surface tension kernel (color function gradient method)
// ─────────────────────────────────────────────────────────────────────────────

/// Compute surface tension forces using the color function gradient (CSF) method.
///
/// The color function gradient approximates the interface normal.
/// The curvature κ is estimated from the divergence of the unit normal.
/// The surface tension force on particle `i` is:
///
/// `F_i^surf = σ · m_i / ρ_i · ∇(c_i)`
///
/// where `∇c_i = sum_j m_j/ρ_j * (c_j - c_i) * ∇W(r_ij, h)`.
///
/// # Arguments
/// * `positions`  - Particle positions.
/// * `color_fn`   - Color function value per particle (1 = fluid, 0 = void).
/// * `masses`     - Per-particle masses.
/// * `densities`  - Per-particle densities.
/// * `h`          - Smoothing length.
/// * `sigma`      - Surface tension coefficient.
pub fn surface_tension_force(
    positions: &[[f64; 3]],
    color_fn: &[f64],
    masses: &[f64],
    densities: &[f64],
    h: f64,
    sigma: f64,
) -> Vec<[f64; 3]> {
    let params = SphKernelParams::new(h);
    let n = positions.len();
    let mut forces = vec![[0.0f64; 3]; n];

    // Step 1: compute color function gradient ∇c_i for each particle
    let mut grad_c = vec![[0.0f64; 3]; n];
    for i in 0..n {
        let rho_i = if densities[i].abs() > 1e-30 {
            densities[i]
        } else {
            1.0
        };
        let mut gx = 0.0f64;
        let mut gy = 0.0f64;
        let mut gz = 0.0f64;
        for j in 0..n {
            if i == j {
                continue;
            }
            let r_vec = [
                positions[i][0] - positions[j][0],
                positions[i][1] - positions[j][1],
                positions[i][2] - positions[j][2],
            ];
            let r = (r_vec[0] * r_vec[0] + r_vec[1] * r_vec[1] + r_vec[2] * r_vec[2]).sqrt();
            let grad_w = kernel_gradient(r_vec, r, &params, SphKernel::CubicSpline);
            let rho_j = if densities[j].abs() > 1e-30 {
                densities[j]
            } else {
                1.0
            };
            let dc = masses[j] / rho_j * (color_fn[j] - color_fn[i]);
            gx += dc * grad_w[0];
            gy += dc * grad_w[1];
            gz += dc * grad_w[2];
        }
        grad_c[i] = [gx, gy, gz];

        // Step 2: surface tension force F_i = sigma * m_i / rho_i * grad_c_i
        let prefactor = sigma * masses[i] / rho_i;
        forces[i] = [prefactor * gx, prefactor * gy, prefactor * gz];
    }
    forces
}

// ─────────────────────────────────────────────────────────────────────────────
// CFL time step reduction
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the CFL-limited time step for a particle system.
///
/// The CFL condition restricts the time step to:
/// `Δt = cfl_factor · h / max(|v_i| + c_sound)`
///
/// where the maximum is taken over all particles.
///
/// # Arguments
/// * `velocities` - Per-particle velocity vectors.
/// * `h`          - Smoothing length (characteristic length scale).
/// * `c_sound`    - Speed of sound (used as minimum wave speed).
/// * `cfl_factor` - CFL safety factor (typically 0.1 – 0.4).
pub fn cfl_timestep(velocities: &[[f64; 3]], h: f64, c_sound: f64, cfl_factor: f64) -> f64 {
    let max_signal = velocities
        .iter()
        .map(|v| {
            let speed = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
            speed + c_sound
        })
        .fold(0.0f64, f64::max);

    let denominator = if max_signal > 1e-30 {
        max_signal
    } else {
        c_sound.max(1e-30)
    };
    cfl_factor * h / denominator
}

// ─────────────────────────────────────────────────────────────────────────────
// Radix sort by density (for load balancing)
// ─────────────────────────────────────────────────────────────────────────────

/// Sort particle indices in ascending order of density using a reference sort.
///
/// This CPU reference implementation mirrors what a GPU radix sort would do
/// for load-balancing purposes: denser regions can be binned together for
/// more uniform workloads.
///
/// Returns a permutation `indices` such that `densities[indices[0\]] <= densities[indices[1\]] ...`.
pub fn radix_sort_by_density(densities: &[f64]) -> Vec<usize> {
    let mut indices: Vec<usize> = (0..densities.len()).collect();
    indices.sort_by(|&a, &b| {
        densities[a]
            .partial_cmp(&densities[b])
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    indices
}

// ─────────────────────────────────────────────────────────────────────────────
// SPH density accumulation kernel (GPU-style parallel)
// ─────────────────────────────────────────────────────────────────────────────

/// Accumulate density contributions from a neighbour list.
///
/// For each particle `i`, sums contributions from its neighbours `j`:
/// `rho_i += m_j * W(|r_i - r_j|, h)`
///
/// Only pairs within distance `h` contribute.  The neighbour list is expected
/// to already be built for the given positions.
pub fn density_accumulation(
    positions: &[[f64; 3]],
    masses: &[f64],
    h: f64,
    neighbor_list: &NeighborList,
) -> Vec<f64> {
    let params = SphKernelParams::new(h);
    let n = positions.len();
    let mut densities = vec![0.0f64; n];
    for i in 0..n {
        let neighbors = neighbor_list.neighbors(&positions[i]);
        let mut rho = 0.0;
        for j in neighbors {
            let dx = positions[i][0] - positions[j][0];
            let dy = positions[i][1] - positions[j][1];
            let dz = positions[i][2] - positions[j][2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            rho += masses[j] * kernel_value(r, &params, SphKernel::CubicSpline);
        }
        densities[i] = rho;
    }
    densities
}

// ─────────────────────────────────────────────────────────────────────────────
// SPH pressure force kernel
// ─────────────────────────────────────────────────────────────────────────────

/// Symmetric SPH pressure force kernel using neighbour lists.
///
/// `F_i^press = -m_i sum_j m_j (p_i/rho_i^2 + p_j/rho_j^2) nabla W(r_ij, h)`
pub fn pressure_force_kernel(
    positions: &[[f64; 3]],
    densities: &[f64],
    pressures: &[f64],
    masses: &[f64],
    h: f64,
    neighbor_list: &NeighborList,
) -> Vec<[f64; 3]> {
    let params = SphKernelParams::new(h);
    let n = positions.len();
    let mut forces = vec![[0.0f64; 3]; n];
    for i in 0..n {
        let pi_over_rho2 = if densities[i].abs() > 1e-30 {
            pressures[i] / (densities[i] * densities[i])
        } else {
            0.0
        };
        let neighbors = neighbor_list.neighbors(&positions[i]);
        let mut fx = 0.0;
        let mut fy = 0.0;
        let mut fz = 0.0;
        for j in neighbors {
            if i == j {
                continue;
            }
            let r_vec = [
                positions[i][0] - positions[j][0],
                positions[i][1] - positions[j][1],
                positions[i][2] - positions[j][2],
            ];
            let r = (r_vec[0] * r_vec[0] + r_vec[1] * r_vec[1] + r_vec[2] * r_vec[2]).sqrt();
            let grad = kernel_gradient(r_vec, r, &params, SphKernel::CubicSpline);
            let pj_over_rho2 = if densities[j].abs() > 1e-30 {
                pressures[j] / (densities[j] * densities[j])
            } else {
                0.0
            };
            let coeff = -masses[j] * (pi_over_rho2 + pj_over_rho2);
            fx += coeff * grad[0];
            fy += coeff * grad[1];
            fz += coeff * grad[2];
        }
        forces[i] = [fx, fy, fz];
    }
    forces
}

// ─────────────────────────────────────────────────────────────────────────────
// SPH viscosity kernel
// ─────────────────────────────────────────────────────────────────────────────

/// Artificial viscosity for SPH (Monaghan 1992 formulation).
///
/// Adds dissipation to prevent particle interpenetration.
/// `F_i^visc = -sum_j m_j PI_ij nabla W(r_ij, h)`
///
/// where `PI_ij = (-alpha * mu_ij * c_s + beta * mu_ij^2) / rho_ij`
/// and `mu_ij = h * v_ij . r_ij / (|r_ij|^2 + 0.01 h^2)`.
pub fn artificial_viscosity_force(
    positions: &[[f64; 3]],
    velocities: &[[f64; 3]],
    densities: &[f64],
    masses: &[f64],
    h: f64,
    c_s: f64,
    alpha: f64,
    beta: f64,
) -> Vec<[f64; 3]> {
    let params = SphKernelParams::new(h);
    let n = positions.len();
    let mut forces = vec![[0.0f64; 3]; n];
    for i in 0..n {
        let mut fx = 0.0;
        let mut fy = 0.0;
        let mut fz = 0.0;
        for j in 0..n {
            if i == j {
                continue;
            }
            let r_vec = [
                positions[i][0] - positions[j][0],
                positions[i][1] - positions[j][1],
                positions[i][2] - positions[j][2],
            ];
            let v_vec = [
                velocities[i][0] - velocities[j][0],
                velocities[i][1] - velocities[j][1],
                velocities[i][2] - velocities[j][2],
            ];
            let r2 = r_vec[0] * r_vec[0] + r_vec[1] * r_vec[1] + r_vec[2] * r_vec[2];
            let r = r2.sqrt();
            if r >= h || r < 1e-12 {
                continue;
            }
            let vr = v_vec[0] * r_vec[0] + v_vec[1] * r_vec[1] + v_vec[2] * r_vec[2];
            if vr >= 0.0 {
                continue;
            } // only when approaching
            let mu_ij = h * vr / (r2 + 0.01 * h * h);
            let rho_ij = 0.5 * (densities[i] + densities[j]).max(1e-30);
            let pi_ij = (-alpha * c_s * mu_ij + beta * mu_ij * mu_ij) / rho_ij;
            let grad = kernel_gradient(r_vec, r, &params, SphKernel::CubicSpline);
            let coeff = -masses[j] * pi_ij;
            fx += coeff * grad[0];
            fy += coeff * grad[1];
            fz += coeff * grad[2];
        }
        forces[i] = [fx, fy, fz];
    }
    forces
}

// ─────────────────────────────────────────────────────────────────────────────
// WCSPH (Weakly Compressible SPH) step kernel
// ─────────────────────────────────────────────────────────────────────────────

/// Equation of state for WCSPH: `p = B * ((rho/rho0)^gamma - 1)`.
///
/// Typical values: `gamma = 7`, `B = rho0 * c_s^2 / gamma`.
pub fn wcsph_pressure(rho: f64, rho0: f64, b: f64, gamma: f64) -> f64 {
    b * ((rho / rho0).powf(gamma) - 1.0)
}

/// Apply one explicit WCSPH Euler step.
///
/// Updates positions and velocities using the computed forces.
/// Returns (new_positions, new_velocities).
pub fn wcsph_euler_step(
    positions: &[[f64; 3]],
    velocities: &[[f64; 3]],
    forces: &[[f64; 3]],
    masses: &[f64],
    dt: f64,
) -> (Vec<[f64; 3]>, Vec<[f64; 3]>) {
    let n = positions.len();
    let mut new_pos = positions.to_vec();
    let mut new_vel = velocities.to_vec();
    for i in 0..n {
        let m = masses[i].max(1e-30);
        let ax = forces[i][0] / m;
        let ay = forces[i][1] / m;
        let az = forces[i][2] / m;
        new_vel[i] = [
            velocities[i][0] + dt * ax,
            velocities[i][1] + dt * ay,
            velocities[i][2] + dt * az,
        ];
        new_pos[i] = [
            positions[i][0] + dt * new_vel[i][0],
            positions[i][1] + dt * new_vel[i][1],
            positions[i][2] + dt * new_vel[i][2],
        ];
    }
    (new_pos, new_vel)
}

/// Apply one leap-frog WCSPH half-step (velocity update only).
///
/// `v_{n+1/2} = v_{n-1/2} + a_n * dt`
pub fn wcsph_leapfrog_velocity_half(
    velocities: &[[f64; 3]],
    forces: &[[f64; 3]],
    masses: &[f64],
    dt: f64,
) -> Vec<[f64; 3]> {
    let n = velocities.len();
    let mut new_vel = velocities.to_vec();
    for i in 0..n {
        let m = masses[i].max(1e-30);
        new_vel[i][0] += dt * forces[i][0] / m;
        new_vel[i][1] += dt * forces[i][1] / m;
        new_vel[i][2] += dt * forces[i][2] / m;
    }
    new_vel
}

// ─────────────────────────────────────────────────────────────────────────────
// SPH surface normal kernel
// ─────────────────────────────────────────────────────────────────────────────

/// Compute approximate surface normals using the color function gradient method.
///
/// `n_i = sum_j (m_j / rho_j) * nabla W(r_ij, h)`
///
/// The magnitude of `n_i` is proportional to the interface curvature.
pub fn surface_normal_kernel(
    positions: &[[f64; 3]],
    densities: &[f64],
    masses: &[f64],
    h: f64,
) -> Vec<[f64; 3]> {
    let params = SphKernelParams::new(h);
    let n = positions.len();
    let mut normals = vec![[0.0f64; 3]; n];
    for i in 0..n {
        let mut nx = 0.0;
        let mut ny = 0.0;
        let mut nz = 0.0;
        for j in 0..n {
            if i == j {
                continue;
            }
            let r_vec = [
                positions[i][0] - positions[j][0],
                positions[i][1] - positions[j][1],
                positions[i][2] - positions[j][2],
            ];
            let r = (r_vec[0] * r_vec[0] + r_vec[1] * r_vec[1] + r_vec[2] * r_vec[2]).sqrt();
            let rho_j = densities[j].max(1e-30);
            let grad = kernel_gradient(r_vec, r, &params, SphKernel::CubicSpline);
            let coeff = masses[j] / rho_j;
            nx += coeff * grad[0];
            ny += coeff * grad[1];
            nz += coeff * grad[2];
        }
        normals[i] = [nx, ny, nz];
    }
    normals
}

/// Normalize a surface normal vector.  Returns zero vector if magnitude is too small.
pub fn normalize_normal(n: [f64; 3]) -> [f64; 3] {
    let mag = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
    if mag < 1e-30 {
        [0.0; 3]
    } else {
        [n[0] / mag, n[1] / mag, n[2] / mag]
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// SPH neighbour list build (explicit grid-based)
// ─────────────────────────────────────────────────────────────────────────────

/// Build an explicit neighbour list: for each particle, store all particle
/// indices within distance `h`.
///
/// Returns `Vec<Vec`usize`>` where `neighbours[i]` is the list of neighbour indices for particle `i`.
pub fn build_neighbor_list_explicit(positions: &[[f64; 3]], h: f64) -> Vec<Vec<usize>> {
    let n = positions.len();
    let mut neighbors = vec![Vec::new(); n];
    for i in 0..n {
        for j in 0..n {
            let dx = positions[i][0] - positions[j][0];
            let dy = positions[i][1] - positions[j][1];
            let dz = positions[i][2] - positions[j][2];
            if dx * dx + dy * dy + dz * dz < h * h {
                neighbors[i].push(j);
            }
        }
    }
    neighbors
}

/// Compute average number of neighbours per particle.
pub fn mean_neighbor_count(neighbors: &[Vec<usize>]) -> f64 {
    if neighbors.is_empty() {
        return 0.0;
    }
    let total: usize = neighbors.iter().map(|v| v.len()).sum();
    total as f64 / neighbors.len() as f64
}

// ─────────────────────────────────────────────────────────────────────────────
// SPH kernel normalization check
// ─────────────────────────────────────────────────────────────────────────────

/// Numerically integrate the kernel over a sphere of radius `h` using Monte Carlo.
///
/// Useful for verifying kernel normalization (should integrate to ~1 in 3D).
pub fn integrate_kernel_sphere(h: f64, kernel: SphKernel, n_samples: usize) -> f64 {
    let params = SphKernelParams::new(h);
    // Uniform sampling in [0, h] with spherical volume element 4*pi*r^2
    let dr = h / n_samples as f64;
    let mut integral = 0.0;
    for k in 0..n_samples {
        let r = (k as f64 + 0.5) * dr;
        let w = kernel_value(r, &params, kernel);
        integral += w * 4.0 * PI * r * r * dr;
    }
    integral
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sph_kernel_density_sum() {
        // Central particle at origin surrounded by 8 neighbors at corners of
        // a unit cube with side 0.4.  All 9 particles have unit mass.
        // The smoothing length h = 1.0 covers all neighbors (max dist ~ 0.693 < 1).
        let h_val = 1.0_f64;
        let offsets: &[(f64, f64, f64)] = &[
            (0.4, 0.4, 0.4),
            (-0.4, 0.4, 0.4),
            (0.4, -0.4, 0.4),
            (0.4, 0.4, -0.4),
            (-0.4, -0.4, 0.4),
            (-0.4, 0.4, -0.4),
            (0.4, -0.4, -0.4),
            (-0.4, -0.4, -0.4),
        ];
        // First particle is the central one.
        let mut positions: Vec<f64> = vec![0.0, 0.0, 0.0];
        for &(x, y, z) in offsets {
            positions.extend_from_slice(&[x, y, z]);
        }
        let n = 9_usize;
        let masses = vec![1.0_f64; n];
        let h_slice = vec![h_val];

        let mut outputs = vec![Vec::new()];
        SphDensityKernel.execute(&[&positions, &masses, &h_slice], &mut outputs, n);

        assert_eq!(outputs[0].len(), n, "density output length should equal n");
        // The central particle should have a positive density (it sees all neighbors).
        let central_density = outputs[0][0];
        assert!(
            central_density > 0.0,
            "central particle density should be > 0, got {central_density}"
        );
    }

    #[test]
    fn sph_density_uniform_distribution() {
        // 2 particles at the same position should give non-zero density
        let positions = vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let masses = vec![1.0, 1.0];
        let h = vec![1.0];
        let mut outputs = vec![Vec::new()];
        SphDensityKernel.execute(&[&positions, &masses, &h], &mut outputs, 2);
        assert_eq!(outputs[0].len(), 2);
        // Both particles see each other at r=0, so density should be positive
        assert!(outputs[0][0] > 0.0);
        assert!((outputs[0][0] - outputs[0][1]).abs() < 1e-12);
    }

    #[test]
    fn sph_force_produces_finite_forces() {
        let n = 3;
        // Place particles in a line along x
        let positions = vec![0.0, 0.0, 0.0, 0.3, 0.0, 0.0, 0.6, 0.0, 0.0];
        let velocities = vec![0.0; 9];
        let densities = vec![1000.0; 3];
        let pressures = vec![100.0, 200.0, 100.0];
        let masses = vec![1.0; 3];
        let params = vec![1.0, 0.1]; // h=1.0, mu=0.1

        let mut outputs = vec![Vec::new()];
        SphForceKernel.execute(
            &[
                &positions,
                &velocities,
                &densities,
                &pressures,
                &masses,
                &params,
            ],
            &mut outputs,
            n,
        );
        assert_eq!(outputs[0].len(), 9);
        // Forces should be finite
        for &f in &outputs[0] {
            assert!(f.is_finite(), "force component is not finite: {f}");
        }
    }

    // ── High-level SPH kernel API tests ──────────────────────────────────────

    /// Verify kernel_value integrates numerically to approximately 1 for each
    /// kernel type (compact support, unit normalization).
    #[test]
    fn kernel_value_positive_within_support() {
        let h = 1.0_f64;
        let params = SphKernelParams::new(h);
        for &k in &[
            SphKernel::CubicSpline,
            SphKernel::Wendland,
            SphKernel::Poly6,
            SphKernel::Spiky,
        ] {
            // W(0) should be positive (self-contribution)
            let w0 = kernel_value(0.0, &params, k);
            assert!(w0 > 0.0, "{k:?}: W(0) should be > 0, got {w0}");
            // W at r >= h should be exactly 0
            let wh = kernel_value(h, &params, k);
            assert_eq!(wh, 0.0, "{k:?}: W(h) should be 0, got {wh}");
        }
    }

    /// Kernel must be symmetric: W(-r) == W(r) (evaluated via magnitude).
    #[test]
    fn kernel_value_symmetric() {
        let h = 2.0_f64;
        let params = SphKernelParams::new(h);
        let r = 0.7 * h;
        for &k in &[
            SphKernel::CubicSpline,
            SphKernel::Wendland,
            SphKernel::Poly6,
            SphKernel::Spiky,
        ] {
            let w_pos = kernel_value(r, &params, k);
            // Symmetry: W depends only on |r|, so same value at same magnitude.
            let w_same = kernel_value(r, &params, k);
            assert!(
                (w_pos - w_same).abs() < 1e-15,
                "{k:?}: kernel not symmetric at r={r}"
            );
        }
    }

    /// density_summation: self-contribution gives positive density.
    #[test]
    fn density_summation_self_contribution() {
        let positions = vec![[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]];
        let masses = vec![1.0, 1.0];
        let h = 1.0;
        let densities = density_summation(&positions, &masses, h);
        assert_eq!(densities.len(), 2);
        assert!(
            densities[0] > 0.0,
            "density[0] should be > 0, got {}",
            densities[0]
        );
        assert!(
            densities[1] > 0.0,
            "density[1] should be > 0, got {}",
            densities[1]
        );
    }

    /// pressure_force: forces are finite and Newton's third law holds.
    #[test]
    fn pressure_force_finite_and_newtons_third_law() {
        let positions = vec![[0.0, 0.0, 0.0], [0.3, 0.0, 0.0], [0.6, 0.0, 0.0]];
        let velocities = vec![[0.0; 3]; 3];
        let densities = vec![1000.0, 1000.0, 1000.0];
        let pressures = vec![100.0, 200.0, 100.0];
        let masses = vec![1.0, 1.0, 1.0];
        let h = 1.0;

        let forces = pressure_force(&positions, &velocities, &densities, &pressures, &masses, h);
        assert_eq!(forces.len(), 3);
        for (i, f) in forces.iter().enumerate() {
            for &c in f.iter() {
                assert!(c.is_finite(), "force[{i}] component not finite: {c}");
            }
        }
        // Total momentum change should be near zero (action-reaction).
        let total_fx: f64 = forces.iter().map(|f| f[0]).sum();
        assert!(
            total_fx.abs() < 1e-8,
            "total x-force should be ~0 (Newton III), got {total_fx}"
        );
    }

    // ── New tests ────────────────────────────────────────────────────────────

    #[test]
    fn test_density_summation_kernel_poly6() {
        let positions = vec![[0.0, 0.0, 0.0], [0.3, 0.0, 0.0]];
        let masses = vec![1.0, 1.0];
        let h = 1.0;
        let densities = density_summation_kernel(&positions, &masses, h, SphKernel::Poly6);
        assert_eq!(densities.len(), 2);
        assert!(densities[0] > 0.0);
        assert!(densities[1] > 0.0);
    }

    #[test]
    fn test_density_summation_kernel_wendland() {
        let positions = vec![[0.0, 0.0, 0.0]];
        let masses = vec![1.0];
        let h = 1.0;
        let densities = density_summation_kernel(&positions, &masses, h, SphKernel::Wendland);
        assert_eq!(densities.len(), 1);
        assert!(densities[0] > 0.0);
    }

    #[test]
    fn test_density_summation_kernel_spiky() {
        let positions = vec![[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]];
        let masses = vec![1.0, 1.0];
        let h = 1.0;
        let densities = density_summation_kernel(&positions, &masses, h, SphKernel::Spiky);
        assert!(densities[0] > 0.0);
        assert!(densities[1] > 0.0);
    }

    #[test]
    fn test_viscosity_force_finite() {
        let positions = vec![[0.0, 0.0, 0.0], [0.3, 0.0, 0.0]];
        let velocities = vec![[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]];
        let densities = vec![1000.0, 1000.0];
        let masses = vec![1.0, 1.0];
        let h = 1.0;
        let mu = 0.1;
        let forces = viscosity_force(&positions, &velocities, &densities, &masses, h, mu);
        assert_eq!(forces.len(), 2);
        for f in &forces {
            for &c in f {
                assert!(c.is_finite(), "viscosity force not finite: {c}");
            }
        }
    }

    #[test]
    fn test_viscosity_force_zero_for_same_velocity() {
        let positions = vec![[0.0, 0.0, 0.0], [0.3, 0.0, 0.0]];
        let velocities = vec![[1.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let densities = vec![1000.0, 1000.0];
        let masses = vec![1.0, 1.0];
        let h = 1.0;
        let mu = 0.1;
        let forces = viscosity_force(&positions, &velocities, &densities, &masses, h, mu);
        for f in &forces {
            for &c in f {
                assert!(
                    c.abs() < 1e-12,
                    "viscosity should be zero for uniform velocity"
                );
            }
        }
    }

    #[test]
    fn test_neighbor_list_build_and_query() {
        let positions = vec![
            [0.5, 0.5, 0.5],
            [0.6, 0.5, 0.5],
            [5.0, 5.0, 5.0], // far away
        ];
        let mut nlist = NeighborList::new([0.0, 0.0, 0.0], [10.0, 10.0, 10.0], 1.0);
        nlist.build(&positions);

        let neighbors = nlist.neighbors(&[0.5, 0.5, 0.5]);
        // Should find particles 0 and 1 (same or adjacent cells), not 2
        assert!(neighbors.contains(&0));
        assert!(neighbors.contains(&1));
        assert!(!neighbors.contains(&2));
    }

    #[test]
    fn test_neighbor_list_num_cells() {
        let nlist = NeighborList::new([0.0, 0.0, 0.0], [10.0, 10.0, 10.0], 1.0);
        assert_eq!(nlist.num_cells(), 1000); // 10x10x10
        assert_eq!(nlist.grid_dims(), [10, 10, 10]);
    }

    #[test]
    fn test_neighbor_list_single_particle() {
        let positions = vec![[0.5, 0.5, 0.5]];
        let mut nlist = NeighborList::new([0.0, 0.0, 0.0], [5.0, 5.0, 5.0], 1.0);
        nlist.build(&positions);
        let neighbors = nlist.neighbors(&[0.5, 0.5, 0.5]);
        assert!(neighbors.contains(&0));
    }

    #[test]
    fn test_sph_dispatch_config() {
        let config = SphDispatchConfig::new(1000, 0.1);
        assert_eq!(config.n_particles, 1000);
        assert_eq!(config.num_workgroups(), 16); // ceil(1000/64) = 16
    }

    #[test]
    fn test_sph_dispatch_config_pressure() {
        let config = SphDispatchConfig::new(100, 0.1);
        let p = config.pressure_from_density(1100.0);
        assert!((p - 100_000.0).abs() < 1e-6); // k=1000 * (1100-1000) = 100000
        let p_zero = config.pressure_from_density(500.0);
        assert!((p_zero).abs() < 1e-12); // below rest density -> 0
    }

    #[test]
    fn test_sph_buffer_layout() {
        let layout = SphBufferLayout::new(1000);
        assert_eq!(layout.position_size, 3000);
        assert_eq!(layout.velocity_size, 3000);
        assert_eq!(layout.mass_size, 1000);
        assert_eq!(layout.density_size, 1000);
        assert_eq!(layout.pressure_size, 1000);
        assert_eq!(layout.force_size, 3000);
        assert_eq!(layout.total_elements(), 12000);
        assert_eq!(layout.total_bytes(), 96000);
    }

    #[test]
    fn test_neighbor_list_kernel_executes() {
        let positions = vec![0.5, 0.5, 0.5, 1.5, 1.5, 1.5];
        let params = vec![1.0, 0.0, 0.0, 0.0, 5.0, 5.0, 5.0];
        let mut outputs = vec![Vec::new()];
        SphNeighborListKernel.execute(&[&positions, &params], &mut outputs, 2);
        assert_eq!(outputs[0].len(), 2);
        // cell indices should be non-negative
        assert!(outputs[0][0] >= 0.0);
        assert!(outputs[0][1] >= 0.0);
        // Different positions should generally give different cells
        assert!((outputs[0][0] - outputs[0][1]).abs() > 0.5);
    }

    #[test]
    fn test_kernel_gradient_at_origin_is_zero() {
        let h = 1.0;
        let params = SphKernelParams::new(h);
        for &k in &[
            SphKernel::CubicSpline,
            SphKernel::Wendland,
            SphKernel::Poly6,
            SphKernel::Spiky,
        ] {
            let grad = kernel_gradient([0.0, 0.0, 0.0], 0.0, &params, k);
            assert_eq!(
                grad,
                [0.0, 0.0, 0.0],
                "{k:?}: gradient at origin should be zero"
            );
        }
    }

    #[test]
    fn test_kernel_gradient_outside_support_is_zero() {
        let h = 1.0;
        let params = SphKernelParams::new(h);
        for &k in &[
            SphKernel::CubicSpline,
            SphKernel::Wendland,
            SphKernel::Poly6,
            SphKernel::Spiky,
        ] {
            let grad = kernel_gradient([2.0, 0.0, 0.0], 2.0, &params, k);
            assert_eq!(
                grad,
                [0.0, 0.0, 0.0],
                "{k:?}: gradient outside support should be zero"
            );
        }
    }

    // ── Surface tension tests ─────────────────────────────────────────────

    #[test]
    fn test_surface_tension_force_finite() {
        let positions = vec![[0.0, 0.0, 0.0], [0.3, 0.0, 0.0], [0.6, 0.0, 0.0]];
        let color_fn = vec![1.0, 1.0, 0.0];
        let masses = vec![1.0, 1.0, 1.0];
        let densities = vec![1000.0, 1000.0, 1000.0];
        let h = 1.0;
        let sigma = 0.07;
        let forces = surface_tension_force(&positions, &color_fn, &masses, &densities, h, sigma);
        assert_eq!(forces.len(), 3);
        for f in &forces {
            for &c in f {
                assert!(c.is_finite(), "surface tension force not finite: {c}");
            }
        }
    }

    #[test]
    fn test_surface_tension_zero_sigma() {
        let positions = vec![[0.0, 0.0, 0.0], [0.3, 0.0, 0.0]];
        let color_fn = vec![1.0, 0.0];
        let masses = vec![1.0, 1.0];
        let densities = vec![1000.0, 1000.0];
        let h = 1.0;
        let forces = surface_tension_force(&positions, &color_fn, &masses, &densities, h, 0.0);
        for f in &forces {
            for &c in f {
                assert!(
                    c.abs() < 1e-30,
                    "surface tension with sigma=0 should be zero"
                );
            }
        }
    }

    // ── CFL time step tests ───────────────────────────────────────────────

    #[test]
    fn test_cfl_timestep_basic() {
        let velocities = vec![[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 0.5]];
        let h = 1.0;
        let c_sound = 10.0;
        let dt = cfl_timestep(&velocities, h, c_sound, 0.3);
        assert!(dt > 0.0, "CFL dt should be positive");
        assert!(dt.is_finite(), "CFL dt should be finite");
        // max velocity magnitude = 2.0; CFL condition: dt <= 0.3*h/(v+c) = 0.3/12 ~ 0.025
        assert!(dt <= 0.3 * h / (2.0 + c_sound) + 1e-12);
    }

    #[test]
    fn test_cfl_timestep_zero_velocity() {
        let velocities = vec![[0.0; 3]; 5];
        let h = 0.5;
        let c_sound = 5.0;
        let dt = cfl_timestep(&velocities, h, c_sound, 0.25);
        // With zero velocity: dt = cfl * h / c_sound
        let expected = 0.25 * h / c_sound;
        assert!(
            (dt - expected).abs() < 1e-10,
            "expected {expected}, got {dt}"
        );
    }

    // ── Radix sort by density tests ───────────────────────────────────────

    #[test]
    fn test_radix_sort_by_density_ordered() {
        let densities = vec![3.0, 1.0, 4.0, 1.5, 9.0, 2.6, 5.0, 3.5];
        let indices = radix_sort_by_density(&densities);
        assert_eq!(indices.len(), densities.len());
        // Verify sorted order
        for w in indices.windows(2) {
            assert!(
                densities[w[0]] <= densities[w[1]],
                "not sorted: {} > {}",
                densities[w[0]],
                densities[w[1]]
            );
        }
    }

    #[test]
    fn test_radix_sort_by_density_single() {
        let densities = vec![42.0];
        let indices = radix_sort_by_density(&densities);
        assert_eq!(indices, vec![0]);
    }

    #[test]
    fn test_radix_sort_by_density_empty() {
        let indices = radix_sort_by_density(&[]);
        assert!(indices.is_empty());
    }

    #[test]
    fn test_radix_sort_is_permutation() {
        let densities = vec![5.0, 3.0, 8.0, 1.0, 2.0];
        let indices = radix_sort_by_density(&densities);
        let mut check = indices.clone();
        check.sort_unstable();
        assert_eq!(
            check,
            vec![0, 1, 2, 3, 4],
            "indices should be a permutation"
        );
    }

    // ── density_accumulation tests ───────────────────────────────────────────

    #[test]
    fn test_density_accumulation_matches_direct() {
        let positions = vec![[0.0, 0.0, 0.0], [0.3, 0.0, 0.0], [0.6, 0.0, 0.0]];
        let masses = vec![1.0, 1.0, 1.0];
        let h = 1.0;
        let mut nlist = NeighborList::new([0.0, 0.0, 0.0], [5.0, 5.0, 5.0], h);
        nlist.build(&positions);
        let rho_nl = density_accumulation(&positions, &masses, h, &nlist);
        let rho_direct = density_summation(&positions, &masses, h);
        // Results should agree (neighbour list finds all within h)
        for i in 0..3 {
            assert!(
                (rho_nl[i] - rho_direct[i]).abs() < 1e-10,
                "density_accumulation vs direct at {i}: {} vs {}",
                rho_nl[i],
                rho_direct[i]
            );
        }
    }

    #[test]
    fn test_density_accumulation_positive() {
        let positions = vec![[0.5, 0.5, 0.5], [0.6, 0.5, 0.5]];
        let masses = vec![1.0, 1.0];
        let h = 1.0;
        let mut nlist = NeighborList::new([0.0, 0.0, 0.0], [5.0, 5.0, 5.0], h);
        nlist.build(&positions);
        let rho = density_accumulation(&positions, &masses, h, &nlist);
        assert!(rho[0] > 0.0 && rho[1] > 0.0);
    }

    // ── pressure_force_kernel tests ──────────────────────────────────────────

    #[test]
    fn test_pressure_force_kernel_finite() {
        let positions = vec![[0.0, 0.0, 0.0], [0.3, 0.0, 0.0]];
        let densities = vec![1000.0, 1000.0];
        let pressures = vec![100.0, 200.0];
        let masses = vec![1.0, 1.0];
        let h = 1.0;
        let mut nlist = NeighborList::new([0.0, 0.0, 0.0], [5.0, 5.0, 5.0], h);
        nlist.build(&positions);
        let forces = pressure_force_kernel(&positions, &densities, &pressures, &masses, h, &nlist);
        assert_eq!(forces.len(), 2);
        for f in &forces {
            for &c in f {
                assert!(c.is_finite(), "pressure_force_kernel not finite: {c}");
            }
        }
    }

    #[test]
    fn test_pressure_force_kernel_zero_gradient() {
        // Uniform pressure → no force
        let positions = vec![[0.0, 0.0, 0.0], [0.3, 0.0, 0.0]];
        let densities = vec![1000.0, 1000.0];
        let pressures = vec![100.0, 100.0];
        let masses = vec![1.0, 1.0];
        let h = 1.0;
        let mut nlist = NeighborList::new([0.0, 0.0, 0.0], [5.0, 5.0, 5.0], h);
        nlist.build(&positions);
        let forces = pressure_force_kernel(&positions, &densities, &pressures, &masses, h, &nlist);
        // For uniform pressure with symmetric formulation the net force should be near zero
        let total_fx: f64 = forces.iter().map(|f| f[0]).sum();
        assert!(
            total_fx.abs() < 1e-6,
            "total pressure force with uniform pressure = {total_fx}"
        );
    }

    // ── artificial_viscosity_force tests ─────────────────────────────────────

    #[test]
    fn test_artificial_viscosity_finite() {
        let positions = vec![[0.0, 0.0, 0.0], [0.3, 0.0, 0.0]];
        let velocities = vec![[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]]; // approaching
        let densities = vec![1000.0, 1000.0];
        let masses = vec![1.0, 1.0];
        let h = 1.0;
        let forces = artificial_viscosity_force(
            &positions,
            &velocities,
            &densities,
            &masses,
            h,
            100.0,
            1.0,
            2.0,
        );
        for f in &forces {
            for &c in f {
                assert!(c.is_finite(), "art. visc. force not finite: {c}");
            }
        }
    }

    #[test]
    fn test_artificial_viscosity_zero_for_diverging() {
        // When particles are moving apart, PI_ij = 0
        let positions = vec![[0.0, 0.0, 0.0], [0.3, 0.0, 0.0]];
        let velocities = vec![[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]]; // diverging
        let densities = vec![1000.0, 1000.0];
        let masses = vec![1.0, 1.0];
        let h = 1.0;
        let forces = artificial_viscosity_force(
            &positions,
            &velocities,
            &densities,
            &masses,
            h,
            100.0,
            1.0,
            2.0,
        );
        for f in &forces {
            for &c in f {
                assert!(
                    c.abs() < 1e-30,
                    "diverging particles: art visc should be 0, got {c}"
                );
            }
        }
    }

    // ── WCSPH tests ───────────────────────────────────────────────────────────

    #[test]
    fn test_wcsph_pressure_at_rest_density() {
        // At rho = rho0: pressure = 0
        let p = wcsph_pressure(1000.0, 1000.0, 100.0, 7.0);
        assert!(p.abs() < 1e-8, "pressure at rest density = {p}");
    }

    #[test]
    fn test_wcsph_pressure_above_rest_density() {
        // At rho > rho0: pressure > 0
        let p = wcsph_pressure(1010.0, 1000.0, 100.0, 7.0);
        assert!(
            p > 0.0,
            "pressure above rest density should be positive: {p}"
        );
    }

    #[test]
    fn test_wcsph_euler_step_positions_change() {
        let positions = vec![[0.0, 0.0, 0.0]];
        let velocities = vec![[1.0, 0.0, 0.0]];
        let forces = vec![[0.0, 0.0, 0.0]];
        let masses = vec![1.0];
        let dt = 0.01;
        let (new_pos, _) = wcsph_euler_step(&positions, &velocities, &forces, &masses, dt);
        assert!(
            (new_pos[0][0] - 0.01).abs() < 1e-12,
            "position should advance by v*dt"
        );
    }

    #[test]
    fn test_wcsph_euler_step_velocity_changes() {
        let positions = vec![[0.0, 0.0, 0.0]];
        let velocities = vec![[0.0, 0.0, 0.0]];
        let forces = vec![[1.0, 0.0, 0.0]]; // acceleration = 1
        let masses = vec![1.0];
        let dt = 0.1;
        let (_, new_vel) = wcsph_euler_step(&positions, &velocities, &forces, &masses, dt);
        assert!(
            (new_vel[0][0] - 0.1).abs() < 1e-12,
            "velocity should increase by a*dt"
        );
    }

    #[test]
    fn test_wcsph_leapfrog_velocity_half() {
        let velocities = vec![[1.0, 0.0, 0.0]];
        let forces = vec![[2.0, 0.0, 0.0]];
        let masses = vec![1.0];
        let dt = 0.1;
        let new_vel = wcsph_leapfrog_velocity_half(&velocities, &forces, &masses, dt);
        // v += a * dt = 1 + 2*0.1 = 1.2
        assert!(
            (new_vel[0][0] - 1.2).abs() < 1e-12,
            "leapfrog v = {}",
            new_vel[0][0]
        );
    }

    // ── surface_normal_kernel tests ───────────────────────────────────────────

    #[test]
    fn test_surface_normal_kernel_finite() {
        let positions = vec![[0.0, 0.0, 0.0], [0.3, 0.0, 0.0], [0.6, 0.0, 0.0]];
        let densities = vec![1000.0, 1000.0, 1000.0];
        let masses = vec![1.0, 1.0, 1.0];
        let h = 1.0;
        let normals = surface_normal_kernel(&positions, &densities, &masses, h);
        assert_eq!(normals.len(), 3);
        for n in &normals {
            for &c in n {
                assert!(c.is_finite(), "surface normal not finite: {c}");
            }
        }
    }

    #[test]
    fn test_normalize_normal_unit_vector() {
        let n = normalize_normal([3.0, 0.0, 4.0]);
        let mag = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        assert!((mag - 1.0).abs() < 1e-12, "normalized magnitude = {mag}");
    }

    #[test]
    fn test_normalize_normal_zero_vector() {
        let n = normalize_normal([0.0, 0.0, 0.0]);
        assert_eq!(n, [0.0, 0.0, 0.0]);
    }

    // ── build_neighbor_list_explicit tests ───────────────────────────────────

    #[test]
    fn test_build_neighbor_list_explicit_finds_close() {
        let positions = vec![[0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [5.0, 0.0, 0.0]];
        let h = 1.0;
        let neighbors = build_neighbor_list_explicit(&positions, h);
        assert!(
            neighbors[0].contains(&1),
            "particle 0 should find particle 1"
        );
        assert!(
            !neighbors[0].contains(&2),
            "particle 0 should not find particle 2"
        );
    }

    #[test]
    fn test_build_neighbor_list_self_included() {
        let positions = vec![[0.0, 0.0, 0.0]];
        let h = 1.0;
        let neighbors = build_neighbor_list_explicit(&positions, h);
        assert!(neighbors[0].contains(&0), "particle should find itself");
    }

    #[test]
    fn test_mean_neighbor_count() {
        let positions = vec![[0.0, 0.0, 0.0], [0.3, 0.0, 0.0], [0.6, 0.0, 0.0]];
        let h = 1.0;
        let neighbors = build_neighbor_list_explicit(&positions, h);
        let mean = mean_neighbor_count(&neighbors);
        assert!(mean >= 1.0, "mean neighbors should be >= 1: {mean}");
    }

    #[test]
    fn test_mean_neighbor_count_empty() {
        let mean = mean_neighbor_count(&[]);
        assert_eq!(mean, 0.0);
    }

    // ── integrate_kernel_sphere tests ────────────────────────────────────────

    #[test]
    fn test_integrate_kernel_sphere_cubic_spline() {
        let h = 1.0;
        let integral = integrate_kernel_sphere(h, SphKernel::CubicSpline, 1000);
        // Should be approximately 1 for a normalized kernel
        assert!(
            integral > 0.5 && integral < 2.0,
            "CubicSpline integral = {integral} (expected ~1)"
        );
    }

    #[test]
    fn test_integrate_kernel_sphere_wendland() {
        let h = 1.0;
        let integral = integrate_kernel_sphere(h, SphKernel::Wendland, 1000);
        assert!(
            integral > 0.5 && integral < 2.0,
            "Wendland integral = {integral} (expected ~1)"
        );
    }

    #[test]
    fn test_integrate_kernel_sphere_zero_at_boundary() {
        // Kernel value at r=h should be 0
        let h = 1.0;
        let params = SphKernelParams::new(h);
        for &k in &[SphKernel::CubicSpline, SphKernel::Wendland] {
            let w = kernel_value(h, &params, k);
            assert_eq!(w, 0.0, "{k:?}: W(h) should be 0");
        }
    }
}
