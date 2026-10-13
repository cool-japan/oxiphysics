// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! GPU-accelerated Smoothed Particle Hydrodynamics (SPH) simulation.
//!
//! This module runs a Weakly Compressible SPH (WCSPH) density–pressure solve on
//! the GPU through the real
//! [`WgpuBackendReal`](crate::compute::wgpu_backend::real::WgpuBackendReal)
//! dispatch path, falling back to a clean CPU implementation when no GPU is
//! available.
//!
//! # Physical model
//!
//! Weakly Compressible SPH (WCSPH) with:
//! - **Density**: ρᵢ = Σⱼ mⱼ W(rᵢⱼ, h)   (cubic-spline W3 kernel, support 2h)
//! - **Pressure**: pᵢ = k (ρᵢ/ρ₀ − 1)     (Tait equation of state)
//! - **Acceleration**: aᵢ = −Σⱼ mⱼ (pᵢ/ρᵢ² + pⱼ/ρⱼ²) ∇W + aᵢ^visc + g
//! - **Viscosity**: aᵢ^visc = ν Σⱼ mⱼ/ρⱼ  (r⃗ᵢⱼ · ∇W / (|r⃗ᵢⱼ|² + ε)) (active when approaching)
//!
//! ## GPU dispatch strategy
//!
//! The GPU path is a real resident pipeline on one long-lived device. Each
//! `step` uploads the current particle state, then dispatches the full kernel
//! sequence without any host round-trip in between:
//!
//! 1. **cell-list**: hash each particle into a spatial-hash bucket
//!    (`sph_cell_list` kernel → `cell_keys[]`, `particle_ids[]`).
//! 2. **sort**: stable radix-sort the `(cell_key, particle_id)` pairs on the GPU
//!    (one pass per 8 key bits, with a multi-level GPU scan of the per-tile digit
//!    counts) so each bucket's particles are contiguous in `sorted_ids[]`.
//! 3. **bucket bounds**: `sph_cell_histogram` (global atomic counts, any bucket
//!    count up to 65536) and a GPU exclusive scan yield `cell_count[]` and
//!    `cell_start[]`.
//! 4. **density**: `sph_density` walks the 27 neighbor buckets per particle,
//!    distance-filtering candidates (collision-safe) and summing W.
//! 5. **force**: `sph_force` computes Tait-EOS pressure and accumulates the
//!    symmetric Monaghan pressure force, artificial viscosity, and gravity.
//! 6. **integrate**: `sph_integrate` advances velocity then position (symplectic).
//! 7. **boundary**: `sph_boundary` clamps to the domain AABB and reflects walls.
//!
//! All buffers (particle state, cell list, sort/scan scratch and one parameter
//! buffer per kernel) are allocated once when the simulation is created, and the
//! shader pipelines are compiled once by the first step, so a step costs only the
//! particle-state upload, the dispatches and the final readback; the GPU buffer
//! count is constant for the lifetime of the simulation. Results are read back
//! into the structure-of-arrays [`SphParticleState`] at the end of `step`.
//!
//! The hash grid spacing is the kernel support (2h), so the 27-cell sweep always
//! covers every neighbor. A non-positive or NaN support (or a particle count the
//! 1-D dispatches cannot address) makes the GPU path inapplicable; the exact CPU
//! solve then runs instead.
//!
//! ## Usage
//!
//! ```
//! use oxiphysics_gpu::sph_gpu::{SphSimulation, SphConfig};
//!
//! let cfg = SphConfig { n_particles: 64, smoothing_h: 0.1, rest_density: 1000.0, ..SphConfig::default() };
//! let mut sim = SphSimulation::new(cfg);
//!
//! // Place particles in a 4×4×4 grid
//! for i in 0..4 { for j in 0..4 { for k in 0..4 {
//!     let idx = i * 16 + j * 4 + k;
//!     sim.state.pos_x[idx] = i as f64 * 0.1;
//!     sim.state.pos_y[idx] = j as f64 * 0.1 + 1.0;
//!     sim.state.pos_z[idx] = k as f64 * 0.1;
//! }}}
//!
//! // Simulate 10 frames at 60 Hz
//! for _ in 0..10 { sim.step(1.0 / 60.0); }
//!
//! // Particles should have moved under gravity
//! assert!(sim.state.pos_y[0] < 1.0 + 0.1,
//!     "particles should fall under gravity");
//! ```

// ── SphConfig ─────────────────────────────────────────────────────────────────

/// Configuration for an SPH simulation.
#[derive(Debug, Clone)]
pub struct SphConfig {
    /// Number of particles.
    pub n_particles: usize,
    /// Smoothing length h (m).  Kernel support radius = 2h.
    pub smoothing_h: f64,
    /// Rest density ρ₀ (kg/m³).
    pub rest_density: f64,
    /// Pressure stiffness constant k (Pa).
    pub pressure_k: f64,
    /// Kinematic viscosity ν (m²/s).
    pub viscosity: f64,
    /// Gravitational acceleration (m/s²), applied in −Y direction.
    pub gravity: f64,
    /// Particle mass (kg).  If 0.0, computed as ρ₀ × (2h)³.
    pub particle_mass: f64,
    /// Simulation domain (AABB) minimum corner.
    pub domain_min: [f64; 3],
    /// Simulation domain maximum corner.
    pub domain_max: [f64; 3],
    /// Boundary restitution coefficient [0, 1].
    pub boundary_restitution: f64,
}

impl Default for SphConfig {
    fn default() -> Self {
        let h = 0.05;
        Self {
            n_particles: 256,
            smoothing_h: h,
            rest_density: 1000.0,
            pressure_k: 100.0,
            viscosity: 0.01,
            gravity: 9.81,
            particle_mass: 0.0, // computed below in new()
            domain_min: [-1.0, 0.0, -1.0],
            domain_max: [1.0, 2.0, 1.0],
            boundary_restitution: 0.3,
        }
    }
}

// ── SphParticleState ──────────────────────────────────────────────────────────

/// Structure-of-Arrays particle state for N SPH particles.
#[derive(Debug)]
pub struct SphParticleState {
    /// Number of particles.
    pub n: usize,
    /// X positions (m).
    pub pos_x: Vec<f64>,
    /// Y positions (m).
    pub pos_y: Vec<f64>,
    /// Z positions (m).
    pub pos_z: Vec<f64>,
    /// X velocities (m/s).
    pub vel_x: Vec<f64>,
    /// Y velocities (m/s).
    pub vel_y: Vec<f64>,
    /// Z velocities (m/s).
    pub vel_z: Vec<f64>,
    /// Density (kg/m³).
    pub density: Vec<f64>,
    /// Pressure (Pa).
    pub pressure: Vec<f64>,
}

impl SphParticleState {
    /// Create a zeroed state for `n` particles.
    pub fn new(n: usize) -> Self {
        Self {
            n,
            pos_x: vec![0.0; n],
            pos_y: vec![0.0; n],
            pos_z: vec![0.0; n],
            vel_x: vec![0.0; n],
            vel_y: vec![0.0; n],
            vel_z: vec![0.0; n],
            density: vec![0.0; n],
            pressure: vec![0.0; n],
        }
    }

    /// Reset velocities to zero.
    pub fn zero_velocities(&mut self) {
        self.vel_x.fill(0.0);
        self.vel_y.fill(0.0);
        self.vel_z.fill(0.0);
    }
}

// ── SPH kernel helper ─────────────────────────────────────────────────────────

/// Cubic-spline kernel W3(r, h), normalized for 3-D with support 2h.
///
/// W(r,h) = σ f(q),  q = r/h,  σ = 1/(π h³),
/// f(q) = 1 − 1.5 q² + 0.75 q³ on [0,1), 0.25 (2−q)³ on [1,2), 0 otherwise.
#[inline]
pub fn cubic_spline_w3(r: f64, h: f64) -> f64 {
    let q = r / h;
    let sigma = 1.0 / (std::f64::consts::PI * h * h * h);
    if q < 1.0 {
        sigma * (1.0 - 1.5 * q * q + 0.75 * q * q * q)
    } else if q < 2.0 {
        let t = 2.0 - q;
        sigma * 0.25 * t * t * t
    } else {
        0.0
    }
}

/// Gradient factor of cubic-spline kernel: dW/dr (for r > 0; multiply by r̂).
///
/// σ₄ = 1/(π h⁴); dW/dr = σ₄ (−3q + 2.25 q²) on [0,1), σ₄ (−0.75 (2−q)²) on [1,2).
#[inline]
pub fn cubic_spline_dw_dr(r: f64, h: f64) -> f64 {
    let q = r / h;
    let sigma4 = 1.0 / (std::f64::consts::PI * h * h * h * h);
    if r < 1e-12 || q >= 2.0 {
        0.0
    } else if q < 1.0 {
        sigma4 * (-3.0 * q + 2.25 * q * q)
    } else {
        let t = 2.0 - q;
        sigma4 * (-0.75 * t * t)
    }
}

// ── GPU backend handle ──────────────────────────────────────────────────────

#[cfg(feature = "wgpu-backend")]
use gpu_pipeline::SphGpuState;

/// GPU-built spatial-hash cell-list, exposed for parity testing.
///
/// Produced by [`gpu_cell_list`]: the particle `(cell_key, index)` pairs are
/// hashed by the `sph_cell_list` kernel, sorted on the GPU, and bucketed via
/// GPU histogram + scan. The bucket layout (`sorted_ids` + `cell_start` +
/// `cell_count`) plus the grid parameters lets a caller reconstruct each
/// particle's neighbor candidate set with the same 27-cell sweep the density
/// and force kernels use — which is exactly what the cell-list parity test
/// checks against a CPU brute-force neighbor search.
#[cfg(feature = "wgpu-backend")]
#[derive(Debug, Clone)]
pub struct GpuCellList {
    /// Particle indices sorted by cell key (length = particle count).
    pub sorted_ids: Vec<u32>,
    /// Exclusive-scan start offset into `sorted_ids` per hash bucket.
    pub cell_start: Vec<u32>,
    /// Particle count per hash bucket.
    pub cell_count: Vec<u32>,
    /// Number of hash buckets.
    pub num_cells: u32,
    /// Grid spacing (= kernel support = 2h).
    pub cell_size: f32,
    /// Domain-min corner used as the hash origin.
    pub origin: [f32; 3],
}

#[cfg(feature = "wgpu-backend")]
impl GpuCellList {
    /// Spatial hash of an integer cell coordinate, matching the WGSL kernel.
    pub fn cell_hash(&self, ix: i32, iy: i32, iz: i32) -> u32 {
        let p1: u32 = 73_856_093;
        let p2: u32 = 19_349_663;
        let p3: u32 = 83_492_791;
        let hx = (ix as u32).wrapping_mul(p1);
        let hy = (iy as u32).wrapping_mul(p2);
        let hz = (iz as u32).wrapping_mul(p3);
        (hx ^ hy ^ hz) % self.num_cells
    }

    /// Reconstruct particle `i`'s neighbor set (within `support`) by sweeping
    /// the 27 adjacent hash buckets and distance-filtering — the same traversal
    /// the GPU density/force kernels perform. `positions` is xyz-interleaved.
    pub fn neighbors_of(&self, i: usize, positions: &[f32], support: f32) -> Vec<u32> {
        let inv_cs = 1.0 / self.cell_size;
        let pix = positions[i * 3];
        let piy = positions[i * 3 + 1];
        let piz = positions[i * 3 + 2];
        let cx = ((pix - self.origin[0]) * inv_cs).floor() as i32;
        let cy = ((piy - self.origin[1]) * inv_cs).floor() as i32;
        let cz = ((piz - self.origin[2]) * inv_cs).floor() as i32;
        let support2 = support * support;

        let mut visited: Vec<u32> = Vec::with_capacity(27);
        let mut out: Vec<u32> = Vec::new();
        for dz in -1..=1 {
            for dy in -1..=1 {
                for dx in -1..=1 {
                    let h = self.cell_hash(cx + dx, cy + dy, cz + dz);
                    if visited.contains(&h) {
                        continue;
                    }
                    visited.push(h);
                    let start = self.cell_start[h as usize] as usize;
                    let count = self.cell_count[h as usize] as usize;
                    for s in 0..count {
                        let j = self.sorted_ids[start + s] as usize;
                        let rx = pix - positions[j * 3];
                        let ry = piy - positions[j * 3 + 1];
                        let rz = piz - positions[j * 3 + 2];
                        if rx * rx + ry * ry + rz * rz < support2 {
                            out.push(j as u32);
                        }
                    }
                }
            }
        }
        out.sort_unstable();
        out
    }
}

/// Build the GPU spatial-hash cell-list for `positions` (xyz-interleaved f32).
///
/// Runs the `sph_cell_list` kernel, the GPU radix sort, and the GPU
/// histogram+scan to produce a [`GpuCellList`]. `smoothing_h` sets the kernel
/// support (= 2h) used as the grid spacing; `domain_min` is the hash origin.
///
/// Returns `None` if no GPU adapter is available or a dispatch fails — callers
/// (tests) treat that as a skip.
#[cfg(feature = "wgpu-backend")]
pub fn gpu_cell_list(
    positions: &[f32],
    smoothing_h: f64,
    domain_min: [f64; 3],
) -> Option<GpuCellList> {
    gpu_pipeline::build_cell_list_standalone(positions, smoothing_h, domain_min)
}

/// SPH simulation that dispatches compute to GPU when available.
pub struct SphSimulation {
    /// Configuration (immutable after construction).
    pub config: SphConfig,
    /// Particle state.
    pub state: SphParticleState,
    /// GPU pipeline state (`None` → CPU fallback).
    #[cfg(feature = "wgpu-backend")]
    gpu: Option<SphGpuState>,
    /// Total elapsed simulation time.
    pub time: f64,
}

impl SphSimulation {
    /// Create a new SPH simulation.
    pub fn new(mut config: SphConfig) -> Self {
        if config.particle_mass == 0.0 {
            let vol = (2.0 * config.smoothing_h).powi(3);
            config.particle_mass = config.rest_density * vol;
        }
        let n = config.n_particles;
        let state = SphParticleState::new(n);

        #[cfg(feature = "wgpu-backend")]
        {
            let gpu = SphGpuState::try_init(n);
            Self {
                config,
                state,
                gpu,
                time: 0.0,
            }
        }
        #[cfg(not(feature = "wgpu-backend"))]
        {
            Self {
                config,
                state,
                time: 0.0,
            }
        }
    }

    /// True if GPU backend is active.
    pub fn has_gpu(&self) -> bool {
        #[cfg(feature = "wgpu-backend")]
        {
            self.gpu.is_some()
        }
        #[cfg(not(feature = "wgpu-backend"))]
        {
            false
        }
    }

    /// Advance the simulation by `dt` seconds.
    ///
    /// Steps:
    /// 1. Density summation (GPU cell-list or CPU brute force)
    /// 2. Pressure update (Tait EOS)
    /// 3. Pressure + viscosity acceleration
    /// 4. Velocity + position integration (symplectic Euler)
    /// 5. Boundary reflection
    pub fn step(&mut self, dt: f64) {
        let n = self.config.n_particles;

        #[cfg(feature = "wgpu-backend")]
        {
            if self.gpu.is_some() && self.try_step_gpu(dt) {
                self.time += dt;
                return;
            }
        }

        self.step_cpu(dt, n);
        self.time += dt;
    }

    // ── GPU step ──────────────────────────────────────────────────────────────

    /// Attempt the full GPU-resident pipeline for one step.
    ///
    /// Returns `false` (so the caller falls back to CPU for this step) when the
    /// GPU path is not applicable: the smoothing length is NaN or not positive,
    /// the particle arrays are shorter than the resident buffers, or a dispatch
    /// or readback failed. The GPU buffers stay resident across calls; particle
    /// state is uploaded at entry and read back at exit.
    #[cfg(feature = "wgpu-backend")]
    fn try_step_gpu(&mut self, dt: f64) -> bool {
        let cfg = self.config.clone();
        let state = &mut self.state;
        match self.gpu.as_mut() {
            Some(gpu) => gpu.step(&cfg, state, dt),
            None => false,
        }
    }

    // ── CPU step ──────────────────────────────────────────────────────────────

    fn step_cpu(&mut self, dt: f64, n: usize) {
        // 1. Density summation
        let h = self.config.smoothing_h;
        let m = self.config.particle_mass;
        let support2 = (2.0 * h) * (2.0 * h);

        for i in 0..n {
            let mut rho = 0.0;
            for j in 0..n {
                let dx = self.state.pos_x[i] - self.state.pos_x[j];
                let dy = self.state.pos_y[i] - self.state.pos_y[j];
                let dz = self.state.pos_z[i] - self.state.pos_z[j];
                let r2 = dx * dx + dy * dy + dz * dz;
                if r2 < support2 {
                    rho += m * cubic_spline_w3(r2.sqrt(), h);
                }
            }
            self.state.density[i] = rho.max(1e-6);
        }

        self.pressure_and_integrate(dt, n);
    }

    fn pressure_and_integrate(&mut self, dt: f64, n: usize) {
        let rho0 = self.config.rest_density;
        let k = self.config.pressure_k;
        let nu = self.config.viscosity;
        let m = self.config.particle_mass;
        let h = self.config.smoothing_h;
        let g = self.config.gravity;
        let support2 = (2.0 * h) * (2.0 * h);

        // 2. Tait EOS: p = k (ρ/ρ₀ − 1)
        for i in 0..n {
            self.state.pressure[i] = k * (self.state.density[i] / rho0 - 1.0);
        }

        // 3. Accelerations (collect then apply to avoid borrow conflict)
        let mut ax = vec![0.0_f64; n];
        let mut ay = vec![-g; n]; // gravity
        let mut az = vec![0.0_f64; n];

        for i in 0..n {
            let pi = self.state.pressure[i];
            let rhi = self.state.density[i];

            for j in 0..n {
                if i == j {
                    continue;
                }
                let dx = self.state.pos_x[i] - self.state.pos_x[j];
                let dy = self.state.pos_y[i] - self.state.pos_y[j];
                let dz = self.state.pos_z[i] - self.state.pos_z[j];
                let r2 = dx * dx + dy * dy + dz * dz;
                if r2 < support2 && r2 > 1e-12 {
                    let r = r2.sqrt();
                    let pj = self.state.pressure[j];
                    let rhj = self.state.density[j];

                    // Pressure term (symmetric)
                    let dw = cubic_spline_dw_dr(r, h);
                    let pf = -m * (pi / (rhi * rhi) + pj / (rhj * rhj)) * dw;
                    ax[i] += pf * dx / r;
                    ay[i] += pf * dy / r;
                    az[i] += pf * dz / r;

                    // Viscosity (Monaghan)
                    let vdotr = (self.state.vel_x[i] - self.state.vel_x[j]) * dx
                        + (self.state.vel_y[i] - self.state.vel_y[j]) * dy
                        + (self.state.vel_z[i] - self.state.vel_z[j]) * dz;
                    if vdotr < 0.0 {
                        let vf = nu * m / rhj * vdotr / (r2 + 0.01 * h * h) * dw / r;
                        ax[i] += vf * dx;
                        ay[i] += vf * dy;
                        az[i] += vf * dz;
                    }
                }
            }
        }

        // 4. Symplectic Euler integration
        for i in 0..n {
            self.state.vel_x[i] += ax[i] * dt;
            self.state.vel_y[i] += ay[i] * dt;
            self.state.vel_z[i] += az[i] * dt;
            self.state.pos_x[i] += self.state.vel_x[i] * dt;
            self.state.pos_y[i] += self.state.vel_y[i] * dt;
            self.state.pos_z[i] += self.state.vel_z[i] * dt;
        }

        // 5. Domain reflection (AABB walls)
        let [xmin, ymin, zmin] = self.config.domain_min;
        let [xmax, ymax, zmax] = self.config.domain_max;
        let e = self.config.boundary_restitution;
        macro_rules! reflect {
            ($pos:expr, $vel:expr, $min:expr, $max:expr) => {
                if $pos < $min {
                    $pos = $min;
                    $vel = $vel.abs() * e;
                }
                if $pos > $max {
                    $pos = $max;
                    $vel = -$vel.abs() * e;
                }
            };
        }
        for i in 0..n {
            reflect!(self.state.pos_x[i], self.state.vel_x[i], xmin, xmax);
            reflect!(self.state.pos_y[i], self.state.vel_y[i], ymin, ymax);
            reflect!(self.state.pos_z[i], self.state.vel_z[i], zmin, zmax);
        }
    }

    /// Compute total kinetic energy (J) across all particles.
    pub fn kinetic_energy(&self) -> f64 {
        let m = self.config.particle_mass;
        let n = self.config.n_particles;
        (0..n)
            .map(|i| {
                let v2 = self.state.vel_x[i].powi(2)
                    + self.state.vel_y[i].powi(2)
                    + self.state.vel_z[i].powi(2);
                0.5 * m * v2
            })
            .sum()
    }

    /// Mean density across all particles.
    pub fn mean_density(&self) -> f64 {
        self.state.density.iter().sum::<f64>() / self.config.n_particles as f64
    }
}

// ── GPU pipeline implementation ─────────────────────────────────────────────

#[cfg(feature = "wgpu-backend")]
mod resident;

#[cfg(feature = "wgpu-backend")]
mod gpu_pipeline {
    use super::resident::ResidentCellList;
    use super::{SphConfig, SphParticleState};
    use crate::compute::wgpu_backend::WgpuBufferHandle;
    use crate::compute::wgpu_backend::real::WgpuBackendReal;
    use crate::kernels_wgsl::{
        SPH_BOUNDARY_WGSL, SPH_DENSITY_GRID_WGSL, SPH_FORCE_GRID_WGSL, SPH_INTEGRATE_WGSL,
    };

    /// Maximum spatial-hash table size (keeps the bucket scan at two levels and
    /// the memory bounded; larger particle counts reuse buckets via hashing).
    const MAX_CELLS: usize = 1 << 16;

    /// Threads per workgroup of the SPH kernels.
    const SPH_WORKGROUP: usize = 64;

    /// Largest particle count the 1-D SPH dispatches can address (65535
    /// workgroups per dimension); bigger problems stay on the CPU path.
    const MAX_GPU_PARTICLES: usize = 65_535 * SPH_WORKGROUP;

    /// Byte size of the `SphGridParams` buffer (12 words, see [`density_params_bytes`]).
    const DENSITY_PARAM_BYTES: u64 = 12 * 4;
    /// Byte size of the `SphForceParams` buffer (16 words, see [`force_params_bytes`]).
    const FORCE_PARAM_BYTES: u64 = 16 * 4;
    /// Byte size of the `IntegrateParams` buffer (4 words, see [`integrate_params_bytes`]).
    const INTEGRATE_PARAM_BYTES: u64 = 4 * 4;
    /// Byte size of the `BoundaryParams` buffer (8 words, see [`boundary_params_bytes`]).
    const BOUNDARY_PARAM_BYTES: u64 = 8 * 4;

    fn ro() -> wgpu::BufferBindingType {
        wgpu::BufferBindingType::Storage { read_only: true }
    }

    fn rw() -> wgpu::BufferBindingType {
        wgpu::BufferBindingType::Storage { read_only: false }
    }

    /// GPU-resident SPH pipeline: one long-lived backend and every buffer of the
    /// step, allocated once in [`SphGpuState::try_init`].
    ///
    /// A step only rewrites buffer contents (particle state, per-kernel
    /// parameters) and dispatches kernels; it neither creates a device, nor
    /// compiles a pipeline after the first step, nor allocates a GPU buffer.
    pub(super) struct SphGpuState {
        backend: WgpuBackendReal,
        n: usize,
        /// Interleaved xyz positions (3·n f32).
        positions: WgpuBufferHandle,
        /// Interleaved xyz velocities (3·n f32).
        velocities: WgpuBufferHandle,
        /// Per-particle density (n f32).
        densities: WgpuBufferHandle,
        /// Interleaved xyz acceleration (3·n f32).
        accel: WgpuBufferHandle,
        /// Spatial-hash cell list: sorted ids, bucket starts and bucket counts.
        cells: ResidentCellList,
        /// `SphGridParams` words of the density kernel.
        density_params: WgpuBufferHandle,
        /// `SphForceParams` words of the force kernel.
        force_params: WgpuBufferHandle,
        /// `IntegrateParams` words of the integration kernel.
        integrate_params: WgpuBufferHandle,
        /// `BoundaryParams` words of the boundary kernel.
        boundary_params: WgpuBufferHandle,
        /// Host staging for the SoA → interleaved position upload.
        host_pos: Vec<f32>,
        /// Host staging for the SoA → interleaved velocity upload.
        host_vel: Vec<f32>,
    }

    impl SphGpuState {
        /// Try to construct the GPU state, allocating all resident buffers.
        ///
        /// Returns `None` if no GPU adapter is available or the particle count is
        /// beyond what the dispatches can address.
        pub(super) fn try_init(n: usize) -> Option<Self> {
            let n_eff = n.max(1);
            if n_eff > MAX_GPU_PARTICLES {
                return None;
            }
            let mut backend = WgpuBackendReal::try_new().ok()?;
            let num_cells = n_eff.next_power_of_two().clamp(1, MAX_CELLS);
            let cells = ResidentCellList::new(&mut backend, n_eff, num_cells).ok()?;

            let positions = backend.create_buffer_storage((n_eff * 3 * 4) as u64);
            let velocities = backend.create_buffer_storage((n_eff * 3 * 4) as u64);
            let densities = backend.create_buffer_storage((n_eff * 4) as u64);
            let accel = backend.create_buffer_storage((n_eff * 3 * 4) as u64);
            let density_params = backend.create_buffer_storage(DENSITY_PARAM_BYTES);
            let force_params = backend.create_buffer_storage(FORCE_PARAM_BYTES);
            let integrate_params = backend.create_buffer_storage(INTEGRATE_PARAM_BYTES);
            let boundary_params = backend.create_buffer_storage(BOUNDARY_PARAM_BYTES);

            Some(Self {
                backend,
                n,
                positions,
                velocities,
                densities,
                accel,
                cells,
                density_params,
                force_params,
                integrate_params,
                boundary_params,
                host_pos: vec![0.0; n * 3],
                host_vel: vec![0.0; n * 3],
            })
        }

        /// Number of GPU buffers currently alive on this state's backend.
        #[cfg(test)]
        pub(super) fn live_buffer_count(&self) -> usize {
            self.backend.live_buffer_count()
        }

        /// Total size in bytes of the live GPU buffers of this state's backend.
        #[cfg(test)]
        pub(super) fn live_buffer_bytes(&self) -> u64 {
            self.backend.live_buffer_bytes()
        }

        /// Number of compute pipelines this state's backend has compiled so far.
        #[cfg(test)]
        pub(super) fn compiled_pipeline_count(&self) -> usize {
            self.backend.compiled_pipeline_count()
        }

        /// Run one full GPU step. Returns `false` if the GPU path is not
        /// applicable for this configuration (caller falls back to CPU).
        pub(super) fn step(
            &mut self,
            cfg: &SphConfig,
            state: &mut SphParticleState,
            dt: f64,
        ) -> bool {
            let n = self.n;
            if n == 0 {
                return true;
            }
            if state.pos_x.len() < n
                || state.pos_y.len() < n
                || state.pos_z.len() < n
                || state.vel_x.len() < n
                || state.vel_y.len() < n
                || state.vel_z.len() < n
                || state.density.len() < n
                || state.pressure.len() < n
            {
                return false;
            }

            // The 27-cell neighbor sweep requires the support radius to fit in
            // one cell; otherwise the GPU result would miss neighbors. The cell
            // size is set to the support (= 2h) so neighbors always lie in the
            // 27 cells. A non-positive (or NaN) support means the GPU path is
            // inapplicable.
            let support = 2.0 * cfg.smoothing_h as f32;
            if support.is_nan() || support <= 0.0 {
                return false;
            }

            // Upload current particle state (SoA → interleaved f32) through the
            // reused host staging vectors.
            for i in 0..n {
                self.host_pos[i * 3] = state.pos_x[i] as f32;
                self.host_pos[i * 3 + 1] = state.pos_y[i] as f32;
                self.host_pos[i * 3 + 2] = state.pos_z[i] as f32;
                self.host_vel[i * 3] = state.vel_x[i] as f32;
                self.host_vel[i * 3 + 1] = state.vel_y[i] as f32;
                self.host_vel[i * 3 + 2] = state.vel_z[i] as f32;
            }
            self.backend
                .queue_write_buffer_f32(&self.positions, &self.host_pos);
            self.backend
                .queue_write_buffer_f32(&self.velocities, &self.host_vel);

            // 1–3: rebuild the spatial-hash cell list on the GPU.
            if !self.rebuild_cell_list(cfg) {
                return false;
            }

            // 4: density.
            if !self.dispatch_density(cfg) {
                return false;
            }
            // 5: force (pressure + viscosity + gravity).
            if !self.dispatch_force(cfg) {
                return false;
            }
            // 6: integrate (symplectic Euler).
            if !self.dispatch_integrate(dt) {
                return false;
            }
            // 7: boundary enforcement.
            if !self.dispatch_boundary(cfg) {
                return false;
            }

            // Read back state (interleaved f32 → SoA).
            let (Ok(pos_out), Ok(vel_out), Ok(dens_out)) = (
                self.backend.read_buffer_f32(self.positions),
                self.backend.read_buffer_f32(self.velocities),
                self.backend.read_buffer_f32(self.densities),
            ) else {
                return false;
            };
            if pos_out.len() < n * 3 || vel_out.len() < n * 3 || dens_out.len() < n {
                return false;
            }
            let rho0 = cfg.rest_density;
            let k = cfg.pressure_k;
            for i in 0..n {
                state.pos_x[i] = pos_out[i * 3] as f64;
                state.pos_y[i] = pos_out[i * 3 + 1] as f64;
                state.pos_z[i] = pos_out[i * 3 + 2] as f64;
                state.vel_x[i] = vel_out[i * 3] as f64;
                state.vel_y[i] = vel_out[i * 3 + 1] as f64;
                state.vel_z[i] = vel_out[i * 3 + 2] as f64;
                state.density[i] = dens_out[i] as f64;
                state.pressure[i] = k * (state.density[i] / rho0 - 1.0);
            }
            true
        }

        /// Common per-step scalars (particle count, bucket count, h, support).
        fn step_scalars(&self, cfg: &SphConfig) -> StepScalars {
            let h = cfg.smoothing_h as f32;
            StepScalars {
                n: self.n as u32,
                num_cells: self.cells.num_cells() as u32,
                h,
                support: 2.0 * h,
            }
        }

        /// Workgroup counts of the 1-D per-particle SPH dispatches, or `None` if
        /// the particle count cannot be covered by any dispatch shape.
        fn particle_workgroups(&self) -> Option<[u32; 3]> {
            WgpuBackendReal::dispatch_count_for(self.n, SPH_WORKGROUP as u32).ok()
        }

        /// Hash kernel → radix sort → bucket histogram → bucket scan, leaving
        /// `sorted_ids`, `cell_start` and `cell_count` on the GPU. The positions
        /// buffer is already uploaded by [`Self::step`]. Returns `false` on
        /// failure.
        fn rebuild_cell_list(&self, cfg: &SphConfig) -> bool {
            let support = 2.0 * cfg.smoothing_h as f32;
            let origin = [
                cfg.domain_min[0] as f32,
                cfg.domain_min[1] as f32,
                cfg.domain_min[2] as f32,
            ];
            self.cells
                .build(&self.backend, self.positions, support, origin)
                .is_ok()
        }

        fn dispatch_density(&self, cfg: &SphConfig) -> bool {
            let Some(workgroups) = self.particle_workgroups() else {
                return false;
            };
            let scalars = self.step_scalars(cfg);
            self.backend
                .queue_write_buffer_raw(&self.density_params, &density_params_bytes(cfg, &scalars));
            self.backend
                .dispatch_wgsl(
                    SPH_DENSITY_GRID_WGSL,
                    "sph_density",
                    &[
                        (self.positions, ro()),
                        (self.cells.sorted_ids(), ro()),
                        (self.cells.cell_start(), ro()),
                        (self.cells.cell_count(), ro()),
                        (self.densities, rw()),
                        (self.density_params, ro()),
                    ],
                    workgroups,
                )
                .is_ok()
        }

        fn dispatch_force(&self, cfg: &SphConfig) -> bool {
            let Some(workgroups) = self.particle_workgroups() else {
                return false;
            };
            let scalars = self.step_scalars(cfg);
            self.backend
                .queue_write_buffer_raw(&self.force_params, &force_params_bytes(cfg, &scalars));
            self.backend
                .dispatch_wgsl(
                    SPH_FORCE_GRID_WGSL,
                    "sph_force",
                    &[
                        (self.positions, ro()),
                        (self.velocities, ro()),
                        (self.densities, ro()),
                        (self.cells.sorted_ids(), ro()),
                        (self.cells.cell_start(), ro()),
                        (self.cells.cell_count(), ro()),
                        (self.accel, rw()),
                        (self.force_params, ro()),
                    ],
                    workgroups,
                )
                .is_ok()
        }

        fn dispatch_integrate(&self, dt: f64) -> bool {
            let Some(workgroups) = self.particle_workgroups() else {
                return false;
            };
            self.backend.queue_write_buffer_raw(
                &self.integrate_params,
                &integrate_params_bytes(self.n as u32, dt as f32),
            );
            self.backend
                .dispatch_wgsl(
                    SPH_INTEGRATE_WGSL,
                    "sph_integrate",
                    &[
                        (self.positions, rw()),
                        (self.velocities, rw()),
                        (self.accel, ro()),
                        (self.integrate_params, ro()),
                    ],
                    workgroups,
                )
                .is_ok()
        }

        fn dispatch_boundary(&self, cfg: &SphConfig) -> bool {
            let Some(workgroups) = self.particle_workgroups() else {
                return false;
            };
            self.backend.queue_write_buffer_raw(
                &self.boundary_params,
                &boundary_params_bytes(cfg, self.n as u32),
            );
            self.backend
                .dispatch_wgsl(
                    SPH_BOUNDARY_WGSL,
                    "sph_boundary",
                    &[
                        (self.positions, rw()),
                        (self.velocities, rw()),
                        (self.boundary_params, ro()),
                    ],
                    workgroups,
                )
                .is_ok()
        }
    }

    /// Standalone GPU cell-list build for parity testing (see [`super::gpu_cell_list`]).
    ///
    /// Allocates a throwaway backend and runs the same resident cell-list code the
    /// stepping pipeline uses (hash kernel, GPU radix sort, GPU bucket histogram
    /// and scan), then reads the layout back. Returns `None` if no adapter is
    /// present or any dispatch fails.
    pub(super) fn build_cell_list_standalone(
        positions: &[f32],
        smoothing_h: f64,
        domain_min: [f64; 3],
    ) -> Option<super::GpuCellList> {
        let n = positions.len() / 3;
        if n == 0 || n > MAX_GPU_PARTICLES {
            return None;
        }
        let support = 2.0 * smoothing_h as f32;
        if support.is_nan() || support <= 0.0 {
            return None;
        }
        let num_cells = n.next_power_of_two().clamp(1, MAX_CELLS);

        let mut backend = WgpuBackendReal::try_new().ok()?;
        let cells = ResidentCellList::new(&mut backend, n, num_cells).ok()?;
        let pos_buf = backend.create_buffer_storage((n * 3 * 4) as u64);
        backend.queue_write_buffer_f32(&pos_buf, &positions[..n * 3]);
        let origin = [
            domain_min[0] as f32,
            domain_min[1] as f32,
            domain_min[2] as f32,
        ];
        cells.build(&backend, pos_buf, support, origin).ok()?;
        let host = cells.read_back(&backend)?;

        Some(super::GpuCellList {
            sorted_ids: host.sorted_ids,
            cell_start: host.cell_start,
            cell_count: host.cell_count,
            num_cells: num_cells as u32,
            cell_size: support,
            origin,
        })
    }

    // ── Params packing (raw words written into read-only storage buffers) ──────

    /// Little-endian word accumulator for packing kernel `params` structs.
    ///
    /// All SPH kernels bind their `params` as a plain `array`-style storage
    /// buffer (matching the all-storage dispatch layout), so each struct is just
    /// a tight little-endian sequence of `u32`/`f32` words. Trailing padding
    /// words keep the 16-byte alignment WGSL expects for struct-in-storage.
    struct ParamWords(Vec<u8>);

    impl ParamWords {
        fn new() -> Self {
            Self(Vec::with_capacity(64))
        }
        fn u32(mut self, x: u32) -> Self {
            self.0.extend_from_slice(&x.to_le_bytes());
            self
        }
        fn f32(mut self, x: f32) -> Self {
            self.0.extend_from_slice(&x.to_le_bytes());
            self
        }
        fn vec3(self, v: [f64; 3]) -> Self {
            self.f32(v[0] as f32).f32(v[1] as f32).f32(v[2] as f32)
        }
        fn pad(mut self, words: usize) -> Self {
            for _ in 0..words {
                self.0.extend_from_slice(&0u32.to_le_bytes());
            }
            self
        }
        fn build(self) -> Vec<u8> {
            self.0
        }
    }

    /// Common derived scalars passed to the per-kernel packers.
    struct StepScalars {
        n: u32,
        num_cells: u32,
        h: f32,
        support: f32,
    }

    /// `SphGridParams`: n, num_cells, h, support, mass, cell_size, origin.xyz,
    /// + 3 pad = 12 words.
    fn density_params_bytes(cfg: &SphConfig, s: &StepScalars) -> Vec<u8> {
        ParamWords::new()
            .u32(s.n)
            .u32(s.num_cells)
            .f32(s.h)
            .f32(s.support)
            .f32(cfg.particle_mass as f32)
            .f32(s.support) // cell_size = support
            .vec3(cfg.domain_min)
            .pad(3)
            .build()
    }

    /// `SphForceParams`: n, num_cells, h, support, mass, cell_size, rho0, k, nu,
    /// gravity, origin.xyz, + 3 pad = 16 words.
    fn force_params_bytes(cfg: &SphConfig, s: &StepScalars) -> Vec<u8> {
        ParamWords::new()
            .u32(s.n)
            .u32(s.num_cells)
            .f32(s.h)
            .f32(s.support)
            .f32(cfg.particle_mass as f32)
            .f32(s.support) // cell_size = support
            .f32(cfg.rest_density as f32)
            .f32(cfg.pressure_k as f32)
            .f32(cfg.viscosity as f32)
            .f32(cfg.gravity as f32)
            .vec3(cfg.domain_min)
            .pad(3)
            .build()
    }

    /// `IntegrateParams`: n, dt, + 2 pad = 4 words.
    fn integrate_params_bytes(n: u32, dt: f32) -> Vec<u8> {
        ParamWords::new().u32(n).f32(dt).pad(2).build()
    }

    /// `BoundaryParams`: n, restitution, min.xyz, max.xyz = 8 words.
    fn boundary_params_bytes(cfg: &SphConfig, n: u32) -> Vec<u8> {
        ParamWords::new()
            .u32(n)
            .f32(cfg.boundary_restitution as f32)
            .vec3(cfg.domain_min)
            .vec3(cfg.domain_max)
            .build()
    }

    // ── tests ──────────────────────────────────────────────────────────────────

    #[cfg(test)]
    mod tests {
        use super::*;
        use crate::sph_gpu::{GpuCellList, SphSimulation};

        /// Deterministic LCG in `[0, 1)`.
        fn lcg_unit(state: &mut u32) -> f32 {
            *state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
            (*state >> 8) as f32 / (1u32 << 24) as f32
        }

        fn adapter_available() -> bool {
            if WgpuBackendReal::try_new().is_err() {
                eprintln!("SKIPPED: no GPU adapter available");
                return false;
            }
            true
        }

        /// A block of particles on a jittered lattice inside the default domain.
        fn lattice_particles(sim: &mut SphSimulation, per_axis: [usize; 3], spacing: f64) {
            let mut rng = 0x00BA_DCAB_u32;
            let [nx, ny, nz] = per_axis;
            for ix in 0..nx {
                for iy in 0..ny {
                    for iz in 0..nz {
                        let idx = (ix * ny + iy) * nz + iz;
                        let jitter = |rng: &mut u32| f64::from(lcg_unit(rng) - 0.5) * 0.1 * spacing;
                        sim.state.pos_x[idx] = -0.5 + ix as f64 * spacing + jitter(&mut rng);
                        sim.state.pos_y[idx] = 0.3 + iy as f64 * spacing + jitter(&mut rng);
                        sim.state.pos_z[idx] = -0.5 + iz as f64 * spacing + jitter(&mut rng);
                    }
                }
            }
        }

        #[test]
        fn param_packers_match_the_declared_buffer_sizes() {
            let cfg = SphConfig::default();
            let scalars = StepScalars {
                n: 8,
                num_cells: 8,
                h: 0.05,
                support: 0.1,
            };
            assert_eq!(
                density_params_bytes(&cfg, &scalars).len() as u64,
                DENSITY_PARAM_BYTES
            );
            assert_eq!(
                force_params_bytes(&cfg, &scalars).len() as u64,
                FORCE_PARAM_BYTES
            );
            assert_eq!(
                integrate_params_bytes(8, 0.01).len() as u64,
                INTEGRATE_PARAM_BYTES
            );
            assert_eq!(
                boundary_params_bytes(&cfg, 8).len() as u64,
                BOUNDARY_PARAM_BYTES
            );
        }

        /// F2095 regression: the old pipeline created 7 buffers per step and kept
        /// them forever. The buffer table must now be constant over 1000 steps.
        #[test]
        fn gpu_steps_do_not_grow_the_buffer_table() {
            if !adapter_available() {
                return;
            }
            let n = 64usize;
            let cfg = SphConfig {
                n_particles: n,
                smoothing_h: 0.1,
                domain_min: [-1.0, 0.0, -1.0],
                domain_max: [1.0, 2.0, 1.0],
                ..SphConfig::default()
            };
            let mut sim = SphSimulation::new(cfg.clone());
            lattice_particles(&mut sim, [4, 4, 4], 0.12);
            let Some(gpu) = sim.gpu.as_mut() else {
                return;
            };

            assert!(gpu.step(&cfg, &mut sim.state, 1.0 / 240.0));
            let count_after_first = gpu.live_buffer_count();
            let bytes_after_first = gpu.live_buffer_bytes();
            assert!(count_after_first > 0);
            for step in 1..1000 {
                assert!(
                    gpu.step(&cfg, &mut sim.state, 1.0 / 240.0),
                    "GPU step {step} reported failure"
                );
            }
            assert_eq!(
                gpu.live_buffer_count(),
                count_after_first,
                "buffer count grew over 1000 steps"
            );
            assert_eq!(gpu.live_buffer_bytes(), bytes_after_first);
            assert!(sim.state.pos_y.iter().all(|y| y.is_finite()));
        }

        /// F2096 regression: the old step created three wgpu devices (the sort, the
        /// bucket histogram and the scan each called `WgpuBackendReal::try_new`), so
        /// every fresh backend recompiled its pipelines on every step.  The resident
        /// pipeline owns one device: the first step fills its pipeline cache and the
        /// following steps must not add to it.
        #[test]
        fn gpu_steps_compile_every_pipeline_exactly_once() {
            if !adapter_available() {
                return;
            }
            let n = 600usize; // three radix tiles -> a two-level scan hierarchy
            let cfg = SphConfig {
                n_particles: n,
                smoothing_h: 0.1,
                ..SphConfig::default()
            };
            let mut sim = SphSimulation::new(cfg.clone());
            lattice_particles(&mut sim, [10, 10, 6], 0.12);
            let Some(gpu) = sim.gpu.as_mut() else {
                return;
            };

            assert_eq!(
                gpu.compiled_pipeline_count(),
                0,
                "nothing is compiled before the first step"
            );
            assert!(gpu.step(&cfg, &mut sim.state, 1.0 / 480.0));
            // sph_cell_list, radix_histogram, scan_block, add_block_offsets,
            // radix_scatter_pairs, sph_cell_histogram, sph_density, sph_force,
            // sph_integrate and sph_boundary.
            let compiled = gpu.compiled_pipeline_count();
            assert_eq!(compiled, 10, "the first step compiles the ten step kernels");
            for step in 1..=20 {
                assert!(
                    gpu.step(&cfg, &mut sim.state, 1.0 / 480.0),
                    "GPU step {step} reported failure"
                );
                assert_eq!(
                    gpu.compiled_pipeline_count(),
                    compiled,
                    "step {step} compiled another pipeline"
                );
            }
        }

        /// One step on more than 256 buckets (the case the old GPU histogram could
        /// not handle) must agree with the CPU reference.
        #[test]
        fn gpu_step_matches_cpu_reference_beyond_256_buckets() {
            if !adapter_available() {
                return;
            }
            let (nx, ny, nz) = (12usize, 10usize, 10usize);
            let n = nx * ny * nz; // 1200 particles -> 2048 buckets
            let h = 0.05_f64;
            let spacing = 0.6 * 2.0 * h;
            let make_config = || SphConfig {
                n_particles: n,
                smoothing_h: h,
                rest_density: 1000.0,
                pressure_k: 50.0,
                viscosity: 0.02,
                gravity: 9.81,
                particle_mass: 0.0,
                domain_min: [-1.0, 0.0, -1.0],
                domain_max: [1.0, 2.0, 1.0],
                boundary_restitution: 0.2,
            };
            let mut gpu_sim = SphSimulation::new(make_config());
            lattice_particles(&mut gpu_sim, [nx, ny, nz], spacing);
            if !gpu_sim.has_gpu() {
                return;
            }
            let mut cpu_sim = SphSimulation::new(make_config());
            lattice_particles(&mut cpu_sim, [nx, ny, nz], spacing);
            cpu_sim.gpu = None;

            let dt = 1.0 / 480.0;
            gpu_sim.step(dt);
            cpu_sim.step(dt);

            let mut sq = 0.0;
            for (g, c) in gpu_sim.state.density.iter().zip(&cpu_sim.state.density) {
                sq += (g - c) * (g - c);
            }
            let l2 = (sq / n as f64).sqrt() / 1000.0;
            assert!(
                l2 < 1e-2,
                "density L2 {l2} between GPU and CPU is too large"
            );
            for i in 0..n {
                let dy = (gpu_sim.state.pos_y[i] - cpu_sim.state.pos_y[i]).abs();
                assert!(dy < 1e-3, "particle {i}: y differs by {dy}");
            }
        }

        /// Multi-tile / multi-block / > 256-bucket cell list: for every particle
        /// the GPU layout must reproduce the brute-force neighbor set. Particles
        /// sit near cell centers so f32 floor rounding cannot differ between host
        /// and shader.
        #[test]
        fn cell_list_neighbors_match_brute_force_at_scale() {
            if !adapter_available() {
                return;
            }
            let smoothing_h = 0.05_f64;
            let support = 2.0 * smoothing_h as f32;
            let cells_per_axis = 12usize;
            let mut rng = 0x00DD_BA11_u32;
            let mut positions: Vec<f32> = Vec::new();
            for cx in 0..cells_per_axis {
                for cy in 0..cells_per_axis {
                    for cz in 0..cells_per_axis {
                        let count = 1 + (lcg_unit(&mut rng) * 2.0) as usize;
                        for _ in 0..count {
                            for c in [cx, cy, cz] {
                                let jitter = (lcg_unit(&mut rng) - 0.5) * 0.6;
                                positions.push((c as f32 + 0.5 + jitter) * support);
                            }
                        }
                    }
                }
            }
            let n = positions.len() / 3;
            assert!(n > 1500, "need a few thousand particles, got {n}");

            let Some(list): Option<GpuCellList> =
                build_cell_list_standalone(&positions, smoothing_h, [0.0; 3])
            else {
                return;
            };
            assert!(list.num_cells > 256, "must exercise more than 256 buckets");
            assert_eq!(list.cell_count.iter().sum::<u32>(), n as u32);
            let mut seen = list.sorted_ids.clone();
            seen.sort_unstable();
            assert_eq!(seen, (0..n as u32).collect::<Vec<_>>());

            let support2 = support * support;
            for i in 0..n {
                let mut expected: Vec<u32> = Vec::new();
                for j in 0..n {
                    let rx = positions[i * 3] - positions[j * 3];
                    let ry = positions[i * 3 + 1] - positions[j * 3 + 1];
                    let rz = positions[i * 3 + 2] - positions[j * 3 + 2];
                    if rx * rx + ry * ry + rz * rz < support2 {
                        expected.push(j as u32);
                    }
                }
                assert_eq!(
                    list.neighbors_of(i, &positions, support),
                    expected,
                    "neighbor set of particle {i} differs from brute force"
                );
            }
        }
    }
}

// ── tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cubic_spline_w3_normalisation() {
        // W(0, h) should be positive; W(2h, h) = 0 (beyond kernel support)
        let h = 0.1;
        assert!(cubic_spline_w3(0.0, h) > 0.0);
        assert_eq!(cubic_spline_w3(2.0 * h, h), 0.0);
        assert_eq!(cubic_spline_w3(2.1 * h, h), 0.0);
    }

    #[test]
    fn test_cubic_spline_dw_dr() {
        let h = 0.1;
        // Gradient at r=0 should be 0 (symmetric kernel)
        assert_eq!(cubic_spline_dw_dr(0.0, h), 0.0);
        // Gradient at r > 2h should be 0
        assert_eq!(cubic_spline_dw_dr(3.0 * h, h), 0.0);
    }

    #[test]
    fn test_sph_construction() {
        let sim = SphSimulation::new(SphConfig {
            n_particles: 8,
            ..SphConfig::default()
        });
        assert_eq!(sim.state.n, 8);
        assert!(sim.config.particle_mass > 0.0);
    }

    #[test]
    fn test_sph_step_falls_under_gravity() {
        let mut sim = SphSimulation::new(SphConfig {
            n_particles: 4,
            smoothing_h: 0.2,
            gravity: 9.81,
            domain_min: [-5., 0., -5.],
            domain_max: [5., 10., 5.],
            ..SphConfig::default()
        });
        // Place particles high up
        for i in 0..4 {
            sim.state.pos_y[i] = 5.0;
        }

        let dt = 1.0 / 60.0;
        for _ in 0..10 {
            sim.step(dt);
        }

        // All particles should have moved down
        for i in 0..4 {
            assert!(
                sim.state.pos_y[i] < 5.0,
                "particle {} should fall, y={}",
                i,
                sim.state.pos_y[i]
            );
        }
    }

    #[test]
    fn test_sph_boundary_reflection() {
        let mut sim = SphSimulation::new(SphConfig {
            n_particles: 1,
            smoothing_h: 0.2,
            gravity: 0.0, // No gravity so we control bounce
            domain_min: [0., 0., 0.],
            domain_max: [1., 1., 1.],
            boundary_restitution: 1.0,
            ..SphConfig::default()
        });
        sim.state.pos_y[0] = 0.5;
        sim.state.vel_y[0] = -10.0; // Moving down fast

        for _ in 0..10 {
            sim.step(0.01);
        }

        // Particle should stay within domain
        assert!(sim.state.pos_y[0] >= 0.0);
        assert!(sim.state.pos_y[0] <= 1.0);
    }

    #[test]
    fn test_sph_kinetic_energy() {
        let mut sim = SphSimulation::new(SphConfig {
            n_particles: 4,
            ..SphConfig::default()
        });
        for i in 0..4 {
            sim.state.vel_y[i] = 1.0;
        }
        let ke = sim.kinetic_energy();
        assert!(ke > 0.0, "KE should be positive");
    }
}
