// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! GPU-accelerated Smoothed Particle Hydrodynamics (SPH) simulation.
//!
//! This module runs a Weakly Compressible SPH (WCSPH) density–pressure solve on
//! the GPU through the real [`WgpuBackendReal`] dispatch path, falling back to a
//! clean CPU implementation when no GPU is available.
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
//! The GPU path is a real resident pipeline. Each `step` uploads the current
//! particle state, then dispatches the full kernel sequence:
//!
//! 1. **cell-list**: hash each particle into a spatial-hash bucket
//!    (`sph_cell_list` kernel → `cell_keys[]`, `particle_ids[]`).
//! 2. **sort**: stable radix-sort the `(cell_key, particle_id)` pairs on the GPU
//!    (`gpu_radix::radix_sort_pairs_gpu`) so each bucket's particles are
//!    contiguous in `sorted_ids[]`.
//! 3. **bucket bounds**: `histogram_u32` + `exclusive_scan_u32` (GPU primitives)
//!    yield `cell_count[]` and `cell_start[]`.
//! 4. **density**: `sph_density` walks the 27 neighbor buckets per particle,
//!    distance-filtering candidates (collision-safe) and summing W.
//! 5. **force**: `sph_force` computes Tait-EOS pressure and accumulates the
//!    symmetric Monaghan pressure force, artificial viscosity, and gravity.
//! 6. **integrate**: `sph_integrate` advances velocity then position (symplectic).
//! 7. **boundary**: `sph_boundary` clamps to the domain AABB and reflects walls.
//!
//! The positions/velocities/densities/acceleration buffers stay GPU-resident
//! across all sub-steps of a `step`; only the spatial-hash rebuild round-trips
//! through the reused GPU radix-sort/scan primitives (mirroring the shipped
//! hybrid LBVH's "GPU sort → host orchestration" structure). Results are read
//! back into the structure-of-arrays [`SphParticleState`] at the end of `step`.
//!
//! When the kernel support radius exceeds the cell size the neighbor sweep would
//! be incomplete, so the GPU path validates `support ≤ cell_size` and otherwise
//! falls back to the exact CPU solve for that step.
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
    /// GPU path is not applicable — e.g. the cell size cannot cover the kernel
    /// support, or a dispatch/readback failed. The GPU buffers stay resident
    /// across calls; particle state is uploaded at entry and read back at exit.
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
mod gpu_pipeline {
    use super::{SphConfig, SphParticleState};
    use crate::compute::wgpu_backend::WgpuBufferHandle;
    use crate::compute::wgpu_backend::real::WgpuBackendReal;
    use crate::kernels_wgsl::{
        SPH_BOUNDARY_WGSL, SPH_CELL_LIST_WGSL, SPH_DENSITY_GRID_WGSL, SPH_FORCE_GRID_WGSL,
        SPH_INTEGRATE_WGSL,
    };

    /// Maximum spatial-hash table size (keeps histogram bins ≤ device limits and
    /// memory bounded; larger particle counts reuse buckets via hashing).
    const MAX_CELLS: usize = 1 << 16;

    fn ro() -> wgpu::BufferBindingType {
        wgpu::BufferBindingType::Storage { read_only: true }
    }

    fn rw() -> wgpu::BufferBindingType {
        wgpu::BufferBindingType::Storage { read_only: false }
    }

    /// GPU-resident SPH pipeline buffers and backend.
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
        /// Sorted particle index list (n u32), produced by the radix sort.
        sorted_ids: WgpuBufferHandle,
        /// Per-bucket start offset (num_cells u32).
        cell_start: WgpuBufferHandle,
        /// Per-bucket particle count (num_cells u32).
        cell_count: WgpuBufferHandle,
        /// Number of hash buckets.
        num_cells: usize,
    }

    impl SphGpuState {
        /// Try to construct the GPU state, allocating all resident buffers.
        ///
        /// Returns `None` if no GPU adapter is available.
        pub(super) fn try_init(n: usize) -> Option<Self> {
            let mut backend = WgpuBackendReal::try_new().ok()?;
            let n_eff = n.max(1);
            let num_cells = n_eff.next_power_of_two().clamp(1, MAX_CELLS);

            let positions = backend.create_buffer_storage((n_eff * 3 * 4) as u64);
            let velocities = backend.create_buffer_storage((n_eff * 3 * 4) as u64);
            let densities = backend.create_buffer_storage((n_eff * 4) as u64);
            let accel = backend.create_buffer_storage((n_eff * 3 * 4) as u64);
            let sorted_ids = backend.create_buffer_storage((n_eff * 4) as u64);
            let cell_start = backend.create_buffer_storage((num_cells * 4) as u64);
            let cell_count = backend.create_buffer_storage((num_cells * 4) as u64);

            Some(Self {
                backend,
                n,
                positions,
                velocities,
                densities,
                accel,
                sorted_ids,
                cell_start,
                cell_count,
                num_cells,
            })
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

            // The 27-cell neighbor sweep requires the support radius to fit in
            // one cell; otherwise the GPU result would miss neighbors. The cell
            // size is set to the support (= 2h) so neighbors always lie in the
            // 27 cells. A non-positive support means the GPU path is inapplicable.
            let support = 2.0 * cfg.smoothing_h as f32;
            if support <= 0.0 {
                return false;
            }

            // Upload current particle state (SoA → interleaved f32).
            let mut pos = vec![0.0f32; n * 3];
            let mut vel = vec![0.0f32; n * 3];
            for i in 0..n {
                pos[i * 3] = state.pos_x[i] as f32;
                pos[i * 3 + 1] = state.pos_y[i] as f32;
                pos[i * 3 + 2] = state.pos_z[i] as f32;
                vel[i * 3] = state.vel_x[i] as f32;
                vel[i * 3 + 1] = state.vel_y[i] as f32;
                vel[i * 3 + 2] = state.vel_z[i] as f32;
            }
            self.backend
                .queue_write_buffer_f32(&self.positions, &pos);
            self.backend
                .queue_write_buffer_f32(&self.velocities, &vel);

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
            if !self.dispatch_integrate(cfg, dt) {
                return false;
            }
            // 7: boundary enforcement.
            if !self.dispatch_boundary(cfg) {
                return false;
            }

            // Read back state (interleaved f32 → SoA).
            let pos_out = self.backend.read_buffer_f32(self.positions);
            let vel_out = self.backend.read_buffer_f32(self.velocities);
            let dens_out = self.backend.read_buffer_f32(self.densities);
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
                num_cells: self.num_cells as u32,
                h,
                support: 2.0 * h,
            }
        }

        /// CELL_LIST → radix sort → histogram → scan, writing `sorted_ids`,
        /// `cell_start`, `cell_count` on the GPU. The positions buffer is already
        /// uploaded by [`Self::step`]. Returns `false` on failure.
        fn rebuild_cell_list(&mut self, cfg: &SphConfig) -> bool {
            let n = self.n;

            // Scratch key/payload buffers.
            let cell_keys = self.backend.create_buffer_storage((n * 4) as u64);
            let particle_ids = self.backend.create_buffer_storage((n * 4) as u64);

            let scalars = self.step_scalars(cfg);
            let params = cell_params_bytes(cfg, &scalars);
            let params_buf = self.backend.create_buffer_storage(params.len() as u64);
            self.backend.queue_write_buffer_raw(&params_buf, &params);

            let wg = WgpuBackendReal::dispatch_count_for(n, 64);
            if self
                .backend
                .dispatch_wgsl(
                    SPH_CELL_LIST_WGSL,
                    "sph_cell_list",
                    &[
                        (self.positions, ro()),
                        (cell_keys, rw()),
                        (particle_ids, rw()),
                        (params_buf, ro()),
                    ],
                    wg,
                )
                .is_err()
            {
                return false;
            }

            // Read keys+ids, sort the pairs on the GPU, derive bucket bounds.
            let keys = self.backend.read_buffer_u32(cell_keys);
            let ids = self.backend.read_buffer_u32(particle_ids);
            if keys.len() < n || ids.len() < n {
                return false;
            }
            let keys = &keys[..n];
            let ids = &ids[..n];

            let (sorted_keys, sorted_ids) = crate::gpu_radix::radix_sort_pairs_gpu(keys, ids);
            let counts = crate::gpu_primitives::histogram_u32(&sorted_keys, self.num_cells);
            let starts = crate::gpu_primitives::exclusive_scan_u32(&counts);
            if sorted_ids.len() < n || counts.len() < self.num_cells || starts.len() < self.num_cells
            {
                return false;
            }

            self.backend
                .queue_write_buffer_raw(&self.sorted_ids, bytemuck::cast_slice(&sorted_ids[..n]));
            self.backend.queue_write_buffer_raw(
                &self.cell_start,
                bytemuck::cast_slice(&starts[..self.num_cells]),
            );
            self.backend.queue_write_buffer_raw(
                &self.cell_count,
                bytemuck::cast_slice(&counts[..self.num_cells]),
            );
            true
        }

        fn dispatch_density(&mut self, cfg: &SphConfig) -> bool {
            let n = self.n;
            let scalars = self.step_scalars(cfg);
            let params = density_params_bytes(cfg, &scalars);
            let params_buf = self.backend.create_buffer_storage(params.len() as u64);
            self.backend.queue_write_buffer_raw(&params_buf, &params);

            let wg = WgpuBackendReal::dispatch_count_for(n, 64);
            self.backend
                .dispatch_wgsl(
                    SPH_DENSITY_GRID_WGSL,
                    "sph_density",
                    &[
                        (self.positions, ro()),
                        (self.sorted_ids, ro()),
                        (self.cell_start, ro()),
                        (self.cell_count, ro()),
                        (self.densities, rw()),
                        (params_buf, ro()),
                    ],
                    wg,
                )
                .is_ok()
        }

        fn dispatch_force(&mut self, cfg: &SphConfig) -> bool {
            let n = self.n;
            let scalars = self.step_scalars(cfg);
            let params = force_params_bytes(cfg, &scalars);
            let params_buf = self.backend.create_buffer_storage(params.len() as u64);
            self.backend.queue_write_buffer_raw(&params_buf, &params);

            let wg = WgpuBackendReal::dispatch_count_for(n, 64);
            self.backend
                .dispatch_wgsl(
                    SPH_FORCE_GRID_WGSL,
                    "sph_force",
                    &[
                        (self.positions, ro()),
                        (self.velocities, ro()),
                        (self.densities, ro()),
                        (self.sorted_ids, ro()),
                        (self.cell_start, ro()),
                        (self.cell_count, ro()),
                        (self.accel, rw()),
                        (params_buf, ro()),
                    ],
                    wg,
                )
                .is_ok()
        }

        fn dispatch_integrate(&mut self, _cfg: &SphConfig, dt: f64) -> bool {
            let n = self.n;
            let params = integrate_params_bytes(n as u32, dt as f32);
            let params_buf = self.backend.create_buffer_storage(params.len() as u64);
            self.backend.queue_write_buffer_raw(&params_buf, &params);

            let wg = WgpuBackendReal::dispatch_count_for(n, 64);
            self.backend
                .dispatch_wgsl(
                    SPH_INTEGRATE_WGSL,
                    "sph_integrate",
                    &[
                        (self.positions, rw()),
                        (self.velocities, rw()),
                        (self.accel, ro()),
                        (params_buf, ro()),
                    ],
                    wg,
                )
                .is_ok()
        }

        fn dispatch_boundary(&mut self, cfg: &SphConfig) -> bool {
            let n = self.n;
            let params = boundary_params_bytes(cfg, n as u32);
            let params_buf = self.backend.create_buffer_storage(params.len() as u64);
            self.backend.queue_write_buffer_raw(&params_buf, &params);

            let wg = WgpuBackendReal::dispatch_count_for(n, 64);
            self.backend
                .dispatch_wgsl(
                    SPH_BOUNDARY_WGSL,
                    "sph_boundary",
                    &[
                        (self.positions, rw()),
                        (self.velocities, rw()),
                        (params_buf, ro()),
                    ],
                    wg,
                )
                .is_ok()
        }
    }

    /// Standalone GPU cell-list build for parity testing (see [`super::gpu_cell_list`]).
    ///
    /// Allocates a throwaway backend, hashes the particles with the cell-list
    /// kernel, sorts the pairs and buckets them, and returns the GPU-produced
    /// layout. Returns `None` if no adapter is present or any dispatch fails.
    pub(super) fn build_cell_list_standalone(
        positions: &[f32],
        smoothing_h: f64,
        domain_min: [f64; 3],
    ) -> Option<super::GpuCellList> {
        let n = positions.len() / 3;
        if n == 0 {
            return None;
        }
        let support = 2.0 * smoothing_h as f32;
        if support <= 0.0 {
            return None;
        }
        let num_cells = n.next_power_of_two().clamp(1, MAX_CELLS);

        let mut backend = WgpuBackendReal::try_new().ok()?;
        let pos_buf = backend.create_buffer_storage((n * 3 * 4) as u64);
        let cell_keys = backend.create_buffer_storage((n * 4) as u64);
        let particle_ids = backend.create_buffer_storage((n * 4) as u64);
        backend.queue_write_buffer_f32(&pos_buf, positions);

        let params = ParamWords::new()
            .u32(n as u32)
            .u32(num_cells as u32)
            .f32(support)
            .f32(domain_min[0] as f32)
            .f32(domain_min[1] as f32)
            .f32(domain_min[2] as f32)
            .pad(2)
            .build();
        let params_buf = backend.create_buffer_storage(params.len() as u64);
        backend.queue_write_buffer_raw(&params_buf, &params);

        let wg = WgpuBackendReal::dispatch_count_for(n, 64);
        backend
            .dispatch_wgsl(
                SPH_CELL_LIST_WGSL,
                "sph_cell_list",
                &[
                    (pos_buf, ro()),
                    (cell_keys, rw()),
                    (particle_ids, rw()),
                    (params_buf, ro()),
                ],
                wg,
            )
            .ok()?;

        let keys = backend.read_buffer_u32(cell_keys);
        let ids = backend.read_buffer_u32(particle_ids);
        if keys.len() < n || ids.len() < n {
            return None;
        }
        let (sorted_keys, sorted_ids) = crate::gpu_radix::radix_sort_pairs_gpu(&keys[..n], &ids[..n]);
        let cell_count = crate::gpu_primitives::histogram_u32(&sorted_keys, num_cells);
        let cell_start = crate::gpu_primitives::exclusive_scan_u32(&cell_count);
        if sorted_ids.len() < n || cell_count.len() < num_cells || cell_start.len() < num_cells {
            return None;
        }

        Some(super::GpuCellList {
            sorted_ids: sorted_ids[..n].to_vec(),
            cell_start: cell_start[..num_cells].to_vec(),
            cell_count: cell_count[..num_cells].to_vec(),
            num_cells: num_cells as u32,
            cell_size: support,
            origin: [
                domain_min[0] as f32,
                domain_min[1] as f32,
                domain_min[2] as f32,
            ],
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

    /// `CellParams`: n, num_cells, cell_size, origin.xyz, + 2 pad = 8 words.
    fn cell_params_bytes(cfg: &SphConfig, s: &StepScalars) -> Vec<u8> {
        ParamWords::new()
            .u32(s.n)
            .u32(s.num_cells)
            .f32(s.support) // cell_size = support
            .vec3(cfg.domain_min)
            .pad(2)
            .build()
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
