// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! GPU-accelerated Lattice Boltzmann Method (LBM) fluid simulation.
//!
//! Implements the D3Q19 single-relaxation-time (SRT) Bhatnagar–Gross–Krook
//! (BGK) LBM on a regular Cartesian grid.  GPU dispatches are made via the
//! `oxiphysics-gpu` compute backend; a reference CPU implementation is used as
//! the fallback.
//!
//! # Physical model
//!
//! The D3Q19 BGK equation:
//!
//! f_i(x + e_i Δt, t + Δt) = f_i(x,t) − (1/τ) [f_i − f_i^eq]
//!
//! where the equilibrium distribution is:
//!
//! f_i^eq = wᵢ ρ [1 + (eᵢ·u)/cs² + (eᵢ·u)²/(2cs⁴) − |u|²/(2cs²)]
//!
//! and cs² = 1/3 (in lattice units, Δx=Δt=1).
//!
//! The relaxation time τ relates to kinematic viscosity: ν = cs²(τ − 0.5)Δt.
//!
//! ## Forcing
//!
//! Body forces use the Guo, Zheng & Shi (2002) scheme: the equilibrium is
//! evaluated at the shifted velocity u = (Σ fᵢ eᵢ + F/2)/ρ and the source term
//! Sᵢ = (1 − 1/(2τ)) wᵢ [3(eᵢ − u) + 9(eᵢ·u) eᵢ]·F is added after collision.
//!
//! ## Boundaries
//!
//! * Solid cells (see [`LbmSimulation::set_solid`]) are no-slip walls handled
//!   by halfway bounce-back; the wall sits halfway between the last fluid node
//!   and the first solid node.
//! * Solid cells on the lid plane `y = ny − 1` move with the lid velocity
//!   ([`LbmSimulation::set_lid_velocity`]) using Ladd's moving-wall
//!   bounce-back fᵢ̄(x, t+1) = fᵢ*(x, t) − 6 wᵢ ρ (eᵢ·u_lid).
//! * Directions without a solid neighbour wrap periodically.
//!
//! The CPU reference ([`LbmSimulation`] with the GPU disabled) and the WGSL
//! kernel (`shaders/lbm_bgk_d3q19.wgsl`) implement exactly the same
//! collide-then-push algorithm; the GPU path stores populations in `f32`.
//!
//! ## Grid layout
//!
//! Flat SoA: one `Vec<f64>` per velocity direction (19 arrays of N×M×P cells).
//! Allows GPU kernels to process each direction slice in parallel.
//!
//! ## Usage
//!
//! ```
//! use oxiphysics_gpu::lbm_gpu::{LbmSimulation, LbmConfig};
//!
//! let cfg = LbmConfig { nx: 8, ny: 8, nz: 8, tau: 0.6, ..LbmConfig::default() };
//! let mut sim = LbmSimulation::new(cfg);
//!
//! // Drive the lid (top layer, y = ny - 1) at u = 0.1 in X
//! sim.set_lid_velocity(0.1, 0.0, 0.0);
//!
//! for _ in 0..20 { sim.step(); }
//!
//! // The lid drags the fluid along +X: the mean flow is positive in X and
//! // the layer just below the lid moves faster than the layer above the floor.
//! let (ux, _, _) = sim.mean_velocity();
//! assert!(ux > 0.0, "mean ux = {ux}");
//! let (_, u_top) = sim.cell_macro(4, 6, 4);
//! let (_, u_bottom) = sim.cell_macro(4, 1, 4);
//! assert!(u_top[0] > u_bottom[0], "top {} vs bottom {}", u_top[0], u_bottom[0]);
//! ```

use crate::GpuError;
use crate::compute_preference::ComputePreference;
#[cfg(feature = "wgpu-backend")]
use {crate::compute::WgpuBufferHandle, crate::compute::wgpu_backend::real::WgpuBackendReal, wgpu};

/// WGSL source of the D3Q19 BGK kernel dispatched by the GPU path.
///
/// Bindings: `0` params (`[nx, ny, nz, ω, Fx, Fy, Fz, u_lid.x, u_lid.y,
/// u_lid.z]` as `u32` words, floats bit-cast), `1` `f_in`, `2` `f_out`
/// (q-major `f32`), `3` solid mask (`u32` per cell).
pub const LBM_BGK_D3Q19_WGSL: &str = include_str!("shaders/lbm_bgk_d3q19.wgsl");

/// Number of `u32` words in the GPU params buffer.
#[cfg(feature = "wgpu-backend")]
const LBM_PARAM_WORDS: usize = 10;

/// Workgroup edge length of the WGSL kernel (`@workgroup_size(4, 4, 4)`).
#[cfg(feature = "wgpu-backend")]
const LBM_WORKGROUP_EDGE: u32 = 4;

// ── D3Q19 velocity set ────────────────────────────────────────────────────────

/// D3Q19 velocity directions (ex, ey, ez) × 19.
pub const D3Q19_EX: [i32; 19] = [0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 0, 0, 0, 0];
pub const D3Q19_EY: [i32; 19] = [0, 0, 0, 1, -1, 0, 0, 1, 1, -1, -1, 0, 0, 0, 0, 1, -1, 1, -1];
pub const D3Q19_EZ: [i32; 19] = [0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, 1, -1, -1, 1, 1, -1, -1];

/// D3Q19 weights wᵢ.
pub const D3Q19_W: [f64; 19] = [
    1.0 / 3.0, // rest
    1.0 / 18.0,
    1.0 / 18.0,
    1.0 / 18.0,
    1.0 / 18.0,
    1.0 / 18.0,
    1.0 / 18.0, // face
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0, // edge
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
];

/// Opposite direction index for bounce-back: opp\[i\] is the index j such that e_j = −e_i.
///
/// Directions and their opposites (all components are negated):
///  0: ( 0, 0, 0) → 0   1: (+1, 0, 0) → 2   2: (−1, 0, 0) → 1
///  3: ( 0,+1, 0) → 4   4: ( 0,−1, 0) → 3   5: ( 0, 0,+1) → 6
///  6: ( 0, 0,−1) → 5   7: (+1,+1, 0) → 10  8: (−1,+1, 0) → 9
///  9: (+1,−1, 0) → 8  10: (−1,−1, 0) → 7  11: (+1, 0,+1) → 14
/// 12: (−1, 0,+1) → 13 13: (+1, 0,−1) → 12 14: (−1, 0,−1) → 11
/// 15: ( 0,+1,+1) → 18 16: ( 0,−1,+1) → 17 17: ( 0,+1,−1) → 16
/// 18: ( 0,−1,−1) → 15
pub const D3Q19_OPP: [usize; 19] = [
    0, 2, 1, 4, 3, 6, 5, 10, 9, 8, 7, 14, 13, 12, 11, 18, 17, 16, 15,
];

// ── LbmConfig ─────────────────────────────────────────────────────────────────

/// Configuration for an LBM simulation.
#[derive(Debug, Clone)]
pub struct LbmConfig {
    /// Grid size in X direction (cells).
    pub nx: usize,
    /// Grid size in Y direction (cells).
    pub ny: usize,
    /// Grid size in Z direction (cells).
    pub nz: usize,
    /// BGK relaxation time τ (lattice units, τ > 0.5 for stability).
    pub tau: f64,
    /// Initial density (ρ₀, lattice units, default 1.0).
    pub rho0: f64,
    /// Body force in X (lattice units per step²).
    pub force_x: f64,
    /// Body force in Y.
    pub force_y: f64,
    /// Body force in Z.
    pub force_z: f64,
    /// Which compute device may run the steps (default
    /// [`ComputePreference::Auto`]).
    ///
    /// * `Auto` — GPU when an adapter is available, CPU reference otherwise.
    /// * `CpuOnly` — never initialise the GPU path (the same as calling
    ///   [`LbmSimulation::disable_gpu`] right after construction).
    /// * `GpuRequired` — [`LbmSimulation::try_step`] reports an error instead
    ///   of running a step on the CPU.
    pub compute: ComputePreference,
}

impl Default for LbmConfig {
    fn default() -> Self {
        Self {
            nx: 16,
            ny: 16,
            nz: 16,
            tau: 0.6,
            rho0: 1.0,
            force_x: 0.0,
            force_y: 0.0,
            force_z: 0.0,
            compute: ComputePreference::Auto,
        }
    }
}

impl LbmConfig {
    /// Kinematic viscosity ν = cs²(τ − 0.5) in lattice units (cs² = 1/3).
    pub fn viscosity(&self) -> f64 {
        (1.0 / 3.0) * (self.tau - 0.5)
    }

    /// Total number of cells.
    pub fn n_cells(&self) -> usize {
        self.nx * self.ny * self.nz
    }
}

// ── LbmSimulation ─────────────────────────────────────────────────────────────

/// GPU resources owned by `LbmSimulation` when the `wgpu-backend` feature is active.
///
/// Holds ping-pong f32 buffers for f_in / f_out, the params buffer and the
/// solid mask.  Upload happens at lazy init (and again whenever the CPU-side
/// geometry, parameters or populations are changed through the setters);
/// readback happens lazily when macroscopic quantities are requested after
/// one or more GPU steps.
#[cfg(feature = "wgpu-backend")]
struct LbmGpuState {
    backend: WgpuBackendReal,
    /// Buffer A — alternates between input and output roles.
    f_buf_a: WgpuBufferHandle,
    /// Buffer B — alternates between input and output roles.
    f_buf_b: WgpuBufferHandle,
    /// Params buffer: `LBM_PARAM_WORDS` × u32 (see [`LBM_BGK_D3Q19_WGSL`]).
    params_buf: WgpuBufferHandle,
    /// Solid mask: one u32 per cell (1 = solid).
    solid_buf: WgpuBufferHandle,
    /// When `false`, A is the current input; when `true`, B is current input.
    b_is_input: bool,
    /// Copy of the params words last written to `params_buf`.
    uploaded_params: [u32; LBM_PARAM_WORDS],
    /// Workgroup counts `[ceil(nx/4), ceil(ny/4), ceil(nz/4)]`.
    workgroups: [u32; 3],
}

#[cfg(feature = "wgpu-backend")]
impl LbmGpuState {
    /// Try to create GPU state, uploading populations, solid mask and
    /// parameters from `sim`.
    ///
    /// Returns `Err` when the grid cannot be represented on the GPU or when no
    /// compatible adapter is available.
    fn try_new(sim: &LbmSimulation) -> Result<Self, GpuError> {
        let nc = sim.config.n_cells();
        if nc == 0 {
            return Err(GpuError::BackendInit("LBM grid has no cells".to_string()));
        }
        let params = sim.gpu_params()?;
        let workgroups = [
            params[0].div_ceil(LBM_WORKGROUP_EDGE),
            params[1].div_ceil(LBM_WORKGROUP_EDGE),
            params[2].div_ceil(LBM_WORKGROUP_EDGE),
        ];
        let f_bytes = 19_usize
            .checked_mul(nc)
            .and_then(|v| v.checked_mul(std::mem::size_of::<f32>()))
            .and_then(|v| u64::try_from(v).ok())
            .ok_or_else(|| {
                GpuError::BackendInit(format!("LBM population buffer for {nc} cells overflows"))
            })?;
        let solid_bytes = nc
            .checked_mul(std::mem::size_of::<u32>())
            .and_then(|v| u64::try_from(v).ok())
            .ok_or_else(|| {
                GpuError::BackendInit(format!("LBM solid mask for {nc} cells overflows"))
            })?;
        let params_bytes = (LBM_PARAM_WORDS * std::mem::size_of::<u32>()) as u64;

        let mut backend =
            WgpuBackendReal::try_new().map_err(|e| GpuError::BackendInit(e.to_string()))?;

        let f_buf_a = backend.create_buffer_storage(f_bytes);
        let f_buf_b = backend.create_buffer_storage(f_bytes);
        let params_buf = backend.create_buffer_storage(params_bytes);
        let solid_buf = backend.create_buffer_storage(solid_bytes);

        backend.queue_write_buffer_raw(&params_buf, bytemuck::cast_slice(&params));
        backend.queue_write_buffer_raw(&solid_buf, bytemuck::cast_slice(&solid_words(&sim.solid)));
        backend.queue_write_buffer_f32(&f_buf_a, &flatten_soa_to_f32(&sim.f, nc));

        Ok(Self {
            backend,
            f_buf_a,
            f_buf_b,
            params_buf,
            solid_buf,
            b_is_input: false,
            uploaded_params: params,
            workgroups,
        })
    }

    /// Current input buffer (the one last written by the shader).
    fn input_buf(&self) -> WgpuBufferHandle {
        if self.b_is_input {
            self.f_buf_b
        } else {
            self.f_buf_a
        }
    }

    /// Current output buffer (the one the shader will write next).
    fn output_buf(&self) -> WgpuBufferHandle {
        if self.b_is_input {
            self.f_buf_a
        } else {
            self.f_buf_b
        }
    }

    /// Overwrite the current input populations with the CPU-side `f`.
    fn upload_populations(&self, f: &[Vec<f64>], nc: usize) {
        let input = self.input_buf();
        self.backend
            .queue_write_buffer_f32(&input, &flatten_soa_to_f32(f, nc));
    }

    /// Overwrite the GPU solid mask.
    fn upload_solid(&self, solid: &[bool]) {
        self.backend
            .queue_write_buffer_raw(&self.solid_buf, bytemuck::cast_slice(&solid_words(solid)));
    }

    /// Re-upload the params words when they differ from the last upload.
    fn upload_params_if_changed(&mut self, params: [u32; LBM_PARAM_WORDS]) {
        if params != self.uploaded_params {
            self.backend
                .queue_write_buffer_raw(&self.params_buf, bytemuck::cast_slice(&params));
            self.uploaded_params = params;
        }
    }

    /// Run one collide-and-stream kernel dispatch and swap the ping-pong
    /// buffers.  On error the buffers are left unswapped, so the input buffer
    /// still holds the last valid state.
    fn dispatch(&mut self) -> Result<(), GpuError> {
        let input = self.input_buf();
        let output = self.output_buf();
        self.backend
            .dispatch_wgsl(
                LBM_BGK_D3Q19_WGSL,
                "main",
                &[
                    (
                        self.params_buf,
                        wgpu::BufferBindingType::Storage { read_only: true },
                    ),
                    (input, wgpu::BufferBindingType::Storage { read_only: true }),
                    (
                        output,
                        wgpu::BufferBindingType::Storage { read_only: false },
                    ),
                    (
                        self.solid_buf,
                        wgpu::BufferBindingType::Storage { read_only: true },
                    ),
                ],
                self.workgroups,
            )
            .map_err(|e| GpuError::ShaderDispatch(e.to_string()))?;
        self.b_is_input = !self.b_is_input;
        Ok(())
    }

    /// Read the current populations back into SoA `f64` form.
    fn read_populations(&self, nc: usize) -> Result<Vec<Vec<f64>>, GpuError> {
        let flat = self
            .backend
            .read_buffer_f32(self.input_buf())
            .map_err(|e| GpuError::ReadBack(e.to_string()))?;
        unflatten_f32_to_soa(&flat, nc)
    }
}

/// Availability of the GPU path of an [`LbmSimulation`].
#[cfg(feature = "wgpu-backend")]
enum GpuSlot {
    /// GPU initialisation has not been attempted yet (lazy, on first step).
    Untried,
    /// GPU resources are live; steps run on the GPU.
    Ready(Box<LbmGpuState>),
    /// Initialisation or a later dispatch/readback failed.  The failure is
    /// remembered so the (expensive) adapter request is never retried and all
    /// further steps run on the CPU.
    Unavailable(GpuError),
    /// The GPU path was switched off explicitly with
    /// [`LbmSimulation::disable_gpu`].
    Disabled,
}

/// Convert the boolean solid mask into the `u32` words used by the kernel.
#[cfg(feature = "wgpu-backend")]
fn solid_words(solid: &[bool]) -> Vec<u32> {
    solid.iter().map(|&s| u32::from(s)).collect()
}

/// Flatten SoA distribution functions to a flat f32 buffer using q-major layout.
///
/// q-major: `flat[q * nc + cell] = f[q][cell]`
///
/// This matches the WGSL indexing `q * (nx*ny*nz) + cell` with
/// `cell = x + nx*(y + ny*z)`.
#[cfg(feature = "wgpu-backend")]
fn flatten_soa_to_f32(f: &[Vec<f64>], nc: usize) -> Vec<f32> {
    let mut flat = Vec::with_capacity(19 * nc);
    for dir in f.iter().take(19) {
        for &val in dir.iter().take(nc) {
            flat.push(val as f32);
        }
    }
    flat
}

/// Unflatten a q-major f32 buffer back into SoA f64 format.
///
/// Inverse of `flatten_soa_to_f32`: `f[q][cell] = flat[q * nc + cell]`.
/// A buffer of the wrong length (e.g. an empty vector from a failed
/// readback) is reported as [`GpuError::ReadBack`] instead of being padded
/// with zeros.
#[cfg(feature = "wgpu-backend")]
fn unflatten_f32_to_soa(flat: &[f32], nc: usize) -> Result<Vec<Vec<f64>>, GpuError> {
    if flat.len() != 19 * nc {
        return Err(GpuError::ReadBack(format!(
            "expected {} population values, got {}",
            19 * nc,
            flat.len()
        )));
    }
    Ok(flat
        .chunks_exact(nc.max(1))
        .take(19)
        .map(|dir| dir.iter().map(|&v| f64::from(v)).collect())
        .collect())
}

/// D3Q19 BGK LBM simulation.
///
/// Population arrays are indexed `f[dir][cell_index]` where
/// `cell_index = x + nx * (y + ny * z)`.
///
/// With the `wgpu-backend` feature (default) [`step`](Self::step) runs the
/// WGSL kernel on the first available GPU adapter and falls back to the CPU
/// reference when none is available; both paths implement the same physics
/// (Guo forcing, halfway bounce-back walls, Ladd moving lid).
pub struct LbmSimulation {
    /// Configuration.
    ///
    /// `tau` and `force_*` may be changed between steps; the GPU parameters
    /// are re-uploaded automatically.  The grid dimensions must not change.
    /// Setting `compute` to [`ComputePreference::CpuOnly`] switches the GPU
    /// path off for good at the next step (see
    /// [`disable_gpu`](Self::disable_gpu)).
    pub config: LbmConfig,
    /// Distribution functions f_i (19 × N_cells).
    ///
    /// Reading is always safe after an observation method
    /// ([`mean_velocity`](Self::mean_velocity), [`populations`](Self::populations),
    /// …) has synchronised the GPU state.  Mutate the populations through
    /// [`populations_mut`](Self::populations_mut), which also schedules the
    /// re-upload to the GPU; writes made directly to this field while the GPU
    /// path is active are not seen by the GPU.
    pub f: Vec<Vec<f64>>,
    /// Temporary buffer for streaming step.
    f_tmp: Vec<Vec<f64>>,
    /// Solid (no-slip) mask: true = bounce-back wall.  Changed through
    /// [`set_solid`](Self::set_solid) / [`set_solid_mask`](Self::set_solid_mask)
    /// so that the GPU copy stays in sync.
    solid: Vec<bool>,
    /// Lid velocity (X component).
    pub lid_vel_x: f64,
    /// Lid velocity (Y component).
    pub lid_vel_y: f64,
    /// Lid velocity (Z component).
    pub lid_vel_z: f64,
    /// Total steps executed.
    pub step_count: u64,
    /// GPU availability / resources (lazily initialised on the first step).
    #[cfg(feature = "wgpu-backend")]
    gpu: GpuSlot,
    /// Set to `true` after a GPU step; cleared when `f` is synchronised from GPU.
    ///
    /// Observation methods (`mean_velocity`, `mean_density`, etc.) call
    /// `sync_from_gpu()` when this flag is set.
    #[cfg(feature = "wgpu-backend")]
    gpu_dirty: bool,
    /// The CPU-side populations were modified and must be uploaded before the
    /// next GPU step.
    #[cfg(feature = "wgpu-backend")]
    f_upload_pending: bool,
    /// The solid mask was modified and must be uploaded before the next GPU step.
    #[cfg(feature = "wgpu-backend")]
    solid_upload_pending: bool,
    /// Number of GPU initialisation attempts (test instrumentation).
    #[cfg(all(test, feature = "wgpu-backend"))]
    gpu_init_attempts: u32,
}

impl LbmSimulation {
    /// Create a new LBM simulation, initialised at rest with density ρ₀.
    pub fn new(config: LbmConfig) -> Self {
        let nc = config.n_cells();
        let rho0 = config.rho0;

        // Initialise all populations at equilibrium for zero velocity
        let mut f = Vec::with_capacity(19);
        let mut f_tmp = Vec::with_capacity(19);
        for &wi in D3Q19_W.iter() {
            f.push(vec![wi * rho0; nc]);
            f_tmp.push(vec![0.0; nc]);
        }

        // Default solid geometry: bounce-back walls on the Y and Z faces only.
        // The X direction is left open so that the periodic streaming (rem_euclid)
        // acts as a true periodic boundary condition there.  This is the canonical
        // Poiseuille-flow / body-force-driven-channel setup: walls perpendicular to
        // Y and Z provide the no-slip surfaces, while X is the flow direction with
        // periodic inflow/outflow.  The solid layer at y = ny − 1 doubles as the
        // moving lid when a lid velocity is set.
        let mut solid = vec![false; nc];
        let (nx, ny, nz) = (config.nx, config.ny, config.nz);
        for z in 0..nz {
            for y in 0..ny {
                for x in 0..nx {
                    let idx = x + nx * (y + ny * z);
                    if y == 0 || y + 1 == ny || z == 0 || z + 1 == nz {
                        solid[idx] = true;
                    }
                }
            }
        }

        #[cfg(feature = "wgpu-backend")]
        let gpu = if config.compute.allows_gpu() {
            GpuSlot::Untried
        } else {
            GpuSlot::Disabled
        };
        Self {
            config,
            f,
            f_tmp,
            solid,
            lid_vel_x: 0.0,
            lid_vel_y: 0.0,
            lid_vel_z: 0.0,
            step_count: 0,
            #[cfg(feature = "wgpu-backend")]
            gpu,
            #[cfg(feature = "wgpu-backend")]
            gpu_dirty: false,
            #[cfg(feature = "wgpu-backend")]
            f_upload_pending: false,
            #[cfg(feature = "wgpu-backend")]
            solid_upload_pending: false,
            #[cfg(all(test, feature = "wgpu-backend"))]
            gpu_init_attempts: 0,
        }
    }

    /// Set the lid (the solid layer on the top face, y = ny − 1) moving at
    /// (ux, uy, uz).
    ///
    /// The lid is applied with Ladd's moving-wall bounce-back to every solid
    /// cell on the plane y = ny − 1.  Tangential components (`ux`, `uz`)
    /// conserve mass exactly; a non-zero wall-normal component `uy` makes the
    /// lid act as a porous wall that injects (uy < 0) or removes (uy > 0) mass.
    pub fn set_lid_velocity(&mut self, ux: f64, uy: f64, uz: f64) {
        self.lid_vel_x = ux;
        self.lid_vel_y = uy;
        self.lid_vel_z = uz;
    }

    /// Solid mask (`true` = bounce-back wall), indexed by
    /// `x + nx * (y + ny * z)`.
    pub fn solid(&self) -> &[bool] {
        &self.solid
    }

    /// Mark cell `(x, y, z)` as solid (`true`) or fluid (`false`).
    ///
    /// # Errors
    ///
    /// Returns [`GpuError::InvalidHandle`] carrying the flat index when the
    /// coordinates lie outside the grid.
    pub fn set_solid(
        &mut self,
        x: usize,
        y: usize,
        z: usize,
        is_solid: bool,
    ) -> Result<(), GpuError> {
        let (nx, ny, nz) = (self.config.nx, self.config.ny, self.config.nz);
        if x >= nx || y >= ny || z >= nz {
            return Err(GpuError::InvalidHandle(x + nx * (y + ny * z)));
        }
        let idx = x + nx * (y + ny * z);
        if let Some(cell) = self.solid.get_mut(idx)
            && *cell != is_solid
        {
            *cell = is_solid;
            self.mark_solid_changed();
        }
        Ok(())
    }

    /// Replace the whole solid mask.
    ///
    /// # Errors
    ///
    /// Returns [`GpuError::InvalidHandle`] carrying the supplied length when
    /// `mask.len()` differs from the number of cells.
    pub fn set_solid_mask(&mut self, mask: Vec<bool>) -> Result<(), GpuError> {
        if mask.len() != self.config.n_cells() {
            return Err(GpuError::InvalidHandle(mask.len()));
        }
        self.solid = mask;
        self.mark_solid_changed();
        Ok(())
    }

    /// Record that the solid mask changed so the GPU copy is refreshed.
    fn mark_solid_changed(&mut self) {
        #[cfg(feature = "wgpu-backend")]
        {
            self.solid_upload_pending = true;
        }
    }

    /// Current populations (synchronised from the GPU when necessary).
    pub fn populations(&mut self) -> &[Vec<f64>] {
        #[cfg(feature = "wgpu-backend")]
        self.sync_from_gpu();
        &self.f
    }

    /// Mutable access to the populations.
    ///
    /// The GPU state is synchronised first, and the modified populations are
    /// uploaded again before the next GPU step.
    pub fn populations_mut(&mut self) -> &mut [Vec<f64>] {
        #[cfg(feature = "wgpu-backend")]
        {
            self.sync_from_gpu();
            self.f_upload_pending = true;
        }
        &mut self.f
    }

    /// True if a real GPU adapter was successfully initialised.
    ///
    /// Returns `false` before the first call to `step()` (GPU state is lazy;
    /// use [`ensure_gpu`](Self::ensure_gpu) to initialise eagerly).
    pub fn has_gpu(&self) -> bool {
        #[cfg(feature = "wgpu-backend")]
        {
            matches!(self.gpu, GpuSlot::Ready(_))
        }
        #[cfg(not(feature = "wgpu-backend"))]
        {
            false
        }
    }

    /// Initialise the GPU path now instead of on the first step.
    ///
    /// Returns `true` when the GPU is active.  A failed attempt is remembered
    /// (see [`gpu_error`](Self::gpu_error)) and never retried.  A simulation
    /// whose [`config.compute`](LbmConfig::compute) is
    /// [`ComputePreference::CpuOnly`] never initialises the GPU.
    pub fn ensure_gpu(&mut self) -> bool {
        if !self.config.compute.allows_gpu() {
            self.disable_gpu();
        }
        #[cfg(feature = "wgpu-backend")]
        {
            if matches!(self.gpu, GpuSlot::Untried) {
                self.try_init_gpu();
            }
        }
        self.has_gpu()
    }

    /// Switch the GPU path off; all further steps run on the CPU reference.
    ///
    /// Pending GPU results are synchronised into [`f`](Self::f) first.  The
    /// switch is permanent and sets [`config.compute`](LbmConfig::compute) to
    /// [`ComputePreference::CpuOnly`] so the configuration stays consistent.
    pub fn disable_gpu(&mut self) {
        self.config.compute = ComputePreference::CpuOnly;
        #[cfg(feature = "wgpu-backend")]
        {
            self.sync_from_gpu();
            if !matches!(self.gpu, GpuSlot::Unavailable(_)) {
                self.gpu = GpuSlot::Disabled;
            }
            self.gpu_dirty = false;
        }
    }

    /// Why the GPU path is unavailable, if it failed.
    ///
    /// `Some` after GPU initialisation failed (e.g. no adapter) or a later
    /// dispatch / readback failed; the simulation then runs on the CPU.
    /// `None` while the GPU is untried, active, explicitly disabled, or when
    /// the crate is built without the `wgpu-backend` feature.
    pub fn gpu_error(&self) -> Option<&GpuError> {
        #[cfg(feature = "wgpu-backend")]
        {
            match &self.gpu {
                GpuSlot::Unavailable(e) => Some(e),
                _ => None,
            }
        }
        #[cfg(not(feature = "wgpu-backend"))]
        {
            None
        }
    }

    /// Advance one LBM step: BGK collision (Guo forcing) + streaming with
    /// bounce-back walls and the moving lid.
    ///
    /// When the `wgpu-backend` feature is enabled the GPU path is used when an
    /// adapter is available and the CPU reference otherwise.  If the GPU fails
    /// mid-run and its state cannot even be read back (device lost), the step
    /// is taken on the CPU from the most recently synchronised populations; the
    /// failure is reported by [`gpu_error`](Self::gpu_error) and by
    /// [`try_step`](Self::try_step).
    ///
    /// With [`ComputePreference::GpuRequired`] a step the GPU cannot take is
    /// *not* replaced by a CPU step: the simulation is left unchanged (use
    /// [`try_step`](Self::try_step) to obtain the error).
    pub fn step(&mut self) {
        if self.try_step().is_err() && self.config.compute.allows_cpu() {
            self.step_cpu();
            self.step_count += 1;
        }
    }

    /// Advance one LBM step, reporting GPU failures that must not be hidden.
    ///
    /// With [`ComputePreference::Auto`] a GPU that is missing, or whose
    /// dispatch fails while its populations can still be read back, is not an
    /// error: the step is taken on the CPU.  [`ComputePreference::CpuOnly`]
    /// always takes the CPU step.
    ///
    /// # Errors
    ///
    /// * [`GpuError::ReadBack`] when a GPU dispatch failed and the GPU
    ///   populations could not be recovered.  No step is taken in that case;
    ///   the GPU path is marked unavailable and [`f`](Self::f) holds the most
    ///   recently synchronised populations.
    /// * With [`ComputePreference::GpuRequired`], the initialisation or
    ///   dispatch error whenever the GPU cannot take the step
    ///   ([`GpuError::BackendInit`] when no adapter is available or the GPU
    ///   path was disabled, the dispatch error otherwise).  The populations
    ///   stay valid and the step counter is not advanced.
    pub fn try_step(&mut self) -> Result<(), GpuError> {
        if !self.config.compute.allows_gpu() {
            self.disable_gpu();
        }
        #[cfg(feature = "wgpu-backend")]
        {
            if !self.try_step_gpu()? {
                if !self.config.compute.allows_cpu() {
                    return Err(self.gpu_required_error());
                }
                self.step_cpu();
            }
        }
        #[cfg(not(feature = "wgpu-backend"))]
        {
            if !self.config.compute.allows_cpu() {
                return Err(GpuError::BackendInit(
                    "compute preference is GpuRequired but the crate was built without the \
                     `wgpu-backend` feature"
                        .to_string(),
                ));
            }
            self.step_cpu();
        }
        self.step_count += 1;
        Ok(())
    }

    /// The error reported when [`ComputePreference::GpuRequired`] cannot be
    /// honoured because the GPU did not take the step.
    #[cfg(feature = "wgpu-backend")]
    fn gpu_required_error(&self) -> GpuError {
        match &self.gpu {
            GpuSlot::Unavailable(e) => e.clone(),
            GpuSlot::Disabled => GpuError::BackendInit(
                "compute preference is GpuRequired but the GPU path was disabled".to_string(),
            ),
            GpuSlot::Untried | GpuSlot::Ready(_) => GpuError::BackendInit(
                "compute preference is GpuRequired but the GPU did not take the step".to_string(),
            ),
        }
    }

    // ── GPU step ──────────────────────────────────────────────────────────────

    /// GPU parameter words `[nx, ny, nz, ω, Fx, Fy, Fz, u_lid]` (floats bit-cast).
    #[cfg(feature = "wgpu-backend")]
    fn gpu_params(&self) -> Result<[u32; LBM_PARAM_WORDS], GpuError> {
        let dim = |n: usize| {
            u32::try_from(n).map_err(|_| {
                GpuError::BackendInit(format!("LBM grid dimension {n} exceeds u32 range"))
            })
        };
        let omega = (1.0 / self.config.tau) as f32;
        Ok([
            dim(self.config.nx)?,
            dim(self.config.ny)?,
            dim(self.config.nz)?,
            omega.to_bits(),
            (self.config.force_x as f32).to_bits(),
            (self.config.force_y as f32).to_bits(),
            (self.config.force_z as f32).to_bits(),
            (self.lid_vel_x as f32).to_bits(),
            (self.lid_vel_y as f32).to_bits(),
            (self.lid_vel_z as f32).to_bits(),
        ])
    }

    /// Attempt to create the GPU state exactly once; failures are remembered.
    #[cfg(feature = "wgpu-backend")]
    fn try_init_gpu(&mut self) {
        #[cfg(test)]
        {
            self.gpu_init_attempts += 1;
        }
        self.gpu = match LbmGpuState::try_new(self) {
            Ok(state) => {
                self.f_upload_pending = false;
                self.solid_upload_pending = false;
                GpuSlot::Ready(Box::new(state))
            }
            Err(e) => GpuSlot::Unavailable(e),
        };
    }

    /// Run one step on the GPU.
    ///
    /// Returns `Ok(true)` when the GPU took the step, `Ok(false)` when the
    /// caller must take it on the CPU (GPU missing, disabled, or a dispatch
    /// failed after the populations were recovered), and `Err` when the GPU
    /// failed and its populations could not be recovered.
    #[cfg(feature = "wgpu-backend")]
    fn try_step_gpu(&mut self) -> Result<bool, GpuError> {
        if matches!(self.gpu, GpuSlot::Untried) {
            self.try_init_gpu();
        }
        let params = self.gpu_params();
        let nc = self.config.n_cells();
        let GpuSlot::Ready(state) = &mut self.gpu else {
            return Ok(false);
        };
        let dispatched = params.and_then(|params| {
            if self.f_upload_pending {
                state.upload_populations(&self.f, nc);
            }
            if self.solid_upload_pending {
                state.upload_solid(&self.solid);
            }
            state.upload_params_if_changed(params);
            state.dispatch()
        });
        match dispatched {
            Ok(()) => {
                self.f_upload_pending = false;
                self.solid_upload_pending = false;
                self.gpu_dirty = true;
                Ok(true)
            }
            Err(dispatch_err) => {
                // Recover the last valid GPU state (the unswapped input
                // buffer) before handing over to the CPU for good.
                let recovered = if self.gpu_dirty {
                    state.read_populations(nc).map(Some)
                } else {
                    Ok(None)
                };
                match recovered {
                    Ok(populations) => {
                        if let Some(f) = populations {
                            self.f = f;
                        }
                        self.gpu_dirty = false;
                        self.gpu = GpuSlot::Unavailable(dispatch_err);
                        Ok(false)
                    }
                    Err(read_err) => {
                        self.gpu_dirty = false;
                        self.gpu = GpuSlot::Unavailable(GpuError::ReadBack(format!(
                            "{dispatch_err}; recovering the GPU populations failed: {read_err}"
                        )));
                        Err(GpuError::ReadBack(read_err.to_string()))
                    }
                }
            }
        }
    }

    /// Synchronise CPU SoA populations from the current GPU input buffer.
    ///
    /// Called lazily by observation methods when `gpu_dirty` is `true`.  A
    /// failed readback marks the GPU path unavailable (see
    /// [`gpu_error`](Self::gpu_error)); `f` then keeps the most recently
    /// synchronised populations.
    #[cfg(feature = "wgpu-backend")]
    fn sync_from_gpu(&mut self) {
        if !self.gpu_dirty {
            return;
        }
        self.gpu_dirty = false;
        let nc = self.config.n_cells();
        let GpuSlot::Ready(state) = &self.gpu else {
            return;
        };
        match state.read_populations(nc) {
            Ok(f) => self.f = f,
            Err(e) => self.gpu = GpuSlot::Unavailable(e),
        }
    }

    // ── CPU step ──────────────────────────────────────────────────────────────

    /// CPU reference step (f64): BGK collision with Guo forcing, then push
    /// streaming with halfway bounce-back and the Ladd moving lid.
    ///
    /// This is the exact algorithm of the WGSL kernel.
    fn step_cpu(&mut self) {
        let nc = self.config.n_cells();
        let nx = self.config.nx;
        let ny = self.config.ny;
        let nz = self.config.nz;
        let omega = 1.0 / self.config.tau;
        let force = [
            self.config.force_x,
            self.config.force_y,
            self.config.force_z,
        ];
        let lid = [self.lid_vel_x, self.lid_vel_y, self.lid_vel_z];
        let guo = 1.0 - 0.5 * omega;

        // ── Collision (BGK + Guo forcing) ────────────────────────────────────
        // `cell_rho` keeps the local density for the moving-wall term.
        let mut cell_rho = vec![0.0_f64; nc];
        for (cell, rho_out) in cell_rho.iter_mut().enumerate() {
            if self.solid[cell] {
                continue;
            }

            let mut rho = 0.0_f64;
            let mut mom = [0.0_f64; 3];
            for i in 0..19 {
                let fi = self.f[i][cell];
                rho += fi;
                mom[0] += fi * D3Q19_EX[i] as f64;
                mom[1] += fi * D3Q19_EY[i] as f64;
                mom[2] += fi * D3Q19_EZ[i] as f64;
            }
            *rho_out = rho;
            // Degenerate (empty) cells skip the collision and only stream.
            if rho < 1e-10 {
                continue;
            }
            // Guo (2002): equilibrium at u = (Σ fᵢ eᵢ + F/2) / ρ ...
            let ux = (mom[0] + 0.5 * force[0]) / rho;
            let uy = (mom[1] + 0.5 * force[1]) / rho;
            let uz = (mom[2] + 0.5 * force[2]) / rho;
            let u2 = ux * ux + uy * uy + uz * uz;

            for i in 0..19 {
                let (ex, ey, ez) = (D3Q19_EX[i] as f64, D3Q19_EY[i] as f64, D3Q19_EZ[i] as f64);
                let eu = ex * ux + ey * uy + ez * uz;
                let feq = D3Q19_W[i] * rho * (1.0 + 3.0 * eu + 4.5 * eu * eu - 1.5 * u2);
                // ... plus the source Sᵢ = (1 − ω/2) wᵢ [3(eᵢ − u) + 9(eᵢ·u) eᵢ]·F.
                let e_minus_u_f =
                    (ex - ux) * force[0] + (ey - uy) * force[1] + (ez - uz) * force[2];
                let e_f = ex * force[0] + ey * force[1] + ez * force[2];
                let source = guo * D3Q19_W[i] * (3.0 * e_minus_u_f + 9.0 * eu * e_f);
                self.f[i][cell] += omega * (feq - self.f[i][cell]) + source;
            }
        }

        // ── Streaming (push) with halfway bounce-back ────────────────────────
        // Solid cells carry their populations through unchanged.  Every fluid
        // slot of `f_tmp` is written exactly once below (streaming is a
        // bijection per direction, and a bounce-back slot is only reachable
        // from a solid cell, which never streams).
        for i in 0..19 {
            for cell in 0..nc {
                self.f_tmp[i][cell] = if self.solid[cell] {
                    self.f[i][cell]
                } else {
                    0.0
                };
            }
        }

        let (nx2, ny2, nz2) = (nx as i64, ny as i64, nz as i64);
        for z in 0..nz {
            for y in 0..ny {
                for x in 0..nx {
                    let src = x + nx * (y + ny * z);
                    if self.solid[src] {
                        continue;
                    }
                    for i in 0..19 {
                        let xd = (x as i64 + D3Q19_EX[i] as i64).rem_euclid(nx2) as usize;
                        let yd = (y as i64 + D3Q19_EY[i] as i64).rem_euclid(ny2) as usize;
                        let zd = (z as i64 + D3Q19_EZ[i] as i64).rem_euclid(nz2) as usize;
                        let dst = xd + nx * (yd + ny * zd);

                        if self.solid[dst] {
                            let mut value = self.f[i][src];
                            if yd + 1 == ny {
                                // Ladd moving-wall bounce-back on the lid plane:
                                // f_opp(x, t+1) = f*_i(x, t) − 6 wᵢ ρ (eᵢ·u_lid).
                                let e_u = D3Q19_EX[i] as f64 * lid[0]
                                    + D3Q19_EY[i] as f64 * lid[1]
                                    + D3Q19_EZ[i] as f64 * lid[2];
                                value -= 6.0 * D3Q19_W[i] * cell_rho[src] * e_u;
                            }
                            self.f_tmp[D3Q19_OPP[i]][src] = value;
                        } else {
                            self.f_tmp[i][dst] = self.f[i][src];
                        }
                    }
                }
            }
        }

        std::mem::swap(&mut self.f, &mut self.f_tmp);
    }

    // ── Macro quantities ──────────────────────────────────────────────────────

    /// Density and physical velocity of the flat cell index `cell`.
    ///
    /// For fluid cells the velocity includes the Guo half-force correction
    /// u = (Σ fᵢ eᵢ + F/2)/ρ; solid cells (which feel no body force) report
    /// Σ fᵢ eᵢ / ρ.  Out-of-range indices yield zero density and velocity.
    fn macro_at(&self, cell: usize) -> (f64, [f64; 3]) {
        let mut rho = 0.0_f64;
        let mut mom = [0.0_f64; 3];
        for (i, dir) in self.f.iter().enumerate().take(19) {
            let Some(&fi) = dir.get(cell) else {
                return (0.0, [0.0; 3]);
            };
            rho += fi;
            mom[0] += fi * D3Q19_EX[i] as f64;
            mom[1] += fi * D3Q19_EY[i] as f64;
            mom[2] += fi * D3Q19_EZ[i] as f64;
        }
        if rho <= 1e-10 {
            return (rho, [0.0; 3]);
        }
        if !self.solid.get(cell).copied().unwrap_or(true) {
            mom[0] += 0.5 * self.config.force_x;
            mom[1] += 0.5 * self.config.force_y;
            mom[2] += 0.5 * self.config.force_z;
        }
        (rho, [mom[0] / rho, mom[1] / rho, mom[2] / rho])
    }

    /// Density and velocity at cell `(x, y, z)`.
    ///
    /// The velocity is the physical Guo-scheme velocity
    /// `(Σ fᵢ eᵢ + F/2) / ρ` for fluid cells.  Coordinates outside the grid
    /// yield zero density and velocity.
    pub fn cell_macro(&mut self, x: usize, y: usize, z: usize) -> (f64, [f64; 3]) {
        #[cfg(feature = "wgpu-backend")]
        self.sync_from_gpu();
        if x >= self.config.nx || y >= self.config.ny || z >= self.config.nz {
            return (0.0, [0.0; 3]);
        }
        self.macro_at(x + self.config.nx * (y + self.config.ny * z))
    }

    /// Per-cell density ρ = Σᵢ fᵢ (length `nx * ny * nz`).
    pub fn density_field(&mut self) -> Vec<f64> {
        #[cfg(feature = "wgpu-backend")]
        self.sync_from_gpu();
        (0..self.config.n_cells())
            .map(|cell| self.macro_at(cell).0)
            .collect()
    }

    /// Per-cell physical velocity (length `nx * ny * nz`), see
    /// [`cell_macro`](Self::cell_macro).
    pub fn velocity_field(&mut self) -> Vec<[f64; 3]> {
        #[cfg(feature = "wgpu-backend")]
        self.sync_from_gpu();
        (0..self.config.n_cells())
            .map(|cell| self.macro_at(cell).1)
            .collect()
    }

    /// Mean velocity (ux, uy, uz) averaged over all fluid cells.
    ///
    /// When body forces are active the physical (observable) velocity in the Guo
    /// forcing scheme is `u_phys = (Σ f_i e_i + F/2) / ρ`, so `config.force_*`
    /// contributes a half-force correction here.
    ///
    /// If GPU steps have been taken, the CPU-side populations are synchronised
    /// from the GPU before computing the mean.
    pub fn mean_velocity(&mut self) -> (f64, f64, f64) {
        #[cfg(feature = "wgpu-backend")]
        self.sync_from_gpu();
        let mut sum = [0.0_f64; 3];
        let mut count = 0_usize;
        for cell in 0..self.config.n_cells() {
            if self.solid[cell] {
                continue;
            }
            let (_, u) = self.macro_at(cell);
            sum[0] += u[0];
            sum[1] += u[1];
            sum[2] += u[2];
            count += 1;
        }
        if count == 0 {
            return (0.0, 0.0, 0.0);
        }
        let n = count as f64;
        (sum[0] / n, sum[1] / n, sum[2] / n)
    }

    /// Mean density over all fluid cells.
    ///
    /// If GPU steps have been taken, the CPU-side populations are synchronised
    /// from the GPU before computing the mean.
    pub fn mean_density(&mut self) -> f64 {
        #[cfg(feature = "wgpu-backend")]
        self.sync_from_gpu();
        let nc = self.config.n_cells();
        let (sum, count) = (0..nc)
            .filter(|&i| !self.solid[i])
            .map(|i| (0..19_usize).map(|d| self.f[d][i]).sum::<f64>())
            .fold((0.0_f64, 0_usize), |(s, c), rho| (s + rho, c + 1));
        if count == 0 { 0.0 } else { sum / count as f64 }
    }

    /// Maximum physical velocity magnitude across all fluid cells.
    ///
    /// If GPU steps have been taken, the CPU-side populations are synchronised
    /// from the GPU before computing the maximum.
    pub fn max_velocity_magnitude(&mut self) -> f64 {
        #[cfg(feature = "wgpu-backend")]
        self.sync_from_gpu();
        (0..self.config.n_cells())
            .filter(|&cell| !self.solid[cell])
            .map(|cell| {
                let (_, u) = self.macro_at(cell);
                (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]).sqrt()
            })
            .fold(0.0_f64, f64::max)
    }
}

// ── LbmGpuSolver ──────────────────────────────────────────────────────────────

/// GPU-accelerated D3Q19 BGK LBM solver with an `f32`-centric interface.
///
/// A thin wrapper around [`LbmSimulation`]: [`new`](Self::new) initialises
/// the GPU path eagerly (dispatching `lbm_bgk_d3q19.wgsl` with ping-pong
/// buffers) and falls back to the CPU reference when no adapter is
/// available; [`new_cpu`](Self::new_cpu) always uses the CPU reference.  Both
/// paths simulate the same physics: the default channel geometry (bounce-back
/// walls on the Y and Z faces, periodic in X), the optional moving lid on the
/// top face and Guo body forcing.
pub struct LbmGpuSolver {
    /// Number of cells in X.
    pub nx: u32,
    /// Number of cells in Y.
    pub ny: u32,
    /// Number of cells in Z.
    pub nz: u32,
    /// BGK relaxation frequency ω = 1/τ.
    pub omega: f32,
    /// Number of steps executed so far.
    pub step_count: u64,
    /// Underlying simulation (owns the GPU resources when active).
    sim: LbmSimulation,
}

impl LbmGpuSolver {
    /// Build the configuration shared by both constructors.
    fn config(nx: u32, ny: u32, nz: u32, omega: f32, compute: ComputePreference) -> LbmConfig {
        LbmConfig {
            nx: nx as usize,
            ny: ny as usize,
            nz: nz as usize,
            tau: if omega > 0.0 {
                1.0 / f64::from(omega)
            } else {
                0.6
            },
            compute,
            ..LbmConfig::default()
        }
    }

    /// Create a new solver.  Attempts to use a real GPU adapter; if none is
    /// available, falls back to the CPU reference (the reason is available
    /// from [`gpu_error`](Self::gpu_error)).
    ///
    /// The simulation keeps [`ComputePreference::Auto`] while it runs on the
    /// GPU.  When the GPU cannot be initialised the CPU reference is final, so
    /// [`ComputePreference::CpuOnly`] is pinned and nothing probes wgpu again.
    ///
    /// `omega = 1/tau` (e.g. 1.5 for τ = 2/3, giving ν = 1/18 in lattice units).
    pub fn new(nx: u32, ny: u32, nz: u32, omega: f32) -> Self {
        let mut sim = LbmSimulation::new(Self::config(nx, ny, nz, omega, ComputePreference::Auto));
        if !sim.ensure_gpu() {
            // The failure stays observable through `gpu_error`.
            sim.disable_gpu();
        }
        Self {
            nx,
            ny,
            nz,
            omega,
            step_count: 0,
            sim,
        }
    }

    /// Create directly with the CPU reference (useful for testing without GPU).
    ///
    /// The simulation is built with [`ComputePreference::CpuOnly`]: no GPU
    /// adapter is ever requested and every observation reads the CPU state.
    pub fn new_cpu(nx: u32, ny: u32, nz: u32, omega: f32) -> Self {
        let sim = LbmSimulation::new(Self::config(nx, ny, nz, omega, ComputePreference::CpuOnly));
        Self {
            nx,
            ny,
            nz,
            omega,
            step_count: 0,
            sim,
        }
    }

    /// Returns `true` if a real GPU backend is active.
    pub fn is_gpu(&self) -> bool {
        self.sim.has_gpu()
    }

    /// Why the GPU path is unavailable, if it failed (see
    /// [`LbmSimulation::gpu_error`]).
    pub fn gpu_error(&self) -> Option<&GpuError> {
        self.sim.gpu_error()
    }

    /// Set the lid (top face, y = ny − 1) velocity.
    pub fn set_lid_velocity(&mut self, ux: f64, uy: f64, uz: f64) {
        self.sim.set_lid_velocity(ux, uy, uz);
    }

    /// The underlying simulation.
    pub fn simulation(&self) -> &LbmSimulation {
        &self.sim
    }

    /// Mutable access to the underlying simulation (geometry, forcing, …).
    pub fn simulation_mut(&mut self) -> &mut LbmSimulation {
        &mut self.sim
    }

    /// Advance one BGK collision + streaming step.
    ///
    /// # Errors
    ///
    /// Returns [`GpuError::ReadBack`] when a GPU dispatch failed and the GPU
    /// populations could not be recovered (see [`LbmSimulation::try_step`]).
    /// A dispatch failure whose state can be recovered continues on the CPU
    /// and is reported through [`gpu_error`](Self::gpu_error).
    pub fn step(&mut self) -> Result<(), GpuError> {
        self.sim.try_step()?;
        self.step_count += 1;
        Ok(())
    }

    /// Per-cell density ρ = Σᵢ fᵢ (length `nx * ny * nz`), read back from the
    /// GPU when it is active.
    pub fn read_density(&mut self) -> Vec<f32> {
        self.sim
            .density_field()
            .into_iter()
            .map(|rho| rho as f32)
            .collect()
    }

    /// Per-cell velocity `[ux, uy, uz]` (length `nx * ny * nz`), read back
    /// from the GPU when it is active.
    pub fn read_velocity(&mut self) -> Vec<[f32; 3]> {
        self.sim
            .velocity_field()
            .into_iter()
            .map(|u| [u[0] as f32, u[1] as f32, u[2] as f32])
            .collect()
    }
}

#[cfg(test)]
mod tests;
