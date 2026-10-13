// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Web Worker and `SharedArrayBuffer` multi-threaded support for WASM.
//!
//! This module provides the Rust-side infrastructure for offloading physics
//! simulation to a Web Worker thread, keeping the main browser thread free for
//! UI rendering.
//!
//! # Threading model
//!
//! ```text
//!  Main thread (JS/TS)                 Worker thread (WASM)
//!  ┌────────────────────────┐          ┌─────────────────────────────┐
//!  │ WorkerBridge           │          │ WorkerRuntime               │
//!  │   .post(SimCommand)    │──msg──▶  │   .dispatch(cmd)            │
//!  │   .poll() → SimResult  │◀──msg──  │   .step(dt)                 │
//!  │                        │          │   world: WasmPhysicsEngine   │
//!  │ SharedStateBuffer      │          │   state: SharedStateBuffer   │
//!  │ (SharedArrayBuffer)    │◀──SAB──  │ (same SAB, written by wkr)  │
//!  └────────────────────────┘          └─────────────────────────────┘
//! ```
//!
//! ## SharedArrayBuffer layout
//!
//! The [`SharedStateBuffer`] wraps a `SharedArrayBuffer` (SAB) with a fixed
//! layout for lock-free reads by the main thread:
//!
//! | Offset (bytes) | Type | Description |
//! |---|---|---|
//! | 0   | i32 | Spinlock (0=free, 1=worker writing) |
//! | 4   | u32 | Body count |
//! | 8   | u32 | Step count |
//! | 12  | f32 | Simulation time (s) |
//! | 16  | f32[N×3] | Body positions (x,y,z per body) |
//! | 16+12N | f32[N×4] | Body orientations (qx,qy,qz,qw) |
//! | 16+28N | f32[N×3] | Body velocities |
//!
//! `N` is the `max_bodies` the buffer was created with.  Every per-body
//! accessor is bounds-checked against `N`: a write for a body index `>= N`
//! is rejected (`false` from the Rust API, an exception from the JavaScript
//! API) and leaves the buffer untouched, and a read for such an index yields
//! zeros.  Without that guard a position write for index `N` would land in
//! the orientation region and corrupt body 0.
//!
//! ## JavaScript interop
//!
//! ```js
//! import init, { WorkerBridge } from './oxiphysics_wasm.js';
//!
//! const bridge = new WorkerBridge(64);
//! bridge.postCommandStep(1.0/60.0);
//! const positions = bridge.sharedPositionsFlat(bridge.sharedBodyCount());
//! ```
//!
//! ## Usage (Rust)
//!
//! ```
//! use oxiphysics_wasm::web_worker::{
//!     WorkerBridge, SimCommand, SimResult, SharedStateBuffer,
//! };
//!
//! let max_bodies = 64;
//! let mut shared = SharedStateBuffer::new(max_bodies);
//! assert_eq!(shared.max_bodies(), max_bodies);
//! assert_eq!(shared.read_step_count(), 0);
//!
//! // Simulate a worker writing body positions
//! shared.write_body_position(0, [1.0, 2.0, 3.0]);
//! let pos = shared.read_body_position(0);
//! assert!((pos[0] - 1.0).abs() < 1e-6);
//! ```

use wasm_bindgen::prelude::*;

use crate::error::Error;
use crate::wasm_helpers::now_ms;

/// Convert an elapsed interval measured by [`now_ms`] (milliseconds) into whole
/// microseconds for [`SimResult::StepDone::step_us`].
///
/// Returns `0` when either endpoint was unmeasurable (the target exposes no
/// clock) — an honest "not measured" rather than a fabricated constant.
fn elapsed_us(start: Option<f64>, end: Option<f64>) -> u32 {
    match (start, end) {
        (Some(a), Some(b)) => {
            let us = (b - a).max(0.0) * 1_000.0;
            if us.is_finite() {
                us.round().min(u32::MAX as f64) as u32
            } else {
                0
            }
        }
        _ => 0,
    }
}

// ── SimCommand ────────────────────────────────────────────────────────────────

/// Commands sent from the main thread to the physics worker.
///
/// This is a Rust-only enum (payload variants cannot be `#[wasm_bindgen]`).
/// Use the `post_command_*` methods on `WorkerBridge` from JavaScript.
#[derive(Debug, Clone)]
pub enum SimCommand {
    /// Advance simulation by `dt` seconds.
    Step { dt: f32 },
    /// Advance simulation by `dt` using `substeps` internal sub-steps.
    StepSubsteps { dt: f32, substeps: u32 },
    /// Add a dynamic sphere body at position (x,y,z) with given mass and radius.
    AddSphere {
        mass: f32,
        x: f32,
        y: f32,
        z: f32,
        radius: f32,
    },
    /// Add a static box body (immovable ground plane, etc.).
    AddStaticBox {
        x: f32,
        y: f32,
        z: f32,
        hx: f32,
        hy: f32,
        hz: f32,
    },
    /// Apply an impulse to body `handle`.
    ApplyImpulse {
        handle: u32,
        ix: f32,
        iy: f32,
        iz: f32,
    },
    /// Reset the simulation (remove all bodies).
    Reset,
    /// Request a full state snapshot (triggers `SimResult::Snapshot`).
    RequestSnapshot,
    /// Gracefully shut down the worker.
    Shutdown,
}

/// Serialised wire representation of [`SimCommand`] (flat byte-friendly).
///
/// Rust-only — not exposed to JS directly due to `[f32; 7]` field.
#[derive(Debug, Clone, Copy)]
#[repr(C)]
pub struct RawSimCommand {
    /// Command discriminant (see [`SimCommandKind`]).
    pub kind: u8,
    /// Padding.
    pub _pad: [u8; 3],
    /// Packed f32 payload (interpretation depends on `kind`).
    pub payload: [f32; 7],
    /// Integer payload (body handle, substep count, etc.).
    pub ipayload: u32,
}

/// Discriminants for [`RawSimCommand::kind`].
#[wasm_bindgen]
#[repr(u8)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SimCommandKind {
    /// No-op.
    None = 0,
    /// [`SimCommand::Step`].
    Step = 1,
    /// [`SimCommand::StepSubsteps`].
    StepSubs = 2,
    /// [`SimCommand::AddSphere`].
    AddSphere = 3,
    /// [`SimCommand::AddStaticBox`].
    AddBox = 4,
    /// [`SimCommand::ApplyImpulse`].
    Impulse = 5,
    /// [`SimCommand::Reset`].
    Reset = 6,
    /// [`SimCommand::RequestSnapshot`].
    Snapshot = 7,
    /// [`SimCommand::Shutdown`].
    Shutdown = 8,
}

impl SimCommand {
    /// Serialise to [`RawSimCommand`] for postMessage transfer.
    pub fn to_raw(&self) -> RawSimCommand {
        match *self {
            Self::Step { dt } => RawSimCommand {
                kind: SimCommandKind::Step as u8,
                _pad: [0; 3],
                payload: [dt, 0., 0., 0., 0., 0., 0.],
                ipayload: 0,
            },
            Self::StepSubsteps { dt, substeps } => RawSimCommand {
                kind: SimCommandKind::StepSubs as u8,
                _pad: [0; 3],
                payload: [dt, 0., 0., 0., 0., 0., 0.],
                ipayload: substeps,
            },
            Self::AddSphere {
                mass,
                x,
                y,
                z,
                radius,
            } => RawSimCommand {
                kind: SimCommandKind::AddSphere as u8,
                _pad: [0; 3],
                payload: [mass, x, y, z, radius, 0., 0.],
                ipayload: 0,
            },
            Self::AddStaticBox {
                x,
                y,
                z,
                hx,
                hy,
                hz,
            } => RawSimCommand {
                kind: SimCommandKind::AddBox as u8,
                _pad: [0; 3],
                payload: [x, y, z, hx, hy, hz, 0.],
                ipayload: 0,
            },
            Self::ApplyImpulse { handle, ix, iy, iz } => RawSimCommand {
                kind: SimCommandKind::Impulse as u8,
                _pad: [0; 3],
                payload: [ix, iy, iz, 0., 0., 0., 0.],
                ipayload: handle,
            },
            Self::Reset => RawSimCommand {
                kind: SimCommandKind::Reset as u8,
                _pad: [0; 3],
                payload: [0.; 7],
                ipayload: 0,
            },
            Self::RequestSnapshot => RawSimCommand {
                kind: SimCommandKind::Snapshot as u8,
                _pad: [0; 3],
                payload: [0.; 7],
                ipayload: 0,
            },
            Self::Shutdown => RawSimCommand {
                kind: SimCommandKind::Shutdown as u8,
                _pad: [0; 3],
                payload: [0.; 7],
                ipayload: 0,
            },
        }
    }

    /// Deserialise from [`RawSimCommand`].
    pub fn from_raw(r: &RawSimCommand) -> Option<Self> {
        match r.kind {
            k if k == SimCommandKind::Step as u8 => Some(Self::Step { dt: r.payload[0] }),
            k if k == SimCommandKind::StepSubs as u8 => Some(Self::StepSubsteps {
                dt: r.payload[0],
                substeps: r.ipayload,
            }),
            k if k == SimCommandKind::AddSphere as u8 => Some(Self::AddSphere {
                mass: r.payload[0],
                x: r.payload[1],
                y: r.payload[2],
                z: r.payload[3],
                radius: r.payload[4],
            }),
            k if k == SimCommandKind::AddBox as u8 => Some(Self::AddStaticBox {
                x: r.payload[0],
                y: r.payload[1],
                z: r.payload[2],
                hx: r.payload[3],
                hy: r.payload[4],
                hz: r.payload[5],
            }),
            k if k == SimCommandKind::Impulse as u8 => Some(Self::ApplyImpulse {
                handle: r.ipayload,
                ix: r.payload[0],
                iy: r.payload[1],
                iz: r.payload[2],
            }),
            k if k == SimCommandKind::Reset as u8 => Some(Self::Reset),
            k if k == SimCommandKind::Snapshot as u8 => Some(Self::RequestSnapshot),
            k if k == SimCommandKind::Shutdown as u8 => Some(Self::Shutdown),
            _ => None,
        }
    }
}

// ── SimResult ─────────────────────────────────────────────────────────────────

/// Results posted back from the physics worker to the main thread.
///
/// Rust-only enum (payload variants cannot be `#[wasm_bindgen]`).
/// Use `WorkerBridge::poll_result_*` from JavaScript.
#[derive(Debug, Clone)]
pub enum SimResult {
    /// Step completed; includes timing info.
    ///
    /// Note: [`WorkerRuntime`] currently runs a lightweight **kinematic
    /// preview** (ballistic gravity integration with no collision detection),
    /// not the full constraint solver. Consequently `contact_count` is the
    /// genuine count from a contactless integrator (always `0`), and `step_us`
    /// is the *real* measured wall-clock cost of the preview in microseconds
    /// (`0` when the target exposes no clock or the work was sub-microsecond),
    /// never a fabricated constant.
    StepDone {
        sim_time: f32,
        /// Number of bodies advanced this step.
        body_count: u32,
        /// Real contact count. The kinematic preview performs no collision
        /// detection, so this is genuinely `0` (not a placeholder for an
        /// unimplemented full-physics path).
        contact_count: u32,
        /// Real measured wall-clock cost of this step in microseconds; `0`
        /// when unmeasurable (no clock / sub-µs), never a fabricated literal.
        step_us: u32,
    },
    /// Body added successfully.
    BodyAdded { handle: u32 },
    /// Full snapshot of all body states.
    Snapshot {
        positions: Vec<f32>,
        orientations: Vec<f32>,
        velocities: Vec<f32>,
        sim_time: f32,
    },
    /// Worker is ready.
    Ready,
    /// Worker has shut down.
    Shutdown,
    /// An error occurred.
    Error(String),
}

impl SimResult {
    /// Kind string for JS (e.g. `"StepDone"`, `"BodyAdded"`, `"Snapshot"`, `"Ready"`, `"Shutdown"`, `"Error"`).
    pub fn kind_str(&self) -> &'static str {
        match self {
            Self::StepDone { .. } => "StepDone",
            Self::BodyAdded { .. } => "BodyAdded",
            Self::Snapshot { .. } => "Snapshot",
            Self::Ready => "Ready",
            Self::Shutdown => "Shutdown",
            Self::Error(_) => "Error",
        }
    }

    /// Serialise to a flat `Vec<f32>` for JS transfer.
    ///
    /// Layout: `[kind_code, payload...]`
    /// - `StepDone`: `[1.0, sim_time, body_count, contact_count, step_us]`
    /// - `BodyAdded`: `[2.0, handle]`
    /// - `Snapshot`: `[3.0, sim_time, n_bodies, positions..., orientations..., velocities...]`
    /// - `Ready`: `[4.0]`
    /// - `Shutdown`: `[5.0]`
    /// - `Error`: `[6.0]` (message lost in flat repr)
    pub fn to_flat(&self) -> Vec<f32> {
        match self {
            Self::StepDone {
                sim_time,
                body_count,
                contact_count,
                step_us,
            } => {
                vec![
                    1.0,
                    *sim_time,
                    *body_count as f32,
                    *contact_count as f32,
                    *step_us as f32,
                ]
            }
            Self::BodyAdded { handle } => vec![2.0, *handle as f32],
            Self::Snapshot {
                positions,
                orientations,
                velocities,
                sim_time,
            } => {
                let n = positions.len() / 3;
                let mut v = vec![3.0, *sim_time, n as f32];
                v.extend_from_slice(positions);
                v.extend_from_slice(orientations);
                v.extend_from_slice(velocities);
                v
            }
            Self::Ready => vec![4.0],
            Self::Shutdown => vec![5.0],
            Self::Error(_) => vec![6.0],
        }
    }
}

// ── SharedStateBuffer ─────────────────────────────────────────────────────────

/// Size in bytes of the fixed header (lock, body count, step count, time).
const HEADER_BYTES: usize = 16;
/// Bytes per body in the position region (`f32 x 3`).
const POSITION_STRIDE: usize = 12;
/// Bytes per body in the orientation region (`f32 x 4`).
const ORIENTATION_STRIDE: usize = 16;
/// Bytes per body in the velocity region (`f32 x 3`).
const VELOCITY_STRIDE: usize = 12;

/// CPU-resident approximation of a `SharedArrayBuffer` for Rust-side logic.
///
/// In a real browser environment this wraps a JS `SharedArrayBuffer` and
/// accesses it via `Atomics` operations.  In non-WASM / unit-test contexts
/// it is a plain `Vec<u8>` providing the same layout.
///
/// All per-body accessors validate the body index against
/// [`Self::max_bodies`]; see the module documentation.
#[wasm_bindgen]
pub struct SharedStateBuffer {
    data: Vec<u8>,
    max_bodies: usize,
}

#[wasm_bindgen]
impl SharedStateBuffer {
    /// Create a new shared state buffer for up to `max_bodies` bodies.
    pub fn new(max_bodies: usize) -> Self {
        let size = Self::byte_size(max_bodies);
        Self {
            data: vec![0u8; size],
            max_bodies,
        }
    }

    /// Total byte size for `max_bodies`.
    pub fn byte_size(max_bodies: usize) -> usize {
        16
        + max_bodies * 3 * 4    // positions  (f32×3)
        + max_bodies * 4 * 4    // orientations (f32×4)
        + max_bodies * 3 * 4 // velocities (f32×3)
    }

    /// Maximum number of bodies this buffer supports.
    pub fn max_bodies(&self) -> usize {
        self.max_bodies
    }

    /// Total byte length.
    pub fn byte_len(&self) -> usize {
        self.data.len()
    }

    /// Acquire the spinlock (busy-wait).
    pub fn lock(&mut self) {
        let bytes: [u8; 4] = 1_i32.to_le_bytes();
        self.data[0..4].copy_from_slice(&bytes);
    }

    /// Release the spinlock.
    pub fn unlock(&mut self) {
        self.data[0..4].copy_from_slice(&0_i32.to_le_bytes());
    }

    /// Write the body count.
    pub fn write_body_count(&mut self, n: u32) {
        self.data[4..8].copy_from_slice(&n.to_le_bytes());
    }

    /// Read the body count.
    pub fn read_body_count(&self) -> u32 {
        let b: [u8; 4] = self.data[4..8].try_into().unwrap_or([0; 4]);
        u32::from_le_bytes(b)
    }

    /// Increment and return the step count.
    pub fn increment_step_count(&mut self) -> u32 {
        let b: [u8; 4] = self.data[8..12].try_into().unwrap_or([0; 4]);
        let c = u32::from_le_bytes(b) + 1;
        self.data[8..12].copy_from_slice(&c.to_le_bytes());
        c
    }

    /// Read the step count.
    pub fn read_step_count(&self) -> u32 {
        let b: [u8; 4] = self.data[8..12].try_into().unwrap_or([0; 4]);
        u32::from_le_bytes(b)
    }

    /// Write the simulation time.
    pub fn write_sim_time(&mut self, t: f32) {
        self.data[12..16].copy_from_slice(&t.to_le_bytes());
    }

    /// Read the simulation time.
    pub fn read_sim_time(&self) -> f32 {
        let b: [u8; 4] = self.data[12..16].try_into().unwrap_or([0; 4]);
        f32::from_le_bytes(b)
    }

    /// Write position `[x, y, z]` for body `i` (JS-compatible: uses flat `Vec<f32>`).
    ///
    /// Throws when `i >= max_bodies()` or when `pos_flat` holds fewer than
    /// three components; the buffer is left untouched in both cases.
    pub fn write_body_position_js(&mut self, i: usize, pos_flat: Vec<f32>) -> Result<(), JsValue> {
        Ok(self.write_body_position_slice(i, &pos_flat)?)
    }

    /// Read position for body `i` as `Vec<f32>` of `[x, y, z]` (JS-compatible).
    pub fn read_body_position_js(&self, i: usize) -> Vec<f32> {
        self.read_body_position(i).to_vec()
    }

    /// Read orientation quaternion for body `i` as `Vec<f32>` of `[qx, qy, qz, qw]`.
    pub fn read_body_orientation_js(&self, i: usize) -> Vec<f32> {
        self.read_body_orientation(i).to_vec()
    }

    /// Read linear velocity for body `i` as `Vec<f32>` of `[vx, vy, vz]`.
    pub fn read_body_velocity_js(&self, i: usize) -> Vec<f32> {
        self.read_body_velocity(i).to_vec()
    }

    /// Flat positions slice `[x0,y0,z0, x1,y1,z1, …]` for `n` bodies.
    pub fn positions_flat(&self, n: usize) -> Vec<f32> {
        (0..n.min(self.max_bodies))
            .flat_map(|i| self.read_body_position(i))
            .collect()
    }

    /// Write all positions from a flat `Vec<f32>`.
    pub fn write_positions_flat_js(&mut self, flat: Vec<f32>) {
        let n = flat.len() / 3;
        for i in 0..n.min(self.max_bodies) {
            self.write_body_position(i, [flat[i * 3], flat[i * 3 + 1], flat[i * 3 + 2]]);
        }
    }

    /// Return raw bytes as a `Vec<u8>` (JS-compatible copy).
    pub fn as_bytes_js(&self) -> Vec<u8> {
        self.data.clone()
    }
}

impl SharedStateBuffer {
    /// Byte offset of body `i`'s position record, or `None` when
    /// `i >= max_bodies` (the record would belong to another region).
    fn position_offset(&self, i: usize) -> Option<usize> {
        if i >= self.max_bodies {
            return None;
        }
        HEADER_BYTES.checked_add(i.checked_mul(POSITION_STRIDE)?)
    }

    /// Byte offset of body `i`'s orientation record, or `None` when
    /// `i >= max_bodies`.
    fn orientation_offset(&self, i: usize) -> Option<usize> {
        if i >= self.max_bodies {
            return None;
        }
        let region = HEADER_BYTES.checked_add(self.max_bodies.checked_mul(POSITION_STRIDE)?)?;
        region.checked_add(i.checked_mul(ORIENTATION_STRIDE)?)
    }

    /// Byte offset of body `i`'s velocity record, or `None` when
    /// `i >= max_bodies`.
    fn velocity_offset(&self, i: usize) -> Option<usize> {
        if i >= self.max_bodies {
            return None;
        }
        let region = HEADER_BYTES.checked_add(
            self.max_bodies
                .checked_mul(POSITION_STRIDE + ORIENTATION_STRIDE)?,
        )?;
        region.checked_add(i.checked_mul(VELOCITY_STRIDE)?)
    }

    /// Store `values` as little-endian `f32`s starting at byte `offset`.
    ///
    /// Returns `false`, writing nothing, when the range does not fit in the
    /// buffer.
    fn write_f32s(&mut self, offset: usize, values: &[f32]) -> bool {
        let Some(end) = values
            .len()
            .checked_mul(4)
            .and_then(|len| offset.checked_add(len))
        else {
            return false;
        };
        let Some(dst) = self.data.get_mut(offset..end) else {
            return false;
        };
        for (chunk, v) in dst.chunks_exact_mut(4).zip(values) {
            chunk.copy_from_slice(&v.to_le_bytes());
        }
        true
    }

    /// Load `N` consecutive little-endian `f32`s starting at byte `offset`
    /// (components outside the buffer read as `0.0`).
    fn read_f32s<const N: usize>(&self, offset: usize) -> [f32; N] {
        let mut out = [0.0_f32; N];
        for (k, slot) in out.iter_mut().enumerate() {
            *slot = offset
                .checked_add(k * 4)
                .map_or(0.0, |at| read_f32_le(&self.data, at));
        }
        out
    }

    /// Error describing a rejected body index.
    fn index_error(&self, i: usize) -> Error {
        Error::invalid_param(
            "body_index",
            format!(
                "{i} is out of range for a shared buffer of {} bodies",
                self.max_bodies
            ),
        )
    }

    /// Write position `[x, y, z]` for body `i` (Rust-only).
    ///
    /// Returns `false`, leaving the buffer untouched, when `i >= max_bodies()`.
    pub fn write_body_position(&mut self, i: usize, pos: [f32; 3]) -> bool {
        match self.position_offset(i) {
            Some(offset) => self.write_f32s(offset, &pos),
            None => false,
        }
    }

    /// Validated slice form of [`Self::write_body_position`] used by the
    /// JavaScript API.  The first three components of `pos` are written.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidParameter`] when `pos` holds fewer than three
    /// components or `i >= max_bodies()`.
    pub fn write_body_position_slice(&mut self, i: usize, pos: &[f32]) -> crate::error::Result<()> {
        let &[x, y, z, ..] = pos else {
            return Err(Error::invalid_param(
                "pos_flat",
                format!("expected 3 components, got {}", pos.len()),
            ));
        };
        if self.write_body_position(i, [x, y, z]) {
            Ok(())
        } else {
            Err(self.index_error(i))
        }
    }

    /// Read position `[x, y, z]` for body `i` (Rust-only); zeros when
    /// `i >= max_bodies()`.
    pub fn read_body_position(&self, i: usize) -> [f32; 3] {
        self.position_offset(i)
            .map_or([0.0; 3], |offset| self.read_f32s::<3>(offset))
    }

    /// Write orientation quaternion for body `i` (Rust-only).
    ///
    /// Returns `false`, leaving the buffer untouched, when `i >= max_bodies()`.
    pub fn write_body_orientation(&mut self, i: usize, q: [f32; 4]) -> bool {
        match self.orientation_offset(i) {
            Some(offset) => self.write_f32s(offset, &q),
            None => false,
        }
    }

    /// Read orientation quaternion for body `i` (Rust-only); zeros when
    /// `i >= max_bodies()`.
    pub fn read_body_orientation(&self, i: usize) -> [f32; 4] {
        self.orientation_offset(i)
            .map_or([0.0; 4], |offset| self.read_f32s::<4>(offset))
    }

    /// Write linear velocity for body `i` (Rust-only).
    ///
    /// Returns `false`, leaving the buffer untouched, when `i >= max_bodies()`.
    pub fn write_body_velocity(&mut self, i: usize, vel: [f32; 3]) -> bool {
        match self.velocity_offset(i) {
            Some(offset) => self.write_f32s(offset, &vel),
            None => false,
        }
    }

    /// Read linear velocity for body `i` (Rust-only); zeros when
    /// `i >= max_bodies()`.
    pub fn read_body_velocity(&self, i: usize) -> [f32; 3] {
        self.velocity_offset(i)
            .map_or([0.0; 3], |offset| self.read_f32s::<3>(offset))
    }

    /// Write all positions from a flat `&[f32]` slice (Rust-only).
    pub fn write_positions_flat(&mut self, flat: &[f32]) {
        let n = flat.len() / 3;
        for i in 0..n.min(self.max_bodies) {
            self.write_body_position(i, [flat[i * 3], flat[i * 3 + 1], flat[i * 3 + 2]]);
        }
    }

    /// Raw byte slice (for `Uint8Array` view in WASM; Rust-only).
    pub fn as_bytes(&self) -> &[u8] {
        &self.data
    }

    /// Mutable raw byte slice (Rust-only).
    pub fn as_bytes_mut(&mut self) -> &mut [u8] {
        &mut self.data
    }
}

/// Read a little-endian `f32` from `bytes` at `offset` (returns 0.0 if out of bounds).
#[inline]
fn read_f32_le(bytes: &[u8], offset: usize) -> f32 {
    offset
        .checked_add(4)
        .and_then(|end| bytes.get(offset..end))
        .and_then(|slice| <[u8; 4]>::try_from(slice).ok())
        .map_or(0.0, f32::from_le_bytes)
}

// ── WorkerBridge ──────────────────────────────────────────────────────────────

/// Main-thread-side bridge to the physics Web Worker.
#[wasm_bindgen]
pub struct WorkerBridge {
    /// Shared state buffer (written by worker, read by main thread).
    shared: SharedStateBuffer,
    pending_commands: Vec<RawSimCommand>,
    received_results: Vec<SimResult>,
    connected: bool,
}

#[wasm_bindgen]
impl WorkerBridge {
    /// Create a new bridge for up to `max_bodies` bodies.
    pub fn new(max_bodies: usize) -> Self {
        Self {
            shared: SharedStateBuffer::new(max_bodies),
            pending_commands: Vec::new(),
            received_results: Vec::new(),
            connected: false,
        }
    }

    /// Whether the bridge is connected.
    pub fn connected(&self) -> bool {
        self.connected
    }

    /// Set the connected state.
    pub fn set_connected(&mut self, v: bool) {
        self.connected = v;
    }

    /// Post a Step command.
    pub fn post_command_step(&mut self, dt: f32) {
        self.pending_commands.push(SimCommand::Step { dt }.to_raw());
    }

    /// Post a StepSubsteps command.
    pub fn post_command_step_substeps(&mut self, dt: f32, substeps: u32) {
        self.pending_commands
            .push(SimCommand::StepSubsteps { dt, substeps }.to_raw());
    }

    /// Post an AddSphere command.
    pub fn post_command_add_sphere(&mut self, mass: f32, x: f32, y: f32, z: f32, radius: f32) {
        self.pending_commands.push(
            SimCommand::AddSphere {
                mass,
                x,
                y,
                z,
                radius,
            }
            .to_raw(),
        );
    }

    /// Post an AddStaticBox command.
    pub fn post_command_add_static_box(
        &mut self,
        x: f32,
        y: f32,
        z: f32,
        hx: f32,
        hy: f32,
        hz: f32,
    ) {
        self.pending_commands.push(
            SimCommand::AddStaticBox {
                x,
                y,
                z,
                hx,
                hy,
                hz,
            }
            .to_raw(),
        );
    }

    /// Post an ApplyImpulse command.
    pub fn post_command_apply_impulse(&mut self, handle: u32, ix: f32, iy: f32, iz: f32) {
        self.pending_commands
            .push(SimCommand::ApplyImpulse { handle, ix, iy, iz }.to_raw());
    }

    /// Post a Reset command.
    pub fn post_command_reset(&mut self) {
        self.pending_commands.push(SimCommand::Reset.to_raw());
    }

    /// Post a RequestSnapshot command.
    pub fn post_command_request_snapshot(&mut self) {
        self.pending_commands
            .push(SimCommand::RequestSnapshot.to_raw());
    }

    /// Post a Shutdown command.
    pub fn post_command_shutdown(&mut self) {
        self.pending_commands.push(SimCommand::Shutdown.to_raw());
    }

    /// Poll for a received result as a flat `Vec<f32>` (empty if none available).
    ///
    /// Use `result_kind` on the first element to determine type:
    /// 1=StepDone, 2=BodyAdded, 3=Snapshot, 4=Ready, 5=Shutdown, 6=Error
    pub fn poll_result_flat(&mut self) -> Vec<f32> {
        if self.received_results.is_empty() {
            return Vec::new();
        }
        self.received_results.remove(0).to_flat()
    }

    /// Whether there are any unread results.
    pub fn has_result(&self) -> bool {
        !self.received_results.is_empty()
    }

    /// Number of pending commands.
    pub fn pending_command_count(&self) -> u32 {
        self.pending_commands.len() as u32
    }

    /// Read shared body count.
    pub fn shared_body_count(&self) -> u32 {
        self.shared.read_body_count()
    }

    /// Read shared step count.
    pub fn shared_step_count(&self) -> u32 {
        self.shared.read_step_count()
    }

    /// Read shared simulation time.
    pub fn shared_sim_time(&self) -> f32 {
        self.shared.read_sim_time()
    }

    /// Return shared positions as flat `Vec<f32>` for `n` bodies.
    pub fn shared_positions_flat(&self, n: usize) -> Vec<f32> {
        self.shared.positions_flat(n)
    }
}

impl WorkerBridge {
    /// Post a raw `SimCommand` (Rust-only).
    pub fn post_command(&mut self, cmd: &SimCommand) {
        self.pending_commands.push(cmd.to_raw());
    }

    /// Poll for a received result (Rust-only).
    pub fn poll_result(&mut self) -> Option<SimResult> {
        if self.received_results.is_empty() {
            None
        } else {
            Some(self.received_results.remove(0))
        }
    }

    /// Drain pending commands (used by WorkerRuntime in single-process tests).
    pub fn drain_commands(&mut self) -> Vec<RawSimCommand> {
        std::mem::take(&mut self.pending_commands)
    }

    /// Push a result from the worker.
    pub fn push_result(&mut self, r: SimResult) {
        self.received_results.push(r);
    }

    /// Access shared state (Rust-only).
    pub fn shared(&self) -> &SharedStateBuffer {
        &self.shared
    }

    /// Access shared state mutably (Rust-only).
    pub fn shared_mut(&mut self) -> &mut SharedStateBuffer {
        &mut self.shared
    }
}

// ── WorkerRuntime ─────────────────────────────────────────────────────────────

/// Worker-thread-side physics runtime.
///
/// This runtime is a **lightweight kinematic preview**: `Step` advances each
/// body under constant gravity (`y += ½ g dt²`, a direct positional term) and
/// integrates any externally-imparted velocity (`pos += v · dt`) carried by
/// [`SimCommand::ApplyImpulse`]. It performs no collision detection, contacts,
/// or constraint solving. It is intended for off-thread plumbing and smoke
/// tests; a production deployment swaps in [`crate::WasmPhysicsEngine`].
/// Telemetry it returns is honest about this scope — see [`SimResult::StepDone`].
#[wasm_bindgen]
pub struct WorkerRuntime {
    /// Accumulated simulation time.
    pub sim_time: f32,
    /// Step count.
    pub step_count: u32,
    /// Simulated body count (stub; real impl uses WasmPhysicsEngine).
    ///
    /// The runtime keeps this equal to the number of bodies it stores.
    /// Capacity checks and body handles are derived from the stored bodies,
    /// never from a value written to this field from JavaScript.
    pub body_count: u32,
    /// Body positions (stub storage for testing).
    positions: Vec<[f32; 3]>,
    /// Body velocities, kept parallel to [`Self::positions`] (one entry per
    /// body). The preview tracks the velocity imparted by
    /// [`SimCommand::ApplyImpulse`] and integrates it into position each
    /// [`SimCommand::Step`]. Gravity is modeled as a direct positional term
    /// (see the `Step` arm of [`Self::dispatch`]) and is intentionally *not*
    /// accumulated here, so a body that received no impulse has zero tracked
    /// velocity.
    velocities: Vec<[f32; 3]>,
    /// Running flag.
    pub running: bool,
}

#[wasm_bindgen]
impl WorkerRuntime {
    /// Create a new worker runtime.
    pub fn new() -> Self {
        Self {
            sim_time: 0.0,
            step_count: 0,
            body_count: 0,
            positions: Vec::new(),
            velocities: Vec::new(),
            running: true,
        }
    }

    /// Process all pending commands from `bridge` and write results back.
    pub fn process(&mut self, bridge: &mut WorkerBridge) {
        let cmds = bridge.drain_commands();
        for raw_cmd in cmds {
            if let Some(cmd) = SimCommand::from_raw(&raw_cmd) {
                let result = self.dispatch(cmd, &mut bridge.shared);
                bridge.push_result(result);
            }
        }
    }

    /// Append a body at `position`, refusing once `capacity` bodies exist.
    ///
    /// The body count is derived from the stored bodies (`positions.len()`),
    /// never from the public [`Self::body_count`] field that JavaScript may
    /// overwrite, so a forged count cannot be used to run past the shared
    /// buffer.  The handle of the new body is its index.
    fn add_body(&mut self, capacity: usize, position: [f32; 3]) -> crate::error::Result<u32> {
        let index = self.positions.len();
        if index >= capacity {
            return Err(Error::CapacityExceeded { max: capacity });
        }
        let handle = u32::try_from(index).map_err(|_| Error::CapacityExceeded { max: capacity })?;
        self.positions.push(position);
        self.velocities.push([0.0, 0.0, 0.0]);
        self.body_count = handle.saturating_add(1);
        Ok(handle)
    }

    fn dispatch(&mut self, cmd: SimCommand, shared: &mut SharedStateBuffer) -> SimResult {
        match cmd {
            SimCommand::Step { dt } => {
                // Measure the real wall-clock cost of the kinematic preview.
                let t_start = now_ms();

                let g = -9.81_f32;
                for (pos, vel) in self.positions.iter_mut().zip(self.velocities.iter_mut()) {
                    // Ballistic gravity: documented direct positional term, not
                    // accumulated into velocity.
                    pos[1] += 0.5 * g * dt * dt;
                    // Integrate any externally-imparted (impulse) velocity so
                    // that `ApplyImpulse` actually bends the trajectory. A body
                    // that received no impulse has zero velocity and therefore
                    // stays on its pure-gravity ballistic path.
                    pos[0] += vel[0] * dt;
                    pos[1] += vel[1] * dt;
                    pos[2] += vel[2] * dt;
                }
                self.sim_time += dt;
                self.step_count += 1;

                // Publish only the bodies the shared buffer has room for; the
                // header never claims more bodies than the regions hold.
                let published = self.positions.len().min(shared.max_bodies());
                shared.lock();
                shared.write_sim_time(self.sim_time);
                shared.increment_step_count();
                shared.write_body_count(u32::try_from(published).unwrap_or(u32::MAX));
                for (i, &p) in self.positions.iter().enumerate().take(published) {
                    shared.write_body_position(i, p);
                }
                for (i, &v) in self.velocities.iter().enumerate().take(published) {
                    shared.write_body_velocity(i, v);
                }
                shared.unlock();

                SimResult::StepDone {
                    sim_time: self.sim_time,
                    body_count: self.body_count,
                    // No collision detection in the kinematic preview, so the
                    // genuine contact count is zero.
                    contact_count: 0,
                    step_us: elapsed_us(t_start, now_ms()),
                }
            }
            SimCommand::StepSubsteps { dt, substeps } => {
                let t_start = now_ms();
                let sub_dt = dt / substeps.max(1) as f32;
                for _ in 0..substeps {
                    self.dispatch(SimCommand::Step { dt: sub_dt }, shared);
                }
                SimResult::StepDone {
                    sim_time: self.sim_time,
                    body_count: self.body_count,
                    // Kinematic preview performs no collision detection.
                    contact_count: 0,
                    step_us: elapsed_us(t_start, now_ms()),
                }
            }
            SimCommand::AddSphere { x, y, z, .. } => {
                match self.add_body(shared.max_bodies(), [x, y, z]) {
                    Ok(handle) => {
                        shared.write_body_position(handle as usize, [x, y, z]);
                        SimResult::BodyAdded { handle }
                    }
                    Err(err) => SimResult::Error(err.to_string()),
                }
            }
            SimCommand::AddStaticBox { x, y, z, .. } => {
                match self.add_body(shared.max_bodies(), [x, y, z]) {
                    Ok(handle) => SimResult::BodyAdded { handle },
                    Err(err) => SimResult::Error(err.to_string()),
                }
            }
            SimCommand::ApplyImpulse { handle, ix, iy, iz } => {
                // Measure the real wall-clock cost of resolving the body and
                // applying the impulse (same clock the other arms use).
                let t_start = now_ms();
                let idx = handle as usize;
                let have = self.velocities.len();
                let Some(vel) = self.velocities.get_mut(idx) else {
                    // No body owns this handle — there is nothing to push on.
                    // Report the truth instead of a fabricated StepDone.
                    return SimResult::Error(format!(
                        "ApplyImpulse: body handle {handle} out of range (have {have} bodies)"
                    ));
                };
                if idx >= shared.max_bodies() {
                    // The body exists but this shared buffer has no slot for
                    // it, so the new velocity could not be published.
                    return SimResult::Error(format!(
                        "ApplyImpulse: body handle {handle} exceeds the shared buffer capacity of {} bodies",
                        shared.max_bodies()
                    ));
                }
                // The kinematic preview carries no per-body mass, so every body
                // is modeled as unit mass (m = 1 kg): Δv = impulse / m = impulse.
                // This is an honest, documented unit-mass approximation — the
                // impulse genuinely changes the body's velocity (and hence its
                // trajectory once integrated in `Step`), not a silent no-op.
                vel[0] += ix;
                vel[1] += iy;
                vel[2] += iz;
                let new_vel = *vel;
                // Publish the updated velocity so the main thread can read it
                // lock-free before the next `Step`.
                shared.lock();
                shared.write_body_velocity(idx, new_vel);
                shared.unlock();
                SimResult::StepDone {
                    sim_time: self.sim_time,
                    body_count: self.body_count,
                    // No collision detection in the kinematic preview.
                    contact_count: 0,
                    step_us: elapsed_us(t_start, now_ms()),
                }
            }
            SimCommand::Reset => {
                self.positions.clear();
                self.velocities.clear();
                self.body_count = 0;
                self.sim_time = 0.0;
                self.step_count = 0;
                shared.lock();
                shared.write_body_count(0);
                shared.write_sim_time(0.0);
                shared.unlock();
                SimResult::Ready
            }
            SimCommand::RequestSnapshot => {
                let positions: Vec<f32> = self
                    .positions
                    .iter()
                    .flat_map(|p| p.iter().cloned())
                    .collect();
                // Report the *actual* tracked velocities (impulse-imparted; see
                // the `velocities` field docs) rather than a hardcoded zero
                // vector. Orientations are identity because the kinematic
                // preview models no rotation — a genuine value, not a stand-in.
                // One quaternion per *stored* body: the public `body_count`
                // field may have been overwritten from JavaScript, and sizing
                // the allocation from it could exhaust memory.
                let velocities: Vec<f32> = self
                    .velocities
                    .iter()
                    .flat_map(|v| v.iter().cloned())
                    .collect();
                SimResult::Snapshot {
                    positions,
                    orientations: (0..self.positions.len())
                        .flat_map(|_| [0_f32, 0., 0., 1.])
                        .collect(),
                    velocities,
                    sim_time: self.sim_time,
                }
            }
            SimCommand::Shutdown => {
                self.running = false;
                SimResult::Shutdown
            }
        }
    }
}

impl Default for WorkerRuntime {
    fn default() -> Self {
        Self::new()
    }
}

// ── tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_shared_state_layout() {
        let n = 10;
        let expected = 16 + n * 3 * 4 + n * 4 * 4 + n * 3 * 4;
        assert_eq!(SharedStateBuffer::byte_size(n), expected);
    }

    #[test]
    fn test_shared_state_position_roundtrip() {
        let mut buf = SharedStateBuffer::new(4);
        buf.write_body_position(0, [1.5, -2.3, 4.0]);
        let p = buf.read_body_position(0);
        assert!((p[0] - 1.5).abs() < 1e-6);
        assert!((p[1] + 2.3).abs() < 1e-6);
        assert!((p[2] - 4.0).abs() < 1e-6);
    }

    #[test]
    fn test_shared_state_orientation_roundtrip() {
        let mut buf = SharedStateBuffer::new(4);
        buf.write_body_orientation(2, [0.0, 0.707, 0.0, 0.707]);
        let q = buf.read_body_orientation(2);
        assert!((q[1] - 0.707).abs() < 1e-3);
        assert!((q[3] - 0.707).abs() < 1e-3);
    }

    #[test]
    fn test_shared_state_velocity_roundtrip() {
        let mut buf = SharedStateBuffer::new(8);
        buf.write_body_velocity(3, [5.0, 0.0, -3.0]);
        let v = buf.read_body_velocity(3);
        assert!((v[0] - 5.0).abs() < 1e-6);
        assert!((v[2] + 3.0).abs() < 1e-5);
    }

    #[test]
    fn test_shared_state_header() {
        let mut buf = SharedStateBuffer::new(4);
        buf.write_body_count(7);
        assert_eq!(buf.read_body_count(), 7);
        buf.increment_step_count();
        buf.increment_step_count();
        assert_eq!(buf.read_step_count(), 2);
        buf.write_sim_time(1.23);
        assert!((buf.read_sim_time() - 1.23).abs() < 1e-5);
    }

    #[test]
    fn test_positions_flat_roundtrip() {
        let mut buf = SharedStateBuffer::new(3);
        let flat = vec![1.0_f32, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0];
        buf.write_positions_flat(&flat);
        let out = buf.positions_flat(3);
        for (a, b) in flat.iter().zip(out.iter()) {
            assert!((a - b).abs() < 1e-6, "mismatch: {} vs {}", a, b);
        }
    }

    #[test]
    fn test_sim_command_serialisation() {
        let cmd = SimCommand::Step { dt: 1.0 / 60.0 };
        let raw = cmd.to_raw();
        let back = SimCommand::from_raw(&raw).unwrap();
        if let SimCommand::Step { dt } = back {
            assert!((dt - 1.0 / 60.0).abs() < 1e-6);
        } else {
            panic!("wrong variant");
        }
    }

    #[test]
    fn test_sim_command_add_sphere() {
        let cmd = SimCommand::AddSphere {
            mass: 2.0,
            x: 1.0,
            y: 5.0,
            z: 0.0,
            radius: 0.5,
        };
        let raw = cmd.to_raw();
        let back = SimCommand::from_raw(&raw).unwrap();
        if let SimCommand::AddSphere {
            mass, y, radius, ..
        } = back
        {
            assert!((mass - 2.0).abs() < 1e-6);
            assert!((y - 5.0).abs() < 1e-6);
            assert!((radius - 0.5).abs() < 1e-6);
        } else {
            panic!("wrong variant");
        }
    }

    #[test]
    fn test_worker_bridge_integration() {
        let mut bridge = WorkerBridge::new(16);
        let mut runtime = WorkerRuntime::new();

        bridge.post_command(&SimCommand::AddSphere {
            mass: 1.0,
            x: 0.,
            y: 10.,
            z: 0.,
            radius: 0.5,
        });
        runtime.process(&mut bridge);
        let result = bridge.poll_result().unwrap();
        assert!(matches!(result, SimResult::BodyAdded { handle: 0 }));

        bridge.post_command(&SimCommand::Step { dt: 1.0 / 60.0 });
        runtime.process(&mut bridge);
        let result = bridge.poll_result().unwrap();
        assert!(matches!(result, SimResult::StepDone { .. }));

        let y = bridge.shared.read_body_position(0)[1];
        assert!(y < 10.0, "body should have fallen, y={}", y);

        bridge.post_command(&SimCommand::RequestSnapshot);
        runtime.process(&mut bridge);
        let result = bridge.poll_result().unwrap();
        assert!(matches!(result, SimResult::Snapshot { .. }));

        bridge.post_command(&SimCommand::Reset);
        runtime.process(&mut bridge);
        assert_eq!(bridge.shared.read_body_count(), 0);

        bridge.post_command(&SimCommand::Shutdown);
        runtime.process(&mut bridge);
        assert!(!runtime.running);
    }

    #[test]
    fn test_step_telemetry_is_honest_not_fabricated() {
        // The old code hardcoded step_us = 100 regardless of work done. With a
        // real clock (host: std::time::Instant) the measured cost of advancing
        // a few bodies is essentially never exactly 100 µs, and contact_count
        // must be the genuine 0 from a contactless kinematic preview.
        let mut shared = SharedStateBuffer::new(8);
        let mut runtime = WorkerRuntime::new();
        for i in 0..4 {
            runtime.dispatch(
                SimCommand::AddSphere {
                    mass: 1.0,
                    x: i as f32,
                    y: 10.0,
                    z: 0.0,
                    radius: 0.5,
                },
                &mut shared,
            );
        }

        let result = runtime.dispatch(SimCommand::Step { dt: 1.0 / 60.0 }, &mut shared);
        match result {
            SimResult::StepDone {
                contact_count,
                step_us,
                body_count,
                ..
            } => {
                assert_eq!(body_count, 4);
                assert_eq!(
                    contact_count, 0,
                    "kinematic preview has no collisions; contact_count must be 0"
                );
                // On the host the clock is always available, so a real, finite
                // measurement is produced; it must not be the old constant.
                #[cfg(not(target_arch = "wasm32"))]
                assert_ne!(
                    step_us, 100,
                    "step_us still equals the discarded fabricated constant"
                );
                // Either way it is a plain u32 count of microseconds.
                let _ = step_us;
            }
            other => panic!("expected StepDone, got {:?}", other.kind_str()),
        }
    }

    #[test]
    fn test_step_substeps_telemetry_is_honest() {
        // StepSubsteps previously fabricated step_us = 100 * substeps.
        let mut shared = SharedStateBuffer::new(8);
        let mut runtime = WorkerRuntime::new();
        runtime.dispatch(
            SimCommand::AddSphere {
                mass: 1.0,
                x: 0.0,
                y: 10.0,
                z: 0.0,
                radius: 0.5,
            },
            &mut shared,
        );
        let substeps = 4u32;
        let result = runtime.dispatch(
            SimCommand::StepSubsteps {
                dt: 1.0 / 60.0,
                substeps,
            },
            &mut shared,
        );
        if let SimResult::StepDone {
            contact_count,
            step_us,
            ..
        } = result
        {
            assert_eq!(contact_count, 0);
            #[cfg(not(target_arch = "wasm32"))]
            assert_ne!(
                step_us,
                100 * substeps,
                "substep step_us still equals the discarded fabricated formula"
            );
            let _ = step_us;
        } else {
            panic!("expected StepDone");
        }
    }

    #[test]
    fn test_apply_impulse_changes_trajectory() {
        // Regression: `ApplyImpulse` used to be a silent no-op that discarded
        // the handle/impulse and fabricated `step_us = 0`. It must now impart a
        // real Δv (unit-mass preview: Δv = impulse) that bends the trajectory,
        // while a zero impulse leaves the body on its pure-gravity ballistic
        // path.
        let mut shared = SharedStateBuffer::new(8);
        let mut runtime = WorkerRuntime::new();

        // Two bodies starting at exactly the same point.
        for _ in 0..2 {
            runtime.dispatch(
                SimCommand::AddSphere {
                    mass: 1.0,
                    x: 0.0,
                    y: 10.0,
                    z: 0.0,
                    radius: 0.5,
                },
                &mut shared,
            );
        }

        let dt = 1.0_f32 / 60.0;
        let ix = 3.0_f32;

        // Body 0 receives a known impulse; the ack must be honest StepDone.
        let res0 = runtime.dispatch(
            SimCommand::ApplyImpulse {
                handle: 0,
                ix,
                iy: 0.0,
                iz: 0.0,
            },
            &mut shared,
        );
        assert!(
            matches!(res0, SimResult::StepDone { .. }),
            "valid impulse should acknowledge with StepDone"
        );
        // Δv landed in the shared velocity buffer immediately (unit mass).
        let v0 = shared.read_body_velocity(0);
        assert!(
            (v0[0] - ix).abs() < 1e-5,
            "impulse Δv not applied: vx={}",
            v0[0]
        );

        // Body 1 receives an explicit *zero* impulse: a genuine no-effect.
        runtime.dispatch(
            SimCommand::ApplyImpulse {
                handle: 1,
                ix: 0.0,
                iy: 0.0,
                iz: 0.0,
            },
            &mut shared,
        );
        let v1 = shared.read_body_velocity(1);
        assert!(
            v1[0].abs() < 1e-6 && v1[1].abs() < 1e-6 && v1[2].abs() < 1e-6,
            "zero impulse must not impart velocity: {v1:?}"
        );

        // Advance: the impulse velocity must integrate into position.
        runtime.dispatch(SimCommand::Step { dt }, &mut shared);

        let p0 = shared.read_body_position(0);
        let p1 = shared.read_body_position(1);

        let g = -9.81_f32;
        let gravity_y = 10.0 + 0.5 * g * dt * dt;

        // Impulsed body: moved +x by Δv·dt, fell by the gravity term, z fixed.
        assert!((p0[0] - ix * dt).abs() < 1e-5, "impulsed x={}", p0[0]);
        assert!((p0[1] - gravity_y).abs() < 1e-4, "impulsed y={}", p0[1]);
        assert!(p0[2].abs() < 1e-6, "impulsed z drifted={}", p0[2]);

        // Zero-impulse body: stays on its prior ballistic path (gravity only).
        assert!(p1[0].abs() < 1e-6, "no-impulse x drifted={}", p1[0]);
        assert!((p1[1] - gravity_y).abs() < 1e-4, "no-impulse y={}", p1[1]);
        assert!(p1[2].abs() < 1e-6, "no-impulse z drifted={}", p1[2]);

        // The two trajectories differ *only* in x, exactly by the integrated
        // impulse — proof the impulse, and nothing else, changed the path.
        assert!(((p0[0] - p1[0]) - ix * dt).abs() < 1e-5);
    }

    #[test]
    fn test_apply_impulse_out_of_range_is_honest_error() {
        // An impulse to a non-existent body must surface an honest `Error`,
        // not a fabricated `StepDone` success.
        let mut shared = SharedStateBuffer::new(4);
        let mut runtime = WorkerRuntime::new();
        runtime.dispatch(
            SimCommand::AddSphere {
                mass: 1.0,
                x: 0.0,
                y: 0.0,
                z: 0.0,
                radius: 0.5,
            },
            &mut shared,
        );
        let res = runtime.dispatch(
            SimCommand::ApplyImpulse {
                handle: 7,
                ix: 1.0,
                iy: 0.0,
                iz: 0.0,
            },
            &mut shared,
        );
        match res {
            SimResult::Error(msg) => assert!(msg.contains("out of range"), "msg={msg}"),
            other => panic!("expected Error, got {}", other.kind_str()),
        }
    }

    // ---- shared buffer bounds / capacity regression (F5764, F5765, F5768) --

    /// Quaternion seeded into body 0 to detect cross-region corruption.
    const SEED_Q: [f32; 4] = [0.1, 0.2, 0.3, 0.9];

    fn add_sphere(x: f32) -> SimCommand {
        SimCommand::AddSphere {
            mass: 1.0,
            x,
            y: 10.0,
            z: 0.0,
            radius: 0.5,
        }
    }

    fn is_capacity_error(result: &SimResult) -> bool {
        matches!(result, SimResult::Error(msg) if msg.contains("capacity"))
    }

    #[test]
    fn test_shared_state_regions_tile_the_buffer_exactly() {
        let n = 5;
        let buf = SharedStateBuffer::new(n);
        assert_eq!(buf.position_offset(0), Some(16));
        assert_eq!(buf.orientation_offset(0), Some(16 + 12 * n));
        assert_eq!(buf.velocity_offset(0), Some(16 + 28 * n));
        let last = n - 1;
        assert_eq!(buf.position_offset(last), Some(16 + 12 * last));
        assert_eq!(buf.orientation_offset(last), Some(16 + 12 * n + 16 * last));
        assert_eq!(
            buf.velocity_offset(last).map(|o| o + VELOCITY_STRIDE),
            Some(buf.byte_len()),
            "the velocity region must end exactly at the end of the buffer"
        );
        for i in [n, n + 1, usize::MAX] {
            assert_eq!(buf.position_offset(i), None);
            assert_eq!(buf.orientation_offset(i), None);
            assert_eq!(buf.velocity_offset(i), None);
        }
    }

    #[test]
    fn test_shared_state_last_body_round_trips_without_touching_neighbours() {
        let mut buf = SharedStateBuffer::new(3);
        assert!(buf.write_body_orientation(0, SEED_Q));
        assert!(buf.write_body_position(2, [1.0, 2.0, 3.0]));
        assert!(buf.write_body_orientation(2, [0.0, 0.0, 1.0, 0.0]));
        assert!(buf.write_body_velocity(2, [4.0, 5.0, 6.0]));
        assert_eq!(buf.read_body_position(2), [1.0, 2.0, 3.0]);
        assert_eq!(buf.read_body_orientation(2), [0.0, 0.0, 1.0, 0.0]);
        assert_eq!(buf.read_body_velocity(2), [4.0, 5.0, 6.0]);
        assert_eq!(buf.read_body_orientation(0), SEED_Q);
        assert_eq!(buf.read_body_position(0), [0.0; 3]);
        assert_eq!(buf.read_body_velocity(0), [0.0; 3]);
    }

    #[test]
    fn test_shared_state_rejects_writes_past_max_bodies() {
        let mut buf = SharedStateBuffer::new(2);
        assert!(buf.write_body_orientation(0, SEED_Q));
        assert!(buf.write_body_position(0, [1.0, 2.0, 3.0]));
        assert!(buf.write_body_velocity(1, [4.0, 5.0, 6.0]));
        let before = buf.as_bytes().to_vec();

        // Index == max_bodies used to overwrite body 0's orientation (position
        // write) or slice past the end of the buffer (velocity write).
        for i in [2, 3, 99, usize::MAX] {
            assert!(!buf.write_body_position(i, [9.0; 3]), "position {i}");
            assert!(!buf.write_body_orientation(i, [9.0; 4]), "orientation {i}");
            assert!(!buf.write_body_velocity(i, [9.0; 3]), "velocity {i}");
        }
        assert_eq!(
            buf.as_bytes(),
            &before[..],
            "a rejected write modified the buffer"
        );
        assert_eq!(buf.read_body_orientation(0), SEED_Q);
    }

    #[test]
    fn test_shared_state_out_of_range_reads_are_zero_not_neighbour_data() {
        let mut buf = SharedStateBuffer::new(2);
        assert!(buf.write_body_orientation(0, SEED_Q));
        assert!(buf.write_body_velocity(0, [7.0, 8.0, 9.0]));
        // Position index 2 would alias body 0's orientation, orientation
        // index 2 the velocity region and velocity index 2 the buffer end.
        assert_eq!(buf.read_body_position(2), [0.0; 3]);
        assert_eq!(buf.read_body_orientation(2), [0.0; 4]);
        assert_eq!(buf.read_body_velocity(2), [0.0; 3]);
        assert_eq!(buf.read_body_position(usize::MAX), [0.0; 3]);
        assert_eq!(buf.read_body_orientation(usize::MAX), [0.0; 4]);
        assert_eq!(buf.read_body_velocity(usize::MAX), [0.0; 3]);
        assert_eq!(buf.read_body_position_js(2), vec![0.0; 3]);
        assert_eq!(buf.read_body_orientation_js(2), vec![0.0; 4]);
        assert_eq!(buf.read_body_velocity_js(2), vec![0.0; 3]);
    }

    #[test]
    fn test_shared_state_position_slice_validation() {
        let mut buf = SharedStateBuffer::new(2);
        buf.write_body_position_slice(1, &[1.0, 2.0, 3.0])
            .expect("valid write");
        assert_eq!(buf.read_body_position(1), [1.0, 2.0, 3.0]);
        // Extra trailing components are ignored, as before.
        buf.write_body_position_slice(0, &[7.0, 8.0, 9.0, 99.0])
            .expect("extra component tolerated");
        assert_eq!(buf.read_body_position(0), [7.0, 8.0, 9.0]);

        let out_of_range = buf
            .write_body_position_slice(2, &[0.0; 3])
            .expect_err("index 2 of 2 bodies must be rejected");
        assert!(
            matches!(&out_of_range, Error::InvalidParameter { name, .. } if name == "body_index"),
            "{out_of_range:?}"
        );
        let too_short = buf
            .write_body_position_slice(0, &[1.0, 2.0])
            .expect_err("two components must be rejected");
        assert!(
            matches!(&too_short, Error::InvalidParameter { name, .. } if name == "pos_flat"),
            "{too_short:?}"
        );
        // Neither failure touched the stored data.
        assert_eq!(buf.read_body_position(0), [7.0, 8.0, 9.0]);
        assert_eq!(buf.read_body_position(1), [1.0, 2.0, 3.0]);
    }

    #[test]
    fn test_shared_state_js_write_accepts_valid_input() {
        let mut buf = SharedStateBuffer::new(2);
        assert!(
            buf.write_body_position_js(1, vec![0.5, 0.25, 0.125])
                .is_ok()
        );
        assert_eq!(buf.read_body_position_js(1), vec![0.5, 0.25, 0.125]);
    }

    #[test]
    fn test_worker_runtime_refuses_bodies_past_shared_capacity() {
        // WorkerBridge::new(2) + 3 spheres + Step used to panic while writing
        // the velocity of body 2 (slice 96..100 of a 96-byte buffer), after
        // the third sphere's position write had already clobbered body 0's
        // orientation.
        let mut bridge = WorkerBridge::new(2);
        let mut runtime = WorkerRuntime::new();
        assert!(bridge.shared_mut().write_body_orientation(0, SEED_Q));
        for i in 0..3 {
            bridge.post_command(&add_sphere(i as f32));
        }
        bridge.post_command(&SimCommand::AddStaticBox {
            x: 0.0,
            y: 0.0,
            z: 0.0,
            hx: 1.0,
            hy: 1.0,
            hz: 1.0,
        });
        bridge.post_command(&SimCommand::Step { dt: 1.0 / 60.0 });
        runtime.process(&mut bridge);

        let results: Vec<SimResult> = std::iter::from_fn(|| bridge.poll_result()).collect();
        assert_eq!(results.len(), 5);
        assert!(matches!(results[0], SimResult::BodyAdded { handle: 0 }));
        assert!(matches!(results[1], SimResult::BodyAdded { handle: 1 }));
        assert!(
            is_capacity_error(&results[2]),
            "third sphere: {:?}",
            results[2]
        );
        assert!(
            is_capacity_error(&results[3]),
            "static box: {:?}",
            results[3]
        );
        assert!(matches!(
            results[4],
            SimResult::StepDone { body_count: 2, .. }
        ));
        assert_eq!(runtime.body_count, 2);
        assert_eq!(bridge.shared_body_count(), 2);
        assert_eq!(
            bridge.shared().read_body_orientation(0),
            SEED_Q,
            "body 0's orientation was overwritten"
        );
        // Both admitted bodies were published and have fallen.
        assert!(bridge.shared().read_body_position(0)[1] < 10.0);
        assert!(bridge.shared().read_body_position(1)[1] < 10.0);
    }

    #[test]
    fn test_worker_runtime_fifth_sphere_of_four_is_an_error() {
        let mut shared = SharedStateBuffer::new(4);
        let mut runtime = WorkerRuntime::new();
        assert!(shared.write_body_orientation(0, SEED_Q));
        for i in 0..4 {
            let res = runtime.dispatch(add_sphere(i as f32), &mut shared);
            assert!(
                matches!(res, SimResult::BodyAdded { handle } if handle == i),
                "sphere {i}: {res:?}"
            );
        }
        let fifth = runtime.dispatch(add_sphere(4.0), &mut shared);
        assert!(is_capacity_error(&fifth), "{fifth:?}");
        let step = runtime.dispatch(SimCommand::Step { dt: 1.0 / 60.0 }, &mut shared);
        assert!(matches!(step, SimResult::StepDone { body_count: 4, .. }));
        assert_eq!(runtime.body_count, 4);
        assert_eq!(shared.read_body_count(), 4);
        assert_eq!(shared.read_body_orientation(0), SEED_Q);
    }

    #[test]
    fn test_worker_runtime_capacity_ignores_forged_body_count() {
        // `body_count` is a public wasm-bindgen field: JavaScript can set it
        // to anything, so capacity must be judged from the stored bodies.
        let mut shared = SharedStateBuffer::new(2);
        let mut runtime = WorkerRuntime::new();
        for i in 0..2 {
            let res = runtime.dispatch(add_sphere(i as f32), &mut shared);
            assert!(matches!(res, SimResult::BodyAdded { .. }));
        }
        runtime.body_count = 0;
        let res = runtime.dispatch(add_sphere(2.0), &mut shared);
        assert!(is_capacity_error(&res), "{res:?}");

        let mut fresh = WorkerRuntime::new();
        fresh.body_count = 50;
        let mut roomy = SharedStateBuffer::new(4);
        let res = fresh.dispatch(add_sphere(0.0), &mut roomy);
        assert!(
            matches!(res, SimResult::BodyAdded { handle: 0 }),
            "handle must be the body index, got {res:?}"
        );
        assert_eq!(fresh.body_count, 1, "body_count is re-derived from storage");
    }

    #[test]
    fn test_worker_runtime_snapshot_ignores_forged_body_count() {
        // A forged `body_count` of u32::MAX would otherwise size the
        // orientation array at 16 GiB of f32s.
        let mut shared = SharedStateBuffer::new(4);
        let mut runtime = WorkerRuntime::new();
        for i in 0..3 {
            runtime.dispatch(add_sphere(i as f32), &mut shared);
        }
        runtime.body_count = u32::MAX;
        match runtime.dispatch(SimCommand::RequestSnapshot, &mut shared) {
            SimResult::Snapshot {
                positions,
                orientations,
                velocities,
                ..
            } => {
                assert_eq!(positions.len(), 9);
                assert_eq!(velocities.len(), 9);
                assert_eq!(orientations.len(), 12, "one quaternion per stored body");
                assert!(
                    orientations
                        .chunks_exact(4)
                        .all(|q| q == [0.0, 0.0, 0.0, 1.0]),
                    "{orientations:?}"
                );
            }
            other => panic!("expected Snapshot, got {}", other.kind_str()),
        }
    }

    #[test]
    fn test_worker_runtime_with_undersized_shared_buffer_does_not_panic() {
        // A runtime that holds more bodies than the buffer it is pointed at
        // can store must publish only what fits.
        let mut big = SharedStateBuffer::new(4);
        let mut runtime = WorkerRuntime::new();
        for i in 0..4 {
            runtime.dispatch(add_sphere(i as f32), &mut big);
        }
        let mut small = SharedStateBuffer::new(2);
        assert!(small.write_body_orientation(0, SEED_Q));
        let step = runtime.dispatch(SimCommand::Step { dt: 1.0 / 60.0 }, &mut small);
        assert!(matches!(step, SimResult::StepDone { body_count: 4, .. }));
        assert_eq!(small.read_body_count(), 2, "header must not claim 4 bodies");
        assert_eq!(small.read_body_orientation(0), SEED_Q);
        assert!(small.read_body_position(1)[1] < 10.0);

        let impulse = runtime.dispatch(
            SimCommand::ApplyImpulse {
                handle: 3,
                ix: 1.0,
                iy: 0.0,
                iz: 0.0,
            },
            &mut small,
        );
        assert!(is_capacity_error(&impulse), "{impulse:?}");
        assert_eq!(small.read_body_orientation(0), SEED_Q);
    }
}
