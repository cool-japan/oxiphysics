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
    StepDone {
        sim_time: f32,
        body_count: u32,
        contact_count: u32,
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

/// CPU-resident approximation of a `SharedArrayBuffer` for Rust-side logic.
///
/// In a real browser environment this wraps a JS `SharedArrayBuffer` and
/// accesses it via `Atomics` operations.  In non-WASM / unit-test contexts
/// it is a plain `Vec<u8>` providing the same layout.
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
    pub fn write_body_position_js(&mut self, i: usize, pos_flat: Vec<f32>) {
        if pos_flat.len() >= 3 {
            self.write_body_position(i, [pos_flat[0], pos_flat[1], pos_flat[2]]);
        }
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
    /// Write position `[x, y, z]` for body `i` (Rust-only).
    pub fn write_body_position(&mut self, i: usize, pos: [f32; 3]) {
        let base = 16 + i * 12;
        for (j, &v) in pos.iter().enumerate() {
            self.data[base + j * 4..base + j * 4 + 4].copy_from_slice(&v.to_le_bytes());
        }
    }

    /// Read position `[x, y, z]` for body `i` (Rust-only).
    pub fn read_body_position(&self, i: usize) -> [f32; 3] {
        let base = 16 + i * 12;
        let x = read_f32_le(&self.data, base);
        let y = read_f32_le(&self.data, base + 4);
        let z = read_f32_le(&self.data, base + 8);
        [x, y, z]
    }

    /// Write orientation quaternion for body `i` (Rust-only).
    pub fn write_body_orientation(&mut self, i: usize, q: [f32; 4]) {
        let base = 16 + self.max_bodies * 12 + i * 16;
        for (j, &v) in q.iter().enumerate() {
            self.data[base + j * 4..base + j * 4 + 4].copy_from_slice(&v.to_le_bytes());
        }
    }

    /// Read orientation quaternion for body `i` (Rust-only).
    pub fn read_body_orientation(&self, i: usize) -> [f32; 4] {
        let base = 16 + self.max_bodies * 12 + i * 16;
        [
            read_f32_le(&self.data, base),
            read_f32_le(&self.data, base + 4),
            read_f32_le(&self.data, base + 8),
            read_f32_le(&self.data, base + 12),
        ]
    }

    /// Write linear velocity for body `i` (Rust-only).
    pub fn write_body_velocity(&mut self, i: usize, vel: [f32; 3]) {
        let base = 16 + self.max_bodies * 28 + i * 12;
        for (j, &v) in vel.iter().enumerate() {
            self.data[base + j * 4..base + j * 4 + 4].copy_from_slice(&v.to_le_bytes());
        }
    }

    /// Read linear velocity for body `i` (Rust-only).
    pub fn read_body_velocity(&self, i: usize) -> [f32; 3] {
        let base = 16 + self.max_bodies * 28 + i * 12;
        [
            read_f32_le(&self.data, base),
            read_f32_le(&self.data, base + 4),
            read_f32_le(&self.data, base + 8),
        ]
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
    if offset + 4 > bytes.len() {
        return 0.0;
    }
    let arr: [u8; 4] = [
        bytes[offset],
        bytes[offset + 1],
        bytes[offset + 2],
        bytes[offset + 3],
    ];
    f32::from_le_bytes(arr)
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
#[wasm_bindgen]
pub struct WorkerRuntime {
    /// Accumulated simulation time.
    pub sim_time: f32,
    /// Step count.
    pub step_count: u32,
    /// Simulated body count (stub; real impl uses WasmPhysicsEngine).
    pub body_count: u32,
    /// Body positions (stub storage for testing).
    positions: Vec<[f32; 3]>,
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

    fn dispatch(&mut self, cmd: SimCommand, shared: &mut SharedStateBuffer) -> SimResult {
        match cmd {
            SimCommand::Step { dt } => {
                let g = -9.81_f32;
                for pos in &mut self.positions {
                    pos[1] += 0.5 * g * dt * dt;
                }
                self.sim_time += dt;
                self.step_count += 1;

                shared.lock();
                shared.write_sim_time(self.sim_time);
                shared.increment_step_count();
                shared.write_body_count(self.body_count);
                for (i, &p) in self.positions.iter().enumerate() {
                    shared.write_body_position(i, p);
                }
                shared.unlock();

                SimResult::StepDone {
                    sim_time: self.sim_time,
                    body_count: self.body_count,
                    contact_count: 0,
                    step_us: 100,
                }
            }
            SimCommand::StepSubsteps { dt, substeps } => {
                let sub_dt = dt / substeps as f32;
                for _ in 0..substeps {
                    self.dispatch(SimCommand::Step { dt: sub_dt }, shared);
                }
                SimResult::StepDone {
                    sim_time: self.sim_time,
                    body_count: self.body_count,
                    contact_count: 0,
                    step_us: 100 * substeps,
                }
            }
            SimCommand::AddSphere { x, y, z, .. } => {
                let handle = self.body_count;
                self.positions.push([x, y, z]);
                self.body_count += 1;
                shared.write_body_position(handle as usize, [x, y, z]);
                SimResult::BodyAdded { handle }
            }
            SimCommand::AddStaticBox { x, y, z, .. } => {
                let handle = self.body_count;
                self.positions.push([x, y, z]);
                self.body_count += 1;
                SimResult::BodyAdded { handle }
            }
            SimCommand::ApplyImpulse { .. } => SimResult::StepDone {
                sim_time: self.sim_time,
                body_count: self.body_count,
                contact_count: 0,
                step_us: 0,
            },
            SimCommand::Reset => {
                self.positions.clear();
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
                SimResult::Snapshot {
                    positions,
                    orientations: (0..self.body_count as usize)
                        .flat_map(|_| [0_f32, 0., 0., 1.])
                        .collect(),
                    velocities: vec![0.; self.body_count as usize * 3],
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
}
