// Copyright (C) 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! WebAssembly bridge for the rollback / deterministic lockstep system.
//!
//! Exposes a JSON-oriented surface around rollback frames, desync detection,
//! and a tick-based snapshot ring buffer suitable for use across the WASM boundary.

use serde::{Deserialize, Serialize};
use std::collections::VecDeque;
use wasm_bindgen::prelude::*;

use crate::wasm_helpers::{err_to_jsvalue, to_js_value};

// ---------------------------------------------------------------------------
// WasmPeerId
// ---------------------------------------------------------------------------

/// Identifies a network peer. [`LOCAL_PEER_ID`] is used for local inputs.
pub type WasmPeerId = u64;

/// Reserved peer ID for locally submitted inputs.
pub const LOCAL_PEER_ID: WasmPeerId = u64::MAX;

// ---------------------------------------------------------------------------
// WasmFrameRecord — one tick's saved state
// ---------------------------------------------------------------------------

/// A single tick frame: snapshot bytes + local input hash.
///
/// `tick`, `snapshot_hash` (`u64`), `snapshot` (`Vec<u8>`), and
/// `local_input_json` (`String`) are private; JavaScript reads them through
/// accessor methods on the `#[wasm_bindgen]` impl block.
#[wasm_bindgen]
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct WasmFrameRecord {
    /// Tick number — private; JS reads as `f64`.
    pub(crate) tick: u64,
    /// Raw snapshot bytes — private; JS reads via `get_snapshot()`.
    pub(crate) snapshot: Vec<u8>,
    /// FNV-1a hash of the snapshot — private; JS reads as `f64`.
    pub(crate) snapshot_hash: u64,
    /// Local input JSON — private; JS reads via `get_local_input_json()`.
    pub(crate) local_input_json: String,
}

impl WasmFrameRecord {
    /// Compute FNV-1a 64-bit hash of byte slice.
    pub fn hash_bytes(data: &[u8]) -> u64 {
        let mut hash: u64 = 14_695_981_039_346_656_037;
        for &b in data {
            hash ^= b as u64;
            hash = hash.wrapping_mul(1_099_511_628_211);
        }
        hash
    }
}

#[wasm_bindgen]
impl WasmFrameRecord {
    /// Tick number as `f64` (`u64` unsupported by wasm-bindgen).
    #[wasm_bindgen(js_name = "tick")]
    pub fn tick_js(&self) -> f64 {
        self.tick as f64
    }

    /// FNV-1a hash of the snapshot as `f64` (`u64` unsupported).
    ///
    /// Note: only ~53 bits are exactly representable in `f64`; round-trip
    /// through JS may lose precision.  Use `snapshot_hash_lo`/`_hi` for an
    /// exact split.
    #[wasm_bindgen(js_name = "snapshot_hash")]
    pub fn snapshot_hash_js(&self) -> f64 {
        self.snapshot_hash as f64
    }

    /// Low 32 bits of the snapshot hash.
    #[wasm_bindgen(js_name = "snapshot_hash_lo")]
    pub fn snapshot_hash_lo_js(&self) -> u32 {
        (self.snapshot_hash & 0xffff_ffff) as u32
    }

    /// High 32 bits of the snapshot hash.
    #[wasm_bindgen(js_name = "snapshot_hash_hi")]
    pub fn snapshot_hash_hi_js(&self) -> u32 {
        (self.snapshot_hash >> 32) as u32
    }

    /// Raw snapshot bytes as `Uint8Array`.
    #[wasm_bindgen(js_name = "get_snapshot")]
    pub fn get_snapshot_js(&self) -> Vec<u8> {
        self.snapshot.clone()
    }

    /// Local input JSON string.
    #[wasm_bindgen(js_name = "get_local_input_json")]
    pub fn get_local_input_json_js(&self) -> String {
        self.local_input_json.clone()
    }

    /// Compute the FNV-1a 64-bit hash of the supplied byte slice as `f64`.
    #[wasm_bindgen(js_name = "compute_hash")]
    pub fn compute_hash_js(data: &[u8]) -> f64 {
        WasmFrameRecord::hash_bytes(data) as f64
    }
}

// ---------------------------------------------------------------------------
// WasmDesyncReport
// ---------------------------------------------------------------------------

/// Emitted when a hash mismatch is detected between local and remote peers.
///
/// All fields are `u64` and therefore private; JavaScript reads them as
/// `f64` accessors on the `#[wasm_bindgen]` impl block.
#[wasm_bindgen]
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct WasmDesyncReport {
    /// Tick at which the desync was detected — private.
    pub(crate) tick: u64,
    /// Local snapshot hash — private.
    pub(crate) local_hash: u64,
    /// Remote (authoritative) snapshot hash — private.
    pub(crate) remote_hash: u64,
}

#[wasm_bindgen]
impl WasmDesyncReport {
    /// Desync tick as `f64`.
    #[wasm_bindgen(js_name = "tick")]
    pub fn tick_js(&self) -> f64 {
        self.tick as f64
    }

    /// Local snapshot hash as `f64`.
    #[wasm_bindgen(js_name = "local_hash")]
    pub fn local_hash_js(&self) -> f64 {
        self.local_hash as f64
    }

    /// Remote snapshot hash as `f64`.
    #[wasm_bindgen(js_name = "remote_hash")]
    pub fn remote_hash_js(&self) -> f64 {
        self.remote_hash as f64
    }
}

// ---------------------------------------------------------------------------
// WasmRollbackError
// ---------------------------------------------------------------------------

/// Error variants for rollback operations.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum WasmRollbackError {
    /// Requested tick is not in the buffer.
    TickNotInBuffer(u64),
    /// Buffer capacity must be > 0.
    ZeroCapacity,
}

impl std::fmt::Display for WasmRollbackError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            WasmRollbackError::TickNotInBuffer(t) => {
                write!(f, "tick {t} is not in the rollback buffer")
            }
            WasmRollbackError::ZeroCapacity => write!(f, "buffer capacity must be > 0"),
        }
    }
}

// ---------------------------------------------------------------------------
// WasmRollbackBuffer — ring buffer of frames
// ---------------------------------------------------------------------------

/// Ring buffer of recent [`WasmFrameRecord`]s.
#[wasm_bindgen]
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct WasmRollbackBuffer {
    frames: VecDeque<WasmFrameRecord>,
    capacity: usize,
}

impl WasmRollbackBuffer {
    /// Create a new buffer with the given `capacity` (number of ticks retained).
    pub fn new(capacity: usize) -> Result<Self, WasmRollbackError> {
        if capacity == 0 {
            return Err(WasmRollbackError::ZeroCapacity);
        }
        Ok(WasmRollbackBuffer {
            frames: VecDeque::with_capacity(capacity),
            capacity,
        })
    }

    /// Push a frame, evicting the oldest if at capacity.
    pub fn push(&mut self, frame: WasmFrameRecord) {
        if self.frames.len() == self.capacity {
            self.frames.pop_front();
        }
        self.frames.push_back(frame);
    }

    /// Retrieve a frame by tick, or `None`.
    pub fn get(&self, tick: u64) -> Option<&WasmFrameRecord> {
        self.frames.iter().find(|f| f.tick == tick)
    }

    /// Number of frames currently stored.
    pub fn len(&self) -> usize {
        self.frames.len()
    }

    /// Returns `true` if the buffer is empty.
    pub fn is_empty(&self) -> bool {
        self.frames.is_empty()
    }

    /// Tick of the oldest retained frame, or `None`.
    pub fn oldest_tick(&self) -> Option<u64> {
        self.frames.front().map(|f| f.tick)
    }

    /// Tick of the newest retained frame, or `None`.
    pub fn newest_tick(&self) -> Option<u64> {
        self.frames.back().map(|f| f.tick)
    }

    /// Serialise the buffer as a JSON string.
    pub fn to_json(&self) -> String {
        serde_json::to_string(self).unwrap_or_default()
    }
}

#[wasm_bindgen]
impl WasmRollbackBuffer {
    /// Create a new buffer with the given `capacity` (JS constructor).
    ///
    /// # Errors
    ///
    /// Returns a `JsValue` error string when `capacity == 0`.
    #[wasm_bindgen(constructor)]
    pub fn new_js(capacity: u32) -> std::result::Result<WasmRollbackBuffer, JsValue> {
        WasmRollbackBuffer::new(capacity as usize).map_err(err_to_jsvalue)
    }

    /// Push a pre-built `WasmFrameRecord`.
    #[wasm_bindgen(js_name = "push")]
    pub fn push_js(&mut self, frame: WasmFrameRecord) {
        self.push(frame);
    }

    /// Retrieve a frame by tick number, returning a fresh `WasmFrameRecord`
    /// or `None` when the tick is not in the buffer.
    ///
    /// `tick` is `f64` because `u64` is unsupported by wasm-bindgen.
    #[wasm_bindgen(js_name = "get")]
    pub fn get_js(&self, tick: f64) -> Option<WasmFrameRecord> {
        self.get(tick as u64).cloned()
    }

    /// Number of frames currently stored.
    #[wasm_bindgen(js_name = "len")]
    pub fn len_js(&self) -> u32 {
        self.len() as u32
    }

    /// `true` when the buffer is empty.
    #[wasm_bindgen(js_name = "is_empty")]
    pub fn is_empty_js(&self) -> bool {
        self.is_empty()
    }

    /// Maximum number of retained frames.
    #[wasm_bindgen(js_name = "capacity")]
    pub fn capacity_js(&self) -> u32 {
        self.capacity as u32
    }

    /// Tick of the oldest retained frame, or `None`.
    #[wasm_bindgen(js_name = "oldest_tick")]
    pub fn oldest_tick_js(&self) -> Option<f64> {
        self.oldest_tick().map(|t| t as f64)
    }

    /// Tick of the newest retained frame, or `None`.
    #[wasm_bindgen(js_name = "newest_tick")]
    pub fn newest_tick_js(&self) -> Option<f64> {
        self.newest_tick().map(|t| t as f64)
    }

    /// Serialise the buffer as a JSON string.
    #[wasm_bindgen(js_name = "to_json")]
    pub fn to_json_js(&self) -> String {
        self.to_json()
    }
}

// ---------------------------------------------------------------------------
// WasmRollbackSession
// ---------------------------------------------------------------------------

/// WASM rollback session that manages a tick counter, frame buffer,
/// and simple desync detection.
///
/// Fields are private; JavaScript reads them through accessor methods on
/// the `#[wasm_bindgen]` impl block.
#[wasm_bindgen]
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct WasmRollbackSession {
    /// Monotonically-increasing tick counter — private; JS reads as `f64`.
    pub(crate) current_tick: u64,
    /// Underlying frame ring buffer — private; JS clones via `get_buffer()`.
    pub(crate) buffer: WasmRollbackBuffer,
}

impl WasmRollbackSession {
    /// Create a new session with the given buffer `capacity`.
    pub fn new(capacity: usize) -> Result<Self, WasmRollbackError> {
        Ok(WasmRollbackSession {
            current_tick: 0,
            buffer: WasmRollbackBuffer::new(capacity)?,
        })
    }

    /// Record a new frame from the given snapshot bytes and input JSON.
    ///
    /// Advances `current_tick` by 1.
    pub fn record_frame(&mut self, snapshot: Vec<u8>, local_input_json: String) {
        let hash = WasmFrameRecord::hash_bytes(&snapshot);
        let frame = WasmFrameRecord {
            tick: self.current_tick,
            snapshot,
            snapshot_hash: hash,
            local_input_json,
        };
        self.buffer.push(frame);
        self.current_tick += 1;
    }

    /// Check for a desync at `tick` against the provided `remote_hash`.
    ///
    /// Returns `Some(WasmDesyncReport)` if a mismatch is detected,
    /// `None` if the tick is not in the buffer or hashes match.
    pub fn check_desync(&self, tick: u64, remote_hash: u64) -> Option<WasmDesyncReport> {
        let frame = self.buffer.get(tick)?;
        if frame.snapshot_hash != remote_hash {
            Some(WasmDesyncReport {
                tick,
                local_hash: frame.snapshot_hash,
                remote_hash,
            })
        } else {
            None
        }
    }

    /// List ticks currently retained in the buffer.
    pub fn retained_ticks(&self) -> Vec<u64> {
        self.buffer.frames.iter().map(|f| f.tick).collect()
    }

    /// Serialise the session as a JSON string.
    pub fn to_json(&self) -> String {
        serde_json::to_string(self).unwrap_or_default()
    }
}

// ---------------------------------------------------------------------------
// wasm-bindgen JavaScript API
// ---------------------------------------------------------------------------

#[wasm_bindgen]
impl WasmRollbackSession {
    /// Create a new session (JS constructor).
    ///
    /// # Errors
    ///
    /// Returns a `JsValue` error string when `capacity == 0`.
    #[wasm_bindgen(constructor)]
    pub fn new_js(capacity: u32) -> std::result::Result<WasmRollbackSession, JsValue> {
        WasmRollbackSession::new(capacity as usize).map_err(err_to_jsvalue)
    }

    /// Current tick as `f64` (`u64` unsupported by wasm-bindgen).
    #[wasm_bindgen(js_name = "current_tick")]
    pub fn current_tick_js(&self) -> f64 {
        self.current_tick as f64
    }

    /// Snapshot of the underlying frame ring buffer.
    #[wasm_bindgen(js_name = "get_buffer")]
    pub fn get_buffer_js(&self) -> WasmRollbackBuffer {
        self.buffer.clone()
    }

    /// Record a new frame and advance `current_tick`.
    #[wasm_bindgen(js_name = "record_frame")]
    pub fn record_frame_js(&mut self, snapshot: Vec<u8>, local_input_json: String) {
        self.record_frame(snapshot, local_input_json);
    }

    /// Look up a frame by tick (`f64` because `u64` is unsupported).
    #[wasm_bindgen(js_name = "get_frame")]
    pub fn get_frame_js(&self, tick: f64) -> Option<WasmFrameRecord> {
        self.buffer.get(tick as u64).cloned()
    }

    /// Detect a desync at the given tick by comparing against `remote_hash`.
    ///
    /// `tick` and `remote_hash` are `f64` (the latter is exact only for
    /// values that fit in 53 bits).
    #[wasm_bindgen(js_name = "check_desync")]
    pub fn check_desync_js(&self, tick: f64, remote_hash: f64) -> Option<WasmDesyncReport> {
        self.check_desync(tick as u64, remote_hash as u64)
    }

    /// Detect a desync by passing the remote hash as low/high `u32` halves
    /// for full 64-bit precision.
    #[wasm_bindgen(js_name = "check_desync_u64_parts")]
    pub fn check_desync_u64_parts_js(
        &self,
        tick: f64,
        remote_hash_lo: u32,
        remote_hash_hi: u32,
    ) -> Option<WasmDesyncReport> {
        let remote_hash = (u64::from(remote_hash_hi) << 32) | u64::from(remote_hash_lo);
        self.check_desync(tick as u64, remote_hash)
    }

    /// Ticks currently retained in the buffer (as `Float64Array`).
    #[wasm_bindgen(js_name = "retained_ticks")]
    pub fn retained_ticks_js(&self) -> Vec<f64> {
        self.retained_ticks()
            .into_iter()
            .map(|t| t as f64)
            .collect()
    }

    /// Serialise the session as a `JsValue` (plain JS object).
    ///
    /// # Errors
    ///
    /// Returns a `JsValue` error string when serialisation fails.
    #[wasm_bindgen(js_name = "to_js_value")]
    pub fn to_js_value_js(&self) -> std::result::Result<JsValue, JsValue> {
        to_js_value(self)
    }

    /// Serialise the session as a JSON string.
    #[wasm_bindgen(js_name = "to_json")]
    pub fn to_json_js(&self) -> String {
        self.to_json()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rollback_bridge_instantiation() {
        let session = WasmRollbackSession::new(16).expect("valid capacity");
        let json = session.to_json();
        assert!(!json.is_empty(), "to_json should return non-empty string");
    }

    #[test]
    fn test_rollback_bridge_record_and_retrieve() {
        let mut session = WasmRollbackSession::new(8).expect("valid capacity");
        let snap: Vec<u8> = 42_i64.to_le_bytes().to_vec();
        session.record_frame(snap.clone(), "\"move\"".to_string());
        assert_eq!(session.current_tick, 1);
        let frame = session.buffer.get(0).expect("tick 0 should exist");
        assert_eq!(frame.snapshot, snap);
        assert_eq!(frame.tick, 0);
    }

    #[test]
    fn test_rollback_bridge_desync_detection() {
        let mut session = WasmRollbackSession::new(8).expect("valid capacity");
        let snap: Vec<u8> = 5_i64.to_le_bytes().to_vec();
        session.record_frame(snap, "{}".to_string());
        let true_snap: Vec<u8> = 99_i64.to_le_bytes().to_vec();
        let true_hash = WasmFrameRecord::hash_bytes(&true_snap);
        let report = session.check_desync(0, true_hash);
        assert!(
            report.is_some(),
            "different snapshots should trigger desync"
        );
        let r = report.expect("should have report");
        assert_eq!(r.tick, 0);
        assert_ne!(r.local_hash, r.remote_hash);
    }

    #[test]
    fn test_rollback_bridge_no_desync() {
        let mut session = WasmRollbackSession::new(8).expect("valid capacity");
        let snap: Vec<u8> = 42_i64.to_le_bytes().to_vec();
        let hash = WasmFrameRecord::hash_bytes(&snap);
        session.record_frame(snap, "{}".to_string());
        assert!(
            session.check_desync(0, hash).is_none(),
            "matching hash should not trigger desync"
        );
    }

    #[test]
    fn test_rollback_bridge_buffer_capacity() {
        let mut session = WasmRollbackSession::new(3).expect("valid capacity");
        for i in 0..10_u8 {
            let snap = vec![i];
            session.record_frame(snap, "{}".to_string());
        }
        assert_eq!(
            session.buffer.len(),
            3,
            "buffer should hold at most 3 frames"
        );
        let oldest = session.buffer.oldest_tick().expect("has frames");
        assert!(oldest > 0, "oldest tick should be > 0 after 10 frames");
    }

    #[test]
    fn test_rollback_bridge_zero_capacity_error() {
        assert_eq!(
            WasmRollbackSession::new(0),
            Err(WasmRollbackError::ZeroCapacity)
        );
    }
}
