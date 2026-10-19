// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Python bindings for the oxiphysics rollback module.
//!
//! Exposes a byte-command rollback buffer (`RollbackBuffer<Vec<u8>>`) and
//! desync detection to Python.

use oxiphysics::rollback::{Frame, RollbackBuffer};
use pyo3::prelude::*;
use std::collections::HashMap;

// ─────────────────────────────────────────────────────────────────────────────
// PyRollbackBuffer
// ─────────────────────────────────────────────────────────────────────────────

/// Rolling ring-buffer of recent frames for deterministic rollback/replay.
///
/// Commands are passed as raw bytes (`bytes` or `list[int]` in Python).
/// Use `push_frame` to record a snapshot + inputs for a tick, `get_frame_json`
/// to inspect a stored frame, and `check_desync_json` to detect hash mismatches.
#[pyclass(name = "RollbackBuffer")]
pub struct PyRollbackBuffer {
    inner: RollbackBuffer<Vec<u8>>,
}

#[pymethods]
impl PyRollbackBuffer {
    /// Create a new rollback buffer with `capacity` retained frames.
    #[new]
    pub fn new(capacity: usize) -> Self {
        Self {
            inner: RollbackBuffer::new(capacity),
        }
    }

    /// Push a frame at the given `tick`.
    ///
    /// `snapshot` — raw snapshot bytes.
    /// `inputs`   — mapping from peer id (u64) to command bytes.
    pub fn push_frame(&mut self, tick: u64, snapshot: Vec<u8>, inputs: HashMap<u64, Vec<u8>>) {
        let frame = Frame {
            tick,
            snapshot,
            inputs: inputs.into_iter().collect(),
        };
        self.inner.push(frame);
    }

    /// Return the stored frame at `tick` as JSON, or `None` if not present.
    ///
    /// JSON shape: `{"tick": u64, "snapshot": [u8,...], "inputs": {peer_id: [u8,...]}}`
    pub fn get_frame_json(&self, tick: u64) -> PyResult<Option<String>> {
        match self.inner.get(tick) {
            None => Ok(None),
            Some(frame) => serde_json::to_string(frame)
                .map(Some)
                .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string())),
        }
    }

    /// Number of frames currently in the buffer.
    pub fn len(&self) -> usize {
        self.inner.len()
    }

    /// `true` when no frames are stored.
    pub fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }

    /// The tick of the newest (most recently pushed) frame.
    pub fn head_tick(&self) -> u64 {
        self.inner.head_tick
    }

    /// The tick of the oldest frame, or `None` if the buffer is empty.
    pub fn tail_tick(&self) -> Option<u64> {
        self.inner.tail_tick()
    }

    /// Check for a desync at `tick` by comparing the local snapshot hash
    /// against a remote hash.
    ///
    /// Returns a JSON `DesyncReport` if a mismatch is found, `None` if
    /// hashes match or the tick is not in the buffer.
    ///
    /// The hash is the standard 64-bit FNV-1a (offset basis
    /// `0xcbf29ce484222325`, prime `0x100000001b3`) over the raw snapshot
    /// bytes -- the same function the native `RollbackWorld` and the WASM
    /// bindings use.
    pub fn check_desync_json(&self, tick: u64, remote_hash: u64) -> PyResult<Option<String>> {
        match self.inner.check_desync(tick, remote_hash) {
            None => Ok(None),
            Some(report) => serde_json::to_string(&report)
                .map(Some)
                .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string())),
        }
    }

    /// Compute the 64-bit FNV-1a hash of the snapshot at `tick`, or `None`.
    ///
    /// Identical to `oxiphysics_core::fnv1a_64` over the stored snapshot bytes
    /// and to `RollbackWorld::snapshot_hash` on the native side.
    pub fn snapshot_hash(&self, tick: u64) -> Option<u64> {
        self.inner.snapshot_hash(tick)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Module registration
// ─────────────────────────────────────────────────────────────────────────────

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyRollbackBuffer>()?;
    // Expose LOCAL_PEER_ID constant so Python callers can identify local inputs.
    m.add(
        "ROLLBACK_LOCAL_PEER_ID",
        oxiphysics::rollback::LOCAL_PEER_ID,
    )?;
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use oxiphysics_core::fnv1a_64;

    #[test]
    fn test_rollback_buffer_instantiation() {
        let buf = PyRollbackBuffer::new(8);
        assert!(buf.is_empty());
        assert_eq!(buf.len(), 0);
    }

    #[test]
    fn test_rollback_push_and_get() {
        let mut buf = PyRollbackBuffer::new(8);
        let snapshot = vec![1u8, 2, 3, 4];
        let mut inputs: HashMap<u64, Vec<u8>> = HashMap::new();
        inputs.insert(1, vec![10u8, 20]);
        buf.push_frame(1, snapshot, inputs);
        assert_eq!(buf.len(), 1);
        assert_eq!(buf.head_tick(), 1);
        let json = buf.get_frame_json(1).expect("get_frame_json failed");
        assert!(json.is_some());
        assert!(json.expect("should have json").contains("tick"));
    }

    #[test]
    fn test_rollback_capacity_eviction() {
        let mut buf = PyRollbackBuffer::new(3);
        for i in 0u64..5 {
            buf.push_frame(i, vec![i as u8], HashMap::new());
        }
        // Only last 3 frames retained
        assert_eq!(buf.len(), 3);
        assert!(buf.get_frame_json(0).expect("ok").is_none());
        assert!(buf.get_frame_json(4).expect("ok").is_some());
    }

    #[test]
    fn test_rollback_desync_detection() {
        let mut buf = PyRollbackBuffer::new(8);
        let snapshot = vec![0xABu8, 0xCD];
        buf.push_frame(1, snapshot.clone(), HashMap::new());

        let local_hash = fnv1a_64(&snapshot);
        // Same hash → no desync
        let no_desync = buf.check_desync_json(1, local_hash).expect("check ok");
        assert!(no_desync.is_none());

        // Different hash → desync
        let desync = buf
            .check_desync_json(1, local_hash ^ 0xDEAD)
            .expect("check ok");
        assert!(desync.is_some());
        let json = desync.expect("should have desync report");
        assert!(json.contains("local_hash"));
    }

    #[test]
    fn test_snapshot_hash_is_stable() {
        let data = b"hello rollback";
        let h1 = fnv1a_64(data);
        let h2 = fnv1a_64(data);
        assert_eq!(h1, h2);
    }

    /// Regression: the binding used to multiply by `0x100000003d811`, which is
    /// not the FNV prime, so its hashes never matched a peer following the
    /// documented FNV-1a-64 (such as the WASM bridge).
    #[test]
    fn test_snapshot_hash_is_standard_fnv1a_64() {
        let mut buf = PyRollbackBuffer::new(8);
        buf.push_frame(1, b"a".to_vec(), HashMap::new());
        buf.push_frame(2, Vec::new(), HashMap::new());
        buf.push_frame(3, 5_i64.to_le_bytes().to_vec(), HashMap::new());

        assert_eq!(buf.snapshot_hash(1), Some(0xaf63_dc4c_8601_ec8c));
        assert_eq!(buf.snapshot_hash(2), Some(0xcbf2_9ce4_8422_2325));
        assert_eq!(buf.snapshot_hash(3), Some(0x0de2_1504_f16d_c720));
        assert_eq!(buf.snapshot_hash(4), None);
        // The pre-fix value for b"a".
        assert_ne!(buf.snapshot_hash(1), Some(0x396e_c70f_fdd1_b784));
    }

    #[test]
    fn test_remote_peer_with_standard_hash_sees_no_desync() {
        let mut buf = PyRollbackBuffer::new(8);
        buf.push_frame(7, b"foobar".to_vec(), HashMap::new());

        // A remote peer hashing the same bytes with standard FNV-1a-64.
        let remote_hash = 0x8594_4171_f739_67e8_u64;
        assert!(
            buf.check_desync_json(7, remote_hash)
                .expect("check ok")
                .is_none()
        );

        let report = buf
            .check_desync_json(7, remote_hash ^ 1)
            .expect("check ok")
            .expect("one flipped bit must be reported");
        assert!(report.contains("\"tick\":7"));
        assert!(report.contains(&format!("\"local_hash\":{remote_hash}")));
    }

    #[test]
    fn test_local_peer_id_exported() {
        use oxiphysics::rollback::PeerId;
        // PeerId type alias is u64; LOCAL_PEER_ID == u64::MAX
        let lid: PeerId = oxiphysics::rollback::LOCAL_PEER_ID;
        assert_eq!(lid, u64::MAX);
    }
}
