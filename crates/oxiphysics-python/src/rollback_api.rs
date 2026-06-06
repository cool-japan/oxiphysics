// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Python bindings for the oxiphysics rollback module.
//!
//! Exposes a byte-command rollback buffer (`RollbackBuffer<Vec<u8>>`) and
//! desync detection to Python.

use oxiphysics::rollback::{DesyncReport, Frame, RollbackBuffer};
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
    /// The hash is FNV-1a 64-bit over the raw snapshot bytes.
    pub fn check_desync_json(&self, tick: u64, remote_hash: u64) -> PyResult<Option<String>> {
        let report = self.inner.get(tick).and_then(|frame| {
            let local_hash = fnv1a_hash(&frame.snapshot);
            if local_hash != remote_hash {
                Some(DesyncReport {
                    tick,
                    local_hash,
                    remote_hash,
                })
            } else {
                None
            }
        });
        match report {
            None => Ok(None),
            Some(r) => serde_json::to_string(&r)
                .map(Some)
                .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string())),
        }
    }

    /// Compute the FNV-1a 64-bit hash of the snapshot at `tick`, or `None`.
    pub fn snapshot_hash(&self, tick: u64) -> Option<u64> {
        self.inner.get(tick).map(|f| fnv1a_hash(&f.snapshot))
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Helper: FNV-1a 64-bit hash (mirrors oxiphysics internals)
// ─────────────────────────────────────────────────────────────────────────────

fn fnv1a_hash(data: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in data {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000003d811_u64);
    }
    h
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

        let local_hash = fnv1a_hash(&snapshot);
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
        let h1 = fnv1a_hash(data);
        let h2 = fnv1a_hash(data);
        assert_eq!(h1, h2);
    }

    #[test]
    fn test_local_peer_id_exported() {
        use oxiphysics::rollback::PeerId;
        // PeerId type alias is u64; LOCAL_PEER_ID == u64::MAX
        let lid: PeerId = oxiphysics::rollback::LOCAL_PEER_ID;
        assert_eq!(lid, u64::MAX);
    }
}
