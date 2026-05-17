// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Snapshot-based rollback and deterministic lockstep input buffer.
//!
//! This module provides client-prediction / authoritative-server rollback for
//! game and real-time physics simulations. The key idea:
//!
//! - Every tick, the simulation captures a snapshot.
//! - When late input arrives (or a correction from the server), we resimulate
//!   from the confirmed tick forward.
//! - Desync is detected by comparing snapshot hashes between peers.
//!
//! **The module is NOT coupled to any specific physics engine.** The user
//! provides the world type, a command type, and three function pointers:
//!
//! | Function | Signature | Purpose |
//! |----------|-----------|---------|
//! | `advance_fn` | `fn(&mut World, f64, &[Cmd])` | Step simulation one tick |
//! | `snapshot_fn` | `fn(&World) -> Vec<u8>` | Serialise world state |
//! | `restore_fn` | `fn(&mut World, &[u8])` | Restore world from bytes |
//!
//! ## Types
//!
//! - `Frame` — one tick's saved state (snapshot + inputs).
//! - `RollbackBuffer` — ring buffer of recent frames.
//! - `RollbackWorld` — orchestrates stepping, rollback, and desync detection.
//! - `DesyncReport` — emitted when a hash mismatch is detected.
//! - `RollbackError` — error variants for resimulation failures.
//!
//! ## Example
//!
//! ```rust,no_run
//! use oxiphysics::rollback::RollbackWorld;
//!
//! struct Sim { state: i64 }
//!
//! fn advance(w: &mut Sim, _dt: f64, cmds: &[i32]) {
//!     w.state += cmds.iter().map(|c| *c as i64).sum::<i64>();
//! }
//! fn snapshot(w: &Sim) -> Vec<u8> { w.state.to_le_bytes().to_vec() }
//! fn restore(w: &mut Sim, snap: &[u8]) {
//!     if snap.len() == 8 {
//!         let arr: [u8; 8] = snap.try_into().unwrap_or([0; 8]);
//!         w.state = i64::from_le_bytes(arr);
//!     }
//! }
//!
//! let mut rw = RollbackWorld::new(Sim { state: 0 }, 64, advance, snapshot, restore);
//! rw.step(1.0 / 60.0, &[5_i32]);
//! assert_eq!(rw.world.state, 5);
//! ```

#![allow(missing_docs)]
#![allow(dead_code)]

use serde::de::DeserializeOwned;
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, VecDeque};

// ---------------------------------------------------------------------------
// PeerId
// ---------------------------------------------------------------------------

/// Identifies a network peer. `LOCAL_PEER_ID` is used for local inputs.
pub type PeerId = u64;

/// Reserved peer ID for local commands submitted via [`RollbackWorld::step`].
pub const LOCAL_PEER_ID: PeerId = u64::MAX;

// ---------------------------------------------------------------------------
// Frame
// ---------------------------------------------------------------------------

/// One tick's saved state: snapshot bytes plus all inputs for that tick.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(bound(serialize = "Cmd: Serialize", deserialize = "Cmd: DeserializeOwned"))]
pub struct Frame<Cmd>
where
    Cmd: Serialize + DeserializeOwned,
{
    /// Tick number this frame corresponds to.
    pub tick: u64,
    /// Raw snapshot bytes from the user-provided snapshot function.
    pub snapshot: Vec<u8>,
    /// Inputs keyed by peer. Local inputs are stored under [`LOCAL_PEER_ID`].
    pub inputs: HashMap<PeerId, Cmd>,
}

// ---------------------------------------------------------------------------
// RollbackBuffer
// ---------------------------------------------------------------------------

/// Rolling ring-buffer of recent [`Frame`]s.
///
/// Once the buffer reaches `capacity`, the oldest frame is dropped whenever a
/// new frame is pushed.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(bound(serialize = "Cmd: Serialize", deserialize = "Cmd: DeserializeOwned"))]
pub struct RollbackBuffer<Cmd>
where
    Cmd: Serialize + DeserializeOwned,
{
    frames: VecDeque<Frame<Cmd>>,
    /// The tick number of the newest frame (the "head").
    pub head_tick: u64,
    /// Maximum number of frames retained.
    pub capacity: usize,
}

impl<Cmd: Clone + Serialize + DeserializeOwned> RollbackBuffer<Cmd> {
    /// Create a new, empty buffer with the given capacity.
    pub fn new(capacity: usize) -> Self {
        Self {
            frames: VecDeque::with_capacity(capacity),
            head_tick: 0,
            capacity,
        }
    }

    /// Push a new frame. Drops the oldest frame if the buffer is at capacity.
    pub fn push(&mut self, frame: Frame<Cmd>) {
        if self.frames.len() >= self.capacity {
            self.frames.pop_front();
        }
        self.head_tick = frame.tick;
        self.frames.push_back(frame);
    }

    /// Return a reference to the frame at `tick`, or `None` if not in buffer.
    pub fn get(&self, tick: u64) -> Option<&Frame<Cmd>> {
        self.frames.iter().find(|f| f.tick == tick)
    }

    /// Return a mutable reference to the frame at `tick`, or `None`.
    pub fn get_mut(&mut self, tick: u64) -> Option<&mut Frame<Cmd>> {
        self.frames.iter_mut().find(|f| f.tick == tick)
    }

    /// The tick of the oldest frame, or `None` if the buffer is empty.
    pub fn tail_tick(&self) -> Option<u64> {
        self.frames.front().map(|f| f.tick)
    }

    /// Number of frames currently stored.
    pub fn len(&self) -> usize {
        self.frames.len()
    }

    /// `true` when no frames are stored.
    pub fn is_empty(&self) -> bool {
        self.frames.is_empty()
    }
}

// ---------------------------------------------------------------------------
// DesyncReport
// ---------------------------------------------------------------------------

/// Emitted by [`RollbackWorld::check_desync`] when a hash mismatch is found.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct DesyncReport {
    /// The tick at which the desync was detected.
    pub tick: u64,
    /// FNV-1a hash of the local snapshot at that tick.
    pub local_hash: u64,
    /// The remote hash that was compared against.
    pub remote_hash: u64,
}

// ---------------------------------------------------------------------------
// RollbackError
// ---------------------------------------------------------------------------

/// Error variants returned by [`RollbackWorld::resimulate_from`].
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum RollbackError {
    /// The requested tick is not present in the rollback buffer.
    TickNotInBuffer(u64),
    /// The snapshot for the tick immediately before `from_tick` is unavailable.
    PriorSnapshotUnavailable,
}

// ---------------------------------------------------------------------------
// Snapshot hashing (FNV-1a 64-bit)
// ---------------------------------------------------------------------------

/// Compute an FNV-1a 64-bit hash of a snapshot byte slice.
fn hash_snapshot(snap: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in snap {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000003d811_u64);
    }
    h
}

// ---------------------------------------------------------------------------
// RollbackWorld
// ---------------------------------------------------------------------------

/// The full rollback world — generic over the simulation world type `World`
/// and the command type `Cmd`.
///
/// `RollbackWorld` is intentionally **not** `Serialize`/`Deserialize` because
/// it holds function pointers, which are not serialisable.
pub struct RollbackWorld<World, Cmd>
where
    Cmd: Clone + PartialEq + Serialize + DeserializeOwned,
{
    /// The live simulation world.
    pub world: World,
    /// Rolling ring-buffer of past frames.
    pub buffer: RollbackBuffer<Cmd>,
    /// User-provided simulation advance function.
    advance_fn: fn(&mut World, f64, &[Cmd]),
    /// User-provided snapshot function (world → bytes).
    snapshot_fn: fn(&World) -> Vec<u8>,
    /// User-provided restore function (bytes → world).
    restore_fn: fn(&mut World, &[u8]),
    /// The current tick (matches the newest frame in the buffer).
    pub current_tick: u64,
    /// Buffered future inputs: tick → (peer → cmd).
    pending_inputs: HashMap<u64, HashMap<PeerId, Cmd>>,
}

impl<World, Cmd> RollbackWorld<World, Cmd>
where
    Cmd: Clone + PartialEq + Serialize + DeserializeOwned,
{
    /// Create a new `RollbackWorld`.
    ///
    /// An initial frame at tick 0 is captured immediately so that
    /// [`resimulate_from`](Self::resimulate_from) can always restore the state
    /// prior to tick 1.
    pub fn new(
        world: World,
        capacity: usize,
        advance_fn: fn(&mut World, f64, &[Cmd]),
        snapshot_fn: fn(&World) -> Vec<u8>,
        restore_fn: fn(&mut World, &[u8]),
    ) -> Self {
        let initial_snapshot = snapshot_fn(&world);
        let mut buffer = RollbackBuffer::new(capacity);
        buffer.push(Frame {
            tick: 0,
            snapshot: initial_snapshot,
            inputs: HashMap::new(),
        });

        Self {
            world,
            buffer,
            advance_fn,
            snapshot_fn,
            restore_fn,
            current_tick: 0,
            pending_inputs: HashMap::new(),
        }
    }

    /// Record an input from `peer` for `tick`.
    ///
    /// - If `tick` is already in the buffer, the frame's input map is updated.
    /// - If `tick` is in the future, the input is stored in a pending map and
    ///   merged into the frame when that tick is stepped.
    pub fn record_input(&mut self, peer: PeerId, tick: u64, cmd: Cmd) {
        if tick <= self.current_tick {
            // Update an existing frame in the buffer.
            if let Some(frame) = self.buffer.get_mut(tick) {
                frame.inputs.insert(peer, cmd);
            }
        } else {
            // Buffer for a future tick.
            self.pending_inputs
                .entry(tick)
                .or_default()
                .insert(peer, cmd);
        }
    }

    /// Advance by one tick.
    ///
    /// Steps:
    /// 1. Collect all inputs for `current_tick + 1` (local cmds + any pending
    ///    remote inputs).
    /// 2. Call `advance_fn` with a deterministic (PeerId-sorted) input slice.
    /// 3. Capture a snapshot.
    /// 4. Push the frame to the buffer.
    /// 5. Increment `current_tick`.
    pub fn step(&mut self, dt: f64, local_cmds: &[Cmd]) {
        let next_tick = self.current_tick + 1;

        // Start with any buffered future inputs for this tick.
        let mut inputs: HashMap<PeerId, Cmd> =
            self.pending_inputs.remove(&next_tick).unwrap_or_default();

        // Merge local commands under LOCAL_PEER_ID.
        // Only the first local command is stored (single-cmd API).
        if let Some(cmd) = local_cmds.first() {
            inputs.insert(LOCAL_PEER_ID, cmd.clone());
        }

        // Build a deterministic ordered slice (sort by PeerId).
        let ordered_cmds = Self::inputs_to_sorted_slice(&inputs);

        // Advance the world.
        (self.advance_fn)(&mut self.world, dt, &ordered_cmds);

        // Snapshot after advancing.
        let snapshot = (self.snapshot_fn)(&self.world);

        // Store the inputs so resimulate can replay them.
        let frame = Frame {
            tick: next_tick,
            snapshot,
            inputs,
        };
        self.buffer.push(frame);
        self.current_tick = next_tick;
    }

    /// Resimulate from `from_tick` (inclusive) to `current_tick` (inclusive)
    /// using the stored inputs in the buffer.
    ///
    /// # Process
    /// 1. Restore the snapshot at `from_tick - 1` (the state BEFORE `from_tick`).
    /// 2. For each tick `t` in `from_tick ..= current_tick`:
    ///    a. Collect stored inputs (sorted by PeerId for determinism).
    ///    b. Call `advance_fn`.
    ///    c. Update `frame.snapshot` with the new snapshot.
    ///
    /// # Errors
    /// - [`RollbackError::TickNotInBuffer`] if `from_tick` is not in the buffer.
    /// - [`RollbackError::PriorSnapshotUnavailable`] if `from_tick == 0` (no
    ///   predecessor) or if the frame at `from_tick - 1` is missing.
    pub fn resimulate_from(&mut self, from_tick: u64, dt: f64) -> Result<(), RollbackError> {
        // Guard: from_tick must be in the buffer.
        if self.buffer.get(from_tick).is_none() {
            return Err(RollbackError::TickNotInBuffer(from_tick));
        }

        // Restore prior snapshot.
        if from_tick == 0 {
            return Err(RollbackError::PriorSnapshotUnavailable);
        }
        let prior_tick = from_tick - 1;
        let prior_snapshot = self
            .buffer
            .get(prior_tick)
            .map(|f| f.snapshot.clone())
            .ok_or(RollbackError::PriorSnapshotUnavailable)?;

        (self.restore_fn)(&mut self.world, &prior_snapshot);

        // Re-advance each tick and refresh its snapshot.
        let end_tick = self.current_tick;
        for t in from_tick..=end_tick {
            // Collect inputs for tick t, sorted by PeerId.
            let ordered_cmds: Vec<Cmd> = self
                .buffer
                .get(t)
                .map(|f| Self::inputs_to_sorted_slice(&f.inputs))
                .unwrap_or_default();

            (self.advance_fn)(&mut self.world, dt, &ordered_cmds);

            // Update the stored snapshot.
            let new_snapshot = (self.snapshot_fn)(&self.world);
            if let Some(frame) = self.buffer.get_mut(t) {
                frame.snapshot = new_snapshot;
            }
        }

        Ok(())
    }

    /// Compare the local snapshot hash at `tick` to `remote_hash`.
    ///
    /// Returns `Some(DesyncReport)` on a mismatch, or `None` if the hashes
    /// match or the tick is not in the buffer.
    pub fn check_desync(&self, tick: u64, remote_hash: u64) -> Option<DesyncReport> {
        let frame = self.buffer.get(tick)?;
        let local_hash = hash_snapshot(&frame.snapshot);
        if local_hash != remote_hash {
            Some(DesyncReport {
                tick,
                local_hash,
                remote_hash,
            })
        } else {
            None
        }
    }

    // -----------------------------------------------------------------------
    // Private helpers
    // -----------------------------------------------------------------------

    /// Collect a `HashMap<PeerId, Cmd>` into a `Vec<Cmd>` sorted by `PeerId`.
    fn inputs_to_sorted_slice(inputs: &HashMap<PeerId, Cmd>) -> Vec<Cmd> {
        let mut pairs: Vec<(PeerId, &Cmd)> = inputs.iter().map(|(&k, v)| (k, v)).collect();
        pairs.sort_by_key(|(peer, _)| *peer);
        pairs.into_iter().map(|(_, cmd)| cmd.clone()).collect()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // -----------------------------------------------------------------------
    // Toy world
    // -----------------------------------------------------------------------

    #[derive(Debug, Default)]
    struct ToyWorld {
        state: i64,
        history: Vec<i64>,
    }

    fn advance(w: &mut ToyWorld, _dt: f64, cmds: &[i32]) {
        let delta: i32 = cmds.iter().sum();
        w.state += delta as i64;
        w.history.push(w.state);
    }

    fn snapshot(w: &ToyWorld) -> Vec<u8> {
        w.state.to_le_bytes().to_vec()
    }

    fn restore(w: &mut ToyWorld, snap: &[u8]) {
        if snap.len() == 8 {
            w.state = i64::from_le_bytes(snap.try_into().unwrap_or([0u8; 8]));
        }
    }

    fn make_world(capacity: usize) -> RollbackWorld<ToyWorld, i32> {
        RollbackWorld::new(ToyWorld::default(), capacity, advance, snapshot, restore)
    }

    // -----------------------------------------------------------------------
    // Test 1: Late input triggers resimulate
    // -----------------------------------------------------------------------

    #[test]
    fn test_late_input_triggers_resimulate() {
        let mut rw = make_world(16);

        // Tick 1: local cmd = 1 → state becomes 1.
        rw.step(1.0 / 60.0, &[1_i32]);
        assert_eq!(rw.world.state, 1, "after tick 1 state should be 1");
        assert_eq!(rw.current_tick, 1);

        // Late input from peer 2 for tick 1: cmd = 10.
        rw.record_input(2, 1, 10_i32);

        // Resimulate from tick 1 — should replay local(1) + peer2(10) = 11.
        rw.resimulate_from(1, 1.0 / 60.0)
            .expect("resimulate should succeed");

        assert_eq!(rw.world.state, 11, "after resim state should be 11");
    }

    // -----------------------------------------------------------------------
    // Test 2: Identical input → hash equality
    // -----------------------------------------------------------------------

    #[test]
    fn test_identical_input_hash_equality() {
        let cmds: &[i32] = &[3, 7, -2, 1, 5, 0, 4, 8, -3, 2];

        let mut rw1 = make_world(32);
        let mut rw2 = make_world(32);

        let dt = 1.0 / 60.0;
        for &c in cmds {
            rw1.step(dt, &[c]);
            rw2.step(dt, &[c]);
        }

        let tick = rw1.current_tick;
        let snap1 = rw1.buffer.get(tick).map(|f| f.snapshot.clone()).unwrap();
        let snap2 = rw2.buffer.get(tick).map(|f| f.snapshot.clone()).unwrap();

        assert_eq!(
            hash_snapshot(&snap1),
            hash_snapshot(&snap2),
            "identical simulations should produce identical snapshot hashes"
        );
    }

    // -----------------------------------------------------------------------
    // Test 3: Desync detection
    // -----------------------------------------------------------------------

    #[test]
    fn test_desync_detection() {
        let mut rw = make_world(16);
        rw.step(1.0 / 60.0, &[5_i32]);

        let tick = rw.current_tick;

        // Tamper: flip one byte in the stored snapshot.
        {
            let frame = rw.buffer.get_mut(tick).unwrap();
            if let Some(b) = frame.snapshot.first_mut() {
                *b ^= 0xFF;
            }
        }

        // The tampered local hash should differ from the "true" hash.
        let true_state: i64 = 5;
        let true_snap = true_state.to_le_bytes().to_vec();
        let true_hash = hash_snapshot(&true_snap);

        let report = rw.check_desync(tick, true_hash);
        assert!(
            report.is_some(),
            "check_desync should detect the tampered snapshot"
        );

        let r = report.unwrap();
        assert_eq!(r.tick, tick);
        assert_eq!(r.remote_hash, true_hash);
        assert_ne!(r.local_hash, true_hash);
    }

    // -----------------------------------------------------------------------
    // Test 4: Buffer capacity
    // -----------------------------------------------------------------------

    #[test]
    fn test_buffer_capacity() {
        let mut rw = make_world(3);
        let dt = 1.0 / 60.0;

        for _ in 0..10 {
            rw.step(dt, &[1_i32]);
        }

        // capacity = 3, so exactly 3 frames retained.
        assert_eq!(rw.buffer.len(), 3, "buffer should hold at most 3 frames");

        // Oldest frame should not be tick 0 (it got evicted).
        let oldest = rw.buffer.tail_tick().unwrap();
        assert!(oldest > 0, "oldest tick should be > 0 after 10 steps");
    }

    // -----------------------------------------------------------------------
    // Test 5: Serde Frame<i32>
    // -----------------------------------------------------------------------

    #[test]
    fn test_serde_frame() {
        let mut inputs = HashMap::new();
        inputs.insert(1u64, 42_i32);
        inputs.insert(2u64, -7_i32);

        let frame: Frame<i32> = Frame {
            tick: 99,
            snapshot: vec![0xDE, 0xAD, 0xBE, 0xEF],
            inputs,
        };

        let json = serde_json::to_string(&frame).expect("serialize frame");
        let restored: Frame<i32> = serde_json::from_str(&json).expect("deserialize frame");

        assert_eq!(restored.tick, 99);
        assert_eq!(restored.inputs.get(&1u64), Some(&42_i32));
        assert_eq!(restored.inputs.get(&2u64), Some(&-7_i32));
        assert_eq!(restored.snapshot, vec![0xDE, 0xAD, 0xBE, 0xEF]);
    }

    // -----------------------------------------------------------------------
    // Test 6: resimulate_from on missing tick → TickNotInBuffer
    // -----------------------------------------------------------------------

    #[test]
    fn test_resimulate_from_missing_tick() {
        let mut rw = make_world(4);
        rw.step(1.0 / 60.0, &[1_i32]);

        // Tick 999 is not in the buffer.
        let result = rw.resimulate_from(999, 1.0 / 60.0);
        assert_eq!(result, Err(RollbackError::TickNotInBuffer(999)));
    }
}
