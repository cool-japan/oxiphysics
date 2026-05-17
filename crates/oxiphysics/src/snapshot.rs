// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! World-state snapshots with delta tracking and a ring-buffer manager.
//!
//! `WorldSnapshot` captures a complete set of body states at a given
//! simulation step.  `SnapshotManager` stores a bounded history of
//! snapshots and supports step-indexed lookup and inter-snapshot delta
//! computation.
//!
//! ## Typical workflow
//!
//! ```rust
//! use oxiphysics::snapshot::{BodySnapshot, WorldSnapshot, SnapshotManager};
//!
//! let mut manager = SnapshotManager::new(10);   // keep last 10 snapshots
//!
//! // After each simulation step…
//! let snap = WorldSnapshot::new(1, 1.0 / 60.0)
//!     .add_body(BodySnapshot::at_rest(0, [0.0, 0.0, 0.0]));
//! manager.push(snap);
//!
//! assert_eq!(manager.len(), 1);
//! assert!(manager.latest().is_some());
//! ```

#![allow(missing_docs)]
#![allow(dead_code)]

use serde::{Deserialize, Serialize};

// ============================================================================
// Helper — monotonic nanoseconds
// ============================================================================

fn monotonic_ns() -> u64 {
    use std::time::{SystemTime, UNIX_EPOCH};
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_nanos() as u64)
        .unwrap_or(0)
}

// ============================================================================
// BodySnapshot
// ============================================================================

/// State of a single body captured at a snapshot instant.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BodySnapshot {
    /// Application-assigned body index.
    pub index: usize,
    /// Optional human-readable label.
    #[serde(default)]
    pub label: String,
    /// World-space position [x, y, z] in metres.
    pub position: [f64; 3],
    /// World-space orientation as unit quaternion [qx, qy, qz, qw].
    pub rotation: [f64; 4],
    /// Linear velocity in m/s.
    #[serde(default)]
    pub velocity: [f64; 3],
    /// Angular velocity in rad/s.
    #[serde(default)]
    pub angular_velocity: [f64; 3],
    /// Whether the body is in the sleeping (deactivated) state.
    #[serde(default)]
    pub is_sleeping: bool,
    /// Whether the body is kinematic / static.
    #[serde(default)]
    pub is_static: bool,
}

impl BodySnapshot {
    /// Create a snapshot for a body at rest.
    pub fn at_rest(index: usize, position: [f64; 3]) -> Self {
        Self {
            index,
            label: String::new(),
            position,
            rotation: [0.0, 0.0, 0.0, 1.0],
            velocity: [0.0; 3],
            angular_velocity: [0.0; 3],
            is_sleeping: false,
            is_static: false,
        }
    }

    /// Create with full state.
    pub fn new(
        index: usize,
        position: [f64; 3],
        rotation: [f64; 4],
        velocity: [f64; 3],
        angular_velocity: [f64; 3],
    ) -> Self {
        Self {
            index,
            label: String::new(),
            position,
            rotation,
            velocity,
            angular_velocity,
            is_sleeping: false,
            is_static: false,
        }
    }

    /// Builder: attach a label.
    pub fn with_label(mut self, label: impl Into<String>) -> Self {
        self.label = label.into();
        self
    }

    /// Builder: mark sleeping.
    pub fn sleeping(mut self) -> Self {
        self.is_sleeping = true;
        self
    }

    /// Builder: mark static.
    pub fn as_static(mut self) -> Self {
        self.is_static = true;
        self
    }

    /// Linear speed (magnitude of velocity) in m/s.
    pub fn speed(&self) -> f64 {
        let v = &self.velocity;
        (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
    }

    /// Angular speed (magnitude of angular velocity) in rad/s.
    pub fn angular_speed(&self) -> f64 {
        let w = &self.angular_velocity;
        (w[0] * w[0] + w[1] * w[1] + w[2] * w[2]).sqrt()
    }

    /// Approximate kinetic energy (requires mass as input) in J.
    pub fn kinetic_energy(&self, mass: f64, inertia: f64) -> f64 {
        0.5 * mass * self.speed().powi(2) + 0.5 * inertia * self.angular_speed().powi(2)
    }
}

// ============================================================================
// BodyDelta
// ============================================================================

/// Change in a body's state between two snapshots.
#[derive(Debug, Clone, PartialEq)]
pub struct BodyDelta {
    /// Body index.
    pub index: usize,
    /// Change in position [dx, dy, dz].
    pub position_delta: [f64; 3],
    /// Distance moved (magnitude of position_delta).
    pub distance_moved: f64,
    /// Linear speed in the *later* snapshot.
    pub speed: f64,
    /// Angular speed in the *later* snapshot.
    pub angular_speed: f64,
}

impl BodyDelta {
    /// Returns `true` if this body barely moved.
    pub fn is_still(&self, threshold: f64) -> bool {
        self.distance_moved < threshold && self.speed < threshold
    }
}

// ============================================================================
// WorldSnapshot
// ============================================================================

/// Complete state of all bodies at a given simulation step.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct WorldSnapshot {
    /// Simulation step index.
    pub step: u64,
    /// Simulated time at this snapshot (seconds).
    pub sim_time: f64,
    /// Optional human-readable label.
    #[serde(default)]
    pub label: String,
    /// All body states.
    pub bodies: Vec<BodySnapshot>,
    /// Monotonic wall-clock timestamp in nanoseconds (for ordering).
    #[serde(default)]
    pub created_at_ns: u64,
}

impl WorldSnapshot {
    /// Create an empty snapshot.
    pub fn new(step: u64, sim_time: f64) -> Self {
        Self {
            step,
            sim_time,
            label: String::new(),
            bodies: Vec::new(),
            created_at_ns: monotonic_ns(),
        }
    }

    /// Builder: add a body state.
    pub fn add_body(mut self, body: BodySnapshot) -> Self {
        self.bodies.push(body);
        self
    }

    /// Builder: attach a label.
    pub fn with_label(mut self, label: impl Into<String>) -> Self {
        self.label = label.into();
        self
    }

    /// Add multiple bodies at once.
    pub fn extend_bodies(&mut self, bodies: impl IntoIterator<Item = BodySnapshot>) {
        self.bodies.extend(bodies);
    }

    /// Number of body entries.
    pub fn body_count(&self) -> usize {
        self.bodies.len()
    }

    /// Look up a body snapshot by application index.
    pub fn body_at_index(&self, index: usize) -> Option<&BodySnapshot> {
        self.bodies.iter().find(|b| b.index == index)
    }

    /// All sleeping bodies.
    pub fn sleeping_bodies(&self) -> Vec<&BodySnapshot> {
        self.bodies.iter().filter(|b| b.is_sleeping).collect()
    }

    /// All dynamic (non-static) bodies.
    pub fn dynamic_bodies(&self) -> Vec<&BodySnapshot> {
        self.bodies.iter().filter(|b| !b.is_static).collect()
    }

    /// Compute per-body deltas between this snapshot and an earlier one.
    ///
    /// Bodies present in `earlier` but not here are omitted.  New bodies
    /// (present only in `self`) are also omitted — this only tracks bodies
    /// present in *both* snapshots.
    pub fn body_delta(&self, earlier: &WorldSnapshot) -> Vec<BodyDelta> {
        self.bodies
            .iter()
            .filter_map(|body| {
                let prev = earlier.body_at_index(body.index)?;
                let dp = [
                    body.position[0] - prev.position[0],
                    body.position[1] - prev.position[1],
                    body.position[2] - prev.position[2],
                ];
                let distance_moved = (dp[0] * dp[0] + dp[1] * dp[1] + dp[2] * dp[2]).sqrt();
                Some(BodyDelta {
                    index: body.index,
                    position_delta: dp,
                    distance_moved,
                    speed: body.speed(),
                    angular_speed: body.angular_speed(),
                })
            })
            .collect()
    }

    /// Total kinetic energy estimate (needs mass/inertia per body — uses 1.0 as unit mass).
    pub fn total_kinetic_energy_unit_mass(&self) -> f64 {
        self.bodies
            .iter()
            .filter(|b| !b.is_static)
            .map(|b| b.kinetic_energy(1.0, 1.0))
            .sum()
    }

    /// Deserialize from JSON.
    pub fn from_json(s: &str) -> Result<Self, String> {
        serde_json::from_str(s).map_err(|e| e.to_string())
    }

    /// Serialize to pretty-printed JSON.
    pub fn to_json(&self) -> Result<String, String> {
        serde_json::to_string_pretty(self).map_err(|e| e.to_string())
    }

    /// Serialize to compact JSON.
    pub fn to_json_compact(&self) -> Result<String, String> {
        serde_json::to_string(self).map_err(|e| e.to_string())
    }

    /// One-line summary.
    pub fn summary(&self) -> String {
        format!(
            "WorldSnapshot {{ step={}, sim_time={:.4}s, bodies={}, sleeping={} }}",
            self.step,
            self.sim_time,
            self.body_count(),
            self.sleeping_bodies().len(),
        )
    }
}

impl std::fmt::Display for WorldSnapshot {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.summary())
    }
}

// ============================================================================
// SnapshotDiff
// ============================================================================

/// Aggregated diff between two world snapshots.
#[derive(Debug, Clone)]
pub struct SnapshotDiff {
    /// Step distance between the two snapshots.
    pub step_delta: u64,
    /// Simulated time elapsed between them.
    pub sim_time_delta: f64,
    /// Per-body deltas.
    pub body_deltas: Vec<BodyDelta>,
    /// Number of bodies that changed position by more than `threshold`.
    pub moving_body_count: usize,
}

impl SnapshotDiff {
    fn compute(later: &WorldSnapshot, earlier: &WorldSnapshot, still_threshold: f64) -> Self {
        let body_deltas = later.body_delta(earlier);
        let moving_body_count = body_deltas
            .iter()
            .filter(|d| !d.is_still(still_threshold))
            .count();
        Self {
            step_delta: later.step.saturating_sub(earlier.step),
            sim_time_delta: later.sim_time - earlier.sim_time,
            body_deltas,
            moving_body_count,
        }
    }
}

// ============================================================================
// SnapshotManager
// ============================================================================

/// A bounded ring-buffer of [`WorldSnapshot`]s, ordered by step index.
///
/// When the capacity is exceeded the oldest snapshot is dropped to make room.
///
/// # Example
///
/// ```rust
/// use oxiphysics::snapshot::{BodySnapshot, WorldSnapshot, SnapshotManager};
///
/// let mut manager = SnapshotManager::new(5);
/// for step in 0..8_u64 {
///     let snap = WorldSnapshot::new(step, step as f64 / 60.0)
///         .add_body(BodySnapshot::at_rest(0, [0.0, step as f64, 0.0]));
///     manager.push(snap);
/// }
/// // Only the last 5 are kept
/// assert_eq!(manager.len(), 5);
/// assert_eq!(manager.latest().unwrap().step, 7);
/// ```
#[derive(Debug)]
pub struct SnapshotManager {
    snapshots: Vec<WorldSnapshot>,
    max_count: usize,
}

impl SnapshotManager {
    /// Create a manager that keeps at most `max_count` snapshots.
    pub fn new(max_count: usize) -> Self {
        Self {
            snapshots: Vec::with_capacity(max_count),
            max_count,
        }
    }

    // ------------------------------------------------------------------
    // Insertion
    // ------------------------------------------------------------------

    /// Push a new snapshot.  If capacity is exceeded the oldest is dropped.
    pub fn push(&mut self, snapshot: WorldSnapshot) {
        if self.snapshots.len() >= self.max_count {
            self.snapshots.remove(0);
        }
        self.snapshots.push(snapshot);
    }

    // ------------------------------------------------------------------
    // Lookup
    // ------------------------------------------------------------------

    /// The most recently pushed snapshot.
    pub fn latest(&self) -> Option<&WorldSnapshot> {
        self.snapshots.last()
    }

    /// Exact lookup by step index.
    pub fn at_step(&self, step: u64) -> Option<&WorldSnapshot> {
        self.snapshots.iter().find(|s| s.step == step)
    }

    /// The snapshot whose step is closest to `step` (ties: prefer earlier).
    pub fn closest_to_step(&self, step: u64) -> Option<&WorldSnapshot> {
        self.snapshots
            .iter()
            .min_by_key(|s| s.step.max(step) - s.step.min(step))
    }

    /// The earliest stored snapshot.
    pub fn earliest(&self) -> Option<&WorldSnapshot> {
        self.snapshots.first()
    }

    // ------------------------------------------------------------------
    // Pruning
    // ------------------------------------------------------------------

    /// Drop all snapshots with step < `min_step`.  Returns count removed.
    pub fn prune_older_than(&mut self, min_step: u64) -> usize {
        let before = self.snapshots.len();
        self.snapshots.retain(|s| s.step >= min_step);
        before - self.snapshots.len()
    }

    /// Drop all snapshots older than `keep_last` seconds of sim time.
    pub fn prune_old_sim_time(&mut self, current_sim_time: f64, keep_last: f64) -> usize {
        let cutoff = current_sim_time - keep_last;
        let before = self.snapshots.len();
        self.snapshots.retain(|s| s.sim_time >= cutoff);
        before - self.snapshots.len()
    }

    /// Remove all snapshots.
    pub fn clear(&mut self) {
        self.snapshots.clear();
    }

    // ------------------------------------------------------------------
    // Introspection
    // ------------------------------------------------------------------

    /// Number of snapshots currently held.
    pub fn len(&self) -> usize {
        self.snapshots.len()
    }

    /// Returns `true` if no snapshots are stored.
    pub fn is_empty(&self) -> bool {
        self.snapshots.is_empty()
    }

    /// The configured maximum capacity.
    pub fn max_count(&self) -> usize {
        self.max_count
    }

    /// Sorted list of step indices.
    pub fn steps(&self) -> Vec<u64> {
        let mut steps: Vec<u64> = self.snapshots.iter().map(|s| s.step).collect();
        steps.sort_unstable();
        steps
    }

    /// Iterate over all snapshots (oldest first).
    pub fn iter(&self) -> impl Iterator<Item = &WorldSnapshot> {
        self.snapshots.iter()
    }

    // ------------------------------------------------------------------
    // Diff
    // ------------------------------------------------------------------

    /// Compute the diff between the latest and the snapshot immediately before it.
    ///
    /// Returns `None` if fewer than 2 snapshots are stored.
    pub fn last_diff(&self, still_threshold: f64) -> Option<SnapshotDiff> {
        if self.snapshots.len() < 2 {
            return None;
        }
        let n = self.snapshots.len();
        Some(SnapshotDiff::compute(
            &self.snapshots[n - 1],
            &self.snapshots[n - 2],
            still_threshold,
        ))
    }

    /// Diff between two specific step indices.
    ///
    /// Returns `None` if either step is not found, or if `later_step ≤ earlier_step`.
    pub fn diff_between(
        &self,
        earlier_step: u64,
        later_step: u64,
        still_threshold: f64,
    ) -> Option<SnapshotDiff> {
        if later_step <= earlier_step {
            return None;
        }
        let earlier = self.at_step(earlier_step)?;
        let later = self.at_step(later_step)?;
        Some(SnapshotDiff::compute(later, earlier, still_threshold))
    }

    // ------------------------------------------------------------------
    // Summary
    // ------------------------------------------------------------------

    /// Multi-line diagnostic string.
    pub fn summary(&self) -> String {
        let step_range = if self.snapshots.is_empty() {
            "empty".to_string()
        } else {
            format!(
                "steps {}..={}",
                self.snapshots
                    .first()
                    .expect("non-empty checked above")
                    .step,
                self.snapshots.last().expect("non-empty checked above").step
            )
        };
        format!(
            "SnapshotManager {{ capacity={}, stored={}, {} }}",
            self.max_count,
            self.len(),
            step_range
        )
    }
}

impl std::fmt::Display for SnapshotManager {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.summary())
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    fn make_snap(step: u64) -> WorldSnapshot {
        WorldSnapshot::new(step, step as f64 * (1.0 / 60.0))
            .add_body(BodySnapshot::at_rest(0, [0.0, step as f64 * 0.01, 0.0]))
    }

    #[test]
    fn ring_buffer_evicts_oldest() {
        let mut mgr = SnapshotManager::new(3);
        for i in 0..5_u64 {
            mgr.push(make_snap(i));
        }
        assert_eq!(mgr.len(), 3);
        assert_eq!(mgr.earliest().unwrap().step, 2);
        assert_eq!(mgr.latest().unwrap().step, 4);
    }

    #[test]
    fn at_step_lookup() {
        let mut mgr = SnapshotManager::new(10);
        for i in 0..5_u64 {
            mgr.push(make_snap(i));
        }
        assert!(mgr.at_step(2).is_some());
        assert!(mgr.at_step(99).is_none());
    }

    #[test]
    fn closest_to_step() {
        let mut mgr = SnapshotManager::new(10);
        mgr.push(make_snap(0));
        mgr.push(make_snap(10));
        mgr.push(make_snap(20));
        let s = mgr.closest_to_step(12).unwrap();
        assert_eq!(s.step, 10);
    }

    #[test]
    fn prune_older_than() {
        let mut mgr = SnapshotManager::new(10);
        for i in 0..6_u64 {
            mgr.push(make_snap(i));
        }
        let removed = mgr.prune_older_than(3);
        assert_eq!(removed, 3);
        assert_eq!(mgr.len(), 3);
        assert_eq!(mgr.earliest().unwrap().step, 3);
    }

    #[test]
    fn body_delta() {
        let snap_a = WorldSnapshot::new(0, 0.0).add_body(BodySnapshot::new(
            0,
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
            [1.0, 0.0, 0.0],
            [0.0; 3],
        ));
        let snap_b = WorldSnapshot::new(1, 1.0).add_body(BodySnapshot::new(
            0,
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
            [1.0, 0.0, 0.0],
            [0.0; 3],
        ));

        let deltas = snap_b.body_delta(&snap_a);
        assert_eq!(deltas.len(), 1);
        assert!((deltas[0].distance_moved - 1.0).abs() < 1e-6);
    }

    #[test]
    fn last_diff() {
        let mut mgr = SnapshotManager::new(10);
        mgr.push(WorldSnapshot::new(0, 0.0).add_body(BodySnapshot::new(
            0,
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
            [0.0; 3],
            [0.0; 3],
        )));
        mgr.push(WorldSnapshot::new(1, 1.0).add_body(BodySnapshot::new(
            0,
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
            [1.0, 0.0, 0.0],
            [0.0; 3],
        )));
        let diff = mgr.last_diff(0.01).unwrap();
        assert_eq!(diff.step_delta, 1);
        assert_eq!(diff.moving_body_count, 1);
    }

    #[test]
    fn json_round_trip() {
        let snap = WorldSnapshot::new(42, 0.7)
            .with_label("checkpoint")
            .add_body(BodySnapshot::at_rest(0, [1.0, 2.0, 3.0]));

        let json = snap.to_json().unwrap();
        let restored = WorldSnapshot::from_json(&json).unwrap();
        assert_eq!(restored.step, 42);
        assert_eq!(restored.label, "checkpoint");
        assert_eq!(restored.bodies[0].position, [1.0, 2.0, 3.0]);
    }

    #[test]
    fn snapshot_display() {
        let snap = WorldSnapshot::new(5, 0.083);
        let s = snap.to_string();
        assert!(s.contains("step=5"));
    }
}
