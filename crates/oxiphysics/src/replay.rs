// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Deterministic simulation replay: record all external inputs and replay them.
//!
//! The replay system has two components:
//!
//! - **`SimRecorder`** — wraps a simulation, records every external command
//!   applied to it (forces, impulses, velocity overrides, position resets),
//!   and produces a `ReplayRecord` on completion.
//!
//! - **`SimReplayer`** — consumes a `ReplayRecord` and, step by step,
//!   yields the `ReplayStep`s that were originally recorded.  The caller is
//!   responsible for re-applying the commands to a fresh simulation.
//!
//! ## Guarantees
//!
//! Replay is **deterministically reproducible** when:
//! 1. The simulation is seeded with the same random seed.
//! 2. All external inputs are applied via the recorder.
//! 3. The simulation uses the same time steps.
//!
//! ## Storage format
//!
//! Records can be serialised to and from a compact JSON string via
//! `ReplayRecord::to_json` / `ReplayRecord::from_json`, suitable for
//! logging, debugging, and regression testing.
//!
//! ## Example
//!
//! ```rust
//! use oxiphysics::replay::{SimRecorder, SimReplayer, ReplayCommand};
//!
//! let mut rec = SimRecorder::new(42);
//!
//! // Simulate 3 steps, recording the inputs
//! for step in 0..3 {
//!     rec.begin_step(1.0 / 60.0);
//!     rec.record_force(0, [0.0, 100.0, 0.0]);
//!     rec.record_impulse(1, [1.0, 0.0, 0.0]);
//!     rec.end_step();
//! }
//! let record = rec.finish();
//! assert_eq!(record.step_count(), 3);
//!
//! // Replay
//! let mut replayer = SimReplayer::new(record);
//! while let Some(step) = replayer.advance() {
//!     for cmd in &step.commands {
//!         let _ = cmd; // re-apply each command to the simulation
//!     }
//! }
//! assert!(replayer.is_done());
//! ```

#![allow(missing_docs)]
#![allow(dead_code)]

use serde::{Deserialize, Serialize};

// ============================================================================
// ReplayCommand
// ============================================================================

/// A single external input applied to the simulation during one step.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum ReplayCommand {
    /// Apply a continuous force \[N\] to a body for this step.
    ApplyForce {
        /// Zero-based index of the target body.
        body_index: usize,
        /// Force vector in world space \[N\].
        force: [f64; 3],
    },

    /// Apply an instantaneous impulse \[N·s\] to a body.
    ApplyImpulse {
        body_index: usize,
        /// Impulse vector in world space \[N·s\].
        impulse: [f64; 3],
    },

    /// Apply a torque \[N·m\] to a body for this step.
    ApplyTorque {
        body_index: usize,
        /// Torque vector in world space \[N·m\].
        torque: [f64; 3],
    },

    /// Directly override a body's linear velocity.
    SetVelocity {
        body_index: usize,
        velocity: [f64; 3],
    },

    /// Directly override a body's angular velocity.
    SetAngularVelocity {
        body_index: usize,
        angular_velocity: [f64; 3],
    },

    /// Teleport a body to a new position (used for kinematic controllers).
    SetPosition {
        body_index: usize,
        position: [f64; 3],
    },

    /// Wake a sleeping body.
    Wake { body_index: usize },
}

impl ReplayCommand {
    /// Return the index of the body this command targets.
    pub fn body_index(&self) -> usize {
        match self {
            ReplayCommand::ApplyForce { body_index, .. } => *body_index,
            ReplayCommand::ApplyImpulse { body_index, .. } => *body_index,
            ReplayCommand::ApplyTorque { body_index, .. } => *body_index,
            ReplayCommand::SetVelocity { body_index, .. } => *body_index,
            ReplayCommand::SetAngularVelocity { body_index, .. } => *body_index,
            ReplayCommand::SetPosition { body_index, .. } => *body_index,
            ReplayCommand::Wake { body_index } => *body_index,
        }
    }

    /// Human-readable command kind name.
    pub fn kind_name(&self) -> &'static str {
        match self {
            ReplayCommand::ApplyForce { .. } => "ApplyForce",
            ReplayCommand::ApplyImpulse { .. } => "ApplyImpulse",
            ReplayCommand::ApplyTorque { .. } => "ApplyTorque",
            ReplayCommand::SetVelocity { .. } => "SetVelocity",
            ReplayCommand::SetAngularVelocity { .. } => "SetAngularVelocity",
            ReplayCommand::SetPosition { .. } => "SetPosition",
            ReplayCommand::Wake { .. } => "Wake",
        }
    }
}

// ============================================================================
// ReplayStep
// ============================================================================

/// All commands recorded during a single simulation step.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ReplayStep {
    /// Sequential step counter (0-based).
    pub step_index: u64,
    /// Duration of this step \[s\].
    pub dt: f64,
    /// Commands applied during this step, in the order they were recorded.
    pub commands: Vec<ReplayCommand>,
}

impl ReplayStep {
    /// Total number of body-targeting commands in this step.
    pub fn command_count(&self) -> usize {
        self.commands.len()
    }

    /// Returns all forces recorded in this step as `(body_index, force)` pairs.
    pub fn forces(&self) -> impl Iterator<Item = (usize, [f64; 3])> + '_ {
        self.commands.iter().filter_map(|c| {
            if let ReplayCommand::ApplyForce { body_index, force } = c {
                Some((*body_index, *force))
            } else {
                None
            }
        })
    }
}

// ============================================================================
// ReplayRecord
// ============================================================================

/// A complete recording of all external inputs over a simulation run.
///
/// Can be serialised to JSON with [`to_json`][ReplayRecord::to_json] and
/// deserialised with [`from_json`][ReplayRecord::from_json].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ReplayRecord {
    /// Format version for future compatibility.
    pub version: u32,
    /// Random seed used when the simulation was first created.
    pub seed: u64,
    /// Human-readable label (e.g. "crash test scenario 3").
    pub label: String,
    /// All recorded steps in order.
    pub steps: Vec<ReplayStep>,
}

impl ReplayRecord {
    /// Create a new empty record.
    pub fn new(seed: u64) -> Self {
        Self {
            version: 1,
            seed,
            label: String::new(),
            steps: Vec::new(),
        }
    }

    /// Create a new empty record with a label.
    pub fn with_label(seed: u64, label: impl Into<String>) -> Self {
        Self {
            version: 1,
            seed,
            label: label.into(),
            steps: Vec::new(),
        }
    }

    /// Total number of recorded steps.
    pub fn step_count(&self) -> usize {
        self.steps.len()
    }

    /// Total simulated time (sum of all `dt` values).
    pub fn total_sim_time(&self) -> f64 {
        self.steps.iter().map(|s| s.dt).sum()
    }

    /// Total number of commands across all steps.
    pub fn total_commands(&self) -> usize {
        self.steps.iter().map(|s| s.commands.len()).sum()
    }

    /// Returns the step at the given index, or `None` if out of bounds.
    pub fn step(&self, index: usize) -> Option<&ReplayStep> {
        self.steps.get(index)
    }

    /// Serialise this record to a JSON string.
    pub fn to_json(&self) -> Result<String, String> {
        serde_json::to_string_pretty(self).map_err(|e| e.to_string())
    }

    /// Deserialise a record from a JSON string produced by [`Self::to_json`].
    pub fn from_json(s: &str) -> Result<Self, String> {
        serde_json::from_str(s).map_err(|e| e.to_string())
    }

    /// Serialise to compact (non-pretty) JSON.
    pub fn to_json_compact(&self) -> Result<String, String> {
        serde_json::to_string(self).map_err(|e| e.to_string())
    }

    /// Build a summary of this record.
    pub fn summary(&self) -> String {
        format!(
            "ReplayRecord {{ version={}, seed={}, label={:?}, steps={}, \
             total_time={:.4}s, total_commands={} }}",
            self.version,
            self.seed,
            self.label,
            self.step_count(),
            self.total_sim_time(),
            self.total_commands()
        )
    }
}

// ============================================================================
// SimRecorder
// ============================================================================

/// Records external inputs applied to a simulation.
///
/// Call [`begin_step`][SimRecorder::begin_step] before the simulation step,
/// record commands via the `record_*` methods, then call
/// [`end_step`][SimRecorder::end_step].  When the simulation is done call
/// [`finish`][SimRecorder::finish] to take ownership of the [`ReplayRecord`].
#[derive(Debug)]
pub struct SimRecorder {
    record: ReplayRecord,
    current: Option<ReplayStep>,
}

impl SimRecorder {
    /// Create a new recorder.  `seed` should match the simulation's RNG seed.
    pub fn new(seed: u64) -> Self {
        Self {
            record: ReplayRecord::new(seed),
            current: None,
        }
    }

    /// Create a new recorder with a label.
    pub fn with_label(seed: u64, label: impl Into<String>) -> Self {
        Self {
            record: ReplayRecord::with_label(seed, label),
            current: None,
        }
    }

    /// Mark the start of a new step with time step `dt`.
    ///
    /// Must be followed by zero or more `record_*` calls and then
    /// [`end_step`][Self::end_step].
    pub fn begin_step(&mut self, dt: f64) {
        debug_assert!(self.current.is_none(), "begin_step called without end_step");
        let step_index = self.record.steps.len() as u64;
        self.current = Some(ReplayStep {
            step_index,
            dt,
            commands: Vec::new(),
        });
    }

    fn push(&mut self, cmd: ReplayCommand) {
        if let Some(s) = &mut self.current {
            s.commands.push(cmd);
        }
    }

    /// Record a force applied to `body_index`.
    pub fn record_force(&mut self, body_index: usize, force: [f64; 3]) {
        self.push(ReplayCommand::ApplyForce { body_index, force });
    }

    /// Record an impulse applied to `body_index`.
    pub fn record_impulse(&mut self, body_index: usize, impulse: [f64; 3]) {
        self.push(ReplayCommand::ApplyImpulse {
            body_index,
            impulse,
        });
    }

    /// Record a torque applied to `body_index`.
    pub fn record_torque(&mut self, body_index: usize, torque: [f64; 3]) {
        self.push(ReplayCommand::ApplyTorque { body_index, torque });
    }

    /// Record an explicit velocity override on `body_index`.
    pub fn record_set_velocity(&mut self, body_index: usize, velocity: [f64; 3]) {
        self.push(ReplayCommand::SetVelocity {
            body_index,
            velocity,
        });
    }

    /// Record an explicit angular-velocity override on `body_index`.
    pub fn record_set_angular_velocity(&mut self, body_index: usize, av: [f64; 3]) {
        self.push(ReplayCommand::SetAngularVelocity {
            body_index,
            angular_velocity: av,
        });
    }

    /// Record a teleport (position override) of `body_index`.
    pub fn record_set_position(&mut self, body_index: usize, position: [f64; 3]) {
        self.push(ReplayCommand::SetPosition {
            body_index,
            position,
        });
    }

    /// Record a wake command for a sleeping body.
    pub fn record_wake(&mut self, body_index: usize) {
        self.push(ReplayCommand::Wake { body_index });
    }

    /// Commit the current step to the record.
    ///
    /// Must be called after [`begin_step`][Self::begin_step].
    pub fn end_step(&mut self) {
        if let Some(step) = self.current.take() {
            self.record.steps.push(step);
        }
    }

    /// Finish recording and return the completed [`ReplayRecord`].
    ///
    /// If a step was begun but not ended, it is committed automatically.
    pub fn finish(mut self) -> ReplayRecord {
        if self.current.is_some() {
            self.end_step();
        }
        self.record
    }

    /// Return the number of steps recorded so far.
    pub fn recorded_steps(&self) -> usize {
        self.record.steps.len()
    }

    /// Return a reference to the in-progress record (for inspection).
    pub fn record(&self) -> &ReplayRecord {
        &self.record
    }
}

// ============================================================================
// SimReplayer
// ============================================================================

/// Iterates through a [`ReplayRecord`] step by step.
///
/// Call [`advance`][SimReplayer::advance] once per simulation step to get the
/// [`ReplayStep`] containing all commands that were recorded for that step.
/// Apply the commands to your simulation (same order as recorded) to reproduce
/// the original run.
#[derive(Debug)]
pub struct SimReplayer {
    record: ReplayRecord,
    cursor: usize,
}

impl SimReplayer {
    /// Create a replayer for the given record.
    pub fn new(record: ReplayRecord) -> Self {
        Self { record, cursor: 0 }
    }

    /// Returns `true` when all steps have been yielded.
    pub fn is_done(&self) -> bool {
        self.cursor >= self.record.steps.len()
    }

    /// Return the current step without advancing.
    pub fn current_step(&self) -> Option<&ReplayStep> {
        self.record.steps.get(self.cursor)
    }

    /// Advance to the next step, returning a reference to it (or `None` if done).
    pub fn advance(&mut self) -> Option<&ReplayStep> {
        let step = self.record.steps.get(self.cursor)?;
        self.cursor += 1;
        Some(step)
    }

    /// Number of steps not yet yielded.
    pub fn steps_remaining(&self) -> usize {
        self.record.steps.len().saturating_sub(self.cursor)
    }

    /// Total steps in the record.
    pub fn total_steps(&self) -> usize {
        self.record.steps.len()
    }

    /// Reset the cursor to the beginning.
    pub fn reset(&mut self) {
        self.cursor = 0;
    }

    /// Seek to a specific step index.  Returns `false` if out of bounds.
    pub fn seek(&mut self, step: usize) -> bool {
        if step <= self.record.steps.len() {
            self.cursor = step;
            true
        } else {
            false
        }
    }

    /// Consume the replayer and return the underlying record.
    pub fn into_record(self) -> ReplayRecord {
        self.record
    }

    /// Return a reference to the underlying record.
    pub fn record(&self) -> &ReplayRecord {
        &self.record
    }

    /// Return the current cursor position (0-based step index).
    pub fn cursor(&self) -> usize {
        self.cursor
    }
}

impl Iterator for SimReplayer {
    type Item = ReplayStep;

    fn next(&mut self) -> Option<Self::Item> {
        if self.cursor < self.record.steps.len() {
            let step = self.record.steps[self.cursor].clone();
            self.cursor += 1;
            Some(step)
        } else {
            None
        }
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    fn make_record() -> ReplayRecord {
        let mut rec = SimRecorder::new(12345);
        for i in 0..5u64 {
            rec.begin_step(1.0 / 60.0);
            rec.record_force(0, [0.0, f64::from(i as u32) * 10.0, 0.0]);
            if i % 2 == 0 {
                rec.record_impulse(1, [1.0, 0.0, 0.0]);
            }
            rec.end_step();
        }
        rec.finish()
    }

    #[test]
    fn record_step_count() {
        let r = make_record();
        assert_eq!(r.step_count(), 5);
    }

    #[test]
    fn record_total_commands() {
        let r = make_record();
        // Steps 0,2,4 have 2 commands; steps 1,3 have 1 → 3*2 + 2*1 = 8
        assert_eq!(r.total_commands(), 8);
    }

    #[test]
    fn record_total_sim_time() {
        let r = make_record();
        let expected = 5.0 / 60.0;
        assert!((r.total_sim_time() - expected).abs() < 1e-14);
    }

    #[test]
    fn replayer_all_steps() {
        let r = make_record();
        let mut rep = SimReplayer::new(r);
        let mut count = 0;
        while let Some(_s) = rep.advance() {
            count += 1;
        }
        assert_eq!(count, 5);
        assert!(rep.is_done());
    }

    #[test]
    fn replayer_reset() {
        let r = make_record();
        let mut rep = SimReplayer::new(r);
        assert!(rep.advance().is_some());
        rep.reset();
        assert_eq!(rep.cursor(), 0);
        assert!(!rep.is_done());
    }

    #[test]
    fn replayer_seek() {
        let r = make_record();
        let mut rep = SimReplayer::new(r);
        assert!(rep.seek(3));
        assert_eq!(rep.cursor(), 3);
        assert_eq!(rep.steps_remaining(), 2);
    }

    #[test]
    fn replayer_as_iterator() {
        let r = make_record();
        let steps: Vec<_> = SimReplayer::new(r).collect();
        assert_eq!(steps.len(), 5);
        assert_eq!(steps[0].step_index, 0);
        assert_eq!(steps[4].step_index, 4);
    }

    #[test]
    fn json_round_trip() {
        let r = make_record();
        let json = r.to_json().expect("serialise");
        let r2 = ReplayRecord::from_json(&json).expect("deserialise");
        assert_eq!(r, r2);
    }

    #[test]
    fn json_compact_round_trip() {
        let r = make_record();
        let json = r.to_json_compact().expect("compact serialise");
        let r2 = ReplayRecord::from_json(&json).expect("deserialise");
        assert_eq!(r, r2);
    }

    #[test]
    fn summary_contains_step_count() {
        let r = make_record();
        let s = r.summary();
        assert!(s.contains("steps=5"));
    }

    #[test]
    fn command_kind_names() {
        let cmds = [
            ReplayCommand::ApplyForce {
                body_index: 0,
                force: [0.0; 3],
            },
            ReplayCommand::ApplyImpulse {
                body_index: 0,
                impulse: [0.0; 3],
            },
            ReplayCommand::SetPosition {
                body_index: 0,
                position: [0.0; 3],
            },
            ReplayCommand::Wake { body_index: 0 },
        ];
        assert_eq!(cmds[0].kind_name(), "ApplyForce");
        assert_eq!(cmds[1].kind_name(), "ApplyImpulse");
        assert_eq!(cmds[2].kind_name(), "SetPosition");
        assert_eq!(cmds[3].kind_name(), "Wake");
    }

    #[test]
    fn step_forces_iterator() {
        let mut rec = SimRecorder::new(0);
        rec.begin_step(0.01);
        rec.record_force(0, [1.0, 0.0, 0.0]);
        rec.record_force(1, [0.0, 2.0, 0.0]);
        rec.record_impulse(0, [0.0, 0.0, 1.0]); // not a force
        rec.end_step();
        let rec = rec.finish();
        let forces: Vec<_> = rec.steps[0].forces().collect();
        assert_eq!(forces.len(), 2);
        assert_eq!(forces[0], (0, [1.0, 0.0, 0.0]));
    }

    #[test]
    fn recorder_finish_auto_end_step() {
        let mut rec = SimRecorder::new(0);
        rec.begin_step(0.01);
        rec.record_force(0, [0.0; 3]);
        // Deliberately skip end_step — finish() should commit it
        let r = rec.finish();
        assert_eq!(r.step_count(), 1);
    }
}
