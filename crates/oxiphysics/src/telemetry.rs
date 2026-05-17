// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Physics simulation telemetry — per-step stats, rolling averages, and CSV
//! export.
//!
//! ## Types
//!
//! - `PhysicsStats` — a snapshot of per-step metrics (body count, contacts,
//!   kinetic energy, solver time, …).
//! - `PhysicsAverages` — windowed averages of the continuous metrics (f64
//!   fields only; not averaging `step` or count fields).
//! - `TelemetrySession` — ring-buffer window of `PhysicsStats` with
//!   helpers for peaks, trend analysis, and CSV export.
//!
//! ## Example
//!
//! ```rust,no_run
//! use oxiphysics::telemetry::{PhysicsStats, TelemetrySession};
//!
//! let mut session = TelemetrySession::new(60);
//!
//! session.push(PhysicsStats {
//!     step: 1, dt: 0.016,
//!     body_count: 100, sleeping_count: 20, contact_count: 50,
//!     island_count: 5, solve_iterations: 8, broad_phase_pairs: 200,
//!     kinetic_energy: 42.0, elapsed_ms: 2.1,
//! });
//!
//! let avg = session.average().unwrap();
//! assert!(avg.kinetic_energy > 0.0);
//!
//! println!("{}", session.to_csv());
//! ```

#![allow(missing_docs)]
#![allow(dead_code)]

use serde::{Deserialize, Serialize};
use std::collections::VecDeque;

// ---------------------------------------------------------------------------
// PhysicsStats
// ---------------------------------------------------------------------------

/// Per-step physics simulation metrics.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PhysicsStats {
    /// Simulation step index (monotonically increasing).
    pub step: u64,
    /// Physics time-step duration (seconds).
    pub dt: f64,
    /// Total number of active bodies.
    pub body_count: usize,
    /// Number of bodies currently sleeping.
    pub sleeping_count: usize,
    /// Number of active contact pairs.
    pub contact_count: usize,
    /// Number of simulation islands.
    pub island_count: usize,
    /// Constraint solver iterations executed this step.
    pub solve_iterations: usize,
    /// Broad-phase collision candidate pairs.
    pub broad_phase_pairs: usize,
    /// Total kinetic energy of all bodies (J).
    pub kinetic_energy: f64,
    /// Wall-clock time consumed by this step (milliseconds).
    pub elapsed_ms: f64,
}

// ---------------------------------------------------------------------------
// PhysicsAverages
// ---------------------------------------------------------------------------

/// Windowed averages of continuous per-step metrics.
///
/// Only f64-valued fields are averaged; integer counts and the step index are
/// excluded.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PhysicsAverages {
    /// Mean physics time-step (s).
    pub dt: f64,
    /// Mean kinetic energy (J).
    pub kinetic_energy: f64,
    /// Mean wall-clock step duration (ms).
    pub elapsed_ms: f64,
    /// Number of samples used to compute these averages.
    pub sample_count: usize,
}

// ---------------------------------------------------------------------------
// TelemetrySession
// ---------------------------------------------------------------------------

/// Rolling telemetry window for a running simulation.
///
/// Maintains the last `window_size` [`PhysicsStats`] entries in a
/// `VecDeque`.  When the window is full, the oldest entry is dropped
/// automatically.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TelemetrySession {
    window: VecDeque<PhysicsStats>,
    /// Maximum number of retained samples.
    pub window_size: usize,
    /// Total steps pushed since the session was created.
    pub total_steps: u64,
}

impl TelemetrySession {
    /// Create a new session with the given rolling-window size.
    ///
    /// Panics if `window_size == 0`.
    pub fn new(window_size: usize) -> Self {
        assert!(window_size > 0, "window_size must be positive");
        Self {
            window: VecDeque::with_capacity(window_size),
            window_size,
            total_steps: 0,
        }
    }

    // -----------------------------------------------------------------------
    // Mutation
    // -----------------------------------------------------------------------

    /// Append a stats record.  Drops the oldest entry when the window is full.
    pub fn push(&mut self, stats: PhysicsStats) {
        if self.window.len() == self.window_size {
            self.window.pop_front();
        }
        self.window.push_back(stats);
        self.total_steps += 1;
    }

    /// Remove all entries.
    pub fn clear(&mut self) {
        self.window.clear();
    }

    // -----------------------------------------------------------------------
    // Queries — latest
    // -----------------------------------------------------------------------

    /// Most recent stats record, or `None` if the session is empty.
    pub fn latest(&self) -> Option<&PhysicsStats> {
        self.window.back()
    }

    /// Number of samples currently in the window.
    pub fn len(&self) -> usize {
        self.window.len()
    }

    /// Return `true` if no samples have been pushed.
    pub fn is_empty(&self) -> bool {
        self.window.is_empty()
    }

    /// Iterate over samples from oldest to newest.
    pub fn iter(&self) -> impl Iterator<Item = &PhysicsStats> {
        self.window.iter()
    }

    // -----------------------------------------------------------------------
    // Aggregate statistics
    // -----------------------------------------------------------------------

    /// Windowed averages over the current window.
    ///
    /// Returns `None` if the session is empty.
    pub fn average(&self) -> Option<PhysicsAverages> {
        let n = self.window.len();
        if n == 0 {
            return None;
        }
        let n_f = n as f64;
        let mut sum_dt = 0.0_f64;
        let mut sum_ke = 0.0_f64;
        let mut sum_ms = 0.0_f64;
        for s in &self.window {
            sum_dt += s.dt;
            sum_ke += s.kinetic_energy;
            sum_ms += s.elapsed_ms;
        }
        Some(PhysicsAverages {
            dt: sum_dt / n_f,
            kinetic_energy: sum_ke / n_f,
            elapsed_ms: sum_ms / n_f,
            sample_count: n,
        })
    }

    /// Peak kinetic energy within the current window.
    ///
    /// Returns `None` if the session is empty.
    pub fn peak_kinetic_energy(&self) -> Option<f64> {
        self.window
            .iter()
            .map(|s| s.kinetic_energy)
            .reduce(f64::max)
    }

    /// Peak elapsed time (ms) within the current window.
    ///
    /// Returns `None` if the session is empty.
    pub fn peak_elapsed_ms(&self) -> Option<f64> {
        self.window.iter().map(|s| s.elapsed_ms).reduce(f64::max)
    }

    /// Mean elapsed time (ms) over the current window.
    ///
    /// Returns `None` if the session is empty.
    pub fn avg_elapsed_ms(&self) -> Option<f64> {
        self.average().map(|a| a.elapsed_ms)
    }

    /// Linear regression slope of kinetic energy over the current window
    /// (J/step).
    ///
    /// A positive slope means energy is growing; negative means dissipating.
    /// Returns `0.0` if fewer than two samples are present.
    pub fn energy_trend(&self) -> f64 {
        let n = self.window.len();
        if n < 2 {
            return 0.0;
        }
        let n_f = n as f64;
        // x = step index 0..n-1, y = kinetic_energy
        let sum_x: f64 = (0..n).map(|i| i as f64).sum();
        let sum_y: f64 = self.window.iter().map(|s| s.kinetic_energy).sum();
        let sum_xx: f64 = (0..n).map(|i| (i as f64) * (i as f64)).sum();
        let sum_xy: f64 = self
            .window
            .iter()
            .enumerate()
            .map(|(i, s)| (i as f64) * s.kinetic_energy)
            .sum();

        let denom = n_f * sum_xx - sum_x * sum_x;
        if denom.abs() < 1e-12 {
            return 0.0;
        }
        (n_f * sum_xy - sum_x * sum_y) / denom
    }

    // -----------------------------------------------------------------------
    // Export
    // -----------------------------------------------------------------------

    /// Render the current window as a CSV string with a header row.
    ///
    /// Column order: `step,dt,body_count,sleeping_count,contact_count,
    /// island_count,solve_iterations,broad_phase_pairs,kinetic_energy,elapsed_ms`
    pub fn to_csv(&self) -> String {
        let mut out = String::from(
            "step,dt,body_count,sleeping_count,contact_count,\
             island_count,solve_iterations,broad_phase_pairs,kinetic_energy,elapsed_ms\n",
        );
        for s in &self.window {
            out.push_str(&format!(
                "{},{},{},{},{},{},{},{},{},{}\n",
                s.step,
                s.dt,
                s.body_count,
                s.sleeping_count,
                s.contact_count,
                s.island_count,
                s.solve_iterations,
                s.broad_phase_pairs,
                s.kinetic_energy,
                s.elapsed_ms,
            ));
        }
        out
    }
}
