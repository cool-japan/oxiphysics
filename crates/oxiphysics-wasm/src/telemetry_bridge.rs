// Copyright (C) 2026 COOLJAPAN OU (Team Kitasan)
// SPDX-License-Identifier: Apache-2.0

//! WebAssembly bridge for physics simulation telemetry.
//!
//! Wraps per-step physics statistics and a rolling window session with
//! a JSON-oriented surface suitable for use across the WASM boundary.

use serde::{Deserialize, Serialize};
use std::collections::VecDeque;
use wasm_bindgen::prelude::*;

use crate::wasm_helpers::to_js_value;

// ---------------------------------------------------------------------------
// WasmPhysicsStats
// ---------------------------------------------------------------------------

/// Per-step physics simulation metrics.
///
/// `usize` count fields and `u64` step are private; JavaScript reads them
/// through `u32` / `f64` accessors on the `#[wasm_bindgen]` impl block.
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmPhysicsStats {
    /// Simulation step index — private; JS reads as `f64` via `step()`.
    pub(crate) step: u64,
    /// Physics time-step duration (seconds).
    pub dt: f64,
    /// Total active bodies — private; JS reads as `u32`.
    pub(crate) body_count: usize,
    /// Sleeping bodies — private; JS reads as `u32`.
    pub(crate) sleeping_count: usize,
    /// Active contact pairs — private; JS reads as `u32`.
    pub(crate) contact_count: usize,
    /// Simulation islands — private; JS reads as `u32`.
    pub(crate) island_count: usize,
    /// Constraint solver iterations — private; JS reads as `u32`.
    pub(crate) solve_iterations: usize,
    /// Broad-phase candidate pairs — private; JS reads as `u32`.
    pub(crate) broad_phase_pairs: usize,
    /// Total kinetic energy of all bodies (J).
    pub kinetic_energy: f64,
    /// Wall-clock time consumed by this step (milliseconds).
    pub elapsed_ms: f64,
}

impl WasmPhysicsStats {
    /// Create a zeroed stats record for `step` and `dt`.
    pub fn new(step: u64, dt: f64) -> Self {
        WasmPhysicsStats {
            step,
            dt,
            body_count: 0,
            sleeping_count: 0,
            contact_count: 0,
            island_count: 0,
            solve_iterations: 0,
            broad_phase_pairs: 0,
            kinetic_energy: 0.0,
            elapsed_ms: 0.0,
        }
    }

    /// Serialise this record as a JSON string.
    pub fn to_json(&self) -> String {
        serde_json::to_string(self).unwrap_or_default()
    }
}

#[wasm_bindgen]
impl WasmPhysicsStats {
    /// Create a zeroed stats record (JS constructor).
    ///
    /// `step` is accepted as `f64` because `u64` is unsupported by
    /// wasm-bindgen.
    #[wasm_bindgen(constructor)]
    pub fn new_js(step: f64, dt: f64) -> WasmPhysicsStats {
        WasmPhysicsStats::new(step as u64, dt)
    }

    /// Simulation step index as `f64` (`u64` unsupported by wasm-bindgen).
    #[wasm_bindgen(js_name = "step")]
    pub fn step_js(&self) -> f64 {
        self.step as f64
    }

    /// Set simulation step index.
    #[wasm_bindgen(js_name = "set_step")]
    pub fn set_step_js(&mut self, step: f64) {
        self.step = step as u64;
    }

    /// Active body count as `u32`.
    #[wasm_bindgen(js_name = "body_count")]
    pub fn body_count_js(&self) -> u32 {
        self.body_count as u32
    }

    /// Set active body count.
    #[wasm_bindgen(js_name = "set_body_count")]
    pub fn set_body_count_js(&mut self, count: u32) {
        self.body_count = count as usize;
    }

    /// Sleeping body count as `u32`.
    #[wasm_bindgen(js_name = "sleeping_count")]
    pub fn sleeping_count_js(&self) -> u32 {
        self.sleeping_count as u32
    }

    /// Set sleeping body count.
    #[wasm_bindgen(js_name = "set_sleeping_count")]
    pub fn set_sleeping_count_js(&mut self, count: u32) {
        self.sleeping_count = count as usize;
    }

    /// Active contact pair count as `u32`.
    #[wasm_bindgen(js_name = "contact_count")]
    pub fn contact_count_js(&self) -> u32 {
        self.contact_count as u32
    }

    /// Set contact count.
    #[wasm_bindgen(js_name = "set_contact_count")]
    pub fn set_contact_count_js(&mut self, count: u32) {
        self.contact_count = count as usize;
    }

    /// Simulation island count as `u32`.
    #[wasm_bindgen(js_name = "island_count")]
    pub fn island_count_js(&self) -> u32 {
        self.island_count as u32
    }

    /// Set island count.
    #[wasm_bindgen(js_name = "set_island_count")]
    pub fn set_island_count_js(&mut self, count: u32) {
        self.island_count = count as usize;
    }

    /// Solver iterations performed this step as `u32`.
    #[wasm_bindgen(js_name = "solve_iterations")]
    pub fn solve_iterations_js(&self) -> u32 {
        self.solve_iterations as u32
    }

    /// Set solver iterations.
    #[wasm_bindgen(js_name = "set_solve_iterations")]
    pub fn set_solve_iterations_js(&mut self, n: u32) {
        self.solve_iterations = n as usize;
    }

    /// Broad-phase candidate pair count as `u32`.
    #[wasm_bindgen(js_name = "broad_phase_pairs")]
    pub fn broad_phase_pairs_js(&self) -> u32 {
        self.broad_phase_pairs as u32
    }

    /// Set broad-phase pair count.
    #[wasm_bindgen(js_name = "set_broad_phase_pairs")]
    pub fn set_broad_phase_pairs_js(&mut self, n: u32) {
        self.broad_phase_pairs = n as usize;
    }

    /// Serialise this record as a JSON string.
    #[wasm_bindgen(js_name = "to_json")]
    pub fn to_json_js(&self) -> String {
        self.to_json()
    }
}

// ---------------------------------------------------------------------------
// WasmPhysicsAverages
// ---------------------------------------------------------------------------

/// Windowed averages of continuous per-step metrics.
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmPhysicsAverages {
    /// Average physics time step (seconds).
    pub dt: f64,
    /// Average kinetic energy (J).
    pub kinetic_energy: f64,
    /// Average wall-clock elapsed time (milliseconds).
    pub elapsed_ms: f64,
    /// Number of samples averaged — private; JS reads as `u32`.
    pub(crate) sample_count: usize,
}

#[wasm_bindgen]
impl WasmPhysicsAverages {
    /// Number of samples averaged.
    #[wasm_bindgen(js_name = "sample_count")]
    pub fn sample_count_js(&self) -> u32 {
        self.sample_count as u32
    }
}

// ---------------------------------------------------------------------------
// WasmTelemetrySession
// ---------------------------------------------------------------------------

/// WASM wrapper for a rolling-window telemetry session.
///
/// Maintains the last `window_size` [`WasmPhysicsStats`] entries.
/// When the window is full the oldest entry is dropped automatically.
///
/// `window_size: usize` and `total_steps: u64` are private; JavaScript reads
/// them via `u32` / `f64` accessors on the `#[wasm_bindgen]` impl block.
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmTelemetrySession {
    window: VecDeque<WasmPhysicsStats>,
    /// Maximum retained samples — private.
    pub(crate) window_size: usize,
    /// Total steps pushed since creation — private.
    pub(crate) total_steps: u64,
}

impl WasmTelemetrySession {
    /// Create a new session with the given rolling-window size (must be > 0).
    pub fn new(window_size: usize) -> Self {
        let window_size = window_size.max(1);
        WasmTelemetrySession {
            window: VecDeque::with_capacity(window_size),
            window_size,
            total_steps: 0,
        }
    }

    /// Append a stats record.
    pub fn push(&mut self, stats: WasmPhysicsStats) {
        if self.window.len() == self.window_size {
            self.window.pop_front();
        }
        self.window.push_back(stats);
        self.total_steps += 1;
    }

    /// Clear all entries.
    pub fn clear(&mut self) {
        self.window.clear();
    }

    /// Most recent stats record, or `None` if empty.
    pub fn latest(&self) -> Option<&WasmPhysicsStats> {
        self.window.back()
    }

    /// Number of samples currently in the window.
    pub fn len(&self) -> usize {
        self.window.len()
    }

    /// Returns `true` if no samples have been pushed.
    pub fn is_empty(&self) -> bool {
        self.window.is_empty()
    }

    /// Windowed averages over the current window, or `None` if empty.
    pub fn average(&self) -> Option<WasmPhysicsAverages> {
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
        Some(WasmPhysicsAverages {
            dt: sum_dt / n_f,
            kinetic_energy: sum_ke / n_f,
            elapsed_ms: sum_ms / n_f,
            sample_count: n,
        })
    }

    /// Peak kinetic energy within the current window.
    pub fn peak_kinetic_energy(&self) -> Option<f64> {
        self.window
            .iter()
            .map(|s| s.kinetic_energy)
            .reduce(f64::max)
    }

    /// Peak elapsed time (ms) within the current window.
    pub fn peak_elapsed_ms(&self) -> Option<f64> {
        self.window.iter().map(|s| s.elapsed_ms).reduce(f64::max)
    }

    /// Linear regression slope of kinetic energy over the current window (J/step).
    pub fn energy_trend(&self) -> f64 {
        let n = self.window.len();
        if n < 2 {
            return 0.0;
        }
        let n_f = n as f64;
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

    /// Render the current window as a CSV string with a header row.
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

    /// Serialise the session as a JSON string.
    pub fn to_json(&self) -> String {
        serde_json::to_string(self).unwrap_or_default()
    }
}

// ---------------------------------------------------------------------------
// wasm-bindgen JavaScript API
// ---------------------------------------------------------------------------

#[wasm_bindgen]
impl WasmTelemetrySession {
    /// Create a new session (JS constructor).
    #[wasm_bindgen(constructor)]
    pub fn new_js(window_size: u32) -> WasmTelemetrySession {
        WasmTelemetrySession::new(window_size as usize)
    }

    /// Append a stats record (consumes the JS-side handle).
    #[wasm_bindgen(js_name = "push")]
    pub fn push_js(&mut self, stats: WasmPhysicsStats) {
        self.push(stats);
    }

    /// Append a stats record from a `WasmPhysicsStats` object.
    ///
    /// Aliased as `push_record` on the JS side; prefer constructing a
    /// `WasmPhysicsStats` directly rather than passing ten scalar arguments.
    #[wasm_bindgen(js_name = "push_record")]
    pub fn push_record_js(&mut self, stats: WasmPhysicsStats) {
        self.push(stats);
    }

    /// Drop all retained samples.
    #[wasm_bindgen(js_name = "clear")]
    pub fn clear_js(&mut self) {
        self.clear();
    }

    /// Return the most recent stats record as a fresh `WasmPhysicsStats`,
    /// or `None` when empty.
    #[wasm_bindgen(js_name = "latest")]
    pub fn latest_js(&self) -> Option<WasmPhysicsStats> {
        self.latest().cloned()
    }

    /// Number of samples currently retained.
    #[wasm_bindgen(js_name = "len")]
    pub fn len_js(&self) -> u32 {
        self.len() as u32
    }

    /// `true` when no samples have been pushed.
    #[wasm_bindgen(js_name = "is_empty")]
    pub fn is_empty_js(&self) -> bool {
        self.is_empty()
    }

    /// Maximum retained samples.
    #[wasm_bindgen(js_name = "window_size")]
    pub fn window_size_js(&self) -> u32 {
        self.window_size as u32
    }

    /// Total samples pushed since creation as `f64` (`u64` unsupported).
    #[wasm_bindgen(js_name = "total_steps")]
    pub fn total_steps_js(&self) -> f64 {
        self.total_steps as f64
    }

    /// Windowed averages, or `None` when the window is empty.
    #[wasm_bindgen(js_name = "average")]
    pub fn average_js(&self) -> Option<WasmPhysicsAverages> {
        self.average()
    }

    /// Peak kinetic energy across the window, or `None` when empty.
    #[wasm_bindgen(js_name = "peak_kinetic_energy")]
    pub fn peak_kinetic_energy_js(&self) -> Option<f64> {
        self.peak_kinetic_energy()
    }

    /// Peak wall-clock elapsed time across the window, or `None` when empty.
    #[wasm_bindgen(js_name = "peak_elapsed_ms")]
    pub fn peak_elapsed_ms_js(&self) -> Option<f64> {
        self.peak_elapsed_ms()
    }

    /// Linear regression slope of kinetic energy across the window (J/step).
    #[wasm_bindgen(js_name = "energy_trend")]
    pub fn energy_trend_js(&self) -> f64 {
        self.energy_trend()
    }

    /// Render the current window as CSV (header + data rows).
    #[wasm_bindgen(js_name = "to_csv")]
    pub fn to_csv_js(&self) -> String {
        self.to_csv()
    }

    /// Serialise the session as a `JsValue` (plain JS object).
    ///
    /// # Errors
    ///
    /// Returns a `JsValue` error string if serialisation fails.
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

    fn sample(step: u64, ke: f64) -> WasmPhysicsStats {
        WasmPhysicsStats {
            step,
            dt: 0.016,
            body_count: 10,
            sleeping_count: 2,
            contact_count: 5,
            island_count: 1,
            solve_iterations: 8,
            broad_phase_pairs: 20,
            kinetic_energy: ke,
            elapsed_ms: 1.5,
        }
    }

    #[test]
    fn test_telemetry_bridge_instantiation() {
        let session = WasmTelemetrySession::new(60);
        let json = session.to_json();
        assert!(!json.is_empty(), "to_json should return non-empty string");
    }

    #[test]
    fn test_telemetry_bridge_push_and_average() {
        let mut session = WasmTelemetrySession::new(60);
        session.push(sample(1, 42.0));
        session.push(sample(2, 58.0));
        let avg = session.average().expect("should have average");
        assert!((avg.kinetic_energy - 50.0).abs() < 1e-9);
        assert_eq!(avg.sample_count, 2);
    }

    #[test]
    fn test_telemetry_bridge_window_eviction() {
        let mut session = WasmTelemetrySession::new(3);
        for i in 0..10 {
            session.push(sample(i, i as f64));
        }
        assert_eq!(session.len(), 3);
        assert_eq!(session.latest().expect("has latest").step, 9);
    }

    #[test]
    fn test_telemetry_bridge_peak_ke() {
        let mut session = WasmTelemetrySession::new(10);
        session.push(sample(0, 10.0));
        session.push(sample(1, 99.0));
        session.push(sample(2, 5.0));
        assert!((session.peak_kinetic_energy().expect("has peak") - 99.0).abs() < 1e-9);
    }

    #[test]
    fn test_telemetry_bridge_csv_format() {
        let mut session = WasmTelemetrySession::new(60);
        session.push(sample(0, 1.0));
        let csv = session.to_csv();
        let lines: Vec<&str> = csv.lines().collect();
        assert!(lines.len() >= 2, "CSV should have header + at least 1 row");
        let header = lines[0];
        assert!(header.contains("step"), "header should contain 'step'");
        assert!(
            header.contains("kinetic_energy"),
            "header should contain 'kinetic_energy'"
        );
    }

    #[test]
    fn test_telemetry_bridge_energy_trend() {
        let mut session = WasmTelemetrySession::new(60);
        for i in 0..5 {
            session.push(sample(i, i as f64 * 2.0)); // increasing
        }
        let trend = session.energy_trend();
        assert!(
            trend > 0.0,
            "energy trend should be positive for increasing KE"
        );
    }
}
