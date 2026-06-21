// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Simulation state, stats, and per-step result types.

use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

use crate::wasm_helpers::err_to_jsvalue;

// ---------------------------------------------------------------------------
// SimulationState
// ---------------------------------------------------------------------------

/// The current execution state of the simulation.
#[wasm_bindgen]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum SimulationState {
    /// Simulation is actively advancing every frame.
    Running,
    /// Simulation is paused; no integration is performed.
    #[default]
    Paused,
    /// Simulation advances exactly one timestep and then pauses.
    Stepping,
    /// Simulation is playing back recorded frames in reverse.
    Rewinding,
    /// Simulation is running and every frame is being recorded.
    Recording,
}

impl SimulationState {
    /// Returns `true` if integration should occur this frame.
    pub fn is_active(&self) -> bool {
        matches!(
            self,
            SimulationState::Running | SimulationState::Stepping | SimulationState::Recording
        )
    }

    /// Returns `true` if the simulation is recording frames.
    pub fn is_recording(&self) -> bool {
        *self == SimulationState::Recording
    }

    /// Human-readable label for the state.
    pub fn label(&self) -> &'static str {
        match self {
            SimulationState::Running => "running",
            SimulationState::Paused => "paused",
            SimulationState::Stepping => "stepping",
            SimulationState::Rewinding => "rewinding",
            SimulationState::Recording => "recording",
        }
    }
}

/// Free helpers exposed to JavaScript for [`SimulationState`] (enums in
/// `#[wasm_bindgen]` cannot have inherent methods taking `&self`).
#[wasm_bindgen]
pub fn sim_state_is_active(state: SimulationState) -> bool {
    state.is_active()
}

/// JS-friendly wrapper for [`SimulationState::is_recording`].
#[wasm_bindgen]
pub fn sim_state_is_recording(state: SimulationState) -> bool {
    state.is_recording()
}

/// JS-friendly wrapper returning the [`SimulationState`] human-readable label.
#[wasm_bindgen]
pub fn sim_state_label(state: SimulationState) -> String {
    state.label().to_string()
}

// ---------------------------------------------------------------------------
// SimulationStats
// ---------------------------------------------------------------------------

/// Aggregate statistics collected over the lifetime of the simulation.
#[wasm_bindgen]
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct SimulationStats {
    /// Total number of integration steps completed (skipped in JS — use
    /// `step_count_f64()` to obtain it as a JavaScript `Number`).
    #[wasm_bindgen(skip)]
    pub step_count: u64,
    /// Total simulated time in seconds.
    pub elapsed_time: f64,
    /// Current number of rigid bodies in the scene.
    pub body_count: u32,
    /// Current number of active constraints.
    pub constraint_count: u32,
    /// Number of simulation islands (disconnected groups) this step.
    pub island_count: u32,
    /// Number of active collision pairs detected this step.
    pub collision_pairs: u32,
    /// Peak memory used by the physics system in bytes (skipped in JS — use
    /// `peak_memory_bytes_f64()` to obtain it as a JavaScript `Number`).
    #[wasm_bindgen(skip)]
    pub peak_memory_bytes: u64,
    /// Number of bodies currently sleeping.
    pub sleeping_bodies: u32,
    /// Number of CCD sub-steps taken this step.
    pub ccd_substeps: u32,
}

impl SimulationStats {
    /// Create a zeroed `SimulationStats`.
    pub fn new() -> Self {
        Self::default()
    }

    /// Advance step count and elapsed time.
    pub fn advance(&mut self, dt: f64) {
        self.step_count += 1;
        self.elapsed_time += dt;
    }

    /// Reset all counters while keeping lifetime step_count/elapsed_time.
    pub fn reset_frame_counters(&mut self) {
        self.collision_pairs = 0;
        self.island_count = 0;
        self.ccd_substeps = 0;
    }
}

#[wasm_bindgen]
impl SimulationStats {
    /// Construct a zeroed [`SimulationStats`] from JavaScript.
    #[wasm_bindgen(constructor)]
    pub fn new_js() -> SimulationStats {
        SimulationStats::default()
    }

    /// Total number of integration steps completed, as a JS-friendly `f64`.
    #[wasm_bindgen(js_name = "step_count_f64")]
    pub fn step_count_f64(&self) -> f64 {
        self.step_count as f64
    }

    /// Peak memory used (bytes), as a JS-friendly `f64`.
    #[wasm_bindgen(js_name = "peak_memory_bytes_f64")]
    pub fn peak_memory_bytes_f64(&self) -> f64 {
        self.peak_memory_bytes as f64
    }

    /// Advance step count and elapsed time by one tick of size `dt`.
    #[wasm_bindgen(js_name = "advance")]
    pub fn advance_js(&mut self, dt: f64) {
        self.advance(dt);
    }

    /// Reset frame counters but keep lifetime totals.
    #[wasm_bindgen(js_name = "reset_frame_counters")]
    pub fn reset_frame_counters_js(&mut self) {
        self.reset_frame_counters();
    }

    /// Serialise to a JSON string for transfer to JavaScript.
    #[wasm_bindgen(js_name = "to_json")]
    pub fn to_json_js(&self) -> Result<String, JsValue> {
        serde_json::to_string(self).map_err(err_to_jsvalue)
    }
}

// ---------------------------------------------------------------------------
// StepResult
// ---------------------------------------------------------------------------

/// Per-step diagnostic result returned by `SimulationController::step`.
#[wasm_bindgen]
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct StepResult {
    /// Number of new contact manifolds generated this step.
    pub new_contacts: u32,
    /// Number of contact manifolds resolved (separated) this step.
    pub resolved_contacts: u32,
    /// Number of bodies that transitioned to sleep.
    pub sleeping_bodies: u32,
    /// Number of active simulation islands this step.
    pub active_islands: u32,
    /// Wall-clock time spent in physics this step (ms).
    ///
    /// Measured with a real monotonic clock (`std::time::Instant` natively,
    /// `performance.now()` on wasm). `0.0` when the step did nothing (paused /
    /// no substeps). `f64::NAN` is an honest sentinel meaning the target
    /// exposes no clock — measure on the JS side via `performance.now()`.
    pub perf_ms: f64,
    /// Whether any constraint was violated beyond the tolerance.
    pub constraint_violation: bool,
}

impl StepResult {
    /// Create a zeroed `StepResult`.
    pub fn new() -> Self {
        Self::default()
    }

    /// Return `true` if any performance-relevant event occurred.
    pub fn has_events(&self) -> bool {
        self.new_contacts > 0 || self.sleeping_bodies > 0 || self.constraint_violation
    }
}

#[wasm_bindgen]
impl StepResult {
    /// Construct a zeroed [`StepResult`] from JavaScript.
    #[wasm_bindgen(constructor)]
    pub fn new_js() -> StepResult {
        StepResult::default()
    }

    /// JS-friendly wrapper for [`StepResult::has_events`].
    #[wasm_bindgen(js_name = "has_events")]
    pub fn has_events_js(&self) -> bool {
        self.has_events()
    }

    /// Serialise to a JSON string for transfer to JavaScript.
    #[wasm_bindgen(js_name = "to_json")]
    pub fn to_json_js(&self) -> Result<String, JsValue> {
        serde_json::to_string(self).map_err(err_to_jsvalue)
    }
}
