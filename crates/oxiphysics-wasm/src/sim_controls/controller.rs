// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! SimulationController — drives the main simulation loop.

use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

use super::config::SimulationConfig;
use super::state::{SimulationState, SimulationStats, StepResult};

// ---------------------------------------------------------------------------
// SimulationController
// ---------------------------------------------------------------------------

/// High-level controller that owns a `SimulationConfig` and drives the loop.
///
/// The controller is annotated with `#[wasm_bindgen]` so it can be created and
/// driven directly from JavaScript. Fields holding non-primitive types are
/// marked `#[wasm_bindgen(skip)]`; JS-friendly getters/setters are provided
/// in the dedicated [`#[wasm_bindgen] impl`](#impl-SimulationController-1)
/// block below.
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SimulationController {
    /// Active configuration.
    #[wasm_bindgen(skip)]
    pub config: SimulationConfig,
    /// Current execution state.
    pub state: SimulationState,
    /// Accumulated statistics.
    #[wasm_bindgen(skip)]
    pub stats: SimulationStats,
    /// Accumulated time not yet consumed by fixed steps.
    pub accumulator: f64,
}

impl SimulationController {
    /// Create a new controller with the given config. Initial state is `Paused`.
    pub fn new(config: SimulationConfig) -> Self {
        SimulationController {
            config,
            state: SimulationState::Paused,
            stats: SimulationStats::new(),
            accumulator: 0.0,
        }
    }

    /// Advance the simulation by `frame_time` seconds (wall-clock delta).
    ///
    /// Returns a `StepResult` summarising this frame's activity.
    pub fn step(&mut self, frame_time: f64) -> StepResult {
        let mut result = StepResult::new();

        if !self.state.is_active() {
            return result;
        }

        self.accumulator += frame_time;
        let dt = self.config.timestep;
        let mut substeps = 0u32;

        while self.accumulator >= dt && substeps < self.config.max_substeps {
            self.stats.advance(dt);
            self.accumulator -= dt;
            substeps += 1;
            result.active_islands = self.stats.island_count;
        }

        // If we were in Stepping mode, pause after one logical frame.
        if self.state == SimulationState::Stepping {
            self.state = SimulationState::Paused;
        }

        result.perf_ms = frame_time * 1_000.0 * 0.01; // mock
        result
    }

    /// Pause the simulation.
    pub fn pause(&mut self) {
        self.state = SimulationState::Paused;
    }

    /// Resume from paused.
    pub fn resume(&mut self) {
        self.state = SimulationState::Running;
    }

    /// Advance exactly one fixed-timestep and pause.
    pub fn step_once(&mut self) {
        self.state = SimulationState::Stepping;
        self.accumulator += self.config.timestep;
    }

    /// Begin recording frames.
    pub fn start_recording(&mut self) {
        self.state = SimulationState::Recording;
    }

    /// Begin rewinding recorded frames.
    pub fn start_rewind(&mut self) {
        self.state = SimulationState::Rewinding;
    }

    /// Reset simulation to t=0.
    pub fn reset(&mut self) {
        self.stats = SimulationStats::new();
        self.accumulator = 0.0;
        self.state = SimulationState::Paused;
    }

    /// Get a reference to the current statistics.
    pub fn get_stats(&self) -> &SimulationStats {
        &self.stats
    }

    /// Update body and constraint counts in stats.
    pub fn update_counts(&mut self, bodies: u32, constraints: u32) {
        self.stats.body_count = bodies;
        self.stats.constraint_count = constraints;
    }
}

// ---------------------------------------------------------------------------
// SimulationController — wasm-bindgen impl
// ---------------------------------------------------------------------------

#[wasm_bindgen]
impl SimulationController {
    /// Construct a controller with default configuration values.
    #[wasm_bindgen(constructor)]
    pub fn new_js() -> SimulationController {
        SimulationController::new(SimulationConfig::default())
    }

    /// Replace the active configuration with a clone of `config`.
    #[wasm_bindgen(js_name = "set_config")]
    pub fn set_config_js(&mut self, config: &SimulationConfig) {
        self.config = config.clone();
    }

    /// Return a clone of the active configuration.
    #[wasm_bindgen(js_name = "get_config")]
    pub fn get_config_js(&self) -> SimulationConfig {
        self.config.clone()
    }

    /// Return a clone of the current statistics.
    #[wasm_bindgen(js_name = "get_stats")]
    pub fn get_stats_js(&self) -> SimulationStats {
        self.stats.clone()
    }

    /// Advance the simulation by `frame_time` seconds and return a `StepResult`.
    #[wasm_bindgen(js_name = "step")]
    pub fn step_js(&mut self, frame_time: f64) -> StepResult {
        self.step(frame_time)
    }

    /// Pause the simulation.
    #[wasm_bindgen(js_name = "pause")]
    pub fn pause_js(&mut self) {
        self.pause();
    }

    /// Resume from paused.
    #[wasm_bindgen(js_name = "resume")]
    pub fn resume_js(&mut self) {
        self.resume();
    }

    /// Advance exactly one fixed-timestep and pause afterwards.
    #[wasm_bindgen(js_name = "step_once")]
    pub fn step_once_js(&mut self) {
        self.step_once();
    }

    /// Begin recording frames.
    #[wasm_bindgen(js_name = "start_recording")]
    pub fn start_recording_js(&mut self) {
        self.start_recording();
    }

    /// Begin rewinding recorded frames.
    #[wasm_bindgen(js_name = "start_rewind")]
    pub fn start_rewind_js(&mut self) {
        self.start_rewind();
    }

    /// Reset the simulation to t=0.
    #[wasm_bindgen(js_name = "reset")]
    pub fn reset_js(&mut self) {
        self.reset();
    }

    /// Update body and constraint counts in the stats.
    #[wasm_bindgen(js_name = "update_counts")]
    pub fn update_counts_js(&mut self, bodies: u32, constraints: u32) {
        self.update_counts(bodies, constraints);
    }
}
