// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Physics performance profiling types.

use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

use crate::wasm_helpers::err_to_jsvalue;

// ---------------------------------------------------------------------------
// PhysicsProfilerSample
// ---------------------------------------------------------------------------

/// Per-step timing breakdown for the physics pipeline.
#[wasm_bindgen]
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct PhysicsProfilerSample {
    /// Time spent in broadphase collision detection (ms).
    pub broadphase_ms: f64,
    /// Time spent in narrowphase collision detection (ms).
    pub narrowphase_ms: f64,
    /// Time spent in constraint solver (ms).
    pub solver_ms: f64,
    /// Time spent in integration (ms).
    pub integration_ms: f64,
    /// Total physics time this step (ms).
    pub total_ms: f64,
    /// Step index at which this sample was taken (skipped in JS — use
    /// `step_index_f64()`).
    #[wasm_bindgen(skip)]
    pub step_index: u64,
}

impl PhysicsProfilerSample {
    /// Compute total from sub-phases.
    pub fn compute_total(&mut self) {
        self.total_ms =
            self.broadphase_ms + self.narrowphase_ms + self.solver_ms + self.integration_ms;
    }
}

#[wasm_bindgen]
impl PhysicsProfilerSample {
    /// Construct a zero-valued sample.
    #[wasm_bindgen(constructor)]
    pub fn new_js() -> PhysicsProfilerSample {
        PhysicsProfilerSample::default()
    }

    /// Step index as a JS-friendly `f64`.
    #[wasm_bindgen(js_name = "step_index_f64")]
    pub fn step_index_f64(&self) -> f64 {
        self.step_index as f64
    }

    /// Recompute `total_ms` from the per-phase breakdown.
    #[wasm_bindgen(js_name = "compute_total")]
    pub fn compute_total_js(&mut self) {
        self.compute_total();
    }

    /// Serialise to a JSON string for transfer to JavaScript.
    #[wasm_bindgen(js_name = "to_json")]
    pub fn to_json_js(&self) -> Result<String, JsValue> {
        serde_json::to_string(self).map_err(err_to_jsvalue)
    }
}

// ---------------------------------------------------------------------------
// PhysicsProfiler
// ---------------------------------------------------------------------------

/// Profiler that accumulates timing history over multiple steps.
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PhysicsProfiler {
    /// Ring buffer of per-step samples.
    #[wasm_bindgen(skip)]
    pub history: Vec<PhysicsProfilerSample>,
    /// Capacity of the history ring buffer (skipped in JS — use
    /// `history_capacity_js`).
    #[wasm_bindgen(skip)]
    pub history_capacity: usize,
    /// Running average total_ms.
    pub avg_total_ms: f64,
    /// Peak total_ms seen.
    pub peak_total_ms: f64,
}

impl PhysicsProfiler {
    /// Create a new profiler with a history ring buffer of the given size.
    pub fn new(history_capacity: usize) -> Self {
        PhysicsProfiler {
            history: Vec::with_capacity(history_capacity),
            history_capacity,
            avg_total_ms: 0.0,
            peak_total_ms: 0.0,
        }
    }

    /// Record a new sample.
    pub fn record(&mut self, mut sample: PhysicsProfilerSample) {
        sample.compute_total();
        if sample.total_ms > self.peak_total_ms {
            self.peak_total_ms = sample.total_ms;
        }
        // Update running average.
        let n = (self.history.len() + 1) as f64;
        self.avg_total_ms = (self.avg_total_ms * (n - 1.0) + sample.total_ms) / n;
        if self.history.len() >= self.history_capacity {
            self.history.remove(0);
        }
        self.history.push(sample);
    }

    /// Return the most recent sample, if any.
    pub fn latest(&self) -> Option<&PhysicsProfilerSample> {
        self.history.last()
    }

    /// Reset all history.
    pub fn reset(&mut self) {
        self.history.clear();
        self.avg_total_ms = 0.0;
        self.peak_total_ms = 0.0;
    }

    /// Return the number of recorded samples.
    pub fn sample_count(&self) -> usize {
        self.history.len()
    }
}

#[wasm_bindgen]
impl PhysicsProfiler {
    /// Construct a profiler with the given history ring-buffer capacity.
    #[wasm_bindgen(constructor)]
    pub fn new_js(history_capacity: u32) -> PhysicsProfiler {
        PhysicsProfiler::new(history_capacity as usize)
    }

    /// History ring-buffer capacity.
    #[wasm_bindgen(js_name = "history_capacity")]
    pub fn history_capacity_js(&self) -> u32 {
        self.history_capacity as u32
    }

    /// Number of samples currently stored in the history ring.
    #[wasm_bindgen(js_name = "sample_count")]
    pub fn sample_count_js(&self) -> u32 {
        self.sample_count() as u32
    }

    /// Record a new sample (the sample's `total_ms` is recomputed first).
    #[wasm_bindgen(js_name = "record")]
    pub fn record_js(&mut self, sample: &PhysicsProfilerSample) {
        self.record(sample.clone());
    }

    /// Return a clone of the most recent sample, or `None` when empty.
    #[wasm_bindgen(js_name = "latest")]
    pub fn latest_js(&self) -> Option<PhysicsProfilerSample> {
        self.latest().cloned()
    }

    /// Reset the profiler (clear history and zero peak/average).
    #[wasm_bindgen(js_name = "reset")]
    pub fn reset_js(&mut self) {
        self.reset();
    }

    /// Serialise the entire profiler state to a JSON string.
    #[wasm_bindgen(js_name = "to_json")]
    pub fn to_json_js(&self) -> Result<String, JsValue> {
        serde_json::to_string(self).map_err(err_to_jsvalue)
    }
}
