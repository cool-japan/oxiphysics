// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Performance HUD — rolling-average frame timings and CSV export.

use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

use crate::wasm_helpers::to_js_value;

/// Per-frame timing record.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FrameTiming {
    /// Frame index.
    pub frame: u64,
    /// Frame duration (ms).
    pub ms: f64,
}

/// Rolling-average performance HUD.
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct WasmPerformanceHud {
    /// History of frame timings.
    #[wasm_bindgen(skip)]
    pub history: Vec<FrameTiming>,
    /// Capacity limit.
    pub capacity: usize,
    /// Current frame counter.
    #[wasm_bindgen(skip)]
    pub frame: u64,
}

impl WasmPerformanceHud {
    /// Create with capacity.
    pub fn new(capacity: usize) -> Self {
        Self {
            capacity,
            ..Default::default()
        }
    }

    /// Record a frame timing.
    pub fn record_frame(&mut self, ms: f64) {
        if self.capacity > 0 && self.history.len() >= self.capacity {
            self.history.remove(0);
        }
        self.history.push(FrameTiming {
            frame: self.frame,
            ms,
        });
        self.frame += 1;
    }

    /// Rolling average frame time.
    pub fn rolling_average(&self) -> f64 {
        if self.history.is_empty() {
            return 0.0;
        }
        self.history.iter().map(|f| f.ms).sum::<f64>() / self.history.len() as f64
    }

    /// Peak frame time.
    pub fn peak_ms(&self) -> f64 {
        self.history.iter().map(|f| f.ms).fold(0.0_f64, f64::max)
    }

    /// Export CSV string.
    pub fn export_csv(&self) -> String {
        let mut out = String::from("frame,ms\n");
        for f in &self.history {
            out.push_str(&format!("{},{:.4}\n", f.frame, f.ms));
        }
        out
    }
}

#[wasm_bindgen]
impl WasmPerformanceHud {
    /// Create with capacity (JS-compatible constructor).
    #[wasm_bindgen(constructor)]
    pub fn wasm_new(capacity: usize) -> WasmPerformanceHud {
        WasmPerformanceHud::new(capacity)
    }

    /// Record a frame timing of `ms` milliseconds.
    #[wasm_bindgen(js_name = "record_frame")]
    pub fn record_frame_js(&mut self, ms: f64) {
        self.record_frame(ms);
    }

    /// Rolling average frame time (ms).
    #[wasm_bindgen(js_name = "rolling_average")]
    pub fn rolling_average_js(&self) -> f64 {
        self.rolling_average()
    }

    /// Peak (maximum) recorded frame time (ms).
    #[wasm_bindgen(js_name = "peak_ms")]
    pub fn peak_ms_js(&self) -> f64 {
        self.peak_ms()
    }

    /// Number of recorded frames.
    #[wasm_bindgen(js_name = "history_len")]
    pub fn history_len_js(&self) -> usize {
        self.history.len()
    }

    /// Total frames recorded so far (cast to `f64` so JavaScript can read it
    /// without a `BigInt`).
    #[wasm_bindgen(js_name = "frame_count")]
    pub fn frame_count_js(&self) -> f64 {
        self.frame as f64
    }

    /// Export the recorded history as a CSV string.
    #[wasm_bindgen(js_name = "export_csv")]
    pub fn export_csv_js(&self) -> String {
        self.export_csv()
    }

    /// Serialise the HUD state to a `JsValue` object.
    #[wasm_bindgen(js_name = "to_js_value")]
    pub fn to_js_value_js(&self) -> Result<JsValue, JsValue> {
        to_js_value(self)
    }
}
