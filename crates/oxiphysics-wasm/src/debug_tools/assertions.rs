// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Runtime physics assertions and their results.

use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

use crate::wasm_helpers::to_js_value;

/// Result of a physics assertion check.
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AssertResult {
    /// Whether the assertion passed.
    pub passed: bool,
    /// Human-readable message.
    #[wasm_bindgen(skip)]
    pub message: String,
}

impl AssertResult {
    pub(super) fn pass(msg: impl Into<String>) -> Self {
        Self {
            passed: true,
            message: msg.into(),
        }
    }
    pub(super) fn fail(msg: impl Into<String>) -> Self {
        Self {
            passed: false,
            message: msg.into(),
        }
    }
}

#[wasm_bindgen]
impl AssertResult {
    /// Human-readable message describing the assertion.
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.message.clone()
    }

    /// Serialise to a `JsValue` object.
    #[wasm_bindgen(js_name = "to_js_value")]
    pub fn to_js_value_js(&self) -> Result<JsValue, JsValue> {
        to_js_value(self)
    }
}

/// Runtime physics assertion checker.
#[wasm_bindgen]
#[derive(Debug, Clone, Default)]
pub struct WasmPhysicsAssert {
    /// History of all assertion results.
    #[wasm_bindgen(skip)]
    pub results: Vec<AssertResult>,
}

impl WasmPhysicsAssert {
    /// Create a new asserter.
    pub fn new() -> Self {
        Self::default()
    }

    /// Assert that energy is below `max_energy`.
    pub fn assert_energy_lt(&mut self, energy: f64, max_energy: f64) -> bool {
        if energy < max_energy {
            self.results.push(AssertResult::pass(format!(
                "energy {energy:.4} < {max_energy:.4} OK"
            )));
            true
        } else {
            self.results.push(AssertResult::fail(format!(
                "energy {energy:.4} >= {max_energy:.4} FAIL"
            )));
            false
        }
    }

    /// Assert that penetration depth is below epsilon.
    pub fn assert_penetration_lt(&mut self, depth: f64, eps: f64) -> bool {
        if depth < eps {
            self.results.push(AssertResult::pass(format!(
                "depth {depth:.6} < {eps:.6} OK"
            )));
            true
        } else {
            self.results.push(AssertResult::fail(format!(
                "depth {depth:.6} >= {eps:.6} FAIL"
            )));
            false
        }
    }

    /// Assert that velocity magnitude is below `max_vel`.
    pub fn assert_velocity_lt(&mut self, vel_mag: f64, max_vel: f64) -> bool {
        if vel_mag < max_vel {
            self.results.push(AssertResult::pass(format!(
                "vel {vel_mag:.4} < {max_vel:.4} OK"
            )));
            true
        } else {
            self.results.push(AssertResult::fail(format!(
                "vel {vel_mag:.4} >= {max_vel:.4} FAIL"
            )));
            false
        }
    }

    /// Number of failing assertions.
    pub fn failure_count(&self) -> usize {
        self.results.iter().filter(|r| !r.passed).count()
    }

    /// True if all assertions passed.
    pub fn all_passed(&self) -> bool {
        self.failure_count() == 0
    }

    /// Clear results.
    pub fn clear(&mut self) {
        self.results.clear();
    }
}

#[wasm_bindgen]
impl WasmPhysicsAssert {
    /// Create a new asserter (JS-compatible constructor).
    #[wasm_bindgen(constructor)]
    pub fn wasm_new() -> WasmPhysicsAssert {
        WasmPhysicsAssert::new()
    }

    /// Assert that energy is below `max_energy`. Returns whether the
    /// assertion passed.
    #[wasm_bindgen(js_name = "assert_energy_lt")]
    pub fn assert_energy_lt_js(&mut self, energy: f64, max_energy: f64) -> bool {
        self.assert_energy_lt(energy, max_energy)
    }

    /// Assert that penetration depth is below `eps`.
    #[wasm_bindgen(js_name = "assert_penetration_lt")]
    pub fn assert_penetration_lt_js(&mut self, depth: f64, eps: f64) -> bool {
        self.assert_penetration_lt(depth, eps)
    }

    /// Assert that velocity magnitude is below `max_vel`.
    #[wasm_bindgen(js_name = "assert_velocity_lt")]
    pub fn assert_velocity_lt_js(&mut self, vel_mag: f64, max_vel: f64) -> bool {
        self.assert_velocity_lt(vel_mag, max_vel)
    }

    /// Number of failing assertions.
    #[wasm_bindgen(js_name = "failure_count")]
    pub fn failure_count_js(&self) -> usize {
        self.failure_count()
    }

    /// Total number of recorded assertions.
    #[wasm_bindgen(js_name = "total_count")]
    pub fn total_count_js(&self) -> usize {
        self.results.len()
    }

    /// True if all recorded assertions passed.
    #[wasm_bindgen(js_name = "all_passed")]
    pub fn all_passed_js(&self) -> bool {
        self.all_passed()
    }

    /// Clear all recorded assertions.
    #[wasm_bindgen(js_name = "clear")]
    pub fn clear_js(&mut self) {
        self.clear();
    }

    /// Export all recorded assertions as a `JsValue` array.
    #[wasm_bindgen(js_name = "results_to_js_value")]
    pub fn results_to_js_value_js(&self) -> Result<JsValue, JsValue> {
        to_js_value(&self.results)
    }
}
