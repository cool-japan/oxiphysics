// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Test harness — expected-state comparisons and baseline diffs.

use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

use crate::wasm_helpers::{err_to_jsvalue, to_js_value};

/// Expected state for a single body.
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ExpectedBodyState {
    /// Body id.
    pub id: u32,
    /// Expected position per axis (private; use accessors).
    #[wasm_bindgen(skip)]
    pub position: [f64; 3],
    /// Position tolerance.
    pub position_tol: f64,
    /// Expected velocity (private; use accessors).
    #[wasm_bindgen(skip)]
    pub velocity: [f64; 3],
    /// Velocity tolerance.
    pub velocity_tol: f64,
}

#[wasm_bindgen]
impl ExpectedBodyState {
    /// Construct an expected-state record from individual `f64` components.
    #[wasm_bindgen(constructor)]
    pub fn wasm_new(
        id: u32,
        px: f64,
        py: f64,
        pz: f64,
        position_tol: f64,
        vx: f64,
        vy: f64,
        vz: f64,
        velocity_tol: f64,
    ) -> ExpectedBodyState {
        ExpectedBodyState {
            id,
            position: [px, py, pz],
            position_tol,
            velocity: [vx, vy, vz],
            velocity_tol,
        }
    }

    /// Expected position as `[x, y, z]`.
    #[wasm_bindgen(js_name = "get_position")]
    pub fn get_position_js(&self) -> Vec<f64> {
        self.position.to_vec()
    }

    /// Expected velocity as `[vx, vy, vz]`.
    #[wasm_bindgen(js_name = "get_velocity")]
    pub fn get_velocity_js(&self) -> Vec<f64> {
        self.velocity.to_vec()
    }

    /// Serialise to a `JsValue` object.
    #[wasm_bindgen(js_name = "to_js_value")]
    pub fn to_js_value_js(&self) -> Result<JsValue, JsValue> {
        to_js_value(self)
    }
}

/// Result of a test scenario.
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct TestResult {
    /// Whether the test passed.
    pub passed: bool,
    /// Messages from assertions.
    #[wasm_bindgen(skip)]
    pub messages: Vec<String>,
}

#[wasm_bindgen]
impl TestResult {
    /// Number of recorded messages.
    #[wasm_bindgen(js_name = "message_count")]
    pub fn message_count_js(&self) -> usize {
        self.messages.len()
    }

    /// Get the message at index `idx`, or an empty string if out of range.
    #[wasm_bindgen(js_name = "get_message")]
    pub fn get_message_js(&self, idx: usize) -> String {
        self.messages.get(idx).cloned().unwrap_or_default()
    }

    /// Concatenate all messages with newlines.
    #[wasm_bindgen(js_name = "messages_joined")]
    pub fn messages_joined_js(&self) -> String {
        self.messages.join("\n")
    }

    /// Serialise to a `JsValue` object.
    #[wasm_bindgen(js_name = "to_js_value")]
    pub fn to_js_value_js(&self) -> Result<JsValue, JsValue> {
        to_js_value(self)
    }
}

impl TestResult {
    fn pass() -> Self {
        Self {
            passed: true,
            messages: Vec::new(),
        }
    }

    fn fail(messages: Vec<String>) -> Self {
        Self {
            passed: false,
            messages,
        }
    }
}

/// Runs a physics scenario and asserts expected state.
#[wasm_bindgen]
#[derive(Debug, Clone, Default)]
pub struct WasmTestHarness {
    /// Baseline states (name -> JSON snapshot).
    baselines: Vec<(String, String)>,
    /// Last test result.
    #[wasm_bindgen(skip)]
    pub last_result: Option<TestResult>,
}

impl WasmTestHarness {
    /// Create a new harness.
    pub fn new() -> Self {
        Self::default()
    }

    /// Save a named baseline state.
    pub fn save_baseline(&mut self, name: impl Into<String>, snapshot_json: impl Into<String>) {
        let n = name.into();
        self.baselines.retain(|(bn, _)| bn != &n);
        self.baselines.push((n, snapshot_json.into()));
    }

    /// Load a baseline by name.
    pub fn load_baseline(&self, name: &str) -> Option<&str> {
        self.baselines
            .iter()
            .find_map(|(n, s)| if n == name { Some(s.as_str()) } else { None })
    }

    /// Assert that a set of actual body states matches expectations.
    /// Each actual state is `(id, position, velocity)`.
    pub fn assert_expected(
        &mut self,
        actual: &[(u32, [f64; 3], [f64; 3])],
        expected: &[ExpectedBodyState],
    ) -> bool {
        let mut messages = Vec::new();
        for exp in expected {
            if let Some(&(_id, pos, vel)) = actual.iter().find(|&&(id, _, _)| id == exp.id) {
                // Check position
                for (i, (got, want)) in pos.iter().zip(exp.position.iter()).enumerate() {
                    let diff = (got - want).abs();
                    if diff > exp.position_tol {
                        messages.push(format!(
                            "Body {} pos[{i}] diff={diff:.6} > tol={:.6}",
                            exp.id, exp.position_tol
                        ));
                    }
                }
                // Check velocity
                for (i, (got, want)) in vel.iter().zip(exp.velocity.iter()).enumerate() {
                    let diff = (got - want).abs();
                    if diff > exp.velocity_tol {
                        messages.push(format!(
                            "Body {} vel[{i}] diff={diff:.6} > tol={:.6}",
                            exp.id, exp.velocity_tol
                        ));
                    }
                }
            } else {
                messages.push(format!("Body {} not found in actual state", exp.id));
            }
        }
        let passed = messages.is_empty();
        self.last_result = Some(if passed {
            TestResult::pass()
        } else {
            TestResult::fail(messages)
        });
        passed
    }

    /// Diff current state against a saved baseline.
    pub fn diff_against_baseline(&self, name: &str, current_json: &str) -> Option<String> {
        let baseline = self.load_baseline(name)?;
        if baseline == current_json {
            Some("MATCH".to_string())
        } else {
            Some(format!(
                "DIFFER\nBaseline: {baseline}\nCurrent: {current_json}"
            ))
        }
    }
}

#[wasm_bindgen]
impl WasmTestHarness {
    /// Create a new harness (JS-compatible constructor).
    #[wasm_bindgen(constructor)]
    pub fn wasm_new() -> WasmTestHarness {
        WasmTestHarness::new()
    }

    /// Save a named baseline JSON snapshot.
    #[wasm_bindgen(js_name = "save_baseline")]
    pub fn save_baseline_js(&mut self, name: String, snapshot_json: String) {
        self.save_baseline(name, snapshot_json);
    }

    /// True when a baseline with the given name exists.
    #[wasm_bindgen(js_name = "has_baseline")]
    pub fn has_baseline_js(&self, name: String) -> bool {
        self.load_baseline(&name).is_some()
    }

    /// Load a baseline by name. Returns the JSON snapshot or an empty string
    /// if the baseline is not registered. Use [`WasmTestHarness::has_baseline_js`]
    /// to disambiguate.
    #[wasm_bindgen(js_name = "load_baseline")]
    pub fn load_baseline_js(&self, name: String) -> String {
        self.load_baseline(&name).unwrap_or("").to_string()
    }

    /// Assert that a set of actual body states matches the expected list.
    ///
    /// `actual` is a flat `Float64Array` packed as `[id, x, y, z, vx, vy, vz]`
    /// per body. `expected` is a `JsValue` array of `ExpectedBodyState`-like
    /// records (matching the JSON shape produced by `to_js_value`).
    ///
    /// # Errors
    ///
    /// Returns an error if `actual.len()` is not a multiple of 7 or if the
    /// `expected` value cannot be deserialised.
    #[wasm_bindgen(js_name = "assert_expected")]
    pub fn assert_expected_js(
        &mut self,
        actual: &[f64],
        expected: JsValue,
    ) -> Result<bool, JsValue> {
        if !actual.len().is_multiple_of(7) {
            return Err(JsValue::from_str(&format!(
                "actual length {} is not a multiple of 7 (id, x, y, z, vx, vy, vz per body)",
                actual.len()
            )));
        }
        let actual_records: Vec<(u32, [f64; 3], [f64; 3])> = actual
            .chunks_exact(7)
            .map(|c| (c[0] as u32, [c[1], c[2], c[3]], [c[4], c[5], c[6]]))
            .collect();
        let expected_records: Vec<ExpectedBodyState> =
            serde_wasm_bindgen::from_value(expected).map_err(err_to_jsvalue)?;
        Ok(self.assert_expected(&actual_records, &expected_records))
    }

    /// True when [`WasmTestHarness::assert_expected_js`] has been called.
    #[wasm_bindgen(js_name = "has_last_result")]
    pub fn has_last_result_js(&self) -> bool {
        self.last_result.is_some()
    }

    /// Get a clone of the most recent test result, or `None`/`undefined` if
    /// no test has been run yet.
    #[wasm_bindgen(js_name = "get_last_result")]
    pub fn get_last_result_js(&self) -> Option<TestResult> {
        self.last_result.clone()
    }

    /// Diff current state against a saved baseline. Returns `"MATCH"`,
    /// a `"DIFFER\n..."` description, or an empty string if the baseline does
    /// not exist.
    #[wasm_bindgen(js_name = "diff_against_baseline")]
    pub fn diff_against_baseline_js(&self, name: String, current_json: String) -> String {
        self.diff_against_baseline(&name, &current_json)
            .unwrap_or_default()
    }
}
