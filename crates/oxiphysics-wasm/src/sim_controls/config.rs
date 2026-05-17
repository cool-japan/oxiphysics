// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Simulation configuration types.

use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

use crate::physics_config::BroadphaseType;
use crate::wasm_helpers::{err_to_jsvalue, to_js_value};

// ---------------------------------------------------------------------------
// SimulationConfig
// ---------------------------------------------------------------------------

/// Configuration parameters for the simulation.
///
/// Controls the integration timestep, gravity vector, solver iterations,
/// sleeping threshold, and broadphase algorithm.
///
/// Exposed to JavaScript as `SimControlsConfig` to avoid colliding with the
/// engine-level `SimulationConfig` defined in [`crate::types`].
#[wasm_bindgen(js_name = "SimControlsConfig")]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SimulationConfig {
    /// Fixed timestep in seconds (e.g. 1/60).
    pub timestep: f64,
    /// Gravity vector \[x, y, z\] in m/s² (skipped from JS — use `gravity_x/y/z` getters/setters).
    #[wasm_bindgen(skip)]
    pub gravity: [f64; 3],
    /// Number of constraint solver iterations per step.
    pub iterations: u32,
    /// Linear velocity threshold below which a body may sleep (m/s).
    pub sleeping_threshold_linear: f64,
    /// Angular velocity threshold below which a body may sleep (rad/s).
    pub sleeping_threshold_angular: f64,
    /// Broadphase algorithm to use.
    pub broadphase: BroadphaseType,
    /// Maximum number of sub-steps per call to `step`.
    pub max_substeps: u32,
    /// Whether to allow bodies to sleep.
    pub sleeping_enabled: bool,
    /// CCD (continuous collision detection) enabled.
    pub ccd_enabled: bool,
}

impl Default for SimulationConfig {
    fn default() -> Self {
        SimulationConfig {
            timestep: 1.0 / 60.0,
            gravity: [0.0, -9.81, 0.0],
            iterations: 10,
            sleeping_threshold_linear: 0.01,
            sleeping_threshold_angular: 0.01,
            broadphase: BroadphaseType::Bvh,
            max_substeps: 4,
            sleeping_enabled: true,
            ccd_enabled: false,
        }
    }
}

impl SimulationConfig {
    /// Create a new `SimulationConfig` with the given timestep and gravity.
    pub fn new(timestep: f64, gravity: [f64; 3]) -> Self {
        SimulationConfig {
            timestep,
            gravity,
            ..Default::default()
        }
    }

    /// Set the number of solver iterations.
    pub fn with_iterations(mut self, iterations: u32) -> Self {
        self.iterations = iterations;
        self
    }

    /// Set sleeping thresholds (linear and angular).
    pub fn with_sleeping(mut self, linear: f64, angular: f64) -> Self {
        self.sleeping_threshold_linear = linear;
        self.sleeping_threshold_angular = angular;
        self
    }

    /// Set the broadphase algorithm.
    pub fn with_broadphase(mut self, broadphase: BroadphaseType) -> Self {
        self.broadphase = broadphase;
        self
    }

    /// Enable or disable CCD.
    pub fn with_ccd(mut self, enabled: bool) -> Self {
        self.ccd_enabled = enabled;
        self
    }

    /// Validate configuration values.
    pub fn validate(&self) -> Result<(), String> {
        if self.timestep <= 0.0 {
            return Err("timestep must be positive".to_string());
        }
        if self.iterations == 0 {
            return Err("iterations must be at least 1".to_string());
        }
        if self.max_substeps == 0 {
            return Err("max_substeps must be at least 1".to_string());
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// SimulationConfig — wasm-bindgen impl
// ---------------------------------------------------------------------------

#[wasm_bindgen(js_class = "SimControlsConfig")]
impl SimulationConfig {
    /// Construct a new config with the given timestep and gravity components.
    #[wasm_bindgen(constructor)]
    pub fn new_js(timestep: f64, gx: f64, gy: f64, gz: f64) -> SimulationConfig {
        SimulationConfig::new(timestep, [gx, gy, gz])
    }

    /// Create a config using all default values.
    #[wasm_bindgen(js_name = "default_config")]
    pub fn default_config_js() -> SimulationConfig {
        SimulationConfig::default()
    }

    /// Gravity X component (m/s²).
    #[wasm_bindgen(getter, js_name = "gravity_x")]
    pub fn gravity_x_js(&self) -> f64 {
        self.gravity[0]
    }

    /// Gravity Y component (m/s²).
    #[wasm_bindgen(getter, js_name = "gravity_y")]
    pub fn gravity_y_js(&self) -> f64 {
        self.gravity[1]
    }

    /// Gravity Z component (m/s²).
    #[wasm_bindgen(getter, js_name = "gravity_z")]
    pub fn gravity_z_js(&self) -> f64 {
        self.gravity[2]
    }

    /// Set all three gravity components at once.
    #[wasm_bindgen(js_name = "set_gravity")]
    pub fn set_gravity_js(&mut self, gx: f64, gy: f64, gz: f64) {
        self.gravity = [gx, gy, gz];
    }

    /// Return the gravity vector as a flat `[x, y, z]` array.
    #[wasm_bindgen(js_name = "gravity")]
    pub fn gravity_js(&self) -> Vec<f64> {
        self.gravity.to_vec()
    }

    /// Validate configuration values; returns an empty string on success or an
    /// error description string on failure.
    #[wasm_bindgen(js_name = "validate")]
    pub fn validate_js(&self) -> String {
        match self.validate() {
            Ok(()) => String::new(),
            Err(e) => e,
        }
    }

    /// Serialise to a JSON string for transfer to JavaScript.
    #[wasm_bindgen(js_name = "to_json")]
    pub fn to_json_js(&self) -> Result<String, JsValue> {
        serde_json::to_string(self).map_err(err_to_jsvalue)
    }

    /// Serialise to a `JsValue` object (structured-clone-friendly).
    #[wasm_bindgen(js_name = "to_js_value")]
    pub fn to_js_value_js(&self) -> Result<JsValue, JsValue> {
        to_js_value(self)
    }
}
