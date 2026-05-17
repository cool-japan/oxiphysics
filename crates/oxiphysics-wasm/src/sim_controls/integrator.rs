// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Time integration schemes and configuration.

use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

use crate::wasm_helpers::err_to_jsvalue;

// ---------------------------------------------------------------------------
// IntegratorType
// ---------------------------------------------------------------------------

/// Integration scheme for advancing body state.
#[wasm_bindgen]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum IntegratorType {
    /// Forward (explicit) Euler.
    Euler,
    /// Semi-implicit Euler (symplectic Euler).
    #[default]
    SemiImplicit,
    /// Velocity Verlet.
    Verlet,
    /// Fourth-order Runge-Kutta.
    Rk4,
}

impl IntegratorType {
    /// Human-readable name.
    pub fn name(&self) -> &'static str {
        match self {
            IntegratorType::Euler => "euler",
            IntegratorType::SemiImplicit => "semi_implicit",
            IntegratorType::Verlet => "verlet",
            IntegratorType::Rk4 => "rk4",
        }
    }

    /// Whether the integrator is symplectic (energy-conserving in limit).
    pub fn is_symplectic(&self) -> bool {
        matches!(self, IntegratorType::SemiImplicit | IntegratorType::Verlet)
    }
}

/// Return the human-readable scheme name for an [`IntegratorType`].
///
/// Free helpers are required because enums in `#[wasm_bindgen]` cannot have
/// inherent methods taking `&self`.
#[wasm_bindgen]
pub fn integrator_name(t: IntegratorType) -> String {
    t.name().to_string()
}

/// JS-friendly wrapper for [`IntegratorType::is_symplectic`].
#[wasm_bindgen]
pub fn integrator_is_symplectic(t: IntegratorType) -> bool {
    t.is_symplectic()
}

// ---------------------------------------------------------------------------
// TimeIntegrator
// ---------------------------------------------------------------------------

/// Configuration for the time integrator.
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TimeIntegrator {
    /// Selected integration scheme.
    pub integrator_type: IntegratorType,
    /// Number of sub-steps per timestep.
    pub sub_steps: u32,
    /// Damping coefficient applied during integration.
    pub damping: f64,
    /// Angular damping coefficient.
    pub angular_damping: f64,
}

impl Default for TimeIntegrator {
    fn default() -> Self {
        TimeIntegrator {
            integrator_type: IntegratorType::SemiImplicit,
            sub_steps: 1,
            damping: 0.0,
            angular_damping: 0.0,
        }
    }
}

impl TimeIntegrator {
    /// Create with a specific scheme.
    pub fn new(integrator_type: IntegratorType) -> Self {
        TimeIntegrator {
            integrator_type,
            ..Default::default()
        }
    }

    /// Set sub-step count.
    pub fn with_substeps(mut self, n: u32) -> Self {
        self.sub_steps = n.max(1);
        self
    }

    /// Set linear damping.
    pub fn with_damping(mut self, damping: f64) -> Self {
        self.damping = damping.clamp(0.0, 1.0);
        self
    }

    /// Integrate a 1D position/velocity pair using the configured scheme.
    ///
    /// Returns `(new_position, new_velocity)`.
    pub fn integrate_1d(&self, pos: f64, vel: f64, acc: f64, dt: f64) -> (f64, f64) {
        let sub_dt = dt / self.sub_steps as f64;
        let mut p = pos;
        let mut v = vel;
        for _ in 0..self.sub_steps {
            (p, v) = match self.integrator_type {
                IntegratorType::Euler => (p + v * sub_dt, v + acc * sub_dt),
                IntegratorType::SemiImplicit => {
                    let v2 = v + acc * sub_dt;
                    (p + v2 * sub_dt, v2)
                }
                IntegratorType::Verlet => {
                    let p2 = p + v * sub_dt + 0.5 * acc * sub_dt * sub_dt;
                    let v2 = v + acc * sub_dt;
                    (p2, v2)
                }
                IntegratorType::Rk4 => {
                    // Simplified RK4 for constant acceleration.
                    let k1v = acc;
                    let k2v = acc;
                    let k3v = acc;
                    let k4v = acc;
                    let k1p = v;
                    let k2p = v + 0.5 * sub_dt * k1v;
                    let k3p = v + 0.5 * sub_dt * k2v;
                    let k4p = v + sub_dt * k3v;
                    let v2 = v + (sub_dt / 6.0) * (k1v + 2.0 * k2v + 2.0 * k3v + k4v);
                    let p2 = p + (sub_dt / 6.0) * (k1p + 2.0 * k2p + 2.0 * k3p + k4p);
                    (p2, v2)
                }
            };
            v *= 1.0 - self.damping * sub_dt;
        }
        (p, v)
    }
}

#[wasm_bindgen]
impl TimeIntegrator {
    /// Construct an integrator using the given scheme.
    #[wasm_bindgen(constructor)]
    pub fn new_js(integrator_type: IntegratorType) -> TimeIntegrator {
        TimeIntegrator::new(integrator_type)
    }

    /// Construct an integrator from default values (semi-implicit Euler).
    #[wasm_bindgen(js_name = "default_integrator")]
    pub fn default_integrator_js() -> TimeIntegrator {
        TimeIntegrator::default()
    }

    /// Set the sub-step count (clamped to a minimum of 1).
    #[wasm_bindgen(js_name = "set_substeps")]
    pub fn set_substeps_js(&mut self, n: u32) {
        self.sub_steps = n.max(1);
    }

    /// Set the linear damping coefficient (clamped to `[0, 1]`).
    #[wasm_bindgen(js_name = "set_damping")]
    pub fn set_damping_js(&mut self, damping: f64) {
        self.damping = damping.clamp(0.0, 1.0);
    }

    /// Set the angular damping coefficient.
    #[wasm_bindgen(js_name = "set_angular_damping")]
    pub fn set_angular_damping_js(&mut self, angular_damping: f64) {
        self.angular_damping = angular_damping;
    }

    /// Integrate a 1D position/velocity pair using the configured scheme.
    ///
    /// Returns a flat `[new_position, new_velocity]` `Float64Array` for use
    /// from JavaScript.
    #[wasm_bindgen(js_name = "integrate_1d")]
    pub fn integrate_1d_js(&self, pos: f64, vel: f64, acc: f64, dt: f64) -> Vec<f64> {
        let (p, v) = self.integrate_1d(pos, vel, acc, dt);
        vec![p, v]
    }

    /// Serialise to a JSON string.
    #[wasm_bindgen(js_name = "to_json")]
    pub fn to_json_js(&self) -> Result<String, JsValue> {
        serde_json::to_string(self).map_err(err_to_jsvalue)
    }
}
