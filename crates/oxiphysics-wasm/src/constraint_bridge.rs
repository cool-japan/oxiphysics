// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! WebAssembly constraint system bridge.
//!
//! Provides Rust types for joints, contacts, motors, the constraint solver,
//! islands, CCD, and ragdoll construction, fully exposed to JavaScript via
//! wasm-bindgen.

use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

use crate::wasm_helpers::{err_to_jsvalue, to_js_value};

// ---------------------------------------------------------------------------
// WasmJointHandle (re-exported from simulation_api)
// ---------------------------------------------------------------------------

pub use crate::simulation_api::WasmJointHandle;

// ---------------------------------------------------------------------------
// WasmJointConfig
// ---------------------------------------------------------------------------

/// Configuration for a joint connecting two bodies.
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmJointConfig {
    /// Joint type: `"revolute"`, `"prismatic"`, `"ball"`, `"fixed"`, or `"spring"`.
    #[wasm_bindgen(skip)]
    pub joint_type: String,
    /// Identifier of body A (index or handle).
    pub body_a: u64,
    /// Identifier of body B (index or handle).
    pub body_b: u64,
    /// Anchor point on body A in local space \[x, y, z\].
    #[wasm_bindgen(skip)]
    pub anchor_a: [f64; 3],
    /// Anchor point on body B in local space \[x, y, z\].
    #[wasm_bindgen(skip)]
    pub anchor_b: [f64; 3],
    /// Joint axis in body A local space \[x, y, z\] (used by revolute/prismatic).
    #[wasm_bindgen(skip)]
    pub axis_a: [f64; 3],
    /// Joint axis in body B local space \[x, y, z\].
    #[wasm_bindgen(skip)]
    pub axis_b: [f64; 3],
    /// Lower limit (angle in radians for revolute, distance in metres for prismatic).
    pub lower_limit: f64,
    /// Upper limit.
    pub upper_limit: f64,
    /// Whether limits are enabled.
    pub limits_enabled: bool,
    /// Motor enabled flag.
    pub motor_enabled: bool,
    /// Maximum motor force/torque (N or N·m).
    pub motor_max_force: f64,
    /// Motor target velocity (rad/s or m/s).
    pub motor_target_velocity: f64,
    /// Spring stiffness (N/m) – used when `joint_type == "spring"`.
    pub spring_stiffness: f64,
    /// Spring rest length (m).
    pub spring_rest_length: f64,
    /// Spring damping (N·s/m).
    pub spring_damping: f64,
}

impl Default for WasmJointConfig {
    fn default() -> Self {
        WasmJointConfig {
            joint_type: "fixed".to_string(),
            body_a: 0,
            body_b: 1,
            anchor_a: [0.0; 3],
            anchor_b: [0.0; 3],
            axis_a: [0.0, 1.0, 0.0],
            axis_b: [0.0, 1.0, 0.0],
            lower_limit: -std::f64::consts::PI,
            upper_limit: std::f64::consts::PI,
            limits_enabled: false,
            motor_enabled: false,
            motor_max_force: 100.0,
            motor_target_velocity: 0.0,
            spring_stiffness: 1000.0,
            spring_rest_length: 1.0,
            spring_damping: 10.0,
        }
    }
}

impl WasmJointConfig {
    /// Create a revolute joint.
    pub fn revolute(body_a: u64, body_b: u64, anchor: [f64; 3], axis: [f64; 3]) -> Self {
        WasmJointConfig {
            joint_type: "revolute".to_string(),
            body_a,
            body_b,
            anchor_a: anchor,
            anchor_b: anchor,
            axis_a: axis,
            axis_b: axis,
            ..Default::default()
        }
    }

    /// Create a prismatic joint.
    pub fn prismatic(body_a: u64, body_b: u64, anchor: [f64; 3], axis: [f64; 3]) -> Self {
        WasmJointConfig {
            joint_type: "prismatic".to_string(),
            body_a,
            body_b,
            anchor_a: anchor,
            anchor_b: anchor,
            axis_a: axis,
            axis_b: axis,
            ..Default::default()
        }
    }

    /// Create a ball-and-socket joint.
    pub fn ball(body_a: u64, body_b: u64, anchor: [f64; 3]) -> Self {
        WasmJointConfig {
            joint_type: "ball".to_string(),
            body_a,
            body_b,
            anchor_a: anchor,
            anchor_b: anchor,
            ..Default::default()
        }
    }

    /// Create a spring joint.
    pub fn spring(
        body_a: u64,
        body_b: u64,
        anchor_a: [f64; 3],
        anchor_b: [f64; 3],
        stiffness: f64,
        damping: f64,
    ) -> Self {
        WasmJointConfig {
            joint_type: "spring".to_string(),
            body_a,
            body_b,
            anchor_a,
            anchor_b,
            spring_stiffness: stiffness,
            spring_damping: damping,
            ..Default::default()
        }
    }

    /// Returns `true` if the joint type is recognised.
    pub fn is_valid_type(&self) -> bool {
        matches!(
            self.joint_type.as_str(),
            "revolute" | "prismatic" | "ball" | "fixed" | "spring"
        )
    }
}

#[wasm_bindgen]
impl WasmJointConfig {
    /// Construct a default fixed joint between body 0 and body 1.
    #[wasm_bindgen(constructor)]
    pub fn wasm_new() -> WasmJointConfig {
        WasmJointConfig::default()
    }

    /// Joint type string getter.
    #[wasm_bindgen(getter, js_name = "joint_type")]
    pub fn joint_type_js(&self) -> String {
        self.joint_type.clone()
    }

    /// Joint type string setter.
    #[wasm_bindgen(setter, js_name = "joint_type")]
    pub fn set_joint_type_js(&mut self, joint_type: String) {
        self.joint_type = joint_type;
    }

    /// Anchor on body A as `[x, y, z]`.
    #[wasm_bindgen(js_name = "get_anchor_a")]
    pub fn get_anchor_a(&self) -> Vec<f64> {
        self.anchor_a.to_vec()
    }

    /// Set anchor on body A.
    #[wasm_bindgen(js_name = "set_anchor_a")]
    pub fn set_anchor_a(&mut self, x: f64, y: f64, z: f64) {
        self.anchor_a = [x, y, z];
    }

    /// Anchor on body B as `[x, y, z]`.
    #[wasm_bindgen(js_name = "get_anchor_b")]
    pub fn get_anchor_b(&self) -> Vec<f64> {
        self.anchor_b.to_vec()
    }

    /// Set anchor on body B.
    #[wasm_bindgen(js_name = "set_anchor_b")]
    pub fn set_anchor_b(&mut self, x: f64, y: f64, z: f64) {
        self.anchor_b = [x, y, z];
    }

    /// Axis on body A as `[x, y, z]`.
    #[wasm_bindgen(js_name = "get_axis_a")]
    pub fn get_axis_a(&self) -> Vec<f64> {
        self.axis_a.to_vec()
    }

    /// Set axis on body A.
    #[wasm_bindgen(js_name = "set_axis_a")]
    pub fn set_axis_a(&mut self, x: f64, y: f64, z: f64) {
        self.axis_a = [x, y, z];
    }

    /// Axis on body B as `[x, y, z]`.
    #[wasm_bindgen(js_name = "get_axis_b")]
    pub fn get_axis_b(&self) -> Vec<f64> {
        self.axis_b.to_vec()
    }

    /// Set axis on body B.
    #[wasm_bindgen(js_name = "set_axis_b")]
    pub fn set_axis_b(&mut self, x: f64, y: f64, z: f64) {
        self.axis_b = [x, y, z];
    }

    /// Returns `true` if the joint type string is recognised.
    #[wasm_bindgen(js_name = "is_valid_type")]
    pub fn is_valid_type_js(&self) -> bool {
        self.is_valid_type()
    }

    /// Serialise to a JSON string for consumption in JavaScript.
    #[wasm_bindgen(js_name = "to_json")]
    pub fn to_json_js(&self) -> Result<String, JsValue> {
        serde_json::to_string(self).map_err(err_to_jsvalue)
    }

    /// Serialise to a structured `JsValue` object.
    #[wasm_bindgen(js_name = "to_js_value")]
    pub fn to_js_value_js(&self) -> Result<JsValue, JsValue> {
        to_js_value(self)
    }
}

/// Construct a revolute joint, returned as a structured `JsValue`.
#[wasm_bindgen(js_name = "joint_config_revolute")]
pub fn joint_config_revolute_js(
    body_a: u64,
    body_b: u64,
    anchor_x: f64,
    anchor_y: f64,
    anchor_z: f64,
    axis_x: f64,
    axis_y: f64,
    axis_z: f64,
) -> WasmJointConfig {
    WasmJointConfig::revolute(
        body_a,
        body_b,
        [anchor_x, anchor_y, anchor_z],
        [axis_x, axis_y, axis_z],
    )
}

/// Construct a prismatic joint.
#[wasm_bindgen(js_name = "joint_config_prismatic")]
pub fn joint_config_prismatic_js(
    body_a: u64,
    body_b: u64,
    anchor_x: f64,
    anchor_y: f64,
    anchor_z: f64,
    axis_x: f64,
    axis_y: f64,
    axis_z: f64,
) -> WasmJointConfig {
    WasmJointConfig::prismatic(
        body_a,
        body_b,
        [anchor_x, anchor_y, anchor_z],
        [axis_x, axis_y, axis_z],
    )
}

/// Construct a ball-and-socket joint.
#[wasm_bindgen(js_name = "joint_config_ball")]
pub fn joint_config_ball_js(
    body_a: u64,
    body_b: u64,
    anchor_x: f64,
    anchor_y: f64,
    anchor_z: f64,
) -> WasmJointConfig {
    WasmJointConfig::ball(body_a, body_b, [anchor_x, anchor_y, anchor_z])
}

/// Construct a spring joint.
#[wasm_bindgen(js_name = "joint_config_spring")]
pub fn joint_config_spring_js(
    body_a: u64,
    body_b: u64,
    anchor_a_x: f64,
    anchor_a_y: f64,
    anchor_a_z: f64,
    anchor_b_x: f64,
    anchor_b_y: f64,
    anchor_b_z: f64,
    stiffness: f64,
    damping: f64,
) -> WasmJointConfig {
    WasmJointConfig::spring(
        body_a,
        body_b,
        [anchor_a_x, anchor_a_y, anchor_a_z],
        [anchor_b_x, anchor_b_y, anchor_b_z],
        stiffness,
        damping,
    )
}

// ---------------------------------------------------------------------------
// WasmJointState
// ---------------------------------------------------------------------------

/// Runtime state of a joint queried from the solver.
#[wasm_bindgen]
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct WasmJointState {
    /// Current angle (revolute) or translation (prismatic) in rad or m.
    pub position: f64,
    /// Current velocity in rad/s or m/s.
    pub velocity: f64,
    /// Applied motor torque/force this step (N·m or N).
    pub motor_force: f64,
    /// Constraint reaction force magnitude (N).
    pub constraint_force: f64,
    /// Constraint reaction torque magnitude (N·m).
    pub constraint_torque: f64,
    /// Whether the joint limit is currently active.
    pub limit_active: bool,
    /// Current constraint violation (m or rad).
    pub violation: f64,
}

impl WasmJointState {
    /// Create a zeroed state.
    pub fn new() -> Self {
        Self::default()
    }

    /// Returns `true` if constraint force exceeds `threshold`.
    pub fn is_overloaded(&self, threshold: f64) -> bool {
        self.constraint_force > threshold
    }
}

#[wasm_bindgen]
impl WasmJointState {
    /// Construct a zeroed joint state.
    #[wasm_bindgen(constructor)]
    pub fn wasm_new() -> WasmJointState {
        WasmJointState::default()
    }

    /// Returns `true` if constraint force exceeds `threshold`.
    #[wasm_bindgen(js_name = "is_overloaded")]
    pub fn is_overloaded_js(&self, threshold: f64) -> bool {
        self.is_overloaded(threshold)
    }

    /// Serialise to a structured `JsValue` object.
    #[wasm_bindgen(js_name = "to_js_value")]
    pub fn to_js_value_js(&self) -> Result<JsValue, JsValue> {
        to_js_value(self)
    }
}

// ---------------------------------------------------------------------------
// WasmContactConfig
// ---------------------------------------------------------------------------

/// Per-contact material parameters.
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmContactConfig {
    /// Coefficient of restitution \[0, 1\].
    pub restitution: f64,
    /// Static friction coefficient.
    pub friction_static: f64,
    /// Dynamic (kinetic) friction coefficient.
    pub friction_dynamic: f64,
    /// Rolling friction coefficient.
    pub rolling_friction: f64,
    /// Spinning friction coefficient.
    pub spinning_friction: f64,
    /// Contact compliance (ERP-like parameter).
    pub compliance: f64,
    /// Contact damping (CFM-like parameter).
    pub damping: f64,
}

impl Default for WasmContactConfig {
    fn default() -> Self {
        WasmContactConfig {
            restitution: 0.3,
            friction_static: 0.5,
            friction_dynamic: 0.4,
            rolling_friction: 0.01,
            spinning_friction: 0.005,
            compliance: 0.0,
            damping: 0.0,
        }
    }
}

impl WasmContactConfig {
    /// Create a perfectly rigid, frictionless contact.
    pub fn frictionless() -> Self {
        WasmContactConfig {
            friction_static: 0.0,
            friction_dynamic: 0.0,
            rolling_friction: 0.0,
            spinning_friction: 0.0,
            ..Default::default()
        }
    }

    /// Create a high-friction rubber-like contact.
    pub fn rubber() -> Self {
        WasmContactConfig {
            restitution: 0.5,
            friction_static: 1.0,
            friction_dynamic: 0.8,
            rolling_friction: 0.1,
            spinning_friction: 0.05,
            ..Default::default()
        }
    }

    /// Validate that coefficients are in sensible ranges.
    pub fn validate(&self) -> Result<(), String> {
        if !(0.0..=1.0).contains(&self.restitution) {
            return Err("restitution must be in [0,1]".to_string());
        }
        if self.friction_static < 0.0 {
            return Err("friction_static must be >= 0".to_string());
        }
        Ok(())
    }
}

#[wasm_bindgen]
impl WasmContactConfig {
    /// Construct a default contact configuration.
    #[wasm_bindgen(constructor)]
    pub fn wasm_new() -> WasmContactConfig {
        WasmContactConfig::default()
    }

    /// Construct a perfectly rigid, frictionless contact.
    #[wasm_bindgen(js_name = "frictionless")]
    pub fn frictionless_js() -> WasmContactConfig {
        WasmContactConfig::frictionless()
    }

    /// Construct a high-friction rubber-like contact.
    #[wasm_bindgen(js_name = "rubber")]
    pub fn rubber_js() -> WasmContactConfig {
        WasmContactConfig::rubber()
    }

    /// Validate the configuration; returns the empty string on success.
    #[wasm_bindgen(js_name = "validate")]
    pub fn validate_js(&self) -> String {
        match self.validate() {
            Ok(()) => String::new(),
            Err(e) => e,
        }
    }

    /// Serialise to a structured `JsValue` object.
    #[wasm_bindgen(js_name = "to_js_value")]
    pub fn to_js_value_js(&self) -> Result<JsValue, JsValue> {
        to_js_value(self)
    }
}

// ---------------------------------------------------------------------------
// WasmMotorTarget
// ---------------------------------------------------------------------------

/// Desired motor command for a joint.
///
/// This enum carries payload data and therefore cannot be exposed as a
/// `#[wasm_bindgen]` enum directly. Use the helper free functions
/// `motor_target_*` to build instances on the JavaScript side.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum WasmMotorTarget {
    /// Drive to a target velocity (rad/s or m/s).
    TargetVelocity(f64),
    /// Drive to a target position (rad or m).
    TargetPosition(f64),
    /// Apply a constant torque/force cap (N·m or N).
    MaxTorque(f64),
}

impl WasmMotorTarget {
    /// Return the numeric value regardless of variant.
    pub fn value(&self) -> f64 {
        match self {
            WasmMotorTarget::TargetVelocity(v) => *v,
            WasmMotorTarget::TargetPosition(p) => *p,
            WasmMotorTarget::MaxTorque(t) => *t,
        }
    }

    /// Returns `true` if this is a velocity target.
    pub fn is_velocity(&self) -> bool {
        matches!(self, WasmMotorTarget::TargetVelocity(_))
    }

    /// Returns `true` if this is a position target.
    pub fn is_position(&self) -> bool {
        matches!(self, WasmMotorTarget::TargetPosition(_))
    }
}

/// Build a `TargetVelocity` motor command as a structured `JsValue`.
#[wasm_bindgen(js_name = "motor_target_velocity")]
pub fn motor_target_velocity_js(value: f64) -> Result<JsValue, JsValue> {
    to_js_value(&WasmMotorTarget::TargetVelocity(value))
}

/// Build a `TargetPosition` motor command as a structured `JsValue`.
#[wasm_bindgen(js_name = "motor_target_position")]
pub fn motor_target_position_js(value: f64) -> Result<JsValue, JsValue> {
    to_js_value(&WasmMotorTarget::TargetPosition(value))
}

/// Build a `MaxTorque` motor command as a structured `JsValue`.
#[wasm_bindgen(js_name = "motor_target_max_torque")]
pub fn motor_target_max_torque_js(value: f64) -> Result<JsValue, JsValue> {
    to_js_value(&WasmMotorTarget::MaxTorque(value))
}

/// Numeric value of a JSON-encoded `WasmMotorTarget`.
#[wasm_bindgen(js_name = "motor_target_value")]
pub fn motor_target_value_js(target_json: &str) -> Result<f64, JsValue> {
    let target: WasmMotorTarget = serde_json::from_str(target_json).map_err(err_to_jsvalue)?;
    Ok(target.value())
}

/// Returns `true` when the JSON-encoded `WasmMotorTarget` is a velocity command.
#[wasm_bindgen(js_name = "motor_target_is_velocity")]
pub fn motor_target_is_velocity_js(target_json: &str) -> Result<bool, JsValue> {
    let target: WasmMotorTarget = serde_json::from_str(target_json).map_err(err_to_jsvalue)?;
    Ok(target.is_velocity())
}

/// Returns `true` when the JSON-encoded `WasmMotorTarget` is a position command.
#[wasm_bindgen(js_name = "motor_target_is_position")]
pub fn motor_target_is_position_js(target_json: &str) -> Result<bool, JsValue> {
    let target: WasmMotorTarget = serde_json::from_str(target_json).map_err(err_to_jsvalue)?;
    Ok(target.is_position())
}

// ---------------------------------------------------------------------------
// WasmConstraintSolver
// ---------------------------------------------------------------------------

/// Simple constraint solver managing a collection of joints.
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmConstraintSolver {
    /// Solver iteration count per step.
    pub iterations: u32,
    /// Whether to warm-start from the previous step's lambdas.
    pub warm_start: bool,
    /// All joints currently managed by this solver.
    pub(crate) joints: Vec<(WasmJointHandle, WasmJointConfig)>,
    /// Per-joint state from the last solve.
    pub(crate) joint_states: Vec<WasmJointState>,
    /// Next handle to assign.
    next_handle: u32,
    /// Solver error from the last step.
    pub last_error: f64,
    /// Number of iterations actually used in the last solve.
    pub iterations_used: u32,
}

impl WasmConstraintSolver {
    /// Create a new solver with `iterations` iterations and optional warm-starting.
    pub fn new(iterations: u32, warm_start: bool) -> Self {
        WasmConstraintSolver {
            iterations,
            warm_start,
            joints: Vec::new(),
            joint_states: Vec::new(),
            next_handle: 1,
            last_error: 0.0,
            iterations_used: 0,
        }
    }

    /// Add a joint and return its handle.
    pub fn add_joint(&mut self, config: WasmJointConfig) -> WasmJointHandle {
        let handle = WasmJointHandle::new(self.next_handle);
        self.next_handle += 1;
        self.joints.push((handle, config));
        self.joint_states.push(WasmJointState::new());
        handle
    }

    /// Remove a joint by handle. Returns `true` if found and removed.
    pub fn remove_joint(&mut self, handle: WasmJointHandle) -> bool {
        if let Some(pos) = self.joints.iter().position(|(h, _)| *h == handle) {
            self.joints.remove(pos);
            self.joint_states.remove(pos);
            true
        } else {
            false
        }
    }

    /// Perform one Gauss-Seidel Baumgarte constraint solver pass.
    ///
    /// For each joint, computes a velocity correction impulse based on the
    /// current velocity violation and positional bias (Baumgarte stabilisation),
    /// clamps to the feasible impulse range, and applies the correction.
    pub fn solve(&mut self, dt: f64) {
        self.iterations_used = self.iterations;
        self.last_error = 0.0;
        let safe_dt = dt.max(1e-12);
        let beta = 0.1_f64;
        let regularization = 1e-4_f64;
        // Unit effective mass for the WASM bridge (no rigid-body state).
        let effective_mass = 1.0_f64;
        let inv_m = 1.0 / effective_mass;

        for ((_handle, cfg), state) in self.joints.iter().zip(self.joint_states.iter_mut()) {
            // Integrate position from velocity.
            state.position += state.velocity * safe_dt;

            // Velocity violation: motor joints drive toward target velocity,
            // all other joint types drive velocity toward zero.
            let v_err = if cfg.motor_enabled {
                state.velocity - cfg.motor_target_velocity
            } else {
                state.velocity
            };

            // Baumgarte position stabilisation bias.
            let bias = beta / safe_dt * state.violation;

            // Constraint impulse increment.
            let delta_lambda = -(v_err + bias) / (inv_m + regularization);

            // Clamp to the feasible impulse range.
            let max_force = cfg.motor_max_force.max(0.0);
            let clamped = delta_lambda.clamp(-max_force, max_force);

            // Apply correction.
            state.velocity += clamped * inv_m;
            state.motor_force = clamped;

            // Update violation: reduce residual toward zero.
            state.violation = (state.violation - clamped.abs() * safe_dt * 0.1).max(0.0);
            self.last_error += state.violation;
        }
    }

    /// Get a reference to the state of a joint, if it exists.
    pub fn get_joint_state(&self, handle: WasmJointHandle) -> Option<&WasmJointState> {
        self.joints
            .iter()
            .position(|(h, _)| *h == handle)
            .map(|i| &self.joint_states[i])
    }

    /// Set the motor target for a joint.
    pub fn set_motor_target(&mut self, handle: WasmJointHandle, target: WasmMotorTarget) {
        if let Some(pos) = self.joints.iter().position(|(h, _)| *h == handle) {
            let cfg = &mut self.joints[pos].1;
            match target {
                WasmMotorTarget::TargetVelocity(v) => cfg.motor_target_velocity = v,
                WasmMotorTarget::TargetPosition(p) => cfg.motor_target_velocity = p,
                WasmMotorTarget::MaxTorque(t) => cfg.motor_max_force = t,
            }
        }
    }

    /// Number of joints currently managed.
    pub fn joint_count(&self) -> usize {
        self.joints.len()
    }
}

#[wasm_bindgen]
impl WasmConstraintSolver {
    /// Construct a solver with the given iteration count and warm-start flag.
    #[wasm_bindgen(constructor)]
    pub fn wasm_new(iterations: u32, warm_start: bool) -> WasmConstraintSolver {
        WasmConstraintSolver::new(iterations, warm_start)
    }

    /// Add a joint and return the assigned raw handle id.
    #[wasm_bindgen(js_name = "add_joint")]
    pub fn add_joint_js(&mut self, config: WasmJointConfig) -> u32 {
        self.add_joint(config).raw()
    }

    /// Remove a joint by raw handle id; returns `true` when removed.
    #[wasm_bindgen(js_name = "remove_joint")]
    pub fn remove_joint_js(&mut self, handle_raw: u32) -> bool {
        self.remove_joint(WasmJointHandle::new(handle_raw))
    }

    /// Step the constraint solver forward by `dt` seconds.
    #[wasm_bindgen(js_name = "solve")]
    pub fn solve_js(&mut self, dt: f64) {
        self.solve(dt);
    }

    /// Joint state for the given raw handle as a structured `JsValue`.
    ///
    /// Returns `JsValue::NULL` when the handle is unknown.
    #[wasm_bindgen(js_name = "get_joint_state")]
    pub fn get_joint_state_js(&self, handle_raw: u32) -> Result<JsValue, JsValue> {
        match self.get_joint_state(WasmJointHandle::new(handle_raw)) {
            Some(state) => to_js_value(state),
            None => Ok(JsValue::NULL),
        }
    }

    /// Set the motor target velocity for a joint by raw handle.
    #[wasm_bindgen(js_name = "set_motor_target_velocity")]
    pub fn set_motor_target_velocity_js(&mut self, handle_raw: u32, value: f64) {
        self.set_motor_target(
            WasmJointHandle::new(handle_raw),
            WasmMotorTarget::TargetVelocity(value),
        );
    }

    /// Set the motor target position for a joint by raw handle.
    #[wasm_bindgen(js_name = "set_motor_target_position")]
    pub fn set_motor_target_position_js(&mut self, handle_raw: u32, value: f64) {
        self.set_motor_target(
            WasmJointHandle::new(handle_raw),
            WasmMotorTarget::TargetPosition(value),
        );
    }

    /// Cap the motor torque/force for a joint by raw handle.
    #[wasm_bindgen(js_name = "set_motor_max_torque")]
    pub fn set_motor_max_torque_js(&mut self, handle_raw: u32, value: f64) {
        self.set_motor_target(
            WasmJointHandle::new(handle_raw),
            WasmMotorTarget::MaxTorque(value),
        );
    }

    /// Number of joints currently managed.
    #[wasm_bindgen(js_name = "joint_count")]
    pub fn joint_count_js(&self) -> usize {
        self.joint_count()
    }

    /// Number of (handle, config) pairs currently stored — JS-friendly count.
    #[wasm_bindgen(js_name = "joints_count")]
    pub fn joints_count(&self) -> usize {
        self.joints.len()
    }

    /// Raw handle id of the joint at index `i`, if any.
    #[wasm_bindgen(js_name = "joint_handle_at")]
    pub fn joint_handle_at(&self, i: usize) -> Option<u32> {
        self.joints.get(i).map(|(h, _)| h.raw())
    }

    /// Cloned configuration of the joint at index `i`, if any.
    #[wasm_bindgen(js_name = "joint_config_at")]
    pub fn joint_config_at(&self, i: usize) -> Option<WasmJointConfig> {
        self.joints.get(i).map(|(_, cfg)| cfg.clone())
    }

    /// Cloned state of the joint at index `i`, if any.
    #[wasm_bindgen(js_name = "joint_state_at")]
    pub fn joint_state_at(&self, i: usize) -> Option<WasmJointState> {
        self.joint_states.get(i).cloned()
    }

    /// Serialise the entire solver state to a JSON string.
    #[wasm_bindgen(js_name = "to_json")]
    pub fn to_json_js(&self) -> Result<String, JsValue> {
        serde_json::to_string(self).map_err(err_to_jsvalue)
    }
}

// ---------------------------------------------------------------------------
// WasmConstraintDebug
// ---------------------------------------------------------------------------

/// Debug information for a single constraint.
#[wasm_bindgen]
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct WasmConstraintDebug {
    /// Handle of the joint.
    pub handle: u64,
    /// Constraint impulse (lambda) values from last step.
    #[wasm_bindgen(skip)]
    pub lambdas: Vec<f64>,
    /// Constraint force vector magnitude.
    pub force_magnitude: f64,
    /// Current constraint violation.
    pub violation: f64,
    /// Number of solver iterations used.
    pub iterations_used: u32,
    /// Whether the constraint is active (non-sleeping).
    pub is_active: bool,
}

impl WasmConstraintDebug {
    /// Create debug info for a handle.
    pub fn new(handle: u64) -> Self {
        WasmConstraintDebug {
            handle,
            is_active: true,
            ..Default::default()
        }
    }
}

#[wasm_bindgen]
impl WasmConstraintDebug {
    /// Construct debug info for a joint handle.
    #[wasm_bindgen(constructor)]
    pub fn wasm_new(handle: u64) -> WasmConstraintDebug {
        WasmConstraintDebug::new(handle)
    }

    /// Lambda values from the last solve.
    #[wasm_bindgen(js_name = "get_lambdas")]
    pub fn get_lambdas(&self) -> Vec<f64> {
        self.lambdas.clone()
    }

    /// Replace the lambda values.
    #[wasm_bindgen(js_name = "set_lambdas")]
    pub fn set_lambdas(&mut self, lambdas: Vec<f64>) {
        self.lambdas = lambdas;
    }

    /// Number of stored lambda values.
    #[wasm_bindgen(js_name = "lambdas_count")]
    pub fn lambdas_count(&self) -> usize {
        self.lambdas.len()
    }

    /// Serialise to a structured `JsValue` object.
    #[wasm_bindgen(js_name = "to_js_value")]
    pub fn to_js_value_js(&self) -> Result<JsValue, JsValue> {
        to_js_value(self)
    }
}

// ---------------------------------------------------------------------------
// WasmIslandManager
// ---------------------------------------------------------------------------

/// Statistics and configuration for the simulation island manager.
#[wasm_bindgen]
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct WasmIslandManager {
    /// Number of active simulation islands.
    pub island_count: u32,
    /// Number of bodies per island (sorted by size descending).
    #[wasm_bindgen(skip)]
    pub bodies_per_island: Vec<u32>,
    /// Linear velocity threshold for sleeping (m/s).
    pub sleep_threshold_linear: f64,
    /// Angular velocity threshold for sleeping (rad/s).
    pub sleep_threshold_angular: f64,
    /// Number of consecutive frames below threshold before sleeping.
    pub sleep_delay_frames: u32,
    /// Log of recent wake events (island indices).
    #[wasm_bindgen(skip)]
    pub wake_events: Vec<u32>,
    /// Log of recent sleep events (island indices).
    #[wasm_bindgen(skip)]
    pub sleep_events: Vec<u32>,
}

impl WasmIslandManager {
    /// Create a new island manager with default thresholds.
    pub fn new() -> Self {
        WasmIslandManager {
            sleep_threshold_linear: 0.01,
            sleep_threshold_angular: 0.01,
            sleep_delay_frames: 60,
            ..Default::default()
        }
    }

    /// Register a wake event for island `id`.
    pub fn record_wake(&mut self, island_id: u32) {
        self.wake_events.push(island_id);
    }

    /// Register a sleep event for island `id`.
    pub fn record_sleep(&mut self, island_id: u32) {
        self.sleep_events.push(island_id);
    }

    /// Total bodies across all islands.
    pub fn total_bodies(&self) -> u32 {
        self.bodies_per_island.iter().sum()
    }

    /// Largest island body count.
    pub fn largest_island(&self) -> u32 {
        self.bodies_per_island.iter().copied().max().unwrap_or(0)
    }

    /// Clear event logs.
    pub fn flush_events(&mut self) {
        self.wake_events.clear();
        self.sleep_events.clear();
    }
}

#[wasm_bindgen]
impl WasmIslandManager {
    /// Construct a manager with default sleep thresholds.
    #[wasm_bindgen(constructor)]
    pub fn wasm_new() -> WasmIslandManager {
        WasmIslandManager::new()
    }

    /// Record a wake event for the given island id.
    #[wasm_bindgen(js_name = "record_wake")]
    pub fn record_wake_js(&mut self, island_id: u32) {
        self.record_wake(island_id);
    }

    /// Record a sleep event for the given island id.
    #[wasm_bindgen(js_name = "record_sleep")]
    pub fn record_sleep_js(&mut self, island_id: u32) {
        self.record_sleep(island_id);
    }

    /// Total bodies across all islands.
    #[wasm_bindgen(js_name = "total_bodies")]
    pub fn total_bodies_js(&self) -> u32 {
        self.total_bodies()
    }

    /// Largest island body count.
    #[wasm_bindgen(js_name = "largest_island")]
    pub fn largest_island_js(&self) -> u32 {
        self.largest_island()
    }

    /// Clear event logs.
    #[wasm_bindgen(js_name = "flush_events")]
    pub fn flush_events_js(&mut self) {
        self.flush_events();
    }

    /// Snapshot of the bodies-per-island vector.
    #[wasm_bindgen(js_name = "get_bodies_per_island")]
    pub fn get_bodies_per_island(&self) -> Vec<u32> {
        self.bodies_per_island.clone()
    }

    /// Replace the bodies-per-island counts.
    #[wasm_bindgen(js_name = "set_bodies_per_island")]
    pub fn set_bodies_per_island(&mut self, counts: Vec<u32>) {
        self.bodies_per_island = counts;
    }

    /// Snapshot of pending wake events.
    #[wasm_bindgen(js_name = "get_wake_events")]
    pub fn get_wake_events(&self) -> Vec<u32> {
        self.wake_events.clone()
    }

    /// Snapshot of pending sleep events.
    #[wasm_bindgen(js_name = "get_sleep_events")]
    pub fn get_sleep_events(&self) -> Vec<u32> {
        self.sleep_events.clone()
    }

    /// Serialise to a structured `JsValue` object.
    #[wasm_bindgen(js_name = "to_js_value")]
    pub fn to_js_value_js(&self) -> Result<JsValue, JsValue> {
        to_js_value(self)
    }
}

// ---------------------------------------------------------------------------
// WasmCcdConfig
// ---------------------------------------------------------------------------

/// Continuous collision detection configuration.
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmCcdConfig {
    /// Whether CCD is globally enabled.
    pub enabled: bool,
    /// Maximum time-of-impact value searched \[0, 1\].
    pub max_toi: f64,
    /// Maximum number of CCD sub-steps per timestep.
    pub substeps: u32,
    /// Speculative contact margin (m).
    pub speculative_margin: f64,
    /// Minimum relative velocity to trigger CCD (m/s).
    pub min_velocity_threshold: f64,
    /// Whether to use motion clamping mode.
    pub motion_clamping: bool,
}

impl Default for WasmCcdConfig {
    fn default() -> Self {
        WasmCcdConfig {
            enabled: false,
            max_toi: 1.0,
            substeps: 4,
            speculative_margin: 0.001,
            min_velocity_threshold: 1.0,
            motion_clamping: true,
        }
    }
}

impl WasmCcdConfig {
    /// Create a CCD config with CCD enabled.
    pub fn enabled() -> Self {
        WasmCcdConfig {
            enabled: true,
            ..Default::default()
        }
    }

    /// Validate the configuration.
    pub fn validate(&self) -> Result<(), String> {
        if !(0.0..=1.0).contains(&self.max_toi) {
            return Err("max_toi must be in [0,1]".to_string());
        }
        if self.substeps == 0 {
            return Err("substeps must be >= 1".to_string());
        }
        Ok(())
    }
}

#[wasm_bindgen]
impl WasmCcdConfig {
    /// Construct a default (disabled) CCD configuration.
    #[wasm_bindgen(constructor)]
    pub fn wasm_new() -> WasmCcdConfig {
        WasmCcdConfig::default()
    }

    /// Construct an enabled CCD configuration with default parameters.
    #[wasm_bindgen(js_name = "create_enabled")]
    pub fn create_enabled() -> WasmCcdConfig {
        WasmCcdConfig::enabled()
    }

    /// Validate the configuration; returns the empty string on success.
    #[wasm_bindgen(js_name = "validate")]
    pub fn validate_js(&self) -> String {
        match self.validate() {
            Ok(()) => String::new(),
            Err(e) => e,
        }
    }

    /// Serialise to a structured `JsValue` object.
    #[wasm_bindgen(js_name = "to_js_value")]
    pub fn to_js_value_js(&self) -> Result<JsValue, JsValue> {
        to_js_value(self)
    }
}

// ---------------------------------------------------------------------------
// WasmRagdollBuilder
// ---------------------------------------------------------------------------

/// Description of a single bone in a ragdoll.
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmBone {
    /// Unique name for this bone.
    #[wasm_bindgen(skip)]
    pub name: String,
    /// Bone length in metres.
    pub length: f64,
    /// Bone mass in kilograms.
    pub mass: f64,
    /// Radius of the capsule collider.
    pub radius: f64,
    /// Index of the parent bone (-1 for root).
    pub parent_index: i32,
    /// Joint connecting this bone to its parent.
    #[wasm_bindgen(skip)]
    pub joint: Option<WasmJointConfig>,
}

impl WasmBone {
    /// Create a new bone.
    pub fn new(name: String, length: f64, mass: f64) -> Self {
        WasmBone {
            name,
            length,
            mass,
            radius: length * 0.1,
            parent_index: -1,
            joint: None,
        }
    }
}

#[wasm_bindgen]
impl WasmBone {
    /// Construct a new bone with the given name, length (m) and mass (kg).
    #[wasm_bindgen(constructor)]
    pub fn wasm_new(name: String, length: f64, mass: f64) -> WasmBone {
        WasmBone::new(name, length, mass)
    }

    /// Bone name.
    #[wasm_bindgen(getter)]
    pub fn name(&self) -> String {
        self.name.clone()
    }

    /// Set the bone name.
    #[wasm_bindgen(setter)]
    pub fn set_name(&mut self, name: String) {
        self.name = name;
    }

    /// Returns `true` when this bone has a joint connecting it to its parent.
    #[wasm_bindgen(js_name = "has_joint")]
    pub fn has_joint(&self) -> bool {
        self.joint.is_some()
    }

    /// Cloned joint configuration, if any.
    #[wasm_bindgen(js_name = "get_joint")]
    pub fn get_joint(&self) -> Option<WasmJointConfig> {
        self.joint.clone()
    }

    /// Replace the joint configuration. Pass `None` to clear (use
    /// `clear_joint` from JS).
    pub fn set_joint(&mut self, joint: WasmJointConfig) {
        self.joint = Some(joint);
    }

    /// Remove the joint configuration.
    #[wasm_bindgen(js_name = "clear_joint")]
    pub fn clear_joint(&mut self) {
        self.joint = None;
    }

    /// Serialise to a structured `JsValue` object.
    #[wasm_bindgen(js_name = "to_js_value")]
    pub fn to_js_value_js(&self) -> Result<JsValue, JsValue> {
        to_js_value(self)
    }
}

/// Builder for constructing a ragdoll hierarchy.
#[wasm_bindgen]
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct WasmRagdollBuilder {
    /// Bones added so far.
    pub(crate) bones: Vec<WasmBone>,
}

impl WasmRagdollBuilder {
    /// Create an empty builder.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a bone and return its index.
    pub fn add_bone(&mut self, name: String, length: f64, mass: f64) -> usize {
        let idx = self.bones.len();
        self.bones.push(WasmBone::new(name, length, mass));
        idx
    }

    /// Connect bone at `child_index` to bone at `parent_index` with the given joint.
    ///
    /// Returns `Err` if either index is out of bounds.
    pub fn connect_joint(
        &mut self,
        child_index: usize,
        parent_index: usize,
        joint: WasmJointConfig,
    ) -> Result<(), String> {
        if child_index >= self.bones.len() {
            return Err(format!("child_index {child_index} out of bounds"));
        }
        if parent_index >= self.bones.len() {
            return Err(format!("parent_index {parent_index} out of bounds"));
        }
        self.bones[child_index].parent_index = parent_index as i32;
        self.bones[child_index].joint = Some(joint);
        Ok(())
    }

    /// Build the ragdoll, returning a handle per bone joint (root has no joint handle).
    ///
    /// Handles are assigned sequentially starting from 1.
    pub fn build(&self) -> Vec<WasmJointHandle> {
        self.bones
            .iter()
            .enumerate()
            .filter_map(|(i, bone)| {
                if bone.joint.is_some() {
                    Some(WasmJointHandle::new(i as u32 + 1))
                } else {
                    None
                }
            })
            .collect()
    }

    /// Number of bones.
    pub fn bone_count(&self) -> usize {
        self.bones.len()
    }

    /// Find a bone by name.
    pub fn find_bone(&self, name: &str) -> Option<usize> {
        self.bones.iter().position(|b| b.name == name)
    }
}

#[wasm_bindgen]
impl WasmRagdollBuilder {
    /// Construct an empty ragdoll builder.
    #[wasm_bindgen(constructor)]
    pub fn wasm_new() -> WasmRagdollBuilder {
        WasmRagdollBuilder::new()
    }

    /// Add a bone and return its index.
    #[wasm_bindgen(js_name = "add_bone")]
    pub fn add_bone_js(&mut self, name: String, length: f64, mass: f64) -> usize {
        self.add_bone(name, length, mass)
    }

    /// Connect bone `child_index` to bone `parent_index` with the given joint.
    ///
    /// Returns the empty string on success, or the error message on failure.
    #[wasm_bindgen(js_name = "connect_joint")]
    pub fn connect_joint_js(
        &mut self,
        child_index: usize,
        parent_index: usize,
        joint: WasmJointConfig,
    ) -> String {
        match self.connect_joint(child_index, parent_index, joint) {
            Ok(()) => String::new(),
            Err(e) => e,
        }
    }

    /// Build the ragdoll, returning the raw handle id for every bone joint.
    #[wasm_bindgen(js_name = "build")]
    pub fn build_js(&self) -> Vec<u32> {
        self.build().into_iter().map(|h| h.raw()).collect()
    }

    /// Number of bones.
    #[wasm_bindgen(js_name = "bone_count")]
    pub fn bone_count_js(&self) -> usize {
        self.bone_count()
    }

    /// Find a bone by name; returns the index or `usize::MAX` if not found.
    #[wasm_bindgen(js_name = "find_bone")]
    pub fn find_bone_js(&self, name: &str) -> Option<usize> {
        self.find_bone(name)
    }

    /// Cloned bone at index `i`, if any.
    #[wasm_bindgen(js_name = "bone_at")]
    pub fn bone_at(&self, i: usize) -> Option<WasmBone> {
        self.bones.get(i).cloned()
    }

    /// Number of bones (synonym of `bone_count`).
    #[wasm_bindgen(js_name = "bones_count")]
    pub fn bones_count(&self) -> usize {
        self.bones.len()
    }

    /// Serialise the entire builder state to a JSON string.
    #[wasm_bindgen(js_name = "to_json")]
    pub fn to_json_js(&self) -> Result<String, JsValue> {
        serde_json::to_string(self).map_err(err_to_jsvalue)
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(not(target_arch = "wasm32"))]
#[cfg(test)]
mod tests {
    use super::*;

    // --- WasmJointHandle ---

    #[test]
    fn test_handle_new_and_raw() {
        let h = WasmJointHandle::new(42);
        assert_eq!(h.raw(), 42);
    }

    #[test]
    fn test_handle_invalid() {
        let h = WasmJointHandle::invalid();
        assert!(h.is_invalid());
    }

    #[test]
    fn test_handle_valid_is_not_invalid() {
        let h = WasmJointHandle::new(1);
        assert!(!h.is_invalid());
    }

    // --- WasmJointConfig ---

    #[test]
    fn test_joint_config_revolute() {
        let cfg = WasmJointConfig::revolute(0, 1, [0.0; 3], [0.0, 1.0, 0.0]);
        assert_eq!(cfg.joint_type, "revolute");
        assert!(cfg.is_valid_type());
    }

    #[test]
    fn test_joint_config_prismatic() {
        let cfg = WasmJointConfig::prismatic(0, 1, [0.0; 3], [1.0, 0.0, 0.0]);
        assert_eq!(cfg.joint_type, "prismatic");
        assert!(cfg.is_valid_type());
    }

    #[test]
    fn test_joint_config_ball() {
        let cfg = WasmJointConfig::ball(0, 1, [0.5, 0.0, 0.0]);
        assert_eq!(cfg.joint_type, "ball");
        assert!(cfg.is_valid_type());
    }

    #[test]
    fn test_joint_config_spring() {
        let cfg = WasmJointConfig::spring(0, 1, [0.0; 3], [1.0; 3], 500.0, 5.0);
        assert_eq!(cfg.joint_type, "spring");
        assert!((cfg.spring_stiffness - 500.0).abs() < 1e-10);
    }

    #[test]
    fn test_joint_config_invalid_type() {
        let cfg = WasmJointConfig {
            joint_type: "unknown".to_string(),
            ..Default::default()
        };
        assert!(!cfg.is_valid_type());
    }

    #[test]
    fn test_joint_config_serialization() {
        let cfg = WasmJointConfig::revolute(0, 1, [0.0; 3], [0.0, 1.0, 0.0]);
        let json = serde_json::to_string(&cfg).expect("serialize joint config");
        let cfg2: WasmJointConfig = serde_json::from_str(&json).expect("deserialize joint config");
        assert_eq!(cfg2.joint_type, "revolute");
    }

    // --- WasmJointState ---

    #[test]
    fn test_joint_state_overloaded() {
        let mut s = WasmJointState::new();
        s.constraint_force = 1000.0;
        assert!(s.is_overloaded(500.0));
        assert!(!s.is_overloaded(2000.0));
    }

    // --- WasmContactConfig ---

    #[test]
    fn test_contact_frictionless() {
        let cfg = WasmContactConfig::frictionless();
        assert_eq!(cfg.friction_static, 0.0);
        assert_eq!(cfg.friction_dynamic, 0.0);
    }

    #[test]
    fn test_contact_rubber() {
        let cfg = WasmContactConfig::rubber();
        assert!(cfg.friction_static >= 1.0);
    }

    #[test]
    fn test_contact_validate_ok() {
        let cfg = WasmContactConfig::default();
        assert!(cfg.validate().is_ok());
    }

    #[test]
    fn test_contact_validate_bad_restitution() {
        let cfg = WasmContactConfig {
            restitution: 1.5,
            ..Default::default()
        };
        assert!(cfg.validate().is_err());
    }

    // --- WasmMotorTarget ---

    #[test]
    fn test_motor_target_velocity() {
        let t = WasmMotorTarget::TargetVelocity(2.72);
        assert!(t.is_velocity());
        assert!((t.value() - 2.72).abs() < 1e-10);
    }

    #[test]
    fn test_motor_target_position() {
        let t = WasmMotorTarget::TargetPosition(1.0);
        assert!(t.is_position());
    }

    #[test]
    fn test_motor_target_max_torque() {
        let t = WasmMotorTarget::MaxTorque(50.0);
        assert!(!t.is_velocity());
        assert!(!t.is_position());
        assert!((t.value() - 50.0).abs() < 1e-10);
    }

    // --- WasmConstraintSolver ---

    #[test]
    fn test_solver_add_remove_joint() {
        let mut solver = WasmConstraintSolver::new(10, true);
        let cfg = WasmJointConfig::revolute(0, 1, [0.0; 3], [0.0, 1.0, 0.0]);
        let h = solver.add_joint(cfg);
        assert_eq!(solver.joint_count(), 1);
        assert!(solver.remove_joint(h));
        assert_eq!(solver.joint_count(), 0);
    }

    #[test]
    fn test_solver_remove_nonexistent() {
        let mut solver = WasmConstraintSolver::new(10, false);
        assert!(!solver.remove_joint(WasmJointHandle::new(99)));
    }

    #[test]
    fn test_solver_get_joint_state() {
        let mut solver = WasmConstraintSolver::new(10, true);
        let h = solver.add_joint(WasmJointConfig::default());
        let state = solver.get_joint_state(h);
        assert!(state.is_some());
    }

    #[test]
    fn test_solver_solve_advances_position() {
        let mut solver = WasmConstraintSolver::new(10, true);
        let h = solver.add_joint(WasmJointConfig::default());
        if let Some(pos) = solver.joints.iter().position(|(hh, _)| *hh == h) {
            solver.joint_states[pos].velocity = 1.0;
        }
        solver.solve(0.1);
        let state = solver
            .get_joint_state(h)
            .expect("joint state must exist after add_joint");
        assert!((state.position - 0.1).abs() < 1e-10);
    }

    #[test]
    fn test_solver_set_motor_target() {
        let mut solver = WasmConstraintSolver::new(10, true);
        let h = solver.add_joint(WasmJointConfig::default());
        solver.set_motor_target(h, WasmMotorTarget::TargetVelocity(2.0));
        let (_, cfg) = solver
            .joints
            .iter()
            .find(|(hh, _)| *hh == h)
            .expect("joint must exist after add_joint");
        assert!((cfg.motor_target_velocity - 2.0).abs() < 1e-10);
    }

    // --- WasmIslandManager ---

    #[test]
    fn test_island_total_bodies() {
        let mut mgr = WasmIslandManager::new();
        mgr.bodies_per_island = vec![10, 5, 3];
        assert_eq!(mgr.total_bodies(), 18);
    }

    #[test]
    fn test_island_largest() {
        let mut mgr = WasmIslandManager::new();
        mgr.bodies_per_island = vec![4, 12, 7];
        assert_eq!(mgr.largest_island(), 12);
    }

    #[test]
    fn test_island_events() {
        let mut mgr = WasmIslandManager::new();
        mgr.record_wake(0);
        mgr.record_sleep(1);
        assert_eq!(mgr.wake_events.len(), 1);
        mgr.flush_events();
        assert!(mgr.wake_events.is_empty());
    }

    // --- WasmCcdConfig ---

    #[test]
    fn test_ccd_default_disabled() {
        let cfg = WasmCcdConfig::default();
        assert!(!cfg.enabled);
    }

    #[test]
    fn test_ccd_enabled_constructor() {
        let cfg = WasmCcdConfig::enabled();
        assert!(cfg.enabled);
    }

    #[test]
    fn test_ccd_validate_ok() {
        let cfg = WasmCcdConfig::default();
        assert!(cfg.validate().is_ok());
    }

    #[test]
    fn test_ccd_validate_bad_max_toi() {
        let cfg = WasmCcdConfig {
            max_toi: 2.0,
            ..Default::default()
        };
        assert!(cfg.validate().is_err());
    }

    // --- WasmRagdollBuilder ---

    #[test]
    fn test_ragdoll_add_bones() {
        let mut rb = WasmRagdollBuilder::new();
        let _torso = rb.add_bone("torso".to_string(), 0.5, 10.0);
        let _head = rb.add_bone("head".to_string(), 0.25, 3.0);
        assert_eq!(rb.bone_count(), 2);
    }

    #[test]
    fn test_ragdoll_connect_joint() {
        let mut rb = WasmRagdollBuilder::new();
        let torso = rb.add_bone("torso".to_string(), 0.5, 10.0);
        let head = rb.add_bone("head".to_string(), 0.25, 3.0);
        let joint = WasmJointConfig::ball(torso as u64, head as u64, [0.0, 0.25, 0.0]);
        assert!(rb.connect_joint(head, torso, joint).is_ok());
        assert_eq!(rb.bones[head].parent_index, torso as i32);
    }

    #[test]
    fn test_ragdoll_connect_out_of_bounds() {
        let mut rb = WasmRagdollBuilder::new();
        let _ = rb.add_bone("root".to_string(), 1.0, 5.0);
        let joint = WasmJointConfig::default();
        assert!(rb.connect_joint(99, 0, joint).is_err());
    }

    #[test]
    fn test_ragdoll_build_returns_handles() {
        let mut rb = WasmRagdollBuilder::new();
        let torso = rb.add_bone("torso".to_string(), 0.5, 10.0);
        let head = rb.add_bone("head".to_string(), 0.25, 3.0);
        let joint = WasmJointConfig::ball(torso as u64, head as u64, [0.0, 0.25, 0.0]);
        rb.connect_joint(head, torso, joint)
            .expect("connect_joint must succeed for valid indices");
        let handles = rb.build();
        assert_eq!(handles.len(), 1);
    }

    #[test]
    fn test_ragdoll_find_bone() {
        let mut rb = WasmRagdollBuilder::new();
        rb.add_bone("spine".to_string(), 0.4, 8.0);
        assert_eq!(rb.find_bone("spine"), Some(0));
        assert!(rb.find_bone("missing").is_none());
    }

    // --- WasmConstraintDebug ---

    #[test]
    fn test_constraint_debug_new() {
        let dbg = WasmConstraintDebug::new(5);
        assert_eq!(dbg.handle, 5);
        assert!(dbg.is_active);
    }

    #[test]
    fn test_constraint_debug_serialization() {
        let dbg = WasmConstraintDebug::new(3);
        let json = serde_json::to_string(&dbg).expect("serialize constraint debug");
        let dbg2: WasmConstraintDebug =
            serde_json::from_str(&json).expect("deserialize constraint debug");
        assert_eq!(dbg2.handle, 3);
    }

    // --- H5: Gauss-Seidel Baumgarte solver ---

    #[test]
    fn solver_reduces_violation() {
        let mut solver = WasmConstraintSolver::new(1, false);
        let cfg = WasmJointConfig {
            motor_enabled: false,
            motor_max_force: 1000.0,
            ..Default::default()
        };
        let h = solver.add_joint(cfg);

        // Manually set a non-zero violation.
        if let Some(pos) = solver.joints.iter().position(|(hh, _)| *hh == h) {
            solver.joint_states[pos].violation = 0.5;
            solver.joint_states[pos].velocity = 1.0;
        }

        solver.solve(0.01);

        // After the corrective pass, last_error must be less than the initial violation.
        assert!(
            solver.last_error < 0.5,
            "last_error ({}) should be less than initial violation 0.5",
            solver.last_error
        );
    }
}
