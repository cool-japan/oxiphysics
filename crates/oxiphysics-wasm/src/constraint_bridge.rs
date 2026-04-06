// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! WebAssembly constraint system bridge.
//!
//! Provides pure-Rust wrappers around constraint/joint data structures
//! intended for serialisation to and from JavaScript via wasm-bindgen.
//! No wasm-bindgen annotations are placed here; the structs are plain Rust.

#![allow(dead_code)]

use serde::{Deserialize, Serialize};

// ---------------------------------------------------------------------------
// WasmJointHandle
// ---------------------------------------------------------------------------

/// Opaque 64-bit handle identifying a joint in the constraint solver.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct WasmJointHandle(pub u64);

impl WasmJointHandle {
    /// Create a handle from a raw `u64`.
    pub fn new(id: u64) -> Self {
        WasmJointHandle(id)
    }

    /// Return the raw identifier.
    pub fn raw(&self) -> u64 {
        self.0
    }

    /// Sentinel representing an invalid/null handle.
    pub fn invalid() -> Self {
        WasmJointHandle(u64::MAX)
    }

    /// Returns `true` if this handle is the invalid sentinel.
    pub fn is_invalid(&self) -> bool {
        self.0 == u64::MAX
    }
}

// ---------------------------------------------------------------------------
// WasmJointConfig
// ---------------------------------------------------------------------------

/// Configuration for a joint connecting two bodies.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmJointConfig {
    /// Joint type: `"revolute"`, `"prismatic"`, `"ball"`, `"fixed"`, or `"spring"`.
    pub joint_type: String,
    /// Identifier of body A (index or handle).
    pub body_a: u64,
    /// Identifier of body B (index or handle).
    pub body_b: u64,
    /// Anchor point on body A in local space \[x, y, z\].
    pub anchor_a: [f64; 3],
    /// Anchor point on body B in local space \[x, y, z\].
    pub anchor_b: [f64; 3],
    /// Joint axis in body A local space \[x, y, z\] (used by revolute/prismatic).
    pub axis_a: [f64; 3],
    /// Joint axis in body B local space \[x, y, z\].
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

// ---------------------------------------------------------------------------
// WasmJointState
// ---------------------------------------------------------------------------

/// Runtime state of a joint queried from the solver.
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

// ---------------------------------------------------------------------------
// WasmContactConfig
// ---------------------------------------------------------------------------

/// Per-contact material parameters.
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

// ---------------------------------------------------------------------------
// WasmMotorTarget
// ---------------------------------------------------------------------------

/// Desired motor command for a joint.
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

// ---------------------------------------------------------------------------
// WasmConstraintSolver
// ---------------------------------------------------------------------------

/// Simple constraint solver managing a collection of joints.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmConstraintSolver {
    /// Solver iteration count per step.
    pub iterations: u32,
    /// Whether to warm-start from the previous step's lambdas.
    pub warm_start: bool,
    /// All joints currently managed by this solver.
    pub joints: Vec<(WasmJointHandle, WasmJointConfig)>,
    /// Per-joint state from the last solve.
    pub joint_states: Vec<WasmJointState>,
    /// Next handle to assign.
    next_handle: u64,
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

    /// Perform one solver pass (stub — updates violation bookkeeping).
    pub fn solve(&mut self, dt: f64) {
        self.iterations_used = self.iterations;
        self.last_error = 0.0;
        for state in &mut self.joint_states {
            // Integrate position from velocity (stub).
            state.position += state.velocity * dt;
            state.violation = state.position.abs() * 0.001;
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

// ---------------------------------------------------------------------------
// WasmConstraintDebug
// ---------------------------------------------------------------------------

/// Debug information for a single constraint.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct WasmConstraintDebug {
    /// Handle of the joint.
    pub handle: u64,
    /// Constraint impulse (lambda) values from last step.
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

// ---------------------------------------------------------------------------
// WasmIslandManager
// ---------------------------------------------------------------------------

/// Statistics and configuration for the simulation island manager.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct WasmIslandManager {
    /// Number of active simulation islands.
    pub island_count: u32,
    /// Number of bodies per island (sorted by size descending).
    pub bodies_per_island: Vec<u32>,
    /// Linear velocity threshold for sleeping (m/s).
    pub sleep_threshold_linear: f64,
    /// Angular velocity threshold for sleeping (rad/s).
    pub sleep_threshold_angular: f64,
    /// Number of consecutive frames below threshold before sleeping.
    pub sleep_delay_frames: u32,
    /// Log of recent wake events (island indices).
    pub wake_events: Vec<u32>,
    /// Log of recent sleep events (island indices).
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

// ---------------------------------------------------------------------------
// WasmCcdConfig
// ---------------------------------------------------------------------------

/// Continuous collision detection configuration.
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

// ---------------------------------------------------------------------------
// WasmRagdollBuilder
// ---------------------------------------------------------------------------

/// Description of a single bone in a ragdoll.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmBone {
    /// Unique name for this bone.
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
    pub joint: Option<WasmJointConfig>,
}

impl WasmBone {
    /// Create a new bone.
    pub fn new(name: impl Into<String>, length: f64, mass: f64) -> Self {
        WasmBone {
            name: name.into(),
            length,
            mass,
            radius: length * 0.1,
            parent_index: -1,
            joint: None,
        }
    }
}

/// Builder for constructing a ragdoll hierarchy.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct WasmRagdollBuilder {
    /// Bones added so far.
    pub bones: Vec<WasmBone>,
}

impl WasmRagdollBuilder {
    /// Create an empty builder.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a bone and return its index.
    pub fn add_bone(&mut self, name: impl Into<String>, length: f64, mass: f64) -> usize {
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
                    Some(WasmJointHandle::new(i as u64 + 1))
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
        let json = serde_json::to_string(&cfg).unwrap();
        let cfg2: WasmJointConfig = serde_json::from_str(&json).unwrap();
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
        let state = solver.get_joint_state(h).unwrap();
        assert!((state.position - 0.1).abs() < 1e-10);
    }

    #[test]
    fn test_solver_set_motor_target() {
        let mut solver = WasmConstraintSolver::new(10, true);
        let h = solver.add_joint(WasmJointConfig::default());
        solver.set_motor_target(h, WasmMotorTarget::TargetVelocity(2.0));
        let (_, cfg) = solver.joints.iter().find(|(hh, _)| *hh == h).unwrap();
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
        let _torso = rb.add_bone("torso", 0.5, 10.0);
        let _head = rb.add_bone("head", 0.25, 3.0);
        assert_eq!(rb.bone_count(), 2);
    }

    #[test]
    fn test_ragdoll_connect_joint() {
        let mut rb = WasmRagdollBuilder::new();
        let torso = rb.add_bone("torso", 0.5, 10.0);
        let head = rb.add_bone("head", 0.25, 3.0);
        let joint = WasmJointConfig::ball(torso as u64, head as u64, [0.0, 0.25, 0.0]);
        assert!(rb.connect_joint(head, torso, joint).is_ok());
        assert_eq!(rb.bones[head].parent_index, torso as i32);
    }

    #[test]
    fn test_ragdoll_connect_out_of_bounds() {
        let mut rb = WasmRagdollBuilder::new();
        let _ = rb.add_bone("root", 1.0, 5.0);
        let joint = WasmJointConfig::default();
        assert!(rb.connect_joint(99, 0, joint).is_err());
    }

    #[test]
    fn test_ragdoll_build_returns_handles() {
        let mut rb = WasmRagdollBuilder::new();
        let torso = rb.add_bone("torso", 0.5, 10.0);
        let head = rb.add_bone("head", 0.25, 3.0);
        let joint = WasmJointConfig::ball(torso as u64, head as u64, [0.0, 0.25, 0.0]);
        rb.connect_joint(head, torso, joint).unwrap();
        let handles = rb.build();
        assert_eq!(handles.len(), 1);
    }

    #[test]
    fn test_ragdoll_find_bone() {
        let mut rb = WasmRagdollBuilder::new();
        rb.add_bone("spine", 0.4, 8.0);
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
        let json = serde_json::to_string(&dbg).unwrap();
        let dbg2: WasmConstraintDebug = serde_json::from_str(&json).unwrap();
        assert_eq!(dbg2.handle, 3);
    }
}
