//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

use super::types::{BroadPhaseAlgorithm, IntegrationMethod, WasmColliderShape, WasmJointHandle};
use crate::wasm_helpers::to_js_value;

/// A collider attached to a rigid body.
#[derive(Debug, Clone)]
pub(super) struct WasmCollider {
    pub(super) handle: WasmColliderHandle,
    pub(super) body: Option<WasmRigidBodyHandle>,
    pub(super) shape: WasmColliderShape,
    pub(super) friction: f64,
    pub(super) restitution: f64,
    pub(super) is_sensor: bool,
}
impl WasmCollider {
    pub(super) fn new(handle: WasmColliderHandle, shape: WasmColliderShape) -> Self {
        Self {
            handle,
            body: None,
            shape,
            friction: 0.5,
            restitution: 0.3,
            is_sensor: false,
        }
    }
}
/// State of a rigid body for serialization and queries.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmBodyState {
    /// Handle of the body.
    pub handle: WasmRigidBodyHandle,
    /// World-space position.
    pub position: [f64; 3],
    /// Orientation as quaternion (x, y, z, w).
    pub rotation: [f64; 4],
    /// Linear velocity.
    pub linear_vel: [f64; 3],
    /// Angular velocity.
    pub angular_vel: [f64; 3],
    /// Whether the body is currently sleeping.
    pub sleeping: bool,
    /// Body type.
    pub body_type: WasmBodyType,
}
/// Type of a rigid body.
#[wasm_bindgen]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum WasmBodyType {
    /// Fully simulated (mass > 0).
    Dynamic,
    /// Immovable (infinite mass).
    Static,
    /// Driven externally (velocity-controlled).
    Kinematic,
}
/// Handle to a rigid body inside WasmWorld.
#[wasm_bindgen]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct WasmRigidBodyHandle(pub u32);
#[wasm_bindgen]
impl WasmRigidBodyHandle {
    /// Construct a rigid body handle from a raw `u32` id.
    #[wasm_bindgen(constructor)]
    pub fn wasm_new(id: u32) -> Self {
        WasmRigidBodyHandle(id)
    }
    /// Return the raw identifier.
    #[wasm_bindgen(getter)]
    pub fn id(&self) -> u32 {
        self.0
    }
}
/// Joint/constraint type.
#[wasm_bindgen]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum WasmJointType {
    /// Ball-and-socket joint (3 DOF rotation).
    Ball,
    /// Fixed joint (0 DOF).
    Fixed,
    /// Revolute joint (1 DOF rotation).
    Revolute,
    /// Prismatic joint (1 DOF translation).
    Prismatic,
    /// Generic 6-DOF joint.
    Generic6Dof,
    /// Distance constraint.
    Distance,
    /// Spring-damper constraint.
    Spring,
}
/// Configuration for a WASM simulation world.
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmSimulationConfig {
    /// Gravity vector (m/s²).
    #[wasm_bindgen(skip)]
    pub gravity: [f64; 3],
    /// Fixed simulation time step (s).
    pub dt: f64,
    /// Number of velocity solver iterations.
    pub velocity_iters: u32,
    /// Number of position solver iterations.
    pub position_iters: u32,
    /// Broad-phase algorithm.
    pub broad_phase: BroadPhaseAlgorithm,
    /// Integration method.
    pub integration: IntegrationMethod,
    /// Enable continuous collision detection.
    pub ccd_enabled: bool,
    /// Maximum linear velocity (m/s) before clamping.
    pub max_linear_vel: f64,
    /// Maximum angular velocity (rad/s) before clamping.
    pub max_angular_vel: f64,
    /// Linear damping applied to all bodies.
    pub linear_damping: f64,
    /// Angular damping applied to all bodies.
    pub angular_damping: f64,
    /// Allow sleeping of inactive bodies.
    pub allow_sleeping: bool,
    /// Linear speed below which body enters sleep (m/s).
    pub sleep_linear_threshold: f64,
    /// Angular speed below which body enters sleep (rad/s).
    pub sleep_angular_threshold: f64,
    /// Time before sleep activates (s).
    pub sleep_time_threshold: f64,
}
impl WasmSimulationConfig {
    /// Create a config with Earth-standard gravity.
    pub fn earth() -> Self {
        Self::default()
    }
    /// Create a config with Moon-level gravity (1/6 of Earth).
    pub fn moon() -> Self {
        Self {
            gravity: [0.0, -1.625, 0.0],
            ..Default::default()
        }
    }
    /// Create a zero-gravity (space) configuration.
    pub fn zero_gravity() -> Self {
        Self {
            gravity: [0.0; 3],
            ..Default::default()
        }
    }
    /// Create a config with custom gravity.
    pub fn with_gravity(gx: f64, gy: f64, gz: f64) -> Self {
        Self {
            gravity: [gx, gy, gz],
            ..Default::default()
        }
    }
}
#[wasm_bindgen]
impl WasmSimulationConfig {
    /// Construct a default config (Earth gravity, 60 Hz).
    #[wasm_bindgen(constructor)]
    pub fn wasm_new() -> Self {
        Self::default()
    }
    /// Earth-standard gravity preset.
    #[wasm_bindgen(js_name = "earth")]
    pub fn earth_js() -> Self {
        Self::earth()
    }
    /// Lunar-gravity preset (1/6 Earth).
    #[wasm_bindgen(js_name = "moon")]
    pub fn moon_js() -> Self {
        Self::moon()
    }
    /// Zero-gravity preset (deep space).
    #[wasm_bindgen(js_name = "zeroGravity")]
    pub fn zero_gravity_js() -> Self {
        Self::zero_gravity()
    }
    /// Construct a config with explicit gravity components.
    #[wasm_bindgen(js_name = "withGravity")]
    pub fn with_gravity_js(gx: f64, gy: f64, gz: f64) -> Self {
        Self::with_gravity(gx, gy, gz)
    }
    /// Get gravity as `[x, y, z]`.
    #[wasm_bindgen(js_name = "getGravity")]
    pub fn get_gravity_js(&self) -> Vec<f64> {
        self.gravity.to_vec()
    }
    /// Set gravity components.
    #[wasm_bindgen(js_name = "setGravity")]
    pub fn set_gravity_js(&mut self, gx: f64, gy: f64, gz: f64) {
        self.gravity = [gx, gy, gz];
    }
    /// Serialise the config to a `JsValue` object.
    #[wasm_bindgen(js_name = "toJsValue")]
    pub fn to_js_value_js(&self) -> Result<JsValue, JsValue> {
        to_js_value(self)
    }
}
/// A joint connecting two rigid bodies.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmJoint {
    /// Joint handle.
    pub handle: WasmJointHandle,
    /// First body.
    pub body_a: WasmRigidBodyHandle,
    /// Second body.
    pub body_b: WasmRigidBodyHandle,
    /// Joint type.
    pub joint_type: WasmJointType,
    /// Anchor in body A local space.
    pub anchor_a: [f64; 3],
    /// Anchor in body B local space.
    pub anchor_b: [f64; 3],
    /// Frame in body A local space (quaternion).
    pub frame_a: [f64; 4],
    /// Frame in body B local space (quaternion).
    pub frame_b: [f64; 4],
    /// Motor target velocity.
    pub motor_velocity: f64,
    /// Motor max force.
    pub motor_max_force: f64,
    /// Spring stiffness.
    pub spring_stiffness: f64,
    /// Spring damping.
    pub spring_damping: f64,
    /// Lower limit (angle or distance).
    pub lower_limit: f64,
    /// Upper limit (angle or distance).
    pub upper_limit: f64,
    /// Whether limits are enabled.
    pub limits_enabled: bool,
}
impl WasmJoint {
    /// Create a ball joint between two bodies.
    pub fn ball(
        handle: WasmJointHandle,
        body_a: WasmRigidBodyHandle,
        body_b: WasmRigidBodyHandle,
        anchor_a: [f64; 3],
        anchor_b: [f64; 3],
    ) -> Self {
        Self {
            handle,
            body_a,
            body_b,
            joint_type: WasmJointType::Ball,
            anchor_a,
            anchor_b,
            frame_a: [0.0, 0.0, 0.0, 1.0],
            frame_b: [0.0, 0.0, 0.0, 1.0],
            motor_velocity: 0.0,
            motor_max_force: f64::INFINITY,
            spring_stiffness: 0.0,
            spring_damping: 0.0,
            lower_limit: f64::NEG_INFINITY,
            upper_limit: f64::INFINITY,
            limits_enabled: false,
        }
    }
    /// Create a distance constraint.
    pub fn distance(
        handle: WasmJointHandle,
        body_a: WasmRigidBodyHandle,
        body_b: WasmRigidBodyHandle,
        target_dist: f64,
    ) -> Self {
        Self {
            handle,
            body_a,
            body_b,
            joint_type: WasmJointType::Distance,
            anchor_a: [0.0; 3],
            anchor_b: [0.0; 3],
            frame_a: [0.0, 0.0, 0.0, 1.0],
            frame_b: [0.0, 0.0, 0.0, 1.0],
            motor_velocity: 0.0,
            motor_max_force: f64::INFINITY,
            spring_stiffness: 0.0,
            spring_damping: 0.0,
            lower_limit: target_dist,
            upper_limit: target_dist,
            limits_enabled: true,
        }
    }
}
/// Handle to a collider inside WasmWorld.
#[wasm_bindgen]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct WasmColliderHandle(pub u32);
#[wasm_bindgen]
impl WasmColliderHandle {
    /// Construct a collider handle from a raw `u32` id.
    #[wasm_bindgen(constructor)]
    pub fn wasm_new(id: u32) -> Self {
        WasmColliderHandle(id)
    }
    /// Return the raw identifier.
    #[wasm_bindgen(getter)]
    pub fn id(&self) -> u32 {
        self.0
    }
}
/// Result of an overlap query (AABB or shape).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmOverlapResult {
    /// List of collider handles that overlap the query shape.
    pub colliders: Vec<WasmColliderHandle>,
    /// Associated rigid body handles.
    pub bodies: Vec<Option<WasmRigidBodyHandle>>,
}
impl WasmOverlapResult {
    /// Create an empty result.
    pub fn empty() -> Self {
        Self {
            colliders: Vec::new(),
            bodies: Vec::new(),
        }
    }
    /// Number of overlapping shapes.
    pub fn count(&self) -> usize {
        self.colliders.len()
    }
}
