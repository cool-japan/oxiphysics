//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

use super::functions::{quaternion_multiply, quaternion_normalize};
use super::types_3::{SerializedBody, SerializedCollider, WasmWorld};
use super::types_4::{WasmBodyState, WasmBodyType, WasmColliderHandle, WasmRigidBodyHandle};
use crate::wasm_helpers::to_js_value;

/// Handle to a joint (constraint) inside WasmWorld.
#[wasm_bindgen]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct WasmJointHandle(pub u32);
impl WasmJointHandle {
    /// Create a handle from a raw `u32`.
    pub fn new(id: u32) -> Self {
        WasmJointHandle(id)
    }
    /// Return the raw identifier.
    pub fn raw(&self) -> u32 {
        self.0
    }
    /// Sentinel representing an invalid/null handle.
    pub fn invalid() -> Self {
        WasmJointHandle(u32::MAX)
    }
    /// Returns `true` if this handle is the invalid sentinel.
    pub fn is_invalid(&self) -> bool {
        self.0 == u32::MAX
    }
}
#[wasm_bindgen]
impl WasmJointHandle {
    /// Construct a joint handle from a raw `u32` id.
    #[wasm_bindgen(constructor)]
    pub fn wasm_new(id: u32) -> Self {
        WasmJointHandle(id)
    }
    /// Return the raw identifier.
    #[wasm_bindgen(getter)]
    pub fn id(&self) -> u32 {
        self.0
    }
    /// Sentinel representing an invalid/null handle.
    #[wasm_bindgen(js_name = "invalid")]
    pub fn invalid_js() -> Self {
        WasmJointHandle(u32::MAX)
    }
    /// Returns `true` if this handle is the invalid sentinel.
    #[wasm_bindgen(js_name = "isInvalid")]
    pub fn is_invalid_js(&self) -> bool {
        self.0 == u32::MAX
    }
}
/// Broad-phase algorithm selection.
#[wasm_bindgen]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum BroadPhaseAlgorithm {
    /// Simple axis-aligned bounding box sweep and prune.
    Sap,
    /// BVH-based broad-phase.
    Bvh,
    /// Grid-based spatial hash.
    SpatialHash,
}
/// A contact event with detailed contact information.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmContactEvent {
    /// First collider.
    pub collider_a: WasmColliderHandle,
    /// Second collider.
    pub collider_b: WasmColliderHandle,
    /// Contact normal (from B to A).
    pub normal: [f64; 3],
    /// Contact point in world space.
    pub contact_point: [f64; 3],
    /// Penetration depth (positive = overlap).
    pub depth: f64,
    /// Impulse applied to resolve the contact.
    pub impulse: f64,
    /// Friction impulse applied.
    pub friction_impulse: [f64; 3],
}
/// Result of a single raycast query.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmRaycastResult {
    /// Whether the ray hit anything.
    pub hit: bool,
    /// Hit collider.
    pub collider: Option<WasmColliderHandle>,
    /// Associated rigid body.
    pub body: Option<WasmRigidBodyHandle>,
    /// Hit point in world space.
    pub point: [f64; 3],
    /// Surface normal at hit point.
    pub normal: [f64; 3],
    /// Distance from ray origin to hit.
    pub toi: f64,
    /// Feature ID (sub-shape element hit, e.g. triangle index).
    pub feature_id: u32,
}
impl WasmRaycastResult {
    /// Create a "no hit" result.
    pub fn miss() -> Self {
        Self {
            hit: false,
            collider: None,
            body: None,
            point: [0.0; 3],
            normal: [0.0; 3],
            toi: f64::INFINITY,
            feature_id: u32::MAX,
        }
    }
}
/// Collider shape types.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub enum WasmColliderShape {
    /// Sphere with radius.
    Sphere { radius: f64 },
    /// Box with half-extents.
    Box { half_extents: [f64; 3] },
    /// Capsule with half-height and radius.
    Capsule { half_height: f64, radius: f64 },
    /// Infinite plane: normal (nx, ny, nz) and offset d.
    Plane { normal: [f64; 3], offset: f64 },
    /// Cylinder with half-height and radius.
    Cylinder { half_height: f64, radius: f64 },
    /// Cone with half-height and base radius.
    Cone { half_height: f64, radius: f64 },
    /// Convex hull (list of points).
    ConvexHull { points: Vec<[f64; 3]> },
}
impl WasmColliderShape {
    /// Compute an approximate bounding sphere radius.
    pub fn bounding_radius(&self) -> f64 {
        match self {
            Self::Sphere { radius } => *radius,
            Self::Box { half_extents } => (half_extents[0] * half_extents[0]
                + half_extents[1] * half_extents[1]
                + half_extents[2] * half_extents[2])
                .sqrt(),
            Self::Capsule {
                half_height,
                radius,
            } => half_height + radius,
            Self::Plane { .. } => f64::INFINITY,
            Self::Cylinder {
                half_height,
                radius,
            } => (half_height * half_height + radius * radius).sqrt(),
            Self::Cone {
                half_height,
                radius,
            } => (half_height * half_height + radius * radius).sqrt(),
            Self::ConvexHull { points } => points
                .iter()
                .map(|p| (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt())
                .fold(0.0f64, f64::max),
        }
    }
}
/// Full scene snapshot for save/restore.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmSceneSnapshot {
    /// Simulation time.
    pub time: f64,
    /// Step count.
    pub step_count: u64,
    /// Serialized bodies.
    pub bodies: Vec<SerializedBody>,
    /// Serialized colliders.
    pub colliders: Vec<SerializedCollider>,
    /// Gravity.
    pub gravity: [f64; 3],
    /// Configuration dt.
    pub dt: f64,
}
/// WASM scene serializer: converts WasmWorld state to/from a snapshot.
#[wasm_bindgen]
pub struct WasmSceneSerializer;
impl WasmSceneSerializer {
    /// Serialize the given world to a snapshot.
    pub fn serialize(world: &WasmWorld) -> WasmSceneSnapshot {
        let bodies: Vec<SerializedBody> = world
            .bodies
            .values()
            .map(|b| SerializedBody {
                id: b.handle.0,
                body_type: match b.body_type {
                    WasmBodyType::Dynamic => "dynamic".to_string(),
                    WasmBodyType::Static => "static".to_string(),
                    WasmBodyType::Kinematic => "kinematic".to_string(),
                },
                position: b.position,
                rotation: b.rotation,
                linear_vel: b.linear_vel,
                angular_vel: b.angular_vel,
                mass: b.mass,
                sleeping: b.sleeping,
                user_data: b.user_data,
            })
            .collect();
        let colliders: Vec<SerializedCollider> = world
            .colliders
            .values()
            .map(|c| {
                let (shape_type, shape_params) = match &c.shape {
                    WasmColliderShape::Sphere { radius } => ("sphere".to_string(), vec![*radius]),
                    WasmColliderShape::Box { half_extents } => {
                        ("box".to_string(), half_extents.to_vec())
                    }
                    WasmColliderShape::Capsule {
                        half_height,
                        radius,
                    } => ("capsule".to_string(), vec![*half_height, *radius]),
                    WasmColliderShape::Plane { normal, offset } => (
                        "plane".to_string(),
                        vec![normal[0], normal[1], normal[2], *offset],
                    ),
                    WasmColliderShape::Cylinder {
                        half_height,
                        radius,
                    } => ("cylinder".to_string(), vec![*half_height, *radius]),
                    WasmColliderShape::Cone {
                        half_height,
                        radius,
                    } => ("cone".to_string(), vec![*half_height, *radius]),
                    WasmColliderShape::ConvexHull { .. } => ("convex_hull".to_string(), vec![]),
                };
                SerializedCollider {
                    id: c.handle.0,
                    body_id: c.body.map(|h| h.0),
                    shape_type,
                    shape_params,
                    friction: c.friction,
                    restitution: c.restitution,
                    is_sensor: c.is_sensor,
                }
            })
            .collect();
        WasmSceneSnapshot {
            time: world.time,
            step_count: world.step_count,
            bodies,
            colliders,
            gravity: world.config.gravity,
            dt: world.config.dt,
        }
    }
    /// Build a simple JSON string from the snapshot (no serde dependency).
    pub fn to_json_string(snapshot: &WasmSceneSnapshot) -> String {
        let mut out = String::new();
        out.push_str("{\n");
        out.push_str(&format!("  \"time\": {},\n", snapshot.time));
        out.push_str(&format!("  \"step_count\": {},\n", snapshot.step_count));
        out.push_str(&format!(
            "  \"gravity\": [{},{},{}],\n",
            snapshot.gravity[0], snapshot.gravity[1], snapshot.gravity[2]
        ));
        out.push_str(&format!("  \"dt\": {},\n", snapshot.dt));
        out.push_str(&format!("  \"body_count\": {},\n", snapshot.bodies.len()));
        out.push_str(&format!(
            "  \"collider_count\": {}\n",
            snapshot.colliders.len()
        ));
        out.push('}');
        out
    }
}
#[wasm_bindgen]
impl WasmSceneSerializer {
    /// Serialize the world's state to a `JsValue` snapshot object.
    #[wasm_bindgen(js_name = "serialize")]
    pub fn serialize_js(world: &WasmWorld) -> Result<JsValue, JsValue> {
        let snapshot = Self::serialize(world);
        to_js_value(&snapshot)
    }
    /// Serialize the world's state to a JSON string summary.
    #[wasm_bindgen(js_name = "toJsonString")]
    pub fn to_json_string_js(world: &WasmWorld) -> String {
        let snapshot = Self::serialize(world);
        Self::to_json_string(&snapshot)
    }
}
/// Integration method for rigid body dynamics.
#[wasm_bindgen]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum IntegrationMethod {
    /// Semi-implicit Euler (first-order, fast).
    SemiImplicitEuler,
    /// Runge-Kutta 4 (higher accuracy).
    Rk4,
    /// Velocity Verlet.
    Verlet,
}
/// Internal rigid body data.
#[derive(Debug, Clone)]
pub(super) struct WasmRigidBody {
    pub(super) handle: WasmRigidBodyHandle,
    pub(super) body_type: WasmBodyType,
    pub(super) position: [f64; 3],
    pub(super) rotation: [f64; 4],
    pub(super) linear_vel: [f64; 3],
    pub(super) angular_vel: [f64; 3],
    pub(super) mass: f64,
    pub(super) inv_mass: f64,
    pub(super) inv_inertia: [f64; 3],
    pub(super) force_accum: [f64; 3],
    pub(super) torque_accum: [f64; 3],
    pub(super) linear_damping: f64,
    pub(super) angular_damping: f64,
    pub(super) sleeping: bool,
    pub(super) sleep_time: f64,
    pub(super) gravity_scale: f64,
    pub(super) user_data: u64,
}
impl WasmRigidBody {
    pub(super) fn new_dynamic(handle: WasmRigidBodyHandle, position: [f64; 3], mass: f64) -> Self {
        let inv_mass = if mass > 1e-15 { 1.0 / mass } else { 0.0 };
        let inertia_diag = mass * 0.4;
        let inv_inertia = [if inertia_diag > 1e-15 {
            1.0 / inertia_diag
        } else {
            0.0
        }; 3];
        Self {
            handle,
            body_type: WasmBodyType::Dynamic,
            position,
            rotation: [0.0, 0.0, 0.0, 1.0],
            linear_vel: [0.0; 3],
            angular_vel: [0.0; 3],
            mass,
            inv_mass,
            inv_inertia,
            force_accum: [0.0; 3],
            torque_accum: [0.0; 3],
            linear_damping: 0.0,
            angular_damping: 0.0,
            sleeping: false,
            sleep_time: 0.0,
            gravity_scale: 1.0,
            user_data: 0,
        }
    }
    pub(super) fn new_static(handle: WasmRigidBodyHandle, position: [f64; 3]) -> Self {
        Self {
            handle,
            body_type: WasmBodyType::Static,
            position,
            rotation: [0.0, 0.0, 0.0, 1.0],
            linear_vel: [0.0; 3],
            angular_vel: [0.0; 3],
            mass: 0.0,
            inv_mass: 0.0,
            inv_inertia: [0.0; 3],
            force_accum: [0.0; 3],
            torque_accum: [0.0; 3],
            linear_damping: 0.0,
            angular_damping: 0.0,
            sleeping: true,
            sleep_time: f64::INFINITY,
            gravity_scale: 0.0,
            user_data: 0,
        }
    }
    pub(super) fn integrate(&mut self, dt: f64, gravity: [f64; 3]) {
        if self.body_type != WasmBodyType::Dynamic || self.sleeping {
            self.force_accum = [0.0; 3];
            self.torque_accum = [0.0; 3];
            return;
        }
        let grav_force = [
            gravity[0] * self.mass * self.gravity_scale,
            gravity[1] * self.mass * self.gravity_scale,
            gravity[2] * self.mass * self.gravity_scale,
        ];
        let total_force = [
            self.force_accum[0] + grav_force[0],
            self.force_accum[1] + grav_force[1],
            self.force_accum[2] + grav_force[2],
        ];
        let accel = [
            total_force[0] * self.inv_mass,
            total_force[1] * self.inv_mass,
            total_force[2] * self.inv_mass,
        ];
        self.linear_vel[0] += accel[0] * dt;
        self.linear_vel[1] += accel[1] * dt;
        self.linear_vel[2] += accel[2] * dt;
        let lin_damp = (1.0 - self.linear_damping * dt).max(0.0);
        self.linear_vel[0] *= lin_damp;
        self.linear_vel[1] *= lin_damp;
        self.linear_vel[2] *= lin_damp;
        let speed2 = self.linear_vel[0] * self.linear_vel[0]
            + self.linear_vel[1] * self.linear_vel[1]
            + self.linear_vel[2] * self.linear_vel[2];
        let max_v = 200.0f64;
        if speed2 > max_v * max_v {
            let scale = max_v / speed2.sqrt();
            self.linear_vel[0] *= scale;
            self.linear_vel[1] *= scale;
            self.linear_vel[2] *= scale;
        }
        self.position[0] += self.linear_vel[0] * dt;
        self.position[1] += self.linear_vel[1] * dt;
        self.position[2] += self.linear_vel[2] * dt;
        let alpha = [
            self.torque_accum[0] * self.inv_inertia[0],
            self.torque_accum[1] * self.inv_inertia[1],
            self.torque_accum[2] * self.inv_inertia[2],
        ];
        let ang_damp = (1.0 - self.angular_damping * dt).max(0.0);
        self.angular_vel[0] = (self.angular_vel[0] + alpha[0] * dt) * ang_damp;
        self.angular_vel[1] = (self.angular_vel[1] + alpha[1] * dt) * ang_damp;
        self.angular_vel[2] = (self.angular_vel[2] + alpha[2] * dt) * ang_damp;
        let omega = self.angular_vel;
        let angle = (omega[0] * omega[0] + omega[1] * omega[1] + omega[2] * omega[2]).sqrt();
        if angle > 1e-10 {
            let half_angle = angle * dt * 0.5;
            let sin_h = half_angle.sin();
            let cos_h = half_angle.cos();
            let ax = omega[0] / angle;
            let ay = omega[1] / angle;
            let az = omega[2] / angle;
            let dq = [ax * sin_h, ay * sin_h, az * sin_h, cos_h];
            let q = self.rotation;
            self.rotation = quaternion_multiply(dq, q);
            self.rotation = quaternion_normalize(self.rotation);
        }
        self.force_accum = [0.0; 3];
        self.torque_accum = [0.0; 3];
    }
    pub(super) fn state(&self) -> WasmBodyState {
        WasmBodyState {
            handle: self.handle,
            position: self.position,
            rotation: self.rotation,
            linear_vel: self.linear_vel,
            angular_vel: self.angular_vel,
            sleeping: self.sleeping,
            body_type: self.body_type,
        }
    }
}
