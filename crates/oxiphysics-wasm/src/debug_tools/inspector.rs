// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Physics inspector — body / joint / island snapshots.
//!
//! [`WasmPhysicsInspector`] keeps a mock registry of bodies and joints whose
//! state can be inspected as JSON snapshots, ready for transferring to a
//! JavaScript debugger UI.

use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

use super::vec3_from_slice;

/// A body state snapshot returned by the inspector.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BodyInspectResult {
    /// Body handle.
    pub id: u32,
    /// Position `[x, y, z]`.
    pub position: [f64; 3],
    /// Velocity `[vx, vy, vz]`.
    pub velocity: [f64; 3],
    /// Angular velocity `[wx, wy, wz]`.
    pub ang_vel: [f64; 3],
    /// Kinetic energy (J).
    pub kinetic_energy: f64,
    /// Whether sleeping.
    pub sleeping: bool,
    /// Whether static.
    pub is_static: bool,
}

/// A joint state snapshot.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct JointInspectResult {
    /// Joint handle.
    pub id: u32,
    /// Body A handle.
    pub body_a: u32,
    /// Body B handle.
    pub body_b: u32,
    /// Joint type name.
    pub joint_type: String,
    /// Reaction force `[fx, fy, fz]`.
    pub reaction_force: [f64; 3],
    /// Motor enabled.
    pub motor_enabled: bool,
}

/// Inspects simulation objects and returns JSON-serialisable snapshots.
#[wasm_bindgen]
#[derive(Debug, Clone, Default)]
pub struct WasmPhysicsInspector {
    /// Mock body registry indexed by id.
    bodies: Vec<InspectorBodyRecord>,
    /// Mock joint registry indexed by id.
    joints: Vec<InspectorJointRecord>,
}

/// Internal record describing a registered body inside the inspector.
#[derive(Debug, Clone)]
struct InspectorBodyRecord {
    id: u32,
    position: [f64; 3],
    velocity: [f64; 3],
    ang_vel: [f64; 3],
    sleeping: bool,
    is_static: bool,
}

/// Internal record describing a registered joint inside the inspector.
#[derive(Debug, Clone)]
struct InspectorJointRecord {
    id: u32,
    body_a: u32,
    body_b: u32,
    type_name: String,
}

impl WasmPhysicsInspector {
    /// Create a new inspector.
    pub fn new() -> Self {
        Self::default()
    }

    /// Register a body for inspection.
    pub fn register_body(
        &mut self,
        id: u32,
        position: [f64; 3],
        velocity: [f64; 3],
        ang_vel: [f64; 3],
        sleeping: bool,
        is_static: bool,
    ) {
        self.bodies.retain(|b| b.id != id);
        self.bodies.push(InspectorBodyRecord {
            id,
            position,
            velocity,
            ang_vel,
            sleeping,
            is_static,
        });
    }

    /// Register a joint for inspection.
    pub fn register_joint(
        &mut self,
        id: u32,
        body_a: u32,
        body_b: u32,
        type_name: impl Into<String>,
    ) {
        self.joints.retain(|j| j.id != id);
        self.joints.push(InspectorJointRecord {
            id,
            body_a,
            body_b,
            type_name: type_name.into(),
        });
    }

    /// Inspect a body; returns JSON string or `None` if not found.
    pub fn inspect_body(&self, id: u32) -> Option<String> {
        let body = self.bodies.iter().find(|b| b.id == id)?;
        let v = body.velocity;
        let v2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
        let ke = 0.5 * v2; // assume mass = 1 for mock
        let result = BodyInspectResult {
            id: body.id,
            position: body.position,
            velocity: body.velocity,
            ang_vel: body.ang_vel,
            kinetic_energy: ke,
            sleeping: body.sleeping,
            is_static: body.is_static,
        };
        serde_json::to_string_pretty(&result).ok()
    }

    /// Inspect a joint; returns JSON or `None`.
    pub fn inspect_joint(&self, id: u32) -> Option<String> {
        let joint = self.joints.iter().find(|j| j.id == id)?;
        let result = JointInspectResult {
            id: joint.id,
            body_a: joint.body_a,
            body_b: joint.body_b,
            joint_type: joint.type_name.clone(),
            reaction_force: [0.0; 3],
            motor_enabled: false,
        };
        serde_json::to_string_pretty(&result).ok()
    }

    /// Inspect an island (mock: returns body count JSON).
    pub fn inspect_island(&self, island_id: u32) -> String {
        format!(
            "{{\"island_id\":{},\"body_count\":{}}}",
            island_id,
            self.bodies.len()
        )
    }
}

#[wasm_bindgen]
impl WasmPhysicsInspector {
    /// Create a new inspector (JS-compatible constructor).
    #[wasm_bindgen(constructor)]
    pub fn wasm_new() -> WasmPhysicsInspector {
        WasmPhysicsInspector::new()
    }

    /// Register a body. `position`, `velocity` and `ang_vel` are flat
    /// `[x, y, z]` arrays passed as `Float64Array` from JavaScript.
    ///
    /// # Errors
    ///
    /// Returns an error if any of the input arrays is not of length 3.
    #[wasm_bindgen(js_name = "register_body")]
    pub fn register_body_js(
        &mut self,
        id: u32,
        position: &[f64],
        velocity: &[f64],
        ang_vel: &[f64],
        sleeping: bool,
        is_static: bool,
    ) -> Result<(), JsValue> {
        let pos = vec3_from_slice(position, "position")?;
        let vel = vec3_from_slice(velocity, "velocity")?;
        let ang = vec3_from_slice(ang_vel, "ang_vel")?;
        self.register_body(id, pos, vel, ang, sleeping, is_static);
        Ok(())
    }

    /// Register a joint between two bodies.
    #[wasm_bindgen(js_name = "register_joint")]
    pub fn register_joint_js(&mut self, id: u32, body_a: u32, body_b: u32, type_name: String) {
        self.register_joint(id, body_a, body_b, type_name);
    }

    /// True when a body with the given id has been registered.
    #[wasm_bindgen(js_name = "has_body")]
    pub fn has_body_js(&self, id: u32) -> bool {
        self.bodies.iter().any(|b| b.id == id)
    }

    /// Inspect a body; returns the JSON snapshot, or an empty string if the
    /// body is not registered. Use [`WasmPhysicsInspector::has_body_js`] to
    /// distinguish "missing" from "empty" results.
    #[wasm_bindgen(js_name = "inspect_body")]
    pub fn inspect_body_js(&self, id: u32) -> String {
        self.inspect_body(id).unwrap_or_default()
    }

    /// True when a joint with the given id has been registered.
    #[wasm_bindgen(js_name = "has_joint")]
    pub fn has_joint_js(&self, id: u32) -> bool {
        self.joints.iter().any(|j| j.id == id)
    }

    /// Inspect a joint; returns the JSON snapshot, or an empty string if the
    /// joint is not registered.
    #[wasm_bindgen(js_name = "inspect_joint")]
    pub fn inspect_joint_js(&self, id: u32) -> String {
        self.inspect_joint(id).unwrap_or_default()
    }

    /// Inspect an island (mock: returns a JSON object with the island id and
    /// total registered body count).
    #[wasm_bindgen(js_name = "inspect_island")]
    pub fn inspect_island_js(&self, island_id: u32) -> String {
        self.inspect_island(island_id)
    }
}
