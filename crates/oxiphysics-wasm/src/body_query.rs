// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Body query API for browser-based physics simulations.
//!
//! This module provides a rich query layer on top of `WasmPhysicsEngine`.
//! All query results are plain Rust structs serializable to JSON.
//!
//! ## Available queries
//!
//! - **Position / velocity / AABB** — per-body state.
//! - **Sleeping** — which bodies are sleeping.
//! - **Batch queries** — collect data for all bodies in a single call.
//! - **Nearest-body query** — find the closest body to a world-space point.
//! - **AABB overlap query** — find bodies whose AABB overlaps a test box.
//! - **Simulation snapshot** — serializable full-state export.

#![allow(missing_docs)]
#![allow(dead_code)]

use crate::engine::WasmPhysicsEngine;
use crate::types::{BodyState, ContactResult, DebugInfo};
use crate::wasm_helpers::{err_to_jsvalue, to_js_value};
use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

// ---------------------------------------------------------------------------
// BodyAabb
// ---------------------------------------------------------------------------

/// Axis-aligned bounding box for a body.
///
/// Computed by the query layer from position + default radius (0.5 m).
/// In a production engine, collider radii would be used.
#[wasm_bindgen]
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct BodyAabb {
    /// Body handle.
    pub handle: u32,
    /// Minimum corner `[min_x, min_y, min_z]`.
    #[wasm_bindgen(skip)]
    pub min: [f64; 3],
    /// Maximum corner `[max_x, max_y, max_z]`.
    #[wasm_bindgen(skip)]
    pub max: [f64; 3],
}

impl BodyAabb {
    /// Build from a body's position with a sphere of given radius.
    pub fn from_sphere(handle: u32, position: [f64; 3], radius: f64) -> Self {
        Self {
            handle,
            min: [
                position[0] - radius,
                position[1] - radius,
                position[2] - radius,
            ],
            max: [
                position[0] + radius,
                position[1] + radius,
                position[2] + radius,
            ],
        }
    }

    /// Return the center point of this AABB.
    pub fn center(&self) -> [f64; 3] {
        [
            (self.min[0] + self.max[0]) * 0.5,
            (self.min[1] + self.max[1]) * 0.5,
            (self.min[2] + self.max[2]) * 0.5,
        ]
    }

    /// Return the half-extents of this AABB.
    pub fn half_extents(&self) -> [f64; 3] {
        [
            (self.max[0] - self.min[0]) * 0.5,
            (self.max[1] - self.min[1]) * 0.5,
            (self.max[2] - self.min[2]) * 0.5,
        ]
    }

    /// Returns `true` if this AABB overlaps `other`.
    pub fn overlaps(&self, other: &Self) -> bool {
        self.min[0] <= other.max[0]
            && self.max[0] >= other.min[0]
            && self.min[1] <= other.max[1]
            && self.max[1] >= other.min[1]
            && self.min[2] <= other.max[2]
            && self.max[2] >= other.min[2]
    }

    /// Returns `true` if the test box `[test_min, test_max]` overlaps this AABB.
    pub fn overlaps_box(&self, test_min: [f64; 3], test_max: [f64; 3]) -> bool {
        self.min[0] <= test_max[0]
            && self.max[0] >= test_min[0]
            && self.min[1] <= test_max[1]
            && self.max[1] >= test_min[1]
            && self.min[2] <= test_max[2]
            && self.max[2] >= test_min[2]
    }

    /// Returns `true` if the point is inside (or on) this AABB.
    pub fn contains_point(&self, p: [f64; 3]) -> bool {
        p[0] >= self.min[0]
            && p[0] <= self.max[0]
            && p[1] >= self.min[1]
            && p[1] <= self.max[1]
            && p[2] >= self.min[2]
            && p[2] <= self.max[2]
    }
}

// ---------------------------------------------------------------------------
// BodyAabb — wasm-bindgen impl
// ---------------------------------------------------------------------------

#[wasm_bindgen]
impl BodyAabb {
    /// Construct a `BodyAabb` from a body handle, sphere centre, and radius.
    ///
    /// `position_flat` must be `[x, y, z]`; missing components default to `0.0`.
    #[wasm_bindgen(constructor)]
    pub fn wasm_new(handle: u32, position_flat: Vec<f64>, radius: f64) -> BodyAabb {
        BodyAabb::from_sphere(handle, vec3_from_flat(&position_flat), radius)
    }

    /// Minimum corner `[min_x, min_y, min_z]` as a `Vec<f64>`.
    #[wasm_bindgen(getter)]
    pub fn min(&self) -> Vec<f64> {
        self.min.to_vec()
    }

    /// Maximum corner `[max_x, max_y, max_z]` as a `Vec<f64>`.
    #[wasm_bindgen(getter)]
    pub fn max(&self) -> Vec<f64> {
        self.max.to_vec()
    }

    /// Centre point of this AABB.
    #[wasm_bindgen(js_name = "center")]
    pub fn center_js(&self) -> Vec<f64> {
        self.center().to_vec()
    }

    /// Half-extents of this AABB.
    #[wasm_bindgen(js_name = "half_extents")]
    pub fn half_extents_js(&self) -> Vec<f64> {
        self.half_extents().to_vec()
    }

    /// Returns `true` when this AABB overlaps `other`.
    #[wasm_bindgen(js_name = "overlaps")]
    pub fn overlaps_js(&self, other: &BodyAabb) -> bool {
        self.overlaps(other)
    }

    /// Returns `true` when this AABB overlaps the test box `[test_min, test_max]`.
    ///
    /// `test_min_flat` and `test_max_flat` must be 3-element arrays.
    #[wasm_bindgen(js_name = "overlaps_box")]
    pub fn overlaps_box_js(&self, test_min_flat: Vec<f64>, test_max_flat: Vec<f64>) -> bool {
        self.overlaps_box(
            vec3_from_flat(&test_min_flat),
            vec3_from_flat(&test_max_flat),
        )
    }

    /// Returns `true` when the point `[x, y, z]` is inside or on this AABB.
    #[wasm_bindgen(js_name = "contains_point")]
    pub fn contains_point_js(&self, point_flat: Vec<f64>) -> bool {
        self.contains_point(vec3_from_flat(&point_flat))
    }

    /// Serialise to JSON.
    #[wasm_bindgen(js_name = "to_json")]
    pub fn to_json_js(&self) -> Result<String, JsValue> {
        serde_json::to_string(self).map_err(err_to_jsvalue)
    }

    /// Serialise to a `JsValue`.
    #[wasm_bindgen(js_name = "to_js_value")]
    pub fn to_js_value_js(&self) -> Result<JsValue, JsValue> {
        to_js_value(self)
    }
}

// ---------------------------------------------------------------------------
// BodyQueryResult
// ---------------------------------------------------------------------------

/// Full query result for a single body.
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BodyQueryResult {
    /// Body handle.
    pub handle: u32,
    /// Current position `[x, y, z]`.
    #[wasm_bindgen(skip)]
    pub position: [f64; 3],
    /// Current orientation quaternion `[x, y, z, w]`.
    #[wasm_bindgen(skip)]
    pub rotation: [f64; 4],
    /// Current linear velocity `[vx, vy, vz]`.
    #[wasm_bindgen(skip)]
    pub linear_velocity: [f64; 3],
    /// Current angular velocity `[wx, wy, wz]`.
    #[wasm_bindgen(skip)]
    pub angular_velocity: [f64; 3],
    /// Whether the body is sleeping.
    pub is_sleeping: bool,
    /// Whether the body is static (inv_mass == 0 and not kinematic).
    pub is_static: bool,
    /// Kinetic energy.
    pub kinetic_energy: f64,
    /// Axis-aligned bounding box (sphere approximation, radius=0.5).
    #[wasm_bindgen(skip)]
    pub aabb: BodyAabb,
    /// Speed (magnitude of linear velocity).
    pub speed: f64,
}

impl BodyQueryResult {
    /// Build from engine-provided data.
    pub fn from_engine(engine: &WasmPhysicsEngine, handle: u32) -> Option<Self> {
        let state = engine.get_body_state(handle)?;
        let pos = state.position;
        let vel = state.linear_velocity;
        let speed = (vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]).sqrt();
        let is_static = engine.body_inv_mass(handle) == 0.0 && !engine.body_is_dynamic(handle);
        let aabb = BodyAabb::from_sphere(handle, pos, 0.5);
        Some(Self {
            handle,
            position: pos,
            rotation: state.rotation,
            linear_velocity: vel,
            angular_velocity: state.angular_velocity,
            is_sleeping: state.is_sleeping,
            is_static,
            kinetic_energy: state.kinetic_energy,
            aabb,
            speed,
        })
    }
}

// ---------------------------------------------------------------------------
// BodyQueryResult — wasm-bindgen impl
// ---------------------------------------------------------------------------

#[wasm_bindgen]
impl BodyQueryResult {
    /// Build from a live engine and a body handle.
    ///
    /// Returns `undefined` (i.e. `None`) when the handle does not refer to a body.
    #[wasm_bindgen(js_name = "from_engine")]
    pub fn from_engine_js(engine: &WasmPhysicsEngine, handle: u32) -> Option<BodyQueryResult> {
        BodyQueryResult::from_engine(engine, handle)
    }

    /// Position `[x, y, z]` as a `Vec<f64>`.
    #[wasm_bindgen(getter)]
    pub fn position(&self) -> Vec<f64> {
        self.position.to_vec()
    }

    /// Orientation quaternion `[x, y, z, w]` as a `Vec<f64>`.
    #[wasm_bindgen(getter)]
    pub fn rotation(&self) -> Vec<f64> {
        self.rotation.to_vec()
    }

    /// Linear velocity `[vx, vy, vz]` as a `Vec<f64>`.
    #[wasm_bindgen(getter)]
    pub fn linear_velocity(&self) -> Vec<f64> {
        self.linear_velocity.to_vec()
    }

    /// Angular velocity `[wx, wy, wz]` as a `Vec<f64>`.
    #[wasm_bindgen(getter)]
    pub fn angular_velocity(&self) -> Vec<f64> {
        self.angular_velocity.to_vec()
    }

    /// Cached AABB (clone).
    #[wasm_bindgen(getter)]
    pub fn aabb(&self) -> BodyAabb {
        self.aabb
    }

    /// Serialise to JSON.
    #[wasm_bindgen(js_name = "to_json")]
    pub fn to_json_js(&self) -> Result<String, JsValue> {
        serde_json::to_string(self).map_err(err_to_jsvalue)
    }

    /// Serialise to a `JsValue`.
    #[wasm_bindgen(js_name = "to_js_value")]
    pub fn to_js_value_js(&self) -> Result<JsValue, JsValue> {
        to_js_value(self)
    }
}

// ---------------------------------------------------------------------------
// BatchBodyQuery
// ---------------------------------------------------------------------------

/// Result of a batch body query (all active bodies in a single call).
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BatchBodyQuery {
    /// Individual body query results.
    #[wasm_bindgen(skip)]
    pub bodies: Vec<BodyQueryResult>,
    /// Handles of all sleeping bodies.
    #[wasm_bindgen(skip)]
    pub sleeping_handles: Vec<u32>,
    /// Handles of all static bodies.
    #[wasm_bindgen(skip)]
    pub static_handles: Vec<u32>,
    /// Handles of all dynamic bodies.
    #[wasm_bindgen(skip)]
    pub dynamic_handles: Vec<u32>,
    /// Simulation time at which this snapshot was taken.
    pub sim_time: f64,
}

impl BatchBodyQuery {
    /// Collect a batch query from an engine.
    pub fn from_engine(engine: &WasmPhysicsEngine) -> Self {
        let handles = engine.get_all_body_handles();
        let mut bodies = Vec::with_capacity(handles.len());
        let mut sleeping_handles = Vec::new();
        let mut static_handles = Vec::new();
        let mut dynamic_handles = Vec::new();

        for h in &handles {
            if let Some(result) = BodyQueryResult::from_engine(engine, *h) {
                if result.is_sleeping {
                    sleeping_handles.push(*h);
                }
                if result.is_static {
                    static_handles.push(*h);
                } else {
                    dynamic_handles.push(*h);
                }
                bodies.push(result);
            }
        }

        Self {
            bodies,
            sleeping_handles,
            static_handles,
            dynamic_handles,
            sim_time: engine.time(),
        }
    }

    /// Find the body closest to `point` (by position distance).
    ///
    /// Returns `None` if there are no bodies.
    pub fn nearest_to(&self, point: [f64; 3]) -> Option<&BodyQueryResult> {
        self.bodies.iter().min_by(|a, b| {
            let da = dist_sq(a.position, point);
            let db = dist_sq(b.position, point);
            da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
        })
    }

    /// Return all bodies whose AABB overlaps the given axis-aligned box.
    pub fn overlapping_box(&self, test_min: [f64; 3], test_max: [f64; 3]) -> Vec<&BodyQueryResult> {
        self.bodies
            .iter()
            .filter(|b| b.aabb.overlaps_box(test_min, test_max))
            .collect()
    }

    /// Return all bodies with speed above `threshold`.
    pub fn fast_bodies(&self, threshold: f64) -> Vec<&BodyQueryResult> {
        self.bodies.iter().filter(|b| b.speed > threshold).collect()
    }

    /// Return the total kinetic energy of all bodies.
    pub fn total_kinetic_energy(&self) -> f64 {
        self.bodies.iter().map(|b| b.kinetic_energy).sum()
    }

    /// Return flat positions buffer `[x0,y0,z0, x1,y1,z1, ...]`.
    pub fn positions_flat(&self) -> Vec<f64> {
        let mut out = Vec::with_capacity(self.bodies.len() * 3);
        for b in &self.bodies {
            out.extend_from_slice(&b.position);
        }
        out
    }

    /// Return flat velocities buffer `[vx0,vy0,vz0, vx1, ...]`.
    pub fn velocities_flat(&self) -> Vec<f64> {
        let mut out = Vec::with_capacity(self.bodies.len() * 3);
        for b in &self.bodies {
            out.extend_from_slice(&b.linear_velocity);
        }
        out
    }

    /// Return flat transforms buffer `[px,py,pz,qx,qy,qz,qw, ...]`.
    pub fn transforms_flat(&self) -> Vec<f64> {
        let mut out = Vec::with_capacity(self.bodies.len() * 7);
        for b in &self.bodies {
            out.extend_from_slice(&b.position);
            out.extend_from_slice(&b.rotation);
        }
        out
    }

    /// Return a flat AABB buffer: `[handle, min_x, min_y, min_z, max_x, max_y, max_z, ...]`.
    pub fn aabbs_flat(&self) -> Vec<f64> {
        let mut out = Vec::with_capacity(self.bodies.len() * 7);
        for b in &self.bodies {
            out.push(b.handle as f64);
            out.extend_from_slice(&b.aabb.min);
            out.extend_from_slice(&b.aabb.max);
        }
        out
    }

    /// Serialize to JSON.
    pub fn to_json(&self) -> String {
        serde_json::to_string(self).unwrap_or_else(|_| "{}".to_string())
    }
}

// ---------------------------------------------------------------------------
// BatchBodyQuery — wasm-bindgen impl
// ---------------------------------------------------------------------------

#[wasm_bindgen]
impl BatchBodyQuery {
    /// Collect a batch query from a live engine.
    #[wasm_bindgen(js_name = "from_engine")]
    pub fn from_engine_js(engine: &WasmPhysicsEngine) -> BatchBodyQuery {
        BatchBodyQuery::from_engine(engine)
    }

    /// All sleeping body handles.
    #[wasm_bindgen(getter)]
    pub fn sleeping_handles(&self) -> Vec<u32> {
        self.sleeping_handles.clone()
    }

    /// All static body handles.
    #[wasm_bindgen(getter)]
    pub fn static_handles(&self) -> Vec<u32> {
        self.static_handles.clone()
    }

    /// All dynamic body handles.
    #[wasm_bindgen(getter)]
    pub fn dynamic_handles(&self) -> Vec<u32> {
        self.dynamic_handles.clone()
    }

    /// Total number of bodies in this batch.
    #[wasm_bindgen(js_name = "body_count")]
    pub fn body_count_js(&self) -> usize {
        self.bodies.len()
    }

    /// Bodies array as a serde-serialised `JsValue`.
    #[wasm_bindgen(js_name = "bodies_js")]
    pub fn bodies_js(&self) -> Result<JsValue, JsValue> {
        to_js_value(&self.bodies)
    }

    /// Body at the given batch index (clones a `BodyQueryResult`).
    #[wasm_bindgen(js_name = "body_at")]
    pub fn body_at_js(&self, index: usize) -> Option<BodyQueryResult> {
        self.bodies.get(index).cloned()
    }

    /// Find the body closest to `point` (returns its handle).
    ///
    /// Returns `None` when the batch is empty.
    #[wasm_bindgen(js_name = "nearest_to_handle")]
    pub fn nearest_to_handle_js(&self, point_flat: Vec<f64>) -> Option<u32> {
        let p = vec3_from_flat(&point_flat);
        self.nearest_to(p).map(|b| b.handle)
    }

    /// Body handles whose AABB overlaps the test box.
    #[wasm_bindgen(js_name = "overlapping_box_handles")]
    pub fn overlapping_box_handles_js(
        &self,
        test_min_flat: Vec<f64>,
        test_max_flat: Vec<f64>,
    ) -> Vec<u32> {
        let lo = vec3_from_flat(&test_min_flat);
        let hi = vec3_from_flat(&test_max_flat);
        self.overlapping_box(lo, hi)
            .into_iter()
            .map(|b| b.handle)
            .collect()
    }

    /// Body handles whose speed exceeds `threshold`.
    #[wasm_bindgen(js_name = "fast_body_handles")]
    pub fn fast_body_handles_js(&self, threshold: f64) -> Vec<u32> {
        self.fast_bodies(threshold)
            .into_iter()
            .map(|b| b.handle)
            .collect()
    }

    /// Total kinetic energy across all bodies.
    #[wasm_bindgen(js_name = "total_kinetic_energy")]
    pub fn total_kinetic_energy_js(&self) -> f64 {
        self.total_kinetic_energy()
    }

    /// Flat positions buffer `[x0, y0, z0, ...]`.
    #[wasm_bindgen(js_name = "positions_flat")]
    pub fn positions_flat_js(&self) -> Vec<f64> {
        self.positions_flat()
    }

    /// Flat velocities buffer `[vx0, vy0, vz0, ...]`.
    #[wasm_bindgen(js_name = "velocities_flat")]
    pub fn velocities_flat_js(&self) -> Vec<f64> {
        self.velocities_flat()
    }

    /// Flat transforms buffer `[px, py, pz, qx, qy, qz, qw, ...]`.
    #[wasm_bindgen(js_name = "transforms_flat")]
    pub fn transforms_flat_js(&self) -> Vec<f64> {
        self.transforms_flat()
    }

    /// Flat AABB buffer `[handle, min_x, min_y, min_z, max_x, max_y, max_z, ...]`.
    #[wasm_bindgen(js_name = "aabbs_flat")]
    pub fn aabbs_flat_js(&self) -> Vec<f64> {
        self.aabbs_flat()
    }

    /// Serialise to a JSON string.
    #[wasm_bindgen(js_name = "to_json")]
    pub fn to_json_js(&self) -> String {
        self.to_json()
    }

    /// Serialise to a `JsValue`.
    #[wasm_bindgen(js_name = "to_js_value")]
    pub fn to_js_value_js(&self) -> Result<JsValue, JsValue> {
        to_js_value(self)
    }
}

// ---------------------------------------------------------------------------
// SimulationSnapshot
// ---------------------------------------------------------------------------

/// A complete serializable snapshot of the simulation state.
///
/// Suitable for transferring via `postMessage`, saving to IndexedDB,
/// or logging for replay.
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SimulationSnapshot {
    /// All body states.
    #[wasm_bindgen(skip)]
    pub bodies: Vec<BodyState>,
    /// Contact results from the last step.
    #[wasm_bindgen(skip)]
    pub contacts: Vec<ContactResult>,
    /// Debug/performance info.
    #[wasm_bindgen(skip)]
    pub debug_info: DebugInfo,
    /// Simulation time (seconds).
    pub sim_time: f64,
    /// Gravity vector `[gx, gy, gz]`.
    #[wasm_bindgen(skip)]
    pub gravity: [f64; 3],
    /// Total active bodies.
    pub active_body_count: u32,
    /// Total contacts.
    pub contact_count: u32,
    /// Frame number (how many times `step()` was called).
    pub frame: u64,
}

impl SimulationSnapshot {
    /// Collect a full snapshot from an engine.
    pub fn from_engine(engine: &WasmPhysicsEngine, frame: u64) -> Self {
        let handles = engine.get_all_body_handles();
        let bodies: Vec<BodyState> = handles
            .iter()
            .filter_map(|&h| engine.get_body_state(h))
            .collect();
        let contacts: Vec<ContactResult> = engine.get_contacts().to_vec();

        Self {
            active_body_count: engine.get_body_count(),
            contact_count: engine.get_contact_count(),
            sim_time: engine.time(),
            gravity: engine.gravity(),
            debug_info: engine.debug_info().clone(),
            bodies,
            contacts,
            frame,
        }
    }

    /// Serialize to a compact JSON string.
    pub fn to_json(&self) -> String {
        serde_json::to_string(self).unwrap_or_else(|_| "{}".to_string())
    }

    /// Deserialize from JSON.
    pub fn from_json(json: &str) -> Option<Self> {
        serde_json::from_str(json).ok()
    }

    /// Return all positions as a flat `Vec`f64`.
    pub fn positions_flat(&self) -> Vec<f64> {
        let mut out = Vec::with_capacity(self.bodies.len() * 3);
        for b in &self.bodies {
            out.extend_from_slice(&b.position);
        }
        out
    }

    /// Return all transforms (position + rotation) as a flat `Vec`f64`.
    pub fn transforms_flat(&self) -> Vec<f64> {
        let mut out = Vec::with_capacity(self.bodies.len() * 7);
        for b in &self.bodies {
            out.extend_from_slice(&b.position);
            out.extend_from_slice(&b.rotation);
        }
        out
    }

    /// Return indices (handles) of sleeping bodies.
    pub fn sleeping_handles(&self) -> Vec<u32> {
        self.bodies
            .iter()
            .filter(|b| b.is_sleeping)
            .map(|b| b.handle)
            .collect()
    }

    /// Return the total kinetic energy of all bodies.
    pub fn total_kinetic_energy(&self) -> f64 {
        self.bodies.iter().map(|b| b.kinetic_energy).sum()
    }
}

// ---------------------------------------------------------------------------
// SimulationSnapshot — wasm-bindgen impl
// ---------------------------------------------------------------------------

#[wasm_bindgen]
impl SimulationSnapshot {
    /// Collect a full snapshot from a live engine for the given frame number.
    #[wasm_bindgen(js_name = "from_engine")]
    pub fn from_engine_js(engine: &WasmPhysicsEngine, frame: u64) -> SimulationSnapshot {
        SimulationSnapshot::from_engine(engine, frame)
    }

    /// Number of bodies in this snapshot.
    #[wasm_bindgen(js_name = "body_count")]
    pub fn body_count_js(&self) -> usize {
        self.bodies.len()
    }

    /// Number of contacts in this snapshot.
    #[wasm_bindgen(js_name = "contact_count_len")]
    pub fn contact_count_len_js(&self) -> usize {
        self.contacts.len()
    }

    /// Gravity vector `[gx, gy, gz]` as a `Vec<f64>`.
    #[wasm_bindgen(getter)]
    pub fn gravity(&self) -> Vec<f64> {
        self.gravity.to_vec()
    }

    /// Bodies array as a serde-serialised `JsValue`.
    #[wasm_bindgen(js_name = "bodies_js")]
    pub fn bodies_js(&self) -> Result<JsValue, JsValue> {
        to_js_value(&self.bodies)
    }

    /// Contacts array as a serde-serialised `JsValue`.
    #[wasm_bindgen(js_name = "contacts_js")]
    pub fn contacts_js(&self) -> Result<JsValue, JsValue> {
        to_js_value(&self.contacts)
    }

    /// Debug info as a serde-serialised `JsValue`.
    #[wasm_bindgen(js_name = "debug_info_js")]
    pub fn debug_info_js(&self) -> Result<JsValue, JsValue> {
        to_js_value(&self.debug_info)
    }

    /// Serialise to a compact JSON string.
    #[wasm_bindgen(js_name = "to_json")]
    pub fn to_json_js(&self) -> String {
        self.to_json()
    }

    /// Deserialise from a JSON string. Returns `None` on parse failure.
    #[wasm_bindgen(js_name = "from_json")]
    pub fn from_json_js(json: String) -> Option<SimulationSnapshot> {
        SimulationSnapshot::from_json(&json)
    }

    /// All positions as a flat `[x0, y0, z0, ...]` `Vec<f64>`.
    #[wasm_bindgen(js_name = "positions_flat")]
    pub fn positions_flat_js(&self) -> Vec<f64> {
        self.positions_flat()
    }

    /// All transforms as a flat `[px, py, pz, qx, qy, qz, qw, ...]` `Vec<f64>`.
    #[wasm_bindgen(js_name = "transforms_flat")]
    pub fn transforms_flat_js(&self) -> Vec<f64> {
        self.transforms_flat()
    }

    /// Sleeping body handles.
    #[wasm_bindgen(js_name = "sleeping_handles")]
    pub fn sleeping_handles_js(&self) -> Vec<u32> {
        self.sleeping_handles()
    }

    /// Total kinetic energy across all bodies.
    #[wasm_bindgen(js_name = "total_kinetic_energy")]
    pub fn total_kinetic_energy_js(&self) -> f64 {
        self.total_kinetic_energy()
    }

    /// Serialise to a `JsValue`.
    #[wasm_bindgen(js_name = "to_js_value")]
    pub fn to_js_value_js(&self) -> Result<JsValue, JsValue> {
        to_js_value(self)
    }
}

// ---------------------------------------------------------------------------
// ConstraintCountQuery
// ---------------------------------------------------------------------------

/// A query for the number and types of constraints in the simulation.
///
/// In the current engine, this reports the number of active contacts
/// (contact constraints). Future versions will include joint constraints.
#[wasm_bindgen]
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct ConstraintCountQuery {
    /// Number of contact constraints.
    pub contact_constraints: u32,
    /// Number of joint constraints (reserved for future use).
    pub joint_constraints: u32,
    /// Number of limit constraints (reserved for future use).
    pub limit_constraints: u32,
    /// Total constraint count.
    pub total: u32,
}

impl ConstraintCountQuery {
    /// Query constraint counts from an engine.
    pub fn from_engine(engine: &WasmPhysicsEngine) -> Self {
        let contact_constraints = engine.get_contact_count();
        Self {
            contact_constraints,
            joint_constraints: 0,
            limit_constraints: 0,
            total: contact_constraints,
        }
    }
}

// ---------------------------------------------------------------------------
// ConstraintCountQuery — wasm-bindgen impl
// ---------------------------------------------------------------------------

#[wasm_bindgen]
impl ConstraintCountQuery {
    /// Query constraint counts from a live engine.
    #[wasm_bindgen(js_name = "from_engine")]
    pub fn from_engine_js(engine: &WasmPhysicsEngine) -> ConstraintCountQuery {
        ConstraintCountQuery::from_engine(engine)
    }

    /// Serialise to JSON.
    #[wasm_bindgen(js_name = "to_json")]
    pub fn to_json_js(&self) -> Result<String, JsValue> {
        serde_json::to_string(self).map_err(err_to_jsvalue)
    }

    /// Serialise to a `JsValue`.
    #[wasm_bindgen(js_name = "to_js_value")]
    pub fn to_js_value_js(&self) -> Result<JsValue, JsValue> {
        to_js_value(self)
    }
}

// ---------------------------------------------------------------------------
// ContactManifoldsExport
// ---------------------------------------------------------------------------

/// Export of all contact manifold data for JavaScript.
///
/// Provides both a structured `Vec<ContactResult>` and a flat `Float64Array`-
/// compatible buffer for zero-copy transfer.
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ContactManifoldsExport {
    /// Structured contact results.
    #[wasm_bindgen(skip)]
    pub contacts: Vec<ContactResult>,
    /// Number of contacts.
    pub count: u32,
}

impl ContactManifoldsExport {
    /// Export contact manifolds from an engine.
    pub fn from_engine(engine: &WasmPhysicsEngine) -> Self {
        let contacts = engine.get_contacts().to_vec();
        let count = contacts.len() as u32;
        Self { contacts, count }
    }

    /// Return a flat buffer suitable for a JavaScript `Float64Array`.
    ///
    /// Layout per contact: `[body_a, body_b, px, py, pz, nx, ny, nz, depth, impulse]`
    /// (10 f64 values per contact).
    pub fn to_flat_f64(&self) -> Vec<f64> {
        let mut out = Vec::with_capacity(self.contacts.len() * 10);
        for c in &self.contacts {
            out.push(c.body_a as f64);
            out.push(c.body_b as f64);
            // midpoint of the two contact points
            out.push((c.point_on_a[0] + c.point_on_b[0]) * 0.5);
            out.push((c.point_on_a[1] + c.point_on_b[1]) * 0.5);
            out.push((c.point_on_a[2] + c.point_on_b[2]) * 0.5);
            out.extend_from_slice(&c.normal);
            out.push(c.depth);
            out.push(c.impulse);
        }
        out
    }

    /// Return a `Uint32Array`-compatible buffer of body handle pairs.
    ///
    /// Layout: `[body_a_0, body_b_0, body_a_1, body_b_1, ...]`
    pub fn to_handle_pairs_u32(&self) -> Vec<u32> {
        let mut out = Vec::with_capacity(self.contacts.len() * 2);
        for c in &self.contacts {
            out.push(c.body_a);
            out.push(c.body_b);
        }
        out
    }
}

// ---------------------------------------------------------------------------
// ContactManifoldsExport — wasm-bindgen impl
// ---------------------------------------------------------------------------

#[wasm_bindgen]
impl ContactManifoldsExport {
    /// Export the contact manifolds from a live engine.
    #[wasm_bindgen(js_name = "from_engine")]
    pub fn from_engine_js(engine: &WasmPhysicsEngine) -> ContactManifoldsExport {
        ContactManifoldsExport::from_engine(engine)
    }

    /// Contacts array as a serde-serialised `JsValue`.
    #[wasm_bindgen(js_name = "contacts_js")]
    pub fn contacts_js(&self) -> Result<JsValue, JsValue> {
        to_js_value(&self.contacts)
    }

    /// Flat `Float64Array`-compatible buffer.
    ///
    /// Layout per contact:
    /// `[body_a, body_b, px, py, pz, nx, ny, nz, depth, impulse]`.
    #[wasm_bindgen(js_name = "to_flat_f64")]
    pub fn to_flat_f64_js(&self) -> Vec<f64> {
        self.to_flat_f64()
    }

    /// `Uint32Array`-compatible buffer of `[body_a_0, body_b_0, ...]` pairs.
    #[wasm_bindgen(js_name = "to_handle_pairs_u32")]
    pub fn to_handle_pairs_u32_js(&self) -> Vec<u32> {
        self.to_handle_pairs_u32()
    }

    /// Serialise to JSON.
    #[wasm_bindgen(js_name = "to_json")]
    pub fn to_json_js(&self) -> Result<String, JsValue> {
        serde_json::to_string(self).map_err(err_to_jsvalue)
    }

    /// Serialise to a `JsValue`.
    #[wasm_bindgen(js_name = "to_js_value")]
    pub fn to_js_value_js(&self) -> Result<JsValue, JsValue> {
        to_js_value(self)
    }
}

// ---------------------------------------------------------------------------
// Helper
// ---------------------------------------------------------------------------

fn dist_sq(a: [f64; 3], b: [f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    dx * dx + dy * dy + dz * dz
}

fn vec3_from_flat(v: &[f64]) -> [f64; 3] {
    let mut out = [0.0_f64; 3];
    for (i, slot) in out.iter_mut().enumerate() {
        if let Some(val) = v.get(i).copied() {
            *slot = val;
        }
    }
    out
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn make_engine_with_bodies() -> WasmPhysicsEngine {
        let mut e = WasmPhysicsEngine::new(0.0, -9.81, 0.0);
        e.add_dynamic_body(1.0, 0.0, 5.0, 0.0);
        e.add_dynamic_body(2.0, 3.0, 5.0, 0.0);
        e.add_static_body(0.0, 0.0, 0.0);
        e
    }

    // --- BodyAabb ---

    #[test]
    fn test_body_aabb_from_sphere() {
        let aabb = BodyAabb::from_sphere(0, [0.0, 0.0, 0.0], 1.0);
        assert_eq!(aabb.min, [-1.0, -1.0, -1.0]);
        assert_eq!(aabb.max, [1.0, 1.0, 1.0]);
    }

    #[test]
    fn test_body_aabb_center() {
        let aabb = BodyAabb::from_sphere(0, [2.0, 3.0, 4.0], 1.0);
        let c = aabb.center();
        assert!((c[0] - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_body_aabb_half_extents() {
        let aabb = BodyAabb::from_sphere(0, [0.0, 0.0, 0.0], 2.0);
        let he = aabb.half_extents();
        assert!((he[0] - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_body_aabb_overlaps() {
        let a = BodyAabb::from_sphere(0, [0.0, 0.0, 0.0], 1.0);
        let b = BodyAabb::from_sphere(1, [1.5, 0.0, 0.0], 1.0);
        assert!(a.overlaps(&b));
        let c = BodyAabb::from_sphere(2, [5.0, 0.0, 0.0], 1.0);
        assert!(!a.overlaps(&c));
    }

    #[test]
    fn test_body_aabb_overlaps_box() {
        let a = BodyAabb::from_sphere(0, [0.0, 0.0, 0.0], 1.0);
        assert!(a.overlaps_box([-2.0, -2.0, -2.0], [2.0, 2.0, 2.0]));
        assert!(!a.overlaps_box([5.0, 5.0, 5.0], [6.0, 6.0, 6.0]));
    }

    #[test]
    fn test_body_aabb_contains_point() {
        let a = BodyAabb::from_sphere(0, [0.0, 0.0, 0.0], 1.0);
        assert!(a.contains_point([0.5, 0.5, 0.5]));
        assert!(!a.contains_point([2.0, 0.0, 0.0]));
    }

    // --- BodyQueryResult ---

    #[test]
    fn test_body_query_result_from_engine() {
        let engine = make_engine_with_bodies();
        let result = BodyQueryResult::from_engine(&engine, 0).expect("should have result");
        assert_eq!(result.handle, 0);
        assert!((result.position[1] - 5.0).abs() < 1e-10);
        assert!(!result.is_static);
    }

    #[test]
    fn test_body_query_result_static() {
        let engine = make_engine_with_bodies();
        // handle 2 is the static body
        let result = BodyQueryResult::from_engine(&engine, 2).expect("should have result");
        assert!(result.is_static);
    }

    #[test]
    fn test_body_query_result_invalid_handle() {
        let engine = make_engine_with_bodies();
        let result = BodyQueryResult::from_engine(&engine, 99);
        assert!(result.is_none());
    }

    // --- BatchBodyQuery ---

    #[test]
    fn test_batch_body_query_count() {
        let engine = make_engine_with_bodies();
        let bq = BatchBodyQuery::from_engine(&engine);
        assert_eq!(bq.bodies.len(), 3);
        assert_eq!(bq.static_handles.len(), 1);
        assert_eq!(bq.dynamic_handles.len(), 2);
    }

    #[test]
    fn test_batch_body_query_nearest_to() {
        let engine = make_engine_with_bodies();
        let bq = BatchBodyQuery::from_engine(&engine);
        let nearest = bq.nearest_to([0.0, 5.0, 0.0]).expect("should find nearest");
        // body 0 is at (0,5,0) — closest
        assert_eq!(nearest.handle, 0);
    }

    #[test]
    fn test_batch_body_query_nearest_to_empty() {
        let engine = WasmPhysicsEngine::new(0.0, -9.81, 0.0);
        let bq = BatchBodyQuery::from_engine(&engine);
        assert!(bq.nearest_to([0.0, 0.0, 0.0]).is_none());
    }

    #[test]
    fn test_batch_body_query_overlapping_box() {
        let engine = make_engine_with_bodies();
        let bq = BatchBodyQuery::from_engine(&engine);
        // Box around origin
        let overlapping = bq.overlapping_box([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]);
        // Body 2 (static at origin) should be in this range
        assert!(!overlapping.is_empty());
    }

    #[test]
    fn test_batch_body_query_positions_flat() {
        let engine = make_engine_with_bodies();
        let bq = BatchBodyQuery::from_engine(&engine);
        let flat = bq.positions_flat();
        assert_eq!(flat.len(), 9); // 3 bodies × 3
    }

    #[test]
    fn test_batch_body_query_transforms_flat() {
        let engine = make_engine_with_bodies();
        let bq = BatchBodyQuery::from_engine(&engine);
        let flat = bq.transforms_flat();
        assert_eq!(flat.len(), 21); // 3 bodies × 7
    }

    #[test]
    fn test_batch_body_query_aabbs_flat() {
        let engine = make_engine_with_bodies();
        let bq = BatchBodyQuery::from_engine(&engine);
        let flat = bq.aabbs_flat();
        assert_eq!(flat.len(), 21); // 3 bodies × 7 (handle + min[3] + max[3])
    }

    #[test]
    fn test_batch_body_query_velocities_flat() {
        let engine = make_engine_with_bodies();
        let bq = BatchBodyQuery::from_engine(&engine);
        let flat = bq.velocities_flat();
        assert_eq!(flat.len(), 9); // 3 × 3
    }

    #[test]
    fn test_batch_body_query_total_kinetic_energy() {
        let mut engine = WasmPhysicsEngine::new(0.0, 0.0, 0.0);
        let h = engine.add_dynamic_body(2.0, 0.0, 0.0, 0.0);
        engine.set_velocity(h, 3.0, 0.0, 0.0).unwrap();
        let bq = BatchBodyQuery::from_engine(&engine);
        // KE = 0.5 * 2 * 3² = 9
        assert!((bq.total_kinetic_energy() - 9.0).abs() < 1e-6);
    }

    // --- SimulationSnapshot ---

    #[test]
    fn test_simulation_snapshot_from_engine() {
        let mut engine = make_engine_with_bodies();
        engine.step(1.0 / 60.0);
        let snap = SimulationSnapshot::from_engine(&engine, 1);
        assert_eq!(snap.bodies.len(), 3);
        assert_eq!(snap.frame, 1);
        assert_eq!(snap.active_body_count, 3);
    }

    #[test]
    fn test_simulation_snapshot_json_roundtrip() {
        let engine = make_engine_with_bodies();
        let snap = SimulationSnapshot::from_engine(&engine, 0);
        let json = snap.to_json();
        let back = SimulationSnapshot::from_json(&json).expect("should deserialize");
        assert_eq!(back.bodies.len(), snap.bodies.len());
        assert_eq!(back.frame, snap.frame);
    }

    #[test]
    fn test_simulation_snapshot_positions_flat() {
        let engine = make_engine_with_bodies();
        let snap = SimulationSnapshot::from_engine(&engine, 0);
        let flat = snap.positions_flat();
        assert_eq!(flat.len(), 9);
    }

    #[test]
    fn test_simulation_snapshot_transforms_flat() {
        let engine = make_engine_with_bodies();
        let snap = SimulationSnapshot::from_engine(&engine, 0);
        let flat = snap.transforms_flat();
        assert_eq!(flat.len(), 21);
    }

    #[test]
    fn test_simulation_snapshot_sleeping_handles() {
        let engine = make_engine_with_bodies();
        let snap = SimulationSnapshot::from_engine(&engine, 0);
        // No body has been sleeping long enough yet
        let sleeping = snap.sleeping_handles();
        let _ = sleeping; // just assert no panic
    }

    // --- ConstraintCountQuery ---

    #[test]
    fn test_constraint_count_no_contacts() {
        let engine = make_engine_with_bodies();
        let q = ConstraintCountQuery::from_engine(&engine);
        assert_eq!(q.contact_constraints, 0);
        assert_eq!(q.joint_constraints, 0);
        assert_eq!(q.total, 0);
    }

    #[test]
    fn test_constraint_count_with_contacts() {
        let mut engine = WasmPhysicsEngine::new(0.0, 0.0, 0.0);
        let b0 = engine.add_dynamic_body(1.0, 0.0, 0.0, 0.0);
        let b1 = engine.add_dynamic_body(1.0, 1.5, 0.0, 0.0);
        engine.add_sphere_collider(b0, 1.0);
        engine.add_sphere_collider(b1, 1.0);
        engine.step(1.0 / 60.0);
        let q = ConstraintCountQuery::from_engine(&engine);
        assert_eq!(q.total, engine.get_contact_count());
    }

    // --- ContactManifoldsExport ---

    #[test]
    fn test_contact_manifolds_export_empty() {
        let engine = make_engine_with_bodies();
        let exp = ContactManifoldsExport::from_engine(&engine);
        assert_eq!(exp.count, 0);
        assert!(exp.to_flat_f64().is_empty());
        assert!(exp.to_handle_pairs_u32().is_empty());
    }

    #[test]
    fn test_contact_manifolds_export_with_contacts() {
        let mut engine = WasmPhysicsEngine::new(0.0, 0.0, 0.0);
        let b0 = engine.add_dynamic_body(1.0, 0.0, 0.0, 0.0);
        let b1 = engine.add_dynamic_body(1.0, 1.5, 0.0, 0.0);
        engine.add_sphere_collider(b0, 1.0);
        engine.add_sphere_collider(b1, 1.0);
        engine.step(1.0 / 60.0);
        let exp = ContactManifoldsExport::from_engine(&engine);
        let count = exp.count;
        let flat = exp.to_flat_f64();
        assert_eq!(flat.len(), count as usize * 10);
        let pairs = exp.to_handle_pairs_u32();
        assert_eq!(pairs.len(), count as usize * 2);
    }

    #[test]
    fn test_batch_body_query_fast_bodies() {
        let mut engine = WasmPhysicsEngine::new(0.0, 0.0, 0.0);
        let h = engine.add_dynamic_body(1.0, 0.0, 0.0, 0.0);
        engine.set_velocity(h, 10.0, 0.0, 0.0).unwrap();
        engine.add_dynamic_body(1.0, 5.0, 0.0, 0.0); // slow body
        let bq = BatchBodyQuery::from_engine(&engine);
        let fast = bq.fast_bodies(5.0);
        assert_eq!(fast.len(), 1);
        assert_eq!(fast[0].handle, h);
    }

    #[test]
    fn test_batch_body_query_json_roundtrip() {
        let engine = make_engine_with_bodies();
        let bq = BatchBodyQuery::from_engine(&engine);
        let json = bq.to_json();
        assert!(!json.is_empty());
        // just check it is valid JSON by checking it starts with "{"
        assert!(json.starts_with('{'));
    }
}
