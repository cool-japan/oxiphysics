//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::engine::WasmPhysicsEngine;
use crate::js_api::engine_collect_contact_events;
use crate::types::{BodyState, ContactResult, DebugInfo, SimulationConfig};
use crate::wasm_helpers::to_js_value;
use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

/// A contact event ready for dispatch to JavaScript event listeners.
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct JsContactEvent {
    /// Type of event: "contact_begin" or "contact_end".
    #[wasm_bindgen(skip)]
    pub event_type: String,
    /// Handle of the first body.
    pub body_a: u32,
    /// Handle of the second body.
    pub body_b: u32,
    /// World-space contact point `[x, y, z]`.
    #[wasm_bindgen(skip)]
    pub contact_point: [f64; 3],
    /// Contact normal `[nx, ny, nz]`.
    #[wasm_bindgen(skip)]
    pub normal: [f64; 3],
    /// Penetration depth.
    pub depth: f64,
    /// Impulse magnitude applied.
    pub impulse: f64,
}
impl JsContactEvent {
    /// Create a "contact_begin" event from a `ContactResult`.
    pub fn begin(cr: &ContactResult) -> Self {
        Self {
            event_type: "contact_begin".to_string(),
            body_a: cr.body_a,
            body_b: cr.body_b,
            contact_point: cr.contact_midpoint(),
            normal: cr.normal,
            depth: cr.depth,
            impulse: cr.impulse,
        }
    }
    /// Serialize to JSON.
    pub fn to_json(&self) -> String {
        serde_json::to_string(self).unwrap_or_else(|_| "{}".to_string())
    }
}
/// An ordered queue of `JsPhysicsEvent` objects accumulated over simulation frames.
///
/// Intended to be flushed once per render frame and dispatched to JS event listeners
/// via `postMessage`.
#[wasm_bindgen]
#[derive(Debug, Default)]
pub struct JsEventQueue {
    events: Vec<JsPhysicsEvent>,
    max_events: usize,
}
impl JsEventQueue {
    /// Create a new event queue with an optional event cap.
    /// `max_events = 0` means unlimited.
    #[allow(dead_code)]
    pub fn new(max_events: usize) -> Self {
        Self {
            events: Vec::new(),
            max_events,
        }
    }
    /// Push an event onto the queue (drops oldest event if at capacity).
    #[allow(dead_code)]
    pub fn push(&mut self, event: JsPhysicsEvent) {
        if self.max_events > 0 && self.events.len() >= self.max_events {
            self.events.remove(0);
        }
        self.events.push(event);
    }
    /// Drain all events and return them.
    #[allow(dead_code)]
    pub fn drain(&mut self) -> Vec<JsPhysicsEvent> {
        std::mem::take(&mut self.events)
    }
    /// Number of queued events.
    #[allow(dead_code)]
    pub fn len(&self) -> usize {
        self.events.len()
    }
    /// Whether the queue is empty.
    #[allow(dead_code)]
    pub fn is_empty(&self) -> bool {
        self.events.is_empty()
    }
    /// Serialize all queued events to a JSON array without draining.
    #[allow(dead_code)]
    pub fn to_json(&self) -> String {
        serde_json::to_string(&self.events).unwrap_or_else(|_| "[]".to_string())
    }
    /// Drain all events and serialize to JSON in one operation.
    #[allow(dead_code)]
    pub fn drain_to_json(&mut self) -> String {
        let events = self.drain();
        serde_json::to_string(&events).unwrap_or_else(|_| "[]".to_string())
    }
    /// Push a step_complete event automatically.
    #[allow(dead_code)]
    pub fn push_step_complete(&mut self, time: f64) {
        self.push(JsPhysicsEvent::step_complete(time));
    }
    /// Push all contact events from the engine into the queue.
    #[allow(dead_code)]
    pub fn collect_contacts_from_engine(&mut self, engine: &WasmPhysicsEngine) {
        for event in engine_collect_contact_events(engine) {
            self.push(event);
        }
    }
}
#[wasm_bindgen]
impl JsEventQueue {
    /// Create a new event queue (JS constructor; `max_events = 0` means unlimited).
    #[wasm_bindgen(constructor)]
    pub fn new_js(max_events: u32) -> JsEventQueue {
        JsEventQueue::new(max_events as usize)
    }
    /// Push a `JsPhysicsEvent` serialized as a `JsValue` onto the queue (JS-accessible).
    ///
    /// The JS caller should pass a `JsPhysicsEvent` object. Returns an error if
    /// deserialization fails.
    #[wasm_bindgen(js_name = "push")]
    pub fn push_js(&mut self, event_js: JsValue) -> Result<(), JsValue> {
        let event: JsPhysicsEvent = serde_wasm_bindgen::from_value(event_js)
            .map_err(|e| JsValue::from_str(&e.to_string()))?;
        self.push(event);
        Ok(())
    }
    /// Drain all events and return them as a JSON string (JS-accessible).
    #[wasm_bindgen(js_name = "drainToJson")]
    pub fn drain_to_json_js(&mut self) -> String {
        self.drain_to_json()
    }
    /// Serialize queued events to a JSON string without draining (JS-accessible).
    #[wasm_bindgen(js_name = "toJson")]
    pub fn to_json_js(&self) -> String {
        self.to_json()
    }
    /// Number of queued events (JS getter).
    #[wasm_bindgen(getter, js_name = "length")]
    pub fn length_js(&self) -> u32 {
        self.events.len() as u32
    }
    /// Push all contact events from the engine (JS-accessible).
    #[wasm_bindgen(js_name = "collectContactsFromEngine")]
    pub fn collect_contacts_from_engine_js(&mut self, engine: &WasmPhysicsEngine) {
        self.collect_contacts_from_engine(engine);
    }
    /// Push a step_complete event (JS-accessible).
    #[wasm_bindgen(js_name = "pushStepComplete")]
    pub fn push_step_complete_js(&mut self, time: f64) {
        self.push_step_complete(time);
    }
}
/// A descriptor for a single body and its collider, for batch creation.
///
/// This is designed to be passed from JavaScript as a JSON array when
/// setting up a scene with many bodies at once.
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct JsBodyDesc {
    /// Mass in kg (0 = static).
    pub mass: f64,
    /// Initial position `[x, y, z]`.
    #[wasm_bindgen(skip)]
    pub position: [f64; 3],
    /// Initial velocity `[vx, vy, vz]`.
    #[wasm_bindgen(skip)]
    pub velocity: [f64; 3],
    /// Collider shape: "sphere", "box", "plane", "capsule".
    #[wasm_bindgen(skip)]
    pub shape: String,
    /// Radius (for sphere/capsule).
    pub radius: f64,
    /// Half-extents `[hx, hy, hz]` (for box).
    #[wasm_bindgen(skip)]
    pub half_extents: [f64; 3],
    /// Height (for capsule).
    pub height: f64,
    /// Restitution \[0, 1\].
    pub restitution: f64,
    /// Friction coefficient.
    pub friction: f64,
    /// User-supplied tag for identifying the body on the JS side.
    #[wasm_bindgen(skip)]
    pub tag: String,
}
impl JsBodyDesc {
    /// Create a sphere body descriptor.
    pub fn sphere(mass: f64, x: f64, y: f64, z: f64, radius: f64) -> Self {
        Self {
            mass,
            position: [x, y, z],
            shape: "sphere".to_string(),
            radius,
            ..Default::default()
        }
    }
    /// Create a box body descriptor.
    pub fn cuboid(mass: f64, x: f64, y: f64, z: f64, hx: f64, hy: f64, hz: f64) -> Self {
        Self {
            mass,
            position: [x, y, z],
            shape: "box".to_string(),
            half_extents: [hx, hy, hz],
            ..Default::default()
        }
    }
    /// Create a static plane descriptor.
    pub fn static_plane(nx: f64, ny: f64, nz: f64, offset: f64) -> Self {
        Self {
            mass: 0.0,
            position: [0.0, offset, 0.0],
            shape: "plane".to_string(),
            half_extents: [nx, ny, nz],
            ..Default::default()
        }
    }
}
/// A serializable snapshot of the simulation state for web worker messaging.
///
/// This can be posted from a WASM worker to the main thread using `postMessage`,
/// or stored in `localStorage` for session restore.
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct JsSimulationState {
    /// Accumulated simulation time (seconds).
    pub time: f64,
    /// Number of active bodies.
    pub num_bodies: u32,
    /// All active body positions as a flat array `[x0, y0, z0, x1, y1, z1, ...]`.
    #[wasm_bindgen(skip)]
    pub positions: Vec<f64>,
    /// Gravity vector `[gx, gy, gz]`.
    #[wasm_bindgen(skip)]
    pub gravity: [f64; 3],
    /// Number of contacts detected in the last step.
    pub contact_count: u32,
    /// Full body states for all active bodies.
    #[wasm_bindgen(skip)]
    pub body_states: Vec<BodyState>,
    /// Debug information from the last step.
    #[wasm_bindgen(skip)]
    pub debug_info: DebugInfo,
}
impl JsSimulationState {
    /// Serialize to a JSON string.
    pub fn to_json(&self) -> String {
        serde_json::to_string(self).unwrap_or_else(|_| "{}".to_string())
    }
    /// Deserialize from a JSON string.
    pub fn from_json(json: &str) -> Option<Self> {
        serde_json::from_str(json).ok()
    }
}
#[wasm_bindgen]
impl JsSimulationState {
    /// Positions as a flat `Vec<f64>` (JS getter).
    #[wasm_bindgen(getter, js_name = "positions")]
    pub fn positions_js(&self) -> Vec<f64> {
        self.positions.clone()
    }
    /// Gravity as `Vec<f64>` `[gx, gy, gz]` (JS getter).
    #[wasm_bindgen(getter, js_name = "gravity")]
    pub fn gravity_js(&self) -> Vec<f64> {
        self.gravity.to_vec()
    }
    /// Serialize to JSON (JS-accessible).
    #[wasm_bindgen(js_name = "toJson")]
    pub fn to_json_js(&self) -> String {
        self.to_json()
    }
    /// Deserialize from JSON (JS static factory).
    #[wasm_bindgen(js_name = "fromJson")]
    pub fn from_json_js(json: &str) -> Option<JsSimulationState> {
        JsSimulationState::from_json(json)
    }
}
/// A record describing a physics event ready to be dispatched to JS.
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct JsPhysicsEvent {
    /// Event category: "contact", "sleep", "wake", "step_complete".
    #[wasm_bindgen(skip)]
    pub category: String,
    /// Event sub-type (e.g. "begin", "end", "body_slept").
    #[wasm_bindgen(skip)]
    pub sub_type: String,
    /// Primary body handle (may be `u32::MAX` if not applicable).
    pub body_a: u32,
    /// Secondary body handle (may be `u32::MAX` if not applicable).
    pub body_b: u32,
    /// Extra floating-point payload (e.g. impulse magnitude, speed).
    pub value: f64,
    /// Simulation time when the event occurred.
    pub time: f64,
}
impl JsPhysicsEvent {
    /// Create a "step_complete" event.
    #[allow(dead_code)]
    pub fn step_complete(time: f64) -> Self {
        Self {
            category: "step".to_string(),
            sub_type: "step_complete".to_string(),
            body_a: u32::MAX,
            body_b: u32::MAX,
            value: time,
            time,
        }
    }
    /// Create a contact-begin event.
    #[allow(dead_code)]
    pub fn contact_begin(body_a: u32, body_b: u32, impulse: f64, time: f64) -> Self {
        Self {
            category: "contact".to_string(),
            sub_type: "begin".to_string(),
            body_a,
            body_b,
            value: impulse,
            time,
        }
    }
    /// Create a contact-end event.
    #[allow(dead_code)]
    pub fn contact_end(body_a: u32, body_b: u32, time: f64) -> Self {
        Self {
            category: "contact".to_string(),
            sub_type: "end".to_string(),
            body_a,
            body_b,
            value: 0.0,
            time,
        }
    }
    /// Create a body-slept event.
    #[allow(dead_code)]
    pub fn body_slept(handle: u32, time: f64) -> Self {
        Self {
            category: "sleep".to_string(),
            sub_type: "body_slept".to_string(),
            body_a: handle,
            body_b: u32::MAX,
            value: 0.0,
            time,
        }
    }
    /// Create a body-woke event.
    #[allow(dead_code)]
    pub fn body_woke(handle: u32, time: f64) -> Self {
        Self {
            category: "wake".to_string(),
            sub_type: "body_woke".to_string(),
            body_a: handle,
            body_b: u32::MAX,
            value: 0.0,
            time,
        }
    }
    /// Serialize to JSON.
    #[allow(dead_code)]
    pub fn to_json(&self) -> String {
        serde_json::to_string(self).unwrap_or_else(|_| "{}".to_string())
    }
}
/// Memory statistics for the WASM heap (estimated).
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmMemoryStats {
    /// Estimated bytes used by all body states (as u32 for WASM compatibility).
    pub body_state_bytes: u32,
    /// Estimated bytes used by contact data (as u32 for WASM compatibility).
    pub contact_bytes: u32,
    /// Number of active bodies (as u32 for WASM compatibility).
    pub body_count: u32,
    /// Number of active contacts (as u32 for WASM compatibility).
    pub contact_count: u32,
}
/// A simple arena that pre-allocates a contiguous `Vec<f64>` scratch buffer
/// for zero-allocation JS↔WASM data transfers each frame.
#[wasm_bindgen]
#[derive(Debug)]
pub struct WasmScratchBuffer {
    buffer: Vec<f64>,
    capacity: usize,
}
impl WasmScratchBuffer {
    /// Allocate a scratch buffer for `capacity` f64 values.
    #[allow(dead_code)]
    pub fn new(capacity: usize) -> Self {
        Self {
            buffer: Vec::with_capacity(capacity),
            capacity,
        }
    }
    /// Fill the buffer with the current body positions.
    /// Returns a slice into the internal buffer.
    #[allow(dead_code)]
    pub fn fill_positions<'a>(&'a mut self, engine: &WasmPhysicsEngine) -> &'a [f64] {
        self.buffer.clear();
        self.buffer.extend_from_slice(&engine.get_all_positions());
        &self.buffer
    }
    /// Fill the buffer with all transforms (pos + quat, 7 per body).
    #[allow(dead_code)]
    pub fn fill_transforms<'a>(&'a mut self, engine: &WasmPhysicsEngine) -> &'a [f64] {
        self.buffer.clear();
        self.buffer.extend_from_slice(&engine.get_all_transforms());
        &self.buffer
    }
    /// Return the current buffer contents without refilling.
    #[allow(dead_code)]
    pub fn as_slice(&self) -> &[f64] {
        &self.buffer
    }
    /// Current number of values in the buffer.
    #[allow(dead_code)]
    pub fn len(&self) -> usize {
        self.buffer.len()
    }
    /// Whether the buffer is empty.
    #[allow(dead_code)]
    pub fn is_empty(&self) -> bool {
        self.buffer.is_empty()
    }
    /// Capacity of the pre-allocated buffer.
    #[allow(dead_code)]
    pub fn capacity(&self) -> usize {
        self.capacity
    }
    /// Reset the buffer without deallocating.
    #[allow(dead_code)]
    pub fn clear(&mut self) {
        self.buffer.clear();
    }
}
#[wasm_bindgen]
impl WasmScratchBuffer {
    /// Create a scratch buffer (JS constructor).
    #[wasm_bindgen(constructor)]
    pub fn new_js(capacity: u32) -> WasmScratchBuffer {
        WasmScratchBuffer::new(capacity as usize)
    }
    /// Fill with body positions and return a cloned `Vec<f64>` (JS-compatible).
    #[wasm_bindgen(js_name = "fillPositions")]
    pub fn fill_positions_js(&mut self, engine: &WasmPhysicsEngine) -> Vec<f64> {
        self.buffer.clear();
        self.buffer.extend_from_slice(&engine.get_all_positions());
        self.buffer.clone()
    }
    /// Fill with body transforms and return a cloned `Vec<f64>` (JS-compatible).
    #[wasm_bindgen(js_name = "fillTransforms")]
    pub fn fill_transforms_js(&mut self, engine: &WasmPhysicsEngine) -> Vec<f64> {
        self.buffer.clear();
        self.buffer.extend_from_slice(&engine.get_all_transforms());
        self.buffer.clone()
    }
    /// Current buffer contents as a cloned `Vec<f64>` (JS-compatible).
    #[wasm_bindgen(js_name = "toVec")]
    pub fn to_vec_js(&self) -> Vec<f64> {
        self.buffer.clone()
    }
    /// Number of values in the buffer (JS getter).
    #[wasm_bindgen(getter, js_name = "length")]
    pub fn length_js(&self) -> u32 {
        self.buffer.len() as u32
    }
    /// Buffer capacity (JS getter).
    #[wasm_bindgen(getter, js_name = "capacity")]
    pub fn capacity_js(&self) -> u32 {
        self.capacity as u32
    }
    /// Clear buffer contents (JS-accessible).
    #[wasm_bindgen(js_name = "clear")]
    pub fn clear_js(&mut self) {
        self.buffer.clear();
    }
}
/// Extended serialized state including velocity histograms and per-body energy.
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct JsExtendedState {
    /// Base simulation state (skipped for wasm_bindgen — access via `baseJson()`).
    #[wasm_bindgen(skip)]
    pub base: JsSimulationState,
    /// Per-body kinetic energy (skipped for wasm_bindgen — access via `kineticEnergies()`).
    #[wasm_bindgen(skip)]
    pub kinetic_energies: Vec<f64>,
    /// Bounding box (skipped for wasm_bindgen — access via `boundingBox()`).
    #[wasm_bindgen(skip)]
    pub bounding_box: [f64; 6],
    /// Center of mass (skipped for wasm_bindgen — access via `centerOfMass()`).
    #[wasm_bindgen(skip)]
    pub center_of_mass: [f64; 3],
    /// Total number of contacts.
    pub contact_count: u32,
}
impl JsExtendedState {
    /// Serialize to JSON.
    #[allow(dead_code)]
    pub fn to_json(&self) -> String {
        serde_json::to_string(self).unwrap_or_else(|_| "{}".to_string())
    }
    /// Deserialize from JSON.
    #[allow(dead_code)]
    pub fn from_json(json: &str) -> Option<Self> {
        serde_json::from_str(json).ok()
    }
}
#[wasm_bindgen]
impl JsExtendedState {
    /// Serialize to JSON (JS-accessible).
    #[wasm_bindgen(js_name = "toJson")]
    pub fn to_json_js(&self) -> String {
        self.to_json()
    }
    /// Deserialize from JSON (JS static factory).
    #[wasm_bindgen(js_name = "fromJson")]
    pub fn from_json_js(json: &str) -> Option<JsExtendedState> {
        JsExtendedState::from_json(json)
    }
    /// Kinetic energies as `Vec<f64>` (JS getter).
    #[wasm_bindgen(getter, js_name = "kineticEnergies")]
    pub fn kinetic_energies_js(&self) -> Vec<f64> {
        self.kinetic_energies.clone()
    }
    /// Bounding box as `Vec<f64>` `[minX, minY, minZ, maxX, maxY, maxZ]` (JS getter).
    #[wasm_bindgen(getter, js_name = "boundingBox")]
    pub fn bounding_box_js(&self) -> Vec<f64> {
        self.bounding_box.to_vec()
    }
    /// Center of mass as `Vec<f64>` `[x, y, z]` (JS getter).
    #[wasm_bindgen(getter, js_name = "centerOfMass")]
    pub fn center_of_mass_js(&self) -> Vec<f64> {
        self.center_of_mass.to_vec()
    }
    /// Full state as a `JsValue` object.
    #[wasm_bindgen(js_name = "toJsValue")]
    pub fn to_js_value_js(&self) -> Result<JsValue, JsValue> {
        to_js_value(self)
    }
}
/// Builder for `WasmPhysicsEngine` configuration.
///
/// Serializable as JSON for storage or transfer via `postMessage`.
///
/// # Example
///
/// ```no_run
/// use oxiphysics_wasm::js_api::JsPhysicsConfig;
///
/// let engine = JsPhysicsConfig::new()
///     .with_gravity(0.0, -9.81, 0.0)
///     .with_fixed_dt(1.0 / 60.0)
///     .build();
/// assert_eq!(engine.get_body_count(), 0);
/// ```
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct JsPhysicsConfig {
    /// Gravity X component (m/s²).
    pub gravity_x: f64,
    /// Gravity Y component (m/s²). Default: -9.81.
    pub gravity_y: f64,
    /// Gravity Z component (m/s²).
    pub gravity_z: f64,
    /// Fixed integration time step in seconds. Default: 1/60.
    pub fixed_dt: f64,
    /// Maximum substeps per `step()` call. Default: 4.
    pub max_substeps: u32,
    /// Number of constraint solver iterations. Default: 8.
    pub solver_iterations: u32,
    /// Enable continuous collision detection. Default: false.
    pub ccd_enabled: bool,
    /// Enable body sleeping. Default: true.
    pub sleeping_enabled: bool,
    /// Linear sleep velocity threshold (m/s). Default: 0.01.
    pub linear_sleep_threshold: f64,
    /// Angular sleep velocity threshold (rad/s). Default: 0.01.
    pub angular_sleep_threshold: f64,
}
impl JsPhysicsConfig {
    /// Create a new config with default values.
    pub fn new() -> Self {
        Self::default()
    }
    /// Set the gravity vector.
    pub fn with_gravity(mut self, x: f64, y: f64, z: f64) -> Self {
        self.gravity_x = x;
        self.gravity_y = y;
        self.gravity_z = z;
        self
    }
    /// Set the fixed time step (seconds).
    pub fn with_fixed_dt(mut self, dt: f64) -> Self {
        self.fixed_dt = dt;
        self
    }
    /// Set the maximum substeps per frame.
    pub fn with_max_substeps(mut self, max: u32) -> Self {
        self.max_substeps = max;
        self
    }
    /// Set the solver iteration count.
    pub fn with_solver_iterations(mut self, n: u32) -> Self {
        self.solver_iterations = n;
        self
    }
    /// Enable or disable CCD.
    pub fn with_ccd(mut self, enabled: bool) -> Self {
        self.ccd_enabled = enabled;
        self
    }
    /// Enable or disable body sleeping.
    pub fn with_sleeping(mut self, enabled: bool) -> Self {
        self.sleeping_enabled = enabled;
        self
    }
    /// Convert this config to a `SimulationConfig`.
    pub fn to_sim_config(&self) -> SimulationConfig {
        SimulationConfig {
            gravity: [self.gravity_x, self.gravity_y, self.gravity_z],
            fixed_dt: self.fixed_dt,
            max_substeps: self.max_substeps,
            solver_iterations: self.solver_iterations,
            ccd_enabled: self.ccd_enabled,
            sleeping_enabled: self.sleeping_enabled,
            linear_sleep_threshold: self.linear_sleep_threshold,
            angular_sleep_threshold: self.angular_sleep_threshold,
            ..Default::default()
        }
    }
    /// Build a `WasmPhysicsEngine` from this config.
    pub fn build(&self) -> WasmPhysicsEngine {
        WasmPhysicsEngine::from_config(self.to_sim_config())
    }
    /// Serialize this config to a JSON string.
    pub fn to_json(&self) -> String {
        serde_json::to_string(self).unwrap_or_else(|_| "{}".to_string())
    }
    /// Deserialize a `JsPhysicsConfig` from a JSON string.
    pub fn from_json(json: &str) -> Option<Self> {
        serde_json::from_str(json).ok()
    }
}
#[wasm_bindgen]
impl JsPhysicsConfig {
    /// Create a new JsPhysicsConfig with default values (JS constructor).
    #[wasm_bindgen(constructor)]
    pub fn new_js() -> JsPhysicsConfig {
        JsPhysicsConfig::new()
    }
    /// Build a WasmPhysicsEngine from this config (JS-accessible).
    #[wasm_bindgen(js_name = "build")]
    pub fn build_js(&self) -> WasmPhysicsEngine {
        self.build()
    }
    /// Serialize to JSON (JS-accessible).
    #[wasm_bindgen(js_name = "toJson")]
    pub fn to_json_js(&self) -> String {
        self.to_json()
    }
    /// Deserialize from JSON (JS static factory).
    #[wasm_bindgen(js_name = "fromJson")]
    pub fn from_json_js(json: &str) -> Option<JsPhysicsConfig> {
        JsPhysicsConfig::from_json(json)
    }
}
