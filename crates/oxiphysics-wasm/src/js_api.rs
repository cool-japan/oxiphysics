// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Flat JavaScript-interop API for the OxiPhysics WASM engine.
//!
//! This module provides two layers:
//!
//! 1. **`JsPhysicsConfig`** — a builder for engine configuration.
//! 2. **`JsSimulationState`** — a serializable snapshot for web worker messaging.
//! 3. **`JsBodyDesc`** — a serializable body descriptor for batch creation.
//! 4. **Free functions** — `create_engine`, `destroy_engine`, `engine_step`, etc.
//!    that mirror the wasm-bindgen-style flat API.
//!
//! All types use `serde` for JSON round-tripping and are designed to be
//! transferred through the structured-clone algorithm or `postMessage`.

#![allow(missing_docs)]

use crate::engine::WasmPhysicsEngine;
use crate::types::{
    BodyState, ColliderConfig, ContactResult, DebugInfo, RigidBodyConfig, SimulationConfig,
    Vec3Wasm,
};
use serde::{Deserialize, Serialize};

// ---------------------------------------------------------------------------
// JsPhysicsConfig
// ---------------------------------------------------------------------------

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

impl Default for JsPhysicsConfig {
    fn default() -> Self {
        let sim = SimulationConfig::default();
        Self {
            gravity_x: sim.gravity[0],
            gravity_y: sim.gravity[1],
            gravity_z: sim.gravity[2],
            fixed_dt: sim.fixed_dt,
            max_substeps: sim.max_substeps,
            solver_iterations: sim.solver_iterations,
            ccd_enabled: sim.ccd_enabled,
            sleeping_enabled: sim.sleeping_enabled,
            linear_sleep_threshold: sim.linear_sleep_threshold,
            angular_sleep_threshold: sim.angular_sleep_threshold,
        }
    }
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

// ---------------------------------------------------------------------------
// JsSimulationState
// ---------------------------------------------------------------------------

/// A serializable snapshot of the simulation state for web worker messaging.
///
/// This can be posted from a WASM worker to the main thread using `postMessage`,
/// or stored in `localStorage` for session restore.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct JsSimulationState {
    /// Accumulated simulation time (seconds).
    pub time: f64,
    /// Number of active bodies.
    pub num_bodies: u32,
    /// All active body positions as a flat array `[x0, y0, z0, x1, y1, z1, ...]`.
    pub positions: Vec<f64>,
    /// Gravity vector `[gx, gy, gz]`.
    pub gravity: [f64; 3],
    /// Number of contacts detected in the last step.
    pub contact_count: u32,
    /// Full body states for all active bodies.
    pub body_states: Vec<BodyState>,
    /// Debug information from the last step.
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

/// Serialize the current engine state to a JSON string for web worker messaging.
pub fn serialize_state(engine: &WasmPhysicsEngine) -> String {
    let body_states: Vec<BodyState> = engine
        .get_all_body_handles()
        .into_iter()
        .filter_map(|h| engine.get_body_state(h))
        .collect();

    let state = JsSimulationState {
        time: engine.time(),
        num_bodies: engine.get_body_count(),
        positions: engine.get_all_positions(),
        gravity: engine.gravity(),
        contact_count: engine.get_contact_count(),
        body_states,
        debug_info: engine.debug_info().clone(),
    };
    serde_json::to_string(&state).unwrap_or_else(|_| "{}".to_string())
}

/// Deserialize a `JsSimulationState` from a JSON string.
pub fn deserialize_state(json: &str) -> Option<JsSimulationState> {
    serde_json::from_str(json).ok()
}

// ---------------------------------------------------------------------------
// JsBodyDesc
// ---------------------------------------------------------------------------

/// A descriptor for a single body and its collider, for batch creation.
///
/// This is designed to be passed from JavaScript as a JSON array when
/// setting up a scene with many bodies at once.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct JsBodyDesc {
    /// Mass in kg (0 = static).
    pub mass: f64,
    /// Initial position `[x, y, z]`.
    pub position: [f64; 3],
    /// Initial velocity `[vx, vy, vz]`.
    pub velocity: [f64; 3],
    /// Collider shape: "sphere", "box", "plane", "capsule".
    pub shape: String,
    /// Radius (for sphere/capsule).
    pub radius: f64,
    /// Half-extents `[hx, hy, hz]` (for box).
    pub half_extents: [f64; 3],
    /// Height (for capsule).
    pub height: f64,
    /// Restitution \[0, 1\].
    pub restitution: f64,
    /// Friction coefficient.
    pub friction: f64,
    /// User-supplied tag for identifying the body on the JS side.
    pub tag: String,
}

impl Default for JsBodyDesc {
    fn default() -> Self {
        Self {
            mass: 1.0,
            position: [0.0; 3],
            velocity: [0.0; 3],
            shape: "sphere".to_string(),
            radius: 0.5,
            half_extents: [0.5; 3],
            height: 1.0,
            restitution: 0.3,
            friction: 0.5,
            tag: String::new(),
        }
    }
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
            half_extents: [nx, ny, nz], // re-use field as normal
            ..Default::default()
        }
    }
}

// ---------------------------------------------------------------------------
// JsContactEvent
// ---------------------------------------------------------------------------

/// A contact event ready for dispatch to JavaScript event listeners.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct JsContactEvent {
    /// Type of event: "contact_begin" or "contact_end".
    pub event_type: String,
    /// Handle of the first body.
    pub body_a: u32,
    /// Handle of the second body.
    pub body_b: u32,
    /// World-space contact point `[x, y, z]`.
    pub contact_point: [f64; 3],
    /// Contact normal `[nx, ny, nz]`.
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

// ---------------------------------------------------------------------------
// Flat API functions (future wasm-bindgen targets)
// ---------------------------------------------------------------------------

/// Create a new physics engine with the given gravity and return it.
///
/// This is the primary entry point for JavaScript usage.
pub fn create_engine(gravity_x: f64, gravity_y: f64, gravity_z: f64) -> WasmPhysicsEngine {
    WasmPhysicsEngine::new(gravity_x, gravity_y, gravity_z)
}

/// Create a physics engine from a JSON configuration string.
///
/// Falls back to default gravity `(0, -9.81, 0)` if parsing fails.
pub fn create_engine_from_json(config_json: &str) -> WasmPhysicsEngine {
    if let Some(cfg) = JsPhysicsConfig::from_json(config_json) {
        cfg.build()
    } else {
        WasmPhysicsEngine::new(0.0, -9.81, 0.0)
    }
}

/// Advance a mutable engine reference by `dt` seconds.
///
/// This is the primary simulation loop function; call it every animation frame.
pub fn engine_step(engine: &mut WasmPhysicsEngine, dt: f64) {
    engine.step(dt);
}

/// Add a sphere body to the engine.
///
/// Returns the body handle.
pub fn engine_add_sphere(
    engine: &mut WasmPhysicsEngine,
    mass: f64,
    x: f64,
    y: f64,
    z: f64,
    radius: f64,
) -> u32 {
    let cfg = RigidBodyConfig::new()
        .with_mass(mass)
        .with_position(x, y, z);
    let handle = engine.add_rigid_body(&cfg);
    engine.add_sphere_collider(handle, radius);
    handle
}

/// Add a box body to the engine.
///
/// Returns the body handle.
#[allow(clippy::too_many_arguments)]
pub fn engine_add_box(
    engine: &mut WasmPhysicsEngine,
    mass: f64,
    x: f64,
    y: f64,
    z: f64,
    hx: f64,
    hy: f64,
    hz: f64,
) -> u32 {
    let cfg = RigidBodyConfig::new()
        .with_mass(mass)
        .with_position(x, y, z);
    let handle = engine.add_rigid_body(&cfg);
    engine.add_box_collider(handle, hx, hy, hz);
    handle
}

/// Add a static infinite plane to the engine.
///
/// The plane normal `(nx, ny, nz)` should be normalized.
/// Returns the body handle.
pub fn engine_add_static_plane(
    engine: &mut WasmPhysicsEngine,
    nx: f64,
    ny: f64,
    nz: f64,
    offset: f64,
) -> u32 {
    let cfg = RigidBodyConfig::static_body();
    let handle = engine.add_rigid_body(&cfg);
    engine.add_plane_collider(handle, nx, ny, nz, offset);
    handle
}

/// Change the gravity of a running engine.
pub fn engine_set_gravity(engine: &mut WasmPhysicsEngine, gx: f64, gy: f64, gz: f64) {
    engine.set_gravity(gx, gy, gz);
}

/// Get the position of a body as `[x, y, z]`.
pub fn engine_get_position(engine: &WasmPhysicsEngine, handle: u32) -> [f64; 3] {
    engine.get_position(handle)
}

/// Get the velocity of a body as `[vx, vy, vz]`.
pub fn engine_get_velocity(engine: &WasmPhysicsEngine, handle: u32) -> [f64; 3] {
    engine.get_velocity(handle)
}

/// Get all body positions as a flat `Vec`f64`.
pub fn engine_get_all_positions(engine: &WasmPhysicsEngine) -> Vec<f64> {
    engine.get_all_positions()
}

/// Get all body transforms (position + quaternion) as a flat `Vec`f64`.
pub fn engine_get_all_transforms(engine: &WasmPhysicsEngine) -> Vec<f64> {
    engine.get_all_transforms()
}

/// Get the number of active bodies.
pub fn engine_get_body_count(engine: &WasmPhysicsEngine) -> u32 {
    engine.get_body_count()
}

/// Get the number of contacts from the last step.
pub fn engine_get_contact_count(engine: &WasmPhysicsEngine) -> u32 {
    engine.get_contact_count()
}

/// Reset the engine (remove all bodies, reset time).
pub fn engine_reset(engine: &mut WasmPhysicsEngine) {
    engine.reset();
}

/// Serialize the engine state to JSON.
pub fn engine_serialize(engine: &WasmPhysicsEngine) -> String {
    serialize_state(engine)
}

/// Apply an impulse to a body by handle.
pub fn engine_apply_impulse(
    engine: &mut WasmPhysicsEngine,
    handle: u32,
    ix: f64,
    iy: f64,
    iz: f64,
) -> bool {
    engine.apply_impulse(handle, ix, iy, iz).is_ok()
}

/// Apply a force to a body for the next step.
pub fn engine_apply_force(
    engine: &mut WasmPhysicsEngine,
    handle: u32,
    fx: f64,
    fy: f64,
    fz: f64,
) -> bool {
    engine.apply_force(handle, fx, fy, fz).is_ok()
}

/// Set the linear velocity of a body.
pub fn engine_set_velocity(
    engine: &mut WasmPhysicsEngine,
    handle: u32,
    vx: f64,
    vy: f64,
    vz: f64,
) -> bool {
    engine.set_velocity(handle, vx, vy, vz).is_ok()
}

/// Remove a body from the engine.
pub fn engine_remove_body(engine: &mut WasmPhysicsEngine, handle: u32) -> bool {
    engine.remove_body(handle).is_ok()
}

/// Cast a ray and return the result as a JSON string.
#[allow(clippy::too_many_arguments)]
pub fn engine_raycast_json(
    engine: &WasmPhysicsEngine,
    ox: f64,
    oy: f64,
    oz: f64,
    dx: f64,
    dy: f64,
    dz: f64,
    max_dist: f64,
) -> String {
    let result = engine.raycast(ox, oy, oz, dx, dy, dz, max_dist);
    serde_json::to_string(&result).unwrap_or_else(|_| "{}".to_string())
}

/// Get contact events from the last step as a JSON string.
pub fn engine_get_contact_events_json(engine: &WasmPhysicsEngine) -> String {
    let events: Vec<JsContactEvent> = engine
        .get_contacts()
        .iter()
        .map(JsContactEvent::begin)
        .collect();
    serde_json::to_string(&events).unwrap_or_else(|_| "[]".to_string())
}

/// Batch-add bodies from a JSON array of `JsBodyDesc` objects.
///
/// Returns the number of bodies successfully added.
pub fn engine_add_bodies_from_json(engine: &mut WasmPhysicsEngine, json: &str) -> u32 {
    let descs: Vec<JsBodyDesc> = match serde_json::from_str(json) {
        Ok(v) => v,
        Err(_) => return 0,
    };
    let mut count = 0u32;
    for desc in &descs {
        let cfg = RigidBodyConfig::new()
            .with_mass(desc.mass)
            .with_position(desc.position[0], desc.position[1], desc.position[2])
            .with_linear_velocity(desc.velocity[0], desc.velocity[1], desc.velocity[2])
            .with_restitution(desc.restitution)
            .with_friction(desc.friction);
        let cfg = if desc.mass == 0.0 {
            RigidBodyConfig::static_body()
                .with_position(desc.position[0], desc.position[1], desc.position[2])
                .with_restitution(desc.restitution)
                .with_friction(desc.friction)
        } else {
            cfg
        };
        let handle = engine.add_rigid_body(&cfg);

        match desc.shape.as_str() {
            "sphere" => {
                engine.add_sphere_collider(handle, desc.radius);
            }
            "box" => {
                engine.add_box_collider(
                    handle,
                    desc.half_extents[0],
                    desc.half_extents[1],
                    desc.half_extents[2],
                );
            }
            "capsule" => {
                engine.add_capsule_collider(handle, desc.radius, desc.height);
            }
            "plane" => {
                let n = desc.half_extents; // plane normal stored in half_extents
                let norm = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt().max(1e-15);
                engine.add_plane_collider(handle, n[0] / norm, n[1] / norm, n[2] / norm, 0.0);
            }
            _ => {
                // default: sphere
                engine.add_sphere_collider(handle, desc.radius);
            }
        }
        count += 1;
    }
    count
}

// ---------------------------------------------------------------------------
// Utility: compute a simple hash of the engine state (for dirty-checking)
// ---------------------------------------------------------------------------

/// Compute a lightweight hash of all body positions (for dirty-checking from JS).
///
/// The returned `u64` changes whenever any body moves. It is not
/// cryptographically secure and should only be used for change detection.
pub fn engine_state_hash(engine: &WasmPhysicsEngine) -> u64 {
    let positions = engine.get_all_positions();
    let mut hash: u64 = 0xcbf2_9ce4_8422_2325; // FNV offset basis
    for &v in &positions {
        let bits = v.to_bits();
        hash ^= bits;
        hash = hash.wrapping_mul(0x0000_0100_0000_01b3); // FNV prime
    }
    hash
}

// ---------------------------------------------------------------------------
// Utility: scene setup helpers
// ---------------------------------------------------------------------------

/// Set up a standard demo scene: ground plane + N falling spheres.
///
/// Returns a `Vec`u32` of sphere body handles.
pub fn setup_falling_spheres_scene(
    engine: &mut WasmPhysicsEngine,
    count: u32,
    radius: f64,
    start_height: f64,
) -> Vec<u32> {
    // Ground plane
    engine_add_static_plane(engine, 0.0, 1.0, 0.0, 0.0);

    let mut handles = Vec::with_capacity(count as usize);
    let cols = (count as f64).sqrt().ceil() as u32;
    for i in 0..count {
        let row = i / cols;
        let col = i % cols;
        let x = col as f64 * (radius * 2.5) - (cols as f64 * radius * 1.25);
        let y = start_height + row as f64 * (radius * 3.0);
        let z = 0.0;
        let h = engine_add_sphere(engine, 1.0, x, y, z, radius);
        handles.push(h);
    }
    handles
}

/// Compute the bounding box of all bodies as `\[min_x, min_y, min_z, max_x, max_y, max_z\]`.
pub fn engine_bounding_box(engine: &WasmPhysicsEngine) -> [f64; 6] {
    let positions = engine.get_all_positions();
    if positions.is_empty() {
        return [0.0; 6];
    }
    let mut min = [f64::INFINITY; 3];
    let mut max = [f64::NEG_INFINITY; 3];
    for chunk in positions.chunks(3) {
        for i in 0..3 {
            min[i] = min[i].min(chunk[i]);
            max[i] = max[i].max(chunk[i]);
        }
    }
    [min[0], min[1], min[2], max[0], max[1], max[2]]
}

/// Compute the center of mass of all active dynamic bodies.
pub fn engine_center_of_mass(engine: &WasmPhysicsEngine) -> Vec3Wasm {
    let positions = engine.get_all_positions();
    if positions.is_empty() {
        return Vec3Wasm::zero();
    }
    let n = (positions.len() / 3) as f64;
    let mut sum = Vec3Wasm::zero();
    for chunk in positions.chunks(3) {
        sum.x += chunk[0];
        sum.y += chunk[1];
        sum.z += chunk[2];
    }
    Vec3Wasm::new(sum.x / n, sum.y / n, sum.z / n)
}

/// Add a collider to an existing body, specified via JSON (a `ColliderConfig`).
///
/// Returns the collider handle, or `u32::MAX` on failure.
pub fn engine_add_collider_json(
    engine: &mut WasmPhysicsEngine,
    body_handle: u32,
    collider_json: &str,
) -> u32 {
    match serde_json::from_str::<ColliderConfig>(collider_json) {
        Ok(cfg) => engine.add_collider(body_handle, &cfg),
        Err(_) => u32::MAX,
    }
}

// ---------------------------------------------------------------------------
// WebGL rendering helpers
// ---------------------------------------------------------------------------

/// A flat buffer of per-body transform matrices (column-major 4×4, 16 floats each)
/// ready for upload to a WebGL uniform array or instanced transform buffer.
///
/// Each body contributes 16 f32 values; the layout is column-major as expected
/// by `gl.uniformMatrix4fv`.
#[allow(dead_code)]
pub fn engine_get_webgl_transform_matrices(engine: &WasmPhysicsEngine) -> Vec<f32> {
    let transforms = engine.get_all_transforms(); // [x,y,z, qx,qy,qz,qw] per body
    let n = transforms.len() / 7;
    let mut matrices = Vec::with_capacity(n * 16);
    for i in 0..n {
        let base = i * 7;
        let (px, py, pz) = (transforms[base], transforms[base + 1], transforms[base + 2]);
        let (qx, qy, qz, qw) = (
            transforms[base + 3] as f32,
            transforms[base + 4] as f32,
            transforms[base + 5] as f32,
            transforms[base + 6] as f32,
        );
        // Build 4×4 rotation-translation matrix from quaternion (column-major)
        let x2 = qx * qx;
        let y2 = qy * qy;
        let z2 = qz * qz;
        let xy = qx * qy;
        let xz = qx * qz;
        let yz = qy * qz;
        let wx = qw * qx;
        let wy = qw * qy;
        let wz = qw * qz;
        // Column 0
        matrices.push(1.0 - 2.0 * (y2 + z2));
        matrices.push(2.0 * (xy + wz));
        matrices.push(2.0 * (xz - wy));
        matrices.push(0.0);
        // Column 1
        matrices.push(2.0 * (xy - wz));
        matrices.push(1.0 - 2.0 * (x2 + z2));
        matrices.push(2.0 * (yz + wx));
        matrices.push(0.0);
        // Column 2
        matrices.push(2.0 * (xz + wy));
        matrices.push(2.0 * (yz - wx));
        matrices.push(1.0 - 2.0 * (x2 + y2));
        matrices.push(0.0);
        // Column 3 (translation)
        matrices.push(px as f32);
        matrices.push(py as f32);
        matrices.push(pz as f32);
        matrices.push(1.0);
    }
    matrices
}

/// Produce a flat `Vec`f32` of per-body colors `[r, g, b, a]` for WebGL instanced rendering.
///
/// Active bodies → `active_color`, sleeping bodies → `sleep_color`.
#[allow(dead_code)]
pub fn engine_get_body_colors(
    engine: &WasmPhysicsEngine,
    active_color: [f32; 4],
    sleep_color: [f32; 4],
) -> Vec<f32> {
    let handles = engine.get_all_body_handles();
    let mut colors = Vec::with_capacity(handles.len() * 4);
    for h in &handles {
        let is_sleeping = engine
            .get_body_state(*h)
            .map(|s| s.is_sleeping)
            .unwrap_or(false);
        let c = if is_sleeping {
            sleep_color
        } else {
            active_color
        };
        colors.extend_from_slice(&c);
    }
    colors
}

/// Build a JSON string describing the WebGL scene for use with a JS renderer.
///
/// The returned object has:
/// - `"transforms"`: flat f32 array (16 per body, column-major 4×4)
/// - `"colors"`: flat f32 array (4 per body, RGBA)
/// - `"body_count"`: number of bodies
#[allow(dead_code)]
pub fn engine_get_webgl_scene_json(engine: &WasmPhysicsEngine) -> String {
    let matrices = engine_get_webgl_transform_matrices(engine);
    let colors = engine_get_body_colors(engine, [0.2, 0.6, 1.0, 1.0], [0.5, 0.5, 0.5, 1.0]);
    serde_json::json!({
        "transforms": matrices,
        "colors": colors,
        "body_count": engine.get_body_count(),
    })
    .to_string()
}

// ---------------------------------------------------------------------------
// JS physics callback system
// ---------------------------------------------------------------------------

/// A record describing a physics event ready to be dispatched to JS.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct JsPhysicsEvent {
    /// Event category: "contact", "sleep", "wake", "step_complete".
    pub category: String,
    /// Event sub-type (e.g. "begin", "end", "body_slept").
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

/// Collect all contact events from the engine as `JsPhysicsEvent` objects.
#[allow(dead_code)]
pub fn engine_collect_contact_events(engine: &WasmPhysicsEngine) -> Vec<JsPhysicsEvent> {
    engine
        .get_contacts()
        .iter()
        .map(|c| JsPhysicsEvent::contact_begin(c.body_a, c.body_b, c.impulse, engine.time()))
        .collect()
}

/// Serialize all contact events from the engine to a JSON array string.
#[allow(dead_code)]
pub fn engine_collect_contact_events_json(engine: &WasmPhysicsEngine) -> String {
    let events = engine_collect_contact_events(engine);
    serde_json::to_string(&events).unwrap_or_else(|_| "[]".to_string())
}

// ---------------------------------------------------------------------------
// WASM memory management utilities
// ---------------------------------------------------------------------------

/// A simple arena that pre-allocates a contiguous `Vec`f64` scratch buffer
/// for zero-allocation JS↔WASM data transfers each frame.
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

/// Memory statistics for the WASM heap (estimated).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmMemoryStats {
    /// Estimated bytes used by all body states.
    pub body_state_bytes: usize,
    /// Estimated bytes used by contact data.
    pub contact_bytes: usize,
    /// Number of active bodies.
    pub body_count: usize,
    /// Number of active contacts.
    pub contact_count: usize,
}

/// Compute estimated memory usage of the engine.
#[allow(dead_code)]
pub fn engine_memory_stats(engine: &WasmPhysicsEngine) -> WasmMemoryStats {
    let body_count = engine.get_body_count() as usize;
    let contact_count = engine.get_contact_count() as usize;
    // Each body state: ~7 f64s for position + quaternion = 56 bytes,
    // plus velocity, angular velocity, mass etc. ≈ 200 bytes total.
    let body_state_bytes = body_count * 200;
    // Each contact: ~5 f64s (point, normal, depth, impulse) = 40 bytes.
    let contact_bytes = contact_count * 40;
    WasmMemoryStats {
        body_state_bytes,
        contact_bytes,
        body_count,
        contact_count,
    }
}

/// Serialize memory statistics to JSON.
#[allow(dead_code)]
pub fn engine_memory_stats_json(engine: &WasmPhysicsEngine) -> String {
    let stats = engine_memory_stats(engine);
    serde_json::to_string(&stats).unwrap_or_else(|_| "{}".to_string())
}

// ---------------------------------------------------------------------------
// JS event system
// ---------------------------------------------------------------------------

/// An ordered queue of `JsPhysicsEvent` objects accumulated over simulation frames.
///
/// Intended to be flushed once per render frame and dispatched to JS event listeners
/// via `postMessage`.
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

// ---------------------------------------------------------------------------
// Physics state serialization to JSON (extended)
// ---------------------------------------------------------------------------

/// Extended serialized state including velocity histograms and per-body energy.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct JsExtendedState {
    /// Base simulation state.
    pub base: JsSimulationState,
    /// Per-body kinetic energy (0.5 * m * v²), one entry per body (mass=1 proxy).
    pub kinetic_energies: Vec<f64>,
    /// Bounding box: `\[min_x, min_y, min_z, max_x, max_y, max_z\]`.
    pub bounding_box: [f64; 6],
    /// Center of mass of all dynamic bodies.
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

/// Build an extended state snapshot from the engine.
#[allow(dead_code)]
pub fn build_extended_state(engine: &WasmPhysicsEngine) -> JsExtendedState {
    let body_states: Vec<crate::types::BodyState> = engine
        .get_all_body_handles()
        .into_iter()
        .filter_map(|h| engine.get_body_state(h))
        .collect();

    let kinetic_energies: Vec<f64> = body_states
        .iter()
        .map(|s| {
            let v = &s.linear_velocity;
            0.5 * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2])
        })
        .collect();

    let bb = engine_bounding_box(engine);
    let com = engine_center_of_mass(engine);

    let base = JsSimulationState {
        time: engine.time(),
        num_bodies: engine.get_body_count(),
        positions: engine.get_all_positions(),
        gravity: engine.gravity(),
        contact_count: engine.get_contact_count(),
        body_states,
        debug_info: engine.debug_info().clone(),
    };

    JsExtendedState {
        base,
        kinetic_energies,
        bounding_box: bb,
        center_of_mass: [com.x, com.y, com.z],
        contact_count: engine.get_contact_count(),
    }
}

/// Serialize the extended state to JSON.
#[allow(dead_code)]
pub fn serialize_extended_state(engine: &WasmPhysicsEngine) -> String {
    build_extended_state(engine).to_json()
}

/// Compute total kinetic energy of all dynamic bodies (mass=1 proxy).
#[allow(dead_code)]
pub fn engine_total_kinetic_energy(engine: &WasmPhysicsEngine) -> f64 {
    engine
        .get_all_body_handles()
        .into_iter()
        .filter_map(|h| engine.get_body_state(h))
        .filter(|s| !s.is_sleeping)
        .map(|s| {
            let v = &s.linear_velocity;
            0.5 * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2])
        })
        .sum()
}

/// Serialize a compact per-body velocity magnitude array to JSON.
#[allow(dead_code)]
pub fn engine_get_speed_array_json(engine: &WasmPhysicsEngine) -> String {
    let speeds: Vec<f64> = engine
        .get_all_body_handles()
        .into_iter()
        .filter_map(|h| engine.get_body_state(h))
        .map(|s| {
            let v = &s.linear_velocity;
            (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
        })
        .collect();
    serde_json::to_string(&speeds).unwrap_or_else(|_| "[]".to_string())
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_config_builder() {
        let config = JsPhysicsConfig::new()
            .with_gravity(0.0, -20.0, 0.0)
            .with_fixed_dt(1.0 / 120.0)
            .with_max_substeps(8)
            .with_solver_iterations(12);

        assert!((config.gravity_y + 20.0).abs() < 1e-10);
        assert!((config.fixed_dt - 1.0 / 120.0).abs() < 1e-10);
        assert_eq!(config.max_substeps, 8);
        assert_eq!(config.solver_iterations, 12);
    }

    #[test]
    fn test_config_json_roundtrip() {
        let config = JsPhysicsConfig::new().with_gravity(0.0, -5.0, 0.0);
        let json = config.to_json();
        let restored = JsPhysicsConfig::from_json(&json).expect("should parse");
        assert!((restored.gravity_y + 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_serialize_state() {
        let mut engine = WasmPhysicsEngine::new(0.0, -9.81, 0.0);
        engine.add_dynamic_body(1.0, 1.0, 2.0, 3.0);

        let json = serialize_state(&engine);
        let state = deserialize_state(&json).expect("should deserialize");
        assert_eq!(state.num_bodies, 1);
        assert_eq!(state.positions.len(), 3);
        assert!((state.positions[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_create_engine() {
        let engine = create_engine(0.0, -9.81, 0.0);
        assert_eq!(engine.get_body_count(), 0);
    }

    #[test]
    fn test_create_engine_from_json() {
        let json = r#"{"gravity_x":0.0,"gravity_y":-1.62,"gravity_z":0.0,"fixed_dt":0.016,"max_substeps":4,"solver_iterations":8,"ccd_enabled":false,"sleeping_enabled":true,"linear_sleep_threshold":0.01,"angular_sleep_threshold":0.01}"#;
        let engine = create_engine_from_json(json);
        assert!((engine.gravity()[1] + 1.62).abs() < 1e-10);
    }

    #[test]
    fn test_engine_add_sphere() {
        let mut engine = create_engine(0.0, -9.81, 0.0);
        let h = engine_add_sphere(&mut engine, 1.0, 0.0, 10.0, 0.0, 0.5);
        assert_eq!(engine_get_body_count(&engine), 1);
        let pos = engine_get_position(&engine, h);
        assert!((pos[1] - 10.0).abs() < 1e-10);
    }

    #[test]
    fn test_engine_add_box() {
        let mut engine = create_engine(0.0, -9.81, 0.0);
        let h = engine_add_box(&mut engine, 2.0, 1.0, 2.0, 3.0, 0.5, 0.5, 0.5);
        let pos = engine_get_position(&engine, h);
        assert!((pos[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_engine_add_static_plane() {
        let mut engine = create_engine(0.0, -9.81, 0.0);
        let _h = engine_add_static_plane(&mut engine, 0.0, 1.0, 0.0, 0.0);
        assert_eq!(engine_get_body_count(&engine), 1);
    }

    #[test]
    fn test_engine_step_and_get_position() {
        let mut engine = create_engine(0.0, -10.0, 0.0);
        let h = engine_add_sphere(&mut engine, 1.0, 0.0, 100.0, 0.0, 0.5);
        engine_step(&mut engine, 1.0);
        let pos = engine_get_position(&engine, h);
        assert!(pos[1] < 100.0, "body should have fallen");
    }

    #[test]
    fn test_engine_set_gravity() {
        let mut engine = create_engine(0.0, -9.81, 0.0);
        engine_set_gravity(&mut engine, 0.0, -1.62, 0.0);
        assert!((engine.gravity()[1] + 1.62).abs() < 1e-10);
    }

    #[test]
    fn test_engine_apply_impulse() {
        let mut engine = create_engine(0.0, 0.0, 0.0);
        let h = engine_add_sphere(&mut engine, 2.0, 0.0, 0.0, 0.0, 0.5);
        let ok = engine_apply_impulse(&mut engine, h, 10.0, 0.0, 0.0);
        assert!(ok);
        let vel = engine_get_velocity(&engine, h);
        assert!((vel[0] - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_engine_apply_force() {
        let mut engine = create_engine(0.0, 0.0, 0.0);
        let h = engine_add_sphere(&mut engine, 1.0, 0.0, 0.0, 0.0, 0.5);
        engine.set_body_linear_damping(h, 0.0).unwrap();
        let ok = engine_apply_force(&mut engine, h, 10.0, 0.0, 0.0);
        assert!(ok);
        engine.set_fixed_dt(1.0);
        engine_step(&mut engine, 1.0);
        let vel = engine_get_velocity(&engine, h);
        assert!((vel[0] - 10.0).abs() < 1e-6);
    }

    #[test]
    fn test_engine_remove_body() {
        let mut engine = create_engine(0.0, -9.81, 0.0);
        let h = engine_add_sphere(&mut engine, 1.0, 0.0, 0.0, 0.0, 0.5);
        assert!(engine_remove_body(&mut engine, h));
        assert!(!engine_remove_body(&mut engine, h)); // already removed
    }

    #[test]
    fn test_engine_set_velocity() {
        let mut engine = create_engine(0.0, 0.0, 0.0);
        let h = engine_add_sphere(&mut engine, 1.0, 0.0, 0.0, 0.0, 0.5);
        let ok = engine_set_velocity(&mut engine, h, 5.0, 0.0, 0.0);
        assert!(ok);
        let vel = engine_get_velocity(&engine, h);
        assert!((vel[0] - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_engine_reset() {
        let mut engine = create_engine(0.0, -9.81, 0.0);
        engine_add_sphere(&mut engine, 1.0, 0.0, 0.0, 0.0, 0.5);
        engine_step(&mut engine, 0.1);
        engine_reset(&mut engine);
        assert_eq!(engine_get_body_count(&engine), 0);
    }

    #[test]
    fn test_engine_serialize() {
        let mut engine = create_engine(0.0, -9.81, 0.0);
        engine_add_sphere(&mut engine, 1.0, 1.0, 2.0, 3.0, 0.5);
        let json = engine_serialize(&engine);
        assert!(json.contains("num_bodies"));
    }

    #[test]
    fn test_state_hash_changes() {
        let mut engine = create_engine(0.0, -10.0, 0.0);
        engine_add_sphere(&mut engine, 1.0, 0.0, 100.0, 0.0, 0.5);
        let h0 = engine_state_hash(&engine);
        engine_step(&mut engine, 1.0);
        let h1 = engine_state_hash(&engine);
        assert_ne!(h0, h1, "hash should change after step");
    }

    #[test]
    fn test_setup_falling_spheres_scene() {
        let mut engine = create_engine(0.0, -9.81, 0.0);
        let handles = setup_falling_spheres_scene(&mut engine, 4, 0.5, 10.0);
        assert_eq!(handles.len(), 4);
        // 4 spheres + 1 ground plane
        assert_eq!(engine_get_body_count(&engine), 5);
    }

    #[test]
    fn test_bounding_box() {
        let mut engine = create_engine(0.0, 0.0, 0.0);
        engine_add_sphere(&mut engine, 1.0, -5.0, 0.0, 0.0, 0.5);
        engine_add_sphere(&mut engine, 1.0, 5.0, 10.0, 0.0, 0.5);
        let bb = engine_bounding_box(&engine);
        assert!((bb[0] + 5.0).abs() < 1e-10); // min_x
        assert!((bb[3] - 5.0).abs() < 1e-10); // max_x
        assert!((bb[4] - 10.0).abs() < 1e-10); // max_y
    }

    #[test]
    fn test_center_of_mass() {
        let mut engine = create_engine(0.0, 0.0, 0.0);
        engine_add_sphere(&mut engine, 1.0, -2.0, 0.0, 0.0, 0.5);
        engine_add_sphere(&mut engine, 1.0, 2.0, 0.0, 0.0, 0.5);
        let com = engine_center_of_mass(&engine);
        assert!(com.x.abs() < 1e-10, "COM.x should be 0, got {}", com.x);
    }

    #[test]
    fn test_batch_add_from_json() {
        let mut engine = create_engine(0.0, -9.81, 0.0);
        let json = r#"[
            {"mass":1.0,"position":[0.0,5.0,0.0],"velocity":[0.0,0.0,0.0],"shape":"sphere","radius":0.5,"half_extents":[0.5,0.5,0.5],"height":1.0,"restitution":0.3,"friction":0.5,"tag":"ball1"},
            {"mass":1.0,"position":[2.0,5.0,0.0],"velocity":[0.0,0.0,0.0],"shape":"box","radius":0.5,"half_extents":[0.5,0.5,0.5],"height":1.0,"restitution":0.3,"friction":0.5,"tag":"box1"}
        ]"#;
        let count = engine_add_bodies_from_json(&mut engine, json);
        assert_eq!(count, 2);
        assert_eq!(engine_get_body_count(&engine), 2);
    }

    #[test]
    fn test_contact_events_json() {
        let mut engine = create_engine(0.0, 0.0, 0.0);
        let b0 = engine_add_sphere(&mut engine, 1.0, 0.0, 0.0, 0.0, 1.0);
        let b1 = engine_add_sphere(&mut engine, 1.0, 1.5, 0.0, 0.0, 1.0);
        engine_set_velocity(&mut engine, b0, 2.0, 0.0, 0.0);
        engine_set_velocity(&mut engine, b1, -2.0, 0.0, 0.0);
        engine_step(&mut engine, 1.0 / 60.0);
        let json = engine_get_contact_events_json(&engine);
        // Should be a JSON array
        assert!(json.starts_with('['));
    }

    #[test]
    fn test_raycast_json() {
        let mut engine = create_engine(0.0, 0.0, 0.0);
        let b = engine_add_sphere(&mut engine, 1.0, 10.0, 0.0, 0.0, 1.0);
        let _ = b;
        let json = engine_raycast_json(&engine, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 100.0);
        assert!(json.contains("hit"));
    }

    #[test]
    fn test_get_all_transforms() {
        let mut engine = create_engine(0.0, 0.0, 0.0);
        engine_add_sphere(&mut engine, 1.0, 1.0, 2.0, 3.0, 0.5);
        let t = engine_get_all_transforms(&engine);
        assert_eq!(t.len(), 7); // 3 position + 4 quaternion
    }

    #[test]
    fn test_js_body_desc_constructors() {
        let s = JsBodyDesc::sphere(1.0, 0.0, 5.0, 0.0, 0.5);
        assert_eq!(s.shape, "sphere");
        assert!((s.radius - 0.5).abs() < 1e-10);

        let b = JsBodyDesc::cuboid(2.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        assert_eq!(b.shape, "box");

        let p = JsBodyDesc::static_plane(0.0, 1.0, 0.0, 0.0);
        assert_eq!(p.mass, 0.0);
    }

    #[test]
    fn test_js_contact_event_begin() {
        let cr = ContactResult {
            body_a: 0,
            body_b: 1,
            point_on_a: [1.0, 0.0, 0.0],
            point_on_b: [0.0, 0.0, 0.0],
            normal: [1.0, 0.0, 0.0],
            depth: 0.1,
            relative_velocity: -1.0,
            impulse: 2.0,
            is_new: true,
            friction_impulse: 0.5,
        };
        let ev = JsContactEvent::begin(&cr);
        assert_eq!(ev.event_type, "contact_begin");
        assert_eq!(ev.body_a, 0);
        let json = ev.to_json();
        assert!(json.contains("contact_begin"));
    }

    // ── WebGL rendering helpers ──────────────────────────────────────────

    #[test]
    fn test_webgl_transform_matrices_empty() {
        let engine = create_engine(0.0, 0.0, 0.0);
        let mats = engine_get_webgl_transform_matrices(&engine);
        assert!(mats.is_empty());
    }

    #[test]
    fn test_webgl_transform_matrices_one_body() {
        let mut engine = create_engine(0.0, 0.0, 0.0);
        engine_add_sphere(&mut engine, 1.0, 1.0, 2.0, 3.0, 0.5);
        let mats = engine_get_webgl_transform_matrices(&engine);
        // One body → 16 floats
        assert_eq!(mats.len(), 16);
        // Translation column (column 3): indices 12, 13, 14 = x, y, z
        assert!(
            (mats[12] - 1.0).abs() < 1e-5,
            "translation x should be 1.0, got {}",
            mats[12]
        );
        assert!(
            (mats[13] - 2.0).abs() < 1e-5,
            "translation y should be 2.0, got {}",
            mats[13]
        );
        assert!(
            (mats[14] - 3.0).abs() < 1e-5,
            "translation z should be 3.0, got {}",
            mats[14]
        );
        assert!((mats[15] - 1.0).abs() < 1e-5, "w should be 1.0");
    }

    #[test]
    fn test_webgl_transform_matrices_two_bodies() {
        let mut engine = create_engine(0.0, 0.0, 0.0);
        engine_add_sphere(&mut engine, 1.0, 0.0, 0.0, 0.0, 0.5);
        engine_add_sphere(&mut engine, 1.0, 5.0, 0.0, 0.0, 0.5);
        let mats = engine_get_webgl_transform_matrices(&engine);
        assert_eq!(mats.len(), 32);
    }

    #[test]
    fn test_engine_get_body_colors() {
        let mut engine = create_engine(0.0, 0.0, 0.0);
        engine_add_sphere(&mut engine, 1.0, 0.0, 0.0, 0.0, 0.5);
        engine_add_sphere(&mut engine, 1.0, 5.0, 0.0, 0.0, 0.5);
        let active_color = [1.0f32, 0.0, 0.0, 1.0];
        let sleep_color = [0.0f32, 0.0, 1.0, 1.0];
        let colors = engine_get_body_colors(&engine, active_color, sleep_color);
        // 2 bodies × 4 channels
        assert_eq!(colors.len(), 8);
    }

    #[test]
    fn test_engine_get_webgl_scene_json() {
        let mut engine = create_engine(0.0, 0.0, 0.0);
        engine_add_sphere(&mut engine, 1.0, 0.0, 0.0, 0.0, 0.5);
        let json = engine_get_webgl_scene_json(&engine);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(v["body_count"], 1);
        assert!(v["transforms"].is_array());
        assert!(v["colors"].is_array());
    }

    // ── JS physics callback system ───────────────────────────────────────

    #[test]
    fn test_js_physics_event_step_complete() {
        let ev = JsPhysicsEvent::step_complete(1.5);
        assert_eq!(ev.category, "step");
        assert_eq!(ev.sub_type, "step_complete");
        assert!((ev.time - 1.5).abs() < 1e-10);
        let json = ev.to_json();
        assert!(json.contains("step_complete"));
    }

    #[test]
    fn test_js_physics_event_contact() {
        let ev = JsPhysicsEvent::contact_begin(0, 1, 5.0, 0.016);
        assert_eq!(ev.category, "contact");
        assert_eq!(ev.body_a, 0);
        assert_eq!(ev.body_b, 1);
        assert!((ev.value - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_js_physics_event_body_slept() {
        let ev = JsPhysicsEvent::body_slept(42, 2.0);
        assert_eq!(ev.category, "sleep");
        assert_eq!(ev.body_a, 42);
        assert_eq!(ev.body_b, u32::MAX);
    }

    #[test]
    fn test_js_physics_event_body_woke() {
        let ev = JsPhysicsEvent::body_woke(7, 3.0);
        assert_eq!(ev.category, "wake");
        assert_eq!(ev.sub_type, "body_woke");
    }

    #[test]
    fn test_engine_collect_contact_events_json() {
        let engine = create_engine(0.0, 0.0, 0.0);
        let json = engine_collect_contact_events_json(&engine);
        assert!(json.starts_with('['));
    }

    // ── WASM memory management utilities ────────────────────────────────

    #[test]
    fn test_wasm_scratch_buffer_fill_positions() {
        let mut engine = create_engine(0.0, 0.0, 0.0);
        engine_add_sphere(&mut engine, 1.0, 1.0, 2.0, 3.0, 0.5);
        let mut scratch = WasmScratchBuffer::new(64);
        let slice = scratch.fill_positions(&engine);
        assert_eq!(slice.len(), 3);
        assert!((slice[0] - 1.0).abs() < 1e-10);
        assert!((slice[1] - 2.0).abs() < 1e-10);
        assert!((slice[2] - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_wasm_scratch_buffer_fill_transforms() {
        let mut engine = create_engine(0.0, 0.0, 0.0);
        engine_add_sphere(&mut engine, 1.0, 0.0, 0.0, 0.0, 0.5);
        let mut scratch = WasmScratchBuffer::new(64);
        let slice = scratch.fill_transforms(&engine);
        // 1 body × 7 values
        assert_eq!(slice.len(), 7);
    }

    #[test]
    fn test_wasm_scratch_buffer_clear() {
        let mut engine = create_engine(0.0, 0.0, 0.0);
        engine_add_sphere(&mut engine, 1.0, 0.0, 0.0, 0.0, 0.5);
        let mut scratch = WasmScratchBuffer::new(16);
        scratch.fill_positions(&engine);
        assert!(!scratch.is_empty());
        scratch.clear();
        assert!(scratch.is_empty());
    }

    #[test]
    fn test_wasm_scratch_buffer_capacity() {
        let scratch = WasmScratchBuffer::new(128);
        assert_eq!(scratch.capacity(), 128);
    }

    #[test]
    fn test_engine_memory_stats() {
        let mut engine = create_engine(0.0, 0.0, 0.0);
        engine_add_sphere(&mut engine, 1.0, 0.0, 0.0, 0.0, 0.5);
        engine_add_sphere(&mut engine, 1.0, 5.0, 0.0, 0.0, 0.5);
        let stats = engine_memory_stats(&engine);
        assert_eq!(stats.body_count, 2);
        assert!(stats.body_state_bytes > 0);
    }

    #[test]
    fn test_engine_memory_stats_json() {
        let mut engine = create_engine(0.0, 0.0, 0.0);
        engine_add_sphere(&mut engine, 1.0, 0.0, 0.0, 0.0, 0.5);
        let json = engine_memory_stats_json(&engine);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(v["body_count"], 1);
    }

    // ── JS event system ──────────────────────────────────────────────────

    #[test]
    fn test_js_event_queue_push_drain() {
        let mut queue = JsEventQueue::new(0);
        queue.push(JsPhysicsEvent::step_complete(0.5));
        queue.push(JsPhysicsEvent::step_complete(1.0));
        assert_eq!(queue.len(), 2);
        let drained = queue.drain();
        assert_eq!(drained.len(), 2);
        assert!(queue.is_empty());
    }

    #[test]
    fn test_js_event_queue_max_events() {
        let mut queue = JsEventQueue::new(3);
        for i in 0..5 {
            queue.push(JsPhysicsEvent::step_complete(i as f64));
        }
        assert_eq!(queue.len(), 3, "queue should be capped at 3");
    }

    #[test]
    fn test_js_event_queue_to_json() {
        let mut queue = JsEventQueue::new(0);
        queue.push(JsPhysicsEvent::step_complete(1.0));
        let json = queue.to_json();
        assert!(json.starts_with('['));
        assert!(json.contains("step_complete"));
    }

    #[test]
    fn test_js_event_queue_drain_to_json() {
        let mut queue = JsEventQueue::new(0);
        queue.push(JsPhysicsEvent::body_slept(5, 1.0));
        let json = queue.drain_to_json();
        assert!(json.contains("body_slept"));
        assert!(queue.is_empty());
    }

    #[test]
    fn test_js_event_queue_push_step_complete() {
        let mut queue = JsEventQueue::new(0);
        queue.push_step_complete(2.5);
        assert_eq!(queue.len(), 1);
        let evs = queue.drain();
        assert!((evs[0].time - 2.5).abs() < 1e-10);
    }

    #[test]
    fn test_js_event_queue_collect_contacts() {
        let mut engine = create_engine(0.0, 0.0, 0.0);
        engine_add_sphere(&mut engine, 1.0, 0.0, 0.0, 0.0, 1.0);
        engine_add_sphere(&mut engine, 1.0, 1.5, 0.0, 0.0, 1.0);
        engine_step(&mut engine, 1.0 / 60.0);
        let mut queue = JsEventQueue::new(0);
        queue.collect_contacts_from_engine(&engine);
        // Just ensure it doesn't panic; number of contacts depends on engine state
        let _json = queue.to_json();
    }

    // ── Extended state serialization ─────────────────────────────────────

    #[test]
    fn test_build_extended_state() {
        let mut engine = create_engine(0.0, 0.0, 0.0);
        engine_add_sphere(&mut engine, 1.0, 1.0, 0.0, 0.0, 0.5);
        engine_add_sphere(&mut engine, 1.0, -1.0, 0.0, 0.0, 0.5);
        let state = build_extended_state(&engine);
        assert_eq!(state.base.num_bodies, 2);
        assert_eq!(state.kinetic_energies.len(), 2);
    }

    #[test]
    fn test_serialize_extended_state_json() {
        let mut engine = create_engine(0.0, 0.0, 0.0);
        engine_add_sphere(&mut engine, 1.0, 0.0, 5.0, 0.0, 0.5);
        let json = serialize_extended_state(&engine);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert!(v["base"].is_object());
        assert!(v["kinetic_energies"].is_array());
        assert!(v["bounding_box"].is_array());
    }

    #[test]
    fn test_extended_state_roundtrip() {
        let mut engine = create_engine(0.0, 0.0, 0.0);
        engine_add_sphere(&mut engine, 1.0, 3.0, 0.0, 0.0, 0.5);
        let state = build_extended_state(&engine);
        let json = state.to_json();
        let back = JsExtendedState::from_json(&json).expect("should deserialize");
        assert_eq!(back.base.num_bodies, 1);
        assert!((back.bounding_box[0] - 3.0).abs() < 1e-9);
    }

    #[test]
    fn test_engine_total_kinetic_energy_at_rest() {
        let mut engine = create_engine(0.0, 0.0, 0.0);
        engine_add_sphere(&mut engine, 1.0, 0.0, 0.0, 0.0, 0.5);
        let ke = engine_total_kinetic_energy(&engine);
        assert!(
            ke.abs() < 1e-10,
            "resting body should have ~zero KE, got {ke}"
        );
    }

    #[test]
    fn test_engine_total_kinetic_energy_moving() {
        let mut engine = create_engine(0.0, 0.0, 0.0);
        let h = engine_add_sphere(&mut engine, 1.0, 0.0, 0.0, 0.0, 0.5);
        engine_set_velocity(&mut engine, h, 3.0, 4.0, 0.0);
        let ke = engine_total_kinetic_energy(&engine);
        // KE = 0.5 * (9 + 16) = 12.5
        assert!((ke - 12.5).abs() < 1e-6, "KE should be 12.5, got {ke}");
    }

    #[test]
    fn test_engine_get_speed_array_json() {
        let mut engine = create_engine(0.0, 0.0, 0.0);
        let h = engine_add_sphere(&mut engine, 1.0, 0.0, 0.0, 0.0, 0.5);
        engine_set_velocity(&mut engine, h, 3.0, 4.0, 0.0);
        let json = engine_get_speed_array_json(&engine);
        let speeds: Vec<f64> = serde_json::from_str(&json).unwrap();
        assert_eq!(speeds.len(), 1);
        assert!(
            (speeds[0] - 5.0).abs() < 1e-6,
            "speed should be 5.0, got {}",
            speeds[0]
        );
    }
}
