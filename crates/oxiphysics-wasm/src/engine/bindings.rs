// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! `WasmBindings` — flat wasm-bindgen-style API wrapper + `WasmContactList`.

use wasm_bindgen::prelude::*;

use super::WasmPhysicsEngine;
use super::wasm_types::WasmTransform;

// ===========================================================================
// WasmBindings — flat wasm-bindgen-style API wrapper
// ===========================================================================

/// A flat, JS-idiomatic wrapper around [`WasmPhysicsEngine`].
///
/// All method signatures use primitive types and `Vec<f64>` / `Vec<u32>`
/// instead of Rust structs, making them trivially mappable to wasm-bindgen
/// exported functions.
///
/// # Example
///
/// ```no_run
/// use oxiphysics_wasm::engine::WasmBindings;
///
/// let mut api = WasmBindings::new(0.0, -9.81, 0.0);
///
/// // Add a 1 kg body at height 10
/// let id = api.add_rigid_body(0.0, 10.0, 0.0, 1.0);
/// api.add_collider_sphere(id, 0.5);
///
/// api.step(1.0 / 60.0);
///
/// let pos = api.get_position(id);
/// assert!(pos[1] < 10.0, "body should fall");
/// assert_eq!(api.body_count(), 1);
/// ```
#[wasm_bindgen]
#[derive(Debug)]
pub struct WasmBindings {
    engine: WasmPhysicsEngine,
}

#[wasm_bindgen]
impl WasmBindings {
    // -----------------------------------------------------------------------
    // Construction
    // -----------------------------------------------------------------------

    /// Create a new engine with the given gravity components.
    pub fn new(gx: f64, gy: f64, gz: f64) -> Self {
        Self {
            engine: WasmPhysicsEngine::new(gx, gy, gz),
        }
    }

    // -----------------------------------------------------------------------
    // Body management
    // -----------------------------------------------------------------------

    /// Add a rigid body at `(x, y, z)` with the given `mass` (kg).
    ///
    /// Returns the stable body handle (`u32`).  Pass `mass = 0.0` for a
    /// static (immovable) body.
    pub fn add_rigid_body(&mut self, x: f64, y: f64, z: f64, mass: f64) -> u32 {
        if mass <= 0.0 {
            self.engine.add_static_body(x, y, z)
        } else {
            self.engine.add_dynamic_body(mass, x, y, z)
        }
    }

    /// Set linear velocity of body `id`.
    ///
    /// Silently ignores invalid handles (JS callers should not need to
    /// propagate physics errors in hot loops).
    pub fn set_velocity(&mut self, id: u32, vx: f64, vy: f64, vz: f64) {
        let _ = self.engine.set_velocity(id, vx, vy, vz);
    }

    /// Get position of body `id` as a flat `[x, y, z]` `Vec<f64>`.
    ///
    /// In a wasm-bindgen build this would return a `Float64Array` view.
    pub fn get_position(&self, id: u32) -> Vec<f64> {
        self.engine.get_position(id).to_vec()
    }

    /// Advance the simulation by `dt` seconds.
    pub fn step(&mut self, dt: f64) {
        self.engine.step(dt);
    }

    /// Return all active body positions as a flat `Vec<f64>`:
    /// `[x0, y0, z0, x1, y1, z1, ...]`.
    ///
    /// Suitable for wrapping in a `Float64Array` on the JS side.
    pub fn get_all_positions(&self) -> Vec<f64> {
        self.engine.get_all_positions()
    }

    /// Set the gravity vector.
    pub fn set_gravity(&mut self, gx: f64, gy: f64, gz: f64) {
        self.engine.set_gravity(gx, gy, gz);
    }

    /// Return the number of active bodies.
    pub fn body_count(&self) -> u32 {
        self.engine.get_body_count()
    }

    // -----------------------------------------------------------------------
    // Collider management
    // -----------------------------------------------------------------------

    /// Add a sphere collider to `body_id` with the given `radius`.
    ///
    /// Returns the collider handle.
    pub fn add_collider_sphere(&mut self, body_id: u32, radius: f64) -> u32 {
        self.engine.add_sphere_collider(body_id, radius)
    }

    /// Add a box collider to `body_id` with half-extents `(hx, hy, hz)`.
    ///
    /// Returns the collider handle.
    pub fn add_collider_box(&mut self, body_id: u32, hx: f64, hy: f64, hz: f64) -> u32 {
        self.engine.add_box_collider(body_id, hx, hy, hz)
    }

    // -----------------------------------------------------------------------
    // Contacts
    // -----------------------------------------------------------------------

    /// Return all contacts from the last `step()` call as a flat `Vec<f64>`.
    ///
    /// Each contact is encoded as 12 consecutive floats:
    /// `[body_a, body_b, nx, ny, nz, depth, rel_vel, impulse, is_new,
    ///   friction_impulse, point_ax, point_ay]`
    ///
    /// (body handles are transmitted as `f64` for uniformity.)
    pub fn get_contacts_flat(&self) -> Vec<f64> {
        let contacts = self.engine.get_contacts();
        let mut out = Vec::with_capacity(contacts.len() * 12);
        for c in contacts {
            out.push(c.body_a as f64);
            out.push(c.body_b as f64);
            out.push(c.normal[0]);
            out.push(c.normal[1]);
            out.push(c.normal[2]);
            out.push(c.depth);
            out.push(c.relative_velocity);
            out.push(c.impulse);
            out.push(if c.is_new { 1.0 } else { 0.0 });
            out.push(c.friction_impulse);
            out.push(c.point_on_a[0]);
            out.push(c.point_on_a[1]);
        }
        out
    }

    /// Return contact info as a `Vec` of JSON strings, one per contact.
    ///
    /// On the JS side each string can be parsed with `JSON.parse`.
    pub fn get_contacts_json(&self) -> Vec<String> {
        self.engine
            .get_contacts()
            .iter()
            .map(|c| {
                format!(
                    r#"{{"bodyA":{},"bodyB":{},"normal":[{},{},{}],"depth":{},"impulse":{}}}"#,
                    c.body_a, c.body_b, c.normal[0], c.normal[1], c.normal[2], c.depth, c.impulse,
                )
            })
            .collect()
    }

    /// Return the number of contacts detected in the last step.
    pub fn contact_count(&self) -> u32 {
        self.engine.get_contact_count()
    }

    // -----------------------------------------------------------------------
    // WebGPU compute placeholder
    // -----------------------------------------------------------------------

    /// WebGPU compute shader placeholder.
    ///
    /// The mock result contains:
    /// - slot 0: body count (as f64)
    /// - slot 1: simulation time
    /// - remaining: zeros
    ///
    /// # Example
    ///
    /// ```no_run
    /// use oxiphysics_wasm::engine::WasmBindings;
    ///
    /// let mut api = WasmBindings::new(0.0, -9.81, 0.0);
    /// api.add_rigid_body(0.0, 1.0, 0.0, 1.0);
    /// api.step(1.0 / 60.0);
    ///
    /// let result = api.webgpu_compute_placeholder(64);
    /// assert_eq!(result.len(), 64);
    /// assert!((result[0] - 1.0).abs() < 1e-10); // body count
    /// ```
    pub fn webgpu_compute_placeholder(&self, workgroup_size: u32) -> Vec<f64> {
        let n = (workgroup_size as usize).max(2);
        let mut out = vec![0.0f64; n];
        out[0] = self.engine.get_body_count() as f64;
        out[1] = self.engine.time();
        out
    }

    // -----------------------------------------------------------------------
    // Transform helpers
    // -----------------------------------------------------------------------

    /// Get the full transform of body `id` as a `WasmTransform`.
    ///
    /// Returns the identity transform if the handle is invalid.
    pub fn get_transform(&self, id: u32) -> WasmTransform {
        let pos = self.engine.get_position(id);
        let rot = self.engine.get_rotation(id);
        WasmTransform::new(pos[0], pos[1], pos[2], rot[0], rot[1], rot[2], rot[3])
    }

    /// Get all transforms as a flat `Vec<f64>` of 7-element blocks
    /// `[px, py, pz, qx, qy, qz, qw, ...]`.
    pub fn get_all_transforms_flat(&self) -> Vec<f64> {
        self.engine.get_all_transforms()
    }

    // -----------------------------------------------------------------------
    // Gravity query
    // -----------------------------------------------------------------------

    /// Return the current gravity as `Vec<f64>` of `[gx, gy, gz]`.
    pub fn gravity(&self) -> Vec<f64> {
        self.engine.gravity().to_vec()
    }

    // -----------------------------------------------------------------------
    // Serialization helpers
    // -----------------------------------------------------------------------

    /// Serialize the current simulation state to a compact JSON string.
    pub fn state_to_json(&self) -> String {
        let bodies = self.engine.get_all_body_handles();
        let mut entries = Vec::with_capacity(bodies.len());
        for h in &bodies {
            let p = self.engine.get_position(*h);
            let v = self.engine.get_velocity(*h);
            entries.push(format!(
                r#"{{"id":{},"pos":[{},{},{}],"vel":[{},{},{}]}}"#,
                h, p[0], p[1], p[2], v[0], v[1], v[2]
            ));
        }
        format!(
            r#"{{"time":{},"gravity":[{},{},{}],"bodies":[{}]}}"#,
            self.engine.time(),
            self.engine.gravity()[0],
            self.engine.gravity()[1],
            self.engine.gravity()[2],
            entries.join(",")
        )
    }

    // -----------------------------------------------------------------------
    // Reset helpers
    // -----------------------------------------------------------------------

    /// Reset the simulation (remove all bodies / colliders / contacts).
    pub fn reset(&mut self) {
        self.engine.reset();
    }

    /// Apply a force to body `id` (silently ignores invalid handles).
    pub fn apply_force(&mut self, id: u32, fx: f64, fy: f64, fz: f64) {
        let _ = self.engine.apply_force(id, fx, fy, fz);
    }

    /// Apply an impulse to body `id`.
    pub fn apply_impulse(&mut self, id: u32, ix: f64, iy: f64, iz: f64) {
        let _ = self.engine.apply_impulse(id, ix, iy, iz);
    }

    /// Return contacts as flat data: 8 floats per contact
    /// `[body_a, body_b, nx, ny, nz, depth, impulse, is_new]`.
    pub fn get_contacts_structured_flat(&self) -> Vec<f64> {
        WasmContactList::from_contacts(self.engine.get_contacts()).to_flat()
    }

    /// Return contacts as a JSON array string.
    pub fn get_contacts_structured_json(&self) -> String {
        WasmContactList::from_contacts(self.engine.get_contacts()).to_json()
    }
}

impl WasmBindings {
    /// Return contacts as a `WasmContactList` (Rust-only).
    pub fn get_contacts_structured(&self) -> WasmContactList {
        WasmContactList::from_contacts(self.engine.get_contacts())
    }
}

// ===========================================================================
// WasmContactList — structured contact result list
// ===========================================================================

/// A structured list of contact results from a simulation step.
///
/// Provides both array-indexed access and JSON serialization.
#[wasm_bindgen]
#[derive(Debug)]
pub struct WasmContactList {
    items: Vec<ContactInfoEntry>,
}

/// A single contact event returned from `WasmBindings::get_contacts_structured`.
#[derive(Debug, Clone)]
pub struct ContactInfoEntry {
    /// Handle of body A.
    pub body_a: u32,
    /// Handle of body B.
    pub body_b: u32,
    /// Contact normal `[nx, ny, nz]`.
    pub normal: [f64; 3],
    /// Penetration depth (m).
    pub depth: f64,
    /// Normal impulse magnitude.
    pub impulse: f64,
    /// Whether this is a new contact (not present in the previous step).
    pub is_new: bool,
}

impl ContactInfoEntry {
    /// Serialize to JSON.
    pub fn to_json(&self) -> String {
        format!(
            r#"{{"bodyA":{},"bodyB":{},"normal":[{},{},{}],"depth":{},"impulse":{},"isNew":{}}}"#,
            self.body_a,
            self.body_b,
            self.normal[0],
            self.normal[1],
            self.normal[2],
            self.depth,
            self.impulse,
            self.is_new,
        )
    }
}

#[wasm_bindgen]
impl WasmContactList {
    /// Number of contacts.
    pub fn len(&self) -> u32 {
        self.items.len() as u32
    }

    /// Whether the list is empty.
    pub fn is_empty(&self) -> bool {
        self.items.is_empty()
    }

    /// Get body A handle for contact at index `i` (returns `u32::MAX` if out of bounds).
    pub fn get_body_a(&self, i: u32) -> u32 {
        self.items.get(i as usize).map_or(u32::MAX, |c| c.body_a)
    }

    /// Get body B handle for contact at index `i` (returns `u32::MAX` if out of bounds).
    pub fn get_body_b(&self, i: u32) -> u32 {
        self.items.get(i as usize).map_or(u32::MAX, |c| c.body_b)
    }

    /// Get contact normal as `Vec<f64>` of `[nx, ny, nz]` for contact `i`.
    pub fn get_normal(&self, i: u32) -> Vec<f64> {
        self.items
            .get(i as usize)
            .map_or_else(|| vec![0.0; 3], |c| c.normal.to_vec())
    }

    /// Get penetration depth for contact `i` (returns `0.0` if out of bounds).
    pub fn get_depth(&self, i: u32) -> f64 {
        self.items.get(i as usize).map_or(0.0, |c| c.depth)
    }

    /// Get impulse magnitude for contact `i`.
    pub fn get_impulse(&self, i: u32) -> f64 {
        self.items.get(i as usize).map_or(0.0, |c| c.impulse)
    }

    /// Whether contact `i` is new.
    pub fn get_is_new(&self, i: u32) -> bool {
        self.items.get(i as usize).is_some_and(|c| c.is_new)
    }

    /// Serialize the whole list to a JSON array string.
    pub fn to_json(&self) -> String {
        let entries: Vec<String> = self.items.iter().map(|c| c.to_json()).collect();
        format!("[{}]", entries.join(","))
    }

    /// Flatten to a `Vec<f64>` with 8 floats per contact:
    /// `[body_a, body_b, nx, ny, nz, depth, impulse, is_new]`.
    pub fn to_flat(&self) -> Vec<f64> {
        let mut out = Vec::with_capacity(self.items.len() * 8);
        for c in &self.items {
            out.push(c.body_a as f64);
            out.push(c.body_b as f64);
            out.push(c.normal[0]);
            out.push(c.normal[1]);
            out.push(c.normal[2]);
            out.push(c.depth);
            out.push(c.impulse);
            out.push(if c.is_new { 1.0 } else { 0.0 });
        }
        out
    }
}

impl WasmContactList {
    /// Build from the raw `ContactResult` slice from the engine (Rust-only).
    pub fn from_contacts(contacts: &[crate::types::ContactResult]) -> Self {
        let items = contacts
            .iter()
            .map(|c| ContactInfoEntry {
                body_a: c.body_a,
                body_b: c.body_b,
                normal: c.normal,
                depth: c.depth,
                impulse: c.impulse,
                is_new: c.is_new,
            })
            .collect();
        Self { items }
    }

    /// Get contact at index `i` (Rust-only; use per-field getters from JS).
    pub fn get(&self, i: usize) -> Option<&ContactInfoEntry> {
        self.items.get(i)
    }
}
