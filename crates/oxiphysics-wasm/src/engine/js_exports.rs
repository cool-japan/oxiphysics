// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! JavaScript-facing exports of [`WasmPhysicsEngine`].
//!
//! The methods in this module are thin `#[wasm_bindgen]` wrappers around the
//! Rust API of the engine.  Every operation that can fail (creating bodies or
//! colliders, stepping with an invalid time step, ...) returns a
//! `Result<_, JsValue>`, which wasm-bindgen turns into a thrown JavaScript
//! exception carrying the error message.

use wasm_bindgen::prelude::*;

use super::WasmPhysicsEngine;
use super::handles::decode_handle;

/// Result type of the exported methods: errors surface as JS exceptions.
type JsResult<T> = std::result::Result<T, JsValue>;

// ============================================================================
// wasm-bindgen JavaScript API
// ============================================================================
//
// This impl block exports `WasmPhysicsEngine` as a JavaScript class.
// All methods use primitive types or `Vec<f64>` (returned as `Float64Array`)
// so no JS/Rust type-conversion boilerplate is needed on the JS side.
//
// Naming convention: Rust method names end in `_js`; `js_name` attributes
// expose them under the clean, idiomatic JavaScript names.
//
// ## JavaScript usage
//
// ```js
// const engine = new WasmPhysicsEngine(0.0, -9.81, 0.0);
// const ball  = engine.add_dynamic_body(1.0, 0.0, 10.0, 0.0);
// engine.add_sphere_collider(ball, 0.5);
//
// const ground = engine.add_static_body(0.0, 0.0, 0.0);
// engine.add_plane_collider(ground, 0.0, 1.0, 0.0, 0.0);
//
// // step() consumes real time in fixed substeps (1/60 s by default, at most
// // 4 per call): feed it the frame time, not a large "simulate N seconds".
// for (let i = 0; i < 60; i++) engine.step(1 / 60);
//
// const pos = engine.get_position(ball);  // Float64Array [x, y, z]
// console.log(`y = ${pos[1].toFixed(3)}`);
//
// // Bulk positions (good for instanced mesh rendering)
// const all = engine.get_all_positions();  // Float64Array [x0,y0,z0, x1,y1,z1, ...]
// ```
//
// ## Errors
//
// Methods that create bodies or colliders (`add_*`) and `step` throw a
// JavaScript exception when given invalid input (non-positive mass or
// dimensions, NaN, a stale body handle, a negative time step, ...).
// Handles are opaque integers: do not do arithmetic on them or use them as
// array indices; use `handle_index` for that.

#[wasm_bindgen]
impl WasmPhysicsEngine {
    // -----------------------------------------------------------------------
    // Construction
    // -----------------------------------------------------------------------

    /// Create an engine with the given gravity components.
    ///
    /// Equivalent to `new WasmPhysicsEngine(0.0, -9.81, 0.0)` in JavaScript.
    #[wasm_bindgen(constructor)]
    pub fn wasm_new(gx: f64, gy: f64, gz: f64) -> WasmPhysicsEngine {
        WasmPhysicsEngine::new(gx, gy, gz)
    }

    // -----------------------------------------------------------------------
    // Body management
    // -----------------------------------------------------------------------

    /// Add a dynamic (physics-simulated) body at `(x, y, z)`.
    ///
    /// Returns the body handle (an opaque `u32` used in subsequent calls).
    /// Throws if `mass` is not finite and positive or the position is not
    /// finite.
    #[wasm_bindgen(js_name = "add_dynamic_body")]
    pub fn add_dynamic_body_js(&mut self, mass: f64, x: f64, y: f64, z: f64) -> JsResult<u32> {
        Ok(self.add_dynamic_body(mass, x, y, z)?)
    }

    /// Add a static (immovable) body at `(x, y, z)`.
    ///
    /// Throws if the position is not finite.
    #[wasm_bindgen(js_name = "add_static_body")]
    pub fn add_static_body_js(&mut self, x: f64, y: f64, z: f64) -> JsResult<u32> {
        Ok(self.add_static_body(x, y, z)?)
    }

    /// Remove a body and its colliders. Returns `true` on success and `false`
    /// when the handle is not (or no longer) valid.
    #[wasm_bindgen(js_name = "remove_body")]
    pub fn remove_body_js(&mut self, handle: u32) -> bool {
        self.remove_body(handle).is_ok()
    }

    /// `true` if `handle` refers to a live body (stale handles of removed
    /// bodies return `false` even after their slot has been reused).
    #[wasm_bindgen(js_name = "contains_body")]
    pub fn contains_body_js(&self, handle: u32) -> bool {
        self.contains_body(handle)
    }

    /// Slot index encoded in a handle (a small integer, usable as an array
    /// index in renderer-side tables; it is reused after a body is removed).
    #[wasm_bindgen(js_name = "handle_index")]
    pub fn handle_index_js(handle: u32) -> u32 {
        decode_handle(handle).0 as u32
    }

    /// Generation encoded in a handle (how often its slot was recycled,
    /// modulo 4096).
    #[wasm_bindgen(js_name = "handle_generation")]
    pub fn handle_generation_js(handle: u32) -> u32 {
        decode_handle(handle).1
    }

    /// Number of active (non-removed) bodies.
    #[wasm_bindgen(js_name = "get_body_count")]
    pub fn get_body_count_js(&self) -> u32 {
        self.get_body_count()
    }

    // -----------------------------------------------------------------------
    // Colliders
    // -----------------------------------------------------------------------

    /// Attach a sphere collider of `radius` to `body`. Returns collider handle.
    ///
    /// Throws if `body` is not a live body or `radius` is not finite and
    /// positive.
    #[wasm_bindgen(js_name = "add_sphere_collider")]
    pub fn add_sphere_collider_js(&mut self, body: u32, radius: f64) -> JsResult<u32> {
        Ok(self.add_sphere_collider(body, radius)?)
    }

    /// Attach an oriented box collider (half-extents along the body's local
    /// axes, so the box rotates with its body) to `body`. Returns collider
    /// handle.
    ///
    /// Throws if `body` is not a live body or a half-extent is not finite and
    /// positive.
    #[wasm_bindgen(js_name = "add_box_collider")]
    pub fn add_box_collider_js(&mut self, body: u32, hx: f64, hy: f64, hz: f64) -> JsResult<u32> {
        Ok(self.add_box_collider(body, hx, hy, hz)?)
    }

    /// Attach a capsule collider (axis along the body's local Y axis; `height`
    /// is the length of the cylindrical part) to `body`. Returns collider
    /// handle.
    ///
    /// Throws if `body` is not a live body or `radius`/`height` are not finite
    /// and positive.
    #[wasm_bindgen(js_name = "add_capsule_collider")]
    pub fn add_capsule_collider_js(
        &mut self,
        body: u32,
        radius: f64,
        height: f64,
    ) -> JsResult<u32> {
        Ok(self.add_capsule_collider(body, radius, height)?)
    }

    /// Attach a cylinder collider (axis along the body's local Y axis) to
    /// `body`. Returns collider handle.
    ///
    /// Throws if `body` is not a live body or `radius`/`height` are not finite
    /// and positive.
    #[wasm_bindgen(js_name = "add_cylinder_collider")]
    pub fn add_cylinder_collider_js(
        &mut self,
        body: u32,
        radius: f64,
        height: f64,
    ) -> JsResult<u32> {
        Ok(self.add_cylinder_collider(body, radius, height)?)
    }

    /// Attach an infinite plane collider (`n . x = offset`, normal need not be
    /// unit length). Returns collider handle.
    ///
    /// Throws if `body` is not a live body or the normal is zero or not finite.
    #[wasm_bindgen(js_name = "add_plane_collider")]
    pub fn add_plane_collider_js(
        &mut self,
        body: u32,
        nx: f64,
        ny: f64,
        nz: f64,
        offset: f64,
    ) -> JsResult<u32> {
        Ok(self.add_plane_collider(body, nx, ny, nz, offset)?)
    }

    // -----------------------------------------------------------------------
    // Simulation
    // -----------------------------------------------------------------------

    /// Advance the simulation by `dt` seconds of real time.
    ///
    /// The time is consumed in fixed substeps (see `set_fixed_dt`, 1/60 s by
    /// default) at most `max_substeps` per call; anything beyond that is
    /// dropped. Throws if `dt` is negative, `NaN` or infinite.
    #[wasm_bindgen(js_name = "step")]
    pub fn step_js(&mut self, dt: f64) -> JsResult<()> {
        Ok(self.step(dt)?)
    }

    /// Set the fixed substep length in seconds. Throws unless it is finite
    /// and positive.
    #[wasm_bindgen(js_name = "set_fixed_dt")]
    pub fn set_fixed_dt_js(&mut self, dt: f64) -> JsResult<()> {
        Ok(self.set_fixed_dt(dt)?)
    }

    /// Fraction in `[0, 1)` of a fixed substep that is still pending after the
    /// last `step()`; use it to interpolate rendering between physics states.
    #[wasm_bindgen(js_name = "interpolation_alpha")]
    pub fn interpolation_alpha_js(&self) -> f64 {
        self.interpolation_alpha()
    }

    /// Number of fixed substeps the last `step()` call executed (0 when the
    /// elapsed time did not fill a whole substep).
    #[wasm_bindgen(js_name = "substeps_performed")]
    pub fn substeps_performed_js(&self) -> u32 {
        self.debug_info().substeps_performed
    }

    /// Seconds of simulated time the last `step()` call discarded because it
    /// would have needed more than `max_substeps` substeps (the simulation
    /// clock, `time()`, only advances by the time that was integrated).
    #[wasm_bindgen(js_name = "dropped_time")]
    pub fn dropped_time_js(&self) -> f64 {
        self.debug_info().dropped_time
    }

    /// Remove all bodies, colliders and contacts; reset simulation time.
    #[wasm_bindgen(js_name = "reset")]
    pub fn reset_js(&mut self) {
        self.reset();
    }

    // -----------------------------------------------------------------------
    // Body state queries — returns Float64Array via Vec<f64>
    // -----------------------------------------------------------------------

    /// Position of body `handle` as `[x, y, z]`.
    ///
    /// Returns `[0, 0, 0]` for invalid or stale handles (check `contains_body`).
    #[wasm_bindgen(js_name = "get_position")]
    pub fn get_position_js(&self, handle: u32) -> Vec<f64> {
        self.get_position(handle).to_vec()
    }

    /// Orientation quaternion of body `handle` as `[qx, qy, qz, qw]`.
    #[wasm_bindgen(js_name = "get_rotation")]
    pub fn get_rotation_js(&self, handle: u32) -> Vec<f64> {
        self.get_rotation(handle).to_vec()
    }

    /// Linear velocity of body `handle` as `[vx, vy, vz]`.
    #[wasm_bindgen(js_name = "get_velocity")]
    pub fn get_velocity_js(&self, handle: u32) -> Vec<f64> {
        self.get_velocity(handle).to_vec()
    }

    /// Angular velocity of body `handle` as `[wx, wy, wz]` (rad/s).
    #[wasm_bindgen(js_name = "get_angular_velocity")]
    pub fn get_angular_velocity_js(&self, handle: u32) -> Vec<f64> {
        self.get_angular_velocity(handle).to_vec()
    }

    // -----------------------------------------------------------------------
    // Body state mutations
    // -----------------------------------------------------------------------

    /// Set linear velocity of body `handle`. Returns `true` on success and
    /// `false` (leaving the body untouched) for a dead handle or a non-finite
    /// component.
    #[wasm_bindgen(js_name = "set_velocity")]
    pub fn set_velocity_js(&mut self, handle: u32, vx: f64, vy: f64, vz: f64) -> bool {
        self.set_velocity(handle, vx, vy, vz).is_ok()
    }

    /// Set angular velocity of body `handle` in rad/s. Returns `true` on
    /// success and `false` for a dead handle or a non-finite component.
    #[wasm_bindgen(js_name = "set_angular_velocity")]
    pub fn set_angular_velocity_js(&mut self, handle: u32, wx: f64, wy: f64, wz: f64) -> bool {
        self.set_angular_velocity(handle, wx, wy, wz).is_ok()
    }

    /// Teleport body `handle` to `(x, y, z)`. Returns `true` on success and
    /// `false` for a dead handle or a non-finite coordinate.
    #[wasm_bindgen(js_name = "set_position")]
    pub fn set_position_js(&mut self, handle: u32, x: f64, y: f64, z: f64) -> bool {
        self.set_position(handle, x, y, z).is_ok()
    }

    // -----------------------------------------------------------------------
    // Force / impulse
    // -----------------------------------------------------------------------

    /// Apply a world-space force to body `handle` for the current step.
    /// Returns `true` on success and `false` for a dead handle or a
    /// non-finite component.
    #[wasm_bindgen(js_name = "apply_force")]
    pub fn apply_force_js(&mut self, handle: u32, fx: f64, fy: f64, fz: f64) -> bool {
        self.apply_force(handle, fx, fy, fz).is_ok()
    }

    /// Apply a world-space torque to body `handle` for the current step.
    /// Returns `true` on success and `false` for a dead handle or a
    /// non-finite component.
    #[wasm_bindgen(js_name = "apply_torque")]
    pub fn apply_torque_js(&mut self, handle: u32, tx: f64, ty: f64, tz: f64) -> bool {
        self.apply_torque(handle, tx, ty, tz).is_ok()
    }

    /// Apply an instantaneous linear impulse. Returns `true` on success and
    /// `false` for a dead handle or a non-finite component.
    #[wasm_bindgen(js_name = "apply_impulse")]
    pub fn apply_impulse_js(&mut self, handle: u32, ix: f64, iy: f64, iz: f64) -> bool {
        self.apply_impulse(handle, ix, iy, iz).is_ok()
    }

    /// Apply an instantaneous impulse `(ix, iy, iz)` at the world-space offset
    /// `(rx, ry, rz)` from the body centre; the body spins according to its
    /// inertia tensor. Returns `true` on success and `false` for a dead handle
    /// or a non-finite value.
    #[wasm_bindgen(js_name = "apply_impulse_at_point")]
    pub fn apply_impulse_at_point_js(
        &mut self,
        handle: u32,
        ix: f64,
        iy: f64,
        iz: f64,
        rx: f64,
        ry: f64,
        rz: f64,
    ) -> bool {
        self.apply_impulse_at_point(handle, ix, iy, iz, rx, ry, rz)
            .is_ok()
    }

    // -----------------------------------------------------------------------
    // Gravity
    // -----------------------------------------------------------------------

    /// Set the global gravity vector.
    #[wasm_bindgen(js_name = "set_gravity")]
    pub fn set_gravity_js(&mut self, gx: f64, gy: f64, gz: f64) {
        self.set_gravity(gx, gy, gz);
    }

    /// Current gravity as `[gx, gy, gz]`.
    #[wasm_bindgen(js_name = "get_gravity")]
    pub fn get_gravity_js(&self) -> Vec<f64> {
        self.gravity().to_vec()
    }

    // -----------------------------------------------------------------------
    // Statistics
    // -----------------------------------------------------------------------

    /// Total accumulated simulation time in seconds.
    #[wasm_bindgen(js_name = "time")]
    pub fn time_js(&self) -> f64 {
        self.time()
    }

    /// Number of contacts detected in the last `step()`.
    #[wasm_bindgen(js_name = "get_contact_count")]
    pub fn get_contact_count_js(&self) -> u32 {
        self.get_contact_count()
    }

    // -----------------------------------------------------------------------
    // Bulk queries — efficient for rendering
    // -----------------------------------------------------------------------

    /// All active body handles as a `Uint32Array`.
    #[wasm_bindgen(js_name = "get_all_body_handles")]
    pub fn get_all_body_handles_js(&self) -> Vec<u32> {
        self.get_all_body_handles()
    }

    /// All active body positions as a flat `Float64Array` `[x0,y0,z0, x1,y1,z1, ...]`.
    ///
    /// Suitable for instanced mesh updates in Three.js / Babylon.js:
    ///
    /// ```js
    /// const pos = engine.get_all_positions();
    /// for (let i = 0; i < pos.length / 3; i++) {
    ///   mesh.setPositionAt(i, new THREE.Vector3(pos[i*3], pos[i*3+1], pos[i*3+2]));
    /// }
    /// ```
    #[wasm_bindgen(js_name = "get_all_positions")]
    pub fn get_all_positions_js(&self) -> Vec<f64> {
        self.get_all_positions()
    }

    /// All active body transforms as flat `Float64Array` `[x,y,z,qx,qy,qz,qw, ...]`
    /// (7 values per body).
    #[wasm_bindgen(js_name = "get_all_transforms")]
    pub fn get_all_transforms_js(&self) -> Vec<f64> {
        self.get_all_transforms()
    }

    /// Contacts from the last step encoded as a flat `Float64Array`.
    ///
    /// Layout per contact (10 values):
    /// `[body_a, body_b, nx, ny, nz, depth, impulse, px_a, py_a, pz_a]`
    ///
    /// ```js
    /// const c = engine.get_contacts_flat();
    /// for (let i = 0; i < c.length; i += 10) {
    ///   const bodyA = c[i + 0];
    ///   const depth  = c[i + 5];
    /// }
    /// ```
    #[wasm_bindgen(js_name = "get_contacts_flat")]
    pub fn get_contacts_flat_js(&self) -> Vec<f64> {
        let contacts = self.get_contacts();
        let mut out = Vec::with_capacity(contacts.len() * 10);
        for c in contacts {
            out.push(c.body_a as f64);
            out.push(c.body_b as f64);
            out.push(c.normal[0]);
            out.push(c.normal[1]);
            out.push(c.normal[2]);
            out.push(c.depth);
            out.push(c.impulse);
            out.push(c.point_on_a[0]);
            out.push(c.point_on_a[1]);
            out.push(c.point_on_a[2]);
        }
        out
    }

    // -----------------------------------------------------------------------
    // Raycasting
    // -----------------------------------------------------------------------

    /// Cast a ray; returns `[]` (no hit) or `[handle, distance]` (nearest hit).
    ///
    /// ```js
    /// const hit = engine.raycast(ox, oy, oz, dx, dy, dz, maxDist);
    /// if (hit.length > 0) {
    ///   const handle   = hit[0];
    ///   const distance = hit[1];
    /// }
    /// ```
    #[wasm_bindgen(js_name = "raycast")]
    pub fn raycast_js(
        &self,
        ox: f64,
        oy: f64,
        oz: f64,
        dx: f64,
        dy: f64,
        dz: f64,
        max_distance: f64,
    ) -> Vec<f64> {
        let result = self.raycast(ox, oy, oz, dx, dy, dz, max_distance);
        if result.hit {
            vec![result.body_handle as f64, result.distance]
        } else {
            vec![]
        }
    }
}
