// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! `WasmEngine` — extended engine wrapper with performance metrics, debug draw,
//! ray-cast, simulation control, queries, and serialisation helpers.

#![allow(missing_docs)]

use super::WasmPhysicsEngine;
use crate::types::RaycastResult;

// ===========================================================================
// Performance metrics, debug draw, and ray-cast extensions
// ===========================================================================

/// Performance metrics snapshot from the engine.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct PerformanceMetrics {
    /// Estimated frames per second (1 / last_step_time_s).
    pub fps: f64,
    /// Wall-clock time of the last simulation step (seconds).
    pub step_time_s: f64,
    /// Number of active bodies.
    pub body_count: u32,
    /// Number of active colliders.
    pub collider_count: u32,
    /// Number of contacts in the last step.
    pub contact_count: u32,
    /// Total accumulated simulation time.
    pub sim_time: f64,
    /// Integration time of the last step (microseconds).
    pub integration_time_us: u64,
    /// Collision detection time (microseconds).
    pub collision_time_us: u64,
    /// Solver time (microseconds).
    pub solver_time_us: u64,
}

/// Debug visualization flags for the WASM engine.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct DebugDrawFlags {
    /// Draw AABB wireframes for all bodies.
    pub show_aabbs: bool,
    /// Draw velocity vectors.
    pub show_velocities: bool,
    /// Draw contact points.
    pub show_contacts: bool,
    /// Draw center of mass indicators.
    pub show_centers: bool,
    /// Draw body labels (handle IDs).
    pub show_labels: bool,
    /// Global draw scale.
    pub draw_scale: f64,
}

impl DebugDrawFlags {
    /// Create with all flags disabled.
    pub fn disabled() -> Self {
        Self {
            show_aabbs: false,
            show_velocities: false,
            show_contacts: false,
            show_centers: false,
            show_labels: false,
            draw_scale: 1.0,
        }
    }

    /// Create with all flags enabled.
    pub fn all_enabled() -> Self {
        Self {
            show_aabbs: true,
            show_velocities: true,
            show_contacts: true,
            show_centers: true,
            show_labels: true,
            draw_scale: 1.0,
        }
    }

    /// Returns `true` if any flag is enabled.
    pub fn any_enabled(&self) -> bool {
        self.show_aabbs
            || self.show_velocities
            || self.show_contacts
            || self.show_centers
            || self.show_labels
    }
}

/// Extended engine wrapper exposing performance metrics, debug draw, and ray-cast.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct WasmEngine {
    /// The underlying WASM physics engine.
    engine: WasmPhysicsEngine,
    /// Current debug draw flags.
    debug_flags: DebugDrawFlags,
    /// Simulated last-step wall-clock duration (seconds).
    /// In native tests we store a fixed value; in real WASM this would be measured.
    last_step_time_s: f64,
}

impl WasmEngine {
    /// Create a new `WasmEngine` wrapping a fresh `WasmPhysicsEngine`.
    pub fn new(gx: f64, gy: f64, gz: f64) -> Self {
        Self {
            engine: WasmPhysicsEngine::new(gx, gy, gz),
            debug_flags: DebugDrawFlags::disabled(),
            last_step_time_s: 1.0 / 60.0,
        }
    }

    /// Advance the simulation and record step timing.
    pub fn step(&mut self, dt: f64) {
        self.engine.step(dt);
        self.last_step_time_s = dt.max(1e-6);
    }

    /// Add a dynamic body at `(x, y, z)` with the given mass.
    pub fn add_dynamic_body(&mut self, mass: f64, x: f64, y: f64, z: f64) -> u32 {
        self.engine.add_dynamic_body(mass, x, y, z)
    }

    /// Add a sphere collider to `body_id`.
    pub fn add_sphere_collider(&mut self, body_id: u32, radius: f64) -> u32 {
        self.engine.add_sphere_collider(body_id, radius)
    }

    // -----------------------------------------------------------------------
    // Performance metrics
    // -----------------------------------------------------------------------

    /// Return a snapshot of current engine performance metrics.
    pub fn get_performance_metrics(&self) -> PerformanceMetrics {
        let di = self.engine.debug_info();
        let body_count = self.engine.get_body_count();
        let contact_count = self.engine.get_contacts().len() as u32;
        // Count active colliders by querying the engine's active collider count
        let collider_count = self.engine.get_active_collider_count();
        let step_time_s = self.last_step_time_s.max(1e-9);
        PerformanceMetrics {
            fps: 1.0 / step_time_s,
            step_time_s,
            body_count,
            collider_count,
            contact_count,
            sim_time: self.engine.time(),
            integration_time_us: di.integration_time_us,
            collision_time_us: di.collision_time_us,
            solver_time_us: di.solver_time_us,
        }
    }

    // -----------------------------------------------------------------------
    // Debug draw
    // -----------------------------------------------------------------------

    /// Set the debug visualization flags for this engine.
    pub fn set_debug_draw(&mut self, flags: DebugDrawFlags) {
        self.debug_flags = flags;
    }

    /// Return the current debug draw flags.
    pub fn debug_draw_flags(&self) -> &DebugDrawFlags {
        &self.debug_flags
    }

    /// Toggle a specific debug flag by name.
    ///
    /// Recognized names: `"aabbs"`, `"velocities"`, `"contacts"`, `"centers"`, `"labels"`.
    /// Unknown names are ignored.
    pub fn toggle_debug_flag(&mut self, flag: &str, enabled: bool) {
        match flag {
            "aabbs" => self.debug_flags.show_aabbs = enabled,
            "velocities" => self.debug_flags.show_velocities = enabled,
            "contacts" => self.debug_flags.show_contacts = enabled,
            "centers" => self.debug_flags.show_centers = enabled,
            "labels" => self.debug_flags.show_labels = enabled,
            _ => {}
        }
    }

    // -----------------------------------------------------------------------
    // Ray cast
    // -----------------------------------------------------------------------

    /// Cast a ray from `origin` in `direction` and return the closest hit.
    ///
    /// Delegates to `WasmPhysicsEngine::raycast`. `max_dist` limits the search.
    #[allow(clippy::too_many_arguments)]
    pub fn compute_ray_cast(
        &self,
        ox: f64,
        oy: f64,
        oz: f64,
        dx: f64,
        dy: f64,
        dz: f64,
        max_dist: f64,
    ) -> RaycastResult {
        self.engine.raycast(ox, oy, oz, dx, dy, dz, max_dist)
    }

    /// Cast a ray and return a flat representation `[hit, body_handle, px, py, pz, nx, ny, nz, dist]`.
    ///
    /// `hit` is 1.0 if something was hit, 0.0 otherwise.
    #[allow(clippy::too_many_arguments)]
    pub fn compute_ray_cast_flat(
        &self,
        ox: f64,
        oy: f64,
        oz: f64,
        dx: f64,
        dy: f64,
        dz: f64,
        max_dist: f64,
    ) -> Vec<f64> {
        let r = self.compute_ray_cast(ox, oy, oz, dx, dy, dz, max_dist);
        vec![
            if r.hit { 1.0 } else { 0.0 },
            r.body_handle as f64,
            r.point[0],
            r.point[1],
            r.point[2],
            r.normal[0],
            r.normal[1],
            r.normal[2],
            r.distance,
        ]
    }
}

// ===========================================================================
// Extended simulation control (JS-callable)
// ===========================================================================

/// Extended body snapshot returned by `get_body_snapshot`.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct BodySnapshot {
    /// Body handle.
    pub handle: u32,
    /// World-space position `[x, y, z]`.
    pub position: [f64; 3],
    /// Orientation quaternion `[x, y, z, w]`.
    pub rotation: [f64; 4],
    /// Linear velocity `[vx, vy, vz]`.
    pub velocity: [f64; 3],
    /// Angular velocity `[wx, wy, wz]`.
    pub angular_velocity: [f64; 3],
    /// Whether the body is asleep.
    pub sleeping: bool,
    /// Kinetic energy (translational + rotational approximation).
    pub kinetic_energy: f64,
    /// Speed (magnitude of linear velocity).
    pub speed: f64,
}

impl BodySnapshot {
    /// Build from the engine.
    fn from_engine(engine: &WasmPhysicsEngine, handle: u32) -> Option<Self> {
        let state = engine.get_body_state(handle)?;
        let v = state.linear_velocity;
        let speed = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        let ke = state.kinetic_energy;
        Some(Self {
            handle,
            position: state.position,
            rotation: state.rotation,
            velocity: state.linear_velocity,
            angular_velocity: state.angular_velocity,
            sleeping: state.is_sleeping,
            kinetic_energy: ke,
            speed,
        })
    }
}

/// Simulation summary returned by `simulation_summary`.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct SimulationSummary {
    /// Total elapsed simulation time.
    pub elapsed_time: f64,
    /// Number of active bodies.
    pub active_body_count: u32,
    /// Number of contacts in the last step.
    pub contact_count: u32,
    /// Total kinetic energy of all dynamic bodies.
    pub total_kinetic_energy: f64,
    /// Average speed across all dynamic bodies.
    pub average_speed: f64,
    /// Maximum speed of any body.
    pub max_speed: f64,
    /// Gravity vector.
    pub gravity: [f64; 3],
}

/// AABB query result.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct AabbQuery {
    /// Handles of bodies whose positions fall within the AABB.
    pub handles: Vec<u32>,
    /// Number of matching bodies.
    pub count: usize,
}

/// Sphere overlap query result.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct SphereOverlapResult {
    /// Handles of bodies whose positions fall within the sphere.
    pub handles: Vec<u32>,
    /// Distances from the query center to each matching body.
    pub distances: Vec<f64>,
}

/// JSON-serialisable flat body state for WASM boundary crossing.
///
/// Flat layout:  `[handle, px, py, pz, rx, ry, rz, rw, vx, vy, vz, wx, wy, wz, speed]`
pub const BODY_STATE_FLAT_LEN: usize = 15;

impl WasmEngine {
    // ------------------------------------------------------------------
    // Simulation control
    // ------------------------------------------------------------------

    /// Step the simulation `n` times with the given `dt`.
    ///
    /// Equivalent to calling `step(dt)` in a loop but avoids repeated
    /// JS->WASM boundary crossings.
    pub fn step_n(&mut self, n: u32, dt: f64) {
        for _ in 0..n {
            self.step(dt);
        }
    }

    /// Reset the engine to an empty state and set new gravity.
    pub fn reset_with_new_gravity(&mut self, gx: f64, gy: f64, gz: f64) {
        self.engine.reset_with_gravity(gx, gy, gz);
    }

    /// Change gravity without resetting the simulation.
    pub fn set_gravity(&mut self, gx: f64, gy: f64, gz: f64) {
        self.engine.set_gravity(gx, gy, gz);
    }

    /// Return the current gravity vector as `[gx, gy, gz]`.
    pub fn get_gravity(&self) -> [f64; 3] {
        self.engine.gravity()
    }

    /// Add a dynamic body and attach a sphere collider in one call.
    ///
    /// Returns the body handle.
    pub fn add_dynamic_sphere(&mut self, mass: f64, x: f64, y: f64, z: f64, radius: f64) -> u32 {
        let handle = self.engine.add_dynamic_body(mass, x, y, z);
        self.engine.add_sphere_collider(handle, radius);
        handle
    }

    /// Add a static body and attach a plane collider in one call.
    ///
    /// Returns the body handle.
    #[allow(clippy::too_many_arguments)]
    pub fn add_static_plane(
        &mut self,
        x: f64,
        y: f64,
        z: f64,
        nx: f64,
        ny: f64,
        nz: f64,
        offset: f64,
    ) -> u32 {
        let handle = self.engine.add_static_body(x, y, z);
        self.engine.add_plane_collider(handle, nx, ny, nz, offset);
        handle
    }

    /// Remove a body by handle; returns `true` if it existed.
    pub fn remove_body(&mut self, handle: u32) -> bool {
        self.engine.remove_body(handle).is_ok()
    }

    /// Set the position of a body.  Returns `true` on success.
    pub fn set_body_position(&mut self, handle: u32, x: f64, y: f64, z: f64) -> bool {
        self.engine.set_position(handle, x, y, z).is_ok()
    }

    /// Set the linear velocity of a body.  Returns `true` on success.
    pub fn set_body_velocity(&mut self, handle: u32, vx: f64, vy: f64, vz: f64) -> bool {
        self.engine.set_velocity(handle, vx, vy, vz).is_ok()
    }

    /// Apply a force to a body for the current step.  Returns `true` on success.
    pub fn apply_body_force(&mut self, handle: u32, fx: f64, fy: f64, fz: f64) -> bool {
        self.engine.apply_force(handle, fx, fy, fz).is_ok()
    }

    /// Apply an impulse to a body.  Returns `true` on success.
    pub fn apply_body_impulse(&mut self, handle: u32, ix: f64, iy: f64, iz: f64) -> bool {
        self.engine.apply_impulse(handle, ix, iy, iz).is_ok()
    }

    /// Apply a torque to a body for the current step.  Returns `true` on success.
    pub fn apply_body_torque(&mut self, handle: u32, tx: f64, ty: f64, tz: f64) -> bool {
        self.engine.apply_torque(handle, tx, ty, tz).is_ok()
    }

    /// Set the linear damping coefficient of a body.
    pub fn set_linear_damping(&mut self, handle: u32, damping: f64) -> bool {
        self.engine.set_body_linear_damping(handle, damping).is_ok()
    }

    /// Set the angular damping coefficient of a body.
    pub fn set_angular_damping(&mut self, handle: u32, damping: f64) -> bool {
        self.engine
            .set_body_angular_damping(handle, damping)
            .is_ok()
    }

    // ------------------------------------------------------------------
    // Physics query functions
    // ------------------------------------------------------------------

    /// Return a `BodySnapshot` for the given handle, or `None` if invalid.
    pub fn get_body_snapshot(&self, handle: u32) -> Option<BodySnapshot> {
        BodySnapshot::from_engine(&self.engine, handle)
    }

    /// Return all body handles as a `Vec`u32`.
    pub fn get_all_handles(&self) -> Vec<u32> {
        self.engine.get_all_body_handles()
    }

    /// Return all body positions as flat `\[x0,y0,z0, x1,y1,z1, ...\]`.
    pub fn get_all_positions_flat(&self) -> Vec<f64> {
        self.engine.get_all_positions()
    }

    /// Return all body transforms as flat `\[px,py,pz, rx,ry,rz,rw, ...\]` (7 per body).
    pub fn get_all_transforms_flat(&self) -> Vec<f64> {
        self.engine.get_all_transforms()
    }

    /// Return position of a specific body as `\[x, y, z\]`.
    pub fn get_body_position(&self, handle: u32) -> [f64; 3] {
        self.engine.get_position(handle)
    }

    /// Return linear velocity of a specific body as `\[vx, vy, vz\]`.
    pub fn get_body_velocity(&self, handle: u32) -> [f64; 3] {
        self.engine.get_velocity(handle)
    }

    /// Return speed (velocity magnitude) of a specific body.
    pub fn get_body_speed(&self, handle: u32) -> f64 {
        let v = self.engine.get_velocity(handle);
        (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
    }

    /// Return the total number of active bodies.
    pub fn body_count(&self) -> u32 {
        self.engine.get_body_count()
    }

    /// Return a complete `SimulationSummary` for the current state.
    pub fn simulation_summary(&self) -> SimulationSummary {
        let handles = self.engine.get_all_body_handles();
        let mut total_ke = 0.0f64;
        let mut total_speed = 0.0f64;
        let mut max_speed = 0.0f64;
        let mut dynamic_count = 0u32;

        for h in &handles {
            if let Some(state) = self.engine.get_body_state(*h) {
                let v = state.linear_velocity;
                let spd = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
                // Use pre-computed kinetic energy from BodyState
                if state.kinetic_energy > 0.0 || spd > 0.0 {
                    total_ke += state.kinetic_energy;
                    total_speed += spd;
                    if spd > max_speed {
                        max_speed = spd;
                    }
                    dynamic_count += 1;
                }
            }
        }

        let avg_speed = if dynamic_count > 0 {
            total_speed / dynamic_count as f64
        } else {
            0.0
        };

        SimulationSummary {
            elapsed_time: self.engine.time(),
            active_body_count: self.engine.get_body_count(),
            contact_count: self.engine.get_contact_count(),
            total_kinetic_energy: total_ke,
            average_speed: avg_speed,
            max_speed,
            gravity: self.engine.gravity(),
        }
    }

    /// Return bodies within an axis-aligned bounding box.
    ///
    /// `aabb` is `\[min_x, min_y, min_z, max_x, max_y, max_z\]`.
    pub fn query_aabb(&self, aabb: [f64; 6]) -> AabbQuery {
        let handles: Vec<u32> = self
            .engine
            .get_all_body_handles()
            .into_iter()
            .filter(|&h| {
                let p = self.engine.get_position(h);
                p[0] >= aabb[0]
                    && p[0] <= aabb[3]
                    && p[1] >= aabb[1]
                    && p[1] <= aabb[4]
                    && p[2] >= aabb[2]
                    && p[2] <= aabb[5]
            })
            .collect();
        let count = handles.len();
        AabbQuery { handles, count }
    }

    /// Return bodies within a sphere centered at `(cx, cy, cz)` with `radius`.
    pub fn query_sphere_overlap(
        &self,
        cx: f64,
        cy: f64,
        cz: f64,
        radius: f64,
    ) -> SphereOverlapResult {
        let r2 = radius * radius;
        let mut handles = Vec::new();
        let mut distances = Vec::new();
        for h in self.engine.get_all_body_handles() {
            let p = self.engine.get_position(h);
            let dx = p[0] - cx;
            let dy = p[1] - cy;
            let dz = p[2] - cz;
            let d2 = dx * dx + dy * dy + dz * dz;
            if d2 <= r2 {
                handles.push(h);
                distances.push(d2.sqrt());
            }
        }
        SphereOverlapResult { handles, distances }
    }

    /// Return the handle of the body closest to the given point, or `None` if empty.
    pub fn find_closest_body(&self, x: f64, y: f64, z: f64) -> Option<u32> {
        let handles = self.engine.get_all_body_handles();
        if handles.is_empty() {
            return None;
        }
        let mut best_handle = handles[0];
        let mut best_d2 = f64::MAX;
        for h in handles {
            let p = self.engine.get_position(h);
            let dx = p[0] - x;
            let dy = p[1] - y;
            let dz = p[2] - z;
            let d2 = dx * dx + dy * dy + dz * dz;
            if d2 < best_d2 {
                best_d2 = d2;
                best_handle = h;
            }
        }
        Some(best_handle)
    }

    /// Return all current contacts as a flat array.
    ///
    /// Each contact is encoded as:
    /// `\[body_a, body_b, px, py, pz, nx, ny, nz, depth\]` (9 elements).
    pub fn get_contacts_flat_extended(&self) -> Vec<f64> {
        let contacts = self.engine.get_contacts();
        let mut out = Vec::with_capacity(contacts.len() * 9);
        for c in contacts {
            // Use point_on_a as representative contact point
            let cp = c.contact_midpoint();
            out.push(c.body_a as f64);
            out.push(c.body_b as f64);
            out.push(cp[0]);
            out.push(cp[1]);
            out.push(cp[2]);
            out.push(c.normal[0]);
            out.push(c.normal[1]);
            out.push(c.normal[2]);
            out.push(c.depth);
        }
        out
    }

    /// Return the number of active contacts.
    pub fn contact_count(&self) -> u32 {
        self.engine.get_contact_count()
    }

    // ------------------------------------------------------------------
    // Serialisation / deserialisation helpers
    // ------------------------------------------------------------------

    /// Serialise a body's state to a flat array of `BODY_STATE_FLAT_LEN` elements.
    ///
    /// Layout: `\[handle, px, py, pz, rx, ry, rz, rw, vx, vy, vz, wx, wy, wz, speed\]`
    ///
    /// Returns an array filled with `f64::NAN` if the handle is invalid.
    pub fn serialise_body_flat(&self, handle: u32) -> [f64; BODY_STATE_FLAT_LEN] {
        match BodySnapshot::from_engine(&self.engine, handle) {
            Some(s) => [
                handle as f64,
                s.position[0],
                s.position[1],
                s.position[2],
                s.rotation[0],
                s.rotation[1],
                s.rotation[2],
                s.rotation[3],
                s.velocity[0],
                s.velocity[1],
                s.velocity[2],
                s.angular_velocity[0],
                s.angular_velocity[1],
                s.angular_velocity[2],
                s.speed,
            ],
            None => [f64::NAN; BODY_STATE_FLAT_LEN],
        }
    }

    /// Serialise all body states to a flat array.
    ///
    /// Each body occupies `BODY_STATE_FLAT_LEN` consecutive elements.
    pub fn serialise_all_bodies_flat(&self) -> Vec<f64> {
        let handles = self.engine.get_all_body_handles();
        let mut out = Vec::with_capacity(handles.len() * BODY_STATE_FLAT_LEN);
        for h in handles {
            out.extend_from_slice(&self.serialise_body_flat(h));
        }
        out
    }

    /// Deserialise a flat body record and update that body's position and
    /// velocity.
    ///
    /// `data` must be at least `BODY_STATE_FLAT_LEN` elements long.
    /// Returns `true` if the update succeeded.
    pub fn deserialise_body_flat(&mut self, data: &[f64]) -> bool {
        if data.len() < BODY_STATE_FLAT_LEN {
            return false;
        }
        let handle = data[0] as u32;
        let ok_pos = self
            .engine
            .set_position(handle, data[1], data[2], data[3])
            .is_ok();
        let ok_vel = self
            .engine
            .set_velocity(handle, data[8], data[9], data[10])
            .is_ok();
        ok_pos && ok_vel
    }

    /// Serialise the entire simulation state to a JSON-compatible `String`.
    ///
    /// Format: `{"time": f64, "bodies": \[flat_array, ...\], "gravity": \[gx,gy,gz\]}`
    pub fn state_to_json_extended(&self) -> String {
        let handles = self.engine.get_all_body_handles();
        let g = self.engine.gravity();
        let t = self.engine.time();

        let bodies_json: Vec<String> = handles
            .iter()
            .map(|&h| {
                let s = self.serialise_body_flat(h);
                format!(
                    "{{\"handle\":{},\"pos\":[{},{},{}],\"vel\":[{},{},{}],\"speed\":{}}}",
                    s[0] as u32, s[1], s[2], s[3], s[8], s[9], s[10], s[14]
                )
            })
            .collect();

        format!(
            "{{\"time\":{},\"gravity\":[{},{},{}],\"bodies\":[{}]}}",
            t,
            g[0],
            g[1],
            g[2],
            bodies_json.join(",")
        )
    }

    /// Serialise contact list to a JSON string.
    pub fn contacts_to_json(&self) -> String {
        let contacts = self.engine.get_contacts();
        let entries: Vec<String> = contacts
            .iter()
            .map(|c| {
                format!(
                    "{{\"a\":{},\"b\":{},\"depth\":{},\"n\":[{},{},{}]}}",
                    c.body_a, c.body_b, c.depth, c.normal[0], c.normal[1], c.normal[2]
                )
            })
            .collect();
        format!("[{}]", entries.join(","))
    }

    /// Return total elapsed simulation time.
    pub fn elapsed_time(&self) -> f64 {
        self.engine.time()
    }

    /// Return `true` if there are no active bodies.
    pub fn is_empty(&self) -> bool {
        self.engine.get_body_count() == 0
    }
}
