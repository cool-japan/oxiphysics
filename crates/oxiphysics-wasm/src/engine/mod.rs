// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Core WASM physics engine (`WasmPhysicsEngine`).
//!
//! This module provides a self-contained, pure-Rust physics engine designed
//! for use at the WASM boundary. It does **not** depend on nalgebra directly;
//! all math uses flat `[f64; N]` arrays or the types defined in `crate::types`.
//!
//! ## Architecture
//!
//! The engine maintains:
//! - A `Vec`BodyEntry` of body states (position, velocity, forces, etc.)
//! - A `Vec<ColliderEntry>` of collider shapes attached to bodies
//! - A `Vec<ContactResult>` list refreshed each step
//! - A `Vec`u32` free-list for handle reuse
//!
//! Integration uses symplectic Euler (velocity Verlet variant) which is
//! stable for small time steps and straightforward to implement.
//!
//! ## Example
//!
//! ```no_run
//! use oxiphysics_wasm::WasmPhysicsEngine;
//!
//! let mut engine = WasmPhysicsEngine::new(0.0, -9.81, 0.0);
//! let body = engine.add_dynamic_body(1.0, 0.0, 10.0, 0.0);
//! engine.step(1.0 / 60.0);
//! let pos = engine.get_position(body);
//! assert!(pos[1] < 10.0); // body fell due to gravity
//! ```

#![allow(missing_docs)]

// Sub-modules
mod bindings;
mod constraint;
mod float64_view;
mod lbm;
mod sph;
mod vehicle;
mod wasm_engine;
mod wasm_types;

#[cfg(test)]
mod tests;

// Re-export everything from submodules
pub use bindings::{ContactInfoEntry, WasmBindings, WasmContactList};
pub use constraint::{ContactData, WasmConstraint, WasmConstraintSolver};
pub use float64_view::Float64View;
pub use lbm::{WasmLbmConfig, WasmLbmSim};
pub use sph::WasmSphSim;
pub use vehicle::{VehicleState, WasmVehicleSim};
pub use wasm_engine::{
    AabbQuery, BODY_STATE_FLAT_LEN, BodySnapshot, DebugDrawFlags, PerformanceMetrics,
    SimulationSummary, SphereOverlapResult, WasmEngine,
};
pub use wasm_types::{
    TS_CONTACT_DEF, TS_DEFINITIONS, TS_TRANSFORM_DEF, TS_VEC3_DEF, WASM_PAGE_SIZE, WasmTransform,
    WasmVec3, error_to_js_string, free_f64_buffer, free_u8_buffer, free_u32_buffer, result_to_js,
};

use crate::error::{Error, Result};
use crate::types::{
    BodyState, ColliderConfig, ColliderShapeType, ContactResult, DebugInfo, RaycastResult,
    RigidBodyConfig, SimulationConfig, Vec3Wasm,
};

// ---------------------------------------------------------------------------
// Internal body storage
// ---------------------------------------------------------------------------

/// Internal body state stored in the engine.
#[derive(Debug, Clone)]
struct BodyEntry {
    /// Current position.
    position: [f64; 3],
    /// Current orientation quaternion `[x, y, z, w]`.
    rotation: [f64; 4],
    /// Current linear velocity.
    velocity: [f64; 3],
    /// Current angular velocity.
    angular_velocity: [f64; 3],
    /// Accumulated force for this step.
    force: [f64; 3],
    /// Accumulated torque for this step.
    torque: [f64; 3],
    /// Mass (0 for static).
    mass: f64,
    /// Inverse mass (0 for static).
    inv_mass: f64,
    /// Linear damping.
    linear_damping: f64,
    /// Angular damping.
    angular_damping: f64,
    /// Restitution.
    restitution: f64,
    /// Friction coefficient.
    friction: f64,
    /// Gravity scale multiplier.
    gravity_scale: f64,
    /// Whether this entry is occupied (vs. a freed slot).
    active: bool,
    /// Whether this body is static.
    is_static: bool,
    /// Whether this body is kinematic.
    is_kinematic: bool,
    /// Whether this body is a sensor.
    #[allow(dead_code)]
    is_sensor: bool,
    /// Whether the body is currently sleeping.
    sleeping: bool,
    /// Time the body has been below sleep thresholds (seconds).
    sleep_timer: f64,
    /// Generation counter for safe handle reuse.
    generation: u32,
}

impl BodyEntry {
    fn from_config(cfg: &RigidBodyConfig, generation: u32) -> Self {
        let inv_mass = if cfg.mass > 0.0 { 1.0 / cfg.mass } else { 0.0 };
        Self {
            position: cfg.position,
            rotation: cfg.rotation,
            velocity: cfg.linear_velocity,
            angular_velocity: cfg.angular_velocity,
            force: [0.0; 3],
            torque: [0.0; 3],
            mass: cfg.mass,
            inv_mass,
            linear_damping: cfg.linear_damping,
            angular_damping: cfg.angular_damping,
            restitution: cfg.restitution,
            friction: cfg.friction,
            gravity_scale: cfg.gravity_scale,
            active: true,
            is_static: cfg.is_static(),
            is_kinematic: cfg.is_kinematic,
            is_sensor: cfg.is_sensor,
            sleeping: false,
            sleep_timer: 0.0,
            generation,
        }
    }

    /// Returns `true` if the body is movable (dynamic, not static or kinematic).
    fn is_dynamic(&self) -> bool {
        !self.is_static && !self.is_kinematic && self.mass > 0.0
    }

    /// Kinetic energy (translational only).
    fn kinetic_energy(&self) -> f64 {
        if self.mass <= 0.0 {
            return 0.0;
        }
        let v2 = self.velocity[0] * self.velocity[0]
            + self.velocity[1] * self.velocity[1]
            + self.velocity[2] * self.velocity[2];
        0.5 * self.mass * v2
    }
}

// ---------------------------------------------------------------------------
// Internal collider storage
// ---------------------------------------------------------------------------

/// Internal collider state stored in the engine.
#[derive(Debug, Clone)]
struct ColliderEntry {
    /// Handle of the owning body.
    body_index: u32,
    /// Shape configuration.
    config: ColliderConfig,
    /// Whether this slot is occupied.
    active: bool,
    /// Generation counter.
    generation: u32,
}

// ---------------------------------------------------------------------------
// Contact manifold
// ---------------------------------------------------------------------------

/// Simple axis-aligned sphere-sphere contact detection result.
#[derive(Debug, Clone)]
struct ManifoldPair {
    body_a: u32,
    body_b: u32,
    normal: [f64; 3],
    depth: f64,
    point_a: [f64; 3],
    point_b: [f64; 3],
}

// ---------------------------------------------------------------------------
// WasmPhysicsEngine
// ---------------------------------------------------------------------------

/// A self-contained physics engine designed for WASM/JavaScript interop.
///
/// All methods use primitive types (`f64`, `u32`) and flat `[f64; N]` arrays
/// to minimize overhead at the WASM boundary. No nalgebra types are exposed.
///
/// ## Lifecycle
///
/// 1. Create with [`WasmPhysicsEngine::new`] or [`WasmPhysicsEngine::from_config`].
/// 2. Add bodies with [`WasmPhysicsEngine::add_rigid_body`], [`WasmPhysicsEngine::add_dynamic_body`], [`WasmPhysicsEngine::add_static_body`].
/// 3. Add colliders with [`WasmPhysicsEngine::add_sphere_collider`], [`WasmPhysicsEngine::add_box_collider`].
/// 4. Advance simulation with [`WasmPhysicsEngine::step`].
/// 5. Query state with [`WasmPhysicsEngine::get_position`], [`WasmPhysicsEngine::get_velocity`], [`WasmPhysicsEngine::get_contacts`], etc.
/// 6. Remove bodies with [`WasmPhysicsEngine::remove_body`] or reset everything with [`WasmPhysicsEngine::reset`].
#[derive(Debug, Clone)]
pub struct WasmPhysicsEngine {
    /// Gravity components `[gx, gy, gz]`.
    gravity: [f64; 3],
    /// All body slots (may contain freed entries).
    bodies: Vec<BodyEntry>,
    /// Free list of body indices available for reuse.
    body_free_list: Vec<usize>,
    /// All collider slots (may contain freed entries).
    colliders: Vec<ColliderEntry>,
    /// Free list of collider indices available for reuse.
    collider_free_list: Vec<usize>,
    /// Contact results from the last call to `step()`.
    contacts: Vec<ContactResult>,
    /// Total accumulated simulation time (seconds).
    time: f64,
    /// Simulation configuration (time step, solver settings, etc.).
    config: SimulationConfig,
    /// Debug/performance information from the last step.
    debug_info: DebugInfo,
}

impl WasmPhysicsEngine {
    // -----------------------------------------------------------------------
    // Construction
    // -----------------------------------------------------------------------

    /// Create a new engine with the given gravity components.
    ///
    /// # Example
    ///
    /// ```no_run
    /// use oxiphysics_wasm::WasmPhysicsEngine;
    /// let engine = WasmPhysicsEngine::new(0.0, -9.81, 0.0);
    /// assert_eq!(engine.get_body_count(), 0);
    /// ```
    pub fn new(gx: f64, gy: f64, gz: f64) -> Self {
        let config = SimulationConfig {
            gravity: [gx, gy, gz],
            ..Default::default()
        };
        Self {
            gravity: [gx, gy, gz],
            bodies: Vec::new(),
            body_free_list: Vec::new(),
            colliders: Vec::new(),
            collider_free_list: Vec::new(),
            contacts: Vec::new(),
            time: 0.0,
            config,
            debug_info: DebugInfo::default(),
        }
    }

    /// Create an engine from a `SimulationConfig`.
    pub fn from_config(config: SimulationConfig) -> Self {
        let gravity = config.gravity;
        Self {
            gravity,
            bodies: Vec::new(),
            body_free_list: Vec::new(),
            colliders: Vec::new(),
            collider_free_list: Vec::new(),
            contacts: Vec::new(),
            time: 0.0,
            config,
            debug_info: DebugInfo::default(),
        }
    }

    // -----------------------------------------------------------------------
    // Body management
    // -----------------------------------------------------------------------

    /// Add a rigid body using a [`RigidBodyConfig`] and return its handle.
    pub fn add_rigid_body(&mut self, config: &RigidBodyConfig) -> u32 {
        if let Some(idx) = self.body_free_list.pop() {
            let old_gen = self.bodies[idx].generation;
            self.bodies[idx] = BodyEntry::from_config(config, old_gen + 1);
            idx as u32
        } else {
            let idx = self.bodies.len();
            self.bodies.push(BodyEntry::from_config(config, 0));
            idx as u32
        }
    }

    /// Add a dynamic (movable) body and return its handle.
    pub fn add_dynamic_body(&mut self, mass: f64, x: f64, y: f64, z: f64) -> u32 {
        let cfg = RigidBodyConfig::new()
            .with_mass(mass)
            .with_position(x, y, z);
        self.add_rigid_body(&cfg)
    }

    /// Add a static (immovable) body and return its handle.
    pub fn add_static_body(&mut self, x: f64, y: f64, z: f64) -> u32 {
        let cfg = RigidBodyConfig::static_body().with_position(x, y, z);
        self.add_rigid_body(&cfg)
    }

    /// Remove the body with the given handle.
    pub fn remove_body(&mut self, handle: u32) -> Result<()> {
        let idx = handle as usize;
        match self.bodies.get_mut(idx) {
            Some(body) if body.active => {
                body.active = false;
                for col in &mut self.colliders {
                    if col.active && col.body_index == handle {
                        col.active = false;
                    }
                }
                self.body_free_list.push(idx);
                Ok(())
            }
            Some(_) => Err(Error::InvalidHandle(handle)),
            None => Err(Error::InvalidHandle(handle)),
        }
    }

    /// Return the number of **active** bodies (not counting freed slots).
    pub fn get_body_count(&self) -> u32 {
        self.bodies.iter().filter(|b| b.active).count() as u32
    }

    /// Return the number of **active** colliders.
    pub fn get_active_collider_count(&self) -> u32 {
        self.colliders.iter().filter(|c| c.active).count() as u32
    }

    // -----------------------------------------------------------------------
    // Velocity & force control
    // -----------------------------------------------------------------------

    /// Set the linear velocity of a body.
    pub fn set_velocity(&mut self, handle: u32, vx: f64, vy: f64, vz: f64) -> Result<()> {
        self.get_body_mut(handle).map(|b| {
            b.velocity = [vx, vy, vz];
            b.sleeping = false;
            b.sleep_timer = 0.0;
        })
    }

    /// Set the angular velocity of a body (radians/second).
    pub fn set_angular_velocity(&mut self, handle: u32, wx: f64, wy: f64, wz: f64) -> Result<()> {
        self.get_body_mut(handle).map(|b| {
            b.angular_velocity = [wx, wy, wz];
            b.sleeping = false;
            b.sleep_timer = 0.0;
        })
    }

    /// Set the position of a body directly (teleport).
    pub fn set_position(&mut self, handle: u32, x: f64, y: f64, z: f64) -> Result<()> {
        self.get_body_mut(handle).map(|b| {
            b.position = [x, y, z];
        })
    }

    /// Apply a world-space force to a body for the current step.
    pub fn apply_force(&mut self, handle: u32, fx: f64, fy: f64, fz: f64) -> Result<()> {
        self.get_body_mut(handle).map(|b| {
            b.force[0] += fx;
            b.force[1] += fy;
            b.force[2] += fz;
            b.sleeping = false;
            b.sleep_timer = 0.0;
        })
    }

    /// Apply a world-space torque to a body for the current step.
    pub fn apply_torque(&mut self, handle: u32, tx: f64, ty: f64, tz: f64) -> Result<()> {
        self.get_body_mut(handle).map(|b| {
            b.torque[0] += tx;
            b.torque[1] += ty;
            b.torque[2] += tz;
            b.sleeping = false;
        })
    }

    /// Apply an instantaneous linear impulse to a body.
    pub fn apply_impulse(&mut self, handle: u32, ix: f64, iy: f64, iz: f64) -> Result<()> {
        self.get_body_mut(handle).map(|b| {
            if b.inv_mass > 0.0 {
                b.velocity[0] += ix * b.inv_mass;
                b.velocity[1] += iy * b.inv_mass;
                b.velocity[2] += iz * b.inv_mass;
                b.sleeping = false;
                b.sleep_timer = 0.0;
            }
        })
    }

    /// Apply an impulse at a point relative to the body center (generates torque).
    #[allow(clippy::too_many_arguments)]
    pub fn apply_impulse_at_point(
        &mut self,
        handle: u32,
        ix: f64,
        iy: f64,
        iz: f64,
        rx: f64,
        ry: f64,
        rz: f64,
    ) -> Result<()> {
        self.get_body_mut(handle).map(|b| {
            if b.inv_mass > 0.0 {
                b.velocity[0] += ix * b.inv_mass;
                b.velocity[1] += iy * b.inv_mass;
                b.velocity[2] += iz * b.inv_mass;
                b.angular_velocity[0] += ry * iz - rz * iy;
                b.angular_velocity[1] += rz * ix - rx * iz;
                b.angular_velocity[2] += rx * iy - ry * ix;
                b.sleeping = false;
                b.sleep_timer = 0.0;
            }
        })
    }

    // -----------------------------------------------------------------------
    // Queries
    // -----------------------------------------------------------------------

    /// Get the position of a body as `[x, y, z]`.
    pub fn get_position(&self, handle: u32) -> [f64; 3] {
        self.get_body(handle)
            .map(|b| b.position)
            .unwrap_or([0.0; 3])
    }

    /// Get the orientation quaternion `[x, y, z, w]` of a body.
    pub fn get_rotation(&self, handle: u32) -> [f64; 4] {
        self.get_body(handle)
            .map(|b| b.rotation)
            .unwrap_or([0.0, 0.0, 0.0, 1.0])
    }

    /// Get the linear velocity of a body as `[vx, vy, vz]`.
    pub fn get_velocity(&self, handle: u32) -> [f64; 3] {
        self.get_body(handle)
            .map(|b| b.velocity)
            .unwrap_or([0.0; 3])
    }

    /// Get the angular velocity of a body as `[wx, wy, wz]` (rad/s).
    pub fn get_angular_velocity(&self, handle: u32) -> [f64; 3] {
        self.get_body(handle)
            .map(|b| b.angular_velocity)
            .unwrap_or([0.0; 3])
    }

    /// Get the full runtime state of a body.
    pub fn get_body_state(&self, handle: u32) -> Option<BodyState> {
        self.get_body(handle).map(|b| BodyState {
            handle,
            position: b.position,
            rotation: b.rotation,
            linear_velocity: b.velocity,
            angular_velocity: b.angular_velocity,
            is_sleeping: b.sleeping,
            is_active: b.active,
            kinetic_energy: b.kinetic_energy(),
        })
    }

    /// Return all active body handles as a `Vec`u32`.
    pub fn get_all_body_handles(&self) -> Vec<u32> {
        self.bodies
            .iter()
            .enumerate()
            .filter(|(_, b)| b.active)
            .map(|(i, _)| i as u32)
            .collect()
    }

    /// Return all body positions as a flat array `\[x0, y0, z0, x1, y1, z1, ...\]`.
    pub fn get_all_positions(&self) -> Vec<f64> {
        let mut result = Vec::with_capacity(self.bodies.len() * 3);
        for body in &self.bodies {
            if body.active {
                result.extend_from_slice(&body.position);
            }
        }
        result
    }

    /// Return all body transforms as a flat array of `\[f64; 7\]` values.
    pub fn get_all_transforms(&self) -> Vec<f64> {
        let mut result = Vec::with_capacity(self.bodies.len() * 7);
        for body in &self.bodies {
            if body.active {
                result.extend_from_slice(&body.position);
                result.extend_from_slice(&body.rotation);
            }
        }
        result
    }

    /// Return the contact results from the last `step()` call.
    pub fn get_contacts(&self) -> &[ContactResult] {
        &self.contacts
    }

    /// Return the number of contacts detected in the last `step()`.
    pub fn get_contact_count(&self) -> u32 {
        self.contacts.len() as u32
    }

    /// Return contacts involving a specific body.
    pub fn get_contacts_for_body(&self, handle: u32) -> Vec<ContactResult> {
        self.contacts
            .iter()
            .filter(|c| c.body_a == handle || c.body_b == handle)
            .cloned()
            .collect()
    }

    /// Return the accumulated simulation time in seconds.
    pub fn time(&self) -> f64 {
        self.time
    }

    /// Return the current gravity vector `\[gx, gy, gz\]`.
    pub fn gravity(&self) -> [f64; 3] {
        self.gravity
    }

    /// Return the number of **active** bodies (legacy alias for `get_body_count`).
    pub fn get_num_bodies(&self) -> u32 {
        self.get_body_count()
    }

    /// Return debug information from the last `step()` call.
    pub fn debug_info(&self) -> &DebugInfo {
        &self.debug_info
    }

    // -----------------------------------------------------------------------
    // Collider management
    // -----------------------------------------------------------------------

    /// Add a sphere collider to a body and return the collider handle.
    pub fn add_sphere_collider(&mut self, body: u32, radius: f64) -> u32 {
        let cfg = ColliderConfig::sphere(radius);
        self.add_collider(body, &cfg)
    }

    /// Add a box collider to a body and return the collider handle.
    pub fn add_box_collider(&mut self, body: u32, hx: f64, hy: f64, hz: f64) -> u32 {
        let cfg = ColliderConfig::cuboid(hx, hy, hz);
        self.add_collider(body, &cfg)
    }

    /// Add a capsule collider to a body and return the collider handle.
    pub fn add_capsule_collider(&mut self, body: u32, radius: f64, height: f64) -> u32 {
        let cfg = ColliderConfig::capsule(radius, height);
        self.add_collider(body, &cfg)
    }

    /// Add a static plane collider (infinite ground plane).
    pub fn add_plane_collider(&mut self, body: u32, nx: f64, ny: f64, nz: f64, offset: f64) -> u32 {
        let cfg = ColliderConfig::plane(nx, ny, nz, offset);
        self.add_collider(body, &cfg)
    }

    /// Add a collider using a full `ColliderConfig` and return the collider handle.
    pub fn add_collider(&mut self, body_index: u32, config: &ColliderConfig) -> u32 {
        if let Some(idx) = self.collider_free_list.pop() {
            let old_gen = self.colliders[idx].generation;
            self.colliders[idx] = ColliderEntry {
                body_index,
                config: config.clone(),
                active: true,
                generation: old_gen + 1,
            };
            idx as u32
        } else {
            let idx = self.colliders.len();
            self.colliders.push(ColliderEntry {
                body_index,
                config: config.clone(),
                active: true,
                generation: 0,
            });
            idx as u32
        }
    }

    /// Remove the collider with the given handle.
    pub fn remove_collider(&mut self, handle: u32) -> Result<()> {
        let idx = handle as usize;
        match self.colliders.get_mut(idx) {
            Some(col) if col.active => {
                col.active = false;
                self.collider_free_list.push(idx);
                Ok(())
            }
            Some(_) => Err(Error::InvalidHandle(handle)),
            None => Err(Error::InvalidHandle(handle)),
        }
    }

    // -----------------------------------------------------------------------
    // Configuration
    // -----------------------------------------------------------------------

    /// Change the gravity vector at runtime.
    pub fn set_gravity(&mut self, gx: f64, gy: f64, gz: f64) {
        self.gravity = [gx, gy, gz];
        self.config.gravity = [gx, gy, gz];
    }

    /// Return a reference to the current simulation config.
    pub fn config(&self) -> &SimulationConfig {
        &self.config
    }

    /// Set the fixed time step used by `step()`.
    pub fn set_fixed_dt(&mut self, dt: f64) {
        self.config.fixed_dt = dt;
    }

    /// Set the linear damping on a body (0 = no damping).
    pub fn set_body_linear_damping(
        &mut self,
        handle: u32,
        damping: f64,
    ) -> crate::error::Result<()> {
        self.get_body_mut(handle)
            .map(|b| b.linear_damping = damping)
    }

    /// Set the angular damping on a body (0 = no damping).
    pub fn set_body_angular_damping(
        &mut self,
        handle: u32,
        damping: f64,
    ) -> crate::error::Result<()> {
        self.get_body_mut(handle)
            .map(|b| b.angular_damping = damping)
    }

    // -----------------------------------------------------------------------
    // Raycasting
    // -----------------------------------------------------------------------

    /// Cast a ray from `origin` in `direction` and return the closest hit.
    #[allow(clippy::too_many_arguments)]
    pub fn raycast(
        &self,
        ox: f64,
        oy: f64,
        oz: f64,
        dx: f64,
        dy: f64,
        dz: f64,
        max_distance: f64,
    ) -> RaycastResult {
        let origin = Vec3Wasm::new(ox, oy, oz);
        let dir = Vec3Wasm::new(dx, dy, dz).normalized();

        let mut best = RaycastResult::no_hit();

        for col in &self.colliders {
            if !col.active {
                continue;
            }
            let body_idx = col.body_index as usize;
            let body = match self.bodies.get(body_idx) {
                Some(b) if b.active => b,
                _ => continue,
            };

            if col.config.shape_type == ColliderShapeType::Sphere {
                let center = Vec3Wasm::from_array(body.position);
                let oc = origin.sub(&center);
                let a = dir.dot(&dir);
                let b = 2.0 * oc.dot(&dir);
                let c = oc.dot(&oc) - col.config.radius * col.config.radius;
                let discriminant = b * b - 4.0 * a * c;

                if discriminant >= 0.0 {
                    let t = (-b - discriminant.sqrt()) / (2.0 * a);
                    if t > 1e-6 && t < max_distance && t < best.distance {
                        let hit_point = origin.add(&dir.scale(t));
                        let normal = hit_point.sub(&center).normalized();
                        best = RaycastResult {
                            hit: true,
                            body_handle: col.body_index,
                            point: hit_point.to_array(),
                            normal: normal.to_array(),
                            distance: t,
                        };
                    }
                }
            }
        }

        best
    }

    // -----------------------------------------------------------------------
    // Simulation step
    // -----------------------------------------------------------------------

    /// Step the simulation forward by `dt` seconds.
    pub fn step(&mut self, dt: f64) {
        if dt <= 0.0 {
            return;
        }

        let substep_dt = self.config.fixed_dt;
        let max_substeps = self.config.max_substeps;
        let mut remaining = dt;
        let mut substeps = 0u32;

        self.contacts.clear();

        while remaining > 1e-12 && substeps < max_substeps {
            let sub_dt = remaining.min(substep_dt);
            self.integrate(sub_dt);
            self.detect_and_resolve_contacts(sub_dt);
            remaining -= sub_dt;
            substeps += 1;
        }

        self.time += dt;
        self.debug_info.solver_iterations_performed = substeps;
    }

    // -----------------------------------------------------------------------
    // Reset
    // -----------------------------------------------------------------------

    /// Remove all bodies, colliders and contacts, resetting to an empty state.
    pub fn reset(&mut self) {
        self.bodies.clear();
        self.body_free_list.clear();
        self.colliders.clear();
        self.collider_free_list.clear();
        self.contacts.clear();
        self.time = 0.0;
        self.debug_info = DebugInfo::default();
    }

    /// Reset and change the gravity at the same time.
    pub fn reset_with_gravity(&mut self, gx: f64, gy: f64, gz: f64) {
        self.reset();
        self.set_gravity(gx, gy, gz);
    }

    // -----------------------------------------------------------------------
    // Private helpers
    // -----------------------------------------------------------------------

    fn get_body(&self, handle: u32) -> Option<&BodyEntry> {
        self.bodies.get(handle as usize).filter(|b| b.active)
    }

    /// Get inverse mass for a body handle (0 for static/inactive bodies).
    pub(crate) fn body_inv_mass(&self, handle: u32) -> f64 {
        self.get_body(handle).map(|b| b.inv_mass).unwrap_or(0.0)
    }

    /// Whether a body is dynamic (movable).
    pub(crate) fn body_is_dynamic(&self, handle: u32) -> bool {
        self.get_body(handle)
            .map(|b| b.is_dynamic())
            .unwrap_or(false)
    }

    fn get_body_mut(&mut self, handle: u32) -> Result<&mut BodyEntry> {
        let idx = handle as usize;
        self.bodies
            .get_mut(idx)
            .filter(|b| b.active)
            .ok_or(Error::InvalidHandle(handle))
    }

    /// Integrate velocities and positions (symplectic Euler).
    fn integrate(&mut self, dt: f64) {
        let gravity = self.gravity;
        let lin_sleep = self.config.linear_sleep_threshold;
        let ang_sleep = self.config.angular_sleep_threshold;
        let time_sleep = self.config.time_before_sleep;
        let sleeping_enabled = self.config.sleeping_enabled;

        for body in &mut self.bodies {
            if !body.active || !body.is_dynamic() {
                continue;
            }
            if body.sleeping {
                let f2 = body.force[0] * body.force[0]
                    + body.force[1] * body.force[1]
                    + body.force[2] * body.force[2];
                if f2 > 1e-20 {
                    body.sleeping = false;
                    body.sleep_timer = 0.0;
                } else {
                    body.force = [0.0; 3];
                    body.torque = [0.0; 3];
                    continue;
                }
            }

            // 1. Apply gravity scaled by gravity_scale
            let gs = body.gravity_scale;
            let ax = gravity[0] * gs + body.force[0] * body.inv_mass;
            let ay = gravity[1] * gs + body.force[1] * body.inv_mass;
            let az = gravity[2] * gs + body.force[2] * body.inv_mass;

            // 2. Update velocity (symplectic Euler: update v first)
            body.velocity[0] += ax * dt;
            body.velocity[1] += ay * dt;
            body.velocity[2] += az * dt;

            // 3. Apply damping
            let ld = (1.0 - body.linear_damping * dt).max(0.0);
            let ad = (1.0 - body.angular_damping * dt).max(0.0);
            body.velocity[0] *= ld;
            body.velocity[1] *= ld;
            body.velocity[2] *= ld;
            body.angular_velocity[0] *= ad;
            body.angular_velocity[1] *= ad;
            body.angular_velocity[2] *= ad;

            // 4. Update position
            body.position[0] += body.velocity[0] * dt;
            body.position[1] += body.velocity[1] * dt;
            body.position[2] += body.velocity[2] * dt;

            // 5. Integrate angular velocity -> update orientation (small angle approx)
            let wx = body.angular_velocity[0];
            let wy = body.angular_velocity[1];
            let wz = body.angular_velocity[2];
            let half_dt = 0.5 * dt;
            let q = body.rotation;
            let dqx = half_dt * (q[3] * wx - q[2] * wy + q[1] * wz);
            let dqy = half_dt * (q[3] * wy + q[2] * wx - q[0] * wz);
            let dqz = half_dt * (q[3] * wz - q[1] * wx + q[0] * wy);
            let dqw = half_dt * (-q[0] * wx - q[1] * wy - q[2] * wz);
            let nx = q[0] + dqx;
            let ny = q[1] + dqy;
            let nz = q[2] + dqz;
            let nw = q[3] + dqw;
            let qnorm = (nx * nx + ny * ny + nz * nz + nw * nw).sqrt().max(1e-15);
            body.rotation = [nx / qnorm, ny / qnorm, nz / qnorm, nw / qnorm];

            // 6. Update sleep timer
            if sleeping_enabled {
                let v2 = body.velocity[0] * body.velocity[0]
                    + body.velocity[1] * body.velocity[1]
                    + body.velocity[2] * body.velocity[2];
                let w2 = body.angular_velocity[0] * body.angular_velocity[0]
                    + body.angular_velocity[1] * body.angular_velocity[1]
                    + body.angular_velocity[2] * body.angular_velocity[2];
                if v2 < lin_sleep * lin_sleep && w2 < ang_sleep * ang_sleep {
                    body.sleep_timer += dt;
                    if body.sleep_timer > time_sleep {
                        body.sleeping = true;
                    }
                } else {
                    body.sleep_timer = 0.0;
                }
            }

            // 7. Clear accumulated forces
            body.force = [0.0; 3];
            body.torque = [0.0; 3];
        }
    }

    /// Detect sphere-sphere and sphere-plane contacts and resolve with impulses.
    fn detect_and_resolve_contacts(&mut self, dt: f64) {
        let manifolds = self.build_manifolds();
        self.debug_info.narrowphase_tests = manifolds.len() as u32;

        for manifold in &manifolds {
            let contact = self.resolve_manifold(manifold, dt);
            self.contacts.push(contact);
        }

        self.debug_info.active_contacts = self.contacts.len() as u32;
    }

    /// Build contact manifolds by testing all collider pairs.
    fn build_manifolds(&self) -> Vec<ManifoldPair> {
        let mut manifolds = Vec::new();

        let n = self.colliders.len();
        for i in 0..n {
            if !self.colliders[i].active {
                continue;
            }
            let ci = &self.colliders[i];
            let bi_idx = ci.body_index as usize;
            let body_i = match self.bodies.get(bi_idx) {
                Some(b) if b.active => b,
                _ => continue,
            };
            let pos_i = body_i.position;

            for j in (i + 1)..n {
                if !self.colliders[j].active {
                    continue;
                }
                let cj = &self.colliders[j];
                if ci.body_index == cj.body_index {
                    continue;
                }
                let bj_idx = cj.body_index as usize;
                let body_j = match self.bodies.get(bj_idx) {
                    Some(b) if b.active => b,
                    _ => continue,
                };
                let pos_j = body_j.position;

                match (&ci.config.shape_type, &cj.config.shape_type) {
                    (ColliderShapeType::Sphere, ColliderShapeType::Sphere) => {
                        let dx = pos_j[0] - pos_i[0];
                        let dy = pos_j[1] - pos_i[1];
                        let dz = pos_j[2] - pos_i[2];
                        let dist2 = dx * dx + dy * dy + dz * dz;
                        let sum_r = ci.config.radius + cj.config.radius;
                        if dist2 < sum_r * sum_r {
                            let dist = dist2.sqrt().max(1e-12);
                            let depth = sum_r - dist;
                            let nx = dx / dist;
                            let ny = dy / dist;
                            let nz = dz / dist;
                            let pa = [
                                pos_i[0] + nx * ci.config.radius,
                                pos_i[1] + ny * ci.config.radius,
                                pos_i[2] + nz * ci.config.radius,
                            ];
                            let pb = [
                                pos_j[0] - nx * cj.config.radius,
                                pos_j[1] - ny * cj.config.radius,
                                pos_j[2] - nz * cj.config.radius,
                            ];
                            manifolds.push(ManifoldPair {
                                body_a: ci.body_index,
                                body_b: cj.body_index,
                                normal: [nx, ny, nz],
                                depth,
                                point_a: pa,
                                point_b: pb,
                            });
                        }
                    }
                    (ColliderShapeType::Sphere, ColliderShapeType::Plane)
                    | (ColliderShapeType::Plane, ColliderShapeType::Sphere) => {
                        let (sphere_body_pos, sphere_r, plane_cfg, sphere_body_idx, plane_body_idx) =
                            if ci.config.shape_type == ColliderShapeType::Sphere {
                                (
                                    pos_i,
                                    ci.config.radius,
                                    &cj.config,
                                    ci.body_index,
                                    cj.body_index,
                                )
                            } else {
                                (
                                    pos_j,
                                    cj.config.radius,
                                    &ci.config,
                                    cj.body_index,
                                    ci.body_index,
                                )
                            };
                        let pn = plane_cfg.plane_normal;
                        let offset = plane_cfg.plane_offset;
                        let signed_dist = pn[0] * sphere_body_pos[0]
                            + pn[1] * sphere_body_pos[1]
                            + pn[2] * sphere_body_pos[2]
                            - offset;
                        let depth = sphere_r - signed_dist;
                        if depth > 0.0 {
                            let pa = [
                                sphere_body_pos[0] - pn[0] * sphere_r,
                                sphere_body_pos[1] - pn[1] * sphere_r,
                                sphere_body_pos[2] - pn[2] * sphere_r,
                            ];
                            let pb = [
                                sphere_body_pos[0] - pn[0] * signed_dist,
                                sphere_body_pos[1] - pn[1] * signed_dist,
                                sphere_body_pos[2] - pn[2] * signed_dist,
                            ];
                            manifolds.push(ManifoldPair {
                                body_a: sphere_body_idx,
                                body_b: plane_body_idx,
                                normal: pn,
                                depth,
                                point_a: pa,
                                point_b: pb,
                            });
                        }
                    }
                    _ => {}
                }
            }
        }

        manifolds
    }

    /// Resolve a single contact manifold using impulse-based collision response.
    fn resolve_manifold(&mut self, m: &ManifoldPair, _dt: f64) -> ContactResult {
        let mut cr = ContactResult {
            body_a: m.body_a,
            body_b: m.body_b,
            point_on_a: m.point_a,
            point_on_b: m.point_b,
            normal: m.normal,
            depth: m.depth,
            relative_velocity: 0.0,
            impulse: 0.0,
            is_new: true,
            friction_impulse: 0.0,
        };

        let ia = m.body_a as usize;
        let ib = m.body_b as usize;

        let inv_mass_a = self
            .bodies
            .get(ia)
            .map(|b| if b.is_dynamic() { b.inv_mass } else { 0.0 })
            .unwrap_or(0.0);
        let inv_mass_b = self
            .bodies
            .get(ib)
            .map(|b| if b.is_dynamic() { b.inv_mass } else { 0.0 })
            .unwrap_or(0.0);
        let total_inv = inv_mass_a + inv_mass_b;
        if total_inv == 0.0 {
            return cr;
        }

        let vel_a = self.bodies.get(ia).map(|b| b.velocity).unwrap_or([0.0; 3]);
        let vel_b = self.bodies.get(ib).map(|b| b.velocity).unwrap_or([0.0; 3]);
        let rest_a = self.bodies.get(ia).map(|b| b.restitution).unwrap_or(0.0);
        let rest_b = self.bodies.get(ib).map(|b| b.restitution).unwrap_or(0.0);
        let fric_a = self.bodies.get(ia).map(|b| b.friction).unwrap_or(0.5);
        let fric_b = self.bodies.get(ib).map(|b| b.friction).unwrap_or(0.5);

        let nx = m.normal[0];
        let ny = m.normal[1];
        let nz = m.normal[2];

        let rel_vx = vel_a[0] - vel_b[0];
        let rel_vy = vel_a[1] - vel_b[1];
        let rel_vz = vel_a[2] - vel_b[2];
        let rel_vn = rel_vx * nx + rel_vy * ny + rel_vz * nz;

        cr.relative_velocity = rel_vn;

        if rel_vn > 0.0 {
            return cr;
        }

        let e = (rest_a * rest_b).sqrt();
        let j = -(1.0 + e) * rel_vn / total_inv;
        cr.impulse = j;

        if let Some(body) = self.bodies.get_mut(ia)
            && body.is_dynamic()
        {
            body.velocity[0] += j * inv_mass_a * nx;
            body.velocity[1] += j * inv_mass_a * ny;
            body.velocity[2] += j * inv_mass_a * nz;
        }
        if let Some(body) = self.bodies.get_mut(ib)
            && body.is_dynamic()
        {
            body.velocity[0] -= j * inv_mass_b * nx;
            body.velocity[1] -= j * inv_mass_b * ny;
            body.velocity[2] -= j * inv_mass_b * nz;
        }

        // Friction impulse (Coulomb model)
        let mu = (fric_a * fric_b).sqrt();
        let vel_a2 = self.bodies.get(ia).map(|b| b.velocity).unwrap_or([0.0; 3]);
        let vel_b2 = self.bodies.get(ib).map(|b| b.velocity).unwrap_or([0.0; 3]);
        let rel_tx = vel_a2[0] - vel_b2[0];
        let rel_ty = vel_a2[1] - vel_b2[1];
        let rel_tz = vel_a2[2] - vel_b2[2];
        let rel_tn = rel_tx * nx + rel_ty * ny + rel_tz * nz;
        let tx = rel_tx - rel_tn * nx;
        let ty = rel_ty - rel_tn * ny;
        let tz = rel_tz - rel_tn * nz;
        let t_len = (tx * tx + ty * ty + tz * tz).sqrt();
        if t_len > 1e-10 {
            let tx = tx / t_len;
            let ty = ty / t_len;
            let tz = tz / t_len;
            let jt_raw = -(rel_tx * tx + rel_ty * ty + rel_tz * tz) / total_inv;
            let jt = jt_raw.clamp(-mu * j.abs(), mu * j.abs());
            cr.friction_impulse = jt;
            if let Some(body) = self.bodies.get_mut(ia)
                && body.is_dynamic()
            {
                body.velocity[0] += jt * inv_mass_a * tx;
                body.velocity[1] += jt * inv_mass_a * ty;
                body.velocity[2] += jt * inv_mass_a * tz;
            }
            if let Some(body) = self.bodies.get_mut(ib)
                && body.is_dynamic()
            {
                body.velocity[0] -= jt * inv_mass_b * tx;
                body.velocity[1] -= jt * inv_mass_b * ty;
                body.velocity[2] -= jt * inv_mass_b * tz;
            }
        }

        // Positional correction (Baumgarte)
        let slop = self.config.contact_slop;
        let corr_mag = (m.depth - slop).max(0.0) * 0.2 / total_inv;
        if corr_mag > 1e-10 {
            if let Some(body) = self.bodies.get_mut(ia)
                && body.is_dynamic()
            {
                body.position[0] += corr_mag * inv_mass_a * nx;
                body.position[1] += corr_mag * inv_mass_a * ny;
                body.position[2] += corr_mag * inv_mass_a * nz;
            }
            if let Some(body) = self.bodies.get_mut(ib)
                && body.is_dynamic()
            {
                body.position[0] -= corr_mag * inv_mass_b * nx;
                body.position[1] -= corr_mag * inv_mass_b * ny;
                body.position[2] -= corr_mag * inv_mass_b * nz;
            }
        }

        cr
    }
}
