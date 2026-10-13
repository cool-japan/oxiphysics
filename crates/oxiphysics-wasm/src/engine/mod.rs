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
//! - A `Vec<BodyEntry>` of body states (position, velocity, forces, inertia, ...)
//! - A `Vec<ColliderEntry>` of collider shapes attached to bodies
//! - A `Vec<ContactResult>` list refreshed each step
//! - Free lists of recycled body and collider slots
//!
//! ## Handles
//!
//! Bodies and colliders are addressed by opaque `u32` handles.  A handle packs
//! the slot index (low 20 bits) together with a generation counter (high 12
//! bits) that is bumped every time the slot is recycled, so a handle that was
//! returned for a body that has since been removed can never be used to
//! control whatever body reuses its slot: the stale handle is rejected with
//! [`Error::InvalidHandle`].  On an engine that never removes anything the
//! `n`-th body simply has handle `n`.  At most [`MAX_SLOTS`] bodies and
//! [`MAX_SLOTS`] colliders can exist at once ([`Error::CapacityExceeded`]).
//!
//! ## Time stepping
//!
//! [`WasmPhysicsEngine::step`] adds the elapsed time to a fixed-step
//! accumulator and then runs as many substeps of exactly
//! [`SimulationConfig::fixed_dt`] as fit (at most
//! [`SimulationConfig::max_substeps`]); [`WasmPhysicsEngine::time`] advances by
//! the time that was really integrated.  The left-over fraction of a step is
//! exposed through [`WasmPhysicsEngine::interpolation_alpha`] for render
//! interpolation.  If the cap is hit the surplus time is dropped (and reported
//! in [`DebugInfo::dropped_time`](crate::types::DebugInfo::dropped_time))
//! instead of letting the simulation fall further and further behind.
//!
//! ## Dynamics
//!
//! Each substep integrates with symplectic Euler.  Rotation follows Euler's
//! equations for a rigid body: external torques are mapped through the
//! world-space inverse inertia tensor (computed from the attached colliders)
//! and a gyroscopic term is applied implicitly, so non-spherical bodies
//! precess correctly and remain stable at large spin rates.  Contacts are
//! resolved with linear impulses (normal, Coulomb friction) plus positional
//! correction; contact impulses do not exert torque.
//!
//! ## Collision shapes
//!
//! Spheres, oriented boxes, capsules, cylinders and planes collide with each
//! other (the `narrowphase` module documents the algorithm used per pair:
//! closed-form tests for sphere, capsule and plane pairs, SAT with face
//! clipping for boxes, GJK/EPA for the remaining convex pairs).  Ray queries
//! support the same shapes.  Convex-hull colliders carry no vertex data in
//! [`ColliderConfig`] and are rejected.
//!
//! ## Example
//!
//! ```no_run
//! use oxiphysics_wasm::WasmPhysicsEngine;
//!
//! # fn main() -> Result<(), oxiphysics_wasm::Error> {
//! let mut engine = WasmPhysicsEngine::new(0.0, -9.81, 0.0);
//! let body = engine.add_dynamic_body(1.0, 0.0, 10.0, 0.0)?;
//! engine.step(1.0 / 60.0)?;
//! let pos = engine.get_position(body);
//! assert!(pos[1] < 10.0); // body fell due to gravity
//! # Ok(())
//! # }
//! ```

// Sub-modules
mod bindings;
mod constraint;
mod float64_view;
mod geom;
mod gjk;
mod handles;
mod inertia;
mod js_exports;
mod lbm;
mod narrowphase;
mod raycast;
mod shapes;
mod sph;
mod validation;
mod vehicle;
mod wasm_engine;
mod wasm_types;

#[cfg(test)]
mod tests;
#[cfg(test)]
mod tests_collision;
#[cfg(test)]
mod tests_dynamics;
#[cfg(test)]
mod tests_handles;
#[cfg(test)]
mod tests_raycast;
#[cfg(test)]
mod tests_validation;

// Re-export everything from submodules
pub use bindings::{ContactInfoEntry, WasmBindings, WasmContactList};
pub use constraint::{ContactData, EngineConstraintSolver, WasmConstraint};
pub use float64_view::Float64View;
pub use handles::MAX_SLOTS;
pub use lbm::{EngineWasmLbmConfig, WasmLbmSim};
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

use wasm_bindgen::prelude::*;

use self::geom::{
    MAT3_ZERO, Mat3, Vec3, add, cross, dot, length, length_sq, madd, mat_add, mat_inverse, mat_mul,
    mat_scale, mat_skew, mat_t_vec, mat_transpose, mat_vec, normalized, quat_from_axis_angle,
    quat_mul, quat_normalized, quat_to_mat3, scale, sub,
};
use self::handles::{decode_handle, encode_handle, next_generation};
use self::inertia::{ShapePart, compound_inertia};
use self::narrowphase::{bounds_overlap, collide};
use self::raycast::{point_at, ray_shape};
use self::shapes::{Shape, WorldShape};
use self::validation::{
    non_negative_param, require_finite3, validate_body_config, validate_collider_config,
    validate_dynamic_mass,
};
use crate::error::{Error, Result};
use crate::types::{
    BodyState, ColliderConfig, ContactResult, DebugInfo, RaycastResult, RigidBodyConfig,
    SimulationConfig,
};

/// Relative slack of the fixed-step accumulator test.
///
/// Time that sums to one ulp short of `fixed_dt` (for example `0.3 - 0.1 - 0.1`
/// against `0.1`) still triggers the substep instead of being deferred to the
/// next call.
const STEP_EPSILON: f64 = 1e-9;

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
    /// Current angular velocity (world frame, rad/s).
    angular_velocity: [f64; 3],
    /// Force accumulated for the next `step()` (world frame).
    force: [f64; 3],
    /// Torque accumulated for the next `step()` (world frame).
    torque: [f64; 3],
    /// Mass (0 for static).
    mass: f64,
    /// Inverse mass (0 for static).
    inv_mass: f64,
    /// Rotational inertia tensor in the body frame (zero until the body has a
    /// collider that encloses a volume).
    inertia_local: Mat3,
    /// Inverse of `inertia_local` (zero when the body has no inertia, in
    /// which case torques and angular impulses have no effect).
    inv_inertia_local: Mat3,
    /// Whether `inertia_local` / `inv_inertia_local` are valid.
    has_inertia: bool,
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
    /// Whether the body is currently sleeping.
    sleeping: bool,
    /// Time the body has been below sleep thresholds (seconds).
    sleep_timer: f64,
    /// Generation of the slot; part of every handle that refers to this body.
    generation: u32,
}

impl BodyEntry {
    /// Build a body from an already validated configuration.
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
            inertia_local: MAT3_ZERO,
            inv_inertia_local: MAT3_ZERO,
            has_inertia: false,
            linear_damping: cfg.linear_damping,
            angular_damping: cfg.angular_damping,
            restitution: cfg.restitution,
            friction: cfg.friction,
            gravity_scale: cfg.gravity_scale,
            active: true,
            is_static: cfg.is_static(),
            is_kinematic: cfg.is_kinematic,
            sleeping: false,
            sleep_timer: 0.0,
            generation,
        }
    }

    /// Returns `true` if the body is movable (dynamic, not static or kinematic).
    fn is_dynamic(&self) -> bool {
        !self.is_static && !self.is_kinematic && self.mass > 0.0
    }

    /// Kinetic energy: the translational part `m v^2 / 2` plus the rotational
    /// part `w . (I w) / 2` (evaluated in the body frame; zero for a body
    /// without rotational inertia).
    fn kinetic_energy(&self) -> f64 {
        if self.mass <= 0.0 {
            return 0.0;
        }
        let translational = 0.5 * self.mass * length_sq(self.velocity);
        let rotational = if self.has_inertia {
            let w_body = mat_t_vec(&quat_to_mat3(self.rotation), self.angular_velocity);
            0.5 * dot(w_body, mat_vec(&self.inertia_local, w_body))
        } else {
            0.0
        };
        translational + rotational
    }

    /// Inverse inertia tensor in world coordinates, `R * I_body^-1 * R^T`.
    fn world_inv_inertia(&self) -> Mat3 {
        let r = quat_to_mat3(self.rotation);
        mat_mul(&mat_mul(&r, &self.inv_inertia_local), &mat_transpose(&r))
    }

    /// Change the angular velocity over `dt` by the applied torque and the
    /// gyroscopic term.
    ///
    /// The angular velocity is expressed in the world frame, so Euler's
    /// equations `I w' + w x (I w) = tau` read, in the body frame,
    /// `w_b' = I_b^-1 (R^T tau - w_b x (I_b w_b))`.  The torque is applied
    /// through the world-space inverse inertia; the gyroscopic term is solved
    /// implicitly with a single Newton step (unconditionally stable, exact
    /// for isotropic inertia where the term vanishes).
    fn accelerate_rotation(&mut self, dt: f64) {
        if self.has_inertia {
            let r = quat_to_mat3(self.rotation);
            let inv_world = mat_mul(&mat_mul(&r, &self.inv_inertia_local), &mat_transpose(&r));
            self.angular_velocity =
                madd(self.angular_velocity, mat_vec(&inv_world, self.torque), dt);

            let w_body = mat_t_vec(&r, self.angular_velocity);
            let i_body = &self.inertia_local;
            let l_body = mat_vec(i_body, w_body);
            // Residual of the implicit gyroscopic equation at the current w,
            // and its Jacobian  I + dt ([w]x I - [I w]x).
            let residual = scale(cross(w_body, l_body), dt);
            let jacobian = mat_add(
                i_body,
                &mat_scale(
                    &mat_add(
                        &mat_mul(&mat_skew(w_body), i_body),
                        &mat_scale(&mat_skew(l_body), -1.0),
                    ),
                    dt,
                ),
            );
            if let Some(jacobian_inv) = mat_inverse(&jacobian) {
                let delta_body = mat_vec(&jacobian_inv, residual);
                self.angular_velocity = sub(self.angular_velocity, mat_vec(&r, delta_body));
            }
        }
    }

    /// Advance the orientation by the world-frame angular velocity.
    ///
    /// The exact rotation `exp(dt * w)` is composed on the left of the current
    /// orientation (`q' = dq * q`), the convention for a world-space `w`.
    fn integrate_orientation(&mut self, dt: f64) {
        let speed = length(self.angular_velocity);
        let angle = speed * dt;
        if !angle.is_finite() || angle <= 0.0 {
            return;
        }
        let dq = if angle < 1e-9 {
            // First-order expansion keeps the tiny-angle case free of 0/0.
            let half_dt = 0.5 * dt;
            [
                self.angular_velocity[0] * half_dt,
                self.angular_velocity[1] * half_dt,
                self.angular_velocity[2] * half_dt,
                1.0,
            ]
        } else {
            quat_from_axis_angle(scale(self.angular_velocity, 1.0 / speed), angle)
        };
        if let Some(q) = quat_normalized(quat_mul(dq, self.rotation)) {
            self.rotation = q;
        }
    }
}

// ---------------------------------------------------------------------------
// Internal collider storage
// ---------------------------------------------------------------------------

/// Internal collider state stored in the engine.
#[derive(Debug, Clone)]
struct ColliderEntry {
    /// Handle of the owning body (generation included).
    body_handle: u32,
    /// Validated shape configuration (unit plane normal and quaternion).
    config: ColliderConfig,
    /// Local geometry derived from `config`.
    shape: Shape,
    /// Rotation of the collider relative to its body.
    local_rotation: Mat3,
    /// Whether this slot is occupied.
    active: bool,
    /// Generation of the slot; part of every handle that refers to it.
    generation: u32,
}

// ---------------------------------------------------------------------------
// Contact manifold
// ---------------------------------------------------------------------------

/// A contact found by the narrow phase and waiting to be resolved.
///
/// `body_a` / `body_b` are body handles; `normal` is the unit normal pointing
/// from body B toward body A (the convention of
/// [`ContactResult`](crate::types::ContactResult)).
#[derive(Debug, Clone)]
struct ManifoldPair {
    body_a: u32,
    body_b: u32,
    normal: [f64; 3],
    depth: f64,
    point_a: [f64; 3],
    point_b: [f64; 3],
}

/// Counters of one collision pass (one fixed substep).
#[derive(Debug, Clone, Copy, Default)]
struct PassStats {
    /// Collider pairs on different bodies that were examined.
    candidate_pairs: u32,
    /// Pairs whose bounding spheres overlap and that went through the narrow
    /// phase.
    narrowphase_tests: u32,
    /// Contacts found and resolved.
    contacts: u32,
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
/// 1. Create with `new WasmPhysicsEngine(gx, gy, gz)` (JavaScript) or
///    [`WasmPhysicsEngine::new`] (Rust).
/// 2. Add bodies with `add_dynamic_body`, `add_static_body`.
/// 3. Add colliders with `add_sphere_collider`, `add_box_collider`, etc.
/// 4. Advance simulation with `step(dt)`.
/// 5. Query state with `get_position(handle)`, `get_velocity(handle)`, etc.
/// 6. Remove bodies with `remove_body(handle)` or reset with `reset()`.
///
/// Every operation that creates something (`add_*`) validates its input and
/// reports problems as an [`Error`] (an exception in JavaScript) instead of
/// storing nonsensical state.
#[wasm_bindgen]
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
    /// Elapsed time handed to `step()` that has not been consumed by a fixed
    /// substep yet (always below `fixed_dt` between calls).
    accumulator: f64,
    /// Maximum number of slots in each of the body and collider tables
    /// ([`MAX_SLOTS`]; lowered only by tests).
    slot_limit: usize,
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
        Self::with_valid_config(config)
    }

    /// Create an engine from a `SimulationConfig`.
    ///
    /// # Errors
    ///
    /// Returns the error of [`SimulationConfig::validate`] when the
    /// configuration is unusable (non-positive or non-finite `fixed_dt`, zero
    /// `max_substeps`, non-finite gravity, negative thresholds, ...).
    pub fn from_config(config: SimulationConfig) -> Result<Self> {
        config.validate()?;
        Ok(Self::with_valid_config(config))
    }

    /// Build an engine around a configuration that is known to be valid.
    fn with_valid_config(config: SimulationConfig) -> Self {
        Self {
            gravity: config.gravity,
            bodies: Vec::new(),
            body_free_list: Vec::new(),
            colliders: Vec::new(),
            collider_free_list: Vec::new(),
            contacts: Vec::new(),
            time: 0.0,
            accumulator: 0.0,
            slot_limit: MAX_SLOTS,
            config,
            debug_info: DebugInfo::default(),
        }
    }

    // -----------------------------------------------------------------------
    // Body management
    // -----------------------------------------------------------------------

    /// Add a rigid body using a [`RigidBodyConfig`] and return its handle.
    ///
    /// A mass of exactly `0` creates a static body (or a kinematic one when
    /// `is_kinematic` is set); any positive mass creates a dynamic body whose
    /// rotational inertia follows from the colliders attached later.
    ///
    /// # Errors
    ///
    /// * [`Error::InvalidMass`] - the mass is negative or not finite;
    /// * [`Error::InvalidRestitution`] / [`Error::InvalidFriction`] - a
    ///   material coefficient is out of range;
    /// * [`Error::InvalidParameter`] - a position, velocity, damping value or
    ///   the orientation is not finite (or the quaternion has zero length);
    /// * [`Error::CapacityExceeded`] - all [`MAX_SLOTS`] body slots are in use.
    pub fn add_rigid_body(&mut self, config: &RigidBodyConfig) -> Result<u32> {
        let config = validate_body_config(config)?;
        if let Some(idx) = self.body_free_list.pop() {
            let slot = self
                .bodies
                .get_mut(idx)
                .ok_or_else(|| Error::internal("body free list references a missing slot"))?;
            let generation = next_generation(slot.generation);
            *slot = BodyEntry::from_config(&config, generation);
            return Ok(encode_handle(idx, generation));
        }
        if self.bodies.len() >= self.slot_limit {
            return Err(Error::CapacityExceeded {
                max: self.slot_limit,
            });
        }
        let idx = self.bodies.len();
        self.bodies.push(BodyEntry::from_config(&config, 0));
        Ok(encode_handle(idx, 0))
    }

    /// Add a dynamic (movable) body and return its handle.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidMass`] unless `mass` is finite and strictly positive
    /// (a zero mass would silently create a static body - use
    /// [`add_static_body`](Self::add_static_body) for that), plus the errors
    /// of [`add_rigid_body`](Self::add_rigid_body).
    pub fn add_dynamic_body(&mut self, mass: f64, x: f64, y: f64, z: f64) -> Result<u32> {
        validate_dynamic_mass(mass)?;
        let cfg = RigidBodyConfig::new()
            .with_mass(mass)
            .with_position(x, y, z);
        self.add_rigid_body(&cfg)
    }

    /// Add a static (immovable) body and return its handle.
    ///
    /// # Errors
    ///
    /// The errors of [`add_rigid_body`](Self::add_rigid_body).
    pub fn add_static_body(&mut self, x: f64, y: f64, z: f64) -> Result<u32> {
        let cfg = RigidBodyConfig::static_body().with_position(x, y, z);
        self.add_rigid_body(&cfg)
    }

    /// Remove the body with the given handle together with all its colliders.
    ///
    /// The slot (and the slots of its colliders) is recycled with a new
    /// generation, so every handle that referred to the removed body or its
    /// colliders becomes invalid.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidHandle`] if the handle does not refer to a live body
    /// (never allocated, already removed, or stale).
    pub fn remove_body(&mut self, handle: u32) -> Result<()> {
        let idx = self.body_slot(handle).ok_or(Error::InvalidHandle(handle))?;
        if let Some(body) = self.bodies.get_mut(idx) {
            body.active = false;
        }
        for (collider_idx, collider) in self.colliders.iter_mut().enumerate() {
            if collider.active && collider.body_handle == handle {
                collider.active = false;
                self.collider_free_list.push(collider_idx);
            }
        }
        self.body_free_list.push(idx);
        Ok(())
    }

    /// `true` if `handle` refers to a live body of this engine.
    pub fn contains_body(&self, handle: u32) -> bool {
        self.body_slot(handle).is_some()
    }

    /// Slot index of the body a handle refers to (`None` for dead or stale
    /// handles).
    fn body_slot(&self, handle: u32) -> Option<usize> {
        let (idx, generation) = decode_handle(handle);
        self.bodies
            .get(idx)
            .filter(|b| b.active && b.generation == generation)
            .map(|_| idx)
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
    ///
    /// # Errors
    ///
    /// [`Error::InvalidHandle`] for a dead or stale handle,
    /// [`Error::InvalidParameter`] if a component is not finite (the body is
    /// left untouched).
    pub fn set_velocity(&mut self, handle: u32, vx: f64, vy: f64, vz: f64) -> Result<()> {
        let body = self.get_body_mut(handle)?;
        require_finite3("velocity", [vx, vy, vz])?;
        body.velocity = [vx, vy, vz];
        body.sleeping = false;
        body.sleep_timer = 0.0;
        Ok(())
    }

    /// Set the angular velocity of a body (radians/second).
    ///
    /// # Errors
    ///
    /// As [`set_velocity`](Self::set_velocity).
    pub fn set_angular_velocity(&mut self, handle: u32, wx: f64, wy: f64, wz: f64) -> Result<()> {
        let body = self.get_body_mut(handle)?;
        require_finite3("angular_velocity", [wx, wy, wz])?;
        body.angular_velocity = [wx, wy, wz];
        body.sleeping = false;
        body.sleep_timer = 0.0;
        Ok(())
    }

    /// Set the position of a body directly (teleport).
    ///
    /// # Errors
    ///
    /// As [`set_velocity`](Self::set_velocity).
    pub fn set_position(&mut self, handle: u32, x: f64, y: f64, z: f64) -> Result<()> {
        let body = self.get_body_mut(handle)?;
        require_finite3("position", [x, y, z])?;
        body.position = [x, y, z];
        Ok(())
    }

    /// Apply a world-space force to a body.
    ///
    /// Forces accumulate and act during **every** fixed substep of the next
    /// [`step`](Self::step) call, after which they are cleared (so the force
    /// has to be re-applied each frame).
    ///
    /// # Errors
    ///
    /// As [`set_velocity`](Self::set_velocity).
    pub fn apply_force(&mut self, handle: u32, fx: f64, fy: f64, fz: f64) -> Result<()> {
        let body = self.get_body_mut(handle)?;
        require_finite3("force", [fx, fy, fz])?;
        body.force = add(body.force, [fx, fy, fz]);
        body.sleeping = false;
        body.sleep_timer = 0.0;
        Ok(())
    }

    /// Apply a world-space torque to a body.
    ///
    /// Like [`apply_force`](Self::apply_force), the torque acts during every
    /// fixed substep of the next [`step`](Self::step) call and is then
    /// cleared.  It accelerates the body through its world-space inverse
    /// inertia tensor, which is derived from the shapes of the attached
    /// colliders; a body without a collider that encloses a volume (for
    /// example one that only carries planes) has infinite rotational inertia
    /// and is not spun up by torques.
    ///
    /// # Errors
    ///
    /// As [`set_velocity`](Self::set_velocity).
    pub fn apply_torque(&mut self, handle: u32, tx: f64, ty: f64, tz: f64) -> Result<()> {
        let body = self.get_body_mut(handle)?;
        require_finite3("torque", [tx, ty, tz])?;
        body.torque = add(body.torque, [tx, ty, tz]);
        body.sleeping = false;
        body.sleep_timer = 0.0;
        Ok(())
    }

    /// Apply an instantaneous linear impulse to a body.
    ///
    /// # Errors
    ///
    /// As [`set_velocity`](Self::set_velocity).
    pub fn apply_impulse(&mut self, handle: u32, ix: f64, iy: f64, iz: f64) -> Result<()> {
        let body = self.get_body_mut(handle)?;
        require_finite3("impulse", [ix, iy, iz])?;
        if body.inv_mass > 0.0 {
            body.velocity = madd(body.velocity, [ix, iy, iz], body.inv_mass);
            body.sleeping = false;
            body.sleep_timer = 0.0;
        }
        Ok(())
    }

    /// Apply an impulse `(ix, iy, iz)` at the world-space offset `(rx, ry, rz)`
    /// from the body centre.
    ///
    /// The linear velocity changes by `J / m` and the angular velocity by
    /// `I^-1 (r x J)` with the world-space inverse inertia tensor `I^-1`, so
    /// an off-centre impulse spins the body according to its mass
    /// distribution.
    ///
    /// # Errors
    ///
    /// As [`set_velocity`](Self::set_velocity).
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
        let body = self.get_body_mut(handle)?;
        require_finite3("impulse", [ix, iy, iz])?;
        require_finite3("offset", [rx, ry, rz])?;
        if body.inv_mass > 0.0 {
            let impulse = [ix, iy, iz];
            body.velocity = madd(body.velocity, impulse, body.inv_mass);
            let angular_impulse = cross([rx, ry, rz], impulse);
            let delta_w = mat_vec(&body.world_inv_inertia(), angular_impulse);
            body.angular_velocity = add(body.angular_velocity, delta_w);
            body.sleeping = false;
            body.sleep_timer = 0.0;
        }
        Ok(())
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

    /// Return the handles of all active bodies as a `Vec<u32>`.
    pub fn get_all_body_handles(&self) -> Vec<u32> {
        self.bodies
            .iter()
            .enumerate()
            .filter(|(_, b)| b.active)
            .map(|(i, b)| encode_handle(i, b.generation))
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

    /// Return the contact results from the last `step()` call that ran at
    /// least one substep.
    ///
    /// Every fixed substep contributes its own contacts, so a `step()` that
    /// catches up with several substeps lists a pair that stays in contact
    /// once per substep; a call that is too short to fill a substep leaves the
    /// previous list untouched.
    pub fn get_contacts(&self) -> &[ContactResult] {
        &self.contacts
    }

    /// Return the number of entries of [`get_contacts`](Self::get_contacts).
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
    ///
    /// # Errors
    ///
    /// See [`add_collider`](Self::add_collider).
    pub fn add_sphere_collider(&mut self, body: u32, radius: f64) -> Result<u32> {
        let cfg = ColliderConfig::sphere(radius);
        self.add_collider(body, &cfg)
    }

    /// Add an oriented box collider (half-extents along the body's local axes)
    /// to a body and return the collider handle.
    ///
    /// # Errors
    ///
    /// See [`add_collider`](Self::add_collider).
    pub fn add_box_collider(&mut self, body: u32, hx: f64, hy: f64, hz: f64) -> Result<u32> {
        let cfg = ColliderConfig::cuboid(hx, hy, hz);
        self.add_collider(body, &cfg)
    }

    /// Add a capsule collider to a body and return the collider handle.
    ///
    /// The capsule is symmetric about the body's local Y axis; `height` is the
    /// length of the cylindrical part (the hemispherical caps come on top).
    ///
    /// # Errors
    ///
    /// See [`add_collider`](Self::add_collider).
    pub fn add_capsule_collider(&mut self, body: u32, radius: f64, height: f64) -> Result<u32> {
        let cfg = ColliderConfig::capsule(radius, height);
        self.add_collider(body, &cfg)
    }

    /// Add a cylinder collider to a body and return the collider handle.
    ///
    /// The cylinder is symmetric about the body's local Y axis.
    ///
    /// # Errors
    ///
    /// See [`add_collider`](Self::add_collider).
    pub fn add_cylinder_collider(&mut self, body: u32, radius: f64, height: f64) -> Result<u32> {
        let cfg = ColliderConfig::cylinder(radius, height);
        self.add_collider(body, &cfg)
    }

    /// Add a plane collider (infinite half-space boundary `n . x = offset`).
    ///
    /// The normal is normalised (and `offset` rescaled accordingly) and, like
    /// every collider, the plane is expressed in the frame of its body.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidParameter`] for a zero or non-finite `plane_normal`,
    /// plus the errors of [`add_collider`](Self::add_collider).
    pub fn add_plane_collider(
        &mut self,
        body: u32,
        nx: f64,
        ny: f64,
        nz: f64,
        offset: f64,
    ) -> Result<u32> {
        let cfg = ColliderConfig::plane(nx, ny, nz, offset);
        self.add_collider(body, &cfg)
    }

    /// Add a collider using a full `ColliderConfig` and return the collider handle.
    ///
    /// The configuration is validated and a normalised copy is stored (unit
    /// plane normal, unit local rotation).  The rotational inertia of the
    /// owning body is recomputed from all its colliders.
    ///
    /// # Errors
    ///
    /// * [`Error::BodyNotFound`] - `body_handle` does not refer to a live body;
    /// * [`Error::InvalidDimension`] - a radius, half-extent or height is not
    ///   finite and strictly positive;
    /// * [`Error::InvalidParameter`] - a zero/non-finite plane normal, a
    ///   non-finite local offset, a zero local rotation or an invalid density;
    /// * [`Error::UnsupportedShape`] - a convex hull (its vertices are not part
    ///   of [`ColliderConfig`]);
    /// * [`Error::CapacityExceeded`] - all [`MAX_SLOTS`] collider slots are in
    ///   use.
    pub fn add_collider(&mut self, body_handle: u32, config: &ColliderConfig) -> Result<u32> {
        if self.body_slot(body_handle).is_none() {
            return Err(Error::BodyNotFound(body_handle));
        }
        let config = validate_collider_config(config)?;
        let shape = Shape::from_config(&config).ok_or_else(|| {
            Error::UnsupportedShape(format!(
                "{:?} colliders are not supported",
                config.shape_type
            ))
        })?;
        let local_rotation = quat_to_mat3(config.local_rotation);
        let make_entry = |generation: u32| ColliderEntry {
            body_handle,
            config: config.clone(),
            shape,
            local_rotation,
            active: true,
            generation,
        };
        let handle = if let Some(idx) = self.collider_free_list.pop() {
            let slot = self
                .colliders
                .get_mut(idx)
                .ok_or_else(|| Error::internal("collider free list references a missing slot"))?;
            let generation = next_generation(slot.generation);
            *slot = make_entry(generation);
            encode_handle(idx, generation)
        } else {
            if self.colliders.len() >= self.slot_limit {
                return Err(Error::CapacityExceeded {
                    max: self.slot_limit,
                });
            }
            let idx = self.colliders.len();
            self.colliders.push(make_entry(0));
            encode_handle(idx, 0)
        };
        self.refresh_body_inertia(body_handle);
        Ok(handle)
    }

    /// Add a body together with one collider and return both handles
    /// `(body, collider)`.
    ///
    /// The operation is atomic: both configurations are validated before
    /// anything is created, and if the collider cannot be attached after all
    /// (for example because the collider table is full) the freshly created
    /// body is removed again, so a failure never leaves a collider-less body
    /// behind.
    ///
    /// # Errors
    ///
    /// The errors of [`add_rigid_body`](Self::add_rigid_body) and
    /// [`add_collider`](Self::add_collider) (with the exception of
    /// [`Error::BodyNotFound`], which cannot occur).
    pub fn add_body_with_collider(
        &mut self,
        body: &RigidBodyConfig,
        collider: &ColliderConfig,
    ) -> Result<(u32, u32)> {
        validate_body_config(body)?;
        validate_collider_config(collider)?;
        let body_handle = self.add_rigid_body(body)?;
        match self.add_collider(body_handle, collider) {
            Ok(collider_handle) => Ok((body_handle, collider_handle)),
            Err(error) => {
                self.remove_body(body_handle).map_err(|rollback| {
                    Error::internal(format!(
                        "could not roll back body {body_handle} after `{error}`: {rollback}"
                    ))
                })?;
                Err(error)
            }
        }
    }

    /// Remove the collider with the given handle.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidHandle`] if the handle does not refer to a live
    /// collider (never allocated, already removed, stale, or removed together
    /// with its body).
    pub fn remove_collider(&mut self, handle: u32) -> Result<()> {
        let (idx, generation) = decode_handle(handle);
        let owner = match self.colliders.get_mut(idx) {
            Some(col) if col.active && col.generation == generation => {
                col.active = false;
                col.body_handle
            }
            _ => return Err(Error::InvalidHandle(handle)),
        };
        self.collider_free_list.push(idx);
        self.refresh_body_inertia(owner);
        Ok(())
    }

    /// Recompute the rotational inertia of a body from its active colliders.
    fn refresh_body_inertia(&mut self, body_handle: u32) {
        let Some(slot) = self.body_slot(body_handle) else {
            return;
        };
        let parts: Vec<ShapePart> = self
            .colliders
            .iter()
            .filter(|c| c.active && c.body_handle == body_handle)
            .map(|c| ShapePart {
                shape: c.shape,
                position: c.config.local_position,
                rotation: c.local_rotation,
            })
            .collect();
        let Some(body) = self.bodies.get_mut(slot) else {
            return;
        };
        let tensor = if body.is_dynamic() {
            compound_inertia(body.mass, &parts)
        } else {
            None
        };
        match tensor.and_then(|t| mat_inverse(&t).map(|inverse| (t, inverse))) {
            Some((inertia, inverse)) => {
                body.inertia_local = inertia;
                body.inv_inertia_local = inverse;
                body.has_inertia = true;
            }
            None => {
                body.inertia_local = MAT3_ZERO;
                body.inv_inertia_local = MAT3_ZERO;
                body.has_inertia = false;
            }
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
    ///
    /// # Errors
    ///
    /// [`Error::InvalidTimeStep`] unless `dt` is finite and strictly positive
    /// (a negative step would integrate backwards).
    pub fn set_fixed_dt(&mut self, dt: f64) -> Result<()> {
        if !(dt.is_finite() && dt > 0.0) {
            return Err(Error::InvalidTimeStep(dt));
        }
        self.config.fixed_dt = dt;
        Ok(())
    }

    /// Set the linear damping on a body (0 = no damping).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidHandle`] for a dead or stale handle,
    /// [`Error::InvalidParameter`] unless `damping` is finite and `>= 0` (a
    /// negative damping would make the body speed up on its own).
    pub fn set_body_linear_damping(&mut self, handle: u32, damping: f64) -> Result<()> {
        let body = self.get_body_mut(handle)?;
        non_negative_param("linear_damping", damping)?;
        body.linear_damping = damping;
        Ok(())
    }

    /// Set the angular damping on a body (0 = no damping).
    ///
    /// # Errors
    ///
    /// As [`set_body_linear_damping`](Self::set_body_linear_damping).
    pub fn set_body_angular_damping(&mut self, handle: u32, damping: f64) -> Result<()> {
        let body = self.get_body_mut(handle)?;
        non_negative_param("angular_damping", damping)?;
        body.angular_damping = damping;
        Ok(())
    }

    // -----------------------------------------------------------------------
    // Raycasting
    // -----------------------------------------------------------------------

    /// Cast a ray from `origin` in `direction` and return the closest hit.
    ///
    /// Every collider shape is intersected exactly: spheres, oriented boxes,
    /// capsules, cylinders and (one-sided) planes.  A ray that starts inside a
    /// solid collider reports the surface it leaves through; hits farther than
    /// `max_distance` are ignored.  `direction` need not be normalised; a
    /// zero or non-finite direction, a non-finite origin or a non-positive
    /// `max_distance` yield no hit.
    ///
    /// The result's `normal` is the outward surface normal and `distance` is
    /// measured along the normalised direction.
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
        let origin = [ox, oy, oz];
        let direction = [dx, dy, dz];
        let mut best = RaycastResult::no_hit();
        if !origin.iter().all(|v| v.is_finite())
            || !direction.iter().all(|v| v.is_finite())
            || max_distance.is_nan()
            || max_distance <= 0.0
        {
            return best;
        }
        let Some(dir) = normalized(direction) else {
            return best;
        };

        for col in &self.colliders {
            if !col.active {
                continue;
            }
            let Some(body) = self.get_body(col.body_handle) else {
                continue;
            };
            let placed = Self::collider_world(col, body);
            if let Some(hit) = ray_shape(origin, dir, &placed)
                && hit.distance <= max_distance
                && hit.distance < best.distance
            {
                best = RaycastResult {
                    hit: true,
                    body_handle: col.body_handle,
                    point: point_at(origin, dir, hit.distance),
                    normal: hit.normal,
                    distance: hit.distance,
                };
            }
        }

        best
    }

    // -----------------------------------------------------------------------
    // Simulation step
    // -----------------------------------------------------------------------

    /// Advance the simulation by `dt` seconds of real time.
    ///
    /// `dt` is added to a fixed-step accumulator and the engine runs as many
    /// substeps of exactly [`SimulationConfig::fixed_dt`] as the accumulator
    /// allows, at most [`SimulationConfig::max_substeps`] per call; the
    /// simulation clock ([`time`](Self::time)) advances by the time that was
    /// actually integrated.  Time that does not fill a whole substep stays in
    /// the accumulator (see [`interpolation_alpha`](Self::interpolation_alpha)),
    /// so a call with `dt < fixed_dt` may integrate nothing at all and the
    /// contacts reported by [`get_contacts`](Self::get_contacts) are then those
    /// of the previous step.  If the substep cap is reached the remaining time
    /// is discarded and reported in
    /// [`DebugInfo::dropped_time`](crate::types::DebugInfo::dropped_time)
    /// rather than being carried into later frames.
    ///
    /// Forces and torques applied before this call act during every substep of
    /// it and are cleared afterwards.  A `dt` of exactly zero is a no-op.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidTimeStep`] if `dt` is negative, `NaN` or infinite (the
    /// engine state is left untouched).
    pub fn step(&mut self, dt: f64) -> Result<()> {
        if !dt.is_finite() || dt < 0.0 {
            return Err(Error::InvalidTimeStep(dt));
        }
        if dt == 0.0 {
            return Ok(());
        }

        let fixed_dt = self.config.fixed_dt;
        let max_substeps = self.config.max_substeps;
        let threshold = fixed_dt * (1.0 - STEP_EPSILON);

        self.accumulator += dt;
        if self.accumulator >= threshold {
            // At least one substep will run, so this call defines the new
            // contact list.
            self.contacts.clear();
        }

        let mut substeps = 0u32;
        let mut sweeps = 0u32;
        let mut candidate_pairs = 0u32;
        let mut narrowphase_tests = 0u32;
        while self.accumulator >= threshold && substeps < max_substeps {
            self.integrate(fixed_dt);
            let pass = self.detect_and_resolve_contacts(fixed_dt);
            candidate_pairs += pass.candidate_pairs;
            narrowphase_tests += pass.narrowphase_tests;
            if pass.contacts > 0 {
                sweeps += 1;
            }
            self.accumulator = (self.accumulator - fixed_dt).max(0.0);
            self.time += fixed_dt;
            substeps += 1;
        }

        let mut dropped = 0.0;
        if self.accumulator >= threshold {
            dropped = self.accumulator;
            self.accumulator = 0.0;
        }

        for body in &mut self.bodies {
            body.force = [0.0; 3];
            body.torque = [0.0; 3];
        }

        self.debug_info.substeps_performed = substeps;
        self.debug_info.dropped_time = dropped;
        if substeps > 0 {
            self.debug_info.broadphase_pairs = candidate_pairs;
            self.debug_info.narrowphase_tests = narrowphase_tests;
            self.debug_info.active_contacts = self.contacts.len() as u32;
            self.debug_info.solver_iterations_performed = sweeps;
        }
        Ok(())
    }

    /// Fraction of a fixed substep that is still pending in the accumulator,
    /// in `[0, 1)`.
    ///
    /// Renderers can blend the previous and the current state with this value
    /// (`state = prev + alpha * (current - prev)`) to display smooth motion
    /// when the render rate is not a multiple of the physics rate.
    pub fn interpolation_alpha(&self) -> f64 {
        (self.accumulator / self.config.fixed_dt).clamp(0.0, 1.0)
    }

    // -----------------------------------------------------------------------
    // Reset
    // -----------------------------------------------------------------------

    /// Remove all bodies, colliders and contacts, resetting to an empty state.
    ///
    /// The simulation clock and the fixed-step accumulator return to zero.
    /// The slot tables keep their generation counters (every slot simply goes
    /// back to the free lists), so a handle obtained before the reset is
    /// rejected afterwards instead of silently addressing a new body.
    pub fn reset(&mut self) {
        for body in &mut self.bodies {
            body.active = false;
        }
        for collider in &mut self.colliders {
            collider.active = false;
        }
        // Free lists are popped from the back: reuse the lowest slots first.
        self.body_free_list = (0..self.bodies.len()).rev().collect();
        self.collider_free_list = (0..self.colliders.len()).rev().collect();
        self.contacts.clear();
        self.time = 0.0;
        self.accumulator = 0.0;
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

    /// Live body a handle refers to (`None` for dead or stale handles).
    fn get_body(&self, handle: u32) -> Option<&BodyEntry> {
        let (idx, generation) = decode_handle(handle);
        self.bodies
            .get(idx)
            .filter(|b| b.active && b.generation == generation)
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
        let (idx, generation) = decode_handle(handle);
        self.bodies
            .get_mut(idx)
            .filter(|b| b.active && b.generation == generation)
            .ok_or(Error::InvalidHandle(handle))
    }

    /// Place a collider into the world using the pose of its body.
    fn collider_world(collider: &ColliderEntry, body: &BodyEntry) -> WorldShape {
        let rotation = quat_to_mat3(body.rotation);
        WorldShape::place(
            collider.shape,
            body.position,
            &rotation,
            collider.config.local_position,
            &collider.local_rotation,
        )
    }

    /// Integrate velocities, orientations and positions over one fixed
    /// substep (symplectic Euler).
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
                if length_sq(body.force) > 1e-20 || length_sq(body.torque) > 1e-20 {
                    body.sleeping = false;
                    body.sleep_timer = 0.0;
                } else {
                    continue;
                }
            }

            // 1. Gravity (scaled by gravity_scale) and the applied force give
            //    the linear acceleration.
            let gs = body.gravity_scale;
            let acceleration = madd(
                [gravity[0] * gs, gravity[1] * gs, gravity[2] * gs],
                body.force,
                body.inv_mass,
            );

            // 2. Update velocity (symplectic Euler: update v first)
            body.velocity = madd(body.velocity, acceleration, dt);

            // 3. Applied torque and gyroscopic acceleration change the angular
            //    velocity (before damping).
            body.accelerate_rotation(dt);

            // 4. Apply damping
            let ld = (1.0 - body.linear_damping * dt).max(0.0);
            let ad = (1.0 - body.angular_damping * dt).max(0.0);
            body.velocity = scale(body.velocity, ld);
            body.angular_velocity = scale(body.angular_velocity, ad);

            // 5. Update position
            body.position = madd(body.position, body.velocity, dt);

            // 6. Update orientation from the world-frame angular velocity
            body.integrate_orientation(dt);

            // 7. Update sleep timer
            if sleeping_enabled {
                let v2 = length_sq(body.velocity);
                let w2 = length_sq(body.angular_velocity);
                if v2 < lin_sleep * lin_sleep && w2 < ang_sleep * ang_sleep {
                    body.sleep_timer += dt;
                    if body.sleep_timer > time_sleep {
                        body.sleeping = true;
                    }
                } else {
                    body.sleep_timer = 0.0;
                }
            }
        }
    }

    /// Run the narrow phase on all collider pairs and resolve every contact
    /// found with impulses (normal, friction) and positional correction.
    ///
    /// Contacts are appended to `self.contacts`.  Contact impulses act on the
    /// linear velocity only; they do not exert torque.
    fn detect_and_resolve_contacts(&mut self, dt: f64) -> PassStats {
        let (manifolds, mut stats) = self.build_manifolds();
        for manifold in &manifolds {
            let contact = self.resolve_manifold(manifold, dt);
            self.contacts.push(contact);
        }
        stats.contacts = manifolds.len() as u32;
        stats
    }

    /// Build contact manifolds by testing all pairs of colliders that belong
    /// to different bodies.
    ///
    /// Each pair is ordered by shape rank so that, for instance, the sphere of
    /// a sphere-plane pair is always body A; colliders of equal rank keep
    /// their creation order.  The returned normals point from body B toward
    /// body A.
    fn build_manifolds(&self) -> (Vec<ManifoldPair>, PassStats) {
        let placed: Vec<(u32, WorldShape)> = self
            .colliders
            .iter()
            .filter(|c| c.active)
            .filter_map(|c| {
                let body = self.get_body(c.body_handle)?;
                Some((c.body_handle, Self::collider_world(c, body)))
            })
            .collect();

        let mut stats = PassStats::default();
        let mut manifolds = Vec::new();
        for (i, (handle_i, shape_i)) in placed.iter().enumerate() {
            for (handle_j, shape_j) in placed.iter().skip(i + 1) {
                if handle_i == handle_j {
                    continue;
                }
                stats.candidate_pairs += 1;
                if !bounds_overlap(shape_i, shape_j) {
                    continue;
                }
                stats.narrowphase_tests += 1;
                let ((body_a, shape_a), (body_b, shape_b)) =
                    if shape_i.shape.rank() <= shape_j.shape.rank() {
                        ((*handle_i, shape_i), (*handle_j, shape_j))
                    } else {
                        ((*handle_j, shape_j), (*handle_i, shape_i))
                    };
                if let Some(contact) = collide(shape_a, shape_b) {
                    manifolds.push(ManifoldPair {
                        body_a,
                        body_b,
                        normal: contact.normal,
                        depth: contact.depth,
                        point_a: contact.point_a,
                        point_b: contact.point_b,
                    });
                }
            }
        }
        (manifolds, stats)
    }

    /// `(inverse mass used for contacts, velocity, restitution, friction)` of
    /// the body in a slot; only dynamic bodies respond to contacts.
    fn contact_properties(&self, slot: usize) -> (f64, [f64; 3], f64, f64) {
        self.bodies
            .get(slot)
            .map_or((0.0, [0.0; 3], 0.0, 0.5), |b| {
                let inv_mass = if b.is_dynamic() { b.inv_mass } else { 0.0 };
                (inv_mass, b.velocity, b.restitution, b.friction)
            })
    }

    /// Change the linear velocity and position of a dynamic body (contact
    /// response); other bodies are left untouched.
    fn nudge_body(&mut self, slot: usize, velocity_delta: Vec3, position_delta: Vec3) {
        if let Some(body) = self.bodies.get_mut(slot)
            && body.is_dynamic()
        {
            body.velocity = add(body.velocity, velocity_delta);
            body.position = add(body.position, position_delta);
        }
    }

    /// Resolve a single contact manifold using impulse-based collision response.
    ///
    /// * the normal impulse removes the approach velocity along the contact
    ///   normal (with restitution `sqrt(e_a * e_b)`);
    /// * a Coulomb friction impulse (`mu = sqrt(mu_a * mu_b)`) acts along the
    ///   tangential relative velocity;
    /// * a Baumgarte-style positional correction pushes the bodies apart by
    ///   20 % of the penetration beyond `contact_slop` per substep - also for
    ///   pairs that are already separating, so an overlap can never persist.
    ///
    /// The normal points from body B toward body A, so approaching bodies have
    /// a negative relative normal velocity and a positive impulse pushes A
    /// along `+normal` and B along `-normal`.
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

        let (Some(ia), Some(ib)) = (self.body_slot(m.body_a), self.body_slot(m.body_b)) else {
            return cr;
        };
        let (inv_mass_a, vel_a, rest_a, fric_a) = self.contact_properties(ia);
        let (inv_mass_b, vel_b, rest_b, fric_b) = self.contact_properties(ib);
        let total_inv = inv_mass_a + inv_mass_b;
        if total_inv == 0.0 {
            return cr;
        }

        let n = m.normal;
        let rel_vn = dot(sub(vel_a, vel_b), n);
        cr.relative_velocity = rel_vn;

        if rel_vn <= 0.0 {
            // Normal impulse (approaching or resting contact).
            let e = (rest_a * rest_b).sqrt();
            let j = -(1.0 + e) * rel_vn / total_inv;
            cr.impulse = j;
            self.nudge_body(ia, scale(n, j * inv_mass_a), [0.0; 3]);
            self.nudge_body(ib, scale(n, -j * inv_mass_b), [0.0; 3]);

            // Friction impulse (Coulomb model) on the post-impulse velocities.
            let mu = (fric_a * fric_b).sqrt();
            let (_, vel_a2, _, _) = self.contact_properties(ia);
            let (_, vel_b2, _, _) = self.contact_properties(ib);
            let rel2 = sub(vel_a2, vel_b2);
            let tangential = sub(rel2, scale(n, dot(rel2, n)));
            let t_len = length(tangential);
            if t_len > 1e-10 {
                let t_dir = scale(tangential, 1.0 / t_len);
                let jt_raw = -dot(rel2, t_dir) / total_inv;
                let jt = jt_raw.clamp(-mu * j.abs(), mu * j.abs());
                cr.friction_impulse = jt;
                self.nudge_body(ia, scale(t_dir, jt * inv_mass_a), [0.0; 3]);
                self.nudge_body(ib, scale(t_dir, -jt * inv_mass_b), [0.0; 3]);
            }
        }

        // Positional correction (Baumgarte).
        let slop = self.config.contact_slop;
        let corr_mag = (m.depth - slop).max(0.0) * 0.2 / total_inv;
        if corr_mag > 1e-10 {
            self.nudge_body(ia, [0.0; 3], scale(n, corr_mag * inv_mass_a));
            self.nudge_body(ib, [0.0; 3], scale(n, -corr_mag * inv_mass_b));
        }

        cr
    }
}
