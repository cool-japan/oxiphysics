// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! High-level physics world API designed for Python consumers.
//!
//! `PyPhysicsWorld` provides the primary interface: add/remove bodies, apply
//! forces and impulses, step the simulation, and query state. Handles are
//! `u32` integers for easy FFI transmission. All state uses plain arrays and
//! primitive types; no nalgebra types are exposed through the public API.

// Sub-modules
mod constraints;
mod extensions;
mod fem;
mod geometry;
mod lbm;
mod materials;
mod math_helpers;
mod py_bindings;
mod sph;
mod stats;

#[cfg(test)]
mod tests;

// Re-exports
pub use constraints::{ConstraintType, PyConstraint};
pub use extensions::{ContactPair, InertiaTensor, PyRigidBody};
pub use fem::{FemBarElement, PyFemAssembly};
pub use geometry::{PyAabb, PyConvexHull, PySphere};
pub use lbm::{PyLbmConfig, PyLbmGrid};
pub use materials::{MaterialClass, PyMaterial};
pub use math_helpers::{
    array_to_vec3, quat_conjugate, quat_from_axis_angle, quat_mul, quat_normalize, quat_rotate_vec,
    vec3_to_array,
};
pub use py_bindings::{
    MaterialProperties, PyFemBinding, PyLbmBinding, PyMdBinding, PySphBinding, py_p_wave_speed,
    py_query_material, py_s_wave_speed,
};
pub use sph::{PySphConfig, PySphSim};
pub use stats::SimStats;

use crate::types::{
    PyColliderShape, PyContactResult, PyRigidBodyConfig, PyRigidBodyDesc, PySimConfig, PyVec3,
};
use pyo3::prelude::*;
use serde::{Deserialize, Serialize};

// ---------------------------------------------------------------------------
// Internal body state
// ---------------------------------------------------------------------------

/// Full internal state of a single rigid body.
#[derive(Debug, Clone, Serialize, Deserialize)]
struct InternalBody {
    /// Position `[x, y, z]`.
    position: [f64; 3],
    /// Linear velocity `[vx, vy, vz]`.
    velocity: [f64; 3],
    /// Orientation quaternion `[x, y, z, w]`.
    orientation: [f64; 4],
    /// Angular velocity `[wx, wy, wz]`.
    angular_velocity: [f64; 3],
    /// Mass in kilograms.
    mass: f64,
    /// Whether this body is static (immovable).
    is_static: bool,
    /// Whether this body is kinematic (moved manually).
    is_kinematic: bool,
    /// Whether this body is currently sleeping.
    is_sleeping: bool,
    /// Whether this body participates in sleeping.
    can_sleep: bool,
    /// Friction coefficient.
    friction: f64,
    /// Restitution coefficient.
    restitution: f64,
    /// Linear damping.
    linear_damping: f64,
    /// Angular damping.
    angular_damping: f64,
    /// Accumulated force for this step `[fx, fy, fz]`.
    accumulated_force: [f64; 3],
    /// Accumulated torque for this step `[tx, ty, tz]`.
    accumulated_torque: [f64; 3],
    /// Velocity-below-threshold counter (for sleep detection).
    sleep_timer: f64,
    /// Attached collider shapes.
    shapes: Vec<PyColliderShape>,
    /// Optional tag.
    tag: Option<String>,
}

impl InternalBody {
    fn from_config(config: &PyRigidBodyConfig) -> Self {
        Self {
            position: config.position,
            velocity: config.velocity,
            orientation: config.orientation,
            angular_velocity: config.angular_velocity,
            mass: config.mass,
            is_static: config.is_static,
            is_kinematic: config.is_kinematic,
            is_sleeping: false,
            can_sleep: config.can_sleep,
            friction: config.friction,
            restitution: config.restitution,
            linear_damping: config.linear_damping,
            angular_damping: config.angular_damping,
            accumulated_force: [0.0; 3],
            accumulated_torque: [0.0; 3],
            sleep_timer: 0.0,
            shapes: config.shapes.clone(),
            tag: config.tag.clone(),
        }
    }

    /// Effective inverse mass (0 for static or infinite-mass).
    fn inv_mass(&self) -> f64 {
        if self.is_static || self.is_kinematic || self.mass <= 0.0 {
            0.0
        } else {
            1.0 / self.mass
        }
    }

    /// Whether this body is mobile (can be accelerated by forces).
    fn is_dynamic(&self) -> bool {
        !self.is_static && !self.is_kinematic
    }
}

// ---------------------------------------------------------------------------
// Slot-based storage with generation counters
// ---------------------------------------------------------------------------

#[derive(Debug, Clone)]
struct Slot {
    body: Option<InternalBody>,
    generation: u32,
}

// ---------------------------------------------------------------------------
// PyPhysicsWorld
// ---------------------------------------------------------------------------

/// A self-contained physics simulation world.
///
/// Provides Python-friendly methods using integer handles (`u32`) and plain
/// array types. Internally uses a slot-based arena with generation counters
/// so that stale handles are detected as absent bodies.
#[pyclass(from_py_object)]
#[derive(Debug, Clone)]
pub struct PyPhysicsWorld {
    /// Body storage slots.
    slots: Vec<Slot>,
    /// Free-list of slot indices available for reuse.
    free_list: Vec<u32>,
    /// Global simulation configuration.
    config: PySimConfig,
    /// Accumulated simulation time.
    time: f64,
    /// Contacts detected in the most recent step.
    contacts: Vec<PyContactResult>,
    /// Active constraints.
    constraints: Vec<PyConstraint>,
}

#[pymethods]
impl PyPhysicsWorld {
    // -----------------------------------------------------------------------
    // Construction
    // -----------------------------------------------------------------------

    /// Create a new physics world with the given gravity vector.
    ///
    /// By default uses Earth gravity `(0, -9.81, 0)`.
    #[new]
    #[pyo3(signature = (gx = 0.0, gy = -9.81, gz = 0.0))]
    pub fn py_new(gx: f64, gy: f64, gz: f64) -> Self {
        Self::new_with_config(PySimConfig {
            gravity: [gx, gy, gz],
            ..PySimConfig::default()
        })
    }

    // -----------------------------------------------------------------------
    // Body management
    // -----------------------------------------------------------------------

    /// Remove a body by handle. Returns `true` if the body existed.
    pub fn remove_body(&mut self, handle: u32) -> bool {
        if let Some(slot) = self.slots.get_mut(handle as usize)
            && slot.body.is_some()
        {
            slot.body = None;
            slot.generation = slot.generation.wrapping_add(1);
            self.free_list.push(handle);
            return true;
        }
        false
    }

    // -----------------------------------------------------------------------
    // Querying body state
    // -----------------------------------------------------------------------

    /// Get the position of a body, or `None` if the handle is invalid.
    pub fn get_position(&self, handle: u32) -> Option<[f64; 3]> {
        self.get_body(handle).map(|b| b.position)
    }

    /// Get the linear velocity of a body, or `None` if the handle is invalid.
    pub fn get_velocity(&self, handle: u32) -> Option<[f64; 3]> {
        self.get_body(handle).map(|b| b.velocity)
    }

    /// Get the orientation quaternion `[x, y, z, w]`, or `None`.
    pub fn get_orientation(&self, handle: u32) -> Option<[f64; 4]> {
        self.get_body(handle).map(|b| b.orientation)
    }

    /// Get the angular velocity, or `None`.
    pub fn get_angular_velocity(&self, handle: u32) -> Option<[f64; 3]> {
        self.get_body(handle).map(|b| b.angular_velocity)
    }

    /// Whether the body with `handle` is currently sleeping.
    pub fn is_sleeping(&self, handle: u32) -> bool {
        self.get_body(handle)
            .map(|b| b.is_sleeping)
            .unwrap_or(false)
    }

    /// Get the tag of a body, or `None`.
    pub fn get_tag(&self, handle: u32) -> Option<String> {
        self.get_body(handle).and_then(|b| b.tag.clone())
    }

    // -----------------------------------------------------------------------
    // Mutating body state
    // -----------------------------------------------------------------------

    /// Set the position of a body.
    pub fn set_position(&mut self, handle: u32, pos: [f64; 3]) {
        if let Some(body) = self.get_body_mut(handle) {
            body.position = pos;
            body.is_sleeping = false;
        }
    }

    /// Set the linear velocity of a body.
    pub fn set_velocity(&mut self, handle: u32, vel: [f64; 3]) {
        if let Some(body) = self.get_body_mut(handle) {
            body.velocity = vel;
            body.is_sleeping = false;
        }
    }

    /// Set the orientation (quaternion `[x, y, z, w]`) of a body.
    pub fn set_orientation(&mut self, handle: u32, orientation: [f64; 4]) {
        if let Some(body) = self.get_body_mut(handle) {
            body.orientation = normalize_quat(orientation);
            body.is_sleeping = false;
        }
    }

    /// Set the angular velocity of a body.
    pub fn set_angular_velocity(&mut self, handle: u32, omega: [f64; 3]) {
        if let Some(body) = self.get_body_mut(handle) {
            body.angular_velocity = omega;
            body.is_sleeping = false;
        }
    }

    // -----------------------------------------------------------------------
    // Force / impulse API
    // -----------------------------------------------------------------------

    /// Apply a force `[fx, fy, fz]` to a body, optionally at a world-space point.
    ///
    /// If `point` is `None`, the force is applied at the center of mass.
    /// Force accumulates until the next `step()`.
    pub fn apply_force(&mut self, handle: u32, force: [f64; 3], point: Option<[f64; 3]>) {
        if let Some(body) = self.get_body_mut(handle) {
            if !body.is_dynamic() {
                return;
            }
            body.accumulated_force[0] += force[0];
            body.accumulated_force[1] += force[1];
            body.accumulated_force[2] += force[2];

            if let Some(pt) = point {
                // Torque = (pt - center) x force
                let r = [
                    pt[0] - body.position[0],
                    pt[1] - body.position[1],
                    pt[2] - body.position[2],
                ];
                let torque = cross3(r, force);
                body.accumulated_torque[0] += torque[0];
                body.accumulated_torque[1] += torque[1];
                body.accumulated_torque[2] += torque[2];
            }
            body.is_sleeping = false;
        }
    }

    /// Apply an impulse `[ix, iy, iz]` directly to a body's velocity.
    ///
    /// An impulse is a instantaneous change in momentum: `Δv = impulse / mass`.
    /// If `point` is provided, an angular impulse is also applied.
    pub fn apply_impulse(&mut self, handle: u32, impulse: [f64; 3], point: Option<[f64; 3]>) {
        if let Some(body) = self.get_body_mut(handle) {
            if !body.is_dynamic() {
                return;
            }
            let inv_m = body.inv_mass();
            body.velocity[0] += impulse[0] * inv_m;
            body.velocity[1] += impulse[1] * inv_m;
            body.velocity[2] += impulse[2] * inv_m;

            if let Some(pt) = point {
                // Angular impulse approximation using scalar inertia
                let r = [
                    pt[0] - body.position[0],
                    pt[1] - body.position[1],
                    pt[2] - body.position[2],
                ];
                let ang_impulse = cross3(r, impulse);
                // Approximate moment of inertia: 0.4 * m * r^2 (solid sphere)
                // Better approximation from shape if available
                let approx_inertia = approximate_inertia(body);
                let inv_i = if approx_inertia > 1e-12 {
                    1.0 / approx_inertia
                } else {
                    0.0
                };
                body.angular_velocity[0] += ang_impulse[0] * inv_i;
                body.angular_velocity[1] += ang_impulse[1] * inv_i;
                body.angular_velocity[2] += ang_impulse[2] * inv_i;
            }
            body.is_sleeping = false;
        }
    }

    /// Apply a torque `[tx, ty, tz]` to a body (accumulates until next step).
    pub fn apply_torque(&mut self, handle: u32, torque: [f64; 3]) {
        if let Some(body) = self.get_body_mut(handle) {
            if !body.is_dynamic() {
                return;
            }
            body.accumulated_torque[0] += torque[0];
            body.accumulated_torque[1] += torque[1];
            body.accumulated_torque[2] += torque[2];
            body.is_sleeping = false;
        }
    }

    /// Wake up a sleeping body.
    pub fn wake_body(&mut self, handle: u32) {
        if let Some(body) = self.get_body_mut(handle) {
            body.is_sleeping = false;
            body.sleep_timer = 0.0;
        }
    }

    /// Put a body to sleep.
    pub fn sleep_body(&mut self, handle: u32) {
        if let Some(body) = self.get_body_mut(handle) {
            body.is_sleeping = true;
            body.velocity = [0.0; 3];
            body.angular_velocity = [0.0; 3];
        }
    }

    // -----------------------------------------------------------------------
    // Global configuration
    // -----------------------------------------------------------------------

    /// Set the gravity vector.
    pub fn set_gravity(&mut self, g: [f64; 3]) {
        self.config.gravity = g;
    }

    /// Get the current gravity vector.
    pub fn gravity(&self) -> [f64; 3] {
        self.config.gravity
    }

    /// Get the gravity and key config values as `[gx, gy, gz]`.
    pub fn gravity_vec(&self) -> [f64; 3] {
        self.config.gravity
    }

    // -----------------------------------------------------------------------
    // Statistics
    // -----------------------------------------------------------------------

    /// Total number of allocated body slots (includes removed ones until reuse).
    pub fn body_count(&self) -> usize {
        self.slots.iter().filter(|s| s.body.is_some()).count()
    }

    /// Number of bodies currently sleeping.
    pub fn sleeping_count(&self) -> usize {
        self.slots
            .iter()
            .filter_map(|s| s.body.as_ref())
            .filter(|b| b.is_sleeping)
            .count()
    }

    /// Accumulated simulation time in seconds.
    pub fn time(&self) -> f64 {
        self.time
    }

    /// Total number of contacts from the most recent step.
    pub fn get_contact_count(&self) -> usize {
        self.contacts.len()
    }

    // -----------------------------------------------------------------------
    // Positions / velocities bulk query
    // -----------------------------------------------------------------------

    /// Return all body positions as a flat `Vec<[f64;3]>` in handle order (skipping removed).
    pub fn all_positions(&self) -> Vec<[f64; 3]> {
        self.slots
            .iter()
            .filter_map(|s| s.body.as_ref().map(|b| b.position))
            .collect()
    }

    /// Return all body velocities in handle order (skipping removed).
    pub fn all_velocities(&self) -> Vec<[f64; 3]> {
        self.slots
            .iter()
            .filter_map(|s| s.body.as_ref().map(|b| b.velocity))
            .collect()
    }

    /// Return all active (non-removed) handles.
    pub fn active_handles(&self) -> Vec<u32> {
        self.slots
            .iter()
            .enumerate()
            .filter(|(_, s)| s.body.is_some())
            .map(|(i, _)| i as u32)
            .collect()
    }

    // -----------------------------------------------------------------------
    // Simulation step
    // -----------------------------------------------------------------------

    /// Advance the simulation by `dt` seconds.
    ///
    /// Uses symplectic Euler integration:
    /// 1. Accumulate gravity + applied forces.
    /// 2. Update velocity from forces.
    /// 3. Apply damping.
    /// 4. Update position from velocity.
    /// 5. Integrate orientation from angular velocity.
    /// 6. Run a simple sphere-sphere narrow-phase to detect contacts.
    /// 7. Resolve contacts with impulse-based collision response.
    /// 8. Update sleep state.
    /// 9. Clear accumulated forces.
    pub fn step(&mut self, dt: f64) {
        let g = self.config.gravity;

        // --- Integrate velocities and positions ---
        for slot in &mut self.slots {
            let body = match slot.body.as_mut() {
                Some(b) => b,
                None => continue,
            };
            if body.is_static || body.is_kinematic || body.is_sleeping {
                // Still clear forces for static/sleeping bodies
                body.accumulated_force = [0.0; 3];
                body.accumulated_torque = [0.0; 3];
                continue;
            }

            let inv_m = body.inv_mass();

            // Gravity
            let total_force = [
                body.accumulated_force[0] + g[0] * body.mass,
                body.accumulated_force[1] + g[1] * body.mass,
                body.accumulated_force[2] + g[2] * body.mass,
            ];

            // Velocity update (symplectic Euler)
            body.velocity[0] += total_force[0] * inv_m * dt;
            body.velocity[1] += total_force[1] * inv_m * dt;
            body.velocity[2] += total_force[2] * inv_m * dt;

            // Linear damping
            let lin_damp = (1.0 - body.linear_damping * dt).max(0.0);
            body.velocity[0] *= lin_damp;
            body.velocity[1] *= lin_damp;
            body.velocity[2] *= lin_damp;

            // Position update
            body.position[0] += body.velocity[0] * dt;
            body.position[1] += body.velocity[1] * dt;
            body.position[2] += body.velocity[2] * dt;

            // Angular velocity update from torque
            let approx_i = approximate_inertia(body);
            let inv_i = if approx_i > 1e-12 {
                1.0 / approx_i
            } else {
                0.0
            };
            body.angular_velocity[0] += body.accumulated_torque[0] * inv_i * dt;
            body.angular_velocity[1] += body.accumulated_torque[1] * inv_i * dt;
            body.angular_velocity[2] += body.accumulated_torque[2] * inv_i * dt;

            // Angular damping
            let ang_damp = (1.0 - body.angular_damping * dt).max(0.0);
            body.angular_velocity[0] *= ang_damp;
            body.angular_velocity[1] *= ang_damp;
            body.angular_velocity[2] *= ang_damp;

            // Integrate orientation via quaternion derivative
            body.orientation = integrate_orientation(body.orientation, body.angular_velocity, dt);

            // Clear accumulated forces
            body.accumulated_force = [0.0; 3];
            body.accumulated_torque = [0.0; 3];
        }

        // --- Narrow-phase collision detection + response ---
        self.contacts.clear();
        self.detect_and_resolve_contacts(dt);

        // --- Constraint resolution ---
        self.resolve_constraints();

        // --- Sleep detection ---
        if self.config.sleep_enabled {
            self.update_sleep(dt);
        }

        self.time += dt;
    }

    /// Reset the world: remove all bodies, clear contacts, constraints, reset time.
    pub fn reset(&mut self) {
        self.slots.clear();
        self.free_list.clear();
        self.contacts.clear();
        self.constraints.clear();
        self.time = 0.0;
    }

    // -----------------------------------------------------------------------
    // Substep simulation control
    // -----------------------------------------------------------------------

    /// Advance the simulation by `dt` seconds using `substeps` equal sub-steps.
    ///
    /// This is useful when CCD is needed for fast-moving objects, or when
    /// a larger integration step must be subdivided for stability.
    ///
    /// # Example
    /// ```ignore
    /// world.step_substeps(1.0 / 60.0, 4); // 4 sub-steps of 1/240 s each
    /// ```
    pub fn step_substeps(&mut self, dt: f64, substeps: u32) {
        let substeps = substeps.max(1);
        let sub_dt = dt / substeps as f64;
        for _ in 0..substeps {
            self.step(sub_dt);
        }
    }

    // -----------------------------------------------------------------------
    // Query functions
    // -----------------------------------------------------------------------

    /// Return the handles of all bodies whose centre-of-mass lies within `aabb`.
    ///
    /// `aabb` is given as `[xmin, ymin, zmin, xmax, ymax, zmax]`.
    pub fn bodies_in_aabb(&self, aabb: [f64; 6]) -> Vec<u32> {
        let [xmin, ymin, zmin, xmax, ymax, zmax] = aabb;
        self.slots
            .iter()
            .enumerate()
            .filter_map(|(i, slot)| {
                let body = slot.body.as_ref()?;
                let p = body.position;
                if p[0] >= xmin
                    && p[0] <= xmax
                    && p[1] >= ymin
                    && p[1] <= ymax
                    && p[2] >= zmin
                    && p[2] <= zmax
                {
                    Some(i as u32)
                } else {
                    None
                }
            })
            .collect()
    }

    /// Cast a ray from `origin` in direction `dir` and return the handle and
    /// hit distance of the nearest body whose sphere collider is struck.
    ///
    /// Returns `None` if no body is hit within `max_dist`.
    ///
    /// Direction `dir` need not be normalized; it is normalised internally.
    pub fn raycast(&self, origin: [f64; 3], dir: [f64; 3], max_dist: f64) -> Option<(u32, f64)> {
        let len = (dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]).sqrt();
        if len < 1e-12 {
            return None;
        }
        let d = [dir[0] / len, dir[1] / len, dir[2] / len];

        let mut best_handle: Option<u32> = None;
        let mut best_t = max_dist;

        for (i, slot) in self.slots.iter().enumerate() {
            let body = match slot.body.as_ref() {
                Some(b) => b,
                None => continue,
            };
            let radius = first_sphere_radius(&body.shapes);
            if radius <= 0.0 {
                continue;
            }
            // Ray-sphere intersection
            let oc = [
                origin[0] - body.position[0],
                origin[1] - body.position[1],
                origin[2] - body.position[2],
            ];
            let b_coeff = oc[0] * d[0] + oc[1] * d[1] + oc[2] * d[2];
            let c_coeff = oc[0] * oc[0] + oc[1] * oc[1] + oc[2] * oc[2] - radius * radius;
            let discriminant = b_coeff * b_coeff - c_coeff;
            if discriminant < 0.0 {
                continue;
            }
            let sqrt_disc = discriminant.sqrt();
            let t = -b_coeff - sqrt_disc;
            let t = if t >= 0.0 { t } else { -b_coeff + sqrt_disc };
            if t >= 0.0 && t < best_t {
                best_t = t;
                best_handle = Some(i as u32);
            }
        }

        best_handle.map(|h| (h, best_t))
    }

    /// Return the total number of contacts from the most recent step.
    pub fn contact_count(&self) -> usize {
        self.contacts.len()
    }
}

// ---------------------------------------------------------------------------
// Private helpers and Rust-internal constructors (NOT exposed to Python)
// ---------------------------------------------------------------------------

impl PyPhysicsWorld {
    /// Create a world from an explicit `PySimConfig` (Rust-only constructor).
    ///
    /// Rust tests and the `py_new` Python constructor both delegate here.
    pub fn new(config: PySimConfig) -> Self {
        Self::new_with_config(config)
    }

    /// Get the current simulation configuration (Rust-only access).
    pub fn config(&self) -> &PySimConfig {
        &self.config
    }

    /// Contacts from the most recent simulation step (Rust-only access).
    pub fn get_contacts(&self) -> Vec<PyContactResult> {
        self.contacts.clone()
    }

    /// Add a body from a legacy `PyRigidBodyDesc` (Rust-only, for backward compat).
    pub fn add_body_legacy(&mut self, desc: &PyRigidBodyDesc) -> u32 {
        let config = if desc.is_static {
            PyRigidBodyConfig::static_body([desc.position.x, desc.position.y, desc.position.z])
        } else {
            PyRigidBodyConfig::dynamic(
                desc.mass,
                [desc.position.x, desc.position.y, desc.position.z],
            )
        };
        self.add_rigid_body(config)
    }

    /// Get the position of a body by slot index as `PyVec3`, or `None` if absent.
    pub fn get_body_position(&self, idx: usize) -> Option<PyVec3> {
        let body = self.slots.get(idx)?.body.as_ref()?;
        Some(PyVec3::new(
            body.position[0],
            body.position[1],
            body.position[2],
        ))
    }

    /// Create a world from an explicit `PySimConfig`.
    ///
    /// This is a Rust-only constructor; Python uses the `#[new]` constructor
    /// which takes primitive gravity components.
    pub fn new_with_config(config: PySimConfig) -> Self {
        Self {
            slots: Vec::new(),
            free_list: Vec::new(),
            config,
            time: 0.0,
            contacts: Vec::new(),
            constraints: Vec::new(),
        }
    }

    /// Add a body described by `config` and return its integer handle.
    pub fn add_rigid_body(&mut self, config: PyRigidBodyConfig) -> u32 {
        let body = InternalBody::from_config(&config);
        if let Some(idx) = self.free_list.pop() {
            let slot = &mut self.slots[idx as usize];
            slot.body = Some(body);
            slot.generation = slot.generation.wrapping_add(1);
            idx
        } else {
            let idx = self.slots.len() as u32;
            self.slots.push(Slot {
                body: Some(body),
                generation: 0,
            });
            idx
        }
    }

    fn get_body(&self, handle: u32) -> Option<&InternalBody> {
        self.slots.get(handle as usize)?.body.as_ref()
    }

    fn get_body_mut(&mut self, handle: u32) -> Option<&mut InternalBody> {
        self.slots.get_mut(handle as usize)?.body.as_mut()
    }

    /// Very simple sphere-sphere broad/narrow phase.
    ///
    /// For each pair of bodies that both have sphere shapes, compute penetration
    /// depth. If penetrating, record a contact and apply an impulse-based
    /// resolution (one iteration).
    fn detect_and_resolve_contacts(&mut self, dt: f64) {
        let n = self.slots.len();
        for i in 0..n {
            for j in (i + 1)..n {
                // Collect needed data without borrow issues
                let (pos_a, vel_a, mass_a, static_a, sleeping_a, friction_a, rest_a, radius_a) = {
                    let body = match self.slots[i].body.as_ref() {
                        Some(b) => b,
                        None => continue,
                    };
                    let radius = first_sphere_radius(&body.shapes);
                    if radius <= 0.0 {
                        continue;
                    }
                    (
                        body.position,
                        body.velocity,
                        body.mass,
                        body.is_static,
                        body.is_sleeping,
                        body.friction,
                        body.restitution,
                        radius,
                    )
                };

                let (pos_b, vel_b, mass_b, static_b, sleeping_b, friction_b, rest_b, radius_b) = {
                    let body = match self.slots[j].body.as_ref() {
                        Some(b) => b,
                        None => continue,
                    };
                    let radius = first_sphere_radius(&body.shapes);
                    if radius <= 0.0 {
                        continue;
                    }
                    (
                        body.position,
                        body.velocity,
                        body.mass,
                        body.is_static,
                        body.is_sleeping,
                        body.friction,
                        body.restitution,
                        radius,
                    )
                };

                // Distance between centers
                let diff = [
                    pos_b[0] - pos_a[0],
                    pos_b[1] - pos_a[1],
                    pos_b[2] - pos_a[2],
                ];
                let dist_sq = diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2];
                let min_dist = radius_a + radius_b;
                if dist_sq >= min_dist * min_dist {
                    continue;
                }

                let dist = dist_sq.sqrt();
                let depth = min_dist - dist;
                let normal = if dist > 1e-12 {
                    [diff[0] / dist, diff[1] / dist, diff[2] / dist]
                } else {
                    [0.0, 1.0, 0.0]
                };

                // Contact point (midpoint on surface of A towards B)
                let cp = [
                    pos_a[0] + normal[0] * radius_a,
                    pos_a[1] + normal[1] * radius_a,
                    pos_a[2] + normal[2] * radius_a,
                ];

                // Combined coefficients
                let combined_friction = (friction_a + friction_b) * 0.5;
                let combined_rest = (rest_a + rest_b) * 0.5;

                // Relative velocity along normal
                let rel_vel = [
                    vel_a[0] - vel_b[0],
                    vel_a[1] - vel_b[1],
                    vel_a[2] - vel_b[2],
                ];
                let rel_vel_normal =
                    rel_vel[0] * normal[0] + rel_vel[1] * normal[1] + rel_vel[2] * normal[2];

                // Impulse magnitude (only if approaching)
                let impulse_mag = if rel_vel_normal < 0.0 {
                    let inv_m_a = if static_a || mass_a <= 0.0 {
                        0.0
                    } else {
                        1.0 / mass_a
                    };
                    let inv_m_b = if static_b || mass_b <= 0.0 {
                        0.0
                    } else {
                        1.0 / mass_b
                    };
                    let denom = inv_m_a + inv_m_b;
                    if denom > 1e-12 {
                        -(1.0 + combined_rest) * rel_vel_normal / denom
                    } else {
                        0.0
                    }
                } else {
                    0.0
                };

                // Baumgarte positional correction
                let correction_mag = (depth * self.config.baumgarte_factor) / dt.max(1e-6);

                let contact = PyContactResult {
                    body_a: i as u32,
                    body_b: j as u32,
                    contact_point: cp,
                    normal,
                    depth,
                    friction: combined_friction,
                    restitution: combined_rest,
                    impulse: impulse_mag,
                };
                self.contacts.push(contact);

                // Apply impulse to body A (not static/sleeping)
                if !static_a && !sleeping_a && mass_a > 0.0 {
                    let inv_m_a = 1.0 / mass_a;
                    if let Some(ba) = self.slots[i].body.as_mut() {
                        ba.velocity[0] += impulse_mag * normal[0] * inv_m_a;
                        ba.velocity[1] += impulse_mag * normal[1] * inv_m_a;
                        ba.velocity[2] += impulse_mag * normal[2] * inv_m_a;
                        // Positional correction
                        ba.position[0] -= normal[0] * correction_mag * inv_m_a * dt;
                        ba.position[1] -= normal[1] * correction_mag * inv_m_a * dt;
                        ba.position[2] -= normal[2] * correction_mag * inv_m_a * dt;
                    }
                }
                // Apply impulse to body B (opposite direction)
                if !static_b && !sleeping_b && mass_b > 0.0 {
                    let inv_m_b = 1.0 / mass_b;
                    if let Some(bb) = self.slots[j].body.as_mut() {
                        bb.velocity[0] -= impulse_mag * normal[0] * inv_m_b;
                        bb.velocity[1] -= impulse_mag * normal[1] * inv_m_b;
                        bb.velocity[2] -= impulse_mag * normal[2] * inv_m_b;
                        // Positional correction
                        bb.position[0] += normal[0] * correction_mag * inv_m_b * dt;
                        bb.position[1] += normal[1] * correction_mag * inv_m_b * dt;
                        bb.position[2] += normal[2] * correction_mag * inv_m_b * dt;
                    }
                }
            }
        }
    }

    fn update_sleep(&mut self, dt: f64) {
        let lin_thresh = self.config.linear_sleep_threshold;
        let ang_thresh = self.config.angular_sleep_threshold;
        let time_thresh = self.config.time_before_sleep;

        for slot in &mut self.slots {
            let body = match slot.body.as_mut() {
                Some(b) => b,
                None => continue,
            };
            if !body.can_sleep || body.is_static || body.is_kinematic {
                continue;
            }
            if body.is_sleeping {
                continue;
            }

            let lin_speed = (body.velocity[0] * body.velocity[0]
                + body.velocity[1] * body.velocity[1]
                + body.velocity[2] * body.velocity[2])
                .sqrt();
            let ang_speed = (body.angular_velocity[0] * body.angular_velocity[0]
                + body.angular_velocity[1] * body.angular_velocity[1]
                + body.angular_velocity[2] * body.angular_velocity[2])
                .sqrt();

            if lin_speed < lin_thresh && ang_speed < ang_thresh {
                body.sleep_timer += dt;
                if body.sleep_timer >= time_thresh {
                    body.is_sleeping = true;
                    body.velocity = [0.0; 3];
                    body.angular_velocity = [0.0; 3];
                }
            } else {
                body.sleep_timer = 0.0;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Pure helper functions (no `self` borrow issues)
// ---------------------------------------------------------------------------

fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn normalize_quat(q: [f64; 4]) -> [f64; 4] {
    let norm = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]).sqrt();
    if norm < 1e-12 {
        [0.0, 0.0, 0.0, 1.0]
    } else {
        [q[0] / norm, q[1] / norm, q[2] / norm, q[3] / norm]
    }
}

/// Integrate an orientation quaternion by angular velocity over dt.
///
/// Uses the quaternion derivative: dq/dt = 0.5 * \[wx,wy,wz,0\] * q
fn integrate_orientation(q: [f64; 4], omega: [f64; 3], dt: f64) -> [f64; 4] {
    let (qx, qy, qz, qw) = (q[0], q[1], q[2], q[3]);
    let (wx, wy, wz) = (omega[0], omega[1], omega[2]);
    let half_dt = 0.5 * dt;
    let dqx = half_dt * (qw * wx + qy * wz - qz * wy);
    let dqy = half_dt * (qw * wy - qx * wz + qz * wx);
    let dqz = half_dt * (qw * wz + qx * wy - qy * wx);
    let dqw = half_dt * (-qx * wx - qy * wy - qz * wz);
    normalize_quat([qx + dqx, qy + dqy, qz + dqz, qw + dqw])
}

/// Compute an approximate scalar moment of inertia for a body.
///
/// Uses the first sphere shape if present; otherwise falls back to a point-mass
/// approximation.
fn approximate_inertia(body: &InternalBody) -> f64 {
    if body.mass <= 0.0 {
        return 0.0;
    }
    for shape in &body.shapes {
        if let PyColliderShape::Sphere { radius } = shape {
            // Solid sphere: I = 2/5 * m * r^2
            return 0.4 * body.mass * radius * radius;
        }
        if let PyColliderShape::Box { half_extents } = shape {
            // Box: I_xx = 1/12 * m * (hy^2 + hz^2), use average
            let hx = half_extents[0];
            let hy = half_extents[1];
            let hz = half_extents[2];
            let ixx = (1.0 / 3.0) * body.mass * (hy * hy + hz * hz);
            let iyy = (1.0 / 3.0) * body.mass * (hx * hx + hz * hz);
            let izz = (1.0 / 3.0) * body.mass * (hx * hx + hy * hy);
            return (ixx + iyy + izz) / 3.0;
        }
    }
    // Point mass fallback: use a unit sphere assumption
    0.4 * body.mass
}

/// Return the radius of the first sphere shape in the list, or 0.0.
fn first_sphere_radius(shapes: &[PyColliderShape]) -> f64 {
    for shape in shapes {
        if let PyColliderShape::Sphere { radius } = shape {
            return *radius;
        }
    }
    0.0
}

/// Register all `world_api` classes and functions into a Python sub-module named `"world"`.
pub fn register_world_module(
    parent: &pyo3::Bound<'_, pyo3::types::PyModule>,
) -> pyo3::PyResult<()> {
    use pyo3::types::PyModuleMethods;
    let child = pyo3::types::PyModule::new(parent.py(), "world")?;

    // Core world class
    child.add_class::<PyPhysicsWorld>()?;

    // Constraint types
    child.add_class::<constraints::ConstraintType>()?;
    child.add_class::<constraints::PyConstraint>()?;

    // Extension types
    child.add_class::<extensions::ContactPair>()?;
    child.add_class::<extensions::InertiaTensor>()?;
    child.add_class::<extensions::PyRigidBody>()?;

    // Geometry types
    child.add_class::<geometry::PyAabb>()?;
    child.add_class::<geometry::PySphere>()?;
    child.add_class::<geometry::PyConvexHull>()?;

    // Material types
    child.add_class::<materials::MaterialClass>()?;
    child.add_class::<materials::PyMaterial>()?;

    // LBM types
    child.add_class::<lbm::PyLbmConfig>()?;
    child.add_class::<lbm::PyLbmGrid>()?;

    // SPH types
    child.add_class::<sph::PySphConfig>()?;
    child.add_class::<sph::PySphSim>()?;

    // FEM types
    child.add_class::<fem::FemBarElement>()?;
    child.add_class::<fem::PyFemAssembly>()?;

    // Binding wrappers
    child.add_class::<py_bindings::PyMdBinding>()?;
    child.add_class::<py_bindings::PyLbmBinding>()?;
    child.add_class::<py_bindings::PyFemBinding>()?;
    child.add_class::<py_bindings::PySphBinding>()?;
    child.add_class::<py_bindings::MaterialProperties>()?;

    // Statistics
    child.add_class::<stats::SimStats>()?;

    // PyFunctions
    child.add_function(pyo3::wrap_pyfunction!(
        py_bindings::py_query_material,
        &child
    )?)?;
    child.add_function(pyo3::wrap_pyfunction!(
        py_bindings::py_p_wave_speed,
        &child
    )?)?;
    child.add_function(pyo3::wrap_pyfunction!(
        py_bindings::py_s_wave_speed,
        &child
    )?)?;

    child.add_function(pyo3::wrap_pyfunction!(
        math_helpers::quat_normalize,
        &child
    )?)?;
    child.add_function(pyo3::wrap_pyfunction!(math_helpers::quat_mul, &child)?)?;
    child.add_function(pyo3::wrap_pyfunction!(
        math_helpers::quat_conjugate,
        &child
    )?)?;
    child.add_function(pyo3::wrap_pyfunction!(
        math_helpers::quat_rotate_vec,
        &child
    )?)?;
    child.add_function(pyo3::wrap_pyfunction!(
        math_helpers::quat_from_axis_angle,
        &child
    )?)?;

    parent.add_submodule(&child)?;
    Ok(())
}
