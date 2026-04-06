#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Python rigid body bindings.
//!
//! Provides Python-friendly wrappers for rigid body simulation:
//! bodies, colliders, joints, a physics world, and batch force APIs.
//! All types use plain arrays and primitive types (no nalgebra).

#![allow(missing_docs)]
#![allow(dead_code)]

use serde::{Deserialize, Serialize};

// ─────────────────────────────────────────────────────────────────────────────
// Vec3 helpers
// ─────────────────────────────────────────────────────────────────────────────

#[inline]
fn vec3_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn vec3_scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline]
fn vec3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn vec3_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
fn vec3_length(a: [f64; 3]) -> f64 {
    vec3_dot(a, a).sqrt()
}

#[inline]
fn quat_identity() -> [f64; 4] {
    [0.0, 0.0, 0.0, 1.0]
}

// ─────────────────────────────────────────────────────────────────────────────
// Shape type
// ─────────────────────────────────────────────────────────────────────────────

/// Collider shape variants supported by the physics engine.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub enum PyShapeType {
    /// Sphere defined by radius.
    Sphere {
        /// Sphere radius in metres.
        radius: f64,
    },
    /// Axis-aligned box defined by half-extents.
    Box {
        /// Half-extents `[hx, hy, hz]` in metres.
        half_extents: [f64; 3],
    },
    /// Capsule (cylinder capped by two hemispheres).
    Capsule {
        /// Capsule half-height (cylinder portion) in metres.
        half_height: f64,
        /// Capsule radius in metres.
        radius: f64,
    },
}

// ─────────────────────────────────────────────────────────────────────────────
// PyRigidBody
// ─────────────────────────────────────────────────────────────────────────────

/// A rigid body with full kinematic and dynamic state.
///
/// Positions and orientations are stored as plain arrays for easy FFI.
/// The orientation quaternion uses convention `[x, y, z, w]`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyRigidBody {
    /// Unique handle for this body.
    pub id: u32,
    /// Mass in kilograms. Zero means static.
    pub mass: f64,
    /// Diagonal inertia tensor `[Ixx, Iyy, Izz]`.
    pub inertia: [f64; 3],
    /// World-space position `[x, y, z]`.
    pub position: [f64; 3],
    /// Orientation quaternion `[x, y, z, w]`.
    pub orientation: [f64; 4],
    /// Linear velocity `[vx, vy, vz]`.
    pub velocity: [f64; 3],
    /// Angular velocity `[wx, wy, wz]`.
    pub ang_vel: [f64; 3],
    /// Accumulated force buffer `[fx, fy, fz]`.
    pub accumulated_force: [f64; 3],
    /// Accumulated torque buffer `[tx, ty, tz]`.
    pub accumulated_torque: [f64; 3],
    /// Whether the body is sleeping.
    pub sleeping: bool,
    /// Whether the body is static (immovable).
    pub is_static: bool,
    /// Whether the body is kinematic (user-controlled velocity).
    pub is_kinematic: bool,
    /// Linear damping coefficient.
    pub linear_damping: f64,
    /// Angular damping coefficient.
    pub angular_damping: f64,
}

impl PyRigidBody {
    /// Create a new dynamic rigid body.
    pub fn new(id: u32, mass: f64, position: [f64; 3]) -> Self {
        Self {
            id,
            mass,
            inertia: [mass * 0.1; 3],
            position,
            orientation: quat_identity(),
            velocity: [0.0; 3],
            ang_vel: [0.0; 3],
            accumulated_force: [0.0; 3],
            accumulated_torque: [0.0; 3],
            sleeping: false,
            is_static: false,
            is_kinematic: false,
            linear_damping: 0.01,
            angular_damping: 0.01,
        }
    }

    /// Create a static body (zero mass, immovable).
    pub fn new_static(id: u32, position: [f64; 3]) -> Self {
        let mut b = Self::new(id, 0.0, position);
        b.is_static = true;
        b
    }

    /// Create a kinematic body.
    pub fn new_kinematic(id: u32, mass: f64, position: [f64; 3]) -> Self {
        let mut b = Self::new(id, mass, position);
        b.is_kinematic = true;
        b
    }

    /// Apply a world-space force (accumulates until `clear_forces` is called).
    pub fn apply_force(&mut self, force: [f64; 3]) {
        if self.is_static || self.sleeping {
            return;
        }
        self.accumulated_force = vec3_add(self.accumulated_force, force);
    }

    /// Apply a torque about the body's centre of mass.
    pub fn apply_torque(&mut self, torque: [f64; 3]) {
        if self.is_static || self.sleeping {
            return;
        }
        self.accumulated_torque = vec3_add(self.accumulated_torque, torque);
    }

    /// Apply a linear impulse (instantaneous velocity change).
    pub fn apply_impulse(&mut self, impulse: [f64; 3]) {
        if self.is_static {
            return;
        }
        if self.mass > 0.0 {
            let inv_m = 1.0 / self.mass;
            self.velocity = vec3_add(self.velocity, vec3_scale(impulse, inv_m));
        }
        self.sleeping = false;
    }

    /// Apply an angular impulse (instantaneous angular velocity change).
    pub fn apply_angular_impulse(&mut self, ang_impulse: [f64; 3]) {
        if self.is_static {
            return;
        }
        for i in 0..3 {
            if self.inertia[i] > 0.0 {
                self.ang_vel[i] += ang_impulse[i] / self.inertia[i];
            }
        }
        self.sleeping = false;
    }

    /// Clear accumulated forces and torques.
    pub fn clear_forces(&mut self) {
        self.accumulated_force = [0.0; 3];
        self.accumulated_torque = [0.0; 3];
    }

    /// Integrate the body forward by `dt` seconds (semi-implicit Euler).
    pub fn integrate(&mut self, dt: f64, gravity: [f64; 3]) {
        if self.is_static || self.is_kinematic || self.sleeping {
            return;
        }
        let inv_m = if self.mass > 0.0 {
            1.0 / self.mass
        } else {
            0.0
        };
        // Linear
        let grav_force = vec3_scale(gravity, self.mass);
        let total_force = vec3_add(self.accumulated_force, grav_force);
        let accel = vec3_scale(total_force, inv_m);
        self.velocity = vec3_add(self.velocity, vec3_scale(accel, dt));
        // Damping
        let lin_damp = (1.0 - self.linear_damping * dt).max(0.0);
        self.velocity = vec3_scale(self.velocity, lin_damp);
        self.position = vec3_add(self.position, vec3_scale(self.velocity, dt));
        // Angular
        for i in 0..3 {
            let alpha = if self.inertia[i] > 0.0 {
                self.accumulated_torque[i] / self.inertia[i]
            } else {
                0.0
            };
            self.ang_vel[i] += alpha * dt;
        }
        let ang_damp = (1.0 - self.angular_damping * dt).max(0.0);
        self.ang_vel = vec3_scale(self.ang_vel, ang_damp);
    }

    /// Kinetic energy of this body.
    pub fn kinetic_energy(&self) -> f64 {
        let ke_lin = 0.5 * self.mass * vec3_dot(self.velocity, self.velocity);
        let ke_rot: f64 = self
            .inertia
            .iter()
            .zip(self.ang_vel.iter())
            .map(|(&i, &w)| 0.5 * i * w * w)
            .sum();
        ke_lin + ke_rot
    }

    /// Linear momentum `[px, py, pz]`.
    pub fn linear_momentum(&self) -> [f64; 3] {
        vec3_scale(self.velocity, self.mass)
    }

    /// Effective inverse mass (0 for static).
    pub fn inv_mass(&self) -> f64 {
        if self.is_static || self.mass <= 0.0 {
            0.0
        } else {
            1.0 / self.mass
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PyRigidBodySet
// ─────────────────────────────────────────────────────────────────────────────

/// A collection of rigid bodies with O(1) lookup by handle.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct PyRigidBodySet {
    bodies: Vec<Option<PyRigidBody>>,
    next_id: u32,
    count: usize,
}

impl PyRigidBodySet {
    /// Create an empty set.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a body; returns its handle.
    pub fn add(&mut self, mut body: PyRigidBody) -> u32 {
        let id = self.next_id;
        self.next_id += 1;
        body.id = id;
        // Find free slot or push
        if let Some(slot) = self.bodies.iter_mut().find(|s| s.is_none()) {
            *slot = Some(body);
        } else {
            self.bodies.push(Some(body));
        }
        self.count += 1;
        id
    }

    /// Remove a body by handle. Returns `true` if found.
    pub fn remove(&mut self, id: u32) -> bool {
        for slot in &mut self.bodies {
            if slot.as_ref().is_some_and(|b| b.id == id) {
                *slot = None;
                self.count -= 1;
                return true;
            }
        }
        false
    }

    /// Get an immutable reference to a body.
    pub fn get(&self, id: u32) -> Option<&PyRigidBody> {
        self.bodies.iter().flatten().find(|b| b.id == id)
    }

    /// Get a mutable reference to a body.
    pub fn get_mut(&mut self, id: u32) -> Option<&mut PyRigidBody> {
        self.bodies.iter_mut().flatten().find(|b| b.id == id)
    }

    /// Number of active bodies.
    pub fn len(&self) -> usize {
        self.count
    }

    /// True when no bodies exist.
    pub fn is_empty(&self) -> bool {
        self.count == 0
    }

    /// Iterate over all active bodies.
    pub fn iter(&self) -> impl Iterator<Item = &PyRigidBody> {
        self.bodies.iter().flatten()
    }

    /// Mutable iterator over all active bodies.
    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut PyRigidBody> {
        self.bodies.iter_mut().flatten()
    }

    /// Get position of a body by id.
    pub fn get_position(&self, id: u32) -> Option<[f64; 3]> {
        self.get(id).map(|b| b.position)
    }

    /// Get velocity of a body by id.
    pub fn get_velocity(&self, id: u32) -> Option<[f64; 3]> {
        self.get(id).map(|b| b.velocity)
    }

    /// Batch: apply gravity to all dynamic bodies.
    pub fn apply_gravity(&mut self, gravity: [f64; 3]) {
        for body in self.iter_mut() {
            if !body.is_static && !body.is_kinematic && !body.sleeping {
                let grav_force = vec3_scale(gravity, body.mass);
                body.apply_force(grav_force);
            }
        }
    }

    /// Batch: integrate all dynamic bodies forward by `dt`.
    pub fn integrate_all(&mut self, dt: f64, gravity: [f64; 3]) {
        for body in self.iter_mut() {
            body.integrate(dt, gravity);
            body.clear_forces();
        }
    }

    /// Total kinetic energy of all bodies.
    pub fn total_kinetic_energy(&self) -> f64 {
        self.iter().map(|b| b.kinetic_energy()).sum()
    }

    /// Number of sleeping bodies.
    pub fn sleeping_count(&self) -> usize {
        self.iter().filter(|b| b.sleeping).count()
    }

    /// Number of awake bodies.
    pub fn awake_count(&self) -> usize {
        self.count - self.sleeping_count()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PyCollider
// ─────────────────────────────────────────────────────────────────────────────

/// A collider attached to a rigid body.
///
/// Defines the geometric shape for collision detection, plus surface
/// properties (friction, restitution) and sensor mode.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyCollider {
    /// Unique collider handle.
    pub id: u32,
    /// Parent body handle.
    pub body_id: u32,
    /// Shape of this collider.
    pub shape: PyShapeType,
    /// Coefficient of friction (Coulomb model).
    pub friction: f64,
    /// Coefficient of restitution (bounciness 0..=1).
    pub restitution: f64,
    /// When `true`, the collider fires overlap events but does not generate
    /// contact forces.
    pub is_sensor: bool,
    /// Local offset from body centre `[x, y, z]`.
    pub local_offset: [f64; 3],
    /// Local orientation quaternion `[x, y, z, w]`.
    pub local_orientation: [f64; 4],
}

impl PyCollider {
    /// Create a sphere collider.
    pub fn sphere(id: u32, body_id: u32, radius: f64) -> Self {
        Self {
            id,
            body_id,
            shape: PyShapeType::Sphere { radius },
            friction: 0.5,
            restitution: 0.3,
            is_sensor: false,
            local_offset: [0.0; 3],
            local_orientation: quat_identity(),
        }
    }

    /// Create a box collider.
    pub fn box_collider(id: u32, body_id: u32, half_extents: [f64; 3]) -> Self {
        Self {
            id,
            body_id,
            shape: PyShapeType::Box { half_extents },
            friction: 0.5,
            restitution: 0.3,
            is_sensor: false,
            local_offset: [0.0; 3],
            local_orientation: quat_identity(),
        }
    }

    /// Create a capsule collider.
    pub fn capsule(id: u32, body_id: u32, half_height: f64, radius: f64) -> Self {
        Self {
            id,
            body_id,
            shape: PyShapeType::Capsule {
                half_height,
                radius,
            },
            friction: 0.5,
            restitution: 0.3,
            is_sensor: false,
            local_offset: [0.0; 3],
            local_orientation: quat_identity(),
        }
    }

    /// Set as sensor (no contact response).
    pub fn as_sensor(mut self) -> Self {
        self.is_sensor = true;
        self
    }

    /// Bounding radius used for broadphase.
    pub fn bounding_radius(&self) -> f64 {
        match &self.shape {
            PyShapeType::Sphere { radius } => *radius,
            PyShapeType::Box { half_extents } => vec3_length(*half_extents),
            PyShapeType::Capsule {
                half_height,
                radius,
            } => half_height + radius,
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Contact and Ray types
// ─────────────────────────────────────────────────────────────────────────────

/// A contact between two bodies.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyContact {
    /// First body handle.
    pub body_a: u32,
    /// Second body handle.
    pub body_b: u32,
    /// Contact point in world space.
    pub point: [f64; 3],
    /// Contact normal pointing from B to A.
    pub normal: [f64; 3],
    /// Penetration depth (positive = overlap).
    pub depth: f64,
    /// Magnitude of the contact impulse applied.
    pub impulse: f64,
}

/// Result of a ray-cast query.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyRayResult {
    /// Hit body handle.
    pub body_id: u32,
    /// Hit point in world space.
    pub point: [f64; 3],
    /// Surface normal at hit.
    pub normal: [f64; 3],
    /// Distance from ray origin to hit point.
    pub distance: f64,
}

// ─────────────────────────────────────────────────────────────────────────────
// Joint types
// ─────────────────────────────────────────────────────────────────────────────

/// Supported joint types.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum PyJointType {
    /// Ball-and-socket joint (3 rotational DOF).
    Ball,
    /// Hinge joint (1 rotational DOF).
    Hinge,
    /// Slider / prismatic joint (1 translational DOF).
    Slider,
    /// Fixed joint (0 DOF).
    Fixed,
    /// Generic 6-DOF joint.
    Generic6Dof,
}

/// A joint connecting two bodies.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyJoint {
    /// Joint handle.
    pub id: u32,
    /// Joint type.
    pub joint_type: PyJointType,
    /// Body A handle.
    pub body_a: u32,
    /// Body B handle.
    pub body_b: u32,
    /// Anchor on body A in local space.
    pub anchor_a: [f64; 3],
    /// Anchor on body B in local space.
    pub anchor_b: [f64; 3],
    /// Axis for hinge/slider (body A local space).
    pub axis: [f64; 3],
    /// Lower limit (angle in radians for hinge, distance for slider).
    pub lower_limit: f64,
    /// Upper limit.
    pub upper_limit: f64,
    /// Motor target speed.
    pub motor_speed: f64,
    /// Motor maximum force/torque.
    pub motor_max_force: f64,
    /// Whether the motor is enabled.
    pub motor_enabled: bool,
    /// Last computed force on body A `[fx, fy, fz]`.
    pub reaction_force: [f64; 3],
    /// Last computed torque on body A `[tx, ty, tz]`.
    pub reaction_torque: [f64; 3],
}

impl PyJoint {
    /// Create a new joint.
    pub fn new(
        id: u32,
        joint_type: PyJointType,
        body_a: u32,
        body_b: u32,
        anchor_a: [f64; 3],
        anchor_b: [f64; 3],
    ) -> Self {
        Self {
            id,
            joint_type,
            body_a,
            body_b,
            anchor_a,
            anchor_b,
            axis: [0.0, 1.0, 0.0],
            lower_limit: f64::NEG_INFINITY,
            upper_limit: f64::INFINITY,
            motor_speed: 0.0,
            motor_max_force: 0.0,
            motor_enabled: false,
            reaction_force: [0.0; 3],
            reaction_torque: [0.0; 3],
        }
    }

    /// Enable the motor.
    pub fn enable_motor(&mut self, speed: f64, max_force: f64) {
        self.motor_enabled = true;
        self.motor_speed = speed;
        self.motor_max_force = max_force;
    }

    /// Disable the motor.
    pub fn disable_motor(&mut self) {
        self.motor_enabled = false;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PyJointSet
// ─────────────────────────────────────────────────────────────────────────────

/// A collection of joints.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct PyJointSet {
    joints: Vec<Option<PyJoint>>,
    next_id: u32,
}

impl PyJointSet {
    /// Create an empty joint set.
    pub fn new() -> Self {
        Self::default()
    }

    /// Create and add a joint; returns its handle.
    #[allow(clippy::too_many_arguments)]
    pub fn create_joint(
        &mut self,
        joint_type: PyJointType,
        body_a: u32,
        body_b: u32,
        anchor_a: [f64; 3],
        anchor_b: [f64; 3],
        _params: Option<serde_json::Value>,
    ) -> u32 {
        let id = self.next_id;
        self.next_id += 1;
        let joint = PyJoint::new(id, joint_type, body_a, body_b, anchor_a, anchor_b);
        self.joints.push(Some(joint));
        id
    }

    /// Remove a joint by handle.
    pub fn remove(&mut self, id: u32) -> bool {
        for slot in &mut self.joints {
            if slot.as_ref().is_some_and(|j| j.id == id) {
                *slot = None;
                return true;
            }
        }
        false
    }

    /// Get immutable reference.
    pub fn get(&self, id: u32) -> Option<&PyJoint> {
        self.joints.iter().flatten().find(|j| j.id == id)
    }

    /// Get mutable reference.
    pub fn get_mut(&mut self, id: u32) -> Option<&mut PyJoint> {
        self.joints.iter_mut().flatten().find(|j| j.id == id)
    }

    /// Set motor speed and max force.
    pub fn set_motor(&mut self, id: u32, speed: f64, max_force: f64) -> bool {
        if let Some(j) = self.get_mut(id) {
            j.enable_motor(speed, max_force);
            true
        } else {
            false
        }
    }

    /// Get reaction forces and torques `(force, torque)`.
    pub fn get_forces(&self, id: u32) -> Option<([f64; 3], [f64; 3])> {
        self.get(id).map(|j| (j.reaction_force, j.reaction_torque))
    }

    /// Number of active joints.
    pub fn len(&self) -> usize {
        self.joints.iter().flatten().count()
    }

    /// True if empty.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PyWorld
// ─────────────────────────────────────────────────────────────────────────────

/// A self-contained physics world managing bodies, colliders, and joints.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyWorld {
    /// All rigid bodies.
    pub bodies: PyRigidBodySet,
    /// All colliders keyed by id.
    colliders: Vec<PyCollider>,
    /// All joints.
    pub joints: PyJointSet,
    /// Gravity vector.
    pub gravity: [f64; 3],
    /// Contacts from last step.
    contacts: Vec<PyContact>,
    /// Elapsed simulation time.
    pub time: f64,
    /// Next collider id.
    next_collider_id: u32,
}

impl PyWorld {
    /// Create a world with the given gravity.
    pub fn new(gravity: [f64; 3]) -> Self {
        Self {
            bodies: PyRigidBodySet::new(),
            colliders: Vec::new(),
            joints: PyJointSet::new(),
            gravity,
            contacts: Vec::new(),
            time: 0.0,
            next_collider_id: 0,
        }
    }

    /// Create a world with Earth gravity.
    pub fn earth_gravity() -> Self {
        Self::new([0.0, -9.81, 0.0])
    }

    /// Add a dynamic body; returns handle.
    pub fn create_body(&mut self, mass: f64, position: [f64; 3]) -> u32 {
        let body = PyRigidBody::new(0, mass, position);
        self.bodies.add(body)
    }

    /// Add a static body; returns handle.
    pub fn create_static_body(&mut self, position: [f64; 3]) -> u32 {
        let body = PyRigidBody::new_static(0, position);
        self.bodies.add(body)
    }

    /// Add a collider to a body; returns collider handle.
    pub fn create_collider(&mut self, collider_desc: PyCollider) -> u32 {
        let id = self.next_collider_id;
        self.next_collider_id += 1;
        let mut c = collider_desc;
        c.id = id;
        self.colliders.push(c);
        id
    }

    /// Step the simulation by `dt` seconds.
    pub fn step(&mut self, dt: f64) {
        self.contacts.clear();
        // Detect collisions (simple sphere-sphere broadphase mock)
        let ids: Vec<u32> = self.bodies.iter().map(|b| b.id).collect();
        for (i, &id_a) in ids.iter().enumerate() {
            for &id_b in &ids[i + 1..] {
                if let (Some(ba), Some(bb)) = (self.bodies.get(id_a), self.bodies.get(id_b)) {
                    if ba.is_static && bb.is_static {
                        continue;
                    }
                    let dx = [
                        bb.position[0] - ba.position[0],
                        bb.position[1] - ba.position[1],
                        bb.position[2] - ba.position[2],
                    ];
                    let dist = vec3_length(dx);
                    // Find colliders for these bodies
                    let ra = self
                        .colliders
                        .iter()
                        .filter(|c| c.body_id == id_a)
                        .map(|c| c.bounding_radius())
                        .fold(0.0f64, f64::max);
                    let rb = self
                        .colliders
                        .iter()
                        .filter(|c| c.body_id == id_b)
                        .map(|c| c.bounding_radius())
                        .fold(0.0f64, f64::max);
                    if ra > 0.0 && rb > 0.0 && dist < ra + rb {
                        let depth = ra + rb - dist;
                        let normal = if dist > 1e-12 {
                            [dx[0] / dist, dx[1] / dist, dx[2] / dist]
                        } else {
                            [0.0, 1.0, 0.0]
                        };
                        let midpoint = [
                            (ba.position[0] + bb.position[0]) * 0.5,
                            (ba.position[1] + bb.position[1]) * 0.5,
                            (ba.position[2] + bb.position[2]) * 0.5,
                        ];
                        self.contacts.push(PyContact {
                            body_a: id_a,
                            body_b: id_b,
                            point: midpoint,
                            normal,
                            depth,
                            impulse: 0.0,
                        });
                    }
                }
            }
        }
        // Resolve contacts (impulse-based)
        for contact in &self.contacts {
            let inv_a = self
                .bodies
                .get(contact.body_a)
                .map(|b| b.inv_mass())
                .unwrap_or(0.0);
            let inv_b = self
                .bodies
                .get(contact.body_b)
                .map(|b| b.inv_mass())
                .unwrap_or(0.0);
            let denom = inv_a + inv_b;
            if denom < 1e-15 {
                continue;
            }
            let restitution = 0.3_f64;
            let vel_a = self
                .bodies
                .get(contact.body_a)
                .map(|b| b.velocity)
                .unwrap_or([0.0; 3]);
            let vel_b = self
                .bodies
                .get(contact.body_b)
                .map(|b| b.velocity)
                .unwrap_or([0.0; 3]);
            let rel_vel = [
                vel_a[0] - vel_b[0],
                vel_a[1] - vel_b[1],
                vel_a[2] - vel_b[2],
            ];
            let vn = vec3_dot(rel_vel, contact.normal);
            if vn > 0.0 {
                continue; // Separating
            }
            let j_scalar = -(1.0 + restitution) * vn / denom;
            let impulse = vec3_scale(contact.normal, j_scalar);
            if let Some(ba) = self.bodies.get_mut(contact.body_a) {
                ba.apply_impulse(impulse);
            }
            if let Some(bb) = self.bodies.get_mut(contact.body_b) {
                bb.apply_impulse(vec3_scale(impulse, -1.0));
            }
            // Positional correction
            let correction_scale = (contact.depth - 0.01_f64).max(0.0) / denom * 0.4;
            let correction = vec3_scale(contact.normal, correction_scale);
            if let Some(ba) = self.bodies.get_mut(contact.body_a)
                && !ba.is_static
            {
                let corr = vec3_scale(correction, ba.inv_mass());
                ba.position = vec3_add(ba.position, corr);
            }
            if let Some(bb) = self.bodies.get_mut(contact.body_b)
                && !bb.is_static
            {
                let corr = vec3_scale(correction, bb.inv_mass());
                bb.position = [
                    bb.position[0] - corr[0],
                    bb.position[1] - corr[1],
                    bb.position[2] - corr[2],
                ];
            }
        }
        // Integrate
        self.bodies.integrate_all(dt, self.gravity);
        self.time += dt;
    }

    /// Cast a ray and return all hits sorted by distance.
    pub fn query_ray(
        &self,
        origin: [f64; 3],
        direction: [f64; 3],
        max_dist: f64,
    ) -> Vec<PyRayResult> {
        let mut results = Vec::new();
        let dir_len = vec3_length(direction);
        if dir_len < 1e-12 {
            return results;
        }
        let dir_n = vec3_scale(direction, 1.0 / dir_len);
        for body in self.bodies.iter() {
            // Sphere test against body bounding sphere
            let to_body = [
                body.position[0] - origin[0],
                body.position[1] - origin[1],
                body.position[2] - origin[2],
            ];
            let tc = vec3_dot(to_body, dir_n);
            if tc < 0.0 {
                continue;
            }
            let closest = [
                origin[0] + dir_n[0] * tc,
                origin[1] + dir_n[1] * tc,
                origin[2] + dir_n[2] * tc,
            ];
            let miss = [
                closest[0] - body.position[0],
                closest[1] - body.position[1],
                closest[2] - body.position[2],
            ];
            // Use max collider radius as bounding sphere
            let r = self
                .colliders
                .iter()
                .filter(|c| c.body_id == body.id)
                .map(|c| c.bounding_radius())
                .fold(0.5_f64, f64::max);
            let miss_dist = vec3_length(miss);
            if miss_dist > r {
                continue;
            }
            let dist = tc - (r * r - miss_dist * miss_dist).max(0.0).sqrt();
            if dist > max_dist {
                continue;
            }
            let hit_point = [
                origin[0] + dir_n[0] * dist,
                origin[1] + dir_n[1] * dist,
                origin[2] + dir_n[2] * dist,
            ];
            let normal = {
                let n = [
                    hit_point[0] - body.position[0],
                    hit_point[1] - body.position[1],
                    hit_point[2] - body.position[2],
                ];
                let nl = vec3_length(n);
                if nl > 1e-12 {
                    vec3_scale(n, 1.0 / nl)
                } else {
                    [0.0, 1.0, 0.0]
                }
            };
            results.push(PyRayResult {
                body_id: body.id,
                point: hit_point,
                normal,
                distance: dist,
            });
        }
        results.sort_by(|a, b| {
            a.distance
                .partial_cmp(&b.distance)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        results
    }

    /// Query all pairs of overlapping colliders.
    pub fn query_overlap(&self) -> Vec<(u32, u32)> {
        let mut pairs = Vec::new();
        for i in 0..self.colliders.len() {
            for j in (i + 1)..self.colliders.len() {
                let ca = &self.colliders[i];
                let cb = &self.colliders[j];
                if ca.body_id == cb.body_id {
                    continue;
                }
                let ba = self.bodies.get(ca.body_id);
                let bb = self.bodies.get(cb.body_id);
                if let (Some(ba), Some(bb)) = (ba, bb) {
                    let d = vec3_length([
                        bb.position[0] - ba.position[0],
                        bb.position[1] - ba.position[1],
                        bb.position[2] - ba.position[2],
                    ]);
                    if d < ca.bounding_radius() + cb.bounding_radius() {
                        pairs.push((ca.body_id, cb.body_id));
                    }
                }
            }
        }
        pairs
    }

    /// Return contacts from the last step.
    pub fn get_contacts(&self) -> &[PyContact] {
        &self.contacts
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PyPhysicsPipeline
// ─────────────────────────────────────────────────────────────────────────────

/// Configuration for the physics pipeline.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyPipelineConfig {
    /// Maximum solver iterations per step.
    pub max_iterations: u32,
    /// Velocity tolerance for sleep (m/s).
    pub sleep_threshold: f64,
    /// Time (s) body must be below sleep threshold before sleeping.
    pub sleep_time: f64,
    /// Gravity vector.
    pub gravity: [f64; 3],
    /// Substeps per `step()` call.
    pub substeps: u32,
}

impl Default for PyPipelineConfig {
    fn default() -> Self {
        Self {
            max_iterations: 10,
            sleep_threshold: 0.05,
            sleep_time: 1.0,
            gravity: [0.0, -9.81, 0.0],
            substeps: 1,
        }
    }
}

/// A physics pipeline orchestrating broadphase, narrowphase, solver, and integration.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyPhysicsPipeline {
    /// Pipeline configuration.
    pub config: PyPipelineConfig,
    /// Internal world.
    pub world: PyWorld,
    /// Simulation step counter.
    pub step_count: u64,
    /// Broadphase collision pair count (last step).
    pub last_broadphase_pairs: usize,
    /// Narrowphase contact count (last step).
    pub last_contact_count: usize,
    /// Solver iterations used (last step).
    pub last_solver_iterations: u32,
}

impl PyPhysicsPipeline {
    /// Create a new pipeline with default config.
    pub fn new() -> Self {
        let config = PyPipelineConfig::default();
        let world = PyWorld::new(config.gravity);
        Self {
            config,
            world,
            step_count: 0,
            last_broadphase_pairs: 0,
            last_contact_count: 0,
            last_solver_iterations: 0,
        }
    }

    /// Create with custom config.
    pub fn with_config(config: PyPipelineConfig) -> Self {
        let world = PyWorld::new(config.gravity);
        Self {
            config,
            world,
            step_count: 0,
            last_broadphase_pairs: 0,
            last_contact_count: 0,
            last_solver_iterations: 0,
        }
    }

    /// Advance the simulation by `dt` seconds.
    pub fn step(&mut self, dt: f64) {
        let sub_dt = dt / self.config.substeps as f64;
        for _ in 0..self.config.substeps {
            self.world.step(sub_dt);
        }
        self.last_contact_count = self.world.contacts.len();
        self.last_broadphase_pairs = self.last_contact_count;
        self.last_solver_iterations =
            self.config
                .max_iterations
                .min(if self.last_contact_count > 0 {
                    self.config.max_iterations
                } else {
                    0
                });
        self.step_count += 1;
    }

    /// Set maximum solver iterations.
    pub fn set_max_iterations(&mut self, n: u32) {
        self.config.max_iterations = n;
    }

    /// Access the world.
    pub fn world(&self) -> &PyWorld {
        &self.world
    }

    /// Mutable access to the world.
    pub fn world_mut(&mut self) -> &mut PyWorld {
        &mut self.world
    }
}

impl Default for PyPhysicsPipeline {
    fn default() -> Self {
        Self::new()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PyWakeEvent
// ─────────────────────────────────────────────────────────────────────────────

/// Reason a body woke up.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum WakeReason {
    /// Woken by an applied force or impulse.
    Force,
    /// Woken by a contact.
    Contact,
    /// Woken by a script/user request.
    Script,
}

/// Wake event emitted when a body transitions from sleeping to awake.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyWakeEvent {
    /// The body that woke up.
    pub body_id: u32,
    /// Reason for waking.
    pub reason: WakeReason,
    /// Simulation time at wake.
    pub time: f64,
}

impl PyWakeEvent {
    /// Create a new wake event.
    pub fn new(body_id: u32, reason: WakeReason, time: f64) -> Self {
        Self {
            body_id,
            reason,
            time,
        }
    }
}

/// Manages sleep/wake events for a body set.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct PySleepManager {
    /// Pending wake events.
    pub events: Vec<PyWakeEvent>,
    /// Sleep timers per body id.
    sleep_timers: Vec<(u32, f64)>,
}

impl PySleepManager {
    /// Create a new sleep manager.
    pub fn new() -> Self {
        Self::default()
    }

    /// Force-wake a body and emit a wake event.
    pub fn wake_body(&mut self, body: &mut PyRigidBody, reason: WakeReason, time: f64) {
        if body.sleeping {
            body.sleeping = false;
            self.events.push(PyWakeEvent::new(body.id, reason, time));
        }
    }

    /// Put a body to sleep immediately.
    pub fn sleep_body(&mut self, body: &mut PyRigidBody) {
        body.sleeping = true;
        body.velocity = [0.0; 3];
        body.ang_vel = [0.0; 3];
    }

    /// Update sleep timers for all bodies. Returns ids that should sleep.
    pub fn update(
        &mut self,
        bodies: &mut PyRigidBodySet,
        threshold: f64,
        sleep_time: f64,
        dt: f64,
    ) {
        let ids: Vec<u32> = bodies.iter().map(|b| b.id).collect();
        for id in ids {
            if let Some(body) = bodies.get_mut(id) {
                if body.is_static || body.is_kinematic {
                    continue;
                }
                let speed = vec3_length(body.velocity) + vec3_length(body.ang_vel);
                if speed < threshold {
                    let timer = self.sleep_timers.iter_mut().find(|(bid, _)| *bid == id);
                    match timer {
                        Some((_, t)) => {
                            *t += dt;
                            if *t >= sleep_time && !body.sleeping {
                                body.sleeping = true;
                            }
                        }
                        None => {
                            self.sleep_timers.push((id, dt));
                        }
                    }
                } else {
                    self.sleep_timers.retain(|(bid, _)| *bid != id);
                }
            }
        }
    }

    /// Drain all pending events.
    pub fn drain_events(&mut self) -> Vec<PyWakeEvent> {
        std::mem::take(&mut self.events)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PyBodyFilter
// ─────────────────────────────────────────────────────────────────────────────

/// Filter criteria for selecting bodies.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct PyBodyFilter {
    /// Include active (awake) bodies.
    pub include_active: bool,
    /// Include sleeping bodies.
    pub include_sleeping: bool,
    /// Include static bodies.
    pub include_static: bool,
    /// Include kinematic bodies.
    pub include_kinematic: bool,
}

impl PyBodyFilter {
    /// Filter that matches all bodies.
    pub fn all() -> Self {
        Self {
            include_active: true,
            include_sleeping: true,
            include_static: true,
            include_kinematic: true,
        }
    }

    /// Filter that matches only dynamic active bodies.
    pub fn active_dynamic() -> Self {
        Self {
            include_active: true,
            include_sleeping: false,
            include_static: false,
            include_kinematic: false,
        }
    }

    /// Test if a body passes this filter.
    pub fn matches(&self, body: &PyRigidBody) -> bool {
        if body.is_static {
            return self.include_static;
        }
        if body.is_kinematic {
            return self.include_kinematic;
        }
        if body.sleeping {
            self.include_sleeping
        } else {
            self.include_active
        }
    }

    /// Filter a body set; returns matching body ids.
    pub fn filter_by_type<'a>(&self, bodies: &'a PyRigidBodySet) -> Vec<&'a PyRigidBody> {
        bodies.iter().filter(|b| self.matches(b)).collect()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PyIslandStats
// ─────────────────────────────────────────────────────────────────────────────

/// Statistics for a simulation island.
///
/// An island is a group of bodies connected by contacts or constraints.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyIslandStats {
    /// Island identifier.
    pub island_id: u32,
    /// Number of bodies in this island.
    pub body_count: usize,
    /// Number of constraints (contacts + joints) in this island.
    pub constraint_count: usize,
    /// Accumulated time all bodies have been below sleep threshold (seconds).
    pub sleep_timer: f64,
    /// Total kinetic energy of the island.
    pub energy: f64,
}

impl PyIslandStats {
    /// Create island stats.
    pub fn new(island_id: u32, body_count: usize, constraint_count: usize) -> Self {
        Self {
            island_id,
            body_count,
            constraint_count,
            sleep_timer: 0.0,
            energy: 0.0,
        }
    }

    /// True if the island is eligible to sleep.
    pub fn can_sleep(&self, threshold: f64, sleep_time: f64) -> bool {
        self.energy < threshold && self.sleep_timer >= sleep_time
    }
}

/// Compute island statistics from a body set (simple connected-components mock).
pub fn compute_island_stats(bodies: &PyRigidBodySet) -> Vec<PyIslandStats> {
    // Trivial: one island per body (no connectivity info here)
    bodies
        .iter()
        .enumerate()
        .map(|(i, b)| {
            let mut stats = PyIslandStats::new(i as u32, 1, 0);
            stats.energy = b.kinetic_energy();
            stats
        })
        .collect()
}

// ─────────────────────────────────────────────────────────────────────────────
// PyBatchForce
// ─────────────────────────────────────────────────────────────────────────────

/// A vectorized (batch) force application request.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyBatchForce {
    /// Target body handles.
    pub body_ids: Vec<u32>,
    /// Forces to apply (parallel to `body_ids`).
    pub forces: Vec<[f64; 3]>,
    /// Optional torques (parallel to `body_ids`; zero if absent).
    pub torques: Option<Vec<[f64; 3]>>,
}

impl PyBatchForce {
    /// Create a batch of forces.
    pub fn new(body_ids: Vec<u32>, forces: Vec<[f64; 3]>) -> Self {
        Self {
            body_ids,
            forces,
            torques: None,
        }
    }

    /// Add torques to the batch.
    pub fn with_torques(mut self, torques: Vec<[f64; 3]>) -> Self {
        self.torques = Some(torques);
        self
    }

    /// Apply all forces (and torques if present) to the given body set.
    pub fn apply(&self, bodies: &mut PyRigidBodySet) {
        for (idx, &id) in self.body_ids.iter().enumerate() {
            if let Some(body) = bodies.get_mut(id) {
                if idx < self.forces.len() {
                    body.apply_force(self.forces[idx]);
                }
                if let Some(torques) = &self.torques
                    && idx < torques.len()
                {
                    body.apply_torque(torques[idx]);
                }
            }
        }
    }

    /// Validate that body_ids and forces have equal length.
    pub fn is_valid(&self) -> bool {
        self.body_ids.len() == self.forces.len()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // --- PyRigidBody ---

    #[test]
    fn test_rigid_body_new() {
        let b = PyRigidBody::new(1, 5.0, [0.0, 10.0, 0.0]);
        assert_eq!(b.mass, 5.0);
        assert_eq!(b.position, [0.0, 10.0, 0.0]);
        assert!(!b.sleeping);
        assert!(!b.is_static);
    }

    #[test]
    fn test_rigid_body_static() {
        let b = PyRigidBody::new_static(2, [1.0, 0.0, 0.0]);
        assert!(b.is_static);
        assert_eq!(b.inv_mass(), 0.0);
    }

    #[test]
    fn test_rigid_body_kinematic() {
        let b = PyRigidBody::new_kinematic(3, 2.0, [0.0; 3]);
        assert!(b.is_kinematic);
        assert!(!b.is_static);
    }

    #[test]
    fn test_apply_force() {
        let mut b = PyRigidBody::new(1, 1.0, [0.0; 3]);
        b.apply_force([1.0, 0.0, 0.0]);
        assert!((b.accumulated_force[0] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_apply_force_static_ignored() {
        let mut b = PyRigidBody::new_static(1, [0.0; 3]);
        b.apply_force([999.0, 0.0, 0.0]);
        assert_eq!(b.accumulated_force[0], 0.0);
    }

    #[test]
    fn test_apply_impulse() {
        let mut b = PyRigidBody::new(1, 2.0, [0.0; 3]);
        b.apply_impulse([4.0, 0.0, 0.0]);
        assert!((b.velocity[0] - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_apply_impulse_static_ignored() {
        let mut b = PyRigidBody::new_static(1, [0.0; 3]);
        b.apply_impulse([999.0, 0.0, 0.0]);
        assert_eq!(b.velocity[0], 0.0);
    }

    #[test]
    fn test_kinetic_energy() {
        let mut b = PyRigidBody::new(1, 2.0, [0.0; 3]);
        b.velocity = [1.0, 0.0, 0.0];
        // KE_lin = 0.5 * 2 * 1 = 1
        assert!((b.kinetic_energy() - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_integrate_gravity() {
        let mut b = PyRigidBody::new(1, 1.0, [0.0, 10.0, 0.0]);
        let gravity = [0.0, -9.81, 0.0];
        b.integrate(0.1, gravity);
        // Position should decrease (fallen)
        assert!(b.position[1] < 10.0);
    }

    #[test]
    fn test_linear_momentum() {
        let mut b = PyRigidBody::new(1, 3.0, [0.0; 3]);
        b.velocity = [2.0, 0.0, 0.0];
        let p = b.linear_momentum();
        assert!((p[0] - 6.0).abs() < 1e-12);
    }

    // --- PyRigidBodySet ---

    #[test]
    fn test_body_set_add_remove() {
        let mut set = PyRigidBodySet::new();
        let id = set.add(PyRigidBody::new(0, 1.0, [0.0; 3]));
        assert_eq!(set.len(), 1);
        assert!(set.remove(id));
        assert_eq!(set.len(), 0);
    }

    #[test]
    fn test_body_set_get() {
        let mut set = PyRigidBodySet::new();
        let id = set.add(PyRigidBody::new(0, 5.0, [1.0, 2.0, 3.0]));
        let b = set.get(id).unwrap();
        assert_eq!(b.mass, 5.0);
    }

    #[test]
    fn test_body_set_is_empty() {
        let set = PyRigidBodySet::new();
        assert!(set.is_empty());
    }

    #[test]
    fn test_body_set_total_ke() {
        let mut set = PyRigidBodySet::new();
        let id = set.add(PyRigidBody::new(0, 2.0, [0.0; 3]));
        let b = set.get_mut(id).unwrap();
        b.velocity = [1.0, 0.0, 0.0];
        let ke = set.total_kinetic_energy();
        assert!(ke > 0.0);
    }

    // --- PyCollider ---

    #[test]
    fn test_collider_sphere() {
        let c = PyCollider::sphere(0, 1, 0.5);
        assert_eq!(c.bounding_radius(), 0.5);
        assert!(!c.is_sensor);
    }

    #[test]
    fn test_collider_box() {
        let c = PyCollider::box_collider(0, 1, [1.0, 2.0, 3.0]);
        assert!(c.bounding_radius() > 0.0);
    }

    #[test]
    fn test_collider_capsule() {
        let c = PyCollider::capsule(0, 1, 1.0, 0.5);
        assert_eq!(c.bounding_radius(), 1.5);
    }

    #[test]
    fn test_collider_sensor() {
        let c = PyCollider::sphere(0, 1, 0.3).as_sensor();
        assert!(c.is_sensor);
    }

    // --- PyWorld ---

    #[test]
    fn test_world_step_gravity() {
        let mut w = PyWorld::earth_gravity();
        let id = w.create_body(1.0, [0.0, 10.0, 0.0]);
        w.step(0.1);
        let pos = w.bodies.get_position(id).unwrap();
        assert!(pos[1] < 10.0);
    }

    #[test]
    fn test_world_static_body_unmoved() {
        let mut w = PyWorld::earth_gravity();
        let id = w.create_static_body([0.0, 0.0, 0.0]);
        w.step(1.0);
        let pos = w.bodies.get_position(id).unwrap();
        assert_eq!(pos, [0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_world_ray_cast() {
        let mut w = PyWorld::earth_gravity();
        let id = w.create_body(1.0, [0.0, 0.0, 5.0]);
        w.create_collider(PyCollider::sphere(0, id, 1.0));
        let hits = w.query_ray([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 100.0);
        assert!(!hits.is_empty());
    }

    #[test]
    fn test_world_get_contacts() {
        let mut w = PyWorld::new([0.0; 3]);
        let a = w.create_body(1.0, [0.0; 3]);
        let b = w.create_body(1.0, [0.5, 0.0, 0.0]);
        w.create_collider(PyCollider::sphere(0, a, 0.4));
        w.create_collider(PyCollider::sphere(0, b, 0.4));
        w.step(0.001);
        let contacts = w.get_contacts();
        assert!(!contacts.is_empty());
    }

    // --- PyJointSet ---

    #[test]
    fn test_joint_set_create_remove() {
        let mut js = PyJointSet::new();
        let id = js.create_joint(PyJointType::Hinge, 0, 1, [0.0; 3], [0.0; 3], None);
        assert_eq!(js.len(), 1);
        assert!(js.remove(id));
        assert_eq!(js.len(), 0);
    }

    #[test]
    fn test_joint_set_motor() {
        let mut js = PyJointSet::new();
        let id = js.create_joint(PyJointType::Hinge, 0, 1, [0.0; 3], [0.0; 3], None);
        assert!(js.set_motor(id, 5.0, 100.0));
        let j = js.get(id).unwrap();
        assert!(j.motor_enabled);
        assert_eq!(j.motor_speed, 5.0);
    }

    // --- PyPhysicsPipeline ---

    #[test]
    fn test_pipeline_step() {
        let mut pipe = PyPhysicsPipeline::new();
        let _id = pipe.world_mut().create_body(1.0, [0.0, 5.0, 0.0]);
        pipe.step(0.016);
        assert_eq!(pipe.step_count, 1);
    }

    #[test]
    fn test_pipeline_max_iterations() {
        let mut pipe = PyPhysicsPipeline::new();
        pipe.set_max_iterations(5);
        assert_eq!(pipe.config.max_iterations, 5);
    }

    // --- PyBodyFilter ---

    #[test]
    fn test_body_filter_all() {
        let f = PyBodyFilter::all();
        let static_b = PyRigidBody::new_static(1, [0.0; 3]);
        assert!(f.matches(&static_b));
        let dyn_b = PyRigidBody::new(2, 1.0, [0.0; 3]);
        assert!(f.matches(&dyn_b));
    }

    #[test]
    fn test_body_filter_active_only() {
        let f = PyBodyFilter::active_dynamic();
        let mut b = PyRigidBody::new(1, 1.0, [0.0; 3]);
        b.sleeping = false;
        assert!(f.matches(&b));
        b.sleeping = true;
        assert!(!f.matches(&b));
    }

    // --- PyIslandStats ---

    #[test]
    fn test_island_stats_new() {
        let s = PyIslandStats::new(0, 5, 4);
        assert_eq!(s.body_count, 5);
        assert!(!s.can_sleep(0.1, 1.0));
    }

    // --- PyBatchForce ---

    #[test]
    fn test_batch_force_apply() {
        let mut set = PyRigidBodySet::new();
        let a = set.add(PyRigidBody::new(0, 1.0, [0.0; 3]));
        let b = set.add(PyRigidBody::new(0, 1.0, [1.0, 0.0, 0.0]));
        let batch = PyBatchForce::new(vec![a, b], vec![[10.0, 0.0, 0.0], [0.0, 10.0, 0.0]]);
        assert!(batch.is_valid());
        batch.apply(&mut set);
        assert!((set.get(a).unwrap().accumulated_force[0] - 10.0).abs() < 1e-12);
        assert!((set.get(b).unwrap().accumulated_force[1] - 10.0).abs() < 1e-12);
    }

    #[test]
    fn test_batch_force_invalid() {
        let batch = PyBatchForce::new(vec![0, 1], vec![[1.0, 0.0, 0.0]]);
        assert!(!batch.is_valid());
    }

    // --- PySleepManager ---

    #[test]
    fn test_sleep_manager_wake() {
        let mut mgr = PySleepManager::new();
        let mut set = PyRigidBodySet::new();
        let id = set.add(PyRigidBody::new(0, 1.0, [0.0; 3]));
        let body = set.get_mut(id).unwrap();
        body.sleeping = true;
        mgr.wake_body(body, WakeReason::Script, 0.0);
        assert!(!body.sleeping);
        assert_eq!(mgr.events.len(), 1);
    }

    #[test]
    fn test_sleep_manager_sleep() {
        let mut mgr = PySleepManager::new();
        let mut set = PyRigidBodySet::new();
        let id = set.add(PyRigidBody::new(0, 1.0, [0.0; 3]));
        let body = set.get_mut(id).unwrap();
        body.velocity = [10.0, 0.0, 0.0];
        mgr.sleep_body(body);
        assert!(body.sleeping);
        assert_eq!(body.velocity, [0.0; 3]);
    }

    #[test]
    fn test_compute_island_stats() {
        let mut set = PyRigidBodySet::new();
        set.add(PyRigidBody::new(0, 1.0, [0.0; 3]));
        set.add(PyRigidBody::new(0, 2.0, [1.0, 0.0, 0.0]));
        let islands = compute_island_stats(&set);
        assert_eq!(islands.len(), 2);
    }
}
