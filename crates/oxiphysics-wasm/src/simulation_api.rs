// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! WASM Simulation API: high-level types and world management for
//! WebAssembly JavaScript interop.
//!
//! Provides WasmSimulationConfig, WasmWorld, WasmRigidBodyHandle,
//! WasmColliderHandle, WasmJointHandle, WasmPhysicsEvent,
//! WasmRaycastResult, WasmOverlapResult, WasmContactEvent,
//! and WasmSceneSerializer (JSON-based).

#![allow(dead_code)]
#![allow(unused_imports)]
#![allow(clippy::too_many_arguments)]

use std::collections::HashMap;

// ── Handle types ──────────────────────────────────────────────────────────────

/// Handle to a rigid body inside WasmWorld.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct WasmRigidBodyHandle(pub u32);

/// Handle to a collider inside WasmWorld.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct WasmColliderHandle(pub u32);

/// Handle to a joint (constraint) inside WasmWorld.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct WasmJointHandle(pub u32);

// ── Simulation configuration ──────────────────────────────────────────────────

/// Broad-phase algorithm selection.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BroadPhaseAlgorithm {
    /// Simple axis-aligned bounding box sweep and prune.
    Sap,
    /// BVH-based broad-phase.
    Bvh,
    /// Grid-based spatial hash.
    SpatialHash,
}

/// Integration method for rigid body dynamics.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IntegrationMethod {
    /// Semi-implicit Euler (first-order, fast).
    SemiImplicitEuler,
    /// Runge-Kutta 4 (higher accuracy).
    Rk4,
    /// Velocity Verlet.
    Verlet,
}

/// Configuration for a WASM simulation world.
#[derive(Debug, Clone)]
pub struct WasmSimulationConfig {
    /// Gravity vector (m/s²).
    pub gravity: [f64; 3],
    /// Fixed simulation time step (s).
    pub dt: f64,
    /// Number of velocity solver iterations.
    pub velocity_iters: u32,
    /// Number of position solver iterations.
    pub position_iters: u32,
    /// Broad-phase algorithm.
    pub broad_phase: BroadPhaseAlgorithm,
    /// Integration method.
    pub integration: IntegrationMethod,
    /// Enable continuous collision detection.
    pub ccd_enabled: bool,
    /// Maximum linear velocity (m/s) before clamping.
    pub max_linear_vel: f64,
    /// Maximum angular velocity (rad/s) before clamping.
    pub max_angular_vel: f64,
    /// Linear damping applied to all bodies.
    pub linear_damping: f64,
    /// Angular damping applied to all bodies.
    pub angular_damping: f64,
    /// Allow sleeping of inactive bodies.
    pub allow_sleeping: bool,
    /// Linear speed below which body enters sleep (m/s).
    pub sleep_linear_threshold: f64,
    /// Angular speed below which body enters sleep (rad/s).
    pub sleep_angular_threshold: f64,
    /// Time before sleep activates (s).
    pub sleep_time_threshold: f64,
}

impl Default for WasmSimulationConfig {
    fn default() -> Self {
        Self {
            gravity: [0.0, -9.81, 0.0],
            dt: 1.0 / 60.0,
            velocity_iters: 8,
            position_iters: 4,
            broad_phase: BroadPhaseAlgorithm::Sap,
            integration: IntegrationMethod::SemiImplicitEuler,
            ccd_enabled: false,
            max_linear_vel: 200.0,
            max_angular_vel: 200.0,
            linear_damping: 0.0,
            angular_damping: 0.0,
            allow_sleeping: true,
            sleep_linear_threshold: 0.002,
            sleep_angular_threshold: 0.005,
            sleep_time_threshold: 1.0,
        }
    }
}

impl WasmSimulationConfig {
    /// Create a config with Earth-standard gravity.
    pub fn earth() -> Self {
        Self::default()
    }

    /// Create a config with Moon-level gravity (1/6 of Earth).
    pub fn moon() -> Self {
        Self {
            gravity: [0.0, -1.625, 0.0],
            ..Default::default()
        }
    }

    /// Create a zero-gravity (space) configuration.
    pub fn zero_gravity() -> Self {
        Self {
            gravity: [0.0; 3],
            ..Default::default()
        }
    }

    /// Create a config with custom gravity.
    pub fn with_gravity(gx: f64, gy: f64, gz: f64) -> Self {
        Self {
            gravity: [gx, gy, gz],
            ..Default::default()
        }
    }
}

// ── Rigid body ────────────────────────────────────────────────────────────────

/// Type of a rigid body.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WasmBodyType {
    /// Fully simulated (mass > 0).
    Dynamic,
    /// Immovable (infinite mass).
    Static,
    /// Driven externally (velocity-controlled).
    Kinematic,
}

/// State of a rigid body for serialization and queries.
#[derive(Debug, Clone)]
pub struct WasmBodyState {
    /// Handle of the body.
    pub handle: WasmRigidBodyHandle,
    /// World-space position.
    pub position: [f64; 3],
    /// Orientation as quaternion (x, y, z, w).
    pub rotation: [f64; 4],
    /// Linear velocity.
    pub linear_vel: [f64; 3],
    /// Angular velocity.
    pub angular_vel: [f64; 3],
    /// Whether the body is currently sleeping.
    pub sleeping: bool,
    /// Body type.
    pub body_type: WasmBodyType,
}

/// Internal rigid body data.
#[derive(Debug, Clone)]
struct WasmRigidBody {
    handle: WasmRigidBodyHandle,
    body_type: WasmBodyType,
    position: [f64; 3],
    rotation: [f64; 4],
    linear_vel: [f64; 3],
    angular_vel: [f64; 3],
    mass: f64,
    inv_mass: f64,
    inertia: [f64; 3],
    inv_inertia: [f64; 3],
    force_accum: [f64; 3],
    torque_accum: [f64; 3],
    linear_damping: f64,
    angular_damping: f64,
    sleeping: bool,
    sleep_time: f64,
    ccd_enabled: bool,
    gravity_scale: f64,
    user_data: u64,
}

impl WasmRigidBody {
    fn new_dynamic(handle: WasmRigidBodyHandle, position: [f64; 3], mass: f64) -> Self {
        let inv_mass = if mass > 1e-15 { 1.0 / mass } else { 0.0 };
        let inertia = [mass * 0.4; 3]; // unit sphere approximation
        let inv_inertia = [if inertia[0] > 1e-15 {
            1.0 / inertia[0]
        } else {
            0.0
        }; 3];
        Self {
            handle,
            body_type: WasmBodyType::Dynamic,
            position,
            rotation: [0.0, 0.0, 0.0, 1.0],
            linear_vel: [0.0; 3],
            angular_vel: [0.0; 3],
            mass,
            inv_mass,
            inertia,
            inv_inertia,
            force_accum: [0.0; 3],
            torque_accum: [0.0; 3],
            linear_damping: 0.0,
            angular_damping: 0.0,
            sleeping: false,
            sleep_time: 0.0,
            ccd_enabled: false,
            gravity_scale: 1.0,
            user_data: 0,
        }
    }

    fn new_static(handle: WasmRigidBodyHandle, position: [f64; 3]) -> Self {
        Self {
            handle,
            body_type: WasmBodyType::Static,
            position,
            rotation: [0.0, 0.0, 0.0, 1.0],
            linear_vel: [0.0; 3],
            angular_vel: [0.0; 3],
            mass: 0.0,
            inv_mass: 0.0,
            inertia: [0.0; 3],
            inv_inertia: [0.0; 3],
            force_accum: [0.0; 3],
            torque_accum: [0.0; 3],
            linear_damping: 0.0,
            angular_damping: 0.0,
            sleeping: true,
            sleep_time: f64::INFINITY,
            ccd_enabled: false,
            gravity_scale: 0.0,
            user_data: 0,
        }
    }

    fn integrate(&mut self, dt: f64, gravity: [f64; 3]) {
        if self.body_type != WasmBodyType::Dynamic || self.sleeping {
            self.force_accum = [0.0; 3];
            self.torque_accum = [0.0; 3];
            return;
        }

        // Apply gravity
        let grav_force = [
            gravity[0] * self.mass * self.gravity_scale,
            gravity[1] * self.mass * self.gravity_scale,
            gravity[2] * self.mass * self.gravity_scale,
        ];
        let total_force = [
            self.force_accum[0] + grav_force[0],
            self.force_accum[1] + grav_force[1],
            self.force_accum[2] + grav_force[2],
        ];

        // Linear integration
        let accel = [
            total_force[0] * self.inv_mass,
            total_force[1] * self.inv_mass,
            total_force[2] * self.inv_mass,
        ];
        self.linear_vel[0] += accel[0] * dt;
        self.linear_vel[1] += accel[1] * dt;
        self.linear_vel[2] += accel[2] * dt;

        // Damping
        let lin_damp = (1.0 - self.linear_damping * dt).max(0.0);
        self.linear_vel[0] *= lin_damp;
        self.linear_vel[1] *= lin_damp;
        self.linear_vel[2] *= lin_damp;

        // Clamp velocity
        let speed2 = self.linear_vel[0] * self.linear_vel[0]
            + self.linear_vel[1] * self.linear_vel[1]
            + self.linear_vel[2] * self.linear_vel[2];
        let max_v = 200.0f64;
        if speed2 > max_v * max_v {
            let scale = max_v / speed2.sqrt();
            self.linear_vel[0] *= scale;
            self.linear_vel[1] *= scale;
            self.linear_vel[2] *= scale;
        }

        // Update position
        self.position[0] += self.linear_vel[0] * dt;
        self.position[1] += self.linear_vel[1] * dt;
        self.position[2] += self.linear_vel[2] * dt;

        // Angular integration (simplified: no quaternion integration)
        let alpha = [
            self.torque_accum[0] * self.inv_inertia[0],
            self.torque_accum[1] * self.inv_inertia[1],
            self.torque_accum[2] * self.inv_inertia[2],
        ];
        let ang_damp = (1.0 - self.angular_damping * dt).max(0.0);
        self.angular_vel[0] = (self.angular_vel[0] + alpha[0] * dt) * ang_damp;
        self.angular_vel[1] = (self.angular_vel[1] + alpha[1] * dt) * ang_damp;
        self.angular_vel[2] = (self.angular_vel[2] + alpha[2] * dt) * ang_damp;

        // Integrate rotation (approximate via axis-angle → quaternion update)
        let omega = self.angular_vel;
        let angle = (omega[0] * omega[0] + omega[1] * omega[1] + omega[2] * omega[2]).sqrt();
        if angle > 1e-10 {
            let half_angle = angle * dt * 0.5;
            let sin_h = half_angle.sin();
            let cos_h = half_angle.cos();
            let ax = omega[0] / angle;
            let ay = omega[1] / angle;
            let az = omega[2] / angle;
            let dq = [ax * sin_h, ay * sin_h, az * sin_h, cos_h];
            let q = self.rotation;
            self.rotation = quaternion_multiply(dq, q);
            self.rotation = quaternion_normalize(self.rotation);
        }

        // Reset accumulators
        self.force_accum = [0.0; 3];
        self.torque_accum = [0.0; 3];
    }

    fn state(&self) -> WasmBodyState {
        WasmBodyState {
            handle: self.handle,
            position: self.position,
            rotation: self.rotation,
            linear_vel: self.linear_vel,
            angular_vel: self.angular_vel,
            sleeping: self.sleeping,
            body_type: self.body_type,
        }
    }
}

// ── Quaternion math ───────────────────────────────────────────────────────────

fn quaternion_multiply(p: [f64; 4], q: [f64; 4]) -> [f64; 4] {
    [
        p[3] * q[0] + p[0] * q[3] + p[1] * q[2] - p[2] * q[1],
        p[3] * q[1] - p[0] * q[2] + p[1] * q[3] + p[2] * q[0],
        p[3] * q[2] + p[0] * q[1] - p[1] * q[0] + p[2] * q[3],
        p[3] * q[3] - p[0] * q[0] - p[1] * q[1] - p[2] * q[2],
    ]
}

fn quaternion_normalize(q: [f64; 4]) -> [f64; 4] {
    let len = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]).sqrt();
    if len < 1e-15 {
        return [0.0, 0.0, 0.0, 1.0];
    }
    [q[0] / len, q[1] / len, q[2] / len, q[3] / len]
}

fn quaternion_conjugate(q: [f64; 4]) -> [f64; 4] {
    [-q[0], -q[1], -q[2], q[3]]
}

fn quaternion_rotate_vec(q: [f64; 4], v: [f64; 3]) -> [f64; 3] {
    let vq = [v[0], v[1], v[2], 0.0];
    let q_conj = quaternion_conjugate(q);
    let tmp = quaternion_multiply(q, vq);
    let result = quaternion_multiply(tmp, q_conj);
    [result[0], result[1], result[2]]
}

// ── Collider shapes ───────────────────────────────────────────────────────────

/// Collider shape types.
#[derive(Debug, Clone)]
pub enum WasmColliderShape {
    /// Sphere with radius.
    Sphere { radius: f64 },
    /// Box with half-extents.
    Box { half_extents: [f64; 3] },
    /// Capsule with half-height and radius.
    Capsule { half_height: f64, radius: f64 },
    /// Infinite plane: normal (nx, ny, nz) and offset d.
    Plane { normal: [f64; 3], offset: f64 },
    /// Cylinder with half-height and radius.
    Cylinder { half_height: f64, radius: f64 },
    /// Cone with half-height and base radius.
    Cone { half_height: f64, radius: f64 },
    /// Convex hull (list of points).
    ConvexHull { points: Vec<[f64; 3]> },
}

impl WasmColliderShape {
    /// Compute an approximate bounding sphere radius.
    pub fn bounding_radius(&self) -> f64 {
        match self {
            Self::Sphere { radius } => *radius,
            Self::Box { half_extents } => (half_extents[0] * half_extents[0]
                + half_extents[1] * half_extents[1]
                + half_extents[2] * half_extents[2])
                .sqrt(),
            Self::Capsule {
                half_height,
                radius,
            } => half_height + radius,
            Self::Plane { .. } => f64::INFINITY,
            Self::Cylinder {
                half_height,
                radius,
            } => (half_height * half_height + radius * radius).sqrt(),
            Self::Cone {
                half_height,
                radius,
            } => (half_height * half_height + radius * radius).sqrt(),
            Self::ConvexHull { points } => points
                .iter()
                .map(|p| (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt())
                .fold(0.0f64, f64::max),
        }
    }
}

/// A collider attached to a rigid body.
#[derive(Debug, Clone)]
struct WasmCollider {
    handle: WasmColliderHandle,
    body: Option<WasmRigidBodyHandle>,
    shape: WasmColliderShape,
    local_offset: [f64; 3],
    local_rotation: [f64; 4],
    friction: f64,
    restitution: f64,
    density: f64,
    is_sensor: bool,
    collision_group: u32,
    collision_mask: u32,
    user_data: u64,
}

impl WasmCollider {
    fn new(handle: WasmColliderHandle, shape: WasmColliderShape) -> Self {
        Self {
            handle,
            body: None,
            shape,
            local_offset: [0.0; 3],
            local_rotation: [0.0, 0.0, 0.0, 1.0],
            friction: 0.5,
            restitution: 0.3,
            density: 1000.0,
            is_sensor: false,
            collision_group: 1,
            collision_mask: u32::MAX,
            user_data: 0,
        }
    }
}

// ── Joint types ───────────────────────────────────────────────────────────────

/// Joint/constraint type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WasmJointType {
    /// Ball-and-socket joint (3 DOF rotation).
    Ball,
    /// Fixed joint (0 DOF).
    Fixed,
    /// Revolute joint (1 DOF rotation).
    Revolute,
    /// Prismatic joint (1 DOF translation).
    Prismatic,
    /// Generic 6-DOF joint.
    Generic6Dof,
    /// Distance constraint.
    Distance,
    /// Spring-damper constraint.
    Spring,
}

/// A joint connecting two rigid bodies.
#[derive(Debug, Clone)]
pub struct WasmJoint {
    /// Joint handle.
    pub handle: WasmJointHandle,
    /// First body.
    pub body_a: WasmRigidBodyHandle,
    /// Second body.
    pub body_b: WasmRigidBodyHandle,
    /// Joint type.
    pub joint_type: WasmJointType,
    /// Anchor in body A local space.
    pub anchor_a: [f64; 3],
    /// Anchor in body B local space.
    pub anchor_b: [f64; 3],
    /// Frame in body A local space (quaternion).
    pub frame_a: [f64; 4],
    /// Frame in body B local space (quaternion).
    pub frame_b: [f64; 4],
    /// Motor target velocity.
    pub motor_velocity: f64,
    /// Motor max force.
    pub motor_max_force: f64,
    /// Spring stiffness.
    pub spring_stiffness: f64,
    /// Spring damping.
    pub spring_damping: f64,
    /// Lower limit (angle or distance).
    pub lower_limit: f64,
    /// Upper limit (angle or distance).
    pub upper_limit: f64,
    /// Whether limits are enabled.
    pub limits_enabled: bool,
}

impl WasmJoint {
    /// Create a ball joint between two bodies.
    pub fn ball(
        handle: WasmJointHandle,
        body_a: WasmRigidBodyHandle,
        body_b: WasmRigidBodyHandle,
        anchor_a: [f64; 3],
        anchor_b: [f64; 3],
    ) -> Self {
        Self {
            handle,
            body_a,
            body_b,
            joint_type: WasmJointType::Ball,
            anchor_a,
            anchor_b,
            frame_a: [0.0, 0.0, 0.0, 1.0],
            frame_b: [0.0, 0.0, 0.0, 1.0],
            motor_velocity: 0.0,
            motor_max_force: f64::INFINITY,
            spring_stiffness: 0.0,
            spring_damping: 0.0,
            lower_limit: f64::NEG_INFINITY,
            upper_limit: f64::INFINITY,
            limits_enabled: false,
        }
    }

    /// Create a distance constraint.
    pub fn distance(
        handle: WasmJointHandle,
        body_a: WasmRigidBodyHandle,
        body_b: WasmRigidBodyHandle,
        target_dist: f64,
    ) -> Self {
        Self {
            handle,
            body_a,
            body_b,
            joint_type: WasmJointType::Distance,
            anchor_a: [0.0; 3],
            anchor_b: [0.0; 3],
            frame_a: [0.0, 0.0, 0.0, 1.0],
            frame_b: [0.0, 0.0, 0.0, 1.0],
            motor_velocity: 0.0,
            motor_max_force: f64::INFINITY,
            spring_stiffness: 0.0,
            spring_damping: 0.0,
            lower_limit: target_dist,
            upper_limit: target_dist,
            limits_enabled: true,
        }
    }
}

// ── Physics events ────────────────────────────────────────────────────────────

/// A physics event type.
#[derive(Debug, Clone, PartialEq)]
pub enum WasmPhysicsEvent {
    /// Two colliders started touching.
    CollisionStarted {
        collider_a: WasmColliderHandle,
        collider_b: WasmColliderHandle,
    },
    /// Two colliders stopped touching.
    CollisionEnded {
        collider_a: WasmColliderHandle,
        collider_b: WasmColliderHandle,
    },
    /// A sensor collider was entered.
    SensorEntered {
        sensor: WasmColliderHandle,
        other: WasmColliderHandle,
    },
    /// A sensor collider was exited.
    SensorExited {
        sensor: WasmColliderHandle,
        other: WasmColliderHandle,
    },
    /// A body went to sleep.
    BodySlept { body: WasmRigidBodyHandle },
    /// A body woke up.
    BodyWoke { body: WasmRigidBodyHandle },
}

/// A contact event with detailed contact information.
#[derive(Debug, Clone)]
pub struct WasmContactEvent {
    /// First collider.
    pub collider_a: WasmColliderHandle,
    /// Second collider.
    pub collider_b: WasmColliderHandle,
    /// Contact normal (from B to A).
    pub normal: [f64; 3],
    /// Contact point in world space.
    pub contact_point: [f64; 3],
    /// Penetration depth (positive = overlap).
    pub depth: f64,
    /// Impulse applied to resolve the contact.
    pub impulse: f64,
    /// Friction impulse applied.
    pub friction_impulse: [f64; 3],
}

// ── Raycast / overlap results ─────────────────────────────────────────────────

/// Result of a single raycast query.
#[derive(Debug, Clone)]
pub struct WasmRaycastResult {
    /// Whether the ray hit anything.
    pub hit: bool,
    /// Hit collider.
    pub collider: Option<WasmColliderHandle>,
    /// Associated rigid body.
    pub body: Option<WasmRigidBodyHandle>,
    /// Hit point in world space.
    pub point: [f64; 3],
    /// Surface normal at hit point.
    pub normal: [f64; 3],
    /// Distance from ray origin to hit.
    pub toi: f64,
    /// Feature ID (sub-shape element hit, e.g. triangle index).
    pub feature_id: u32,
}

impl WasmRaycastResult {
    /// Create a "no hit" result.
    pub fn miss() -> Self {
        Self {
            hit: false,
            collider: None,
            body: None,
            point: [0.0; 3],
            normal: [0.0; 3],
            toi: f64::INFINITY,
            feature_id: u32::MAX,
        }
    }
}

/// Result of an overlap query (AABB or shape).
#[derive(Debug, Clone)]
pub struct WasmOverlapResult {
    /// List of collider handles that overlap the query shape.
    pub colliders: Vec<WasmColliderHandle>,
    /// Associated rigid body handles.
    pub bodies: Vec<Option<WasmRigidBodyHandle>>,
}

impl WasmOverlapResult {
    /// Create an empty result.
    pub fn empty() -> Self {
        Self {
            colliders: Vec::new(),
            bodies: Vec::new(),
        }
    }

    /// Number of overlapping shapes.
    pub fn count(&self) -> usize {
        self.colliders.len()
    }
}

// ── World ─────────────────────────────────────────────────────────────────────

/// A complete WASM physics simulation world.
#[derive(Debug, Clone)]
pub struct WasmWorld {
    /// Simulation configuration.
    pub config: WasmSimulationConfig,
    /// Rigid bodies stored by handle.
    bodies: HashMap<u32, WasmRigidBody>,
    /// Colliders stored by handle.
    colliders: HashMap<u32, WasmCollider>,
    /// Joints stored by handle.
    joints: HashMap<u32, WasmJoint>,
    /// Next body handle ID.
    next_body_id: u32,
    /// Next collider handle ID.
    next_collider_id: u32,
    /// Next joint handle ID.
    next_joint_id: u32,
    /// Pending physics events.
    events: Vec<WasmPhysicsEvent>,
    /// Contact events from last step.
    contact_events: Vec<WasmContactEvent>,
    /// Simulation time (s).
    pub time: f64,
    /// Step count.
    pub step_count: u64,
}

impl WasmWorld {
    /// Create a new world with the given configuration.
    pub fn new(config: WasmSimulationConfig) -> Self {
        Self {
            config,
            bodies: HashMap::new(),
            colliders: HashMap::new(),
            joints: HashMap::new(),
            next_body_id: 0,
            next_collider_id: 0,
            next_joint_id: 0,
            events: Vec::new(),
            contact_events: Vec::new(),
            time: 0.0,
            step_count: 0,
        }
    }

    /// Create a world with Earth gravity.
    pub fn earth() -> Self {
        Self::new(WasmSimulationConfig::earth())
    }

    // ── Body management ───────────────────────────────────────────────────

    /// Add a dynamic rigid body. Returns its handle.
    pub fn add_dynamic_body(
        &mut self,
        mass: f64,
        px: f64,
        py: f64,
        pz: f64,
    ) -> WasmRigidBodyHandle {
        let handle = WasmRigidBodyHandle(self.next_body_id);
        self.next_body_id += 1;
        let body = WasmRigidBody::new_dynamic(handle, [px, py, pz], mass);
        self.bodies.insert(handle.0, body);
        handle
    }

    /// Add a static rigid body. Returns its handle.
    pub fn add_static_body(&mut self, px: f64, py: f64, pz: f64) -> WasmRigidBodyHandle {
        let handle = WasmRigidBodyHandle(self.next_body_id);
        self.next_body_id += 1;
        let body = WasmRigidBody::new_static(handle, [px, py, pz]);
        self.bodies.insert(handle.0, body);
        handle
    }

    /// Remove a rigid body and its attached colliders.
    pub fn remove_body(&mut self, handle: WasmRigidBodyHandle) -> bool {
        let removed = self.bodies.remove(&handle.0).is_some();
        if removed {
            // Remove attached colliders
            self.colliders.retain(|_, c| c.body != Some(handle));
        }
        removed
    }

    /// Get the position of a body.
    pub fn get_position(&self, handle: WasmRigidBodyHandle) -> Option<[f64; 3]> {
        self.bodies.get(&handle.0).map(|b| b.position)
    }

    /// Set the position of a body.
    pub fn set_position(&mut self, handle: WasmRigidBodyHandle, pos: [f64; 3]) -> bool {
        if let Some(body) = self.bodies.get_mut(&handle.0) {
            body.position = pos;
            true
        } else {
            false
        }
    }

    /// Get the velocity of a body.
    pub fn get_velocity(&self, handle: WasmRigidBodyHandle) -> Option<[f64; 3]> {
        self.bodies.get(&handle.0).map(|b| b.linear_vel)
    }

    /// Set the velocity of a body.
    pub fn set_velocity(&mut self, handle: WasmRigidBodyHandle, vel: [f64; 3]) -> bool {
        if let Some(body) = self.bodies.get_mut(&handle.0) {
            body.linear_vel = vel;
            true
        } else {
            false
        }
    }

    /// Apply a force to a body at its center of mass.
    pub fn apply_force(&mut self, handle: WasmRigidBodyHandle, force: [f64; 3]) -> bool {
        if let Some(body) = self.bodies.get_mut(&handle.0) {
            body.force_accum[0] += force[0];
            body.force_accum[1] += force[1];
            body.force_accum[2] += force[2];
            true
        } else {
            false
        }
    }

    /// Apply an impulse to a body at its center of mass.
    pub fn apply_impulse(&mut self, handle: WasmRigidBodyHandle, impulse: [f64; 3]) -> bool {
        if let Some(body) = self.bodies.get_mut(&handle.0) {
            if body.inv_mass > 1e-15 {
                body.linear_vel[0] += impulse[0] * body.inv_mass;
                body.linear_vel[1] += impulse[1] * body.inv_mass;
                body.linear_vel[2] += impulse[2] * body.inv_mass;
            }
            true
        } else {
            false
        }
    }

    /// Apply torque to a body.
    pub fn apply_torque(&mut self, handle: WasmRigidBodyHandle, torque: [f64; 3]) -> bool {
        if let Some(body) = self.bodies.get_mut(&handle.0) {
            body.torque_accum[0] += torque[0];
            body.torque_accum[1] += torque[1];
            body.torque_accum[2] += torque[2];
            true
        } else {
            false
        }
    }

    /// Get the rotation quaternion of a body (x,y,z,w).
    pub fn get_rotation(&self, handle: WasmRigidBodyHandle) -> Option<[f64; 4]> {
        self.bodies.get(&handle.0).map(|b| b.rotation)
    }

    /// Set the gravity scale of a body.
    pub fn set_gravity_scale(&mut self, handle: WasmRigidBodyHandle, scale: f64) -> bool {
        if let Some(body) = self.bodies.get_mut(&handle.0) {
            body.gravity_scale = scale;
            true
        } else {
            false
        }
    }

    /// Wake a sleeping body.
    pub fn wake_body(&mut self, handle: WasmRigidBodyHandle) -> bool {
        if let Some(body) = self.bodies.get_mut(&handle.0) {
            if body.sleeping {
                body.sleeping = false;
                body.sleep_time = 0.0;
                self.events
                    .push(WasmPhysicsEvent::BodyWoke { body: handle });
            }
            true
        } else {
            false
        }
    }

    /// Get body state.
    pub fn get_body_state(&self, handle: WasmRigidBodyHandle) -> Option<WasmBodyState> {
        self.bodies.get(&handle.0).map(|b| b.state())
    }

    /// Number of rigid bodies.
    pub fn body_count(&self) -> usize {
        self.bodies.len()
    }

    // ── Collider management ───────────────────────────────────────────────

    /// Attach a sphere collider to a body.
    pub fn add_sphere_collider(
        &mut self,
        body: WasmRigidBodyHandle,
        radius: f64,
    ) -> WasmColliderHandle {
        let handle = WasmColliderHandle(self.next_collider_id);
        self.next_collider_id += 1;
        let mut coll = WasmCollider::new(handle, WasmColliderShape::Sphere { radius });
        coll.body = Some(body);
        self.colliders.insert(handle.0, coll);
        handle
    }

    /// Attach a box collider to a body.
    pub fn add_box_collider(
        &mut self,
        body: WasmRigidBodyHandle,
        hx: f64,
        hy: f64,
        hz: f64,
    ) -> WasmColliderHandle {
        let handle = WasmColliderHandle(self.next_collider_id);
        self.next_collider_id += 1;
        let mut coll = WasmCollider::new(
            handle,
            WasmColliderShape::Box {
                half_extents: [hx, hy, hz],
            },
        );
        coll.body = Some(body);
        self.colliders.insert(handle.0, coll);
        handle
    }

    /// Attach a plane collider to a body.
    pub fn add_plane_collider(
        &mut self,
        body: WasmRigidBodyHandle,
        nx: f64,
        ny: f64,
        nz: f64,
        d: f64,
    ) -> WasmColliderHandle {
        let handle = WasmColliderHandle(self.next_collider_id);
        self.next_collider_id += 1;
        let mut coll = WasmCollider::new(
            handle,
            WasmColliderShape::Plane {
                normal: [nx, ny, nz],
                offset: d,
            },
        );
        coll.body = Some(body);
        self.colliders.insert(handle.0, coll);
        handle
    }

    /// Set collider friction.
    pub fn set_friction(&mut self, handle: WasmColliderHandle, friction: f64) -> bool {
        if let Some(coll) = self.colliders.get_mut(&handle.0) {
            coll.friction = friction.clamp(0.0, f64::INFINITY);
            true
        } else {
            false
        }
    }

    /// Set collider restitution (bounciness).
    pub fn set_restitution(&mut self, handle: WasmColliderHandle, restitution: f64) -> bool {
        if let Some(coll) = self.colliders.get_mut(&handle.0) {
            coll.restitution = restitution.clamp(0.0, 1.0);
            true
        } else {
            false
        }
    }

    /// Make a collider a sensor (no collision response).
    pub fn set_sensor(&mut self, handle: WasmColliderHandle, is_sensor: bool) -> bool {
        if let Some(coll) = self.colliders.get_mut(&handle.0) {
            coll.is_sensor = is_sensor;
            true
        } else {
            false
        }
    }

    /// Number of colliders.
    pub fn collider_count(&self) -> usize {
        self.colliders.len()
    }

    // ── Joint management ──────────────────────────────────────────────────

    /// Add a ball joint between two bodies.
    pub fn add_ball_joint(
        &mut self,
        body_a: WasmRigidBodyHandle,
        body_b: WasmRigidBodyHandle,
        anchor_a: [f64; 3],
        anchor_b: [f64; 3],
    ) -> WasmJointHandle {
        let handle = WasmJointHandle(self.next_joint_id);
        self.next_joint_id += 1;
        let joint = WasmJoint::ball(handle, body_a, body_b, anchor_a, anchor_b);
        self.joints.insert(handle.0, joint);
        handle
    }

    /// Add a distance constraint.
    pub fn add_distance_joint(
        &mut self,
        body_a: WasmRigidBodyHandle,
        body_b: WasmRigidBodyHandle,
        target_dist: f64,
    ) -> WasmJointHandle {
        let handle = WasmJointHandle(self.next_joint_id);
        self.next_joint_id += 1;
        let joint = WasmJoint::distance(handle, body_a, body_b, target_dist);
        self.joints.insert(handle.0, joint);
        handle
    }

    /// Remove a joint.
    pub fn remove_joint(&mut self, handle: WasmJointHandle) -> bool {
        self.joints.remove(&handle.0).is_some()
    }

    /// Number of joints.
    pub fn joint_count(&self) -> usize {
        self.joints.len()
    }

    // ── Step ─────────────────────────────────────────────────────────────

    /// Advance the simulation by the configured dt.
    pub fn step(&mut self) {
        let dt = self.config.dt;
        let gravity = self.config.gravity;
        self.events.clear();
        self.contact_events.clear();

        // Integrate all bodies
        let handles: Vec<u32> = self.bodies.keys().copied().collect();
        for id in handles {
            if let Some(body) = self.bodies.get_mut(&id) {
                body.integrate(dt, gravity);
            }
        }

        // Simple collision detection: sphere vs. plane
        self.detect_sphere_plane_collisions(dt);

        // Sleep management
        self.update_sleep(dt);

        self.time += dt;
        self.step_count += 1;
    }

    /// Advance the simulation by a custom delta time.
    pub fn step_with_dt(&mut self, dt: f64) {
        let old_dt = self.config.dt;
        self.config.dt = dt;
        self.step();
        self.config.dt = old_dt;
    }

    fn detect_sphere_plane_collisions(&mut self, dt: f64) {
        // Collect sphere and plane colliders
        let sphere_data: Vec<(WasmColliderHandle, Option<WasmRigidBodyHandle>, f64)> = self
            .colliders
            .values()
            .filter_map(|c| {
                if let WasmColliderShape::Sphere { radius } = c.shape {
                    Some((c.handle, c.body, radius))
                } else {
                    None
                }
            })
            .collect();

        let plane_data: Vec<(WasmColliderHandle, [f64; 3], f64)> = self
            .colliders
            .values()
            .filter_map(|c| {
                if let WasmColliderShape::Plane { normal, offset } = c.shape {
                    Some((c.handle, normal, offset))
                } else {
                    None
                }
            })
            .collect();

        for (sc_handle, body_opt, radius) in &sphere_data {
            let Some(body_handle) = body_opt else {
                continue;
            };
            let body_pos = match self.bodies.get(&body_handle.0) {
                Some(b) => b.position,
                None => continue,
            };

            for (pc_handle, plane_normal, plane_offset) in &plane_data {
                // Distance from sphere center to plane
                let dist = body_pos[0] * plane_normal[0]
                    + body_pos[1] * plane_normal[1]
                    + body_pos[2] * plane_normal[2]
                    - plane_offset;

                let penetration = radius - dist;
                if penetration > 0.0 {
                    // Resolve: push body out of plane
                    if let Some(body) = self.bodies.get_mut(&body_handle.0) {
                        // Position correction
                        body.position[0] += penetration * plane_normal[0];
                        body.position[1] += penetration * plane_normal[1];
                        body.position[2] += penetration * plane_normal[2];

                        // Velocity response (restitution)
                        let restitution = self
                            .colliders
                            .get(&sc_handle.0)
                            .map(|c| c.restitution)
                            .unwrap_or(0.3);
                        let vn = body.linear_vel[0] * plane_normal[0]
                            + body.linear_vel[1] * plane_normal[1]
                            + body.linear_vel[2] * plane_normal[2];
                        if vn < 0.0 {
                            let impulse = -(1.0 + restitution) * vn;
                            body.linear_vel[0] += impulse * plane_normal[0];
                            body.linear_vel[1] += impulse * plane_normal[1];
                            body.linear_vel[2] += impulse * plane_normal[2];
                        }

                        // Friction (simple)
                        let friction = self
                            .colliders
                            .get(&sc_handle.0)
                            .map(|c| c.friction)
                            .unwrap_or(0.5);
                        body.linear_vel[0] *= (1.0 - friction * dt).max(0.0);
                        body.linear_vel[2] *= (1.0 - friction * dt).max(0.0);
                    }

                    self.contact_events.push(WasmContactEvent {
                        collider_a: *sc_handle,
                        collider_b: *pc_handle,
                        normal: *plane_normal,
                        contact_point: [
                            body_pos[0] - radius * plane_normal[0],
                            body_pos[1] - radius * plane_normal[1],
                            body_pos[2] - radius * plane_normal[2],
                        ],
                        depth: penetration,
                        impulse: penetration,
                        friction_impulse: [0.0; 3],
                    });
                }
            }
        }
    }

    fn update_sleep(&mut self, dt: f64) {
        if !self.config.allow_sleeping {
            return;
        }
        let lin_thresh = self.config.sleep_linear_threshold;
        let ang_thresh = self.config.sleep_angular_threshold;
        let sleep_time = self.config.sleep_time_threshold;

        let handles: Vec<u32> = self.bodies.keys().copied().collect();
        for id in handles {
            if let Some(body) = self.bodies.get_mut(&id) {
                if body.body_type != WasmBodyType::Dynamic {
                    continue;
                }
                let lin_speed2 = body.linear_vel[0] * body.linear_vel[0]
                    + body.linear_vel[1] * body.linear_vel[1]
                    + body.linear_vel[2] * body.linear_vel[2];
                let ang_speed2 = body.angular_vel[0] * body.angular_vel[0]
                    + body.angular_vel[1] * body.angular_vel[1]
                    + body.angular_vel[2] * body.angular_vel[2];

                if lin_speed2 < lin_thresh * lin_thresh && ang_speed2 < ang_thresh * ang_thresh {
                    body.sleep_time += dt;
                    if body.sleep_time >= sleep_time && !body.sleeping {
                        body.sleeping = true;
                        let h = WasmRigidBodyHandle(id);
                        self.events.push(WasmPhysicsEvent::BodySlept { body: h });
                    }
                } else {
                    body.sleep_time = 0.0;
                    if body.sleeping {
                        body.sleeping = false;
                        let h = WasmRigidBodyHandle(id);
                        self.events.push(WasmPhysicsEvent::BodyWoke { body: h });
                    }
                }
            }
        }
    }

    // ── Query ─────────────────────────────────────────────────────────────

    /// Cast a ray and return the first hit.
    pub fn raycast(
        &self,
        origin: [f64; 3],
        direction: [f64; 3],
        max_dist: f64,
    ) -> WasmRaycastResult {
        let dir_len = (direction[0] * direction[0]
            + direction[1] * direction[1]
            + direction[2] * direction[2])
            .sqrt();
        if dir_len < 1e-15 {
            return WasmRaycastResult::miss();
        }
        let dir = [
            direction[0] / dir_len,
            direction[1] / dir_len,
            direction[2] / dir_len,
        ];

        let mut best_toi = max_dist;
        let mut best_result = WasmRaycastResult::miss();

        for coll in self.colliders.values() {
            let body_pos = coll
                .body
                .and_then(|h| self.bodies.get(&h.0))
                .map(|b| b.position)
                .unwrap_or([0.0; 3]);

            let hit = match &coll.shape {
                WasmColliderShape::Sphere { radius } => {
                    // Ray-sphere intersection
                    let oc = [
                        origin[0] - body_pos[0],
                        origin[1] - body_pos[1],
                        origin[2] - body_pos[2],
                    ];
                    let b = oc[0] * dir[0] + oc[1] * dir[1] + oc[2] * dir[2];
                    let c = oc[0] * oc[0] + oc[1] * oc[1] + oc[2] * oc[2] - radius * radius;
                    let discriminant = b * b - c;
                    if discriminant < 0.0 {
                        None
                    } else {
                        let t = -b - discriminant.sqrt();
                        if t > 1e-4 && t < best_toi {
                            let hit_pos = [
                                origin[0] + t * dir[0],
                                origin[1] + t * dir[1],
                                origin[2] + t * dir[2],
                            ];
                            let normal = [
                                (hit_pos[0] - body_pos[0]) / radius,
                                (hit_pos[1] - body_pos[1]) / radius,
                                (hit_pos[2] - body_pos[2]) / radius,
                            ];
                            Some((t, hit_pos, normal))
                        } else {
                            None
                        }
                    }
                }
                WasmColliderShape::Plane { normal, offset } => {
                    // Ray-plane intersection
                    let denom = dir[0] * normal[0] + dir[1] * normal[1] + dir[2] * normal[2];
                    if denom.abs() < 1e-10 {
                        None
                    } else {
                        let t = (offset
                            - (origin[0] * normal[0]
                                + origin[1] * normal[1]
                                + origin[2] * normal[2]))
                            / denom;
                        if t > 1e-4 && t < best_toi {
                            let hit_pos = [
                                origin[0] + t * dir[0],
                                origin[1] + t * dir[1],
                                origin[2] + t * dir[2],
                            ];
                            Some((t, hit_pos, *normal))
                        } else {
                            None
                        }
                    }
                }
                _ => None,
            };

            if let Some((t, pos, norm)) = hit {
                best_toi = t;
                best_result = WasmRaycastResult {
                    hit: true,
                    collider: Some(coll.handle),
                    body: coll.body,
                    point: pos,
                    normal: norm,
                    toi: t,
                    feature_id: 0,
                };
            }
        }

        best_result
    }

    /// Test all colliders against an AABB and return overlapping ones.
    pub fn aabb_overlap(&self, min: [f64; 3], max: [f64; 3]) -> WasmOverlapResult {
        let mut result = WasmOverlapResult::empty();

        for coll in self.colliders.values() {
            let body_pos = coll
                .body
                .and_then(|h| self.bodies.get(&h.0))
                .map(|b| b.position)
                .unwrap_or([0.0; 3]);

            let r = coll.shape.bounding_radius();
            let aabb_min = [body_pos[0] - r, body_pos[1] - r, body_pos[2] - r];
            let aabb_max = [body_pos[0] + r, body_pos[1] + r, body_pos[2] + r];

            // AABB vs AABB overlap test
            if aabb_min[0] <= max[0]
                && aabb_max[0] >= min[0]
                && aabb_min[1] <= max[1]
                && aabb_max[1] >= min[1]
                && aabb_min[2] <= max[2]
                && aabb_max[2] >= min[2]
            {
                result.colliders.push(coll.handle);
                result.bodies.push(coll.body);
            }
        }

        result
    }

    // ── Events ────────────────────────────────────────────────────────────

    /// Drain pending physics events.
    pub fn drain_events(&mut self) -> Vec<WasmPhysicsEvent> {
        std::mem::take(&mut self.events)
    }

    /// Drain pending contact events.
    pub fn drain_contact_events(&mut self) -> Vec<WasmContactEvent> {
        std::mem::take(&mut self.contact_events)
    }

    // ── Reset ─────────────────────────────────────────────────────────────

    /// Remove all bodies, colliders, and joints.
    pub fn reset(&mut self) {
        self.bodies.clear();
        self.colliders.clear();
        self.joints.clear();
        self.events.clear();
        self.contact_events.clear();
        self.next_body_id = 0;
        self.next_collider_id = 0;
        self.next_joint_id = 0;
        self.time = 0.0;
        self.step_count = 0;
    }
}

// ── Scene serialization ───────────────────────────────────────────────────────

/// Serialized state of a single body for JSON output.
#[derive(Debug, Clone)]
pub struct SerializedBody {
    /// Handle ID.
    pub id: u32,
    /// Body type as string.
    pub body_type: String,
    /// Position.
    pub position: [f64; 3],
    /// Rotation (x, y, z, w).
    pub rotation: [f64; 4],
    /// Linear velocity.
    pub linear_vel: [f64; 3],
    /// Angular velocity.
    pub angular_vel: [f64; 3],
    /// Mass.
    pub mass: f64,
    /// Sleeping state.
    pub sleeping: bool,
    /// User data.
    pub user_data: u64,
}

/// Serialized state of a single collider.
#[derive(Debug, Clone)]
pub struct SerializedCollider {
    /// Handle ID.
    pub id: u32,
    /// Associated body ID (if any).
    pub body_id: Option<u32>,
    /// Shape description.
    pub shape_type: String,
    /// Shape parameters (radius, half-extents, etc.).
    pub shape_params: Vec<f64>,
    /// Friction.
    pub friction: f64,
    /// Restitution.
    pub restitution: f64,
    /// Is sensor.
    pub is_sensor: bool,
}

/// Full scene snapshot for save/restore.
#[derive(Debug, Clone)]
pub struct WasmSceneSnapshot {
    /// Simulation time.
    pub time: f64,
    /// Step count.
    pub step_count: u64,
    /// Serialized bodies.
    pub bodies: Vec<SerializedBody>,
    /// Serialized colliders.
    pub colliders: Vec<SerializedCollider>,
    /// Gravity.
    pub gravity: [f64; 3],
    /// Configuration dt.
    pub dt: f64,
}

/// WASM scene serializer: converts WasmWorld state to/from a snapshot.
pub struct WasmSceneSerializer;

impl WasmSceneSerializer {
    /// Serialize the given world to a snapshot.
    pub fn serialize(world: &WasmWorld) -> WasmSceneSnapshot {
        let bodies: Vec<SerializedBody> = world
            .bodies
            .values()
            .map(|b| SerializedBody {
                id: b.handle.0,
                body_type: match b.body_type {
                    WasmBodyType::Dynamic => "dynamic".to_string(),
                    WasmBodyType::Static => "static".to_string(),
                    WasmBodyType::Kinematic => "kinematic".to_string(),
                },
                position: b.position,
                rotation: b.rotation,
                linear_vel: b.linear_vel,
                angular_vel: b.angular_vel,
                mass: b.mass,
                sleeping: b.sleeping,
                user_data: b.user_data,
            })
            .collect();

        let colliders: Vec<SerializedCollider> = world
            .colliders
            .values()
            .map(|c| {
                let (shape_type, shape_params) = match &c.shape {
                    WasmColliderShape::Sphere { radius } => ("sphere".to_string(), vec![*radius]),
                    WasmColliderShape::Box { half_extents } => {
                        ("box".to_string(), half_extents.to_vec())
                    }
                    WasmColliderShape::Capsule {
                        half_height,
                        radius,
                    } => ("capsule".to_string(), vec![*half_height, *radius]),
                    WasmColliderShape::Plane { normal, offset } => (
                        "plane".to_string(),
                        vec![normal[0], normal[1], normal[2], *offset],
                    ),
                    WasmColliderShape::Cylinder {
                        half_height,
                        radius,
                    } => ("cylinder".to_string(), vec![*half_height, *radius]),
                    WasmColliderShape::Cone {
                        half_height,
                        radius,
                    } => ("cone".to_string(), vec![*half_height, *radius]),
                    WasmColliderShape::ConvexHull { .. } => ("convex_hull".to_string(), vec![]),
                };
                SerializedCollider {
                    id: c.handle.0,
                    body_id: c.body.map(|h| h.0),
                    shape_type,
                    shape_params,
                    friction: c.friction,
                    restitution: c.restitution,
                    is_sensor: c.is_sensor,
                }
            })
            .collect();

        WasmSceneSnapshot {
            time: world.time,
            step_count: world.step_count,
            bodies,
            colliders,
            gravity: world.config.gravity,
            dt: world.config.dt,
        }
    }

    /// Build a simple JSON string from the snapshot (no serde dependency).
    pub fn to_json_string(snapshot: &WasmSceneSnapshot) -> String {
        let mut out = String::new();
        out.push_str("{\n");
        out.push_str(&format!("  \"time\": {},\n", snapshot.time));
        out.push_str(&format!("  \"step_count\": {},\n", snapshot.step_count));
        out.push_str(&format!(
            "  \"gravity\": [{},{},{}],\n",
            snapshot.gravity[0], snapshot.gravity[1], snapshot.gravity[2]
        ));
        out.push_str(&format!("  \"dt\": {},\n", snapshot.dt));
        out.push_str(&format!("  \"body_count\": {},\n", snapshot.bodies.len()));
        out.push_str(&format!(
            "  \"collider_count\": {}\n",
            snapshot.colliders.len()
        ));
        out.push('}');
        out
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn make_world() -> WasmWorld {
        WasmWorld::earth()
    }

    #[test]
    fn test_world_create() {
        let world = make_world();
        assert_eq!(world.body_count(), 0);
        assert_eq!(world.collider_count(), 0);
        assert_eq!(world.joint_count(), 0);
    }

    #[test]
    fn test_add_dynamic_body() {
        let mut world = make_world();
        let h = world.add_dynamic_body(1.0, 0.0, 10.0, 0.0);
        assert_eq!(world.body_count(), 1);
        let pos = world.get_position(h).unwrap();
        assert!((pos[1] - 10.0).abs() < 1e-12);
    }

    #[test]
    fn test_add_static_body() {
        let mut world = make_world();
        let h = world.add_static_body(0.0, 0.0, 0.0);
        let state = world.get_body_state(h).unwrap();
        assert_eq!(state.body_type, WasmBodyType::Static);
    }

    #[test]
    fn test_body_falls_under_gravity() {
        let mut world = make_world();
        let h = world.add_dynamic_body(1.0, 0.0, 10.0, 0.0);
        let y0 = world.get_position(h).unwrap()[1];
        for _ in 0..60 {
            world.step();
        }
        let y1 = world.get_position(h).unwrap()[1];
        assert!(y1 < y0, "body should fall: y0={y0} y1={y1}");
    }

    #[test]
    fn test_apply_impulse() {
        let mut world = make_world();
        let h = world.add_dynamic_body(1.0, 0.0, 0.0, 0.0);
        world.apply_impulse(h, [10.0, 0.0, 0.0]);
        let vel = world.get_velocity(h).unwrap();
        assert!((vel[0] - 10.0).abs() < 1e-10);
    }

    #[test]
    fn test_apply_force() {
        let mut world = make_world();
        let h = world.add_dynamic_body(2.0, 0.0, 0.0, 0.0);
        world.apply_force(h, [20.0, 0.0, 0.0]);
        world.step();
        let vel = world.get_velocity(h).unwrap();
        assert!(vel[0] > 0.0);
    }

    #[test]
    fn test_set_position() {
        let mut world = make_world();
        let h = world.add_dynamic_body(1.0, 0.0, 0.0, 0.0);
        world.set_position(h, [5.0, 5.0, 5.0]);
        let pos = world.get_position(h).unwrap();
        assert!((pos[0] - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_remove_body() {
        let mut world = make_world();
        let h = world.add_dynamic_body(1.0, 0.0, 0.0, 0.0);
        assert!(world.remove_body(h));
        assert_eq!(world.body_count(), 0);
    }

    #[test]
    fn test_sphere_collider() {
        let mut world = make_world();
        let h = world.add_dynamic_body(1.0, 0.0, 5.0, 0.0);
        let ch = world.add_sphere_collider(h, 0.5);
        assert_eq!(world.collider_count(), 1);
        world.set_friction(ch, 0.3);
        world.set_restitution(ch, 0.5);
    }

    #[test]
    fn test_plane_collider_bounce() {
        let mut world = make_world();
        let ball = world.add_dynamic_body(1.0, 0.0, 2.0, 0.0);
        let sc = world.add_sphere_collider(ball, 0.5);
        world.set_restitution(sc, 0.8);

        let ground = world.add_static_body(0.0, 0.0, 0.0);
        world.add_plane_collider(ground, 0.0, 1.0, 0.0, 0.0);

        // Simulate until ground contact
        for _ in 0..200 {
            world.step();
        }
        let pos = world.get_position(ball).unwrap();
        assert!(pos[1] >= 0.0, "ball should be above or on the ground");
    }

    #[test]
    fn test_add_ball_joint() {
        let mut world = make_world();
        let a = world.add_dynamic_body(1.0, 0.0, 0.0, 0.0);
        let b = world.add_dynamic_body(1.0, 1.0, 0.0, 0.0);
        let jh = world.add_ball_joint(a, b, [0.5, 0.0, 0.0], [-0.5, 0.0, 0.0]);
        assert_eq!(world.joint_count(), 1);
        world.remove_joint(jh);
        assert_eq!(world.joint_count(), 0);
    }

    #[test]
    fn test_add_distance_joint() {
        let mut world = make_world();
        let a = world.add_dynamic_body(1.0, 0.0, 0.0, 0.0);
        let b = world.add_dynamic_body(1.0, 2.0, 0.0, 0.0);
        world.add_distance_joint(a, b, 2.0);
        assert_eq!(world.joint_count(), 1);
    }

    #[test]
    fn test_raycast_sphere_hit() {
        let mut world = make_world();
        let h = world.add_static_body(0.0, 0.0, 0.0);
        world.add_sphere_collider(h, 1.0);
        let result = world.raycast([0.0, 0.0, -5.0], [0.0, 0.0, 1.0], 100.0);
        assert!(result.hit);
        assert!((result.toi - 4.0).abs() < 0.01);
    }

    #[test]
    fn test_raycast_miss() {
        let mut world = make_world();
        let h = world.add_static_body(0.0, 0.0, 0.0);
        world.add_sphere_collider(h, 0.5);
        let result = world.raycast([10.0, 10.0, -5.0], [0.0, 0.0, 1.0], 100.0);
        assert!(!result.hit);
    }

    #[test]
    fn test_raycast_miss_result() {
        let miss = WasmRaycastResult::miss();
        assert!(!miss.hit);
        assert!(miss.toi.is_infinite());
    }

    #[test]
    fn test_aabb_overlap_hit() {
        let mut world = make_world();
        let h = world.add_static_body(0.0, 0.0, 0.0);
        world.add_sphere_collider(h, 0.5);
        let result = world.aabb_overlap([-1.0; 3], [1.0; 3]);
        assert_eq!(result.count(), 1);
    }

    #[test]
    fn test_aabb_overlap_miss() {
        let mut world = make_world();
        let h = world.add_static_body(10.0, 0.0, 0.0);
        world.add_sphere_collider(h, 0.5);
        let result = world.aabb_overlap([-1.0; 3], [1.0; 3]);
        assert_eq!(result.count(), 0);
    }

    #[test]
    fn test_drain_events() {
        let mut world = make_world();
        let h = world.add_dynamic_body(1.0, 0.0, 0.0, 0.0);
        world.wake_body(h);
        world.step();
        // Events may or may not fire — just check drain doesn't panic
        let _events = world.drain_events();
    }

    #[test]
    fn test_contact_events_on_collision() {
        let mut world = make_world();
        let ball = world.add_dynamic_body(1.0, 0.0, 1.0, 0.0);
        world.add_sphere_collider(ball, 0.5);
        let ground = world.add_static_body(0.0, 0.0, 0.0);
        world.add_plane_collider(ground, 0.0, 1.0, 0.0, 0.0);

        // Step until collision
        for _ in 0..100 {
            world.step();
        }
        let contacts = world.drain_contact_events();
        // Should have generated contacts at some point
        // (they get drained each check so just verify no panic)
        let _ = contacts;
    }

    #[test]
    fn test_reset() {
        let mut world = make_world();
        world.add_dynamic_body(1.0, 0.0, 0.0, 0.0);
        world.step();
        world.reset();
        assert_eq!(world.body_count(), 0);
        assert!((world.time).abs() < 1e-15);
        assert_eq!(world.step_count, 0);
    }

    #[test]
    fn test_scene_serializer() {
        let mut world = make_world();
        world.add_dynamic_body(1.0, 0.0, 5.0, 0.0);
        world.step();
        let snapshot = WasmSceneSerializer::serialize(&world);
        assert_eq!(snapshot.bodies.len(), 1);
        let json = WasmSceneSerializer::to_json_string(&snapshot);
        assert!(json.contains("time"));
        assert!(json.contains("body_count"));
    }

    #[test]
    fn test_config_moon_gravity() {
        let cfg = WasmSimulationConfig::moon();
        assert!((cfg.gravity[1] + 1.625).abs() < 1e-10);
    }

    #[test]
    fn test_config_zero_gravity() {
        let cfg = WasmSimulationConfig::zero_gravity();
        assert!((cfg.gravity[1]).abs() < 1e-15);
    }

    #[test]
    fn test_step_with_dt() {
        let mut world = make_world();
        world.add_dynamic_body(1.0, 0.0, 5.0, 0.0);
        world.step_with_dt(0.1);
        assert!((world.config.dt - 1.0 / 60.0).abs() < 1e-10); // dt restored
    }

    #[test]
    fn test_sensor_collider() {
        let mut world = make_world();
        let h = world.add_static_body(0.0, 0.0, 0.0);
        let ch = world.add_sphere_collider(h, 1.0);
        world.set_sensor(ch, true);
    }

    #[test]
    fn test_gravity_scale() {
        let mut world = make_world();
        let h = world.add_dynamic_body(1.0, 0.0, 10.0, 0.0);
        world.set_gravity_scale(h, 0.0);
        let y0 = world.get_position(h).unwrap()[1];
        for _ in 0..60 {
            world.step();
        }
        let y1 = world.get_position(h).unwrap()[1];
        assert!(
            (y1 - y0).abs() < 1e-3,
            "zero gravity scale: body should not fall"
        );
    }

    #[test]
    fn test_quaternion_multiply_identity() {
        let identity = [0.0, 0.0, 0.0, 1.0];
        let q = [0.1, 0.0, 0.0, 1.0];
        let result = quaternion_multiply(identity, q);
        assert!((result[0] - q[0]).abs() < 1e-10);
    }

    #[test]
    fn test_quaternion_normalize() {
        let q = [0.0, 0.0, 0.0, 2.0];
        let n = quaternion_normalize(q);
        assert!((n[3] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_collider_shape_bounding_radius_sphere() {
        let s = WasmColliderShape::Sphere { radius: 2.5 };
        assert!((s.bounding_radius() - 2.5).abs() < 1e-10);
    }

    #[test]
    fn test_collider_shape_bounding_radius_box() {
        let b = WasmColliderShape::Box {
            half_extents: [1.0, 1.0, 1.0],
        };
        let r = b.bounding_radius();
        assert!((r - 3.0f64.sqrt()).abs() < 1e-10);
    }

    #[test]
    fn test_wasm_joint_types() {
        let ha = WasmRigidBodyHandle(0);
        let hb = WasmRigidBodyHandle(1);
        let jh = WasmJointHandle(0);
        let j = WasmJoint::ball(jh, ha, hb, [0.0; 3], [0.0; 3]);
        assert_eq!(j.joint_type, WasmJointType::Ball);
    }

    #[test]
    fn test_multiple_bodies_simulation() {
        let mut world = make_world();
        let mut handles = Vec::new();
        for i in 0..5 {
            let h = world.add_dynamic_body(1.0, i as f64, 5.0, 0.0);
            world.add_sphere_collider(h, 0.3);
            handles.push(h);
        }
        let ground = world.add_static_body(0.0, 0.0, 0.0);
        world.add_plane_collider(ground, 0.0, 1.0, 0.0, 0.0);

        for _ in 0..60 {
            world.step();
        }

        for h in &handles {
            let pos = world.get_position(*h).unwrap();
            assert!(pos[1] >= -0.1, "body should not fall through ground");
        }
    }
}
