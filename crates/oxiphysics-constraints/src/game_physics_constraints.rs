// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Game-physics oriented constraint systems.
//!
//! Implements breakable joints, ragdoll hierarchies, motor constraints,
//! spring/rope/wheel constraints, character controllers, trigger volumes,
//! composite constraint chains, and a sequential-impulse game solver.

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Helper math — [f64; 3] arrays only, no nalgebra
// ---------------------------------------------------------------------------

/// Subtract two 3-vectors (a − b).
#[inline]
fn v3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Scale a 3-vector.
#[inline]
fn v3_scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

/// Dot product of two 3-vectors.
#[inline]
fn v3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Euclidean length of a 3-vector.
#[inline]
fn v3_len(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

/// Normalize a 3-vector, returning zero if near-degenerate.
#[inline]
fn v3_norm(a: [f64; 3]) -> [f64; 3] {
    let l = v3_len(a);
    if l < 1e-15 {
        [0.0; 3]
    } else {
        v3_scale(a, 1.0 / l)
    }
}

/// Clamp `x` to `[lo, hi]`.
#[inline]
fn clampf(x: f64, lo: f64, hi: f64) -> f64 {
    x.max(lo).min(hi)
}

// ---------------------------------------------------------------------------
// BreakableJoint
// ---------------------------------------------------------------------------

/// State of a breakable joint.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum JointState {
    /// Joint is intact and applying constraints normally.
    Intact,
    /// Joint has fractured and is no longer constraining.
    Broken,
}

/// A joint that fractures when the accumulated impulse exceeds a threshold.
///
/// When `accumulated_impulse > break_threshold`, the joint transitions from
/// `Intact` to `Broken` and no longer contributes any constraint force.
#[derive(Debug, Clone)]
pub struct BreakableJoint {
    /// Index of body A.
    pub body_a: usize,
    /// Index of body B.
    pub body_b: usize,
    /// Anchor point relative to body A (local space).
    pub anchor_a: [f64; 3],
    /// Anchor point relative to body B (local space).
    pub anchor_b: [f64; 3],
    /// Maximum impulse the joint can withstand before breaking.
    pub break_threshold: f64,
    /// Accumulated impulse magnitude since last reset.
    pub accumulated_impulse: f64,
    /// Current state of the joint.
    pub state: JointState,
    /// Impulse magnitude at the moment of fracture (zero if intact).
    pub break_impulse: f64,
    /// Rest length of the joint anchor offset.
    pub rest_length: f64,
    /// Stiffness coefficient for compliance model.
    pub stiffness: f64,
    /// Damping coefficient.
    pub damping: f64,
}

impl BreakableJoint {
    /// Create a new intact breakable joint.
    pub fn new(
        body_a: usize,
        body_b: usize,
        anchor_a: [f64; 3],
        anchor_b: [f64; 3],
        break_threshold: f64,
        stiffness: f64,
        damping: f64,
    ) -> Self {
        let rest_length = v3_len(v3_sub(anchor_b, anchor_a));
        Self {
            body_a,
            body_b,
            anchor_a,
            anchor_b,
            break_threshold,
            accumulated_impulse: 0.0,
            state: JointState::Intact,
            break_impulse: 0.0,
            rest_length,
            stiffness,
            damping,
        }
    }

    /// Apply a single impulse step. Checks the break condition.
    ///
    /// Returns the impulse applied this step (0 if broken).
    pub fn apply_impulse(&mut self, impulse: f64) -> f64 {
        if self.state == JointState::Broken {
            return 0.0;
        }
        self.accumulated_impulse += impulse.abs();
        if self.accumulated_impulse >= self.break_threshold {
            self.break_impulse = self.accumulated_impulse;
            self.state = JointState::Broken;
            return 0.0;
        }
        impulse
    }

    /// Reset accumulated impulse for a new time step.
    pub fn reset_accumulated(&mut self) {
        self.accumulated_impulse = 0.0;
    }

    /// Returns `true` if the joint is broken.
    pub fn is_broken(&self) -> bool {
        self.state == JointState::Broken
    }

    /// Compute the spring force magnitude for a given separation distance.
    ///
    /// Uses a damped Hooke model: `F = stiffness * (dist - rest) + damping * vel_along`.
    pub fn spring_force(&self, dist: f64, vel_along: f64) -> f64 {
        if self.state == JointState::Broken {
            return 0.0;
        }
        self.stiffness * (dist - self.rest_length) + self.damping * vel_along
    }
}

// ---------------------------------------------------------------------------
// RagdollConstraints
// ---------------------------------------------------------------------------

/// A single bone segment in a ragdoll hierarchy.
#[derive(Debug, Clone)]
pub struct RagdollBone {
    /// Human-readable bone name (e.g., "left_forearm").
    pub name: String,
    /// Parent bone index (None for root).
    pub parent: Option<usize>,
    /// Rest pose orientation as a unit quaternion `[x, y, z, w]`.
    pub rest_orientation: [f64; 4],
    /// Current world-space position.
    pub position: [f64; 3],
    /// Current orientation as a unit quaternion.
    pub orientation: [f64; 4],
    /// Half-angle swing limits `[cone_half_x, cone_half_y]` in radians.
    pub swing_limits: [f64; 2],
    /// Twist limit half-angle in radians.
    pub twist_limit: f64,
    /// Bone mass in kilograms.
    pub mass: f64,
    /// Half-extents of the bone's collision capsule.
    pub half_length: f64,
    /// Bone radius for capsule shape.
    pub radius: f64,
}

impl RagdollBone {
    /// Create a new ragdoll bone with default physics parameters.
    pub fn new(
        name: impl Into<String>,
        parent: Option<usize>,
        mass: f64,
        half_length: f64,
        radius: f64,
    ) -> Self {
        Self {
            name: name.into(),
            parent,
            rest_orientation: [0.0, 0.0, 0.0, 1.0],
            position: [0.0; 3],
            orientation: [0.0, 0.0, 0.0, 1.0],
            swing_limits: [PI / 4.0, PI / 4.0],
            twist_limit: PI / 6.0,
            mass,
            half_length,
            radius,
        }
    }

    /// Set anatomical swing limits.
    pub fn with_swing_limits(mut self, cone_x: f64, cone_y: f64) -> Self {
        self.swing_limits = [cone_x, cone_y];
        self
    }

    /// Set twist limit.
    pub fn with_twist_limit(mut self, twist: f64) -> Self {
        self.twist_limit = twist;
        self
    }
}

/// A ragdoll composed of multiple bones with hierarchical joint constraints.
#[derive(Debug, Clone)]
pub struct RagdollConstraints {
    /// Ordered list of bones (index 0 = root).
    pub bones: Vec<RagdollBone>,
    /// Gravity acceleration vector applied during fall prediction.
    pub gravity: [f64; 3],
    /// Whether the ragdoll is currently active (simulating).
    pub active: bool,
    /// Blend weight between animation pose and physics (0 = full anim, 1 = full physics).
    pub physics_blend: f64,
}

impl RagdollConstraints {
    /// Create an empty ragdoll.
    pub fn new(gravity: [f64; 3]) -> Self {
        Self {
            bones: Vec::new(),
            gravity,
            active: false,
            physics_blend: 1.0,
        }
    }

    /// Add a bone and return its index.
    pub fn add_bone(&mut self, bone: RagdollBone) -> usize {
        let idx = self.bones.len();
        self.bones.push(bone);
        idx
    }

    /// Activate ragdoll physics.
    pub fn activate(&mut self) {
        self.active = true;
    }

    /// Deactivate (return control to animation).
    pub fn deactivate(&mut self) {
        self.active = false;
    }

    /// Predict the landing position of the root bone after `time` seconds,
    /// assuming constant `gravity` and initial velocity `vel`.
    pub fn predict_fall_position(
        &self,
        initial_pos: [f64; 3],
        vel: [f64; 3],
        time: f64,
    ) -> [f64; 3] {
        // p(t) = p0 + v*t + 0.5*g*t^2
        let half_t2 = 0.5 * time * time;
        [
            initial_pos[0] + vel[0] * time + self.gravity[0] * half_t2,
            initial_pos[1] + vel[1] * time + self.gravity[1] * half_t2,
            initial_pos[2] + vel[2] * time + self.gravity[2] * half_t2,
        ]
    }

    /// Check if a joint angle violates swing limits for a given bone.
    ///
    /// `swing_angle` is the half-cone angle in radians; returns the clamped value.
    pub fn clamp_swing(&self, bone_idx: usize, swing_angle: f64) -> f64 {
        if bone_idx >= self.bones.len() {
            return swing_angle;
        }
        let limit = self.bones[bone_idx].swing_limits[0].min(self.bones[bone_idx].swing_limits[1]);
        clampf(swing_angle, -limit, limit)
    }

    /// Total mass of the ragdoll (sum of all bone masses).
    pub fn total_mass(&self) -> f64 {
        self.bones.iter().map(|b| b.mass).sum()
    }

    /// Compute the centre of mass position (mass-weighted average of bone positions).
    pub fn centre_of_mass(&self) -> [f64; 3] {
        let total = self.total_mass();
        if total < 1e-15 {
            return [0.0; 3];
        }
        let mut com = [0.0f64; 3];
        for bone in &self.bones {
            let w = bone.mass / total;
            com[0] += bone.position[0] * w;
            com[1] += bone.position[1] * w;
            com[2] += bone.position[2] * w;
        }
        com
    }
}

// ---------------------------------------------------------------------------
// MotorConstraintGame
// ---------------------------------------------------------------------------

/// Control mode for a game motor constraint.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MotorMode {
    /// Target angular velocity (rad/s).
    VelocityControl,
    /// Target angle (rad).
    PositionControl,
    /// Target RPM (revolutions per minute).
    RpmControl,
    /// Motor disabled (free-spinning).
    Disabled,
}

/// A servo/motor constraint for game physics (wheels, doors, hinges, etc.).
#[derive(Debug, Clone)]
pub struct MotorConstraintGame {
    /// Body index for the driven body.
    pub body: usize,
    /// Rotation axis in world space.
    pub axis: [f64; 3],
    /// Current control mode.
    pub mode: MotorMode,
    /// Target value (angle \[rad\], velocity \[rad/s\], or RPM depending on mode).
    pub target: f64,
    /// Current angle in radians.
    pub current_angle: f64,
    /// Current angular velocity in rad/s.
    pub current_velocity: f64,
    /// Maximum torque the motor can apply (N·m).
    pub max_torque: f64,
    /// Maximum angular velocity (rad/s) — speed limit.
    pub max_velocity: f64,
    /// PD proportional gain.
    pub kp: f64,
    /// PD derivative gain.
    pub kd: f64,
    /// Accumulated torque for warm-starting.
    pub accumulated_torque: f64,
}

impl MotorConstraintGame {
    /// Create a new motor constraint in velocity-control mode.
    pub fn new(body: usize, axis: [f64; 3], max_torque: f64, max_velocity: f64) -> Self {
        Self {
            body,
            axis: v3_norm(axis),
            mode: MotorMode::VelocityControl,
            target: 0.0,
            current_angle: 0.0,
            current_velocity: 0.0,
            max_torque,
            max_velocity,
            kp: 100.0,
            kd: 10.0,
            accumulated_torque: 0.0,
        }
    }

    /// Switch to RPM control and set target RPM.
    pub fn set_rpm(&mut self, rpm: f64) {
        self.mode = MotorMode::RpmControl;
        self.target = rpm;
    }

    /// Convert target RPM to rad/s.
    pub fn target_rad_per_sec(&self) -> f64 {
        match self.mode {
            MotorMode::RpmControl => self.target * (2.0 * PI / 60.0),
            MotorMode::VelocityControl => self.target,
            MotorMode::PositionControl => 0.0,
            MotorMode::Disabled => 0.0,
        }
    }

    /// Compute the torque to apply for the current time step.
    pub fn compute_torque(&self, dt: f64) -> f64 {
        if self.mode == MotorMode::Disabled {
            return 0.0;
        }
        let torque = match self.mode {
            MotorMode::VelocityControl | MotorMode::RpmControl => {
                let target_vel = self.target_rad_per_sec();
                let vel_err = target_vel - self.current_velocity;
                self.kp * vel_err
            }
            MotorMode::PositionControl => {
                let angle_err = self.target - self.current_angle;
                let vel_err = -self.current_velocity;
                self.kp * angle_err + self.kd * vel_err
            }
            MotorMode::Disabled => 0.0,
        };
        let _ = dt;
        clampf(torque, -self.max_torque, self.max_torque)
    }

    /// Integrate the motor state forward by `dt` seconds under the computed torque,
    /// assuming a rigid body with rotational inertia `inertia`.
    pub fn integrate(&mut self, inertia: f64, dt: f64) {
        let torque = self.compute_torque(dt);
        let alpha = torque / inertia.max(1e-15);
        self.current_velocity = clampf(
            self.current_velocity + alpha * dt,
            -self.max_velocity,
            self.max_velocity,
        );
        self.current_angle += self.current_velocity * dt;
        self.accumulated_torque += torque.abs();
    }
}

// ---------------------------------------------------------------------------
// SpringJoint
// ---------------------------------------------------------------------------

/// A damped spring-damper joint constraint.
///
/// Models a 1-D spring between two anchor points. Uses the standard
/// second-order oscillator formulation with natural frequency `ω_n` and
/// damping ratio `ζ`.
#[derive(Debug, Clone)]
pub struct SpringJoint {
    /// Index of body A.
    pub body_a: usize,
    /// Index of body B.
    pub body_b: usize,
    /// Local-space attachment point on body A.
    pub anchor_a: [f64; 3],
    /// Local-space attachment point on body B.
    pub anchor_b: [f64; 3],
    /// Natural (rest) length of the spring.
    pub rest_length: f64,
    /// Spring stiffness constant `k` (N/m).
    pub stiffness: f64,
    /// Damping coefficient `c` (N·s/m).
    pub damping: f64,
    /// Natural angular frequency `ω_n = sqrt(k/m)`.
    pub natural_frequency: f64,
    /// Damping ratio `ζ = c / (2 * sqrt(k * m))`.
    pub damping_ratio: f64,
    /// Equilibrium position offset (world space).
    pub equilibrium_offset: [f64; 3],
    /// Accumulated spring impulse for warm-starting.
    pub accumulated_impulse: f64,
    /// Whether the spring can only pull (not push — rubber-band spring).
    pub tension_only: bool,
}

impl SpringJoint {
    /// Create a spring joint with explicit stiffness and damping.
    pub fn new(
        body_a: usize,
        body_b: usize,
        anchor_a: [f64; 3],
        anchor_b: [f64; 3],
        rest_length: f64,
        stiffness: f64,
        damping: f64,
        mass: f64,
    ) -> Self {
        let natural_frequency = if mass > 1e-15 {
            (stiffness / mass).sqrt()
        } else {
            0.0
        };
        let critical_damping = 2.0 * (stiffness * mass).sqrt();
        let damping_ratio = if critical_damping > 1e-15 {
            damping / critical_damping
        } else {
            0.0
        };
        Self {
            body_a,
            body_b,
            anchor_a,
            anchor_b,
            rest_length,
            stiffness,
            damping,
            natural_frequency,
            damping_ratio,
            equilibrium_offset: [0.0; 3],
            accumulated_impulse: 0.0,
            tension_only: false,
        }
    }

    /// Enable tension-only (rubber-band) mode.
    pub fn set_tension_only(&mut self, enabled: bool) {
        self.tension_only = enabled;
    }

    /// Compute the spring force for current separation distance and velocity.
    pub fn compute_force(&self, dist: f64, vel: f64) -> f64 {
        let extension = dist - self.rest_length;
        if self.tension_only && extension < 0.0 {
            return 0.0;
        }
        -(self.stiffness * extension + self.damping * vel)
    }

    /// Period of oscillation: `T = 2π / ω_n`.
    pub fn period(&self) -> f64 {
        if self.natural_frequency < 1e-15 {
            f64::INFINITY
        } else {
            2.0 * PI / self.natural_frequency
        }
    }

    /// Return `true` if the system is over-damped (`ζ > 1`).
    pub fn is_overdamped(&self) -> bool {
        self.damping_ratio > 1.0
    }

    /// Return `true` if the system is critically damped (`ζ ≈ 1`).
    pub fn is_critically_damped(&self) -> bool {
        (self.damping_ratio - 1.0).abs() < 0.01
    }

    /// Compute the decay envelope amplitude at time `t` for underdamped systems.
    pub fn decay_envelope(&self, t: f64) -> f64 {
        (-self.damping_ratio * self.natural_frequency * t).exp()
    }
}

// ---------------------------------------------------------------------------
// RopeConstraint
// ---------------------------------------------------------------------------

/// A single node in a rope simulation.
#[derive(Debug, Clone)]
pub struct RopeNode {
    /// Current world-space position.
    pub position: [f64; 3],
    /// Current velocity.
    pub velocity: [f64; 3],
    /// Node mass in kilograms.
    pub mass: f64,
    /// Whether this node is pinned (static).
    pub pinned: bool,
}

impl RopeNode {
    /// Create a new rope node at the given position.
    pub fn new(position: [f64; 3], mass: f64) -> Self {
        Self {
            position,
            velocity: [0.0; 3],
            mass,
            pinned: false,
        }
    }

    /// Pin this node (make it static).
    pub fn pin(&mut self) {
        self.pinned = true;
    }
}

/// An n-segment rope with positional constraints, sag, and damping.
#[derive(Debug, Clone)]
pub struct RopeConstraint {
    /// Ordered rope nodes from start to end.
    pub nodes: Vec<RopeNode>,
    /// Rest segment length (uniform for all segments).
    pub segment_rest_length: f64,
    /// Stiffness factor for position constraint projection (XPBD compliance).
    pub stiffness: f64,
    /// Velocity damping factor per segment.
    pub damping: f64,
    /// Gravity applied to each un-pinned node.
    pub gravity: [f64; 3],
    /// Extra sag factor (0 = tight rope, 1 = very saggy).
    pub sag_factor: f64,
    /// Number of solver iterations per step.
    pub solver_iterations: usize,
}

impl RopeConstraint {
    /// Create a straight rope from `start` to `end` with `num_segments` segments.
    pub fn new(
        start: [f64; 3],
        end: [f64; 3],
        num_segments: usize,
        mass_per_node: f64,
        stiffness: f64,
        damping: f64,
        gravity: [f64; 3],
    ) -> Self {
        let n = num_segments + 1;
        let mut nodes = Vec::with_capacity(n);
        let dir = v3_sub(end, start);
        let total_len = v3_len(dir);
        let segment_rest_length = if num_segments > 0 {
            total_len / num_segments as f64
        } else {
            total_len
        };
        for i in 0..n {
            let t = if num_segments > 0 {
                i as f64 / num_segments as f64
            } else {
                0.0
            };
            let pos = [
                start[0] + dir[0] * t,
                start[1] + dir[1] * t,
                start[2] + dir[2] * t,
            ];
            nodes.push(RopeNode::new(pos, mass_per_node));
        }
        // Pin first and last nodes by default
        if let Some(first) = nodes.first_mut() {
            first.pin();
        }
        if let Some(last) = nodes.last_mut() {
            last.pin();
        }
        Self {
            nodes,
            segment_rest_length,
            stiffness,
            damping,
            gravity,
            sag_factor: 0.0,
            solver_iterations: 4,
        }
    }

    /// Number of segments in the rope.
    pub fn num_segments(&self) -> usize {
        self.nodes.len().saturating_sub(1)
    }

    /// Compute the total rope length (sum of current segment lengths).
    pub fn current_length(&self) -> f64 {
        self.nodes
            .windows(2)
            .map(|w| v3_len(v3_sub(w[1].position, w[0].position)))
            .sum()
    }

    /// Apply gravity to all un-pinned nodes.
    pub fn apply_gravity(&mut self, dt: f64) {
        for node in &mut self.nodes {
            if node.pinned {
                continue;
            }
            node.velocity[0] += self.gravity[0] * dt;
            node.velocity[1] += self.gravity[1] * dt;
            node.velocity[2] += self.gravity[2] * dt;
        }
    }

    /// Integrate node positions by one time step.
    pub fn integrate(&mut self, dt: f64) {
        for node in &mut self.nodes {
            if node.pinned {
                continue;
            }
            node.position[0] += node.velocity[0] * dt;
            node.position[1] += node.velocity[1] * dt;
            node.position[2] += node.velocity[2] * dt;
        }
    }

    /// Project all segment length constraints (one Gauss-Seidel pass).
    pub fn project_constraints(&mut self) {
        let rest = self.segment_rest_length * (1.0 + self.sag_factor);
        let n = self.nodes.len();
        for _ in 0..self.solver_iterations {
            for i in 0..n.saturating_sub(1) {
                let pa = self.nodes[i].position;
                let pb = self.nodes[i + 1].position;
                let diff = v3_sub(pb, pa);
                let dist = v3_len(diff);
                if dist < 1e-15 {
                    continue;
                }
                let correction = (dist - rest) / dist;
                let inv_ma = if self.nodes[i].pinned {
                    0.0
                } else {
                    1.0 / self.nodes[i].mass
                };
                let inv_mb = if self.nodes[i + 1].pinned {
                    0.0
                } else {
                    1.0 / self.nodes[i + 1].mass
                };
                let total_inv = inv_ma + inv_mb;
                if total_inv < 1e-15 {
                    continue;
                }
                let scale = self.stiffness * correction / total_inv;
                if !self.nodes[i].pinned {
                    let wa = inv_ma * scale;
                    self.nodes[i].position[0] += diff[0] * wa;
                    self.nodes[i].position[1] += diff[1] * wa;
                    self.nodes[i].position[2] += diff[2] * wa;
                }
                if !self.nodes[i + 1].pinned {
                    let wb = inv_mb * scale;
                    self.nodes[i + 1].position[0] -= diff[0] * wb;
                    self.nodes[i + 1].position[1] -= diff[1] * wb;
                    self.nodes[i + 1].position[2] -= diff[2] * wb;
                }
            }
        }
    }

    /// Step the rope simulation forward by `dt` seconds.
    pub fn step(&mut self, dt: f64) {
        self.apply_gravity(dt);
        self.integrate(dt);
        self.project_constraints();
    }
}

// ---------------------------------------------------------------------------
// WheelConstraint
// ---------------------------------------------------------------------------

/// Ground contact state for a wheel.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WheelContact {
    /// Wheel is in contact with the ground.
    Grounded,
    /// Wheel is airborne.
    Airborne,
}

/// Constraint that models a vehicle wheel on an axle.
#[derive(Debug, Clone)]
pub struct WheelConstraint {
    /// Index of the wheel body.
    pub body: usize,
    /// Spin axis direction (local space).
    pub spin_axis: [f64; 3],
    /// Current angular velocity of the wheel (rad/s).
    pub spin_velocity: f64,
    /// Maximum suspension travel (m).
    pub suspension_travel: f64,
    /// Current suspension compression (0 = unloaded, 1 = fully compressed).
    pub suspension_compression: f64,
    /// Suspension spring stiffness (N/m).
    pub suspension_stiffness: f64,
    /// Suspension damping (N·s/m).
    pub suspension_damping: f64,
    /// Wheel radius (m).
    pub radius: f64,
    /// Camber angle (rad, positive = top outward).
    pub camber: f64,
    /// Toe angle (rad, positive = inward).
    pub toe: f64,
    /// Current ground contact state.
    pub contact: WheelContact,
    /// Brake torque being applied (N·m).
    pub brake_torque: f64,
    /// Drive torque being applied (N·m).
    pub drive_torque: f64,
}

impl WheelConstraint {
    /// Create a wheel constraint with default parameters.
    pub fn new(body: usize, radius: f64, suspension_travel: f64) -> Self {
        Self {
            body,
            spin_axis: [1.0, 0.0, 0.0],
            spin_velocity: 0.0,
            suspension_travel,
            suspension_compression: 0.0,
            suspension_stiffness: 30000.0,
            suspension_damping: 2000.0,
            radius,
            camber: 0.0,
            toe: 0.0,
            contact: WheelContact::Airborne,
            brake_torque: 0.0,
            drive_torque: 0.0,
        }
    }

    /// Compute the suspension force from current compression state.
    ///
    /// Uses a damped spring model. `compression_velocity` is the rate of change
    /// of `suspension_compression` (positive = compressing).
    pub fn suspension_force(&self, compression_velocity: f64) -> f64 {
        if self.contact == WheelContact::Airborne {
            return 0.0;
        }
        let spring_f =
            self.suspension_stiffness * self.suspension_compression * self.suspension_travel;
        let damp_f = self.suspension_damping * compression_velocity;
        (spring_f + damp_f).max(0.0)
    }

    /// Apply brake torque, returning the resulting angular deceleration (rad/s²)
    /// assuming a wheel with rotational inertia `inertia` (kg·m²).
    pub fn apply_brake(&self, inertia: f64) -> f64 {
        if inertia < 1e-15 {
            return 0.0;
        }
        let sign = -self.spin_velocity.signum();
        let decel = sign * self.brake_torque / inertia;
        // Clamp so we don't reverse direction
        if decel.abs() * 1.0 > self.spin_velocity.abs() {
            -self.spin_velocity
        } else {
            decel
        }
    }

    /// Update wheel spin from drive torque, brake, and wheel inertia.
    pub fn integrate_spin(&mut self, inertia: f64, dt: f64) {
        let drive_alpha = if inertia > 1e-15 {
            self.drive_torque / inertia
        } else {
            0.0
        };
        let brake_alpha = self.apply_brake(inertia);
        self.spin_velocity += (drive_alpha + brake_alpha) * dt;
    }

    /// Speed of the contact patch (m/s) based on spin velocity and radius.
    pub fn contact_patch_speed(&self) -> f64 {
        self.spin_velocity * self.radius
    }
}

// ---------------------------------------------------------------------------
// CharacterConstraint
// ---------------------------------------------------------------------------

/// A character controller constraint (ground contact, slope, step-up, gravity override).
#[derive(Debug, Clone)]
pub struct CharacterConstraint {
    /// Body index of the character.
    pub body: usize,
    /// Maximum slope angle (radians) the character can walk on.
    pub max_slope_angle: f64,
    /// Maximum step height (m) the character can step up.
    pub max_step_height: f64,
    /// Whether the character is currently grounded.
    pub grounded: bool,
    /// Current ground normal (world space).
    pub ground_normal: [f64; 3],
    /// Custom gravity override (None = use world gravity).
    pub gravity_override: Option<[f64; 3]>,
    /// Character capsule height.
    pub capsule_height: f64,
    /// Character capsule radius.
    pub capsule_radius: f64,
    /// Current linear velocity.
    pub velocity: [f64; 3],
    /// Desired movement direction (normalized world space).
    pub move_direction: [f64; 3],
    /// Move speed (m/s).
    pub move_speed: f64,
    /// Jump impulse magnitude (m/s).
    pub jump_velocity: f64,
    /// Whether a jump was requested this frame.
    pub jump_requested: bool,
}

impl CharacterConstraint {
    /// Create a new character constraint with typical defaults.
    pub fn new(body: usize, capsule_height: f64, capsule_radius: f64) -> Self {
        Self {
            body,
            max_slope_angle: PI / 4.0,
            max_step_height: 0.4,
            grounded: false,
            ground_normal: [0.0, 1.0, 0.0],
            gravity_override: None,
            capsule_height,
            capsule_radius,
            velocity: [0.0; 3],
            move_direction: [0.0; 3],
            move_speed: 5.0,
            jump_velocity: 6.0,
            jump_requested: false,
        }
    }

    /// Returns `true` if the current ground slope is walkable.
    pub fn is_slope_walkable(&self) -> bool {
        let up = [0.0f64, 1.0, 0.0];
        let cos_angle = v3_dot(self.ground_normal, up);
        let angle = cos_angle.clamp(-1.0, 1.0).acos();
        angle <= self.max_slope_angle
    }

    /// Determine whether the character can step up over an obstacle of height `h`.
    pub fn can_step_up(&self, h: f64) -> bool {
        h <= self.max_step_height
    }

    /// Apply a jump if jump was requested and character is grounded.
    /// Returns the vertical impulse applied (0 if not jumping).
    pub fn try_jump(&mut self) -> f64 {
        if self.grounded && self.jump_requested {
            self.jump_requested = false;
            self.velocity[1] = self.jump_velocity;
            self.jump_velocity
        } else {
            self.jump_requested = false;
            0.0
        }
    }

    /// Project velocity onto the walkable surface plane (remove component into ground).
    pub fn project_velocity_on_ground(&self) -> [f64; 3] {
        let n = self.ground_normal;
        let vn = v3_dot(self.velocity, n);
        if vn >= 0.0 {
            return self.velocity;
        }
        v3_sub(self.velocity, v3_scale(n, vn))
    }
}

// ---------------------------------------------------------------------------
// TriggerVolume
// ---------------------------------------------------------------------------

/// Shape type for a trigger volume.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum TriggerShape {
    /// Axis-aligned bounding box defined by centre and half-extents.
    Aabb {
        /// Centre of the AABB.
        centre: [f64; 3],
        /// Half-extents (width/2, height/2, depth/2).
        half_extents: [f64; 3],
    },
    /// Sphere defined by centre and radius.
    Sphere {
        /// Centre of the sphere.
        centre: [f64; 3],
        /// Radius of the sphere.
        radius: f64,
    },
}

impl TriggerShape {
    /// Test whether a world-space point `p` is inside this trigger.
    pub fn contains(&self, p: [f64; 3]) -> bool {
        match *self {
            TriggerShape::Aabb {
                centre,
                half_extents,
            } => {
                (p[0] - centre[0]).abs() <= half_extents[0]
                    && (p[1] - centre[1]).abs() <= half_extents[1]
                    && (p[2] - centre[2]).abs() <= half_extents[2]
            }
            TriggerShape::Sphere { centre, radius } => v3_len(v3_sub(p, centre)) <= radius,
        }
    }
}

/// An event generated by a trigger volume.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TriggerEvent {
    /// A body entered the trigger this frame.
    Enter,
    /// A body is inside the trigger.
    Stay,
    /// A body exited the trigger this frame.
    Exit,
}

/// A trigger volume that detects body overlap and fires enter/exit events.
#[derive(Debug, Clone)]
pub struct TriggerVolume {
    /// Trigger identifier.
    pub id: u64,
    /// Shape defining the trigger region.
    pub shape: TriggerShape,
    /// Layer mask — only bodies with matching `layer_bit` are detected.
    pub filter_mask: u32,
    /// Set of body indices currently inside the volume.
    pub inside_set: Vec<usize>,
    /// Events generated during the last `update()` call.
    pub events: Vec<(usize, TriggerEvent)>,
    /// Whether the trigger is enabled.
    pub enabled: bool,
}

impl TriggerVolume {
    /// Create a new trigger volume.
    pub fn new(id: u64, shape: TriggerShape, filter_mask: u32) -> Self {
        Self {
            id,
            shape,
            filter_mask,
            inside_set: Vec::new(),
            events: Vec::new(),
            enabled: true,
        }
    }

    /// Update overlap state given the current positions of candidate bodies.
    ///
    /// `bodies` is a slice of `(body_index, position, layer_bit)` tuples.
    /// Returns a slice of the generated events.
    pub fn update(&mut self, bodies: &[(usize, [f64; 3], u32)]) -> &[(usize, TriggerEvent)] {
        self.events.clear();
        if !self.enabled {
            return &self.events;
        }

        let mut new_inside: Vec<usize> = Vec::new();
        for &(idx, pos, layer) in bodies {
            if (layer & self.filter_mask) == 0 {
                continue;
            }
            if self.shape.contains(pos) {
                new_inside.push(idx);
                if !self.inside_set.contains(&idx) {
                    self.events.push((idx, TriggerEvent::Enter));
                } else {
                    self.events.push((idx, TriggerEvent::Stay));
                }
            }
        }
        // Check for exits
        for &old in &self.inside_set {
            if !new_inside.contains(&old) {
                self.events.push((old, TriggerEvent::Exit));
            }
        }
        self.inside_set = new_inside;
        &self.events
    }
}

// ---------------------------------------------------------------------------
// CompositeConstraint
// ---------------------------------------------------------------------------

/// Priority level for constraint solving order.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum ConstraintPriority {
    /// Solved first — hard constraints (non-penetration, joints).
    Critical = 0,
    /// Solved second — gameplay constraints (triggers, motors).
    High = 1,
    /// Solved last — soft/optional constraints.
    Low = 2,
}

/// A descriptor for one sub-constraint in a composite.
#[derive(Debug, Clone)]
pub struct SubConstraintRef {
    /// Human-readable identifier.
    pub name: String,
    /// Solve priority.
    pub priority: ConstraintPriority,
    /// Whether this sub-constraint is currently active.
    pub active: bool,
    /// Warm-start accumulated impulse.
    pub warm_impulse: f64,
    /// Compliance (softness): 0 = hard, larger = softer.
    pub compliance: f64,
}

impl SubConstraintRef {
    /// Create a new sub-constraint reference.
    pub fn new(name: impl Into<String>, priority: ConstraintPriority, compliance: f64) -> Self {
        Self {
            name: name.into(),
            priority,
            active: true,
            warm_impulse: 0.0,
            compliance,
        }
    }
}

/// A composite constraint that chains multiple sub-constraints with priority solving.
#[derive(Debug, Clone)]
pub struct CompositeConstraint {
    /// All sub-constraints in this composite.
    pub constraints: Vec<SubConstraintRef>,
    /// Whether the entire group is active.
    pub group_active: bool,
    /// Number of inner solver iterations for this composite.
    pub iterations: usize,
}

impl CompositeConstraint {
    /// Create an empty composite constraint.
    pub fn new(iterations: usize) -> Self {
        Self {
            constraints: Vec::new(),
            group_active: true,
            iterations,
        }
    }

    /// Add a sub-constraint.
    pub fn add(&mut self, sub: SubConstraintRef) {
        self.constraints.push(sub);
    }

    /// Return constraints sorted by priority (Critical first).
    pub fn sorted_by_priority(&self) -> Vec<&SubConstraintRef> {
        let mut refs: Vec<&SubConstraintRef> = self.constraints.iter().collect();
        refs.sort_by_key(|c| c.priority);
        refs
    }

    /// Activate all sub-constraints.
    pub fn activate_all(&mut self) {
        for c in &mut self.constraints {
            c.active = true;
        }
    }

    /// Deactivate all sub-constraints.
    pub fn deactivate_all(&mut self) {
        for c in &mut self.constraints {
            c.active = false;
        }
    }

    /// Count active sub-constraints.
    pub fn active_count(&self) -> usize {
        self.constraints.iter().filter(|c| c.active).count()
    }
}

// ---------------------------------------------------------------------------
// ConstraintSolver — sequential impulse game solver
// ---------------------------------------------------------------------------

/// State of a single 1-DOF constraint row for the sequential impulse solver.
#[derive(Debug, Clone)]
pub struct ConstraintRow {
    /// Jacobian linear component for body A.
    pub j_lin_a: [f64; 3],
    /// Jacobian angular component for body A.
    pub j_ang_a: [f64; 3],
    /// Jacobian linear component for body B.
    pub j_lin_b: [f64; 3],
    /// Jacobian angular component for body B.
    pub j_ang_b: [f64; 3],
    /// Bias velocity (Baumgarte + restitution).
    pub bias: f64,
    /// Effective mass `1/K`.
    pub eff_mass: f64,
    /// Lower bound for the clamped impulse.
    pub lo: f64,
    /// Upper bound for the clamped impulse.
    pub hi: f64,
    /// Accumulated (warm-start) impulse.
    pub accumulated: f64,
}

impl ConstraintRow {
    /// Create a new constraint row with given Jacobian entries.
    pub fn new(
        j_lin_a: [f64; 3],
        j_ang_a: [f64; 3],
        j_lin_b: [f64; 3],
        j_ang_b: [f64; 3],
        bias: f64,
        eff_mass: f64,
        lo: f64,
        hi: f64,
    ) -> Self {
        Self {
            j_lin_a,
            j_ang_a,
            j_lin_b,
            j_ang_b,
            bias,
            eff_mass,
            lo,
            hi,
            accumulated: 0.0,
        }
    }

    /// Compute the velocity error Cdot for given body velocities.
    pub fn cdot(&self, v_a: [f64; 3], w_a: [f64; 3], v_b: [f64; 3], w_b: [f64; 3]) -> f64 {
        v3_dot(self.j_lin_a, v_a)
            + v3_dot(self.j_ang_a, w_a)
            + v3_dot(self.j_lin_b, v_b)
            + v3_dot(self.j_ang_b, w_b)
    }

    /// Perform one impulse iteration. Returns the delta impulse applied.
    pub fn solve(&mut self, cdot: f64) -> f64 {
        let raw = -self.eff_mass * (cdot + self.bias);
        let old = self.accumulated;
        self.accumulated = clampf(old + raw, self.lo, self.hi);
        self.accumulated - old
    }
}

/// Body velocity state for the sequential impulse solver.
#[derive(Debug, Clone, Copy)]
pub struct SolverBody {
    /// Linear velocity.
    pub vel: [f64; 3],
    /// Angular velocity.
    pub omega: [f64; 3],
    /// Inverse mass (0 = static).
    pub inv_mass: f64,
    /// World-space inverse inertia diagonal (simplified diagonal tensor).
    pub inv_inertia: [f64; 3],
    /// Whether this body is sleeping.
    pub sleeping: bool,
}

impl SolverBody {
    /// Create a dynamic solver body.
    pub fn new(vel: [f64; 3], omega: [f64; 3], inv_mass: f64, inv_inertia: [f64; 3]) -> Self {
        Self {
            vel,
            omega,
            inv_mass,
            inv_inertia,
            sleeping: false,
        }
    }

    /// Create a static (kinematic) body.
    pub fn static_body() -> Self {
        Self {
            vel: [0.0; 3],
            omega: [0.0; 3],
            inv_mass: 0.0,
            inv_inertia: [0.0; 3],
            sleeping: false,
        }
    }

    /// Apply a delta impulse in direction `j` to this body.
    pub fn apply_delta(&mut self, j_lin: [f64; 3], j_ang: [f64; 3], delta: f64, sign: f64) {
        if self.inv_mass < 1e-30 {
            return;
        }
        self.vel[0] += sign * delta * j_lin[0] * self.inv_mass;
        self.vel[1] += sign * delta * j_lin[1] * self.inv_mass;
        self.vel[2] += sign * delta * j_lin[2] * self.inv_mass;
        self.omega[0] += sign * delta * j_ang[0] * self.inv_inertia[0];
        self.omega[1] += sign * delta * j_ang[1] * self.inv_inertia[1];
        self.omega[2] += sign * delta * j_ang[2] * self.inv_inertia[2];
    }

    /// Return the kinetic energy of this body.
    pub fn kinetic_energy(&self) -> f64 {
        let lin_ke = if self.inv_mass > 1e-30 {
            0.5 / self.inv_mass * v3_dot(self.vel, self.vel)
        } else {
            0.0
        };
        let ang_ke = 0.5
            * (self.omega[0] * self.omega[0] / self.inv_inertia[0].max(1e-30)
                + self.omega[1] * self.omega[1] / self.inv_inertia[1].max(1e-30)
                + self.omega[2] * self.omega[2] / self.inv_inertia[2].max(1e-30));
        lin_ke + ang_ke
    }
}

/// A constraint pair binding two solver bodies with a list of constraint rows.
#[derive(Debug, Clone)]
pub struct ConstraintPair {
    /// Index of body A in the bodies array.
    pub body_a: usize,
    /// Index of body B in the bodies array.
    pub body_b: usize,
    /// Constraint rows (one per DOF).
    pub rows: Vec<ConstraintRow>,
}

impl ConstraintPair {
    /// Create a constraint pair.
    pub fn new(body_a: usize, body_b: usize) -> Self {
        Self {
            body_a,
            body_b,
            rows: Vec::new(),
        }
    }

    /// Add a constraint row.
    pub fn add_row(&mut self, row: ConstraintRow) {
        self.rows.push(row);
    }
}

/// Sequential-impulse game constraint solver with warm-starting and sleeping.
#[derive(Debug, Clone)]
pub struct ConstraintSolver {
    /// Number of solver iterations per step.
    pub iterations: usize,
    /// Warm-start scale factor (typically 0.8–1.0).
    pub warm_start_factor: f64,
    /// Sleep threshold: bodies with KE below this are candidates for sleep.
    pub sleep_threshold: f64,
    /// Number of consecutive steps below threshold before sleeping.
    pub sleep_steps: usize,
    /// Internal sleep-step counter per body (indexed same as bodies slice).
    sleep_counters: Vec<usize>,
}

impl ConstraintSolver {
    /// Create a new solver with default parameters.
    pub fn new(iterations: usize) -> Self {
        Self {
            iterations,
            warm_start_factor: 0.85,
            sleep_threshold: 0.01,
            sleep_steps: 60,
            sleep_counters: Vec::new(),
        }
    }

    /// Warm-start all constraint rows by applying the stored accumulated impulses.
    pub fn warm_start(&self, pairs: &[ConstraintPair], bodies: &mut [SolverBody]) {
        for pair in pairs {
            for row in &pair.rows {
                let delta = row.accumulated * self.warm_start_factor;
                if delta.abs() < 1e-30 {
                    continue;
                }
                if pair.body_a < bodies.len() {
                    bodies[pair.body_a].apply_delta(row.j_lin_a, row.j_ang_a, delta, 1.0);
                }
                if pair.body_b < bodies.len() {
                    bodies[pair.body_b].apply_delta(row.j_lin_b, row.j_ang_b, delta, 1.0);
                }
            }
        }
    }

    /// Run the sequential impulse solve loop.
    pub fn solve(&self, pairs: &mut [ConstraintPair], bodies: &mut [SolverBody]) {
        for _ in 0..self.iterations {
            for pair in pairs.iter_mut() {
                let ba = pair.body_a;
                let bb = pair.body_b;
                for row in &mut pair.rows {
                    let (va, wa) = if ba < bodies.len() {
                        (bodies[ba].vel, bodies[ba].omega)
                    } else {
                        ([0.0; 3], [0.0; 3])
                    };
                    let (vb, wb) = if bb < bodies.len() {
                        (bodies[bb].vel, bodies[bb].omega)
                    } else {
                        ([0.0; 3], [0.0; 3])
                    };
                    let cdot = row.cdot(va, wa, vb, wb);
                    let delta = row.solve(cdot);
                    if delta.abs() < 1e-30 {
                        continue;
                    }
                    if ba < bodies.len() {
                        let b = &mut bodies[ba];
                        b.apply_delta(row.j_lin_a, row.j_ang_a, delta, 1.0);
                    }
                    if bb < bodies.len() {
                        let b = &mut bodies[bb];
                        b.apply_delta(row.j_lin_b, row.j_ang_b, delta, 1.0);
                    }
                }
            }
        }
    }

    /// Update sleep state for each body. Bodies sleeping for `sleep_steps`
    /// consecutive steps are flagged `sleeping = true`.
    pub fn update_sleep(&mut self, bodies: &mut [SolverBody]) {
        if self.sleep_counters.len() < bodies.len() {
            self.sleep_counters.resize(bodies.len(), 0);
        }
        for (i, body) in bodies.iter_mut().enumerate() {
            if body.inv_mass < 1e-30 {
                // Static bodies never sleep
                self.sleep_counters[i] = 0;
                continue;
            }
            if body.kinetic_energy() < self.sleep_threshold {
                self.sleep_counters[i] += 1;
                if self.sleep_counters[i] >= self.sleep_steps {
                    body.sleeping = true;
                }
            } else {
                self.sleep_counters[i] = 0;
                body.sleeping = false;
            }
        }
    }

    /// Wake all bodies (reset sleep counters).
    pub fn wake_all(&mut self, bodies: &mut [SolverBody]) {
        for body in bodies.iter_mut() {
            body.sleeping = false;
        }
        for c in &mut self.sleep_counters {
            *c = 0;
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- BreakableJoint ---

    #[test]
    fn test_breakable_joint_intact_below_threshold() {
        let mut j = BreakableJoint::new(0, 1, [0.0; 3], [1.0, 0.0, 0.0], 10.0, 100.0, 5.0);
        let applied = j.apply_impulse(4.0);
        assert!((applied - 4.0).abs() < 1e-12);
        assert_eq!(j.state, JointState::Intact);
    }

    #[test]
    fn test_breakable_joint_breaks_at_threshold() {
        let mut j = BreakableJoint::new(0, 1, [0.0; 3], [1.0, 0.0, 0.0], 5.0, 100.0, 5.0);
        j.apply_impulse(3.0);
        let applied = j.apply_impulse(3.0); // total = 6 > 5 → break
        assert_eq!(applied, 0.0);
        assert_eq!(j.state, JointState::Broken);
        assert!(j.break_impulse >= 5.0);
    }

    #[test]
    fn test_breakable_joint_broken_returns_zero() {
        let mut j = BreakableJoint::new(0, 1, [0.0; 3], [1.0, 0.0, 0.0], 1.0, 100.0, 5.0);
        j.apply_impulse(2.0);
        let a1 = j.apply_impulse(10.0);
        let a2 = j.apply_impulse(10.0);
        assert_eq!(a1, 0.0);
        assert_eq!(a2, 0.0);
    }

    #[test]
    fn test_breakable_joint_spring_force() {
        let j = BreakableJoint::new(0, 1, [0.0; 3], [2.0, 0.0, 0.0], 100.0, 50.0, 3.0);
        // dist=2.0, rest_length=2.0 → extension=0, force = 3.0 * vel_along
        let f = j.spring_force(2.0, 1.0);
        assert!((f - 3.0).abs() < 1e-12, "f={f}");
    }

    #[test]
    fn test_breakable_joint_spring_force_broken() {
        let mut j = BreakableJoint::new(0, 1, [0.0; 3], [1.0, 0.0, 0.0], 1.0, 50.0, 3.0);
        j.apply_impulse(5.0);
        let f = j.spring_force(1.5, 1.0);
        assert_eq!(f, 0.0);
    }

    // --- RagdollConstraints ---

    #[test]
    fn test_ragdoll_total_mass() {
        let mut rag = RagdollConstraints::new([0.0, -9.81, 0.0]);
        rag.add_bone(RagdollBone::new("hip", None, 5.0, 0.2, 0.05));
        rag.add_bone(RagdollBone::new("spine", Some(0), 3.0, 0.25, 0.04));
        assert!((rag.total_mass() - 8.0).abs() < 1e-12);
    }

    #[test]
    fn test_ragdoll_predict_fall_position() {
        let rag = RagdollConstraints::new([0.0, -9.81, 0.0]);
        let pos = rag.predict_fall_position([0.0, 10.0, 0.0], [0.0, 0.0, 0.0], 1.0);
        let expected_y = 10.0 + 0.5 * (-9.81) * 1.0;
        assert!((pos[1] - expected_y).abs() < 1e-6, "pos_y={}", pos[1]);
    }

    #[test]
    fn test_ragdoll_centre_of_mass() {
        let mut rag = RagdollConstraints::new([0.0, -9.81, 0.0]);
        let mut b0 = RagdollBone::new("a", None, 1.0, 0.1, 0.05);
        b0.position = [0.0, 0.0, 0.0];
        let mut b1 = RagdollBone::new("b", Some(0), 1.0, 0.1, 0.05);
        b1.position = [2.0, 0.0, 0.0];
        rag.add_bone(b0);
        rag.add_bone(b1);
        let com = rag.centre_of_mass();
        assert!((com[0] - 1.0).abs() < 1e-12, "com_x={}", com[0]);
    }

    #[test]
    fn test_ragdoll_clamp_swing() {
        let mut rag = RagdollConstraints::new([0.0, -9.81, 0.0]);
        let bone = RagdollBone::new("hip", None, 5.0, 0.2, 0.05).with_swing_limits(0.5, 0.5);
        rag.add_bone(bone);
        let clamped = rag.clamp_swing(0, 1.0);
        assert!((clamped - 0.5).abs() < 1e-12);
    }

    // --- MotorConstraintGame ---

    #[test]
    fn test_motor_rpm_to_rad_per_sec() {
        let mut m = MotorConstraintGame::new(0, [0.0, 0.0, 1.0], 100.0, 50.0);
        m.set_rpm(60.0);
        let rad = m.target_rad_per_sec();
        assert!((rad - 2.0 * PI).abs() < 1e-8, "rad={rad}");
    }

    #[test]
    fn test_motor_velocity_control_torque() {
        let mut m = MotorConstraintGame::new(0, [0.0, 0.0, 1.0], 500.0, 100.0);
        m.mode = MotorMode::VelocityControl;
        m.target = 10.0;
        m.current_velocity = 0.0;
        m.kp = 50.0;
        let t = m.compute_torque(0.01);
        assert!((t - 500.0).abs() < 1e-6, "t={t}"); // clamped at max_torque
    }

    #[test]
    fn test_motor_disabled_no_torque() {
        let mut m = MotorConstraintGame::new(0, [0.0, 0.0, 1.0], 100.0, 50.0);
        m.mode = MotorMode::Disabled;
        assert_eq!(m.compute_torque(0.01), 0.0);
    }

    #[test]
    fn test_motor_integrate_position() {
        let mut m = MotorConstraintGame::new(0, [0.0, 0.0, 1.0], 1000.0, 100.0);
        m.mode = MotorMode::VelocityControl;
        m.target = 1.0;
        m.kp = 10.0;
        m.current_velocity = 0.0;
        let inertia = 1.0;
        let dt = 0.1;
        m.integrate(inertia, dt);
        // After one step velocity should have increased
        assert!(m.current_velocity > 0.0);
    }

    // --- SpringJoint ---

    #[test]
    fn test_spring_force_at_rest() {
        let s = SpringJoint::new(0, 1, [0.0; 3], [1.0, 0.0, 0.0], 1.0, 100.0, 5.0, 1.0);
        let f = s.compute_force(1.0, 0.0);
        assert!(f.abs() < 1e-12, "f={f}");
    }

    #[test]
    fn test_spring_force_extended() {
        let s = SpringJoint::new(0, 1, [0.0; 3], [1.0, 0.0, 0.0], 1.0, 100.0, 0.0, 1.0);
        // Extended by 0.5 m, no damping
        let f = s.compute_force(1.5, 0.0);
        // force = -(100 * 0.5) = -50
        assert!((f + 50.0).abs() < 1e-10, "f={f}");
    }

    #[test]
    fn test_spring_natural_frequency() {
        // k=100, m=1 → ω_n = 10 rad/s
        let s = SpringJoint::new(0, 1, [0.0; 3], [1.0, 0.0, 0.0], 1.0, 100.0, 0.0, 1.0);
        assert!(
            (s.natural_frequency - 10.0).abs() < 1e-8,
            "ωn={}",
            s.natural_frequency
        );
    }

    #[test]
    fn test_spring_period() {
        let s = SpringJoint::new(0, 1, [0.0; 3], [1.0, 0.0, 0.0], 1.0, 100.0, 0.0, 1.0);
        let expected = 2.0 * PI / 10.0;
        assert!((s.period() - expected).abs() < 1e-8);
    }

    #[test]
    fn test_spring_overdamped() {
        // critical_damping = 2*sqrt(100*1) = 20; damping = 30 → ζ = 1.5
        let s = SpringJoint::new(0, 1, [0.0; 3], [1.0, 0.0, 0.0], 1.0, 100.0, 30.0, 1.0);
        assert!(s.is_overdamped());
    }

    #[test]
    fn test_spring_tension_only() {
        let mut s = SpringJoint::new(0, 1, [0.0; 3], [2.0, 0.0, 0.0], 2.0, 100.0, 0.0, 1.0);
        s.set_tension_only(true);
        // dist < rest_length → compression → no force
        let f = s.compute_force(1.5, 0.0);
        assert_eq!(f, 0.0);
    }

    // --- RopeConstraint ---

    #[test]
    fn test_rope_num_segments() {
        let rope = RopeConstraint::new(
            [0.0; 3],
            [10.0, 0.0, 0.0],
            5,
            0.2,
            1.0,
            0.1,
            [0.0, -9.81, 0.0],
        );
        assert_eq!(rope.num_segments(), 5);
    }

    #[test]
    fn test_rope_initial_length() {
        let rope = RopeConstraint::new(
            [0.0; 3],
            [10.0, 0.0, 0.0],
            4,
            0.2,
            1.0,
            0.1,
            [0.0, -9.81, 0.0],
        );
        let len = rope.current_length();
        assert!((len - 10.0).abs() < 1e-8, "len={len}");
    }

    #[test]
    fn test_rope_step_gravity() {
        let mut rope = RopeConstraint::new(
            [0.0; 3],
            [0.0, 0.0, 0.0],
            3,
            1.0,
            1.0,
            0.1,
            [0.0, -9.81, 0.0],
        );
        // Unpin middle node and step
        rope.nodes[1].pinned = false;
        rope.nodes[2].pinned = false;
        let y_before = rope.nodes[1].position[1];
        rope.step(0.016);
        let y_after = rope.nodes[1].position[1];
        assert!(y_after < y_before + 1e-3, "gravity should pull node down");
    }

    // --- WheelConstraint ---

    #[test]
    fn test_wheel_contact_patch_speed() {
        let mut w = WheelConstraint::new(0, 0.3, 0.2);
        w.spin_velocity = 10.0;
        let speed = w.contact_patch_speed();
        assert!((speed - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_wheel_suspension_force_airborne() {
        let w = WheelConstraint::new(0, 0.3, 0.2);
        assert_eq!(w.suspension_force(0.0), 0.0);
    }

    #[test]
    fn test_wheel_suspension_force_grounded() {
        let mut w = WheelConstraint::new(0, 0.3, 0.2);
        w.contact = WheelContact::Grounded;
        w.suspension_compression = 0.5;
        let f = w.suspension_force(0.0);
        // 30000 * 0.5 * 0.2 = 3000
        assert!((f - 3000.0).abs() < 1e-6, "f={f}");
    }

    #[test]
    fn test_wheel_integrate_spin() {
        let mut w = WheelConstraint::new(0, 0.3, 0.2);
        w.drive_torque = 100.0;
        w.integrate_spin(2.0, 0.1); // alpha = 50 rad/s²
        assert!(w.spin_velocity > 0.0);
    }

    // --- CharacterConstraint ---

    #[test]
    fn test_character_slope_walkable_flat() {
        let c = CharacterConstraint::new(0, 1.8, 0.4);
        assert!(c.is_slope_walkable());
    }

    #[test]
    fn test_character_slope_steep() {
        let mut c = CharacterConstraint::new(0, 1.8, 0.4);
        // Steep slope — normal pointing almost horizontally
        c.ground_normal = v3_norm([0.1, 0.1, 0.0]);
        assert!(!c.is_slope_walkable());
    }

    #[test]
    fn test_character_step_up() {
        let c = CharacterConstraint::new(0, 1.8, 0.4);
        assert!(c.can_step_up(0.3));
        assert!(!c.can_step_up(0.5));
    }

    #[test]
    fn test_character_jump() {
        let mut c = CharacterConstraint::new(0, 1.8, 0.4);
        c.grounded = true;
        c.jump_requested = true;
        let impulse = c.try_jump();
        assert!((impulse - c.jump_velocity).abs() < 1e-12);
        assert!(c.velocity[1] > 0.0);
    }

    #[test]
    fn test_character_no_jump_airborne() {
        let mut c = CharacterConstraint::new(0, 1.8, 0.4);
        c.grounded = false;
        c.jump_requested = true;
        let impulse = c.try_jump();
        assert_eq!(impulse, 0.0);
    }

    // --- TriggerVolume ---

    #[test]
    fn test_trigger_aabb_enter() {
        let shape = TriggerShape::Aabb {
            centre: [0.0; 3],
            half_extents: [1.0; 3],
        };
        let mut vol = TriggerVolume::new(1, shape, 0xFFFF_FFFF);
        let events = vol.update(&[(0, [0.5, 0.5, 0.5], 0x0001)]);
        assert_eq!(events.len(), 1);
        assert_eq!(events[0], (0, TriggerEvent::Enter));
    }

    #[test]
    fn test_trigger_sphere_enter_exit() {
        let shape = TriggerShape::Sphere {
            centre: [0.0; 3],
            radius: 1.0,
        };
        let mut vol = TriggerVolume::new(2, shape, 0xFF);
        // Step 1: enter
        vol.update(&[(5, [0.0, 0.0, 0.5], 0x01)]);
        // Step 2: exit
        let events = vol.update(&[(5, [2.0, 0.0, 0.0], 0x01)]);
        let has_exit = events.iter().any(|(_, e)| *e == TriggerEvent::Exit);
        assert!(has_exit);
    }

    #[test]
    fn test_trigger_filter_mask() {
        let shape = TriggerShape::Aabb {
            centre: [0.0; 3],
            half_extents: [2.0; 3],
        };
        let mut vol = TriggerVolume::new(3, shape, 0x01); // only layer bit 0
        // Body with layer 0x02 should be ignored
        let events = vol.update(&[(0, [0.0; 3], 0x02)]);
        assert_eq!(events.len(), 0);
    }

    // --- CompositeConstraint ---

    #[test]
    fn test_composite_sorting() {
        let mut comp = CompositeConstraint::new(4);
        comp.add(SubConstraintRef::new("low", ConstraintPriority::Low, 0.0));
        comp.add(SubConstraintRef::new(
            "critical",
            ConstraintPriority::Critical,
            0.0,
        ));
        comp.add(SubConstraintRef::new("high", ConstraintPriority::High, 0.0));
        let sorted = comp.sorted_by_priority();
        assert_eq!(sorted[0].name, "critical");
        assert_eq!(sorted[1].name, "high");
        assert_eq!(sorted[2].name, "low");
    }

    #[test]
    fn test_composite_activate_deactivate() {
        let mut comp = CompositeConstraint::new(2);
        comp.add(SubConstraintRef::new("a", ConstraintPriority::High, 0.0));
        comp.add(SubConstraintRef::new("b", ConstraintPriority::High, 0.0));
        comp.deactivate_all();
        assert_eq!(comp.active_count(), 0);
        comp.activate_all();
        assert_eq!(comp.active_count(), 2);
    }

    // --- ConstraintSolver ---

    #[test]
    fn test_solver_single_contact_constraint() {
        // Two bodies approaching each other; contact normal = (1,0,0)
        let mut bodies = vec![
            SolverBody::new([-1.0, 0.0, 0.0], [0.0; 3], 1.0, [1.0; 3]),
            SolverBody::new([1.0, 0.0, 0.0], [0.0; 3], 1.0, [1.0; 3]),
        ];
        let n = [1.0f64, 0.0, 0.0];
        let row = ConstraintRow::new(
            v3_scale(n, -1.0), // j_lin_a = -n
            [0.0; 3],
            n, // j_lin_b = +n
            [0.0; 3],
            0.0,
            0.5, // eff_mass = 0.5
            0.0,
            f64::INFINITY,
        );
        let mut pair = ConstraintPair::new(0, 1);
        pair.add_row(row);
        let mut pairs = vec![pair];
        let solver = ConstraintSolver::new(10);
        solver.solve(&mut pairs, &mut bodies);
        // After solving, relative normal velocity should be >= 0 (no longer approaching)
        let rel_vn = bodies[1].vel[0] - bodies[0].vel[0];
        assert!(rel_vn >= -1e-6, "rel_vn={rel_vn}");
    }

    #[test]
    fn test_solver_warm_start_applies_impulse() {
        let mut bodies = vec![
            SolverBody::new([0.0; 3], [0.0; 3], 1.0, [1.0; 3]),
            SolverBody::new([0.0; 3], [0.0; 3], 1.0, [1.0; 3]),
        ];
        let mut row = ConstraintRow::new(
            [-1.0, 0.0, 0.0],
            [0.0; 3],
            [1.0, 0.0, 0.0],
            [0.0; 3],
            0.0,
            0.5,
            f64::NEG_INFINITY,
            f64::INFINITY,
        );
        row.accumulated = 2.0;
        let mut pair = ConstraintPair::new(0, 1);
        pair.add_row(row);
        let pairs = vec![pair];
        let solver = ConstraintSolver::new(4);
        solver.warm_start(&pairs, &mut bodies);
        // Body A should have gained velocity in -x direction
        assert!(bodies[0].vel[0] < 0.0, "vel_a_x={}", bodies[0].vel[0]);
    }

    #[test]
    fn test_solver_sleep_transition() {
        let mut solver = ConstraintSolver::new(4);
        let mut bodies = vec![SolverBody::new([0.0; 3], [0.0; 3], 1.0, [1.0; 3])];
        // Kinetic energy is zero → should accumulate sleep counter
        for _ in 0..65 {
            solver.update_sleep(&mut bodies);
        }
        assert!(bodies[0].sleeping, "body should be sleeping");
    }

    #[test]
    fn test_solver_wake_all() {
        let mut solver = ConstraintSolver::new(4);
        let mut bodies = vec![SolverBody::new([0.0; 3], [0.0; 3], 1.0, [1.0; 3])];
        for _ in 0..65 {
            solver.update_sleep(&mut bodies);
        }
        assert!(bodies[0].sleeping);
        solver.wake_all(&mut bodies);
        assert!(!bodies[0].sleeping);
    }
}
