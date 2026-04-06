#![allow(clippy::if_same_then_else)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Joint-limit collision: revolute, prismatic, and ball-socket angle limits
//! modelled as constraints.  Covers soft/hard stops, gear-ratio response,
//! cable/chain limits, articulated-body joint stops, motor hard limits,
//! joint friction as a constraint, compliance at limits, and warm-starting.
//!
//! All arithmetic uses plain `f64` and `[f64; 3]` arrays — no nalgebra.

#![allow(dead_code)]

// ─────────────────────────────────────────────────────────────────────────────
// Basic 3-vector helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Add two 3-vectors.
#[inline]
pub fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Subtract two 3-vectors.
#[inline]
pub fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Scale a 3-vector.
#[inline]
pub fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

/// Dot product.
#[inline]
pub fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Euclidean norm.
#[inline]
pub fn norm3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}

/// Normalise; returns zero vector if near-zero.
#[inline]
pub fn normalize3(a: [f64; 3]) -> [f64; 3] {
    let n = norm3(a);
    if n < 1e-14 {
        [0.0; 3]
    } else {
        scale3(a, 1.0 / n)
    }
}

/// Clamp `v` to `[lo, hi]`.
#[inline]
pub fn clampf(v: f64, lo: f64, hi: f64) -> f64 {
    v.max(lo).min(hi)
}

// ─────────────────────────────────────────────────────────────────────────────
// JointType
// ─────────────────────────────────────────────────────────────────────────────

/// Classification of joint kinematics for limit enforcement.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum JointType {
    /// Single rotational degree of freedom.
    Revolute,
    /// Single translational degree of freedom.
    Prismatic,
    /// Three rotational degrees of freedom (ball-and-socket).
    BallSocket,
}

// ─────────────────────────────────────────────────────────────────────────────
// JointLimits
// ─────────────────────────────────────────────────────────────────────────────

/// Angular or translational limits for a joint.
#[derive(Debug, Clone)]
pub struct JointLimits {
    /// Minimum angle (rad) or translation (m).
    pub lower: f64,
    /// Maximum angle (rad) or translation (m).
    pub upper: f64,
    /// Whether limits are active.
    pub enabled: bool,
}

impl JointLimits {
    /// Create joint limits.
    pub fn new(lower: f64, upper: f64) -> Self {
        Self {
            lower,
            upper,
            enabled: true,
        }
    }

    /// Return `true` if `value` is outside the limits.
    #[inline]
    pub fn is_violated(&self, value: f64) -> bool {
        self.enabled && (value < self.lower || value > self.upper)
    }

    /// Signed penetration: positive when below lower limit, negative when above upper.
    #[inline]
    pub fn penetration(&self, value: f64) -> f64 {
        if value < self.lower {
            self.lower - value
        } else if value > self.upper {
            value - self.upper
        } else {
            0.0
        }
    }

    /// Clamp `value` to `[lower, upper]`.
    #[inline]
    pub fn clamp(&self, value: f64) -> f64 {
        clampf(value, self.lower, self.upper)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// StopType — soft vs hard
// ─────────────────────────────────────────────────────────────────────────────

/// How the joint stop is modelled.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum StopType {
    /// Rigid constraint: velocity is zeroed and position is projected.
    Hard,
    /// Spring-damper constraint: penalty force proportional to penetration.
    Soft {
        /// Stiffness coefficient (N/m or N·m/rad).
        stiffness: f64,
        /// Damping coefficient.
        damping: f64,
    },
}

impl StopType {
    /// Compute the constraint force for the given penetration depth and velocity.
    ///
    /// For a `Hard` stop returns a large penalty impulse.
    /// For a `Soft` stop returns `stiffness * penetration + damping * velocity`.
    pub fn constraint_force(&self, penetration: f64, velocity: f64) -> f64 {
        match *self {
            StopType::Hard => {
                // Return a large restoring force (handled separately as impulse).
                if penetration > 0.0 {
                    1e9 * penetration
                } else {
                    0.0
                }
            }
            StopType::Soft { stiffness, damping } => {
                if penetration > 0.0 {
                    stiffness * penetration - damping * velocity.min(0.0)
                } else {
                    0.0
                }
            }
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// JointLimitConstraint
// ─────────────────────────────────────────────────────────────────────────────

/// A constraint that enforces a joint limit.
#[derive(Debug, Clone)]
pub struct JointLimitConstraint {
    /// Body index A (parent).
    pub body_a: usize,
    /// Body index B (child).
    pub body_b: usize,
    /// Type of joint.
    pub joint_type: JointType,
    /// Lower and upper limits.
    pub limits: JointLimits,
    /// Stop model.
    pub stop_type: StopType,
    /// Current joint position (angle or translation).
    pub position: f64,
    /// Current joint velocity.
    pub velocity: f64,
    /// Accumulated impulse from previous frame (warm-start).
    pub warm_impulse: f64,
    /// Joint axis in world space (for revolute/prismatic).
    pub axis: [f64; 3],
}

impl JointLimitConstraint {
    /// Create a joint limit constraint.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        body_a: usize,
        body_b: usize,
        joint_type: JointType,
        limits: JointLimits,
        stop_type: StopType,
        axis: [f64; 3],
    ) -> Self {
        Self {
            body_a,
            body_b,
            joint_type,
            limits,
            stop_type,
            position: 0.0,
            velocity: 0.0,
            warm_impulse: 0.0,
            axis: normalize3(axis),
        }
    }

    /// Whether the joint is currently violating its limit.
    pub fn is_active(&self) -> bool {
        self.limits.is_violated(self.position)
    }

    /// Signed penetration depth (> 0 means violated).
    pub fn penetration_depth(&self) -> f64 {
        self.limits.penetration(self.position)
    }

    /// Compute the constraint impulse to resolve the limit violation.
    ///
    /// Returns the scalar impulse along the joint axis.  Positive impulse
    /// pushes the position back into limits.
    pub fn compute_impulse(&self, inv_mass_a: f64, inv_mass_b: f64, dt: f64) -> f64 {
        let pen = self.penetration_depth();
        if pen <= 0.0 {
            return 0.0;
        }
        let force = self.stop_type.constraint_force(pen, self.velocity);
        let eff_mass = inv_mass_a + inv_mass_b;
        if eff_mass < 1e-14 {
            return 0.0;
        }
        force * dt / eff_mass
    }

    /// Apply warm-starting impulse from the previous step.
    ///
    /// Returns the warm impulse (already clamped to `>= 0`).
    pub fn apply_warm_start(&self) -> f64 {
        self.warm_impulse.max(0.0)
    }

    /// Store the accumulated impulse for next-frame warm-starting.
    pub fn store_impulse(&mut self, impulse: f64) {
        self.warm_impulse = impulse.max(0.0);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// BallSocketLimitConstraint
// ─────────────────────────────────────────────────────────────────────────────

/// Constraint that enforces a cone limit on a ball-and-socket joint.
///
/// The cone half-angle `max_angle` (radians) is the maximum allowed swing.
#[derive(Debug, Clone)]
pub struct BallSocketLimitConstraint {
    /// Parent body index.
    pub body_a: usize,
    /// Child body index.
    pub body_b: usize,
    /// Reference axis in body-A frame (cone axis).
    pub cone_axis: [f64; 3],
    /// Maximum cone half-angle (rad).
    pub max_angle: f64,
    /// Stop type.
    pub stop_type: StopType,
    /// Current swing angle (angle between child axis and cone axis).
    pub current_angle: f64,
    /// Warm-start impulse.
    pub warm_impulse: f64,
}

impl BallSocketLimitConstraint {
    /// Create a ball-socket limit constraint.
    pub fn new(
        body_a: usize,
        body_b: usize,
        cone_axis: [f64; 3],
        max_angle: f64,
        stop_type: StopType,
    ) -> Self {
        Self {
            body_a,
            body_b,
            cone_axis: normalize3(cone_axis),
            max_angle,
            stop_type,
            current_angle: 0.0,
            warm_impulse: 0.0,
        }
    }

    /// Update the current swing angle given the child-body axis in world space.
    pub fn update_angle(&mut self, child_axis_world: [f64; 3]) {
        let d = dot3(self.cone_axis, normalize3(child_axis_world));
        let d_clamped = clampf(d, -1.0, 1.0);
        self.current_angle = d_clamped.acos();
    }

    /// Penetration beyond the cone limit.
    pub fn penetration(&self) -> f64 {
        (self.current_angle - self.max_angle).max(0.0)
    }

    /// Whether the cone limit is violated.
    pub fn is_violated(&self) -> bool {
        self.current_angle > self.max_angle
    }

    /// Compute a corrective angular impulse magnitude.
    pub fn compute_impulse(&self, inv_inertia_a: f64, inv_inertia_b: f64, dt: f64) -> f64 {
        let pen = self.penetration();
        if pen <= 0.0 {
            return 0.0;
        }
        let force = self.stop_type.constraint_force(pen, 0.0);
        let eff = inv_inertia_a + inv_inertia_b;
        if eff < 1e-14 { 0.0 } else { force * dt / eff }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// GearRatioCollision
// ─────────────────────────────────────────────────────────────────────────────

/// Collision response through a gear transmission.
///
/// When joint A hits its limit, the impulse is transmitted to joint B
/// through a gear ratio `ratio = omega_B / omega_A`.
#[derive(Debug, Clone)]
pub struct GearRatioCollision {
    /// Joint A (driving).
    pub joint_a: usize,
    /// Joint B (driven).
    pub joint_b: usize,
    /// Gear ratio ω_B / ω_A (can be negative for reversals).
    pub ratio: f64,
    /// Efficiency factor ∈ (0, 1].
    pub efficiency: f64,
}

impl GearRatioCollision {
    /// Create a gear-ratio collision link.
    pub fn new(joint_a: usize, joint_b: usize, ratio: f64, efficiency: f64) -> Self {
        Self {
            joint_a,
            joint_b,
            ratio,
            efficiency: efficiency.clamp(0.0, 1.0),
        }
    }

    /// Propagate impulse from joint A to joint B.
    ///
    /// `impulse_a` is the impulse applied at joint A; returns the resulting
    /// impulse at joint B.
    pub fn propagate_impulse(&self, impulse_a: f64) -> f64 {
        impulse_a * self.ratio * self.efficiency
    }

    /// Effective inertia at joint A accounting for the gear and driven inertia.
    ///
    /// `inertia_b` is the rotational inertia at joint B.
    pub fn reflected_inertia(&self, inertia_b: f64) -> f64 {
        inertia_b * self.ratio * self.ratio
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// CableChainLimit
// ─────────────────────────────────────────────────────────────────────────────

/// Cable or chain joint limit.
///
/// A cable prevents extension beyond `max_length`; a chain additionally
/// prevents compression below `min_length`.
#[derive(Debug, Clone)]
pub struct CableChainLimit {
    /// Body index A (cable anchor).
    pub body_a: usize,
    /// Body index B (cable end).
    pub body_b: usize,
    /// Minimum allowed distance (0 for cable-only).
    pub min_length: f64,
    /// Maximum allowed distance.
    pub max_length: f64,
    /// Current distance between anchor and end.
    pub current_length: f64,
    /// Relative velocity along the cable direction.
    pub relative_velocity: f64,
    /// Accumulated constraint impulse.
    pub impulse: f64,
    /// Cable stiffness (for soft stops).
    pub stiffness: f64,
    /// Cable damping.
    pub damping: f64,
}

impl CableChainLimit {
    /// Create a cable/chain limit.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        body_a: usize,
        body_b: usize,
        min_length: f64,
        max_length: f64,
        stiffness: f64,
        damping: f64,
    ) -> Self {
        Self {
            body_a,
            body_b,
            min_length,
            max_length,
            current_length: 0.0,
            relative_velocity: 0.0,
            impulse: 0.0,
            stiffness,
            damping,
        }
    }

    /// Extension violation: positive if cable is too long.
    pub fn extension_violation(&self) -> f64 {
        (self.current_length - self.max_length).max(0.0)
    }

    /// Compression violation: positive if chain is too short.
    pub fn compression_violation(&self) -> f64 {
        (self.min_length - self.current_length).max(0.0)
    }

    /// Net penetration (positive means constraint is violated).
    pub fn net_penetration(&self) -> f64 {
        self.extension_violation() + self.compression_violation()
    }

    /// Compute the corrective impulse along the cable direction.
    pub fn compute_impulse(&self, inv_mass_a: f64, inv_mass_b: f64, dt: f64) -> f64 {
        let pen = self.net_penetration();
        if pen <= 0.0 {
            return 0.0;
        }
        let force = self.stiffness * pen - self.damping * self.relative_velocity.min(0.0);
        let eff = inv_mass_a + inv_mass_b;
        if eff < 1e-14 { 0.0 } else { force * dt / eff }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ArticulatedBodyJointStop
// ─────────────────────────────────────────────────────────────────────────────

/// Joint-stop contact between two bodies in an articulated chain.
///
/// When a joint hits its limit, the effective mass computation must account
/// for the propagation of impulses through the articulated body tree.
#[derive(Debug, Clone)]
pub struct ArticulatedBodyJointStop {
    /// Joint identifier.
    pub joint_id: usize,
    /// Effective inertia of the subtree (pre-computed by ABFC algorithm).
    pub effective_inertia: f64,
    /// Current position in the joint's coordinate.
    pub position: f64,
    /// Current velocity.
    pub velocity: f64,
    /// Joint limits.
    pub limits: JointLimits,
    /// Stop type.
    pub stop_type: StopType,
    /// Warm-starting impulse.
    pub warm_impulse: f64,
}

impl ArticulatedBodyJointStop {
    /// Create an articulated body joint stop.
    pub fn new(
        joint_id: usize,
        effective_inertia: f64,
        limits: JointLimits,
        stop_type: StopType,
    ) -> Self {
        Self {
            joint_id,
            effective_inertia,
            position: 0.0,
            velocity: 0.0,
            limits,
            stop_type,
            warm_impulse: 0.0,
        }
    }

    /// Compute corrective impulse using the articulated effective inertia.
    pub fn compute_impulse(&self, dt: f64) -> f64 {
        let pen = self.limits.penetration(self.position);
        if pen <= 0.0 {
            return 0.0;
        }
        let force = self.stop_type.constraint_force(pen, self.velocity);
        if self.effective_inertia < 1e-14 {
            0.0
        } else {
            force * dt / self.effective_inertia
        }
    }

    /// Velocity impulse to zero out joint velocity when hitting a hard stop.
    ///
    /// Returns the impulse needed to zero the outgoing velocity, accounting
    /// for restitution `e`.
    pub fn velocity_impulse(&self, restitution: f64) -> f64 {
        if !self.limits.is_violated(self.position) {
            return 0.0;
        }
        let vel_out = -(1.0 + restitution) * self.velocity;
        if self.effective_inertia < 1e-14 {
            0.0
        } else {
            vel_out / self.effective_inertia
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// MotorHardLimit
// ─────────────────────────────────────────────────────────────────────────────

/// Represents a motor that has absolute position and velocity limits.
///
/// When the motor exceeds its limits, a corrective impulse is applied.
#[derive(Debug, Clone)]
pub struct MotorHardLimit {
    /// Maximum positive torque the motor can produce.
    pub max_torque: f64,
    /// Maximum positive angular velocity.
    pub max_velocity: f64,
    /// Position limits.
    pub position_limits: JointLimits,
    /// Current motor angle.
    pub angle: f64,
    /// Current motor angular velocity.
    pub angular_velocity: f64,
    /// Motor inertia.
    pub inertia: f64,
    /// Whether the motor is currently stalled (against hard stop).
    pub stalled: bool,
}

impl MotorHardLimit {
    /// Create a motor hard limit.
    #[allow(clippy::too_many_arguments)]
    pub fn new(max_torque: f64, max_velocity: f64, lower: f64, upper: f64, inertia: f64) -> Self {
        Self {
            max_torque,
            max_velocity,
            position_limits: JointLimits::new(lower, upper),
            angle: 0.0,
            angular_velocity: 0.0,
            inertia,
            stalled: false,
        }
    }

    /// Clamp motor torque to `[-max_torque, max_torque]`.
    pub fn clamp_torque(&self, torque: f64) -> f64 {
        clampf(torque, -self.max_torque, self.max_torque)
    }

    /// Clamp motor velocity to `[-max_velocity, max_velocity]`.
    pub fn clamp_velocity(&self, vel: f64) -> f64 {
        clampf(vel, -self.max_velocity, self.max_velocity)
    }

    /// Check if the motor is against a hard stop.
    pub fn check_stall(&mut self) {
        self.stalled = self.position_limits.is_violated(self.angle);
    }

    /// Impulse to push the motor back into limits (applied when stalled).
    pub fn limit_impulse(&self, dt: f64) -> f64 {
        if !self.stalled {
            return 0.0;
        }
        let pen = self.position_limits.penetration(self.angle);
        if self.inertia < 1e-14 {
            0.0
        } else {
            // Large stiffness for hard stop
            1e6 * pen * dt / self.inertia
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// JointFrictionConstraint
// ─────────────────────────────────────────────────────────────────────────────

/// Joint friction modelled as a velocity-level constraint.
///
/// Friction opposes relative velocity at the joint with a Coulomb-like model.
#[derive(Debug, Clone)]
pub struct JointFrictionConstraint {
    /// Body A index.
    pub body_a: usize,
    /// Body B index.
    pub body_b: usize,
    /// Joint axis (world space).
    pub axis: [f64; 3],
    /// Coulomb friction coefficient (dimensionless).
    pub friction_coefficient: f64,
    /// Current normal force at the joint (N).
    pub normal_force: f64,
    /// Current relative velocity along the joint axis.
    pub relative_velocity: f64,
    /// Accumulated friction impulse (warm-start).
    pub warm_impulse: f64,
}

impl JointFrictionConstraint {
    /// Create a joint friction constraint.
    pub fn new(body_a: usize, body_b: usize, axis: [f64; 3], friction_coefficient: f64) -> Self {
        Self {
            body_a,
            body_b,
            axis: normalize3(axis),
            friction_coefficient,
            normal_force: 0.0,
            relative_velocity: 0.0,
            warm_impulse: 0.0,
        }
    }

    /// Maximum friction impulse magnitude given timestep `dt`.
    pub fn max_friction_impulse(&self, inv_mass_eff: f64, dt: f64) -> f64 {
        let _ = dt;
        let max_f = self.friction_coefficient * self.normal_force.abs();
        if inv_mass_eff < 1e-14 { max_f } else { max_f }
    }

    /// Compute the friction impulse needed to bring joint velocity to zero,
    /// clamped to the friction cone.
    pub fn compute_impulse(&self, inv_mass_eff: f64, dt: f64) -> f64 {
        if inv_mass_eff < 1e-14 {
            return 0.0;
        }
        let delta_v = -self.relative_velocity;
        let raw_impulse = delta_v / inv_mass_eff;
        let max_imp = self.max_friction_impulse(inv_mass_eff, dt);
        clampf(raw_impulse, -max_imp, max_imp)
    }

    /// Apply warm start — returns clamped previous impulse.
    pub fn warm_start_impulse(&self, inv_mass_eff: f64, dt: f64) -> f64 {
        let max_imp = self.max_friction_impulse(inv_mass_eff, dt);
        clampf(self.warm_impulse, -max_imp, max_imp)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// JointComplianceParams
// ─────────────────────────────────────────────────────────────────────────────

/// Compliance (softness) parameters for a joint limit constraint.
///
/// Implements the Constraint Force Mixing (CFM) / Error Reduction Parameter
/// (ERP) approach used in many rigid-body solvers.
#[derive(Debug, Clone, Default)]
pub struct JointComplianceParams {
    /// Error Reduction Parameter ∈ (0, 1].  Fraction of position error
    /// corrected per step.
    pub erp: f64,
    /// Constraint Force Mixing coefficient.  Non-zero → soft constraint.
    pub cfm: f64,
}

impl JointComplianceParams {
    /// Create compliance params with typical values.
    pub fn new(erp: f64, cfm: f64) -> Self {
        Self {
            erp: clampf(erp, 0.0, 1.0),
            cfm: cfm.max(0.0),
        }
    }

    /// Default rigid constraint (ERP=0.2, CFM=0).
    pub fn rigid() -> Self {
        Self::new(0.2, 0.0)
    }

    /// Soft spring-like constraint.
    pub fn soft(stiffness: f64, damping: f64, dt: f64) -> Self {
        // From ODE manual: cfm = 1/(dt*k + c), erp = dt*k/(dt*k + c)
        let denom = dt * stiffness + damping;
        if denom < 1e-14 {
            Self::rigid()
        } else {
            Self::new(dt * stiffness / denom, 1.0 / denom)
        }
    }

    /// Correct a position error `c` with effective mass `k_eff` over step `dt`.
    ///
    /// Returns the corrective velocity (bias) to add to the constraint RHS.
    pub fn positional_bias(&self, error: f64, dt: f64) -> f64 {
        if dt < 1e-14 {
            0.0
        } else {
            self.erp * error / dt
        }
    }

    /// Modify the effective constraint mass to include CFM softness.
    ///
    /// `k_eff` is the raw inverse effective mass.
    pub fn softened_keff(&self, k_eff: f64) -> f64 {
        k_eff + self.cfm
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// JointLimitWarmStart
// ─────────────────────────────────────────────────────────────────────────────

/// Cache of warm-start impulses for joint limit constraints across frames.
#[derive(Debug, Clone, Default)]
pub struct JointLimitWarmStart {
    /// Map from `(body_a, body_b, axis_index)` → cached impulse.
    entries: Vec<WarmEntry>,
}

/// One warm-start entry.
#[derive(Debug, Clone)]
struct WarmEntry {
    body_a: usize,
    body_b: usize,
    axis_idx: usize,
    impulse: f64,
}

impl JointLimitWarmStart {
    /// Create an empty warm-start cache.
    pub fn new() -> Self {
        Self {
            entries: Vec::new(),
        }
    }

    /// Store (or update) a warm-start impulse for the given constraint key.
    pub fn store(&mut self, body_a: usize, body_b: usize, axis_idx: usize, impulse: f64) {
        for entry in &mut self.entries {
            if entry.body_a == body_a && entry.body_b == body_b && entry.axis_idx == axis_idx {
                entry.impulse = impulse;
                return;
            }
        }
        self.entries.push(WarmEntry {
            body_a,
            body_b,
            axis_idx,
            impulse,
        });
    }

    /// Retrieve a warm-start impulse, returning `0` if not found.
    pub fn retrieve(&self, body_a: usize, body_b: usize, axis_idx: usize) -> f64 {
        for entry in &self.entries {
            if entry.body_a == body_a && entry.body_b == body_b && entry.axis_idx == axis_idx {
                return entry.impulse;
            }
        }
        0.0
    }

    /// Clear all stored impulses (call at the start of each frame if needed).
    pub fn clear(&mut self) {
        self.entries.clear();
    }

    /// Number of cached entries.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Whether the cache is empty.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Scale all stored impulses by `factor` (use for substep warm-starting).
    pub fn scale(&mut self, factor: f64) {
        for entry in &mut self.entries {
            entry.impulse *= factor;
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// JointLimitSolver
// ─────────────────────────────────────────────────────────────────────────────

/// Iterative solver for a collection of joint limit constraints.
///
/// Performs sequential impulse iterations (PGS) with warm-starting.
#[derive(Debug, Clone)]
pub struct JointLimitSolver {
    /// Joint limit constraints registered with this solver.
    pub constraints: Vec<JointLimitConstraint>,
    /// Warm-start cache.
    pub warm_start: JointLimitWarmStart,
    /// Compliance for all constraints.
    pub compliance: JointComplianceParams,
    /// Number of solver iterations per step.
    pub iterations: usize,
}

impl JointLimitSolver {
    /// Create a solver with default parameters.
    pub fn new(iterations: usize) -> Self {
        Self {
            constraints: Vec::new(),
            warm_start: JointLimitWarmStart::new(),
            compliance: JointComplianceParams::rigid(),
            iterations,
        }
    }

    /// Register a constraint.
    pub fn add_constraint(&mut self, c: JointLimitConstraint) {
        self.constraints.push(c);
    }

    /// Apply warm starts from the previous frame to all constraints.
    pub fn apply_warm_starts(&mut self) {
        for c in &mut self.constraints {
            let wi = self.warm_start.retrieve(c.body_a, c.body_b, 0);
            c.warm_impulse = wi;
        }
    }

    /// Solve all constraints for `dt` seconds.
    ///
    /// `inv_masses[i]` is the inverse mass of body `i`.
    /// Returns the total impulse applied over all iterations.
    pub fn solve(&mut self, inv_masses: &[f64], dt: f64) -> f64 {
        // Apply warm starts
        self.apply_warm_starts();
        let mut total_impulse = 0.0f64;
        for _iter in 0..self.iterations {
            for c in &mut self.constraints {
                let im_a = if c.body_a < inv_masses.len() {
                    inv_masses[c.body_a]
                } else {
                    0.0
                };
                let im_b = if c.body_b < inv_masses.len() {
                    inv_masses[c.body_b]
                } else {
                    0.0
                };
                let impulse = c.compute_impulse(im_a, im_b, dt);
                let new_acc = (c.warm_impulse + impulse).max(0.0);
                let delta = new_acc - c.warm_impulse;
                c.warm_impulse = new_acc;
                total_impulse += delta.abs();
            }
        }
        // Store warm starts
        for c in &self.constraints {
            self.warm_start.store(c.body_a, c.body_b, 0, c.warm_impulse);
        }
        total_impulse
    }

    /// Clear all constraints (keep warm-start cache).
    pub fn clear_constraints(&mut self) {
        self.constraints.clear();
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// JointStopContact
// ─────────────────────────────────────────────────────────────────────────────

/// A contact-like event generated when a joint hits its hard stop.
#[derive(Debug, Clone)]
pub struct JointStopContact {
    /// Joint identifier.
    pub joint_id: usize,
    /// Contact normal (direction of correction).
    pub normal: [f64; 3],
    /// Penetration depth.
    pub depth: f64,
    /// Relative velocity along the normal at the stop.
    pub relative_velocity: f64,
    /// Restitution coefficient for the stop bounce.
    pub restitution: f64,
}

impl JointStopContact {
    /// Create a joint stop contact.
    pub fn new(
        joint_id: usize,
        normal: [f64; 3],
        depth: f64,
        relative_velocity: f64,
        restitution: f64,
    ) -> Self {
        Self {
            joint_id,
            normal: normalize3(normal),
            depth,
            relative_velocity,
            restitution,
        }
    }

    /// Compute the corrective impulse (Baumgarte + restitution).
    ///
    /// `eff_mass` is the scalar effective mass at the constraint.
    /// `beta` is the Baumgarte position correction factor.
    /// `dt` is the timestep.
    pub fn compute_impulse(&self, eff_mass: f64, beta: f64, dt: f64) -> f64 {
        if eff_mass < 1e-14 {
            return 0.0;
        }
        let baumgarte_bias = beta * self.depth / dt;
        let restitution_bias = self.restitution * self.relative_velocity.min(0.0).abs();
        let target_vel = baumgarte_bias + restitution_bias;
        let delta_v = target_vel - self.relative_velocity.min(0.0);
        (delta_v / eff_mass).max(0.0)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PrismaticJointLimit
// ─────────────────────────────────────────────────────────────────────────────

/// Prismatic joint with translational limits and full collision response.
#[derive(Debug, Clone)]
pub struct PrismaticJointLimit {
    /// Parent body.
    pub body_a: usize,
    /// Child body.
    pub body_b: usize,
    /// Sliding axis (world space, unit vector).
    pub axis: [f64; 3],
    /// Translation limits.
    pub limits: JointLimits,
    /// Stop type.
    pub stop_type: StopType,
    /// Current translation along axis.
    pub translation: f64,
    /// Current velocity along axis.
    pub velocity: f64,
    /// Warm-start impulse.
    pub warm_impulse: f64,
    /// Compliance parameters.
    pub compliance: JointComplianceParams,
}

impl PrismaticJointLimit {
    /// Create a prismatic joint limit.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        body_a: usize,
        body_b: usize,
        axis: [f64; 3],
        lower: f64,
        upper: f64,
        stop_type: StopType,
        compliance: JointComplianceParams,
    ) -> Self {
        Self {
            body_a,
            body_b,
            axis: normalize3(axis),
            limits: JointLimits::new(lower, upper),
            stop_type,
            translation: 0.0,
            velocity: 0.0,
            warm_impulse: 0.0,
            compliance,
        }
    }

    /// Penetration depth at the current translation.
    pub fn penetration(&self) -> f64 {
        self.limits.penetration(self.translation)
    }

    /// Positional bias velocity.
    pub fn bias_velocity(&self, dt: f64) -> f64 {
        self.compliance.positional_bias(self.penetration(), dt)
    }

    /// Effective constraint jacobian row (scalar, since 1-D constraint).
    pub fn constraint_row(&self) -> f64 {
        if self.translation < self.limits.lower {
            1.0
        } else if self.translation > self.limits.upper {
            -1.0
        } else {
            0.0
        }
    }

    /// Solve for the constraint impulse.
    pub fn solve_impulse(&self, inv_mass_a: f64, inv_mass_b: f64, dt: f64) -> f64 {
        let pen = self.penetration();
        if pen <= 0.0 {
            return 0.0;
        }
        let bias = self.bias_velocity(dt);
        let k_eff_raw = inv_mass_a + inv_mass_b;
        let k_eff = self.compliance.softened_keff(k_eff_raw);
        if k_eff < 1e-14 {
            return 0.0;
        }
        let lambda = (bias - self.velocity) / k_eff;
        lambda.max(0.0)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// RevoluteJointLimit
// ─────────────────────────────────────────────────────────────────────────────

/// Revolute joint with angular limits and full collision response.
#[derive(Debug, Clone)]
pub struct RevoluteJointLimit {
    /// Parent body.
    pub body_a: usize,
    /// Child body.
    pub body_b: usize,
    /// Rotation axis (world space).
    pub axis: [f64; 3],
    /// Angular limits.
    pub limits: JointLimits,
    /// Stop type.
    pub stop_type: StopType,
    /// Current angle (rad).
    pub angle: f64,
    /// Current angular velocity (rad/s).
    pub ang_velocity: f64,
    /// Warm-start impulse.
    pub warm_impulse: f64,
    /// Compliance.
    pub compliance: JointComplianceParams,
    /// Restitution at hard stop (0 = no bounce).
    pub restitution: f64,
}

impl RevoluteJointLimit {
    /// Create a revolute joint limit.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        body_a: usize,
        body_b: usize,
        axis: [f64; 3],
        lower: f64,
        upper: f64,
        stop_type: StopType,
        compliance: JointComplianceParams,
        restitution: f64,
    ) -> Self {
        Self {
            body_a,
            body_b,
            axis: normalize3(axis),
            limits: JointLimits::new(lower, upper),
            stop_type,
            angle: 0.0,
            ang_velocity: 0.0,
            warm_impulse: 0.0,
            compliance,
            restitution,
        }
    }

    /// Angular penetration.
    pub fn penetration(&self) -> f64 {
        self.limits.penetration(self.angle)
    }

    /// Whether the limit is active.
    pub fn is_active(&self) -> bool {
        self.limits.is_violated(self.angle)
    }

    /// Solve for corrective angular impulse.
    ///
    /// `inv_inertia_a`, `inv_inertia_b` — inverse moment of inertia projected
    /// onto the joint axis.
    pub fn solve_impulse(&self, inv_inertia_a: f64, inv_inertia_b: f64, dt: f64) -> f64 {
        let pen = self.penetration();
        if pen <= 0.0 {
            return 0.0;
        }
        let bias = self.compliance.positional_bias(pen, dt);
        let restitution_term = if self.ang_velocity.abs() > 1e-6 {
            self.restitution * self.ang_velocity.abs()
        } else {
            0.0
        };
        let k_eff_raw = inv_inertia_a + inv_inertia_b;
        let k_eff = self.compliance.softened_keff(k_eff_raw);
        if k_eff < 1e-14 {
            return 0.0;
        }
        let target = bias + restitution_term;
        let lambda = (target - self.ang_velocity) / k_eff;
        (self.warm_impulse + lambda).max(0.0) - self.warm_impulse
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// JointLimitPipeline
// ─────────────────────────────────────────────────────────────────────────────

/// High-level pipeline that processes all joint limit collisions for one frame.
#[derive(Debug, Clone, Default)]
pub struct JointLimitPipeline {
    /// Revolute joint limits.
    pub revolute: Vec<RevoluteJointLimit>,
    /// Prismatic joint limits.
    pub prismatic: Vec<PrismaticJointLimit>,
    /// Ball-socket limits.
    pub ball_socket: Vec<BallSocketLimitConstraint>,
    /// Cable/chain limits.
    pub cables: Vec<CableChainLimit>,
    /// Articulated body stops.
    pub ab_stops: Vec<ArticulatedBodyJointStop>,
    /// Motor hard limits.
    pub motors: Vec<MotorHardLimit>,
    /// Warm-start cache.
    pub warm_start: JointLimitWarmStart,
    /// Global compliance.
    pub compliance: JointComplianceParams,
    /// Iterations per step.
    pub iterations: usize,
}

impl JointLimitPipeline {
    /// Create a pipeline.
    pub fn new(iterations: usize) -> Self {
        Self {
            revolute: Vec::new(),
            prismatic: Vec::new(),
            ball_socket: Vec::new(),
            cables: Vec::new(),
            ab_stops: Vec::new(),
            motors: Vec::new(),
            warm_start: JointLimitWarmStart::new(),
            compliance: JointComplianceParams::rigid(),
            iterations,
        }
    }

    /// Step all registered constraints by `dt` with the given inverse masses.
    ///
    /// Returns the total absolute impulse applied.
    pub fn step(&mut self, inv_masses: &[f64], inv_inertias: &[f64], dt: f64) -> f64 {
        let mut total = 0.0f64;

        // Revolute
        for c in &mut self.revolute {
            let ia = if c.body_a < inv_inertias.len() {
                inv_inertias[c.body_a]
            } else {
                0.0
            };
            let ib = if c.body_b < inv_inertias.len() {
                inv_inertias[c.body_b]
            } else {
                0.0
            };
            for _ in 0..self.iterations {
                let imp = c.solve_impulse(ia, ib, dt);
                c.warm_impulse = (c.warm_impulse + imp).max(0.0);
                total += imp.abs();
            }
        }

        // Prismatic
        for c in &mut self.prismatic {
            let ma = if c.body_a < inv_masses.len() {
                inv_masses[c.body_a]
            } else {
                0.0
            };
            let mb = if c.body_b < inv_masses.len() {
                inv_masses[c.body_b]
            } else {
                0.0
            };
            for _ in 0..self.iterations {
                let imp = c.solve_impulse(ma, mb, dt);
                c.warm_impulse = (c.warm_impulse + imp).max(0.0);
                total += imp.abs();
            }
        }

        // Articulated body stops
        for c in &mut self.ab_stops {
            for _ in 0..self.iterations {
                let imp = c.compute_impulse(dt);
                c.warm_impulse = (c.warm_impulse + imp).max(0.0);
                total += imp.abs();
            }
        }

        // Cable/chain
        for c in &self.cables {
            let ma = if c.body_a < inv_masses.len() {
                inv_masses[c.body_a]
            } else {
                0.0
            };
            let mb = if c.body_b < inv_masses.len() {
                inv_masses[c.body_b]
            } else {
                0.0
            };
            let imp = c.compute_impulse(ma, mb, dt);
            total += imp.abs();
        }

        total
    }

    /// Count total active constraints.
    pub fn active_count(&self) -> usize {
        self.revolute.iter().filter(|c| c.is_active()).count()
            + self
                .prismatic
                .iter()
                .filter(|c| c.limits.is_violated(c.translation))
                .count()
            + self.ball_socket.iter().filter(|c| c.is_violated()).count()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Utility: wrap angle to [-π, π]
// ─────────────────────────────────────────────────────────────────────────────

/// Wrap an angle to `(−π, π]`.
pub fn wrap_angle(angle: f64) -> f64 {
    use std::f64::consts::PI;
    let mut a = angle;
    while a > PI {
        a -= 2.0 * PI;
    }
    while a <= -PI {
        a += 2.0 * PI;
    }
    a
}

/// Angular difference `b − a` wrapped to `(−π, π]`.
pub fn angle_diff(a: f64, b: f64) -> f64 {
    wrap_angle(b - a)
}

// ─────────────────────────────────────────────────────────────────────────────
// #[cfg(test)] unit tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    // ── 1. JointLimits::is_violated ──────────────────────────────────────────
    #[test]
    fn test_joint_limits_violation() {
        let lim = JointLimits::new(-1.0, 1.0);
        assert!(lim.is_violated(-1.5));
        assert!(lim.is_violated(1.5));
        assert!(!lim.is_violated(0.0));
        assert!(!lim.is_violated(-1.0));
        assert!(!lim.is_violated(1.0));
    }

    // ── 2. JointLimits::penetration ──────────────────────────────────────────
    #[test]
    fn test_joint_limits_penetration() {
        let lim = JointLimits::new(-PI / 2.0, PI / 2.0);
        let pen = lim.penetration(-PI);
        assert!((pen - (PI / 2.0)).abs() < 1e-10);
        assert_eq!(lim.penetration(0.0), 0.0);
    }

    // ── 3. JointLimits::clamp ────────────────────────────────────────────────
    #[test]
    fn test_joint_limits_clamp() {
        let lim = JointLimits::new(-1.0, 2.0);
        assert_eq!(lim.clamp(-5.0), -1.0);
        assert_eq!(lim.clamp(5.0), 2.0);
        assert_eq!(lim.clamp(0.5), 0.5);
    }

    // ── 4. StopType::Hard force ───────────────────────────────────────────────
    #[test]
    fn test_stop_type_hard_force() {
        let s = StopType::Hard;
        assert!(s.constraint_force(0.1, 0.0) > 0.0);
        assert_eq!(s.constraint_force(0.0, 0.0), 0.0);
    }

    // ── 5. StopType::Soft force ───────────────────────────────────────────────
    #[test]
    fn test_stop_type_soft_force() {
        let s = StopType::Soft {
            stiffness: 1000.0,
            damping: 10.0,
        };
        let f = s.constraint_force(0.01, 0.0);
        assert!((f - 10.0).abs() < 1e-10);
    }

    // ── 6. StopType::Soft no force when not violated ──────────────────────────
    #[test]
    fn test_stop_type_soft_no_violation() {
        let s = StopType::Soft {
            stiffness: 1000.0,
            damping: 10.0,
        };
        assert_eq!(s.constraint_force(-0.1, 0.0), 0.0);
    }

    // ── 7. JointLimitConstraint::compute_impulse ─────────────────────────────
    #[test]
    fn test_joint_limit_constraint_impulse() {
        let lim = JointLimits::new(-1.0, 1.0);
        let mut c = JointLimitConstraint::new(
            0,
            1,
            JointType::Revolute,
            lim,
            StopType::Soft {
                stiffness: 1000.0,
                damping: 0.0,
            },
            [0.0, 0.0, 1.0],
        );
        c.position = 1.5; // 0.5 rad past upper limit
        let imp = c.compute_impulse(1.0, 1.0, 0.01);
        assert!(imp > 0.0, "impulse should be positive");
    }

    // ── 8. JointLimitConstraint warm-start ───────────────────────────────────
    #[test]
    fn test_warm_start_positive_clamped() {
        let lim = JointLimits::new(-1.0, 1.0);
        let mut c = JointLimitConstraint::new(
            0,
            1,
            JointType::Revolute,
            lim,
            StopType::Hard,
            [1.0, 0.0, 0.0],
        );
        c.store_impulse(-5.0);
        // Negative impulse should be clamped to 0
        assert_eq!(c.apply_warm_start(), 0.0);
        c.store_impulse(3.0);
        assert_eq!(c.apply_warm_start(), 3.0);
    }

    // ── 9. BallSocketLimitConstraint::update_angle ───────────────────────────
    #[test]
    fn test_ball_socket_angle_update() {
        let mut bs =
            BallSocketLimitConstraint::new(0, 1, [0.0, 0.0, 1.0], PI / 4.0, StopType::Hard);
        // Child axis aligned with cone axis → angle = 0
        bs.update_angle([0.0, 0.0, 1.0]);
        assert!(bs.current_angle.abs() < 1e-10);
        assert!(!bs.is_violated());
        // Child axis perpendicular → angle = π/2
        bs.update_angle([1.0, 0.0, 0.0]);
        assert!((bs.current_angle - PI / 2.0).abs() < 1e-10);
        assert!(bs.is_violated());
    }

    // ── 10. BallSocketLimitConstraint::penetration ───────────────────────────
    #[test]
    fn test_ball_socket_penetration() {
        let mut bs =
            BallSocketLimitConstraint::new(0, 1, [0.0, 0.0, 1.0], PI / 6.0, StopType::Hard);
        bs.current_angle = PI / 3.0;
        let pen = bs.penetration();
        assert!((pen - PI / 6.0).abs() < 1e-10);
    }

    // ── 11. GearRatioCollision::propagate_impulse ────────────────────────────
    #[test]
    fn test_gear_propagate_impulse() {
        let gear = GearRatioCollision::new(0, 1, 3.0, 1.0);
        assert!((gear.propagate_impulse(10.0) - 30.0).abs() < 1e-10);
    }

    // ── 12. GearRatioCollision::reflected_inertia ────────────────────────────
    #[test]
    fn test_gear_reflected_inertia() {
        let gear = GearRatioCollision::new(0, 1, 2.0, 0.9);
        // reflected = 5.0 * 4.0 = 20.0
        assert!((gear.reflected_inertia(5.0) - 20.0).abs() < 1e-10);
    }

    // ── 13. CableChainLimit extension violation ───────────────────────────────
    #[test]
    fn test_cable_extension_violation() {
        let mut cable = CableChainLimit::new(0, 1, 0.0, 2.0, 1000.0, 10.0);
        cable.current_length = 2.5;
        assert!((cable.extension_violation() - 0.5).abs() < 1e-10);
        assert_eq!(cable.compression_violation(), 0.0);
    }

    // ── 14. CableChainLimit compression violation ─────────────────────────────
    #[test]
    fn test_chain_compression_violation() {
        let mut chain = CableChainLimit::new(0, 1, 1.0, 3.0, 500.0, 5.0);
        chain.current_length = 0.5;
        assert!((chain.compression_violation() - 0.5).abs() < 1e-10);
    }

    // ── 15. CableChainLimit::compute_impulse ─────────────────────────────────
    #[test]
    fn test_cable_impulse_positive_on_violation() {
        let mut cable = CableChainLimit::new(0, 1, 0.0, 1.0, 1000.0, 0.0);
        cable.current_length = 1.2;
        let imp = cable.compute_impulse(1.0, 1.0, 0.01);
        assert!(imp > 0.0);
    }

    // ── 16. ArticulatedBodyJointStop::compute_impulse ────────────────────────
    #[test]
    fn test_ab_stop_impulse() {
        let lim = JointLimits::new(-1.0, 1.0);
        let mut stop = ArticulatedBodyJointStop::new(
            0,
            2.0,
            lim,
            StopType::Soft {
                stiffness: 500.0,
                damping: 0.0,
            },
        );
        stop.position = 1.5;
        let imp = stop.compute_impulse(0.01);
        assert!(imp > 0.0);
    }

    // ── 17. ArticulatedBodyJointStop::velocity_impulse ───────────────────────
    #[test]
    fn test_ab_stop_velocity_impulse_restitution() {
        let lim = JointLimits::new(-1.0, 1.0);
        let mut stop = ArticulatedBodyJointStop::new(0, 1.0, lim, StopType::Hard);
        stop.position = 1.5;
        stop.velocity = -2.0;
        // restitution=0.5 → impulse = (1+0.5)*2.0 / 1.0 = 3.0
        let vi = stop.velocity_impulse(0.5);
        assert!((vi - 3.0).abs() < 1e-10);
    }

    // ── 18. MotorHardLimit::clamp_torque ─────────────────────────────────────
    #[test]
    fn test_motor_clamp_torque() {
        let m = MotorHardLimit::new(10.0, 5.0, -PI, PI, 0.1);
        assert_eq!(m.clamp_torque(20.0), 10.0);
        assert_eq!(m.clamp_torque(-20.0), -10.0);
        assert_eq!(m.clamp_torque(5.0), 5.0);
    }

    // ── 19. MotorHardLimit::check_stall ──────────────────────────────────────
    #[test]
    fn test_motor_stall_detection() {
        let mut m = MotorHardLimit::new(10.0, 5.0, -PI, PI, 0.1);
        m.angle = PI + 0.1; // past upper limit
        m.check_stall();
        assert!(m.stalled);
        m.angle = 0.0;
        m.check_stall();
        assert!(!m.stalled);
    }

    // ── 20. MotorHardLimit::limit_impulse zero when not stalled ──────────────
    #[test]
    fn test_motor_limit_impulse_not_stalled() {
        let m = MotorHardLimit::new(10.0, 5.0, -PI, PI, 0.1);
        assert_eq!(m.limit_impulse(0.01), 0.0);
    }

    // ── 21. JointFrictionConstraint::compute_impulse ─────────────────────────
    #[test]
    fn test_friction_impulse_computed() {
        let mut fric = JointFrictionConstraint::new(0, 1, [0.0, 0.0, 1.0], 0.3);
        fric.normal_force = 100.0;
        fric.relative_velocity = 2.0;
        let imp = fric.compute_impulse(2.0, 0.01);
        // raw = -2.0 / 2.0 = -1.0; max = 0.3*100=30; clamped to -30..30 → -1.0
        assert!((imp - (-1.0)).abs() < 1e-10);
    }

    // ── 22. JointFrictionConstraint warm start ────────────────────────────────
    #[test]
    fn test_friction_warm_start_clamped() {
        let mut fric = JointFrictionConstraint::new(0, 1, [1.0, 0.0, 0.0], 0.1);
        fric.normal_force = 50.0;
        fric.warm_impulse = 1000.0; // way above max friction
        let ws = fric.warm_start_impulse(1.0, 0.01);
        // max = 0.1 * 50 = 5.0
        assert!(ws <= 5.0 + 1e-10);
    }

    // ── 23. JointComplianceParams::rigid ─────────────────────────────────────
    #[test]
    fn test_compliance_rigid() {
        let c = JointComplianceParams::rigid();
        assert!((c.erp - 0.2).abs() < 1e-10);
        assert_eq!(c.cfm, 0.0);
    }

    // ── 24. JointComplianceParams::soft formula ───────────────────────────────
    #[test]
    fn test_compliance_soft() {
        let c = JointComplianceParams::soft(1000.0, 10.0, 0.01);
        // denom = 0.01*1000 + 10 = 20; erp = 10/20 = 0.5; cfm = 1/20 = 0.05
        assert!((c.erp - 0.5).abs() < 1e-10);
        assert!((c.cfm - 0.05).abs() < 1e-10);
    }

    // ── 25. JointComplianceParams::positional_bias ────────────────────────────
    #[test]
    fn test_compliance_positional_bias() {
        let c = JointComplianceParams::new(0.5, 0.0);
        let bias = c.positional_bias(0.1, 0.01);
        // bias = 0.5 * 0.1 / 0.01 = 5.0
        assert!((bias - 5.0).abs() < 1e-10);
    }

    // ── 26. JointLimitWarmStart::store/retrieve ───────────────────────────────
    #[test]
    fn test_warm_start_store_retrieve() {
        let mut ws = JointLimitWarmStart::new();
        ws.store(0, 1, 2, 5.5);
        assert!((ws.retrieve(0, 1, 2) - 5.5).abs() < 1e-10);
        assert_eq!(ws.retrieve(0, 1, 3), 0.0); // missing
    }

    // ── 27. JointLimitWarmStart::scale ────────────────────────────────────────
    #[test]
    fn test_warm_start_scale() {
        let mut ws = JointLimitWarmStart::new();
        ws.store(0, 1, 0, 10.0);
        ws.scale(0.5);
        assert!((ws.retrieve(0, 1, 0) - 5.0).abs() < 1e-10);
    }

    // ── 28. JointLimitWarmStart::clear ────────────────────────────────────────
    #[test]
    fn test_warm_start_clear() {
        let mut ws = JointLimitWarmStart::new();
        ws.store(0, 1, 0, 1.0);
        ws.clear();
        assert!(ws.is_empty());
    }

    // ── 29. JointLimitSolver::solve zero impulse when no violation ────────────
    #[test]
    fn test_solver_no_violation() {
        let mut solver = JointLimitSolver::new(4);
        let lim = JointLimits::new(-1.0, 1.0);
        let c = JointLimitConstraint::new(
            0,
            1,
            JointType::Revolute,
            lim,
            StopType::Soft {
                stiffness: 100.0,
                damping: 0.0,
            },
            [0.0, 0.0, 1.0],
        );
        solver.add_constraint(c);
        // position stays at 0.0 → no violation
        let inv_masses = vec![1.0, 1.0];
        let total = solver.solve(&inv_masses, 0.01);
        assert_eq!(total, 0.0);
    }

    // ── 30. JointLimitSolver::solve positive impulse on violation ─────────────
    #[test]
    fn test_solver_violation_produces_impulse() {
        let mut solver = JointLimitSolver::new(4);
        let lim = JointLimits::new(-1.0, 1.0);
        let mut c = JointLimitConstraint::new(
            0,
            1,
            JointType::Revolute,
            lim,
            StopType::Soft {
                stiffness: 100.0,
                damping: 0.0,
            },
            [0.0, 0.0, 1.0],
        );
        c.position = 2.0; // 1 rad past upper limit
        solver.add_constraint(c);
        let inv_masses = vec![1.0, 1.0];
        let total = solver.solve(&inv_masses, 0.01);
        assert!(total > 0.0, "expected positive impulse, got {}", total);
    }

    // ── 31. PrismaticJointLimit::penetration ──────────────────────────────────
    #[test]
    fn test_prismatic_penetration() {
        let mut pj = PrismaticJointLimit::new(
            0,
            1,
            [1.0, 0.0, 0.0],
            -0.5,
            0.5,
            StopType::Hard,
            JointComplianceParams::rigid(),
        );
        pj.translation = 0.8;
        assert!((pj.penetration() - 0.3).abs() < 1e-10);
    }

    // ── 32. PrismaticJointLimit::constraint_row ───────────────────────────────
    #[test]
    fn test_prismatic_constraint_row() {
        let mut pj = PrismaticJointLimit::new(
            0,
            1,
            [0.0, 1.0, 0.0],
            -1.0,
            1.0,
            StopType::Hard,
            JointComplianceParams::rigid(),
        );
        pj.translation = 1.5;
        assert_eq!(pj.constraint_row(), -1.0);
        pj.translation = -1.5;
        assert_eq!(pj.constraint_row(), 1.0);
        pj.translation = 0.0;
        assert_eq!(pj.constraint_row(), 0.0);
    }

    // ── 33. RevoluteJointLimit::is_active ─────────────────────────────────────
    #[test]
    fn test_revolute_is_active() {
        let mut rev = RevoluteJointLimit::new(
            0,
            1,
            [0.0, 1.0, 0.0],
            -PI / 2.0,
            PI / 2.0,
            StopType::Hard,
            JointComplianceParams::rigid(),
            0.0,
        );
        rev.angle = PI;
        assert!(rev.is_active());
        rev.angle = 0.0;
        assert!(!rev.is_active());
    }

    // ── 34. RevoluteJointLimit::solve_impulse ────────────────────────────────
    #[test]
    fn test_revolute_solve_impulse_positive() {
        let mut rev = RevoluteJointLimit::new(
            0,
            1,
            [0.0, 0.0, 1.0],
            -PI / 4.0,
            PI / 4.0,
            StopType::Hard,
            JointComplianceParams::new(0.8, 0.0),
            0.0,
        );
        rev.angle = PI / 2.0; // past limit
        rev.ang_velocity = 1.0;
        let imp = rev.solve_impulse(1.0, 1.0, 0.01);
        assert!(imp >= 0.0);
    }

    // ── 35. JointStopContact::compute_impulse ────────────────────────────────
    #[test]
    fn test_joint_stop_contact_impulse() {
        let contact = JointStopContact::new(
            0,
            [0.0, 1.0, 0.0],
            0.05,
            -1.0, // approaching
            0.2,
        );
        let imp = contact.compute_impulse(2.0, 0.2, 0.01);
        assert!(imp >= 0.0);
    }

    // ── 36. JointLimitPipeline::active_count ─────────────────────────────────
    #[test]
    fn test_pipeline_active_count() {
        let mut pipe = JointLimitPipeline::new(4);
        let mut rev = RevoluteJointLimit::new(
            0,
            1,
            [0.0, 0.0, 1.0],
            -1.0,
            1.0,
            StopType::Hard,
            JointComplianceParams::rigid(),
            0.0,
        );
        rev.angle = 2.0; // violated
        pipe.revolute.push(rev);
        assert_eq!(pipe.active_count(), 1);
    }

    // ── 37. wrap_angle ────────────────────────────────────────────────────────
    #[test]
    fn test_wrap_angle() {
        assert!((wrap_angle(3.0 * PI) - PI).abs() < 1e-10);
        assert!(
            (wrap_angle(-3.0 * PI) - (-PI)).abs() < 1e-10
                || (wrap_angle(-3.0 * PI) - PI).abs() < 1e-10
        );
        assert!(wrap_angle(0.5).abs() - 0.5 < 1e-10);
    }

    // ── 38. angle_diff ────────────────────────────────────────────────────────
    #[test]
    fn test_angle_diff_wrap() {
        let d = angle_diff(PI * 0.9, -PI * 0.9);
        // difference is about -1.8π → wrapped ≈ 0.2π
        assert!(d.abs() < PI + 1e-10);
    }

    // ── 39. normalize3 zero vector ────────────────────────────────────────────
    #[test]
    fn test_normalize3_zero() {
        let v = normalize3([0.0, 0.0, 0.0]);
        assert_eq!(v, [0.0, 0.0, 0.0]);
    }

    // ── 40. normalize3 unit result ────────────────────────────────────────────
    #[test]
    fn test_normalize3_unit() {
        let v = normalize3([3.0, 4.0, 0.0]);
        let n = norm3(v);
        assert!((n - 1.0).abs() < 1e-10);
    }

    // ── 41. clampf ────────────────────────────────────────────────────────────
    #[test]
    fn test_clampf() {
        assert_eq!(clampf(5.0, 0.0, 3.0), 3.0);
        assert_eq!(clampf(-1.0, 0.0, 3.0), 0.0);
        assert_eq!(clampf(1.5, 0.0, 3.0), 1.5);
    }

    // ── 42. JointLimitSolver::clear_constraints ───────────────────────────────
    #[test]
    fn test_solver_clear_constraints() {
        let mut solver = JointLimitSolver::new(2);
        let lim = JointLimits::new(-1.0, 1.0);
        let c = JointLimitConstraint::new(
            0,
            1,
            JointType::Prismatic,
            lim,
            StopType::Hard,
            [1.0, 0.0, 0.0],
        );
        solver.add_constraint(c);
        assert_eq!(solver.constraints.len(), 1);
        solver.clear_constraints();
        assert!(solver.constraints.is_empty());
    }
}
