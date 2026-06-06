// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Penalty / soft constraint methods for rigid body simulation.
//!
//! This module implements soft constraints using the penalty method:
//! constraint violations generate restoring forces proportional to
//! the magnitude of the violation.
//!
//! # Overview
//!
//! - [`SoftConstraint`] — trait for computing penalty forces/torques
//! - [`SpringConstraint`] — spring attached from body to a fixed anchor
//! - [`DistanceConstraint`] — maintain a target distance between two bodies
//! - [`SoftJoint`] — position/orientation target with separate stiffness
//! - [`ConstraintSolver`] — aggregate solver applying all constraints
//! - [`PenaltyMethod`] — helper for violation-to-force conversion
//! - [`RigidBodyState`] — minimal state used by this module (no nalgebra)

// ---------------------------------------------------------------------------
// 3-D vector helpers (plain arrays, no nalgebra)
// ---------------------------------------------------------------------------

#[inline]
fn vec3_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn vec3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
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
fn vec3_norm(a: [f64; 3]) -> f64 {
    vec3_dot(a, a).sqrt()
}

#[inline]
fn vec3_normalize(a: [f64; 3]) -> [f64; 3] {
    let n = vec3_norm(a);
    if n > 1e-12 {
        vec3_scale(a, 1.0 / n)
    } else {
        [0.0; 3]
    }
}

#[inline]
fn vec3_neg(a: [f64; 3]) -> [f64; 3] {
    [-a[0], -a[1], -a[2]]
}

#[cfg(test)]
#[inline]
fn vec3_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

// ---------------------------------------------------------------------------
// Minimal rigid-body state (no nalgebra)
// ---------------------------------------------------------------------------

/// Minimal rigid-body state used by soft constraint calculations.
///
/// Contains position, velocity, angular velocity, and mass.
/// Orientation is stored as a quaternion `[w, x, y, z]`.
#[derive(Debug, Clone)]
pub struct RigidBodyState {
    /// World-space position of the centre of mass.
    pub position: [f64; 3],
    /// Linear velocity in world space.
    pub velocity: [f64; 3],
    /// Angular velocity in world space (rad/s).
    pub angular_velocity: [f64; 3],
    /// Orientation quaternion `[w, x, y, z]` (unit quaternion).
    pub orientation: [f64; 4],
    /// Body mass in kg.
    pub mass: f64,
}

impl RigidBodyState {
    /// Construct a state at the origin with zero velocity and unit mass.
    pub fn new() -> Self {
        Self {
            position: [0.0; 3],
            velocity: [0.0; 3],
            angular_velocity: [0.0; 3],
            orientation: [1.0, 0.0, 0.0, 0.0],
            mass: 1.0,
        }
    }

    /// Construct a state with given position and mass.
    pub fn with_position(position: [f64; 3], mass: f64) -> Self {
        Self {
            position,
            mass,
            ..Self::new()
        }
    }
}

impl Default for RigidBodyState {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// SoftConstraint trait
// ---------------------------------------------------------------------------

/// A soft (penalty-based) constraint that contributes forces and torques to a
/// rigid body.
///
/// Implementations compute a restoring force/torque proportional to the
/// constraint violation (penalty method).
pub trait SoftConstraint: Send + Sync {
    /// Compute the world-space force to apply to the body (N).
    fn compute_force(&self, state: &RigidBodyState) -> [f64; 3];

    /// Compute the world-space torque to apply to the body (N·m).
    fn compute_torque(&self, state: &RigidBodyState) -> [f64; 3];

    /// Combined call returning `(force, torque)`.
    fn compute(&self, state: &RigidBodyState) -> ([f64; 3], [f64; 3]) {
        (self.compute_force(state), self.compute_torque(state))
    }
}

// ---------------------------------------------------------------------------
// PenaltyMethod helper
// ---------------------------------------------------------------------------

/// Helper for converting constraint violations into penalty forces.
///
/// Computes a restoring force proportional to violation magnitude:
/// `F = -stiffness * violation - damping * velocity_along_violation`.
#[derive(Debug, Clone)]
pub struct PenaltyMethod {
    /// Spring stiffness coefficient (N/m).
    pub stiffness: f64,
    /// Damping coefficient (N·s/m).
    pub damping: f64,
}

impl PenaltyMethod {
    /// Create a new `PenaltyMethod` with given stiffness and damping.
    pub fn new(stiffness: f64, damping: f64) -> Self {
        Self { stiffness, damping }
    }

    /// Compute penalty force given a scalar violation and relative velocity
    /// along the constraint direction.
    ///
    /// `violation` is positive when the constraint is violated (stretched/compressed).
    /// `relative_velocity` is the component of velocity in the constraint direction.
    pub fn compute_scalar_force(&self, violation: f64, relative_velocity: f64) -> f64 {
        -self.stiffness * violation - self.damping * relative_velocity
    }

    /// Compute penalty force vector given a 3-D violation vector and
    /// the body's velocity.
    pub fn compute_vector_force(&self, violation: [f64; 3], velocity: [f64; 3]) -> [f64; 3] {
        let v_dot = vec3_dot(violation, velocity);
        let violation_len = vec3_norm(violation);
        if violation_len < 1e-12 {
            return [0.0; 3];
        }
        let spring = vec3_scale(violation, -self.stiffness);
        let damp_scalar = -self.damping * v_dot / violation_len;
        let dir = vec3_normalize(violation);
        let damp = vec3_scale(dir, damp_scalar);
        vec3_add(spring, damp)
    }
}

// ---------------------------------------------------------------------------
// SpringConstraint
// ---------------------------------------------------------------------------

/// A spring constraint that connects a body's centre of mass to a fixed world
/// anchor point.
///
/// Applies Hooke's law with optional damping:
/// `F = -k * (d - L0) * n̂ - c * (v · n̂) * n̂`
/// where `d` is the current distance, `L0` the rest length, and `n̂` the
/// unit direction from anchor to body.
#[derive(Debug, Clone)]
pub struct SpringConstraint {
    /// Fixed anchor position in world space.
    pub anchor: [f64; 3],
    /// Natural (rest) length of the spring (m).
    pub rest_length: f64,
    /// Spring stiffness (N/m).
    pub stiffness: f64,
    /// Damping coefficient (N·s/m).
    pub damping: f64,
}

impl SpringConstraint {
    /// Create a new `SpringConstraint`.
    pub fn new(anchor: [f64; 3], rest_length: f64, stiffness: f64, damping: f64) -> Self {
        Self {
            anchor,
            rest_length,
            stiffness,
            damping,
        }
    }
}

impl SoftConstraint for SpringConstraint {
    fn compute_force(&self, state: &RigidBodyState) -> [f64; 3] {
        let delta = vec3_sub(state.position, self.anchor);
        let dist = vec3_norm(delta);
        if dist < 1e-12 {
            return [0.0; 3];
        }
        let dir = vec3_scale(delta, 1.0 / dist);
        let extension = dist - self.rest_length;
        let v_along = vec3_dot(state.velocity, dir);
        let f_mag = -self.stiffness * extension - self.damping * v_along;
        vec3_scale(dir, f_mag)
    }

    fn compute_torque(&self, _state: &RigidBodyState) -> [f64; 3] {
        [0.0; 3]
    }
}

// ---------------------------------------------------------------------------
// DistanceConstraint
// ---------------------------------------------------------------------------

/// A soft distance constraint between two bodies.
///
/// Applies equal and opposite penalty forces to maintain `target_distance`
/// between the two centres of mass.
#[derive(Debug, Clone)]
pub struct DistanceConstraint {
    /// State of body A.
    pub body_a: RigidBodyState,
    /// State of body B.
    pub body_b: RigidBodyState,
    /// Target separation distance (m).
    pub target_distance: f64,
    /// Penalty stiffness (N/m).
    pub stiffness: f64,
}

impl DistanceConstraint {
    /// Create a new `DistanceConstraint`.
    pub fn new(
        body_a: RigidBodyState,
        body_b: RigidBodyState,
        target_distance: f64,
        stiffness: f64,
    ) -> Self {
        Self {
            body_a,
            body_b,
            target_distance,
            stiffness,
        }
    }

    /// Compute violation: actual distance minus target distance.
    pub fn violation(&self) -> f64 {
        let delta = vec3_sub(self.body_a.position, self.body_b.position);
        vec3_norm(delta) - self.target_distance
    }

    /// Compute the force applied to body A (body B receives the negation).
    pub fn force_on_a(&self) -> [f64; 3] {
        let delta = vec3_sub(self.body_a.position, self.body_b.position);
        let dist = vec3_norm(delta);
        if dist < 1e-12 {
            return [0.0; 3];
        }
        let violation = dist - self.target_distance;
        let dir = vec3_scale(delta, 1.0 / dist);
        vec3_scale(dir, -self.stiffness * violation)
    }

    /// Compute the force applied to body B.
    pub fn force_on_b(&self) -> [f64; 3] {
        vec3_neg(self.force_on_a())
    }
}

impl SoftConstraint for DistanceConstraint {
    fn compute_force(&self, _state: &RigidBodyState) -> [f64; 3] {
        self.force_on_a()
    }

    fn compute_torque(&self, _state: &RigidBodyState) -> [f64; 3] {
        [0.0; 3]
    }
}

// ---------------------------------------------------------------------------
// SoftJoint
// ---------------------------------------------------------------------------

/// A soft joint that drives a body toward a target position and orientation.
///
/// Applies separate linear and angular penalty forces.
#[derive(Debug, Clone)]
pub struct SoftJoint {
    /// Target world-space position.
    pub position_target: [f64; 3],
    /// Target orientation as unit quaternion `[w, x, y, z]`.
    pub orientation_target: [f64; 4],
    /// Stiffness for position error (N/m).
    pub linear_stiffness: f64,
    /// Stiffness for orientation error (N·m/rad).
    pub angular_stiffness: f64,
    /// Linear damping coefficient.
    pub linear_damping: f64,
    /// Angular damping coefficient.
    pub angular_damping: f64,
}

impl SoftJoint {
    /// Create a new `SoftJoint` with given targets and stiffnesses.
    pub fn new(
        position_target: [f64; 3],
        orientation_target: [f64; 4],
        linear_stiffness: f64,
        angular_stiffness: f64,
    ) -> Self {
        Self {
            position_target,
            orientation_target,
            linear_stiffness,
            angular_stiffness,
            linear_damping: 0.0,
            angular_damping: 0.0,
        }
    }

    /// Set damping values and return `self` for chaining.
    pub fn with_damping(mut self, linear_damping: f64, angular_damping: f64) -> Self {
        self.linear_damping = linear_damping;
        self.angular_damping = angular_damping;
        self
    }

    /// Compute the linear position error vector (target - current).
    pub fn position_error(&self, state: &RigidBodyState) -> [f64; 3] {
        vec3_sub(self.position_target, state.position)
    }

    /// Compute the angular error axis-angle from the current orientation to
    /// the target.  Returns the 3-D angular error vector (rad).
    pub fn orientation_error(&self, state: &RigidBodyState) -> [f64; 3] {
        // Compute relative quaternion: q_rel = q_target * q_current^{-1}
        // For unit quaternions: q^{-1} = [w, -x, -y, -z]
        let [qw, qx, qy, qz] = state.orientation;
        let [tw, tx, ty, tz] = self.orientation_target;
        // q_rel = q_target * conj(q_current)
        let rw = tw * qw + tx * qx + ty * qy + tz * qz;
        let rx = tx * qw - tw * qx + tz * qy - ty * qz;
        let ry = ty * qw - tz * qx - tw * qy + tx * qz;
        let rz = tz * qw + ty * qx - tx * qy - tw * qz;
        // Convert to axis-angle: angle = 2 * acos(rw), axis = [rx,ry,rz]/sin(angle/2)
        let half_angle = rw.clamp(-1.0, 1.0).acos();
        let sin_half = half_angle.sin();
        if sin_half.abs() < 1e-10 {
            [0.0; 3]
        } else {
            let angle = 2.0 * half_angle;
            let k = angle / sin_half;
            [rx * k, ry * k, rz * k]
        }
    }
}

impl SoftConstraint for SoftJoint {
    fn compute_force(&self, state: &RigidBodyState) -> [f64; 3] {
        let err = self.position_error(state);
        let spring = vec3_scale(err, self.linear_stiffness);
        let damp = vec3_scale(state.velocity, -self.linear_damping);
        vec3_add(spring, damp)
    }

    fn compute_torque(&self, state: &RigidBodyState) -> [f64; 3] {
        let err = self.orientation_error(state);
        let spring = vec3_scale(err, self.angular_stiffness);
        let damp = vec3_scale(state.angular_velocity, -self.angular_damping);
        vec3_add(spring, damp)
    }
}

// ---------------------------------------------------------------------------
// ConstraintSolver
// ---------------------------------------------------------------------------

/// Aggregate solver that applies all registered soft constraints to a body.
///
/// Each constraint in the list contributes penalty forces and torques which
/// are summed and returned.
pub struct ConstraintSolver {
    /// Collection of soft constraints to solve.
    pub constraints: Vec<Box<dyn SoftConstraint>>,
}

impl ConstraintSolver {
    /// Create an empty solver.
    pub fn new() -> Self {
        Self {
            constraints: Vec::new(),
        }
    }

    /// Add a constraint to the solver.
    pub fn add<C: SoftConstraint + 'static>(&mut self, constraint: C) {
        self.constraints.push(Box::new(constraint));
    }

    /// Apply all constraints to `state`, returning `(total_force, total_torque)`.
    pub fn solve(&self, state: &RigidBodyState) -> ([f64; 3], [f64; 3]) {
        let mut force = [0.0_f64; 3];
        let mut torque = [0.0_f64; 3];
        for c in &self.constraints {
            let (f, t) = c.compute(state);
            force = vec3_add(force, f);
            torque = vec3_add(torque, t);
        }
        (force, torque)
    }

    /// Return the number of registered constraints.
    pub fn len(&self) -> usize {
        self.constraints.len()
    }

    /// Return `true` if no constraints are registered.
    pub fn is_empty(&self) -> bool {
        self.constraints.is_empty()
    }
}

impl Default for ConstraintSolver {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ── helpers ─────────────────────────────────────────────────────────────

    fn approx_eq(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    fn approx_eq3(a: [f64; 3], b: [f64; 3], tol: f64) -> bool {
        approx_eq(a[0], b[0], tol) && approx_eq(a[1], b[1], tol) && approx_eq(a[2], b[2], tol)
    }

    fn state_at(pos: [f64; 3]) -> RigidBodyState {
        RigidBodyState::with_position(pos, 1.0)
    }

    // ── vec3 helpers ────────────────────────────────────────────────────────

    #[test]
    fn vec3_add_basic() {
        assert_eq!(vec3_add([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]), [5.0, 7.0, 9.0]);
    }

    #[test]
    fn vec3_sub_basic() {
        assert_eq!(vec3_sub([5.0, 5.0, 5.0], [2.0, 3.0, 1.0]), [3.0, 2.0, 4.0]);
    }

    #[test]
    fn vec3_norm_unit_x() {
        assert!(approx_eq(vec3_norm([1.0, 0.0, 0.0]), 1.0, 1e-12));
    }

    #[test]
    fn vec3_norm_pythagoras() {
        assert!(approx_eq(vec3_norm([3.0, 4.0, 0.0]), 5.0, 1e-10));
    }

    #[test]
    fn vec3_normalize_gives_unit() {
        let n = vec3_normalize([3.0, 4.0, 0.0]);
        assert!(approx_eq(vec3_norm(n), 1.0, 1e-10));
    }

    #[test]
    fn vec3_normalize_zero_returns_zero() {
        assert_eq!(vec3_normalize([0.0; 3]), [0.0; 3]);
    }

    #[test]
    fn vec3_cross_basic() {
        let c = vec3_cross([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert!(approx_eq3(c, [0.0, 0.0, 1.0], 1e-12));
    }

    // ── RigidBodyState ──────────────────────────────────────────────────────

    #[test]
    fn rigid_body_state_default_position() {
        let s = RigidBodyState::new();
        assert_eq!(s.position, [0.0; 3]);
    }

    #[test]
    fn rigid_body_state_with_position() {
        let s = RigidBodyState::with_position([1.0, 2.0, 3.0], 5.0);
        assert_eq!(s.position, [1.0, 2.0, 3.0]);
        assert!(approx_eq(s.mass, 5.0, 1e-12));
    }

    // ── PenaltyMethod ────────────────────────────────────────────────────────

    #[test]
    fn penalty_method_scalar_no_damping() {
        let pm = PenaltyMethod::new(100.0, 0.0);
        // violation = 0.5 m, v = 0 → F = -100 * 0.5 = -50
        let f = pm.compute_scalar_force(0.5, 0.0);
        assert!(approx_eq(f, -50.0, 1e-10));
    }

    #[test]
    fn penalty_method_scalar_with_damping() {
        let pm = PenaltyMethod::new(100.0, 10.0);
        // violation = 0.1 m, velocity = 2.0 → F = -10 - 20 = -30
        let f = pm.compute_scalar_force(0.1, 2.0);
        assert!(approx_eq(f, -30.0, 1e-10));
    }

    #[test]
    fn penalty_method_vector_force_along_x() {
        let pm = PenaltyMethod::new(100.0, 0.0);
        let violation = [0.1, 0.0, 0.0]; // 0.1 m in X
        let velocity = [0.0; 3];
        let f = pm.compute_vector_force(violation, velocity);
        // F = -k * violation = [-10, 0, 0]
        assert!(approx_eq(f[0], -10.0, 1e-10));
        assert!(approx_eq(f[1], 0.0, 1e-10));
        assert!(approx_eq(f[2], 0.0, 1e-10));
    }

    #[test]
    fn penalty_method_vector_zero_violation() {
        let pm = PenaltyMethod::new(100.0, 10.0);
        let f = pm.compute_vector_force([0.0; 3], [1.0, 0.0, 0.0]);
        assert_eq!(f, [0.0; 3]);
    }

    // ── SpringConstraint ─────────────────────────────────────────────────────

    #[test]
    fn spring_at_rest_length_produces_zero_force() {
        let spring = SpringConstraint::new([0.0; 3], 1.0, 500.0, 0.0);
        let state = state_at([1.0, 0.0, 0.0]); // exactly 1 m from anchor
        let f = spring.compute_force(&state);
        assert!(approx_eq3(f, [0.0; 3], 1e-10));
    }

    #[test]
    fn spring_stretched_produces_restoring_force() {
        let spring = SpringConstraint::new([0.0; 3], 1.0, 100.0, 0.0);
        // Body at 2.0 m, rest length 1.0 → extension = 1.0 → F = -100 N in +X
        let state = state_at([2.0, 0.0, 0.0]);
        let f = spring.compute_force(&state);
        assert!(approx_eq(f[0], -100.0, 1e-8));
        assert!(approx_eq(f[1], 0.0, 1e-10));
    }

    #[test]
    fn spring_compressed_produces_repulsive_force() {
        let spring = SpringConstraint::new([0.0; 3], 2.0, 100.0, 0.0);
        // Body at 1.0 m, rest length 2.0 → extension = -1.0 → F = +100 N in +X
        let state = state_at([1.0, 0.0, 0.0]);
        let f = spring.compute_force(&state);
        assert!(approx_eq(f[0], 100.0, 1e-8));
    }

    #[test]
    fn spring_torque_is_zero() {
        let spring = SpringConstraint::new([0.0; 3], 1.0, 100.0, 0.0);
        let state = state_at([2.0, 0.0, 0.0]);
        let t = spring.compute_torque(&state);
        assert_eq!(t, [0.0; 3]);
    }

    #[test]
    fn spring_damping_reduces_velocity_force() {
        let spring = SpringConstraint::new([0.0; 3], 1.0, 100.0, 20.0);
        let mut state = state_at([2.0, 0.0, 0.0]);
        // Body moving away from anchor at 1.0 m/s
        state.velocity = [1.0, 0.0, 0.0];
        let f = spring.compute_force(&state);
        // spring: -100 N, damping: -20 N → total = -120 N
        assert!(approx_eq(f[0], -120.0, 1e-8));
    }

    #[test]
    fn spring_zero_distance_no_panic() {
        let spring = SpringConstraint::new([0.0; 3], 1.0, 100.0, 10.0);
        let state = state_at([0.0; 3]); // at anchor
        let f = spring.compute_force(&state);
        assert_eq!(f, [0.0; 3]);
    }

    #[test]
    fn spring_diagonal_anchor() {
        let spring = SpringConstraint::new([1.0, 1.0, 0.0], 0.0, 100.0, 0.0);
        // Body at origin, anchor at (1,1,0), rest_length=0 → dist=sqrt(2)
        // extension = sqrt(2) → force toward anchor = -100 * sqrt(2) * dir
        let state = state_at([0.0; 3]);
        let f = spring.compute_force(&state);
        // dir from anchor to body is [-1,-1,0]/sqrt(2)
        let expected = 100.0 * 2.0_f64.sqrt();
        let fmag = vec3_norm(f);
        assert!(approx_eq(fmag, expected, 1e-6));
    }

    // ── DistanceConstraint ───────────────────────────────────────────────────

    #[test]
    fn distance_constraint_at_target_no_force() {
        let a = RigidBodyState::with_position([0.0, 0.0, 0.0], 1.0);
        let b = RigidBodyState::with_position([1.0, 0.0, 0.0], 1.0);
        let dc = DistanceConstraint::new(a, b, 1.0, 500.0);
        let fa = dc.force_on_a();
        assert!(approx_eq3(fa, [0.0; 3], 1e-10));
    }

    #[test]
    fn distance_constraint_violation_positive_when_far() {
        let a = RigidBodyState::with_position([0.0, 0.0, 0.0], 1.0);
        let b = RigidBodyState::with_position([3.0, 0.0, 0.0], 1.0);
        let dc = DistanceConstraint::new(a, b, 1.0, 100.0);
        assert!(approx_eq(dc.violation(), 2.0, 1e-10));
    }

    #[test]
    fn distance_constraint_forces_are_equal_opposite() {
        let a = RigidBodyState::with_position([0.0, 0.0, 0.0], 1.0);
        let b = RigidBodyState::with_position([2.0, 0.0, 0.0], 1.0);
        let dc = DistanceConstraint::new(a, b, 1.0, 100.0);
        let fa = dc.force_on_a();
        let fb = dc.force_on_b();
        assert!(approx_eq(fa[0], -fb[0], 1e-10));
        assert!(approx_eq(fa[1], -fb[1], 1e-10));
        assert!(approx_eq(fa[2], -fb[2], 1e-10));
    }

    #[test]
    fn distance_constraint_force_attracts_when_too_far() {
        let a = RigidBodyState::with_position([0.0, 0.0, 0.0], 1.0);
        let b = RigidBodyState::with_position([3.0, 0.0, 0.0], 1.0);
        // target = 1.0, actual = 3.0, violation = 2.0
        // force_on_a should pull toward b → positive X
        let dc = DistanceConstraint::new(a, b, 1.0, 100.0);
        let fa = dc.force_on_a();
        // violation = 2.0, dir from b to a is -X → force_on_a = -k*violation*dir = -100*2*(-1) = +200
        // But force_on_a uses delta = a - b = -3x, dir = -x, violation = 2 → F = -k*v*dir = -100*2*(-1,0,0) = (+200,0,0)
        // Wait: delta = a-b = (-3,0,0), dist=3, violation=3-1=2, dir=(-1,0,0), F=-k*v*dir = -100*2*(-1,0,0) = (200,0,0)
        // Hmm: that repels. Let me recheck the sign in force_on_a:
        // force_on_a = scale(dir, -k*violation) = scale((-1,0,0), -100*2) = (200,0,0)
        // That pushes a away from b which is wrong for attraction. But the implementation is correct:
        // when bodies are too far apart, violation > 0, dir points from b to a (outward),
        // and F = -k*violation * dir is inward (toward b). Let me re-check:
        // delta = a.pos - b.pos = -3x → dist=3, dir = delta/dist = -x
        // F = scale(dir, -k*violation) = scale(-x, -100*2) = +200x  ← pushes a in +x (away from b at +3x)
        // That seems wrong. But re-reading: a is at 0, b is at +3. delta = 0-3 = -3x. dir=-x.
        // -k*v*dir = -100*2*(-x) = +200x pushes a in +x direction (AWAY from b). This is wrong physically
        // but the test should just verify the implementation's sign.
        // Actually for a distance constraint keeping bodies at target, when dist>target the bodies should attract.
        // Let's check: we want force on a to pull it toward b (+x direction): so +200x is correct here!
        // a is at 0, b is at +3. Force on a should be +x (toward b). F=(200,0,0) → correct.
        assert!(fa[0] > 0.0, "force on a should attract toward b: fa={fa:?}");
    }

    #[test]
    fn distance_constraint_soft_constraint_trait() {
        let a = RigidBodyState::with_position([0.0, 0.0, 0.0], 1.0);
        let b = RigidBodyState::with_position([2.0, 0.0, 0.0], 1.0);
        let dc = DistanceConstraint::new(a.clone(), b, 1.0, 100.0);
        let f = dc.compute_force(&a);
        assert!(f[0].is_finite());
    }

    // ── SoftJoint ────────────────────────────────────────────────────────────

    #[test]
    fn soft_joint_at_target_zero_force() {
        let joint = SoftJoint::new([1.0, 2.0, 3.0], [1.0, 0.0, 0.0, 0.0], 500.0, 100.0);
        let state = RigidBodyState::with_position([1.0, 2.0, 3.0], 1.0);
        let f = joint.compute_force(&state);
        assert!(approx_eq3(f, [0.0; 3], 1e-10));
    }

    #[test]
    fn soft_joint_displaced_produces_force() {
        let joint = SoftJoint::new([0.0; 3], [1.0, 0.0, 0.0, 0.0], 100.0, 0.0);
        let state = RigidBodyState::with_position([1.0, 0.0, 0.0], 1.0);
        let f = joint.compute_force(&state);
        // error = target - pos = -1 x → force = 100 * -1 = -100 x
        assert!(approx_eq(f[0], -100.0, 1e-10));
    }

    #[test]
    fn soft_joint_identity_orientation_zero_torque() {
        let joint = SoftJoint::new([0.0; 3], [1.0, 0.0, 0.0, 0.0], 100.0, 100.0);
        let state = RigidBodyState::new(); // identity orientation
        let t = joint.compute_torque(&state);
        assert!(approx_eq3(t, [0.0; 3], 1e-8));
    }

    #[test]
    fn soft_joint_damping_opposes_velocity() {
        let joint =
            SoftJoint::new([0.0; 3], [1.0, 0.0, 0.0, 0.0], 0.0, 0.0).with_damping(10.0, 0.0);
        let mut state = RigidBodyState::with_position([0.0; 3], 1.0);
        state.velocity = [5.0, 0.0, 0.0];
        let f = joint.compute_force(&state);
        // damp = -10 * 5 = -50 in X
        assert!(approx_eq(f[0], -50.0, 1e-10));
    }

    #[test]
    fn soft_joint_angular_damping_opposes_spin() {
        let joint = SoftJoint::new([0.0; 3], [1.0, 0.0, 0.0, 0.0], 0.0, 0.0).with_damping(0.0, 5.0);
        let mut state = RigidBodyState::new();
        state.angular_velocity = [0.0, 2.0, 0.0];
        let t = joint.compute_torque(&state);
        // angular damp = -5 * 2 = -10 in Y
        assert!(approx_eq(t[1], -10.0, 1e-10));
    }

    #[test]
    fn soft_joint_position_error_correct() {
        let joint = SoftJoint::new([3.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0], 100.0, 0.0);
        let state = state_at([1.0, 0.0, 0.0]);
        let err = joint.position_error(&state);
        assert!(approx_eq(err[0], 2.0, 1e-10));
    }

    // ── ConstraintSolver ─────────────────────────────────────────────────────

    #[test]
    fn solver_empty_gives_zero() {
        let solver = ConstraintSolver::new();
        let state = RigidBodyState::new();
        let (f, t) = solver.solve(&state);
        assert_eq!(f, [0.0; 3]);
        assert_eq!(t, [0.0; 3]);
    }

    #[test]
    fn solver_single_spring() {
        let mut solver = ConstraintSolver::new();
        solver.add(SpringConstraint::new([0.0; 3], 0.0, 100.0, 0.0));
        let state = state_at([1.0, 0.0, 0.0]);
        let (f, _) = solver.solve(&state);
        assert!(approx_eq(f[0], -100.0, 1e-8));
    }

    #[test]
    fn solver_multiple_constraints_sum() {
        let mut solver = ConstraintSolver::new();
        // Two springs from origin, both with stiffness 100 and rest 0
        solver.add(SpringConstraint::new([0.0; 3], 0.0, 100.0, 0.0));
        solver.add(SpringConstraint::new([0.0; 3], 0.0, 100.0, 0.0));
        let state = state_at([1.0, 0.0, 0.0]);
        let (f, _) = solver.solve(&state);
        assert!(approx_eq(f[0], -200.0, 1e-8));
    }

    #[test]
    fn solver_len_and_is_empty() {
        let mut solver = ConstraintSolver::new();
        assert!(solver.is_empty());
        solver.add(SpringConstraint::new([0.0; 3], 1.0, 100.0, 0.0));
        assert_eq!(solver.len(), 1);
        assert!(!solver.is_empty());
    }

    #[test]
    fn solver_spring_plus_joint() {
        let mut solver = ConstraintSolver::new();
        solver.add(SpringConstraint::new([0.0; 3], 0.0, 50.0, 0.0));
        solver.add(SoftJoint::new([0.0; 3], [1.0, 0.0, 0.0, 0.0], 50.0, 0.0));
        let state = state_at([1.0, 0.0, 0.0]);
        let (f, _) = solver.solve(&state);
        // spring: -50, joint: 50*(0-1) = -50 → total = -100
        assert!(approx_eq(f[0], -100.0, 1e-8));
    }

    // ── Integration / energy tests ───────────────────────────────────────────

    #[test]
    fn spring_energy_proportional_to_stiffness() {
        // Two springs with different stiffness; force should scale proportionally
        let s1 = SpringConstraint::new([0.0; 3], 0.0, 100.0, 0.0);
        let s2 = SpringConstraint::new([0.0; 3], 0.0, 200.0, 0.0);
        let state = state_at([1.0, 0.0, 0.0]);
        let f1 = s1.compute_force(&state);
        let f2 = s2.compute_force(&state);
        assert!(approx_eq(f2[0] / f1[0], 2.0, 1e-10));
    }

    #[test]
    fn soft_joint_force_scales_with_stiffness() {
        let j1 = SoftJoint::new([0.0; 3], [1.0, 0.0, 0.0, 0.0], 100.0, 0.0);
        let j2 = SoftJoint::new([0.0; 3], [1.0, 0.0, 0.0, 0.0], 300.0, 0.0);
        let state = state_at([1.0, 0.0, 0.0]);
        let f1 = j1.compute_force(&state);
        let f2 = j2.compute_force(&state);
        assert!(approx_eq(f2[0] / f1[0], 3.0, 1e-10));
    }

    #[test]
    fn spring_force_3d() {
        let spring = SpringConstraint::new([0.0; 3], 0.0, 100.0, 0.0);
        let state = state_at([3.0, 4.0, 0.0]); // 5 m from origin
        let f = spring.compute_force(&state);
        let mag = vec3_norm(f);
        // F = 100 * 5 = 500 N toward origin
        assert!(approx_eq(mag, 500.0, 1e-6));
    }

    #[test]
    fn distance_constraint_zero_stiffness_no_force() {
        let a = RigidBodyState::with_position([0.0, 0.0, 0.0], 1.0);
        let b = RigidBodyState::with_position([5.0, 0.0, 0.0], 1.0);
        let dc = DistanceConstraint::new(a, b, 1.0, 0.0);
        let fa = dc.force_on_a();
        assert!(approx_eq3(fa, [0.0; 3], 1e-12));
    }

    #[test]
    fn penalty_method_symmetry() {
        let pm = PenaltyMethod::new(50.0, 0.0);
        // Equal and opposite violations should give equal and opposite forces
        let f_pos = pm.compute_scalar_force(1.0, 0.0);
        let f_neg = pm.compute_scalar_force(-1.0, 0.0);
        assert!(approx_eq(f_pos, -f_neg, 1e-10));
    }

    #[test]
    fn rigid_body_state_default_trait() {
        let s: RigidBodyState = Default::default();
        assert_eq!(s.position, [0.0; 3]);
        assert!(approx_eq(s.mass, 1.0, 1e-12));
    }

    #[test]
    fn constraint_solver_default_trait() {
        let solver: ConstraintSolver = Default::default();
        assert!(solver.is_empty());
    }
}
