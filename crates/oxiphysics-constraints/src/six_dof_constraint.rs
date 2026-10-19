// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Six-degree-of-freedom (6-DOF) constraint.
//!
//! Provides per-axis control over all 6 rigid-body DOFs:
//! - 3 linear (translation along X, Y, Z)
//! - 3 angular (rotation about X, Y, Z)
//!
//! Each axis can be `Free`, `Limited` (with min/max bounds), or `Locked`.
//! An optional motor can drive any axis toward a target velocity or position.

use oxiphysics_core::BodyHandle;
use oxiphysics_core::math::{Mat3, Quat, Real, Vec3};
use oxiphysics_rigid::RigidBodySet;

use crate::traits::Constraint;
use crate::{Error, Result};

// ─────────────────────────────────────────────────────────────────────────────
// MotorConfig
// ─────────────────────────────────────────────────────────────────────────────

/// Configuration for a motor on one DOF axis.
#[derive(Debug, Clone)]
pub struct MotorConfig {
    /// Target velocity for this axis (m/s or rad/s): the rate of change of the
    /// joint coordinate, i.e. the velocity of body B relative to body A along
    /// (about) the axis.  A positive target moves the coordinate towards its
    /// upper limit.
    pub target_velocity: f64,
    /// Maximum force/torque for this axis (N or N·m); the motor impulse of one
    /// step is limited to `max_force·dt`.  Must not be negative.
    pub max_force: f64,
    /// Simple proportional gain for velocity error.
    pub gain: f64,
}

impl MotorConfig {
    /// Create a new motor configuration.
    pub fn new(target_velocity: f64, max_force: f64, gain: f64) -> Self {
        Self {
            target_velocity,
            max_force,
            gain,
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// AxisConfig
// ─────────────────────────────────────────────────────────────────────────────

/// Per-axis configuration for a 6-DOF constraint.
#[derive(Debug, Clone)]
pub struct AxisConfig {
    /// If `true`, this axis is rigidly locked: its joint coordinate is held at
    /// zero.  A motor and a limit on a locked axis are ignored.
    pub locked: bool,
    /// Optional `(min, max)` limit on the joint coordinate of this axis
    /// (metres for a linear axis, radians for an angular one).
    /// Ignored when `locked == true`.
    ///
    /// The bounds are used in ascending order; a pair containing NaN is
    /// ignored.  Prefer [`SixDofConstraint::set_linear_limit`] /
    /// [`SixDofConstraint::set_angular_limit`], which reject such pairs.
    pub limited: Option<(f64, f64)>,
    /// Optional motor configuration.
    pub motor: Option<MotorConfig>,
}

impl AxisConfig {
    /// A completely free (unconstrained) axis.
    pub fn free() -> Self {
        Self {
            locked: false,
            limited: None,
            motor: None,
        }
    }

    /// A locked axis (zero relative motion).
    pub fn locked() -> Self {
        Self {
            locked: true,
            limited: None,
            motor: None,
        }
    }

    /// A limited axis with the given `[min, max]` bounds.
    pub fn limited(min: f64, max: f64) -> Self {
        Self {
            locked: false,
            limited: Some((min, max)),
            motor: None,
        }
    }

    /// Attach a motor to this axis.
    pub fn with_motor(mut self, motor: MotorConfig) -> Self {
        self.motor = Some(motor);
        self
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// SixDofConstraint
// ─────────────────────────────────────────────────────────────────────────────

/// A general 6-DOF constraint between two rigid bodies.
///
/// The constraint is expressed in the local frame of body A.  Per-axis
/// `AxisConfig` entries control which degrees of freedom are free, limited,
/// or locked.
///
/// # Axis convention
/// - `linear[0]` → translation along local X
/// - `linear[1]` → translation along local Y
/// - `linear[2]` → translation along local Z
/// - `angular[0]` → rotation about local X
/// - `angular[1]` → rotation about local Y
/// - `angular[2]` → rotation about local Z
///
/// # Joint coordinates
///
/// Every degree of freedom has a joint coordinate that describes body B
/// relative to body A, along (or about) the corresponding local axis of A:
///
/// - a linear coordinate is the displacement `(anchor_b − anchor_a)·axis`
///   (metres), zero when the two anchors coincide along the axis,
/// - an angular coordinate is the component, along the axis, of the rotation
///   vector of `rot_a⁻¹·rot_b·ref⁻¹` (radians), where `ref` is the *reference
///   rotation*: the relative orientation `rot_a⁻¹·rot_b` the joint was in when
///   it was first prepared (or the one set with
///   [`with_reference_rotation`](SixDofConstraint::with_reference_rotation)).
///   All angular coordinates are therefore zero in the initial configuration.
///
/// A locked axis keeps its coordinate at zero, a limited axis keeps it inside
/// `[min, max]` and a motor drives its rate of change (the relative velocity
/// of B with respect to A along the axis) to the motor's target velocity.
///
/// # Limits
///
/// A limit is a one-sided velocity constraint: the coordinate may approach a
/// bound by at most the remaining distance in one step (so it stops exactly at
/// the bound instead of overshooting), and beyond a bound it may not move
/// further out.  A violated limit is brought back by the position pass
/// ([`solve_position`](Constraint::solve_position)), which moves the bodies
/// without giving them a velocity.  A lock, in contrast, is a two-sided
/// constraint that also uses a Baumgarte velocity.
///
/// # Impulse sign convention
///
/// A positive row impulse pushes body B along `+axis` (or turns it about
/// `+axis`) and body A the opposite way, so it increases the coordinate.  The
/// accumulated impulse of a lower limit is therefore `≥ 0` and that of an upper
/// limit `≤ 0`.
#[derive(Debug, Clone)]
pub struct SixDofConstraint {
    /// Handle of body A.
    pub body_a: BodyHandle,
    /// Handle of body B.
    pub body_b: BodyHandle,
    /// Anchor point on body A in local-A space.
    pub anchor_a: Vec3,
    /// Anchor point on body B in local-B space.
    pub anchor_b: Vec3,
    /// Per-axis config for the three linear DOFs.
    pub linear: [AxisConfig; 3],
    /// Per-axis config for the three angular DOFs.
    pub angular: [AxisConfig; 3],

    // ── Cached solver data ────────────────────────────────────────────────────
    /// World-space offset from body A CoM to anchor (updated each `prepare`).
    r_a: Vec3,
    /// World-space offset from body B CoM to anchor (updated each `prepare`).
    r_b: Vec3,
    /// Relative orientation `rot_a⁻¹·rot_b` the angular coordinates are
    /// measured from (captured on the first `prepare`).
    reference_rotation: Option<Quat>,
    /// Solver rows of the six DOFs: three linear, then three angular.
    rows: [DofRow; 6],
    /// Joint coordinates measured by the latest `prepare`.
    coordinates: [Real; 6],
    /// Size of the velocity update of the latest velocity solve.
    last_residual: Real,
}

/// Baumgarte factor for positional correction.
const BAUMGARTE_6DOF: f64 = 0.2;
/// Largest linear correction (m) one position pass applies to a single row.
const MAX_LINEAR_CORRECTION: f64 = 0.2;
/// Largest angular correction (rad, about 8°) one position pass applies to a
/// single row.
const MAX_ANGULAR_CORRECTION: f64 = 0.14;
/// Errors below this size are not corrected by the position pass.
const MIN_POSITION_ERROR: f64 = 1e-12;
/// Smallest inverse effective mass a row needs to be solved.
const MIN_INVERSE_EFFECTIVE_MASS: f64 = 1e-12;
/// Smallest SOR factor accepted by [`Constraint::solve_velocity_relaxed`].
const MIN_SOR_OMEGA: f64 = 0.01;
/// Largest SOR factor accepted by [`Constraint::solve_velocity_relaxed`].
const MAX_SOR_OMEGA: f64 = 1.99;

/// One one-sided limit row of a DOF.
#[derive(Debug, Clone, Copy, Default)]
struct LimitRow {
    /// Whether the limit exists and is finite.
    enabled: bool,
    /// Velocity bound the row enforces: the coordinate rate must stay above
    /// (lower limit) or below (upper limit) this value.
    target: Real,
    /// Accumulated impulse of the row (`≥ 0` for a lower, `≤ 0` for an upper
    /// limit).
    lambda: Real,
}

/// Solver data of one degree of freedom.
#[derive(Debug, Clone, Copy)]
struct DofRow {
    /// World-space direction of the DOF (a local axis of body A).
    axis: Vec3,
    /// Effective mass `1 / K` of the DOF, `0` if no row is active or `K` vanishes.
    eff_mass: Real,
    /// Whether the DOF is locked.
    locked: bool,
    /// Rate the lock row drives the coordinate to (Baumgarte).
    lock_target: Real,
    /// Accumulated impulse of the lock row.
    lock_lambda: Real,
    /// Lower limit row.
    lower: LimitRow,
    /// Upper limit row.
    upper: LimitRow,
    /// Whether a motor drives the DOF.
    motor: bool,
    /// Rate the motor drives the coordinate to.
    motor_target: Real,
    /// Largest motor force or torque (non-negative).
    motor_max_force: Real,
    /// Accumulated impulse of the motor row.
    motor_lambda: Real,
}

impl DofRow {
    /// A DOF without any active row.
    fn idle() -> Self {
        Self {
            axis: Vec3::zeros(),
            eff_mass: 0.0,
            locked: false,
            lock_target: 0.0,
            lock_lambda: 0.0,
            lower: LimitRow::default(),
            upper: LimitRow::default(),
            motor: false,
            motor_target: 0.0,
            motor_max_force: 0.0,
            motor_lambda: 0.0,
        }
    }
}

/// State of one body that the constraint needs.
#[derive(Debug, Clone, Copy)]
struct Pose {
    position: Vec3,
    rotation: Quat,
    inv_mass: Real,
    inv_inertia: Mat3,
}

impl Pose {
    /// Read the pose of `handle`, `None` if the body does not exist.
    fn read(bodies: &RigidBodySet, handle: BodyHandle) -> Option<Self> {
        let body = bodies.get(handle)?;
        Some(Self {
            position: body.transform.position,
            rotation: body.transform.rotation,
            inv_mass: body.inverse_mass,
            inv_inertia: body.world_inverse_inertia,
        })
    }
}

/// The geometry shared by all rows of one velocity solve.
#[derive(Debug, Clone, Copy)]
struct Frame {
    body_a: BodyHandle,
    body_b: BodyHandle,
    r_a: Vec3,
    r_b: Vec3,
}

/// A rotation as the vector `angle·axis` along the shortest path.
fn rotation_vector(rotation: &Quat) -> Vec3 {
    let (w, v) = (rotation.w, Vec3::new(rotation.i, rotation.j, rotation.k));
    // q and -q are the same rotation: take the one with w >= 0 (angle <= π).
    let (w, v) = if w < 0.0 { (-w, -v) } else { (w, v) };
    let sin_half = v.norm();
    if sin_half < 1e-12 {
        return v * 2.0;
    }
    v * (2.0 * sin_half.atan2(w) / sin_half)
}

/// Project `x` onto `[lo, hi]`; never panics on inverted or NaN bounds.
fn project(x: Real, lo: Real, hi: Real) -> Real {
    x.max(lo).min(hi)
}

/// `(min, max)` of a limit pair in ascending order, `None` if a bound is NaN.
fn ordered_limits(limited: Option<(f64, f64)>) -> Option<(f64, f64)> {
    let (a, b) = limited?;
    if a.is_nan() || b.is_nan() {
        return None;
    }
    Some((a.min(b), a.max(b)))
}

/// The magnitude of a force limit (a NaN limit counts as 0), so the clamp of a
/// motor impulse always has `-limit <= limit`.
fn sanitized_force(max_force: f64) -> f64 {
    let magnitude = max_force.abs();
    if magnitude.is_nan() { 0.0 } else { magnitude }
}

/// Validate the bounds of a limit.
fn check_limits(min: f64, max: f64) -> Result<()> {
    if min.is_nan() || max.is_nan() || min > max {
        return Err(Error::joint_limits(min, max));
    }
    Ok(())
}

/// Validate a motor configuration.
fn check_motor(motor: &MotorConfig) -> Result<()> {
    if motor.target_velocity.is_nan() {
        return Err(Error::controller("motor target velocity is NaN"));
    }
    if motor.max_force.is_nan() || motor.max_force < 0.0 {
        return Err(Error::controller(format!(
            "motor max_force must be non-negative, got {}",
            motor.max_force
        )));
    }
    Ok(())
}

/// Validate an axis index.
fn check_axis(axis: usize) -> Result<()> {
    if axis >= 3 {
        return Err(Error::config("axis", axis as f64, 0.0, 2.0));
    }
    Ok(())
}

impl SixDofConstraint {
    /// Create a new 6-DOF constraint with all axes free.
    pub fn new(body_a: BodyHandle, body_b: BodyHandle) -> Self {
        Self {
            body_a,
            body_b,
            anchor_a: Vec3::zeros(),
            anchor_b: Vec3::zeros(),
            linear: std::array::from_fn(|_| AxisConfig::free()),
            angular: std::array::from_fn(|_| AxisConfig::free()),
            r_a: Vec3::zeros(),
            r_b: Vec3::zeros(),
            reference_rotation: None,
            rows: [DofRow::idle(); 6],
            coordinates: [0.0; 6],
            last_residual: 0.0,
        }
    }

    /// Set the anchor points in local body-space.
    pub fn with_anchors(mut self, anchor_a: Vec3, anchor_b: Vec3) -> Self {
        self.anchor_a = anchor_a;
        self.anchor_b = anchor_b;
        self
    }

    /// Fix the relative orientation `rot_a⁻¹·rot_b` that the angular
    /// coordinates (locks, limits, motors) are measured from.
    ///
    /// Without it the joint uses the relative orientation of its bodies at the
    /// first [`prepare`](Constraint::prepare), so the initial configuration has
    /// all angular coordinates at zero.  Pass the identity to measure angles
    /// from the configuration in which the two body frames are aligned.
    pub fn with_reference_rotation(mut self, reference: Quat) -> Self {
        self.reference_rotation = Some(reference);
        self
    }

    /// The reference rotation of the angular coordinates, `None` until it is
    /// set or captured by the first `prepare`.
    pub fn reference_rotation(&self) -> Option<Quat> {
        self.reference_rotation
    }

    /// The joint coordinates measured by the latest `prepare`: the three linear
    /// displacements in metres followed by the three angles in radians, each
    /// along (about) the local axis of body A.  All zero before the first
    /// `prepare` and while a body of the constraint is missing.
    pub fn coordinates(&self) -> [f64; 6] {
        self.coordinates
    }

    /// Lock the linear DOF along the given axis (0=X, 1=Y, 2=Z).
    ///
    /// Panics if `axis >= 3`.
    pub fn lock_linear(&mut self, axis: usize) {
        assert!(axis < 3, "axis must be 0, 1, or 2");
        self.linear[axis] = AxisConfig::locked();
    }

    /// Lock the angular DOF about the given axis (0=X, 1=Y, 2=Z).
    ///
    /// Panics if `axis >= 3`.
    pub fn lock_angular(&mut self, axis: usize) {
        assert!(axis < 3, "axis must be 0, 1, or 2");
        self.angular[axis] = AxisConfig::locked();
    }

    /// Set a linear limit `[min, max]` (metres of the joint coordinate) for the
    /// given axis.
    ///
    /// Fails with [`Error::InvalidJointLimits`] if `min > max` or a bound is
    /// NaN, and with [`Error::ConfigurationError`] if `axis >= 3`.
    pub fn set_linear_limit(&mut self, axis: usize, min: f64, max: f64) -> Result<()> {
        check_axis(axis)?;
        check_limits(min, max)?;
        self.linear[axis] = AxisConfig::limited(min, max);
        Ok(())
    }

    /// Set an angular limit `[min, max]` (radians of the joint coordinate) for
    /// the given axis.
    ///
    /// Fails with [`Error::InvalidJointLimits`] if `min > max` or a bound is
    /// NaN, and with [`Error::ConfigurationError`] if `axis >= 3`.
    pub fn set_angular_limit(&mut self, axis: usize, min: f64, max: f64) -> Result<()> {
        check_axis(axis)?;
        check_limits(min, max)?;
        self.angular[axis] = AxisConfig::limited(min, max);
        Ok(())
    }

    /// Attach a motor to a linear axis.
    ///
    /// Fails with [`Error::ControllerError`] for a NaN target velocity or a
    /// negative or NaN `max_force`, and with [`Error::ConfigurationError`] if
    /// `axis >= 3`.
    pub fn set_linear_motor(&mut self, axis: usize, motor: MotorConfig) -> Result<()> {
        check_axis(axis)?;
        check_motor(&motor)?;
        self.linear[axis].motor = Some(motor);
        Ok(())
    }

    /// Attach a motor to an angular axis.
    ///
    /// Fails with [`Error::ControllerError`] for a NaN target velocity or a
    /// negative or NaN `max_force`, and with [`Error::ConfigurationError`] if
    /// `axis >= 3`.
    pub fn set_angular_motor(&mut self, axis: usize, motor: MotorConfig) -> Result<()> {
        check_axis(axis)?;
        check_motor(&motor)?;
        self.angular[axis].motor = Some(motor);
        Ok(())
    }

    /// Create a fully-locked constraint (weld joint equivalent).
    pub fn weld(body_a: BodyHandle, body_b: BodyHandle) -> Self {
        let mut c = Self::new(body_a, body_b);
        for i in 0..3 {
            c.lock_linear(i);
            c.lock_angular(i);
        }
        c
    }

    // ── Jacobian helpers ──────────────────────────────────────────────────────

    /// Inverse effective mass of a linear row along `axis`:
    ///   k = inv_m_a + inv_m_b + (r_a×d)·I_a^{-1}·(r_a×d) + (r_b×d)·I_b^{-1}·(r_b×d)
    fn linear_inverse_mass(a: &Pose, b: &Pose, r_a: &Vec3, r_b: &Vec3, axis: &Vec3) -> f64 {
        let rn_a = r_a.cross(axis);
        let rn_b = r_b.cross(axis);
        a.inv_mass
            + b.inv_mass
            + rn_a.dot(&(a.inv_inertia * rn_a))
            + rn_b.dot(&(b.inv_inertia * rn_b))
    }

    /// Inverse effective mass of an angular row about `axis`:
    ///   k = d·I_a^{-1}·d + d·I_b^{-1}·d
    fn angular_inverse_mass(a: &Pose, b: &Pose, axis: &Vec3) -> f64 {
        axis.dot(&(a.inv_inertia * axis)) + axis.dot(&(b.inv_inertia * axis))
    }

    /// The configuration of DOF `dof` (0..3 linear, 3..6 angular).
    fn config(&self, dof: usize) -> &AxisConfig {
        if dof < 3 {
            &self.linear[dof]
        } else {
            &self.angular[dof - 3]
        }
    }

    /// The three world-space axes (local axes of body A) and the six joint
    /// coordinates for the given poses.  Captures the reference rotation if it
    /// is not known yet.
    fn measure(&mut self, a: &Pose, b: &Pose) -> ([Vec3; 3], [Real; 6]) {
        let axes = [
            a.rotation * Vec3::new(1.0, 0.0, 0.0),
            a.rotation * Vec3::new(0.0, 1.0, 0.0),
            a.rotation * Vec3::new(0.0, 0.0, 1.0),
        ];
        let r_a = a.rotation * self.anchor_a;
        let r_b = b.rotation * self.anchor_b;
        let offset = (b.position + r_b) - (a.position + r_a);

        // Rotation of B relative to A, in A's frame, beyond the reference.
        let relative = a.rotation.inverse() * b.rotation;
        let reference = *self.reference_rotation.get_or_insert(relative);
        let angles = rotation_vector(&(relative * reference.inverse()));

        let coordinates = [
            offset.dot(&axes[0]),
            offset.dot(&axes[1]),
            offset.dot(&axes[2]),
            angles.x,
            angles.y,
            angles.z,
        ];
        (axes, coordinates)
    }

    /// Residual position error of DOF `dof` for its current coordinate: the
    /// coordinate itself for a locked DOF, the distance to the violated limit
    /// (zero inside the limits) for a limited one.
    fn position_error(&self, dof: usize, coordinate: Real) -> Real {
        let config = self.config(dof);
        if config.locked {
            return coordinate;
        }
        let Some((min, max)) = ordered_limits(config.limited) else {
            return 0.0;
        };
        if coordinate < min {
            coordinate - min
        } else if coordinate > max {
            coordinate - max
        } else {
            0.0
        }
    }

    /// Translate and rotate `handle` by a position-level impulse `linear` at
    /// lever arm `r` (zero for a pure torque) and an angular impulse `angular`.
    fn move_body(
        bodies: &mut RigidBodySet,
        handle: BodyHandle,
        linear: &Vec3,
        r: &Vec3,
        angular: &Vec3,
        sign: Real,
    ) {
        let Some(body) = bodies.get_mut(handle) else {
            return;
        };
        body.transform.position += linear * (sign * body.inverse_mass);
        let turn = body.world_inverse_inertia * (r.cross(linear) + angular) * sign;
        if turn.norm_squared() > 0.0 {
            body.transform.rotation = Quat::from_scaled_axis(turn) * body.transform.rotation;
            body.update_world_inertia();
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Row solving
// ─────────────────────────────────────────────────────────────────────────────

/// Rate of change of the coordinate of DOF `dof`: the velocity of B relative to
/// A along (about) `axis`.
fn coordinate_rate(bodies: &RigidBodySet, frame: &Frame, dof: usize, axis: &Vec3) -> Option<Real> {
    let a = bodies.get(frame.body_a)?;
    let b = bodies.get(frame.body_b)?;
    Some(if dof < 3 {
        let point_a = a.velocity + a.angular_velocity.cross(&frame.r_a);
        let point_b = b.velocity + b.angular_velocity.cross(&frame.r_b);
        (point_b - point_a).dot(axis)
    } else {
        (b.angular_velocity - a.angular_velocity).dot(axis)
    })
}

/// Apply the row impulse `delta` of DOF `dof`: `+delta` on body B and `-delta`
/// on body A, along (about) `axis`.
fn apply_row_impulse(
    bodies: &mut RigidBodySet,
    frame: &Frame,
    dof: usize,
    axis: &Vec3,
    delta: Real,
) {
    let impulse = axis * delta;
    if let Some(b) = bodies.get_mut(frame.body_b) {
        if dof < 3 {
            b.velocity += impulse * b.inverse_mass;
            b.angular_velocity += b.world_inverse_inertia * frame.r_b.cross(&impulse);
        } else {
            b.angular_velocity += b.world_inverse_inertia * impulse;
        }
    }
    if let Some(a) = bodies.get_mut(frame.body_a) {
        if dof < 3 {
            a.velocity -= impulse * a.inverse_mass;
            a.angular_velocity -= a.world_inverse_inertia * frame.r_a.cross(&impulse);
        } else {
            a.angular_velocity -= a.world_inverse_inertia * impulse;
        }
    }
}

/// What a row does to the coordinate rate of its DOF.
#[derive(Debug, Clone, Copy)]
struct RowGoal {
    /// Index of the DOF (0..3 linear, 3..6 angular).
    dof: usize,
    /// World-space direction of the DOF.
    axis: Vec3,
    /// Effective mass of the DOF.
    eff_mass: Real,
    /// Coordinate rate the row drives towards.
    target: Real,
    /// Bounds of the accumulated impulse.
    bounds: (Real, Real),
}

/// One Gauss-Seidel update of a row: `λ ← project(λ + ω·m·(target − rate))`,
/// applied to the bodies.  Returns the size of the velocity change the update
/// made to the coordinate rate (`|Δλ|·K`).
fn solve_row(
    bodies: &mut RigidBodySet,
    frame: &Frame,
    goal: &RowGoal,
    lambda: &mut Real,
    omega: Real,
) -> Real {
    let Some(rate) = coordinate_rate(bodies, frame, goal.dof, &goal.axis) else {
        return 0.0;
    };
    let old = *lambda;
    let new = project(
        old + omega * goal.eff_mass * (goal.target - rate),
        goal.bounds.0,
        goal.bounds.1,
    );
    let applied = new - old;
    *lambda = new;
    if applied != 0.0 {
        apply_row_impulse(bodies, frame, goal.dof, &goal.axis, applied);
    }
    (applied / goal.eff_mass).abs()
}

impl SixDofConstraint {
    /// One sweep over all rows with the SOR factor `omega`; stores and returns
    /// the size of the velocity update.
    fn solve_rows(&mut self, bodies: &mut RigidBodySet, dt: Real, omega: Real) -> Real {
        let frame = Frame {
            body_a: self.body_a,
            body_b: self.body_b,
            r_a: self.r_a,
            r_b: self.r_b,
        };
        let step = if dt > 0.0 && dt.is_finite() { dt } else { 0.0 };
        let mut sum_sq = 0.0;
        for dof in 0..6 {
            let row = self.rows[dof];
            if row.eff_mass <= 0.0 {
                continue;
            }
            let goal = |target: Real, bounds: (Real, Real)| RowGoal {
                dof,
                axis: row.axis,
                eff_mass: row.eff_mass,
                target,
                bounds,
            };
            if row.locked {
                let lock = goal(row.lock_target, (f64::NEG_INFINITY, f64::INFINITY));
                let lambda = &mut self.rows[dof].lock_lambda;
                sum_sq += solve_row(bodies, &frame, &lock, lambda, omega).powi(2);
                continue;
            }
            if row.motor {
                let limit = row.motor_max_force * step;
                let motor = goal(row.motor_target, (-limit, limit));
                let lambda = &mut self.rows[dof].motor_lambda;
                sum_sq += solve_row(bodies, &frame, &motor, lambda, omega).powi(2);
            }
            if row.lower.enabled {
                let lower = goal(row.lower.target, (0.0, f64::INFINITY));
                let lambda = &mut self.rows[dof].lower.lambda;
                sum_sq += solve_row(bodies, &frame, &lower, lambda, omega).powi(2);
            }
            if row.upper.enabled {
                let upper = goal(row.upper.target, (f64::NEG_INFINITY, 0.0));
                let lambda = &mut self.rows[dof].upper.lambda;
                sum_sq += solve_row(bodies, &frame, &upper, lambda, omega).powi(2);
            }
        }
        self.last_residual = sum_sq.sqrt();
        self.last_residual
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Constraint impl
// ─────────────────────────────────────────────────────────────────────────────

impl Constraint for SixDofConstraint {
    fn prepare(&mut self, bodies: &RigidBodySet, dt: f64) {
        self.last_residual = 0.0;
        self.rows = [DofRow::idle(); 6];
        self.coordinates = [0.0; 6];
        let (Some(a), Some(b)) = (
            Pose::read(bodies, self.body_a),
            Pose::read(bodies, self.body_b),
        ) else {
            return;
        };

        // Anchor offsets in world space
        self.r_a = a.rotation * self.anchor_a;
        self.r_b = b.rotation * self.anchor_b;
        let (axes, coordinates) = self.measure(&a, &b);
        self.coordinates = coordinates;

        // A paused (dt = 0) or invalid step has no rate to correct at: locks
        // stay at zero rate, limits are left alone.
        let inv_dt = if dt > 1e-12 && dt.is_finite() {
            1.0 / dt
        } else {
            0.0
        };

        for dof in 0..6 {
            let config = self.config(dof).clone();
            let limits = if config.locked {
                None
            } else {
                ordered_limits(config.limited)
            };
            let motor = if config.locked {
                None
            } else {
                config.motor.as_ref()
            };
            if !config.locked && limits.is_none() && motor.is_none() {
                continue;
            }

            let (axis, coordinate) = (axes[dof % 3], coordinates[dof]);
            let inverse_mass = if dof < 3 {
                Self::linear_inverse_mass(&a, &b, &self.r_a, &self.r_b, &axis)
            } else {
                Self::angular_inverse_mass(&a, &b, &axis)
            };
            let mut row = DofRow::idle();
            row.axis = axis;
            row.eff_mass = if inverse_mass > MIN_INVERSE_EFFECTIVE_MASS {
                1.0 / inverse_mass
            } else {
                0.0
            };

            if config.locked {
                row.locked = true;
                // Drive the coordinate to zero with the Baumgarte rate.
                row.lock_target = -BAUMGARTE_6DOF * coordinate * inv_dt;
            }
            if let Some((min, max)) = limits.filter(|_| inv_dt > 0.0) {
                // Speculative one-sided rows: the coordinate may approach a
                // limit by at most the remaining gap in one step, and once it
                // is beyond a limit it may not move further out.  The violation
                // itself is removed by the position pass; pushing it back with
                // a Baumgarte velocity would leave that velocity in the bodies
                // and carry them across the allowed range.
                let lower_gap = coordinate - min;
                row.lower = LimitRow {
                    enabled: min.is_finite(),
                    target: if lower_gap >= 0.0 {
                        -lower_gap * inv_dt
                    } else {
                        0.0
                    },
                    lambda: 0.0,
                };
                let upper_gap = max - coordinate;
                row.upper = LimitRow {
                    enabled: max.is_finite(),
                    target: if upper_gap >= 0.0 {
                        upper_gap * inv_dt
                    } else {
                        0.0
                    },
                    lambda: 0.0,
                };
            }
            if let Some(m) = motor {
                row.motor = true;
                row.motor_target = if m.target_velocity.is_nan() {
                    0.0
                } else {
                    m.target_velocity
                };
                row.motor_max_force = sanitized_force(m.max_force);
            }
            self.rows[dof] = row;
        }
    }

    fn solve_velocity(&mut self, bodies: &mut RigidBodySet, dt: f64) {
        self.solve_rows(bodies, dt, 1.0);
    }

    fn solve_velocity_relaxed(&mut self, bodies: &mut RigidBodySet, dt: f64, omega: f64) -> f64 {
        let omega = if omega.is_nan() {
            1.0
        } else {
            omega.clamp(MIN_SOR_OMEGA, MAX_SOR_OMEGA)
        };
        self.solve_rows(bodies, dt, omega)
    }

    fn residual(&self) -> f64 {
        self.last_residual
    }

    fn solve_position(&mut self, bodies: &mut RigidBodySet, _dt: f64) {
        for dof in 0..6 {
            let config = self.config(dof);
            let corrects = config.locked || ordered_limits(config.limited).is_some();
            if !corrects {
                continue;
            }
            // Re-measure for every row: the previous rows moved the bodies.
            let (Some(a), Some(b)) = (
                Pose::read(bodies, self.body_a),
                Pose::read(bodies, self.body_b),
            ) else {
                return;
            };
            let (axes, coordinates) = self.measure(&a, &b);
            let error = self.position_error(dof, coordinates[dof]);
            if error.abs() < MIN_POSITION_ERROR {
                continue;
            }
            let axis = axes[dof % 3];
            let max_correction = if dof < 3 {
                MAX_LINEAR_CORRECTION
            } else {
                MAX_ANGULAR_CORRECTION
            };
            // Change of the coordinate the pass asks for (opposite to the error).
            let wanted = project(-BAUMGARTE_6DOF * error, -max_correction, max_correction);

            if dof < 3 {
                let r_a = a.rotation * self.anchor_a;
                let r_b = b.rotation * self.anchor_b;
                let inverse_mass = Self::linear_inverse_mass(&a, &b, &r_a, &r_b, &axis);
                if inverse_mass < MIN_INVERSE_EFFECTIVE_MASS {
                    continue;
                }
                // A position impulse along +axis moves B along +axis and A the
                // other way, which increases the coordinate.
                let impulse = axis * (wanted / inverse_mass);
                Self::move_body(bodies, self.body_b, &impulse, &r_b, &Vec3::zeros(), 1.0);
                Self::move_body(bodies, self.body_a, &impulse, &r_a, &Vec3::zeros(), -1.0);
            } else {
                let inverse_mass = Self::angular_inverse_mass(&a, &b, &axis);
                if inverse_mass < MIN_INVERSE_EFFECTIVE_MASS {
                    continue;
                }
                let torque = axis * (wanted / inverse_mass);
                let nowhere = Vec3::zeros();
                Self::move_body(bodies, self.body_b, &nowhere, &nowhere, &torque, 1.0);
                Self::move_body(bodies, self.body_a, &nowhere, &nowhere, &torque, -1.0);
            }
        }
    }

    fn body_handles(&self) -> Vec<BodyHandle> {
        vec![self.body_a, self.body_b]
    }

    /// The number of locked or limited axes (motors do not constrain).
    fn dof_count(&self) -> usize {
        (0..6)
            .filter(|&dof| {
                let config = self.config(dof);
                config.locked || ordered_limits(config.limited).is_some()
            })
            .count()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Jacobian matrix helper (6×12 matrix, rows = constraints, cols = velocities)
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the 6×12 Jacobian for a fully-locked 6-DOF constraint.
///
/// The Jacobian maps `[v_a, ω_a, v_b, ω_b]` (12 entries) to the 6 constraint
/// velocities (3 linear + 3 angular).
///
/// Row layout:
/// - Rows 0-2: linear constraints along X, Y, Z
/// - Rows 3-5: angular constraints about X, Y, Z
///
/// Column layout per body: `[v.x, v.y, v.z, ω.x, ω.y, ω.z]`
pub fn compute_jacobian_6dof(r_a: Vec3, r_b: Vec3, axes: &[Vec3; 3]) -> [[f64; 12]; 6] {
    let mut j = [[0.0f64; 12]; 6];

    // Linear rows
    for (row, ax) in axes.iter().enumerate() {
        // Body A linear part: +ax
        j[row][0] = ax.x;
        j[row][1] = ax.y;
        j[row][2] = ax.z;
        // Body A angular part: +(r_a × ax)
        let ra_x_ax = r_a.cross(ax);
        j[row][3] = ra_x_ax.x;
        j[row][4] = ra_x_ax.y;
        j[row][5] = ra_x_ax.z;
        // Body B linear part: -ax
        j[row][6] = -ax.x;
        j[row][7] = -ax.y;
        j[row][8] = -ax.z;
        // Body B angular part: -(r_b × ax)
        let rb_x_ax = r_b.cross(ax);
        j[row][9] = -rb_x_ax.x;
        j[row][10] = -rb_x_ax.y;
        j[row][11] = -rb_x_ax.z;
    }

    // Angular rows
    for (i, ax) in axes.iter().enumerate() {
        let row = 3 + i;
        // Body A angular: +ax
        j[row][3] = ax.x;
        j[row][4] = ax.y;
        j[row][5] = ax.z;
        // Body B angular: -ax
        j[row][9] = -ax.x;
        j[row][10] = -ax.y;
        j[row][11] = -ax.z;
    }

    j
}

// ─────────────────────────────────────────────────────────────────────────────
// SpringDamperConfig
// ─────────────────────────────────────────────────────────────────────────────

/// Spring-damper configuration for a single DOF axis.
///
/// The force law is:
///   `F = -k * (x - x_rest) - c * v`
///
/// where `k` is the stiffness, `c` is the damping, `x` is the current
/// displacement, `x_rest` is the natural length and `v` is the velocity.
#[derive(Debug, Clone)]
pub struct SpringDamperConfig {
    /// Spring stiffness \[N/m or N·m/rad\].
    pub stiffness: f64,
    /// Damping coefficient \[N·s/m or N·m·s/rad\].
    pub damping: f64,
    /// Natural (rest) length / angle \[m or rad\].
    pub rest_length: f64,
    /// Pre-load force/torque \[N or N·m\].
    pub preload: f64,
}

impl SpringDamperConfig {
    /// Create a new spring-damper configuration.
    pub fn new(stiffness: f64, damping: f64, rest_length: f64, preload: f64) -> Self {
        Self {
            stiffness,
            damping,
            rest_length,
            preload,
        }
    }

    /// Compute the force/torque: `F = -k*(x - x_rest) - c*v + F_preload`.
    pub fn force(&self, position: f64, velocity: f64) -> f64 {
        -self.stiffness * (position - self.rest_length) - self.damping * velocity + self.preload
    }

    /// Compute the damping ratio: `zeta = c / (2 * sqrt(k * m))`.
    /// Requires the effective mass `m`.
    pub fn damping_ratio(&self, effective_mass: f64) -> f64 {
        if self.stiffness < 1e-30 || effective_mass < 1e-30 {
            return 0.0;
        }
        self.damping / (2.0 * (self.stiffness * effective_mass).sqrt())
    }

    /// Natural frequency: `omega_n = sqrt(k / m)`.
    pub fn natural_frequency(&self, effective_mass: f64) -> f64 {
        if effective_mass < 1e-30 {
            return 0.0;
        }
        (self.stiffness / effective_mass).sqrt()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// VelocityDriveConfig
// ─────────────────────────────────────────────────────────────────────────────

/// Velocity drive configuration for a single DOF axis.
///
/// A PD-style velocity controller that computes the constraint impulse
/// needed to achieve the target velocity.
#[derive(Debug, Clone)]
pub struct VelocityDriveConfig {
    /// Target velocity \[m/s or rad/s\].
    pub target_velocity: f64,
    /// Maximum force/torque magnitude \[N or N·m\].
    pub max_force: f64,
    /// Proportional gain for velocity error.
    pub kp: f64,
}

impl VelocityDriveConfig {
    /// Create a simple linear velocity drive.
    pub fn new_linear(target_velocity: f64, max_force: f64) -> Self {
        Self {
            target_velocity,
            max_force,
            kp: 1.0,
        }
    }

    /// Create a velocity drive with explicit gain.
    pub fn new(target_velocity: f64, max_force: f64, kp: f64) -> Self {
        Self {
            target_velocity,
            max_force,
            kp,
        }
    }

    /// Compute the constraint impulse: `imp = kp * (v_target - v_actual) * dt`.
    ///
    /// Clamped to `[-max_force * dt, max_force * dt]`.
    pub fn compute_impulse(&self, effective_mass: f64, actual_velocity: f64, dt: f64) -> f64 {
        let vel_error = self.target_velocity - actual_velocity;
        let imp = effective_mass * self.kp * vel_error;
        let limit = self.max_force * dt;
        imp.clamp(-limit, limit)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PositionDriveConfig
// ─────────────────────────────────────────────────────────────────────────────

/// Position drive configuration for a single DOF axis.
///
/// A PD (proportional-derivative) controller that drives the position to
/// a target while damping velocity.
#[derive(Debug, Clone)]
pub struct PositionDriveConfig {
    /// Target position \[m or rad\].
    pub target_position: f64,
    /// Proportional gain (stiffness-like) \[N/m or N·m/rad\].
    pub kp: f64,
    /// Derivative gain (damping-like) \[N·s/m or N·m·s/rad\].
    pub kd: f64,
    /// Maximum force/torque magnitude \[N or N·m\].
    pub max_force: f64,
}

impl PositionDriveConfig {
    /// Create a new position drive.
    pub fn new(target_position: f64, kp: f64, kd: f64, max_force: f64) -> Self {
        Self {
            target_position,
            kp,
            kd,
            max_force,
        }
    }

    /// Compute the drive force: `F = kp * (x_target - x) - kd * v`.
    ///
    /// Clamped to `[-max_force, max_force]`.
    pub fn compute_force(&self, position: f64, velocity: f64) -> f64 {
        let f = self.kp * (self.target_position - position) - self.kd * velocity;
        f.clamp(-self.max_force, self.max_force)
    }

    /// Error from target position.
    pub fn position_error(&self, position: f64) -> f64 {
        self.target_position - position
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Compliance and restitution helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Apply constraint compliance (softness) to the effective mass.
///
/// CFM (Constraint Force Mixing): modifies the diagonal of the system matrix
/// to allow slight violation.  The effective mass becomes:
///
///   `m_eff_c = m_eff / (1 + cfm * m_eff / dt)`
///
/// where `cfm` is the compliance parameter.
pub fn apply_compliance(effective_mass: f64, cfm: f64, dt: f64) -> f64 {
    let denom = 1.0 + cfm * effective_mass / dt.max(1e-30);
    effective_mass / denom
}

/// Compute the restitution impulse for a collision along a constraint axis.
///
/// `J_restitution = -(1 + e) * v_approach * effective_mass`
///
/// where `e` is the coefficient of restitution \[0, 1\] and `v_approach < 0`.
pub fn restitution_impulse(effective_mass: f64, v_approach: f64, restitution: f64) -> f64 {
    if v_approach >= 0.0 {
        return 0.0;
    }
    -(1.0 + restitution) * v_approach * effective_mass
}

/// Compute the post-collision velocity of body A using the impulse-momentum law.
///
/// `v_a_new = v_a + J * inv_mass_a`
pub fn post_collision_velocity(v_a: f64, impulse: f64, inv_mass_a: f64) -> f64 {
    v_a + impulse * inv_mass_a
}

/// Coefficient of restitution from pre/post velocities.
///
/// `e = -(v_a_post - v_b_post) / (v_a_pre - v_b_pre)`
pub fn coefficient_of_restitution(v_a_pre: f64, v_b_pre: f64, v_a_post: f64, v_b_post: f64) -> f64 {
    let dv_pre = v_a_pre - v_b_pre;
    if dv_pre.abs() < 1e-30 {
        return 0.0;
    }
    -(v_a_post - v_b_post) / dv_pre
}

// ─────────────────────────────────────────────────────────────────────────────
// 6-DOF limit utilities
// ─────────────────────────────────────────────────────────────────────────────

/// Check whether a DOF value is within its limit.
///
/// Returns `true` if `value` is within `[min, max]`.
pub fn within_limit(value: f64, min: f64, max: f64) -> bool {
    value >= min && value <= max
}

/// Compute the limit violation (penetration) for a constrained DOF.
///
/// Returns `> 0` if above the max, `< 0` if below the min, `0` if within.
pub fn limit_violation(value: f64, min: f64, max: f64) -> f64 {
    if value < min {
        value - min
    } else if value > max {
        value - max
    } else {
        0.0
    }
}

/// Clamp an accumulated impulse to its limit range and return the applied delta.
pub fn clamp_accumulated_impulse(old: f64, delta: f64, min: f64, max: f64) -> (f64, f64) {
    let new = (old + delta).clamp(min, max);
    let applied = new - old;
    (new, applied)
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests;
