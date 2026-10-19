// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Shared constraint math for the joint solvers.
//!
//! The joints in this module solve their rows exactly instead of with a scalar
//! approximation.  This file collects the building blocks they share:
//!
//! * the 3x3 effective-mass matrix of a point-to-point constraint
//!   ([`point_constraint_k`]) and the exact solves built on it,
//! * projected solves of one or two rows inside a subspace
//!   ([`constrained_impulse`]) for hinges, sliders and universal joints,
//! * the inertia-weighted orientation correction used by the position pass
//!   ([`apply_orientation_correction`]),
//! * the hinge angle and axis misalignment of a revolute joint,
//! * the soft (implicit spring-damper) velocity row ([`SoftRow`]) shared by the
//!   spring, breakable and XPBD joints.
//!
//! Conventions: an impulse `P` handed to [`apply_pair_impulse`] adds `+P` to
//! body A and `-P` to body B, so the relative anchor velocity `v_A - v_B`
//! changes by `K·P`.  An impulse handed to [`apply_angular_impulse`] changes the
//! relative angular velocity `ω_A - ω_B` by `(I_A⁻¹ + I_B⁻¹)·P`.

use std::f64::consts::{PI, TAU};

use oxiphysics_core::BodyHandle;
use oxiphysics_core::math::{Mat3, Quat, Real, Vec3};
use oxiphysics_rigid::RigidBodySet;

use super::functions::{apply_angular_impulse, apply_pair_impulse, read_body, small_rotation_quat};
use super::types::BodyProps;

/// Relative size below which a pivot or eigenvalue is treated as zero.
const PIVOT_TOLERANCE: Real = 1e-12;

/// Skew-symmetric (cross product) matrix `[v]×`, so that `[v]× w = v × w`.
pub(super) fn skew(v: &Vec3) -> Mat3 {
    Mat3::new(0.0, -v.z, v.y, v.z, 0.0, -v.x, -v.y, v.x, 0.0)
}

/// Effective-mass matrix `K` of a point constraint with lever arms `r_a`, `r_b`.
///
/// `K = (1/m_A + 1/m_B)·I − [r_A]×·I_A⁻¹·[r_A]× − [r_B]×·I_B⁻¹·[r_B]×`
///
/// Applying the impulse `P` (see the module conventions) changes the relative
/// anchor velocity by `K·P`; the exact impulse that cancels a relative velocity
/// `dv` is therefore `-K⁻¹·dv`.
pub(super) fn point_constraint_k(a: &BodyProps, b: &BodyProps, r_a: &Vec3, r_b: &Vec3) -> Mat3 {
    let skew_a = skew(r_a);
    let skew_b = skew(r_b);
    Mat3::identity() * (a.inv_mass + b.inv_mass)
        - skew_a * a.inv_inertia * skew_a
        - skew_b * b.inv_inertia * skew_b
}

/// Inverse of a symmetric positive semi-definite matrix.
///
/// Well-conditioned matrices are inverted directly.  A singular or nearly
/// singular matrix (for example two bodies that cannot translate and an anchor
/// that lies on the rotation axis) gets its Moore-Penrose pseudo-inverse, so
/// the result never contains the huge entries a plain inverse of a numerically
/// singular matrix would; the all-zero matrix maps to the all-zero matrix.
pub(super) fn stable_inverse(k: &Mat3) -> Mat3 {
    let scale = k.norm();
    if !scale.is_finite() || scale <= 0.0 {
        return Mat3::zeros();
    }
    if k.determinant().abs() > PIVOT_TOLERANCE * scale * scale * scale
        && let Some(inverse) = k.try_inverse()
    {
        return inverse;
    }
    k.pseudo_inverse(PIVOT_TOLERANCE * scale)
        .unwrap_or_else(|_| Mat3::zeros())
}

/// Inverse effective mass `K = J·M⁻¹·Jᵀ` of the linear row along unit `dir`
/// between the anchors `r_a` and `r_b`.
pub(super) fn inverse_effective_mass(
    a: &BodyProps,
    b: &BodyProps,
    r_a: &Vec3,
    r_b: &Vec3,
    dir: &Vec3,
) -> Real {
    let torque_arm_a = r_a.cross(dir);
    let torque_arm_b = r_b.cross(dir);
    a.inv_mass
        + b.inv_mass
        + torque_arm_a.dot(&(a.inv_inertia * torque_arm_a))
        + torque_arm_b.dot(&(b.inv_inertia * torque_arm_b))
}

/// Two unit vectors perpendicular to the unit vector `axis` and to each other.
pub(super) fn perpendicular_basis(axis: &Vec3) -> (Vec3, Vec3) {
    let reference = if axis.x.abs() < 0.9 {
        Vec3::new(1.0, 0.0, 0.0)
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    let first = axis.cross(&reference).normalize();
    let second = axis.cross(&first);
    (first, second)
}

/// Solve the symmetric 2x2 system `[[a, b], [b, d]]·x = rhs` with a
/// pseudo-inverse: directions whose eigenvalue is negligible are left alone.
fn solve_symmetric_2x2(a: Real, b: Real, d: Real, rhs: [Real; 2]) -> [Real; 2] {
    let mean = 0.5 * (a + d);
    let spread = (0.25 * (a - d) * (a - d) + b * b).sqrt();
    let largest = mean + spread;
    if !largest.is_finite() || largest <= 0.0 {
        return [0.0, 0.0];
    }
    let smallest = mean - spread;
    // Unit eigenvector of the larger eigenvalue.
    let (ex, ey) = if b.abs() > 1e-300 {
        let from_first_row = (largest - d, b);
        let from_second_row = (b, largest - a);
        let first_norm = from_first_row.0.hypot(from_first_row.1);
        let second_norm = from_second_row.0.hypot(from_second_row.1);
        let (x, y, norm) = if first_norm >= second_norm {
            (from_first_row.0, from_first_row.1, first_norm)
        } else {
            (from_second_row.0, from_second_row.1, second_norm)
        };
        (x / norm, y / norm)
    } else if a >= d {
        (1.0, 0.0)
    } else {
        (0.0, 1.0)
    };
    let cutoff = PIVOT_TOLERANCE * largest;
    let mut solution = [0.0, 0.0];
    let along_largest = (ex * rhs[0] + ey * rhs[1]) / largest;
    solution[0] += along_largest * ex;
    solution[1] += along_largest * ey;
    if smallest > cutoff {
        let (fx, fy) = (-ey, ex);
        let along_smallest = (fx * rhs[0] + fy * rhs[1]) / smallest;
        solution[0] += along_smallest * fx;
        solution[1] += along_smallest * fy;
    }
    solution
}

/// Impulse inside the span of `rows` that drives the rows of `rate` to zero.
///
/// `rate` is the relative velocity (anchor velocity or angular velocity) and
/// `k` the matching effective-mass matrix, so an impulse `P` changes the rate
/// by `K·P`.  The returned `P = Σ λ_i·rows[i]` solves
/// `rows[i]·(rate + K·P) = 0` exactly, which is the textbook result for a
/// constraint that only acts along its own rows (the hinge or slider force is
/// perpendicular to the free axis, so it does no work along it).
///
/// One and two rows are solved in the projected `1x1` / `2x2` system; three or
/// more rows are assumed to span the whole space and use the full `K⁻¹`.
pub(super) fn constrained_impulse(k: &Mat3, rows: &[Vec3], rate: &Vec3) -> Vec3 {
    // A degenerate row (for example the normalised zero vector of a joint with
    // a zero axis) must not turn into a NaN impulse.
    let finite = |v: &Vec3| v.iter().all(|c| c.is_finite());
    if !rows.iter().all(finite) || !finite(rate) || !k.iter().all(|c| c.is_finite()) {
        return Vec3::zeros();
    }
    match rows {
        [] => Vec3::zeros(),
        [row] => {
            let stiffness = row.dot(&(k * row));
            if stiffness > PIVOT_TOLERANCE * k.norm() {
                row * (-row.dot(rate) / stiffness)
            } else {
                Vec3::zeros()
            }
        }
        [first, second] => {
            let k_first = k * first;
            let k_second = k * second;
            let multipliers = solve_symmetric_2x2(
                first.dot(&k_first),
                first.dot(&k_second),
                second.dot(&k_second),
                [-first.dot(rate), -second.dot(rate)],
            );
            first * multipliers[0] + second * multipliers[1]
        }
        _ => -(stable_inverse(k) * rate),
    }
}

/// Cancel the relative velocity of the two anchors in all three directions.
///
/// Returns the size of the velocity change that was applied to the anchors.
pub(super) fn solve_point_velocity(
    bodies: &mut RigidBodySet,
    handle_a: BodyHandle,
    handle_b: BodyHandle,
    r_a: &Vec3,
    r_b: &Vec3,
) -> Real {
    let (Some(a), Some(b)) = (read_body(bodies, handle_a), read_body(bodies, handle_b)) else {
        return 0.0;
    };
    let relative_velocity = (a.vel + a.ang_vel.cross(r_a)) - (b.vel + b.ang_vel.cross(r_b));
    let k_inverse = stable_inverse(&point_constraint_k(&a, &b, r_a, r_b));
    let impulse = k_inverse * (-relative_velocity);
    apply_pair_impulse(bodies, handle_a, handle_b, impulse, r_a, r_b);
    relative_velocity.norm()
}

/// Cancel the relative velocity of the two anchors along the two directions
/// `first` and `second` only (the third direction stays free).
///
/// Returns the size of the relative velocity that was removed.
pub(super) fn solve_point_velocity_in_plane(
    bodies: &mut RigidBodySet,
    handle_a: BodyHandle,
    handle_b: BodyHandle,
    r_a: &Vec3,
    r_b: &Vec3,
    first: &Vec3,
    second: &Vec3,
) -> Real {
    let (Some(a), Some(b)) = (read_body(bodies, handle_a), read_body(bodies, handle_b)) else {
        return 0.0;
    };
    let relative_velocity = (a.vel + a.ang_vel.cross(r_a)) - (b.vel + b.ang_vel.cross(r_b));
    let k = point_constraint_k(&a, &b, r_a, r_b);
    let impulse = constrained_impulse(&k, &[*first, *second], &relative_velocity);
    apply_pair_impulse(bodies, handle_a, handle_b, impulse, r_a, r_b);
    first
        .dot(&relative_velocity)
        .hypot(second.dot(&relative_velocity))
}

/// Cancel the relative angular velocity along the given `rows` (one or two
/// directions; three rows lock the relative rotation completely).
///
/// Returns the size of the relative angular velocity that was removed along
/// the rows.
pub(super) fn solve_angular_rows(
    bodies: &mut RigidBodySet,
    handle_a: BodyHandle,
    handle_b: BodyHandle,
    rows: &[Vec3],
) -> Real {
    let (Some(a), Some(b)) = (read_body(bodies, handle_a), read_body(bodies, handle_b)) else {
        return 0.0;
    };
    let relative_rate = a.ang_vel - b.ang_vel;
    let k = a.inv_inertia + b.inv_inertia;
    let impulse = constrained_impulse(&k, rows, &relative_rate);
    apply_angular_impulse(bodies, handle_a, handle_b, impulse);
    rows.iter()
        .map(|row| row.dot(&relative_rate).powi(2))
        .sum::<Real>()
        .sqrt()
}

/// Rotate a body by the world-frame rotation vector `rotation`, keeping its
/// world-space inertia in step with the new orientation.
fn rotate_body(bodies: &mut RigidBodySet, handle: BodyHandle, rotation: &Vec3) {
    let size_squared = rotation.norm_squared();
    if !size_squared.is_finite() || size_squared < 1e-30 {
        return;
    }
    if let Some(body) = bodies.get_mut(handle) {
        body.transform.rotation = small_rotation_quat(rotation) * body.transform.rotation;
        body.update_world_inertia();
    }
}

/// Position-level orientation correction.
///
/// `error_world` is the world-frame rotation vector (axis times angle) that
/// takes the orientation B *should* have, relative to A, to the orientation it
/// *has*.  The correction removes the fraction `beta` of it: with
/// `K = I_A⁻¹ + I_B⁻¹` the shared multiplier is `λ = K⁻¹·beta·error_world`, and
/// body A rotates by `+I_A⁻¹·λ` while body B rotates by `−I_B⁻¹·λ`, so the
/// relative rotation changes by exactly `−beta·error_world` whatever the
/// inertia split.  Bodies without rotational freedom (static, kinematic or an
/// all-zero inverse inertia) are never rotated and the other body absorbs the
/// whole correction; if neither can rotate nothing happens.
pub(super) fn apply_orientation_correction(
    bodies: &mut RigidBodySet,
    handle_a: BodyHandle,
    handle_b: BodyHandle,
    error_world: Vec3,
    beta: Real,
) {
    let error_size_squared = error_world.norm_squared();
    if !error_size_squared.is_finite() || error_size_squared < 1e-30 {
        return;
    }
    let (Some(a), Some(b)) = (read_body(bodies, handle_a), read_body(bodies, handle_b)) else {
        return;
    };
    let k = a.inv_inertia + b.inv_inertia;
    let multiplier = stable_inverse(&k) * (error_world * beta);
    let rotation_a = a.inv_inertia * multiplier;
    let rotation_b = -(b.inv_inertia * multiplier);
    rotate_body(bodies, handle_a, &rotation_a);
    rotate_body(bodies, handle_b, &rotation_b);
}

/// Wrap an angle into `(-π, π]`.
pub(super) fn wrap_to_pi(angle: Real) -> Real {
    let wrapped = angle.rem_euclid(TAU);
    if wrapped > PI { wrapped - TAU } else { wrapped }
}

/// Hinge angle of body B relative to body A about `local_axis` (a unit vector
/// in A's frame), measured from the `reference` relative rotation `qA⁻¹·qB` at
/// which the hinge angle is zero.  Positive angles are right-handed about the
/// axis; the result lies in `(-π, π]`.
///
/// The angle is the twist of `(qA⁻¹·qB)·reference⁻¹` about the axis (swing-twist
/// decomposition), so it is exact for an ideal hinge and ignores any small
/// off-axis misalignment.
pub(super) fn hinge_angle(
    rotation_a: &Quat,
    rotation_b: &Quat,
    reference: &Quat,
    local_axis: &Vec3,
) -> Real {
    let joint_rotation = (rotation_a.inverse() * rotation_b) * reference.inverse();
    let quaternion = joint_rotation.quaternion();
    let imaginary = Vec3::new(quaternion.i, quaternion.j, quaternion.k);
    wrap_to_pi(2.0 * imaginary.dot(local_axis).atan2(quaternion.w))
}

/// World-frame rotation vector that takes body A's hinge axis onto body B's.
///
/// The hinge axis is `local_axis` in A's frame and `reference⁻¹·local_axis` in
/// B's frame (the two coincide at the reference configuration).  The result is
/// zero for an aligned hinge and has the direction `a × b` and the magnitude of
/// the angle between the axes otherwise, so it plugs into
/// [`apply_orientation_correction`] as `error_world`.
pub(super) fn hinge_alignment_error(
    rotation_a: &Quat,
    rotation_b: &Quat,
    reference: &Quat,
    local_axis: &Vec3,
) -> Vec3 {
    let axis_a = rotation_a * local_axis;
    let axis_b = rotation_b * (reference.inverse() * local_axis);
    let cross = axis_a.cross(&axis_b);
    let sine = cross.norm();
    let cosine = axis_a.dot(&axis_b);
    if sine > 1e-12 {
        cross * (sine.atan2(cosine) / sine)
    } else if cosine < 0.0 {
        perpendicular_basis(&axis_a).0 * PI
    } else {
        Vec3::zeros()
    }
}

/// Implicit spring-damper coefficients over a step `dt`.
///
/// Returns `(γ, bias_scale)` for a spring of stiffness `stiffness` (N/m) and
/// damping `damping` (N·s/m): `γ = 1 / (dt·(c + dt·k))` is the softness of the
/// row and `bias = bias_scale · C` the velocity bias for a stretch `C`
/// (`bias_scale = dt·k·γ = k / (c + dt·k)`).  Negative parameters count as zero.
/// Returns `None` when the step is not positive or the spring has neither
/// stiffness nor damping, so there is no row to solve.
pub(super) fn spring_softness(stiffness: Real, damping: Real, dt: Real) -> Option<(Real, Real)> {
    let stiffness = stiffness.max(0.0);
    let damping = damping.max(0.0);
    let denominator = damping + dt * stiffness;
    if !dt.is_finite() || dt <= 0.0 || !denominator.is_finite() || denominator <= 1e-12 {
        return None;
    }
    let gamma = 1.0 / (dt * denominator);
    Some((gamma, dt * stiffness * gamma))
}

/// A soft (compliant) velocity row along a fixed direction between two anchors.
///
/// The row solves the implicit-Euler spring-damper update as a regularised
/// constraint: with the relative speed `u = (v_B − v_A)·dir` of the anchors, the
/// inverse effective mass `K`, the softness `γ`, the bias `b` and the
/// accumulated *tension* impulse `λ` (positive pulls the anchors together), one
/// Gauss-Seidel update is
///
/// ```text
/// Δλ = (u + b − γ·λ) / (K + γ),   λ += Δλ,   A += Δλ·dir,  B −= Δλ·dir
/// ```
///
/// Because `λ` is accumulated and `γ` regularises against it, repeated updates
/// converge to the same solution however many solver iterations run, unlike an
/// explicit spring impulse that would be applied again on every iteration.
#[derive(Debug, Clone, Copy)]
pub(super) struct SoftRow {
    /// Unit direction from anchor A towards anchor B.
    dir: Vec3,
    /// Inverse effective mass `K` of the row.
    inv_eff_mass: Real,
    /// Softness `γ` (zero for a rigid row).
    gamma: Real,
    /// Velocity bias `b` (`C·k/(c + dt·k)` for a spring with stretch `C`).
    bias: Real,
    /// Accumulated tension impulse `λ`.
    lambda: Real,
    /// Size of the velocity change of the last update.
    last_update: Real,
    /// Whether the row takes part in the solve.
    enabled: bool,
}

impl SoftRow {
    /// A row that does nothing (a spring without stiffness and damping, a
    /// non-positive step, or a joint whose bodies are missing).
    pub(super) fn inactive() -> Self {
        Self {
            dir: Vec3::new(1.0, 0.0, 0.0),
            inv_eff_mass: 0.0,
            gamma: 0.0,
            bias: 0.0,
            lambda: 0.0,
            last_update: 0.0,
            enabled: false,
        }
    }

    /// Build a row along `dir` between the anchors `r_a` and `r_b` with
    /// softness `gamma` and bias `bias`, with a zeroed accumulator.
    pub(super) fn new(
        a: &BodyProps,
        b: &BodyProps,
        r_a: &Vec3,
        r_b: &Vec3,
        dir: Vec3,
        gamma: Real,
        bias: Real,
    ) -> Self {
        Self {
            dir,
            inv_eff_mass: inverse_effective_mass(a, b, r_a, r_b, &dir),
            gamma,
            bias,
            lambda: 0.0,
            last_update: 0.0,
            enabled: true,
        }
    }

    /// Accumulated tension impulse of the current step.
    pub(super) fn lambda(&self) -> Real {
        self.lambda
    }

    /// Velocity change made by the most recent update.
    pub(super) fn last_update(&self) -> Real {
        self.last_update
    }

    /// Run one Gauss-Seidel update, scaled by the relaxation factor `omega`
    /// (`1.0` is the plain update), and return the size of the velocity change.
    pub(super) fn iterate(
        &mut self,
        bodies: &mut RigidBodySet,
        handle_a: BodyHandle,
        handle_b: BodyHandle,
        r_a: &Vec3,
        r_b: &Vec3,
        omega: Real,
    ) -> Real {
        self.last_update = 0.0;
        if !self.enabled {
            return 0.0;
        }
        let (Some(a), Some(b)) = (read_body(bodies, handle_a), read_body(bodies, handle_b)) else {
            return 0.0;
        };
        let denominator = self.inv_eff_mass + self.gamma;
        if !denominator.is_finite() || denominator <= 1e-30 {
            return 0.0;
        }
        let speed =
            ((b.vel + b.ang_vel.cross(r_b)) - (a.vel + a.ang_vel.cross(r_a))).dot(&self.dir);
        let relaxation = if omega > 0.0 && omega.is_finite() {
            omega
        } else {
            1.0
        };
        let delta = relaxation * (speed + self.bias - self.gamma * self.lambda) / denominator;
        self.lambda += delta;
        apply_pair_impulse(bodies, handle_a, handle_b, self.dir * delta, r_a, r_b);
        self.last_update = (delta * self.inv_eff_mass).abs();
        self.last_update
    }
}
