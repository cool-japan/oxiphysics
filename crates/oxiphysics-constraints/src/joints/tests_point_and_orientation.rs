// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the point-constraint effective mass and the
//! orientation correction of the position pass.

use std::f64::consts::FRAC_PI_2;

use oxiphysics_core::BodyHandle;
use oxiphysics_core::math::{Mat3, Quat, Real, Vec3};
use oxiphysics_rigid::{RigidBody, RigidBodySet};

use super::constraint_math::{
    apply_orientation_correction, constrained_impulse, hinge_alignment_error, hinge_angle,
    point_constraint_k, spring_softness, stable_inverse, wrap_to_pi,
};
use super::functions::{apply_pair_impulse, read_body};
use super::test_support::{
    anchor_relative_velocity, body_with_inertia, linear_momentum, relative_rotation_angle,
    rotate_body, rotation_angle,
};
use super::{
    BallJoint, ConeJoint, CylindricalJoint, FixedJoint, PrismaticJoint, RevoluteJoint,
    SphericalJoint, UniversalJoint, WeldJoint,
};
use crate::traits::Constraint;

const DT: Real = 1.0 / 60.0;

fn insert(bodies: &mut RigidBodySet, body: RigidBody) -> BodyHandle {
    bodies.insert(body)
}

/// Two dynamic bodies with different masses and anisotropic inertia, both
/// moving and spinning, whose anchors `(0.5, 0.2, 0)` / `(-0.5, 0.2, 0)` meet
/// at the world point `(0.5, 0.2, 0)`.
fn offset_anchor_scene() -> (RigidBodySet, BodyHandle, BodyHandle, Vec3, Vec3) {
    let mut bodies = RigidBodySet::new();
    let mut body_a = body_with_inertia(2.0, Vec3::new(0.4, 0.9, 1.3));
    body_a.angular_velocity = Vec3::new(0.3, -0.2, 0.9);
    let mut body_b = body_with_inertia(0.5, Vec3::new(0.2, 0.1, 0.3));
    body_b.transform.position = Vec3::new(1.0, 0.0, 0.0);
    body_b.velocity = Vec3::new(1.0, -2.0, 0.5);
    body_b.angular_velocity = Vec3::new(-0.7, 0.4, 0.1);
    let handle_a = insert(&mut bodies, body_a);
    let handle_b = insert(&mut bodies, body_b);
    (
        bodies,
        handle_a,
        handle_b,
        Vec3::new(0.5, 0.2, 0.0),
        Vec3::new(-0.5, 0.2, 0.0),
    )
}

// ---------------------------------------------------------------------------
// F0384: point constraints use the full 3x3 effective mass
// ---------------------------------------------------------------------------

#[test]
fn point_constraint_k_predicts_the_anchor_velocity_response_to_an_impulse() {
    let (mut bodies, handle_a, handle_b, local_a, local_b) = offset_anchor_scene();
    // Tilt both bodies so the world-space inertia is not diagonal.
    if let Some(body) = bodies.get_mut(handle_a) {
        rotate_body(body, Vec3::new(1.0, 1.0, 0.0), 0.7);
    }
    if let Some(body) = bodies.get_mut(handle_b) {
        rotate_body(body, Vec3::new(0.0, 1.0, 1.0), -0.4);
    }
    let a = read_body(&bodies, handle_a).expect("body A");
    let b = read_body(&bodies, handle_b).expect("body B");
    let r_a = a.rotation * local_a;
    let r_b = b.rotation * local_b;
    let k = point_constraint_k(&a, &b, &r_a, &r_b);

    let impulse = Vec3::new(0.3, -1.1, 0.7);
    let before = anchor_relative_velocity(&bodies, handle_a, handle_b, local_a, local_b);
    apply_pair_impulse(&mut bodies, handle_a, handle_b, impulse, &r_a, &r_b);
    let after = anchor_relative_velocity(&bodies, handle_a, handle_b, local_a, local_b);

    let predicted = k * impulse;
    assert!(
        ((after - before) - predicted).norm() < 1e-12,
        "K must map an impulse to the change of the relative anchor velocity: \
         actual {:?} vs K*P {:?}",
        after - before,
        predicted
    );
    assert!(
        (k - k.transpose()).norm() < 1e-12,
        "the effective-mass matrix is symmetric"
    );
}

#[test]
fn ball_joint_with_offset_anchors_cancels_anchor_velocity_in_one_iteration() {
    let (mut bodies, handle_a, handle_b, local_a, local_b) = offset_anchor_scene();
    let momentum_before = linear_momentum(&bodies);
    let mut joint = BallJoint::new(handle_a, handle_b, local_a, local_b);
    joint.prepare(&bodies, DT);
    let initial = anchor_relative_velocity(&bodies, handle_a, handle_b, local_a, local_b);
    assert!(
        initial.norm() > 1.0,
        "the scene must start with a real error"
    );
    joint.solve_velocity(&mut bodies, DT);
    let residual = anchor_relative_velocity(&bodies, handle_a, handle_b, local_a, local_b);
    assert!(
        residual.norm() < 1e-9,
        "one iteration must cancel the anchor velocity exactly, residual {residual:?}"
    );
    let momentum_after = linear_momentum(&bodies);
    assert!(
        (momentum_after - momentum_before).norm() < 1e-12,
        "the joint impulse is internal: linear momentum is conserved"
    );
}

#[test]
fn spherical_joint_adds_inertia_with_the_right_units() {
    // The old matrix added 1/(kg m^2) to 1/kg; with a heavy body of tiny
    // inertia the result was far from the exact solution.
    let (mut bodies, handle_a, handle_b, local_a, local_b) = offset_anchor_scene();
    let mut joint = SphericalJoint::new(
        handle_a,
        handle_b,
        local_a,
        local_b,
        Vec3::new(0.0, 0.0, 1.0),
    );
    joint.prepare(&bodies, DT);
    joint.solve_velocity(&mut bodies, DT);
    let residual = anchor_relative_velocity(&bodies, handle_a, handle_b, local_a, local_b);
    assert!(residual.norm() < 1e-9, "residual {residual:?}");
}

#[test]
fn point_joints_converge_to_zero_anchor_velocity() {
    // The cone joint only has the point rows, so one iteration is exact.
    let (mut bodies, handle_a, handle_b, local_a, local_b) = offset_anchor_scene();
    let mut cone = ConeJoint::new(
        handle_a,
        handle_b,
        local_a,
        local_b,
        Vec3::new(1.0, 0.0, 0.0),
        1.0,
    );
    cone.prepare(&bodies, DT);
    cone.solve_velocity(&mut bodies, DT);
    let residual = anchor_relative_velocity(&bodies, handle_a, handle_b, local_a, local_b);
    assert!(residual.norm() < 1e-9, "cone residual {residual:?}");

    // Revolute and universal joints solve the point rows first and then lock
    // rotation rows that move the anchors again, so they converge over the
    // solver iterations.
    let (mut bodies, handle_a, handle_b, local_a, local_b) = offset_anchor_scene();
    let mut revolute = RevoluteJoint::new(handle_a, handle_b, local_a, local_b, Vec3::z());
    revolute.prepare(&bodies, DT);
    for _ in 0..400 {
        revolute.solve_velocity(&mut bodies, DT);
    }
    let residual = anchor_relative_velocity(&bodies, handle_a, handle_b, local_a, local_b);
    assert!(residual.norm() < 1e-6, "revolute residual {residual:?}");

    let (mut bodies, handle_a, handle_b, local_a, local_b) = offset_anchor_scene();
    let mut universal =
        UniversalJoint::new(handle_a, handle_b, local_a, local_b, Vec3::x(), Vec3::z());
    universal.prepare(&bodies, DT);
    for _ in 0..400 {
        universal.solve_velocity(&mut bodies, DT);
    }
    let residual = anchor_relative_velocity(&bodies, handle_a, handle_b, local_a, local_b);
    assert!(residual.norm() < 1e-6, "universal residual {residual:?}");
}

#[test]
fn prismatic_joint_locks_perpendicular_motion_and_rotation_but_keeps_sliding() {
    let mut bodies = RigidBodySet::new();
    let handle_a = insert(
        &mut bodies,
        body_with_inertia(1.0, Vec3::new(0.5, 0.8, 1.1)),
    );
    let mut body_b = body_with_inertia(3.0, Vec3::new(0.3, 0.6, 0.2));
    body_b.transform.position = Vec3::new(1.0, 0.0, 0.0);
    body_b.velocity = Vec3::new(3.0, 4.0, -2.0);
    body_b.angular_velocity = Vec3::new(0.5, -0.5, 1.0);
    let handle_b = insert(&mut bodies, body_b);
    // Anchors off the centres so the lever arms matter.
    let local_a = Vec3::new(0.5, 0.3, 0.0);
    let local_b = Vec3::new(-0.5, 0.3, 0.0);
    let mut joint = PrismaticJoint::new(
        handle_a,
        handle_b,
        local_a,
        local_b,
        Vec3::new(1.0, 0.0, 0.0),
    );
    joint.prepare(&bodies, DT);
    for _ in 0..300 {
        joint.solve_velocity(&mut bodies, DT);
    }
    let relative = anchor_relative_velocity(&bodies, handle_a, handle_b, local_a, local_b);
    assert!(
        relative.y.abs() < 1e-9 && relative.z.abs() < 1e-9,
        "perpendicular anchor velocity must vanish, got {relative:?}"
    );
    assert!(
        relative.x.abs() > 0.1,
        "sliding along the axis stays free, got {relative:?}"
    );
    let spin = bodies.get(handle_a).unwrap().angular_velocity
        - bodies.get(handle_b).unwrap().angular_velocity;
    assert!(
        spin.norm() < 1e-9,
        "relative spin must be locked, got {spin:?}"
    );
}

#[test]
fn cylindrical_joint_keeps_slide_and_spin_about_its_axis_free() {
    let mut bodies = RigidBodySet::new();
    let handle_a = insert(
        &mut bodies,
        body_with_inertia(1.0, Vec3::new(0.5, 0.5, 0.5)),
    );
    let mut body_b = body_with_inertia(2.0, Vec3::new(0.5, 0.5, 0.5));
    body_b.transform.position = Vec3::new(1.0, 0.0, 0.0);
    body_b.velocity = Vec3::new(2.0, 3.0, 1.0);
    body_b.angular_velocity = Vec3::new(4.0, 1.0, -1.0);
    let handle_b = insert(&mut bodies, body_b);
    let mut joint = CylindricalJoint::new(
        handle_a,
        handle_b,
        Vec3::zeros(),
        Vec3::zeros(),
        Vec3::new(1.0, 0.0, 0.0),
    );
    joint.prepare(&bodies, DT);
    for _ in 0..100 {
        joint.solve_velocity(&mut bodies, DT);
    }
    let a = bodies.get(handle_a).unwrap();
    let b = bodies.get(handle_b).unwrap();
    let relative_velocity = a.velocity - b.velocity;
    assert!(
        relative_velocity.y.abs() < 1e-9 && relative_velocity.z.abs() < 1e-9,
        "perpendicular velocity must vanish, got {relative_velocity:?}"
    );
    let relative_spin = a.angular_velocity - b.angular_velocity;
    assert!(
        relative_spin.y.abs() < 1e-9 && relative_spin.z.abs() < 1e-9,
        "perpendicular spin must vanish, got {relative_spin:?}"
    );
    assert!(
        (relative_velocity.x + 2.0).abs() < 1e-9,
        "sliding along the axis is untouched, got {relative_velocity:?}"
    );
    assert!(
        (relative_spin.x + 4.0).abs() < 1e-9,
        "spin about the axis is untouched, got {relative_spin:?}"
    );
}

#[test]
fn point_rows_do_nothing_when_neither_body_can_move() {
    let mut bodies = RigidBodySet::new();
    let handle_a = insert(&mut bodies, RigidBody::new_static());
    let mut body_b = RigidBody::new_static();
    body_b.transform.position = Vec3::new(1.0, 0.0, 0.0);
    let handle_b = insert(&mut bodies, body_b);
    let mut joint = BallJoint::new(
        handle_a,
        handle_b,
        Vec3::new(0.5, 0.0, 0.0),
        Vec3::new(-0.5, 0.0, 0.0),
    );
    joint.prepare(&bodies, DT);
    joint.solve_velocity(&mut bodies, DT);
    for handle in [handle_a, handle_b] {
        let body = bodies.get(handle).unwrap();
        assert_eq!(body.velocity, Vec3::zeros());
        assert_eq!(body.angular_velocity, Vec3::zeros());
    }
}

#[test]
fn stable_inverse_handles_singular_and_regular_matrices() {
    let regular = Mat3::from_diagonal(&Vec3::new(2.0, 4.0, 8.0));
    let inverse = stable_inverse(&regular);
    assert!((inverse * regular - Mat3::identity()).norm() < 1e-12);

    // Rank two: no direction along z.
    let rank_two = Mat3::from_diagonal(&Vec3::new(2.0, 4.0, 0.0));
    let pseudo = stable_inverse(&rank_two);
    assert!((pseudo - Mat3::from_diagonal(&Vec3::new(0.5, 0.25, 0.0))).norm() < 1e-12);

    // Numerically singular (determinant ~ 1e-30): must not blow up.
    let nearly_singular = Mat3::from_diagonal(&Vec3::new(1.0, 1.0, 1e-30));
    assert!(stable_inverse(&nearly_singular).norm() < 10.0);

    assert_eq!(stable_inverse(&Mat3::zeros()), Mat3::zeros());
}

#[test]
fn constrained_impulse_removes_exactly_the_rows_it_is_given() {
    // A non-diagonal positive definite K.
    let k = Mat3::new(2.0, 0.5, 0.1, 0.5, 1.5, -0.3, 0.1, -0.3, 1.0);
    let rate = Vec3::new(0.7, -1.2, 0.4);
    let rows = [Vec3::new(0.0, 1.0, 0.0), Vec3::new(0.0, 0.0, 1.0)];
    let impulse = constrained_impulse(&k, &rows, &rate);
    let after = rate + k * impulse;
    assert!(after.y.abs() < 1e-12 && after.z.abs() < 1e-12, "{after:?}");
    assert!(
        impulse.x.abs() < 1e-12,
        "the impulse stays inside the row span"
    );

    let single = constrained_impulse(&k, &[Vec3::new(1.0, 0.0, 0.0)], &rate);
    assert!((rate + k * single).x.abs() < 1e-12);
    assert!(single.y.abs() < 1e-12 && single.z.abs() < 1e-12);

    // Rank-deficient pair of rows: the duplicated direction is solved once.
    let duplicated = constrained_impulse(&k, &[rows[0], rows[0]], &rate);
    assert!((rate + k * duplicated).y.abs() < 1e-12);
}

// ---------------------------------------------------------------------------
// F0381: orientation correction has the right sign and units
// ---------------------------------------------------------------------------

/// Static body A at the origin and dynamic body B (unit inertia) rotated by
/// `angle` about `axis`.
fn pre_rotated_pair(axis: Vec3, angle: Real) -> (RigidBodySet, BodyHandle, BodyHandle) {
    let mut bodies = RigidBodySet::new();
    let handle_a = insert(&mut bodies, RigidBody::new_static());
    let mut body_b = RigidBody::new(1.0);
    rotate_body(&mut body_b, axis, angle);
    let handle_b = insert(&mut bodies, body_b);
    (bodies, handle_a, handle_b)
}

/// Joints stop correcting below a 1e-6 dead-band, so the error may plateau
/// there; above it every iteration must shrink it.
const DEAD_BAND: Real = 2e-6;

fn assert_monotone_convergence(errors: &[Real], final_bound: Real) {
    for pair in errors.windows(2) {
        assert!(
            pair[1] < pair[0] || pair[0] < DEAD_BAND,
            "the error must shrink every iteration: {} -> {}",
            pair[0],
            pair[1]
        );
    }
    let last = errors.last().copied().unwrap_or(Real::MAX);
    assert!(last < final_bound, "final error {last} >= {final_bound}");
}

#[test]
fn fixed_joint_removes_a_pre_rotation_monotonically_without_rotating_the_static_body() {
    let (mut bodies, handle_a, handle_b) = pre_rotated_pair(Vec3::new(0.0, 0.0, 1.0), 0.3);
    let mut joint = FixedJoint::new(
        handle_a,
        handle_b,
        Vec3::zeros(),
        Vec3::zeros(),
        Quat::identity(),
    );
    let mut errors = vec![relative_rotation_angle(&bodies, handle_a, handle_b)];
    for _ in 0..100 {
        joint.solve_position(&mut bodies, DT);
        errors.push(relative_rotation_angle(&bodies, handle_a, handle_b));
    }
    // A static anchor takes no share: each pass removes exactly 10 % (the
    // Baumgarte factor) of the error.
    assert!(
        (errors[1] - 0.3 * 0.9).abs() < 1e-12,
        "first pass must leave 90 % of 0.3 rad, got {}",
        errors[1]
    );
    assert_monotone_convergence(&errors, 1e-3);
    assert_eq!(
        rotation_angle(&bodies, handle_a),
        0.0,
        "the static body must not rotate"
    );
}

#[test]
fn fixed_joint_shares_the_orientation_correction_by_inverse_inertia() {
    // A has inertia 1 (inverse 1), B has inertia 3 (inverse 1/3): K = 4/3 and
    // A takes three quarters of the correction, B one quarter.
    let mut bodies = RigidBodySet::new();
    let handle_a = insert(
        &mut bodies,
        body_with_inertia(1.0, Vec3::new(1.0, 1.0, 1.0)),
    );
    let mut body_b = body_with_inertia(1.0, Vec3::new(3.0, 3.0, 3.0));
    rotate_body(&mut body_b, Vec3::new(0.0, 0.0, 1.0), 0.3);
    let handle_b = insert(&mut bodies, body_b);
    let mut joint = FixedJoint::new(
        handle_a,
        handle_b,
        Vec3::zeros(),
        Vec3::zeros(),
        Quat::identity(),
    );
    joint.solve_position(&mut bodies, DT);
    let beta = 0.1;
    let removed = beta * 0.3;
    let angle_a = bodies
        .get(handle_a)
        .unwrap()
        .transform
        .rotation
        .scaled_axis();
    let angle_b = bodies
        .get(handle_b)
        .unwrap()
        .transform
        .rotation
        .scaled_axis();
    assert!(
        (angle_a.z - 0.75 * removed).abs() < 1e-12,
        "A turns towards B by 3/4 of the correction, got {}",
        angle_a.z
    );
    assert!(
        (angle_b.z - (0.3 - 0.25 * removed)).abs() < 1e-12,
        "B turns back by 1/4 of the correction, got {}",
        angle_b.z
    );
}

#[test]
fn fixed_joint_respects_a_non_identity_target_rotation() {
    let target = Quat::from_axis_angle(&Vec3::y_axis(), 0.5);
    let mut bodies = RigidBodySet::new();
    let handle_a = insert(&mut bodies, RigidBody::new_static());
    let mut body_b = RigidBody::new(1.0);
    // B sits at the target orientation plus an extra 0.2 rad about x.
    body_b.transform.rotation = Quat::from_axis_angle(&Vec3::x_axis(), 0.2) * target;
    body_b.update_world_inertia();
    let handle_b = insert(&mut bodies, body_b);
    let mut joint = FixedJoint::new(handle_a, handle_b, Vec3::zeros(), Vec3::zeros(), target);
    let deviation = |bodies: &RigidBodySet| {
        let b = bodies.get(handle_b).unwrap();
        (target.inverse() * b.transform.rotation).angle()
    };
    let mut errors = vec![deviation(&bodies)];
    for _ in 0..120 {
        joint.solve_position(&mut bodies, DT);
        errors.push(deviation(&bodies));
    }
    assert_monotone_convergence(&errors, 1e-4);
}

#[test]
fn weld_joint_removes_a_pre_rotation_monotonically_without_rotating_the_static_body() {
    let desired = Quat::from_axis_angle(&Vec3::z_axis(), 0.2);
    let mut bodies = RigidBodySet::new();
    let handle_a = insert(&mut bodies, RigidBody::new_static());
    let mut body_b = RigidBody::new(1.0);
    body_b.transform.rotation = Quat::from_axis_angle(&Vec3::x_axis(), 0.25) * desired;
    body_b.update_world_inertia();
    let handle_b = insert(&mut bodies, body_b);
    let mut joint = WeldJoint::new(handle_a, handle_b, Vec3::zeros(), desired);
    let deviation = |bodies: &RigidBodySet| {
        let b = bodies.get(handle_b).unwrap();
        (desired.inverse() * b.transform.rotation).angle()
    };
    let mut errors = vec![deviation(&bodies)];
    for _ in 0..100 {
        joint.solve_position(&mut bodies, DT);
        errors.push(deviation(&bodies));
    }
    assert_monotone_convergence(&errors, 1e-3);
    assert_eq!(rotation_angle(&bodies, handle_a), 0.0);
}

#[test]
fn cone_joint_pulls_the_axis_back_to_the_cone_boundary_without_crossing_it() {
    let half_angle = 0.3;
    let (mut bodies, handle_a, handle_b) = pre_rotated_pair(Vec3::new(0.0, 1.0, 0.0), 0.8);
    let mut joint = ConeJoint::new(
        handle_a,
        handle_b,
        Vec3::zeros(),
        Vec3::zeros(),
        Vec3::new(1.0, 0.0, 0.0),
        half_angle,
    );
    let swing = |bodies: &RigidBodySet| joint.current_swing_angle(bodies);
    let mut previous = swing(&bodies);
    assert!((previous - 0.8).abs() < 1e-12);
    for _ in 0..200 {
        joint.solve_position(&mut bodies, DT);
        let now = joint.current_swing_angle(&bodies);
        assert!(
            now <= previous + 1e-15,
            "swing must not grow: {previous} -> {now}"
        );
        assert!(
            now >= half_angle - 1e-9,
            "the correction must not overshoot the boundary, swing {now}"
        );
        previous = now;
    }
    assert!(
        previous - half_angle < 1e-5,
        "swing must end at the boundary, got {previous}"
    );
    assert_eq!(rotation_angle(&bodies, handle_a), 0.0);
}

#[test]
fn cone_joint_leaves_an_axis_inside_the_cone_alone() {
    let (mut bodies, handle_a, handle_b) = pre_rotated_pair(Vec3::new(0.0, 1.0, 0.0), 0.2);
    let mut joint = ConeJoint::new(
        handle_a,
        handle_b,
        Vec3::zeros(),
        Vec3::zeros(),
        Vec3::new(1.0, 0.0, 0.0),
        0.5,
    );
    joint.solve_position(&mut bodies, DT);
    assert!((joint.current_swing_angle(&bodies) - 0.2).abs() < 1e-12);
}

#[test]
fn universal_joint_restores_perpendicular_axes_from_either_side() {
    for tilt in [0.2, -0.2] {
        // Axis A = x (body A), axis B = z (body B); tilting B about y moves its
        // axis towards or away from x, so the axes stop being perpendicular.
        let (mut bodies, handle_a, handle_b) = pre_rotated_pair(Vec3::new(0.0, 1.0, 0.0), tilt);
        let mut joint = UniversalJoint::new(
            handle_a,
            handle_b,
            Vec3::zeros(),
            Vec3::zeros(),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        );
        let dot = |bodies: &RigidBodySet| {
            let a = bodies.get(handle_a).unwrap();
            let b = bodies.get(handle_b).unwrap();
            (a.transform.rotation * Vec3::new(1.0, 0.0, 0.0))
                .dot(&(b.transform.rotation * Vec3::new(0.0, 0.0, 1.0)))
                .abs()
        };
        let mut errors = vec![dot(&bodies)];
        assert!((errors[0] - tilt.abs().sin()).abs() < 1e-12);
        for _ in 0..120 {
            joint.solve_position(&mut bodies, DT);
            errors.push(dot(&bodies));
        }
        assert_monotone_convergence(&errors, 1e-4);
        assert_eq!(rotation_angle(&bodies, handle_a), 0.0);
    }
}

#[test]
fn orientation_correction_moves_neither_body_when_both_are_static() {
    let mut bodies = RigidBodySet::new();
    let handle_a = insert(&mut bodies, RigidBody::new_static());
    let mut body_b = RigidBody::new_static();
    body_b.transform.rotation = Quat::from_axis_angle(&Vec3::z_axis(), 0.4);
    let handle_b = insert(&mut bodies, body_b);
    apply_orientation_correction(
        &mut bodies,
        handle_a,
        handle_b,
        Vec3::new(0.0, 0.0, 0.4),
        0.5,
    );
    assert_eq!(rotation_angle(&bodies, handle_a), 0.0);
    assert!((rotation_angle(&bodies, handle_b) - 0.4).abs() < 1e-15);
}

#[test]
fn orientation_correction_ignores_non_finite_errors() {
    let (mut bodies, handle_a, handle_b) = pre_rotated_pair(Vec3::new(0.0, 0.0, 1.0), 0.3);
    apply_orientation_correction(
        &mut bodies,
        handle_a,
        handle_b,
        Vec3::new(Real::NAN, 0.0, 0.0),
        0.5,
    );
    apply_orientation_correction(
        &mut bodies,
        handle_a,
        handle_b,
        Vec3::new(Real::INFINITY, 0.0, 0.0),
        0.5,
    );
    let rotation = bodies.get(handle_b).unwrap().transform.rotation;
    assert!(rotation.coords.iter().all(|c| c.is_finite()));
    assert!((rotation.angle() - 0.3).abs() < 1e-15);
}

#[test]
fn a_joint_with_a_degenerate_axis_never_corrupts_the_body_rotation() {
    let (mut bodies, handle_a, handle_b) = pre_rotated_pair(Vec3::new(0.0, 0.0, 1.0), 0.3);
    let mut joint = RevoluteJoint::new(
        handle_a,
        handle_b,
        Vec3::zeros(),
        Vec3::zeros(),
        Vec3::zeros(),
    );
    joint.prepare(&bodies, DT);
    joint.solve_velocity(&mut bodies, DT);
    joint.solve_position(&mut bodies, DT);
    let body = bodies.get(handle_b).unwrap();
    assert!(body.transform.rotation.coords.iter().all(|c| c.is_finite()));
    assert!(body.angular_velocity.iter().all(|c| c.is_finite()));
    assert!((body.transform.rotation.angle() - 0.3).abs() < 1e-15);
}

// ---------------------------------------------------------------------------
// Hinge angle helpers and spring coefficients
// ---------------------------------------------------------------------------

#[test]
fn hinge_angle_is_the_signed_twist_about_the_axis() {
    let axis = Vec3::new(0.0, 0.0, 1.0);
    let identity = Quat::identity();
    for angle in [-3.0, -1.2, -0.01, 0.0, 0.4, 2.0, 3.1] {
        let rotated = Quat::from_axis_angle(&Vec3::z_axis(), angle);
        let measured = hinge_angle(&identity, &rotated, &identity, &axis);
        assert!(
            (measured - angle).abs() < 1e-12,
            "angle {angle} measured as {measured}"
        );
    }
    // A rotation about another axis has no twist about z.
    let swing = Quat::from_axis_angle(&Vec3::x_axis(), 0.7);
    assert!(hinge_angle(&identity, &swing, &identity, &axis).abs() < 1e-12);
    // The reference rotation moves the zero point.
    let reference = Quat::from_axis_angle(&Vec3::z_axis(), 0.5);
    let rotated = Quat::from_axis_angle(&Vec3::z_axis(), 0.9);
    assert!((hinge_angle(&identity, &rotated, &reference, &axis) - 0.4).abs() < 1e-12);
}

#[test]
fn hinge_angle_is_relative_to_body_a() {
    let axis = Vec3::new(0.0, 1.0, 0.0);
    let body_a = Quat::from_axis_angle(&Vec3::y_axis(), 0.6);
    let body_b = Quat::from_axis_angle(&Vec3::y_axis(), 1.0);
    let measured = hinge_angle(&body_a, &body_b, &Quat::identity(), &axis);
    assert!((measured - 0.4).abs() < 1e-12, "got {measured}");
}

#[test]
fn hinge_alignment_error_points_from_axis_a_to_axis_b() {
    let identity = Quat::identity();
    let axis = Vec3::new(0.0, 0.0, 1.0);
    assert!(hinge_alignment_error(&identity, &identity, &identity, &axis).norm() < 1e-15);
    // Tilting B about +x by 0.3 takes its z axis towards -y, a rotation about +x.
    let tilted = Quat::from_axis_angle(&Vec3::x_axis(), 0.3);
    let error = hinge_alignment_error(&identity, &tilted, &identity, &axis);
    assert!(
        (error - Vec3::new(0.3, 0.0, 0.0)).norm() < 1e-12,
        "{error:?}"
    );
    // Anti-parallel axes still give a half-turn about some perpendicular axis.
    let flipped = Quat::from_axis_angle(&Vec3::x_axis(), std::f64::consts::PI);
    let half_turn = hinge_alignment_error(&identity, &flipped, &identity, &axis);
    assert!((half_turn.norm() - std::f64::consts::PI).abs() < 1e-9);
}

#[test]
fn wrap_to_pi_maps_into_the_principal_range() {
    assert!((wrap_to_pi(0.0)).abs() < 1e-15);
    assert!((wrap_to_pi(3.0 * FRAC_PI_2) + FRAC_PI_2).abs() < 1e-12);
    assert!((wrap_to_pi(-3.0 * FRAC_PI_2) - FRAC_PI_2).abs() < 1e-12);
    assert!((wrap_to_pi(std::f64::consts::PI) - std::f64::consts::PI).abs() < 1e-12);
    assert!((wrap_to_pi(7.0) - (7.0 - std::f64::consts::TAU)).abs() < 1e-12);
}

#[test]
fn spring_softness_follows_the_implicit_euler_coefficients() {
    let (gamma, bias_scale) = spring_softness(100.0, 0.0, DT).expect("a spring has a row");
    assert!((gamma - 1.0 / (DT * DT * 100.0)).abs() < 1e-9);
    // bias_scale = k / (c + dt k) = 1 / dt for an undamped spring: the whole
    // stretch is removed within one step.
    assert!((bias_scale - 1.0 / DT).abs() < 1e-9);

    let (gamma, bias_scale) = spring_softness(100.0, 20.0, DT).expect("a damped spring has a row");
    let expected_gamma = 1.0 / (DT * (20.0 + DT * 100.0));
    assert!((gamma - expected_gamma).abs() < 1e-9);
    assert!((bias_scale - 100.0 / (20.0 + DT * 100.0)).abs() < 1e-9);

    assert!(spring_softness(0.0, 0.0, DT).is_none());
    assert!(spring_softness(100.0, 5.0, 0.0).is_none());
    assert!(spring_softness(100.0, 5.0, -DT).is_none());
    // Negative parameters count as zero.
    assert!(spring_softness(-1.0, -1.0, DT).is_none());
}

#[test]
fn rotated_pair_helper_places_body_b_at_the_requested_angle() {
    let (bodies, handle_a, handle_b) = pre_rotated_pair(Vec3::new(1.0, 2.0, 3.0), 0.9);
    assert!((relative_rotation_angle(&bodies, handle_a, handle_b) - 0.9).abs() < 1e-12);
}
