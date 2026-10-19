// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the revolute joint: angle limits, the reference
//! rotation, axis alignment and the exactness of the hinge rows.

use oxiphysics_core::BodyHandle;
use oxiphysics_core::math::{Quat, Real, Vec3};
use oxiphysics_rigid::{RigidBody, RigidBodySet};

use super::RevoluteJoint;
use super::test_support::{body_with_inertia, rotate_body, rotation_angle, step_joint};
use crate::pgs_solver::PgsSolver;
use crate::traits::Constraint;

const DT: Real = 1.0 / 60.0;
const LIMIT: Real = 0.5;

/// Static body A and dynamic unit-mass, unit-inertia body B, both at the
/// origin so that zero anchors coincide.
fn hinge_pair() -> (RigidBodySet, BodyHandle, BodyHandle) {
    let mut bodies = RigidBodySet::new();
    let handle_a = bodies.insert(RigidBody::new_static());
    let handle_b = bodies.insert(RigidBody::new(1.0));
    (bodies, handle_a, handle_b)
}

fn limited_hinge(handle_a: BodyHandle, handle_b: BodyHandle) -> RevoluteJoint {
    RevoluteJoint::new(handle_a, handle_b, Vec3::zeros(), Vec3::zeros(), Vec3::z())
        .with_limits(-LIMIT, LIMIT)
}

#[test]
fn constant_torque_cannot_drive_the_hinge_past_its_limits() {
    for torque in [5.0, -5.0] {
        let (mut bodies, handle_a, handle_b) = hinge_pair();
        let mut joint = limited_hinge(handle_a, handle_b);
        let mut extreme = 0.0_f64;
        for _ in 0..600 {
            if let Some(body) = bodies.get_mut(handle_b) {
                body.apply_torque(Vec3::new(0.0, 0.0, torque));
            }
            step_joint(&mut joint, &mut bodies, &Vec3::zeros(), DT, 8, 4);
            let angle = joint.hinge_angle(&bodies);
            if angle.abs() > extreme.abs() {
                extreme = angle;
            }
        }
        assert!(
            extreme.abs() < LIMIT + 1e-6,
            "torque {torque}: the hinge reached {extreme} rad, past the limit {LIMIT}"
        );
        assert!(
            extreme.abs() > LIMIT - 1e-6,
            "torque {torque}: the hinge should reach its limit, only got {extreme}"
        );
        let spin = bodies.get(handle_b).unwrap().angular_velocity.z;
        assert!(
            spin.abs() < 1e-6,
            "the body must rest on the limit, spin {spin}"
        );
        let pushed_back = if torque > 0.0 {
            joint.upper_limit_impulse()
        } else {
            joint.lower_limit_impulse()
        };
        assert!(
            (pushed_back - torque.abs() * DT).abs() < 1e-6,
            "the limit must carry the torque impulse {}, got {pushed_back}",
            torque.abs() * DT
        );
        assert_eq!(
            if torque > 0.0 {
                joint.lower_limit_impulse()
            } else {
                joint.upper_limit_impulse()
            },
            0.0,
            "the opposite limit stays inactive"
        );
    }
}

#[test]
fn limits_work_through_the_pgs_solver() {
    let (mut bodies, handle_a, handle_b) = hinge_pair();
    let mut constraints: Vec<Box<dyn Constraint>> =
        vec![Box::new(limited_hinge(handle_a, handle_b))];
    let solver = PgsSolver::new(8, 4);
    let mut extreme = 0.0_f64;
    for _ in 0..400 {
        if let Some(body) = bodies.get_mut(handle_b) {
            body.apply_torque(Vec3::new(0.0, 0.0, 8.0));
            body.integrate_forces(DT, &Vec3::zeros());
        }
        solver.solve(&mut constraints, &mut bodies, DT);
        if let Some(body) = bodies.get_mut(handle_b) {
            body.integrate_velocity(DT);
            extreme = extreme.max(body.transform.rotation.scaled_axis().z);
        }
    }
    assert!(
        extreme < LIMIT + 1e-6 && extreme > LIMIT - 1e-6,
        "hinge angle ended at {extreme}, expected the limit {LIMIT}"
    );
}

#[test]
fn a_fast_swing_stops_at_the_limit_instead_of_overshooting_it() {
    // 20 rad/s would carry the hinge 0.33 rad per step, i.e. far past 0.5 rad
    // in two steps; the speculative row hands the joint exactly to the limit.
    let (mut bodies, handle_a, handle_b) = hinge_pair();
    if let Some(body) = bodies.get_mut(handle_b) {
        body.angular_velocity = Vec3::new(0.0, 0.0, 20.0);
    }
    let mut joint = limited_hinge(handle_a, handle_b);
    for step in 0..6 {
        step_joint(&mut joint, &mut bodies, &Vec3::zeros(), DT, 8, 4);
        let angle = joint.hinge_angle(&bodies);
        assert!(
            angle <= LIMIT + 1e-9,
            "step {step}: hinge angle {angle} overshoots the limit"
        );
    }
    assert!((joint.hinge_angle(&bodies) - LIMIT).abs() < 1e-9);
    assert!(bodies.get(handle_b).unwrap().angular_velocity.z.abs() < 1e-9);
}

#[test]
fn limits_never_pull_a_hinge_that_moves_freely_inside_its_range() {
    let (mut bodies, handle_a, handle_b) = hinge_pair();
    if let Some(body) = bodies.get_mut(handle_b) {
        body.angular_velocity = Vec3::new(0.0, 0.0, 1.0);
    }
    let mut joint = limited_hinge(handle_a, handle_b);
    joint.prepare(&bodies, DT);
    joint.solve_velocity(&mut bodies, DT);
    assert_eq!(joint.lower_limit_impulse(), 0.0);
    assert_eq!(joint.upper_limit_impulse(), 0.0);
    let spin = bodies.get(handle_b).unwrap().angular_velocity.z;
    assert!(
        (spin - 1.0).abs() < 1e-12,
        "free spin must be untouched, got {spin}"
    );
}

#[test]
fn a_joint_without_limits_spins_freely() {
    let (mut bodies, handle_a, handle_b) = hinge_pair();
    if let Some(body) = bodies.get_mut(handle_b) {
        body.angular_velocity = Vec3::new(0.0, 0.0, 10.0);
    }
    let mut joint = RevoluteJoint::new(handle_a, handle_b, Vec3::zeros(), Vec3::zeros(), Vec3::z());
    for _ in 0..120 {
        step_joint(&mut joint, &mut bodies, &Vec3::zeros(), DT, 8, 4);
    }
    let spin = bodies.get(handle_b).unwrap().angular_velocity.z;
    assert!((spin - 10.0).abs() < 1e-9, "spin {spin}");
}

#[test]
fn the_position_pass_pushes_a_violation_back_to_the_limit_without_overshoot() {
    let (mut bodies, handle_a, handle_b) = hinge_pair();
    if let Some(body) = bodies.get_mut(handle_b) {
        rotate_body(body, Vec3::z(), 0.7);
    }
    let mut joint = limited_hinge(handle_a, handle_b);
    let mut previous = joint.hinge_angle(&bodies);
    assert!((previous - 0.7).abs() < 1e-12);
    for _ in 0..200 {
        joint.solve_position(&mut bodies, DT);
        let angle = joint.hinge_angle(&bodies);
        assert!(angle <= previous + 1e-15, "{previous} -> {angle}");
        assert!(angle >= LIMIT - 1e-9, "overshoot below the limit: {angle}");
        previous = angle;
    }
    assert!(previous - LIMIT < 1e-6, "ended at {previous}");
    assert_eq!(rotation_angle(&bodies, handle_a), 0.0);
}

#[test]
fn the_position_pass_also_fixes_a_violation_of_the_lower_limit() {
    let (mut bodies, handle_a, handle_b) = hinge_pair();
    if let Some(body) = bodies.get_mut(handle_b) {
        rotate_body(body, Vec3::z(), -0.9);
    }
    let mut joint = limited_hinge(handle_a, handle_b);
    for _ in 0..200 {
        joint.solve_position(&mut bodies, DT);
    }
    let angle = joint.hinge_angle(&bodies);
    assert!(
        (angle + LIMIT).abs() < 1e-6,
        "ended at {angle}, expected {}",
        -LIMIT
    );
}

#[test]
fn the_reference_rotation_defines_the_zero_angle_of_the_limits() {
    let offset = 0.4;
    let build = || {
        let (mut bodies, handle_a, handle_b) = hinge_pair();
        if let Some(body) = bodies.get_mut(handle_b) {
            rotate_body(body, Vec3::z(), offset);
        }
        (bodies, handle_a, handle_b)
    };

    // Without a reference the 0.4 rad is a violation of the +-0.1 rad limits.
    let (mut bodies, handle_a, handle_b) = build();
    let mut plain = RevoluteJoint::new(handle_a, handle_b, Vec3::zeros(), Vec3::zeros(), Vec3::z())
        .with_limits(-0.1, 0.1);
    for _ in 0..100 {
        plain.solve_position(&mut bodies, DT);
    }
    assert!(plain.hinge_angle(&bodies) < 0.2);

    // With the reference at the current pose nothing is violated and nothing moves.
    let (mut bodies, handle_a, handle_b) = build();
    let reference = Quat::from_axis_angle(&Vec3::z_axis(), offset);
    let mut referenced =
        RevoluteJoint::new(handle_a, handle_b, Vec3::zeros(), Vec3::zeros(), Vec3::z())
            .with_limits(-0.1, 0.1)
            .with_reference_rotation(reference);
    assert!(referenced.hinge_angle(&bodies).abs() < 1e-12);
    for _ in 0..100 {
        referenced.solve_position(&mut bodies, DT);
    }
    assert!((rotation_angle(&bodies, handle_b) - offset).abs() < 1e-12);

    // `capture_reference` reads the same pose off the bodies.
    let (bodies, handle_a, handle_b) = build();
    let mut captured =
        RevoluteJoint::new(handle_a, handle_b, Vec3::zeros(), Vec3::zeros(), Vec3::z())
            .with_limits(-0.1, 0.1);
    captured.capture_reference(&bodies);
    assert!(captured.hinge_angle(&bodies).abs() < 1e-12);
    assert!((captured.reference_rotation.angle() - offset).abs() < 1e-12);
}

#[test]
fn hinge_angle_reports_the_signed_rotation_about_the_axis() {
    for angle in [-1.3, -0.2, 0.0, 0.7, 2.5] {
        let (mut bodies, handle_a, handle_b) = hinge_pair();
        if let Some(body) = bodies.get_mut(handle_b) {
            rotate_body(body, Vec3::z(), angle);
        }
        let joint = RevoluteJoint::new(handle_a, handle_b, Vec3::zeros(), Vec3::zeros(), Vec3::z());
        assert!(
            (joint.hinge_angle(&bodies) - angle).abs() < 1e-12,
            "angle {angle}"
        );
    }
}

#[test]
fn the_position_pass_realigns_a_tilted_hinge_axis() {
    let (mut bodies, handle_a, handle_b) = hinge_pair();
    if let Some(body) = bodies.get_mut(handle_b) {
        rotate_body(body, Vec3::x(), 0.2);
    }
    let mut joint = RevoluteJoint::new(handle_a, handle_b, Vec3::zeros(), Vec3::zeros(), Vec3::z());
    let tilt = |bodies: &RigidBodySet| {
        let b = bodies.get(handle_b).unwrap();
        let tilted_axis = b.transform.rotation * Vec3::z();
        tilted_axis
            .cross(&Vec3::z())
            .norm()
            .atan2(tilted_axis.dot(&Vec3::z()))
    };
    let mut previous = tilt(&bodies);
    assert!((previous - 0.2).abs() < 1e-12);
    for _ in 0..100 {
        joint.solve_position(&mut bodies, DT);
        let now = tilt(&bodies);
        assert!(now < previous, "tilt must shrink: {previous} -> {now}");
        previous = now;
    }
    assert!(previous < 1e-3, "tilt left over: {previous}");
    assert_eq!(rotation_angle(&bodies, handle_a), 0.0);
}

#[test]
fn alignment_correction_shares_the_work_between_two_dynamic_bodies() {
    let mut bodies = RigidBodySet::new();
    let handle_a = bodies.insert(body_with_inertia(1.0, Vec3::new(1.0, 1.0, 1.0)));
    let mut body_b = body_with_inertia(1.0, Vec3::new(1.0, 1.0, 1.0));
    rotate_body(&mut body_b, Vec3::x(), 0.2);
    let handle_b = bodies.insert(body_b);
    let mut joint = RevoluteJoint::new(handle_a, handle_b, Vec3::zeros(), Vec3::zeros(), Vec3::z());
    joint.solve_position(&mut bodies, DT);
    // Equal inertia: A and B each take half of the 10 % correction.
    assert!((rotation_angle(&bodies, handle_a) - 0.5 * 0.1 * 0.2).abs() < 1e-12);
    assert!((rotation_angle(&bodies, handle_b) - (0.2 - 0.5 * 0.1 * 0.2)).abs() < 1e-12);
}

#[test]
fn hinge_rows_exert_no_torque_about_the_free_axis() {
    // B has an anisotropic inertia tensor and is tilted, so its world-space
    // inertia is not diagonal.  The locked rows only produce torque
    // perpendicular to the hinge axis: the angular momentum about the axis
    // must not change, whatever the inertia coupling does to the spin.
    let mut bodies = RigidBodySet::new();
    let handle_a = bodies.insert(RigidBody::new_static());
    let mut body_b = body_with_inertia(1.0, Vec3::new(1.0, 2.0, 3.0));
    rotate_body(&mut body_b, Vec3::x(), 0.5);
    body_b.angular_velocity = Vec3::new(0.5, -0.3, 4.0);
    let handle_b = bodies.insert(body_b);

    let axial_momentum = |bodies: &RigidBodySet| {
        let b = bodies.get(handle_b).unwrap();
        let inertia = b
            .world_inverse_inertia
            .try_inverse()
            .expect("the inertia tensor is invertible");
        (inertia * b.angular_velocity).z
    };
    let before = axial_momentum(&bodies);

    let mut joint = RevoluteJoint::new(handle_a, handle_b, Vec3::zeros(), Vec3::zeros(), Vec3::z());
    joint.prepare(&bodies, DT);
    joint.solve_velocity(&mut bodies, DT);

    let after = axial_momentum(&bodies);
    assert!(
        (after - before).abs() < 1e-12,
        "axial angular momentum changed from {before} to {after}"
    );
    // The two rows perpendicular to the axis are locked exactly: what is left
    // of the spin lies along the hinge axis.
    let final_spin = bodies.get(handle_b).unwrap().angular_velocity;
    let (first, second) = super::constraint_math::perpendicular_basis(&Vec3::z());
    assert!(first.dot(&final_spin).abs() < 1e-12, "{final_spin:?}");
    assert!(second.dot(&final_spin).abs() < 1e-12, "{final_spin:?}");
}

#[test]
fn limits_hold_when_both_bodies_are_dynamic() {
    let mut bodies = RigidBodySet::new();
    let handle_a = bodies.insert(body_with_inertia(1.0, Vec3::new(1.0, 1.0, 1.0)));
    let handle_b = bodies.insert(body_with_inertia(1.0, Vec3::new(1.0, 1.0, 1.0)));
    let mut joint = limited_hinge(handle_a, handle_b);
    let mut extreme = 0.0_f64;
    for _ in 0..600 {
        if let Some(body) = bodies.get_mut(handle_b) {
            body.apply_torque(Vec3::new(0.0, 0.0, 3.0));
        }
        step_joint(&mut joint, &mut bodies, &Vec3::zeros(), DT, 12, 4);
        extreme = extreme.max(joint.hinge_angle(&bodies));
    }
    assert!(
        extreme < LIMIT + 1e-4,
        "relative angle {extreme} exceeded the limit {LIMIT}"
    );
    assert!(extreme > LIMIT - 1e-3, "the hinge should reach its limit");
}
