// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

use oxiphysics_articulated::{
    aba::aba,
    body::RigidBody,
    joint::RevoluteJoint,
    joint_limits::{JointLimit, JointLimitSet, aba_with_limits},
    model::ArticulatedModel,
    spatial::{SpatialInertia, SpatialTransform},
};
use std::f64::consts::PI;

fn make_pendulum() -> ArticulatedModel {
    let g = 9.81_f64;
    let mut model = ArticulatedModel::new([0.0, 0.0, -g]);
    let inertia = SpatialInertia::from_com(1.0, [0.0, 0.0, -1.0], [[0.0; 3]; 3]);
    let body = RigidBody::new("pendulum", inertia, None, SpatialTransform::IDENTITY);
    model.add_body(body, Box::new(RevoluteJoint::new([1.0, 0.0, 0.0])));
    model
}

#[test]
fn test_limit_set_no_limit_unchanged_dynamics() {
    let model = make_pendulum();
    let q = vec![PI / 6.0];
    let qd = vec![0.0];
    let tau = vec![0.0];
    let limits = JointLimitSet::for_model(&model);
    let qdd_no_limit = aba(&model, &q, &qd, &tau);
    let qdd_with_limit = aba_with_limits(&model, &q, &qd, &tau, &limits);
    let diff = (qdd_no_limit[0] - qdd_with_limit[0]).abs();
    assert!(
        diff < 1e-15,
        "no limits should not change qdd, diff={diff:.2e}"
    );
}

#[test]
fn test_pendulum_with_soft_limit_at_pi_over_2() {
    let model = make_pendulum();
    let hi = PI / 2.0;
    let mut limits = JointLimitSet::for_model(&model);
    limits.set(
        0,
        0,
        JointLimit {
            lo: -hi,
            hi,
            stiffness: 1e4,
            damping: 50.0,
        },
    );

    // Release from q = PI/3, qd = 0 under gravity; integrate 1000 × 1ms
    let mut q = vec![PI / 3.0];
    let mut qd = vec![0.0];
    let dt = 0.001_f64;
    let tau = vec![0.0];
    let mut max_q = q[0].abs();

    for _ in 0..1000 {
        let qdd = aba_with_limits(&model, &q, &qd, &tau, &limits);
        qd[0] += qdd[0] * dt;
        q[0] += qd[0] * dt;
        if q[0].abs() > max_q {
            max_q = q[0].abs();
        }
    }

    assert!(
        max_q <= hi + 0.02,
        "pendulum should not exceed limit+0.02 rad, max_q={max_q:.4} hi={hi:.4}"
    );
}

#[test]
fn test_pendulum_reverses_direction_at_limit() {
    let model = make_pendulum();
    let hi = PI / 2.0;
    let mut limits = JointLimitSet::for_model(&model);
    limits.set(
        0,
        0,
        JointLimit {
            lo: -hi,
            hi,
            stiffness: 1e4,
            damping: 10.0,
        },
    );

    // Start well above hi with positive velocity (moving into the limit)
    let mut q = vec![hi + 0.05];
    let mut qd = vec![0.5];
    let dt = 0.001_f64;
    let tau = vec![0.0];
    let mut sign_changes = 0;
    let mut last_sign = (qd[0] as f64).signum();

    for _ in 0..500 {
        let qdd = aba_with_limits(&model, &q, &qd, &tau, &limits);
        qd[0] += qdd[0] * dt;
        q[0] += qd[0] * dt;
        let sign = qd[0].signum();
        if sign != last_sign {
            sign_changes += 1;
        }
        last_sign = sign;
    }

    assert!(
        sign_changes >= 1,
        "velocity should reverse direction at least once"
    );
}

#[test]
fn test_limit_set_indexing_two_body_model() {
    // 2-body model: body 0 has 1 DOF (revolute), body 1 has 1 DOF (revolute)
    let g = 9.81_f64;
    let mut model = ArticulatedModel::new([0.0, 0.0, -g]);
    let inertia = SpatialInertia::from_com(1.0, [0.0, 0.0, -0.5], [[0.0; 3]; 3]);
    let b0 = RigidBody::new("link0", inertia, None, SpatialTransform::IDENTITY);
    model.add_body(b0, Box::new(RevoluteJoint::new([1.0, 0.0, 0.0])));
    let b1 = RigidBody::new(
        "link1",
        inertia,
        Some(0),
        SpatialTransform::from_translation([0.0, 0.0, -0.5]),
    );
    model.add_body(b1, Box::new(RevoluteJoint::new([0.0, 1.0, 0.0])));

    let mut limits = JointLimitSet::for_model(&model);
    // Set limit only on body 1, DOF 0 (global dof index = 1)
    limits.set(
        1,
        0,
        JointLimit {
            lo: -0.5,
            hi: 0.5,
            stiffness: 1000.0,
            damping: 10.0,
        },
    );

    // q = [0.0, 1.0] — body 1's joint is well outside its limit
    let q = vec![0.0, 1.0];
    let qd = vec![0.0, 0.0];
    let mut tau = vec![0.0, 0.0];
    limits.apply(&model, &q, &qd, &mut tau);

    // tau[0] should be unchanged (no limit on body 0)
    assert_eq!(tau[0], 0.0, "tau[0] should not be affected");
    // tau[1] should be negative (pushing back below hi=0.5)
    assert!(
        tau[1] < -400.0,
        "tau[1] should be a strong negative penalty, got {}",
        tau[1]
    );
}
