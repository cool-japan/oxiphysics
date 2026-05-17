// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

use oxiphysics_articulated::{
    aba::aba,
    body::RigidBody,
    joint::{HelicalJoint, Joint, RevoluteJoint},
    model::ArticulatedModel,
    spatial::{SpatialInertia, SpatialTransform},
};
use std::f64::consts::PI;

#[test]
fn test_helical_pitch_zero_reduces_to_revolute() {
    let helical = HelicalJoint::new([0.0, 0.0, 1.0], 0.0);
    let revolute = RevoluteJoint::new([0.0, 0.0, 1.0]);
    for k in 0..16 {
        let angle = k as f64 * PI / 8.0;
        let xh = helical.transform(&[angle]);
        let xr = revolute.transform(&[angle]);
        for i in 0..3 {
            for j in 0..3 {
                let diff = (xh.rot[i][j] - xr.rot[i][j]).abs();
                assert!(
                    diff < 1e-14,
                    "rot[{i}][{j}] at angle={angle:.3}: helical={:.8} revolute={:.8}",
                    xh.rot[i][j],
                    xr.rot[i][j]
                );
            }
        }
        for i in 0..3 {
            assert!(
                xh.trans[i].abs() < 1e-14,
                "trans[{i}] at angle={angle:.3}: {:.2e}",
                xh.trans[i]
            );
        }
    }
}

#[test]
fn test_helical_large_pitch_dominantly_translates() {
    let j = HelicalJoint::new([1.0, 0.0, 0.0], 10.0);
    let x = j.transform(&[1.0]);
    let expected_tx = 10.0;
    assert!(
        (x.trans[0] - expected_tx).abs() < 1e-12,
        "tx should be {expected_tx}, got {}",
        x.trans[0]
    );
}

#[test]
fn test_helical_motion_subspace_coupled() {
    let j = HelicalJoint::new([0.0, 0.0, 1.0], 0.5);
    let s = j.motion_subspace();
    assert_eq!(s.len(), 1);
    // angular = [0,0,1], linear = [0,0,0.5]
    for i in 0..3 {
        assert!(
            (s[0].angular[i] - [0.0, 0.0, 1.0][i]).abs() < 1e-15,
            "angular[{i}]"
        );
        assert!(
            (s[0].linear[i] - [0.0, 0.0, 0.5][i]).abs() < 1e-15,
            "linear[{i}]"
        );
    }
}

#[test]
fn test_helical_from_turn_pitch_units() {
    // 5 mm / turn = 0.005 m / turn
    let j = HelicalJoint::from_turn_pitch([0.0, 0.0, 1.0], 0.005);
    let x = j.transform(&[2.0 * PI]); // one full turn
    assert!(
        (x.trans[2] - 0.005).abs() < 1e-12,
        "one turn should advance 5mm, got {}",
        x.trans[2]
    );
}

#[test]
fn test_helical_dynamics_coupling() {
    // A body attached via a helical joint along z.
    // Under a torque about z, we get both rotation AND translation.
    let mut model = ArticulatedModel::new([0.0; 3]);
    let inertia = SpatialInertia::from_com(
        1.0,
        [0.0; 3],
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
    );
    let body = RigidBody::new("screw_body", inertia, None, SpatialTransform::IDENTITY);
    model.add_body(body, Box::new(HelicalJoint::new([0.0, 0.0, 1.0], 0.1)));
    let q = vec![0.0];
    let qd = vec![0.0];
    let tau = vec![1.0]; // 1 N·m about z
    let qdd = aba(&model, &q, &qd, &tau);
    assert_eq!(qdd.len(), 1);
    assert!(
        qdd[0] > 0.0,
        "positive torque should give positive qdd, got {}",
        qdd[0]
    );
    // After one step dt=0.1: q ≈ qdd * dt^2 / 2 → translation = pitch * q ≈ nonzero
    let dt = 0.1;
    let q_new = q[0] + qd[0] * dt + 0.5 * qdd[0] * dt * dt;
    let x_new = HelicalJoint::new([0.0, 0.0, 1.0], 0.1).transform(&[q_new]);
    assert!(
        x_new.trans[2].abs() > 1e-6,
        "should have nonzero z-translation after helical motion"
    );
}
