// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

use oxiphysics_articulated::{
    aba::aba,
    body::RigidBody,
    joint::{Joint, RevoluteJoint, UniversalJoint},
    model::ArticulatedModel,
    spatial::{SpatialInertia, SpatialTransform},
};
use std::f64::consts::PI;

fn assert_mat3_approx(a: [[f64; 3]; 3], b: [[f64; 3]; 3], tol: f64, msg: &str) {
    for i in 0..3 {
        for j in 0..3 {
            let diff = (a[i][j] - b[i][j]).abs();
            assert!(
                diff < tol,
                "{}: a[{}][{}]={:.2e} b[{}][{}]={:.2e} diff={:.2e}",
                msg,
                i,
                j,
                a[i][j],
                i,
                j,
                b[i][j],
                diff
            );
        }
    }
}

#[test]
fn test_universal_joint_dof_and_transform_at_zero() {
    let j = UniversalJoint::new([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
    assert_eq!(j.dof(), 2);
    let x = j.transform(&[0.0, 0.0]);
    assert_mat3_approx(
        x.rot,
        SpatialTransform::IDENTITY.rot,
        1e-14,
        "zero-q transform should be identity",
    );
    for component in &x.trans {
        assert!(component.abs() < 1e-14, "translation should be zero");
    }
}

#[test]
fn test_universal_locked_first_axis_reduces_to_revolute() {
    let universal = UniversalJoint::new([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
    let revolute = RevoluteJoint::new([0.0, 1.0, 0.0]);
    for k in 0..8 {
        let angle = k as f64 * PI / 4.0;
        let xu = universal.transform(&[0.0, angle]);
        let xr = revolute.transform(&[angle]);
        assert_mat3_approx(xu.rot, xr.rot, 1e-12, &format!("angle={:.3}", angle));
    }
}

#[test]
fn test_universal_motion_subspace_q_dependence() {
    let j = UniversalJoint::new([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
    let s0 = j.motion_subspace(); // at q=[0,0]
    let sq = j.motion_subspace_at(&[0.3, 0.7]); // at non-zero q
    // column 0 should differ between q-independent and q-dependent
    let col0_s0 = s0[0].angular;
    let col0_sq = sq[0].angular;
    let diff = ((col0_s0[0] - col0_sq[0]).powi(2)
        + (col0_s0[1] - col0_sq[1]).powi(2)
        + (col0_s0[2] - col0_sq[2]).powi(2))
    .sqrt();
    assert!(
        diff > 0.05,
        "motion_subspace_at should differ from motion_subspace() at non-zero q, got diff={:.4}",
        diff
    );
}

#[test]
fn test_universal_dynamics_gravity_sign() {
    // 1-body universal pendulum, mass 1 kg at (0, 0, -1) below the joint.
    let g = 9.81_f64;
    let mut model = ArticulatedModel::new([0.0, 0.0, -g]);
    let inertia = SpatialInertia::from_com(1.0, [0.0, 0.0, -1.0], [[0.0; 3]; 3]);
    let body = RigidBody::new("link1", inertia, None, SpatialTransform::IDENTITY);
    model.add_body(
        body,
        Box::new(UniversalJoint::new([1.0, 0.0, 0.0], [0.0, 1.0, 0.0])),
    );
    let q = vec![0.0, 0.0];
    let qd = vec![0.0, 0.0];
    let tau = vec![0.0, 0.0];
    let qdd = aba(&model, &q, &qd, &tau);
    assert_eq!(qdd.len(), 2);
    // At q=[0,0], COM is directly below. With mass at (0,0,-1) and gravity in -z,
    // the torques about both axes are zero (symmetric configuration). So qdd~[0,0].
    for (k, &val) in qdd.iter().enumerate() {
        assert!(val.is_finite(), "qdd[{}] must be finite", k);
    }
    // Now test with first joint rotated so mass is displaced:
    let q2 = vec![PI / 6.0, 0.0]; // 30 degrees about x
    let qdd2 = aba(&model, &q2, &qd, &tau);
    assert!(qdd2[0].is_finite() && qdd2[1].is_finite());
    // With mass displaced, at least one qdd should be non-zero
    let any_nonzero = qdd2[0].abs() > 1e-6 || qdd2[1].abs() > 1e-6;
    assert!(
        any_nonzero,
        "qdd should be non-zero when joint is at non-equilibrium, got {:?}",
        qdd2
    );
}

#[test]
fn test_universal_rnea_aba_round_trip() {
    use oxiphysics_articulated::rnea::rnea;
    let g = 9.81_f64;
    let mut model = ArticulatedModel::new([0.0, 0.0, -g]);
    let inertia = SpatialInertia::from_com(
        2.0,
        [0.0, 0.0, -0.5],
        [[0.1, 0.0, 0.0], [0.0, 0.1, 0.0], [0.0, 0.0, 0.05]],
    );
    let body = RigidBody::new("link1", inertia, None, SpatialTransform::IDENTITY);
    model.add_body(
        body,
        Box::new(UniversalJoint::new([1.0, 0.0, 0.0], [0.0, 1.0, 0.0])),
    );
    let q = vec![0.2, -0.3];
    let qd = vec![0.5, -0.1];
    let qdd_target = vec![1.0, -2.0];
    // RNEA: given qdd_target, compute τ
    let tau = rnea(&model, &q, &qd, &qdd_target);
    // ABA: given τ, recover qdd
    let qdd_recovered = aba(&model, &q, &qd, &tau);
    assert_eq!(qdd_recovered.len(), 2);
    for (k, (&target, &recovered)) in qdd_target.iter().zip(qdd_recovered.iter()).enumerate() {
        let diff = (recovered - target).abs();
        assert!(
            diff < 1e-8,
            "qdd[{}] round-trip: target={:.6} recovered={:.6} diff={:.2e}",
            k,
            target,
            recovered,
            diff
        );
    }
}
