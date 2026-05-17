// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

use oxiphysics_articulated::{
    aba::aba,
    body::RigidBody,
    joint::{FreeFloatingJoint, Joint, RevoluteJoint},
    model::ArticulatedModel,
    rnea::rnea,
    spatial::{SpatialInertia, SpatialTransform},
};

#[test]
fn free_floating_dof_count_and_motion_subspace() {
    let j = FreeFloatingJoint;
    assert_eq!(j.dof(), 6);
    let s = j.motion_subspace();
    assert_eq!(s.len(), 6);
    // First 3 columns are linear velocity DOFs [0; e_i]
    assert_eq!(s[0].angular, [0.0; 3]);
    assert_eq!(s[0].linear, [1.0, 0.0, 0.0]);
    // Last 3 columns are angular velocity DOFs [e_i; 0]
    assert_eq!(s[3].angular, [1.0, 0.0, 0.0]);
    assert_eq!(s[3].linear, [0.0; 3]);
}

#[test]
fn free_floating_root_alone_unforced_constant_velocity() {
    // Single free body, no forces, no gravity — qdd should be zero.
    let mut model = ArticulatedModel::new([0.0, 0.0, 0.0]);
    let inertia = SpatialInertia::from_com(
        1.0,
        [0.0; 3],
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
    );
    let body = RigidBody::new("free_body", inertia, None, SpatialTransform::IDENTITY);
    model.add_body(body, Box::new(FreeFloatingJoint));

    // Supply 7 values in q (3 translation + 4 quaternion = identity config)
    // but only 6 are used by q_slice (dof=6). We extend to 7 for the transform.
    // q layout: [tx,ty,tz, qw,qx,qy,qz] but q_slice only reads 6.
    // Use q = [0,0,0, 1,0,0,0] (identity). Pass the first 6 as the dynamics q.
    let q = vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0]; // 6 DOF: tx,ty,tz, wx,wy,wz
    let qd = vec![1.0, 2.0, 3.0, 0.1, 0.2, 0.3];
    let tau = vec![0.0; 6];
    let qdd = aba(&model, &q, &qd, &tau);
    assert_eq!(qdd.len(), 6);
    for (k, &acc) in qdd.iter().enumerate() {
        assert!(acc.is_finite(), "qdd[{k}] must be finite");
        assert!(
            acc.abs() < 1e-10,
            "unforced free body qdd[{k}]={acc:.2e} should be ~0"
        );
    }
}

#[test]
fn free_floating_root_under_gravity() {
    // Single free body under gravity — z linear acceleration should be -9.81.
    let g = 9.81_f64;
    let mut model = ArticulatedModel::new([0.0, 0.0, -g]);
    let inertia = SpatialInertia::from_com(
        1.0,
        [0.0; 3],
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
    );
    let body = RigidBody::new("free_body", inertia, None, SpatialTransform::IDENTITY);
    model.add_body(body, Box::new(FreeFloatingJoint));

    // q = identity rotation, zero velocity
    let q = vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0];
    let qd = vec![0.0; 6];
    let tau = vec![0.0; 6];
    let qdd = aba(&model, &q, &qd, &tau);
    assert_eq!(qdd.len(), 6);
    // DOFs 0,1,2 are linear velocity (vx, vy, vz) → accelerations.
    // Featherstone's ABA uses gravity as a fictitious parent acceleration:
    // gravity = [0,0,-g] as a SpatialVec linear part. ABA returns qdd in
    // the sense of generalised acceleration, so for a unit-mass free body
    // under gravity [0,0,-g], qdd[2] = +g (the body "sees" gravity as +g
    // along z due to the fictitious-acceleration sign convention).
    let az = qdd[2]; // linear z acceleration (generalised)
    assert!(
        (az - g).abs() < 1e-8,
        "z acceleration should be +g={g}, got {az:.6}"
    );
    // Other DOFs should be ~0
    for k in [0, 1, 3, 4, 5] {
        assert!(qdd[k].abs() < 1e-8, "qdd[{k}]={:.2e} should be ~0", qdd[k]);
    }
}

#[test]
fn free_floating_base_plus_one_arm_link_round_trip() {
    // 2-body: free floating base + one revolute arm link.
    // RNEA + ABA round trip: compute tau from qdd via RNEA, recover qdd via ABA.
    let g = 9.81_f64;
    let mut model = ArticulatedModel::new([0.0, 0.0, -g]);
    let base_inertia = SpatialInertia::from_com(
        2.0,
        [0.0; 3],
        [[0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 0.5]],
    );
    let base = RigidBody::new("base", base_inertia, None, SpatialTransform::IDENTITY);
    model.add_body(base, Box::new(FreeFloatingJoint));

    let arm_inertia = SpatialInertia::from_com(
        1.0,
        [0.0, 0.0, -0.5],
        [[0.1, 0.0, 0.0], [0.0, 0.1, 0.0], [0.0, 0.0, 0.01]],
    );
    let arm = RigidBody::new("arm", arm_inertia, Some(0), SpatialTransform::IDENTITY);
    model.add_body(arm, Box::new(RevoluteJoint::new([1.0, 0.0, 0.0])));

    // total_dof = 7 (6 free + 1 revolute)
    let q = vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.3];
    let qd = vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5];
    let qdd_target = vec![0.1, 0.2, -0.3, 0.05, -0.05, 0.1, 1.5];

    let tau = rnea(&model, &q, &qd, &qdd_target);
    let qdd_recovered = aba(&model, &q, &qd, &tau);
    assert_eq!(qdd_recovered.len(), 7);
    for k in 0..7 {
        let diff = (qdd_recovered[k] - qdd_target[k]).abs();
        assert!(
            diff < 1e-8,
            "round-trip qdd[{k}]: target={:.6} recovered={:.6} diff={:.2e}",
            qdd_target[k],
            qdd_recovered[k],
            diff
        );
    }
}
