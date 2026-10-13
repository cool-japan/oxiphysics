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

    // q layout (nq = 7): [tx,ty,tz, qw,qx,qy,qz]; identity pose.
    let q = vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0];
    // Velocity (nv = 6): body twist [v; ω]. With the identity inertia about the
    // COM at the origin, ω × (I·ω) = 0 and ω × v is cancelled by the bias
    // term, so the body moves with constant twist.
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
    let q = vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0];
    let qd = vec![0.0; 6];
    let tau = vec![0.0; 6];
    let qdd = aba(&model, &q, &qd, &tau);
    assert_eq!(qdd.len(), 6);
    // DOFs 0,1,2 are the body-frame linear acceleration. A free body with
    // identity orientation under gravity [0,0,-g] falls: qdd_z = -g.
    let az = qdd[2];
    assert!(
        (az + g).abs() < 1e-8,
        "z acceleration should be -g={}, got {az:.6}",
        -g
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

    // total_dof = 7 (6 free + 1 revolute), total_nq = 8 (7 free + 1 revolute)
    assert_eq!(model.total_dof(), 7);
    assert_eq!(model.total_nq(), 8);
    let q = vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.3];
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

// ─── Orientation and SE(3) integration ───────────────────────────────────────

#[test]
fn free_floating_nq_and_neutral() {
    assert_eq!(FreeFloatingJoint.nq(), 7);
    assert_eq!(
        FreeFloatingJoint.neutral(),
        vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0]
    );
}

#[test]
fn free_floating_rotated_body_falls_in_body_frame() {
    // Orientation: +90° about x. World gravity −z expressed in the body frame
    // is Rᵀ·[0,0,−g] = [0,−g,0], so the body-frame linear acceleration is
    // [0, −g, 0] and the angular acceleration is zero.
    let g = 9.81_f64;
    let mut model = ArticulatedModel::new([0.0, 0.0, -g]);
    let inertia = SpatialInertia::from_com(
        1.7,
        [0.0; 3],
        [[0.4, 0.0, 0.0], [0.0, 0.3, 0.0], [0.0, 0.0, 0.2]],
    );
    model.add_body(
        RigidBody::new("free", inertia, None, SpatialTransform::IDENTITY),
        Box::new(FreeFloatingJoint),
    );
    let h = std::f64::consts::FRAC_PI_4;
    let q = vec![1.0, 2.0, 3.0, h.cos(), h.sin(), 0.0, 0.0];
    let qdd = aba(&model, &q, &[0.0; 6], &[0.0; 6]);
    let expected = [0.0, -g, 0.0, 0.0, 0.0, 0.0];
    for k in 0..6 {
        assert!(
            (qdd[k] - expected[k]).abs() < 1e-12,
            "qdd[{k}]={} expected {}",
            qdd[k],
            expected[k]
        );
    }
}

#[test]
fn free_floating_integrate_difference_round_trip() {
    let j = FreeFloatingJoint;
    let n = (0.7_f64 * 0.7 + 0.2 * 0.2 + 0.5 * 0.5 + 0.4 * 0.4).sqrt();
    let q0 = [0.3, -0.4, 1.2, 0.7 / n, 0.2 / n, -0.5 / n, 0.4 / n];
    let v = [0.5, -1.0, 0.25, 0.9, -0.3, 1.7];
    let dt = 0.8;
    let mut q1 = [0.0; 7];
    j.integrate(&q0, &v, dt, &mut q1);
    let mut back = [0.0; 6];
    j.difference(&q0, &q1, &mut back);
    for k in 0..6 {
        assert!(
            (back[k] - dt * v[k]).abs() < 1e-12,
            "difference[{k}]={} expected {}",
            back[k],
            dt * v[k]
        );
    }
}

#[test]
fn free_floating_integrate_matches_constant_twist_flow() {
    // Integrate the same constant body twist in many small steps and in one
    // step: the SE(3) exponential is exact, so both must agree.
    let j = FreeFloatingJoint;
    let q0 = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0];
    let v = [1.0, 0.0, 0.0, 0.0, 0.0, 2.0]; // forward + yaw: a circular arc
    let mut q_one = [0.0; 7];
    j.integrate(&q0, &v, 1.0, &mut q_one);
    let mut q_many = q0;
    for _ in 0..1000 {
        let mut next = [0.0; 7];
        j.integrate(&q_many, &v, 1e-3, &mut next);
        q_many = next;
    }
    for k in 0..7 {
        assert!(
            (q_one[k] - q_many[k]).abs() < 1e-10,
            "component {k}: one-step {} vs many-step {}",
            q_one[k],
            q_many[k]
        );
    }
    // Circular arc of radius v/ω = 0.5 swept through 2 rad.
    let (r, th) = (0.5_f64, 2.0_f64);
    assert!((q_one[0] - r * th.sin()).abs() < 1e-12);
    assert!((q_one[1] - r * (1.0 - th.cos())).abs() < 1e-12);
}

#[test]
fn free_floating_rotated_round_trip_with_arm() {
    let g = 9.81_f64;
    let mut model = ArticulatedModel::new([0.0, 0.0, -g]);
    model.add_body(
        RigidBody::new(
            "base",
            SpatialInertia::from_com(
                2.0,
                [0.1, 0.0, 0.05],
                [[0.5, 0.02, 0.0], [0.02, 0.4, 0.0], [0.0, 0.0, 0.3]],
            ),
            None,
            SpatialTransform::IDENTITY,
        ),
        Box::new(FreeFloatingJoint),
    );
    model.add_body(
        RigidBody::new(
            "arm",
            SpatialInertia::from_com(
                1.0,
                [0.0, 0.0, -0.5],
                [[0.1, 0.0, 0.0], [0.0, 0.1, 0.0], [0.0, 0.0, 0.01]],
            ),
            Some(0),
            SpatialTransform::from_translation([0.2, 0.0, 0.0]),
        ),
        Box::new(RevoluteJoint::new([1.0, 0.0, 0.0])),
    );
    let h = 0.4_f64;
    let q = vec![
        0.5,
        -0.2,
        1.0,
        h.cos(),
        0.0,
        h.sin() * 0.6,
        h.sin() * 0.8,
        0.3,
    ];
    let qd = vec![0.3, -0.1, 0.2, 0.4, -0.6, 0.2, 0.5];
    let qdd_target = vec![0.1, 0.2, -0.3, 0.05, -0.05, 0.1, 1.5];
    let tau = rnea(&model, &q, &qd, &qdd_target);
    let qdd = aba(&model, &q, &qd, &tau);
    for k in 0..7 {
        assert!(
            (qdd[k] - qdd_target[k]).abs() < 1e-9,
            "qdd[{k}]={} target={}",
            qdd[k],
            qdd_target[k]
        );
    }
}
