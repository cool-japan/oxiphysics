// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

use oxiphysics_articulated::{
    aba::aba,
    body::RigidBody,
    joint::{Joint, RevoluteJoint, SphericalJoint},
    model::ArticulatedModel,
    rnea::rnea,
    rnea_derivatives::rnea_derivatives,
    spatial::{SpatialInertia, SpatialTransform},
};
use std::f64::consts::PI;

#[test]
fn test_spherical_dof_count() {
    assert_eq!(SphericalJoint.dof(), 3);
}

#[test]
fn test_spherical_identity_quaternion_yields_identity_transform() {
    let x = SphericalJoint.transform(&[1.0, 0.0, 0.0, 0.0]);
    let id = SpatialTransform::IDENTITY;
    for i in 0..3 {
        for j in 0..3 {
            assert!(
                (x.rot[i][j] - id.rot[i][j]).abs() < 1e-15,
                "rot[{i}][{j}]: got {:.2e}",
                x.rot[i][j] - id.rot[i][j]
            );
        }
    }
}

#[test]
fn test_spherical_90deg_about_x_quaternion() {
    // q = [cos(π/4), sin(π/4), 0, 0] = 90 deg about x
    let half = PI / 4.0;
    let q = [half.cos(), half.sin(), 0.0, 0.0];
    let x = SphericalJoint.transform(&q);
    // Rotation about x by 90 deg: e_y → e_z, e_z → -e_y
    assert!((x.rot[0][0] - 1.0).abs() < 1e-12, "R[0][0] should be 1");
    assert!(
        x.rot[1][2].abs() < 1e-12 || (x.rot[1][2] - (-1.0)).abs() < 1e-12,
        "unexpected R[1][2]={}",
        x.rot[1][2]
    );
}

#[test]
fn test_spherical_motion_subspace_constant() {
    let s = SphericalJoint.motion_subspace();
    assert_eq!(s.len(), 3);
    // Column 0: angular=[1,0,0], linear=[0,0,0]
    assert_eq!(s[0].angular, [1.0, 0.0, 0.0]);
    assert_eq!(s[0].linear, [0.0; 3]);
    // Column 1: angular=[0,1,0], linear=[0,0,0]
    assert_eq!(s[1].angular, [0.0, 1.0, 0.0]);
    assert_eq!(s[1].linear, [0.0; 3]);
    // Column 2: angular=[0,0,1], linear=[0,0,0]
    assert_eq!(s[2].angular, [0.0, 0.0, 1.0]);
    assert_eq!(s[2].linear, [0.0; 3]);
}

#[test]
fn test_spherical_dynamics_no_gravity_constant_velocity() {
    let mut model = ArticulatedModel::new([0.0; 3]);
    let inertia = SpatialInertia::from_com(
        1.0,
        [0.0; 3],
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
    );
    let body = RigidBody::new("sphere", inertia, None, SpatialTransform::IDENTITY);
    model.add_body(body, Box::new(SphericalJoint));
    // q = identity quaternion (nq = 4), qd = body-frame ω (nv = 3)
    assert_eq!(model.total_nq(), 4);
    assert_eq!(model.total_dof(), 3);
    let q = vec![1.0, 0.0, 0.0, 0.0];
    let qd = vec![1.0, 0.0, 0.0];
    let tau = vec![0.0; 3];
    let qdd = aba(&model, &q, &qd, &tau);
    assert_eq!(qdd.len(), 3);
    for (k, &acc) in qdd.iter().enumerate() {
        assert!(
            acc.abs() < 1e-10,
            "unforced sphere qdd[{k}]={acc:.2e} should be ~0"
        );
    }
}

#[test]
fn test_spherical_pendulum_equilibrium() {
    // Body hanging straight down: COM at (0, 0, -0.5), gravity in -z.
    // With identity quaternion, COM is directly below → zero net torque → qdd ≈ 0.
    let g = 9.81_f64;
    let mut model = ArticulatedModel::new([0.0, 0.0, -g]);
    // Use non-degenerate inertia to avoid singular mass matrix in ABA.
    let inertia = SpatialInertia::from_com(
        1.0,
        [0.0, 0.0, -0.5],
        [[0.1, 0.0, 0.0], [0.0, 0.1, 0.0], [0.0, 0.0, 0.1]],
    );
    let body = RigidBody::new("pendulum", inertia, None, SpatialTransform::IDENTITY);
    model.add_body(body, Box::new(SphericalJoint));
    let q = vec![1.0, 0.0, 0.0, 0.0]; // identity quaternion
    let qd = vec![0.0; 3];
    let tau = vec![0.0; 3];
    let qdd = aba(&model, &q, &qd, &tau);
    // At equilibrium orientation, all angular accelerations should be ~0
    for (k, &acc) in qdd.iter().enumerate() {
        assert!(acc.abs() < 1e-12, "qdd[{k}]={acc:.3e} should be ~0");
    }
}

// ─── Quaternion configuration actually reaches the algorithms ────────────────

fn spherical_pendulum(m: f64, l: f64, g: f64) -> ArticulatedModel {
    let mut model = ArticulatedModel::new([0.0, 0.0, -g]);
    let inertia = SpatialInertia::from_com(m, [0.0, 0.0, -l], [[0.0; 3]; 3]);
    model.add_body(
        RigidBody::new("bob", inertia, None, SpatialTransform::IDENTITY),
        Box::new(SphericalJoint),
    );
    model
}

#[test]
fn test_spherical_nq_and_neutral() {
    assert_eq!(SphericalJoint.nq(), 4);
    assert_eq!(SphericalJoint.neutral(), vec![1.0, 0.0, 0.0, 0.0]);
    let model = spherical_pendulum(1.0, 1.0, 9.81);
    assert_eq!(model.total_nq(), 4);
    assert_eq!(model.total_dof(), 3);
    assert_eq!(model.neutral_configuration(), vec![1.0, 0.0, 0.0, 0.0]);
}

#[test]
fn test_spherical_pendulum_tilted_90deg_gravity_torque() {
    // Tilt 90° about x: the bob at local [0,0,−l] moves to world [0,+l,0].
    // Gravity −z then needs a holding torque of +m·g·l about x, exactly like a
    // revolute-x pendulum at q = π/2.
    let (m, l, g) = (1.3, 0.8, 9.81);
    let model = spherical_pendulum(m, l, g);
    let half = PI / 4.0;
    let q = [half.cos(), half.sin(), 0.0, 0.0];
    let tau = rnea(&model, &q, &[0.0; 3], &[0.0; 3]);
    assert!(
        (tau[0] - m * g * l).abs() < 1e-12,
        "tau_x={} expected {}",
        tau[0],
        m * g * l
    );
    assert!(tau[1].abs() < 1e-12 && tau[2].abs() < 1e-12, "tau={tau:?}");

    let mut revolute = ArticulatedModel::new([0.0, 0.0, -g]);
    revolute.add_body(
        RigidBody::new(
            "bob",
            SpatialInertia::from_com(m, [0.0, 0.0, -l], [[0.0; 3]; 3]),
            None,
            SpatialTransform::IDENTITY,
        ),
        Box::new(RevoluteJoint::new([1.0, 0.0, 0.0])),
    );
    let tau_rev = rnea(&revolute, &[PI / 2.0], &[0.0], &[0.0]);
    assert!((tau[0] - tau_rev[0]).abs() < 1e-12);

    // Released from rest, it swings back: θ̈_x = −g/l.
    let qdd = aba(&model, &q, &[0.0; 3], &[0.0; 3]);
    assert!((qdd[0] + g / l).abs() < 1e-10, "qdd={qdd:?}");
}

#[test]
fn test_spherical_rnea_aba_round_trip_general_state() {
    let g = 9.81;
    let mut model = ArticulatedModel::new([0.0, -g, 0.0]);
    let inertia = SpatialInertia::from_com(
        2.0,
        [0.1, -0.4, 0.2],
        [[0.2, 0.01, 0.0], [0.01, 0.15, -0.02], [0.0, -0.02, 0.1]],
    );
    model.add_body(
        RigidBody::new("ball", inertia, None, SpatialTransform::IDENTITY),
        Box::new(SphericalJoint),
    );
    let n = (0.3_f64 * 0.3 + 0.5 * 0.5 + 0.2 * 0.2 + 0.8 * 0.8).sqrt();
    let q = [0.8 / n, 0.3 / n, -0.5 / n, 0.2 / n];
    let qd = [0.4, -1.2, 0.7];
    let qdd_target = [1.0, -0.5, 2.0];
    let tau = rnea(&model, &q, &qd, &qdd_target);
    let qdd = aba(&model, &q, &qd, &tau);
    for k in 0..3 {
        assert!(
            (qdd[k] - qdd_target[k]).abs() < 1e-9,
            "qdd[{k}]={} target={}",
            qdd[k],
            qdd_target[k]
        );
    }
}

#[test]
fn test_spherical_integrate_difference_round_trip() {
    let j = SphericalJoint;
    let q0 = {
        let n = (0.9_f64 * 0.9 + 0.1 * 0.1 + 0.3 * 0.3 + 0.2 * 0.2).sqrt();
        [0.9 / n, 0.1 / n, -0.3 / n, 0.2 / n]
    };
    let v = [0.7, -1.4, 0.4];
    let mut q1 = [0.0; 4];
    j.integrate(&q0, &v, 0.5, &mut q1);
    let norm = q1.iter().map(|x| x * x).sum::<f64>().sqrt();
    assert!(
        (norm - 1.0).abs() < 1e-14,
        "integrated quaternion stays unit"
    );
    let mut back = [0.0; 3];
    j.difference(&q0, &q1, &mut back);
    for k in 0..3 {
        assert!(
            (back[k] - 0.5 * v[k]).abs() < 1e-12,
            "difference[{k}]={} expected {}",
            back[k],
            0.5 * v[k]
        );
    }
    // Body-frame semantics: R1 = R0 · exp(v·dt).
    let r0 = j.transform(&q0).rot;
    let r1 = j.transform(&q1).rot;
    let phi = [0.5 * v[0], 0.5 * v[1], 0.5 * v[2]];
    let theta = (phi[0] * phi[0] + phi[1] * phi[1] + phi[2] * phi[2]).sqrt();
    let r_delta = SpatialTransform::from_axis_angle(phi, theta).rot;
    for i in 0..3 {
        for c in 0..3 {
            let expected: f64 = (0..3).map(|k| r0[i][k] * r_delta[k][c]).sum();
            assert!((r1[i][c] - expected).abs() < 1e-12, "R1[{i}][{c}]");
        }
    }
}

#[test]
fn test_spherical_rnea_derivative_in_tangent_space() {
    // τ(R) = −m·c × (Rᵀ·g) for a static spherical pendulum, so along a body
    // tangent perturbation R·exp(δ): ∂τ/∂δ = −m·[c]×·[u]× with u = Rᵀ·g.
    let (m, g) = (1.5, 9.81);
    let c = [0.1, -0.2, -0.7];
    let mut model = ArticulatedModel::new([0.0, 0.0, -g]);
    model.add_body(
        RigidBody::new(
            "bob",
            SpatialInertia::from_com(m, c, [[0.0; 3]; 3]),
            None,
            SpatialTransform::IDENTITY,
        ),
        Box::new(SphericalJoint),
    );
    let half = 0.6_f64;
    let axis = [0.0, 0.6, 0.8];
    let q = [
        half.cos(),
        half.sin() * axis[0],
        half.sin() * axis[1],
        half.sin() * axis[2],
    ];
    let d = rnea_derivatives(&model, &q, &[0.0; 3], &[0.0; 3]);
    assert_eq!(d.dtau_dq.len(), 3);
    assert_eq!(d.dtau_dq[0].len(), 3);

    let r = SphericalJoint.transform(&q).rot;
    let g_world = [0.0, 0.0, -g];
    let u: Vec<f64> = (0..3)
        .map(|i| (0..3).map(|k| r[k][i] * g_world[k]).sum())
        .collect();
    let skew = |a: &[f64]| [[0.0, -a[2], a[1]], [a[2], 0.0, -a[0]], [-a[1], a[0], 0.0]];
    let sc = skew(&c);
    let su = skew(&u);
    for (i, sc_i) in sc.iter().enumerate() {
        for jcol in 0..3 {
            let expected: f64 = -m
                * sc_i
                    .iter()
                    .zip(su.iter())
                    .map(|(s_ik, su_k)| s_ik * su_k[jcol])
                    .sum::<f64>();
            assert!(
                (d.dtau_dq[i][jcol] - expected).abs() < 1e-6,
                "dtau_dq[{i}][{jcol}]={} expected {}",
                d.dtau_dq[i][jcol],
                expected
            );
        }
    }
}
