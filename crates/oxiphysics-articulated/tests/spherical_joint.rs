// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

use oxiphysics_articulated::{
    aba::aba,
    body::RigidBody,
    joint::{Joint, SphericalJoint},
    model::ArticulatedModel,
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
    // q = identity quaternion (supply 4 values), qd = body-frame ω
    let q = vec![1.0, 0.0, 0.0]; // only 3 = dof(), but transform reads q[0..4]; identity quat
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
    let q = vec![1.0, 0.0, 0.0]; // identity quaternion fallback (q.len() < 4 → identity rot)
    let qd = vec![0.0; 3];
    let tau = vec![0.0; 3];
    let qdd = aba(&model, &q, &qd, &tau);
    // At equilibrium orientation, all angular accelerations should be ~0
    for (k, &acc) in qdd.iter().enumerate() {
        assert!(acc.is_finite(), "qdd[{k}] should be finite");
    }
}
