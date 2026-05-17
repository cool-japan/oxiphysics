// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

use oxiphysics_articulated::{
    body::RigidBody,
    centroidal::compute_cmm,
    joint::RevoluteJoint,
    model::ArticulatedModel,
    spatial::{SpatialInertia, SpatialTransform},
};

#[test]
fn test_cmm_dimensions() {
    let mut model = ArticulatedModel::new([0.0; 3]);
    let inertia = SpatialInertia::from_com(
        1.0,
        [0.0; 3],
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
    );
    for i in 0..4 {
        let parent = if i == 0 { None } else { Some(i - 1) };
        let body = RigidBody::new(
            format!("link{i}"),
            inertia,
            parent,
            SpatialTransform::IDENTITY,
        );
        model.add_body(body, Box::new(RevoluteJoint::new([0.0, 0.0, 1.0])));
    }
    let q = vec![0.0; 4];
    let a_g = compute_cmm(&mut model, &q);
    assert_eq!(a_g.len(), 6, "A_G must have 6 rows");
    assert_eq!(a_g[0].len(), 4, "A_G must have n columns");
}

#[test]
fn test_cmm_nonzero_for_revolute_chain() {
    // A 2-link revolute arm should have non-zero CMM entries
    let mut model = ArticulatedModel::new([0.0; 3]);
    let i0 = SpatialInertia::from_com(
        1.0,
        [0.0, 0.0, -0.5],
        [[0.1, 0.0, 0.0], [0.0, 0.1, 0.0], [0.0, 0.0, 0.05]],
    );
    let b0 = RigidBody::new("b0", i0, None, SpatialTransform::IDENTITY);
    model.add_body(b0, Box::new(RevoluteJoint::new([1.0, 0.0, 0.0])));
    let i1 = SpatialInertia::from_com(
        1.0,
        [0.0, 0.0, -0.5],
        [[0.1, 0.0, 0.0], [0.0, 0.1, 0.0], [0.0, 0.0, 0.05]],
    );
    let b1 = RigidBody::new(
        "b1",
        i1,
        Some(0),
        SpatialTransform::from_translation([0.0, 0.0, -1.0]),
    );
    model.add_body(b1, Box::new(RevoluteJoint::new([0.0, 1.0, 0.0])));
    let q = vec![0.2, -0.3];
    let a_g = compute_cmm(&mut model, &q);
    // At least some entries should be non-zero
    let any_nonzero = a_g.iter().any(|row| row.iter().any(|&v| v.abs() > 1e-10));
    assert!(
        any_nonzero,
        "CMM should have non-zero entries for a physical system"
    );
}
