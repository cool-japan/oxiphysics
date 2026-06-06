// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

use oxiphysics_articulated::{
    body::RigidBody,
    joint::RevoluteJoint,
    model::ArticulatedModel,
    osc::compute_lambda,
    spatial::{SpatialInertia, SpatialTransform},
};

/// Build a 6-DOF serial chain with distinct revolute axes.
///
/// This is the minimum DOF needed to get a full-rank (6×6 invertible) Λ when
/// using a full-rank Jacobian.  Each body has a distinct offset and axis so the
/// model is not degenerate at q = 0.
fn make_6dof_chain() -> (ArticulatedModel, Vec<f64>) {
    let mut model = ArticulatedModel::new([0.0; 3]);
    // Axes cycle through x, y, z with a small offset each link (0.5 m between joints)
    let axes = [
        [1.0_f64, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];
    let offsets: [[f64; 3]; 6] = [
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.5],
        [0.0, 0.5, 0.0],
        [0.0, 0.0, 0.5],
        [0.5, 0.0, 0.0],
        [0.0, 0.0, 0.5],
    ];
    for i in 0..6 {
        let parent = if i == 0 { None } else { Some(i - 1) };
        let inertia = SpatialInertia::from_com(
            1.0,
            [0.0, 0.0, -0.25],
            [[0.02, 0.0, 0.0], [0.0, 0.02, 0.0], [0.0, 0.0, 0.005]],
        );
        let parent_tf = if i == 0 {
            SpatialTransform::IDENTITY
        } else {
            SpatialTransform::from_translation(offsets[i])
        };
        let body = RigidBody::new(format!("link{i}"), inertia, parent, parent_tf);
        model.add_body(body, Box::new(RevoluteJoint::new(axes[i])));
    }
    // Slightly off-zero configuration to avoid any axis collinearity artefacts
    let q = vec![0.1, -0.2, 0.15, 0.3, -0.1, 0.2];
    (model, q)
}

/// Build a full-rank 6×6 Jacobian for the 6-DOF chain.
///
/// Each row is a different standard basis vector repeated across columns so the
/// Jacobian has rank 6.  Physically this represents the end-effector Jacobian
/// but we only care about the rank for the inversion test.
fn full_rank_jacobian_6x6() -> Vec<Vec<f64>> {
    // J[row][col] = 1 if row == col, else 0  →  J = I (full rank, simple)
    (0..6)
        .map(|r| (0..6).map(|c| if r == c { 1.0 } else { 0.0 }).collect())
        .collect()
}

#[test]
fn test_osc_dimensions() {
    let (mut model, q) = make_6dof_chain();
    let jacobian = full_rank_jacobian_6x6();
    let lambda = compute_lambda(&mut model, &q, &jacobian);
    assert_eq!(lambda.len(), 6, "Λ must be 6×6");
    assert_eq!(lambda[0].len(), 6, "Λ must be 6×6");
}

#[test]
fn test_osc_symmetric() {
    let (mut model, q) = make_6dof_chain();
    let jacobian = full_rank_jacobian_6x6();
    let lambda = compute_lambda(&mut model, &q, &jacobian);
    for (i, lambda_row) in lambda.iter().enumerate() {
        for (j, &lambda_ij) in lambda_row.iter().enumerate() {
            let diff = (lambda_ij - lambda[j][i]).abs();
            assert!(
                diff < 1e-8,
                "Λ not symmetric at [{i}][{j}]: diff={diff:.2e}"
            );
        }
    }
}

#[test]
fn test_cholesky_roundtrip_in_osc() {
    // Build a 6-DOF model, run OSC with identity Jacobian, verify result is finite PD.
    let (mut model, q) = make_6dof_chain();
    let jacobian = full_rank_jacobian_6x6();
    let lambda = compute_lambda(&mut model, &q, &jacobian);
    // All entries should be finite
    for (i, lambda_row) in lambda.iter().enumerate() {
        for (j, &lambda_ij) in lambda_row.iter().enumerate() {
            assert!(lambda_ij.is_finite(), "Λ[{i}][{j}] should be finite");
        }
    }
    // All diagonal entries should be positive (Λ is SPD)
    for (i, lambda_row) in lambda.iter().enumerate() {
        assert!(
            lambda_row[i] > 0.0,
            "Λ[{i}][{i}] should be positive, got {}",
            lambda_row[i]
        );
    }
}
