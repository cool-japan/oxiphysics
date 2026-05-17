// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

#![allow(clippy::needless_range_loop)]

use oxiphysics_articulated::{
    body::RigidBody,
    crba::compute_mass_matrix_crba,
    joint::RevoluteJoint,
    model::ArticulatedModel,
    rnea::rnea,
    spatial::{SpatialInertia, SpatialTransform},
};

fn make_2link_arm() -> (ArticulatedModel, Vec<f64>) {
    let mut model = ArticulatedModel::new([0.0; 3]);
    let i0 = SpatialInertia::from_com(
        2.0,
        [0.0, 0.0, -0.5],
        [[0.1, 0.0, 0.0], [0.0, 0.1, 0.0], [0.0, 0.0, 0.05]],
    );
    let b0 = RigidBody::new("link0", i0, None, SpatialTransform::IDENTITY);
    model.add_body(b0, Box::new(RevoluteJoint::new([1.0, 0.0, 0.0])));
    let i1 = SpatialInertia::from_com(
        1.0,
        [0.0, 0.0, -0.4],
        [[0.05, 0.0, 0.0], [0.0, 0.05, 0.0], [0.0, 0.0, 0.01]],
    );
    let b1 = RigidBody::new(
        "link1",
        i1,
        Some(0),
        SpatialTransform::from_translation([0.0, 0.0, -1.0]),
    );
    model.add_body(b1, Box::new(RevoluteJoint::new([0.0, 1.0, 0.0])));
    let q = vec![0.3, -0.5];
    (model, q)
}

#[test]
fn test_crba_vs_rnea_2link_arm() {
    let (mut model, q) = make_2link_arm();
    let m_crba = compute_mass_matrix_crba(&mut model, &q);
    let n = 2;

    // Compute M via n RNEA calls: M[:, k] = RNEA(q, 0, e_k) (no gravity, no velocity)
    let mut m_rnea = vec![vec![0.0f64; n]; n];
    let qd = vec![0.0; n];
    for k in 0..n {
        let mut qdd = vec![0.0; n];
        qdd[k] = 1.0;
        let tau = rnea(&model, &q, &qd, &qdd);
        for i in 0..n {
            m_rnea[i][k] = tau[i];
        }
    }
    // Zero gravity for RNEA (gravity torques cancel out in the difference)
    // Note: gravity was set to zero in make_2link_arm, so tau = M · qdd
    for i in 0..n {
        for j in 0..n {
            let diff = (m_crba[i][j] - m_rnea[i][j]).abs();
            assert!(
                diff < 1e-8,
                "M_crba[{i}][{j}]={:.6} M_rnea[{i}][{j}]={:.6} diff={:.2e}",
                m_crba[i][j],
                m_rnea[i][j],
                diff
            );
        }
    }
}

#[test]
fn test_crba_symmetric() {
    let (mut model, q) = make_2link_arm();
    let m = compute_mass_matrix_crba(&mut model, &q);
    let n = m.len();
    for i in 0..n {
        for j in 0..n {
            let diff = (m[i][j] - m[j][i]).abs();
            assert!(
                diff < 1e-10,
                "M not symmetric: M[{i}][{j}]={} M[{j}][{i}]={} diff={:.2e}",
                m[i][j],
                m[j][i],
                diff
            );
        }
    }
}

#[test]
fn test_crba_positive_definite() {
    let (mut model, q) = make_2link_arm();
    let m = compute_mass_matrix_crba(&mut model, &q);
    // Check PD via Cholesky: all diagonal entries of L must be > 0
    let n = m.len();
    let mut l = m.clone();
    for i in 0..n {
        for j in 0..=i {
            let mut sum = l[i][j];
            for k in 0..j {
                sum -= l[i][k] * l[j][k];
            }
            if i == j {
                assert!(
                    sum > 1e-14,
                    "M is not positive-definite at pivot ({i},{i}): sum={sum:.2e}"
                );
                l[i][j] = sum.sqrt();
            } else {
                l[i][j] = sum / l[j][j];
            }
        }
    }
}
