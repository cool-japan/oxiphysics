// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

use oxiphysics_articulated::{
    body::RigidBody,
    crba::compute_mass_matrix_crba,
    joint::{FreeFloatingJoint, RevoluteJoint, SphericalJoint, UniversalJoint},
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
    let (model, q) = make_2link_arm();
    let m_crba = compute_mass_matrix_crba(&model, &q);
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
    let (model, q) = make_2link_arm();
    let m = compute_mass_matrix_crba(&model, &q);
    let _n = m.len();
    for (i, m_row) in m.iter().enumerate() {
        for (j, &m_ij) in m_row.iter().enumerate() {
            let diff = (m_ij - m[j][i]).abs();
            assert!(
                diff < 1e-10,
                "M not symmetric: M[{i}][{j}]={} M[{j}][{i}]={} diff={:.2e}",
                m_ij,
                m[j][i],
                diff
            );
        }
    }
}

#[test]
fn test_crba_positive_definite() {
    let (model, q) = make_2link_arm();
    let m = compute_mass_matrix_crba(&model, &q);
    // Check PD via Cholesky: all diagonal entries of L must be > 0
    let n = m.len();
    let mut l = m.clone();
    for i in 0..n {
        for j in 0..=i {
            let mut sum = l[i][j];
            for (k, &l_jk) in l[j][..j].iter().enumerate() {
                sum -= l[i][k] * l_jk;
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

// ─── Multi-DoF joints: intra-joint off-diagonal entries ──────────────────────

/// Assert CRBA's `M(q)` equals the RNEA-derived matrix (`M[:, k] = RNEA(q, 0, e_k)`
/// with zero gravity) entry by entry.
fn assert_crba_matches_rnea(model: &ArticulatedModel, q: &[f64], tol: f64) {
    let n = model.total_dof();
    let m_crba = compute_mass_matrix_crba(model, q);
    let qd = vec![0.0; n];
    for k in 0..n {
        let mut qdd = vec![0.0; n];
        qdd[k] = 1.0;
        let col = rnea(model, q, &qd, &qdd);
        for i in 0..n {
            assert!(
                (m_crba[i][k] - col[i]).abs() < tol,
                "M_crba[{i}][{k}]={:.12} vs RNEA {:.12}",
                m_crba[i][k],
                col[i]
            );
        }
    }
}

fn full_inertia() -> [[f64; 3]; 3] {
    [[0.31, 0.04, -0.02], [0.04, 0.22, 0.05], [-0.02, 0.05, 0.17]]
}

#[test]
fn test_crba_universal_off_diagonal_block() {
    let mut model = ArticulatedModel::new([0.0; 3]);
    let inertia = SpatialInertia::from_com(1.4, [0.3, -0.2, -0.5], full_inertia());
    model.add_body(
        RigidBody::new("u", inertia, None, SpatialTransform::IDENTITY),
        Box::new(UniversalJoint::new([1.0, 0.0, 0.0], [0.0, 1.0, 0.0])),
    );
    let q = vec![0.4, -0.9];
    let m = compute_mass_matrix_crba(&model, &q);
    assert!(
        m[0][1].abs() > 1e-3,
        "universal joint must have a non-zero intra-joint coupling, got {}",
        m[0][1]
    );
    assert_crba_matches_rnea(&model, &q, 1e-10);
}

#[test]
fn test_crba_free_floating_off_diagonal_block() {
    let mut model = ArticulatedModel::new([0.0; 3]);
    let inertia = SpatialInertia::from_com(2.2, [0.25, 0.1, -0.4], full_inertia());
    model.add_body(
        RigidBody::new("base", inertia, None, SpatialTransform::IDENTITY),
        Box::new(FreeFloatingJoint),
    );
    // Non-identity orientation: unit quaternion about a skew axis.
    let (half, axis) = (0.35_f64, [0.48_f64, -0.6, 0.64]);
    let q = vec![
        0.3,
        -1.2,
        0.8,
        half.cos(),
        half.sin() * axis[0],
        half.sin() * axis[1],
        half.sin() * axis[2],
    ];
    let m = compute_mass_matrix_crba(&model, &q);
    // Linear–angular coupling m·[c]× is non-zero for an off-origin COM.
    assert!(
        m[0][4].abs() > 1e-3,
        "M[0][4]={} should be non-zero",
        m[0][4]
    );
    assert_crba_matches_rnea(&model, &q, 1e-10);
}

#[test]
fn test_crba_mixed_tree_with_spherical_and_free_base() {
    let mut model = ArticulatedModel::new([0.0; 3]);
    model.add_body(
        RigidBody::new(
            "base",
            SpatialInertia::from_com(3.0, [0.05, 0.0, 0.1], full_inertia()),
            None,
            SpatialTransform::IDENTITY,
        ),
        Box::new(FreeFloatingJoint),
    );
    model.add_body(
        RigidBody::new(
            "shoulder",
            SpatialInertia::from_com(1.1, [0.0, 0.0, -0.3], full_inertia()),
            Some(0),
            SpatialTransform::from_translation([0.2, 0.1, 0.0]),
        ),
        Box::new(SphericalJoint),
    );
    model.add_body(
        RigidBody::new(
            "wrist",
            SpatialInertia::from_com(0.6, [0.1, 0.05, -0.2], full_inertia()),
            Some(1),
            SpatialTransform::from_translation([0.0, 0.0, -0.6]),
        ),
        Box::new(UniversalJoint::new([0.0, 1.0, 0.0], [1.0, 0.0, 0.0])),
    );
    assert_eq!(model.total_dof(), 6 + 3 + 2);
    assert_eq!(model.total_nq(), 7 + 4 + 2);
    let (h1, h2) = (0.2_f64, 0.45_f64);
    let q = vec![
        0.1,
        0.2,
        0.3,
        h1.cos(),
        0.0,
        h1.sin(),
        0.0,
        h2.cos(),
        h2.sin() * 0.6,
        0.0,
        h2.sin() * 0.8,
        0.3,
        -0.5,
    ];
    assert_crba_matches_rnea(&model, &q, 1e-10);
}
