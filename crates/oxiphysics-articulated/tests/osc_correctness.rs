// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

use oxiphysics_articulated::{
    ArticulatedError,
    body::RigidBody,
    crba::compute_mass_matrix_crba,
    ik::geometric_jacobian,
    joint::RevoluteJoint,
    model::ArticulatedModel,
    osc::{compute_lambda, compute_lambda_damped},
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
    let (model, q) = make_6dof_chain();
    let jacobian = full_rank_jacobian_6x6();
    let lambda = compute_lambda(&model, &q, &jacobian).expect("full-rank 6-DoF task must invert");
    assert_eq!(lambda.len(), 6, "Λ must be 6×6");
    assert_eq!(lambda[0].len(), 6, "Λ must be 6×6");
}

#[test]
fn test_osc_symmetric() {
    let (model, q) = make_6dof_chain();
    let jacobian = full_rank_jacobian_6x6();
    let lambda = compute_lambda(&model, &q, &jacobian).expect("full-rank 6-DoF task must invert");
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
    let (model, q) = make_6dof_chain();
    let jacobian = full_rank_jacobian_6x6();
    let lambda = compute_lambda(&model, &q, &jacobian).expect("full-rank 6-DoF task must invert");
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

// ─── Rank-deficient tasks and error handling ─────────────────────────────────

/// 3-DoF planar arm (all joints about z), links 0.5 m along x.
fn make_planar_3dof() -> (ArticulatedModel, Vec<f64>) {
    let mut model = ArticulatedModel::new([0.0, 0.0, -9.81]);
    for i in 0..3 {
        let parent = if i == 0 { None } else { Some(i - 1) };
        let x_t = if i == 0 {
            SpatialTransform::IDENTITY
        } else {
            SpatialTransform::from_translation([0.5, 0.0, 0.0])
        };
        let inertia = SpatialInertia::from_com(
            1.0,
            [0.25, 0.0, 0.0],
            [[0.001, 0.0, 0.0], [0.0, 0.02, 0.0], [0.0, 0.0, 0.02]],
        );
        model.add_body(
            RigidBody::new(format!("l{i}"), inertia, parent, x_t),
            Box::new(RevoluteJoint::new([0.0, 0.0, 1.0])),
        );
    }
    (model, vec![0.3, -0.6, 0.9])
}

/// The 6×3 end-effector Jacobian as row vectors.
fn jacobian_rows(model: &ArticulatedModel, q: &[f64], ee: usize) -> Vec<Vec<f64>> {
    let n = model.total_dof();
    let flat = geometric_jacobian(model, q, ee);
    (0..6).map(|r| flat[r * n..(r + 1) * n].to_vec()).collect()
}

#[test]
fn test_osc_three_dof_task_is_singular_not_panic() {
    let (model, q) = make_planar_3dof();
    let jac = jacobian_rows(&model, &q, 2);
    let result = compute_lambda(&model, &q, &jac);
    assert_eq!(result, Err(ArticulatedError::SingularTaskInertia));
}

#[test]
fn test_osc_damped_three_dof_task_is_regularised() {
    let (model, q) = make_planar_3dof();
    let jac = jacobian_rows(&model, &q, 2);
    let damping = 0.05;
    let lambda = compute_lambda_damped(&model, &q, &jac, damping)
        .expect("damped task inertia must be invertible");
    // Verify Λ·(J·M⁻¹·Jᵀ + λ²I) = I using an independent dense solve of M.
    let m = compute_mass_matrix_crba(&model, &q);
    let n = model.total_dof();
    let minv_jt: Vec<Vec<f64>> = jac.iter().map(|row| solve_dense(&m, row)).collect();
    let mut a = [[0.0f64; 6]; 6];
    for (r, (a_row, jac_r)) in a.iter_mut().zip(jac.iter()).enumerate() {
        for (a_rc, minv_c) in a_row.iter_mut().zip(minv_jt.iter()) {
            *a_rc = (0..n).map(|k| jac_r[k] * minv_c[k]).sum::<f64>();
        }
        a_row[r] += damping * damping;
    }
    for (r, lambda_r) in lambda.iter().enumerate() {
        for c in 0..6 {
            let prod: f64 = lambda_r
                .iter()
                .zip(a.iter())
                .map(|(l, a_k)| l * a_k[c])
                .sum();
            let expected = if r == c { 1.0 } else { 0.0 };
            assert!((prod - expected).abs() < 1e-8, "(Λ·A)[{r}][{c}]={prod:.3e}");
        }
    }
    // Λ is symmetric positive-definite.
    for (r, lambda_r) in lambda.iter().enumerate() {
        assert!(lambda_r[r] > 0.0, "Λ[{r}][{r}]={}", lambda_r[r]);
        for (c, lambda_c) in lambda.iter().enumerate() {
            assert!((lambda_r[c] - lambda_c[r]).abs() < 1e-8);
        }
    }
}

#[test]
fn test_osc_full_rank_matches_inverse_of_task_inertia() {
    let (model, q) = make_6dof_chain();
    let jac = jacobian_rows(&model, &q, 5);
    let lambda = compute_lambda(&model, &q, &jac).expect("6-DoF arm away from singularity");
    let lambda_zero_damping =
        compute_lambda_damped(&model, &q, &jac, 0.0).expect("λ = 0 equals the undamped case");
    let m = compute_mass_matrix_crba(&model, &q);
    let n = model.total_dof();
    let minv_jt: Vec<Vec<f64>> = jac.iter().map(|row| solve_dense(&m, row)).collect();
    for r in 0..6 {
        for c in 0..6 {
            assert!((lambda[r][c] - lambda_zero_damping[r][c]).abs() < 1e-12);
            let prod: f64 = (0..6)
                .map(|k| lambda[r][k] * (0..n).map(|i| jac[k][i] * minv_jt[c][i]).sum::<f64>())
                .sum();
            let expected = if r == c { 1.0 } else { 0.0 };
            assert!((prod - expected).abs() < 1e-7, "(Λ·A)[{r}][{c}]={prod:.3e}");
        }
    }
}

#[test]
fn test_osc_dimension_and_parameter_errors() {
    let (model, q) = make_planar_3dof();
    let jac = jacobian_rows(&model, &q, 2);
    let five_rows: Vec<Vec<f64>> = jac[..5].to_vec();
    assert_eq!(
        compute_lambda(&model, &q, &five_rows),
        Err(ArticulatedError::DimensionMismatch {
            what: "jacobian rows",
            expected: 6,
            got: 5
        })
    );
    let mut short_row = jac.clone();
    short_row[2].pop();
    assert!(matches!(
        compute_lambda(&model, &q, &short_row),
        Err(ArticulatedError::DimensionMismatch {
            what: "jacobian row length",
            ..
        })
    ));
    assert!(matches!(
        compute_lambda(&model, &q[..2], &jac),
        Err(ArticulatedError::DimensionMismatch { what: "q", .. })
    ));
    assert!(matches!(
        compute_lambda_damped(&model, &q, &jac, -1.0),
        Err(ArticulatedError::InvalidParameter {
            what: "damping",
            ..
        })
    ));
}

/// Dense Gaussian elimination with partial pivoting: solve `m·x = b`.
fn solve_dense(m: &[Vec<f64>], b: &[f64]) -> Vec<f64> {
    let n = b.len();
    let mut a: Vec<Vec<f64>> = m
        .iter()
        .zip(b.iter())
        .map(|(row, &bi)| {
            let mut r = row.clone();
            r.push(bi);
            r
        })
        .collect();
    for col in 0..n {
        let piv = (col..n)
            .max_by(|&x, &y| a[x][col].abs().total_cmp(&a[y][col].abs()))
            .expect("non-empty pivot range");
        a.swap(col, piv);
        let pivot_row = a[col].clone();
        for row in a.iter_mut().skip(col + 1) {
            let f = row[col] / pivot_row[col];
            for (entry, p) in row.iter_mut().zip(pivot_row.iter()).skip(col) {
                *entry -= f * p;
            }
        }
    }
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let s: f64 = (i + 1..n).map(|k| a[i][k] * x[k]).sum();
        x[i] = (a[i][n] - s) / a[i][i];
    }
    x
}
