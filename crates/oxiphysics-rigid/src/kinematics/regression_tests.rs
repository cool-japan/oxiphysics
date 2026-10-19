// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the kinematics module: Jacobian manipulability
//! (any number of joints, no unit mixing, wrist singularities) and the
//! Pfaffian constraint projection / Frobenius integrability test.

use super::types::*;
use super::types_2::*;
use std::f64::consts::{FRAC_PI_2, PI};

// ── helpers ─────────────────────────────────────────────────────────────

/// Determinant by Gaussian elimination with partial pivoting.
fn determinant(mut m: Vec<Vec<f64>>) -> f64 {
    let n = m.len();
    let mut det = 1.0;
    for col in 0..n {
        let Some(pivot) = (col..n).max_by(|&i, &j| m[i][col].abs().total_cmp(&m[j][col].abs()))
        else {
            return 0.0;
        };
        if m[pivot][col] == 0.0 {
            return 0.0;
        }
        if pivot != col {
            m.swap(pivot, col);
            det = -det;
        }
        det *= m[col][col];
        for row in (col + 1)..n {
            let factor = m[row][col] / m[col][col];
            let pivot_row = m[col].clone();
            for (target, source) in m[row].iter_mut().zip(pivot_row.iter()).skip(col) {
                *target -= factor * source;
            }
        }
    }
    det
}

/// `J · Jᵀ` restricted to the given rows of the 6 x n Jacobian.
fn j_jt(jac: &Jacobian, rows: &[usize]) -> Vec<Vec<f64>> {
    rows.iter()
        .map(|&r| {
            rows.iter()
                .map(|&c| jac.columns.iter().map(|col| col[r] * col[c]).sum())
                .collect()
        })
        .collect()
}

/// Franka-like 7-axis arm (modified-DH-style table used as standard DH here);
/// only the first `n` joints are used.
fn seven_axis_arm(n: usize, q: &[f64]) -> SerialManipulator {
    let table = [
        (0.0, 0.333, -FRAC_PI_2),
        (0.0, 0.0, FRAC_PI_2),
        (0.0825, 0.316, FRAC_PI_2),
        (-0.0825, 0.0, -FRAC_PI_2),
        (0.0, 0.384, FRAC_PI_2),
        (0.088, 0.0, FRAC_PI_2),
        (0.0, 0.107, 0.0),
    ];
    let dh = table
        .iter()
        .take(n)
        .map(|&(a, d, alpha)| DhParams::revolute(alpha, a, d, 0.0))
        .collect();
    let mut arm = SerialManipulator::new(dh);
    arm.set_q(q);
    arm
}

const Q7: [f64; 7] = [0.3, -0.5, 0.2, -1.4, 0.4, 1.1, -0.3];

/// PUMA-560-style anthropomorphic arm with a spherical wrist (standard DH
/// table); only the first `n` joints are used.  The end point of the full arm
/// is the wrist centre, so joints 4-6 add no linear velocity there.
fn puma_style_arm(n: usize, q: &[f64]) -> SerialManipulator {
    let table = [
        (FRAC_PI_2, 0.0, 0.0),
        (0.0, 0.4318, 0.0),
        (-FRAC_PI_2, 0.0203, 0.15),
        (FRAC_PI_2, 0.0, 0.4318),
        (-FRAC_PI_2, 0.0, 0.0),
        (0.0, 0.0, 0.0),
    ];
    let dh = table
        .iter()
        .take(n)
        .map(|&(alpha, a, d)| DhParams::revolute(alpha, a, d, 0.0))
        .collect();
    let mut arm = SerialManipulator::new(dh);
    arm.set_q(q);
    arm
}

// ── manipulability ──────────────────────────────────────────────────────

#[test]
fn planar_2r_manipulability_is_l1_l2_abs_sin_q2() {
    for &(l1, l2, q1, q2) in &[
        (1.0, 1.0, 0.4, 0.9),
        (0.7, 1.3, -1.2, 2.1),
        (2.0, 0.5, 0.0, -0.6),
        (1.5, 1.5, 3.0, FRAC_PI_2),
    ] {
        let mut arm = SerialManipulator::planar_2r(l1, l2);
        arm.set_q(&[q1, q2]);
        let jac = Jacobian::compute(&arm);
        let expected = l1 * l2 * q2.sin().abs();
        let w = jac.manipulability();
        assert!(
            (w - expected).abs() < 1e-12 * (1.0 + expected),
            "w = {w}, expected {expected} for l = ({l1}, {l2}), q = ({q1}, {q2})"
        );
        assert!((jac.translational_manipulability() - expected).abs() < 1e-12 * (1.0 + expected));
    }
}

#[test]
fn planar_2r_is_singular_when_stretched_or_folded() {
    for q2 in [0.0, PI] {
        let mut arm = SerialManipulator::planar_2r(1.0, 0.8);
        arm.set_q(&[0.5, q2]);
        let jac = Jacobian::compute(&arm);
        assert!(
            jac.manipulability() < 1e-12,
            "q2 = {q2}: {}",
            jac.manipulability()
        );
        assert!(jac.is_singular(1e-9));
    }
    let mut arm = SerialManipulator::planar_2r(1.0, 0.8);
    arm.set_q(&[0.5, 1.0]);
    assert!(!Jacobian::compute(&arm).is_singular(1e-3));
}

#[test]
fn seven_joint_arm_manipulability_is_sqrt_det_j_jt() {
    let arm = seven_axis_arm(7, &Q7);
    let jac = Jacobian::compute(&arm);
    assert_eq!(jac.n_joints, 7);
    let w = jac.manipulability();
    let all_rows = [0usize, 1, 2, 3, 4, 5];
    let expected = determinant(j_jt(&jac, &all_rows)).sqrt();
    assert!(w > 1e-4, "regular 7-axis configuration, w = {w}");
    assert!(
        (w - expected).abs() <= 1e-9 * expected,
        "w = {w} vs sqrt(det(J J^T)) = {expected}"
    );
    assert!((jac.manipulability_full() - w).abs() <= 1e-12 * w);
    // The same index through the explicit 6x6 product.
    let jjt = jac.jjt();
    let from_jjt = determinant(jjt.iter().map(|row| row.to_vec()).collect()).sqrt();
    assert!((from_jjt - w).abs() <= 1e-9 * w);
}

#[test]
fn six_joint_arm_manipulability_is_abs_det_j() {
    let arm = seven_axis_arm(6, &Q7[..6]);
    let jac = Jacobian::compute(&arm);
    let square: Vec<Vec<f64>> = (0..6)
        .map(|r| (0..6).map(|c| jac.columns[c][r]).collect())
        .collect();
    let expected = determinant(square).abs();
    let w = jac.manipulability();
    assert!(
        expected > 1e-6,
        "regular configuration, |det J| = {expected}"
    );
    assert!(
        (w - expected).abs() <= 1e-9 * expected,
        "w = {w} vs {expected}"
    );
}

#[test]
fn four_and_five_joint_arms_use_translational_rows_only() {
    for n in [4usize, 5] {
        let arm = seven_axis_arm(n, &Q7[..n]);
        let mut jac = Jacobian::compute(&arm);
        let w = jac.manipulability();
        let expected = determinant(j_jt(&jac, &[0, 1, 2])).sqrt();
        assert!(expected > 1e-4, "n = {n}: translational index {expected}");
        assert!(
            (w - expected).abs() <= 1e-9 * expected,
            "n = {n}: w = {w} vs sqrt(det(Jv Jv^T)) = {expected}"
        );
        // Angular rows must not influence the index (no unit mixing).
        for col in jac.columns.iter_mut() {
            col[3] = 1.0e3;
            col[4] = -2.0e3;
            col[5] = 3.0e3;
        }
        assert_eq!(jac.manipulability(), w, "n = {n}: angular rows leaked in");
    }
}

#[test]
fn three_joint_arm_manipulability_is_abs_det_of_the_linear_block() {
    // Shoulder (vertical), shoulder pitch and elbow of an anthropomorphic arm.
    let dh = vec![
        DhParams::revolute(FRAC_PI_2, 0.0, 0.0, 0.0),
        DhParams::revolute(0.0, 0.5, 0.0, 0.0),
        DhParams::revolute(0.0, 0.4, 0.0, 0.0),
    ];
    let mut arm = SerialManipulator::new(dh);
    arm.set_q(&[0.3, 0.6, 0.9]);
    let jac = Jacobian::compute(&arm);
    let square: Vec<Vec<f64>> = (0..3)
        .map(|r| (0..3).map(|c| jac.columns[c][r]).collect())
        .collect();
    let expected = determinant(square).abs();
    assert!(
        expected > 1e-3,
        "regular configuration, |det Jv| = {expected}"
    );
    assert!((jac.manipulability() - expected).abs() <= 1e-12 * expected);
    // Closed form for this arm: l1 l2 sin(q3) (l1 cos q2 + l2 cos(q2 + q3)).
    let (q2, q3) = (0.6_f64, 0.9_f64);
    let closed_form = 0.5 * 0.4 * q3.sin() * (0.5 * q2.cos() + 0.4 * (q2 + q3).cos());
    assert!(
        (jac.manipulability() - closed_form.abs()).abs() < 1e-12,
        "{} vs {closed_form}",
        jac.manipulability()
    );
    // Stretched elbow: the arm loses a translational degree of freedom.
    arm.set_q(&[0.3, 0.6, 0.0]);
    assert!(Jacobian::compute(&arm).manipulability() < 1e-12);
}

#[test]
fn six_axis_wrist_singularity_is_detected_by_the_full_measure() {
    // Joints 4 and 6 become collinear at q5 = 0.
    let arm = puma_style_arm(6, &[0.3, -0.4, 0.5, 0.2, 0.0, 0.1]);
    let jac = Jacobian::compute(&arm);
    assert_eq!(jac.n_joints, 6);
    assert!(jac.manipulability() < 1e-12, "w = {}", jac.manipulability());
    assert!(jac.is_singular(1e-9));
    // The wrist centre is still fully mobile: a translational-only index
    // would have missed this singularity.
    assert!(
        jac.translational_manipulability() > 1e-3,
        "translational index {}",
        jac.translational_manipulability()
    );

    let arm = puma_style_arm(6, &[0.3, -0.4, 0.5, 0.2, 0.7, 0.1]);
    let jac = Jacobian::compute(&arm);
    assert!(jac.manipulability() > 1e-4, "w = {}", jac.manipulability());
    assert!(!jac.is_singular(1e-6));
    // For six joints the index is |det J|.
    let square: Vec<Vec<f64>> = (0..6)
        .map(|r| (0..6).map(|c| jac.columns[c][r]).collect())
        .collect();
    let expected = determinant(square).abs();
    assert!((jac.manipulability() - expected).abs() <= 1e-9 * expected);
}

#[test]
fn malformed_jacobians_report_zero() {
    let empty = Jacobian {
        columns: Vec::new(),
        n_joints: 0,
    };
    assert_eq!(empty.manipulability(), 0.0);
    assert_eq!(empty.translational_manipulability(), 0.0);
    assert_eq!(empty.manipulability_full(), 0.0);

    let arm = seven_axis_arm(7, &Q7);
    let valid = Jacobian::compute(&arm);
    assert!(valid.manipulability() > 0.0);

    let mut nan = valid.clone();
    nan.columns[2][1] = f64::NAN;
    assert_eq!(
        nan.manipulability(),
        0.0,
        "NaN entry is treated as singular"
    );

    let mut infinite = valid.clone();
    infinite.columns[0][4] = f64::INFINITY;
    assert_eq!(infinite.manipulability(), 0.0);

    let mut short = valid.clone();
    short.columns.pop();
    assert_eq!(short.manipulability(), 0.0, "column count != n_joints");
    let mut long = valid.clone();
    long.n_joints = 9;
    assert_eq!(long.manipulability(), 0.0);
    // The 6x6 product must not index out of range either.
    let _ = long.jjt();
}

// ── Pfaffian constraint projection ──────────────────────────────────────

fn norm(v: &[f64]) -> f64 {
    v.iter().map(|x| x * x).sum::<f64>().sqrt()
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

/// A velocity satisfying both bicycle constraints.
fn bicycle_admissible(
    theta: f64,
    phi: f64,
    wheelbase: f64,
    speed: f64,
    steer_rate: f64,
) -> Vec<f64> {
    vec![
        speed * theta.cos(),
        speed * theta.sin(),
        speed * phi.tan() / wheelbase,
        steer_rate,
    ]
}

#[test]
fn bicycle_projection_satisfies_both_constraints() {
    let (theta, phi, wheelbase) = (0.35, 0.4, 2.5);
    let constraint = PfaffianConstraint::bicycle(theta, phi, wheelbase);
    for v in [
        vec![1.0, 2.0, 3.0, 4.0],
        vec![-0.5, 0.25, 7.0, -3.0],
        vec![10.0, -10.0, 0.1, 0.0],
        vec![0.3, 0.0, 0.0, 1.0],
    ] {
        assert!(
            !constraint.is_satisfied(&v, 1e-6),
            "the raw test velocity {v:?} must violate the constraints"
        );
        let p = constraint.project_velocity(&v);
        assert_eq!(p.len(), 4);
        assert!(constraint.is_satisfied(&p, 1e-10), "{v:?} -> {p:?}");
        // Idempotent.
        let pp = constraint.project_velocity(&p);
        for (x, y) in p.iter().zip(&pp) {
            assert!((x - y).abs() < 1e-12, "projection must be idempotent");
        }
        // Orthogonal decomposition: v = p + r with r perpendicular to p and
        // the residual lies in the row space (its own projection vanishes).
        let r: Vec<f64> = v.iter().zip(&p).map(|(a, b)| a - b).collect();
        assert!(dot(&p, &r).abs() < 1e-10, "residual not orthogonal");
        assert!(norm(&constraint.project_velocity(&r)) < 1e-10);
        // The projection is not the zero vector for generic input.
        assert!(norm(&p) > 1e-3);
    }
}

#[test]
fn bicycle_projection_keeps_admissible_velocities_and_is_closest() {
    let (theta, phi, wheelbase) = (-0.8, 0.25, 1.7);
    let constraint = PfaffianConstraint::bicycle(theta, phi, wheelbase);
    let admissible = bicycle_admissible(theta, phi, wheelbase, 1.3, -0.4);
    assert!(constraint.is_satisfied(&admissible, 1e-12));
    let kept = constraint.project_velocity(&admissible);
    for (x, y) in kept.iter().zip(&admissible) {
        assert!((x - y).abs() < 1e-12, "{kept:?} vs {admissible:?}");
    }

    let v = [2.0, -1.0, 0.5, 3.0];
    let projected = constraint.project_velocity(&v);
    let best = norm(
        &v.iter()
            .zip(&projected)
            .map(|(a, b)| a - b)
            .collect::<Vec<_>>(),
    );
    for (speed, steer) in [(0.0, 0.0), (1.0, 1.0), (-2.0, 0.5), (3.0, -2.0), (0.7, 4.0)] {
        let other = bicycle_admissible(theta, phi, wheelbase, speed, steer);
        let distance = norm(&v.iter().zip(&other).map(|(a, b)| a - b).collect::<Vec<_>>());
        assert!(
            best <= distance + 1e-12,
            "projection must be the closest admissible velocity: {best} > {distance}"
        );
    }
}

#[test]
fn single_constraint_projection_matches_closed_form() {
    let theta = 0.9;
    let constraint = PfaffianConstraint::unicycle(theta);
    let v = [1.0, -2.0, 0.7];
    let a = [theta.sin(), -theta.cos(), 0.0];
    let lambda = dot(&a, &v) / dot(&a, &a);
    let expected: Vec<f64> = v.iter().zip(&a).map(|(x, y)| x - lambda * y).collect();
    let p = constraint.project_velocity(&v);
    for (x, y) in p.iter().zip(&expected) {
        assert!((x - y).abs() < 1e-14, "{p:?} vs {expected:?}");
    }
}

#[test]
fn rank_deficient_constraint_sets_are_projected_like_their_row_space() {
    let v = [3.0, 1.0, -2.0];
    let single = PfaffianConstraint {
        a_matrix: vec![vec![1.0, -1.0, 0.0]],
        n_constraints: 1,
        n_dof: 3,
    };
    let reference = single.project_velocity(&v);
    // a x - a y = 0 projection: (v_x - v_y)/2 removed along (1, -1, 0).
    assert!((reference[0] - 2.0).abs() < 1e-14 && (reference[1] - 2.0).abs() < 1e-14);
    assert!((reference[2] + 2.0).abs() < 1e-14);

    // Duplicated row, a scaled copy, a zero row and a linear combination.
    let redundant = PfaffianConstraint {
        a_matrix: vec![
            vec![1.0, -1.0, 0.0],
            vec![1.0, -1.0, 0.0],
            vec![-3.0, 3.0, 0.0],
            vec![0.0, 0.0, 0.0],
            vec![2.0, -2.0, 0.0],
        ],
        n_constraints: 5,
        n_dof: 3,
    };
    let p = redundant.project_velocity(&v);
    for (x, y) in p.iter().zip(&reference) {
        assert!((x - y).abs() < 1e-12, "{p:?} vs {reference:?}");
    }
    assert!(redundant.is_satisfied(&p, 1e-12));

    // Bicycle with zero wheelbase and zero steering: both rows coincide.
    let degenerate = PfaffianConstraint::bicycle(0.4, 0.0, 0.0);
    let w = [1.0, 2.0, 3.0, 4.0];
    let p = degenerate.project_velocity(&w);
    assert!(degenerate.is_satisfied(&p, 1e-12));
    let unicycle_like = PfaffianConstraint::unicycle(0.4);
    let expected = unicycle_like.project_velocity(&[1.0, 2.0, 3.0]);
    for k in 0..3 {
        assert!((p[k] - expected[k]).abs() < 1e-12);
    }
    assert!((p[3] - 4.0).abs() < 1e-12, "free coordinate untouched");
}

#[test]
fn projection_without_constraints_is_the_identity_and_odd_lengths_do_not_panic() {
    let free = PfaffianConstraint {
        a_matrix: Vec::new(),
        n_constraints: 0,
        n_dof: 3,
    };
    assert_eq!(free.project_velocity(&[1.0, 2.0, 3.0]), vec![1.0, 2.0, 3.0]);
    assert!(free.project_velocity(&[]).is_empty());

    let unicycle = PfaffianConstraint::unicycle(0.5);
    // Longer than n_dof: the extra coordinate is unconstrained.
    let long = unicycle.project_velocity(&[1.0, 0.0, 0.0, 9.0]);
    assert_eq!(long.len(), 4);
    assert!((long[3] - 9.0).abs() < 1e-15);
    // Shorter than n_dof: rows are truncated.
    let short = unicycle.project_velocity(&[1.0, 0.0]);
    assert_eq!(short.len(), 2);
    assert!(unicycle.is_satisfied(&[short[0], short[1], 0.0], 1e-12));
}

// ── Frobenius integrability ─────────────────────────────────────────────

const STEP: f64 = 1e-6;
const TOL: f64 = 1e-5;

#[test]
fn unicycle_and_bicycle_are_nonholonomic() {
    let q = [0.2, -0.1, 0.3];
    assert!(!PfaffianConstraint::is_integrable_at(
        |q| PfaffianConstraint::unicycle(q[2]).a_matrix,
        &q,
        STEP,
        TOL
    ));
    let q = [0.0, 0.0, 0.3, 0.4];
    assert!(!PfaffianConstraint::is_integrable_at(
        |q| PfaffianConstraint::bicycle(q[2], q[3], 2.0).a_matrix,
        &q,
        STEP,
        TOL
    ));
}

#[test]
fn constant_coefficient_constraints_are_integrable() {
    let q = [0.1, 0.2, 0.3, 0.4];
    let constant = |_: &[f64]| vec![vec![1.0, -1.0, 0.0, 0.5], vec![0.0, 2.0, 1.0, 0.0]];
    assert!(PfaffianConstraint::is_integrable_at(
        constant, &q, STEP, TOL
    ));
}

#[test]
fn exact_forms_are_integrable() {
    // omega = d(x y z) = (y z, x z, x y)
    let q = [1.0, 2.0, 3.0];
    let exact = |q: &[f64]| vec![vec![q[1] * q[2], q[0] * q[2], q[0] * q[1]]];
    assert!(PfaffianConstraint::is_integrable_at(exact, &q, STEP, TOL));
}

#[test]
fn integrable_but_not_closed_forms_are_recognised() {
    // omega = z dx - x dz = z^2 d(x/z) is integrable (integrating factor
    // 1/z^2) although d(omega) != 0, i.e. the derivative-symmetry condition
    // da_i/dq_j = da_j/dq_i fails: da_x/dz = 1 but da_z/dx = -1.
    let form = |q: &[f64]| vec![vec![q[2], 0.0, -q[0]]];
    let q = [1.0, 0.5, 2.0];
    let h = 1e-4;
    let da_x_dz =
        (form(&[q[0], q[1], q[2] + h])[0][0] - form(&[q[0], q[1], q[2] - h])[0][0]) / (2.0 * h);
    let da_z_dx =
        (form(&[q[0] + h, q[1], q[2]])[0][2] - form(&[q[0] - h, q[1], q[2]])[0][2]) / (2.0 * h);
    assert!((da_x_dz - 1.0).abs() < 1e-6 && (da_z_dx + 1.0).abs() < 1e-6);
    assert!(PfaffianConstraint::is_integrable_at(form, &q, STEP, TOL));

    // omega = z dx is integrable too (integrating factor 1/z).
    let z_dx = |q: &[f64]| vec![vec![q[2], 0.0, 0.0]];
    assert!(PfaffianConstraint::is_integrable_at(z_dx, &q, STEP, TOL));
}

#[test]
fn contact_form_is_not_integrable() {
    // omega = dz - y dx is the standard non-integrable (contact) distribution.
    let contact = |q: &[f64]| vec![vec![-q[1], 0.0, 1.0]];
    assert!(!PfaffianConstraint::is_integrable_at(
        contact,
        &[0.3, 0.7, -0.2],
        STEP,
        TOL
    ));
}

#[test]
fn a_single_constraint_in_two_dimensions_is_always_integrable() {
    let q = [0.4, 0.9];
    let form = |q: &[f64]| vec![vec![q[1].sin(), -q[0].cos()]];
    assert!(PfaffianConstraint::is_integrable_at(form, &q, STEP, TOL));
}
