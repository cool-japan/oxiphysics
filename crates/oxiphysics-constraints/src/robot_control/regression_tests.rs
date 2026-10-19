//! Regression tests for the robot safety monitor and the Jacobian
//! manipulability index.

use super::types::{DhParams, Jacobian, RobotSafetyMonitor, SafetyLevel, SerialManipulator};
use std::f64::consts::FRAC_PI_2;

// ── RobotSafetyMonitor ──────────────────────────────────────────────────

/// Two joints, positions in [-1, 1] and [-2, 2], |velocity| <= 5, |torque| <= 20.
fn monitor() -> RobotSafetyMonitor {
    RobotSafetyMonitor::new(
        vec![[-1.0, 1.0], [-2.0, 2.0]],
        vec![5.0, 5.0],
        vec![20.0, 20.0],
    )
}

const REST: [f64; 2] = [0.0, 0.0];

#[test]
fn position_inside_the_warning_band_is_a_warning() {
    let mut mon = monitor();
    // The default warning band is the outer 10 % of the range: above 0.8 for joint 0.
    let level = mon.check(&[0.95, 0.0], &REST, &REST).clone();
    assert_eq!(level, SafetyLevel::Warning);
    assert_eq!(mon.reason, "Joint 0 near position limit");
    // The lower edge of the other joint's range warns as well.
    let level = mon.check(&[0.0, -1.9], &REST, &REST).clone();
    assert_eq!(level, SafetyLevel::Warning);
    assert_eq!(mon.reason, "Joint 1 near position limit");
    // Just inside the band is normal.
    assert_eq!(*mon.check(&[0.79, 1.5], &REST, &REST), SafetyLevel::Normal);
}

#[test]
fn velocity_and_torque_above_the_warning_fraction_warn() {
    let mut mon = monitor();
    assert_eq!(*mon.check(&REST, &[4.6, 0.0], &REST), SafetyLevel::Warning);
    assert_eq!(mon.reason, "Joint 0 near velocity limit");
    assert_eq!(*mon.check(&REST, &[-4.4, 0.0], &REST), SafetyLevel::Normal);
    assert_eq!(
        *mon.check(&REST, &REST, &[0.0, -19.0]),
        SafetyLevel::Warning
    );
    assert_eq!(mon.reason, "Joint 1 near torque limit");
    assert_eq!(
        *mon.check(&REST, &REST, &[17.0, -17.0]),
        SafetyLevel::Normal
    );
}

#[test]
fn a_warning_does_not_outlive_its_sample() {
    let mut mon = monitor();
    assert_eq!(*mon.check(&[0.95, 0.0], &REST, &REST), SafetyLevel::Warning);
    assert_eq!(*mon.check(&REST, &REST, &REST), SafetyLevel::Normal);
    assert!(mon.reason.is_empty(), "stale reason: {}", mon.reason);
    assert_eq!(mon.level, SafetyLevel::Normal);
}

#[test]
fn every_warning_of_a_sample_is_reported() {
    let mut mon = monitor();
    let level = mon.check(&[0.95, 0.0], &[0.0, 4.8], &[19.5, 0.0]).clone();
    assert_eq!(level, SafetyLevel::Warning);
    assert!(mon.reason.contains("Joint 0 near position limit"));
    assert!(mon.reason.contains("Joint 1 near velocity limit"));
    assert!(mon.reason.contains("Joint 0 near torque limit"));
}

#[test]
fn limit_violations_are_emergency_stops() {
    let mut mon = monitor();
    assert_eq!(
        *mon.check(&[1.2, 0.0], &REST, &REST),
        SafetyLevel::EmergencyStop
    );
    assert!(mon.reason.contains("position"), "{}", mon.reason);
    mon.reset();
    assert_eq!(
        *mon.check(&REST, &[0.0, -5.5], &REST),
        SafetyLevel::EmergencyStop
    );
    assert!(mon.reason.contains("velocity"), "{}", mon.reason);
    mon.reset();
    assert_eq!(
        *mon.check(&REST, &REST, &[21.0, 0.0]),
        SafetyLevel::EmergencyStop
    );
    assert!(mon.reason.contains("torque"), "{}", mon.reason);
}

#[test]
fn a_violation_outranks_a_warning_in_the_same_sample() {
    let mut mon = monitor();
    // Joint 0 only warns, joint 1 exceeds its velocity limit.
    let level = mon.check(&[0.95, 0.0], &[0.0, 6.0], &REST).clone();
    assert_eq!(level, SafetyLevel::EmergencyStop);
    assert!(mon.reason.contains("Joint 1 velocity"), "{}", mon.reason);
}

#[test]
fn non_finite_samples_are_emergency_stops() {
    let mut mon = monitor();
    for (q, qd, tau, what) in [
        ([f64::NAN, 0.0], REST, REST, "position"),
        ([0.0, f64::INFINITY], REST, REST, "position"),
        (REST, [0.0, f64::NAN], REST, "velocity"),
        (REST, [f64::NEG_INFINITY, 0.0], REST, "velocity"),
        (REST, REST, [f64::NAN, 0.0], "torque"),
    ] {
        mon.reset();
        let level = mon.check(&q, &qd, &tau).clone();
        assert_eq!(
            level,
            SafetyLevel::EmergencyStop,
            "{what} sample {q:?} {qd:?} {tau:?}"
        );
        assert!(
            mon.reason.starts_with("non-finite state"),
            "reason for {what}: {}",
            mon.reason
        );
        assert!(
            mon.reason.contains(what),
            "reason for {what}: {}",
            mon.reason
        );
    }
}

#[test]
fn non_finite_values_beyond_the_configured_joints_are_still_detected() {
    let mut mon = monitor();
    // Three joints reported although only two have limits.
    let level = mon.check(&[0.0, 0.0, f64::NAN], &REST, &REST).clone();
    assert_eq!(level, SafetyLevel::EmergencyStop);
    assert!(mon.reason.contains("joint 2"), "{}", mon.reason);
}

#[test]
fn an_emergency_stop_latches_until_reset() {
    let mut mon = monitor();
    assert_eq!(
        *mon.check(&[2.0, 0.0], &REST, &REST),
        SafetyLevel::EmergencyStop
    );
    // A perfectly normal sample does not clear the stop ...
    assert_eq!(*mon.check(&REST, &REST, &REST), SafetyLevel::EmergencyStop);
    // ... and neither does a warning-level one.
    assert_eq!(
        *mon.check(&[0.95, 0.0], &REST, &REST),
        SafetyLevel::EmergencyStop
    );
    mon.reset();
    assert_eq!(mon.level, SafetyLevel::Normal);
    assert!(mon.reason.is_empty());
    assert_eq!(*mon.check(&REST, &REST, &REST), SafetyLevel::Normal);
}

// ── Jacobian::manipulability ────────────────────────────────────────────

/// Serial chain from `(a, d, alpha)` rows (theta offset zero) at `angles`.
fn chain(rows: &[(f64, f64, f64)], angles: &[f64]) -> SerialManipulator {
    let links = rows
        .iter()
        .map(|&(a, d, alpha)| DhParams::new(a, d, alpha, 0.0))
        .collect();
    let mut robot = SerialManipulator::new(links);
    robot.set_joint_angles(angles);
    robot
}

/// Determinant of a row-major `n x n` matrix by LU with partial pivoting.
fn determinant(a: &[f64], n: usize) -> f64 {
    let mut lu = a.to_vec();
    let mut det = 1.0;
    for col in 0..n {
        let mut pivot = col;
        for row in (col + 1)..n {
            if lu[row * n + col].abs() > lu[pivot * n + col].abs() {
                pivot = row;
            }
        }
        if lu[pivot * n + col] == 0.0 {
            return 0.0;
        }
        if pivot != col {
            for k in 0..n {
                lu.swap(col * n + k, pivot * n + k);
            }
            det = -det;
        }
        det *= lu[col * n + col];
        for row in (col + 1)..n {
            let f = lu[row * n + col] / lu[col * n + col];
            for k in col..n {
                lu[row * n + k] -= f * lu[col * n + k];
            }
        }
    }
    det
}

/// `J Jᵀ` (6 x 6) of the row-major `6 x n` Jacobian.
fn jjt(j: &[f64], n: usize) -> Vec<f64> {
    let mut g = vec![0.0; 36];
    for r in 0..6 {
        for c in 0..6 {
            g[r * 6 + c] = (0..n).map(|k| j[r * n + k] * j[c * n + k]).sum();
        }
    }
    g
}

/// `Jᵀ J` (n x n) of the row-major `6 x n` Jacobian.
fn jtj(j: &[f64], n: usize) -> Vec<f64> {
    let mut g = vec![0.0; n * n];
    for r in 0..n {
        for c in 0..n {
            g[r * n + c] = (0..6).map(|k| j[k * n + r] * j[k * n + c]).sum();
        }
    }
    g
}

fn seven_dof() -> SerialManipulator {
    let h = FRAC_PI_2;
    chain(
        &[
            (0.0, 0.333, -h),
            (0.0, 0.0, h),
            (0.0825, 0.316, h),
            (-0.0825, 0.0, -h),
            (0.0, 0.384, h),
            (0.088, 0.0, h),
            (0.0, 0.107, 0.0),
        ],
        &[0.3, -0.5, 0.2, -1.4, 0.4, 1.1, -0.3],
    )
}

#[test]
fn redundant_arm_manipulability_is_sqrt_det_j_jt() {
    let robot = seven_dof();
    let jac = Jacobian::compute(&robot);
    let w = jac.manipulability();
    let expected = determinant(&jjt(&jac.j, 7), 6).sqrt();
    assert!(
        (w - expected).abs() <= 1e-9 * expected,
        "manipulability {w} vs sqrt(det(J J^T)) = {expected}"
    );
    assert!(
        (w - 0.060892469298682916).abs() < 1e-9,
        "reference value {w}"
    );
    // The product of the diagonal of J J^T (the previous implementation) is
    // only a Hadamard upper bound and is off by almost an order of magnitude.
    let g = jjt(&jac.j, 7);
    let diagonal_product: f64 = (0..6).map(|i| g[i * 6 + i]).product();
    assert!(
        diagonal_product.sqrt() > 5.0 * w,
        "diagonal product {} should be far above the true index {w}",
        diagonal_product.sqrt()
    );
}

#[test]
fn six_joint_arm_manipulability_is_the_absolute_jacobian_determinant() {
    let h = FRAC_PI_2;
    let robot = chain(
        &[
            (0.0, 0.333, -h),
            (0.0, 0.0, h),
            (0.0825, 0.316, h),
            (-0.0825, 0.0, -h),
            (0.0, 0.384, h),
            (0.088, 0.0, h),
        ],
        &[0.3, -0.5, 0.2, -1.4, 0.4, 1.1],
    );
    let jac = Jacobian::compute(&robot);
    let expected = determinant(&jac.j, 6).abs();
    let w = jac.manipulability();
    assert!(
        (w - expected).abs() <= 1e-9 * expected,
        "manipulability {w} vs |det J| = {expected}"
    );
    assert!(
        (w - 0.017287307177442274).abs() < 1e-9,
        "reference value {w}"
    );
}

#[test]
fn arm_with_fewer_than_six_joints_uses_sqrt_det_jt_j() {
    let h = FRAC_PI_2;
    let robot = chain(
        &[(0.0, 0.3, -h), (0.4, 0.0, 0.0), (0.35, 0.0, 0.0)],
        &[0.4, -0.7, 1.1],
    );
    let jac = Jacobian::compute(&robot);
    let expected = determinant(&jtj(&jac.j, 3), 3).sqrt();
    let w = jac.manipulability();
    assert!(
        (w - expected).abs() <= 1e-9 * expected,
        "manipulability {w} vs sqrt(det(J^T J)) = {expected}"
    );
    // det(J J^T) of a 6 x 3 Jacobian is identically zero, so the index must
    // not be computed from it.
    let g = jjt(&jac.j, 3);
    assert!(determinant(&g, 6).abs() < 1e-12);
    assert!(w > 0.1, "regular configuration, got {w}");
}

#[test]
fn coincident_joint_axes_make_the_index_vanish() {
    let h = FRAC_PI_2;
    // Joint 2 repeats joint 1 (same axis through the same point): two equal
    // columns, rank five.
    let singular = chain(
        &[
            (0.0, 0.3, -h),
            (0.0, 0.0, 0.0),
            (0.4, 0.0, 0.0),
            (0.35, 0.0, h),
            (0.0, 0.2, -h),
            (0.0, 0.1, 0.0),
        ],
        &[0.4, -0.7, 1.1, 0.3, -0.9, 0.5],
    );
    let jac = Jacobian::compute(&singular);
    let w = jac.manipulability();
    assert!(w < 1e-9, "singular chain must report ~0, got {w}");
    // The previous implementation reported a healthy 1.19 for this chain.
    let g = jjt(&jac.j, 6);
    let diagonal_product: f64 = (0..6).map(|i| g[i * 6 + i]).product();
    assert!(diagonal_product.sqrt() > 1.0);
    // Replacing the repeated axis by a perpendicular one restores the index.
    let regular = chain(
        &[
            (0.0, 0.3, -h),
            (0.0, 0.0, h),
            (0.4, 0.0, 0.0),
            (0.35, 0.0, h),
            (0.0, 0.2, -h),
            (0.0, 0.1, 0.0),
        ],
        &[0.4, -0.7, 1.1, 0.3, -0.9, 0.5],
    );
    let w_regular = Jacobian::compute(&regular).manipulability();
    assert!(
        (w_regular - 0.010277328623643936).abs() < 1e-9,
        "{w_regular}"
    );
}

#[test]
fn degenerate_jacobians_report_zero_manipulability() {
    let empty = Jacobian {
        dh: Vec::new(),
        joint_angles: Vec::new(),
        j: Vec::new(),
    };
    assert_eq!(empty.manipulability(), 0.0);

    let mut jac = Jacobian::compute(&seven_dof());
    let valid = jac.manipulability();
    assert!(valid > 0.0);
    jac.j[3] = f64::NAN;
    assert_eq!(
        jac.manipulability(),
        0.0,
        "NaN Jacobian is treated as singular"
    );
    jac.j[3] = f64::INFINITY;
    assert_eq!(
        jac.manipulability(),
        0.0,
        "infinite Jacobian is treated as singular"
    );

    let mut short = Jacobian::compute(&seven_dof());
    short.j.pop();
    assert_eq!(short.manipulability(), 0.0, "wrong-sized Jacobian");
}
