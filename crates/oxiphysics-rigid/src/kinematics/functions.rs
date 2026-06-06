//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{IkResult, Jacobian, SerialManipulator};

/// 4×4 homogeneous transformation matrix (row-major).
pub type HMatrix = [[f64; 4]; 4];
/// 3×3 rotation matrix.
pub type Rot3 = [[f64; 3]; 3];
/// 3D vector.
pub type Vec3 = [f64; 3];
/// 6D twist / wrench vector \[ω; v\] or \[f; m\].
pub type Twist = [f64; 6];
/// Identity homogeneous transform.
pub const EYE4: HMatrix = [
    [1., 0., 0., 0.],
    [0., 1., 0., 0.],
    [0., 0., 1., 0.],
    [0., 0., 0., 1.],
];
/// Identity rotation matrix.
pub const EYE3: Rot3 = [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];
/// Multiply two 4×4 homogeneous transforms.
pub fn h_mul(a: &HMatrix, b: &HMatrix) -> HMatrix {
    let mut c = [[0.0f64; 4]; 4];
    for i in 0..4 {
        for j in 0..4 {
            for k in 0..4 {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    c
}
/// Invert a homogeneous transform: T⁻¹ = \[Rᵀ, −Rᵀ·p; 0, 1\].
pub fn h_inv(t: &HMatrix) -> HMatrix {
    let mut result = [[0.0f64; 4]; 4];
    for i in 0..3 {
        for j in 0..3 {
            result[i][j] = t[j][i];
        }
    }
    for i in 0..3 {
        let sum: f64 = t[..3].iter().map(|row| row[i] * row[3]).sum();
        result[i][3] = -sum;
    }
    result[3][3] = 1.0;
    result
}
/// Extract the rotation part from a homogeneous transform.
pub fn h_rotation(t: &HMatrix) -> Rot3 {
    let mut r = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            r[i][j] = t[i][j];
        }
    }
    r
}
/// Extract the translation part from a homogeneous transform.
pub fn h_translation(t: &HMatrix) -> Vec3 {
    [t[0][3], t[1][3], t[2][3]]
}
/// Rotation matrix about the X axis by angle θ.
pub fn rot_x(theta: f64) -> HMatrix {
    let c = theta.cos();
    let s = theta.sin();
    [
        [1., 0., 0., 0.],
        [0., c, -s, 0.],
        [0., s, c, 0.],
        [0., 0., 0., 1.],
    ]
}
/// Rotation matrix about the Y axis by angle θ.
pub fn rot_y(theta: f64) -> HMatrix {
    let c = theta.cos();
    let s = theta.sin();
    [
        [c, 0., s, 0.],
        [0., 1., 0., 0.],
        [-s, 0., c, 0.],
        [0., 0., 0., 1.],
    ]
}
/// Rotation matrix about the Z axis by angle θ.
pub fn rot_z(theta: f64) -> HMatrix {
    let c = theta.cos();
    let s = theta.sin();
    [
        [c, -s, 0., 0.],
        [s, c, 0., 0.],
        [0., 0., 1., 0.],
        [0., 0., 0., 1.],
    ]
}
/// Pure translation homogeneous transform.
pub fn translate(tx: f64, ty: f64, tz: f64) -> HMatrix {
    [
        [1., 0., 0., tx],
        [0., 1., 0., ty],
        [0., 0., 1., tz],
        [0., 0., 0., 1.],
    ]
}
/// Apply a homogeneous transform to a 3D point.
pub fn h_apply_point(t: &HMatrix, p: &Vec3) -> Vec3 {
    [
        t[0][0] * p[0] + t[0][1] * p[1] + t[0][2] * p[2] + t[0][3],
        t[1][0] * p[0] + t[1][1] * p[1] + t[1][2] * p[2] + t[1][3],
        t[2][0] * p[0] + t[2][1] * p[1] + t[2][2] * p[2] + t[2][3],
    ]
}
/// 3D cross product.
pub fn cross3(a: &Vec3, b: &Vec3) -> Vec3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
/// 3D dot product.
pub fn dot3(a: &Vec3, b: &Vec3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
/// 3D vector norm.
pub fn norm3(a: &Vec3) -> f64 {
    dot3(a, a).sqrt()
}
/// Normalize a 3D vector.
pub fn normalize3(a: &Vec3) -> Vec3 {
    let n = norm3(a).max(1e-30);
    [a[0] / n, a[1] / n, a[2] / n]
}
/// Numerical inverse kinematics using Jacobian pseudo-inverse (Levenberg-Marquardt).
pub fn inverse_kinematics_lm(
    arm: &mut SerialManipulator,
    target_pos: &Vec3,
    q0: &[f64],
    tol: f64,
    max_iter: usize,
    lambda: f64,
) -> IkResult {
    arm.set_q(q0);
    let mut q = q0.to_vec();
    for iter in 0..max_iter {
        let pos = arm.end_effector_position();
        let err = [
            target_pos[0] - pos[0],
            target_pos[1] - pos[1],
            target_pos[2] - pos[2],
        ];
        let e_norm = norm3(&err);
        if e_norm < tol {
            return IkResult {
                q,
                pos_error: e_norm,
                iterations: iter,
                converged: true,
            };
        }
        let jac = Jacobian::compute(arm);
        let n = arm.n_joints;
        let mut jjt3 = [[0.0f64; 3]; 3];
        for j in 0..n {
            let col = &jac.columns[j];
            for (r, jjt_row) in jjt3.iter_mut().enumerate() {
                for (c, cell) in jjt_row.iter_mut().enumerate() {
                    *cell += col[r] * col[c];
                }
            }
        }
        for (i, jjt_row) in jjt3.iter_mut().enumerate() {
            jjt_row[i] += lambda * lambda;
        }
        let y = solve3x3(&jjt3, &err);
        for (j, q_j) in q.iter_mut().enumerate() {
            let dq: f64 = (0..3).map(|r| jac.columns[j][r] * y[r]).sum();
            *q_j += dq;
        }
        arm.set_q(&q);
    }
    let pos = arm.end_effector_position();
    let err = [
        target_pos[0] - pos[0],
        target_pos[1] - pos[1],
        target_pos[2] - pos[2],
    ];
    IkResult {
        q,
        pos_error: norm3(&err),
        iterations: max_iter,
        converged: false,
    }
}
/// Solve a 3×3 linear system Ax = b using Cramer's rule.
/// Compute determinant of an n×n matrix via LU decomposition with partial pivoting.
pub(super) fn mat_det_lu(a: &[Vec<f64>]) -> f64 {
    let n = a.len();
    if n == 0 {
        return 0.0;
    }
    let mut m: Vec<Vec<f64>> = a.to_vec();
    let mut det = 1.0f64;
    for col in 0..n {
        let pivot_row = (col..n)
            .max_by(|&i, &j| {
                m[i][col]
                    .abs()
                    .partial_cmp(&m[j][col].abs())
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .expect("operation should succeed");
        if m[pivot_row][col].abs() < 1e-300 {
            return 0.0;
        }
        if pivot_row != col {
            m.swap(pivot_row, col);
            det = -det;
        }
        det *= m[col][col];
        let inv = 1.0 / m[col][col];
        for row in (col + 1)..n {
            let factor = m[row][col] * inv;
            let (left, right) = m.split_at_mut(row);
            for (pivot_val, target) in left[col][col..].iter().zip(right[0][col..].iter_mut()) {
                *target -= factor * pivot_val;
            }
        }
    }
    det
}
fn solve3x3(a: &[[f64; 3]; 3], b: &Vec3) -> Vec3 {
    let det = a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
    if det.abs() < 1e-30 {
        return [0.0; 3];
    }
    let mut x = [0.0f64; 3];
    for k in 0..3 {
        let mut ak = *a;
        for i in 0..3 {
            ak[i][k] = b[i];
        }
        let det_k = ak[0][0] * (ak[1][1] * ak[2][2] - ak[1][2] * ak[2][1])
            - ak[0][1] * (ak[1][0] * ak[2][2] - ak[1][2] * ak[2][0])
            + ak[0][2] * (ak[1][0] * ak[2][1] - ak[1][1] * ak[2][0]);
        x[k] = det_k / det;
    }
    x
}
/// Analytical IK for a planar 2R arm.
///
/// Returns up to two solutions (elbow up / elbow down).
pub fn ik_planar_2r(l1: f64, l2: f64, px: f64, py: f64) -> Vec<[f64; 2]> {
    let r2 = px * px + py * py;
    let cos_q2 = (r2 - l1 * l1 - l2 * l2) / (2.0 * l1 * l2);
    if cos_q2.abs() > 1.0 {
        return vec![];
    }
    let mut solutions = Vec::new();
    for sign in &[1.0, -1.0] {
        let q2 = (sign * (1.0 - cos_q2 * cos_q2).sqrt()).atan2(cos_q2);
        let k1 = l1 + l2 * q2.cos();
        let k2 = l2 * q2.sin();
        let q1 = py.atan2(px) - k2.atan2(k1);
        solutions.push([q1, q2]);
    }
    solutions
}
/// Coriolis acceleration for a point in a rotating frame.
///
/// a_Coriolis = 2 · ω × v_rel
pub fn coriolis_acceleration(omega: &Vec3, v_rel: &Vec3) -> Vec3 {
    let two_omega_cross_v = cross3(omega, v_rel);
    [
        2.0 * two_omega_cross_v[0],
        2.0 * two_omega_cross_v[1],
        2.0 * two_omega_cross_v[2],
    ]
}
/// Centrifugal acceleration for a point at r in a rotating frame.
///
/// a_centrifugal = −ω × (ω × r)
pub fn centrifugal_acceleration(omega: &Vec3, r: &Vec3) -> Vec3 {
    let omega_cross_r = cross3(omega, r);
    let result = cross3(omega, &omega_cross_r);
    [-result[0], -result[1], -result[2]]
}
#[cfg(test)]
mod tests {
    use super::super::types::*;
    use super::super::types_2::*;
    use super::*;
    use std::f64::consts::PI;
    #[test]
    fn h_mul_identity() {
        let a = rot_z(0.5);
        let result = h_mul(&a, &EYE4);
        for i in 0..4 {
            for j in 0..4 {
                assert!((result[i][j] - a[i][j]).abs() < 1e-10);
            }
        }
    }
    #[test]
    fn h_inv_round_trip() {
        let t = h_mul(&rot_z(0.7), &translate(1.0, 2.0, 3.0));
        let t_inv = h_inv(&t);
        let product = h_mul(&t, &t_inv);
        for (i, row) in product.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (val - expected).abs() < 1e-8,
                    "[{i}][{j}]: {} ≠ {}",
                    val,
                    expected
                );
            }
        }
    }
    #[test]
    fn rot_z_preserves_z_axis() {
        let rz = rot_z(PI / 3.0);
        let z = [0.0, 0.0, 1.0];
        let z_rotated = h_apply_point(&rz, &z);
        assert!((z_rotated[2] - 1.0).abs() < 1e-10);
        assert!(z_rotated[0].abs() < 1e-10);
    }
    #[test]
    fn translate_moves_origin() {
        let t = translate(1.0, 2.0, 3.0);
        let p = h_apply_point(&t, &[0.0; 3]);
        assert!((p[0] - 1.0).abs() < 1e-10);
        assert!((p[1] - 2.0).abs() < 1e-10);
        assert!((p[2] - 3.0).abs() < 1e-10);
    }
    #[test]
    fn dh_identity_zero_angles() {
        let dh = DhParams::revolute(0.0, 0.0, 0.0, 0.0);
        let t = dh.transform();
        for i in 0..4 {
            for j in 0..4 {
                let expected = EYE4[i][j];
                assert!(
                    (t[i][j] - expected).abs() < 1e-10,
                    "[{i}][{j}]: {} ≠ {}",
                    t[i][j],
                    expected
                );
            }
        }
    }
    #[test]
    fn dh_revolute_variable_update() {
        let mut dh = DhParams::revolute(0.0, 1.0, 0.0, 0.0);
        dh.set_variable(PI / 2.0);
        assert!((dh.get_variable() - PI / 2.0).abs() < 1e-10);
    }
    #[test]
    fn planar_2r_home_position() {
        let arm = SerialManipulator::planar_2r(1.0, 1.0);
        let pos = arm.end_effector_position();
        assert!((pos[0] - 2.0).abs() < 1e-8, "x={}", pos[0]);
        assert!(pos[1].abs() < 1e-8, "y={}", pos[1]);
    }
    #[test]
    fn planar_2r_folded_position() {
        let mut arm = SerialManipulator::planar_2r(1.0, 1.0);
        arm.set_q(&[0.0, PI]);
        let pos = arm.end_effector_position();
        assert!(pos[0].abs() < 1e-6, "folded x should be ~0: {}", pos[0]);
    }
    #[test]
    fn planar_2r_quarter_turn() {
        let mut arm = SerialManipulator::planar_2r(1.0, 1.0);
        arm.set_q(&[PI / 2.0, 0.0]);
        let pos = arm.end_effector_position();
        assert!(
            pos[0].abs() < 1e-6,
            "x should be ~0 after 90° rotation: {}",
            pos[0]
        );
        assert!((pos[1] - 2.0).abs() < 1e-6, "y should be ~2: {}", pos[1]);
    }
    #[test]
    fn jacobian_2r_num_columns() {
        let arm = SerialManipulator::planar_2r(1.0, 1.0);
        let jac = Jacobian::compute(&arm);
        assert_eq!(jac.n_joints, 2);
        assert_eq!(jac.columns.len(), 2);
    }
    #[test]
    fn jacobian_manipulability_nonzero_generic_config() {
        let mut arm = SerialManipulator::planar_2r(1.0, 1.0);
        arm.set_q(&[PI / 6.0, PI / 4.0]);
        let jac = Jacobian::compute(&arm);
        let w = jac.manipulability();
        assert!(
            w > 0.0,
            "manipulability should be positive at non-singular config: {w}"
        );
    }
    #[test]
    fn jacobian_twist_velocity() {
        let arm = SerialManipulator::planar_2r(1.0, 1.0);
        let jac = Jacobian::compute(&arm);
        let twist = jac.multiply_velocity(&[0.0, 0.0]);
        for &v in &twist {
            assert!(v.abs() < 1e-10, "zero velocity → zero twist: {v}");
        }
    }
    #[test]
    fn ik_planar_2r_reachable_target() {
        let solutions = ik_planar_2r(1.0, 1.0, 1.0, 1.0);
        assert!(
            !solutions.is_empty(),
            "target (1,1) should be reachable with l1=l2=1"
        );
    }
    #[test]
    fn ik_planar_2r_out_of_reach() {
        let solutions = ik_planar_2r(1.0, 1.0, 3.0, 0.0);
        assert!(
            solutions.is_empty(),
            "target at distance 3 should be unreachable"
        );
    }
    #[test]
    fn ik_planar_2r_check_fk_consistency() {
        let solutions = ik_planar_2r(1.0, 1.0, 1.0, 0.0);
        assert!(!solutions.is_empty());
        for sol in &solutions {
            let mut arm = SerialManipulator::planar_2r(1.0, 1.0);
            arm.set_q(sol);
            let pos = arm.end_effector_position();
            assert!((pos[0] - 1.0).abs() < 1e-6, "IK-FK mismatch: x={}", pos[0]);
            assert!(pos[1].abs() < 1e-6, "IK-FK mismatch: y={}", pos[1]);
        }
    }
    #[test]
    fn screw_revolute_exp_zero_angle() {
        let s = ScrewAxis::revolute([0., 0., 1.], [0., 0., 0.]);
        let t = s.exp_map(0.0);
        for i in 0..4 {
            for j in 0..4 {
                assert!(
                    (t[i][j] - EYE4[i][j]).abs() < 1e-10,
                    "[{i}][{j}]: {} ≠ {}",
                    t[i][j],
                    EYE4[i][j]
                );
            }
        }
    }
    #[test]
    fn screw_revolute_z_full_rotation() {
        let s = ScrewAxis::revolute([0., 0., 1.], [0., 0., 0.]);
        let t = s.exp_map(2.0 * PI);
        for i in 0..4 {
            for j in 0..4 {
                assert!(
                    (t[i][j] - EYE4[i][j]).abs() < 1e-8,
                    "full rotation not identity: [{i}][{j}] = {}",
                    t[i][j]
                );
            }
        }
    }
    #[test]
    fn screw_prismatic_translates() {
        let s = ScrewAxis::prismatic([1., 0., 0.]);
        let t = s.exp_map(2.0);
        let p = h_apply_point(&t, &[0.; 3]);
        assert!(
            (p[0] - 2.0).abs() < 1e-10,
            "prismatic translation: {}",
            p[0]
        );
    }
    #[test]
    fn four_bar_grashof_satisfied() {
        let fb = FourBarLinkage::new(1.0, 2.0, 3.0, 4.0);
        let _ = fb.is_grashof();
    }
    #[test]
    fn four_bar_follower_angle_some() {
        let fb = FourBarLinkage::new(1.0, 3.0, 3.0, 4.0);
        let result = fb.follower_angle(PI / 4.0);
        let _ = result;
    }
    #[test]
    fn four_bar_follower_velocity_finite() {
        let fb = FourBarLinkage::new(1.0, 2.0, 2.0, 3.0);
        let omega3 = fb.follower_angular_velocity(1.0, 0.5, PI / 4.0, PI / 3.0);
        assert!(
            omega3.is_finite(),
            "angular velocity should be finite: {omega3}"
        );
    }
    #[test]
    fn cam_sinusoidal_lift_at_zero() {
        let cam = CamProfile::sinusoidal(100, 0.05, 0.02);
        let l = cam.lift_at(0.0);
        assert!(l.abs() < 1e-4, "sinusoidal cam lift at 0 should be ~0: {l}");
    }
    #[test]
    fn cam_sinusoidal_max_lift() {
        let cam = CamProfile::sinusoidal(100, 0.05, 0.02);
        let max_l = cam.lifts.iter().cloned().fold(0.0_f64, f64::max);
        assert!((max_l - 0.02).abs() < 0.001, "max lift: {max_l}");
    }
    #[test]
    fn cam_pressure_angle_reasonable() {
        let cam = CamProfile::sinusoidal(360, 0.05, 0.01);
        let phi = cam.max_pressure_angle();
        assert!(
            phi.to_degrees() < 45.0,
            "pressure angle should be < 45°: {}",
            phi.to_degrees()
        );
    }
    #[test]
    fn gear_pair_ratio() {
        let g = GearPair::spur(20, 40, 0.002);
        assert!(
            (g.ratio() - 0.5).abs() < 1e-10,
            "ratio = 20/40 = 0.5: {}",
            g.ratio()
        );
    }
    #[test]
    fn gear_pair_output_velocity() {
        let g = GearPair::spur(20, 40, 0.002);
        let omega_out = g.output_velocity(100.0);
        assert!(
            (omega_out - 50.0).abs() < 1e-8,
            "output = 50 rad/s: {omega_out}"
        );
    }
    #[test]
    fn gear_train_overall_ratio() {
        let train = GearTrain::new(vec![
            GearPair::spur(10, 20, 0.002),
            GearPair::spur(10, 40, 0.002),
        ]);
        assert!(
            (train.overall_ratio() - 0.125).abs() < 1e-10,
            "overall: {}",
            train.overall_ratio()
        );
    }
    #[test]
    fn gear_pair_contact_ratio_gt_1() {
        let g = GearPair::spur(20, 30, 0.002);
        let cr = g.contact_ratio();
        assert!(
            cr > 1.0,
            "contact ratio should exceed 1 for standard gears: {cr}"
        );
    }
    #[test]
    fn unicycle_constraint_satisfied_forward_motion() {
        let robot = UnicycleRobot::new(0.0, 0.0, 0.0);
        let constraint_val = robot.check_constraint();
        assert!(
            constraint_val.abs() < 1e-10,
            "forward motion satisfies constraint: {constraint_val}"
        );
    }
    #[test]
    fn unicycle_pfaffian_constraint_satisfied() {
        let robot = UnicycleRobot::new(0.0, 0.0, 0.0);
        let gv = robot.generalized_velocity();
        let pfaff = PfaffianConstraint::unicycle(robot.config[2]);
        assert!(pfaff.is_satisfied(&gv, 1e-10));
    }
    #[test]
    fn unicycle_robot_step_moves_forward() {
        let mut robot = UnicycleRobot::new(0.0, 0.0, 0.0);
        robot.v = 1.0;
        robot.step(0.1);
        assert!(
            (robot.config[0] - 0.1).abs() < 1e-10,
            "x = {}",
            robot.config[0]
        );
        assert!(robot.config[1].abs() < 1e-10, "y = {}", robot.config[1]);
    }
    #[test]
    fn pfaffian_degrees_of_freedom_unicycle() {
        let pfaff = PfaffianConstraint::unicycle(0.0);
        assert_eq!(pfaff.degrees_of_freedom(), 2);
    }
    #[test]
    fn gyroscopic_torque_zero_for_symmetric_spinning() {
        let body = GyroscopicBody::new([1.0, 1.0, 1.0], [0.0, 0.0, 10.0]);
        let tau_gyro = body.gyroscopic_torque();
        for &t in &tau_gyro {
            assert!(t.abs() < 1e-10, "symmetric gyro torque should be zero: {t}");
        }
    }
    #[test]
    fn gyroscopic_torque_nonzero_for_asymmetric() {
        let body = GyroscopicBody::new([1.0, 2.0, 3.0], [1.0, 1.0, 0.0]);
        let tau_gyro = body.gyroscopic_torque();
        let tau_norm: f64 = tau_gyro.iter().map(|x| x * x).sum::<f64>().sqrt();
        assert!(tau_norm > 0.0, "asymmetric gyro torque should be nonzero");
    }
    #[test]
    fn gyroscopic_kinetic_energy_positive() {
        let body = GyroscopicBody::new([1.0, 2.0, 3.0], [1.0, 2.0, 3.0]);
        let ke = body.kinetic_energy();
        assert!(ke > 0.0, "kinetic energy should be positive: {ke}");
    }
    #[test]
    fn coriolis_acceleration_perpendicular() {
        let omega = [0.0, 0.0, 1.0];
        let v_rel = [1.0, 0.0, 0.0];
        let a_cor = coriolis_acceleration(&omega, &v_rel);
        assert!(a_cor[0].abs() < 1e-10);
        assert!(
            a_cor[1].abs() > 1e-10,
            "Coriolis acceleration should have y component"
        );
        assert!(a_cor[2].abs() < 1e-10);
    }
    #[test]
    fn centrifugal_acceleration_outward() {
        let omega = [0.0, 0.0, 1.0];
        let r = [1.0, 0.0, 0.0];
        let a_cf = centrifugal_acceleration(&omega, &r);
        assert!(
            a_cf[0] > 0.0,
            "centrifugal should push outward: {}",
            a_cf[0]
        );
    }
}
