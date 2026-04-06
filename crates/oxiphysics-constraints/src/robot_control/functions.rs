//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
use std::f64::consts::PI;

use super::types::SerialManipulator;

/// Compute the 3-D cross product of arrays `a` and `b`.
pub(super) fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
/// Dot product of two 3-vectors.
pub(super) fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
/// Euclidean norm of a 3-vector.
pub(super) fn norm3(v: [f64; 3]) -> f64 {
    dot3(v, v).sqrt()
}
/// Normalise a 3-vector, returning `[0,0,0]` for near-zero input.
pub(super) fn normalize3(v: [f64; 3]) -> [f64; 3] {
    let n = norm3(v);
    if n < 1e-15 {
        [0.0; 3]
    } else {
        [v[0] / n, v[1] / n, v[2] / n]
    }
}
/// Multiply a 4×4 column-major matrix (stored row-wise as `[f64; 16]`) by another.
///
/// Layout: element `(row, col)` lives at index `row * 4 + col`.
pub fn mat4_mul(a: [f64; 16], b: [f64; 16]) -> [f64; 16] {
    let mut out = [0.0f64; 16];
    for row in 0..4 {
        for col in 0..4 {
            let mut s = 0.0;
            for k in 0..4 {
                s += a[row * 4 + k] * b[k * 4 + col];
            }
            out[row * 4 + col] = s;
        }
    }
    out
}
/// Identity 4×4 matrix.
pub(super) fn mat4_identity() -> [f64; 16] {
    let mut m = [0.0f64; 16];
    m[0] = 1.0;
    m[5] = 1.0;
    m[10] = 1.0;
    m[15] = 1.0;
    m
}
/// Extract the 3×3 rotation submatrix from a 4×4 transform as `[[f64;3\];3]`.
pub(super) fn mat4_rotation(t: [f64; 16]) -> [[f64; 3]; 3] {
    [[t[0], t[1], t[2]], [t[4], t[5], t[6]], [t[8], t[9], t[10]]]
}
/// Extract the translation column from a 4×4 transform.
pub(super) fn mat4_translation(t: [f64; 16]) -> [f64; 3] {
    [t[3], t[7], t[11]]
}
/// Build the 4×4 DH homogeneous transform for the given parameters.
///
/// Row-major layout: index `i*4 + j` gives element `(row=i, col=j)`.
pub fn dh_matrix(a: f64, d: f64, alpha: f64, theta: f64) -> [f64; 16] {
    let ct = theta.cos();
    let st = theta.sin();
    let ca = alpha.cos();
    let sa = alpha.sin();
    [
        ct,
        -st * ca,
        st * sa,
        a * ct,
        st,
        ct * ca,
        -ct * sa,
        a * st,
        0.0,
        sa,
        ca,
        d,
        0.0,
        0.0,
        0.0,
        1.0,
    ]
}
/// Convert a 3×3 rotation matrix to ZYX Euler angles `[roll, pitch, yaw]` (radians).
pub fn rotation_to_euler(r: [[f64; 3]; 3]) -> [f64; 3] {
    let pitch = (-r[2][0]).asin().clamp(-PI / 2.0, PI / 2.0);
    let cos_pitch = pitch.cos();
    if cos_pitch.abs() < 1e-6 {
        let roll = r[0][1].atan2(r[1][1]);
        [roll, pitch, 0.0]
    } else {
        let yaw = r[1][0].atan2(r[0][0]);
        let roll = r[2][1].atan2(r[2][2]);
        [roll, pitch, yaw]
    }
}
/// Build a 3×3 rotation matrix from ZYX Euler angles `[roll, pitch, yaw]` (radians).
pub fn euler_to_rotation(roll: f64, pitch: f64, yaw: f64) -> [[f64; 3]; 3] {
    let cr = roll.cos();
    let sr = roll.sin();
    let cp = pitch.cos();
    let sp = pitch.sin();
    let cy = yaw.cos();
    let sy = yaw.sin();
    [
        [cy * cp, cy * sp * sr - sy * cr, cy * sp * cr + sy * sr],
        [sy * cp, sy * sp * sr + cy * cr, sy * sp * cr - cy * sr],
        [-sp, cp * sr, cp * cr],
    ]
}
/// Compute the 6×6 adjoint (spatial transform) matrix from a 4×4 homogeneous transform.
///
/// `Ad(T) = [[R, 0\], [p×R, R]]` in column-major vector form stored as `Vec`f64` of length 36.
pub fn adjoint_matrix(t: [f64; 16]) -> Vec<f64> {
    let r = mat4_rotation(t);
    let p = mat4_translation(t);
    let px = [[0.0, -p[2], p[1]], [p[2], 0.0, -p[0]], [-p[1], p[0], 0.0]];
    let mut pxr = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            pxr[i][j] = px[i][0] * r[0][j] + px[i][1] * r[1][j] + px[i][2] * r[2][j];
        }
    }
    let mut ad = vec![0.0f64; 36];
    for i in 0..3 {
        for j in 0..3 {
            ad[i * 6 + j] = r[i][j];
        }
    }
    for i in 0..3 {
        for j in 0..3 {
            ad[(i + 3) * 6 + (j + 3)] = r[i][j];
        }
    }
    for i in 0..3 {
        for j in 0..3 {
            ad[(i + 3) * 6 + j] = pxr[i][j];
        }
    }
    ad
}
/// Gauss-Jordan inversion of a 6×6 matrix (row-major flat `\[f64; 36\]`-compatible `Vec`).
pub(super) fn gauss_jordan_6x6(m: &[f64]) -> Vec<f64> {
    let n = 6usize;
    let mut aug = vec![0.0f64; n * n * 2];
    for i in 0..n {
        for j in 0..n {
            aug[i * (n * 2) + j] = m[i * n + j];
        }
        aug[i * (n * 2) + (n + i)] = 1.0;
    }
    for col in 0..n {
        let mut max_row = col;
        let mut max_val = aug[col * (n * 2) + col].abs();
        for row in (col + 1)..n {
            let v = aug[row * (n * 2) + col].abs();
            if v > max_val {
                max_val = v;
                max_row = row;
            }
        }
        if max_val < 1e-14 {
            continue;
        }
        for j in 0..(n * 2) {
            aug.swap(col * (n * 2) + j, max_row * (n * 2) + j);
        }
        let pivot = aug[col * (n * 2) + col];
        for j in 0..(n * 2) {
            aug[col * (n * 2) + j] /= pivot;
        }
        for row in 0..n {
            if row == col {
                continue;
            }
            let factor = aug[row * (n * 2) + col];
            for j in 0..(n * 2) {
                let val = aug[col * (n * 2) + j];
                aug[row * (n * 2) + j] -= factor * val;
            }
        }
    }
    let mut inv = vec![0.0f64; n * n];
    for i in 0..n {
        for j in 0..n {
            inv[i * n + j] = aug[i * (n * 2) + (n + j)];
        }
    }
    inv
}
/// Damped pseudoinverse of a 3×n matrix (row-major).
pub(super) fn damped_pseudoinverse_3xn(j: &[f64], n: usize, lambda: f64) -> Vec<f64> {
    let mut jjt = [[0.0f64; 3]; 3];
    for r in 0..3 {
        for c in 0..3 {
            let mut s = 0.0;
            for k in 0..n {
                s += j[r * n + k] * j[c * n + k];
            }
            jjt[r][c] = s;
        }
        jjt[r][r] += lambda * lambda;
    }
    let inv = invert3x3(jjt);
    let mut jp = vec![0.0f64; n * 3];
    for i in 0..n {
        for c in 0..3 {
            let mut s = 0.0;
            for r in 0..3 {
                s += j[r * n + i] * inv[r][c];
            }
            jp[i * 3 + c] = s;
        }
    }
    jp
}
/// Invert a 3×3 matrix.  Returns identity on near-singular inputs.
pub(super) fn invert3x3(m: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    if det.abs() < 1e-15 {
        return [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    }
    let inv_det = 1.0 / det;
    [
        [
            (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * inv_det,
            (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * inv_det,
            (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * inv_det,
        ],
        [
            (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * inv_det,
            (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * inv_det,
            (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * inv_det,
        ],
        [
            (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * inv_det,
            (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * inv_det,
            (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * inv_det,
        ],
    ]
}
/// Gradient of the joint-limit avoidance potential.
pub(super) fn null_space_joint_limit_gradient(robot: &SerialManipulator) -> Vec<f64> {
    let n = robot.dof();
    let mut g = vec![0.0f64; n];
    for i in 0..n {
        let (lo, hi) = robot.joint_limits[i];
        let mid = (lo + hi) * 0.5;
        let range = (hi - lo).max(1e-10);
        g[i] = -(robot.joint_angles[i] - mid) / (range * range);
    }
    g
}
/// Project a null-space gradient into the actual null-space of J (3×n).
pub(super) fn project_null_space(jp: &[f64], j: &[f64], n: usize, g: &[f64], w: f64) -> Vec<f64> {
    let mut jpj = vec![0.0f64; n * n];
    for i in 0..n {
        for k in 0..n {
            let mut s = 0.0;
            for r in 0..3 {
                s += jp[i * 3 + r] * j[r * n + k];
            }
            jpj[i * n + k] = s;
        }
    }
    let mut out = vec![0.0f64; n];
    for i in 0..n {
        let mut s = g[i];
        for k in 0..n {
            s -= jpj[i * n + k] * g[k];
        }
        out[i] = w * s;
    }
    out
}
/// Build a 3×3 rotation matrix from axis-angle representation.
pub fn rotation_from_axis_angle(axis: [f64; 3], angle: f64) -> [[f64; 3]; 3] {
    let c = angle.cos();
    let s = angle.sin();
    let t = 1.0 - c;
    let [x, y, z] = axis;
    [
        [t * x * x + c, t * x * y - s * z, t * x * z + s * y],
        [t * x * y + s * z, t * y * y + c, t * y * z - s * x],
        [t * x * z - s * y, t * y * z + s * x, t * z * z + c],
    ]
}
/// Multiply 3×3 matrix by 3-vector.
pub(super) fn mat3_vec_mul(m: [[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::robot_control::AdaptiveRobotControl;

    use crate::robot_control::CartesianTrajectory;
    use crate::robot_control::CollisionAvoidance;

    use crate::robot_control::DhParams;

    use crate::robot_control::ImpedanceControl;
    use crate::robot_control::InverseKinematics;
    use crate::robot_control::Jacobian;

    use crate::robot_control::TorqueFeedforward;

    use crate::robot_control::TrajectoryPlanning;

    fn make_2dof_robot() -> SerialManipulator {
        SerialManipulator::new(vec![
            DhParams::new(1.0, 0.0, 0.0, 0.0),
            DhParams::new(1.0, 0.0, 0.0, 0.0),
        ])
    }
    #[test]
    fn test_dh_identity_at_zero() {
        let m = dh_matrix(0.0, 0.0, 0.0, 0.0);
        assert!((m[0] - 1.0).abs() < 1e-10);
        assert!((m[5] - 1.0).abs() < 1e-10);
        assert!((m[10] - 1.0).abs() < 1e-10);
        assert!((m[15] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_dh_translation_x() {
        let m = dh_matrix(1.0, 0.0, 0.0, 0.0);
        assert!((m[3] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_dh_rotation_theta() {
        let m = dh_matrix(0.0, 0.0, 0.0, PI / 2.0);
        assert!(m[0].abs() < 1e-10);
        assert!((m[4] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_mat4_mul_identity() {
        let i = mat4_identity();
        let m = dh_matrix(1.0, 2.0, 0.3, 0.5);
        let r = mat4_mul(i, m);
        for k in 0..16 {
            assert!((r[k] - m[k]).abs() < 1e-10, "k={k}: {} vs {}", r[k], m[k]);
        }
    }
    #[test]
    fn test_serial_fk_zero_angles() {
        let robot = make_2dof_robot();
        let pos = robot.end_effector_position();
        assert!((pos[0] - 2.0).abs() < 1e-9, "x={}", pos[0]);
        assert!(pos[1].abs() < 1e-9);
        assert!(pos[2].abs() < 1e-9);
    }
    #[test]
    fn test_serial_fk_rotated_joint() {
        let mut robot = make_2dof_robot();
        robot.set_joint_angles(&[PI / 2.0, 0.0]);
        let pos = robot.end_effector_position();
        assert!(pos[0].abs() < 1e-9, "x={}", pos[0]);
        assert!((pos[1] - 2.0).abs() < 1e-9, "y={}", pos[1]);
    }
    #[test]
    fn test_joint_angle_clamping() {
        let mut robot = make_2dof_robot();
        robot.joint_limits[0] = (-1.0, 1.0);
        robot.set_joint_angles(&[5.0, 0.0]);
        assert!((robot.joint_angles[0] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_jacobian_size() {
        let robot = make_2dof_robot();
        let jac = Jacobian::compute(&robot);
        assert_eq!(jac.j.len(), 6 * robot.dof());
    }
    #[test]
    fn test_jacobian_manipulability_nonneg() {
        let robot = make_2dof_robot();
        let jac = Jacobian::compute(&robot);
        assert!(jac.manipulability() >= 0.0);
    }
    #[test]
    fn test_damped_pseudoinverse_size() {
        let robot = make_2dof_robot();
        let jac = Jacobian::compute(&robot);
        let jp = jac.damped_pseudoinverse(0.01);
        assert_eq!(jp.len(), robot.dof() * 6);
    }
    #[test]
    fn test_ik_solver_converges() {
        let mut robot = make_2dof_robot();
        robot.set_joint_angles(&[0.3, -0.5]);
        let target = [1.5, 0.5, 0.0];
        let ik = InverseKinematics::new();
        let converged = ik.solve_position(&mut robot, target);
        let pos = robot.end_effector_position();
        let err = ((pos[0] - target[0]).powi(2)
            + (pos[1] - target[1]).powi(2)
            + (pos[2] - target[2]).powi(2))
        .sqrt();
        assert!(converged || err < 0.1, "err={err}");
    }
    #[test]
    fn test_cubic_trajectory_endpoints() {
        let q0 = vec![0.0, 0.0];
        let qf = vec![1.0, -1.0];
        let traj = TrajectoryPlanning::new(q0.clone(), qf.clone(), 1.0);
        let start = traj.cubic_at(0.0);
        let end = traj.cubic_at(1.0);
        for (s, &qs) in start.iter().zip(q0.iter()) {
            assert!((s - qs).abs() < 1e-10);
        }
        for (e, &qg) in end.iter().zip(qf.iter()) {
            assert!((e - qg).abs() < 1e-10);
        }
    }
    #[test]
    fn test_quintic_trajectory_endpoints() {
        let q0 = vec![0.0];
        let qf = vec![2.0];
        let traj = TrajectoryPlanning::new(q0.clone(), qf.clone(), 2.0);
        assert!((traj.quintic_at(0.0)[0]).abs() < 1e-10);
        assert!((traj.quintic_at(2.0)[0] - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_cubic_velocity_zero_at_endpoints() {
        let traj = TrajectoryPlanning::new(vec![0.0], vec![1.0], 1.0);
        let v0 = traj.cubic_velocity_at(0.0);
        let vf = traj.cubic_velocity_at(1.0);
        assert!(v0[0].abs() < 1e-10, "v0={}", v0[0]);
        assert!(vf[0].abs() < 1e-10, "vf={}", vf[0]);
    }
    #[test]
    fn test_cartesian_linear_endpoints() {
        let p0 = [0.0, 0.0, 0.0];
        let pf = [1.0, 2.0, 3.0];
        let traj = CartesianTrajectory::linear(p0, pf, [0.0; 3], [0.0; 3], 1.0);
        let (p_start, _) = traj.evaluate(0.0);
        let (p_end, _) = traj.evaluate(1.0);
        for k in 0..3 {
            assert!((p_start[k] - p0[k]).abs() < 1e-10);
            assert!((p_end[k] - pf[k]).abs() < 1e-10);
        }
    }
    #[test]
    fn test_gravity_torques_size() {
        let robot = make_2dof_robot();
        let ff = TorqueFeedforward::new(vec![1.0, 1.0], vec![0.5, 0.5], [0.0, -9.81, 0.0]);
        let taus = ff.gravity_torques(&robot);
        assert_eq!(taus.len(), robot.dof());
    }
    #[test]
    fn test_impedance_wrench_zero_at_desired() {
        let mut ctrl = ImpedanceControl::new([100.0; 3], [10.0; 3]);
        ctrl.x_desired = [1.0, 2.0, 3.0];
        ctrl.xdot_desired = [0.0; 3];
        let w = ctrl.wrench([1.0, 2.0, 3.0], [0.0; 3]);
        for c in w {
            assert!(c.abs() < 1e-10);
        }
    }
    #[test]
    fn test_impedance_wrench_proportional_to_stiffness() {
        let ctrl1 = ImpedanceControl::new([100.0; 3], [0.0; 3]);
        let ctrl2 = ImpedanceControl::new([200.0; 3], [0.0; 3]);
        let w1 = ctrl1.wrench([0.0; 3], [0.0; 3]);
        let w2 = ctrl2.wrench([0.0; 3], [0.0; 3]);
        for k in 0..3 {
            assert!((w2[k] - 2.0 * w1[k]).abs() < 1e-10, "k={k}");
        }
    }
    #[test]
    fn test_repulsive_potential_near_obstacle() {
        let mut ca = CollisionAvoidance::new(1.0, 1.0);
        ca.add_obstacle([0.0; 3]);
        let u = ca.repulsive_potential([0.1, 0.0, 0.0]);
        assert!(u > 0.0, "u={u}");
    }
    #[test]
    fn test_repulsive_potential_outside_influence() {
        let mut ca = CollisionAvoidance::new(1.0, 1.0);
        ca.add_obstacle([0.0; 3]);
        let u = ca.repulsive_potential([2.0, 0.0, 0.0]);
        assert!(u.abs() < 1e-10, "u={u}");
    }
    #[test]
    fn test_attractive_force_direction() {
        let ca = CollisionAvoidance::new(1.0, 1.0);
        let q = [0.0; 3];
        let goal = [1.0, 0.0, 0.0];
        let f = ca.attractive_force(q, goal, 1.0);
        assert!(f[0] > 0.0, "Should point toward goal");
    }
    #[test]
    fn test_euler_rotation_roundtrip() {
        let roll = 0.3;
        let pitch = 0.1;
        let yaw = -0.5;
        let r = euler_to_rotation(roll, pitch, yaw);
        let euler = rotation_to_euler(r);
        assert!((euler[0] - roll).abs() < 1e-10, "roll={}", euler[0]);
        assert!((euler[1] - pitch).abs() < 1e-10, "pitch={}", euler[1]);
        assert!((euler[2] - yaw).abs() < 1e-10, "yaw={}", euler[2]);
    }
    #[test]
    fn test_adjoint_matrix_size() {
        let t = mat4_identity();
        let ad = adjoint_matrix(t);
        assert_eq!(ad.len(), 36);
    }
    #[test]
    fn test_rotation_from_axis_angle_zero() {
        let r = rotation_from_axis_angle([0.0, 0.0, 1.0], 0.0);
        assert!((r[0][0] - 1.0).abs() < 1e-10);
        assert!((r[1][1] - 1.0).abs() < 1e-10);
        assert!((r[2][2] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_adaptive_controller_torque_size() {
        let mut ctrl = AdaptiveRobotControl::new(3, 1.0, 10.0, 5.0);
        let q = vec![0.0; 3];
        let qdot = vec![0.0; 3];
        let q_d = vec![1.0; 3];
        let qdot_d = vec![0.0; 3];
        let qdotdot_d = vec![0.0; 3];
        let taus = ctrl.update(&q, &qdot, &q_d, &qdot_d, &qdotdot_d, 0.01);
        assert_eq!(taus.len(), 3);
    }
    #[test]
    fn test_trapezoidal_profile_integration() {
        let traj = TrajectoryPlanning::new(vec![0.0], vec![1.0], 1.0);
        let dt = 0.001;
        let mut integral = 0.0;
        let mut t = 0.0;
        while t < 1.0 {
            integral += traj.trapezoidal_profile(t, 0.2) * dt;
            t += dt;
        }
        assert!((integral - 1.0).abs() < 0.01, "integral={integral}");
    }
    #[test]
    fn test_inertia_torques_scale() {
        let robot = make_2dof_robot();
        let ff = TorqueFeedforward::new(vec![1.0, 1.0], vec![0.5, 0.5], [0.0, -9.81, 0.0]);
        let t1 = ff.inertia_torques(&robot, &[1.0, 1.0]);
        let t2 = ff.inertia_torques(&robot, &[2.0, 2.0]);
        for k in 0..2 {
            assert!((t2[k] - 2.0 * t1[k]).abs() < 1e-10, "k={k}");
        }
    }
    #[test]
    fn test_dh_params_forward_kinematics() {
        let dh = DhParams::new(0.5, 0.3, 0.2, 0.7);
        let m1 = dh.forward_kinematics();
        let m2 = dh_matrix(0.5, 0.3, 0.2, 0.7);
        for k in 0..16 {
            assert!((m1[k] - m2[k]).abs() < 1e-12);
        }
    }
}
#[allow(dead_code)]
pub(super) fn dot_n(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}
#[allow(dead_code)]
pub(super) fn norm_n(v: &[f64]) -> f64 {
    dot_n(v, v).sqrt()
}
#[allow(dead_code)]
pub(super) fn scale_vec(v: &[f64], s: f64) -> Vec<f64> {
    v.iter().map(|x| x * s).collect()
}
#[allow(dead_code)]
pub(super) fn add_vec(a: &[f64], b: &[f64]) -> Vec<f64> {
    a.iter().zip(b.iter()).map(|(x, y)| x + y).collect()
}
#[allow(dead_code)]
pub(super) fn sub_vec(a: &[f64], b: &[f64]) -> Vec<f64> {
    a.iter().zip(b.iter()).map(|(x, y)| x - y).collect()
}
#[allow(dead_code)]
pub(super) fn clamp_vec(v: &[f64], lo: f64, hi: f64) -> Vec<f64> {
    v.iter().map(|x| x.max(lo).min(hi)).collect()
}
/// Simple Gauss-Seidel solver for small dense systems A*x = b.
#[allow(dead_code)]
pub(super) fn gauss_seidel_solve(a: &[Vec<f64>], b: &[f64], max_iter: usize) -> Vec<f64> {
    let n = b.len();
    let mut x = vec![0.0; n];
    for _ in 0..max_iter {
        for i in 0..n {
            let mut sum = b[i];
            for j in 0..n {
                if j != i {
                    sum -= a[i][j] * x[j];
                }
            }
            if a[i][i].abs() > 1e-15 {
                x[i] = sum / a[i][i];
            }
        }
    }
    x
}
/// Invert a dense square matrix (n×n, row-major) using Gauss-Jordan.
///
/// Returns the inverse as a flat Vec of length n×n, or identity on failure.
#[allow(dead_code)]
pub(super) fn invert_dense_matrix(a: &[f64], n: usize) -> Vec<f64> {
    let mut aug = vec![0.0; n * 2 * n];
    for i in 0..n {
        for j in 0..n {
            aug[i * 2 * n + j] = a[i * n + j];
        }
        aug[i * 2 * n + n + i] = 1.0;
    }
    for col in 0..n {
        let mut best = col;
        for row in (col + 1)..n {
            if aug[row * 2 * n + col].abs() > aug[best * 2 * n + col].abs() {
                best = row;
            }
        }
        for k in 0..2 * n {
            aug.swap(col * 2 * n + k, best * 2 * n + k);
        }
        let pivot = aug[col * 2 * n + col];
        if pivot.abs() < 1e-15 {
            let mut id = vec![0.0; n * n];
            for i in 0..n {
                id[i * n + i] = 1.0;
            }
            return id;
        }
        for k in 0..2 * n {
            aug[col * 2 * n + k] /= pivot;
        }
        for row in 0..n {
            if row != col {
                let factor = aug[row * 2 * n + col];
                for k in 0..2 * n {
                    let val = aug[col * 2 * n + k];
                    aug[row * 2 * n + k] -= factor * val;
                }
            }
        }
    }
    let mut inv = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            inv[i * n + j] = aug[i * 2 * n + n + j];
        }
    }
    inv
}
#[cfg(test)]
mod expanded_robot_tests {
    use super::*;

    use crate::robot_control::BroydenJacobianUpdate;
    use crate::robot_control::CartesianImpedanceController;

    use crate::robot_control::ContactPhase;
    use crate::robot_control::ContactTransitionDetector;

    use crate::robot_control::FtSensorData;
    use crate::robot_control::FtSensorFilter;

    use crate::robot_control::JointMpc;
    use crate::robot_control::JointTrajectory;
    use crate::robot_control::NullSpaceProjector;
    use crate::robot_control::PassivityBasedController;
    use crate::robot_control::PayloadIdentifier;
    use crate::robot_control::QuinticTrajectory;
    use crate::robot_control::RneaLink;
    use crate::robot_control::RneaSolver;
    use crate::robot_control::RobotSafetyMonitor;
    use crate::robot_control::SafetyLevel;
    use crate::robot_control::SingularityRobustController;
    use crate::robot_control::StiffnessRenderer;

    use crate::robot_control::TorqueLimits;

    use crate::robot_control::VelocityLimits;
    use crate::robot_control::WorkspaceAnalysis;
    #[test]
    fn test_torque_limits_clamp() {
        let lim = TorqueLimits::new(vec![10.0, 5.0]);
        let clamped = lim.clamp(&[15.0, -8.0]);
        assert!((clamped[0] - 10.0).abs() < 1e-10);
        assert!((clamped[1] + 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_torque_limits_scale() {
        let lim = TorqueLimits::uniform(2, 10.0);
        let scaled = lim.scale_to_satisfy(&[20.0, 5.0]);
        assert!((scaled[0] - 10.0).abs() < 1e-10);
        assert!((scaled[1] - 2.5).abs() < 1e-10);
    }
    #[test]
    fn test_torque_limits_satisfied() {
        let lim = TorqueLimits::uniform(2, 10.0);
        assert!(lim.satisfied(&[5.0, -5.0], 1e-6));
        assert!(!lim.satisfied(&[12.0, 5.0], 1e-6));
    }
    #[test]
    fn test_velocity_limits_clamp() {
        let lim = VelocityLimits::uniform(3, 2.0);
        let clamped = lim.clamp(&[3.0, 1.0, -5.0]);
        assert!((clamped[0] - 2.0).abs() < 1e-10);
        assert!((clamped[1] - 1.0).abs() < 1e-10);
        assert!((clamped[2] + 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_rnea_link_axis_inertia() {
        let link = RneaLink::cylindrical(1.0, 0.5, 0.05);
        let j = link.axis_inertia();
        assert!((j - 0.5 * 1.0 * 0.05 * 0.05).abs() < 1e-10, "j={j}");
    }
    #[test]
    fn test_rnea_compute_torques() {
        let links = vec![
            RneaLink::cylindrical(1.0, 0.4, 0.04),
            RneaLink::cylindrical(0.5, 0.3, 0.03),
        ];
        let solver = RneaSolver::new(links, [0.0, -9.81, 0.0]);
        let tau = solver.compute_torques(&[0.0, 0.5], &[0.1, 0.1], &[0.2, 0.2]);
        assert_eq!(tau.len(), 2);
    }
    #[test]
    fn test_gauss_seidel_identity() {
        let a = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let b = vec![3.0, 5.0];
        let x = gauss_seidel_solve(&a, &b, 10);
        assert!((x[0] - 3.0).abs() < 1e-6, "x[0]={}", x[0]);
        assert!((x[1] - 5.0).abs() < 1e-6, "x[1]={}", x[1]);
    }
    #[test]
    fn test_quintic_boundary() {
        let traj = QuinticTrajectory::new(0.0, 1.0, 1.0);
        assert!((traj.position(0.0) - 0.0).abs() < 1e-10);
        assert!((traj.position(1.0) - 1.0).abs() < 1e-10);
        assert!(traj.velocity(0.0).abs() < 1e-10);
        assert!(traj.velocity(1.0).abs() < 1e-10);
        assert!(traj.acceleration(0.0).abs() < 1e-10);
        assert!(traj.acceleration(1.0).abs() < 1e-10);
    }
    #[test]
    fn test_quintic_midpoint() {
        let traj = QuinticTrajectory::new(0.0, 2.0, 2.0);
        assert!((traj.position(1.0) - 1.0).abs() < 1e-6);
    }
    #[test]
    fn test_joint_trajectory_positions() {
        let traj = JointTrajectory::new(&[0.0, 0.0], &[1.0, 2.0], 1.0);
        let p0 = traj.position(0.0);
        let p1 = traj.position(1.0);
        assert!((p0[0]).abs() < 1e-10);
        assert!((p1[0] - 1.0).abs() < 1e-10);
        assert!((p1[1] - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_workspace_analysis() {
        let positions = vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.5]];
        let manip = vec![0.5, 0.3, 0.8];
        let ws = WorkspaceAnalysis::from_samples(&positions, &manip);
        assert!(ws.max_reach >= 1.0);
        assert!(ws.n_samples == 3);
        assert!(ws.contains([0.5, 0.5, 0.0]));
    }
    #[test]
    fn test_passivity_damping_torque() {
        let ctrl = PassivityBasedController::new(vec![2.0, 3.0]);
        let tau = ctrl.damping_torque(&[1.0, -1.0]);
        assert!((tau[0] + 2.0).abs() < 1e-10, "tau[0]={}", tau[0]);
        assert!((tau[1] - 3.0).abs() < 1e-10, "tau[1]={}", tau[1]);
    }
    #[test]
    fn test_passivity_condition() {
        let ctrl = PassivityBasedController::new(vec![1.0, 1.0]);
        assert!(ctrl.is_passive(&[1.0, 1.0]));
    }
    #[test]
    fn test_cartesian_impedance_wrench() {
        let k = [100.0; 6];
        let b = [10.0; 6];
        let m = [1.0; 6];
        let ctrl = CartesianImpedanceController::new(k, b, m);
        let pos_err = [0.01, 0.0, 0.0, 0.0, 0.0, 0.0];
        let vel_err = [0.0; 6];
        let acc_des = [0.0; 6];
        let w = ctrl.compute_wrench(&pos_err, &vel_err, &acc_des);
        assert!((w[0] - 1.0).abs() < 1e-10, "w[0]={}", w[0]);
    }
    #[test]
    fn test_impedance_damping_ratio() {
        let k = [100.0; 6];
        let m = [1.0; 6];
        let ctrl = CartesianImpedanceController::critically_damped(k, m);
        let zeta = ctrl.damping_ratio(0);
        assert!((zeta - 1.0).abs() < 1e-6, "zeta={zeta}");
    }
    #[test]
    fn test_null_space_projector_no_panic() {
        let jac = vec![vec![1.0, 0.0, 0.0], vec![0.0, 1.0, 0.0]];
        let proj = NullSpaceProjector::new(jac, 0.01);
        let result = proj.project(&[1.0, 1.0, 1.0]);
        assert_eq!(result.len(), 3);
    }
    #[test]
    fn test_joint_limit_gradient() {
        let q = &[0.5];
        let q_min = &[0.0];
        let q_max = &[1.0];
        let grad = NullSpaceProjector::joint_limit_gradient(q, q_min, q_max);
        assert!(grad[0].abs() < 1e-10);
    }
    #[test]
    fn test_invert_identity() {
        let id = [1.0_f64, 0.0, 0.0, 1.0];
        let inv = invert_dense_matrix(&id, 2);
        assert!((inv[0] - 1.0).abs() < 1e-10);
        assert!((inv[3] - 1.0).abs() < 1e-10);
        assert!(inv[1].abs() < 1e-10);
    }
    #[test]
    fn test_broyden_update_no_nan() {
        let jac = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let mut bj = BroydenJacobianUpdate::new(jac, 0.99);
        bj.update(&[0.01, 0.0], &[0.01, 0.0]);
        let flat = bj.flat();
        assert!(flat.iter().all(|v| v.is_finite()), "Jacobian has NaN/inf");
    }
    #[test]
    fn test_singularity_robust_damping() {
        let ctrl = SingularityRobustController::new(0.1, 0.5, 0.05);
        let d_far = ctrl.adaptive_damping(1.0);
        assert!(d_far.abs() < 1e-10, "d_far={d_far}");
        let d_near = ctrl.adaptive_damping(0.0);
        assert!((d_near - 0.5).abs() < 1e-10, "d_near={d_near}");
    }
    #[test]
    fn test_mpc_rollout_length() {
        let mpc = JointMpc::new(5, 0.01, 100.0, 1.0, 0.1, 10.0, 0.1);
        let q_ref: Vec<f64> = vec![0.5; 20];
        let traj = mpc.rollout(0.0, 0.0, &q_ref, 10);
        assert_eq!(traj.len(), 11);
    }
    #[test]
    fn test_ft_sensor_force_magnitude() {
        let ft = FtSensorData::from_wrench([3.0, 4.0, 0.0, 0.0, 0.0, 0.0]);
        assert!((ft.force_magnitude() - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_ft_filter_convergence() {
        let mut filter = FtSensorFilter::new(10.0, 1000.0);
        let raw = [1.0; 6];
        for _ in 0..1000 {
            filter.update(raw);
        }
        let state = filter.filtered();
        assert!((state[0] - 1.0).abs() < 0.01, "state[0]={}", state[0]);
    }
    #[test]
    fn test_contact_detector_transitions() {
        let mut det = ContactTransitionDetector::new(1.0, 5.0, 20.0, 0.2);
        det.debounce_count = 1;
        det.update(0.5);
        assert_eq!(*det.update(0.5), ContactPhase::FreeMotion);
        for _ in 0..5 {
            det.update(2.0);
        }
        assert_eq!(det.phase, ContactPhase::LightContact);
        for _ in 0..5 {
            det.update(8.0);
        }
        assert_eq!(det.phase, ContactPhase::StableContact);
    }
    #[test]
    fn test_payload_identifier_zero_payload() {
        let mut ident = PayloadIdentifier::new(1000.0, 0.99);
        let reg = [[0.0_f64; 4]; 3];
        let force = [0.0_f64; 3];
        for _ in 0..10 {
            ident.update(&reg, &force);
        }
        assert!(ident.mass().abs() < 100.0);
    }
    #[test]
    fn test_safety_monitor_normal() {
        let mut mon = RobotSafetyMonitor::new(
            vec![[-1.0, 1.0], [-1.0, 1.0]],
            vec![5.0, 5.0],
            vec![20.0, 20.0],
        );
        let level = mon.check(&[0.0, 0.0], &[0.1, 0.1], &[1.0, 1.0]);
        assert_eq!(*level, SafetyLevel::Normal);
    }
    #[test]
    fn test_safety_monitor_estop() {
        let mut mon = RobotSafetyMonitor::new(vec![[-1.0, 1.0]], vec![5.0], vec![20.0]);
        let level = mon.check(&[2.0], &[0.0], &[0.0]);
        assert_eq!(*level, SafetyLevel::EmergencyStop);
    }
    #[test]
    fn test_stiffness_renderer_contact() {
        let wall = StiffnessRenderer::new(0.0, 1.0, 1000.0, 50.0);
        assert!(!wall.in_contact(-0.01));
        assert!(wall.in_contact(0.01));
        let f = wall.force(0.01, 0.0);
        assert!(f < 0.0, "Contact force should push back, f={f}");
    }
    #[test]
    fn test_stiffness_renderer_no_contact() {
        let wall = StiffnessRenderer::new(0.0, 1.0, 1000.0, 50.0);
        let f = wall.force(-0.1, 0.0);
        assert!(f.abs() < 1e-10, "No force outside wall, f={f}");
    }
    #[test]
    fn test_joint_trajectory_sample() {
        let traj = JointTrajectory::new(&[0.0], &[1.0], 2.0);
        let samples = traj.sample(10);
        assert_eq!(samples.len(), 11);
    }
}
