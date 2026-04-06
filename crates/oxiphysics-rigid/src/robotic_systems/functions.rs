//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
#[inline]
pub(super) fn vec3_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
#[inline]
pub(super) fn vec3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
#[inline]
pub(super) fn vec3_scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
#[inline]
pub(super) fn vec3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
#[inline]
pub(super) fn vec3_norm(a: [f64; 3]) -> f64 {
    vec3_dot(a, a).sqrt()
}
#[inline]
pub(super) fn vec3_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
/// Multiply two 4×4 homogeneous matrices (row-major).
pub(super) fn mat4_mul(a: &[[f64; 4]; 4], b: &[[f64; 4]; 4]) -> [[f64; 4]; 4] {
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
/// 4×4 identity matrix.
pub(super) fn mat4_eye() -> [[f64; 4]; 4] {
    let mut m = [[0.0f64; 4]; 4];
    m[0][0] = 1.0;
    m[1][1] = 1.0;
    m[2][2] = 1.0;
    m[3][3] = 1.0;
    m
}
/// Build DH homogeneous transform: T = Rot_z(theta) * Trans_z(d) * Trans_x(a) * Rot_x(alpha).
pub(super) fn dh_transform(theta: f64, d: f64, a: f64, alpha: f64) -> [[f64; 4]; 4] {
    let ct = theta.cos();
    let st = theta.sin();
    let ca = alpha.cos();
    let sa = alpha.sin();
    [
        [ct, -st * ca, st * sa, a * ct],
        [st, ct * ca, -ct * sa, a * st],
        [0.0, sa, ca, d],
        [0.0, 0.0, 0.0, 1.0],
    ]
}
/// Multiply an NxM matrix (stored as `Vec<Vec`f64`>`) by a vector.
pub(super) fn mat_vec_mul(m: &[Vec<f64>], v: &[f64]) -> Vec<f64> {
    let rows = m.len();
    let mut result = vec![0.0; rows];
    for (i, row) in m.iter().enumerate() {
        for (j, &vj) in v.iter().enumerate() {
            if j < row.len() {
                result[i] += row[j] * vj;
            }
        }
    }
    result
}
/// Transpose of an NxM matrix.
pub(super) fn mat_transpose(m: &[Vec<f64>]) -> Vec<Vec<f64>> {
    if m.is_empty() {
        return Vec::new();
    }
    let rows = m.len();
    let cols = m[0].len();
    let mut t = vec![vec![0.0; rows]; cols];
    for i in 0..rows {
        for j in 0..cols {
            t[j][i] = m[i][j];
        }
    }
    t
}
/// Matrix multiply A (n×k) × B (k×m) → C (n×m).
pub(super) fn mat_mul_dyn(a: &[Vec<f64>], b: &[Vec<f64>]) -> Vec<Vec<f64>> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let n = a.len();
    let k = b.len();
    let m = b[0].len();
    let mut c = vec![vec![0.0; m]; n];
    for i in 0..n {
        for j in 0..m {
            for l in 0..k {
                if l < a[i].len() && j < b[l].len() {
                    c[i][j] += a[i][l] * b[l][j];
                }
            }
        }
    }
    c
}
/// Add scalar * identity to a square matrix (for damped LS).
pub(super) fn mat_add_lambda_eye(m: &mut [Vec<f64>], lambda: f64) {
    for (i, row) in m.iter_mut().enumerate() {
        if i < row.len() {
            row[i] += lambda;
        }
    }
}
/// Solve a small NxN linear system Ax = b using Gaussian elimination.
/// Returns None if singular (det ≈ 0).
pub(super) fn solve_linear(a: &[Vec<f64>], b: &[f64]) -> Option<Vec<f64>> {
    let n = a.len();
    if n == 0 {
        return Some(Vec::new());
    }
    let mut aug: Vec<Vec<f64>> = a
        .iter()
        .zip(b.iter())
        .map(|(row, &rhs)| {
            let mut r = row.clone();
            r.push(rhs);
            r
        })
        .collect();
    for col in 0..n {
        let pivot = (col..n).max_by(|&i, &j| {
            aug[i][col]
                .abs()
                .partial_cmp(&aug[j][col].abs())
                .unwrap_or(std::cmp::Ordering::Equal)
        })?;
        aug.swap(col, pivot);
        let diag = aug[col][col];
        if diag.abs() < 1e-14 {
            return None;
        }
        for k in col..=n {
            aug[col][k] /= diag;
        }
        for r in 0..n {
            if r != col {
                let factor = aug[r][col];
                for k in col..=n {
                    let val = aug[col][k];
                    aug[r][k] -= factor * val;
                }
            }
        }
    }
    Some((0..n).map(|i| aug[i][n]).collect())
}
/// Compute determinant of a 3×3 matrix stored as `Vec<Vec`f64`>`.
#[allow(dead_code)]
pub(super) fn det3x3(m: &[Vec<f64>]) -> f64 {
    if m.len() < 3 || m[0].len() < 3 {
        return 0.0;
    }
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}
/// Compute determinant of an n×n matrix via Gaussian elimination.
pub(super) fn det_nxn(m: &[Vec<f64>]) -> f64 {
    let n = m.len();
    if n == 0 {
        return 1.0;
    }
    if n == 1 {
        return if m[0].is_empty() { 0.0 } else { m[0][0] };
    }
    let mut a: Vec<Vec<f64>> = m.to_vec();
    let mut det = 1.0f64;
    for col in 0..n {
        let pivot_row = (col..n).max_by(|&i, &j| {
            a[i].get(col)
                .unwrap_or(&0.0)
                .abs()
                .partial_cmp(&a[j].get(col).unwrap_or(&0.0).abs())
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        if let Some(pivot) = pivot_row
            && pivot != col
        {
            a.swap(col, pivot);
            det = -det;
        }
        let diag = *a[col].get(col).unwrap_or(&0.0);
        if diag.abs() < 1e-14 {
            return 0.0;
        }
        det *= diag;
        for k in col..n {
            a[col][k] /= diag;
        }
        for r in 0..n {
            if r != col {
                let factor = *a[r].get(col).unwrap_or(&0.0);
                if factor.abs() < 1e-300 {
                    continue;
                }
                for k in col..n {
                    let val = a[col][k];
                    a[r][k] -= factor * val;
                }
            }
        }
    }
    det
}
/// Select active (non-zero) position rows from the 6×n Jacobian matrix.
///
/// Returns at most `n` rows that have meaningful (non-zero) norms, so the
/// resulting system is well-conditioned.
#[allow(dead_code)]
pub(super) fn active_position_jacobian(j_mat: &[Vec<f64>], n: usize) -> Vec<Vec<f64>> {
    let pos_rows: Vec<Vec<f64>> = j_mat[..3.min(j_mat.len())]
        .iter()
        .filter(|row| row.iter().map(|x| x * x).sum::<f64>().sqrt() > 1e-12)
        .cloned()
        .collect();
    pos_rows[..n.min(pos_rows.len())].to_vec()
}
/// Extract the error components corresponding to the active Jacobian rows.
#[allow(dead_code)]
pub(super) fn active_error(err: &[f64; 3], active_j: &[Vec<f64>], full_j: &[Vec<f64>]) -> Vec<f64> {
    let pos_rows: Vec<(usize, Vec<f64>)> = full_j[..3.min(full_j.len())]
        .iter()
        .enumerate()
        .filter(|(_, row)| row.iter().map(|x| x * x).sum::<f64>().sqrt() > 1e-12)
        .map(|(i, row)| (i, row.clone()))
        .collect();
    let n_active = active_j.len();
    pos_rows[..n_active]
        .iter()
        .map(|(i, _)| if *i < 3 { err[*i] } else { 0.0 })
        .collect()
}
/// Compute a tangent vector orthogonal to `n`.
pub(super) fn tangent_vec(n: [f64; 3]) -> [f64; 3] {
    let t = if n[0].abs() < 0.9 {
        vec3_cross(n, [1.0, 0.0, 0.0])
    } else {
        vec3_cross(n, [0.0, 1.0, 0.0])
    };
    let len = vec3_norm(t);
    if len > 1e-10 {
        vec3_scale(t, 1.0 / len)
    } else {
        [0.0, 1.0, 0.0]
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::robotic_systems::CollisionAvoidance;
    use crate::robotic_systems::ContactPoint;
    use crate::robotic_systems::DhParam;
    use crate::robotic_systems::ForceControl;
    use crate::robotic_systems::GraspPlanning;
    use crate::robotic_systems::GraspStability;
    use crate::robotic_systems::ImpedanceParams;
    use crate::robotic_systems::InverseKinematics;
    use crate::robotic_systems::JacobianMatrix;
    use crate::robotic_systems::JointType;
    use crate::robotic_systems::KinematicChain;
    use crate::robotic_systems::LinkInertia;
    use crate::robotic_systems::RobotDynamics;
    use crate::robotic_systems::SerialManipulator;
    use crate::robotic_systems::SphereObstacle;
    use crate::robotic_systems::TrajectoryPlanning;
    use crate::robotic_systems::Waypoint;
    use crate::robotic_systems::WorkspaceAnalysis;
    use std::f64::consts::PI;
    /// Build a 2-DOF planar RR robot (two revolute joints of unit length).
    fn make_rr_robot(q0: f64, q1: f64) -> KinematicChain {
        let params = vec![
            DhParam::new(1.0, 0.0, 0.0, 0.0),
            DhParam::new(1.0, 0.0, 0.0, 0.0),
        ];
        let types = vec![JointType::Revolute, JointType::Revolute];
        KinematicChain::new(params, types, vec![q0, q1])
    }
    /// Build a 3-DOF robot.
    fn make_rrr_robot(q: Vec<f64>) -> KinematicChain {
        let params = vec![
            DhParam::new(1.0, 0.0, 0.0, 0.0),
            DhParam::new(1.0, 0.0, 0.0, 0.0),
            DhParam::new(0.5, 0.0, 0.0, 0.0),
        ];
        let types = vec![
            JointType::Revolute,
            JointType::Revolute,
            JointType::Revolute,
        ];
        KinematicChain::new(params, types, q)
    }
    #[test]
    fn dh_identity_at_zero() {
        let t = dh_transform(0.0, 0.0, 0.0, 0.0);
        assert!((t[0][0] - 1.0).abs() < 1e-10);
        assert!((t[1][1] - 1.0).abs() < 1e-10);
        assert!((t[2][2] - 1.0).abs() < 1e-10);
        assert!((t[3][3] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn dh_pure_rotation_z() {
        let t = dh_transform(PI / 2.0, 0.0, 0.0, 0.0);
        assert!((t[0][0]).abs() < 1e-10);
        assert!((t[0][1] + 1.0).abs() < 1e-10);
        assert!((t[1][0] - 1.0).abs() < 1e-10);
        assert!((t[1][1]).abs() < 1e-10);
    }
    #[test]
    fn dh_pure_translation_a() {
        let t = dh_transform(0.0, 0.0, 1.0, 0.0);
        assert!((t[0][3] - 1.0).abs() < 1e-10);
        assert!((t[1][3]).abs() < 1e-10);
        assert!((t[2][3]).abs() < 1e-10);
    }
    #[test]
    fn dh_pure_translation_d() {
        let t = dh_transform(0.0, 2.0, 0.0, 0.0);
        assert!((t[0][3]).abs() < 1e-10);
        assert!((t[1][3]).abs() < 1e-10);
        assert!((t[2][3] - 2.0).abs() < 1e-10);
    }
    #[test]
    fn rr_robot_zero_config_tip_at_2_0() {
        let chain = make_rr_robot(0.0, 0.0);
        let pos = chain.end_effector_position();
        assert!((pos[0] - 2.0).abs() < 1e-10, "x={}", pos[0]);
        assert!((pos[1]).abs() < 1e-10, "y={}", pos[1]);
    }
    #[test]
    fn rr_robot_q0_pi_half_tip_y() {
        let chain = make_rr_robot(PI / 2.0, 0.0);
        let pos = chain.end_effector_position();
        assert!(pos[0].abs() < 1e-9, "x={}", pos[0]);
        assert!((pos[1] - 2.0).abs() < 1e-9, "y={}", pos[1]);
    }
    #[test]
    fn rr_robot_folded_q1_pi_tip_at_origin() {
        let chain = make_rr_robot(0.0, PI);
        let pos = chain.end_effector_position();
        assert!(pos[0].abs() < 1e-9, "x={}", pos[0]);
    }
    #[test]
    fn forward_kinematics_length_equals_dof() {
        let chain = make_rrr_robot(vec![0.0, 0.0, 0.0]);
        let transforms = chain.forward_kinematics();
        assert_eq!(transforms.len(), 3);
    }
    #[test]
    fn dof_returns_correct_count() {
        let chain = make_rrr_robot(vec![0.0, 0.1, 0.2]);
        assert_eq!(chain.dof(), 3);
    }
    #[test]
    fn jacobian_size_is_6n() {
        let chain = make_rr_robot(0.0, 0.0);
        let j = chain.jacobian();
        assert_eq!(j.len(), 12);
    }
    #[test]
    fn jacobian_matrix_has_6_rows() {
        let chain = make_rrr_robot(vec![0.0, 0.0, 0.0]);
        let j = chain.jacobian_matrix();
        assert_eq!(j.len(), 6);
        assert_eq!(j[0].len(), 3);
    }
    #[test]
    fn jacobian_nonzero_at_standard_config() {
        let chain = make_rr_robot(PI / 4.0, PI / 4.0);
        let j = chain.jacobian();
        let has_nonzero = j.iter().any(|&v| v.abs() > 1e-10);
        assert!(has_nonzero, "Jacobian should be nonzero");
    }
    #[test]
    fn manipulability_zero_at_singularity() {
        let chain = make_rr_robot(0.0, 0.0);
        let m = chain.manipulability();
        assert!(m < 1e-8, "manipulability at singularity: {m}");
    }
    #[test]
    fn manipulability_nonzero_at_generic() {
        let chain = make_rr_robot(PI / 4.0, PI / 3.0);
        let m = chain.manipulability();
        assert!(m > 0.0, "manipulability should be positive");
    }
    #[test]
    fn clamp_joints_enforces_limits() {
        let mut chain = make_rr_robot(0.0, 0.0);
        chain.q[0] = 5.0;
        chain.clamp_joints();
        assert!(chain.q[0] <= PI);
    }
    #[test]
    fn within_limits_fails_if_violated() {
        let mut chain = make_rr_robot(0.0, 0.0);
        chain.q[1] = 4.0;
        assert!(!chain.within_limits());
    }
    #[test]
    fn within_limits_passes_at_zero() {
        let chain = make_rr_robot(0.0, 0.0);
        assert!(chain.within_limits());
    }
    #[test]
    fn ik_pseudoinverse_reaches_reachable_target() {
        let mut chain = make_rr_robot(0.1, 0.1);
        let ik = InverseKinematics::with_params(500, 1e-3, 0.1, 0.05);
        let target = [1.0, 0.5, 0.0];
        let (_q, converged) = ik.solve_pseudoinverse(&mut chain, target);
        let pos = chain.end_effector_position();
        let err = ((pos[0] - target[0]).powi(2) + (pos[1] - target[1]).powi(2)).sqrt();
        assert!(
            converged || err < 0.1,
            "IK should converge or be close: err={err}, converged={converged}"
        );
    }
    #[test]
    fn ik_damped_ls_reaches_reachable_target() {
        let mut chain = make_rr_robot(0.1, 0.1);
        let ik = InverseKinematics::new();
        let target = [1.0, 0.5, 0.0];
        let (_q, converged) = ik.solve_damped_ls(&mut chain, target);
        assert!(converged, "Damped LS IK should converge");
    }
    #[test]
    fn ik_solution_reduces_error() {
        let mut chain = make_rr_robot(0.1, 0.1);
        let ik = InverseKinematics::new();
        let target = [1.2, 0.3, 0.0];
        let initial_pos = chain.end_effector_position();
        let initial_err = vec3_norm(vec3_sub(target, initial_pos));
        ik.solve_damped_ls(&mut chain, target);
        let final_pos = chain.end_effector_position();
        let final_err = vec3_norm(vec3_sub(target, final_pos));
        assert!(final_err < initial_err, "IK should reduce error");
    }
    #[test]
    fn ik_null_space_projection_size() {
        let chain = make_rrr_robot(vec![0.2, 0.3, 0.1]);
        let ik = InverseKinematics::new();
        let grad = vec![0.1, -0.1, 0.05];
        let proj = ik.null_space_projection(&chain, &grad);
        assert_eq!(
            proj.len(),
            3,
            "Null space projection should have DOF entries"
        );
    }
    fn make_dynamics(n: usize) -> RobotDynamics {
        let inertias = (0..n)
            .map(|_| LinkInertia::point_mass(1.0, [0.5, 0.0, 0.0]))
            .collect();
        RobotDynamics::new(inertias, [0.0, -9.81, 0.0])
    }
    #[test]
    fn gravity_torques_size_equals_dof() {
        let chain = make_rr_robot(0.0, 0.0);
        let dyn_model = make_dynamics(2);
        let tau = dyn_model.gravity_torques(&chain);
        assert_eq!(tau.len(), 2);
    }
    #[test]
    fn gravity_torques_nonzero_at_zero_config() {
        let chain = make_rr_robot(0.0, 0.0);
        let dyn_model = make_dynamics(2);
        let tau = dyn_model.gravity_torques(&chain);
        let nonzero = tau.iter().any(|&t| t.abs() > 1e-6);
        assert!(
            nonzero,
            "gravity torques should be nonzero at q=[0,0]: {:?}",
            tau
        );
        assert!(tau[0].abs() > 1e-6, "tau[0] should be nonzero: {}", tau[0]);
    }
    #[test]
    fn diagonal_mass_matrix_positive() {
        let chain = make_rr_robot(0.2, 0.3);
        let dyn_model = make_dynamics(2);
        let m = dyn_model.diagonal_mass_matrix(&chain);
        assert!(
            m.iter().all(|&v| v >= 0.0),
            "mass matrix diagonal should be non-negative"
        );
    }
    #[test]
    fn coriolis_zero_when_qdot_zero() {
        let chain = make_rr_robot(0.0, 0.0);
        let dyn_model = make_dynamics(2);
        let c = dyn_model.coriolis_torques(&chain, &[0.0, 0.0]);
        assert!(c.iter().all(|&v| v.abs() < 1e-10));
    }
    #[test]
    fn link_inertia_box_positive_values() {
        let l = LinkInertia::box_link(2.0, 0.1, 0.2, 0.3);
        assert!(l.inertia[0][0] > 0.0);
        assert!(l.inertia[1][1] > 0.0);
        assert!(l.inertia[2][2] > 0.0);
    }
    #[test]
    fn linear_interp_at_start() {
        let wps = vec![Waypoint::new(vec![0.0], 0.0), Waypoint::new(vec![1.0], 1.0)];
        let tp = TrajectoryPlanning::new(wps);
        let q = tp.linear_interp(0.0).unwrap();
        assert!((q[0] - 0.0).abs() < 1e-10);
    }
    #[test]
    fn linear_interp_at_end() {
        let wps = vec![Waypoint::new(vec![0.0], 0.0), Waypoint::new(vec![1.0], 1.0)];
        let tp = TrajectoryPlanning::new(wps);
        let q = tp.linear_interp(1.0).unwrap();
        assert!((q[0] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn linear_interp_midpoint() {
        let wps = vec![Waypoint::new(vec![0.0], 0.0), Waypoint::new(vec![2.0], 2.0)];
        let tp = TrajectoryPlanning::new(wps);
        let q = tp.linear_interp(1.0).unwrap();
        assert!((q[0] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn cubic_poly_boundary_conditions() {
        let (p0, v0) = TrajectoryPlanning::cubic_poly(0.0, 1.0, 0.0);
        let (pf, vf) = TrajectoryPlanning::cubic_poly(0.0, 1.0, 1.0);
        assert!((p0 - 0.0).abs() < 1e-10);
        assert!((pf - 1.0).abs() < 1e-10);
        assert!(v0.abs() < 1e-10, "initial velocity should be zero");
        assert!(vf.abs() < 1e-10, "final velocity should be zero");
    }
    #[test]
    fn quintic_poly_boundary_conditions() {
        let (p0, v0, a0) = TrajectoryPlanning::quintic_poly(0.0, 1.0, 0.0);
        let (pf, vf, af) = TrajectoryPlanning::quintic_poly(0.0, 1.0, 1.0);
        assert!((p0 - 0.0).abs() < 1e-10);
        assert!((pf - 1.0).abs() < 1e-10);
        assert!(v0.abs() < 1e-10);
        assert!(vf.abs() < 1e-10);
        assert!(a0.abs() < 1e-10);
        assert!(af.abs() < 1e-10);
    }
    #[test]
    fn trapezoidal_profile_at_endpoints() {
        let q0 = TrajectoryPlanning::trapezoidal_profile(0.0, 2.0, 0.0, 2.0, 0.5).unwrap();
        let qf = TrajectoryPlanning::trapezoidal_profile(0.0, 2.0, 2.0, 2.0, 0.5).unwrap();
        assert!((q0 - 0.0).abs() < 1e-10);
        assert!((qf - 2.0).abs() < 1e-8);
    }
    #[test]
    fn trapezoidal_profile_none_for_invalid_ta() {
        let result = TrajectoryPlanning::trapezoidal_profile(0.0, 1.0, 0.5, 1.0, 0.6);
        assert!(result.is_none());
    }
    #[test]
    fn minimum_jerk_boundary_conditions() {
        let (p0, v0, a0) = TrajectoryPlanning::minimum_jerk(0.0, 1.0, 0.0);
        let (pf, vf, af) = TrajectoryPlanning::minimum_jerk(0.0, 1.0, 1.0);
        assert!((p0 - 0.0).abs() < 1e-10);
        assert!((pf - 1.0).abs() < 1e-10);
        assert!(v0.abs() < 1e-10);
        assert!(vf.abs() < 1e-10);
        assert!(a0.abs() < 1e-10);
        assert!(af.abs() < 1e-10);
    }
    #[test]
    fn cubic_spline_midpoint() {
        let wps = vec![Waypoint::new(vec![0.0], 0.0), Waypoint::new(vec![2.0], 2.0)];
        let tp = TrajectoryPlanning::new(wps);
        let q = tp.cubic_spline(1.0).unwrap();
        assert!((q[0] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn arc_length_positive_for_moving_trajectory() {
        let wps = vec![
            Waypoint::new(vec![0.0, 0.0], 0.0),
            Waypoint::new(vec![1.0, 1.0], 1.0),
        ];
        let tp = TrajectoryPlanning::new(wps);
        let len = tp.arc_length(20);
        assert!(len > 0.0, "arc length should be positive");
    }
    #[test]
    fn no_collision_when_no_obstacles() {
        let ca = CollisionAvoidance::new(vec![], 0.1, 0.05);
        let chain = make_rr_robot(PI / 4.0, PI / 4.0);
        assert!(!ca.check_obstacle_collision(&chain));
    }
    #[test]
    fn collision_detected_with_obstacle_at_tip() {
        let chain = make_rr_robot(0.0, 0.0);
        let tip = chain.end_effector_position();
        let obs = SphereObstacle::new(tip, 0.5);
        let ca = CollisionAvoidance::new(vec![obs], 0.1, 0.05);
        assert!(ca.check_obstacle_collision(&chain));
    }
    #[test]
    fn no_self_collision_for_extended_robot() {
        let ca = CollisionAvoidance::new(vec![], 0.05, 0.01);
        let chain = make_rr_robot(0.0, 0.0);
        assert!(!ca.check_self_collision(&chain));
    }
    #[test]
    fn sphere_obstacle_signed_distance() {
        let obs = SphereObstacle::new([0.0; 3], 1.0);
        assert!((obs.signed_distance([2.0, 0.0, 0.0]) - 1.0).abs() < 1e-10);
        assert!((obs.signed_distance([0.5, 0.0, 0.0]) + 0.5).abs() < 1e-10);
    }
    #[test]
    fn min_obstacle_distance_positive_when_clear() {
        let obs = SphereObstacle::new([100.0, 0.0, 0.0], 0.1);
        let ca = CollisionAvoidance::new(vec![obs], 0.1, 0.05);
        let chain = make_rr_robot(0.0, 0.0);
        let d = ca.min_obstacle_distance(&chain);
        assert!(
            d > 0.0,
            "distance should be positive when far from obstacle"
        );
    }
    #[test]
    fn impedance_natural_frequency() {
        let p = ImpedanceParams::new(1.0, 2.0, 4.0);
        assert!((p.natural_frequency() - 2.0).abs() < 1e-10);
    }
    #[test]
    fn impedance_damping_ratio() {
        let p = ImpedanceParams::new(1.0, 2.0, 1.0);
        assert!((p.damping_ratio() - 1.0).abs() < 1e-10);
    }
    #[test]
    fn impedance_critically_damped() {
        let p = ImpedanceParams::new(1.0, 2.0, 1.0);
        assert!(p.is_critically_damped());
    }
    #[test]
    fn force_control_acceleration_zero_at_equilibrium() {
        let p = ImpedanceParams::new(1.0, 1.0, 1.0);
        let fc = ForceControl::new(p, [5.0, 0.0, 0.0]);
        let acc = fc.compute_acceleration([0.0; 3], [0.0; 3], [5.0, 0.0, 0.0]);
        assert!(
            acc[0].abs() < 1e-10,
            "acceleration should be zero at equilibrium"
        );
    }
    #[test]
    fn admittance_velocity_nonzero_for_force_error() {
        let p = ImpedanceParams::new(1.0, 10.0, 1.0);
        let fc = ForceControl::new(p, [0.0; 3]);
        let vel = fc.admittance_velocity([5.0, 0.0, 0.0]);
        assert!(
            vel[0].abs() > 0.0,
            "velocity should be nonzero for force error"
        );
    }
    #[test]
    fn pi_force_control_output() {
        let p = ImpedanceParams::new(1.0, 1.0, 1.0);
        let fc = ForceControl::new(p, [10.0, 0.0, 0.0]);
        let dx = fc.pi_force_control([5.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.1, 0.01);
        assert!((dx[0] - 0.51).abs() < 1e-10, "dx={}", dx[0]);
    }
    #[test]
    fn force_closure_two_opposing_contacts() {
        let c1 = ContactPoint::new([1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], 0.5);
        let c2 = ContactPoint::new([-1.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.5);
        let gp = GraspPlanning::new(vec![c1, c2], [0.0; 3]);
        assert!(gp.has_force_closure());
    }
    #[test]
    fn no_force_closure_single_contact() {
        let c = ContactPoint::new([1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], 0.5);
        let gp = GraspPlanning::new(vec![c], [0.0; 3]);
        assert!(!gp.has_force_closure());
    }
    #[test]
    fn grasp_matrix_size() {
        let c1 = ContactPoint::new([1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], 0.5);
        let c2 = ContactPoint::new([-1.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.5);
        let gp = GraspPlanning::new(vec![c1, c2], [0.0; 3]);
        let g = gp.grasp_matrix();
        assert_eq!(g.len(), 36);
    }
    #[test]
    fn contact_point_normal_is_unit() {
        let c = ContactPoint::new([0.0; 3], [3.0, 4.0, 0.0], 0.5);
        let n = vec3_norm(c.normal);
        assert!((n - 1.0).abs() < 1e-10, "normal should be unit: {n}");
    }
    #[test]
    fn wrench_basis_not_empty() {
        let c = ContactPoint::new([1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], 0.5);
        let wb = c.wrench_basis([0.0; 3]);
        assert!(!wb.is_empty(), "wrench basis should not be empty");
    }
    #[test]
    fn cone_half_angle_for_steel() {
        let c = ContactPoint::new([0.0; 3], [0.0, 1.0, 0.0], 0.3);
        let angle_deg = c.cone_half_angle().to_degrees();
        assert!((angle_deg - 16.7).abs() < 0.2, "angle={angle_deg}");
    }
    #[test]
    fn min_norm_forces_returns_correct_size() {
        let c1 = ContactPoint::new([1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], 0.5);
        let c2 = ContactPoint::new([-1.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.5);
        let gp = GraspPlanning::new(vec![c1, c2], [0.0; 3]);
        let f = gp.min_norm_forces([0.0; 6]);
        assert_eq!(f.len(), 6, "should have 3 * 2 contact force components");
    }
    #[test]
    fn quality_q1_positive_for_good_grasp() {
        let c1 = ContactPoint::new([1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], 0.5);
        let c2 = ContactPoint::new([-1.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.5);
        let c3 = ContactPoint::new([0.0, 1.0, 0.0], [0.0, -1.0, 0.0], 0.5);
        let gp = GraspPlanning::new(vec![c1, c2, c3], [0.0; 3]);
        let q = gp.quality_q1();
        assert!(q >= 0.0, "grasp quality should be non-negative: {q}");
    }
    #[test]
    fn mat4_mul_identity() {
        let eye = mat4_eye();
        let result = mat4_mul(&eye, &eye);
        for i in 0..4 {
            for j in 0..4 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((result[i][j] - expected).abs() < 1e-10);
            }
        }
    }
    #[test]
    fn vec3_cross_orthogonal() {
        let a = [1.0, 0.0, 0.0];
        let b = [0.0, 1.0, 0.0];
        let c = vec3_cross(a, b);
        assert!((c[2] - 1.0).abs() < 1e-10);
        assert!(c[0].abs() < 1e-10);
        assert!(c[1].abs() < 1e-10);
    }
    #[test]
    fn solve_linear_identity_system() {
        let a = vec![
            vec![1.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
            vec![0.0, 0.0, 1.0],
        ];
        let b = vec![3.0, 2.0, 1.0];
        let x = solve_linear(&a, &b).unwrap();
        assert!((x[0] - 3.0).abs() < 1e-10);
        assert!((x[1] - 2.0).abs() < 1e-10);
        assert!((x[2] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn det3x3_identity_is_one() {
        let eye = vec![
            vec![1.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
            vec![0.0, 0.0, 1.0],
        ];
        assert!((det3x3(&eye) - 1.0).abs() < 1e-10);
    }
    #[test]
    fn prismatic_joint_translates_along_z() {
        let params = vec![DhParam::new(0.0, 0.0, 0.0, 0.0)];
        let types = vec![JointType::Prismatic];
        let mut chain = KinematicChain::new(params, types, vec![0.0]);
        chain.q[0] = 1.5;
        let pos = chain.end_effector_position();
        assert!((pos[2] - 1.5).abs() < 1e-10, "z={}", pos[2]);
    }
    #[test]
    fn ik_default_solver_can_be_created() {
        let ik = InverseKinematics::default();
        assert_eq!(ik.max_iter, 200);
        assert!(ik.tol > 0.0);
    }
    #[test]
    fn trajectory_empty_returns_none() {
        let tp = TrajectoryPlanning::new(vec![]);
        assert!(tp.linear_interp(0.5).is_none());
        assert!(tp.cubic_spline(0.5).is_none());
    }
    #[test]
    fn forces_within_cones_normal_force() {
        let c = ContactPoint::new([0.0; 3], [0.0, 0.0, 1.0], 0.5);
        let gp = GraspPlanning::new(vec![c], [0.0; 3]);
        let forces = [0.0, 0.0, 1.0];
        assert!(gp.forces_within_cones(&forces));
    }
    #[test]
    fn forces_outside_cones_detected() {
        let c = ContactPoint::new([0.0; 3], [0.0, 0.0, 1.0], 0.1);
        let gp = GraspPlanning::new(vec![c], [0.0; 3]);
        let forces = [1.0, 0.0, 0.1];
        assert!(!gp.forces_within_cones(&forces));
    }
    #[test]
    fn serial_manipulator_fk_matches_kinematic_chain() {
        let sm = SerialManipulator::puma_like();
        let pos = sm.chain.end_effector_position();
        let r = (pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]).sqrt();
        assert!(r > 0.0, "end-effector should not be at origin: r={r}");
    }
    #[test]
    fn serial_manipulator_dof_matches_params() {
        let sm = SerialManipulator::puma_like();
        assert_eq!(sm.chain.dof(), 6);
    }
    #[test]
    fn serial_manipulator_tip_position_finite() {
        let sm = SerialManipulator::puma_like();
        let pos = sm.chain.end_effector_position();
        assert!(
            pos.iter().all(|v| v.is_finite()),
            "tip must be finite: {pos:?}"
        );
    }
    #[test]
    fn jacobian_matrix_struct_rank_2dof_singular() {
        let chain = make_rr_robot(0.0, 0.0);
        let jm = JacobianMatrix::from_chain(&chain);
        assert!(jm.rank(1e-8) <= 1, "rank at singularity should be ≤1");
    }
    #[test]
    fn jacobian_matrix_struct_rank_2dof_generic() {
        let chain = make_rr_robot(PI / 4.0, PI / 3.0);
        let jm = JacobianMatrix::from_chain(&chain);
        assert!(jm.rank(1e-8) >= 2, "rank at generic config should be ≥2");
    }
    #[test]
    fn jacobian_matrix_condition_number_large_at_singularity() {
        let chain = make_rr_robot(0.0, 0.0);
        let jm = JacobianMatrix::from_chain(&chain);
        let cond = jm.condition_number();
        assert!(
            cond > 10.0 || cond == 0.0,
            "condition number at singularity: {cond}"
        );
    }
    #[test]
    fn workspace_analysis_samples_positive_count() {
        let wa = WorkspaceAnalysis::new(100);
        let chain = make_rr_robot(0.0, 0.0);
        let samples = wa.sample_workspace(&chain);
        assert!(!samples.is_empty(), "workspace samples should be non-empty");
    }
    #[test]
    fn workspace_analysis_bounding_box_valid() {
        let wa = WorkspaceAnalysis::new(200);
        let chain = make_rr_robot(0.0, 0.0);
        let (lo, hi) = wa.bounding_box(&chain);
        for i in 0..3 {
            assert!(lo[i] <= hi[i], "lo[{i}]={} > hi[{i}]={}", lo[i], hi[i]);
        }
    }
    #[test]
    fn workspace_analysis_singular_configs_detected() {
        let wa = WorkspaceAnalysis::new(50);
        let chain = make_rr_robot(0.0, 0.0);
        let singulars = wa.detect_singularities(&chain, 1e-6);
        assert!(
            !singulars.is_empty(),
            "should detect at least one singularity"
        );
    }
    #[test]
    fn workspace_radius_2dof_approx() {
        let wa = WorkspaceAnalysis::new(300);
        let chain = make_rr_robot(0.0, 0.0);
        let r = wa.max_reach(&chain);
        assert!((r - 2.0).abs() < 0.05, "max reach should be ≈2.0, got {r}");
    }
    #[test]
    fn grasp_stability_force_closure_opposing_contacts() {
        let c1 = ContactPoint::new([1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], 0.5);
        let c2 = ContactPoint::new([-1.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.5);
        let gs = GraspStability::new(vec![c1, c2], [0.0; 3]);
        assert!(
            gs.force_closure(),
            "opposing contacts should give force closure"
        );
    }
    #[test]
    fn grasp_stability_no_closure_single_contact() {
        let c = ContactPoint::new([1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], 0.5);
        let gs = GraspStability::new(vec![c], [0.0; 3]);
        assert!(
            !gs.force_closure(),
            "single contact cannot give force closure"
        );
    }
    #[test]
    fn grasp_stability_wrench_space_not_empty() {
        let c1 = ContactPoint::new([1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], 0.5);
        let c2 = ContactPoint::new([-1.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.5);
        let gs = GraspStability::new(vec![c1, c2], [0.0; 3]);
        let ws = gs.wrench_space();
        assert!(!ws.is_empty(), "wrench space should not be empty");
    }
    #[test]
    fn grasp_stability_stability_index_nonneg() {
        let c1 = ContactPoint::new([0.0, 1.0, 0.0], [0.0, -1.0, 0.0], 0.4);
        let c2 = ContactPoint::new([0.0, -1.0, 0.0], [0.0, 1.0, 0.0], 0.4);
        let gs = GraspStability::new(vec![c1, c2], [0.0; 3]);
        let idx = gs.stability_index();
        assert!(idx >= 0.0, "stability index should be non-negative: {idx}");
    }
}
