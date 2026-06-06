//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

use super::types::{
    AccelerationAnalysis, DHParam, FourBarLinkage, Joint, JointType, SerialChain, SingularityInfo,
    VelocityAnalysis, WorkspaceResult,
};

/// 3D vector as a plain array (no nalgebra dependency).
pub type Vec3 = [f64; 3];
/// 4×4 homogeneous transformation matrix (row-major).
pub type HMatrix = [[f64; 4]; 4];
/// 3×3 rotation matrix (row-major).
pub type Rot3 = [[f64; 3]; 3];
/// 6D twist / wrench vector `[ωx, ωy, ωz, vx, vy, vz]`.
pub type Twist6 = [f64; 6];
/// Identity 4×4 homogeneous transform.
pub const EYE4: HMatrix = [
    [1.0, 0.0, 0.0, 0.0],
    [0.0, 1.0, 0.0, 0.0],
    [0.0, 0.0, 1.0, 0.0],
    [0.0, 0.0, 0.0, 1.0],
];
/// Identity 3×3 rotation.
pub const EYE3: Rot3 = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
/// Subtract `b` from `a`.
pub(super) fn v3_sub(a: &Vec3, b: &Vec3) -> Vec3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
/// Dot product of two 3-vectors.
pub(super) fn v3_dot(a: &Vec3, b: &Vec3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
/// Cross product of two 3-vectors.
pub(super) fn v3_cross(a: &Vec3, b: &Vec3) -> Vec3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
/// Euclidean norm.
pub(super) fn v3_norm(v: &Vec3) -> f64 {
    v3_dot(v, v).sqrt()
}
/// Scale a 3-vector.
#[cfg(test)]
fn v3_scale(v: &Vec3, s: f64) -> Vec3 {
    [v[0] * s, v[1] * s, v[2] * s]
}
/// Multiply two 4×4 matrices.
pub(super) fn h_mul(a: &HMatrix, b: &HMatrix) -> HMatrix {
    let mut c = [[0.0_f64; 4]; 4];
    for i in 0..4 {
        for j in 0..4 {
            for k in 0..4 {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    c
}
/// Transform a 3D point by a 4×4 homogeneous matrix.
pub(super) fn h_transform_point(m: &HMatrix, p: &Vec3) -> Vec3 {
    [
        m[0][0] * p[0] + m[0][1] * p[1] + m[0][2] * p[2] + m[0][3],
        m[1][0] * p[0] + m[1][1] * p[1] + m[1][2] * p[2] + m[1][3],
        m[2][0] * p[0] + m[2][1] * p[1] + m[2][2] * p[2] + m[2][3],
    ]
}
/// Extract the translation part of a 4×4 homogeneous matrix.
pub(super) fn h_translation(m: &HMatrix) -> Vec3 {
    [m[0][3], m[1][3], m[2][3]]
}
/// Extract the 3×3 rotation part.
#[cfg(test)]
fn h_rotation(m: &HMatrix) -> Rot3 {
    [
        [m[0][0], m[0][1], m[0][2]],
        [m[1][0], m[1][1], m[1][2]],
        [m[2][0], m[2][1], m[2][2]],
    ]
}
/// Rotate a 3D vector by a 3×3 rotation matrix.
#[cfg(test)]
fn rot_mul_vec(r: &Rot3, v: &Vec3) -> Vec3 {
    [
        r[0][0] * v[0] + r[0][1] * v[1] + r[0][2] * v[2],
        r[1][0] * v[0] + r[1][1] * v[1] + r[1][2] * v[2],
        r[2][0] * v[0] + r[2][1] * v[1] + r[2][2] * v[2],
    ]
}
/// Transpose a 3×3 matrix (rotation inverse for orthonormal matrices).
#[cfg(test)]
fn rot_transpose(r: &Rot3) -> Rot3 {
    [
        [r[0][0], r[1][0], r[2][0]],
        [r[0][1], r[1][1], r[2][1]],
        [r[0][2], r[1][2], r[2][2]],
    ]
}
/// Invert a 4×4 homogeneous transform (assuming orthonormal rotation + translation).
#[cfg(test)]
fn h_inv(m: &HMatrix) -> HMatrix {
    let rt = rot_transpose(&h_rotation(m));
    let t = h_translation(m);
    let inv_t = rot_mul_vec(&rt, &[-t[0], -t[1], -t[2]]);
    [
        [rt[0][0], rt[0][1], rt[0][2], inv_t[0]],
        [rt[1][0], rt[1][1], rt[1][2], inv_t[1]],
        [rt[2][0], rt[2][1], rt[2][2], inv_t[2]],
        [0.0, 0.0, 0.0, 1.0],
    ]
}
/// Build a rotation matrix about the Z axis.
pub(super) fn rot_z(angle: f64) -> HMatrix {
    let (s, c) = angle.sin_cos();
    [
        [c, -s, 0.0, 0.0],
        [s, c, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ]
}
/// Build a rotation matrix about the X axis.
pub(super) fn rot_x(angle: f64) -> HMatrix {
    let (s, c) = angle.sin_cos();
    [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, c, -s, 0.0],
        [0.0, s, c, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ]
}
/// Build a translation matrix.
pub(super) fn trans(dx: f64, dy: f64, dz: f64) -> HMatrix {
    [
        [1.0, 0.0, 0.0, dx],
        [0.0, 1.0, 0.0, dy],
        [0.0, 0.0, 1.0, dz],
        [0.0, 0.0, 0.0, 1.0],
    ]
}
/// Determinant of a 3×3 matrix.
pub(super) fn det3x3(m: &[[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}
/// Solve a 3×3 linear system Ax = b using Cramer's rule.
pub(super) fn solve_3x3(a: &[[f64; 3]; 3], b: &Vec3) -> Vec3 {
    let d = det3x3(a);
    if d.abs() < 1e-30 {
        return [0.0; 3];
    }
    let inv_d = 1.0 / d;
    let x0 = det3x3(&[
        [b[0], a[0][1], a[0][2]],
        [b[1], a[1][1], a[1][2]],
        [b[2], a[2][1], a[2][2]],
    ]) * inv_d;
    let x1 = det3x3(&[
        [a[0][0], b[0], a[0][2]],
        [a[1][0], b[1], a[1][2]],
        [a[2][0], b[2], a[2][2]],
    ]) * inv_d;
    let x2 = det3x3(&[
        [a[0][0], a[0][1], b[0]],
        [a[1][0], a[1][1], b[1]],
        [a[2][0], a[2][1], b[2]],
    ]) * inv_d;
    [x0, x1, x2]
}
/// Sample the reachable workspace of a serial chain by uniformly scanning
/// each joint through `samples_per_joint` steps within its limits.
///
/// Returns a `WorkspaceResult` with all sampled end-effector positions.
pub fn workspace_analysis(chain: &SerialChain, samples_per_joint: usize) -> WorkspaceResult {
    let n = chain.dof();
    if n == 0 {
        let p = chain.end_effector_position();
        return WorkspaceResult {
            points: vec![p],
            bb_min: p,
            bb_max: p,
            approx_volume: 0.0,
        };
    }
    let samples = samples_per_joint.max(2);
    let total = samples.pow(n as u32);
    let mut points = Vec::with_capacity(total.min(500_000));
    let mut bb_min = [f64::INFINITY; 3];
    let mut bb_max = [f64::NEG_INFINITY; 3];
    let limits: Vec<(f64, f64)> = chain
        .joints
        .iter()
        .filter(|j| j.joint_type != JointType::Fixed)
        .map(|j| (j.q_min, j.q_max))
        .collect();
    let mut indices = vec![0usize; n];
    let mut test_chain = chain.clone();
    for _ in 0..total {
        let q_vals: Vec<f64> = indices
            .iter()
            .enumerate()
            .map(|(i, &idx)| {
                let (lo, hi) = limits[i];
                if samples <= 1 {
                    lo
                } else {
                    lo + (hi - lo) * (idx as f64) / ((samples - 1) as f64)
                }
            })
            .collect();
        test_chain.set_q(&q_vals);
        let p = test_chain.end_effector_position();
        for k in 0..3 {
            if p[k] < bb_min[k] {
                bb_min[k] = p[k];
            }
            if p[k] > bb_max[k] {
                bb_max[k] = p[k];
            }
        }
        points.push(p);
        let mut carry = true;
        for d in (0..n).rev() {
            if carry {
                indices[d] += 1;
                if indices[d] >= samples {
                    indices[d] = 0;
                } else {
                    carry = false;
                }
            }
        }
        if carry {
            break;
        }
    }
    let volume = (bb_max[0] - bb_min[0]) * (bb_max[1] - bb_min[1]) * (bb_max[2] - bb_min[2]);
    WorkspaceResult {
        points,
        bb_min,
        bb_max,
        approx_volume: volume.abs(),
    }
}
/// Analyse singularity at the current chain configuration.
pub fn singularity_analysis(chain: &SerialChain, threshold: f64) -> SingularityInfo {
    let manip = chain.manipulability();
    let jac = chain.jacobian();
    let n = jac.len();
    let mut jjt = [[0.0_f64; 3]; 3];
    for col in &jac {
        for r in 0..3 {
            for c in 0..3 {
                jjt[r][c] += col[r + 3] * col[c + 3];
            }
        }
    }
    let trace = jjt[0][0] + jjt[1][1] + jjt[2][2];
    let det = det3x3(&jjt);
    let condition = if det.abs() > 1e-30 && n >= 3 {
        let ratio = trace * trace * trace / det.abs();
        ratio.abs().cbrt()
    } else {
        f64::INFINITY
    };
    SingularityInfo {
        manipulability: manip,
        condition_number: condition,
        is_singular: manip < threshold,
    }
}
/// Perform velocity analysis on a serial chain.
pub fn velocity_analysis(chain: &SerialChain) -> VelocityAnalysis {
    let lv = chain.end_effector_velocity();
    let av = chain.end_effector_angular_velocity();
    let speed = v3_norm(&lv);
    VelocityAnalysis {
        linear_velocity: lv,
        angular_velocity: av,
        speed,
    }
}
/// Perform acceleration analysis on a serial chain (J * ddq term only).
pub fn acceleration_analysis(chain: &SerialChain) -> AccelerationAnalysis {
    let acc = chain.end_effector_acceleration_linear_term();
    let mag = v3_norm(&acc);
    AccelerationAnalysis {
        linear_acceleration: acc,
        magnitude: mag,
    }
}
/// Build a 2-DOF planar RR (revolute-revolute) serial chain.
///
/// `l1`, `l2` — link lengths. Both joints rotate about the Z axis.
pub fn build_planar_rr(l1: f64, l2: f64) -> SerialChain {
    let j1 = Joint::revolute(DHParam::new(0.0, 0.0, l1, 0.0), -PI, PI);
    let j2 = Joint::revolute(DHParam::new(0.0, 0.0, l2, 0.0), -PI, PI);
    SerialChain::new(vec![j1, j2])
}
/// Build a 3-DOF planar RRR chain.
pub fn build_planar_rrr(l1: f64, l2: f64, l3: f64) -> SerialChain {
    let j1 = Joint::revolute(DHParam::new(0.0, 0.0, l1, 0.0), -PI, PI);
    let j2 = Joint::revolute(DHParam::new(0.0, 0.0, l2, 0.0), -PI, PI);
    let j3 = Joint::revolute(DHParam::new(0.0, 0.0, l3, 0.0), -PI, PI);
    SerialChain::new(vec![j1, j2, j3])
}
/// Build a 3-DOF spatial RPR (revolute-prismatic-revolute) chain.
pub fn build_rpr(d_range: (f64, f64)) -> SerialChain {
    let j1 = Joint::revolute(DHParam::new(0.0, 0.0, 0.0, PI / 2.0), -PI, PI);
    let j2 = Joint::prismatic(DHParam::new(0.0, 0.0, 0.0, -PI / 2.0), d_range.0, d_range.1);
    let j3 = Joint::revolute(DHParam::new(0.0, 0.0, 0.0, 0.0), -PI, PI);
    SerialChain::new(vec![j1, j2, j3])
}
/// Build a standard 6-DOF industrial robot arm (PUMA-like DH parameters).
pub fn build_6dof_arm(d1: f64, a2: f64, a3: f64, d4: f64, _d5: f64, d6: f64) -> SerialChain {
    let joints = vec![
        Joint::revolute(DHParam::new(0.0, d1, 0.0, PI / 2.0), -PI, PI),
        Joint::revolute(DHParam::new(0.0, 0.0, a2, 0.0), -PI, PI),
        Joint::revolute(DHParam::new(0.0, 0.0, a3, PI / 2.0), -PI, PI),
        Joint::revolute(DHParam::new(0.0, d4, 0.0, -PI / 2.0), -PI, PI),
        Joint::revolute(DHParam::new(0.0, 0.0, 0.0, PI / 2.0), -PI, PI),
        Joint::revolute(DHParam::new(0.0, d6, 0.0, 0.0), -PI, PI),
    ];
    SerialChain::new(joints)
}
/// Build a SCARA robot (selective compliance assembly robot arm).
pub fn build_scara(l1: f64, l2: f64, d_range: (f64, f64)) -> SerialChain {
    let joints = vec![
        Joint::revolute(DHParam::new(0.0, 0.0, l1, 0.0), -PI, PI),
        Joint::revolute(DHParam::new(0.0, 0.0, l2, PI), -PI, PI),
        Joint::prismatic(DHParam::new(0.0, 0.0, 0.0, 0.0), d_range.0, d_range.1),
        Joint::revolute(DHParam::new(0.0, 0.0, 0.0, 0.0), -PI, PI),
    ];
    SerialChain::new(joints)
}
/// Trace the coupler curve of a four-bar linkage over a full input rotation.
///
/// `coupler_point` is specified in the local frame of the coupler link
/// as `(fraction_along_coupler, perpendicular_offset)`.
///
/// Returns a list of `[x, y]` points.
pub fn trace_coupler_curve(
    linkage: &FourBarLinkage,
    coupler_point: (f64, f64),
    num_samples: usize,
    mode: i32,
) -> Vec<[f64; 2]> {
    let mut curve = Vec::with_capacity(num_samples);
    for i in 0..num_samples {
        let theta2 = 2.0 * PI * (i as f64) / (num_samples as f64);
        if let Some(result) = linkage.position_analysis(theta2, mode) {
            let (frac, perp) = coupler_point;
            let dx = result.point_c[0] - result.point_b[0];
            let dy = result.point_c[1] - result.point_b[1];
            let len = (dx * dx + dy * dy).sqrt();
            if len > 1e-15 {
                let ux = dx / len;
                let uy = dy / len;
                let px = result.point_b[0] + frac * dx + perp * (-uy);
                let py = result.point_b[1] + frac * dy + perp * ux;
                curve.push([px, py]);
            }
        }
    }
    curve
}
/// Compute the degree of freedom of a planar mechanism using the Gruebler
/// equation: `F = 3(n-1) - 2*j1 - j2`, where `n` = number of links,
/// `j1` = full joints (1 DOF removed), `j2` = half joints.
pub fn gruebler_planar(num_links: usize, full_joints: usize, half_joints: usize) -> i32 {
    3 * (num_links as i32 - 1) - 2 * full_joints as i32 - half_joints as i32
}
/// Compute the mobility (DOF) of a spatial mechanism using the Kutzbach
/// criterion: `F = 6(n-1) - Σ(6-fi)`, where `fi` is the DOF of joint i.
pub fn kutzbach_spatial(num_links: usize, joint_dofs: &[usize]) -> i32 {
    let n = num_links as i32;
    let constraint_sum: i32 = joint_dofs.iter().map(|&fi| 6 - fi as i32).sum();
    6 * (n - 1) - constraint_sum
}
/// Convert degrees to radians.
pub fn deg2rad(deg: f64) -> f64 {
    deg * PI / 180.0
}
/// Convert radians to degrees.
pub fn rad2deg(rad: f64) -> f64 {
    rad * 180.0 / PI
}
/// Normalize an angle to the range `[-π, π]`.
pub fn normalize_angle(angle: f64) -> f64 {
    let mut a = angle % (2.0 * PI);
    if a > PI {
        a -= 2.0 * PI;
    }
    if a < -PI {
        a += 2.0 * PI;
    }
    a
}
/// Linear interpolation between two joint configurations.
pub fn lerp_config(q0: &[f64], q1: &[f64], t: f64) -> Vec<f64> {
    q0.iter()
        .zip(q1.iter())
        .map(|(&a, &b)| a + t * (b - a))
        .collect()
}
/// Compute the distance between two 3D points.
pub fn point_distance(a: &Vec3, b: &Vec3) -> f64 {
    v3_norm(&v3_sub(a, b))
}
/// Check if a target point is within the reachable workspace of a serial chain.
///
/// This is a conservative check using the sum-of-link-lengths criterion.
pub fn is_reachable(chain: &SerialChain, target: &Vec3) -> bool {
    let base_pos = h_translation(&chain.base);
    let dist = point_distance(&base_pos, target);
    let max_reach: f64 = chain
        .joints
        .iter()
        .map(|j| j.dh.a.abs() + j.dh.d.abs())
        .sum();
    dist <= max_reach * 1.01
}
/// Compute the maximum reach distance of a serial chain from its base.
pub fn max_reach(chain: &SerialChain) -> f64 {
    chain
        .joints
        .iter()
        .map(|j| j.dh.a.abs() + j.dh.d.abs())
        .sum()
}
#[cfg(test)]
mod tests {
    use super::super::types::*;
    use super::*;
    const TOL: f64 = 1e-8;
    #[test]
    fn test_dh_identity() {
        let dh = DHParam::new(0.0, 0.0, 0.0, 0.0);
        let t = dh.to_transform();
        for (i, row) in t.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (val - expected).abs() < TOL,
                    "t[{i}][{j}] = {}, expected {expected}",
                    val
                );
            }
        }
    }
    #[test]
    fn test_dh_pure_rotation_z() {
        let dh = DHParam::new(PI / 2.0, 0.0, 0.0, 0.0);
        let t = dh.to_transform();
        assert!((t[0][0] - 0.0).abs() < TOL);
        assert!((t[0][1] - (-1.0)).abs() < TOL);
        assert!((t[1][0] - 1.0).abs() < TOL);
        assert!((t[1][1] - 0.0).abs() < TOL);
    }
    #[test]
    fn test_dh_translation_a() {
        let dh = DHParam::new(0.0, 0.0, 5.0, 0.0);
        let t = dh.to_transform();
        let pos = h_translation(&t);
        assert!((pos[0] - 5.0).abs() < TOL, "x = {}", pos[0]);
        assert!(pos[1].abs() < TOL);
        assert!(pos[2].abs() < TOL);
    }
    #[test]
    fn test_dh_translation_d() {
        let dh = DHParam::new(0.0, 3.0, 0.0, 0.0);
        let t = dh.to_transform();
        let pos = h_translation(&t);
        assert!(pos[0].abs() < TOL);
        assert!(pos[1].abs() < TOL);
        assert!((pos[2] - 3.0).abs() < TOL, "z = {}", pos[2]);
    }
    #[test]
    fn test_joint_clamp() {
        let mut j = Joint::revolute(DHParam::new(0.0, 0.0, 1.0, 0.0), -1.0, 1.0);
        j.q = 5.0;
        j.clamp_q();
        assert!((j.q - 1.0).abs() < TOL);
        j.q = -5.0;
        j.clamp_q();
        assert!((j.q - (-1.0)).abs() < TOL);
    }
    #[test]
    fn test_joint_effective_dh_revolute() {
        let j = Joint {
            joint_type: JointType::Revolute,
            dh: DHParam::new(0.5, 0.0, 1.0, 0.0),
            q: 0.3,
            q_min: -PI,
            q_max: PI,
            dq: 0.0,
            ddq: 0.0,
            parent_body: 0,
            child_body: 1,
            lower_limit: -PI,
            upper_limit: PI,
        };
        let edh = j.effective_dh();
        assert!((edh.theta - 0.8).abs() < TOL);
    }
    #[test]
    fn test_joint_effective_dh_prismatic() {
        let j = Joint {
            joint_type: JointType::Prismatic,
            dh: DHParam::new(0.0, 1.0, 0.0, 0.0),
            q: 2.0,
            q_min: 0.0,
            q_max: 5.0,
            dq: 0.0,
            ddq: 0.0,
            parent_body: 0,
            child_body: 1,
            lower_limit: 0.0,
            upper_limit: 5.0,
        };
        let edh = j.effective_dh();
        assert!((edh.d - 3.0).abs() < TOL);
    }
    #[test]
    fn test_planar_rr_fk_zero_config() {
        let chain = build_planar_rr(1.0, 1.0);
        let pos = chain.end_effector_position();
        assert!((pos[0] - 2.0).abs() < TOL, "x = {}", pos[0]);
        assert!(pos[1].abs() < TOL);
        assert!(pos[2].abs() < TOL);
    }
    #[test]
    fn test_planar_rr_fk_90_degrees() {
        let mut chain = build_planar_rr(1.0, 1.0);
        chain.set_q(&[PI / 2.0, 0.0]);
        let pos = chain.end_effector_position();
        assert!(pos[0].abs() < TOL, "x = {}", pos[0]);
        assert!((pos[1] - 2.0).abs() < TOL, "y = {}", pos[1]);
    }
    #[test]
    fn test_planar_rr_fk_folded() {
        let mut chain = build_planar_rr(1.0, 1.0);
        chain.set_q(&[0.0, PI]);
        let pos = chain.end_effector_position();
        assert!(pos[0].abs() < TOL, "x = {}", pos[0]);
        assert!(pos[1].abs() < TOL, "y = {}", pos[1]);
    }
    #[test]
    fn test_serial_chain_dof() {
        let chain = build_planar_rrr(1.0, 1.0, 1.0);
        assert_eq!(chain.dof(), 3);
    }
    #[test]
    fn test_serial_chain_set_get_q() {
        let mut chain = build_planar_rr(1.0, 1.0);
        chain.set_q(&[0.5, 1.0]);
        let q = chain.get_q();
        assert_eq!(q.len(), 2);
        assert!((q[0] - 0.5).abs() < TOL);
        assert!((q[1] - 1.0).abs() < TOL);
    }
    #[test]
    fn test_jacobian_column_count() {
        let chain = build_planar_rrr(1.0, 1.0, 1.0);
        let jac = chain.jacobian();
        assert_eq!(jac.len(), 3, "3 DOF → 3 Jacobian columns");
    }
    #[test]
    fn test_jacobian_velocity_consistency() {
        let mut chain = build_planar_rr(1.0, 1.0);
        chain.set_q(&[0.3, 0.5]);
        chain.set_dq(&[1.0, 0.0]);
        let v = chain.end_effector_velocity();
        let eps = 1e-7;
        let q0 = chain.get_q();
        let p0 = chain.end_effector_position();
        chain.set_q(&[q0[0] + eps, q0[1]]);
        let p1 = chain.end_effector_position();
        chain.set_q(&q0);
        let v_num = v3_scale(&v3_sub(&p1, &p0), 1.0 / eps);
        for k in 0..3 {
            assert!(
                (v[k] - v_num[k]).abs() < 1e-4,
                "v[{k}] = {}, v_num = {}",
                v[k],
                v_num[k]
            );
        }
    }
    #[test]
    fn test_manipulability_nonzero_at_generic_config() {
        let mut chain = build_planar_rr(1.0, 1.0);
        chain.set_q(&[0.5, 0.8]);
        let m = chain.manipulability();
        assert!(m > 0.0, "manipulability = {m}");
    }
    #[test]
    fn test_singular_at_full_extension() {
        let mut chain = build_planar_rr(1.0, 1.0);
        chain.set_q(&[0.0, 0.0]);
        let info = singularity_analysis(&chain, 0.01);
        assert!(
            info.is_singular,
            "fully extended RR should be singular, manip = {}",
            info.manipulability
        );
    }
    #[test]
    fn test_ik_reachable_target() {
        let mut chain = build_planar_rr(1.0, 1.0);
        chain.set_q(&[0.5, 0.5]);
        let target = [1.0, 1.0, 0.0];
        let ok = chain.inverse_kinematics(&target, 500, 1e-3, 0.001);
        assert!(ok, "IK should converge for reachable target");
        let p = chain.end_effector_position();
        let err = v3_norm(&v3_sub(&p, &target));
        assert!(err < 1e-3, "IK error = {err}");
    }
    #[test]
    fn test_ik_unreachable_target() {
        let mut chain = build_planar_rr(1.0, 1.0);
        let target = [10.0, 10.0, 0.0];
        let ok = chain.inverse_kinematics(&target, 50, 1e-4, 0.01);
        assert!(!ok, "IK should not converge for unreachable target");
    }
    #[test]
    fn test_workspace_rr_bounding_box() {
        let chain = build_planar_rr(1.0, 1.0);
        let ws = workspace_analysis(&chain, 10);
        assert!(ws.bb_max[0] > 1.5, "bb_max x = {}", ws.bb_max[0]);
        assert!(ws.bb_min[0] < -1.5, "bb_min x = {}", ws.bb_min[0]);
        let area = (ws.bb_max[0] - ws.bb_min[0]) * (ws.bb_max[1] - ws.bb_min[1]);
        assert!(area > 0.0, "2D workspace area = {area}");
    }
    #[test]
    fn test_workspace_has_points() {
        let chain = build_planar_rr(1.0, 1.0);
        let ws = workspace_analysis(&chain, 5);
        assert_eq!(ws.points.len(), 25);
    }
    #[test]
    fn test_fourbar_grashof() {
        let fb = FourBarLinkage::new(4.0, 2.0, 3.0, 3.5);
        assert!(fb.is_grashof(), "Should be Grashof: 2+4 <= 3+3.5");
    }
    #[test]
    fn test_fourbar_not_grashof() {
        let fb = FourBarLinkage::new(4.0, 1.0, 1.0, 1.0);
        assert!(!fb.is_grashof(), "Should not be Grashof");
    }
    #[test]
    fn test_fourbar_position_analysis() {
        let fb = FourBarLinkage::new(4.0, 2.0, 3.0, 3.5);
        let result = fb.position_analysis(0.5, 1);
        assert!(result.is_some(), "Should find a valid configuration");
        let r = result.unwrap();
        let rb = (r.point_b[0] * r.point_b[0] + r.point_b[1] * r.point_b[1]).sqrt();
        assert!((rb - 2.0).abs() < TOL, "B radius = {rb}");
    }
    #[test]
    fn test_fourbar_velocity_analysis() {
        let fb = FourBarLinkage::new(4.0, 2.0, 3.0, 3.5);
        let result = fb.velocity_analysis(0.5, 1, 1.0);
        assert!(result.is_some());
        let (omega3, omega4) = result.unwrap();
        assert!(omega3.is_finite(), "omega3 = {omega3}");
        assert!(omega4.is_finite(), "omega4 = {omega4}");
    }
    #[test]
    fn test_fourbar_transmission_angle() {
        let fb = FourBarLinkage::new(4.0, 2.0, 3.0, 3.5);
        let mu = fb.transmission_angle(0.5, 1);
        assert!(mu.is_some());
        let mu = mu.unwrap();
        assert!((0.0..=PI).contains(&mu), "mu = {mu}");
    }
    #[test]
    fn test_crank_slider_tdc() {
        let cs = CrankSlider::new(1.0, 3.0);
        let p_tdc = cs.slider_position(0.0);
        assert!((p_tdc - 4.0).abs() < TOL, "TDC = {p_tdc}");
    }
    #[test]
    fn test_crank_slider_bdc() {
        let cs = CrankSlider::new(1.0, 3.0);
        let p_bdc = cs.slider_position(PI);
        assert!((p_bdc - 2.0).abs() < TOL, "BDC = {p_bdc}");
    }
    #[test]
    fn test_crank_slider_stroke() {
        let cs = CrankSlider::new(1.0, 3.0);
        let stroke = cs.stroke();
        assert!((stroke - 2.0).abs() < TOL, "stroke = {stroke}");
    }
    #[test]
    fn test_crank_slider_rod_angle_at_tdc() {
        let cs = CrankSlider::new(1.0, 3.0);
        let phi = cs.rod_angle(0.0);
        assert!(phi.abs() < TOL, "rod angle at TDC = {phi}");
    }
    #[test]
    fn test_geneva_advance_angle() {
        let gm = GenevaMechanism::new(4, 10.0);
        let advance = gm.advance_angle();
        assert!((advance - PI / 2.0).abs() < TOL, "advance = {advance}");
    }
    #[test]
    fn test_geneva_six_slot_advance() {
        let gm = GenevaMechanism::new(6, 10.0);
        let advance = gm.advance_angle();
        assert!(
            (advance - PI / 3.0).abs() < TOL,
            "6-slot advance = {advance}"
        );
    }
    #[test]
    fn test_geneva_max_velocity_ratio() {
        let gm = GenevaMechanism::new(4, 10.0);
        let vr = gm.max_velocity_ratio();
        let expected = gm.pin_radius
            / (gm.center_distance * gm.center_distance + gm.pin_radius * gm.pin_radius).sqrt();
        assert!(
            (vr - expected).abs() < 1e-6,
            "max vr = {vr}, expected {expected}"
        );
    }
    #[test]
    fn test_geneva_slot_depth_positive() {
        let gm = GenevaMechanism::new(6, 10.0);
        let sd = gm.slot_depth();
        assert!(sd > 0.0, "slot depth = {sd}");
    }
    #[test]
    fn test_gruebler_four_bar() {
        let f = gruebler_planar(4, 4, 0);
        assert_eq!(f, 1);
    }
    #[test]
    fn test_kutzbach_spatial_6dof() {
        let f = kutzbach_spatial(7, &[1, 1, 1, 1, 1, 1]);
        assert_eq!(f, 6);
    }
    #[test]
    fn test_normalize_angle() {
        assert!((normalize_angle(3.0 * PI) - PI).abs() < TOL);
        assert!((normalize_angle(-3.0 * PI) - (-PI)).abs() < TOL);
        assert!((normalize_angle(0.5) - 0.5).abs() < TOL);
    }
    #[test]
    fn test_deg2rad_rad2deg_roundtrip() {
        let deg = 45.0;
        let rad = deg2rad(deg);
        let back = rad2deg(rad);
        assert!((back - deg).abs() < TOL);
    }
    #[test]
    fn test_lerp_config() {
        let q0 = vec![0.0, 0.0];
        let q1 = vec![1.0, 2.0];
        let q_mid = lerp_config(&q0, &q1, 0.5);
        assert!((q_mid[0] - 0.5).abs() < TOL);
        assert!((q_mid[1] - 1.0).abs() < TOL);
    }
    #[test]
    fn test_is_reachable() {
        let chain = build_planar_rr(1.0, 1.0);
        assert!(is_reachable(&chain, &[1.5, 0.0, 0.0]));
        assert!(!is_reachable(&chain, &[10.0, 0.0, 0.0]));
    }
    #[test]
    fn test_max_reach() {
        let chain = build_planar_rr(1.0, 1.5);
        assert!((max_reach(&chain) - 2.5).abs() < TOL);
    }
    #[test]
    fn test_parallel_chain_loop_closure() {
        let leg1 = build_planar_rr(1.0, 1.0);
        let leg2 = build_planar_rr(1.0, 1.0);
        let pc = ParallelChain::new(vec![leg1, leg2], vec![[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]]);
        let err = pc.loop_closure_error_norm(&EYE4);
        assert!(err > 0.0);
    }
    #[test]
    fn test_h_inv_roundtrip() {
        let dh = DHParam::new(0.5, 1.0, 2.0, 0.3);
        let t = dh.to_transform();
        let ti = h_inv(&t);
        let product = h_mul(&t, &ti);
        for (i, row) in product.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (val - expected).abs() < 1e-10,
                    "product[{i}][{j}] = {}",
                    val
                );
            }
        }
    }
    #[test]
    fn test_velocity_analysis_zero_config() {
        let mut chain = build_planar_rr(1.0, 1.0);
        chain.set_dq(&[1.0, 0.0]);
        let va = velocity_analysis(&chain);
        assert!(
            (va.linear_velocity[1] - 2.0).abs() < 1e-6,
            "vy = {}",
            va.linear_velocity[1]
        );
        assert!(va.speed > 0.0);
    }
    #[test]
    fn test_coupler_curve_has_points() {
        let fb = FourBarLinkage::new(4.0, 2.0, 3.0, 3.5);
        let curve = trace_coupler_curve(&fb, (0.5, 0.0), 100, 1);
        assert!(!curve.is_empty(), "coupler curve should have points");
    }
    #[test]
    fn test_build_6dof_arm() {
        let arm = build_6dof_arm(0.5, 1.0, 0.5, 0.5, 0.0, 0.1);
        assert_eq!(arm.dof(), 6);
        let pos = arm.end_effector_position();
        for (k, &p) in pos.iter().enumerate() {
            assert!(p.is_finite(), "pos[{k}] = {}", p);
        }
    }
    #[test]
    fn test_build_scara() {
        let scara = build_scara(1.0, 1.0, (0.0, 0.5));
        assert_eq!(scara.dof(), 4);
    }
}
