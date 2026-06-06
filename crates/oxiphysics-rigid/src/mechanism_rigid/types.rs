//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

use super::functions::{EYE4, HMatrix, Twist6, Vec3};
use super::functions::{
    det3x3, h_mul, h_transform_point, h_translation, rot_x, rot_z, solve_3x3, trans, v3_cross,
    v3_dot, v3_norm, v3_sub,
};

/// Result of a workspace sampling analysis.
#[derive(Debug, Clone)]
pub struct WorkspaceResult {
    /// Reachable points in 3D.
    pub points: Vec<Vec3>,
    /// Bounding box minimum corner.
    pub bb_min: Vec3,
    /// Bounding box maximum corner.
    pub bb_max: Vec3,
    /// Approximate volume (from bounding box).
    pub approx_volume: f64,
}
/// Full acceleration analysis result (first-order approximation).
#[derive(Debug, Clone)]
pub struct AccelerationAnalysis {
    /// Linear acceleration from J * ddq term.
    pub linear_acceleration: Vec3,
    /// Magnitude of linear acceleration.
    pub magnitude: f64,
}
/// A parallel kinematic mechanism consisting of multiple serial chains
/// sharing a common base and end-effector platform.
#[derive(Debug, Clone)]
pub struct ParallelChain {
    /// The individual serial chains (legs).
    pub legs: Vec<SerialChain>,
    /// Platform attachment points in the platform frame.
    pub platform_points: Vec<Vec3>,
}
impl ParallelChain {
    /// Create a new parallel mechanism.
    pub fn new(legs: Vec<SerialChain>, platform_points: Vec<Vec3>) -> Self {
        Self {
            legs,
            platform_points,
        }
    }
    /// Number of legs.
    pub fn num_legs(&self) -> usize {
        self.legs.len()
    }
    /// Total degrees of freedom across all legs.
    pub fn total_dof(&self) -> usize {
        self.legs.iter().map(|l| l.dof()).sum()
    }
    /// Compute the position error for each leg: difference between the leg
    /// end-effector and the corresponding platform attachment point after
    /// the platform is placed at `platform_pose`.
    pub fn loop_closure_errors(&self, platform_pose: &HMatrix) -> Vec<Vec3> {
        self.legs
            .iter()
            .zip(self.platform_points.iter())
            .map(|(leg, pp)| {
                let leg_ee = leg.end_effector_position();
                let target = h_transform_point(platform_pose, pp);
                v3_sub(&target, &leg_ee)
            })
            .collect()
    }
    /// Sum of squared loop-closure position errors.
    pub fn loop_closure_error_norm(&self, platform_pose: &HMatrix) -> f64 {
        self.loop_closure_errors(platform_pose)
            .iter()
            .map(|e| v3_dot(e, e))
            .sum::<f64>()
            .sqrt()
    }
}
/// Result of four-bar linkage position analysis.
#[derive(Debug, Clone)]
pub struct FourBarResult {
    /// Input crank angle θ₂ (radians).
    pub theta2: f64,
    /// Coupler angle θ₃ (radians).
    pub theta3: f64,
    /// Output rocker angle θ₄ (radians).
    pub theta4: f64,
    /// Position of the coupler joint B (crank-coupler pivot).
    pub point_b: [f64; 2],
    /// Position of the coupler joint C (coupler-rocker pivot).
    pub point_c: [f64; 2],
    /// Assembly mode (+1 or -1).
    pub mode: i32,
}
/// A Geneva mechanism (intermittent-motion device).
///
/// The driving wheel has a single pin that engages slots on the driven
/// (Geneva) wheel, producing intermittent rotation.
#[derive(Debug, Clone)]
pub struct GenevaMechanism {
    /// Number of slots on the Geneva wheel.
    pub num_slots: usize,
    /// Centre distance between driving and Geneva wheel axes.
    pub center_distance: f64,
    /// Radius of the pin circle on the driving wheel.
    pub pin_radius: f64,
}
impl GenevaMechanism {
    /// Create a new Geneva mechanism.
    ///
    /// For a standard Geneva mechanism, `pin_radius = center_distance * sin(π / num_slots)`.
    pub fn new(num_slots: usize, center_distance: f64) -> Self {
        let pin_radius = center_distance * (PI / num_slots as f64).sin();
        Self {
            num_slots,
            center_distance,
            pin_radius,
        }
    }
    /// Create a Geneva mechanism with a custom pin radius.
    pub fn with_pin_radius(num_slots: usize, center_distance: f64, pin_radius: f64) -> Self {
        Self {
            num_slots,
            center_distance,
            pin_radius,
        }
    }
    /// Half-angle of the driving crank during which the Geneva wheel moves (radians).
    pub fn driving_half_angle(&self) -> f64 {
        (self.pin_radius / self.center_distance)
            .clamp(-1.0, 1.0)
            .asin()
    }
    /// Angular advance of the Geneva wheel per cycle (radians).
    pub fn advance_angle(&self) -> f64 {
        2.0 * PI / self.num_slots as f64
    }
    /// Ratio of motion time to total cycle time.
    pub fn motion_time_ratio(&self) -> f64 {
        let half_angle = self.driving_half_angle();
        2.0 * half_angle / (2.0 * PI)
    }
    /// Angular velocity ratio of Geneva wheel to driving wheel at the engagement
    /// point where the driving angle is `phi` (measured from the centre-line).
    ///
    /// The kinematic relationship is:
    /// `ω_g / ω_d = r_p * cos(φ) / sqrt(d² + r_p² - 2 d r_p sin(φ))`
    /// where `d` = centre distance, `r_p` = pin radius.
    pub fn velocity_ratio(&self, phi: f64) -> f64 {
        let d = self.center_distance;
        let rp = self.pin_radius;
        let denom_sq = d * d + rp * rp - 2.0 * d * rp * phi.sin();
        if denom_sq < 1e-30 {
            return 0.0;
        }
        rp * phi.cos() / denom_sq.sqrt()
    }
    /// Maximum angular velocity ratio (occurs when the pin is on the
    /// centre-line, φ = 0).
    pub fn max_velocity_ratio(&self) -> f64 {
        self.velocity_ratio(0.0)
    }
    /// Angular acceleration ratio at driving angle `phi`.
    ///
    /// Numerically approximated.
    pub fn acceleration_ratio(&self, phi: f64, _omega_drive: f64) -> f64 {
        let eps = 1e-7;
        let vr_plus = self.velocity_ratio(phi + eps);
        let vr_minus = self.velocity_ratio(phi - eps);
        (vr_plus - vr_minus) / (2.0 * eps)
    }
    /// Slot depth required on the Geneva wheel.
    pub fn slot_depth(&self) -> f64 {
        let d = self.center_distance;
        let rp = self.pin_radius;
        let r_inner = d - rp;
        let r_outer = (d * d + rp * rp).sqrt();
        r_outer - r_inner
    }
    /// Geneva wheel outer radius (distance from centre to slot bottom).
    pub fn geneva_radius(&self) -> f64 {
        let d = self.center_distance;
        let rp = self.pin_radius;
        (d * d + rp * rp).sqrt()
    }
}
/// Types of mechanical joints.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum JointType {
    /// Revolute joint — rotation about the Z axis.
    Revolute,
    /// Prismatic joint — translation along the Z axis.
    Prismatic,
    /// Fixed / welded joint — no degree of freedom.
    Fixed,
}
/// Result of a singularity analysis at a given configuration.
#[derive(Debug, Clone)]
pub struct SingularityInfo {
    /// Manipulability index w = sqrt(det(J J^T)).
    pub manipulability: f64,
    /// Condition number of the linear Jacobian (ratio of largest to smallest
    /// singular value, approximated).
    pub condition_number: f64,
    /// Whether the configuration is singular (manipulability below threshold).
    pub is_singular: bool,
}
/// Standard Denavit-Hartenberg parameters for a single joint/link.
///
/// Convention: `T_i = Rot_z(theta) * Trans_z(d) * Trans_x(a) * Rot_x(alpha)`.
#[derive(Debug, Clone, Copy)]
pub struct DHParam {
    /// Joint angle θ (radians). For revolute joints this is the variable.
    pub theta: f64,
    /// Link offset d along the Z axis. For prismatic joints this is the variable.
    pub d: f64,
    /// Link length a (distance along X axis).
    pub a: f64,
    /// Link twist α (rotation about X axis, radians).
    pub alpha: f64,
}
impl DHParam {
    /// Create a new DH parameter set.
    pub fn new(theta: f64, d: f64, a: f64, alpha: f64) -> Self {
        Self { theta, d, a, alpha }
    }
    /// Compute the 4×4 homogeneous transform for this DH parameter set.
    pub fn to_transform(&self) -> HMatrix {
        let t1 = rot_z(self.theta);
        let t2 = trans(0.0, 0.0, self.d);
        let t3 = trans(self.a, 0.0, 0.0);
        let t4 = rot_x(self.alpha);
        h_mul(&h_mul(&h_mul(&t1, &t2), &t3), &t4)
    }
}
/// Full velocity analysis result.
#[derive(Debug, Clone)]
pub struct VelocityAnalysis {
    /// End-effector linear velocity \[vx, vy, vz\].
    pub linear_velocity: Vec3,
    /// End-effector angular velocity \[ωx, ωy, ωz\].
    pub angular_velocity: Vec3,
    /// Speed magnitude.
    pub speed: f64,
}
/// A serial (open) kinematic chain.
#[derive(Debug, Clone)]
pub struct SerialChain {
    /// Base transform (world-to-base).
    pub base: HMatrix,
    /// Joints in order from base to end-effector.
    pub joints: Vec<Joint>,
}
impl SerialChain {
    /// Create a new serial chain with an identity base.
    pub fn new(joints: Vec<Joint>) -> Self {
        Self { base: EYE4, joints }
    }
    /// Create a serial chain with a custom base transform.
    pub fn with_base(base: HMatrix, joints: Vec<Joint>) -> Self {
        Self { base, joints }
    }
    /// Number of joints.
    pub fn num_joints(&self) -> usize {
        self.joints.len()
    }
    /// Number of degrees of freedom (excludes fixed joints).
    pub fn dof(&self) -> usize {
        self.joints
            .iter()
            .filter(|j| j.joint_type != JointType::Fixed)
            .count()
    }
    /// Set all joint variables from a slice.
    pub fn set_q(&mut self, q: &[f64]) {
        let mut idx = 0;
        for joint in &mut self.joints {
            if joint.joint_type != JointType::Fixed && idx < q.len() {
                joint.q = q[idx];
                idx += 1;
            }
        }
    }
    /// Get all joint variables as a vector.
    pub fn get_q(&self) -> Vec<f64> {
        self.joints
            .iter()
            .filter(|j| j.joint_type != JointType::Fixed)
            .map(|j| j.q)
            .collect()
    }
    /// Set all joint velocities from a slice.
    pub fn set_dq(&mut self, dq: &[f64]) {
        let mut idx = 0;
        for joint in &mut self.joints {
            if joint.joint_type != JointType::Fixed && idx < dq.len() {
                joint.dq = dq[idx];
                idx += 1;
            }
        }
    }
    /// Get all joint velocities as a vector.
    pub fn get_dq(&self) -> Vec<f64> {
        self.joints
            .iter()
            .filter(|j| j.joint_type != JointType::Fixed)
            .map(|j| j.dq)
            .collect()
    }
    /// Set all joint accelerations from a slice.
    pub fn set_ddq(&mut self, ddq: &[f64]) {
        let mut idx = 0;
        for joint in &mut self.joints {
            if joint.joint_type != JointType::Fixed && idx < ddq.len() {
                joint.ddq = ddq[idx];
                idx += 1;
            }
        }
    }
    /// Compute all intermediate transforms: `T[i]` = base * T_0 * T_1 * ... * T_i.
    pub fn forward_transforms(&self) -> Vec<HMatrix> {
        let mut transforms = Vec::with_capacity(self.joints.len() + 1);
        let mut current = self.base;
        transforms.push(current);
        for joint in &self.joints {
            current = h_mul(&current, &joint.transform());
            transforms.push(current);
        }
        transforms
    }
    /// Forward kinematics — returns the end-effector pose as a 4×4 matrix.
    pub fn forward_kinematics(&self) -> HMatrix {
        let transforms = self.forward_transforms();
        *transforms.last().unwrap_or(&self.base)
    }
    /// End-effector position (translation part of the FK result).
    pub fn end_effector_position(&self) -> Vec3 {
        h_translation(&self.forward_kinematics())
    }
    /// Compute the geometric Jacobian (6×n) at the current configuration.
    ///
    /// Returns a `Vec` of 6-element columns: `[Jω; Jv]` for each DOF.
    /// For revolute joints: `Jω_i = z_i`, `Jv_i = z_i × (p_e - p_i)`.
    /// For prismatic joints: `Jω_i = 0`, `Jv_i = z_i`.
    pub fn jacobian(&self) -> Vec<Twist6> {
        let transforms = self.forward_transforms();
        let p_ee = h_translation(transforms.last().expect("collection should not be empty"));
        let mut cols = Vec::with_capacity(self.dof());
        for (i, joint) in self.joints.iter().enumerate() {
            if joint.joint_type == JointType::Fixed {
                continue;
            }
            let ti = &transforms[i];
            let z_i = [ti[0][2], ti[1][2], ti[2][2]];
            let p_i = h_translation(ti);
            match joint.joint_type {
                JointType::Revolute => {
                    let r = v3_sub(&p_ee, &p_i);
                    let jv = v3_cross(&z_i, &r);
                    cols.push([z_i[0], z_i[1], z_i[2], jv[0], jv[1], jv[2]]);
                }
                JointType::Prismatic => {
                    cols.push([0.0, 0.0, 0.0, z_i[0], z_i[1], z_i[2]]);
                }
                JointType::Fixed => unreachable!(),
            }
        }
        cols
    }
    /// Compute end-effector linear velocity given current joint velocities.
    ///
    /// Returns `[vx, vy, vz]`.
    pub fn end_effector_velocity(&self) -> Vec3 {
        let jac = self.jacobian();
        let dq = self.get_dq();
        let mut v = [0.0; 3];
        for (col, &qi) in jac.iter().zip(dq.iter()) {
            v[0] += col[3] * qi;
            v[1] += col[4] * qi;
            v[2] += col[5] * qi;
        }
        v
    }
    /// Compute end-effector angular velocity given current joint velocities.
    ///
    /// Returns `[ωx, ωy, ωz]`.
    pub fn end_effector_angular_velocity(&self) -> Vec3 {
        let jac = self.jacobian();
        let dq = self.get_dq();
        let mut omega = [0.0; 3];
        for (col, &qi) in jac.iter().zip(dq.iter()) {
            omega[0] += col[0] * qi;
            omega[1] += col[1] * qi;
            omega[2] += col[2] * qi;
        }
        omega
    }
    /// Compute end-effector linear acceleration (approximate, ignoring
    /// Coriolis/centrifugal terms — first-order Jacobian-based).
    ///
    /// `a = J * ddq + dJ * dq`  (this returns only the `J * ddq` term).
    pub fn end_effector_acceleration_linear_term(&self) -> Vec3 {
        let jac = self.jacobian();
        let mut acc = [0.0; 3];
        let mut idx = 0;
        for joint in &self.joints {
            if joint.joint_type == JointType::Fixed {
                continue;
            }
            if idx < jac.len() {
                acc[0] += jac[idx][3] * joint.ddq;
                acc[1] += jac[idx][4] * joint.ddq;
                acc[2] += jac[idx][5] * joint.ddq;
            }
            idx += 1;
        }
        acc
    }
    /// Manipulability index for the linear part of the Jacobian.
    ///
    /// For chains with motion in all 3 dimensions uses `sqrt(det(J J^T))` (3×3).
    /// For planar chains (z-column all zero) uses the 2×2 sub-Jacobian in XY.
    /// A value near zero indicates proximity to a singularity.
    pub fn manipulability(&self) -> f64 {
        let jac = self.jacobian();
        let n = jac.len();
        if n == 0 {
            return 0.0;
        }
        let has_z = jac.iter().any(|col| col[5].abs() > 1e-12);
        if has_z {
            let mut jjt = [[0.0_f64; 3]; 3];
            for col in &jac {
                for r in 0..3 {
                    for c in 0..3 {
                        jjt[r][c] += col[r + 3] * col[c + 3];
                    }
                }
            }
            let det = det3x3(&jjt);
            if det < 0.0 { 0.0 } else { det.sqrt() }
        } else {
            let mut jjt = [[0.0_f64; 2]; 2];
            for col in &jac {
                for r in 0..2 {
                    for c in 0..2 {
                        jjt[r][c] += col[r + 3] * col[c + 3];
                    }
                }
            }
            let det = jjt[0][0] * jjt[1][1] - jjt[0][1] * jjt[1][0];
            if det < 0.0 { 0.0 } else { det.sqrt() }
        }
    }
    /// Check if the current configuration is near a singularity.
    ///
    /// Returns `true` when the manipulability index is below `threshold`.
    pub fn is_singular(&self, threshold: f64) -> bool {
        self.manipulability() < threshold
    }
    /// Simple iterative inverse kinematics using the Jacobian pseudo-inverse
    /// (damped least-squares / Levenberg-Marquardt).
    ///
    /// `target` — desired end-effector position `[x, y, z]`.
    /// `max_iter` — maximum iterations.
    /// `tol` — position error tolerance.
    /// `damping` — damping factor lambda for singularity robustness.
    ///
    /// Automatically detects planar (2D) vs. spatial (3D) chains.
    /// Returns `true` if converged.
    pub fn inverse_kinematics(
        &mut self,
        target: &Vec3,
        max_iter: usize,
        tol: f64,
        damping: f64,
    ) -> bool {
        for _iter in 0..max_iter {
            let p = self.end_effector_position();
            let err = v3_sub(target, &p);
            if v3_norm(&err) < tol {
                return true;
            }
            let jac = self.jacobian();
            let n = jac.len();
            if n == 0 {
                return false;
            }
            let has_z = jac.iter().any(|col| col[5].abs() > 1e-12);
            let dq_vals = if has_z {
                let mut jjt = [[0.0_f64; 3]; 3];
                for col in &jac {
                    for r in 0..3 {
                        for c in 0..3 {
                            jjt[r][c] += col[r + 3] * col[c + 3];
                        }
                    }
                }
                for (i, jjt_row) in jjt.iter_mut().enumerate() {
                    jjt_row[i] += damping * damping;
                }
                let y = solve_3x3(&jjt, &err);
                let mut dq = vec![0.0_f64; n];
                for (j, col) in jac.iter().enumerate() {
                    for k in 0..3 {
                        dq[j] += col[k + 3] * y[k];
                    }
                }
                dq
            } else {
                let err2 = [err[0], err[1]];
                let mut jjt = [[0.0_f64; 2]; 2];
                for col in &jac {
                    for r in 0..2 {
                        for c in 0..2 {
                            jjt[r][c] += col[r + 3] * col[c + 3];
                        }
                    }
                }
                for (i, jjt_row) in jjt.iter_mut().enumerate() {
                    jjt_row[i] += damping * damping;
                }
                let det = jjt[0][0] * jjt[1][1] - jjt[0][1] * jjt[1][0];
                let y = if det.abs() < 1e-30 {
                    [0.0; 2]
                } else {
                    let inv_d = 1.0 / det;
                    [
                        (jjt[1][1] * err2[0] - jjt[0][1] * err2[1]) * inv_d,
                        (-jjt[1][0] * err2[0] + jjt[0][0] * err2[1]) * inv_d,
                    ]
                };
                let mut dq = vec![0.0_f64; n];
                for (j, col) in jac.iter().enumerate() {
                    for k in 0..2 {
                        dq[j] += col[k + 3] * y[k];
                    }
                }
                dq
            };
            let mut idx = 0;
            for joint in &mut self.joints {
                if joint.joint_type == JointType::Fixed {
                    continue;
                }
                if idx < n {
                    joint.q += dq_vals[idx];
                    joint.clamp_q();
                }
                idx += 1;
            }
        }
        let p = self.end_effector_position();
        let err = v3_sub(target, &p);
        v3_norm(&err) < tol
    }
}
/// A single joint in a kinematic chain.
#[derive(Debug, Clone)]
pub struct Joint {
    /// Joint type.
    pub joint_type: JointType,
    /// DH parameters (base values; θ or d may be overridden by joint variable).
    pub dh: DHParam,
    /// Current joint variable value (angle for revolute, displacement for prismatic).
    pub q: f64,
    /// Joint lower limit.
    pub q_min: f64,
    /// Joint upper limit.
    pub q_max: f64,
    /// Joint velocity.
    pub dq: f64,
    /// Joint acceleration.
    pub ddq: f64,
    /// Index of the parent body.
    pub parent_body: usize,
    /// Index of the child body.
    pub child_body: usize,
    /// Joint lower limit (alias).
    pub lower_limit: f64,
    /// Joint upper limit (alias).
    pub upper_limit: f64,
}
impl Joint {
    /// Create a revolute joint with the given DH parameters and limits.
    pub fn revolute(dh: DHParam, q_min: f64, q_max: f64) -> Self {
        Self {
            joint_type: JointType::Revolute,
            dh,
            q: 0.0,
            q_min,
            q_max,
            dq: 0.0,
            ddq: 0.0,
            parent_body: 0,
            child_body: 0,
            lower_limit: q_min,
            upper_limit: q_max,
        }
    }
    /// Create a prismatic joint with the given DH parameters and limits.
    pub fn prismatic(dh: DHParam, q_min: f64, q_max: f64) -> Self {
        Self {
            joint_type: JointType::Prismatic,
            dh,
            q: 0.0,
            q_min,
            q_max,
            dq: 0.0,
            ddq: 0.0,
            parent_body: 0,
            child_body: 0,
            lower_limit: q_min,
            upper_limit: q_max,
        }
    }
    /// Create a fixed joint.
    pub fn fixed(dh: DHParam) -> Self {
        Self {
            joint_type: JointType::Fixed,
            dh,
            q: 0.0,
            q_min: 0.0,
            q_max: 0.0,
            dq: 0.0,
            ddq: 0.0,
            parent_body: 0,
            child_body: 0,
            lower_limit: 0.0,
            upper_limit: 0.0,
        }
    }
    /// Create a revolute joint connecting two bodies with limits.
    pub fn new_revolute(parent_body: usize, child_body: usize, q_min: f64, q_max: f64) -> Self {
        Self {
            joint_type: JointType::Revolute,
            dh: DHParam::new(0.0, 0.0, 0.0, 0.0),
            q: 0.0,
            q_min,
            q_max,
            dq: 0.0,
            ddq: 0.0,
            parent_body,
            child_body,
            lower_limit: q_min,
            upper_limit: q_max,
        }
    }
    /// Create a prismatic joint connecting two bodies with limits.
    pub fn new_prismatic(parent_body: usize, child_body: usize, q_min: f64, q_max: f64) -> Self {
        Self {
            joint_type: JointType::Prismatic,
            dh: DHParam::new(0.0, 0.0, 0.0, 0.0),
            q: 0.0,
            q_min,
            q_max,
            dq: 0.0,
            ddq: 0.0,
            parent_body,
            child_body,
            lower_limit: q_min,
            upper_limit: q_max,
        }
    }
    /// Create a fixed joint connecting two bodies.
    pub fn new_fixed(parent_body: usize, child_body: usize) -> Self {
        Self {
            joint_type: JointType::Fixed,
            dh: DHParam::new(0.0, 0.0, 0.0, 0.0),
            q: 0.0,
            q_min: 0.0,
            q_max: 0.0,
            dq: 0.0,
            ddq: 0.0,
            parent_body,
            child_body,
            lower_limit: 0.0,
            upper_limit: 0.0,
        }
    }
    /// Number of degrees of freedom for this joint.
    pub fn num_dof(&self) -> usize {
        match self.joint_type {
            JointType::Revolute | JointType::Prismatic => 1,
            JointType::Fixed => 0,
        }
    }
    /// Clamp a position value to this joint's limits.
    pub fn clamp_position(&self, value: f64) -> f64 {
        value.clamp(self.q_min, self.q_max)
    }
    /// Clamp the joint variable to its limits.
    pub fn clamp_q(&mut self) {
        if self.q < self.q_min {
            self.q = self.q_min;
        }
        if self.q > self.q_max {
            self.q = self.q_max;
        }
    }
    /// Effective DH parameters (with the joint variable applied).
    pub fn effective_dh(&self) -> DHParam {
        let mut dh = self.dh;
        match self.joint_type {
            JointType::Revolute => dh.theta += self.q,
            JointType::Prismatic => dh.d += self.q,
            JointType::Fixed => {}
        }
        dh
    }
    /// The 4×4 transform for this joint at its current configuration.
    pub fn transform(&self) -> HMatrix {
        self.effective_dh().to_transform()
    }
}
/// A planar four-bar linkage defined by its four link lengths.
///
/// Links: ground (`a`), input crank (`b`), coupler (`c`), output rocker (`d`).
/// The ground link is along the X axis from the origin to `(a, 0)`.
#[derive(Debug, Clone)]
pub struct FourBarLinkage {
    /// Ground link length.
    pub a: f64,
    /// Input crank length.
    pub b: f64,
    /// Coupler length.
    pub c: f64,
    /// Output rocker length.
    pub d: f64,
}
impl FourBarLinkage {
    /// Create a new four-bar linkage.
    pub fn new(a: f64, b: f64, c: f64, d: f64) -> Self {
        Self { a, b, c, d }
    }
    /// Grashof condition: the sum of the shortest and longest links must be
    /// less than or equal to the sum of the other two for at least one link
    /// to make a full rotation.
    pub fn is_grashof(&self) -> bool {
        let mut links = [self.a, self.b, self.c, self.d];
        links.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
        links[0] + links[3] <= links[1] + links[2]
    }
    /// Position analysis for a given input angle θ₂ (radians).
    ///
    /// `mode` selects the assembly mode: +1 (open) or -1 (crossed).
    /// Returns `None` if the configuration is impossible.
    pub fn position_analysis(&self, theta2: f64, mode: i32) -> Option<FourBarResult> {
        let (s2, c2) = theta2.sin_cos();
        let bx = self.b * c2;
        let by = self.b * s2;
        let dx = self.a - bx;
        let dy = -by;
        let diag = (dx * dx + dy * dy).sqrt();
        if diag < 1e-15 {
            return None;
        }
        let phi = dy.atan2(dx);
        let cos_gamma = (diag * diag + self.d * self.d - self.c * self.c) / (2.0 * diag * self.d);
        if cos_gamma.abs() > 1.0 + 1e-10 {
            return None;
        }
        let cos_gamma = cos_gamma.clamp(-1.0, 1.0);
        let gamma = cos_gamma.acos();
        let sign = if mode >= 0 { 1.0 } else { -1.0 };
        let theta4 = phi + sign * gamma;
        let cx = self.a + self.d * theta4.cos();
        let cy = self.d * theta4.sin();
        let theta3 = (cy - by).atan2(cx - bx);
        Some(FourBarResult {
            theta2,
            theta3,
            theta4,
            point_b: [bx, by],
            point_c: [cx, cy],
            mode: if mode >= 0 { 1 } else { -1 },
        })
    }
    /// Velocity analysis for the four-bar linkage.
    ///
    /// Given input angular velocity `omega2`, returns `(omega3, omega4)`.
    pub fn velocity_analysis(&self, theta2: f64, mode: i32, omega2: f64) -> Option<(f64, f64)> {
        let result = self.position_analysis(theta2, mode)?;
        let s2 = theta2.sin();
        let c2 = theta2.cos();
        let s3 = result.theta3.sin();
        let c3 = result.theta3.cos();
        let s4 = result.theta4.sin();
        let c4 = result.theta4.cos();
        let det = self.c * s3 * self.d * c4 - self.c * c3 * self.d * s4;
        if det.abs() < 1e-15 {
            return None;
        }
        let omega3 = (self.b * omega2 * (s2 * self.d * c4 - c2 * self.d * s4)) / det;
        let omega4 = (self.b * omega2 * (s2 * self.c * c3 - c2 * self.c * s3)) / det;
        Some((omega3, omega4))
    }
    /// Compute the transmission angle (angle at the coupler-rocker joint).
    /// A good mechanism has a transmission angle between 40° and 140°.
    pub fn transmission_angle(&self, theta2: f64, mode: i32) -> Option<f64> {
        let result = self.position_analysis(theta2, mode)?;
        let mu = (result.theta4 - result.theta3).abs();
        let mu = mu % PI;
        Some(mu)
    }
    /// Mechanical advantage: ratio of output torque to input torque.
    pub fn mechanical_advantage(&self, theta2: f64, mode: i32) -> Option<f64> {
        let (omega3, omega4) = self.velocity_analysis(theta2, mode, 1.0)?;
        let _omega3 = omega3;
        if omega4.abs() < 1e-15 {
            return None;
        }
        Some(1.0 / omega4)
    }
}
/// A planar crank-slider mechanism.
///
/// The crank rotates about the origin, the slider moves along the X axis.
#[derive(Debug, Clone)]
pub struct CrankSlider {
    /// Crank length.
    pub r: f64,
    /// Connecting rod length.
    pub l: f64,
    /// Offset of the slider from the crank axis (eccentricity), default 0.
    pub offset: f64,
}
impl CrankSlider {
    /// Create a new crank-slider mechanism.
    pub fn new(r: f64, l: f64) -> Self {
        Self { r, l, offset: 0.0 }
    }
    /// Create a crank-slider with eccentricity.
    pub fn with_offset(r: f64, l: f64, offset: f64) -> Self {
        Self { r, l, offset }
    }
    /// Slider position for a given crank angle θ (radians).
    pub fn slider_position(&self, theta: f64) -> f64 {
        let y_crank = self.r * theta.sin() - self.offset;
        let x_crank = self.r * theta.cos();
        let arg = (self.l * self.l - y_crank * y_crank).max(0.0);
        x_crank + arg.sqrt()
    }
    /// Slider velocity for a given crank angle and angular velocity.
    pub fn slider_velocity(&self, theta: f64, omega: f64) -> f64 {
        let (s, c) = theta.sin_cos();
        let y_crank = self.r * s - self.offset;
        let denom = (self.l * self.l - y_crank * y_crank).max(1e-30).sqrt();
        -self.r * s * omega + (self.r * c * omega * y_crank) / denom
    }
    /// Slider acceleration for a given crank angle, angular velocity,
    /// and angular acceleration.
    pub fn slider_acceleration(&self, theta: f64, omega: f64, alpha: f64) -> f64 {
        let eps = 1e-7;
        let v1 = self.slider_velocity(theta + eps, omega);
        let v0 = self.slider_velocity(theta - eps, omega);
        let dvdtheta = (v1 - v0) / (2.0 * eps);
        let v = self.slider_velocity(theta, omega);
        if omega.abs() < 1e-15 {
            self.slider_velocity(theta, 1.0) * alpha
        } else {
            dvdtheta * omega + (v / omega) * alpha
        }
    }
    /// Crank pin position \[x, y\].
    pub fn crank_pin(&self, theta: f64) -> [f64; 2] {
        [self.r * theta.cos(), self.r * theta.sin()]
    }
    /// Connecting rod angle (angle between rod and horizontal).
    pub fn rod_angle(&self, theta: f64) -> f64 {
        let y_crank = self.r * theta.sin() - self.offset;
        let sin_phi = y_crank / self.l;
        sin_phi.clamp(-1.0, 1.0).asin()
    }
    /// Stroke length (difference between TDC and BDC slider positions).
    pub fn stroke(&self) -> f64 {
        let p_tdc = self.slider_position(0.0);
        let p_bdc = self.slider_position(PI);
        (p_tdc - p_bdc).abs()
    }
}
