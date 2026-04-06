#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Mechanism analysis for rigid body systems.
//!
//! Provides kinematic chain modelling (serial and parallel), Denavit-Hartenberg
//! parameterisation, forward and inverse kinematics, Jacobian computation,
//! workspace analysis, singularity detection, velocity/acceleration analysis,
//! four-bar linkage, crank-slider mechanism, and Geneva mechanism analysis.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ══════════════════════════════════════════════════════════════════════════════
// § 1  PRIMITIVE TYPES
// ══════════════════════════════════════════════════════════════════════════════

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

// ══════════════════════════════════════════════════════════════════════════════
// § 2  VECTOR / MATRIX HELPERS
// ══════════════════════════════════════════════════════════════════════════════

/// Add two 3-vectors.
fn v3_add(a: &Vec3, b: &Vec3) -> Vec3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Subtract `b` from `a`.
fn v3_sub(a: &Vec3, b: &Vec3) -> Vec3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Dot product of two 3-vectors.
fn v3_dot(a: &Vec3, b: &Vec3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Cross product of two 3-vectors.
fn v3_cross(a: &Vec3, b: &Vec3) -> Vec3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Euclidean norm.
fn v3_norm(v: &Vec3) -> f64 {
    v3_dot(v, v).sqrt()
}

/// Scale a 3-vector.
fn v3_scale(v: &Vec3, s: f64) -> Vec3 {
    [v[0] * s, v[1] * s, v[2] * s]
}

/// Normalize a 3-vector. Returns zero vector if length is near zero.
fn v3_normalize(v: &Vec3) -> Vec3 {
    let n = v3_norm(v);
    if n < 1e-15 {
        [0.0; 3]
    } else {
        v3_scale(v, 1.0 / n)
    }
}

/// Multiply two 4×4 matrices.
fn h_mul(a: &HMatrix, b: &HMatrix) -> HMatrix {
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
fn h_transform_point(m: &HMatrix, p: &Vec3) -> Vec3 {
    [
        m[0][0] * p[0] + m[0][1] * p[1] + m[0][2] * p[2] + m[0][3],
        m[1][0] * p[0] + m[1][1] * p[1] + m[1][2] * p[2] + m[1][3],
        m[2][0] * p[0] + m[2][1] * p[1] + m[2][2] * p[2] + m[2][3],
    ]
}

/// Extract the translation part of a 4×4 homogeneous matrix.
fn h_translation(m: &HMatrix) -> Vec3 {
    [m[0][3], m[1][3], m[2][3]]
}

/// Extract the 3×3 rotation part.
fn h_rotation(m: &HMatrix) -> Rot3 {
    [
        [m[0][0], m[0][1], m[0][2]],
        [m[1][0], m[1][1], m[1][2]],
        [m[2][0], m[2][1], m[2][2]],
    ]
}

/// Rotate a 3D vector by a 3×3 rotation matrix.
fn rot_mul_vec(r: &Rot3, v: &Vec3) -> Vec3 {
    [
        r[0][0] * v[0] + r[0][1] * v[1] + r[0][2] * v[2],
        r[1][0] * v[0] + r[1][1] * v[1] + r[1][2] * v[2],
        r[2][0] * v[0] + r[2][1] * v[1] + r[2][2] * v[2],
    ]
}

/// Transpose a 3×3 matrix (rotation inverse for orthonormal matrices).
fn rot_transpose(r: &Rot3) -> Rot3 {
    [
        [r[0][0], r[1][0], r[2][0]],
        [r[0][1], r[1][1], r[2][1]],
        [r[0][2], r[1][2], r[2][2]],
    ]
}

/// Invert a 4×4 homogeneous transform (assuming orthonormal rotation + translation).
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
fn rot_z(angle: f64) -> HMatrix {
    let (s, c) = angle.sin_cos();
    [
        [c, -s, 0.0, 0.0],
        [s, c, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ]
}

/// Build a rotation matrix about the X axis.
fn rot_x(angle: f64) -> HMatrix {
    let (s, c) = angle.sin_cos();
    [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, c, -s, 0.0],
        [0.0, s, c, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ]
}

/// Build a translation matrix.
fn trans(dx: f64, dy: f64, dz: f64) -> HMatrix {
    [
        [1.0, 0.0, 0.0, dx],
        [0.0, 1.0, 0.0, dy],
        [0.0, 0.0, 1.0, dz],
        [0.0, 0.0, 0.0, 1.0],
    ]
}

// ══════════════════════════════════════════════════════════════════════════════
// § 3  DENAVIT-HARTENBERG PARAMETERS
// ══════════════════════════════════════════════════════════════════════════════

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

// ══════════════════════════════════════════════════════════════════════════════
// § 4  JOINT TYPES
// ══════════════════════════════════════════════════════════════════════════════

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

// ══════════════════════════════════════════════════════════════════════════════
// § 5  SERIAL KINEMATIC CHAIN
// ══════════════════════════════════════════════════════════════════════════════

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
            let ti = &transforms[i]; // transform up to joint i (before applying)
            let z_i = [ti[0][2], ti[1][2], ti[2][2]]; // Z-axis of frame i
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
        // Check if any column has a non-negligible z linear component
        let has_z = jac.iter().any(|col| col[5].abs() > 1e-12);

        if has_z {
            // Full 3D: 3×3 JJT
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
            // Planar: use 2×2 JJT (XY only, rows 3 and 4)
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

            // Detect if chain is planar (no z linear motion)
            let has_z = jac.iter().any(|col| col[5].abs() > 1e-12);

            let dq_vals = if has_z {
                // Full 3D damped least-squares
                let mut jjt = [[0.0_f64; 3]; 3];
                for col in &jac {
                    for r in 0..3 {
                        for c in 0..3 {
                            jjt[r][c] += col[r + 3] * col[c + 3];
                        }
                    }
                }
                for i in 0..3 {
                    jjt[i][i] += damping * damping;
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
                // 2D (XY) damped least-squares
                let err2 = [err[0], err[1]];
                let mut jjt = [[0.0_f64; 2]; 2];
                for col in &jac {
                    for r in 0..2 {
                        for c in 0..2 {
                            jjt[r][c] += col[r + 3] * col[c + 3];
                        }
                    }
                }
                for i in 0..2 {
                    jjt[i][i] += damping * damping;
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

            // Update joint variables
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

/// Determinant of a 3×3 matrix.
fn det3x3(m: &[[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

/// Solve a 3×3 linear system Ax = b using Cramer's rule.
fn solve_3x3(a: &[[f64; 3]; 3], b: &Vec3) -> Vec3 {
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

// ══════════════════════════════════════════════════════════════════════════════
// § 6  PARALLEL KINEMATIC CHAIN
// ══════════════════════════════════════════════════════════════════════════════

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

// ══════════════════════════════════════════════════════════════════════════════
// § 7  WORKSPACE ANALYSIS
// ══════════════════════════════════════════════════════════════════════════════

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

    // Collect joint limits
    let limits: Vec<(f64, f64)> = chain
        .joints
        .iter()
        .filter(|j| j.joint_type != JointType::Fixed)
        .map(|j| (j.q_min, j.q_max))
        .collect();

    let mut indices = vec![0usize; n];
    let mut test_chain = chain.clone();

    for _ in 0..total {
        // Set joint values from indices
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

        // Increment multi-index
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

// ══════════════════════════════════════════════════════════════════════════════
// § 8  SINGULARITY DETECTION
// ══════════════════════════════════════════════════════════════════════════════

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

/// Analyse singularity at the current chain configuration.
pub fn singularity_analysis(chain: &SerialChain, threshold: f64) -> SingularityInfo {
    let manip = chain.manipulability();
    let jac = chain.jacobian();
    let n = jac.len();

    // Approximate condition number from J_lin * J_lin^T eigenvalues
    let mut jjt = [[0.0_f64; 3]; 3];
    for col in &jac {
        for r in 0..3 {
            for c in 0..3 {
                jjt[r][c] += col[r + 3] * col[c + 3];
            }
        }
    }

    // Use trace and determinant to get eigenvalue bounds
    let trace = jjt[0][0] + jjt[1][1] + jjt[2][2];
    let det = det3x3(&jjt);

    let condition = if det.abs() > 1e-30 && n >= 3 {
        // Rough estimate: max_eig / min_eig ≈ trace / (det / (trace^2/4))
        // Use a simpler heuristic: trace^3 / det provides a rough condition number.
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

// ══════════════════════════════════════════════════════════════════════════════
// § 9  VELOCITY & ACCELERATION ANALYSIS
// ══════════════════════════════════════════════════════════════════════════════

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

/// Full acceleration analysis result (first-order approximation).
#[derive(Debug, Clone)]
pub struct AccelerationAnalysis {
    /// Linear acceleration from J * ddq term.
    pub linear_acceleration: Vec3,
    /// Magnitude of linear acceleration.
    pub magnitude: f64,
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

// ══════════════════════════════════════════════════════════════════════════════
// § 10  FOUR-BAR LINKAGE
// ══════════════════════════════════════════════════════════════════════════════

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

        // Use the law of cosines approach
        let dx = self.a - bx;
        let dy = -by;
        let diag = (dx * dx + dy * dy).sqrt();

        if diag < 1e-15 {
            return None;
        }

        // Angle of the diagonal from B to D (fixed pivot at (a, 0))
        let phi = dy.atan2(dx);

        // Angle at D using law of cosines: c^2 = diag^2 + d^2 - 2*diag*d*cos(gamma)
        let cos_gamma = (diag * diag + self.d * self.d - self.c * self.c) / (2.0 * diag * self.d);

        if cos_gamma.abs() > 1.0 + 1e-10 {
            return None; // impossible configuration
        }
        let cos_gamma = cos_gamma.clamp(-1.0, 1.0);
        let gamma = cos_gamma.acos();

        let sign = if mode >= 0 { 1.0 } else { -1.0 };
        let theta4 = phi + sign * gamma;

        // Point C
        let cx = self.a + self.d * theta4.cos();
        let cy = self.d * theta4.sin();

        // Coupler angle
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

        // Velocity equations from the loop closure derivative:
        // -b*s2*ω2 - c*s3*ω3 + d*s4*ω4 = 0
        //  b*c2*ω2 + c*c3*ω3 - d*c4*ω4 = 0
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
        // Normalize to [0, π]
        let mu = mu % PI;
        Some(mu)
    }

    /// Mechanical advantage: ratio of output torque to input torque.
    pub fn mechanical_advantage(&self, theta2: f64, mode: i32) -> Option<f64> {
        let (omega3, omega4) = self.velocity_analysis(theta2, mode, 1.0)?;
        let _omega3 = omega3; // coupler angular velocity (unused for MA calculation)
        if omega4.abs() < 1e-15 {
            return None;
        }
        // MA = ω2 / ω4 (with ω2=1.0)
        Some(1.0 / omega4)
    }
}

// ══════════════════════════════════════════════════════════════════════════════
// § 11  CRANK-SLIDER MECHANISM
// ══════════════════════════════════════════════════════════════════════════════

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
        // Numerical differentiation of velocity for simplicity
        let eps = 1e-7;
        let v1 = self.slider_velocity(theta + eps, omega);
        let v0 = self.slider_velocity(theta - eps, omega);
        let dvdtheta = (v1 - v0) / (2.0 * eps);
        // a = dv/dt = (dv/dθ) * ω + (∂v/∂ω) * α
        // The (∂v/∂ω) part: v is linear in ω, so ∂v/∂ω = v/ω (when ω ≠ 0)
        let v = self.slider_velocity(theta, omega);
        if omega.abs() < 1e-15 {
            // When ω=0, only α contribution matters
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

// ══════════════════════════════════════════════════════════════════════════════
// § 12  GENEVA MECHANISM
// ══════════════════════════════════════════════════════════════════════════════

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
        // At full engagement the pin is at distance d - rp from the Geneva center,
        // and at max extension it's at sqrt(d^2 + rp^2).
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

// ══════════════════════════════════════════════════════════════════════════════
// § 13  MECHANISM BUILDER HELPERS
// ══════════════════════════════════════════════════════════════════════════════

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

// ══════════════════════════════════════════════════════════════════════════════
// § 14  COUPLER CURVE TRACING
// ══════════════════════════════════════════════════════════════════════════════

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

// ══════════════════════════════════════════════════════════════════════════════
// § 15  ADJACENCY MATRIX / GRAPH REPRESENTATION
// ══════════════════════════════════════════════════════════════════════════════

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

// ══════════════════════════════════════════════════════════════════════════════
// § 16  MISCELLANEOUS UTILITIES
// ══════════════════════════════════════════════════════════════════════════════

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
    dist <= max_reach * 1.01 // small tolerance
}

/// Compute the maximum reach distance of a serial chain from its base.
pub fn max_reach(chain: &SerialChain) -> f64 {
    chain
        .joints
        .iter()
        .map(|j| j.dh.a.abs() + j.dh.d.abs())
        .sum()
}

// ══════════════════════════════════════════════════════════════════════════════
// § 17  TESTS
// ══════════════════════════════════════════════════════════════════════════════

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1e-8;

    // ── DH transform tests ──────────────────────────────────────────────

    #[test]
    fn test_dh_identity() {
        let dh = DHParam::new(0.0, 0.0, 0.0, 0.0);
        let t = dh.to_transform();
        for i in 0..4 {
            for j in 0..4 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (t[i][j] - expected).abs() < TOL,
                    "t[{i}][{j}] = {}, expected {expected}",
                    t[i][j]
                );
            }
        }
    }

    #[test]
    fn test_dh_pure_rotation_z() {
        let dh = DHParam::new(PI / 2.0, 0.0, 0.0, 0.0);
        let t = dh.to_transform();
        // Should rotate 90° about Z: x -> y, y -> -x
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

    // ── Joint tests ─────────────────────────────────────────────────────

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

    // ── Serial chain FK tests ───────────────────────────────────────────

    #[test]
    fn test_planar_rr_fk_zero_config() {
        let chain = build_planar_rr(1.0, 1.0);
        let pos = chain.end_effector_position();
        // Both joints at 0 → EE at (2, 0, 0)
        assert!((pos[0] - 2.0).abs() < TOL, "x = {}", pos[0]);
        assert!(pos[1].abs() < TOL);
        assert!(pos[2].abs() < TOL);
    }

    #[test]
    fn test_planar_rr_fk_90_degrees() {
        let mut chain = build_planar_rr(1.0, 1.0);
        chain.set_q(&[PI / 2.0, 0.0]);
        let pos = chain.end_effector_position();
        // Joint 1 at 90°: EE at (0, 2, 0)
        assert!(pos[0].abs() < TOL, "x = {}", pos[0]);
        assert!((pos[1] - 2.0).abs() < TOL, "y = {}", pos[1]);
    }

    #[test]
    fn test_planar_rr_fk_folded() {
        let mut chain = build_planar_rr(1.0, 1.0);
        chain.set_q(&[0.0, PI]);
        let pos = chain.end_effector_position();
        // Second link folds back → EE at (0, 0, 0)
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

    // ── Jacobian tests ──────────────────────────────────────────────────

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
        // Numerical differentiation check
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

    // ── Manipulability / singularity tests ──────────────────────────────

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
        chain.set_q(&[0.0, 0.0]); // fully extended along X
        let info = singularity_analysis(&chain, 0.01);
        assert!(
            info.is_singular,
            "fully extended RR should be singular, manip = {}",
            info.manipulability
        );
    }

    // ── Inverse kinematics tests ────────────────────────────────────────

    #[test]
    fn test_ik_reachable_target() {
        let mut chain = build_planar_rr(1.0, 1.0);
        // Start from a non-singular configuration
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
        let target = [10.0, 10.0, 0.0]; // way out of reach
        let ok = chain.inverse_kinematics(&target, 50, 1e-4, 0.01);
        assert!(!ok, "IK should not converge for unreachable target");
    }

    // ── Workspace tests ─────────────────────────────────────────────────

    #[test]
    fn test_workspace_rr_bounding_box() {
        let chain = build_planar_rr(1.0, 1.0);
        let ws = workspace_analysis(&chain, 10);
        // Max reach = 2.0, so the bounding box should span about [-2, 2] in x and y
        assert!(ws.bb_max[0] > 1.5, "bb_max x = {}", ws.bb_max[0]);
        assert!(ws.bb_min[0] < -1.5, "bb_min x = {}", ws.bb_min[0]);
        // Planar chain: z extent is 0, so approx_volume is 0. Check 2D area instead.
        let area = (ws.bb_max[0] - ws.bb_min[0]) * (ws.bb_max[1] - ws.bb_min[1]);
        assert!(area > 0.0, "2D workspace area = {area}");
    }

    #[test]
    fn test_workspace_has_points() {
        let chain = build_planar_rr(1.0, 1.0);
        let ws = workspace_analysis(&chain, 5);
        assert_eq!(ws.points.len(), 25); // 5^2
    }

    // ── Four-bar linkage tests ──────────────────────────────────────────

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
        // Verify point B is on the crank circle
        let rb = (r.point_b[0] * r.point_b[0] + r.point_b[1] * r.point_b[1]).sqrt();
        assert!((rb - 2.0).abs() < TOL, "B radius = {rb}");
    }

    #[test]
    fn test_fourbar_velocity_analysis() {
        let fb = FourBarLinkage::new(4.0, 2.0, 3.0, 3.5);
        let result = fb.velocity_analysis(0.5, 1, 1.0);
        assert!(result.is_some());
        let (omega3, omega4) = result.unwrap();
        // Just check they are finite
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

    // ── Crank-slider tests ──────────────────────────────────────────────

    #[test]
    fn test_crank_slider_tdc() {
        let cs = CrankSlider::new(1.0, 3.0);
        let p_tdc = cs.slider_position(0.0);
        // At θ=0: x_crank = r, slider at r + l = 4.0
        assert!((p_tdc - 4.0).abs() < TOL, "TDC = {p_tdc}");
    }

    #[test]
    fn test_crank_slider_bdc() {
        let cs = CrankSlider::new(1.0, 3.0);
        let p_bdc = cs.slider_position(PI);
        // At θ=π: x_crank = -r, slider at -r + l = 2.0
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

    // ── Geneva mechanism tests ──────────────────────────────────────────

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
        // At φ=0 (sin(0)=0, cos(0)=1):
        // vr = rp * cos(0) / sqrt(d² + rp² - 2*d*rp*sin(0))
        //    = rp / sqrt(d² + rp²)
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

    // ── Gruebler / Kutzbach tests ───────────────────────────────────────

    #[test]
    fn test_gruebler_four_bar() {
        // 4-bar: 4 links, 4 full joints → F = 3*(4-1) - 2*4 = 9 - 8 = 1
        let f = gruebler_planar(4, 4, 0);
        assert_eq!(f, 1);
    }

    #[test]
    fn test_kutzbach_spatial_6dof() {
        // 7 links, 6 revolute joints (each 1 DOF)
        // F = 6*(7-1) - 6*5 = 36 - 30 = 6
        let f = kutzbach_spatial(7, &[1, 1, 1, 1, 1, 1]);
        assert_eq!(f, 6);
    }

    // ── Utility tests ───────────────────────────────────────────────────

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
        // With zero joint angles, both legs reach (2,0,0).
        // Platform point 0 at (0,0,0) → error (0,0,0)-(2,0,0) = large
        assert!(err > 0.0);
    }

    #[test]
    fn test_h_inv_roundtrip() {
        let dh = DHParam::new(0.5, 1.0, 2.0, 0.3);
        let t = dh.to_transform();
        let ti = h_inv(&t);
        let product = h_mul(&t, &ti);
        for i in 0..4 {
            for j in 0..4 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (product[i][j] - expected).abs() < 1e-10,
                    "product[{i}][{j}] = {}",
                    product[i][j]
                );
            }
        }
    }

    #[test]
    fn test_velocity_analysis_zero_config() {
        let mut chain = build_planar_rr(1.0, 1.0);
        chain.set_dq(&[1.0, 0.0]);
        let va = velocity_analysis(&chain);
        // At q=(0,0), joint 1 rotates about Z. EE at (2,0,0).
        // v = ω × r = (0,0,1) × (2,0,0) = (0,2,0)
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
        // Just check it returns a finite position
        for k in 0..3 {
            assert!(pos[k].is_finite(), "pos[{k}] = {}", pos[k]);
        }
    }

    #[test]
    fn test_build_scara() {
        let scara = build_scara(1.0, 1.0, (0.0, 0.5));
        assert_eq!(scara.dof(), 4);
    }
}
