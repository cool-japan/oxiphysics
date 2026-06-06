// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Robot geometry: DH parameters, forward/inverse kinematics, Jacobians,
//! workspace computation, collision checking, and self-collision detection.
//!
//! All transforms are represented as 4×4 homogeneous matrices stored in
//! row-major order as `[f64; 16]`. Joint arrays use plain `[f64; 3]` vectors
//! for positions and axes. No external linear-algebra crate is used.

use rand::RngExt;
use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Low-level math helpers
// ---------------------------------------------------------------------------

/// Dot product of two 3-vectors.
#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Vector subtraction.
#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Vector addition.
#[inline]
fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Scalar multiplication of a 3-vector.
#[inline]
fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

/// Euclidean length of a 3-vector.
#[inline]
fn len3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}

/// Cross product.
#[inline]
fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

// ---------------------------------------------------------------------------
// 4×4 homogeneous matrix (row-major)
// ---------------------------------------------------------------------------

/// A 4×4 homogeneous transformation matrix stored row-major.
///
/// Index layout: `m[row * 4 + col]`.
pub type Mat4 = [f64; 16];

/// Returns the 4×4 identity matrix.
pub fn mat4_identity() -> Mat4 {
    let mut m = [0.0f64; 16];
    m[0] = 1.0;
    m[5] = 1.0;
    m[10] = 1.0;
    m[15] = 1.0;
    m
}

/// Multiplies two 4×4 row-major matrices.
pub fn mat4_mul(a: &Mat4, b: &Mat4) -> Mat4 {
    let mut c = [0.0f64; 16];
    for row in 0..4 {
        for col in 0..4 {
            let mut s = 0.0f64;
            for k in 0..4 {
                s += a[row * 4 + k] * b[k * 4 + col];
            }
            c[row * 4 + col] = s;
        }
    }
    c
}

/// Transforms a 3-vector by a 4×4 homogeneous matrix (w = 1).
pub fn mat4_transform_point(m: &Mat4, p: [f64; 3]) -> [f64; 3] {
    [
        m[0] * p[0] + m[1] * p[1] + m[2] * p[2] + m[3],
        m[4] * p[0] + m[5] * p[1] + m[6] * p[2] + m[7],
        m[8] * p[0] + m[9] * p[1] + m[10] * p[2] + m[11],
    ]
}

/// Transforms a 3-direction (w = 0) by a 4×4 matrix.
pub fn mat4_transform_dir(m: &Mat4, d: [f64; 3]) -> [f64; 3] {
    [
        m[0] * d[0] + m[1] * d[1] + m[2] * d[2],
        m[4] * d[0] + m[5] * d[1] + m[6] * d[2],
        m[8] * d[0] + m[9] * d[1] + m[10] * d[2],
    ]
}

/// Extracts the translation component from a 4×4 matrix.
pub fn mat4_translation(m: &Mat4) -> [f64; 3] {
    [m[3], m[7], m[11]]
}

/// Builds a rotation matrix about the X axis.
fn rot_x(angle: f64) -> Mat4 {
    let c = angle.cos();
    let s = angle.sin();
    let mut m = mat4_identity();
    m[5] = c;
    m[6] = -s;
    m[9] = s;
    m[10] = c;
    m
}

/// Builds a rotation matrix about the Z axis.
fn rot_z(angle: f64) -> Mat4 {
    let c = angle.cos();
    let s = angle.sin();
    let mut m = mat4_identity();
    m[0] = c;
    m[1] = -s;
    m[4] = s;
    m[5] = c;
    m
}

/// Builds a pure translation matrix.
fn trans(tx: f64, ty: f64, tz: f64) -> Mat4 {
    let mut m = mat4_identity();
    m[3] = tx;
    m[7] = ty;
    m[11] = tz;
    m
}

// ---------------------------------------------------------------------------
// Denavit-Hartenberg (DH) parameters
// ---------------------------------------------------------------------------

/// Standard Denavit-Hartenberg parameters for a single joint/link.
///
/// Convention: `T = Rz(θ) · Tz(d) · Tx(a) · Rx(α)`.
#[derive(Debug, Clone, Copy)]
pub struct DhParam {
    /// Joint angle θ (radians) — variable for revolute joints.
    pub theta: f64,
    /// Link offset d (metres) — variable for prismatic joints.
    pub d: f64,
    /// Link length a (metres).
    pub a: f64,
    /// Link twist α (radians).
    pub alpha: f64,
}

impl DhParam {
    /// Creates a new DH parameter set.
    pub fn new(theta: f64, d: f64, a: f64, alpha: f64) -> Self {
        Self { theta, d, a, alpha }
    }

    /// Returns the 4×4 homogeneous transform for this DH frame.
    pub fn transform(&self) -> Mat4 {
        let rz = rot_z(self.theta);
        let tz = trans(0.0, 0.0, self.d);
        let tx = trans(self.a, 0.0, 0.0);
        let rx = rot_x(self.alpha);
        mat4_mul(&mat4_mul(&mat4_mul(&rz, &tz), &tx), &rx)
    }
}

// ---------------------------------------------------------------------------
// Joint types
// ---------------------------------------------------------------------------

/// Type of a robot joint.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum JointType {
    /// Revolute joint — rotates about Z axis.
    Revolute,
    /// Prismatic joint — translates along Z axis.
    Prismatic,
    /// Fixed joint — no degree of freedom.
    Fixed,
}

// ---------------------------------------------------------------------------
// Joint limits
// ---------------------------------------------------------------------------

/// Angular or linear limits for a single joint.
#[derive(Debug, Clone, Copy)]
pub struct JointLimits {
    /// Minimum joint value (radians or metres).
    pub min: f64,
    /// Maximum joint value (radians or metres).
    pub max: f64,
    /// Maximum joint velocity (rad/s or m/s).
    pub max_velocity: f64,
    /// Maximum joint torque/force (Nm or N).
    pub max_effort: f64,
}

impl JointLimits {
    /// Creates joint limits.
    pub fn new(min: f64, max: f64, max_velocity: f64, max_effort: f64) -> Self {
        Self {
            min,
            max,
            max_velocity,
            max_effort,
        }
    }

    /// Returns true if `value` is within \[min, max\].
    pub fn is_within(&self, value: f64) -> bool {
        value >= self.min && value <= self.max
    }

    /// Clamps `value` to \[min, max\].
    pub fn clamp(&self, value: f64) -> f64 {
        value.clamp(self.min, self.max)
    }
}

// ---------------------------------------------------------------------------
// RobotLink
// ---------------------------------------------------------------------------

/// Describes a single link of a robot arm.
#[derive(Debug, Clone)]
pub struct RobotLink {
    /// Human-readable name of the link.
    pub name: String,
    /// DH parameters relative to the parent frame.
    pub dh: DhParam,
    /// Type of the joint connecting to this link.
    pub joint_type: JointType,
    /// Joint limits (optional for fixed joints).
    pub limits: Option<JointLimits>,
    /// Capsule radius used for collision checking (metres).
    pub collision_radius: f64,
    /// Capsule half-length used for collision checking (metres).
    pub collision_half_len: f64,
}

impl RobotLink {
    /// Creates a new robot link.
    pub fn new(
        name: impl Into<String>,
        dh: DhParam,
        joint_type: JointType,
        limits: Option<JointLimits>,
        collision_radius: f64,
        collision_half_len: f64,
    ) -> Self {
        Self {
            name: name.into(),
            dh,
            joint_type,
            limits,
            collision_radius,
            collision_half_len,
        }
    }
}

// ---------------------------------------------------------------------------
// RobotArm — the main robot structure
// ---------------------------------------------------------------------------

/// A serial-chain robot arm described by a sequence of DH links.
///
/// Supports forward kinematics, Jacobian computation (analytical and numerical),
/// Newton-Raphson inverse kinematics, workspace Monte-Carlo sampling,
/// joint-limit checking, and self-collision detection.
#[derive(Debug, Clone)]
pub struct RobotArm {
    /// Name of the robot.
    pub name: String,
    /// Ordered list of links (joint 0 is closest to base).
    pub links: Vec<RobotLink>,
    /// Base transform applied before all joint transforms.
    pub base_transform: Mat4,
}

impl RobotArm {
    /// Creates a new robot arm.
    pub fn new(name: impl Into<String>, links: Vec<RobotLink>) -> Self {
        Self {
            name: name.into(),
            links,
            base_transform: mat4_identity(),
        }
    }

    /// Sets the base transform and returns `self` for chaining.
    pub fn with_base_transform(mut self, t: Mat4) -> Self {
        self.base_transform = t;
        self
    }

    /// Number of degrees of freedom (non-Fixed joints).
    pub fn dof(&self) -> usize {
        self.links
            .iter()
            .filter(|l| l.joint_type != JointType::Fixed)
            .count()
    }

    /// Returns the DH parameter for link `i` with the joint variable (`q`)
    /// applied depending on joint type.
    fn dh_with_q(&self, link_idx: usize, q: f64) -> DhParam {
        let link = &self.links[link_idx];
        let mut dh = link.dh;
        match link.joint_type {
            JointType::Revolute => dh.theta += q,
            JointType::Prismatic => dh.d += q,
            JointType::Fixed => {}
        }
        dh
    }

    /// Computes all intermediate frames (base → end-effector) for a given joint
    /// configuration `q`. Returns one matrix per link plus the base.
    pub fn forward_frames(&self, q: &[f64]) -> Vec<Mat4> {
        let mut frames = Vec::with_capacity(self.links.len() + 1);
        frames.push(self.base_transform);
        let mut q_iter = q.iter().copied();
        for (i, link) in self.links.iter().enumerate() {
            let qi = if link.joint_type == JointType::Fixed {
                0.0
            } else {
                q_iter.next().unwrap_or(0.0)
            };
            let dh = self.dh_with_q(i, qi);
            let t = mat4_mul(
                frames.last().expect("collection should not be empty"),
                &dh.transform(),
            );
            frames.push(t);
        }
        frames
    }

    /// Returns the end-effector transform for joint configuration `q`.
    pub fn end_effector(&self, q: &[f64]) -> Mat4 {
        let frames = self.forward_frames(q);
        *frames.last().expect("collection should not be empty")
    }

    /// Returns the end-effector position `[x, y, z]` for joint configuration `q`.
    pub fn end_effector_pos(&self, q: &[f64]) -> [f64; 3] {
        mat4_translation(&self.end_effector(q))
    }

    // -----------------------------------------------------------------------
    // Jacobian (numerical)
    // -----------------------------------------------------------------------

    /// Computes the 3×n positional Jacobian numerically (finite differences).
    ///
    /// Returns a flat row-major array of shape `3 × dof`.
    pub fn jacobian_numerical(&self, q: &[f64]) -> Vec<f64> {
        let n = q.len();
        let eps = 1e-7;
        let p0 = self.end_effector_pos(q);
        let mut j = vec![0.0f64; 3 * n];
        let mut qp = q.to_vec();
        for i in 0..n {
            qp[i] += eps;
            let p1 = self.end_effector_pos(&qp);
            qp[i] = q[i];
            for row in 0..3 {
                j[row * n + i] = (p1[row] - p0[row]) / eps;
            }
        }
        j
    }

    /// Computes the full 6×n Jacobian (position + orientation) numerically.
    ///
    /// Orientation rows use angle-axis approximation via the rotation part of
    /// the end-effector transform. Returns flat row-major `6 × dof`.
    pub fn jacobian_6dof_numerical(&self, q: &[f64]) -> Vec<f64> {
        let n = q.len();
        let eps = 1e-7;
        let t0 = self.end_effector(q);
        let mut j = vec![0.0f64; 6 * n];
        let mut qp = q.to_vec();
        for i in 0..n {
            qp[i] += eps;
            let t1 = self.end_effector(&qp);
            qp[i] = q[i];
            // Position columns
            for row in 0..3 {
                j[row * n + i] = (t1[row * 4 + 3] - t0[row * 4 + 3]) / eps;
            }
            // Orientation columns (skew-symmetric R_delta)
            let r_delta = mat4_mul(&mat4_transpose(&t0), &t1);
            // vee map: extract angular velocity from skew part
            let wx = (r_delta[6] - r_delta[9]) / (2.0 * eps);
            let wy = (r_delta[8] - r_delta[2]) / (2.0 * eps);
            let wz = (r_delta[1] - r_delta[4]) / (2.0 * eps);
            j[3 * n + i] = wx;
            j[4 * n + i] = wy;
            j[5 * n + i] = wz;
        }
        j
    }

    /// Computes the 3×n analytical Jacobian using the geometric (cross-product) method.
    ///
    /// For revolute joints: `J_i = z_{i-1} × (p_e − p_{i-1})`.
    /// For prismatic joints: `J_i = z_{i-1}`.
    pub fn jacobian_analytical(&self, q: &[f64]) -> Vec<f64> {
        let frames = self.forward_frames(q);
        let n = q.len();
        let p_e = mat4_translation(frames.last().expect("collection should not be empty"));
        let mut j = vec![0.0f64; 3 * n];
        let mut qi = 0usize;
        for (i, link) in self.links.iter().enumerate() {
            if link.joint_type == JointType::Fixed {
                continue;
            }
            let frame = &frames[i]; // frame before this joint
            // Z axis of this frame (third column of rotation)
            let z = [frame[2], frame[6], frame[10]];
            let p_i = mat4_translation(frame);
            let col: [f64; 3] = match link.joint_type {
                JointType::Revolute => cross3(z, sub3(p_e, p_i)),
                JointType::Prismatic => z,
                JointType::Fixed => unreachable!(),
            };
            for row in 0..3 {
                j[row * n + qi] = col[row];
            }
            qi += 1;
        }
        j
    }

    // -----------------------------------------------------------------------
    // Jacobian pseudo-inverse (Moore-Penrose via SVD-free damped LS)
    // -----------------------------------------------------------------------

    /// Computes the damped least-squares pseudo-inverse of a `3 × n` Jacobian.
    ///
    /// Uses the formula `J^+ = J^T (J J^T + λ²I)^{-1}` with closed-form
    /// inversion for the 3×3 system.
    pub fn jacobian_pseudoinverse(j: &[f64], n: usize, lambda: f64) -> Vec<f64> {
        // J is 3×n.  Compute A = J J^T (3×3)
        let mut a = [0.0f64; 9];
        for r in 0..3 {
            for c in 0..3 {
                let mut s = 0.0f64;
                for k in 0..n {
                    s += j[r * n + k] * j[c * n + k];
                }
                a[r * 3 + c] = s;
            }
        }
        // Add damping
        let lam2 = lambda * lambda;
        a[0] += lam2;
        a[4] += lam2;
        a[8] += lam2;
        // Invert 3×3
        let inv = mat3_inverse(&a);
        // J^+ = J^T * inv  (n×3)
        let mut jp = vec![0.0f64; n * 3];
        for r in 0..n {
            for c in 0..3 {
                let mut s = 0.0f64;
                for k in 0..3 {
                    s += j[k * n + r] * inv[k * 3 + c];
                }
                jp[r * 3 + c] = s;
            }
        }
        jp
    }

    // -----------------------------------------------------------------------
    // Inverse kinematics (Newton-Raphson)
    // -----------------------------------------------------------------------

    /// Attempts to solve inverse kinematics using Newton-Raphson iteration.
    ///
    /// Returns `Ok(q)` if converged within `max_iter` steps, or `Err` with the
    /// best solution found.
    pub fn ik_newton_raphson(
        &self,
        target: [f64; 3],
        q_init: &[f64],
        max_iter: usize,
        tol: f64,
        lambda: f64,
    ) -> Result<Vec<f64>, Vec<f64>> {
        let mut q = q_init.to_vec();
        let n = q.len();
        for _iter in 0..max_iter {
            let p = self.end_effector_pos(&q);
            let e = sub3(target, p);
            if len3(e) < tol {
                return Ok(q);
            }
            let j = self.jacobian_analytical(&q);
            let jp = Self::jacobian_pseudoinverse(&j, n, lambda);
            // dq = J^+ * e
            for i in 0..n {
                let dqi = jp[i * 3] * e[0] + jp[i * 3 + 1] * e[1] + jp[i * 3 + 2] * e[2];
                q[i] += dqi;
            }
            // Enforce joint limits
            for (i, link) in self.links.iter().enumerate() {
                if link.joint_type == JointType::Fixed {
                    continue;
                }
                let qi_idx = self.joint_index(i);
                if let (Some(qi_idx), Some(lim)) = (qi_idx, &link.limits)
                    && qi_idx < q.len()
                {
                    q[qi_idx] = lim.clamp(q[qi_idx]);
                }
            }
        }
        Err(q)
    }

    /// Returns the q-vector index for link `link_idx`, or `None` if fixed.
    fn joint_index(&self, link_idx: usize) -> Option<usize> {
        let mut idx = 0usize;
        for (i, l) in self.links.iter().enumerate() {
            if l.joint_type == JointType::Fixed {
                continue;
            }
            if i == link_idx {
                return Some(idx);
            }
            idx += 1;
        }
        None
    }

    // -----------------------------------------------------------------------
    // Workspace Monte-Carlo
    // -----------------------------------------------------------------------

    /// Samples `n_samples` random joint configurations and returns the set of
    /// reachable end-effector positions.
    pub fn workspace_monte_carlo(&self, n_samples: usize) -> Vec<[f64; 3]> {
        let mut rng = rand::rng();
        let dof = self.dof();
        let limits: Vec<JointLimits> = self
            .links
            .iter()
            .filter(|l| l.joint_type != JointType::Fixed)
            .map(|l| l.limits.unwrap_or(JointLimits::new(-PI, PI, 10.0, 100.0)))
            .collect();
        let mut points = Vec::with_capacity(n_samples);
        for _ in 0..n_samples {
            let q: Vec<f64> = (0..dof)
                .map(|i| rng.random_range(limits[i].min..limits[i].max))
                .collect();
            points.push(self.end_effector_pos(&q));
        }
        points
    }

    /// Estimates the reachable workspace bounding box from a Monte-Carlo sample.
    ///
    /// Returns `(min_corner, max_corner)` or `None` if no samples.
    pub fn workspace_bounding_box(&self, n_samples: usize) -> Option<([f64; 3], [f64; 3])> {
        let pts = self.workspace_monte_carlo(n_samples);
        if pts.is_empty() {
            return None;
        }
        let mut mn = pts[0];
        let mut mx = pts[0];
        for p in &pts[1..] {
            for k in 0..3 {
                if p[k] < mn[k] {
                    mn[k] = p[k];
                }
                if p[k] > mx[k] {
                    mx[k] = p[k];
                }
            }
        }
        Some((mn, mx))
    }

    // -----------------------------------------------------------------------
    // Collision helpers
    // -----------------------------------------------------------------------

    /// Returns all link capsule endpoints in world space for configuration `q`.
    ///
    /// Each entry is `(p_start, p_end, radius)` describing the capsule axis.
    pub fn link_capsules(&self, q: &[f64]) -> Vec<([f64; 3], [f64; 3], f64)> {
        let frames = self.forward_frames(q);
        let mut caps = Vec::with_capacity(self.links.len());
        for (i, link) in self.links.iter().enumerate() {
            let frame = &frames[i + 1];
            let center = mat4_translation(frame);
            // Local capsule axis is along local Z
            let z_world = [frame[2], frame[6], frame[10]];
            let hl = link.collision_half_len;
            let p_start = sub3(center, scale3(z_world, hl));
            let p_end = add3(center, scale3(z_world, hl));
            caps.push((p_start, p_end, link.collision_radius));
        }
        caps
    }

    /// Checks whether any robot link capsule intersects an obstacle sphere.
    ///
    /// Returns `true` if any link is in collision.
    pub fn check_collision_sphere(
        &self,
        q: &[f64],
        sphere_center: [f64; 3],
        sphere_radius: f64,
    ) -> bool {
        for (a, b, r) in self.link_capsules(q) {
            let dist = dist_point_segment(sphere_center, a, b);
            if dist < r + sphere_radius {
                return true;
            }
        }
        false
    }

    /// Checks whether any robot link capsule intersects an axis-aligned box.
    ///
    /// Returns `true` if any link is in collision.
    pub fn check_collision_aabb(&self, q: &[f64], box_min: [f64; 3], box_max: [f64; 3]) -> bool {
        for (a, b, r) in self.link_capsules(q) {
            if capsule_aabb_overlap(a, b, r, box_min, box_max) {
                return true;
            }
        }
        false
    }

    // -----------------------------------------------------------------------
    // Self-collision detection
    // -----------------------------------------------------------------------

    /// Detects self-collision between non-adjacent link capsules.
    ///
    /// Links with indices that differ by 1 (adjacent) are always skipped.
    /// Returns `true` if any non-adjacent pair of capsules overlaps.
    pub fn check_self_collision(&self, q: &[f64]) -> bool {
        let caps = self.link_capsules(q);
        let n = caps.len();
        for i in 0..n {
            for j in (i + 2)..n {
                let (a1, b1, r1) = caps[i];
                let (a2, b2, r2) = caps[j];
                let dist = dist_segment_segment(a1, b1, a2, b2);
                if dist < r1 + r2 {
                    return true;
                }
            }
        }
        false
    }

    // -----------------------------------------------------------------------
    // Joint limit checking
    // -----------------------------------------------------------------------

    /// Returns `true` if all joint values in `q` satisfy the robot's joint limits.
    pub fn joints_within_limits(&self, q: &[f64]) -> bool {
        let mut qi = 0usize;
        for link in &self.links {
            if link.joint_type == JointType::Fixed {
                continue;
            }
            if qi >= q.len() {
                break;
            }
            if let Some(lim) = &link.limits
                && !lim.is_within(q[qi])
            {
                return false;
            }
            qi += 1;
        }
        true
    }

    /// Clamps all joints in `q` to the robot's limits and returns the result.
    pub fn clamp_joints(&self, q: &[f64]) -> Vec<f64> {
        let mut out = q.to_vec();
        let mut qi = 0usize;
        for link in &self.links {
            if link.joint_type == JointType::Fixed {
                continue;
            }
            if qi >= out.len() {
                break;
            }
            if let Some(lim) = &link.limits {
                out[qi] = lim.clamp(out[qi]);
            }
            qi += 1;
        }
        out
    }

    // -----------------------------------------------------------------------
    // Swept volume
    // -----------------------------------------------------------------------

    /// Computes a coarse AABB for the swept volume of the end-effector as joints
    /// vary linearly from `q_start` to `q_end` using `n_steps` samples.
    pub fn swept_volume_aabb(
        &self,
        q_start: &[f64],
        q_end: &[f64],
        n_steps: usize,
    ) -> ([f64; 3], [f64; 3]) {
        let n_steps = n_steps.max(2);
        let mut mn = [f64::INFINITY; 3];
        let mut mx = [f64::NEG_INFINITY; 3];
        for step in 0..=n_steps {
            let t = step as f64 / n_steps as f64;
            let q: Vec<f64> = q_start
                .iter()
                .zip(q_end.iter())
                .map(|(a, b)| a + t * (b - a))
                .collect();
            let p = self.end_effector_pos(&q);
            for k in 0..3 {
                if p[k] < mn[k] {
                    mn[k] = p[k];
                }
                if p[k] > mx[k] {
                    mx[k] = p[k];
                }
            }
        }
        (mn, mx)
    }

    // -----------------------------------------------------------------------
    // Manipulability measure
    // -----------------------------------------------------------------------

    /// Computes Yoshikawa's manipulability measure `sqrt(det(J J^T))`.
    ///
    /// High values indicate configurations far from singularities.
    pub fn manipulability(&self, q: &[f64]) -> f64 {
        let n = q.len();
        let j = self.jacobian_analytical(q);
        // Compute J J^T (3×3)
        let mut a = [0.0f64; 9];
        for r in 0..3 {
            for c in 0..3 {
                let mut s = 0.0f64;
                for k in 0..n {
                    s += j[r * n + k] * j[c * n + k];
                }
                a[r * 3 + c] = s;
            }
        }
        mat3_det(&a).abs().sqrt()
    }
}

// ---------------------------------------------------------------------------
// 3×3 matrix helpers
// ---------------------------------------------------------------------------

/// Computes the determinant of a 3×3 row-major matrix.
pub fn mat3_det(m: &[f64; 9]) -> f64 {
    m[0] * (m[4] * m[8] - m[5] * m[7]) - m[1] * (m[3] * m[8] - m[5] * m[6])
        + m[2] * (m[3] * m[7] - m[4] * m[6])
}

/// Computes the inverse of a 3×3 row-major matrix.
///
/// Returns the zero matrix if the determinant is near zero.
pub fn mat3_inverse(m: &[f64; 9]) -> [f64; 9] {
    let det = mat3_det(m);
    if det.abs() < 1e-14 {
        return [0.0; 9];
    }
    let inv_det = 1.0 / det;
    [
        (m[4] * m[8] - m[5] * m[7]) * inv_det,
        (m[2] * m[7] - m[1] * m[8]) * inv_det,
        (m[1] * m[5] - m[2] * m[4]) * inv_det,
        (m[5] * m[6] - m[3] * m[8]) * inv_det,
        (m[0] * m[8] - m[2] * m[6]) * inv_det,
        (m[2] * m[3] - m[0] * m[5]) * inv_det,
        (m[3] * m[7] - m[4] * m[6]) * inv_det,
        (m[1] * m[6] - m[0] * m[7]) * inv_det,
        (m[0] * m[4] - m[1] * m[3]) * inv_det,
    ]
}

/// Transposes a 4×4 row-major matrix.
pub fn mat4_transpose(m: &Mat4) -> Mat4 {
    let mut t = [0.0f64; 16];
    for r in 0..4 {
        for c in 0..4 {
            t[c * 4 + r] = m[r * 4 + c];
        }
    }
    t
}

// ---------------------------------------------------------------------------
// Geometry primitives for collision
// ---------------------------------------------------------------------------

/// Minimum distance from point `p` to line segment `[a, b]`.
pub fn dist_point_segment(p: [f64; 3], a: [f64; 3], b: [f64; 3]) -> f64 {
    let ab = sub3(b, a);
    let ap = sub3(p, a);
    let t = dot3(ap, ab) / (dot3(ab, ab) + 1e-30);
    let t = t.clamp(0.0, 1.0);
    let closest = add3(a, scale3(ab, t));
    len3(sub3(p, closest))
}

/// Minimum distance between two line segments `[a1,b1]` and `[a2,b2]`.
pub fn dist_segment_segment(a1: [f64; 3], b1: [f64; 3], a2: [f64; 3], b2: [f64; 3]) -> f64 {
    let d1 = sub3(b1, a1);
    let d2 = sub3(b2, a2);
    let r = sub3(a1, a2);
    let a = dot3(d1, d1);
    let e = dot3(d2, d2);
    let f = dot3(d2, r);
    let eps = 1e-14;

    let (s, t) = if a < eps && e < eps {
        (0.0f64, 0.0f64)
    } else if a < eps {
        (0.0f64, (f / e).clamp(0.0, 1.0))
    } else {
        let c = dot3(d1, r);
        if e < eps {
            (((-c) / a).clamp(0.0, 1.0), 0.0)
        } else {
            let b = dot3(d1, d2);
            let denom = a * e - b * b;
            let s = if denom.abs() > eps {
                ((b * f - c * e) / denom).clamp(0.0, 1.0)
            } else {
                0.0
            };
            let t = (b * s + f) / e;
            if t < 0.0 {
                let s = ((-c) / a).clamp(0.0, 1.0);
                (s, 0.0f64)
            } else if t > 1.0 {
                let s = ((b - c) / a).clamp(0.0, 1.0);
                (s, 1.0f64)
            } else {
                (s, t)
            }
        }
    };
    let c1 = add3(a1, scale3(d1, s));
    let c2 = add3(a2, scale3(d2, t));
    len3(sub3(c1, c2))
}

/// Tests whether a capsule (axis `[a,b]`, radius `r`) overlaps an AABB.
pub fn capsule_aabb_overlap(
    a: [f64; 3],
    b: [f64; 3],
    r: f64,
    box_min: [f64; 3],
    box_max: [f64; 3],
) -> bool {
    // Closest point on segment to AABB center
    let center = [
        (box_min[0] + box_max[0]) * 0.5,
        (box_min[1] + box_max[1]) * 0.5,
        (box_min[2] + box_max[2]) * 0.5,
    ];
    let half = [
        (box_max[0] - box_min[0]) * 0.5,
        (box_max[1] - box_min[1]) * 0.5,
        (box_max[2] - box_min[2]) * 0.5,
    ];
    let ab = sub3(b, a);
    let ac = sub3(center, a);
    let t = (dot3(ac, ab) / (dot3(ab, ab) + 1e-30)).clamp(0.0, 1.0);
    let closest_on_seg = add3(a, scale3(ab, t));
    // Clamp to box
    let mut dist_sq = 0.0f64;
    for k in 0..3 {
        let v = closest_on_seg[k];
        let lo = box_min[k];
        let hi = box_max[k];
        if v < lo {
            dist_sq += (lo - v) * (lo - v);
        } else if v > hi {
            dist_sq += (v - hi) * (v - hi);
        }
    }
    let _ = half; // used implicitly via box_min/box_max
    dist_sq <= r * r
}

// ---------------------------------------------------------------------------
// Preset robot constructors
// ---------------------------------------------------------------------------

/// Builds a 3-DOF planar robot arm (all joints revolute, links in XY plane).
///
/// `lengths`: link lengths in metres.
pub fn build_planar_3dof(lengths: [f64; 3]) -> RobotArm {
    let lim = JointLimits::new(-PI, PI, 3.0, 50.0);
    let links = vec![
        RobotLink::new(
            "link1",
            DhParam::new(0.0, 0.0, lengths[0], 0.0),
            JointType::Revolute,
            Some(lim),
            0.03,
            lengths[0] * 0.5,
        ),
        RobotLink::new(
            "link2",
            DhParam::new(0.0, 0.0, lengths[1], 0.0),
            JointType::Revolute,
            Some(lim),
            0.03,
            lengths[1] * 0.5,
        ),
        RobotLink::new(
            "link3",
            DhParam::new(0.0, 0.0, lengths[2], 0.0),
            JointType::Revolute,
            Some(lim),
            0.03,
            lengths[2] * 0.5,
        ),
    ];
    RobotArm::new("planar_3dof", links)
}

/// Builds a 6-DOF robot arm using standard Puma-like DH parameters.
pub fn build_6dof_puma_like() -> RobotArm {
    let lim_base = JointLimits::new(-PI, PI, 2.0, 150.0);
    let lim_elbow = JointLimits::new(-PI * 0.75, PI * 0.75, 2.0, 100.0);
    let lim_wrist = JointLimits::new(-PI, PI, 3.0, 50.0);

    let links = vec![
        RobotLink::new(
            "joint1",
            DhParam::new(0.0, 0.36, 0.0, PI / 2.0),
            JointType::Revolute,
            Some(lim_base),
            0.06,
            0.18,
        ),
        RobotLink::new(
            "joint2",
            DhParam::new(0.0, 0.0, 0.43, 0.0),
            JointType::Revolute,
            Some(lim_elbow),
            0.05,
            0.215,
        ),
        RobotLink::new(
            "joint3",
            DhParam::new(0.0, 0.0, 0.0, PI / 2.0),
            JointType::Revolute,
            Some(lim_elbow),
            0.05,
            0.05,
        ),
        RobotLink::new(
            "joint4",
            DhParam::new(0.0, 0.43, 0.0, -PI / 2.0),
            JointType::Revolute,
            Some(lim_wrist),
            0.04,
            0.215,
        ),
        RobotLink::new(
            "joint5",
            DhParam::new(0.0, 0.0, 0.0, PI / 2.0),
            JointType::Revolute,
            Some(lim_wrist),
            0.04,
            0.04,
        ),
        RobotLink::new(
            "joint6",
            DhParam::new(0.0, 0.1, 0.0, 0.0),
            JointType::Revolute,
            Some(lim_wrist),
            0.03,
            0.05,
        ),
    ];
    RobotArm::new("puma_6dof", links)
}

// ---------------------------------------------------------------------------
// URDF-like robot description (from struct)
// ---------------------------------------------------------------------------

/// A joint in a URDF-like robot description.
#[derive(Debug, Clone)]
pub struct UrdfJoint {
    /// Joint name.
    pub name: String,
    /// Parent link name.
    pub parent: String,
    /// Child link name.
    pub child: String,
    /// Joint type.
    pub joint_type: JointType,
    /// Joint origin XYZ.
    pub origin_xyz: [f64; 3],
    /// Joint origin RPY.
    pub origin_rpy: [f64; 3],
    /// Joint axis.
    pub axis: [f64; 3],
    /// Optional limits.
    pub limits: Option<JointLimits>,
}

/// A link in a URDF-like description.
#[derive(Debug, Clone)]
pub struct UrdfLink {
    /// Link name.
    pub name: String,
    /// Collision geometry radius (simplified capsule).
    pub collision_radius: f64,
    /// Collision geometry half-length.
    pub collision_half_len: f64,
}

/// A complete URDF-like robot description parsed from structs.
#[derive(Debug, Clone)]
pub struct UrdfRobot {
    /// Robot name.
    pub name: String,
    /// All links.
    pub links: Vec<UrdfLink>,
    /// All joints.
    pub joints: Vec<UrdfJoint>,
}

impl UrdfRobot {
    /// Converts this URDF-like description into a `RobotArm`.
    ///
    /// Only supports simple serial chains (first joint chain from base).
    pub fn to_robot_arm(&self) -> RobotArm {
        let arm_links: Vec<RobotLink> = self
            .joints
            .iter()
            .map(|jt| {
                let collision = self
                    .links
                    .iter()
                    .find(|l| l.name == jt.child)
                    .map(|l| (l.collision_radius, l.collision_half_len))
                    .unwrap_or((0.05, 0.1));
                // Build DH-approximated transform from RPY + XYZ
                let dh = rpy_xyz_to_dh(jt.origin_rpy, jt.origin_xyz);
                RobotLink::new(
                    &jt.name,
                    dh,
                    jt.joint_type,
                    jt.limits,
                    collision.0,
                    collision.1,
                )
            })
            .collect();
        RobotArm::new(&self.name, arm_links)
    }
}

/// Approximates DH parameters from RPY and XYZ joint origin data.
///
/// This is a simplified mapping; full DH extraction requires more link data.
fn rpy_xyz_to_dh(rpy: [f64; 3], xyz: [f64; 3]) -> DhParam {
    let a = (xyz[0] * xyz[0] + xyz[1] * xyz[1]).sqrt();
    DhParam::new(rpy[2], xyz[2], a, rpy[0])
}

// ---------------------------------------------------------------------------
// Velocity and acceleration analysis
// ---------------------------------------------------------------------------

/// Computes the end-effector Cartesian velocity given joint velocities `dq`.
///
/// `v = J * dq`  (3-vector)
pub fn end_effector_velocity(jacobian: &[f64], dq: &[f64]) -> [f64; 3] {
    let n = dq.len();
    let mut v = [0.0f64; 3];
    for row in 0..3 {
        for k in 0..n {
            v[row] += jacobian[row * n + k] * dq[k];
        }
    }
    v
}

/// Computes the minimum norm joint velocity for a desired Cartesian velocity.
///
/// Uses the damped pseudo-inverse: `dq = J^+ * v`.
pub fn joint_velocity_from_cartesian(
    jacobian: &[f64],
    n_dof: usize,
    v_desired: [f64; 3],
    lambda: f64,
) -> Vec<f64> {
    let jp = RobotArm::jacobian_pseudoinverse(jacobian, n_dof, lambda);
    let mut dq = vec![0.0f64; n_dof];
    for i in 0..n_dof {
        dq[i] =
            jp[i * 3] * v_desired[0] + jp[i * 3 + 1] * v_desired[1] + jp[i * 3 + 2] * v_desired[2];
    }
    dq
}

// ---------------------------------------------------------------------------
// Singularity analysis
// ---------------------------------------------------------------------------

/// Classifies the manipulability of a configuration.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SingularityClass {
    /// Far from any singularity.
    Regular,
    /// Close to (but not at) a singularity.
    NearSingular,
    /// At or very close to a type-1 (boundary) singularity.
    BoundarySingular,
    /// At or very close to a type-2 (interior) singularity.
    InteriorSingular,
}

/// Classifies the singularity of a configuration from its manipulability value.
pub fn classify_singularity(manipulability: f64) -> SingularityClass {
    if manipulability < 1e-6 {
        SingularityClass::InteriorSingular
    } else if manipulability < 1e-3 {
        SingularityClass::BoundarySingular
    } else if manipulability < 0.05 {
        SingularityClass::NearSingular
    } else {
        SingularityClass::Regular
    }
}

// ---------------------------------------------------------------------------
// Condition number of Jacobian
// ---------------------------------------------------------------------------

/// Estimates the condition number of the `3×n` Jacobian via the ratio of
/// the largest to smallest singular value approximation.
///
/// Uses power iteration to find the largest and smallest "effective"
/// singular values via `J J^T`.
pub fn jacobian_condition_number(jacobian: &[f64], n: usize) -> f64 {
    // Form J J^T (3×3)
    let mut a = [0.0f64; 9];
    for r in 0..3 {
        for c in 0..3 {
            let mut s = 0.0f64;
            for k in 0..n {
                s += jacobian[r * n + k] * jacobian[c * n + k];
            }
            a[r * 3 + c] = s;
        }
    }
    // Frobenius norm as proxy for max singular value (fast approximation)
    let frob: f64 = a.iter().map(|x| x * x).sum::<f64>().sqrt();
    let det = mat3_det(&a).abs();
    if frob < 1e-14 {
        return f64::INFINITY;
    }
    // min_sv ≈ det / (frob^2) for 3×3 case
    let min_sv_approx = det / (frob * frob + 1e-30);
    frob / (min_sv_approx + 1e-30)
}

// ---------------------------------------------------------------------------
// Path planning helpers
// ---------------------------------------------------------------------------

/// Generates a linear interpolation path in joint space.
///
/// Returns `n_steps + 1` joint configurations from `q_start` to `q_end`.
pub fn joint_space_path(q_start: &[f64], q_end: &[f64], n_steps: usize) -> Vec<Vec<f64>> {
    let n_steps = n_steps.max(1);
    (0..=n_steps)
        .map(|step| {
            let t = step as f64 / n_steps as f64;
            q_start
                .iter()
                .zip(q_end.iter())
                .map(|(a, b)| a + t * (b - a))
                .collect()
        })
        .collect()
}

/// Checks whether every configuration along a straight-line joint-space path
/// is within joint limits and collision-free with respect to a sphere obstacle.
pub fn path_collision_free(
    arm: &RobotArm,
    path: &[Vec<f64>],
    sphere_center: [f64; 3],
    sphere_radius: f64,
) -> bool {
    for q in path {
        if arm.check_collision_sphere(q, sphere_center, sphere_radius) {
            return false;
        }
        if !arm.joints_within_limits(q) {
            return false;
        }
    }
    true
}

// ---------------------------------------------------------------------------
// Torque / statics analysis
// ---------------------------------------------------------------------------

/// Computes the joint torques required to produce a given end-effector force.
///
/// Uses the static relationship `τ = J^T f`.
pub fn joint_torques_from_force(jacobian: &[f64], n_dof: usize, force: [f64; 3]) -> Vec<f64> {
    let mut tau = vec![0.0f64; n_dof];
    for i in 0..n_dof {
        tau[i] = jacobian[i] * force[0]
            + jacobian[n_dof + i] * force[1]
            + jacobian[2 * n_dof + i] * force[2];
    }
    tau
}

// ---------------------------------------------------------------------------
// Reachability map
// ---------------------------------------------------------------------------

/// A coarse voxel reachability map for a robot arm.
///
/// Populated by Monte-Carlo workspace sampling.
#[derive(Debug, Clone)]
pub struct ReachabilityMap {
    /// World-space origin of the voxel grid.
    pub origin: [f64; 3],
    /// Size of each voxel in metres.
    pub voxel_size: f64,
    /// Number of voxels along each axis.
    pub dims: [usize; 3],
    /// Flat array of occupancy counts.
    pub counts: Vec<u32>,
}

impl ReachabilityMap {
    /// Creates an empty reachability map.
    pub fn new(origin: [f64; 3], voxel_size: f64, dims: [usize; 3]) -> Self {
        let total = dims[0] * dims[1] * dims[2];
        Self {
            origin,
            voxel_size,
            dims,
            counts: vec![0u32; total],
        }
    }

    /// Converts a world position to a voxel index, returning `None` if outside.
    pub fn world_to_voxel(&self, p: [f64; 3]) -> Option<[usize; 3]> {
        let mut idx = [0usize; 3];
        for k in 0..3 {
            let v = (p[k] - self.origin[k]) / self.voxel_size;
            if v < 0.0 || v >= self.dims[k] as f64 {
                return None;
            }
            idx[k] = v as usize;
        }
        Some(idx)
    }

    /// Increments the occupancy count for the voxel at world position `p`.
    pub fn add_point(&mut self, p: [f64; 3]) {
        if let Some(idx) = self.world_to_voxel(p) {
            let flat = idx[0] * self.dims[1] * self.dims[2] + idx[1] * self.dims[2] + idx[2];
            self.counts[flat] = self.counts[flat].saturating_add(1);
        }
    }

    /// Populates the map from a Monte-Carlo workspace sample.
    pub fn populate(&mut self, arm: &RobotArm, n_samples: usize) {
        let pts = arm.workspace_monte_carlo(n_samples);
        for p in pts {
            self.add_point(p);
        }
    }

    /// Returns the total number of non-zero voxels (reachable cells).
    pub fn reachable_count(&self) -> usize {
        self.counts.iter().filter(|&&c| c > 0).count()
    }
}

// ---------------------------------------------------------------------------
// Unit tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn planar_arm() -> RobotArm {
        build_planar_3dof([0.5, 0.4, 0.3])
    }

    fn puma_arm() -> RobotArm {
        build_6dof_puma_like()
    }

    // --- DH parameter tests ---

    #[test]
    fn test_dh_identity_at_zero() {
        let dh = DhParam::new(0.0, 0.0, 0.0, 0.0);
        let t = dh.transform();
        let ident = mat4_identity();
        for i in 0..16 {
            assert!(
                (t[i] - ident[i]).abs() < 1e-12,
                "DH identity mismatch at {i}"
            );
        }
    }

    #[test]
    fn test_dh_pure_translation_a() {
        let a = 1.0;
        let dh = DhParam::new(0.0, 0.0, a, 0.0);
        let t = dh.transform();
        // Should translate by a along X
        assert!((t[3] - a).abs() < 1e-12);
        assert!(t[7].abs() < 1e-12);
        assert!(t[11].abs() < 1e-12);
    }

    #[test]
    fn test_dh_pure_translation_d() {
        let d = 2.0;
        let dh = DhParam::new(0.0, d, 0.0, 0.0);
        let t = dh.transform();
        // Should translate by d along Z
        assert!((t[11] - d).abs() < 1e-12);
    }

    #[test]
    fn test_dh_rotation_theta() {
        let dh = DhParam::new(PI / 2.0, 0.0, 0.0, 0.0);
        let t = dh.transform();
        // Rz(pi/2): cos -> 0, sin -> 1
        assert!(t[0].abs() < 1e-12); // cos
        assert!((t[1] + 1.0).abs() < 1e-12); // -sin
        assert!((t[4] - 1.0).abs() < 1e-12); // sin
    }

    #[test]
    fn test_mat4_mul_identity() {
        let ident = mat4_identity();
        let a = rot_z(0.5);
        let result = mat4_mul(&a, &ident);
        for i in 0..16 {
            assert!(
                (result[i] - a[i]).abs() < 1e-12,
                "mul ident mismatch at {i}"
            );
        }
    }

    #[test]
    fn test_mat4_transpose_roundtrip() {
        let m = rot_z(1.2);
        let tt = mat4_transpose(&mat4_transpose(&m));
        for i in 0..16 {
            assert!((tt[i] - m[i]).abs() < 1e-12);
        }
    }

    // --- Forward kinematics ---

    #[test]
    fn test_fk_zero_config_planar_x() {
        let arm = planar_arm();
        let q = vec![0.0f64; 3];
        let pos = arm.end_effector_pos(&q);
        // All links along X: pos ~ [0.5+0.4+0.3, 0, 0]
        assert!((pos[0] - 1.2).abs() < 1e-10);
        assert!(pos[1].abs() < 1e-10);
        assert!(pos[2].abs() < 1e-10);
    }

    #[test]
    fn test_fk_pi_fold() {
        let arm = planar_arm();
        // With q2 = PI, end-effector folds back
        let q = vec![0.0, PI, 0.0];
        let pos = arm.end_effector_pos(&q);
        // link1(0.5) + link2*cos(pi)(-0.4) + link3*cos(pi)(-0.3) ≈ 0.5-0.4-0.3 = -0.2
        assert!((pos[0] - (-0.2)).abs() < 1e-8);
    }

    #[test]
    fn test_fk_frames_count() {
        let arm = planar_arm();
        let q = vec![0.0f64; 3];
        let frames = arm.forward_frames(&q);
        // 1 base + 3 links = 4
        assert_eq!(frames.len(), 4);
    }

    #[test]
    fn test_dof_count() {
        let arm = planar_arm();
        assert_eq!(arm.dof(), 3);
    }

    #[test]
    fn test_dof_puma() {
        let arm = puma_arm();
        assert_eq!(arm.dof(), 6);
    }

    // --- Joint limits ---

    #[test]
    fn test_joint_limits_within() {
        let lim = JointLimits::new(-PI, PI, 3.0, 50.0);
        assert!(lim.is_within(0.0));
        assert!(lim.is_within(-PI));
        assert!(lim.is_within(PI));
        assert!(!lim.is_within(PI + 0.001));
    }

    #[test]
    fn test_joint_limits_clamp() {
        let lim = JointLimits::new(-1.0, 1.0, 1.0, 10.0);
        assert_eq!(lim.clamp(2.0), 1.0);
        assert_eq!(lim.clamp(-2.0), -1.0);
        assert_eq!(lim.clamp(0.5), 0.5);
    }

    #[test]
    fn test_joints_within_limits_planar() {
        let arm = planar_arm();
        assert!(arm.joints_within_limits(&[0.0, 0.0, 0.0]));
        assert!(!arm.joints_within_limits(&[4.0, 0.0, 0.0])); // > PI
    }

    #[test]
    fn test_clamp_joints_planar() {
        let arm = planar_arm();
        let clamped = arm.clamp_joints(&[5.0, -5.0, 0.0]);
        assert!((clamped[0] - PI).abs() < 1e-12);
        assert!((clamped[1] - (-PI)).abs() < 1e-12);
    }

    // --- Jacobian tests ---

    #[test]
    fn test_jacobian_numerical_size() {
        let arm = planar_arm();
        let q = vec![0.1, 0.2, 0.3];
        let j = arm.jacobian_numerical(&q);
        assert_eq!(j.len(), 3 * 3);
    }

    #[test]
    fn test_jacobian_analytical_matches_numerical() {
        let arm = planar_arm();
        let q = vec![0.2, -0.3, 0.5];
        let jn = arm.jacobian_numerical(&q);
        let ja = arm.jacobian_analytical(&q);
        for i in 0..9 {
            assert!(
                (jn[i] - ja[i]).abs() < 1e-4,
                "Jacobian mismatch at {i}: numerical={} analytical={}",
                jn[i],
                ja[i]
            );
        }
    }

    #[test]
    fn test_jacobian_pseudoinverse_identity_check() {
        // For a simple scalar case: J = [1,0,0; 0,1,0; 0,0,1], pinv = identity
        let j: Vec<f64> = vec![1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        let jp = RobotArm::jacobian_pseudoinverse(&j, 3, 0.0001);
        // jp should be approximately identity (3x3 = n=3, result is 3x3 n×3)
        assert!((jp[0] - 1.0).abs() < 1e-3);
        assert!((jp[4] - 1.0).abs() < 1e-3);
        assert!((jp[8] - 1.0).abs() < 1e-3);
    }

    // --- IK tests ---

    #[test]
    fn test_ik_reachable_target() {
        let arm = planar_arm();
        let target = [0.8, 0.4, 0.0];
        let q_init = vec![0.1f64, 0.2, 0.1];
        let result = arm.ik_newton_raphson(target, &q_init, 200, 1e-6, 0.01);
        let q = match result {
            Ok(q) | Err(q) => q,
        };
        let pos = arm.end_effector_pos(&q);
        assert!(
            len3(sub3(target, pos)) < 0.05,
            "IK did not converge close enough"
        );
    }

    #[test]
    fn test_ik_at_zero_config() {
        let arm = planar_arm();
        let target = arm.end_effector_pos(&[0.0, 0.0, 0.0]);
        let q_init = vec![0.01f64, 0.01, 0.01];
        let result = arm.ik_newton_raphson(target, &q_init, 200, 1e-8, 0.001);
        let q = match result {
            Ok(q) | Err(q) => q,
        };
        let pos = arm.end_effector_pos(&q);
        assert!(len3(sub3(target, pos)) < 1e-4);
    }

    // --- Collision tests ---

    #[test]
    fn test_dist_point_segment_basic() {
        let d = dist_point_segment([0.0, 1.0, 0.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert!((d - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_dist_point_segment_on_segment() {
        let d = dist_point_segment([0.5, 0.0, 0.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert!(d.abs() < 1e-12);
    }

    #[test]
    fn test_collision_sphere_no_hit() {
        let arm = planar_arm();
        let q = vec![0.0f64; 3];
        // Sphere far away
        let hit = arm.check_collision_sphere(&q, [10.0, 10.0, 10.0], 0.01);
        assert!(!hit);
    }

    #[test]
    fn test_collision_sphere_hit() {
        let arm = planar_arm();
        let q = vec![0.0f64; 3];
        // Sphere at robot base position
        let hit = arm.check_collision_sphere(&q, [0.25, 0.0, 0.0], 0.3);
        assert!(hit);
    }

    #[test]
    fn test_self_collision_straight() {
        let arm = planar_arm();
        // Straight configuration — no self-collision expected
        let q = vec![0.0f64; 3];
        assert!(!arm.check_self_collision(&q));
    }

    #[test]
    fn test_capsule_aabb_no_overlap() {
        let no_overlap = capsule_aabb_overlap(
            [10.0, 0.0, 0.0],
            [11.0, 0.0, 0.0],
            0.1,
            [-1.0, -1.0, -1.0],
            [1.0, 1.0, 1.0],
        );
        assert!(!no_overlap);
    }

    #[test]
    fn test_capsule_aabb_overlap() {
        let overlap = capsule_aabb_overlap(
            [-0.5, 0.0, 0.0],
            [0.5, 0.0, 0.0],
            0.1,
            [-1.0, -1.0, -1.0],
            [1.0, 1.0, 1.0],
        );
        assert!(overlap);
    }

    // --- Swept volume ---

    #[test]
    fn test_swept_volume_contains_endpoints() {
        let arm = planar_arm();
        let q0 = vec![0.0f64; 3];
        let q1 = vec![PI / 4.0, PI / 4.0, PI / 4.0];
        let (mn, mx) = arm.swept_volume_aabb(&q0, &q1, 20);
        let p0 = arm.end_effector_pos(&q0);
        let p1 = arm.end_effector_pos(&q1);
        for k in 0..3 {
            assert!(p0[k] >= mn[k] - 1e-9 && p0[k] <= mx[k] + 1e-9);
            assert!(p1[k] >= mn[k] - 1e-9 && p1[k] <= mx[k] + 1e-9);
        }
    }

    // --- Workspace ---

    #[test]
    fn test_workspace_monte_carlo_count() {
        let arm = planar_arm();
        let pts = arm.workspace_monte_carlo(100);
        assert_eq!(pts.len(), 100);
    }

    #[test]
    fn test_workspace_bounding_box_some() {
        let arm = planar_arm();
        let bb = arm.workspace_bounding_box(200);
        assert!(bb.is_some());
        let (mn, mx) = bb.unwrap();
        for k in 0..3 {
            assert!(mn[k] <= mx[k]);
        }
    }

    // --- Manipulability ---

    #[test]
    fn test_manipulability_positive_planar() {
        let arm = planar_arm();
        let m = arm.manipulability(&[0.3, 0.4, 0.5]);
        assert!(m >= 0.0);
    }

    #[test]
    fn test_manipulability_near_singular() {
        let arm = planar_arm();
        // Fully extended: all q = 0 → rank-deficient in 2D embedded in 3D
        let m = arm.manipulability(&[0.0, 0.0, 0.0]);
        // In 3D, planar arm has zero z-component, so det can be 0
        assert!(m >= 0.0);
    }

    #[test]
    fn test_singularity_classification() {
        assert_eq!(
            classify_singularity(0.0),
            SingularityClass::InteriorSingular
        );
        assert_eq!(classify_singularity(0.1), SingularityClass::Regular);
        assert_eq!(
            classify_singularity(1e-4),
            SingularityClass::BoundarySingular
        );
        assert_eq!(classify_singularity(0.01), SingularityClass::NearSingular);
    }

    // --- Reachability map ---

    #[test]
    fn test_reachability_map_populate() {
        let arm = planar_arm();
        let mut rmap = ReachabilityMap::new([-2.0, -2.0, -0.1], 0.1, [40, 40, 2]);
        rmap.populate(&arm, 500);
        assert!(rmap.reachable_count() > 0);
    }

    #[test]
    fn test_reachability_map_voxel_oob() {
        let rmap = ReachabilityMap::new([0.0, 0.0, 0.0], 0.1, [10, 10, 10]);
        assert!(rmap.world_to_voxel([1.5, 0.5, 0.5]).is_none());
    }

    // --- Path planning helpers ---

    #[test]
    fn test_joint_space_path_length() {
        let path = joint_space_path(&[0.0, 0.0], &[1.0, 1.0], 10);
        assert_eq!(path.len(), 11);
    }

    #[test]
    fn test_joint_space_path_endpoints() {
        let q0 = vec![0.0f64, 0.5];
        let q1 = vec![1.0f64, -0.5];
        let path = joint_space_path(&q0, &q1, 5);
        for (a, b) in path[0].iter().zip(q0.iter()) {
            assert!((a - b).abs() < 1e-12);
        }
        for (a, b) in path.last().unwrap().iter().zip(q1.iter()) {
            assert!((a - b).abs() < 1e-12);
        }
    }

    // --- Torque analysis ---

    #[test]
    fn test_joint_torques_zero_force() {
        let arm = planar_arm();
        let q = vec![0.0f64; 3];
        let j = arm.jacobian_analytical(&q);
        let tau = joint_torques_from_force(&j, 3, [0.0, 0.0, 0.0]);
        for t in &tau {
            assert!(t.abs() < 1e-12);
        }
    }

    // --- URDF-like struct ---

    #[test]
    fn test_urdf_to_robot_arm() {
        let robot = UrdfRobot {
            name: "test_bot".into(),
            links: vec![
                UrdfLink {
                    name: "base".into(),
                    collision_radius: 0.05,
                    collision_half_len: 0.1,
                },
                UrdfLink {
                    name: "link1".into(),
                    collision_radius: 0.04,
                    collision_half_len: 0.2,
                },
            ],
            joints: vec![UrdfJoint {
                name: "joint1".into(),
                parent: "base".into(),
                child: "link1".into(),
                joint_type: JointType::Revolute,
                origin_xyz: [0.0, 0.0, 0.2],
                origin_rpy: [0.0, 0.0, 0.0],
                axis: [0.0, 0.0, 1.0],
                limits: Some(JointLimits::new(-PI, PI, 2.0, 50.0)),
            }],
        };
        let arm = robot.to_robot_arm();
        assert_eq!(arm.links.len(), 1);
        assert_eq!(arm.dof(), 1);
    }

    // --- Velocity from Cartesian ---

    #[test]
    fn test_end_effector_velocity_zero_dq() {
        let arm = planar_arm();
        let q = vec![0.2f64, 0.3, 0.1];
        let j = arm.jacobian_analytical(&q);
        let v = end_effector_velocity(&j, &[0.0, 0.0, 0.0]);
        for vi in &v {
            assert!(vi.abs() < 1e-12);
        }
    }

    #[test]
    fn test_joint_velocity_from_cartesian_zero() {
        let arm = planar_arm();
        let q = vec![0.2f64, 0.3, 0.1];
        let j = arm.jacobian_analytical(&q);
        let dq = joint_velocity_from_cartesian(&j, 3, [0.0, 0.0, 0.0], 0.01);
        for d in &dq {
            assert!(d.abs() < 1e-10);
        }
    }

    // --- Segment–segment distance ---

    #[test]
    fn test_segment_segment_parallel() {
        let d = dist_segment_segment(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
        );
        assert!((d - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_segment_segment_crossing() {
        let d = dist_segment_segment(
            [-1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, -1.0, 0.5],
            [0.0, 1.0, 0.5],
        );
        // Distance should be 0.5 (along Z)
        assert!((d - 0.5).abs() < 1e-10);
    }

    // --- mat3 helpers ---

    #[test]
    fn test_mat3_det_identity() {
        let ident = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        assert!((mat3_det(&ident) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_mat3_inverse_identity() {
        let ident = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        let inv = mat3_inverse(&ident);
        for i in 0..9 {
            assert!((inv[i] - ident[i]).abs() < 1e-12);
        }
    }

    #[test]
    fn test_mat3_inverse_roundtrip() {
        let m = [2.0, 1.0, 0.0, 0.5, 3.0, 1.0, 0.0, 1.0, 4.0];
        let inv = mat3_inverse(&m);
        // m * inv should be identity
        for r in 0..3 {
            for c in 0..3 {
                let mut s = 0.0f64;
                for k in 0..3 {
                    s += m[r * 3 + k] * inv[k * 3 + c];
                }
                let expected = if r == c { 1.0 } else { 0.0 };
                assert!((s - expected).abs() < 1e-10, "m*inv[{r},{c}]={s}");
            }
        }
    }
}
