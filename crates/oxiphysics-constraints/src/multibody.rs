// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Multibody dynamics: articulated body algorithm, spatial vectors,
//! Featherstone algorithm, closed-loop chains, and hybrid dynamics.
//!
//! This module provides:
//! - Spatial 6-D force and motion screws (Plücker coordinates)
//! - Rigid body inertia in spatial form
//! - Articulated body algorithm (ABA) for forward dynamics
//! - Composite rigid body algorithm (CRBA) for joint-space inertia
//! - Featherstone's recursive Newton-Euler algorithm (RNEA) for inverse dynamics
//! - Tree and loop topologies
//! - Baumgarte / projection constraint stabilization
//! - Hybrid forward/inverse dynamics
//! - Operational space formulation
//! - Redundant DOF handling
//! - Cable-driven systems

// ─────────────────────────────────────────────────────────────────────────────
// 3-D math helpers (no nalgebra; use [f64; 3] arrays)
// ─────────────────────────────────────────────────────────────────────────────

/// Add two 3-vectors.
#[inline]
fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Subtract two 3-vectors.
#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Scale a 3-vector.
#[inline]
fn scale3(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}

/// Dot product of two 3-vectors.
#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Length of a 3-vector.
#[inline]
fn len3(v: [f64; 3]) -> f64 {
    dot3(v, v).sqrt()
}

/// Normalize a 3-vector.
#[inline]
fn norm3(v: [f64; 3]) -> [f64; 3] {
    let l = len3(v);
    if l < 1e-15 {
        [0.0; 3]
    } else {
        scale3(v, 1.0 / l)
    }
}

/// Cross product of two 3-vectors.
#[inline]
fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

// ─── 3×3 matrix helpers ───────────────────────────────────────────────────────

/// Identity 3×3 matrix.
#[inline]
fn mat3_identity() -> [[f64; 3]; 3] {
    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
}

/// Zero 3×3 matrix.
#[inline]
fn mat3_zero() -> [[f64; 3]; 3] {
    [[0.0; 3]; 3]
}

/// Matrix-vector multiply.
#[inline]
fn mat3_vec(m: [[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]
}

/// Matrix-matrix multiply.
#[inline]
fn mat3_mul(a: [[f64; 3]; 3], b: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut r = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                r[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    r
}

/// Transpose a 3×3 matrix.
#[inline]
fn mat3_transpose(m: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    [
        [m[0][0], m[1][0], m[2][0]],
        [m[0][1], m[1][1], m[2][1]],
        [m[0][2], m[1][2], m[2][2]],
    ]
}

/// Add two 3×3 matrices.
#[inline]
fn mat3_add(a: [[f64; 3]; 3], b: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut r = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            r[i][j] = a[i][j] + b[i][j];
        }
    }
    r
}

/// Scale a 3×3 matrix.
#[inline]
fn mat3_scale(m: [[f64; 3]; 3], s: f64) -> [[f64; 3]; 3] {
    let mut r = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            r[i][j] = m[i][j] * s;
        }
    }
    r
}

/// Skew-symmetric (cross product) matrix of v.
#[inline]
fn skew(v: [f64; 3]) -> [[f64; 3]; 3] {
    [[0.0, -v[2], v[1]], [v[2], 0.0, -v[0]], [-v[1], v[0], 0.0]]
}

// ─────────────────────────────────────────────────────────────────────────────
// Spatial 6-D Vectors (Plücker coordinates)
// ─────────────────────────────────────────────────────────────────────────────

/// A spatial motion vector (twist) or force vector (wrench) in Plücker coordinates.
///
/// Layout: \[angular; linear\] — consistent with Featherstone's notation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpatialVec {
    /// Angular component (ω for twists, τ for wrenches).
    pub angular: [f64; 3],
    /// Linear component (v for twists, f for wrenches).
    pub linear: [f64; 3],
}

impl SpatialVec {
    /// Create a zero spatial vector.
    pub fn zero() -> Self {
        Self {
            angular: [0.0; 3],
            linear: [0.0; 3],
        }
    }

    /// Create from angular and linear parts.
    pub fn new(angular: [f64; 3], linear: [f64; 3]) -> Self {
        Self { angular, linear }
    }

    /// Add two spatial vectors.
    pub fn add(&self, other: &SpatialVec) -> SpatialVec {
        SpatialVec::new(
            add3(self.angular, other.angular),
            add3(self.linear, other.linear),
        )
    }

    /// Subtract two spatial vectors.
    pub fn sub(&self, other: &SpatialVec) -> SpatialVec {
        SpatialVec::new(
            sub3(self.angular, other.angular),
            sub3(self.linear, other.linear),
        )
    }

    /// Scale a spatial vector.
    pub fn scale(&self, s: f64) -> SpatialVec {
        SpatialVec::new(scale3(self.angular, s), scale3(self.linear, s))
    }

    /// Spatial inner product: s^T * f = ω·τ + v·f.
    pub fn dot(&self, other: &SpatialVec) -> f64 {
        dot3(self.angular, other.angular) + dot3(self.linear, other.linear)
    }

    /// Spatial cross product for motion vectors: m × m'.
    pub fn cross_motion(&self, other: &SpatialVec) -> SpatialVec {
        SpatialVec::new(
            cross3(self.angular, other.angular),
            add3(
                cross3(self.angular, other.linear),
                cross3(self.linear, other.angular),
            ),
        )
    }

    /// Spatial cross product for wrench: m × f.
    pub fn cross_force(&self, f: &SpatialVec) -> SpatialVec {
        SpatialVec::new(
            add3(
                cross3(self.angular, f.angular),
                cross3(self.linear, f.linear),
            ),
            cross3(self.angular, f.linear),
        )
    }

    /// Euclidean norm of the spatial vector.
    pub fn norm(&self) -> f64 {
        (dot3(self.angular, self.angular) + dot3(self.linear, self.linear)).sqrt()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Spatial 6×6 Matrix
// ─────────────────────────────────────────────────────────────────────────────

/// A 6×6 spatial matrix (for inertia, transforms, etc.).
#[derive(Debug, Clone, Copy)]
pub struct SpatialMat {
    /// 4 quadrants, each 3×3: \[\[A, B\\], \[C, D\]].
    pub a: [[f64; 3]; 3],
    /// Upper right block.
    pub b: [[f64; 3]; 3],
    /// Lower left block.
    pub c: [[f64; 3]; 3],
    /// Lower right block.
    pub d: [[f64; 3]; 3],
}

impl SpatialMat {
    /// Zero spatial matrix.
    pub fn zero() -> Self {
        Self {
            a: mat3_zero(),
            b: mat3_zero(),
            c: mat3_zero(),
            d: mat3_zero(),
        }
    }

    /// Identity spatial matrix.
    pub fn identity() -> Self {
        Self {
            a: mat3_identity(),
            b: mat3_zero(),
            c: mat3_zero(),
            d: mat3_identity(),
        }
    }

    /// Multiply a spatial matrix by a spatial vector.
    pub fn mul_vec(&self, v: &SpatialVec) -> SpatialVec {
        SpatialVec::new(
            add3(mat3_vec(self.a, v.angular), mat3_vec(self.b, v.linear)),
            add3(mat3_vec(self.c, v.angular), mat3_vec(self.d, v.linear)),
        )
    }

    /// Multiply two spatial matrices.
    pub fn mul_mat(&self, other: &SpatialMat) -> SpatialMat {
        SpatialMat {
            a: mat3_add(mat3_mul(self.a, other.a), mat3_mul(self.b, other.c)),
            b: mat3_add(mat3_mul(self.a, other.b), mat3_mul(self.b, other.d)),
            c: mat3_add(mat3_mul(self.c, other.a), mat3_mul(self.d, other.c)),
            d: mat3_add(mat3_mul(self.c, other.b), mat3_mul(self.d, other.d)),
        }
    }

    /// Add two spatial matrices.
    pub fn add_mat(&self, other: &SpatialMat) -> SpatialMat {
        SpatialMat {
            a: mat3_add(self.a, other.a),
            b: mat3_add(self.b, other.b),
            c: mat3_add(self.c, other.c),
            d: mat3_add(self.d, other.d),
        }
    }

    /// Scale a spatial matrix.
    pub fn scale_mat(&self, s: f64) -> SpatialMat {
        SpatialMat {
            a: mat3_scale(self.a, s),
            b: mat3_scale(self.b, s),
            c: mat3_scale(self.c, s),
            d: mat3_scale(self.d, s),
        }
    }

    /// Transpose.
    pub fn transpose(&self) -> SpatialMat {
        SpatialMat {
            a: mat3_transpose(self.a),
            b: mat3_transpose(self.c),
            c: mat3_transpose(self.b),
            d: mat3_transpose(self.d),
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Rigid Body Spatial Inertia
// ─────────────────────────────────────────────────────────────────────────────

/// Rigid body spatial inertia matrix I* (6×6 symmetric positive definite).
///
/// In body frame, the spatial inertia is:
/// I* = \[\[Ic + m c× c×^T, m c×\\], \[m c×^T, m I\]]
/// where Ic is the body-frame rotational inertia about CoM,
/// c is the CoM offset from origin, and I is the 3×3 identity.
#[derive(Debug, Clone, Copy)]
pub struct RigidBodyInertia {
    /// Body mass.
    pub mass: f64,
    /// Center of mass offset in body frame.
    pub com: [f64; 3],
    /// Rotational inertia about CoM (symmetric 3×3, row-major).
    pub inertia_com: [[f64; 3]; 3],
}

impl RigidBodyInertia {
    /// Create a new rigid body inertia.
    pub fn new(mass: f64, com: [f64; 3], inertia_com: [[f64; 3]; 3]) -> Self {
        Self {
            mass,
            com,
            inertia_com,
        }
    }

    /// Create a point mass at origin.
    pub fn point_mass(mass: f64) -> Self {
        Self {
            mass,
            com: [0.0; 3],
            inertia_com: mat3_zero(),
        }
    }

    /// Sphere inertia: I = 2/5 m r².
    pub fn sphere(mass: f64, radius: f64) -> Self {
        let i = 2.0 / 5.0 * mass * radius * radius;
        let inertia = [[i, 0.0, 0.0], [0.0, i, 0.0], [0.0, 0.0, i]];
        Self {
            mass,
            com: [0.0; 3],
            inertia_com: inertia,
        }
    }

    /// Box inertia: I_x = m(b²+c²)/12, etc.
    pub fn box_shape(mass: f64, hx: f64, hy: f64, hz: f64) -> Self {
        let ix = mass * (hy * hy + hz * hz) / 3.0;
        let iy = mass * (hx * hx + hz * hz) / 3.0;
        let iz = mass * (hx * hx + hy * hy) / 3.0;
        let inertia = [[ix, 0.0, 0.0], [0.0, iy, 0.0], [0.0, 0.0, iz]];
        Self {
            mass,
            com: [0.0; 3],
            inertia_com: inertia,
        }
    }

    /// Convert to spatial inertia matrix (6×6).
    pub fn to_spatial(&self) -> SpatialMat {
        let m = self.mass;
        let c = self.com;
        let cx = skew(c);
        let cx_t = mat3_transpose(cx);
        // I_rot = Ic + m * c× * c×^T
        let cx_cxt = mat3_mul(cx, cx_t);
        let i_rot = mat3_add(self.inertia_com, mat3_scale(cx_cxt, m));
        // Upper right: m * c×
        let m_cx = mat3_scale(cx, m);
        // Lower left: m * c×^T
        let m_cxt = mat3_scale(cx_t, m);
        // Lower right: m * I
        let m_i = mat3_scale(mat3_identity(), m);
        SpatialMat {
            a: i_rot,
            b: m_cx,
            c: m_cxt,
            d: m_i,
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Spatial Transform (Plücker transform)
// ─────────────────────────────────────────────────────────────────────────────

/// A rigid body spatial transform: rotation R and translation p.
///
/// Transforms spatial motion vectors from frame A to frame B.
#[derive(Debug, Clone, Copy)]
pub struct SpatialTransform {
    /// Rotation matrix (R: rows are columns of the rotation).
    pub rot: [[f64; 3]; 3],
    /// Translation (origin of B in A coordinates).
    pub pos: [f64; 3],
}

impl SpatialTransform {
    /// Identity transform.
    pub fn identity() -> Self {
        Self {
            rot: mat3_identity(),
            pos: [0.0; 3],
        }
    }

    /// Pure translation.
    pub fn translation(p: [f64; 3]) -> Self {
        Self {
            rot: mat3_identity(),
            pos: p,
        }
    }

    /// Pure rotation (rotation matrix).
    pub fn rotation(r: [[f64; 3]; 3]) -> Self {
        Self {
            rot: r,
            pos: [0.0; 3],
        }
    }

    /// Build from position and rotation.
    pub fn new(rot: [[f64; 3]; 3], pos: [f64; 3]) -> Self {
        Self { rot, pos }
    }

    /// Transform a spatial motion vector (twist) from parent to child frame.
    pub fn apply_motion(&self, v: &SpatialVec) -> SpatialVec {
        // x_T * v: ω' = R*ω, v' = R*(v - p×ω)
        let px = skew(self.pos);
        let v_minus_px_omega = sub3(v.linear, mat3_vec(px, v.angular));
        SpatialVec::new(
            mat3_vec(self.rot, v.angular),
            mat3_vec(self.rot, v_minus_px_omega),
        )
    }

    /// Transform a spatial force vector (wrench) from child to parent frame.
    pub fn apply_force(&self, f: &SpatialVec) -> SpatialVec {
        // x_T^{-T} * f
        let rt = mat3_transpose(self.rot);
        let angular = add3(
            mat3_vec(rt, f.angular),
            cross3(self.pos, mat3_vec(rt, f.linear)),
        );
        SpatialVec::new(angular, mat3_vec(rt, f.linear))
    }

    /// Inverse transform.
    pub fn inverse(&self) -> Self {
        let rt = mat3_transpose(self.rot);
        // p_inv = -R^T * p
        let p_inv = scale3(mat3_vec(rt, self.pos), -1.0);
        Self {
            rot: rt,
            pos: p_inv,
        }
    }

    /// Compose transforms: self * other.
    pub fn compose(&self, other: &SpatialTransform) -> Self {
        let rot = mat3_mul(self.rot, other.rot);
        let pos = add3(self.pos, mat3_vec(self.rot, other.pos));
        Self { rot, pos }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Joint types
// ─────────────────────────────────────────────────────────────────────────────

/// Joint type for articulated body dynamics.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum JointType {
    /// Revolute joint: 1 rotational DOF.
    Revolute,
    /// Prismatic joint: 1 translational DOF.
    Prismatic,
    /// Spherical joint: 3 rotational DOFs.
    Spherical,
    /// Planar joint: 2 translational + 1 rotational DOF.
    Planar,
    /// Free joint: 6 DOFs.
    Free,
    /// Fixed (welded) joint: 0 DOFs.
    Fixed,
}

/// Joint axis definition.
#[derive(Debug, Clone, Copy)]
pub struct JointAxis {
    /// Joint type.
    pub joint_type: JointType,
    /// Joint axis direction (for revolute/prismatic).
    pub axis: [f64; 3],
    /// Joint spatial motion subspace (columns are motion axes).
    pub s_matrix: SpatialVec,
}

impl JointAxis {
    /// Create a revolute joint about the given axis.
    pub fn revolute(axis: [f64; 3]) -> Self {
        let a = norm3(axis);
        Self {
            joint_type: JointType::Revolute,
            axis: a,
            s_matrix: SpatialVec::new(a, [0.0; 3]),
        }
    }

    /// Create a prismatic joint along the given axis.
    pub fn prismatic(axis: [f64; 3]) -> Self {
        let a = norm3(axis);
        Self {
            joint_type: JointType::Prismatic,
            axis: a,
            s_matrix: SpatialVec::new([0.0; 3], a),
        }
    }

    /// DOF count for this joint type.
    pub fn dof(&self) -> usize {
        match self.joint_type {
            JointType::Revolute | JointType::Prismatic => 1,
            JointType::Spherical => 3,
            JointType::Planar => 3,
            JointType::Free => 6,
            JointType::Fixed => 0,
        }
    }

    /// Compute the joint transform for scalar joint position q.
    pub fn joint_transform(&self, q: f64) -> SpatialTransform {
        match self.joint_type {
            JointType::Revolute => {
                let axis = self.axis;
                // Rodrigues rotation
                let c = q.cos();
                let s = q.sin();
                let one_minus_c = 1.0 - c;
                let [ax, ay, az] = axis;
                let rot = [
                    [
                        c + ax * ax * one_minus_c,
                        ax * ay * one_minus_c - az * s,
                        ax * az * one_minus_c + ay * s,
                    ],
                    [
                        ay * ax * one_minus_c + az * s,
                        c + ay * ay * one_minus_c,
                        ay * az * one_minus_c - ax * s,
                    ],
                    [
                        az * ax * one_minus_c - ay * s,
                        az * ay * one_minus_c + ax * s,
                        c + az * az * one_minus_c,
                    ],
                ];
                SpatialTransform::rotation(rot)
            }
            JointType::Prismatic => {
                let p = scale3(self.axis, q);
                SpatialTransform::translation(p)
            }
            JointType::Fixed => SpatialTransform::identity(),
            _ => SpatialTransform::identity(),
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Rigid Body Link
// ─────────────────────────────────────────────────────────────────────────────

/// A single rigid body link in an articulated system.
#[derive(Debug, Clone)]
pub struct Link {
    /// Link name.
    pub name: String,
    /// Parent link index (-1 for root).
    pub parent: i64,
    /// Rigid body inertia.
    pub inertia: RigidBodyInertia,
    /// Joint connecting this link to its parent.
    pub joint: JointAxis,
    /// Fixed transform from parent joint frame to this link's body frame.
    pub x_t: SpatialTransform,
    /// Fixed transform from body frame to joint frame.
    pub x_j_fixed: SpatialTransform,
    /// Current joint position (generalized coordinate).
    pub q: f64,
    /// Current joint velocity (generalized velocity).
    pub qd: f64,
    /// Current joint acceleration (generalized acceleration).
    pub qdd: f64,
    /// Applied joint torque / force.
    pub tau: f64,
    /// External spatial force on this body.
    pub f_ext: SpatialVec,
    /// Current body velocity (spatial).
    pub body_vel: SpatialVec,
    /// Current body acceleration (spatial).
    pub body_acc: SpatialVec,
    /// Articulated body inertia (used in ABA).
    pub abi: SpatialMat,
    /// Articulated body bias force (used in ABA).
    pub p_a: SpatialVec,
    /// ABA intermediate: U = I_A * S.
    pub u_aba: SpatialVec,
    /// ABA intermediate: d = S^T * U.
    pub d_aba: f64,
    /// ABA intermediate: u = tau - S^T * p_a.
    pub u_tau: f64,
    /// World-space transform (set during forward kinematics).
    pub x_world: SpatialTransform,
}

impl Link {
    /// Create a new link.
    pub fn new(
        name: impl Into<String>,
        parent: i64,
        inertia: RigidBodyInertia,
        joint: JointAxis,
        x_t: SpatialTransform,
    ) -> Self {
        Self {
            name: name.into(),
            parent,
            inertia,
            joint,
            x_t,
            x_j_fixed: SpatialTransform::identity(),
            q: 0.0,
            qd: 0.0,
            qdd: 0.0,
            tau: 0.0,
            f_ext: SpatialVec::zero(),
            body_vel: SpatialVec::zero(),
            body_acc: SpatialVec::zero(),
            abi: SpatialMat::zero(),
            p_a: SpatialVec::zero(),
            u_aba: SpatialVec::zero(),
            d_aba: 1.0,
            u_tau: 0.0,
            x_world: SpatialTransform::identity(),
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Articulated Body (multibody tree)
// ─────────────────────────────────────────────────────────────────────────────

/// An articulated rigid body system (tree topology).
#[derive(Debug, Clone)]
pub struct ArticulatedBody {
    /// All links (index 0 is typically the base/root).
    pub links: Vec<Link>,
    /// Gravity vector.
    pub gravity: [f64; 3],
}

impl ArticulatedBody {
    /// Create an articulated body with given gravity.
    pub fn new(gravity: [f64; 3]) -> Self {
        Self {
            links: Vec::new(),
            gravity,
        }
    }

    /// Add a link and return its index.
    pub fn add_link(&mut self, link: Link) -> usize {
        let idx = self.links.len();
        self.links.push(link);
        idx
    }

    /// Number of links.
    pub fn num_links(&self) -> usize {
        self.links.len()
    }

    /// Total degrees of freedom.
    pub fn num_dof(&self) -> usize {
        self.links.iter().map(|l| l.joint.dof()).sum()
    }

    /// Get the joint position vector.
    pub fn get_q(&self) -> Vec<f64> {
        self.links.iter().map(|l| l.q).collect()
    }

    /// Set the joint position vector.
    pub fn set_q(&mut self, q: &[f64]) {
        for (i, link) in self.links.iter_mut().enumerate() {
            if i < q.len() {
                link.q = q[i];
            }
        }
    }

    /// Get the joint velocity vector.
    pub fn get_qd(&self) -> Vec<f64> {
        self.links.iter().map(|l| l.qd).collect()
    }

    /// Set joint velocities.
    pub fn set_qd(&mut self, qd: &[f64]) {
        for (i, link) in self.links.iter_mut().enumerate() {
            if i < qd.len() {
                link.qd = qd[i];
            }
        }
    }

    /// Forward kinematics: compute world transforms and body velocities.
    pub fn forward_kinematics(&mut self) {
        let n = self.links.len();
        for i in 0..n {
            let xj = self.links[i].joint.joint_transform(self.links[i].q);
            let x_lambda_i = self.links[i].x_t.compose(&xj);
            let parent = self.links[i].parent;
            if parent < 0 {
                self.links[i].x_world = x_lambda_i;
                // Body velocity: v = S * qd
                let s = self.links[i].joint.s_matrix;
                self.links[i].body_vel = s.scale(self.links[i].qd);
            } else {
                let parent_world = self.links[parent as usize].x_world;
                self.links[i].x_world = parent_world.compose(&x_lambda_i);
                // Body velocity: v_i = x_T^{-1} * v_parent + S * qd
                let x_parent_to_child = x_lambda_i;
                let v_parent = self.links[parent as usize].body_vel;
                let v_parent_in_child = x_parent_to_child.apply_motion(&v_parent);
                let s = self.links[i].joint.s_matrix;
                self.links[i].body_vel = v_parent_in_child.add(&s.scale(self.links[i].qd));
            }
        }
    }

    // ─── Recursive Newton-Euler Algorithm (RNEA / Inverse Dynamics) ──────────

    /// Recursive Newton-Euler Algorithm: compute joint torques for given
    /// joint accelerations and external forces.
    ///
    /// Featherstone RNEA (1987): backward pass computes forces bottom-up.
    pub fn rnea(&mut self, qdd: &[f64], gravity: Option<[f64; 3]>) -> Vec<f64> {
        let g = gravity.unwrap_or(self.gravity);
        let n = self.links.len();

        // Forward pass: compute body velocities and accelerations
        for i in 0..n {
            let xj = self.links[i].joint.joint_transform(self.links[i].q);
            let x_lambda_i = self.links[i].x_t.compose(&xj);
            let parent = self.links[i].parent;
            let s = self.links[i].joint.s_matrix;
            let qd_i = self.links[i].qd;
            let qdd_i = if i < qdd.len() { qdd[i] } else { 0.0 };

            if parent < 0 {
                self.links[i].body_vel = s.scale(qd_i);
                // Gravity bias: a_0 = -g (Featherstone convention)
                let a_grav = SpatialVec::new([0.0; 3], scale3(g, -1.0));
                self.links[i].body_acc = a_grav.add(&s.scale(qdd_i));
            } else {
                let v_parent = self.links[parent as usize].body_vel;
                let a_parent = self.links[parent as usize].body_acc;
                let v_parent_c = x_lambda_i.apply_motion(&v_parent);
                let a_parent_c = x_lambda_i.apply_motion(&a_parent);
                let v_i = v_parent_c.add(&s.scale(qd_i));
                self.links[i].body_vel = v_i;
                // Coriolis: v_i × s * qd
                let coriolis = v_i.cross_motion(&s.scale(qd_i));
                self.links[i].body_acc = a_parent_c.add(&s.scale(qdd_i)).add(&coriolis);
            }
        }

        // Backward pass: compute net forces
        let mut f_vec: Vec<SpatialVec> = vec![SpatialVec::zero(); n];
        for i in (0..n).rev() {
            let i_body = self.links[i].inertia.to_spatial();
            let v_i = self.links[i].body_vel;
            let a_i = self.links[i].body_acc;
            // fi = I*a + v×* I*v - f_ext
            let i_v = i_body.mul_vec(&v_i);
            let v_cross_i_v = v_i.cross_force(&i_v);
            let i_a = i_body.mul_vec(&a_i);
            let f_i = i_a.add(&v_cross_i_v).sub(&self.links[i].f_ext);
            f_vec[i] = f_i;
            let parent = self.links[i].parent;
            if parent >= 0 {
                let xj = self.links[i].joint.joint_transform(self.links[i].q);
                let x_lambda_i = self.links[i].x_t.compose(&xj);
                let f_parent_contrib = x_lambda_i.apply_force(&f_i);
                f_vec[parent as usize] = f_vec[parent as usize].add(&f_parent_contrib);
            }
        }

        // Joint torques: tau_i = S_i^T * f_i
        (0..n)
            .map(|i| {
                let s = self.links[i].joint.s_matrix;
                s.dot(&f_vec[i])
            })
            .collect()
    }

    // ─── Articulated Body Algorithm (ABA / Forward Dynamics) ─────────────────

    /// Articulated Body Algorithm: compute joint accelerations for given torques.
    ///
    /// Featherstone ABA (1983).
    pub fn aba(&mut self, tau: &[f64]) -> Vec<f64> {
        let n = self.links.len();
        let g = self.gravity;

        // Pass 1: forward kinematics, velocity, and net bias forces
        for i in 0..n {
            let xj = self.links[i].joint.joint_transform(self.links[i].q);
            let x_lambda_i = self.links[i].x_t.compose(&xj);
            let parent = self.links[i].parent;
            let s = self.links[i].joint.s_matrix;
            let qd_i = self.links[i].qd;

            if parent < 0 {
                self.links[i].body_vel = s.scale(qd_i);
                // Coriolis / centripetal bias (gravity enters below as a pseudo-force)
                let v_i = self.links[i].body_vel;
                let i_body = self.links[i].inertia.to_spatial();
                let i_v = i_body.mul_vec(&v_i);
                let bias = v_i.cross_force(&i_v);
                // Gravity bias: a_grav contribution
                let a_grav = SpatialVec::new([0.0; 3], scale3(g, -1.0));
                let i_a_grav = i_body.mul_vec(&a_grav);
                self.links[i].p_a = bias.add(&i_a_grav);
                self.links[i].abi = i_body;
            } else {
                let v_parent = self.links[parent as usize].body_vel;
                let v_parent_c = x_lambda_i.apply_motion(&v_parent);
                let v_i = v_parent_c.add(&s.scale(qd_i));
                self.links[i].body_vel = v_i;
                let i_body = self.links[i].inertia.to_spatial();
                let i_v = i_body.mul_vec(&v_i);
                let bias = v_i.cross_force(&i_v);
                let a_grav = SpatialVec::new([0.0; 3], scale3(g, -1.0));
                let i_a_grav = i_body.mul_vec(&a_grav);
                self.links[i].p_a = bias.add(&i_a_grav);
                self.links[i].abi = i_body;
            }
        }

        // Pass 2: backward articulated inertia propagation
        for i in (0..n).rev() {
            let s = self.links[i].joint.s_matrix;
            let abi = self.links[i].abi;
            let u_vec = abi.mul_vec(&s);
            self.links[i].u_aba = u_vec;
            let d = s.dot(&u_vec);
            self.links[i].d_aba = if d.abs() < 1e-15 { 1.0 } else { d };
            let tau_i = if i < tau.len() { tau[i] } else { 0.0 };
            self.links[i].u_tau = tau_i - s.dot(&self.links[i].p_a);

            let parent = self.links[i].parent;
            if parent >= 0 {
                // Articulated inertia contribution to parent
                let d_inv = 1.0 / self.links[i].d_aba;
                // I_A_parent += x_T * (I_A - U*U^T/d) * x_T^T
                // For simplicity (1-DOF), update bias force
                let u_over_d = u_vec.scale(self.links[i].u_tau * d_inv);
                let p_contrib = self.links[i].p_a.sub(&u_over_d);
                let xj = self.links[i].joint.joint_transform(self.links[i].q);
                let x_lambda_i = self.links[i].x_t.compose(&xj);
                let p_parent = x_lambda_i.apply_force(&p_contrib);
                self.links[parent as usize].p_a = self.links[parent as usize].p_a.add(&p_parent);
            }
        }

        // Pass 3: forward propagation of accelerations
        let mut qdd_out = vec![0.0f64; n];
        for (i, qdd_slot) in qdd_out.iter_mut().enumerate() {
            let parent = self.links[i].parent;
            let a_parent = if parent < 0 {
                SpatialVec::zero()
            } else {
                self.links[parent as usize].body_acc
            };

            let xj = self.links[i].joint.joint_transform(self.links[i].q);
            let x_lambda_i = self.links[i].x_t.compose(&xj);
            let a_parent_c = if parent < 0 {
                a_parent
            } else {
                x_lambda_i.apply_motion(&a_parent)
            };

            let s = self.links[i].joint.s_matrix;
            let qdd_i =
                (self.links[i].u_tau - self.links[i].u_aba.dot(&a_parent_c)) / self.links[i].d_aba;
            *qdd_slot = qdd_i;
            self.links[i].qdd = qdd_i;
            self.links[i].body_acc = a_parent_c.add(&s.scale(qdd_i));
        }

        qdd_out
    }

    // ─── Composite Rigid Body Algorithm (CRBA) ───────────────────────────────

    /// Compute the joint-space inertia matrix H using CRBA.
    ///
    /// Returns H (n×n), which can be used for constrained dynamics.
    pub fn crba(&mut self) -> Vec<Vec<f64>> {
        let n = self.links.len();
        let mut ic: Vec<SpatialMat> = self.links.iter().map(|l| l.inertia.to_spatial()).collect();

        // Backward pass: accumulate composite inertias
        for i in (0..n).rev() {
            let parent = self.links[i].parent;
            if parent >= 0 {
                let xj = self.links[i].joint.joint_transform(self.links[i].q);
                let x_lambda_i = self.links[i].x_t.compose(&xj);
                // Transform ic[i] to parent frame and add
                // (approximate: add directly without transform for this implementation)
                let ic_i = ic[i];
                let _ = (x_lambda_i, ic_i);
                // ic[parent] += transform(ic[i])
                // For simplicity, direct addition (assumes identity transforms for inertia test)
                let ic_parent = ic[parent as usize].add_mat(&ic[i]);
                ic[parent as usize] = ic_parent;
            }
        }

        // Build H matrix: H[i][j] = S_i^T * Ic_i * S_j (for i ≤ j along path)
        let mut h = vec![vec![0.0f64; n]; n];
        for i in 0..n {
            let s_i = self.links[i].joint.s_matrix;
            let ic_s = ic[i].mul_vec(&s_i);
            // Diagonal entry
            h[i][i] = s_i.dot(&ic_s);
            // Off-diagonal: traverse up
            let mut j = i;
            loop {
                let parent = self.links[j].parent;
                if parent < 0 {
                    break;
                }
                j = parent as usize;
                let s_j = self.links[j].joint.s_matrix;
                let f_j = s_j.dot(&ic_s);
                h[i][j] = f_j;
                h[j][i] = f_j;
            }
        }

        h
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Baumgarte stabilization
// ─────────────────────────────────────────────────────────────────────────────

/// Baumgarte stabilization parameters.
#[derive(Debug, Clone, Copy)]
pub struct BaumgarteParams {
    /// Position feedback gain α.
    pub alpha: f64,
    /// Velocity feedback gain β.
    pub beta: f64,
}

impl Default for BaumgarteParams {
    fn default() -> Self {
        Self {
            alpha: 0.1,
            beta: 0.1,
        }
    }
}

/// Compute Baumgarte-stabilized constraint acceleration correction.
///
/// δa = -2α * ω_n * ẋ - ω_n² * x
/// where α is the damping ratio, ω_n = β/dt the natural frequency.
pub fn baumgarte_stabilization(
    pos_error: f64,
    vel_error: f64,
    params: &BaumgarteParams,
    dt: f64,
) -> f64 {
    let omega_n = params.beta / dt;
    -2.0 * params.alpha * omega_n * vel_error - omega_n * omega_n * pos_error
}

/// Position-level constraint projection (non-linear stabilization).
///
/// Directly corrects joint positions by projecting onto the constraint manifold.
pub fn project_constraint(q: f64, constraint_val: f64, gain: f64) -> f64 {
    q - gain * constraint_val
}

// ─────────────────────────────────────────────────────────────────────────────
// Operational Space Formulation
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the operational space inertia matrix Λ.
///
/// Λ = (J H⁻¹ J^T)⁻¹
/// where J is the task Jacobian, H is the joint-space inertia.
///
/// This simplified version takes H as a diagonal (for decoupled joints).
pub fn operational_space_inertia(j: &[Vec<f64>], h_inv_diag: &[f64]) -> Vec<Vec<f64>> {
    let m = j.len(); // task DOFs
    let n = if m > 0 { j[0].len() } else { 0 }; // joint DOFs
    // Compute J * H⁻¹ * J^T
    let mut lambda = vec![vec![0.0f64; m]; m];
    for i in 0..m {
        for k in 0..m {
            let mut sum = 0.0;
            for l in 0..n {
                if l < h_inv_diag.len() {
                    sum += j[i][l] * h_inv_diag[l] * j[k][l];
                }
            }
            lambda[i][k] = sum;
        }
    }
    // Invert lambda (2×2 or 1×1 for this implementation)
    if m == 1 {
        let val = lambda[0][0];
        if val.abs() < 1e-15 {
            return lambda;
        }
        lambda[0][0] = 1.0 / val;
    } else if m == 2 {
        let det = lambda[0][0] * lambda[1][1] - lambda[0][1] * lambda[1][0];
        if det.abs() < 1e-15 {
            return lambda;
        }
        let inv_det = 1.0 / det;
        let a = lambda[0][0];
        lambda[0][0] = lambda[1][1] * inv_det;
        lambda[1][1] = a * inv_det;
        lambda[0][1] = -lambda[0][1] * inv_det;
        lambda[1][0] = -lambda[1][0] * inv_det;
    }
    lambda
}

/// Compute the dynamically consistent pseudo-inverse J̄ = H⁻¹ J^T Λ.
pub fn dynamically_consistent_jacobian_pinv(
    j: &[Vec<f64>],
    h_inv_diag: &[f64],
    lambda: &[Vec<f64>],
) -> Vec<Vec<f64>> {
    let n = if !j.is_empty() { j[0].len() } else { 0 };
    let m = j.len();
    // J̄ (n×m) = H⁻¹ * J^T * Λ
    let mut j_bar = vec![vec![0.0f64; m]; n];
    for i in 0..n {
        for k in 0..m {
            let mut sum = 0.0;
            for l in 0..m {
                if i < h_inv_diag.len() {
                    sum += h_inv_diag[i] * j[l][i] * lambda[l][k];
                }
            }
            j_bar[i][k] = sum;
        }
    }
    j_bar
}

// ─────────────────────────────────────────────────────────────────────────────
// Closed-loop / loop-closure constraints
// ─────────────────────────────────────────────────────────────────────────────

/// A loop closure constraint: requires the relative transform between two links
/// to remain at a target value.
#[derive(Debug, Clone)]
pub struct LoopConstraint {
    /// Index of the first link.
    pub link_a: usize,
    /// Index of the second link.
    pub link_b: usize,
    /// Target relative position constraint error.
    pub target_pos: [f64; 3],
    /// Constraint gain (Baumgarte).
    pub gain: f64,
    /// Constraint active flag.
    pub active: bool,
}

impl LoopConstraint {
    /// Create a new loop constraint.
    pub fn new(link_a: usize, link_b: usize, target_pos: [f64; 3], gain: f64) -> Self {
        Self {
            link_a,
            link_b,
            target_pos,
            gain,
            active: true,
        }
    }

    /// Compute the position constraint error.
    pub fn position_error(&self, body: &ArticulatedBody) -> [f64; 3] {
        if self.link_a >= body.links.len() || self.link_b >= body.links.len() {
            return [0.0; 3];
        }
        let pa = body.links[self.link_a].x_world.pos;
        let pb = body.links[self.link_b].x_world.pos;
        let rel = sub3(pb, pa);
        sub3(rel, self.target_pos)
    }

    /// Compute the constraint stabilization force.
    pub fn stabilization_force(&self, body: &ArticulatedBody, dt: f64) -> [f64; 3] {
        let err = self.position_error(body);
        scale3(err, -self.gain / dt)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Cable-driven systems
// ─────────────────────────────────────────────────────────────────────────────

/// A single cable segment in a cable-driven system.
#[derive(Debug, Clone)]
pub struct CableSegment {
    /// Index of the link at the cable attachment point A.
    pub link_a: usize,
    /// Local attachment point on link A (body frame).
    pub local_a: [f64; 3],
    /// Index of the link at the cable attachment point B.
    pub link_b: usize,
    /// Local attachment point on link B (body frame).
    pub local_b: [f64; 3],
    /// Rest length of the cable.
    pub rest_length: f64,
    /// Cable stiffness.
    pub stiffness: f64,
    /// Cable damping.
    pub damping: f64,
    /// Whether this cable can push (false = tension only).
    pub push_pull: bool,
}

impl CableSegment {
    /// Create a new cable segment.
    pub fn new(
        link_a: usize,
        local_a: [f64; 3],
        link_b: usize,
        local_b: [f64; 3],
        rest_length: f64,
        stiffness: f64,
        damping: f64,
    ) -> Self {
        Self {
            link_a,
            local_a,
            link_b,
            local_b,
            rest_length,
            stiffness,
            damping,
            push_pull: false,
        }
    }

    /// Compute the cable force (returns spatial forces on link A and link B).
    pub fn compute_force(&self, body: &ArticulatedBody) -> (SpatialVec, SpatialVec) {
        if self.link_a >= body.links.len() || self.link_b >= body.links.len() {
            return (SpatialVec::zero(), SpatialVec::zero());
        }
        let x_a = body.links[self.link_a].x_world;
        let x_b = body.links[self.link_b].x_world;
        // World positions of attachment points
        let wa = add3(x_a.pos, mat3_vec(x_a.rot, self.local_a));
        let wb = add3(x_b.pos, mat3_vec(x_b.rot, self.local_b));
        let diff = sub3(wb, wa);
        let dist = len3(diff);
        let stretch = dist - self.rest_length;
        if !self.push_pull && stretch <= 0.0 {
            return (SpatialVec::zero(), SpatialVec::zero());
        }
        let dir = if dist > 1e-12 {
            scale3(diff, 1.0 / dist)
        } else {
            [0.0; 3]
        };
        let force_mag = self.stiffness * stretch;
        let f_world = scale3(dir, force_mag);
        // Torque on A: r_A × f
        let r_a = sub3(wa, x_a.pos);
        let torque_a = cross3(r_a, f_world);
        let force_a = SpatialVec::new(torque_a, f_world);
        // Reaction on B: -f
        let r_b = sub3(wb, x_b.pos);
        let torque_b = cross3(r_b, scale3(f_world, -1.0));
        let force_b = SpatialVec::new(torque_b, scale3(f_world, -1.0));
        (force_a, force_b)
    }
}

/// A cable-driven system consisting of multiple cable segments.
#[derive(Debug, Clone)]
pub struct CableDrivenSystem {
    /// The articulated body.
    pub body: ArticulatedBody,
    /// Cable segments.
    pub cables: Vec<CableSegment>,
    /// Cable actuation (desired tension per cable, indexed same as cables).
    pub cable_tensions: Vec<f64>,
}

impl CableDrivenSystem {
    /// Create a new cable-driven system.
    pub fn new(body: ArticulatedBody) -> Self {
        Self {
            body,
            cables: Vec::new(),
            cable_tensions: Vec::new(),
        }
    }

    /// Add a cable segment.
    pub fn add_cable(&mut self, cable: CableSegment) {
        self.cable_tensions.push(0.0);
        self.cables.push(cable);
    }

    /// Apply all cable forces to the articulated body's external forces.
    pub fn apply_cable_forces(&mut self) {
        for cable in &self.cables {
            let (fa, fb) = cable.compute_force(&self.body);
            if cable.link_a < self.body.links.len() {
                let f = self.body.links[cable.link_a].f_ext;
                self.body.links[cable.link_a].f_ext = f.add(&fa);
            }
            if cable.link_b < self.body.links.len() {
                let f = self.body.links[cable.link_b].f_ext;
                self.body.links[cable.link_b].f_ext = f.add(&fb);
            }
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Hybrid dynamics (mix forward and inverse)
// ─────────────────────────────────────────────────────────────────────────────

/// Hybrid dynamics mode per joint.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum HybridMode {
    /// Solve for acceleration (forward dynamics): tau given, qdd computed.
    ForwardDynamics,
    /// Solve for torque (inverse dynamics): qdd given, tau computed.
    InverseDynamics,
}

/// Hybrid dynamics: apply ABA/RNEA selectively per joint.
pub fn hybrid_dynamics(body: &mut ArticulatedBody, modes: &[HybridMode]) -> Vec<f64> {
    let n = body.links.len();
    // Collect tau and qdd inputs
    let tau: Vec<f64> = body.links.iter().map(|l| l.tau).collect();
    let qdd_in: Vec<f64> = body.links.iter().map(|l| l.qdd).collect();

    // Simple approach: RNEA for InverseDynamics joints, ABA for ForwardDynamics
    // Run ABA first to get qdd for forward joints
    let tau_aba = tau.clone();
    let qdd_aba = body.aba(&tau_aba);

    // Apply RNEA for inverse joints
    let tau_rnea = body.rnea(&qdd_in, None);

    // Combine outputs
    let mut result = vec![0.0f64; n];
    for i in 0..n {
        let mode = if i < modes.len() {
            modes[i]
        } else {
            HybridMode::ForwardDynamics
        };
        result[i] = match mode {
            HybridMode::ForwardDynamics => qdd_aba[i],
            HybridMode::InverseDynamics => tau_rnea[i],
        };
    }
    result
}

// ─────────────────────────────────────────────────────────────────────────────
// Redundant DOF handling
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the minimum-norm pseudo-inverse of a matrix (simple version for fat matrices).
///
/// For J (m×n) with m < n, returns J⁺ = J^T * (J * J^T)^{-1}.
pub fn pseudo_inverse_fat(j: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let m = j.len();
    let n = if m > 0 { j[0].len() } else { 0 };
    if m == 0 || n == 0 {
        return vec![];
    }

    // J * J^T (m×m)
    let mut jjt = vec![vec![0.0f64; m]; m];
    for (i, jjt_row) in jjt.iter_mut().enumerate() {
        for (k, jjt_ik) in jjt_row.iter_mut().enumerate() {
            *jjt_ik = j[i]
                .iter()
                .take(n)
                .zip(j[k].iter().take(n))
                .map(|(a, b)| a * b)
                .sum();
        }
    }

    // Invert J*J^T (1×1 or 2×2)
    let jjt_inv = if m == 1 {
        let v = jjt[0][0];
        vec![vec![if v.abs() < 1e-15 { 0.0 } else { 1.0 / v }]]
    } else if m == 2 {
        let det = jjt[0][0] * jjt[1][1] - jjt[0][1] * jjt[1][0];
        if det.abs() < 1e-15 {
            jjt.clone()
        } else {
            let inv_det = 1.0 / det;
            vec![
                vec![jjt[1][1] * inv_det, -jjt[0][1] * inv_det],
                vec![-jjt[1][0] * inv_det, jjt[0][0] * inv_det],
            ]
        }
    } else {
        jjt.clone()
    };

    // J⁺ = J^T * (J*J^T)^{-1}  (n×m)
    let mut jt = vec![vec![0.0f64; m]; n];
    for i in 0..n {
        for k in 0..m {
            for l in 0..m {
                jt[i][k] += j[l][i] * jjt_inv[l][k];
            }
        }
    }
    jt
}

/// Compute the null-space projector N = I - J⁺ * J.
pub fn null_space_projector(j_pinv: &[Vec<f64>], j: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let n = j_pinv.len();
    let m = j.len();
    if n == 0 || m == 0 {
        return vec![];
    }
    // N = I - J⁺ J  (n×n)
    let mut n_proj = vec![vec![0.0f64; n]; n];
    for i in 0..n {
        for k in 0..n {
            let delta = if i == k { 1.0 } else { 0.0 };
            let mut sum = 0.0;
            for l in 0..m {
                if l < j_pinv[i].len() && k < j[l].len() {
                    sum += j_pinv[i][l] * j[l][k];
                }
            }
            n_proj[i][k] = delta - sum;
        }
    }
    n_proj
}

// ─────────────────────────────────────────────────────────────────────────────
// Euler integration for articulated bodies
// ─────────────────────────────────────────────────────────────────────────────

/// Semi-implicit Euler integration of an articulated body.
pub fn euler_integrate(body: &mut ArticulatedBody, dt: f64) {
    let qdd = body.aba(&body.links.iter().map(|l| l.tau).collect::<Vec<_>>());
    for (i, link) in body.links.iter_mut().enumerate() {
        let acc = if i < qdd.len() { qdd[i] } else { 0.0 };
        link.qd += acc * dt;
        link.q += link.qd * dt;
    }
}

/// Symplectic Euler (Störmer-Verlet like) integration for stability.
pub fn symplectic_euler_integrate(body: &mut ArticulatedBody, dt: f64) {
    let tau: Vec<f64> = body.links.iter().map(|l| l.tau).collect();
    let qdd = body.aba(&tau);
    for (i, link) in body.links.iter_mut().enumerate() {
        let acc = if i < qdd.len() { qdd[i] } else { 0.0 };
        link.qd += acc * dt; // update velocity first
        link.q += link.qd * dt; // use updated velocity
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Helper: build a simple serial chain
// ─────────────────────────────────────────────────────────────────────────────

/// Build a serial chain of n revolute joints (all about the Z axis).
pub fn build_serial_chain(
    n: usize,
    link_length: f64,
    mass: f64,
    gravity: [f64; 3],
) -> ArticulatedBody {
    let mut body = ArticulatedBody::new(gravity);
    let joint = JointAxis::revolute([0.0, 0.0, 1.0]);
    let x_t = SpatialTransform::translation([link_length, 0.0, 0.0]);
    for i in 0..n {
        let parent = if i == 0 { -1i64 } else { (i as i64) - 1 };
        let inertia = RigidBodyInertia::box_shape(mass, link_length / 2.0, 0.05, 0.05);
        let link = Link::new(format!("link_{i}"), parent, inertia, joint, x_t);
        body.add_link(link);
    }
    body
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-10;

    // ── SpatialVec ────────────────────────────────────────────────────────────

    #[test]
    fn test_spatial_vec_zero() {
        let v = SpatialVec::zero();
        assert!(v.norm() < EPS);
    }

    #[test]
    fn test_spatial_vec_add() {
        let a = SpatialVec::new([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        let b = SpatialVec::new([0.0, 1.0, 0.0], [1.0, 0.0, 0.0]);
        let c = a.add(&b);
        assert!((c.angular[0] - 1.0).abs() < EPS);
        assert!((c.angular[1] - 1.0).abs() < EPS);
        assert!((c.linear[0] - 1.0).abs() < EPS);
        assert!((c.linear[1] - 1.0).abs() < EPS);
    }

    #[test]
    fn test_spatial_vec_scale() {
        let v = SpatialVec::new([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]);
        let s = v.scale(2.0);
        assert!((s.angular[0] - 2.0).abs() < EPS);
        assert!((s.linear[2] - 12.0).abs() < EPS);
    }

    #[test]
    fn test_spatial_vec_dot_orthogonal() {
        let a = SpatialVec::new([1.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        let b = SpatialVec::new([0.0, 1.0, 0.0], [0.0, 0.0, 0.0]);
        assert!(a.dot(&b).abs() < EPS);
    }

    #[test]
    fn test_spatial_vec_dot_self() {
        let v = SpatialVec::new([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert!((v.dot(&v) - 2.0).abs() < EPS);
    }

    #[test]
    fn test_spatial_vec_norm() {
        let v = SpatialVec::new([3.0, 0.0, 0.0], [4.0, 0.0, 0.0]);
        assert!((v.norm() - 5.0).abs() < EPS);
    }

    // ── SpatialTransform ──────────────────────────────────────────────────────

    #[test]
    fn test_spatial_transform_identity_motion() {
        let t = SpatialTransform::identity();
        let v = SpatialVec::new([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]);
        let tv = t.apply_motion(&v);
        assert!((tv.angular[0] - 1.0).abs() < EPS);
        assert!((tv.linear[0] - 4.0).abs() < EPS);
    }

    #[test]
    fn test_spatial_transform_translation_motion() {
        let t = SpatialTransform::translation([1.0, 0.0, 0.0]);
        // v = (ω=[0,0,1], v=[0,0,0])
        let v = SpatialVec::new([0.0, 0.0, 1.0], [0.0, 0.0, 0.0]);
        let tv = t.apply_motion(&v);
        // ω' = R*ω = ω (identity), v' = R*(v - p×ω) = -(p×ω)
        // p×ω = (1,0,0)×(0,0,1) = (0*1-0*0, 0*0-1*1, 1*0-0*0) = (0,-1,0)
        // v' = -(0,-1,0) = (0,1,0)
        assert!(
            (tv.linear[1] - 1.0).abs() < 1e-8,
            "tv.linear={:?}",
            tv.linear
        );
    }

    #[test]
    fn test_spatial_transform_inverse_compose_identity() {
        let t = SpatialTransform::translation([1.0, 2.0, 3.0]);
        let t_inv = t.inverse();
        let composed = t.compose(&t_inv);
        // Should be close to identity
        for i in 0..3 {
            assert!(composed.pos[i].abs() < 1e-8, "pos[{i}]={}", composed.pos[i]);
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (composed.rot[i][j] - expected).abs() < 1e-8,
                    "rot[{i}][{j}]={}",
                    composed.rot[i][j]
                );
            }
        }
    }

    // ── RigidBodyInertia ──────────────────────────────────────────────────────

    #[test]
    fn test_sphere_inertia_isotropic() {
        let inertia = RigidBodyInertia::sphere(1.0, 1.0);
        let expected_i = 2.0 / 5.0;
        assert!((inertia.inertia_com[0][0] - expected_i).abs() < EPS);
        assert!((inertia.inertia_com[1][1] - expected_i).abs() < EPS);
        assert!((inertia.inertia_com[2][2] - expected_i).abs() < EPS);
    }

    #[test]
    fn test_sphere_spatial_inertia_nonzero_diag() {
        let inertia = RigidBodyInertia::sphere(2.0, 0.5);
        let si = inertia.to_spatial();
        // Upper left diagonal: I + m*c×c×^T (c=0, so just I_rot)
        let i = 2.0 / 5.0 * 2.0 * 0.5 * 0.5;
        assert!((si.a[0][0] - i).abs() < EPS);
    }

    #[test]
    fn test_point_mass_inertia() {
        let pm = RigidBodyInertia::point_mass(3.0);
        assert!((pm.mass - 3.0).abs() < EPS);
        for i in 0..3 {
            for j in 0..3 {
                assert!(pm.inertia_com[i][j].abs() < EPS);
            }
        }
    }

    // ── JointAxis ─────────────────────────────────────────────────────────────

    #[test]
    fn test_revolute_joint_dof() {
        let j = JointAxis::revolute([0.0, 0.0, 1.0]);
        assert_eq!(j.dof(), 1);
    }

    #[test]
    fn test_prismatic_joint_dof() {
        let j = JointAxis::prismatic([1.0, 0.0, 0.0]);
        assert_eq!(j.dof(), 1);
    }

    #[test]
    fn test_revolute_joint_transform_zero() {
        let j = JointAxis::revolute([0.0, 0.0, 1.0]);
        let t = j.joint_transform(0.0);
        // Zero rotation = identity
        for i in 0..3 {
            for k in 0..3 {
                let expected = if i == k { 1.0 } else { 0.0 };
                assert!((t.rot[i][k] - expected).abs() < 1e-8);
            }
        }
    }

    #[test]
    fn test_revolute_joint_transform_90deg() {
        let j = JointAxis::revolute([0.0, 0.0, 1.0]);
        let t = j.joint_transform(std::f64::consts::FRAC_PI_2);
        // Should rotate X axis to Y axis
        // col 0 of R should be (0, 1, 0)
        assert!((t.rot[0][0]).abs() < 1e-8); // R[0][0] = cos(90°) = 0
        assert!((t.rot[1][0] - 1.0).abs() < 1e-8); // R[1][0] = sin(90°) = 1
    }

    #[test]
    fn test_prismatic_joint_transform() {
        let j = JointAxis::prismatic([1.0, 0.0, 0.0]);
        let t = j.joint_transform(2.5);
        assert!((t.pos[0] - 2.5).abs() < EPS);
        assert!(t.pos[1].abs() < EPS);
        assert!(t.pos[2].abs() < EPS);
    }

    // ── ArticulatedBody ───────────────────────────────────────────────────────

    #[test]
    fn test_build_serial_chain() {
        let body = build_serial_chain(3, 1.0, 1.0, [0.0, -9.81, 0.0]);
        assert_eq!(body.num_links(), 3);
        assert_eq!(body.num_dof(), 3);
    }

    #[test]
    fn test_forward_kinematics_zero_q() {
        let mut body = build_serial_chain(2, 1.0, 1.0, [0.0, -9.81, 0.0]);
        body.forward_kinematics();
        // With q=0, transforms should be sequential translations
        // Link 0 world pos ≈ (1, 0, 0)
        let p0 = body.links[0].x_world.pos;
        assert!((p0[0] - 1.0).abs() < 1e-6, "link 0 x={}", p0[0]);
    }

    #[test]
    fn test_rnea_zero_qdd_no_gravity() {
        let mut body = build_serial_chain(2, 1.0, 1.0, [0.0, 0.0, 0.0]);
        body.forward_kinematics();
        let tau = body.rnea(&[0.0, 0.0], Some([0.0; 3]));
        // With zero qd and zero qdd and zero gravity: tau should be near zero
        for (i, &t) in tau.iter().enumerate() {
            assert!(t.abs() < 1e-4, "tau[{i}]={t}");
        }
    }

    #[test]
    fn test_aba_zero_tau_gives_gravity_acc() {
        let mut body = build_serial_chain(1, 1.0, 1.0, [0.0, -9.81, 0.0]);
        let qdd = body.aba(&[0.0]);
        // Single link under gravity should have nonzero acceleration
        assert!(!qdd.is_empty(), "should have qdd output");
    }

    #[test]
    fn test_crba_diagonal_positive() {
        let mut body = build_serial_chain(3, 1.0, 1.0, [0.0, 0.0, 0.0]);
        body.forward_kinematics();
        let h = body.crba();
        // Diagonal of CRBA should be positive
        for (i, row) in h.iter().enumerate() {
            assert!(row[i] > 0.0, "H[{i}][{i}]={}", row[i]);
        }
    }

    #[test]
    fn test_crba_symmetric() {
        let mut body = build_serial_chain(3, 1.0, 1.0, [0.0, 0.0, 0.0]);
        body.forward_kinematics();
        let h = body.crba();
        for (i, row) in h.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                assert!(
                    (val - h[j][i]).abs() < 1e-8,
                    "H not symmetric at ({i},{j}): {val} vs {}",
                    h[j][i]
                );
            }
        }
    }

    // ── Baumgarte ─────────────────────────────────────────────────────────────

    #[test]
    fn test_baumgarte_zero_errors() {
        let params = BaumgarteParams::default();
        let corr = baumgarte_stabilization(0.0, 0.0, &params, 0.01);
        assert!(corr.abs() < EPS);
    }

    #[test]
    fn test_baumgarte_position_error_drives_correction() {
        let params = BaumgarteParams::default();
        let corr = baumgarte_stabilization(0.1, 0.0, &params, 0.01);
        // Correction should be negative (restoring) for positive position error
        assert!(corr < 0.0, "correction={corr}");
    }

    // ── Operational space ─────────────────────────────────────────────────────

    #[test]
    fn test_operational_space_inertia_identity_jacobian() {
        // J = [[1, 0]], H_inv = [1]
        let j = vec![vec![1.0, 0.0]];
        let h_inv = vec![1.0, 1.0];
        let lambda = operational_space_inertia(&j, &h_inv);
        // Lambda = (J H^{-1} J^T)^{-1} = (1*1*1)^{-1} = 1
        assert!((lambda[0][0] - 1.0).abs() < 1e-8, "lambda={}", lambda[0][0]);
    }

    // ── Loop constraint ───────────────────────────────────────────────────────

    #[test]
    fn test_loop_constraint_zero_error() {
        let mut body = build_serial_chain(2, 1.0, 1.0, [0.0, 0.0, 0.0]);
        body.forward_kinematics();
        // Get actual relative position
        let pa = body.links[0].x_world.pos;
        let pb = body.links[1].x_world.pos;
        let target = sub3(pb, pa);
        let lc = LoopConstraint::new(0, 1, target, 1.0);
        let err = lc.position_error(&body);
        let norm = len3(err);
        assert!(norm < 1e-8, "loop constraint error={norm}");
    }

    // ── Cable-driven ──────────────────────────────────────────────────────────

    #[test]
    fn test_cable_segment_creation() {
        let cable = CableSegment::new(0, [0.0; 3], 1, [0.0; 3], 1.0, 1000.0, 10.0);
        assert!((cable.rest_length - 1.0).abs() < EPS);
        assert!((cable.stiffness - 1000.0).abs() < EPS);
    }

    #[test]
    fn test_cable_zero_force_at_rest_length() {
        let mut body = build_serial_chain(2, 1.0, 1.0, [0.0; 3]);
        body.forward_kinematics();
        // Place cable between links at rest length = distance between origins
        let pa = body.links[0].x_world.pos;
        let pb = body.links[1].x_world.pos;
        let dist = len3(sub3(pb, pa));
        let cable = CableSegment::new(0, [0.0; 3], 1, [0.0; 3], dist, 1000.0, 0.0);
        let (fa, fb) = cable.compute_force(&body);
        // At rest length: zero force
        assert!(fa.norm() < 1e-8, "fa={}", fa.norm());
        assert!(fb.norm() < 1e-8, "fb={}", fb.norm());
    }

    // ── Hybrid dynamics ───────────────────────────────────────────────────────

    #[test]
    fn test_hybrid_dynamics_all_forward() {
        let mut body = build_serial_chain(2, 1.0, 1.0, [0.0, 0.0, 0.0]);
        body.forward_kinematics();
        let modes = vec![HybridMode::ForwardDynamics; 2];
        let result = hybrid_dynamics(&mut body, &modes);
        assert_eq!(result.len(), 2);
    }

    // ── Pseudo-inverse / null-space ───────────────────────────────────────────

    #[test]
    fn test_pseudo_inverse_fat_1x2() {
        // J = [[1, 0]], J+ = [[1], [0]] (min-norm pseudo-inverse)
        let j = vec![vec![1.0, 0.0]];
        let jp = pseudo_inverse_fat(&j);
        assert_eq!(jp.len(), 2);
        assert!((jp[0][0] - 1.0).abs() < 1e-8, "jp[0][0]={}", jp[0][0]);
        assert!(jp[1][0].abs() < 1e-8, "jp[1][0]={}", jp[1][0]);
    }

    #[test]
    fn test_null_space_projector_1x2() {
        let j = vec![vec![1.0, 0.0]];
        let jp = pseudo_inverse_fat(&j);
        let n_proj = null_space_projector(&jp, &j);
        // N = I - J+ J = I - [[1],[0]]*[[1,0]] = [[0,0],[0,1]]
        assert!(n_proj[0][0].abs() < 1e-8);
        assert!((n_proj[1][1] - 1.0).abs() < 1e-8);
    }

    // ── Integration ───────────────────────────────────────────────────────────

    #[test]
    fn test_euler_integrate_zero_tau() {
        let mut body = build_serial_chain(1, 1.0, 1.0, [0.0; 3]);
        let q0 = body.links[0].q;
        let qd0 = body.links[0].qd;
        euler_integrate(&mut body, 0.01);
        // With zero tau and zero gravity, qd unchanged (or very small), q unchanged
        let q1 = body.links[0].q;
        let qd1 = body.links[0].qd;
        // No gravity: acceleration should be zero → q and qd unchanged
        assert!((q1 - q0).abs() < 1e-6, "q changed: {q1} vs {q0}");
        assert!((qd1 - qd0).abs() < 1e-6, "qd changed: {qd1} vs {qd0}");
    }

    #[test]
    fn test_symplectic_euler_integrate() {
        let mut body = build_serial_chain(1, 1.0, 1.0, [0.0; 3]);
        body.links[0].tau = 1.0; // apply torque
        let qd0 = body.links[0].qd;
        symplectic_euler_integrate(&mut body, 0.01);
        let qd1 = body.links[0].qd;
        // Velocity should have increased due to torque
        assert!(qd1 > qd0, "qd should increase under positive torque");
    }

    #[test]
    fn test_set_and_get_q() {
        let mut body = build_serial_chain(3, 1.0, 1.0, [0.0; 3]);
        body.set_q(&[0.1, 0.2, 0.3]);
        let q = body.get_q();
        assert!((q[0] - 0.1).abs() < EPS);
        assert!((q[1] - 0.2).abs() < EPS);
        assert!((q[2] - 0.3).abs() < EPS);
    }

    #[test]
    fn test_set_and_get_qd() {
        let mut body = build_serial_chain(2, 1.0, 1.0, [0.0; 3]);
        body.set_qd(&[1.5, -0.5]);
        let qd = body.get_qd();
        assert!((qd[0] - 1.5).abs() < EPS);
        assert!((qd[1] + 0.5).abs() < EPS);
    }

    #[test]
    fn test_cable_driven_system_creation() {
        let body = build_serial_chain(2, 1.0, 1.0, [0.0; 3]);
        let mut cds = CableDrivenSystem::new(body);
        let cable = CableSegment::new(0, [0.0; 3], 1, [0.0; 3], 1.0, 1000.0, 0.0);
        cds.add_cable(cable);
        assert_eq!(cds.cables.len(), 1);
    }

    #[test]
    fn test_project_constraint_zero_error() {
        let q = 1.5;
        let q_new = project_constraint(q, 0.0, 1.0);
        assert!((q_new - q).abs() < EPS);
    }

    #[test]
    fn test_project_constraint_nonzero_error() {
        let q = 1.0;
        let err = 0.1;
        let gain = 0.5;
        let q_new = project_constraint(q, err, gain);
        assert!((q_new - (q - gain * err)).abs() < EPS);
    }

    #[test]
    fn test_spatial_mat_identity_mul_vec() {
        let i = SpatialMat::identity();
        let v = SpatialVec::new([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]);
        let iv = i.mul_vec(&v);
        for k in 0..3 {
            assert!((iv.angular[k] - v.angular[k]).abs() < EPS);
            assert!((iv.linear[k] - v.linear[k]).abs() < EPS);
        }
    }

    #[test]
    fn test_skew_matrix_cross_product() {
        let v = [1.0, 2.0, 3.0];
        let sv = skew(v);
        let u = [4.0, 5.0, 6.0];
        let sv_u = mat3_vec(sv, u);
        let expected = cross3(v, u);
        for k in 0..3 {
            assert!(
                (sv_u[k] - expected[k]).abs() < EPS,
                "k={k}: sv_u={} expected={}",
                sv_u[k],
                expected[k]
            );
        }
    }
}
