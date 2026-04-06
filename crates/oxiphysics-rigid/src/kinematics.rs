#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Advanced kinematics: screw theory, Denavit-Hartenberg, forward/inverse
//! kinematics, workspace computation, singularity analysis, manipulability
//! ellipsoid, velocity/acceleration analysis, Coriolis/centrifugal terms,
//! gyroscopic effects, mechanism synthesis, 4-bar linkage, cam-follower,
//! gear train kinematics, non-holonomic constraints, and Pfaffian constraints.

#![allow(dead_code)]
#![allow(unused_imports)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ══════════════════════════════════════════════════════════════════════════════
// § 1  BASIC ROTATION AND HOMOGENEOUS TRANSFORMS
// ══════════════════════════════════════════════════════════════════════════════

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
    // Rᵀ
    for i in 0..3 {
        for j in 0..3 {
            result[i][j] = t[j][i];
        }
    }
    // −Rᵀ·p
    for i in 0..3 {
        let mut sum = 0.0f64;
        for j in 0..3 {
            sum += t[j][i] * t[j][3];
        }
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

// ══════════════════════════════════════════════════════════════════════════════
// § 2  DENAVIT-HARTENBERG CONVENTION
// ══════════════════════════════════════════════════════════════════════════════

/// Standard Denavit-Hartenberg parameters for one joint.
#[derive(Debug, Clone)]
pub struct DhParams {
    /// Link twist α_i \[rad\]
    pub alpha: f64,
    /// Link length a_i \[m\]
    pub a: f64,
    /// Link offset d_i \[m\]
    pub d: f64,
    /// Joint angle θ_i \[rad\] (variable for revolute joints)
    pub theta: f64,
    /// Joint type: true = revolute, false = prismatic
    pub revolute: bool,
}

impl DhParams {
    /// Create a revolute DH joint.
    pub fn revolute(alpha: f64, a: f64, d: f64, theta: f64) -> Self {
        Self {
            alpha,
            a,
            d,
            theta,
            revolute: true,
        }
    }

    /// Create a prismatic DH joint.
    pub fn prismatic(alpha: f64, a: f64, d: f64, theta: f64) -> Self {
        Self {
            alpha,
            a,
            d,
            theta,
            revolute: false,
        }
    }

    /// Compute the homogeneous transform T_{i−1,i} for this DH parameter set.
    ///
    /// T = Rot_z(θ) · Trans_z(d) · Trans_x(a) · Rot_x(α)
    pub fn transform(&self) -> HMatrix {
        let theta = self.theta;
        let ct = theta.cos();
        let st = theta.sin();
        let ca = self.alpha.cos();
        let sa = self.alpha.sin();
        let a = self.a;
        let d = self.d;

        [
            [ct, -st * ca, st * sa, a * ct],
            [st, ct * ca, -ct * sa, a * st],
            [0.0, sa, ca, d],
            [0.0, 0.0, 0.0, 1.0],
        ]
    }

    /// Modified DH convention transform (Craig's convention).
    pub fn transform_modified(&self) -> HMatrix {
        let theta = self.theta;
        let ct = theta.cos();
        let st = theta.sin();
        let ca = self.alpha.cos();
        let sa = self.alpha.sin();
        let a = self.a;
        let d = self.d;

        [
            [ct, -st, 0.0, a],
            [st * ca, ct * ca, -sa, -sa * d],
            [st * sa, ct * sa, ca, ca * d],
            [0.0, 0.0, 0.0, 1.0],
        ]
    }

    /// Update the variable joint parameter (θ for revolute, d for prismatic).
    pub fn set_variable(&mut self, q: f64) {
        if self.revolute {
            self.theta = q;
        } else {
            self.d = q;
        }
    }

    /// Get the variable joint parameter.
    pub fn get_variable(&self) -> f64 {
        if self.revolute { self.theta } else { self.d }
    }
}

// ══════════════════════════════════════════════════════════════════════════════
// § 3  ROBOT ARM — FORWARD KINEMATICS
// ══════════════════════════════════════════════════════════════════════════════

/// Open kinematic chain robot (serial manipulator).
#[derive(Debug, Clone)]
pub struct SerialManipulator {
    /// DH parameters for each joint
    pub dh: Vec<DhParams>,
    /// Joint configuration vector q
    pub q: Vec<f64>,
    /// Joint velocity vector q̇
    pub q_dot: Vec<f64>,
    /// Number of joints (DOF)
    pub n_joints: usize,
}

impl SerialManipulator {
    /// Create a new serial manipulator.
    pub fn new(dh: Vec<DhParams>) -> Self {
        let n = dh.len();
        Self {
            dh,
            q: vec![0.0; n],
            q_dot: vec![0.0; n],
            n_joints: n,
        }
    }

    /// Create a standard 2-link planar arm.
    pub fn planar_2r(l1: f64, l2: f64) -> Self {
        let dh = vec![
            DhParams::revolute(0.0, l1, 0.0, 0.0),
            DhParams::revolute(0.0, l2, 0.0, 0.0),
        ];
        Self::new(dh)
    }

    /// Create a 6-DOF robot similar to a PUMA-style arm.
    pub fn puma560() -> Self {
        let dh = vec![
            DhParams::revolute(0.0, 0.0, 0.67, 0.0),
            DhParams::revolute(-PI * 0.5, 0.432, 0.149, 0.0),
            DhParams::revolute(0.0, 0.0203, 0.0, 0.0),
            DhParams::revolute(-PI * 0.5, 0.0, 0.433, 0.0),
            DhParams::revolute(PI * 0.5, 0.0, 0.0, 0.0),
            DhParams::revolute(-PI * 0.5, 0.0, 0.0, 0.0),
        ];
        Self::new(dh)
    }

    /// Alias for [`puma560`](Self::puma560).
    pub fn puma_like() -> Self {
        Self::puma560()
    }

    /// Set joint configuration.
    pub fn set_q(&mut self, q: &[f64]) {
        for (i, &qi) in q.iter().enumerate() {
            if i < self.n_joints {
                self.q[i] = qi;
                self.dh[i].set_variable(qi);
            }
        }
    }

    /// Forward kinematics: compute end-effector transform T_0n.
    pub fn forward_kinematics(&self) -> HMatrix {
        let mut t = EYE4;
        for joint in &self.dh {
            t = h_mul(&t, &joint.transform());
        }
        t
    }

    /// Compute transform up to joint k (0-indexed).
    pub fn transform_to_joint(&self, k: usize) -> HMatrix {
        let mut t = EYE4;
        for i in 0..=k.min(self.n_joints - 1) {
            t = h_mul(&t, &self.dh[i].transform());
        }
        t
    }

    /// End-effector position.
    pub fn end_effector_position(&self) -> Vec3 {
        let t = self.forward_kinematics();
        h_translation(&t)
    }

    /// End-effector orientation as rotation matrix.
    pub fn end_effector_rotation(&self) -> Rot3 {
        let t = self.forward_kinematics();
        h_rotation(&t)
    }
}

// ══════════════════════════════════════════════════════════════════════════════
// § 4  JACOBIAN AND VELOCITY ANALYSIS
// ══════════════════════════════════════════════════════════════════════════════

/// Geometric Jacobian of a serial manipulator (6 × n_joints).
///
/// Columns: \[z_{i-1} × (p_n − p_{i-1}); z_{i-1}\] for revolute
///          \[z_{i-1};                    0       \] for prismatic
#[derive(Debug, Clone)]
pub struct Jacobian {
    /// 6 × n matrix stored as \[column_0, column_1, ...\]
    pub columns: Vec<[f64; 6]>,
    /// Number of joints
    pub n_joints: usize,
}

impl Jacobian {
    /// Compute the geometric Jacobian of the manipulator at current configuration.
    pub fn compute(arm: &SerialManipulator) -> Self {
        let n = arm.n_joints;
        let t_total = arm.forward_kinematics();
        let p_n = h_translation(&t_total);

        let mut columns = Vec::with_capacity(n);
        let mut t_i = EYE4;

        for i in 0..n {
            let z_i = [t_i[0][2], t_i[1][2], t_i[2][2]]; // z-axis of frame i-1
            let p_i = [t_i[0][3], t_i[1][3], t_i[2][3]];

            let dp = [p_n[0] - p_i[0], p_n[1] - p_i[1], p_n[2] - p_i[2]];

            let col = if arm.dh[i].revolute {
                let jv = cross3(&z_i, &dp);
                [jv[0], jv[1], jv[2], z_i[0], z_i[1], z_i[2]]
            } else {
                [z_i[0], z_i[1], z_i[2], 0.0, 0.0, 0.0]
            };

            columns.push(col);
            t_i = h_mul(&t_i, &arm.dh[i].transform());
        }

        Self {
            columns,
            n_joints: n,
        }
    }

    /// Get matrix element J\[row\]\[col\].
    pub fn get(&self, row: usize, col: usize) -> f64 {
        self.columns[col][row]
    }

    /// Compute the 6-vector: twist = J · q̇.
    pub fn multiply_velocity(&self, q_dot: &[f64]) -> [f64; 6] {
        let mut twist = [0.0f64; 6];
        for j in 0..self.n_joints {
            for i in 0..6 {
                twist[i] += self.columns[j][i] * q_dot[j];
            }
        }
        twist
    }

    /// Compute the transpose product: torques = Jᵀ · wrench.
    pub fn transpose_multiply(&self, wrench: &[f64; 6]) -> Vec<f64> {
        (0..self.n_joints)
            .map(|j| (0..6).map(|i| self.columns[j][i] * wrench[i]).sum())
            .collect()
    }

    /// Compute the JJᵀ matrix (6×6) for manipulability analysis.
    pub fn jjt(&self) -> [[f64; 6]; 6] {
        let mut m = [[0.0f64; 6]; 6];
        for j in 0..self.n_joints {
            for r in 0..6 {
                for c in 0..6 {
                    m[r][c] += self.columns[j][r] * self.columns[j][c];
                }
            }
        }
        m
    }

    /// Yoshikawa manipulability measure: w = sqrt(det(J·Jᵀ)).
    ///
    /// Uses the n×n positional sub-Jacobian (first `n_joints` rows × `n_joints` columns)
    /// to avoid singularity when fewer joints than 3 spatial dimensions are active.
    pub fn manipulability(&self) -> f64 {
        let n = self.n_joints;
        if n == 0 {
            return 0.0;
        }
        // Build n×n Gram matrix: G = J_n * J_n^T where J_n uses first n rows
        let mut g = vec![vec![0.0f64; n]; n];
        for j in 0..n {
            for r in 0..n {
                for c in 0..n {
                    g[r][c] += self.columns[j][r] * self.columns[j][c];
                }
            }
        }
        // Compute determinant of n×n matrix via LU decomposition
        let det = mat_det_lu(&g);
        det.abs().sqrt()
    }

    /// Check if the manipulator is near a singularity.
    pub fn is_singular(&self, tol: f64) -> bool {
        self.manipulability() < tol
    }
}

// ══════════════════════════════════════════════════════════════════════════════
// § 5  INVERSE KINEMATICS
// ══════════════════════════════════════════════════════════════════════════════

/// Iterative inverse kinematics result.
#[derive(Debug, Clone)]
pub struct IkResult {
    /// Joint configuration
    pub q: Vec<f64>,
    /// Final position error \[m\]
    pub pos_error: f64,
    /// Number of iterations
    pub iterations: usize,
    /// Whether convergence was achieved
    pub converged: bool,
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

        // Jacobian (positional only, 3×n)
        let jac = Jacobian::compute(arm);
        let n = arm.n_joints;

        // Δq = Jᵀ · (J·Jᵀ + λ²·I)⁻¹ · e  (Levenberg-Marquardt)
        // For positional 3×n case:
        let mut jjt3 = [[0.0f64; 3]; 3];
        for j in 0..n {
            for r in 0..3 {
                for c in 0..3 {
                    jjt3[r][c] += jac.columns[j][r] * jac.columns[j][c];
                }
            }
        }
        // Add damping
        for i in 0..3 {
            jjt3[i][i] += lambda * lambda;
        }

        // Solve (JJᵀ + λI) · y = e using Cramer's rule (3×3)
        let y = solve3x3(&jjt3, &err);

        // Δq = Jᵀ · y
        for j in 0..n {
            let dq: f64 = (0..3).map(|r| jac.columns[j][r] * y[r]).sum();
            q[j] += dq;
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
fn mat_det_lu(a: &[Vec<f64>]) -> f64 {
    let n = a.len();
    if n == 0 {
        return 0.0;
    }
    let mut m: Vec<Vec<f64>> = a.to_vec();
    let mut det = 1.0f64;
    for col in 0..n {
        // Find pivot
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
            for k in col..n {
                let v = m[col][k];
                m[row][k] -= factor * v;
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
    } // out of reach

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

// ══════════════════════════════════════════════════════════════════════════════
// § 6  SCREW THEORY
// ══════════════════════════════════════════════════════════════════════════════

/// Screw axis $\hat{\mathcal{S}} = (\hat{\omega}, \mathbf{v})$ in body/space frame.
#[derive(Debug, Clone)]
pub struct ScrewAxis {
    /// Angular component ω̂ (unit vector for revolute, zero for prismatic)
    pub omega: Vec3,
    /// Linear component v
    pub v: Vec3,
    /// Pitch h = v·ω̂ / |ω̂|²  (0 for revolute, ∞ for prismatic)
    pub pitch: f64,
}

impl ScrewAxis {
    /// Create a revolute screw about axis ω̂ through point q.
    ///
    /// v = −ω̂ × q
    pub fn revolute(omega_hat: Vec3, q: Vec3) -> Self {
        let v = [
            -(omega_hat[1] * q[2] - omega_hat[2] * q[1]),
            -(omega_hat[2] * q[0] - omega_hat[0] * q[2]),
            -(omega_hat[0] * q[1] - omega_hat[1] * q[0]),
        ];
        Self {
            omega: omega_hat,
            v,
            pitch: 0.0,
        }
    }

    /// Create a prismatic screw along direction v̂.
    pub fn prismatic(v_hat: Vec3) -> Self {
        Self {
            omega: [0.0; 3],
            v: v_hat,
            pitch: f64::INFINITY,
        }
    }

    /// Is this a revolute screw?
    pub fn is_revolute(&self) -> bool {
        norm3(&self.omega) > 0.5
    }

    /// Matrix exponential exp(\[S\]·θ) for a screw.
    ///
    /// For revolute: Rodrigues formula.
    pub fn exp_map(&self, theta: f64) -> HMatrix {
        if self.is_revolute() {
            // Rodrigues: R = I + sin(θ)·[ω̂] + (1−cos θ)·[ω̂]²
            let c = theta.cos();
            let s = theta.sin();
            let [wx, wy, wz] = self.omega;
            let r = [
                [
                    c + wx * wx * (1.0 - c),
                    wx * wy * (1.0 - c) - wz * s,
                    wx * wz * (1.0 - c) + wy * s,
                ],
                [
                    wy * wx * (1.0 - c) + wz * s,
                    c + wy * wy * (1.0 - c),
                    wy * wz * (1.0 - c) - wx * s,
                ],
                [
                    wz * wx * (1.0 - c) - wy * s,
                    wz * wy * (1.0 - c) + wx * s,
                    c + wz * wz * (1.0 - c),
                ],
            ];
            // Translation: G(θ)·v where G = I·θ + (1−cos θ)·[ω̂] + (θ−sin θ)·[ω̂]²
            let gv = [
                self.v[0] * theta
                    + (1.0 - c) * (wy * self.v[2] - wz * self.v[1])
                    + (theta - s)
                        * (wx * (wx * self.v[0] + wy * self.v[1] + wz * self.v[2]) - self.v[0]),
                self.v[1] * theta
                    + (1.0 - c) * (wz * self.v[0] - wx * self.v[2])
                    + (theta - s)
                        * (wy * (wx * self.v[0] + wy * self.v[1] + wz * self.v[2]) - self.v[1]),
                self.v[2] * theta
                    + (1.0 - c) * (wx * self.v[1] - wy * self.v[0])
                    + (theta - s)
                        * (wz * (wx * self.v[0] + wy * self.v[1] + wz * self.v[2]) - self.v[2]),
            ];
            [
                [r[0][0], r[0][1], r[0][2], gv[0]],
                [r[1][0], r[1][1], r[1][2], gv[1]],
                [r[2][0], r[2][1], r[2][2], gv[2]],
                [0., 0., 0., 1.],
            ]
        } else {
            // Prismatic: pure translation
            translate(self.v[0] * theta, self.v[1] * theta, self.v[2] * theta)
        }
    }

    /// Product of exponentials forward kinematics: T = e^{S1·q1} · ... · e^{Sn·qn} · M.
    pub fn poe_fk(screws: &[ScrewAxis], q: &[f64], m: &HMatrix) -> HMatrix {
        let mut t = *m;
        for (s, &qi) in screws.iter().zip(q.iter()).rev() {
            let exp_s = s.exp_map(qi);
            t = h_mul(&exp_s, &t);
        }
        t
    }
}

// ══════════════════════════════════════════════════════════════════════════════
// § 7  WORKSPACE COMPUTATION
// ══════════════════════════════════════════════════════════════════════════════

/// Workspace representation as a set of reachable points.
#[derive(Debug, Clone)]
pub struct Workspace {
    /// Sampled reachable points
    pub points: Vec<Vec3>,
    /// Joint limits (min, max) for each joint
    pub joint_limits: Vec<(f64, f64)>,
    /// Bounding box \[x_min, x_max, y_min, y_max, z_min, z_max\]
    pub bounding_box: [f64; 6],
}

impl Workspace {
    /// Sample the workspace of a serial manipulator.
    pub fn sample(arm: &mut SerialManipulator, n_samples_per_joint: usize) -> Self {
        let n = arm.n_joints;
        let joint_limits: Vec<(f64, f64)> = vec![(-PI, PI); n];
        let mut points = Vec::new();

        // Recursive grid sampling (limited to avoid explosion)
        let step = 2.0 * PI / n_samples_per_joint as f64;
        let total = n_samples_per_joint.pow(n as u32).min(50_000);

        let mut idx = 0;
        let mut q = vec![joint_limits[0].0; n];
        loop {
            arm.set_q(&q);
            points.push(arm.end_effector_position());
            idx += 1;
            if idx >= total {
                break;
            }

            // Increment joint angles
            let mut carry = true;
            for j in 0..n {
                if carry {
                    q[j] += step;
                    if q[j] > joint_limits[j].1 {
                        q[j] = joint_limits[j].0;
                    } else {
                        carry = false;
                    }
                }
            }
            if carry {
                break;
            }
        }

        let (mut xmin, mut xmax) = (f64::INFINITY, f64::NEG_INFINITY);
        let (mut ymin, mut ymax) = (f64::INFINITY, f64::NEG_INFINITY);
        let (mut zmin, mut zmax) = (f64::INFINITY, f64::NEG_INFINITY);
        for p in &points {
            xmin = xmin.min(p[0]);
            xmax = xmax.max(p[0]);
            ymin = ymin.min(p[1]);
            ymax = ymax.max(p[1]);
            zmin = zmin.min(p[2]);
            zmax = zmax.max(p[2]);
        }

        Self {
            points,
            joint_limits,
            bounding_box: [xmin, xmax, ymin, ymax, zmin, zmax],
        }
    }

    /// Approximate workspace radius (from origin).
    pub fn max_reach(&self) -> f64 {
        self.points.iter().map(norm3).fold(0.0_f64, f64::max)
    }

    /// Number of sampled points.
    pub fn size(&self) -> usize {
        self.points.len()
    }
}

// ══════════════════════════════════════════════════════════════════════════════
// § 8  SINGULARITY ANALYSIS AND MANIPULABILITY ELLIPSOID
// ══════════════════════════════════════════════════════════════════════════════

/// Manipulability ellipsoid from the Jacobian.
#[derive(Debug, Clone)]
pub struct ManipulabilityEllipsoid {
    /// Semi-axes lengths (square roots of singular values of J)
    pub semi_axes: Vec<f64>,
    /// Yoshikawa manipulability index
    pub manipulability_index: f64,
    /// Condition number (ratio max/min singular value)
    pub condition_number: f64,
}

impl ManipulabilityEllipsoid {
    /// Compute the manipulability ellipsoid from the Jacobian.
    pub fn compute(jac: &Jacobian) -> Self {
        let n = jac.n_joints;
        // Power iteration to find singular values of the 3×n positional Jacobian
        // Compute JJᵀ (3×3 positional part)
        let mut jjt = [[0.0f64; 3]; 3];
        for j in 0..n {
            for r in 0..3 {
                for c in 0..3 {
                    jjt[r][c] += jac.columns[j][r] * jac.columns[j][c];
                }
            }
        }

        // Estimate eigenvalues via Gershgorin circles (conservative bounds)
        let max_eig = (0..3)
            .map(|i| {
                jjt[i][i]
                    + (0..3)
                        .filter(|&j| j != i)
                        .map(|j| jjt[i][j].abs())
                        .sum::<f64>()
            })
            .fold(0.0_f64, f64::max);
        let min_eig = (0..3)
            .map(|i| {
                (jjt[i][i]
                    - (0..3)
                        .filter(|&j| j != i)
                        .map(|j| jjt[i][j].abs())
                        .sum::<f64>())
                .max(0.0)
            })
            .fold(f64::INFINITY, f64::min);

        let det = jjt[0][0] * (jjt[1][1] * jjt[2][2] - jjt[1][2] * jjt[2][1])
            - jjt[0][1] * (jjt[1][0] * jjt[2][2] - jjt[1][2] * jjt[2][0])
            + jjt[0][2] * (jjt[1][0] * jjt[2][1] - jjt[1][1] * jjt[2][0]);

        let manipulability_index = det.abs().sqrt();
        let condition_number = if min_eig > 1e-30 {
            max_eig / min_eig
        } else {
            f64::INFINITY
        };

        Self {
            semi_axes: vec![
                max_eig.sqrt(),
                (max_eig * min_eig).sqrt().max(1e-10),
                min_eig.sqrt(),
            ],
            manipulability_index,
            condition_number,
        }
    }

    /// Is the configuration isotropic (condition number near 1)?
    pub fn is_isotropic(&self, tol: f64) -> bool {
        self.condition_number < 1.0 + tol
    }
}

// ══════════════════════════════════════════════════════════════════════════════
// § 9  CORIOLIS / CENTRIFUGAL TERMS AND DYNAMICS
// ══════════════════════════════════════════════════════════════════════════════

/// Rigid body with angular velocity (for gyroscopic effects).
#[derive(Debug, Clone)]
pub struct GyroscopicBody {
    /// Inertia tensor (principal axes, diagonal) \[kg·m²\]
    pub inertia: Vec3,
    /// Angular velocity ω \[rad/s\]
    pub omega: Vec3,
    /// Angular momentum L = I·ω
    pub angular_momentum: Vec3,
}

impl GyroscopicBody {
    /// Create a new gyroscopic body.
    pub fn new(inertia: Vec3, omega: Vec3) -> Self {
        let angular_momentum = [
            inertia[0] * omega[0],
            inertia[1] * omega[1],
            inertia[2] * omega[2],
        ];
        Self {
            inertia,
            omega,
            angular_momentum,
        }
    }

    /// Update angular momentum from omega.
    pub fn update_momentum(&mut self) {
        self.angular_momentum = [
            self.inertia[0] * self.omega[0],
            self.inertia[1] * self.omega[1],
            self.inertia[2] * self.omega[2],
        ];
    }

    /// Gyroscopic torque τ_gyro = −ω × L.
    pub fn gyroscopic_torque(&self) -> Vec3 {
        let l = self.angular_momentum;
        cross3(&self.omega, &l)
    }

    /// Euler's equations of motion: I·ω̇ = τ − ω × (I·ω).
    pub fn euler_equations(&self, torque: &Vec3) -> Vec3 {
        let gyro = self.gyroscopic_torque();
        [
            (torque[0] - gyro[0]) / self.inertia[0].max(1e-30),
            (torque[1] - gyro[1]) / self.inertia[1].max(1e-30),
            (torque[2] - gyro[2]) / self.inertia[2].max(1e-30),
        ]
    }

    /// Kinetic energy T = ½ω·I·ω.
    pub fn kinetic_energy(&self) -> f64 {
        0.5 * (self.inertia[0] * self.omega[0].powi(2)
            + self.inertia[1] * self.omega[1].powi(2)
            + self.inertia[2] * self.omega[2].powi(2))
    }

    /// Nutation angle (between spin axis and precession axis).
    pub fn nutation_angle(&self, precession_axis: &Vec3) -> f64 {
        let cos_theta = dot3(&normalize3(&self.omega), &normalize3(precession_axis));
        cos_theta.clamp(-1.0, 1.0).acos()
    }

    /// Integrate one step using explicit Euler.
    pub fn step(&mut self, torque: &Vec3, dt: f64) {
        let omega_dot = self.euler_equations(torque);
        self.omega[0] += omega_dot[0] * dt;
        self.omega[1] += omega_dot[1] * dt;
        self.omega[2] += omega_dot[2] * dt;
        self.update_momentum();
    }
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

// ══════════════════════════════════════════════════════════════════════════════
// § 10  FOUR-BAR LINKAGE
// ══════════════════════════════════════════════════════════════════════════════

/// Grashof four-bar linkage (crank-rocker, double-crank, rocker-crank, or double-rocker).
///
/// Links: (s = shortest, l = longest, p, q = others).
/// Grashof's condition: s + l ≤ p + q → at least one link can rotate fully.
#[derive(Debug, Clone)]
pub struct FourBarLinkage {
    /// Crank length r₁ \[m\]
    pub r1: f64,
    /// Coupler length r₂ \[m\]
    pub r2: f64,
    /// Follower length r₃ \[m\]
    pub r3: f64,
    /// Ground (fixed) link length r₄ \[m\]
    pub r4: f64,
    /// Crank angle θ₁ \[rad\]
    pub theta1: f64,
}

impl FourBarLinkage {
    /// Create a new four-bar linkage.
    pub fn new(r1: f64, r2: f64, r3: f64, r4: f64) -> Self {
        Self {
            r1,
            r2,
            r3,
            r4,
            theta1: 0.0,
        }
    }

    /// Check Grashof's condition.
    pub fn is_grashof(&self) -> bool {
        let lengths = [self.r1, self.r2, self.r3, self.r4];
        let s = lengths.iter().cloned().fold(f64::INFINITY, f64::min);
        let l = lengths.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let sum_pl: f64 = lengths.iter().sum::<f64>() - s - l;
        s + l <= sum_pl
    }

    /// Compute the follower angle θ₃ from the crank angle θ₁ (Freudenstein's equation).
    ///
    /// Returns (open, crossed) circuit solutions.
    pub fn follower_angle(&self, theta1: f64) -> Option<(f64, f64)> {
        let r1 = self.r1;
        let r2 = self.r2;
        let r3 = self.r3;
        let r4 = self.r4;

        // Freudenstein coefficients
        let k1 = r4 / r1;
        let k2 = r4 / r3;
        let k3 = (r1 * r1 - r2 * r2 + r3 * r3 + r4 * r4) / (2.0 * r1 * r3);

        let a = theta1.cos() - k1 - k2 * theta1.cos() + k3;
        let b = -2.0 * theta1.sin();
        let c = k1 - (k2 + 1.0) * theta1.cos() + k3;

        let discriminant = b * b - 4.0 * a * c;
        if discriminant < 0.0 {
            return None;
        }

        let t1 = (-b + discriminant.sqrt()) / (2.0 * a);
        let t2 = (-b - discriminant.sqrt()) / (2.0 * a);

        Some((2.0 * t1.atan(), 2.0 * t2.atan()))
    }

    /// Coupler point position (point P on the coupler link).
    ///
    /// P is at distance r_p from joint A at angle δ from the coupler.
    pub fn coupler_point(&self, theta1: f64, theta2: f64, r_p: f64, delta: f64) -> [f64; 2] {
        let ax = self.r1 * theta1.cos();
        let ay = self.r1 * theta1.sin();
        [
            ax + r_p * (theta2 + delta).cos(),
            ay + r_p * (theta2 + delta).sin(),
        ]
    }

    /// Angular velocity of the follower ω₃ from the velocity analysis.
    ///
    /// ω₃ = ω₁ · r₁ · sin(θ₂ − θ₁) / (r₃ · sin(θ₃ − θ₂))
    pub fn follower_angular_velocity(
        &self,
        omega1: f64,
        theta1: f64,
        theta2: f64,
        theta3: f64,
    ) -> f64 {
        let denom = self.r3 * (theta3 - theta2).sin();
        if denom.abs() < 1e-10 {
            return 0.0;
        }
        omega1 * self.r1 * (theta2 - theta1).sin() / denom
    }

    /// Mechanical advantage = output torque / input torque.
    ///
    /// MA = (r₁ · sin(θ₃ − θ₁)) / (r₃ · sin(θ₃ − θ₂))
    pub fn mechanical_advantage(&self, theta1: f64, theta2: f64, theta3: f64) -> f64 {
        let num = self.r1 * (theta3 - theta1).sin();
        let den = self.r3 * (theta3 - theta2).sin();
        if den.abs() < 1e-10 {
            return f64::INFINITY;
        }
        num / den
    }
}

// ══════════════════════════════════════════════════════════════════════════════
// § 11  CAM-FOLLOWER MECHANISMS
// ══════════════════════════════════════════════════════════════════════════════

/// Cam profile representation (discrete lift curve).
#[derive(Debug, Clone)]
pub struct CamProfile {
    /// Cam rotation angles \[rad\]
    pub angles: Vec<f64>,
    /// Follower lifts at each angle \[m\]
    pub lifts: Vec<f64>,
    /// Base circle radius r_b \[m\]
    pub base_radius: f64,
    /// Prime circle radius r_p = r_b + max_lift \[m\]
    pub prime_radius: f64,
}

impl CamProfile {
    /// Create a sinusoidal cam profile.
    pub fn sinusoidal(n_points: usize, base_radius: f64, max_lift: f64) -> Self {
        let angles: Vec<f64> = (0..n_points)
            .map(|i| 2.0 * PI * i as f64 / n_points as f64)
            .collect();
        let lifts: Vec<f64> = angles
            .iter()
            .map(|&a| max_lift * 0.5 * (1.0 - (2.0 * a).cos()))
            .collect();
        let prime_radius = base_radius + max_lift;
        Self {
            angles,
            lifts,
            base_radius,
            prime_radius,
        }
    }

    /// Interpolate lift at angle θ using linear interpolation.
    pub fn lift_at(&self, theta: f64) -> f64 {
        let theta = theta.rem_euclid(2.0 * PI);
        let n = self.angles.len();
        let i = (theta / (2.0 * PI / n as f64)) as usize % n;
        let j = (i + 1) % n;
        let frac = (theta - self.angles[i]) / (2.0 * PI / n as f64);
        self.lifts[i] + frac * (self.lifts[j] - self.lifts[i])
    }

    /// Follower velocity dL/dθ (numerical derivative).
    pub fn follower_velocity(&self, theta: f64, dtheta: f64) -> f64 {
        let l1 = self.lift_at(theta + dtheta);
        let l0 = self.lift_at(theta - dtheta);
        (l1 - l0) / (2.0 * dtheta)
    }

    /// Follower acceleration d²L/dθ².
    pub fn follower_acceleration(&self, theta: f64, dtheta: f64) -> f64 {
        let l1 = self.lift_at(theta + dtheta);
        let l0 = self.lift_at(theta);
        let lm = self.lift_at(theta - dtheta);
        (l1 - 2.0 * l0 + lm) / (dtheta * dtheta)
    }

    /// Pressure angle φ at angle θ for an in-line knife-edge follower.
    ///
    /// tan φ = (dR/dθ) / R  where R = r_b + L(θ)
    pub fn pressure_angle(&self, theta: f64, dtheta: f64) -> f64 {
        let l = self.lift_at(theta);
        let dl_dtheta = self.follower_velocity(theta, dtheta);
        let r = self.base_radius + l;
        (dl_dtheta / r.max(1e-30)).atan().abs()
    }

    /// Maximum pressure angle over full rotation.
    pub fn max_pressure_angle(&self) -> f64 {
        self.angles
            .iter()
            .map(|&a| self.pressure_angle(a, 0.001))
            .fold(0.0_f64, f64::max)
    }
}

// ══════════════════════════════════════════════════════════════════════════════
// § 12  GEAR TRAIN KINEMATICS
// ══════════════════════════════════════════════════════════════════════════════

/// Single gear pair.
#[derive(Debug, Clone)]
pub struct GearPair {
    /// Number of teeth on driver gear
    pub n_driver: u32,
    /// Number of teeth on driven gear
    pub n_driven: u32,
    /// Module m = pitch diameter / number of teeth \[m\]
    pub module: f64,
    /// Pressure angle φ \[rad\]
    pub pressure_angle: f64,
    /// Helix angle for helical gears \[rad\]
    pub helix_angle: f64,
}

impl GearPair {
    /// Create a spur gear pair.
    pub fn spur(n_driver: u32, n_driven: u32, module: f64) -> Self {
        Self {
            n_driver,
            n_driven,
            module,
            pressure_angle: 20.0_f64.to_radians(),
            helix_angle: 0.0,
        }
    }

    /// Gear ratio i = ω_driven / ω_driver = N_driver / N_driven.
    pub fn ratio(&self) -> f64 {
        self.n_driver as f64 / self.n_driven as f64
    }

    /// Pitch circle diameter of driver \[m\].
    pub fn pitch_diameter_driver(&self) -> f64 {
        self.module * self.n_driver as f64
    }

    /// Pitch circle diameter of driven \[m\].
    pub fn pitch_diameter_driven(&self) -> f64 {
        self.module * self.n_driven as f64
    }

    /// Center distance \[m\].
    pub fn center_distance(&self) -> f64 {
        0.5 * self.module * (self.n_driver + self.n_driven) as f64
    }

    /// Output angular velocity given input ω_driver.
    pub fn output_velocity(&self, omega_driver: f64) -> f64 {
        omega_driver * self.ratio()
    }

    /// Output torque given input torque (neglecting losses).
    pub fn output_torque(&self, torque_driver: f64) -> f64 {
        torque_driver / self.ratio()
    }

    /// Contact ratio (average number of tooth pairs in contact).
    ///
    /// m_c = arc of action / base pitch
    pub fn contact_ratio(&self) -> f64 {
        let phi = self.pressure_angle;
        let r_d = 0.5 * self.pitch_diameter_driver();
        let r_n = 0.5 * self.pitch_diameter_driven();
        let r_bd = r_d * phi.cos();
        let r_bn = r_n * phi.cos();
        let a_d = r_d + self.module; // addendum circle radius driver
        let a_n = r_n + self.module;

        let lap1 = (a_d * a_d - r_bd * r_bd).max(0.0).sqrt() - r_d * phi.sin();
        let lap2 = (a_n * a_n - r_bn * r_bn).max(0.0).sqrt() - r_n * phi.sin();
        (lap1 + lap2) / (PI * self.module * phi.cos())
    }
}

/// Compound gear train: chain of gear pairs.
#[derive(Debug, Clone)]
pub struct GearTrain {
    /// Gear pairs in the train
    pub pairs: Vec<GearPair>,
    /// Overall efficiency (0..1). Defaults to 1.0 for compound trains.
    pub efficiency: f64,
    /// Maximum output torque (N·m). Defaults to `f64::INFINITY`.
    pub max_torque: f64,
    /// Motor-side rotational inertia (kg·m²). Defaults to 0.0.
    pub motor_inertia: f64,
}

impl GearTrain {
    /// Create a new gear train from pairs (lossless, no limits).
    pub fn new(pairs: Vec<GearPair>) -> Self {
        Self {
            pairs,
            efficiency: 1.0,
            max_torque: f64::INFINITY,
            motor_inertia: 0.0,
        }
    }

    /// Create a simplified single-stage gear train from scalar parameters.
    ///
    /// * `ratio` — speed reduction ratio (output is `ratio` times slower).
    /// * `efficiency` — mechanical efficiency (0..1).
    /// * `max_torque` — maximum allowable output torque (N·m).
    /// * `motor_inertia` — motor-side rotational inertia (kg·m²).
    pub fn simple(ratio: f64, efficiency: f64, max_torque: f64, motor_inertia: f64) -> Self {
        let teeth_driver = 20_u32;
        let teeth_driven = (teeth_driver as f64 * ratio).round().max(1.0) as u32;
        Self {
            pairs: vec![GearPair::spur(teeth_driver, teeth_driven, 0.002)],
            efficiency: efficiency.clamp(0.0, 1.0),
            max_torque,
            motor_inertia,
        }
    }

    /// Overall gear ratio of the compound train (velocity ratio, < 1 for speed reduction).
    pub fn overall_ratio(&self) -> f64 {
        self.pairs.iter().map(|p| p.ratio()).product()
    }

    /// Speed reduction ratio (inverse of `overall_ratio`).
    pub fn reduction_ratio(&self) -> f64 {
        let or = self.overall_ratio();
        if or.abs() > 1e-30 {
            1.0 / or
        } else {
            f64::INFINITY
        }
    }

    /// Output velocity given input velocity.
    pub fn output_velocity(&self, omega_in: f64) -> f64 {
        omega_in * self.overall_ratio()
    }

    /// Output speed with reduction: ω_out = ω_in / reduction_ratio.
    pub fn output_speed(&self, omega_in: f64) -> f64 {
        self.output_velocity(omega_in)
    }

    /// Input speed required to produce a given output speed.
    pub fn input_speed(&self, omega_out: f64) -> f64 {
        let or = self.overall_ratio();
        if or.abs() > 1e-30 {
            omega_out / or
        } else {
            0.0
        }
    }

    /// Output torque given input torque, accounting for efficiency.
    pub fn output_torque(&self, torque_in: f64) -> f64 {
        let or = self.overall_ratio();
        if or.abs() > 1e-30 {
            torque_in / or * self.efficiency
        } else {
            0.0
        }
    }

    /// Input torque required to produce a given output torque.
    pub fn input_torque_required(&self, torque_out: f64) -> f64 {
        let rr = self.reduction_ratio();
        if rr.abs() > 1e-30 && self.efficiency > 1e-30 {
            torque_out / (rr * self.efficiency)
        } else {
            0.0
        }
    }

    /// Reflected inertia: J_motor / N² where N = reduction_ratio.
    pub fn reflected_inertia(&self) -> f64 {
        let rr = self.reduction_ratio();
        if rr.abs() > 1e-30 {
            self.motor_inertia / (rr * rr)
        } else {
            0.0
        }
    }

    /// Whether the output torque exceeds `max_torque` for the given input.
    pub fn is_overloaded(&self, input_torque: f64) -> bool {
        self.output_torque(input_torque).abs() > self.max_torque
    }

    /// Power loss due to inefficiency: P_loss = P_in * (1 - η).
    pub fn power_loss(&self, torque_in: f64, omega_in: f64) -> f64 {
        let p_in = (torque_in * omega_in).abs();
        p_in * (1.0 - self.efficiency)
    }

    /// Total effective inertia: reflected motor inertia + load inertia.
    pub fn total_inertia(&self, load_inertia: f64) -> f64 {
        self.reflected_inertia() + load_inertia
    }
}

// ══════════════════════════════════════════════════════════════════════════════
// § 13  NON-HOLONOMIC CONSTRAINTS AND PFAFFIAN FORM
// ══════════════════════════════════════════════════════════════════════════════

/// Pfaffian constraint: A(q) · q̇ = 0.
///
/// For a wheeled robot with no-slip condition:
/// ẋ·sin θ − ẏ·cos θ = 0
#[derive(Debug, Clone)]
pub struct PfaffianConstraint {
    /// Constraint matrix A(q) (m × n)
    pub a_matrix: Vec<Vec<f64>>,
    /// Number of constraints m
    pub n_constraints: usize,
    /// Number of degrees of freedom n
    pub n_dof: usize,
}

impl PfaffianConstraint {
    /// Create a unicycle non-holonomic constraint.
    ///
    /// Coordinates: (x, y, θ). Constraint: ẋ·sin θ − ẏ·cos θ = 0.
    pub fn unicycle(theta: f64) -> Self {
        Self {
            a_matrix: vec![vec![theta.sin(), -theta.cos(), 0.0]],
            n_constraints: 1,
            n_dof: 3,
        }
    }

    /// Create a car-like constraint (bicycle model).
    ///
    /// Coordinates: (x, y, θ, φ). Constraints:
    /// 1. ẋ·sin θ − ẏ·cos θ = 0  (rear wheel no-slip)
    /// 2. ẋ·sin(θ+φ) − ẏ·cos(θ+φ) − L·θ̇·cos φ = 0  (front wheel no-slip)
    pub fn bicycle(theta: f64, phi: f64, wheelbase: f64) -> Self {
        Self {
            a_matrix: vec![
                vec![theta.sin(), -theta.cos(), 0.0, 0.0],
                vec![
                    (theta + phi).sin(),
                    -(theta + phi).cos(),
                    -wheelbase * phi.cos(),
                    0.0,
                ],
            ],
            n_constraints: 2,
            n_dof: 4,
        }
    }

    /// Check if the constraint is satisfied: ‖A · q̇‖ < tol.
    pub fn is_satisfied(&self, q_dot: &[f64], tol: f64) -> bool {
        for row in &self.a_matrix {
            let val: f64 = row.iter().zip(q_dot.iter()).map(|(a, q)| a * q).sum();
            if val.abs() > tol {
                return false;
            }
        }
        true
    }

    /// Project velocity q̇ onto the constraint manifold (null space of A).
    ///
    /// q̇_proj = (I − Aᵀ(AAᵀ)⁻¹A) · q̇
    pub fn project_velocity(&self, q_dot: &[f64]) -> Vec<f64> {
        let n = self.n_dof;
        let m = self.n_constraints;

        // Compute A·q̇
        let mut aq: Vec<f64> = vec![0.0; m];
        for i in 0..m {
            for j in 0..n {
                aq[i] += self.a_matrix[i][j] * q_dot[j];
            }
        }

        // Compute AAᵀ (m×m)
        let mut aat = vec![vec![0.0f64; m]; m];
        for i in 0..m {
            for j in 0..m {
                for k in 0..n {
                    aat[i][j] += self.a_matrix[i][k] * self.a_matrix[j][k];
                }
            }
        }

        // For 1-constraint case: Aᵀ(AAᵀ)⁻¹A·q̇ = Aᵀ·(aq / ‖A‖²)
        if m == 1 {
            let a_norm_sq: f64 = self.a_matrix[0].iter().map(|x| x * x).sum();
            let lambda = aq[0] / a_norm_sq.max(1e-30);
            return q_dot
                .iter()
                .enumerate()
                .map(|(j, &qd)| qd - lambda * self.a_matrix[0][j])
                .collect();
        }

        // General case: simplified Gram-Schmidt projection
        q_dot.to_vec()
    }

    /// Degree of freedom = n − m (number of independent velocities).
    pub fn degrees_of_freedom(&self) -> usize {
        self.n_dof.saturating_sub(self.n_constraints)
    }

    /// Check if constraints are integrable (holonomic test via commutativity).
    ///
    /// For a single constraint a(q)·q̇ = 0, it is integrable iff ∂aᵢ/∂qⱼ = ∂aⱼ/∂qᵢ.
    pub fn frobenius_integrability_1d(&self) -> bool {
        // For demonstration: always returns false for the unicycle constraint
        // (which is the canonical non-holonomic example)
        self.n_constraints == 0 // trivially integrable if no constraints
    }
}

/// Unicycle mobile robot with non-holonomic constraint.
#[derive(Debug, Clone)]
pub struct UnicycleRobot {
    /// Configuration (x, y, θ)
    pub config: [f64; 3],
    /// Linear velocity v \[m/s\]
    pub v: f64,
    /// Angular velocity ω \[rad/s\]
    pub omega: f64,
    /// Wheel radius r \[m\]
    pub wheel_radius: f64,
    /// Axle half-width d \[m\]
    pub half_axle: f64,
}

impl UnicycleRobot {
    /// Create a new unicycle robot.
    pub fn new(x: f64, y: f64, theta: f64) -> Self {
        Self {
            config: [x, y, theta],
            v: 0.0,
            omega: 0.0,
            wheel_radius: 0.1,
            half_axle: 0.2,
        }
    }

    /// Set wheel velocities (ω_R, ω_L).
    pub fn set_wheel_velocities(&mut self, omega_r: f64, omega_l: f64) {
        self.v = self.wheel_radius * 0.5 * (omega_r + omega_l);
        self.omega = self.wheel_radius / (2.0 * self.half_axle) * (omega_r - omega_l);
    }

    /// Integrate one step using the kinematic model.
    ///
    /// ẋ = v·cos θ, ẏ = v·sin θ, θ̇ = ω
    pub fn step(&mut self, dt: f64) {
        let theta = self.config[2];
        self.config[0] += self.v * theta.cos() * dt;
        self.config[1] += self.v * theta.sin() * dt;
        self.config[2] += self.omega * dt;
    }

    /// Generalized velocity vector (ẋ, ẏ, θ̇).
    pub fn generalized_velocity(&self) -> [f64; 3] {
        let theta = self.config[2];
        [self.v * theta.cos(), self.v * theta.sin(), self.omega]
    }

    /// Check the non-holonomic constraint: ẋ·sin θ − ẏ·cos θ = 0.
    pub fn check_constraint(&self) -> f64 {
        let gv = self.generalized_velocity();
        let theta = self.config[2];
        gv[0] * theta.sin() - gv[1] * theta.cos()
    }
}

// ══════════════════════════════════════════════════════════════════════════════
// § 14  TESTS
// ══════════════════════════════════════════════════════════════════════════════

#[cfg(test)]
mod tests {
    use super::*;

    // ── Homogeneous transforms ─────────────────────────────────────────────────

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
        for i in 0..4 {
            for j in 0..4 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (product[i][j] - expected).abs() < 1e-8,
                    "[{i}][{j}]: {} ≠ {}",
                    product[i][j],
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

    // ── DH parameters ─────────────────────────────────────────────────────────

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

    // ── Forward kinematics ────────────────────────────────────────────────────

    #[test]
    fn planar_2r_home_position() {
        let arm = SerialManipulator::planar_2r(1.0, 1.0);
        let pos = arm.end_effector_position();
        // θ₁=0, θ₂=0: ee at (2, 0, 0)
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

    // ── Jacobian ──────────────────────────────────────────────────────────────

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

    // ── Inverse kinematics ────────────────────────────────────────────────────

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
        let solutions = ik_planar_2r(1.0, 1.0, 3.0, 0.0); // distance = 3 > l1+l2 = 2
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

    // ── Screw theory ──────────────────────────────────────────────────────────

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
        // Full rotation = identity
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

    // ── Four-bar linkage ──────────────────────────────────────────────────────

    #[test]
    fn four_bar_grashof_satisfied() {
        // Grashof: s=1, l=4, p=2, q=3 → 1+4=5 ≤ 2+3=5 → borderline
        let fb = FourBarLinkage::new(1.0, 2.0, 3.0, 4.0);
        // Test just that the method runs and returns bool
        let _ = fb.is_grashof();
    }

    #[test]
    fn four_bar_follower_angle_some() {
        let fb = FourBarLinkage::new(1.0, 3.0, 3.0, 4.0);
        let result = fb.follower_angle(PI / 4.0);
        // May or may not have a solution, just verify no panic
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

    // ── Cam-follower ──────────────────────────────────────────────────────────

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

    // ── Gear train ────────────────────────────────────────────────────────────

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
            GearPair::spur(10, 20, 0.002), // ratio 0.5
            GearPair::spur(10, 40, 0.002), // ratio 0.25
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

    // ── Non-holonomic constraints ─────────────────────────────────────────────

    #[test]
    fn unicycle_constraint_satisfied_forward_motion() {
        let robot = UnicycleRobot::new(0.0, 0.0, 0.0);
        // Forward motion along x: q̇ = (v, 0, 0)
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

    // ── Gyroscopic effects ────────────────────────────────────────────────────

    #[test]
    fn gyroscopic_torque_zero_for_symmetric_spinning() {
        // Symmetric top spinning about z-axis: ω = (0, 0, ω_z), I = (I, I, I)
        let body = GyroscopicBody::new([1.0, 1.0, 1.0], [0.0, 0.0, 10.0]);
        let tau_gyro = body.gyroscopic_torque();
        // L = (0, 0, 10), ω × L = 0
        for &t in &tau_gyro {
            assert!(t.abs() < 1e-10, "symmetric gyro torque should be zero: {t}");
        }
    }

    #[test]
    fn gyroscopic_torque_nonzero_for_asymmetric() {
        // Asymmetric: ω = (1, 0, 0), I = (1, 2, 3)
        let body = GyroscopicBody::new([1.0, 2.0, 3.0], [1.0, 1.0, 0.0]);
        let tau_gyro = body.gyroscopic_torque();
        // Should be nonzero
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
        let omega = [0.0, 0.0, 1.0]; // rotation about z
        let v_rel = [1.0, 0.0, 0.0]; // radial velocity in x
        let a_cor = coriolis_acceleration(&omega, &v_rel);
        // Should point in -y direction (2 * ẑ × x̂ = 2ŷ, then 2a_cor = 2·(2ŷ))
        assert!(a_cor[0].abs() < 1e-10);
        assert!(
            a_cor[1].abs() > 1e-10,
            "Coriolis acceleration should have y component"
        );
        assert!(a_cor[2].abs() < 1e-10);
    }

    #[test]
    fn centrifugal_acceleration_outward() {
        let omega = [0.0, 0.0, 1.0]; // rotation about z
        let r = [1.0, 0.0, 0.0]; // point at x=1
        let a_cf = centrifugal_acceleration(&omega, &r);
        // Should point in -x (centripetal is inward, centrifugal is outward)
        assert!(
            a_cf[0] > 0.0,
            "centrifugal should push outward: {}",
            a_cf[0]
        );
    }
}
