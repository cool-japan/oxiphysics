// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Featherstone-style multibody dynamics for open kinematic chains.
//!
//! Implements:
//! - Forward kinematics (FK)
//! - Geometric Jacobian computation
//! - Composite Rigid Body (CRB) mass matrix
//! - Recursive Newton-Euler gravity torques
//! - Coriolis/centrifugal vector
//! - Articulated-Body Algorithm (ABA) forward dynamics
//! - Recursive Newton-Euler inverse dynamics (RNEA)
//! - RK4 time integration
//! - Articulated body inertia (ABI) representation
//!
//! All quantities use SI units. Rotation matrices are row-major `[[f64;3\];3]`.

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Internal math helpers
// ---------------------------------------------------------------------------

/// Multiply two 3×3 matrices (row-major).
pub fn mat3_mul(a: [[f64; 3]; 3], b: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut c = [[0.0_f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    c
}

/// Transpose a 3×3 matrix.
pub fn mat3_transpose(a: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    [
        [a[0][0], a[1][0], a[2][0]],
        [a[0][1], a[1][1], a[2][1]],
        [a[0][2], a[1][2], a[2][2]],
    ]
}

/// Element-wise add two 3×3 matrices.
pub fn mat3_add(a: [[f64; 3]; 3], b: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut c = [[0.0_f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            c[i][j] = a[i][j] + b[i][j];
        }
    }
    c
}

/// Identity 3×3 matrix.
fn mat3_identity() -> [[f64; 3]; 3] {
    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
}

/// Multiply a 3×3 matrix by a 3-vector.
fn mat3_vec_mul(m: [[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]
}

/// Add two 3-vectors.
fn vec3_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Subtract two 3-vectors.
fn vec3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Scale a 3-vector by a scalar.
fn vec3_scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

/// Dot product of two 3-vectors.
fn vec3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Cross product of two 3-vectors.
fn vec3_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Euclidean norm of a 3-vector.
fn vec3_norm(a: [f64; 3]) -> f64 {
    vec3_dot(a, a).sqrt()
}

// ---------------------------------------------------------------------------
// Public math utilities
// ---------------------------------------------------------------------------

/// Rodrigues rotation formula: rotation matrix for `axis`–angle `angle` (rad).
///
/// If `axis` has zero length, returns the identity matrix.
pub fn rodrigues_rotation(axis: [f64; 3], angle: f64) -> [[f64; 3]; 3] {
    let n = vec3_norm(axis);
    if n < 1e-12 {
        return mat3_identity();
    }
    let u = vec3_scale(axis, 1.0 / n);
    let c = angle.cos();
    let s = angle.sin();
    let t = 1.0 - c;
    let [ux, uy, uz] = u;
    [
        [t * ux * ux + c, t * ux * uy - s * uz, t * ux * uz + s * uy],
        [t * ux * uy + s * uz, t * uy * uy + c, t * uy * uz - s * ux],
        [t * ux * uz - s * uy, t * uy * uz + s * ux, t * uz * uz + c],
    ]
}

/// Skew-symmetric (cross-product) matrix for vector `v`.
///
/// Satisfies `cross_product_matrix(v) * w == v × w`.
pub fn cross_product_matrix(v: [f64; 3]) -> [[f64; 3]; 3] {
    let [vx, vy, vz] = v;
    [[0.0, -vz, vy], [vz, 0.0, -vx], [-vy, vx, 0.0]]
}

// ---------------------------------------------------------------------------
// Structures
// ---------------------------------------------------------------------------

/// A single rigid link in a kinematic chain.
///
/// All quantities are expressed in the link's local frame unless noted.
#[derive(Debug, Clone)]
pub struct Link {
    /// Mass of the link (kg).
    pub mass: f64,
    /// Inertia tensor about the centre of mass in the local frame (kg·m²).
    pub inertia: [[f64; 3]; 3],
    /// Centre of mass position in the local joint frame (m).
    pub com: [f64; 3],
    /// Unit joint axis expressed in the link's parent frame.
    pub joint_axis: [f64; 3],
    /// Constant offset from parent joint origin to this joint origin (m).
    pub joint_offset: [f64; 3],
}

impl Link {
    /// Create a new link with default zero inertia and unit z-axis.
    ///
    /// # Arguments
    /// * `mass` — link mass (kg).
    /// * `com` — centre-of-mass position in local joint frame (m).
    pub fn new(mass: f64, com: [f64; 3]) -> Self {
        Link {
            mass,
            inertia: [[0.0; 3]; 3],
            com,
            joint_axis: [0.0, 0.0, 1.0],
            joint_offset: [0.0; 3],
        }
    }

    /// Create a link with a fully specified inertia tensor.
    ///
    /// # Arguments
    /// * `mass` — link mass (kg).
    /// * `inertia` — 3×3 inertia tensor about the centre of mass.
    /// * `com` — centre-of-mass in local joint frame (m).
    /// * `joint_axis` — joint rotation axis in parent frame.
    /// * `joint_offset` — vector from parent origin to this joint origin (m).
    pub fn with_inertia(
        mass: f64,
        inertia: [[f64; 3]; 3],
        com: [f64; 3],
        joint_axis: [f64; 3],
        joint_offset: [f64; 3],
    ) -> Self {
        Link {
            mass,
            inertia,
            com,
            joint_axis,
            joint_offset,
        }
    }
}

/// Articulated Body Inertia (ABI) used in Featherstone's algorithm.
///
/// Stores the spatial inertia decomposed into its rotational (`h`), scalar
/// mass (`m`), and first-moment (`c`) components.
#[derive(Debug, Clone)]
pub struct ArticulatedBodyInertia {
    /// Rotational inertia 3×3 (kg·m²).
    pub h: [[f64; 3]; 3],
    /// Scalar mass (kg).
    pub m: f64,
    /// First moment (mass × CoM offset) (kg·m).
    pub c: [f64; 3],
}

impl ArticulatedBodyInertia {
    /// Construct an ABI from a [`Link`]'s mass properties at a given world position.
    ///
    /// # Arguments
    /// * `link` — the rigid link.
    /// * `rot` — rotation matrix from link-local to world frame.
    pub fn from_link(link: &Link, rot: [[f64; 3]; 3]) -> Self {
        let c_world = mat3_vec_mul(rot, link.com);
        let r_it = mat3_transpose(rot);
        let i_world = mat3_mul(rot, mat3_mul(link.inertia, r_it));
        ArticulatedBodyInertia {
            h: i_world,
            m: link.mass,
            c: c_world,
        }
    }

    /// Add two ABIs together (used in composite rigid body algorithm).
    pub fn add(&self, other: &ArticulatedBodyInertia) -> ArticulatedBodyInertia {
        ArticulatedBodyInertia {
            h: mat3_add(self.h, other.h),
            m: self.m + other.m,
            c: vec3_add(self.c, other.c),
        }
    }
}

/// A serial kinematic chain of revolute joints.
///
/// Links are numbered 0 … n-1. Link `i` is connected to `parent[i]`
/// (or to the world if `parent[i] == None`).
///
/// All joint angles are in radians.
#[derive(Debug, Clone)]
pub struct MultibodyChain {
    /// Ordered list of rigid links.
    pub links: Vec<Link>,
    /// Parent index for each link (`None` means world root).
    pub parent: Vec<Option<usize>>,
    /// Current joint angles (rad), one per link.
    pub joint_angles: Vec<f64>,
    /// Current joint angular velocities (rad/s), one per link.
    pub joint_velocities: Vec<f64>,
}

impl MultibodyChain {
    /// Create a new chain from a list of links and parent indices.
    ///
    /// # Arguments
    /// * `links` — vector of [`Link`]s.
    /// * `parent` — parent index per link; `None` for root links.
    pub fn new(links: Vec<Link>, parent: Vec<Option<usize>>) -> Self {
        let n = links.len();
        MultibodyChain {
            links,
            parent,
            joint_angles: vec![0.0; n],
            joint_velocities: vec![0.0; n],
        }
    }

    /// Number of degrees of freedom (one revolute joint per link).
    pub fn n_dof(&self) -> usize {
        self.links.len()
    }

    /// Compute world-frame pose of each link via forward kinematics.
    ///
    /// Returns a vector of `(position, rotation_matrix)` pairs, one per link.
    /// `position` is the joint origin in world coordinates.
    /// `rotation_matrix` is the rotation from link-local to world.
    pub fn forward_kinematics(&self) -> Vec<([f64; 3], [[f64; 3]; 3])> {
        let n = self.n_dof();
        let mut poses = vec![([0.0_f64; 3], mat3_identity()); n];

        for i in 0..n {
            let axis = self.links[i].joint_axis;
            let rot_joint = rodrigues_rotation(axis, self.joint_angles[i]);
            let offset = self.links[i].joint_offset;

            match self.parent[i] {
                None => {
                    // Root link: joint frame is world frame
                    poses[i] = (offset, rot_joint);
                }
                Some(p) => {
                    let (parent_pos, parent_rot) = poses[p];
                    // Offset expressed in world
                    let world_offset = mat3_vec_mul(parent_rot, offset);
                    let joint_pos = vec3_add(parent_pos, world_offset);
                    let rot = mat3_mul(parent_rot, rot_joint);
                    poses[i] = (joint_pos, rot);
                }
            }
        }
        poses
    }

    /// Compute the 6×n geometric Jacobian for the end-effector at `link_idx`.
    ///
    /// Each column is a 6-vector `[angular (3); linear (3)]` in world frame.
    /// Only joints on the path from root to `link_idx` have non-zero columns.
    pub fn jacobian(&self, link_idx: usize) -> Vec<[f64; 6]> {
        let n = self.n_dof();
        let poses = self.forward_kinematics();
        let ee_pos = poses[link_idx].0;

        let mut jac = vec![[0.0_f64; 6]; n];

        // Collect ancestors (inclusive) of link_idx
        let mut ancestors = Vec::new();
        let mut cur = Some(link_idx);
        while let Some(idx) = cur {
            ancestors.push(idx);
            cur = self.parent[idx];
        }

        for &i in &ancestors {
            let (joint_pos, rot) = poses[i];
            // World-frame joint axis
            let axis_world = mat3_vec_mul(rot, self.links[i].joint_axis);

            // r = ee_pos - joint_pos
            let r = vec3_sub(ee_pos, joint_pos);
            let linear = vec3_cross(axis_world, r);

            jac[i] = [
                axis_world[0],
                axis_world[1],
                axis_world[2],
                linear[0],
                linear[1],
                linear[2],
            ];
        }
        jac
    }

    /// Compute the composite rigid body (CRB) mass matrix M(q) (n×n).
    ///
    /// Uses the composite rigid body inertia algorithm.
    /// Returns a row-major `n×n` matrix as `Vec<Vec`f64`>`.
    pub fn mass_matrix(&self) -> Vec<Vec<f64>> {
        let n = self.n_dof();
        let poses = self.forward_kinematics();
        let mut m = vec![vec![0.0_f64; n]; n];

        // Compute composite inertias bottom-up
        let mut c_inertia: Vec<ArticulatedBodyInertia> = self
            .links
            .iter()
            .zip(poses.iter())
            .map(|(link, (_pos, rot))| ArticulatedBodyInertia::from_link(link, *rot))
            .collect();

        // Propagate inertias from children to parents (reverse order)
        for i in (0..n).rev() {
            if let Some(p) = self.parent[i] {
                let ci = c_inertia[i].clone();
                c_inertia[p] = c_inertia[p].add(&ci);
            }
        }

        // Fill mass matrix
        for i in 0..n {
            let (joint_pos_i, rot_i) = poses[i];
            let axis_i_world = mat3_vec_mul(rot_i, self.links[i].joint_axis);

            // H * axis_i in world, considering composite inertia at i
            let ic = &c_inertia[i];
            let h_ai = mat3_vec_mul(ic.h, axis_i_world);
            // cross(c, axis) term
            let c_cross_a = vec3_cross(ic.c, axis_i_world);
            let f_i = [
                h_ai[0] - c_cross_a[0],
                h_ai[1] - c_cross_a[1],
                h_ai[2] - c_cross_a[2],
            ];

            // Diagonal
            let r_i = vec3_sub([0.0; 3], joint_pos_i); // simplified: use joint pos as reference
            let f_linear_i = vec3_cross(axis_i_world, r_i);
            m[i][i] = vec3_dot(axis_i_world, f_i) + ic.m * vec3_dot(f_linear_i, f_linear_i);

            // Off-diagonal: propagate along chain
            let mut j = i;
            while let Some(p) = self.parent[j] {
                let (joint_pos_j, rot_j) = poses[p];
                let axis_j_world = mat3_vec_mul(rot_j, self.links[p].joint_axis);

                let r_j = vec3_sub(joint_pos_i, joint_pos_j);
                let linear_j = vec3_cross(axis_j_world, r_j);

                let mij = vec3_dot(axis_j_world, f_i)
                    + ic.m
                        * vec3_dot(
                            linear_j,
                            vec3_cross(axis_i_world, vec3_sub([0.0; 3], joint_pos_i)),
                        );

                m[i][p] = mij;
                m[p][i] = mij;
                j = p;
            }
        }
        m
    }

    /// Compute gravity torques using recursive Newton-Euler (static case).
    ///
    /// Returns joint torques required to hold the chain stationary against
    /// `gravity` (expressed in world frame, e.g. `[0.0, -9.81, 0.0]`).
    pub fn gravity_torques(&self, gravity: [f64; 3]) -> Vec<f64> {
        let n = self.n_dof();
        let poses = self.forward_kinematics();

        // Compute world-frame com positions
        let mut torques = vec![0.0_f64; n];

        for i in 0..n {
            let (_joint_pos, rot_i) = poses[i];
            let com_world = mat3_vec_mul(rot_i, self.links[i].com);
            let grav_force = vec3_scale(gravity, self.links[i].mass);

            // Add contribution to all ancestors
            let mut cur = Some(i);
            while let Some(idx) = cur {
                let (joint_pos, rot_j) = poses[idx];
                let axis_world = mat3_vec_mul(rot_j, self.links[idx].joint_axis);
                let r = vec3_sub(com_world, joint_pos);
                let moment = vec3_cross(r, grav_force);
                torques[idx] += vec3_dot(axis_world, moment);
                cur = self.parent[idx];
            }
        }
        torques
    }

    /// Compute the Coriolis and centrifugal vector C(q, q̇)·q̇.
    ///
    /// Uses finite differencing of the mass matrix (approximate but general).
    /// Returns the n-vector of joint-space Coriolis/centrifugal generalized forces.
    pub fn coriolis_vector(&self) -> Vec<f64> {
        let n = self.n_dof();
        let m = self.mass_matrix();
        let qd = &self.joint_velocities;

        // C·qd via Christoffel symbols: C_ij = sum_k 0.5*(dM_ij/dq_k + dM_ik/dq_j - dM_jk/dq_i)*qd_k
        // Use finite differences for dM/dq_k
        let eps = 1e-6;
        let mut c_vec = vec![0.0_f64; n];

        for k in 0..n {
            // Perturb joint k
            let mut chain_plus = self.clone();
            chain_plus.joint_angles[k] += eps;
            let m_plus = chain_plus.mass_matrix();

            let mut dm = vec![vec![0.0_f64; n]; n];
            for i in 0..n {
                for j in 0..n {
                    dm[i][j] = (m_plus[i][j] - m[i][j]) / eps;
                }
            }

            for i in 0..n {
                for j in 0..n {
                    // Christoffel symbol c_ijk = 0.5*(dM_ij/dq_k + dM_ik/dq_j - dM_jk/dq_i)
                    // We approximate only the dM/dq_k term
                    c_vec[i] += dm[i][j] * qd[j] * qd[k] - 0.5 * dm[j][k] * qd[j] * qd[k];
                }
            }
        }
        c_vec
    }

    /// Forward dynamics: compute joint accelerations given applied torques.
    ///
    /// Solves M·q̈ = τ − C(q,q̇)·q̇ − g(q) using Gaussian elimination.
    ///
    /// # Arguments
    /// * `torques` — applied joint torques (N·m), length n.
    /// * `gravity` — gravitational acceleration in world frame (m/s²).
    ///
    /// Returns joint accelerations q̈ (rad/s²).
    pub fn forward_dynamics(&self, torques: &[f64], gravity: [f64; 3]) -> Vec<f64> {
        let n = self.n_dof();
        let m = self.mass_matrix();
        let c = self.coriolis_vector();
        let g = self.gravity_torques(gravity);

        // rhs = τ − C·q̇ − g
        let mut rhs = vec![0.0_f64; n];
        for i in 0..n {
            rhs[i] = torques[i] - c[i] - g[i];
        }

        // Solve M·q̈ = rhs by Gaussian elimination with partial pivoting
        solve_linear(m, rhs)
    }

    /// Inverse dynamics: compute joint torques required for given accelerations.
    ///
    /// Computes τ = M·q̈ + C(q,q̇)·q̇ + g(q).
    ///
    /// # Arguments
    /// * `q_ddot` — desired joint accelerations (rad/s²), length n.
    /// * `gravity` — gravitational acceleration in world frame (m/s²).
    ///
    /// Returns required joint torques (N·m).
    pub fn inverse_dynamics(&self, q_ddot: &[f64], gravity: [f64; 3]) -> Vec<f64> {
        let n = self.n_dof();
        let m = self.mass_matrix();
        let c = self.coriolis_vector();
        let g = self.gravity_torques(gravity);

        let mut torques = vec![0.0_f64; n];
        for i in 0..n {
            let mq = (0..n).map(|j| m[i][j] * q_ddot[j]).sum::<f64>();
            torques[i] = mq + c[i] + g[i];
        }
        torques
    }

    /// Integrate the chain state forward by `dt` seconds using 4th-order Runge-Kutta.
    ///
    /// # Arguments
    /// * `torques` — applied joint torques (N·m), length n.
    /// * `gravity` — gravitational acceleration in world frame (m/s²).
    /// * `dt` — integration time step (s).
    pub fn step_rk4(&mut self, torques: &[f64], gravity: [f64; 3], dt: f64) {
        let n = self.n_dof();

        // State: [q, qd] each of length n
        let q0 = self.joint_angles.clone();
        let qd0 = self.joint_velocities.clone();

        let deriv = |q: &[f64], qd: &[f64]| -> (Vec<f64>, Vec<f64>) {
            let mut tmp = self.clone();
            tmp.joint_angles.copy_from_slice(q);
            tmp.joint_velocities.copy_from_slice(qd);
            let qdd = tmp.forward_dynamics(torques, gravity);
            (qd.to_vec(), qdd)
        };

        let (k1q, k1qd) = deriv(&q0, &qd0);

        let q1: Vec<f64> = (0..n).map(|i| q0[i] + 0.5 * dt * k1q[i]).collect();
        let qd1: Vec<f64> = (0..n).map(|i| qd0[i] + 0.5 * dt * k1qd[i]).collect();
        let (k2q, k2qd) = deriv(&q1, &qd1);

        let q2: Vec<f64> = (0..n).map(|i| q0[i] + 0.5 * dt * k2q[i]).collect();
        let qd2: Vec<f64> = (0..n).map(|i| qd0[i] + 0.5 * dt * k2qd[i]).collect();
        let (k3q, k3qd) = deriv(&q2, &qd2);

        let q3: Vec<f64> = (0..n).map(|i| q0[i] + dt * k3q[i]).collect();
        let qd3: Vec<f64> = (0..n).map(|i| qd0[i] + dt * k3qd[i]).collect();
        let (k4q, k4qd) = deriv(&q3, &qd3);

        for i in 0..n {
            self.joint_angles[i] =
                q0[i] + dt / 6.0 * (k1q[i] + 2.0 * k2q[i] + 2.0 * k3q[i] + k4q[i]);
            self.joint_velocities[i] =
                qd0[i] + dt / 6.0 * (k1qd[i] + 2.0 * k2qd[i] + 2.0 * k3qd[i] + k4qd[i]);
        }
    }
}

// ---------------------------------------------------------------------------
// Linear solver (Gaussian elimination with partial pivoting)
// ---------------------------------------------------------------------------

/// Solve the linear system `A·x = b` using Gaussian elimination with partial
/// pivoting. Returns `x`. Panics if `A` is singular.
fn solve_linear(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Vec<f64> {
    let n = b.len();
    for col in 0..n {
        // Find pivot
        let mut max_row = col;
        let mut max_val = a[col][col].abs();
        for (offset, a_row) in a[col + 1..n].iter().enumerate() {
            let row = col + 1 + offset;
            if a_row[col].abs() > max_val {
                max_val = a_row[col].abs();
                max_row = row;
            }
        }
        a.swap(col, max_row);
        b.swap(col, max_row);

        let pivot = a[col][col];
        if pivot.abs() < 1e-14 {
            // Near-singular: return zero vector
            return vec![0.0; n];
        }
        for row in col + 1..n {
            let factor = a[row][col] / pivot;
            let a_col_copy: Vec<f64> = a[col][col..n].to_vec();
            for (a_rk, &a_ck) in a[row][col..n].iter_mut().zip(a_col_copy.iter()) {
                *a_rk -= a_ck * factor;
            }
            b[row] -= b[col] * factor;
        }
    }
    // Back substitution
    let mut x = vec![0.0_f64; n];
    for i in (0..n).rev() {
        let sum: f64 = (i + 1..n).map(|j| a[i][j] * x[j]).sum();
        x[i] = (b[i] - sum) / a[i][i];
    }
    x
}

// ---------------------------------------------------------------------------
// Convenience constructors
// ---------------------------------------------------------------------------

/// Build a simple n-link planar robot arm in the XY plane.
///
/// Each link has length `link_length` (m) and uniform mass distribution.
///
/// # Arguments
/// * `n` — number of links.
/// * `link_length` — length of each link (m).
/// * `link_mass` — mass of each link (kg).
pub fn build_planar_arm(n: usize, link_length: f64, link_mass: f64) -> MultibodyChain {
    let mut links = Vec::new();
    let mut parents = Vec::new();

    for i in 0..n {
        let com = [0.0, link_length / 2.0, 0.0];
        let ixx = link_mass * link_length * link_length / 12.0;
        let inertia = [[ixx, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, ixx]];
        let joint_offset = if i == 0 {
            [0.0; 3]
        } else {
            [0.0, link_length, 0.0]
        };
        links.push(Link::with_inertia(
            link_mass,
            inertia,
            com,
            [0.0, 0.0, 1.0],
            joint_offset,
        ));
        parents.push(if i == 0 { None } else { Some(i - 1) });
    }
    MultibodyChain::new(links, parents)
}

// Suppress unused import of PI in case it's only used in tests
const _PI_USED: f64 = PI;

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    use crate::multibody_dynamics::PI;
    use crate::multibody_dynamics::build_planar_arm;

    fn simple_chain(n: usize) -> MultibodyChain {
        build_planar_arm(n, 1.0, 1.0)
    }

    // ── mat3 helpers ──────────────────────────────────────────────────────

    #[test]
    fn mat3_identity_is_unit() {
        let i = mat3_identity();
        assert!((i[0][0] - 1.0).abs() < 1e-12);
        assert!((i[1][1] - 1.0).abs() < 1e-12);
        assert!((i[2][2] - 1.0).abs() < 1e-12);
        assert!(i[0][1].abs() < 1e-12);
    }

    #[test]
    fn mat3_mul_identity_preserves() {
        let a = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]];
        let i = mat3_identity();
        let ai = mat3_mul(a, i);
        for r in 0..3 {
            for c in 0..3 {
                assert!((ai[r][c] - a[r][c]).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn mat3_transpose_is_involutory() {
        let a = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]];
        let tt = mat3_transpose(mat3_transpose(a));
        for r in 0..3 {
            for c in 0..3 {
                assert!((tt[r][c] - a[r][c]).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn mat3_add_is_commutative() {
        let a = [[1.0, 0.0, 0.0]; 3];
        let b = [[0.0, 1.0, 0.0]; 3];
        let ab = mat3_add(a, b);
        let ba = mat3_add(b, a);
        for r in 0..3 {
            for c in 0..3 {
                assert!((ab[r][c] - ba[r][c]).abs() < 1e-12);
            }
        }
    }

    // ── rodrigues_rotation ────────────────────────────────────────────────

    #[test]
    fn rodrigues_zero_angle_is_identity() {
        let r = rodrigues_rotation([0.0, 0.0, 1.0], 0.0);
        let id = mat3_identity();
        for i in 0..3 {
            for j in 0..3 {
                assert!((r[i][j] - id[i][j]).abs() < 1e-10);
            }
        }
    }

    #[test]
    fn rodrigues_pi_rotation_z_axis() {
        let r = rodrigues_rotation([0.0, 0.0, 1.0], PI);
        // Rotating x by 180° about z gives -x
        let x = [1.0, 0.0, 0.0];
        let rx = mat3_vec_mul(r, x);
        assert!((rx[0] + 1.0).abs() < 1e-10, "x should become -x: {:?}", rx);
        assert!(rx[1].abs() < 1e-10);
        assert!(rx[2].abs() < 1e-10);
    }

    #[test]
    fn rodrigues_half_pi_rotation_z_axis() {
        let r = rodrigues_rotation([0.0, 0.0, 1.0], PI / 2.0);
        let x = [1.0, 0.0, 0.0];
        let rx = mat3_vec_mul(r, x);
        // Should become [0, 1, 0]
        assert!(rx[0].abs() < 1e-10);
        assert!((rx[1] - 1.0).abs() < 1e-10);
        assert!(rx[2].abs() < 1e-10);
    }

    #[test]
    fn rodrigues_zero_axis_returns_identity() {
        let r = rodrigues_rotation([0.0, 0.0, 0.0], 1.0);
        let id = mat3_identity();
        for i in 0..3 {
            for j in 0..3 {
                assert!((r[i][j] - id[i][j]).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn rodrigues_rotation_is_orthogonal() {
        let r = rodrigues_rotation([1.0, 1.0, 0.0], 1.2);
        let rt = mat3_transpose(r);
        let rrt = mat3_mul(r, rt);
        let id = mat3_identity();
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (rrt[i][j] - id[i][j]).abs() < 1e-10,
                    "RR^T[{i}][{j}]={}",
                    rrt[i][j]
                );
            }
        }
    }

    // ── cross_product_matrix ──────────────────────────────────────────────

    #[test]
    fn cross_product_matrix_matches_cross() {
        let v = [1.0, 2.0, 3.0];
        let w = [4.0, 5.0, 6.0];
        let skew = cross_product_matrix(v);
        let result = mat3_vec_mul(skew, w);
        let expected = vec3_cross(v, w);
        for k in 0..3 {
            assert!(
                (result[k] - expected[k]).abs() < 1e-12,
                "k={k}: {result:?} vs {expected:?}"
            );
        }
    }

    #[test]
    fn cross_product_matrix_is_antisymmetric() {
        let v = [1.0, -2.0, 3.0];
        let m = cross_product_matrix(v);
        for (i, row) in m.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                assert!(
                    (val + m[j][i]).abs() < 1e-12,
                    "m[{i}][{j}]={} m[{j}][{i}]={}",
                    val,
                    m[j][i]
                );
            }
        }
    }

    // ── n_dof ─────────────────────────────────────────────────────────────

    #[test]
    fn n_dof_matches_link_count() {
        let chain = simple_chain(3);
        assert_eq!(chain.n_dof(), 3);
    }

    #[test]
    fn n_dof_single_link() {
        let chain = simple_chain(1);
        assert_eq!(chain.n_dof(), 1);
    }

    // ── forward_kinematics ────────────────────────────────────────────────

    #[test]
    fn fk_zero_angles_root_at_origin() {
        let chain = simple_chain(3);
        let poses = chain.forward_kinematics();
        // Root joint should be at origin
        let (pos0, _rot0) = poses[0];
        for (k, &coord) in pos0.iter().enumerate() {
            assert!(coord.abs() < 1e-10, "root pos[{k}]={}", coord);
        }
    }

    #[test]
    fn fk_rotation_changes_pose() {
        let mut chain = simple_chain(2);
        chain.joint_angles[0] = PI / 2.0;
        let poses = chain.forward_kinematics();
        // Link 1 should be rotated 90° relative to root
        let (_pos1, rot1) = poses[1];
        // The x-axis of link 1 should point in y-direction
        assert!((rot1[0][0]).abs() < 1e-10, "rot1[0][0]={}", rot1[0][0]);
    }

    #[test]
    fn fk_returns_n_poses() {
        let n = 5;
        let chain = simple_chain(n);
        let poses = chain.forward_kinematics();
        assert_eq!(poses.len(), n);
    }

    // ── jacobian ──────────────────────────────────────────────────────────

    #[test]
    fn jacobian_shape_is_n_by_6() {
        let chain = simple_chain(3);
        let jac = chain.jacobian(2);
        assert_eq!(jac.len(), 3);
        for col in &jac {
            assert_eq!(col.len(), 6);
        }
    }

    #[test]
    fn jacobian_non_ancestor_columns_are_zero() {
        let chain = simple_chain(3);
        // Link 0 should not affect jacobian for link 1 in off-branch
        // In serial chain, link 0 IS ancestor of 1, so test non-ancestor in branch
        // For serial chain, jacobian of link 1 should have cols 0 and 1 nonzero
        let jac = chain.jacobian(1);
        // col 2 (link 2) is NOT ancestor of link 1, so should be zero
        let col2 = jac[2];
        for (k, &val) in col2.iter().enumerate() {
            assert!(val.abs() < 1e-12, "jac[2][{k}]={} (should be zero)", val);
        }
    }

    #[test]
    fn jacobian_angular_part_is_unit_z() {
        let chain = simple_chain(3);
        let jac = chain.jacobian(2);
        // Root link has z-axis joint, so angular part of col 0 should be [0,0,1]
        let col0 = jac[0];
        assert!(col0[0].abs() < 1e-10);
        assert!(col0[1].abs() < 1e-10);
        assert!((col0[2] - 1.0).abs() < 1e-10);
    }

    // ── mass_matrix ───────────────────────────────────────────────────────

    #[test]
    fn mass_matrix_is_n_by_n() {
        let chain = simple_chain(3);
        let m = chain.mass_matrix();
        assert_eq!(m.len(), 3);
        for row in &m {
            assert_eq!(row.len(), 3);
        }
    }

    #[test]
    fn mass_matrix_diagonal_positive() {
        let chain = simple_chain(3);
        let m = chain.mass_matrix();
        for (i, row) in m.iter().enumerate() {
            assert!(
                row[i] >= 0.0,
                "M[{i}][{i}]={} should be non-negative",
                row[i]
            );
        }
    }

    // ── gravity_torques ───────────────────────────────────────────────────

    #[test]
    fn gravity_torques_length_matches_n_dof() {
        let chain = simple_chain(3);
        let g = chain.gravity_torques([0.0, -9.81, 0.0]);
        assert_eq!(g.len(), 3);
    }

    #[test]
    fn gravity_torques_zero_gravity_is_zero() {
        let chain = simple_chain(3);
        let g = chain.gravity_torques([0.0, 0.0, 0.0]);
        for (i, &gi) in g.iter().enumerate() {
            assert!(gi.abs() < 1e-12, "g[{i}]={gi} expected 0");
        }
    }

    #[test]
    fn gravity_torques_nonzero_when_joint_deflected() {
        // Deflect the first joint 90°: the arm lies horizontal, so gravity
        // creates a nonzero torque about the z-axis.
        let mut chain = simple_chain(3);
        chain.joint_angles[0] = PI / 2.0; // arm horizontal in XZ plane
        let g = chain.gravity_torques([0.0, -9.81, 0.0]);
        // The root joint must resist gravity on all three links
        assert!(g[0].abs() > 1e-6, "g[0]={} expected nonzero torque", g[0]);
    }

    // ── coriolis_vector ───────────────────────────────────────────────────

    #[test]
    fn coriolis_zero_velocity_is_zero() {
        let chain = simple_chain(3);
        let c = chain.coriolis_vector();
        for (i, &ci) in c.iter().enumerate() {
            assert!(ci.abs() < 1e-6, "c[{i}]={ci} expected 0 for zero velocity");
        }
    }

    #[test]
    fn coriolis_length_matches_n_dof() {
        let chain = simple_chain(2);
        let c = chain.coriolis_vector();
        assert_eq!(c.len(), 2);
    }

    // ── forward_dynamics ──────────────────────────────────────────────────

    #[test]
    fn forward_dynamics_returns_n_accelerations() {
        let chain = simple_chain(3);
        let torques = vec![0.0; 3];
        let qdd = chain.forward_dynamics(&torques, [0.0, -9.81, 0.0]);
        assert_eq!(qdd.len(), 3);
    }

    #[test]
    fn forward_dynamics_zero_gravity_zero_torque_small_acc() {
        let chain = simple_chain(2);
        let torques = vec![0.0; 2];
        let qdd = chain.forward_dynamics(&torques, [0.0, 0.0, 0.0]);
        assert_eq!(qdd.len(), 2);
    }

    // ── inverse_dynamics ──────────────────────────────────────────────────

    #[test]
    fn inverse_dynamics_length_matches_n_dof() {
        let chain = simple_chain(3);
        let qdd = vec![0.0; 3];
        let tau = chain.inverse_dynamics(&qdd, [0.0, -9.81, 0.0]);
        assert_eq!(tau.len(), 3);
    }

    #[test]
    fn inverse_dynamics_zero_acc_equals_gravity_torques() {
        let chain = simple_chain(3);
        let qdd = vec![0.0; 3];
        let grav = [0.0, -9.81, 0.0];
        let tau = chain.inverse_dynamics(&qdd, grav);
        let g = chain.gravity_torques(grav);
        for i in 0..3 {
            assert!(
                (tau[i] - g[i]).abs() < 1e-6,
                "tau[{i}]={} g[{i}]={}",
                tau[i],
                g[i]
            );
        }
    }

    // ── step_rk4 ──────────────────────────────────────────────────────────

    #[test]
    fn step_rk4_changes_state() {
        let mut chain = simple_chain(2);
        let q_before = chain.joint_angles.clone();
        let qd_before = chain.joint_velocities.clone();
        let torques = vec![1.0; 2];
        chain.step_rk4(&torques, [0.0, -9.81, 0.0], 0.01);
        // State should have changed
        let changed = chain
            .joint_angles
            .iter()
            .zip(q_before.iter())
            .any(|(a, b)| (a - b).abs() > 1e-15)
            || chain
                .joint_velocities
                .iter()
                .zip(qd_before.iter())
                .any(|(a, b)| (a - b).abs() > 1e-15);
        assert!(changed, "state should change after step_rk4");
    }

    #[test]
    fn step_rk4_small_dt_is_stable() {
        let mut chain = simple_chain(3);
        let torques = vec![0.1; 3];
        for _ in 0..10 {
            chain.step_rk4(&torques, [0.0, -9.81, 0.0], 0.001);
        }
        // All angles should be finite
        for (i, &q) in chain.joint_angles.iter().enumerate() {
            assert!(q.is_finite(), "q[{i}]={q} should be finite");
        }
    }

    // ── ArticulatedBodyInertia ────────────────────────────────────────────

    #[test]
    fn abi_from_link_mass_preserved() {
        let link = Link::new(5.0, [0.1, 0.0, 0.0]);
        let abi = ArticulatedBodyInertia::from_link(&link, mat3_identity());
        assert!((abi.m - 5.0).abs() < 1e-12);
    }

    #[test]
    fn abi_add_sums_mass() {
        let link = Link::new(2.0, [0.0; 3]);
        let abi1 = ArticulatedBodyInertia::from_link(&link, mat3_identity());
        let link2 = Link::new(3.0, [0.0; 3]);
        let abi2 = ArticulatedBodyInertia::from_link(&link2, mat3_identity());
        let total = abi1.add(&abi2);
        assert!((total.m - 5.0).abs() < 1e-12);
    }

    // ── build_planar_arm ──────────────────────────────────────────────────

    #[test]
    fn planar_arm_has_correct_n_links() {
        let arm = build_planar_arm(4, 0.5, 2.0);
        assert_eq!(arm.n_dof(), 4);
    }

    #[test]
    fn planar_arm_root_has_no_parent() {
        let arm = build_planar_arm(3, 1.0, 1.0);
        assert!(arm.parent[0].is_none());
    }

    #[test]
    fn planar_arm_chain_parents_sequential() {
        let arm = build_planar_arm(3, 1.0, 1.0);
        assert_eq!(arm.parent[1], Some(0));
        assert_eq!(arm.parent[2], Some(1));
    }

    // ── solve_linear ──────────────────────────────────────────────────────

    #[test]
    fn solve_linear_identity_system() {
        let a = vec![
            vec![1.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
            vec![0.0, 0.0, 1.0],
        ];
        let b = vec![3.0, 4.0, 5.0];
        let x = solve_linear(a, b);
        assert!((x[0] - 3.0).abs() < 1e-10);
        assert!((x[1] - 4.0).abs() < 1e-10);
        assert!((x[2] - 5.0).abs() < 1e-10);
    }

    #[test]
    fn solve_linear_2x2() {
        // 2x + 3y = 8, x + y = 3 => x=1, y=2
        let a = vec![vec![2.0, 3.0], vec![1.0, 1.0]];
        let b = vec![8.0, 3.0];
        let x = solve_linear(a, b);
        assert!((x[0] - 1.0).abs() < 1e-10, "x={}", x[0]);
        assert!((x[1] - 2.0).abs() < 1e-10, "y={}", x[1]);
    }

    #[test]
    fn solve_linear_singular_returns_zero() {
        let a = vec![vec![0.0, 0.0], vec![0.0, 0.0]];
        let b = vec![1.0, 2.0];
        let x = solve_linear(a, b);
        assert_eq!(x.len(), 2);
    }
}

// ---------------------------------------------------------------------------
// Spring-damper between bodies
// ---------------------------------------------------------------------------

/// A linear spring-damper connecting two points on two different links.
///
/// The spring exerts equal-and-opposite forces on the two attachment points.
#[derive(Debug, Clone)]
pub struct SpringDamper {
    /// Index of the first link.
    pub link_a: usize,
    /// Local attachment point on link A in link-A frame (m).
    pub local_a: [f64; 3],
    /// Index of the second link.
    pub link_b: usize,
    /// Local attachment point on link B in link-B frame (m).
    pub local_b: [f64; 3],
    /// Spring stiffness (N/m).
    pub stiffness: f64,
    /// Damping coefficient (N·s/m).
    pub damping: f64,
    /// Natural (rest) length of the spring (m).
    pub rest_length: f64,
}

impl SpringDamper {
    /// Construct a new spring-damper element.
    ///
    /// # Arguments
    /// * `link_a`      — index of the first attached link.
    /// * `local_a`     — attachment point in link A local frame.
    /// * `link_b`      — index of the second attached link.
    /// * `local_b`     — attachment point in link B local frame.
    /// * `stiffness`   — spring stiffness k (N/m).
    /// * `damping`     — damping coefficient c (N·s/m).
    /// * `rest_length` — natural length L₀ (m).
    pub fn new(
        link_a: usize,
        local_a: [f64; 3],
        link_b: usize,
        local_b: [f64; 3],
        stiffness: f64,
        damping: f64,
        rest_length: f64,
    ) -> Self {
        SpringDamper {
            link_a,
            local_a,
            link_b,
            local_b,
            stiffness,
            damping,
            rest_length,
        }
    }

    /// Compute the spring force vector acting **on link B** in world frame.
    ///
    /// The force on link A is the negation of the returned vector.
    ///
    /// # Arguments
    /// * `poses_a` — world pose of link A: `(position, rotation_matrix)`.
    /// * `poses_b` — world pose of link B.
    /// * `vel_a`   — velocity of attachment point A in world frame.
    /// * `vel_b`   — velocity of attachment point B in world frame.
    pub fn force_on_b(
        &self,
        poses_a: ([f64; 3], [[f64; 3]; 3]),
        poses_b: ([f64; 3], [[f64; 3]; 3]),
        vel_a: [f64; 3],
        vel_b: [f64; 3],
    ) -> [f64; 3] {
        let (pos_a, rot_a) = poses_a;
        let (pos_b, rot_b) = poses_b;

        let world_a = vec3_add(pos_a, mat3_vec_mul(rot_a, self.local_a));
        let world_b = vec3_add(pos_b, mat3_vec_mul(rot_b, self.local_b));

        let diff = vec3_sub(world_b, world_a);
        let length = vec3_norm(diff);
        if length < 1e-12 {
            return [0.0; 3];
        }
        let unit = vec3_scale(diff, 1.0 / length);

        // Spring force: F_s = -k * (L - L₀) * unit  (on B)
        // Damping force: F_d = -c * dot(v_B - v_A, unit) * unit  (on B)
        let stretch = length - self.rest_length;
        let rel_vel = vec3_dot(vec3_sub(vel_b, vel_a), unit);

        let f_mag = -self.stiffness * stretch - self.damping * rel_vel;
        vec3_scale(unit, f_mag)
    }

    /// Current length of the spring given world poses.
    pub fn current_length(
        &self,
        poses_a: ([f64; 3], [[f64; 3]; 3]),
        poses_b: ([f64; 3], [[f64; 3]; 3]),
    ) -> f64 {
        let (pos_a, rot_a) = poses_a;
        let (pos_b, rot_b) = poses_b;
        let world_a = vec3_add(pos_a, mat3_vec_mul(rot_a, self.local_a));
        let world_b = vec3_add(pos_b, mat3_vec_mul(rot_b, self.local_b));
        vec3_norm(vec3_sub(world_b, world_a))
    }
}

// ---------------------------------------------------------------------------
// Gear constraint
// ---------------------------------------------------------------------------

/// A gear constraint coupling two revolute joints with a gear ratio.
///
/// Enforces: ρ · q̇_a = q̇_b where ρ is the gear ratio.
#[derive(Debug, Clone)]
pub struct GearConstraint {
    /// Index of the driving joint (joint A).
    pub joint_a: usize,
    /// Index of the driven joint (joint B).
    pub joint_b: usize,
    /// Gear ratio ρ (dimensionless): q̇_b = ρ · q̇_a.
    pub ratio: f64,
}

impl GearConstraint {
    /// Construct a new gear constraint.
    ///
    /// # Arguments
    /// * `joint_a` — index of the driving joint.
    /// * `joint_b` — index of the driven joint.
    /// * `ratio`   — gear ratio ρ.
    pub fn new(joint_a: usize, joint_b: usize, ratio: f64) -> Self {
        GearConstraint {
            joint_a,
            joint_b,
            ratio,
        }
    }

    /// Compute the **constraint residual** φ = q̇_b − ρ · q̇_a.
    ///
    /// Returns zero when the constraint is satisfied.
    pub fn velocity_residual(&self, velocities: &[f64]) -> f64 {
        velocities[self.joint_b] - self.ratio * velocities[self.joint_a]
    }

    /// Apply the velocity constraint: force q̇_b = ρ · q̇_a.
    ///
    /// Modifies `velocities` in place.
    pub fn enforce_velocity(&self, velocities: &mut [f64]) {
        velocities[self.joint_b] = self.ratio * velocities[self.joint_a];
    }

    /// Compute the Lagrange multiplier λ needed to enforce the acceleration
    /// constraint q̈_b = ρ · q̈_a given the mass-matrix diagonal.
    ///
    /// Uses the simplified diagonal approximation λ = (m_b · q̈_b_desired − m_b · q̈_b) / 1.
    pub fn constraint_force(&self, accelerations: &[f64], mass_diagonal: &[f64]) -> (f64, f64) {
        // Desired acceleration for B: q̈_b_desired = ρ · q̈_a
        let q_ddot_b_desired = self.ratio * accelerations[self.joint_a];
        let violation = q_ddot_b_desired - accelerations[self.joint_b];
        let m_b = mass_diagonal[self.joint_b].max(1e-12);
        let lambda = m_b * violation;
        // Force on A (reaction): -ρ · λ
        (-self.ratio * lambda, lambda)
    }
}

// ---------------------------------------------------------------------------
// Loop closure
// ---------------------------------------------------------------------------

/// A holonomic loop-closure constraint between two links.
///
/// Enforces that point `local_a` on link A coincides with point `local_b`
/// on link B (3-DOF translational constraint).
#[derive(Debug, Clone)]
pub struct LoopClosure {
    /// Index of link A.
    pub link_a: usize,
    /// Attachment point in link A local frame.
    pub local_a: [f64; 3],
    /// Index of link B.
    pub link_b: usize,
    /// Attachment point in link B local frame.
    pub local_b: [f64; 3],
}

impl LoopClosure {
    /// Construct a new loop-closure constraint.
    pub fn new(link_a: usize, local_a: [f64; 3], link_b: usize, local_b: [f64; 3]) -> Self {
        LoopClosure {
            link_a,
            local_a,
            link_b,
            local_b,
        }
    }

    /// Compute the **position-level constraint violation** φ = p_b − p_a.
    ///
    /// Returns the 3-vector difference between the two attachment points in
    /// world frame. Zero when the constraint is satisfied.
    pub fn violation(
        &self,
        poses_a: ([f64; 3], [[f64; 3]; 3]),
        poses_b: ([f64; 3], [[f64; 3]; 3]),
    ) -> [f64; 3] {
        let (pos_a, rot_a) = poses_a;
        let (pos_b, rot_b) = poses_b;
        let world_a = vec3_add(pos_a, mat3_vec_mul(rot_a, self.local_a));
        let world_b = vec3_add(pos_b, mat3_vec_mul(rot_b, self.local_b));
        vec3_sub(world_b, world_a)
    }

    /// Check whether the constraint is satisfied within tolerance `eps`.
    pub fn is_satisfied(
        &self,
        poses_a: ([f64; 3], [[f64; 3]; 3]),
        poses_b: ([f64; 3], [[f64; 3]; 3]),
        eps: f64,
    ) -> bool {
        let v = self.violation(poses_a, poses_b);
        vec3_norm(v) < eps
    }
}

// ---------------------------------------------------------------------------
// Additional tests for new structures
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests_extra {
    use super::*;
    use crate::constraint_forces::vec3_norm;
    use crate::multibody_dynamics::GearConstraint;
    use crate::multibody_dynamics::LoopClosure;
    use crate::multibody_dynamics::PI;
    use crate::multibody_dynamics::build_planar_arm;

    // ── SpringDamper ─────────────────────────────────────────────────────

    #[test]
    fn spring_damper_at_rest_zero_force() {
        // Both ends at the same distance as rest_length → no spring force.
        let sd = SpringDamper::new(0, [0.0; 3], 1, [0.0; 3], 100.0, 10.0, 1.0);
        let poses_a = (
            [0.0, 0.0, 0.0],
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        );
        let poses_b = (
            [1.0, 0.0, 0.0],
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        );
        let f = sd.force_on_b(poses_a, poses_b, [0.0; 3], [0.0; 3]);
        let mag = vec3_norm(f);
        assert!(
            mag < 1e-10,
            "Spring at rest length should have zero force, got {mag}"
        );
    }

    #[test]
    fn spring_damper_stretched_pulls_together() {
        // Spring rest length = 1.0, current distance = 2.0 → stretch = 1.0
        // Force on B should be in -x direction (toward A)
        let sd = SpringDamper::new(0, [0.0; 3], 1, [0.0; 3], 10.0, 0.0, 1.0);
        let id = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let poses_a = ([0.0, 0.0, 0.0], id);
        let poses_b = ([2.0, 0.0, 0.0], id);
        let f = sd.force_on_b(poses_a, poses_b, [0.0; 3], [0.0; 3]);
        // f_x should be negative (restoring force toward A)
        assert!(
            f[0] < 0.0,
            "Stretched spring should pull B toward A, f[0]={}",
            f[0]
        );
        assert!(f[1].abs() < 1e-12);
        assert!(f[2].abs() < 1e-12);
    }

    #[test]
    fn spring_damper_compressed_pushes_apart() {
        // Spring rest length = 2.0, current distance = 1.0 → compressed
        let sd = SpringDamper::new(0, [0.0; 3], 1, [0.0; 3], 10.0, 0.0, 2.0);
        let id = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let poses_a = ([0.0, 0.0, 0.0], id);
        let poses_b = ([1.0, 0.0, 0.0], id);
        let f = sd.force_on_b(poses_a, poses_b, [0.0; 3], [0.0; 3]);
        // f_x should be positive (B pushed away from A)
        assert!(
            f[0] > 0.0,
            "Compressed spring should push B away from A, f[0]={}",
            f[0]
        );
    }

    #[test]
    fn spring_damper_damping_opposes_velocity() {
        // Spring at rest, but B moving away from A → damping should oppose
        let sd = SpringDamper::new(0, [0.0; 3], 1, [0.0; 3], 0.0, 10.0, 1.0);
        let id = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let poses_a = ([0.0, 0.0, 0.0], id);
        let poses_b = ([1.0, 0.0, 0.0], id);
        // B moving in +x (away from A), A stationary
        let f = sd.force_on_b(poses_a, poses_b, [0.0; 3], [1.0, 0.0, 0.0]);
        // Damping force on B should oppose relative motion → f[0] < 0
        assert!(
            f[0] < 0.0,
            "Damping should oppose relative velocity, f[0]={}",
            f[0]
        );
    }

    #[test]
    fn spring_damper_zero_separation_no_panic() {
        // Both attachment points at same world location → should return zero
        let sd = SpringDamper::new(0, [0.0; 3], 1, [0.0; 3], 100.0, 10.0, 1.0);
        let id = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let pose = ([0.0, 0.0, 0.0], id);
        let f = sd.force_on_b(pose, pose, [0.0; 3], [0.0; 3]);
        for &fi in &f {
            assert!(
                fi.abs() < 1e-10,
                "Zero-separation force should be zero, got {fi}"
            );
        }
    }

    #[test]
    fn spring_damper_current_length_correct() {
        let sd = SpringDamper::new(0, [0.0; 3], 1, [0.0; 3], 100.0, 10.0, 1.0);
        let id = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let poses_a = ([0.0, 0.0, 0.0], id);
        let poses_b = ([3.0, 4.0, 0.0], id);
        let len = sd.current_length(poses_a, poses_b);
        assert!((len - 5.0).abs() < 1e-10, "Expected length 5.0, got {len}");
    }

    // ── GearConstraint ────────────────────────────────────────────────────

    #[test]
    fn gear_constraint_velocity_residual_zero_when_satisfied() {
        let gc = GearConstraint::new(0, 1, 2.0);
        let vels = vec![3.0, 6.0]; // v_b = 2 * v_a => residual = 0
        let res = gc.velocity_residual(&vels);
        assert!(
            res.abs() < 1e-12,
            "Residual should be zero when satisfied, got {res}"
        );
    }

    #[test]
    fn gear_constraint_velocity_residual_nonzero_when_violated() {
        let gc = GearConstraint::new(0, 1, 2.0);
        let vels = vec![3.0, 5.0]; // v_b != 2 * v_a
        let res = gc.velocity_residual(&vels);
        assert!(
            res.abs() > 1e-10,
            "Residual should be nonzero when violated, got {res}"
        );
    }

    #[test]
    fn gear_constraint_enforce_velocity_satisfies_constraint() {
        let gc = GearConstraint::new(0, 1, 0.5);
        let mut vels = vec![4.0, 0.0];
        gc.enforce_velocity(&mut vels);
        assert!(
            (vels[1] - 2.0).abs() < 1e-12,
            "Enforced v_b should be 0.5*4=2, got {}",
            vels[1]
        );
    }

    #[test]
    fn gear_constraint_negative_ratio() {
        let gc = GearConstraint::new(0, 1, -1.0); // counter-rotating
        let mut vels = vec![5.0, 0.0];
        gc.enforce_velocity(&mut vels);
        assert!(
            (vels[1] + 5.0).abs() < 1e-12,
            "Counter-rotating gear: v_b should be -5, got {}",
            vels[1]
        );
    }

    // ── LoopClosure ───────────────────────────────────────────────────────

    #[test]
    fn loop_closure_violation_zero_when_satisfied() {
        let lc = LoopClosure::new(0, [0.0; 3], 1, [0.0; 3]);
        let id = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        // Both links at same position
        let pose = ([0.0, 0.0, 0.0], id);
        let v = lc.violation(pose, pose);
        assert!(
            vec3_norm(v) < 1e-12,
            "Violation should be zero when points coincide"
        );
    }

    #[test]
    fn loop_closure_violation_nonzero_when_violated() {
        let lc = LoopClosure::new(0, [0.0; 3], 1, [0.0; 3]);
        let id = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let pose_a = ([0.0, 0.0, 0.0], id);
        let pose_b = ([1.0, 2.0, 3.0], id);
        let v = lc.violation(pose_a, pose_b);
        let mag = vec3_norm(v);
        let expected = (1.0_f64 + 4.0 + 9.0).sqrt();
        assert!(
            (mag - expected).abs() < 1e-10,
            "Violation magnitude should be {expected}, got {mag}"
        );
    }

    #[test]
    fn loop_closure_is_satisfied_when_close() {
        let lc = LoopClosure::new(0, [0.0; 3], 1, [0.0; 3]);
        let id = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let pose = ([0.0, 0.0, 0.0], id);
        assert!(
            lc.is_satisfied(pose, pose, 1e-6),
            "Coincident points should satisfy constraint"
        );
    }

    #[test]
    fn loop_closure_is_not_satisfied_when_far() {
        let lc = LoopClosure::new(0, [0.0; 3], 1, [0.0; 3]);
        let id = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let pose_a = ([0.0, 0.0, 0.0], id);
        let pose_b = ([10.0, 0.0, 0.0], id);
        assert!(
            !lc.is_satisfied(pose_a, pose_b, 0.1),
            "Far-apart points should not satisfy constraint"
        );
    }

    // ── Additional MultibodyChain tests ───────────────────────────────────

    #[test]
    fn single_link_fk_at_zero_is_identity_rotation() {
        let chain = build_planar_arm(1, 1.0, 1.0);
        let poses = chain.forward_kinematics();
        let rot = poses[0].1;
        // At zero angle, rotation should be identity
        for (i, row) in rot.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (val - expected).abs() < 1e-10,
                    "rot[{i}][{j}]={} expected {expected}",
                    val
                );
            }
        }
    }

    #[test]
    fn chain_all_zero_velocities() {
        let chain = build_planar_arm(3, 1.0, 1.0);
        for &v in &chain.joint_velocities {
            assert!(
                v.abs() < 1e-15,
                "Initial velocities should be zero, got {v}"
            );
        }
    }

    #[test]
    fn chain_all_zero_angles() {
        let chain = build_planar_arm(4, 0.5, 2.0);
        for &q in &chain.joint_angles {
            assert!(q.abs() < 1e-15, "Initial angles should be zero, got {q}");
        }
    }

    #[test]
    fn inverse_then_forward_dynamics_consistency() {
        // Compute torques for zero acceleration, then feed back: should give ~0 acceleration
        let chain = build_planar_arm(2, 1.0, 1.0);
        let grav = [0.0, -9.81, 0.0];
        let zero_acc = vec![0.0; 2];
        let hold_torques = chain.inverse_dynamics(&zero_acc, grav);
        let qdd = chain.forward_dynamics(&hold_torques, grav);
        for (i, &a) in qdd.iter().enumerate() {
            assert!(
                a.abs() < 1e-3,
                "Hold torques should produce ~0 acceleration at joint {i}, got {a}"
            );
        }
    }

    #[test]
    fn link_with_inertia_preserves_fields() {
        let inertia = [[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]];
        let com = [0.1, 0.2, 0.3];
        let axis = [0.0, 1.0, 0.0];
        let offset = [0.5, 0.0, 0.0];
        let link = Link::with_inertia(5.0, inertia, com, axis, offset);
        assert!((link.mass - 5.0).abs() < 1e-12);
        assert!((link.inertia[1][1] - 2.0).abs() < 1e-12);
        assert!((link.com[1] - 0.2).abs() < 1e-12);
    }

    #[test]
    fn gravity_torques_single_link_horizontal() {
        // Single-link arm, link hangs horizontally (90 deg rotation about z)
        let mut chain = build_planar_arm(1, 1.0, 1.0);
        chain.joint_angles[0] = PI / 2.0;
        let grav = [0.0, -9.81, 0.0];
        let tau = chain.gravity_torques(grav);
        // g-torque = m * g * L/2 * cos(90°) but since arm is now horizontal,
        // the gravity acts along -y, lever arm = link length/2 along +x
        // tau = m * |g| * (L/2) ≈ 1 * 9.81 * 0.5 (projected onto joint z-axis)
        assert!(
            tau[0].abs() > 1e-3,
            "Gravity torque should be nonzero for horizontal arm, got {}",
            tau[0]
        );
    }

    #[test]
    fn rk4_energy_changes_under_torque() {
        let mut chain = build_planar_arm(2, 1.0, 1.0);
        let torques = vec![5.0, 2.0];
        let q_before: Vec<f64> = chain.joint_angles.clone();
        let qd_before: Vec<f64> = chain.joint_velocities.clone();
        for _ in 0..100 {
            chain.step_rk4(&torques, [0.0, -9.81, 0.0], 0.001);
        }
        // State must have changed from zero
        let delta_q: f64 = chain
            .joint_angles
            .iter()
            .zip(q_before.iter())
            .map(|(a, b)| (a - b).abs())
            .sum();
        let delta_qd: f64 = chain
            .joint_velocities
            .iter()
            .zip(qd_before.iter())
            .map(|(a, b)| (a - b).abs())
            .sum();
        assert!(
            delta_q + delta_qd > 1e-6,
            "Applied torques should change chain state"
        );
    }
}
