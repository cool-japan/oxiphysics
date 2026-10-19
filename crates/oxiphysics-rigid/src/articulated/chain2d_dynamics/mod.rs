// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Rigid-body dynamics of [`ArticulatedChain2D`].
//!
//! # Model
//!
//! The chain is a tree of rigid rods moving in the xy-plane.
//!
//! * Link `i` is a rod of length `l_i` and mass `m_i` whose centre of mass lies
//!   at its midpoint. [`LinkBody::inertia`](super::LinkBody::inertia) is its
//!   moment of inertia about the joint axis at the proximal end (`m l^2 / 3`
//!   for the uniform rod created by
//!   [`add_link`](ArticulatedChain2D::add_link)); about the centre of mass it
//!   is `I_c = max(I - m (l / 2)^2, 0)`.
//! * Joint `i` is a revolute joint about the z axis at the joint origin `o_i`:
//!   the distal end of the parent link, or the world origin for a link without
//!   a parent. `q_i` is [`joint_angle`](super::ArticulatedLink2D::joint_angle),
//!   the angle of link `i` relative to its parent, so link `i` points along the
//!   absolute angle `theta_i = theta_parent + q_i`. The velocity `q'_i` is
//!   [`joint_velocity`](super::ArticulatedLink2D::joint_velocity). Joint types
//!   other than revolute are not distinguished, matching
//!   [`forward_kinematics`](ArticulatedChain2D::forward_kinematics).
//! * The base frame is the world frame (base at the origin with zero angle) and
//!   `gravity` is the acceleration vector `[gx, gy]` in it. A link whose
//!   `parent` index does not refer to an earlier link is attached to the base.
//!
//! # Equations of motion
//!
//! `H(q) q'' + C(q, q') + G(q) = tau`, with
//!
//! * `H_kl = sum_j ( m_j (c_j - o_k) . (c_j - o_l) + I_c,j )`, the sum running
//!   over the links `j` distal to *both* joints `k` and `l` (links of the
//!   subtrees of `k` and `l`; `c_j` is the centre of mass of link `j`). `H`
//!   depends on the configuration and couples every joint with the joints on
//!   its path to the root;
//! * `G_k = - sum_j m_j [(c_j - o_k) x g]_z` over the links distal to joint
//!   `k`: the holding torque against gravity with the true lever arm of every
//!   distal link;
//! * `C` the Coriolis and centrifugal terms, the inverse dynamics at zero
//!   acceleration and without gravity.
//!
//! The inverse dynamics is the recursive Newton-Euler algorithm (`O(n)`); `H`
//! is assembled from the definition and the forward dynamics solves
//! `H q'' = tau - C - G`.

use nalgebra::{DMatrix, DVector};

use super::types::ArticulatedChain2D;
use crate::error::{Error, Result};

#[cfg(test)]
mod tests;

/// Entry `i` of a state/torque slice, or 0 when the slice is too short.
#[inline]
fn coordinate(values: &[f64], i: usize) -> f64 {
    values.get(i).copied().unwrap_or(0.0)
}

/// Difference of two planar vectors.
#[inline]
fn sub2(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [a[0] - b[0], a[1] - b[1]]
}

/// Dot product of two planar vectors.
#[inline]
fn dot2(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}

/// z component of the cross product of two planar vectors.
#[inline]
fn cross2(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[1] - a[1] * b[0]
}

/// Acceleration of a point at offset `r` from a reference point of a body that
/// rotates with angular velocity `omega` and angular acceleration `alpha`:
/// `a_ref + alpha z x r - omega^2 r`.
#[inline]
fn point_acceleration(a_ref: [f64; 2], alpha: f64, omega: f64, r: [f64; 2]) -> [f64; 2] {
    [
        a_ref[0] - alpha * r[1] - omega * omega * r[0],
        a_ref[1] + alpha * r[0] - omega * omega * r[1],
    ]
}

/// Planar kinematics of the chain in its current configuration.
struct TreePose {
    /// Parent of every link; `None` for links attached to the base (including
    /// links whose stored parent index does not refer to an earlier link).
    parent: Vec<Option<usize>>,
    /// Joint origin of every link.
    origin: Vec<[f64; 2]>,
    /// Centre of mass of every link.
    com: Vec<[f64; 2]>,
    /// Distal end of every link.
    tip: Vec<[f64; 2]>,
}

impl TreePose {
    /// Vector from the joint origin of link `i` to its centre of mass.
    #[inline]
    fn com_offset(&self, i: usize) -> [f64; 2] {
        sub2(self.com[i], self.origin[i])
    }

    /// Vector from the joint origin of link `i` to its distal end.
    #[inline]
    fn link_vector(&self, i: usize) -> [f64; 2] {
        sub2(self.tip[i], self.origin[i])
    }
}

impl ArticulatedChain2D {
    /// Joint origins, centres of mass and tips for the current joint angles.
    fn tree_pose(&self) -> TreePose {
        let n = self.links.len();
        let mut pose = TreePose {
            parent: Vec::with_capacity(n),
            origin: Vec::with_capacity(n),
            com: Vec::with_capacity(n),
            tip: Vec::with_capacity(n),
        };
        let mut angle: Vec<f64> = Vec::with_capacity(n);
        for (i, link) in self.links.iter().enumerate() {
            let parent = link.parent.filter(|&p| p < i);
            let (origin, base_angle) = match parent {
                Some(p) => (pose.tip[p], angle[p]),
                None => ([0.0, 0.0], 0.0),
            };
            let absolute_angle = base_angle + link.joint_angle;
            let (sin_a, cos_a) = absolute_angle.sin_cos();
            let length = link.body.length;
            pose.parent.push(parent);
            pose.origin.push(origin);
            pose.com.push([
                origin[0] + 0.5 * length * cos_a,
                origin[1] + 0.5 * length * sin_a,
            ]);
            pose.tip
                .push([origin[0] + length * cos_a, origin[1] + length * sin_a]);
            angle.push(absolute_angle);
        }
        pose
    }

    /// Moment of inertia of link `i` about its centre of mass.
    fn com_inertia(&self, i: usize) -> f64 {
        let body = &self.links[i].body;
        let half_length = 0.5 * body.length;
        (body.inertia - body.mass * half_length * half_length).max(0.0)
    }

    /// Joint velocities of all links.
    fn joint_velocities(&self) -> Vec<f64> {
        self.links.iter().map(|l| l.joint_velocity).collect()
    }

    /// Recursive Newton-Euler algorithm: the joint torques that produce the
    /// accelerations `qdd` in the state `(q, qd)` of `pose` under `gravity`.
    ///
    /// Forward pass (roots to leaves): absolute angular velocity and
    /// acceleration of every link and the acceleration of every joint origin and
    /// centre of mass, with gravity folded in as a fictitious base acceleration
    /// `-g`. Backward pass (leaves to roots): the force `f_i = m_i a_ci + sum
    /// f_k` transmitted through joint `i` and the moment about its origin
    /// `N_i = I_c alpha_i + r_ci x m_i a_ci + sum (N_k + L_i x f_k)` over the
    /// child joints `k` (`L_i` is the link vector); the joint torque is `N_i`.
    fn recursive_newton_euler(
        &self,
        pose: &TreePose,
        qd: &[f64],
        qdd: &[f64],
        gravity: [f64; 2],
    ) -> Vec<f64> {
        let n = self.links.len();
        let mut omega = vec![0.0_f64; n];
        let mut alpha = vec![0.0_f64; n];
        let mut origin_acceleration = vec![[0.0_f64; 2]; n];
        let mut force = vec![[0.0_f64; 2]; n];
        let mut torque = vec![0.0_f64; n];

        for i in 0..n {
            let (omega_parent, alpha_parent, origin_acc) = match pose.parent[i] {
                Some(p) => (
                    omega[p],
                    alpha[p],
                    point_acceleration(
                        origin_acceleration[p],
                        alpha[p],
                        omega[p],
                        pose.link_vector(p),
                    ),
                ),
                None => (0.0, 0.0, [-gravity[0], -gravity[1]]),
            };
            omega[i] = omega_parent + coordinate(qd, i);
            alpha[i] = alpha_parent + coordinate(qdd, i);
            origin_acceleration[i] = origin_acc;
            let r_com = pose.com_offset(i);
            let com_acc = point_acceleration(origin_acc, alpha[i], omega[i], r_com);
            let mass = self.links[i].body.mass;
            force[i] = [mass * com_acc[0], mass * com_acc[1]];
            torque[i] = self.com_inertia(i) * alpha[i] + cross2(r_com, force[i]);
        }

        for i in (0..n).rev() {
            if let Some(p) = pose.parent[i] {
                let transmitted = force[i];
                torque[p] += torque[i] + cross2(pose.link_vector(p), transmitted);
                force[p][0] += transmitted[0];
                force[p][1] += transmitted[1];
            }
        }
        torque
    }

    /// Joint-space inertia matrix `H(q)`, row-major `n x n`.
    ///
    /// `H_kl = sum_j ( m_j (c_j - o_k) . (c_j - o_l) + I_c,j )` over the links
    /// `j` distal to both joints `k` and `l`. The matrix is symmetric, positive
    /// semi-definite (positive definite when every link carries mass or
    /// inertia) and depends on the configuration. Joints of different branches
    /// of the tree do not couple.
    pub fn mass_matrix(&self) -> Vec<f64> {
        let n = self.links.len();
        let pose = self.tree_pose();
        let mut h = vec![0.0_f64; n * n];
        let mut path: Vec<usize> = Vec::new();
        for j in 0..n {
            path.clear();
            let mut cursor = Some(j);
            while let Some(k) = cursor {
                path.push(k);
                cursor = pose.parent[k];
            }
            let mass = self.links[j].body.mass;
            let inertia = self.com_inertia(j);
            for &k in &path {
                let r_k = sub2(pose.com[j], pose.origin[k]);
                for &l in &path {
                    let r_l = sub2(pose.com[j], pose.origin[l]);
                    h[k * n + l] += mass * dot2(r_k, r_l) + inertia;
                }
            }
        }
        h
    }

    /// Gravity torques `G(q)`: the joint torques that hold the chain still
    /// against `gravity` (`[gx, gy]`, in the world frame).
    ///
    /// `G_k = - sum_j m_j [(c_j - o_k) x g]_z` over the links distal to joint
    /// `k`. For `g = (0, -9.81)` a link pointing along +x needs a positive
    /// (counter-clockwise) holding torque `m g l / 2`.
    pub fn gravity_torques(&self, gravity: [f64; 2]) -> Vec<f64> {
        let pose = self.tree_pose();
        self.recursive_newton_euler(&pose, &[], &[], gravity)
    }

    /// Coriolis and centrifugal joint torques `C(q, q')` of the current state:
    /// the inverse dynamics at zero joint acceleration and without gravity.
    pub fn coriolis_forces(&self) -> Vec<f64> {
        let pose = self.tree_pose();
        self.recursive_newton_euler(&pose, &self.joint_velocities(), &[], [0.0, 0.0])
    }

    /// Newton-Euler inverse dynamics: the joint torques
    /// `H(q) q'' + C(q, q') + G(q)` that produce the joint accelerations
    /// `joint_accelerations` in the current state under `gravity`
    /// (`[gx, gy]`, in the world frame with the chain base at the origin).
    ///
    /// The state is the `joint_angle` and `joint_velocity` of the links; missing
    /// entries of `joint_accelerations` count as zero, extra entries are
    /// ignored. The result has one torque per link.
    pub fn inverse_dynamics(&self, joint_accelerations: &[f64], gravity: [f64; 2]) -> Vec<f64> {
        let pose = self.tree_pose();
        self.recursive_newton_euler(
            &pose,
            &self.joint_velocities(),
            joint_accelerations,
            gravity,
        )
    }

    /// The dynamics equation `H q'' = tau - C - G` as the matrix `H` (row-major)
    /// and the right-hand side.
    fn equations_of_motion(&self, torques: &[f64], gravity: [f64; 2]) -> (Vec<f64>, DVector<f64>) {
        let n = self.links.len();
        let pose = self.tree_pose();
        let bias = self.recursive_newton_euler(&pose, &self.joint_velocities(), &[], gravity);
        let rhs = DVector::from_iterator(n, (0..n).map(|i| coordinate(torques, i) - bias[i]));
        (self.mass_matrix(), rhs)
    }

    /// Forward dynamics: the joint accelerations `q''` produced by the joint
    /// torques `torques` in the current state under `gravity`, from
    /// `H(q) q'' = tau - C(q, q') - G(q)` solved with a Cholesky factorisation
    /// of the joint-space inertia matrix. Missing torques count as zero.
    ///
    /// Unlike a chain of independent joints, the joints are coupled: a torque on
    /// one joint accelerates every joint it shares a branch with.
    ///
    /// # Errors
    ///
    /// Returns [`Error::ForwardDynamicsFailure`] when `H` is not positive
    /// definite (for example a link without mass and inertia, or non-finite
    /// state) or the resulting accelerations are not finite.
    pub fn forward_dynamics(&self, torques: &[f64], gravity: [f64; 2]) -> Result<Vec<f64>> {
        let n = self.links.len();
        if n == 0 {
            return Ok(Vec::new());
        }
        let (h, rhs) = self.equations_of_motion(torques, gravity);
        let factorisation = DMatrix::from_row_slice(n, n, &h)
            .cholesky()
            .ok_or_else(|| {
                let body_index = (0..n)
                    .find(|&i| !(h[i * n + i].is_finite() && h[i * n + i] > 0.0))
                    .unwrap_or(0);
                Error::forward_dynamics_failure(
                    body_index,
                    "joint-space inertia matrix is not positive definite",
                )
            })?;
        let qdd = factorisation.solve(&rhs);
        match qdd.iter().position(|value| !value.is_finite()) {
            Some(body_index) => Err(Error::forward_dynamics_failure(
                body_index,
                "joint acceleration is not finite",
            )),
            None => Ok(qdd.iter().copied().collect()),
        }
    }

    /// Joint accelerations for the integration step: the solution of
    /// [`forward_dynamics`](Self::forward_dynamics) when `H` is positive
    /// definite; otherwise the minimum-norm least-squares solution, in which
    /// joints that carry no inertia cannot be accelerated and keep zero
    /// acceleration. Non-finite input gives NaN accelerations.
    fn integration_accelerations(&self, torques: &[f64], gravity: [f64; 2]) -> Vec<f64> {
        let n = self.links.len();
        if n == 0 {
            return Vec::new();
        }
        let (h, rhs) = self.equations_of_motion(torques, gravity);
        if !h.iter().all(|v| v.is_finite()) || !rhs.iter().all(|v| v.is_finite()) {
            return vec![f64::NAN; n];
        }
        let matrix = DMatrix::from_row_slice(n, n, &h);
        if let Some(factorisation) = matrix.clone().cholesky() {
            return factorisation.solve(&rhs).iter().copied().collect();
        }
        let solution = matrix
            .try_svd(true, true, f64::EPSILON, 1000)
            .and_then(|svd| {
                let largest = svd.singular_values.iter().copied().fold(0.0_f64, f64::max);
                svd.solve(&rhs, largest * n as f64 * f64::EPSILON).ok()
            });
        match solution {
            Some(qdd) => qdd.iter().copied().collect(),
            None => vec![0.0; n],
        }
    }

    /// Integrate joint velocities and angles by `dt` under the applied joint
    /// torques, without gravity.
    ///
    /// Solves `H(q) q'' = tau - C(q, q')` for the joint accelerations (so the
    /// joints are dynamically coupled, see
    /// [`forward_dynamics`](Self::forward_dynamics)) and then takes a
    /// semi-implicit Euler step: `q'_new = q' + q'' dt`, `q_new = q + q'_new dt`.
    /// Missing torques count as zero. Joints without any inertia cannot be
    /// accelerated and stay locked.
    pub fn apply_torques(&mut self, torques: &[f64], dt: f64) {
        self.apply_torques_with_gravity(torques, [0.0, 0.0], dt);
    }

    /// Like [`apply_torques`](Self::apply_torques) with gravity (`[gx, gy]`, in
    /// the world frame with the chain base at the origin) acting on the links:
    /// the joint accelerations solve `H q'' = tau - C - G`.
    pub fn apply_torques_with_gravity(&mut self, torques: &[f64], gravity: [f64; 2], dt: f64) {
        let accelerations = self.integration_accelerations(torques, gravity);
        for (link, acceleration) in self.links.iter_mut().zip(accelerations) {
            link.joint_velocity += acceleration * dt;
            link.joint_angle += link.joint_velocity * dt;
        }
    }

    /// Total kinetic energy of the chain, `1/2 q'^T H(q) q'`: the translation
    /// of the centre of mass of every link and its rotation,
    /// `sum_j ( m_j |v_cj|^2 / 2 + I_c,j omega_j^2 / 2 )`.
    pub fn total_kinetic_energy(&self) -> f64 {
        let n = self.links.len();
        let pose = self.tree_pose();
        let mut omega = vec![0.0_f64; n];
        let mut origin_velocity = vec![[0.0_f64; 2]; n];
        let mut energy = 0.0_f64;
        for i in 0..n {
            let (omega_parent, origin_vel) = match pose.parent[i] {
                Some(p) => {
                    let l = pose.link_vector(p);
                    (
                        omega[p],
                        [
                            origin_velocity[p][0] - omega[p] * l[1],
                            origin_velocity[p][1] + omega[p] * l[0],
                        ],
                    )
                }
                None => (0.0, [0.0, 0.0]),
            };
            omega[i] = omega_parent + self.links[i].joint_velocity;
            origin_velocity[i] = origin_vel;
            let r = pose.com_offset(i);
            let com_velocity = [
                origin_vel[0] - omega[i] * r[1],
                origin_vel[1] + omega[i] * r[0],
            ];
            energy += 0.5 * self.links[i].body.mass * dot2(com_velocity, com_velocity)
                + 0.5 * self.com_inertia(i) * omega[i] * omega[i];
        }
        energy
    }
}
