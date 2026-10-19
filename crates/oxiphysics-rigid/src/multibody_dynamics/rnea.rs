// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Recursive Newton-Euler algorithm (RNEA) for [`MultibodyChain`].
//!
//! The pass works entirely in the world frame.  Every revolute joint `i` has
//! the world-frame axis `aᵢ` through its joint origin `pᵢ`; the joint origin
//! is a point fixed in the parent link.

use super::{
    MultibodyChain, mat3_mul, mat3_transpose, mat3_vec_mul, unit_vector, vec3_add, vec3_cross,
    vec3_dot, vec3_scale, vec3_sub,
};

impl MultibodyChain {
    /// Recursive Newton-Euler inverse dynamics in the world frame.
    ///
    /// Returns the joint torques `τ` that realise the joint accelerations
    /// `q_ddot` for the current configuration, the joint velocities `q_dot`
    /// and a base frame that accelerates with `base_acceleration`:
    ///
    /// * `base_acceleration = -gravity` gives the full inverse dynamics
    ///   `τ = M(q)·q̈ + C(q, q̇)·q̇ + g(q)` (a uniform gravity field `g` is
    ///   equivalent to accelerating the base by `-g`);
    /// * `q_ddot` empty and `base_acceleration = [0; 3]` gives the Coriolis
    ///   and centrifugal vector `C(q, q̇)·q̇`;
    /// * `q_dot` empty, `q_ddot` empty and `base_acceleration = -gravity`
    ///   gives the gravity (holding) torques `g(q)`.
    ///
    /// Entries missing from `q_dot` or `q_ddot` are treated as zero.  The
    /// chain must be stored parents-first (`parent[i] < i`), as required by
    /// [`MultibodyChain::forward_kinematics`].
    pub fn rnea(&self, q_dot: &[f64], q_ddot: &[f64], base_acceleration: [f64; 3]) -> Vec<f64> {
        let n = self.n_dof();
        let poses = self.forward_kinematics();
        let rate = |values: &[f64], i: usize| values.get(i).copied().unwrap_or(0.0);

        let mut axes = vec![[0.0_f64; 3]; n];
        let mut omega = vec![[0.0_f64; 3]; n];
        let mut alpha = vec![[0.0_f64; 3]; n];
        let mut origin_accel = vec![[0.0_f64; 3]; n];
        // Net inertial force of each link, and its moment about the joint origin.
        let mut force = vec![[0.0_f64; 3]; n];
        let mut moment = vec![[0.0_f64; 3]; n];

        // Forward pass: velocities and accelerations, root to leaves.
        for i in 0..n {
            let (joint_pos, rot) = poses[i];
            let link = &self.links[i];
            let axis = unit_vector(mat3_vec_mul(rot, link.joint_axis));
            axes[i] = axis;

            let (omega_parent, alpha_parent, accel_origin) = match self.parent[i] {
                None => ([0.0; 3], [0.0; 3], base_acceleration),
                Some(p) => {
                    let r = vec3_sub(joint_pos, poses[p].0);
                    let accel = vec3_add(
                        origin_accel[p],
                        vec3_add(
                            vec3_cross(alpha[p], r),
                            vec3_cross(omega[p], vec3_cross(omega[p], r)),
                        ),
                    );
                    (omega[p], alpha[p], accel)
                }
            };
            // The axis is fixed in the parent link, so it rotates with `omega_parent`.
            let relative_rate = vec3_scale(axis, rate(q_dot, i));
            omega[i] = vec3_add(omega_parent, relative_rate);
            alpha[i] = vec3_add(
                vec3_add(alpha_parent, vec3_scale(axis, rate(q_ddot, i))),
                vec3_cross(omega_parent, relative_rate),
            );
            origin_accel[i] = accel_origin;

            // Centre-of-mass acceleration and Newton-Euler equations.
            let r_com = mat3_vec_mul(rot, link.com);
            let accel_com = vec3_add(
                accel_origin,
                vec3_add(
                    vec3_cross(alpha[i], r_com),
                    vec3_cross(omega[i], vec3_cross(omega[i], r_com)),
                ),
            );
            let inertia = mat3_mul(rot, mat3_mul(link.inertia, mat3_transpose(rot)));
            let link_force = vec3_scale(accel_com, link.mass);
            let moment_com = vec3_add(
                mat3_vec_mul(inertia, alpha[i]),
                vec3_cross(omega[i], mat3_vec_mul(inertia, omega[i])),
            );
            force[i] = link_force;
            moment[i] = vec3_add(moment_com, vec3_cross(r_com, link_force));
        }

        // Backward pass: accumulate the wrenches transmitted through the joints.
        let mut torques = vec![0.0_f64; n];
        for i in (0..n).rev() {
            torques[i] = vec3_dot(axes[i], moment[i]);
            if let Some(p) = self.parent[i] {
                let r = vec3_sub(poses[i].0, poses[p].0);
                let child_force = force[i];
                let child_moment = moment[i];
                force[p] = vec3_add(force[p], child_force);
                moment[p] = vec3_add(
                    moment[p],
                    vec3_add(child_moment, vec3_cross(r, child_force)),
                );
            }
        }
        torques
    }
}
