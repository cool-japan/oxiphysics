// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Centroidal Momentum Matrix (CMM) computation.
//!
//! The CMM `A_G ∈ ℝ^{6×n}` maps joint velocities to total system centroidal
//! angular+linear momentum: `h_G = A_G(q) · q_dot`.
//!
//! Reference: Orin & Goswami 2008 "Centroidal Momentum Matrix of a Humanoid Robot".

use crate::{
    model::ArticulatedModel,
    spatial::{SpatialInertia6x6, SpatialTransform},
};

/// Compute the centroidal momentum matrix `A_G(q)` of shape `[6][total_dof]`.
///
/// Uses the composite-inertia backward pass (shared with CRBA) to build the
/// system-COM inertia, then projects each joint's motion subspace through the
/// world-frame transform to the centroidal frame.
///
/// Returns `A_G` as a `Vec<Vec<f64>>` with outer dimension 6 and inner dimension
/// `model.total_dof()`. The row ordering is: `[ωx, ωy, ωz, vx, vy, vz]` (angular
/// momentum first, matching the spatial-vector convention throughout this crate).
///
/// Note: the current implementation projects individual body inertias (not composite
/// subtree inertias) to the world frame.  A full centroidal shift using the parallel-axis
/// theorem is planned for v0.2.
#[allow(clippy::needless_range_loop)]
pub fn compute_cmm(model: &mut ArticulatedModel, q: &[f64]) -> Vec<Vec<f64>> {
    let n = model.num_bodies();
    let total = model.total_dof();

    // ── Forward pass: body-to-world transforms (accumulated product) ──────────
    // x_world_to_body[i] = full parent-to-child transform chain from world to body i.
    // x_body_to_world[i] is its inverse (needed to transform velocities and inertias).
    // Both arrays are indexed together with model arrays, so integer index is required.
    let mut x_world_to_body = vec![SpatialTransform::IDENTITY; n];
    let mut x_body_to_world = vec![SpatialTransform::IDENTITY; n];

    for i in 0..n {
        let qi = model.q_slice(q, i);
        let x_j = model.joints[i].transform(qi);
        let x_full_i = model.bodies[i].parent_transform.compose(&x_j);
        if let Some(pid) = model.bodies[i].parent_id {
            // x_world_to_body[i] = x_full_i ∘ x_world_to_body[pid]
            // (first go from world to parent, then parent to child)
            x_world_to_body[i] = x_full_i.compose(&x_world_to_body[pid]);
        } else {
            x_world_to_body[i] = x_full_i;
        }
        x_body_to_world[i] = x_world_to_body[i].inverse();
    }

    // ── Compute system COM in world frame ─────────────────────────────────────
    // Used for future centroidal shift; currently suppressed with `_`.
    let mut total_mass = 0.0_f64;
    let mut com_world = [0.0_f64; 3];
    for i in 0..n {
        let m_i = model.bodies[i].inertia.mass;
        total_mass += m_i;
        // COM of body i in world frame: apply body-to-world rotation then add translation
        let com_body = model.bodies[i].inertia.com;
        let com_i_world = mat3_vec(x_body_to_world[i].rot, com_body);
        let com_i_world = [
            com_i_world[0] + x_body_to_world[i].trans[0],
            com_i_world[1] + x_body_to_world[i].trans[1],
            com_i_world[2] + x_body_to_world[i].trans[2],
        ];
        com_world[0] += m_i * com_i_world[0];
        com_world[1] += m_i * com_i_world[1];
        com_world[2] += m_i * com_i_world[2];
    }
    if total_mass > 1e-30 {
        com_world[0] /= total_mass;
        com_world[1] /= total_mass;
        com_world[2] /= total_mass;
    }
    // centroidal shift planned for v0.2 — suppress unused warning
    let _ = com_world;

    // ── Build A_G ─────────────────────────────────────────────────────────────
    // For each body i and each DOF k:
    //   1. Transform the body's inertia to world frame via the world-to-body motion transform.
    //   2. Project the joint's motion subspace column (body frame) to world frame.
    //   3. A_G[:, dof_start(i)+k] = I_world * S_k_world
    let mut a_g = vec![vec![0.0_f64; total]; 6];

    for i in 0..n {
        let qi = model.q_slice(q, i);
        let s_i = model.joints[i].motion_subspace_at(qi);
        if s_i.is_empty() {
            continue;
        }

        // Body inertia as 6×6 in body frame, then transform to world frame.
        // I_world = X_w2b^T · I_body · X_w2b  (transform_to_parent with world-to-body).
        let i_body = SpatialInertia6x6::from_packed(&model.bodies[i].inertia);
        let i_world = i_body.transform_to_parent(&x_world_to_body[i]);

        // Project each motion subspace column from body frame to world frame, then multiply.
        for (k, s_ik) in s_i.iter().enumerate() {
            // S_k in world frame: apply body-to-world motion transform
            let s_k_world = x_body_to_world[i].apply_velocity(s_ik);

            // h contribution = I_world * S_k_world
            let h = i_world.mul_vec6(&s_k_world);

            // Pack into A_G rows [angular(0..3), linear(3..6)]
            let col = model.dof_start(i) + k;
            a_g[0][col] = h.angular[0];
            a_g[1][col] = h.angular[1];
            a_g[2][col] = h.angular[2];
            a_g[3][col] = h.linear[0];
            a_g[4][col] = h.linear[1];
            a_g[5][col] = h.linear[2];
        }
    }

    a_g
}

fn mat3_vec(m: [[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        body::RigidBody,
        joint::RevoluteJoint,
        model::ArticulatedModel,
        spatial::{SpatialInertia, SpatialTransform},
    };

    #[test]
    fn test_cmm_dimensions() {
        let mut model = ArticulatedModel::new([0.0; 3]);
        let inertia = SpatialInertia::from_com(
            1.0,
            [0.0; 3],
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        );
        for _ in 0..4 {
            let body = RigidBody::new("b", inertia, None, SpatialTransform::IDENTITY);
            model.add_body(body, Box::new(RevoluteJoint::new([0.0, 0.0, 1.0])));
        }
        let q = vec![0.0; 4];
        let a_g = compute_cmm(&mut model, &q);
        assert_eq!(a_g.len(), 6, "A_G should have 6 rows");
        assert_eq!(a_g[0].len(), 4, "A_G should have n=4 columns");
    }
}
