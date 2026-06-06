// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Composite Rigid Body Algorithm (CRBA) for O(n²) mass-matrix computation.
//!
//! Reference: Featherstone 2008 "Rigid Body Dynamics Algorithms", Algorithm 6.2.

use crate::{
    model::ArticulatedModel,
    spatial::{SpatialInertia6x6, SpatialTransform},
};

/// Compute the joint-space mass matrix `M(q)` using the Composite Rigid Body Algorithm.
///
/// Returns an `n×n` symmetric positive-definite matrix where `n = model.total_dof()`.
/// Complexity: O(n²) — faster than the `n` RNEA calls that the naive method requires.
///
/// # Algorithm (Featherstone §6.2)
///
/// 1. **Forward pass**: compute the full parent-to-child transform `X_i(q)` for each body.
/// 2. **Composite-inertia backward pass**: accumulate `I^c_i = I_i + Σ_children X_c^T · I^c_c · X_c`.
/// 3. **Mass-matrix assembly**: for each body `i` and DOF `k`, propagate the force column
///    `F_k = I^c_i · S_k` upward through the kinematic chain.
pub fn compute_mass_matrix_crba(model: &mut ArticulatedModel, q: &[f64]) -> Vec<Vec<f64>> {
    let n = model.num_bodies();
    let total = model.total_dof();

    // ── Pass 1: joint transforms ─────────────────────────────────────────────
    // Index i required to access model.bodies[i], model.joints[i], and x_full[i] together.
    let mut x_full = vec![SpatialTransform::IDENTITY; n];
    for (i, x_full_i) in x_full.iter_mut().enumerate() {
        let qi = model.q_slice(q, i);
        let x_j = model.joints[i].transform(qi);
        *x_full_i = model.bodies[i].parent_transform.compose(&x_j);
    }

    // ── Pass 2: composite-inertia backward pass ───────────────────────────────
    // Index i required to access model.bodies[i], i_c[i], and i_c[pid] together.
    let mut i_c: Vec<SpatialInertia6x6> = (0..n)
        .map(|i| SpatialInertia6x6::from_packed(&model.bodies[i].inertia))
        .collect();

    for i in (0..n).rev() {
        if let Some(pid) = model.bodies[i].parent_id {
            // Transform I^c_i to parent frame and accumulate
            let i_c_in_parent = i_c[i].transform_to_parent(&x_full[i]);
            i_c[pid] = i_c[pid].add(&i_c_in_parent);
        }
    }

    // ── Pass 3: mass-matrix assembly ─────────────────────────────────────────
    // Index i required to access model.bodies[i], model.joints[i], i_c[i] together.
    let mut m = vec![vec![0.0f64; total]; total];

    for (i, i_c_i) in i_c.iter().enumerate() {
        let qi = model.q_slice(q, i);
        let s_i = model.joints[i].motion_subspace_at(qi);
        if s_i.is_empty() {
            continue;
        }

        // For each DOF k of body i, propagate force F_k upward
        for (k, s_ik) in s_i.iter().enumerate() {
            let mut f_k = i_c_i.mul_vec6(s_ik);
            let row = model.dof_start(i) + k;

            // Diagonal block: M[row][row] = S_i[k]^T · F_k
            m[row][row] = s_ik.dot(&f_k);

            // Walk up the kinematic chain
            let mut j = i;
            while let Some(pid) = model.bodies[j].parent_id {
                // Transform F_k to parent frame
                f_k = x_full[j].apply_force(&f_k);
                j = pid;

                let qj = model.q_slice(q, j);
                let s_j = model.joints[j].motion_subspace_at(qj);
                for (l, s_jl) in s_j.iter().enumerate() {
                    let col = model.dof_start(j) + l;
                    let val = s_jl.dot(&f_k);
                    m[row][col] = val;
                    m[col][row] = val; // symmetric
                }
            }
        }
    }

    m
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
    fn test_crba_single_revolute_matches_inertia() {
        // 1-body revolute about z, unit mass at (1, 0, 0).
        // I_z = m * r^2 = 1 * 1^2 = 1.0
        let mut model = ArticulatedModel::new([0.0; 3]);
        let inertia = SpatialInertia::from_com(1.0, [1.0, 0.0, 0.0], [[0.0; 3]; 3]);
        let body = RigidBody::new("b", inertia, None, SpatialTransform::IDENTITY);
        model.add_body(body, Box::new(RevoluteJoint::new([0.0, 0.0, 1.0])));
        let q = vec![0.0];
        let m = compute_mass_matrix_crba(&mut model, &q);
        assert_eq!(m.len(), 1);
        let diff = (m[0][0] - 1.0).abs();
        assert!(
            diff < 1e-10,
            "M[0][0] should be ~1.0, got {} diff={:.2e}",
            m[0][0],
            diff
        );
    }
}
