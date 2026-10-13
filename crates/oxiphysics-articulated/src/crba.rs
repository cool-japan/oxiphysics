// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Composite Rigid Body Algorithm (CRBA) for O(n²) mass-matrix computation.
//!
//! Reference: Featherstone 2008 "Rigid Body Dynamics Algorithms", Algorithm 6.2.

use crate::{
    model::ArticulatedModel,
    spatial::{SpatialInertia6x6, SpatialTransform},
};

/// Per-body parent←child poses `X_T ∘ X_J(q_i)` for configuration `q`.
pub(crate) fn joint_poses(model: &ArticulatedModel, q: &[f64]) -> Vec<SpatialTransform> {
    (0..model.num_bodies())
        .map(|i| {
            let qi = model.q_slice(q, i);
            let x_j = model.joints[i].transform(qi);
            model.bodies[i].parent_transform.compose(&x_j)
        })
        .collect()
}

/// Composite (subtree) spatial inertias `I^c_i`, each expressed in body `i`'s
/// own frame (Featherstone §6.2, composite-inertia backward pass).
///
/// `x_full[i]` must be the pose of body `i` relative to its parent (as returned
/// by [`joint_poses`]).
pub(crate) fn composite_inertias(
    model: &ArticulatedModel,
    x_full: &[SpatialTransform],
) -> Vec<SpatialInertia6x6> {
    let n = model.num_bodies();
    let mut i_c: Vec<SpatialInertia6x6> = (0..n)
        .map(|i| SpatialInertia6x6::from_packed(&model.bodies[i].inertia))
        .collect();
    for i in (0..n).rev() {
        if let Some(pid) = model.bodies[i].parent_id {
            let i_c_in_parent = i_c[i].transform_to_parent(&x_full[i]);
            i_c[pid] = i_c[pid].add(&i_c_in_parent);
        }
    }
    i_c
}

/// Compute the joint-space mass matrix `M(q)` using the Composite Rigid Body Algorithm.
///
/// Returns an `n×n` symmetric positive-definite matrix where
/// `n = model.total_dof()` (velocity dimension). `q` is the configuration
/// vector of length [`ArticulatedModel::total_nq`].
/// Complexity: O(n²) — faster than the `n` RNEA calls that the naive method requires.
///
/// # Algorithm (Featherstone §6.2)
///
/// 1. **Forward pass**: compute the pose `X_T ∘ X_J(q_i)` of each body in its parent.
/// 2. **Composite-inertia backward pass**: accumulate `I^c_i = I_i + Σ_children X_cᵀ · I^c_c · X_c`.
/// 3. **Mass-matrix assembly**: for each body `i` and DOF `k`, form the force
///    column `F_k = I^c_i · S_k`; the block of body `i`'s own joint is
///    `M[i_k][i_l] = S_lᵀ · F_k` (all `l`, including the off-diagonal entries of
///    multi-DoF joints), then `F_k` is carried up the chain to fill the
///    ancestor blocks.
///
/// # Panics
/// Panics if `q.len() != model.total_nq()`.
pub fn compute_mass_matrix_crba(model: &ArticulatedModel, q: &[f64]) -> Vec<Vec<f64>> {
    assert_eq!(
        q.len(),
        model.total_nq(),
        "q length mismatch: got {}, expected {}",
        q.len(),
        model.total_nq()
    );
    let total = model.total_dof();

    // ── Pass 1: joint transforms ─────────────────────────────────────────────
    let x_full = joint_poses(model, q);

    // ── Pass 2: composite-inertia backward pass ───────────────────────────────
    let i_c = composite_inertias(model, &x_full);

    // ── Pass 3: mass-matrix assembly ─────────────────────────────────────────
    let mut m = vec![vec![0.0f64; total]; total];

    for (i, i_c_i) in i_c.iter().enumerate() {
        let qi = model.q_slice(q, i);
        let s_i = model.joints[i].motion_subspace_at(qi);
        if s_i.is_empty() {
            continue;
        }
        let start_i = model.dof_start(i);

        // For each DOF k of body i, propagate force F_k upward
        for (k, s_ik) in s_i.iter().enumerate() {
            let mut f_k = i_c_i.mul_vec6(s_ik);
            let row = start_i + k;

            // Own-joint block: M[row][start_i + l] = S_i[l]ᵀ · F_k for every
            // column l of this joint (diagonal and intra-joint off-diagonal).
            for (l, s_il) in s_i.iter().enumerate() {
                m[row][start_i + l] = s_il.dot(&f_k);
            }

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
        let m = compute_mass_matrix_crba(&model, &q);
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
