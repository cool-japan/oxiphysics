// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Recursive Newton-Euler Algorithm (RNEA) for inverse dynamics.
//!
//! Given joint positions `q`, velocities `q̇`, and accelerations `q̈`,
//! computes the joint torques `τ` required to produce that motion.
//!
//! Reference: Featherstone 2008 "Rigid Body Dynamics Algorithms", Algorithm 5.1.
//!
//! # Algorithm overview
//!
//! **Pass 1 (outward — root to leaves):** Propagate kinematics.
//! For each body `i`:
//! 1. Compute `X_J = joint.transform(q_i)` and `v_J = S_i · q̇_i`.
//! 2. Compose the full transform: `X = X_T ∘ X_J` (fixed offset then joint motion).
//! 3. Body velocity: `v[i] = X · v[parent] + v_J`.
//! 4. Body acceleration: `a[i] = X · a[parent] + S_i · q̈_i + v[i] × v_J + c_J`.
//!
//! For the root body the parent velocity is zero and the parent acceleration is
//! `a_gravity = −g_vec` (fictitious gravity, expressed as upward acceleration on
//! the root's virtual parent so that it flows outward correctly).
//!
//! **Pass 2 (inward — leaves to root):** Propagate forces.
//! For each body `i` (reverse order):
//! 1. Net force: `f[i] = I[i] · a[i] + v[i] ×* (I[i] · v[i])`.
//! 2. Joint torque: `τ[i] = S_i^T · f[i]`.
//! 3. Propagate to parent: `f[parent] += X^{force-T} · f[i]`.

use crate::{
    model::ArticulatedModel,
    spatial::{SpatialTransform, SpatialVec},
};

/// Compute joint torques via RNEA (inverse dynamics).
///
/// # Arguments
/// - `model`: Articulated-body model (kinematic tree, topological order).
/// - `q`: Joint positions; length must equal `model.total_dof()`.
/// - `q_dot`: Joint velocities; same length as `q`.
/// - `q_ddot`: Joint accelerations; same length as `q`.
///
/// # Returns
/// Joint torques/forces `τ` (length = `model.total_dof()`).
///
/// # Panics
/// Panics if `q`, `q_dot`, or `q_ddot` are not the right length.
pub fn rnea(model: &ArticulatedModel, q: &[f64], q_dot: &[f64], q_ddot: &[f64]) -> Vec<f64> {
    let n = model.num_bodies();
    assert_eq!(
        q.len(),
        model.total_dof(),
        "q length mismatch: got {}, expected {}",
        q.len(),
        model.total_dof()
    );
    assert_eq!(q_dot.len(), model.total_dof(), "q_dot length mismatch");
    assert_eq!(q_ddot.len(), model.total_dof(), "q_ddot length mismatch");

    // Spatial velocity of each body (in body frame)
    let mut v = vec![SpatialVec::ZERO; n];
    // Spatial acceleration of each body (in body frame)
    let mut a = vec![SpatialVec::ZERO; n];
    // Net spatial force on each body (in body frame)
    let mut f = vec![SpatialVec::ZERO; n];
    // Total joint-to-parent transform for each body (X_T ∘ X_J)
    let mut x_total = vec![SpatialTransform::IDENTITY; n];

    // ── Pass 1: Outward kinematics ────────────────────────────────────────────
    //
    // Gravity is treated as a fictitious acceleration of the root's parent.
    // a_parent_of_root = gravity (not negative, because a_root gets X*gravity
    // which carries the sign correctly when gravity=[0,-g,0]).
    // Featherstone convention: a_root = a_gravity propagated outward.
    // We set a_parent_root = gravity so the root sees it correctly.
    let a_gravity = model.gravity; // [0; g_vec], e.g. [0; [0,-9.81,0]]

    for i in 0..n {
        let body = &model.bodies[i];
        let joint = &model.joints[i];
        let qi = model.q_slice(q, i);
        let q_dot_i = model.q_slice(q_dot, i);
        let q_ddot_i = model.q_slice(q_ddot, i);

        // Joint transform X_J and fixed transform X_T
        let x_j = joint.transform(qi);
        let x_t = &body.parent_transform;
        // Full transform: child ← parent (compose X_T then X_J)
        let x_full = x_t.compose(&x_j);
        x_total[i] = x_full;

        // Joint velocity: v_J = Σ_k S_k * q_dot_k
        let s = joint.motion_subspace_at(qi);
        let v_j = subspace_velocity(&s, q_dot_i);

        // Coriolis term c_J = Ṡ * q_dot
        let c_j = joint.coriolis(qi, q_dot_i);

        // Get parent velocity and acceleration (or gravity for root's parent)
        let (v_parent, a_parent) = if let Some(pid) = body.parent_id {
            (v[pid], a[pid])
        } else {
            // Root body: parent is the world. Velocity = 0, acceleration = gravity.
            (SpatialVec::ZERO, a_gravity)
        };

        // Transform parent velocity to child frame
        let v_parent_in_child = x_full.apply_velocity(&v_parent);
        let a_parent_in_child = x_full.apply_velocity(&a_parent);

        // Body velocity: v[i] = X * v_parent + v_J
        v[i] = v_parent_in_child + v_j;

        // Joint acceleration: S * q_ddot_i
        let s_q_ddot = subspace_velocity(&s, q_ddot_i);

        // Body acceleration: a[i] = X*a_parent + S*q_ddot + v[i] × v_J + c_J
        // Note: v[i] × v_J is the Coriolis/centrifugal from velocity propagation
        let v_cross_vj = v[i].cross_motion(&v_j);
        a[i] = a_parent_in_child + s_q_ddot + v_cross_vj + c_j;
    }

    // ── Pass 2: Inward force propagation ─────────────────────────────────────
    let mut tau = vec![0.0f64; model.total_dof()];

    for i in (0..n).rev() {
        let body = &model.bodies[i];
        let joint = &model.joints[i];
        let qi = model.q_slice(q, i);
        let s = joint.motion_subspace_at(qi);

        // Inertia-weighted acceleration: I[i] * a[i]
        let ia = body.inertia.mul_vec(&a[i]);

        // Centrifugal/gyroscopic: v[i] ×* (I[i] * v[i])
        let iv = body.inertia.mul_vec(&v[i]);
        let v_cross_f = v[i].cross_force(&iv);

        // Net force on body i: f[i] = I*a + v ×* (I*v)
        f[i] = f[i] + ia + v_cross_f;

        // Joint generalised force: τ_i = S_i^T * f[i]
        let start = model.dof_start(i);
        for (k, sk) in s.iter().enumerate() {
            tau[start + k] = sk.dot(&f[i]);
        }

        // Propagate force to parent: f[parent] += X^force * f[i]
        if let Some(pid) = body.parent_id {
            let force_in_parent = x_total[i].apply_force(&f[i]);
            f[pid] = f[pid] + force_in_parent;
        }
    }

    tau
}

/// Compute `Σ_k S_k * v_k` from the motion subspace columns and scalar velocities.
fn subspace_velocity(s: &[SpatialVec], v: &[f64]) -> SpatialVec {
    let mut result = SpatialVec::ZERO;
    for (sk, &vk) in s.iter().zip(v.iter()) {
        result = result + sk.scale(vk);
    }
    result
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

    fn approx_eq(a: f64, b: f64, eps: f64) -> bool {
        (a - b).abs() < eps
    }

    /// Build a simple pendulum: revolute about z at origin, point mass at [0,-L,0].
    fn build_pendulum(m: f64, l: f64, g: f64) -> ArticulatedModel {
        let mut model = ArticulatedModel::new([0.0, -g, 0.0]);

        // Revolute joint about z at origin, no fixed offset
        let joint = Box::new(RevoluteJoint::new([0.0, 0.0, 1.0]));

        // Point mass at [0, -l, 0] in the body frame (hanging position at q=0)
        // I_com = 0 for point mass; I_origin = m*(l^2 I - c*c^T)
        let com = [0.0, -l, 0.0];
        let inertia = SpatialInertia::from_com(m, com, [[0.0; 3]; 3]);

        let body = RigidBody::new("link1", inertia, None, SpatialTransform::IDENTITY);
        model.add_body(body, joint);
        model
    }

    #[test]
    fn test_rnea_zero_config() {
        // At q=0 (hanging down), with zero velocity and acceleration,
        // the gravitational torque about z should be zero (gravity acts through pivot).
        let (m, l, g) = (1.0, 1.0, 9.81);
        let model = build_pendulum(m, l, g);

        let tau = rnea(&model, &[0.0], &[0.0], &[0.0]);
        assert!(
            approx_eq(tau[0], 0.0, 1e-8),
            "At q=0 (hanging), tau should be 0, got {:.8}",
            tau[0]
        );
    }

    #[test]
    fn test_rnea_horizontal() {
        // At q=π/2 (horizontal — mass at [l, 0, 0] after rotation), zero vel/acc:
        // Gravity [0,-g,0] creates torque = m*g*l about z (positive rotation
        // of body about z moves COM from [0,-l,0] to [l,0,0]).
        // τ = -m*g*l (restoring torque opposing positive q direction).
        let (m, l, g) = (1.0, 1.0, 9.81);
        let model = build_pendulum(m, l, g);

        let tau = rnea(&model, &[std::f64::consts::FRAC_PI_2], &[0.0], &[0.0]);
        // Expected: torque = m * g * l in the direction opposing positive rotation
        // The RNEA returns the required torque for the given motion. At q=π/2 with
        // zero acceleration, the required torque must cancel gravity → τ = m*g*l.
        let expected = m * g * l;
        assert!(
            approx_eq(tau[0], expected, 1e-6),
            "At q=π/2, tau should be {expected:.6}, got {:.6}",
            tau[0]
        );
    }
}
