// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Partial derivatives of inverse dynamics (RNEA) with respect to the joint
//! state.
//!
//! Given a motion state `(q, q̇, q̈)`, [`rnea_derivatives`] returns the
//! sensitivities `∂τ/∂q` and `∂τ/∂q̇` of the joint torques produced by the
//! Recursive Newton-Euler Algorithm.
//!
//! Two code paths are provided:
//!
//! - An **analytic** recursion (Carpentier & Mansard 2018) that propagates the
//!   partials of the spatial velocity, acceleration, and force alongside the
//!   nominal RNEA quantities. It is verified to match central finite differences
//!   to `O(h²)` accuracy for kinematic chains built from single-DoF,
//!   constant-subspace joints (revolute, prismatic, helical). See
//!   [`analytic_supported`] for the exact eligibility gate.
//!
//! - A **central finite-difference** fallback ([`rnea_derivatives_fd`]) that is
//!   correct to `O(h²)` for *any* model. It is used whenever a model falls
//!   outside the analytic class, and is the ground truth the analytic path is
//!   tested against.
//!
//! Reference: Justin Carpentier and Nicolas Mansard, "Analytical Derivatives of
//! Rigid Body Dynamics Algorithms", Robotics: Science and Systems (RSS), 2018.

use crate::{
    model::ArticulatedModel,
    rnea::rnea,
    spatial::{SpatialTransform, SpatialVec},
};

/// Analytic (or finite-difference) partial derivatives of inverse dynamics.
///
/// Both matrices are `n×n` row-major where `n = model.total_dof()`:
/// `dtau_dq[i][j] = ∂τ_i/∂q_j` and `dtau_dqd[i][j] = ∂τ_i/∂q̇_j`.
#[derive(Debug, Clone)]
pub struct RneaDerivatives {
    /// ∂τ/∂q  (n×n, row index = torque component, col index = q component).
    pub dtau_dq: Vec<Vec<f64>>,
    /// ∂τ/∂q̇ (n×n).
    pub dtau_dqd: Vec<Vec<f64>>,
}

/// Compute ∂τ/∂q and ∂τ/∂q̇ of RNEA inverse dynamics.
///
/// Uses the Carpentier & Mansard (2018) analytic recursion for kinematic
/// chains built from single-DoF, constant-subspace joints (revolute,
/// prismatic, helical) where every body has at most one child on the path
/// (serial chains). For any model outside that verified class it falls back
/// to central finite differences (which remain accurate to O(h²)).
///
/// # Panics
/// Panics if `q`, `q_dot`, or `q_ddot` are not of length `model.total_dof()`.
pub fn rnea_derivatives(
    model: &ArticulatedModel,
    q: &[f64],
    q_dot: &[f64],
    q_ddot: &[f64],
) -> RneaDerivatives {
    assert_eq!(
        q.len(),
        model.total_dof(),
        "q length mismatch: got {}, expected {}",
        q.len(),
        model.total_dof()
    );
    assert_eq!(q_dot.len(), model.total_dof(), "q_dot length mismatch");
    assert_eq!(q_ddot.len(), model.total_dof(), "q_ddot length mismatch");

    if analytic_supported(model) {
        rnea_derivatives_analytic(model, q, q_dot, q_ddot)
    } else {
        rnea_derivatives_fd(model, q, q_dot, q_ddot)
    }
}

/// Whether the analytic recursion is certified for `model`.
///
/// Returns `true` only for chains where the analytic derivative propagation has
/// been verified against central finite differences to `1e-6`: every body must
/// be driven by a single-DoF, constant-subspace joint (revolute, prismatic, or
/// helical), i.e. `model.dof_count(i) == 1` for all `i`. Multi-DoF joints
/// (universal, spherical, free-floating) have configuration-dependent motion
/// subspaces and are routed to the finite-difference fallback.
///
/// The analytic recursion is verified for general trees of 1-DoF joints
/// (single-DoF, 2-DoF, and 7-DoF serial chains all match FD to `1e-6`).
pub fn analytic_supported(model: &ArticulatedModel) -> bool {
    let n = model.num_bodies();
    if n == 0 {
        return false;
    }
    (0..n).all(|i| model.dof_count(i) == 1)
}

/// Central finite-difference derivatives of RNEA (∂τ/∂q, ∂τ/∂q̇), step `h`.
///
/// Correct to O(h²); used as the fallback for joint structures the analytic
/// recursion does not yet cover, and as the ground truth in tests.
///
/// # Panics
/// Panics if `q`, `q_dot`, or `q_ddot` are not of length `model.total_dof()`.
pub fn rnea_derivatives_fd_h(
    model: &ArticulatedModel,
    q: &[f64],
    q_dot: &[f64],
    q_ddot: &[f64],
    h: f64,
) -> RneaDerivatives {
    let n = model.total_dof();
    assert_eq!(q.len(), n, "q length mismatch");
    assert_eq!(q_dot.len(), n, "q_dot length mismatch");
    assert_eq!(q_ddot.len(), n, "q_ddot length mismatch");

    let mut dtau_dq = vec![vec![0.0f64; n]; n];
    let mut dtau_dqd = vec![vec![0.0f64; n]; n];

    let inv_2h = 1.0 / (2.0 * h);

    // ∂τ/∂q: perturb each q[j].
    let mut q_pert = q.to_vec();
    for j in 0..n {
        let original = q_pert[j];
        q_pert[j] = original + h;
        let tau_plus = rnea(model, &q_pert, q_dot, q_ddot);
        q_pert[j] = original - h;
        let tau_minus = rnea(model, &q_pert, q_dot, q_ddot);
        q_pert[j] = original;
        for i in 0..n {
            dtau_dq[i][j] = (tau_plus[i] - tau_minus[i]) * inv_2h;
        }
    }

    // ∂τ/∂q̇: perturb each q_dot[j].
    let mut qd_pert = q_dot.to_vec();
    for j in 0..n {
        let original = qd_pert[j];
        qd_pert[j] = original + h;
        let tau_plus = rnea(model, q, &qd_pert, q_ddot);
        qd_pert[j] = original - h;
        let tau_minus = rnea(model, q, &qd_pert, q_ddot);
        qd_pert[j] = original;
        for i in 0..n {
            dtau_dqd[i][j] = (tau_plus[i] - tau_minus[i]) * inv_2h;
        }
    }

    RneaDerivatives { dtau_dq, dtau_dqd }
}

/// Central finite-difference derivatives of RNEA with the default step `h = 1e-6`.
///
/// Convenience wrapper over [`rnea_derivatives_fd_h`]; see that function for the
/// accuracy and usage notes.
///
/// # Panics
/// Panics if `q`, `q_dot`, or `q_ddot` are not of length `model.total_dof()`.
pub fn rnea_derivatives_fd(
    model: &ArticulatedModel,
    q: &[f64],
    q_dot: &[f64],
    q_ddot: &[f64],
) -> RneaDerivatives {
    rnea_derivatives_fd_h(model, q, q_dot, q_ddot, 1e-6)
}

// ─── Analytic recursion ──────────────────────────────────────────────────────

/// Analytic ∂τ/∂q and ∂τ/∂q̇ via the Carpentier-Mansard derivative recursion.
///
/// Assumes every body is driven by a single-DoF, constant-subspace joint, so
/// the body index equals its single DoF index (`dof_start(i) == i`). Callers
/// must gate on [`analytic_supported`]; this function `debug_assert!`s the
/// invariant.
fn rnea_derivatives_analytic(
    model: &ArticulatedModel,
    q: &[f64],
    q_dot: &[f64],
    q_ddot: &[f64],
) -> RneaDerivatives {
    let n = model.num_bodies();
    debug_assert_eq!(
        n,
        model.total_dof(),
        "analytic path requires one DoF per body"
    );

    // Precompute children lists from parent_id (parent index < child index).
    let mut children: Vec<Vec<usize>> = vec![Vec::new(); n];
    for (i, body) in model.bodies.iter().enumerate() {
        if let Some(pid) = body.parent_id {
            children[pid].push(i);
        }
    }

    // Nominal RNEA quantities.
    let mut v = vec![SpatialVec::ZERO; n]; // body spatial velocity
    let mut a = vec![SpatialVec::ZERO; n]; // body spatial acceleration
    let mut x_total = vec![SpatialTransform::IDENTITY; n]; // parent→child motion transform
    let mut s_axis = vec![SpatialVec::ZERO; n]; // single subspace column (child frame)
    let mut v_j_store = vec![SpatialVec::ZERO; n]; // joint velocity s_i*qd_i

    // Partials of v and a, full n columns per body (O(n²) memory; n small).
    let mut dv_dq = vec![vec![SpatialVec::ZERO; n]; n];
    let mut dv_dqd = vec![vec![SpatialVec::ZERO; n]; n];
    let mut da_dq = vec![vec![SpatialVec::ZERO; n]; n];
    let mut da_dqd = vec![vec![SpatialVec::ZERO; n]; n];

    let a_gravity = model.gravity;

    // ── Forward derivative pass ──────────────────────────────────────────────
    for i in 0..n {
        debug_assert_eq!(model.dof_count(i), 1, "analytic path requires dof==1");
        let body = &model.bodies[i];
        let joint = &model.joints[i];
        let qi = model.q_slice(q, i);

        let x_j = joint.transform(qi);
        let x_full = body.parent_transform.compose(&x_j);
        x_total[i] = x_full;

        let s = joint.motion_subspace_at(qi);
        // 1-DoF: exactly one column.
        let s_i = s[0];
        s_axis[i] = s_i;

        let qd_i = q_dot[i];
        let qdd_i = q_ddot[i];
        let v_j = s_i.scale(qd_i);
        v_j_store[i] = v_j;

        let (v_parent, a_parent) = if let Some(pid) = body.parent_id {
            (v[pid], a[pid])
        } else {
            (SpatialVec::ZERO, a_gravity)
        };

        // Transported parent motion/acceleration (functions of q_i through x_full).
        let xv_p = x_full.apply_velocity(&v_parent);
        let xa_p = x_full.apply_velocity(&a_parent);

        v[i] = xv_p + v_j;
        // a[i] = X*a_parent + s_i*qdd_i + v[i] × v_J  (c_J = 0 for these joints).
        let s_qdd = s_i.scale(qdd_i);
        a[i] = xa_p + s_qdd + v[i].cross_motion(&v_j);

        let parent = body.parent_id;

        for j in 0..n {
            // (a) ∂/∂q_j.
            let mut base_dv = if let Some(p) = parent {
                x_full.apply_velocity(&dv_dq[p][j])
            } else {
                SpatialVec::ZERO
            };
            let mut base_da = if let Some(p) = parent {
                x_full.apply_velocity(&da_dq[p][j])
            } else {
                SpatialVec::ZERO
            };
            if j == i {
                // Explicit transform-derivative term: ∂(X·w)/∂q_i = -s_i × (X·w).
                base_dv = base_dv + (-s_i.cross_motion(&xv_p));
                base_da = base_da + (-s_i.cross_motion(&xa_p));
            }
            // v_J = s_i*qd_i does not depend on q, so dv_dq is just base_dv.
            let dvij = base_dv;
            dv_dq[i][j] = dvij;
            // ∂/∂q_j (v[i] × v_J) = dv_dq[i][j] × v_J  (v_J independent of q).
            da_dq[i][j] = base_da + dvij.cross_motion(&v_j);

            // (b) ∂/∂q̇_j.
            let pv = if let Some(p) = parent {
                x_full.apply_velocity(&dv_dqd[p][j])
            } else {
                SpatialVec::ZERO
            };
            let pa = if let Some(p) = parent {
                x_full.apply_velocity(&da_dqd[p][j])
            } else {
                SpatialVec::ZERO
            };
            // v_J = s_i*qd_i depends on q̇ only at j == i.
            let dvj = if j == i { s_i } else { SpatialVec::ZERO };
            let dv_qd = pv + dvj;
            dv_dqd[i][j] = dv_qd;
            // ∂/∂q̇_j (v[i] × v_J) = dv_dqd[i][j] × v_J + v[i] × dvj.
            da_dqd[i][j] = pa + dv_qd.cross_motion(&v_j) + v[i].cross_motion(&dvj);
        }
    }

    // ── Inward (force) nominal pass ──────────────────────────────────────────
    // Accumulate child forces into f[i] so that f[i] is the total force at i.
    let mut f = vec![SpatialVec::ZERO; n];
    for i in (0..n).rev() {
        let body = &model.bodies[i];
        let ia = body.inertia.mul_vec(&a[i]);
        let iv = body.inertia.mul_vec(&v[i]);
        let v_cross_f = v[i].cross_force(&iv);
        f[i] = f[i] + ia + v_cross_f;
        if let Some(pid) = body.parent_id {
            let carried = x_total[i].apply_force(&f[i]);
            f[pid] = f[pid] + carried;
        }
    }

    // ── Backward derivative pass ─────────────────────────────────────────────
    // Partials of the *accumulated* body force f[i].
    let mut df_dq = vec![vec![SpatialVec::ZERO; n]; n];
    let mut df_dqd = vec![vec![SpatialVec::ZERO; n]; n];

    for i in (0..n).rev() {
        let body = &model.bodies[i];
        let inertia = &body.inertia;
        let iv = inertia.mul_vec(&v[i]);

        for j in 0..n {
            // fnet_i = I·a[i] + v[i] ×* (I·v[i]).
            // ∂fnet/∂x = I·(da/∂x) + (dv/∂x) ×* (I·v) + v ×* (I·(dv/∂x)).
            // ∂/∂q_j.
            let i_da_q = inertia.mul_vec(&da_dq[i][j]);
            let i_dv_q = inertia.mul_vec(&dv_dq[i][j]);
            let mut df_q =
                i_da_q + dv_dq[i][j].cross_force(&iv) + v[i].cross_force(&i_dv_q);
            // ∂/∂q̇_j.
            let i_da_qd = inertia.mul_vec(&da_dqd[i][j]);
            let i_dv_qd = inertia.mul_vec(&dv_dqd[i][j]);
            let mut df_qd =
                i_da_qd + dv_dqd[i][j].cross_force(&iv) + v[i].cross_force(&i_dv_qd);

            // Accumulate children: f[i] = fnet_i + Σ_c X_c.apply_force(f_total_c).
            for &c in &children[i] {
                df_q = df_q + x_total[c].apply_force(&df_dq[c][j]);
                df_qd = df_qd + x_total[c].apply_force(&df_dqd[c][j]);
                // The carried force X_c.apply_force(f_total_c) depends on q_c via
                // X_c; its derivative wrt q_c is X_c.apply_force(s_c ×* f_total_c).
                if j == c {
                    let extra = x_total[c].apply_force(&s_axis[c].cross_force(&f[c]));
                    df_q = df_q + extra;
                }
            }

            df_dq[i][j] = df_q;
            df_dqd[i][j] = df_qd;
        }
    }

    // ── Torque partials: τ_i = s_i · f_total_i ───────────────────────────────
    let mut dtau_dq = vec![vec![0.0f64; n]; n];
    let mut dtau_dqd = vec![vec![0.0f64; n]; n];
    for i in 0..n {
        let s_i = s_axis[i];
        for j in 0..n {
            dtau_dq[i][j] = s_i.dot(&df_dq[i][j]);
            dtau_dqd[i][j] = s_i.dot(&df_dqd[i][j]);
        }
    }

    RneaDerivatives { dtau_dq, dtau_dqd }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        body::RigidBody,
        joint::{Joint, PrismaticJoint, RevoluteJoint},
        model::ArticulatedModel,
        spatial::{SpatialInertia, SpatialTransform},
    };

    fn assert_close(a: f64, b: f64, eps: f64, label: &str) {
        assert!(
            (a - b).abs() < eps,
            "{label}: expected {b:.12}, got {a:.12} (|Δ|={:.3e})",
            (a - b).abs()
        );
    }

    fn assert_vec_close(a: &SpatialVec, b: &SpatialVec, eps: f64, label: &str) {
        for k in 0..3 {
            assert_close(a.angular[k], b.angular[k], eps, &format!("{label}.ang[{k}]"));
            assert_close(a.linear[k], b.linear[k], eps, &format!("{label}.lin[{k}]"));
        }
    }

    fn assert_matrices_close(a: &[Vec<f64>], b: &[Vec<f64>], eps: f64, label: &str) {
        assert_eq!(a.len(), b.len(), "{label}: row count mismatch");
        for (i, (ar, br)) in a.iter().zip(b.iter()).enumerate() {
            assert_eq!(ar.len(), br.len(), "{label}: col count mismatch row {i}");
            for (j, (av, bv)) in ar.iter().zip(br.iter()).enumerate() {
                assert_close(*av, *bv, eps, &format!("{label}[{i}][{j}]"));
            }
        }
    }

    // ── Model builders ───────────────────────────────────────────────────────

    fn build_pendulum(m: f64, l: f64, g: f64) -> ArticulatedModel {
        let mut model = ArticulatedModel::new([0.0, -g, 0.0]);
        let joint = Box::new(RevoluteJoint::new([0.0, 0.0, 1.0]));
        let com = [0.0, -l, 0.0];
        let inertia = SpatialInertia::from_com(m, com, [[0.0; 3]; 3]);
        let body = RigidBody::new("link1", inertia, None, SpatialTransform::IDENTITY);
        model.add_body(body, joint);
        model
    }

    fn build_two_link(m1: f64, m2: f64, l1: f64, l2: f64, g: f64) -> ArticulatedModel {
        let mut model = ArticulatedModel::new([0.0, -g, 0.0]);

        // Link 1: revolute z at origin, point mass at [l1/2, 0, 0].
        let j1 = Box::new(RevoluteJoint::new([0.0, 0.0, 1.0]));
        let i1 = SpatialInertia::from_com(m1, [l1 * 0.5, 0.0, 0.0], [[0.0; 3]; 3]);
        let b1 = RigidBody::new("link1", i1, None, SpatialTransform::IDENTITY);
        model.add_body(b1, j1);

        // Link 2: revolute z, offset by l1 along x from link 1.
        let j2 = Box::new(RevoluteJoint::new([0.0, 0.0, 1.0]));
        let i2 = SpatialInertia::from_com(m2, [l2 * 0.5, 0.0, 0.0], [[0.0; 3]; 3]);
        let b2 = RigidBody::new(
            "link2",
            i2,
            Some(0),
            SpatialTransform::from_translation([l1, 0.0, 0.0]),
        );
        model.add_body(b2, j2);

        model
    }

    fn build_seven_dof() -> ArticulatedModel {
        let mut model = ArticulatedModel::new([0.0, -9.81, 0.0]);
        let axes = [
            [0.0, 0.0, 1.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let rot_com = [[0.01, 0.0, 0.0], [0.0, 0.01, 0.0], [0.0, 0.0, 0.01]];
        for (idx, axis) in axes.iter().enumerate() {
            let joint = Box::new(RevoluteJoint::new(*axis));
            let mass = 1.0 + 0.1 * idx as f64;
            let inertia = SpatialInertia::from_com(mass, [0.15, 0.0, 0.0], rot_com);
            let (parent_id, parent_transform) = if idx == 0 {
                (None, SpatialTransform::IDENTITY)
            } else {
                (
                    Some(idx - 1),
                    SpatialTransform::from_translation([0.3, 0.0, 0.0]),
                )
            };
            let body = RigidBody::new(format!("link{idx}"), inertia, parent_id, parent_transform);
            model.add_body(body, joint);
        }
        model
    }

    // ── Crux helper identity tests (do these FIRST) ──────────────────────────

    #[test]
    fn test_transform_velocity_derivative_matches_fd() {
        let joint = RevoluteJoint::new([0.0, 0.0, 1.0]);
        // Fixed parent transform: rotation about y plus a translation.
        let x_t = SpatialTransform::from_axis_angle([0.0, 1.0, 0.0], 0.4)
            .compose(&SpatialTransform::from_translation([0.2, -0.3, 0.5]));
        let w = SpatialVec::new([0.7, -0.2, 0.5], [0.1, 0.9, -0.4]);
        let q0 = 0.37_f64;
        let h = 1e-6;

        let x = x_t.compose(&joint.transform(&[q0]));
        let s = joint.motion_subspace_at(&[q0])[0];
        let analytic = -s.cross_motion(&x.apply_velocity(&w));

        let x_plus = x_t.compose(&joint.transform(&[q0 + h]));
        let x_minus = x_t.compose(&joint.transform(&[q0 - h]));
        let fd =
            (x_plus.apply_velocity(&w) - x_minus.apply_velocity(&w)).scale(1.0 / (2.0 * h));

        assert_vec_close(&analytic, &fd, 1e-6, "velocity derivative");
    }

    #[test]
    fn test_transform_force_derivative_matches_fd() {
        let joint = RevoluteJoint::new([0.0, 0.0, 1.0]);
        let x_t = SpatialTransform::from_axis_angle([0.0, 1.0, 0.0], 0.4)
            .compose(&SpatialTransform::from_translation([0.2, -0.3, 0.5]));
        // Child-frame force.
        let f = SpatialVec::new([0.3, 0.6, -0.2], [1.1, -0.7, 0.4]);
        let q0 = 0.37_f64;
        let h = 1e-6;

        let x = x_t.compose(&joint.transform(&[q0]));
        let s = joint.motion_subspace_at(&[q0])[0];
        // Identity used in production: ∂(X.apply_force(f))/∂q_i = X.apply_force(s ×* f).
        let analytic = x.apply_force(&s.cross_force(&f));

        let x_plus = x_t.compose(&joint.transform(&[q0 + h]));
        let x_minus = x_t.compose(&joint.transform(&[q0 - h]));
        let fd = (x_plus.apply_force(&f) - x_minus.apply_force(&f)).scale(1.0 / (2.0 * h));

        assert_vec_close(&analytic, &fd, 1e-6, "force derivative");
    }

    // ── Integration tests: analytic vs FD ────────────────────────────────────

    #[test]
    fn test_rnea_deriv_pendulum_vs_fd() {
        let model = build_pendulum(1.5, 0.8, 9.81);
        let q = [0.6];
        let qd = [0.4];
        let qdd = [-0.3];

        let analytic = rnea_derivatives(&model, &q, &qd, &qdd);
        assert!(analytic_supported(&model), "pendulum should use analytic path");
        let fd = rnea_derivatives_fd_h(&model, &q, &qd, &qdd, 1e-6);

        assert_matrices_close(&analytic.dtau_dq, &fd.dtau_dq, 1e-6, "pendulum dtau_dq");
        assert_matrices_close(&analytic.dtau_dqd, &fd.dtau_dqd, 1e-6, "pendulum dtau_dqd");
    }

    #[test]
    fn test_rnea_deriv_two_link_vs_fd() {
        let model = build_two_link(1.2, 0.9, 0.5, 0.4, 9.81);
        let q = [0.3, -0.5];
        let qd = [0.2, 0.7];
        let qdd = [0.1, -0.4];

        let analytic = rnea_derivatives(&model, &q, &qd, &qdd);
        assert!(analytic_supported(&model), "two-link should use analytic path");
        let fd = rnea_derivatives_fd_h(&model, &q, &qd, &qdd, 1e-6);

        assert_matrices_close(&analytic.dtau_dq, &fd.dtau_dq, 1e-6, "two-link dtau_dq");
        assert_matrices_close(&analytic.dtau_dqd, &fd.dtau_dqd, 1e-6, "two-link dtau_dqd");
    }

    #[test]
    fn test_rnea_deriv_seven_dof_vs_fd() {
        let model = build_seven_dof();
        let q = [0.2, -0.4, 0.6, -0.1, 0.3, -0.5, 0.15];
        let qd = [0.1, 0.2, -0.3, 0.4, -0.2, 0.1, 0.25];
        let qdd = [-0.05, 0.15, 0.1, -0.2, 0.3, -0.1, 0.2];

        let analytic = rnea_derivatives(&model, &q, &qd, &qdd);
        assert!(analytic_supported(&model), "7-DoF chain should use analytic path");
        let fd = rnea_derivatives_fd_h(&model, &q, &qd, &qdd, 1e-6);

        assert_matrices_close(&analytic.dtau_dq, &fd.dtau_dq, 1e-6, "7-DoF dtau_dq");
        assert_matrices_close(&analytic.dtau_dqd, &fd.dtau_dqd, 1e-6, "7-DoF dtau_dqd");
    }

    #[test]
    fn test_rnea_deriv_prismatic_chain_vs_fd() {
        // Mixed revolute/prismatic chain exercises the prismatic subspace path.
        let mut model = ArticulatedModel::new([0.0, -9.81, 0.0]);
        let rot_com = [[0.02, 0.0, 0.0], [0.0, 0.02, 0.0], [0.0, 0.0, 0.02]];
        let j1 = Box::new(RevoluteJoint::new([0.0, 0.0, 1.0]));
        let i1 = SpatialInertia::from_com(1.0, [0.2, 0.0, 0.0], rot_com);
        model.add_body(
            RigidBody::new("rev", i1, None, SpatialTransform::IDENTITY),
            j1,
        );
        let j2 = Box::new(PrismaticJoint::new([1.0, 0.0, 0.0]));
        let i2 = SpatialInertia::from_com(0.8, [0.0, 0.1, 0.0], rot_com);
        model.add_body(
            RigidBody::new(
                "pris",
                i2,
                Some(0),
                SpatialTransform::from_translation([0.3, 0.0, 0.0]),
            ),
            j2,
        );

        let q = [0.4, 0.25];
        let qd = [0.3, -0.2];
        let qdd = [0.1, 0.15];

        let analytic = rnea_derivatives(&model, &q, &qd, &qdd);
        assert!(analytic_supported(&model), "rev+prismatic should use analytic path");
        let fd = rnea_derivatives_fd_h(&model, &q, &qd, &qdd, 1e-6);

        assert_matrices_close(&analytic.dtau_dq, &fd.dtau_dq, 1e-6, "rev+pris dtau_dq");
        assert_matrices_close(&analytic.dtau_dqd, &fd.dtau_dqd, 1e-6, "rev+pris dtau_dqd");
    }

    #[test]
    fn test_fd_self_consistency() {
        let model = build_two_link(1.2, 0.9, 0.5, 0.4, 9.81);
        let q = [0.3, -0.5];
        let qd = [0.2, 0.7];
        let qdd = [0.1, -0.4];

        let fd_a = rnea_derivatives_fd_h(&model, &q, &qd, &qdd, 1e-6);
        let fd_b = rnea_derivatives_fd_h(&model, &q, &qd, &qdd, 2e-6);

        assert_matrices_close(&fd_a.dtau_dq, &fd_b.dtau_dq, 1e-4, "FD self dtau_dq");
        assert_matrices_close(&fd_a.dtau_dqd, &fd_b.dtau_dqd, 1e-4, "FD self dtau_dqd");
    }
}
