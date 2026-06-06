// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Articulated Body Algorithm (ABA) for forward dynamics.
//!
//! Given joint positions `q`, velocities `q̇`, and applied torques `τ`,
//! computes joint accelerations `q̈` that would result from those inputs.
//!
//! Reference: Featherstone 2008 "Rigid Body Dynamics Algorithms", Algorithm 7.4.
//!
//! # Algorithm overview
//!
//! **Pass 1 (outward):** Compute velocity terms and transform bodies.
//! For each body `i` (root to leaves):
//! - Compute `X_J`, compose full transform `X = X_T ∘ X_J`.
//! - Propagate velocity: `v[i] = X · v[parent] + S · q̇_i`.
//! - Compute bias acceleration: `c[i] = v[i] × (S · q̇_i)  +  c_J`.
//!
//! **Pass 2 (inward):** Compute articulated inertias and bias forces.
//! For each body `i` (leaves to root):
//! - Initialise articulated inertia: `I_A[i] = I[i]`.
//! - For each child `j`: back-project child's articulated inertia to parent.
//! - Compute Articulated inertia Schur complement update:
//!   ```text
//!   U = I_A[child] · S
//!   D = S^T · U     (scalar DOF-×-DOF matrix, scalar for 1-DOF joints)
//!   I_A[parent] += X^* · (I_A[child] − U · D^{−1} · U^T) · X^{*T}
//!   p_A[parent] += X^* · (p_A[child] + I_A[child] · c_child + U · D^{−1} · (τ_child − S^T · p_A[child]))
//!   ```
//!
//! **Pass 3 (outward):** Compute accelerations.
//! For each body `i` (root to leaves):
//! - `a_hat = X · a[parent] + c[i]`
//! - `q̈_i = D^{−1} · (τ_i − S^T · (p_A[i] + I_A[i] · a_hat))`
//! - `a[i] = a_hat + S · q̈_i`

use crate::{
    model::ArticulatedModel,
    spatial::{SpatialInertia6x6, SpatialTransform, SpatialVec},
};

// ─── Internal per-body workspace ────────────────────────────────────────────

/// Per-body workspace for ABA intermediate values.
struct AbaBody {
    /// Full joint-to-parent transform (X_T ∘ X_J).
    x_full: SpatialTransform,
    /// Body spatial velocity (in body frame).
    v: SpatialVec,
    /// Bias acceleration c[i] = v × v_J + c_J (no parent contribution yet).
    c: SpatialVec,
    /// Articulated inertia (6×6, starts as body's own inertia, then updated).
    i_a: SpatialInertia6x6,
    /// Articulated bias force (6-vector, starts as v ×* (I·v), then updated).
    p_a: SpatialVec,
    /// U = I_A · S (one column per DOF; here we assume max 6 DOF per joint).
    u: Vec<SpatialVec>,
    /// D = S^T · U (DOF×DOF matrix, stored as flat row-major vec).
    d: Vec<f64>,
    /// D^{−1} (inverse of D; for 1-DOF joints this is a scalar).
    d_inv: Vec<f64>,
    /// τ_i − S^T · p_A (bias-corrected torque).
    tau_bias: Vec<f64>,
    /// Body acceleration (accumulated in pass 3).
    a: SpatialVec,
}

impl AbaBody {
    fn new(n_dof: usize) -> Self {
        Self {
            x_full: SpatialTransform::IDENTITY,
            v: SpatialVec::ZERO,
            c: SpatialVec::ZERO,
            i_a: SpatialInertia6x6 {
                data: [[0.0; 6]; 6],
            },
            p_a: SpatialVec::ZERO,
            u: vec![SpatialVec::ZERO; n_dof],
            d: vec![0.0; n_dof * n_dof],
            d_inv: vec![0.0; n_dof * n_dof],
            tau_bias: vec![0.0; n_dof],
            a: SpatialVec::ZERO,
        }
    }
}

/// Compute joint accelerations via ABA (forward dynamics).
///
/// # Arguments
/// - `model`: Articulated-body model (kinematic tree, topological order).
/// - `q`: Joint positions; length must equal `model.total_dof()`.
/// - `q_dot`: Joint velocities; same length as `q`.
/// - `tau`: Applied joint torques; same length as `q`.
///
/// # Returns
/// Joint accelerations `q̈` (length = `model.total_dof()`).
pub fn aba(model: &ArticulatedModel, q: &[f64], q_dot: &[f64], tau: &[f64]) -> Vec<f64> {
    let n = model.num_bodies();
    assert_eq!(
        q.len(),
        model.total_dof(),
        "q length mismatch: got {}, expected {}",
        q.len(),
        model.total_dof()
    );
    assert_eq!(q_dot.len(), model.total_dof(), "q_dot length mismatch");
    assert_eq!(tau.len(), model.total_dof(), "tau length mismatch");

    // Allocate workspace (one entry per body)
    let mut wb: Vec<AbaBody> = (0..n)
        .map(|i| AbaBody::new(model.joints[i].dof()))
        .collect();

    let a_gravity = model.gravity;

    // ── Pass 1: Outward — velocities and bias accelerations ──────────────────

    for i in 0..n {
        let body = &model.bodies[i];
        let joint = &model.joints[i];
        let qi = model.q_slice(q, i);
        let q_dot_i = model.q_slice(q_dot, i);

        let x_j = joint.transform(qi);
        let x_t = &body.parent_transform;
        let x_full = x_t.compose(&x_j);
        wb[i].x_full = x_full;

        let s = joint.motion_subspace_at(qi);

        // Joint velocity v_J = S · q̇_i
        let v_j = subspace_velocity(&s, q_dot_i);

        // Coriolis c_J
        let c_j = joint.coriolis(qi, q_dot_i);

        let v_parent = if let Some(pid) = body.parent_id {
            wb[pid].v
        } else {
            SpatialVec::ZERO
        };

        // Body velocity: v[i] = X * v_parent + v_J
        let v_parent_in_child = x_full.apply_velocity(&v_parent);
        wb[i].v = v_parent_in_child + v_j;

        // Bias acceleration: c[i] = v[i] × v_J + c_J
        let v_cross_vj = wb[i].v.cross_motion(&v_j);
        wb[i].c = v_cross_vj + c_j;

        // Initialise articulated inertia with body's own inertia
        wb[i].i_a = SpatialInertia6x6::from_packed(&body.inertia);

        // Initialise bias force: p_A[i] = v[i] ×* (I[i] * v[i])
        let iv = body.inertia.mul_vec(&wb[i].v);
        wb[i].p_a = wb[i].v.cross_force(&iv);
    }

    // ── Pass 2: Inward — articulated inertias and bias forces ────────────────

    for i in (0..n).rev() {
        let body = &model.bodies[i];
        let joint = &model.joints[i];
        let tau_i = model.q_slice(tau, i);
        let qi = model.q_slice(q, i);
        let s = joint.motion_subspace_at(qi);
        let dof = joint.dof();

        if dof > 0 {
            // U_k = I_A[i] · S_k  (for each column k of motion subspace)
            let u: Vec<SpatialVec> = s.iter().map(|sk| wb[i].i_a.mul_vec6(sk)).collect();

            // D[k,l] = S_k^T · U_l  (DOF×DOF matrix)
            let mut d = vec![0.0f64; dof * dof];
            for k in 0..dof {
                for l in 0..dof {
                    d[k * dof + l] = s[k].dot(&u[l]);
                }
            }

            // D^{-1}: for 1-DOF this is 1/d[0,0]
            let d_inv = invert_dof_matrix(&d, dof);

            // τ_bias_k = τ_k − S_k^T · p_A[i]
            let tau_bias: Vec<f64> = (0..dof).map(|k| tau_i[k] - s[k].dot(&wb[i].p_a)).collect();

            // Store for pass 3
            wb[i].u = u.clone();
            wb[i].d = d.clone();
            wb[i].d_inv = d_inv.clone();
            wb[i].tau_bias = tau_bias.clone();

            if let Some(pid) = body.parent_id {
                // Back-project articulated inertia to parent
                // I_A_hat = I_A[i] - Σ_{k,l} U_k · D_inv[k,l] · U_l^T
                let mut i_a_hat = wb[i].i_a;
                for k in 0..dof {
                    for l in 0..dof {
                        let scale = d_inv[k * dof + l];
                        // rank-1 subtraction via outer product U_k ⊗ U_l
                        i_a_hat = i_a_hat.sub_rank1_pair(&u[k], &u[l], scale);
                    }
                }

                // Transform I_A_hat to parent frame and add.
                // x_full is the parent→child motion transform (X_pc).
                // Formula: I_parent = X_pc^T * I_hat * X_pc.
                // transform_to_parent takes X_pc and computes X^T * I * X.
                let i_a_parent = i_a_hat.transform_to_parent(&wb[i].x_full);
                wb[pid].i_a = wb[pid].i_a.add(&i_a_parent);

                // Back-project bias force to parent
                // p_a_contribution = p_A[i] + I_A_hat · c[i] + U · D^{-1} · tau_bias
                // Must use the Schur complement i_a_hat (not full I_A):
                // Featherstone eq. 7.42: pA_parent += X^*(pA + Ia*c + U*D^-1*u)
                // where Ia = I_A - U·D^{-1}·U^T (same Schur complement used for inertia).
                let ia_c = i_a_hat.mul_vec6(&wb[i].c);
                let mut u_d_inv_tau = SpatialVec::ZERO;
                for k in 0..dof {
                    let mut coeff = 0.0;
                    for l in 0..dof {
                        coeff += d_inv[k * dof + l] * tau_bias[l];
                    }
                    u_d_inv_tau = u_d_inv_tau + u[k].scale(coeff);
                }
                let p_contribution = wb[i].p_a + ia_c + u_d_inv_tau;
                let p_in_parent = wb[i].x_full.apply_force(&p_contribution);
                wb[pid].p_a = wb[pid].p_a + p_in_parent;
            }
        } else {
            // Fixed joint: propagate inertia and bias directly (no DOF reduction)
            if let Some(pid) = body.parent_id {
                let i_a_parent = wb[i].i_a.transform_to_parent(&wb[i].x_full);
                wb[pid].i_a = wb[pid].i_a.add(&i_a_parent);

                let ia_c = wb[i].i_a.mul_vec6(&wb[i].c);
                let p_contribution = wb[i].p_a + ia_c;
                let p_in_parent = wb[i].x_full.apply_force(&p_contribution);
                wb[pid].p_a = wb[pid].p_a + p_in_parent;
            }
        }
    }

    // ── Pass 3: Outward — accelerations ──────────────────────────────────────

    let mut q_ddot = vec![0.0f64; model.total_dof()];

    for i in 0..n {
        let body = &model.bodies[i];
        let joint = &model.joints[i];
        let qi = model.q_slice(q, i);
        let s = joint.motion_subspace_at(qi);
        let dof = joint.dof();

        // a_hat = X * a_parent + c[i]
        let a_parent = if let Some(pid) = body.parent_id {
            wb[pid].a
        } else {
            a_gravity // gravity as fictitious parent acceleration
        };
        let a_parent_in_child = wb[i].x_full.apply_velocity(&a_parent);
        let a_hat = a_parent_in_child + wb[i].c;

        if dof > 0 {
            // Featherstone eq. 7.56: q̈_i = D_i^{-1} · (u_i − U_i^T · a_hat)
            //
            // where u_i = τ_i − S_i^T · p_A[i]  (= tau_bias, stored in Pass 2)
            // and U_i = I_A[i] · S_i              (= wb[i].u, stored in Pass 2)
            //
            // Note: tau_bias already contains the −S^T·p_A term; do NOT
            // subtract S^T·p_A again here or the bias force is double-counted.
            let ia_a_hat = wb[i].i_a.mul_vec6(&a_hat);

            let mut rhs = vec![0.0f64; dof];
            for l in 0..dof {
                // rhs = u_i - U_i^T · a_hat = tau_bias - S^T·(I_A·a_hat)
                rhs[l] = wb[i].tau_bias[l] - s[l].dot(&ia_a_hat);
            }

            let start = model.dof_start(i);
            let mut q_ddot_i = vec![0.0f64; dof];
            for k in 0..dof {
                for (l, rhs_l) in rhs.iter().enumerate() {
                    q_ddot_i[k] += wb[i].d_inv[k * dof + l] * rhs_l;
                }
                q_ddot[start + k] = q_ddot_i[k];
            }

            // a[i] = a_hat + S · q̈_i
            let s_q_ddot = subspace_velocity(&s, &q_ddot_i);
            wb[i].a = a_hat + s_q_ddot;
        } else {
            // Fixed joint: acceleration propagates without DOF contribution
            wb[i].a = a_hat;
        }
    }

    q_ddot
}

// ─── Helpers ────────────────────────────────────────────────────────────────

fn subspace_velocity(s: &[SpatialVec], v: &[f64]) -> SpatialVec {
    let mut result = SpatialVec::ZERO;
    for (sk, &vk) in s.iter().zip(v.iter()) {
        result = result + sk.scale(vk);
    }
    result
}

/// Invert a small DOF×DOF matrix (row-major).
///
/// For 1×1: trivial reciprocal.
/// For 2×2: analytical 2×2 inverse.
/// For larger: Gauss-Jordan elimination.
fn invert_dof_matrix(d: &[f64], dof: usize) -> Vec<f64> {
    match dof {
        0 => vec![],
        1 => {
            let d00 = d[0];
            let inv = if d00.abs() > 1e-30 { 1.0 / d00 } else { 0.0 };
            vec![inv]
        }
        2 => {
            let det = d[0] * d[3] - d[1] * d[2];
            let det_inv = if det.abs() > 1e-30 { 1.0 / det } else { 0.0 };
            vec![
                d[3] * det_inv,
                -d[1] * det_inv,
                -d[2] * det_inv,
                d[0] * det_inv,
            ]
        }
        n => {
            // Gauss-Jordan for n×n
            let mut a = vec![0.0f64; n * n];
            let mut inv = vec![0.0f64; n * n];
            a.copy_from_slice(d);
            // Initialise identity
            for i in 0..n {
                inv[i * n + i] = 1.0;
            }
            for col in 0..n {
                // Find pivot
                let mut pivot_row = col;
                let mut pivot_val = a[col * n + col].abs();
                for row in (col + 1)..n {
                    let v = a[row * n + col].abs();
                    if v > pivot_val {
                        pivot_val = v;
                        pivot_row = row;
                    }
                }
                if pivot_row != col {
                    for j in 0..n {
                        a.swap(col * n + j, pivot_row * n + j);
                        inv.swap(col * n + j, pivot_row * n + j);
                    }
                }
                let diag = a[col * n + col];
                if diag.abs() < 1e-30 {
                    continue;
                }
                let scale = 1.0 / diag;
                for j in 0..n {
                    a[col * n + j] *= scale;
                    inv[col * n + j] *= scale;
                }
                for row in 0..n {
                    if row == col {
                        continue;
                    }
                    let factor = a[row * n + col];
                    for j in 0..n {
                        let a_col_j = a[col * n + j];
                        let inv_col_j = inv[col * n + j];
                        a[row * n + j] -= factor * a_col_j;
                        inv[row * n + j] -= factor * inv_col_j;
                    }
                }
            }
            inv
        }
    }
}

// ─── Extension for SpatialInertia6x6 needed by ABA ───────────────────────────

impl SpatialInertia6x6 {
    /// Subtract outer product of two *different* spatial vectors (for ABA):
    /// `I -= u ⊗ w * scale`
    pub(crate) fn sub_rank1_pair(&self, u: &SpatialVec, w: &SpatialVec, scale: f64) -> Self {
        let uu = [
            u.angular[0],
            u.angular[1],
            u.angular[2],
            u.linear[0],
            u.linear[1],
            u.linear[2],
        ];
        let ww = [
            w.angular[0],
            w.angular[1],
            w.angular[2],
            w.linear[0],
            w.linear[1],
            w.linear[2],
        ];
        let mut d = self.data;
        for i in 0..6 {
            for j in 0..6 {
                d[i][j] -= uu[i] * ww[j] * scale;
            }
        }
        Self { data: d }
    }
}

// ─── Tests ────────────────────────────────────────────────────────────────────

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

    fn build_pendulum(m: f64, l: f64, g: f64) -> ArticulatedModel {
        let mut model = ArticulatedModel::new([0.0, -g, 0.0]);
        let joint = Box::new(RevoluteJoint::new([0.0, 0.0, 1.0]));
        let com = [0.0, -l, 0.0];
        let inertia = SpatialInertia::from_com(m, com, [[0.0; 3]; 3]);
        let body = RigidBody::new("link1", inertia, None, SpatialTransform::IDENTITY);
        model.add_body(body, joint);
        model
    }

    #[test]
    fn test_aba_pendulum_horizontal() {
        // Single pendulum at q=π/2 (horizontal), zero velocity, no torque.
        // Expected q_ddot = -g/L (for a point mass: θ̈ = -g/L·sin(π/2) = -g/L).
        let (m, l, g) = (1.0, 1.0, 9.81);
        let model = build_pendulum(m, l, g);

        let q_ddot = aba(&model, &[std::f64::consts::FRAC_PI_2], &[0.0], &[0.0]);
        let expected = -g / l;
        assert!(
            approx_eq(q_ddot[0], expected, 1e-6),
            "Pendulum ABA at horizontal: expected q_ddot={expected:.6}, got {:.6}",
            q_ddot[0]
        );
    }

    #[test]
    fn test_aba_hanging_no_torque() {
        // At q=0 (hanging), no torque. Expected q_ddot=0 (stable equilibrium).
        let (m, l, g) = (1.0, 1.0, 9.81);
        let model = build_pendulum(m, l, g);

        let q_ddot = aba(&model, &[0.0], &[0.0], &[0.0]);
        assert!(
            approx_eq(q_ddot[0], 0.0, 1e-8),
            "At q=0 (hanging), q_ddot should be 0, got {:.8}",
            q_ddot[0]
        );
    }
}
