// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

use oxiphysics_articulated::{
    aba::aba,
    body::RigidBody,
    joint::{Joint, RevoluteJoint, UniversalJoint},
    model::ArticulatedModel,
    spatial::{SpatialInertia, SpatialTransform},
};
use std::f64::consts::PI;

fn assert_mat3_approx(a: [[f64; 3]; 3], b: [[f64; 3]; 3], tol: f64, msg: &str) {
    for i in 0..3 {
        for j in 0..3 {
            let diff = (a[i][j] - b[i][j]).abs();
            assert!(
                diff < tol,
                "{}: a[{}][{}]={:.2e} b[{}][{}]={:.2e} diff={:.2e}",
                msg,
                i,
                j,
                a[i][j],
                i,
                j,
                b[i][j],
                diff
            );
        }
    }
}

#[test]
fn test_universal_joint_dof_and_transform_at_zero() {
    let j = UniversalJoint::new([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
    assert_eq!(j.dof(), 2);
    let x = j.transform(&[0.0, 0.0]);
    assert_mat3_approx(
        x.rot,
        SpatialTransform::IDENTITY.rot,
        1e-14,
        "zero-q transform should be identity",
    );
    for component in &x.trans {
        assert!(component.abs() < 1e-14, "translation should be zero");
    }
}

#[test]
fn test_universal_locked_first_axis_reduces_to_revolute() {
    let universal = UniversalJoint::new([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
    let revolute = RevoluteJoint::new([0.0, 1.0, 0.0]);
    for k in 0..8 {
        let angle = k as f64 * PI / 4.0;
        let xu = universal.transform(&[0.0, angle]);
        let xr = revolute.transform(&[angle]);
        assert_mat3_approx(xu.rot, xr.rot, 1e-12, &format!("angle={:.3}", angle));
    }
}

#[test]
fn test_universal_motion_subspace_q_dependence() {
    let j = UniversalJoint::new([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
    let s0 = j.motion_subspace(); // at q=[0,0]
    let sq = j.motion_subspace_at(&[0.3, 0.7]); // at non-zero q
    // column 0 should differ between q-independent and q-dependent
    let col0_s0 = s0[0].angular;
    let col0_sq = sq[0].angular;
    let diff = ((col0_s0[0] - col0_sq[0]).powi(2)
        + (col0_s0[1] - col0_sq[1]).powi(2)
        + (col0_s0[2] - col0_sq[2]).powi(2))
    .sqrt();
    assert!(
        diff > 0.05,
        "motion_subspace_at should differ from motion_subspace() at non-zero q, got diff={:.4}",
        diff
    );
}

#[test]
fn test_universal_dynamics_gravity_sign() {
    // 1-body universal pendulum, mass 1 kg at (0, 0, -1) below the joint.
    let g = 9.81_f64;
    let mut model = ArticulatedModel::new([0.0, 0.0, -g]);
    let inertia = SpatialInertia::from_com(1.0, [0.0, 0.0, -1.0], [[0.0; 3]; 3]);
    let body = RigidBody::new("link1", inertia, None, SpatialTransform::IDENTITY);
    model.add_body(
        body,
        Box::new(UniversalJoint::new([1.0, 0.0, 0.0], [0.0, 1.0, 0.0])),
    );
    let q = vec![0.0, 0.0];
    let qd = vec![0.0, 0.0];
    let tau = vec![0.0, 0.0];
    let qdd = aba(&model, &q, &qd, &tau);
    assert_eq!(qdd.len(), 2);
    // At q=[0,0], COM is directly below. With mass at (0,0,-1) and gravity in -z,
    // the torques about both axes are zero (symmetric configuration). So qdd~[0,0].
    for (k, &val) in qdd.iter().enumerate() {
        assert!(val.is_finite(), "qdd[{}] must be finite", k);
    }
    // Now test with first joint rotated so mass is displaced:
    let q2 = vec![PI / 6.0, 0.0]; // 30 degrees about x
    let qdd2 = aba(&model, &q2, &qd, &tau);
    assert!(qdd2[0].is_finite() && qdd2[1].is_finite());
    // With mass displaced, at least one qdd should be non-zero
    let any_nonzero = qdd2[0].abs() > 1e-6 || qdd2[1].abs() > 1e-6;
    assert!(
        any_nonzero,
        "qdd should be non-zero when joint is at non-equilibrium, got {:?}",
        qdd2
    );
}

#[test]
fn test_universal_rnea_aba_round_trip() {
    use oxiphysics_articulated::rnea::rnea;
    let g = 9.81_f64;
    let mut model = ArticulatedModel::new([0.0, 0.0, -g]);
    let inertia = SpatialInertia::from_com(
        2.0,
        [0.0, 0.0, -0.5],
        [[0.1, 0.0, 0.0], [0.0, 0.1, 0.0], [0.0, 0.0, 0.05]],
    );
    let body = RigidBody::new("link1", inertia, None, SpatialTransform::IDENTITY);
    model.add_body(
        body,
        Box::new(UniversalJoint::new([1.0, 0.0, 0.0], [0.0, 1.0, 0.0])),
    );
    let q = vec![0.2, -0.3];
    let qd = vec![0.5, -0.1];
    let qdd_target = vec![1.0, -2.0];
    // RNEA: given qdd_target, compute τ
    let tau = rnea(&model, &q, &qd, &qdd_target);
    // ABA: given τ, recover qdd
    let qdd_recovered = aba(&model, &q, &qd, &tau);
    assert_eq!(qdd_recovered.len(), 2);
    for (k, (&target, &recovered)) in qdd_target.iter().zip(qdd_recovered.iter()).enumerate() {
        let diff = (recovered - target).abs();
        assert!(
            diff < 1e-8,
            "qdd[{}] round-trip: target={:.6} recovered={:.6} diff={:.2e}",
            k,
            target,
            recovered,
            diff
        );
    }
}

// ─── Kinematic consistency of transform, S(q) and c_J ────────────────────────

fn mat3_mul(a: [[f64; 3]; 3], b: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut out = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                out[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    out
}

fn mat3_t(a: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut out = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            out[i][j] = a[j][i];
        }
    }
    out
}

/// Body-frame angular velocity `vee(Rᵀ·Ṙ)` of the joint pose along `q(t) = q + t·qd`.
fn body_omega_fd(j: &UniversalJoint, q: [f64; 2], qd: [f64; 2], h: f64) -> [f64; 3] {
    let r0 = j.transform(&q).rot;
    let rp = j.transform(&[q[0] + h * qd[0], q[1] + h * qd[1]]).rot;
    let rm = j.transform(&[q[0] - h * qd[0], q[1] - h * qd[1]]).rot;
    let mut dr = [[0.0; 3]; 3];
    for i in 0..3 {
        for k in 0..3 {
            dr[i][k] = (rp[i][k] - rm[i][k]) / (2.0 * h);
        }
    }
    let w = mat3_mul(mat3_t(r0), dr);
    [
        0.5 * (w[2][1] - w[1][2]),
        0.5 * (w[0][2] - w[2][0]),
        0.5 * (w[1][0] - w[0][1]),
    ]
}

#[test]
fn test_universal_motion_subspace_matches_transform_derivative() {
    let j = UniversalJoint::new([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
    let cases = [
        ([0.4, -0.9], [0.7, -1.1]),
        ([-1.3, 0.6], [0.2, 0.9]),
        ([2.1, 1.4], [-1.5, 0.3]),
    ];
    for (q, qd) in cases {
        let omega_fd = body_omega_fd(&j, q, qd, 1e-6);
        let s = j.motion_subspace_at(&q);
        let omega_s = [
            s[0].angular[0] * qd[0] + s[1].angular[0] * qd[1],
            s[0].angular[1] * qd[0] + s[1].angular[1] * qd[1],
            s[0].angular[2] * qd[0] + s[1].angular[2] * qd[1],
        ];
        for k in 0..3 {
            assert!(
                (omega_fd[k] - omega_s[k]).abs() < 1e-7,
                "q={q:?} qd={qd:?}: ω_fd[{k}]={:.10} vs S·qd={:.10}",
                omega_fd[k],
                omega_s[k]
            );
        }
        for col in &s {
            assert_eq!(col.linear, [0.0; 3], "universal S has no linear part");
        }
    }
}

#[test]
fn test_universal_coriolis_matches_motion_subspace_derivative() {
    let j = UniversalJoint::new([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
    let q = [0.4, -0.9];
    let qd = [0.7, -1.1];
    let h = 1e-6;
    let sp = j.motion_subspace_at(&[q[0] + h * qd[0], q[1] + h * qd[1]]);
    let sm = j.motion_subspace_at(&[q[0] - h * qd[0], q[1] - h * qd[1]]);
    // Ṡ·q̇ by central differences.
    let mut sdot_qd = [0.0; 3];
    for (col, &qdk) in qd.iter().enumerate() {
        for (k, acc) in sdot_qd.iter_mut().enumerate() {
            *acc += (sp[col].angular[k] - sm[col].angular[k]) / (2.0 * h) * qdk;
        }
    }
    let c = j.coriolis(&q, &qd);
    for (k, fd) in sdot_qd.iter().enumerate() {
        assert!(
            (c.angular[k] - fd).abs() < 1e-7,
            "c_J[{k}]={:.10} vs FD Ṡq̇={:.10}",
            c.angular[k],
            fd
        );
    }
    assert_eq!(c.linear, [0.0; 3]);
    // The bias term is genuinely non-zero here.
    let norm = (c.angular[0].powi(2) + c.angular[1].powi(2) + c.angular[2].powi(2)).sqrt();
    assert!(norm > 0.1, "c_J should be non-zero, |c_J|={norm:.3e}");
}

/// Universal joint ≡ revolute(axis1) → massless intermediate link → revolute(axis2).
#[test]
fn test_universal_dynamics_equal_two_revolute_chain() {
    use oxiphysics_articulated::rnea::rnea;
    let g = 9.81_f64;
    let i_com = [[0.12, 0.01, -0.02], [0.01, 0.08, 0.03], [-0.02, 0.03, 0.05]];
    let inertia = SpatialInertia::from_com(1.7, [0.2, -0.1, -0.6], i_com);

    let mut universal = ArticulatedModel::new([0.0, 0.0, -g]);
    universal.add_body(
        RigidBody::new("link", inertia, None, SpatialTransform::IDENTITY),
        Box::new(UniversalJoint::new([1.0, 0.0, 0.0], [0.0, 1.0, 0.0])),
    );

    let mut chain = ArticulatedModel::new([0.0, 0.0, -g]);
    let massless = SpatialInertia::from_com(0.0, [0.0; 3], [[0.0; 3]; 3]);
    chain.add_body(
        RigidBody::new("cross", massless, None, SpatialTransform::IDENTITY),
        Box::new(RevoluteJoint::new([1.0, 0.0, 0.0])),
    );
    chain.add_body(
        RigidBody::new("link", inertia, Some(0), SpatialTransform::IDENTITY),
        Box::new(RevoluteJoint::new([0.0, 1.0, 0.0])),
    );

    let q = [0.4, -0.9];
    let qd = [0.7, -1.1];
    let qdd = [0.3, 1.2];
    let tau_u = rnea(&universal, &q, &qd, &qdd);
    let tau_c = rnea(&chain, &q, &qd, &qdd);
    for k in 0..2 {
        assert!(
            (tau_u[k] - tau_c[k]).abs() < 1e-10,
            "tau[{k}]: universal={:.12} chain={:.12}",
            tau_u[k],
            tau_c[k]
        );
    }
    // Forward dynamics agree as well.
    let tau = [0.5, -0.25];
    let qdd_u = aba(&universal, &q, &qd, &tau);
    let qdd_c = aba(&chain, &q, &qd, &tau);
    for k in 0..2 {
        assert!(
            (qdd_u[k] - qdd_c[k]).abs() < 1e-9,
            "qdd[{k}]: universal={:.12} chain={:.12}",
            qdd_u[k],
            qdd_c[k]
        );
    }
}
