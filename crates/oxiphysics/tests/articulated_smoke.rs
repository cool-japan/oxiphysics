// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Integration smoke test: 3-link revolute pendulum via ABA forward dynamics.
//!
//! Verifies that the `oxiphysics::articulated` re-export is wired correctly and
//! that one ABA step under gravity produces finite, physically non-trivial
//! joint accelerations.

use oxiphysics::articulated::{
    aba::aba,
    body::RigidBody,
    joint::RevoluteJoint,
    model::ArticulatedModel,
    spatial::{SpatialInertia, SpatialTransform},
};

/// Construct a 3-link planar revolute pendulum, run one ABA step, check results.
///
/// Chain geometry (all links horizontal in Y direction):
///   - link 0 (root): attached to world; joint axis = Z
///   - link 1: offset [0, -L, 0] from link 0; joint axis = Z
///   - link 2: offset [0, -L, 0] from link 1; joint axis = Z
///
/// Initial configuration: q = [π/2, 0, 0] — first link horizontal, others straight.
/// This gives a non-zero gravity moment on link 0, so q_ddot[0] != 0.
#[test]
fn articulated_smoke_3link_pendulum() {
    let g = 9.81_f64;
    let l = 1.0_f64; // link length [m]
    let m = 1.0_f64; // link mass [kg]

    // Build model with gravity pointing -Y
    let mut model = ArticulatedModel::new([0.0, -g, 0.0]);

    // Inertia for each link: point mass at COM = [0, -l, 0] in link frame,
    // zero rotational inertia about COM (thin rod approximation — dominated by
    // the parallel-axis term m*l^2 anyway, and from_com handles that).
    let make_inertia = || SpatialInertia::from_com(m, [0.0, -l, 0.0], [[0.0; 3]; 3]);

    // link 0 — root; no parent
    let body0 = RigidBody::new("link0", make_inertia(), None, SpatialTransform::IDENTITY);
    let joint0 = Box::new(RevoluteJoint::new([0.0, 0.0, 1.0]));
    model.add_body(body0, joint0);

    // link 1 — child of link 0; fixed offset [0, -l, 0] from link 0 origin
    let x_t1 = SpatialTransform::from_translation([0.0, -l, 0.0]);
    let body1 = RigidBody::new("link1", make_inertia(), Some(0), x_t1);
    let joint1 = Box::new(RevoluteJoint::new([0.0, 0.0, 1.0]));
    model.add_body(body1, joint1);

    // link 2 — child of link 1; same offset
    let x_t2 = SpatialTransform::from_translation([0.0, -l, 0.0]);
    let body2 = RigidBody::new("link2", make_inertia(), Some(1), x_t2);
    let joint2 = Box::new(RevoluteJoint::new([0.0, 0.0, 1.0]));
    model.add_body(body2, joint2);

    assert_eq!(
        model.total_dof(),
        3,
        "3-link revolute chain must have 3 DOF"
    );

    // Initial state: first link at 90° (horizontal), others straight (0°)
    let q = [std::f64::consts::FRAC_PI_2, 0.0, 0.0];
    let q_dot = [0.0_f64; 3]; // at rest
    let tau = [0.0_f64; 3]; // no applied torques

    // Run forward dynamics (ABA)
    let q_ddot = aba(&model, &q, &q_dot, &tau);

    // ── Assertions ────────────────────────────────────────────────────────────

    assert_eq!(
        q_ddot.len(),
        3,
        "ABA must return exactly 3 accelerations for a 3-DOF model"
    );

    // All accelerations must be finite (no NaN / Inf)
    for (i, &acc) in q_ddot.iter().enumerate() {
        assert!(
            acc.is_finite(),
            "joint acceleration q_ddot[{i}] must be finite, got {acc}"
        );
    }

    // At q = [π/2, 0, 0] the chain is not at equilibrium: gravity creates a
    // non-zero torque on link 0, so at least one acceleration must be non-zero.
    let any_nonzero = q_ddot.iter().any(|&a| a.abs() > 1e-10);
    assert!(
        any_nonzero,
        "at least one joint acceleration must be non-zero under gravity; got {q_ddot:?}"
    );

    // Physical sanity: q_ddot[0] for a horizontal revolute link under gravity
    // should be negative (gravity pulls the link down, causing clockwise rotation).
    // For a single-link pendulum q_ddot = -g/l ≈ -9.81 when horizontal.
    // With 3 links the coupling changes the magnitude but the sign should hold.
    assert!(
        q_ddot[0] < 0.0,
        "root link acceleration should be negative (gravity pulls the chain down), \
         got q_ddot[0] = {}",
        q_ddot[0]
    );
}
