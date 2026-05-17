// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Integration tests for articulated-body dynamics (RNEA + ABA).
//!
//! Tests validate:
//! 1. Single-pendulum RNEA against analytical gravitational torque.
//! 2. Single-pendulum ABA against analytical gravitational acceleration.
//! 3. RNEA/ABA inverse consistency on a 2-DOF robot arm.

use oxiphysics_articulated::{
    aba::aba,
    body::RigidBody,
    joint::RevoluteJoint,
    model::ArticulatedModel,
    rnea::rnea,
    spatial::{SpatialInertia, SpatialTransform},
};
use std::f64::consts::FRAC_PI_2;

fn approx_eq(a: f64, b: f64, eps: f64, label: &str) {
    assert!(
        (a - b).abs() < eps,
        "{label}: expected {b:.8}, got {a:.8}, diff={:.2e}",
        (a - b).abs()
    );
}

/// Build a single-link revolute pendulum.
///
/// - Revolute joint about z-axis at origin.
/// - Point mass `m` at `[0, -l, 0]` in body frame (hanging at q=0).
/// - Gravity [0, -g, 0].
fn build_pendulum(m: f64, l: f64, g: f64) -> ArticulatedModel {
    let mut model = ArticulatedModel::new([0.0, -g, 0.0]);
    let joint = Box::new(RevoluteJoint::new([0.0, 0.0, 1.0]));
    // COM is at [0, -l, 0] in body frame; rotational inertia of point mass = 0 about COM
    let inertia = SpatialInertia::from_com(m, [0.0, -l, 0.0], [[0.0; 3]; 3]);
    let body = RigidBody::new("link1", inertia, None, SpatialTransform::IDENTITY);
    model.add_body(body, joint);
    model
}

/// Build a 2-DOF planar robot arm:
/// - Link 1: revolute about z at origin, COM at [l1, 0, 0] in link-1 frame.
/// - Link 2: revolute about z at tip of link 1 ([l1, 0, 0]), COM at [l2, 0, 0] in link-2 frame.
fn build_two_link_arm(m1: f64, l1: f64, m2: f64, l2: f64, g: f64) -> ArticulatedModel {
    let mut model = ArticulatedModel::new([0.0, -g, 0.0]);

    // Link 1: revolute about z at origin
    let j1 = Box::new(RevoluteJoint::new([0.0, 0.0, 1.0]));
    // COM at [l1, 0, 0] in link-1 frame (extends in +x direction)
    let inertia1 = SpatialInertia::from_com(m1, [l1, 0.0, 0.0], [[0.0; 3]; 3]);
    let body1 = RigidBody::new("link1", inertia1, None, SpatialTransform::IDENTITY);
    let idx1 = model.add_body(body1, j1);

    // Link 2: revolute about z at [2*l1, 0, 0] in link-1 frame
    // The fixed offset (X_T) puts the joint at the tip of link 1
    let x_t2 = SpatialTransform::from_translation([2.0 * l1, 0.0, 0.0]);
    let j2 = Box::new(RevoluteJoint::new([0.0, 0.0, 1.0]));
    // COM at [l2, 0, 0] in link-2 frame (relative to joint 2 origin)
    let inertia2 = SpatialInertia::from_com(m2, [l2, 0.0, 0.0], [[0.0; 3]; 3]);
    let body2 = RigidBody::new("link2", inertia2, Some(idx1), x_t2);
    model.add_body(body2, j2);

    model
}

// ─── Test 1: Single pendulum RNEA ────────────────────────────────────────────

#[test]
fn test_single_pendulum_rnea_hanging() {
    // At q=0 (mass hanging directly below pivot), gravitational torque = 0.
    let (m, l, g) = (1.0, 1.0, 9.81);
    let model = build_pendulum(m, l, g);

    let tau = rnea(&model, &[0.0], &[0.0], &[0.0]);
    approx_eq(tau[0], 0.0, 1e-8, "RNEA at q=0: tau should be 0");
}

#[test]
fn test_single_pendulum_rnea_horizontal() {
    // At q=π/2 (mass at [l, 0, 0] after rotation), with zero velocity and acceleration:
    // Gravity [0,-g,0] acts on mass at [l,0,0] → torque = m·g·l about z-axis.
    // RNEA returns the torque needed to keep that motion state, so τ = m·g·l.
    let (m, l, g) = (1.0, 1.0, 9.81);
    let model = build_pendulum(m, l, g);

    let tau = rnea(&model, &[FRAC_PI_2], &[0.0], &[0.0]);
    let expected = m * g * l;
    approx_eq(tau[0], expected, 1e-6, "RNEA at q=π/2: tau should be m*g*l");
}

// ─── Test 2: Single pendulum ABA ─────────────────────────────────────────────

#[test]
fn test_single_pendulum_aba_gravity() {
    // At q=π/2 (horizontal), zero velocity, no applied torque:
    // Expected q̈ = -g/l (for a point mass pendulum: θ̈ = −(g/l)·sin(π/2) = −g/l).
    let (m, l, g) = (1.0, 1.0, 9.81);
    let model = build_pendulum(m, l, g);

    let q_ddot = aba(&model, &[FRAC_PI_2], &[0.0], &[0.0]);
    let expected = -g / l;
    approx_eq(
        q_ddot[0],
        expected,
        1e-6,
        "ABA pendulum at horizontal: q̈ should be -g/L",
    );
}

#[test]
fn test_single_pendulum_aba_hanging() {
    // At q=0 (hanging), no torque → q̈ = 0 (equilibrium).
    let (m, l, g) = (1.0, 1.0, 9.81);
    let model = build_pendulum(m, l, g);

    let q_ddot = aba(&model, &[0.0], &[0.0], &[0.0]);
    approx_eq(q_ddot[0], 0.0, 1e-8, "ABA at q=0: q̈ should be 0");
}

// ─── Test 3: RNEA/ABA inverse consistency ────────────────────────────────────

#[test]
fn test_rnea_aba_inverse_consistency_pendulum() {
    // D'Alembert principle: RNEA(q, q̇, q̈) = τ, then ABA(q, q̇, τ) = q̈.
    // Must hold for any motion state.
    let (m, l, g) = (2.0, 0.8, 9.81);
    let model = build_pendulum(m, l, g);

    let q = [0.4];
    let q_dot = [0.3];
    let q_ddot_in = [1.2];

    let tau = rnea(&model, &q, &q_dot, &q_ddot_in);
    let q_ddot_out = aba(&model, &q, &q_dot, &tau);

    approx_eq(
        q_ddot_out[0],
        q_ddot_in[0],
        1e-8,
        "RNEA/ABA inverse consistency (pendulum)",
    );
}

#[test]
fn test_rnea_aba_inverse_consistency_two_link() {
    // Two-link arm: RNEA then ABA must be exact inverses.
    let (m1, l1, m2, l2, g) = (2.0, 0.5, 1.0, 0.4, 9.81);
    let model = build_two_link_arm(m1, l1, m2, l2, g);

    let q = [0.3, -0.7];
    let q_dot = [0.5, -0.3];
    let q_ddot_in = [0.1, -0.2];

    let tau = rnea(&model, &q, &q_dot, &q_ddot_in);
    let q_ddot_out = aba(&model, &q, &q_dot, &tau);

    for (k, (a_in, a_out)) in q_ddot_in.iter().zip(q_ddot_out.iter()).enumerate() {
        approx_eq(
            *a_out,
            *a_in,
            1e-8,
            &format!("RNEA/ABA inverse consistency (2-link), joint {k}"),
        );
    }
}
