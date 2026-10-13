// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the Plücker-transform convention (`rot` is the
//! child→parent pose) and the gravity sign (`a_0 = −g`), checked against
//! closed-form Lagrangian dynamics.

use oxiphysics_articulated::{
    aba::aba,
    body::RigidBody,
    crba::compute_mass_matrix_crba,
    joint::{PrismaticJoint, RevoluteJoint},
    model::ArticulatedModel,
    rnea::rnea,
    spatial::{SpatialInertia, SpatialTransform},
};

fn assert_close(a: f64, b: f64, eps: f64, label: &str) {
    assert!(
        (a - b).abs() < eps,
        "{label}: expected {b:.12}, got {a:.12}, diff={:.3e}",
        (a - b).abs()
    );
}

/// Planar 2-link arm in the x–y plane (both joints about +z), gravity −y.
///
/// Link i has mass `m_i`, COM at distance `lc_i` along its local x axis, and
/// rotational inertia `i_i` about the COM's z axis. Link 2's joint sits at
/// `[l1, 0, 0]` in link 1's frame.
struct PlanarArm {
    m1: f64,
    m2: f64,
    l1: f64,
    lc1: f64,
    lc2: f64,
    i1: f64,
    i2: f64,
    g: f64,
}

impl PlanarArm {
    fn model(&self) -> ArticulatedModel {
        let mut model = ArticulatedModel::new([0.0, -self.g, 0.0]);
        let inertia1 = SpatialInertia::from_com(
            self.m1,
            [self.lc1, 0.0, 0.0],
            [[0.3, 0.0, 0.0], [0.0, 0.2, 0.0], [0.0, 0.0, self.i1]],
        );
        model.add_body(
            RigidBody::new("link1", inertia1, None, SpatialTransform::IDENTITY),
            Box::new(RevoluteJoint::new([0.0, 0.0, 1.0])),
        );
        let inertia2 = SpatialInertia::from_com(
            self.m2,
            [self.lc2, 0.0, 0.0],
            [[0.1, 0.0, 0.0], [0.0, 0.15, 0.0], [0.0, 0.0, self.i2]],
        );
        model.add_body(
            RigidBody::new(
                "link2",
                inertia2,
                Some(0),
                SpatialTransform::from_translation([self.l1, 0.0, 0.0]),
            ),
            Box::new(RevoluteJoint::new([0.0, 0.0, 1.0])),
        );
        model
    }

    /// Closed-form inverse dynamics τ = M(q)·q̈ + C(q, q̇)·q̇ + G(q).
    fn lagrangian_tau(&self, q: [f64; 2], qd: [f64; 2], qdd: [f64; 2]) -> [f64; 2] {
        let (c2, s2) = (q[1].cos(), q[1].sin());
        let m11 = self.m1 * self.lc1 * self.lc1
            + self.i1
            + self.m2 * (self.l1 * self.l1 + self.lc2 * self.lc2 + 2.0 * self.l1 * self.lc2 * c2)
            + self.i2;
        let m12 = self.m2 * (self.lc2 * self.lc2 + self.l1 * self.lc2 * c2) + self.i2;
        let m22 = self.m2 * self.lc2 * self.lc2 + self.i2;
        let h = self.m2 * self.l1 * self.lc2 * s2;
        let g1 = (self.m1 * self.lc1 + self.m2 * self.l1) * self.g * q[0].cos()
            + self.m2 * self.lc2 * self.g * (q[0] + q[1]).cos();
        let g2 = self.m2 * self.lc2 * self.g * (q[0] + q[1]).cos();
        [
            m11 * qdd[0] + m12 * qdd[1] - h * (2.0 * qd[0] * qd[1] + qd[1] * qd[1]) + g1,
            m12 * qdd[0] + m22 * qdd[1] + h * qd[0] * qd[0] + g2,
        ]
    }
}

fn arm() -> PlanarArm {
    PlanarArm {
        m1: 1.7,
        m2: 0.9,
        l1: 0.8,
        lc1: 0.35,
        lc2: 0.45,
        i1: 0.04,
        i2: 0.02,
        g: 9.81,
    }
}

#[test]
fn rnea_two_link_matches_closed_form_lagrangian() {
    let arm = arm();
    let model = arm.model();
    let cases = [
        ([0.3, -0.7], [1.2, -0.8], [0.5, -1.1]),
        ([-1.1, 0.9], [-0.6, 2.0], [0.0, 0.0]),
        ([0.0, 0.0], [0.0, 0.0], [0.0, 0.0]),
        ([2.4, -2.0], [3.0, 1.5], [-2.0, 0.7]),
    ];
    for (q, qd, qdd) in cases {
        let tau = rnea(&model, &q, &qd, &qdd);
        let expected = arm.lagrangian_tau(q, qd, qdd);
        assert_close(tau[0], expected[0], 1e-9, "tau[0]");
        assert_close(tau[1], expected[1], 1e-9, "tau[1]");
    }
}

#[test]
fn aba_two_link_inverts_closed_form_lagrangian() {
    let arm = arm();
    let model = arm.model();
    let q = [0.4, 1.1];
    let qd = [-0.9, 0.6];
    let qdd = [1.3, -0.4];
    let tau = arm.lagrangian_tau(q, qd, qdd);
    let qdd_aba = aba(&model, &q, &qd, &tau);
    assert_close(qdd_aba[0], qdd[0], 1e-9, "qdd[0]");
    assert_close(qdd_aba[1], qdd[1], 1e-9, "qdd[1]");
}

#[test]
fn crba_two_link_matches_closed_form_mass_matrix() {
    let arm = arm();
    let model = arm.model();
    let q = [0.2, -1.3];
    let m = compute_mass_matrix_crba(&model, &q);
    // Columns of M from the closed form with q̇ = 0 and gravity removed.
    let g_only = arm.lagrangian_tau(q, [0.0; 2], [0.0; 2]);
    for k in 0..2 {
        let mut e = [0.0; 2];
        e[k] = 1.0;
        let col = arm.lagrangian_tau(q, [0.0; 2], e);
        for i in 0..2 {
            assert_close(m[i][k], col[i] - g_only[i], 1e-10, "M entry");
        }
    }
}

#[test]
fn prismatic_slider_falls_with_minus_g() {
    // A free slider along world y under gravity [0, -g, 0] must accelerate
    // downward: q̈ = −g.
    let g = 9.81;
    let mut model = ArticulatedModel::new([0.0, -g, 0.0]);
    let inertia = SpatialInertia::from_com(
        2.5,
        [0.1, 0.2, -0.3],
        [[0.1, 0.0, 0.0], [0.0, 0.1, 0.0], [0.0, 0.0, 0.1]],
    );
    model.add_body(
        RigidBody::new("slider", inertia, None, SpatialTransform::IDENTITY),
        Box::new(PrismaticJoint::new([0.0, 1.0, 0.0])),
    );
    let qdd = aba(&model, &[0.7], &[0.0], &[0.0]);
    assert_close(qdd[0], -g, 1e-12, "slider q̈");
    // Holding force is +m·g (pushing up).
    let tau = rnea(&model, &[0.7], &[0.0], &[0.0]);
    assert_close(tau[0], 2.5 * g, 1e-12, "slider holding force");
}

#[test]
fn horizontal_link_holding_torque_is_plus_mgl() {
    // Revolute about z, COM at [l, 0, 0], gravity −y: at q = 0 gravity pulls
    // the COM down (negative torque), so the holding torque is +m·g·l.
    let (m, l, g) = (1.3, 0.75, 9.81);
    let mut model = ArticulatedModel::new([0.0, -g, 0.0]);
    let inertia = SpatialInertia::from_com(m, [l, 0.0, 0.0], [[0.0; 3]; 3]);
    model.add_body(
        RigidBody::new("link", inertia, None, SpatialTransform::IDENTITY),
        Box::new(RevoluteJoint::new([0.0, 0.0, 1.0])),
    );
    let tau = rnea(&model, &[0.0], &[0.0], &[0.0]);
    assert_close(tau[0], m * g * l, 1e-12, "holding torque");
    let qdd = aba(&model, &[0.0], &[0.0], &[0.0]);
    assert_close(qdd[0], -g / l, 1e-12, "free fall angular acceleration");
}

#[test]
fn rotated_fixed_offset_uses_pose_convention() {
    // The fixed transform is the pose of the joint frame: rotated +90° about
    // z, so the child's local +x axis points along world +y. A COM at local
    // [l, 0, 0] therefore sits straight above the pivot (zero gravity torque
    // about z), while a COM at local [0, −l, 0] sits at world [+l, 0, 0] and
    // needs a holding torque of +m·g·l.
    let (m, l, g) = (2.0, 0.5, 9.81);
    let x_t = SpatialTransform::from_axis_angle([0.0, 0.0, 1.0], std::f64::consts::FRAC_PI_2);

    let mut upright = ArticulatedModel::new([0.0, -g, 0.0]);
    upright.add_body(
        RigidBody::new(
            "upright",
            SpatialInertia::from_com(m, [l, 0.0, 0.0], [[0.0; 3]; 3]),
            None,
            x_t,
        ),
        Box::new(RevoluteJoint::new([0.0, 0.0, 1.0])),
    );
    let tau_up = rnea(&upright, &[0.0], &[0.0], &[0.0]);
    assert_close(tau_up[0], 0.0, 1e-12, "upright torque");

    let mut sideways = ArticulatedModel::new([0.0, -g, 0.0]);
    sideways.add_body(
        RigidBody::new(
            "sideways",
            SpatialInertia::from_com(m, [0.0, -l, 0.0], [[0.0; 3]; 3]),
            None,
            x_t,
        ),
        Box::new(RevoluteJoint::new([0.0, 0.0, 1.0])),
    );
    let tau_side = rnea(&sideways, &[0.0], &[0.0], &[0.0]);
    assert_close(tau_side[0], m * g * l, 1e-12, "sideways torque");
}
