// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tests for the world-frame convention of the angular velocity.
//!
//! `Body::integrate_positions` used to rotate a body by the angular velocity
//! as if it were expressed in the body frame, while impulses, contacts and the
//! constraint solver all treat it as a world-frame vector. A body that had
//! already turned then spun about `R * w` instead of `w`, which is harmless
//! for planar motion but makes any constrained 3-D rotation unstable.

use super::{Body, Constraint, ConstraintKind, PhysicsWorld};
use crate::math::{Quat, Vec3};
use crate::types::{BodyHandle, TimeStep};
use std::f64::consts::FRAC_PI_2;

fn undamped_body(position: [f64; 3], mass: f64, inertia: f64) -> Body {
    let mut body = Body::new_dynamic(position, mass);
    body.inertia = inertia;
    body.linear_damping = 0.0;
    body.angular_damping = 0.0;
    body
}

fn ball_socket(
    body_a: BodyHandle,
    body_b: Option<BodyHandle>,
    anchor_a: [f64; 3],
    anchor_b: [f64; 3],
) -> Constraint {
    Constraint {
        body_a,
        body_b,
        kind: ConstraintKind::BallSocket,
        anchor_a,
        anchor_b,
        target: 0.0,
        compliance: 0.0,
        active: true,
    }
}

fn world_anchor(world: &PhysicsWorld, handle: BodyHandle, local: [f64; 3]) -> Vec3 {
    world
        .body(handle)
        .map(|b| {
            b.transform
                .transform_point(&Vec3::new(local[0], local[1], local[2]))
        })
        .unwrap_or_else(Vec3::zeros)
}

#[test]
fn a_rotated_body_spins_about_the_world_axis_of_its_angular_velocity() {
    // Pre-rotated by 90 degrees about Z and spun with w = (1, 0, 0) for pi/2
    // seconds: the extra rotation is about the *world* X axis, so the result
    // is Rx(90) * Rz(90); the body-frame reading gives Rz(90) * Rx(90).
    let dt = 1.0e-3;
    let steps = (FRAC_PI_2 / dt).round() as usize;
    let swept_angle = steps as f64 * dt;
    let mut body = Body::new_dynamic([0.0; 3], 1.0);
    body.transform.rotation = Quat::from_axis_angle(&Vec3::z_axis(), FRAC_PI_2);
    body.velocity.angular = [1.0, 0.0, 0.0];
    for _ in 0..steps {
        body.integrate_positions(dt);
    }
    let world_frame = Quat::from_axis_angle(&Vec3::x_axis(), swept_angle)
        * Quat::from_axis_angle(&Vec3::z_axis(), FRAC_PI_2);
    let body_frame = Quat::from_axis_angle(&Vec3::z_axis(), FRAC_PI_2)
        * Quat::from_axis_angle(&Vec3::x_axis(), swept_angle);
    let rotation = body.transform.rotation;
    assert!(
        rotation.angle_to(&world_frame) < 1e-5,
        "off the world-frame result by {}",
        rotation.angle_to(&world_frame)
    );
    assert!(
        rotation.angle_to(&body_frame) > 1.0,
        "must not match the body-frame reading"
    );
    // The orientation stays a unit quaternion.
    assert!((rotation.quaternion().norm() - 1.0).abs() < 1e-12);
}

#[test]
fn a_body_without_prior_rotation_is_unaffected_by_the_frame_convention() {
    let dt = 1.0e-3;
    let steps = 500;
    let mut body = Body::new_dynamic([0.0; 3], 1.0);
    body.velocity.angular = [0.0, 0.0, 2.0];
    for _ in 0..steps {
        body.integrate_positions(dt);
    }
    let expected = Quat::from_axis_angle(&Vec3::z_axis(), 2.0 * steps as f64 * dt);
    assert!(body.transform.rotation.angle_to(&expected) < 1e-5);
}

#[test]
fn offset_ball_socket_with_a_general_spin_stays_pinned_and_never_gains_energy() {
    // A spherical pendulum: the body's centre is 1 m from the pivot and it is
    // spun about all three axes at once. With the body-frame reading of the
    // angular velocity this blows up (energy grows by nine orders of magnitude
    // and the angular velocity reaches tens of thousands of rad/s).
    let dt = 1.0 / 240.0;
    let mut world = PhysicsWorld::new();
    let mut body = undamped_body([1.0, 0.0, 0.0], 1.0, 0.5);
    body.velocity.angular = [0.7, 1.1, -0.4];
    let handle = world.add_body(body);
    world.add_constraint(ball_socket(handle, None, [-1.0, 0.0, 0.0], [0.0, 0.0, 0.0]));

    let energy = |world: &PhysicsWorld| {
        let height = world.body(handle).map_or(0.0, |b| b.transform.position.y);
        world.total_kinetic_energy() + 9.81 * height
    };
    let initial_energy = energy(&world);
    let mut worst_gap = 0.0_f64;
    let mut highest_energy = initial_energy;
    let mut fastest_spin = 0.0_f64;
    for _ in 0..(6.0 / dt) as usize {
        world.step(&TimeStep::new(dt));
        worst_gap = worst_gap.max(world_anchor(&world, handle, [-1.0, 0.0, 0.0]).norm());
        highest_energy = highest_energy.max(energy(&world));
        let spin = world
            .body(handle)
            .map_or(0.0, |b| Vec3::from(b.velocity.angular).norm());
        fastest_spin = fastest_spin.max(spin);
    }
    assert!(worst_gap < 1e-4, "pivot drifted by {worst_gap}");
    assert!(
        highest_energy <= initial_energy + 1e-6,
        "energy rose from {initial_energy} to {highest_energy}"
    );
    assert!(fastest_spin < 10.0, "spin reached {fastest_spin} rad/s");
}

#[test]
fn spinning_two_body_ball_socket_conserves_momentum_in_three_dimensions() {
    let dt = 1.0 / 120.0;
    let mut world = PhysicsWorld::new();
    world.config.gravity = Vec3::zeros();
    let mut light = undamped_body([0.0, 0.0, 0.0], 1.0, 0.5);
    light.velocity.linear = [0.0, 1.0, 0.2];
    light.velocity.angular = [0.5, -0.3, 0.8];
    let mut heavy = undamped_body([2.0, 0.0, 0.0], 3.0, 2.0);
    heavy.velocity.angular = [-0.2, 0.4, 0.1];
    let a = world.add_body(light);
    let b = world.add_body(heavy);
    world.add_constraint(ball_socket(a, Some(b), [1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]));

    let momentum = world.total_momentum();
    let angular_momentum = world.total_angular_momentum();
    let angular_scale = Vec3::from(angular_momentum).norm();
    let initial_energy = world.total_kinetic_energy();
    let mut worst_gap = 0.0_f64;
    let mut worst_angular_drift = 0.0_f64;
    for _ in 0..(5.0 / dt) as usize {
        world.step(&TimeStep::new(dt));
        let gap = (world_anchor(&world, a, [1.0, 0.0, 0.0])
            - world_anchor(&world, b, [-1.0, 0.0, 0.0]))
        .norm();
        worst_gap = worst_gap.max(gap);
        let drift =
            (Vec3::from(world.total_angular_momentum()) - Vec3::from(angular_momentum)).norm();
        worst_angular_drift = worst_angular_drift.max(drift);
        assert!(world.total_kinetic_energy() <= initial_energy + 1e-9);
    }
    assert!(worst_gap < 1e-5, "anchors drifted by {worst_gap}");
    let final_momentum = world.total_momentum();
    for k in 0..3 {
        assert!(
            (final_momentum[k] - momentum[k]).abs() < 1e-12,
            "linear momentum component {k} changed"
        );
    }
    assert!(
        worst_angular_drift < 0.02 * angular_scale,
        "angular momentum drifted by {worst_angular_drift} of {angular_scale}"
    );
}
