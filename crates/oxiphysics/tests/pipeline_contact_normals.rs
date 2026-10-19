// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Pipeline-level regression test for the narrow-phase contact-normal
//! convention.
//!
//! The contact solver applies an impulse only while the bodies approach along
//! the contact normal, which presumes a normal pointing from body B toward body
//! A. Two dynamic spheres (or capsules) that meet head-on must therefore
//! rebound, whichever of them is stored first in the collider set.

use oxiphysics::core::math::{Quat, Vec3};
use oxiphysics::core::{PhysicsConfig, Transform};
use oxiphysics::geometry::{Capsule, Shape, Sphere};
use oxiphysics::pipeline::PhysicsPipeline;
use oxiphysics::rigid::{Collider, ColliderSet, RigidBody, RigidBodySet};
use std::f64::consts::FRAC_PI_2;
use std::sync::Arc;

/// Zero gravity and bodies that never fall asleep.
fn free_space_config() -> PhysicsConfig {
    PhysicsConfig {
        gravity: Vec3::zeros(),
        time_before_sleep: 1.0e9,
        linear_sleep_threshold: 0.0,
        angular_sleep_threshold: 0.0,
        ..PhysicsConfig::default()
    }
}

/// A body of unit mass at `position` with the given orientation and velocity,
/// without damping.
fn free_body(position: Vec3, rotation: Quat, velocity: Vec3) -> RigidBody {
    let mut body = RigidBody::new(1.0);
    body.transform = Transform::new(position, rotation);
    body.velocity = velocity;
    body.linear_damping = 0.0;
    body.angular_damping = 0.0;
    body
}

/// Velocities after one step of two elastic bodies that overlap by 0.02 and
/// move toward each other at 1 m/s each. The left body moves along +x. When
/// `right_first` is set the right body is inserted first, which swaps the
/// roles of bodies A and B in the narrow phase. The right body is turned by
/// `right_rotation`.
fn collide_head_on(
    make_shape: impl Fn() -> Arc<dyn Shape>,
    right_rotation: Quat,
    spacing: f64,
    right_first: bool,
) -> (Vec3, Vec3) {
    let left = free_body(
        Vec3::new(-0.5 * spacing, 0.0, 0.0),
        Quat::identity(),
        Vec3::new(1.0, 0.0, 0.0),
    );
    let right = free_body(
        Vec3::new(0.5 * spacing, 0.0, 0.0),
        right_rotation,
        Vec3::new(-1.0, 0.0, 0.0),
    );
    let mut bodies = RigidBodySet::new();
    let mut colliders = ColliderSet::new();
    let order = if right_first {
        [right, left]
    } else {
        [left, right]
    };
    let mut handles = Vec::new();
    for body in order {
        let handle = bodies.insert(body);
        colliders.insert(
            Collider::new(make_shape())
                .with_body(handle)
                .with_friction(0.0)
                .with_restitution(1.0),
        );
        handles.push(handle);
    }
    let mut pipeline = PhysicsPipeline::with_config(free_space_config());
    pipeline.step(1.0 / 240.0, &mut bodies, &colliders);
    let velocity = |index: usize| {
        handles
            .get(index)
            .and_then(|&h| bodies.get(h))
            .map(|b| b.velocity)
            .unwrap_or_default()
    };
    // Report (left body, right body) whatever the insertion order.
    if right_first {
        (velocity(1), velocity(0))
    } else {
        (velocity(0), velocity(1))
    }
}

fn assert_rebound(left: Vec3, right: Vec3, what: &str) {
    assert!(
        right.x - left.x > 1.5,
        "{what}: the bodies must separate after the impact, got left {left:?} right {right:?}"
    );
    assert!(left.x < 0.0 && right.x > 0.0, "{what}: velocities reverse");
    assert!(
        (left.x + right.x).abs() < 1e-6,
        "{what}: equal masses conserve momentum, got {} and {}",
        left.x,
        right.x
    );
    assert!(
        left.y.abs() < 1e-9 && left.z.abs() < 1e-9,
        "{what}: {left:?}"
    );
    assert!(
        right.y.abs() < 1e-9 && right.z.abs() < 1e-9,
        "{what}: {right:?}"
    );
}

#[test]
fn head_on_spheres_rebound_in_either_collider_order() {
    for right_first in [false, true] {
        let (left, right) = collide_head_on(
            || -> Arc<dyn Shape> { Arc::new(Sphere::new(0.5)) },
            Quat::identity(),
            0.98,
            right_first,
        );
        assert_rebound(left, right, &format!("spheres, right first: {right_first}"));
    }
}

#[test]
fn crossed_capsules_rebound_in_either_collider_order() {
    // The left capsule stands along y, the right one lies along z (a quarter
    // turn about x). Their axes cross on the line between the centres, so the
    // contact acts through both centres of mass and the exchange is purely
    // linear (parallel capsules would put the contact at an end point of the
    // axes, off the line of centres, and part of the impulse would spin them).
    for right_first in [false, true] {
        let (left, right) = collide_head_on(
            || -> Arc<dyn Shape> { Arc::new(Capsule::new(0.3, 0.5)) },
            Quat::from_axis_angle(&Vec3::x_axis(), FRAC_PI_2),
            0.58,
            right_first,
        );
        assert_rebound(
            left,
            right,
            &format!("capsules, right first: {right_first}"),
        );
    }
}

#[test]
fn separated_spheres_do_not_interact() {
    let (left, right) = collide_head_on(
        || -> Arc<dyn Shape> { Arc::new(Sphere::new(0.5)) },
        Quat::identity(),
        1.05,
        false,
    );
    assert_eq!(left, Vec3::new(1.0, 0.0, 0.0));
    assert_eq!(right, Vec3::new(-1.0, 0.0, 0.0));
}
