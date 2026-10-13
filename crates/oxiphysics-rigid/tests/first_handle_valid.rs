// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//! Regression tests for the null-handle sentinel.
//!
//! The first body and the first collider of every rigid set are issued as
//! `(index 0, generation 0)`. `BodyHandle::null()` and `ColliderHandle::null()`
//! used to be `(0, 0)` as well, so the very first body reported
//! `is_valid() == false` and compared equal to the null handle.

use oxiphysics_core::{BodyHandle, ColliderHandle};
use oxiphysics_geometry::Sphere;
use oxiphysics_rigid::world::PhysicsWorld;
use oxiphysics_rigid::{Collider, ColliderSet, RigidBody, RigidBodySet};
use std::sync::Arc;

#[test]
fn first_body_of_a_rigid_set_is_valid_and_distinct_from_null() {
    let mut set = RigidBodySet::new();
    let first = set.insert(RigidBody::new(1.0));
    assert_eq!((first.index, first.generation), (0, 0));
    assert!(first.is_valid(), "the first body must be a valid handle");
    assert_ne!(first, BodyHandle::null());
    assert!(set.get(first).is_some());

    assert!(!BodyHandle::null().is_valid());
    assert!(set.get(BodyHandle::null()).is_none());

    // A recycled slot gets a new generation and stays valid.
    assert!(set.remove(first).is_some());
    let recycled = set.insert(RigidBody::new(2.0));
    assert_eq!(recycled.index, first.index);
    assert!(recycled.is_valid());
    assert!(set.get(first).is_none(), "the old handle is stale");
}

#[test]
fn first_collider_of_a_collider_set_is_valid_and_distinct_from_null() {
    let mut set = ColliderSet::new();
    let first = set.insert(Collider::new(Arc::new(Sphere::new(0.5))));
    assert_eq!((first.index, first.generation), (0, 0));
    assert!(
        first.is_valid(),
        "the first collider must be a valid handle"
    );
    assert_ne!(first, ColliderHandle::null());
    assert!(set.get(first).is_some());

    assert!(!ColliderHandle::null().is_valid());
    assert!(set.get(ColliderHandle::null()).is_none());
}

#[test]
fn world_hands_out_valid_handles_from_the_very_first_body_and_collider() {
    let mut world = PhysicsWorld::new([0.0, -9.81, 0.0], 1.0 / 60.0);
    let body = world.add_rigid_body(RigidBody::new(1.0));
    let collider = world.add_collider(Collider::new(Arc::new(Sphere::new(0.5))), Some(body));
    assert!(body.is_valid());
    assert!(collider.is_valid());
    assert!(world.get_body(body).is_some());
    assert!(world.get_body(BodyHandle::null()).is_none());
}
