// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for `NarrowPhaseDispatcher::generate_contacts`.
//!
//! Every manifold must carry normals pointing FROM body B TOWARD body A (see
//! `Contact::normal`), witness points on the respective surfaces, and box-box
//! pairs must come with the full face-clip manifold instead of one averaged
//! point.

use oxiphysics_collision::narrowphase::dispatch::{box_box_manifold, box_box_manifold_clip};
use oxiphysics_collision::{CollisionPair, Contact, NarrowPhaseDispatcher};
use oxiphysics_core::Transform;
use oxiphysics_core::math::Vec3;
use oxiphysics_geometry::{BoxShape, Capsule, Shape, Sphere};

const PAIR: CollisionPair = CollisionPair { a: 0, b: 1 };

fn contacts(a: &dyn Shape, ta: &Transform, b: &dyn Shape, tb: &Transform) -> Vec<Contact> {
    NarrowPhaseDispatcher::generate_contacts(a, ta, b, tb, PAIR)
        .map(|m| m.contacts)
        .unwrap_or_default()
}

fn at(x: f64, y: f64, z: f64) -> Transform {
    Transform::from_position(Vec3::new(x, y, z))
}

#[test]
fn box_on_box_yields_a_four_point_manifold_with_b_to_a_normals() {
    let upper = BoxShape::new(Vec3::new(0.5, 0.5, 0.5));
    let lower = BoxShape::new(Vec3::new(1.0, 0.5, 1.0));
    let found = contacts(&upper, &at(0.0, 0.99, 0.0), &lower, &at(0.0, 0.0, 0.0));
    assert_eq!(found.len(), 4, "a face-on-face patch has four corners");
    for c in &found {
        assert!(
            (c.normal - Vec3::new(0.0, 1.0, 0.0)).norm() < 1e-12,
            "{:?}",
            c.normal
        );
        assert!((c.depth - 0.01).abs() < 1e-9, "depth {}", c.depth);
        assert!(
            (c.point_a.y - 0.49).abs() < 1e-9,
            "A's bottom face: {:?}",
            c.point_a
        );
        assert!(
            (c.point_b.y - 0.5).abs() < 1e-9,
            "B's top face: {:?}",
            c.point_b
        );
        assert!((c.point_a.x.abs() - 0.5).abs() < 1e-9 && (c.point_a.z.abs() - 0.5).abs() < 1e-9);
    }
}

#[test]
fn box_box_normal_flips_with_the_argument_order() {
    let cube = BoxShape::new(Vec3::new(1.0, 1.0, 1.0));
    let origin = at(0.0, 0.0, 0.0);
    let right = at(1.5, 0.0, 0.0);
    // B on the +x side of A: the normal from B towards A is -x.
    let ab = contacts(&cube, &origin, &cube, &right);
    assert!(!ab.is_empty());
    for c in &ab {
        assert!(c.normal.x < -0.99, "normal {:?}", c.normal);
        assert!((c.depth - 0.5).abs() < 1e-9);
        assert!(
            (c.point_a.x - 1.0).abs() < 1e-9,
            "A's +x face: {:?}",
            c.point_a
        );
        assert!(
            (c.point_b.x - 0.5).abs() < 1e-9,
            "B's -x face: {:?}",
            c.point_b
        );
    }
    let ba = contacts(&cube, &right, &cube, &origin);
    assert_eq!(ab.len(), ba.len());
    for c in &ba {
        assert!(c.normal.x > 0.99, "normal {:?}", c.normal);
    }
}

#[test]
fn box_box_manifold_clip_is_the_full_manifold() {
    let upper = BoxShape::new(Vec3::new(0.5, 0.5, 0.5));
    let lower = BoxShape::new(Vec3::new(1.0, 0.5, 1.0));
    let ta = at(0.0, 0.99, 0.0);
    let tb = at(0.0, 0.0, 0.0);
    let clip = box_box_manifold_clip(&upper, &ta, &lower, &tb, PAIR).expect("overlap");
    let full = box_box_manifold(&upper, &ta, &lower, &tb, PAIR)
        .manifold
        .expect("overlap");
    assert_eq!(clip.contacts.len(), 4);
    assert_eq!(clip.contacts.len(), full.contacts.len());
    for (c, f) in clip.contacts.iter().zip(&full.contacts) {
        assert!((c.point_a - f.point_a).norm() < 1e-12);
        assert!((c.point_b - f.point_b).norm() < 1e-12);
        assert!((c.normal - f.normal).norm() < 1e-12);
        assert!((c.depth - f.depth).abs() < 1e-12);
    }
    // Separated boxes produce no manifold.
    assert!(box_box_manifold_clip(&upper, &at(0.0, 1.2, 0.0), &lower, &tb, PAIR).is_none());
}

#[test]
fn sphere_sphere_and_capsule_capsule_normals_point_from_b_to_a() {
    let sphere = Sphere::new(0.5);
    let origin = at(0.0, 0.0, 0.0);
    let right = at(0.9, 0.0, 0.0);
    let found = contacts(&sphere, &right, &sphere, &origin);
    assert_eq!(found.len(), 1);
    assert!((found[0].normal - Vec3::new(1.0, 0.0, 0.0)).norm() < 1e-12);
    assert!((found[0].depth - 0.1).abs() < 1e-12);
    assert!((found[0].point_a.x - 0.4).abs() < 1e-12, "A's facing point");
    assert!((found[0].point_b.x - 0.5).abs() < 1e-12, "B's facing point");
    let swapped = contacts(&sphere, &origin, &sphere, &right);
    assert!((swapped[0].normal - Vec3::new(-1.0, 0.0, 0.0)).norm() < 1e-12);

    let capsule = Capsule::new(0.5, 1.0);
    let right = at(0.8, 0.0, 0.0);
    let found = contacts(&capsule, &right, &capsule, &origin);
    assert_eq!(found.len(), 1);
    assert!((found[0].normal - Vec3::new(1.0, 0.0, 0.0)).norm() < 1e-9);
    assert!((found[0].depth - 0.2).abs() < 1e-9);
    let swapped = contacts(&capsule, &origin, &capsule, &right);
    assert!((swapped[0].normal - Vec3::new(-1.0, 0.0, 0.0)).norm() < 1e-9);
}

#[test]
fn sphere_box_normals_point_from_b_to_a_in_both_orders() {
    let sphere = Sphere::new(1.0);
    let cube = BoxShape::new(Vec3::new(1.0, 1.0, 1.0));
    let above = at(0.0, 1.5, 0.0);
    let origin = at(0.0, 0.0, 0.0);
    // Sphere (A) above the box (B): from the box towards the sphere is +y.
    let found = contacts(&sphere, &above, &cube, &origin);
    assert_eq!(found.len(), 1);
    assert!((found[0].normal - Vec3::new(0.0, 1.0, 0.0)).norm() < 1e-12);
    assert!((found[0].depth - 0.5).abs() < 1e-12);
    assert!((found[0].point_a.y - 0.5).abs() < 1e-12, "sphere surface");
    assert!((found[0].point_b.y - 1.0).abs() < 1e-12, "box face");
    // Box (A) below the sphere (B): from the sphere towards the box is -y.
    let swapped = contacts(&cube, &origin, &sphere, &above);
    assert_eq!(swapped.len(), 1);
    assert!((swapped[0].normal - Vec3::new(0.0, -1.0, 0.0)).norm() < 1e-12);
    assert!((swapped[0].point_a.y - 1.0).abs() < 1e-12, "box face");
    assert!((swapped[0].point_b.y - 0.5).abs() < 1e-12, "sphere surface");
}

#[test]
fn sphere_centre_inside_a_box_is_pushed_out_through_the_nearest_face() {
    let sphere = Sphere::new(0.1);
    let cube = BoxShape::new(Vec3::new(1.0, 1.0, 1.0));
    // 0.2 below the +y face.
    let found = contacts(&sphere, &at(0.0, 0.8, 0.0), &cube, &at(0.0, 0.0, 0.0));
    assert_eq!(found.len(), 1);
    assert!((found[0].normal - Vec3::new(0.0, 1.0, 0.0)).norm() < 1e-12);
    assert!(
        (found[0].depth - 0.3).abs() < 1e-12,
        "depth {}",
        found[0].depth
    );
    assert!((found[0].point_b.y - 1.0).abs() < 1e-12);
    assert!((found[0].point_a.y - 0.7).abs() < 1e-12);
}

/// A unit-mass box dropped on a static slab, resolved with the impulse scheme
/// of the physics pipeline (`rel_vel = vel_a - vel_b`, skip when the contact is
/// separating, push A out along the normal): it has to come to rest on top.
#[test]
fn box_resting_on_a_static_box_stays_at_rest_for_120_steps() {
    let dt = 1.0 / 60.0;
    let falling = BoxShape::new(Vec3::new(0.5, 0.5, 0.5));
    let slab = BoxShape::new(Vec3::new(5.0, 0.5, 5.0));
    let slab_t = at(0.0, -0.5, 0.0);
    let mut position = Vec3::new(0.0, 0.52, 0.0);
    let mut velocity = Vec3::zeros();
    for step in 0..120 {
        velocity.y -= 9.81 * dt;
        let ta = Transform::from_position(position);
        if let Some(manifold) = NarrowPhaseDispatcher::generate_contacts(
            &falling as &dyn Shape,
            &ta,
            &slab as &dyn Shape,
            &slab_t,
            PAIR,
        ) {
            let mut deepest: Option<&Contact> = None;
            for c in &manifold.contacts {
                let normal_velocity = velocity.dot(&c.normal);
                if normal_velocity < 0.0 {
                    velocity -= c.normal * normal_velocity;
                }
                if deepest.is_none_or(|d| c.depth > d.depth) {
                    deepest = Some(c);
                }
            }
            // One positional correction per manifold, along the deepest contact.
            if let Some(c) = deepest {
                position += c.normal * c.depth;
            }
        }
        position += velocity * dt;
        assert!(
            position.y > 0.45,
            "step {step}: the box sank through the slab to y = {}",
            position.y
        );
    }
    assert!(
        (position.y - 0.5).abs() < 0.02,
        "the box rests on the slab, y = {}",
        position.y
    );
    assert!(velocity.norm() < 0.2, "velocity {velocity:?}");
}
