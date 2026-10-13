// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tests of the pipeline's collision detection: rotation-aware world bounds,
//! shape classification, sensors and pair filtering.

use super::*;
use crate::pipeline::ccd::predicted_transform;
use crate::pipeline::collision::{classify_shape, collect_entries, world_aabb};
use oxiphysics_collision::narrowphase::dispatch::ShapeType;
use oxiphysics_core::math::Quat;
use oxiphysics_geometry::{Capsule, Cone, ConvexHull, Cylinder};
use std::f64::consts::{FRAC_PI_2, FRAC_PI_4};

/// Two dynamic boxes (half-extent 1) rotated 45 degrees about Y, 2.2 apart on X.
fn rotated_box_pair(rotation: Quat) -> Scene {
    let mut scene = Scene::new();
    for x in [0.0, 2.2] {
        let mut body = RigidBody::new(1.0);
        body.transform = Transform::new(Vec3::new(x, 0.0, 0.0), rotation);
        scene.add(body, cuboid(1.0, 1.0, 1.0), 0.5, 0.0);
    }
    scene
}

#[test]
fn rotated_boxes_are_found_by_the_broadphase_and_make_a_manifold() {
    // Seen along the boxes' own axes the centres are 2.2 * cos(45 deg) = 1.56
    // apart on both horizontal axes, so the boxes interpenetrate by 0.44. Their
    // *local* AABBs (x in [-1, 1] and [1.2, 3.2]) do not touch, so a broadphase
    // that ignores rotation misses the pair.
    let scene = rotated_box_pair(Quat::from_axis_angle(&Vec3::y_axis(), FRAC_PI_4));
    let pipeline = PhysicsPipeline::new();
    let entries = collect_entries(&scene.colliders);
    let detection = pipeline.detect_collisions(&entries, &scene.bodies);

    assert_eq!(detection.broadphase_pairs, 1, "world AABBs overlap");
    assert_eq!(detection.contacts.len(), 1, "the boxes interpenetrate");
    let depth = detection.contacts[0].manifold.max_depth();
    assert!(
        (0.3..0.6).contains(&depth),
        "expected a penetration of about 0.44, got {depth}"
    );
}

#[test]
fn unrotated_boxes_at_the_same_distance_do_not_touch() {
    let scene = rotated_box_pair(Quat::identity());
    let pipeline = PhysicsPipeline::new();
    let entries = collect_entries(&scene.colliders);
    let detection = pipeline.detect_collisions(&entries, &scene.bodies);

    // Centres 2.2 apart, half-extents 1 + 1: a 0.2 gap.
    assert_eq!(detection.broadphase_pairs, 0);
    assert!(detection.contacts.is_empty());
}

#[test]
fn upright_plank_collides_with_a_sphere_at_its_top_end() {
    // A 4 m long plank lying along X, turned upright by a quarter turn about Z.
    let mut scene = Scene::new();
    let mut plank = RigidBody::new_static();
    plank.transform = Transform::new(
        Vec3::zeros(),
        Quat::from_axis_angle(&Vec3::z_axis(), FRAC_PI_2),
    );
    scene.add(plank, cuboid(2.0, 0.1, 0.1), 0.5, 0.0);
    // The sphere dips 0.1 m into the plank's top end (y = 2).
    scene.dynamic(1.0, Vec3::new(0.0, 2.2, 0.0), sphere(0.3), 0.5, 0.0);

    let pipeline = PhysicsPipeline::new();
    let entries = collect_entries(&scene.colliders);
    let detection = pipeline.detect_collisions(&entries, &scene.bodies);

    assert_eq!(detection.contacts.len(), 1);
    let contact = &detection.contacts[0].manifold.contacts[0];
    assert!(contact.normal.y.abs() > 0.99, "normal {:?}", contact.normal);
    assert!(
        (contact.depth - 0.1).abs() < 0.02,
        "depth {} should be about 0.1",
        contact.depth
    );
}

#[test]
fn world_aabb_uses_the_collider_offset_and_body_rotation() {
    let mut scene = Scene::new();
    let mut body = RigidBody::new(1.0);
    body.transform = Transform::new(
        Vec3::new(1.0, 0.0, 0.0),
        Quat::from_axis_angle(&Vec3::z_axis(), FRAC_PI_2),
    );
    let body_handle = scene.bodies.insert(body);
    // The collider sits 2 m along the body's local X axis.
    scene.colliders.insert(
        Collider::new(cuboid(0.5, 0.25, 0.25))
            .with_body(body_handle)
            .with_local_transform(Transform::from_position(Vec3::new(2.0, 0.0, 0.0))),
    );
    let entries = collect_entries(&scene.colliders);
    let world = entries[0].world_transform(&scene.bodies);
    let aabb = world_aabb(entries[0].collider, &world);

    // The quarter turn about Z carries local +X to world +Y: the collider
    // centre is at (1, 2, 0) and the box is 0.5 high-extent along Y.
    assert!((aabb.center() - Vec3::new(1.0, 2.0, 0.0)).norm() < 1e-12);
    assert!((aabb.half_extents() - Vec3::new(0.25, 0.5, 0.25)).norm() < 1e-12);
}

#[test]
fn predicted_transform_advances_position_and_orientation() {
    let mut body = RigidBody::new(1.0);
    body.transform = Transform::from_position(Vec3::new(1.0, 2.0, 3.0));
    body.velocity = Vec3::new(2.0, 0.0, -4.0);
    body.angular_velocity = Vec3::new(0.0, 0.0, std::f64::consts::PI);

    // Half a second: translate by v/2 and turn a quarter circle about Z.
    let end = predicted_transform(&body, 0.5);
    assert!((end.position - Vec3::new(2.0, 2.0, 1.0)).norm() < 1e-12);
    let expected = Quat::from_axis_angle(&Vec3::z_axis(), FRAC_PI_2);
    assert!(
        end.rotation.angle_to(&expected) < 1e-12,
        "orientation was not advanced: {:?}",
        end.rotation
    );

    // The prediction matches what the integrator does in one step.
    let mut stepped = body.clone();
    stepped.integrate_velocity(0.5);
    assert!((stepped.transform.position - end.position).norm() < 1e-12);
    assert!(stepped.transform.rotation.angle_to(&end.rotation) < 1e-12);
}

#[test]
fn predicted_transform_keeps_static_sleeping_and_scripted_bodies_in_place() {
    let mut fixed = RigidBody::new_static();
    fixed.velocity = Vec3::new(5.0, 0.0, 0.0);
    assert!((predicted_transform(&fixed, 1.0).position - fixed.transform.position).norm() < 1e-15);

    let mut asleep = RigidBody::new(1.0);
    asleep.velocity = Vec3::new(5.0, 0.0, 0.0);
    asleep.state = BodyState::Sleeping;
    assert!(
        (predicted_transform(&asleep, 1.0).position - asleep.transform.position).norm() < 1e-15
    );

    let mut scripted = RigidBody::new_kinematic();
    scripted.velocity = Vec3::new(5.0, 0.0, 0.0);
    scripted.kinematic_target = Some(oxiphysics_rigid::KinematicTarget {
        target_position: Vec3::new(9.0, 0.0, 0.0),
        target_rotation: Quat::identity(),
    });
    assert!(
        (predicted_transform(&scripted, 1.0).position - scripted.transform.position).norm() < 1e-15
    );
}

#[test]
fn shapes_are_classified_for_the_narrowphase() {
    assert_eq!(classify_shape(&Sphere::new(1.0)), ShapeType::Sphere);
    assert_eq!(
        classify_shape(&BoxShape::new(Vec3::new(1.0, 1.0, 1.0))),
        ShapeType::Box
    );
    assert_eq!(classify_shape(&Capsule::new(0.25, 0.5)), ShapeType::Capsule);
    assert_eq!(
        classify_shape(&Cylinder::new(0.5, 0.5)),
        ShapeType::Cylinder
    );
    assert_eq!(classify_shape(&Cone::new(0.5, 0.5)), ShapeType::Cone);
    let hull = ConvexHull::new(vec![
        Vec3::new(0.0, 0.0, 0.0),
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
    ]);
    assert_eq!(classify_shape(&hull), ShapeType::ConvexHull);
}

#[test]
fn sphere_capsule_pair_has_a_b_to_a_normal_in_both_insertion_orders() {
    // A sphere resting on top of a capsule: the contact normal points from
    // B toward A whichever of the two was inserted first.
    for capsule_first in [true, false] {
        let mut scene = Scene::new();
        let mut capsule_body = RigidBody::new_static();
        capsule_body.transform = Transform::from_position(Vec3::zeros());
        let ball_body = |scene: &mut Scene| {
            scene.dynamic(1.0, Vec3::new(0.0, 1.2, 0.0), sphere(0.5), 0.5, 0.0);
        };
        if capsule_first {
            scene.add(capsule_body, Arc::new(Capsule::new(0.25, 0.5)), 0.5, 0.0);
            ball_body(&mut scene);
        } else {
            ball_body(&mut scene);
            scene.add(capsule_body, Arc::new(Capsule::new(0.25, 0.5)), 0.5, 0.0);
        }
        let pipeline = PhysicsPipeline::new();
        let entries = collect_entries(&scene.colliders);
        let detection = pipeline.detect_collisions(&entries, &scene.bodies);
        assert_eq!(detection.contacts.len(), 1, "capsule_first={capsule_first}");
        let record = &detection.contacts[0];
        let normal = record.manifold.contacts[0].normal;
        // The ball is "A" exactly when it was inserted first (index 0).
        let ball_is_a = !capsule_first;
        let expected_y = if ball_is_a { 1.0 } else { -1.0 };
        assert!(
            normal.y * expected_y > 0.9,
            "capsule_first={capsule_first}: normal {normal:?} should point from B to A"
        );
    }
}

// ---------------------------------------------------------------------------
// Sensors
// ---------------------------------------------------------------------------

/// A static sensor box spanning y in [0, 2] and a ball dropped from y = 3.
fn sensor_scene() -> (Scene, BodyHandle, ColliderHandle, ColliderHandle) {
    let mut scene = Scene::new();
    let mut sensor_body = RigidBody::new_static();
    sensor_body.transform = Transform::from_position(Vec3::new(0.0, 1.0, 0.0));
    let sensor_body = scene.bodies.insert(sensor_body);
    let sensor = scene.colliders.insert(
        Collider::new(cuboid(2.0, 1.0, 2.0))
            .with_body(sensor_body)
            .as_sensor(),
    );
    let (ball, ball_collider) =
        scene.dynamic(1.0, Vec3::new(0.0, 3.0, 0.0), sphere(0.25), 0.5, 0.0);
    (scene, ball, sensor, ball_collider)
}

fn run_sensor_fall(ccd_enabled: bool) {
    let config = PhysicsConfig {
        ccd_enabled,
        ..awake_config()
    };
    let mut pipeline = PhysicsPipeline::with_config(config);
    let (mut scene, ball, sensor, ball_collider) = sensor_scene();

    let g = 9.81;
    let mut overlap_steps = Vec::new();
    for n in 1..=90_u32 {
        scene.step(&mut pipeline);

        // A sensor exerts no force: the ball falls exactly as in free fall
        // (semi-implicit Euler: v_n = -g * n * dt, no damping).
        let expected_vy = -g * DT * f64::from(n);
        let vy = scene.velocity(ball).y;
        assert!(
            (vy - expected_vy).abs() < 1e-9,
            "ccd={ccd_enabled} step {n}: vy {vy} should be free fall {expected_vy}"
        );
        assert_eq!(scene.velocity(ball).x, 0.0);
        assert_eq!(scene.velocity(ball).z, 0.0);

        let overlaps = pipeline.last_sensor_overlaps();
        if !overlaps.is_empty() {
            overlap_steps.push(n);
            assert_eq!(overlaps, &[(sensor, ball_collider)]);
        }
        assert_eq!(
            pipeline.stats.last_contact_manifolds, 0,
            "sensor pairs are never solved"
        );
    }

    // The ball is inside the sensor from t = 0.39 s (step ~24) to 0.82 s
    // (step ~49), and nowhere else.
    let first = *overlap_steps
        .first()
        .expect("the ball must reach the sensor");
    let last = *overlap_steps
        .last()
        .expect("the ball must reach the sensor");
    assert!((22..=26).contains(&first), "first overlap at step {first}");
    assert!((47..=51).contains(&last), "last overlap at step {last}");
    assert_eq!(
        overlap_steps.len(),
        (last - first + 1) as usize,
        "the overlap is reported on every step in between"
    );
}

#[test]
fn sensor_does_not_stop_a_falling_sphere_with_ccd_off() {
    run_sensor_fall(false);
}

#[test]
fn sensor_does_not_stop_a_falling_sphere_with_ccd_on() {
    run_sensor_fall(true);
}

#[test]
fn sensor_overlaps_are_cleared_when_nothing_overlaps_any_more() {
    let mut pipeline = PhysicsPipeline::with_config(awake_config());
    let (mut scene, ball, _sensor, _ball_collider) = sensor_scene();
    // Put the ball inside the sensor: one overlap is reported ...
    if let Some(body) = scene.bodies.get_mut(ball) {
        body.transform.position = Vec3::new(0.0, 1.0, 0.0);
        body.gravity_scale = 0.0;
    }
    scene.step(&mut pipeline);
    assert_eq!(pipeline.last_sensor_overlaps().len(), 1);
    // ... and moving it out empties the list again.
    if let Some(body) = scene.bodies.get_mut(ball) {
        body.transform.position = Vec3::new(0.0, 50.0, 0.0);
    }
    scene.step(&mut pipeline);
    assert!(pipeline.last_sensor_overlaps().is_empty());
}

#[test]
fn sensor_with_a_solid_floor_does_not_disturb_normal_contacts() {
    let mut pipeline = PhysicsPipeline::with_config(awake_config());
    let (mut scene, ball, _sensor, _ball_collider) = sensor_scene();
    scene.static_box(
        Vec3::new(0.0, -5.5, 0.0),
        Vec3::new(10.0, 0.5, 10.0),
        0.5,
        0.0,
    );
    for _ in 0..240 {
        scene.step(&mut pipeline);
    }
    // The ball fell through the sensor and came to rest on the floor (top at
    // y = -5) with its centre one radius above it.
    let p = scene.position(ball);
    assert!((p.y - (-4.75)).abs() < 0.02, "ball rests at y = {}", p.y);
    assert!(scene.velocity(ball).norm() < 1e-3);
}

#[test]
fn static_sensor_overlapping_static_geometry_is_not_reported() {
    let mut pipeline = PhysicsPipeline::with_config(awake_config());
    let mut scene = Scene::new();
    let mut sensor_body = RigidBody::new_static();
    sensor_body.transform = Transform::from_position(Vec3::zeros());
    let sensor_body = scene.bodies.insert(sensor_body);
    scene.colliders.insert(
        Collider::new(cuboid(1.0, 1.0, 1.0))
            .with_body(sensor_body)
            .as_sensor(),
    );
    scene.static_box(Vec3::zeros(), Vec3::new(1.0, 1.0, 1.0), 0.5, 0.0);

    scene.step(&mut pipeline);
    assert!(
        pipeline.last_sensor_overlaps().is_empty(),
        "nothing can enter a static-static overlap"
    );
}

#[test]
fn colliders_of_one_body_never_collide_with_each_other() {
    let mut pipeline = PhysicsPipeline::with_config(awake_config());
    let mut scene = Scene::new();
    let mut body = RigidBody::new(2.0);
    body.gravity_scale = 0.0;
    let body_handle = scene.bodies.insert(body);
    for x in [-0.1, 0.1] {
        scene.colliders.insert(
            Collider::new(sphere(0.5))
                .with_body(body_handle)
                .with_local_transform(Transform::from_position(Vec3::new(x, 0.0, 0.0))),
        );
    }
    for _ in 0..10 {
        scene.step(&mut pipeline);
    }
    assert_eq!(pipeline.stats.last_contact_manifolds, 0);
    assert_eq!(scene.velocity(body_handle), Vec3::zeros());
    assert_eq!(scene.position(body_handle), Vec3::zeros());
}

#[test]
fn collider_whose_body_was_removed_acts_as_world_geometry() {
    let mut pipeline = PhysicsPipeline::with_config(awake_config());
    let mut scene = Scene::new();
    // A floor collider whose body is gone: it keeps its local transform only.
    let stale_body = scene.bodies.insert(RigidBody::new_static());
    scene.colliders.insert(
        Collider::new(cuboid(10.0, 0.5, 10.0))
            .with_body(stale_body)
            .with_local_transform(Transform::from_position(Vec3::new(0.0, -0.5, 0.0))),
    );
    scene.bodies.remove(stale_body);
    let (ball, _) = scene.dynamic(1.0, Vec3::new(0.0, 0.6, 0.0), sphere(0.5), 0.5, 0.0);

    for _ in 0..120 {
        scene.step(&mut pipeline);
    }
    let p = scene.position(ball);
    assert!(
        (p.y - 0.5).abs() < 0.02,
        "ball rests on the floor: y = {}",
        p.y
    );
}
