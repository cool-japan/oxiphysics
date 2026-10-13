// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tests: ray casting against every collider shape.
//!
//! The ray caster used to know spheres only and to miss a sphere whenever the
//! ray started inside it.  These tests cover spheres, oriented boxes,
//! capsules, cylinders and planes, rays that start inside a solid, the
//! `max_distance` cut-off and invalid input.

#![cfg(test)]

use std::f64::consts::{FRAC_1_SQRT_2, FRAC_PI_2, FRAC_PI_4, SQRT_2};

use super::geom::quat_from_axis_angle;
use super::*;

const IDENTITY: [f64; 4] = [0.0, 0.0, 0.0, 1.0];

fn engine() -> WasmPhysicsEngine {
    WasmPhysicsEngine::new(0.0, 0.0, 0.0)
}

/// Add a static body at `position`, turned by `rotation`, carrying `collider`.
fn add_shape(
    engine: &mut WasmPhysicsEngine,
    position: [f64; 3],
    rotation: [f64; 4],
    collider: ColliderConfig,
) -> u32 {
    let body = engine
        .add_rigid_body(&RigidBodyConfig {
            rotation,
            ..RigidBodyConfig::static_body().with_position(position[0], position[1], position[2])
        })
        .expect("static body");
    engine.add_collider(body, &collider).expect("collider");
    body
}

fn assert_hit(result: &RaycastResult, body: u32, distance: f64, point: [f64; 3], normal: [f64; 3]) {
    assert!(result.hit, "expected a hit, got {result:?}");
    assert_eq!(result.body_handle, body, "body handle");
    assert!(
        (result.distance - distance).abs() < 1e-9,
        "distance: expected {distance}, got {}",
        result.distance
    );
    for axis in 0..3 {
        assert!(
            (result.point[axis] - point[axis]).abs() < 1e-9,
            "point[{axis}]: expected {}, got {}",
            point[axis],
            result.point[axis]
        );
        assert!(
            (result.normal[axis] - normal[axis]).abs() < 1e-9,
            "normal[{axis}]: expected {}, got {}",
            normal[axis],
            result.normal[axis]
        );
    }
}

fn cast(
    engine: &WasmPhysicsEngine,
    origin: [f64; 3],
    direction: [f64; 3],
    max_distance: f64,
) -> RaycastResult {
    engine.raycast(
        origin[0],
        origin[1],
        origin[2],
        direction[0],
        direction[1],
        direction[2],
        max_distance,
    )
}

// ---------------------------------------------------------------------------
// Spheres
// ---------------------------------------------------------------------------

#[test]
fn ray_hits_the_near_side_of_a_sphere() {
    let mut e = engine();
    let body = add_shape(
        &mut e,
        [10.0, 0.0, 0.0],
        IDENTITY,
        ColliderConfig::sphere(1.0),
    );
    let hit = cast(&e, [0.0; 3], [1.0, 0.0, 0.0], 100.0);
    assert_hit(&hit, body, 9.0, [9.0, 0.0, 0.0], [-1.0, 0.0, 0.0]);
}

#[test]
fn ray_starting_inside_a_sphere_reports_the_exit_surface() {
    let mut e = engine();
    let body = add_shape(&mut e, [0.0; 3], IDENTITY, ColliderConfig::sphere(2.0));
    // From the centre.
    let hit = cast(&e, [0.0; 3], [1.0, 0.0, 0.0], 100.0);
    assert_hit(&hit, body, 2.0, [2.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
    // Off centre, towards the near wall: the ray leaves through the far side.
    let hit = cast(&e, [0.5, 0.0, 0.0], [-1.0, 0.0, 0.0], 100.0);
    assert_hit(&hit, body, 2.5, [-2.0, 0.0, 0.0], [-1.0, 0.0, 0.0]);
    // The exit obeys `max_distance` as well.
    assert!(!cast(&e, [0.0; 3], [1.0, 0.0, 0.0], 1.9).hit);
}

#[test]
fn ray_on_the_surface_moving_inward_reports_the_exit_and_moving_outward_misses() {
    let mut e = engine();
    let body = add_shape(&mut e, [0.0; 3], IDENTITY, ColliderConfig::sphere(1.0));
    let inward = cast(&e, [1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], 10.0);
    assert_hit(&inward, body, 2.0, [-1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]);
    assert!(!cast(&e, [1.0, 0.0, 0.0], [1.0, 0.0, 0.0], 10.0).hit);
}

#[test]
fn sphere_behind_the_ray_origin_is_not_hit() {
    let mut e = engine();
    add_shape(
        &mut e,
        [-5.0, 0.0, 0.0],
        IDENTITY,
        ColliderConfig::sphere(1.0),
    );
    assert!(!cast(&e, [0.0; 3], [1.0, 0.0, 0.0], 100.0).hit);
}

#[test]
fn ray_passing_beside_a_sphere_misses() {
    let mut e = engine();
    add_shape(
        &mut e,
        [5.0, 0.0, 0.0],
        IDENTITY,
        ColliderConfig::sphere(1.0),
    );
    assert!(!cast(&e, [0.0, 1.5, 0.0], [1.0, 0.0, 0.0], 100.0).hit);
    assert!(cast(&e, [0.0, 0.9, 0.0], [1.0, 0.0, 0.0], 100.0).hit);
}

// ---------------------------------------------------------------------------
// Boxes
// ---------------------------------------------------------------------------

#[test]
fn ray_hits_an_axis_aligned_box_face() {
    let mut e = engine();
    let body = add_shape(
        &mut e,
        [5.0, 0.0, 0.0],
        IDENTITY,
        ColliderConfig::cuboid(1.0, 1.0, 1.0),
    );
    let hit = cast(&e, [0.0; 3], [1.0, 0.0, 0.0], 100.0);
    assert_hit(&hit, body, 4.0, [4.0, 0.0, 0.0], [-1.0, 0.0, 0.0]);
}

#[test]
fn ray_hits_every_face_of_a_box() {
    let mut e = engine();
    let body = add_shape(
        &mut e,
        [0.0; 3],
        IDENTITY,
        ColliderConfig::cuboid(1.0, 2.0, 3.0),
    );
    let cases: [([f64; 3], [f64; 3], f64); 6] = [
        ([10.0, 0.5, 0.5], [-1.0, 0.0, 0.0], 9.0),
        ([-10.0, 0.5, 0.5], [1.0, 0.0, 0.0], 9.0),
        ([0.5, 10.0, 0.5], [0.0, -1.0, 0.0], 8.0),
        ([0.5, -10.0, 0.5], [0.0, 1.0, 0.0], 8.0),
        ([0.5, 0.5, 10.0], [0.0, 0.0, -1.0], 7.0),
        ([0.5, 0.5, -10.0], [0.0, 0.0, 1.0], 7.0),
    ];
    for (origin, direction, distance) in cases {
        let hit = cast(&e, origin, direction, 100.0);
        let normal = [-direction[0], -direction[1], -direction[2]];
        let point = [
            origin[0] + direction[0] * distance,
            origin[1] + direction[1] * distance,
            origin[2] + direction[2] * distance,
        ];
        assert_hit(&hit, body, distance, point, normal);
    }
}

#[test]
fn ray_hits_a_rotated_box_on_its_tilted_face() {
    // A unit box turned by 45 degrees about z is a diamond with its corner
    // pointing at the origin; the ray at y = 0.5 meets the upper-left face
    // `-(x - 5) + y = sqrt(2)`.
    let mut e = engine();
    let rotation = quat_from_axis_angle([0.0, 0.0, 1.0], FRAC_PI_4);
    let body = add_shape(
        &mut e,
        [5.0, 0.0, 0.0],
        rotation,
        ColliderConfig::cuboid(1.0, 1.0, 1.0),
    );
    let hit = cast(&e, [0.0, 0.5, 0.0], [1.0, 0.0, 0.0], 100.0);
    let x = 5.0 + 0.5 - SQRT_2;
    assert_hit(
        &hit,
        body,
        x,
        [x, 0.5, 0.0],
        [-FRAC_1_SQRT_2, FRAC_1_SQRT_2, 0.0],
    );
    // The lower half of the diamond has the mirrored normal.
    let hit = cast(&e, [0.0, -0.5, 0.0], [1.0, 0.0, 0.0], 100.0);
    assert_hit(
        &hit,
        body,
        x,
        [x, -0.5, 0.0],
        [-FRAC_1_SQRT_2, -FRAC_1_SQRT_2, 0.0],
    );
}

#[test]
fn box_collider_follows_the_rotation_of_its_body() {
    // A long thin box lying along x becomes tall after a quarter turn.
    let mut e = engine();
    let rotation = quat_from_axis_angle([0.0, 0.0, 1.0], FRAC_PI_2);
    let body = add_shape(
        &mut e,
        [5.0, 0.0, 0.0],
        rotation,
        ColliderConfig::cuboid(2.0, 0.25, 0.25),
    );
    // Unturned the box would span x in [3, 7]; turned it spans x in
    // [4.75, 5.25] and y in [-2, 2].
    let hit = cast(&e, [0.0, 1.0, 0.0], [1.0, 0.0, 0.0], 100.0);
    assert_hit(&hit, body, 4.75, [4.75, 1.0, 0.0], [-1.0, 0.0, 0.0]);
    assert!(cast(&e, [0.0, 1.9, 0.0], [1.0, 0.0, 0.0], 100.0).hit);
    assert!(!cast(&e, [0.0, 2.1, 0.0], [1.0, 0.0, 0.0], 100.0).hit);
}

#[test]
fn ray_inside_a_box_reports_the_exit_face() {
    let mut e = engine();
    let body = add_shape(
        &mut e,
        [0.0; 3],
        IDENTITY,
        ColliderConfig::cuboid(1.0, 1.0, 1.0),
    );
    let hit = cast(&e, [0.0; 3], [0.0, 1.0, 0.0], 100.0);
    assert_hit(&hit, body, 1.0, [0.0, 1.0, 0.0], [0.0, 1.0, 0.0]);
}

#[test]
fn ray_parallel_to_a_box_slab_outside_it_misses() {
    let mut e = engine();
    add_shape(
        &mut e,
        [0.0; 3],
        IDENTITY,
        ColliderConfig::cuboid(1.0, 1.0, 1.0),
    );
    assert!(!cast(&e, [-5.0, 1.5, 0.0], [1.0, 0.0, 0.0], 100.0).hit);
    assert!(cast(&e, [-5.0, 0.99, 0.0], [1.0, 0.0, 0.0], 100.0).hit);
}

// ---------------------------------------------------------------------------
// Capsules
// ---------------------------------------------------------------------------

#[test]
fn ray_hits_the_side_of_a_capsule() {
    // Radius 0.5 and a 2 m cylindrical part along y: the caps are centred at
    // y = +-1.
    let mut e = engine();
    let body = add_shape(
        &mut e,
        [5.0, 0.0, 0.0],
        IDENTITY,
        ColliderConfig::capsule(0.5, 2.0),
    );
    let hit = cast(&e, [0.0, 0.3, 0.0], [1.0, 0.0, 0.0], 100.0);
    assert_hit(&hit, body, 4.5, [4.5, 0.3, 0.0], [-1.0, 0.0, 0.0]);
}

#[test]
fn ray_hits_the_rounded_cap_of_a_capsule() {
    let mut e = engine();
    let body = add_shape(
        &mut e,
        [5.0, 0.0, 0.0],
        IDENTITY,
        ColliderConfig::capsule(0.5, 2.0),
    );
    // At y = 1.3 the ray crosses the upper cap sphere (centre (5, 1, 0)) 0.3
    // above its centre: dx = sqrt(0.25 - 0.09) = 0.4.
    let hit = cast(&e, [0.0, 1.3, 0.0], [1.0, 0.0, 0.0], 100.0);
    assert_hit(&hit, body, 4.6, [4.6, 1.3, 0.0], [-0.8, 0.6, 0.0]);
    // Straight down onto the top of the capsule.
    let hit = cast(&e, [5.0, 5.0, 0.0], [0.0, -1.0, 0.0], 100.0);
    assert_hit(&hit, body, 3.5, [5.0, 1.5, 0.0], [0.0, 1.0, 0.0]);
    // Above the cap: miss.
    assert!(!cast(&e, [0.0, 1.6, 0.0], [1.0, 0.0, 0.0], 100.0).hit);
}

#[test]
fn capsule_axis_follows_the_rotation_of_its_body() {
    // Turned by 90 degrees about z the axis lies along x: its end cap points
    // at the origin and reaches x = 5 - 1 - 0.5 = 3.5.
    let mut e = engine();
    let rotation = quat_from_axis_angle([0.0, 0.0, 1.0], FRAC_PI_2);
    let body = add_shape(
        &mut e,
        [5.0, 0.0, 0.0],
        rotation,
        ColliderConfig::capsule(0.5, 2.0),
    );
    let hit = cast(&e, [0.0; 3], [1.0, 0.0, 0.0], 100.0);
    assert_hit(&hit, body, 3.5, [3.5, 0.0, 0.0], [-1.0, 0.0, 0.0]);
}

#[test]
fn ray_inside_a_capsule_reports_the_exit_surface() {
    let mut e = engine();
    let body = add_shape(
        &mut e,
        [0.0; 3],
        IDENTITY,
        ColliderConfig::capsule(0.5, 2.0),
    );
    let hit = cast(&e, [0.0; 3], [0.0, 1.0, 0.0], 100.0);
    assert_hit(&hit, body, 1.5, [0.0, 1.5, 0.0], [0.0, 1.0, 0.0]);
    let hit = cast(&e, [0.0; 3], [1.0, 0.0, 0.0], 100.0);
    assert_hit(&hit, body, 0.5, [0.5, 0.0, 0.0], [1.0, 0.0, 0.0]);
}

// ---------------------------------------------------------------------------
// Cylinders
// ---------------------------------------------------------------------------

#[test]
fn ray_hits_the_side_and_the_caps_of_a_cylinder() {
    let mut e = engine();
    let body = add_shape(
        &mut e,
        [5.0, 0.0, 0.0],
        IDENTITY,
        ColliderConfig::cylinder(0.5, 2.0),
    );
    let side = cast(&e, [0.0; 3], [1.0, 0.0, 0.0], 100.0);
    assert_hit(&side, body, 4.5, [4.5, 0.0, 0.0], [-1.0, 0.0, 0.0]);

    let top = cast(&e, [5.0, 5.0, 0.0], [0.0, -1.0, 0.0], 100.0);
    assert_hit(&top, body, 4.0, [5.0, 1.0, 0.0], [0.0, 1.0, 0.0]);

    let bottom = cast(&e, [5.0, -5.0, 0.0], [0.0, 1.0, 0.0], 100.0);
    assert_hit(&bottom, body, 4.0, [5.0, -1.0, 0.0], [0.0, -1.0, 0.0]);
}

#[test]
fn cylinder_has_sharp_edges_unlike_a_capsule() {
    let mut e = engine();
    add_shape(
        &mut e,
        [5.0, 0.0, 0.0],
        IDENTITY,
        ColliderConfig::cylinder(0.5, 2.0),
    );
    // Just above the cap plane: no hit; just below: a side hit.
    assert!(!cast(&e, [0.0, 1.05, 0.0], [1.0, 0.0, 0.0], 100.0).hit);
    assert!(cast(&e, [0.0, 0.95, 0.0], [1.0, 0.0, 0.0], 100.0).hit);
}

#[test]
fn ray_entering_the_top_cap_off_centre_hits_the_cap() {
    let mut e = engine();
    let body = add_shape(
        &mut e,
        [0.0; 3],
        IDENTITY,
        ColliderConfig::cylinder(1.0, 2.0),
    );
    let hit = cast(&e, [-0.5, 3.0, 0.0], [0.0, -1.0, 0.0], 100.0);
    assert_hit(&hit, body, 2.0, [-0.5, 1.0, 0.0], [0.0, 1.0, 0.0]);
}

#[test]
fn ray_inside_a_cylinder_reports_the_exit_surface() {
    let mut e = engine();
    let body = add_shape(
        &mut e,
        [0.0; 3],
        IDENTITY,
        ColliderConfig::cylinder(1.0, 2.0),
    );
    let hit = cast(&e, [0.0; 3], [0.0, -1.0, 0.0], 100.0);
    assert_hit(&hit, body, 1.0, [0.0, -1.0, 0.0], [0.0, -1.0, 0.0]);
    let hit = cast(&e, [0.0; 3], [0.0, 0.0, 1.0], 100.0);
    assert_hit(&hit, body, 1.0, [0.0, 0.0, 1.0], [0.0, 0.0, 1.0]);
}

#[test]
fn ray_along_the_cylinder_axis_outside_the_radius_misses() {
    let mut e = engine();
    add_shape(
        &mut e,
        [0.0; 3],
        IDENTITY,
        ColliderConfig::cylinder(1.0, 2.0),
    );
    assert!(!cast(&e, [1.5, -5.0, 0.0], [0.0, 1.0, 0.0], 100.0).hit);
    assert!(cast(&e, [0.9, -5.0, 0.0], [0.0, 1.0, 0.0], 100.0).hit);
}

// ---------------------------------------------------------------------------
// Planes
// ---------------------------------------------------------------------------

#[test]
fn ray_hits_a_plane_from_the_free_side() {
    let mut e = engine();
    let ground = add_shape(
        &mut e,
        [0.0; 3],
        IDENTITY,
        ColliderConfig::plane(0.0, 1.0, 0.0, 0.0),
    );
    let down = cast(&e, [0.0, 5.0, 0.0], [0.0, -1.0, 0.0], 100.0);
    assert_hit(&down, ground, 5.0, [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]);

    let oblique = cast(&e, [0.0, 5.0, 0.0], [1.0, -1.0, 0.0], 100.0);
    assert_hit(
        &oblique,
        ground,
        5.0 * SQRT_2,
        [5.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
    );
}

#[test]
fn planes_are_one_sided() {
    let mut e = engine();
    add_shape(
        &mut e,
        [0.0; 3],
        IDENTITY,
        ColliderConfig::plane(0.0, 1.0, 0.0, 0.0),
    );
    // From below, travelling with the normal: no hit.
    assert!(!cast(&e, [0.0, -5.0, 0.0], [0.0, 1.0, 0.0], 100.0).hit);
    // Parallel to the plane: no hit.
    assert!(!cast(&e, [0.0, 5.0, 0.0], [1.0, 0.0, 0.0], 100.0).hit);
    // Pointing away from the plane: no hit.
    assert!(!cast(&e, [0.0, 5.0, 0.0], [0.0, 1.0, 0.0], 100.0).hit);
}

#[test]
fn plane_offset_translation_and_rotation_are_applied() {
    let mut e = engine();
    // Offset: the plane y = 2.
    let offset_plane = add_shape(
        &mut e,
        [0.0; 3],
        IDENTITY,
        ColliderConfig::plane(0.0, 1.0, 0.0, 2.0),
    );
    let hit = cast(&e, [0.0, 5.0, 0.0], [0.0, -1.0, 0.0], 100.0);
    assert_hit(&hit, offset_plane, 3.0, [0.0, 2.0, 0.0], [0.0, 1.0, 0.0]);

    // Translation: the plane through (0, 3, 0) of a second body is hit first.
    let mut e = engine();
    let raised = add_shape(
        &mut e,
        [0.0, 3.0, 0.0],
        IDENTITY,
        ColliderConfig::plane(0.0, 1.0, 0.0, 0.0),
    );
    let hit = cast(&e, [0.0, 5.0, 0.0], [0.0, -1.0, 0.0], 100.0);
    assert_hit(&hit, raised, 2.0, [0.0, 3.0, 0.0], [0.0, 1.0, 0.0]);

    // Rotation: turned by 90 degrees about z the normal points along -x.
    let mut e = engine();
    let rotation = quat_from_axis_angle([0.0, 0.0, 1.0], FRAC_PI_2);
    let wall = add_shape(
        &mut e,
        [0.0; 3],
        rotation,
        ColliderConfig::plane(0.0, 1.0, 0.0, 0.0),
    );
    let hit = cast(&e, [-5.0, 0.0, 0.0], [1.0, 0.0, 0.0], 100.0);
    assert_hit(&hit, wall, 5.0, [0.0, 0.0, 0.0], [-1.0, 0.0, 0.0]);
    assert!(!cast(&e, [5.0, 0.0, 0.0], [-1.0, 0.0, 0.0], 100.0).hit);
}

// ---------------------------------------------------------------------------
// Several colliders, limits and invalid input
// ---------------------------------------------------------------------------

#[test]
fn nearest_hit_wins_across_shapes() {
    let mut e = engine();
    add_shape(
        &mut e,
        [8.0, 0.0, 0.0],
        IDENTITY,
        ColliderConfig::sphere(1.0),
    );
    let near_box = add_shape(
        &mut e,
        [3.0, 0.0, 0.0],
        IDENTITY,
        ColliderConfig::cuboid(0.5, 0.5, 0.5),
    );
    let behind = add_shape(
        &mut e,
        [-4.0, 0.0, 0.0],
        IDENTITY,
        ColliderConfig::sphere(1.0),
    );
    let hit = cast(&e, [0.0; 3], [1.0, 0.0, 0.0], 100.0);
    assert_hit(&hit, near_box, 2.5, [2.5, 0.0, 0.0], [-1.0, 0.0, 0.0]);
    // From the other side the sphere behind is nearest.
    let hit = cast(&e, [0.0; 3], [-1.0, 0.0, 0.0], 100.0);
    assert_hit(&hit, behind, 3.0, [-3.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
}

#[test]
fn compound_body_reports_the_body_of_the_collider_that_was_hit() {
    let mut e = engine();
    let body = e.add_static_body(0.0, 0.0, 0.0).expect("body");
    e.add_collider(
        body,
        &ColliderConfig {
            local_position: [4.0, 0.0, 0.0],
            ..ColliderConfig::sphere(0.5)
        },
    )
    .expect("sphere collider");
    e.add_collider(
        body,
        &ColliderConfig {
            local_position: [8.0, 0.0, 0.0],
            ..ColliderConfig::cuboid(0.5, 0.5, 0.5)
        },
    )
    .expect("box collider");
    let hit = cast(&e, [0.0; 3], [1.0, 0.0, 0.0], 100.0);
    assert_hit(&hit, body, 3.5, [3.5, 0.0, 0.0], [-1.0, 0.0, 0.0]);
    let hit = cast(&e, [6.0, 0.0, 0.0], [1.0, 0.0, 0.0], 100.0);
    assert_hit(&hit, body, 1.5, [7.5, 0.0, 0.0], [-1.0, 0.0, 0.0]);
}

#[test]
fn collider_local_pose_is_part_of_the_ray_test() {
    let mut e = engine();
    let body = e.add_static_body(0.0, 0.0, 0.0).expect("body");
    // A capsule offset by 5 along x and turned so that its axis lies along x.
    e.add_collider(
        body,
        &ColliderConfig {
            local_position: [5.0, 0.0, 0.0],
            local_rotation: quat_from_axis_angle([0.0, 0.0, 1.0], FRAC_PI_2),
            ..ColliderConfig::capsule(0.5, 2.0)
        },
    )
    .expect("capsule collider");
    let hit = cast(&e, [0.0; 3], [1.0, 0.0, 0.0], 100.0);
    assert_hit(&hit, body, 3.5, [3.5, 0.0, 0.0], [-1.0, 0.0, 0.0]);
}

#[test]
fn max_distance_cuts_off_hits() {
    let mut e = engine();
    add_shape(
        &mut e,
        [5.0, 0.0, 0.0],
        IDENTITY,
        ColliderConfig::cuboid(1.0, 1.0, 1.0),
    );
    assert!(!cast(&e, [0.0; 3], [1.0, 0.0, 0.0], 3.9).hit);
    assert!(cast(&e, [0.0; 3], [1.0, 0.0, 0.0], 4.1).hit);
    assert!(cast(&e, [0.0; 3], [1.0, 0.0, 0.0], f64::INFINITY).hit);
}

#[test]
fn ray_direction_need_not_be_normalised() {
    let mut e = engine();
    let body = add_shape(
        &mut e,
        [5.0, 0.0, 0.0],
        IDENTITY,
        ColliderConfig::cuboid(1.0, 1.0, 1.0),
    );
    let hit = cast(&e, [0.0; 3], [3.0, 0.0, 0.0], 100.0);
    assert_hit(&hit, body, 4.0, [4.0, 0.0, 0.0], [-1.0, 0.0, 0.0]);
    let hit = cast(&e, [0.0; 3], [1e-3, 0.0, 0.0], 100.0);
    assert_hit(&hit, body, 4.0, [4.0, 0.0, 0.0], [-1.0, 0.0, 0.0]);
}

#[test]
fn invalid_rays_report_no_hit() {
    let mut e = engine();
    add_shape(
        &mut e,
        [5.0, 0.0, 0.0],
        IDENTITY,
        ColliderConfig::sphere(1.0),
    );
    assert!(!cast(&e, [0.0; 3], [0.0; 3], 100.0).hit, "zero direction");
    assert!(!cast(&e, [f64::NAN, 0.0, 0.0], [1.0, 0.0, 0.0], 100.0).hit);
    assert!(!cast(&e, [0.0; 3], [f64::NAN, 0.0, 0.0], 100.0).hit);
    assert!(!cast(&e, [0.0; 3], [f64::INFINITY, 0.0, 0.0], 100.0).hit);
    assert!(!cast(&e, [0.0; 3], [1.0, 0.0, 0.0], 0.0).hit);
    assert!(!cast(&e, [0.0; 3], [1.0, 0.0, 0.0], -5.0).hit);
    assert!(!cast(&e, [0.0; 3], [1.0, 0.0, 0.0], f64::NAN).hit);
    // The same ray with valid parameters does hit.
    assert!(cast(&e, [0.0; 3], [1.0, 0.0, 0.0], 100.0).hit);
}

#[test]
fn no_hit_result_is_the_documented_sentinel() {
    let e = engine();
    let miss = cast(&e, [0.0; 3], [1.0, 0.0, 0.0], 100.0);
    assert!(!miss.hit);
    assert_eq!(miss.body_handle, u32::MAX);
    assert!(miss.distance.is_infinite());
}

#[test]
fn removed_bodies_are_not_hit() {
    let mut e = engine();
    let body = add_shape(
        &mut e,
        [5.0, 0.0, 0.0],
        IDENTITY,
        ColliderConfig::sphere(1.0),
    );
    assert!(cast(&e, [0.0; 3], [1.0, 0.0, 0.0], 100.0).hit);
    e.remove_body(body).expect("remove body");
    assert!(!cast(&e, [0.0; 3], [1.0, 0.0, 0.0], 100.0).hit);
}

#[test]
fn removing_a_collider_removes_it_from_ray_tests() {
    let mut e = engine();
    let body = e.add_static_body(5.0, 0.0, 0.0).expect("body");
    let collider = e.add_sphere_collider(body, 1.0).expect("collider");
    assert!(cast(&e, [0.0; 3], [1.0, 0.0, 0.0], 100.0).hit);
    e.remove_collider(collider).expect("remove collider");
    assert!(!cast(&e, [0.0; 3], [1.0, 0.0, 0.0], 100.0).hit);
}

#[test]
fn rays_follow_turning_bodies() {
    // A long box along x spins a quarter turn about z in one second (60
    // substeps) and becomes tall.
    let mut e = WasmPhysicsEngine::from_config(SimulationConfig {
        gravity: [0.0; 3],
        max_substeps: 100,
        ..SimulationConfig::default()
    })
    .expect("engine");
    let body = e
        .add_rigid_body(&RigidBodyConfig {
            angular_velocity: [0.0, 0.0, FRAC_PI_2],
            linear_damping: 0.0,
            angular_damping: 0.0,
            ..RigidBodyConfig::new().with_position(5.0, 0.0, 0.0)
        })
        .expect("spinning body");
    e.add_box_collider(body, 2.0, 0.25, 0.25)
        .expect("box collider");

    // Before the turn the ray at y = 1 passes above the thin box.
    assert!(!cast(&e, [0.0, 1.0, 0.0], [1.0, 0.0, 0.0], 100.0).hit);
    e.step(1.0).expect("one second");
    let hit = cast(&e, [0.0, 1.0, 0.0], [1.0, 0.0, 0.0], 100.0);
    assert_hit(&hit, body, 4.75, [4.75, 1.0, 0.0], [-1.0, 0.0, 0.0]);
}

#[test]
fn rays_follow_moving_bodies() {
    let mut e = WasmPhysicsEngine::new(0.0, 0.0, 0.0);
    let body = e.add_dynamic_body(1.0, 5.0, 0.0, 0.0).expect("body");
    e.add_sphere_collider(body, 1.0).expect("collider");
    e.set_body_linear_damping(body, 0.0).expect("damping");
    e.set_velocity(body, 0.0, 6.0, 0.0).expect("velocity");
    assert!(cast(&e, [0.0; 3], [1.0, 0.0, 0.0], 100.0).hit);
    // After 0.5 s (30 substeps would exceed the default cap, so step in
    // frames) the sphere has moved to y = 3 and the ray passes below it.
    for _ in 0..30 {
        e.step(1.0 / 60.0).expect("frame");
    }
    assert!(!cast(&e, [0.0; 3], [1.0, 0.0, 0.0], 100.0).hit);
    let hit = cast(&e, [0.0, 3.0, 0.0], [1.0, 0.0, 0.0], 100.0);
    assert!(
        (hit.distance - 4.0).abs() < 1e-6,
        "sphere now at y = 3: {hit:?}"
    );
}

#[test]
fn raycast_js_returns_handle_and_distance_or_an_empty_array() {
    let mut e = engine();
    let body = add_shape(
        &mut e,
        [5.0, 0.0, 0.0],
        IDENTITY,
        ColliderConfig::sphere(1.0),
    );
    let hit = e.raycast_js(0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 100.0);
    assert_eq!(hit.len(), 2);
    assert_eq!(hit[0], f64::from(body));
    assert!((hit[1] - 4.0).abs() < 1e-9);
    assert!(e.raycast_js(0.0, 5.0, 0.0, 1.0, 0.0, 0.0, 100.0).is_empty());
}
