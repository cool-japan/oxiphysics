// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tests: narrow phase and contact response.
//!
//! The sphere-sphere contact normal used to point from A toward B, the
//! opposite of the documented (and plane) convention, so approaching spheres
//! received no impulse and were pushed *into* each other; boxes, capsules and
//! cylinders did not collide at all.  The tests below pin the contact
//! convention (`normal` points from B toward A), check the contact geometry of
//! every shape pair analytically and with randomised invariants, and run the
//! scenarios from the audit through the full engine.

#![cfg(test)]

use std::f64::consts::{FRAC_1_SQRT_2, FRAC_PI_2, FRAC_PI_4, SQRT_2};

use super::geom::{Vec3, add, dot, length, madd, quat_from_axis_angle, quat_to_mat3, scale, sub};
use super::narrowphase::{Contact, collide};
use super::shapes::{Shape, WorldShape};
use super::*;

const IDENTITY: [f64; 4] = [0.0, 0.0, 0.0, 1.0];
const DT: f64 = 1.0 / 60.0;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn placed(shape: Shape, position: Vec3, rotation: [f64; 4]) -> WorldShape {
    WorldShape::place(
        shape,
        position,
        &quat_to_mat3(rotation),
        [0.0; 3],
        &quat_to_mat3(IDENTITY),
    )
}

fn sphere(radius: f64) -> Shape {
    Shape::Sphere { radius }
}

fn capsule(radius: f64, half_height: f64) -> Shape {
    Shape::Capsule {
        radius,
        half_height,
    }
}

fn cuboid(hx: f64, hy: f64, hz: f64) -> Shape {
    Shape::Box {
        half_extents: [hx, hy, hz],
    }
}

fn cylinder(radius: f64, half_height: f64) -> Shape {
    Shape::Cylinder {
        radius,
        half_height,
    }
}

fn plane(normal: Vec3, offset: f64) -> Shape {
    Shape::Plane { normal, offset }
}

fn assert_close(actual: f64, expected: f64, tolerance: f64, what: &str) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{what}: expected {expected}, got {actual} (tolerance {tolerance})"
    );
}

fn assert_vec_close(actual: Vec3, expected: Vec3, tolerance: f64, what: &str) {
    for axis in 0..3 {
        assert_close(
            actual[axis],
            expected[axis],
            tolerance,
            &format!("{what}[{axis}]"),
        );
    }
}

fn expect_contact(a: &WorldShape, b: &WorldShape) -> Contact {
    collide(a, b).unwrap_or_else(|| panic!("expected a contact between {a:?} and {b:?}"))
}

// ---------------------------------------------------------------------------
// Sphere pairs (F5665)
// ---------------------------------------------------------------------------

#[test]
fn sphere_sphere_normal_points_from_b_toward_a() {
    let a = placed(sphere(1.0), [0.0, 0.0, 0.0], IDENTITY);
    let b = placed(sphere(1.0), [1.5, 0.0, 0.0], IDENTITY);
    let contact = expect_contact(&a, &b);
    // A sits at the smaller x, so A must be pushed toward -x.
    assert_vec_close(contact.normal, [-1.0, 0.0, 0.0], 1e-12, "normal");
    assert_close(contact.depth, 0.5, 1e-12, "depth");
    assert_vec_close(contact.point_a, [1.0, 0.0, 0.0], 1e-12, "point on A");
    assert_vec_close(contact.point_b, [0.5, 0.0, 0.0], 1e-12, "point on B");
    // Swapping the roles flips the normal and keeps the depth.
    let mirrored = expect_contact(&b, &a);
    assert_vec_close(mirrored.normal, [1.0, 0.0, 0.0], 1e-12, "mirrored normal");
    assert_close(mirrored.depth, 0.5, 1e-12, "mirrored depth");
}

#[test]
fn sphere_sphere_touching_or_apart_has_no_contact() {
    let a = placed(sphere(1.0), [0.0; 3], IDENTITY);
    assert!(collide(&a, &placed(sphere(1.0), [2.0, 0.0, 0.0], IDENTITY)).is_none());
    assert!(collide(&a, &placed(sphere(1.0), [2.5, 0.0, 0.0], IDENTITY)).is_none());
    assert!(collide(&a, &placed(sphere(1.0), [1.999, 0.0, 0.0], IDENTITY)).is_some());
}

#[test]
fn concentric_spheres_still_produce_a_unit_normal() {
    let a = placed(sphere(1.0), [0.0; 3], IDENTITY);
    let b = placed(sphere(0.5), [0.0; 3], IDENTITY);
    let contact = expect_contact(&a, &b);
    assert_close(length(contact.normal), 1.0, 1e-12, "normal length");
    assert_close(contact.depth, 1.5, 1e-12, "depth");
}

#[test]
fn sphere_capsule_contacts_the_side_and_the_cap() {
    let cap = placed(capsule(0.5, 1.0), [0.0; 3], IDENTITY);
    // Beside the cylindrical part.
    let side = expect_contact(&placed(sphere(0.5), [0.8, 0.3, 0.0], IDENTITY), &cap);
    assert_vec_close(side.normal, [1.0, 0.0, 0.0], 1e-12, "side normal");
    assert_close(side.depth, 0.2, 1e-12, "side depth");
    // Above the end cap.
    let top = expect_contact(&placed(sphere(0.5), [0.0, 1.8, 0.0], IDENTITY), &cap);
    assert_vec_close(top.normal, [0.0, 1.0, 0.0], 1e-12, "cap normal");
    assert_close(top.depth, 0.2, 1e-12, "cap depth");
    // Diagonally beyond the cap.
    let diagonal = expect_contact(&placed(sphere(0.5), [0.6, 1.6, 0.0], IDENTITY), &cap);
    let expected = [0.6 / 0.6f64.hypot(0.6), 0.6 / 0.6f64.hypot(0.6), 0.0];
    assert_vec_close(diagonal.normal, expected, 1e-12, "diagonal normal");
    assert_close(
        diagonal.depth,
        1.0 - 0.6f64.hypot(0.6),
        1e-12,
        "diagonal depth",
    );
    assert!(collide(&placed(sphere(0.5), [1.2, 0.3, 0.0], IDENTITY), &cap).is_none());
}

#[test]
fn sphere_box_face_edge_and_inside() {
    let cube = placed(cuboid(1.0, 1.0, 1.0), [0.0; 3], IDENTITY);
    // Over a face.
    let face = expect_contact(&placed(sphere(0.5), [0.0, 1.2, 0.0], IDENTITY), &cube);
    assert_vec_close(face.normal, [0.0, 1.0, 0.0], 1e-12, "face normal");
    assert_close(face.depth, 0.3, 1e-12, "face depth");
    assert_vec_close(face.point_a, [0.0, 0.7, 0.0], 1e-12, "point on the sphere");
    assert_vec_close(face.point_b, [0.0, 1.0, 0.0], 1e-12, "point on the box");

    // Over an edge: the closest point is the edge point (1, 1, 0).
    let edge = expect_contact(&placed(sphere(0.5), [1.2, 1.2, 0.0], IDENTITY), &cube);
    assert_vec_close(
        edge.normal,
        [FRAC_1_SQRT_2, FRAC_1_SQRT_2, 0.0],
        1e-12,
        "edge normal",
    );
    assert_close(edge.depth, 0.5 - 0.2 * SQRT_2, 1e-12, "edge depth");

    // Centre inside the box: leave through the nearest face (+x, gap 0.1).
    let inside = expect_contact(&placed(sphere(0.3), [0.9, 0.2, 0.0], IDENTITY), &cube);
    assert_vec_close(inside.normal, [1.0, 0.0, 0.0], 1e-12, "inside normal");
    assert_close(inside.depth, 0.4, 1e-12, "inside depth");
    assert_vec_close(
        inside.point_b,
        [1.0, 0.2, 0.0],
        1e-12,
        "exit point on the box",
    );

    assert!(collide(&placed(sphere(0.5), [0.0, 1.6, 0.0], IDENTITY), &cube).is_none());
}

#[test]
fn sphere_against_a_rotated_box_vertex() {
    // The unit box turned by 45 degrees about z has a vertex at (sqrt(2), 0).
    let diamond = placed(
        cuboid(1.0, 1.0, 1.0),
        [0.0; 3],
        quat_from_axis_angle([0.0, 0.0, 1.0], FRAC_PI_4),
    );
    let contact = expect_contact(&placed(sphere(0.5), [1.8, 0.0, 0.0], IDENTITY), &diamond);
    assert_vec_close(contact.normal, [1.0, 0.0, 0.0], 1e-9, "normal");
    assert_close(contact.depth, 0.5 - (1.8 - SQRT_2), 1e-9, "depth");
    assert_vec_close(contact.point_b, [SQRT_2, 0.0, 0.0], 1e-9, "vertex");
}

#[test]
fn sphere_cylinder_uses_the_convex_solver() {
    let can = placed(cylinder(1.0, 1.0), [0.0; 3], IDENTITY);
    // Beside the curved wall.
    let side = expect_contact(&placed(sphere(0.5), [1.2, 0.0, 0.0], IDENTITY), &can);
    assert_vec_close(side.normal, [1.0, 0.0, 0.0], 1e-9, "wall normal");
    assert_close(side.depth, 0.3, 1e-9, "wall depth");
    // On top of the flat cap.
    let top = expect_contact(&placed(sphere(0.5), [0.2, 1.3, 0.0], IDENTITY), &can);
    assert_vec_close(top.normal, [0.0, 1.0, 0.0], 1e-9, "cap normal");
    assert_close(top.depth, 0.2, 1e-9, "cap depth");
    // Centre inside the cylinder: leave through the wall (gap 0.2).
    let inside = expect_contact(&placed(sphere(0.3), [0.8, 0.0, 0.0], IDENTITY), &can);
    assert_vec_close(inside.normal, [1.0, 0.0, 0.0], 1e-3, "inside normal");
    assert_close(inside.depth, 0.5, 1e-3, "inside depth");
    assert!(collide(&placed(sphere(0.5), [1.6, 0.0, 0.0], IDENTITY), &can).is_none());
}

// ---------------------------------------------------------------------------
// Capsules, boxes and cylinders against each other
// ---------------------------------------------------------------------------

#[test]
fn parallel_and_crossing_capsules() {
    let upright = placed(capsule(0.5, 1.0), [0.0; 3], IDENTITY);
    // Parallel, side by side.
    let parallel = expect_contact(
        &placed(capsule(0.5, 1.0), [0.8, 0.5, 0.0], IDENTITY),
        &upright,
    );
    assert_vec_close(parallel.normal, [1.0, 0.0, 0.0], 1e-12, "parallel normal");
    assert_close(parallel.depth, 0.2, 1e-12, "parallel depth");

    // Crossing: B lies along x, 0.8 behind A in z.
    let lying = placed(
        capsule(0.5, 1.0),
        [0.0, 0.0, 0.8],
        quat_from_axis_angle([0.0, 0.0, 1.0], FRAC_PI_2),
    );
    let crossing = expect_contact(&upright, &lying);
    assert_vec_close(crossing.normal, [0.0, 0.0, -1.0], 1e-9, "crossing normal");
    assert_close(crossing.depth, 0.2, 1e-9, "crossing depth");

    // Coincident axes still give a unit normal.
    let same = expect_contact(&upright, &placed(capsule(0.5, 1.0), [0.0; 3], IDENTITY));
    assert_close(length(same.normal), 1.0, 1e-12, "degenerate normal length");
    assert_close(same.depth, 1.0, 1e-12, "degenerate depth");
}

#[test]
fn aligned_boxes_overlap_along_the_shortest_axis() {
    let a = placed(cuboid(1.0, 1.0, 1.0), [0.0; 3], IDENTITY);
    let b = placed(cuboid(1.0, 1.0, 1.0), [1.5, 0.2, 0.0], IDENTITY);
    let contact = expect_contact(&a, &b);
    assert_vec_close(
        contact.normal,
        [-1.0, 0.0, 0.0],
        1e-12,
        "normal from B to A",
    );
    assert_close(contact.depth, 0.5, 1e-12, "depth");
    // A's face x = 1 and B's face x = 0.5 interpenetrate; both contact points
    // lie in the overlap of the two faces (y in [-0.8, 1], z in [-1, 1]).
    assert_close(contact.point_a[0], 1.0, 1e-12, "point on A's face");
    assert_close(contact.point_b[0], 0.5, 1e-12, "point on B's face");
    for point in [contact.point_a, contact.point_b] {
        assert!(
            (-0.8..=1.0).contains(&point[1]) && point[2].abs() <= 1.0,
            "{point:?}"
        );
    }
}

#[test]
fn rotated_box_against_a_box_uses_the_tilted_face_axis() {
    // B is the unit box turned by 45 degrees about y: its extent along x is
    // sqrt(2), so the centres 2.2 apart overlap by 1 + sqrt(2) - 2.2 along x.
    let a = placed(cuboid(1.0, 1.0, 1.0), [0.0; 3], IDENTITY);
    let b = placed(
        cuboid(1.0, 1.0, 1.0),
        [2.2, 0.0, 0.0],
        quat_from_axis_angle([0.0, 1.0, 0.0], FRAC_PI_4),
    );
    let contact = expect_contact(&a, &b);
    assert_vec_close(contact.normal, [-1.0, 0.0, 0.0], 1e-9, "normal");
    assert_close(contact.depth, 1.0 + SQRT_2 - 2.2, 1e-9, "depth");
    assert!(
        collide(
            &a,
            &placed(
                cuboid(1.0, 1.0, 1.0),
                [2.5, 0.0, 0.0],
                quat_from_axis_angle([0.0, 1.0, 0.0], FRAC_PI_4)
            )
        )
        .is_none()
    );
}

#[test]
fn box_resting_on_the_edge_of_another_box() {
    // A tilted box balanced with an edge over the top face of a flat box.
    let floor = placed(cuboid(2.0, 0.5, 2.0), [0.0; 3], IDENTITY);
    let tilted = placed(
        cuboid(0.5, 0.5, 0.5),
        [0.0, 0.5 + SQRT_2 * 0.5 - 0.1, 0.0],
        quat_from_axis_angle([0.0, 0.0, 1.0], FRAC_PI_4),
    );
    let contact = expect_contact(&tilted, &floor);
    assert_vec_close(contact.normal, [0.0, 1.0, 0.0], 1e-9, "normal");
    assert_close(contact.depth, 0.1, 1e-9, "depth");
}

#[test]
fn capsule_against_a_box_uses_the_convex_solver() {
    let cube = placed(cuboid(1.0, 1.0, 1.0), [0.0; 3], IDENTITY);
    let upright = placed(capsule(0.5, 1.0), [1.3, 0.0, 0.0], IDENTITY);
    let contact = expect_contact(&upright, &cube);
    // Capsule (rank lower) is A; the box face x = 1 is 0.3 from the axis.
    assert_vec_close(contact.normal, [1.0, 0.0, 0.0], 1e-9, "normal");
    assert_close(contact.depth, 0.2, 1e-9, "depth");
    let mirrored = expect_contact(&cube, &upright);
    assert_vec_close(mirrored.normal, [-1.0, 0.0, 0.0], 1e-9, "mirrored normal");
    assert_close(mirrored.depth, 0.2, 1e-9, "mirrored depth");
}

#[test]
fn overlapping_cylinders_separate_along_the_shortest_direction() {
    let a = placed(cylinder(1.0, 1.0), [0.0; 3], IDENTITY);
    let b = placed(cylinder(1.0, 1.0), [1.5, 0.0, 0.0], IDENTITY);
    let contact = expect_contact(&a, &b);
    assert_vec_close(contact.normal, [-1.0, 0.0, 0.0], 1e-3, "normal");
    assert_close(contact.depth, 0.5, 1e-3, "depth");
    // Stacked cylinders overlap along the axis.
    let above = placed(cylinder(1.0, 1.0), [0.1, 1.7, 0.0], IDENTITY);
    let stacked = expect_contact(&above, &a);
    assert_vec_close(stacked.normal, [0.0, 1.0, 0.0], 1e-3, "stacked normal");
    assert_close(stacked.depth, 0.3, 1e-3, "stacked depth");
}

#[test]
fn box_cylinder_and_capsule_cylinder_contacts() {
    let can = placed(cylinder(1.0, 1.0), [0.0; 3], IDENTITY);
    let cube = placed(cuboid(0.5, 0.5, 0.5), [1.4, 0.0, 0.0], IDENTITY);
    // Box (rank 2) is A, cylinder (rank 3) is B.
    let contact = expect_contact(&cube, &can);
    assert_vec_close(contact.normal, [1.0, 0.0, 0.0], 1e-6, "box-cylinder normal");
    assert_close(contact.depth, 0.1, 1e-6, "box-cylinder depth");

    let pill = placed(capsule(0.25, 0.5), [0.0, 1.6, 0.0], IDENTITY);
    let contact = expect_contact(&pill, &can);
    assert_vec_close(
        contact.normal,
        [0.0, 1.0, 0.0],
        1e-6,
        "capsule-cylinder normal",
    );
    // The cylinder top is at y = 1, the capsule reaches down to 1.6 - 0.5 - 0.25.
    assert_close(
        contact.depth,
        1.0 - (1.6 - 0.5 - 0.25),
        1e-6,
        "capsule-cylinder depth",
    );
}

// ---------------------------------------------------------------------------
// Planes
// ---------------------------------------------------------------------------

#[test]
fn sphere_plane_normal_is_the_plane_normal() {
    let ground = placed(plane([0.0, 1.0, 0.0], 0.0), [0.0; 3], IDENTITY);
    let contact = expect_contact(&placed(sphere(0.5), [3.0, 0.3, 1.0], IDENTITY), &ground);
    assert_vec_close(contact.normal, [0.0, 1.0, 0.0], 1e-12, "normal");
    assert_close(contact.depth, 0.2, 1e-12, "depth");
    assert_vec_close(
        contact.point_a,
        [3.0, -0.2, 1.0],
        1e-12,
        "point on the sphere",
    );
    assert_vec_close(
        contact.point_b,
        [3.0, 0.0, 1.0],
        1e-12,
        "point on the plane",
    );
    // The plane is the second shape however the arguments are ordered.
    let mirrored = expect_contact(&ground, &placed(sphere(0.5), [3.0, 0.3, 1.0], IDENTITY));
    assert_vec_close(mirrored.normal, [0.0, -1.0, 0.0], 1e-12, "mirrored normal");
    assert!(collide(&placed(sphere(0.5), [0.0, 0.6, 0.0], IDENTITY), &ground).is_none());
}

#[test]
fn plane_offset_and_tilt_are_honoured() {
    // Plane y = 2.
    let raised = placed(plane([0.0, 1.0, 0.0], 2.0), [0.0; 3], IDENTITY);
    let contact = expect_contact(&placed(sphere(0.5), [0.0, 2.3, 0.0], IDENTITY), &raised);
    assert_close(contact.depth, 0.2, 1e-12, "raised depth");

    // Tilted plane through the origin with normal (1, 1, 0) / sqrt(2).
    let n = [FRAC_1_SQRT_2, FRAC_1_SQRT_2, 0.0];
    let tilted = placed(plane(n, 0.0), [0.0; 3], IDENTITY);
    let contact = expect_contact(&placed(sphere(0.5), [0.2, 0.2, 0.0], IDENTITY), &tilted);
    assert_vec_close(contact.normal, n, 1e-12, "tilted normal");
    assert_close(
        contact.depth,
        0.5 - 0.4 * FRAC_1_SQRT_2,
        1e-12,
        "tilted depth",
    );

    // A plane attached to a translated body moves with it.
    let moved = placed(plane([0.0, 1.0, 0.0], 0.0), [0.0, 3.0, 0.0], IDENTITY);
    let contact = expect_contact(&placed(sphere(0.5), [0.0, 3.2, 0.0], IDENTITY), &moved);
    assert_close(contact.depth, 0.3, 1e-12, "moved depth");
}

#[test]
fn box_plane_flat_and_tilted() {
    let ground = placed(plane([0.0, 1.0, 0.0], 0.0), [0.0; 3], IDENTITY);
    let flat = expect_contact(
        &placed(cuboid(0.5, 0.5, 0.5), [0.0, 0.4, 0.0], IDENTITY),
        &ground,
    );
    assert_vec_close(flat.normal, [0.0, 1.0, 0.0], 1e-12, "flat normal");
    assert_close(flat.depth, 0.1, 1e-12, "flat depth");
    assert_vec_close(flat.point_a, [0.0, -0.1, 0.0], 1e-12, "flat contact point");

    // Turned by 45 degrees about z the lowest edge is at y - sqrt(2) / 2.
    let tilted = expect_contact(
        &placed(
            cuboid(0.5, 0.5, 0.5),
            [0.0, 0.5, 0.0],
            quat_from_axis_angle([0.0, 0.0, 1.0], FRAC_PI_4),
        ),
        &ground,
    );
    assert_close(tilted.depth, SQRT_2 / 2.0 - 0.5, 1e-12, "tilted depth");
    assert_vec_close(
        tilted.point_a,
        [0.0, 0.5 - SQRT_2 / 2.0, 0.0],
        1e-12,
        "edge midpoint",
    );

    // Rolled about x: the lowest edge runs along x and lies at y - sqrt(2) / 2.
    let rolled = placed(
        cuboid(0.5, 0.5, 0.5),
        [0.0, 0.6, 0.0],
        quat_from_axis_angle([1.0, 0.0, 0.0], FRAC_PI_4),
    );
    let contact = expect_contact(&rolled, &ground);
    assert_close(contact.depth, SQRT_2 / 2.0 - 0.6, 1e-12, "rolled depth");
    assert_vec_close(contact.normal, [0.0, 1.0, 0.0], 1e-12, "normal");
    assert_close(contact.point_a[1], 0.6 - SQRT_2 / 2.0, 1e-12, "edge height");
}

#[test]
fn capsule_plane_lying_and_standing() {
    let ground = placed(plane([0.0, 1.0, 0.0], 0.0), [0.0; 3], IDENTITY);
    let lying = expect_contact(
        &placed(
            capsule(0.5, 1.0),
            [0.0, 0.3, 0.0],
            quat_from_axis_angle([0.0, 0.0, 1.0], FRAC_PI_2),
        ),
        &ground,
    );
    assert_close(lying.depth, 0.2, 1e-12, "lying depth");
    assert_vec_close(lying.point_a, [0.0, -0.2, 0.0], 1e-9, "centre of pressure");

    let standing = expect_contact(
        &placed(capsule(0.5, 1.0), [0.0, 1.2, 0.0], IDENTITY),
        &ground,
    );
    assert_close(standing.depth, 0.3, 1e-12, "standing depth");
    assert_vec_close(standing.point_a, [0.0, -0.3, 0.0], 1e-12, "lowest point");
}

#[test]
fn cylinder_plane_standing_lying_and_tilted() {
    let ground = placed(plane([0.0, 1.0, 0.0], 0.0), [0.0; 3], IDENTITY);
    let standing = expect_contact(
        &placed(cylinder(0.5, 1.0), [0.0, 0.9, 0.0], IDENTITY),
        &ground,
    );
    assert_close(standing.depth, 0.1, 1e-12, "standing depth");
    assert_vec_close(standing.point_a, [0.0, -0.1, 0.0], 1e-12, "cap centre");

    let lying = expect_contact(
        &placed(
            cylinder(0.5, 1.0),
            [0.0, 0.4, 0.0],
            quat_from_axis_angle([0.0, 0.0, 1.0], FRAC_PI_2),
        ),
        &ground,
    );
    assert_close(lying.depth, 0.1, 1e-9, "lying depth");

    let tilted = expect_contact(
        &placed(
            cylinder(0.5, 1.0),
            [0.0, 1.0, 0.0],
            quat_from_axis_angle([0.0, 0.0, 1.0], FRAC_PI_4),
        ),
        &ground,
    );
    // Lowest point: centre - axis * hh * cos45 - radial * r * sin45 ...
    let expected_low = 1.0 - (1.0 * FRAC_1_SQRT_2 + 0.5 * FRAC_1_SQRT_2);
    assert_close(tilted.depth, -expected_low, 1e-9, "tilted depth");
}

#[test]
fn plane_plane_pairs_never_collide() {
    let a = placed(plane([0.0, 1.0, 0.0], 0.0), [0.0; 3], IDENTITY);
    let b = placed(plane([1.0, 0.0, 0.0], 0.0), [0.0; 3], IDENTITY);
    assert!(collide(&a, &b).is_none());
}

// ---------------------------------------------------------------------------
// Randomised invariants over every shape pair
// ---------------------------------------------------------------------------

/// Small deterministic pseudo random generator (64-bit LCG).
struct Lcg(u64);

impl Lcg {
    fn next_unit(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }

    fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.next_unit()
    }

    fn vector(&mut self, extent: f64) -> Vec3 {
        [
            self.range(-extent, extent),
            self.range(-extent, extent),
            self.range(-extent, extent),
        ]
    }

    fn rotation(&mut self) -> [f64; 4] {
        loop {
            let q = [
                self.range(-1.0, 1.0),
                self.range(-1.0, 1.0),
                self.range(-1.0, 1.0),
                self.range(-1.0, 1.0),
            ];
            let n = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]).sqrt();
            if n > 0.2 && n <= 1.0 {
                return [q[0] / n, q[1] / n, q[2] / n, q[3] / n];
            }
        }
    }

    fn shape(&mut self, kind: usize) -> Shape {
        match kind {
            0 => sphere(self.range(0.2, 1.0)),
            1 => capsule(self.range(0.15, 0.6), self.range(0.2, 1.0)),
            2 => cuboid(
                self.range(0.2, 1.0),
                self.range(0.2, 1.0),
                self.range(0.2, 1.0),
            ),
            _ => cylinder(self.range(0.2, 0.8), self.range(0.2, 1.0)),
        }
    }
}

const KIND_NAMES: [&str; 4] = ["sphere", "capsule", "box", "cylinder"];

/// Translate a placed shape by `offset`.
fn translated(shape: &WorldShape, offset: Vec3) -> WorldShape {
    let mut moved = *shape;
    moved.position = add(shape.position, offset);
    moved
}

/// Invariants every contact must satisfy.
///
/// * `normal` is a unit vector and `depth` is positive;
/// * `normal` points from B toward A: moving A along `+normal` by the depth
///   (plus a small allowance) separates the pair;
/// * the depth is not grossly over-estimated: moving A by 80 % of it leaves
///   the pair in contact;
/// * exchanging A and B mirrors the contact.
fn check_contact_invariants(
    a: &WorldShape,
    b: &WorldShape,
    contact: &Contact,
    tolerance: f64,
    label: &str,
) {
    assert_close(
        length(contact.normal),
        1.0,
        1e-9,
        &format!("{label}: normal length"),
    );
    assert!(contact.depth > 0.0, "{label}: depth must be positive");
    let margin = contact.depth * tolerance + 1e-7;

    let separated = translated(a, scale(contact.normal, contact.depth + margin));
    if let Some(remaining) = collide(&separated, b) {
        panic!(
            "{label}: moving A along the normal by the depth {} did not separate the pair \
             (remaining depth {}, normal {:?})\n a = {a:?}\n b = {b:?}",
            contact.depth, remaining.depth, contact.normal
        );
    }
    let partly = translated(a, scale(contact.normal, 0.8 * contact.depth));
    assert!(
        collide(&partly, b).is_some(),
        "{label}: depth {} over-estimated, 80 % of it already separates the pair\n a = {a:?}\n b = {b:?}",
        contact.depth
    );

    let mirrored =
        collide(b, a).unwrap_or_else(|| panic!("{label}: the mirrored pair has no contact"));
    assert!(
        (mirrored.depth - contact.depth).abs() <= tolerance * contact.depth + 1e-9,
        "{label}: mirrored depth {} vs {}",
        mirrored.depth,
        contact.depth
    );
    assert!(
        dot(mirrored.normal, contact.normal) < -1.0 + 1e-3 + tolerance,
        "{label}: mirrored normal {:?} vs {:?}",
        mirrored.normal,
        contact.normal
    );
}

#[test]
fn contact_invariants_hold_for_every_solid_pair() {
    let mut rng = Lcg(0x5EED_1234_ABCD_0001);
    for (kind_a, name_a) in KIND_NAMES.iter().enumerate() {
        for (kind_b, name_b) in KIND_NAMES.iter().enumerate() {
            let analytic = !(kind_a == 3
                || kind_b == 3
                || (kind_a == 1 && kind_b == 2)
                || (kind_a == 2 && kind_b == 1));
            let tolerance = if analytic { 1e-9 } else { 2e-3 };
            let mut found = 0;
            for trial in 0..400 {
                let a = placed(rng.shape(kind_a), rng.vector(1.0), rng.rotation());
                let b = placed(rng.shape(kind_b), [0.0; 3], rng.rotation());
                let Some(contact) = collide(&a, &b) else {
                    continue;
                };
                found += 1;
                let label = format!("{name_a} vs {name_b} (trial {trial})");
                check_contact_invariants(&a, &b, &contact, tolerance, &label);
            }
            assert!(
                found >= 40,
                "{name_a} vs {name_b}: only {found} of 400 random poses touched"
            );
        }
    }
}

/// `true` when `p` lies inside `shape` shrunk by `margin` (a negative margin
/// grows the shape).
fn contains(shape: &WorldShape, p: Vec3, margin: f64) -> bool {
    let local = shape.to_local_point(p);
    match shape.shape {
        Shape::Sphere { radius } => length(local) <= radius - margin,
        Shape::Capsule {
            radius,
            half_height,
        } => {
            let axis_point = [0.0, local[1].clamp(-half_height, half_height), 0.0];
            length(sub(local, axis_point)) <= radius - margin
        }
        Shape::Box { half_extents } => (0..3).all(|i| local[i].abs() <= half_extents[i] - margin),
        Shape::Cylinder {
            radius,
            half_height,
        } => local[0].hypot(local[2]) <= radius - margin && local[1].abs() <= half_height - margin,
        Shape::Plane { .. } => false,
    }
}

#[test]
fn overlapping_solids_always_produce_a_contact() {
    // Independent oracle: a sampled point that lies inside both shapes proves
    // that they overlap, so the narrow phase must report a contact.
    let mut rng = Lcg(0x0BAD_5EED_7777_0004);
    for (kind_a, name_a) in KIND_NAMES.iter().enumerate() {
        for (kind_b, name_b) in KIND_NAMES.iter().enumerate() {
            let mut witnessed = 0;
            for trial in 0..300 {
                let a = placed(rng.shape(kind_a), rng.vector(1.2), rng.rotation());
                let b = placed(rng.shape(kind_b), [0.0; 3], rng.rotation());
                let witness = (0..400).any(|_| {
                    let p = add(a.position, rng.vector(a.bound_radius));
                    contains(&a, p, 1e-3) && contains(&b, p, 1e-3)
                });
                if witness {
                    witnessed += 1;
                    assert!(
                        collide(&a, &b).is_some(),
                        "{name_a} vs {name_b} (trial {trial}): the shapes share a point but no contact was found\n a = {a:?}\n b = {b:?}"
                    );
                }
            }
            assert!(
                witnessed >= 20,
                "{name_a} vs {name_b}: only {witnessed} overlapping poses were generated"
            );
        }
    }
}

#[test]
fn contact_invariants_hold_for_every_shape_against_a_plane() {
    let mut rng = Lcg(0x5EED_9876_FEDC_0002);
    for (kind, name) in KIND_NAMES.iter().enumerate() {
        let mut found = 0;
        for trial in 0..400 {
            let normal = loop {
                let v = rng.vector(1.0);
                let len = length(v);
                if len > 0.2 {
                    break scale(v, 1.0 / len);
                }
            };
            let ground = placed(plane(normal, rng.range(-0.3, 0.3)), [0.0; 3], IDENTITY);
            let a = placed(rng.shape(kind), rng.vector(1.0), rng.rotation());
            let Some(contact) = collide(&a, &ground) else {
                continue;
            };
            found += 1;
            let label = format!("{name} vs plane (trial {trial})");
            assert_vec_close(contact.normal, normal, 1e-12, &format!("{label}: normal"));
            // Plane contacts are exact: A moves out by exactly the depth.
            let lifted = translated(&a, scale(normal, contact.depth + 1e-9));
            assert!(
                collide(&lifted, &ground).is_none(),
                "{label}: still touching"
            );
            let partly = translated(&a, scale(normal, 0.99 * contact.depth));
            assert!(
                collide(&partly, &ground).is_some(),
                "{label}: depth too large"
            );
            // The contact point lies on A's side of the plane, the matching
            // point on B lies on the plane.
            let gap = dot(sub(contact.point_a, contact.point_b), normal);
            assert!(gap <= 1e-9, "{label}: point A is above point B by {gap}");
        }
        assert!(found >= 100, "{name}: only {found} contacts");
    }
}

#[test]
fn single_point_contacts_satisfy_the_documented_point_relation() {
    // For single-point contacts `(point_a - point_b) . normal == -depth`.
    let a = placed(sphere(0.7), [0.3, 0.2, -0.1], IDENTITY);
    let b = placed(sphere(0.9), [1.0, 0.1, 0.2], IDENTITY);
    let contact = expect_contact(&a, &b);
    assert_close(
        dot(sub(contact.point_a, contact.point_b), contact.normal),
        -contact.depth,
        1e-12,
        "sphere-sphere",
    );
    // `quat_from_axis_angle` expects a unit axis.
    let axis = [0.3, 0.5, 0.2];
    let axis = scale(axis, 1.0 / length(axis));
    let cube = placed(
        cuboid(1.0, 1.0, 1.0),
        [0.0; 3],
        quat_from_axis_angle(axis, 0.7),
    );
    let ball = placed(sphere(0.6), [1.1, 0.4, 0.2], IDENTITY);
    let contact = expect_contact(&ball, &cube);
    assert_close(
        dot(sub(contact.point_a, contact.point_b), contact.normal),
        -contact.depth,
        1e-12,
        "sphere-box",
    );
    // A's deepest point lies `depth` inside B, i.e. against the normal.
    let shifted = madd(contact.point_b, contact.normal, -contact.depth);
    assert_vec_close(shifted, contact.point_a, 1e-12, "point_b - n * depth");
}

// ---------------------------------------------------------------------------
// Whole-engine scenarios
// ---------------------------------------------------------------------------

fn space_engine() -> WasmPhysicsEngine {
    WasmPhysicsEngine::new(0.0, 0.0, 0.0)
}

/// Add a dynamic sphere without damping, so that the velocities a test sets
/// are exactly the velocities the contact solver sees.
fn add_sphere_body(engine: &mut WasmPhysicsEngine, mass: f64, position: Vec3, radius: f64) -> u32 {
    let body = engine
        .add_dynamic_body(mass, position[0], position[1], position[2])
        .expect("dynamic body");
    engine
        .set_body_linear_damping(body, 0.0)
        .expect("linear damping");
    engine
        .set_body_angular_damping(body, 0.0)
        .expect("angular damping");
    engine
        .add_sphere_collider(body, radius)
        .expect("sphere collider");
    body
}

#[test]
fn approaching_spheres_receive_an_impulse_and_separate() {
    // Unit spheres at x = 0 and x = 1.5 (0.5 overlap) moving towards each
    // other at 2 m/s each.
    let mut engine = space_engine();
    let a = add_sphere_body(&mut engine, 1.0, [0.0; 3], 1.0);
    let b = add_sphere_body(&mut engine, 1.0, [1.5, 0.0, 0.0], 1.0);
    engine.set_velocity(a, 2.0, 0.0, 0.0).expect("velocity a");
    engine.set_velocity(b, -2.0, 0.0, 0.0).expect("velocity b");
    engine.step(DT).expect("step");

    let va = engine.get_velocity(a);
    let vb = engine.get_velocity(b);
    assert!(va[0] <= 0.0, "A must be pushed back, vx_a = {}", va[0]);
    assert!(vb[0] >= 0.0, "B must be pushed back, vx_b = {}", vb[0]);
    // Restitution 0.3: the closing speed of 4 m/s becomes a 1.2 m/s recession.
    assert_close(va[0], -0.6, 1e-12, "vx_a");
    assert_close(vb[0], 0.6, 1e-12, "vx_b");

    let contacts = engine.get_contacts();
    assert_eq!(contacts.len(), 1);
    let contact = &contacts[0];
    assert_eq!((contact.body_a, contact.body_b), (a, b));
    assert_vec_close(contact.normal, [-1.0, 0.0, 0.0], 1e-12, "reported normal");
    assert!(
        contact.relative_velocity < 0.0,
        "the bodies were approaching"
    );
    assert!(contact.impulse > 0.0);
    assert!(contact.depth > 0.0);
}

#[test]
fn spheres_at_rest_with_overlap_are_pushed_apart() {
    let mut engine = space_engine();
    let a = add_sphere_body(&mut engine, 1.0, [0.0; 3], 1.0);
    let b = add_sphere_body(&mut engine, 1.0, [1.5, 0.0, 0.0], 1.0);
    let slop = engine.config().contact_slop;

    let mut previous = 1.5;
    for _ in 0..60 {
        engine.step(DT).expect("step");
        let distance = engine.get_position(b)[0] - engine.get_position(a)[0];
        assert!(
            distance >= previous - 1e-12,
            "the bodies must never move closer"
        );
        previous = distance;
    }
    // Each substep removes 20 % of the penetration beyond the slop, so after
    // 60 steps only (0.5 - slop) * 0.8^60 < 1e-6 of it is left.
    assert!(
        previous >= 2.0 - slop - 1e-5,
        "after 60 steps the overlap is down to the slop, distance = {previous}"
    );
    assert!(engine.get_position(a)[0] < 0.0 && engine.get_position(b)[0] > 1.5);
    // Positional correction does not add speed.
    assert_eq!(engine.get_velocity(a), [0.0; 3]);
    assert_eq!(engine.get_velocity(b), [0.0; 3]);
}

#[test]
fn collision_conserves_momentum_and_the_light_body_gives_way() {
    let mut engine = space_engine();
    let heavy = add_sphere_body(&mut engine, 3.0, [0.0; 3], 1.0);
    let light = add_sphere_body(&mut engine, 1.0, [1.8, 0.0, 0.0], 1.0);
    engine.set_velocity(heavy, 1.0, 0.0, 0.0).expect("velocity");
    engine
        .set_velocity(light, -1.0, 0.0, 0.0)
        .expect("velocity");
    let momentum_before = 3.0 * 1.0 - 1.0;
    engine.step(DT).expect("step");
    let momentum_after = 3.0 * engine.get_velocity(heavy)[0] + engine.get_velocity(light)[0];
    assert_close(momentum_after, momentum_before, 1e-12, "momentum");
    // Closing speed 2 m/s becomes a 0.6 m/s recession (e = 0.3).
    let closing = engine.get_velocity(heavy)[0] - engine.get_velocity(light)[0];
    assert_close(closing, -0.6, 1e-12, "post-collision closing speed");
    assert!(
        engine.get_velocity(light)[0] > 0.0,
        "the light body is thrown back"
    );
}

#[test]
fn spheres_moving_apart_are_not_slowed_down() {
    let mut engine = space_engine();
    let a = add_sphere_body(&mut engine, 1.0, [0.0; 3], 1.0);
    let b = add_sphere_body(&mut engine, 1.0, [1.9, 0.0, 0.0], 1.0);
    engine.set_velocity(a, -1.0, 0.0, 0.0).expect("velocity");
    engine.set_velocity(b, 1.0, 0.0, 0.0).expect("velocity");
    engine.step(DT).expect("step");
    // They overlap but are already separating: no impulse, only position
    // correction.
    assert_eq!(engine.get_velocity(a), [-1.0, 0.0, 0.0]);
    assert_eq!(engine.get_velocity(b), [1.0, 0.0, 0.0]);
    assert!(engine.get_contacts()[0].impulse == 0.0);
}

#[test]
fn static_bodies_are_not_moved_by_contacts() {
    let mut engine = space_engine();
    let wall = engine.add_static_body(1.5, 0.0, 0.0).expect("wall");
    engine
        .add_sphere_collider(wall, 1.0)
        .expect("wall collider");
    let ball = add_sphere_body(&mut engine, 1.0, [0.0; 3], 1.0);
    engine.set_velocity(ball, 2.0, 0.0, 0.0).expect("velocity");
    engine.step(DT).expect("step");
    assert_eq!(engine.get_position(wall), [1.5, 0.0, 0.0]);
    assert_eq!(engine.get_velocity(wall), [0.0; 3]);
    assert!(
        engine.get_velocity(ball)[0] < 0.0,
        "the ball bounces off the wall"
    );
}

#[test]
fn contact_between_a_compound_collider_and_a_sphere_uses_the_local_offset() {
    let mut engine = space_engine();
    let body = engine
        .add_dynamic_body(1.0, 0.0, 0.0, 0.0)
        .expect("compound body");
    for x in [-3.0, 3.0] {
        engine
            .add_collider(
                body,
                &ColliderConfig {
                    local_position: [x, 0.0, 0.0],
                    ..ColliderConfig::sphere(0.5)
                },
            )
            .expect("local sphere");
    }
    let ball = add_sphere_body(&mut engine, 1.0, [3.8, 0.0, 0.0], 0.5);
    engine.step(DT).expect("step");
    let contacts = engine.get_contacts();
    assert_eq!(contacts.len(), 1, "only the right-hand sphere is touched");
    assert_eq!(contacts[0].body_a, body);
    assert_eq!(contacts[0].body_b, ball);
    assert_vec_close(contacts[0].normal, [-1.0, 0.0, 0.0], 1e-9, "normal");
    assert_close(contacts[0].depth, 0.2, 1e-9, "depth");
}

#[test]
fn colliders_of_one_body_do_not_collide_with_each_other() {
    let mut engine = space_engine();
    let body = engine.add_dynamic_body(1.0, 0.0, 0.0, 0.0).expect("body");
    engine.add_sphere_collider(body, 1.0).expect("collider a");
    engine.add_sphere_collider(body, 1.0).expect("collider b");
    engine.step(DT).expect("step");
    assert_eq!(engine.get_contact_count(), 0);
}

#[test]
fn contacts_of_a_removed_body_disappear() {
    let mut engine = space_engine();
    let a = add_sphere_body(&mut engine, 1.0, [0.0; 3], 1.0);
    let b = add_sphere_body(&mut engine, 1.0, [1.5, 0.0, 0.0], 1.0);
    engine.step(DT).expect("step");
    assert_eq!(engine.get_contact_count(), 1);
    engine.remove_body(b).expect("remove");
    engine.step(DT).expect("step");
    assert_eq!(engine.get_contact_count(), 0);
    assert!(engine.contains_body(a));
}

/// Drop `collider` (on a dynamic body of 1 kg) from `height` onto a ground
/// plane through the origin and return the body after 3 s.
fn settle_on_ground(
    collider: ColliderConfig,
    rotation: [f64; 4],
    height: f64,
) -> (WasmPhysicsEngine, u32) {
    let mut engine = WasmPhysicsEngine::new(0.0, -9.81, 0.0);
    let ground = engine.add_static_body(0.0, 0.0, 0.0).expect("ground");
    engine
        .add_plane_collider(ground, 0.0, 1.0, 0.0, 0.0)
        .expect("ground plane");
    let body = engine
        .add_rigid_body(&RigidBodyConfig {
            rotation,
            ..RigidBodyConfig::new().with_position(0.0, height, 0.0)
        })
        .expect("falling body");
    engine.add_collider(body, &collider).expect("collider");
    for _ in 0..180 {
        engine.step(DT).expect("step");
    }
    (engine, body)
}

/// A body resting under gravity sinks slightly into its support: gravity adds
/// about 2 mm of penetration per substep, the positional correction removes
/// 20 % of what exceeds the slop, so the penetration settles around 1.5 cm.
const MAX_RESTING_SINK: f64 = 0.02;

fn assert_resting_at(engine: &WasmPhysicsEngine, body: u32, height: f64) {
    let position = engine.get_position(body);
    assert!(
        position[1] <= height + 1e-9 && position[1] >= height - MAX_RESTING_SINK,
        "expected the body to rest at y = {height}, got {}",
        position[1]
    );
    assert!(
        engine.get_velocity(body)[1].abs() < 0.2,
        "vy = {}",
        engine.get_velocity(body)[1]
    );
}

#[test]
fn box_dropped_on_a_plane_rests_at_its_half_extent() {
    let (engine, body) = settle_on_ground(ColliderConfig::cuboid(0.5, 0.5, 0.5), IDENTITY, 2.0);
    assert_resting_at(&engine, body, 0.5);
    assert!(engine.get_contact_count() >= 1);
    let contact = &engine.get_contacts()[0];
    assert_eq!(contact.body_a, body);
    assert_vec_close(contact.normal, [0.0, 1.0, 0.0], 1e-12, "normal");
    // Contacts exert no torque, so the box does not tumble.
    let q = engine.get_rotation(body);
    assert_vec_close([q[0], q[1], q[2]], [0.0; 3], 1e-12, "rotation axis part");
    assert_close(q[3].abs(), 1.0, 1e-12, "rotation w");
}

#[test]
fn every_solid_comes_to_rest_on_the_ground_at_its_support_height() {
    let cases: [(&str, ColliderConfig, f64); 4] = [
        ("sphere", ColliderConfig::sphere(0.6), 0.6),
        ("box", ColliderConfig::cuboid(0.4, 0.3, 0.5), 0.3),
        ("capsule", ColliderConfig::capsule(0.25, 1.0), 0.75),
        ("cylinder", ColliderConfig::cylinder(0.5, 1.2), 0.6),
    ];
    for (name, collider, rest_height) in cases {
        let (engine, body) = settle_on_ground(collider, IDENTITY, 2.5);
        let position = engine.get_position(body);
        assert!(
            position[1] <= rest_height + 1e-9 && position[1] >= rest_height - MAX_RESTING_SINK,
            "{name}: expected to rest at y = {rest_height}, got {}",
            position[1]
        );
    }
}

#[test]
fn tilted_box_dropped_on_a_plane_stops_at_the_lowest_corner() {
    // Turned 45 degrees about z, the lowest edge is sqrt(2) / 2 below the
    // centre; contacts exert no torque, so the box keeps its attitude.
    let rotation = quat_from_axis_angle([0.0, 0.0, 1.0], FRAC_PI_4);
    let (engine, body) = settle_on_ground(ColliderConfig::cuboid(0.5, 0.5, 0.5), rotation, 3.0);
    assert_resting_at(&engine, body, SQRT_2 / 2.0);
}

#[test]
fn sphere_dropped_on_a_box_rests_on_top() {
    let mut engine = WasmPhysicsEngine::new(0.0, -9.81, 0.0);
    let base = engine.add_static_body(0.0, 0.0, 0.0).expect("base");
    engine
        .add_box_collider(base, 2.0, 1.0, 2.0)
        .expect("base collider");
    let ball = add_sphere_body(&mut engine, 1.0, [0.5, 3.0, 0.2], 0.5);
    for _ in 0..180 {
        engine.step(DT).expect("step");
    }
    assert_resting_at(&engine, ball, 1.5);
    assert_close(engine.get_position(ball)[0], 0.5, 1e-9, "x");
}

#[test]
fn box_dropped_on_a_box_rests_on_top() {
    let mut engine = WasmPhysicsEngine::new(0.0, -9.81, 0.0);
    let base = engine.add_static_body(0.0, 0.0, 0.0).expect("base");
    engine
        .add_box_collider(base, 2.0, 1.0, 2.0)
        .expect("base collider");
    let crate_body = engine.add_dynamic_body(2.0, 0.3, 3.0, -0.2).expect("crate");
    engine
        .add_box_collider(crate_body, 0.5, 0.4, 0.6)
        .expect("crate collider");
    for _ in 0..180 {
        engine.step(DT).expect("step");
    }
    assert_resting_at(&engine, crate_body, 1.4);
}

#[test]
fn capsule_and_sphere_dropped_on_a_cylinder_rest_on_top() {
    for collider in [
        ColliderConfig::capsule(0.25, 1.0),
        ColliderConfig::sphere(0.5),
    ] {
        let mut engine = WasmPhysicsEngine::new(0.0, -9.81, 0.0);
        let base = engine.add_static_body(0.0, 0.0, 0.0).expect("base");
        engine
            .add_cylinder_collider(base, 1.5, 2.0)
            .expect("base collider");
        let body = engine.add_dynamic_body(1.0, 0.3, 3.0, 0.0).expect("body");
        let rest = match collider.shape_type {
            crate::types::ColliderShapeType::Capsule => 1.0 + 0.5 + 0.25,
            _ => 1.0 + 0.5,
        };
        engine.add_collider(body, &collider).expect("collider");
        for _ in 0..240 {
            engine.step(DT).expect("step");
        }
        assert_resting_at(&engine, body, rest);
    }
}

#[test]
fn bodies_do_not_tunnel_through_a_thin_box_at_moderate_speed() {
    let mut engine = WasmPhysicsEngine::new(0.0, 0.0, 0.0);
    let wall = engine.add_static_body(0.0, 0.0, 0.0).expect("wall");
    engine
        .add_box_collider(wall, 0.5, 5.0, 5.0)
        .expect("wall collider");
    let ball = add_sphere_body(&mut engine, 1.0, [-3.0, 0.0, 0.0], 0.25);
    engine.set_velocity(ball, 15.0, 0.0, 0.0).expect("velocity");
    for _ in 0..120 {
        engine.step(DT).expect("step");
    }
    assert!(engine.get_position(ball)[0] < 0.0, "the ball bounced back");
    assert!(engine.get_velocity(ball)[0] < 0.0);
}

#[test]
fn engine_contact_normals_always_point_from_b_toward_a() {
    // A small pile: every reported contact must push body A away from B.
    let mut engine = WasmPhysicsEngine::new(0.0, -9.81, 0.0);
    let ground = engine.add_static_body(0.0, 0.0, 0.0).expect("ground");
    engine
        .add_plane_collider(ground, 0.0, 1.0, 0.0, 0.0)
        .expect("ground plane");
    let mut rng = Lcg(0x00C0_FFEE_0000_0003);
    for _ in 0..12 {
        let position = [
            rng.range(-1.5, 1.5),
            rng.range(0.4, 4.0),
            rng.range(-1.5, 1.5),
        ];
        let body = engine
            .add_dynamic_body(1.0, position[0], position[1], position[2])
            .expect("body");
        match (rng.next_unit() * 3.0) as u32 {
            0 => engine.add_sphere_collider(body, 0.4).expect("sphere"),
            1 => engine.add_box_collider(body, 0.3, 0.3, 0.3).expect("box"),
            _ => engine
                .add_capsule_collider(body, 0.2, 0.5)
                .expect("capsule"),
        };
    }
    let mut checked = 0;
    for _ in 0..120 {
        engine.step(DT).expect("step");
        for contact in engine.get_contacts() {
            checked += 1;
            let separation = sub(
                engine.get_position(contact.body_a),
                engine.get_position(contact.body_b),
            );
            if contact.body_b == ground {
                assert_vec_close(contact.normal, [0.0, 1.0, 0.0], 1e-12, "ground normal");
            } else {
                assert!(
                    dot(separation, contact.normal) > 0.0,
                    "contact normal {:?} does not point from B to A (offset {:?})",
                    contact.normal,
                    separation
                );
            }
        }
    }
    assert!(checked > 20, "the pile produced only {checked} contacts");
}

// ---------------------------------------------------------------------------
// Plane colliders live in the frame of their body
// ---------------------------------------------------------------------------

#[test]
fn plane_colliders_are_expressed_in_the_frame_of_their_body() {
    // The plane y = 0 of a static body at y = 5 is a floor at y = 5, and a
    // local offset of 1 raises it by one more metre.  (The body position used
    // to be ignored, which left such a floor at y = 0.)
    for (body_y, offset, floor_y) in [(5.0, 0.0, 5.0), (4.0, 1.0, 5.0), (-2.0, 0.0, -2.0)] {
        let mut engine = WasmPhysicsEngine::new(0.0, -9.81, 0.0);
        let floor = engine
            .add_static_body(3.0, body_y, -1.0)
            .expect("floor body");
        engine
            .add_plane_collider(floor, 0.0, 1.0, 0.0, offset)
            .expect("floor plane");
        let ball = add_sphere_body(&mut engine, 1.0, [3.0, floor_y + 3.0, -1.0], 0.5);
        for _ in 0..240 {
            engine.step(DT).expect("step");
        }
        assert_resting_at(&engine, ball, floor_y + 0.5);
    }
}

#[test]
fn plane_colliders_turn_with_their_body() {
    // Turned by 90 degrees about z, the plane y = 0 faces -x: a wall whose
    // free side is x < 0.
    let mut engine = WasmPhysicsEngine::new(0.0, 0.0, 0.0);
    let wall = engine
        .add_rigid_body(&RigidBodyConfig {
            rotation: quat_from_axis_angle([0.0, 0.0, 1.0], FRAC_PI_2),
            ..RigidBodyConfig::static_body()
        })
        .expect("wall body");
    engine
        .add_plane_collider(wall, 0.0, 1.0, 0.0, 0.0)
        .expect("wall plane");
    // The ball reaches 0.2 into the wall and is pushed back toward -x.
    let ball = add_sphere_body(&mut engine, 1.0, [-0.3, 0.0, 0.0], 0.5);
    engine.step(DT).expect("step");
    let contacts = engine.get_contacts();
    assert_eq!(contacts.len(), 1);
    assert_eq!((contacts[0].body_a, contacts[0].body_b), (ball, wall));
    assert_vec_close(contacts[0].normal, [-1.0, 0.0, 0.0], 1e-12, "normal");
    assert_close(contacts[0].depth, 0.2, 1e-12, "depth");
    assert!(
        engine.get_position(ball)[0] < -0.3,
        "pushed out of the wall"
    );
    assert_eq!(engine.get_position(wall), [0.0; 3]);
}

// ---------------------------------------------------------------------------
// Soak test: random piles of every shape
// ---------------------------------------------------------------------------

/// Add a floor and four walls (an 8 m x 8 m container, open at the top) as
/// plane colliders of one static body.
fn add_container(engine: &mut WasmPhysicsEngine) {
    let container = engine.add_static_body(0.0, 0.0, 0.0).expect("container");
    for (normal, offset) in [
        ([0.0, 1.0, 0.0], 0.0),
        ([1.0, 0.0, 0.0], -4.0),
        ([-1.0, 0.0, 0.0], -4.0),
        ([0.0, 0.0, 1.0], -4.0),
        ([0.0, 0.0, -1.0], -4.0),
    ] {
        engine
            .add_plane_collider(container, normal[0], normal[1], normal[2], offset)
            .expect("container plane");
    }
}

/// A random collider of the given kind (`0` sphere, `1` capsule, `2` box,
/// `3` cylinder) of about half a metre in size.
fn random_collider(rng: &mut Lcg, kind: usize) -> ColliderConfig {
    match kind {
        0 => ColliderConfig::sphere(rng.range(0.2, 0.5)),
        1 => ColliderConfig::capsule(rng.range(0.15, 0.35), rng.range(0.3, 0.8)),
        2 => ColliderConfig::cuboid(
            rng.range(0.2, 0.5),
            rng.range(0.2, 0.5),
            rng.range(0.2, 0.5),
        ),
        _ => ColliderConfig::cylinder(rng.range(0.2, 0.4), rng.range(0.3, 0.8)),
    }
}

/// Drop `count` randomly shaped, randomly turned bodies of the given `kinds`
/// into the container for 15 s and check that the simulation never produces
/// a non-finite value or an explosion, that the bodies really collide with
/// each other, and that nothing leaves the container.
fn soak_pile(seed: u64, kinds: &[usize], count: usize) {
    let mut rng = Lcg(seed);
    let mut engine = WasmPhysicsEngine::new(0.0, -9.81, 0.0);
    add_container(&mut engine);
    let mut bodies = Vec::new();
    for index in 0..count {
        let kind = kinds[index % kinds.len()];
        let body = engine
            .add_rigid_body(&RigidBodyConfig {
                rotation: rng.rotation(),
                linear_velocity: rng.vector(2.0),
                ..RigidBodyConfig::new()
                    .with_mass(rng.range(0.5, 3.0))
                    .with_position(
                        rng.range(-3.0, 3.0),
                        1.0 + 0.7 * index as f64 + rng.range(0.0, 0.3),
                        rng.range(-3.0, 3.0),
                    )
            })
            .expect("pile body");
        engine
            .add_collider(body, &random_collider(&mut rng, kind))
            .expect("pile collider");
        bodies.push(body);
    }

    let mut pair_contacts = 0usize;
    for step in 0..900 {
        engine.step(DT).expect("step");
        pair_contacts += engine
            .get_contacts()
            .iter()
            .filter(|contact| bodies.contains(&contact.body_a) && bodies.contains(&contact.body_b))
            .count();
        for &body in &bodies {
            let position = engine.get_position(body);
            let velocity = engine.get_velocity(body);
            let rotation = engine.get_rotation(body);
            let angular = engine.get_angular_velocity(body);
            assert!(
                position
                    .iter()
                    .chain(velocity.iter())
                    .chain(rotation.iter())
                    .chain(angular.iter())
                    .all(|v| v.is_finite()),
                "seed {seed:#x}, step {step}: body {body} left the finite range \
                 (position {position:?}, velocity {velocity:?})"
            );
            assert!(
                length(velocity) < 60.0,
                "seed {seed:#x}, step {step}: body {body} exploded, velocity {velocity:?}"
            );
        }
    }

    assert!(
        pair_contacts > 100,
        "seed {seed:#x}: the bodies hardly touched each other ({pair_contacts} pair contacts)"
    );

    let mut total_speed = 0.0;
    for &body in &bodies {
        let position = engine.get_position(body);
        assert!(
            position[1] > 0.0 && position[0].abs() < 4.0 && position[2].abs() < 4.0,
            "seed {seed:#x}: body {body} left the container, position {position:?}"
        );
        total_speed += length(engine.get_velocity(body));
    }
    assert!(
        total_speed / (bodies.len() as f64) < 2.0,
        "seed {seed:#x}: the pile has not calmed down, mean speed {}",
        total_speed / (bodies.len() as f64)
    );
}

#[test]
fn mixed_shape_piles_stay_finite_and_inside_the_container() {
    for seed in [0x51_u64, 0xA7, 0x1234_5678] {
        soak_pile(seed, &[0, 1, 2, 3], 24);
    }
}

#[test]
fn piles_of_round_and_convex_solids_exercise_the_convex_solver() {
    // Capsule-box, capsule-cylinder, box-cylinder and cylinder-cylinder pairs
    // go through GJK/EPA; none of them may ever produce a non-finite value.
    for seed in [0xC0FFEE_u64, 0xBEEF] {
        soak_pile(seed, &[1, 2, 3], 21);
    }
}
