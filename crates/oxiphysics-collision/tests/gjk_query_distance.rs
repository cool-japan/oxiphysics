// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the converged distance of separated GJK queries.
//!
//! `Gjk::query` used to stop at the first separating axis and report the
//! distance and witness points of the terminating simplex, an upper bound that
//! is generally not the closest pair of the shapes. Conservative advancement
//! (`time_of_impact`) took steps of that overestimated length, and the contact
//! normal and witnesses it returned were built from the unconverged pair.

use oxiphysics_collision::narrowphase::gjk::{GjkResult, gjk_distance_query};
use oxiphysics_collision::{Gjk, time_of_impact};
use oxiphysics_core::Transform;
use oxiphysics_core::math::Vec3;
use oxiphysics_geometry::{BoxShape, Capsule, Sphere};

fn at(x: f64, y: f64, z: f64) -> Transform {
    Transform::from_position(Vec3::new(x, y, z))
}

fn separated(result: GjkResult) -> (f64, Vec3, Vec3) {
    match result {
        GjkResult::Separated {
            point_a,
            point_b,
            distance,
        } => (distance, point_a, point_b),
        GjkResult::Intersecting(_) => panic!("the shapes are separated"),
    }
}

#[test]
fn unit_spheres_five_apart_report_the_gap_and_the_facing_poles() {
    let sphere = Sphere::new(1.0);
    let (distance, point_a, point_b) = separated(Gjk::query(
        &sphere,
        &at(0.0, 0.0, 0.0),
        &sphere,
        &at(5.0, 0.0, 0.0),
    ));
    assert!((distance - 3.0).abs() < 1e-9, "distance {distance}");
    assert!(
        (point_a - Vec3::new(1.0, 0.0, 0.0)).norm() < 1e-9,
        "{point_a:?}"
    );
    assert!(
        (point_b - Vec3::new(4.0, 0.0, 0.0)).norm() < 1e-9,
        "{point_b:?}"
    );
}

#[test]
fn staggered_unit_boxes_are_two_apart() {
    // Boxes of half-extent 1 at the origin and at (4, 0.5, 0): the x faces are
    // two apart and the y/z ranges of the facing sides overlap.
    let cube = BoxShape::new(Vec3::new(1.0, 1.0, 1.0));
    let (distance, point_a, point_b) = separated(Gjk::query(
        &cube,
        &at(0.0, 0.0, 0.0),
        &cube,
        &at(4.0, 0.5, 0.0),
    ));
    assert!((distance - 2.0).abs() < 1e-9, "distance {distance}");
    assert!((point_a.x - 1.0).abs() < 1e-9 && (point_b.x - 3.0).abs() < 1e-9);
    assert!(((point_a - point_b).norm() - distance).abs() < 1e-9);
    // The same pair through the detailed query.
    let detailed = gjk_distance_query(&cube, &at(0.0, 0.0, 0.0), &cube, &at(4.0, 0.5, 0.0));
    assert!((detailed.distance - 2.0).abs() < 1e-9);
    assert!((detailed.closest_a.x - 1.0).abs() < 1e-9);
    assert!((detailed.closest_b.x - 3.0).abs() < 1e-9);
}

#[test]
fn a_tilted_capsule_over_a_slab_reports_the_distance_to_its_lower_cap() {
    // The capsule axis is turned 30 degrees about z, so its lower cap centre sits
    // at (1.3 + 0.8 sin 30, 2 - 0.8 cos 30, -0.7); the nearest point of the
    // capsule is the bottom of that cap sphere, straight above the slab's top.
    let slab = BoxShape::new(Vec3::new(5.0, 0.5, 5.0));
    let capsule = Capsule::new(0.4, 0.8);
    let tilt = std::f64::consts::FRAC_PI_6;
    let pose = Transform::from_axis_angle(Vec3::new(1.3, 2.0, -0.7), Vec3::z(), tilt);
    let (distance, point_a, point_b) =
        separated(Gjk::query(&slab, &at(0.0, 0.0, 0.0), &capsule, &pose));
    let cap_x = 1.3 + 0.8 * tilt.sin();
    let lowest = 2.0 - 0.8 * tilt.cos() - 0.4;
    assert!(
        (distance - (lowest - 0.5)).abs() < 1e-8,
        "distance {distance} vs {}",
        lowest - 0.5
    );
    assert!((point_a.y - 0.5).abs() < 1e-8, "{point_a:?}");
    assert!((point_b.y - lowest).abs() < 1e-7, "{point_b:?}");
    // On a curved cap the witnesses are only fixed to about the square root of
    // the distance tolerance along the surface.
    assert!((point_a.x - cap_x).abs() < 1e-3 && (point_a.z + 0.7).abs() < 1e-3);
    assert!((point_b.x - cap_x).abs() < 1e-3 && (point_b.z + 0.7).abs() < 1e-3);
}

#[test]
fn time_of_impact_advances_by_the_true_gap_and_reports_the_face_normal() {
    // A box flies at a second box whose centre is offset in y, so the first
    // support point of the search is far from the closest one. They touch when
    // the faces meet at x = 4, half way through the step.
    let cube = BoxShape::new(Vec3::new(1.0, 1.0, 1.0));
    let start = at(0.0, 0.0, 0.0);
    let end = at(6.0, 0.0, 0.0);
    let target = at(5.0, 0.5, 0.0);
    let hit = time_of_impact(&cube, &start, &end, &cube, &target, &target)
        .expect("the box reaches the target within the step");
    assert!((hit.toi - 0.5).abs() < 1e-5, "toi {}", hit.toi);
    // The centre-to-centre direction (5, 0.5, 0) is 0.1 away from the true normal.
    assert!(
        (hit.normal - Vec3::new(1.0, 0.0, 0.0)).norm() < 1e-2,
        "normal {:?}",
        hit.normal
    );
    // Both witnesses sit on the meeting plane x = 4 and on the shared part of
    // the faces (y in [-0.5, 1], z in [-1, 1]), not at the body origins.
    assert!((hit.witness_a.x - 4.0).abs() < 1e-5, "{:?}", hit.witness_a);
    assert!((hit.witness_b.x - 4.0).abs() < 1e-5, "{:?}", hit.witness_b);
    assert!(hit.witness_b.y >= -0.5 - 1e-6 && hit.witness_b.y <= 1.0 + 1e-6);
    assert!(hit.witness_b.z.abs() <= 1.0 + 1e-6);
}

#[test]
fn a_fast_sphere_does_not_tunnel_through_a_thin_plate() {
    // The plate is 0.1 thick; the sphere moves 10 units in one step, offset in y
    // from the plate centre.
    let plate = BoxShape::new(Vec3::new(0.05, 2.0, 2.0));
    let sphere = Sphere::new(0.25);
    let start = at(-5.0, 0.3, 0.2);
    let end = at(5.0, 0.3, 0.2);
    let rest = at(0.0, 0.0, 0.0);
    let hit = time_of_impact(&sphere, &start, &end, &plate, &rest, &rest)
        .expect("the sphere hits the plate");
    // The sphere's front reaches the plate face x = -0.05 when its centre is at
    // x = -0.3, i.e. after 4.7 of the 10 units.
    assert!((hit.toi - 0.47).abs() < 1e-4, "toi {}", hit.toi);
    assert!(hit.normal.x > 0.99, "normal {:?}", hit.normal);
}
