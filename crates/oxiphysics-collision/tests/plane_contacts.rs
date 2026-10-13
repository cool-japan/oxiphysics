// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for shape-vs-plane contacts of the `ShapeKind` routing layer.
//!
//! A plane `n . x = offset` bounds the solid half-space `n . x <= offset`.
//! Contacts follow the `NarrowPhaseContact` convention: the normal points from
//! shape A towards shape B, `point_a` lies on A and `point_b` on B.

use oxiphysics_collision::{ShapeKind, shape_shape_contact};

const TOL: f64 = 1e-12;

fn ground() -> ShapeKind {
    ShapeKind::Plane {
        normal: [0.0, 1.0, 0.0],
        offset: 0.0,
    }
}

fn unit_box_at(y: f64) -> ShapeKind {
    ShapeKind::Box {
        center: [0.0, y, 0.0],
        half_extents: [0.5, 0.5, 0.5],
    }
}

fn assert_vec(actual: [f64; 3], expected: [f64; 3], what: &str) {
    let gap = (0..3)
        .map(|i| (actual[i] - expected[i]).powi(2))
        .sum::<f64>()
        .sqrt();
    assert!(gap < TOL, "{what}: {actual:?} differs from {expected:?}");
}

#[test]
fn box_straddling_the_plane_penetrates_by_its_lowest_corners() {
    // Bottom face at y = -0.1.
    let c = shape_shape_contact(&unit_box_at(0.4), &ground()).expect("the box dips into the plane");
    assert!((c.depth - 0.1).abs() < TOL, "depth {}", c.depth);
    // A is the box, B the plane: from A towards B is straight down.
    assert_vec(c.normal, [0.0, -1.0, 0.0], "normal");
    // The four lowest corners tie; the witness is their centroid.
    assert_vec(c.point_a, [0.0, -0.1, 0.0], "point_a");
    assert_vec(c.point_b, [0.0, 0.0, 0.0], "point_b");
    // With the plane first the contact is mirrored.
    let m = shape_shape_contact(&ground(), &unit_box_at(0.4)).expect("mirrored contact");
    assert_vec(m.normal, [0.0, 1.0, 0.0], "mirrored normal");
    assert_vec(m.point_a, c.point_b, "mirrored point_a");
    assert_vec(m.point_b, c.point_a, "mirrored point_b");
    assert!((m.depth - c.depth).abs() < TOL);
}

#[test]
fn box_above_the_plane_has_no_contact_and_a_touching_box_has_zero_depth() {
    assert!(shape_shape_contact(&unit_box_at(0.6), &ground()).is_none());
    assert!(shape_shape_contact(&ground(), &unit_box_at(0.6)).is_none());
    let touching = shape_shape_contact(&unit_box_at(0.5), &ground()).expect("touching contact");
    assert!(touching.depth.abs() < TOL, "depth {}", touching.depth);
}

#[test]
fn the_plane_is_unbounded() {
    let far = ShapeKind::Box {
        center: [1.0e6, 0.4, -2.0e6],
        half_extents: [0.5, 0.5, 0.5],
    };
    let c = shape_shape_contact(&far, &ground()).expect("the plane extends to infinity");
    assert!((c.depth - 0.1).abs() < 1e-9, "depth {}", c.depth);
}

#[test]
fn rotated_cube_touches_the_plane_with_its_lowest_vertex() {
    // A cube tilted so that none of its faces is horizontal, given as a convex
    // vertex cloud. Its edges along `u` are horizontal, so the two lowest
    // vertices tie in height.
    let half = 0.5_f64;
    let centre = [0.0, 0.5, 0.0];
    let u = [1.0 / 2.0_f64.sqrt(), 0.0, -1.0 / 2.0_f64.sqrt()];
    let w = [
        1.0 / 6.0_f64.sqrt(),
        (2.0_f64 / 3.0).sqrt(),
        1.0 / 6.0_f64.sqrt(),
    ];
    let v = [
        -1.0 / 3.0_f64.sqrt(),
        1.0 / 3.0_f64.sqrt(),
        -1.0 / 3.0_f64.sqrt(),
    ];
    let mut vertices = Vec::new();
    for sx in [-half, half] {
        for sy in [-half, half] {
            for sz in [-half, half] {
                vertices.push([
                    centre[0] + sx * u[0] + sy * w[0] + sz * v[0],
                    centre[1] + sx * u[1] + sy * w[1] + sz * v[1],
                    centre[2] + sx * u[2] + sy * w[2] + sz * v[2],
                ]);
            }
        }
    }
    let lowest = vertices.iter().map(|p| p[1]).fold(f64::INFINITY, f64::min);
    assert!(lowest < -0.19, "the cube must dip into the plane: {lowest}");
    let raised = ShapeKind::Convex {
        vertices: vertices
            .iter()
            .map(|p| [p[0], p[1] - lowest + 0.01, p[2]])
            .collect(),
    };
    let convex = ShapeKind::Convex { vertices };
    let c = shape_shape_contact(&convex, &ground()).expect("the vertex dips into the plane");
    assert!(
        (c.depth + lowest).abs() < 1e-9,
        "depth {} vs lowest {lowest}",
        c.depth
    );
    assert!(
        (c.point_a[1] - lowest).abs() < 1e-9,
        "point_a {:?}",
        c.point_a
    );
    assert!(c.point_b[1].abs() < 1e-9);
    assert_vec(c.normal, [0.0, -1.0, 0.0], "normal");
    // The same cloud lifted 0.01 clear of the plane.
    assert!(shape_shape_contact(&raised, &ground()).is_none());
}

#[test]
fn capsule_contacts_the_plane_with_its_lowest_end_cap() {
    // Upright capsule: lowest surface point is p0 - radius.
    let upright = ShapeKind::Capsule {
        p0: [0.0, 0.3, 0.0],
        p1: [0.0, 1.3, 0.0],
        radius: 0.5,
    };
    let c = shape_shape_contact(&upright, &ground()).expect("capsule dips into the plane");
    assert!((c.depth - 0.2).abs() < TOL, "depth {}", c.depth);
    assert_vec(c.normal, [0.0, -1.0, 0.0], "normal");
    assert_vec(c.point_a, [0.0, -0.2, 0.0], "point_a");
    assert_vec(c.point_b, [0.0, 0.0, 0.0], "point_b");
    // Lying capsule: both end caps tie, the witness is their centroid.
    let lying = ShapeKind::Capsule {
        p0: [-1.0, 0.4, 0.0],
        p1: [1.0, 0.4, 0.0],
        radius: 0.5,
    };
    let c = shape_shape_contact(&lying, &ground()).expect("lying capsule dips into the plane");
    assert!((c.depth - 0.1).abs() < TOL, "depth {}", c.depth);
    assert_vec(c.point_a, [0.0, -0.1, 0.0], "point_a");
    // Tilted capsule: the lower end decides.
    let tilted = ShapeKind::Capsule {
        p0: [0.0, 0.3, 0.0],
        p1: [1.0, 1.3, 0.0],
        radius: 0.5,
    };
    let c = shape_shape_contact(&tilted, &ground()).expect("tilted capsule dips into the plane");
    assert!((c.depth - 0.2).abs() < TOL, "depth {}", c.depth);
    assert_vec(c.point_a, [0.0, -0.2, 0.0], "point_a");
    // Reversed argument order and a capsule that is clear of the plane.
    let m = shape_shape_contact(&ground(), &upright).expect("mirrored contact");
    assert_vec(m.normal, [0.0, 1.0, 0.0], "mirrored normal");
    let high = ShapeKind::Capsule {
        p0: [0.0, 0.6, 0.0],
        p1: [0.0, 1.6, 0.0],
        radius: 0.5,
    };
    assert!(shape_shape_contact(&high, &ground()).is_none());
}

#[test]
fn sphere_plane_normal_points_from_the_sphere_to_the_plane() {
    let sphere = ShapeKind::Sphere {
        center: [0.0, 0.3, 0.0],
        radius: 0.5,
    };
    let c = shape_shape_contact(&sphere, &ground()).expect("sphere dips into the plane");
    assert!((c.depth - 0.2).abs() < TOL, "depth {}", c.depth);
    assert_vec(c.normal, [0.0, -1.0, 0.0], "normal");
    assert_vec(c.point_a, [0.0, -0.2, 0.0], "point_a");
    assert_vec(c.point_b, [0.0, 0.0, 0.0], "point_b");
    let m = shape_shape_contact(&ground(), &sphere).expect("mirrored contact");
    assert_vec(m.normal, [0.0, 1.0, 0.0], "mirrored normal");
    // The box, capsule and sphere routines agree on the sign for a sphere-sized
    // box resting on the same plane.
    let cube = shape_shape_contact(&unit_box_at(0.3), &ground()).expect("box contact");
    assert_vec(cube.normal, c.normal, "box and sphere normals agree");
}

#[test]
fn tilted_and_offset_planes_are_normalised() {
    // Plane through y = 1 given by 2y = 2.
    let plane = ShapeKind::Plane {
        normal: [0.0, 2.0, 0.0],
        offset: 2.0,
    };
    let c = shape_shape_contact(&unit_box_at(1.3), &plane).expect("box dips into the plane");
    assert!((c.depth - 0.2).abs() < TOL, "depth {}", c.depth);
    assert_vec(c.point_b, [0.0, 1.0, 0.0], "point_b");
    // A 45-degree plane x + y = 0: the box corners at (-0.9, -0.9, +-0.5) dip
    // deepest.
    let slanted = ShapeKind::Plane {
        normal: [1.0, 1.0, 0.0],
        offset: 0.0,
    };
    let b = ShapeKind::Box {
        center: [-0.4, -0.4, 0.0],
        half_extents: [0.5, 0.5, 0.5],
    };
    let c = shape_shape_contact(&b, &slanted).expect("box crosses the slanted plane");
    // Signed distance of the lowest corner (-0.9, -0.9): -1.8 / sqrt(2).
    assert!(
        (c.depth - 1.8 / 2.0_f64.sqrt()).abs() < 1e-9,
        "depth {}",
        c.depth
    );
}

#[test]
fn two_planes_have_no_contact() {
    assert!(shape_shape_contact(&ground(), &ground()).is_none());
}

#[test]
fn batch_narrow_phase_resolves_plane_pairs() {
    use oxiphysics_collision::BatchNarrowPhase;
    let shapes = [unit_box_at(0.4), ground(), unit_box_at(2.0)];
    let batch = BatchNarrowPhase::new();
    let results = batch.run(&shapes, &[(0, 1), (2, 1), (1, 0), (1, 1)]);
    assert!(results[0].is_some(), "box dipping into the plane");
    assert!(results[1].is_none(), "box clear of the plane");
    assert!(results[2].is_some(), "plane first, box second");
    assert!(results[3].is_none(), "a plane never touches itself");
    let compact = batch.run_compact(&shapes, &[(0, 1), (2, 1)]);
    assert_eq!(compact.len(), 1);
    assert_eq!(compact[0].0, (0, 1));
    assert!((compact[0].1.depth - 0.1).abs() < TOL);
}

#[test]
fn plane_support_point_depends_on_the_direction() {
    let plane = ShapeKind::Plane {
        normal: [0.0, 1.0, 0.0],
        offset: 2.0,
    };
    // Towards the outside of the half-space the support point stays on the plane.
    let up = plane.support([0.0, 1.0, 0.0]);
    assert!((up[1] - 2.0).abs() < TOL, "{up:?}");
    // Into the half-space it runs to the far side of the stand-in slab.
    let down = plane.support([0.0, -1.0, 0.0]);
    assert!(down[1] < -1.0e8, "{down:?}");
    // Along the plane it slides towards the requested side.
    let slid = plane.support([1.0, 1.0, 0.0]);
    assert!(slid[0] > 1.0e8 && (slid[1] - 2.0).abs() < TOL, "{slid:?}");
    let slid_back = plane.support([-1.0, 1.0, 0.0]);
    assert!(slid_back[0] < -1.0e8, "{slid_back:?}");
}
