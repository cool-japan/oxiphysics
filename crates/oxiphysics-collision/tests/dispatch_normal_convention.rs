// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the crate-wide contact-normal convention.
//!
//! Every manifold produced by [`NarrowPhaseDispatcher::dispatch`] must carry a
//! normal pointing FROM body B TOWARD body A (see [`Contact::normal`]), with
//! `point_a` on A, `point_b` on B and a positive penetration depth equal to the
//! true minimum overlap. Each primitive pair is stacked with A directly above B
//! and an overlap of 0.05 m, and checked in BOTH argument orders — specialised
//! closures (sphere/sphere, capsule/capsule, ...), the GJK/EPA fallback and the
//! convex-hull face-clip manifold (which is seeded by the EPA normal) alike.

use oxiphysics_collision::narrowphase::dispatch::ShapeType;
use oxiphysics_collision::narrowphase::epa::Epa;
use oxiphysics_collision::narrowphase::gjk::{Gjk, GjkResult, mpr_contact};
use oxiphysics_collision::types::CollisionPair;
use oxiphysics_collision::{Contact, NarrowPhaseDispatcher};
use oxiphysics_core::Transform;
use oxiphysics_core::math::Vec3;
use oxiphysics_geometry::{BoxShape, Capsule, Cone, ConvexHull, Cylinder, Shape, Sphere};

/// Expected overlap of every stacked configuration below.
const OVERLAP: f64 = 0.05;

/// A unit cube (half-extent 0.5) expressed as a convex hull centred at the origin.
fn cube_hull() -> ConvexHull {
    let mut v = Vec::new();
    for sx in [-0.5, 0.5] {
        for sy in [-0.5, 0.5] {
            for sz in [-0.5, 0.5] {
                v.push(Vec3::new(sx, sy, sz));
            }
        }
    }
    ConvexHull::new(v)
}

/// Half-height (extent along +y) of each probe shape.
fn half_y(kind: ShapeType) -> f64 {
    match kind {
        ShapeType::Capsule => 0.75,
        _ => 0.5,
    }
}

fn make(kind: ShapeType) -> Box<dyn Shape> {
    match kind {
        ShapeType::Sphere => Box::new(Sphere::new(0.5)),
        ShapeType::Box => Box::new(BoxShape::new(Vec3::new(0.5, 0.5, 0.5))),
        ShapeType::Capsule => Box::new(Capsule::new(0.25, 0.5)),
        ShapeType::Cylinder => Box::new(Cylinder::new(0.5, 0.5)),
        ShapeType::Cone => Box::new(Cone::new(0.5, 0.5)),
        _ => Box::new(cube_hull()),
    }
}

/// Dispatch A (above) against B (below) and return the manifold contacts.
fn stacked_contacts(kind_a: ShapeType, kind_b: ShapeType) -> Vec<Contact> {
    let dispatcher = NarrowPhaseDispatcher::default();
    let shape_a = make(kind_a);
    let shape_b = make(kind_b);
    let y_a = half_y(kind_a) + half_y(kind_b) - OVERLAP;
    let t_a = Transform::from_position(Vec3::new(0.0, y_a, 0.0));
    let t_b = Transform::from_position(Vec3::zeros());
    let result = dispatcher.dispatch(
        shape_a.as_ref(),
        kind_a,
        &t_a,
        shape_b.as_ref(),
        kind_b,
        &t_b,
        CollisionPair::new(0, 1),
    );
    let manifold = result
        .manifold
        .unwrap_or_else(|| panic!("{kind_a:?} above {kind_b:?} must overlap"));
    assert!(
        !manifold.contacts.is_empty(),
        "{kind_a:?}/{kind_b:?}: empty manifold"
    );
    manifold.contacts
}

fn assert_b_to_a(kind_a: ShapeType, kind_b: ShapeType) {
    for c in stacked_contacts(kind_a, kind_b) {
        assert!(
            c.normal.y > 0.99,
            "{kind_a:?} above {kind_b:?}: normal must point from B up to A, got {:?}",
            c.normal
        );
        assert!(
            (c.depth - OVERLAP).abs() < 5e-3,
            "{kind_a:?} above {kind_b:?}: depth {} != overlap {OVERLAP}",
            c.depth
        );
        // A sits above B, so the witness on A lies at or below the witness on B.
        assert!(
            c.point_a.y <= c.point_b.y + 1e-6,
            "{kind_a:?}/{kind_b:?}: point_a {:?} must lie below point_b {:?}",
            c.point_a,
            c.point_b
        );
    }
}

const KINDS: [ShapeType; 6] = [
    ShapeType::Sphere,
    ShapeType::Box,
    ShapeType::Capsule,
    ShapeType::Cylinder,
    ShapeType::Cone,
    ShapeType::ConvexHull,
];

#[test]
fn dispatch_normals_point_from_b_to_a_for_every_pair_and_order() {
    for &a in &KINDS {
        for &b in &KINDS {
            assert_b_to_a(a, b);
        }
    }
}

#[test]
fn hull_box_manifold_has_true_depth_and_four_points() {
    // The face-clip manifold is seeded by the EPA normal; with a wrong seed sign
    // it reported depth 1.95 and witness points on the far faces.
    let contacts = stacked_contacts(ShapeType::ConvexHull, ShapeType::Box);
    assert_eq!(contacts.len(), 4, "face-on-face hull/box gives 4 corners");
    for c in &contacts {
        assert!((c.point_a.y - 0.45).abs() < 1e-6, "A bottom face at 0.45");
        assert!((c.point_b.y - 0.5).abs() < 1e-6, "B top face at 0.5");
    }
}

#[test]
fn epa_and_mpr_report_b_to_a_normals() {
    // Box A centred at x = 0.9 overlaps box B at the origin by 0.1 along x.
    let a = BoxShape::new(Vec3::new(0.5, 0.5, 0.5));
    let b = BoxShape::new(Vec3::new(0.5, 0.5, 0.5));
    let t_a = Transform::from_position(Vec3::new(0.9, 0.0, 0.0));
    let t_b = Transform::from_position(Vec3::zeros());
    let simplex = match Gjk::query(&a, &t_a, &b, &t_b) {
        GjkResult::Intersecting(s) => s,
        GjkResult::Separated { .. } => panic!("boxes overlap"),
    };
    let epa = Epa::penetration_depth(&a, &t_a, &b, &t_b, &simplex).expect("EPA contact");
    assert!(
        epa.normal.x > 0.99,
        "EPA normal {:?} must be +x",
        epa.normal
    );
    assert!((epa.depth - 0.1).abs() < 1e-3, "EPA depth {}", epa.depth);
    let mpr = mpr_contact(&a, &t_a, &b, &t_b).expect("MPR contact");
    assert!(
        mpr.normal.x > 0.99,
        "MPR normal {:?} must be +x",
        mpr.normal
    );
    assert!((mpr.depth - 0.1).abs() < 1e-3, "MPR depth {}", mpr.depth);
}

#[test]
fn gjk_query_separated_witnesses_lie_on_the_shapes() {
    // Unit box 0.1 m above a wide ground slab: the witness points must lie on
    // the facing surfaces and their difference must be the separation vector.
    let a = BoxShape::new(Vec3::new(0.5, 0.5, 0.5));
    let b = BoxShape::new(Vec3::new(20.0, 0.5, 20.0));
    let t_a = Transform::from_position(Vec3::new(0.0, 0.6, 0.0));
    let t_b = Transform::from_position(Vec3::new(0.0, -0.5, 0.0));
    match Gjk::query(&a, &t_a, &b, &t_b) {
        GjkResult::Separated {
            point_a,
            point_b,
            distance,
        } => {
            assert!((distance - 0.1).abs() < 1e-9, "distance {distance}");
            assert!((point_a.y - 0.1).abs() < 1e-9, "A bottom face: {point_a:?}");
            assert!(point_b.y.abs() < 1e-9, "B top face: {point_b:?}");
            assert!(
                ((point_a - point_b).norm() - distance).abs() < 1e-9,
                "witness gap must equal the distance"
            );
        }
        GjkResult::Intersecting(_) => panic!("shapes are 0.1 m apart"),
    }
}
