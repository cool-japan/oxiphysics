// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for [`dispatch_concave_vs_convex`].
//!
//! The previous implementation ignored the convex shape: it intersected every
//! convex piece of the mesh with itself, reported a fabricated depth of `0.1`
//! with the piece centroid as the contact point, culled pieces with a hard-coded
//! `2.0` window around the shape's centre and numbered the manifolds with
//! made-up body ids (`b + piece * 100`). These tests exercise the real
//! piece-versus-shape dispatch.

use oxiphysics_collision::NarrowPhaseDispatcher;
use oxiphysics_collision::narrowphase::dispatch::{
    ConcaveMesh, ConvexPiece, ShapeType, dispatch_concave_vs_convex,
};
use oxiphysics_collision::types::{CollisionPair, ContactManifold, MAX_CONTACTS};
use oxiphysics_core::Transform;
use oxiphysics_core::math::Vec3;
use oxiphysics_geometry::{BoxShape, Shape, Sphere};

/// Axis-aligned box as a convex piece with its eight corners.
fn box_piece(min: [f64; 3], max: [f64; 3]) -> ConvexPiece {
    let mut corners = Vec::with_capacity(8);
    for x in [min[0], max[0]] {
        for y in [min[1], max[1]] {
            for z in [min[2], max[2]] {
                corners.push([x, y, z]);
            }
        }
    }
    ConvexPiece::new(corners)
}

/// Two-piece slab. Its top face is the plane `y = 0` over `x in [-2, 2]`,
/// `z in [-1, 1]`, with the seam between the pieces at `x = 0`.
fn slab() -> ConcaveMesh {
    ConcaveMesh::new(vec![
        box_piece([-2.0, -1.0, -1.0], [0.0, 0.0, 1.0]),
        box_piece([0.0, -1.0, -1.0], [2.0, 0.0, 1.0]),
    ])
}

fn identity() -> Transform {
    Transform::from_position(Vec3::zeros())
}

fn dispatch(
    mesh: &ConcaveMesh,
    mesh_transform: &Transform,
    shape: &dyn Shape,
    kind: ShapeType,
    shape_transform: &Transform,
    pair: CollisionPair,
) -> Vec<ContactManifold> {
    let dispatcher = NarrowPhaseDispatcher::default();
    dispatch_concave_vs_convex(
        mesh,
        mesh_transform,
        shape,
        kind,
        shape_transform,
        pair,
        &dispatcher,
    )
}

fn sphere_at(mesh: &ConcaveMesh, center: [f64; 3], radius: f64) -> Vec<ContactManifold> {
    let sphere = Sphere::new(radius);
    let transform = Transform::from_position(Vec3::new(center[0], center[1], center[2]));
    dispatch(
        mesh,
        &identity(),
        &sphere,
        ShapeType::Sphere,
        &transform,
        CollisionPair::new(0, 1),
    )
}

#[test]
fn sphere_sunk_into_a_piece_reports_the_geometric_depth() {
    // Sphere radius 0.5 centred 0.3 above the top face: its bottom is 0.2 below.
    let manifolds = sphere_at(&slab(), [-1.0, 0.3, 0.0], 0.5);
    assert_eq!(manifolds.len(), 1, "one manifold for the whole mesh");
    let manifold = &manifolds[0];
    assert!(!manifold.contacts.is_empty());
    for contact in &manifold.contacts {
        assert!(
            (contact.depth - 0.2).abs() < 1e-3,
            "depth {}",
            contact.depth
        );
        // Normal from the sphere (B) toward the mesh (A): straight down.
        assert!(
            (contact.normal - Vec3::new(0.0, -1.0, 0.0)).norm() < 1e-3,
            "normal {:?}",
            contact.normal
        );
        // The witness on the mesh lies on its top face, the one on the sphere at its bottom.
        assert!(
            contact.point_a.y.abs() < 1e-3,
            "point_a {:?}",
            contact.point_a
        );
        assert!(
            (contact.point_b.y + 0.2).abs() < 1e-3,
            "point_b {:?}",
            contact.point_b
        );
        assert!((contact.point_a.x + 1.0).abs() < 1e-2);
    }
}

#[test]
fn sphere_above_the_slab_has_no_contact() {
    // Bottom of the sphere is 0.5 above the face.
    assert!(sphere_at(&slab(), [-1.0, 1.0, 0.0], 0.5).is_empty());
}

#[test]
fn sphere_inside_the_old_fabricated_window_but_clear_of_the_mesh_has_no_contact() {
    // The old code culled with a fixed 2.0 window around the centre, then tested
    // the piece against itself and reported a contact of depth 0.1.
    assert!(sphere_at(&slab(), [-1.0, 1.9, 0.0], 0.1).is_empty());
    assert!(sphere_at(&slab(), [-1.0, 0.0, 3.5], 1.0).is_empty());
}

#[test]
fn sphere_far_from_every_piece_is_culled() {
    assert!(sphere_at(&slab(), [100.0, 0.0, 0.0], 0.5).is_empty());
    assert!(sphere_at(&slab(), [0.0, -50.0, 0.0], 0.5).is_empty());
}

#[test]
fn sphere_on_the_seam_yields_one_manifold_with_the_base_pair() {
    let manifolds = sphere_at(&slab(), [0.0, 0.3, 0.0], 0.5);
    assert_eq!(manifolds.len(), 1);
    let manifold = &manifolds[0];
    assert!(!manifold.contacts.is_empty() && manifold.contacts.len() <= MAX_CONTACTS);
    assert_eq!(manifold.pair, CollisionPair::new(0, 1));
    for contact in &manifold.contacts {
        assert!(
            (contact.depth - 0.2).abs() < 1e-3,
            "depth {}",
            contact.depth
        );
    }
}

#[test]
fn manifold_keeps_the_caller_pair_instead_of_inventing_body_ids() {
    let sphere = Sphere::new(0.5);
    let transform = Transform::from_position(Vec3::new(1.0, 0.3, 0.0));
    let manifolds = dispatch(
        &slab(),
        &identity(),
        &sphere,
        ShapeType::Sphere,
        &transform,
        CollisionPair::new(7, 9),
    );
    assert_eq!(manifolds.len(), 1);
    assert_eq!(manifolds[0].pair, CollisionPair::new(7, 9));
}

#[test]
fn box_straddling_the_seam_gets_contacts_on_both_pieces() {
    let shape = BoxShape::new(Vec3::new(0.5, 0.5, 0.5));
    // Bottom face 0.1 below the slab top, spanning x in [-0.5, 0.5].
    let transform = Transform::from_position(Vec3::new(0.0, 0.4, 0.0));
    let manifolds = dispatch(
        &slab(),
        &identity(),
        &shape,
        ShapeType::Box,
        &transform,
        CollisionPair::new(3, 4),
    );
    assert_eq!(manifolds.len(), 1);
    let manifold = &manifolds[0];
    assert_eq!(manifold.pair, CollisionPair::new(3, 4));
    assert!(manifold.contacts.len() <= MAX_CONTACTS);
    assert!(
        manifold.contacts.len() >= 3,
        "{} contacts",
        manifold.contacts.len()
    );
    for contact in &manifold.contacts {
        assert!(
            (contact.depth - 0.1).abs() < 1e-6,
            "depth {}",
            contact.depth
        );
        assert!(
            (contact.normal - Vec3::new(0.0, -1.0, 0.0)).norm() < 1e-6,
            "normal {:?}",
            contact.normal
        );
    }
    // The reduction keeps the footprint: points on both sides of the seam.
    assert!(manifold.contacts.iter().any(|c| c.point_a.x < -0.2));
    assert!(manifold.contacts.iter().any(|c| c.point_a.x > 0.2));
    let z_min = manifold
        .contacts
        .iter()
        .map(|c| c.point_a.z)
        .fold(f64::INFINITY, f64::min);
    let z_max = manifold
        .contacts
        .iter()
        .map(|c| c.point_a.z)
        .fold(f64::NEG_INFINITY, f64::max);
    assert!(z_max - z_min > 0.5, "contacts span z in [{z_min}, {z_max}]");
}

#[test]
fn box_resting_on_one_piece_only_uses_that_piece() {
    let shape = BoxShape::new(Vec3::new(0.25, 0.5, 0.25));
    let transform = Transform::from_position(Vec3::new(1.0, 0.45, 0.0));
    let manifolds = dispatch(
        &slab(),
        &identity(),
        &shape,
        ShapeType::Box,
        &transform,
        CollisionPair::new(0, 1),
    );
    assert_eq!(manifolds.len(), 1);
    for contact in &manifolds[0].contacts {
        assert!(
            (contact.depth - 0.05).abs() < 1e-6,
            "depth {}",
            contact.depth
        );
        assert!(
            contact.point_a.x > 0.5 && contact.point_a.x < 1.5,
            "{:?}",
            contact.point_a
        );
    }
}

#[test]
fn mesh_transform_moves_and_rotates_the_pieces() {
    // Rotate the slab a quarter turn about z and move it to x = 10: its local +y
    // (the top face normal) points along world -x and the top face is at x = 10.
    let mesh_transform = Transform::from_axis_angle(
        Vec3::new(10.0, 0.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
        std::f64::consts::FRAC_PI_2,
    );
    // Local point (-1, 0.3, 0) lies 0.3 above the top face: world (9.7, -1, 0).
    let sphere = Sphere::new(0.5);
    let sphere_transform = Transform::from_position(Vec3::new(9.7, -1.0, 0.0));
    let manifolds = dispatch(
        &slab(),
        &mesh_transform,
        &sphere,
        ShapeType::Sphere,
        &sphere_transform,
        CollisionPair::new(0, 1),
    );
    assert_eq!(manifolds.len(), 1);
    for contact in &manifolds[0].contacts {
        assert!(
            (contact.depth - 0.2).abs() < 1e-3,
            "depth {}",
            contact.depth
        );
        // From the sphere toward the slab: world +x.
        assert!(
            (contact.normal - Vec3::new(1.0, 0.0, 0.0)).norm() < 1e-3,
            "normal {:?}",
            contact.normal
        );
        assert!(
            (contact.point_a.x - 10.0).abs() < 1e-3,
            "point_a {:?}",
            contact.point_a
        );
    }
    // The same sphere placed at the un-transformed location collides with nothing.
    let away = Transform::from_position(Vec3::new(-1.0, 0.3, 0.0));
    assert!(
        dispatch(
            &slab(),
            &mesh_transform,
            &sphere,
            ShapeType::Sphere,
            &away,
            CollisionPair::new(0, 1)
        )
        .is_empty()
    );
}

#[test]
fn rotated_convex_shape_is_culled_with_its_rotated_bounds() {
    // A long thin box standing upright over the slab, then lying down: only the
    // orientation changes which pieces it reaches.
    let shape = BoxShape::new(Vec3::new(0.1, 1.5, 0.1));
    let standing = Transform::from_position(Vec3::new(-1.0, 1.4, 0.0));
    let hit = dispatch(
        &slab(),
        &identity(),
        &shape,
        ShapeType::Box,
        &standing,
        CollisionPair::new(0, 1),
    );
    assert_eq!(hit.len(), 1, "box bottom at -0.1 penetrates the slab");

    let lying = Transform::from_axis_angle(
        Vec3::new(-1.0, 1.4, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
        std::f64::consts::FRAC_PI_2,
    );
    let miss = dispatch(
        &slab(),
        &identity(),
        &shape,
        ShapeType::Box,
        &lying,
        CollisionPair::new(0, 1),
    );
    assert!(
        miss.is_empty(),
        "lying box (centre 1.4, half thickness 0.1) is clear of the slab"
    );
}

#[test]
fn empty_pieces_and_empty_meshes_are_ignored() {
    let sphere = Sphere::new(0.5);
    let transform = Transform::from_position(Vec3::new(1.0, 0.3, 0.0));
    let pair = CollisionPair::new(0, 1);

    let empty = ConcaveMesh::new(Vec::new());
    assert!(
        dispatch(
            &empty,
            &identity(),
            &sphere,
            ShapeType::Sphere,
            &transform,
            pair
        )
        .is_empty()
    );

    let mut pieces = vec![ConvexPiece::new(Vec::new())];
    pieces.extend(slab().pieces);
    let with_empty_piece = ConcaveMesh::new(pieces);
    let manifolds = dispatch(
        &with_empty_piece,
        &identity(),
        &sphere,
        ShapeType::Sphere,
        &transform,
        pair,
    );
    assert_eq!(manifolds.len(), 1);
}

#[test]
fn flat_piece_against_a_box_does_not_panic() {
    // A zero-thickness quad in the z = 0 plane.
    let quad = ConvexPiece::new(vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],
    ]);
    let mesh = ConcaveMesh::new(vec![quad]);
    let shape = BoxShape::new(Vec3::new(0.5, 0.5, 0.5));
    let transform = Transform::from_position(Vec3::new(0.5, 0.5, 0.2));
    let manifolds = dispatch(
        &mesh,
        &identity(),
        &shape,
        ShapeType::Box,
        &transform,
        CollisionPair::new(0, 1),
    );
    for manifold in &manifolds {
        for contact in &manifold.contacts {
            assert!(contact.depth.is_finite() && contact.depth >= 0.0);
            assert!(contact.normal.iter().all(|c| c.is_finite()));
            assert!(contact.point_a.iter().all(|c| c.is_finite()));
            assert!(contact.point_b.iter().all(|c| c.is_finite()));
        }
    }
}

#[test]
fn many_pieces_are_reduced_to_the_manifold_capacity() {
    // A row of ten unit cubes; a wide box lies across all of them.
    let pieces: Vec<ConvexPiece> = (0..10)
        .map(|k| {
            let x = k as f64;
            box_piece([x, -1.0, -1.0], [x + 1.0, 0.0, 1.0])
        })
        .collect();
    let mesh = ConcaveMesh::new(pieces);
    let shape = BoxShape::new(Vec3::new(5.0, 0.5, 0.5));
    let transform = Transform::from_position(Vec3::new(5.0, 0.4, 0.0));
    let manifolds = dispatch(
        &mesh,
        &identity(),
        &shape,
        ShapeType::Box,
        &transform,
        CollisionPair::new(1, 2),
    );
    assert_eq!(manifolds.len(), 1);
    let contacts = &manifolds[0].contacts;
    assert!(contacts.len() <= MAX_CONTACTS);
    let x_min = contacts
        .iter()
        .map(|c| c.point_a.x)
        .fold(f64::INFINITY, f64::min);
    let x_max = contacts
        .iter()
        .map(|c| c.point_a.x)
        .fold(f64::NEG_INFINITY, f64::max);
    assert!(x_max - x_min > 8.0, "contacts span x in [{x_min}, {x_max}]");
}
