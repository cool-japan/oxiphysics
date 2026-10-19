// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for `build_box_manifold` through its public paths.
//!
//! The function used to offset the reference corners along the face normal,
//! clip in an arbitrary 2-D frame and unproject onto the plane through the box
//! centre, so every contact got the half-extent as its depth. These tests pin
//! the real face-clipping manifold: contact points on the overlap patch, the
//! geometric depth of each point, and the B-to-A normal convention.

use oxiphysics_collision::narrowphase::build_box_manifold;
use oxiphysics_collision::narrowphase::manifold::{
    ContactManifold, build_box_manifold as via_manifold_module, reduce_manifold,
};

type Matrix4 = [[f64; 4]; 4];

fn pose(centre: [f64; 3]) -> Matrix4 {
    [
        [1.0, 0.0, 0.0, centre[0]],
        [0.0, 1.0, 0.0, centre[1]],
        [0.0, 0.0, 1.0, centre[2]],
        [0.0, 0.0, 0.0, 1.0],
    ]
}

/// Rotation by `angle` about the vertical axis, row-major.
fn yawed(centre: [f64; 3], angle: f64) -> Matrix4 {
    let (s, c) = angle.sin_cos();
    [
        [c, 0.0, s, centre[0]],
        [0.0, 1.0, 0.0, centre[1]],
        [-s, 0.0, c, centre[2]],
        [0.0, 0.0, 0.0, 1.0],
    ]
}

const HALF: [f64; 3] = [0.5, 0.5, 0.5];

fn stacked(upper_y: f64, depth: f64) -> ContactManifold {
    build_box_manifold(
        HALF,
        pose([0.0, upper_y, 0.0]),
        HALF,
        pose([0.0; 3]),
        [0.0, 1.0, 0.0],
        depth,
    )
}

#[test]
fn stacked_unit_cubes_with_a_tenth_of_overlap_give_four_points_of_that_depth() {
    let manifold = stacked(0.9, 0.1);
    assert_eq!(manifold.normal, [0.0, 1.0, 0.0]);
    assert_eq!(manifold.contacts.len(), 4, "{:?}", manifold.contacts);
    for contact in &manifold.contacts {
        // Not the half extent (0.5) of the old constant-depth construction.
        assert!(
            (contact.depth - 0.1).abs() < 1e-12,
            "depth {}",
            contact.depth
        );
        assert_eq!(contact.normal, [0.0, 1.0, 0.0]);
        assert!(
            (contact.point[1] - 0.45).abs() < 1e-12,
            "{:?}",
            contact.point
        );
        assert!((contact.point[0].abs() - 0.5).abs() < 1e-12);
        assert!((contact.point[2].abs() - 0.5).abs() < 1e-12);
    }
}

#[test]
fn the_function_is_reachable_through_every_public_path() {
    let direct = stacked(0.9, 0.1);
    let through_module = via_manifold_module(
        HALF,
        pose([0.0, 0.9, 0.0]),
        HALF,
        pose([0.0; 3]),
        [0.0, 1.0, 0.0],
        0.1,
    );
    assert_eq!(direct.contacts.len(), through_module.contacts.len());
    for (a, b) in direct.contacts.iter().zip(&through_module.contacts) {
        assert_eq!(a.point, b.point);
        assert_eq!(a.depth, b.depth);
    }
}

#[test]
fn hovering_cubes_have_no_contacts_and_touching_cubes_have_zero_depth() {
    assert!(stacked(1.1, 0.0).contacts.is_empty());
    let touching = stacked(1.0, 0.0);
    assert_eq!(touching.contacts.len(), 4);
    assert!(touching.contacts.iter().all(|c| c.depth == 0.0));
}

#[test]
fn a_turned_cube_gives_an_octagon_that_reduces_to_four_points() {
    let mut manifold = build_box_manifold(
        HALF,
        yawed([0.0, 0.9, 0.0], std::f64::consts::FRAC_PI_4),
        HALF,
        pose([0.0; 3]),
        [0.0, 1.0, 0.0],
        0.1,
    );
    assert_eq!(manifold.contacts.len(), 8, "{:?}", manifold.contacts);
    for contact in &manifold.contacts {
        assert!((contact.depth - 0.1).abs() < 1e-9);
    }
    reduce_manifold(&mut manifold, 4);
    assert_eq!(manifold.contacts.len(), 4);
}

#[test]
fn the_normal_convention_is_from_b_to_a() {
    // A above B and the normal pointing up gives the overlap patch between
    // B's top face and A's bottom face; flipping both the boxes and the normal
    // gives the same patch.
    let upper = build_box_manifold(
        HALF,
        pose([0.0, 0.9, 0.0]),
        HALF,
        pose([0.0; 3]),
        [0.0, 1.0, 0.0],
        0.1,
    );
    let lower = build_box_manifold(
        HALF,
        pose([0.0; 3]),
        HALF,
        pose([0.0, 0.9, 0.0]),
        [0.0, -1.0, 0.0],
        0.1,
    );
    assert_eq!(upper.contacts.len(), 4);
    assert_eq!(lower.contacts.len(), 4);
    for contact in &upper.contacts {
        assert!(
            lower
                .contacts
                .iter()
                .any(|other| (0..3).all(|i| (other.point[i] - contact.point[i]).abs() < 1e-12)),
            "{:?} missing from {:?}",
            contact.point,
            lower.contacts
        );
    }
    assert!(lower.contacts.iter().all(|c| c.normal == [0.0, -1.0, 0.0]));
}
