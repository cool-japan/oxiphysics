// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Black-box regression tests for the broad-phase pair queries, the box-box and
//! persistent contact manifolds, convex decomposition and the GJK time of
//! impact.
//!
//! Everything goes through the public paths of the crate, so the tests also
//! guard the module layout (`contact_generation::BoxBoxManifold`,
//! `compound_shapes::ConvexDecomposition`, ...).

use oxiphysics_collision::compound_shapes::{
    CollisionMeshMode, ConvexDecomposition, ImportedCollisionMesh, PrimitiveShape,
};
use oxiphysics_collision::contact_generation::{BoxBoxManifold, ContactPoint, PersistentManifold};
use oxiphysics_collision::gjk_enhanced::{EnhBox, EnhSphere, gjk_time_of_impact};
use oxiphysics_collision::{
    BroadPhase, BruteForceBroadPhase, BvhAabb, BvhBroadphase, CollisionPair, DynamicBvh,
};
use oxiphysics_core::{Aabb, Vec3};

const IDENTITY: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

fn cube(centre: [f64; 3], half: f64) -> Aabb {
    Aabb::new(
        Vec3::new(centre[0] - half, centre[1] - half, centre[2] - half),
        Vec3::new(centre[0] + half, centre[1] + half, centre[2] + half),
    )
}

/// A 6 x 6 x 6 lattice of boxes with side 1 spaced 0.9 apart: every box
/// overlaps its 26 neighbours and no one else.
fn lattice() -> Vec<Aabb> {
    let mut boxes = Vec::new();
    for i in 0..6 {
        for j in 0..6 {
            for k in 0..6 {
                boxes.push(cube([i as f64 * 0.9, j as f64 * 0.9, k as f64 * 0.9], 0.5));
            }
        }
    }
    boxes
}

#[test]
fn unbuilt_bvh_broadphase_finds_the_overlap() {
    let aabbs = [
        cube([0.0, 0.0, 0.0], 1.0),
        cube([1.5, 0.0, 0.0], 1.0),
        cube([9.0, 9.0, 9.0], 1.0),
    ];
    let pairs = BvhBroadphase::default().find_pairs(&aabbs);
    assert_eq!(pairs, vec![CollisionPair::new(0, 1)]);
}

#[test]
fn bvh_broadphase_matches_brute_force_on_a_lattice() {
    let aabbs = lattice();
    let expected = BruteForceBroadPhase.find_pairs(&aabbs);
    assert!(expected.len() > 1000, "only {} pairs", expected.len());
    let mut bvh = BvhBroadphase::new();
    bvh.build(&aabbs);
    assert_eq!(bvh.find_pairs(&aabbs), expected);
    // A hierarchy that was never built answers the same query.
    assert_eq!(BvhBroadphase::new().find_pairs(&aabbs), expected);
}

#[test]
fn dynamic_bvh_pairs_match_brute_force() {
    // Boxes of side 2 on a 1.5 grid: neighbours (including diagonal ones)
    // overlap by 0.5, second neighbours are 1.0 apart, so the result does not
    // depend on the fattening of the stored leaf boxes.
    let mut tree = DynamicBvh::new();
    let mut boxes = Vec::new();
    for data in 0..60u32 {
        let x = f64::from(data % 6) * 1.5;
        let y = f64::from((data / 6) % 5) * 1.5;
        let z = f64::from(data / 30) * 1.5;
        let aabb = BvhAabb::new(Vec3::new(x, y, z), Vec3::new(x + 2.0, y + 2.0, z + 2.0));
        tree.insert(aabb, data);
        boxes.push(aabb);
    }
    let mut expected = Vec::new();
    for (i, a) in boxes.iter().enumerate() {
        for (j, b) in boxes.iter().enumerate().skip(i + 1) {
            if a.intersects(b) {
                expected.push((i as u32, j as u32));
            }
        }
    }
    assert!(expected.len() > 100, "only {} pairs", expected.len());
    assert_eq!(tree.query_overlapping_pairs(), expected);
}

#[test]
fn box_box_manifold_reports_the_least_penetrating_face() {
    let manifold = BoxBoxManifold::new(0.01);
    let half = [1.0; 3];

    // A gap of 1 along x: no contact (the one-sided test reported four).
    let separated = manifold.compute([0.0; 3], half, IDENTITY, [-3.0, 0.0, 0.0], half, IDENTITY);
    assert!(separated.is_empty(), "got {separated:?}");

    // An overlap of 0.5 along x: four contacts, normal from A toward B.
    let contacts = manifold.compute([0.0; 3], half, IDENTITY, [1.5, 0.0, 0.0], half, IDENTITY);
    assert_eq!(contacts.len(), 4, "got {contacts:?}");
    for contact in &contacts {
        assert!((contact.normal[0] - 1.0).abs() < 1e-9, "{contact:?}");
        assert!(contact.normal[1].abs() < 1e-9 && contact.normal[2].abs() < 1e-9);
        assert!((contact.depth - 0.5).abs() < 1e-9, "{contact:?}");
    }
}

#[test]
fn persistent_manifold_drops_ghost_contacts() {
    let mut manifold = PersistentManifold::new(1, 2, 0.1);
    let contact = ContactPoint::new([0.0; 3], [0.0; 3], [0.0, 0.0, 1.0], 0.1);
    manifold.update(&[contact], 1);
    assert_eq!(manifold.num_points(), 1);

    // The bodies separate: the point survives `max_lifetime` updates and no more.
    let lifetime = manifold.max_lifetime as u64;
    for frame in 2..=(1 + lifetime) {
        manifold.update(&[], frame);
        assert_eq!(manifold.num_points(), 1, "dropped early at frame {frame}");
    }
    manifold.update(&[], 2 + lifetime);
    assert_eq!(manifold.num_points(), 0, "ghost contact never expired");
}

/// Prism over the L-shaped polygon (0,0) (2,0) (2,1) (1,1) (1,2) (0,2), from
/// `z = 0` to `z = 1`, with outward-facing triangles. Volume 3; the volume of
/// its convex hull is 3.5.
fn l_prism() -> (Vec<[f64; 3]>, Vec<[usize; 3]>) {
    let polygon = [
        [0.0, 0.0],
        [2.0, 0.0],
        [2.0, 1.0],
        [1.0, 1.0],
        [1.0, 2.0],
        [0.0, 2.0],
    ];
    let n = polygon.len();
    let mut vertices: Vec<[f64; 3]> = polygon.iter().map(|p| [p[0], p[1], 0.0]).collect();
    vertices.extend(polygon.iter().map(|p| [p[0], p[1], 1.0]));
    let mut indices = Vec::new();
    // A fan from the first corner triangulates this particular polygon.
    for k in 1..n - 1 {
        indices.push([n, n + k, n + k + 1]);
        indices.push([0, k + 1, k]);
    }
    for i in 0..n {
        let j = (i + 1) % n;
        indices.push([i, j, n + j]);
        indices.push([i, n + j, n + i]);
    }
    (vertices, indices)
}

#[test]
fn l_shaped_prism_decomposes_into_two_convex_parts() {
    let (vertices, indices) = l_prism();
    let mut decomposition = ConvexDecomposition::new(8, 0.05, 4);
    decomposition.decompose(&vertices, &indices);
    assert_eq!(decomposition.num_hulls(), 2, "{:?}", decomposition.hulls);
    let volume: f64 = decomposition
        .hulls
        .iter()
        .map(|hull| PrimitiveShape::ConvexHull(hull.clone()).volume())
        .sum();
    assert!((volume - 3.0).abs() < 1e-6, "total hull volume {volume}");
    // Without the split the single hull would be 3.5.
    let whole = PrimitiveShape::ConvexHull(vertices.clone()).volume();
    assert!((whole - 3.5).abs() < 1e-9, "hull volume {whole}");

    let mesh =
        ImportedCollisionMesh::new(vertices, indices, CollisionMeshMode::ConvexDecomposition);
    assert_eq!(mesh.decomposed_hulls.len(), 2);
}

#[test]
fn imported_mesh_hull_keeps_every_cube_corner() {
    let mut vertices = Vec::new();
    for sx in [-1.0, 1.0] {
        for sy in [-1.0, 1.0] {
            for sz in [-1.0, 1.0] {
                vertices.push([sx, sy, sz]);
            }
        }
    }
    vertices.push([0.0, 0.0, 0.0]);
    let mesh = ImportedCollisionMesh::new(vertices, Vec::new(), CollisionMeshMode::ConvexHull);
    // Six directional extremes would keep at most six of the eight corners.
    assert_eq!(mesh.convex_hull.len(), 8);
    assert_eq!(mesh.support([1.0, 1.0, 1.0]), [1.0, 1.0, 1.0]);
    let volume = PrimitiveShape::ConvexHull(mesh.convex_hull.clone()).volume();
    assert!((volume - 8.0).abs() < 1e-9, "hull volume {volume}");
}

#[test]
fn fast_sphere_does_not_tunnel_through_a_box() {
    // Separated at both ends of the sweep, but the path crosses the box.
    let sphere = EnhSphere::new([-5.0, 0.0, 0.0], 0.1);
    let wall = EnhBox::new([0.0; 3], [0.5; 3]);
    let result = gjk_time_of_impact(&sphere, &wall, [10.0, 0.0, 0.0], 64, 1e-3);
    assert!(result.hit, "tunnelled through the box: {result:?}");
    assert!(result.toi <= 0.44 + 1e-9, "toi {} is late", result.toi);
    assert!((result.toi - 0.44).abs() < 2e-3, "toi {}", result.toi);
}
