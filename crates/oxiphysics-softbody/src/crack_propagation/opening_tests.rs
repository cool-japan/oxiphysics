// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tests for the topological crack opening of [`SoftBodyMesh`].
//!
//! The previous crack propagation only split edges at their midpoints, which
//! left `a - mid - b` connected with the triangles still referencing `a - b`:
//! the mesh never separated, so no new fragment could ever appear.

use super::*;
use std::collections::HashSet;

/// `nx` x `ny` nodes at `(i h, j h, 0)` with unit masses, horizontal, vertical and
/// "/" diagonal edges and two counter-clockwise triangles per cell.
///
/// Cell `(i, j)` owns the triangles `2 (j (nx - 1) + i)` (lower right) and
/// `2 (j (nx - 1) + i) + 1` (upper left).
fn grid(nx: usize, ny: usize, h: f64) -> SoftBodyMesh {
    let mut mesh = SoftBodyMesh::new();
    for j in 0..ny {
        for i in 0..nx {
            mesh.add_node([i as f64 * h, j as f64 * h, 0.0], 1.0);
        }
    }
    let id = |i: usize, j: usize| j * nx + i;
    for j in 0..ny {
        for i in 0..nx {
            if i + 1 < nx {
                mesh.add_edge(id(i, j), id(i + 1, j));
            }
            if j + 1 < ny {
                mesh.add_edge(id(i, j), id(i, j + 1));
            }
            if i + 1 < nx && j + 1 < ny {
                mesh.add_edge(id(i, j), id(i + 1, j + 1));
                mesh.add_triangle(id(i, j), id(i + 1, j), id(i + 1, j + 1));
                mesh.add_triangle(id(i, j), id(i + 1, j + 1), id(i, j + 1));
            }
        }
    }
    mesh
}

fn total_mass(mesh: &SoftBodyMesh) -> f64 {
    mesh.nodes.iter().map(|n| n.mass).sum()
}

/// Fragment index of every node.
fn fragment_ids(mesh: &SoftBodyMesh) -> Vec<usize> {
    let mut ids = vec![usize::MAX; mesh.node_count()];
    for fragment in identify_fragments(mesh) {
        for &node in &fragment.node_indices {
            ids[node] = fragment.id;
        }
    }
    ids
}

/// Structural invariants of a conforming mesh: valid indices, every triangle edge
/// present exactly once in the edge list, all three corners in one fragment.
fn assert_conforming(mesh: &SoftBodyMesh) {
    let n = mesh.node_count();
    let mut listed: HashSet<(usize, usize)> = HashSet::new();
    for &[a, b] in &mesh.edges {
        assert!(a < n && b < n, "edge ({a}, {b}) out of range");
        assert!(
            listed.insert((a.min(b), a.max(b))),
            "edge ({a}, {b}) is listed twice"
        );
    }
    let fragments = fragment_ids(mesh);
    for (t, &[a, b, c]) in mesh.triangles.iter().enumerate() {
        assert!(a < n && b < n && c < n, "triangle {t} out of range");
        for (x, y) in [(a, b), (b, c), (c, a)] {
            assert!(
                listed.contains(&(x.min(y), x.max(y))),
                "edge ({x}, {y}) of triangle {t} is missing from the edge list"
            );
        }
        assert!(
            fragments[a] == fragments[b] && fragments[b] == fragments[c],
            "triangle {t} straddles two fragments"
        );
    }
}

#[test]
fn crack_across_the_strip_leaves_two_fragments() {
    let mut mesh = grid(5, 5, 1.0);
    assert_eq!(identify_fragments(&mesh).len(), 1);
    let nodes_before = mesh.node_count();
    // Vertical crack at x = 1.5 reaching beyond the top and bottom boundary.
    let severed = open_crack_segment(&mut mesh, [1.5, -1.0, 0.0], [1.5, 5.0, 0.0]);
    // 5 horizontal + 4 diagonal edges connect the columns x = 1 and x = 2.
    assert_eq!(severed, 9);
    assert_eq!(identify_fragments(&mesh).len(), 2);
    assert!(
        mesh.node_count() > nodes_before,
        "nodes along the crack split"
    );
    assert_conforming(&mesh);

    // The two fragments are the two sides of the crack.
    let ids = fragment_ids(&mesh);
    let left = ids[0];
    for (i, node) in mesh.nodes.iter().enumerate().take(25) {
        let x = node.position[0];
        if x <= 1.0 {
            assert_eq!(
                ids[i], left,
                "node {i} at x = {x} belongs to the left piece"
            );
        }
        if x >= 2.0 {
            assert_ne!(
                ids[i], left,
                "node {i} at x = {x} belongs to the right piece"
            );
        }
    }
}

#[test]
fn mass_positions_and_velocities_survive_the_cut() {
    let mut mesh = grid(5, 5, 0.5);
    for (i, node) in mesh.nodes.iter_mut().enumerate() {
        node.mass = 1.0 + 0.1 * i as f64;
        node.velocity = [i as f64, -(i as f64), 0.25];
    }
    let mass_before = total_mass(&mesh);
    let original: Vec<([f64; 3], [f64; 3], f64)> = mesh
        .nodes
        .iter()
        .map(|n| (n.position, n.velocity, n.mass))
        .collect();
    let opening = open_crack(&mut mesh, [0.75, -1.0, 0.0], [0.75, 3.0, 0.0], None);
    assert!(opening.severed_edges > 0);
    assert!(!opening.duplicated_nodes.is_empty());
    assert!(
        (total_mass(&mesh) - mass_before).abs() < 1e-12 * mass_before,
        "splitting nodes must not create or destroy mass"
    );
    for (i, &(position, velocity, mass)) in original.iter().enumerate() {
        assert_eq!(mesh.nodes[i].position, position, "node {i} moved");
        assert_eq!(
            mesh.nodes[i].velocity, velocity,
            "node {i} changed velocity"
        );
        assert!(mesh.nodes[i].mass <= mass + 1e-15, "node {i} gained mass");
        assert!(mesh.nodes[i].mass > 0.0, "node {i} lost all its mass");
    }
    for &(source, duplicate) in &opening.duplicated_nodes {
        assert!(duplicate >= original.len());
        assert_eq!(mesh.nodes[duplicate].position, original[source].0);
        assert_eq!(mesh.nodes[duplicate].velocity, original[source].1);
        assert!(mesh.nodes[duplicate].mass > 0.0);
        assert_eq!(mesh.nodes[duplicate].id, duplicate);
    }
}

#[test]
fn partial_crack_opens_a_slit_and_keeps_the_body_connected() {
    let mut mesh = grid(5, 5, 1.0);
    let nodes_before = mesh.node_count();
    // Horizontal crack at y = 2.5 from outside the left boundary to x = 2.2.
    let opening = open_crack(&mut mesh, [-1.0, 2.5, 0.0], [2.2, 2.5, 0.0], None);
    // Vertical edges at x = 0, 1, 2 and diagonals at x = 0.5, 1.5 are crossed.
    assert_eq!(opening.severed_edges, 5);
    assert_eq!(
        identify_fragments(&mesh).len(),
        1,
        "a crack that ends inside the body does not separate it"
    );
    assert!(mesh.node_count() > nodes_before);
    assert_eq!(
        mesh.node_count(),
        nodes_before + opening.duplicated_nodes.len()
    );
    assert_conforming(&mesh);

    // Along the slit the lower and the upper triangle of each cut cell share no
    // node any more; beyond the tip the cell is still joined.
    let cell_triangles = |i: usize| (2 * (2 * 4 + i), 2 * (2 * 4 + i) + 1);
    for i in 0..2 {
        let (lower, upper) = cell_triangles(i);
        for node in mesh.triangles[lower] {
            assert!(
                !mesh.triangles[upper].contains(&node),
                "cell {i}: the crack faces still share node {node}"
            );
        }
    }
    let (lower, upper) = cell_triangles(3);
    assert!(
        mesh.triangles[lower]
            .iter()
            .any(|n| mesh.triangles[upper].contains(n)),
        "cell 3 lies beyond the tip and must stay intact"
    );
}

#[test]
fn edges_without_triangles_are_removed_when_crossed() {
    // A chain along x with no triangles at all.
    let mut mesh = SoftBodyMesh::new();
    for i in 0..5 {
        mesh.add_node([i as f64, 0.0, 0.0], 1.0);
    }
    for i in 0..4 {
        mesh.add_edge(i, i + 1);
    }
    let severed = open_crack_segment(&mut mesh, [2.5, -1.0, 0.0], [2.5, 1.0, 0.0]);
    assert_eq!(severed, 1);
    assert_eq!(mesh.edges, vec![[0, 1], [1, 2], [3, 4]]);
    assert_eq!(mesh.node_count(), 5, "graph meshes need no duplicate nodes");
    assert_eq!(identify_fragments(&mesh).len(), 2);
}

#[test]
fn opening_the_same_segment_twice_changes_nothing_more() {
    let mut mesh = grid(5, 5, 1.0);
    open_crack_segment(&mut mesh, [-1.0, 2.5, 0.0], [2.2, 2.5, 0.0]);
    let snapshot = mesh.clone();
    let again = open_crack(&mut mesh, [-1.0, 2.5, 0.0], [2.2, 2.5, 0.0], None);
    assert_eq!(again, CrackOpening::default());
    assert_eq!(mesh.edges, snapshot.edges);
    assert_eq!(mesh.triangles, snapshot.triangles);
    assert_eq!(mesh.node_count(), snapshot.node_count());
}

#[test]
fn a_crack_in_two_steps_opens_like_one_long_crack() {
    let mut once = grid(5, 5, 1.0);
    open_crack_segment(&mut once, [-1.0, 2.5, 0.0], [6.0, 2.5, 0.0]);
    assert_eq!(identify_fragments(&once).len(), 2);
    assert_conforming(&once);

    let mut twice = grid(5, 5, 1.0);
    open_crack_segment(&mut twice, [-1.0, 2.5, 0.0], [2.0, 2.5, 0.0]);
    assert_eq!(identify_fragments(&twice).len(), 1);
    open_crack_segment(&mut twice, [2.0, 2.5, 0.0], [6.0, 2.5, 0.0]);
    assert_eq!(identify_fragments(&twice).len(), 2);
    assert_conforming(&twice);
    assert_eq!(twice.node_count(), once.node_count());
    assert_eq!(twice.edges.len(), once.edges.len());
}

#[test]
fn direction_of_travel_does_not_matter() {
    for (from, to) in [
        ([1.5, -1.0, 0.0], [1.5, 5.0, 0.0]),
        ([1.5, 5.0, 0.0], [1.5, -1.0, 0.0]),
    ] {
        let mut mesh = grid(5, 5, 1.0);
        let severed = open_crack_segment(&mut mesh, from, to);
        assert_eq!(severed, 9);
        assert_eq!(identify_fragments(&mesh).len(), 2);
        assert_conforming(&mesh);
    }
}

#[test]
fn diagonal_crack_cuts_off_a_corner() {
    let mut mesh = grid(5, 5, 1.0);
    // The line x + y = 1.5 separates (0,0), (1,0), (0,1) from the rest.  Both
    // triangles of the corner cell lie on the near side although the crack crosses
    // their shared diagonal: they must stay joined to each other.
    let severed = open_crack_segment(&mut mesh, [-1.0, 2.5, 0.0], [2.5, -1.0, 0.0]);
    assert!(severed > 0);
    assert_eq!(identify_fragments(&mesh).len(), 2);
    assert_conforming(&mesh);
    let ids = fragment_ids(&mesh);
    let corner = ids[mesh.triangles[0][0]];
    for (t, tri) in mesh.triangles.iter().enumerate().take(2) {
        for &node in tri {
            assert_eq!(
                ids[node], corner,
                "corner triangle {t} left the corner piece"
            );
        }
    }
    let shared = mesh.triangles[0]
        .iter()
        .filter(|n| mesh.triangles[1].contains(n))
        .count();
    assert!(
        shared >= 2,
        "the corner triangles keep sharing their diagonal"
    );
    for (t, tri) in mesh.triangles.iter().enumerate().skip(2) {
        for &node in tri {
            assert_ne!(ids[node], corner, "triangle {t} must be on the far side");
        }
    }
}

#[test]
fn cracks_that_miss_the_mesh_or_are_degenerate_leave_it_untouched() {
    let reference = grid(4, 4, 1.0);
    let attempts: [([f64; 3], [f64; 3], Option<[f64; 3]>); 5] = [
        // Entirely outside the body.
        ([10.0, 0.0, 0.0], [10.0, 5.0, 0.0], None),
        // Zero length.
        ([1.5, 1.5, 0.0], [1.5, 1.5, 0.0], None),
        // Not finite.
        ([f64::NAN, 0.0, 0.0], [1.0, 1.0, 0.0], None),
        // A crack running along the surface normal cannot cut a surface.
        ([1.5, 1.5, -1.0], [1.5, 1.5, 1.0], Some([0.0, 0.0, 1.0])),
        // Ends before reaching the first edge.
        ([1.5, 1.4, 0.0], [1.5, 1.45, 0.0], None),
    ];
    for (from, to, normal) in attempts {
        let mut mesh = reference.clone();
        let opening = open_crack(&mut mesh, from, to, normal);
        assert_eq!(opening, CrackOpening::default(), "{from:?} -> {to:?}");
        assert_eq!(mesh.edges, reference.edges);
        assert_eq!(mesh.triangles, reference.triangles);
        assert_eq!(mesh.node_count(), reference.node_count());
    }
}

#[test]
fn explicit_surface_normal_selects_the_cutting_plane() {
    // Chain along x at y = 0, z = 0.
    let build = || {
        let mut mesh = SoftBodyMesh::new();
        for i in 0..4 {
            mesh.add_node([i as f64, 0.0, 0.0], 1.0);
        }
        for i in 0..3 {
            mesh.add_edge(i, i + 1);
        }
        mesh
    };
    // A crack along z through x = 1.5: with the default (y) surface normal it cuts
    // the plane x = 1.5 and severs edge (1, 2).
    let mut cut = build();
    assert_eq!(
        open_crack_segment(&mut cut, [1.5, 0.0, -1.0], [1.5, 0.0, 1.0]),
        1
    );
    // With the surface normal along x the cutting plane is y = 0, which contains
    // the whole chain: nothing is crossed.
    let mut grazing = build();
    let opening = open_crack(
        &mut grazing,
        [1.5, 0.0, -1.0],
        [1.5, 0.0, 1.0],
        Some([1.0, 0.0, 0.0]),
    );
    assert_eq!(opening, CrackOpening::default());
    assert_eq!(grazing.edges.len(), 3);
}

#[test]
fn propagating_crack_cuts_a_triangle_mesh() {
    let mut mesh = grid(5, 5, 1.0);
    // Starts below the strip and travels 9.75 m in +y at x = 1.5.
    let mut crack = CrackPath::new([1.5, -1.0, 0.0], [0.0, 1.0, 0.0]);
    let propagated = propagate_crack(&mut mesh, &mut crack, &CrackCriteria::glass(), 2.0, 0.01);
    assert!(propagated);
    assert_eq!(identify_fragments(&mesh).len(), 2);
    assert_conforming(&mesh);
    // Below the toughness the mesh stays whole.
    let mut intact = grid(5, 5, 1.0);
    let mut other = CrackPath::new([1.5, -1.0, 0.0], [0.0, 1.0, 0.0]);
    assert!(!propagate_crack(
        &mut intact,
        &mut other,
        &CrackCriteria::glass(),
        0.3,
        0.01
    ));
    assert_eq!(identify_fragments(&intact).len(), 1);
}

#[test]
fn dynamic_fracture_step_cuts_the_mesh_and_tracks_energy() {
    let mut mesh = grid(5, 5, 1.0);
    let mut crack = CrackPath::new([1.5, -1.0, 0.0], [0.0, 1.0, 0.0]);
    let mut fracture = DynamicFracture::new(CrackCriteria::glass(), 1.0);
    // G = K^2 / E exceeds G_c = 1: K = 3e5 Pa sqrt(m) at E = 70 GPa.
    let advance = fracture.step(&mut crack, &mut mesh, 3.0e5, 70.0e9, 0.001);
    assert!(advance > 0.0);
    assert!(fracture.total_energy_dissipated() > 0.0);
    // The tip travelled 1.5 km/s * 1 ms = 1.5 m: it enters the strip at y = 0
    // and stops inside it, so the body is slit but not separated.
    let tip = crack
        .current_tip()
        .map(|t| t.position[1])
        .unwrap_or(f64::NAN);
    assert!(tip > 0.0 && tip < 4.0, "tip at y = {tip}");
    assert_eq!(identify_fragments(&mesh).len(), 1);
    assert!(
        mesh.node_count() > 25,
        "the slit must have duplicated nodes"
    );
    assert_conforming(&mesh);
}

#[test]
fn repeated_cracks_keep_a_larger_mesh_consistent() {
    let mut mesh = grid(9, 9, 1.0);
    let mass = total_mass(&mesh);
    let cracks = [
        ([-1.0, 3.3, 0.0], [4.2, 3.3, 0.0]),
        ([2.5, -1.0, 0.0], [2.5, 6.1, 0.0]),
        ([9.0, 7.7, 0.0], [3.1, 4.9, 0.0]),
        ([6.4, 9.0, 0.0], [6.4, 1.2, 0.0]),
    ];
    let mut previous_fragments = 1;
    for (from, to) in cracks {
        open_crack_segment(&mut mesh, from, to);
        assert_conforming(&mesh);
        assert!((total_mass(&mesh) - mass).abs() < 1e-9 * mass);
        let fragments = identify_fragments(&mesh).len();
        assert!(fragments >= previous_fragments, "cracks never heal");
        previous_fragments = fragments;
    }
    assert!(previous_fragments >= 2, "the cracks must separate the body");
}
