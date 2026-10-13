// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the empty-mesh BVH and the grid-accelerated cloth
//! self-collision.
//!
//! * `TriMeshBvh::build` on an empty triangle list recursed on the same empty
//!   slice until the stack overflowed (an uncatchable abort).
//! * `cloth_self_collision` scanned the whole triangle list for every vertex
//!   pair (`O(V^2 T)`); it now uses a spatial hash and must keep the exact
//!   behaviour of the brute-force scan, position corrections included.

use super::*;

// ── empty meshes (F0150) ─────────────────────────────────────────────────────

fn unit_triangle_at(x: f64) -> Triangle {
    Triangle::new([x, 0.0, 0.0], [x + 1.0, 0.0, 0.0], [x, 1.0, 0.0])
}

#[test]
fn empty_bvh_has_no_nodes_and_the_sentinel_root() {
    let bvh = TriMeshBvh::build(Vec::new());
    assert_eq!(bvh.num_triangles(), 0);
    assert_eq!(bvh.num_nodes(), 0);
    assert_eq!(bvh.root, usize::MAX);
}

#[test]
fn building_a_node_range_without_indices_returns_the_sentinel() {
    // The recursion guard on its own: an empty index range has no node.
    let tris = vec![unit_triangle_at(0.0)];
    let mut nodes = Vec::new();
    let root = TriMeshBvh::build_recursive(&tris, &[], &mut nodes);
    assert_eq!(root, usize::MAX);
    assert!(nodes.is_empty());
}

#[test]
fn queries_on_an_empty_bvh_find_nothing() {
    let mut bvh = TriMeshBvh::build(Vec::new());
    assert!(
        bvh.ray_cast([0.0, 0.0, 5.0], [0.0, 0.0, -1.0], 100.0)
            .is_none()
    );
    assert!(
        bvh.query_aabb(&MeshAabb::new([-10.0; 3], [10.0; 3]))
            .is_empty()
    );
    bvh.refit();
    assert_eq!(bvh.num_nodes(), 0);
}

#[test]
fn mesh_broad_phase_with_an_empty_side_is_empty() {
    let empty = TriMeshBvh::build(Vec::new());
    let full = TriMeshBvh::build(vec![unit_triangle_at(0.0), unit_triangle_at(2.0)]);
    assert!(mesh_broad_phase(&empty, &full).is_empty());
    assert!(mesh_broad_phase(&full, &empty).is_empty());
    assert!(mesh_broad_phase(&empty, &TriMeshBvh::build(Vec::new())).is_empty());
    assert!(mesh_mesh_collision(&empty, &full).is_empty());
    assert!(mesh_mesh_collision(&full, &empty).is_empty());
}

#[test]
fn mesh_level_entry_points_accept_an_empty_bvh() {
    let empty = TriMeshBvh::build(Vec::new());
    let full = TriMeshBvh::build(vec![unit_triangle_at(0.0)]);
    assert_eq!(penetration_depth_estimate(&empty, &full, 4), 0.0);
    assert_eq!(penetration_depth_estimate(&full, &empty, 4), 0.0);
    assert_eq!(
        generate_mesh_contact_manifold(&empty, &full).num_contacts(),
        0
    );
    let swept = swept_triangle_vs_mesh(&unit_triangle_at(0.0), [0.0, 0.0, 1.0], &empty);
    assert!(swept.is_none());
    assert!(mesh_narrow_phase(&empty, &full, &[]).is_empty());
}

#[test]
fn deformable_mesh_without_faces_builds_and_updates() {
    let verts = vec![[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    let mut dm = DeformableMesh::build(verts, Vec::new());
    assert_eq!(dm.num_faces(), 0);
    assert_eq!(dm.bvh.num_nodes(), 0);
    dm.update_vertex(1, [4.0, 0.0, 0.0]);
    dm.apply_displacement(&[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]);
    assert_eq!(dm.vertices[1], [4.0, 1.0, 0.0]);
    assert!(
        dm.bvh
            .ray_cast([0.0, 0.0, 5.0], [0.0, 0.0, -1.0], 100.0)
            .is_none()
    );
}

#[test]
fn non_empty_bvh_is_unchanged_by_the_empty_guard() {
    let tris: Vec<Triangle> = (0..9).map(|i| unit_triangle_at(i as f64 * 2.0)).collect();
    let bvh = TriMeshBvh::build(tris);
    assert_eq!(bvh.num_triangles(), 9);
    // A full binary tree over 9 leaves has 17 nodes.
    assert_eq!(bvh.num_nodes(), 17);
    assert!(bvh.root < bvh.num_nodes());
    let hit = bvh.ray_cast([6.2, 0.2, 3.0], [0.0, 0.0, -1.0], 10.0);
    assert_eq!(hit.map(|(_, tri)| tri), Some(3));
}

// ── cloth self-collision (F0151) ─────────────────────────────────────────────

/// The original brute-force algorithm, kept verbatim as the reference.
fn cloth_self_collision_reference(
    vertices: &mut [ClothVertex],
    triangles: &[[usize; 3]],
    thickness: f64,
) -> usize {
    let mut n_collisions = 0;
    let n = vertices.len();
    for i in 0..n {
        for j in (i + 1)..n {
            let shared = triangles.iter().any(|t| {
                (t[0] == i || t[1] == i || t[2] == i) && (t[0] == j || t[1] == j || t[2] == j)
            });
            if shared {
                continue;
            }
            let diff = vec3_sub(vertices[i].pos, vertices[j].pos);
            let dist = vec3_len(diff);
            if dist < thickness && dist > 1e-15 {
                let push = vec3_scale(vec3_norm(diff), (thickness - dist) * 0.5);
                if !vertices[i].pinned {
                    vertices[i].pos = vec3_add(vertices[i].pos, push);
                }
                if !vertices[j].pinned {
                    vertices[j].pos = vec3_sub(vertices[j].pos, push);
                }
                n_collisions += 1;
            }
        }
    }
    n_collisions
}

/// Deterministic xorshift generator.
struct Rng(u64);

impl Rng {
    fn unit(&mut self) -> f64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }

    fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.unit()
    }

    fn index(&mut self, n: usize) -> usize {
        ((self.unit() * n as f64) as usize).min(n.saturating_sub(1))
    }
}

fn cloth_from_positions(positions: &[[f64; 3]]) -> Vec<ClothVertex> {
    positions
        .iter()
        .map(|&p| ClothVertex::new(p, 1.0))
        .collect()
}

fn assert_same_cloth(a: &[ClothVertex], b: &[ClothVertex]) {
    assert_eq!(a.len(), b.len());
    for (i, (va, vb)) in a.iter().zip(b.iter()).enumerate() {
        assert_eq!(va.pos, vb.pos, "vertex {i} moved differently");
    }
}

#[test]
fn grid_collision_matches_brute_force_on_random_clouds() {
    let mut rng = Rng(0xC107_45E1_FC01_11E5);
    let mut total_collisions = 0;
    for case in 0..250 {
        let n = 2 + rng.index(70);
        let extent = rng.range(0.5, 2.0);
        let thickness = rng.range(0.02, 0.6) * extent;
        let positions: Vec<[f64; 3]> = (0..n)
            .map(|_| {
                [
                    rng.range(0.0, extent),
                    rng.range(0.0, extent),
                    rng.range(0.0, extent),
                ]
            })
            .collect();
        // Random triangles: shared vertex pairs are skipped by the collision pass.
        let tri_count = rng.index(n + 1);
        let triangles: Vec<[usize; 3]> = (0..tri_count)
            .map(|_| [rng.index(n), rng.index(n), rng.index(n)])
            .collect();

        let mut expected = cloth_from_positions(&positions);
        let mut actual = cloth_from_positions(&positions);
        for (e, a) in expected.iter_mut().zip(actual.iter_mut()) {
            if rng.unit() < 0.15 {
                e.pinned = true;
                a.pinned = true;
            }
        }

        let n_expected = cloth_self_collision_reference(&mut expected, &triangles, thickness);
        let n_actual = cloth_self_collision(&mut actual, &triangles, thickness);
        assert_eq!(n_actual, n_expected, "case {case}: collision count");
        assert_same_cloth(&actual, &expected);
        total_collisions += n_expected;
    }
    // The clouds are dense enough that corrections cascade through many pairs.
    assert!(
        total_collisions > 1500,
        "only {total_collisions} collisions exercised"
    );
}

#[test]
fn grid_collision_matches_brute_force_with_duplicates_and_bad_triangles() {
    // Coincident vertices (distance 0 is skipped), out-of-range triangle corners,
    // repeated triangle corners and a NaN position.
    let positions = [
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.01, 0.0, 0.0],
        [0.02, 0.01, 0.0],
        [f64::NAN, 0.0, 0.0],
        [0.03, 0.0, 0.01],
        [5.0, 5.0, 5.0],
    ];
    let triangles = [[0, 2, 99], [3, 3, 5], [100, 101, 102], [1, 1, 1]];
    let mut expected = cloth_from_positions(&positions);
    let mut actual = cloth_from_positions(&positions);
    let n_expected = cloth_self_collision_reference(&mut expected, &triangles, 0.05);
    let n_actual = cloth_self_collision(&mut actual, &triangles, 0.05);
    assert_eq!(n_actual, n_expected);
    assert!(n_actual > 0);
    for (a, e) in actual.iter().zip(expected.iter()) {
        // NaN != NaN, so compare bit patterns.
        for (x, y) in a.pos.iter().zip(e.pos.iter()) {
            assert_eq!(x.to_bits(), y.to_bits());
        }
    }
}

#[test]
fn grid_collision_matches_brute_force_far_from_the_origin() {
    let mut rng = Rng(0xFA12_0123_4567_89AB);
    let offset = [1.0e6, -3.0e6, 2.0e5];
    let positions: Vec<[f64; 3]> = (0..60)
        .map(|_| {
            [
                offset[0] + rng.range(0.0, 1.0),
                offset[1] + rng.range(0.0, 1.0),
                offset[2] + rng.range(0.0, 1.0),
            ]
        })
        .collect();
    let mut expected = cloth_from_positions(&positions);
    let mut actual = cloth_from_positions(&positions);
    let n_expected = cloth_self_collision_reference(&mut expected, &[], 0.2);
    let n_actual = cloth_self_collision(&mut actual, &[], 0.2);
    assert_eq!(n_actual, n_expected);
    assert!(n_actual > 0);
    assert_same_cloth(&actual, &expected);
}

#[test]
fn isolated_close_pairs_on_a_lattice_are_counted_once_each() {
    // A 12 x 12 x 1 lattice with spacing 1 and thickness 0.1, plus three
    // isolated vertices that each sit 0.02 from one lattice vertex.
    let mut positions = Vec::new();
    for ix in 0..12 {
        for iy in 0..12 {
            positions.push([ix as f64, iy as f64, 0.0]);
        }
    }
    let lattice = positions.len();
    positions.push([3.02, 4.0, 0.0]);
    positions.push([7.0, 8.0, 0.02]);
    positions.push([10.0, 1.02, 0.0]);
    let mut verts = cloth_from_positions(&positions);
    let n = cloth_self_collision(&mut verts, &[], 0.1);
    assert_eq!(n, 3);
    // Only the six vertices of the three pairs moved, apart by the thickness.
    let moved: Vec<usize> = (0..positions.len())
        .filter(|&i| verts[i].pos != positions[i])
        .collect();
    assert_eq!(moved.len(), 6, "moved vertices: {moved:?}");
    for (lattice_vertex, extra) in [
        (3 * 12 + 4, lattice),
        (7 * 12 + 8, lattice + 1),
        (10 * 12 + 1, lattice + 2),
    ] {
        let gap = vec3_len(vec3_sub(verts[lattice_vertex].pos, verts[extra].pos));
        assert!((gap - 0.1).abs() < 1e-12, "pair separated to {gap}");
    }
}

#[test]
fn vertices_sharing_a_triangle_never_collide() {
    let positions = [
        [0.0, 0.0, 0.0],
        [0.001, 0.0, 0.0],
        [0.0, 0.001, 0.0],
        [0.5, 0.0, 0.0],
    ];
    let triangles = [[0, 1, 2]];
    let mut verts = cloth_from_positions(&positions);
    let n = cloth_self_collision(&mut verts, &triangles, 0.05);
    assert_eq!(n, 0);
    assert_same_cloth(&verts, &cloth_from_positions(&positions));
}

#[test]
fn pinned_vertices_do_not_move_and_the_free_one_takes_half_the_correction() {
    let mut verts = cloth_from_positions(&[[0.0, 0.0, 0.0], [0.02, 0.0, 0.0]]);
    verts[0].pinned = true;
    let n = cloth_self_collision(&mut verts, &[], 0.1);
    assert_eq!(n, 1);
    assert_eq!(verts[0].pos, [0.0, 0.0, 0.0]);
    // push = (0.1 - 0.02) / 2 along -x for the second vertex.
    assert!((verts[1].pos[0] - 0.06).abs() < 1e-15, "{:?}", verts[1].pos);
}

#[test]
fn degenerate_thickness_produces_no_collisions() {
    let positions = [[0.0; 3], [0.001, 0.0, 0.0], [0.002, 0.0, 0.0]];
    for thickness in [0.0, -1.0, 1e-16, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let mut verts = cloth_from_positions(&positions);
        assert_eq!(
            cloth_self_collision(&mut verts, &[], thickness),
            0,
            "thickness {thickness}"
        );
        assert_same_cloth(&verts, &cloth_from_positions(&positions));
    }
    let mut empty: Vec<ClothVertex> = Vec::new();
    assert_eq!(cloth_self_collision(&mut empty, &[], 0.1), 0);
    let mut single = cloth_from_positions(&[[1.0, 2.0, 3.0]]);
    assert_eq!(cloth_self_collision(&mut single, &[], 0.1), 0);
}

#[test]
fn ten_thousand_vertex_cloth_finishes_and_finds_the_injected_pairs() {
    // 100 x 100 sheet, spacing 0.01, structural triangles; the O(V^2 T) scan
    // needed about 1e12 operations for this input.
    let side = 100;
    let spacing = 0.01;
    let mut positions = Vec::with_capacity(side * side);
    for iy in 0..side {
        for ix in 0..side {
            positions.push([ix as f64 * spacing, iy as f64 * spacing, 0.0]);
        }
    }
    let mut triangles = Vec::new();
    for iy in 0..side - 1 {
        for ix in 0..side - 1 {
            let a = iy * side + ix;
            let (b, c, d) = (a + 1, a + side, a + side + 1);
            triangles.push([a, b, c]);
            triangles.push([b, d, c]);
        }
    }
    // Three extra vertices hover 0.003 above sheet vertices, i.e. inside the
    // collision thickness of exactly one non-adjacent vertex each (the sheet
    // neighbours are sqrt(0.01^2 + 0.003^2) = 0.0104 away).
    let thickness = 0.005;
    let mut verts = cloth_from_positions(&positions);
    let mut pairs = Vec::new();
    for sheet in [10 * side + 10, 40 * side + 60, 80 * side + 25] {
        let mut partner = positions[sheet];
        partner[2] = 0.003;
        verts.push(ClothVertex::new(partner, 1.0));
        pairs.push((sheet, verts.len() - 1));
    }
    let originals = cloth_from_positions(&positions);
    let n = cloth_self_collision(&mut verts, &triangles, thickness);
    assert_eq!(n, pairs.len());
    for &(sheet, partner) in &pairs {
        let gap = vec3_len(vec3_sub(verts[sheet].pos, verts[partner].pos));
        assert!((gap - thickness).abs() < 1e-12, "pair separated to {gap}");
    }
    // Nothing else moved.
    let moved = originals
        .iter()
        .zip(verts.iter())
        .filter(|(before, after)| before.pos != after.pos)
        .count();
    assert_eq!(moved, pairs.len());
}
