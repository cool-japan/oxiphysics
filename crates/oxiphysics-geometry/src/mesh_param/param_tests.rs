// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the Tutte, harmonic-map and LSCM parameterizations.

use super::lscm::{NormalEquations, assemble_normal_equations, conformal_row, farthest_pair_pins};
use super::*;
use crate::error::Error;
use std::f64::consts::{PI, TAU};

// ── Test meshes ──────────────────────────────────────────────────────────────

/// Planar `n x n` vertex grid with spacing `h`, counter-clockwise in the xy plane;
/// vertex `(i, j)` has id `j * n + i`.
fn grid_mesh(n: usize, h: f64) -> ParamTriMesh {
    let mut mesh = ParamTriMesh::new();
    for j in 0..n {
        for i in 0..n {
            mesh.add_vertex([i as f64 * h, j as f64 * h, 0.0], [0.0, 0.0]);
        }
    }
    let id = |i: usize, j: usize| j * n + i;
    for j in 0..n - 1 {
        for i in 0..n - 1 {
            mesh.add_triangle(id(i, j), id(i + 1, j), id(i + 1, j + 1));
            mesh.add_triangle(id(i, j), id(i + 1, j + 1), id(i, j + 1));
        }
    }
    mesh
}

/// Mesh with the triangles of the cells `(i, j)` selected by `drop` removed.
fn grid_without(n: usize, h: f64, drop: impl Fn(usize, usize) -> bool) -> ParamTriMesh {
    let full = grid_mesh(n, h);
    let mut mesh = ParamTriMesh::new();
    for p in &full.positions {
        mesh.add_vertex(*p, [0.0, 0.0]);
    }
    for (t, tri) in full.triangles.iter().enumerate() {
        let cell = t / 2;
        if !drop(cell % (n - 1), cell / (n - 1)) {
            mesh.add_triangle(tri[0], tri[1], tri[2]);
        }
    }
    mesh
}

/// Octahedron with outward winding: a closed mesh.
fn octahedron() -> ParamTriMesh {
    let mut mesh = ParamTriMesh::new();
    for p in [
        [1.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, -1.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 0.0, -1.0],
    ] {
        mesh.add_vertex(p, [0.0, 0.0]);
    }
    for t in [
        [0, 2, 4],
        [2, 1, 4],
        [1, 3, 4],
        [3, 0, 4],
        [2, 0, 5],
        [1, 2, 5],
        [3, 1, 5],
        [0, 3, 5],
    ] {
        mesh.add_triangle(t[0], t[1], t[2]);
    }
    mesh
}

/// Facets around a half cylinder of radius 1 (outward winding): `segments` facets
/// around the axis and `rows` quads along it, every quad split into two triangles.
/// A strip of planar facets is exactly developable.
fn cylinder_strip(segments: usize, rows: usize, height: f64) -> ParamTriMesh {
    let mut mesh = ParamTriMesh::new();
    for j in 0..=rows {
        for i in 0..=segments {
            let theta = PI * i as f64 / segments as f64;
            mesh.add_vertex([theta.cos(), theta.sin(), height * j as f64], [0.0, 0.0]);
        }
    }
    let id = |i: usize, j: usize| j * (segments + 1) + i;
    for j in 0..rows {
        for i in 0..segments {
            mesh.add_triangle(id(i, j), id(i + 1, j), id(i + 1, j + 1));
            mesh.add_triangle(id(i, j), id(i + 1, j + 1), id(i, j + 1));
        }
    }
    mesh
}

/// Spherical cap of angular radius `half_angle` around +z: a centre vertex and
/// `rings` rings of `6 * ring` vertices, triangulated like a hexagonal disc.
fn spherical_cap(rings: usize, half_angle: f64) -> ParamTriMesh {
    let mut mesh = ParamTriMesh::new();
    mesh.add_vertex([0.0, 0.0, 1.0], [0.0, 0.0]);
    let mut ring_start = vec![0usize];
    for ring in 1..=rings {
        ring_start.push(mesh.vertex_count());
        let polar = half_angle * ring as f64 / rings as f64;
        for k in 0..6 * ring {
            let azimuth = TAU * k as f64 / (6 * ring) as f64;
            mesh.add_vertex(
                [
                    polar.sin() * azimuth.cos(),
                    polar.sin() * azimuth.sin(),
                    polar.cos(),
                ],
                [0.0, 0.0],
            );
        }
    }
    // Centre fan.
    for k in 0..6 {
        mesh.add_triangle(0, 1 + k, 1 + (k + 1) % 6);
    }
    // Strips between consecutive rings.
    for ring in 1..rings {
        let (inner, outer) = (ring_start[ring], ring_start[ring + 1]);
        let (inner_count, outer_count) = (6 * ring, 6 * (ring + 1));
        let (mut a, mut b) = (0usize, 0usize);
        while a < inner_count || b < outer_count {
            let advance_inner = if a >= inner_count {
                false
            } else if b >= outer_count {
                true
            } else {
                // Advance the ring whose next vertex has the smaller azimuth.
                (a + 1) as f64 / inner_count as f64 <= (b + 1) as f64 / outer_count as f64
            };
            let (ia, ib) = (inner + a % inner_count, outer + b % outer_count);
            if advance_inner {
                let next = inner + (a + 1) % inner_count;
                mesh.add_triangle(ia, ib, next);
                a += 1;
            } else {
                let next = outer + (b + 1) % outer_count;
                mesh.add_triangle(ia, ib, next);
                b += 1;
            }
        }
    }
    mesh
}

/// Rotate (about a fixed axis) and translate the positions of `mesh`.
fn rigidly_moved(mesh: &ParamTriMesh) -> ParamTriMesh {
    let (c, s) = (0.8f64, 0.6f64);
    let mut moved = mesh.clone();
    for p in &mut moved.positions {
        // Rotation about the axis (1, 2, 2) / 3 by the angle with cos = c, sin = s.
        let k = [1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0];
        let dot = k[0] * p[0] + k[1] * p[1] + k[2] * p[2];
        let cross = [
            k[1] * p[2] - k[2] * p[1],
            k[2] * p[0] - k[0] * p[2],
            k[0] * p[1] - k[1] * p[0],
        ];
        for a in 0..3 {
            p[a] = p[a] * c + cross[a] * s + k[a] * dot * (1.0 - c);
        }
        p[0] += 3.0;
        p[1] -= 2.0;
        p[2] += 7.5;
    }
    moved
}

fn signed_uv_areas(uvs: &[[f64; 2]], triangles: &[[usize; 3]]) -> Vec<f64> {
    triangles
        .iter()
        .map(|t| {
            let (a, b, c) = (uvs[t[0]], uvs[t[1]], uvs[t[2]]);
            0.5 * ((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]))
        })
        .collect()
}

fn dist2(a: [f64; 2], b: [f64; 2]) -> f64 {
    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2)).sqrt()
}

fn dist3(a: [f64; 3], b: [f64; 3]) -> f64 {
    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
}

fn is_invalid_mesh<T>(result: Result<T>) -> bool {
    matches!(result, Err(Error::InvalidMesh { .. }))
}

// ── Tutte and harmonic maps ─────────────────────────────────────────────────

#[test]
fn grid_boundaries_are_placed_in_loop_order_without_inverted_triangles() {
    for n in [3usize, 5, 9] {
        let mesh = grid_mesh(n, 1.0);
        for (name, uvs) in [
            ("tutte", TutteParameterization::compute(&mesh).unwrap()),
            (
                "harmonic",
                HarmonicMapParameterization::compute(&mesh).unwrap(),
            ),
        ] {
            assert_eq!(uvs.len(), n * n);
            assert_eq!(
                uv_overlap_check(&uvs, &mesh.triangles),
                0,
                "{name} {n}x{n}: inverted triangles"
            );
            let areas = signed_uv_areas(&uvs, &mesh.triangles);
            assert!(
                areas.iter().all(|&a| a > 1e-9),
                "{name} {n}x{n}: every UV triangle needs positive area"
            );
            // Total UV area never exceeds the disc and the boundary is a convex
            // polygon inscribed in it.
            let total: f64 = areas.iter().sum();
            assert!(total < PI && total > 1.0, "{name}: area {total}");
        }
    }
}

#[test]
fn boundary_lies_on_the_circle_and_interior_inside_it() {
    let mesh = grid_mesh(6, 0.5);
    let loop_vertices = boundary_loop(&mesh.triangles).unwrap();
    for uvs in [
        TutteParameterization::compute(&mesh).unwrap(),
        HarmonicMapParameterization::compute(&mesh).unwrap(),
    ] {
        for (v, uv) in uvs.iter().enumerate() {
            let radius = dist2(*uv, [0.0, 0.0]);
            if loop_vertices.contains(&v) {
                assert!((radius - 1.0).abs() < 1e-12, "vertex {v}: {radius}");
            } else {
                assert!(radius < 1.0 - 1e-6, "interior vertex {v}: {radius}");
            }
        }
    }
}

#[test]
fn boundary_angles_follow_arc_length_in_loop_order() {
    // 4 x 2 vertex grid with a long and a short direction.
    let mut mesh = ParamTriMesh::new();
    let xs = [0.0, 1.0, 3.0, 6.0];
    for y in [0.0, 1.0] {
        for x in xs {
            mesh.add_vertex([x, y, 0.0], [0.0, 0.0]);
        }
    }
    for i in 0..3 {
        mesh.add_triangle(i, i + 1, 4 + i + 1);
        mesh.add_triangle(i, 4 + i + 1, 4 + i);
    }
    let uvs = TutteParameterization::compute(&mesh).unwrap();
    // No interior vertices: the whole result is the arc-length placement.
    let loop_vertices = boundary_loop(&mesh.triangles).unwrap();
    assert_eq!(loop_vertices, vec![0, 1, 2, 3, 7, 6, 5, 4]);
    let lengths = [1.0, 2.0, 3.0, 1.0, 3.0, 2.0, 1.0, 1.0];
    let perimeter: f64 = lengths.iter().sum();
    let mut expected_angle = 0.0;
    for (k, &v) in loop_vertices.iter().enumerate() {
        let angle = uvs[v][1].atan2(uvs[v][0]).rem_euclid(TAU);
        assert!(
            (angle - expected_angle).abs() < 1e-12
                || (angle - expected_angle - TAU).abs() < 1e-12
                || (angle + TAU - expected_angle).abs() < 1e-12,
            "vertex {v}: angle {angle}, expected {expected_angle}"
        );
        expected_angle += TAU * lengths[k] / perimeter;
    }
    assert_eq!(uv_overlap_check(&uvs, &mesh.triangles), 0);
}

#[test]
fn map_boundary_loop_to_circle_uses_chord_length() {
    // A 2 x 1 rectangle: edge lengths 2, 1, 2, 1 -> angles 0, 120, 180, 300 degrees.
    let positions = [
        [0.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [2.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],
        [9.0, 9.0, 9.0],
    ];
    let uvs = TutteParameterization::map_boundary_loop_to_circle(&positions, &[0, 1, 2, 3], 5);
    for (v, degrees) in [(0usize, 0.0f64), (1, 120.0), (2, 180.0), (3, 300.0)] {
        let angle = degrees.to_radians();
        assert!((uvs[v][0] - angle.cos()).abs() < 1e-12, "vertex {v}");
        assert!((uvs[v][1] - angle.sin()).abs() < 1e-12, "vertex {v}");
    }
    assert_eq!(uvs[4], [0.0, 0.0]);
    // Coincident points fall back to equal spacing; out-of-range ids are skipped.
    let coincident = [[1.0, 1.0, 1.0]; 4];
    let fallback =
        TutteParameterization::map_boundary_loop_to_circle(&coincident, &[0, 1, 2, 3], 4);
    assert!((fallback[1][0]).abs() < 1e-12 && (fallback[1][1] - 1.0).abs() < 1e-12);
    let skipped = TutteParameterization::map_boundary_loop_to_circle(&coincident, &[0, 9], 4);
    assert_eq!(skipped.len(), 4);
    assert!(
        TutteParameterization::map_boundary_loop_to_circle(&[], &[], 2)
            .iter()
            .all(|uv| *uv == [0.0, 0.0])
    );
}

#[test]
fn equally_spaced_boundary_placement_follows_the_given_order_and_skips_bad_ids() {
    // Order 0, 2, 1 puts vertex 2 at a third and vertex 1 at two thirds of the circle.
    let uvs = TutteParameterization::map_boundary_to_circle(&[0, 2, 1], 3);
    for (v, third) in [(0usize, 0.0f64), (2, 1.0), (1, 2.0)] {
        let angle = TAU * third / 3.0;
        assert!((uvs[v][0] - angle.cos()).abs() < 1e-12, "vertex {v}");
        assert!((uvs[v][1] - angle.sin()).abs() < 1e-12, "vertex {v}");
    }
    // An id beyond `n_verts` is skipped instead of panicking.
    let skipped = TutteParameterization::map_boundary_to_circle(&[0, 7, 1], 3);
    assert_eq!(skipped.len(), 3);
    assert_eq!(skipped[2], [0.0, 0.0]);
    assert!((skipped[1][0] - (TAU * 2.0 / 3.0).cos()).abs() < 1e-12);
}

#[test]
fn centre_of_the_symmetric_3x3_grid_maps_to_the_origin() {
    let mesh = grid_mesh(3, 1.0);
    for uvs in [
        TutteParameterization::compute(&mesh).unwrap(),
        HarmonicMapParameterization::compute(&mesh).unwrap(),
    ] {
        assert!(dist2(uvs[4], [0.0, 0.0]) < 1e-9, "centre at {:?}", uvs[4]);
    }
}

#[test]
fn tutte_interior_vertices_are_the_barycentre_of_their_neighbours() {
    let mesh = spherical_cap(4, 0.9);
    let uvs = TutteParameterization::compute(&mesh).unwrap();
    let boundary = boundary_loop(&mesh.triangles).unwrap();
    for v in 0..mesh.vertex_count() {
        if boundary.contains(&v) {
            continue;
        }
        let neighbours = vertex_neighbors(v, &mesh.triangles);
        let mean = neighbours.iter().fold([0.0, 0.0], |acc, &w| {
            [acc[0] + uvs[w][0], acc[1] + uvs[w][1]]
        });
        let k = neighbours.len() as f64;
        assert!(
            dist2(uvs[v], [mean[0] / k, mean[1] / k]) < 1e-9,
            "vertex {v}"
        );
    }
    assert_eq!(uv_overlap_check(&uvs, &mesh.triangles), 0);
}

/// Cotangent weight of the edge `(i, j)`: the sum of the cotangents of the angles
/// opposite to it, computed independently of the implementation.
fn reference_cotangent_weight(mesh: &ParamTriMesh, i: usize, j: usize) -> f64 {
    let mut weight = 0.0;
    for tri in &mesh.triangles {
        if tri.contains(&i) && tri.contains(&j) {
            let k = tri[0] + tri[1] + tri[2] - i - j;
            let u = [
                mesh.positions[i][0] - mesh.positions[k][0],
                mesh.positions[i][1] - mesh.positions[k][1],
                mesh.positions[i][2] - mesh.positions[k][2],
            ];
            let v = [
                mesh.positions[j][0] - mesh.positions[k][0],
                mesh.positions[j][1] - mesh.positions[k][1],
                mesh.positions[j][2] - mesh.positions[k][2],
            ];
            let cross = [
                u[1] * v[2] - u[2] * v[1],
                u[2] * v[0] - u[0] * v[2],
                u[0] * v[1] - u[1] * v[0],
            ];
            let sin = (cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]).sqrt();
            weight += (u[0] * v[0] + u[1] * v[1] + u[2] * v[2]) / sin;
        }
    }
    weight.max(1e-8)
}

#[test]
fn harmonic_interior_vertices_satisfy_the_cotangent_laplace_equation() {
    let mesh = spherical_cap(4, 0.9);
    let uvs = HarmonicMapParameterization::compute(&mesh).unwrap();
    let boundary = boundary_loop(&mesh.triangles).unwrap();
    for v in 0..mesh.vertex_count() {
        if boundary.contains(&v) {
            continue;
        }
        let mut weighted = [0.0, 0.0];
        let mut total = 0.0;
        for w in vertex_neighbors(v, &mesh.triangles) {
            let weight = reference_cotangent_weight(&mesh, v, w);
            weighted[0] += weight * uvs[w][0];
            weighted[1] += weight * uvs[w][1];
            total += weight;
        }
        let target = [weighted[0] / total, weighted[1] / total];
        assert!(
            dist2(uvs[v], target) < 1e-8,
            "vertex {v}: {:?} vs {target:?}",
            uvs[v]
        );
    }
    assert_eq!(uv_overlap_check(&uvs, &mesh.triangles), 0);
}

#[test]
fn harmonic_weights_differ_from_uniform_weights_on_irregular_meshes() {
    // A fan around an off-centre interior vertex: cotangent weights are not 1.
    let mut mesh = ParamTriMesh::new();
    mesh.add_vertex([0.35, 0.1, 0.0], [0.0, 0.0]);
    for k in 0..6 {
        let angle = TAU * k as f64 / 6.0;
        let radius = if k % 2 == 0 { 1.0 } else { 2.0 };
        mesh.add_vertex(
            [radius * angle.cos(), radius * angle.sin(), 0.0],
            [0.0, 0.0],
        );
    }
    for k in 0..6 {
        mesh.add_triangle(0, 1 + k, 1 + (k + 1) % 6);
    }
    let tutte = TutteParameterization::compute(&mesh).unwrap();
    let harmonic = HarmonicMapParameterization::compute(&mesh).unwrap();
    assert!(
        dist2(tutte[0], harmonic[0]) > 1e-3,
        "{:?} vs {:?}",
        tutte[0],
        harmonic[0]
    );
    assert_eq!(uv_overlap_check(&tutte, &mesh.triangles), 0);
    assert_eq!(uv_overlap_check(&harmonic, &mesh.triangles), 0);
}

#[test]
fn concave_patches_embed_without_folds() {
    // L-shaped region: the upper-right 2 x 2 cells of a 5 x 5 vertex grid are gone,
    // so the boundary has a reflex corner.
    let mesh = grid_without(5, 1.0, |i, j| i >= 2 && j >= 2);
    for uvs in [
        TutteParameterization::compute(&mesh).unwrap(),
        HarmonicMapParameterization::compute(&mesh).unwrap(),
    ] {
        assert_eq!(uv_overlap_check(&uvs, &mesh.triangles), 0);
        // Vertices of the removed cells that no triangle uses stay at the origin.
        assert_eq!(uvs[4 * 5 + 4], [0.0, 0.0]);
        assert!(
            signed_uv_areas(&uvs, &mesh.triangles)
                .iter()
                .all(|&a| a > 0.0)
        );
    }
}

#[test]
fn vertex_numbering_does_not_matter() {
    let mesh = grid_mesh(5, 1.0);
    // A fixed relabelling of the vertices: new id = (7 v + 3) mod 25.
    let n = mesh.vertex_count();
    let permutation: Vec<usize> = (0..n).map(|v| (v * 7 + 3) % n).collect();
    let mut sorted = permutation.clone();
    sorted.sort_unstable();
    assert_eq!(sorted, (0..n).collect::<Vec<_>>());
    let mut shuffled = ParamTriMesh::new();
    shuffled.positions = vec![[0.0; 3]; n];
    shuffled.uvs = vec![[0.0; 2]; n];
    shuffled.normals = vec![[0.0; 3]; n];
    for (v, &target) in permutation.iter().enumerate() {
        shuffled.positions[target] = mesh.positions[v];
    }
    for t in &mesh.triangles {
        shuffled.add_triangle(permutation[t[0]], permutation[t[1]], permutation[t[2]]);
    }
    type Solver = fn(&ParamTriMesh) -> Result<Vec<[f64; 2]>>;
    let solvers: [Solver; 2] = [
        TutteParameterization::compute,
        HarmonicMapParameterization::compute,
    ];
    for solve in solvers {
        let original = solve(&mesh).unwrap();
        let relabelled = solve(&shuffled).unwrap();
        assert_eq!(uv_overlap_check(&relabelled, &shuffled.triangles), 0);
        // The start of the boundary loop depends on the ids, which rotates the
        // circle; the distance of every vertex from the centre does not.
        for v in 0..n {
            let r0 = dist2(original[v], [0.0; 2]);
            let r1 = dist2(relabelled[permutation[v]], [0.0; 2]);
            assert!((r0 - r1).abs() < 1e-9, "vertex {v}: {r0} vs {r1}");
        }
    }
}

#[test]
fn patches_with_holes_keep_the_hole_inside_the_disc() {
    // 5 x 5 grid without its centre cells: an annulus.
    let mesh = grid_without(5, 1.0, |i, j| (1..=2).contains(&i) && (1..=2).contains(&j));
    let loops = boundary_loops(&mesh.triangles);
    assert_eq!(loops.len(), 2);
    let uvs = TutteParameterization::compute(&mesh).unwrap();
    for &v in &loops[0] {
        assert!((dist2(uvs[v], [0.0; 2]) - 1.0).abs() < 1e-12);
    }
    for &v in &loops[1] {
        assert!(dist2(uvs[v], [0.0; 2]) < 1.0 - 1e-6, "hole vertex {v}");
        assert!(uvs[v][0].is_finite() && uvs[v][1].is_finite());
    }
}

#[test]
fn invalid_meshes_are_reported_instead_of_returning_zeros() {
    let closed = octahedron();
    assert!(is_invalid_mesh(TutteParameterization::compute(&closed)));
    assert!(is_invalid_mesh(HarmonicMapParameterization::compute(
        &closed
    )));
    assert!(is_invalid_mesh(LscmParameterization::compute(&closed)));

    // Two disconnected triangles.
    let mut pieces = ParamTriMesh::new();
    for k in 0..2 {
        let o = 5.0 * k as f64;
        pieces.add_vertex([o, 0.0, 0.0], [0.0; 2]);
        pieces.add_vertex([o + 1.0, 0.0, 0.0], [0.0; 2]);
        pieces.add_vertex([o, 1.0, 0.0], [0.0; 2]);
        pieces.add_triangle(3 * k, 3 * k + 1, 3 * k + 2);
    }
    assert!(is_invalid_mesh(TutteParameterization::compute(&pieces)));
    assert!(is_invalid_mesh(HarmonicMapParameterization::compute(
        &pieces
    )));
    assert!(is_invalid_mesh(LscmParameterization::compute(&pieces)));

    // Vertices without triangles.
    let mut points = ParamTriMesh::new();
    points.add_vertex([0.0; 3], [0.0; 2]);
    assert!(is_invalid_mesh(TutteParameterization::compute(&points)));
    assert!(is_invalid_mesh(LscmParameterization::compute(&points)));

    // A triangle that refers to a missing vertex.
    let mut dangling = grid_mesh(2, 1.0);
    dangling.add_triangle(0, 1, 40);
    for result in [
        TutteParameterization::compute(&dangling),
        HarmonicMapParameterization::compute(&dangling),
        LscmParameterization::compute(&dangling),
    ] {
        assert!(matches!(
            result,
            Err(Error::IndexOutOfBounds { index: 40, len: 4 })
        ));
    }

    // The empty mesh has nothing to parameterize and is not an error.
    let empty = ParamTriMesh::new();
    assert!(TutteParameterization::compute(&empty).unwrap().is_empty());
    assert!(
        HarmonicMapParameterization::compute(&empty)
            .unwrap()
            .is_empty()
    );
    assert!(LscmParameterization::compute(&empty).unwrap().is_empty());
}

#[test]
fn unused_vertices_keep_the_origin() {
    let mut mesh = grid_mesh(4, 1.0);
    let spare = mesh.add_vertex([9.0, 9.0, 9.0], [0.5, 0.5]);
    for uvs in [
        TutteParameterization::compute(&mesh).unwrap(),
        HarmonicMapParameterization::compute(&mesh).unwrap(),
        LscmParameterization::compute(&mesh).unwrap(),
    ] {
        assert_eq!(uvs[spare], [0.0, 0.0]);
        assert_eq!(uv_overlap_check(&uvs, &mesh.triangles), 0);
    }
}

#[test]
fn a_large_grid_is_solved_quickly_and_accurately() {
    let mesh = grid_mesh(41, 0.25);
    for uvs in [
        TutteParameterization::compute(&mesh).unwrap(),
        HarmonicMapParameterization::compute(&mesh).unwrap(),
    ] {
        assert_eq!(uv_overlap_check(&uvs, &mesh.triangles), 0);
        let total: f64 = signed_uv_areas(&uvs, &mesh.triangles).iter().sum();
        assert!(total > 2.5 && total < PI, "area {total}");
    }
}

// ── LSCM ─────────────────────────────────────────────────────────────────────

#[test]
fn planar_grid_is_unfolded_isometrically() {
    for (n, h) in [(5usize, 1.0f64), (5, 0.37), (9, 2.0)] {
        let mesh = grid_mesh(n, h);
        let uvs = LscmParameterization::compute(&mesh).unwrap();
        let energy = LscmParameterization::conformal_energy(&mesh.positions, &uvs, &mesh.triangles);
        assert!(energy < 1e-8, "{n}x{n} h={h}: conformal energy {energy}");
        let areas = signed_uv_areas(&uvs, &mesh.triangles);
        assert!(
            areas.iter().all(|&a| a > 0.0),
            "every UV triangle needs positive area"
        );
        assert_eq!(uv_overlap_check(&uvs, &mesh.triangles), 0);
        // Pins at their true distance: the layout reproduces the mesh lengths.
        for t in &mesh.triangles {
            for k in 0..3 {
                let (a, b) = (t[k], t[(k + 1) % 3]);
                let expected = dist3(mesh.positions[a], mesh.positions[b]);
                assert!(
                    (dist2(uvs[a], uvs[b]) - expected).abs() < 1e-6 * (1.0 + expected),
                    "edge {a}-{b}"
                );
            }
        }
        // The farthest pair of boundary vertices are opposite corners.
        assert_eq!(uvs[0], [0.0, 0.0]);
        let diagonal = dist3(mesh.positions[0], mesh.positions[n * n - 1]);
        assert!((uvs[n * n - 1][0] - diagonal).abs() < 1e-9 && uvs[n * n - 1][1].abs() < 1e-9);
    }
}

#[test]
fn lscm_does_not_depend_on_the_pose_of_the_mesh() {
    let mesh = grid_mesh(6, 0.8);
    let moved = rigidly_moved(&mesh);
    let reference = LscmParameterization::compute(&mesh).unwrap();
    let uvs = LscmParameterization::compute(&moved).unwrap();
    for (a, b) in reference.iter().zip(uvs.iter()) {
        assert!(dist2(*a, *b) < 1e-8, "{a:?} vs {b:?}");
    }
    let energy = LscmParameterization::conformal_energy(&moved.positions, &uvs, &moved.triangles);
    assert!(energy < 1e-8);
}

#[test]
fn developable_cylinder_strip_unfolds_with_zero_energy_and_uniform_scale() {
    let mesh = cylinder_strip(6, 3, 0.7);
    let uvs = LscmParameterization::compute(&mesh).unwrap();
    let energy = LscmParameterization::conformal_energy(&mesh.positions, &uvs, &mesh.triangles);
    assert!(energy < 1e-8, "conformal energy {energy}");
    assert!(
        signed_uv_areas(&uvs, &mesh.triangles)
            .iter()
            .all(|&a| a > 0.0)
    );
    // Conformal and (as the strip is developable) rigid up to one global scale:
    // every edge shrinks by the same factor.
    let mut ratios = Vec::new();
    for t in &mesh.triangles {
        for k in 0..3 {
            let (a, b) = (t[k], t[(k + 1) % 3]);
            ratios.push(dist2(uvs[a], uvs[b]) / dist3(mesh.positions[a], mesh.positions[b]));
        }
    }
    let first = ratios[0];
    assert!(first > 0.0 && first <= 1.0 + 1e-9, "scale {first}");
    for r in &ratios {
        assert!((r - first).abs() < 1e-6, "{r} vs {first}");
    }
    let distortion = mean_angle_distortion(&mesh.positions, &uvs, &mesh.triangles);
    assert!(distortion < 1e-6, "angle distortion {distortion}");
}

#[test]
fn curved_patch_gets_a_valid_low_distortion_layout() {
    let mesh = spherical_cap(5, 0.9);
    let uvs = LscmParameterization::compute(&mesh).unwrap();
    assert!(uvs.iter().all(|uv| uv[0].is_finite() && uv[1].is_finite()));
    assert_eq!(uv_overlap_check(&uvs, &mesh.triangles), 0);
    // A sphere is not developable, but a cap this small is nearly flat.
    let distortion = mean_angle_distortion(&mesh.positions, &uvs, &mesh.triangles);
    assert!(distortion < 0.05, "mean angle distortion {distortion}");
    let energy = LscmParameterization::conformal_energy(&mesh.positions, &uvs, &mesh.triangles);
    assert!(energy > 0.0 && energy < 0.05, "energy {energy}");
}

#[test]
fn explicit_pins_are_honoured() {
    let mesh = grid_mesh(5, 1.0);
    let pins = [(12usize, [2.0, 3.0]), (3, [2.0, 7.0])];
    let uvs = lscm_solve(&mesh.positions, &mesh.triangles, pins).unwrap();
    assert_eq!(uvs[12], [2.0, 3.0]);
    assert_eq!(uvs[3], [2.0, 7.0]);
    // The pins fix a rotation by 90 degrees about the first pin and a scaling:
    // the distance of 4 between the pins is mapped from the 3-D distance
    // |p12 - p3| = sqrt(1 + 4) -> the map is the similarity that does this.
    let energy = LscmParameterization::conformal_energy(&mesh.positions, &uvs, &mesh.triangles);
    assert!(energy < 1e-8, "energy {energy}");
    let scale = 4.0 / dist3(mesh.positions[12], mesh.positions[3]);
    for t in &mesh.triangles {
        let (a, b) = (t[0], t[1]);
        let ratio = dist2(uvs[a], uvs[b]) / dist3(mesh.positions[a], mesh.positions[b]);
        assert!((ratio - scale).abs() < 1e-6, "{ratio} vs {scale}");
    }
    assert!(
        signed_uv_areas(&uvs, &mesh.triangles)
            .iter()
            .all(|&a| a > 0.0)
    );
}

#[test]
fn invalid_lscm_inputs_are_rejected() {
    let mesh = grid_mesh(3, 1.0);
    let ok = [(0usize, [0.0, 0.0]), (8, [1.0, 0.0])];
    assert!(lscm_solve(&mesh.positions, &mesh.triangles, ok).is_ok());
    // Equal pins.
    assert!(is_invalid_mesh(lscm_solve(
        &mesh.positions,
        &mesh.triangles,
        [(1, [0.0, 0.0]), (1, [1.0, 0.0])]
    )));
    // Pin outside the mesh / not used by any triangle.
    assert!(matches!(
        lscm_solve(
            &mesh.positions,
            &mesh.triangles,
            [(0, [0.0; 2]), (99, [1.0, 0.0])]
        ),
        Err(Error::IndexOutOfBounds { index: 99, len: 9 })
    ));
    let mut spare = mesh.clone();
    let unused = spare.add_vertex([5.0, 5.0, 0.0], [0.0; 2]);
    assert!(is_invalid_mesh(lscm_solve(
        &spare.positions,
        &spare.triangles,
        [(0, [0.0; 2]), (unused, [1.0, 0.0])]
    )));
    // Non-finite coordinates.
    let mut broken = mesh.clone();
    broken.positions[4][2] = f64::NAN;
    assert!(matches!(
        lscm_solve(&broken.positions, &broken.triangles, ok),
        Err(Error::DegenerateGeometry { .. })
    ));
    // Only degenerate triangles.
    let flat = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
    assert!(matches!(
        lscm_solve(&flat, &[[0, 1, 2]], [(0, [0.0; 2]), (2, [1.0, 0.0])]),
        Err(Error::DegenerateGeometry { .. })
    ));
}

#[test]
fn degenerate_triangles_are_ignored_when_the_rest_determines_the_map() {
    let mut mesh = grid_mesh(5, 1.0);
    // A zero-area sliver on three collinear interior vertices (11, 12, 13).
    mesh.add_triangle(11, 12, 13);
    let uvs = LscmParameterization::compute(&mesh).unwrap();
    let regular = &mesh.triangles[..mesh.triangles.len() - 1];
    let energy = LscmParameterization::conformal_energy(&mesh.positions, &uvs, regular);
    assert!(energy < 1e-8, "energy {energy}");
    assert!(signed_uv_areas(&uvs, regular).iter().all(|&a| a > 0.0));

    // A vertex that only degenerate triangles use cannot be placed.
    let mut dangling = grid_mesh(3, 1.0);
    let a = dangling.add_vertex([-1.0, 0.0, 0.0], [0.0; 2]);
    let b = dangling.add_vertex([-2.0, 0.0, 0.0], [0.0; 2]);
    dangling.add_triangle(0, a, b);
    assert!(matches!(
        lscm_solve(
            &dangling.positions,
            &dangling.triangles,
            [(0, [0.0, 0.0]), (8, [1.0, 0.0])]
        ),
        Err(Error::DegenerateGeometry { .. })
    ));
}

#[test]
fn patches_joined_only_by_a_degenerate_triangle_are_rejected() {
    // Two 3 x 3 grids four units apart; the collinear triangle (2, 9, 10) is the only
    // link, so the second grid would be free to move, turn and scale.
    let left = grid_mesh(3, 1.0);
    let mut mesh = ParamTriMesh::new();
    for p in &left.positions {
        mesh.add_vertex(*p, [0.0; 2]);
    }
    for p in &left.positions {
        mesh.add_vertex([p[0] + 4.0, p[1], p[2]], [0.0; 2]);
    }
    for t in &left.triangles {
        mesh.add_triangle(t[0], t[1], t[2]);
        mesh.add_triangle(t[0] + 9, t[1] + 9, t[2] + 9);
    }
    mesh.add_triangle(2, 9, 10);
    assert!(matches!(
        lscm_solve(
            &mesh.positions,
            &mesh.triangles,
            [(0, [0.0, 0.0]), (8, [1.0, 0.0])]
        ),
        Err(Error::DegenerateGeometry { .. })
    ));
    assert!(matches!(
        LscmParameterization::compute(&mesh),
        Err(Error::DegenerateGeometry { .. })
    ));
    // A pin that only a degenerate triangle uses cannot fix anything either.
    let mut pinned_by_sliver = grid_mesh(3, 1.0);
    let stray = pinned_by_sliver.add_vertex([3.0, 0.0, 0.0], [0.0; 2]);
    pinned_by_sliver.add_triangle(1, 2, stray);
    assert!(matches!(
        lscm_solve(
            &pinned_by_sliver.positions,
            &pinned_by_sliver.triangles,
            [(0, [0.0, 0.0]), (stray, [3.0, 0.0])]
        ),
        Err(Error::DegenerateGeometry { .. })
    ));
}

#[test]
fn closed_meshes_can_still_be_solved_with_explicit_pins() {
    let mesh = octahedron();
    let uvs = lscm_solve(
        &mesh.positions,
        &mesh.triangles,
        [(0, [0.0, 0.0]), (1, [2.0, 0.0])],
    )
    .unwrap();
    assert!(uvs.iter().all(|uv| uv[0].is_finite() && uv[1].is_finite()));
    assert_eq!(uvs[0], [0.0, 0.0]);
    assert_eq!(uvs[1], [2.0, 0.0]);
}

#[test]
fn farthest_boundary_pair_is_pinned_at_the_true_distance() {
    let mesh = grid_mesh(5, 1.0);
    let boundary: Vec<usize> = boundary_loops(&mesh.triangles)
        .into_iter()
        .flatten()
        .collect();
    let pins = farthest_pair_pins(&mesh.positions, boundary).unwrap();
    // Ties between the two diagonals resolve to the smaller ids.
    assert_eq!(pins[0].0, 0);
    assert_eq!(pins[1].0, 24);
    assert_eq!(pins[0].1, [0.0, 0.0]);
    assert!((pins[1].1[0] - 32.0_f64.sqrt()).abs() < 1e-12);
    assert_eq!(pins[1].1[1], 0.0);
    // Fewer than two candidates and invalid ids are errors.
    assert!(matches!(
        farthest_pair_pins(&mesh.positions, vec![3, 3]),
        Err(Error::TooFewPoints {
            required: 2,
            actual: 1
        })
    ));
    assert!(matches!(
        farthest_pair_pins(&mesh.positions, vec![0, 25]),
        Err(Error::IndexOutOfBounds { index: 25, len: 25 })
    ));
    // Coincident candidates still give two distinct pins with unit separation.
    let coincident = [[1.0; 3]; 3];
    let pins = farthest_pair_pins(&coincident, vec![0, 1, 2]).unwrap();
    assert_eq!((pins[0].0, pins[1].0), (0, 1));
    assert_eq!(pins[1].1, [1.0, 0.0]);
}

#[test]
fn small_candidate_sets_get_the_exact_diameter_pair() {
    // A double sweep from vertex 0 would pick the pair (0, 1) at squared distance
    // 133; the true diameter is the pair (2, 3) at squared distance 164.
    let positions = [
        [9.0, 9.0, 10.0],
        [3.0, 5.0, 1.0],
        [8.0, 1.0, 9.0],
        [0.0, 9.0, 3.0],
    ];
    let pins = farthest_pair_pins(&positions, vec![3, 1, 0, 2, 2]).unwrap();
    assert_eq!((pins[0].0, pins[1].0), (2, 3));
    assert_eq!(pins[0].1, [0.0, 0.0]);
    assert!((pins[1].1[0] - 164.0_f64.sqrt()).abs() < 1e-12);
    // The result does not depend on the order in which the candidates are listed.
    let reordered = farthest_pair_pins(&positions, vec![0, 2, 3, 1]).unwrap();
    assert_eq!((reordered[0].0, reordered[1].0), (2, 3));
}

/// A centre vertex (id 0) with `ring` vertices equally spaced on the unit circle
/// around it (ids `1..=ring`), moved by a rigid motion: every ring vertex is
/// equally far from the centre and every antipodal ring pair is equally far apart,
/// up to rounding noise.
fn ring_positions(ring: usize, moved: bool) -> Vec<[f64; 3]> {
    let mut mesh = ParamTriMesh::new();
    mesh.add_vertex([0.0, 0.0, 0.0], [0.0; 2]);
    for k in 0..ring {
        let angle = TAU * k as f64 / ring as f64;
        mesh.add_vertex([angle.cos(), angle.sin(), 0.0], [0.0; 2]);
    }
    if moved {
        mesh = rigidly_moved(&mesh);
    }
    mesh.positions
}

#[test]
fn exactly_tied_pairs_resolve_to_the_smallest_ids_whatever_the_pose() {
    // 40 equally spaced vertices: the 20 antipodal pairs are tied for the diameter.
    for moved in [false, true] {
        let positions = ring_positions(40, moved);
        let candidates: Vec<usize> = (1..=40).collect();
        let pins = farthest_pair_pins(&positions, candidates).unwrap();
        assert_eq!(
            (pins[0].0, pins[1].0),
            (1, 21),
            "moved = {moved}: rounding noise must not pick another diameter"
        );
        assert!((pins[1].1[0] - 2.0).abs() < 1e-12);
    }
    // The same for the corners of a square, whose two diagonals tie.
    let square = [
        [0.0, 0.0, 0.0],
        [0.3, 0.0, 0.0],
        [0.3, 0.3, 0.0],
        [0.0, 0.3, 0.0],
    ];
    let mut mesh = ParamTriMesh::new();
    for p in square {
        mesh.add_vertex(p, [0.0; 2]);
    }
    for positions in [mesh.positions.clone(), rigidly_moved(&mesh).positions] {
        let pins = farthest_pair_pins(&positions, vec![0, 1, 2, 3]).unwrap();
        assert_eq!((pins[0].0, pins[1].0), (0, 2));
    }
}

#[test]
fn tied_sweep_candidates_resolve_to_the_smallest_ids_whatever_the_pose() {
    // More candidates than the exact search takes: every ring vertex is equally far
    // from the centre, so the first sweep is one big tie.
    let ring = 4200;
    let candidates: Vec<usize> = (0..=ring).collect();
    let reference = farthest_pair_pins(&ring_positions(ring, false), candidates.clone()).unwrap();
    // Farthest from the centre: ring vertex 1 (smallest id); farthest from that:
    // its antipode, ring vertex 1 + ring / 2.
    assert_eq!((reference[0].0, reference[1].0), (1, 1 + ring / 2));
    let moved = farthest_pair_pins(&ring_positions(ring, true), candidates).unwrap();
    assert_eq!((moved[0].0, moved[1].0), (1, 1 + ring / 2));
    assert!((moved[1].1[0] - 2.0).abs() < 1e-9);
}

#[test]
fn large_candidate_sets_use_the_double_sweep_and_find_the_ends_of_a_strip() {
    // More candidates than the exact search takes: points on a long strip, listed
    // in a scrambled order.  The double sweep finds the two ends.
    let count = 5000usize;
    let positions: Vec<[f64; 3]> = (0..count)
        .map(|k| [0.5 * k as f64, 0.01 * ((k * 7919) % 13) as f64, 0.0])
        .collect();
    let mut candidates: Vec<usize> = (0..count).collect();
    candidates.reverse();
    let pins = farthest_pair_pins(&positions, candidates).unwrap();
    assert_eq!((pins[0].0, pins[1].0), (0, count - 1));
    assert_eq!(pins[0].1, [0.0, 0.0]);
    assert!((pins[1].1[0] - dist3(positions[0], positions[count - 1])).abs() < 1e-9);
}

/// Dense reference for the LSCM normal equations, assembled from the two real rows
/// of every triangle (`Re` and `Im` of the complex conformality condition).
fn dense_reference(
    mesh: &ParamTriMesh,
    system: &NormalEquations,
    pins: [(usize, [f64; 2]); 2],
) -> (Vec<Vec<f64>>, Vec<f64>) {
    let free = &system.free_vertices;
    let size = 2 * free.len();
    let dof = |v: usize, axis: usize| free.iter().position(|&f| f == v).map(|s| 2 * s + axis);
    let pinned_value = |v: usize, axis: usize| {
        pins.iter()
            .find(|&&(p, _)| p == v)
            .map_or(0.0, |&(_, uv)| uv[axis])
    };
    let mut matrix = vec![vec![0.0; size]; size];
    let mut rhs = vec![0.0; size];
    for tri in &mesh.triangles {
        let Some(c) = conformal_row(
            mesh.positions[tri[0]],
            mesh.positions[tri[1]],
            mesh.positions[tri[2]],
        ) else {
            continue;
        };
        // Real rows over [u0, v0, u1, v1, u2, v2].
        let re = [c[0].0, -c[0].1, c[1].0, -c[1].1, c[2].0, -c[2].1];
        let im = [c[0].1, c[0].0, c[1].1, c[1].0, c[2].1, c[2].0];
        for row in [re, im] {
            for a in 0..6 {
                let (va, axis_a) = (tri[a / 2], a % 2);
                let Some(da) = dof(va, axis_a) else { continue };
                for b in 0..6 {
                    let (vb, axis_b) = (tri[b / 2], b % 2);
                    match dof(vb, axis_b) {
                        Some(db) => matrix[da][db] += row[a] * row[b],
                        None => rhs[da] -= row[a] * row[b] * pinned_value(vb, axis_b),
                    }
                }
            }
        }
    }
    (matrix, rhs)
}

#[test]
fn normal_equations_match_a_dense_assembly_from_the_real_rows() {
    // An irregular mesh: the spherical cap has no symmetric structure to hide bugs.
    let mesh = spherical_cap(3, 1.2);
    let pins = [(5usize, [0.3, -0.2]), (11, [1.4, 0.9])];
    let system = assemble_normal_equations(&mesh.positions, &mesh.triangles, pins).unwrap();
    let (matrix, rhs) = dense_reference(&mesh, &system, pins);
    let sparse = system.matrix.to_dense();
    assert_eq!(sparse.len(), matrix.len());
    for (r, row) in matrix.iter().enumerate() {
        for (c, &expected) in row.iter().enumerate() {
            assert!(
                (sparse[r][c] - expected).abs() < 1e-12,
                "entry ({r}, {c}): {} vs {expected}",
                sparse[r][c]
            );
        }
        assert!((system.rhs[r] - rhs[r]).abs() < 1e-12, "rhs {r}");
    }
    // Symmetric and positive semi-definite.
    assert!(system.matrix.is_symmetric(1e-12));
    for seed in 1..6 {
        let x: Vec<f64> = (0..matrix.len())
            .map(|k| ((k * 31 + seed * 17) % 23) as f64 - 11.0)
            .collect();
        let mx = system.matrix.matvec(&x);
        let energy: f64 = x.iter().zip(mx.iter()).map(|(a, b)| a * b).sum();
        assert!(energy >= -1e-9, "x^T M x = {energy}");
    }
}

#[test]
fn conformal_energy_vanishes_exactly_for_similarities() {
    let mesh = grid_mesh(4, 1.0);
    let rotation = |p: [f64; 3], scale: f64| -> [f64; 2] {
        let (c, s) = (0.6f64, 0.8f64);
        [
            scale * (c * p[0] - s * p[1]) + 5.0,
            scale * (s * p[0] + c * p[1]) - 1.0,
        ]
    };
    let similar: Vec<[f64; 2]> = mesh.positions.iter().map(|&p| rotation(p, 2.5)).collect();
    let energy = LscmParameterization::conformal_energy(&mesh.positions, &similar, &mesh.triangles);
    assert!(energy < 1e-12, "{energy}");
    let sheared: Vec<[f64; 2]> = mesh
        .positions
        .iter()
        .map(|p| [p[0] + 0.5 * p[1], p[1]])
        .collect();
    let energy = LscmParameterization::conformal_energy(&mesh.positions, &sheared, &mesh.triangles);
    assert!(energy > 0.1, "{energy}");
    // A reflection is conformal as a map of the plane, so it is not penalised.
    let mirrored: Vec<[f64; 2]> = mesh.positions.iter().map(|p| [-p[0], p[1]]).collect();
    let energy =
        LscmParameterization::conformal_energy(&mesh.positions, &mirrored, &mesh.triangles);
    assert!(energy < 1e-12, "{energy}");
}

#[test]
fn single_triangles_are_handled_by_every_parameterization() {
    // A scalene triangle tilted in space: nothing but the boundary exists.
    let mut mesh = ParamTriMesh::new();
    mesh.add_vertex([0.0, 0.0, 0.0], [0.0; 2]);
    mesh.add_vertex([3.0, 0.0, 1.0], [0.0; 2]);
    mesh.add_vertex([0.5, 2.0, 2.0], [0.0; 2]);
    mesh.add_triangle(0, 1, 2);

    for uvs in [
        TutteParameterization::compute(&mesh).unwrap(),
        HarmonicMapParameterization::compute(&mesh).unwrap(),
    ] {
        for uv in &uvs {
            assert!((dist2(*uv, [0.0; 2]) - 1.0).abs() < 1e-12);
        }
        assert_eq!(uv_overlap_check(&uvs, &mesh.triangles), 0);
    }

    // LSCM reproduces the triangle exactly (it is trivially developable).
    let uvs = LscmParameterization::compute(&mesh).unwrap();
    for (a, b) in [(0usize, 1usize), (1, 2), (2, 0)] {
        let expected = dist3(mesh.positions[a], mesh.positions[b]);
        assert!(
            (dist2(uvs[a], uvs[b]) - expected).abs() < 1e-9,
            "edge {a}-{b}"
        );
    }
    assert!(signed_uv_areas(&uvs, &mesh.triangles)[0] > 0.0);
    let energy = LscmParameterization::conformal_energy(&mesh.positions, &uvs, &mesh.triangles);
    assert!(energy < 1e-12, "energy {energy}");
}

#[test]
fn lscm_scales_to_larger_meshes() {
    let mesh = grid_mesh(33, 0.1);
    let uvs = LscmParameterization::compute(&mesh).unwrap();
    let energy = LscmParameterization::conformal_energy(&mesh.positions, &uvs, &mesh.triangles);
    assert!(energy < 1e-8, "energy {energy}");
    assert_eq!(uv_overlap_check(&uvs, &mesh.triangles), 0);
}
