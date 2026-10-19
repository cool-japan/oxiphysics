//! Regression tests for the butterfly subdivision stencils.
//!
//! An earlier version dropped the wing term of the stencil and subtracted a
//! sixteenth of every ring neighbour instead, so the weights summed to `7/8` (and
//! to `5/4` on a tetrahedron): the output shrank towards the origin and was not
//! invariant under translation.

use super::*;
use crate::test_meshes::{flat_grid, icosphere};
use std::collections::HashSet;

const TOLERANCE: f64 = 1e-12;

const IDENTITY: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

fn mesh_of((vertices, faces): (Vec<[f64; 3]>, Vec<[usize; 3]>)) -> SubdivMesh {
    SubdivMesh::new(vertices, faces)
}

/// Regular tetrahedron inscribed in the cube `[-1, 1]^3`, outward oriented.
fn tetrahedron() -> SubdivMesh {
    SubdivMesh::new(
        vec![
            [1.0, 1.0, 1.0],
            [-1.0, -1.0, 1.0],
            [-1.0, 1.0, -1.0],
            [1.0, -1.0, -1.0],
        ],
        vec![[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]],
    )
}

fn assert_close(actual: [f64; 3], expected: [f64; 3], what: &str) {
    for k in 0..3 {
        assert!(
            (actual[k] - expected[k]).abs() < TOLERANCE,
            "{what}: component {k} is {} but {} was expected ({actual:?} vs {expected:?})",
            actual[k],
            expected[k]
        );
    }
}

/// Whether `vertices` contains a point within `tolerance` of `target`.
fn contains_point(vertices: &[[f64; 3]], target: [f64; 3], tolerance: f64) -> bool {
    vertices.iter().any(|p| len3(sub3(*p, target)) < tolerance)
}

/// Applies the affine map `x -> M x + t` to every vertex.
fn affine(mesh: &SubdivMesh, m: [[f64; 3]; 3], t: [f64; 3]) -> SubdivMesh {
    let vertices = mesh
        .vertices
        .iter()
        .map(|p| {
            let mut q = t;
            for (row, value) in q.iter_mut().enumerate() {
                for (col, coordinate) in p.iter().enumerate() {
                    *value += m[row][col] * coordinate;
                }
            }
            q
        })
        .collect();
    SubdivMesh::new(vertices, mesh.faces.clone())
}

/// A closed, consistently oriented mesh uses every directed edge once and each
/// directed edge has its opposite.
fn assert_closed_and_oriented(mesh: &SubdivMesh) {
    let mut directed: HashSet<(usize, usize)> = HashSet::new();
    for face in &mesh.faces {
        for k in 0..3 {
            let edge = (face[k], face[(k + 1) % 3]);
            assert!(
                directed.insert(edge),
                "directed edge {edge:?} appears twice"
            );
        }
    }
    for &(a, b) in &directed {
        assert!(
            directed.contains(&(b, a)),
            "edge ({a}, {b}) has no opposite half-edge"
        );
    }
}

// ── Affine invariance ────────────────────────────────────────────────────────

fn assert_commutes_with_translation(mesh: &SubdivMesh, shift: [f64; 3]) {
    let moved = affine(mesh, IDENTITY, shift);
    let refined = butterfly_subdivision(mesh);
    let refined_moved = butterfly_subdivision(&moved);
    assert_eq!(refined.faces, refined_moved.faces);
    for (index, (p, q)) in refined
        .vertices
        .iter()
        .zip(&refined_moved.vertices)
        .enumerate()
    {
        assert_close(
            *q,
            [p[0] + shift[0], p[1] + shift[1], p[2] + shift[2]],
            &format!("vertex {index}"),
        );
    }
}

#[test]
fn translating_the_input_translates_the_output_on_an_icosphere() {
    assert_commutes_with_translation(&mesh_of(icosphere(1)), [3.5, -2.25, 10.0]);
}

#[test]
fn translating_the_input_translates_the_output_on_a_tetrahedron() {
    assert_commutes_with_translation(&tetrahedron(), [100.0, -50.0, 25.0]);
}

#[test]
fn general_affine_maps_commute_with_the_subdivision() {
    // Non-orthogonal linear part: a shear combined with anisotropic scaling.
    let m = [[2.0, 0.5, -0.3], [0.1, 0.7, 0.4], [-0.6, 0.2, 1.5]];
    let t = [1.0, -2.0, 0.5];
    for mesh in [mesh_of(icosphere(1)), tetrahedron()] {
        let mapped_then_refined = butterfly_subdivision(&affine(&mesh, m, t));
        let refined_then_mapped = affine(&butterfly_subdivision(&mesh), m, t);
        assert_eq!(mapped_then_refined.faces, refined_then_mapped.faces);
        for (index, (p, q)) in mapped_then_refined
            .vertices
            .iter()
            .zip(&refined_then_mapped.vertices)
            .enumerate()
        {
            assert_close(*p, *q, &format!("vertex {index}"));
        }
    }
}

#[test]
fn a_sphere_does_not_shrink_towards_the_origin() {
    // New vertices of an interpolating scheme on an icosphere lie close to the
    // unit sphere; a stencil whose weights sum to 7/8 would put them at a radius
    // below 0.9.
    let mesh = mesh_of(icosphere(2));
    let refined = butterfly_subdivision(&mesh);
    for (index, p) in refined.vertices.iter().enumerate().skip(mesh.n_vertices()) {
        let radius = len3(*p);
        assert!(
            radius > 0.98 && radius < 1.02,
            "new vertex {index} has radius {radius}"
        );
    }
}

// ── Interpolation and topology ───────────────────────────────────────────────

#[test]
fn original_vertices_are_interpolated_exactly() {
    for mesh in [
        mesh_of(icosphere(1)),
        tetrahedron(),
        mesh_of(flat_grid(3, 3, 1.0)),
    ] {
        let refined = butterfly_subdivision(&mesh);
        assert_eq!(&refined.vertices[..mesh.n_vertices()], &mesh.vertices[..]);
    }
}

#[test]
fn closed_meshes_stay_closed_and_consistently_oriented() {
    for mesh in [mesh_of(icosphere(0)), mesh_of(icosphere(1)), tetrahedron()] {
        assert_closed_and_oriented(&mesh);
        let refined = butterfly_subdivision(&mesh);
        assert_eq!(refined.n_faces(), 4 * mesh.n_faces());
        let edges = 3 * mesh.n_faces() / 2;
        assert_eq!(refined.n_vertices(), mesh.n_vertices() + edges);
        assert_closed_and_oriented(&refined);
    }
}

// ── Regular stencil ──────────────────────────────────────────────────────────

/// Key of an edge midpoint of the unit-spacing flat grid (multiples of one half).
fn midpoint_key(p: [f64; 3]) -> (i64, i64) {
    ((p[0] * 2.0).round() as i64, (p[1] * 2.0).round() as i64)
}

fn cubic_height(x: f64, y: f64) -> f64 {
    0.3 * x * x * x - 0.7 * x * x * y + 0.2 * x * y * y + 0.5 * y * y * y + 0.9 * x * x
        - 0.4 * x * y
        + 0.6 * y * y
        + 0.8 * x
        - 0.3 * y
        + 1.2
}

#[test]
fn flat_regular_grid_edges_land_on_their_midpoints() {
    let n = 8;
    let mesh = mesh_of(flat_grid(n, n, 1.0));
    let refined = butterfly_subdivision(&mesh);
    let mut checked = 0;
    for p in &refined.vertices[mesh.n_vertices()..] {
        assert!(p[2].abs() < TOLERANCE, "vertex {p:?} left the plane");
        let (kx, ky) = midpoint_key(*p);
        let on_boundary = kx == 0 || ky == 0 || kx == 2 * n as i64 || ky == 2 * n as i64;
        if on_boundary {
            // The four-point boundary rule overshoots at the corners of the
            // rectangle, so only interior edge points are exact midpoints.
            continue;
        }
        // The stencil is point-symmetric on a regular lattice: the new vertex is
        // exactly the midpoint of its edge.
        assert!(
            (p[0] * 2.0 - kx as f64).abs() < TOLERANCE
                && (p[1] * 2.0 - ky as f64).abs() < TOLERANCE,
            "vertex {p:?} is not on a half-integer lattice point"
        );
        checked += 1;
    }
    // The grid has n (n + 1) horizontal, n (n + 1) vertical and n^2 diagonal
    // edges, 4 n of which lie on the boundary.
    assert_eq!(checked, 3 * n * n - 2 * n);
}

#[test]
fn regular_stencil_reproduces_cubic_polynomials_in_the_interior() {
    let n = 9;
    let (mut vertices, faces) = flat_grid(n, n, 1.0);
    for v in &mut vertices {
        v[2] = cubic_height(v[0], v[1]);
    }
    let mesh = SubdivMesh::new(vertices, faces);
    let refined = butterfly_subdivision(&mesh);
    let mut checked = 0;
    for p in &refined.vertices[mesh.n_vertices()..] {
        // Midpoints with all eight stencil vertices inside the grid.
        if p[0] >= 2.0 && p[0] <= n as f64 - 2.0 && p[1] >= 2.0 && p[1] <= n as f64 - 2.0 {
            let expected = cubic_height(p[0], p[1]);
            assert!(
                (p[2] - expected).abs() < 1e-9,
                "height at ({}, {}) is {} instead of {expected}",
                p[0],
                p[1],
                p[2]
            );
            checked += 1;
        }
    }
    assert!(checked > 60, "only {checked} stencils were checked");
}

#[test]
fn regular_stencil_uses_the_documented_wing_vertices_for_every_edge_direction() {
    let n = 8;
    let (mut vertices, faces) = flat_grid(n, n, 1.0);
    // Irregular heights so that every vertex of the stencil matters.
    for (index, v) in vertices.iter_mut().enumerate() {
        v[2] = ((index * 37 + 11) % 17) as f64 * 0.25;
    }
    let mesh = SubdivMesh::new(vertices, faces);
    let refined = butterfly_subdivision(&mesh);
    let height = |i: i64, j: i64| mesh.vertices[(j as usize) * (n + 1) + i as usize][2];

    let mut by_midpoint: HashMap<(i64, i64), f64> = HashMap::new();
    for p in &refined.vertices[mesh.n_vertices()..] {
        by_midpoint.insert(midpoint_key(*p), p[2]);
    }

    // Lattice offsets of the stencil `a, b, c, d, w1, w2, w3, w4` of one edge of
    // each direction, relative to the lower-left end point `a`.
    let stencils: [[(i64, i64); 8]; 3] = [
        // horizontal edge (0,0)-(1,0)
        [
            (0, 0),
            (1, 0),
            (1, 1),
            (0, -1),
            (2, 1),
            (0, 1),
            (-1, -1),
            (1, -1),
        ],
        // vertical edge (0,0)-(0,1)
        [
            (0, 0),
            (0, 1),
            (1, 1),
            (-1, 0),
            (1, 2),
            (1, 0),
            (-1, -1),
            (-1, 1),
        ],
        // diagonal edge (0,0)-(1,1)
        [
            (0, 0),
            (1, 1),
            (1, 0),
            (0, 1),
            (2, 1),
            (0, -1),
            (-1, 0),
            (1, 2),
        ],
    ];
    let mut checked = 0;
    for stencil in &stencils {
        for i in 2..(n as i64 - 2) {
            for j in 2..(n as i64 - 2) {
                let at = |k: usize| height(i + stencil[k].0, j + stencil[k].1);
                let expected = 0.5 * (at(0) + at(1)) + 0.125 * (at(2) + at(3))
                    - 0.0625 * (at(4) + at(5) + at(6) + at(7));
                let key = (
                    2 * i + stencil[0].0 + stencil[1].0,
                    2 * j + stencil[0].1 + stencil[1].1,
                );
                let actual = by_midpoint[&key];
                assert!(
                    (actual - expected).abs() < TOLERANCE,
                    "stencil {stencil:?} at ({i}, {j}): {actual} vs {expected}"
                );
                checked += 1;
            }
        }
    }
    assert_eq!(checked, 3 * 4 * 4);
}

// ── Extraordinary vertices ───────────────────────────────────────────────────

#[test]
fn tetrahedron_edge_points_use_the_valence_three_stencil() {
    let mesh = tetrahedron();
    let refined = butterfly_subdivision(&mesh);
    // Both end points have valence 3, so the result is the average of the two
    // stencils 3/4 a + 5/12 b - 1/12 (c + d) and 3/4 b + 5/12 a - 1/12 (c + d),
    // i.e. 7/12 (a + b) - 1/12 (c + d).
    let mut seen = HashSet::new();
    for face in &mesh.faces {
        for k in 0..3 {
            let (a, b) = (face[k], face[(k + 1) % 3]);
            if !seen.insert((a.min(b), a.max(b))) {
                continue;
            }
            let others: Vec<usize> = (0..4).filter(|&v| v != a && v != b).collect();
            let ends = scale3(add3(mesh.vertices[a], mesh.vertices[b]), 7.0 / 12.0);
            let opposite = scale3(
                add3(mesh.vertices[others[0]], mesh.vertices[others[1]]),
                -1.0 / 12.0,
            );
            let expected = add3(ends, opposite);
            assert!(
                contains_point(&refined.vertices[mesh.n_vertices()..], expected, TOLERANCE),
                "no refined vertex at {expected:?} for edge ({a}, {b})"
            );
        }
    }
    assert_eq!(seen.len(), 6);
    // The edge midpoints of this tetrahedron have radius 1; the new vertices are
    // pushed outwards (to 4/3) towards the circumscribed sphere (radius sqrt(3))
    // instead of collapsing onto the midpoints.
    for p in &refined.vertices[mesh.n_vertices()..] {
        let radius = len3(*p);
        assert!((radius - 4.0 / 3.0).abs() < TOLERANCE, "radius {radius}");
    }
}

/// Edge point of the unit icosahedron predicted by the modified butterfly rule
/// for the end point `a` (valence 5, ring ordered by the angle around the radial
/// direction and starting at the other end point `b`).
fn icosahedron_stencil_from(vertices: &[[f64; 3]], a: usize, b: usize) -> [f64; 3] {
    let normal = vertices[a];
    let mut ring: Vec<usize> = (0..vertices.len())
        .filter(|&q| q != a && len3(sub3(vertices[q], vertices[a])) < 1.2)
        .collect();
    assert_eq!(ring.len(), 5);
    let tangent = {
        let first = sub3(vertices[b], vertices[a]);
        let along = dot3(first, normal);
        let t = sub3(first, scale3(normal, along));
        scale3(t, 1.0 / len3(t))
    };
    let bitangent = [
        normal[1] * tangent[2] - normal[2] * tangent[1],
        normal[2] * tangent[0] - normal[0] * tangent[2],
        normal[0] * tangent[1] - normal[1] * tangent[0],
    ];
    let angle_of = |q: usize| {
        let d = sub3(vertices[q], vertices[a]);
        let angle = dot3(d, bitangent).atan2(dot3(d, tangent));
        if angle < -1e-9 {
            angle + 2.0 * std::f64::consts::PI
        } else {
            angle
        }
    };
    ring.sort_by(|&p, &q| angle_of(p).total_cmp(&angle_of(q)));
    assert_eq!(ring[0], b);
    let mut point = scale3(vertices[a], 0.75);
    for (j, &q) in ring.iter().enumerate() {
        let angle = 2.0 * std::f64::consts::PI * j as f64 / 5.0;
        let weight = (0.25 + angle.cos() + 0.5 * (2.0 * angle).cos()) / 5.0;
        point = add3(point, scale3(vertices[q], weight));
    }
    point
}

#[test]
fn icosahedron_edge_points_follow_the_valence_five_modified_stencil() {
    let mesh = mesh_of(icosphere(0));
    let refined = butterfly_subdivision(&mesh);
    let mut seen = HashSet::new();
    for face in &mesh.faces {
        for k in 0..3 {
            let (a, b) = (face[k], face[(k + 1) % 3]);
            if !seen.insert((a.min(b), a.max(b))) {
                continue;
            }
            // Both end points are extraordinary: the two stencils are averaged.
            let from_a = icosahedron_stencil_from(&mesh.vertices, a, b);
            let from_b = icosahedron_stencil_from(&mesh.vertices, b, a);
            let expected = scale3(add3(from_a, from_b), 0.5);
            assert!(
                contains_point(&refined.vertices[mesh.n_vertices()..], expected, 1e-10),
                "no refined vertex at {expected:?} for edge ({a}, {b})"
            );
            // The edge midpoints of the icosahedron have radius 0.8507; the
            // stencil moves the new vertices outwards, closer to the sphere.
            let radius = len3(expected);
            assert!(radius > 0.86 && radius < 1.0, "radius {radius}");
        }
    }
    assert_eq!(seen.len(), 30);
}

/// A planar fan of `valence` triangles around the interior vertex 0 (the
/// centre); the other vertices lie on a circle of the given radius at equal
/// angles and are boundary vertices.
fn fan(valence: usize, radius: f64, centre: [f64; 3]) -> SubdivMesh {
    let mut vertices = vec![centre];
    for j in 0..valence {
        let angle = 2.0 * std::f64::consts::PI * j as f64 / valence as f64;
        vertices.push([
            centre[0] + radius * angle.cos(),
            centre[1] + radius * angle.sin(),
            centre[2],
        ]);
    }
    let faces = (0..valence)
        .map(|j| [0, 1 + j, 1 + (j + 1) % valence])
        .collect();
    SubdivMesh::new(vertices, faces)
}

#[test]
fn extraordinary_stencils_have_linear_precision_on_a_planar_fan() {
    // The edge from the centre to the first ring vertex touches one extraordinary
    // vertex (the centre) and one boundary vertex, so exactly the valence-K
    // stencil 3/4 c + sum s_j q_j is used. On an equiangular planar fan it must
    // reproduce the midpoint of the edge: the weights have to sum to one and
    // their first moment sum s_j cos(2 pi j / K) has to be 1/2. This checks the
    // weights without restating them. (Valence 6 has no extraordinary stencil and
    // falls back to the midpoint, as the missing wing vertices leave no choice.)
    for valence in [3, 4, 5, 6, 7, 8, 9, 12] {
        for (radius, centre) in [(1.0, [0.0; 3]), (2.5, [0.3, -0.2, 0.1])] {
            let mesh = fan(valence, radius, centre);
            let refined = butterfly_subdivision(&mesh);
            // The first face is (0, 1, 2); its first edge (0, 1) gets the first
            // new vertex.
            let expected = midpoint3(mesh.vertices[0], mesh.vertices[1]);
            assert_close(
                refined.vertices[mesh.n_vertices()],
                expected,
                &format!("valence {valence}, radius {radius}"),
            );
        }
    }
}

#[test]
fn octahedron_edge_points_use_the_valence_four_stencil() {
    // Every vertex has valence 4 and its ring alternates between the other end
    // point b of the edge and vertices that carry weight 0 or -1/8:
    // 3/4 a + 3/8 b + 0 - 1/8 (-b) + 0 = 3/4 a + 1/2 b. Averaging the stencils of
    // both end points gives 5/8 (a + b) for every edge.
    let vertices = vec![
        [1.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, -1.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 0.0, -1.0],
    ];
    let faces = vec![
        [0, 2, 4],
        [2, 1, 4],
        [1, 3, 4],
        [3, 0, 4],
        [2, 0, 5],
        [1, 2, 5],
        [3, 1, 5],
        [0, 3, 5],
    ];
    let mesh = SubdivMesh::new(vertices, faces);
    assert_closed_and_oriented(&mesh);
    let refined = butterfly_subdivision(&mesh);
    let new_vertices = &refined.vertices[mesh.n_vertices()..];
    assert_eq!(new_vertices.len(), 12);
    let mut seen = HashSet::new();
    for face in &mesh.faces {
        for k in 0..3 {
            let (a, b) = (face[k], face[(k + 1) % 3]);
            if !seen.insert((a.min(b), a.max(b))) {
                continue;
            }
            let expected = scale3(add3(mesh.vertices[a], mesh.vertices[b]), 5.0 / 8.0);
            assert!(
                contains_point(new_vertices, expected, TOLERANCE),
                "no refined vertex at {expected:?} for edge ({a}, {b})"
            );
        }
    }
    assert_eq!(seen.len(), 12);
}

// ── Boundaries ───────────────────────────────────────────────────────────────

#[test]
fn boundary_edges_use_the_four_point_rule() {
    let (mut vertices, faces) = flat_grid(4, 1, 1.0);
    // Bottom boundary heights 0, 0, 1, 0, 0 at x = 0..4.
    vertices[2][2] = 1.0;
    let mesh = SubdivMesh::new(vertices, faces);
    let refined = butterfly_subdivision(&mesh);
    let new_vertices = &refined.vertices[mesh.n_vertices()..];
    // Edge (1,0)-(2,0): 9/16 (z1 + z2) - 1/16 (z0 + z3) = 9/16, and the plain
    // midpoint x = 1.5 because the spacing is uniform.
    assert!(contains_point(
        new_vertices,
        [1.5, 0.0, 9.0 / 16.0],
        TOLERANCE
    ));
    // Edge (2,0)-(3,0) mirrors it.
    assert!(contains_point(
        new_vertices,
        [2.5, 0.0, 9.0 / 16.0],
        TOLERANCE
    ));
}

#[test]
fn a_single_triangle_uses_its_closed_three_gon_boundary() {
    let mesh = SubdivMesh::new(
        vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        vec![[0, 1, 2]],
    );
    let refined = butterfly_subdivision(&mesh);
    assert_eq!(refined.n_faces(), 4);
    assert_eq!(refined.n_vertices(), 6);
    // The four-point rule on the closed boundary 3-gon reads
    // 9/16 (a + b) - 1/16 (c + c), whose weights sum to one.
    assert_close(
        refined.vertices[3],
        [9.0 / 16.0, -2.0 / 16.0, 0.0],
        "edge (0, 1)",
    );
    for p in &refined.vertices {
        assert!(p.iter().all(|c| c.is_finite()));
    }
}

#[test]
fn edges_without_a_stencil_fall_back_to_the_midpoint() {
    // Three triangles share the edge (0, 1): the edge is non-manifold.
    let mesh = SubdivMesh::new(
        vec![
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [1.0, -1.0, 0.0],
            [1.0, 0.0, 1.0],
        ],
        vec![[0, 1, 2], [1, 0, 3], [0, 1, 4]],
    );
    let refined = butterfly_subdivision(&mesh);
    assert_eq!(refined.n_faces(), 12);
    assert!(contains_point(
        &refined.vertices[mesh.n_vertices()..],
        [1.0, 0.0, 0.0],
        TOLERANCE
    ));
    for face in &refined.faces {
        assert!(face.iter().all(|&v| v < refined.n_vertices()));
    }
    for p in &refined.vertices {
        assert!(p.iter().all(|c| c.is_finite()));
    }
}

#[test]
fn interior_edge_next_to_the_boundary_without_wings_uses_the_midpoint() {
    // Two triangles sharing one edge: every other edge is a boundary edge and no
    // wing vertex exists for the shared edge.
    let mesh = SubdivMesh::new(
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, 1.0, 0.3],
            [1.5, 1.0, -0.2],
        ],
        vec![[0, 1, 2], [1, 3, 2]],
    );
    let refined = butterfly_subdivision(&mesh);
    // Shared edge (1, 2).
    let midpoint = midpoint3(mesh.vertices[1], mesh.vertices[2]);
    assert!(contains_point(
        &refined.vertices[mesh.n_vertices()..],
        midpoint,
        TOLERANCE
    ));
}

// ── Orientation independence ─────────────────────────────────────────────────

#[test]
fn reversing_every_face_gives_the_same_vertex_set() {
    let mesh = mesh_of(icosphere(1));
    let reversed = SubdivMesh::new(
        mesh.vertices.clone(),
        mesh.faces.iter().map(|&[a, b, c]| [a, c, b]).collect(),
    );
    let refined = butterfly_subdivision(&mesh);
    let refined_reversed = butterfly_subdivision(&reversed);
    assert_eq!(refined.n_vertices(), refined_reversed.n_vertices());
    for (index, p) in refined.vertices.iter().enumerate() {
        assert!(
            contains_point(&refined_reversed.vertices, *p, 1e-12),
            "vertex {index} = {p:?} is missing after reversing the faces"
        );
    }
}

#[test]
fn butterfly_is_deterministic() {
    let mesh = mesh_of(icosphere(1));
    let first = butterfly_subdivision(&mesh);
    let second = butterfly_subdivision(&mesh);
    assert_eq!(first.vertices, second.vertices);
    assert_eq!(first.faces, second.faces);
}
