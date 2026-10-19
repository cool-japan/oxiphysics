//! Regression tests for Catmull–Clark subdivision.
//!
//! An earlier version applied the interior vertex rule to every vertex, which
//! mapped the corner `(0, 0)` of a unit quad onto the opposite corner `(1, 1)` and
//! pulled open boundaries off their boundary curve, and it numbered the edge
//! points in the iteration order of a `HashMap`, so two runs on the same mesh
//! produced different vertex arrays.

use super::*;
use std::collections::HashSet;

const TOLERANCE: f64 = 1e-12;

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

fn contains_point(vertices: &[[f64; 3]], target: [f64; 3]) -> bool {
    vertices.iter().any(|p| len3(sub3(*p, target)) < TOLERANCE)
}

fn unit_quad() -> QuadMesh {
    QuadMesh::new(
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ],
        vec![[0, 1, 2, 3]],
    )
}

/// `nx x ny` grid of quads (vertex `(i, j)` has index `j * (nx + 1) + i`) with
/// irregular heights, counter-clockwise seen from above.
fn bumpy_grid(nx: usize, ny: usize) -> QuadMesh {
    let mut vertices = Vec::new();
    for j in 0..=ny {
        for i in 0..=nx {
            let index = j * (nx + 1) + i;
            vertices.push([i as f64, j as f64, ((index * 37 + 11) % 17) as f64 * 0.25]);
        }
    }
    let at = |i: usize, j: usize| j * (nx + 1) + i;
    let mut faces = Vec::new();
    for j in 0..ny {
        for i in 0..nx {
            faces.push([at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1)]);
        }
    }
    QuadMesh::new(vertices, faces)
}

/// Cube `[-1, 1]^3` with outward counter-clockwise quads. Vertex `i` has the
/// coordinates `(+-1, +-1, +-1)` with the bits of `i` selecting the signs of
/// `x`, `y`, `z`.
fn cube() -> QuadMesh {
    let vertices = (0..8)
        .map(|i: usize| {
            [
                if i & 1 == 0 { -1.0 } else { 1.0 },
                if i & 2 == 0 { -1.0 } else { 1.0 },
                if i & 4 == 0 { -1.0 } else { 1.0 },
            ]
        })
        .collect();
    let faces = vec![
        [0, 2, 3, 1], // z = -1
        [4, 5, 7, 6], // z = +1
        [0, 1, 5, 4], // y = -1
        [2, 6, 7, 3], // y = +1
        [0, 4, 6, 2], // x = -1
        [1, 3, 7, 5], // x = +1
    ];
    QuadMesh::new(vertices, faces)
}

// ── Corners and boundaries ───────────────────────────────────────────────────

#[test]
fn single_quad_keeps_its_corners() {
    let mesh = unit_quad();
    let refined = catmull_clark_subdivision(&mesh);
    assert_eq!(refined.n_vertices(), 9);
    assert_eq!(refined.n_faces(), 4);
    // The interior rule would have mapped the corner (0, 0) to (1, 1) and
    // exchanged opposite corners.
    for (index, corner) in mesh.vertices.iter().enumerate() {
        assert_close(refined.vertices[index], *corner, &format!("corner {index}"));
    }
}

#[test]
fn single_quad_has_the_expected_face_point_edge_points_and_faces() {
    let refined = catmull_clark_subdivision(&unit_quad());
    // Face point, then the edge points in the order in which the sides are met.
    assert_close(refined.vertices[4], [0.5, 0.5, 0.0], "face point");
    assert_close(refined.vertices[5], [0.5, 0.0, 0.0], "edge (0, 1)");
    assert_close(refined.vertices[6], [1.0, 0.5, 0.0], "edge (1, 2)");
    assert_close(refined.vertices[7], [0.5, 1.0, 0.0], "edge (2, 3)");
    assert_close(refined.vertices[8], [0.0, 0.5, 0.0], "edge (3, 0)");
    assert_eq!(
        refined.faces,
        vec![[0, 5, 4, 8], [1, 6, 4, 5], [2, 7, 4, 6], [3, 8, 4, 7]]
    );
}

#[test]
fn boundary_vertices_follow_the_cubic_bspline_boundary_rule() {
    let mesh = bumpy_grid(4, 4);
    let refined = catmull_clark_subdivision(&mesh);
    let at = |i: usize, j: usize| j * 5 + i;
    let p = |i: usize, j: usize| mesh.vertices[at(i, j)];
    let boundary_rule = |centre: [f64; 3], first: [f64; 3], second: [f64; 3]| {
        add3(scale3(centre, 0.75), scale3(add3(first, second), 0.125))
    };
    // Bottom, top, left and right boundary vertices (not corners).
    for k in 1..4 {
        assert_close(
            refined.vertices[at(k, 0)],
            boundary_rule(p(k, 0), p(k - 1, 0), p(k + 1, 0)),
            &format!("bottom {k}"),
        );
        assert_close(
            refined.vertices[at(k, 4)],
            boundary_rule(p(k, 4), p(k - 1, 4), p(k + 1, 4)),
            &format!("top {k}"),
        );
        assert_close(
            refined.vertices[at(0, k)],
            boundary_rule(p(0, k), p(0, k - 1), p(0, k + 1)),
            &format!("left {k}"),
        );
        assert_close(
            refined.vertices[at(4, k)],
            boundary_rule(p(4, k), p(4, k - 1), p(4, k + 1)),
            &format!("right {k}"),
        );
    }
    // The four corners of the grid stay where they are.
    for (i, j) in [(0, 0), (4, 0), (0, 4), (4, 4)] {
        assert_close(
            refined.vertices[at(i, j)],
            p(i, j),
            &format!("corner ({i}, {j})"),
        );
    }
}

#[test]
fn regular_interior_vertices_use_the_bicubic_bspline_mask() {
    let mesh = bumpy_grid(4, 4);
    let refined = catmull_clark_subdivision(&mesh);
    let at = |i: usize, j: usize| j * 5 + i;
    let p = |i: usize, j: usize| mesh.vertices[at(i, j)];
    for j in 1..4 {
        for i in 1..4 {
            // (36 P + 6 sum of the 4 edge neighbours + sum of the 4 diagonals) / 64
            let mut expected = scale3(p(i, j), 36.0 / 64.0);
            for (di, dj) in [(-1_i32, 0_i32), (1, 0), (0, -1), (0, 1)] {
                let q = p((i as i32 + di) as usize, (j as i32 + dj) as usize);
                expected = add3(expected, scale3(q, 6.0 / 64.0));
            }
            for (di, dj) in [(-1_i32, -1_i32), (1, -1), (-1, 1), (1, 1)] {
                let q = p((i as i32 + di) as usize, (j as i32 + dj) as usize);
                expected = add3(expected, scale3(q, 1.0 / 64.0));
            }
            assert_close(refined.vertices[at(i, j)], expected, &format!("({i}, {j})"));
        }
    }
}

#[test]
fn interior_and_boundary_edge_points_follow_their_rules() {
    let mesh = bumpy_grid(4, 4);
    let refined = catmull_clark_subdivision(&mesh);
    let at = |i: usize, j: usize| j * 5 + i;
    let p = |i: usize, j: usize| mesh.vertices[at(i, j)];
    let face_point = |i: usize, j: usize| {
        let sum = [p(i, j), p(i + 1, j), p(i + 1, j + 1), p(i, j + 1)]
            .iter()
            .fold([0.0; 3], |acc, &q| add3(acc, q));
        scale3(sum, 0.25)
    };
    // Interior edge (1,1)-(2,1) lies between the cells (1,0) and (1,1).
    let interior = scale3(
        add3(
            add3(p(1, 1), p(2, 1)),
            add3(face_point(1, 0), face_point(1, 1)),
        ),
        0.25,
    );
    assert!(
        contains_point(&refined.vertices, interior),
        "interior edge point"
    );
    // Boundary edge (1,0)-(2,0) uses the midpoint.
    let boundary = midpoint3(p(1, 0), p(2, 0));
    assert!(
        contains_point(&refined.vertices, boundary),
        "boundary edge point"
    );
    // Face points are the face centroids.
    let faces_start = mesh.n_vertices();
    assert_close(
        refined.vertices[faces_start + 5],
        face_point(1, 1),
        "face point of cell (1, 1)",
    );
}

// ── Closed meshes ────────────────────────────────────────────────────────────

#[test]
fn cube_corners_move_to_five_ninths() {
    let mesh = cube();
    let refined = catmull_clark_subdivision(&mesh);
    // Valence 3: F = (1/3, 1/3, 1/3), R = (2/3, 2/3, 2/3), (F + 2 R) / 3 = 5/9.
    for (index, corner) in mesh.vertices.iter().enumerate() {
        let expected = scale3(*corner, 5.0 / 9.0);
        assert_close(
            refined.vertices[index],
            expected,
            &format!("corner {index}"),
        );
    }
    // Face points are the face centres.
    let centres = [
        [0.0, 0.0, -1.0],
        [0.0, 0.0, 1.0],
        [0.0, -1.0, 0.0],
        [0.0, 1.0, 0.0],
        [-1.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
    ];
    for (index, centre) in centres.iter().enumerate() {
        assert_close(
            refined.vertices[8 + index],
            *centre,
            &format!("face {index}"),
        );
    }
    // Edge points (a + b + f1 + f2) / 4: two coordinates +-3/4 and one zero.
    for point in &refined.vertices[14..] {
        let zeros = point.iter().filter(|c| c.abs() < TOLERANCE).count();
        let three_quarters = point
            .iter()
            .filter(|c| (c.abs() - 0.75).abs() < TOLERANCE)
            .count();
        assert_eq!((zeros, three_quarters), (1, 2), "edge point {point:?}");
    }
}

#[test]
fn closed_input_gives_a_closed_consistently_oriented_output() {
    let refined = catmull_clark_subdivision(&cube());
    assert_eq!(refined.n_vertices(), 8 + 6 + 12);
    assert_eq!(refined.n_faces(), 24);
    let mut directed: HashSet<(usize, usize)> = HashSet::new();
    for face in &refined.faces {
        for k in 0..4 {
            let edge = (face[k], face[(k + 1) % 4]);
            assert!(
                directed.insert(edge),
                "directed edge {edge:?} appears twice"
            );
        }
    }
    for &(a, b) in &directed {
        assert!(
            directed.contains(&(b, a)),
            "edge ({a}, {b}) has no opposite"
        );
    }
    // Euler characteristic of a sphere: V - E + F = 2.
    let edges = directed.len() / 2;
    assert_eq!(refined.n_vertices() + refined.n_faces() - edges, 2);
}

#[test]
fn repeated_subdivision_of_the_cube_shrinks_towards_the_limit_surface() {
    let mut mesh = cube();
    let mut previous_extent = f64::INFINITY;
    for _ in 0..3 {
        mesh = catmull_clark_subdivision(&mesh);
        let extent = mesh
            .vertices
            .iter()
            .map(|p| p[0].abs().max(p[1].abs()).max(p[2].abs()))
            .fold(0.0, f64::max);
        assert!(extent < previous_extent, "extent {extent}");
        assert!(extent > 0.5);
        previous_extent = extent;
    }
    // Eight corner vertices of the original cube stay on the diagonals by symmetry.
    for p in &mesh.vertices[..8] {
        assert!((p[0].abs() - p[1].abs()).abs() < 1e-12);
        assert!((p[1].abs() - p[2].abs()).abs() < 1e-12);
    }
}

// ── Determinism ──────────────────────────────────────────────────────────────

#[test]
fn two_runs_give_identical_output() {
    for mesh in [cube(), bumpy_grid(4, 3), unit_quad()] {
        let reference = catmull_clark_subdivision(&mesh);
        for _ in 0..20 {
            let again = catmull_clark_subdivision(&mesh);
            assert_eq!(again.vertices, reference.vertices);
            assert_eq!(again.faces, reference.faces);
        }
    }
}

#[test]
fn edge_points_are_numbered_in_face_order() {
    // Two quads sharing the edge (1, 4): the shared edge is numbered when the
    // first quad is traversed.
    let mesh = QuadMesh::new(
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
            [2.0, 1.0, 0.0],
        ],
        vec![[0, 1, 4, 3], [1, 2, 5, 4]],
    );
    let refined = catmull_clark_subdivision(&mesh);
    let edge_start = 6 + 2;
    // First quad: (0,1), (1,4), (4,3), (3,0); second quad adds (1,2), (2,5), (5,4).
    let expected_midpoints = [
        [0.5, 0.0, 0.0],
        [1.0, 0.5, 0.0], // shared interior edge: (a + b + f1 + f2) / 4 with f = (0.5,0.5), (1.5,0.5)
        [0.5, 1.0, 0.0],
        [0.0, 0.5, 0.0],
        [1.5, 0.0, 0.0],
        [2.0, 0.5, 0.0],
        [1.5, 1.0, 0.0],
    ];
    assert_eq!(refined.n_vertices(), edge_start + 7);
    for (offset, expected) in expected_midpoints.iter().enumerate() {
        assert_close(
            refined.vertices[edge_start + offset],
            *expected,
            &format!("edge point {offset}"),
        );
    }
}

// ── Unusual input ────────────────────────────────────────────────────────────

#[test]
fn a_non_manifold_edge_uses_midpoints_and_fixed_end_points() {
    // Three quads share the edge (0, 1).
    let mesh = QuadMesh::new(
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, -1.0, 0.0],
            [0.0, -1.0, 0.0],
            [1.0, 0.0, 1.0],
            [0.0, 0.0, 1.0],
        ],
        vec![[0, 1, 2, 3], [1, 0, 5, 4], [0, 1, 6, 7]],
    );
    let refined = catmull_clark_subdivision(&mesh);
    assert_eq!(refined.n_faces(), 12);
    assert_close(refined.vertices[0], mesh.vertices[0], "end point 0");
    assert_close(refined.vertices[1], mesh.vertices[1], "end point 1");
    assert!(contains_point(&refined.vertices, [0.5, 0.0, 0.0]));
    for p in &refined.vertices {
        assert!(p.iter().all(|c| c.is_finite()));
    }
}

#[test]
fn vertices_that_belong_to_no_face_are_kept() {
    let mut mesh = unit_quad();
    mesh.vertices.push([5.0, 5.0, 5.0]);
    let refined = catmull_clark_subdivision(&mesh);
    assert_close(refined.vertices[4], [5.0, 5.0, 5.0], "isolated vertex");
    // The remaining vertices follow the single-quad layout shifted by one slot.
    assert_close(refined.vertices[5], [0.5, 0.5, 0.0], "face point");
    assert_eq!(refined.n_vertices(), 5 + 1 + 4);
}

#[test]
fn empty_mesh_stays_empty() {
    let refined = catmull_clark_subdivision(&QuadMesh::new(Vec::new(), Vec::new()));
    assert_eq!(refined.n_vertices(), 0);
    assert_eq!(refined.n_faces(), 0);
}
