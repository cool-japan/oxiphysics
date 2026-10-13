// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Unit tests of the convex-decomposition helpers: vector helpers, the
//! half-edge convex hull, support functions, Minkowski sums, GJK distance,
//! containment, SAT intersection, hull inertia and hull erosion.

use super::*;

fn approx_eq(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() < tol
}

/// Build a simple octahedron point set (convex by construction).
fn octahedron_points() -> Vec<Vec3> {
    vec![
        [1.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, -1.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 0.0, -1.0],
    ]
}

/// Build a cube point set.
fn cube_points(half: f64) -> Vec<Vec3> {
    let h = half;
    vec![
        [-h, -h, -h],
        [h, -h, -h],
        [-h, h, -h],
        [h, h, -h],
        [-h, -h, h],
        [h, -h, h],
        [-h, h, h],
        [h, h, h],
    ]
}

// ── Vec3 helpers ──

#[test]
fn test_dot_product() {
    assert!(approx_eq(
        dot([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]),
        32.0,
        1e-12
    ));
}

#[test]
fn test_cross_product_orthogonal() {
    let c = cross([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
    assert!(approx_eq(c[0], 0.0, 1e-12));
    assert!(approx_eq(c[1], 0.0, 1e-12));
    assert!(approx_eq(c[2], 1.0, 1e-12));
}

#[test]
fn test_normalize_unit_length() {
    let v = normalize([3.0, 4.0, 0.0]);
    assert!(approx_eq(norm(v), 1.0, 1e-12), "norm={:.6}", norm(v));
}

#[test]
fn test_norm_zero_vector() {
    let v = normalize([0.0, 0.0, 0.0]);
    assert!(approx_eq(norm(v), 0.0, 1e-12));
}

// ── ConvexHull3d ──

#[test]
fn test_hull_builds_from_octahedron() {
    let pts = octahedron_points();
    let hull = ConvexHull3d::build(&pts);
    assert!(!hull.triangles.is_empty(), "hull should have triangles");
}

#[test]
fn test_hull_all_points_inside_or_on() {
    let pts = octahedron_points();
    let hull = ConvexHull3d::build(&pts);
    assert!(
        hull.all_points_inside_or_on(&pts, 1e-9),
        "all input points should be on or inside hull"
    );
}

#[test]
fn test_hull_volume_positive() {
    let pts = octahedron_points();
    let hull = ConvexHull3d::build(&pts);
    let vol = hull.volume();
    assert!(vol > 0.0, "hull volume should be positive: {:.6}", vol);
}

#[test]
fn test_hull_volume_cube() {
    // Cube with side 2: volume 8.  The first four points of
    // `cube_points` are coplanar, which used to break the seed
    // tetrahedron and produce a garbage hull.
    let pts = cube_points(1.0);
    let hull = ConvexHull3d::build(&pts);
    let vol = hull.volume();
    assert!(approx_eq(vol, 8.0, 1e-9), "cube hull volume = {vol}");
    assert_eq!(hull.triangles.len(), 12);
    assert!(hull.all_points_inside_or_on(&pts, 1e-9));
}

/// Deterministic permutation of `items` (xorshift Fisher-Yates).
fn permuted<T: Clone>(items: &[T], seed: u64) -> Vec<T> {
    let mut state = seed | 1;
    let mut out = items.to_vec();
    for i in (1..out.len()).rev() {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        out.swap(i, (state % (i as u64 + 1)) as usize);
    }
    out
}

/// The hull surface must be a closed, consistently wound 2-manifold with
/// outward facing facets.
fn assert_closed_outward_surface(hull: &ConvexHull3d, label: &str) {
    for (i, edge) in hull.half_edges.iter().enumerate() {
        assert_ne!(edge.twin, usize::MAX, "{label}: half-edge {i} has no twin");
        let twin = &hull.half_edges[edge.twin];
        assert_eq!(twin.twin, i, "{label}: twins are not mutual");
        assert_eq!(
            twin.origin, hull.half_edges[edge.next].origin,
            "{label}: twin does not run backwards"
        );
    }
    let mut used: Vec<usize> = hull.triangles.iter().flatten().copied().collect();
    used.sort_unstable();
    used.dedup();
    let (v, e, f) = (
        used.len() as i64,
        (hull.half_edges.len() / 2) as i64,
        hull.triangles.len() as i64,
    );
    assert_eq!(
        v - e + f,
        2,
        "{label}: Euler characteristic V={v} E={e} F={f}"
    );
    let inside = hull.centroid();
    for (fi, face) in hull.faces.iter().enumerate() {
        let corner = hull.vertices[hull.triangles[fi][0]];
        assert!(
            dot(face.normal, sub(inside, corner)) < 0.0,
            "{label}: facet {fi} looks inward"
        );
    }
}

#[test]
fn test_hull_of_cube_and_octahedron_is_exact_for_every_point_order() {
    for seed in 0..24 {
        let cube = permuted(&cube_points(1.0), seed);
        let hull = ConvexHull3d::build(&cube);
        let label = format!("cube, order {seed}");
        assert_closed_outward_surface(&hull, &label);
        assert_eq!(hull.triangles.len(), 12, "{label}");
        assert!(approx_eq(hull.volume(), 8.0, 1e-9), "{label}");
        assert!(hull.all_points_inside_or_on(&cube, 1e-9), "{label}");
        let octahedron = permuted(&octahedron_points(), seed);
        let hull = ConvexHull3d::build(&octahedron);
        let label = format!("octahedron, order {seed}");
        assert_closed_outward_surface(&hull, &label);
        assert_eq!(hull.triangles.len(), 8, "{label}");
        assert!(approx_eq(hull.volume(), 4.0 / 3.0, 1e-12), "{label}");
    }
}

#[test]
fn test_hull_ignores_interior_face_and_edge_points() {
    let mut pts = cube_points(1.0);
    pts.push([0.0, 0.0, 0.0]); // interior
    pts.push([0.3, -0.2, 0.9]); // interior
    pts.push([0.0, 0.0, 1.0]); // centre of a face
    pts.push([1.0, 1.0, 0.0]); // middle of an edge
    for seed in 0..12 {
        let order = permuted(&pts, seed);
        let hull = ConvexHull3d::build(&order);
        let label = format!("cube with extra points, order {seed}");
        assert_closed_outward_surface(&hull, &label);
        assert!(approx_eq(hull.volume(), 8.0, 1e-9), "{label}");
        assert!(hull.all_points_inside_or_on(&order, 1e-9), "{label}");
    }
}

#[test]
fn test_hull_is_scale_independent() {
    for &factor in &[1e-9, 1e-4, 1.0, 1e6] {
        let pts: Vec<Vec3> = cube_points(1.0).iter().map(|&p| scale(p, factor)).collect();
        let hull = ConvexHull3d::build(&pts);
        assert_eq!(hull.triangles.len(), 12, "scale {factor}");
        let want = 8.0 * factor.powi(3);
        assert!(
            (hull.volume() - want).abs() <= 1e-9 * want,
            "scale {factor}: volume {} vs {want}",
            hull.volume()
        );
        assert_closed_outward_surface(&hull, &format!("scale {factor}"));
    }
}

#[test]
fn test_hull_of_degenerate_point_sets_has_no_facets() {
    let line: Vec<Vec3> = (0..6)
        .map(|i| [i as f64, 2.0 * i as f64, -(i as f64)])
        .collect();
    let cases: Vec<(&str, Vec<Vec3>)> = vec![
        ("empty", Vec::new()),
        ("single point", vec![[1.0, 2.0, 3.0]]),
        (
            "three points",
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        ),
        (
            "coplanar square",
            vec![
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [1.0, 1.0, 0.0],
                [0.0, 1.0, 0.0],
            ],
        ),
        ("collinear", line),
        ("coincident", vec![[1.0, 2.0, 3.0]; 5]),
        (
            "non-finite",
            vec![
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, f64::NAN],
            ],
        ),
    ];
    for (label, pts) in cases {
        // Used to panic for fewer than four points.
        let hull = ConvexHull3d::build(&pts);
        assert!(hull.triangles.is_empty(), "{label}");
        assert!(hull.faces.is_empty(), "{label}");
        assert!(hull.half_edges.is_empty(), "{label}");
        assert_eq!(hull.volume(), 0.0, "{label}");
        assert_eq!(hull.vertices.len(), pts.len(), "{label}");
        assert!(hull.all_points_inside_or_on(&pts, 1e-9), "{label}");
        assert!(hull.all_points_on_boundary(&pts, 1e-9), "{label}");
    }
}

#[test]
fn test_convex_part_detects_buried_vertices() {
    assert!(ConvexPart::new(cube_points(1.0)).is_convex(1e-9));
    // A vertex strictly inside the cube is not in convex position.
    let mut buried = cube_points(1.0);
    buried.push([0.1, 0.2, -0.3]);
    assert!(!ConvexPart::new(buried).is_convex(1e-9));
    // A vertex on a face or an edge still is.
    let mut on_boundary = cube_points(1.0);
    on_boundary.push([0.0, 0.0, 1.0]);
    on_boundary.push([1.0, 1.0, 0.0]);
    assert!(ConvexPart::new(on_boundary).is_convex(1e-9));
    // Fewer than four vertices, or vertices that span no volume, count
    // as convex.
    assert!(ConvexPart::new(vec![[0.0; 3], [1.0, 0.0, 0.0]]).is_convex(1e-9));
    let flat = vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.5, 0.5, 0.0],
    ];
    assert!(ConvexPart::new(flat).is_convex(1e-9));
}

#[test]
fn test_hull_centroid_near_origin_for_octahedron() {
    let pts = octahedron_points();
    let hull = ConvexHull3d::build(&pts);
    let c = hull.centroid();
    assert!(
        approx_eq(norm(c), 0.0, 1e-12),
        "centroid of symmetric octahedron should be near origin: {:.6}",
        norm(c)
    );
}

#[test]
fn test_support_function_maximum() {
    let pts = octahedron_points();
    let hull = ConvexHull3d::build(&pts);
    let dir = normalize([1.0, 0.0, 0.0]);
    let h = hull.support(dir);
    // Max dot product along X is 1.0 for unit octahedron
    assert!(
        approx_eq(h, 1.0, 1e-9),
        "support along x should be 1.0: {:.6}",
        h
    );
}

#[test]
fn test_support_point_achieves_max() {
    let pts = cube_points(2.0);
    let hull = ConvexHull3d::build(&pts);
    let dir = [1.0, 1.0, 1.0_f64];
    let sp = hull.support_point(normalize(dir));
    let h = hull.support(normalize(dir));
    let achieved = dot(sp, normalize(dir));
    assert!(
        approx_eq(achieved, h, 1e-9),
        "support point should achieve max: {:.6} vs {:.6}",
        achieved,
        h
    );
}

// ── SupportFunction ──

/// An empty vertex set has no support vertex: the origin is returned
/// instead of panicking.
#[test]
fn test_support_point_of_an_empty_vertex_set_is_the_origin() {
    let direction = [1.0, -2.0, 3.0];
    let function = SupportFunction::new(Vec::new());
    assert_eq!(function.support_point(direction), [0.0; 3]);
    assert_eq!(function.support(direction), f64::NEG_INFINITY);
    let hull = ConvexHull3d {
        vertices: Vec::new(),
        half_edges: Vec::new(),
        faces: Vec::new(),
        triangles: Vec::new(),
    };
    assert_eq!(hull.support_point(direction), [0.0; 3]);
    // A non-empty set still reports its extreme vertex.
    let cube = SupportFunction::new(cube_points(1.0));
    assert_eq!(cube.support_point([1.0, 1.0, 1.0]), [1.0, 1.0, 1.0]);
}

#[test]
fn test_support_fn_positive_along_positive_dir() {
    let sf = SupportFunction::new(cube_points(1.0));
    let h = sf.support([1.0, 0.0, 0.0]);
    assert!(
        approx_eq(h, 1.0, 1e-9),
        "support along +x should be 1.0: {:.6}",
        h
    );
}

#[test]
fn test_support_fn_symmetric() {
    let sf = SupportFunction::new(cube_points(1.0));
    let h_pos = sf.support([0.0, 1.0, 0.0]);
    let h_neg = sf.support([0.0, -1.0, 0.0]);
    assert!(
        approx_eq(h_pos, h_neg, 1e-9),
        "cube support ±y should be equal"
    );
}

// ── MinkowskiSum ──

#[test]
fn test_minkowski_sum_contains_centroid_of_a() {
    let a = SupportFunction::new(cube_points(1.0));
    let b = SupportFunction::new(octahedron_points());
    let ms = MinkowskiSum::new(a, b);
    // centroid of A is origin; it should be in A⊕B
    assert!(
        ms.contains_point([0.0; 3], 1e-9),
        "centroid should be inside Minkowski sum"
    );
}

#[test]
fn test_minkowski_sum_support_additive() {
    let a = SupportFunction::new(cube_points(1.0));
    let b = SupportFunction::new(cube_points(0.5));
    let ms = MinkowskiSum::new(a.clone(), b.clone());
    let dir = normalize([1.0, 1.0, 0.0]);
    let h_ms = ms.support(dir);
    let h_sum = a.support(dir) + b.support(dir);
    assert!(
        approx_eq(h_ms, h_sum, 1e-9),
        "Minkowski sum support should be additive"
    );
}

// ── GjkDistance ──

fn translated(points: &[Vec3], offset: Vec3) -> Vec<Vec3> {
    points.iter().map(|&v| add(v, offset)).collect()
}

/// Distance between the cube `[-1, 1]³` and the same cube moved by `offset`.
fn cube_pair_distance(offset: Vec3) -> f64 {
    let a = SupportFunction::new(cube_points(1.0));
    let b = SupportFunction::new(translated(&cube_points(1.0), offset));
    GjkDistance::default().distance(&a, &b)
}

#[test]
fn test_gjk_distance_cubes_face_to_face() {
    // Gap of 1 along x.
    let d = cube_pair_distance([3.0, 0.0, 0.0]);
    assert!(approx_eq(d, 1.0, 1e-9), "face-to-face distance = {d}");
    // Gap of 3 along y (the configuration that used to return 50.75).
    let d = cube_pair_distance([0.0, 5.0, 0.0]);
    assert!(approx_eq(d, 3.0, 1e-9), "face-to-face distance = {d}");
    let d = cube_pair_distance([0.0, 0.0, -4.5]);
    assert!(approx_eq(d, 2.5, 1e-9), "face-to-face distance = {d}");
}

#[test]
fn test_gjk_distance_cubes_edge_to_edge_and_vertex_to_vertex() {
    // Edge to edge: gaps of 2 along x and y.
    let d = cube_pair_distance([4.0, 4.0, 0.0]);
    assert!(
        approx_eq(d, 2.0 * 2.0_f64.sqrt(), 1e-9),
        "edge-to-edge distance = {d}"
    );
    // Vertex to vertex: gaps of 1 along every axis.
    let d = cube_pair_distance([3.0, 3.0, 3.0]);
    assert!(
        approx_eq(d, 3.0_f64.sqrt(), 1e-9),
        "vertex-to-vertex distance = {d}"
    );
    // Mixed offset with negative components and a face-region contact.
    let d = cube_pair_distance([-3.5, 0.5, -0.25]);
    assert!(approx_eq(d, 1.5, 1e-9), "face-to-face distance = {d}");
}

#[test]
fn test_gjk_distance_is_symmetric_and_translation_invariant() {
    let a_pts = octahedron_points();
    let b_pts = translated(&cube_points(0.5), [2.0, 1.5, -1.0]);
    let a = SupportFunction::new(a_pts.clone());
    let b = SupportFunction::new(b_pts.clone());
    let gjk = GjkDistance::default();
    let d_ab = gjk.distance(&a, &b);
    let d_ba = gjk.distance(&b, &a);
    assert!(d_ab > 0.0);
    assert!(approx_eq(d_ab, d_ba, 1e-12), "{d_ab} vs {d_ba}");
    // Moving both bodies by the same large vector must not change it.
    let shift = [1.0e3, -2.0e3, 5.0e2];
    let a_far = SupportFunction::new(translated(&a_pts, shift));
    let b_far = SupportFunction::new(translated(&b_pts, shift));
    let d_far = gjk.distance(&a_far, &b_far);
    assert!(approx_eq(d_ab, d_far, 1e-8), "{d_ab} vs {d_far}");
}

#[test]
fn test_gjk_distance_overlapping_returns_exactly_zero() {
    // Nested cubes: the origin is deep inside the Minkowski difference.
    let a = SupportFunction::new(cube_points(1.0));
    let b = SupportFunction::new(cube_points(0.5));
    let gjk = GjkDistance::default();
    assert_eq!(gjk.distance(&a, &b), 0.0);
    assert_eq!(gjk.distance(&b, &a), 0.0);
    // Partially overlapping cubes, off-axis.
    assert_eq!(cube_pair_distance([1.0, 0.5, -0.25]), 0.0);
    assert_eq!(cube_pair_distance([0.3, -1.2, 1.9]), 0.0);
}

#[test]
fn test_gjk_distance_touching_is_zero() {
    // Face contact, edge contact and vertex contact.
    assert!(cube_pair_distance([2.0, 0.0, 0.0]) <= 1e-9);
    assert!(cube_pair_distance([2.0, 2.0, 0.0]) <= 1e-9);
    assert!(cube_pair_distance([2.0, 2.0, 2.0]) <= 1e-9);
    assert!(cube_pair_distance([2.0, 0.7, -0.4]) <= 1e-9);
}

#[test]
fn test_gjk_distance_degenerate_shapes() {
    let gjk = GjkDistance::default();
    // A single point against a cube.
    let point = SupportFunction::new(vec![[0.0, 0.0, 0.0]]);
    let cube = SupportFunction::new(translated(&cube_points(1.0), [4.0, 0.0, 0.0]));
    assert!(approx_eq(gjk.distance(&point, &cube), 3.0, 1e-12));
    assert!(approx_eq(gjk.distance(&cube, &point), 3.0, 1e-12));
    // Two skew perpendicular segments at height 2.
    let seg_a = SupportFunction::new(vec![[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]]);
    let seg_b = SupportFunction::new(vec![[0.0, -1.0, 2.0], [0.0, 1.0, 2.0]]);
    assert!(approx_eq(gjk.distance(&seg_a, &seg_b), 2.0, 1e-12));
    // Two parallel unit squares (planar vertex sets) one apart.
    let square = |z: f64| {
        SupportFunction::new(vec![
            [-1.0, -1.0, z],
            [1.0, -1.0, z],
            [1.0, 1.0, z],
            [-1.0, 1.0, z],
        ])
    };
    assert!(approx_eq(
        gjk.distance(&square(0.0), &square(1.0)),
        1.0,
        1e-12
    ));
    // Coplanar squares side by side with a gap of 1.5.
    let shifted = SupportFunction::new(vec![
        [2.5, -1.0, 0.0],
        [4.5, -1.0, 0.0],
        [4.5, 1.0, 0.0],
        [2.5, 1.0, 0.0],
    ]);
    assert!(approx_eq(gjk.distance(&square(0.0), &shifted), 1.5, 1e-12));
    // Identical points.
    assert_eq!(gjk.distance(&point, &point), 0.0);
}

#[test]
fn test_gjk_distance_empty_and_non_finite_inputs() {
    let gjk = GjkDistance::default();
    let empty = SupportFunction::new(Vec::new());
    let cube = SupportFunction::new(cube_points(1.0));
    assert_eq!(gjk.distance(&empty, &cube), f64::INFINITY);
    assert_eq!(gjk.distance(&cube, &empty), f64::INFINITY);
    assert_eq!(gjk.distance(&empty, &empty), f64::INFINITY);
    let bad = SupportFunction::new(vec![[0.0, f64::NAN, 0.0], [1.0, 0.0, 0.0]]);
    assert!(gjk.distance(&bad, &cube).is_nan());
    let inf = SupportFunction::new(vec![[f64::INFINITY, 0.0, 0.0]]);
    assert!(gjk.distance(&cube, &inf).is_nan());
}

/// Exact distance between two disjoint convex polytopes from the closest
/// feature pairs of their hull triangulations (vertex-face and edge-edge).
fn feature_pair_distance(a: &[Vec3], b: &[Vec3]) -> f64 {
    use crate::convex_hull::ConvexHull3D;
    let hull_a = ConvexHull3D::build(a).expect("hull of A");
    let hull_b = ConvexHull3D::build(b).expect("hull of B");
    let mut best = f64::INFINITY;
    for (hull, other) in [(&hull_a, &hull_b), (&hull_b, &hull_a)] {
        for &p in &hull.vertices {
            for f in &other.faces {
                let (q0, q1, q2) = (
                    other.vertices[f[0]],
                    other.vertices[f[1]],
                    other.vertices[f[2]],
                );
                best = best.min(norm(sub(
                    p,
                    acd_mesh::closest_point_on_triangle(p, q0, q1, q2),
                )));
            }
        }
    }
    let edges = |h: &ConvexHull3D| -> Vec<(Vec3, Vec3)> {
        h.faces
            .iter()
            .flat_map(|f| (0..3).map(move |k| (f[k], f[(k + 1) % 3])))
            .map(|(i, j)| (h.vertices[i], h.vertices[j]))
            .collect()
    };
    for (p1, q1) in edges(&hull_a) {
        for (p2, q2) in edges(&hull_b) {
            best = best.min(segment_segment_distance(p1, q1, p2, q2));
        }
    }
    best
}

/// Distance between segments `p1q1` and `p2q2` (Ericson 5.1.9).
fn segment_segment_distance(p1: Vec3, q1: Vec3, p2: Vec3, q2: Vec3) -> f64 {
    let d1 = sub(q1, p1);
    let d2 = sub(q2, p2);
    let r = sub(p1, p2);
    let a = dot(d1, d1);
    let e = dot(d2, d2);
    let f = dot(d2, r);
    let (s, t);
    if a <= 1e-30 && e <= 1e-30 {
        s = 0.0;
        t = 0.0;
    } else if a <= 1e-30 {
        s = 0.0;
        t = (f / e).clamp(0.0, 1.0);
    } else {
        let c = dot(d1, r);
        if e <= 1e-30 {
            t = 0.0;
            s = (-c / a).clamp(0.0, 1.0);
        } else {
            let b = dot(d1, d2);
            let denom = a * e - b * b;
            let mut s_val = if denom > 1e-30 {
                ((b * f - c * e) / denom).clamp(0.0, 1.0)
            } else {
                0.0
            };
            let mut t_val = (b * s_val + f) / e;
            if t_val < 0.0 {
                t_val = 0.0;
                s_val = (-c / a).clamp(0.0, 1.0);
            } else if t_val > 1.0 {
                t_val = 1.0;
                s_val = ((b - c) / a).clamp(0.0, 1.0);
            }
            s = s_val;
            t = t_val;
        }
    }
    let c1 = add(p1, scale(d1, s));
    let c2 = add(p2, scale(d2, t));
    norm(sub(c1, c2))
}

#[test]
fn test_gjk_distance_matches_feature_pair_oracle_on_random_polytopes() {
    // Deterministic xorshift generator.
    let mut state = 0x9E37_79B9_7F4A_7C15_u64;
    let mut next = move || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state >> 11) as f64 / (1_u64 << 53) as f64
    };
    let gjk = GjkDistance::default();
    for case in 0..60 {
        let n_a = 5 + (next() * 8.0) as usize;
        let n_b = 5 + (next() * 8.0) as usize;
        let a: Vec<Vec3> = (0..n_a)
            .map(|_| [next() * 2.0 - 1.0, next() * 2.0 - 1.0, next() * 2.0 - 1.0])
            .collect();
        // Direction on the sphere and a separation of 3.7..6.5, which
        // keeps the two polytopes (radius <= sqrt(3)) disjoint.
        let dir = normalize([next() - 0.5, next() - 0.5, next() - 0.5]);
        let sep = 3.7 + 2.8 * next();
        let off = scale(dir, sep);
        let b: Vec<Vec3> = (0..n_b)
            .map(|_| {
                add(
                    [next() * 2.0 - 1.0, next() * 2.0 - 1.0, next() * 2.0 - 1.0],
                    off,
                )
            })
            .collect();
        let expected = feature_pair_distance(&a, &b);
        let got = gjk.distance(&SupportFunction::new(a), &SupportFunction::new(b));
        assert!(
            (got - expected).abs() <= 1e-9 * expected.max(1.0),
            "case {case}: gjk {got} oracle {expected}"
        );
    }
}

#[test]
fn test_gjk_distance_many_vertex_spheres() {
    // Fibonacci-lattice point sets approximating unit spheres.
    let sphere = |n: usize, centre: Vec3| -> Vec<Vec3> {
        let golden = std::f64::consts::PI * (3.0 - 5.0_f64.sqrt());
        (0..n)
            .map(|i| {
                let y = 1.0 - 2.0 * (i as f64 + 0.5) / n as f64;
                let r = (1.0 - y * y).sqrt();
                let phi = golden * i as f64;
                add(centre, [r * phi.cos(), y, r * phi.sin()])
            })
            .collect()
    };
    let a = sphere(400, [0.0, 0.0, 0.0]);
    let b = sphere(400, [0.0, 4.0, 0.0]);
    let expected = feature_pair_distance(&a, &b);
    let got = GjkDistance::default().distance(&SupportFunction::new(a), &SupportFunction::new(b));
    assert!(
        (got - expected).abs() <= 1e-9,
        "gjk {got} oracle {expected}"
    );
    // The faceted spheres are slightly smaller than the true ones, so the
    // distance is a little more than the analytic gap of 2.
    assert!(got > 2.0 && got < 2.05, "sphere gap = {got}");
}

// ── ConvexContainment ──

#[test]
fn test_centroid_always_inside_hull() {
    let pts = octahedron_points();
    let hull = ConvexHull3d::build(&pts);
    let c = hull.centroid();
    let cc = ConvexContainment::new(hull);
    assert!(cc.contains_point(c, 1e-9), "centroid should be inside hull");
}

#[test]
fn test_far_point_outside_hull() {
    let pts = octahedron_points();
    let hull = ConvexHull3d::build(&pts);
    let cc = ConvexContainment::new(hull);
    assert!(
        !cc.contains_point([10.0, 0.0, 0.0], 1e-9),
        "far point should be outside hull"
    );
}

#[test]
fn test_aabb_inside_hull() {
    let pts = cube_points(2.0);
    let hull = ConvexHull3d::build(&pts);
    let cc = ConvexContainment::new(hull);
    // A small AABB inside the cube
    let inside = cc.contains_aabb([-0.5, -0.5, -0.5], [0.5, 0.5, 0.5]);
    assert!(inside, "small AABB should be inside larger cube hull");
}

// ── ConvexIntersection ──

#[test]
fn test_sat_overlapping_hulls_intersect() {
    let pts_a = cube_points(1.0);
    let pts_b = cube_points(0.5); // inside a
    let hull_a = ConvexHull3d::build(&pts_a);
    let hull_b = ConvexHull3d::build(&pts_b);
    assert!(
        ConvexIntersection::intersects(&hull_a, &hull_b),
        "nested hulls should intersect"
    );
}

#[test]
fn test_sat_separated_hulls_no_intersect() {
    let pts_a = cube_points(0.4);
    let pts_b: Vec<Vec3> = cube_points(0.4)
        .iter()
        .map(|&v| add(v, [5.0, 0.0, 0.0]))
        .collect();
    let hull_a = ConvexHull3d::build(&pts_a);
    let hull_b = ConvexHull3d::build(&pts_b);
    assert!(
        !ConvexIntersection::intersects(&hull_a, &hull_b),
        "separated hulls should not intersect"
    );
}

// ── InertiaFromHull ──

#[test]
fn test_inertia_diagonal_positive() {
    let pts = cube_points(1.0);
    let hull = ConvexHull3d::build(&pts);
    let inertia = InertiaFromHull::compute(&hull, 1000.0);
    for (i, row) in inertia.iter().enumerate() {
        assert!(
            row[i] > 0.0,
            "diagonal inertia[{}][{}] should be positive: {:.6}",
            i,
            i,
            row[i]
        );
    }
}

#[test]
fn test_inertia_symmetric() {
    let pts = octahedron_points();
    let hull = ConvexHull3d::build(&pts);
    let inertia = InertiaFromHull::compute(&hull, 1000.0);
    for (i, row_i) in inertia.iter().enumerate() {
        for (j, &val_ij) in row_i.iter().enumerate() {
            assert!(
                approx_eq(val_ij, inertia[j][i], 1e-9),
                "inertia tensor should be symmetric at ({},{}): {:.6} vs {:.6}",
                i,
                j,
                val_ij,
                inertia[j][i]
            );
        }
    }
}

// ── HullShrink ──

#[test]
fn test_shrink_reduces_support() {
    let pts = cube_points(1.0);
    let hull = ConvexHull3d::build(&pts);
    let shrunk = HullShrink::shrink(&hull, 0.1);
    let sf_orig = SupportFunction::new(pts.clone());
    let sf_shrunk = SupportFunction::new(shrunk);
    let dir = normalize([1.0, 0.0, 0.0]);
    assert!(
        sf_shrunk.support(dir) < sf_orig.support(dir),
        "shrunk hull should have smaller support"
    );
}
