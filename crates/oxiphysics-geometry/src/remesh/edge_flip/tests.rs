//! Regression tests for the edge-flip engine.
//!
//! The previous flip loops iterated a prebuilt edge map while rewriting the
//! triangles, so later edges used stale opposite vertices; they wrote the flipped
//! triangles with the wrong winding, never checked that the new edge was free or
//! that the quad was convex, and evaluated the in-circle test without regard to
//! the orientation of the triangle.
//!
//! Every test asserts that a flip actually happened (or, for the negative cases,
//! pairs the refusal with a nearly identical configuration that does flip), so a
//! policy that never fires cannot pass.

use super::*;
use crate::remesh::{isotropic_remesh, mesh_aspect_ratio_avg, mesh_quality_min_angle};
use crate::test_meshes::{flat_grid, icosphere};
use crate::triangle_mesh::TriangleMesh;
use oxiphysics_core::math::Vec3;
use std::collections::{BTreeMap, BTreeSet};

type Verts = Vec<[f64; 3]>;
type Tris = Vec<[usize; 3]>;

/// Deterministic uniform samples in `[0, 1)` from a 64-bit LCG.
struct Lcg(u64);

impl Lcg {
    fn next_unit(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }

    fn next_signed(&mut self) -> f64 {
        self.next_unit() * 2.0 - 1.0
    }
}

fn no_lock(_: usize, _: usize) -> bool {
    false
}

/// Policy that accepts every candidate: exercises the topology bookkeeping of
/// the engine independently of any geometric criterion.
struct AlwaysFlip;

impl FlipPolicy for AlwaysFlip {
    fn accepts(&self, _flip: &FlipCandidate) -> bool {
        true
    }
}

// ── Mesh invariants ──────────────────────────────────────────────────────────

fn directed_edges(tris: &[[usize; 3]]) -> Vec<(usize, usize)> {
    tris.iter()
        .flat_map(|t| (0..3).map(move |k| (t[k], t[(k + 1) % 3])))
        .collect()
}

fn undirected_edges(tris: &[[usize; 3]]) -> BTreeSet<(usize, usize)> {
    directed_edges(tris)
        .into_iter()
        .map(|(a, b)| (a.min(b), a.max(b)))
        .collect()
}

/// No directed edge is used twice (so every edge has at most two triangles, with
/// opposite directions) and no triangle is duplicated or degenerate.
fn assert_manifold_and_consistent(tris: &[[usize; 3]]) {
    let mut seen: BTreeSet<(usize, usize)> = BTreeSet::new();
    for edge in directed_edges(tris) {
        assert!(seen.insert(edge), "directed edge {edge:?} is used twice");
    }
    let mut faces: BTreeSet<[usize; 3]> = BTreeSet::new();
    for t in tris {
        assert!(
            t[0] != t[1] && t[1] != t[2] && t[0] != t[2],
            "degenerate triangle {t:?}"
        );
        let mut sorted = *t;
        sorted.sort_unstable();
        assert!(faces.insert(sorted), "triangle {sorted:?} appears twice");
    }
}

/// Every directed edge has its opposite (a closed surface).
fn assert_closed(tris: &[[usize; 3]]) {
    let edges: BTreeSet<(usize, usize)> = directed_edges(tris).into_iter().collect();
    for &(a, b) in &edges {
        assert!(edges.contains(&(b, a)), "edge ({a}, {b}) has no opposite");
    }
}

/// Undirected edges used by exactly one triangle.
fn boundary_edges(tris: &[[usize; 3]]) -> BTreeSet<(usize, usize)> {
    let mut count: BTreeMap<(usize, usize), usize> = BTreeMap::new();
    for (a, b) in directed_edges(tris) {
        *count.entry((a.min(b), a.max(b))).or_default() += 1;
    }
    count
        .into_iter()
        .filter(|&(_, uses)| uses == 1)
        .map(|(edge, _)| edge)
        .collect()
}

/// Every triangle of a closed surface around the origin faces away from it.
fn assert_outward(verts: &[[f64; 3]], tris: &[[usize; 3]]) {
    for t in tris {
        let [a, b, c] = [verts[t[0]], verts[t[1]], verts[t[2]]];
        let normal = cross3(sub3(b, a), sub3(c, a));
        let centroid = [
            (a[0] + b[0] + c[0]) / 3.0,
            (a[1] + b[1] + c[1]) / 3.0,
            (a[2] + b[2] + c[2]) / 3.0,
        ];
        assert!(
            dot3(normal, centroid) > 0.0,
            "triangle {t:?} faces the centre of the sphere"
        );
    }
}

/// Sum of the unsigned triangle areas.
fn total_area(verts: &[[f64; 3]], tris: &[[usize; 3]]) -> f64 {
    tris.iter()
        .map(|t| {
            let (a, b, c) = (verts[t[0]], verts[t[1]], verts[t[2]]);
            0.5 * len3(cross3(sub3(b, a), sub3(c, a)))
        })
        .sum()
}

/// z component of the triangle normal (positive: counter-clockwise in x-y).
fn normal_z(verts: &[[f64; 3]], t: [usize; 3]) -> f64 {
    let (a, b, c) = (verts[t[0]], verts[t[1]], verts[t[2]]);
    cross3(sub3(b, a), sub3(c, a))[2]
}

/// Valence (number of incident edges) of every vertex.
fn valences(vertex_count: usize, tris: &[[usize; 3]]) -> Vec<usize> {
    let mut valence = vec![0; vertex_count];
    for (a, b) in undirected_edges(tris) {
        valence[a] += 1;
        valence[b] += 1;
    }
    valence
}

fn valence_deviation_sum(vertex_count: usize, tris: &[[usize; 3]]) -> i64 {
    valences(vertex_count, tris)
        .iter()
        .map(|&v| (v as i64 - 6).pow(2))
        .sum()
}

/// Flips the edge `{a, b}` of a consistently oriented mesh without any checks.
/// Independent of the engine; used to build meshes with a known bad edge.
fn flip_manually(tris: &mut [[usize; 3]], a: usize, b: usize) {
    let ti = tris
        .iter()
        .position(|t| (0..3).any(|k| t[k] == a && t[(k + 1) % 3] == b))
        .expect("a triangle holds a -> b");
    let tj = tris
        .iter()
        .position(|t| (0..3).any(|k| t[k] == b && t[(k + 1) % 3] == a))
        .expect("a triangle holds b -> a");
    let opposite = |t: [usize; 3]| {
        *t.iter()
            .find(|&&v| v != a && v != b)
            .expect("a triangle has a third vertex")
    };
    let (c, d) = (opposite(tris[ti]), opposite(tris[tj]));
    tris[ti] = [c, a, d];
    tris[tj] = [d, b, c];
}

// ── Two-triangle quads ───────────────────────────────────────────────────────

/// Thin rhombus with a long diagonal `a-b` (indices 0-1) and a short one `c-d`
/// (indices 2-3) in the plane `y = 0`, i.e. in the X-Z plane used by the Delaunay
/// flips, where `z` plays the role of the second coordinate.
fn rhombus_xz() -> Verts {
    vec![
        [0.0, 0.0, 0.0],  // a
        [4.0, 0.0, 0.0],  // b
        [2.0, 0.0, 1.0],  // c
        [2.0, 0.0, -1.0], // d
    ]
}

/// A thin rhombus in the plane `z = 0` (counter-clockwise seen from +z).
fn rhombus_xy() -> Verts {
    vec![
        [0.0, 0.0, 0.0],
        [4.0, 0.0, 0.0],
        [2.0, 0.5, 0.0],
        [2.0, -0.5, 0.0],
    ]
}

/// Triangles `(a, b, c)` and `(b, a, d)` with the long diagonal `a-b`.
const LONG_DIAGONAL: [[usize; 3]; 2] = [[0, 1, 2], [1, 0, 3]];
/// The same quad after the flip: `(c, a, d)` and `(d, b, c)`.
const SHORT_DIAGONAL: [[usize; 3]; 2] = [[2, 0, 3], [3, 1, 2]];

fn reversed(tris: &[[usize; 3]]) -> Tris {
    tris.iter().map(|&[a, b, c]| [a, c, b]).collect()
}

#[test]
fn quality_flip_writes_the_new_triangles_with_the_old_winding() {
    let verts = rhombus_xy();
    let mut tris = LONG_DIAGONAL.to_vec();
    let flips = flip_edges_for_quality(&mut tris, &verts);
    assert_eq!(flips, 1);
    assert_eq!(tris, SHORT_DIAGONAL.to_vec());
    // Counter-clockwise (+z) before, counter-clockwise after.
    for &t in &LONG_DIAGONAL {
        assert!(normal_z(&verts, t) > 0.0);
    }
    for &t in &tris {
        assert!(normal_z(&verts, t) > 0.0, "{t:?} was inverted");
    }
    // The flip improved the worst angle (14 degrees -> 28 degrees).
    assert!(mesh_quality_min_angle(&verts, &tris) > mesh_quality_min_angle(&verts, &LONG_DIAGONAL));
    // A second call has nothing left to do.
    assert_eq!(flip_edges_for_quality(&mut tris, &verts), 0);
}

#[test]
fn quality_flip_keeps_clockwise_triangles_clockwise() {
    let verts = rhombus_xy();
    let mut tris = reversed(&LONG_DIAGONAL);
    assert_eq!(flip_edges_for_quality(&mut tris, &verts), 1);
    // The engine names the first triangle's directed edge b -> a, so the quad is
    // rewritten as (c, b, d) and (d, a, c).
    assert_eq!(tris, vec![[2, 1, 3], [3, 0, 2]]);
    for &t in &tris {
        assert!(normal_z(&verts, t) < 0.0, "{t:?} was inverted");
    }
    assert_eq!(flip_edges_for_quality(&mut tris, &verts), 0);
}

#[test]
fn quality_flip_leaves_the_good_diagonal_alone() {
    let verts = rhombus_xy();
    for mut tris in [SHORT_DIAGONAL.to_vec(), reversed(&SHORT_DIAGONAL)] {
        let before = tris.clone();
        assert_eq!(flip_edges_for_quality(&mut tris, &verts), 0);
        assert_eq!(tris, before);
    }
}

#[test]
fn delaunay_flip_fixes_the_bad_diagonal_for_both_windings() {
    let verts = rhombus_xz();
    let mut ccw = LONG_DIAGONAL.to_vec();
    assert_eq!(delaunay_edge_flip(&mut ccw, &verts), 1);
    assert_eq!(ccw, SHORT_DIAGONAL.to_vec());

    // The old in-circle test ignored the orientation, so for the clockwise
    // variant it neither flipped this edge nor kept the Delaunay one.
    let mut cw = reversed(&LONG_DIAGONAL);
    assert_eq!(delaunay_edge_flip(&mut cw, &verts), 1);
    assert_eq!(cw, vec![[2, 1, 3], [3, 0, 2]]);
}

#[test]
fn delaunay_flip_keeps_an_already_delaunay_pair_for_both_windings() {
    let verts = rhombus_xz();
    for mut tris in [SHORT_DIAGONAL.to_vec(), reversed(&SHORT_DIAGONAL)] {
        let before = tris.clone();
        assert_eq!(delaunay_edge_flip(&mut tris, &verts), 0);
        assert_eq!(tris, before);
    }
}

#[test]
fn delaunay_flip_ignores_the_y_coordinate() {
    let mut verts = rhombus_xz();
    for (index, v) in verts.iter_mut().enumerate() {
        v[1] = 10.0 * index as f64 - 3.0;
    }
    let mut tris = LONG_DIAGONAL.to_vec();
    assert_eq!(delaunay_edge_flip(&mut tris, &verts), 1);
    assert_eq!(tris, SHORT_DIAGONAL.to_vec());
}

#[test]
fn delaunay_flip_refuses_folded_triangle_pairs() {
    // d lies inside the circumcircle of (a, b, c) but on the same side of `a-b`
    // as c: the triangles overlap in the projection and the quad has no
    // flipped triangulation.
    let verts: Verts = vec![
        [0.0, 0.0, 0.0],
        [4.0, 0.0, 0.0],
        [2.0, 0.0, 1.0],
        [2.0, 0.0, 0.5],
    ];
    assert!(point_in_circumcircle(
        verts[0], verts[1], verts[2], verts[3]
    ));
    let mut tris = LONG_DIAGONAL.to_vec();
    assert_eq!(delaunay_edge_flip(&mut tris, &verts), 0);
    assert_eq!(tris, LONG_DIAGONAL.to_vec());
}

#[test]
fn delaunay_flip_ignores_meshes_that_are_degenerate_in_the_projection() {
    // All points in the plane z = 0 project onto a line in the X-Z plane.
    let verts = rhombus_xy();
    let mut tris = LONG_DIAGONAL.to_vec();
    assert_eq!(delaunay_edge_flip(&mut tris, &verts), 0);
    assert_eq!(tris, LONG_DIAGONAL.to_vec());
}

// ── In-circle predicate ──────────────────────────────────────────────────────

#[test]
fn circumcircle_test_does_not_depend_on_the_winding() {
    let a = [0.0, 0.0, 0.0];
    let b = [4.0, 0.0, 0.0];
    let c = [2.0, 0.0, 1.0];
    let inside = [2.0, 0.0, -1.0];
    let outside = [2.0, 0.0, -4.5];
    for (p, q, r) in [
        (a, b, c),
        (b, c, a),
        (c, a, b),
        (a, c, b),
        (c, b, a),
        (b, a, c),
    ] {
        assert!(point_in_circumcircle(p, q, r, inside), "{p:?} {q:?} {r:?}");
        assert!(
            !point_in_circumcircle(p, q, r, outside),
            "{p:?} {q:?} {r:?}"
        );
        // The vertices themselves are on the circle.
        assert!(!point_in_circumcircle(p, q, r, a));
        assert!(!point_in_circumcircle(p, q, r, b));
        assert!(!point_in_circumcircle(p, q, r, c));
    }
}

#[test]
fn circumcircle_test_treats_cocircular_points_as_outside() {
    // The corners of a unit square are cocircular.
    let a = [0.0, 0.0, 0.0];
    let b = [1.0, 0.0, 0.0];
    let c = [1.0, 0.0, 1.0];
    let d = [0.0, 0.0, 1.0];
    assert!(!point_in_circumcircle(a, b, c, d));
    assert!(!point_in_circumcircle(a, c, b, d));
    // Same with a large offset, where rounding errors matter.
    let shift = |p: [f64; 3]| [p[0] + 1.0e6, p[1], p[2] - 3.0e5];
    assert!(!point_in_circumcircle(
        shift(a),
        shift(b),
        shift(c),
        shift(d)
    ));
}

#[test]
fn circumcircle_test_rejects_degenerate_triangles() {
    let a = [0.0, 0.0, 0.0];
    let b = [1.0, 0.0, 0.0];
    let c = [2.0, 0.0, 0.0];
    assert!(!point_in_circumcircle(a, b, c, [1.0, 0.0, 0.1]));
    assert!(!point_in_circumcircle(a, a, a, [0.0, 0.0, 0.0]));
    assert!(!point_in_circumcircle(
        [f64::NAN, 0.0, 0.0],
        b,
        [0.0, 0.0, 1.0],
        [0.2, 0.0, 0.2]
    ));
}

// ── Refusals ─────────────────────────────────────────────────────────────────

/// Closed tetrahedron whose faces are `(a, b, c)`, `(b, a, d)`, `(c, b, d)` and
/// `(a, c, d)` over a thin rhombus: flipping `a-b` would create `c-d`, which is
/// already an edge.
fn thin_tetrahedron() -> (Verts, Tris) {
    let verts: Verts = vec![
        [0.0, 0.0, 0.0],
        [4.0, 0.0, 0.0],
        [2.0, 0.5, 0.3],
        [2.0, -0.5, 0.3],
    ];
    (verts, vec![[0, 1, 2], [1, 0, 3], [2, 1, 3], [0, 2, 3]])
}

#[test]
fn a_flip_that_would_duplicate_an_existing_edge_is_refused() {
    let (verts, tris) = thin_tetrahedron();
    assert_closed(&tris);
    // The same two triangles without the closing faces do flip, so the refusal
    // below is caused by the existing edge alone.
    let mut open = tris[..2].to_vec();
    assert_eq!(flip_edges_for_quality(&mut open, &verts), 1);

    let mut quality = tris.clone();
    assert_eq!(flip_edges_for_quality(&mut quality, &verts), 0);
    assert_eq!(quality, tris);

    let mut valence = tris.clone();
    assert_eq!(flip_for_valence(&mut valence, &verts, no_lock), 0);
    assert_eq!(valence, tris);
}

#[test]
fn delaunay_flip_does_not_duplicate_an_existing_edge() {
    // In the X-Z plane the rhombus is flippable; closing it with faces that
    // already contain `c-d` must block the flip.
    let verts = rhombus_xz();
    let mut open = LONG_DIAGONAL.to_vec();
    assert_eq!(delaunay_edge_flip(&mut open, &verts), 1);

    let mut closed: Tris = vec![[0, 1, 2], [1, 0, 3], [2, 1, 3], [0, 2, 3]];
    assert_closed(&closed);
    let before = closed.clone();
    assert_eq!(delaunay_edge_flip(&mut closed, &verts), 0);
    assert_eq!(closed, before);
}

#[test]
fn edges_with_inconsistent_winding_are_not_flipped() {
    let verts = rhombus_xy();
    // Both triangles traverse a -> b.
    let mut tris: Tris = vec![[0, 1, 2], [0, 1, 3]];
    let before = tris.clone();
    assert_eq!(flip_edges_for_quality(&mut tris, &verts), 0);
    assert_eq!(flip_for_valence(&mut tris, &verts, no_lock), 0);
    assert_eq!(delaunay_edge_flip(&mut tris, &rhombus_xz()), 0);
    assert_eq!(tris, before);
}

#[test]
fn non_manifold_edges_are_not_flipped() {
    let mut verts = rhombus_xy();
    verts.push([2.0, 0.0, 1.0]);
    // Three triangles share the edge a-b.
    let mut tris: Tris = vec![[0, 1, 2], [1, 0, 3], [0, 1, 4]];
    let before = tris.clone();
    assert_eq!(flip_edges_for_quality(&mut tris, &verts), 0);
    assert_eq!(tris, before);
}

#[test]
fn triangle_pairs_with_the_same_third_vertex_are_not_flipped() {
    let verts = rhombus_xy();
    // (a, b, c) and (b, a, c): a two-triangle pillow.
    let mut tris: Tris = vec![[0, 1, 2], [1, 0, 2]];
    let before = tris.clone();
    assert_eq!(flip_edges_for_quality(&mut tris, &verts), 0);
    assert_eq!(tris, before);
}

#[test]
fn non_convex_quads_are_not_flipped() {
    // Vertex a is reflex: the segment c-d misses a-b, and the flipped triangle
    // (c, a, d) would be inverted. Its worst angle (12.5 degrees) is still larger
    // than the current one (12.1 degrees), so only the convexity check refuses.
    let verts: Verts = vec![
        [0.0, 0.0, 0.0],
        [4.0, 0.0, 0.0],
        [-0.2, 0.9, 0.0],
        [-0.2, -0.9, 0.0],
    ];
    let mut tris = LONG_DIAGONAL.to_vec();
    assert_eq!(flip_edges_for_quality(&mut tris, &verts), 0);
    assert_eq!(tris, LONG_DIAGONAL.to_vec());
    let before = mesh_quality_min_angle(&verts, &tris);
    let flipped = mesh_quality_min_angle(&verts, &SHORT_DIAGONAL);
    assert!(
        flipped > before,
        "the flip must look attractive: {flipped} vs {before}"
    );
}

#[test]
fn invalid_indices_and_degenerate_triangles_are_ignored() {
    let verts = rhombus_xy();
    let mut tris: Tris = vec![[0, 1, 2], [1, 0, 3], [0, 1, 99], [2, 2, 3]];
    let flips = flip_edges_for_quality(&mut tris, &verts);
    assert_eq!(flips, 1);
    assert_eq!(&tris[..2], &SHORT_DIAGONAL[..]);
    assert_eq!(tris[2], [0, 1, 99]);
    assert_eq!(tris[3], [2, 2, 3]);
}

#[test]
fn empty_inputs_do_nothing() {
    let verts = rhombus_xy();
    let mut none: Tris = Vec::new();
    assert_eq!(flip_edges_for_quality(&mut none, &verts), 0);
    assert_eq!(delaunay_edge_flip(&mut none, &verts), 0);
    assert_eq!(flip_for_valence(&mut none, &verts, no_lock), 0);
    let mut tris = LONG_DIAGONAL.to_vec();
    assert_eq!(flip_edges_for_quality(&mut tris, &[]), 0);
    assert_eq!(tris, LONG_DIAGONAL.to_vec());
}

// ── Engine bookkeeping ───────────────────────────────────────────────────────

#[test]
fn passes_defer_edges_that_touch_a_flipped_triangle() {
    // With a policy that accepts everything, any stale bookkeeping would show up
    // as duplicated faces, reused directed edges or changed boundaries.
    let (verts, initial) = flat_grid(4, 4, 1.0);
    let boundary = boundary_edges(&initial);
    let mut tris = initial.clone();
    let mut first_pass = 0;
    for pass in 0..6 {
        let flips = flip_pass(&mut tris, verts.len(), &mut AlwaysFlip);
        if pass == 0 {
            first_pass = flips;
        }
        // Each flip rewrites two distinct triangles, and a rewritten triangle is
        // not touched again in the same pass.
        assert!(2 * flips <= tris.len(), "pass {pass}: {flips} flips");
        assert_manifold_and_consistent(&tris);
        assert_eq!(tris.len(), initial.len());
        assert_eq!(boundary_edges(&tris), boundary, "pass {pass}");
        assert_eq!(
            valences(verts.len(), &tris).iter().sum::<usize>(),
            2 * undirected_edges(&tris).len()
        );
    }
    assert!(first_pass > 4, "only {first_pass} flips in the first pass");
}

#[test]
fn engine_stops_when_the_policy_accepts_nothing() {
    struct Never;
    impl FlipPolicy for Never {
        fn accepts(&self, _flip: &FlipCandidate) -> bool {
            false
        }
    }
    let (verts, initial) = flat_grid(3, 3, 1.0);
    let mut tris = initial.clone();
    assert_eq!(flip_edges(&mut tris, verts.len(), &mut Never), 0);
    assert_eq!(tris, initial);
}

#[test]
fn valence_criterion_targets_four_on_the_boundary_and_six_inside() {
    let verts = rhombus_xy();
    let flip = FlipCandidate {
        a: 0,
        b: 1,
        c: 2,
        d: 3,
        tri_ab: 0,
        tri_ba: 1,
    };
    // Valences 5, 7, 5, 5 become 4, 6, 6, 6 after the flip. Vertex a is on the
    // boundary (ideal valence 4) in the first policy and in the interior
    // (ideal valence 6) in the second.
    let policy = |a_on_boundary: bool| ValencePolicy {
        verts: &verts,
        locked: no_lock,
        degree: vec![5, 7, 5, 5],
        on_boundary: vec![a_on_boundary, false, false, false],
    };
    assert!(policy(true).accepts(&flip), "deviation 4 -> 0 must flip");
    assert!(
        !policy(false).accepts(&flip),
        "deviation 4 -> 4 must not flip"
    );
    // A vertex never drops below three edges.
    let thin = ValencePolicy {
        verts: &verts,
        locked: no_lock,
        degree: vec![3, 9, 4, 4],
        on_boundary: vec![false; 4],
    };
    assert!(!thin.accepts(&flip));
}

// ── Delaunay on a point cloud ────────────────────────────────────────────────

/// Jittered `n x n` cell grid in the X-Z plane, counter-clockwise in `(x, z)`,
/// with a random diagonal in every cell.
fn jittered_grid_xz(n: usize, jitter: f64, point_seed: u64, diagonal_seed: u64) -> (Verts, Tris) {
    let mut points = Lcg(point_seed);
    let mut verts = Vec::new();
    for j in 0..=n {
        for i in 0..=n {
            let dx = if i == 0 || i == n {
                0.0
            } else {
                jitter * points.next_signed()
            };
            let dz = if j == 0 || j == n {
                0.0
            } else {
                jitter * points.next_signed()
            };
            verts.push([i as f64 + dx, 0.0, j as f64 + dz]);
        }
    }
    let at = |i: usize, j: usize| j * (n + 1) + i;
    let mut diagonals = Lcg(diagonal_seed);
    let mut tris = Vec::new();
    for j in 0..n {
        for i in 0..n {
            let (p00, p10, p11, p01) = (at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1));
            if diagonals.next_unit() < 0.5 {
                tris.push([p00, p10, p11]);
                tris.push([p00, p11, p01]);
            } else {
                tris.push([p00, p10, p01]);
                tris.push([p10, p11, p01]);
            }
        }
    }
    (verts, tris)
}

/// Circumcentre and squared radius of a triangle in the X-Z projection.
fn circumcircle_xz(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> ([f64; 2], f64) {
    let (ax, az) = (a[0], a[2]);
    let (bx, bz) = (b[0] - ax, b[2] - az);
    let (cx, cz) = (c[0] - ax, c[2] - az);
    let d = 2.0 * (bx * cz - bz * cx);
    let ux = (cz * (bx * bx + bz * bz) - bz * (cx * cx + cz * cz)) / d;
    let uz = (bx * (cx * cx + cz * cz) - cx * (bx * bx + bz * bz)) / d;
    ([ax + ux, az + uz], ux * ux + uz * uz)
}

/// Every interior edge is locally Delaunay according to explicit circumcircles.
fn assert_locally_delaunay(verts: &[[f64; 3]], tris: &[[usize; 3]]) {
    let mut by_edge: BTreeMap<(usize, usize), Vec<(usize, usize)>> = BTreeMap::new();
    for (ti, t) in tris.iter().enumerate() {
        for k in 0..3 {
            let (a, b, c) = (t[k], t[(k + 1) % 3], t[(k + 2) % 3]);
            by_edge
                .entry((a.min(b), a.max(b)))
                .or_default()
                .push((ti, c));
        }
    }
    let mut interior = 0;
    for (&(a, b), uses) in &by_edge {
        let [(first, _), (_, d)] = uses[..] else {
            continue;
        };
        interior += 1;
        let t = tris[first];
        let (centre, radius_sq) = circumcircle_xz(verts[t[0]], verts[t[1]], verts[t[2]]);
        let dist_sq = (verts[d][0] - centre[0]).powi(2) + (verts[d][2] - centre[1]).powi(2);
        assert!(
            dist_sq >= radius_sq * (1.0 - 1e-9),
            "edge ({a}, {b}) is not Delaunay: vertex {d} is inside the circle of {t:?}"
        );
    }
    assert!(interior > 0);
}

fn assert_counter_clockwise_xz(verts: &[[f64; 3]], tris: &[[usize; 3]]) {
    for &t in tris {
        assert!(
            orient_xz(verts[t[0]], verts[t[1]], verts[t[2]]) > 0.0,
            "{t:?} lost its orientation"
        );
    }
}

#[test]
fn delaunay_flip_on_a_jittered_grid_reaches_the_delaunay_triangulation() {
    let (verts, initial) = jittered_grid_xz(8, 0.3, 7, 1);
    assert_counter_clockwise_xz(&verts, &initial);
    let mut tris = initial.clone();
    let boundary = boundary_edges(&tris);
    let area = total_area(&verts, &tris);

    let flips = delaunay_edge_flip(&mut tris, &verts);
    assert!(flips > 10, "only {flips} flips on a random triangulation");

    assert_manifold_and_consistent(&tris);
    assert_eq!(tris.len(), initial.len());
    assert_eq!(boundary_edges(&tris), boundary);
    assert!((total_area(&verts, &tris) - area).abs() < 1e-9);
    assert_counter_clockwise_xz(&verts, &tris);
    assert_locally_delaunay(&verts, &tris);
    // Nothing is left to flip.
    assert_eq!(delaunay_edge_flip(&mut tris, &verts), 0);
}

#[test]
fn delaunay_triangulation_of_points_in_general_position_is_unique() {
    // Two different initial triangulations of the same points end in the same
    // set of edges.
    let (verts, mut first) = jittered_grid_xz(7, 0.3, 21, 1);
    let (_, mut second) = jittered_grid_xz(7, 0.3, 21, 2);
    assert_ne!(undirected_edges(&first), undirected_edges(&second));
    assert!(delaunay_edge_flip(&mut first, &verts) > 0);
    assert!(delaunay_edge_flip(&mut second, &verts) > 0);
    assert_eq!(undirected_edges(&first), undirected_edges(&second));
}

#[test]
fn cocircular_lattice_points_are_left_alone() {
    // flat_grid lies in z = 0; swap y and z to put it in the X-Z plane.
    let (flat, tris) = flat_grid(6, 6, 1.0);
    let verts: Verts = flat.iter().map(|v| [v[0], 0.0, v[1]]).collect();
    assert_counter_clockwise_xz(&verts, &tris);
    let mut work = tris.clone();
    assert_eq!(delaunay_edge_flip(&mut work, &verts), 0);
    assert_eq!(work, tris);

    // The same grid in its own plane (z = 0) is degenerate in the X-Z
    // projection and must not be touched either.
    let mut flat_work = tris.clone();
    assert_eq!(delaunay_edge_flip(&mut flat_work, &flat), 0);
    assert_eq!(flat_work, tris);
}

#[test]
fn moving_one_lattice_point_creates_a_delaunay_violation_that_is_repaired() {
    let (flat, tris) = flat_grid(6, 6, 1.0);
    let mut verts: Verts = flat.iter().map(|v| [v[0], 0.0, v[1]]).collect();
    // Push the vertex (3, 3) towards (4, 4).
    let moved = 3 * 7 + 3;
    verts[moved] = [3.45, 0.0, 3.45];
    let mut work = tris.clone();
    let flips = delaunay_edge_flip(&mut work, &verts);
    assert!(flips > 0);
    assert_manifold_and_consistent(&work);
    assert_locally_delaunay(&verts, &work);
    assert_counter_clockwise_xz(&verts, &work);
}

// ── Quality flips on surfaces ────────────────────────────────────────────────

/// Unit icosphere with every coordinate displaced by a deterministic offset.
fn perturbed_icosphere(subdivisions: usize, amount: f64, seed: u64) -> (Verts, Tris) {
    let (mut verts, tris) = icosphere(subdivisions);
    let mut rng = Lcg(seed);
    for v in &mut verts {
        for coordinate in v.iter_mut() {
            *coordinate += amount * rng.next_signed();
        }
    }
    (verts, tris)
}

#[test]
fn quality_flips_keep_a_closed_sphere_closed_and_outward() {
    // Displacements of up to 0.09 per coordinate (a third of the edge length)
    // make a number of diagonals worth flipping while every triangle still faces
    // outward.
    let (verts, initial) = perturbed_icosphere(2, 0.09, 12);
    assert_outward(&verts, &initial);
    let mut tris = initial.clone();
    let quality_before = mesh_quality_min_angle(&verts, &tris);
    let aspect_before = mesh_aspect_ratio_avg(&verts, &tris);

    let flips = flip_edges_for_quality(&mut tris, &verts);
    assert!(flips > 10, "only {flips} flips on a perturbed icosphere");

    assert_manifold_and_consistent(&tris);
    assert_closed(&tris);
    assert_outward(&verts, &tris);
    assert_eq!(tris.len(), initial.len());
    // The smallest angle never gets worse, and the triangles get rounder.
    assert!(mesh_quality_min_angle(&verts, &tris) >= quality_before);
    assert!(
        mesh_aspect_ratio_avg(&verts, &tris) < aspect_before,
        "{} vs {aspect_before}",
        mesh_aspect_ratio_avg(&verts, &tris)
    );
    assert_eq!(flip_edges_for_quality(&mut tris, &verts), 0);
}

#[test]
fn quality_flips_on_a_sheared_strip_survive_neighbouring_flips() {
    // A strip of parallelogram cells, each split along its long diagonal. All the
    // diagonals flip, and the slanted edges between the cells touch the flipped
    // triangles.
    let cells = 12;
    let shear = 1.5;
    let mut verts: Verts = Vec::new();
    for i in 0..=cells {
        verts.push([i as f64, 0.0, 0.0]);
    }
    for i in 0..=cells {
        verts.push([i as f64 + shear, 1.0, 0.0]);
    }
    let top = |i: usize| cells + 1 + i;
    let mut tris: Tris = Vec::new();
    for i in 0..cells {
        // Long diagonal from the lower left to the upper right corner.
        tris.push([i, i + 1, top(i + 1)]);
        tris.push([i, top(i + 1), top(i)]);
    }
    let initial = tris.clone();
    let boundary = boundary_edges(&tris);
    let area = total_area(&verts, &tris);
    for &t in &initial {
        assert!(normal_z(&verts, t) > 0.0);
    }

    let flips = flip_edges_for_quality(&mut tris, &verts);
    assert!(
        flips >= cells,
        "only {flips} flips for {cells} bad diagonals"
    );

    assert_manifold_and_consistent(&tris);
    assert_eq!(tris.len(), initial.len());
    assert_eq!(boundary_edges(&tris), boundary);
    assert!((total_area(&verts, &tris) - area).abs() < 1e-9);
    for &t in &tris {
        assert!(normal_z(&verts, t) > 0.0, "{t:?} was inverted");
    }
    assert!(mesh_quality_min_angle(&verts, &tris) > mesh_quality_min_angle(&verts, &initial));
    assert_eq!(flip_edges_for_quality(&mut tris, &verts), 0);
}

// ── Valence flips ────────────────────────────────────────────────────────────

#[test]
fn valence_flips_repair_manually_flipped_edges_on_an_icosphere() {
    let (verts, regular) = icosphere(1);
    let regular_deviation = valence_deviation_sum(verts.len(), &regular);
    let valence = valences(verts.len(), &regular);

    // Flip edges between valence-6 vertices (greedily, with disjoint quads) to
    // create valence-5 and valence-7 pairs, using code that shares nothing with
    // the engine.
    let mut tris = regular.clone();
    let mut used: BTreeSet<usize> = BTreeSet::new();
    let mut flipped_edges = 0;
    for (a, b) in undirected_edges(&regular) {
        if flipped_edges == 4 || valence[a] != 6 || valence[b] != 6 {
            continue;
        }
        let quad: BTreeSet<usize> = regular
            .iter()
            .filter(|t| t.contains(&a) && t.contains(&b))
            .flat_map(|t| t.iter().copied())
            .collect();
        if quad.iter().any(|v| used.contains(v)) {
            continue;
        }
        used.extend(quad);
        flip_manually(&mut tris, a, b);
        flipped_edges += 1;
    }
    assert_eq!(flipped_edges, 4);
    let damaged_deviation = valence_deviation_sum(verts.len(), &tris);
    assert!(damaged_deviation > regular_deviation);
    assert_outward(&verts, &tris);

    let flips = flip_for_valence(&mut tris, &verts, no_lock);
    assert!(flips >= 1, "no valence flip happened");
    let repaired_deviation = valence_deviation_sum(verts.len(), &tris);
    assert!(
        repaired_deviation < damaged_deviation,
        "{repaired_deviation} vs {damaged_deviation}"
    );

    assert_manifold_and_consistent(&tris);
    assert_closed(&tris);
    assert_outward(&verts, &tris);
    assert_eq!(tris.len(), regular.len());
    assert_eq!(flip_for_valence(&mut tris, &verts, no_lock), 0);
}

#[test]
fn valence_flips_leave_a_regular_icosphere_alone() {
    let (verts, tris) = icosphere(2);
    let mut work = tris.clone();
    assert_eq!(flip_for_valence(&mut work, &verts, no_lock), 0);
    assert_eq!(work, tris);
}

#[test]
fn locked_edges_are_never_flipped_by_the_valence_engine() {
    // An icosahedron (all valences 5) with one manually flipped edge: the only
    // improving flip undoes it.
    let (verts, regular) = icosphere(0);
    let mut damaged = regular.clone();
    flip_manually(&mut damaged, 0, 11);
    let new_edge = *undirected_edges(&damaged)
        .difference(&undirected_edges(&regular))
        .next()
        .expect("the manual flip created an edge");
    assert_eq!(valence_deviation_sum(verts.len(), &damaged), 16);

    // Locking the new edge (given in either argument order) blocks the repair.
    let lock = |u: usize, v: usize| (u.min(v), u.max(v)) == new_edge;
    let mut locked = damaged.clone();
    assert_eq!(flip_for_valence(&mut locked, &verts, lock), 0);
    assert_eq!(locked, damaged);

    // Locking everything blocks everything.
    let mut all_locked = damaged.clone();
    assert_eq!(flip_for_valence(&mut all_locked, &verts, |_, _| true), 0);
    assert_eq!(all_locked, damaged);

    // Without a lock the icosahedron is restored.
    let mut free = damaged.clone();
    assert!(flip_for_valence(&mut free, &verts, no_lock) >= 1);
    assert_eq!(valence_deviation_sum(verts.len(), &free), 12);
    assert_eq!(undirected_edges(&free), undirected_edges(&regular));
    assert_closed(&free);
    assert_outward(&verts, &free);
}

// ── Remeshing loops ──────────────────────────────────────────────────────────

fn sphere_mesh(subdivisions: usize) -> (TriangleMesh, Tris) {
    let (verts, tris) = icosphere(subdivisions);
    let mesh = TriangleMesh::new(
        verts.iter().map(|v| Vec3::new(v[0], v[1], v[2])).collect(),
        tris.clone(),
    );
    (mesh, tris)
}

#[test]
fn isotropic_remeshing_of_a_sphere_keeps_its_triangles_outward() {
    // With a target length close to the edge length nothing is split or
    // collapsed, so the loop reduces to flips, smoothing and projection and has
    // to leave an outward oriented sphere.
    let (mesh, tris) = sphere_mesh(2);
    let result = isotropic_remesh(&mesh, 0.28, 2);
    let verts: Verts = result.vertices.iter().map(|v| [v.x, v.y, v.z]).collect();
    assert_eq!(result.indices.len(), tris.len());
    assert_outward(&verts, &result.indices);
}
