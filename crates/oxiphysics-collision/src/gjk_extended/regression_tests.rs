// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the polytope expansion, the 2-D GJK distance and the
//! GJK ray cast of [`crate::gjk_extended`].
//!
//! * `ExpandingPolytope::expand` never found a horizon edge (the owning face is
//!   itself in the visible set), so every expansion deleted the visible faces
//!   and replaced them with nothing.
//! * `Gjk2D::distance` never updated its search direction and returned the
//!   nearest simplex vertex, or `0.0` after the iteration limit.
//! * `GjkRaycast::cast` advanced along the ray with a fixed support direction
//!   and decided hits with the support point of an unrelated `+z` direction.

use super::*;
use crate::gjk_epa::simplex_closest::closest_point_to_origin;
use std::collections::BTreeMap;

const TOL: f64 = 1e-6;

fn assert_vec_close(actual: [f64; 3], expected: [f64; 3], tol: f64) {
    assert!(
        len3(sub3(actual, expected)) <= tol,
        "expected {expected:?}, got {actual:?} (tolerance {tol})"
    );
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

    fn vec(&mut self, lo: f64, hi: f64) -> [f64; 3] {
        [self.range(lo, hi), self.range(lo, hi), self.range(lo, hi)]
    }
}

// ── ExpandingPolytope::expand ────────────────────────────────────────────────

/// Eight corners of an axis-aligned cube with the given centre and half extent.
fn cube(center: [f64; 3], half: f64) -> Vec<[f64; 3]> {
    let mut corners = Vec::with_capacity(8);
    for sx in [-1.0, 1.0] {
        for sy in [-1.0, 1.0] {
            for sz in [-1.0, 1.0] {
                corners.push([
                    center[0] + sx * half,
                    center[1] + sy * half,
                    center[2] + sz * half,
                ]);
            }
        }
    }
    corners
}

/// Tetrahedral GJK simplex of `A - B` built from four Minkowski supports.
fn tetrahedron_of(mds: &MinkowskiDiffSupport) -> SimplexManager {
    let mut simplex = SimplexManager::new();
    for dir in [
        [1.0, 1.0, 1.0],
        [-1.0, -1.0, 1.0],
        [-1.0, 1.0, -1.0],
        [1.0, -1.0, -1.0],
    ] {
        let (v, sa, sb) = mds.support(dir);
        simplex.add(v, sa, sb);
    }
    simplex
}

/// Number of faces that use each undirected edge of the polytope.
fn edge_uses(polytope: &ExpandingPolytope) -> BTreeMap<(usize, usize), usize> {
    let mut uses = BTreeMap::new();
    for face in &polytope.faces {
        for k in 0..3 {
            let a = face.indices[k];
            let b = face.indices[(k + 1) % 3];
            *uses.entry((a.min(b), a.max(b))).or_insert(0) += 1;
        }
    }
    uses
}

/// A closed triangulated surface: every edge is shared by exactly two faces and
/// Euler's formula `V - E + F = 2` holds on the vertices that are used.
fn assert_closed_surface(polytope: &ExpandingPolytope) {
    let uses = edge_uses(polytope);
    for (edge, count) in &uses {
        assert_eq!(*count, 2, "edge {edge:?} is used by {count} faces");
    }
    let mut used_vertices: Vec<usize> = polytope
        .faces
        .iter()
        .flat_map(|f| f.indices.iter().copied())
        .collect();
    used_vertices.sort_unstable();
    used_vertices.dedup();
    let v = used_vertices.len() as isize;
    let e = uses.len() as isize;
    let f = polytope.faces.len() as isize;
    assert_eq!(v - e + f, 2, "V={v} E={e} F={f} violates Euler's formula");
}

#[test]
fn expand_replaces_the_visible_faces_with_a_cone_over_the_horizon() {
    let a = cube([0.0; 3], 0.5);
    let b = cube([0.7, 0.0, 0.0], 0.5);
    let mds = MinkowskiDiffSupport::new(a, b);
    let simplex = tetrahedron_of(&mds);
    let mut polytope = ExpandingPolytope::from_simplex(&simplex);
    assert_eq!(polytope.faces.len(), 4);
    assert_closed_surface(&polytope);

    // The corner (0.3, 1, -1) of the Minkowski difference lies outside the
    // initial tetrahedron and is visible from exactly one of its faces.
    let (v, sa, sb) = mds.support([1.0, 1.0, -1.0]);
    assert_vec_close(v, [0.3, 1.0, -1.0], 1e-12);
    let before = polytope.faces.len();
    polytope.expand(v, sa, sb);

    let new_index = polytope.vertices.len() - 1;
    assert!(
        !polytope.faces.is_empty(),
        "the polytope must not be emptied"
    );
    assert!(
        polytope
            .faces
            .iter()
            .any(|f| f.indices.contains(&new_index)),
        "the new vertex must be connected to the horizon"
    );
    // One visible face (3 horizon edges) is replaced by three faces: a cone over
    // the horizon adds two faces to a convex polytope.
    assert_eq!(polytope.faces.len(), before + 2);
    assert_closed_surface(&polytope);
    // Every face keeps an outward normal and a non-negative distance.
    for face in &polytope.faces {
        assert!(face.dist >= 0.0);
        assert!((len3(face.normal) - 1.0).abs() < 1e-12);
    }
}

#[test]
fn expand_with_an_interior_point_leaves_the_faces_alone() {
    let mds = MinkowskiDiffSupport::new(cube([0.0; 3], 0.5), cube([0.7, 0.0, 0.0], 0.5));
    let simplex = tetrahedron_of(&mds);
    let mut polytope = ExpandingPolytope::from_simplex(&simplex);
    let faces_before = polytope.faces.len();
    // The centroid of the tetrahedron is visible from no face.
    polytope.expand([-0.7, 0.0, 0.0], [0.0; 3], [0.0; 3]);
    assert_eq!(polytope.faces.len(), faces_before);
    assert_closed_surface(&polytope);
}

#[test]
fn expand_loop_stays_closed_and_converges_on_the_penetration() {
    // Two unit cubes (half extent 0.5) overlapping by 0.3 along x.
    let mds = MinkowskiDiffSupport::new(cube([0.0; 3], 0.5), cube([0.7, 0.0, 0.0], 0.5));
    let simplex = tetrahedron_of(&mds);
    let mut polytope = ExpandingPolytope::from_simplex(&simplex);

    let mut converged = None;
    for _ in 0..64 {
        assert!(
            !polytope.faces.is_empty(),
            "faces vanished during expansion"
        );
        let face = polytope.closest_face().cloned().expect("a closest face");
        let (v, sa, sb) = mds.support(face.normal);
        if dot3(face.normal, v) - face.dist < 1e-9 {
            converged = Some(face);
            break;
        }
        polytope.expand(v, sa, sb);
        assert_closed_surface(&polytope);
    }
    let face = converged.expect("EPA must converge on two overlapping cubes");
    assert!((face.dist - 0.3).abs() < 1e-9, "depth {}", face.dist);
    assert!(
        (face.normal[0].abs() - 1.0).abs() < 1e-9,
        "normal {:?}",
        face.normal
    );
    assert!(face.normal[1].abs() < 1e-9 && face.normal[2].abs() < 1e-9);
}

#[test]
fn epa_expander_reports_depth_and_axis_for_overlapping_cubes() {
    let mds = MinkowskiDiffSupport::new(cube([0.0; 3], 0.5), cube([0.7, 0.0, 0.0], 0.5));
    let simplex = tetrahedron_of(&mds);
    let mut polytope = ExpandingPolytope::from_simplex(&simplex);
    let expander = EpaExpander::new(64, 1e-9);
    let result = expander.expand(&mut polytope, |d| mds.support(d));
    assert!((result.depth - 0.3).abs() < TOL, "depth {}", result.depth);
    assert!(
        (result.normal[0].abs() - 1.0).abs() < TOL,
        "normal {:?}",
        result.normal
    );
    assert!(result.normal[1].abs() < TOL && result.normal[2].abs() < TOL);
    assert!(!polytope.faces.is_empty());
}

#[test]
fn epa_expander_finds_the_shallowest_axis() {
    // Overlap 0.2 along y but 0.95 along x and z: the minimum translation is along y.
    let mds = MinkowskiDiffSupport::new(cube([0.0; 3], 0.5), cube([0.05, 0.8, -0.05], 0.5));
    let simplex = tetrahedron_of(&mds);
    let mut polytope = ExpandingPolytope::from_simplex(&simplex);
    let result = EpaExpander::new(64, 1e-9).expand(&mut polytope, |d| mds.support(d));
    assert!((result.depth - 0.2).abs() < TOL, "depth {}", result.depth);
    assert!(
        (result.normal[1].abs() - 1.0).abs() < TOL,
        "normal {:?}",
        result.normal
    );
}

// ── Gjk2D::distance ──────────────────────────────────────────────────────────

fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Vec<[f64; 2]> {
    vec![[x0, y0], [x1, y0], [x1, y1], [x0, y1]]
}

#[test]
fn gjk2d_distance_between_offset_squares_is_the_horizontal_gap() {
    // F0133: the old code returned the nearest vertex, 1.118.
    let a = rect(0.0, 0.0, 1.0, 1.0);
    let b = rect(2.0, 0.5, 3.0, 1.5);
    let d = Gjk2D::new(64, 1e-9).distance(&a, &b);
    assert!((d - 1.0).abs() < 1e-8, "distance {d}");
}

#[test]
fn gjk2d_distance_of_tall_polygons_is_not_reported_as_contact() {
    // F0133: the old code ran out of iterations and returned 0.0.
    let a = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 10.0], [0.0, 10.0]];
    let b = vec![[2.0, 10.0], [3.0, 10.0], [3.0, 0.0], [2.0, 0.0]];
    let d = Gjk2D::new(64, 1e-9).distance(&a, &b);
    assert!((d - 1.0).abs() < 1e-8, "distance {d}");
}

#[test]
fn gjk2d_distance_between_separated_squares_on_the_axis() {
    let a = vec![[-2.0_f64, 0.0], [-1.0, 0.0], [-1.0, 1.0], [-2.0, 1.0]];
    let b = vec![[1.0_f64, 0.0], [2.0, 0.0], [2.0, 1.0], [1.0, 1.0]];
    let d = Gjk2D::new(64, 1e-9).distance(&a, &b);
    assert!((d - 2.0).abs() < 1e-8, "distance {d}");
}

#[test]
fn gjk2d_distance_between_corner_facing_squares_is_the_diagonal() {
    let a = rect(0.0, 0.0, 1.0, 1.0);
    let b = rect(2.0, 3.0, 3.0, 4.0);
    let d = Gjk2D::new(64, 1e-9).distance(&a, &b);
    assert!((d - 5.0_f64.sqrt()).abs() < 1e-8, "distance {d}");
}

#[test]
fn gjk2d_distance_to_an_edge_interior() {
    // Triangle apex pointing down at the middle of a square's top edge.
    let square = rect(0.0, 0.0, 2.0, 2.0);
    let triangle = vec![[1.0, 3.5], [0.0, 5.0], [2.0, 5.0]];
    let d = Gjk2D::new(64, 1e-9).distance(&square, &triangle);
    assert!((d - 1.5).abs() < 1e-8, "distance {d}");
}

#[test]
fn gjk2d_overlapping_and_touching_polygons_have_zero_distance() {
    let a = rect(0.0, 0.0, 2.0, 2.0);
    let g = Gjk2D::new(64, 1e-9);
    assert_eq!(
        g.distance(&a, &rect(1.0, 1.0, 3.0, 3.0)),
        0.0,
        "overlapping"
    );
    assert_eq!(g.distance(&a, &rect(0.5, 0.5, 1.5, 1.5)), 0.0, "contained");
    assert_eq!(g.distance(&a, &a), 0.0, "identical");
    let touching = g.distance(&a, &rect(2.0, 0.0, 4.0, 2.0));
    assert!(touching < 1e-9, "edge contact reported distance {touching}");
}

#[test]
fn gjk2d_distance_is_symmetric_and_handles_point_polygons() {
    let hexagon: Vec<[f64; 2]> = (0..6)
        .map(|k| {
            let angle = k as f64 * std::f64::consts::FRAC_PI_3;
            [3.0 + angle.cos(), 2.0 + angle.sin()]
        })
        .collect();
    let square = rect(-1.0, -1.0, 1.0, 1.0);
    let g = Gjk2D::new(64, 1e-9);
    let ab = g.distance(&square, &hexagon);
    let ba = g.distance(&hexagon, &square);
    assert!((ab - ba).abs() < 1e-8, "{ab} vs {ba}");
    // The square corner (1, 1) is nearest to the hexagon edge from (2, 2) to
    // (2.5, 1.134): the foot point is (2.183, 1.683), at distance (1 + sqrt 3) / 2.
    assert!(
        (ab - (1.0 + 3.0_f64.sqrt()) / 2.0).abs() < 1e-8,
        "distance {ab}"
    );
    assert!((ab - polygon_distance_by_enumeration(&square, &hexagon)).abs() < 1e-8);

    // A single point against a square: plain point-box distance.
    let point = vec![[4.0, 5.0]];
    let d = g.distance(&point, &square);
    assert!((d - 5.0).abs() < 1e-9, "distance {d}");

    assert_eq!(g.distance(&[], &square), f64::MAX);
    assert_eq!(g.distance(&square, &[]), f64::MAX);
}

#[test]
fn gjk2d_accuracy_is_relative_so_tiny_polygons_are_measured_exactly() {
    // The stop test compares the remaining gap with the distance itself, so the
    // answer must not depend on the units the polygons are expressed in.
    for scale in [1e-9, 1e-4, 1.0, 1e6] {
        let a = rect(0.0, 0.0, scale, scale);
        let b = rect(3.0 * scale, 0.5 * scale, 4.0 * scale, 1.5 * scale);
        let d = Gjk2D::new(64, 1e-8).distance(&a, &b);
        assert!(
            (d - 2.0 * scale).abs() <= 1e-8 * 2.0 * scale,
            "scale {scale}: distance {d}"
        );
    }
}

#[test]
fn gjk2d_iteration_limit_never_reports_false_contact() {
    // With a single iteration the answer is only an upper bound, but a pair of
    // separated polygons must never be reported as touching.
    let a = rect(0.0, 0.0, 1.0, 1.0);
    let b = rect(2.0, 0.5, 3.0, 1.5);
    let d = Gjk2D::new(1, 1e-9).distance(&a, &b);
    assert!(d >= 1.0 - 1e-12 && d.is_finite(), "distance {d}");
    let d = Gjk2D::new(0, 1e-9).distance(&a, &b);
    assert!(d >= 1.0 - 1e-12 && d.is_finite(), "distance {d}");
}

/// Distance from point `p` to segment `a`-`b` in the plane.
fn point_segment_distance_2d(p: [f64; 2], a: [f64; 2], b: [f64; 2]) -> f64 {
    let ab = [b[0] - a[0], b[1] - a[1]];
    let ap = [p[0] - a[0], p[1] - a[1]];
    let len_sq = ab[0] * ab[0] + ab[1] * ab[1];
    let t = if len_sq > 0.0 {
        ((ap[0] * ab[0] + ap[1] * ab[1]) / len_sq).clamp(0.0, 1.0)
    } else {
        0.0
    };
    let dx = p[0] - (a[0] + t * ab[0]);
    let dy = p[1] - (a[1] + t * ab[1]);
    (dx * dx + dy * dy).sqrt()
}

/// Distance between two disjoint convex polygons by checking every
/// vertex-edge pair in both directions.
fn polygon_distance_by_enumeration(a: &[[f64; 2]], b: &[[f64; 2]]) -> f64 {
    let mut best = f64::MAX;
    for (poly, other) in [(a, b), (b, a)] {
        for p in poly {
            for k in 0..other.len() {
                let d = point_segment_distance_2d(*p, other[k], other[(k + 1) % other.len()]);
                best = best.min(d);
            }
        }
    }
    best
}

/// Separating-axis test for convex counter-clockwise polygons.
fn polygons_overlap(a: &[[f64; 2]], b: &[[f64; 2]]) -> bool {
    for poly in [a, b] {
        for k in 0..poly.len() {
            let p = poly[k];
            let q = poly[(k + 1) % poly.len()];
            let axis = [q[1] - p[1], -(q[0] - p[0])];
            let project = |pts: &[[f64; 2]]| {
                pts.iter()
                    .map(|v| v[0] * axis[0] + v[1] * axis[1])
                    .fold((f64::MAX, f64::MIN), |(lo, hi), x| (lo.min(x), hi.max(x)))
            };
            let (a_lo, a_hi) = project(a);
            let (b_lo, b_hi) = project(b);
            if a_hi < b_lo || b_hi < a_lo {
                return false;
            }
        }
    }
    true
}

/// Regular `n`-gon with the given centre, circumradius and rotation.
fn regular_polygon(center: [f64; 2], radius: f64, n: usize, phase: f64) -> Vec<[f64; 2]> {
    (0..n)
        .map(|k| {
            let angle = phase + k as f64 * std::f64::consts::TAU / n as f64;
            [
                center[0] + radius * angle.cos(),
                center[1] + radius * angle.sin(),
            ]
        })
        .collect()
}

#[test]
fn gjk2d_matches_vertex_edge_enumeration_on_random_convex_polygons() {
    let mut rng = Rng(0xFEED_FACE_CAFE_BEEF);
    let g = Gjk2D::new(128, 1e-10);
    let (mut disjoint, mut overlapping) = (0, 0);
    for _ in 0..400 {
        let a = regular_polygon(
            [rng.range(-4.0, 4.0), rng.range(-4.0, 4.0)],
            rng.range(0.3, 2.0),
            3 + (rng.unit() * 6.0) as usize,
            rng.range(0.0, std::f64::consts::TAU),
        );
        let b = regular_polygon(
            [rng.range(-4.0, 4.0), rng.range(-4.0, 4.0)],
            rng.range(0.3, 2.0),
            3 + (rng.unit() * 6.0) as usize,
            rng.range(0.0, std::f64::consts::TAU),
        );
        let d = g.distance(&a, &b);
        if polygons_overlap(&a, &b) {
            assert_eq!(d, 0.0, "overlapping polygons must have zero distance");
            overlapping += 1;
        } else {
            let expected = polygon_distance_by_enumeration(&a, &b);
            assert!((d - expected).abs() < 1e-7, "distance {d} vs {expected}");
            disjoint += 1;
        }
    }
    assert!(
        disjoint > 50 && overlapping > 50,
        "{disjoint} disjoint, {overlapping} overlapping"
    );
}

#[test]
fn gjk2d_matches_the_axis_aligned_box_gap_on_random_rectangles() {
    let mut rng = Rng(0x1234_5678_9ABC_DEF1);
    let g = Gjk2D::new(64, 1e-10);
    for _ in 0..300 {
        let (ax, ay) = (rng.range(-5.0, 5.0), rng.range(-5.0, 5.0));
        let (aw, ah) = (rng.range(0.2, 2.0), rng.range(0.2, 2.0));
        let (bx, by) = (rng.range(-5.0, 5.0), rng.range(-5.0, 5.0));
        let (bw, bh) = (rng.range(0.2, 2.0), rng.range(0.2, 2.0));
        let a = rect(ax, ay, ax + aw, ay + ah);
        let b = rect(bx, by, bx + bw, by + bh);
        let gap_x = (bx - (ax + aw)).max(ax - (bx + bw)).max(0.0);
        let gap_y = (by - (ay + ah)).max(ay - (by + bh)).max(0.0);
        let expected = (gap_x * gap_x + gap_y * gap_y).sqrt();
        let d = g.distance(&a, &b);
        assert!((d - expected).abs() < 1e-7, "distance {d} vs {expected}");
    }
}

// ── GjkRaycast::cast ─────────────────────────────────────────────────────────

fn sphere_support(center: [f64; 3], radius: f64) -> impl Fn([f64; 3]) -> [f64; 3] {
    move |d| add3(center, scale3(norm3(d), radius))
}

fn box_support(center: [f64; 3], half: [f64; 3]) -> impl Fn([f64; 3]) -> [f64; 3] {
    move |d| {
        [
            center[0] + half[0] * d[0].signum(),
            center[1] + half[1] * d[1].signum(),
            center[2] + half[2] * d[2].signum(),
        ]
    }
}

fn hull_support(vertices: Vec<[f64; 3]>) -> impl Fn([f64; 3]) -> [f64; 3] {
    move |d| MinkowskiDiffSupport::support_point(&vertices, d)
}

fn caster() -> GjkRaycast {
    GjkRaycast::new(128, 1e-9)
}

#[test]
fn ray_hits_a_sphere_head_on() {
    let r = caster().cast(
        [10.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0],
        100.0,
        sphere_support([0.0; 3], 1.0),
    );
    let t = r.t.expect("the ray must hit");
    assert!((t - 9.0).abs() < TOL, "t={t}");
    assert_vec_close(r.point, [1.0, 0.0, 0.0], TOL);
    assert_vec_close(r.normal, [1.0, 0.0, 0.0], TOL);
}

#[test]
fn ray_hits_a_sphere_straight_ahead_along_z() {
    // F0135: unit sphere at (0,0,5), ray from the origin along +z. The old code
    // measured the support point against an unrelated +z direction and missed.
    let r = caster().cast(
        [0.0; 3],
        [0.0, 0.0, 1.0],
        100.0,
        sphere_support([0.0, 0.0, 5.0], 1.0),
    );
    let t = r.t.expect("the ray must hit");
    assert!((t - 4.0).abs() < TOL, "t={t}");
    assert_vec_close(r.point, [0.0, 0.0, 4.0], TOL);
    assert_vec_close(r.normal, [0.0, 0.0, -1.0], TOL);
}

#[test]
fn ray_direction_is_normalised_so_t_is_a_distance() {
    let r = caster().cast(
        [0.0; 3],
        [0.0, 0.0, 7.5],
        100.0,
        sphere_support([0.0, 0.0, 5.0], 1.0),
    );
    assert!((r.t.expect("hit") - 4.0).abs() < TOL);
}

#[test]
fn oblique_ray_hits_a_sphere_at_the_analytic_point() {
    let r = caster().cast(
        [-5.0, -5.0, 0.0],
        [2.0, 2.0, 0.0],
        100.0,
        sphere_support([0.0; 3], 1.0),
    );
    let t = r.t.expect("the ray must hit");
    assert!((t - (50.0_f64.sqrt() - 1.0)).abs() < TOL, "t={t}");
    let s = std::f64::consts::FRAC_1_SQRT_2;
    assert_vec_close(r.point, [-s, -s, 0.0], TOL);
    assert_vec_close(r.normal, [-s, -s, 0.0], TOL);
}

#[test]
fn off_centre_ray_hits_a_sphere_with_the_matching_normal() {
    // Ray at height 0.6 hits the unit sphere at x = 0.8: normal (0.8, 0.6, 0).
    let r = caster().cast(
        [10.0, 0.6, 0.0],
        [-1.0, 0.0, 0.0],
        100.0,
        sphere_support([0.0; 3], 1.0),
    );
    let t = r.t.expect("the ray must hit");
    assert!((t - 9.2).abs() < TOL, "t={t}");
    assert_vec_close(r.point, [0.8, 0.6, 0.0], TOL);
    // On a curved shape the normal is the direction of the last supporting plane
    // and is accurate to about the square root of the tolerance.
    assert_vec_close(r.normal, [0.8, 0.6, 0.0], 1e-4);
}

#[test]
fn ray_hits_a_box_face_head_on() {
    // With signum support the same corner is returned again on a later
    // iteration; that repetition must not be mistaken for convergence.
    let r = caster().cast(
        [5.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0],
        100.0,
        box_support([0.0; 3], [1.0; 3]),
    );
    let t = r.t.expect("the ray must hit");
    assert!((t - 4.0).abs() < TOL, "t={t}");
    assert_vec_close(r.normal, [1.0, 0.0, 0.0], TOL);
    assert_vec_close(r.point, [1.0, 0.0, 0.0], TOL);
}

#[test]
fn ray_hits_a_box_face_off_centre() {
    let r = caster().cast(
        [5.0, 0.3, 0.2],
        [-1.0, 0.0, 0.0],
        100.0,
        box_support([0.0; 3], [1.0; 3]),
    );
    let t = r.t.expect("the ray must hit");
    assert!((t - 4.0).abs() < TOL, "t={t}");
    assert_vec_close(r.point, [1.0, 0.3, 0.2], TOL);
    assert_vec_close(r.normal, [1.0, 0.0, 0.0], TOL);
}

#[test]
fn diagonal_ray_enters_a_box_through_the_nearest_face() {
    // From (3, 3, 0.5) towards (-1, -1.2, 0): reaches the x = 1 plane first.
    let dir = [-1.0, -1.2, 0.0];
    let r = caster().cast([3.0, 3.0, 0.5], dir, 100.0, box_support([0.0; 3], [1.0; 3]));
    let t = r.t.expect("the ray must hit");
    let len = len3(dir);
    // x reaches 1 after parameter 2 / 1 of the un-normalised direction.
    assert!((t - 2.0 * len).abs() < TOL, "t={t}");
    assert_vec_close(r.point, [1.0, 3.0 - 2.4, 0.5], TOL);
    assert_vec_close(r.normal, [1.0, 0.0, 0.0], TOL);
}

#[test]
fn ray_hits_a_tetrahedron_through_its_slanted_face() {
    let verts = vec![
        [0.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [0.0, 2.0, 0.0],
        [0.0, 0.0, 2.0],
    ];
    let r = caster().cast(
        [0.5, 0.5, 5.0],
        [0.0, 0.0, -1.0],
        100.0,
        hull_support(verts),
    );
    let t = r.t.expect("the ray must hit");
    // The slanted face is x + y + z = 2, so z = 1 above (0.5, 0.5).
    assert!((t - 4.0).abs() < TOL, "t={t}");
    let s = 1.0 / 3.0_f64.sqrt();
    assert_vec_close(r.normal, [s, s, s], TOL);
}

#[test]
fn point_shape_is_hit_at_its_location() {
    let r = caster().cast([10.0, 0.0, 0.0], [-1.0, 0.0, 0.0], 20.0, |_d| [0.0; 3]);
    assert!((r.t.expect("the ray passes through the point") - 10.0).abs() < TOL);
    assert_vec_close(r.point, [0.0; 3], TOL);
    assert_vec_close(r.normal, [1.0, 0.0, 0.0], TOL);
}

#[test]
fn ray_that_passes_beside_a_sphere_misses() {
    let c = caster();
    assert!(
        c.cast(
            [10.0, 1.5, 0.0],
            [-1.0, 0.0, 0.0],
            100.0,
            sphere_support([0.0; 3], 1.0)
        )
        .t
        .is_none()
    );
    assert!(
        c.cast(
            [10.0, 0.0, 2.0],
            [-1.0, 0.0, 0.0],
            100.0,
            sphere_support([0.0; 3], 1.0)
        )
        .t
        .is_none()
    );
}

#[test]
fn ray_that_passes_beside_a_box_misses() {
    let c = caster();
    assert!(
        c.cast(
            [5.0, 1.2, 0.0],
            [-1.0, 0.0, 0.0],
            100.0,
            box_support([0.0; 3], [1.0; 3])
        )
        .t
        .is_none()
    );
    assert!(
        c.cast(
            [5.0, 0.0, -1.5],
            [-1.0, 0.0, 0.0],
            100.0,
            box_support([0.0; 3], [1.0; 3])
        )
        .t
        .is_none()
    );
}

#[test]
fn nearly_grazing_ray_still_hits_the_sphere() {
    // Height 0.99 on the unit sphere: x = sqrt(1 - 0.9801).
    let r = caster().cast(
        [10.0, 0.99, 0.0],
        [-1.0, 0.0, 0.0],
        100.0,
        sphere_support([0.0; 3], 1.0),
    );
    let t = r.t.expect("the ray must hit");
    let x_hit = (1.0_f64 - 0.99 * 0.99).sqrt();
    assert!((t - (10.0 - x_hit)).abs() < 1e-5, "t={t}");
    assert_vec_close(r.normal, [x_hit, 0.99, 0.0], 1e-5);
}

#[test]
fn ray_pointing_away_from_the_shape_misses() {
    let r = caster().cast(
        [10.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        100.0,
        sphere_support([0.0; 3], 1.0),
    );
    assert!(r.t.is_none());
}

#[test]
fn hit_beyond_max_t_is_a_miss() {
    let support = sphere_support([0.0; 3], 1.0);
    let c = caster();
    assert!(
        c.cast([10.0, 0.0, 0.0], [-1.0, 0.0, 0.0], 5.0, &support)
            .t
            .is_none()
    );
    let hit = c.cast([10.0, 0.0, 0.0], [-1.0, 0.0, 0.0], 9.5, &support);
    assert!((hit.t.expect("within reach") - 9.0).abs() < TOL);
}

#[test]
fn ray_starting_inside_hits_at_zero_with_the_reversed_ray_as_normal() {
    let r = caster().cast(
        [0.2, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        100.0,
        sphere_support([0.0; 3], 1.0),
    );
    assert_eq!(r.t, Some(0.0));
    assert_vec_close(r.point, [0.2, 0.0, 0.0], 1e-12);
    assert_vec_close(r.normal, [0.0, -1.0, 0.0], 1e-12);
}

#[test]
fn degenerate_rays_miss() {
    let support = sphere_support([0.0; 3], 1.0);
    let c = caster();
    assert!(
        c.cast([10.0, 0.0, 0.0], [0.0; 3], 100.0, &support)
            .t
            .is_none()
    );
    assert!(
        c.cast([10.0, 0.0, 0.0], [-1.0, 0.0, 0.0], -1.0, &support)
            .t
            .is_none()
    );
    assert!(
        c.cast([10.0, 0.0, 0.0], [-1.0, 0.0, 0.0], f64::NAN, &support)
            .t
            .is_none()
    );
    assert!(
        c.cast([10.0, 0.0, 0.0], [f64::NAN, 0.0, 0.0], 100.0, &support)
            .t
            .is_none()
    );
}

/// Analytic ray-sphere entry distance (`None` for a miss or a hit behind the origin).
fn sphere_entry(
    origin: [f64; 3],
    dir: [f64; 3],
    center: [f64; 3],
    radius: f64,
) -> Option<(f64, f64)> {
    let d = norm3(dir);
    let oc = sub3(origin, center);
    let b = dot3(d, oc);
    let c = dot3(oc, oc) - radius * radius;
    let disc = b * b - c;
    if disc < 0.0 {
        return None;
    }
    let t = -b - disc.sqrt();
    if t < 0.0 { None } else { Some((t, disc)) }
}

/// Analytic ray-box entry distance by the slab method.
fn box_entry(origin: [f64; 3], dir: [f64; 3], half: [f64; 3]) -> Option<f64> {
    let d = norm3(dir);
    let mut t_min = f64::NEG_INFINITY;
    let mut t_max = f64::INFINITY;
    for k in 0..3 {
        if d[k].abs() < 1e-12 {
            if origin[k].abs() > half[k] {
                return None;
            }
        } else {
            let (t0, t1) = ((-half[k] - origin[k]) / d[k], (half[k] - origin[k]) / d[k]);
            t_min = t_min.max(t0.min(t1));
            t_max = t_max.min(t0.max(t1));
        }
    }
    if t_min > t_max || t_max < 0.0 || t_min < 0.0 {
        None
    } else {
        Some(t_min)
    }
}

#[test]
fn random_rays_agree_with_analytic_sphere_and_box_intersections() {
    let mut rng = Rng(0x0DDB_1A5E_5BAD_5EED);
    let c = caster();
    let (mut sphere_hits, mut box_hits, mut misses) = (0, 0, 0);
    for _ in 0..400 {
        let origin = rng.vec(-6.0, 6.0);
        let target = rng.vec(-1.5, 1.5);
        let dir = sub3(target, origin);
        if len3(dir) < 0.5 {
            continue;
        }

        // Unit sphere at the origin.
        match sphere_entry(origin, dir, [0.0; 3], 1.0) {
            Some((t, disc)) if disc > 1e-2 && t > 1e-2 => {
                let r = c.cast(origin, dir, 1000.0, sphere_support([0.0; 3], 1.0));
                let got = r.t.expect("analytic hit must be found");
                assert!((got - t).abs() < 1e-6, "sphere t={got} vs {t}");
                assert!(
                    (len3(r.point) - 1.0).abs() < 1e-6,
                    "hit point {:?} is not on the sphere",
                    r.point
                );
                assert!((len3(r.normal) - 1.0).abs() < 1e-9);
                sphere_hits += 1;
            }
            None => {
                let oc = origin;
                let d = norm3(dir);
                let b = dot3(d, oc);
                let miss_decisively =
                    b * b - (dot3(oc, oc) - 1.0) < -1e-2 || (b > 0.0 && dot3(oc, oc) > 1.0);
                if miss_decisively {
                    assert!(
                        c.cast(origin, dir, 1000.0, sphere_support([0.0; 3], 1.0))
                            .t
                            .is_none()
                    );
                    misses += 1;
                }
            }
            _ => {}
        }

        // Box of half extents (1, 0.5, 0.8) at the origin.
        let half = [1.0, 0.5, 0.8];
        if let Some(t) = box_entry(origin, dir, half)
            && t > 1e-2
        {
            let r = c.cast(origin, dir, 1000.0, box_support([0.0; 3], half));
            let got = r.t.expect("analytic box hit must be found");
            assert!((got - t).abs() < 1e-6, "box t={got} vs {t}");
            assert!((len3(r.normal) - 1.0).abs() < 1e-9);
            box_hits += 1;
        }
    }
    assert!(
        sphere_hits > 20 && box_hits > 20 && misses > 20,
        "{sphere_hits} {box_hits} {misses}"
    );
}

// ── randomised comparison with exact answers for convex polyhedra ────────────

/// Largest separating-axis gap of the convex hulls of two point sets together
/// with the unit axis that attains it.
///
/// The candidates are the normals of all triangles spanned by three points of
/// either set and the cross products of all point-pair directions, a superset
/// of the face normals and edge pairs of the hulls. The gap is positive exactly
/// when the hulls are disjoint; for overlapping hulls `-gap` is the exact
/// penetration depth and the axis the direction of the minimum translation.
fn best_separation(a: &[[f64; 3]], b: &[[f64; 3]]) -> (f64, [f64; 3]) {
    let mut axes = Vec::new();
    for cloud in [a, b] {
        for i in 0..cloud.len() {
            for j in (i + 1)..cloud.len() {
                for k in (j + 1)..cloud.len() {
                    axes.push(cross3(sub3(cloud[j], cloud[i]), sub3(cloud[k], cloud[i])));
                }
            }
        }
    }
    for i in 0..a.len() {
        for j in (i + 1)..a.len() {
            for k in 0..b.len() {
                for l in (k + 1)..b.len() {
                    axes.push(cross3(sub3(a[j], a[i]), sub3(b[l], b[k])));
                }
            }
        }
    }
    let interval = |cloud: &[[f64; 3]], u: [f64; 3]| {
        cloud
            .iter()
            .map(|p| dot3(*p, u))
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), x| {
                (lo.min(x), hi.max(x))
            })
    };
    let mut best = (f64::NEG_INFINITY, [0.0; 3]);
    for axis in axes {
        if len3(axis) < 1e-9 {
            continue;
        }
        let u = norm3(axis);
        let (a_lo, a_hi) = interval(a, u);
        let (b_lo, b_hi) = interval(b, u);
        let gap = (b_lo - a_hi).max(a_lo - b_hi);
        if gap > best.0 {
            best = (gap, u);
        }
    }
    best
}

/// A tetrahedron of Minkowski-difference supports (along random directions)
/// that encloses the origin, the starting simplex EPA needs.
fn enclosing_tetrahedron(mds: &MinkowskiDiffSupport, rng: &mut Rng) -> Option<SimplexManager> {
    for _ in 0..200 {
        let supports: Vec<_> = (0..4).map(|_| mds.support(rng.vec(-1.0, 1.0))).collect();
        let points = [supports[0].0, supports[1].0, supports[2].0, supports[3].0];
        if closest_point_to_origin(&points).encloses_origin {
            let mut simplex = SimplexManager::new();
            for (v, sa, sb) in supports {
                simplex.add(v, sa, sb);
            }
            return Some(simplex);
        }
    }
    None
}

#[test]
fn epa_matches_the_exact_penetration_depth_of_random_convex_polyhedra() {
    // Random polyhedra expose multi-face visible regions (horizons of four or
    // more edges), which the axis-aligned cube cases never produce.
    let mut rng = Rng(0x0A11_CE5E_ED51_7E57);
    let mut verified = 0;
    for case in 0..300 {
        let count_a = 5 + (rng.unit() * 4.0) as usize;
        let count_b = 5 + (rng.unit() * 4.0) as usize;
        let offset = rng.vec(-0.6, 0.6);
        let cloud_a: Vec<[f64; 3]> = (0..count_a).map(|_| rng.vec(-1.0, 1.0)).collect();
        let cloud_b: Vec<[f64; 3]> = (0..count_b)
            .map(|_| add3(offset, rng.vec(-1.0, 1.0)))
            .collect();
        let (gap, axis) = best_separation(&cloud_a, &cloud_b);
        if gap > -0.05 {
            continue; // disjoint, or too shallow to start EPA robustly
        }
        let mds = MinkowskiDiffSupport::new(cloud_a, cloud_b);
        let Some(simplex) = enclosing_tetrahedron(&mds, &mut rng) else {
            continue;
        };
        let mut polytope = ExpandingPolytope::from_simplex(&simplex);
        if polytope.closest_face().is_none_or(|face| face.dist < 1e-3) {
            continue;
        }
        let result = EpaExpander::new(256, 1e-10).expand(&mut polytope, |d| mds.support(d));
        assert!(
            (result.depth + gap).abs() < 1e-6,
            "case {case}: depth {} vs exact {}",
            result.depth,
            -gap
        );
        assert!(
            dot3(result.normal, axis).abs() > 1.0 - 1e-6,
            "case {case}: normal {:?} vs axis {axis:?}",
            result.normal
        );
        assert!(!polytope.faces.is_empty());
        assert_closed_surface(&polytope);
        verified += 1;
    }
    assert!(
        verified > 60,
        "only {verified} configurations were decisive"
    );
}

/// Möller-Trumbore: distance along the unit direction `dir` to the triangle
/// `tri`, accepting barycentric coordinates down to `-slack`. Also returns the
/// triangle normal facing against the ray.
fn ray_triangle_hit(
    origin: [f64; 3],
    dir: [f64; 3],
    tri: [[f64; 3]; 3],
    slack: f64,
) -> Option<(f64, [f64; 3])> {
    let e1 = sub3(tri[1], tri[0]);
    let e2 = sub3(tri[2], tri[0]);
    let p = cross3(dir, e2);
    let det = dot3(e1, p);
    if det.abs() < 1e-12 {
        return None;
    }
    let inv = 1.0 / det;
    let s = sub3(origin, tri[0]);
    let u = dot3(s, p) * inv;
    let q = cross3(s, e1);
    let v = dot3(dir, q) * inv;
    if u < -slack || v < -slack || u + v > 1.0 + slack {
        return None;
    }
    let t = dot3(e2, q) * inv;
    if t < 0.0 {
        return None;
    }
    let normal = norm3(cross3(e1, e2));
    let facing = if dot3(normal, dir) > 0.0 {
        neg3(normal)
    } else {
        normal
    };
    Some((t, facing))
}

/// First hit of a ray that starts outside the convex hull of `cloud`, found by
/// testing the ray against every triangle spanned by three of the points. A
/// triangle through the interior is always entered after the boundary, so the
/// nearest hit of all of them is the boundary hit.
fn first_hull_hit(
    origin: [f64; 3],
    dir: [f64; 3],
    cloud: &[[f64; 3]],
    slack: f64,
) -> Option<(f64, [f64; 3])> {
    let mut best: Option<(f64, [f64; 3])> = None;
    for i in 0..cloud.len() {
        for j in (i + 1)..cloud.len() {
            for k in (j + 1)..cloud.len() {
                let tri = [cloud[i], cloud[j], cloud[k]];
                if let Some(hit) = ray_triangle_hit(origin, dir, tri, slack)
                    && best.is_none_or(|(t, _)| hit.0 < t)
                {
                    best = Some(hit);
                }
            }
        }
    }
    best
}

#[test]
fn random_rays_agree_with_triangle_enumeration_on_random_polyhedra() {
    let mut rng = Rng(0x7A11_5EED_C0DE_0042);
    let c = caster();
    let (mut hits, mut misses) = (0, 0);
    for case in 0..600 {
        let count = 5 + (rng.unit() * 4.0) as usize;
        let cloud: Vec<[f64; 3]> = (0..count).map(|_| rng.vec(-1.0, 1.0)).collect();
        let direction_to_origin = rng.vec(-1.0, 1.0);
        if len3(direction_to_origin) < 0.1 {
            continue;
        }
        let origin = scale3(norm3(direction_to_origin), 5.0);
        let dir = sub3(rng.vec(-1.2, 1.2), origin);
        let ray = norm3(dir);
        // Decisive only when widening and shrinking every triangle by the same
        // margin gives the same answer (no grazing hit of an edge or vertex).
        let loose = first_hull_hit(origin, ray, &cloud, 1e-3);
        let tight = first_hull_hit(origin, ray, &cloud, -1e-3);
        let expected = match (loose, tight) {
            (None, None) => None,
            (Some(l), Some(t)) if (l.0 - t.0).abs() < 1e-9 => Some(t),
            _ => continue,
        };
        let result = c.cast(origin, dir, 100.0, hull_support(cloud.clone()));
        match expected {
            None => {
                assert!(
                    result.t.is_none(),
                    "case {case}: the ray misses the hull but t={:?}",
                    result.t
                );
                misses += 1;
            }
            Some((t, normal)) => {
                let got = result
                    .t
                    .unwrap_or_else(|| panic!("case {case}: the ray hits at t={t}"));
                assert!((got - t).abs() < 1e-6, "case {case}: t={got} vs {t}");
                assert_vec_close(result.point, add3(origin, scale3(ray, got)), 1e-9);
                assert_vec_close(result.normal, normal, 1e-6);
                hits += 1;
            }
        }
    }
    assert!(
        hits > 100 && misses > 40,
        "{hits} hits and {misses} misses were decisive"
    );
}
