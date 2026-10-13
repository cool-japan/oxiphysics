// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Closest point of a small simplex to the origin (Johnson's sub-algorithm).
//!
//! GJK needs, in every iteration, the point of the convex hull of at most four
//! vertices that is nearest to the origin together with the barycentric weights
//! that express it. Vertices whose weight is exactly zero do not support the
//! closest point and can be dropped from the simplex.
//!
//! The solver is shared by the GJK distance query ([`super::distance`]), the
//! 2-D GJK and the GJK ray cast of [`crate::gjk_extended`].
//!
//! # Degenerate input
//!
//! GJK regularly produces duplicate, collinear or coplanar vertices (a box
//! support function returns corners, so repeated points are the norm). Every
//! case is handled explicitly instead of dividing by a vanishing quantity:
//!
//! * a segment with coincident end points returns that point,
//! * a zero-area (or numerically flat) triangle returns the best of its three
//!   edges,
//! * a near-flat tetrahedron never reports that it encloses the origin; the
//!   best of its four faces is returned instead.

use super::{add3, cross3, dot3, len3, neg3, scale3, sub3};

/// A triangle whose squared double area is below this fraction of the fourth
/// power of its longest edge is treated as degenerate (edges closer than about
/// `1e-8` radians to parallel).
const FLAT_TRIANGLE_REL: f64 = 1e-16;

/// A tetrahedron whose triple product is below this fraction of the product of
/// its three edge lengths at one vertex is treated as flat.
const FLAT_TETRAHEDRON_REL: f64 = 1e-10;

/// Closest point of a simplex to the origin together with its barycentric
/// description.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct SimplexClosest {
    /// Point of the simplex nearest to the origin.
    pub point: [f64; 3],
    /// Barycentric weight of every input vertex, in input order.
    ///
    /// The weights are non-negative and sum to one. Vertices that do not
    /// support the closest point have a weight of exactly `0.0`; entries
    /// beyond the number of input vertices are `0.0`.
    pub weights: [f64; 4],
    /// `true` when the origin lies inside (or on the boundary of) a
    /// non-degenerate tetrahedron. Only four input vertices can enclose the
    /// origin.
    pub encloses_origin: bool,
}

impl SimplexClosest {
    /// Closest point on a feature spanned by `weights` (no enclosure).
    fn on_feature(point: [f64; 3], weights: [f64; 4]) -> Self {
        Self {
            point,
            weights,
            encloses_origin: false,
        }
    }
}

/// Closest point of the convex hull of `points` to the origin.
///
/// At most the first four points are used. An empty slice yields the origin
/// with all-zero weights.
pub(crate) fn closest_point_to_origin(points: &[[f64; 3]]) -> SimplexClosest {
    match points {
        [] => SimplexClosest::on_feature([0.0; 3], [0.0; 4]),
        [a] => SimplexClosest::on_feature(*a, [1.0, 0.0, 0.0, 0.0]),
        [a, b] => {
            let (point, w) = segment_closest(*a, *b);
            SimplexClosest::on_feature(point, [w[0], w[1], 0.0, 0.0])
        }
        [a, b, c] => {
            let (point, w) = triangle_closest(*a, *b, *c);
            SimplexClosest::on_feature(point, [w[0], w[1], w[2], 0.0])
        }
        [a, b, c, d, ..] => tetrahedron_closest(*a, *b, *c, *d),
    }
}

/// Closest point on segment `ab` to the origin and its weights `[wa, wb]`.
///
/// Coincident end points return `a`.
fn segment_closest(a: [f64; 3], b: [f64; 3]) -> ([f64; 3], [f64; 2]) {
    let ab = sub3(b, a);
    let length_sq = dot3(ab, ab);
    // (origin - a) . ab : positive when the origin projects beyond `a` towards `b`.
    let along = -dot3(a, ab);
    if along <= 0.0 {
        return (a, [1.0, 0.0]);
    }
    if along >= length_sq {
        return (b, [0.0, 1.0]);
    }
    let t = along / length_sq;
    (add3(a, scale3(ab, t)), [1.0 - t, t])
}

/// Closest point on triangle `abc` to the origin and its weights.
///
/// The origin is projected onto the plane of the triangle; when the projection
/// falls inside, it is the answer. Otherwise (or when the triangle is
/// degenerate) the nearest point lies on an edge.
fn triangle_closest(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    let ab = sub3(b, a);
    let ac = sub3(c, a);
    let bc = sub3(c, b);
    let normal = cross3(ab, ac);
    let normal_sq = dot3(normal, normal);
    let longest_sq = dot3(ab, ab).max(dot3(ac, ac)).max(dot3(bc, bc));
    let flat_limit = FLAT_TRIANGLE_REL * longest_sq * longest_sq;
    if normal_sq.is_finite() && normal_sq > flat_limit {
        // Normal equations of the projection of `origin - a` onto span{ab, ac},
        // translated to `a` so that the arithmetic stays well conditioned for
        // triangles far from the origin.
        let ao = neg3(a);
        let d00 = dot3(ab, ab);
        let d01 = dot3(ab, ac);
        let d11 = dot3(ac, ac);
        let r0 = dot3(ab, ao);
        let r1 = dot3(ac, ao);
        let s = (d11 * r0 - d01 * r1) / normal_sq;
        let t = (d00 * r1 - d01 * r0) / normal_sq;
        let u = 1.0 - s - t;
        if s >= 0.0 && t >= 0.0 && u >= 0.0 {
            let point = add3(a, add3(scale3(ab, s), scale3(ac, t)));
            return (point, [u, s, t]);
        }
    }
    best_edge_of_triangle(a, b, c)
}

/// Nearest point to the origin on the boundary of triangle `abc`.
fn best_edge_of_triangle(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    let (p_ab, w_ab) = segment_closest(a, b);
    let (p_ac, w_ac) = segment_closest(a, c);
    let (p_bc, w_bc) = segment_closest(b, c);
    let mut best = (p_ab, [w_ab[0], w_ab[1], 0.0]);
    let mut best_sq = dot3(p_ab, p_ab);
    let sq_ac = dot3(p_ac, p_ac);
    if sq_ac < best_sq {
        best_sq = sq_ac;
        best = (p_ac, [w_ac[0], 0.0, w_ac[1]]);
    }
    let sq_bc = dot3(p_bc, p_bc);
    if sq_bc < best_sq {
        best = (p_bc, [0.0, w_bc[0], w_bc[1]]);
    }
    best
}

/// Closest point of tetrahedron `abcd` to the origin.
fn tetrahedron_closest(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> SimplexClosest {
    let ab = sub3(b, a);
    let ac = sub3(c, a);
    let ad = sub3(d, a);
    let ac_x_ad = cross3(ac, ad);
    let det = dot3(ab, ac_x_ad);
    let edge_product = len3(ab) * len3(ac) * len3(ad);
    if det.is_finite() && det.abs() > FLAT_TETRAHEDRON_REL * edge_product {
        // Barycentric coordinates of the origin (Cramer's rule on
        // `origin - a = lb ab + lc ac + ld ad`).
        let ao = neg3(a);
        let lb = dot3(ao, ac_x_ad) / det;
        let lc = dot3(ab, cross3(ao, ad)) / det;
        let ld = dot3(ab, cross3(ac, ao)) / det;
        let la = 1.0 - lb - lc - ld;
        if la >= 0.0 && lb >= 0.0 && lc >= 0.0 && ld >= 0.0 {
            let sum = la + lb + lc + ld;
            return SimplexClosest {
                point: [0.0; 3],
                weights: [la / sum, lb / sum, lc / sum, ld / sum],
                encloses_origin: true,
            };
        }
    }
    // The origin is outside (or the tetrahedron is flat): the nearest point
    // lies on the union of the four faces.
    let vertices = [a, b, c, d];
    let faces = [[0_usize, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]];
    let mut best_point = a;
    let mut best_weights = [1.0, 0.0, 0.0, 0.0];
    let mut best_sq = f64::INFINITY;
    for face in faces {
        let (point, w) = triangle_closest(vertices[face[0]], vertices[face[1]], vertices[face[2]]);
        let dist_sq = dot3(point, point);
        if dist_sq < best_sq {
            best_sq = dist_sq;
            best_point = point;
            best_weights = [0.0; 4];
            for (slot, &vertex_index) in face.iter().enumerate() {
                best_weights[vertex_index] = w[slot];
            }
        }
    }
    SimplexClosest::on_feature(best_point, best_weights)
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-12;

    fn weight_sum(c: &SimplexClosest) -> f64 {
        c.weights.iter().sum()
    }

    /// Reconstruct the closest point from the barycentric weights.
    fn blend(points: &[[f64; 3]], weights: &[f64; 4]) -> [f64; 3] {
        let mut p = [0.0; 3];
        for (v, &w) in points.iter().zip(weights.iter()) {
            p = add3(p, scale3(*v, w));
        }
        p
    }

    fn assert_close(a: [f64; 3], b: [f64; 3], tol: f64) {
        assert!(
            len3(sub3(a, b)) <= tol,
            "expected {b:?}, got {a:?} (tolerance {tol})"
        );
    }

    // ── degenerate inputs (verified before any caller relies on the solver) ───

    #[test]
    fn coincident_segment_returns_the_shared_point() {
        let p = [2.0, -1.0, 0.5];
        let c = closest_point_to_origin(&[p, p]);
        assert_eq!(c.point, p);
        assert_eq!(c.weights, [1.0, 0.0, 0.0, 0.0]);
        assert!(!c.encloses_origin);
    }

    #[test]
    fn zero_area_triangle_with_coincident_vertices_uses_the_remaining_edge() {
        // a == b: the triangle is the segment a-c. The origin projects inside it.
        let a = [1.0, 1.0, 0.0];
        let c = [-1.0, 1.0, 0.0];
        let r = closest_point_to_origin(&[a, a, c]);
        assert_close(r.point, [0.0, 1.0, 0.0], EPS);
        assert!(r.point.iter().all(|x| x.is_finite()));
        assert!((weight_sum(&r) - 1.0).abs() < EPS);
        assert!(r.weights.iter().all(|w| *w >= 0.0));
    }

    #[test]
    fn collinear_triangle_returns_the_best_edge() {
        // Three collinear points on the line y = 2 (x = -3, 1, 5): nearest is (0,2,0)
        // on the edge between the first two.
        let pts = [[-3.0, 2.0, 0.0], [1.0, 2.0, 0.0], [5.0, 2.0, 0.0]];
        let r = closest_point_to_origin(&pts);
        assert_close(r.point, [0.0, 2.0, 0.0], EPS);
        assert!(r.weights.iter().all(|w| w.is_finite() && *w >= 0.0));
        assert!((weight_sum(&r) - 1.0).abs() < EPS);
        assert_close(blend(&pts, &r.weights), r.point, 1e-12);
        // The far vertex is not needed.
        assert_eq!(r.weights[2], 0.0);
    }

    #[test]
    fn all_three_vertices_identical() {
        let p = [0.0, 3.0, 4.0];
        let r = closest_point_to_origin(&[p, p, p]);
        assert_close(r.point, p, EPS);
        assert!((weight_sum(&r) - 1.0).abs() < EPS);
    }

    #[test]
    fn sliver_triangle_stays_finite_and_optimal() {
        // Needle-thin triangle (width 1e-12 of its length): fall back to edges.
        let pts = [[1.0, 1.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 1.0 + 1e-13, 0.0]];
        let r = closest_point_to_origin(&pts);
        assert!(r.point.iter().all(|x| x.is_finite()));
        assert!(r.weights.iter().all(|w| w.is_finite() && *w >= 0.0));
        assert!((r.point[1] - 1.0).abs() < 1e-9, "point {:?}", r.point);
        assert!(r.point[0].abs() < 1e-9, "point {:?}", r.point);
    }

    #[test]
    fn flat_tetrahedron_never_reports_enclosure() {
        // Four coplanar points (z = 0) whose hull contains the origin in-plane.
        let pts = [
            [1.0, 0.0, 0.0],
            [-1.0, 1.0, 0.0],
            [-1.0, -1.0, 0.0],
            [0.5, 0.5, 0.0],
        ];
        let r = closest_point_to_origin(&pts);
        assert!(!r.encloses_origin);
        assert!(
            len3(r.point) < 1e-12,
            "origin lies in the flat hull: {:?}",
            r.point
        );
        assert!((weight_sum(&r) - 1.0).abs() < EPS);
        assert_close(blend(&pts, &r.weights), r.point, 1e-12);
    }

    #[test]
    fn flat_tetrahedron_off_plane_returns_plane_distance() {
        let pts = [
            [1.0, 0.0, 2.0],
            [-1.0, 1.0, 2.0],
            [-1.0, -1.0, 2.0],
            [0.0, 0.0, 2.0],
        ];
        let r = closest_point_to_origin(&pts);
        assert!(!r.encloses_origin);
        assert_close(r.point, [0.0, 0.0, 2.0], 1e-12);
    }

    #[test]
    fn tetrahedron_with_duplicated_vertex_is_flat() {
        let a = [1.0, 1.0, 1.0];
        let r = closest_point_to_origin(&[a, a, [2.0, 1.0, 1.0], [1.0, 2.0, 1.0]]);
        assert!(!r.encloses_origin);
        assert!(r.point.iter().all(|x| x.is_finite()));
        assert_close(r.point, a, 1e-12);
    }

    #[test]
    fn empty_input_is_the_origin() {
        let r = closest_point_to_origin(&[]);
        assert_eq!(r.point, [0.0; 3]);
        assert_eq!(r.weights, [0.0; 4]);
    }

    // ── regular cases ────────────────────────────────────────────────────────

    #[test]
    fn single_vertex() {
        let r = closest_point_to_origin(&[[3.0, 4.0, 0.0]]);
        assert_eq!(r.point, [3.0, 4.0, 0.0]);
        assert_eq!(r.weights, [1.0, 0.0, 0.0, 0.0]);
    }

    #[test]
    fn segment_interior_vertex_and_origin_on_segment() {
        let r = closest_point_to_origin(&[[-1.0, 2.0, 0.0], [3.0, 2.0, 0.0]]);
        assert_close(r.point, [0.0, 2.0, 0.0], EPS);
        assert!((r.weights[0] - 0.75).abs() < EPS && (r.weights[1] - 0.25).abs() < EPS);

        let r = closest_point_to_origin(&[[2.0, 0.0, 0.0], [4.0, 0.0, 0.0]]);
        assert_eq!(r.point, [2.0, 0.0, 0.0]);
        assert_eq!(r.weights, [1.0, 0.0, 0.0, 0.0]);

        let r = closest_point_to_origin(&[[-4.0, 0.0, 0.0], [-2.0, 0.0, 0.0]]);
        assert_eq!(r.point, [-2.0, 0.0, 0.0]);
        assert_eq!(r.weights, [0.0, 1.0, 0.0, 0.0]);

        let r = closest_point_to_origin(&[[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]]);
        assert!(len3(r.point) < EPS);
        assert!((r.weights[0] - 0.5).abs() < EPS && (r.weights[1] - 0.5).abs() < EPS);
    }

    #[test]
    fn triangle_interior_projection() {
        let pts = [[-1.0, -1.0, 3.0], [2.0, -1.0, 3.0], [-1.0, 2.0, 3.0]];
        let r = closest_point_to_origin(&pts);
        assert_close(r.point, [0.0, 0.0, 3.0], EPS);
        assert!(r.weights[..3].iter().all(|w| *w > 0.0));
        assert!((weight_sum(&r) - 1.0).abs() < EPS);
        assert_close(blend(&pts, &r.weights), r.point, 1e-12);
    }

    #[test]
    fn triangle_vertex_and_edge_regions() {
        // Vertex region: nearest to `a`.
        let pts = [[1.0, 0.0, 0.0], [3.0, 0.0, 0.0], [1.0, 2.0, 0.0]];
        let r = closest_point_to_origin(&pts);
        assert_eq!(r.point, [1.0, 0.0, 0.0]);
        assert_eq!(r.weights, [1.0, 0.0, 0.0, 0.0]);
        // Edge region: the origin projects beyond edge a-b.
        let pts = [[-2.0, 1.0, 0.0], [2.0, 1.0, 0.0], [0.0, 5.0, 0.0]];
        let r = closest_point_to_origin(&pts);
        assert_close(r.point, [0.0, 1.0, 0.0], EPS);
        assert_eq!(r.weights[2], 0.0);
        assert!((r.weights[0] - 0.5).abs() < EPS);
    }

    #[test]
    fn triangle_containing_the_origin_in_plane() {
        let pts = [[-1.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 1.0, 0.0]];
        let r = closest_point_to_origin(&pts);
        assert!(len3(r.point) < EPS);
        assert!(r.weights[..3].iter().all(|w| *w > 0.0));
    }

    #[test]
    fn tetrahedron_containing_the_origin() {
        let pts = [
            [1.0, 1.0, 1.0],
            [-1.0, -1.0, 1.0],
            [-1.0, 1.0, -1.0],
            [1.0, -1.0, -1.0],
        ];
        let r = closest_point_to_origin(&pts);
        assert!(r.encloses_origin);
        assert_eq!(r.point, [0.0; 3]);
        assert!((weight_sum(&r) - 1.0).abs() < EPS);
        assert!(r.weights.iter().all(|w| *w >= 0.0));
        assert_close(blend(&pts, &r.weights), [0.0; 3], 1e-12);
    }

    #[test]
    fn tetrahedron_outside_picks_the_nearest_face() {
        // Regular-ish tetrahedron shifted along +x so the origin is outside.
        let pts = [
            [3.0, 1.0, 1.0],
            [1.0, -1.0, 1.0],
            [1.0, 1.0, -1.0],
            [3.0, -1.0, -1.0],
        ];
        let r = closest_point_to_origin(&pts);
        assert!(!r.encloses_origin);
        // The vertex (1,-1,1)/(1,1,-1) edge at x = 1 is the nearest feature of
        // the solid: the closest point lies on the segment between them.
        assert_close(r.point, [1.0, 0.0, 0.0], 1e-12);
        assert_close(blend(&pts, &r.weights), r.point, 1e-12);
        assert_eq!(r.weights[0], 0.0);
        assert_eq!(r.weights[3], 0.0);
    }

    // ── randomised optimality checks ─────────────────────────────────────────

    /// Deterministic xorshift generator in `[-1, 1)`.
    struct Rng(u64);

    impl Rng {
        fn uniform(&mut self) -> f64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            ((self.0 >> 11) as f64 / (1u64 << 53) as f64) * 2.0 - 1.0
        }

        fn point(&mut self, scale: f64, offset: [f64; 3]) -> [f64; 3] {
            [
                offset[0] + scale * self.uniform(),
                offset[1] + scale * self.uniform(),
                offset[2] + scale * self.uniform(),
            ]
        }
    }

    #[test]
    fn random_simplices_satisfy_the_optimality_conditions() {
        let mut rng = Rng(0x9E37_79B9_7F4A_7C15);
        for case in 0..4000 {
            let n = 1 + case % 4;
            let offset = rng.point(3.0, [0.0; 3]);
            let pts: Vec<[f64; 3]> = (0..n).map(|_| rng.point(2.0, offset)).collect();
            let r = closest_point_to_origin(&pts);

            // Valid barycentric description.
            assert!(
                r.weights.iter().all(|w| w.is_finite() && *w >= 0.0),
                "{r:?}"
            );
            assert!((weight_sum(&r) - 1.0).abs() < 1e-9, "{r:?}");
            assert!(r.weights[n..].iter().all(|w| *w == 0.0), "{r:?}");
            assert_close(blend(&pts, &r.weights), r.point, 1e-9);

            // KKT condition for the nearest point of a convex hull to the origin:
            // v . (p_i - v) >= 0 for every vertex.
            let v = r.point;
            let v_sq = dot3(v, v);
            for p in &pts {
                assert!(
                    dot3(v, *p) >= v_sq - 1e-9,
                    "vertex {p:?} lies on the origin side of the supporting plane at {v:?}"
                );
            }
            // No random convex combination is closer than the reported point.
            for _ in 0..8 {
                let mut weights = [0.0; 4];
                let mut total = 0.0;
                for w in weights.iter_mut().take(n) {
                    *w = rng.uniform().abs() + 1e-3;
                    total += *w;
                }
                for w in weights.iter_mut().take(n) {
                    *w /= total;
                }
                let q = blend(&pts, &weights);
                assert!(dot3(q, q) >= v_sq - 1e-9, "{q:?} beats {v:?}");
            }
            if r.encloses_origin {
                assert_eq!(n, 4);
                assert_eq!(v, [0.0; 3]);
            }
        }
    }
}
