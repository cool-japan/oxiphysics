// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Closest-point projection onto the convex hull of a vertex set.
//!
//! The vertices are first turned into a proper triangulated hull with
//! QuickHull; the projection is then the closest of the closest points on the
//! hull facets. Vertex sets that do not span three dimensions (coplanar,
//! collinear or coincident points) have a flat hull; they are handled by
//! dedicated 2-D / 1-D / 0-D paths.

use super::functions::{
    closest_point_on_segment, closest_point_on_triangle_local, dot, scale, sub,
};
use oxiphysics_core::math::Vec3;
use oxiphysics_geometry::quickhull::ConvexHull3DVec;

/// Relative tolerance (to the diagonal of the bounding box) below which a
/// vertex set is considered to lie on a line or a plane. It matches the
/// tolerance of the QuickHull construction, so a set the hull builder rejects
/// as flat is also flat here.
const FLATNESS_TOLERANCE: f64 = 1e-10;

/// Project a point onto the surface of the convex hull of `vertices`.
///
/// Returns the point of the hull boundary that is closest to `point`:
///
/// * a point outside the hull is mapped to the closest point of the hull
///   (a facet, edge or vertex point);
/// * a point inside the hull is mapped to the closest point of the boundary,
///   i.e. it is pushed out through the nearest facet.
///
/// The hull is built with QuickHull, so `vertices` may be given in any order
/// and may contain points that lie inside the hull (they are ignored), as well
/// as duplicates. Vertices with a non-finite coordinate are ignored.
///
/// Vertex sets that do not span three dimensions have a flat hull, which is
/// treated as the shape itself:
///
/// * coplanar points - the filled convex polygon (a point above or below the
///   polygon is projected onto its plane, a point beside it onto its boundary);
/// * collinear points - the segment between the two extreme points;
/// * one distinct point - that point.
///
/// An empty vertex set has no hull; `point` is returned unchanged.
///
/// # Examples
/// ```
/// use oxiphysics_collision::query::project_point_convex;
/// let cube = [
///     [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0],
///     [0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [0.0, 1.0, 1.0], [1.0, 1.0, 1.0],
/// ];
/// let p = project_point_convex([2.0, 0.5, 0.5], &cube);
/// assert!((p[0] - 1.0).abs() < 1e-12 && (p[1] - 0.5).abs() < 1e-12);
/// ```
pub fn project_point_convex(point: [f64; 3], vertices: &[[f64; 3]]) -> [f64; 3] {
    let cloud: Vec<[f64; 3]> = vertices
        .iter()
        .copied()
        .filter(|v| v.iter().all(|c| c.is_finite()))
        .collect();
    match cloud.as_slice() {
        [] => point,
        [only] => *only,
        _ => {
            let points: Vec<Vec3> = cloud.iter().map(|v| Vec3::new(v[0], v[1], v[2])).collect();
            match ConvexHull3DVec::build(&points) {
                Some(hull) => closest_on_hull_facets(point, &hull),
                None => project_onto_flat_hull(point, &cloud),
            }
        }
    }
}

fn to_array(v: &Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}

/// Closest point to `point` over all triangular facets of `hull`.
fn closest_on_hull_facets(point: [f64; 3], hull: &ConvexHull3DVec) -> [f64; 3] {
    let mut best = hull.vertices.first().map(to_array).unwrap_or(point);
    let mut best_dist_sq = f64::INFINITY;
    for face in &hull.faces {
        let (Some(a), Some(b), Some(c)) = (
            hull.vertices.get(face[0]),
            hull.vertices.get(face[1]),
            hull.vertices.get(face[2]),
        ) else {
            continue;
        };
        let candidate =
            closest_point_on_triangle_local(point, to_array(a), to_array(b), to_array(c));
        if !candidate.iter().all(|c| c.is_finite()) {
            continue;
        }
        let d = sub(point, candidate);
        let dist_sq = dot(d, d);
        if dist_sq < best_dist_sq {
            best_dist_sq = dist_sq;
            best = candidate;
        }
    }
    best
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn length(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

/// Index of the point of `cloud` farthest from `from`, with its distance.
fn farthest_from(cloud: &[[f64; 3]], from: [f64; 3]) -> (usize, f64) {
    let mut best_index = 0;
    let mut best_dist = 0.0;
    for (i, v) in cloud.iter().enumerate() {
        let d = length(sub(*v, from));
        if d > best_dist {
            best_dist = d;
            best_index = i;
        }
    }
    (best_index, best_dist)
}

/// Projection onto the hull of a vertex set that does not span three
/// dimensions. `cloud` holds at least two finite points.
fn project_onto_flat_hull(point: [f64; 3], cloud: &[[f64; 3]]) -> [f64; 3] {
    let Some(&origin) = cloud.first() else {
        return point;
    };
    let mut lowest = origin;
    let mut highest = origin;
    for v in cloud {
        for k in 0..3 {
            lowest[k] = lowest[k].min(v[k]);
            highest[k] = highest[k].max(v[k]);
        }
    }
    let tolerance = FLATNESS_TOLERANCE * length(sub(highest, lowest));

    let (far_index, far_dist) = farthest_from(cloud, origin);
    if far_dist <= tolerance {
        return origin;
    }
    let far = cloud[far_index];
    let line_dir = sub(far, origin);
    let line_len = length(line_dir);
    let unit_dir = scale(line_dir, 1.0 / line_len);

    // The point farthest from the line through `origin` and `far` decides
    // between a segment and a polygon.
    let mut third_index = 0;
    let mut third_dist = 0.0;
    for (i, v) in cloud.iter().enumerate() {
        let rel = sub(*v, origin);
        let along = dot(rel, unit_dir);
        let perp = sub(rel, scale(unit_dir, along));
        let d = length(perp);
        if d > third_dist {
            third_dist = d;
            third_index = i;
        }
    }
    if third_dist <= tolerance {
        return project_onto_segment_hull(point, cloud, origin, unit_dir);
    }

    // Coplanar points: express them in an orthonormal basis of the plane and
    // keep the corners of their 2-D hull.
    let normal = cross(line_dir, sub(cloud[third_index], origin));
    let in_plane = cross(normal, unit_dir);
    let in_plane_len = length(in_plane);
    let basis_v = scale(in_plane, 1.0 / in_plane_len);
    let coords: Vec<[f64; 2]> = cloud
        .iter()
        .map(|v| {
            let rel = sub(*v, origin);
            [dot(rel, unit_dir), dot(rel, basis_v)]
        })
        .collect();
    let corners = convex_hull_2d(&coords);
    if corners.len() < 3 {
        return project_onto_segment_hull(point, cloud, origin, unit_dir);
    }

    // Fan triangulation of the convex polygon; the closest of the closest
    // points on its triangles is the closest point of the polygon.
    let mut best = cloud[corners[0]];
    let mut best_dist_sq = f64::INFINITY;
    for pair in corners[1..].windows(2) {
        let candidate = closest_point_on_triangle_local(
            point,
            cloud[corners[0]],
            cloud[pair[0]],
            cloud[pair[1]],
        );
        if !candidate.iter().all(|c| c.is_finite()) {
            continue;
        }
        let d = sub(point, candidate);
        let dist_sq = dot(d, d);
        if dist_sq < best_dist_sq {
            best_dist_sq = dist_sq;
            best = candidate;
        }
    }
    best
}

/// Projection onto the segment between the two extreme points of collinear
/// points lying along `unit_dir` through `origin`.
fn project_onto_segment_hull(
    point: [f64; 3],
    cloud: &[[f64; 3]],
    origin: [f64; 3],
    unit_dir: [f64; 3],
) -> [f64; 3] {
    let mut low = origin;
    let mut high = origin;
    let mut low_t = 0.0;
    let mut high_t = 0.0;
    for v in cloud {
        let t = dot(sub(*v, origin), unit_dir);
        if t < low_t {
            low_t = t;
            low = *v;
        }
        if t > high_t {
            high_t = t;
            high = *v;
        }
    }
    let (closest, _) = closest_point_on_segment(point, low, high);
    closest
}

/// Indices of the corners of the convex hull of `coords`, counter-clockwise
/// (Andrew's monotone chain). Fewer than three indices are returned when the
/// points do not span two dimensions.
fn convex_hull_2d(coords: &[[f64; 2]]) -> Vec<usize> {
    let mut order: Vec<usize> = (0..coords.len()).collect();
    order.sort_by(|&a, &b| {
        coords[a][0]
            .total_cmp(&coords[b][0])
            .then(coords[a][1].total_cmp(&coords[b][1]))
    });
    order.dedup_by(|a, b| coords[*a] == coords[*b]);
    if order.len() < 3 {
        return order;
    }
    let turn = |o: usize, a: usize, b: usize| {
        (coords[a][0] - coords[o][0]) * (coords[b][1] - coords[o][1])
            - (coords[a][1] - coords[o][1]) * (coords[b][0] - coords[o][0])
    };
    let mut lower = half_hull(order.iter().copied(), &turn);
    let mut upper = half_hull(order.iter().rev().copied(), &turn);
    lower.pop();
    upper.pop();
    lower.extend(upper);
    lower
}

/// One chain of the monotone-chain hull: the points of `indices` (sorted by
/// increasing or decreasing coordinate) that make only left turns, where
/// `turn(o, a, b)` is the cross product of `o -> a` and `o -> b`.
fn half_hull<I>(indices: I, turn: &impl Fn(usize, usize, usize) -> f64) -> Vec<usize>
where
    I: Iterator<Item = usize>,
{
    let mut hull: Vec<usize> = Vec::new();
    for i in indices {
        while hull.len() >= 2 && turn(hull[hull.len() - 2], hull[hull.len() - 1], i) <= 0.0 {
            hull.pop();
        }
        hull.push(i);
    }
    hull
}

#[cfg(test)]
mod tests {
    use super::super::functions::add;
    use super::*;

    fn centroid(points: &[[f64; 3]]) -> [f64; 3] {
        let mut sum = [0.0; 3];
        for p in points {
            sum = add(sum, *p);
        }
        let n = points.len() as f64;
        [sum[0] / n, sum[1] / n, sum[2] / n]
    }

    fn unit_cube_binary_order() -> Vec<[f64; 3]> {
        let mut corners = Vec::new();
        for i in 0..8u32 {
            corners.push([
                f64::from(i & 1),
                f64::from((i >> 1) & 1),
                f64::from((i >> 2) & 1),
            ]);
        }
        corners
    }

    fn assert_close(actual: [f64; 3], expected: [f64; 3], what: &str) {
        let d = length(sub(actual, expected));
        assert!(
            d < 1e-9,
            "{what}: got {actual:?}, expected {expected:?} (distance {d})"
        );
    }

    /// Deterministic generator in `[0, 1)`.
    struct Lcg(u64);

    impl Lcg {
        fn next_unit(&mut self) -> f64 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (self.0 >> 11) as f64 / (1u64 << 53) as f64
        }

        fn next_below(&mut self, bound: usize) -> usize {
            ((self.next_unit() * bound as f64) as usize).min(bound - 1)
        }
    }

    #[test]
    fn cube_face_vertex_and_edge_regions_in_binary_corner_order() {
        let cube = unit_cube_binary_order();
        assert_close(
            project_point_convex([2.0, 0.5, 0.5], &cube),
            [1.0, 0.5, 0.5],
            "face region",
        );
        assert_close(
            project_point_convex([2.0, 2.0, 2.0], &cube),
            [1.0, 1.0, 1.0],
            "vertex region",
        );
        assert_close(
            project_point_convex([2.0, 2.0, 0.5], &cube),
            [1.0, 1.0, 0.5],
            "edge region",
        );
        assert_close(
            project_point_convex([-3.0, 0.25, 4.0], &cube),
            [0.0, 0.25, 1.0],
            "edge/face region on the opposite side",
        );
    }

    #[test]
    fn interior_points_are_pushed_out_through_the_nearest_face() {
        let cube = unit_cube_binary_order();
        assert_close(
            project_point_convex([0.5, 0.5, 0.9], &cube),
            [0.5, 0.5, 1.0],
            "top face",
        );
        assert_close(
            project_point_convex([0.2, 0.5, 0.5], &cube),
            [0.0, 0.5, 0.5],
            "x = 0 face",
        );
        assert_close(
            project_point_convex([0.5, 0.15, 0.45], &cube),
            [0.5, 0.0, 0.45],
            "y = 0 face",
        );
    }

    #[test]
    fn points_on_the_surface_are_fixed() {
        let cube = unit_cube_binary_order();
        for p in [[1.0, 0.3, 0.6], [0.0, 0.0, 0.0], [0.5, 1.0, 0.5]] {
            assert_close(project_point_convex(p, &cube), p, "surface point");
        }
    }

    #[test]
    fn vertex_order_and_interior_vertices_do_not_matter() {
        let cube = unit_cube_binary_order();
        let queries = [
            [2.0, 0.5, 0.5],
            [0.5, 0.5, 0.9],
            [-1.0, -1.0, 0.5],
            [0.3, 1.7, 1.2],
        ];
        let expected: Vec<[f64; 3]> = queries
            .iter()
            .map(|q| project_point_convex(*q, &cube))
            .collect();
        let mut rng = Lcg(0x0123_4567_89AB_CDEF);
        for _ in 0..20 {
            let mut shuffled = cube.clone();
            for i in (1..shuffled.len()).rev() {
                let j = rng.next_below(i + 1);
                shuffled.swap(i, j);
            }
            // Interior points and duplicates of corners.
            shuffled.push([0.5, 0.5, 0.5]);
            shuffled.push([0.25, 0.6, 0.4]);
            shuffled.push(cube[3]);
            for (q, want) in queries.iter().zip(expected.iter()) {
                assert_close(project_point_convex(*q, &shuffled), *want, "shuffled cube");
            }
        }
    }

    #[test]
    fn matches_the_analytic_projection_onto_a_random_axis_aligned_box() {
        let mut rng = Lcg(0xDEAD_BEEF_1234_5678);
        for case in 0..60 {
            let lo = [
                -2.0 * rng.next_unit(),
                -rng.next_unit(),
                -3.0 * rng.next_unit(),
            ];
            let hi = [
                0.5 + 2.0 * rng.next_unit(),
                0.5 + rng.next_unit(),
                0.5 + 3.0 * rng.next_unit(),
            ];
            let mut corners = Vec::new();
            for i in 0..8usize {
                corners.push([
                    if i & 1 == 0 { lo[0] } else { hi[0] },
                    if i & 2 == 0 { lo[1] } else { hi[1] },
                    if i & 4 == 0 { lo[2] } else { hi[2] },
                ]);
            }
            let q = [
                lo[0] - 1.5 + (hi[0] - lo[0] + 3.0) * rng.next_unit(),
                lo[1] - 1.5 + (hi[1] - lo[1] + 3.0) * rng.next_unit(),
                lo[2] - 1.5 + (hi[2] - lo[2] + 3.0) * rng.next_unit(),
            ];
            let inside = (0..3).all(|k| q[k] >= lo[k] && q[k] <= hi[k]);
            let mut expected = [0.0; 3];
            if inside {
                // Push out through the face with the smallest gap.
                let mut best_gap = f64::INFINITY;
                expected = q;
                for k in 0..3 {
                    for (bound, gap) in [(lo[k], q[k] - lo[k]), (hi[k], hi[k] - q[k])] {
                        if gap < best_gap {
                            best_gap = gap;
                            expected = q;
                            expected[k] = bound;
                        }
                    }
                }
            } else {
                for k in 0..3 {
                    expected[k] = q[k].clamp(lo[k], hi[k]);
                }
            }
            let got = project_point_convex(q, &corners);
            let d = length(sub(got, expected));
            assert!(
                d < 1e-9,
                "case {case}: query {q:?} box {lo:?}..{hi:?}: got {got:?}, expected {expected:?}"
            );
        }
    }

    /// Distance from `point` to the hull of `cloud` as the minimum over the
    /// triangles of all vertex triples (every such triangle lies inside the
    /// hull and every facet is a union of them, so this is exact for points
    /// outside the hull).
    fn brute_force_outside_distance(point: [f64; 3], cloud: &[[f64; 3]]) -> f64 {
        let mut best = f64::INFINITY;
        for i in 0..cloud.len() {
            for j in (i + 1)..cloud.len() {
                for k in (j + 1)..cloud.len() {
                    let cp = closest_point_on_triangle_local(point, cloud[i], cloud[j], cloud[k]);
                    best = best.min(length(sub(point, cp)));
                }
            }
        }
        best
    }

    #[test]
    fn random_clouds_satisfy_the_optimality_conditions_of_the_projection() {
        let mut rng = Lcg(0x0BAD_C0DE_5EED_0001);
        for case in 0..40 {
            let n = 5 + rng.next_below(10);
            let cloud: Vec<[f64; 3]> = (0..n)
                .map(|_| {
                    [
                        4.0 * rng.next_unit() - 2.0,
                        3.0 * rng.next_unit() - 1.0,
                        5.0 * rng.next_unit() - 2.5,
                    ]
                })
                .collect();
            let hull_center = centroid(&cloud);
            // A query far outside the cloud, in a random direction.
            let dir = [
                rng.next_unit() - 0.5,
                rng.next_unit() - 0.5,
                rng.next_unit() - 0.5,
            ];
            let dir_len = length(dir);
            if dir_len < 1e-3 {
                continue;
            }
            let q = add(
                hull_center,
                [
                    dir[0] / dir_len * 7.0,
                    dir[1] / dir_len * 7.0,
                    dir[2] / dir_len * 7.0,
                ],
            );
            let p = project_point_convex(q, &cloud);
            // Variational inequality: (q - p) . (v - p) <= 0 for every vertex v.
            let toward_query = sub(q, p);
            for v in &cloud {
                let slack = dot(toward_query, sub(*v, p));
                assert!(
                    slack <= 1e-9,
                    "case {case}: vertex {v:?} lies beyond the supporting plane at {p:?} (slack {slack})"
                );
            }
            let exact = brute_force_outside_distance(q, &cloud);
            let got = length(sub(q, p));
            assert!(
                (got - exact).abs() < 1e-9,
                "case {case}: distance {got} differs from brute force {exact}"
            );
        }
    }

    #[test]
    fn coplanar_vertices_project_onto_the_filled_polygon() {
        // The square from the diamond test: (+-1, 0, 0), (0, 0, +-1).
        let diamond = [
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
            [-1.0, 0.0, 0.0],
            [0.0, 0.0, -1.0],
        ];
        assert_close(
            project_point_convex([3.0, 0.0, 0.0], &diamond),
            [1.0, 0.0, 0.0],
            "outside, vertex region",
        );
        assert_close(
            project_point_convex([0.2, 0.7, 0.1], &diamond),
            [0.2, 0.0, 0.1],
            "above the polygon",
        );
        assert_close(
            project_point_convex([1.0, 0.5, 1.0], &diamond),
            [0.5, 0.0, 0.5],
            "beside the polygon, edge region",
        );
        // Interior corner points of a larger point set do not change the polygon.
        let mut with_interior = diamond.to_vec();
        with_interior.push([0.1, 0.0, 0.1]);
        with_interior.push([0.0, 0.0, 0.0]);
        assert_close(
            project_point_convex([0.2, -0.4, 0.1], &with_interior),
            [0.2, 0.0, 0.1],
            "below the polygon",
        );
    }

    #[test]
    fn collinear_vertices_project_onto_the_extreme_segment() {
        let line = [
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
        ];
        assert_close(
            project_point_convex([5.0, 1.0, 0.0], &line),
            [2.0, 0.0, 0.0],
            "beyond the end",
        );
        assert_close(
            project_point_convex([-4.0, 1.0, 3.0], &line),
            [-1.0, 0.0, 0.0],
            "beyond the start",
        );
        assert_close(
            project_point_convex([0.7, 2.0, -1.0], &line),
            [0.7, 0.0, 0.0],
            "beside the segment",
        );
    }

    #[test]
    fn tiny_clouds_and_degenerate_inputs() {
        assert_eq!(project_point_convex([1.0, 2.0, 3.0], &[]), [1.0, 2.0, 3.0]);
        assert_eq!(
            project_point_convex([5.0, 5.0, 5.0], &[[1.0, 2.0, 3.0]]),
            [1.0, 2.0, 3.0]
        );
        // Two distinct points form a segment.
        let segment = [[0.0, 0.0, 0.0], [2.0, 2.0, 0.0]];
        assert_close(
            project_point_convex([2.0, 0.0, 0.0], &segment),
            [1.0, 1.0, 0.0],
            "segment",
        );
        // Many copies of one point are one point.
        let same = [[1.0, -1.0, 0.5]; 6];
        assert_close(
            project_point_convex([9.0, 9.0, 9.0], &same),
            [1.0, -1.0, 0.5],
            "coincident",
        );
        // Non-finite vertices are ignored.
        let mut with_nan = unit_cube_binary_order();
        with_nan.push([f64::NAN, 0.0, 0.0]);
        with_nan.push([0.0, f64::INFINITY, 0.0]);
        assert_close(
            project_point_convex([2.0, 0.5, 0.5], &with_nan),
            [1.0, 0.5, 0.5],
            "non-finite vertices",
        );
    }

    #[test]
    fn tetrahedron_projection_matches_the_closest_facet() {
        let tet = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        // Beyond the slanted face x + y + z = 1: the foot point is the facet
        // centroid direction.
        let q = [1.0, 1.0, 1.0];
        let third = 1.0 / 3.0;
        assert_close(
            project_point_convex(q, &tet),
            [third, third, third],
            "slanted face",
        );
        // The point (0.25, 0.25, 0.25) is 0.25 from the three coordinate faces
        // but only 0.25 / sqrt(3) from the slanted face, whose foot point is
        // the face centroid.
        assert_close(
            project_point_convex([0.25, 0.25, 0.25], &tet),
            [third, third, third],
            "interior point leaves through the slanted face",
        );
    }
}
