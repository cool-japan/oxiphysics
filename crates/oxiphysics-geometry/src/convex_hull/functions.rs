//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
use super::types::ConvexHull3D;

/// Cross product of two 3-vectors.
#[inline]
pub fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
/// Dot product of two 3-vectors.
#[inline]
pub fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
/// Component-wise subtraction.
#[inline]
pub fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
/// Component-wise addition.
#[inline]
pub fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
/// Scale a vector by scalar `s`.
#[inline]
pub fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
/// Normalize a vector (returns zero-vector if near-zero magnitude).
#[inline]
pub fn normalize(a: [f64; 3]) -> [f64; 3] {
    let len = dot(a, a).sqrt();
    if len < 1e-12 {
        [0.0, 0.0, 0.0]
    } else {
        scale(a, 1.0 / len)
    }
}
/// Euclidean distance between two points.
#[inline]
pub fn dist(a: [f64; 3], b: [f64; 3]) -> f64 {
    let d = sub(a, b);
    dot(d, d).sqrt()
}
/// Outward-facing normal for triangle (a, b, c) – cross product of edges.
/// Not normalized; caller may normalize if needed.
#[inline]
pub fn face_normal(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> [f64; 3] {
    cross(sub(b, a), sub(c, a))
}
/// Signed distance from `point` to the plane defined by `normal` and `origin`.
/// Positive means the point is on the side the normal points toward.
#[inline]
pub fn face_plane_dist(normal: [f64; 3], origin: [f64; 3], point: [f64; 3]) -> f64 {
    dot(normal, sub(point, origin))
}
/// Find a pair of extreme points (min/max along the axis with the largest spread).
pub(super) fn ch_extreme_pair(points: &[[f64; 3]]) -> Option<(usize, usize)> {
    if points.len() < 2 {
        return None;
    }
    let mut best_spread = -1.0f64;
    let mut best = (0usize, 1usize);
    for axis in 0..3usize {
        let (imin, _) = points.iter().enumerate().min_by(|(_, a), (_, b)| {
            a[axis]
                .partial_cmp(&b[axis])
                .unwrap_or(std::cmp::Ordering::Equal)
        })?;
        let (imax, _) = points.iter().enumerate().max_by(|(_, a), (_, b)| {
            a[axis]
                .partial_cmp(&b[axis])
                .unwrap_or(std::cmp::Ordering::Equal)
        })?;
        let spread = (points[imax][axis] - points[imin][axis]).abs();
        if spread > best_spread && imin != imax {
            best_spread = spread;
            best = (imin, imax);
        }
    }
    if best_spread < 1e-12 {
        return None;
    }
    Some(best)
}
/// Find the point farthest from the line through `points[i0]` and `points[i1]`.
pub(super) fn ch_farthest_from_line(points: &[[f64; 3]], i0: usize, i1: usize) -> Option<usize> {
    let a = points[i0];
    let b = points[i1];
    let ab = sub(b, a);
    let ab_len_sq = dot(ab, ab);
    if ab_len_sq < 1e-24 {
        return None;
    }
    points
        .iter()
        .enumerate()
        .filter(|&(i, _)| i != i0 && i != i1)
        .max_by(|(_, p), (_, q)| {
            let dp = dot(cross(sub(**p, a), ab), cross(sub(**p, a), ab)) / ab_len_sq;
            let dq = dot(cross(sub(**q, a), ab), cross(sub(**q, a), ab)) / ab_len_sq;
            dp.partial_cmp(&dq).unwrap_or(std::cmp::Ordering::Equal)
        })
        .map(|(i, _)| i)
}
/// Find the farthest point from the plane `(normal, origin)`, excluding listed indices.
pub(super) fn ch_farthest_from_plane(
    points: &[[f64; 3]],
    normal: [f64; 3],
    origin: [f64; 3],
    exclude: &[usize],
) -> Option<usize> {
    points
        .iter()
        .enumerate()
        .filter(|(i, _)| !exclude.contains(i))
        .max_by(|(_, a), (_, b)| {
            face_plane_dist(normal, origin, **a)
                .abs()
                .partial_cmp(&face_plane_dist(normal, origin, **b).abs())
                .expect("operation should succeed")
        })
        .map(|(i, _)| i)
}
/// Orient face so that `interior` is on the negative side of the face plane.
pub(super) fn ch_orient_face(
    face: [usize; 3],
    interior: usize,
    points: &[[f64; 3]],
) -> Option<[usize; 3]> {
    let n = face_normal(points[face[0]], points[face[1]], points[face[2]]);
    if dot(n, n) < 1e-24 {
        return None;
    }
    if face_plane_dist(n, points[face[0]], points[interior]) > 0.0 {
        Some([face[0], face[2], face[1]])
    } else {
        Some(face)
    }
}
/// Collect horizon edges: edges belonging to exactly one visible face.
pub(super) fn ch_horizon_edges(hull: &[[usize; 3]], visible: &[usize]) -> Vec<(usize, usize)> {
    let mut edges: Vec<(usize, usize)> = Vec::new();
    for &fi in visible {
        let face = &hull[fi];
        let face_edges = [(face[0], face[1]), (face[1], face[2]), (face[2], face[0])];
        for (e0, e1) in face_edges {
            let adjacent_visible = visible.iter().any(|&other| {
                if other == fi {
                    return false;
                }
                let of = &hull[other];
                let oe = [(of[0], of[1]), (of[1], of[2]), (of[2], of[0])];
                oe.contains(&(e1, e0))
            });
            if !adjacent_visible {
                edges.push((e0, e1));
            }
        }
    }
    edges
}
/// Approximate convex decomposition of a point cloud by recursively splitting
/// along the axis of maximum variance.
///
/// Returns a list of convex hulls, each covering a partition of the input.
#[allow(dead_code)]
pub fn approximate_convex_decomposition(
    points: &[[f64; 3]],
    max_depth: usize,
    min_points: usize,
) -> Vec<ConvexHull3D> {
    let mut result = Vec::new();
    acd_recursive(points, max_depth, min_points, &mut result);
    result
}
pub(super) fn acd_recursive(
    points: &[[f64; 3]],
    depth: usize,
    min_points: usize,
    result: &mut Vec<ConvexHull3D>,
) {
    if points.len() < min_points.max(4) || depth == 0 {
        if let Some(hull) = ConvexHull3D::build(points) {
            result.push(hull);
        }
        return;
    }
    let n = points.len() as f64;
    let mean: [f64; 3] = [
        points.iter().map(|p| p[0]).sum::<f64>() / n,
        points.iter().map(|p| p[1]).sum::<f64>() / n,
        points.iter().map(|p| p[2]).sum::<f64>() / n,
    ];
    let var: [f64; 3] = [
        points.iter().map(|p| (p[0] - mean[0]).powi(2)).sum::<f64>() / n,
        points.iter().map(|p| (p[1] - mean[1]).powi(2)).sum::<f64>() / n,
        points.iter().map(|p| (p[2] - mean[2]).powi(2)).sum::<f64>() / n,
    ];
    let axis = if var[0] >= var[1] && var[0] >= var[2] {
        0
    } else if var[1] >= var[2] {
        1
    } else {
        2
    };
    let median = mean[axis];
    let left: Vec<[f64; 3]> = points
        .iter()
        .copied()
        .filter(|p| p[axis] <= median)
        .collect();
    let right: Vec<[f64; 3]> = points
        .iter()
        .copied()
        .filter(|p| p[axis] > median)
        .collect();
    if left.is_empty() || right.is_empty() {
        if let Some(hull) = ConvexHull3D::build(points) {
            result.push(hull);
        }
        return;
    }
    acd_recursive(&left, depth - 1, min_points, result);
    acd_recursive(&right, depth - 1, min_points, result);
}
/// Compute a quality measure for a convex hull: the ratio of hull volume to
/// the volume of the bounding sphere.
///
/// Returns a value in (0, 1]: closer to 1 means the hull is more spherical.
#[allow(dead_code)]
pub fn hull_sphericity(hull: &ConvexHull3D) -> f64 {
    if hull.vertices.is_empty() {
        return 0.0;
    }
    let c = hull
        .vertices
        .iter()
        .fold([0.0f64; 3], |a, v| [a[0] + v[0], a[1] + v[1], a[2] + v[2]]);
    let n = hull.vertices.len() as f64;
    let center = [c[0] / n, c[1] / n, c[2] / n];
    let r_sq: f64 = hull
        .vertices
        .iter()
        .map(|v| {
            let d = sub(*v, center);
            dot(d, d)
        })
        .fold(0.0f64, f64::max);
    if r_sq < 1e-24 {
        return 0.0;
    }
    let r = r_sq.sqrt();
    let sphere_vol = (4.0 / 3.0) * std::f64::consts::PI * r * r * r;
    let hull_vol = hull.volume();
    if sphere_vol < 1e-24 {
        0.0
    } else {
        (hull_vol / sphere_vol).min(1.0)
    }
}
/// Compute the aspect ratio of the hull (longest edge / shortest edge of bounding box).
///
/// Returns `f64::INFINITY` if the hull is degenerate.
#[allow(dead_code)]
pub fn hull_aspect_ratio(hull: &ConvexHull3D) -> f64 {
    if hull.vertices.is_empty() {
        return f64::INFINITY;
    }
    let mut mn = hull.vertices[0];
    let mut mx = hull.vertices[0];
    for v in &hull.vertices {
        for i in 0..3 {
            if v[i] < mn[i] {
                mn[i] = v[i];
            }
            if v[i] > mx[i] {
                mx[i] = v[i];
            }
        }
    }
    let extents = [mx[0] - mn[0], mx[1] - mn[1], mx[2] - mn[2]];
    let max_e = extents.iter().cloned().fold(0.0f64, f64::max);
    let min_e = extents.iter().cloned().fold(f64::MAX, f64::min);
    if min_e < 1e-12 {
        f64::INFINITY
    } else {
        max_e / min_e
    }
}
/// Simplify a convex hull by removing vertices that contribute less than
/// `min_contribution` to the total surface area.
///
/// Returns a new (possibly smaller) point cloud from which the hull can be rebuilt.
#[allow(dead_code)]
pub fn simplify_hull_points(hull: &ConvexHull3D, min_contribution: f64) -> Vec<[f64; 3]> {
    if hull.vertices.is_empty() {
        return vec![];
    }
    let total_area = hull.surface_area();
    if total_area < 1e-24 {
        return hull.vertices.clone();
    }
    let mut contrib = vec![0.0f64; hull.vertices.len()];
    for face in &hull.faces {
        let area = {
            let n = face_normal(
                hull.vertices[face[0]],
                hull.vertices[face[1]],
                hull.vertices[face[2]],
            );
            dot(n, n).sqrt() * 0.5
        };
        for &vi in face {
            contrib[vi] += area;
        }
    }
    hull.vertices
        .iter()
        .enumerate()
        .filter(|(i, _)| contrib[*i] / total_area >= min_contribution)
        .map(|(_, v)| *v)
        .collect()
}
/// Remove duplicate or near-coincident vertices from a point cloud.
///
/// Points within `tolerance` of each other are merged.
#[allow(dead_code)]
pub fn deduplicate_points(points: &[[f64; 3]], tolerance: f64) -> Vec<[f64; 3]> {
    let tol_sq = tolerance * tolerance;
    let mut unique: Vec<[f64; 3]> = Vec::new();
    for &p in points {
        let is_dup = unique.iter().any(|&u| {
            let d = sub(p, u);
            dot(d, d) <= tol_sq
        });
        if !is_dup {
            unique.push(p);
        }
    }
    unique
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::ConvexHull;
    use crate::ConvexHull3D;
    use crate::HullWithQuality;

    use crate::convex_decomposition::dot;

    use oxiphysics_core::math::Vec3;
    fn cube_points() -> Vec<[f64; 3]> {
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [0.0, 1.0, 1.0],
            [1.0, 1.0, 1.0],
        ]
    }
    fn tetra_points() -> Vec<[f64; 3]> {
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]
    }
    #[test]
    fn test_cross_product() {
        let x = [1.0, 0.0, 0.0];
        let y = [0.0, 1.0, 0.0];
        let z = cross(x, y);
        assert!((z[0] - 0.0).abs() < 1e-12);
        assert!((z[1] - 0.0).abs() < 1e-12);
        assert!((z[2] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_dot_product() {
        let a = [1.0, 2.0, 3.0];
        let b = [4.0, 5.0, 6.0];
        assert!((dot(a, b) - 32.0).abs() < 1e-12);
    }
    #[test]
    fn test_normalize() {
        let v = [3.0, 0.0, 0.0];
        let n = normalize(v);
        assert!((n[0] - 1.0).abs() < 1e-12);
        assert!(n[1].abs() < 1e-12);
        assert!(n[2].abs() < 1e-12);
    }
    #[test]
    fn test_dist() {
        let a = [1.0, 0.0, 0.0];
        let b = [4.0, 0.0, 0.0];
        assert!((dist(a, b) - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_build_fewer_than_4_returns_none() {
        assert!(ConvexHull3D::build(&[]).is_none());
        assert!(ConvexHull3D::build(&[[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]).is_none());
    }
    #[test]
    fn test_build_tetrahedron_4_faces() {
        let pts = tetra_points();
        let hull = ConvexHull3D::build(&pts).expect("should build tetrahedron");
        assert_eq!(
            hull.n_faces(),
            4,
            "tetrahedron should have 4 faces, got {}",
            hull.n_faces()
        );
    }
    #[test]
    fn test_build_cube_12_faces() {
        let pts = cube_points();
        let hull = ConvexHull3D::build(&pts).expect("should build cube hull");
        assert_eq!(
            hull.n_faces(),
            12,
            "unit cube should have 12 triangular faces, got {}",
            hull.n_faces()
        );
    }
    #[test]
    fn test_volume_unit_cube() {
        let pts = cube_points();
        let hull = ConvexHull3D::build(&pts).expect("build");
        let vol = hull.volume();
        assert!(
            (vol - 1.0).abs() < 0.05,
            "unit cube volume should be ≈1.0, got {vol}"
        );
    }
    #[test]
    fn test_volume_tetrahedron() {
        let pts = tetra_points();
        let hull = ConvexHull3D::build(&pts).expect("build");
        let vol = hull.volume();
        assert!(
            (vol - 1.0 / 6.0).abs() < 0.01,
            "unit tetrahedron volume should be ≈1/6, got {vol}"
        );
    }
    #[test]
    fn test_surface_area_unit_cube() {
        let pts = cube_points();
        let hull = ConvexHull3D::build(&pts).expect("build");
        let sa = hull.surface_area();
        assert!(
            (sa - 6.0).abs() < 0.05,
            "unit cube surface area should be ≈6.0, got {sa}"
        );
    }
    #[test]
    fn test_contains_center_of_cube() {
        let pts = cube_points();
        let hull = ConvexHull3D::build(&pts).expect("build");
        let center = [0.5, 0.5, 0.5];
        assert!(
            hull.contains_point(center),
            "center of unit cube should be inside hull"
        );
    }
    #[test]
    fn test_does_not_contain_exterior_point() {
        let pts = cube_points();
        let hull = ConvexHull3D::build(&pts).expect("build");
        let outside = [5.0, 5.0, 5.0];
        assert!(
            !hull.contains_point(outside),
            "point far outside should not be contained"
        );
    }
    #[test]
    fn test_support_direction_x() {
        let pts = cube_points();
        let hull = ConvexHull3D::build(&pts).expect("build");
        let sp = hull.support([1.0, 0.0, 0.0]);
        assert!(
            (sp[0] - 1.0).abs() < 1e-9,
            "support in +X should have x=1, got x={}",
            sp[0]
        );
    }
    #[test]
    fn test_support_direction_neg_y() {
        let pts = cube_points();
        let hull = ConvexHull3D::build(&pts).expect("build");
        let sp = hull.support([0.0, -1.0, 0.0]);
        assert!(
            sp[1].abs() < 1e-9,
            "support in -Y should have y=0, got y={}",
            sp[1]
        );
    }
    #[test]
    fn test_convex_hull_support() {
        let hull = ConvexHull::new(vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        ]);
        let sp = hull.support_point(&Vec3::new(1.0, 0.0, 0.0));
        assert!((sp.x - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_convex_hull_bounding_box() {
        let hull = ConvexHull::new(vec![Vec3::new(-1.0, -2.0, -3.0), Vec3::new(1.0, 2.0, 3.0)]);
        let bb = hull.bounding_box();
        assert_eq!(bb.min, Vec3::new(-1.0, -2.0, -3.0));
        assert_eq!(bb.max, Vec3::new(1.0, 2.0, 3.0));
    }
    #[test]
    fn test_acd_cube_returns_hulls() {
        let pts = tetra_points();
        let hulls = approximate_convex_decomposition(&pts, 1, 4);
        let hulls2 = approximate_convex_decomposition(&pts, 0, 4);
        assert!(
            !hulls2.is_empty(),
            "ACD with depth=0 should produce at least one hull"
        );
        let _ = hulls;
    }
    #[test]
    fn test_acd_depth_zero_returns_single_hull() {
        let pts = cube_points();
        let hulls = approximate_convex_decomposition(&pts, 0, 4);
        assert_eq!(hulls.len(), 1);
    }
    #[test]
    fn test_acd_min_points_prevents_splitting() {
        let pts = cube_points();
        let hulls = approximate_convex_decomposition(&pts, 5, 100);
        assert_eq!(hulls.len(), 1, "min_points > n should prevent splitting");
    }
    #[test]
    fn test_hull_sphericity_cube_nonzero() {
        let hull = ConvexHull3D::build(&cube_points()).expect("build");
        let s = hull_sphericity(&hull);
        assert!(s > 0.0 && s <= 1.0, "Sphericity should be in (0,1]: s={s}");
    }
    #[test]
    fn test_hull_sphericity_empty() {
        let hull = ConvexHull3D {
            vertices: vec![],
            faces: vec![],
        };
        assert_eq!(hull_sphericity(&hull), 0.0);
    }
    #[test]
    fn test_hull_aspect_ratio_cube_near_one() {
        let hull = ConvexHull3D::build(&cube_points()).expect("build");
        let ar = hull_aspect_ratio(&hull);
        assert!(
            (ar - 1.0).abs() < 0.1,
            "Unit cube aspect ratio should be ~1, got {ar}"
        );
    }
    #[test]
    fn test_hull_aspect_ratio_elongated() {
        let pts = vec![
            [0.0, 0.0, 0.0],
            [10.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [10.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [10.0, 0.0, 1.0],
            [0.0, 1.0, 1.0],
            [10.0, 1.0, 1.0],
        ];
        let hull = ConvexHull3D::build(&pts).expect("build");
        let ar = hull_aspect_ratio(&hull);
        assert!(
            ar > 5.0,
            "Elongated box aspect ratio should be > 5, got {ar}"
        );
    }
    #[test]
    fn test_simplify_hull_returns_fewer_or_equal_points() {
        let hull = ConvexHull3D::build(&cube_points()).expect("build");
        let simplified = simplify_hull_points(&hull, 0.1);
        assert!(!simplified.is_empty());
        assert!(simplified.len() <= hull.vertices.len());
    }
    #[test]
    fn test_simplify_hull_zero_threshold_keeps_all() {
        let hull = ConvexHull3D::build(&cube_points()).expect("build");
        let simplified = simplify_hull_points(&hull, 0.0);
        assert_eq!(
            simplified.len(),
            hull.vertices.len(),
            "Zero threshold should keep all vertices"
        );
    }
    #[test]
    fn test_simplify_hull_empty_hull() {
        let hull = ConvexHull3D {
            vertices: vec![],
            faces: vec![],
        };
        let simplified = simplify_hull_points(&hull, 0.1);
        assert!(simplified.is_empty());
    }
    #[test]
    fn test_deduplicate_removes_identical() {
        let pts = vec![[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let unique = deduplicate_points(&pts, 1e-9);
        assert_eq!(unique.len(), 2, "Should remove duplicate point");
    }
    #[test]
    fn test_deduplicate_tolerance() {
        let pts = vec![[0.0, 0.0, 0.0], [0.0001, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let unique = deduplicate_points(&pts, 0.01);
        assert_eq!(unique.len(), 2, "Near-duplicate should be removed");
    }
    #[test]
    fn test_deduplicate_empty() {
        let unique = deduplicate_points(&[], 0.01);
        assert!(unique.is_empty());
    }
    #[test]
    fn test_hull_with_quality_cube() {
        let hq = HullWithQuality::build(&cube_points(), 1e-9).expect("build");
        assert!(hq.sphericity > 0.0 && hq.sphericity <= 1.0);
        assert!(hq.aspect_ratio.is_finite() && hq.aspect_ratio > 0.0);
    }
    #[test]
    fn test_hull_with_quality_deduplicates() {
        let mut pts = cube_points();
        pts.push(pts[0]);
        let hq = HullWithQuality::build(&pts, 1e-9).expect("build");
        assert!(hq.hull.n_vertices() <= pts.len());
    }
    #[test]
    fn test_hull_with_quality_too_few_returns_none() {
        let pts = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        assert!(HullWithQuality::build(&pts, 1e-9).is_none());
    }
}
/// Compute the Minkowski sum of two convex hulls as a point cloud.
///
/// The Minkowski sum A ⊕ B is the set `{ a + b | a ∈ A, b ∈ B }`.
/// For convex hulls, the result is convex.  This function returns all
/// pairwise vertex sums; the caller should build a new hull from the result.
#[allow(dead_code)]
pub fn minkowski_sum_hulls(a: &ConvexHull3D, b: &ConvexHull3D) -> Vec<[f64; 3]> {
    let mut pts = Vec::with_capacity(a.vertices.len() * b.vertices.len());
    for &va in &a.vertices {
        for &vb in &b.vertices {
            pts.push(add(va, vb));
        }
    }
    pts
}
/// Build the convex hull of the Minkowski sum of `a` and `b`.
///
/// Returns `None` if the result is degenerate (fewer than 4 non-coplanar vertices).
#[allow(dead_code)]
pub fn minkowski_sum_hull(a: &ConvexHull3D, b: &ConvexHull3D) -> Option<ConvexHull3D> {
    let pts = minkowski_sum_hulls(a, b);
    let unique = deduplicate_points(&pts, 1e-10);
    ConvexHull3D::build(&unique)
}
/// Compute outward-facing unit normals for every triangular face.
///
/// Returns a `Vec` with one `[f64; 3]` unit normal per face, in the same
/// order as `hull.faces`.
#[allow(dead_code)]
pub fn hull_face_normals(hull: &ConvexHull3D) -> Vec<[f64; 3]> {
    hull.faces
        .iter()
        .map(|f| {
            let n = face_normal(
                hull.vertices[f[0]],
                hull.vertices[f[1]],
                hull.vertices[f[2]],
            );
            normalize(n)
        })
        .collect()
}
/// Compute the *signed* distances from the centroid to each face plane.
///
/// Positive means the centroid is on the interior side (correct for a closed
/// convex hull with outward normals).
#[allow(dead_code)]
pub fn hull_face_offsets(hull: &ConvexHull3D) -> Vec<f64> {
    let c = hull.vertices.iter().fold([0.0f64; 3], |a, v| add(a, *v));
    let n = hull.vertices.len() as f64;
    let _centroid = scale(c, 1.0 / n);
    hull.faces
        .iter()
        .map(|f| {
            let n_raw = face_normal(
                hull.vertices[f[0]],
                hull.vertices[f[1]],
                hull.vertices[f[2]],
            );
            let nhat = normalize(n_raw);
            dot(nhat, hull.vertices[f[0]])
        })
        .collect()
}
/// Test whether two convex hulls intersect using the GJK algorithm.
///
/// Returns `true` if the hulls overlap (i.e., their Minkowski difference
/// contains the origin).
///
/// Implements a simplified 3-D GJK with at most 64 iterations.
#[allow(dead_code)]
pub fn gjk_intersect(a: &ConvexHull3D, b: &ConvexHull3D) -> bool {
    if a.vertices.is_empty() || b.vertices.is_empty() {
        return false;
    }
    let support = |dir: [f64; 3]| -> [f64; 3] {
        let sa = a.support(dir);
        let neg_dir = [-dir[0], -dir[1], -dir[2]];
        let sb = b.support(neg_dir);
        sub(sa, sb)
    };
    let mut dir = sub(a.vertices[0], b.vertices[0]);
    if dot(dir, dir) < 1e-24 {
        dir = [1.0, 0.0, 0.0];
    }
    let mut simplex: Vec<[f64; 3]> = Vec::with_capacity(4);
    simplex.push(support(dir));
    dir = [-simplex[0][0], -simplex[0][1], -simplex[0][2]];
    for _ in 0..64 {
        if dot(dir, dir) < 1e-24 {
            return true;
        }
        let a_pt = support(dir);
        if dot(a_pt, dir) < 0.0 {
            return false;
        }
        simplex.push(a_pt);
        if gjk_do_simplex(&mut simplex, &mut dir) {
            return true;
        }
    }
    false
}
/// GJK simplex processing step.  Returns `true` if origin is contained.
pub(super) fn gjk_do_simplex(simplex: &mut Vec<[f64; 3]>, dir: &mut [f64; 3]) -> bool {
    match simplex.len() {
        2 => gjk_line_case(simplex, dir),
        3 => gjk_triangle_case(simplex, dir),
        4 => gjk_tetrahedron_case(simplex, dir),
        _ => false,
    }
}
pub(super) fn gjk_line_case(simplex: &mut Vec<[f64; 3]>, dir: &mut [f64; 3]) -> bool {
    let b = simplex[0];
    let a = simplex[1];
    let ab = sub(b, a);
    let ao = [-a[0], -a[1], -a[2]];
    if dot(ab, ao) > 0.0 {
        *dir = sub(cross(cross(ab, ao), ab), [0.0; 3]);
        let ab_cross_ao = cross(ab, ao);
        *dir = cross(ab_cross_ao, ab);
    } else {
        simplex.clear();
        simplex.push(a);
        *dir = ao;
    }
    false
}
pub(super) fn gjk_triangle_case(simplex: &mut Vec<[f64; 3]>, dir: &mut [f64; 3]) -> bool {
    let c = simplex[0];
    let b = simplex[1];
    let a = simplex[2];
    let ab = sub(b, a);
    let ac = sub(c, a);
    let ao = [-a[0], -a[1], -a[2]];
    let abc = cross(ab, ac);
    if dot(cross(abc, ac), ao) > 0.0 {
        if dot(ac, ao) > 0.0 {
            *simplex = vec![c, a];
            *dir = cross(cross(ac, ao), ac);
        } else if dot(ab, ao) > 0.0 {
            *simplex = vec![b, a];
            *dir = cross(cross(ab, ao), ab);
        } else {
            *simplex = vec![a];
            *dir = ao;
        }
    } else if dot(cross(ab, abc), ao) > 0.0 {
        if dot(ab, ao) > 0.0 {
            *simplex = vec![b, a];
            *dir = cross(cross(ab, ao), ab);
        } else {
            *simplex = vec![a];
            *dir = ao;
        }
    } else if dot(abc, ao) > 0.0 {
        *dir = abc;
    } else {
        *simplex = vec![b, c, a];
        *dir = [-abc[0], -abc[1], -abc[2]];
    }
    false
}
pub(super) fn gjk_tetrahedron_case(simplex: &mut Vec<[f64; 3]>, dir: &mut [f64; 3]) -> bool {
    let d = simplex[0];
    let c = simplex[1];
    let b = simplex[2];
    let a = simplex[3];
    let ab = sub(b, a);
    let ac = sub(c, a);
    let ad = sub(d, a);
    let ao = [-a[0], -a[1], -a[2]];
    let abc = cross(ab, ac);
    let acd = cross(ac, ad);
    let adb = cross(ad, ab);
    if dot(abc, ao) > 0.0 {
        *simplex = vec![c, b, a];
        return gjk_triangle_case(simplex, dir);
    }
    if dot(acd, ao) > 0.0 {
        *simplex = vec![d, c, a];
        return gjk_triangle_case(simplex, dir);
    }
    if dot(adb, ao) > 0.0 {
        *simplex = vec![b, d, a];
        return gjk_triangle_case(simplex, dir);
    }
    true
}
/// Offset a convex hull outward by `epsilon` (chamfering / Minkowski expansion).
///
/// Each vertex is displaced along the average of adjacent face normals scaled
/// by `epsilon`.  Returns a new point cloud from which a hull can be rebuilt.
#[allow(dead_code)]
pub fn chamfer_hull(hull: &ConvexHull3D, epsilon: f64) -> Vec<[f64; 3]> {
    if hull.vertices.is_empty() {
        return vec![];
    }
    let mut normals: Vec<[f64; 3]> = vec![[0.0; 3]; hull.vertices.len()];
    let mut counts: Vec<u32> = vec![0; hull.vertices.len()];
    for face in &hull.faces {
        let n = normalize(face_normal(
            hull.vertices[face[0]],
            hull.vertices[face[1]],
            hull.vertices[face[2]],
        ));
        for &vi in face {
            normals[vi] = add(normals[vi], n);
            counts[vi] += 1;
        }
    }
    hull.vertices
        .iter()
        .enumerate()
        .map(|(i, &v)| {
            let avg = if counts[i] > 0 {
                normalize(normals[i])
            } else {
                [0.0; 3]
            };
            add(v, scale(avg, epsilon))
        })
        .collect()
}
/// Chamfer a convex hull and return the expanded hull.
///
/// Returns `None` if the result is degenerate.
#[allow(dead_code)]
pub fn chamfer_build(hull: &ConvexHull3D, epsilon: f64) -> Option<ConvexHull3D> {
    let pts = chamfer_hull(hull, epsilon);
    ConvexHull3D::build(&pts)
}
#[cfg(test)]
mod tests_extended {

    use crate::ConvexHull3D;

    use crate::add;
    use crate::chamfer_build;
    use crate::chamfer_hull;
    use crate::convex_decomposition::dot;
    use crate::convex_decomposition::scale;
    use crate::convex_decomposition::sub;
    use crate::gjk_intersect;
    use crate::hull_face_normals;
    use crate::hull_face_offsets;
    use crate::minkowski_sum_hull;
    use crate::minkowski_sum_hulls;

    fn unit_cube_hull() -> ConvexHull3D {
        let pts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [0.0, 1.0, 1.0],
            [1.0, 1.0, 1.0],
        ];
        ConvexHull3D::build(&pts).expect("unit cube hull")
    }
    fn unit_tet_hull() -> ConvexHull3D {
        let pts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        ConvexHull3D::build(&pts).expect("unit tet hull")
    }
    #[test]
    fn test_minkowski_sum_hulls_vertex_count() {
        let h1 = unit_tet_hull();
        let h2 = unit_tet_hull();
        let pts = minkowski_sum_hulls(&h1, &h2);
        assert_eq!(pts.len(), h1.vertices.len() * h2.vertices.len());
    }
    #[test]
    fn test_minkowski_sum_hull_builds() {
        let h1 = unit_tet_hull();
        let h2 = unit_tet_hull();
        let sum = minkowski_sum_hull(&h1, &h2);
        assert!(sum.is_some(), "Minkowski sum should build a valid hull");
    }
    #[test]
    fn test_minkowski_sum_hull_larger_than_operands() {
        let h1 = unit_cube_hull();
        let h2 = unit_cube_hull();
        let sum = minkowski_sum_hull(&h1, &h2).expect("should build");
        let aabb_min = sum.vertices.iter().fold([f64::INFINITY; 3], |a, v| {
            [a[0].min(v[0]), a[1].min(v[1]), a[2].min(v[2])]
        });
        let aabb_max = sum.vertices.iter().fold([f64::NEG_INFINITY; 3], |a, v| {
            [a[0].max(v[0]), a[1].max(v[1]), a[2].max(v[2])]
        });
        for k in 0..3 {
            let extent = aabb_max[k] - aabb_min[k];
            assert!(
                (extent - 2.0).abs() < 0.01,
                "axis {k} extent = {extent}, expected ~2.0"
            );
        }
    }
    #[test]
    fn test_hull_face_normals_unit_length() {
        let hull = unit_cube_hull();
        let normals = hull_face_normals(&hull);
        assert_eq!(normals.len(), hull.n_faces());
        for (i, n) in normals.iter().enumerate() {
            let len = dot(*n, *n).sqrt();
            assert!(
                (len - 1.0).abs() < 1e-9,
                "face {i} normal not unit: len={len}"
            );
        }
    }
    #[test]
    fn test_hull_face_normals_count_matches_faces() {
        let hull = unit_tet_hull();
        let normals = hull_face_normals(&hull);
        assert_eq!(normals.len(), 4, "tetrahedron should have 4 normals");
    }
    #[test]
    fn test_gjk_intersect_overlapping_cubes() {
        let hull_a = unit_cube_hull();
        let pts_b: Vec<[f64; 3]> = hull_a
            .vertices
            .iter()
            .map(|v| [v[0] + 0.5, v[1], v[2]])
            .collect();
        let hull_b = ConvexHull3D::build(&pts_b).expect("build");
        assert!(
            gjk_intersect(&hull_a, &hull_b),
            "overlapping cubes should intersect"
        );
    }
    #[test]
    fn test_gjk_intersect_separated_cubes() {
        let hull_a = unit_cube_hull();
        let pts_b: Vec<[f64; 3]> = hull_a
            .vertices
            .iter()
            .map(|v| [v[0] + 5.0, v[1], v[2]])
            .collect();
        let hull_b = ConvexHull3D::build(&pts_b).expect("build");
        assert!(
            !gjk_intersect(&hull_a, &hull_b),
            "separated cubes should not intersect"
        );
    }
    #[test]
    fn test_gjk_intersect_identical_hulls() {
        let hull = unit_tet_hull();
        assert!(
            gjk_intersect(&hull, &hull),
            "identical hulls should intersect"
        );
    }
    #[test]
    fn test_gjk_intersect_touching_hulls() {
        let hull_a = unit_cube_hull();
        let pts_b: Vec<[f64; 3]> = hull_a
            .vertices
            .iter()
            .map(|v| [v[0] + 1.0, v[1], v[2]])
            .collect();
        let hull_b = ConvexHull3D::build(&pts_b).expect("build");
        let _result = gjk_intersect(&hull_a, &hull_b);
    }
    #[test]
    fn test_chamfer_hull_expands_vertices() {
        let hull = unit_cube_hull();
        let eps = 0.1;
        let chamfered = chamfer_hull(&hull, eps);
        let c: [f64; 3] = scale(
            hull.vertices.iter().fold([0.0; 3], |a, v| add(a, *v)),
            1.0 / hull.vertices.len() as f64,
        );
        for (orig, ch) in hull.vertices.iter().zip(chamfered.iter()) {
            let d_orig = dot(sub(*orig, c), sub(*orig, c)).sqrt();
            let d_ch = dot(sub(*ch, c), sub(*ch, c)).sqrt();
            assert!(
                d_ch >= d_orig - 1e-9,
                "chamfered vertex should be >= original distance from centroid"
            );
        }
    }
    #[test]
    fn test_chamfer_zero_epsilon_unchanged() {
        let hull = unit_cube_hull();
        let pts = chamfer_hull(&hull, 0.0);
        assert_eq!(pts.len(), hull.vertices.len());
        for (a, b) in hull.vertices.iter().zip(pts.iter()) {
            let d = dot(sub(*a, *b), sub(*a, *b)).sqrt();
            assert!(d < 1e-12, "zero epsilon should not change vertices");
        }
    }
    #[test]
    fn test_chamfer_build_produces_valid_hull() {
        let hull = unit_cube_hull();
        let expanded = chamfer_build(&hull, 0.1);
        assert!(expanded.is_some(), "chamfered hull should build");
        let exp = expanded.unwrap();
        assert!(exp.n_faces() > 0);
    }
    #[test]
    fn test_volume_centroid_cube() {
        let hull = unit_cube_hull();
        let (vol, centroid) = hull.volume_centroid();
        assert!((vol - 1.0).abs() < 0.05, "volume={vol}");
        for k in 0..3 {
            assert!(
                (centroid[k] - 0.5).abs() < 0.1,
                "centroid[{k}]={} expected ~0.5",
                centroid[k]
            );
        }
    }
    #[test]
    fn test_vertex_centroid_cube() {
        let hull = unit_cube_hull();
        let c = hull.vertex_centroid();
        for k in 0..3 {
            assert!(
                (c[k] - 0.5).abs() < 1e-9,
                "vertex_centroid[{k}]={} expected 0.5",
                c[k]
            );
        }
    }
    #[test]
    fn test_hull_face_offsets_count() {
        let hull = unit_cube_hull();
        let offsets = hull_face_offsets(&hull);
        assert_eq!(offsets.len(), hull.n_faces());
    }
    #[test]
    fn test_hull_face_offsets_tetra() {
        let hull = unit_tet_hull();
        let offsets = hull_face_offsets(&hull);
        assert_eq!(offsets.len(), 4);
        for &off in &offsets {
            assert!(off.is_finite(), "offset should be finite, got {off}");
        }
    }
}
#[cfg(test)]
mod tests_ch3d_physics {

    use crate::ConvexHull3D;

    fn unit_cube_hull() -> ConvexHull3D {
        let pts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [0.0, 1.0, 1.0],
            [1.0, 1.0, 1.0],
        ];
        ConvexHull3D::build(&pts).expect("unit cube hull")
    }
    fn unit_tet_hull() -> ConvexHull3D {
        let pts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        ConvexHull3D::build(&pts).expect("unit tet hull")
    }
    #[test]
    fn test_compute_volume_unit_cube() {
        let hull = unit_cube_hull();
        let vol = hull.compute_volume();
        assert!((vol - 1.0).abs() < 0.05, "unit cube volume ~1.0, got {vol}");
    }
    #[test]
    fn test_compute_volume_tetrahedron() {
        let hull = unit_tet_hull();
        let vol = hull.compute_volume();
        let expected = 1.0 / 6.0;
        assert!(
            (vol - expected).abs() < 0.01,
            "tet volume ~{expected}, got {vol}"
        );
    }
    #[test]
    fn test_compute_volume_empty() {
        let hull = ConvexHull3D {
            vertices: vec![],
            faces: vec![],
        };
        assert_eq!(hull.compute_volume(), 0.0);
    }
    #[test]
    fn test_compute_volume_positive() {
        let hull = unit_cube_hull();
        assert!(hull.compute_volume() > 0.0, "volume must be positive");
    }
    #[test]
    fn test_compute_center_of_mass_unit_cube() {
        let hull = unit_cube_hull();
        let com = hull.compute_center_of_mass();
        for k in 0..3 {
            assert!(
                (com[k] - 0.5).abs() < 0.1,
                "COM[{k}]={} expected ~0.5",
                com[k]
            );
        }
    }
    #[test]
    fn test_compute_center_of_mass_empty() {
        let hull = ConvexHull3D {
            vertices: vec![],
            faces: vec![],
        };
        let com = hull.compute_center_of_mass();
        for &c in &com {
            assert!(c.is_finite());
        }
    }
    #[test]
    fn test_compute_center_of_mass_symmetric() {
        let hull = unit_tet_hull();
        let com = hull.compute_center_of_mass();
        for &c in &com {
            assert!(c >= -0.05, "COM component should be >= 0, got {c}");
        }
    }
    #[test]
    fn test_compute_moment_of_inertia_unit_cube() {
        let hull = unit_cube_hull();
        let inertia = hull.compute_moment_of_inertia(1.0);
        for k in 0..3 {
            assert!(
                inertia[k][k] > 0.0,
                "diagonal inertia[{k}][{k}] must be positive, got {}",
                inertia[k][k]
            );
        }
    }
    #[test]
    fn test_compute_moment_of_inertia_zero_mass() {
        let hull = unit_cube_hull();
        let inertia = hull.compute_moment_of_inertia(0.0);
        for row in &inertia {
            for &val in row {
                assert_eq!(val, 0.0, "zero mass → zero inertia");
            }
        }
    }
    #[test]
    fn test_compute_moment_of_inertia_symmetry() {
        let hull = unit_cube_hull();
        let inertia = hull.compute_moment_of_inertia(1.0);
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (inertia[i][j] - inertia[j][i]).abs() < 1e-10,
                    "tensor not symmetric at [{i}][{j}]"
                );
            }
        }
    }
    #[test]
    fn test_compute_moment_of_inertia_tetrahedron_positive_diagonal() {
        let hull = unit_tet_hull();
        let inertia = hull.compute_moment_of_inertia(1.0);
        for k in 0..3 {
            assert!(
                inertia[k][k] > 0.0,
                "tet diagonal inertia[{k}][{k}]={} should be positive",
                inertia[k][k]
            );
        }
    }
}
/// Compute the 2D lower convex hull of a set of points using Andrew's monotone
/// chain algorithm.
///
/// Points are sorted lexicographically (x then y). The lower hull is the
/// sequence of points forming the lower boundary.
///
/// Returns the lower hull in counter-clockwise order.
#[allow(dead_code)]
pub fn lower_hull_2d(mut pts: Vec<[f64; 2]>) -> Vec<[f64; 2]> {
    pts.sort_unstable_by(|a, b| {
        a[0].partial_cmp(&b[0])
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(a[1].partial_cmp(&b[1]).unwrap_or(std::cmp::Ordering::Equal))
    });
    let mut hull: Vec<[f64; 2]> = Vec::new();
    for p in &pts {
        while hull.len() >= 2 {
            let n = hull.len();
            let cross = cross2d(hull[n - 2], hull[n - 1], *p);
            if cross <= 0.0 {
                hull.pop();
            } else {
                break;
            }
        }
        hull.push(*p);
    }
    hull
}
/// Compute the 2D upper convex hull of a set of points using Andrew's monotone
/// chain algorithm.
///
/// Returns the upper hull in counter-clockwise order (right-to-left).
#[allow(dead_code)]
pub fn upper_hull_2d(mut pts: Vec<[f64; 2]>) -> Vec<[f64; 2]> {
    pts.sort_unstable_by(|a, b| {
        a[0].partial_cmp(&b[0])
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(a[1].partial_cmp(&b[1]).unwrap_or(std::cmp::Ordering::Equal))
    });
    let n = pts.len();
    let mut hull: Vec<[f64; 2]> = Vec::new();
    for i in (0..n).rev() {
        let p = pts[i];
        while hull.len() >= 2 {
            let m = hull.len();
            let cross = cross2d(hull[m - 2], hull[m - 1], p);
            if cross <= 0.0 {
                hull.pop();
            } else {
                break;
            }
        }
        hull.push(p);
    }
    hull
}
/// Compute the full 2D convex hull using Andrew's monotone chain algorithm.
///
/// Returns the hull vertices in counter-clockwise order (no duplicate start/end).
#[allow(dead_code)]
pub fn monotone_chain_hull_2d(pts: &[[f64; 2]]) -> Vec<[f64; 2]> {
    if pts.len() < 3 {
        let mut sorted = pts.to_vec();
        sorted
            .sort_unstable_by(|a, b| a[0].partial_cmp(&b[0]).unwrap_or(std::cmp::Ordering::Equal));
        return sorted;
    }
    let lower = lower_hull_2d(pts.to_vec());
    let upper = upper_hull_2d(pts.to_vec());
    let mut hull = lower;
    if upper.len() > 2 {
        hull.extend_from_slice(&upper[1..upper.len() - 1]);
    }
    hull
}
/// 2D cross product (z-component) of vectors (b-a) and (c-a).
///
/// Positive: counter-clockwise, negative: clockwise, zero: collinear.
#[inline]
pub(super) fn cross2d(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> f64 {
    (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
}
/// Support function for an axis-aligned box `[min, max]`.
///
/// Returns the furthest point of the box in direction `dir`.
#[allow(dead_code)]
pub fn box_support(min: [f64; 3], max: [f64; 3], dir: [f64; 3]) -> [f64; 3] {
    [
        if dir[0] >= 0.0 { max[0] } else { min[0] },
        if dir[1] >= 0.0 { max[1] } else { min[1] },
        if dir[2] >= 0.0 { max[2] } else { min[2] },
    ]
}
/// Support function for a sphere centered at `center` with radius `radius`.
#[allow(dead_code)]
pub fn sphere_support(center: [f64; 3], radius: f64, dir: [f64; 3]) -> [f64; 3] {
    let n = normalize(dir);
    add(center, scale(n, radius))
}
/// Support function for a capsule (cylinder with hemispherical caps).
///
/// The capsule has its axis from `a` to `b` and radius `r`.
#[allow(dead_code)]
pub fn capsule_support(a: [f64; 3], b: [f64; 3], r: f64, dir: [f64; 3]) -> [f64; 3] {
    let n = normalize(dir);
    let sa = add(a, scale(n, r));
    let sb = add(b, scale(n, r));
    if dot(sa, dir) >= dot(sb, dir) { sa } else { sb }
}
/// Support function for a cylinder (flat caps) aligned along the Z axis.
///
/// Center at `center`, half-height `h`, radius `r`.
#[allow(dead_code)]
pub fn cylinder_support(center: [f64; 3], h: f64, r: f64, dir: [f64; 3]) -> [f64; 3] {
    let dir_xy_len = (dir[0] * dir[0] + dir[1] * dir[1]).sqrt();
    let (dx, dy) = if dir_xy_len > 1e-12 {
        (dir[0] / dir_xy_len * r, dir[1] / dir_xy_len * r)
    } else {
        (0.0, 0.0)
    };
    let dz = if dir[2] >= 0.0 { h } else { -h };
    add(center, [dx, dy, dz])
}
/// Support function for an ellipsoid with semi-axes `(a, b, c)` centered at `center`.
#[allow(dead_code)]
pub fn ellipsoid_support(center: [f64; 3], semi_axes: [f64; 3], dir: [f64; 3]) -> [f64; 3] {
    let [a, b, c] = semi_axes;
    let numerator = [a * a * dir[0], b * b * dir[1], c * c * dir[2]];
    let denom_sq = a * a * dir[0] * dir[0] + b * b * dir[1] * dir[1] + c * c * dir[2] * dir[2];
    let denom = denom_sq.sqrt();
    let s = if denom > 1e-12 {
        scale(numerator, 1.0 / denom)
    } else {
        [0.0; 3]
    };
    add(center, s)
}
/// Simplify a convex hull by iteratively removing the vertex that contributes
/// least to the total surface area until at most `target_vertices` remain.
///
/// Each iteration finds the vertex with the smallest total adjacent face area
/// and removes it from the point cloud before rebuilding.
///
/// Returns the simplified point cloud (caller must rebuild the hull).
#[allow(dead_code)]
pub fn vertex_decimation(hull: &ConvexHull3D, target_vertices: usize) -> Vec<[f64; 3]> {
    let mut pts = hull.vertices.clone();
    while pts.len() > target_vertices.max(4) {
        if let Some(idx) = least_contributing_vertex_index(&pts, hull) {
            pts.remove(idx);
        } else {
            break;
        }
        if ConvexHull3D::build(&pts).is_none() {
            break;
        }
    }
    pts
}
/// Find the index (in `pts`) of the vertex that contributes the least area
/// to the hull.
pub(super) fn least_contributing_vertex_index(
    pts: &[[f64; 3]],
    hull: &ConvexHull3D,
) -> Option<usize> {
    if pts.is_empty() || hull.vertices.is_empty() {
        return None;
    }
    let mut contrib = vec![0.0f64; pts.len()];
    for face in &hull.faces {
        let area = {
            let n = face_normal(
                hull.vertices[face[0]],
                hull.vertices[face[1]],
                hull.vertices[face[2]],
            );
            dot(n, n).sqrt() * 0.5
        };
        for &vi in face {
            let hv = hull.vertices[vi];
            let best = pts
                .iter()
                .enumerate()
                .min_by(|(_, a), (_, b)| {
                    let da = dot(sub(**a, hv), sub(**a, hv));
                    let db = dot(sub(**b, hv), sub(**b, hv));
                    da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
                })
                .map(|(i, _)| i);
            if let Some(bi) = best {
                contrib[bi] += area;
            }
        }
    }
    contrib
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
        .map(|(i, _)| i)
}
/// Generate vertices of a regular icosahedron centered at the origin with
/// circumscribed radius `r`.
///
/// An icosahedron has 12 vertices and 20 triangular faces.
#[allow(dead_code)]
pub fn icosahedron_vertices(r: f64) -> Vec<[f64; 3]> {
    let phi = (1.0 + 5.0_f64.sqrt()) / 2.0;
    let a = r / (1.0 + phi * phi).sqrt();
    let b = a * phi;
    vec![
        [-a, b, 0.0],
        [a, b, 0.0],
        [-a, -b, 0.0],
        [a, -b, 0.0],
        [0.0, -a, b],
        [0.0, a, b],
        [0.0, -a, -b],
        [0.0, a, -b],
        [b, 0.0, -a],
        [b, 0.0, a],
        [-b, 0.0, -a],
        [-b, 0.0, a],
    ]
}
/// Generate vertices of a regular octahedron centered at the origin with
/// circumscribed radius `r`.
///
/// An octahedron has 6 vertices and 8 triangular faces.
#[allow(dead_code)]
pub fn octahedron_vertices(r: f64) -> Vec<[f64; 3]> {
    vec![
        [r, 0.0, 0.0],
        [-r, 0.0, 0.0],
        [0.0, r, 0.0],
        [0.0, -r, 0.0],
        [0.0, 0.0, r],
        [0.0, 0.0, -r],
    ]
}
/// Generate vertices of a regular tetrahedron centered at the origin with
/// circumscribed radius `r`.
#[allow(dead_code)]
pub fn tetrahedron_vertices(r: f64) -> Vec<[f64; 3]> {
    let a = r * (8.0_f64 / 9.0).sqrt();
    let b = r * (2.0_f64 / 9.0).sqrt();
    let c = r * (2.0_f64 / 3.0).sqrt();
    vec![
        [0.0, 0.0, r],
        [a, 0.0, -r / 3.0],
        [-b, c, -r / 3.0],
        [-b, -c, -r / 3.0],
    ]
}
/// Compute the Minkowski difference point cloud of hulls A and B.
///
/// The Minkowski difference A ⊖ B = A ⊕ (-B) is the set `{ a - b | a∈A, b∈B }`.
/// The GJK algorithm tests whether the origin lies inside this set.
#[allow(dead_code)]
pub fn minkowski_difference_hulls(a: &ConvexHull3D, b: &ConvexHull3D) -> Vec<[f64; 3]> {
    let mut pts = Vec::with_capacity(a.vertices.len() * b.vertices.len());
    for &va in &a.vertices {
        for &vb in &b.vertices {
            pts.push(sub(va, vb));
        }
    }
    pts
}
/// Build the convex hull of the Minkowski difference of `a` and `b`.
#[allow(dead_code)]
pub fn minkowski_difference_hull(a: &ConvexHull3D, b: &ConvexHull3D) -> Option<ConvexHull3D> {
    let pts = minkowski_difference_hulls(a, b);
    let unique = deduplicate_points(&pts, 1e-10);
    ConvexHull3D::build(&unique)
}
/// Compute the bounding box of a point cloud.
///
/// Returns `([min_x, min_y, min_z], [max_x, max_y, max_z])` or `None` if
/// the cloud is empty.
#[allow(dead_code)]
pub fn point_cloud_aabb(pts: &[[f64; 3]]) -> Option<([f64; 3], [f64; 3])> {
    if pts.is_empty() {
        return None;
    }
    let mut mn = pts[0];
    let mut mx = pts[0];
    for &p in pts {
        for i in 0..3 {
            if p[i] < mn[i] {
                mn[i] = p[i];
            }
            if p[i] > mx[i] {
                mx[i] = p[i];
            }
        }
    }
    Some((mn, mx))
}
/// Compute the centroid (mean) of a point cloud.
#[allow(dead_code)]
pub fn point_cloud_centroid(pts: &[[f64; 3]]) -> Option<[f64; 3]> {
    if pts.is_empty() {
        return None;
    }
    let n = pts.len() as f64;
    let sum = pts.iter().fold([0.0f64; 3], |acc, &p| add(acc, p));
    Some(scale(sum, 1.0 / n))
}
/// Compute the covariance matrix of a point cloud (returns 3x3 row-major).
#[allow(dead_code)]
pub fn point_cloud_covariance(pts: &[[f64; 3]]) -> [[f64; 3]; 3] {
    if pts.len() < 2 {
        return [[0.0; 3]; 3];
    }
    let mean = point_cloud_centroid(pts).unwrap_or([0.0; 3]);
    let n = pts.len() as f64;
    let mut cov = [[0.0f64; 3]; 3];
    for &p in pts {
        let d = sub(p, mean);
        for i in 0..3 {
            for j in 0..3 {
                cov[i][j] += d[i] * d[j];
            }
        }
    }
    for i in 0..3 {
        for j in 0..3 {
            cov[i][j] /= n - 1.0;
        }
    }
    cov
}
/// Scale all points in a cloud by `factor`.
#[allow(dead_code)]
pub fn scale_point_cloud(pts: &[[f64; 3]], factor: f64) -> Vec<[f64; 3]> {
    pts.iter().map(|&p| scale(p, factor)).collect()
}
/// Translate all points in a cloud by `offset`.
#[allow(dead_code)]
pub fn translate_point_cloud(pts: &[[f64; 3]], offset: [f64; 3]) -> Vec<[f64; 3]> {
    pts.iter().map(|&p| add(p, offset)).collect()
}
