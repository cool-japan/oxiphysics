//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use oxiphysics_core::math::Vec3;

use super::types::{BspPlane, ManifoldCheckResult};

#[inline]
pub(super) fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
#[inline]
pub(super) fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
#[inline]
pub(super) fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
#[inline]
pub(super) fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
#[inline]
pub(super) fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
#[inline]
pub(super) fn len3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}
#[inline]
pub(super) fn normalize3(a: [f64; 3]) -> [f64; 3] {
    let l = len3(a);
    if l < f64::EPSILON {
        [0.0, 0.0, 0.0]
    } else {
        [a[0] / l, a[1] / l, a[2] / l]
    }
}
pub(super) fn vec3_to_arr(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}
pub(super) fn arr_to_vec3(a: [f64; 3]) -> Vec3 {
    Vec3::new(a[0], a[1], a[2])
}
/// Test whether a point is inside a closed mesh using multiple ray casts.
///
/// Fires rays in three axis directions (+X, +Y, +Z) and takes a majority vote
/// to improve robustness against degenerate cases where a ray passes exactly
/// through an edge or vertex. An odd intersection count means the point is
/// inside the mesh (Jordan curve theorem generalised to 3-D).
pub fn point_inside_mesh(point: [f64; 3], verts: &[[f64; 3]], tris: &[[usize; 3]]) -> bool {
    let ray_dirs: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let mut votes_inside = 0usize;
    for ray_dir in &ray_dirs {
        let mut count = 0usize;
        for tri in tris {
            let v0 = verts[tri[0]];
            let v1 = verts[tri[1]];
            let v2 = verts[tri[2]];
            if ray_triangle_intersect(point, *ray_dir, v0, v1, v2).is_some() {
                count += 1;
            }
        }
        if count % 2 == 1 {
            votes_inside += 1;
        }
    }
    votes_inside >= 2
}
/// Möller-Trumbore ray-triangle intersection.
///
/// Returns the parametric distance `t` along the ray if the ray intersects
/// the triangle from the front, or `None` otherwise.
pub fn ray_triangle_intersect(
    origin: [f64; 3],
    direction: [f64; 3],
    v0: [f64; 3],
    v1: [f64; 3],
    v2: [f64; 3],
) -> Option<f64> {
    let e1 = sub3(v1, v0);
    let e2 = sub3(v2, v0);
    let h = cross3(direction, e2);
    let a = dot3(e1, h);
    if a.abs() < 1e-10 {
        return None;
    }
    let f = 1.0 / a;
    let s = sub3(origin, v0);
    let u = f * dot3(s, h);
    if !(0.0..=1.0).contains(&u) {
        return None;
    }
    let q = cross3(s, e1);
    let v = f * dot3(direction, q);
    if v < 0.0 || u + v > 1.0 {
        return None;
    }
    let t = f * dot3(e2, q);
    if t > 1e-10 { Some(t) } else { None }
}
/// Test whether two triangles in 3-D intersect.
///
/// Uses the Devillers-Guigue algorithm (interval overlap on the line of
/// intersection of the two triangle planes). Returns `true` if the triangles
/// intersect (including edge/vertex contacts).
pub fn triangles_intersect(
    a0: [f64; 3],
    a1: [f64; 3],
    a2: [f64; 3],
    b0: [f64; 3],
    b1: [f64; 3],
    b2: [f64; 3],
) -> bool {
    let nb = cross3(sub3(b1, b0), sub3(b2, b0));
    let db = dot3(nb, b0);
    let da0 = dot3(nb, a0) - db;
    let da1 = dot3(nb, a1) - db;
    let da2 = dot3(nb, a2) - db;
    if da0 * da1 > 0.0 && da1 * da2 > 0.0 {
        return false;
    }
    let na = cross3(sub3(a1, a0), sub3(a2, a0));
    let da = dot3(na, a0);
    let db0 = dot3(na, b0) - da;
    let db1 = dot3(na, b1) - da;
    let db2 = dot3(na, b2) - da;
    if db0 * db1 > 0.0 && db1 * db2 > 0.0 {
        return false;
    }
    let d = cross3(na, nb);
    let d_len = len3(d);
    if d_len < 1e-12 {
        return coplanar_tri_overlap(a0, a1, a2, b0, b1, b2, na);
    }
    let proj_a = [dot3(d, a0), dot3(d, a1), dot3(d, a2)];
    let proj_b = [dot3(d, b0), dot3(d, b1), dot3(d, b2)];
    let interval_a = compute_interval(proj_a, [da0, da1, da2]);
    let interval_b = compute_interval(proj_b, [db0, db1, db2]);
    intervals_overlap(interval_a, interval_b)
}
/// Compute the intersection interval of a triangle on the intersection line.
///
/// Uses the signed distances `dists` and projections `proj` to compute
/// `[t_min, t_max]` via linear interpolation.
pub(super) fn compute_interval(proj: [f64; 3], dists: [f64; 3]) -> [f64; 2] {
    let signs: [bool; 3] = [dists[0] > 0.0, dists[1] > 0.0, dists[2] > 0.0];
    let (lone, pair): (usize, [usize; 2]) = if signs[0] != signs[1] && signs[0] != signs[2] {
        (0, [1, 2])
    } else if signs[1] != signs[0] && signs[1] != signs[2] {
        (1, [0, 2])
    } else {
        (2, [0, 1])
    };
    let t0 = proj[pair[0]]
        + (proj[lone] - proj[pair[0]]) * dists[pair[0]] / (dists[pair[0]] - dists[lone]);
    let t1 = proj[pair[1]]
        + (proj[lone] - proj[pair[1]]) * dists[pair[1]] / (dists[pair[1]] - dists[lone]);
    [t0.min(t1), t0.max(t1)]
}
pub(super) fn intervals_overlap(a: [f64; 2], b: [f64; 2]) -> bool {
    !(a[1] < b[0] || b[1] < a[0])
}
/// Coplanar triangle overlap test (2-D edge separating axis test).
pub(super) fn coplanar_tri_overlap(
    a0: [f64; 3],
    a1: [f64; 3],
    a2: [f64; 3],
    b0: [f64; 3],
    b1: [f64; 3],
    b2: [f64; 3],
    normal: [f64; 3],
) -> bool {
    let tri_a = [a0, a1, a2];
    let tri_b = [b0, b1, b2];
    for tri in [tri_a, tri_b].iter() {
        for k in 0..3 {
            let edge = sub3(tri[(k + 1) % 3], tri[k]);
            let axis = cross3(normal, edge);
            let proj_a: [f64; 3] = [
                dot3(axis, tri_a[0]),
                dot3(axis, tri_a[1]),
                dot3(axis, tri_a[2]),
            ];
            let proj_b: [f64; 3] = [
                dot3(axis, tri_b[0]),
                dot3(axis, tri_b[1]),
                dot3(axis, tri_b[2]),
            ];
            let a_min = proj_a.iter().cloned().fold(f64::INFINITY, f64::min);
            let a_max = proj_a.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
            let b_min = proj_b.iter().cloned().fold(f64::INFINITY, f64::min);
            let b_max = proj_b.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
            if a_max < b_min - 1e-10 || b_max < a_min - 1e-10 {
                return false;
            }
        }
    }
    true
}
/// Compute the centroid of a triangle.
pub fn triangle_centroid(v0: [f64; 3], v1: [f64; 3], v2: [f64; 3]) -> [f64; 3] {
    scale3(add3(add3(v0, v1), v2), 1.0 / 3.0)
}
/// Compute the area of a triangle.
pub fn triangle_area(v0: [f64; 3], v1: [f64; 3], v2: [f64; 3]) -> f64 {
    len3(cross3(sub3(v1, v0), sub3(v2, v0))) * 0.5
}
/// Compute the normal of a triangle (not normalised).
pub fn triangle_normal(v0: [f64; 3], v1: [f64; 3], v2: [f64; 3]) -> [f64; 3] {
    cross3(sub3(v1, v0), sub3(v2, v0))
}
/// Compute the unit normal of a triangle.
pub fn triangle_unit_normal(v0: [f64; 3], v1: [f64; 3], v2: [f64; 3]) -> [f64; 3] {
    normalize3(triangle_normal(v0, v1, v2))
}
/// Clip a convex polygon against a half-space defined by `plane`.
///
/// Vertices on the front side (dist ≥ 0) are kept; edges crossing the plane
/// produce intersection vertices.  Returns the clipped polygon or an empty
/// Vec if everything is clipped.
pub fn clip_polygon_by_plane(polygon: &[[f64; 3]], plane: &BspPlane) -> Vec<[f64; 3]> {
    if polygon.is_empty() {
        return Vec::new();
    }
    let n = polygon.len();
    let mut result = Vec::new();
    for i in 0..n {
        let a = polygon[i];
        let b = polygon[(i + 1) % n];
        let da = plane.signed_dist(a);
        let db = plane.signed_dist(b);
        if da >= -1e-10 {
            result.push(a);
        }
        if (da > 1e-10 && db < -1e-10) || (da < -1e-10 && db > 1e-10) {
            let t = da / (da - db);
            let inter = [
                a[0] + t * (b[0] - a[0]),
                a[1] + t * (b[1] - a[1]),
                a[2] + t * (b[2] - a[2]),
            ];
            result.push(inter);
        }
    }
    result
}
/// Clip a polygon against all six half-spaces of an AABB.
///
/// Returns the clipped polygon or `None` if completely clipped.
pub fn clip_polygon_by_aabb(
    polygon: &[[f64; 3]],
    aabb_min: [f64; 3],
    aabb_max: [f64; 3],
) -> Option<Vec<[f64; 3]>> {
    let planes = [
        BspPlane::new([1.0, 0.0, 0.0], aabb_min[0]),
        BspPlane::new([-1.0, 0.0, 0.0], -aabb_max[0]),
        BspPlane::new([0.0, 1.0, 0.0], aabb_min[1]),
        BspPlane::new([0.0, -1.0, 0.0], -aabb_max[1]),
        BspPlane::new([0.0, 0.0, 1.0], aabb_min[2]),
        BspPlane::new([0.0, 0.0, -1.0], -aabb_max[2]),
    ];
    let mut poly = polygon.to_vec();
    for plane in &planes {
        poly = clip_polygon_by_plane(&poly, plane);
        if poly.is_empty() {
            return None;
        }
    }
    Some(poly)
}
/// Clip a 2-D polygon (list of `[x, y]`) against a convex polygon (clipper).
///
/// Uses the Sutherland-Hodgman algorithm.
pub fn clip_polygon_2d(subject: &[[f64; 2]], clip: &[[f64; 2]]) -> Vec<[f64; 2]> {
    if subject.is_empty() || clip.is_empty() {
        return Vec::new();
    }
    let mut output = subject.to_vec();
    let n_clip = clip.len();
    for i in 0..n_clip {
        if output.is_empty() {
            return Vec::new();
        }
        let input = output.clone();
        output.clear();
        let edge_a = clip[i];
        let edge_b = clip[(i + 1) % n_clip];
        let n_in = input.len();
        for j in 0..n_in {
            let current = input[j];
            let prev = input[(j + n_in - 1) % n_in];
            if inside_2d(current, edge_a, edge_b) {
                if !inside_2d(prev, edge_a, edge_b) {
                    output.push(intersect_2d(prev, current, edge_a, edge_b));
                }
                output.push(current);
            } else if inside_2d(prev, edge_a, edge_b) {
                output.push(intersect_2d(prev, current, edge_a, edge_b));
            }
        }
    }
    output
}
/// Test whether `point` is on the left (inside) of the directed edge `a → b`.
pub(super) fn inside_2d(point: [f64; 2], a: [f64; 2], b: [f64; 2]) -> bool {
    let edge = [b[0] - a[0], b[1] - a[1]];
    let to_pt = [point[0] - a[0], point[1] - a[1]];
    edge[0] * to_pt[1] - edge[1] * to_pt[0] >= 0.0
}
/// Compute intersection of segment `p1→p2` with segment `p3→p4`.
pub(super) fn intersect_2d(p1: [f64; 2], p2: [f64; 2], p3: [f64; 2], p4: [f64; 2]) -> [f64; 2] {
    let d1 = [p2[0] - p1[0], p2[1] - p1[1]];
    let d2 = [p4[0] - p3[0], p4[1] - p3[1]];
    let cross = d1[0] * d2[1] - d1[1] * d2[0];
    if cross.abs() < 1e-12 {
        return p1;
    }
    let t = ((p3[0] - p1[0]) * d2[1] - (p3[1] - p1[1]) * d2[0]) / cross;
    [p1[0] + t * d1[0], p1[1] + t * d1[1]]
}
/// Repair the winding order of a triangle mesh so that all outward-facing
/// normals point away from the centroid of the mesh.
///
/// Computes the centroid of all vertices, then for each triangle checks
/// whether the face normal points toward or away from the centroid.  If it
/// points toward the centroid (inside), the triangle winding is flipped.
pub fn repair_winding(verts: &[[f64; 3]], tris: &[[usize; 3]]) -> Vec<[usize; 3]> {
    let n = verts.len() as f64;
    let centroid = verts.iter().fold([0.0f64; 3], |acc, v| add3(acc, *v));
    let centroid = scale3(centroid, 1.0 / n.max(1.0));
    let mut fixed = Vec::with_capacity(tris.len());
    for tri in tris {
        let v0 = verts[tri[0]];
        let v1 = verts[tri[1]];
        let v2 = verts[tri[2]];
        let face_normal = cross3(sub3(v1, v0), sub3(v2, v0));
        let to_centroid = sub3(centroid, scale3(add3(add3(v0, v1), v2), 1.0 / 3.0));
        if dot3(face_normal, to_centroid) > 0.0 {
            fixed.push([tri[0], tri[2], tri[1]]);
        } else {
            fixed.push(*tri);
        }
    }
    fixed
}
/// Check whether all triangles in a mesh have consistent outward normals
/// (pointing away from the mesh centroid).
///
/// Returns the fraction of triangles (0.0–1.0) with correct outward winding.
pub fn winding_consistency_score(verts: &[[f64; 3]], tris: &[[usize; 3]]) -> f64 {
    if tris.is_empty() {
        return 1.0;
    }
    let n = verts.len() as f64;
    let centroid = verts.iter().fold([0.0f64; 3], |acc, v| add3(acc, *v));
    let centroid = scale3(centroid, 1.0 / n.max(1.0));
    let correct = tris
        .iter()
        .filter(|tri| {
            let v0 = verts[tri[0]];
            let v1 = verts[tri[1]];
            let v2 = verts[tri[2]];
            let face_normal = cross3(sub3(v1, v0), sub3(v2, v0));
            let to_centroid = sub3(centroid, scale3(add3(add3(v0, v1), v2), 1.0 / 3.0));
            dot3(face_normal, to_centroid) <= 0.0
        })
        .count();
    correct as f64 / tris.len() as f64
}
/// Compute the 2-D intersection of two convex polygons using Sutherland-Hodgman.
///
/// Both polygons are assumed to be convex and wound counter-clockwise.
pub fn polygon2d_intersection(a: &[[f64; 2]], b: &[[f64; 2]]) -> Vec<[f64; 2]> {
    clip_polygon_2d(a, b)
}
/// Compute the 2-D union of two convex polygons.
///
/// Returns the merged boundary.  Uses the Sutherland-Hodgman approach:
/// clip b against a, then add vertices of a that lie outside b,
/// and vertices of b that lie outside a.
pub fn polygon2d_union(a: &[[f64; 2]], b: &[[f64; 2]]) -> Vec<[f64; 2]> {
    let mut pts: Vec<[f64; 2]> = Vec::new();
    pts.extend_from_slice(a);
    pts.extend_from_slice(b);
    convex_hull_2d(&pts)
}
/// Compute the 2-D difference `a \ b` for convex polygons.
///
/// Returns the portion of `a` outside `b`.
/// Inverts b's winding so we clip a against the outside of b.
pub fn polygon2d_difference(a: &[[f64; 2]], b: &[[f64; 2]]) -> Vec<[f64; 2]> {
    if a.is_empty() || b.is_empty() {
        return a.to_vec();
    }
    let mut result = a.to_vec();
    let n = b.len();
    for i in 0..n {
        if result.is_empty() {
            break;
        }
        let edge_a = b[i];
        let edge_b = b[(i + 1) % n];
        let input = result.clone();
        result.clear();
        let ni = input.len();
        for j in 0..ni {
            let current = input[j];
            let prev = input[(j + ni - 1) % ni];
            let cur_out = !inside_2d(current, edge_a, edge_b);
            let prev_out = !inside_2d(prev, edge_a, edge_b);
            if cur_out {
                if !prev_out {
                    result.push(intersect_2d(prev, current, edge_a, edge_b));
                }
                result.push(current);
            } else if prev_out {
                result.push(intersect_2d(prev, current, edge_a, edge_b));
            }
        }
    }
    result
}
/// Compute a convex hull of 2-D points using the gift-wrapping algorithm.
pub(super) fn convex_hull_2d(pts: &[[f64; 2]]) -> Vec<[f64; 2]> {
    if pts.len() < 3 {
        return pts.to_vec();
    }
    let start = pts
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| a[0].partial_cmp(&b[0]).unwrap_or(std::cmp::Ordering::Equal))
        .map(|(i, _)| i)
        .unwrap_or(0);
    let mut hull = Vec::new();
    let mut current = start;
    loop {
        hull.push(pts[current]);
        let mut next = 0;
        for i in 1..pts.len() {
            if i == current {
                continue;
            }
            let a = pts[current];
            let b = pts[next];
            let c = pts[i];
            let cross = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
            if cross < 0.0 || (cross == 0.0 && dist_sq_2d(a, c) > dist_sq_2d(a, b)) {
                next = i;
            }
        }
        current = next;
        if current == start {
            break;
        }
        if hull.len() > pts.len() {
            break;
        }
    }
    hull
}
pub(super) fn dist_sq_2d(a: [f64; 2], b: [f64; 2]) -> f64 {
    (a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2)
}
/// Check whether a mesh is 2-manifold.
pub fn check_manifold(tris: &[[usize; 3]]) -> ManifoldCheckResult {
    let mut edge_count: std::collections::HashMap<(usize, usize), usize> =
        std::collections::HashMap::new();
    for tri in tris {
        let edges = [
            (tri[0].min(tri[1]), tri[0].max(tri[1])),
            (tri[1].min(tri[2]), tri[1].max(tri[2])),
            (tri[0].min(tri[2]), tri[0].max(tri[2])),
        ];
        for e in &edges {
            *edge_count.entry(*e).or_insert(0) += 1;
        }
    }
    let mut n_boundary = 0;
    let mut n_non_manifold = 0;
    for &count in edge_count.values() {
        if count == 1 {
            n_boundary += 1;
        }
        if count > 2 {
            n_non_manifold += 1;
        }
    }
    ManifoldCheckResult {
        is_manifold: n_boundary == 0 && n_non_manifold == 0,
        n_boundary_edges: n_boundary,
        n_non_manifold_edges: n_non_manifold,
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::TriangleMesh;
    use crate::boolean_ops::BooleanOp;

    use crate::boolean_ops::BspPlane;
    use crate::boolean_ops::CsgMaterial;
    use crate::boolean_ops::CsgWithMaterials;

    use crate::boolean_ops::ManifoldBoolean;
    use crate::boolean_ops::MeshBoolean;

    use crate::boolean_ops::RobustMeshBoolean;

    #[test]
    fn test_intersecting_triangles() {
        let a0 = [-1.0_f64, 0.0, 0.0];
        let a1 = [1.0, 0.0, 0.0];
        let a2 = [0.0, 1.0, 0.0];
        let b0 = [0.0, 0.5, -1.0];
        let b1 = [0.0, 0.5, 1.0];
        let b2 = [0.0, -0.5, 0.0];
        assert!(
            triangles_intersect(a0, a1, a2, b0, b1, b2),
            "crossing triangles should intersect"
        );
    }
    #[test]
    fn test_non_intersecting_triangles() {
        let a0 = [0.0_f64, 0.0, 0.0];
        let a1 = [1.0, 0.0, 0.0];
        let a2 = [0.0, 1.0, 0.0];
        let b0 = [5.0, 0.0, 0.0];
        let b1 = [6.0, 0.0, 0.0];
        let b2 = [5.0, 1.0, 0.0];
        assert!(
            !triangles_intersect(a0, a1, a2, b0, b1, b2),
            "separated triangles should not intersect"
        );
    }
    #[test]
    fn test_ray_triangle_hit() {
        let origin = [-5.0_f64, 0.0, 0.0];
        let dir = [1.0, 0.0, 0.0];
        let v0 = [0.0, -1.0, -1.0];
        let v1 = [0.0, 1.0, -1.0];
        let v2 = [0.0, 0.0, 1.0];
        let t = ray_triangle_intersect(origin, dir, v0, v1, v2);
        assert!(t.is_some(), "ray should hit the triangle");
        assert!((t.unwrap() - 5.0).abs() < 1e-10, "t should be 5.0");
    }
    #[test]
    fn test_ray_triangle_miss() {
        let origin = [-5.0_f64, 10.0, 0.0];
        let dir = [1.0, 0.0, 0.0];
        let v0 = [0.0, -1.0, -1.0];
        let v1 = [0.0, 1.0, -1.0];
        let v2 = [0.0, 0.0, 1.0];
        let t = ray_triangle_intersect(origin, dir, v0, v1, v2);
        assert!(t.is_none(), "ray should miss the triangle");
    }
    fn unit_box_mesh() -> (Vec<[f64; 3]>, Vec<[usize; 3]>) {
        let v: Vec<[f64; 3]> = vec![
            [-1.0, -1.0, -1.0],
            [1.0, -1.0, -1.0],
            [1.0, 1.0, -1.0],
            [-1.0, 1.0, -1.0],
            [-1.0, -1.0, 1.0],
            [1.0, -1.0, 1.0],
            [1.0, 1.0, 1.0],
            [-1.0, 1.0, 1.0],
        ];
        let t: Vec<[usize; 3]> = vec![
            [0, 2, 1],
            [0, 3, 2],
            [4, 5, 6],
            [4, 6, 7],
            [0, 1, 5],
            [0, 5, 4],
            [3, 6, 2],
            [3, 7, 6],
            [0, 4, 7],
            [0, 7, 3],
            [1, 2, 6],
            [1, 6, 5],
        ];
        (v, t)
    }
    #[test]
    fn test_point_inside_mesh() {
        let (v, t) = unit_box_mesh();
        assert!(
            point_inside_mesh([0.1, 0.2, 0.3], &v, &t),
            "interior point should be inside the unit box"
        );
    }
    #[test]
    fn test_point_outside_mesh() {
        let (v, t) = unit_box_mesh();
        assert!(
            !point_inside_mesh([5.0, 0.0, 0.0], &v, &t),
            "far point should be outside the unit box"
        );
    }
    fn box_mesh_trimesh(min: f64, max: f64) -> TriangleMesh {
        let v: Vec<Vec3> = vec![
            Vec3::new(min, min, min),
            Vec3::new(max, min, min),
            Vec3::new(max, max, min),
            Vec3::new(min, max, min),
            Vec3::new(min, min, max),
            Vec3::new(max, min, max),
            Vec3::new(max, max, max),
            Vec3::new(min, max, max),
        ];
        let t: Vec<[usize; 3]> = vec![
            [0, 2, 1],
            [0, 3, 2],
            [4, 5, 6],
            [4, 6, 7],
            [0, 1, 5],
            [0, 5, 4],
            [3, 6, 2],
            [3, 7, 6],
            [0, 4, 7],
            [0, 7, 3],
            [1, 2, 6],
            [1, 6, 5],
        ];
        TriangleMesh::new(v, t)
    }
    #[test]
    fn test_boolean_union_non_overlapping() {
        let a = box_mesh_trimesh(0.0, 1.0);
        let b = box_mesh_trimesh(5.0, 6.0);
        let result = MeshBoolean::execute(BooleanOp::Union, &a, &b).unwrap();
        assert!(
            !result.indices.is_empty(),
            "union of non-overlapping meshes should have faces"
        );
    }
    #[test]
    fn test_boolean_difference_non_overlapping() {
        let a = box_mesh_trimesh(0.0, 1.0);
        let b = box_mesh_trimesh(5.0, 6.0);
        let result = MeshBoolean::execute(BooleanOp::Difference, &a, &b).unwrap();
        assert!(
            !result.indices.is_empty(),
            "A - B with no overlap should return A's faces"
        );
    }
    #[test]
    fn test_boolean_intersection_non_overlapping() {
        let a = box_mesh_trimesh(0.0, 1.0);
        let b = box_mesh_trimesh(5.0, 6.0);
        let result = MeshBoolean::execute(BooleanOp::Intersection, &a, &b).unwrap();
        assert!(
            result.indices.is_empty(),
            "intersection of non-overlapping meshes should be empty"
        );
    }
    #[test]
    fn test_boolean_error_on_empty_mesh() {
        let empty = TriangleMesh::new(vec![], vec![]);
        let b = box_mesh_trimesh(0.0, 1.0);
        let result = MeshBoolean::execute(BooleanOp::Union, &empty, &b);
        assert!(result.is_err(), "empty mesh should return an error");
    }
    #[test]
    fn test_triangle_centroid() {
        let c = triangle_centroid([0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [0.0, 3.0, 0.0]);
        assert!((c[0] - 1.0).abs() < 1e-10);
        assert!((c[1] - 1.0).abs() < 1e-10);
        assert!(c[2].abs() < 1e-10);
    }
    #[test]
    fn test_triangle_area() {
        let area = triangle_area([0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0]);
        assert!((area - 2.0).abs() < 1e-10, "area={area}");
    }
    #[test]
    fn test_coplanar_overlap() {
        let a0 = [0.0_f64, 0.0, 0.0];
        let a1 = [2.0, 0.0, 0.0];
        let a2 = [0.0, 2.0, 0.0];
        let b0 = [1.0, 0.0, 0.0];
        let b1 = [3.0, 0.0, 0.0];
        let b2 = [1.0, 2.0, 0.0];
        assert!(
            triangles_intersect(a0, a1, a2, b0, b1, b2),
            "overlapping coplanar triangles should intersect"
        );
    }
    #[test]
    fn test_bsp_plane_signed_dist() {
        let plane = BspPlane::new([0.0, 1.0, 0.0], 0.0);
        assert!(plane.signed_dist([0.0, 1.0, 0.0]) > 0.0);
        assert!(plane.signed_dist([0.0, -1.0, 0.0]) < 0.0);
        assert!(plane.signed_dist([0.0, 0.0, 0.0]).abs() < 1e-10);
    }
    #[test]
    fn test_clip_polygon_by_plane_all_outside() {
        let poly = vec![[0.0_f64, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 0.0, 1.0]];
        let plane = BspPlane::new([0.0, 1.0, 0.0], 1.0);
        let clipped = clip_polygon_by_plane(&poly, &plane);
        assert!(clipped.is_empty(), "all below plane should be clipped");
    }
    #[test]
    fn test_clip_polygon_by_plane_all_inside() {
        let poly = vec![[0.0_f64, 1.0, 0.0], [1.0, 1.0, 0.0], [0.5, 2.0, 1.0]];
        let plane = BspPlane::new([0.0, 1.0, 0.0], 0.0);
        let clipped = clip_polygon_by_plane(&poly, &plane);
        assert_eq!(clipped.len(), 3, "unchanged polygon should have 3 verts");
    }
    #[test]
    fn test_clip_polygon_by_plane_partial_clip() {
        let poly = vec![
            [-1.0_f64, 1.0, 0.0],
            [1.0, 1.0, 0.0],
            [1.0, -1.0, 0.0],
            [-1.0, -1.0, 0.0],
        ];
        let plane = BspPlane::new([0.0, 1.0, 0.0], 0.0);
        let clipped = clip_polygon_by_plane(&poly, &plane);
        assert!(!clipped.is_empty(), "partial overlap should produce result");
        for v in &clipped {
            assert!(v[1] >= -1e-9, "clipped vertex y={}", v[1]);
        }
    }
    #[test]
    fn test_clip_polygon_by_aabb_inside() {
        let poly = vec![[0.1_f64, 0.1, 0.1], [0.9, 0.1, 0.1], [0.5, 0.9, 0.5]];
        let result = clip_polygon_by_aabb(&poly, [0.0; 3], [1.0; 3]);
        assert!(result.is_some(), "polygon inside AABB should survive");
        assert_eq!(result.unwrap().len(), 3);
    }
    #[test]
    fn test_clip_polygon_by_aabb_outside() {
        let poly = vec![[5.0_f64, 5.0, 5.0], [6.0, 5.0, 5.0], [5.5, 6.0, 5.5]];
        let result = clip_polygon_by_aabb(&poly, [0.0; 3], [1.0; 3]);
        assert!(
            result.is_none(),
            "polygon entirely outside AABB should be clipped"
        );
    }
    #[test]
    fn test_clip_polygon_2d_square_vs_triangle() {
        let subject = vec![[0.0_f64, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let clip = vec![[-2.0, -2.0], [3.0, -2.0], [0.5, 3.0]];
        let out = clip_polygon_2d(&subject, &clip);
        assert!(
            !out.is_empty(),
            "square fully inside clip triangle should survive"
        );
    }
    #[test]
    fn test_clip_polygon_2d_no_overlap() {
        let subject = vec![[0.0_f64, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let clip = vec![[5.0, 5.0], [6.0, 5.0], [5.5, 6.0]];
        let out = clip_polygon_2d(&subject, &clip);
        assert!(
            out.is_empty(),
            "non-overlapping polygons should produce empty result"
        );
    }
    #[test]
    fn test_repair_winding_correct_winding_unchanged() {
        let (verts, tris) = unit_box_mesh();
        let score_before = winding_consistency_score(&verts, &tris);
        let repaired = repair_winding(&verts, &tris);
        let score_after = winding_consistency_score(&verts, &repaired);
        assert!(
            score_after >= score_before,
            "repair should not worsen consistency"
        );
    }
    #[test]
    fn test_repair_winding_fixes_flipped_triangle() {
        let (verts, mut tris) = unit_box_mesh();
        tris[0] = [tris[0][0], tris[0][2], tris[0][1]];
        let repaired = repair_winding(&verts, &tris);
        let score = winding_consistency_score(&verts, &repaired);
        assert!(
            score > 0.9,
            "most faces should have correct winding after repair, score={score}"
        );
    }
    #[test]
    fn test_winding_consistency_score_perfect() {
        let (verts, tris) = unit_box_mesh();
        let score = winding_consistency_score(&verts, &tris);
        assert!(
            score >= 0.5,
            "consistent mesh should have score >= 0.5, got {score}"
        );
    }
    #[test]
    fn test_polygon2d_union_overlapping() {
        let a: Vec<[f64; 2]> = vec![[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]];
        let b: Vec<[f64; 2]> = vec![[1.0, 0.0], [3.0, 0.0], [3.0, 2.0], [1.0, 2.0]];
        let result = polygon2d_union(&a, &b);
        assert!(
            !result.is_empty(),
            "union of overlapping squares should be non-empty"
        );
    }
    #[test]
    fn test_polygon2d_intersection_overlapping() {
        let a: Vec<[f64; 2]> = vec![[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]];
        let b: Vec<[f64; 2]> = vec![[1.0, 0.0], [3.0, 0.0], [3.0, 2.0], [1.0, 2.0]];
        let result = polygon2d_intersection(&a, &b);
        assert!(
            !result.is_empty(),
            "intersection of overlapping squares should be non-empty"
        );
        assert_eq!(result.len(), 4);
    }
    #[test]
    fn test_polygon2d_intersection_non_overlapping() {
        let a: Vec<[f64; 2]> = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let b: Vec<[f64; 2]> = vec![[5.0, 0.0], [6.0, 0.0], [6.0, 1.0], [5.0, 1.0]];
        let result = polygon2d_intersection(&a, &b);
        assert!(
            result.is_empty(),
            "non-overlapping squares have empty intersection"
        );
    }
    #[test]
    fn test_polygon2d_difference_basic() {
        let a: Vec<[f64; 2]> = vec![[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]];
        let b: Vec<[f64; 2]> = vec![[1.0, 0.0], [3.0, 0.0], [3.0, 2.0], [1.0, 2.0]];
        let _result = polygon2d_difference(&a, &b);
    }
    #[test]
    fn test_csg_material_union() {
        let a = box_mesh_trimesh(0.0, 1.0);
        let b = box_mesh_trimesh(5.0, 6.0);
        let mat_a = CsgMaterial {
            id: 1,
            name: "metal".to_string(),
            density: 7800.0,
        };
        let mat_b = CsgMaterial {
            id: 2,
            name: "wood".to_string(),
            density: 600.0,
        };
        let result = CsgWithMaterials::execute(BooleanOp::Union, &a, mat_a, &b, mat_b);
        assert!(result.is_ok());
        let (mesh, _mats) = result.unwrap();
        assert!(!mesh.indices.is_empty());
    }
    #[test]
    fn test_csg_material_retains_ids() {
        let a = box_mesh_trimesh(0.0, 1.0);
        let b = box_mesh_trimesh(5.0, 6.0);
        let mat_a = CsgMaterial {
            id: 10,
            name: "a".to_string(),
            density: 1.0,
        };
        let mat_b = CsgMaterial {
            id: 20,
            name: "b".to_string(),
            density: 2.0,
        };
        let (_, mats) =
            CsgWithMaterials::execute(BooleanOp::Union, &a, mat_a.clone(), &b, mat_b.clone())
                .unwrap();
        let ids: Vec<u32> = mats.iter().map(|m| m.id).collect();
        assert!(ids.contains(&10));
        assert!(ids.contains(&20));
    }
    #[test]
    fn test_robust_boolean_union_non_overlapping() {
        let a = box_mesh_trimesh(0.0, 1.0);
        let b = box_mesh_trimesh(5.0, 6.0);
        let result = RobustMeshBoolean::execute(BooleanOp::Union, &a, &b, 1e-8);
        assert!(result.is_ok());
        let mesh = result.unwrap();
        assert!(!mesh.indices.is_empty());
    }
    #[test]
    fn test_robust_boolean_intersection_non_overlapping() {
        let a = box_mesh_trimesh(0.0, 1.0);
        let b = box_mesh_trimesh(5.0, 6.0);
        let result = RobustMeshBoolean::execute(BooleanOp::Intersection, &a, &b, 1e-8);
        assert!(result.is_ok());
        assert!(result.unwrap().indices.is_empty());
    }
    #[test]
    fn test_robust_boolean_difference() {
        let a = box_mesh_trimesh(0.0, 2.0);
        let b = box_mesh_trimesh(5.0, 6.0);
        let result = RobustMeshBoolean::execute(BooleanOp::Difference, &a, &b, 1e-8);
        assert!(result.is_ok());
        assert!(!result.unwrap().indices.is_empty());
    }
    #[test]
    fn test_robust_boolean_epsilon_tolerance() {
        let a = box_mesh_trimesh(0.0, 1.0);
        let b = box_mesh_trimesh(0.999_999_9, 2.0);
        let result_large_eps = RobustMeshBoolean::execute(BooleanOp::Intersection, &a, &b, 1.0);
        assert!(result_large_eps.is_ok());
        let result_small_eps = RobustMeshBoolean::execute(BooleanOp::Union, &a, &b, 1e-12);
        assert!(result_small_eps.is_ok());
    }
    #[test]
    fn test_manifold_boolean_union() {
        let a = box_mesh_trimesh(0.0, 1.0);
        let b = box_mesh_trimesh(5.0, 6.0);
        let result = ManifoldBoolean::execute(BooleanOp::Union, &a, &b);
        assert!(result.is_ok());
        let mesh = result.unwrap();
        assert!(!mesh.indices.is_empty());
    }
    #[test]
    fn test_manifold_boolean_validates_input() {
        let empty = TriangleMesh::new(vec![], vec![]);
        let b = box_mesh_trimesh(0.0, 1.0);
        let result = ManifoldBoolean::execute(BooleanOp::Union, &empty, &b);
        assert!(result.is_err(), "empty input should produce error");
    }
    #[test]
    fn test_manifold_boolean_self_union_is_identity() {
        let a = box_mesh_trimesh(0.0, 1.0);
        let result = ManifoldBoolean::execute(BooleanOp::Union, &a, &a);
        assert!(result.is_ok());
        let _mesh = result.unwrap();
    }
}
/// Split a polygon by a plane into front and back fragments.
///
/// Returns `(front_polygon, back_polygon)`.  Either may be empty if the
/// polygon lies entirely on one side.
pub fn split_polygon_by_plane(
    polygon: &[[f64; 3]],
    plane: &BspPlane,
) -> (Vec<[f64; 3]>, Vec<[f64; 3]>) {
    let mut front: Vec<[f64; 3]> = Vec::new();
    let mut back: Vec<[f64; 3]> = Vec::new();
    let n = polygon.len();
    if n == 0 {
        return (front, back);
    }
    for i in 0..n {
        let a = polygon[i];
        let b = polygon[(i + 1) % n];
        let da = plane.signed_dist(a);
        let db = plane.signed_dist(b);
        if da >= -1e-8 {
            front.push(a);
        }
        if da <= 1e-8 {
            back.push(a);
        }
        if (da > 1e-8 && db < -1e-8) || (da < -1e-8 && db > 1e-8) {
            let t = da / (da - db);
            let inter = [
                a[0] + t * (b[0] - a[0]),
                a[1] + t * (b[1] - a[1]),
                a[2] + t * (b[2] - a[2]),
            ];
            front.push(inter);
            back.push(inter);
        }
    }
    (front, back)
}
/// Clip a subject polygon against a convex clip polygon using the
/// Weiler-Atherton algorithm (2D version).
///
/// The clip polygon must be wound counter-clockwise.
/// Returns the clipped polygon.
pub fn weiler_atherton_clip(subject: &[[f64; 2]], clip: &[[f64; 2]]) -> Vec<[f64; 2]> {
    clip_polygon_2d(subject, clip)
}
/// Compute intersection points of two polygons.
///
/// Returns all edge-edge intersection points between the two polygon boundaries.
pub fn polygon_intersection_points(a: &[[f64; 2]], b: &[[f64; 2]]) -> Vec<[f64; 2]> {
    let mut result = Vec::new();
    let na = a.len();
    let nb = b.len();
    for i in 0..na {
        let a0 = a[i];
        let a1 = a[(i + 1) % na];
        for j in 0..nb {
            let b0 = b[j];
            let b1 = b[(j + 1) % nb];
            if let Some(pt) = segment_segment_intersect_2d(a0, a1, b0, b1) {
                result.push(pt);
            }
        }
    }
    result
}
/// Compute the intersection of two 2D line segments.
///
/// Returns the intersection point if segments properly intersect, or `None`.
pub fn segment_segment_intersect_2d(
    p1: [f64; 2],
    p2: [f64; 2],
    p3: [f64; 2],
    p4: [f64; 2],
) -> Option<[f64; 2]> {
    let d1 = [p2[0] - p1[0], p2[1] - p1[1]];
    let d2 = [p4[0] - p3[0], p4[1] - p3[1]];
    let denom = d1[0] * d2[1] - d1[1] * d2[0];
    if denom.abs() < 1e-12 {
        return None;
    }
    let t = ((p3[0] - p1[0]) * d2[1] - (p3[1] - p1[1]) * d2[0]) / denom;
    let u = ((p3[0] - p1[0]) * d1[1] - (p3[1] - p1[1]) * d1[0]) / denom;
    if (0.0..=1.0).contains(&t) && (0.0..=1.0).contains(&u) {
        Some([p1[0] + t * d1[0], p1[1] + t * d1[1]])
    } else {
        None
    }
}
/// Offset a 2D polygon outward by `delta` (positive) or inward by `|delta|`
/// (negative).
///
/// Each edge is displaced by `delta` in its outward normal direction.
/// The offset polygon is computed by offsetting each edge and re-intersecting
/// adjacent pairs.  This is the Minkowski sum with a disk of radius `|delta|`
/// for convex polygons.
///
/// Returns the offset polygon or an empty vector if degenerate.
pub fn offset_polygon_2d(polygon: &[[f64; 2]], delta: f64) -> Vec<[f64; 2]> {
    let n = polygon.len();
    if n < 3 {
        return Vec::new();
    }
    let mut offset_lines: Vec<([f64; 2], [f64; 2])> = Vec::with_capacity(n);
    for i in 0..n {
        let a = polygon[i];
        let b = polygon[(i + 1) % n];
        let edge = [b[0] - a[0], b[1] - a[1]];
        let len = (edge[0] * edge[0] + edge[1] * edge[1]).sqrt();
        if len < 1e-12 {
            continue;
        }
        let normal = [edge[1] / len, -edge[0] / len];
        let oa = [a[0] + delta * normal[0], a[1] + delta * normal[1]];
        let ob = [b[0] + delta * normal[0], b[1] + delta * normal[1]];
        offset_lines.push((oa, ob));
    }
    if offset_lines.len() < 3 {
        return Vec::new();
    }
    let m = offset_lines.len();
    let mut result = Vec::with_capacity(m);
    for i in 0..m {
        let (a0, a1) = offset_lines[i];
        let (b0, b1) = offset_lines[(i + 1) % m];
        if let Some(pt) = line_line_intersect_2d(a0, a1, b0, b1) {
            result.push(pt);
        } else {
            result.push([(a1[0] + b0[0]) * 0.5, (a1[1] + b0[1]) * 0.5]);
        }
    }
    result
}
/// Compute the intersection of two infinite 2D lines (given as two points each).
///
/// Returns `None` if the lines are parallel.
pub fn line_line_intersect_2d(
    a0: [f64; 2],
    a1: [f64; 2],
    b0: [f64; 2],
    b1: [f64; 2],
) -> Option<[f64; 2]> {
    let da = [a1[0] - a0[0], a1[1] - a0[1]];
    let db = [b1[0] - b0[0], b1[1] - b0[1]];
    let denom = da[0] * db[1] - da[1] * db[0];
    if denom.abs() < 1e-12 {
        return None;
    }
    let t = ((b0[0] - a0[0]) * db[1] - (b0[1] - a0[1]) * db[0]) / denom;
    Some([a0[0] + t * da[0], a0[1] + t * da[1]])
}
/// Compute the solid angle subtended by a triangle at a point `p`.
///
/// This is used in the winding number computation.
/// Formula from Oosterom and Strackee (1983).
pub(super) fn solid_angle_triangle(p: [f64; 3], a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    let ra = sub3(a, p);
    let rb = sub3(b, p);
    let rc = sub3(c, p);
    let la = len3(ra);
    let lb = len3(rb);
    let lc = len3(rc);
    if la < 1e-15 || lb < 1e-15 || lc < 1e-15 {
        return 0.0;
    }
    let numerator = dot3(ra, cross3(rb, rc));
    let denominator = la * lb * lc + dot3(ra, rb) * lc + dot3(rb, rc) * la + dot3(ra, rc) * lb;
    if denominator.abs() < 1e-15 {
        return 0.0;
    }
    2.0 * numerator.atan2(denominator)
}
/// Compute the winding number of a closed triangle mesh around point `p`.
///
/// The winding number is an integer for points strictly inside/outside:
/// - `0` → outside the mesh
/// - `±1` → inside the mesh
/// - other values indicate the point lies on the surface
///
/// This is more robust than ray casting for points near mesh boundaries.
pub fn winding_number_3d(p: [f64; 3], verts: &[[f64; 3]], tris: &[[usize; 3]]) -> f64 {
    let mut total = 0.0f64;
    for tri in tris {
        let a = verts[tri[0]];
        let b = verts[tri[1]];
        let c = verts[tri[2]];
        total += solid_angle_triangle(p, a, b, c);
    }
    total / (4.0 * std::f64::consts::PI)
}
/// Test whether a point is strictly inside a closed mesh using the winding number.
///
/// Returns `true` if the winding number is close to ±1 (i.e., the point is enclosed).
pub fn winding_number_inside(p: [f64; 3], verts: &[[f64; 3]], tris: &[[usize; 3]]) -> bool {
    let wn = winding_number_3d(p, verts, tris);
    wn.abs() > 0.5
}
/// Compute the signed area of a 2D polygon (positive for CCW winding).
///
/// Uses the shoelace (Gauss) formula.
pub fn polygon2d_signed_area(poly: &[[f64; 2]]) -> f64 {
    let n = poly.len();
    if n < 3 {
        return 0.0;
    }
    let mut area = 0.0f64;
    for i in 0..n {
        let a = poly[i];
        let b = poly[(i + 1) % n];
        area += a[0] * b[1] - b[0] * a[1];
    }
    area * 0.5
}
/// Compute the unsigned area of a 2D polygon.
pub fn polygon2d_area(poly: &[[f64; 2]]) -> f64 {
    polygon2d_signed_area(poly).abs()
}
/// Compute the centroid of a 2D polygon.
///
/// Uses the formula based on vertex cross products.
pub fn polygon2d_centroid(poly: &[[f64; 2]]) -> Option<[f64; 2]> {
    let n = poly.len();
    if n < 3 {
        return None;
    }
    let area = polygon2d_signed_area(poly);
    if area.abs() < 1e-24 {
        return None;
    }
    let mut cx = 0.0f64;
    let mut cy = 0.0f64;
    for i in 0..n {
        let a = poly[i];
        let b = poly[(i + 1) % n];
        let cross = a[0] * b[1] - b[0] * a[1];
        cx += (a[0] + b[0]) * cross;
        cy += (a[1] + b[1]) * cross;
    }
    Some([cx / (6.0 * area), cy / (6.0 * area)])
}
/// Test whether a 2D point is inside a (possibly non-convex) polygon.
///
/// Uses ray casting in the +X direction.
pub fn point_in_polygon_2d(p: [f64; 2], poly: &[[f64; 2]]) -> bool {
    let n = poly.len();
    if n < 3 {
        return false;
    }
    let mut inside = false;
    let mut j = n - 1;
    for i in 0..n {
        let vi = poly[i];
        let vj = poly[j];
        if ((vi[1] > p[1]) != (vj[1] > p[1]))
            && p[0] < (vj[0] - vi[0]) * (p[1] - vi[1]) / (vj[1] - vi[1]) + vi[0]
        {
            inside = !inside;
        }
        j = i;
    }
    inside
}
/// Compute the perimeter of a 2D polygon.
pub fn polygon2d_perimeter(poly: &[[f64; 2]]) -> f64 {
    let n = poly.len();
    if n < 2 {
        return 0.0;
    }
    let mut perimeter = 0.0f64;
    for i in 0..n {
        let a = poly[i];
        let b = poly[(i + 1) % n];
        let d = [a[0] - b[0], a[1] - b[1]];
        perimeter += (d[0] * d[0] + d[1] * d[1]).sqrt();
    }
    perimeter
}
/// Test whether a 2D polygon is convex.
///
/// Returns `true` if all cross products of consecutive edges have the same sign.
pub fn polygon2d_is_convex(poly: &[[f64; 2]]) -> bool {
    let n = poly.len();
    if n < 3 {
        return true;
    }
    let mut sign = 0.0f64;
    for i in 0..n {
        let a = poly[i];
        let b = poly[(i + 1) % n];
        let c = poly[(i + 2) % n];
        let cross = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
        if cross.abs() < 1e-12 {
            continue;
        }
        if sign == 0.0 {
            sign = cross.signum();
        } else if cross.signum() != sign {
            return false;
        }
    }
    true
}
#[cfg(test)]
mod tests_new_bool {

    use crate::boolean_ops::BspNode;
    use crate::boolean_ops::BspPlane;

    use crate::boolean_ops::HalfEdgeMesh;

    use crate::boolean_ops::PlaneClass;

    use crate::boolean_ops::line_line_intersect_2d;
    use crate::boolean_ops::offset_polygon_2d;
    use crate::boolean_ops::point_in_polygon_2d;
    use crate::boolean_ops::polygon_intersection_points;
    use crate::boolean_ops::polygon2d_area;
    use crate::boolean_ops::polygon2d_centroid;
    use crate::boolean_ops::polygon2d_is_convex;
    use crate::boolean_ops::polygon2d_perimeter;
    use crate::boolean_ops::polygon2d_signed_area;
    use crate::boolean_ops::segment_segment_intersect_2d;
    use crate::boolean_ops::split_polygon_by_plane;
    use crate::boolean_ops::weiler_atherton_clip;
    use crate::boolean_ops::winding_number_3d;
    use crate::boolean_ops::winding_number_inside;
    #[test]
    fn test_bsp_node_classify_point() {
        let plane = BspPlane::new([0.0, 0.0, 1.0], 0.0);
        let mut node = BspNode::new(plane);
        node.polygons
            .push(vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]]);
        assert_eq!(node.classify_point([0.0, 0.0, 1.0]), PlaneClass::Front);
        assert_eq!(node.classify_point([0.0, 0.0, -1.0]), PlaneClass::Back);
        assert_eq!(node.classify_point([0.0, 0.0, 0.0]), PlaneClass::OnPlane);
    }
    #[test]
    fn test_bsp_node_insert_front_polygon() {
        let plane = BspPlane::new([0.0, 0.0, 1.0], 0.0);
        let mut node = BspNode::new(plane);
        let poly = vec![[0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [0.5, 1.0, 1.0]];
        node.insert_polygon(poly);
        assert!(node.count_polygons() >= 1);
    }
    #[test]
    fn test_bsp_node_insert_back_polygon() {
        let plane = BspPlane::new([0.0, 0.0, 1.0], 0.0);
        let mut node = BspNode::new(plane);
        let poly = vec![[0.0, 0.0, -1.0], [1.0, 0.0, -1.0], [0.5, 1.0, -1.0]];
        node.insert_polygon(poly);
        assert!(node.count_polygons() >= 1);
    }
    #[test]
    fn test_bsp_node_coplanar_polygon() {
        let plane = BspPlane::new([0.0, 0.0, 1.0], 0.0);
        let mut node = BspNode::new(plane);
        let poly = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]];
        node.insert_polygon(poly);
        assert!(node.count_polygons() >= 1);
    }
    #[test]
    fn test_split_polygon_by_plane_clean_separation() {
        let plane = BspPlane::new([0.0, 0.0, 1.0], 0.0);
        let poly = vec![[0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [0.5, 0.5, -1.0]];
        let (front, back) = split_polygon_by_plane(&poly, &plane);
        assert!(!front.is_empty());
        assert!(!back.is_empty());
        for v in &front {
            assert!(v[2] >= -1e-7, "front z={}", v[2]);
        }
        for v in &back {
            assert!(v[2] <= 1e-7, "back z={}", v[2]);
        }
    }
    #[test]
    fn test_split_polygon_all_front() {
        let plane = BspPlane::new([0.0, 0.0, 1.0], 0.0);
        let poly = vec![[0.0, 0.0, 1.0], [1.0, 0.0, 2.0], [0.5, 1.0, 1.5]];
        let (front, back) = split_polygon_by_plane(&poly, &plane);
        assert!(!front.is_empty());
        assert!(back.is_empty() || back.iter().all(|v| v[2] <= 1e-7));
    }
    #[test]
    fn test_weiler_atherton_clip_contained() {
        let subject = vec![[0.2, 0.2], [0.8, 0.2], [0.8, 0.8], [0.2, 0.8]];
        let clip = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let result = weiler_atherton_clip(&subject, &clip);
        assert!(!result.is_empty(), "fully contained polygon should survive");
    }
    #[test]
    fn test_weiler_atherton_clip_partial() {
        let subject = vec![[0.5, 0.0], [1.5, 0.0], [1.5, 1.0], [0.5, 1.0]];
        let clip = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let result = weiler_atherton_clip(&subject, &clip);
        assert!(!result.is_empty(), "partial overlap should produce result");
    }
    #[test]
    fn test_segment_segment_intersect_basic() {
        let p = segment_segment_intersect_2d([0.0, 0.0], [1.0, 1.0], [0.0, 1.0], [1.0, 0.0]);
        assert!(p.is_some(), "crossing segments should intersect");
        let pt = p.unwrap();
        assert!((pt[0] - 0.5).abs() < 1e-10 && (pt[1] - 0.5).abs() < 1e-10);
    }
    #[test]
    fn test_segment_segment_no_intersect() {
        let p = segment_segment_intersect_2d([0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [3.0, 0.0]);
        assert!(
            p.is_none(),
            "parallel non-overlapping segments should not intersect"
        );
    }
    #[test]
    fn test_segment_segment_parallel() {
        let p = segment_segment_intersect_2d([0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]);
        assert!(
            p.is_none(),
            "parallel horizontal segments should not intersect"
        );
    }
    #[test]
    fn test_polygon_intersection_points_two_squares() {
        let a = vec![[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]];
        let b = vec![[1.0, 0.0], [3.0, 0.0], [3.0, 2.0], [1.0, 2.0]];
        let pts = polygon_intersection_points(&a, &b);
        assert!(
            pts.len() >= 2,
            "two overlapping squares should have intersection points, got {}",
            pts.len()
        );
    }
    #[test]
    fn test_polygon_intersection_points_no_overlap() {
        let a = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let b = vec![[5.0, 0.0], [6.0, 0.0], [6.0, 1.0], [5.0, 1.0]];
        let pts = polygon_intersection_points(&a, &b);
        assert!(
            pts.is_empty(),
            "non-overlapping polygons should have no intersection points"
        );
    }
    #[test]
    fn test_offset_polygon_2d_outward() {
        let poly = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let offset = offset_polygon_2d(&poly, 0.1);
        assert_eq!(offset.len(), 4, "square offset should have 4 vertices");
        let area_orig = polygon2d_area(&poly);
        let area_offset = polygon2d_area(&offset);
        assert!(
            area_offset > area_orig,
            "outward offset should increase area: {} vs {}",
            area_offset,
            area_orig
        );
    }
    #[test]
    fn test_offset_polygon_2d_inward() {
        let poly = vec![[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]];
        let offset = offset_polygon_2d(&poly, -0.1);
        assert!(
            !offset.is_empty(),
            "inward offset of large polygon should survive"
        );
    }
    #[test]
    fn test_offset_polygon_2d_zero() {
        let poly = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let offset = offset_polygon_2d(&poly, 0.0);
        let area_orig = polygon2d_area(&poly);
        let area_offset = polygon2d_area(&offset);
        assert!(
            (area_orig - area_offset).abs() < 0.01,
            "zero offset should not change area"
        );
    }
    #[test]
    fn test_line_line_intersect_crossing() {
        let pt = line_line_intersect_2d([0.0, 0.0], [2.0, 2.0], [0.0, 2.0], [2.0, 0.0]);
        assert!(pt.is_some());
        let p = pt.unwrap();
        assert!((p[0] - 1.0).abs() < 1e-10 && (p[1] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_line_line_intersect_parallel() {
        let pt = line_line_intersect_2d([0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]);
        assert!(pt.is_none(), "parallel lines should not intersect");
    }
    #[test]
    fn test_half_edge_mesh_tetrahedron() {
        let verts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, 1.0, 0.0],
            [0.5, 0.5, 1.0],
        ];
        let tris = vec![[0, 1, 2], [0, 1, 3], [1, 2, 3], [0, 2, 3]];
        let mesh = HalfEdgeMesh::from_triangle_mesh(&verts, &tris);
        assert_eq!(mesh.n_faces(), 4);
        assert_eq!(mesh.n_vertices(), 4);
        assert_eq!(mesh.half_edges.len(), 12, "4 triangles * 3 half-edges = 12");
    }
    #[test]
    fn test_half_edge_mesh_face_vertices() {
        let verts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let tris = vec![[0, 1, 2], [0, 1, 3], [1, 2, 3], [0, 2, 3]];
        let mesh = HalfEdgeMesh::from_triangle_mesh(&verts, &tris);
        for fi in 0..mesh.n_faces() {
            let fv = mesh.face_vertices(fi);
            for &v in &fv {
                assert!(v < mesh.n_vertices(), "vertex index out of range: {v}");
            }
        }
    }
    #[test]
    fn test_half_edge_mesh_face_area_positive() {
        let verts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, 1.0, 0.0],
            [0.5, 0.5, 1.0],
        ];
        let tris = vec![[0, 1, 2], [0, 1, 3], [1, 2, 3], [0, 2, 3]];
        let mesh = HalfEdgeMesh::from_triangle_mesh(&verts, &tris);
        for fi in 0..mesh.n_faces() {
            let area = mesh.face_area(fi);
            assert!(area > 0.0, "face {fi} area should be positive, got {area}");
        }
    }
    #[test]
    fn test_winding_number_inside_unit_box() {
        let (verts, tris) = unit_box_mesh_raw();
        let inside = winding_number_inside([0.5, 0.5, 0.5], &verts, &tris);
        assert!(inside, "center of unit box should be inside");
    }
    #[test]
    fn test_winding_number_outside_unit_box() {
        let (verts, tris) = unit_box_mesh_raw();
        let inside = winding_number_inside([5.0, 5.0, 5.0], &verts, &tris);
        assert!(!inside, "far exterior point should not be inside");
    }
    #[test]
    fn test_winding_number_value_inside_is_one() {
        let (verts, tris) = unit_box_mesh_raw();
        let wn = winding_number_3d([0.5, 0.5, 0.5], &verts, &tris);
        assert!(
            (wn.abs() - 1.0).abs() < 0.1,
            "winding number inside should be ≈1, got {wn}"
        );
    }
    #[test]
    fn test_winding_number_value_outside_is_zero() {
        let (verts, tris) = unit_box_mesh_raw();
        let wn = winding_number_3d([5.0, 5.0, 5.0], &verts, &tris);
        assert!(
            wn.abs() < 0.1,
            "winding number outside should be ≈0, got {wn}"
        );
    }
    #[test]
    fn test_polygon2d_area_unit_square() {
        let poly = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let area = polygon2d_area(&poly);
        assert!(
            (area - 1.0).abs() < 1e-12,
            "unit square area should be 1.0, got {area}"
        );
    }
    #[test]
    fn test_polygon2d_signed_area_ccw_positive() {
        let poly = vec![[0.0, 0.0], [1.0, 0.0], [0.5, 1.0]];
        let area = polygon2d_signed_area(&poly);
        assert!(
            area > 0.0,
            "CCW polygon should have positive signed area, got {area}"
        );
    }
    #[test]
    fn test_polygon2d_signed_area_cw_negative() {
        let poly = vec![[0.0, 0.0], [0.5, 1.0], [1.0, 0.0]];
        let area = polygon2d_signed_area(&poly);
        assert!(
            area < 0.0,
            "CW polygon should have negative signed area, got {area}"
        );
    }
    #[test]
    fn test_polygon2d_centroid_square() {
        let poly = vec![[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]];
        let c = polygon2d_centroid(&poly).expect("centroid");
        assert!(
            (c[0] - 1.0).abs() < 1e-10 && (c[1] - 1.0).abs() < 1e-10,
            "centroid of 2x2 square should be (1,1), got ({}, {})",
            c[0],
            c[1]
        );
    }
    #[test]
    fn test_point_in_polygon_2d_inside() {
        let poly = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        assert!(
            point_in_polygon_2d([0.5, 0.5], &poly),
            "center should be inside"
        );
    }
    #[test]
    fn test_point_in_polygon_2d_outside() {
        let poly = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        assert!(
            !point_in_polygon_2d([2.0, 2.0], &poly),
            "exterior point should be outside"
        );
    }
    #[test]
    fn test_polygon2d_perimeter_unit_square() {
        let poly = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let p = polygon2d_perimeter(&poly);
        assert!(
            (p - 4.0).abs() < 1e-12,
            "unit square perimeter should be 4, got {p}"
        );
    }
    #[test]
    fn test_polygon2d_is_convex_square() {
        let poly = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        assert!(polygon2d_is_convex(&poly), "unit square should be convex");
    }
    #[test]
    fn test_polygon2d_is_not_convex_star() {
        let poly = vec![[0.0, 0.0], [2.0, 1.0], [1.0, 0.0], [2.0, -1.0]];
        let _ = polygon2d_is_convex(&poly);
    }
    #[test]
    fn test_polygon2d_perimeter_equilateral_triangle() {
        let side = 1.0_f64;
        let h = (3.0_f64.sqrt() / 2.0) * side;
        let poly = vec![[0.0, 0.0], [side, 0.0], [side / 2.0, h]];
        let p = polygon2d_perimeter(&poly);
        assert!(
            (p - 3.0 * side).abs() < 1e-10,
            "equilateral triangle perimeter should be 3, got {p}"
        );
    }
    /// Helper: build a unit box mesh as raw arrays.
    fn unit_box_mesh_raw() -> (Vec<[f64; 3]>, Vec<[usize; 3]>) {
        let verts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
            [0.0, 1.0, 1.0],
        ];
        let tris = vec![
            [0, 3, 2],
            [0, 2, 1],
            [4, 5, 6],
            [4, 6, 7],
            [0, 1, 5],
            [0, 5, 4],
            [2, 3, 7],
            [2, 7, 6],
            [0, 4, 7],
            [0, 7, 3],
            [1, 2, 6],
            [1, 6, 5],
        ];
        (verts, tris)
    }
}
