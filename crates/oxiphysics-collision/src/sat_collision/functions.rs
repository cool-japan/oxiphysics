//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{Obb, ObbBvh, ObbCollision};

pub(super) fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
pub(super) fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
pub(super) fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
pub(super) fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
pub(super) fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
pub(super) fn length(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}
pub(super) fn normalize(a: [f64; 3]) -> [f64; 3] {
    let len = length(a);
    if len < 1e-12 {
        [0.0, 0.0, 0.0]
    } else {
        scale(a, 1.0 / len)
    }
}
pub(super) fn rotate(axes: [[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [
        axes[0][0] * v[0] + axes[1][0] * v[1] + axes[2][0] * v[2],
        axes[0][1] * v[0] + axes[1][1] * v[1] + axes[2][1] * v[2],
        axes[0][2] * v[0] + axes[1][2] * v[1] + axes[2][2] * v[2],
    ]
}
#[cfg(test)]
mod tests {
    use super::*;

    use crate::ContactManifoldBuilder;

    use crate::ConvexPolytope;
    use crate::Obb;

    use crate::ObbCollision;
    use crate::ObbTree;

    use crate::sat_collision::length;
    use crate::sat_collision::sub;
    #[test]
    fn test_obb_from_aabb_axes_and_half_extents() {
        let obb = Obb::from_aabb([0.0, 0.0, 0.0], [2.0, 4.0, 6.0]);
        assert_eq!(
            obb.axes,
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        );
        assert!((obb.half_extents[0] - 1.0).abs() < 1e-12);
        assert!((obb.half_extents[1] - 2.0).abs() < 1e-12);
        assert!((obb.half_extents[2] - 3.0).abs() < 1e-12);
        assert!((obb.center[0] - 1.0).abs() < 1e-12);
        assert!((obb.center[1] - 2.0).abs() < 1e-12);
        assert!((obb.center[2] - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_obb_vertices_eight_distinct() {
        let obb = Obb::from_aabb([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]);
        let verts = obb.vertices();
        for i in 0..8 {
            for j in (i + 1)..8 {
                let d = length(sub(verts[i], verts[j]));
                assert!(d > 1e-9, "vertices {} and {} are the same", i, j);
            }
        }
    }
    #[test]
    fn test_obb_obb_sat_overlapping() {
        let a = Obb::from_aabb([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]);
        let b = Obb::from_aabb([-0.5, -0.5, -0.5], [1.5, 1.5, 1.5]);
        let result = ObbCollision::obb_obb_sat(&a, &b);
        assert!(result.is_some(), "overlapping OBBs should return Some");
    }
    #[test]
    fn test_obb_obb_sat_separated() {
        let a = Obb::from_aabb([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]);
        let b = Obb::from_aabb([3.0, 0.0, 0.0], [5.0, 2.0, 2.0]);
        let result = ObbCollision::obb_obb_sat(&a, &b);
        assert!(result.is_none(), "separated OBBs should return None");
    }
    #[test]
    fn test_test_overlap_on_axis_overlapping() {
        let a = Obb::from_aabb([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]);
        let b = Obb::from_aabb([0.5, -1.0, -1.0], [2.5, 1.0, 1.0]);
        let axis = [1.0, 0.0, 0.0];
        let overlap = ObbCollision::test_overlap_on_axis(&a, &b, axis);
        assert!(overlap.is_some(), "should report overlap");
        assert!(
            (overlap.unwrap() - 0.5).abs() < 1e-9,
            "overlap should be 0.5"
        );
    }
    #[test]
    fn test_clip_polygon_by_plane_square() {
        let square = [
            [-1.0f64, -1.0, 0.0],
            [1.0, -1.0, 0.0],
            [1.0, 1.0, 0.0],
            [-1.0, 1.0, 0.0],
        ];
        let clipped = ContactManifoldBuilder::clip_polygon_by_plane(&square, [1.0, 0.0, 0.0], 0.0);
        assert!(!clipped.is_empty(), "clipped polygon should not be empty");
        for p in &clipped {
            assert!(
                p[0] >= -1e-9,
                "all clipped points should have x >= 0, got {}",
                p[0]
            );
        }
    }
    #[test]
    fn test_obb_tree_build_single_triangle() {
        let vertices = [[0.0f64, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]];
        let triangles = [[0usize, 1, 2]];
        let tree = ObbTree::build(&vertices, &triangles);
        assert!(tree.obb.half_extents[0] >= 0.0);
        assert!(tree.obb.half_extents[1] >= 0.0);
        assert!(tree.obb.half_extents[2] >= 0.0);
    }
    #[test]
    fn test_convex_polytope_support() {
        let obb = Obb::from_aabb([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]);
        let poly = ConvexPolytope::from_obb(&obb);
        let identity_transform = (
            [0.0, 0.0, 0.0],
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        );
        let dir = [1.0, 0.0, 0.0];
        let sp = poly.support(dir, identity_transform);
        assert!(
            (sp[0] - 1.0).abs() < 1e-9,
            "support x should be 1.0, got {}",
            sp[0]
        );
        let best = poly
            .vertices
            .iter()
            .map(|&v| dot(v, dir))
            .fold(f64::NEG_INFINITY, f64::max);
        assert!((dot(sp, dir) - best).abs() < 1e-9);
    }
}
/// Compute the closest points between two line segments.
///
/// Returns `(t_a, t_b, point_a, point_b)` where `t` is the parametric distance
/// along each segment.
pub fn closest_points_on_segments(
    p1: [f64; 3],
    d1: [f64; 3],
    p2: [f64; 3],
    d2: [f64; 3],
) -> (f64, f64, [f64; 3], [f64; 3]) {
    let r = sub(p1, p2);
    let a = dot(d1, d1);
    let e = dot(d2, d2);
    let f = dot(d2, r);
    let (t, s);
    if a < 1e-12 && e < 1e-12 {
        t = 0.0;
        s = 0.0;
    } else if a < 1e-12 {
        s = (f / e).clamp(0.0, 1.0);
        t = 0.0;
    } else {
        let c = dot(d1, r);
        if e < 1e-12 {
            t = (-c / a).clamp(0.0, 1.0);
            s = 0.0;
        } else {
            let b = dot(d1, d2);
            let denom = a * e - b * b;
            t = if denom.abs() > 1e-12 {
                ((b * f - c * e) / denom).clamp(0.0, 1.0)
            } else {
                0.0
            };
            s = ((b * t + f) / e).clamp(0.0, 1.0);
        }
    }
    let pt_a = add(p1, scale(d1, t));
    let pt_b = add(p2, scale(d2, s));
    (t, s, pt_a, pt_b)
}
/// Generate edge-edge contact points from a SAT overlap.
///
/// Given the two closest edges (as segment start points and direction vectors),
/// returns the contact point and normal.
pub fn edge_edge_contact(
    p_a: [f64; 3],
    d_a: [f64; 3],
    p_b: [f64; 3],
    d_b: [f64; 3],
) -> ([f64; 3], [f64; 3]) {
    let (_, _, pt_a, pt_b) = closest_points_on_segments(p_a, d_a, p_b, d_b);
    let diff = sub(pt_a, pt_b);
    let len = length(diff);
    let normal = if len > 1e-10 {
        scale(diff, 1.0 / len)
    } else {
        [0.0, 1.0, 0.0]
    };
    let contact_point = scale(add(pt_a, pt_b), 0.5);
    (contact_point, normal)
}
#[cfg(test)]
mod tests_extended {
    use super::*;
    use crate::Capsule;

    use crate::Obb;

    use crate::ObbCapsuleCollision;
    use crate::ObbCollision;

    use crate::ObbTriangleCollision;

    use crate::SatContactPointGenerator;
    use crate::closest_points_on_segments;

    use crate::sat_collision::length;
    use crate::sat_collision::sub;
    #[test]
    fn capsule_axis() {
        let cap = Capsule::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.5);
        let ax = cap.axis();
        assert!((ax[0] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn capsule_closest_axis_point_clamped() {
        let cap = Capsule::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.5);
        let cp = cap.closest_axis_point([3.0, 0.0, 0.0]);
        assert!((cp[0] - 1.0).abs() < 1e-12, "x={}", cp[0]);
    }
    #[test]
    fn capsule_closest_axis_point_inside() {
        let cap = Capsule::new([0.0, 0.0, 0.0], [2.0, 0.0, 0.0], 0.5);
        let cp = cap.closest_axis_point([1.0, 1.0, 0.0]);
        assert!((cp[0] - 1.0).abs() < 1e-10, "x={}", cp[0]);
        assert!(cp[1].abs() < 1e-12, "y={}", cp[1]);
    }
    #[test]
    fn capsule_support_along_axis() {
        let cap = Capsule::new([0.0, 0.0, 0.0], [2.0, 0.0, 0.0], 1.0);
        let sp = cap.support([1.0, 0.0, 0.0]);
        assert!((sp[0] - 3.0).abs() < 1e-10, "support x={}", sp[0]);
    }
    #[test]
    fn obb_capsule_intersecting() {
        let obb = Obb::from_aabb([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]);
        let cap = Capsule::new([0.0, 0.0, -3.0], [0.0, 0.0, 3.0], 0.5);
        let result = ObbCapsuleCollision::test(&obb, &cap);
        assert!(result.is_some(), "capsule through OBB should intersect");
    }
    #[test]
    fn obb_capsule_separated() {
        let obb = Obb::from_aabb([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]);
        let cap = Capsule::new([5.0, 0.0, 0.0], [7.0, 0.0, 0.0], 0.3);
        let result = ObbCapsuleCollision::test(&obb, &cap);
        assert!(
            result.is_none(),
            "separated OBB and capsule should not intersect"
        );
    }
    #[test]
    fn obb_capsule_closest_point_on_obb() {
        let obb = Obb::from_aabb([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]);
        let q = [2.0, 0.0, 0.0];
        let cp = ObbCapsuleCollision::closest_point_on_obb(&obb, q);
        assert!((cp[0] - 1.0).abs() < 1e-9, "closest x={}", cp[0]);
    }
    #[test]
    fn obb_triangle_intersecting() {
        let obb = Obb::from_aabb([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]);
        let tri = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0]];
        let result = ObbTriangleCollision::test(&obb, tri);
        assert!(
            result.is_some(),
            "intersecting triangle should produce contact"
        );
    }
    #[test]
    fn obb_triangle_separated() {
        let obb = Obb::from_aabb([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]);
        let tri = [[5.0, 0.0, 0.0], [7.0, 0.0, 0.0], [6.0, 2.0, 0.0]];
        let result = ObbTriangleCollision::test(&obb, tri);
        assert!(result.is_none(), "separated triangle should return None");
    }
    #[test]
    fn edge_edge_contact_perpendicular() {
        let (cp, n) = edge_edge_contact(
            [-1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 2.0, 0.0],
        );
        assert!(length(cp) < 0.1, "cp={cp:?}");
        let n_len = length(n);
        assert!((n_len - 1.0).abs() < 1e-10, "normal not unit: {n_len}");
    }
    #[test]
    fn closest_points_on_segments_parallel() {
        let (t, s, pa, pb) = closest_points_on_segments(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 0.0, 0.0],
        );
        let _ = (t, s);
        let dist = length(sub(pa, pb));
        assert!((dist - 1.0).abs() < 1e-9, "dist={dist}");
    }
    #[test]
    fn sat_contact_points_obb_obb_face() {
        let a = Obb::from_aabb([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]);
        let b = Obb::from_aabb([0.5, -1.0, -1.0], [2.5, 1.0, 1.0]);
        let contact = ObbCollision::obb_obb_sat(&a, &b).expect("should collide");
        let points = SatContactPointGenerator::from_obb_obb(&a, &b, &contact);
        assert!(
            !points.is_empty(),
            "Should generate at least one contact point"
        );
    }
}
/// Query overlapping pairs between two flat-array OBB BVHs.
#[allow(dead_code)]
pub fn obb_bvh_pair_query(a: &ObbBvh, b: &ObbBvh) -> Vec<(usize, usize)> {
    if a.nodes.is_empty() || b.nodes.is_empty() {
        return Vec::new();
    }
    let mut results = Vec::new();
    obb_bvh_pair_recursive(a, b, 0, 0, &mut results);
    results
}
pub(super) fn obb_bvh_pair_recursive(
    ta: &ObbBvh,
    tb: &ObbBvh,
    ia: usize,
    ib: usize,
    results: &mut Vec<(usize, usize)>,
) {
    let na = &ta.nodes[ia];
    let nb = &tb.nodes[ib];
    if ObbCollision::obb_obb_sat(&na.bounds, &nb.bounds).is_none() {
        return;
    }
    if na.is_leaf() && nb.is_leaf() {
        if let (Some(ga), Some(gb)) = (na.geometry_index, nb.geometry_index) {
            results.push((ga, gb));
        }
        return;
    }
    let descend_a =
        !na.is_leaf() && (nb.is_leaf() || na.bounds.half_extents[0] >= nb.bounds.half_extents[0]);
    if descend_a {
        if na.left < ta.nodes.len() {
            obb_bvh_pair_recursive(ta, tb, na.left, ib, results);
        }
        if na.right < ta.nodes.len() {
            obb_bvh_pair_recursive(ta, tb, na.right, ib, results);
        }
    } else {
        if nb.left < tb.nodes.len() {
            obb_bvh_pair_recursive(ta, tb, ia, nb.left, results);
        }
        if nb.right < tb.nodes.len() {
            obb_bvh_pair_recursive(ta, tb, ia, nb.right, results);
        }
    }
}
/// Project an OBB onto a unit vector and return (center_projection, extent).
#[allow(dead_code)]
pub fn obb_project(obb: &Obb, axis: [f64; 3]) -> (f64, f64) {
    let cp = dot(obb.center, axis);
    let r = obb.half_extents[0] * dot(obb.axes[0], axis).abs()
        + obb.half_extents[1] * dot(obb.axes[1], axis).abs()
        + obb.half_extents[2] * dot(obb.axes[2], axis).abs();
    (cp, r)
}
/// Returns `true` if two OBBs overlap along `axis`.
#[allow(dead_code)]
pub fn obb_overlap_on_axis(a: &Obb, b: &Obb, axis: [f64; 3]) -> bool {
    let len = length(axis);
    if len < 1e-12 {
        return true;
    }
    let ax = scale(axis, 1.0 / len);
    let (ca, ra) = obb_project(a, ax);
    let (cb, rb) = obb_project(b, ax);
    (ca - cb).abs() <= ra + rb
}
#[cfg(test)]
mod extra_tests {

    use crate::ConvexPolyhedraSat;
    use crate::ConvexPolyhedron;

    use crate::Obb;
    use crate::ObbBvh;
    use crate::ObbBvhNode;

    use crate::SatAxisCache;

    use crate::obb_bvh_pair_query;
    use crate::obb_overlap_on_axis;
    use crate::obb_project;
    use crate::sat_collision::add;
    use crate::sat_collision::length;

    #[test]
    fn sat_cache_initially_invalid() {
        let cache = SatAxisCache::new();
        assert!(!cache.valid);
    }
    #[test]
    fn sat_cache_update_and_read() {
        let mut c = SatAxisCache::new();
        c.update([1.0, 0.0, 0.0], 0.5);
        assert!(c.valid);
        assert!((c.overlap - 0.5).abs() < 1e-12);
    }
    #[test]
    fn sat_cache_invalidate() {
        let mut c = SatAxisCache::new();
        c.update([0.0, 1.0, 0.0], 0.3);
        c.invalidate();
        assert!(!c.valid);
    }
    #[test]
    fn sat_cache_check_overlapping_obbs() {
        let a = Obb::from_aabb([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]);
        let b = Obb::from_aabb([0.5, -1.0, -1.0], [2.5, 1.0, 1.0]);
        let mut c = SatAxisCache::new();
        c.update([1.0, 0.0, 0.0], 0.5);
        assert!(c.check(&a, &b).is_some());
    }
    #[test]
    fn sat_cache_check_invalid_returns_none() {
        let a = Obb::from_aabb([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]);
        let b = Obb::from_aabb([0.5, -1.0, -1.0], [2.5, 1.0, 1.0]);
        let c = SatAxisCache::new();
        assert!(c.check(&a, &b).is_none());
    }
    #[test]
    fn sat_cache_default_invalid() {
        let c = SatAxisCache::default();
        assert!(!c.valid);
    }
    fn unit_cube_poly() -> ConvexPolyhedron {
        let vertices = vec![
            [-1.0, -1.0, -1.0],
            [1.0, -1.0, -1.0],
            [-1.0, 1.0, -1.0],
            [1.0, 1.0, -1.0],
            [-1.0, -1.0, 1.0],
            [1.0, -1.0, 1.0],
            [-1.0, 1.0, 1.0],
            [1.0, 1.0, 1.0],
        ];
        let face_normals = vec![
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0],
        ];
        let edge_directions = vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        ConvexPolyhedron::new(vertices, face_normals, edge_directions)
    }
    fn translate_poly(p: &ConvexPolyhedron, t: [f64; 3]) -> ConvexPolyhedron {
        let v: Vec<[f64; 3]> = p.vertices.iter().map(|v| add(*v, t)).collect();
        ConvexPolyhedron::new(v, p.face_normals.clone(), p.edge_directions.clone())
    }
    #[test]
    fn poly_support_max_x() {
        let c = unit_cube_poly();
        let s = c.support([1.0, 0.0, 0.0]);
        assert!((s[0] - 1.0).abs() < 1e-12, "s[0]={}", s[0]);
    }
    #[test]
    fn poly_project_x_axis() {
        let c = unit_cube_poly();
        let (mn, mx) = c.project([1.0, 0.0, 0.0]);
        assert!((mn + 1.0).abs() < 1e-12, "mn={mn}");
        assert!((mx - 1.0).abs() < 1e-12, "mx={mx}");
    }
    #[test]
    fn polyhedra_sat_overlapping() {
        let a = unit_cube_poly();
        let b = translate_poly(&a, [1.5, 0.0, 0.0]);
        assert!(ConvexPolyhedraSat::test(&a, &b).is_some());
    }
    #[test]
    fn polyhedra_sat_separated() {
        let a = unit_cube_poly();
        let b = translate_poly(&a, [5.0, 0.0, 0.0]);
        assert!(ConvexPolyhedraSat::is_separated(&a, &b));
    }
    #[test]
    fn polyhedra_sat_depth_non_negative() {
        let a = unit_cube_poly();
        let b = translate_poly(&a, [1.0, 0.0, 0.0]);
        let c = ConvexPolyhedraSat::test(&a, &b).unwrap();
        assert!(c.depth >= 0.0, "depth={}", c.depth);
    }
    #[test]
    fn polyhedra_sat_normal_unit_length() {
        let a = unit_cube_poly();
        let b = translate_poly(&a, [1.0, 0.5, 0.0]);
        let c = ConvexPolyhedraSat::test(&a, &b).unwrap();
        let len = length(c.normal);
        assert!((len - 1.0).abs() < 1e-9, "len={len}");
    }
    #[test]
    fn polyhedra_sat_empty_edge_list() {
        let a = ConvexPolyhedron::new(
            vec![[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[0.0, 0.0, 1.0]],
            vec![],
        );
        let b = ConvexPolyhedron::new(
            vec![[-1.0, 0.0, 0.5], [1.0, 0.0, 0.5], [0.0, 1.0, 0.5]],
            vec![[0.0, 0.0, 1.0]],
            vec![],
        );
        let _ = ConvexPolyhedraSat::test(&a, &b);
    }
    #[test]
    fn polyhedra_sat_touching_zero_depth() {
        let a = unit_cube_poly();
        let b = translate_poly(&a, [2.0, 0.0, 0.0]);
        if let Some(c) = ConvexPolyhedraSat::test(&a, &b) {
            assert!(c.depth >= 0.0, "touching depth must be >= 0: {}", c.depth);
        }
    }
    #[test]
    fn obb_bvh_empty() {
        let bvh = ObbBvh::build(&[]);
        assert_eq!(bvh.node_count(), 0);
    }
    #[test]
    fn obb_bvh_single_leaf() {
        let obb = Obb::from_aabb([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]);
        let bvh = ObbBvh::build(&[obb]);
        assert_eq!(bvh.leaf_count(), 1);
    }
    #[test]
    fn obb_bvh_multiple_leaves() {
        let obbs: Vec<Obb> = (0..4)
            .map(|i| {
                let x = i as f64 * 2.0;
                Obb::from_aabb([x, 0.0, 0.0], [x + 1.0, 1.0, 1.0])
            })
            .collect();
        let bvh = ObbBvh::build(&obbs);
        assert_eq!(bvh.leaf_count(), 4);
    }
    #[test]
    fn obb_bvh_query_hit() {
        let obbs = vec![
            Obb::from_aabb([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
            Obb::from_aabb([5.0, 0.0, 0.0], [6.0, 1.0, 1.0]),
        ];
        let bvh = ObbBvh::build(&obbs);
        let q = Obb::from_aabb([0.2, 0.2, 0.2], [0.8, 0.8, 0.8]);
        let hits = bvh.query_overlapping(&q);
        assert!(!hits.is_empty());
    }
    #[test]
    fn obb_bvh_query_miss() {
        let obbs = vec![
            Obb::from_aabb([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
            Obb::from_aabb([5.0, 0.0, 0.0], [6.0, 1.0, 1.0]),
        ];
        let bvh = ObbBvh::build(&obbs);
        let q = Obb::from_aabb([20.0, 0.0, 0.0], [21.0, 1.0, 1.0]);
        assert!(bvh.query_overlapping(&q).is_empty());
    }
    #[test]
    fn obb_bvh_pair_query_overlapping() {
        let a = ObbBvh::build(&[Obb::from_aabb([0.0, 0.0, 0.0], [1.0, 1.0, 1.0])]);
        let b = ObbBvh::build(&[Obb::from_aabb([0.5, 0.0, 0.0], [1.5, 1.0, 1.0])]);
        let pairs = obb_bvh_pair_query(&a, &b);
        assert!(!pairs.is_empty());
    }
    #[test]
    fn obb_bvh_pair_query_separated() {
        let a = ObbBvh::build(&[Obb::from_aabb([0.0, 0.0, 0.0], [1.0, 1.0, 1.0])]);
        let b = ObbBvh::build(&[Obb::from_aabb([10.0, 0.0, 0.0], [11.0, 1.0, 1.0])]);
        assert!(obb_bvh_pair_query(&a, &b).is_empty());
    }
    #[test]
    fn obb_bvh_large_build_correct_leaf_count() {
        let obbs: Vec<Obb> = (0..8)
            .map(|i| {
                let x = i as f64 * 2.0;
                Obb::from_aabb([x, 0.0, 0.0], [x + 1.0, 1.0, 1.0])
            })
            .collect();
        let bvh = ObbBvh::build(&obbs);
        assert_eq!(bvh.leaf_count(), 8);
    }
    #[test]
    fn obb_project_unit_cube_x() {
        let obb = Obb::from_aabb([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]);
        let (cp, r) = obb_project(&obb, [1.0, 0.0, 0.0]);
        assert!(cp.abs() < 1e-12, "cp={cp}");
        assert!((r - 1.0).abs() < 1e-12, "r={r}");
    }
    #[test]
    fn obb_overlap_on_axis_overlapping() {
        let a = Obb::from_aabb([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]);
        let b = Obb::from_aabb([0.5, -1.0, -1.0], [2.5, 1.0, 1.0]);
        assert!(obb_overlap_on_axis(&a, &b, [1.0, 0.0, 0.0]));
    }
    #[test]
    fn obb_overlap_on_axis_separated() {
        let a = Obb::from_aabb([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]);
        let b = Obb::from_aabb([3.0, -1.0, -1.0], [5.0, 1.0, 1.0]);
        assert!(!obb_overlap_on_axis(&a, &b, [1.0, 0.0, 0.0]));
    }
    #[test]
    fn obb_overlap_degenerate_zero_axis_returns_true() {
        let a = Obb::from_aabb([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]);
        let b = Obb::from_aabb([3.0, -1.0, -1.0], [5.0, 1.0, 1.0]);
        assert!(obb_overlap_on_axis(&a, &b, [0.0, 0.0, 0.0]));
    }
    #[test]
    fn obb_bvh_node_is_leaf() {
        let n = ObbBvhNode {
            bounds: Obb::from_aabb([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]),
            left: usize::MAX,
            right: usize::MAX,
            geometry_index: Some(0),
        };
        assert!(n.is_leaf());
    }
    #[test]
    fn obb_bvh_node_not_leaf() {
        let n = ObbBvhNode {
            bounds: Obb::from_aabb([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]),
            left: 1,
            right: 2,
            geometry_index: None,
        };
        assert!(!n.is_leaf());
    }
    #[test]
    fn polyhedra_sat_diagonal_overlap() {
        let a = unit_cube_poly();
        let b = translate_poly(&a, [1.0, 1.0, 0.0]);
        assert!(ConvexPolyhedraSat::test(&a, &b).is_some());
    }
    #[test]
    fn polyhedra_sat_z_axis_separation() {
        let a = unit_cube_poly();
        let b = translate_poly(&a, [0.0, 0.0, 5.0]);
        assert!(ConvexPolyhedraSat::is_separated(&a, &b));
    }
    #[test]
    fn polyhedra_sat_all_axes_separation() {
        let a = unit_cube_poly();
        let b = translate_poly(&a, [3.0, 3.0, 3.0]);
        assert!(ConvexPolyhedraSat::is_separated(&a, &b));
    }
}
