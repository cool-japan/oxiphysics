//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{
    CachedSupport, CsoPoint, GjkCache, GjkCacheRegistry, GjkDistanceResult, GjkTermination,
    JohnsonSubResult, NarrowphaseResult, ProximityResult, ShapeDesc,
};

pub(super) fn dot3_arr(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
pub(super) fn sub3_arr(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
pub(super) fn negate3(a: [f64; 3]) -> [f64; 3] {
    [-a[0], -a[1], -a[2]]
}
pub(super) fn scale3_arr(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
pub(super) fn len3_arr(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}
pub(super) fn cross3_arr(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
/// Simplified GJK do-simplex for cached intersection test.
/// Returns `Some(new_direction)` to continue, or `None` if origin is enclosed.
pub(super) fn do_simplex_cached(simplex: &mut Vec<CachedSupport>) -> Option<[f64; 3]> {
    match simplex.len() {
        2 => {
            let a = simplex[1].minkowski_point;
            let b = simplex[0].minkowski_point;
            let ab = sub3_arr(b, a);
            let ao = negate3(a);
            if dot3_arr(ab, ao) > 0.0 {
                let t = cross3_arr(cross3_arr(ab, ao), ab);
                Some(t)
            } else {
                simplex.remove(0);
                Some(ao)
            }
        }
        3 => {
            let a = simplex[2].minkowski_point;
            let b = simplex[1].minkowski_point;
            let c = simplex[0].minkowski_point;
            let ab = sub3_arr(b, a);
            let ac = sub3_arr(c, a);
            let ao = negate3(a);
            let abc = cross3_arr(ab, ac);
            if dot3_arr(cross3_arr(abc, ac), ao) > 0.0 {
                if dot3_arr(ac, ao) > 0.0 {
                    *simplex = vec![simplex[0], simplex[2]];
                    Some(cross3_arr(cross3_arr(ac, ao), ac))
                } else {
                    *simplex = vec![simplex[2]];
                    Some(ao)
                }
            } else if dot3_arr(cross3_arr(ab, abc), ao) > 0.0 {
                if dot3_arr(ab, ao) > 0.0 {
                    *simplex = vec![simplex[1], simplex[2]];
                    Some(cross3_arr(cross3_arr(ab, ao), ab))
                } else {
                    *simplex = vec![simplex[2]];
                    Some(ao)
                }
            } else if dot3_arr(abc, ao) > 0.0 {
                Some(abc)
            } else {
                simplex.swap(0, 1);
                Some(negate3(abc))
            }
        }
        4 => {
            let a = simplex[3].minkowski_point;
            let b = simplex[2].minkowski_point;
            let c = simplex[1].minkowski_point;
            let d = simplex[0].minkowski_point;
            let ab = sub3_arr(b, a);
            let ac = sub3_arr(c, a);
            let ad = sub3_arr(d, a);
            let ao = negate3(a);
            let abc = cross3_arr(ab, ac);
            let acd = cross3_arr(ac, ad);
            let adb = cross3_arr(ad, ab);
            if dot3_arr(abc, ao) > 0.0 {
                *simplex = vec![simplex[1], simplex[2], simplex[3]];
                Some(abc)
            } else if dot3_arr(acd, ao) > 0.0 {
                *simplex = vec![simplex[0], simplex[1], simplex[3]];
                Some(acd)
            } else if dot3_arr(adb, ao) > 0.0 {
                *simplex = vec![simplex[0], simplex[2], simplex[3]];
                Some(adb)
            } else {
                None
            }
        }
        _ => {
            if simplex.is_empty() {
                Some([1.0, 0.0, 0.0])
            } else {
                simplex.last().map(|s| negate3(s.minkowski_point))
            }
        }
    }
}
/// Compute the closest point on the simplex (line or triangle) to the origin.
pub(super) fn closest_point_on_simplex_to_origin(simplex: &[CachedSupport]) -> [f64; 3] {
    match simplex.len() {
        0 => [0.0; 3],
        1 => simplex[0].minkowski_point,
        2 => {
            let a = simplex[0].minkowski_point;
            let b = simplex[1].minkowski_point;
            let ab = sub3_arr(b, a);
            let ao = negate3(a);
            let t = dot3_arr(ao, ab) / dot3_arr(ab, ab);
            let t = t.clamp(0.0, 1.0);
            [a[0] + ab[0] * t, a[1] + ab[1] * t, a[2] + ab[2] * t]
        }
        _ => simplex[0].minkowski_point,
    }
}
#[cfg(test)]
mod tests {

    use crate::narrowphase::gjk_cache::CachedSupport;

    use crate::narrowphase::gjk_cache::GjkCache;

    use crate::narrowphase::gjk_cache::GjkCacheStats;

    use crate::narrowphase::gjk_cache::GjkStats;

    use crate::narrowphase::gjk_cache::GjkWarmStart;

    use crate::narrowphase::gjk_cache::len3_arr;
    use crate::narrowphase::gjk_cache::negate3;
    use oxiphysics_core::math::Vec3;

    use crate::narrowphase::gjk_cache::scale3_arr;
    use crate::narrowphase::gjk_cache::sub3_arr;

    #[test]
    fn test_cache_miss() {
        let cache = GjkWarmStart::new();
        assert!(
            cache.cached_query_hint(1, 2).is_none(),
            "cache miss must return None"
        );
    }
    #[test]
    fn test_cache_hit() {
        let mut cache = GjkWarmStart::new();
        let wa = Vec3::new(0.0, 0.0, 0.0);
        let wb = Vec3::new(3.0, 0.0, 0.0);
        cache.update_cache(1, 2, wa, wb);
        let hint = cache
            .cached_query_hint(1, 2)
            .expect("second call must return cached direction");
        assert!(hint.x > 0.9, "hint must point along +X, got {:?}", hint);
        assert!((hint.norm() - 1.0).abs() < 1e-6, "hint must be unit length");
    }
    #[test]
    fn test_cache_invalidation() {
        let mut cache = GjkWarmStart::new();
        cache.update_cache(10, 20, Vec3::new(1.0, 0.0, 0.0), Vec3::new(4.0, 0.0, 0.0));
        assert!(cache.cached_query_hint(10, 20).is_some());
        cache.invalidate(10, 20);
        assert!(
            cache.cached_query_hint(10, 20).is_none(),
            "after invalidation the hint must be None"
        );
    }
    #[test]
    fn test_stats() {
        let mut stats = GjkStats::new();
        assert_eq!(stats.warm_start_ratio(), 0.0);
        stats.record(false);
        stats.record(true);
        stats.record(true);
        assert_eq!(stats.total_queries, 3);
        assert_eq!(stats.warm_started, 2);
        let ratio = stats.warm_start_ratio();
        assert!(
            (ratio - 2.0 / 3.0).abs() < 1e-10,
            "ratio must be 2/3, got {}",
            ratio
        );
    }
    #[test]
    fn test_multiple_pairs() {
        let mut cache = GjkWarmStart::new();
        cache.update_cache(1, 2, Vec3::new(0.0, 0.0, 0.0), Vec3::new(5.0, 0.0, 0.0));
        cache.update_cache(3, 4, Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 5.0, 0.0));
        cache.update_cache(5, 6, Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 5.0));
        let h12 = cache.cached_query_hint(1, 2).unwrap();
        let h34 = cache.cached_query_hint(3, 4).unwrap();
        let h56 = cache.cached_query_hint(5, 6).unwrap();
        assert!(h12.x > 0.9, "pair (1,2) must point +X, got {:?}", h12);
        assert!(h34.y > 0.9, "pair (3,4) must point +Y, got {:?}", h34);
        assert!(h56.z > 0.9, "pair (5,6) must point +Z, got {:?}", h56);
        assert!(cache.cached_query_hint(99, 100).is_none());
        assert_eq!(cache.len(), 3);
    }
    fn sphere_support(center: [f64; 3], radius: f64, dir: [f64; 3]) -> [f64; 3] {
        let len = len3_arr(dir);
        if len < 1e-10 {
            return center;
        }
        let n = scale3_arr(dir, 1.0 / len);
        [
            center[0] + n[0] * radius,
            center[1] + n[1] * radius,
            center[2] + n[2] * radius,
        ]
    }
    #[test]
    fn test_gjk_cache_intersect_overlapping_spheres() {
        let mut cache = GjkCache::new();
        let c1 = [0.0, 0.0, 0.0];
        let c2 = [1.0, 0.0, 0.0];
        let r = 1.0;
        let result = cache.intersect_cached(|dir| {
            let sa = sphere_support(c1, r, dir);
            let sb = sphere_support(c2, r, negate3(dir));
            CachedSupport::new(sa, sb)
        });
        assert!(result, "overlapping spheres must intersect");
    }
    #[test]
    fn test_gjk_cache_intersect_separated_spheres() {
        let mut cache = GjkCache::new();
        let c1 = [0.0, 0.0, 0.0];
        let c2 = [5.0, 0.0, 0.0];
        let r = 1.0;
        let result = cache.intersect_cached(|dir| {
            let sa = sphere_support(c1, r, dir);
            let sb = sphere_support(c2, r, negate3(dir));
            CachedSupport::new(sa, sb)
        });
        assert!(!result, "separated spheres must not intersect");
    }
    #[test]
    fn test_gjk_cache_warm_start() {
        let mut cache = GjkCache::new();
        let c1 = [0.0, 0.0, 0.0];
        let c2 = [1.0, 0.0, 0.0];
        let r = 1.0;
        let _ = cache.intersect_cached(|dir| {
            let sa = sphere_support(c1, r, dir);
            let sb = sphere_support(c2, r, negate3(dir));
            CachedSupport::new(sa, sb)
        });
        assert!(!cache.simplex_vertices.is_empty());
        let prev = cache.simplex_vertices.clone();
        cache.warm_start_from(&prev);
        let result = cache.intersect_cached(|dir| {
            let sa = sphere_support(c1, r, dir);
            let sb = sphere_support(c2, r, negate3(dir));
            CachedSupport::new(sa, sb)
        });
        assert!(result);
    }
    #[test]
    fn test_gjk_cache_reset() {
        let mut cache = GjkCache::new();
        cache.simplex_vertices.push(CachedSupport::zero());
        cache.last_direction = [0.0, 1.0, 0.0];
        cache.reset();
        assert!(cache.simplex_vertices.is_empty());
        assert_eq!(cache.last_direction, [1.0, 0.0, 0.0]);
    }
    #[test]
    fn test_cached_support_new() {
        let s = CachedSupport::new([3.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert_eq!(s.minkowski_point, [2.0, 0.0, 0.0]);
    }
    #[test]
    fn test_gjk_cache_stats() {
        let mut stats = GjkCacheStats::new();
        stats.record_query(5, true, 10);
        stats.record_query(8, false, 10);
        assert_eq!(stats.queries, 2);
        assert_eq!(stats.iterations, 13);
        assert_eq!(stats.warm_start_saved, 5);
        assert!((stats.avg_iterations() - 6.5).abs() < 1e-10);
    }
    #[test]
    fn test_closest_points_separated() {
        let mut cache = GjkCache::new();
        let c1 = [0.0, 0.0, 0.0];
        let c2 = [5.0, 0.0, 0.0];
        let r = 1.0;
        let result = cache.closest_points_cached(|dir| {
            let sa = sphere_support(c1, r, dir);
            let sb = sphere_support(c2, r, negate3(dir));
            CachedSupport::new(sa, sb)
        });
        assert!(
            result.is_some(),
            "separated spheres should produce closest points"
        );
    }
    #[test]
    fn test_closest_points_intersecting() {
        let mut cache = GjkCache::new();
        let c1 = [0.0, 0.0, 0.0];
        let c2 = [0.5, 0.0, 0.0];
        let r = 1.0;
        let result = cache.closest_points_cached(|dir| {
            let sa = sphere_support(c1, r, dir);
            let sb = sphere_support(c2, r, negate3(dir));
            CachedSupport::new(sa, sb)
        });
        if let Some((pa, pb)) = result {
            let diff = sub3_arr(pa, pb);
            let dist = len3_arr(diff);
            assert!(
                dist.is_finite(),
                "witness distance should be finite, got {}",
                dist
            );
        }
    }
}
/// Compute the closest point on a line segment AB to the origin.
///
/// Returns barycentric weights (u, v) such that closest = u*A + v*B.
pub(super) fn johnson_segment(a: [f64; 3], b: [f64; 3]) -> ([f64; 3], f64, f64) {
    let ab = sub3_arr(b, a);
    let ab_sq = dot3_arr(ab, ab);
    if ab_sq < 1e-14 {
        return (a, 1.0, 0.0);
    }
    let t = -dot3_arr(a, ab) / ab_sq;
    let t = t.clamp(0.0, 1.0);
    let closest = [a[0] + ab[0] * t, a[1] + ab[1] * t, a[2] + ab[2] * t];
    (closest, 1.0 - t, t)
}
/// Compute the closest point on triangle ABC to the origin using Johnson's method.
///
/// Returns (closest, u, v, w) with barycentric coords.
pub(super) fn johnson_triangle(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> ([f64; 3], f64, f64, f64) {
    let ab = sub3_arr(b, a);
    let ac = sub3_arr(c, a);
    let n = cross3_arr(ab, ac);
    let n_sq = dot3_arr(n, n);
    if n_sq < 1e-14 {
        let (cp_ab, u_ab, v_ab) = johnson_segment(a, b);
        let (cp_ac, u_ac, _) = johnson_segment(a, c);
        let (cp_bc, _, v_bc) = johnson_segment(b, c);
        let d_ab = dot3_arr(cp_ab, cp_ab);
        let d_ac = dot3_arr(cp_ac, cp_ac);
        let d_bc = dot3_arr(cp_bc, cp_bc);
        if d_ab <= d_ac && d_ab <= d_bc {
            return (cp_ab, u_ab, v_ab, 0.0);
        } else if d_ac <= d_bc {
            return (cp_ac, u_ac, 0.0, 1.0 - u_ac);
        } else {
            return (cp_bc, 0.0, v_bc, 1.0 - v_bc);
        }
    }
    let ao = negate3(a);
    let dist_to_plane = dot3_arr(ao, n) / n_sq.sqrt();
    let d1 = dot3_arr(ab, ao);
    let d2 = dot3_arr(ac, ao);
    let d3 = dot3_arr(ab, ab);
    let d4 = dot3_arr(ab, ac);
    let d5 = dot3_arr(ac, ac);
    let denom = d3 * d5 - d4 * d4;
    if denom.abs() < 1e-14 {
        return johnson_triangle_fallback(a, b, c);
    }
    let v = (d1 * d5 - d2 * d4) / denom;
    let w = (d2 * d3 - d1 * d4) / denom;
    let u = 1.0 - v - w;
    if u >= 0.0 && v >= 0.0 && w >= 0.0 {
        let closest = [
            a[0] * u + b[0] * v + c[0] * w,
            a[1] * u + b[1] * v + c[1] * w,
            a[2] * u + b[2] * v + c[2] * w,
        ];
        let _ = dist_to_plane;
        return (closest, u, v, w);
    }
    let (cp_ab, u_ab, v_ab) = johnson_segment(a, b);
    let (cp_ac, u_ac, w_ac) = johnson_segment(a, c);
    let (cp_bc, v_bc, w_bc) = johnson_segment(b, c);
    let d_ab = dot3_arr(cp_ab, cp_ab);
    let d_ac = dot3_arr(cp_ac, cp_ac);
    let d_bc = dot3_arr(cp_bc, cp_bc);
    if d_ab <= d_ac && d_ab <= d_bc {
        (cp_ab, u_ab, v_ab, 0.0)
    } else if d_ac <= d_bc {
        (cp_ac, u_ac, 0.0, w_ac)
    } else {
        (cp_bc, 0.0, v_bc, w_bc)
    }
}
/// Fallback for degenerate triangles in johnson_triangle.
pub(super) fn johnson_triangle_fallback(
    a: [f64; 3],
    b: [f64; 3],
    c: [f64; 3],
) -> ([f64; 3], f64, f64, f64) {
    let (cp_ab, u_ab, v_ab) = johnson_segment(a, b);
    let (cp_ac, u_ac, w_ac) = johnson_segment(a, c);
    let d_ab = dot3_arr(cp_ab, cp_ab);
    let d_ac = dot3_arr(cp_ac, cp_ac);
    if d_ab <= d_ac {
        (cp_ab, u_ab, v_ab, 0.0)
    } else {
        (cp_ac, u_ac, 0.0, w_ac)
    }
}
/// Full Johnson sub-algorithm: find closest point on simplex to origin.
///
/// Handles 1-, 2-, 3-, and 4-vertex simplices.  Returns the closest point,
/// barycentric weights, and whether the origin is inside.
pub fn johnson_closest_point(simplex: &[CachedSupport]) -> JohnsonSubResult {
    match simplex.len() {
        0 => JohnsonSubResult {
            closest: [0.0; 3],
            bary: [0.0; 4],
            num_active: 0,
            origin_inside: false,
        },
        1 => {
            let p = simplex[0].minkowski_point;
            JohnsonSubResult {
                closest: p,
                bary: [1.0, 0.0, 0.0, 0.0],
                num_active: 1,
                origin_inside: len3_arr(p) < 1e-10,
            }
        }
        2 => {
            let (closest, u, v) =
                johnson_segment(simplex[0].minkowski_point, simplex[1].minkowski_point);
            JohnsonSubResult {
                closest,
                bary: [u, v, 0.0, 0.0],
                num_active: 2,
                origin_inside: len3_arr(closest) < 1e-10,
            }
        }
        3 => {
            let (closest, u, v, w) = johnson_triangle(
                simplex[0].minkowski_point,
                simplex[1].minkowski_point,
                simplex[2].minkowski_point,
            );
            JohnsonSubResult {
                closest,
                bary: [u, v, w, 0.0],
                num_active: 3,
                origin_inside: len3_arr(closest) < 1e-10,
            }
        }
        _ => {
            let a = simplex[0].minkowski_point;
            let b = simplex[1].minkowski_point;
            let c = simplex[2].minkowski_point;
            let d = simplex[3].minkowski_point;
            if origin_in_tetrahedron(a, b, c, d) {
                return JohnsonSubResult {
                    closest: [0.0; 3],
                    bary: [0.25; 4],
                    num_active: 4,
                    origin_inside: true,
                };
            }
            let faces = [(a, b, c), (a, b, d), (a, c, d), (b, c, d)];
            let bary_indices = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)];
            let mut best_dist = f64::INFINITY;
            let mut best_closest = a;
            let mut best_bary = [1.0_f64, 0.0, 0.0, 0.0];
            for (i, (fa, fb, fc)) in faces.iter().enumerate() {
                let (cp, u, v, w) = johnson_triangle(*fa, *fb, *fc);
                let dist = dot3_arr(cp, cp);
                if dist < best_dist {
                    best_dist = dist;
                    best_closest = cp;
                    let mut bary = [0.0_f64; 4];
                    bary[bary_indices[i].0] = u;
                    bary[bary_indices[i].1] = v;
                    bary[bary_indices[i].2] = w;
                    best_bary = bary;
                }
            }
            JohnsonSubResult {
                closest: best_closest,
                bary: best_bary,
                num_active: 4,
                origin_inside: false,
            }
        }
    }
}
/// Test whether the origin is inside a tetrahedron ABCD.
pub(super) fn origin_in_tetrahedron(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> bool {
    let ab = sub3_arr(b, a);
    let ac = sub3_arr(c, a);
    let ad = sub3_arr(d, a);
    let n_abc = cross3_arr(ab, ac);
    let n_abd = cross3_arr(ab, ad);
    let n_acd = cross3_arr(ac, ad);
    let n_bcd = cross3_arr(sub3_arr(c, b), sub3_arr(d, b));
    let ao = negate3(a);
    let bo = negate3(b);
    let s_abc = dot3_arr(n_abc, ao) * dot3_arr(n_abc, ad) > 0.0;
    let s_abd = dot3_arr(n_abd, ao) * dot3_arr(n_abd, ac) > 0.0;
    let s_acd = dot3_arr(n_acd, ao) * dot3_arr(n_acd, ab) > 0.0;
    let s_bcd = dot3_arr(n_bcd, bo) * dot3_arr(n_bcd, sub3_arr(a, b)) > 0.0;
    s_abc && s_abd && s_acd && s_bcd
}
/// Full GJK distance computation with Johnson's algorithm, warm starting, and
/// configurable termination criteria.
///
/// `support_fn(dir) -> CachedSupport` must return the support point on the
/// Minkowski difference `A ⊖ B` for direction `dir`.
pub fn gjk_distance<F>(
    support_fn: &mut F,
    cache: &mut GjkCache,
    termination: &GjkTermination,
) -> GjkDistanceResult
where
    F: FnMut([f64; 3]) -> CachedSupport,
{
    let warm_started = !cache.simplex_vertices.is_empty();
    let mut simplex: Vec<CachedSupport> = Vec::with_capacity(4);
    let init_dir = if warm_started {
        cache.last_direction
    } else {
        [1.0, 0.0, 0.0]
    };
    let first = support_fn(init_dir);
    simplex.push(first);
    let mut prev_dist = f64::INFINITY;
    for iter in 0..termination.max_iters {
        let result = johnson_closest_point(&simplex);
        if result.origin_inside {
            cache.simplex_vertices = simplex.clone();
            cache.stats.record_query(iter as u64 + 1, warm_started, 10);
            return GjkDistanceResult {
                distance: 0.0,
                point_a: first.point_a,
                point_b: first.point_b,
                iterations: iter + 1,
                warm_started,
                intersecting: true,
            };
        }
        let dist = result.distance();
        if termination.should_terminate(dist, prev_dist, iter) {
            let (pa, pb) = barycentric_witness(&simplex, &result.bary);
            cache.simplex_vertices = simplex;
            cache.last_direction = scale3_arr(negate3(result.closest), 1.0 / (dist + 1e-30));
            cache.stats.record_query(iter as u64 + 1, warm_started, 10);
            return GjkDistanceResult {
                distance: dist,
                point_a: pa,
                point_b: pb,
                iterations: iter + 1,
                warm_started,
                intersecting: false,
            };
        }
        prev_dist = dist;
        let direction = negate3(result.closest);
        let dir_len = len3_arr(direction);
        if dir_len < 1e-14 {
            let (pa, pb) = barycentric_witness(&simplex, &result.bary);
            cache.simplex_vertices = simplex;
            cache.stats.record_query(iter as u64 + 1, warm_started, 10);
            return GjkDistanceResult {
                distance: dist,
                point_a: pa,
                point_b: pb,
                iterations: iter + 1,
                warm_started,
                intersecting: dist < termination.proximity_tol,
            };
        }
        let norm_dir = scale3_arr(direction, 1.0 / dir_len);
        let new_pt = support_fn(norm_dir);
        let v_proj = dot3_arr(new_pt.minkowski_point, norm_dir);
        if v_proj < dist - termination.abs_tol {
            let (pa, pb) = barycentric_witness(&simplex, &result.bary);
            cache.simplex_vertices = simplex;
            cache.last_direction = norm_dir;
            cache.stats.record_query(iter as u64 + 1, warm_started, 10);
            return GjkDistanceResult {
                distance: dist,
                point_a: pa,
                point_b: pb,
                iterations: iter + 1,
                warm_started,
                intersecting: false,
            };
        }
        simplex.push(new_pt);
        prune_simplex_to_active(&mut simplex, &result.bary, result.num_active);
    }
    let result = johnson_closest_point(&simplex);
    let (pa, pb) = barycentric_witness(&simplex, &result.bary);
    cache.simplex_vertices = simplex;
    cache
        .stats
        .record_query(termination.max_iters as u64, warm_started, 10);
    GjkDistanceResult {
        distance: result.distance(),
        point_a: pa,
        point_b: pb,
        iterations: termination.max_iters,
        warm_started,
        intersecting: result.origin_inside,
    }
}
/// Compute witness points on shapes A and B from simplex and barycentric weights.
pub(super) fn barycentric_witness(
    simplex: &[CachedSupport],
    bary: &[f64; 4],
) -> ([f64; 3], [f64; 3]) {
    let n = simplex.len().min(4);
    let mut pa = [0.0_f64; 3];
    let mut pb = [0.0_f64; 3];
    let mut weight_sum = 0.0;
    for (s, &w) in simplex[..n].iter().zip(bary[..n].iter()) {
        pa[0] += w * s.point_a[0];
        pa[1] += w * s.point_a[1];
        pa[2] += w * s.point_a[2];
        pb[0] += w * s.point_b[0];
        pb[1] += w * s.point_b[1];
        pb[2] += w * s.point_b[2];
        weight_sum += w;
    }
    if weight_sum > 1e-14 {
        pa = scale3_arr(pa, 1.0 / weight_sum);
        pb = scale3_arr(pb, 1.0 / weight_sum);
    }
    (pa, pb)
}
/// Prune simplex vertices that have zero barycentric weight.
pub(super) fn prune_simplex_to_active(
    simplex: &mut Vec<CachedSupport>,
    _bary: &[f64; 4],
    num_active: usize,
) {
    if simplex.len() > num_active {
        simplex.truncate(num_active);
    }
    if simplex.len() > 4 {
        simplex.drain(0..simplex.len() - 4);
    }
}
/// Perform a full proximity query using GJK with warm starting.
///
/// Returns the distance and witness points between two convex shapes.
///
/// `support_a(dir)` returns the support of shape A in direction `dir`.
/// `support_b(dir)` returns the support of shape B in direction `dir`.
pub fn gjk_proximity<FA, FB>(
    support_a: &mut FA,
    support_b: &mut FB,
    cache: &mut GjkCache,
    termination: &GjkTermination,
) -> ProximityResult
where
    FA: FnMut([f64; 3]) -> [f64; 3],
    FB: FnMut([f64; 3]) -> [f64; 3],
{
    let mut combined = |dir: [f64; 3]| -> CachedSupport {
        let pa = support_a(dir);
        let neg_dir = negate3(dir);
        let pb = support_b(neg_dir);
        CachedSupport::new(pa, pb)
    };
    let result = gjk_distance(&mut combined, cache, termination);
    let normal = if result.distance > 1e-10 {
        let d = sub3_arr(result.point_a, result.point_b);
        let len = len3_arr(d);
        scale3_arr(d, 1.0 / len)
    } else {
        [0.0, 1.0, 0.0]
    };
    ProximityResult {
        distance: result.distance,
        point_a: result.point_a,
        point_b: result.point_b,
        normal,
        overlapping: result.intersecting,
    }
}
/// Compute a Minkowski-difference support point from two support functions.
///
/// `support_a(dir)` returns the support of shape A in direction `dir`.
/// `support_b(dir)` returns the support of shape B in direction `dir`
/// (the caller should negate `dir` before passing to shape B's own support).
pub fn support_minkowski_diff<FA, FB>(
    mut support_a: FA,
    mut support_b: FB,
    dir: [f64; 3],
) -> CsoPoint
where
    FA: FnMut([f64; 3]) -> [f64; 3],
    FB: FnMut([f64; 3]) -> [f64; 3],
{
    let a = support_a(dir);
    let b = support_b(dir);
    CsoPoint {
        p: sub3_arr(a, b),
        a_support: a,
        b_support: b,
    }
}
#[cfg(test)]
mod tests_extended {
    use super::super::types::SupportCache;
    use super::super::types::WarmStartedGjk;
    use super::*;
    use crate::gjk_enhanced::johnson_segment;
    use crate::gjk_enhanced::johnson_triangle;
    fn sphere_sup(center: [f64; 3], radius: f64, dir: [f64; 3]) -> [f64; 3] {
        let len = len3_arr(dir);
        if len < 1e-10 {
            return center;
        }
        let n = scale3_arr(dir, 1.0 / len);
        [
            center[0] + n[0] * radius,
            center[1] + n[1] * radius,
            center[2] + n[2] * radius,
        ]
    }
    fn aabb_sup(center: [f64; 3], he: [f64; 3], dir: [f64; 3]) -> [f64; 3] {
        let sign = |d: f64| if d >= 0.0 { 1.0_f64 } else { -1.0_f64 };
        [
            center[0] + he[0] * sign(dir[0]),
            center[1] + he[1] * sign(dir[1]),
            center[2] + he[2] * sign(dir[2]),
        ]
    }
    #[test]
    fn test_johnson_segment_midpoint() {
        let a = [-1.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let (u, v, cp) = johnson_segment(a, b);
        assert!(
            len3_arr(cp) < 1e-10,
            "closest should be origin, got {:?}",
            cp
        );
        assert!((u - 0.5).abs() < 1e-10);
        assert!((v - 0.5).abs() < 1e-10);
    }
    #[test]
    fn test_johnson_segment_clamped() {
        let a = [1.0, 0.0, 0.0];
        let b = [3.0, 0.0, 0.0];
        let (u, v, cp) = johnson_segment(a, b);
        assert!((cp[0] - 1.0).abs() < 1e-10, "expected x=1, got {:?}", cp);
        assert!((u - 1.0).abs() < 1e-10);
        assert!((v).abs() < 1e-10);
    }
    #[test]
    fn test_johnson_triangle_interior() {
        let a = [1.0, 0.0, 0.0];
        let b = [-0.5, 0.866, 0.0];
        let c = [-0.5, -0.866, 0.0];
        let (u, v, w, cp) = johnson_triangle(a, b, c);
        assert!(len3_arr(cp) < 0.15, "expected near origin, got {:?}", cp);
        assert!(
            u >= -1e-6 && v >= -1e-6 && w >= -1e-6,
            "bary coords must be non-negative"
        );
        assert!((u + v + w - 1.0).abs() < 1e-6, "bary coords must sum to 1");
    }
    #[test]
    fn test_johnson_closest_point_single_vertex() {
        let s = CachedSupport::new([2.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        let res = johnson_closest_point(&[s]);
        assert!((res.closest[0] - 2.0).abs() < 1e-10);
        assert!(!res.origin_inside);
    }
    #[test]
    fn test_johnson_closest_point_segment_origin_inside() {
        let sa = CachedSupport::new([-1.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        let sb = CachedSupport::new([1.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        let res = johnson_closest_point(&[sa, sb]);
        assert!(
            len3_arr(res.closest) < 1e-10,
            "expected origin, got {:?}",
            res.closest
        );
        assert!(res.origin_inside);
    }
    #[test]
    fn test_origin_in_tetrahedron_inside() {
        let a = [1.0, 0.0, -0.707];
        let b = [-1.0, 0.0, -0.707];
        let c = [0.0, 1.0, 0.707];
        let d = [0.0, -1.0, 0.707];
        assert!(origin_in_tetrahedron(a, b, c, d));
    }
    #[test]
    fn test_origin_in_tetrahedron_outside() {
        let a = [3.0, 0.0, 0.0];
        let b = [4.0, 1.0, 0.0];
        let c = [4.0, -1.0, 0.0];
        let d = [4.0, 0.0, 1.0];
        assert!(!origin_in_tetrahedron(a, b, c, d));
    }
    #[test]
    fn test_support_cache_hit() {
        let mut cache = SupportCache::<8>::new();
        let dir = [1.0, 0.0, 0.0];
        let support = [5.0, 0.0, 0.0];
        cache.insert(dir, support);
        let found = cache.lookup(dir, 0.9999);
        assert!(found.is_some(), "should find cached support");
        assert_eq!(found.unwrap(), support);
        assert!((cache.hit_rate() - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_support_cache_miss_different_direction() {
        let mut cache = SupportCache::<8>::new();
        cache.insert([1.0, 0.0, 0.0], [5.0, 0.0, 0.0]);
        let found = cache.lookup([0.0, 1.0, 0.0], 0.9999);
        assert!(found.is_none(), "perpendicular direction should be a miss");
    }
    #[test]
    fn test_support_cache_eviction() {
        let mut cache = SupportCache::<2>::new();
        cache.insert([1.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        cache.insert([0.0, 1.0, 0.0], [0.0, 1.0, 0.0]);
        cache.insert([0.0, 0.0, 1.0], [0.0, 0.0, 1.0]);
        assert_eq!(cache.len(), 2);
        assert!(cache.lookup([1.0, 0.0, 0.0], 0.9999).is_none());
        assert!(cache.lookup([0.0, 1.0, 0.0], 0.9999).is_some());
        assert!(cache.lookup([0.0, 0.0, 1.0], 0.9999).is_some());
    }
    #[test]
    fn test_support_cache_reset() {
        let mut cache = SupportCache::<8>::new();
        cache.insert([1.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        cache.reset();
        assert!(cache.is_empty());
        assert_eq!(cache.hit_rate(), 0.0);
    }
    #[test]
    fn test_termination_max_iters() {
        let t = GjkTermination {
            max_iters: 5,
            ..GjkTermination::default_criteria()
        };
        assert!(t.should_terminate(1.0, 1.0, 5));
        assert!(!t.should_terminate(1.0, 2.0, 4));
    }
    #[test]
    fn test_termination_proximity() {
        let t = GjkTermination::default_criteria();
        assert!(t.should_terminate(1e-11, 1.0, 0));
        assert!(!t.should_terminate(1.0, 2.0, 0));
    }
    #[test]
    fn test_termination_abs_tol() {
        let t = GjkTermination {
            abs_tol: 1e-4,
            ..GjkTermination::default_criteria()
        };
        assert!(t.should_terminate(1.00005, 1.0, 0));
        assert!(!t.should_terminate(1.001, 1.0, 0));
    }
    #[test]
    fn test_termination_rel_tol() {
        let t = GjkTermination {
            rel_tol: 1e-3,
            abs_tol: 1e-12,
            ..GjkTermination::default_criteria()
        };
        assert!(t.should_terminate(100.05, 100.0, 0));
    }
    #[test]
    fn test_gjk_distance_separated_spheres() {
        let c1 = [0.0, 0.0, 0.0];
        let c2 = [5.0, 0.0, 0.0];
        let r = 1.0;
        let mut cache = GjkCache::new();
        let term = GjkTermination::default_criteria();
        let result = gjk_distance(
            &mut |dir: [f64; 3]| {
                CachedSupport::new(sphere_sup(c1, r, dir), sphere_sup(c2, r, negate3(dir)))
            },
            &mut cache,
            &term,
        );
        assert!(!result.intersecting);
        assert!(
            (result.distance - 3.0).abs() < 0.1,
            "expected ~3.0, got {}",
            result.distance
        );
    }
    #[test]
    fn test_gjk_distance_touching_spheres() {
        let c1 = [0.0, 0.0, 0.0];
        let c2 = [2.0, 0.0, 0.0];
        let r = 1.0;
        let mut cache = GjkCache::new();
        let term = GjkTermination::default_criteria();
        let result = gjk_distance(
            &mut |dir: [f64; 3]| {
                CachedSupport::new(sphere_sup(c1, r, dir), sphere_sup(c2, r, negate3(dir)))
            },
            &mut cache,
            &term,
        );
        assert!(
            result.distance < 0.05,
            "touching spheres: expected ~0, got {}",
            result.distance
        );
    }
    #[test]
    fn test_gjk_proximity_separated_aabbs() {
        let mut cache = GjkCache::new();
        let term = GjkTermination::default_criteria();
        let c1 = [0.0, 0.0, 0.0];
        let he1 = [0.5, 0.5, 0.5];
        let c2 = [3.0, 0.0, 0.0];
        let he2 = [0.5, 0.5, 0.5];
        let result = gjk_proximity(
            &mut |dir| aabb_sup(c1, he1, dir),
            &mut |dir| aabb_sup(c2, he2, dir),
            &mut cache,
            &term,
        );
        assert!(!result.overlapping);
        assert!(
            (result.distance - 2.0).abs() < 0.1,
            "expected distance ≈ 2.0, got {}",
            result.distance
        );
    }
    #[test]
    fn test_gjk_proximity_normal_direction() {
        let mut cache = GjkCache::new();
        let term = GjkTermination::default_criteria();
        let c1 = [0.0, 0.0, 0.0];
        let c2 = [4.0, 0.0, 0.0];
        let r = 0.5;
        let result = gjk_proximity(
            &mut |dir| sphere_sup(c1, r, dir),
            &mut |dir| sphere_sup(c2, r, dir),
            &mut cache,
            &term,
        );
        assert!(!result.overlapping);
        assert!(
            result.normal[0].abs() > 0.7,
            "normal should be along X axis, got {:?}",
            result.normal
        );
    }
    #[test]
    fn test_registry_creates_separate_caches() {
        let mut registry = GjkCacheRegistry::new();
        let c1 = [0.0, 0.0, 0.0];
        let c2 = [5.0, 0.0, 0.0];
        let c3 = [0.0, 5.0, 0.0];
        let r = 0.5;
        let _ = registry.query(1, 2, &mut |dir| sphere_sup(c1, r, dir), &mut |dir| {
            sphere_sup(c2, r, dir)
        });
        let _ = registry.query(1, 3, &mut |dir| sphere_sup(c1, r, dir), &mut |dir| {
            sphere_sup(c3, r, dir)
        });
        assert_eq!(registry.len(), 2);
    }
    #[test]
    fn test_registry_remove_body() {
        let mut registry = GjkCacheRegistry::new();
        let c = [0.0, 0.0, 0.0];
        let r = 0.5;
        let far = [10.0, 0.0, 0.0];
        let _ = registry.query(1, 2, &mut |dir| sphere_sup(c, r, dir), &mut |dir| {
            sphere_sup(far, r, dir)
        });
        let _ = registry.query(1, 3, &mut |dir| sphere_sup(c, r, dir), &mut |dir| {
            sphere_sup(far, r, dir)
        });
        let _ = registry.query(2, 3, &mut |dir| sphere_sup(c, r, dir), &mut |dir| {
            sphere_sup(far, r, dir)
        });
        assert_eq!(registry.len(), 3);
        registry.remove_body(1);
        assert_eq!(registry.len(), 1, "only pair (2,3) should remain");
    }
    #[test]
    fn test_registry_clear() {
        let mut registry = GjkCacheRegistry::new();
        let c = [0.0, 0.0, 0.0];
        let r = 0.5;
        let far = [10.0, 0.0, 0.0];
        let _ = registry.query(1, 2, &mut |dir| sphere_sup(c, r, dir), &mut |dir| {
            sphere_sup(far, r, dir)
        });
        let _ = registry.query(3, 4, &mut |dir| sphere_sup(c, r, dir), &mut |dir| {
            sphere_sup(far, r, dir)
        });
        assert_eq!(registry.len(), 2);
        registry.clear();
        assert!(registry.is_empty());
    }
    #[test]
    fn test_registry_warm_start_reduces_iterations() {
        let mut registry = GjkCacheRegistry::new();
        let c1 = [0.0, 0.0, 0.0];
        let c2 = [4.0, 0.0, 0.0];
        let r = 0.5;
        let r1 = registry.query(1, 2, &mut |dir| sphere_sup(c1, r, dir), &mut |dir| {
            sphere_sup(c2, r, dir)
        });
        let r2 = registry.query(1, 2, &mut |dir| sphere_sup(c1, r, dir), &mut |dir| {
            sphere_sup(c2, r, dir)
        });
        assert!(
            (r1.distance - r2.distance).abs() < 0.5,
            "warm start should give similar result: r1={}, r2={}",
            r1.distance,
            r2.distance
        );
    }
    #[test]
    fn test_aggregate_stats() {
        let mut registry = GjkCacheRegistry::new();
        let c = [0.0, 0.0, 0.0];
        let r = 0.5;
        let far = [10.0, 0.0, 0.0];
        for i in 0..5u64 {
            let _ = registry.query(i, i + 10, &mut |dir| sphere_sup(c, r, dir), &mut |dir| {
                sphere_sup(far, r, dir)
            });
        }
        let stats = registry.aggregate_stats();
        assert!(
            stats.queries >= 5,
            "expected at least 5 queries, got {}",
            stats.queries
        );
        assert!(stats.iterations >= 5);
    }
    #[test]
    fn test_barycentric_witness_single_point() {
        let s = CachedSupport::new([2.0, 3.0, 4.0], [1.0, 1.0, 1.0]);
        let bary = [1.0, 0.0, 0.0, 0.0];
        let (pa, pb) = barycentric_witness(&[s], &bary);
        assert_eq!(pa, [2.0, 3.0, 4.0]);
        assert_eq!(pb, [1.0, 1.0, 1.0]);
    }
    #[test]
    fn test_barycentric_witness_two_points() {
        let s1 = CachedSupport::new([0.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        let s2 = CachedSupport::new([2.0, 0.0, 0.0], [0.0, 2.0, 0.0]);
        let bary = [0.5, 0.5, 0.0, 0.0];
        let (pa, pb) = barycentric_witness(&[s1, s2], &bary);
        assert!((pa[0] - 1.0).abs() < 1e-10);
        assert!((pb[1] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_warm_started_gjk_same_query_twice_uses_cache() {
        let mut wgjk = WarmStartedGjk::new();
        let c1 = [0.0_f64; 3];
        let c2 = [5.0, 0.0, 0.0_f64];
        let r = 1.0_f64;
        let sup_a = |dir: [f64; 3]| sphere_sup(c1, r, dir);
        let sup_b = |dir: [f64; 3]| sphere_sup(c2, r, negate3(dir));
        let hit1 = wgjk.test_overlap(sup_a, sup_b);
        assert!(!hit1, "separated spheres should not overlap");
        let hits_after_first = wgjk.cache_hits;
        let hit2 = wgjk.test_overlap(sup_a, sup_b);
        assert!(!hit2);
        assert!(
            wgjk.cache_hits > hits_after_first,
            "second call should register at least one cache hit"
        );
    }
    #[test]
    fn test_warm_started_gjk_hit_ratio_converges() {
        let mut wgjk = WarmStartedGjk::new();
        let c1 = [0.0_f64; 3];
        let c2 = [5.0, 0.0, 0.0_f64];
        let r = 1.0_f64;
        let sup_a = |dir: [f64; 3]| sphere_sup(c1, r, dir);
        let sup_b = |dir: [f64; 3]| sphere_sup(c2, r, negate3(dir));
        for _ in 0..10 {
            let _ = wgjk.test_overlap(sup_a, sup_b);
        }
        let ratio = wgjk.hit_ratio();
        assert!(
            ratio > 0.5,
            "after 10 repeated queries hit ratio should be high, got {}",
            ratio
        );
    }
    #[test]
    fn test_cso_point_minkowski_diff() {
        let pa = [3.0, 0.0, 0.0_f64];
        let pb = [1.0, 0.0, 0.0_f64];
        let cso = CsoPoint {
            p: sub3_arr(pa, pb),
            a_support: pa,
            b_support: pb,
        };
        assert_eq!(cso.p, [2.0, 0.0, 0.0]);
    }
    #[test]
    fn test_support_minkowski_diff_direction_reuse() {
        let c1 = [0.0_f64; 3];
        let c2 = [1.0, 0.0, 0.0_f64];
        let r = 1.0_f64;
        let dir = [1.0, 0.0, 0.0_f64];
        let cso1 = support_minkowski_diff(
            |d| sphere_sup(c1, r, d),
            |d| sphere_sup(c2, r, negate3(d)),
            dir,
        );
        let cso2 = support_minkowski_diff(
            |d| sphere_sup(c1, r, d),
            |d| sphere_sup(c2, r, negate3(d)),
            dir,
        );
        assert_eq!(cso1.p, cso2.p, "same direction must yield same CSO point");
    }
}
/// Validity threshold: if either body has moved more than this distance since
/// the cache was recorded, the cache is considered stale.
pub(super) const CACHE_VALID_DIST_SQ: f64 = 0.25;
/// Check whether a cached simplex is still valid given the new positions of the
/// bodies' representative points (e.g., centroids).
///
/// Returns `true` if the cache is still usable as a warm-start.
pub fn cache_still_valid(
    cached_pos_a: [f64; 3],
    cached_pos_b: [f64; 3],
    current_pos_a: [f64; 3],
    current_pos_b: [f64; 3],
) -> bool {
    let da = sub3_arr(current_pos_a, cached_pos_a);
    let db = sub3_arr(current_pos_b, cached_pos_b);
    let dist_sq_a = dot3_arr(da, da);
    let dist_sq_b = dot3_arr(db, db);
    dist_sq_a <= CACHE_VALID_DIST_SQ && dist_sq_b <= CACHE_VALID_DIST_SQ
}
/// Sphere–sphere specialised narrowphase (bypasses GJK for speed).
///
/// Returns `None` for the fast path when shapes are clearly separated,
/// `Some(result)` otherwise.
pub fn sphere_sphere_narrowphase(a: &ShapeDesc, b: &ShapeDesc) -> NarrowphaseResult {
    let d = sub3_arr(b.position, a.position);
    let dist = len3_arr(d);
    let radii_sum = a.radius + b.radius;
    let separation = dist - radii_sum;
    let normal = if dist > 1e-10 {
        scale3_arr(d, 1.0 / dist)
    } else {
        [0.0, 1.0, 0.0]
    };
    let point_a = [
        a.position[0] + normal[0] * a.radius,
        a.position[1] + normal[1] * a.radius,
        a.position[2] + normal[2] * a.radius,
    ];
    let point_b = [
        b.position[0] - normal[0] * b.radius,
        b.position[1] - normal[1] * b.radius,
        b.position[2] - normal[2] * b.radius,
    ];
    NarrowphaseResult {
        overlapping: separation < 0.0,
        distance: separation,
        point_a,
        point_b,
        normal,
        iterations: 0,
        warm_started: false,
    }
}
/// Dispatches a GJK proximity query for general convex shapes, with warm starting.
pub fn dispatch_gjk_query<FA, FB>(
    id_a: u64,
    id_b: u64,
    support_a: &mut FA,
    support_b: &mut FB,
    registry: &mut GjkCacheRegistry,
) -> NarrowphaseResult
where
    FA: FnMut([f64; 3]) -> [f64; 3],
    FB: FnMut([f64; 3]) -> [f64; 3],
{
    let result = registry.query(id_a, id_b, support_a, support_b);
    NarrowphaseResult {
        overlapping: result.overlapping,
        distance: result.distance,
        point_a: result.point_a,
        point_b: result.point_b,
        normal: result.normal,
        iterations: 0,
        warm_started: false,
    }
}
/// Compute the distance from the origin to the Minkowski-difference convex hull
/// defined by `simplex`, using Johnson's sub-algorithm.
///
/// Returns `(distance, closest_on_A, closest_on_B)`.
pub fn gjk_distance_from_simplex(simplex: &[CachedSupport]) -> (f64, [f64; 3], [f64; 3]) {
    let res = johnson_closest_point(simplex);
    let (pa, pb) = barycentric_witness(simplex, &res.bary);
    (res.distance(), pa, pb)
}
#[cfg(test)]
mod tests_extra {

    use crate::narrowphase::HitRateStats;
    use crate::narrowphase::NarrowphaseResult;
    use crate::narrowphase::PositionedGjkCache;
    use crate::narrowphase::ShapeDesc;
    use crate::narrowphase::SimplexCache;

    use crate::narrowphase::cache_still_valid;
    use crate::narrowphase::dispatch_gjk_query;
    use crate::narrowphase::gjk_cache::CachedSupport;

    use crate::narrowphase::gjk_cache::GjkCache;
    use crate::narrowphase::gjk_cache::GjkCacheRegistry;
    use crate::narrowphase::gjk_cache::GjkCacheStats;
    use crate::narrowphase::gjk_cache::GjkDistanceResult;
    use crate::narrowphase::gjk_cache::GjkPairCache;
    use crate::narrowphase::gjk_cache::GjkStats;
    use crate::narrowphase::gjk_cache::GjkTermination;
    use crate::narrowphase::gjk_cache::GjkWarmStart;
    use crate::narrowphase::gjk_cache::ShapeType;

    use crate::narrowphase::gjk_cache::cross3_arr;
    use crate::narrowphase::gjk_cache::dot3_arr;
    use crate::narrowphase::gjk_cache::len3_arr;
    use crate::narrowphase::gjk_cache::negate3;
    use oxiphysics_core::math::Vec3;

    use super::gjk_distance_from_simplex;
    use crate::narrowphase::gjk_cache::prune_simplex_to_active;
    use crate::narrowphase::gjk_cache::scale3_arr;
    use crate::narrowphase::gjk_cache::sub3_arr;

    use crate::narrowphase::sphere_sphere_narrowphase;

    fn sphere_sup_local(center: [f64; 3], radius: f64, dir: [f64; 3]) -> [f64; 3] {
        let len = len3_arr(dir);
        if len < 1e-10 {
            return center;
        }
        let n = scale3_arr(dir, 1.0 / len);
        [
            center[0] + n[0] * radius,
            center[1] + n[1] * radius,
            center[2] + n[2] * radius,
        ]
    }
    #[test]
    fn cache_valid_for_small_move() {
        let old_a = [0.0, 0.0, 0.0];
        let old_b = [5.0, 0.0, 0.0];
        let new_a = [0.1, 0.0, 0.0];
        let new_b = [5.1, 0.0, 0.0];
        assert!(cache_still_valid(old_a, old_b, new_a, new_b));
    }
    #[test]
    fn cache_invalid_for_large_move() {
        let old_a = [0.0, 0.0, 0.0];
        let old_b = [5.0, 0.0, 0.0];
        let new_a = [1.0, 0.0, 0.0];
        let new_b = [5.0, 0.0, 0.0];
        assert!(!cache_still_valid(old_a, old_b, new_a, new_b));
    }
    #[test]
    fn cache_invalid_when_b_moves_too_far() {
        let old_a = [0.0, 0.0, 0.0];
        let old_b = [5.0, 0.0, 0.0];
        let new_a = [0.0, 0.0, 0.0];
        let new_b = [7.0, 0.0, 0.0];
        assert!(!cache_still_valid(old_a, old_b, new_a, new_b));
    }
    #[test]
    fn positioned_cache_initially_invalid() {
        let c = PositionedGjkCache::new();
        assert!(!c.is_valid_for([0.0; 3], [1.0; 3]));
    }
    #[test]
    fn positioned_cache_valid_after_update_with_same_positions() {
        let mut c = PositionedGjkCache::new();
        c.update([0.0; 3], [5.0; 3]);
        assert!(c.is_valid_for([0.0; 3], [5.0; 3]));
    }
    #[test]
    fn positioned_cache_invalid_after_large_displacement() {
        let mut c = PositionedGjkCache::new();
        c.update([0.0; 3], [5.0; 3]);
        assert!(!c.is_valid_for([3.0; 3], [5.0; 3]));
    }
    #[test]
    fn positioned_cache_invalidate_clears() {
        let mut c = PositionedGjkCache::new();
        c.update([0.0; 3], [1.0; 3]);
        assert!(c.is_valid_for([0.0; 3], [1.0; 3]));
        c.invalidate();
        assert!(!c.is_valid_for([0.0; 3], [1.0; 3]));
    }
    #[test]
    fn hit_rate_stats_all_valid() {
        let mut s = HitRateStats::new();
        for _ in 0..10 {
            s.record_warm_attempt(true);
        }
        assert_eq!(s.effective_hit_rate(), 1.0);
        assert_eq!(s.staleness_rate(), 0.0);
    }
    #[test]
    fn hit_rate_stats_all_cold() {
        let mut s = HitRateStats::new();
        for _ in 0..5 {
            s.record_cold_start();
        }
        assert_eq!(s.effective_hit_rate(), 0.0);
    }
    #[test]
    fn hit_rate_stats_mixed() {
        let mut s = HitRateStats::new();
        s.record_cold_start();
        s.record_warm_attempt(true);
        s.record_warm_attempt(false);
        assert!((s.effective_hit_rate() - 1.0 / 3.0).abs() < 1e-10);
        assert!((s.staleness_rate() - 0.5).abs() < 1e-10);
    }
    #[test]
    fn hit_rate_stats_zero_queries_returns_zero() {
        let s = HitRateStats::new();
        assert_eq!(s.effective_hit_rate(), 0.0);
        assert_eq!(s.staleness_rate(), 0.0);
    }
    #[test]
    fn sphere_sphere_separated() {
        let a = ShapeDesc {
            id: 1,
            shape_type: ShapeType::Sphere,
            position: [0.0, 0.0, 0.0],
            radius: 1.0,
            half_extents: [0.0; 3],
        };
        let b = ShapeDesc {
            id: 2,
            shape_type: ShapeType::Sphere,
            position: [5.0, 0.0, 0.0],
            radius: 1.0,
            half_extents: [0.0; 3],
        };
        let r = sphere_sphere_narrowphase(&a, &b);
        assert!(!r.overlapping);
        assert!((r.distance - 3.0).abs() < 1e-10);
    }
    #[test]
    fn sphere_sphere_overlapping() {
        let a = ShapeDesc {
            id: 1,
            shape_type: ShapeType::Sphere,
            position: [0.0, 0.0, 0.0],
            radius: 2.0,
            half_extents: [0.0; 3],
        };
        let b = ShapeDesc {
            id: 2,
            shape_type: ShapeType::Sphere,
            position: [1.0, 0.0, 0.0],
            radius: 2.0,
            half_extents: [0.0; 3],
        };
        let r = sphere_sphere_narrowphase(&a, &b);
        assert!(r.overlapping, "overlapping spheres must report overlapping");
        assert!(
            r.distance < 0.0,
            "distance should be negative (penetration)"
        );
    }
    #[test]
    fn sphere_sphere_touching() {
        let a = ShapeDesc {
            id: 1,
            shape_type: ShapeType::Sphere,
            position: [0.0, 0.0, 0.0],
            radius: 1.0,
            half_extents: [0.0; 3],
        };
        let b = ShapeDesc {
            id: 2,
            shape_type: ShapeType::Sphere,
            position: [2.0, 0.0, 0.0],
            radius: 1.0,
            half_extents: [0.0; 3],
        };
        let r = sphere_sphere_narrowphase(&a, &b);
        assert!(
            r.distance.abs() < 1e-10,
            "touching spheres: distance should be ~0"
        );
    }
    #[test]
    fn sphere_sphere_normal_points_from_b_to_a() {
        let a = ShapeDesc {
            id: 1,
            shape_type: ShapeType::Sphere,
            position: [0.0, 0.0, 0.0],
            radius: 0.5,
            half_extents: [0.0; 3],
        };
        let b = ShapeDesc {
            id: 2,
            shape_type: ShapeType::Sphere,
            position: [3.0, 0.0, 0.0],
            radius: 0.5,
            half_extents: [0.0; 3],
        };
        let r = sphere_sphere_narrowphase(&a, &b);
        assert!(
            r.normal[0].abs() > 0.9,
            "normal should be along X, got {:?}",
            r.normal
        );
        let norm_len = len3_arr(r.normal);
        assert!((norm_len - 1.0).abs() < 1e-10, "normal must be unit length");
    }
    #[test]
    fn narrowphase_result_separated_constructor() {
        let pa = [1.0, 0.0, 0.0];
        let pb = [4.0, 0.0, 0.0];
        let r = NarrowphaseResult::separated(pa, pb, 3.0);
        assert!(!r.overlapping);
        assert!((r.distance - 3.0).abs() < 1e-10);
        assert!(
            (len3_arr(r.normal) - 1.0).abs() < 1e-10,
            "normal must be unit length"
        );
    }
    #[test]
    fn gjk_distance_from_single_vertex_simplex() {
        let s = CachedSupport::new([3.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        let (dist, pa, _pb) = gjk_distance_from_simplex(&[s]);
        assert!((dist - 3.0).abs() < 1e-10, "expected dist=3, got {}", dist);
        assert!((pa[0] - 3.0).abs() < 1e-10);
    }
    #[test]
    fn gjk_distance_from_segment_spanning_origin() {
        let sa = CachedSupport::new([-1.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        let sb = CachedSupport::new([1.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        let (dist, _pa, _pb) = gjk_distance_from_simplex(&[sa, sb]);
        assert!(
            dist < 1e-9,
            "origin on segment: expected dist≈0, got {}",
            dist
        );
    }
    #[test]
    fn dispatch_gjk_query_separated_spheres() {
        let c1 = [0.0, 0.0, 0.0_f64];
        let c2 = [5.0, 0.0, 0.0_f64];
        let r = 0.5_f64;
        let mut registry = GjkCacheRegistry::new();
        let result = dispatch_gjk_query(
            1,
            2,
            &mut |dir| sphere_sup_local(c1, r, dir),
            &mut |dir| sphere_sup_local(c2, r, dir),
            &mut registry,
        );
        assert!(!result.overlapping);
        assert!(
            result.distance > 0.0,
            "separated spheres: expected positive distance"
        );
        assert_eq!(registry.len(), 1);
    }
    #[test]
    fn dispatch_gjk_query_second_call_uses_cache() {
        let c1 = [0.0, 0.0, 0.0_f64];
        let c2 = [5.0, 0.0, 0.0_f64];
        let r = 0.5_f64;
        let mut registry = GjkCacheRegistry::new();
        let _r1 = dispatch_gjk_query(
            1,
            2,
            &mut |dir| sphere_sup_local(c1, r, dir),
            &mut |dir| sphere_sup_local(c2, r, dir),
            &mut registry,
        );
        let r2 = dispatch_gjk_query(
            1,
            2,
            &mut |dir| sphere_sup_local(c1, r, dir),
            &mut |dir| sphere_sup_local(c2, r, dir),
            &mut registry,
        );
        assert!(!r2.overlapping);
        assert_eq!(registry.len(), 1);
    }
    #[test]
    fn shape_desc_creation() {
        let desc = ShapeDesc {
            id: 42,
            shape_type: ShapeType::Aabb,
            position: [1.0, 2.0, 3.0],
            radius: 0.0,
            half_extents: [0.5, 0.5, 0.5],
        };
        assert_eq!(desc.id, 42);
        assert_eq!(desc.shape_type, ShapeType::Aabb);
    }
    #[test]
    fn shape_types_are_distinct() {
        assert_ne!(ShapeType::Sphere, ShapeType::Aabb);
        assert_ne!(ShapeType::ConvexHull, ShapeType::Capsule);
    }
    #[test]
    fn utility_functions_cross_product() {
        let x = [1.0_f64, 0.0, 0.0];
        let y = [0.0_f64, 1.0, 0.0];
        let z = cross3_arr(x, y);
        assert!((z[0]).abs() < 1e-10);
        assert!((z[1]).abs() < 1e-10);
        assert!((z[2] - 1.0).abs() < 1e-10, "x × y should be z");
    }
    #[test]
    fn utility_functions_dot_product() {
        let a = [3.0_f64, 4.0, 0.0];
        let b = [1.0_f64, 0.0, 0.0];
        assert!((dot3_arr(a, b) - 3.0).abs() < 1e-10);
    }
    #[test]
    fn utility_functions_length() {
        let a = [3.0_f64, 4.0, 0.0];
        assert!((len3_arr(a) - 5.0).abs() < 1e-10);
    }
    #[test]
    fn utility_functions_negate() {
        let a = [1.0_f64, -2.0, 3.0];
        let n = negate3(a);
        assert_eq!(n, [-1.0, 2.0, -3.0]);
    }
    #[test]
    fn utility_functions_scale() {
        let a = [1.0_f64, 2.0, 3.0];
        let s = scale3_arr(a, 2.0);
        assert_eq!(s, [2.0, 4.0, 6.0]);
    }
    #[test]
    fn utility_functions_sub() {
        let a = [3.0_f64, 2.0, 1.0];
        let b = [1.0_f64, 1.0, 1.0];
        assert_eq!(sub3_arr(a, b), [2.0, 1.0, 0.0]);
    }
    #[test]
    fn warm_start_key_is_symmetric() {
        let mut cache = GjkWarmStart::new();
        let wa = Vec3::new(0.0, 0.0, 0.0);
        let wb = Vec3::new(0.0, 0.0, 5.0);
        cache.update_cache(7, 3, wa, wb);
        let h1 = cache.cached_query_hint(3, 7);
        let h2 = cache.cached_query_hint(7, 3);
        assert!(
            h1.is_some() && h2.is_some(),
            "both orderings should hit the cache"
        );
        let h1 = h1.unwrap();
        let h2 = h2.unwrap();
        assert!((h1.norm() - 1.0).abs() < 1e-6);
        assert!((h2.norm() - 1.0).abs() < 1e-6);
    }
    #[test]
    fn registry_aggregate_stats_increases_with_queries() {
        let mut registry = GjkCacheRegistry::new();
        let c1 = [0.0, 0.0, 0.0_f64];
        let c2 = [4.0, 0.0, 0.0_f64];
        let r = 0.5_f64;
        let stats0 = registry.aggregate_stats();
        assert_eq!(stats0.queries, 0);
        for _ in 0..3 {
            let _ = registry.query(
                10,
                20,
                &mut |dir| sphere_sup_local(c1, r, dir),
                &mut |dir| sphere_sup_local(c2, r, dir),
            );
        }
        let stats3 = registry.aggregate_stats();
        assert!(
            stats3.queries >= 3,
            "expected >= 3 queries, got {}",
            stats3.queries
        );
    }
    #[test]
    fn registry_remove_pair_reduces_len() {
        let mut registry = GjkCacheRegistry::new();
        let c = [0.0, 0.0, 0.0_f64];
        let r = 0.5_f64;
        let far = [10.0, 0.0, 0.0_f64];
        let _ = registry.query(1, 2, &mut |d| sphere_sup_local(c, r, d), &mut |d| {
            sphere_sup_local(far, r, d)
        });
        let _ = registry.query(3, 4, &mut |d| sphere_sup_local(c, r, d), &mut |d| {
            sphere_sup_local(far, r, d)
        });
        assert_eq!(registry.len(), 2);
        registry.remove(1, 2);
        assert_eq!(registry.len(), 1);
    }
    #[test]
    fn gjk_cache_new_is_empty_and_warm_start_enabled() {
        let cache = GjkCache::new();
        assert!(cache.simplex_vertices.is_empty());
        assert_eq!(cache.last_direction, [1.0, 0.0, 0.0]);
        assert!(cache.warm_start);
    }
    #[test]
    fn gjk_cache_default_equals_new() {
        let c1 = GjkCache::new();
        let c2 = GjkCache::default();
        assert_eq!(c1.last_direction, c2.last_direction);
        assert_eq!(c1.warm_start, c2.warm_start);
    }
    #[test]
    fn gjk_cache_stats_zero_queries_avg_zero() {
        let stats = GjkCacheStats::new();
        assert_eq!(stats.avg_iterations(), 0.0);
    }
    #[test]
    fn gjk_cache_stats_multiple_records() {
        let mut stats = GjkCacheStats::new();
        stats.record_query(4, false, 10);
        stats.record_query(6, true, 10);
        stats.record_query(2, true, 10);
        assert_eq!(stats.queries, 3);
        assert_eq!(stats.iterations, 12);
        assert_eq!(stats.warm_start_saved, 4 + 8);
        assert!((stats.avg_iterations() - 4.0).abs() < 1e-10);
    }
    #[test]
    fn cached_support_zero_is_all_zeros() {
        let s = CachedSupport::zero();
        assert_eq!(s.point_a, [0.0; 3]);
        assert_eq!(s.point_b, [0.0; 3]);
        assert_eq!(s.minkowski_point, [0.0; 3]);
    }
    #[test]
    fn cached_support_minkowski_point_computed_correctly() {
        let s = CachedSupport::new([5.0, 3.0, 1.0], [2.0, 1.0, 0.0]);
        assert_eq!(s.minkowski_point, [3.0, 2.0, 1.0]);
    }
    #[test]
    fn gjk_pair_cache_new_is_not_valid() {
        let pc = GjkPairCache::new();
        assert!(!pc.is_valid());
    }
    #[test]
    fn gjk_pair_cache_default_matches_new() {
        let pc = GjkPairCache::default();
        assert!(!pc.is_valid());
        assert_eq!(pc.last_direction, [1.0, 0.0, 0.0]);
    }
    #[test]
    fn gjk_stats_record_warm_false_doesnt_add_warm_started() {
        let mut stats = GjkStats::new();
        stats.record(false);
        stats.record(false);
        assert_eq!(stats.total_queries, 2);
        assert_eq!(stats.warm_started, 0);
        assert_eq!(stats.warm_start_ratio(), 0.0);
    }
    #[test]
    fn gjk_stats_warm_start_ratio_all_warm() {
        let mut stats = GjkStats::new();
        for _ in 0..5 {
            stats.record(true);
        }
        assert_eq!(stats.warm_start_ratio(), 1.0);
    }
    #[test]
    fn simplex_cache_default_is_invalid() {
        let sc = SimplexCache::default();
        assert!(!sc.cache_valid);
        assert!(sc.direction().is_none());
    }
    #[test]
    fn simplex_cache_direction_when_witnesses_coincide() {
        let mut sc = SimplexCache::new();
        sc.witness_a = Vec3::new(1.0, 0.0, 0.0);
        sc.witness_b = Vec3::new(1.0, 0.0, 0.0);
        sc.cache_valid = true;
        assert!(sc.direction().is_none());
    }
    #[test]
    fn simplex_cache_direction_along_z() {
        let mut sc = SimplexCache::new();
        sc.witness_a = Vec3::new(0.0, 0.0, 0.0);
        sc.witness_b = Vec3::new(0.0, 0.0, 7.0);
        sc.cache_valid = true;
        let dir = sc.direction().expect("should have direction");
        assert!((dir.z - 1.0).abs() < 1e-6, "should point +Z, got {:?}", dir);
    }
    #[test]
    fn gjk_warm_start_len_and_is_empty() {
        let mut cache = GjkWarmStart::new();
        assert!(cache.is_empty());
        assert_eq!(cache.len(), 0);
        cache.update_cache(1, 2, Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0));
        assert!(!cache.is_empty());
        assert_eq!(cache.len(), 1);
    }
    #[test]
    fn gjk_warm_start_update_overwrites_previous() {
        let mut cache = GjkWarmStart::new();
        cache.update_cache(1, 2, Vec3::new(0.0, 0.0, 0.0), Vec3::new(3.0, 0.0, 0.0));
        let h1 = cache.cached_query_hint(1, 2).unwrap();
        cache.update_cache(1, 2, Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 4.0, 0.0));
        let h2 = cache.cached_query_hint(1, 2).unwrap();
        assert!(h1.x > 0.9, "h1 should be +X, got {:?}", h1);
        assert!(h2.y > 0.9, "h2 should be +Y, got {:?}", h2);
    }
    #[test]
    fn gjk_termination_loose_more_permissive_than_default() {
        let loose = GjkTermination::loose();
        let default = GjkTermination::default_criteria();
        assert!(loose.max_iters <= default.max_iters);
        assert!(loose.abs_tol > default.abs_tol);
    }
    #[test]
    fn gjk_termination_tight_stricter_than_default() {
        let tight = GjkTermination::tight();
        let default = GjkTermination::default_criteria();
        assert!(tight.max_iters >= default.max_iters);
        assert!(tight.abs_tol < default.abs_tol);
    }
    #[test]
    fn gjk_distance_result_separation_vector() {
        let pa = [3.0, 0.0, 0.0_f64];
        let pb = [1.0, 0.0, 0.0_f64];
        let result = GjkDistanceResult {
            distance: 2.0,
            point_a: pa,
            point_b: pb,
            iterations: 1,
            warm_started: false,
            intersecting: false,
        };
        let sv = result.separation_vector();
        assert_eq!(sv, [2.0, 0.0, 0.0]);
    }
    #[test]
    fn prune_simplex_does_not_grow_beyond_4() {
        let s = CachedSupport::zero();
        let mut simplex = vec![s; 6];
        let bary = [0.25, 0.25, 0.25, 0.25];
        prune_simplex_to_active(&mut simplex, &bary, 4);
        assert!(simplex.len() <= 4);
    }
    #[test]
    fn positioned_cache_update_twice_stays_valid_if_close() {
        let mut c = PositionedGjkCache::new();
        c.update([0.0; 3], [5.0; 3]);
        c.update([0.1; 3], [5.1; 3]);
        assert!(c.is_valid_for([0.1; 3], [5.1; 3]));
    }
    #[test]
    fn hit_rate_stats_100_valid_attempts() {
        let mut s = HitRateStats::new();
        for _ in 0..100 {
            s.record_warm_attempt(true);
        }
        assert_eq!(s.valid_hits, 100);
        assert_eq!(s.stale_misses, 0);
        assert_eq!(s.cold_starts, 0);
        assert_eq!(s.effective_hit_rate(), 1.0);
        assert_eq!(s.staleness_rate(), 0.0);
    }
}
/// Reduce a 4-point simplex to the smallest sub-simplex that contains
/// the closest point to the origin, discarding redundant vertices.
///
/// Returns the reduced simplex and the barycentric weights of the closest point.
pub fn gjk_reduce_simplex(simplex: &[CachedSupport]) -> (Vec<CachedSupport>, [f64; 4]) {
    let result = johnson_closest_point(simplex);
    let n = simplex.len().min(4);
    let mut active: Vec<CachedSupport> = Vec::new();
    let mut active_bary = [0.0_f64; 4];
    let mut k = 0usize;
    for (&s, &w) in simplex[..n].iter().zip(result.bary[..n].iter()) {
        if w > 1e-12 {
            active.push(s);
            active_bary[k] = w;
            k += 1;
        }
    }
    if active.is_empty() && !simplex.is_empty() {
        active.push(simplex[0]);
        active_bary[0] = 1.0;
    }
    (active, active_bary)
}
/// Compute a lower bound on the GJK distance.
///
/// For each support point in the simplex, the dot product of the point with
/// a direction `d` gives a lower bound on `sup(d)`.  This is used to detect
/// early termination.
pub fn gjk_lower_bound(simplex: &[CachedSupport], direction: [f64; 3]) -> f64 {
    simplex.iter().fold(f64::NEG_INFINITY, |acc, s| {
        acc.max(dot3_arr(s.minkowski_point, direction))
    })
}
/// Compute the distance from the origin to the closest vertex in the simplex.
pub fn simplex_vertex_min_dist(simplex: &[CachedSupport]) -> f64 {
    simplex
        .iter()
        .map(|s| len3_arr(s.minkowski_point))
        .fold(f64::INFINITY, f64::min)
}
/// Compute the distance from the origin to the farthest vertex in the simplex.
pub fn simplex_vertex_max_dist(simplex: &[CachedSupport]) -> f64 {
    simplex
        .iter()
        .map(|s| len3_arr(s.minkowski_point))
        .fold(f64::NEG_INFINITY, f64::max)
}
