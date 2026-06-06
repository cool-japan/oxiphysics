//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[cfg(test)]
mod tests_closest_points {

    use super::super::*;

    use crate::gjk_enhanced::EnhBox;

    use crate::gjk_enhanced::EnhSphere;

    #[test]
    fn test_closest_point_on_sphere() {
        let s = EnhSphere::new([0.0; 3], 1.0);
        let pt = closest_point_on_sphere(&s, [3.0, 0.0, 0.0]);
        assert!((pt[0] - 1.0).abs() < 1e-9);
    }
    #[test]
    fn test_closest_point_on_box_inside() {
        let b = EnhBox::new([0.0; 3], [2.0; 3]);
        let pt = closest_point_on_box(&b, [0.5, 0.5, 0.5]);
        assert_eq!(pt, [0.5, 0.5, 0.5]);
    }
    #[test]
    fn test_closest_point_on_box_outside() {
        let b = EnhBox::new([0.0; 3], [1.0; 3]);
        let pt = closest_point_on_box(&b, [5.0, 0.0, 0.0]);
        assert!((pt[0] - 1.0).abs() < 1e-9);
    }
    #[test]
    fn test_closest_point_on_segment_midpoint() {
        let (pt, t) = closest_point_on_segment([-1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        assert!((t - 0.5).abs() < 1e-9);
        assert!((pt[0]).abs() < 1e-9);
    }
    #[test]
    fn test_segment_segment_distance_parallel() {
        let (_, _, dist) = segment_segment_distance(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
        );
        assert!((dist - 1.0).abs() < 1e-9, "dist={dist}");
    }
    #[test]
    fn test_segment_segment_distance_perpendicular() {
        let (_, _, dist) = segment_segment_distance(
            [-1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, -1.0, 1.0],
            [0.0, 1.0, 1.0],
        );
        assert!((dist - 1.0).abs() < 0.01, "dist={dist}");
    }
    #[test]
    fn test_gjk_closest_points_separated() {
        let s1 = EnhSphere::new([0.0; 3], 1.0);
        let s2 = EnhSphere::new([5.0, 0.0, 0.0], 1.0);
        let pts = gjk_closest_points(&s1, &s2);
        assert!(pts.is_some());
        let (a, b) = pts.unwrap();
        assert!((a[0] - 1.0).abs() < 0.05);
        assert!((b[0] - 4.0).abs() < 0.05);
    }
}
#[cfg(test)]
mod tests_support_cache {

    use crate::gjk_enhanced::EnhSphere;
    use crate::gjk_enhanced::SupportCache;

    use crate::gjk_enhanced::gjk_with_support_cache;

    #[test]
    fn test_support_cache_miss_initially() {
        let mut cache = SupportCache::new();
        let r = cache.get([1.0, 0.0, 0.0], 1e-6);
        assert!(r.is_none());
    }
    #[test]
    fn test_support_cache_hit() {
        let mut cache = SupportCache::new();
        let sa = [1.0, 0.0, 0.0];
        let sb = [-1.0, 0.0, 0.0];
        cache.put([1.0, 0.0, 0.0], sa, sb);
        let r = cache.get([1.0, 0.0, 0.0], 1e-6);
        assert!(r.is_some());
        let (ca, cb) = r.unwrap();
        assert_eq!(ca, sa);
        assert_eq!(cb, sb);
    }
    #[test]
    fn test_support_cache_miss_different_dir() {
        let mut cache = SupportCache::new();
        cache.put([1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]);
        let r = cache.get([0.0, 1.0, 0.0], 1e-6);
        assert!(r.is_none());
    }
    #[test]
    fn test_support_cache_hit_rate() {
        let mut cache = SupportCache::new();
        cache.put([1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]);
        let _ = cache.get([1.0, 0.0, 0.0], 1e-6);
        let _ = cache.get([0.0, 1.0, 0.0], 1e-6);
        assert!((cache.hit_rate() - 0.5).abs() < 1e-9);
    }
    #[test]
    fn test_gjk_with_support_cache_separated() {
        let s1 = EnhSphere::new([0.0; 3], 1.0);
        let s2 = EnhSphere::new([5.0, 0.0, 0.0], 1.0);
        let mut cache = SupportCache::new();
        let res = gjk_with_support_cache(&s1, &s2, &mut cache, 64);
        assert!(!res.intersecting);
    }
}
#[cfg(test)]
mod tests_epa_pipeline {

    use crate::gjk_enhanced::EnhSphere;

    use crate::gjk_enhanced::gjk_epa_pipeline;

    use crate::gjk_enhanced::vlen;
    #[test]
    fn test_gjk_epa_pipeline_separated() {
        let s1 = EnhSphere::new([0.0; 3], 1.0);
        let s2 = EnhSphere::new([5.0, 0.0, 0.0], 1.0);
        let res = gjk_epa_pipeline(&s1, &s2, None, 64, 32, 1e-4);
        assert!(!res.in_contact);
        assert!(res.signed_distance > 0.0);
    }
    #[test]
    fn test_gjk_epa_pipeline_intersecting() {
        let s1 = EnhSphere::new([0.0; 3], 2.0);
        let s2 = EnhSphere::new([1.0, 0.0, 0.0], 2.0);
        let res = gjk_epa_pipeline(&s1, &s2, None, 64, 32, 1e-4);
        assert!(res.in_contact);
        assert!(
            res.signed_distance < 0.0,
            "signed_distance={}",
            res.signed_distance
        );
    }
    #[test]
    fn test_gjk_epa_normal_unit_length() {
        let s1 = EnhSphere::new([0.0; 3], 2.0);
        let s2 = EnhSphere::new([1.0, 0.0, 0.0], 2.0);
        let res = gjk_epa_pipeline(&s1, &s2, None, 64, 32, 1e-4);
        let nl = vlen(res.normal);
        assert!((nl - 1.0).abs() < 0.01 || nl < 0.01, "normal len={nl}");
    }
}
#[cfg(test)]
mod tests_toi {
    use super::super::types::*;
    use super::super::*;

    #[test]
    fn test_toi_approaching_spheres() {
        let s1 = EnhSphere::new([0.0, 0.0, 0.0], 0.5);
        let s2 = EnhSphere::new([1.9, 0.0, 0.0], 0.5);
        let res = gjk_time_of_impact(&s1, &s2, [1.0, 0.0, 0.0], 64, 1e-3);
        assert!(res.hit, "approaching spheres should hit");
        assert!(res.toi <= 1.0);
    }
    #[test]
    fn test_toi_diverging_spheres() {
        let s1 = EnhSphere::new([0.0, 0.0, 0.0], 0.5);
        let s2 = EnhSphere::new([5.0, 0.0, 0.0], 0.5);
        let res = gjk_time_of_impact(&s1, &s2, [-1.0, 0.0, 0.0], 32, 1e-3);
        assert!(!res.hit, "diverging spheres should not hit");
    }
    #[test]
    fn test_toi_lateral_motion() {
        let s1 = EnhSphere::new([0.0, 0.0, 0.0], 0.5);
        let s2 = EnhSphere::new([0.0, 5.0, 0.0], 0.5);
        let res = gjk_time_of_impact(&s1, &s2, [1.0, 0.0, 0.0], 32, 1e-3);
        assert!(!res.hit, "lateral motion should not hit");
    }
}
#[cfg(test)]
mod tests_minkowski_sum {

    use crate::compound_shapes::GjkPolytope;
    use crate::compound_shapes::MinkowskiSum;

    fn box_polytope(center: [f64; 3], half: [f64; 3]) -> GjkPolytope {
        let mut verts = Vec::with_capacity(8);
        for &sx in &[-1.0_f64, 1.0] {
            for &sy in &[-1.0_f64, 1.0] {
                for &sz in &[-1.0_f64, 1.0] {
                    verts.push([
                        center[0] + sx * half[0],
                        center[1] + sy * half[1],
                        center[2] + sz * half[2],
                    ]);
                }
            }
        }
        GjkPolytope::new(verts)
    }
    fn sphere_polytope(center: [f64; 3], radius: f64) -> GjkPolytope {
        let mut verts = Vec::new();
        let n = 12;
        for i in 0..n {
            let phi = std::f64::consts::PI * (i as f64 / (n - 1) as f64);
            for j in 0..(2 * n) {
                let theta = 2.0 * std::f64::consts::PI * (j as f64 / (2 * n) as f64);
                verts.push([
                    center[0] + radius * phi.sin() * theta.cos(),
                    center[1] + radius * phi.sin() * theta.sin(),
                    center[2] + radius * phi.cos(),
                ]);
            }
        }
        GjkPolytope::new(verts)
    }
    #[test]
    fn test_minkowski_sum_sphere_box() {
        let bp = box_polytope([0.0; 3], [1.0; 3]);
        let sp = sphere_polytope([0.0; 3], 0.5);
        let ms = MinkowskiSum::new(bp, sp, [0.0; 3]);
        // support_difference computes s_A(d) - s_B(-d)
        // For box at origin with half=[1,1,1] and sphere at origin with r=0.5:
        // s_A([1,0,0]) = [1,1,1] (vertex), s_B([-1,0,0]) ≈ [-0.5,0,0]
        // result ≈ [1,1,1] - [-0.5,0,0] = [1.5,1,1]
        let sup = ms.support_difference([1.0, 0.0, 0.0]);
        assert!(sup[0] > 1.4, "x component got {}", sup[0]);
    }
    #[test]
    fn test_minkowski_sum_origin_inside() {
        let bp = box_polytope([0.0; 3], [1.0; 3]);
        let sp = sphere_polytope([0.0; 3], 0.5);
        let ms = MinkowskiSum::new(bp, sp, [0.0; 3]);
        assert!(
            ms.origin_in_md(),
            "overlapping box and sphere should contain origin"
        );
    }
}
#[cfg(test)]
mod tests_warm_start {

    use crate::gjk_enhanced::SupportPoint;
    use crate::gjk_enhanced::WarmStartCache;

    #[test]
    fn test_warm_start_cache_invalidate() {
        let mut cache = WarmStartCache::new();
        cache.valid = true;
        cache
            .points
            .push(SupportPoint::new([1.0, 0.0, 0.0], [0.0; 3], [0.0; 3]));
        cache.invalidate();
        assert!(!cache.valid);
        assert!(cache.points.is_empty());
    }
    #[test]
    fn test_warm_start_cache_update() {
        let mut cache = WarmStartCache::new();
        let pts = vec![SupportPoint::new(
            [1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
        )];
        cache.update(&pts, [1.0, 0.0, 0.0]);
        assert!(cache.valid);
        assert_eq!(cache.points.len(), 1);
    }
}
#[cfg(test)]
mod tests_simplex {
    use super::super::types::*;

    #[test]
    fn test_simplex_empty_dim() {
        let s = EnhSimplex::new();
        assert_eq!(s.dim(), -1);
    }
    #[test]
    fn test_simplex_add_and_dim() {
        let mut s = EnhSimplex::new();
        s.add(SupportPoint::new([1.0, 0.0, 0.0], [0.0; 3], [0.0; 3]));
        assert_eq!(s.dim(), 0);
        s.add(SupportPoint::new([-1.0, 0.0, 0.0], [0.0; 3], [0.0; 3]));
        assert_eq!(s.dim(), 1);
    }
    #[test]
    fn test_simplex_point_contains_origin() {
        let mut s = EnhSimplex::new();
        s.add(SupportPoint::new([0.0, 0.0, 0.0], [0.0; 3], [0.0; 3]));
        let (_dir, contains) = s.reduce();
        assert!(contains);
    }
    #[test]
    fn test_simplex_line_no_origin() {
        let mut s = EnhSimplex::new();
        s.add(SupportPoint::new([2.0, 0.0, 0.0], [0.0; 3], [0.0; 3]));
        s.add(SupportPoint::new([3.0, 0.0, 0.0], [0.0; 3], [0.0; 3]));
        let (dir, contains) = s.reduce();
        assert!(!contains);
        assert!(dir[0] < 0.0, "dir should point toward origin, got {dir:?}");
    }
}
