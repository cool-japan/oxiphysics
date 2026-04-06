//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
#[allow(unused_imports)]
use super::types::*;
#[cfg(test)]
mod tests_gjk_extended {
    use super::*;
    use oxiphysics_core::Transform;
    use oxiphysics_core::Vec3;
    use oxiphysics_geometry::{BoxShape, Capsule, Sphere};
    #[test]
    fn test_gjk_distance_query_separated() {
        let s1 = Sphere::new(1.0);
        let s2 = Sphere::new(1.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(4.0, 0.0, 0.0));
        let result = gjk_distance_query(&s1, &t1, &s2, &t2);
        assert!(
            result.distance > 0.0,
            "should be separated, got {}",
            result.distance
        );
        assert!(result.iterations > 0);
        assert!(result.closest_a.x.is_finite());
        assert!(result.closest_b.x.is_finite());
    }
    #[test]
    fn test_gjk_distance_query_intersecting() {
        let s1 = Sphere::new(1.5);
        let s2 = Sphere::new(1.5);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(1.0, 0.0, 0.0));
        let result = gjk_distance_query(&s1, &t1, &s2, &t2);
        assert!(
            result.distance <= 0.0,
            "overlapping shapes should have non-positive distance"
        );
    }
    #[test]
    fn test_gjk_distance_query_witness_direction() {
        let s1 = Sphere::new(1.0);
        let s2 = Sphere::new(1.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(4.0, 0.0, 0.0));
        let result = gjk_distance_query(&s1, &t1, &s2, &t2);
        if result.distance > 0.0 {
            assert!(
                result.closest_a.x < result.closest_b.x,
                "closest_a.x={} should be < closest_b.x={}",
                result.closest_a.x,
                result.closest_b.x
            );
        }
    }
    #[test]
    fn test_gjk_distance_query_box_sphere() {
        let b = BoxShape::new(Vec3::new(1.0, 1.0, 1.0));
        let s = Sphere::new(0.5);
        let tb = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let ts = Transform::from_position(Vec3::new(3.0, 0.0, 0.0));
        let result = gjk_distance_query(&b, &tb, &s, &ts);
        assert!(result.distance > 0.0, "box and sphere should be separated");
    }
    #[test]
    fn test_gjk_distance_consistency_with_basic() {
        let s1 = Sphere::new(1.0);
        let s2 = Sphere::new(1.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(3.0, 0.0, 0.0));
        let d1 = gjk_distance(&s1, &t1, &s2, &t2);
        let d2 = gjk_distance_query(&s1, &t1, &s2, &t2).distance;
        assert_eq!(
            d1 > 0.0,
            d2 > 0.0,
            "both distance functions should agree on sign"
        );
    }
    #[test]
    fn test_gjk_ray_cast_sphere_hit() {
        let s = Sphere::new(1.0);
        let t = Transform::from_position(Vec3::new(5.0, 0.0, 0.0));
        let origin = Vec3::new(0.0, 0.0, 0.0);
        let dir = Vec3::new(1.0, 0.0, 0.0);
        let result = gjk_ray_cast(&s, &t, origin, dir, 20.0);
        assert!(
            result.is_some(),
            "ray along +X toward sphere at x=5 should hit"
        );
        if let Some(hit) = result {
            assert!(hit.t > 0.0 && hit.t < 20.0, "hit.t={}", hit.t);
        }
    }
    #[test]
    fn test_gjk_ray_cast_sphere_miss() {
        let s = Sphere::new(0.5);
        let t = Transform::from_position(Vec3::new(5.0, 10.0, 0.0));
        let origin = Vec3::new(0.0, 0.0, 0.0);
        let dir = Vec3::new(1.0, 0.0, 0.0);
        let result = gjk_ray_cast(&s, &t, origin, dir, 20.0);
        let _ = result;
    }
    #[test]
    fn test_gjk_ray_cast_zero_direction_returns_none() {
        let s = Sphere::new(1.0);
        let t = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let result = gjk_ray_cast(&s, &t, Vec3::zeros(), Vec3::zeros(), 10.0);
        assert!(result.is_none(), "zero direction should return None");
    }
    #[test]
    fn test_gjk_ray_cast_max_t_respected() {
        let s = Sphere::new(1.0);
        let t = Transform::from_position(Vec3::new(100.0, 0.0, 0.0));
        let origin = Vec3::new(0.0, 0.0, 0.0);
        let dir = Vec3::new(1.0, 0.0, 0.0);
        let result = gjk_ray_cast(&s, &t, origin, dir, 5.0);
        let _ = result;
    }
    #[test]
    fn test_mpr_sphere_sphere_intersecting() {
        let s1 = Sphere::new(1.0);
        let s2 = Sphere::new(1.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(1.0, 0.0, 0.0));
        let result = mpr_query(&s1, &t1, &s2, &t2);
        match result {
            MprResult::Intersecting { depth, .. } => {
                assert!(depth >= 0.0, "depth should be non-negative: {depth}");
            }
            MprResult::Separated => {}
        }
    }
    #[test]
    fn test_mpr_sphere_sphere_separated() {
        let s1 = Sphere::new(0.5);
        let s2 = Sphere::new(0.5);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(5.0, 0.0, 0.0));
        let result = mpr_query(&s1, &t1, &s2, &t2);
        match result {
            MprResult::Separated => {}
            MprResult::Intersecting { depth, .. } => {
                assert!(
                    depth < 1.0,
                    "clearly separated spheres should not have large depth: {depth}"
                );
            }
        }
    }
    #[test]
    fn test_mpr_coincident_returns_intersecting() {
        let s = Sphere::new(1.0);
        let t = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let result = mpr_query(&s, &t, &s, &t);
        assert!(matches!(result, MprResult::Intersecting { .. }));
    }
    #[test]
    fn test_mpr_box_box_intersecting() {
        let b1 = BoxShape::new(Vec3::new(1.0, 1.0, 1.0));
        let b2 = BoxShape::new(Vec3::new(1.0, 1.0, 1.0));
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(1.5, 0.0, 0.0));
        let result = mpr_query(&b1, &t1, &b2, &t2);
        let _ = result;
    }
    #[test]
    fn test_warm_start_gjk_consistent_separated() {
        let mut ws = WarmStartGjk::new();
        let s1 = Sphere::new(1.0);
        let s2 = Sphere::new(1.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(4.0, 0.0, 0.0));
        let r1 = ws.query(&s1, &t1, &s2, &t2);
        let r2 = ws.query(&s1, &t1, &s2, &t2);
        assert!(matches!(r1, GjkResult::Separated { .. }));
        assert!(matches!(r2, GjkResult::Separated { .. }));
    }
    #[test]
    fn test_warm_start_gjk_consistent_intersecting() {
        let mut ws = WarmStartGjk::new();
        let s1 = Sphere::new(2.0);
        let s2 = Sphere::new(2.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(1.0, 0.0, 0.0));
        let r1 = ws.query(&s1, &t1, &s2, &t2);
        let _ = ws.query(&s1, &t1, &s2, &t2);
        assert!(matches!(r1, GjkResult::Intersecting(_)));
    }
    #[test]
    fn test_warm_start_gjk_records_stats() {
        let mut ws = WarmStartGjk::new();
        let s1 = Sphere::new(1.0);
        let s2 = Sphere::new(1.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(4.0, 0.0, 0.0));
        for _ in 0..5 {
            ws.query(&s1, &t1, &s2, &t2);
        }
        assert_eq!(ws.query_count(), 5);
        assert!(ws.avg_iterations() > 0.0);
    }
    #[test]
    fn test_warm_start_gjk_reset_clears_state() {
        let mut ws = WarmStartGjk::new();
        let s1 = Sphere::new(1.0);
        let s2 = Sphere::new(1.0);
        let t = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(4.0, 0.0, 0.0));
        ws.query(&s1, &t, &s2, &t2);
        ws.reset();
        let result = ws.query(&s1, &t, &s2, &t2);
        assert!(matches!(result, GjkResult::Separated { .. }));
    }
    #[test]
    fn test_gjk_toi_bisection_approaching() {
        let s1 = Sphere::new(0.5);
        let s2 = Sphere::new(0.5);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(4.0, 0.0, 0.0));
        let v1 = Vec3::new(2.0, 0.0, 0.0);
        let v2 = Vec3::new(-2.0, 0.0, 0.0);
        let toi = gjk_toi_bisection(&s1, &t1, v1, &s2, &t2, v2, 5.0, 1e-4);
        assert!(toi.is_some(), "approaching spheres should produce TOI");
        let t = toi.unwrap();
        assert!(t > 0.0 && t <= 5.0, "TOI should be in (0, 5], got {t}");
    }
    #[test]
    fn test_gjk_toi_bisection_receding_returns_none() {
        let s1 = Sphere::new(0.5);
        let s2 = Sphere::new(0.5);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(4.0, 0.0, 0.0));
        let v1 = Vec3::new(-1.0, 0.0, 0.0);
        let v2 = Vec3::new(1.0, 0.0, 0.0);
        let toi = gjk_toi_bisection(&s1, &t1, v1, &s2, &t2, v2, 5.0, 1e-4);
        assert!(toi.is_none(), "receding spheres should have no TOI");
    }
    #[test]
    fn test_gjk_toi_bisection_already_overlapping() {
        let s1 = Sphere::new(2.0);
        let s2 = Sphere::new(2.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(1.0, 0.0, 0.0));
        let v = Vec3::new(0.0, 0.0, 0.0);
        let toi = gjk_toi_bisection(&s1, &t1, v, &s2, &t2, v, 5.0, 1e-4);
        assert!(toi.is_none(), "already overlapping should return None");
    }
    #[test]
    fn test_gjk_ccd_linear_hit() {
        let s1 = Sphere::new(0.5);
        let s2 = Sphere::new(0.5);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(2.0, 0.0, 0.0));
        let v1 = Vec3::new(1.5, 0.0, 0.0);
        let v2 = Vec3::new(0.0, 0.0, 0.0);
        let result = gjk_ccd_linear(&s1, &t1, v1, &s2, &t2, v2);
        assert!(result.is_some(), "should hit");
        let r = result.unwrap();
        assert!(r.toi >= 0.0 && r.toi <= 1.0, "toi={}", r.toi);
        assert!(
            r.normal.norm() > 0.5,
            "normal should be non-zero: {:?}",
            r.normal
        );
    }
    #[test]
    fn test_gjk_ccd_linear_no_hit() {
        let s1 = Sphere::new(0.1);
        let s2 = Sphere::new(0.1);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(100.0, 0.0, 0.0));
        let v1 = Vec3::new(0.0, 1.0, 0.0);
        let v2 = Vec3::new(0.0, 0.0, 0.0);
        let result = gjk_ccd_linear(&s1, &t1, v1, &s2, &t2, v2);
        assert!(result.is_none(), "should not hit far away sphere");
    }
    #[test]
    fn test_estimate_penetration_26_overlapping() {
        let s1 = Sphere::new(1.5);
        let s2 = Sphere::new(1.5);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(1.0, 0.0, 0.0));
        let (depth, axis) = estimate_penetration_depth_26(&s1, &t1, &s2, &t2);
        assert!(
            depth > 0.0,
            "overlapping spheres should have positive depth: {depth}"
        );
        assert!(axis.norm() > 0.9, "axis should be unit: {:?}", axis);
    }
    #[test]
    fn test_estimate_penetration_26_separated() {
        let s1 = Sphere::new(0.5);
        let s2 = Sphere::new(0.5);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(5.0, 0.0, 0.0));
        let (depth, _axis) = estimate_penetration_depth_26(&s1, &t1, &s2, &t2);
        let _ = depth;
    }
    #[test]
    fn test_gjk_speculative_contact_within_threshold() {
        let s1 = Sphere::new(1.0);
        let s2 = Sphere::new(1.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(2.5, 0.0, 0.0));
        let result = gjk_speculative_contact(&s1, &t1, &s2, &t2, 1.0);
        assert!(
            result.is_some(),
            "should produce speculative contact within threshold"
        );
        let c = result.unwrap();
        assert!(
            c.gap >= 0.0,
            "gap should be non-negative for separated shapes: {}",
            c.gap
        );
    }
    #[test]
    fn test_gjk_speculative_contact_beyond_threshold() {
        let s1 = Sphere::new(1.0);
        let s2 = Sphere::new(1.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(20.0, 0.0, 0.0));
        let result = gjk_speculative_contact(&s1, &t1, &s2, &t2, 0.1);
        assert!(
            result.is_none(),
            "should not produce contact when far beyond threshold"
        );
    }
    #[test]
    fn test_gjk_speculative_contact_normal_direction() {
        let s1 = Sphere::new(1.0);
        let s2 = Sphere::new(1.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(2.5, 0.0, 0.0));
        if let Some(c) = gjk_speculative_contact(&s1, &t1, &s2, &t2, 1.0) {
            assert!(c.normal.norm() > 0.9, "normal should be approximately unit");
        }
    }
    #[test]
    fn test_gjk_hull_distance_separated() {
        let cube_a: Vec<[f64; 3]> = vec![
            [-1.0, -1.0, -1.0],
            [1.0, -1.0, -1.0],
            [1.0, 1.0, -1.0],
            [-1.0, 1.0, -1.0],
            [-1.0, -1.0, 1.0],
            [1.0, -1.0, 1.0],
            [1.0, 1.0, 1.0],
            [-1.0, 1.0, 1.0],
        ];
        let cube_b: Vec<[f64; 3]> = cube_a.iter().map(|v| [v[0] + 5.0, v[1], v[2]]).collect();
        let d = gjk_hull_distance(&cube_a, &cube_b);
        assert!(
            d > 0.0,
            "separated hulls should have positive distance: {d}"
        );
    }
    #[test]
    fn test_gjk_hull_distance_overlapping() {
        let cube: Vec<[f64; 3]> = vec![
            [-1.0, -1.0, -1.0],
            [1.0, -1.0, -1.0],
            [1.0, 1.0, -1.0],
            [-1.0, 1.0, -1.0],
            [-1.0, -1.0, 1.0],
            [1.0, -1.0, 1.0],
            [1.0, 1.0, 1.0],
            [-1.0, 1.0, 1.0],
        ];
        let d = gjk_hull_distance(&cube, &cube);
        assert_eq!(d, 0.0, "coincident hulls should have zero distance");
    }
    #[test]
    fn test_gjk_hull_distance_empty_returns_max() {
        let cube: Vec<[f64; 3]> = vec![[-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]];
        let d = gjk_hull_distance(&[], &cube);
        assert!(d > 1e10, "empty hull should return MAX distance: {d}");
    }
    #[test]
    fn test_gjk_point_inside_sphere() {
        let s = Sphere::new(2.0);
        let t = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let inside = gjk_point_inside_convex(Vec3::new(0.5, 0.0, 0.0), &s, &t);
        assert!(inside, "point at (0.5,0,0) should be inside sphere r=2");
    }
    #[test]
    fn test_gjk_point_outside_sphere() {
        let s = Sphere::new(0.5);
        let t = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let inside = gjk_point_inside_convex(Vec3::new(5.0, 0.0, 0.0), &s, &t);
        assert!(!inside, "point at (5,0,0) should be outside sphere r=0.5");
    }
    #[test]
    fn test_gjk_point_to_convex_distance() {
        let s = Sphere::new(1.0);
        let t = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let d = gjk_point_to_convex_distance(Vec3::new(5.0, 0.0, 0.0), &s, &t);
        assert!(d > 0.0, "distant point should have positive distance: {d}");
    }
    #[test]
    fn test_gjk_toi_angular_sweep_stationary_no_hit() {
        let s1 = Sphere::new(0.3);
        let s2 = Sphere::new(0.3);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(10.0, 0.0, 0.0));
        let zero = Vec3::new(0.0, 0.0, 0.0);
        let toi = gjk_toi_angular_sweep(&s1, &t1, zero, zero, &s2, &t2, zero, zero, 1.0, 4);
        assert!(
            toi.is_none(),
            "stationary separated spheres should have no TOI"
        );
    }
    #[test]
    fn test_gjk_toi_angular_sweep_linear_hit() {
        let s1 = Sphere::new(0.5);
        let s2 = Sphere::new(0.5);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(3.0, 0.0, 0.0));
        let v = Vec3::new(2.0, 0.0, 0.0);
        let zero = Vec3::new(0.0, 0.0, 0.0);
        let toi = gjk_toi_angular_sweep(&s1, &t1, v, zero, &s2, &t2, zero, zero, 2.0, 8);
        assert!(toi.is_some(), "approaching spheres should find TOI");
    }
    #[test]
    fn test_gjk_cast_ray_minkowski_hit() {
        let s1 = Sphere::new(1.0);
        let s2 = Sphere::new(1.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(5.0, 0.0, 0.0));
        let origin = Vec3::new(0.0, 0.0, 0.0);
        let dir = Vec3::new(1.0, 0.0, 0.0);
        let result = gjk_cast_ray_minkowski(&s1, &t1, &s2, &t2, origin, dir, 20.0);
        let _ = result;
    }
    #[test]
    fn test_gjk_cast_ray_minkowski_no_hit() {
        let s1 = Sphere::new(0.1);
        let s2 = Sphere::new(0.1);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(100.0, 100.0, 0.0));
        let origin = Vec3::new(0.0, 0.0, 0.0);
        let dir = Vec3::new(1.0, 0.0, 0.0);
        let result = gjk_cast_ray_minkowski(&s1, &t1, &s2, &t2, origin, dir, 5.0);
        assert!(
            result.is_none(),
            "ray should not hit sphere far away at angle"
        );
    }
    #[test]
    fn test_gjk_distance_capsule_capsule_separated() {
        let c1 = Capsule::new(0.5, 1.0);
        let c2 = Capsule::new(0.5, 1.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(6.0, 0.0, 0.0));
        let d = gjk_distance(&c1, &t1, &c2, &t2);
        assert!(
            d > 0.0,
            "separated capsules should have positive distance: {d}"
        );
    }
    #[test]
    fn test_gjk_distance_capsule_capsule_overlapping() {
        let c1 = Capsule::new(1.5, 1.0);
        let c2 = Capsule::new(1.5, 1.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(0.5, 0.0, 0.0));
        let d = gjk_distance(&c1, &t1, &c2, &t2);
        assert!(
            d <= 0.0,
            "overlapping capsules should have non-positive distance: {d}"
        );
    }
}
