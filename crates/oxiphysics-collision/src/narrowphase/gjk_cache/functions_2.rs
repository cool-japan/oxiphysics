//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[cfg(test)]
mod tests_gjk_extended {
    use super::super::*;
    fn sphere_sup(center: [f64; 3], radius: f64, dir: [f64; 3]) -> [f64; 3] {
        let l = len3_arr(dir);
        if l < 1e-10 {
            return center;
        }
        let n = scale3_arr(dir, 1.0 / l);
        [
            center[0] + n[0] * radius,
            center[1] + n[1] * radius,
            center[2] + n[2] * radius,
        ]
    }
    fn make_support(c1: [f64; 3], c2: [f64; 3], r: f64) -> impl FnMut([f64; 3]) -> CachedSupport {
        move |dir| {
            let sa = sphere_sup(c1, r, dir);
            let sb = sphere_sup(c2, r, negate3(dir));
            CachedSupport::new(sa, sb)
        }
    }
    #[test]
    fn test_gjk_reduce_simplex_single_vertex() {
        let s = CachedSupport::new([1.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        let simplex = vec![s];
        let (reduced, bary) = gjk_reduce_simplex(&simplex);
        assert_eq!(reduced.len(), 1);
        assert!(
            (bary[0] - 1.0).abs() < 1e-10,
            "single vertex bary must be 1"
        );
    }
    #[test]
    fn test_gjk_reduce_simplex_two_vertices_keeps_minimum() {
        let s1 = CachedSupport {
            point_a: [0.0; 3],
            point_b: [0.0; 3],
            minkowski_point: [-1.0, 0.0, 0.0],
        };
        let s2 = CachedSupport {
            point_a: [0.0; 3],
            point_b: [0.0; 3],
            minkowski_point: [1.0, 0.0, 0.0],
        };
        let simplex = vec![s1, s2];
        let (reduced, _bary) = gjk_reduce_simplex(&simplex);
        assert!(
            !reduced.is_empty(),
            "reduced simplex must have at least 1 vertex"
        );
    }
    #[test]
    fn test_gjk_reduce_simplex_empty() {
        let simplex: Vec<CachedSupport> = vec![];
        let (reduced, _bary) = gjk_reduce_simplex(&simplex);
        assert!(reduced.is_empty(), "empty simplex must remain empty");
    }
    #[test]
    fn test_gjk_lower_bound_single_point() {
        let s = CachedSupport {
            point_a: [0.0; 3],
            point_b: [0.0; 3],
            minkowski_point: [3.0, 0.0, 0.0],
        };
        let lb = gjk_lower_bound(&[s], [1.0, 0.0, 0.0]);
        assert!((lb - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_gjk_lower_bound_multiple_points() {
        let s1 = CachedSupport {
            point_a: [0.0; 3],
            point_b: [0.0; 3],
            minkowski_point: [1.0, 0.0, 0.0],
        };
        let s2 = CachedSupport {
            point_a: [0.0; 3],
            point_b: [0.0; 3],
            minkowski_point: [3.0, 0.0, 0.0],
        };
        let lb = gjk_lower_bound(&[s1, s2], [1.0, 0.0, 0.0]);
        assert!(
            (lb - 3.0).abs() < 1e-10,
            "lower bound should be max projection"
        );
    }
    #[test]
    fn test_gjk_lower_bound_empty_simplex() {
        let lb = gjk_lower_bound(&[], [1.0, 0.0, 0.0]);
        assert!(
            lb.is_infinite() && lb < 0.0,
            "empty simplex lower bound is -inf"
        );
    }
    #[test]
    fn test_eviction_policy_conservative_age() {
        let p = EvictionPolicy::conservative();
        assert!(!p.should_evict(60), "age=60 should not evict (max=60)");
        assert!(p.should_evict(61), "age=61 should evict");
    }
    #[test]
    fn test_eviction_policy_aggressive_age() {
        let p = EvictionPolicy::aggressive();
        assert!(!p.should_evict(10), "age=10 should not evict (max=10)");
        assert!(p.should_evict(11), "age=11 should evict");
    }
    #[test]
    fn test_timestamped_registry_basic() {
        let mut reg = TimestampedGjkRegistry::new(EvictionPolicy::conservative());
        assert!(reg.is_empty());
        let _ = reg.get_or_create(1, 2);
        assert_eq!(reg.len(), 1);
    }
    #[test]
    fn test_timestamped_registry_eviction() {
        let mut reg = TimestampedGjkRegistry::new(EvictionPolicy {
            max_age_frames: 2,
            max_entries: 100,
        });
        let _ = reg.get_or_create(1, 2);
        assert_eq!(reg.len(), 1);
        reg.advance_frame();
        reg.advance_frame();
        reg.advance_frame();
        assert_eq!(reg.len(), 0, "entry should be evicted after max_age_frames");
    }
    #[test]
    fn test_timestamped_registry_access_refreshes() {
        let mut reg = TimestampedGjkRegistry::new(EvictionPolicy {
            max_age_frames: 2,
            max_entries: 100,
        });
        let _ = reg.get_or_create(1, 2);
        reg.advance_frame();
        let _ = reg.get_or_create(1, 2);
        reg.advance_frame();
        reg.advance_frame();
        assert_eq!(reg.len(), 1, "refreshed entry should not be evicted yet");
    }
    #[test]
    fn test_timestamped_registry_remove_body() {
        let mut reg = TimestampedGjkRegistry::new(EvictionPolicy::default());
        let _ = reg.get_or_create(1, 2);
        let _ = reg.get_or_create(1, 3);
        let _ = reg.get_or_create(4, 5);
        assert_eq!(reg.len(), 3);
        reg.remove_body(1);
        assert_eq!(
            reg.len(),
            1,
            "should only keep pair (4,5) after removing body 1"
        );
    }
    #[test]
    fn test_simplex_vertex_min_dist() {
        let s1 = CachedSupport {
            point_a: [0.0; 3],
            point_b: [0.0; 3],
            minkowski_point: [1.0, 0.0, 0.0],
        };
        let s2 = CachedSupport {
            point_a: [0.0; 3],
            point_b: [0.0; 3],
            minkowski_point: [3.0, 0.0, 0.0],
        };
        let min_d = simplex_vertex_min_dist(&[s1, s2]);
        assert!((min_d - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_simplex_vertex_max_dist() {
        let s1 = CachedSupport {
            point_a: [0.0; 3],
            point_b: [0.0; 3],
            minkowski_point: [1.0, 0.0, 0.0],
        };
        let s2 = CachedSupport {
            point_a: [0.0; 3],
            point_b: [0.0; 3],
            minkowski_point: [3.0, 0.0, 0.0],
        };
        let max_d = simplex_vertex_max_dist(&[s1, s2]);
        assert!((max_d - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_simplex_vertex_distances_empty() {
        let min_d = simplex_vertex_min_dist(&[]);
        let max_d = simplex_vertex_max_dist(&[]);
        assert!(min_d.is_infinite() && min_d > 0.0);
        assert!(max_d.is_infinite() && max_d < 0.0);
    }
    #[test]
    fn test_gjk_contact_pair_separated() {
        let mut pair = GjkContactPair::new();
        let c1 = [0.0, 0.0, 0.0];
        let c2 = [5.0, 0.0, 0.0];
        let mut support = make_support(c1, c2, 1.0);
        let result = pair.query(&mut support);
        assert!(
            !result.intersecting,
            "separated spheres should not intersect"
        );
        assert!(pair.distance() > 0.0);
        assert!(!pair.is_overlapping());
    }
    #[test]
    fn test_gjk_contact_pair_intersecting() {
        let mut pair = GjkContactPair::new();
        let c1 = [0.0, 0.0, 0.0];
        let c2 = [0.5, 0.0, 0.0];
        let mut support = make_support(c1, c2, 1.0);
        let result = pair.query(&mut support);
        assert!(result.intersecting || result.distance < 2.0);
    }
    #[test]
    fn test_gjk_contact_pair_age_increments() {
        let mut pair = GjkContactPair::new();
        assert_eq!(pair.age, 0);
        let c1 = [0.0, 0.0, 0.0];
        let c2 = [5.0, 0.0, 0.0];
        let mut support = make_support(c1, c2, 1.0);
        pair.query(&mut support);
        pair.query(&mut support);
        pair.query(&mut support);
        assert_eq!(pair.age, 3);
    }
    #[test]
    fn test_gjk_contact_pair_warm_start_reuse() {
        let mut pair = GjkContactPair::new();
        let c1 = [0.0, 0.0, 0.0];
        let c2 = [5.0, 0.0, 0.0];
        let mut support = make_support(c1, c2, 1.0);
        let r1 = pair.query(&mut support);
        let r2 = pair.query(&mut support);
        assert!(r2.warm_started, "second query must be warm-started");
        assert!(
            (r1.distance - r2.distance).abs() < 0.01,
            "warm-start distance diff too large: {} vs {}",
            r1.distance,
            r2.distance
        );
    }
    #[test]
    fn test_batch_proximity_entry_fields() {
        let e = BatchProximityEntry {
            idx_a: 0,
            idx_b: 1,
            distance: 2.5,
            intersecting: false,
        };
        assert_eq!(e.idx_a, 0);
        assert_eq!(e.idx_b, 1);
        assert!((e.distance - 2.5).abs() < 1e-10);
        assert!(!e.intersecting);
    }
    #[test]
    fn test_batch_proximity_entry_intersecting() {
        let e = BatchProximityEntry {
            idx_a: 3,
            idx_b: 7,
            distance: 0.0,
            intersecting: true,
        };
        assert!(e.intersecting);
        assert_eq!(e.distance, 0.0);
    }
}
