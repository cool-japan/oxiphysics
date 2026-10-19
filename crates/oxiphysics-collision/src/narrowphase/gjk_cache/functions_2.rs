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
#[cfg(test)]
mod tests_gjk_distance_termination {
    use super::super::*;

    fn sphere_sup(center: [f64; 3], radius: f64, dir: [f64; 3]) -> [f64; 3] {
        let l = len3_arr(dir);
        if l < 1e-12 {
            return [center[0] + radius, center[1], center[2]];
        }
        let n = scale3_arr(dir, radius / l);
        [center[0] + n[0], center[1] + n[1], center[2] + n[2]]
    }

    fn box_sup(center: [f64; 3], half: [f64; 3], dir: [f64; 3]) -> [f64; 3] {
        let sign = |d: f64| if d >= 0.0 { 1.0 } else { -1.0 };
        [
            center[0] + half[0] * sign(dir[0]),
            center[1] + half[1] * sign(dir[1]),
            center[2] + half[2] * sign(dir[2]),
        ]
    }

    /// Accuracy asserted for queries run with the default termination criteria
    /// (`rel_tol = 1e-6`): the stagnation test may stop a few multiples of the
    /// relative tolerance short of the exact gap.
    const DEFAULT_ACCURACY: f64 = 1e-5;

    fn sphere_pair_distance_with(
        c1: [f64; 3],
        c2: [f64; 3],
        r1: f64,
        r2: f64,
        termination: &GjkTermination,
    ) -> GjkDistanceResult {
        let mut cache = GjkCache::new();
        gjk_distance(
            &mut |dir: [f64; 3]| {
                CachedSupport::new(sphere_sup(c1, r1, dir), sphere_sup(c2, r2, negate3(dir)))
            },
            &mut cache,
            termination,
        )
    }

    fn sphere_pair_distance(c1: [f64; 3], c2: [f64; 3], r1: f64, r2: f64) -> GjkDistanceResult {
        sphere_pair_distance_with(c1, c2, r1, r2, &GjkTermination::default_criteria())
    }

    #[test]
    fn spheres_offset_along_y_report_the_true_gap() {
        // The initial search direction (1, 0, 0) is not along the line of centres,
        // so the first support point is far from the closest one.
        let result = sphere_pair_distance([0.0; 3], [0.0, 5.0, 0.0], 1.0, 1.0);
        assert!(!result.intersecting);
        assert!(
            (result.distance - 3.0).abs() < DEFAULT_ACCURACY,
            "expected 3, got {}",
            result.distance
        );
        let gap = sub3_arr(result.point_b, result.point_a);
        assert!(
            (len3_arr(gap) - 3.0).abs() < DEFAULT_ACCURACY && gap[1] > 2.99,
            "witness gap {gap:?}"
        );
        // The distance is an upper bound that tightens with the criteria.
        let tight = sphere_pair_distance_with(
            [0.0; 3],
            [0.0, 5.0, 0.0],
            1.0,
            1.0,
            &GjkTermination::tight(),
        );
        assert!(
            (tight.distance - 3.0).abs() < 1e-8,
            "tight criteria gave {}",
            tight.distance
        );
        assert!(tight.distance <= result.distance + 1e-12);
    }

    #[test]
    fn spheres_offset_along_the_diagonal_report_the_true_gap() {
        let result = sphere_pair_distance([0.0; 3], [3.0, 3.0, 3.0], 1.0, 0.5);
        let expected = 27.0_f64.sqrt() - 1.5;
        assert!(
            (result.distance - expected).abs() < DEFAULT_ACCURACY,
            "expected {expected}, got {}",
            result.distance
        );
        let tight = sphere_pair_distance_with(
            [0.0; 3],
            [3.0, 3.0, 3.0],
            1.0,
            0.5,
            &GjkTermination::tight(),
        );
        assert!(
            (tight.distance - expected).abs() < 1e-8,
            "tight criteria gave {}",
            tight.distance
        );
    }

    #[test]
    fn random_sphere_pairs_match_the_closed_form_distance() {
        let mut state = 0x2545_F491_4F6C_DD1D_u64;
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state >> 11) as f64 / (1u64 << 53) as f64
        };
        for _ in 0..200 {
            let c2 = [
                next() * 12.0 - 6.0,
                next() * 12.0 - 6.0,
                next() * 12.0 - 6.0,
            ];
            let (r1, r2) = (0.2 + next(), 0.2 + next());
            let gap = len3_arr(c2) - r1 - r2;
            if gap < 1e-3 {
                continue;
            }
            let result = sphere_pair_distance([0.0; 3], c2, r1, r2);
            assert!(!result.intersecting);
            assert!(
                (result.distance - gap).abs() < DEFAULT_ACCURACY * gap,
                "gap {gap}, got {}",
                result.distance
            );
            // The reported distance never undercuts the true gap.
            assert!(
                result.distance > gap - 1e-9,
                "gap {gap}, got {}",
                result.distance
            );
        }
    }

    #[test]
    fn separated_boxes_match_the_axis_gaps() {
        let mut cache = GjkCache::new();
        let termination = GjkTermination::default_criteria();
        let half = [0.5, 0.5, 0.5];
        let c2 = [3.0, 4.0, 0.0];
        let result = gjk_distance(
            &mut |dir: [f64; 3]| {
                CachedSupport::new(
                    box_sup([0.0; 3], half, dir),
                    box_sup(c2, half, negate3(dir)),
                )
            },
            &mut cache,
            &termination,
        );
        // Gaps of 2 (x), 3 (y) and 0 (z) between the boxes.
        assert!(!result.intersecting);
        assert!(
            (result.distance - 13.0_f64.sqrt()).abs() < 1e-9,
            "got {}",
            result.distance
        );
    }

    #[test]
    fn overlapping_boxes_are_reported_as_intersecting() {
        let mut cache = GjkCache::new();
        let termination = GjkTermination::default_criteria();
        let half = [1.0, 1.0, 1.0];
        let c2 = [0.7, 0.4, -0.2];
        let result = gjk_distance(
            &mut |dir: [f64; 3]| {
                CachedSupport::new(
                    box_sup([0.0; 3], half, dir),
                    box_sup(c2, half, negate3(dir)),
                )
            },
            &mut cache,
            &termination,
        );
        assert!(result.intersecting, "{result:?}");
        assert_eq!(result.distance, 0.0);
    }

    #[test]
    fn warm_started_query_returns_the_same_distance() {
        let mut cache = GjkCache::new();
        let termination = GjkTermination::default_criteria();
        let c2 = [0.0, 0.0, 4.0];
        let mut support = |dir: [f64; 3]| {
            CachedSupport::new(
                sphere_sup([0.0; 3], 0.75, dir),
                sphere_sup(c2, 0.5, negate3(dir)),
            )
        };
        let cold = gjk_distance(&mut support, &mut cache, &termination);
        let warm = gjk_distance(&mut support, &mut cache, &termination);
        assert!(!cold.warm_started && warm.warm_started);
        assert!(
            (cold.distance - 2.75).abs() < DEFAULT_ACCURACY,
            "cold {}",
            cold.distance
        );
        assert!(
            (warm.distance - 2.75).abs() < DEFAULT_ACCURACY,
            "warm {}",
            warm.distance
        );
        assert!(
            (warm.distance - cold.distance).abs() < DEFAULT_ACCURACY,
            "warm {} vs cold {}",
            warm.distance,
            cold.distance
        );
        assert!(warm.iterations <= cold.iterations);
    }

    #[test]
    fn proximity_and_registry_queries_inherit_the_fixed_distance() {
        let mut registry = GjkCacheRegistry::new();
        let c1 = [0.0; 3];
        let c2 = [0.0, 5.0, 0.0];
        for _ in 0..2 {
            let result = registry.query(1, 2, &mut |dir| sphere_sup(c1, 1.0, dir), &mut |dir| {
                sphere_sup(c2, 1.0, dir)
            });
            assert!(!result.overlapping);
            assert!(
                (result.distance - 3.0).abs() < DEFAULT_ACCURACY,
                "{}",
                result.distance
            );
            assert!(result.normal[1].abs() > 0.999, "{:?}", result.normal);
        }
    }

    #[test]
    fn pruning_keeps_exactly_the_vertices_that_carry_the_closest_point() {
        let make = |x: f64| CachedSupport::new([x, 0.0, 0.0], [0.0; 3]);
        let mut simplex = vec![make(1.0), make(2.0), make(3.0), make(4.0)];
        prune_simplex_to_active(&mut simplex, &[0.5, 0.0, 0.5, 0.0], 2);
        let xs: Vec<f64> = simplex.iter().map(|s| s.minkowski_point[0]).collect();
        assert_eq!(xs, vec![1.0, 3.0]);
        // Numerically vanishing weights never empty the simplex.
        let mut simplex = vec![make(1.0), make(2.0)];
        prune_simplex_to_active(&mut simplex, &[0.0, 0.0, 0.0, 0.0], 2);
        assert_eq!(simplex.len(), 1);
    }

    #[test]
    fn johnson_reports_the_vertices_of_the_closest_feature() {
        // The closest point of this triangle to the origin lies on the edge ab.
        let a = CachedSupport::new([1.0, -1.0, 1.0], [0.0; 3]);
        let b = CachedSupport::new([1.0, 1.0, 1.0], [0.0; 3]);
        let c = CachedSupport::new([4.0, 0.0, 4.0], [0.0; 3]);
        let result = johnson_closest_point(&[a, b, c]);
        assert_eq!(result.num_active, 2);
        assert!(result.bary[2] < 1e-12);
        let edge_only = johnson_closest_point(&[a, b]);
        assert_eq!(edge_only.num_active, 2);
    }

    /// Closest points of two spheres through [`GjkCache::closest_points_cached`].
    fn sphere_closest_points(
        c1: [f64; 3],
        r1: f64,
        c2: [f64; 3],
        r2: f64,
    ) -> Option<([f64; 3], [f64; 3])> {
        let mut cache = GjkCache::new();
        cache.closest_points_cached(|dir| {
            CachedSupport::new(sphere_sup(c1, r1, dir), sphere_sup(c2, r2, negate3(dir)))
        })
    }

    #[test]
    fn closest_points_of_separated_spheres_are_on_the_line_of_centres() {
        let (pa, pb) = sphere_closest_points([0.0; 3], 1.0, [5.0, 0.0, 0.0], 1.0)
            .expect("spheres 3 apart are separated");
        assert!(len3_arr(sub3_arr(pa, [1.0, 0.0, 0.0])) < 1e-9, "{pa:?}");
        assert!(len3_arr(sub3_arr(pb, [4.0, 0.0, 0.0])) < 1e-9, "{pb:?}");
    }

    #[test]
    fn closest_points_converge_off_the_initial_search_direction() {
        // The first search direction (1, 0, 0) is perpendicular to the line of
        // centres, so the first support points are far from the closest ones.
        let (pa, pb) = sphere_closest_points([0.0; 3], 1.0, [0.0, 5.0, 0.0], 1.0)
            .expect("spheres 3 apart are separated");
        assert!((len3_arr(sub3_arr(pb, pa)) - 3.0).abs() < 1e-8);
        assert!(len3_arr(sub3_arr(pa, [0.0, 1.0, 0.0])) < 1e-3, "{pa:?}");
        assert!(len3_arr(sub3_arr(pb, [0.0, 4.0, 0.0])) < 1e-3, "{pb:?}");

        let centre = [3.0, 3.0, 3.0];
        let (pa, pb) = sphere_closest_points([0.0; 3], 1.0, centre, 0.5)
            .expect("spheres on the diagonal are separated");
        let expected = len3_arr(centre) - 1.5;
        assert!((len3_arr(sub3_arr(pb, pa)) - expected).abs() < 1e-8);
    }

    #[test]
    fn overlapping_and_touching_shapes_have_no_closest_pair() {
        // Overlapping spheres: the old inverted termination test reported them
        // as separated.
        assert!(sphere_closest_points([0.0; 3], 1.0, [0.5, 0.0, 0.0], 1.0).is_none());
        assert!(sphere_closest_points([0.0; 3], 1.0, [0.0, 1.5, 0.0], 1.0).is_none());
        // Concentric spheres and exactly touching spheres.
        assert!(sphere_closest_points([0.0; 3], 1.0, [0.0; 3], 0.5).is_none());
        assert!(sphere_closest_points([0.0; 3], 1.0, [2.0, 0.0, 0.0], 1.0).is_none());
        // Overlapping boxes.
        let mut cache = GjkCache::new();
        let half = [1.0, 1.0, 1.0];
        let c2 = [0.7, 0.4, -0.2];
        let boxes = cache.closest_points_cached(|dir| {
            CachedSupport::new(
                box_sup([0.0; 3], half, dir),
                box_sup(c2, half, negate3(dir)),
            )
        });
        assert!(boxes.is_none(), "{boxes:?}");
    }

    #[test]
    fn closest_points_of_point_shapes_are_finite() {
        // Two single points make every support query return the same Minkowski
        // vertex, the degenerate simplex that used to divide 0 by 0.
        let mut cache = GjkCache::new();
        let (pa, pb) = cache
            .closest_points_cached(|_dir| CachedSupport::new([0.0; 3], [3.0, 0.0, 0.0]))
            .expect("distinct points are separated");
        assert_eq!(pa, [0.0; 3]);
        assert_eq!(pb, [3.0, 0.0, 0.0]);
        // Coincident points touch.
        let touching = cache.closest_points_cached(|_dir| CachedSupport::new([1.0; 3], [1.0; 3]));
        assert!(touching.is_none());
    }

    #[test]
    fn non_finite_support_points_have_no_closest_pair() {
        let mut cache = GjkCache::new();
        let result =
            cache.closest_points_cached(|_dir| CachedSupport::new([f64::NAN; 3], [0.0; 3]));
        assert!(result.is_none(), "{result:?}");
        // The garbage simplex must not become the next query's warm start.
        assert!(cache.simplex_vertices.is_empty());
    }

    #[test]
    fn closest_points_keep_the_warm_start_state() {
        let c2 = [0.0, 0.0, 4.0];
        let mut cache = GjkCache::new();
        let mut support = |dir: [f64; 3]| {
            CachedSupport::new(
                sphere_sup([0.0; 3], 0.75, dir),
                sphere_sup(c2, 0.5, negate3(dir)),
            )
        };
        let cold = cache
            .closest_points_cached(&mut support)
            .expect("separated");
        assert!(!cache.simplex_vertices.is_empty());
        // The stored direction points from A towards B.
        assert!(
            cache.last_direction[2] > 0.999,
            "{:?}",
            cache.last_direction
        );
        let warm = cache
            .closest_points_cached(&mut support)
            .expect("separated");
        // Both runs converge to the same distance; on smooth shapes the witnesses
        // are only determined to about the square root of the distance tolerance.
        let cold_gap = len3_arr(sub3_arr(cold.1, cold.0));
        let warm_gap = len3_arr(sub3_arr(warm.1, warm.0));
        assert!((cold_gap - 2.75).abs() < 1e-8, "cold {cold_gap}");
        assert!((warm_gap - 2.75).abs() < 1e-8, "warm {warm_gap}");
        assert!(len3_arr(sub3_arr(warm.0, cold.0)) < 1e-3);
        assert!(len3_arr(sub3_arr(warm.1, cold.1)) < 1e-3);
        assert_eq!(cache.stats.queries, 2);
        assert_eq!(cache.stats.warm_started, 1);
    }

    /// Exact reference shape for the randomised comparison.
    #[derive(Clone, Copy, Debug)]
    enum Reference {
        Sphere { centre: [f64; 3], radius: f64 },
        Cuboid { centre: [f64; 3], half: [f64; 3] },
    }

    impl Reference {
        fn support(&self, dir: [f64; 3]) -> [f64; 3] {
            match *self {
                Reference::Sphere { centre, radius } => sphere_sup(centre, radius, dir),
                Reference::Cuboid { centre, half } => box_sup(centre, half, dir),
            }
        }

        /// Distance from `p` to the shape (0 inside).
        fn distance_to_point(&self, p: [f64; 3]) -> f64 {
            match *self {
                Reference::Sphere { centre, radius } => {
                    (len3_arr(sub3_arr(p, centre)) - radius).max(0.0)
                }
                Reference::Cuboid { centre, half } => {
                    let mut outside = [0.0; 3];
                    for axis in 0..3 {
                        outside[axis] = ((p[axis] - centre[axis]).abs() - half[axis]).max(0.0);
                    }
                    len3_arr(outside)
                }
            }
        }
    }

    /// Exact distance between two reference shapes.
    fn reference_distance(a: Reference, b: Reference) -> f64 {
        match (a, b) {
            (
                Reference::Sphere {
                    centre: ca,
                    radius: ra,
                },
                Reference::Sphere {
                    centre: cb,
                    radius: rb,
                },
            ) => len3_arr(sub3_arr(ca, cb)) - ra - rb,
            (Reference::Sphere { centre, radius }, cuboid @ Reference::Cuboid { .. })
            | (cuboid @ Reference::Cuboid { .. }, Reference::Sphere { centre, radius }) => {
                cuboid.distance_to_point(centre) - radius
            }
            (
                Reference::Cuboid {
                    centre: ca,
                    half: ha,
                },
                Reference::Cuboid {
                    centre: cb,
                    half: hb,
                },
            ) => {
                let mut gap = [0.0; 3];
                for axis in 0..3 {
                    gap[axis] = ((ca[axis] - cb[axis]).abs() - ha[axis] - hb[axis]).max(0.0);
                }
                len3_arr(gap)
            }
        }
    }

    #[test]
    fn random_sphere_and_box_pairs_match_the_exact_distance_and_certificate() {
        let mut state = 0x9E37_79B9_7F4A_7C15_u64;
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state >> 11) as f64 / (1u64 << 53) as f64
        };
        let mut separated = 0;
        for case in 0..240 {
            let mut shape = |kind: usize| {
                let centre = [next() * 8.0 - 4.0, next() * 8.0 - 4.0, next() * 8.0 - 4.0];
                if kind.is_multiple_of(2) {
                    Reference::Sphere {
                        centre,
                        radius: 0.3 + next(),
                    }
                } else {
                    Reference::Cuboid {
                        centre,
                        half: [0.3 + next(), 0.3 + next(), 0.3 + next()],
                    }
                }
            };
            let a = shape(case);
            let b = shape(case / 2);
            let exact = reference_distance(a, b);
            if exact.abs() < 1e-6 {
                continue;
            }
            let mut cache = GjkCache::new();
            let result = cache.closest_points_cached(|dir| {
                CachedSupport::new(a.support(dir), b.support(negate3(dir)))
            });
            if exact < 0.0 {
                assert!(result.is_none(), "overlap {exact} reported as {result:?}");
                continue;
            }
            let (pa, pb) =
                result.unwrap_or_else(|| panic!("separation {exact} reported as overlap"));
            separated += 1;
            let found = len3_arr(sub3_arr(pb, pa));
            assert!(
                (found - exact).abs() < 1e-7 * (1.0 + exact),
                "exact {exact}, found {found} for {a:?} and {b:?}"
            );
            // Optimality certificate: the planes through the witnesses with the
            // closest-point normal support the two shapes from opposite sides.
            let n = scale3_arr(sub3_arr(pb, pa), 1.0 / found);
            let top_a = dot3_arr(n, a.support(n));
            let bottom_b = dot3_arr(n, b.support(negate3(n)));
            assert!(
                (top_a - dot3_arr(n, pa)).abs() < 1e-7,
                "A overshoots its witness"
            );
            assert!(
                (bottom_b - dot3_arr(n, pb)).abs() < 1e-7,
                "B overshoots its witness"
            );
        }
        assert!(separated >= 100, "only {separated} separated pairs");
    }

    /// The point `sum_i bary[i] * simplex[i]`.
    fn weighted_point(simplex: &[CachedSupport], bary: &[f64; 4]) -> [f64; 3] {
        let mut point = [0.0; 3];
        for (vertex, weight) in simplex.iter().zip(bary) {
            for (coordinate, component) in point.iter_mut().zip(vertex.minkowski_point) {
                *coordinate += weight * component;
            }
        }
        point
    }

    #[test]
    fn sliver_triangle_weights_reproduce_the_closest_point_on_the_far_edge() {
        // Vertices `a` and `c` are almost coincident, so the triangle is treated
        // as degenerate and the closest point is found on its edges. It lies on
        // the edge `bc`, whose weights used to come back swapped.
        let a = CachedSupport::new(
            [
                0.117_311_962_788_504_86,
                -0.149_304_520_925_614_97,
                -3.131_343_814_998_244_7,
            ],
            [0.0; 3],
        );
        let b = CachedSupport::new(
            [
                -0.962_663_206_898_577_9,
                -0.149_304_520_897_032_5,
                -3.131_343_814_398_788_6,
            ],
            [0.0; 3],
        );
        let c = CachedSupport::new(
            [
                0.117_283_365_568_119_36,
                -0.149_304_520_897_032_5,
                -3.131_343_814_398_788_6,
            ],
            [0.0; 3],
        );
        let simplex = [a, b, c];
        let result = johnson_closest_point(&simplex);
        assert!(!result.origin_inside);
        assert!(result.bary[0].abs() < 1e-12, "{:?}", result.bary);
        let reproduced = weighted_point(&simplex, &result.bary);
        assert!(
            len3_arr(sub3_arr(reproduced, result.closest)) < 1e-12,
            "weights {:?} give {reproduced:?}, closest is {:?}",
            result.bary,
            result.closest
        );
        // The closest point is where the edge crosses x = 0.
        assert!(result.closest[0].abs() < 1e-12, "{:?}", result.closest);
        let (pa, pb) = barycentric_witness(&simplex, &result.bary);
        assert!(len3_arr(sub3_arr(sub3_arr(pa, pb), result.closest)) < 1e-12);
    }

    #[test]
    fn johnson_weights_reproduce_the_closest_point_on_random_simplices() {
        let mut state = 0xD1B5_4A32_D192_ED03_u64;
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state >> 11) as f64 / (1u64 << 53) as f64
        };
        let mut checked = 0;
        for case in 0..600 {
            let count = 1 + case % 4;
            let mut points: Vec<[f64; 3]> = (0..count)
                .map(|_| {
                    [
                        next() * 6.0 - 3.0,
                        next() * 6.0 - 3.0,
                        next() * 6.0 - 3.0 + 3.5,
                    ]
                })
                .collect();
            // Every third simplex gets slivers: a vertex pushed next to another.
            if case % 3 == 0 && count >= 3 {
                let tiny = 10f64.powf(-3.0 - 3.0 * next());
                points[2] = [
                    points[0][0] + tiny * (next() - 0.5),
                    points[0][1] + tiny * (next() - 0.5),
                    points[0][2] + tiny * (next() - 0.5),
                ];
            }
            let simplex: Vec<CachedSupport> = points
                .iter()
                .map(|p| CachedSupport::new(*p, [0.0; 3]))
                .collect();
            let result = johnson_closest_point(&simplex);
            if result.origin_inside {
                continue;
            }
            checked += 1;
            let active: f64 = result.bary[..count].iter().sum();
            assert!((active - 1.0).abs() < 1e-9, "weights {:?}", result.bary);
            assert!(result.bary[..count].iter().all(|w| *w >= -1e-12));
            assert!(result.bary[count..].iter().all(|w| *w == 0.0));
            let reproduced = weighted_point(&simplex, &result.bary);
            assert!(
                len3_arr(sub3_arr(reproduced, result.closest)) < 1e-9,
                "case {case}: weights {:?} give {reproduced:?}, closest {:?} of {points:?}",
                result.bary,
                result.closest
            );
            // Never farther than a vertex, and never farther than the midpoint of
            // any edge (up to the conditioning of slivers).
            for (i, p) in points.iter().enumerate() {
                assert!(result.distance() <= len3_arr(*p) + 1e-7);
                for q in &points[i + 1..] {
                    let mid = scale3_arr([p[0] + q[0], p[1] + q[1], p[2] + q[2]], 0.5);
                    assert!(result.distance() <= len3_arr(mid) + 1e-7);
                }
            }
        }
        assert!(checked > 500, "only {checked} simplices checked");
    }
}
