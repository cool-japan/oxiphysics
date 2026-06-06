//! SAP tests (part 3)
//!
//! Continuation of functions_2.rs — split to comply with the 2000-line policy.

#[cfg(test)]
mod tests {
    use super::super::functions::*;
    use super::super::types::*;
    use std::collections::HashMap;
    use std::collections::HashSet;
    #[test]
    fn sap_add_object_batch_empty_batch() {
        let mut sap = SweepAndPrune::new();
        sap.add_object_batch(&[]);
        assert_eq!(sap.object_count(), 0, "empty batch should leave SAP empty");
    }
    #[test]
    fn sap_compute_axis_variance_empty() {
        let sap = SweepAndPrune::new();
        assert_eq!(
            sap.compute_axis_variance(),
            [0.0; 3],
            "empty SAP variance is all zeros"
        );
    }
    #[test]
    fn sap_compute_axis_variance_x_dominates() {
        let mut sap = SweepAndPrune::new();
        for i in 0u64..10 {
            let x = i as f64 * 100.0;
            sap.add_object(i, [x, 0.0, 0.0], [x + 1.0, 1.0, 1.0]);
        }
        let var = sap.compute_axis_variance();
        assert!(var[0] > var[1], "X variance should exceed Y variance");
        assert!(var[0] > var[2], "X variance should exceed Z variance");
    }
    #[test]
    fn sap_compute_axis_variance_single_object() {
        let mut sap = SweepAndPrune::new();
        sap.add_object(1, [0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let var = sap.compute_axis_variance();
        assert_eq!(var, [0.0; 3]);
    }
    #[test]
    fn sap_reorder_axes_returns_x_when_x_dominates() {
        let mut sap = SweepAndPrune::new();
        for i in 0u64..5 {
            let x = i as f64 * 50.0;
            sap.add_object(i, [x, 0.0, 0.0], [x + 1.0, 1.0, 1.0]);
        }
        let axis = sap.reorder_axes();
        assert_eq!(axis, 0, "X axis has highest variance, should be chosen");
    }
    #[test]
    fn sap_reorder_axes_sorted_after_reorder() {
        let mut sap = SweepAndPrune::new();
        for i in [5u64, 2, 8, 1, 9] {
            let x = i as f64 * 10.0;
            sap.add_object(i, [x, 0.0, 0.0], [x + 2.0, 1.0, 1.0]);
        }
        sap.reorder_axes();
        for i in 1..sap.sorted_x.len() {
            assert!(
                sap.sorted_x[i - 1].value <= sap.sorted_x[i].value,
                "endpoints must be sorted after reorder_axes"
            );
        }
    }
    #[test]
    fn incremental_sap_compute_axis_variance_y_dominates() {
        let mut sap = IncrementalSap::new();
        for i in 0u32..8 {
            let y = i as f64 * 100.0;
            sap.insert(
                i,
                Aabb3 {
                    min: [0.0, y, 0.0],
                    max: [1.0, y + 1.0, 1.0],
                },
            );
        }
        let var = sap.compute_axis_variance();
        assert!(var[1] > var[0], "Y variance should exceed X");
        assert!(var[1] > var[2], "Y variance should exceed Z");
    }
    #[test]
    fn incremental_sap_reorder_axes_empty() {
        let mut sap = IncrementalSap::new();
        let axis = sap.reorder_axes();
        assert!(axis < 3, "axis must be in 0..3 for empty SAP");
    }
    #[test]
    fn axis_is_sorted_true_for_ascending() {
        let eps = vec![
            SapEndpointU32 {
                value: 1.0,
                body_id: 1,
                is_min: true,
            },
            SapEndpointU32 {
                value: 2.0,
                body_id: 1,
                is_min: false,
            },
            SapEndpointU32 {
                value: 3.0,
                body_id: 2,
                is_min: true,
            },
        ];
        assert!(axis_is_sorted(&eps));
    }
    #[test]
    fn axis_is_sorted_false_for_unsorted() {
        let eps = vec![
            SapEndpointU32 {
                value: 3.0,
                body_id: 1,
                is_min: true,
            },
            SapEndpointU32 {
                value: 1.0,
                body_id: 2,
                is_min: true,
            },
        ];
        assert!(!axis_is_sorted(&eps));
    }
    #[test]
    fn axis_is_sorted_empty_is_sorted() {
        assert!(axis_is_sorted(&[]));
    }
    #[test]
    fn endpoint_range_finds_all_in_window() {
        let eps = vec![
            SapEndpointU32 {
                value: 1.0,
                body_id: 1,
                is_min: true,
            },
            SapEndpointU32 {
                value: 2.0,
                body_id: 1,
                is_min: false,
            },
            SapEndpointU32 {
                value: 3.0,
                body_id: 2,
                is_min: true,
            },
            SapEndpointU32 {
                value: 4.0,
                body_id: 2,
                is_min: false,
            },
        ];
        let (lo, hi) = endpoint_range(&eps, 1.5, 3.5);
        assert_eq!(lo, 1, "lower bound should be index 1");
        assert_eq!(hi, 3, "upper bound should be index 3");
    }
    #[test]
    fn endpoint_range_empty_window() {
        let eps = vec![
            SapEndpointU32 {
                value: 1.0,
                body_id: 1,
                is_min: true,
            },
            SapEndpointU32 {
                value: 2.0,
                body_id: 1,
                is_min: false,
            },
        ];
        let (lo, hi) = endpoint_range(&eps, 5.0, 6.0);
        assert_eq!(lo, hi, "no endpoints in [5,6] — should give empty range");
    }
    #[test]
    fn bipartite_sap_query_cross_set_pairs_only() {
        let mut aabbs = HashMap::new();
        aabbs.insert(
            1u32,
            Aabb3 {
                min: [0.0; 3],
                max: [3.0; 3],
            },
        );
        aabbs.insert(
            2u32,
            Aabb3 {
                min: [1.0; 3],
                max: [4.0; 3],
            },
        );
        aabbs.insert(
            3u32,
            Aabb3 {
                min: [2.0; 3],
                max: [5.0; 3],
            },
        );
        let pairs = bipartite_sap_query(&[1, 2], &[3], &aabbs);
        for &(a, b) in &pairs {
            let cross = (a == 3 || b == 3) && (a != b);
            assert!(cross, "unexpected intra-set pair ({},{})", a, b);
        }
        assert!(!pairs.is_empty(), "should find at least one cross-set pair");
    }
    #[test]
    fn bipartite_sap_query_no_overlap_returns_empty() {
        let mut aabbs = HashMap::new();
        aabbs.insert(
            1u32,
            Aabb3 {
                min: [0.0; 3],
                max: [1.0; 3],
            },
        );
        aabbs.insert(
            2u32,
            Aabb3 {
                min: [10.0; 3],
                max: [11.0; 3],
            },
        );
        let pairs = bipartite_sap_query(&[1], &[2], &aabbs);
        assert!(pairs.is_empty());
    }
    #[test]
    fn sorted_insert_keeps_order() {
        let mut eps = vec![
            SapEndpointU32 {
                value: 1.0,
                body_id: 1,
                is_min: true,
            },
            SapEndpointU32 {
                value: 3.0,
                body_id: 2,
                is_min: true,
            },
        ];
        sorted_insert(
            &mut eps,
            SapEndpointU32 {
                value: 2.0,
                body_id: 3,
                is_min: true,
            },
        );
        assert_eq!(eps.len(), 3);
        assert!(
            axis_is_sorted(&eps),
            "list must remain sorted after sorted_insert"
        );
        assert_eq!(eps[1].value, 2.0, "inserted value should be at index 1");
    }
    #[test]
    fn sorted_insert_at_head() {
        let mut eps = vec![SapEndpointU32 {
            value: 5.0,
            body_id: 1,
            is_min: true,
        }];
        let pos = sorted_insert(
            &mut eps,
            SapEndpointU32 {
                value: 1.0,
                body_id: 2,
                is_min: true,
            },
        );
        assert_eq!(pos, 0, "should be inserted at head");
        assert_eq!(eps[0].value, 1.0);
    }
    #[test]
    fn aabb3_overlaps_true() {
        let a = Aabb3 {
            min: [0.0; 3],
            max: [2.0; 3],
        };
        let b = Aabb3 {
            min: [1.0; 3],
            max: [3.0; 3],
        };
        assert!(aabb3_overlaps(&a, &b));
    }
    #[test]
    fn aabb3_overlaps_false_separated() {
        let a = Aabb3 {
            min: [0.0; 3],
            max: [1.0; 3],
        };
        let b = Aabb3 {
            min: [2.0; 3],
            max: [3.0; 3],
        };
        assert!(!aabb3_overlaps(&a, &b));
    }
    #[test]
    fn aabb3_overlaps_touching_boundary() {
        let a = Aabb3 {
            min: [0.0; 3],
            max: [1.0; 3],
        };
        let b = Aabb3 {
            min: [1.0; 3],
            max: [2.0; 3],
        };
        assert!(
            aabb3_overlaps(&a, &b),
            "touching at boundary should count as overlap"
        );
    }
    #[test]
    fn aabb3_expand_to_include_basic() {
        let mut a = Aabb3 {
            min: [0.0; 3],
            max: [1.0; 3],
        };
        let b = Aabb3 {
            min: [0.5; 3],
            max: [2.0; 3],
        };
        aabb3_expand_to_include(&mut a, &b);
        assert_eq!(a.min, [0.0; 3]);
        assert_eq!(a.max, [2.0; 3]);
    }
    #[test]
    fn aabb3_expand_to_include_contained_is_noop() {
        let mut a = Aabb3 {
            min: [0.0; 3],
            max: [4.0; 3],
        };
        let b = Aabb3 {
            min: [1.0; 3],
            max: [3.0; 3],
        };
        aabb3_expand_to_include(&mut a, &b);
        assert_eq!(a.min, [0.0; 3]);
        assert_eq!(a.max, [4.0; 3]);
    }
    #[test]
    fn aabb3_point_dist_sq_inside_is_zero() {
        let a = Aabb3 {
            min: [0.0; 3],
            max: [2.0; 3],
        };
        assert_eq!(aabb3_point_dist_sq(&a, [1.0, 1.0, 1.0]), 0.0);
    }
    #[test]
    fn aabb3_point_dist_sq_outside() {
        let a = Aabb3 {
            min: [0.0; 3],
            max: [1.0; 3],
        };
        let d = aabb3_point_dist_sq(&a, [2.0, 0.0, 0.0]);
        assert!((d - 1.0).abs() < 1e-10);
    }
    #[test]
    fn isap_any_overlap_true() {
        let mut sap = IncrementalSap::new();
        sap.insert(
            1,
            Aabb3 {
                min: [0.0; 3],
                max: [2.0; 3],
            },
        );
        let q = Aabb3 {
            min: [1.0; 3],
            max: [3.0; 3],
        };
        assert!(sap.any_overlap(&q));
    }
    #[test]
    fn isap_any_overlap_false() {
        let mut sap = IncrementalSap::new();
        sap.insert(
            1,
            Aabb3 {
                min: [0.0; 3],
                max: [1.0; 3],
            },
        );
        let q = Aabb3 {
            min: [5.0; 3],
            max: [6.0; 3],
        };
        assert!(!sap.any_overlap(&q));
    }
    #[test]
    fn isap_query_aabb_returns_overlapping() {
        let mut sap = IncrementalSap::new();
        sap.insert(
            1,
            Aabb3 {
                min: [0.0; 3],
                max: [2.0; 3],
            },
        );
        sap.insert(
            2,
            Aabb3 {
                min: [10.0; 3],
                max: [11.0; 3],
            },
        );
        let hits = sap.query_aabb(&Aabb3 {
            min: [1.0; 3],
            max: [3.0; 3],
        });
        assert!(hits.contains(&1), "body 1 should be in query result");
        assert!(!hits.contains(&2), "body 2 is far away");
    }
    #[test]
    fn isap_query_sphere_sq_basic() {
        let mut sap = IncrementalSap::new();
        sap.insert(
            1,
            Aabb3 {
                min: [0.0; 3],
                max: [1.0; 3],
            },
        );
        sap.insert(
            2,
            Aabb3 {
                min: [10.0; 3],
                max: [11.0; 3],
            },
        );
        let ids = sap.query_sphere_sq([0.5, 0.5, 0.5], 4.0);
        assert!(ids.contains(&1), "body 1 is close enough");
        assert!(!ids.contains(&2), "body 2 is far");
    }
    #[test]
    fn isap_clear_empties_everything() {
        let mut sap = IncrementalSap::new();
        sap.insert(
            1,
            Aabb3 {
                min: [0.0; 3],
                max: [1.0; 3],
            },
        );
        sap.insert(
            2,
            Aabb3 {
                min: [0.5; 3],
                max: [1.5; 3],
            },
        );
        sap.query_pairs();
        sap.clear();
        assert_eq!(sap.body_count(), 0);
        assert!(sap.active_pairs.is_empty());
        assert!(sap.endpoints_x.is_empty());
    }
    #[test]
    fn isap_get_aabb_found_and_missing() {
        let mut sap = IncrementalSap::new();
        sap.insert(
            7,
            Aabb3 {
                min: [1.0; 3],
                max: [2.0; 3],
            },
        );
        assert!(sap.get_aabb(7).is_some());
        assert!(sap.get_aabb(999).is_none());
    }
    #[test]
    fn isap_contains_after_insert_remove() {
        let mut sap = IncrementalSap::new();
        assert!(!sap.contains(1));
        sap.insert(
            1,
            Aabb3 {
                min: [0.0; 3],
                max: [1.0; 3],
            },
        );
        assert!(sap.contains(1));
        sap.remove(1);
        assert!(!sap.contains(1));
    }
    #[test]
    fn isap_merge_from_adds_bodies() {
        let mut sap_a = IncrementalSap::new();
        sap_a.insert(
            1,
            Aabb3 {
                min: [0.0; 3],
                max: [1.0; 3],
            },
        );
        let mut sap_b = IncrementalSap::new();
        sap_b.insert(
            2,
            Aabb3 {
                min: [0.5; 3],
                max: [1.5; 3],
            },
        );
        sap_b.insert(
            3,
            Aabb3 {
                min: [5.0; 3],
                max: [6.0; 3],
            },
        );
        sap_a.merge_from(&sap_b);
        assert_eq!(sap_a.body_count(), 3, "merged sap should have 3 bodies");
        assert!(sap_a.contains(2));
        assert!(sap_a.contains(3));
    }
    #[test]
    fn sap_snapshot_restore_reverts_additions() {
        let mut sap = IncrementalSap::new();
        sap.insert(
            1,
            Aabb3 {
                min: [0.0; 3],
                max: [1.0; 3],
            },
        );
        let snap = SapSnapshot::capture(&sap);
        sap.insert(
            2,
            Aabb3 {
                min: [5.0; 3],
                max: [6.0; 3],
            },
        );
        assert_eq!(sap.body_count(), 2);
        snap.restore(&mut sap);
        assert_eq!(sap.body_count(), 1, "snapshot restore should remove body 2");
        assert!(sap.contains(1));
        assert!(!sap.contains(2));
    }
    #[test]
    fn sap_snapshot_restore_reverts_removals() {
        let mut sap = IncrementalSap::new();
        sap.insert(
            1,
            Aabb3 {
                min: [0.0; 3],
                max: [1.0; 3],
            },
        );
        sap.insert(
            2,
            Aabb3 {
                min: [2.0; 3],
                max: [3.0; 3],
            },
        );
        let snap = SapSnapshot::capture(&sap);
        sap.remove(2);
        assert_eq!(sap.body_count(), 1);
        snap.restore(&mut sap);
        assert_eq!(sap.body_count(), 2, "snapshot restore should re-add body 2");
        assert!(sap.contains(2));
    }
    #[test]
    fn pair_delta_detects_new_and_removed() {
        let prev: HashSet<(u32, u32)> = [(1, 2), (3, 4)].iter().copied().collect();
        let curr: HashSet<(u32, u32)> = [(1, 2), (5, 6)].iter().copied().collect();
        let (new, removed) = pair_delta(&prev, &curr);
        assert_eq!(new, vec![(5, 6)]);
        assert_eq!(removed, vec![(3, 4)]);
    }
    #[test]
    fn pair_delta_empty_both() {
        let prev: HashSet<(u32, u32)> = HashSet::new();
        let curr: HashSet<(u32, u32)> = HashSet::new();
        let (new, removed) = pair_delta(&prev, &curr);
        assert!(new.is_empty());
        assert!(removed.is_empty());
    }
    #[test]
    fn aabb3_clamp_inside_world_unchanged() {
        let a = Aabb3 {
            min: [1.0; 3],
            max: [3.0; 3],
        };
        let c = aabb3_clamp(&a, [0.0; 3], [10.0; 3]);
        assert_eq!(c.min, [1.0; 3]);
        assert_eq!(c.max, [3.0; 3]);
    }
    #[test]
    fn aabb3_clamp_outside_world_clamped() {
        let a = Aabb3 {
            min: [-5.0; 3],
            max: [15.0; 3],
        };
        let c = aabb3_clamp(&a, [0.0; 3], [10.0; 3]);
        assert_eq!(c.min, [0.0; 3]);
        assert_eq!(c.max, [10.0; 3]);
    }
    #[test]
    fn grid_query_point_finds_objects_in_cell() {
        let mut grid = GridBroadphase::new(5.0);
        grid.insert(1, [0.0; 3], [4.0; 3]);
        grid.insert(2, [6.0; 3], [9.0; 3]);
        let hits = grid.query_point([1.0, 1.0, 1.0]);
        assert!(hits.contains(&1), "object 1 should be in cell (0,0,0)");
        assert!(!hits.contains(&2), "object 2 is in a different cell");
    }
    #[test]
    fn grid_total_entries_counts_all_slots() {
        let mut grid = GridBroadphase::new(10.0);
        grid.insert(1, [0.0; 3], [1.0; 3]);
        grid.insert(2, [0.5; 3], [1.5; 3]);
        assert!(grid.total_entries() >= 2, "at least 2 entries");
    }
    #[test]
    fn grid_occupied_cells_non_zero_after_insert() {
        let mut grid = GridBroadphase::new(10.0);
        grid.insert(1, [0.0; 3], [1.0; 3]);
        assert!(grid.occupied_cells() > 0);
    }
    #[test]
    fn grid_remove_clears_object() {
        let mut grid = GridBroadphase::new(10.0);
        grid.insert(1, [0.0; 3], [1.0; 3]);
        grid.insert(2, [0.5; 3], [1.5; 3]);
        assert!(!grid.query_potential_pairs().is_empty());
        grid.remove(2);
        assert!(
            grid.query_potential_pairs().is_empty(),
            "after removing object 2, no pairs"
        );
    }
    #[test]
    fn grid_query_aabb_finds_overlapping_cells() {
        let mut grid = GridBroadphase::new(5.0);
        grid.insert(10, [0.0; 3], [3.0; 3]);
        grid.insert(20, [20.0; 3], [23.0; 3]);
        let hits = grid.query_aabb([0.0; 3], [3.0; 3]);
        assert!(hits.contains(&10), "object 10 should be found");
        assert!(!hits.contains(&20), "object 20 is far away");
    }
    #[test]
    fn sap_axis_body_count_tracks_insertions() {
        let mut axis = SapAxis::new();
        assert_eq!(axis.body_count(), 0);
        axis.insert(1, 0.0, 1.0);
        assert_eq!(axis.body_count(), 1);
        axis.insert(2, 2.0, 3.0);
        assert_eq!(axis.body_count(), 2);
    }
    #[test]
    fn sap_axis_min_max_values_after_insert() {
        let mut axis = SapAxis::new();
        axis.insert(1, 3.0, 7.0);
        axis.insert(2, 1.0, 5.0);
        assert!(
            (axis.min_value().unwrap() - 1.0).abs() < 1e-12,
            "min should be 1.0"
        );
        assert!(
            (axis.max_value().unwrap() - 7.0).abs() < 1e-12,
            "max should be 7.0"
        );
    }
    #[test]
    fn sap_axis_tracked_bodies_returns_all() {
        let mut axis = SapAxis::new();
        axis.insert(5, 0.0, 1.0);
        axis.insert(3, 2.0, 4.0);
        let bodies = axis.tracked_bodies();
        assert_eq!(bodies, vec![3, 5], "tracked bodies should be sorted");
    }
    #[test]
    fn sap_axis_is_sorted_after_insert() {
        let mut axis = SapAxis::new();
        axis.insert(2, 5.0, 9.0);
        axis.insert(1, 1.0, 3.0);
        axis.insert(3, 2.0, 7.0);
        assert!(axis.is_sorted(), "axis should be sorted after insertions");
    }
    #[test]
    fn sweep_window_query_finds_overlapping_bodies() {
        let mut sap = IncrementalSap::new();
        sap.insert(
            1,
            Aabb3 {
                min: [0.0; 3],
                max: [2.0; 3],
            },
        );
        sap.insert(
            2,
            Aabb3 {
                min: [3.0; 3],
                max: [5.0; 3],
            },
        );
        sap.insert(
            3,
            Aabb3 {
                min: [10.0; 3],
                max: [12.0; 3],
            },
        );
        let mut eps = sap.endpoints_x.clone();
        eps.sort_by(|a, b| {
            a.value
                .partial_cmp(&b.value)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        let hits = sweep_window_query(&eps, 1.0, 4.0);
        assert!(hits.contains(&1), "body 1 overlaps [1,4]");
        assert!(hits.contains(&2), "body 2 overlaps [1,4]");
        assert!(!hits.contains(&3), "body 3 starts at 10, outside window");
    }
    #[test]
    fn sweep_window_query_empty_range() {
        let mut sap = IncrementalSap::new();
        sap.insert(
            1,
            Aabb3 {
                min: [0.0; 3],
                max: [1.0; 3],
            },
        );
        let mut eps = sap.endpoints_x.clone();
        eps.sort_by(|a, b| {
            a.value
                .partial_cmp(&b.value)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        let hits = sweep_window_query(&eps, 5.0, 6.0);
        assert!(hits.is_empty(), "no bodies in [5,6]");
    }
}
