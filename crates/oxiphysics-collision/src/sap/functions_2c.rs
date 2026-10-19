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

    // ── reorder_axes must never change which pairs are reported ─────────────

    /// Unit boxes 1 and 2 overlap on Y and Z but are 4 apart on X; box 3 sits
    /// far away on Y so that Y has by far the largest centre variance.
    fn x_separated_scene() -> [(u64, [f64; 3], [f64; 3]); 3] {
        [
            (1, [0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
            (2, [5.0, 0.0, 0.0], [6.0, 1.0, 1.0]),
            (3, [10.0, 200.0, 0.0], [11.0, 201.0, 1.0]),
        ]
    }

    #[test]
    fn sap_reorder_axes_to_y_still_rejects_x_separated_pair() {
        let mut sap = SweepAndPrune::new();
        sap.add_object_batch(&x_separated_scene());
        assert_eq!(sap.sweep_axis(), 0);
        assert!(sap.query_overlapping_pairs().is_empty());
        let axis = sap.reorder_axes();
        assert_eq!(axis, 1, "Y has the largest variance");
        assert_eq!(sap.sweep_axis(), 1);
        // Boxes 1 and 2 overlap on the new sweep axis (Y) and on Z but not on X;
        // the filter must therefore test X, not Y again.
        assert!(
            sap.query_overlapping_pairs().is_empty(),
            "boxes 4 apart on X must not be reported after reorder_axes"
        );
    }

    #[test]
    fn sap_update_and_add_use_the_reordered_sweep_axis() {
        let mut sap = SweepAndPrune::new();
        sap.add_object_batch(&x_separated_scene());
        assert_eq!(sap.reorder_axes(), 1);
        // Move box 2 so that it overlaps box 1 on every axis.
        sap.update_object(2, [0.5, 0.0, 0.0], [1.5, 1.0, 1.0]);
        assert_eq!(sap.query_overlapping_pairs(), vec![(1, 2)]);
        // Every endpoint value must now come from the Y extents.
        for ep in &sap.sorted_x {
            let obj = sap
                .objects
                .iter()
                .find(|o| o.id == ep.object_id)
                .expect("endpoint belongs to a registered object");
            let expected = if ep.is_min { obj.min[1] } else { obj.max[1] };
            assert_eq!(ep.value, expected, "endpoint of object {}", ep.object_id);
        }
        // A late arrival is inserted with Y endpoints as well and found.
        sap.add_object(4, [0.2, 0.2, 0.2], [0.8, 0.8, 0.8]);
        assert_eq!(sap.query_overlapping_pairs(), vec![(1, 2), (1, 4), (2, 4)]);
    }

    #[test]
    fn sap_pairs_match_brute_force_for_every_sweep_axis() {
        // Stretch one axis so that reorder_axes selects it, then compare the
        // reported pairs with the brute-force AABB test.
        for stretch_axis in 0..3usize {
            let mut state = 0x9E37_79B9_7F4A_7C15u64 ^ stretch_axis as u64;
            let mut next = move || {
                state = state
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                (state >> 33) as f64 / (1u64 << 31) as f64
            };
            let mut objects: Vec<(u64, [f64; 3], [f64; 3])> = Vec::new();
            for id in 0..60u64 {
                let mut min = [next() * 10.0, next() * 10.0, next() * 10.0];
                min[stretch_axis] *= 40.0;
                let size = [0.5 + next() * 2.0, 0.5 + next() * 2.0, 0.5 + next() * 2.0];
                let max = [min[0] + size[0], min[1] + size[1], min[2] + size[2]];
                objects.push((id, min, max));
            }
            let mut expected: Vec<(u64, u64)> = Vec::new();
            for (i, a) in objects.iter().enumerate() {
                for b in &objects[i + 1..] {
                    let overlap = (0..3).all(|k| a.1[k] <= b.2[k] && b.1[k] <= a.2[k]);
                    if overlap {
                        expected.push((a.0.min(b.0), a.0.max(b.0)));
                    }
                }
            }
            expected.sort_unstable();
            let mut sap = SweepAndPrune::new();
            sap.add_object_batch(&objects);
            assert_eq!(sap.query_overlapping_pairs(), expected, "sweep on X");
            assert_eq!(sap.reorder_axes(), stretch_axis, "axis chosen by variance");
            assert_eq!(
                sap.query_overlapping_pairs(),
                expected,
                "sweep on axis {stretch_axis}"
            );
        }
    }

    #[test]
    fn incremental_sap_reorder_axes_keeps_every_axis_in_the_pair_test() {
        let mut sap = IncrementalSap::new();
        for (id, min, max) in x_separated_scene() {
            sap.insert(id as u32, Aabb3 { min, max });
        }
        assert!(sap.query_pairs().is_empty());
        let x_before: Vec<f64> = sap.endpoints_x.iter().map(|e| e.value).collect();
        assert_eq!(sap.reorder_axes(), 1, "Y has the largest variance");
        let x_after: Vec<f64> = sap.endpoints_x.iter().map(|e| e.value).collect();
        assert_eq!(
            x_before, x_after,
            "reorder_axes must not overwrite the X endpoint list with Y values"
        );
        assert!(
            sap.query_pairs().is_empty(),
            "boxes 1 and 2 are separated on X and must stay unpaired"
        );
        // Updates after the reorder still write X coordinates into the X list.
        sap.update(
            2,
            Aabb3 {
                min: [0.5, 0.0, 0.0],
                max: [1.5, 1.0, 1.0],
            },
        );
        assert_eq!(sap.query_pairs(), vec![(1, 2)]);
    }

    // ── lazy sorting: add/update only mark the endpoint list dirty ──────────

    fn brute_force_pairs(objects: &[(u64, [f64; 3], [f64; 3])]) -> Vec<(u64, u64)> {
        let mut expected: Vec<(u64, u64)> = Vec::new();
        for (i, a) in objects.iter().enumerate() {
            for b in &objects[i + 1..] {
                if (0..3).all(|k| a.1[k] <= b.2[k] && b.1[k] <= a.2[k]) {
                    expected.push((a.0.min(b.0), a.0.max(b.0)));
                }
            }
        }
        expected.sort_unstable();
        expected
    }

    fn endpoints_ascending(sap: &SweepAndPrune) -> bool {
        sap.sorted_x
            .windows(2)
            .all(|pair| pair[0].value <= pair[1].value)
    }

    #[test]
    fn sap_add_object_does_not_sort_on_every_call_but_the_query_sorts_once() {
        let mut sap = SweepAndPrune::new();
        // Descending positions: an eager per-call sort would keep the list
        // ascending, the lazy implementation leaves it for the query.
        for i in 0u64..6 {
            let x = (5 - i) as f64 * 10.0;
            sap.add_object(i, [x, 0.0, 0.0], [x + 1.0, 1.0, 1.0]);
        }
        assert!(
            !endpoints_ascending(&sap),
            "add_object must only append the endpoints"
        );
        assert!(sap.query_overlapping_pairs().is_empty());
        assert!(
            endpoints_ascending(&sap),
            "the query sorts the endpoint list"
        );
        // A query on a clean list leaves it untouched and sorted.
        let before: Vec<(u64, bool)> = sap
            .sorted_x
            .iter()
            .map(|e| (e.object_id, e.is_min))
            .collect();
        assert!(sap.query_overlapping_pairs().is_empty());
        let after: Vec<(u64, bool)> = sap
            .sorted_x
            .iter()
            .map(|e| (e.object_id, e.is_min))
            .collect();
        assert_eq!(before, after);
    }

    #[test]
    fn sap_many_updates_between_queries_match_brute_force() {
        let mut state = 0xC0FF_EE00_1234_5678u64;
        let mut next = move || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (state >> 33) as f64 / (1u64 << 31) as f64
        };
        let mut objects: Vec<(u64, [f64; 3], [f64; 3])> = Vec::new();
        for id in 0..40u64 {
            let min = [next() * 30.0, next() * 30.0, next() * 30.0];
            let size = [1.0 + next() * 3.0, 1.0 + next() * 3.0, 1.0 + next() * 3.0];
            objects.push((
                id,
                min,
                [min[0] + size[0], min[1] + size[1], min[2] + size[2]],
            ));
        }
        let mut sap = SweepAndPrune::new();
        sap.add_object_batch(&objects);
        assert_eq!(sap.query_overlapping_pairs(), brute_force_pairs(&objects));
        // Several frames: move every object, then query once per frame.
        for frame in 0..5 {
            if frame == 2 {
                sap.reorder_axes();
            }
            for object in objects.iter_mut() {
                let shift = [next() * 4.0 - 2.0, next() * 4.0 - 2.0, next() * 4.0 - 2.0];
                for (k, delta) in shift.iter().enumerate() {
                    object.1[k] += delta;
                    object.2[k] += delta;
                }
                sap.update_object(object.0, object.1, object.2);
            }
            assert_eq!(
                sap.query_overlapping_pairs(),
                brute_force_pairs(&objects),
                "frame {frame}"
            );
            assert!(endpoints_ascending(&sap), "frame {frame}");
        }
    }

    #[test]
    fn sap_touching_boxes_are_reported_in_every_insertion_order() {
        // Box 1 ends where box 2 begins on every axis.
        let one = (1u64, [0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let two = (2u64, [1.0, 0.0, 0.0], [2.0, 1.0, 1.0]);
        for order in [[one, two], [two, one]] {
            let mut sap = SweepAndPrune::new();
            for (id, min, max) in order {
                sap.add_object(id, min, max);
            }
            assert_eq!(sap.query_overlapping_pairs(), vec![(1, 2)]);
        }
        // Separated by the smallest representable gap: no pair.
        let mut sap = SweepAndPrune::new();
        sap.add_object(1, [0.0; 3], [1.0; 3]);
        sap.add_object(2, [1.0 + f64::EPSILON, 0.0, 0.0], [2.0, 1.0, 1.0]);
        assert!(sap.query_overlapping_pairs().is_empty());
    }

    #[test]
    fn sap_negative_zero_end_point_touches_a_positive_zero_start_point() {
        // -0.0 == 0.0 in the interval test, so the boxes touch and must be paired
        // whatever the sign of the zero.
        for (first_max, second_min) in [(-0.0, 0.0), (0.0, -0.0)] {
            let mut sap = SweepAndPrune::new();
            sap.add_object(1, [-1.0, 0.0, 0.0], [first_max, 1.0, 1.0]);
            sap.add_object(2, [second_min, 0.0, 0.0], [1.0, 1.0, 1.0]);
            assert_eq!(
                sap.query_overlapping_pairs(),
                vec![(1, 2)],
                "{first_max:?} vs {second_min:?}"
            );
        }
    }

    #[test]
    fn sap_sort_axis_puts_minimum_endpoints_first_on_equal_values() {
        let mut sap = SweepAndPrune::new();
        sap.add_object(1, [0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        sap.add_object(2, [1.0, 0.0, 0.0], [2.0, 1.0, 1.0]);
        sap.sort_axis();
        let tie: Vec<(u64, bool)> = sap
            .sorted_x
            .iter()
            .filter(|e| e.value == 1.0)
            .map(|e| (e.object_id, e.is_min))
            .collect();
        assert_eq!(tie, vec![(2, true), (1, false)]);
    }

    #[test]
    fn sap_query_repairs_a_list_edited_through_the_public_field() {
        let mut sap = SweepAndPrune::new();
        sap.add_object_batch(&[
            (1, [0.0, 0.0, 0.0], [2.0, 2.0, 2.0]),
            (2, [1.0, 0.0, 0.0], [3.0, 2.0, 2.0]),
            (3, [10.0, 0.0, 0.0], [12.0, 2.0, 2.0]),
        ]);
        assert_eq!(sap.query_overlapping_pairs(), vec![(1, 2)]);
        // Reverse the (clean) list behind the structure's back.
        sap.sorted_x.reverse();
        assert_eq!(sap.query_overlapping_pairs(), vec![(1, 2)]);
        assert!(endpoints_ascending(&sap));
    }

    #[test]
    fn sap_remove_object_keeps_the_remaining_endpoints_sweepable() {
        let mut sap = SweepAndPrune::new();
        for i in 0u64..4 {
            let x = i as f64 * 0.5;
            sap.add_object(i, [x, 0.0, 0.0], [x + 1.0, 1.0, 1.0]);
        }
        sap.remove_object(1);
        let pairs = sap.query_overlapping_pairs();
        assert!(!pairs.iter().any(|&(a, b)| a == 1 || b == 1));
        // [0, 1] touches [1, 2] and [1, 2] overlaps [1.5, 2.5]; [0, 1] and
        // [1.5, 2.5] are apart.
        assert_eq!(pairs, vec![(0, 2), (2, 3)]);
    }
}
