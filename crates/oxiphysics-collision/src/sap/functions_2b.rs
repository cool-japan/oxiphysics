//! SAP tests (part 2)
//!
//! Continuation of functions_2.rs — split to comply with the 2000-line policy.

#[cfg(test)]
mod tests {
    use super::super::functions::*;
    use super::super::types::*;
    #[test]
    fn bubble_sort_already_sorted_no_swaps() {
        let mut eps = vec![
            SapEndpointU32 {
                value: 1.0,
                body_id: 1,
                is_min: true,
            },
            SapEndpointU32 {
                value: 2.0,
                body_id: 2,
                is_min: true,
            },
            SapEndpointU32 {
                value: 3.0,
                body_id: 3,
                is_min: true,
            },
        ];
        let swaps = bubble_sort_endpoint(&mut eps, 1);
        assert_eq!(swaps, 0);
        assert_eq!(eps[1].value, 2.0);
    }
    #[test]
    fn bubble_sort_moves_left() {
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
            SapEndpointU32 {
                value: 2.0,
                body_id: 3,
                is_min: true,
            },
        ];
        bubble_sort_endpoint(&mut eps, 2);
        assert_eq!(eps[0].value, 1.0);
        assert_eq!(eps[1].value, 2.0);
        assert_eq!(eps[2].value, 3.0);
    }
    #[test]
    fn bubble_sort_moves_right() {
        let mut eps = vec![
            SapEndpointU32 {
                value: 2.0,
                body_id: 1,
                is_min: true,
            },
            SapEndpointU32 {
                value: 1.0,
                body_id: 2,
                is_min: true,
            },
            SapEndpointU32 {
                value: 3.0,
                body_id: 3,
                is_min: true,
            },
        ];
        bubble_sort_endpoint(&mut eps, 1);
        assert_eq!(eps[0].value, 1.0);
        assert_eq!(eps[1].value, 2.0);
        assert_eq!(eps[2].value, 3.0);
    }
    #[test]
    fn aabb_union_basic() {
        let a = Aabb3 {
            min: [0.0, 0.0, 0.0],
            max: [1.0, 1.0, 1.0],
        };
        let b = Aabb3 {
            min: [0.5, 0.5, 0.5],
            max: [2.0, 2.0, 2.0],
        };
        let u = aabb_union(&a, &b);
        assert_eq!(u.min, [0.0, 0.0, 0.0]);
        assert_eq!(u.max, [2.0, 2.0, 2.0]);
    }
    #[test]
    fn aabb_intersection_overlapping() {
        let a = Aabb3 {
            min: [0.0, 0.0, 0.0],
            max: [2.0, 2.0, 2.0],
        };
        let b = Aabb3 {
            min: [1.0, 1.0, 1.0],
            max: [3.0, 3.0, 3.0],
        };
        let i = aabb_intersection(&a, &b).expect("should intersect");
        assert_eq!(i.min, [1.0, 1.0, 1.0]);
        assert_eq!(i.max, [2.0, 2.0, 2.0]);
    }
    #[test]
    fn aabb_intersection_separated_returns_none() {
        let a = Aabb3 {
            min: [0.0, 0.0, 0.0],
            max: [1.0, 1.0, 1.0],
        };
        let b = Aabb3 {
            min: [2.0, 2.0, 2.0],
            max: [3.0, 3.0, 3.0],
        };
        assert!(aabb_intersection(&a, &b).is_none());
    }
    #[test]
    fn aabb_surface_area_unit_cube() {
        let a = Aabb3 {
            min: [0.0; 3],
            max: [1.0; 3],
        };
        assert!((aabb_surface_area(&a) - 6.0).abs() < 1e-10);
    }
    #[test]
    fn aabb_volume_unit_cube() {
        let a = Aabb3 {
            min: [0.0; 3],
            max: [1.0; 3],
        };
        assert!((aabb_volume(&a) - 1.0).abs() < 1e-10);
    }
    #[test]
    fn aabb_contains_point_inside() {
        let a = Aabb3 {
            min: [0.0; 3],
            max: [2.0; 3],
        };
        assert!(aabb_contains_point(&a, [1.0, 1.0, 1.0]));
    }
    #[test]
    fn aabb_contains_point_outside() {
        let a = Aabb3 {
            min: [0.0; 3],
            max: [1.0; 3],
        };
        assert!(!aabb_contains_point(&a, [2.0, 0.5, 0.5]));
    }
    #[test]
    fn aabb_contains_aabb_true() {
        let outer = Aabb3 {
            min: [0.0; 3],
            max: [4.0; 3],
        };
        let inner = Aabb3 {
            min: [1.0; 3],
            max: [3.0; 3],
        };
        assert!(aabb_contains_aabb(&outer, &inner));
    }
    #[test]
    fn aabb_contains_aabb_false() {
        let a = Aabb3 {
            min: [0.0; 3],
            max: [2.0; 3],
        };
        let b = Aabb3 {
            min: [1.0; 3],
            max: [4.0; 3],
        };
        assert!(!aabb_contains_aabb(&a, &b));
    }
    #[test]
    fn aabb_pad_increases_size() {
        let a = Aabb3 {
            min: [1.0; 3],
            max: [3.0; 3],
        };
        let padded = aabb_pad(&a, 0.5);
        assert_eq!(padded.min, [0.5; 3]);
        assert_eq!(padded.max, [3.5; 3]);
    }
    #[test]
    fn aabb_center_unit_cube() {
        let a = Aabb3 {
            min: [0.0; 3],
            max: [2.0; 3],
        };
        let c = aabb_center(&a);
        assert_eq!(c, [1.0; 3]);
    }
    #[test]
    fn aabb_half_extents_unit_cube() {
        let a = Aabb3 {
            min: [0.0; 3],
            max: [2.0; 3],
        };
        let he = aabb_half_extents(&a);
        assert_eq!(he, [1.0; 3]);
    }
    #[test]
    fn multi_phase_sap_step_returns_pairs() {
        let mut sap = MultiPhaseSap::new();
        sap.insert(
            1,
            Aabb3 {
                min: [0.0, 0.0, 0.0],
                max: [2.0, 2.0, 2.0],
            },
        );
        sap.insert(
            2,
            Aabb3 {
                min: [1.0, 1.0, 1.0],
                max: [3.0, 3.0, 3.0],
            },
        );
        let result = sap.step();
        assert_eq!(result.pairs.len(), 1);
        assert!(result.new_pairs.contains(&(1, 2)));
        assert!(result.lost_pairs.is_empty());
    }
    #[test]
    fn multi_phase_sap_step_records_stats() {
        let mut sap = MultiPhaseSap::new();
        sap.insert(
            1,
            Aabb3 {
                min: [0.0, 0.0, 0.0],
                max: [2.0, 2.0, 2.0],
            },
        );
        sap.insert(
            2,
            Aabb3 {
                min: [1.0, 1.0, 1.0],
                max: [3.0, 3.0, 3.0],
            },
        );
        let result = sap.step();
        assert_eq!(result.stats.pair_count, 1);
        assert_eq!(result.stats.body_count, 2);
        assert!(result.stats.sweep_count > 0);
    }
    #[test]
    fn multi_phase_sap_remove_fires_lost_pair() {
        let mut sap = MultiPhaseSap::new();
        sap.insert(
            1,
            Aabb3 {
                min: [0.0, 0.0, 0.0],
                max: [2.0, 2.0, 2.0],
            },
        );
        sap.insert(
            2,
            Aabb3 {
                min: [1.0, 1.0, 1.0],
                max: [3.0, 3.0, 3.0],
            },
        );
        let _ = sap.step();
        sap.remove(2);
        let result = sap.step();
        assert!(result.pairs.is_empty());
        assert!(result.lost_pairs.contains(&(1, 2)));
    }
    #[test]
    fn sap_endpoint_stats_counts_correctly() {
        let mut sap = IncrementalSap::new();
        sap.insert_aabb(1, [0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
        sap.insert_aabb(2, [1.0, 1.0, 1.0], [3.0, 3.0, 3.0]);
        let stats = sap_endpoint_stats(&sap);
        assert_eq!(stats.sweep_count, 4);
        assert_eq!(stats.body_count, 2);
        assert_eq!(stats.pair_count, 1);
    }
    #[test]
    fn find_endpoint_helpers() {
        let mut sap = IncrementalSap::new();
        sap.insert_aabb(5, [1.0, 0.0, 0.0], [3.0, 1.0, 1.0]);
        let min_idx = find_min_endpoint(&sap.endpoints_x, 5);
        let max_idx = find_max_endpoint(&sap.endpoints_x, 5);
        assert!(min_idx.is_some(), "min endpoint for body 5 should be found");
        assert!(max_idx.is_some(), "max endpoint for body 5 should be found");
        assert_ne!(
            min_idx.unwrap(),
            max_idx.unwrap(),
            "min and max should be different indices"
        );
    }
    #[test]
    fn sap_empty_gives_no_pairs() {
        let mut sap = SweepAndPrune::new();
        assert!(sap.query_overlapping_pairs().is_empty());
    }
    #[test]
    fn sap_single_object_gives_no_pairs() {
        let mut sap = SweepAndPrune::new();
        sap.add_object(1, [0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert!(sap.query_overlapping_pairs().is_empty());
    }
    #[test]
    fn sap_object_count_after_insert_and_remove() {
        let mut sap = SweepAndPrune::new();
        assert_eq!(sap.object_count(), 0);
        sap.add_object(1, [0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert_eq!(sap.object_count(), 1);
        sap.add_object(2, [2.0, 2.0, 2.0], [3.0, 3.0, 3.0]);
        assert_eq!(sap.object_count(), 2);
        sap.remove_object(1);
        assert_eq!(sap.object_count(), 1);
    }
    #[test]
    fn sap_overlaps_on_axis_partial_overlap() {
        assert!(SweepAndPrune::overlaps_on_axis(0.0, 3.0, 2.0, 5.0));
        assert!(SweepAndPrune::overlaps_on_axis(2.0, 5.0, 0.0, 3.0));
    }
    #[test]
    fn sap_overlaps_on_axis_identical_intervals() {
        assert!(SweepAndPrune::overlaps_on_axis(1.0, 3.0, 1.0, 3.0));
    }
    #[test]
    fn sap_overlaps_on_axis_point_interval() {
        assert!(SweepAndPrune::overlaps_on_axis(1.0, 1.0, 1.0, 2.0));
        assert!(!SweepAndPrune::overlaps_on_axis(1.0, 1.0, 2.0, 3.0));
    }
    #[test]
    fn sap_many_objects_all_separated() {
        let mut sap = SweepAndPrune::new();
        for i in 0..10u64 {
            let x = i as f64 * 10.0;
            sap.add_object(i, [x, 0.0, 0.0], [x + 1.0, 1.0, 1.0]);
        }
        assert!(
            sap.query_overlapping_pairs().is_empty(),
            "all separated objects should give no pairs"
        );
    }
    #[test]
    fn sap_many_objects_all_overlapping() {
        let mut sap = SweepAndPrune::new();
        for i in 0..5u64 {
            sap.add_object(i, [0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
        }
        let pairs = sap.query_overlapping_pairs();
        assert_eq!(pairs.len(), 10, "5 coincident objects should give 10 pairs");
    }
    #[test]
    fn isap_default_is_empty() {
        let sap = IncrementalSap::default();
        assert_eq!(sap.body_count(), 0);
        assert!(sap.active_pairs.is_empty());
    }
    #[test]
    fn isap_insert_updates_body_count() {
        let mut sap = IncrementalSap::new();
        sap.insert(
            1,
            Aabb3 {
                min: [0.0; 3],
                max: [1.0; 3],
            },
        );
        assert_eq!(sap.body_count(), 1);
        sap.insert(
            2,
            Aabb3 {
                min: [0.0; 3],
                max: [1.0; 3],
            },
        );
        assert_eq!(sap.body_count(), 2);
    }
    #[test]
    fn isap_current_pairs_empty_initially() {
        let sap = IncrementalSap::new();
        assert!(sap.current_pairs().is_empty());
    }
    #[test]
    fn isap_current_pairs_updated_after_query() {
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
                min: [1.0; 3],
                max: [3.0; 3],
            },
        );
        sap.query_pairs();
        assert_eq!(sap.current_pairs().len(), 1);
    }
    #[test]
    fn isap_bipartite_empty_sets_returns_empty() {
        let sap = IncrementalSap::new();
        let pairs = sap.bipartite_pairs(&[], &[]);
        assert!(pairs.is_empty());
    }
    #[test]
    fn isap_bipartite_no_overlap_returns_empty() {
        let mut sap = IncrementalSap::new();
        sap.insert_aabb(1, [0.0; 3], [1.0; 3]);
        sap.insert_aabb(2, [10.0; 3], [11.0; 3]);
        let pairs = sap.bipartite_pairs(&[1], &[2]);
        assert!(
            pairs.is_empty(),
            "separated sets should give no bipartite pairs"
        );
    }
    #[test]
    fn aabb3_union_idempotent() {
        let a = Aabb3 {
            min: [1.0, 2.0, 3.0],
            max: [4.0, 5.0, 6.0],
        };
        let u = aabb_union(&a, &a);
        for i in 0..3 {
            assert!((u.min[i] - a.min[i]).abs() < 1e-12);
            assert!((u.max[i] - a.max[i]).abs() < 1e-12);
        }
    }
    #[test]
    fn aabb3_intersection_touching_single_plane() {
        let a = Aabb3 {
            min: [0.0, 0.0, 0.0],
            max: [1.0, 1.0, 1.0],
        };
        let b = Aabb3 {
            min: [1.0, 0.0, 0.0],
            max: [2.0, 1.0, 1.0],
        };
        let i = aabb_intersection(&a, &b);
        assert!(
            i.is_some(),
            "touching at plane should still be an intersection"
        );
        let sect = i.unwrap();
        assert!((sect.min[0] - 1.0).abs() < 1e-12);
        assert!((sect.max[0] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn aabb3_volume_zero_for_degenerate() {
        let a = Aabb3 {
            min: [0.0; 3],
            max: [0.0; 3],
        };
        assert!((aabb_volume(&a)).abs() < 1e-12);
    }
    #[test]
    fn aabb3_surface_area_zero_for_degenerate() {
        let a = Aabb3 {
            min: [0.0; 3],
            max: [0.0; 3],
        };
        assert!((aabb_surface_area(&a)).abs() < 1e-12);
    }
    #[test]
    fn aabb3_pad_zero_margin_unchanged() {
        let a = Aabb3 {
            min: [1.0, 2.0, 3.0],
            max: [4.0, 5.0, 6.0],
        };
        let padded = aabb_pad(&a, 0.0);
        for i in 0..3 {
            assert!((padded.min[i] - a.min[i]).abs() < 1e-12);
            assert!((padded.max[i] - a.max[i]).abs() < 1e-12);
        }
    }
    #[test]
    fn aabb3_half_extents_asymmetric() {
        let a = Aabb3 {
            min: [1.0, 2.0, 3.0],
            max: [3.0, 6.0, 11.0],
        };
        let he = aabb_half_extents(&a);
        assert!((he[0] - 1.0).abs() < 1e-12);
        assert!((he[1] - 2.0).abs() < 1e-12);
        assert!((he[2] - 4.0).abs() < 1e-12);
    }
    #[test]
    fn aabb3_center_asymmetric() {
        let a = Aabb3 {
            min: [0.0, 2.0, 4.0],
            max: [2.0, 4.0, 8.0],
        };
        let c = aabb_center(&a);
        assert!((c[0] - 1.0).abs() < 1e-12);
        assert!((c[1] - 3.0).abs() < 1e-12);
        assert!((c[2] - 6.0).abs() < 1e-12);
    }
    #[test]
    fn aabb3_contains_point_on_boundary() {
        let a = Aabb3 {
            min: [0.0; 3],
            max: [1.0; 3],
        };
        assert!(aabb_contains_point(&a, [0.0, 0.5, 0.5]));
        assert!(aabb_contains_point(&a, [1.0, 0.5, 0.5]));
        assert!(aabb_contains_point(&a, [0.5, 0.0, 0.5]));
        assert!(aabb_contains_point(&a, [0.5, 1.0, 0.5]));
    }
    #[test]
    fn aabb3_contains_aabb_equal() {
        let a = Aabb3 {
            min: [0.0; 3],
            max: [1.0; 3],
        };
        assert!(aabb_contains_aabb(&a, &a), "AABB should contain itself");
    }
    #[test]
    fn aabb3_contains_aabb_partial_overlap() {
        let a = Aabb3 {
            min: [0.0; 3],
            max: [2.0; 3],
        };
        let b = Aabb3 {
            min: [1.0; 3],
            max: [3.0; 3],
        };
        assert!(
            !aabb_contains_aabb(&a, &b),
            "partial overlap should not count as containment"
        );
    }
    #[test]
    fn sap_endpoint_stats_no_pairs() {
        let mut sap = IncrementalSap::new();
        sap.insert_aabb(1, [0.0; 3], [1.0; 3]);
        sap.insert_aabb(2, [10.0; 3], [11.0; 3]);
        let stats = sap_endpoint_stats(&sap);
        assert_eq!(stats.body_count, 2);
        assert_eq!(stats.sweep_count, 4);
        assert_eq!(stats.pair_count, 0);
    }
    #[test]
    fn grid_clear_empties_all_cells() {
        let mut grid = GridBroadphase::new(1.0);
        grid.insert(1, [0.0; 3], [0.5; 3]);
        grid.insert(2, [0.1; 3], [0.4; 3]);
        assert!(!grid.query_potential_pairs().is_empty());
        grid.clear();
        assert!(grid.query_potential_pairs().is_empty());
    }
    #[test]
    fn grid_cell_coord_positive() {
        let grid = GridBroadphase::new(2.0);
        assert_eq!(grid.cell_coord(0.0), 0);
        assert_eq!(grid.cell_coord(1.9), 0);
        assert_eq!(grid.cell_coord(2.0), 1);
        assert_eq!(grid.cell_coord(4.5), 2);
    }
    #[test]
    fn grid_cell_coord_negative() {
        let grid = GridBroadphase::new(2.0);
        assert_eq!(grid.cell_coord(-0.1), -1);
        assert_eq!(grid.cell_coord(-2.0), -1);
        assert_eq!(grid.cell_coord(-2.1), -2);
    }
    #[test]
    fn stat_tracking_sap_default_is_zero() {
        let sap = StatTrackingSap::default();
        assert_eq!(sap.stats.pair_count, 0);
        assert_eq!(sap.stats.body_count, 0);
        assert_eq!(sap.stats.sweep_count, 0);
    }
    #[test]
    fn event_sap_body_count_tracks_insertions_and_removals() {
        let mut sap = EventDrivenSap::new();
        assert_eq!(sap.body_count(), 0);
        sap.insert(
            1,
            Aabb3 {
                min: [0.0; 3],
                max: [1.0; 3],
            },
        );
        assert_eq!(sap.body_count(), 1);
        sap.insert(
            2,
            Aabb3 {
                min: [0.0; 3],
                max: [1.0; 3],
            },
        );
        assert_eq!(sap.body_count(), 2);
        sap.remove(1);
        assert_eq!(sap.body_count(), 1);
    }
    #[test]
    fn event_sap_default_body_count_is_zero() {
        let sap = EventDrivenSap::default();
        assert_eq!(sap.body_count(), 0);
    }
    #[test]
    fn event_sap_no_events_empty() {
        let mut sap = EventDrivenSap::new();
        let (events, pairs) = sap.step_pairs();
        assert!(events.is_empty());
        assert!(pairs.is_empty());
    }
    #[test]
    fn multi_phase_sap_update_aabb() {
        let mut sap = MultiPhaseSap::new();
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
                min: [1.0; 3],
                max: [3.0; 3],
            },
        );
        let r1 = sap.step();
        assert_eq!(r1.pairs.len(), 1);
        sap.update(
            2,
            Aabb3 {
                min: [50.0; 3],
                max: [52.0; 3],
            },
        );
        let r2 = sap.step();
        assert!(r2.pairs.is_empty(), "moved body should not overlap");
    }
    #[test]
    fn multi_phase_sap_default_creates_empty() {
        let mut sap = MultiPhaseSap::default();
        let result = sap.step();
        assert!(result.pairs.is_empty());
        assert!(result.new_pairs.is_empty());
        assert!(result.lost_pairs.is_empty());
    }
    #[test]
    fn isap_five_bodies_chain_overlap() {
        let mut sap = IncrementalSap::new();
        for i in 0..5u32 {
            let x = i as f64 * 1.5;
            sap.insert(
                i,
                Aabb3 {
                    min: [x, 0.0, 0.0],
                    max: [x + 2.0, 1.0, 1.0],
                },
            );
        }
        let pairs = sap.query_pairs();
        assert_eq!(
            pairs.len(),
            4,
            "chain of 5 should give 4 pairs, got {:?}",
            pairs
        );
    }
    #[test]
    fn sap_axis_sorted_after_insert() {
        let mut axis = SapAxis::new();
        axis.insert(3, 5.0, 8.0);
        axis.insert(1, 1.0, 4.0);
        axis.insert(2, 2.0, 6.0);
        for i in 1..axis.endpoints.len() {
            assert!(
                axis.endpoints[i - 1].value <= axis.endpoints[i].value,
                "endpoints must be sorted after insert"
            );
        }
    }
    #[test]
    fn find_min_endpoint_not_found() {
        let eps = vec![SapEndpointU32 {
            value: 1.0,
            body_id: 1,
            is_min: true,
        }];
        assert!(find_min_endpoint(&eps, 99).is_none());
    }
    #[test]
    fn find_max_endpoint_not_found() {
        let eps = vec![SapEndpointU32 {
            value: 1.0,
            body_id: 1,
            is_min: false,
        }];
        assert!(find_max_endpoint(&eps, 99).is_none());
    }
    #[test]
    fn translate_body_does_nothing_for_unknown_id() {
        let mut sap = IncrementalSap::new();
        sap.insert(
            1,
            Aabb3 {
                min: [0.0; 3],
                max: [1.0; 3],
            },
        );
        translate_body(&mut sap, 99, [5.0, 0.0, 0.0]);
        let aabb = sap.aabbs.get(&1).unwrap();
        assert_eq!(aabb.min, [0.0; 3]);
    }
    #[test]
    fn expand_aabb_does_nothing_for_unknown_id() {
        let mut sap = IncrementalSap::new();
        sap.insert(
            1,
            Aabb3 {
                min: [0.0; 3],
                max: [1.0; 3],
            },
        );
        expand_aabb(&mut sap, 99, 5.0);
        let aabb = sap.aabbs.get(&1).unwrap();
        assert_eq!(aabb.max, [1.0; 3]);
    }
    #[test]
    fn sap_add_object_batch_inserts_all() {
        let mut sap = SweepAndPrune::new();
        let objects: Vec<(u64, [f64; 3], [f64; 3])> = (0u64..4)
            .map(|i| {
                (
                    i,
                    [i as f64 * 5.0, 0.0, 0.0],
                    [i as f64 * 5.0 + 1.0, 1.0, 1.0],
                )
            })
            .collect();
        sap.add_object_batch(&objects);
        assert_eq!(sap.object_count(), 4, "all 4 objects should be inserted");
        assert_eq!(sap.sorted_x.len(), 8);
    }
    #[test]
    fn sap_add_object_batch_finds_correct_pairs() {
        let mut sap = SweepAndPrune::new();
        sap.add_object_batch(&[
            (0, [0.0, 0.0, 0.0], [2.0, 2.0, 2.0]),
            (1, [1.0, 0.0, 0.0], [3.0, 2.0, 2.0]),
            (2, [100.0, 0.0, 0.0], [102.0, 2.0, 2.0]),
        ]);
        let pairs = sap.query_overlapping_pairs();
        assert_eq!(pairs.len(), 1, "only (0,1) should overlap, got {:?}", pairs);
        assert!(pairs.contains(&(0, 1)) || pairs.contains(&(1, 0)));
    }
}
