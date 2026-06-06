//! SAP tests (part 1)
//!
//! This file was split from the original functions_2.rs (2175 lines) to comply
//! with the 2000-line policy. See also functions_2b.rs and functions_2c.rs.

#[cfg(test)]
mod tests {
    use super::super::functions::*;
    use super::super::types::*;
    #[test]
    fn sap_two_overlapping_aabbs_gives_one_pair() {
        let mut sap = SweepAndPrune::new();
        sap.add_object(1, [0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
        sap.add_object(2, [1.0, 1.0, 1.0], [3.0, 3.0, 3.0]);
        let pairs = sap.query_overlapping_pairs();
        assert_eq!(pairs.len(), 1);
        assert_eq!(pairs[0], (1, 2));
    }
    #[test]
    fn sap_two_separated_aabbs_gives_no_pairs() {
        let mut sap = SweepAndPrune::new();
        sap.add_object(1, [0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        sap.add_object(2, [5.0, 5.0, 5.0], [6.0, 6.0, 6.0]);
        let pairs = sap.query_overlapping_pairs();
        assert!(pairs.is_empty());
    }
    #[test]
    fn sap_three_aabbs_only_two_overlap() {
        let mut sap = SweepAndPrune::new();
        sap.add_object(1, [0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
        sap.add_object(2, [1.0, 1.0, 1.0], [3.0, 3.0, 3.0]);
        sap.add_object(3, [10.0, 10.0, 10.0], [12.0, 12.0, 12.0]);
        let pairs = sap.query_overlapping_pairs();
        assert_eq!(pairs.len(), 1);
        assert_eq!(pairs[0], (1, 2));
    }
    #[test]
    fn overlaps_on_axis_touching() {
        assert!(SweepAndPrune::overlaps_on_axis(0.0, 1.0, 1.0, 2.0));
    }
    #[test]
    fn overlaps_on_axis_separated() {
        assert!(!SweepAndPrune::overlaps_on_axis(0.0, 1.0, 2.0, 3.0));
    }
    #[test]
    fn overlaps_on_axis_contained() {
        assert!(SweepAndPrune::overlaps_on_axis(0.0, 10.0, 3.0, 7.0));
    }
    #[test]
    fn grid_two_objects_same_cell_gives_one_pair() {
        let mut grid = GridBroadphase::new(10.0);
        grid.insert(1, [1.0, 1.0, 1.0], [2.0, 2.0, 2.0]);
        grid.insert(2, [3.0, 3.0, 3.0], [4.0, 4.0, 4.0]);
        let pairs = grid.query_potential_pairs();
        assert_eq!(pairs.len(), 1);
        assert_eq!(pairs[0], (1, 2));
    }
    #[test]
    fn grid_adjacent_cells_still_potential_pair() {
        let mut grid = GridBroadphase::new(5.0);
        grid.insert(1, [0.0, 0.0, 0.0], [4.0, 4.0, 4.0]);
        grid.insert(2, [5.0, 0.0, 0.0], [9.0, 4.0, 4.0]);
        let pairs = grid.query_potential_pairs();
        assert!(
            pairs.is_empty(),
            "expected no pairs for objects in different cells, got {:?}",
            pairs
        );
        grid.clear();
        grid.insert(10, [3.0, 0.0, 0.0], [6.0, 4.0, 4.0]);
        grid.insert(20, [4.0, 0.0, 0.0], [7.0, 4.0, 4.0]);
        let pairs2 = grid.query_potential_pairs();
        assert!(
            !pairs2.is_empty(),
            "expected at least one pair for objects sharing a cell"
        );
    }
    #[test]
    fn isap_insert_and_query_overlapping() {
        let mut sap = IncrementalSap::new();
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
        let pairs = sap.query_pairs();
        assert_eq!(pairs.len(), 1);
        assert_eq!(pairs[0], (1, 2));
    }
    #[test]
    fn isap_separated_no_pairs() {
        let mut sap = IncrementalSap::new();
        sap.insert(
            1,
            Aabb3 {
                min: [0.0, 0.0, 0.0],
                max: [1.0, 1.0, 1.0],
            },
        );
        sap.insert(
            2,
            Aabb3 {
                min: [5.0, 5.0, 5.0],
                max: [6.0, 6.0, 6.0],
            },
        );
        let pairs = sap.query_pairs();
        assert!(pairs.is_empty());
    }
    #[test]
    fn isap_remove_body() {
        let mut sap = IncrementalSap::new();
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
        sap.remove(2);
        let pairs = sap.query_pairs();
        assert!(pairs.is_empty());
        assert_eq!(sap.body_count(), 1);
    }
    #[test]
    fn isap_update_moves_apart() {
        let mut sap = IncrementalSap::new();
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
        assert_eq!(sap.query_pairs().len(), 1);
        sap.update(
            2,
            Aabb3 {
                min: [10.0, 10.0, 10.0],
                max: [12.0, 12.0, 12.0],
            },
        );
        assert!(sap.query_pairs().is_empty());
    }
    #[test]
    fn isap_multi_axis_filtering() {
        let mut sap = IncrementalSap::new();
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
                min: [1.0, 1.0, 5.0],
                max: [3.0, 3.0, 6.0],
            },
        );
        let pairs = sap.query_pairs();
        assert!(pairs.is_empty(), "should not overlap when Z is separated");
    }
    #[test]
    fn isap_three_bodies_two_pairs() {
        let mut sap = IncrementalSap::new();
        sap.insert(
            1,
            Aabb3 {
                min: [0.0, 0.0, 0.0],
                max: [3.0, 3.0, 3.0],
            },
        );
        sap.insert(
            2,
            Aabb3 {
                min: [2.0, 2.0, 2.0],
                max: [5.0, 5.0, 5.0],
            },
        );
        sap.insert(
            3,
            Aabb3 {
                min: [4.0, 4.0, 4.0],
                max: [7.0, 7.0, 7.0],
            },
        );
        let pairs = sap.query_pairs();
        assert_eq!(pairs.len(), 2);
        assert!(pairs.contains(&(1, 2)));
        assert!(pairs.contains(&(2, 3)));
    }
    #[test]
    fn isap_batch_update() {
        let mut sap = IncrementalSap::new();
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
                min: [10.0, 10.0, 10.0],
                max: [12.0, 12.0, 12.0],
            },
        );
        assert!(sap.query_pairs().is_empty());
        sap.batch_update(&[(
            2,
            Aabb3 {
                min: [1.0, 1.0, 1.0],
                max: [3.0, 3.0, 3.0],
            },
        )]);
        let pairs = sap.query_pairs();
        assert_eq!(pairs.len(), 1);
        assert_eq!(pairs[0], (1, 2));
    }
    #[test]
    fn isap_sort_and_sweep_single_axis() {
        let mut eps = vec![
            SapEndpointU32 {
                value: 0.0,
                body_id: 1,
                is_min: true,
            },
            SapEndpointU32 {
                value: 2.0,
                body_id: 1,
                is_min: false,
            },
            SapEndpointU32 {
                value: 1.0,
                body_id: 2,
                is_min: true,
            },
            SapEndpointU32 {
                value: 3.0,
                body_id: 2,
                is_min: false,
            },
        ];
        let pairs = IncrementalSap::sort_and_sweep_axis(&mut eps);
        assert!(pairs.contains(&(1, 2)));
        assert_eq!(pairs.len(), 1);
    }
    #[test]
    fn stat_sap_tracks_pair_count() {
        let mut sap = StatTrackingSap::new();
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
        let pairs = sap.query_pairs();
        assert_eq!(pairs.len(), 1);
        assert_eq!(sap.stats.pair_count, 1);
        assert_eq!(sap.stats.body_count, 2);
    }
    #[test]
    fn stat_sap_no_pairs_stats() {
        let mut sap = StatTrackingSap::new();
        sap.insert(
            1,
            Aabb3 {
                min: [0.0, 0.0, 0.0],
                max: [1.0, 1.0, 1.0],
            },
        );
        sap.insert(
            2,
            Aabb3 {
                min: [5.0, 5.0, 5.0],
                max: [6.0, 6.0, 6.0],
            },
        );
        let pairs = sap.query_pairs();
        assert!(pairs.is_empty());
        assert_eq!(sap.stats.pair_count, 0);
        assert!(sap.stats.sweep_count > 0);
    }
    #[test]
    fn stat_sap_remove_body() {
        let mut sap = StatTrackingSap::new();
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
        sap.remove(2);
        let pairs = sap.query_pairs();
        assert!(pairs.is_empty());
        assert_eq!(sap.stats.body_count, 1);
    }
    #[test]
    fn translate_body_moves_apart() {
        let mut sap = IncrementalSap::new();
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
        assert_eq!(sap.query_pairs().len(), 1);
        translate_body(&mut sap, 2, [20.0, 0.0, 0.0]);
        assert!(
            sap.query_pairs().is_empty(),
            "bodies should be separated after translation"
        );
    }
    #[test]
    fn expand_aabb_creates_overlap() {
        let mut sap = IncrementalSap::new();
        sap.insert(
            1,
            Aabb3 {
                min: [0.0, 0.0, 0.0],
                max: [1.0, 1.0, 1.0],
            },
        );
        sap.insert(
            2,
            Aabb3 {
                min: [3.0, 3.0, 3.0],
                max: [4.0, 4.0, 4.0],
            },
        );
        assert!(sap.query_pairs().is_empty());
        expand_aabb(&mut sap, 1, 2.5);
        expand_aabb(&mut sap, 2, 2.5);
        let pairs = sap.query_pairs();
        assert_eq!(pairs.len(), 1, "expanding both AABBs should create overlap");
    }
    #[test]
    fn propagate_aabb_update_changes_aabb() {
        let mut sap = IncrementalSap::new();
        sap.insert(
            1,
            Aabb3 {
                min: [0.0, 0.0, 0.0],
                max: [1.0, 1.0, 1.0],
            },
        );
        sap.insert(
            2,
            Aabb3 {
                min: [5.0, 5.0, 5.0],
                max: [6.0, 6.0, 6.0],
            },
        );
        assert!(sap.query_pairs().is_empty());
        propagate_aabb_update(
            &mut sap,
            2,
            Aabb3 {
                min: [0.5, 0.5, 0.5],
                max: [1.5, 1.5, 1.5],
            },
        );
        let pairs = sap.query_pairs();
        assert_eq!(pairs.len(), 1);
    }
    #[test]
    fn grid_three_overlapping_all_pairs() {
        let mut grid = GridBroadphase::new(10.0);
        grid.insert(1, [1.0, 1.0, 1.0], [2.0, 2.0, 2.0]);
        grid.insert(2, [2.0, 2.0, 2.0], [3.0, 3.0, 3.0]);
        grid.insert(3, [3.0, 3.0, 3.0], [4.0, 4.0, 4.0]);
        let pairs = grid.query_potential_pairs();
        assert!(
            pairs.len() >= 3,
            "three objects in same cell should give 3 pairs"
        );
    }
    #[test]
    fn sap_remove_and_re_add() {
        let mut sap = SweepAndPrune::new();
        sap.add_object(1, [0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
        sap.add_object(2, [1.0, 1.0, 1.0], [3.0, 3.0, 3.0]);
        assert_eq!(sap.query_overlapping_pairs().len(), 1);
        sap.remove_object(2);
        assert!(sap.query_overlapping_pairs().is_empty());
        sap.add_object(2, [1.0, 1.0, 1.0], [3.0, 3.0, 3.0]);
        assert_eq!(sap.query_overlapping_pairs().len(), 1);
    }
    #[test]
    fn sap_update_separates_objects() {
        let mut sap = SweepAndPrune::new();
        sap.add_object(1, [0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
        sap.add_object(2, [1.0, 1.0, 1.0], [3.0, 3.0, 3.0]);
        assert_eq!(sap.query_overlapping_pairs().len(), 1);
        sap.update_object(2, [10.0, 10.0, 10.0], [12.0, 12.0, 12.0]);
        assert!(sap.query_overlapping_pairs().is_empty());
    }
    #[test]
    fn sap_axis_insert_and_sweep() {
        let mut axis = SapAxis::new();
        axis.insert(1, 0.0, 2.0);
        axis.insert(2, 1.0, 3.0);
        let pairs = axis.overlapping_pairs();
        assert!(
            pairs.contains(&(1, 2)),
            "overlapping endpoints should form a pair"
        );
    }
    #[test]
    fn sap_axis_non_overlapping() {
        let mut axis = SapAxis::new();
        axis.insert(1, 0.0, 1.0);
        axis.insert(2, 5.0, 6.0);
        assert!(axis.overlapping_pairs().is_empty());
    }
    #[test]
    fn sap_axis_update_no_longer_overlaps() {
        let mut axis = SapAxis::new();
        axis.insert(1, 0.0, 2.0);
        axis.insert(2, 1.0, 3.0);
        assert!(!axis.overlapping_pairs().is_empty());
        axis.update(2, 10.0, 12.0);
        assert!(axis.overlapping_pairs().is_empty());
    }
    #[test]
    fn sap_axis_remove() {
        let mut axis = SapAxis::new();
        axis.insert(1, 0.0, 2.0);
        axis.insert(2, 1.0, 3.0);
        axis.remove(2);
        assert!(axis.overlapping_pairs().is_empty());
        assert_eq!(axis.endpoints.len(), 2);
    }
    #[test]
    fn isap_insert_aabb_two_overlapping() {
        let mut sap = IncrementalSap::new();
        sap.insert_aabb(1, [0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
        sap.insert_aabb(2, [1.0, 1.0, 1.0], [3.0, 3.0, 3.0]);
        assert_eq!(
            sap.active_pairs().len(),
            1,
            "two overlapping AABBs should give one pair"
        );
    }
    #[test]
    fn isap_update_aabb_move_apart() {
        let mut sap = IncrementalSap::new();
        sap.insert_aabb(1, [0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
        sap.insert_aabb(2, [1.0, 1.0, 1.0], [3.0, 3.0, 3.0]);
        assert_eq!(sap.active_pairs().len(), 1);
        sap.update_aabb(2, [20.0, 20.0, 20.0], [22.0, 22.0, 22.0]);
        assert_eq!(
            sap.active_pairs().len(),
            0,
            "after moving apart, no pairs expected"
        );
    }
    #[test]
    fn isap_remove_aabb_clears_pair() {
        let mut sap = IncrementalSap::new();
        sap.insert_aabb(1, [0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
        sap.insert_aabb(2, [1.0, 1.0, 1.0], [3.0, 3.0, 3.0]);
        assert_eq!(sap.active_pairs().len(), 1);
        sap.remove_aabb(2);
        assert_eq!(
            sap.active_pairs().len(),
            0,
            "removing one body should clear its pairs"
        );
    }
    #[test]
    fn isap_bipartite_pairs_cross_set_only() {
        let mut sap = IncrementalSap::new();
        sap.insert_aabb(1, [0.0, 0.0, 0.0], [3.0, 3.0, 3.0]);
        sap.insert_aabb(2, [1.0, 1.0, 1.0], [4.0, 4.0, 4.0]);
        sap.insert_aabb(3, [2.0, 2.0, 2.0], [5.0, 5.0, 5.0]);
        let pairs = sap.bipartite_pairs(&[1, 2], &[3]);
        for &(a, b) in &pairs {
            let in_a = a == 1 || a == 2;
            let in_b = b == 3;
            let in_a2 = a == 3;
            let in_b2 = b == 1 || b == 2;
            assert!(
                (in_a && in_b) || (in_a2 && in_b2),
                "all pairs must be cross-set, got ({}, {})",
                a,
                b
            );
        }
        assert!(!pairs.is_empty(), "should find at least one cross-set pair");
    }
    #[test]
    fn event_sap_begin_event_on_first_overlap() {
        let mut sap = EventDrivenSap::new();
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
        let (events, _pairs) = sap.step_pairs();
        assert!(
            events.contains(&OverlapEvent::Begin(1, 2)),
            "expected Begin(1,2), got {:?}",
            events
        );
    }
    #[test]
    fn event_sap_no_event_for_stable_overlap() {
        let mut sap = EventDrivenSap::new();
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
        let _ = sap.step_pairs();
        let (events2, _) = sap.step_pairs();
        assert!(
            events2.is_empty(),
            "stable overlap should emit no events, got {:?}",
            events2
        );
    }
    #[test]
    fn event_sap_end_event_on_separation() {
        let mut sap = EventDrivenSap::new();
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
        let _ = sap.step_pairs();
        sap.update(
            2,
            Aabb3 {
                min: [20.0, 20.0, 20.0],
                max: [22.0, 22.0, 22.0],
            },
        );
        let (events2, _) = sap.step_pairs();
        assert!(
            events2.contains(&OverlapEvent::End(1, 2)),
            "expected End(1,2), got {:?}",
            events2
        );
    }
    #[test]
    fn event_sap_remove_clears_prev_pair() {
        let mut sap = EventDrivenSap::new();
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
        let _ = sap.step_pairs();
        sap.remove(2);
        let (events, pairs) = sap.step_pairs();
        assert!(pairs.is_empty(), "after removal, no pairs expected");
        assert!(
            events.contains(&OverlapEvent::End(1, 2)),
            "expected End(1,2) after removal, got {:?}",
            events
        );
    }
    #[test]
    fn event_sap_three_bodies_begin_end_sequence() {
        let mut sap = EventDrivenSap::new();
        sap.insert(
            1,
            Aabb3 {
                min: [0.0, 0.0, 0.0],
                max: [3.0, 3.0, 3.0],
            },
        );
        sap.insert(
            2,
            Aabb3 {
                min: [1.0, 1.0, 1.0],
                max: [4.0, 4.0, 4.0],
            },
        );
        sap.insert(
            3,
            Aabb3 {
                min: [2.0, 2.0, 2.0],
                max: [5.0, 5.0, 5.0],
            },
        );
        let (events1, pairs1) = sap.step_pairs();
        assert!(
            pairs1.len() >= 2,
            "three overlapping bodies, at least 2 pairs"
        );
        for ev in &events1 {
            assert!(matches!(ev, OverlapEvent::Begin(_, _)));
        }
        sap.update(
            3,
            Aabb3 {
                min: [50.0, 50.0, 50.0],
                max: [52.0, 52.0, 52.0],
            },
        );
        let (events2, pairs2) = sap.step_pairs();
        assert_eq!(pairs2.len(), 1, "only pair (1,2) should remain");
        let end_count = events2
            .iter()
            .filter(|e| matches!(e, OverlapEvent::End(_, _)))
            .count();
        assert!(end_count >= 1, "at least one End event expected");
    }
}
