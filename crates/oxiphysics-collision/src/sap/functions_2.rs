//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
#[allow(unused_imports)]
use super::types::*;
#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;
    use std::collections::HashSet;
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
