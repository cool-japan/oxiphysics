//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[cfg(test)]
mod tests {
    use super::super::*;
    use crate::BruteForceBroadPhase;
    use crate::BvhBroadphase;
    use crate::CollisionPair;
    use crate::SweepAndPrune;
    use crate::broadphase::BroadphaseProfiler;
    use crate::broadphase::BroadphaseSceneGraph;
    use crate::broadphase::BroadphaseWarmstart;
    use crate::broadphase::DynamicAabbTree;
    use crate::broadphase::Frustum;
    use crate::broadphase::GpuBroadphaseHints;
    use crate::broadphase::ObjectType;
    use oxiphysics_core::Aabb;
    use oxiphysics_core::Real;
    use oxiphysics_core::Vec3;
    fn make_aabb(cx: Real, cy: Real, cz: Real, half: Real) -> Aabb {
        Aabb::new(
            Vec3::new(cx - half, cy - half, cz - half),
            Vec3::new(cx + half, cy + half, cz + half),
        )
    }
    #[test]
    fn test_brute_force() {
        let aabbs = vec![
            make_aabb(0.0, 0.0, 0.0, 1.0),
            make_aabb(1.5, 0.0, 0.0, 1.0),
            make_aabb(10.0, 0.0, 0.0, 1.0),
        ];
        let pairs = BruteForceBroadPhase.find_pairs(&aabbs);
        assert_eq!(pairs.len(), 1);
        assert_eq!(pairs[0].a, 0);
        assert_eq!(pairs[0].b, 1);
    }
    #[test]
    fn test_sweep_and_prune() {
        let aabbs = vec![
            make_aabb(0.0, 0.0, 0.0, 1.0),
            make_aabb(1.5, 0.0, 0.0, 1.0),
            make_aabb(10.0, 0.0, 0.0, 1.0),
        ];
        let sap = SweepAndPrune::x_axis();
        let pairs = sap.find_pairs(&aabbs);
        assert_eq!(pairs.len(), 1);
    }
    #[test]
    fn test_dynamic_bvh() {
        let aabbs = vec![
            make_aabb(0.0, 0.0, 0.0, 1.0),
            make_aabb(1.5, 0.0, 0.0, 1.0),
            make_aabb(10.0, 0.0, 0.0, 1.0),
        ];
        let mut bvh = BvhBroadphase::new();
        bvh.build(&aabbs);
        let pairs = bvh.find_pairs(&aabbs);
        assert_eq!(pairs.len(), 1);
    }
    #[test]
    fn test_bvh_query() {
        let aabbs = vec![
            make_aabb(0.0, 0.0, 0.0, 1.0),
            make_aabb(5.0, 0.0, 0.0, 1.0),
            make_aabb(10.0, 0.0, 0.0, 1.0),
        ];
        let mut bvh = BvhBroadphase::new();
        bvh.build(&aabbs);
        let query_aabb = make_aabb(0.0, 0.0, 0.0, 0.5);
        let hits = bvh.query(&query_aabb);
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0], 0);
    }
    #[test]
    fn test_bvh_ray_query() {
        let aabbs = vec![
            make_aabb(0.0, 0.0, 0.0, 1.0),
            make_aabb(5.0, 0.0, 0.0, 1.0),
            make_aabb(10.0, 0.0, 0.0, 1.0),
        ];
        let mut bvh = BvhBroadphase::new();
        bvh.build(&aabbs);
        let hits = bvh.ray_query(&Vec3::new(-5.0, 0.0, 0.0), &Vec3::new(1.0, 0.0, 0.0), 100.0);
        assert_eq!(hits.len(), 3);
    }
    #[test]
    fn test_sap_empty_input() {
        let sap = SweepAndPrune::x_axis();
        let pairs = sap.find_pairs(&[]);
        assert!(pairs.is_empty());
    }
    #[test]
    fn test_sap_one_aabb() {
        let sap = SweepAndPrune::x_axis();
        let aabbs = vec![make_aabb(0.0, 0.0, 0.0, 1.0)];
        let pairs = sap.find_pairs(&aabbs);
        assert!(pairs.is_empty());
    }
    #[test]
    fn test_sap_overlapping_pair() {
        let sap = SweepAndPrune::x_axis();
        let aabbs = vec![make_aabb(0.0, 0.0, 0.0, 1.0), make_aabb(1.0, 0.0, 0.0, 1.0)];
        let pairs = sap.find_pairs(&aabbs);
        assert_eq!(pairs.len(), 1);
        assert_eq!(pairs[0].a, 0);
        assert_eq!(pairs[0].b, 1);
    }
    #[test]
    fn test_sap_no_overlap() {
        let sap = SweepAndPrune::x_axis();
        let aabbs = vec![make_aabb(0.0, 0.0, 0.0, 1.0), make_aabb(5.0, 0.0, 0.0, 1.0)];
        let pairs = sap.find_pairs(&aabbs);
        assert!(pairs.is_empty());
    }
    /// Helper: sort pairs canonically so SAP and brute-force can be compared.
    fn sorted_pairs(mut pairs: Vec<CollisionPair>) -> Vec<(usize, usize)> {
        let mut v: Vec<(usize, usize)> = pairs
            .drain(..)
            .map(|p| {
                let (a, b) = if p.a < p.b { (p.a, p.b) } else { (p.b, p.a) };
                (a, b)
            })
            .collect();
        v.sort_unstable();
        v.dedup();
        v
    }
    #[test]
    fn test_sap_agrees_with_brute_force() {
        let mut aabbs = Vec::new();
        for i in 0..20_i32 {
            let x = (i % 5) as Real * 1.5;
            let y = (i / 5) as Real * 1.5;
            aabbs.push(make_aabb(x, y, 0.0, 1.0));
        }
        let sap = SweepAndPrune::x_axis();
        let sap_pairs = sorted_pairs(sap.find_pairs(&aabbs));
        let bf_pairs = sorted_pairs(BruteForceBroadPhase.find_pairs(&aabbs));
        assert_eq!(
            sap_pairs, bf_pairs,
            "SAP and brute-force disagree:\n  SAP: {:?}\n  BF:  {:?}",
            sap_pairs, bf_pairs
        );
    }
    #[test]
    fn test_dbvt_insert_query() {
        let aabbs = vec![
            make_aabb(0.0, 0.0, 0.0, 1.0),
            make_aabb(1.5, 0.0, 0.0, 1.0),
            make_aabb(3.0, 0.0, 0.0, 1.0),
            make_aabb(10.0, 0.0, 0.0, 1.0),
            make_aabb(0.5, 0.0, 0.0, 0.3),
        ];
        let mut bvh = BvhBroadphase::new();
        bvh.build(&aabbs);
        let query = Aabb::new(Vec3::new(-0.4, -0.4, -0.4), Vec3::new(0.4, 0.4, 0.4));
        let mut hits = bvh.query(&query);
        hits.sort_unstable();
        assert!(hits.contains(&0), "expected hit index 0, got {:?}", hits);
        assert!(hits.contains(&4), "expected hit index 4, got {:?}", hits);
        assert!(!hits.contains(&1), "unexpected hit index 1, got {:?}", hits);
        assert!(!hits.contains(&3), "unexpected hit index 3, got {:?}", hits);
    }
    #[test]
    fn test_dbvt_remove_reinsert() {
        let mut aabbs = vec![
            make_aabb(0.0, 0.0, 0.0, 1.0),
            make_aabb(1.5, 0.0, 0.0, 1.0),
            make_aabb(10.0, 0.0, 0.0, 1.0),
        ];
        let mut bvh = BvhBroadphase::new();
        bvh.build(&aabbs);
        let pairs_full = sorted_pairs(bvh.find_pairs(&aabbs));
        let aabbs_reduced = aabbs[..2].to_vec();
        bvh.build(&aabbs_reduced);
        let pairs_reduced = sorted_pairs(bvh.find_pairs(&aabbs_reduced));
        assert_eq!(pairs_reduced, vec![(0, 1)]);
        aabbs[2] = make_aabb(10.0, 0.0, 0.0, 1.0);
        bvh.build(&aabbs);
        let pairs_after = sorted_pairs(bvh.find_pairs(&aabbs));
        assert_eq!(
            pairs_after, pairs_full,
            "pairs after reinsert should equal original"
        );
    }
    #[test]
    fn test_dynamic_aabb_tree_insert_query() {
        let mut tree = DynamicAabbTree::new(0.1);
        tree.insert(make_aabb(0.0, 0.0, 0.0, 1.0));
        tree.insert(make_aabb(1.5, 0.0, 0.0, 1.0));
        tree.insert(make_aabb(10.0, 0.0, 0.0, 1.0));
        assert_eq!(tree.len(), 3);
        assert!(!tree.is_empty());
        let pairs = tree.find_pairs();
        assert_eq!(pairs.len(), 1, "Only first two should overlap");
    }
    #[test]
    fn test_dynamic_aabb_tree_update() {
        let mut tree = DynamicAabbTree::new(0.1);
        tree.insert(make_aabb(0.0, 0.0, 0.0, 1.0));
        tree.insert(make_aabb(10.0, 0.0, 0.0, 1.0));
        let pairs = tree.find_pairs();
        assert!(pairs.is_empty(), "Should have no pairs initially");
        tree.update(1, make_aabb(1.5, 0.0, 0.0, 1.0));
        let pairs = tree.find_pairs();
        assert_eq!(pairs.len(), 1, "Should overlap after move");
    }
    #[test]
    fn test_dynamic_aabb_tree_no_rebuild_within_margin() {
        let mut tree = DynamicAabbTree::new(1.0);
        tree.insert(make_aabb(0.0, 0.0, 0.0, 1.0));
        tree.rebuild_if_dirty();
        let needs_rebuild = tree.update(0, make_aabb(0.1, 0.0, 0.0, 1.0));
        assert!(
            !needs_rebuild,
            "Small move within margin should not trigger rebuild"
        );
    }
    #[test]
    fn test_batch_query() {
        let aabbs = vec![
            make_aabb(0.0, 0.0, 0.0, 1.0),
            make_aabb(5.0, 0.0, 0.0, 1.0),
            make_aabb(10.0, 0.0, 0.0, 1.0),
        ];
        let mut bvh = BvhBroadphase::new();
        bvh.build(&aabbs);
        let queries = vec![make_aabb(0.0, 0.0, 0.0, 0.5), make_aabb(5.0, 0.0, 0.0, 0.5)];
        let results = batch_query(&bvh, &queries);
        assert_eq!(results.len(), 2);
        assert!(results[0].contains(&0), "First query should hit object 0");
        assert!(results[1].contains(&1), "Second query should hit object 1");
    }
    #[test]
    fn test_batch_ray_query() {
        let aabbs = vec![make_aabb(0.0, 0.0, 0.0, 1.0), make_aabb(5.0, 0.0, 0.0, 1.0)];
        let mut bvh = BvhBroadphase::new();
        bvh.build(&aabbs);
        let rays = vec![
            (Vec3::new(-5.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0)),
            (Vec3::new(0.0, 5.0, 0.0), Vec3::new(0.0, -1.0, 0.0)),
        ];
        let results = batch_ray_query(&bvh, &rays, 100.0);
        assert_eq!(results.len(), 2);
        assert_eq!(results[0].len(), 2, "Ray along X should hit both objects");
    }
    #[test]
    fn test_frustum_contains_aabb() {
        let frustum = Frustum::from_perspective(std::f64::consts::FRAC_PI_4, 1.0, 0.1, 100.0);
        let visible = make_aabb(0.0, 0.0, -5.0, 1.0);
        assert!(
            frustum.contains_aabb(&visible),
            "Box in front of camera should be visible"
        );
        let behind = make_aabb(0.0, 0.0, 50.0, 1.0);
        assert!(
            !frustum.contains_aabb(&behind),
            "Box behind camera should not be visible"
        );
    }
    #[test]
    fn test_frustum_cull() {
        let frustum = Frustum::from_perspective(std::f64::consts::FRAC_PI_4, 1.0, 0.1, 100.0);
        let aabbs = vec![
            make_aabb(0.0, 0.0, -5.0, 1.0),
            make_aabb(0.0, 0.0, -50.0, 1.0),
            make_aabb(0.0, 0.0, 50.0, 1.0),
        ];
        let mut bvh = BvhBroadphase::new();
        bvh.build(&aabbs);
        let visible = frustum_cull(&bvh, &aabbs, &frustum);
        assert!(!visible.contains(&2), "Box behind camera should be culled");
    }
    #[test]
    fn test_brute_force_stats() {
        let aabbs = vec![
            make_aabb(0.0, 0.0, 0.0, 1.0),
            make_aabb(1.5, 0.0, 0.0, 1.0),
            make_aabb(10.0, 0.0, 0.0, 1.0),
        ];
        let (pairs, stats) = brute_force_with_stats(&aabbs);
        assert_eq!(stats.num_objects, 3);
        assert_eq!(stats.num_pairs, pairs.len());
        assert_eq!(stats.num_tests, 3);
    }
    #[test]
    fn test_sap_stats() {
        let aabbs = vec![make_aabb(0.0, 0.0, 0.0, 1.0), make_aabb(1.5, 0.0, 0.0, 1.0)];
        let (pairs, stats) = sap_with_stats(&aabbs, 0);
        assert_eq!(stats.num_objects, 2);
        assert_eq!(stats.num_pairs, pairs.len());
    }
    #[test]
    fn test_sap_y_axis() {
        let sap = SweepAndPrune::new(1);
        let aabbs = vec![
            make_aabb(0.0, 0.0, 0.0, 1.0),
            make_aabb(0.0, 1.5, 0.0, 1.0),
            make_aabb(0.0, 10.0, 0.0, 1.0),
        ];
        let pairs = sap.find_pairs(&aabbs);
        assert_eq!(pairs.len(), 1, "Only first two should overlap on Y axis");
    }
    #[test]
    fn test_sap_z_axis() {
        let sap = SweepAndPrune::new(2);
        let aabbs = vec![
            make_aabb(0.0, 0.0, 0.0, 1.0),
            make_aabb(0.0, 0.0, 1.5, 1.0),
            make_aabb(0.0, 0.0, 10.0, 1.0),
        ];
        let pairs = sap.find_pairs(&aabbs);
        assert_eq!(pairs.len(), 1, "Only first two should overlap on Z axis");
    }
    #[test]
    fn test_bvh_agrees_with_brute_force() {
        let mut aabbs = Vec::new();
        for i in 0..15_i32 {
            let x = (i % 5) as Real * 1.5;
            let y = (i / 5) as Real * 1.5;
            aabbs.push(make_aabb(x, y, 0.0, 1.0));
        }
        let mut bvh = BvhBroadphase::new();
        bvh.build(&aabbs);
        let bvh_pairs = sorted_pairs(bvh.find_pairs(&aabbs));
        let bf_pairs = sorted_pairs(BruteForceBroadPhase.find_pairs(&aabbs));
        assert_eq!(bvh_pairs, bf_pairs, "BVH and brute-force should agree");
    }
    #[test]
    fn test_parallel_sap_agrees_with_serial() {
        let mut aabbs = Vec::new();
        for i in 0..20_i32 {
            let x = (i % 5) as Real * 1.5;
            let y = (i / 5) as Real * 1.5;
            aabbs.push(make_aabb(x, y, 0.0, 1.0));
        }
        let serial_pairs = sorted_pairs(SweepAndPrune::x_axis().find_pairs(&aabbs));
        let parallel_pairs = sorted_pairs(parallel_sap(&aabbs, 0));
        for p in &serial_pairs {
            assert!(
                parallel_pairs.contains(p),
                "parallel SAP missing pair {:?}",
                p
            );
        }
    }
    #[test]
    fn test_parallel_brute_force_agrees_with_serial() {
        let mut aabbs = Vec::new();
        for i in 0..12_i32 {
            aabbs.push(make_aabb(i as Real * 1.5, 0.0, 0.0, 1.0));
        }
        let serial_pairs = sorted_pairs(BruteForceBroadPhase.find_pairs(&aabbs));
        let parallel_pairs = sorted_pairs(parallel_brute_force(&aabbs));
        assert_eq!(serial_pairs, parallel_pairs);
    }
    #[test]
    fn test_broadphase_profiler_basic() {
        let mut prof = BroadphaseProfiler::new();
        let aabbs = vec![
            make_aabb(0.0, 0.0, 0.0, 1.0),
            make_aabb(1.5, 0.0, 0.0, 1.0),
            make_aabb(10.0, 0.0, 0.0, 1.0),
        ];
        let pairs = prof.profile_sap(&aabbs, 0);
        assert_eq!(pairs.len(), 1);
        assert!(prof.last_elapsed_ns() > 0);
        assert_eq!(prof.call_count(), 1);
    }
    #[test]
    fn test_broadphase_profiler_reset() {
        let mut prof = BroadphaseProfiler::new();
        let aabbs = vec![make_aabb(0.0, 0.0, 0.0, 1.0)];
        prof.profile_sap(&aabbs, 0);
        prof.reset();
        assert_eq!(prof.call_count(), 0);
    }
    #[test]
    fn test_gpu_hints_basic() {
        let aabbs = vec![
            make_aabb(0.0, 0.0, 0.0, 1.0),
            make_aabb(1.5, 0.0, 0.0, 1.0),
            make_aabb(10.0, 0.0, 0.0, 1.0),
        ];
        let hints = GpuBroadphaseHints::from_aabbs(&aabbs);
        assert_eq!(hints.n_objects, 3);
        assert!(hints.scene_min[0] <= -1.0);
        assert!(hints.scene_max[0] >= 11.0);
    }
    #[test]
    fn test_gpu_hints_suggested_cell_size() {
        let aabbs = vec![
            make_aabb(0.0, 0.0, 0.0, 0.5),
            make_aabb(2.0, 0.0, 0.0, 0.5),
            make_aabb(4.0, 0.0, 0.0, 0.5),
        ];
        let hints = GpuBroadphaseHints::from_aabbs(&aabbs);
        let cell_size = hints.suggested_cell_size();
        assert!(cell_size > 0.0);
    }
    #[test]
    fn test_gpu_hints_empty() {
        let hints = GpuBroadphaseHints::from_aabbs(&[]);
        assert_eq!(hints.n_objects, 0);
        assert_eq!(hints.suggested_cell_size(), 1.0);
    }
    #[test]
    fn test_scene_graph_add_node() {
        let mut sg = BroadphaseSceneGraph::new();
        let id = sg.add_node(make_aabb(0.0, 0.0, 0.0, 1.0), None);
        assert_eq!(sg.node_count(), 1);
        assert_eq!(id, 0);
    }
    #[test]
    fn test_scene_graph_parent_child() {
        let mut sg = BroadphaseSceneGraph::new();
        let parent = sg.add_node(make_aabb(0.0, 0.0, 0.0, 5.0), None);
        let child1 = sg.add_node(make_aabb(1.0, 0.0, 0.0, 1.0), Some(parent));
        let child2 = sg.add_node(make_aabb(-1.0, 0.0, 0.0, 1.0), Some(parent));
        let children = sg.children_of(parent);
        assert!(children.contains(&child1));
        assert!(children.contains(&child2));
    }
    #[test]
    fn test_scene_graph_find_pairs() {
        let mut sg = BroadphaseSceneGraph::new();
        sg.add_node(make_aabb(0.0, 0.0, 0.0, 1.0), None);
        sg.add_node(make_aabb(1.5, 0.0, 0.0, 1.0), None);
        sg.add_node(make_aabb(10.0, 0.0, 0.0, 1.0), None);
        let pairs = sg.find_pairs_bvh();
        assert_eq!(pairs.len(), 1);
    }
    #[test]
    fn test_scene_graph_remove_node() {
        let mut sg = BroadphaseSceneGraph::new();
        let id = sg.add_node(make_aabb(0.0, 0.0, 0.0, 1.0), None);
        sg.add_node(make_aabb(1.5, 0.0, 0.0, 1.0), None);
        sg.remove_node(id);
        assert_eq!(sg.node_count(), 1);
    }
    #[test]
    fn test_warmstart_cache_stores_pairs() {
        let mut ws = BroadphaseWarmstart::new();
        let pairs = vec![CollisionPair::new(0, 1), CollisionPair::new(2, 3)];
        ws.update(&pairs);
        assert_eq!(ws.cached_pair_count(), 2);
    }
    #[test]
    fn test_warmstart_filter_active() {
        let mut ws = BroadphaseWarmstart::new();
        let aabbs = vec![
            make_aabb(0.0, 0.0, 0.0, 1.0),
            make_aabb(1.5, 0.0, 0.0, 1.0),
            make_aabb(10.0, 0.0, 0.0, 1.0),
        ];
        let initial = vec![CollisionPair::new(0, 1), CollisionPair::new(1, 2)];
        ws.update(&initial);
        let still_active = ws.filter_still_overlapping(&aabbs);
        assert_eq!(still_active.len(), 1, "only (0,1) should still overlap");
        assert_eq!(still_active[0].a, 0);
        assert_eq!(still_active[0].b, 1);
    }
    #[test]
    fn test_warmstart_merge_with_new() {
        let mut ws = BroadphaseWarmstart::new();
        let aabbs = vec![
            make_aabb(0.0, 0.0, 0.0, 1.0),
            make_aabb(1.5, 0.0, 0.0, 1.0),
            make_aabb(10.0, 0.0, 0.0, 1.0),
        ];
        ws.update(&[CollisionPair::new(0, 1)]);
        let new_pairs = vec![CollisionPair::new(0, 2)];
        let merged = ws.merge_with_new(&new_pairs, &aabbs);
        let merged_sorted = {
            let mut v: Vec<(usize, usize)> = merged.iter().map(|p| (p.a, p.b)).collect();
            v.sort_unstable();
            v
        };
        assert!(merged_sorted.contains(&(0, 1)) || merged_sorted.contains(&(0, 2)));
    }
    #[test]
    fn test_warmstart_clear() {
        let mut ws = BroadphaseWarmstart::new();
        ws.update(&[CollisionPair::new(0, 1)]);
        ws.clear();
        assert_eq!(ws.cached_pair_count(), 0);
    }
    #[test]
    fn test_bvh_quality_empty() {
        let bvh = BvhBroadphase::new();
        assert!(
            bvh.quality_metrics().is_none(),
            "empty BVH should have no metrics"
        );
        assert_eq!(bvh.sah_cost(), 0.0);
    }
    #[test]
    fn test_bvh_quality_single_leaf() {
        let mut bvh = BvhBroadphase::new();
        bvh.build(&[make_aabb(0.0, 0.0, 0.0, 1.0)]);
        let q = bvh
            .quality_metrics()
            .expect("single leaf must have metrics");
        assert_eq!(q.leaf_count, 1);
        assert_eq!(q.internal_count, 0);
        assert_eq!(q.sah_cost, 0.0);
    }
    #[test]
    fn test_bvh_quality_multiple_leaves() {
        let aabbs: Vec<Aabb> = (0..6)
            .map(|i| make_aabb(i as Real * 3.0, 0.0, 0.0, 1.0))
            .collect();
        let mut bvh = BvhBroadphase::new();
        bvh.build(&aabbs);
        let q = bvh
            .quality_metrics()
            .expect("non-empty BVH must have metrics");
        assert_eq!(q.leaf_count, 6);
        assert_eq!(q.internal_count, 5, "6 leaves → 5 internal nodes");
        assert!(q.sah_cost >= 0.0);
        assert!(q.avg_leaf_depth >= 0.0);
    }
    #[test]
    fn test_bvh_rebuild_if_degraded_not_triggered() {
        let aabbs: Vec<Aabb> = (0..4)
            .map(|i| make_aabb(i as Real * 3.0, 0.0, 0.0, 1.0))
            .collect();
        let mut bvh = BvhBroadphase::new();
        bvh.build(&aabbs);
        let rebuilt = bvh.rebuild_if_degraded(&aabbs, 1e12);
        assert!(
            !rebuilt,
            "should not rebuild when cost well below threshold"
        );
    }
    #[test]
    fn test_bvh_rebuild_if_degraded_triggers() {
        let aabbs: Vec<Aabb> = (0..8)
            .map(|i| make_aabb(i as Real * 3.0, 0.0, 0.0, 1.0))
            .collect();
        let mut bvh = BvhBroadphase::new();
        bvh.build(&aabbs);
        if bvh.quality_metrics().is_some_and(|q| q.internal_count > 0) {
            let rebuilt = bvh.rebuild_if_degraded(&aabbs, 0.0);
            assert!(
                rebuilt,
                "zero threshold should always trigger rebuild on non-trivial tree"
            );
        }
    }
    #[test]
    fn test_bvh_validate_empty() {
        let bvh = BvhBroadphase::new();
        assert!(bvh.validate(), "empty BVH should be valid");
    }
    #[test]
    fn test_bvh_validate_after_build() {
        let aabbs: Vec<Aabb> = (0..5)
            .map(|i| make_aabb(i as Real * 2.0, 0.0, 0.0, 0.8))
            .collect();
        let mut bvh = BvhBroadphase::new();
        bvh.build(&aabbs);
        assert!(bvh.validate(), "freshly built BVH must pass validation");
    }
    #[test]
    fn test_ray_query_iterative_hits() {
        let aabbs = vec![
            make_aabb(0.0, 0.0, 0.0, 1.0),
            make_aabb(5.0, 0.0, 0.0, 1.0),
            make_aabb(10.0, 0.0, 0.0, 1.0),
        ];
        let mut bvh = BvhBroadphase::new();
        bvh.build(&aabbs);
        let hits_iter =
            bvh.ray_query_iterative(&Vec3::new(-5.0, 0.0, 0.0), &Vec3::new(1.0, 0.0, 0.0), 100.0);
        let mut hits_rec =
            bvh.ray_query(&Vec3::new(-5.0, 0.0, 0.0), &Vec3::new(1.0, 0.0, 0.0), 100.0);
        let mut hits_iter_sorted = hits_iter.clone();
        hits_iter_sorted.sort_unstable();
        hits_rec.sort_unstable();
        assert_eq!(
            hits_iter_sorted, hits_rec,
            "iterative and recursive ray queries must agree"
        );
    }
    #[test]
    fn test_ray_query_iterative_miss() {
        let aabbs = vec![make_aabb(0.0, 0.0, 0.0, 1.0)];
        let mut bvh = BvhBroadphase::new();
        bvh.build(&aabbs);
        let hits =
            bvh.ray_query_iterative(&Vec3::new(0.0, 10.0, 0.0), &Vec3::new(1.0, 0.0, 0.0), 100.0);
        assert!(hits.is_empty(), "ray above all objects should miss");
    }
    fn open_frustum_planes() -> [(Vec3, Real); 6] {
        [
            (Vec3::new(1.0, 0.0, 0.0), -1e9),
            (Vec3::new(-1.0, 0.0, 0.0), -1e9),
            (Vec3::new(0.0, 1.0, 0.0), -1e9),
            (Vec3::new(0.0, -1.0, 0.0), -1e9),
            (Vec3::new(0.0, 0.0, 1.0), -1e9),
            (Vec3::new(0.0, 0.0, -1.0), -1e9),
        ]
    }
    #[test]
    fn test_frustum_query_iterative_all_visible() {
        let aabbs: Vec<Aabb> = (0..4)
            .map(|i| make_aabb(i as Real, 0.0, 0.0, 0.4))
            .collect();
        let mut bvh = BvhBroadphase::new();
        bvh.build(&aabbs);
        let mut visible = bvh.frustum_query_iterative(&open_frustum_planes());
        visible.sort_unstable();
        assert_eq!(
            visible,
            vec![0, 1, 2, 3],
            "open frustum should see all objects"
        );
    }
    #[test]
    fn test_frustum_query_iterative_none_visible() {
        let aabbs: Vec<Aabb> = (0..3)
            .map(|i| make_aabb(i as Real * 2.0, 0.0, 0.0, 0.5))
            .collect();
        let mut bvh = BvhBroadphase::new();
        bvh.build(&aabbs);
        let planes: [(Vec3, Real); 6] = [
            (Vec3::new(1.0, 0.0, 0.0), 1e9),
            (Vec3::new(-1.0, 0.0, 0.0), -1e9),
            (Vec3::new(0.0, 1.0, 0.0), -1e9),
            (Vec3::new(0.0, -1.0, 0.0), -1e9),
            (Vec3::new(0.0, 0.0, 1.0), -1e9),
            (Vec3::new(0.0, 0.0, -1.0), -1e9),
        ];
        let visible = bvh.frustum_query_iterative(&planes);
        assert!(visible.is_empty(), "closed frustum should cull all objects");
    }
    #[test]
    fn test_bvh_batch_insert_and_query() {
        let aabbs: Vec<Aabb> = (0..5)
            .map(|i| make_aabb(i as Real * 3.0, 0.0, 0.0, 1.0))
            .collect();
        let mut bvh = BvhBroadphase::new();
        bvh.batch_insert(&aabbs);
        assert!(bvh.validate(), "BVH after batch_insert must validate");
        let q = bvh
            .quality_metrics()
            .expect("must have metrics after batch_insert");
        assert_eq!(q.leaf_count, 5);
    }
    #[test]
    fn test_dynamic_tree_batch_insert() {
        let mut tree = DynamicAabbTree::new(0.1);
        let aabbs: Vec<Aabb> = (0..4)
            .map(|i| make_aabb(i as Real * 3.0, 0.0, 0.0, 1.0))
            .collect();
        tree.batch_insert(&aabbs);
        assert_eq!(tree.len(), 4);
        assert!(tree.is_dirty(), "tree should be dirty after batch_insert");
        tree.rebuild_if_dirty();
        assert!(!tree.is_dirty(), "tree should not be dirty after rebuild");
        assert!(tree.validate(), "tree must be valid after rebuild");
    }
    #[test]
    fn test_proximity_pairs_touching() {
        let mut tree = DynamicAabbTree::new(0.1);
        tree.insert(make_aabb(0.0, 0.0, 0.0, 1.0));
        tree.insert(make_aabb(2.1, 0.0, 0.0, 1.0));
        let exact_pairs = tree.find_pairs();
        assert!(
            exact_pairs.is_empty(),
            "boxes should not touch with 0 inflation"
        );
        let prox_pairs = tree.find_proximity_pairs(0.2);
        assert_eq!(
            prox_pairs.len(),
            1,
            "boxes should become proximate with inflation 0.2"
        );
    }
    #[test]
    fn test_proximity_pairs_already_overlapping() {
        let mut tree = DynamicAabbTree::new(0.1);
        tree.insert(make_aabb(0.0, 0.0, 0.0, 1.0));
        tree.insert(make_aabb(1.5, 0.0, 0.0, 1.0));
        let prox_pairs = tree.find_proximity_pairs(0.0);
        assert_eq!(
            prox_pairs.len(),
            1,
            "overlapping boxes should show as proximate pairs"
        );
    }
    #[test]
    fn test_dynamic_tree_mark_dirty() {
        let mut tree = DynamicAabbTree::new(0.1);
        tree.insert(make_aabb(0.0, 0.0, 0.0, 1.0));
        tree.rebuild_if_dirty();
        assert!(!tree.is_dirty(), "should not be dirty after rebuild");
        tree.mark_dirty();
        assert!(tree.is_dirty(), "should be dirty after mark_dirty");
    }
    #[test]
    fn test_dynamic_tree_tree_info() {
        let mut tree = DynamicAabbTree::new(0.1);
        tree.insert(make_aabb(0.0, 0.0, 0.0, 1.0));
        tree.insert(make_aabb(5.0, 0.0, 0.0, 1.0));
        let (count, dirty) = tree.tree_info();
        assert_eq!(count, 2);
        assert!(dirty, "tree should be dirty after inserts");
    }
    #[test]
    fn test_dynamic_tree_clear() {
        let mut tree = DynamicAabbTree::new(0.1);
        for i in 0..4 {
            tree.insert(make_aabb(i as Real * 2.0, 0.0, 0.0, 0.5));
        }
        tree.clear();
        assert_eq!(tree.len(), 0);
        assert!(tree.is_empty());
        assert!(!tree.is_dirty(), "cleared tree should not be dirty");
    }
    #[test]
    fn test_aabb_in_frustum_always_inside() {
        let planes = open_frustum_planes();
        let aabb = make_aabb(0.0, 0.0, 0.0, 1.0);
        assert!(
            aabb_in_frustum(&aabb, &planes),
            "open frustum should contain any AABB"
        );
    }
    #[test]
    fn test_aabb_in_frustum_outside() {
        let planes: [(Vec3, Real); 6] = [
            (Vec3::new(1.0, 0.0, 0.0), 1e9),
            (Vec3::new(-1.0, 0.0, 0.0), -1e9),
            (Vec3::new(0.0, 1.0, 0.0), -1e9),
            (Vec3::new(0.0, -1.0, 0.0), -1e9),
            (Vec3::new(0.0, 0.0, 1.0), -1e9),
            (Vec3::new(0.0, 0.0, -1.0), -1e9),
        ];
        let aabb = make_aabb(0.0, 0.0, 0.0, 1.0);
        assert!(
            !aabb_in_frustum(&aabb, &planes),
            "impossible frustum plane should reject AABB"
        );
    }
    #[test]
    fn test_pair_count_histogram_all_dynamic() {
        let aabbs = vec![make_aabb(0.0, 0.0, 0.0, 1.0), make_aabb(1.5, 0.0, 0.0, 1.0)];
        let types = vec![ObjectType::Dynamic, ObjectType::Dynamic];
        let pairs = BruteForceBroadPhase.find_pairs(&aabbs);
        let hist = compute_pair_count_histogram(&pairs, &types);
        assert_eq!(hist.dynamic_dynamic, 1, "one dynamic-dynamic pair expected");
        assert_eq!(hist.total(), 1);
    }
    #[test]
    fn test_pair_count_histogram_static_dynamic() {
        let aabbs = vec![
            make_aabb(0.0, 0.0, 0.0, 1.0),
            make_aabb(1.5, 0.0, 0.0, 1.0),
            make_aabb(10.0, 0.0, 0.0, 1.0),
        ];
        let types = vec![ObjectType::Static, ObjectType::Dynamic, ObjectType::Dynamic];
        let pairs = BruteForceBroadPhase.find_pairs(&aabbs);
        let hist = compute_pair_count_histogram(&pairs, &types);
        assert_eq!(hist.static_dynamic, 1, "one static-dynamic pair (0,1)");
        assert_eq!(hist.static_static, 0);
        assert_eq!(hist.dynamic_dynamic, 0);
    }
    #[test]
    fn test_pair_count_histogram_empty_pairs() {
        let types = vec![ObjectType::Static, ObjectType::Kinematic];
        let hist = compute_pair_count_histogram(&[], &types);
        assert_eq!(hist.total(), 0, "no pairs => all zeros");
    }
    #[test]
    fn test_pair_count_histogram_kinematic() {
        let aabbs = vec![make_aabb(0.0, 0.0, 0.0, 1.0), make_aabb(1.5, 0.0, 0.0, 1.0)];
        let types = vec![ObjectType::Kinematic, ObjectType::Kinematic];
        let pairs = BruteForceBroadPhase.find_pairs(&aabbs);
        let hist = compute_pair_count_histogram(&pairs, &types);
        assert_eq!(hist.kinematic_kinematic, 1);
    }
    #[test]
    fn test_update_batch_moves_objects() {
        let mut tree = DynamicAabbTree::new(0.1);
        tree.insert(make_aabb(0.0, 0.0, 0.0, 1.0));
        tree.insert(make_aabb(5.0, 0.0, 0.0, 1.0));
        update_batch(&mut tree, &[(0, make_aabb(20.0, 0.0, 0.0, 1.0))]);
        let pairs = tree.find_pairs();
        assert!(pairs.is_empty(), "after moving apart, no pairs expected");
    }
    #[test]
    fn test_update_batch_empty() {
        let mut tree = DynamicAabbTree::new(0.1);
        tree.insert(make_aabb(0.0, 0.0, 0.0, 1.0));
        update_batch(&mut tree, &[]);
        assert_eq!(tree.len(), 1, "empty batch should not change tree");
    }
    #[test]
    fn test_refit_bottom_up_nonempty() {
        let aabbs = vec![
            make_aabb(0.0, 0.0, 0.0, 1.0),
            make_aabb(4.0, 0.0, 0.0, 1.0),
            make_aabb(8.0, 0.0, 0.0, 1.0),
            make_aabb(12.0, 0.0, 0.0, 1.0),
        ];
        let mut bvh = BvhBroadphase::new();
        bvh.build(&aabbs);
        let count = refit_bottom_up(&mut bvh);
        assert!(count > 0, "should have refitted at least one internal node");
        assert!(bvh.validate(), "tree must remain valid after refit");
    }
    #[test]
    fn test_refit_bottom_up_empty() {
        let mut bvh = BvhBroadphase::new();
        let count = refit_bottom_up(&mut bvh);
        assert_eq!(count, 0, "empty BVH has nothing to refit");
    }
    #[test]
    fn test_refit_bottom_up_single_leaf() {
        let mut bvh = BvhBroadphase::new();
        bvh.build(&[make_aabb(0.0, 0.0, 0.0, 1.0)]);
        let count = refit_bottom_up(&mut bvh);
        assert_eq!(count, 0, "single-leaf BVH has no internal nodes");
    }
}
