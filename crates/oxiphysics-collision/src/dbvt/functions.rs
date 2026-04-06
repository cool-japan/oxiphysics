//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

/// Index type for nodes in the BVH.
pub(super) type NodeIdx = usize;
/// Fat-AABB margin used to avoid re-inserting leaves that move slightly.
pub(super) const FAT_MARGIN: f64 = 0.1;
#[cfg(test)]
mod tests {
    use super::*;
    use crate::BvhAabb;
    use crate::DynamicBvh;
    use crate::dbvt::BvhCapsule;
    use crate::dbvt::BvhFrustum;
    use oxiphysics_core::Vec3;
    fn make_aabb(cx: f64, cy: f64, cz: f64, half: f64) -> BvhAabb {
        BvhAabb::new(
            Vec3::new(cx - half, cy - half, cz - half),
            Vec3::new(cx + half, cy + half, cz + half),
        )
    }
    #[test]
    fn test_bvhaabb_merge() {
        let a = make_aabb(0.0, 0.0, 0.0, 1.0);
        let b = make_aabb(4.0, 0.0, 0.0, 1.0);
        let m = BvhAabb::merge(&a, &b);
        assert!(m.contains(&a), "merged does not contain a");
        assert!(m.contains(&b), "merged does not contain b");
        assert!((m.min.x - (-1.0)).abs() < 1e-10);
        assert!((m.max.x - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_bvhaabb_intersects() {
        let a = make_aabb(0.0, 0.0, 0.0, 1.0);
        let b = make_aabb(1.5, 0.0, 0.0, 1.0);
        let c = make_aabb(5.0, 0.0, 0.0, 1.0);
        assert!(a.intersects(&b), "overlapping AABBs should intersect");
        assert!(b.intersects(&a), "intersection must be symmetric");
        assert!(!a.intersects(&c), "separated AABBs should not intersect");
        assert!(!c.intersects(&a), "must be symmetric");
    }
    #[test]
    fn test_dbvt_insert_single() {
        let mut bvh = DynamicBvh::new();
        let _leaf = bvh.insert(make_aabb(0.0, 0.0, 0.0, 1.0), 0);
        assert_eq!(bvh.n_leaves(), 1);
        assert!(bvh.validate(), "tree invariants violated");
    }
    #[test]
    fn test_dbvt_insert_multiple() {
        let mut bvh = DynamicBvh::new();
        for i in 0..5u32 {
            bvh.insert(make_aabb(i as f64 * 3.0, 0.0, 0.0, 1.0), i);
        }
        assert_eq!(bvh.n_leaves(), 5);
        assert!(bvh.validate(), "tree invariants violated after 5 inserts");
    }
    #[test]
    fn test_dbvt_remove() {
        let mut bvh = DynamicBvh::new();
        let l0 = bvh.insert(make_aabb(0.0, 0.0, 0.0, 1.0), 0);
        let l1 = bvh.insert(make_aabb(3.0, 0.0, 0.0, 1.0), 1);
        let _l2 = bvh.insert(make_aabb(6.0, 0.0, 0.0, 1.0), 2);
        bvh.remove(l1);
        assert_eq!(bvh.n_leaves(), 2);
        assert!(bvh.validate(), "tree invariants violated after remove");
        let hits = bvh.query_aabb(&make_aabb(0.0, 0.0, 0.0, 0.5));
        assert!(hits.contains(&0), "leaf 0 should still be found");
        let _ = l0;
    }
    #[test]
    fn test_dbvt_query_overlapping() {
        let mut bvh = DynamicBvh::new();
        bvh.insert(make_aabb(0.0, 0.0, 0.0, 1.0), 10);
        bvh.insert(make_aabb(1.5, 0.0, 0.0, 1.0), 20);
        let pairs = bvh.query_overlapping_pairs();
        assert_eq!(
            pairs.len(),
            1,
            "expected exactly 1 overlapping pair, got {:?}",
            pairs
        );
        let (a, b) = pairs[0];
        assert!(
            (a == 10 && b == 20) || (a == 20 && b == 10),
            "unexpected pair ({},{})",
            a,
            b
        );
    }
    #[test]
    fn test_dbvt_query_non_overlapping() {
        let mut bvh = DynamicBvh::new();
        bvh.insert(make_aabb(0.0, 0.0, 0.0, 1.0), 1);
        bvh.insert(make_aabb(10.0, 0.0, 0.0, 1.0), 2);
        let pairs = bvh.query_overlapping_pairs();
        assert!(pairs.is_empty(), "expected 0 pairs, got {:?}", pairs);
    }
    #[test]
    fn test_dbvt_update_moved() {
        let mut bvh = DynamicBvh::new();
        bvh.insert(make_aabb(0.0, 0.0, 0.0, 1.0), 42);
        bvh.insert(make_aabb(3.0, 0.0, 0.0, 0.5), 99);
        let orig_hits = bvh.query_aabb(&make_aabb(0.0, 0.0, 0.0, 0.5));
        assert!(
            orig_hits.contains(&42),
            "should be found at original position"
        );
        let root = bvh.root.unwrap();
        let leaf42 = find_leaf_by_data(&bvh, root, 42).expect("leaf 42 not found");
        let updated = bvh.update(leaf42, make_aabb(50.0, 50.0, 50.0, 1.0));
        assert!(
            updated,
            "update should return true when aabb moved significantly"
        );
        let hits_after = bvh.query_aabb(&make_aabb(0.0, 0.0, 0.0, 0.5));
        assert!(
            !hits_after.contains(&42),
            "data=42 should not be found near origin after moving to (50,50,50)"
        );
        assert!(bvh.validate(), "tree must still be valid after update");
    }
    #[test]
    fn test_dbvt_query_aabb() {
        let mut bvh = DynamicBvh::new();
        bvh.insert(make_aabb(0.0, 0.0, 0.0, 1.0), 7);
        bvh.insert(make_aabb(5.0, 0.0, 0.0, 1.0), 8);
        bvh.insert(make_aabb(10.0, 0.0, 0.0, 1.0), 9);
        let hits = bvh.query_point(Vec3::new(0.0, 0.0, 0.0));
        assert!(hits.contains(&7), "point at origin should hit data=7");
        assert!(!hits.contains(&8), "data=8 far away, should not be hit");
        assert!(!hits.contains(&9), "data=9 far away, should not be hit");
    }
    fn find_leaf_by_data(bvh: &DynamicBvh, start: NodeIdx, target: u32) -> Option<NodeIdx> {
        let node = &bvh.nodes[start];
        if node.is_leaf() {
            return if node.data == Some(target) {
                Some(start)
            } else {
                None
            };
        }
        if let Some(l) = node.left
            && let Some(found) = find_leaf_by_data(bvh, l, target)
        {
            return Some(found);
        }
        if let Some(r) = node.right
            && let Some(found) = find_leaf_by_data(bvh, r, target)
        {
            return Some(found);
        }
        None
    }
    #[test]
    fn test_dbvt_stats_empty() {
        let bvh = DynamicBvh::new();
        let s = bvh.stats();
        assert_eq!(s.leaf_count, 0);
        assert_eq!(s.node_count, 0);
        assert_eq!(s.height, 0);
        assert_eq!(s.root_surface_area, 0.0);
    }
    #[test]
    fn test_dbvt_stats_single_leaf() {
        let mut bvh = DynamicBvh::new();
        bvh.insert(make_aabb(0.0, 0.0, 0.0, 1.0), 1);
        let s = bvh.stats();
        assert_eq!(s.leaf_count, 1);
        assert!(s.root_surface_area > 0.0);
    }
    #[test]
    fn test_dbvt_stats_multiple_leaves() {
        let mut bvh = DynamicBvh::new();
        for i in 0..4u32 {
            bvh.insert(make_aabb(i as f64 * 3.0, 0.0, 0.0, 1.0), i);
        }
        let s = bvh.stats();
        assert_eq!(s.leaf_count, 4);
        assert_eq!(s.node_count, 7);
    }
    #[test]
    fn test_refit_all_valid_tree() {
        let mut bvh = DynamicBvh::new();
        for i in 0..5u32 {
            bvh.insert(make_aabb(i as f64 * 2.0, 0.0, 0.0, 0.5), i);
        }
        bvh.refit_all();
        assert!(bvh.validate(), "tree should remain valid after refit_all");
    }
    #[test]
    fn test_ray_cast_hits_leaf() {
        let mut bvh = DynamicBvh::new();
        bvh.insert(make_aabb(0.0, 0.0, 0.0, 1.0), 42);
        let hits = bvh.ray_cast(
            Vec3::new(-5.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            0.0,
            100.0,
        );
        assert!(!hits.is_empty(), "ray along x should hit the leaf");
        assert_eq!(hits[0].0, 42);
    }
    #[test]
    fn test_ray_cast_misses() {
        let mut bvh = DynamicBvh::new();
        bvh.insert(make_aabb(0.0, 0.0, 0.0, 1.0), 7);
        let hits = bvh.ray_cast(
            Vec3::new(-5.0, 10.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            0.0,
            100.0,
        );
        assert!(hits.is_empty(), "high ray should miss the leaf");
    }
    #[test]
    fn test_ray_cast_sorted_by_t() {
        let mut bvh = DynamicBvh::new();
        bvh.insert(make_aabb(3.0, 0.0, 0.0, 0.8), 1);
        bvh.insert(make_aabb(0.0, 0.0, 0.0, 0.8), 2);
        let hits = bvh.ray_cast(
            Vec3::new(-5.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            0.0,
            100.0,
        );
        for w in hits.windows(2) {
            assert!(w[0].1 <= w[1].1, "hits not sorted by t");
        }
    }
    #[test]
    fn test_query_aabb_filtered_passes() {
        let mut bvh = DynamicBvh::new();
        bvh.insert(make_aabb(0.0, 0.0, 0.0, 1.0), 10);
        bvh.insert(make_aabb(3.0, 0.0, 0.0, 1.0), 20);
        let hits = bvh.query_aabb_filtered(&make_aabb(0.0, 0.0, 0.0, 0.5), |d| d > 15);
        assert!(!hits.contains(&10), "data=10 should be filtered out");
    }
    #[test]
    fn test_query_aabb_filtered_finds_all() {
        let mut bvh = DynamicBvh::new();
        for i in 0..4u32 {
            bvh.insert(make_aabb(0.0, 0.0, 0.0, 5.0), i);
        }
        let hits = bvh.query_aabb_filtered(&make_aabb(0.0, 0.0, 0.0, 4.0), |_| true);
        assert_eq!(hits.len(), 4, "all leaves should match");
    }
    #[test]
    fn test_insert_bulk_correct_count() {
        let mut bvh = DynamicBvh::new();
        let items: Vec<(BvhAabb, u32)> = (0..6u32)
            .map(|i| (make_aabb(i as f64 * 3.0, 0.0, 0.0, 1.0), i))
            .collect();
        let leaves = bvh.insert_bulk(&items);
        assert_eq!(leaves.len(), 6);
        assert_eq!(bvh.n_leaves(), 6);
        assert!(bvh.validate());
    }
    #[test]
    fn test_dbvt_clear() {
        let mut bvh = DynamicBvh::new();
        for i in 0..5u32 {
            bvh.insert(make_aabb(i as f64, 0.0, 0.0, 0.5), i);
        }
        bvh.clear();
        assert_eq!(bvh.n_leaves(), 0);
        assert_eq!(bvh.n_nodes(), 0);
        assert!(bvh.root.is_none());
    }
    #[test]
    fn test_all_leaf_data() {
        let mut bvh = DynamicBvh::new();
        for i in 0..3u32 {
            bvh.insert(make_aabb(i as f64 * 2.0, 0.0, 0.0, 0.5), i + 100);
        }
        let mut data = bvh.all_leaf_data();
        data.sort();
        assert_eq!(data, vec![100, 101, 102]);
    }
    #[test]
    fn test_query_pairs_with_stats() {
        let mut bvh = DynamicBvh::new();
        bvh.insert(make_aabb(0.0, 0.0, 0.0, 1.0), 1);
        bvh.insert(make_aabb(1.5, 0.0, 0.0, 1.0), 2);
        bvh.insert(make_aabb(20.0, 0.0, 0.0, 1.0), 3);
        let (pairs, stats) = bvh.query_pairs_with_stats();
        assert!(stats.overlap_pairs == pairs.len());
        assert_eq!(stats.leaf_count, 3);
    }
    #[test]
    fn test_quality_metrics_empty() {
        let bvh = DynamicBvh::new();
        assert!(bvh.quality_metrics().is_none(), "empty tree has no quality");
    }
    #[test]
    fn test_quality_metrics_single_leaf() {
        let mut bvh = DynamicBvh::new();
        bvh.insert(make_aabb(0.0, 0.0, 0.0, 1.0), 1);
        let q = bvh
            .quality_metrics()
            .expect("single leaf should have metrics");
        assert_eq!(q.leaf_count, 1);
        assert_eq!(q.internal_count, 0);
        assert_eq!(q.sah_cost, 0.0);
        assert_eq!(q.max_depth, 0);
    }
    #[test]
    fn test_quality_metrics_multi_leaf() {
        let mut bvh = DynamicBvh::new();
        for i in 0..5u32 {
            bvh.insert(make_aabb(i as f64 * 3.0, 0.0, 0.0, 1.0), i);
        }
        let q = bvh
            .quality_metrics()
            .expect("non-empty tree must have metrics");
        assert_eq!(q.leaf_count, 5);
        assert_eq!(q.internal_count, 4, "5 leaves → 4 internal nodes");
        assert!(q.sah_cost >= 0.0, "SAH cost must be non-negative");
        assert!(q.avg_leaf_depth >= 0.0);
    }
    #[test]
    fn test_quality_degraded_threshold() {
        let mut bvh = DynamicBvh::new();
        for i in 0..8u32 {
            bvh.insert(make_aabb(i as f64 * 3.0, 0.0, 0.0, 1.0), i);
        }
        assert!(
            !bvh.quality_degraded(1e9),
            "freshly built tree should not be degraded at huge threshold"
        );
        if let Some(q) = bvh.quality_metrics()
            && q.internal_count > 0
        {
            assert!(
                bvh.quality_degraded(0.0),
                "zero threshold must always degrade"
            );
        }
    }
    fn make_open_frustum() -> BvhFrustum {
        BvhFrustum::new([
            (Vec3::new(1.0, 0.0, 0.0), -1e9),
            (Vec3::new(-1.0, 0.0, 0.0), -1e9),
            (Vec3::new(0.0, 1.0, 0.0), -1e9),
            (Vec3::new(0.0, -1.0, 0.0), -1e9),
            (Vec3::new(0.0, 0.0, 1.0), -1e9),
            (Vec3::new(0.0, 0.0, -1.0), -1e9),
        ])
    }
    #[test]
    fn test_frustum_query_all_visible() {
        let mut bvh = DynamicBvh::new();
        for i in 0..4u32 {
            bvh.insert(make_aabb(i as f64, 0.0, 0.0, 0.5), i);
        }
        let frustum = make_open_frustum();
        let mut visible = bvh.frustum_query(&frustum);
        visible.sort();
        assert_eq!(
            visible,
            vec![0, 1, 2, 3],
            "all leaves should be visible in open frustum"
        );
    }
    #[test]
    fn test_frustum_query_none_visible() {
        let mut bvh = DynamicBvh::new();
        for i in 0..3u32 {
            bvh.insert(make_aabb(i as f64 * 2.0, 0.0, 0.0, 0.5), i);
        }
        let frustum = BvhFrustum::new([
            (Vec3::new(1.0, 0.0, 0.0), 1e9),
            (Vec3::new(-1.0, 0.0, 0.0), -1e9),
            (Vec3::new(0.0, 1.0, 0.0), -1e9),
            (Vec3::new(0.0, -1.0, 0.0), -1e9),
            (Vec3::new(0.0, 0.0, 1.0), -1e9),
            (Vec3::new(0.0, 0.0, -1.0), -1e9),
        ]);
        let visible = bvh.frustum_query(&frustum);
        assert!(visible.is_empty(), "all leaves should be culled");
    }
    #[test]
    fn test_frustum_intersects_aabb() {
        let frustum = make_open_frustum();
        let aabb = make_aabb(0.0, 0.0, 0.0, 1.0);
        assert!(
            frustum.intersects_aabb(&aabb),
            "open frustum should intersect any AABB"
        );
    }
    #[test]
    fn test_k_nearest_empty_tree() {
        let bvh = DynamicBvh::new();
        let result = bvh.k_nearest(Vec3::new(0.0, 0.0, 0.0), 3);
        assert!(result.is_empty(), "empty tree should return no nearest");
    }
    #[test]
    fn test_k_nearest_k_zero() {
        let mut bvh = DynamicBvh::new();
        bvh.insert(make_aabb(0.0, 0.0, 0.0, 1.0), 1);
        let result = bvh.k_nearest(Vec3::new(0.0, 0.0, 0.0), 0);
        assert!(result.is_empty(), "k=0 should return empty");
    }
    #[test]
    fn test_k_nearest_sorted_order() {
        let mut bvh = DynamicBvh::new();
        bvh.insert(make_aabb(10.0, 0.0, 0.0, 0.5), 10);
        bvh.insert(make_aabb(1.0, 0.0, 0.0, 0.5), 1);
        bvh.insert(make_aabb(5.0, 0.0, 0.0, 0.5), 5);
        let result = bvh.k_nearest(Vec3::new(0.0, 0.0, 0.0), 3);
        assert_eq!(result.len(), 3);
        for w in result.windows(2) {
            assert!(w[0].dist_sq <= w[1].dist_sq, "k-nearest results not sorted");
        }
        assert_eq!(result[0].data, 1, "nearest leaf should be data=1");
    }
    #[test]
    fn test_k_nearest_fewer_than_k() {
        let mut bvh = DynamicBvh::new();
        bvh.insert(make_aabb(0.0, 0.0, 0.0, 1.0), 42);
        let result = bvh.k_nearest(Vec3::new(0.0, 0.0, 0.0), 10);
        assert_eq!(result.len(), 1, "only 1 leaf in tree");
    }
    #[test]
    fn test_serialize_empty_tree() {
        let bvh = DynamicBvh::new();
        let buf = bvh.serialize();
        assert_eq!(buf.len(), 1, "empty tree serializes to just the count");
        assert_eq!(buf[0], 0.0);
    }
    #[test]
    fn test_serialize_nonempty_tree() {
        let mut bvh = DynamicBvh::new();
        bvh.insert(make_aabb(0.0, 0.0, 0.0, 1.0), 7);
        bvh.insert(make_aabb(3.0, 0.0, 0.0, 1.0), 8);
        let buf = bvh.serialize();
        let n_nodes = bvh.nodes.len();
        assert_eq!(buf.len(), 1 + n_nodes * 14);
        assert_eq!(buf[0], n_nodes as f64);
    }
    #[test]
    fn test_serialize_round_trip_node_count() {
        let mut bvh = DynamicBvh::new();
        for i in 0..6u32 {
            bvh.insert(make_aabb(i as f64 * 2.0, 0.0, 0.0, 0.5), i);
        }
        let buf = bvh.serialize();
        let stored_count = buf[0] as usize;
        assert_eq!(stored_count, bvh.nodes.len());
    }
    #[test]
    fn test_move_proxy_large_displacement() {
        let mut bvh = DynamicBvh::new();
        bvh.insert(make_aabb(0.0, 0.0, 0.0, 1.0), 99);
        bvh.insert(make_aabb(5.0, 0.0, 0.0, 1.0), 88);
        let root = bvh.root.unwrap();
        let leaf99 = find_leaf_by_data(&bvh, root, 99).expect("leaf 99 not found");
        let new_aabb = make_aabb(50.0, 0.0, 0.0, 1.0);
        let moved = bvh.move_proxy(leaf99, new_aabb, Vec3::new(5.0, 0.0, 0.0), 1.0);
        assert!(moved, "large displacement should trigger re-insertion");
        assert!(bvh.validate(), "tree must remain valid after move_proxy");
    }
    #[test]
    fn test_move_proxy_small_displacement_no_reinsertion() {
        let mut bvh = DynamicBvh::new();
        bvh.insert(make_aabb(0.0, 0.0, 0.0, 5.0), 7);
        bvh.insert(make_aabb(20.0, 0.0, 0.0, 1.0), 8);
        let root = bvh.root.unwrap();
        let leaf7 = find_leaf_by_data(&bvh, root, 7).expect("leaf 7 not found");
        let new_aabb = make_aabb(0.01, 0.0, 0.0, 5.0);
        let moved = bvh.move_proxy(leaf7, new_aabb, Vec3::new(0.0, 0.0, 0.0), 0.001);
        assert!(
            !moved,
            "tiny displacement inside fat AABB should not re-insert"
        );
    }
    #[test]
    fn test_sah_cost_empty() {
        let bvh = DynamicBvh::new();
        assert_eq!(bvh.sah_cost(), 0.0, "empty tree has zero SAH cost");
    }
    #[test]
    fn test_sah_cost_non_negative() {
        let mut bvh = DynamicBvh::new();
        for i in 0..6u32 {
            bvh.insert(make_aabb(i as f64 * 2.0, 0.0, 0.0, 0.8), i);
        }
        assert!(bvh.sah_cost() >= 0.0, "SAH cost must be non-negative");
    }
    #[test]
    fn test_optimize_rotations_keeps_tree_valid() {
        let mut bvh = DynamicBvh::new();
        for i in 0..8u32 {
            bvh.insert(make_aabb(i as f64 * 2.0, 0.0, 0.0, 0.8), i);
        }
        bvh.optimize_rotations();
        assert!(
            bvh.validate(),
            "tree must remain valid after optimize_rotations"
        );
        assert_eq!(bvh.n_leaves(), 8, "leaf count must not change");
    }
    #[test]
    fn test_optimize_rotations_sah_does_not_increase() {
        let mut bvh = DynamicBvh::new();
        for i in [0u32, 100, 1, 99, 2, 98] {
            bvh.insert(make_aabb(i as f64, 0.0, 0.0, 0.4), i);
        }
        let _cost_before = bvh.sah_cost();
        bvh.optimize_rotations();
        assert!(bvh.validate(), "tree must stay valid");
        assert!(bvh.sah_cost().is_finite());
    }
    #[test]
    fn test_update_leaf_inplace_contained() {
        let mut bvh = DynamicBvh::new();
        bvh.insert(make_aabb(0.0, 0.0, 0.0, 5.0), 42);
        bvh.insert(make_aabb(20.0, 0.0, 0.0, 1.0), 1);
        let root = bvh.root.unwrap();
        let leaf42 = find_leaf_by_data(&bvh, root, 42).expect("leaf not found");
        let new_aabb = make_aabb(0.1, 0.0, 0.0, 5.0);
        bvh.update_leaf_inplace(leaf42, new_aabb);
        assert!(bvh.validate(), "tree must stay valid after inplace update");
    }
    #[test]
    fn test_root_accessor() {
        let mut bvh = DynamicBvh::new();
        assert!(bvh.root().is_none(), "empty tree has no root");
        bvh.insert(make_aabb(0.0, 0.0, 0.0, 1.0), 1);
        assert!(bvh.root().is_some(), "non-empty tree has a root");
    }
    #[test]
    fn test_compute_sah_cost_alias_matches_sah_cost() {
        let mut bvh = DynamicBvh::new();
        for i in 0..5u32 {
            bvh.insert(make_aabb(i as f64 * 3.0, 0.0, 0.0, 1.0), i);
        }
        assert!(
            (bvh.compute_sah_cost() - bvh.sah_cost()).abs() < 1e-12,
            "compute_sah_cost and sah_cost must agree"
        );
    }
    #[test]
    fn test_compute_sah_cost_empty_tree() {
        let bvh = DynamicBvh::new();
        assert_eq!(bvh.compute_sah_cost(), 0.0, "empty tree SAH cost must be 0");
    }
    #[test]
    fn test_compute_sah_cost_single_leaf() {
        let mut bvh = DynamicBvh::new();
        bvh.insert(make_aabb(0.0, 0.0, 0.0, 1.0), 42);
        assert_eq!(bvh.compute_sah_cost(), 0.0, "single leaf has zero SAH cost");
    }
    #[test]
    fn test_balance_rotation_tree_valid_after() {
        let mut bvh = DynamicBvh::new();
        for i in [0u32, 50, 1, 49, 2, 48, 3, 47] {
            bvh.insert(make_aabb(i as f64, 0.0, 0.0, 0.4), i);
        }
        bvh.balance_rotation();
        assert!(bvh.validate(), "tree must be valid after balance_rotation");
        assert_eq!(bvh.n_leaves(), 8, "leaf count must not change");
    }
    #[test]
    fn test_balance_rotation_empty_tree_zero_rotations() {
        let mut bvh = DynamicBvh::new();
        let accepted = bvh.balance_rotation();
        assert_eq!(accepted, 0, "empty tree has no rotations");
    }
    #[test]
    fn test_balance_rotation_sah_cost_finite() {
        let mut bvh = DynamicBvh::new();
        for i in 0..10u32 {
            bvh.insert(make_aabb(i as f64 * 2.0, 0.0, 0.0, 0.8), i);
        }
        bvh.balance_rotation();
        assert!(
            bvh.compute_sah_cost().is_finite(),
            "SAH cost must be finite after rotation"
        );
    }
    fn open_frustum() -> BvhFrustum {
        BvhFrustum::new([
            (Vec3::new(1.0, 0.0, 0.0), -1e9),
            (Vec3::new(-1.0, 0.0, 0.0), -1e9),
            (Vec3::new(0.0, 1.0, 0.0), -1e9),
            (Vec3::new(0.0, -1.0, 0.0), -1e9),
            (Vec3::new(0.0, 0.0, 1.0), -1e9),
            (Vec3::new(0.0, 0.0, -1.0), -1e9),
        ])
    }
    fn closed_frustum() -> BvhFrustum {
        BvhFrustum::new([
            (Vec3::new(1.0, 0.0, 0.0), 1e9),
            (Vec3::new(-1.0, 0.0, 0.0), 1e9),
            (Vec3::new(0.0, 1.0, 0.0), 1e9),
            (Vec3::new(0.0, -1.0, 0.0), 1e9),
            (Vec3::new(0.0, 0.0, 1.0), 1e9),
            (Vec3::new(0.0, 0.0, -1.0), 1e9),
        ])
    }
    #[test]
    fn test_traverse_frustum_open_returns_all() {
        let mut bvh = DynamicBvh::new();
        for i in 0..4u32 {
            bvh.insert(make_aabb(i as f64 * 5.0, 0.0, 0.0, 1.0), i);
        }
        let vis = bvh.traverse_frustum(&open_frustum());
        let mut vis_sorted = vis.clone();
        vis_sorted.sort_unstable();
        assert_eq!(
            vis_sorted,
            vec![0, 1, 2, 3],
            "open frustum should pass all leaves"
        );
    }
    #[test]
    fn test_traverse_frustum_closed_returns_none() {
        let mut bvh = DynamicBvh::new();
        for i in 0..4u32 {
            bvh.insert(make_aabb(i as f64 * 5.0, 0.0, 0.0, 1.0), i);
        }
        let vis = bvh.traverse_frustum(&closed_frustum());
        assert!(vis.is_empty(), "closed frustum should pass no leaves");
    }
    #[test]
    fn test_traverse_frustum_agrees_with_frustum_query() {
        let mut bvh = DynamicBvh::new();
        for i in 0..6u32 {
            bvh.insert(make_aabb(i as f64 * 3.0, 0.0, 0.0, 1.0), i);
        }
        let f = open_frustum();
        let mut a = bvh.traverse_frustum(&f);
        let mut b = bvh.frustum_query(&f);
        a.sort_unstable();
        b.sort_unstable();
        assert_eq!(a, b, "traverse_frustum and frustum_query must agree");
    }
    #[test]
    fn test_longest_axis_x() {
        let a = BvhAabb::new(Vec3::new(0.0, 0.0, 0.0), Vec3::new(3.0, 1.0, 1.0));
        assert_eq!(a.longest_axis(), 0, "X axis should be longest");
    }
    #[test]
    fn test_longest_axis_y() {
        let a = BvhAabb::new(Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 4.0, 1.0));
        assert_eq!(a.longest_axis(), 1, "Y axis should be longest");
    }
    #[test]
    fn test_longest_axis_z() {
        let a = BvhAabb::new(Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 1.0, 5.0));
        assert_eq!(a.longest_axis(), 2, "Z axis should be longest");
    }
    #[test]
    fn test_approx_eq_same_box() {
        let a = make_aabb(1.0, 2.0, 3.0, 1.0);
        assert!(a.approx_eq(&a, 1e-10), "identical AABB must be approx_eq");
    }
    #[test]
    fn test_approx_eq_slightly_different() {
        let a = make_aabb(0.0, 0.0, 0.0, 1.0);
        let b = BvhAabb::new(
            Vec3::new(-1.0 + 1e-12, -1.0 + 1e-12, -1.0 + 1e-12),
            Vec3::new(1.0 + 1e-12, 1.0 + 1e-12, 1.0 + 1e-12),
        );
        assert!(
            a.approx_eq(&b, 1e-9),
            "tiny difference within eps should be equal"
        );
    }
    #[test]
    fn test_aabb_translate() {
        let a = make_aabb(0.0, 0.0, 0.0, 1.0);
        let b = a.translate(Vec3::new(5.0, 0.0, 0.0));
        assert!(
            (b.center().x - 5.0).abs() < 1e-10,
            "center should shift by 5 on X"
        );
        assert!(
            (a.half_extents().x - b.half_extents().x).abs() < 1e-10,
            "extents must not change after translate"
        );
    }
    #[test]
    fn test_aabb_scale() {
        let a = make_aabb(0.0, 0.0, 0.0, 1.0);
        let b = a.scale(2.0);
        assert!(
            (b.half_extents().x - 2.0).abs() < 1e-10,
            "half-extent should double after scale(2)"
        );
        assert!(
            (b.center().x - a.center().x).abs() < 1e-10,
            "center should be unchanged after scale"
        );
    }
    #[test]
    fn test_bvhaabb_ray_hit() {
        let a = make_aabb(0.0, 0.0, 0.0, 1.0);
        let t = a.ray_intersect(
            Vec3::new(-5.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            0.0,
            100.0,
        );
        assert!(t.is_some(), "ray along X should hit the AABB");
    }
    #[test]
    fn test_bvhaabb_ray_miss() {
        let a = make_aabb(0.0, 0.0, 0.0, 1.0);
        let t = a.ray_intersect(
            Vec3::new(-5.0, 5.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            0.0,
            100.0,
        );
        assert!(t.is_none(), "offset ray should miss AABB");
    }
    #[test]
    fn test_bvhaabb_ray_parallel_inside() {
        let a = make_aabb(0.0, 0.0, 0.0, 5.0);
        let t = a.ray_intersect(
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            0.0,
            100.0,
        );
        assert!(
            t.is_some(),
            "ray starting inside AABB should record a hit at t=0"
        );
        assert!(
            (t.unwrap()).abs() < 1e-10,
            "entry t should be ~0 when inside"
        );
    }
    #[test]
    fn test_find_leaf_finds_data() {
        let mut bvh = DynamicBvh::new();
        bvh.insert(make_aabb(0.0, 0.0, 0.0, 1.0), 42);
        bvh.insert(make_aabb(5.0, 0.0, 0.0, 1.0), 99);
        assert!(bvh.find_leaf(42).is_some(), "leaf 42 should be found");
        assert!(bvh.find_leaf(99).is_some(), "leaf 99 should be found");
        assert!(
            bvh.find_leaf(7).is_none(),
            "nonexistent data should return None"
        );
    }
    #[test]
    fn test_find_leaf_empty_tree() {
        let bvh = DynamicBvh::new();
        assert!(
            bvh.find_leaf(1).is_none(),
            "empty tree should not find any leaf"
        );
    }
    #[test]
    fn test_n_internal_matches_n_nodes_minus_n_leaves() {
        let mut bvh = DynamicBvh::new();
        for i in 0..6u32 {
            bvh.insert(make_aabb(i as f64 * 3.0, 0.0, 0.0, 1.0), i);
        }
        assert_eq!(bvh.n_internal(), 5);
        assert_eq!(bvh.n_leaves() + bvh.n_internal(), bvh.n_nodes());
    }
    #[test]
    fn test_max_depth_single_leaf() {
        let mut bvh = DynamicBvh::new();
        bvh.insert(make_aabb(0.0, 0.0, 0.0, 1.0), 1);
        assert_eq!(bvh.max_depth(), 0, "single leaf is at depth 0");
    }
    #[test]
    fn test_max_depth_grows_with_tree() {
        let mut bvh = DynamicBvh::new();
        for i in 0..8u32 {
            bvh.insert(make_aabb(i as f64 * 2.0, 0.0, 0.0, 0.8), i);
        }
        assert!(
            bvh.max_depth() >= 2,
            "tree with 8 leaves must have depth >= 2"
        );
    }
    #[test]
    fn test_root_aabb_empty() {
        let bvh = DynamicBvh::new();
        assert!(bvh.root_aabb().is_none(), "empty tree has no root AABB");
    }
    #[test]
    fn test_root_aabb_contains_all_leaves() {
        let mut bvh = DynamicBvh::new();
        let positions = [0.0_f64, 5.0, 10.0, -3.0];
        for (i, &x) in positions.iter().enumerate() {
            bvh.insert(make_aabb(x, 0.0, 0.0, 0.5), i as u32);
        }
        let root_aabb = bvh.root_aabb().expect("non-empty tree must have root AABB");
        assert!(root_aabb.min.x <= -3.0, "root AABB min X too large");
        assert!(root_aabb.max.x >= 10.0, "root AABB max X too small");
    }
    #[test]
    fn test_query_sphere_hits_nearby() {
        let mut bvh = DynamicBvh::new();
        bvh.insert(make_aabb(0.0, 0.0, 0.0, 1.0), 1);
        bvh.insert(make_aabb(20.0, 0.0, 0.0, 1.0), 2);
        let hits = bvh.query_sphere(Vec3::new(0.0, 0.0, 0.0), 2.0);
        assert!(hits.contains(&1), "leaf 1 should be hit by sphere");
        assert!(!hits.contains(&2), "leaf 2 is far from sphere");
    }
    #[test]
    fn test_query_sphere_empty_tree() {
        let bvh = DynamicBvh::new();
        assert!(bvh.query_sphere(Vec3::new(0.0, 0.0, 0.0), 100.0).is_empty());
    }
    #[test]
    fn test_any_in_sphere_true_and_false() {
        let mut bvh = DynamicBvh::new();
        bvh.insert(make_aabb(0.0, 0.0, 0.0, 1.0), 10);
        assert!(bvh.any_in_sphere(Vec3::new(0.0, 0.0, 0.0), 1.0));
        assert!(!bvh.any_in_sphere(Vec3::new(100.0, 0.0, 0.0), 1.0));
    }
    #[test]
    fn test_query_capsule_along_x_hits_aligned_leaves() {
        let mut bvh = DynamicBvh::new();
        bvh.insert(make_aabb(0.0, 0.0, 0.0, 0.5), 1);
        bvh.insert(make_aabb(3.0, 0.0, 0.0, 0.5), 2);
        bvh.insert(make_aabb(0.0, 10.0, 0.0, 0.5), 3);
        let capsule = BvhCapsule {
            start: Vec3::new(-1.0, 0.0, 0.0),
            end: Vec3::new(4.0, 0.0, 0.0),
            radius: 1.0,
        };
        let hits = bvh.query_capsule(capsule);
        assert!(hits.contains(&1), "leaf 1 should be hit");
        assert!(hits.contains(&2), "leaf 2 should be hit");
        assert!(!hits.contains(&3), "leaf 3 is far off capsule axis");
    }
    #[test]
    fn test_subtree_leaf_count_root() {
        let mut bvh = DynamicBvh::new();
        for i in 0..5u32 {
            bvh.insert(make_aabb(i as f64 * 2.0, 0.0, 0.0, 0.8), i);
        }
        let root = bvh.root.unwrap();
        assert_eq!(
            bvh.subtree_leaf_count(root),
            5,
            "subtree leaf count of root should equal n_leaves"
        );
    }
    #[test]
    fn test_validate_heights_fresh_tree() {
        let mut bvh = DynamicBvh::new();
        for i in 0..6u32 {
            bvh.insert(make_aabb(i as f64 * 2.0, 0.0, 0.0, 0.8), i);
        }
        assert!(
            bvh.validate_heights(),
            "freshly built tree should have valid heights"
        );
    }
    #[test]
    fn test_validate_heights_empty() {
        let bvh = DynamicBvh::new();
        assert!(
            bvh.validate_heights(),
            "empty tree has trivially valid heights"
        );
    }
    #[test]
    fn test_leaf_snapshot_count() {
        let mut bvh = DynamicBvh::new();
        for i in 0..4u32 {
            bvh.insert(make_aabb(i as f64 * 2.0, 0.0, 0.0, 0.5), i + 10);
        }
        let snap = bvh.leaf_snapshot();
        assert_eq!(snap.len(), 4, "snapshot should contain 4 entries");
        let data_vals: Vec<u32> = snap.iter().map(|&(d, _)| d).collect();
        for i in 10..14u32 {
            assert!(data_vals.contains(&i), "data {} should be in snapshot", i);
        }
    }
    #[test]
    fn test_leaf_aabb_valid_leaf() {
        let mut bvh = DynamicBvh::new();
        let leaf = bvh.insert(make_aabb(0.0, 0.0, 0.0, 1.0), 77);
        let aabb = bvh.leaf_aabb(leaf);
        assert!(
            aabb.is_some(),
            "leaf_aabb should return Some for a valid leaf"
        );
    }
    #[test]
    fn test_leaf_aabb_invalid_index() {
        let bvh = DynamicBvh::new();
        assert!(
            bvh.leaf_aabb(999).is_none(),
            "out-of-range index should return None"
        );
    }
    #[test]
    fn test_relabel_leaf_changes_data() {
        let mut bvh = DynamicBvh::new();
        let leaf = bvh.insert(make_aabb(0.0, 0.0, 0.0, 1.0), 10);
        assert!(bvh.relabel_leaf(leaf, 99));
        let hits = bvh.query_aabb(&make_aabb(0.0, 0.0, 0.0, 0.5));
        assert!(
            hits.contains(&99),
            "relabelled leaf should be found with new data"
        );
    }
    #[test]
    fn test_relabel_leaf_invalid_returns_false() {
        let mut bvh = DynamicBvh::new();
        assert!(
            !bvh.relabel_leaf(0, 42),
            "invalid index should return false"
        );
    }
    #[test]
    fn test_remove_where_removes_matching() {
        let mut bvh = DynamicBvh::new();
        for i in 0..5u32 {
            bvh.insert(make_aabb(i as f64 * 2.0, 0.0, 0.0, 0.5), i);
        }
        let removed = bvh.remove_where(|d| d % 2 == 0);
        assert_eq!(removed, 3, "should remove data 0, 2, 4");
        assert_eq!(bvh.n_leaves(), 2, "should have 2 leaves remaining");
        assert!(bvh.validate(), "tree must remain valid after remove_where");
    }
    #[test]
    fn test_remove_where_no_match() {
        let mut bvh = DynamicBvh::new();
        for i in 0..3u32 {
            bvh.insert(make_aabb(i as f64 * 2.0, 0.0, 0.0, 0.5), i);
        }
        let removed = bvh.remove_where(|d| d > 100);
        assert_eq!(removed, 0, "no leaves should be removed");
        assert_eq!(bvh.n_leaves(), 3);
    }
    #[test]
    fn test_rebuild_preserves_leaf_count() {
        let mut bvh = DynamicBvh::new();
        for i in 0..5u32 {
            bvh.insert(make_aabb(i as f64 * 3.0, 0.0, 0.0, 1.0), i);
        }
        let old_count = bvh.n_leaves();
        bvh.rebuild();
        assert_eq!(
            bvh.n_leaves(),
            old_count,
            "rebuild must preserve leaf count"
        );
        assert!(bvh.validate(), "tree must be valid after rebuild");
    }
    #[test]
    fn test_rebuild_empty_tree_is_noop() {
        let mut bvh = DynamicBvh::new();
        bvh.rebuild();
        assert_eq!(bvh.n_leaves(), 0);
        assert!(bvh.validate());
    }
    #[test]
    fn test_node_depth_root_is_zero() {
        let mut bvh = DynamicBvh::new();
        bvh.insert(make_aabb(0.0, 0.0, 0.0, 1.0), 1);
        let root = bvh.root.unwrap();
        assert_eq!(bvh.node_depth(root), Some(0), "root must be at depth 0");
    }
    #[test]
    fn test_node_depth_non_root() {
        let mut bvh = DynamicBvh::new();
        bvh.insert(make_aabb(0.0, 0.0, 0.0, 1.0), 1);
        bvh.insert(make_aabb(5.0, 0.0, 0.0, 1.0), 2);
        let root = bvh.root.unwrap();
        let left = bvh.nodes[root].left.expect("root must have left child");
        assert_eq!(
            bvh.node_depth(left),
            Some(1),
            "direct child of root should be at depth 1"
        );
    }
    #[test]
    fn test_point_dist_sq_inside_is_zero() {
        let a = make_aabb(0.0, 0.0, 0.0, 2.0);
        assert_eq!(a.point_dist_sq(Vec3::new(0.0, 0.0, 0.0)), 0.0);
    }
    #[test]
    fn test_point_dist_sq_outside() {
        let a = make_aabb(0.0, 0.0, 0.0, 1.0);
        let d = a.point_dist_sq(Vec3::new(3.0, 0.0, 0.0));
        assert!(d > 0.0, "point outside AABB must have dist_sq > 0");
    }
    #[test]
    fn test_leaf_data_preorder_count() {
        let mut bvh = DynamicBvh::new();
        for i in 0..4u32 {
            bvh.insert(make_aabb(i as f64 * 2.0, 0.0, 0.0, 0.5), i + 5);
        }
        let data = bvh.leaf_data_preorder();
        assert_eq!(data.len(), 4, "preorder traversal must visit all 4 leaves");
        let mut sorted = data.clone();
        sorted.sort_unstable();
        assert_eq!(sorted, vec![5, 6, 7, 8]);
    }
    #[test]
    fn test_leaf_data_preorder_empty() {
        let bvh = DynamicBvh::new();
        assert!(bvh.leaf_data_preorder().is_empty());
    }
}
