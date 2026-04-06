//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::math::Vec3;

use super::types::{GridSpatialIndex, KdNode, KdTree, Octree, OctreeNode, SpatialIndexStats};

/// Compute statistics for an `Octree`.
pub fn octree_stats<T: Clone>(tree: &Octree<T>) -> SpatialIndexStats {
    let mut stats = SpatialIndexStats {
        total_items: tree.len(),
        ..Default::default()
    };
    octree_node_stats(&tree.root, 0, &mut stats);
    if stats.n_leaves > 0 {
        stats.avg_items_per_leaf = stats.total_items as f64 / stats.n_leaves as f64;
    }
    stats
}
pub(super) fn octree_node_stats<T: Clone>(
    node: &OctreeNode<T>,
    depth: usize,
    stats: &mut SpatialIndexStats,
) {
    if depth > stats.max_depth {
        stats.max_depth = depth;
    }
    match node {
        OctreeNode::Leaf { items } => {
            stats.n_leaves += 1;
            let _ = items;
        }
        OctreeNode::Internal { children, .. } => {
            stats.n_internal += 1;
            for child in children.iter() {
                octree_node_stats(child, depth + 1, stats);
            }
        }
    }
}
/// Compute statistics for a `KdTree`.
pub fn kdtree_stats(tree: &KdTree) -> SpatialIndexStats {
    SpatialIndexStats {
        total_items: tree.len(),
        max_depth: kdtree_max_depth(
            &tree.nodes,
            if tree.nodes.is_empty() { None } else { Some(0) },
            0,
        ),
        n_leaves: tree
            .nodes
            .iter()
            .filter(|n| n.left.is_none() && n.right.is_none())
            .count(),
        n_internal: tree
            .nodes
            .iter()
            .filter(|n| n.left.is_some() || n.right.is_some())
            .count(),
        avg_items_per_leaf: 1.0,
    }
}
pub(super) fn kdtree_max_depth(nodes: &[KdNode], idx: Option<usize>, depth: usize) -> usize {
    match idx {
        None => depth,
        Some(i) => {
            let l = kdtree_max_depth(nodes, nodes[i].left, depth + 1);
            let r = kdtree_max_depth(nodes, nodes[i].right, depth + 1);
            l.max(r)
        }
    }
}
/// Spatial join: find all pairs (i, j) from two point sets within distance `r`.
///
/// Returns unique pairs `(i, j)` where point `i` is from set A and `j` from set B.
#[allow(dead_code)]
pub fn spatial_join(points_a: &[Vec3], points_b: &[Vec3], r: f64) -> Vec<(usize, usize)> {
    if points_a.is_empty() || points_b.is_empty() {
        return Vec::new();
    }
    let mut grid = GridSpatialIndex::new(r.max(1e-15));
    for &p in points_b {
        grid.insert(p);
    }
    let mut pairs = Vec::new();
    for (i, &pa) in points_a.iter().enumerate() {
        let neighbors = grid.range_query(pa, r);
        for j in neighbors {
            pairs.push((i, j));
        }
    }
    pairs
}
/// Spatial self-join: find all pairs (i, j) with i < j within distance `r`.
#[allow(dead_code)]
pub fn spatial_self_join(points: &[Vec3], r: f64) -> Vec<(usize, usize)> {
    let n = points.len();
    let mut grid = GridSpatialIndex::new(r.max(1e-15));
    for &p in points {
        grid.insert(p);
    }
    let r_sq = r * r;
    let mut pairs = Vec::new();
    for i in 0..n {
        let neighbors = grid.range_query(points[i], r);
        for j in neighbors {
            if j > i {
                let diff = points[i] - points[j];
                if diff.norm_squared() <= r_sq {
                    pairs.push((i, j));
                }
            }
        }
    }
    pairs
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::Vec3;

    use crate::spatial::FlatRTree;
    use crate::spatial::KdTree;
    use crate::spatial::KdTreeWithDeletion;
    use crate::spatial::LshIndex;

    use crate::spatial::RTree;
    use crate::spatial::RangeTree1D;
    use crate::spatial::SpatialAabb;

    fn unit_octree<T: Clone>(max_items: usize) -> Octree<T> {
        Octree::new(
            SpatialAabb::new(Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 1.0, 1.0)),
            8,
            max_items,
        )
    }
    #[test]
    fn test_octree_insert_query() {
        let mut tree = unit_octree::<u32>(16);
        for i in 0..100u32 {
            let t = i as f64 / 100.0;
            let p = Vec3::new(t, t * t, (1.0 - t) * t);
            tree.insert(p, i);
        }
        assert_eq!(tree.len(), 100);
        let hits = tree.query_sphere(Vec3::new(0.5, 0.25, 0.25), 0.1);
        assert!(!hits.is_empty());
        assert!(hits.len() < 100);
        for (p, _) in &hits {
            let d = (p - Vec3::new(0.5, 0.25, 0.25)).norm();
            assert!(d <= 0.1 + 1e-10, "point outside sphere: d={d}");
        }
    }
    #[test]
    fn test_octree_nearest_neighbor() {
        let mut tree = unit_octree::<usize>(8);
        let pts = [
            Vec3::new(0.1, 0.1, 0.1),
            Vec3::new(0.9, 0.9, 0.9),
            Vec3::new(0.5, 0.5, 0.5),
            Vec3::new(0.2, 0.8, 0.3),
            Vec3::new(0.7, 0.2, 0.6),
            Vec3::new(0.4, 0.6, 0.1),
            Vec3::new(0.3, 0.3, 0.7),
            Vec3::new(0.8, 0.4, 0.2),
            Vec3::new(0.1, 0.7, 0.9),
            Vec3::new(0.6, 0.1, 0.4),
        ];
        for (i, &p) in pts.iter().enumerate() {
            tree.insert(p, i);
        }
        let query = Vec3::new(0.45, 0.45, 0.45);
        let (nearest_pos, _val, dist) = tree.nearest_neighbor(query).unwrap();
        let (bf_pos, bf_dist) = pts
            .iter()
            .map(|&p| (p, (p - query).norm()))
            .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
            .unwrap();
        assert!(
            (dist - bf_dist).abs() < 1e-10,
            "dist mismatch: {dist} vs {bf_dist}"
        );
        assert!((nearest_pos - bf_pos).norm() < 1e-10, "position mismatch");
    }
    #[test]
    fn test_octree_subdivision() {
        let mut tree: Octree<i32> = Octree::new(
            SpatialAabb::new(Vec3::zeros(), Vec3::new(1.0, 1.0, 1.0)),
            8,
            2,
        );
        for i in 0..20i32 {
            let f = i as f64 / 20.0;
            tree.insert(Vec3::new(f, f, f), i);
        }
        assert_eq!(tree.len(), 20);
        let all = tree.query_aabb(&SpatialAabb::new(Vec3::zeros(), Vec3::new(1.0, 1.0, 1.0)));
        assert_eq!(all.len(), 20);
    }
    #[test]
    fn test_kdtree_build_nearest() {
        let pts = vec![
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(2.0, 2.0, 2.0),
            Vec3::new(0.1, 0.1, 0.1),
        ];
        let kd = KdTree::build(pts.clone());
        assert_eq!(kd.len(), 5);
        let query = Vec3::new(0.05, 0.05, 0.05);
        let (idx, d2) = kd.nearest(query).unwrap();
        let expected_idx = pts
            .iter()
            .enumerate()
            .min_by(|a, b| {
                (a.1 - query)
                    .norm_squared()
                    .partial_cmp(&(b.1 - query).norm_squared())
                    .unwrap()
            })
            .map(|(i, _)| i)
            .unwrap();
        assert_eq!(idx, expected_idx, "nearest index mismatch");
        let expected_d2 = (pts[expected_idx] - query).norm_squared();
        assert!((d2 - expected_d2).abs() < 1e-10, "d² mismatch");
    }
    #[test]
    fn test_kdtree_k_nearest() {
        let pts = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(2.0, 0.0, 0.0),
            Vec3::new(3.0, 0.0, 0.0),
            Vec3::new(4.0, 0.0, 0.0),
        ];
        let kd = KdTree::build(pts.clone());
        let query = Vec3::new(1.5, 0.0, 0.0);
        let k3 = kd.k_nearest(query, 3);
        assert_eq!(k3.len(), 3, "expected 3 results");
        for w in k3.windows(2) {
            assert!(w[0].1 <= w[1].1 + 1e-10, "not sorted");
        }
        let first_d2 = k3[0].1;
        assert!(
            first_d2 < 0.26,
            "closest should be within 0.5: d²={first_d2}"
        );
    }
    #[test]
    fn test_kdtree_range() {
        let pts = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.3, 0.0, 0.0),
            Vec3::new(0.6, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(2.0, 0.0, 0.0),
        ];
        let kd = KdTree::build(pts.clone());
        let query = Vec3::new(0.5, 0.0, 0.0);
        let mut within = kd.range_search(query, 0.5);
        within.sort_unstable();
        assert_eq!(
            within.len(),
            4,
            "expected 4 points within radius 0.5, got {}",
            within.len()
        );
        assert!(
            !within.contains(&4),
            "point at distance 1.5 should not be in range"
        );
    }
    #[test]
    fn test_spatial_aabb_operations() {
        let a = SpatialAabb::new(Vec3::zeros(), Vec3::new(2.0, 2.0, 2.0));
        let c = a.center();
        assert!((c - Vec3::new(1.0, 1.0, 1.0)).norm() < 1e-10);
        let he = a.half_extents();
        assert!((he - Vec3::new(1.0, 1.0, 1.0)).norm() < 1e-10);
        assert!(a.contains_point(Vec3::new(1.0, 1.0, 1.0)));
        assert!(a.contains_point(Vec3::new(0.0, 0.0, 0.0)));
        assert!(!a.contains_point(Vec3::new(3.0, 1.0, 1.0)));
        let b = SpatialAabb::new(Vec3::new(1.0, 1.0, 1.0), Vec3::new(3.0, 3.0, 3.0));
        assert!(a.intersects(&b));
        let c_box = SpatialAabb::new(Vec3::new(3.0, 3.0, 3.0), Vec3::new(5.0, 5.0, 5.0));
        assert!(!a.intersects(&c_box));
        let expanded = a.expand(1.0);
        assert!((expanded.min - Vec3::new(-1.0, -1.0, -1.0)).norm() < 1e-10);
        assert!((expanded.max - Vec3::new(3.0, 3.0, 3.0)).norm() < 1e-10);
        assert!(expanded.contains_point(Vec3::new(2.0, 2.0, 2.0)));
        assert!(expanded.contains_point(Vec3::new(-0.5, -0.5, -0.5)));
    }
    fn pt_aabb(x: f64, y: f64, z: f64) -> SpatialAabb {
        SpatialAabb::new(
            Vec3::new(x - 0.05, y - 0.05, z - 0.05),
            Vec3::new(x + 0.05, y + 0.05, z + 0.05),
        )
    }
    #[test]
    fn test_rtree_build_and_len() {
        let items: Vec<(SpatialAabb, u32)> = (0..20u32)
            .map(|i| {
                let f = i as f64 / 20.0;
                (pt_aabb(f, f, f), i)
            })
            .collect();
        let tree = RTree::build(items, 4);
        assert_eq!(tree.len(), 20);
    }
    #[test]
    fn test_rtree_query_finds_overlapping() {
        let items: Vec<(SpatialAabb, u32)> = vec![
            (pt_aabb(0.1, 0.1, 0.1), 0),
            (pt_aabb(0.5, 0.5, 0.5), 1),
            (pt_aabb(0.9, 0.9, 0.9), 2),
        ];
        let tree = RTree::build(items, 4);
        let query = SpatialAabb::new(Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.3, 0.3, 0.3));
        let results = tree.query(&query);
        assert!(!results.is_empty(), "Should find item near (0.1,0.1,0.1)");
        assert!(results.contains(&&0u32), "Should find item 0");
        assert!(!results.contains(&&2u32), "Should NOT find item 2");
    }
    #[test]
    fn test_flat_rtree_empty_original() {
        let tree: RTree<u32> = RTree::build(vec![], 4);
        assert!(tree.is_empty());
        let query = SpatialAabb::new(Vec3::zeros(), Vec3::new(1.0, 1.0, 1.0));
        assert!(tree.query(&query).is_empty());
    }
    #[test]
    fn test_rtree_query_no_overlap() {
        let items = vec![(pt_aabb(0.1, 0.1, 0.1), 0u32)];
        let tree = RTree::build(items, 4);
        let query = SpatialAabb::new(Vec3::new(5.0, 5.0, 5.0), Vec3::new(6.0, 6.0, 6.0));
        assert!(tree.query(&query).is_empty());
    }
    #[test]
    fn test_octree_stats_empty() {
        let tree: Octree<u32> = Octree::new(
            SpatialAabb::new(Vec3::zeros(), Vec3::new(1.0, 1.0, 1.0)),
            4,
            8,
        );
        let stats = octree_stats(&tree);
        assert_eq!(stats.total_items, 0);
        assert_eq!(stats.n_leaves, 1);
        assert_eq!(stats.n_internal, 0);
    }
    #[test]
    fn test_octree_stats_after_inserts() {
        let mut tree: Octree<u32> = Octree::new(
            SpatialAabb::new(Vec3::zeros(), Vec3::new(1.0, 1.0, 1.0)),
            6,
            2,
        );
        for i in 0..12u32 {
            let f = (i as f64 + 0.5) / 12.0;
            tree.insert(Vec3::new(f, f, f), i);
        }
        let stats = octree_stats(&tree);
        assert_eq!(stats.total_items, 12);
        assert!(stats.n_leaves >= 1);
    }
    #[test]
    fn test_kdtree_stats() {
        let pts = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(1.0, 1.0, 1.0),
        ];
        let kd = KdTree::build(pts);
        let stats = kdtree_stats(&kd);
        assert_eq!(stats.total_items, 5);
        assert!(
            stats.max_depth >= 2,
            "Tree of 5 nodes should have depth >= 2"
        );
        assert!(stats.n_leaves >= 1);
    }
    #[test]
    fn test_kdtree_deletion_n_active() {
        let pts = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(2.0, 0.0, 0.0),
        ];
        let mut kd = KdTreeWithDeletion::build(pts);
        assert_eq!(kd.n_active(), 3);
        kd.delete(1);
        assert_eq!(kd.n_active(), 2);
    }
    #[test]
    fn test_kdtree_deletion_nearest_active_skips_deleted() {
        let pts = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.1, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
        ];
        let mut kd = KdTreeWithDeletion::build(pts);
        let query = Vec3::new(0.05, 0.0, 0.0);
        let before = kd.nearest_active(query).map(|(i, _)| i);
        assert_eq!(before, Some(1), "Nearest should be index 1");
        kd.delete(1);
        let after = kd.nearest_active(query).map(|(i, _)| i);
        assert!(after != Some(1), "Deleted index 1 should not be returned");
    }
    #[test]
    fn test_kdtree_deletion_all_deleted_returns_none() {
        let pts = vec![Vec3::new(0.0, 0.0, 0.0)];
        let mut kd = KdTreeWithDeletion::build(pts);
        kd.delete(0);
        let result = kd.nearest_active(Vec3::zeros());
        assert!(result.is_none(), "All deleted: should return None");
    }
    #[test]
    fn test_range_tree_query_basic() {
        let values = vec![1.0, 3.0, 5.0, 7.0, 9.0];
        let tree = RangeTree1D::build(&values);
        let result = tree.range_query(3.0, 7.0);
        let mut sorted = result.clone();
        sorted.sort_unstable();
        assert_eq!(
            sorted,
            vec![1, 2, 3],
            "Range [3,7] should return indices 1,2,3"
        );
    }
    #[test]
    fn test_range_tree_empty_range() {
        let values = vec![1.0, 3.0, 5.0];
        let tree = RangeTree1D::build(&values);
        let result = tree.range_query(10.0, 20.0);
        assert!(result.is_empty(), "No values in [10,20]");
    }
    #[test]
    fn test_range_tree_all_in_range() {
        let values = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let tree = RangeTree1D::build(&values);
        let result = tree.range_query(0.0, 10.0);
        assert_eq!(result.len(), 5, "All 5 values should be in [0,10]");
    }
    #[test]
    fn test_range_tree_min_max() {
        let values = vec![5.0, 1.0, 9.0, 3.0];
        let tree = RangeTree1D::build(&values);
        assert!((tree.min_value().unwrap() - 1.0).abs() < 1e-12);
        assert!((tree.max_value().unwrap() - 9.0).abs() < 1e-12);
    }
    #[test]
    fn test_range_tree_empty() {
        let tree = RangeTree1D::build(&[]);
        assert!(tree.is_empty());
        assert!(tree.min_value().is_none());
        assert!(tree.max_value().is_none());
    }
    #[test]
    fn test_range_tree_single() {
        let tree = RangeTree1D::build(&[42.0]);
        let result = tree.range_query(42.0, 42.0);
        assert_eq!(result, vec![0]);
    }
    #[test]
    fn test_flat_rtree_multiple_overlaps() {
        let mut tree: FlatRTree<u32> = FlatRTree::new(4);
        for i in 0..5_u32 {
            let x = i as f64;
            tree.insert(
                SpatialAabb::new(Vec3::new(x, 0.0, 0.0), Vec3::new(x + 0.5, 0.5, 0.5)),
                i,
            );
        }
        let big_query = SpatialAabb::new(Vec3::new(-1.0, -1.0, -1.0), Vec3::new(10.0, 10.0, 10.0));
        let results = tree.query_overlap(&big_query);
        assert_eq!(results.len(), 5, "Large query should find all 5 items");
    }
    #[test]
    fn test_flat_rtree_nearest_in_new() {
        let mut tree: FlatRTree<usize> = FlatRTree::new(4);
        for i in 0..5_usize {
            let x = i as f64;
            tree.insert(
                SpatialAabb::new(Vec3::new(x, 0.0, 0.0), Vec3::new(x + 0.1, 0.1, 0.1)),
                i,
            );
        }
        let nn = tree.nearest(Vec3::new(2.05, 0.05, 0.05)).unwrap();
        assert_eq!(*nn, 2, "Nearest should be index 2");
    }
    #[test]
    fn test_grid_index_insert_query() {
        let mut grid = GridSpatialIndex::new(1.0);
        for i in 0..5_usize {
            grid.insert(Vec3::new(i as f64, 0.0, 0.0));
        }
        let result = grid.range_query(Vec3::new(2.0, 0.0, 0.0), 1.5);
        assert!(result.contains(&1));
        assert!(result.contains(&2));
        assert!(result.contains(&3));
        assert!(!result.contains(&0), "x=0 is 2.0 away from query");
        assert!(!result.contains(&4), "x=4 is 2.0 away from query");
    }
    #[test]
    fn test_grid_index_empty() {
        let grid = GridSpatialIndex::new(1.0);
        let result = grid.range_query(Vec3::zeros(), 10.0);
        assert!(result.is_empty());
    }
    #[test]
    fn test_grid_index_self_query() {
        let mut grid = GridSpatialIndex::new(1.0);
        let idx = grid.insert(Vec3::new(0.5, 0.5, 0.5));
        let result = grid.range_query(Vec3::new(0.5, 0.5, 0.5), 0.1);
        assert!(result.contains(&idx), "Should find itself");
    }
    #[test]
    fn test_spatial_join_basic() {
        let a = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(10.0, 0.0, 0.0)];
        let b = vec![Vec3::new(0.5, 0.0, 0.0), Vec3::new(10.5, 0.0, 0.0)];
        let pairs = spatial_join(&a, &b, 1.0);
        assert!(pairs.contains(&(0, 0)), "Close pair (0,0) should be found");
        assert!(pairs.contains(&(1, 1)), "Close pair (1,1) should be found");
        assert!(
            !pairs.contains(&(0, 1)),
            "Far pair (0,1) should not be found"
        );
    }
    #[test]
    fn test_spatial_join_empty() {
        let pairs = spatial_join(&[], &[Vec3::zeros()], 1.0);
        assert!(pairs.is_empty());
    }
    #[test]
    fn test_spatial_self_join() {
        let points = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.5, 0.0, 0.0),
            Vec3::new(5.0, 0.0, 0.0),
        ];
        let pairs = spatial_self_join(&points, 1.0);
        assert!(pairs.contains(&(0, 1)), "Close pair should be found");
        let far_pairs: Vec<_> = pairs.iter().filter(|&&(i, j)| i == 0 && j == 2).collect();
        assert!(far_pairs.is_empty(), "Far pair should not be found");
    }
    #[test]
    fn test_lsh_insert_and_nn() {
        let projections = vec![[1.0, 0.0, 0.0_f64], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let mut lsh = LshIndex::new(3, 0.5, &projections);
        for i in 0..10_usize {
            lsh.insert(Vec3::new(i as f64, 0.0, 0.0));
        }
        assert_eq!(lsh.len(), 10);
        let nn = lsh.query_approx_nn(Vec3::new(7.1, 0.0, 0.0)).unwrap();
        assert_eq!(
            lsh.points[nn].x.round() as usize,
            7,
            "NN should be at x=7, got x={}",
            lsh.points[nn].x
        );
    }
    #[test]
    fn test_lsh_empty() {
        let lsh = LshIndex::new(2, 1.0, &[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]);
        assert!(lsh.is_empty());
        assert!(lsh.query_approx_nn(Vec3::zeros()).is_none());
    }
    #[test]
    fn test_lsh_single_point() {
        let projections = vec![[1.0, 0.0, 0.0_f64]];
        let mut lsh = LshIndex::new(1, 1.0, &projections);
        lsh.insert(Vec3::new(3.0, 4.0, 5.0));
        let nn = lsh.query_approx_nn(Vec3::new(3.1, 4.0, 5.0)).unwrap();
        assert_eq!(nn, 0);
    }
}
/// Compute the Morton code (Z-order curve index) for a 3D integer coordinate.
///
/// Interleaves the bits of (x, y, z) to form a single 64-bit Morton code.
/// Each coordinate must be < 2^21 for 63-bit Morton codes.
#[allow(dead_code)]
pub fn morton_encode(x: u32, y: u32, z: u32) -> u64 {
    let x = expand_bits(x as u64);
    let y = expand_bits(y as u64);
    let z = expand_bits(z as u64);
    x | (y << 1) | (z << 2)
}
/// Expand bits by inserting two zero bits between each pair: abcd → a00b00c00d.
pub(super) fn expand_bits(mut v: u64) -> u64 {
    v &= 0x00000000001fffff;
    v = (v | (v << 32)) & 0x001f00000000ffff;
    v = (v | (v << 16)) & 0x001f0000ff0000ff;
    v = (v | (v << 8)) & 0x100f00f00f00f00f;
    v = (v | (v << 4)) & 0x10c30c30c30c30c3;
    v = (v | (v << 2)) & 0x1249249249249249;
    v
}
/// Decode a Morton code back to (x, y, z) integer coordinates.
#[allow(dead_code)]
pub fn morton_decode(code: u64) -> (u32, u32, u32) {
    let x = compact_bits(code);
    let y = compact_bits(code >> 1);
    let z = compact_bits(code >> 2);
    (x as u32, y as u32, z as u32)
}
/// Compact bits: reverse of `expand_bits`.
pub(super) fn compact_bits(mut v: u64) -> u64 {
    v &= 0x1249249249249249;
    v = (v | (v >> 2)) & 0x10c30c30c30c30c3;
    v = (v | (v >> 4)) & 0x100f00f00f00f00f;
    v = (v | (v >> 8)) & 0x001f0000ff0000ff;
    v = (v | (v >> 16)) & 0x001f00000000ffff;
    v = (v | (v >> 32)) & 0x00000000001fffff;
    v
}
/// Batch k-nearest-neighbor query.
///
/// For each query point, returns the k nearest points from the tree.
/// More efficient than calling `k_nearest` repeatedly if many queries share
/// the same tree.
#[allow(dead_code)]
pub fn kd_batch_knn(tree: &KdTree, queries: &[Vec3], k: usize) -> Vec<Vec<(usize, f64)>> {
    queries.iter().map(|&q| tree.k_nearest(q, k)).collect()
}
/// Find all pairs of points within distance `r` between two KD-trees.
///
/// Returns `(index_from_tree_a, index_from_tree_b)` pairs.
#[allow(dead_code)]
pub fn kd_cross_match(tree_a: &KdTree, tree_b: &KdTree, r: f64) -> Vec<(usize, usize)> {
    let mut pairs = Vec::new();
    for (i, &qa) in tree_a.points.iter().enumerate() {
        let near_b = tree_b.range_search(qa, r);
        for j in near_b {
            pairs.push((i, j));
        }
    }
    pairs
}
/// Compute the k-nearest neighbors of `query` from `points` using brute force.
///
/// Useful for verification and small datasets.
#[allow(dead_code)]
pub fn brute_force_knn(points: &[Vec3], query: Vec3, k: usize) -> Vec<(usize, f64)> {
    let mut dists: Vec<(usize, f64)> = points
        .iter()
        .enumerate()
        .map(|(i, &p)| (i, (p - query).norm()))
        .collect();
    dists.sort_unstable_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
    dists.truncate(k);
    dists
}
/// Compute the k-nearest neighbors of each point in `queries` from `points`.
#[allow(dead_code)]
pub fn brute_force_batch_knn(
    points: &[Vec3],
    queries: &[Vec3],
    k: usize,
) -> Vec<Vec<(usize, f64)>> {
    queries
        .iter()
        .map(|&q| brute_force_knn(points, q, k))
        .collect()
}
/// Build a uniform-spacing grid of points in 3D.
///
/// Returns points on a `(nx × ny × nz)` grid with spacing `(dx, dy, dz)` starting at `origin`.
#[allow(dead_code)]
pub fn regular_grid_3d(
    origin: Vec3,
    nx: usize,
    ny: usize,
    nz: usize,
    dx: f64,
    dy: f64,
    dz: f64,
) -> Vec<Vec3> {
    let mut pts = Vec::with_capacity(nx * ny * nz);
    for ix in 0..nx {
        for iy in 0..ny {
            for iz in 0..nz {
                pts.push(Vec3::new(
                    origin.x + ix as f64 * dx,
                    origin.y + iy as f64 * dy,
                    origin.z + iz as f64 * dz,
                ));
            }
        }
    }
    pts
}
#[cfg(test)]
mod tests_new_spatial {

    use crate::Vec3;
    use crate::spatial::BallTree;

    use crate::spatial::KdTree;

    use crate::spatial::MortonSortedIndex;

    use crate::spatial::VoxelGrid;
    use crate::spatial::brute_force_batch_knn;
    use crate::spatial::brute_force_knn;
    use crate::spatial::kd_batch_knn;
    use crate::spatial::kd_cross_match;
    use crate::spatial::morton_decode;
    use crate::spatial::morton_encode;
    use crate::spatial::regular_grid_3d;
    #[test]
    fn test_ball_tree_nearest_single_point() {
        let pts = vec![Vec3::new(1.0, 2.0, 3.0)];
        let tree = BallTree::build(pts);
        let nn = tree.nearest(Vec3::new(0.0, 0.0, 0.0));
        assert!(nn.is_some());
        assert_eq!(nn.unwrap().0, 0);
    }
    #[test]
    fn test_ball_tree_nearest_empty() {
        let tree = BallTree::build(vec![]);
        assert!(tree.nearest(Vec3::zeros()).is_none());
    }
    #[test]
    fn test_ball_tree_nearest_multiple_points() {
        let pts: Vec<Vec3> = (0..10).map(|i| Vec3::new(i as f64, 0.0, 0.0)).collect();
        let tree = BallTree::build(pts);
        let nn = tree.nearest(Vec3::new(7.1, 0.0, 0.0));
        assert!(nn.is_some());
        let (idx, _) = nn.unwrap();
        assert_eq!(idx, 7, "nearest to 7.1 should be at x=7, got idx={idx}");
    }
    #[test]
    fn test_ball_tree_range_query() {
        let pts: Vec<Vec3> = (0..20).map(|i| Vec3::new(i as f64, 0.0, 0.0)).collect();
        let tree = BallTree::build(pts);
        let result = tree.range_query(Vec3::new(10.0, 0.0, 0.0), 2.5);
        assert!(result.contains(&8));
        assert!(result.contains(&10));
        assert!(result.contains(&12));
        assert!(!result.contains(&6), "x=6 is 4 away, should not be found");
    }
    #[test]
    fn test_ball_tree_range_empty() {
        let pts = vec![Vec3::new(100.0, 0.0, 0.0)];
        let tree = BallTree::build(pts);
        let result = tree.range_query(Vec3::zeros(), 1.0);
        assert!(result.is_empty(), "far point should not be in range");
    }
    #[test]
    fn test_ball_tree_len() {
        let pts: Vec<Vec3> = (0..5).map(|i| Vec3::new(i as f64, 0.0, 0.0)).collect();
        let tree = BallTree::build(pts);
        assert_eq!(tree.len(), 5);
        assert!(!tree.is_empty());
    }
    #[test]
    fn test_morton_encode_origin() {
        assert_eq!(morton_encode(0, 0, 0), 0);
    }
    #[test]
    fn test_morton_encode_decode_roundtrip() {
        for &(x, y, z) in &[(1u32, 2u32, 3u32), (7, 15, 0), (100, 200, 50)] {
            let code = morton_encode(x, y, z);
            let (dx, dy, dz) = morton_decode(code);
            assert_eq!(
                (dx, dy, dz),
                (x, y, z),
                "Morton roundtrip failed for ({x},{y},{z}): got ({dx},{dy},{dz})"
            );
        }
    }
    #[test]
    fn test_morton_encode_increasing() {
        let m0 = morton_encode(0, 0, 0);
        let m1 = morton_encode(1, 0, 0);
        let m2 = morton_encode(2, 0, 0);
        assert!(m1 > m0, "Morton(1,0,0) should be > Morton(0,0,0)");
        assert!(m2 > m1, "Morton(2,0,0) should be > Morton(1,0,0)");
    }
    #[test]
    fn test_morton_decode_identity() {
        let (x, y, z) = morton_decode(0);
        assert_eq!((x, y, z), (0, 0, 0));
    }
    #[test]
    fn test_morton_sorted_index_len() {
        let pts: Vec<Vec3> = (0..10).map(|i| Vec3::new(i as f64, 0.0, 0.0)).collect();
        let idx = MortonSortedIndex::build(&pts, 1024);
        assert_eq!(idx.len(), 10);
        assert!(!idx.is_empty());
    }
    #[test]
    fn test_morton_sorted_index_sorted_order() {
        let pts: Vec<Vec3> = (0..8).map(|i| Vec3::new(i as f64, 0.0, 0.0)).collect();
        let idx = MortonSortedIndex::build(&pts, 64);
        let codes: Vec<u64> = idx.sorted.iter().map(|&(c, _)| c).collect();
        for w in codes.windows(2) {
            assert!(
                w[0] <= w[1],
                "Morton codes should be sorted: {} <= {}",
                w[0],
                w[1]
            );
        }
    }
    #[test]
    fn test_morton_sorted_index_empty() {
        let idx = MortonSortedIndex::build(&[], 16);
        assert!(idx.is_empty());
    }
    #[test]
    fn test_morton_sorted_index_all_indices_present() {
        let pts: Vec<Vec3> = (0..5).map(|i| Vec3::new(i as f64, i as f64, 0.0)).collect();
        let idx = MortonSortedIndex::build(&pts, 64);
        let sorted = idx.sorted_indices();
        let mut check = sorted.clone();
        check.sort_unstable();
        assert_eq!(check, vec![0, 1, 2, 3, 4], "all indices should be present");
    }
    #[test]
    fn test_voxel_grid_insert_query() {
        let mut grid = VoxelGrid::new(1.0);
        for i in 0..5 {
            grid.insert(Vec3::new(i as f64, 0.0, 0.0));
        }
        let result = grid.range_query(Vec3::new(2.0, 0.0, 0.0), 1.5);
        assert!(result.contains(&1), "x=1 should be found");
        assert!(result.contains(&2), "x=2 should be found");
        assert!(result.contains(&3), "x=3 should be found");
        assert!(!result.contains(&0), "x=0 is 2 away, should not be found");
    }
    #[test]
    fn test_voxel_grid_nearest_linear() {
        let mut grid = VoxelGrid::new(1.0);
        for i in 0..5 {
            grid.insert(Vec3::new(i as f64, 0.0, 0.0));
        }
        let nn = grid.nearest_linear(Vec3::new(3.9, 0.0, 0.0)).unwrap();
        assert_eq!(nn.0, 4, "nearest to 3.9 should be at x=4, got {}", nn.0);
    }
    #[test]
    fn test_voxel_grid_empty() {
        let grid = VoxelGrid::new(1.0);
        assert!(grid.is_empty());
        let result = grid.range_query(Vec3::zeros(), 10.0);
        assert!(result.is_empty());
        assert!(grid.nearest_linear(Vec3::zeros()).is_none());
    }
    #[test]
    fn test_voxel_grid_centroid_downsampling() {
        let mut grid = VoxelGrid::new(1.0);
        grid.insert(Vec3::new(0.0, 0.0, 0.0));
        grid.insert(Vec3::new(0.5, 0.0, 0.0));
        let centroids = grid.voxel_centroids();
        assert_eq!(
            centroids.len(),
            1,
            "both points in same voxel should give 1 centroid"
        );
        assert!(
            (centroids[0].x - 0.25).abs() < 1e-12,
            "centroid should be at 0.25, got {}",
            centroids[0].x
        );
    }
    #[test]
    fn test_brute_force_knn_basic() {
        let pts: Vec<Vec3> = (0..10).map(|i| Vec3::new(i as f64, 0.0, 0.0)).collect();
        let result = brute_force_knn(&pts, Vec3::new(5.0, 0.0, 0.0), 3);
        assert_eq!(result.len(), 3);
        assert_eq!(result[0].0, 5, "nearest should be at x=5");
    }
    #[test]
    fn test_brute_force_knn_empty() {
        let result = brute_force_knn(&[], Vec3::zeros(), 5);
        assert!(result.is_empty());
    }
    #[test]
    fn test_brute_force_knn_sorted() {
        let pts: Vec<Vec3> = (0..10).map(|i| Vec3::new(i as f64, 0.0, 0.0)).collect();
        let result = brute_force_knn(&pts, Vec3::new(5.0, 0.0, 0.0), 5);
        for w in result.windows(2) {
            assert!(w[0].1 <= w[1].1, "KNN results should be sorted by distance");
        }
    }
    #[test]
    fn test_brute_force_batch_knn() {
        let pts: Vec<Vec3> = (0..5).map(|i| Vec3::new(i as f64, 0.0, 0.0)).collect();
        let queries = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(4.0, 0.0, 0.0)];
        let result = brute_force_batch_knn(&pts, &queries, 2);
        assert_eq!(result.len(), 2);
        assert_eq!(result[0][0].0, 0, "nearest to x=0 should be index 0");
        assert_eq!(result[1][0].0, 4, "nearest to x=4 should be index 4");
    }
    #[test]
    fn test_kd_batch_knn() {
        let pts: Vec<Vec3> = (0..10).map(|i| Vec3::new(i as f64, 0.0, 0.0)).collect();
        let tree = KdTree::build(pts);
        let queries = vec![Vec3::new(3.0, 0.0, 0.0), Vec3::new(8.0, 0.0, 0.0)];
        let results = kd_batch_knn(&tree, &queries, 2);
        assert_eq!(results.len(), 2);
        assert_eq!(results[0].len(), 2);
        assert_eq!(results[1].len(), 2);
    }
    #[test]
    fn test_kd_cross_match_close_pairs() {
        let pts_a: Vec<Vec3> = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(10.0, 0.0, 0.0)];
        let pts_b: Vec<Vec3> = vec![Vec3::new(0.3, 0.0, 0.0), Vec3::new(10.3, 0.0, 0.0)];
        let tree_a = KdTree::build(pts_a);
        let tree_b = KdTree::build(pts_b);
        let pairs = kd_cross_match(&tree_a, &tree_b, 0.5);
        assert!(pairs.contains(&(0, 0)), "close pair (0,0) should be found");
        assert!(pairs.contains(&(1, 1)), "close pair (1,1) should be found");
    }
    #[test]
    fn test_kd_cross_match_no_pairs() {
        let pts_a: Vec<Vec3> = vec![Vec3::new(0.0, 0.0, 0.0)];
        let pts_b: Vec<Vec3> = vec![Vec3::new(100.0, 0.0, 0.0)];
        let tree_a = KdTree::build(pts_a);
        let tree_b = KdTree::build(pts_b);
        let pairs = kd_cross_match(&tree_a, &tree_b, 1.0);
        assert!(pairs.is_empty(), "far-apart points should have no pairs");
    }
    #[test]
    fn test_regular_grid_3d_count() {
        let pts = regular_grid_3d(Vec3::zeros(), 3, 4, 5, 1.0, 1.0, 1.0);
        assert_eq!(pts.len(), 60, "3x4x5 grid should have 60 points");
    }
    #[test]
    fn test_regular_grid_3d_origin() {
        let pts = regular_grid_3d(Vec3::new(1.0, 2.0, 3.0), 2, 2, 2, 1.0, 1.0, 1.0);
        assert_eq!(pts.len(), 8);
        assert!((pts[0].x - 1.0).abs() < 1e-12);
        assert!((pts[0].y - 2.0).abs() < 1e-12);
        assert!((pts[0].z - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_regular_grid_3d_spacing() {
        let pts = regular_grid_3d(Vec3::zeros(), 3, 1, 1, 0.5, 1.0, 1.0);
        assert_eq!(pts.len(), 3);
        assert!(
            (pts[1].x - 0.5).abs() < 1e-12,
            "second point should be at x=0.5"
        );
        assert!(
            (pts[2].x - 1.0).abs() < 1e-12,
            "third point should be at x=1.0"
        );
    }
}
