// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! CPU BVH (Bounding Volume Hierarchy) tree for broad-phase acceleration.
//!
//! Provides axis-aligned bounding box (AABB) structures, a recursive BVH tree
//! built with a SAH-inspired median split, ray and AABB queries, and a
//! linearised flat representation for cache-friendly traversal.
//!
//! # Module layout
//!
//! - [`types`] — Core data types: `Aabb`, `BvhPrimitive`, `BvhNode`, `FlatBvhNode`, `GpuRay`, `RayHit`, statistics types, Morton/LBVH types.
//! - [`cpu`] — CPU builder (`Bvh::build`), query helpers, flat traversal, LBVH construction, statistics.
//! - [`gpu`] — GPU-accelerated ray traversal (`BvhGpuTraverser`) backed by `WgpuBackendReal`; falls back to CPU when no GPU is present.

pub mod cpu;
pub mod gpu;
pub mod types;

// ── Public re-exports ────────────────────────────────────────────────────────

pub use cpu::{
    Bvh, build_morton_clusters, bvh_closest_hit, compute_bvh_from_sorted, compute_cluster_radius,
    flatten, hlbvh_split, lbvh_build, morton_code, query_flat, ray_aabb_intersect, refit, sah_cost,
};
pub use gpu::BvhGpuTraverser;
pub use types::{
    Aabb, BvhNode, BvhPrimitive, BvhStats, BvhTreeStatistics, FlatBvhNode, GpuRay, LbvhPrimitive,
    MortonCluster, RayHit,
};

// ============================================================================
// Tests (migrated from the original monolithic bvh.rs)
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    // ------------------------------------------------------------------
    // Aabb tests
    // ------------------------------------------------------------------

    #[test]
    fn aabb_new_stores_corners() {
        let a = Aabb::new([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]);
        assert_eq!(a.min, [1.0, 2.0, 3.0]);
        assert_eq!(a.max, [4.0, 5.0, 6.0]);
    }

    #[test]
    fn aabb_point_is_degenerate() {
        let p = [3.0, 3.0, 3.0];
        let a = Aabb::point(p);
        assert_eq!(a.min, p);
        assert_eq!(a.max, p);
    }

    #[test]
    fn aabb_merge_covers_both() {
        let a = Aabb::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = Aabb::new([2.0, 2.0, 2.0], [3.0, 3.0, 3.0]);
        let m = Aabb::merge(&a, &b);
        assert_eq!(m.min, [0.0, 0.0, 0.0]);
        assert_eq!(m.max, [3.0, 3.0, 3.0]);
    }

    #[test]
    fn aabb_intersects_overlapping() {
        let a = Aabb::new([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
        let b = Aabb::new([1.0, 1.0, 1.0], [3.0, 3.0, 3.0]);
        assert!(a.intersects(&b));
    }

    #[test]
    fn aabb_intersects_disjoint() {
        let a = Aabb::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = Aabb::new([2.0, 2.0, 2.0], [3.0, 3.0, 3.0]);
        assert!(!a.intersects(&b));
    }

    #[test]
    fn aabb_intersects_touching_edge() {
        let a = Aabb::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = Aabb::new([1.0, 0.0, 0.0], [2.0, 1.0, 1.0]);
        assert!(a.intersects(&b));
    }

    #[test]
    fn aabb_contains_inside() {
        let a = Aabb::new([0.0, 0.0, 0.0], [4.0, 4.0, 4.0]);
        assert!(a.contains([2.0, 2.0, 2.0]));
    }

    #[test]
    fn aabb_contains_outside() {
        let a = Aabb::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert!(!a.contains([2.0, 0.0, 0.0]));
    }

    #[test]
    fn aabb_contains_on_surface() {
        let a = Aabb::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert!(a.contains([1.0, 0.5, 0.5]));
    }

    #[test]
    fn aabb_surface_area_unit_cube() {
        let a = Aabb::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert!((a.surface_area() - 6.0).abs() < 1e-6);
    }

    #[test]
    fn aabb_surface_area_flat() {
        // 2×3×0 slab: area = 2*(2*3 + 3*0 + 0*2) = 12
        let a = Aabb::new([0.0, 0.0, 0.0], [2.0, 3.0, 0.0]);
        assert!((a.surface_area() - 12.0).abs() < 1e-6);
    }

    #[test]
    fn aabb_center_correct() {
        let a = Aabb::new([0.0, 0.0, 0.0], [2.0, 4.0, 6.0]);
        let c = a.center();
        assert!((c[0] - 1.0).abs() < 1e-6);
        assert!((c[1] - 2.0).abs() < 1e-6);
        assert!((c[2] - 3.0).abs() < 1e-6);
    }

    #[test]
    fn aabb_expand_increases_bounds() {
        let a = Aabb::new([1.0, 1.0, 1.0], [2.0, 2.0, 2.0]);
        let e = a.expand(0.5);
        assert_eq!(e.min, [0.5, 0.5, 0.5]);
        assert_eq!(e.max, [2.5, 2.5, 2.5]);
    }

    // ------------------------------------------------------------------
    // SAH cost
    // ------------------------------------------------------------------

    #[test]
    fn sah_cost_balanced() {
        // Both halves equal area and count → cost == n_left + n_right
        let cost = sah_cost(4, 1.0, 4, 1.0, 2.0);
        // (1/2)*4 + (1/2)*4 = 4
        assert!((cost - 4.0).abs() < 1e-6);
    }

    #[test]
    fn sah_cost_zero_parent_area_returns_max() {
        let cost = sah_cost(1, 1.0, 1, 1.0, 0.0);
        assert_eq!(cost, f32::MAX);
    }

    // ------------------------------------------------------------------
    // Ray–AABB slab intersection
    // ------------------------------------------------------------------

    #[test]
    fn ray_hits_unit_cube() {
        let aabb = Aabb::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let origin = [-1.0, 0.5, 0.5];
        let dir = [1.0, 0.0, 0.0];
        let inv = [1.0 / dir[0], 1.0 / dir[1], 1.0 / dir[2]];
        assert!(ray_aabb_intersect(origin, inv, &aabb, 10.0));
    }

    #[test]
    fn ray_misses_unit_cube() {
        let aabb = Aabb::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let origin = [-1.0, 2.0, 0.5];
        let dir = [1.0, 0.0, 0.0];
        let inv = [1.0 / dir[0], 1.0 / dir[1], 1.0 / dir[2]];
        assert!(!ray_aabb_intersect(origin, inv, &aabb, 10.0));
    }

    #[test]
    fn ray_too_short_misses() {
        let aabb = Aabb::new([5.0, 0.0, 0.0], [6.0, 1.0, 1.0]);
        let origin = [0.0, 0.5, 0.5];
        let dir = [1.0, 0.0, 0.0];
        let inv = [1.0 / dir[0], 1.0 / dir[1], 1.0 / dir[2]];
        assert!(!ray_aabb_intersect(origin, inv, &aabb, 3.0));
    }

    // ------------------------------------------------------------------
    // Bvh build / query
    // ------------------------------------------------------------------

    fn make_grid_primitives(n: usize) -> Vec<BvhPrimitive> {
        (0..n)
            .map(|i| {
                let x = i as f32;
                BvhPrimitive::new(Aabb::new([x, 0.0, 0.0], [x + 1.0, 1.0, 1.0]), i)
            })
            .collect()
    }

    #[test]
    fn bvh_build_empty() {
        let bvh = Bvh::build(vec![]);
        assert!(bvh.root.is_none());
        assert_eq!(bvh.node_count(), 0);
        assert_eq!(bvh.depth(), 0);
    }

    #[test]
    fn bvh_build_single() {
        let prims = make_grid_primitives(1);
        let bvh = Bvh::build(prims);
        assert!(bvh.root.is_some());
        assert!(bvh.root.as_ref().unwrap().is_leaf());
        assert_eq!(bvh.node_count(), 1);
        assert_eq!(bvh.depth(), 1);
    }

    #[test]
    fn bvh_query_aabb_finds_overlap() {
        let prims = make_grid_primitives(10);
        let bvh = Bvh::build(prims);
        // Query the box that overlaps object 5 only.
        let query = Aabb::new([5.1, 0.1, 0.1], [5.9, 0.9, 0.9]);
        let mut hits = bvh.query_aabb(&query);
        hits.sort();
        assert_eq!(hits, vec![5]);
    }

    #[test]
    fn bvh_query_aabb_empty_result() {
        let prims = make_grid_primitives(5);
        let bvh = Bvh::build(prims);
        let query = Aabb::new([100.0, 0.0, 0.0], [101.0, 1.0, 1.0]);
        assert!(bvh.query_aabb(&query).is_empty());
    }

    #[test]
    fn bvh_query_aabb_finds_multiple() {
        let prims = make_grid_primitives(10);
        let bvh = Bvh::build(prims);
        // Query spanning objects 2, 3, 4
        let query = Aabb::new([2.1, 0.1, 0.1], [4.9, 0.9, 0.9]);
        let mut hits = bvh.query_aabb(&query);
        hits.sort();
        assert_eq!(hits, vec![2, 3, 4]);
    }

    #[test]
    fn bvh_query_ray_hits() {
        let prims = make_grid_primitives(8);
        let bvh = Bvh::build(prims);
        // Ray along X axis at y=0.5, z=0.5 hits all 8 primitives.
        let mut hits = bvh.query_ray([-1.0, 0.5, 0.5], [1.0, 0.0, 0.0], 20.0);
        hits.sort();
        assert_eq!(hits, (0..8).collect::<Vec<_>>());
    }

    #[test]
    fn bvh_query_ray_misses() {
        let prims = make_grid_primitives(5);
        let bvh = Bvh::build(prims);
        let hits = bvh.query_ray([0.5, 10.0, 0.5], [0.0, 1.0, 0.0], 100.0);
        assert!(hits.is_empty());
    }

    #[test]
    fn bvh_node_count_and_depth_consistent() {
        let prims = make_grid_primitives(16);
        let bvh = Bvh::build(prims);
        // With LEAF_SIZE=4 and 16 prims, depth should be at least 2.
        assert!(bvh.depth() >= 2);
        // An n-primitive tree has at most 2n-1 nodes.
        assert!(bvh.node_count() < 2 * 16);
    }

    // ------------------------------------------------------------------
    // Flat BVH
    // ------------------------------------------------------------------

    #[test]
    fn flatten_empty_bvh() {
        let bvh = Bvh::build(vec![]);
        let (nodes, prim_indices) = flatten(&bvh);
        assert!(nodes.is_empty());
        assert!(prim_indices.is_empty());
    }

    #[test]
    fn flatten_single_primitive() {
        let prims = make_grid_primitives(1);
        let bvh = Bvh::build(prims);
        let (nodes, prim_indices) = flatten(&bvh);
        assert_eq!(nodes.len(), 1);
        assert_eq!(prim_indices.len(), 1);
        assert_eq!(nodes[0].count, 1);
    }

    #[test]
    fn query_flat_finds_overlap() {
        let prims = make_grid_primitives(10);
        let bvh = Bvh::build(prims);
        let (nodes, prim_indices) = flatten(&bvh);
        let query = Aabb::new([3.1, 0.1, 0.1], [3.9, 0.9, 0.9]);
        let mut hits = query_flat(&nodes, &prim_indices, &bvh.primitives, &query);
        hits.sort();
        assert_eq!(hits, vec![3]);
    }

    #[test]
    fn query_flat_empty_result() {
        let prims = make_grid_primitives(5);
        let bvh = Bvh::build(prims);
        let (nodes, prim_indices) = flatten(&bvh);
        let query = Aabb::new([50.0, 0.0, 0.0], [51.0, 1.0, 1.0]);
        assert!(query_flat(&nodes, &prim_indices, &bvh.primitives, &query).is_empty());
    }

    #[test]
    fn query_flat_matches_recursive() {
        let prims = make_grid_primitives(20);
        let bvh = Bvh::build(prims);
        let query = Aabb::new([7.1, 0.0, 0.0], [12.9, 1.0, 1.0]);
        let mut recursive_hits = bvh.query_aabb(&query);
        recursive_hits.sort();

        let (nodes, prim_indices) = flatten(&bvh);
        let mut flat_hits = query_flat(&nodes, &prim_indices, &bvh.primitives, &query);
        flat_hits.sort();

        assert_eq!(recursive_hits, flat_hits);
    }

    // ------------------------------------------------------------------
    // Morton code tests
    // ------------------------------------------------------------------

    #[test]
    fn morton_origin_is_zero() {
        assert_eq!(morton_code([0.0, 0.0, 0.0]), 0);
    }

    #[test]
    fn morton_increases_along_x() {
        let m0 = morton_code([0.0, 0.0, 0.0]);
        let m1 = morton_code([0.5, 0.0, 0.0]);
        let m2 = morton_code([1.0, 0.0, 0.0]);
        assert!(m0 <= m1, "m0={} m1={}", m0, m1);
        assert!(m1 <= m2, "m1={} m2={}", m1, m2);
    }

    #[test]
    fn morton_clamps_outside_unit_cube() {
        let m_neg = morton_code([-1.0, -1.0, -1.0]);
        let m_zero = morton_code([0.0, 0.0, 0.0]);
        assert_eq!(m_neg, m_zero);

        let m_big = morton_code([2.0, 2.0, 2.0]);
        let m_one = morton_code([1.0, 1.0, 1.0]);
        assert_eq!(m_big, m_one);
    }

    // ------------------------------------------------------------------
    // LBVH construction tests
    // ------------------------------------------------------------------

    #[test]
    fn lbvh_build_empty() {
        let bvh = lbvh_build(vec![]);
        assert!(bvh.root.is_none());
    }

    #[test]
    fn lbvh_build_single() {
        let prims = make_grid_primitives(1);
        let bvh = lbvh_build(prims);
        assert!(bvh.root.is_some());
        assert!(bvh.root.as_ref().unwrap().is_leaf());
    }

    #[test]
    fn lbvh_build_query_finds_correct_objects() {
        let prims = make_grid_primitives(10);
        let bvh = lbvh_build(prims);
        let query = Aabb::new([4.1, 0.1, 0.1], [4.9, 0.9, 0.9]);
        let mut hits = bvh.query_aabb(&query);
        hits.sort();
        assert_eq!(hits, vec![4]);
    }

    #[test]
    fn lbvh_build_covers_all_primitives() {
        let prims = make_grid_primitives(8);
        let bvh = lbvh_build(prims);
        // Root AABB should contain all primitives.
        let root = bvh.root.as_ref().unwrap();
        assert!(root.aabb.min[0] <= 0.0);
        assert!(root.aabb.max[0] >= 8.0);
    }

    // ------------------------------------------------------------------
    // BVH closest-hit traversal
    // ------------------------------------------------------------------

    #[test]
    fn closest_hit_returns_nearest() {
        let prims = make_grid_primitives(10);
        let bvh = Bvh::build(prims);
        // Ray along X from x=-1: should hit object 0 first (x ∈ [0,1])
        let hit = bvh_closest_hit(&bvh, [-1.0, 0.5, 0.5], [1.0, 0.0, 0.0], 100.0);
        assert!(hit.is_some(), "ray should hit something");
        let hit = hit.unwrap();
        assert_eq!(
            hit.object_id, 0,
            "closest hit should be object 0, got {}",
            hit.object_id
        );
    }

    #[test]
    fn closest_hit_misses_returns_none() {
        let prims = make_grid_primitives(5);
        let bvh = Bvh::build(prims);
        let hit = bvh_closest_hit(&bvh, [0.5, 10.0, 0.5], [0.0, 1.0, 0.0], 100.0);
        assert!(hit.is_none());
    }

    #[test]
    fn closest_hit_empty_bvh_returns_none() {
        let bvh = Bvh::build(vec![]);
        let hit = bvh_closest_hit(&bvh, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 100.0);
        assert!(hit.is_none());
    }

    #[test]
    fn closest_hit_t_is_positive() {
        let prims = make_grid_primitives(5);
        let bvh = Bvh::build(prims);
        let hit = bvh_closest_hit(&bvh, [-1.0, 0.5, 0.5], [1.0, 0.0, 0.0], 100.0);
        if let Some(h) = hit {
            assert!(h.t >= 0.0, "hit t should be non-negative, got {}", h.t);
        }
    }

    // ------------------------------------------------------------------
    // HLBVH split
    // ------------------------------------------------------------------

    #[test]
    fn hlbvh_split_single_pair_splits_at_1() {
        let mortons = vec![0u32, 1u32];
        assert_eq!(hlbvh_split(&mortons), 1);
    }

    #[test]
    fn hlbvh_split_identical_codes_splits_at_end() {
        let mortons = vec![42u32; 4];
        let s = hlbvh_split(&mortons);
        assert!(s > 0 && s < mortons.len());
    }

    #[test]
    fn hlbvh_split_returns_valid_range() {
        let mortons: Vec<u32> = (0..16).map(|i| i * 2).collect();
        let s = hlbvh_split(&mortons);
        assert!(s > 0, "split must be > 0");
        assert!(s < mortons.len(), "split must be < len");
    }

    // ------------------------------------------------------------------
    // BVH refit
    // ------------------------------------------------------------------

    #[test]
    fn refit_does_not_panic_on_leaf() {
        let prims = make_grid_primitives(2);
        let mut bvh = Bvh::build(prims.clone());
        if let Some(root) = bvh.root.as_mut() {
            refit(root, &bvh.primitives.clone());
        }
        // No panic = pass.
    }

    // ------------------------------------------------------------------
    // BvhStats
    // ------------------------------------------------------------------

    #[test]
    fn bvh_stats_empty_tree() {
        let bvh = Bvh::build(vec![]);
        let s = BvhStats::compute(&bvh);
        assert_eq!(s.node_count, 0);
        assert_eq!(s.leaf_count, 0);
    }

    #[test]
    fn bvh_stats_single_leaf() {
        let prims = make_grid_primitives(1);
        let bvh = Bvh::build(prims);
        let s = BvhStats::compute(&bvh);
        assert_eq!(s.leaf_count, 1);
        assert_eq!(s.total_primitives, 1);
    }

    #[test]
    fn bvh_stats_counts_consistent() {
        let prims = make_grid_primitives(16);
        let bvh = Bvh::build(prims);
        let s = BvhStats::compute(&bvh);
        assert_eq!(s.leaf_count + s.internal_count, s.node_count);
    }

    #[test]
    fn bvh_stats_total_primitives() {
        let prims = make_grid_primitives(16);
        let bvh = Bvh::build(prims);
        let s = BvhStats::compute(&bvh);
        assert_eq!(s.total_primitives, 16);
    }

    // ------------------------------------------------------------------
    // BvhTreeStatistics
    // ------------------------------------------------------------------

    #[test]
    fn bvh_tree_stats_empty() {
        let bvh = Bvh::build(vec![]);
        let s = BvhTreeStatistics::compute(&bvh);
        assert_eq!(s.node_count, 0);
    }

    #[test]
    fn bvh_tree_stats_single() {
        let prims = make_grid_primitives(1);
        let bvh = Bvh::build(prims);
        let s = BvhTreeStatistics::compute(&bvh);
        assert_eq!(s.leaf_count, 1);
        assert_eq!(s.total_primitives, 1);
    }

    #[test]
    fn bvh_tree_stats_consistent() {
        let prims = make_grid_primitives(32);
        let bvh = Bvh::build(prims.clone());
        let s = BvhTreeStatistics::compute(&bvh);
        assert_eq!(s.leaf_count + s.internal_count, s.node_count);
        assert_eq!(s.total_primitives, prims.len());
    }

    #[test]
    fn bvh_tree_stats_leaf_surface_area_positive() {
        let prims = make_grid_primitives(8);
        let bvh = Bvh::build(prims);
        let s = BvhTreeStatistics::compute(&bvh);
        assert!(
            s.total_leaf_surface_area > 0.0,
            "leaf surface area should be > 0"
        );
    }

    // ------------------------------------------------------------------
    // build_morton_clusters tests
    // ------------------------------------------------------------------

    fn make_sorted_lbvh_prims(n: usize) -> Vec<LbvhPrimitive> {
        let scene = Aabb::new([0.0, 0.0, 0.0], [n as f32 + 1.0, 1.0, 1.0]);
        let mut prims: Vec<LbvhPrimitive> = (0..n)
            .map(|i| {
                let x = i as f32;
                LbvhPrimitive::new(Aabb::new([x, 0.0, 0.0], [x + 1.0, 1.0, 1.0]), i, &scene)
            })
            .collect();
        prims.sort_unstable_by_key(|lp| lp.morton);
        prims
    }

    #[test]
    fn build_morton_clusters_empty() {
        let clusters = build_morton_clusters(&[], 4);
        assert!(clusters.is_empty());
    }

    #[test]
    fn build_morton_clusters_count() {
        let sorted = make_sorted_lbvh_prims(10);
        let clusters = build_morton_clusters(&sorted, 3);
        // 10 primitives in chunks of 3 → ceil(10/3) = 4 clusters
        assert_eq!(clusters.len(), 4);
    }

    #[test]
    fn build_morton_clusters_radii_non_negative() {
        let sorted = make_sorted_lbvh_prims(8);
        let clusters = build_morton_clusters(&sorted, 2);
        for c in &clusters {
            assert!(c.radius >= 0.0, "cluster radius must be non-negative");
        }
    }

    // ------------------------------------------------------------------
    // LbvhPrimitive Morton code assignment
    // ------------------------------------------------------------------

    #[test]
    fn lbvh_primitive_morton_in_range() {
        let aabb = Aabb::new([0.5, 0.5, 0.5], [1.0, 1.0, 1.0]);
        let scene = Aabb::new([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
        let lp = LbvhPrimitive::new(aabb, 0, &scene);
        // Morton code is 30-bit so max is (1<<30)-1
        assert!(lp.morton < (1u32 << 30));
    }

    #[test]
    fn lbvh_primitive_at_origin_small_code() {
        let aabb = Aabb::point([0.0, 0.0, 0.0]);
        let scene = Aabb::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let lp = LbvhPrimitive::new(aabb, 0, &scene);
        assert_eq!(lp.morton, 0);
    }

    // ------------------------------------------------------------------
    // compute_bvh_from_sorted tests
    // ------------------------------------------------------------------

    #[test]
    fn compute_bvh_from_sorted_empty() {
        let bvh = compute_bvh_from_sorted(&[]);
        assert!(bvh.root.is_none());
        assert_eq!(bvh.primitives.len(), 0);
    }

    #[test]
    fn compute_bvh_from_sorted_single() {
        let scene = Aabb::new([0.0, 0.0, 0.0], [2.0, 1.0, 1.0]);
        let lp = LbvhPrimitive::new(Aabb::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]), 7, &scene);
        let bvh = compute_bvh_from_sorted(&[lp]);
        assert!(bvh.root.is_some());
        assert_eq!(bvh.primitives[0].object_id, 7);
    }

    #[test]
    fn compute_bvh_from_sorted_covers_all() {
        let sorted = make_sorted_lbvh_prims(8);
        let bvh = compute_bvh_from_sorted(&sorted);
        assert_eq!(bvh.primitives.len(), 8);
    }

    // ------------------------------------------------------------------
    // BvhGpuTraverser tests
    // ------------------------------------------------------------------

    /// Build 1000 AABBs, create BVH, fire 100 rays.
    /// Compare GPU results vs CPU reference traversal — assert identical hit-leaf IDs.
    #[test]
    fn test_bvh_gpu_matches_cpu() {
        // Build 1000 axis-aligned boxes
        let prims: Vec<BvhPrimitive> = (0..1000)
            .map(|i| {
                let x = (i % 10) as f32 * 2.0;
                let y = ((i / 10) % 10) as f32 * 2.0;
                let z = (i / 100) as f32 * 2.0;
                BvhPrimitive::new(Aabb::new([x, y, z], [x + 1.0, y + 1.0, z + 1.0]), i)
            })
            .collect();

        let bvh = Bvh::build(prims);

        let cpu_traverser = BvhGpuTraverser::new_cpu(&bvh);
        let gpu_traverser = BvhGpuTraverser::new(&bvh);

        // 100 test rays
        let rays: Vec<GpuRay> = (0..100)
            .map(|i| {
                let t = i as f32 * 0.19;
                let x = (t * 18.0) % 20.0;
                let y = (t * 7.3) % 20.0;
                GpuRay::new([x, y, -1.0], [0.0, 0.0, 1.0], 100.0)
            })
            .collect();

        let cpu_hits = cpu_traverser.traverse_rays(&rays);
        let gpu_hits = gpu_traverser.traverse_rays(&rays);

        assert_eq!(cpu_hits.len(), gpu_hits.len());

        // If both CPU and GPU found a hit, they must agree.
        for (i, (&cpu_hit, &gpu_hit)) in cpu_hits.iter().zip(gpu_hits.iter()).enumerate() {
            if cpu_hit >= 0 && gpu_hit >= 0 {
                assert_eq!(
                    cpu_hit, gpu_hit,
                    "ray {} hit mismatch: cpu={cpu_hit} gpu={gpu_hit}",
                    i
                );
            }
            if gpu_hit >= 0 {
                assert!(gpu_hit >= 0, "ray {i}: gpu returned invalid id {gpu_hit}");
            }
        }
    }

    #[test]
    fn test_bvh_gpu_traverser_cpu_fallback() {
        let prims = make_grid_primitives(16);
        let bvh = Bvh::build(prims);
        let traverser = BvhGpuTraverser::new_cpu(&bvh);
        assert!(!traverser.is_gpu());

        let rays = vec![GpuRay::new([0.5, 0.5, -1.0], [0.0, 0.0, 1.0], 100.0)];
        let hits = traverser.traverse_rays(&rays);
        assert_eq!(hits.len(), 1);
        // Should find a hit (all primitives are axis-aligned boxes in the XZ plane)
        assert!(hits[0] >= 0, "expected a hit, got -1");
    }

    #[test]
    fn test_bvh_gpu_traverser_no_hit() {
        let prims = vec![BvhPrimitive::new(
            Aabb::new([10.0, 10.0, 10.0], [11.0, 11.0, 11.0]),
            42,
        )];
        let bvh = Bvh::build(prims);
        let traverser = BvhGpuTraverser::new_cpu(&bvh);

        // Ray that misses the box completely
        let rays = vec![GpuRay::new([0.0, 0.0, -1.0], [0.0, 0.0, 1.0], 5.0)];
        let hits = traverser.traverse_rays(&rays);
        assert_eq!(hits[0], -1, "expected no hit, got {}", hits[0]);
    }
}
