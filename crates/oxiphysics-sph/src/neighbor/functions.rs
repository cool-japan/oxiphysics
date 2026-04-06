//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
use super::types::NeighborPairIterator;

/// Cell key in the spatial hash grid.
pub(super) type CellKey = (i32, i32, i32);
/// Cell key in the spatial hash grid (integer coordinates).
pub(super) type CellKey3 = (i32, i32, i32);
/// Collect all unique symmetric neighbor pairs from per-particle neighbor lists.
///
/// Returns a `Vec<(usize, usize)>` where each pair `(i, j)` satisfies `i < j`.
/// Useful for building symmetric SPH force loops that compute each interaction once.
#[allow(dead_code)]
pub fn collect_symmetric_pairs(neighbor_lists: &[Vec<usize>]) -> Vec<(usize, usize)> {
    NeighborPairIterator::new(neighbor_lists).collect()
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::neighbor::types::*;
    use oxiphysics_core::math::Vec3;
    #[test]
    fn spatial_hash_finds_correct_neighbors() {
        let positions = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.1, 0.0, 0.0),
            Vec3::new(0.5, 0.0, 0.0),
            Vec3::new(5.0, 5.0, 5.0),
        ];
        let hash = SpatialHash::build(&positions, 0.5);
        let nbrs = hash.query_neighbors(&positions, &Vec3::new(0.0, 0.0, 0.0), 0.5);
        assert!(nbrs.contains(&0));
        assert!(nbrs.contains(&1));
        assert!(!nbrs.contains(&3));
    }
    #[test]
    fn find_all_neighbors_excludes_self() {
        let positions = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.1, 0.0, 0.0),
            Vec3::new(10.0, 10.0, 10.0),
        ];
        let all_nbrs = SpatialHash::find_all_neighbors(&positions, 0.5);
        assert!(all_nbrs[0].contains(&1));
        assert!(!all_nbrs[0].contains(&0));
        assert!(!all_nbrs[0].contains(&2));
    }
    #[test]
    fn spatial_hash3d_finds_nearby_particles() {
        let positions: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [0.05, 0.0, 0.0], [10.0, 10.0, 10.0]];
        let hash = SpatialHash3D::build(&positions, 0.2);
        let nbrs = hash.find_neighbors([0.0, 0.0, 0.0]);
        assert!(nbrs.contains(&0), "should find self");
        assert!(nbrs.contains(&1), "should find close particle");
        assert!(!nbrs.contains(&2), "should NOT find far particle");
    }
    #[test]
    fn spatial_hash3d_all_neighbors_excludes_self() {
        let positions: Vec<[f64; 3]> = vec![
            [0.0, 0.0, 0.0],
            [0.05, 0.0, 0.0],
            [0.05, 0.05, 0.0],
            [10.0, 10.0, 10.0],
        ];
        let all_nbrs = SpatialHash3D::all_neighbors(&positions, 0.2);
        assert!(all_nbrs[0].contains(&1));
        assert!(all_nbrs[0].contains(&2));
        assert!(!all_nbrs[0].contains(&0), "self must be excluded");
        assert!(!all_nbrs[0].contains(&3));
    }
    #[test]
    fn spatial_hash3d_empty_input() {
        let positions: Vec<[f64; 3]> = vec![];
        let all_nbrs = SpatialHash3D::all_neighbors(&positions, 0.1);
        assert!(all_nbrs.is_empty());
    }
    #[test]
    fn spatial_hash3d_single_particle_no_neighbors() {
        let positions: Vec<[f64; 3]> = vec![[1.0, 2.0, 3.0]];
        let all_nbrs = SpatialHash3D::all_neighbors(&positions, 0.1);
        assert_eq!(all_nbrs[0].len(), 0, "single particle has no neighbors");
    }
    #[test]
    fn verlet_list_build_and_query() {
        let positions: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [0.05, 0.0, 0.0], [10.0, 10.0, 10.0]];
        let vl = VerletList::build(&positions, 0.2, 0.05);
        assert_eq!(vl.len(), 3);
        assert!(vl.neighbors(0).contains(&1));
        assert!(!vl.neighbors(0).contains(&2));
    }
    #[test]
    fn verlet_list_needs_rebuild_false_initially() {
        let positions: Vec<[f64; 3]> = vec![[0.0; 3], [0.1, 0.0, 0.0]];
        let vl = VerletList::build(&positions, 0.2, 0.1);
        assert!(!vl.needs_rebuild(&positions));
    }
    #[test]
    fn verlet_list_needs_rebuild_after_large_move() {
        let positions: Vec<[f64; 3]> = vec![[0.0; 3], [0.1, 0.0, 0.0]];
        let vl = VerletList::build(&positions, 0.2, 0.1);
        let moved = vec![[0.0; 3], [1.0, 0.0, 0.0]];
        assert!(vl.needs_rebuild(&moved));
    }
    #[test]
    fn verlet_list_rebuild_if_needed() {
        let positions: Vec<[f64; 3]> = vec![[0.0; 3], [0.1, 0.0, 0.0]];
        let mut vl = VerletList::build(&positions, 0.2, 0.1);
        let moved = vec![[0.0; 3], [5.0, 0.0, 0.0]];
        let rebuilt = vl.rebuild_if_needed(&moved);
        assert!(rebuilt);
        assert!(!vl.neighbors(0).contains(&1));
    }
    #[test]
    fn linked_cell_list_finds_neighbors() {
        let positions: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [0.05, 0.0, 0.0], [10.0, 10.0, 10.0]];
        let lcl = LinkedCellList::build(&positions, 0.2);
        let nbrs = lcl.query(&positions, [0.0, 0.0, 0.0], 0.2);
        assert!(nbrs.contains(&0));
        assert!(nbrs.contains(&1));
        assert!(!nbrs.contains(&2));
    }
    #[test]
    fn linked_cell_list_all_neighbors_excludes_self() {
        let positions: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [0.05, 0.0, 0.0], [10.0, 10.0, 10.0]];
        let all = LinkedCellList::all_neighbors(&positions, 0.2);
        assert!(all[0].contains(&1));
        assert!(!all[0].contains(&0), "self must be excluded");
        assert!(!all[0].contains(&2));
    }
    #[test]
    fn linked_cell_list_empty_input() {
        let positions: Vec<[f64; 3]> = vec![];
        let all = LinkedCellList::all_neighbors(&positions, 0.1);
        assert!(all.is_empty());
    }
    #[test]
    fn linked_cell_list_matches_brute_force() {
        let h = 0.15_f64;
        let positions: Vec<[f64; 3]> = (0..5_usize)
            .flat_map(|i| (0..5_usize).map(move |j| [i as f64 * 0.1, j as f64 * 0.1, 0.0]))
            .collect();
        let lcl_nbrs = LinkedCellList::all_neighbors(&positions, h);
        let n = positions.len();
        let mut brute: Vec<Vec<usize>> = vec![Vec::new(); n];
        for i in 0..n {
            for j in 0..n {
                if i == j {
                    continue;
                }
                let d2 = (positions[i][0] - positions[j][0]).powi(2)
                    + (positions[i][1] - positions[j][1]).powi(2)
                    + (positions[i][2] - positions[j][2]).powi(2);
                if d2 <= h * h {
                    brute[i].push(j);
                }
            }
            brute[i].sort_unstable();
        }
        for i in 0..n {
            assert_eq!(
                lcl_nbrs[i], brute[i],
                "LinkedCellList neighbors[{i}] differ from brute force"
            );
        }
    }
    #[test]
    fn linked_cell_list_num_cells_positive() {
        let positions: Vec<[f64; 3]> = vec![[0.0; 3], [1.0, 0.0, 0.0]];
        let lcl = LinkedCellList::build(&positions, 0.5);
        assert!(lcl.num_cells() > 0);
    }
    #[test]
    fn flat_list_build_and_query() {
        let positions: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [0.05, 0.0, 0.0], [10.0, 10.0, 10.0]];
        let flat = FlatNeighborList::build(&positions, 0.2);
        assert_eq!(flat.n_particles, 3);
        assert!(flat.neighbors_of(0).contains(&1));
        assert!(!flat.neighbors_of(0).contains(&2));
        assert!(!flat.neighbors_of(0).contains(&0));
    }
    #[test]
    fn flat_list_offsets_length() {
        let positions: Vec<[f64; 3]> = vec![[0.0; 3], [0.1, 0.0, 0.0], [0.2, 0.0, 0.0]];
        let flat = FlatNeighborList::build(&positions, 0.15);
        assert_eq!(flat.offsets.len(), flat.n_particles + 1);
    }
    #[test]
    fn flat_list_avg_neighbors_positive() {
        let positions: Vec<[f64; 3]> = (0..10).map(|i| [i as f64 * 0.05, 0.0, 0.0]).collect();
        let flat = FlatNeighborList::build(&positions, 0.15);
        assert!(flat.avg_neighbors() > 0.0);
    }
    #[test]
    fn flat_list_empty() {
        let flat = FlatNeighborList::build(&[], 0.1);
        assert_eq!(flat.n_particles, 0);
        assert_eq!(flat.total_pairs(), 0);
        assert_eq!(flat.avg_neighbors(), 0.0);
    }
    #[test]
    fn neighbor_stats_empty() {
        let stats = NeighborStats::compute(&[]);
        assert_eq!(stats.n_particles, 0);
        assert_eq!(stats.avg_neighbors, 0.0);
    }
    #[test]
    fn neighbor_stats_uniform_grid() {
        let positions: Vec<[f64; 3]> = (0..5)
            .flat_map(|i| (0..5).map(move |j| [i as f64 * 0.1, j as f64 * 0.1, 0.0]))
            .collect();
        let stats = NeighborStats::from_positions(&positions, 0.15);
        assert_eq!(stats.n_particles, 25);
        assert!(stats.avg_neighbors > 0.0);
        assert!(stats.max_neighbors >= stats.min_neighbors);
        assert!(stats.std_dev >= 0.0);
    }
    #[test]
    fn neighbor_stats_variance_zero_for_equal_counts() {
        let positions: Vec<[f64; 3]> = (0..5).map(|i| [i as f64 * 10.0, 0.0, 0.0]).collect();
        let nls = SpatialHash3D::all_neighbors(&positions, 0.01);
        let stats = NeighborStats::compute(&nls);
        assert_eq!(stats.min_neighbors, 0);
        assert_eq!(stats.max_neighbors, 0);
        assert!(stats.variance.abs() < 1e-14);
    }
    #[test]
    fn octree_finds_nearby_particles() {
        let positions: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [0.05, 0.0, 0.0], [10.0, 10.0, 10.0]];
        let octree = OctreeNeighborSearch::build(&positions, 0.2);
        let nbrs = octree.query([0.0, 0.0, 0.0]);
        assert!(nbrs.contains(&0));
        assert!(nbrs.contains(&1));
        assert!(!nbrs.contains(&2));
    }
    #[test]
    fn octree_all_neighbors_excludes_self() {
        let positions: Vec<[f64; 3]> = vec![
            [0.0, 0.0, 0.0],
            [0.05, 0.0, 0.0],
            [0.05, 0.05, 0.0],
            [10.0, 10.0, 10.0],
        ];
        let octree = OctreeNeighborSearch::build(&positions, 0.2);
        let all = octree.all_neighbors();
        assert!(all[0].contains(&1));
        assert!(!all[0].contains(&0));
        assert!(!all[0].contains(&3));
    }
    #[test]
    fn octree_empty_positions() {
        let octree = OctreeNeighborSearch::build(&[], 0.1);
        let nbrs = octree.query([0.0, 0.0, 0.0]);
        assert!(nbrs.is_empty());
        let all = octree.all_neighbors();
        assert!(all.is_empty());
    }
    #[test]
    fn octree_agrees_with_brute_force() {
        let h = 0.15_f64;
        let positions: Vec<[f64; 3]> = (0..4)
            .flat_map(|i| (0..4).map(move |j| [i as f64 * 0.1, j as f64 * 0.1, 0.0]))
            .collect();
        let octree = OctreeNeighborSearch::build(&positions, h);
        let octree_nbrs = octree.all_neighbors();
        let n = positions.len();
        let mut brute: Vec<Vec<usize>> = vec![Vec::new(); n];
        for i in 0..n {
            for j in 0..n {
                if i == j {
                    continue;
                }
                let d2 = (positions[i][0] - positions[j][0]).powi(2)
                    + (positions[i][1] - positions[j][1]).powi(2)
                    + (positions[i][2] - positions[j][2]).powi(2);
                if d2 <= h * h {
                    brute[i].push(j);
                }
            }
            brute[i].sort_unstable();
        }
        for i in 0..n {
            assert_eq!(
                octree_nbrs[i], brute[i],
                "OctreeNeighborSearch neighbors[{i}] differ from brute force"
            );
        }
    }
    #[test]
    fn verlet_list_displacement_tracking_tight() {
        let positions = vec![[0.0_f64; 3], [0.1, 0.0, 0.0]];
        let vl = VerletList::build(&positions, 0.2, 0.001);
        let moved = vec![[0.001, 0.0, 0.0], [0.1, 0.0, 0.0]];
        assert!(vl.needs_rebuild(&moved));
    }
    #[test]
    fn verlet_list_no_rebuild_within_skin() {
        let positions = vec![[0.0_f64; 3], [0.1, 0.0, 0.0]];
        let vl = VerletList::build(&positions, 0.2, 0.1);
        let moved = vec![[0.04, 0.0, 0.0], [0.1, 0.0, 0.0]];
        assert!(!vl.needs_rebuild(&moved));
    }
    #[test]
    fn spatial_hash3d_correctness_vs_brute_force() {
        let h = 0.15_f64;
        let positions: Vec<[f64; 3]> = (0..5_usize)
            .flat_map(|i| (0..5_usize).map(move |j| [i as f64 * 0.1, j as f64 * 0.1, 0.0]))
            .collect();
        let grid_nbrs = SpatialHash3D::all_neighbors(&positions, h);
        let n = positions.len();
        let mut brute: Vec<Vec<usize>> = vec![Vec::new(); n];
        for i in 0..n {
            for j in 0..n {
                if i == j {
                    continue;
                }
                let d2 = (positions[i][0] - positions[j][0]).powi(2)
                    + (positions[i][1] - positions[j][1]).powi(2)
                    + (positions[i][2] - positions[j][2]).powi(2);
                if d2 <= h * h {
                    brute[i].push(j);
                }
            }
            brute[i].sort_unstable();
        }
        for i in 0..n {
            assert_eq!(
                grid_nbrs[i], brute[i],
                "SpatialHash3D neighbors[{i}] differ from brute force"
            );
        }
    }
    #[test]
    fn pair_iterator_unique_pairs_triangle() {
        let lists: Vec<Vec<usize>> = vec![vec![1, 2], vec![0, 2], vec![0, 1]];
        let pairs: Vec<(usize, usize)> = NeighborPairIterator::new(&lists).collect();
        assert_eq!(pairs.len(), 3, "should produce 3 unique pairs for 3-clique");
        assert!(pairs.contains(&(0, 1)));
        assert!(pairs.contains(&(0, 2)));
        assert!(pairs.contains(&(1, 2)));
    }
    #[test]
    fn pair_iterator_empty_lists() {
        let lists: Vec<Vec<usize>> = vec![vec![], vec![], vec![]];
        let pairs: Vec<(usize, usize)> = NeighborPairIterator::new(&lists).collect();
        assert!(pairs.is_empty());
    }
    #[test]
    fn pair_iterator_no_duplicates() {
        let positions: Vec<[f64; 3]> = (0..10).map(|i| [i as f64 * 0.05, 0.0, 0.0]).collect();
        let lists = SpatialHash3D::all_neighbors(&positions, 0.15);
        let pairs: Vec<(usize, usize)> = NeighborPairIterator::new(&lists).collect();
        let n_pairs = pairs.len();
        let unique: std::collections::HashSet<_> = pairs.into_iter().collect();
        assert_eq!(unique.len(), n_pairs, "all pairs should be unique");
    }
    #[test]
    fn collect_symmetric_pairs_matches_iterator() {
        let lists: Vec<Vec<usize>> = vec![vec![1, 2], vec![0, 2], vec![0, 1]];
        let a: Vec<_> = NeighborPairIterator::new(&lists).collect();
        let b = collect_symmetric_pairs(&lists);
        assert_eq!(a, b);
    }
    #[test]
    fn periodic_hash_finds_neighbor_across_boundary() {
        let box_size = [1.0_f64; 3];
        let positions = vec![[0.02, 0.5, 0.5], [0.98, 0.5, 0.5], [0.5, 0.5, 0.5]];
        let psh = PeriodicSpatialHash::build(&positions, 0.1, box_size);
        let nbrs_0 = psh.query_neighbors(0);
        assert!(
            nbrs_0.contains(&1),
            "PBC should find particle 1 from particle 0"
        );
        assert!(
            !nbrs_0.contains(&2),
            "Centre particle should be out of range"
        );
    }
    #[test]
    fn periodic_hash_min_image_no_wrap_needed() {
        let psh = PeriodicSpatialHash::build(&[[0.0; 3]], 0.1, [1.0; 3]);
        let dr = psh.min_image([0.3, 0.4, 0.5], [0.1, 0.1, 0.1]);
        assert!((dr[0] - 0.2).abs() < 1e-14);
        assert!((dr[1] - 0.3).abs() < 1e-14);
        assert!((dr[2] - 0.4).abs() < 1e-14);
    }
    #[test]
    fn periodic_hash_min_image_wrap_needed() {
        let psh = PeriodicSpatialHash::build(&[[0.0; 3]], 0.2, [1.0; 3]);
        let dr = psh.min_image([0.05, 0.0, 0.0], [0.95, 0.0, 0.0]);
        assert!(
            dr[0].abs() < 0.11,
            "min-image distance should be ~0.10, got |{}|",
            dr[0]
        );
    }
    #[test]
    fn periodic_all_neighbors_symmetric() {
        let positions: Vec<[f64; 3]> = (0..5).map(|i| [i as f64 * 0.2, 0.0, 0.0]).collect();
        let box_size = [1.0_f64, 1.0, 1.0];
        let nls = PeriodicSpatialHash::all_neighbors(&positions, 0.25, box_size);
        for (i, list) in nls.iter().enumerate() {
            for &j in list {
                assert!(
                    nls[j].contains(&i),
                    "PBC neighbor list should be symmetric: {j} should list {i}"
                );
            }
        }
    }
    #[test]
    fn adaptive_search_same_h_matches_fixed() {
        let positions: Vec<[f64; 3]> = (0..5)
            .flat_map(|i| (0..5).map(move |j| [i as f64 * 0.1, j as f64 * 0.1, 0.0]))
            .collect();
        let h = 0.15_f64;
        let h_vec = vec![h; positions.len()];
        let ans = AdaptiveNeighborSearch::build(&positions, &h_vec);
        let fixed_nls = SpatialHash3D::all_neighbors(&positions, h);
        for i in 0..positions.len() {
            let mut adaptive = ans.query(i);
            let mut fixed = fixed_nls[i].clone();
            adaptive.sort_unstable();
            fixed.sort_unstable();
            assert_eq!(
                adaptive, fixed,
                "adaptive search[{i}] should match fixed-h search"
            );
        }
    }
    #[test]
    fn adaptive_search_larger_h_captures_more_neighbors() {
        let positions = vec![[0.0_f64, 0.0, 0.0], [0.1, 0.0, 0.0], [0.3, 0.0, 0.0]];
        let h_small = vec![0.05_f64, 0.05, 0.05];
        let h_large = vec![0.35_f64, 0.35, 0.35];
        let ans_small = AdaptiveNeighborSearch::build(&positions, &h_small);
        let ans_large = AdaptiveNeighborSearch::build(&positions, &h_large);
        assert!(!ans_small.query(0).contains(&2));
        assert!(ans_large.query(0).contains(&2));
    }
    #[test]
    fn adaptive_suggest_h_increases_when_too_few_neighbors() {
        let positions: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [4.0, 0.0, 0.0]];
        let h_vec = vec![0.1_f64; 3];
        let ans = AdaptiveNeighborSearch::build(&positions, &h_vec);
        let new_h = ans.suggest_h(0, 10);
        assert!(new_h > 0.0, "suggested h should be positive");
    }
    #[test]
    fn neighbor_list_cache_no_rebuild_for_small_displacement() {
        let positions: Vec<[f64; 3]> = vec![[0.0; 3], [0.1, 0.0, 0.0]];
        let mut cache = NeighborListCache::new(&positions, 0.2, 0.1);
        let moved = vec![[0.01, 0.0, 0.0], [0.1, 0.0, 0.0]];
        let rebuilt = cache.update(&moved);
        assert!(!rebuilt, "should NOT rebuild for tiny displacement");
        assert_eq!(cache.rebuild_count, 0);
    }
    #[test]
    fn neighbor_list_cache_rebuilds_after_large_displacement() {
        let positions: Vec<[f64; 3]> = vec![[0.0; 3], [0.1, 0.0, 0.0]];
        let mut cache = NeighborListCache::new(&positions, 0.2, 0.1);
        let far = vec![[5.0, 0.0, 0.0], [0.1, 0.0, 0.0]];
        let rebuilt = cache.update(&far);
        assert!(rebuilt, "should rebuild after large displacement");
        assert_eq!(cache.rebuild_count, 1);
        assert!(!cache.neighbors_of(0).contains(&1));
    }
    #[test]
    fn neighbor_list_cache_len_matches_particle_count() {
        let positions: Vec<[f64; 3]> = vec![[0.0; 3], [0.1, 0.0, 0.0], [0.2, 0.0, 0.0]];
        let cache = NeighborListCache::new(&positions, 0.15, 0.05);
        assert_eq!(cache.len(), 3);
        assert!(!cache.is_empty());
    }
    #[test]
    fn morton_code_origin_is_zero() {
        assert_eq!(ZOrderCurve::morton_code(0, 0, 0), 0);
    }
    #[test]
    fn morton_code_axes_are_distinct() {
        let mx = ZOrderCurve::morton_code(1, 0, 0);
        let my = ZOrderCurve::morton_code(0, 1, 0);
        let mz = ZOrderCurve::morton_code(0, 0, 1);
        assert_ne!(mx, my);
        assert_ne!(mx, mz);
        assert_ne!(my, mz);
    }
    #[test]
    fn morton_sort_produces_valid_permutation() {
        let positions: Vec<[f64; 3]> = (0..8)
            .map(|i| [(i % 2) as f64, ((i / 2) % 2) as f64, (i / 4) as f64])
            .collect();
        let lo = [0.0; 3];
        let hi = [2.0; 3];
        let perm = ZOrderCurve::sort_by_morton(&positions, lo, hi);
        assert_eq!(perm.len(), 8);
        let mut sorted = perm.clone();
        sorted.sort_unstable();
        assert_eq!(sorted, (0..8).collect::<Vec<_>>());
    }
    #[test]
    fn hash_grid_stats_empty_input() {
        let stats = HashGridStats::from_positions(&[], 0.1);
        assert_eq!(stats.total_cells, 0);
        assert_eq!(stats.occupied_cells, 0);
        assert_eq!(stats.max_cell_occupancy, 0);
    }
    #[test]
    fn hash_grid_stats_single_cell() {
        let positions: Vec<[f64; 3]> = vec![[0.01; 3], [0.02; 3], [0.03; 3]];
        let stats = HashGridStats::from_positions(&positions, 0.1);
        assert_eq!(stats.occupied_cells, 1);
        assert_eq!(stats.max_cell_occupancy, 3);
        assert!((stats.avg_cell_occupancy - 3.0).abs() < 1e-12);
    }
    #[test]
    fn hash_grid_stats_spread_particles() {
        let positions: Vec<[f64; 3]> = (0..5).map(|i| [i as f64, 0.0, 0.0]).collect();
        let stats = HashGridStats::from_positions(&positions, 0.5);
        assert_eq!(stats.max_cell_occupancy, 1);
        assert!((stats.avg_cell_occupancy - 1.0).abs() < 1e-12);
    }
}
#[cfg(test)]
mod tests_neighbor_new {

    use crate::neighbor::types::*;
    #[test]
    fn skin_manager_no_rebuild_initially() {
        let positions = vec![[0.0_f64; 3], [1.0, 0.0, 0.0]];
        let mgr = SkinManager::new(&positions, 1.0, 0.2);
        assert!(
            !mgr.needs_rebuild(&positions),
            "no displacement → no rebuild"
        );
    }
    #[test]
    fn skin_manager_needs_rebuild_after_large_move() {
        let positions = vec![[0.0_f64; 3], [1.0, 0.0, 0.0]];
        let mgr = SkinManager::new(&positions, 1.0, 0.2);
        let moved = vec![[0.5_f64, 0.0, 0.0], [1.0, 0.0, 0.0]];
        assert!(mgr.needs_rebuild(&moved));
    }
    #[test]
    fn skin_manager_no_rebuild_within_half_skin() {
        let positions = vec![[0.0_f64; 3]];
        let mgr = SkinManager::new(&positions, 1.0, 0.2);
        let moved = vec![[0.09_f64, 0.0, 0.0]];
        assert!(!mgr.needs_rebuild(&moved));
    }
    #[test]
    fn skin_manager_commit_rebuild_increments_count() {
        let positions = vec![[0.0_f64; 3]];
        let mut mgr = SkinManager::new(&positions, 1.0, 0.2);
        let moved = vec![[1.0_f64, 0.0, 0.0]];
        mgr.commit_rebuild(&moved);
        assert_eq!(mgr.rebuild_count, 1);
    }
    #[test]
    fn skin_manager_update_rebuilds_when_needed() {
        let positions = vec![[0.0_f64; 3]];
        let mut mgr = SkinManager::new(&positions, 1.0, 0.2);
        let far = vec![[5.0_f64, 0.0, 0.0]];
        let rebuilt = mgr.update(&far);
        assert!(rebuilt);
        assert_eq!(mgr.rebuild_count, 1);
        assert_eq!(mgr.update_count, 1);
    }
    #[test]
    fn skin_manager_max_displacement_correct() {
        let positions = vec![[0.0_f64; 3], [0.0, 0.0, 0.0]];
        let mgr = SkinManager::new(&positions, 1.0, 0.2);
        let moved = vec![[0.3_f64, 0.4, 0.0], [0.0, 0.0, 0.0]];
        let max_d = mgr.max_displacement(&moved);
        assert!((max_d - 0.5).abs() < 1e-12, "max_d={max_d}");
    }
    #[test]
    fn parallel_builder_matches_linked_cell_list() {
        let positions: Vec<[f64; 3]> = (0..8).map(|i| [i as f64 * 0.1, 0.0, 0.0]).collect();
        let builder = ParallelReadyNeighborBuilder::new(0.15, 0.05);
        let r_total = builder.r_cut + builder.r_skin;
        let nls1 = builder.build_serial(&positions);
        let nls2 = LinkedCellList::all_neighbors(&positions, r_total);
        for i in 0..positions.len() {
            assert_eq!(nls1[i], nls2[i], "mismatch at particle {i}");
        }
    }
    #[test]
    fn parallel_builder_total_pairs_symmetric() {
        let positions: Vec<[f64; 3]> = (0..6).map(|i| [i as f64 * 0.1, 0.0, 0.0]).collect();
        let builder = ParallelReadyNeighborBuilder::new(0.15, 0.0);
        let nls = builder.build_serial(&positions);
        let half_sum = nls.iter().map(|v| v.len()).sum::<usize>() / 2;
        let total = builder.total_pairs(&positions);
        assert_eq!(half_sum, total);
    }
    #[test]
    fn parallel_builder_statistics_sensible() {
        let positions: Vec<[f64; 3]> = (0..10).map(|i| [i as f64 * 0.1, 0.0, 0.0]).collect();
        let builder = ParallelReadyNeighborBuilder::new(0.15, 0.05);
        let stats = builder.statistics(&positions);
        assert_eq!(stats.n_particles, 10);
        assert!(stats.avg_neighbors > 0.0);
    }
    #[test]
    fn pbc_lcl_finds_neighbor_across_x_boundary() {
        let box_size = [1.0_f64, 1.0, 1.0];
        let positions = vec![[0.05, 0.5, 0.5], [0.95, 0.5, 0.5], [0.5, 0.5, 0.5]];
        let nbrs = PbcLinkedCellList::all_neighbors(&positions, 0.15, box_size);
        assert!(
            nbrs[0].contains(&1),
            "PBC: particle 0 should see particle 1"
        );
        assert!(nbrs[1].contains(&0), "PBC: symmetry should hold");
        assert!(
            !nbrs[0].contains(&2),
            "centre particle should not be in range"
        );
    }
    #[test]
    fn pbc_lcl_symmetric_neighbor_lists() {
        let box_size = [2.0_f64, 2.0, 2.0];
        let positions: Vec<[f64; 3]> = (0..8).map(|i| [(i as f64) * 0.25, 0.0, 0.0]).collect();
        let nbrs = PbcLinkedCellList::all_neighbors(&positions, 0.3, box_size);
        for (i, list) in nbrs.iter().enumerate() {
            for &j in list {
                assert!(
                    nbrs[j].contains(&i),
                    "PBC neighbor list asymmetry: {j} doesn't list {i}"
                );
            }
        }
    }
    #[test]
    fn pbc_lcl_num_cells_positive() {
        let positions = vec![[0.0_f64; 3], [0.5, 0.0, 0.0]];
        let lcl = PbcLinkedCellList::build(&positions, 0.3, [1.0; 3]);
        assert!(lcl.num_cells() > 0);
    }
    #[test]
    fn pbc_lcl_wraps_position_into_box() {
        let positions = vec![[1.05, 0.5, 0.5], [0.05, 0.5, 0.5]];
        let box_size = [1.0_f64, 1.0, 1.0];
        let nbrs = PbcLinkedCellList::all_neighbors(&positions, 0.2, box_size);
        assert!(
            nbrs[0].contains(&1),
            "wrapped position should find neighbor"
        );
    }
    #[test]
    fn neighbor_count_stats_empty() {
        let stats = NeighborCountStats::from_neighbor_lists(&[]);
        assert_eq!(stats.min, 0);
        assert_eq!(stats.max, 0);
        assert_eq!(stats.mean, 0.0);
    }
    #[test]
    fn neighbor_count_stats_uniform_grid() {
        let positions: Vec<[f64; 3]> = (0..5)
            .flat_map(|i| (0..5).map(move |j| [i as f64 * 0.1, j as f64 * 0.1, 0.0]))
            .collect();
        let stats = NeighborCountStats::from_positions(&positions, 0.15);
        assert_eq!(stats.counts.len(), 25);
        assert!(stats.mean > 0.0);
        assert!(stats.max >= stats.min);
        assert!(stats.std_dev >= 0.0);
    }
    #[test]
    fn neighbor_count_stats_fraction_with_at_least() {
        let nls = vec![vec![1usize, 2], vec![0usize, 2], vec![0usize]];
        let stats = NeighborCountStats::from_neighbor_lists(&nls);
        assert!((stats.fraction_with_at_least(1) - 1.0).abs() < 1e-10);
        let f2 = stats.fraction_with_at_least(2);
        assert!((f2 - 2.0 / 3.0).abs() < 1e-10, "f2={f2}");
        assert_eq!(stats.fraction_with_at_least(3), 0.0);
    }
    #[test]
    fn neighbor_count_stats_isolated_particles_have_zero_neighbors() {
        let positions: Vec<[f64; 3]> = (0..5).map(|i| [i as f64 * 10.0, 0.0, 0.0]).collect();
        let stats = NeighborCountStats::from_positions(&positions, 0.1);
        assert_eq!(stats.min, 0);
        assert_eq!(stats.max, 0);
        assert_eq!(stats.mean, 0.0);
    }
    #[test]
    fn verlet_list_rebuild_updates_ref_positions() {
        let positions = vec![[0.0_f64; 3], [0.1, 0.0, 0.0]];
        let mut vl = VerletList::build(&positions, 0.2, 0.05);
        let new_positions = vec![[0.5_f64; 3], [0.6, 0.0, 0.0]];
        vl.rebuild(&new_positions);
        for (i, pos) in new_positions.iter().enumerate() {
            for k in 0..3 {
                assert!((vl.ref_positions[i][k] - pos[k]).abs() < 1e-14);
            }
        }
    }
    #[test]
    fn verlet_list_is_empty_for_empty_input() {
        let vl = VerletList::build(&[], 0.2, 0.05);
        assert!(vl.is_empty());
        assert_eq!(vl.len(), 0);
    }
    #[test]
    fn verlet_list_skin_too_small_always_needs_rebuild() {
        let positions = vec![[0.0_f64; 3], [0.1, 0.0, 0.0]];
        let vl = VerletList::build(&positions, 0.2, 0.0);
        let moved = vec![[0.0001_f64; 3], [0.1, 0.0, 0.0]];
        assert!(
            vl.needs_rebuild(&moved),
            "Zero skin should require rebuild on any move"
        );
    }
    #[test]
    fn linked_cell_list_single_particle_no_neighbors() {
        let positions = vec![[5.0_f64, 5.0, 5.0]];
        let all = LinkedCellList::all_neighbors(&positions, 0.5);
        assert_eq!(all.len(), 1);
        assert!(all[0].is_empty(), "Single particle has no neighbors");
    }
    #[test]
    fn linked_cell_list_chain_of_particles() {
        let positions: Vec<[f64; 3]> = (0..10).map(|i| [i as f64 * 0.1, 0.0, 0.0]).collect();
        let all = LinkedCellList::all_neighbors(&positions, 0.15);
        for i in 1..9 {
            assert!(
                all[i].len() >= 2,
                "Interior particle {i} should have ≥2 neighbors"
            );
        }
    }
    #[test]
    fn flat_neighbor_list_from_neighbor_lists_preserves_counts() {
        let nls = vec![vec![1usize, 2], vec![0usize, 2], vec![0usize, 1]];
        let flat = FlatNeighborList::from_neighbor_lists(&nls);
        for i in 0..3 {
            assert_eq!(flat.neighbors_of(i).len(), 2);
        }
    }
    #[test]
    fn flat_neighbor_list_total_pairs_matches() {
        let nls = vec![vec![1usize], vec![0usize, 2], vec![1usize]];
        let flat = FlatNeighborList::from_neighbor_lists(&nls);
        assert_eq!(flat.total_pairs(), 4);
    }
    #[test]
    fn flat_neighbor_list_avg_neighbors_correct() {
        let nls = vec![vec![1usize, 2], vec![0usize], vec![0usize]];
        let flat = FlatNeighborList::from_neighbor_lists(&nls);
        assert!((flat.avg_neighbors() - 4.0 / 3.0).abs() < 1e-12);
    }
    #[test]
    fn neighbor_stats_single_particle_zero_neighbors() {
        let nls = vec![vec![]];
        let stats = NeighborStats::compute(&nls);
        assert_eq!(stats.n_particles, 1);
        assert_eq!(stats.min_neighbors, 0);
        assert_eq!(stats.max_neighbors, 0);
        assert_eq!(stats.avg_neighbors, 0.0);
    }
    #[test]
    fn neighbor_stats_total_pairs_is_sum_of_counts() {
        let nls = vec![vec![1usize, 2], vec![0usize, 2], vec![0usize, 1]];
        let stats = NeighborStats::compute(&nls);
        assert_eq!(stats.total_pairs, 6);
    }
    #[test]
    fn neighbor_stats_variance_non_negative() {
        let positions: Vec<[f64; 3]> = (0..10).map(|i| [i as f64 * 0.1, 0.0, 0.0]).collect();
        let stats = NeighborStats::from_positions(&positions, 0.25);
        assert!(stats.variance >= 0.0);
        assert!(stats.std_dev >= 0.0);
    }
    #[test]
    fn spatial_hash3d_build_and_query_self_included() {
        let positions = vec![[0.5_f64; 3]];
        let hash = SpatialHash3D::build(&positions, 0.1);
        let nbrs = hash.find_neighbors([0.5, 0.5, 0.5]);
        assert!(
            nbrs.contains(&0),
            "Query at particle position should find self"
        );
    }
    #[test]
    fn spatial_hash3d_distance_exactly_h_is_included() {
        let positions = vec![[0.0_f64; 3], [0.1, 0.0, 0.0]];
        let h = 0.1;
        let hash = SpatialHash3D::build(&positions, h);
        let nbrs = hash.find_neighbors([0.0, 0.0, 0.0]);
        assert!(
            nbrs.contains(&1),
            "Particle exactly at distance h should be included"
        );
    }
    #[test]
    fn periodic_hash_all_neighbors_respects_pbc() {
        let box_size = [1.0_f64; 3];
        let positions = vec![[0.01, 0.5, 0.5], [0.99, 0.5, 0.5], [0.5, 0.5, 0.5]];
        let all = PeriodicSpatialHash::all_neighbors(&positions, 0.05, box_size);
        assert!(all[0].contains(&1), "PBC: particle 0 should see particle 1");
        assert!(all[1].contains(&0), "PBC: particle 1 should see particle 0");
        assert!(!all[0].contains(&2));
    }
    #[test]
    fn periodic_hash_min_image_wraps_correctly() {
        let psh = PeriodicSpatialHash::build(&[[0.0; 3]], 0.1, [1.0; 3]);
        let dr = psh.min_image([0.05, 0.0, 0.0], [0.95, 0.0, 0.0]);
        assert!((dr[0] - 0.10).abs() < 1e-12, "min-image dr[0]={}", dr[0]);
    }
    #[test]
    fn adaptive_neighbor_search_finds_particle_within_h_i() {
        let positions = vec![[0.0_f64; 3], [0.09, 0.0, 0.0]];
        let h = vec![0.1_f64, 0.1];
        let search = AdaptiveNeighborSearch::build(&positions, &h);
        let nbrs = search.query(0);
        assert!(
            nbrs.contains(&1),
            "Particle 1 within h_i=0.1 should be found"
        );
    }
    #[test]
    fn adaptive_neighbor_search_excludes_self() {
        let positions = vec![[0.0_f64; 3], [0.05, 0.0, 0.0], [10.0, 0.0, 0.0]];
        let h = vec![0.1; 3];
        let search = AdaptiveNeighborSearch::build(&positions, &h);
        for i in 0..3 {
            let nbrs = search.query(i);
            assert!(
                !nbrs.contains(&i),
                "Self must not appear in neighbor list for particle {i}"
            );
        }
    }
    #[test]
    fn adaptive_neighbor_search_suggest_h_larger_for_few_neighbors() {
        let positions = vec![[0.0_f64; 3], [100.0, 0.0, 0.0]];
        let h = vec![0.1_f64, 0.1];
        let search = AdaptiveNeighborSearch::build(&positions, &h);
        let h_new = search.suggest_h(0, 10);
        assert!(
            (h_new - h[0]).abs() < 1e-14,
            "With 0 neighbors, suggest_h should return h unchanged"
        );
    }
    #[test]
    fn neighbor_list_cache_rebuild_count_increments() {
        let positions = vec![[0.0_f64; 3], [0.1, 0.0, 0.0]];
        let mut cache = NeighborListCache::new(&positions, 0.2, 0.05);
        assert_eq!(cache.rebuild_count, 0);
        let far = vec![[0.0_f64; 3], [5.0, 0.0, 0.0]];
        let rebuilt = cache.update(&far);
        assert!(rebuilt);
        assert_eq!(cache.rebuild_count, 1);
    }
    #[test]
    fn neighbor_list_cache_no_rebuild_when_close() {
        let positions = vec![[0.0_f64; 3], [0.1, 0.0, 0.0]];
        let mut cache = NeighborListCache::new(&positions, 0.2, 0.1);
        let close = vec![[0.01_f64; 3], [0.1, 0.0, 0.0]];
        let rebuilt = cache.update(&close);
        assert!(!rebuilt);
        assert_eq!(cache.rebuild_count, 0);
    }
    #[test]
    fn z_order_curve_origin_is_zero() {
        let code = ZOrderCurve::morton_code(0, 0, 0);
        assert_eq!(code, 0);
    }
    #[test]
    fn z_order_curve_sort_returns_permutation_length() {
        let positions: Vec<[f64; 3]> = (0..8).map(|i| [i as f64 * 0.1, 0.0, 0.0]).collect();
        let lo = [0.0; 3];
        let hi = [1.0; 3];
        let perm = ZOrderCurve::sort_by_morton(&positions, lo, hi);
        assert_eq!(perm.len(), 8);
        let mut sorted = perm.clone();
        sorted.sort_unstable();
        let expected: Vec<usize> = (0..8).collect();
        assert_eq!(sorted, expected);
    }
    #[test]
    fn hash_grid_stats_from_positions_non_empty() {
        let positions: Vec<[f64; 3]> = (0..10).map(|i| [i as f64 * 0.1, 0.0, 0.0]).collect();
        let stats = HashGridStats::from_positions(&positions, 0.2);
        assert!(stats.occupied_cells > 0);
        assert!(stats.max_cell_occupancy > 0);
        assert!(stats.avg_cell_occupancy > 0.0);
    }
    #[test]
    fn hash_grid_stats_from_positions_empty() {
        let stats = HashGridStats::from_positions(&[], 0.1);
        assert_eq!(stats.occupied_cells, 0);
        assert_eq!(stats.max_cell_occupancy, 0);
        assert_eq!(stats.avg_cell_occupancy, 0.0);
    }
}
#[cfg(test)]
mod tests_neighbor_ext {

    use crate::neighbor::types::*;
    #[test]
    fn rebuild_frequency_fast_particles_is_one() {
        let vl = VerletList::build(&[[0.0; 3]], 0.1, 0.02);
        let n = vl.compute_rebuild_frequency(1e9, 1.0);
        assert_eq!(n, 1);
    }
    #[test]
    fn rebuild_frequency_slow_particles_is_large() {
        let vl = VerletList::build(&[[0.0; 3]], 0.1, 0.5);
        let n = vl.compute_rebuild_frequency(0.01, 0.001);
        assert!(
            (24_999..=25_001).contains(&n),
            "Expected ~25000 steps, got {n}",
        );
    }
    #[test]
    fn rebuild_frequency_zero_velocity_is_one() {
        let vl = VerletList::build(&[[0.0; 3]], 0.1, 0.1);
        let n = vl.compute_rebuild_frequency(0.0, 0.001);
        assert_eq!(n, 1);
    }
    #[test]
    fn rebuild_frequency_zero_dt_is_one() {
        let vl = VerletList::build(&[[0.0; 3]], 0.1, 0.1);
        let n = vl.compute_rebuild_frequency(1.0, 0.0);
        assert_eq!(n, 1);
    }
    #[test]
    fn displacement_since_build_stationary_is_zero() {
        let positions = vec![[0.0_f64; 3], [1.0, 0.0, 0.0]];
        let vl = VerletList::build(&positions, 0.2, 0.05);
        let d = vl.compute_displacement_since_build(&positions);
        assert!((d).abs() < 1e-14);
    }
    #[test]
    fn displacement_since_build_one_moved_particle() {
        let ref_positions = vec![[0.0_f64; 3], [1.0, 0.0, 0.0]];
        let vl = VerletList::build(&ref_positions, 0.2, 0.05);
        let new_positions = vec![[0.0_f64; 3], [1.0, 0.03, 0.04]];
        let d = vl.compute_displacement_since_build(&new_positions);
        assert!(
            (d - 0.05).abs() < 1e-12,
            "max displacement should be 0.05, got {d}"
        );
    }
    #[test]
    fn displacement_since_build_empty_is_zero() {
        let vl = VerletList::build(&[], 0.1, 0.05);
        let d = vl.compute_displacement_since_build(&[]);
        assert_eq!(d, 0.0);
    }
    #[test]
    fn occupancy_histogram_empty_list_peak_zero() {
        let lcl = LinkedCellList::build(&[], 0.1);
        let hist = lcl.compute_occupancy_histogram();
        assert_eq!(hist.peak_occupancy, 0);
        assert_eq!(hist.occupied_cells, 0);
    }
    #[test]
    fn occupancy_histogram_single_particle_one_occupied_cell() {
        let positions = vec![[0.5_f64; 3]];
        let lcl = LinkedCellList::build(&positions, 0.1);
        let hist = lcl.compute_occupancy_histogram();
        assert_eq!(hist.occupied_cells, 1);
        assert_eq!(hist.peak_occupancy, 1);
        assert!(hist.bins.len() >= 2);
        assert_eq!(hist.bins[1], 1);
    }
    #[test]
    fn occupancy_histogram_all_particles_same_cell() {
        let positions: Vec<[f64; 3]> = (0..4).map(|_| [0.05, 0.05, 0.05]).collect();
        let lcl = LinkedCellList::build(&positions, 1.0);
        let hist = lcl.compute_occupancy_histogram();
        assert_eq!(hist.peak_occupancy, 4);
        assert_eq!(hist.occupied_cells, 1);
    }
    #[test]
    fn occupancy_histogram_bins_sum_equals_total_cells() {
        let positions: Vec<[f64; 3]> = (0..8).map(|i| [i as f64 * 0.2, 0.0, 0.0]).collect();
        let lcl = LinkedCellList::build(&positions, 0.15);
        let hist = lcl.compute_occupancy_histogram();
        let total_cells = lcl.num_cells();
        let bins_sum: usize = hist.bins.iter().sum();
        assert_eq!(bins_sum, total_cells, "bins must sum to total cells");
    }
}
