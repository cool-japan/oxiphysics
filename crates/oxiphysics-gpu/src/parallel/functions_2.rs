//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
#[allow(unused_imports)]
use super::functions::*;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parallel::LoadBalanceStrategy;
    use crate::parallel::WorkChunker;
    use crate::parallel::WorkGroupConfig;
    use crate::parallel::WorkStealQueue;
    use std::sync::atomic::{AtomicUsize, Ordering};
    #[test]
    fn parallel_for_processes_all_items() {
        let counter = AtomicUsize::new(0);
        let n = 100;
        parallel_for(n, 16, |_i| {
            counter.fetch_add(1, Ordering::Relaxed);
        });
        assert_eq!(counter.load(Ordering::Relaxed), n);
    }
    #[test]
    fn test_parallel_for_produces_correct_results() {
        let n = 64;
        let mut results = vec![0.0f64; n];
        let ptr = results.as_mut_ptr();
        parallel_for(n, 8, |i| unsafe {
            *ptr.add(i) = (i as f64) * (i as f64);
        });
        for (i, &val) in results.iter().enumerate() {
            let expected = (i as f64) * (i as f64);
            assert!(
                (val - expected).abs() < 1e-15,
                "index {i}: expected {expected}, got {val}"
            );
        }
    }
    /// Gaussian-like kernel used in density tests: W(r, h) = exp(-r^2/h^2) / h^3
    fn gauss_kernel(r: f64, h: f64) -> f64 {
        (-(r / h) * (r / h)).exp() / (h * h * h)
    }
    #[test]
    fn test_parallel_density_uniform() {
        let spacing = 0.5_f64;
        let h = 0.6_f64;
        let mut positions: Vec<[f64; 3]> = Vec::new();
        for ix in 0..3_i32 {
            for iy in 0..3_i32 {
                for iz in 0..3_i32 {
                    positions.push([
                        ix as f64 * spacing,
                        iy as f64 * spacing,
                        iz as f64 * spacing,
                    ]);
                }
            }
        }
        let n = positions.len();
        let mass = 1.0_f64;
        let masses = vec![mass; n];
        let densities = parallel_sph_density(&positions, &masses, h, gauss_kernel);
        assert_eq!(densities.len(), n);
        for &rho in &densities {
            assert!(rho > 0.0, "density should be positive, got {rho}");
        }
        let self_contrib = 1.0 / (h * h * h);
        assert!(
            densities[13] >= self_contrib * 0.9,
            "interior density too low: {}",
            densities[13]
        );
    }
    #[test]
    fn test_parallel_lj_repulsion() {
        let sigma = 1.0_f64;
        let epsilon = 1.0_f64;
        let cutoff = 3.0_f64;
        let positions = vec![[0.0, 0.0, 0.0], [0.9 * sigma, 0.0, 0.0]];
        let forces = parallel_lj_forces(&positions, epsilon, sigma, cutoff);
        assert_eq!(forces.len(), 2);
        assert!(
            forces[0][0] < 0.0,
            "expected repulsive force on particle 0 in -x, got {}",
            forces[0][0]
        );
        assert!(
            (forces[0][0] + forces[1][0]).abs() < 1e-12,
            "forces not equal and opposite: {} vs {}",
            forces[0][0],
            forces[1][0]
        );
    }
    #[test]
    fn test_parallel_lj_attraction() {
        let sigma = 1.0_f64;
        let epsilon = 1.0_f64;
        let cutoff = 5.0_f64;
        let r_eq = 2.0_f64.powf(1.0 / 6.0) * sigma;
        let positions = vec![[0.0, 0.0, 0.0], [r_eq, 0.0, 0.0]];
        let forces = parallel_lj_forces(&positions, epsilon, sigma, cutoff);
        assert_eq!(forces.len(), 2);
        for k in 0..3 {
            assert!(
                forces[0][k].abs() < 1e-10,
                "force[0][{k}] should be ~0 at equilibrium, got {}",
                forces[0][k]
            );
        }
    }
    #[test]
    fn test_parallel_verlet() {
        let mut positions = vec![[0.0_f64, 0.0, 0.0]];
        let mut velocities = vec![[0.0_f64, 0.0, 0.0]];
        let forces = vec![[0.0_f64, -9.81, 0.0]];
        let masses = vec![1.0_f64];
        let dt = 0.1_f64;
        parallel_verlet_step(&mut positions, &mut velocities, &forces, &masses, dt);
        let expected_y = 0.5 * (-9.81) * dt * dt;
        let expected_vy = 0.5 * (-9.81) * dt;
        assert!(
            (positions[0][1] - expected_y).abs() < 1e-12,
            "y position: expected {expected_y}, got {}",
            positions[0][1]
        );
        assert!(
            (velocities[0][1] - expected_vy).abs() < 1e-12,
            "vy: expected {expected_vy}, got {}",
            velocities[0][1]
        );
        assert!((positions[0][0]).abs() < 1e-15);
        assert!((positions[0][2]).abs() < 1e-15);
    }
    #[test]
    fn test_parallel_aabb_pairs() {
        let aabbs = vec![
            ([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]),
            ([0.5, 0.5, 0.5], [1.5, 1.5, 1.5]),
            ([0.3, 0.3, 0.3], [1.8, 1.8, 1.8]),
            ([0.6, 0.6, 0.6], [2.5, 2.5, 2.5]),
        ];
        let mut pairs = parallel_aabb_pairs(&aabbs);
        pairs.sort_unstable();
        assert_eq!(
            pairs.len(),
            6,
            "expected 6 overlapping pairs, got {}: {:?}",
            pairs.len(),
            pairs
        );
        for &(a, b) in &pairs {
            assert!(a < b, "pair ({a}, {b}) not in canonical order");
        }
    }
    #[test]
    fn test_work_chunker() {
        let n = 100;
        let chunker = WorkChunker::new(n);
        let chunks = chunker.chunks();
        assert!(!chunks.is_empty(), "should produce at least one chunk");
        let mut covered = vec![false; n];
        for range in &chunks {
            for idx in range.clone() {
                assert!(idx < n, "index {idx} out of bounds");
                assert!(!covered[idx], "index {idx} covered twice");
                covered[idx] = true;
            }
        }
        assert!(
            covered.iter().all(|&c| c),
            "not all indices covered by chunks"
        );
    }
    #[test]
    fn test_work_group_config_new() {
        let cfg = WorkGroupConfig::new(128);
        assert_eq!(cfg.preferred_size, 128);
        assert_eq!(cfg.max_size, 1024);
        assert_eq!(cfg.min_size, 32);
    }
    #[test]
    fn test_work_group_config_zero() {
        let cfg = WorkGroupConfig::new(0);
        assert_eq!(cfg.preferred_size, 1);
    }
    #[test]
    fn test_work_group_optimal_small() {
        let cfg = WorkGroupConfig::new(64);
        let size = cfg.optimal_size(10);
        assert!(size >= cfg.min_size);
        assert!(size <= cfg.max_size);
    }
    #[test]
    fn test_work_group_optimal_large() {
        let cfg = WorkGroupConfig::new(64);
        let size = cfg.optimal_size(1000);
        assert!(size >= cfg.min_size);
        assert!(size <= cfg.max_size);
    }
    #[test]
    fn test_work_group_ranges_cover_all() {
        let cfg = WorkGroupConfig::new(64);
        let total = 200;
        let ranges = cfg.group_ranges(total);
        let mut covered = vec![false; total];
        for range in &ranges {
            for idx in range.clone() {
                assert!(!covered[idx], "index {idx} covered twice");
                covered[idx] = true;
            }
        }
        assert!(covered.iter().all(|&c| c));
    }
    #[test]
    fn test_work_group_num_groups() {
        let cfg = WorkGroupConfig::new(64);
        let ng = cfg.num_groups(256);
        assert!(ng >= 1);
        assert!(ng * cfg.optimal_size(256) >= 256);
    }
    #[test]
    fn test_work_group_cpu_default() {
        let cfg = WorkGroupConfig::cpu_default();
        assert_eq!(cfg.preferred_size, 64);
        assert!(cfg.min_size >= 1);
    }
    #[test]
    fn test_parallel_reduce_sum() {
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        assert!((parallel_reduce_sum(&data) - 15.0).abs() < 1e-10);
    }
    #[test]
    fn test_parallel_reduce_sum_empty() {
        assert!((parallel_reduce_sum(&[]) - 0.0).abs() < 1e-10);
    }
    #[test]
    fn test_parallel_reduce_max() {
        let data = vec![3.0, 1.0, 4.0, 1.0, 5.0, 9.0, 2.0, 6.0];
        assert!((parallel_reduce_max(&data) - 9.0).abs() < 1e-10);
    }
    #[test]
    fn test_parallel_reduce_min() {
        let data = vec![3.0, 1.0, 4.0, 1.0, 5.0, 9.0, 2.0, 6.0];
        assert!((parallel_reduce_min(&data) - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_parallel_dot_product() {
        let a = vec![1.0, 2.0, 3.0];
        let b = vec![4.0, 5.0, 6.0];
        assert!((parallel_dot_product(&a, &b) - 32.0).abs() < 1e-10);
    }
    #[test]
    fn test_parallel_norm2() {
        let data = vec![3.0, 4.0];
        assert!((parallel_norm2(&data) - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_parallel_mean() {
        let data = vec![2.0, 4.0, 6.0, 8.0];
        assert!((parallel_mean(&data) - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_parallel_mean_empty() {
        assert!((parallel_mean(&[]) - 0.0).abs() < 1e-10);
    }
    #[test]
    fn test_parallel_variance() {
        let data = vec![2.0, 4.0, 6.0, 8.0];
        assert!((parallel_variance(&data) - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_parallel_sum_count() {
        let data = vec![10.0, 20.0, 30.0];
        let (sum, count) = parallel_sum_count(&data);
        assert!((sum - 60.0).abs() < 1e-10);
        assert_eq!(count, 3);
    }
    #[test]
    fn test_parallel_reduce_custom_product() {
        let data = vec![2.0, 3.0, 4.0];
        let product = parallel_reduce_custom(&data, 1.0, |a, b| a * b);
        assert!((product - 24.0).abs() < 1e-10);
    }
    #[test]
    fn test_parallel_exclusive_scan() {
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let result = parallel_exclusive_scan(&data);
        assert_eq!(result.len(), 5);
        let expected = [0.0, 1.0, 3.0, 6.0, 10.0];
        for i in 0..5 {
            assert!(
                (result[i] - expected[i]).abs() < 1e-10,
                "exclusive_scan[{i}]: expected {}, got {}",
                expected[i],
                result[i]
            );
        }
    }
    #[test]
    fn test_parallel_exclusive_scan_empty() {
        let result = parallel_exclusive_scan(&[]);
        assert!(result.is_empty());
    }
    #[test]
    fn test_parallel_inclusive_scan() {
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let result = parallel_inclusive_scan(&data);
        let expected = [1.0, 3.0, 6.0, 10.0, 15.0];
        for i in 0..5 {
            assert!(
                (result[i] - expected[i]).abs() < 1e-10,
                "inclusive_scan[{i}]: expected {}, got {}",
                expected[i],
                result[i]
            );
        }
    }
    #[test]
    fn test_parallel_inclusive_scan_large() {
        let n = 1000;
        let data: Vec<f64> = (0..n).map(|i| i as f64).collect();
        let result = parallel_inclusive_scan(&data);
        assert_eq!(result.len(), n);
        let expected_last = (n - 1) as f64 * n as f64 / 2.0;
        assert!(
            (result[n - 1] - expected_last).abs() < 1e-6,
            "last element: expected {expected_last}, got {}",
            result[n - 1]
        );
    }
    #[test]
    fn test_parallel_exclusive_scan_large() {
        let n = 1000;
        let data: Vec<f64> = (0..n).map(|i| i as f64).collect();
        let result = parallel_exclusive_scan(&data);
        assert_eq!(result.len(), n);
        assert!((result[0]).abs() < 1e-10);
        let expected = (n - 1) as f64 * (n - 2) as f64 / 2.0;
        assert!(
            (result[n - 1] - expected).abs() < 1e-6,
            "result[{n}]: expected {expected}, got {}",
            result[n - 1]
        );
    }
    #[test]
    fn test_segmented_exclusive_scan() {
        let data = vec![1.0, 2.0, 3.0, 10.0, 20.0];
        let segs = vec![0, 0, 0, 1, 1];
        let result = segmented_exclusive_scan(&data, &segs);
        let expected = [0.0, 1.0, 3.0, 0.0, 10.0];
        for i in 0..5 {
            assert!(
                (result[i] - expected[i]).abs() < 1e-10,
                "seg_scan[{i}]: expected {}, got {}",
                expected[i],
                result[i]
            );
        }
    }
    #[test]
    fn test_segmented_exclusive_scan_empty() {
        let result = segmented_exclusive_scan(&[], &[]);
        assert!(result.is_empty());
    }
    #[test]
    fn test_parallel_sort_f64() {
        let mut data = vec![5.0, 2.0, 8.0, 1.0, 9.0, 3.0];
        parallel_sort_f64(&mut data);
        assert_eq!(data, vec![1.0, 2.0, 3.0, 5.0, 8.0, 9.0]);
    }
    #[test]
    fn test_parallel_sort_f64_empty() {
        let mut data: Vec<f64> = vec![];
        parallel_sort_f64(&mut data);
        assert!(data.is_empty());
    }
    #[test]
    fn test_parallel_sort_f64_single() {
        let mut data = vec![42.0];
        parallel_sort_f64(&mut data);
        assert_eq!(data, vec![42.0]);
    }
    #[test]
    fn test_parallel_argsort() {
        let data = vec![3.0, 1.0, 4.0, 1.0, 5.0];
        let indices = parallel_argsort(&data);
        assert_eq!(indices.len(), 5);
        for i in 1..indices.len() {
            assert!(data[indices[i]] >= data[indices[i - 1]]);
        }
    }
    #[test]
    fn test_parallel_sort_by_key() {
        let mut items = vec![(3, "c"), (1, "a"), (2, "b")];
        parallel_sort_by_key(&mut items, |item| item.0 as f64);
        assert_eq!(items, vec![(1, "a"), (2, "b"), (3, "c")]);
    }
    #[test]
    fn test_parallel_partition() {
        let data = vec![1, 2, 3, 4, 5, 6, 7, 8];
        let (evens, odds) = parallel_partition(&data, |&x| x % 2 == 0);
        assert_eq!(evens.len(), 4);
        assert_eq!(odds.len(), 4);
        for v in &evens {
            assert_eq!(v % 2, 0);
        }
        for v in &odds {
            assert_eq!(v % 2, 1);
        }
    }
    #[test]
    fn test_parallel_rank() {
        let data = vec![30.0, 10.0, 20.0];
        let ranks = parallel_rank(&data);
        assert_eq!(ranks[0], 2);
        assert_eq!(ranks[1], 0);
        assert_eq!(ranks[2], 1);
    }
    #[test]
    fn test_load_balance_static() {
        let plan = compute_load_balance(100, 4, LoadBalanceStrategy::Static, None);
        assert_eq!(plan.num_workers(), 4);
        let mut covered = [false; 100];
        for range in &plan.ranges {
            for idx in range.clone() {
                covered[idx] = true;
            }
        }
        assert!(covered.iter().all(|&c| c));
    }
    #[test]
    fn test_load_balance_static_imbalance() {
        let plan = compute_load_balance(100, 4, LoadBalanceStrategy::Static, None);
        let ratio = plan.imbalance_ratio();
        assert!(ratio < 1.5, "imbalance too high: {ratio}");
    }
    #[test]
    fn test_load_balance_weighted() {
        let mut weights = vec![1.0; 100];
        for w in weights.iter_mut().take(50) {
            *w = 10.0;
        }
        let plan = compute_load_balance(100, 4, LoadBalanceStrategy::Weighted, Some(&weights));
        assert!(plan.num_workers() >= 1);
        let mut covered = [false; 100];
        for range in &plan.ranges {
            for idx in range.clone() {
                covered[idx] = true;
            }
        }
        assert!(covered.iter().all(|&c| c));
    }
    #[test]
    fn test_load_balance_weighted_better_than_static() {
        let mut weights = vec![1.0; 100];
        for w in weights.iter_mut().take(10) {
            *w = 50.0;
        }
        let static_plan = compute_load_balance(100, 4, LoadBalanceStrategy::Static, Some(&weights));
        let weighted_plan =
            compute_load_balance(100, 4, LoadBalanceStrategy::Weighted, Some(&weights));
        assert!(
            weighted_plan.imbalance_ratio() <= static_plan.imbalance_ratio() + 0.1,
            "weighted imbalance {} should be <= static imbalance {}",
            weighted_plan.imbalance_ratio(),
            static_plan.imbalance_ratio()
        );
    }
    #[test]
    fn test_load_balance_guided() {
        let plan = compute_load_balance(100, 4, LoadBalanceStrategy::Guided, None);
        assert!(plan.num_workers() >= 1);
        let mut covered = [false; 100];
        for range in &plan.ranges {
            for idx in range.clone() {
                covered[idx] = true;
            }
        }
        assert!(covered.iter().all(|&c| c));
    }
    #[test]
    fn test_load_balance_guided_chunk_sizes_decrease() {
        let plan = compute_load_balance(1000, 4, LoadBalanceStrategy::Guided, None);
        if plan.ranges.len() >= 2 {
            let first_len = plan.ranges[0].len();
            let last_len = plan.ranges.last().unwrap().len();
            assert!(
                first_len >= last_len,
                "guided chunks should decrease: first={first_len}, last={last_len}"
            );
        }
    }
    #[test]
    fn test_load_balance_single_worker() {
        let plan = compute_load_balance(50, 1, LoadBalanceStrategy::Static, None);
        assert_eq!(plan.num_workers(), 1);
        assert_eq!(plan.ranges[0], 0..50);
    }
    #[test]
    fn test_load_balance_zero_items() {
        let plan = compute_load_balance(0, 4, LoadBalanceStrategy::Static, None);
        assert!(plan.ranges.is_empty() || plan.ranges.iter().all(|r| r.is_empty()));
    }
    #[test]
    fn test_execute_balanced() {
        let plan = compute_load_balance(100, 4, LoadBalanceStrategy::Static, None);
        let counter = AtomicUsize::new(0);
        execute_balanced(&plan, |_worker, range| {
            counter.fetch_add(range.len(), Ordering::Relaxed);
        });
        assert_eq!(counter.load(Ordering::Relaxed), 100);
    }
    #[test]
    fn test_parallel_map_reduce() {
        let data = vec![1.0, 2.0, 3.0, 4.0];
        let result = parallel_map_reduce(&data, |&x| x * x, 0.0, |a, b| a + b);
        assert!((result - 30.0).abs() < 1e-10);
    }
    #[test]
    fn test_parallel_map_reduce_max() {
        let data = vec![1.0, 5.0, 3.0, 2.0];
        let result = parallel_map_reduce(&data, |&x| x, f64::NEG_INFINITY, f64::max);
        assert!((result - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_parallel_histogram() {
        let data = vec![0.5, 1.5, 2.5, 3.5, 0.1, 1.9, 2.1, 3.9];
        let hist = parallel_histogram(&data, 0.0, 4.0, 4);
        assert_eq!(hist[0], 2);
        assert_eq!(hist[1], 2);
        assert_eq!(hist[2], 2);
        assert_eq!(hist[3], 2);
    }
    #[test]
    fn test_parallel_histogram_empty() {
        let hist = parallel_histogram(&[], 0.0, 10.0, 5);
        assert_eq!(hist, vec![0; 5]);
    }
    #[test]
    fn test_parallel_histogram_boundary() {
        let data = vec![0.0, 10.0];
        let hist = parallel_histogram(&data, 0.0, 10.0, 2);
        assert_eq!(hist[0], 1);
        assert_eq!(hist[1], 1);
    }
    #[test]
    fn test_dist3() {
        let a = [0.0, 0.0, 0.0];
        let b = [3.0, 4.0, 0.0];
        assert!((dist3(a, b) - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_dist3_same_point() {
        let a = [1.0, 2.0, 3.0];
        assert!(dist3(a, a) < 1e-15);
    }
    #[test]
    fn test_stream_compaction_filter_positive() {
        let data = vec![-1.0f64, 2.0, -3.0, 4.0, 0.0, 5.0];
        let (compacted, scatter_map) = stream_compaction(&data, |&v| v > 0.0);
        assert_eq!(compacted, vec![2.0, 4.0, 5.0]);
        assert_eq!(scatter_map, vec![1, 3, 5]);
    }
    #[test]
    fn test_stream_compaction_empty_input() {
        let data: Vec<f64> = Vec::new();
        let (compacted, scatter_map) = stream_compaction(&data, |&v| v > 0.0);
        assert!(compacted.is_empty());
        assert!(scatter_map.is_empty());
    }
    #[test]
    fn test_stream_compaction_all_pass() {
        let data = vec![1.0f64, 2.0, 3.0];
        let (compacted, scatter_map) = stream_compaction(&data, |_| true);
        assert_eq!(compacted, data);
        assert_eq!(scatter_map, vec![0, 1, 2]);
    }
    #[test]
    fn test_stream_compaction_none_pass() {
        let data = vec![1.0f64, 2.0, 3.0];
        let (compacted, scatter_map) = stream_compaction(&data, |_| false);
        assert!(compacted.is_empty());
        assert!(scatter_map.is_empty());
    }
    #[test]
    fn test_parallel_stream_compaction_matches_serial() {
        let data = vec![3.0f64, -1.0, 5.0, -2.0, 7.0];
        let (c_serial, s_serial) = stream_compaction(&data, |&v| v > 0.0);
        let (c_par, s_par) = parallel_stream_compaction(&data, |&v| v > 0.0);
        assert_eq!(c_serial, c_par);
        assert_eq!(s_serial, s_par);
    }
    #[test]
    fn test_segmented_reduce_sum_basic() {
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let segs = vec![0usize, 0, 1, 1, 1, 2];
        let result = segmented_reduce_sum(&data, &segs);
        assert_eq!(result.len(), 3);
        assert!((result[0] - 3.0).abs() < 1e-12, "seg 0: {}", result[0]);
        assert!((result[1] - 12.0).abs() < 1e-12, "seg 1: {}", result[1]);
        assert!((result[2] - 6.0).abs() < 1e-12, "seg 2: {}", result[2]);
    }
    #[test]
    fn test_segmented_reduce_sum_empty() {
        let result = segmented_reduce_sum(&[], &[]);
        assert!(result.is_empty());
    }
    #[test]
    fn test_segmented_reduce_max_basic() {
        let data = vec![1.0, 5.0, 2.0, 4.0, 3.0];
        let segs = vec![0usize, 0, 1, 1, 1];
        let result = segmented_reduce_max(&data, &segs);
        assert!((result[0] - 5.0).abs() < 1e-12);
        assert!((result[1] - 4.0).abs() < 1e-12);
    }
    #[test]
    fn test_segmented_reduce_min_basic() {
        let data = vec![1.0, 5.0, 2.0, 4.0, 3.0];
        let segs = vec![0usize, 0, 1, 1, 1];
        let result = segmented_reduce_min(&data, &segs);
        assert!((result[0] - 1.0).abs() < 1e-12);
        assert!((result[1] - 2.0).abs() < 1e-12);
    }
    #[test]
    fn test_merge_sort_f64_basic() {
        let data = vec![5.0, 3.0, 8.0, 1.0, 9.0, 2.0];
        let sorted = merge_sort_f64(&data);
        let mut expected = data.clone();
        expected.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert_eq!(sorted, expected);
    }
    #[test]
    fn test_merge_sort_f64_empty() {
        assert!(merge_sort_f64(&[]).is_empty());
    }
    #[test]
    fn test_merge_sort_f64_single() {
        assert_eq!(merge_sort_f64(&[42.0]), vec![42.0]);
    }
    #[test]
    fn test_merge_sort_f64_already_sorted() {
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let sorted = merge_sort_f64(&data);
        assert_eq!(sorted, data);
    }
    #[test]
    fn test_merge_sort_argsort_is_permutation() {
        let data = vec![5.0, 3.0, 8.0, 1.0, 9.0];
        let indices = merge_sort_argsort(&data);
        assert_eq!(indices.len(), data.len());
        let mut check = indices.clone();
        check.sort_unstable();
        assert_eq!(check, vec![0, 1, 2, 3, 4]);
        for w in indices.windows(2) {
            assert!(data[w[0]] <= data[w[1]]);
        }
    }
    #[test]
    fn test_merge_sort_does_not_modify_input() {
        let data = vec![3.0, 1.0, 4.0, 1.5, 9.0];
        let original = data.clone();
        let _ = merge_sort_f64(&data);
        assert_eq!(data, original, "input should not be modified");
    }
    #[test]
    fn test_bitonic_sort_power_of_two() {
        let data = vec![8.0, 3.0, 6.0, 1.0, 7.0, 2.0, 5.0, 4.0];
        let sorted = bitonic_sort(&data);
        let mut expected = data.clone();
        expected.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert_eq!(sorted, expected);
    }
    #[test]
    fn test_bitonic_sort_non_power_of_two() {
        let data = vec![5.0, 3.0, 8.0, 1.0, 9.0];
        let sorted = bitonic_sort(&data);
        assert_eq!(sorted.len(), data.len());
        for w in sorted.windows(2) {
            assert!(w[0] <= w[1], "not sorted: {} > {}", w[0], w[1]);
        }
    }
    #[test]
    fn test_bitonic_sort_empty() {
        assert!(bitonic_sort(&[]).is_empty());
    }
    #[test]
    fn test_bitonic_sort_single() {
        assert_eq!(bitonic_sort(&[42.0]), vec![42.0]);
    }
    #[test]
    fn test_bitonic_sort_already_sorted() {
        let data = vec![1.0, 2.0, 3.0, 4.0];
        let sorted = bitonic_sort(&data);
        assert_eq!(sorted, data);
    }
    #[test]
    fn test_bitonic_argsort_is_permutation() {
        let data = vec![5.0, 3.0, 8.0, 1.0, 9.0, 2.0, 4.0, 7.0];
        let indices = bitonic_argsort(&data);
        assert_eq!(indices.len(), data.len());
        let mut check = indices.clone();
        check.sort_unstable();
        assert_eq!(check, vec![0, 1, 2, 3, 4, 5, 6, 7]);
        for w in indices.windows(2) {
            assert!(data[w[0]] <= data[w[1]], "not sorted via indices");
        }
    }
    #[test]
    fn test_bitonic_sort_matches_merge_sort() {
        let data = vec![9.0, 7.0, 5.0, 3.0, 1.0, 2.0, 4.0, 6.0];
        let bitonic_result = bitonic_sort(&data);
        let merge_result = merge_sort_f64(&data);
        assert_eq!(bitonic_result, merge_result);
    }
    #[test]
    fn test_work_steal_queue_push_pop() {
        let mut q: WorkStealQueue<usize> = WorkStealQueue::new();
        q.push(1);
        q.push(2);
        q.push(3);
        assert_eq!(q.pop(), Some(3));
        assert_eq!(q.pop(), Some(2));
        assert_eq!(q.pop(), Some(1));
        assert_eq!(q.pop(), None);
    }
    #[test]
    fn test_work_steal_queue_steal_from_front() {
        let mut q: WorkStealQueue<usize> = WorkStealQueue::new();
        q.push(10);
        q.push(20);
        assert_eq!(q.steal(), Some(10));
        assert_eq!(q.steal(), Some(20));
        assert_eq!(q.steal(), None);
    }
    #[test]
    fn test_work_steal_queue_is_empty() {
        let mut q: WorkStealQueue<i32> = WorkStealQueue::new();
        assert!(q.is_empty());
        q.push(5);
        assert!(!q.is_empty());
        let _ = q.pop();
        assert!(q.is_empty());
    }
    #[test]
    fn test_load_balance_metric_perfect() {
        let loads = vec![10usize, 10, 10, 10];
        let metric = compute_load_balance_metric(&loads);
        assert!(
            (metric - 1.0).abs() < 1e-12,
            "perfect balance should give 1.0"
        );
    }
    #[test]
    fn test_load_balance_metric_imbalanced() {
        let loads = vec![1usize, 1, 1, 100];
        let metric = compute_load_balance_metric(&loads);
        assert!(
            metric < 0.5,
            "heavily imbalanced should give < 0.5, got {metric}"
        );
    }
    #[test]
    fn test_load_balance_metric_empty() {
        let metric = compute_load_balance_metric(&[]);
        assert!((metric - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_merge_sort_parallel_basic() {
        let data = vec![5.0, 1.0, 4.0, 2.0, 8.0, 3.0];
        let sorted = merge_sort_parallel(&data);
        let mut expected = data.clone();
        expected.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert_eq!(sorted, expected);
    }
    #[test]
    fn test_merge_sort_parallel_empty() {
        assert!(merge_sort_parallel(&[]).is_empty());
    }
    #[test]
    fn test_merge_sort_parallel_large() {
        let data: Vec<f64> = (0..512).map(|i| (512 - i) as f64).collect();
        let sorted = merge_sort_parallel(&data);
        for w in sorted.windows(2) {
            assert!(w[0] <= w[1], "not sorted: {} > {}", w[0], w[1]);
        }
        assert_eq!(sorted.len(), 512);
    }
    #[test]
    fn test_merge_two_sorted_basic() {
        let a = [1.0, 3.0, 5.0];
        let b = [2.0, 4.0, 6.0];
        let merged = merge_two_sorted(&a, &b);
        assert_eq!(merged, vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
    }
}
