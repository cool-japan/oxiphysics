//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
#[cfg(test)]
mod tests {
    use super::*;

    use crate::parallel::ParallelFor;
    use crate::parallel::SerialWorkQueue;

    use crate::parallel::ThreadPoolStats;
    use crate::parallel::WorkRange;

    use crate::parallel::WorkStealingPool;
    use crate::parallel::parallel_dot_product;

    use crate::parallel::prefix_sum;

    use std::sync::Arc;
    use std::sync::atomic::{AtomicUsize, Ordering};
    #[test]
    fn test_parallel_map_f64_square() {
        let data = vec![1.0_f64, 2.0, 3.0, 4.0];
        let result = parallel_map_f64(&data, |x| x * x);
        assert_eq!(result, vec![1.0, 4.0, 9.0, 16.0]);
    }
    #[test]
    fn test_parallel_reduce_f64_sum() {
        let data = vec![1.0_f64, 2.0, 3.0, 4.0, 5.0];
        let result = parallel_reduce_f64(&data, 0.0, |a, b| a + b);
        assert!((result - 15.0).abs() < f64::EPSILON);
    }
    #[test]
    fn test_scatter_gather_roundtrip() {
        let data: Vec<i32> = (0..10).collect();
        let parts = scatter_gather(&data, 3);
        assert_eq!(parts.len(), 3);
        assert_eq!(parts[0].len(), 4);
        assert_eq!(parts[1].len(), 3);
        assert_eq!(parts[2].len(), 3);
        let reassembled = gather(parts);
        assert_eq!(reassembled, data);
    }
    #[test]
    fn test_parallel_dot_product() {
        let a = [1.0_f64, 2.0, 3.0];
        let b = [4.0_f64, 5.0, 6.0];
        assert!((parallel_dot_product(&a, &b) - 32.0).abs() < f64::EPSILON);
    }
    #[test]
    fn test_parallel_matrix_vec_multiply_identity() {
        let matrix = vec![vec![1.0_f64, 0.0], vec![0.0_f64, 1.0]];
        let vector = [3.0_f64, 4.0];
        let result = parallel_matrix_vec_multiply(&matrix, &vector);
        assert!((result[0] - 3.0).abs() < f64::EPSILON);
        assert!((result[1] - 4.0).abs() < f64::EPSILON);
    }
    #[test]
    fn test_work_range_chunks() {
        let wr = WorkRange::new(0, 10, 3);
        let chunks = wr.chunks();
        assert_eq!(chunks, vec![(0, 3), (3, 6), (6, 9), (9, 10)]);
        assert_eq!(wr.n_chunks(), 4);
        assert_eq!(wr.total_work(), 10);
    }
    #[test]
    fn test_parallel_map_square() {
        let data: Vec<i32> = (1..=10).collect();
        let result = parallel_map(&data, |x| x * x);
        let expected: Vec<i32> = (1..=10).map(|x| x * x).collect();
        assert_eq!(result, expected);
    }
    #[test]
    fn test_parallel_map_preserves_order() {
        let data: Vec<usize> = (0..100).collect();
        let result = parallel_map(&data, |&x| x * 2);
        let expected: Vec<usize> = (0..100).map(|x| x * 2).collect();
        assert_eq!(result, expected, "parallel_map must preserve order");
    }
    #[test]
    fn test_parallel_map_empty() {
        let data: Vec<i32> = vec![];
        let result = parallel_map(&data, |x| x * x);
        assert!(result.is_empty());
    }
    #[test]
    fn test_parallel_filter_even() {
        let data: Vec<i32> = (0..20).collect();
        let result = parallel_filter(&data, |x| x % 2 == 0);
        let expected: Vec<i32> = (0..20).filter(|x| x % 2 == 0).collect();
        assert_eq!(result, expected);
    }
    #[test]
    fn test_parallel_filter_empty() {
        let data: Vec<i32> = vec![];
        let result = parallel_filter(&data, |_| true);
        assert!(result.is_empty());
    }
    #[test]
    fn test_parallel_filter_none_match() {
        let data: Vec<i32> = (0..10).collect();
        let result = parallel_filter(&data, |_| false);
        assert!(result.is_empty());
    }
    #[test]
    fn test_parallel_for_each_counter() {
        let counter = Arc::new(AtomicUsize::new(0));
        let data: Vec<i32> = (0..50).collect();
        let cc = Arc::clone(&counter);
        parallel_for_each(&data, move |_| {
            cc.fetch_add(1, Ordering::Relaxed);
        });
        assert_eq!(counter.load(Ordering::Relaxed), 50);
    }
    #[test]
    fn test_parallel_reduce_sum() {
        let data: Vec<i64> = (1..=100).collect();
        let result = parallel_reduce(&data, |a, b| a + b, 0_i64);
        assert_eq!(result, 5050);
    }
    #[test]
    fn test_parallel_reduce_empty() {
        let data: Vec<i32> = vec![];
        let result = parallel_reduce(&data, |a, b| a + b, 42_i32);
        assert_eq!(result, 42);
    }
    pub(super) struct SumOp;
    impl ReduceOperator for SumOp {
        type Acc = f64;
        type Item = f64;
        type Result = f64;
        fn identity(&self) -> f64 {
            0.0
        }
        fn fold(&self, acc: f64, item: f64) -> f64 {
            acc + item
        }
        fn combine(&self, left: f64, right: f64) -> f64 {
            left + right
        }
        fn finalize(&self, acc: f64) -> f64 {
            acc
        }
    }
    #[test]
    fn test_reduce_with_op() {
        let data: Vec<f64> = (1..=10).map(|x| x as f64).collect();
        let result = parallel_reduce_with_op(&data, &SumOp);
        assert!((result - 55.0).abs() < 1e-10);
    }
    pub(super) struct MeanOp;
    impl ReduceOperator for MeanOp {
        type Acc = (f64, usize);
        type Item = f64;
        type Result = f64;
        fn identity(&self) -> (f64, usize) {
            (0.0, 0)
        }
        fn fold(&self, acc: (f64, usize), item: f64) -> (f64, usize) {
            (acc.0 + item, acc.1 + 1)
        }
        fn combine(&self, left: (f64, usize), right: (f64, usize)) -> (f64, usize) {
            (left.0 + right.0, left.1 + right.1)
        }
        fn finalize(&self, acc: (f64, usize)) -> f64 {
            if acc.1 == 0 {
                0.0
            } else {
                acc.0 / acc.1 as f64
            }
        }
    }
    #[test]
    fn test_reduce_with_mean_op() {
        let data: Vec<f64> = vec![2.0, 4.0, 6.0, 8.0];
        let result = parallel_reduce_with_op(&data, &MeanOp);
        assert!((result - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_parallel_for_coverage() {
        let counter = Arc::new(AtomicUsize::new(0));
        let counter_clone = Arc::clone(&counter);
        let pf = ParallelFor::with_chunks(4);
        pf.run(0, 100, move |_i| {
            counter_clone.fetch_add(1, Ordering::Relaxed);
        });
        assert_eq!(counter.load(Ordering::Relaxed), 100);
    }
    #[test]
    fn test_parallel_for_run_n() {
        let counter = Arc::new(AtomicUsize::new(0));
        let counter_clone = Arc::clone(&counter);
        let pf = ParallelFor::with_chunks(2);
        pf.run_n(50, move |_| {
            counter_clone.fetch_add(1, Ordering::Relaxed);
        });
        assert_eq!(counter.load(Ordering::Relaxed), 50);
    }
    #[test]
    fn test_parallel_for_empty_range() {
        let counter = Arc::new(AtomicUsize::new(0));
        let cc = Arc::clone(&counter);
        let pf = ParallelFor::new();
        pf.run(5, 5, move |_| {
            cc.fetch_add(1, Ordering::Relaxed);
        });
        assert_eq!(counter.load(Ordering::Relaxed), 0);
    }
    #[test]
    fn test_work_stealing_pool_basic() {
        let pool = WorkStealingPool::new(4);
        let counter = Arc::new(AtomicUsize::new(0));
        for _ in 0..10 {
            let cc = Arc::clone(&counter);
            pool.submit(move || {
                cc.fetch_add(1, Ordering::Relaxed);
            });
        }
        pool.join();
        assert_eq!(counter.load(Ordering::Relaxed), 10);
    }
    #[test]
    fn test_work_stealing_pool_empty() {
        let pool = WorkStealingPool::new(2);
        pool.join();
    }
    #[test]
    fn test_work_stealing_pool_stats() {
        let pool = WorkStealingPool::new(2);
        for _ in 0..5 {
            pool.submit(|| {});
        }
        let stats = pool.stats();
        assert_eq!(stats.tasks_submitted, 5);
        pool.join();
    }
    #[test]
    fn test_thread_pool_stats() {
        let mut stats = ThreadPoolStats::new(4);
        stats.tasks_submitted = 100;
        stats.tasks_completed = 100;
        assert!((stats.tasks_per_worker() - 25.0).abs() < 1e-10);
        assert!((stats.completion_rate() - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_thread_pool_stats_zero_workers() {
        let stats = ThreadPoolStats::new(0);
        assert!(stats.tasks_per_worker().abs() < 1e-10);
    }
    #[test]
    fn test_parallel_merge_sort_basic() {
        let mut data = vec![5.0, 3.0, 1.0, 4.0, 2.0];
        parallel_merge_sort(&mut data);
        assert_eq!(data, vec![1.0, 2.0, 3.0, 4.0, 5.0]);
    }
    #[test]
    fn test_parallel_merge_sort_large() {
        let mut data: Vec<f64> = (0..200).rev().map(|x| x as f64).collect();
        parallel_merge_sort(&mut data);
        let expected: Vec<f64> = (0..200).map(|x| x as f64).collect();
        assert_eq!(data, expected);
    }
    #[test]
    fn test_parallel_merge_sort_empty() {
        let mut data: Vec<f64> = vec![];
        parallel_merge_sort(&mut data);
        assert!(data.is_empty());
    }
    #[test]
    fn test_parallel_merge_sort_single() {
        let mut data = vec![42.0];
        parallel_merge_sort(&mut data);
        assert_eq!(data, vec![42.0]);
    }
    #[test]
    fn test_parallel_merge_sort_already_sorted() {
        let mut data: Vec<f64> = (0..100).map(|x| x as f64).collect();
        parallel_merge_sort(&mut data);
        let expected: Vec<f64> = (0..100).map(|x| x as f64).collect();
        assert_eq!(data, expected);
    }
    #[test]
    fn test_prefix_sum() {
        let data = [1.0, 2.0, 3.0, 4.0];
        let result = prefix_sum(&data);
        assert_eq!(result, vec![1.0, 3.0, 6.0, 10.0]);
    }
    #[test]
    fn test_exclusive_prefix_sum() {
        let data = [1.0, 2.0, 3.0, 4.0];
        let result = exclusive_prefix_sum(&data);
        assert_eq!(result, vec![0.0, 1.0, 3.0, 6.0]);
    }
    #[test]
    fn test_parallel_min() {
        let data = vec![5.0, 2.0, 8.0, 1.0, 9.0];
        assert!((parallel_min(&data) - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_parallel_max() {
        let data = vec![5.0, 2.0, 8.0, 1.0, 9.0];
        assert!((parallel_max(&data) - 9.0).abs() < 1e-10);
    }
    #[test]
    fn test_parallel_min_empty() {
        let data: Vec<f64> = vec![];
        assert_eq!(parallel_min(&data), f64::INFINITY);
    }
    #[test]
    fn test_parallel_max_empty() {
        let data: Vec<f64> = vec![];
        assert_eq!(parallel_max(&data), f64::NEG_INFINITY);
    }
    #[test]
    fn test_merge_sorted() {
        let a = [1.0, 3.0, 5.0];
        let b = [2.0, 4.0, 6.0];
        let merged = merge_sorted(&a, &b);
        assert_eq!(merged, vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
    }
    #[test]
    fn test_merge_sorted_empty() {
        let a: Vec<f64> = vec![];
        let b = vec![1.0, 2.0];
        assert_eq!(merge_sorted(&a, &b), vec![1.0, 2.0]);
        assert_eq!(merge_sorted(&b, &a), vec![1.0, 2.0]);
    }
    #[test]
    fn test_insertion_sort() {
        let mut data = vec![5.0, 1.0, 3.0, 2.0, 4.0];
        insertion_sort(&mut data);
        assert_eq!(data, vec![1.0, 2.0, 3.0, 4.0, 5.0]);
    }
    #[test]
    fn test_available_threads_positive() {
        assert!(available_threads() >= 1, "should report at least 1 thread");
    }
    #[test]
    fn test_suggested_thread_count_at_least_one() {
        assert!(suggested_thread_count() >= 1);
    }
    #[test]
    fn test_chunk_process_sum_chunks() {
        let data: Vec<f64> = (1..=12).map(|x| x as f64).collect();
        let result = chunk_process(&data, 4, |chunk| vec![chunk.iter().sum::<f64>()]);
        assert_eq!(result, vec![10.0, 26.0, 42.0]);
    }
    #[test]
    fn test_chunk_process_empty() {
        let data: Vec<f64> = vec![];
        let result = chunk_process(&data, 4, |c| vec![c.iter().sum::<f64>()]);
        assert!(result.is_empty());
    }
    #[test]
    fn test_chunk_process_single_element_chunks() {
        let data = vec![1.0, 2.0, 3.0];
        let result = chunk_process(&data, 1, |c| vec![c[0] * 2.0]);
        assert_eq!(result, vec![2.0, 4.0, 6.0]);
    }
    #[test]
    fn test_chunk_zip_map_add() {
        let a = [1.0, 2.0, 3.0, 4.0];
        let b = [10.0, 20.0, 30.0, 40.0];
        let result = chunk_zip_map(&a, &b, 2, |x, y| x + y);
        assert_eq!(result, vec![11.0, 22.0, 33.0, 44.0]);
    }
    #[test]
    fn test_chunk_zip_map_multiply() {
        let a = [2.0, 3.0, 4.0];
        let b = [5.0, 6.0, 7.0];
        let result = chunk_zip_map(&a, &b, 8, |x, y| x * y);
        assert_eq!(result, vec![10.0, 18.0, 28.0]);
    }
    #[test]
    fn test_chunk_dot_product_basic() {
        let a = [1.0, 2.0, 3.0];
        let b = [4.0, 5.0, 6.0];
        let result = chunk_dot_product(&a, &b, 2);
        assert!((result - 32.0).abs() < 1e-10);
    }
    #[test]
    fn test_chunk_dot_product_chunk_size_1() {
        let a = [1.0, 2.0, 3.0, 4.0];
        let b = [1.0, 1.0, 1.0, 1.0];
        let result = chunk_dot_product(&a, &b, 1);
        assert!((result - 10.0).abs() < 1e-10);
    }
    #[test]
    fn test_chunk_dot_product_empty() {
        assert!((chunk_dot_product(&[], &[], 4)).abs() < 1e-10);
    }
    #[test]
    fn test_serial_work_queue_push_pop() {
        let mut q: SerialWorkQueue<i32> = SerialWorkQueue::new();
        q.push(1);
        q.push(2);
        q.push(3);
        assert_eq!(q.pop(), Some(1));
        assert_eq!(q.pop(), Some(2));
        assert_eq!(q.len(), 1);
    }
    #[test]
    fn test_serial_work_queue_steal() {
        let mut q: SerialWorkQueue<i32> = SerialWorkQueue::new();
        q.push(1);
        q.push(2);
        q.push(3);
        assert_eq!(q.steal(), Some(3));
        assert_eq!(q.len(), 2);
    }
    #[test]
    fn test_serial_work_queue_empty() {
        let mut q: SerialWorkQueue<i32> = SerialWorkQueue::new();
        assert!(q.is_empty());
        assert_eq!(q.pop(), None);
        assert_eq!(q.steal(), None);
    }
    #[test]
    fn test_serial_work_queue_drain_and_run() {
        let mut q: SerialWorkQueue<i32> = SerialWorkQueue::new();
        for i in 0..5 {
            q.push(i);
        }
        let mut collected = Vec::new();
        q.drain_and_run(|x| collected.push(x));
        assert_eq!(collected, vec![0, 1, 2, 3, 4]);
        assert!(q.is_empty());
    }
    #[test]
    fn test_parallel_sort_basic() {
        let mut data = vec![5.0, 3.0, 1.0, 4.0, 2.0];
        parallel_sort(&mut data);
        assert_eq!(data, vec![1.0, 2.0, 3.0, 4.0, 5.0]);
    }
    #[test]
    fn test_parallel_sort_large() {
        let mut data: Vec<f64> = (0..500).rev().map(|x| x as f64).collect();
        parallel_sort(&mut data);
        let expected: Vec<f64> = (0..500).map(|x| x as f64).collect();
        assert_eq!(data, expected);
    }
    #[test]
    fn test_parallel_sort_empty() {
        let mut data: Vec<f64> = vec![];
        parallel_sort(&mut data);
        assert!(data.is_empty());
    }
    #[test]
    fn test_parallel_sort_already_sorted() {
        let mut data: Vec<f64> = (0..100).map(|x| x as f64).collect();
        parallel_sort(&mut data);
        assert_eq!(data, (0..100).map(|x| x as f64).collect::<Vec<_>>());
    }
    #[test]
    fn test_parallel_sort_duplicates() {
        let mut data = vec![3.0, 1.0, 3.0, 2.0, 1.0];
        parallel_sort(&mut data);
        assert_eq!(data, vec![1.0, 1.0, 2.0, 3.0, 3.0]);
    }
    #[test]
    fn test_sorted_copy_does_not_mutate_input() {
        let original = vec![5.0, 3.0, 1.0, 4.0, 2.0];
        let sorted = sorted_copy(&original);
        assert_eq!(
            original,
            vec![5.0, 3.0, 1.0, 4.0, 2.0],
            "original must be unchanged"
        );
        assert_eq!(sorted, vec![1.0, 2.0, 3.0, 4.0, 5.0]);
    }
    #[test]
    fn test_parallel_histogram_uniform() {
        let data: Vec<f64> = (0..100).map(|i| i as f64).collect();
        let counts = parallel_histogram(&data, 10);
        assert_eq!(counts.len(), 10);
        let total: usize = counts.iter().sum();
        assert_eq!(total, 100);
    }
    #[test]
    fn test_parallel_histogram_empty() {
        assert!(parallel_histogram(&[], 5).is_empty());
    }
    #[test]
    fn test_parallel_histogram_zero_bins() {
        let data = vec![1.0, 2.0, 3.0];
        assert!(parallel_histogram(&data, 0).is_empty());
    }
    #[test]
    fn test_parallel_histogram_single_value() {
        let data = vec![5.0; 10];
        let counts = parallel_histogram(&data, 4);
        let total: usize = counts.iter().sum();
        assert_eq!(total, 10);
    }
}
/// Compute a chunked (SIMD-mimic) dot product of `a` and `b`.
///
/// Unrolls into chunks of 4 elements to mimic 4-wide SIMD lanes.
/// Falls back to scalar for the remaining elements.
#[allow(dead_code)]
pub fn vectorized_dot_product(a: &[f64], b: &[f64]) -> f64 {
    assert_eq!(a.len(), b.len(), "vectorized_dot_product: length mismatch");
    let n = a.len();
    let chunks = n / 4;
    let mut acc0 = 0.0_f64;
    let mut acc1 = 0.0_f64;
    let mut acc2 = 0.0_f64;
    let mut acc3 = 0.0_f64;
    for i in 0..chunks {
        let base = i * 4;
        acc0 += a[base] * b[base];
        acc1 += a[base + 1] * b[base + 1];
        acc2 += a[base + 2] * b[base + 2];
        acc3 += a[base + 3] * b[base + 3];
    }
    let mut total = acc0 + acc1 + acc2 + acc3;
    for i in (chunks * 4)..n {
        total += a[i] * b[i];
    }
    total
}
/// Parallel-style prefix sum using a two-pass approach.
///
/// Phase 1: Compute local prefix sums for each chunk in parallel.
/// Phase 2: Adjust each chunk by the cumulative sum of previous chunks.
///
/// In practice both phases run sequentially here (correctness reference).
#[allow(dead_code)]
pub fn parallel_prefix_scan(data: &[f64], n_chunks: usize) -> Vec<f64> {
    if data.is_empty() {
        return vec![];
    }
    let n = data.len();
    let n_chunks = n_chunks.max(1).min(n);
    let chunk_size = n.div_ceil(n_chunks);
    let mut local: Vec<Vec<f64>> = data
        .chunks(chunk_size)
        .map(|chunk| {
            let mut acc = 0.0;
            chunk
                .iter()
                .map(|&x| {
                    acc += x;
                    acc
                })
                .collect()
        })
        .collect();
    let mut running = 0.0_f64;
    for chunk in &mut local {
        let total = *chunk.last().unwrap_or(&0.0);
        for v in chunk.iter_mut() {
            *v += running;
        }
        running += total;
    }
    local.into_iter().flatten().collect()
}
#[cfg(test)]
mod tests_new_parallel {

    use crate::parallel::ExtendedPoolStats;

    use crate::parallel::SoaVec3;

    use crate::parallel::WorkStealingDeque;

    use crate::parallel::parallel_dot_product;
    use crate::parallel::parallel_prefix_scan;
    use crate::parallel::prefix_sum;
    use crate::parallel::vectorized_dot_product;
    use std::sync::{Arc, Mutex};
    #[test]
    fn test_soa_vec3_push_get() {
        let mut soa = SoaVec3::new();
        soa.push(1.0, 2.0, 3.0);
        soa.push(4.0, 5.0, 6.0);
        assert_eq!(soa.len(), 2);
        assert_eq!(soa.get(0), (1.0, 2.0, 3.0));
        assert_eq!(soa.get(1), (4.0, 5.0, 6.0));
    }
    #[test]
    fn test_soa_vec3_dot_with() {
        let mut soa = SoaVec3::new();
        soa.push(1.0, 0.0, 0.0);
        soa.push(0.0, 1.0, 0.0);
        soa.push(0.0, 0.0, 1.0);
        let dots = soa.dot_with(1.0, 2.0, 3.0);
        assert!((dots[0] - 1.0).abs() < 1e-10);
        assert!((dots[1] - 2.0).abs() < 1e-10);
        assert!((dots[2] - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_soa_vec3_norms_sq() {
        let mut soa = SoaVec3::new();
        soa.push(3.0, 4.0, 0.0);
        let norms = soa.norms_sq();
        assert!(
            (norms[0] - 25.0).abs() < 1e-10,
            "3-4-0 vector norm²=25, got {}",
            norms[0]
        );
    }
    #[test]
    fn test_soa_vec3_empty() {
        let soa = SoaVec3::new();
        assert!(soa.is_empty());
        assert_eq!(soa.len(), 0);
    }
    #[test]
    fn test_soa_vec3_with_capacity() {
        let soa = SoaVec3::with_capacity(100);
        assert!(soa.is_empty());
    }
    #[test]
    fn test_vectorized_dot_product_basic() {
        let a = [1.0, 2.0, 3.0, 4.0, 5.0];
        let b = [1.0, 1.0, 1.0, 1.0, 1.0];
        let result = vectorized_dot_product(&a, &b);
        assert!((result - 15.0).abs() < 1e-10);
    }
    #[test]
    fn test_vectorized_dot_product_orthogonal() {
        let a = [1.0, 0.0, 0.0, 0.0];
        let b = [0.0, 1.0, 0.0, 0.0];
        let result = vectorized_dot_product(&a, &b);
        assert!(
            result.abs() < 1e-10,
            "orthogonal vectors: dot=0, got {result}"
        );
    }
    #[test]
    fn test_vectorized_dot_product_matches_parallel_dot() {
        let a: Vec<f64> = (1..=100).map(|x| x as f64).collect();
        let b: Vec<f64> = (1..=100).map(|x| x as f64).collect();
        let v = vectorized_dot_product(&a, &b);
        let p = parallel_dot_product(&a, &b);
        assert!(
            (v - p).abs() < 1e-6,
            "vectorized and parallel dot should match"
        );
    }
    #[test]
    fn test_vectorized_dot_product_non_multiple_of_4() {
        let a = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0];
        let b = [1.0; 7];
        let result = vectorized_dot_product(&a, &b);
        assert!((result - 28.0).abs() < 1e-10);
    }
    #[test]
    fn test_parallel_prefix_scan_basic() {
        let data = [1.0, 2.0, 3.0, 4.0];
        let result = parallel_prefix_scan(&data, 2);
        assert_eq!(result, vec![1.0, 3.0, 6.0, 10.0]);
    }
    #[test]
    fn test_parallel_prefix_scan_single_chunk() {
        let data = [1.0, 2.0, 3.0, 4.0, 5.0];
        let result = parallel_prefix_scan(&data, 1);
        assert_eq!(result, vec![1.0, 3.0, 6.0, 10.0, 15.0]);
    }
    #[test]
    fn test_parallel_prefix_scan_matches_prefix_sum() {
        let data: Vec<f64> = (1..=10).map(|x| x as f64).collect();
        let scan = parallel_prefix_scan(&data, 4);
        let psum = prefix_sum(&data);
        for (a, b) in scan.iter().zip(psum.iter()) {
            assert!((a - b).abs() < 1e-10, "parallel scan {a} vs prefix_sum {b}");
        }
    }
    #[test]
    fn test_parallel_prefix_scan_empty() {
        let result = parallel_prefix_scan(&[], 4);
        assert!(result.is_empty());
    }
    #[test]
    fn test_extended_pool_stats_throughput() {
        let mut stats = ExtendedPoolStats::new(4);
        stats.base.tasks_completed = 1000;
        stats.total_ns = 1_000_000;
        let tp = stats.throughput_tasks_per_us();
        assert!(
            (tp - 1.0).abs() < 1e-10,
            "1000 tasks / 1000µs = 1 task/µs, got {tp}"
        );
    }
    #[test]
    fn test_extended_pool_stats_zero_ns() {
        let stats = ExtendedPoolStats::new(2);
        assert_eq!(stats.throughput_tasks_per_us(), 0.0);
    }
    #[test]
    fn test_extended_pool_stats_memory_efficiency() {
        let mut stats = ExtendedPoolStats::new(4);
        stats.base.tasks_completed = 100;
        stats.peak_memory_bytes = 1024;
        let eff = stats.memory_efficiency();
        assert!(
            (eff - 100.0).abs() < 1e-10,
            "100 tasks / 1 KB = 100 tasks/KB, got {eff}"
        );
    }
    #[test]
    fn test_work_stealing_deque_push_pop_bottom() {
        let mut deque: WorkStealingDeque<i32> = WorkStealingDeque::new();
        deque.push_bottom(1);
        deque.push_bottom(2);
        deque.push_bottom(3);
        assert_eq!(deque.pop_bottom(), Some(3));
        assert_eq!(deque.pop_bottom(), Some(2));
        assert_eq!(deque.len(), 1);
    }
    #[test]
    fn test_work_stealing_deque_steal_top() {
        let mut deque: WorkStealingDeque<i32> = WorkStealingDeque::new();
        for i in 0..5 {
            deque.push_bottom(i);
        }
        assert_eq!(deque.steal_top(), Some(0));
        assert_eq!(deque.steal_top(), Some(1));
        assert_eq!(deque.steal_count(), 2);
    }
    #[test]
    fn test_work_stealing_deque_empty_returns_none() {
        let mut deque: WorkStealingDeque<i32> = WorkStealingDeque::new();
        assert_eq!(deque.pop_bottom(), None);
        assert_eq!(deque.steal_top(), None);
    }
    #[test]
    fn test_work_stealing_deque_is_empty() {
        let mut deque: WorkStealingDeque<i32> = WorkStealingDeque::new();
        assert!(deque.is_empty());
        deque.push_bottom(42);
        assert!(!deque.is_empty());
    }
    #[test]
    fn test_work_stealing_deque_pop_count() {
        let mut deque: WorkStealingDeque<i32> = WorkStealingDeque::new();
        for i in 0..5 {
            deque.push_bottom(i);
        }
        let _ = deque.pop_bottom();
        let _ = deque.pop_bottom();
        assert_eq!(deque.pop_count(), 2);
    }
    #[test]
    fn test_parallel_prefix_scan_all_ones() {
        let data = vec![1.0; 8];
        let result = parallel_prefix_scan(&data, 4);
        let expected: Vec<f64> = (1..=8).map(|i| i as f64).collect();
        for (a, b) in result.iter().zip(expected.iter()) {
            assert!((a - b).abs() < 1e-10, "scan[i]={a} expected {b}");
        }
    }
    #[test]
    fn test_vectorized_dot_single_element() {
        let a = [3.0];
        let b = [4.0];
        let result = vectorized_dot_product(&a, &b);
        assert!((result - 12.0).abs() < 1e-10);
    }
    #[test]
    fn test_soa_vec3_large_push() {
        let mut soa = SoaVec3::with_capacity(1000);
        for i in 0..1000 {
            soa.push(i as f64, i as f64 * 2.0, i as f64 * 3.0);
        }
        assert_eq!(soa.len(), 1000);
        let norms = soa.norms_sq();
        for (idx, &n) in norms.iter().enumerate() {
            let expected = 14.0 * (idx as f64).powi(2);
            assert!(
                (n - expected).abs() < 1e-6 * (expected.abs() + 1.0),
                "norm²[{idx}]={n}, expected {expected}"
            );
        }
    }
    #[test]
    fn test_extended_pool_stats_zero_memory() {
        let stats = ExtendedPoolStats::new(4);
        assert_eq!(stats.memory_efficiency(), 0.0);
    }
    #[test]
    fn test_work_stealing_simulation() {
        let shared = Arc::new(Mutex::new(WorkStealingDeque::<i32>::new()));
        for i in 0..10 {
            shared
                .lock()
                .unwrap_or_else(|e| e.into_inner())
                .push_bottom(i);
        }
        let mut stolen = Vec::new();
        for _ in 0..3 {
            if let Some(v) = shared.lock().unwrap_or_else(|e| e.into_inner()).steal_top() {
                stolen.push(v);
            }
        }
        assert_eq!(stolen.len(), 3);
        assert_eq!(stolen, vec![0, 1, 2], "steals should come from top (FIFO)");
        assert_eq!(
            shared.lock().unwrap_or_else(|e| e.into_inner()).len(),
            7,
            "7 items should remain"
        );
    }
}
