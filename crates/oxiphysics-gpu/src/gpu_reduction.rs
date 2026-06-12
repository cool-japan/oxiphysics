// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! GPU reduction and scan operations — CPU reference implementation (parity oracle for [`crate::gpu_primitives`]).
//!
//! Provides parallel-style reductions (sum, max, min, dot product),
//! exclusive prefix scan, radix sort (counting sort mock), and histogram.

/// Compute the sum of all elements in a slice (parallel reduction mock).
///
/// Returns `0.0` for an empty slice.
pub fn gpu_sum(data: &[f64]) -> f64 {
    data.iter().copied().sum()
}

/// Compute the maximum element in a slice (parallel reduction mock).
///
/// Returns `f64::NEG_INFINITY` for an empty slice.
pub fn gpu_max(data: &[f64]) -> f64 {
    data.iter().copied().fold(f64::NEG_INFINITY, f64::max)
}

/// Compute the minimum element in a slice (parallel reduction mock).
///
/// Returns `f64::INFINITY` for an empty slice.
pub fn gpu_min(data: &[f64]) -> f64 {
    data.iter().copied().fold(f64::INFINITY, f64::min)
}

/// Compute the dot product of two equal-length slices (parallel mock).
///
/// Panics in debug mode if the slices differ in length.
pub fn gpu_dot(a: &[f64], b: &[f64]) -> f64 {
    debug_assert_eq!(a.len(), b.len(), "gpu_dot: lengths must match");
    a.iter().zip(b.iter()).map(|(&x, &y)| x * y).sum()
}

/// Compute the exclusive prefix scan (prefix sum) of a slice.
///
/// `result[i] = sum(data[0..i])`. The output length equals the input length.
/// Returns an empty `Vec` for an empty input.
pub fn gpu_prefix_sum(data: &[f64]) -> Vec<f64> {
    let mut result = Vec::with_capacity(data.len());
    let mut acc = 0.0_f64;
    for &v in data {
        result.push(acc);
        acc += v;
    }
    result
}

/// Sort a slice of `u64` values using a radix sort (counting sort) mock.
///
/// This is a stable sort over the full 64-bit key space, implemented here
/// as a mock using counting-sort passes over 8 bits at a time (8 passes).
/// Returns a new sorted `Vec`.
pub fn gpu_sort_radix(data: &[u64]) -> Vec<u64> {
    if data.is_empty() {
        return Vec::new();
    }

    let mut src = data.to_vec();
    let mut dst = vec![0u64; data.len()];

    // 8 passes × 8 bits = 64 bits
    for pass in 0..8u64 {
        let shift = pass * 8;
        let mut counts = [0usize; 256];

        for &v in src.iter() {
            let bucket = ((v >> shift) & 0xFF) as usize;
            counts[bucket] += 1;
        }

        // Prefix sum to get output positions
        let mut offsets = [0usize; 256];
        let mut total = 0usize;
        for (i, &c) in counts.iter().enumerate() {
            offsets[i] = total;
            total += c;
        }

        // Scatter
        for &v in src.iter() {
            let bucket = ((v >> shift) & 0xFF) as usize;
            dst[offsets[bucket]] = v;
            offsets[bucket] += 1;
        }

        std::mem::swap(&mut src, &mut dst);
    }

    src
}

/// Compute a histogram of the input data over `num_bins` equal-width bins
/// spanning `[min_val, max_val)`.
///
/// Returns a `Vec`usize` of length `num_bins`. Values outside the range are
/// clamped to the first/last bin.
pub fn gpu_histogram(data: &[f64], min_val: f64, max_val: f64, num_bins: usize) -> Vec<usize> {
    let mut bins = vec![0usize; num_bins];
    if num_bins == 0 {
        return bins;
    }
    let range = max_val - min_val;
    for &v in data {
        let t = if range <= 0.0 {
            0.0
        } else {
            (v - min_val) / range
        };
        let idx = (t * num_bins as f64) as isize;
        let idx = idx.clamp(0, num_bins as isize - 1) as usize;
        bins[idx] += 1;
    }
    bins
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // gpu_sum
    #[test]
    fn test_gpu_sum_basic() {
        assert!((gpu_sum(&[1.0, 2.0, 3.0]) - 6.0).abs() < 1e-12);
    }

    #[test]
    fn test_gpu_sum_empty() {
        assert!((gpu_sum(&[])).abs() < 1e-12);
    }

    #[test]
    fn test_gpu_sum_single() {
        assert!((gpu_sum(&[42.0]) - 42.0).abs() < 1e-12);
    }

    #[test]
    fn test_gpu_sum_negatives() {
        assert!((gpu_sum(&[-1.0, -2.0, 3.0]) - 0.0).abs() < 1e-12);
    }

    #[test]
    fn test_gpu_sum_large() {
        let data: Vec<f64> = (0..1000).map(|i| i as f64).collect();
        let expected = 999.0 * 1000.0 / 2.0;
        assert!((gpu_sum(&data) - expected).abs() < 1e-6);
    }

    // gpu_max
    #[test]
    fn test_gpu_max_basic() {
        assert!((gpu_max(&[1.0, 5.0, 3.0]) - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_gpu_max_empty() {
        assert!(gpu_max(&[]).is_infinite() && gpu_max(&[]) < 0.0);
    }

    #[test]
    fn test_gpu_max_all_negative() {
        assert!((gpu_max(&[-5.0, -3.0, -1.0]) - (-1.0)).abs() < 1e-12);
    }

    #[test]
    fn test_gpu_max_single() {
        assert!((gpu_max(&[7.0]) - 7.0).abs() < 1e-12);
    }

    // gpu_min
    #[test]
    fn test_gpu_min_basic() {
        assert!((gpu_min(&[1.0, 5.0, 3.0]) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_gpu_min_empty() {
        assert!(gpu_min(&[]).is_infinite() && gpu_min(&[]) > 0.0);
    }

    #[test]
    fn test_gpu_min_all_positive() {
        assert!((gpu_min(&[2.0, 4.0, 6.0]) - 2.0).abs() < 1e-12);
    }

    // gpu_dot
    #[test]
    fn test_gpu_dot_basic() {
        let a = [1.0, 2.0, 3.0];
        let b = [4.0, 5.0, 6.0];
        assert!((gpu_dot(&a, &b) - 32.0).abs() < 1e-12);
    }

    #[test]
    fn test_gpu_dot_orthogonal() {
        let a = [1.0, 0.0, 0.0];
        let b = [0.0, 1.0, 0.0];
        assert!(gpu_dot(&a, &b).abs() < 1e-12);
    }

    #[test]
    fn test_gpu_dot_empty() {
        assert!(gpu_dot(&[], &[]).abs() < 1e-12);
    }

    #[test]
    fn test_gpu_dot_self() {
        let a = [3.0, 4.0];
        assert!((gpu_dot(&a, &a) - 25.0).abs() < 1e-12);
    }

    // gpu_prefix_sum
    #[test]
    fn test_gpu_prefix_sum_basic() {
        let result = gpu_prefix_sum(&[1.0, 2.0, 3.0, 4.0]);
        assert_eq!(result, vec![0.0, 1.0, 3.0, 6.0]);
    }

    #[test]
    fn test_gpu_prefix_sum_empty() {
        assert!(gpu_prefix_sum(&[]).is_empty());
    }

    #[test]
    fn test_gpu_prefix_sum_single() {
        let result = gpu_prefix_sum(&[5.0]);
        assert_eq!(result, vec![0.0]);
    }

    #[test]
    fn test_gpu_prefix_sum_all_ones() {
        let data = vec![1.0; 5];
        let result = gpu_prefix_sum(&data);
        assert_eq!(result, vec![0.0, 1.0, 2.0, 3.0, 4.0]);
    }

    #[test]
    fn test_gpu_prefix_sum_length_preserved() {
        let data = vec![2.0; 10];
        let result = gpu_prefix_sum(&data);
        assert_eq!(result.len(), 10);
    }

    // gpu_sort_radix
    #[test]
    fn test_gpu_sort_radix_basic() {
        let data = vec![5u64, 3, 8, 1, 9, 2];
        let sorted = gpu_sort_radix(&data);
        assert_eq!(sorted, vec![1, 2, 3, 5, 8, 9]);
    }

    #[test]
    fn test_gpu_sort_radix_empty() {
        assert!(gpu_sort_radix(&[]).is_empty());
    }

    #[test]
    fn test_gpu_sort_radix_single() {
        assert_eq!(gpu_sort_radix(&[42]), vec![42]);
    }

    #[test]
    fn test_gpu_sort_radix_already_sorted() {
        let data = vec![1u64, 2, 3, 4, 5];
        assert_eq!(gpu_sort_radix(&data), data);
    }

    #[test]
    fn test_gpu_sort_radix_reverse_sorted() {
        let data = vec![5u64, 4, 3, 2, 1];
        let sorted = gpu_sort_radix(&data);
        assert_eq!(sorted, vec![1, 2, 3, 4, 5]);
    }

    #[test]
    fn test_gpu_sort_radix_duplicates() {
        let data = vec![3u64, 1, 3, 2, 1];
        let sorted = gpu_sort_radix(&data);
        assert_eq!(sorted, vec![1, 1, 2, 3, 3]);
    }

    #[test]
    fn test_gpu_sort_radix_large_values() {
        let data = vec![u64::MAX, 0, u64::MAX / 2, 1];
        let sorted = gpu_sort_radix(&data);
        assert_eq!(sorted[0], 0);
        assert_eq!(sorted[3], u64::MAX);
    }

    // gpu_histogram
    #[test]
    fn test_gpu_histogram_basic() {
        let data = vec![0.5, 1.5, 2.5, 3.5];
        let hist = gpu_histogram(&data, 0.0, 4.0, 4);
        assert_eq!(hist, vec![1, 1, 1, 1]);
    }

    #[test]
    fn test_gpu_histogram_empty_data() {
        let hist = gpu_histogram(&[], 0.0, 10.0, 5);
        assert_eq!(hist, vec![0, 0, 0, 0, 0]);
    }

    #[test]
    fn test_gpu_histogram_zero_bins() {
        let hist = gpu_histogram(&[1.0, 2.0], 0.0, 3.0, 0);
        assert!(hist.is_empty());
    }

    #[test]
    fn test_gpu_histogram_all_in_one_bin() {
        let data = vec![1.0, 1.1, 1.2];
        let hist = gpu_histogram(&data, 0.0, 10.0, 10);
        assert_eq!(hist[1], 3);
    }

    #[test]
    fn test_gpu_histogram_out_of_range_clamped() {
        let data = vec![-100.0, 200.0];
        let hist = gpu_histogram(&data, 0.0, 10.0, 5);
        // Both values are clamped: -100 → bin 0, 200 → bin 4
        assert_eq!(hist[0], 1);
        assert_eq!(hist[4], 1);
    }

    #[test]
    fn test_gpu_histogram_bin_count() {
        let data: Vec<f64> = (0..100).map(|i| i as f64).collect();
        let hist = gpu_histogram(&data, 0.0, 100.0, 10);
        assert_eq!(hist.len(), 10);
        // Each bin should have roughly 10 elements
        for &count in &hist {
            assert!((9..=11).contains(&count));
        }
    }

    // Combined / integration tests
    #[test]
    fn test_prefix_sum_then_last_is_total() {
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let ps = gpu_prefix_sum(&data);
        // The total sum should equal last exclusive prefix + last element
        let total = ps.last().unwrap() + data.last().unwrap();
        assert!((total - gpu_sum(&data)).abs() < 1e-12);
    }

    #[test]
    fn test_sort_then_min_max() {
        let data = vec![9u64, 3, 7, 1, 5];
        let sorted = gpu_sort_radix(&data);
        assert_eq!(sorted[0], 1);
        assert_eq!(*sorted.last().unwrap(), 9);
    }
}
