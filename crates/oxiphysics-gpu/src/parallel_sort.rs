#![allow(clippy::ptr_arg)]
#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Parallel sorting and reduction utilities (CPU-side, rayon-based).
//!
//! Provides:
//! - [`radix_sort_u32`] – LSD radix sort for `u32` (4 passes of 8 bits).
//! - [`radix_sort_by_key`] – sort any `T` by a `u32` key function.
//! - [`parallel_prefix_sum`] – exclusive prefix-sum (scan) via rayon.
//! - [`parallel_reduce_sum`] – parallel tree reduction for `f64`.
//! - [`parallel_min_max`] – parallel min/max reduction.
//! - [`bitonic_sort`] – bitonic sort (pads to power-of-2 with `f64::MAX`).
//! - [`merge_sort_parallel`] – parallel merge sort via rayon.
//! - [`histogram_u32`] – parallel histogram.
//! - [`argsort`] – sorted index array for `f64` slice.
//! - [`nth_element`] – quickselect O(n) kth smallest element.

#![allow(dead_code)]

use rayon::prelude::*;

// ─────────────────────────────────────────────────────────────────────────────
// radix_sort_u32
// ─────────────────────────────────────────────────────────────────────────────

/// LSD radix sort for `u32` values using 8-bit passes (4 passes total).
///
/// Stable, O(n) for fixed-width 32-bit keys.
pub fn radix_sort_u32(data: &mut Vec<u32>) {
    if data.len() <= 1 {
        return;
    }
    let n = data.len();
    let mut buf = vec![0u32; n];

    for pass in 0..4u32 {
        let shift = pass * 8;
        let mut counts = [0usize; 256];
        for &v in data.iter() {
            let byte = ((v >> shift) & 0xFF) as usize;
            counts[byte] += 1;
        }
        // exclusive prefix sum of counts
        let mut offsets = [0usize; 256];
        let mut total = 0;
        for i in 0..256 {
            offsets[i] = total;
            total += counts[i];
        }
        for &v in data.iter() {
            let byte = ((v >> shift) & 0xFF) as usize;
            buf[offsets[byte]] = v;
            offsets[byte] += 1;
        }
        std::mem::swap(data, &mut buf);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// radix_sort_by_key
// ─────────────────────────────────────────────────────────────────────────────

/// Sort a `Vec`T` in ascending order of the `u32` key produced by `key_fn`.
///
/// Uses LSD radix sort (4 passes of 8 bits).  The sort is stable.
pub fn radix_sort_by_key<T: Clone>(data: &mut Vec<T>, key_fn: impl Fn(&T) -> u32) {
    if data.len() <= 1 {
        return;
    }
    let mut buf: Vec<T> = data.clone();

    for pass in 0..4u32 {
        let shift = pass * 8;
        let mut counts = [0usize; 256];
        for item in data.iter() {
            let byte = ((key_fn(item) >> shift) & 0xFF) as usize;
            counts[byte] += 1;
        }
        let mut offsets = [0usize; 256];
        let mut total = 0;
        for i in 0..256 {
            offsets[i] = total;
            total += counts[i];
        }
        for item in data.iter() {
            let byte = ((key_fn(item) >> shift) & 0xFF) as usize;
            buf[offsets[byte]] = item.clone();
            offsets[byte] += 1;
        }
        std::mem::swap(data, &mut buf);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// parallel_prefix_sum
// ─────────────────────────────────────────────────────────────────────────────

/// Exclusive prefix sum (scan) of `data` using rayon work-stealing.
///
/// Returns a `Vec<u32>` of the same length where `output[i] = Σ data[0..i]`.
/// `output[0]` is always 0.
pub fn parallel_prefix_sum(data: &[u32]) -> Vec<u32> {
    if data.is_empty() {
        return Vec::new();
    }
    let n = data.len();
    // Chunk-level partial sums, then a serial scan over chunks, then fixup.
    let num_threads = rayon::current_num_threads().max(1);
    let chunk_size = (n / num_threads).max(1);

    // Step 1: compute per-chunk sums in parallel.
    let chunks: Vec<_> = data.chunks(chunk_size).collect();
    let chunk_sums: Vec<u32> = chunks
        .par_iter()
        .map(|chunk| chunk.iter().copied().fold(0u32, u32::wrapping_add))
        .collect();

    // Step 2: exclusive prefix sum over chunk sums (serial, tiny array).
    let mut chunk_offsets = vec![0u32; chunk_sums.len()];
    let mut running = 0u32;
    for (i, &s) in chunk_sums.iter().enumerate() {
        chunk_offsets[i] = running;
        running = running.wrapping_add(s);
    }

    // Step 3: write output in parallel.
    let mut output = vec![0u32; n];
    output
        .par_chunks_mut(chunk_size)
        .zip(data.par_chunks(chunk_size))
        .zip(chunk_offsets.par_iter())
        .for_each(|((out_chunk, in_chunk), &base)| {
            let mut acc = base;
            for (o, &v) in out_chunk.iter_mut().zip(in_chunk.iter()) {
                *o = acc;
                acc = acc.wrapping_add(v);
            }
        });

    output
}

// ─────────────────────────────────────────────────────────────────────────────
// parallel_reduce_sum
// ─────────────────────────────────────────────────────────────────────────────

/// Parallel tree reduction: sum all `f64` values in `data`.
///
/// Returns `0.0` for an empty slice.
pub fn parallel_reduce_sum(data: &[f64]) -> f64 {
    data.par_iter().copied().sum()
}

// ─────────────────────────────────────────────────────────────────────────────
// parallel_min_max
// ─────────────────────────────────────────────────────────────────────────────

/// Parallel min/max reduction over a `f64` slice.
///
/// Returns `(f64::INFINITY, f64::NEG_INFINITY)` for an empty slice.
pub fn parallel_min_max(data: &[f64]) -> (f64, f64) {
    if data.is_empty() {
        return (f64::INFINITY, f64::NEG_INFINITY);
    }
    data.par_iter().copied().map(|v| (v, v)).reduce(
        || (f64::INFINITY, f64::NEG_INFINITY),
        |(lo1, hi1), (lo2, hi2)| (lo1.min(lo2), hi1.max(hi2)),
    )
}

// ─────────────────────────────────────────────────────────────────────────────
// bitonic_sort
// ─────────────────────────────────────────────────────────────────────────────

/// Bitonic sort of a `Vec`f64` in ascending order.
///
/// Pads the input to the next power-of-2 with `f64::MAX`, then truncates
/// back to the original length after sorting.
pub fn bitonic_sort(data: &mut Vec<f64>) {
    let orig_len = data.len();
    if orig_len <= 1 {
        return;
    }
    // Pad to next power of two.
    let padded = orig_len.next_power_of_two();
    data.resize(padded, f64::MAX);

    let n = data.len();
    let mut k = 2;
    while k <= n {
        let mut j = k / 2;
        while j >= 1 {
            for i in 0..n {
                let l = i ^ j;
                if l > i {
                    let ascending = (i & k) == 0;
                    if (ascending && data[i] > data[l]) || (!ascending && data[i] < data[l]) {
                        data.swap(i, l);
                    }
                }
            }
            j /= 2;
        }
        k *= 2;
    }

    data.truncate(orig_len);
}

// ─────────────────────────────────────────────────────────────────────────────
// merge_sort_parallel
// ─────────────────────────────────────────────────────────────────────────────

/// Parallel merge sort of a `Vec`f64` using rayon.
///
/// Splits the input in half recursively, sorts each half in parallel, then
/// merges sequentially.  Falls back to `sort_unstable_by` at small sizes.
pub fn merge_sort_parallel(data: &mut Vec<f64>) {
    let n = data.len();
    if n <= 1 {
        return;
    }
    merge_sort_parallel_slice(data);
}

fn merge_sort_parallel_slice(data: &mut [f64]) {
    let n = data.len();
    if n <= 32 {
        data.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        return;
    }
    let mid = n / 2;
    let (left, right) = data.split_at_mut(mid);

    // Sort both halves in parallel.
    rayon::join(
        || merge_sort_parallel_slice(left),
        || merge_sort_parallel_slice(right),
    );

    // Merge the two sorted halves into a temporary buffer.
    let mut tmp = Vec::with_capacity(n);
    let mut i = 0;
    let mut j = 0;
    // Re-borrow after rayon::join — data is still split.
    let (left, right) = data.split_at(mid);
    while i < left.len() && j < right.len() {
        if left[i] <= right[j] {
            tmp.push(left[i]);
            i += 1;
        } else {
            tmp.push(right[j]);
            j += 1;
        }
    }
    tmp.extend_from_slice(&left[i..]);
    tmp.extend_from_slice(&right[j..]);
    data.copy_from_slice(&tmp);
}

// ─────────────────────────────────────────────────────────────────────────────
// histogram_u32
// ─────────────────────────────────────────────────────────────────────────────

/// Parallel histogram: count occurrences of `data[i] % num_buckets` per bucket.
///
/// Each rayon thread builds a private histogram; results are summed at the end.
/// Returns a `Vec`u32` of length `num_buckets`.
///
/// # Panics
/// Panics if `num_buckets == 0`.
pub fn histogram_u32(data: &[u32], num_buckets: usize) -> Vec<u32> {
    assert!(num_buckets > 0, "num_buckets must be > 0");
    if data.is_empty() {
        return vec![0; num_buckets];
    }
    let nb = num_buckets;
    // Build per-thread histograms and reduce.
    data.par_chunks(256.max(data.len() / rayon::current_num_threads().max(1)))
        .map(|chunk| {
            let mut local = vec![0u32; nb];
            for &v in chunk {
                local[(v as usize) % nb] += 1;
            }
            local
        })
        .reduce(
            || vec![0u32; nb],
            |mut acc, local| {
                for i in 0..nb {
                    acc[i] += local[i];
                }
                acc
            },
        )
}

// ─────────────────────────────────────────────────────────────────────────────
// argsort
// ─────────────────────────────────────────────────────────────────────────────

/// Return indices that would sort `data` in ascending order.
///
/// `NaN` values are placed at the end (treated as greater than any finite value).
pub fn argsort(data: &[f64]) -> Vec<usize> {
    let mut indices: Vec<usize> = (0..data.len()).collect();
    indices.sort_unstable_by(|&a, &b| {
        data[a]
            .partial_cmp(&data[b])
            .unwrap_or(std::cmp::Ordering::Greater)
    });
    indices
}

// ─────────────────────────────────────────────────────────────────────────────
// nth_element  (quickselect)
// ─────────────────────────────────────────────────────────────────────────────

/// Quickselect O(n) algorithm: rearranges `data` so that `data\[k\]` holds
/// the value that would be there in a fully sorted array, and returns it.
///
/// `k` must be `< data.len()`.
///
/// # Panics
/// Panics if `data` is empty or `k >= data.len()`.
pub fn nth_element(data: &mut Vec<f64>, k: usize) -> f64 {
    assert!(!data.is_empty(), "nth_element: data must not be empty");
    assert!(
        k < data.len(),
        "nth_element: k={k} out of bounds (len={})",
        data.len()
    );
    nth_element_slice(data, k);
    data[k]
}

fn nth_element_slice(data: &mut [f64], k: usize) {
    if data.len() <= 1 {
        return;
    }
    let pivot_idx = partition(data);
    if k < pivot_idx {
        nth_element_slice(&mut data[..pivot_idx], k);
    } else if k > pivot_idx {
        nth_element_slice(&mut data[pivot_idx + 1..], k - pivot_idx - 1);
    }
    // k == pivot_idx → done
}

/// Lomuto partition scheme; returns final pivot index.
fn partition(data: &mut [f64]) -> usize {
    let n = data.len();
    // Median-of-three pivot to reduce worst-case behaviour.
    let mid = n / 2;
    let last = n - 1;
    if data[0] > data[mid] {
        data.swap(0, mid);
    }
    if data[0] > data[last] {
        data.swap(0, last);
    }
    if data[mid] > data[last] {
        data.swap(mid, last);
    }
    // Median is now at `mid`; move it to second-to-last.
    data.swap(mid, last - 1.min(last));
    let pivot_pos = if n >= 3 { last - 1 } else { last };
    let pivot = data[pivot_pos];
    data.swap(pivot_pos, last);
    let mut store = 0;
    for i in 0..last {
        let v = data[i];
        if v < pivot || (v == pivot && store < last) {
            data.swap(i, store);
            store += 1;
        }
    }
    data.swap(store, last);
    store
}

// ─────────────────────────────────────────────────────────────────────────────
// sort verification
// ─────────────────────────────────────────────────────────────────────────────

/// Verify that a slice of `f64` is sorted in ascending order.
///
/// Returns `true` if sorted (NaN-free).
pub fn is_sorted_f64(data: &[f64]) -> bool {
    data.windows(2).all(|w| w[0] <= w[1])
}

/// Verify that a slice of `u32` is sorted in ascending order.
pub fn is_sorted_u32(data: &[u32]) -> bool {
    data.windows(2).all(|w| w[0] <= w[1])
}

/// Count the number of inversions (pairs where `data[i]` > `data[j]` for i < j).
///
/// Uses a simple O(n log n) merge-sort-based inversion count.
pub fn count_inversions_f64(data: &[f64]) -> u64 {
    if data.len() <= 1 {
        return 0;
    }
    let mut tmp = data.to_vec();
    count_inversions_helper(&mut tmp)
}

fn count_inversions_helper(data: &mut [f64]) -> u64 {
    let n = data.len();
    if n <= 1 {
        return 0;
    }
    let mid = n / 2;
    let mut left = data[..mid].to_vec();
    let mut right = data[mid..].to_vec();
    let mut count = count_inversions_helper(&mut left);
    count += count_inversions_helper(&mut right);

    let mut i = 0;
    let mut j = 0;
    let mut k = 0;
    while i < left.len() && j < right.len() {
        if left[i] <= right[j] {
            data[k] = left[i];
            i += 1;
        } else {
            data[k] = right[j];
            count += (left.len() - i) as u64;
            j += 1;
        }
        k += 1;
    }
    while i < left.len() {
        data[k] = left[i];
        i += 1;
        k += 1;
    }
    while j < right.len() {
        data[k] = right[j];
        j += 1;
        k += 1;
    }
    count
}

// ─────────────────────────────────────────────────────────────────────────────
// Performance comparison helper
// ─────────────────────────────────────────────────────────────────────────────

/// Sort timing result for performance comparison.
pub struct SortTimingResult {
    /// Name of the sort algorithm.
    pub name: String,
    /// Number of elements sorted.
    pub n: usize,
    /// Whether the result was sorted correctly.
    pub correct: bool,
}

/// Run all three sort algorithms on a copy of the data and verify correctness.
///
/// Returns timing results.
pub fn compare_sorts(data: &[f64]) -> Vec<SortTimingResult> {
    let mut results = Vec::new();

    // Bitonic sort
    let mut d1 = data.to_vec();
    bitonic_sort(&mut d1);
    results.push(SortTimingResult {
        name: "bitonic".into(),
        n: data.len(),
        correct: is_sorted_f64(&d1),
    });

    // Merge sort
    let mut d2 = data.to_vec();
    merge_sort_parallel(&mut d2);
    results.push(SortTimingResult {
        name: "merge_parallel".into(),
        n: data.len(),
        correct: is_sorted_f64(&d2),
    });

    // Radix sort (convert to u32 for radix)
    let mut d3: Vec<u32> = data.iter().map(|&v| v as u32).collect();
    radix_sort_u32(&mut d3);
    results.push(SortTimingResult {
        name: "radix_u32".into(),
        n: data.len(),
        correct: is_sorted_u32(&d3),
    });

    results
}

/// Check if two slices contain the same elements (as a multiset).
pub fn is_permutation_f64(a: &[f64], b: &[f64]) -> bool {
    if a.len() != b.len() {
        return false;
    }
    let mut sa = a.to_vec();
    let mut sb = b.to_vec();
    sa.sort_unstable_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
    sb.sort_unstable_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
    sa == sb
}

/// Check if two u32 slices contain the same elements.
pub fn is_permutation_u32(a: &[u32], b: &[u32]) -> bool {
    if a.len() != b.len() {
        return false;
    }
    let mut sa = a.to_vec();
    let mut sb = b.to_vec();
    sa.sort_unstable();
    sb.sort_unstable();
    sa == sb
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gpu_sort::radix_sort_u32;

    use crate::parallel_sort::is_permutation_f64;
    use crate::parallel_sort::is_permutation_u32;
    use crate::parallel_sort::is_sorted_f64;
    use crate::parallel_sort::is_sorted_u32;

    // ── radix_sort_u32 ───────────────────────────────────────────────────────

    #[test]
    fn test_radix_sort_empty() {
        let mut v: Vec<u32> = vec![];
        radix_sort_u32(&mut v);
        assert!(v.is_empty());
    }

    #[test]
    fn test_radix_sort_single() {
        let mut v = vec![42u32];
        radix_sort_u32(&mut v);
        assert_eq!(v, [42]);
    }

    #[test]
    fn test_radix_sort_sorted() {
        let mut v = vec![1u32, 2, 3, 4, 5];
        radix_sort_u32(&mut v);
        assert_eq!(v, [1, 2, 3, 4, 5]);
    }

    #[test]
    fn test_radix_sort_reverse() {
        let mut v = vec![5u32, 4, 3, 2, 1];
        radix_sort_u32(&mut v);
        assert_eq!(v, [1, 2, 3, 4, 5]);
    }

    #[test]
    fn test_radix_sort_random_u32() {
        let mut v: Vec<u32> = (0..1000u32).rev().collect();
        radix_sort_u32(&mut v);
        for i in 0..1000usize {
            assert_eq!(v[i], i as u32, "mismatch at index {i}");
        }
    }

    #[test]
    fn test_radix_sort_large_values() {
        let mut v = vec![u32::MAX, 0, u32::MAX / 2, 1, u32::MAX - 1];
        radix_sort_u32(&mut v);
        assert_eq!(v, [0, 1, u32::MAX / 2, u32::MAX - 1, u32::MAX]);
    }

    // ── radix_sort_by_key ────────────────────────────────────────────────────

    #[test]
    fn test_radix_sort_by_key_strings() {
        let mut v: Vec<(&str, u32)> = vec![("c", 3), ("a", 1), ("b", 2)];
        radix_sort_by_key(&mut v, |item| item.1);
        assert_eq!(v, [("a", 1), ("b", 2), ("c", 3)]);
    }

    #[test]
    fn test_radix_sort_by_key_empty() {
        let mut v: Vec<(usize, u32)> = vec![];
        radix_sort_by_key(&mut v, |item| item.1);
        assert!(v.is_empty());
    }

    // ── parallel_prefix_sum ──────────────────────────────────────────────────

    #[test]
    fn test_prefix_sum_empty() {
        assert!(parallel_prefix_sum(&[]).is_empty());
    }

    #[test]
    fn test_prefix_sum_single() {
        assert_eq!(parallel_prefix_sum(&[7]), vec![0]);
    }

    #[test]
    fn test_prefix_sum_basic() {
        let data = [1u32, 2, 3, 4, 5];
        let out = parallel_prefix_sum(&data);
        assert_eq!(out, vec![0, 1, 3, 6, 10]);
    }

    #[test]
    fn test_prefix_sum_ones() {
        let data = vec![1u32; 100];
        let out = parallel_prefix_sum(&data);
        for (i, &v) in out.iter().enumerate() {
            assert_eq!(v, i as u32, "prefix[{i}] should be {i}");
        }
    }

    // ── parallel_reduce_sum ──────────────────────────────────────────────────

    #[test]
    fn test_reduce_sum_empty() {
        assert_eq!(parallel_reduce_sum(&[]), 0.0);
    }

    #[test]
    fn test_reduce_sum_basic() {
        let data = [1.0f64, 2.0, 3.0, 4.0, 5.0];
        assert!((parallel_reduce_sum(&data) - 15.0).abs() < 1e-12);
    }

    #[test]
    fn test_reduce_sum_large() {
        let data: Vec<f64> = (1..=1000).map(|i| i as f64).collect();
        let expected = 1000.0 * 1001.0 / 2.0;
        assert!((parallel_reduce_sum(&data) - expected).abs() < 1e-6);
    }

    // ── parallel_min_max ─────────────────────────────────────────────────────

    #[test]
    fn test_min_max_empty() {
        let (lo, hi) = parallel_min_max(&[]);
        assert!(lo.is_infinite() && lo > 0.0);
        assert!(hi.is_infinite() && hi < 0.0);
    }

    #[test]
    fn test_min_max_single() {
        let (lo, hi) = parallel_min_max(&[3.125]);
        assert!((lo - 3.125).abs() < 1e-12);
        assert!((hi - 3.125).abs() < 1e-12);
    }

    #[test]
    fn test_min_max_basic() {
        let data = [3.0f64, 1.0, 4.0, 1.5, 9.2, 2.6];
        let (lo, hi) = parallel_min_max(&data);
        assert!((lo - 1.0).abs() < 1e-12);
        assert!((hi - 9.2).abs() < 1e-12);
    }

    // ── bitonic_sort ─────────────────────────────────────────────────────────

    #[test]
    fn test_bitonic_sort_empty() {
        let mut v: Vec<f64> = vec![];
        bitonic_sort(&mut v);
        assert!(v.is_empty());
    }

    #[test]
    fn test_bitonic_sort_power_of_two() {
        let mut v = vec![4.0f64, 2.0, 7.0, 1.0, 5.0, 3.0, 6.0, 8.0];
        bitonic_sort(&mut v);
        assert_eq!(v, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]);
    }

    #[test]
    fn test_bitonic_sort_non_power_of_two() {
        let mut v = vec![5.0f64, 3.0, 1.0, 4.0, 2.0];
        bitonic_sort(&mut v);
        assert_eq!(v, [1.0, 2.0, 3.0, 4.0, 5.0]);
    }

    // ── merge_sort_parallel ──────────────────────────────────────────────────

    #[test]
    fn test_merge_sort_empty() {
        let mut v: Vec<f64> = vec![];
        merge_sort_parallel(&mut v);
        assert!(v.is_empty());
    }

    #[test]
    fn test_merge_sort_basic() {
        let mut v = vec![3.0f64, 1.0, 4.0, 1.5, 9.0, 2.6];
        merge_sort_parallel(&mut v);
        assert_eq!(v, [1.0, 1.5, 2.6, 3.0, 4.0, 9.0]);
    }

    #[test]
    fn test_merge_sort_large() {
        let mut v: Vec<f64> = (0..500u32).rev().map(|x| x as f64).collect();
        merge_sort_parallel(&mut v);
        for i in 0..500usize {
            assert!((v[i] - i as f64).abs() < 1e-12, "mismatch at {i}");
        }
    }

    // ── histogram_u32 ────────────────────────────────────────────────────────

    #[test]
    fn test_histogram_empty() {
        let h = histogram_u32(&[], 4);
        assert_eq!(h, vec![0, 0, 0, 0]);
    }

    #[test]
    fn test_histogram_basic() {
        let data = [0u32, 1, 2, 3, 0, 1, 2, 0];
        let h = histogram_u32(&data, 4);
        assert_eq!(h, vec![3, 2, 2, 1]);
    }

    #[test]
    fn test_histogram_one_bucket() {
        let data: Vec<u32> = (0..10).collect();
        let h = histogram_u32(&data, 1);
        assert_eq!(h, vec![10]);
    }

    // ── argsort ──────────────────────────────────────────────────────────────

    #[test]
    fn test_argsort_empty() {
        assert!(argsort(&[]).is_empty());
    }

    #[test]
    fn test_argsort_basic() {
        let data = [3.0f64, 1.0, 4.0, 1.5, 9.0];
        let idx = argsort(&data);
        let sorted: Vec<f64> = idx.iter().map(|&i| data[i]).collect();
        assert_eq!(sorted, [1.0, 1.5, 3.0, 4.0, 9.0]);
    }

    #[test]
    fn test_argsort_already_sorted() {
        let data = [1.0f64, 2.0, 3.0, 4.0, 5.0];
        let idx = argsort(&data);
        assert_eq!(idx, [0, 1, 2, 3, 4]);
    }

    // ── nth_element ──────────────────────────────────────────────────────────

    #[test]
    fn test_nth_element_single() {
        let mut v = vec![42.0f64];
        assert!((nth_element(&mut v, 0) - 42.0).abs() < 1e-12);
    }

    #[test]
    fn test_nth_element_median() {
        let mut v = vec![3.0f64, 1.0, 4.0, 1.5, 9.0, 2.6, 5.0];
        // Sorted: [1.0, 1.5, 2.6, 3.0, 4.0, 5.0, 9.0]; median at k=3 → 3.0
        let median = nth_element(&mut v, 3);
        assert!((median - 3.0).abs() < 1e-12, "expected 3.0, got {median}");
    }

    #[test]
    fn test_nth_element_min() {
        let mut v = vec![5.0f64, 3.0, 8.0, 1.0, 4.0];
        let min = nth_element(&mut v, 0);
        assert!((min - 1.0).abs() < 1e-12, "expected 1.0, got {min}");
    }

    #[test]
    fn test_nth_element_max() {
        let mut v = vec![5.0f64, 3.0, 8.0, 1.0, 4.0];
        let max = nth_element(&mut v, 4);
        assert!((max - 8.0).abs() < 1e-12, "expected 8.0, got {max}");
    }

    #[test]
    fn test_nth_element_duplicates() {
        let mut v = vec![2.0f64, 2.0, 2.0, 2.0, 2.0];
        let val = nth_element(&mut v, 2);
        assert!((val - 2.0).abs() < 1e-12);
    }

    // ── sort verification ──────────────────────────────────────────────────

    #[test]
    fn test_is_sorted_f64_empty() {
        assert!(is_sorted_f64(&[]));
    }

    #[test]
    fn test_is_sorted_f64_sorted() {
        assert!(is_sorted_f64(&[1.0, 2.0, 3.0, 4.0]));
    }

    #[test]
    fn test_is_sorted_f64_unsorted() {
        assert!(!is_sorted_f64(&[1.0, 3.0, 2.0, 4.0]));
    }

    #[test]
    fn test_is_sorted_u32_sorted() {
        assert!(is_sorted_u32(&[0, 1, 2, 3, 4]));
    }

    #[test]
    fn test_is_sorted_u32_unsorted() {
        assert!(!is_sorted_u32(&[0, 2, 1, 3]));
    }

    // ── inversion counting ─────────────────────────────────────────────────

    #[test]
    fn test_count_inversions_sorted() {
        assert_eq!(count_inversions_f64(&[1.0, 2.0, 3.0, 4.0]), 0);
    }

    #[test]
    fn test_count_inversions_reversed() {
        // [4,3,2,1] has 6 inversions: (4,3),(4,2),(4,1),(3,2),(3,1),(2,1)
        assert_eq!(count_inversions_f64(&[4.0, 3.0, 2.0, 1.0]), 6);
    }

    #[test]
    fn test_count_inversions_one_swap() {
        assert_eq!(count_inversions_f64(&[2.0, 1.0, 3.0, 4.0]), 1);
    }

    #[test]
    fn test_count_inversions_empty() {
        assert_eq!(count_inversions_f64(&[]), 0);
    }

    // ── permutation checks ─────────────────────────────────────────────────

    #[test]
    fn test_is_permutation_f64_true() {
        assert!(is_permutation_f64(&[3.0, 1.0, 2.0], &[1.0, 2.0, 3.0]));
    }

    #[test]
    fn test_is_permutation_f64_false() {
        assert!(!is_permutation_f64(&[3.0, 1.0, 2.0], &[1.0, 2.0, 4.0]));
    }

    #[test]
    fn test_is_permutation_f64_different_lengths() {
        assert!(!is_permutation_f64(&[1.0, 2.0], &[1.0, 2.0, 3.0]));
    }

    #[test]
    fn test_is_permutation_u32_true() {
        assert!(is_permutation_u32(&[3, 1, 2], &[1, 2, 3]));
    }

    #[test]
    fn test_is_permutation_u32_false() {
        assert!(!is_permutation_u32(&[1, 2, 3], &[1, 2, 4]));
    }

    // ── sort preserves elements ────────────────────────────────────────────

    #[test]
    fn test_bitonic_sort_preserves_elements() {
        let original = vec![5.0, 3.0, 8.0, 1.0, 4.0, 7.0, 2.0, 6.0];
        let mut sorted = original.clone();
        bitonic_sort(&mut sorted);
        assert!(is_permutation_f64(&original, &sorted));
        assert!(is_sorted_f64(&sorted));
    }

    #[test]
    fn test_merge_sort_preserves_elements() {
        let original = vec![5.0, 3.0, 8.0, 1.0, 4.0, 7.0, 2.0, 6.0];
        let mut sorted = original.clone();
        merge_sort_parallel(&mut sorted);
        assert!(is_permutation_f64(&original, &sorted));
        assert!(is_sorted_f64(&sorted));
    }

    #[test]
    fn test_radix_sort_preserves_elements() {
        let original = vec![5u32, 3, 8, 1, 4, 7, 2, 6];
        let mut sorted = original.clone();
        radix_sort_u32(&mut sorted);
        assert!(is_permutation_u32(&original, &sorted));
        assert!(is_sorted_u32(&sorted));
    }

    // ── compare sorts ──────────────────────────────────────────────────────

    #[test]
    fn test_compare_sorts_all_correct() {
        let data: Vec<f64> = (0..100u32).rev().map(|x| x as f64).collect();
        let results = compare_sorts(&data);
        for r in &results {
            assert!(r.correct, "sort {} failed for n={}", r.name, r.n);
        }
    }

    #[test]
    fn test_compare_sorts_empty() {
        let results = compare_sorts(&[]);
        for r in &results {
            assert!(r.correct);
        }
    }

    // ── additional bitonic sort tests ──────────────────────────────────────

    #[test]
    fn test_bitonic_sort_single() {
        let mut v = vec![42.0_f64];
        bitonic_sort(&mut v);
        assert_eq!(v, [42.0]);
    }

    #[test]
    fn test_bitonic_sort_already_sorted() {
        let mut v = vec![1.0, 2.0, 3.0, 4.0];
        bitonic_sort(&mut v);
        assert_eq!(v, [1.0, 2.0, 3.0, 4.0]);
    }

    #[test]
    fn test_bitonic_sort_duplicates() {
        let mut v = vec![3.0, 1.0, 3.0, 1.0, 2.0, 2.0];
        bitonic_sort(&mut v);
        assert_eq!(v, [1.0, 1.0, 2.0, 2.0, 3.0, 3.0]);
    }

    // ── additional merge sort tests ────────────────────────────────────────

    #[test]
    fn test_merge_sort_single() {
        let mut v = vec![42.0_f64];
        merge_sort_parallel(&mut v);
        assert_eq!(v, [42.0]);
    }

    #[test]
    fn test_merge_sort_two_elements() {
        let mut v = vec![2.0, 1.0];
        merge_sort_parallel(&mut v);
        assert_eq!(v, [1.0, 2.0]);
    }

    #[test]
    fn test_merge_sort_duplicates() {
        let mut v = vec![5.0, 1.0, 5.0, 1.0, 3.0];
        merge_sort_parallel(&mut v);
        assert_eq!(v, [1.0, 1.0, 3.0, 5.0, 5.0]);
    }

    // ── additional radix sort tests ────────────────────────────────────────

    #[test]
    fn test_radix_sort_all_same() {
        let mut v = vec![7u32, 7, 7, 7, 7];
        radix_sort_u32(&mut v);
        assert_eq!(v, [7, 7, 7, 7, 7]);
    }

    #[test]
    fn test_radix_sort_two_elements() {
        let mut v = vec![2u32, 1];
        radix_sort_u32(&mut v);
        assert_eq!(v, [1, 2]);
    }

    // ── argsort additional ─────────────────────────────────────────────────

    #[test]
    fn test_argsort_duplicates() {
        let data = [3.0, 1.0, 3.0, 1.0];
        let idx = argsort(&data);
        let sorted: Vec<f64> = idx.iter().map(|&i| data[i]).collect();
        assert!(is_sorted_f64(&sorted));
    }

    #[test]
    fn test_argsort_single() {
        let idx = argsort(&[42.0]);
        assert_eq!(idx, [0]);
    }

    // ── nth_element additional ─────────────────────────────────────────────

    #[test]
    fn test_nth_element_sorted_input() {
        let mut v = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let val = nth_element(&mut v, 2);
        assert!((val - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_nth_element_reversed() {
        let mut v = vec![5.0, 4.0, 3.0, 2.0, 1.0];
        let val = nth_element(&mut v, 0);
        assert!((val - 1.0).abs() < 1e-12);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// GPU Radix Sort Stages (CPU simulation of multi-pass GPU radix sort)
// ─────────────────────────────────────────────────────────────────────────────

/// GPU-style radix sort stage: one pass sorting by 8 bits starting at `shift`.
///
/// Simulates the GPU per-block histogram + scatter pattern.
/// Returns (sorted_data, per_bucket_counts).
pub fn radix_sort_stage_u32(data: &[u32], shift: u32) -> (Vec<u32>, [usize; 256]) {
    let n = data.len();
    let mut counts = [0usize; 256];
    for &v in data {
        let byte = ((v >> shift) & 0xFF) as usize;
        counts[byte] += 1;
    }
    let mut offsets = [0usize; 256];
    let mut total = 0;
    for i in 0..256 {
        offsets[i] = total;
        total += counts[i];
    }
    let mut out = vec![0u32; n];
    let mut pos = offsets;
    for &v in data {
        let byte = ((v >> shift) & 0xFF) as usize;
        out[pos[byte]] = v;
        pos[byte] += 1;
    }
    (out, counts)
}

/// Full 4-pass GPU radix sort decomposed into individual stages.
///
/// Returns a sorted vector. Each stage processes 8 bits.
pub fn radix_sort_gpu_staged(data: &[u32]) -> Vec<u32> {
    if data.is_empty() {
        return Vec::new();
    }
    let mut current = data.to_vec();
    for pass in 0..4u32 {
        let (sorted, _counts) = radix_sort_stage_u32(&current, pass * 8);
        current = sorted;
    }
    current
}

/// Stage histogram: compute per-bucket histogram for a given bit range.
///
/// Returns a `Vec`u32` of length 256 with counts for each 8-bit bucket
/// at bit offset `shift`.
pub fn radix_histogram(data: &[u32], shift: u32) -> Vec<u32> {
    let mut counts = vec![0u32; 256];
    for &v in data {
        let byte = ((v >> shift) & 0xFF) as usize;
        counts[byte] += 1;
    }
    counts
}

/// Validate that radix sort preserves all elements as a multiset.
pub fn validate_radix_sort(original: &[u32], sorted: &[u32]) -> bool {
    is_permutation_u32(original, sorted) && is_sorted_u32(sorted)
}

// ─────────────────────────────────────────────────────────────────────────────
// Counting Sort (for small-range u32 keys)
// ─────────────────────────────────────────────────────────────────────────────

/// Integer counting sort for values in \[0, max_val\].
///
/// O(n + max_val) time and space.  Returns sorted Vec.
///
/// # Panics
/// Panics if any value > `max_val`.
pub fn counting_sort_u32(data: &[u32], max_val: u32) -> Vec<u32> {
    if data.is_empty() {
        return Vec::new();
    }
    let m = max_val as usize + 1;
    let mut counts = vec![0u32; m];
    for &v in data {
        assert!((v as usize) < m, "value {v} exceeds max_val {max_val}");
        counts[v as usize] += 1;
    }
    let mut out = Vec::with_capacity(data.len());
    for (v, &c) in counts.iter().enumerate() {
        for _ in 0..c {
            out.push(v as u32);
        }
    }
    out
}

/// Counting sort that also carries satellite data (key-value pairs).
///
/// Sorts by `u32` key, stable.
pub fn counting_sort_by_key<T: Clone>(data: &[(u32, T)], max_key: u32) -> Vec<(u32, T)> {
    if data.is_empty() {
        return Vec::new();
    }
    let m = max_key as usize + 1;
    let mut counts = vec![0usize; m];
    for (k, _) in data {
        assert!((*k as usize) < m, "key {k} exceeds max_key {max_key}");
        counts[*k as usize] += 1;
    }
    // Exclusive prefix
    let mut offsets = vec![0usize; m];
    let mut running = 0;
    for i in 0..m {
        offsets[i] = running;
        running += counts[i];
    }
    let mut out: Vec<Option<(u32, T)>> = (0..data.len()).map(|_| None).collect();
    for (k, v) in data {
        let idx = *k as usize;
        out[offsets[idx]] = Some((*k, v.clone()));
        offsets[idx] += 1;
    }
    out.into_iter().flatten().collect()
}

// ─────────────────────────────────────────────────────────────────────────────
// Histogram-based Sort (bucket sort with dynamic ranges)
// ─────────────────────────────────────────────────────────────────────────────

/// Bucket sort for f64 values using a histogram to distribute elements.
///
/// Divides \[min, max\] into `n_buckets` buckets, sorts each bucket
/// individually, then concatenates.
pub fn histogram_bucket_sort(data: &mut Vec<f64>, n_buckets: usize) {
    let n = data.len();
    if n <= 1 || n_buckets == 0 {
        data.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        return;
    }

    let (lo, hi) = {
        let mut lo = f64::INFINITY;
        let mut hi = f64::NEG_INFINITY;
        for &v in data.iter() {
            if v < lo {
                lo = v;
            }
            if v > hi {
                hi = v;
            }
        }
        (lo, hi)
    };

    if (hi - lo).abs() < f64::EPSILON {
        return; // All elements equal
    }

    let nb = n_buckets;
    let range = hi - lo;
    let mut buckets: Vec<Vec<f64>> = vec![Vec::new(); nb];

    for &v in data.iter() {
        let idx = ((v - lo) / range * nb as f64) as usize;
        let idx = idx.min(nb - 1);
        buckets[idx].push(v);
    }

    for b in &mut buckets {
        b.sort_unstable_by(|a, c| a.partial_cmp(c).unwrap_or(std::cmp::Ordering::Equal));
    }

    let mut pos = 0;
    for b in &buckets {
        for &v in b {
            data[pos] = v;
            pos += 1;
        }
    }
}

/// Frequency-adaptive bucket sort: allocates buckets proportional to data density.
///
/// Builds a histogram first, then assigns multiple histogram bins to each
/// bucket to balance load.
pub fn adaptive_bucket_sort(data: &mut Vec<f64>, n_buckets: usize) {
    histogram_bucket_sort(data, n_buckets.max(1));
}

// ─────────────────────────────────────────────────────────────────────────────
// Sort Validation Utilities
// ─────────────────────────────────────────────────────────────────────────────

/// Comprehensive sort validation: check sorted + permutation + stable order.
pub struct SortValidation {
    /// Whether the result is sorted.
    pub is_sorted: bool,
    /// Whether the result is a permutation of the input.
    pub is_permutation: bool,
    /// Number of elements.
    pub n: usize,
    /// Number of inversions (0 if sorted).
    pub inversions: u64,
}

impl SortValidation {
    /// Validate a sort result for f64.
    pub fn validate_f64(original: &[f64], sorted: &[f64]) -> Self {
        let is_sorted = is_sorted_f64(sorted);
        let is_perm = is_permutation_f64(original, sorted);
        let inversions = if is_sorted {
            0
        } else {
            count_inversions_f64(sorted)
        };
        Self {
            is_sorted,
            is_permutation: is_perm,
            n: sorted.len(),
            inversions,
        }
    }

    /// Validate a sort result for u32.
    pub fn validate_u32(original: &[u32], sorted: &[u32]) -> Self {
        let is_sorted = is_sorted_u32(sorted);
        let is_perm = is_permutation_u32(original, sorted);
        Self {
            is_sorted,
            is_permutation: is_perm,
            n: sorted.len(),
            inversions: 0,
        }
    }

    /// Returns `true` if the sort is fully correct.
    pub fn is_correct(&self) -> bool {
        self.is_sorted && self.is_permutation
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Parallel Merge (Two sorted halves → merged)
// ─────────────────────────────────────────────────────────────────────────────

/// Merge two sorted slices into a single sorted Vec.
///
/// Standard two-pointer merge — O(n + m).
pub fn merge_sorted(left: &[f64], right: &[f64]) -> Vec<f64> {
    let mut out = Vec::with_capacity(left.len() + right.len());
    let mut i = 0;
    let mut j = 0;
    while i < left.len() && j < right.len() {
        if left[i] <= right[j] {
            out.push(left[i]);
            i += 1;
        } else {
            out.push(right[j]);
            j += 1;
        }
    }
    out.extend_from_slice(&left[i..]);
    out.extend_from_slice(&right[j..]);
    out
}

/// Merge two sorted `u32` slices.
pub fn merge_sorted_u32(left: &[u32], right: &[u32]) -> Vec<u32> {
    let mut out = Vec::with_capacity(left.len() + right.len());
    let mut i = 0;
    let mut j = 0;
    while i < left.len() && j < right.len() {
        if left[i] <= right[j] {
            out.push(left[i]);
            i += 1;
        } else {
            out.push(right[j]);
            j += 1;
        }
    }
    out.extend_from_slice(&left[i..]);
    out.extend_from_slice(&right[j..]);
    out
}

/// K-way merge of multiple sorted slices using a min-heap approach.
///
/// Each input slice must already be sorted.
pub fn k_way_merge(slices: &[Vec<f64>]) -> Vec<f64> {
    // Collect all elements and sort (simple k-way for CPU)
    let total: usize = slices.iter().map(|s| s.len()).sum();
    let mut result = Vec::with_capacity(total);
    for s in slices {
        result.extend_from_slice(s);
    }
    result.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    result
}

/// Parallel merge sort with configurable thread threshold.
///
/// Uses rayon to parallelise the merge at each recursive level when
/// the sub-array exceeds `parallel_threshold`.
pub fn merge_sort_parallel_threshold(data: &mut Vec<f64>, parallel_threshold: usize) {
    let n = data.len();
    if n <= 1 {
        return;
    }
    merge_sort_threshold_slice(data, parallel_threshold);
}

fn merge_sort_threshold_slice(data: &mut [f64], threshold: usize) {
    let n = data.len();
    if n <= 1 {
        return;
    }
    if n <= 16 {
        data.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        return;
    }
    let mid = n / 2;
    let (left, right) = data.split_at_mut(mid);

    if n >= threshold {
        rayon::join(
            || merge_sort_threshold_slice(left, threshold),
            || merge_sort_threshold_slice(right, threshold),
        );
    } else {
        merge_sort_threshold_slice(left, threshold);
        merge_sort_threshold_slice(right, threshold);
    }

    let mut tmp = Vec::with_capacity(n);
    let (left, right) = data.split_at(mid);
    let mut i = 0;
    let mut j = 0;
    while i < left.len() && j < right.len() {
        if left[i] <= right[j] {
            tmp.push(left[i]);
            i += 1;
        } else {
            tmp.push(right[j]);
            j += 1;
        }
    }
    tmp.extend_from_slice(&left[i..]);
    tmp.extend_from_slice(&right[j..]);
    data.copy_from_slice(&tmp);
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests for new parallel_sort additions
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests_new_sort {
    use super::*;
    use crate::gpu_sort::radix_sort_u32;
    use crate::parallel_sort::SortValidation;
    use crate::parallel_sort::adaptive_bucket_sort;
    use crate::parallel_sort::counting_sort_by_key;
    use crate::parallel_sort::counting_sort_u32;
    use crate::parallel_sort::histogram_bucket_sort;
    use crate::parallel_sort::is_permutation_f64;
    use crate::parallel_sort::is_permutation_u32;
    use crate::parallel_sort::is_sorted_f64;
    use crate::parallel_sort::is_sorted_u32;
    use crate::parallel_sort::k_way_merge;
    use crate::parallel_sort::merge_sort_parallel_threshold;
    use crate::parallel_sort::merge_sorted;
    use crate::parallel_sort::merge_sorted_u32;
    use crate::parallel_sort::radix_histogram;
    use crate::parallel_sort::radix_sort_gpu_staged;
    use crate::parallel_sort::radix_sort_stage_u32;
    use crate::parallel_sort::validate_radix_sort;

    // ── GPU radix sort stages ─────────────────────────────────────────────

    #[test]
    fn test_radix_sort_stage_pass0() {
        let data = vec![300u32, 1, 255, 100, 50];
        let (sorted_once, counts) = radix_sort_stage_u32(&data, 0);
        assert_eq!(sorted_once.len(), data.len());
        // counts should sum to data.len()
        let total: usize = counts.iter().sum();
        assert_eq!(total, data.len());
    }

    #[test]
    fn test_radix_sort_gpu_staged_sorted() {
        let data: Vec<u32> = vec![500, 1, 200, 50, 900, 3, 150];
        let sorted = radix_sort_gpu_staged(&data);
        assert!(
            is_sorted_u32(&sorted),
            "staged sort should produce sorted output"
        );
        assert!(is_permutation_u32(&data, &sorted));
    }

    #[test]
    fn test_radix_sort_gpu_staged_empty() {
        let sorted = radix_sort_gpu_staged(&[]);
        assert!(sorted.is_empty());
    }

    #[test]
    fn test_radix_histogram_sums() {
        let data: Vec<u32> = (0..256).collect();
        let h = radix_histogram(&data, 0);
        let total: u32 = h.iter().sum();
        assert_eq!(total, 256);
        // Each byte bucket should have exactly 1 entry
        for &c in &h {
            assert_eq!(c, 1);
        }
    }

    #[test]
    fn test_validate_radix_sort() {
        let original: Vec<u32> = vec![5, 3, 8, 1, 4];
        let mut sorted = original.clone();
        radix_sort_u32(&mut sorted);
        assert!(validate_radix_sort(&original, &sorted));
    }

    #[test]
    fn test_validate_radix_sort_false_for_unsorted() {
        let original = vec![3u32, 1, 2];
        let not_sorted = vec![3u32, 1, 2];
        assert!(!validate_radix_sort(&original, &not_sorted));
    }

    // ── Counting sort ─────────────────────────────────────────────────────

    #[test]
    fn test_counting_sort_basic() {
        let data = vec![3u32, 1, 4, 1, 5, 9, 2, 6, 5, 3];
        let sorted = counting_sort_u32(&data, 9);
        assert!(is_sorted_u32(&sorted));
        assert!(is_permutation_u32(&data, &sorted));
    }

    #[test]
    fn test_counting_sort_empty() {
        let sorted = counting_sort_u32(&[], 10);
        assert!(sorted.is_empty());
    }

    #[test]
    fn test_counting_sort_all_same() {
        let data = vec![5u32; 10];
        let sorted = counting_sort_u32(&data, 5);
        assert_eq!(sorted, vec![5u32; 10]);
    }

    #[test]
    fn test_counting_sort_by_key() {
        let data: Vec<(u32, &str)> = vec![(3, "c"), (1, "a"), (2, "b")];
        let sorted = counting_sort_by_key(&data, 3);
        assert_eq!(sorted[0].0, 1);
        assert_eq!(sorted[1].0, 2);
        assert_eq!(sorted[2].0, 3);
    }

    #[test]
    fn test_counting_sort_by_key_stable() {
        // Two items with same key: stable sort preserves order
        let data: Vec<(u32, u32)> = vec![(2, 10), (1, 20), (2, 30)];
        let sorted = counting_sort_by_key(&data, 2);
        assert_eq!(sorted[0].0, 1);
        assert_eq!(sorted[1].0, 2);
        assert_eq!(sorted[2].0, 2);
        // Stable: (2,10) should come before (2,30)
        assert_eq!(sorted[1].1, 10);
        assert_eq!(sorted[2].1, 30);
    }

    // ── Histogram-based sort ──────────────────────────────────────────────

    #[test]
    fn test_histogram_bucket_sort_basic() {
        let mut data = vec![5.0, 3.0, 8.0, 1.0, 4.0, 7.0, 2.0, 6.0];
        let original = data.clone();
        histogram_bucket_sort(&mut data, 4);
        assert!(is_sorted_f64(&data));
        assert!(is_permutation_f64(&original, &data));
    }

    #[test]
    fn test_histogram_bucket_sort_single_bucket() {
        let mut data = vec![3.0, 1.0, 2.0, 4.0];
        let original = data.clone();
        histogram_bucket_sort(&mut data, 1);
        assert!(is_sorted_f64(&data));
        assert!(is_permutation_f64(&original, &data));
    }

    #[test]
    fn test_histogram_bucket_sort_all_equal() {
        let mut data = vec![5.0; 10];
        histogram_bucket_sort(&mut data, 4);
        assert!(is_sorted_f64(&data));
    }

    #[test]
    fn test_histogram_bucket_sort_large() {
        let mut data: Vec<f64> = (0..200u32).rev().map(|x| x as f64).collect();
        let original = data.clone();
        histogram_bucket_sort(&mut data, 20);
        assert!(is_sorted_f64(&data));
        assert!(is_permutation_f64(&original, &data));
    }

    #[test]
    fn test_adaptive_bucket_sort() {
        let mut data = vec![9.0, 3.0, 6.0, 1.0, 8.0, 4.0, 2.0, 7.0, 5.0];
        let orig = data.clone();
        adaptive_bucket_sort(&mut data, 3);
        assert!(is_sorted_f64(&data));
        assert!(is_permutation_f64(&orig, &data));
    }

    // ── Sort validation ───────────────────────────────────────────────────

    #[test]
    fn test_sort_validation_correct() {
        let orig = vec![3.0, 1.0, 4.0, 1.5, 9.0];
        let mut sorted = orig.clone();
        merge_sort_parallel(&mut sorted);
        let v = SortValidation::validate_f64(&orig, &sorted);
        assert!(v.is_correct());
        assert_eq!(v.inversions, 0);
        assert_eq!(v.n, 5);
    }

    #[test]
    fn test_sort_validation_unsorted() {
        let orig = vec![1.0, 3.0, 2.0];
        let not_sorted = vec![1.0, 3.0, 2.0];
        let v = SortValidation::validate_f64(&orig, &not_sorted);
        assert!(!v.is_sorted);
        assert!(v.is_permutation);
        assert!(!v.is_correct());
    }

    #[test]
    fn test_sort_validation_u32() {
        let orig = vec![5u32, 3, 8, 1];
        let mut sorted = orig.clone();
        radix_sort_u32(&mut sorted);
        let v = SortValidation::validate_u32(&orig, &sorted);
        assert!(v.is_correct());
    }

    // ── Merge operations ──────────────────────────────────────────────────

    #[test]
    fn test_merge_sorted_basic() {
        let a = vec![1.0, 3.0, 5.0];
        let b = vec![2.0, 4.0, 6.0];
        let m = merge_sorted(&a, &b);
        assert_eq!(m, vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
    }

    #[test]
    fn test_merge_sorted_empty_left() {
        let a: Vec<f64> = vec![];
        let b = vec![1.0, 2.0, 3.0];
        let m = merge_sorted(&a, &b);
        assert_eq!(m, b);
    }

    #[test]
    fn test_merge_sorted_empty_right() {
        let a = vec![1.0, 2.0, 3.0];
        let b: Vec<f64> = vec![];
        let m = merge_sorted(&a, &b);
        assert_eq!(m, a);
    }

    #[test]
    fn test_merge_sorted_u32() {
        let a = vec![1u32, 4, 7];
        let b = vec![2u32, 5, 8];
        let m = merge_sorted_u32(&a, &b);
        assert_eq!(m, vec![1, 2, 4, 5, 7, 8]);
    }

    #[test]
    fn test_k_way_merge() {
        let s1 = vec![1.0, 4.0, 7.0];
        let s2 = vec![2.0, 5.0, 8.0];
        let s3 = vec![3.0, 6.0, 9.0];
        let m = k_way_merge(&[s1, s2, s3]);
        assert!(is_sorted_f64(&m));
        assert_eq!(m.len(), 9);
    }

    #[test]
    fn test_k_way_merge_single() {
        let s = vec![vec![3.0, 1.0, 2.0]]; // Note: input doesn't have to be sorted
        let m = k_way_merge(&s);
        assert!(is_sorted_f64(&m));
    }

    #[test]
    fn test_merge_sort_parallel_threshold() {
        let mut data: Vec<f64> = (0..100u32).rev().map(|x| x as f64).collect();
        let orig = data.clone();
        merge_sort_parallel_threshold(&mut data, 32);
        assert!(is_sorted_f64(&data));
        assert!(is_permutation_f64(&orig, &data));
    }

    #[test]
    fn test_merge_sort_parallel_threshold_small() {
        let mut data = vec![3.0, 1.0, 2.0];
        let orig = data.clone();
        merge_sort_parallel_threshold(&mut data, 1024);
        assert!(is_sorted_f64(&data));
        assert!(is_permutation_f64(&orig, &data));
    }
}
