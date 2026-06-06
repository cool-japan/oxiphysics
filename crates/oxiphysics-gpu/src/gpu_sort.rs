// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! GPU-style parallel sorting algorithms (CPU simulation).
//!
//! This module provides GPU-pattern algorithms for sorting and related
//! operations on `f32`/`u32` data:
//! - Bitonic sort (power-of-2 padded)
//! - LSD radix sort for `u32` and `f32`
//! - Exclusive prefix sum (Blelloch scan style)
//! - Parallel histogram
//! - Counting sort
//! - Morton Z-order sort for 3-D point clouds
//! - [`GpuSortBuffer`] — key+value buffer with `sort_pairs`
//! - Parallel merge

// ─────────────────────────────────────────────────────────────────────────────
// Bitonic sort (f32)
// ─────────────────────────────────────────────────────────────────────────────

/// Sort a `Vec`f32` in ascending order using the bitonic sort algorithm.
///
/// If `data.len()` is not a power of two, the vector is padded to the next
/// power of two with `f32::MAX` sentinel values, sorted, then truncated.
pub fn bitonic_sort(data: &mut Vec<f32>) {
    let orig = data.len();
    if orig <= 1 {
        return;
    }
    let padded = next_pow2(orig);
    data.resize(padded, f32::MAX);
    bitonic_sort_slice_f32(data);
    data.truncate(orig);
}

/// Sort a `Vec`T` by the `u32` key returned by `key_fn`.
///
/// Padding uses sentinel key `u32::MAX`; after sorting, only the first
/// `orig_len` elements (those with valid indices) are retained.
pub fn bitonic_sort_by_key<T: Clone>(data: &mut [T], key_fn: impl Fn(&T) -> u32) {
    let orig = data.len();
    if orig <= 1 {
        return;
    }
    let padded = next_pow2(orig);
    // Build (key, original-index) pairs and pad.
    let mut pairs: Vec<(u32, usize)> = data
        .iter()
        .enumerate()
        .map(|(i, v)| (key_fn(v), i))
        .collect();
    pairs.resize(padded, (u32::MAX, usize::MAX));

    // Bitonic sort the pairs.
    let n = pairs.len();
    let mut k = 2;
    while k <= n {
        let mut j = k / 2;
        while j >= 1 {
            for i in 0..n {
                let l = i ^ j;
                if l > i {
                    let ascending = (i & k) == 0;
                    let should_swap = if ascending {
                        pairs[i].0 > pairs[l].0
                    } else {
                        pairs[i].0 < pairs[l].0
                    };
                    if should_swap {
                        pairs.swap(i, l);
                    }
                }
            }
            j /= 2;
        }
        k *= 2;
    }

    // Reconstruct data: take the first `orig` sorted entries that have valid indices.
    let old: Vec<T> = data.to_vec();
    let mut out: Vec<T> = pairs
        .iter()
        .filter(|&&(_, idx)| idx < orig)
        .map(|&(_, idx)| old[idx].clone())
        .collect();
    out.truncate(orig);
    // In case some items were duplicated by the sentinel logic, ensure length.
    for (i, v) in out.into_iter().enumerate().take(orig) {
        data[i] = v;
    }
}

/// Internal: in-place bitonic sort on a slice whose length is a power of two.
fn bitonic_sort_slice_f32(data: &mut [f32]) {
    let n = data.len();
    let mut k = 2;
    while k <= n {
        let mut j = k / 2;
        while j >= 1 {
            for i in 0..n {
                let l = i ^ j;
                if l > i {
                    let ascending = (i & k) == 0;
                    let should_swap = if ascending {
                        data[i] > data[l]
                    } else {
                        data[i] < data[l]
                    };
                    if should_swap {
                        data.swap(i, l);
                    }
                }
            }
            j /= 2;
        }
        k *= 2;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Radix sort
// ─────────────────────────────────────────────────────────────────────────────

/// LSD radix sort for `u32` using 4 passes of 8-bit digits.
///
/// Stable, O(n) complexity for 32-bit keys.
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
            counts[((v >> shift) & 0xFF) as usize] += 1;
        }
        let mut offsets = [0usize; 256];
        let mut total = 0;
        for i in 0..256 {
            offsets[i] = total;
            total += counts[i];
        }
        for &v in data.iter() {
            let b = ((v >> shift) & 0xFF) as usize;
            buf[offsets[b]] = v;
            offsets[b] += 1;
        }
        std::mem::swap(data, &mut buf);
    }
}

/// Sort `f32` values by reinterpreting bits and flipping the sign bit for
/// negatives so that the full IEEE 754 ordering maps to unsigned integer order.
///
/// After sorting, the bits are un-flipped back to valid `f32` values.
pub fn radix_sort_f32(data: &mut [f32]) {
    if data.len() <= 1 {
        return;
    }
    // Map f32 to a sortable u32.
    let mut keys: Vec<u32> = data.iter().map(|&v| f32_to_sort_key(v)).collect();
    radix_sort_u32(&mut keys);
    for (dst, k) in data.iter_mut().zip(keys.iter()) {
        *dst = sort_key_to_f32(*k);
    }
}

/// Convert an `f32` to a sortable `u32` (flip sign bit; for negatives also
/// flip the remaining bits so the order is preserved).
#[inline]
fn f32_to_sort_key(v: f32) -> u32 {
    let bits = v.to_bits();
    if bits >> 31 == 0 {
        bits | 0x8000_0000 // positive: set sign bit
    } else {
        !bits // negative: flip all bits
    }
}

/// Inverse of [`f32_to_sort_key`].
#[inline]
fn sort_key_to_f32(key: u32) -> f32 {
    let bits = if key >> 31 != 0 {
        key & 0x7FFF_FFFF // positive
    } else {
        !key // negative
    };
    f32::from_bits(bits)
}

// ─────────────────────────────────────────────────────────────────────────────
// Prefix sum (Blelloch exclusive scan)
// ─────────────────────────────────────────────────────────────────────────────

/// Exclusive prefix sum (Blelloch scan) for `u32` values.
///
/// Returns a new `Vec`u32` where `result\[i\] = sum(data\[0..i\])`.
/// `result\[0\]` is always `0`.
pub fn prefix_sum(data: &[u32]) -> Vec<u32> {
    let mut result = Vec::with_capacity(data.len());
    let mut acc = 0u32;
    for &v in data {
        result.push(acc);
        acc = acc.wrapping_add(v);
    }
    result
}

// ─────────────────────────────────────────────────────────────────────────────
// Histogram
// ─────────────────────────────────────────────────────────────────────────────

/// Compute a histogram of `data` into `n_bins` equal-width bins over
/// `\[min_val, max_val)`.
///
/// Values outside the range are clamped into the first or last bin.
/// Returns a `Vec`u32` of length `n_bins`.
///
/// # Panics
/// Panics if `n_bins == 0`.
pub fn histogram(data: &[u32], n_bins: usize) -> Vec<u32> {
    assert!(n_bins > 0, "n_bins must be > 0");
    if data.is_empty() {
        return vec![0u32; n_bins];
    }
    let max_val = *data.iter().max().unwrap_or(&0) as u64 + 1;
    let mut bins = vec![0u32; n_bins];
    for &v in data {
        let idx = ((v as u64 * n_bins as u64) / max_val) as usize;
        let idx = idx.min(n_bins - 1);
        bins[idx] += 1;
    }
    bins
}

// ─────────────────────────────────────────────────────────────────────────────
// Counting sort
// ─────────────────────────────────────────────────────────────────────────────

/// Counting sort for `u32` values bounded by `max_val` (inclusive).
///
/// Creates a count array of size `max_val + 1` and reconstructs sorted data
/// from it. O(n + max_val) time and space.
pub fn counting_sort(data: &mut [u32], max_val: u32) {
    if data.len() <= 1 {
        return;
    }
    let size = (max_val as usize).saturating_add(1);
    let mut counts = vec![0u32; size];
    for &v in data.iter() {
        let idx = (v as usize).min(size - 1);
        counts[idx] += 1;
    }
    let mut pos = 0usize;
    for (val, &cnt) in counts.iter().enumerate() {
        for _ in 0..cnt {
            data[pos] = val as u32;
            pos += 1;
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Morton sort (Z-order curve for 3-D points)
// ─────────────────────────────────────────────────────────────────────────────

/// Sort 3-D `f32` points by their Morton Z-order (space-filling curve) code.
///
/// Coordinates are quantised to 10-bit integers before interleaving, which
/// gives 30-bit Morton codes that fit in a `u32`.
pub fn morton_sort_3d(points: &mut [[f32; 3]]) {
    if points.len() <= 1 {
        return;
    }
    // Compute bounding box.
    let mut lo = [f32::INFINITY; 3];
    let mut hi = [f32::NEG_INFINITY; 3];
    for p in points.iter() {
        for d in 0..3 {
            lo[d] = lo[d].min(p[d]);
            hi[d] = hi[d].max(p[d]);
        }
    }
    let scale: Vec<f32> = (0..3)
        .map(|d| {
            let range = hi[d] - lo[d];
            if range > 0.0 { 1023.0 / range } else { 0.0 }
        })
        .collect();

    let mut pairs: Vec<(u32, usize)> = points
        .iter()
        .enumerate()
        .map(|(i, p)| {
            let ix = ((p[0] - lo[0]) * scale[0]) as u32;
            let iy = ((p[1] - lo[1]) * scale[1]) as u32;
            let iz = ((p[2] - lo[2]) * scale[2]) as u32;
            (morton3(ix.min(1023), iy.min(1023), iz.min(1023)), i)
        })
        .collect();

    pairs.sort_unstable_by_key(|&(code, _)| code);

    let old: Vec<[f32; 3]> = points.to_vec();
    for (i, &(_, idx)) in pairs.iter().enumerate() {
        points[i] = old[idx];
    }
}

/// Interleave the lower 10 bits of x, y, z into a 30-bit Morton code.
fn morton3(x: u32, y: u32, z: u32) -> u32 {
    spread_bits(x) | (spread_bits(y) << 1) | (spread_bits(z) << 2)
}

/// Spread the lower 10 bits of `v` into every third bit position.
fn spread_bits(mut v: u32) -> u32 {
    v &= 0x3FF; // keep lower 10 bits
    v = (v | (v << 16)) & 0x030000FF;
    v = (v | (v << 8)) & 0x0300F00F;
    v = (v | (v << 4)) & 0x030C30C3;
    v = (v | (v << 2)) & 0x09249249;
    v
}

// ─────────────────────────────────────────────────────────────────────────────
// GpuSortBuffer
// ─────────────────────────────────────────────────────────────────────────────

/// A buffer abstraction that holds parallel key and value arrays.
///
/// Provides `sort_pairs` to co-sort both arrays by key using radix sort.
#[derive(Debug, Clone)]
pub struct GpuSortBuffer {
    /// Sort keys.
    pub keys: Vec<u32>,
    /// Associated values (same length as `keys`).
    pub values: Vec<u32>,
}

impl GpuSortBuffer {
    /// Create a new `GpuSortBuffer` with the given key and value arrays.
    ///
    /// # Panics
    /// Panics if `keys` and `values` have different lengths.
    pub fn new(keys: Vec<u32>, values: Vec<u32>) -> Self {
        assert_eq!(
            keys.len(),
            values.len(),
            "keys and values must have equal length"
        );
        Self { keys, values }
    }

    /// Create an empty buffer.
    pub fn empty() -> Self {
        Self {
            keys: Vec::new(),
            values: Vec::new(),
        }
    }

    /// Number of key-value pairs.
    pub fn len(&self) -> usize {
        self.keys.len()
    }

    /// Returns `true` if the buffer is empty.
    pub fn is_empty(&self) -> bool {
        self.keys.is_empty()
    }

    /// Sort keys and values together using LSD radix sort, stable by key.
    pub fn sort_pairs(&mut self) {
        if self.len() <= 1 {
            return;
        }
        let n = self.len();
        let mut key_buf = vec![0u32; n];
        let mut val_buf = vec![0u32; n];
        for pass in 0..4u32 {
            let shift = pass * 8;
            let mut counts = [0usize; 256];
            for &k in self.keys.iter() {
                counts[((k >> shift) & 0xFF) as usize] += 1;
            }
            let mut offsets = [0usize; 256];
            let mut total = 0;
            for i in 0..256 {
                offsets[i] = total;
                total += counts[i];
            }
            for (i, &k) in self.keys.iter().enumerate() {
                let b = ((k >> shift) & 0xFF) as usize;
                let dest = offsets[b];
                key_buf[dest] = k;
                val_buf[dest] = self.values[i];
                offsets[b] += 1;
            }
            std::mem::swap(&mut self.keys, &mut key_buf);
            std::mem::swap(&mut self.values, &mut val_buf);
        }
    }

    /// Append a key-value pair to the buffer.
    pub fn push(&mut self, key: u32, value: u32) {
        self.keys.push(key);
        self.values.push(value);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Parallel merge
// ─────────────────────────────────────────────────────────────────────────────

/// Merge two sorted `f32` slices into a new sorted `Vec`f32`.
///
/// Both inputs must already be sorted in non-decreasing order.
/// Uses a standard two-pointer merge (O(n+m)).
pub fn parallel_merge(left: &[f32], right: &[f32]) -> Vec<f32> {
    let mut result = Vec::with_capacity(left.len() + right.len());
    let (mut i, mut j) = (0, 0);
    while i < left.len() && j < right.len() {
        if left[i] <= right[j] {
            result.push(left[i]);
            i += 1;
        } else {
            result.push(right[j]);
            j += 1;
        }
    }
    result.extend_from_slice(&left[i..]);
    result.extend_from_slice(&right[j..]);
    result
}

// ─────────────────────────────────────────────────────────────────────────────
// Internal helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Round `n` up to the next power of two (returns 1 for 0).
fn next_pow2(n: usize) -> usize {
    if n == 0 {
        return 1;
    }
    let mut p = 1usize;
    while p < n {
        p <<= 1;
    }
    p
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── next_pow2 ─────────────────────────────────────────────────────────────

    #[test]
    fn test_next_pow2_zero() {
        assert_eq!(next_pow2(0), 1);
    }

    #[test]
    fn test_next_pow2_one() {
        assert_eq!(next_pow2(1), 1);
    }

    #[test]
    fn test_next_pow2_exact() {
        assert_eq!(next_pow2(8), 8);
    }

    #[test]
    fn test_next_pow2_non_exact() {
        assert_eq!(next_pow2(9), 16);
    }

    // ── bitonic_sort ──────────────────────────────────────────────────────────

    #[test]
    fn test_bitonic_sort_power_of_two() {
        let mut data = vec![4.0f32, 2.0, 7.0, 1.0, 8.0, 3.0, 6.0, 5.0];
        bitonic_sort(&mut data);
        assert_eq!(data, vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]);
    }

    #[test]
    fn test_bitonic_sort_non_power_of_two() {
        let mut data = vec![5.0f32, 3.0, 1.0, 4.0, 2.0];
        bitonic_sort(&mut data);
        assert_eq!(data, vec![1.0, 2.0, 3.0, 4.0, 5.0]);
    }

    #[test]
    fn test_bitonic_sort_empty() {
        let mut data: Vec<f32> = vec![];
        bitonic_sort(&mut data);
        assert!(data.is_empty());
    }

    #[test]
    fn test_bitonic_sort_single() {
        let mut data = vec![42.0f32];
        bitonic_sort(&mut data);
        assert_eq!(data, vec![42.0]);
    }

    #[test]
    fn test_bitonic_sort_already_sorted() {
        let mut data = vec![1.0f32, 2.0, 3.0, 4.0];
        bitonic_sort(&mut data);
        assert_eq!(data, vec![1.0, 2.0, 3.0, 4.0]);
    }

    #[test]
    fn test_bitonic_sort_reverse() {
        let mut data = vec![8.0f32, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0];
        bitonic_sort(&mut data);
        assert_eq!(data, vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]);
    }

    #[test]
    fn test_bitonic_sort_duplicates() {
        let mut data = vec![3.0f32, 1.0, 3.0, 2.0, 1.0];
        bitonic_sort(&mut data);
        assert_eq!(data, vec![1.0, 1.0, 2.0, 3.0, 3.0]);
    }

    #[test]
    fn test_bitonic_sort_large_non_pow2() {
        let mut data: Vec<f32> = (0..100).map(|i| (100 - i) as f32).collect();
        bitonic_sort(&mut data);
        for i in 0..data.len() - 1 {
            assert!(data[i] <= data[i + 1]);
        }
    }

    // ── bitonic_sort_by_key ───────────────────────────────────────────────────

    #[test]
    fn test_bitonic_sort_by_key_u32() {
        let mut data = vec![30u32, 10, 20, 40];
        bitonic_sort_by_key(&mut data, |x| *x);
        assert_eq!(data, vec![10, 20, 30, 40]);
    }

    #[test]
    fn test_bitonic_sort_by_key_empty() {
        let mut data: Vec<u32> = vec![];
        bitonic_sort_by_key(&mut data, |x| *x);
        assert!(data.is_empty());
    }

    #[test]
    fn test_bitonic_sort_by_key_single() {
        let mut data = vec![42u32];
        bitonic_sort_by_key(&mut data, |x| *x);
        assert_eq!(data, vec![42]);
    }

    // ── radix_sort_u32 ────────────────────────────────────────────────────────

    #[test]
    fn test_radix_sort_u32_basic() {
        let mut data = vec![5u32, 3, 8, 1, 9, 2, 7, 4, 6];
        radix_sort_u32(&mut data);
        assert_eq!(data, vec![1, 2, 3, 4, 5, 6, 7, 8, 9]);
    }

    #[test]
    fn test_radix_sort_u32_empty() {
        let mut data: Vec<u32> = vec![];
        radix_sort_u32(&mut data);
        assert!(data.is_empty());
    }

    #[test]
    fn test_radix_sort_u32_single() {
        let mut data = vec![42u32];
        radix_sort_u32(&mut data);
        assert_eq!(data, vec![42]);
    }

    #[test]
    fn test_radix_sort_u32_large_values() {
        let mut data = vec![u32::MAX, 0u32, u32::MAX / 2, 1u32];
        radix_sort_u32(&mut data);
        assert_eq!(data[0], 0);
        assert_eq!(data[3], u32::MAX);
    }

    #[test]
    fn test_radix_sort_u32_duplicates() {
        let mut data = vec![3u32, 1, 3, 2, 1];
        radix_sort_u32(&mut data);
        assert_eq!(data, vec![1, 1, 2, 3, 3]);
    }

    // ── radix_sort_f32 ────────────────────────────────────────────────────────

    #[test]
    fn test_radix_sort_f32_positive_only() {
        let mut data = vec![3.0f32, 1.0, 4.0, 1.5, 0.5];
        radix_sort_f32(&mut data);
        for i in 0..data.len() - 1 {
            assert!(data[i] <= data[i + 1]);
        }
    }

    #[test]
    fn test_radix_sort_f32_with_negatives() {
        let mut data = vec![1.0f32, -2.0, 0.5, -0.5, 3.0, -1.0];
        radix_sort_f32(&mut data);
        for i in 0..data.len() - 1 {
            assert!(
                data[i] <= data[i + 1],
                "not sorted at {i}: {} > {}",
                data[i],
                data[i + 1]
            );
        }
    }

    #[test]
    fn test_radix_sort_f32_empty() {
        let mut data: Vec<f32> = vec![];
        radix_sort_f32(&mut data);
        assert!(data.is_empty());
    }

    #[test]
    fn test_f32_sort_key_roundtrip_positive() {
        let v = 3.125f32;
        assert_eq!(sort_key_to_f32(f32_to_sort_key(v)), v);
    }

    #[test]
    fn test_f32_sort_key_roundtrip_negative() {
        let v = -2.719f32;
        assert_eq!(sort_key_to_f32(f32_to_sort_key(v)), v);
    }

    // ── prefix_sum ────────────────────────────────────────────────────────────

    #[test]
    fn test_prefix_sum_basic() {
        let data = [1u32, 2, 3, 4];
        let result = prefix_sum(&data);
        assert_eq!(result, vec![0, 1, 3, 6]);
    }

    #[test]
    fn test_prefix_sum_empty() {
        let result = prefix_sum(&[]);
        assert!(result.is_empty());
    }

    #[test]
    fn test_prefix_sum_single() {
        let result = prefix_sum(&[5u32]);
        assert_eq!(result, vec![0]);
    }

    #[test]
    fn test_prefix_sum_all_ones() {
        let data = vec![1u32; 5];
        let result = prefix_sum(&data);
        assert_eq!(result, vec![0, 1, 2, 3, 4]);
    }

    // ── histogram ────────────────────────────────────────────────────────────

    #[test]
    fn test_histogram_basic() {
        let data = [0u32, 1, 2, 3, 4, 5, 6, 7, 8, 9];
        let h = histogram(&data, 5);
        assert_eq!(h.len(), 5);
        let total: u32 = h.iter().sum();
        assert_eq!(total, 10);
    }

    #[test]
    fn test_histogram_empty_data() {
        let h = histogram(&[], 4);
        assert_eq!(h, vec![0, 0, 0, 0]);
    }

    #[test]
    fn test_histogram_all_same() {
        let data = vec![5u32; 10];
        let h = histogram(&data, 3);
        assert_eq!(h.iter().sum::<u32>(), 10);
    }

    // ── counting_sort ─────────────────────────────────────────────────────────

    #[test]
    fn test_counting_sort_basic() {
        let mut data = vec![3u32, 1, 4, 1, 5, 9, 2, 6];
        counting_sort(&mut data, 9);
        for i in 0..data.len() - 1 {
            assert!(data[i] <= data[i + 1]);
        }
    }

    #[test]
    fn test_counting_sort_empty() {
        let mut data: Vec<u32> = vec![];
        counting_sort(&mut data, 10);
        assert!(data.is_empty());
    }

    #[test]
    fn test_counting_sort_single() {
        let mut data = vec![7u32];
        counting_sort(&mut data, 10);
        assert_eq!(data, vec![7]);
    }

    #[test]
    fn test_counting_sort_duplicates() {
        let mut data = vec![2u32, 2, 2, 1, 1];
        counting_sort(&mut data, 2);
        assert_eq!(data, vec![1, 1, 2, 2, 2]);
    }

    #[test]
    fn test_counting_sort_all_zero() {
        let mut data = vec![0u32; 5];
        counting_sort(&mut data, 0);
        assert_eq!(data, vec![0u32; 5]);
    }

    // ── morton_sort_3d ────────────────────────────────────────────────────────

    #[test]
    fn test_morton_sort_3d_basic() {
        let mut points = vec![
            [1.0f32, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
        ];
        morton_sort_3d(&mut points);
        // After sorting, [0,0,0] should come first (Morton code 0).
        assert_eq!(points[0], [0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_morton_sort_3d_empty() {
        let mut points: Vec<[f32; 3]> = vec![];
        morton_sort_3d(&mut points);
        assert!(points.is_empty());
    }

    #[test]
    fn test_morton_sort_3d_single() {
        let mut points = vec![[1.0f32, 2.0, 3.0]];
        morton_sort_3d(&mut points);
        assert_eq!(points, vec![[1.0, 2.0, 3.0]]);
    }

    #[test]
    fn test_morton3_origin() {
        assert_eq!(morton3(0, 0, 0), 0);
    }

    #[test]
    fn test_morton3_unit_x() {
        // x=1, y=0, z=0 => bit 0 set in x position => interleaved x is bit 0
        let code = morton3(1, 0, 0);
        assert_ne!(code, 0);
    }

    #[test]
    fn test_spread_bits_zero() {
        assert_eq!(spread_bits(0), 0);
    }

    #[test]
    fn test_spread_bits_one() {
        // bit 0 of input stays as bit 0 of output
        assert_eq!(spread_bits(1) & 1, 1);
    }

    // ── GpuSortBuffer ─────────────────────────────────────────────────────────

    #[test]
    fn test_gpu_sort_buffer_sort_pairs_basic() {
        let keys = vec![3u32, 1, 4, 1, 5, 9, 2, 6];
        let values: Vec<u32> = (0..keys.len() as u32).collect();
        let mut buf = GpuSortBuffer::new(keys, values);
        buf.sort_pairs();
        for i in 0..buf.keys.len() - 1 {
            assert!(buf.keys[i] <= buf.keys[i + 1]);
        }
    }

    #[test]
    fn test_gpu_sort_buffer_empty() {
        let mut buf = GpuSortBuffer::empty();
        buf.sort_pairs();
        assert!(buf.is_empty());
    }

    #[test]
    fn test_gpu_sort_buffer_push() {
        let mut buf = GpuSortBuffer::empty();
        buf.push(5, 100);
        buf.push(2, 200);
        assert_eq!(buf.len(), 2);
        buf.sort_pairs();
        assert_eq!(buf.keys[0], 2);
        assert_eq!(buf.values[0], 200);
    }

    #[test]
    fn test_gpu_sort_buffer_values_follow_keys() {
        let keys = vec![30u32, 10, 20];
        let values = vec![3u32, 1, 2];
        let mut buf = GpuSortBuffer::new(keys, values);
        buf.sort_pairs();
        assert_eq!(buf.keys, vec![10, 20, 30]);
        assert_eq!(buf.values, vec![1, 2, 3]);
    }

    #[test]
    fn test_gpu_sort_buffer_len_is_empty() {
        let buf = GpuSortBuffer::empty();
        assert_eq!(buf.len(), 0);
        assert!(buf.is_empty());
    }

    // ── parallel_merge ────────────────────────────────────────────────────────

    #[test]
    fn test_parallel_merge_basic() {
        let left = vec![1.0f32, 3.0, 5.0];
        let right = vec![2.0f32, 4.0, 6.0];
        let merged = parallel_merge(&left, &right);
        assert_eq!(merged, vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
    }

    #[test]
    fn test_parallel_merge_empty_left() {
        let merged = parallel_merge(&[], &[1.0f32, 2.0]);
        assert_eq!(merged, vec![1.0, 2.0]);
    }

    #[test]
    fn test_parallel_merge_empty_right() {
        let merged = parallel_merge(&[1.0f32, 2.0], &[]);
        assert_eq!(merged, vec![1.0, 2.0]);
    }

    #[test]
    fn test_parallel_merge_both_empty() {
        let merged: Vec<f32> = parallel_merge(&[], &[]);
        assert!(merged.is_empty());
    }

    #[test]
    fn test_parallel_merge_unequal_lengths() {
        let left = vec![1.0f32, 10.0];
        let right = vec![2.0f32, 3.0, 4.0, 5.0];
        let merged = parallel_merge(&left, &right);
        assert_eq!(merged.len(), 6);
        for i in 0..merged.len() - 1 {
            assert!(merged[i] <= merged[i + 1]);
        }
    }

    #[test]
    fn test_parallel_merge_duplicates() {
        let left = vec![1.0f32, 2.0, 2.0];
        let right = vec![2.0f32, 3.0];
        let merged = parallel_merge(&left, &right);
        assert_eq!(merged.len(), 5);
        for i in 0..merged.len() - 1 {
            assert!(merged[i] <= merged[i + 1]);
        }
    }
}
