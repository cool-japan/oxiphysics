// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! GPU-style grid-based reduction and aggregation kernels.
//!
//! Provides tiled parallel reductions, segmented scans, work-group min/max/sum
//! operations, warp-level primitives (simulated on CPU), stream compaction
//! (filter_compact), sparse-to-dense scatter/gather, and occupancy estimation
//! helpers.  All algorithms are CPU-side mocks that mimic GPU execution
//! semantics using Rayon.

use rayon::prelude::*;

// ---------------------------------------------------------------------------
// 1. Tile — a fixed-size work-group window
// ---------------------------------------------------------------------------

/// A single work-group tile holding up to `CAPACITY` f64 values.
#[derive(Debug, Clone)]
pub struct Tile {
    /// Elements in this tile.
    pub data: Vec<f64>,
}

impl Tile {
    /// Create a tile from a slice.
    pub fn from_slice(s: &[f64]) -> Self {
        Self { data: s.to_vec() }
    }

    /// Reduce: sum of all elements.
    pub fn reduce_sum(&self) -> f64 {
        self.data.iter().copied().sum()
    }

    /// Reduce: maximum element.
    pub fn reduce_max(&self) -> f64 {
        self.data.iter().copied().fold(f64::NEG_INFINITY, f64::max)
    }

    /// Reduce: minimum element.
    pub fn reduce_min(&self) -> f64 {
        self.data.iter().copied().fold(f64::INFINITY, f64::min)
    }

    /// Exclusive prefix scan within the tile (in-place).
    pub fn exclusive_scan_inplace(&mut self) {
        let mut acc = 0.0;
        for v in &mut self.data {
            let old = *v;
            *v = acc;
            acc += old;
        }
    }

    /// Inclusive prefix scan within the tile (in-place).
    pub fn inclusive_scan_inplace(&mut self) {
        let mut acc = 0.0;
        for v in &mut self.data {
            acc += *v;
            *v = acc;
        }
    }
}

// ---------------------------------------------------------------------------
// 2. TiledReducer — parallel multi-tile reduction
// ---------------------------------------------------------------------------

/// Parallel tiled reducer.
///
/// Splits the input into tiles of `tile_size`, performs a local reduction
/// on each tile in parallel, then reduces the per-tile results serially.
#[derive(Debug, Clone)]
pub struct TiledReducer {
    /// Number of elements per tile (analogous to GPU work-group size).
    pub tile_size: usize,
}

impl TiledReducer {
    /// Create a new tiled reducer.
    pub fn new(tile_size: usize) -> Self {
        assert!(tile_size > 0, "tile_size must be > 0");
        Self { tile_size }
    }

    /// Compute the global sum using two-level tiled reduction.
    pub fn sum(&self, data: &[f64]) -> f64 {
        if data.is_empty() {
            return 0.0;
        }
        let tile_sums: Vec<f64> = data
            .par_chunks(self.tile_size)
            .map(|chunk| chunk.iter().copied().sum::<f64>())
            .collect();
        tile_sums.iter().copied().sum()
    }

    /// Compute the global maximum.
    pub fn max(&self, data: &[f64]) -> f64 {
        if data.is_empty() {
            return f64::NEG_INFINITY;
        }
        let tile_maxs: Vec<f64> = data
            .par_chunks(self.tile_size)
            .map(|chunk| chunk.iter().copied().fold(f64::NEG_INFINITY, f64::max))
            .collect();
        tile_maxs.iter().copied().fold(f64::NEG_INFINITY, f64::max)
    }

    /// Compute the global minimum.
    pub fn min(&self, data: &[f64]) -> f64 {
        if data.is_empty() {
            return f64::INFINITY;
        }
        let tile_mins: Vec<f64> = data
            .par_chunks(self.tile_size)
            .map(|chunk| chunk.iter().copied().fold(f64::INFINITY, f64::min))
            .collect();
        tile_mins.iter().copied().fold(f64::INFINITY, f64::min)
    }

    /// Compute the global dot product of two equal-length slices.
    pub fn dot(&self, a: &[f64], b: &[f64]) -> f64 {
        assert_eq!(a.len(), b.len(), "dot product requires equal-length inputs");
        a.par_iter()
            .zip(b.par_iter())
            .map(|(&ai, &bi)| ai * bi)
            .sum()
    }

    /// Compute per-tile sums (exposes intermediate tile results).
    pub fn tile_sums(&self, data: &[f64]) -> Vec<f64> {
        data.par_chunks(self.tile_size)
            .map(|chunk| chunk.iter().copied().sum::<f64>())
            .collect()
    }
}

// ---------------------------------------------------------------------------
// 3. SegmentedScan — prefix scan within variable-length segments
// ---------------------------------------------------------------------------

/// Segmented exclusive prefix scan.
///
/// Each segment boundary is marked by a `true` entry in `flags`.  The scan
/// resets to 0 at the start of each segment.
pub fn segmented_exclusive_scan(data: &[f64], flags: &[bool]) -> Vec<f64> {
    assert_eq!(
        data.len(),
        flags.len(),
        "data and flags must be same length"
    );
    let mut result = vec![0.0; data.len()];
    let mut acc = 0.0;
    for i in 0..data.len() {
        if flags[i] {
            acc = 0.0; // new segment
        }
        result[i] = acc;
        acc += data[i];
    }
    result
}

/// Segmented inclusive prefix scan.
pub fn segmented_inclusive_scan(data: &[f64], flags: &[bool]) -> Vec<f64> {
    assert_eq!(data.len(), flags.len());
    let mut result = vec![0.0; data.len()];
    let mut acc = 0.0;
    for i in 0..data.len() {
        if flags[i] {
            acc = 0.0;
        }
        acc += data[i];
        result[i] = acc;
    }
    result
}

/// Segmented reduce: sum within each segment, returns one value per segment.
pub fn segmented_reduce_sum(data: &[f64], flags: &[bool]) -> Vec<f64> {
    assert_eq!(data.len(), flags.len());
    let mut sums: Vec<f64> = Vec::new();
    let mut acc = 0.0;
    for i in 0..data.len() {
        if flags[i] && i > 0 {
            sums.push(acc);
            acc = 0.0;
        }
        acc += data[i];
    }
    sums.push(acc);
    sums
}

// ---------------------------------------------------------------------------
// 4. Stream Compaction (filter_compact)
// ---------------------------------------------------------------------------

/// Stream compaction: collect elements satisfying `predicate` into a new vec.
///
/// Mimics GPU stream compaction (prefix-sum + scatter).
/// The output preserves the relative order of passing elements.
pub fn filter_compact<T, F>(data: &[T], predicate: F) -> Vec<T>
where
    T: Clone + Send + Sync,
    F: Fn(&T) -> bool + Sync,
{
    data.par_iter().filter(|x| predicate(x)).cloned().collect()
}

/// Partition `data` into two groups: (passing, failing) — stable order.
pub fn partition_stable<T, F>(data: &[T], predicate: F) -> (Vec<T>, Vec<T>)
where
    T: Clone,
    F: Fn(&T) -> bool,
{
    let mut pass = Vec::new();
    let mut fail = Vec::new();
    for x in data {
        if predicate(x) {
            pass.push(x.clone());
        } else {
            fail.push(x.clone());
        }
    }
    (pass, fail)
}

// ---------------------------------------------------------------------------
// 5. Scatter / Gather
// ---------------------------------------------------------------------------

/// Scatter: write `src[i]` to `dst[indices[i\]]`.
///
/// Panics if any index is out of bounds.
pub fn scatter(dst: &mut [f64], src: &[f64], indices: &[usize]) {
    assert_eq!(
        src.len(),
        indices.len(),
        "src and indices must have equal length"
    );
    for (&v, &idx) in src.iter().zip(indices.iter()) {
        dst[idx] = v;
    }
}

/// Gather: collect `src[indices[i\]]` into a new vec.
pub fn gather(src: &[f64], indices: &[usize]) -> Vec<f64> {
    indices.iter().map(|&i| src[i]).collect()
}

/// Atomic-add scatter (simulated serially): `dst[idx] += value`.
///
/// In a real GPU kernel this would use `atomicAdd`.
pub fn atomic_scatter_add(dst: &mut [f64], src: &[f64], indices: &[usize]) {
    assert_eq!(src.len(), indices.len());
    for (&v, &idx) in src.iter().zip(indices.iter()) {
        dst[idx] += v;
    }
}

// ---------------------------------------------------------------------------
// 6. Warp-level primitives (simulated on CPU as fixed-width groups)
// ---------------------------------------------------------------------------

/// Simulated warp size: number of lanes in one warp.
pub const WARP_SIZE: usize = 32;

/// Simulate a warp-level broadcast: every lane gets `lane_val[leader]`.
pub fn warp_broadcast(lanes: &[f64], leader: usize) -> Vec<f64> {
    assert!(leader < lanes.len(), "leader lane out of range");
    vec![lanes[leader]; lanes.len()]
}

/// Simulate a warp-level reduce-sum: all lanes get the total sum.
pub fn warp_reduce_sum(lanes: &[f64]) -> Vec<f64> {
    let total: f64 = lanes.iter().copied().sum();
    vec![total; lanes.len()]
}

/// Simulate a warp-level exclusive scan.
pub fn warp_exclusive_scan(lanes: &[f64]) -> Vec<f64> {
    let mut result = vec![0.0; lanes.len()];
    let mut acc = 0.0;
    for (i, &v) in lanes.iter().enumerate() {
        result[i] = acc;
        acc += v;
    }
    result
}

/// Simulate warp vote: `any` — returns true if any lane passes `pred`.
pub fn warp_vote_any<F: Fn(f64) -> bool>(lanes: &[f64], pred: F) -> bool {
    lanes.iter().any(|&v| pred(v))
}

/// Simulate warp vote: `all` — returns true if all lanes pass `pred`.
pub fn warp_vote_all<F: Fn(f64) -> bool>(lanes: &[f64], pred: F) -> bool {
    lanes.iter().all(|&v| pred(v))
}

// ---------------------------------------------------------------------------
// 7. Occupancy estimation helper
// ---------------------------------------------------------------------------

/// Compute the theoretical SM occupancy given resource usage.
///
/// Returns a value in `[0.0, 1.0]` where 1.0 = 100% occupancy.
///
/// # Parameters
/// * `wg_size`         - threads per work-group.
/// * `regs_per_thread` - registers used per thread.
/// * `shared_mem_bytes`- shared memory used per work-group.
/// * `max_wgs_per_sm`  - hardware limit (work-groups per SM).
/// * `max_threads_per_sm` - hardware limit (threads per SM).
/// * `max_regs_per_sm` - hardware limit (total registers per SM).
/// * `max_smem_per_sm` - hardware limit (shared memory bytes per SM).
pub fn estimate_occupancy(
    wg_size: usize,
    regs_per_thread: usize,
    shared_mem_bytes: usize,
    max_wgs_per_sm: usize,
    max_threads_per_sm: usize,
    max_regs_per_sm: usize,
    max_smem_per_sm: usize,
) -> f64 {
    if wg_size == 0 {
        return 0.0;
    }
    // Maximum work-groups limited by each resource.
    let by_threads = max_threads_per_sm / wg_size;
    let by_regs = if regs_per_thread == 0 {
        max_wgs_per_sm
    } else {
        max_regs_per_sm / (regs_per_thread * wg_size)
    };
    let by_smem = max_smem_per_sm
        .checked_div(shared_mem_bytes)
        .unwrap_or(max_wgs_per_sm);
    let actual_wgs = by_threads.min(by_regs).min(by_smem).min(max_wgs_per_sm);
    let active_threads = actual_wgs * wg_size;
    (active_threads as f64 / max_threads_per_sm as f64).min(1.0)
}

// ---------------------------------------------------------------------------
// 8. GridReduceStats — aggregate statistics over a 3-D grid
// ---------------------------------------------------------------------------

/// Aggregate statistics computed over a 3-D grid of f64 values.
#[derive(Debug, Clone)]
pub struct GridReduceStats {
    /// Total number of elements.
    pub count: usize,
    /// Sum of all elements.
    pub sum: f64,
    /// Mean value.
    pub mean: f64,
    /// Variance (population).
    pub variance: f64,
    /// Minimum value.
    pub min: f64,
    /// Maximum value.
    pub max: f64,
}

impl GridReduceStats {
    /// Compute statistics from a flat slice using parallel reductions.
    pub fn compute(data: &[f64]) -> Self {
        let count = data.len();
        if count == 0 {
            return Self {
                count: 0,
                sum: 0.0,
                mean: 0.0,
                variance: 0.0,
                min: 0.0,
                max: 0.0,
            };
        }
        let sum: f64 = data.par_iter().copied().sum();
        let mean = sum / count as f64;
        let variance: f64 = data
            .par_iter()
            .map(|&v| (v - mean) * (v - mean))
            .sum::<f64>()
            / count as f64;
        let min = data.par_iter().copied().reduce(|| f64::INFINITY, f64::min);
        let max = data
            .par_iter()
            .copied()
            .reduce(|| f64::NEG_INFINITY, f64::max);
        Self {
            count,
            sum,
            mean,
            variance,
            min,
            max,
        }
    }

    /// Standard deviation.
    pub fn std_dev(&self) -> f64 {
        self.variance.sqrt()
    }
}

// ---------------------------------------------------------------------------
// 9. Histogram kernel (parallel, fixed-bin-count)
// ---------------------------------------------------------------------------

/// Fixed-bin histogram over a f64 slice.
///
/// Values outside `[lo, hi)` are clamped into the boundary bins.
/// Mimics a GPU atomic histogram with one thread per element.
#[derive(Debug, Clone)]
pub struct Histogram {
    /// Bin counts.
    pub bins: Vec<u64>,
    /// Lower bound of the first bin.
    pub lo: f64,
    /// Upper bound of the last bin (exclusive).
    pub hi: f64,
}

impl Histogram {
    /// Compute a histogram with `n_bins` bins over `[lo, hi)`.
    ///
    /// Panics if `n_bins == 0` or `lo >= hi`.
    pub fn compute(data: &[f64], lo: f64, hi: f64, n_bins: usize) -> Self {
        assert!(n_bins > 0, "n_bins must be > 0");
        assert!(lo < hi, "lo must be < hi");
        let width = hi - lo;
        let mut bins = vec![0u64; n_bins];
        for &v in data {
            let idx = ((v - lo) / width * n_bins as f64) as isize;
            let idx = idx.max(0).min(n_bins as isize - 1) as usize;
            bins[idx] += 1;
        }
        Self { bins, lo, hi }
    }

    /// Total count of elements in all bins.
    pub fn total(&self) -> u64 {
        self.bins.iter().sum()
    }

    /// Centre value of bin `i`.
    pub fn bin_centre(&self, i: usize) -> f64 {
        let bin_width = (self.hi - self.lo) / self.bins.len() as f64;
        self.lo + (i as f64 + 0.5) * bin_width
    }

    /// Index of the most-populated bin (mode).
    pub fn mode_bin(&self) -> usize {
        self.bins
            .iter()
            .enumerate()
            .max_by_key(|&(_, c)| *c)
            .map(|(i, _)| i)
            .unwrap_or(0)
    }

    /// Approximate mean computed from bin centres.
    pub fn approx_mean(&self) -> f64 {
        let total = self.total();
        if total == 0 {
            return 0.0;
        }
        let sum: f64 = self
            .bins
            .iter()
            .enumerate()
            .map(|(i, &c)| self.bin_centre(i) * c as f64)
            .sum();
        sum / total as f64
    }
}

// ---------------------------------------------------------------------------
// 10. L1 / L2 / Linf norms
// ---------------------------------------------------------------------------

/// L1 norm: sum of absolute values.
pub fn norm_l1(data: &[f64]) -> f64 {
    data.par_iter().map(|&v| v.abs()).sum()
}

/// L2 (Euclidean) norm.
pub fn norm_l2(data: &[f64]) -> f64 {
    let sq: f64 = data.par_iter().map(|&v| v * v).sum();
    sq.sqrt()
}

/// L∞ (Chebyshev) norm: maximum absolute value.
pub fn norm_linf(data: &[f64]) -> f64 {
    data.par_iter()
        .map(|&v| v.abs())
        .reduce(|| 0.0_f64, f64::max)
}

/// Squared L2 distance between two equal-length vectors.
pub fn dist_sq_l2(a: &[f64], b: &[f64]) -> f64 {
    assert_eq!(a.len(), b.len());
    a.par_iter()
        .zip(b.par_iter())
        .map(|(&ai, &bi)| (ai - bi) * (ai - bi))
        .sum()
}

/// L2 distance between two equal-length vectors.
pub fn dist_l2(a: &[f64], b: &[f64]) -> f64 {
    dist_sq_l2(a, b).sqrt()
}

// ---------------------------------------------------------------------------
// 11. Covariance matrix (parallel, CPU mock)
// ---------------------------------------------------------------------------

/// Compute the `d×d` (population) covariance matrix for `n` observations of
/// dimension `d`, stored row-major in `data` (shape `n × d`).
///
/// Returns a flat `d*d` vector, row-major.
pub fn covariance_matrix(data: &[f64], n: usize, d: usize) -> Vec<f64> {
    assert_eq!(data.len(), n * d, "data must have n*d elements");
    // Mean of each dimension
    let mut mean = vec![0.0f64; d];
    for row in 0..n {
        for col in 0..d {
            mean[col] += data[row * d + col];
        }
    }
    for m in &mut mean {
        *m /= n as f64;
    }

    // Covariance: C[i][j] = E[(X_i - mean_i)(X_j - mean_j)]
    let mut cov = vec![0.0f64; d * d];
    for row in 0..n {
        for i in 0..d {
            for j in 0..d {
                let xi = data[row * d + i] - mean[i];
                let xj = data[row * d + j] - mean[j];
                cov[i * d + j] += xi * xj;
            }
        }
    }
    for c in &mut cov {
        *c /= n as f64;
    }
    cov
}

/// Extract the diagonal of a `d×d` matrix stored as a flat `d*d` slice.
pub fn matrix_diagonal(mat: &[f64], d: usize) -> Vec<f64> {
    (0..d).map(|i| mat[i * d + i]).collect()
}

// ---------------------------------------------------------------------------
// 12. Dense matrix-vector multiply (CPU mock GPU GEMV)
// ---------------------------------------------------------------------------

/// Compute `y = A * x` where `A` is `m × n` (row-major), `x` has `n`
/// elements, result `y` has `m` elements.
pub fn matvec(a: &[f64], m: usize, n: usize, x: &[f64]) -> Vec<f64> {
    assert_eq!(a.len(), m * n);
    assert_eq!(x.len(), n);
    (0..m)
        .map(|i| {
            a[i * n..(i + 1) * n]
                .iter()
                .zip(x.iter())
                .map(|(&ai, &xi)| ai * xi)
                .sum()
        })
        .collect()
}

/// Compute `C = A * B` where `A` is `m × k` and `B` is `k × n` (all row-major).
/// Returns a flat `m*n` vector.
pub fn matmul(a: &[f64], m: usize, k: usize, b: &[f64], n: usize) -> Vec<f64> {
    assert_eq!(a.len(), m * k);
    assert_eq!(b.len(), k * n);
    let mut c = vec![0.0f64; m * n];
    for i in 0..m {
        for p in 0..k {
            let a_ip = a[i * k + p];
            for j in 0..n {
                c[i * n + j] += a_ip * b[p * n + j];
            }
        }
    }
    c
}

// ---------------------------------------------------------------------------
// 13. Running statistics accumulator (Welford online algorithm)
// ---------------------------------------------------------------------------

/// Welford online statistics: accumulates count, mean, and variance in O(1)
/// per sample.  Suitable for streaming GPU readback values.
#[derive(Debug, Clone, Default)]
pub struct WelfordStats {
    /// Number of samples observed.
    pub count: u64,
    /// Current mean.
    pub mean: f64,
    /// Running M2 (sum of squared deviations from the mean).
    m2: f64,
}

impl WelfordStats {
    /// Feed a new sample.
    pub fn update(&mut self, x: f64) {
        self.count += 1;
        let delta = x - self.mean;
        self.mean += delta / self.count as f64;
        let delta2 = x - self.mean;
        self.m2 += delta * delta2;
    }

    /// Population variance.
    pub fn variance(&self) -> f64 {
        if self.count < 2 {
            return 0.0;
        }
        self.m2 / self.count as f64
    }

    /// Sample variance (Bessel-corrected).
    pub fn sample_variance(&self) -> f64 {
        if self.count < 2 {
            return 0.0;
        }
        self.m2 / (self.count - 1) as f64
    }

    /// Standard deviation (population).
    pub fn std_dev(&self) -> f64 {
        self.variance().sqrt()
    }
}

// ---------------------------------------------------------------------------
// 14. Parallel histogram reduce (multi-bin merge pattern)
// ---------------------------------------------------------------------------

/// Parallel histogram reduce: split `data` into `n_workers` chunks, compute
/// a partial histogram per chunk (in parallel), then merge all partial
/// histograms serially.  Mirrors the GPU pattern of per-work-group private
/// histograms followed by a reduction pass.
pub fn parallel_histogram(
    data: &[f64],
    lo: f64,
    hi: f64,
    n_bins: usize,
    n_workers: usize,
) -> Vec<u64> {
    assert!(n_bins > 0);
    assert!(lo < hi);
    let chunk_size = data.len().div_ceil(n_workers.max(1));
    if chunk_size == 0 {
        return vec![0u64; n_bins];
    }
    let partial: Vec<Vec<u64>> = data
        .par_chunks(chunk_size)
        .map(|chunk| {
            let width = hi - lo;
            let mut bins = vec![0u64; n_bins];
            for &v in chunk {
                let idx = ((v - lo) / width * n_bins as f64) as isize;
                let idx = idx.max(0).min(n_bins as isize - 1) as usize;
                bins[idx] += 1;
            }
            bins
        })
        .collect();

    // Merge
    let mut merged = vec![0u64; n_bins];
    for part in &partial {
        for (m, &p) in merged.iter_mut().zip(part.iter()) {
            *m += p;
        }
    }
    merged
}

// ---------------------------------------------------------------------------
// 15. Prefix-sum on integer counts (used for compaction offsets)
// ---------------------------------------------------------------------------

/// Exclusive prefix sum on a `u64` slice.  Returns a new vec.
pub fn exclusive_scan_u64(data: &[u64]) -> Vec<u64> {
    let mut result = Vec::with_capacity(data.len());
    let mut acc = 0u64;
    for &v in data {
        result.push(acc);
        acc = acc.saturating_add(v);
    }
    result
}

/// Inclusive prefix sum on a `u64` slice.
pub fn inclusive_scan_u64(data: &[u64]) -> Vec<u64> {
    let mut result = Vec::with_capacity(data.len());
    let mut acc = 0u64;
    for &v in data {
        acc = acc.saturating_add(v);
        result.push(acc);
    }
    result
}

// ---------------------------------------------------------------------------
// 16. Tile-based convolution (1-D, CPU mock)
// ---------------------------------------------------------------------------

/// 1-D convolution of `signal` with `kernel` (full output, length = signal+kernel-1).
///
/// This CPU mock mimics the tiled convolution pattern used in GPU compute
/// shaders where each work-group processes one tile with halo elements.
pub fn convolve1d(signal: &[f64], kernel: &[f64]) -> Vec<f64> {
    if signal.is_empty() || kernel.is_empty() {
        return vec![];
    }
    let out_len = signal.len() + kernel.len() - 1;
    let mut out = vec![0.0f64; out_len];
    for (i, &s) in signal.iter().enumerate() {
        for (j, &k) in kernel.iter().enumerate() {
            out[i + j] += s * k;
        }
    }
    out
}

/// 1-D cross-correlation of `signal` with `pattern` (valid region only).
/// Output length = `signal.len() - pattern.len() + 1`.
pub fn correlate1d_valid(signal: &[f64], pattern: &[f64]) -> Vec<f64> {
    if pattern.len() > signal.len() {
        return vec![];
    }
    let out_len = signal.len() - pattern.len() + 1;
    (0..out_len)
        .map(|i| {
            signal[i..i + pattern.len()]
                .iter()
                .zip(pattern.iter())
                .map(|(&s, &p)| s * p)
                .sum()
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod grid_reduce_tests {
    use super::*;
    use crate::grid_reduce::Histogram;

    use crate::grid_reduce::Tile;
    use crate::grid_reduce::TiledReducer;

    use crate::grid_reduce::WelfordStats;

    use crate::grid_reduce::exclusive_scan_u64;

    use crate::grid_reduce::inclusive_scan_u64;

    use crate::grid_reduce::segmented_reduce_sum;

    #[test]
    fn test_tile_reduce_sum() {
        let t = Tile::from_slice(&[1.0, 2.0, 3.0, 4.0]);
        assert!((t.reduce_sum() - 10.0).abs() < 1e-12);
    }

    #[test]
    fn test_tile_exclusive_scan() {
        let mut t = Tile::from_slice(&[1.0, 2.0, 3.0, 4.0]);
        t.exclusive_scan_inplace();
        assert_eq!(t.data, vec![0.0, 1.0, 3.0, 6.0]);
    }

    #[test]
    fn test_tile_inclusive_scan() {
        let mut t = Tile::from_slice(&[1.0, 2.0, 3.0]);
        t.inclusive_scan_inplace();
        assert_eq!(t.data, vec![1.0, 3.0, 6.0]);
    }

    #[test]
    fn test_tiled_reducer_sum() {
        let data: Vec<f64> = (1..=100).map(|i| i as f64).collect();
        let r = TiledReducer::new(16);
        let s = r.sum(&data);
        assert!((s - 5050.0).abs() < 1e-8, "sum 1..100 = 5050, got {s}");
    }

    #[test]
    fn test_tiled_reducer_dot_product() {
        let a = vec![1.0, 2.0, 3.0];
        let b = vec![4.0, 5.0, 6.0];
        let r = TiledReducer::new(8);
        let d = r.dot(&a, &b);
        assert!((d - 32.0).abs() < 1e-12, "dot([1,2,3],[4,5,6]) = 32");
    }

    #[test]
    fn test_segmented_exclusive_scan() {
        let data = [1.0, 2.0, 3.0, 1.0, 2.0];
        let flags = [true, false, false, true, false];
        let out = segmented_exclusive_scan(&data, &flags);
        assert_eq!(out, vec![0.0, 1.0, 3.0, 0.0, 1.0]);
    }

    #[test]
    fn test_segmented_reduce_sum() {
        let data = [1.0, 2.0, 3.0, 10.0, 20.0];
        let flags = [true, false, false, true, false];
        let sums = segmented_reduce_sum(&data, &flags);
        assert_eq!(sums.len(), 2);
        assert!((sums[0] - 6.0).abs() < 1e-12, "first segment sum = 6");
        assert!((sums[1] - 30.0).abs() < 1e-12, "second segment sum = 30");
    }

    #[test]
    fn test_filter_compact() {
        let data = vec![1.0, -2.0, 3.0, -4.0, 5.0];
        let pos: Vec<f64> = filter_compact(&data, |&x| x > 0.0);
        assert_eq!(pos, vec![1.0, 3.0, 5.0]);
    }

    #[test]
    fn test_scatter_gather_roundtrip() {
        let mut dst = vec![0.0; 5];
        let src = vec![10.0, 20.0, 30.0];
        let indices = vec![4, 1, 2];
        scatter(&mut dst, &src, &indices);
        assert!((dst[4] - 10.0).abs() < 1e-12);
        assert!((dst[1] - 20.0).abs() < 1e-12);
        let gathered = gather(&dst, &[4, 1, 2]);
        assert_eq!(gathered, vec![10.0, 20.0, 30.0]);
    }

    #[test]
    fn test_warp_reduce_sum_all_lanes_equal() {
        let lanes = vec![1.0, 2.0, 3.0, 4.0];
        let result = warp_reduce_sum(&lanes);
        assert!(
            result.iter().all(|&v| (v - 10.0).abs() < 1e-12),
            "all lanes should get the total sum"
        );
    }

    #[test]
    fn test_warp_exclusive_scan() {
        let lanes = vec![1.0, 1.0, 1.0, 1.0];
        let out = warp_exclusive_scan(&lanes);
        assert_eq!(out, vec![0.0, 1.0, 2.0, 3.0]);
    }

    #[test]
    fn test_occupancy_estimate_full() {
        // Design: 64 threads, 32 regs, 0 smem. With SM supporting 2048 threads,
        // 64 work-groups limit = 2048/64 = 32. Occupancy = 1.0.
        let occ = estimate_occupancy(64, 32, 0, 32, 2048, 65536, 49152);
        assert!((occ - 1.0).abs() < 1e-9, "should be 100% occupancy");
    }

    #[test]
    fn test_grid_reduce_stats() {
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let stats = GridReduceStats::compute(&data);
        assert_eq!(stats.count, 5);
        assert!((stats.sum - 15.0).abs() < 1e-10);
        assert!((stats.mean - 3.0).abs() < 1e-10);
        assert!((stats.min - 1.0).abs() < 1e-10);
        assert!((stats.max - 5.0).abs() < 1e-10);
        // Variance = E[(X-mean)^2] = (4+1+0+1+4)/5 = 2.0
        assert!((stats.variance - 2.0).abs() < 1e-10);
        assert!((stats.std_dev() - 2.0_f64.sqrt()).abs() < 1e-10);
    }

    #[test]
    fn test_grid_reduce_stats_empty() {
        let stats = GridReduceStats::compute(&[]);
        assert_eq!(stats.count, 0);
        assert!((stats.sum).abs() < 1e-12);
    }

    // ── Histogram tests ──────────────────────────────────────────────────

    #[test]
    fn test_histogram_basic() {
        let data = vec![0.1, 0.5, 0.9, 1.5, 1.9];
        let h = Histogram::compute(&data, 0.0, 2.0, 2);
        // bin 0: [0,1) => 0.1, 0.5, 0.9  count=3
        // bin 1: [1,2) => 1.5, 1.9        count=2
        assert_eq!(h.bins[0], 3);
        assert_eq!(h.bins[1], 2);
        assert_eq!(h.total(), 5);
    }

    #[test]
    fn test_histogram_mode_bin() {
        let data = vec![0.1, 0.2, 0.3, 1.5];
        let h = Histogram::compute(&data, 0.0, 2.0, 2);
        assert_eq!(h.mode_bin(), 0); // bin 0 has 3 elements
    }

    #[test]
    fn test_histogram_bin_centre() {
        let h = Histogram::compute(&[], 0.0, 4.0, 4);
        // each bin width = 1.0, centre of bin 0 = 0.5
        assert!((h.bin_centre(0) - 0.5).abs() < 1e-10);
        assert!((h.bin_centre(3) - 3.5).abs() < 1e-10);
    }

    #[test]
    fn test_histogram_approx_mean() {
        // All data in one bin centred at 0.5
        let data = vec![0.1, 0.2, 0.3, 0.4];
        let h = Histogram::compute(&data, 0.0, 1.0, 1);
        assert!((h.approx_mean() - 0.5).abs() < 1e-10);
    }

    // ── Norm tests ───────────────────────────────────────────────────────

    #[test]
    fn test_norm_l1() {
        let v = vec![1.0, -2.0, 3.0];
        assert!((norm_l1(&v) - 6.0).abs() < 1e-12);
    }

    #[test]
    fn test_norm_l2() {
        let v = vec![3.0, 4.0];
        assert!((norm_l2(&v) - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_norm_linf() {
        let v = vec![1.0, -5.0, 3.0];
        assert!((norm_linf(&v) - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_dist_l2() {
        let a = vec![0.0, 0.0];
        let b = vec![3.0, 4.0];
        assert!((dist_l2(&a, &b) - 5.0).abs() < 1e-12);
    }

    // ── Covariance tests ─────────────────────────────────────────────────

    #[test]
    fn test_covariance_identity_pattern() {
        // Two variables, perfectly correlated: data = [(0,0),(1,1),(2,2)]
        // Cov = [[var_x, cov_xy],[cov_yx, var_y]]
        let data = vec![0.0, 0.0, 1.0, 1.0, 2.0, 2.0];
        let cov = covariance_matrix(&data, 3, 2);
        // variance of x = variance of y = 2/3; covariance = 2/3
        let expected_var = 2.0 / 3.0;
        assert!(
            (cov[0] - expected_var).abs() < 1e-10,
            "cov[0,0] = {}",
            cov[0]
        );
        assert!(
            (cov[1] - expected_var).abs() < 1e-10,
            "cov[0,1] = {}",
            cov[1]
        );
        assert!(
            (cov[3] - expected_var).abs() < 1e-10,
            "cov[1,1] = {}",
            cov[3]
        );
    }

    #[test]
    fn test_matrix_diagonal() {
        let mat = vec![1.0, 2.0, 3.0, 4.0]; // 2×2
        let diag = matrix_diagonal(&mat, 2);
        assert_eq!(diag, vec![1.0, 4.0]);
    }

    // ── GEMV / GEMM tests ────────────────────────────────────────────────

    #[test]
    fn test_matvec_identity() {
        let identity = vec![1.0, 0.0, 0.0, 1.0]; // 2×2 identity
        let x = vec![3.0, 7.0];
        let y = matvec(&identity, 2, 2, &x);
        assert_eq!(y, x);
    }

    #[test]
    fn test_matvec_basic() {
        // A = [[1,2],[3,4]], x = [1,1] => y = [3,7]
        let a = vec![1.0, 2.0, 3.0, 4.0];
        let x = vec![1.0, 1.0];
        let y = matvec(&a, 2, 2, &x);
        assert!((y[0] - 3.0).abs() < 1e-12);
        assert!((y[1] - 7.0).abs() < 1e-12);
    }

    #[test]
    fn test_matmul_2x2() {
        // A = [[1,2],[3,4]], B = [[5,6],[7,8]]
        // C = [[1*5+2*7, 1*6+2*8],[3*5+4*7, 3*6+4*8]] = [[19,22],[43,50]]
        let a = vec![1.0, 2.0, 3.0, 4.0];
        let b = vec![5.0, 6.0, 7.0, 8.0];
        let c = matmul(&a, 2, 2, &b, 2);
        assert!((c[0] - 19.0).abs() < 1e-12);
        assert!((c[1] - 22.0).abs() < 1e-12);
        assert!((c[2] - 43.0).abs() < 1e-12);
        assert!((c[3] - 50.0).abs() < 1e-12);
    }

    // ── WelfordStats tests ───────────────────────────────────────────────

    #[test]
    fn test_welford_mean_and_variance() {
        let mut w = WelfordStats::default();
        for &v in &[2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0] {
            w.update(v);
        }
        assert!((w.mean - 5.0).abs() < 1e-10, "mean = {}", w.mean);
        // population variance = 4.0
        assert!((w.variance() - 4.0).abs() < 1e-10, "var = {}", w.variance());
    }

    #[test]
    fn test_welford_single_sample() {
        let mut w = WelfordStats::default();
        w.update(42.0);
        assert!((w.mean - 42.0).abs() < 1e-12);
        assert!((w.variance()).abs() < 1e-12);
    }

    // ── Parallel histogram tests ─────────────────────────────────────────

    #[test]
    fn test_parallel_histogram_matches_serial() {
        let data: Vec<f64> = (0..200).map(|i| i as f64 / 10.0).collect(); // 0.0 .. 19.9
        let serial = Histogram::compute(&data, 0.0, 20.0, 10);
        let par = parallel_histogram(&data, 0.0, 20.0, 10, 4);
        assert_eq!(
            serial.bins, par,
            "parallel and serial histograms must agree"
        );
    }

    // ── Integer scan tests ───────────────────────────────────────────────

    #[test]
    fn test_exclusive_scan_u64() {
        let data = [1u64, 2, 3, 4];
        let out = exclusive_scan_u64(&data);
        assert_eq!(out, vec![0, 1, 3, 6]);
    }

    #[test]
    fn test_inclusive_scan_u64() {
        let data = [1u64, 2, 3, 4];
        let out = inclusive_scan_u64(&data);
        assert_eq!(out, vec![1, 3, 6, 10]);
    }

    // ── Convolution tests ────────────────────────────────────────────────

    #[test]
    fn test_convolve1d_basic() {
        // [1,2,3] * [0,1,0] = [1,2,3] (identity-ish kernel, padded)
        let sig = vec![1.0, 2.0, 3.0];
        let ker = vec![1.0];
        let out = convolve1d(&sig, &ker);
        assert_eq!(out, sig);
    }

    #[test]
    fn test_convolve1d_box_filter() {
        // Box filter [1/3, 1/3, 1/3] on [1,2,3,4,5]
        let sig = vec![0.0, 6.0, 0.0]; // impulse at centre → box response
        let ker = vec![1.0, 1.0, 1.0];
        let out = convolve1d(&sig, &ker); // len = 5
        // out = [0, 6, 6, 6, 0]
        assert!((out[0]).abs() < 1e-12);
        assert!((out[1] - 6.0).abs() < 1e-12);
        assert!((out[3] - 6.0).abs() < 1e-12);
        assert!((out[4]).abs() < 1e-12);
    }

    #[test]
    fn test_correlate1d_valid() {
        let sig = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let pat = vec![1.0, 0.0, -1.0]; // difference kernel
        let out = correlate1d_valid(&sig, &pat);
        // [1-3, 2-4, 3-5] = [-2, -2, -2]
        assert_eq!(out.len(), 3);
        assert!((out[0] - (1.0 - 3.0)).abs() < 1e-12);
        assert!((out[1] - (2.0 - 4.0)).abs() < 1e-12);
        assert!((out[2] - (3.0 - 5.0)).abs() < 1e-12);
    }
}

// ---------------------------------------------------------------------------
// 17. Blelloch parallel prefix scan (work-efficient)
// ---------------------------------------------------------------------------

/// Blelloch work-efficient parallel prefix scan (exclusive, in-place).
///
/// Implements the classic two-phase up-sweep / down-sweep algorithm.
/// Operates on a power-of-two sized buffer; the input is padded with zeros
/// if its length is not a power of two.
///
/// # Reference
/// Blelloch, G. E. (1990). *Prefix sums and their applications*.
pub fn blelloch_exclusive_scan(data: &[f64]) -> Vec<f64> {
    if data.is_empty() {
        return vec![];
    }
    // Pad to next power of two
    let n = data.len();
    let mut p = 1usize;
    while p < n {
        p <<= 1;
    }
    let mut buf = vec![0.0f64; p];
    buf[..n].copy_from_slice(data);

    // Up-sweep (reduce)
    let mut stride = 1usize;
    while stride < p {
        let step = stride * 2;
        let mut i = step - 1;
        while i < p {
            buf[i] += buf[i - stride];
            i += step;
        }
        stride = step;
    }

    // Clear the last element (identity for +)
    buf[p - 1] = 0.0;

    // Down-sweep
    let mut stride = p / 2;
    while stride >= 1 {
        let step = stride * 2;
        let mut i = step - 1;
        while i < p {
            let t = buf[i - stride];
            buf[i - stride] = buf[i];
            buf[i] += t;
            i += step;
        }
        stride /= 2;
    }

    buf[..n].to_vec()
}

/// Blelloch inclusive scan: built on top of the exclusive scan.
pub fn blelloch_inclusive_scan(data: &[f64]) -> Vec<f64> {
    let excl = blelloch_exclusive_scan(data);
    excl.into_iter()
        .zip(data.iter())
        .map(|(e, &v)| e + v)
        .collect()
}

// ---------------------------------------------------------------------------
// 18. Segmented scan (parallel, Blelloch-style)
// ---------------------------------------------------------------------------

/// Segmented exclusive scan using a parallel Blelloch-style approach.
///
/// `flags[i] == true` marks the start of a new segment.  Within each segment
/// the exclusive prefix sum is computed independently.
pub fn blelloch_segmented_exclusive_scan(data: &[f64], flags: &[bool]) -> Vec<f64> {
    assert_eq!(data.len(), flags.len());
    // We use a simple serial approach that is semantically equivalent to the
    // parallel GPU version with predicate propagation.
    segmented_exclusive_scan(data, flags)
}

/// Segmented reduce: parallel version using Rayon.
///
/// Returns one aggregate per segment.
pub fn parallel_segmented_reduce_sum(data: &[f64], flags: &[bool]) -> Vec<f64> {
    assert_eq!(data.len(), flags.len());
    // Build segment boundaries
    let mut starts = vec![0usize];
    for (i, &flag) in flags.iter().enumerate().skip(1) {
        if flag {
            starts.push(i);
        }
    }
    starts.push(data.len());
    starts
        .windows(2)
        .map(|w| data[w[0]..w[1]].iter().sum())
        .collect()
}

// ---------------------------------------------------------------------------
// 19. Stream compaction with index output
// ---------------------------------------------------------------------------

/// Stream compaction with index tracking.
///
/// Returns `(compacted_values, original_indices)` — the values that pass
/// `predicate` and the indices they came from.
pub fn filter_compact_indexed(
    data: &[f64],
    predicate: impl Fn(f64) -> bool,
) -> (Vec<f64>, Vec<usize>) {
    let mut vals = Vec::new();
    let mut idxs = Vec::new();
    for (i, &v) in data.iter().enumerate() {
        if predicate(v) {
            vals.push(v);
            idxs.push(i);
        }
    }
    (vals, idxs)
}

/// Stream compaction with count: returns (compacted, n_removed).
pub fn filter_compact_counted<T: Clone>(
    data: &[T],
    predicate: impl Fn(&T) -> bool,
) -> (Vec<T>, usize) {
    let compacted: Vec<T> = data.iter().filter(|x| predicate(x)).cloned().collect();
    let n_removed = data.len() - compacted.len();
    (compacted, n_removed)
}

// ---------------------------------------------------------------------------
// 20. Radix sort step (single digit / pass)
// ---------------------------------------------------------------------------

/// Single-pass radix sort step: sort `data` by a single `bit_offset`-wide
/// digit extracted at bit position `bit_pos` with radix `radix` (must be a
/// power of two, e.g. 256 for 8-bit digits).
///
/// Returns a new sorted vec.  `key_fn` maps each element to its sort key.
pub fn radix_sort_pass_u64(data: &[u64], bit_pos: u32, radix: usize) -> Vec<u64> {
    assert!(radix.is_power_of_two(), "radix must be a power of two");
    let mask = (radix - 1) as u64;
    // Count
    let mut counts = vec![0usize; radix];
    for &v in data {
        let digit = ((v >> bit_pos) & mask) as usize;
        counts[digit] += 1;
    }
    // Exclusive prefix sum of counts
    let offsets = exclusive_scan_u64(&counts.iter().map(|&c| c as u64).collect::<Vec<_>>());
    let mut offsets: Vec<usize> = offsets.iter().map(|&o| o as usize).collect();
    // Scatter
    let mut out = vec![0u64; data.len()];
    for &v in data {
        let digit = ((v >> bit_pos) & mask) as usize;
        out[offsets[digit]] = v;
        offsets[digit] += 1;
    }
    out
}

/// Full 64-bit radix sort (8 passes of 8-bit digits).
pub fn radix_sort_u64(data: &[u64]) -> Vec<u64> {
    let mut buf = data.to_vec();
    for pass in 0..8u32 {
        buf = radix_sort_pass_u64(&buf, pass * 8, 256);
    }
    buf
}

/// Radix sort for f64 values (sorts by bit representation, handles sign bit).
///
/// Uses the standard trick of flipping the sign bit (and all bits for negative
/// numbers) so that the radix sort on the bit pattern produces correct order.
pub fn radix_sort_f64(data: &[f64]) -> Vec<f64> {
    let mut keys: Vec<u64> = data
        .iter()
        .map(|&v| {
            let bits = v.to_bits();
            if bits >> 63 == 0 {
                bits | (1u64 << 63) // positive: flip sign bit
            } else {
                !bits // negative: flip all bits
            }
        })
        .collect();
    keys = radix_sort_u64(&keys);
    keys.iter()
        .map(|&bits| {
            let recovered = if bits >> 63 == 1 {
                bits ^ (1u64 << 63) // was positive
            } else {
                !bits // was negative
            };
            f64::from_bits(recovered)
        })
        .collect()
}

// ---------------------------------------------------------------------------
// 21. Work-efficient parallel reduce (tree reduction)
// ---------------------------------------------------------------------------

/// Work-efficient tree reduction: sums `data` using a binary tree pattern.
///
/// This simulates the GPU tree-reduction kernel where each thread handles one
/// element and the active thread count halves each step.
pub fn tree_reduce_sum(data: &[f64]) -> f64 {
    if data.is_empty() {
        return 0.0;
    }
    let mut buf = data.to_vec();
    let mut n = buf.len();
    while n > 1 {
        let half = n / 2;
        for i in 0..half {
            buf[i] += buf[i + half];
        }
        if n % 2 == 1 {
            buf[half - 1] += buf[n - 1];
        }
        n = half;
    }
    buf[0]
}

/// Work-efficient tree reduction for max.
pub fn tree_reduce_max(data: &[f64]) -> f64 {
    if data.is_empty() {
        return f64::NEG_INFINITY;
    }
    let mut buf = data.to_vec();
    let mut n = buf.len();
    while n > 1 {
        let half = n / 2;
        for i in 0..half {
            buf[i] = f64::max(buf[i], buf[i + half]);
        }
        if n % 2 == 1 {
            buf[half - 1] = f64::max(buf[half - 1], buf[n - 1]);
        }
        n = half;
    }
    buf[0]
}

/// Work-efficient tree reduction for min.
pub fn tree_reduce_min(data: &[f64]) -> f64 {
    if data.is_empty() {
        return f64::INFINITY;
    }
    let mut buf = data.to_vec();
    let mut n = buf.len();
    while n > 1 {
        let half = n / 2;
        for i in 0..half {
            buf[i] = f64::min(buf[i], buf[i + half]);
        }
        if n % 2 == 1 {
            buf[half - 1] = f64::min(buf[half - 1], buf[n - 1]);
        }
        n = half;
    }
    buf[0]
}

// ---------------------------------------------------------------------------
// 22. Reduce-then-broadcast (GPU idiom)
// ---------------------------------------------------------------------------

/// Reduce `data` to a scalar and broadcast the result back to all positions.
///
/// Mimics the GPU pattern: reduce in shared memory → broadcast from lane 0.
pub fn reduce_broadcast(data: &[f64]) -> Vec<f64> {
    let total: f64 = data.iter().copied().sum();
    vec![total; data.len()]
}

/// Normalise: divide each element by the total sum.
pub fn normalise_by_sum(data: &[f64]) -> Vec<f64> {
    let s: f64 = data.iter().copied().sum();
    if s.abs() < 1e-30 {
        return data.to_vec();
    }
    data.iter().map(|&v| v / s).collect()
}

// ---------------------------------------------------------------------------
// 23. Multi-level histogram reduce
// ---------------------------------------------------------------------------

/// Two-level histogram: first pass per-tile, second pass merge.
///
/// Returns the merged bin counts.
#[derive(Debug, Clone)]
pub struct TwoLevelHistogram {
    /// Merged histogram bins.
    pub bins: Vec<u64>,
    /// Lower bound of the value range.
    pub lo: f64,
    /// Upper bound of the value range.
    pub hi: f64,
    /// Number of work-groups / tiles used.
    pub n_tiles: usize,
}

impl TwoLevelHistogram {
    /// Compute a two-level histogram.
    pub fn compute(data: &[f64], lo: f64, hi: f64, n_bins: usize, tile_size: usize) -> Self {
        let n_tiles = data.len().div_ceil(tile_size.max(1));
        let bins = parallel_histogram(data, lo, hi, n_bins, n_tiles.max(1));
        Self {
            bins,
            lo,
            hi,
            n_tiles,
        }
    }

    /// Total count of elements.
    pub fn total(&self) -> u64 {
        self.bins.iter().sum()
    }

    /// Compute the approximate median from bin centres.
    pub fn approx_median(&self) -> f64 {
        let total = self.total();
        if total == 0 {
            return (self.lo + self.hi) / 2.0;
        }
        let half = total / 2;
        let n = self.bins.len() as f64;
        let mut acc = 0u64;
        for (i, &c) in self.bins.iter().enumerate() {
            acc += c;
            if acc >= half {
                let bin_width = (self.hi - self.lo) / n;
                return self.lo + (i as f64 + 0.5) * bin_width;
            }
        }
        self.hi
    }
}

// ---------------------------------------------------------------------------
// 24. Running min/max tracker (streaming GPU readback)
// ---------------------------------------------------------------------------

/// Streaming min/max tracker suitable for GPU readback values.
#[derive(Debug, Clone, Default)]
pub struct RunningMinMax {
    /// Current minimum.
    pub min: f64,
    /// Current maximum.
    pub max: f64,
    /// Number of samples observed.
    pub count: u64,
}

impl RunningMinMax {
    /// Create a new tracker.
    pub fn new() -> Self {
        Self {
            min: f64::INFINITY,
            max: f64::NEG_INFINITY,
            count: 0,
        }
    }

    /// Update with a new sample.
    pub fn update(&mut self, v: f64) {
        self.min = f64::min(self.min, v);
        self.max = f64::max(self.max, v);
        self.count += 1;
    }

    /// Update with a batch of samples.
    pub fn update_slice(&mut self, data: &[f64]) {
        for &v in data {
            self.update(v);
        }
    }

    /// Range (max - min).
    pub fn range(&self) -> f64 {
        if self.count == 0 {
            return 0.0;
        }
        self.max - self.min
    }
}

// ---------------------------------------------------------------------------
// 25. Compact scatter (GPU stream compaction output pattern)
// ---------------------------------------------------------------------------

/// Compact scatter: given a predicate mask, scatter `src` elements into a
/// destination buffer at compacted positions.
///
/// Returns the number of elements written.
pub fn compact_scatter(src: &[f64], mask: &[bool], dst: &mut Vec<f64>) -> usize {
    assert_eq!(src.len(), mask.len());
    let before = dst.len();
    for (&v, &keep) in src.iter().zip(mask.iter()) {
        if keep {
            dst.push(v);
        }
    }
    dst.len() - before
}

/// Build a compaction offset table from a boolean mask.
///
/// Returns a vec of length `mask.len()` where `offsets[i]` is the compacted
/// index for element `i`, or `usize::MAX` when `mask[i]` is false.
pub fn compaction_offsets(mask: &[bool]) -> Vec<usize> {
    let mut result = vec![usize::MAX; mask.len()];
    let mut counter = 0usize;
    for (i, &keep) in mask.iter().enumerate() {
        if keep {
            result[i] = counter;
            counter += 1;
        }
    }
    result
}

// ---------------------------------------------------------------------------
// Tests — new additions
// ---------------------------------------------------------------------------

#[cfg(test)]
mod extended_tests {
    use crate::grid_reduce::Histogram;
    use crate::grid_reduce::RunningMinMax;
    use crate::grid_reduce::Tile;
    use crate::grid_reduce::TiledReducer;
    use crate::grid_reduce::TwoLevelHistogram;
    use crate::grid_reduce::WelfordStats;
    use crate::grid_reduce::blelloch_exclusive_scan;
    use crate::grid_reduce::blelloch_inclusive_scan;
    use crate::grid_reduce::compact_scatter;
    use crate::grid_reduce::compaction_offsets;
    use crate::grid_reduce::exclusive_scan_u64;
    use crate::grid_reduce::filter_compact_counted;
    use crate::grid_reduce::filter_compact_indexed;
    use crate::grid_reduce::inclusive_scan_u64;
    use crate::grid_reduce::normalise_by_sum;
    use crate::grid_reduce::parallel_segmented_reduce_sum;
    use crate::grid_reduce::radix_sort_f64;
    use crate::grid_reduce::radix_sort_pass_u64;
    use crate::grid_reduce::radix_sort_u64;
    use crate::grid_reduce::reduce_broadcast;
    use crate::grid_reduce::segmented_reduce_sum;
    use crate::grid_reduce::tree_reduce_max;
    use crate::grid_reduce::tree_reduce_min;
    use crate::grid_reduce::tree_reduce_sum;

    // ── Blelloch scan ────────────────────────────────────────────────────

    #[test]
    fn blelloch_exclusive_scan_matches_serial() {
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
        let serial = {
            let mut r = Vec::new();
            let mut acc = 0.0f64;
            for &v in &data {
                r.push(acc);
                acc += v;
            }
            r
        };
        let blelloch = blelloch_exclusive_scan(&data);
        for (a, b) in serial.iter().zip(blelloch.iter()) {
            assert!((a - b).abs() < 1e-10, "mismatch: serial={a} blelloch={b}");
        }
    }

    #[test]
    fn blelloch_exclusive_scan_non_pow2() {
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0]; // len=5, not power of 2
        let result = blelloch_exclusive_scan(&data);
        assert_eq!(result.len(), 5);
        assert!((result[0] - 0.0).abs() < 1e-10);
        assert!((result[1] - 1.0).abs() < 1e-10);
        assert!((result[2] - 3.0).abs() < 1e-10);
        assert!((result[3] - 6.0).abs() < 1e-10);
        assert!((result[4] - 10.0).abs() < 1e-10);
    }

    #[test]
    fn blelloch_inclusive_scan_correct() {
        let data = vec![1.0, 2.0, 3.0, 4.0];
        let result = blelloch_inclusive_scan(&data);
        assert_eq!(result, vec![1.0, 3.0, 6.0, 10.0]);
    }

    #[test]
    fn blelloch_exclusive_scan_single_element() {
        let result = blelloch_exclusive_scan(&[42.0]);
        assert_eq!(result, vec![0.0]);
    }

    #[test]
    fn blelloch_exclusive_scan_all_zeros() {
        let data = vec![0.0; 8];
        let result = blelloch_exclusive_scan(&data);
        assert!(result.iter().all(|&v| v.abs() < 1e-12));
    }

    // ── Segmented parallel reduce ────────────────────────────────────────

    #[test]
    fn parallel_segmented_reduce_matches_serial() {
        let data = [1.0, 2.0, 3.0, 10.0, 20.0, 30.0];
        let flags = [true, false, false, true, false, false];
        let par = parallel_segmented_reduce_sum(&data, &flags);
        let ser = segmented_reduce_sum(&data, &flags);
        assert_eq!(par, ser);
    }

    #[test]
    fn parallel_segmented_reduce_single_segment() {
        let data = [1.0, 2.0, 3.0];
        let flags = [true, false, false];
        let result = parallel_segmented_reduce_sum(&data, &flags);
        assert_eq!(result.len(), 1);
        assert!((result[0] - 6.0).abs() < 1e-10);
    }

    // ── filter_compact_indexed ───────────────────────────────────────────

    #[test]
    fn filter_compact_indexed_positive() {
        let data = vec![-1.0, 2.0, -3.0, 4.0, 5.0];
        let (vals, idxs) = filter_compact_indexed(&data, |v| v > 0.0);
        assert_eq!(vals, vec![2.0, 4.0, 5.0]);
        assert_eq!(idxs, vec![1, 3, 4]);
    }

    #[test]
    fn filter_compact_indexed_empty_result() {
        let data = vec![-1.0, -2.0, -3.0];
        let (vals, idxs) = filter_compact_indexed(&data, |v| v > 0.0);
        assert!(vals.is_empty());
        assert!(idxs.is_empty());
    }

    #[test]
    fn filter_compact_counted_removes_negatives() {
        let data = vec![1.0, -2.0, 3.0, -4.0, 5.0];
        let (kept, removed) = filter_compact_counted(&data, |v| *v >= 0.0);
        assert_eq!(kept, vec![1.0, 3.0, 5.0]);
        assert_eq!(removed, 2);
    }

    // ── Radix sort ───────────────────────────────────────────────────────

    #[test]
    fn radix_sort_u64_ascending() {
        let mut data = vec![5u64, 3, 8, 1, 9, 2, 7, 4, 6, 0];
        let sorted = radix_sort_u64(&data);
        data.sort_unstable();
        assert_eq!(sorted, data);
    }

    #[test]
    fn radix_sort_u64_empty() {
        let sorted = radix_sort_u64(&[]);
        assert!(sorted.is_empty());
    }

    #[test]
    fn radix_sort_u64_already_sorted() {
        let data = vec![1u64, 2, 3, 4, 5];
        assert_eq!(radix_sort_u64(&data), data);
    }

    #[test]
    fn radix_sort_u64_reverse() {
        let data = vec![5u64, 4, 3, 2, 1];
        let sorted = radix_sort_u64(&data);
        assert_eq!(sorted, vec![1u64, 2, 3, 4, 5]);
    }

    #[test]
    fn radix_sort_f64_positive_values() {
        let data = vec![3.125, 1.41, 2.71, 0.57, 1.73];
        let sorted = radix_sort_f64(&data);
        let mut expected = data.clone();
        expected.sort_by(|a, b| a.partial_cmp(b).unwrap());
        for (a, b) in sorted.iter().zip(expected.iter()) {
            assert!((a - b).abs() < 1e-12, "a={a} b={b}");
        }
    }

    #[test]
    fn radix_sort_pass_u64_single_pass() {
        // Sort by lowest byte
        let data = vec![0x03u64, 0x01, 0x04, 0x01, 0x05];
        let sorted = radix_sort_pass_u64(&data, 0, 256);
        assert_eq!(sorted.len(), data.len());
        // Lowest bytes should be non-decreasing
        for w in sorted.windows(2) {
            assert!(w[0] & 0xFF <= w[1] & 0xFF, "not sorted by low byte");
        }
    }

    // ── Tree reduce ──────────────────────────────────────────────────────

    #[test]
    fn tree_reduce_sum_correct() {
        let data: Vec<f64> = (1..=16).map(|i| i as f64).collect();
        let s = tree_reduce_sum(&data);
        assert!((s - 136.0).abs() < 1e-10, "sum = {s}");
    }

    #[test]
    fn tree_reduce_sum_odd_length() {
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let s = tree_reduce_sum(&data);
        assert!((s - 15.0).abs() < 1e-10, "sum = {s}");
    }

    #[test]
    fn tree_reduce_max_correct() {
        let data = vec![3.0, 1.0, 4.0, 1.0, 5.0, 9.0, 2.0, 6.0];
        assert!((tree_reduce_max(&data) - 9.0).abs() < 1e-12);
    }

    #[test]
    fn tree_reduce_min_correct() {
        let data = vec![3.0, 1.0, 4.0, 1.0, 5.0, 9.0, 2.0, 6.0];
        assert!((tree_reduce_min(&data) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn tree_reduce_empty() {
        assert!((tree_reduce_sum(&[])).abs() < 1e-12);
        assert!(tree_reduce_max(&[]) == f64::NEG_INFINITY);
        assert!(tree_reduce_min(&[]) == f64::INFINITY);
    }

    #[test]
    fn tree_reduce_single() {
        assert!((tree_reduce_sum(&[42.0]) - 42.0).abs() < 1e-12);
        assert!((tree_reduce_max(&[42.0]) - 42.0).abs() < 1e-12);
        assert!((tree_reduce_min(&[42.0]) - 42.0).abs() < 1e-12);
    }

    #[test]
    fn tree_reduce_matches_tiled_reducer() {
        let data: Vec<f64> = (0..100).map(|i| i as f64).collect();
        let tr = TiledReducer::new(16);
        let tiled_sum = tr.sum(&data);
        let tree_sum = tree_reduce_sum(&data);
        assert!(
            (tiled_sum - tree_sum).abs() < 1e-8,
            "tiled={tiled_sum} tree={tree_sum}"
        );
    }

    // ── Reduce + broadcast ───────────────────────────────────────────────

    #[test]
    fn reduce_broadcast_all_equal() {
        let data = vec![1.0, 2.0, 3.0];
        let result = reduce_broadcast(&data);
        assert!(
            result.iter().all(|&v| (v - 6.0).abs() < 1e-12),
            "all should equal 6"
        );
    }

    #[test]
    fn normalise_by_sum_sums_to_one() {
        let data = vec![1.0, 2.0, 3.0, 4.0];
        let normed = normalise_by_sum(&data);
        let s: f64 = normed.iter().sum();
        assert!((s - 1.0).abs() < 1e-10, "sum = {s}");
    }

    #[test]
    fn normalise_by_sum_zero_input_unchanged() {
        let data = vec![0.0, 0.0, 0.0];
        let result = normalise_by_sum(&data);
        assert_eq!(result, data);
    }

    // ── TwoLevelHistogram ────────────────────────────────────────────────

    #[test]
    fn two_level_histogram_total_correct() {
        let data: Vec<f64> = (0..100).map(|i| i as f64 / 10.0).collect();
        let h = TwoLevelHistogram::compute(&data, 0.0, 10.0, 10, 16);
        assert_eq!(h.total(), 100);
    }

    #[test]
    fn two_level_histogram_approx_median() {
        // Uniform [0,10) → median ≈ 5.0
        let data: Vec<f64> = (0..1000).map(|i| i as f64 / 100.0).collect();
        let h = TwoLevelHistogram::compute(&data, 0.0, 10.0, 100, 64);
        let med = h.approx_median();
        assert!((med - 5.0).abs() < 0.2, "approx median = {med}");
    }

    #[test]
    fn two_level_histogram_bins_count_matches() {
        let data = vec![0.5, 1.5, 2.5, 3.5];
        let h = TwoLevelHistogram::compute(&data, 0.0, 4.0, 4, 2);
        assert_eq!(h.total(), 4);
        for &c in &h.bins {
            assert_eq!(c, 1, "each bin should have 1 element");
        }
    }

    // ── RunningMinMax ────────────────────────────────────────────────────

    #[test]
    fn running_min_max_basic() {
        let mut t = RunningMinMax::new();
        t.update_slice(&[3.0, 1.0, 4.0, 1.0, 5.0]);
        assert!((t.min - 1.0).abs() < 1e-12);
        assert!((t.max - 5.0).abs() < 1e-12);
        assert_eq!(t.count, 5);
        assert!((t.range() - 4.0).abs() < 1e-12);
    }

    #[test]
    fn running_min_max_single() {
        let mut t = RunningMinMax::new();
        t.update(42.0);
        assert!((t.min - 42.0).abs() < 1e-12);
        assert!((t.max - 42.0).abs() < 1e-12);
        assert!((t.range()).abs() < 1e-12);
    }

    #[test]
    fn running_min_max_empty_range() {
        let t = RunningMinMax::new();
        assert!((t.range()).abs() < 1e-12);
    }

    // ── compact_scatter / compaction_offsets ─────────────────────────────

    #[test]
    fn compact_scatter_basic() {
        let src = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let mask = vec![true, false, true, false, true];
        let mut dst = Vec::new();
        let n = compact_scatter(&src, &mask, &mut dst);
        assert_eq!(n, 3);
        assert_eq!(dst, vec![1.0, 3.0, 5.0]);
    }

    #[test]
    fn compact_scatter_appends_to_existing() {
        let src = vec![10.0, 20.0];
        let mask = vec![true, true];
        let mut dst = vec![0.0, 0.0];
        compact_scatter(&src, &mask, &mut dst);
        assert_eq!(dst, vec![0.0, 0.0, 10.0, 20.0]);
    }

    #[test]
    fn compaction_offsets_correct() {
        let mask = vec![true, false, true, false, true];
        let offsets = compaction_offsets(&mask);
        assert_eq!(offsets[0], 0);
        assert_eq!(offsets[1], usize::MAX);
        assert_eq!(offsets[2], 1);
        assert_eq!(offsets[3], usize::MAX);
        assert_eq!(offsets[4], 2);
    }

    #[test]
    fn compaction_offsets_all_false() {
        let mask = vec![false; 5];
        let offsets = compaction_offsets(&mask);
        assert!(offsets.iter().all(|&o| o == usize::MAX));
    }

    // ── Additional histogram ──────────────────────────────────────────────

    #[test]
    fn histogram_uniform_distribution() {
        let data: Vec<f64> = (0..10).map(|i| i as f64 + 0.5).collect();
        let h = Histogram::compute(&data, 0.0, 10.0, 10);
        for &c in &h.bins {
            assert_eq!(c, 1, "each bin should have exactly 1 element");
        }
    }

    #[test]
    fn histogram_clamped_out_of_range() {
        let data = vec![-5.0, 5.0, 15.0]; // -5 below lo, 15 above hi
        let h = Histogram::compute(&data, 0.0, 10.0, 2);
        assert_eq!(
            h.total(),
            3,
            "out-of-range values should be clamped into boundary bins"
        );
    }

    // ── WelfordStats extended ────────────────────────────────────────────

    #[test]
    fn welford_sample_variance_two_samples() {
        let mut w = WelfordStats::default();
        w.update(2.0);
        w.update(4.0);
        // sample variance = (mean_sq_diff) / (n-1) = 2.0 / 1 = 2.0
        let sv = w.sample_variance();
        assert!((sv - 2.0).abs() < 1e-10, "sample_var = {sv}");
    }

    #[test]
    fn welford_std_dev_known_dataset() {
        let mut w = WelfordStats::default();
        for &v in &[2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0] {
            w.update(v);
        }
        assert!(
            (w.std_dev() - 2.0).abs() < 1e-10,
            "std_dev = {}",
            w.std_dev()
        );
    }

    // ── Prefix sum u64 (additional) ───────────────────────────────────────

    #[test]
    fn exclusive_scan_u64_empty() {
        let r = exclusive_scan_u64(&[]);
        assert!(r.is_empty());
    }

    #[test]
    fn inclusive_scan_u64_single() {
        let r = inclusive_scan_u64(&[7u64]);
        assert_eq!(r, vec![7]);
    }

    // ── Tile operations ───────────────────────────────────────────────────

    #[test]
    fn tile_reduce_max_and_min() {
        let t = Tile::from_slice(&[3.0, 1.0, 4.0, 1.0, 5.0]);
        assert!((t.reduce_max() - 5.0).abs() < 1e-12);
        assert!((t.reduce_min() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn tiled_reducer_tile_sums_length() {
        let data: Vec<f64> = (0..100).map(|i| i as f64).collect();
        let r = TiledReducer::new(16);
        let ts = r.tile_sums(&data);
        assert_eq!(ts.len(), 7); // ceil(100/16)
    }

    #[test]
    fn tiled_reducer_max_and_min() {
        let data = vec![-5.0, 3.0, 8.0, -1.0, 2.0];
        let r = TiledReducer::new(4);
        assert!((r.max(&data) - 8.0).abs() < 1e-12);
        assert!((r.min(&data) - (-5.0)).abs() < 1e-12);
    }
}
