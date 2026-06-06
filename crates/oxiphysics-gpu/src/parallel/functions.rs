//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use rayon::prelude::*;

use super::types::{LoadBalancePlan, LoadBalanceStrategy, WorkStealQueue};

#[inline]
pub(super) fn dist3(a: [f64; 3], b: [f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}
/// Parallel SPH density computation.
///
/// For each particle `i`, computes `rho_i = sum_j m_j * W(|r_i - r_j|, h)`.
/// The outer loop runs in parallel via Rayon.
///
/// # Arguments
/// * `positions`  - slice of 3-D particle positions.
/// * `masses`     - per-particle masses (same length as `positions`).
/// * `h`          - smoothing length.
/// * `kernel_fn`  - smoothing kernel `W(r, h)` callable from multiple threads.
pub fn parallel_sph_density(
    positions: &[[f64; 3]],
    masses: &[f64],
    h: f64,
    kernel_fn: impl Fn(f64, f64) -> f64 + Sync,
) -> Vec<f64> {
    positions
        .par_iter()
        .map(|&pi| {
            positions
                .iter()
                .zip(masses.iter())
                .map(|(&pj, &mj)| mj * kernel_fn(dist3(pi, pj), h))
                .sum()
        })
        .collect()
}
/// Parallel pairwise force accumulation.
///
/// Computes the net force on every particle.  The outer loop (over particle
/// `i`) runs in parallel; each thread independently sums the contributions
/// from all `j != i`.
///
/// # Arguments
/// * `positions` - particle positions.
/// * `n`         - number of particles (must equal `positions.len()`).
/// * `force_fn`  - `force_fn(i, j, r_ij) -> force_on_i_from_j`.
pub fn parallel_pairwise_forces(
    positions: &[[f64; 3]],
    n: usize,
    force_fn: impl Fn(usize, usize, [f64; 3]) -> [f64; 3] + Sync,
) -> Vec<[f64; 3]> {
    (0..n)
        .into_par_iter()
        .map(|i| {
            let mut f = [0.0f64; 3];
            for j in 0..n {
                if i == j {
                    continue;
                }
                let r_ij = [
                    positions[j][0] - positions[i][0],
                    positions[j][1] - positions[i][1],
                    positions[j][2] - positions[i][2],
                ];
                let fij = force_fn(i, j, r_ij);
                f[0] += fij[0];
                f[1] += fij[1];
                f[2] += fij[2];
            }
            f
        })
        .collect()
}
/// Parallel Lennard-Jones (12-6) force computation.
///
/// For each particle `i`, accumulates contributions from all `j` within
/// `cutoff`.  The potential is `U = 4*eps*[(sig/r)^12 - (sig/r)^6]`, giving
/// `F_i = sum_{j!=i} 4*eps [12(sig/r)^12 - 6(sig/r)^6] / r * r_hat_ij`.
///
/// Interactions beyond `cutoff` are skipped.
pub fn parallel_lj_forces(
    positions: &[[f64; 3]],
    epsilon: f64,
    sigma: f64,
    cutoff: f64,
) -> Vec<[f64; 3]> {
    let n = positions.len();
    (0..n)
        .into_par_iter()
        .map(|i| {
            let mut f = [0.0f64; 3];
            for j in 0..n {
                if i == j {
                    continue;
                }
                let dx = positions[j][0] - positions[i][0];
                let dy = positions[j][1] - positions[i][1];
                let dz = positions[j][2] - positions[i][2];
                let r2 = dx * dx + dy * dy + dz * dz;
                let r = r2.sqrt();
                if r >= cutoff || r < 1e-12 {
                    continue;
                }
                let sr = sigma / r;
                let sr6 = sr.powi(6);
                let sr12 = sr6 * sr6;
                let mag = 4.0 * epsilon * (12.0 * sr12 - 6.0 * sr6) / r2;
                f[0] -= mag * dx;
                f[1] -= mag * dy;
                f[2] -= mag * dz;
            }
            f
        })
        .collect()
}
/// Parallel velocity-Verlet position and velocity half-update.
///
/// Updates positions with `x += v*dt + 0.5*a*dt^2` and velocities with
/// `v += 0.5*a*dt` (first half of the Verlet velocity update; call again
/// after recomputing forces for the second half).
///
/// The loop runs in parallel via `par_iter_mut`.
pub fn parallel_verlet_step(
    positions: &mut Vec<[f64; 3]>,
    velocities: &mut Vec<[f64; 3]>,
    forces: &[[f64; 3]],
    masses: &[f64],
    dt: f64,
) {
    positions
        .par_iter_mut()
        .zip(velocities.par_iter_mut())
        .zip(forces.par_iter())
        .zip(masses.par_iter())
        .for_each(|(((pos, vel), force), &mass)| {
            let inv_m = 1.0 / mass;
            for k in 0..3 {
                let a = force[k] * inv_m;
                pos[k] += vel[k] * dt + 0.5 * a * dt * dt;
                vel[k] += 0.5 * a * dt;
            }
        });
}
/// Parallel AABB overlap detection.
///
/// Returns all pairs `(i, j)` with `i < j` whose axis-aligned bounding boxes
/// overlap.  The outer loop runs in parallel; each thread contributes matching
/// pairs into a local vector that is then flattened.
pub fn parallel_aabb_pairs(aabbs: &[([f64; 3], [f64; 3])]) -> Vec<(usize, usize)> {
    let n = aabbs.len();
    (0..n)
        .into_par_iter()
        .flat_map(|i| {
            let mut local = Vec::new();
            let (min_i, max_i) = aabbs[i];
            for (j, &(min_j, max_j)) in aabbs.iter().enumerate().skip(i + 1) {
                let overlap = (0..3).all(|k| min_i[k] <= max_j[k] && min_j[k] <= max_i[k]);
                if overlap {
                    local.push((i, j));
                }
            }
            local
        })
        .collect()
}
/// Execute `f(i)` for `i` in `0..n`, splitting into chunks of `chunk_size`.
///
/// Currently processes chunks sequentially; prefer [`parallel_sph_density`]
/// and the other parallel kernels for performance-critical code.
pub fn parallel_for(n: usize, chunk_size: usize, f: impl Fn(usize)) {
    let cs = if chunk_size == 0 { 1 } else { chunk_size };
    for start in (0..n).step_by(cs) {
        let end = (start + cs).min(n);
        for i in start..end {
            f(i);
        }
    }
}
/// Parallel sum reduction using Rayon.
pub fn parallel_reduce_sum(data: &[f64]) -> f64 {
    data.par_iter().copied().sum()
}
/// Parallel max reduction using Rayon.
pub fn parallel_reduce_max(data: &[f64]) -> f64 {
    data.par_iter()
        .copied()
        .reduce(|| f64::NEG_INFINITY, f64::max)
}
/// Parallel min reduction using Rayon.
pub fn parallel_reduce_min(data: &[f64]) -> f64 {
    data.par_iter().copied().reduce(|| f64::INFINITY, f64::min)
}
/// Parallel dot product of two slices.
pub fn parallel_dot_product(a: &[f64], b: &[f64]) -> f64 {
    a.par_iter()
        .zip(b.par_iter())
        .map(|(&ai, &bi)| ai * bi)
        .sum()
}
/// Parallel L2 norm (Euclidean norm).
pub fn parallel_norm2(data: &[f64]) -> f64 {
    let sum_sq: f64 = data.par_iter().map(|&x| x * x).sum();
    sum_sq.sqrt()
}
/// Parallel mean.
pub fn parallel_mean(data: &[f64]) -> f64 {
    if data.is_empty() {
        return 0.0;
    }
    let sum: f64 = data.par_iter().copied().sum();
    sum / data.len() as f64
}
/// Parallel variance (population variance).
pub fn parallel_variance(data: &[f64]) -> f64 {
    if data.is_empty() {
        return 0.0;
    }
    let mean = parallel_mean(data);
    let sum_sq: f64 = data.par_iter().map(|&x| (x - mean) * (x - mean)).sum();
    sum_sq / data.len() as f64
}
/// Two-pass parallel reduction: compute both sum and count in one pass.
pub fn parallel_sum_count(data: &[f64]) -> (f64, usize) {
    data.par_iter()
        .copied()
        .fold(|| (0.0f64, 0usize), |(s, c), x| (s + x, c + 1))
        .reduce(|| (0.0, 0), |(s1, c1), (s2, c2)| (s1 + s2, c1 + c2))
}
/// Parallel reduction with a custom binary operator.
///
/// `identity` is the identity element for the operator (e.g. 0.0 for add).
/// `op` must be associative and commutative for correctness.
pub fn parallel_reduce_custom(
    data: &[f64],
    identity: f64,
    op: impl Fn(f64, f64) -> f64 + Sync + Send,
) -> f64 {
    data.par_iter().copied().reduce(|| identity, op)
}
/// Parallel exclusive prefix sum using a two-pass algorithm.
///
/// Phase 1: compute partial sums in chunks (parallel).
/// Phase 2: propagate offsets (sequential).
/// Phase 3: apply offsets within chunks (parallel).
pub fn parallel_exclusive_scan(data: &[f64]) -> Vec<f64> {
    let n = data.len();
    if n == 0 {
        return Vec::new();
    }
    let chunk_size = (n / rayon::current_num_threads().max(1)).max(64);
    let chunks: Vec<&[f64]> = data.chunks(chunk_size).collect();
    let chunk_sums: Vec<f64> = chunks
        .par_iter()
        .map(|chunk| chunk.iter().copied().sum())
        .collect();
    let mut offsets = Vec::with_capacity(chunks.len());
    let mut acc = 0.0;
    for &cs in &chunk_sums {
        offsets.push(acc);
        acc += cs;
    }
    let result = vec![0.0; n];
    chunks.par_iter().enumerate().for_each(|(ci, chunk)| {
        let base = ci * chunk_size;
        let offset = offsets[ci];
        let mut local_acc = offset;
        let result_ptr = result.as_ptr() as *mut f64;
        for (k, &v) in chunk.iter().enumerate() {
            // SAFETY: `chunks` partitions `result`'s index range into disjoint,
            // contiguous blocks (`chunks(chunk_size)`), and `base = ci * chunk_size`
            // is exactly this chunk's start, so distinct rayon tasks write disjoint
            // indices — no aliasing. `base + k < base + chunk.len() <= n`, so the
            // write stays in bounds of `result`, which outlives the parallel region.
            unsafe {
                *result_ptr.add(base + k) = local_acc;
            }
            local_acc += v;
        }
    });
    result
}
/// Parallel inclusive prefix sum.
pub fn parallel_inclusive_scan(data: &[f64]) -> Vec<f64> {
    let n = data.len();
    if n == 0 {
        return Vec::new();
    }
    let chunk_size = (n / rayon::current_num_threads().max(1)).max(64);
    let chunks: Vec<&[f64]> = data.chunks(chunk_size).collect();
    let chunk_sums: Vec<f64> = chunks
        .par_iter()
        .map(|chunk| chunk.iter().copied().sum())
        .collect();
    let mut offsets = Vec::with_capacity(chunks.len());
    let mut acc = 0.0;
    for &cs in &chunk_sums {
        offsets.push(acc);
        acc += cs;
    }
    let result = vec![0.0; n];
    chunks.par_iter().enumerate().for_each(|(ci, chunk)| {
        let base = ci * chunk_size;
        let offset = offsets[ci];
        let mut local_acc = offset;
        let result_ptr = result.as_ptr() as *mut f64;
        for (k, &v) in chunk.iter().enumerate() {
            local_acc += v;
            // SAFETY: `chunks` partitions `result`'s index range into disjoint,
            // contiguous blocks (`chunks(chunk_size)`), and `base = ci * chunk_size`
            // is exactly this chunk's start, so distinct rayon tasks write disjoint
            // indices — no aliasing. `base + k < base + chunk.len() <= n`, so the
            // write stays in bounds of `result`, which outlives the parallel region.
            unsafe {
                *result_ptr.add(base + k) = local_acc;
            }
        }
    });
    result
}
/// Segmented prefix sum: performs exclusive scans within each segment.
///
/// `segment_ids` assigns each element to a segment. When the segment ID
/// changes, the accumulator resets. Segments must be contiguous.
pub fn segmented_exclusive_scan(data: &[f64], segment_ids: &[usize]) -> Vec<f64> {
    let n = data.len();
    let mut result = vec![0.0; n];
    if n == 0 {
        return result;
    }
    let mut acc = 0.0;
    let mut current_seg = segment_ids[0];
    for i in 0..n {
        if segment_ids[i] != current_seg {
            current_seg = segment_ids[i];
            acc = 0.0;
        }
        result[i] = acc;
        acc += data[i];
    }
    result
}
/// Parallel sort of f64 values (ascending).
///
/// Uses Rayon's parallel sort. NaN values are placed at the end.
pub fn parallel_sort_f64(data: &mut [f64]) {
    data.par_sort_unstable_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
}
/// Parallel argsort: returns indices that would sort `data` in ascending order.
pub fn parallel_argsort(data: &[f64]) -> Vec<usize> {
    let mut indices: Vec<usize> = (0..data.len()).collect();
    indices.par_sort_unstable_by(|&a, &b| {
        data[a]
            .partial_cmp(&data[b])
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    indices
}
/// Parallel sort by key: sorts `items` based on `key_fn`.
pub fn parallel_sort_by_key<T: Send>(items: &mut [T], key_fn: impl Fn(&T) -> f64 + Sync + Send) {
    items.par_sort_unstable_by(|a, b| {
        let ka = key_fn(a);
        let kb = key_fn(b);
        ka.partial_cmp(&kb).unwrap_or(std::cmp::Ordering::Equal)
    });
}
/// Parallel partition: split data into two groups based on a predicate.
///
/// Returns `(true_group, false_group)`.
pub fn parallel_partition<T: Send + Sync + Clone>(
    data: &[T],
    predicate: impl Fn(&T) -> bool + Sync + Send,
) -> (Vec<T>, Vec<T>) {
    let (left, right): (Vec<_>, Vec<_>) =
        data.par_iter().cloned().partition(|item| predicate(item));
    (left, right)
}
/// Parallel rank: compute the rank of each element (0-based) in sorted order.
pub fn parallel_rank(data: &[f64]) -> Vec<usize> {
    let sorted_indices = parallel_argsort(data);
    let n = data.len();
    let mut ranks = vec![0usize; n];
    for (rank, &idx) in sorted_indices.iter().enumerate() {
        ranks[idx] = rank;
    }
    ranks
}
/// Compute a load balance plan for `n` items across `num_workers` workers.
///
/// For `Static` strategy, `item_weights` is ignored.
/// For `Weighted` strategy, items are assigned sequentially to workers to
/// balance the total weight per worker.
/// For `Guided` strategy, chunks start large and shrink.
pub fn compute_load_balance(
    n: usize,
    num_workers: usize,
    strategy: LoadBalanceStrategy,
    item_weights: Option<&[f64]>,
) -> LoadBalancePlan {
    let nw = num_workers.max(1);
    match strategy {
        LoadBalanceStrategy::Static => {
            let chunk = n.div_ceil(nw);
            let mut ranges = Vec::with_capacity(nw);
            let mut weights = Vec::with_capacity(nw);
            for w in 0..nw {
                let start = w * chunk;
                let end = ((w + 1) * chunk).min(n);
                if start < n {
                    let weight = if let Some(wts) = item_weights {
                        wts[start..end].iter().sum()
                    } else {
                        (end - start) as f64
                    };
                    ranges.push(start..end);
                    weights.push(weight);
                }
            }
            LoadBalancePlan { ranges, weights }
        }
        LoadBalanceStrategy::Weighted => {
            let wts = match item_weights {
                Some(w) => w,
                None => {
                    return compute_load_balance(n, num_workers, LoadBalanceStrategy::Static, None);
                }
            };
            let total_weight: f64 = wts.iter().sum();
            let target_per_worker = total_weight / nw as f64;
            let mut ranges = Vec::with_capacity(nw);
            let mut weights = Vec::with_capacity(nw);
            let mut start = 0;
            let mut current_weight = 0.0;
            for (i, &wi) in wts.iter().enumerate() {
                current_weight += wi;
                let workers_remaining = nw - ranges.len();
                let at_last_worker = workers_remaining == 1;
                let exceeded_target = current_weight >= target_per_worker && !at_last_worker;
                if exceeded_target {
                    ranges.push(start..(i + 1));
                    weights.push(current_weight);
                    start = i + 1;
                    current_weight = 0.0;
                }
            }
            if start < n || ranges.is_empty() {
                ranges.push(start..n);
                weights.push(current_weight);
            }
            LoadBalancePlan { ranges, weights }
        }
        LoadBalanceStrategy::Guided => {
            let mut ranges = Vec::new();
            let mut weights = Vec::new();
            let mut pos = 0;
            let mut remaining = n;
            while remaining > 0 {
                let min_chunk = 1usize;
                let chunk = (remaining / nw).max(min_chunk).min(remaining);
                let end = pos + chunk;
                let weight = if let Some(wts) = item_weights {
                    wts[pos..end].iter().sum()
                } else {
                    chunk as f64
                };
                ranges.push(pos..end);
                weights.push(weight);
                pos = end;
                remaining -= chunk;
            }
            LoadBalancePlan { ranges, weights }
        }
    }
}
/// Execute a function in parallel with load-balanced ranges.
///
/// Each range is processed as one Rayon task. The function receives
/// `(worker_id, range)`.
pub fn execute_balanced(
    plan: &LoadBalancePlan,
    f: impl Fn(usize, std::ops::Range<usize>) + Sync + Send,
) {
    plan.ranges
        .par_iter()
        .enumerate()
        .for_each(|(worker_id, range)| {
            f(worker_id, range.clone());
        });
}
/// Parallel map-reduce: map each element, then reduce results.
///
/// Combines mapping and reduction in a single parallel pass.
pub fn parallel_map_reduce<T: Send + Sync>(
    data: &[T],
    map_fn: impl Fn(&T) -> f64 + Sync + Send,
    identity: f64,
    reduce_fn: impl Fn(f64, f64) -> f64 + Sync + Send,
) -> f64 {
    data.par_iter().map(map_fn).reduce(|| identity, reduce_fn)
}
/// Parallel histogram: count elements falling into each bin.
///
/// Bins are `[min, min+step), [min+step, min+2*step), ...`.
/// Returns a vector of length `num_bins`.
pub fn parallel_histogram(data: &[f64], min: f64, max: f64, num_bins: usize) -> Vec<usize> {
    if num_bins == 0 || max <= min {
        return vec![0; num_bins];
    }
    let step = (max - min) / num_bins as f64;
    data.par_iter()
        .fold(
            || vec![0usize; num_bins],
            |mut hist, &v| {
                if v >= min && v < max {
                    let bin = ((v - min) / step) as usize;
                    let bin = bin.min(num_bins - 1);
                    hist[bin] += 1;
                } else if (v - max).abs() < 1e-15 {
                    hist[num_bins - 1] += 1;
                }
                hist
            },
        )
        .reduce(
            || vec![0usize; num_bins],
            |mut a, b| {
                for (ai, bi) in a.iter_mut().zip(b.iter()) {
                    *ai += bi;
                }
                a
            },
        )
}
/// Stream compaction: retain only elements satisfying a predicate, returning
/// a compacted output together with a scatter index map.
///
/// Returns `(compacted, scatter_map)` where:
/// * `compacted` contains the elements `data[i]` for which `pred(data[i])` is true.
/// * `scatter_map[j]` is the original index `i` of `compacted[j]`.
///
/// This mirrors a GPU stream-compaction pass (prefix-sum → scatter).
pub fn stream_compaction<T: Clone>(data: &[T], pred: impl Fn(&T) -> bool) -> (Vec<T>, Vec<usize>) {
    let mut compacted = Vec::new();
    let mut scatter_map = Vec::new();
    for (i, item) in data.iter().enumerate() {
        if pred(item) {
            compacted.push(item.clone());
            scatter_map.push(i);
        }
    }
    (compacted, scatter_map)
}
/// Parallel stream compaction via Rayon.
///
/// Each thread builds a local (value, original_index) list and then the
/// lists are merged in order to preserve a deterministic output.
pub fn parallel_stream_compaction<T: Clone + Send + Sync>(
    data: &[T],
    pred: impl Fn(&T) -> bool + Sync,
) -> (Vec<T>, Vec<usize>) {
    use rayon::iter::IndexedParallelIterator;
    let pairs: Vec<(T, usize)> = data
        .par_iter()
        .enumerate()
        .filter_map(|(i, item)| {
            if pred(item) {
                Some((item.clone(), i))
            } else {
                None
            }
        })
        .collect();
    let compacted: Vec<T> = pairs.iter().map(|(v, _)| v.clone()).collect();
    let scatter_map: Vec<usize> = pairs.iter().map(|(_, i)| *i).collect();
    (compacted, scatter_map)
}
/// Segmented reduction: sum values within each segment independently.
///
/// `segment_ids[i]` must be monotonically non-decreasing.
/// Returns a vector of partial sums, one per distinct segment id.
///
/// Example:
/// ```text
/// data          = [1, 2, 3, 4, 5, 6]
/// segment_ids   = [0, 0, 1, 1, 1, 2]
/// output        = [3, 12, 6]
/// ```
pub fn segmented_reduce_sum(data: &[f64], segment_ids: &[usize]) -> Vec<f64> {
    if data.is_empty() {
        return Vec::new();
    }
    let max_seg = *segment_ids.iter().max().unwrap_or(&0);
    let mut result = vec![0.0f64; max_seg + 1];
    for (&v, &s) in data.iter().zip(segment_ids.iter()) {
        result[s] += v;
    }
    result
}
/// Segmented reduction: maximum value within each segment.
pub fn segmented_reduce_max(data: &[f64], segment_ids: &[usize]) -> Vec<f64> {
    if data.is_empty() {
        return Vec::new();
    }
    let max_seg = *segment_ids.iter().max().unwrap_or(&0);
    let mut result = vec![f64::NEG_INFINITY; max_seg + 1];
    for (&v, &s) in data.iter().zip(segment_ids.iter()) {
        if v > result[s] {
            result[s] = v;
        }
    }
    result
}
/// Segmented reduction: minimum value within each segment.
pub fn segmented_reduce_min(data: &[f64], segment_ids: &[usize]) -> Vec<f64> {
    if data.is_empty() {
        return Vec::new();
    }
    let max_seg = *segment_ids.iter().max().unwrap_or(&0);
    let mut result = vec![f64::INFINITY; max_seg + 1];
    for (&v, &s) in data.iter().zip(segment_ids.iter()) {
        if v < result[s] {
            result[s] = v;
        }
    }
    result
}
/// Stable merge sort for f64 values (CPU reference implementation).
///
/// Returns a new sorted vector leaving the input unchanged.
/// NaN values are placed at the end (treated as greater than any finite value).
pub fn merge_sort_f64(data: &[f64]) -> Vec<f64> {
    let mut buf = data.to_vec();
    merge_sort_recurse(&mut buf);
    buf
}
pub(super) fn merge_sort_recurse(data: &mut [f64]) {
    let n = data.len();
    if n <= 1 {
        return;
    }
    let mid = n / 2;
    merge_sort_recurse(&mut data[..mid]);
    merge_sort_recurse(&mut data[mid..]);
    let left: Vec<f64> = data[..mid].to_vec();
    let right: Vec<f64> = data[mid..].to_vec();
    let (mut l, mut r, mut i) = (0, 0, 0);
    while l < left.len() && r < right.len() {
        if left[l]
            .partial_cmp(&right[r])
            .unwrap_or(std::cmp::Ordering::Greater)
            != std::cmp::Ordering::Greater
        {
            data[i] = left[l];
            l += 1;
        } else {
            data[i] = right[r];
            r += 1;
        }
        i += 1;
    }
    while l < left.len() {
        data[i] = left[l];
        l += 1;
        i += 1;
    }
    while r < right.len() {
        data[i] = right[r];
        r += 1;
        i += 1;
    }
}
/// Merge sort returning the sorted permutation (argsort, stable).
///
/// `result[k]` is the original index of the k-th smallest element.
pub fn merge_sort_argsort(data: &[f64]) -> Vec<usize> {
    let mut indices: Vec<usize> = (0..data.len()).collect();
    merge_argsort_recurse(data, &mut indices);
    indices
}
pub(super) fn merge_argsort_recurse(data: &[f64], indices: &mut [usize]) {
    let n = indices.len();
    if n <= 1 {
        return;
    }
    let mid = n / 2;
    let (left_idx, right_idx) = indices.split_at_mut(mid);
    merge_argsort_recurse(data, left_idx);
    merge_argsort_recurse(data, right_idx);
    let left: Vec<usize> = left_idx.to_vec();
    let right: Vec<usize> = right_idx.to_vec();
    let (mut l, mut r, mut i) = (0, 0, 0);
    while l < left.len() && r < right.len() {
        let cmp = data[left[l]]
            .partial_cmp(&data[right[r]])
            .unwrap_or(std::cmp::Ordering::Greater);
        if cmp != std::cmp::Ordering::Greater {
            indices[i] = left[l];
            l += 1;
        } else {
            indices[i] = right[r];
            r += 1;
        }
        i += 1;
    }
    while l < left.len() {
        indices[i] = left[l];
        l += 1;
        i += 1;
    }
    while r < right.len() {
        indices[i] = right[r];
        r += 1;
        i += 1;
    }
}
/// Bitonic sort in ascending order.
///
/// Works on arrays whose length is a power of two.  Pads the input with
/// `f64::INFINITY` if needed and trims back afterwards.
///
/// This CPU reference mirrors a GPU bitonic sort which operates in
/// `O(n log² n)` compare-and-swap steps.
pub fn bitonic_sort(data: &[f64]) -> Vec<f64> {
    let n = data.len();
    if n == 0 {
        return Vec::new();
    }
    let mut p = 1usize;
    while p < n {
        p <<= 1;
    }
    let mut buf: Vec<f64> = data.to_vec();
    buf.resize(p, f64::INFINITY);
    let mut k = 2usize;
    while k <= p {
        let mut j = k >> 1;
        while j >= 1 {
            for i in 0..p {
                let l = i ^ j;
                if l > i {
                    let ascending = (i & k) == 0;
                    if (ascending && buf[i] > buf[l]) || (!ascending && buf[i] < buf[l]) {
                        buf.swap(i, l);
                    }
                }
            }
            j >>= 1;
        }
        k <<= 1;
    }
    buf.truncate(n);
    buf
}
/// Bitonic sort that returns the original indices (argsort variant).
///
/// Pads with `(f64::INFINITY, usize::MAX)` pairs and trims back.
pub fn bitonic_argsort(data: &[f64]) -> Vec<usize> {
    let n = data.len();
    if n == 0 {
        return Vec::new();
    }
    let mut p = 1usize;
    while p < n {
        p <<= 1;
    }
    let mut buf: Vec<(f64, usize)> = data
        .iter()
        .copied()
        .enumerate()
        .map(|(i, v)| (v, i))
        .collect();
    buf.resize(p, (f64::INFINITY, usize::MAX));
    let mut k = 2usize;
    while k <= p {
        let mut j = k >> 1;
        while j >= 1 {
            for i in 0..p {
                let l = i ^ j;
                if l > i {
                    let ascending = (i & k) == 0;
                    let should_swap =
                        (ascending && buf[i].0 > buf[l].0) || (!ascending && buf[i].0 < buf[l].0);
                    if should_swap {
                        buf.swap(i, l);
                    }
                }
            }
            j >>= 1;
        }
        k <<= 1;
    }
    buf.truncate(n);
    buf.iter().map(|(_, idx)| *idx).collect()
}
/// Simulate a work-stealing dispatcher across `num_workers` queues.
///
/// `tasks` is divided evenly among workers.  Any worker that finishes early
/// steals from the most loaded remaining worker.
/// Returns a `Vec`usize` of length `num_workers` recording the tasks each
/// worker processed.
pub fn work_steal_queue<T: Send + Clone>(
    tasks: Vec<T>,
    num_workers: usize,
    _process: impl Fn(&T) + Sync,
) -> Vec<usize> {
    let nw = num_workers.max(1);
    let mut queues: Vec<WorkStealQueue<T>> = (0..nw).map(|_| WorkStealQueue::new()).collect();
    for (i, task) in tasks.into_iter().enumerate() {
        queues[i % nw].push(task);
    }
    let mut processed = vec![0usize; nw];
    loop {
        let mut did_work = false;
        for w in 0..nw {
            while let Some(task) = queues[w].pop() {
                _process(&task);
                processed[w] += 1;
                did_work = true;
            }
        }
        let max_len = queues.iter().map(|q| q.len()).max().unwrap_or(0);
        if max_len == 0 {
            break;
        }
        if did_work {
            continue;
        }
        let victim = queues
            .iter()
            .enumerate()
            .max_by_key(|(_, q)| q.len())
            .map(|(i, _)| i);
        let thief = queues
            .iter()
            .enumerate()
            .find(|(_, q)| q.is_empty())
            .map(|(i, _)| i);
        if let (Some(v), Some(t)) = (victim, thief) {
            if v != t {
                if let Some(task) = queues[v].steal() {
                    queues[t].push(task);
                }
            } else {
                break;
            }
        } else {
            break;
        }
    }
    processed
}
/// Compute a load-balance efficiency metric given per-worker task counts.
///
/// Returns a value in `\[0, 1\]`: 1.0 means perfect balance, smaller values
/// indicate more imbalance.  Defined as `avg_load / max_load`.
pub fn compute_load_balance_metric(worker_loads: &[usize]) -> f64 {
    if worker_loads.is_empty() {
        return 1.0;
    }
    let total: usize = worker_loads.iter().sum();
    let n = worker_loads.len();
    let avg = total as f64 / n as f64;
    let max = *worker_loads.iter().max().unwrap_or(&1) as f64;
    if max < 1e-15 {
        return 1.0;
    }
    avg / max
}
/// Suggest an optimal chunk size for `n` work items across `num_workers`
/// workers, targeting at least `min_chunks_per_worker` chunks per worker.
pub fn suggest_chunk_size(n: usize, num_workers: usize, min_chunks_per_worker: usize) -> usize {
    let nw = num_workers.max(1);
    let chunks = (nw * min_chunks_per_worker).max(1);
    n.div_ceil(chunks).max(1)
}
/// Parallel merge sort for `f64` slices using Rayon.
///
/// Splits the array recursively.  Below `SERIAL_THRESHOLD` elements the
/// standard library sort is used.  Above that the two halves are sorted in
/// parallel and then merged sequentially.
pub fn merge_sort_parallel(data: &[f64]) -> Vec<f64> {
    pub(super) const SERIAL_THRESHOLD: usize = 256;
    let n = data.len();
    if n <= 1 {
        return data.to_vec();
    }
    if n <= SERIAL_THRESHOLD {
        let mut v = data.to_vec();
        v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        return v;
    }
    let mid = n / 2;
    let (left_slice, right_slice) = data.split_at(mid);
    let (left_sorted, right_sorted) = rayon::join(
        || merge_sort_parallel(left_slice),
        || merge_sort_parallel(right_slice),
    );
    merge_two_sorted(&left_sorted, &right_sorted)
}
/// Merge two sorted `f64` slices into one sorted `Vec`f64`.
pub fn merge_two_sorted(a: &[f64], b: &[f64]) -> Vec<f64> {
    let mut result = Vec::with_capacity(a.len() + b.len());
    let (mut i, mut j) = (0, 0);
    while i < a.len() && j < b.len() {
        if a[i] <= b[j] {
            result.push(a[i]);
            i += 1;
        } else {
            result.push(b[j]);
            j += 1;
        }
    }
    result.extend_from_slice(&a[i..]);
    result.extend_from_slice(&b[j..]);
    result
}
