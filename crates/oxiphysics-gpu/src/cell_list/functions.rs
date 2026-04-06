//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
use rayon::prelude::*;

use super::types::GpuCellList;

/// Expand a 10-bit integer by inserting 2 zero bits after each bit.
///
/// `x` must be in `[0, 1023]`.
pub(super) fn expand_bits(mut x: u32) -> u32 {
    x &= 0x000003FF;
    x = (x | (x << 16)) & 0x030000FF;
    x = (x | (x << 8)) & 0x0300F00F;
    x = (x | (x << 4)) & 0x030C30C3;
    x = (x | (x << 2)) & 0x09249249;
    x
}
/// Compute a 30-bit Morton code from 3D integer coordinates.
///
/// Each coordinate must be in `[0, 1023]`.
pub fn morton_encode(x: u32, y: u32, z: u32) -> u32 {
    expand_bits(x) | (expand_bits(y) << 1) | (expand_bits(z) << 2)
}
/// Decode a 30-bit Morton code back to 3D integer coordinates.
pub fn morton_decode(code: u32) -> (u32, u32, u32) {
    (
        compact_bits(code),
        compact_bits(code >> 1),
        compact_bits(code >> 2),
    )
}
/// Compact every third bit (inverse of expand_bits).
pub(super) fn compact_bits(mut x: u32) -> u32 {
    x &= 0x09249249;
    x = (x | (x >> 2)) & 0x030C30C3;
    x = (x | (x >> 4)) & 0x0300F00F;
    x = (x | (x >> 8)) & 0x030000FF;
    x = (x | (x >> 16)) & 0x000003FF;
    x
}
/// Compute Morton codes for a set of positions and sort indices by Morton code.
///
/// `positions` - particle positions
/// `box_min` - minimum corner of the bounding box
/// `box_max` - maximum corner of the bounding box
///
/// Returns `(sorted_indices, morton_codes)` where `sorted_indices[i]` is the
/// original index of the particle that should be in position `i` after sorting.
pub fn morton_sort(
    positions: &[[f64; 3]],
    box_min: [f64; 3],
    box_max: [f64; 3],
) -> (Vec<usize>, Vec<u32>) {
    let range = [
        (box_max[0] - box_min[0]).max(1e-10),
        (box_max[1] - box_min[1]).max(1e-10),
        (box_max[2] - box_min[2]).max(1e-10),
    ];
    let mut codes: Vec<(u32, usize)> = positions
        .par_iter()
        .enumerate()
        .map(|(i, p)| {
            let x = (((p[0] - box_min[0]) / range[0] * 1023.0) as u32).min(1023);
            let y = (((p[1] - box_min[1]) / range[1] * 1023.0) as u32).min(1023);
            let z = (((p[2] - box_min[2]) / range[2] * 1023.0) as u32).min(1023);
            (morton_encode(x, y, z), i)
        })
        .collect();
    codes.sort_by_key(|&(code, _)| code);
    let sorted_indices: Vec<usize> = codes.iter().map(|&(_, idx)| idx).collect();
    let morton_codes: Vec<u32> = codes.iter().map(|&(code, _)| code).collect();
    (sorted_indices, morton_codes)
}
/// Compute an exclusive prefix sum (scan) over `counts`.
pub fn parallel_prefix_sum(counts: &[usize]) -> Vec<usize> {
    let mut out = Vec::with_capacity(counts.len());
    let mut acc = 0usize;
    for &c in counts {
        out.push(acc);
        acc += c;
    }
    out
}
/// Compute the bounding box of a set of positions.
///
/// Returns `(min, max)` corner coordinates.
pub fn compute_bounding_box(positions: &[[f64; 3]]) -> ([f64; 3], [f64; 3]) {
    if positions.is_empty() {
        return ([0.0; 3], [0.0; 3]);
    }
    let mut min = positions[0];
    let mut max = positions[0];
    for p in positions {
        for d in 0..3 {
            if p[d] < min[d] {
                min[d] = p[d];
            }
            if p[d] > max[d] {
                max[d] = p[d];
            }
        }
    }
    (min, max)
}
/// Reorder an array according to a permutation.
///
/// `perm[i]` is the original index of the element that should appear at position `i`.
pub fn reorder_by_permutation<T: Clone>(data: &[T], perm: &[usize]) -> Vec<T> {
    perm.iter().map(|&i| data[i].clone()).collect()
}
/// Mock radix sort: stable sort of `keys` returning `(sorted_keys, sorted_indices)`.
///
/// `sorted_indices[i]` is the original index of the element now at position `i`.
/// This is a CPU reference implementation matching the expected GPU radix sort output.
pub fn radix_sort_mock(keys: &[u32]) -> (Vec<u32>, Vec<usize>) {
    if keys.is_empty() {
        return (vec![], vec![]);
    }
    let mut indexed: Vec<(u32, usize)> = keys.iter().copied().zip(0..).collect();
    indexed.sort_by_key(|&(k, _)| k);
    let sorted_keys: Vec<u32> = indexed.iter().map(|&(k, _)| k).collect();
    let sorted_indices: Vec<usize> = indexed.iter().map(|&(_, i)| i).collect();
    (sorted_keys, sorted_indices)
}
/// Exclusive prefix sum (scan) of `counts`.
///
/// For input `[a, b, c, d]` the output is `[0, a, a+b, a+b+c]`.
/// This mirrors what a parallel GPU prefix-sum kernel would produce.
pub fn gpu_prefix_sum(counts: &[usize]) -> Vec<usize> {
    let mut out = Vec::with_capacity(counts.len());
    let mut running = 0usize;
    for &c in counts {
        out.push(running);
        running += c;
    }
    out
}
/// Count how many particles fall into each cell of a regular grid.
///
/// The grid has `n_cells = [nx, ny, nz]` cells, each of side length `cell_size`.
/// The origin is `[0, 0, 0]`.  Particles outside the grid are clamped to the
/// nearest boundary cell.
///
/// Returns a flat Vec of length `nx * ny * nz`.
pub fn parallel_count_particles(
    positions: &[[f64; 3]],
    n_cells: [usize; 3],
    cell_size: f64,
) -> Vec<usize> {
    let [nx, ny, nz] = n_cells;
    let total = nx * ny * nz;
    let mut counts = vec![0usize; total];
    for p in positions {
        let ix = ((p[0] / cell_size) as isize).clamp(0, nx as isize - 1) as usize;
        let iy = ((p[1] / cell_size) as isize).clamp(0, ny as isize - 1) as usize;
        let iz = ((p[2] / cell_size) as isize).clamp(0, nz as isize - 1) as usize;
        counts[ix + nx * (iy + ny * iz)] += 1;
    }
    counts
}
/// Distribute `n_cells` cells evenly across `n_gpus` GPUs.
///
/// Returns a Vec of non-overlapping contiguous ranges that together cover
/// `0..n_cells`.  Remainder cells are spread among the first GPUs.
pub fn distribute_cells_to_gpus(n_cells: usize, n_gpus: usize) -> Vec<std::ops::Range<usize>> {
    if n_gpus == 0 || n_cells == 0 {
        return vec![];
    }
    let base = n_cells / n_gpus;
    let remainder = n_cells % n_gpus;
    let mut ranges = Vec::with_capacity(n_gpus);
    let mut start = 0;
    for gpu in 0..n_gpus {
        let extra = if gpu < remainder { 1 } else { 0 };
        let end = start + base + extra;
        ranges.push(start..end);
        start = end;
    }
    ranges
}
/// Return all particle pairs `(i, j)` with `i < j` and `dist(i, j) < cutoff`.
///
/// Uses the [`GpuCellList`] for candidate acceleration, filtering by exact
/// distance afterward.  The returned pairs are in undefined order.
pub fn gpu_neighbor_search_kernel(
    cl: &GpuCellList,
    positions: &[[f64; 3]],
    cutoff: f64,
) -> Vec<(usize, usize)> {
    let mut pairs = Vec::new();
    cl.for_each_pair(positions, cutoff, |i, j, _d2| {
        let (a, b) = if i < j { (i, j) } else { (j, i) };
        pairs.push((a, b));
    });
    pairs.sort_unstable();
    pairs.dedup();
    pairs
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::cell_list::CellList;

    use crate::cell_list::GpuCellList;

    use crate::cell_list::SpatialHash;

    #[test]
    fn test_prefix_sum_empty() {
        assert_eq!(parallel_prefix_sum(&[]), Vec::<usize>::new());
    }
    #[test]
    fn test_prefix_sum_basic() {
        let counts = [1usize, 2, 3, 4];
        let result = parallel_prefix_sum(&counts);
        assert_eq!(result, vec![0, 1, 3, 6]);
    }
    #[test]
    fn test_cell_index_clamp() {
        let list = GpuCellList::new([4, 4, 4], 1.0, [4.0, 4.0, 4.0]);
        let idx = list.cell_index([4.5, 4.5, 4.5]);
        assert_eq!(idx, 3 * 4 * 4 + 3 * 4 + 3);
    }
    #[test]
    fn test_total_cells() {
        let list = GpuCellList::new([3, 4, 5], 1.0, [3.0, 4.0, 5.0]);
        assert_eq!(list.total_cells(), 60);
    }
    #[test]
    fn test_build_parallel_counts() {
        let positions: Vec<[f64; 3]> = vec![
            [0.5, 0.5, 0.5],
            [1.5, 0.5, 0.5],
            [0.5, 1.5, 0.5],
            [1.5, 1.5, 0.5],
            [0.5, 0.5, 1.5],
            [1.5, 0.5, 1.5],
            [0.5, 1.5, 1.5],
            [1.5, 1.5, 1.5],
        ];
        let cl = GpuCellList::build_parallel(&positions);
        assert_eq!(cl.sorted_indices.len(), 8);
        for c in 0..cl.total_cells() {
            assert_eq!(cl.cell_counts[c], 1);
        }
    }
    #[test]
    fn test_neighbors_in_radius() {
        let positions: Vec<[f64; 3]> = vec![[0.5, 0.5, 0.5], [0.6, 0.5, 0.5], [5.0, 5.0, 5.0]];
        let cl = GpuCellList::build_parallel(&positions);
        let mut neighbours = cl.neighbors_in_radius(&positions, [0.5, 0.5, 0.5], 0.5);
        neighbours.sort_unstable();
        assert!(neighbours.contains(&0));
        assert!(neighbours.contains(&1));
        assert!(!neighbours.contains(&2));
    }
    #[test]
    fn cell_list_find_neighbors_all_pairs() {
        let positions: Vec<[f64; 3]> = vec![
            [1.0, 1.0, 1.0],
            [1.2, 1.0, 1.0],
            [1.0, 1.3, 1.0],
            [9.0, 9.0, 9.0],
        ];
        let cl = CellList::build(&positions);
        let radius = 0.5;
        let mut neighbours = cl.find_neighbors([1.0, 1.0, 1.0], radius);
        neighbours.sort_unstable();
        assert!(neighbours.contains(&0), "should find self: {neighbours:?}");
        assert!(
            neighbours.contains(&1),
            "should find particle 1: {neighbours:?}"
        );
        assert!(
            neighbours.contains(&2),
            "should find particle 2: {neighbours:?}"
        );
        assert!(
            !neighbours.contains(&3),
            "particle 3 is far: {neighbours:?}"
        );
    }
    #[test]
    fn cell_list_new_compiles() {
        let cl = CellList::new([10.0, 10.0, 10.0], 2.0);
        assert_eq!(cl.inner.total_cells(), 125);
    }
    /// Morton encode/decode should be inverses.
    #[test]
    fn test_morton_roundtrip() {
        let test_cases = [
            (0, 0, 0),
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1),
            (7, 3, 5),
            (1023, 1023, 1023),
            (512, 256, 128),
        ];
        for (x, y, z) in test_cases {
            let code = morton_encode(x, y, z);
            let (dx, dy, dz) = morton_decode(code);
            assert_eq!(dx, x, "x mismatch for ({x},{y},{z})");
            assert_eq!(dy, y, "y mismatch for ({x},{y},{z})");
            assert_eq!(dz, z, "z mismatch for ({x},{y},{z})");
        }
    }
    /// Morton codes should preserve locality: nearby points have nearby codes.
    #[test]
    fn test_morton_locality() {
        let c1 = morton_encode(1, 1, 1);
        let c2 = morton_encode(2, 1, 1);
        let c_far = morton_encode(100, 100, 100);
        let d_near = c1.abs_diff(c2);
        let d_far = c1.abs_diff(c_far);
        assert!(
            d_near < d_far,
            "near distance {d_near} should be less than far {d_far}"
        );
    }
    /// Morton sort should produce a valid permutation.
    #[test]
    fn test_morton_sort_permutation() {
        let positions = vec![[5.0, 5.0, 5.0], [1.0, 1.0, 1.0], [3.0, 3.0, 3.0]];
        let (indices, codes) = morton_sort(&positions, [0.0; 3], [10.0, 10.0, 10.0]);
        assert_eq!(indices.len(), 3);
        assert_eq!(codes.len(), 3);
        for i in 0..codes.len() - 1 {
            assert!(codes[i] <= codes[i + 1], "codes not sorted at {i}");
        }
        let mut sorted = indices.clone();
        sorted.sort();
        assert_eq!(sorted, vec![0, 1, 2]);
    }
    /// Spatial hash should find nearby particles.
    #[test]
    fn test_spatial_hash_query() {
        let positions = vec![[0.5, 0.5, 0.5], [0.6, 0.5, 0.5], [5.0, 5.0, 5.0]];
        let mut hash = SpatialHash::new(64, 1.0);
        hash.build(&positions);
        assert_eq!(hash.len(), 3);
        let mut neighbours = hash.query_radius(&positions, [0.5, 0.5, 0.5], 0.5);
        neighbours.sort_unstable();
        neighbours.dedup();
        assert!(neighbours.contains(&0));
        assert!(neighbours.contains(&1));
        assert!(!neighbours.contains(&2));
    }
    /// Empty spatial hash.
    #[test]
    fn test_spatial_hash_empty() {
        let hash = SpatialHash::new(64, 1.0);
        assert!(hash.is_empty());
        assert_eq!(hash.len(), 0);
    }
    /// Spatial hash clear.
    #[test]
    fn test_spatial_hash_clear() {
        let mut hash = SpatialHash::new(64, 1.0);
        hash.insert(0, [0.5, 0.5, 0.5]);
        assert!(!hash.is_empty());
        hash.clear();
        assert!(hash.is_empty());
    }
    /// Bounding box computation.
    #[test]
    fn test_bounding_box() {
        let positions = vec![[1.0, 2.0, 3.0], [4.0, 0.0, 1.0], [2.0, 5.0, 2.0]];
        let (min, max) = compute_bounding_box(&positions);
        assert_eq!(min, [1.0, 0.0, 1.0]);
        assert_eq!(max, [4.0, 5.0, 3.0]);
    }
    /// Empty bounding box.
    #[test]
    fn test_bounding_box_empty() {
        let (min, max) = compute_bounding_box(&[]);
        assert_eq!(min, [0.0; 3]);
        assert_eq!(max, [0.0; 3]);
    }
    /// Reorder by permutation.
    #[test]
    fn test_reorder() {
        let data = vec![10, 20, 30, 40];
        let perm = vec![3, 1, 0, 2];
        let reordered = reorder_by_permutation(&data, &perm);
        assert_eq!(reordered, vec![40, 20, 10, 30]);
    }
    /// Max cell occupancy.
    #[test]
    fn test_max_cell_occupancy() {
        let positions: Vec<[f64; 3]> = vec![[0.5, 0.5, 0.5], [0.6, 0.5, 0.5], [5.0, 5.0, 5.0]];
        let cl = GpuCellList::build_parallel(&positions);
        let max = cl.max_cell_occupancy();
        assert!(max >= 2, "max occupancy should be at least 2, got {max}");
    }
    /// Non-empty cells count.
    #[test]
    fn test_nonempty_cells() {
        let positions: Vec<[f64; 3]> = vec![[0.5, 0.5, 0.5], [5.0, 5.0, 5.0]];
        let cl = GpuCellList::build_parallel(&positions);
        let ne = cl.num_nonempty_cells();
        assert_eq!(ne, 2, "should have 2 non-empty cells, got {ne}");
    }
    /// for_each_pair should find all pairs within cutoff.
    #[test]
    fn test_for_each_pair() {
        let positions: Vec<[f64; 3]> = vec![[0.5, 0.5, 0.5], [0.6, 0.5, 0.5], [5.0, 5.0, 5.0]];
        let cl = GpuCellList::build_parallel(&positions);
        let mut pairs = Vec::new();
        cl.for_each_pair(&positions, 0.5, |i, j, _d2| {
            pairs.push((i.min(j), i.max(j)));
        });
        pairs.sort();
        pairs.dedup();
        assert!(
            pairs.contains(&(0, 1)),
            "should find pair (0,1), got {pairs:?}"
        );
        assert!(
            !pairs.iter().any(|&(a, b)| a == 2 || b == 2),
            "should not find pairs with particle 2"
        );
    }
    #[test]
    fn test_radix_sort_sorted_output() {
        let keys = vec![5u32, 1, 9, 3, 7, 2];
        let (sorted_keys, sorted_indices) = radix_sort_mock(&keys);
        for i in 0..sorted_keys.len() - 1 {
            assert!(
                sorted_keys[i] <= sorted_keys[i + 1],
                "radix sort not sorted at {i}"
            );
        }
        for &idx in &sorted_indices {
            assert!(idx < keys.len(), "invalid index {idx}");
        }
    }
    #[test]
    fn test_radix_sort_permutation_correct() {
        let keys = vec![30u32, 10, 20];
        let (sorted_keys, sorted_indices) = radix_sort_mock(&keys);
        assert_eq!(sorted_keys[0], 10);
        assert_eq!(sorted_keys[1], 20);
        assert_eq!(sorted_keys[2], 30);
        assert_eq!(sorted_indices[0], 1);
        assert_eq!(sorted_indices[1], 2);
        assert_eq!(sorted_indices[2], 0);
    }
    #[test]
    fn test_radix_sort_empty() {
        let keys: Vec<u32> = vec![];
        let (sk, si) = radix_sort_mock(&keys);
        assert!(sk.is_empty());
        assert!(si.is_empty());
    }
    #[test]
    fn test_radix_sort_all_equal() {
        let keys = vec![7u32; 10];
        let (sorted_keys, sorted_indices) = radix_sort_mock(&keys);
        assert_eq!(sorted_keys.len(), 10);
        assert!(sorted_keys.iter().all(|&k| k == 7));
        assert_eq!(sorted_indices.len(), 10);
    }
    #[test]
    fn test_gpu_prefix_sum_basic() {
        let counts = vec![0usize, 1, 3, 0, 2, 5];
        let result = gpu_prefix_sum(&counts);
        assert_eq!(result, vec![0, 0, 1, 4, 4, 6]);
    }
    #[test]
    fn test_gpu_prefix_sum_all_zeros() {
        let counts = vec![0usize; 5];
        let result = gpu_prefix_sum(&counts);
        assert_eq!(result, vec![0, 0, 0, 0, 0]);
    }
    #[test]
    fn test_gpu_prefix_sum_single() {
        let result = gpu_prefix_sum(&[7usize]);
        assert_eq!(result, vec![0]);
    }
    #[test]
    fn test_parallel_cell_counting() {
        let positions = vec![[0.5, 0.5, 0.5], [0.5, 0.5, 0.5], [2.5, 0.5, 0.5]];
        let n_cells = [4usize, 4, 4];
        let counts = parallel_count_particles(&positions, n_cells, 1.0);
        let total: usize = counts.iter().sum();
        assert_eq!(total, 3, "total count should equal number of particles");
    }
    #[test]
    fn test_parallel_cell_counting_all_in_one_cell() {
        let positions: Vec<[f64; 3]> = vec![[0.1, 0.1, 0.1], [0.2, 0.1, 0.1], [0.1, 0.2, 0.1]];
        let counts = parallel_count_particles(&positions, [4, 4, 4], 1.0);
        let max = counts.iter().cloned().max().unwrap_or(0);
        assert!(max >= 3, "all particles in one cell: max_count={max}");
    }
    #[test]
    fn test_multi_gpu_distribution_two_gpus() {
        let n_cells = 100;
        let n_gpus = 2;
        let ranges = distribute_cells_to_gpus(n_cells, n_gpus);
        assert_eq!(ranges.len(), n_gpus);
        assert_eq!(ranges[0].start, 0);
        assert_eq!(ranges[n_gpus - 1].end, n_cells);
        for i in 0..n_gpus - 1 {
            assert_eq!(ranges[i].end, ranges[i + 1].start, "gap at gpu {i}");
        }
    }
    #[test]
    fn test_multi_gpu_distribution_odd_cells() {
        let ranges = distribute_cells_to_gpus(7, 3);
        assert_eq!(ranges.len(), 3);
        let total: usize = ranges.iter().map(|r| r.end - r.start).sum();
        assert_eq!(total, 7);
    }
    #[test]
    fn test_multi_gpu_distribution_single_gpu() {
        let ranges = distribute_cells_to_gpus(50, 1);
        assert_eq!(ranges.len(), 1);
        assert_eq!(ranges[0].start, 0);
        assert_eq!(ranges[0].end, 50);
    }
    #[test]
    fn test_neighbor_search_kernel_finds_close_pair() {
        let positions: Vec<[f64; 3]> = vec![[1.0, 1.0, 1.0], [1.1, 1.0, 1.0], [5.0, 5.0, 5.0]];
        let cl = GpuCellList::build_parallel(&positions);
        let pairs = gpu_neighbor_search_kernel(&cl, &positions, 0.5);
        assert!(
            pairs.contains(&(0, 1)) || pairs.contains(&(1, 0)),
            "should find pair (0,1), got {pairs:?}"
        );
        assert!(
            !pairs.iter().any(|&(a, b)| a == 2 || b == 2),
            "particle 2 should not appear in pairs"
        );
    }
    #[test]
    fn test_neighbor_search_kernel_no_pairs() {
        let positions: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [20.0, 0.0, 0.0]];
        let cl = GpuCellList::build_parallel(&positions);
        let pairs = gpu_neighbor_search_kernel(&cl, &positions, 0.5);
        assert!(pairs.is_empty(), "well-separated particles → no pairs");
    }
    #[test]
    fn test_spatial_hash_rebuild() {
        let mut hash = SpatialHash::new(128, 1.0);
        let positions1 = vec![[0.5, 0.5, 0.5], [1.5, 0.5, 0.5]];
        hash.build(&positions1);
        assert_eq!(hash.len(), 2);
        let positions2 = vec![[0.1, 0.1, 0.1]];
        hash.build(&positions2);
        assert_eq!(hash.len(), 1, "rebuild should replace old data");
    }
    #[test]
    fn test_spatial_hash_large_number_of_particles() {
        let positions: Vec<[f64; 3]> = (0..200).map(|i| [i as f64 * 0.1, 0.0, 0.0]).collect();
        let mut hash = SpatialHash::new(256, 1.0);
        hash.build(&positions);
        assert_eq!(hash.len(), 200);
    }
}
/// Sort particle indices by Morton code using Rayon for parallel code generation.
///
/// This variant uses Rayon for parallel code generation, then a sequential
/// sort for the final ordering (radix sort on a GPU would be fully parallel).
///
/// Returns `(sorted_indices, sorted_morton_codes)`.
pub fn parallel_morton_sort(
    positions: &[[f64; 3]],
    box_min: [f64; 3],
    box_max: [f64; 3],
) -> (Vec<usize>, Vec<u32>) {
    let range = [
        (box_max[0] - box_min[0]).max(1e-10),
        (box_max[1] - box_min[1]).max(1e-10),
        (box_max[2] - box_min[2]).max(1e-10),
    ];
    let mut code_index_pairs: Vec<(u32, usize)> = positions
        .par_iter()
        .enumerate()
        .map(|(i, p)| {
            let xi = (((p[0] - box_min[0]) / range[0]) * 1023.0) as u32;
            let yi = (((p[1] - box_min[1]) / range[1]) * 1023.0) as u32;
            let zi = (((p[2] - box_min[2]) / range[2]) * 1023.0) as u32;
            let x = xi.min(1023);
            let y = yi.min(1023);
            let z = zi.min(1023);
            (morton_encode(x, y, z), i)
        })
        .collect();
    code_index_pairs.sort_by_key(|&(code, _)| code);
    let sorted_indices: Vec<usize> = code_index_pairs.iter().map(|&(_, i)| i).collect();
    let sorted_codes: Vec<u32> = code_index_pairs.iter().map(|&(c, _)| c).collect();
    (sorted_indices, sorted_codes)
}
/// Compute the 30-bit Morton code for a position relative to a bounding box.
///
/// Normalises `pos` into `[0, 1023]^3` integer coordinates and encodes them.
pub fn position_to_morton(pos: [f64; 3], box_min: [f64; 3], box_max: [f64; 3]) -> u32 {
    let range = [
        (box_max[0] - box_min[0]).max(1e-10),
        (box_max[1] - box_min[1]).max(1e-10),
        (box_max[2] - box_min[2]).max(1e-10),
    ];
    let x = (((pos[0] - box_min[0]) / range[0]) * 1023.0) as u32;
    let y = (((pos[1] - box_min[1]) / range[1]) * 1023.0) as u32;
    let z = (((pos[2] - box_min[2]) / range[2]) * 1023.0) as u32;
    morton_encode(x.min(1023), y.min(1023), z.min(1023))
}
/// Insert particles into an existing cell list without a full rebuild.
///
/// This is a **sequential** incremental insert that updates `cell_counts` and
/// `sorted_indices` in place.  The `cell_starts` offsets remain valid only if
/// the new particles land in cells that already have capacity; otherwise a
/// full rebuild should be used.
///
/// Returns the number of particles successfully inserted.
pub fn insert_particles(cl: &mut GpuCellList, new_positions: &[[f64; 3]]) -> usize {
    let old_n = cl.sorted_indices.len();
    let mut inserted = 0usize;
    for (i, &pos) in new_positions.iter().enumerate() {
        let cell = cl.cell_index(pos);
        cl.sorted_indices.push(old_n + i);
        cl.cell_counts[cell] += 1;
        inserted += 1;
    }
    let new_starts = parallel_prefix_sum(
        &cl.cell_counts
            .iter()
            .map(|&c| c as usize)
            .collect::<Vec<_>>(),
    );
    cl.cell_starts = new_starts.iter().map(|&s| s as i32).collect();
    inserted
}
/// Query all particles within `radius` of `query_pos` from a set of positions,
/// using a pre-built `GpuCellList`.
///
/// This is an alias for `GpuCellList::neighbors_in_radius` exposed at module level
/// for ergonomic access.
pub fn query_neighbors(
    cl: &GpuCellList,
    positions: &[[f64; 3]],
    query_pos: [f64; 3],
    radius: f64,
) -> Vec<usize> {
    cl.neighbors_in_radius(positions, query_pos, radius)
}
#[cfg(test)]
mod extended_cell_tests {
    use crate::cell_list::CellList;
    use crate::cell_list::GhostCellManager;
    use crate::cell_list::GpuCellList;
    use crate::cell_list::GridResizer;
    use crate::cell_list::OccupancyStats;

    use crate::cell_list::insert_particles;
    use crate::cell_list::parallel_morton_sort;
    use crate::cell_list::position_to_morton;
    use crate::cell_list::query_neighbors;
    #[test]
    fn test_occupancy_stats_uniform() {
        let positions: Vec<[f64; 3]> = vec![
            [0.5, 0.5, 0.5],
            [1.5, 0.5, 0.5],
            [0.5, 1.5, 0.5],
            [1.5, 1.5, 0.5],
            [0.5, 0.5, 1.5],
            [1.5, 0.5, 1.5],
            [0.5, 1.5, 1.5],
            [1.5, 1.5, 1.5],
        ];
        let cl = GpuCellList::build_parallel(&positions);
        let stats = OccupancyStats::compute(&cl);
        assert_eq!(stats.total_particles, 8);
        assert_eq!(stats.max_occupancy, 1);
        assert!(stats.is_perfectly_spread());
    }
    #[test]
    fn test_occupancy_stats_clustered() {
        let positions: Vec<[f64; 3]> = vec![
            [0.1, 0.1, 0.1],
            [0.2, 0.1, 0.1],
            [0.1, 0.2, 0.1],
            [10.0, 10.0, 10.0],
        ];
        let cl = GpuCellList::build_parallel(&positions);
        let stats = OccupancyStats::compute(&cl);
        assert_eq!(stats.total_particles, 4);
        assert!(
            stats.max_occupancy >= 2,
            "clustered particles should share a cell"
        );
        assert_eq!(stats.nonempty_cells, 2);
    }
    #[test]
    fn test_occupancy_stats_load_imbalance_uniform() {
        let positions: Vec<[f64; 3]> = vec![[0.5, 0.5, 0.5], [1.5, 0.5, 0.5]];
        let cl = GpuCellList::build_parallel(&positions);
        let stats = OccupancyStats::compute(&cl);
        assert!(
            (stats.load_imbalance - 1.0).abs() < 1e-10,
            "load_imbalance = {}",
            stats.load_imbalance
        );
    }
    #[test]
    fn test_occupancy_stats_completely_unbalanced() {
        let positions: Vec<[f64; 3]> = vec![[0.1, 0.1, 0.1], [0.2, 0.2, 0.2], [0.3, 0.1, 0.1]];
        let cl = GpuCellList::build_parallel(&positions);
        let stats = OccupancyStats::compute(&cl);
        assert!(stats.is_completely_unbalanced() || stats.max_occupancy >= 2);
    }
    #[test]
    fn test_grid_resizer_initial_build() {
        let mut resizer = GridResizer::new(1.0, 0.5);
        let positions = vec![[1.0, 1.0, 1.0], [3.0, 3.0, 3.0]];
        resizer.update(&positions);
        assert!(resizer.get().is_some(), "cell list should be built");
    }
    #[test]
    fn test_grid_resizer_no_resize_needed() {
        let mut resizer = GridResizer::new(1.0, 1.0);
        let positions = vec![[2.0, 2.0, 2.0]];
        resizer.update(&positions);
        let needs = resizer.needs_resize(&positions);
        assert!(!needs, "same positions should not need resize");
    }
    #[test]
    fn test_grid_resizer_escaping_particle() {
        let mut resizer = GridResizer::new(1.0, 0.5);
        let positions = vec![[1.0, 1.0, 1.0]];
        resizer.rebuild(&positions);
        let new_positions = vec![[100.0, 100.0, 100.0]];
        assert!(
            resizer.needs_resize(&new_positions),
            "escaped particle should trigger resize"
        );
    }
    #[test]
    fn test_grid_resizer_empty_positions() {
        let mut resizer = GridResizer::new(1.0, 0.5);
        resizer.rebuild(&[]);
        assert!(
            resizer.get().is_some(),
            "empty rebuild should produce valid list"
        );
    }
    #[test]
    fn test_ghost_manager_no_ghosts_interior() {
        let mut mgr = GhostCellManager::new([10.0, 10.0, 10.0], 1.0);
        let positions = vec![[5.0, 5.0, 5.0], [6.0, 6.0, 6.0]];
        mgr.build_ghosts(&positions);
        assert_eq!(mgr.num_ghosts(), 0, "interior particles need no ghosts");
    }
    #[test]
    fn test_ghost_manager_near_one_face() {
        let mut mgr = GhostCellManager::new([10.0, 10.0, 10.0], 1.5);
        let positions = vec![[0.5, 5.0, 5.0]];
        mgr.build_ghosts(&positions);
        assert_eq!(mgr.num_ghosts(), 1, "should create 1 ghost on +x side");
        assert!(
            (mgr.ghost_positions[0][0] - 10.5).abs() < 1e-10,
            "ghost x = {}",
            mgr.ghost_positions[0][0]
        );
    }
    #[test]
    fn test_ghost_manager_near_two_faces() {
        let mut mgr = GhostCellManager::new([10.0, 10.0, 10.0], 1.5);
        let positions = vec![[0.5, 0.5, 5.0]];
        mgr.build_ghosts(&positions);
        assert_eq!(mgr.num_ghosts(), 2, "particle near two faces → 2 ghosts");
    }
    #[test]
    fn test_ghost_manager_near_corner() {
        let mut mgr = GhostCellManager::new([10.0, 10.0, 10.0], 1.5);
        let positions = vec![[0.5, 0.5, 0.5]];
        mgr.build_ghosts(&positions);
        assert_eq!(mgr.num_ghosts(), 3, "corner particle → 3 primary ghosts");
    }
    #[test]
    fn test_ghost_manager_map_to_real() {
        let mut mgr = GhostCellManager::new([10.0, 10.0, 10.0], 1.0);
        let positions = vec![[0.5, 5.0, 5.0], [9.5, 5.0, 5.0]];
        mgr.build_ghosts(&positions);
        for &ri in &mgr.ghost_to_real {
            assert!(ri < positions.len(), "real index {ri} out of range");
        }
    }
    #[test]
    fn test_minimum_image_convention() {
        let mgr = GhostCellManager::new([10.0, 10.0, 10.0], 1.0);
        let d = mgr.minimum_image([9.0, 0.0, 0.0]);
        assert!((d[0] - (-1.0)).abs() < 1e-10, "min image x = {}", d[0]);
        assert!(d[1].abs() < 1e-12);
        assert!(d[2].abs() < 1e-12);
    }
    #[test]
    fn test_wrap_position_basic() {
        let mgr = GhostCellManager::new([10.0, 10.0, 10.0], 1.0);
        let p = mgr.wrap_position([11.5, -0.5, 10.0]);
        assert!((p[0] - 1.5).abs() < 1e-10, "wrapped x = {}", p[0]);
        assert!((p[1] - 9.5).abs() < 1e-10, "wrapped y = {}", p[1]);
        assert!(p[2].abs() < 1e-10, "wrapped z = {}", p[2]);
    }
    #[test]
    fn test_wrap_all_in_place() {
        let mgr = GhostCellManager::new([5.0, 5.0, 5.0], 0.5);
        let mut positions = vec![[6.0, 7.0, 0.0], [-1.0, 2.5, 11.0]];
        mgr.wrap_all(&mut positions);
        for p in &positions {
            for k in 0..3 {
                assert!(
                    p[k] >= 0.0 && p[k] < 5.0,
                    "wrapped coord out of range: {}",
                    p[k]
                );
            }
        }
    }
    #[test]
    fn test_parallel_morton_sort_sorted_codes() {
        let positions = vec![
            [3.0, 3.0, 3.0],
            [1.0, 1.0, 1.0],
            [7.0, 7.0, 7.0],
            [5.0, 5.0, 5.0],
        ];
        let (_idx, codes) = parallel_morton_sort(&positions, [0.0; 3], [10.0; 3]);
        for i in 0..codes.len() - 1 {
            assert!(
                codes[i] <= codes[i + 1],
                "parallel morton sort codes not sorted at {i}"
            );
        }
    }
    #[test]
    fn test_parallel_morton_sort_valid_permutation() {
        let positions: Vec<[f64; 3]> = (0..10).map(|i| [i as f64, 0.0, 0.0]).collect();
        let (idx, codes) = parallel_morton_sort(&positions, [0.0; 3], [10.0, 1.0, 1.0]);
        assert_eq!(idx.len(), 10);
        assert_eq!(codes.len(), 10);
        let mut sorted_idx = idx.clone();
        sorted_idx.sort_unstable();
        assert_eq!(sorted_idx, (0..10).collect::<Vec<_>>());
    }
    #[test]
    fn test_position_to_morton_corner() {
        let code = position_to_morton([0.0, 0.0, 0.0], [0.0; 3], [1.0; 3]);
        assert_eq!(code, 0, "corner should give Morton code 0");
    }
    #[test]
    fn test_position_to_morton_different_positions() {
        let p1 = position_to_morton([1.0, 0.0, 0.0], [0.0; 3], [10.0; 3]);
        let p2 = position_to_morton([0.0, 1.0, 0.0], [0.0; 3], [10.0; 3]);
        let p3 = position_to_morton([5.0, 5.0, 5.0], [0.0; 3], [10.0; 3]);
        assert_ne!(p1, p3);
        assert_ne!(p2, p3);
    }
    #[test]
    fn test_insert_particles_increases_count() {
        let mut cl = GpuCellList::build_parallel(&[[1.0, 1.0, 1.0]]);
        let original_len = cl.sorted_indices.len();
        let new_particles = vec![[2.0, 2.0, 2.0], [3.0, 3.0, 3.0]];
        let inserted = insert_particles(&mut cl, &new_particles);
        assert_eq!(inserted, 2);
        assert_eq!(cl.sorted_indices.len(), original_len + 2);
    }
    #[test]
    fn test_insert_particles_empty_grid() {
        let mut cl = GpuCellList::new([4, 4, 4], 1.0, [4.0, 4.0, 4.0]);
        let positions = vec![[0.5, 0.5, 0.5], [1.5, 0.5, 0.5]];
        let inserted = insert_particles(&mut cl, &positions);
        assert_eq!(inserted, 2);
    }
    #[test]
    fn test_query_neighbors_finds_close_particle() {
        let positions = vec![[1.0, 1.0, 1.0], [1.1, 1.0, 1.0], [9.0, 9.0, 9.0]];
        let cl = GpuCellList::build_parallel(&positions);
        let mut neighbours = query_neighbors(&cl, &positions, [1.0, 1.0, 1.0], 0.5);
        neighbours.sort_unstable();
        assert!(neighbours.contains(&0), "should find self");
        assert!(neighbours.contains(&1), "should find nearby particle");
        assert!(!neighbours.contains(&2), "should not find far particle");
    }
    #[test]
    fn test_query_neighbors_empty_result() {
        let positions = vec![[0.0, 0.0, 0.0], [100.0, 100.0, 100.0]];
        let cl = GpuCellList::build_parallel(&positions);
        let neighbours = query_neighbors(&cl, &positions, [50.0, 50.0, 50.0], 0.1);
        assert!(neighbours.is_empty(), "no particles near middle of box");
    }
    #[test]
    fn test_verlet_list_close_pair_found() {
        let positions = vec![[0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [10.0, 10.0, 10.0]];
        let cl = CellList::build(&positions);
        let pairs = cl.build_neighbor_list_verlet(1.0, 0.2);
        let has_01 = pairs.contains(&(0, 1));
        assert!(has_01, "pair (0,1) must be in Verlet list");
    }
    #[test]
    fn test_verlet_list_far_pair_excluded() {
        let positions = vec![[0.0, 0.0, 0.0], [20.0, 20.0, 20.0]];
        let cl = CellList::build(&positions);
        let pairs = cl.build_neighbor_list_verlet(1.0, 0.2);
        assert!(pairs.is_empty(), "far pair must not appear in Verlet list");
    }
    #[test]
    fn test_verlet_list_no_self_pairs() {
        let positions = vec![[1.0, 1.0, 1.0], [1.1, 1.0, 1.0], [1.2, 1.0, 1.0]];
        let cl = CellList::build(&positions);
        let pairs = cl.build_neighbor_list_verlet(1.0, 0.5);
        for &(i, j) in &pairs {
            assert_ne!(i, j, "self-pair found");
        }
    }
    #[test]
    fn test_verlet_list_pairs_ordered() {
        let positions: Vec<[f64; 3]> = (0..5).map(|i| [i as f64 * 0.3, 0.0, 0.0]).collect();
        let cl = CellList::build(&positions);
        let pairs = cl.build_neighbor_list_verlet(1.0, 0.1);
        for &(i, j) in &pairs {
            assert!(i < j, "Verlet pair must have i < j");
        }
    }
    #[test]
    fn test_update_incremental_no_move() {
        let positions = vec![[1.0, 1.0, 1.0], [3.0, 3.0, 3.0]];
        let mut cl = CellList::build(&positions);
        let relocated = cl.update_incremental(&positions, &positions, 0.1);
        assert_eq!(relocated, 0, "no particle moved");
    }
    #[test]
    fn test_update_incremental_large_move_counted() {
        let old = vec![[1.0, 1.0, 1.0], [3.0, 3.0, 3.0]];
        let new_pos = vec![[1.0, 1.0, 1.0], [6.0, 6.0, 6.0]];
        let mut cl = CellList::build(&old);
        let relocated = cl.update_incremental(&new_pos, &old, 0.5);
        assert!(relocated >= 1, "at least one particle relocated");
    }
    #[test]
    fn test_update_incremental_threshold_respected() {
        let old = vec![[0.0, 0.0, 0.0], [5.0, 5.0, 5.0]];
        let new_pos = vec![[0.05, 0.0, 0.0], [8.0, 8.0, 8.0]];
        let mut cl = CellList::build(&old);
        let relocated = cl.update_incremental(&new_pos, &old, 1.0);
        assert_eq!(relocated, 1);
    }
    #[test]
    fn test_pair_density_single_bin() {
        let positions = vec![[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]];
        let cl = CellList::build(&positions);
        let hist = cl.compute_pair_density(2.0, 1.0);
        assert!(hist[0] >= 1, "pair must appear in bin 0");
    }
    #[test]
    fn test_pair_density_no_pairs_beyond_max_r() {
        let positions = vec![[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]];
        let cl = CellList::build(&positions);
        let hist = cl.compute_pair_density(2.0, 0.5);
        let total: usize = hist.iter().sum();
        assert_eq!(total, 0, "pair beyond max_r should not be counted");
    }
    #[test]
    fn test_pair_density_histogram_length() {
        let positions = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let cl = CellList::build(&positions);
        let hist = cl.compute_pair_density(5.0, 1.0);
        assert_eq!(hist.len(), 5, "histogram length = ceil(max_r/dr)");
    }
}
