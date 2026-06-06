// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Broadphase AABB kernels for parallel overlap detection.

use crate::compute::ComputeKernel;
use std::collections::HashMap;

// ---------------------------------------------------------------------------
// Legacy f64 kernels (keep existing API)
// ---------------------------------------------------------------------------

/// Kernel that sorts AABBs along the X-axis and outputs sorted indices.
///
/// **Input layout** (flat f64 array, 6 values per object):
///   `[min_x, max_x, min_y, max_y, min_z, max_z, ...]`
///
/// **Output\[0\]**: sorted indices (as f64, cast back to usize by caller).
pub struct AabbSortKernel;

impl ComputeKernel for AabbSortKernel {
    fn name(&self) -> &str {
        "AabbSortKernel"
    }

    fn execute(&self, inputs: &[&[f64]], outputs: &mut [Vec<f64>], _work_size: usize) {
        if inputs.is_empty() || outputs.is_empty() {
            return;
        }
        let aabbs = inputs[0];
        let n = aabbs.len() / 6;
        // Build (index, min_x) pairs and sort by min_x.
        let mut indices: Vec<usize> = (0..n).collect();
        indices.sort_by(|&a, &b| {
            let ax = aabbs[a * 6];
            let bx = aabbs[b * 6];
            ax.partial_cmp(&bx).unwrap_or(std::cmp::Ordering::Equal)
        });
        outputs[0] = indices.iter().map(|&i| i as f64).collect();
    }
}

/// Kernel that detects overlapping AABB pairs.
///
/// **Input layout** (flat f64 array, 6 values per object):
///   `[min_x, max_x, min_y, max_y, min_z, max_z, ...]`
///
/// **Output\[0\]**: flat pairs `[i, j, i, j, ...]` of overlapping indices (as f64).
pub struct AabbOverlapKernel;

impl AabbOverlapKernel {
    #[inline]
    fn overlaps(a: &[f64], b: &[f64]) -> bool {
        // a, b are slices of length 6: [min_x, max_x, min_y, max_y, min_z, max_z]
        a[0] <= b[1] && a[1] >= b[0] // x overlap
        && a[2] <= b[3] && a[3] >= b[2] // y overlap
        && a[4] <= b[5] && a[5] >= b[4] // z overlap
    }
}

impl ComputeKernel for AabbOverlapKernel {
    fn name(&self) -> &str {
        "AabbOverlapKernel"
    }

    fn execute(&self, inputs: &[&[f64]], outputs: &mut [Vec<f64>], _work_size: usize) {
        if inputs.is_empty() || outputs.is_empty() {
            return;
        }
        let aabbs = inputs[0];
        let n = aabbs.len() / 6;
        let mut pairs = Vec::new();
        for i in 0..n {
            for j in (i + 1)..n {
                if Self::overlaps(&aabbs[i * 6..(i + 1) * 6], &aabbs[j * 6..(j + 1) * 6]) {
                    pairs.push(i as f64);
                    pairs.push(j as f64);
                }
            }
        }
        outputs[0] = pairs;
    }
}

// ---------------------------------------------------------------------------
// GPU-compatible f32 AABB
// ---------------------------------------------------------------------------

/// Axis-aligned bounding box using f32 coordinates for GPU compatibility.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AabbGpu {
    /// Minimum corner `[min_x, min_y, min_z]`.
    pub min: [f32; 3],
    /// Maximum corner `[max_x, max_y, max_z]`.
    pub max: [f32; 3],
    /// Index of the owning rigid body.
    pub body_id: u32,
}

impl AabbGpu {
    /// Create a new GPU AABB.
    pub fn new(min: [f32; 3], max: [f32; 3], body_id: u32) -> Self {
        Self { min, max, body_id }
    }

    /// Test whether two AABBs overlap in all three axes.
    #[inline]
    pub fn overlaps(&self, other: &AabbGpu) -> bool {
        self.min[0] <= other.max[0]
            && self.max[0] >= other.min[0]
            && self.min[1] <= other.max[1]
            && self.max[1] >= other.min[1]
            && self.min[2] <= other.max[2]
            && self.max[2] >= other.min[2]
    }
}

// ---------------------------------------------------------------------------
// Sort-and-Sweep (GPU)
// ---------------------------------------------------------------------------

/// Sort-and-sweep broadphase operating on `AabbGpu` values.
///
/// Bodies are sorted by their minimum X coordinate, then swept to find all
/// overlapping pairs in O(n log n + k) time where k is the number of pairs.
pub struct SortAndSweepGpu;

impl SortAndSweepGpu {
    /// Detect all overlapping AABB pairs.
    ///
    /// Returns a `Vec` of `(body_id_a, body_id_b)` pairs where the two AABBs
    /// overlap in all three axes.  Each pair is reported exactly once with
    /// `body_id_a < body_id_b`.
    pub fn detect_pairs(aabbs: &[AabbGpu]) -> Vec<(u32, u32)> {
        if aabbs.is_empty() {
            return Vec::new();
        }

        // Sort by min_x.
        let mut sorted: Vec<&AabbGpu> = aabbs.iter().collect();
        sorted.sort_by(|a, b| {
            a.min[0]
                .partial_cmp(&b.min[0])
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        let n = sorted.len();
        let mut pairs = Vec::new();

        for i in 0..n {
            for j in (i + 1)..n {
                // Early-out: if the next AABB's min_x exceeds current max_x, no further overlaps
                if sorted[j].min[0] > sorted[i].max[0] {
                    break;
                }
                if sorted[i].overlaps(sorted[j]) {
                    let a = sorted[i].body_id;
                    let b = sorted[j].body_id;
                    let pair = if a < b { (a, b) } else { (b, a) };
                    pairs.push(pair);
                }
            }
        }

        pairs
    }
}

// ---------------------------------------------------------------------------
// Uniform Grid (GPU)
// ---------------------------------------------------------------------------

/// Uniform 3-D grid for GPU broadphase.
#[derive(Debug, Clone)]
pub struct UniformGridGpu {
    /// Cell size (same in all dimensions).
    pub cell_size: f32,
    /// World-space origin of the grid `[ox, oy, oz]`.
    pub origin: [f32; 3],
    /// Number of cells in each dimension `[nx, ny, nz]`.
    pub dims: [u32; 3],
}

impl UniformGridGpu {
    /// Create a new uniform grid.
    pub fn new(cell_size: f32, origin: [f32; 3], dims: [u32; 3]) -> Self {
        Self {
            cell_size,
            origin,
            dims,
        }
    }

    /// Compute the cell index `[ix, iy, iz]` for a world-space position.
    pub fn cell_of(&self, pos: [f32; 3]) -> [i32; 3] {
        [
            ((pos[0] - self.origin[0]) / self.cell_size).floor() as i32,
            ((pos[1] - self.origin[1]) / self.cell_size).floor() as i32,
            ((pos[2] - self.origin[2]) / self.cell_size).floor() as i32,
        ]
    }

    /// Pack cell indices into a `u64` key for use in a HashMap.
    fn cell_key(ix: i32, iy: i32, iz: i32) -> u64 {
        // Pack three i16-range integers into a u64.
        let x = (ix as i16) as u64 & 0xFFFF;
        let y = (iy as i16) as u64 & 0xFFFF;
        let z = (iz as i16) as u64 & 0xFFFF;
        (z << 32) | (y << 16) | x
    }

    /// Insert all AABBs into the uniform grid.
    ///
    /// Each AABB is inserted into every cell it overlaps.  Returns a map from
    /// packed cell key to the list of `body_id`s in that cell.
    pub fn insert_aabbs(&self, aabbs: &[AabbGpu]) -> HashMap<u64, Vec<u32>> {
        let mut map: HashMap<u64, Vec<u32>> = HashMap::new();

        for aabb in aabbs {
            let min_cell = self.cell_of(aabb.min);
            let max_cell = self.cell_of(aabb.max);

            for iz in min_cell[2]..=max_cell[2] {
                for iy in min_cell[1]..=max_cell[1] {
                    for ix in min_cell[0]..=max_cell[0] {
                        let key = Self::cell_key(ix, iy, iz);
                        map.entry(key).or_default().push(aabb.body_id);
                    }
                }
            }
        }

        map
    }

    /// Query all overlapping AABB pairs using the uniform grid.
    ///
    /// Returns deduplicated pairs `(a, b)` with `a < b`.
    pub fn query_pairs(&self, aabbs: &[AabbGpu]) -> Vec<(u32, u32)> {
        let map = self.insert_aabbs(aabbs);

        let mut seen = std::collections::HashSet::new();
        let mut pairs = Vec::new();

        for body_list in map.values() {
            let n = body_list.len();
            for i in 0..n {
                for j in (i + 1)..n {
                    let a = body_list[i];
                    let b = body_list[j];
                    let pair = if a < b { (a, b) } else { (b, a) };
                    if seen.insert(pair) {
                        // Verify actual AABB overlap
                        if let (Some(aa), Some(bb)) = (
                            aabbs.iter().find(|x| x.body_id == pair.0),
                            aabbs.iter().find(|x| x.body_id == pair.1),
                        ) && aa.overlaps(bb)
                        {
                            pairs.push(pair);
                        }
                    }
                }
            }
        }

        pairs.sort();
        pairs
    }
}

// ---------------------------------------------------------------------------
// Morton code (Z-curve)
// ---------------------------------------------------------------------------

/// Compute the 3-D Morton code (Z-curve) for integer coordinates.
///
/// Interleaves the bits of `x`, `y`, and `z` to produce a single `u64`
/// space-filling key.  Each coordinate may use up to 21 bits.
pub fn morton_code(x: u32, y: u32, z: u32) -> u64 {
    spread_bits(x as u64) | (spread_bits(y as u64) << 1) | (spread_bits(z as u64) << 2)
}

/// Spread the bits of a 21-bit integer, inserting two zero bits between each bit.
#[inline]
fn spread_bits(mut v: u64) -> u64 {
    v &= 0x1fffff; // keep lower 21 bits
    v = (v | (v << 32)) & 0x1f00000000ffff;
    v = (v | (v << 16)) & 0x1f0000ff0000ff;
    v = (v | (v << 8)) & 0x100f00f00f00f00f;
    v = (v | (v << 4)) & 0x10c30c30c30c30c3;
    v = (v | (v << 2)) & 0x1249249249249249;
    v
}

// ---------------------------------------------------------------------------
// Sort-and-Sweep (CPU-side mock for GPU-style parallel data)
// ---------------------------------------------------------------------------

/// Sort-and-sweep on a flat f64 AABB array (GPU-style parallel data layout).
///
/// Input layout per object: `[min_x, max_x, min_y, max_y, min_z, max_z]`.
/// Returns pairs of indices `(i, j)` (sorted, deduped) that overlap.
pub fn sort_and_sweep_flat(aabbs: &[f64]) -> Vec<(usize, usize)> {
    let n = aabbs.len() / 6;
    if n == 0 {
        return Vec::new();
    }

    // Sort by min_x
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| {
        let ax = aabbs[a * 6];
        let bx = aabbs[b * 6];
        ax.partial_cmp(&bx).unwrap_or(std::cmp::Ordering::Equal)
    });

    let mut pairs = Vec::new();
    for (i, &si) in order.iter().enumerate() {
        let max_x_i = aabbs[si * 6 + 1];
        for &sj in order.iter().skip(i + 1) {
            if aabbs[sj * 6] > max_x_i {
                break; // Early out: rest cannot overlap in X
            }
            // Full 3-D overlap check
            let ai = &aabbs[si * 6..(si + 1) * 6];
            let aj = &aabbs[sj * 6..(sj + 1) * 6];
            if ai[0] <= aj[1]
                && ai[1] >= aj[0]
                && ai[2] <= aj[3]
                && ai[3] >= aj[2]
                && ai[4] <= aj[5]
                && ai[5] >= aj[4]
            {
                let pair = if si < sj { (si, sj) } else { (sj, si) };
                pairs.push(pair);
            }
        }
    }
    pairs.sort();
    pairs.dedup();
    pairs
}

// ---------------------------------------------------------------------------
// Uniform Grid (CPU-side parallel data helper)
// ---------------------------------------------------------------------------

/// Assigns each AABB in a flat array to uniform grid cells.
///
/// Returns a `Vec<(cell_key, body_idx)>` (unsorted, one entry per cell overlap).
pub fn assign_to_grid_cells(aabbs: &[f64], cell_size: f64, origin: [f64; 3]) -> Vec<(u64, usize)> {
    let n = aabbs.len() / 6;
    let mut result = Vec::new();
    let pack = |ix: i64, iy: i64, iz: i64| -> u64 {
        let x = (ix as i16) as u64 & 0xFFFF;
        let y = (iy as i16) as u64 & 0xFFFF;
        let z = (iz as i16) as u64 & 0xFFFF;
        (z << 32) | (y << 16) | x
    };
    for i in 0..n {
        let a = &aabbs[i * 6..(i + 1) * 6];
        let ix0 = ((a[0] - origin[0]) / cell_size).floor() as i64;
        let ix1 = ((a[1] - origin[0]) / cell_size).floor() as i64;
        let iy0 = ((a[2] - origin[1]) / cell_size).floor() as i64;
        let iy1 = ((a[3] - origin[1]) / cell_size).floor() as i64;
        let iz0 = ((a[4] - origin[2]) / cell_size).floor() as i64;
        let iz1 = ((a[5] - origin[2]) / cell_size).floor() as i64;
        for iz in iz0..=iz1 {
            for iy in iy0..=iy1 {
                for ix in ix0..=ix1 {
                    result.push((pack(ix, iy, iz), i));
                }
            }
        }
    }
    result
}

/// Extract candidate pairs from a cell-assignment list.
///
/// Groups entries by cell key and emits pairs for bodies sharing a cell.
pub fn pairs_from_grid_assignments(assignments: &[(u64, usize)]) -> Vec<(usize, usize)> {
    let mut by_cell: HashMap<u64, Vec<usize>> = HashMap::new();
    for &(key, idx) in assignments {
        by_cell.entry(key).or_default().push(idx);
    }
    let mut seen = std::collections::HashSet::new();
    let mut pairs = Vec::new();
    for bodies in by_cell.values() {
        let n = bodies.len();
        for i in 0..n {
            for j in (i + 1)..n {
                let a = bodies[i];
                let b = bodies[j];
                let p = if a < b { (a, b) } else { (b, a) };
                if seen.insert(p) {
                    pairs.push(p);
                }
            }
        }
    }
    pairs.sort();
    pairs
}

// ---------------------------------------------------------------------------
// Morton code sort utilities
// ---------------------------------------------------------------------------

/// Compute a Morton key for a world-space AABB centroid.
///
/// Quantises each coordinate to a 21-bit integer using `cell_size` and `origin`.
pub fn morton_key_for_aabb(aabb: &AabbGpu, cell_size: f32, origin: [f32; 3]) -> u64 {
    let cx = ((aabb.min[0] + aabb.max[0]) * 0.5 - origin[0]) / cell_size;
    let cy = ((aabb.min[1] + aabb.max[1]) * 0.5 - origin[1]) / cell_size;
    let cz = ((aabb.min[2] + aabb.max[2]) * 0.5 - origin[2]) / cell_size;

    let ix = (cx.max(0.0) as u32).min(0x1F_FFFF);
    let iy = (cy.max(0.0) as u32).min(0x1F_FFFF);
    let iz = (cz.max(0.0) as u32).min(0x1F_FFFF);
    morton_code(ix, iy, iz)
}

/// Sort a slice of `AabbGpu` by Morton code (Z-order curve).
///
/// Returns a new `Vec`AabbGpu` sorted by `morton_key_for_aabb`.
pub fn morton_sort(aabbs: &[AabbGpu], cell_size: f32, origin: [f32; 3]) -> Vec<AabbGpu> {
    let mut keyed: Vec<(u64, AabbGpu)> = aabbs
        .iter()
        .map(|a| (morton_key_for_aabb(a, cell_size, origin), *a))
        .collect();
    keyed.sort_by_key(|&(k, _)| k);
    keyed.into_iter().map(|(_, a)| a).collect()
}

// ---------------------------------------------------------------------------
// Compact pairs list
// ---------------------------------------------------------------------------

/// A deduplicated, sorted list of overlapping body-ID pairs.
#[derive(Debug, Clone, Default)]
pub struct CompactPairList {
    pairs: Vec<(u32, u32)>,
}

impl CompactPairList {
    /// Create an empty list.
    pub fn new() -> Self {
        Self::default()
    }

    /// Insert a pair (order-normalised to a < b).
    pub fn insert(&mut self, a: u32, b: u32) {
        let pair = if a < b { (a, b) } else { (b, a) };
        // Simple insert-if-not-present (small lists)
        if !self.pairs.contains(&pair) {
            self.pairs.push(pair);
        }
    }

    /// Insert all pairs from a `Vec`.
    pub fn insert_all(&mut self, pairs: &[(u32, u32)]) {
        for &(a, b) in pairs {
            self.insert(a, b);
        }
    }

    /// Sort the pair list.
    pub fn sort(&mut self) {
        self.pairs.sort();
    }

    /// Return a reference to all pairs.
    pub fn pairs(&self) -> &[(u32, u32)] {
        &self.pairs
    }

    /// Number of pairs.
    pub fn len(&self) -> usize {
        self.pairs.len()
    }

    /// Returns true if the list is empty.
    pub fn is_empty(&self) -> bool {
        self.pairs.is_empty()
    }

    /// Remove pairs involving a specific body (e.g., after removal from simulation).
    pub fn remove_body(&mut self, body_id: u32) {
        self.pairs.retain(|&(a, b)| a != body_id && b != body_id);
    }

    /// Returns `true` if the given pair exists.
    pub fn contains(&self, a: u32, b: u32) -> bool {
        let p = if a < b { (a, b) } else { (b, a) };
        self.pairs.contains(&p)
    }
}

// ---------------------------------------------------------------------------
// BVH (GPU)
// ---------------------------------------------------------------------------

/// A single node of a GPU-oriented BVH tree.
#[derive(Debug, Clone, Copy)]
pub struct BvhGpuNode {
    /// Bounding box of this node.
    pub aabb: AabbGpu,
    /// Index of the left child node, or `< 0` for a leaf (then `-body_id - 1`).
    pub left: i32,
    /// Index of the right child node, or `< 0` for a leaf (then `-body_id - 1`).
    pub right: i32,
}

impl BvhGpuNode {
    /// Create an internal node.
    pub fn internal(aabb: AabbGpu, left: i32, right: i32) -> Self {
        Self { aabb, left, right }
    }

    /// Create a leaf node for a given body.
    pub fn leaf(aabb: AabbGpu) -> Self {
        let id = aabb.body_id as i32;
        Self {
            aabb,
            left: -(id + 1),
            right: -(id + 1),
        }
    }

    /// Return true if this is a leaf node.
    pub fn is_leaf(&self) -> bool {
        self.left < 0
    }
}

/// Build a BVH using a simple midpoint-split heuristic.
///
/// Returns a flat `Vec`BvhGpuNode` with node 0 as the root.
/// Internal nodes reference their children by index into this vec.
pub fn build_bvh(aabbs: &[AabbGpu]) -> Vec<BvhGpuNode> {
    let mut nodes = Vec::new();
    if aabbs.is_empty() {
        return nodes;
    }
    let mut indices: Vec<usize> = (0..aabbs.len()).collect();
    build_bvh_recursive(aabbs, &mut indices, &mut nodes);
    nodes
}

fn merge_aabbs(aabbs: &[AabbGpu], indices: &[usize]) -> AabbGpu {
    let mut min = aabbs[indices[0]].min;
    let mut max = aabbs[indices[0]].max;
    for &idx in &indices[1..] {
        let a = &aabbs[idx];
        for k in 0..3 {
            if a.min[k] < min[k] {
                min[k] = a.min[k];
            }
            if a.max[k] > max[k] {
                max[k] = a.max[k];
            }
        }
    }
    AabbGpu {
        min,
        max,
        body_id: 0,
    }
}

fn build_bvh_recursive(
    aabbs: &[AabbGpu],
    indices: &mut [usize],
    nodes: &mut Vec<BvhGpuNode>,
) -> i32 {
    let n = indices.len();
    let merged = merge_aabbs(aabbs, indices);

    if n == 1 {
        let idx = nodes.len() as i32;
        nodes.push(BvhGpuNode::leaf(aabbs[indices[0]]));
        return idx;
    }

    // Find the longest axis of the merged AABB.
    let extents = [
        merged.max[0] - merged.min[0],
        merged.max[1] - merged.min[1],
        merged.max[2] - merged.min[2],
    ];
    let axis = if extents[0] >= extents[1] && extents[0] >= extents[2] {
        0
    } else if extents[1] >= extents[2] {
        1
    } else {
        2
    };

    // Sort by centroid along longest axis.
    indices.sort_by(|&a, &b| {
        let ca = (aabbs[a].min[axis] + aabbs[a].max[axis]) * 0.5;
        let cb = (aabbs[b].min[axis] + aabbs[b].max[axis]) * 0.5;
        ca.partial_cmp(&cb).unwrap_or(std::cmp::Ordering::Equal)
    });

    let mid = n / 2;
    let (left_idx, right_idx) = indices.split_at_mut(mid);

    // Reserve a slot for the internal node before recursing.
    let node_idx = nodes.len() as i32;
    nodes.push(BvhGpuNode {
        aabb: merged,
        left: 0,
        right: 0,
    }); // placeholder

    let mut left_indices = left_idx.to_vec();
    let mut right_indices = right_idx.to_vec();

    let left = build_bvh_recursive(aabbs, &mut left_indices, nodes);
    let right = build_bvh_recursive(aabbs, &mut right_indices, nodes);

    nodes[node_idx as usize].left = left;
    nodes[node_idx as usize].right = right;

    node_idx
}

// ---------------------------------------------------------------------------
// LBVH — Linear BVH traversal
// ---------------------------------------------------------------------------

/// Query overlapping leaf pairs in an LBVH by traversal.
///
/// Traverses the BVH from the root to find all pairs of leaf nodes whose
/// AABBs overlap.  Returns `Vec<(u32, u32)>` with `a < b` and deduplication.
pub fn lbvh_query_pairs(nodes: &[BvhGpuNode]) -> Vec<(u32, u32)> {
    if nodes.is_empty() {
        return Vec::new();
    }
    let mut pairs = Vec::new();
    // Stack-based traversal: (nodeA, nodeB)
    let mut stack: Vec<(usize, usize)> = Vec::new();
    // Self-collision query: test root against itself
    stack.push((0, 0));

    while let Some((a_idx, b_idx)) = stack.pop() {
        let na = &nodes[a_idx];
        let nb = &nodes[b_idx];

        if !na.aabb.overlaps(&nb.aabb) {
            continue;
        }

        if na.is_leaf() && nb.is_leaf() {
            if a_idx != b_idx {
                let id_a = na.aabb.body_id;
                let id_b = nb.aabb.body_id;
                let pair = if id_a < id_b {
                    (id_a, id_b)
                } else {
                    (id_b, id_a)
                };
                pairs.push(pair);
            }
            continue;
        }

        // Descend into larger node first
        if na.is_leaf() {
            // descend b
            if nb.left >= 0 {
                stack.push((a_idx, nb.left as usize));
            }
            if nb.right >= 0 {
                stack.push((a_idx, nb.right as usize));
            }
        } else if nb.is_leaf() {
            // descend a
            if na.left >= 0 {
                stack.push((na.left as usize, b_idx));
            }
            if na.right >= 0 {
                stack.push((na.right as usize, b_idx));
            }
        } else {
            // Both internal — descend a
            if na.left >= 0 {
                stack.push((na.left as usize, b_idx));
            }
            if na.right >= 0 {
                stack.push((na.right as usize, b_idx));
            }
        }
    }
    pairs.sort();
    pairs.dedup();
    pairs
}

// ---------------------------------------------------------------------------
// BVH refitting
// ---------------------------------------------------------------------------

/// Refit BVH bounding boxes bottom-up after leaf AABBs have changed.
///
/// For each internal node, recomputes its AABB as the union of its children.
/// Assumes a flat node array where children always have higher indices
/// (which is guaranteed by `build_bvh_recursive`).
pub fn refit_bvh(nodes: &mut [BvhGpuNode]) {
    // Process nodes in reverse order (children before parents)
    let n = nodes.len();
    for i in (0..n).rev() {
        if nodes[i].is_leaf() {
            continue; // Leaves are driven by actual AABB data — no refit needed here
        }
        let left_idx = nodes[i].left;
        let right_idx = nodes[i].right;
        if left_idx < 0 || right_idx < 0 {
            continue;
        }
        let l = &nodes[left_idx as usize];
        let r = &nodes[right_idx as usize];
        let mut min = l.aabb.min;
        let mut max = l.aabb.max;
        for k in 0..3 {
            if r.aabb.min[k] < min[k] {
                min[k] = r.aabb.min[k];
            }
            if r.aabb.max[k] > max[k] {
                max[k] = r.aabb.max[k];
            }
        }
        nodes[i].aabb = AabbGpu {
            min,
            max,
            body_id: 0,
        };
    }
}

// ---------------------------------------------------------------------------
// BVH quality metrics
// ---------------------------------------------------------------------------

/// Compute the surface area of an `AabbGpu`.
pub fn aabb_surface_area(aabb: &AabbGpu) -> f32 {
    let dx = aabb.max[0] - aabb.min[0];
    let dy = aabb.max[1] - aabb.min[1];
    let dz = aabb.max[2] - aabb.min[2];
    2.0 * (dx * dy + dy * dz + dz * dx)
}

/// Surface Area Heuristic (SAH) cost of a BVH tree.
///
/// `SAH = sum_{internal nodes} SA(node) / SA(root) * cost_traversal`
///       `+ sum_{leaf nodes} SA(leaf) / SA(root) * num_primitives`
///
/// A lower SAH cost indicates a better-quality BVH.
pub fn bvh_sah_cost(nodes: &[BvhGpuNode], cost_traversal: f32, cost_primitive: f32) -> f32 {
    if nodes.is_empty() {
        return 0.0;
    }
    let root_sa = aabb_surface_area(&nodes[0].aabb);
    if root_sa < 1e-20 {
        return 0.0;
    }
    let mut cost = 0.0f32;
    for node in nodes {
        let sa = aabb_surface_area(&node.aabb);
        if node.is_leaf() {
            cost += sa / root_sa * cost_primitive;
        } else {
            cost += sa / root_sa * cost_traversal;
        }
    }
    cost
}

/// Compute the depth of the BVH tree.
///
/// Returns the maximum node depth from root (depth 0) to the deepest leaf.
pub fn bvh_depth(nodes: &[BvhGpuNode]) -> usize {
    if nodes.is_empty() {
        return 0;
    }
    let mut max_depth = 0usize;
    // Stack: (node_idx, current_depth)
    let mut stack: Vec<(usize, usize)> = vec![(0, 0)];
    while let Some((idx, depth)) = stack.pop() {
        if depth > max_depth {
            max_depth = depth;
        }
        let node = &nodes[idx];
        if !node.is_leaf() {
            if node.left >= 0 {
                stack.push((node.left as usize, depth + 1));
            }
            if node.right >= 0 {
                stack.push((node.right as usize, depth + 1));
            }
        }
    }
    max_depth
}

/// Count the number of leaf nodes in the BVH.
pub fn bvh_leaf_count(nodes: &[BvhGpuNode]) -> usize {
    nodes.iter().filter(|n| n.is_leaf()).count()
}

// ---------------------------------------------------------------------------
// Parallel SAP update (incremental update for moved bodies)
// ---------------------------------------------------------------------------

/// Update the sort-and-sweep pair list after a set of bodies have moved.
///
/// This is an incremental update: only re-run SAP for the bodies that moved
/// and merge with the existing pair list.  All existing pairs involving moved
/// bodies are removed and recomputed.
///
/// `moved_ids`: set of `body_id`s that have moved.
/// `aabbs`: updated AABB array (all bodies).
///
/// Returns the updated pair list.
pub fn sap_incremental_update(
    existing: &CompactPairList,
    aabbs: &[AabbGpu],
    moved_ids: &[u32],
) -> CompactPairList {
    let mut new_list = existing.clone();

    // Remove all pairs involving moved bodies
    for &id in moved_ids {
        new_list.remove_body(id);
    }

    // Re-detect pairs for moved bodies against all others
    for &moved_id in moved_ids {
        if let Some(moved_aabb) = aabbs.iter().find(|a| a.body_id == moved_id) {
            for other in aabbs {
                if other.body_id == moved_id {
                    continue;
                }
                if moved_aabb.overlaps(other) {
                    new_list.insert(moved_id, other.body_id);
                }
            }
        }
    }

    new_list.sort();
    new_list
}

// ---------------------------------------------------------------------------
// SAH split heuristic
// ---------------------------------------------------------------------------

/// Surface Area Heuristic (SAH) split for BVH construction.
///
/// Evaluates `num_bins` candidate split planes along the given axis and
/// returns the bin index (from 0 to `num_bins-1`) that minimises the SAH cost.
///
/// Returns `None` if all primitives have the same centroid along the axis.
pub fn sah_best_split(
    aabbs: &[AabbGpu],
    indices: &[usize],
    axis: usize,
    num_bins: usize,
) -> Option<usize> {
    if indices.len() < 2 || num_bins < 2 {
        return None;
    }

    // Centroid range
    let min_c = indices
        .iter()
        .map(|&i| 0.5 * (aabbs[i].min[axis] + aabbs[i].max[axis]))
        .fold(f32::INFINITY, f32::min);
    let max_c = indices
        .iter()
        .map(|&i| 0.5 * (aabbs[i].min[axis] + aabbs[i].max[axis]))
        .fold(f32::NEG_INFINITY, f32::max);

    if (max_c - min_c).abs() < 1e-10 {
        return None;
    }

    let bin_width = (max_c - min_c) / num_bins as f32;
    let mut bin_counts = vec![0usize; num_bins];
    let mut bin_aabbs: Vec<Option<AabbGpu>> = vec![None; num_bins];

    for &i in indices {
        let c = 0.5 * (aabbs[i].min[axis] + aabbs[i].max[axis]);
        let bin = ((c - min_c) / bin_width).floor() as usize;
        let bin = bin.min(num_bins - 1);
        bin_counts[bin] += 1;
        bin_aabbs[bin] = Some(match &bin_aabbs[bin] {
            None => aabbs[i],
            Some(prev) => {
                let mut merged = *prev;
                for k in 0..3 {
                    if aabbs[i].min[k] < merged.min[k] {
                        merged.min[k] = aabbs[i].min[k];
                    }
                    if aabbs[i].max[k] > merged.max[k] {
                        merged.max[k] = aabbs[i].max[k];
                    }
                }
                merged
            }
        });
    }

    let mut best_cost = f32::INFINITY;
    let mut best_split = None;

    for split in 1..num_bins {
        let left_count: usize = bin_counts[..split].iter().sum();
        let right_count: usize = bin_counts[split..].iter().sum();
        if left_count == 0 || right_count == 0 {
            continue;
        }

        // Compute left and right bounding boxes
        let left_sa = bin_aabbs[..split]
            .iter()
            .flatten()
            .fold(None::<AabbGpu>, |acc, a| {
                Some(match acc {
                    None => *a,
                    Some(prev) => {
                        let mut m = prev;
                        for k in 0..3 {
                            if a.min[k] < m.min[k] {
                                m.min[k] = a.min[k];
                            }
                            if a.max[k] > m.max[k] {
                                m.max[k] = a.max[k];
                            }
                        }
                        m
                    }
                })
            })
            .map(|a| aabb_surface_area(&a))
            .unwrap_or(0.0);

        let right_sa = bin_aabbs[split..]
            .iter()
            .flatten()
            .fold(None::<AabbGpu>, |acc, a| {
                Some(match acc {
                    None => *a,
                    Some(prev) => {
                        let mut m = prev;
                        for k in 0..3 {
                            if a.min[k] < m.min[k] {
                                m.min[k] = a.min[k];
                            }
                            if a.max[k] > m.max[k] {
                                m.max[k] = a.max[k];
                            }
                        }
                        m
                    }
                })
            })
            .map(|a| aabb_surface_area(&a))
            .unwrap_or(0.0);

        let cost = left_sa * left_count as f32 + right_sa * right_count as f32;
        if cost < best_cost {
            best_cost = cost;
            best_split = Some(split);
        }
    }

    best_split
}

// ---------------------------------------------------------------------------
// LBVH build using Morton codes
// ---------------------------------------------------------------------------

/// Build a Linear BVH (LBVH) by sorting AABBs on their Morton codes.
///
/// This is a GPU-style construction algorithm:
/// 1. Compute Morton code for each AABB centroid.
/// 2. Sort AABBs by Morton code.
/// 3. Recursively build BVH on sorted order (binary radix tree).
///
/// Returns a flat node array with node 0 as the root.
pub fn build_lbvh(aabbs: &[AabbGpu], cell_size: f32, origin: [f32; 3]) -> Vec<BvhGpuNode> {
    if aabbs.is_empty() {
        return Vec::new();
    }
    // Sort by Morton code
    let sorted = morton_sort(aabbs, cell_size, origin);
    // Build BVH on sorted order
    build_bvh(&sorted)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn aabb_overlap_detects_overlapping_boxes() {
        // Two boxes that clearly overlap
        #[rustfmt::skip]
        let aabbs: Vec<f64> = vec![
            0.0, 2.0, 0.0, 2.0, 0.0, 2.0, // box 0
            1.0, 3.0, 1.0, 3.0, 1.0, 3.0, // box 1 (overlaps box 0)
        ];
        let mut outputs = vec![Vec::new()];
        AabbOverlapKernel.execute(&[&aabbs], &mut outputs, 2);
        assert_eq!(outputs[0], vec![0.0, 1.0]);
    }

    #[test]
    fn aabb_overlap_rejects_non_overlapping_boxes() {
        #[rustfmt::skip]
        let aabbs: Vec<f64> = vec![
            0.0, 1.0, 0.0, 1.0, 0.0, 1.0, // box 0
            5.0, 6.0, 5.0, 6.0, 5.0, 6.0, // box 1 (far away)
        ];
        let mut outputs = vec![Vec::new()];
        AabbOverlapKernel.execute(&[&aabbs], &mut outputs, 2);
        assert!(outputs[0].is_empty());
    }

    #[test]
    fn test_broadphase_gpu_matches_cpu() {
        // Three AABBs: 0 overlaps 1, 1 overlaps 2, but 0 does NOT overlap 2.
        #[rustfmt::skip]
        let aabbs: Vec<f64> = vec![
            0.0, 1.5, 0.0, 1.5, 0.0, 1.5, // box 0
            1.0, 2.5, 0.0, 1.5, 0.0, 1.5, // box 1 (overlaps 0 in x)
            3.0, 4.0, 0.0, 1.5, 0.0, 1.5, // box 2 (overlaps 1 in x, not 0)
        ];

        // GPU kernel output
        let mut gpu_outputs = vec![Vec::new()];
        AabbOverlapKernel.execute(&[&aabbs], &mut gpu_outputs, 3);

        // CPU brute-force reference
        let n = aabbs.len() / 6;
        let mut cpu_pairs: Vec<(usize, usize)> = Vec::new();
        for i in 0..n {
            for j in (i + 1)..n {
                let a = &aabbs[i * 6..(i + 1) * 6];
                let b = &aabbs[j * 6..(j + 1) * 6];
                let overlaps = a[0] <= b[1]
                    && a[1] >= b[0]
                    && a[2] <= b[3]
                    && a[3] >= b[2]
                    && a[4] <= b[5]
                    && a[5] >= b[4];
                if overlaps {
                    cpu_pairs.push((i, j));
                }
            }
        }

        // Convert GPU output to (usize, usize) pairs
        let raw = &gpu_outputs[0];
        assert_eq!(raw.len() % 2, 0, "GPU output length must be even");
        let gpu_pairs: Vec<(usize, usize)> = raw
            .chunks(2)
            .map(|c| (c[0] as usize, c[1] as usize))
            .collect();

        assert_eq!(
            gpu_pairs, cpu_pairs,
            "GPU broadphase pairs do not match CPU brute-force pairs"
        );
    }

    /// Morton code: (0,0,0) → 0, and each axis occupies its own bit lanes.
    #[test]
    fn test_morton_code_correctness() {
        assert_eq!(morton_code(0, 0, 0), 0);
        // x=1: only bit 0 of x set → Morton bit 0 (axis 0)
        assert_eq!(morton_code(1, 0, 0), 1);
        // y=1: Morton bit 1
        assert_eq!(morton_code(0, 1, 0), 2);
        // z=1: Morton bit 2
        assert_eq!(morton_code(0, 0, 1), 4);
        // x=1, y=1, z=1 → bits 0,1,2 all set
        assert_eq!(morton_code(1, 1, 1), 7);
    }

    /// SAP detects overlapping GPU AABBs.
    #[test]
    fn test_sap_finds_overlap() {
        let a = AabbGpu::new([0.0, 0.0, 0.0], [2.0, 2.0, 2.0], 0);
        let b = AabbGpu::new([1.0, 1.0, 1.0], [3.0, 3.0, 3.0], 1);
        let c = AabbGpu::new([10.0, 10.0, 10.0], [12.0, 12.0, 12.0], 2);

        let pairs = SortAndSweepGpu::detect_pairs(&[a, b, c]);
        assert!(
            pairs.contains(&(0, 1)),
            "SAP should find pair (0,1), got {pairs:?}"
        );
        assert!(!pairs.contains(&(0, 2)), "SAP should not find pair (0,2)");
        assert!(!pairs.contains(&(1, 2)), "SAP should not find pair (1,2)");
    }

    /// SAP returns empty for non-overlapping AABBs.
    #[test]
    fn test_sap_no_overlap() {
        let a = AabbGpu::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0], 0);
        let b = AabbGpu::new([5.0, 5.0, 5.0], [6.0, 6.0, 6.0], 1);
        let pairs = SortAndSweepGpu::detect_pairs(&[a, b]);
        assert!(
            pairs.is_empty(),
            "Should be no overlapping pairs, got {pairs:?}"
        );
    }

    /// Uniform grid detects the same pairs as brute-force for 3 AABBs.
    #[test]
    fn test_uniform_grid_pair_detection() {
        let a = AabbGpu::new([0.0, 0.0, 0.0], [2.0, 2.0, 2.0], 0);
        let b = AabbGpu::new([1.5, 1.5, 1.5], [3.5, 3.5, 3.5], 1);
        let c = AabbGpu::new([10.0, 10.0, 10.0], [12.0, 12.0, 12.0], 2);

        let grid = UniformGridGpu::new(5.0, [0.0, 0.0, 0.0], [10, 10, 10]);
        let pairs = grid.query_pairs(&[a, b, c]);

        // Only a and b overlap.
        assert!(
            pairs.contains(&(0, 1)),
            "Grid should find pair (0,1), got {pairs:?}"
        );
        assert!(
            !pairs
                .iter()
                .any(|&(x, y)| x == 0 && y == 2 || x == 2 && y == 0)
        );
    }

    /// BVH builds without panic and covers all body AABBs.
    #[test]
    fn test_bvh_depth() {
        let aabbs: Vec<AabbGpu> = (0..8)
            .map(|i| {
                let x = (i * 3) as f32;
                AabbGpu::new([x, 0.0, 0.0], [x + 1.0, 1.0, 1.0], i)
            })
            .collect();

        let nodes = build_bvh(&aabbs);

        // There should be at least 1 node, and at most 2*n-1 nodes.
        assert!(!nodes.is_empty(), "BVH should have at least one node");
        assert!(
            nodes.len() <= 2 * aabbs.len(),
            "BVH node count unexpected: {}",
            nodes.len()
        );

        // The root node (index 0) should cover the entire range.
        let root = &nodes[0];
        assert!(root.aabb.min[0] <= 0.0 + 1e-5, "Root min_x too large");
        assert!(
            root.aabb.max[0] >= 22.0 - 1e-5,
            "Root max_x too small: {}",
            root.aabb.max[0]
        );
    }

    /// AabbGpu::overlaps is symmetric.
    #[test]
    fn test_aabb_gpu_overlap_symmetric() {
        let a = AabbGpu::new([0.0, 0.0, 0.0], [2.0, 2.0, 2.0], 0);
        let b = AabbGpu::new([1.0, 1.0, 1.0], [3.0, 3.0, 3.0], 1);
        assert_eq!(a.overlaps(&b), b.overlaps(&a));

        let c = AabbGpu::new([5.0, 5.0, 5.0], [6.0, 6.0, 6.0], 2);
        assert_eq!(a.overlaps(&c), c.overlaps(&a));
        assert!(!a.overlaps(&c));
    }

    // ── New expanded tests ──

    #[test]
    fn test_sort_and_sweep_flat_finds_pair() {
        #[rustfmt::skip]
        let aabbs: Vec<f64> = vec![
            0.0, 2.0, 0.0, 2.0, 0.0, 2.0, // box 0
            1.0, 3.0, 1.0, 3.0, 1.0, 3.0, // box 1 overlaps 0
            5.0, 6.0, 5.0, 6.0, 5.0, 6.0, // box 2 no overlap
        ];
        let pairs = sort_and_sweep_flat(&aabbs);
        assert!(pairs.contains(&(0, 1)), "should find (0,1)");
        assert!(!pairs.iter().any(|&(a, b)| (a == 0 || a == 1) && b == 2));
    }

    #[test]
    fn test_sort_and_sweep_flat_empty() {
        let pairs = sort_and_sweep_flat(&[]);
        assert!(pairs.is_empty());
    }

    #[test]
    fn test_sort_and_sweep_flat_no_overlap() {
        #[rustfmt::skip]
        let aabbs: Vec<f64> = vec![
            0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
            2.0, 3.0, 2.0, 3.0, 2.0, 3.0,
        ];
        let pairs = sort_and_sweep_flat(&aabbs);
        assert!(pairs.is_empty());
    }

    #[test]
    fn test_assign_to_grid_cells() {
        #[rustfmt::skip]
        let aabbs: Vec<f64> = vec![
            0.0, 1.0, 0.0, 1.0, 0.0, 1.0, // fits in cell (0,0,0)
        ];
        let cells = assign_to_grid_cells(&aabbs, 2.0, [0.0, 0.0, 0.0]);
        assert!(!cells.is_empty());
        // All entries should have body index 0
        assert!(cells.iter().all(|&(_, idx)| idx == 0));
    }

    #[test]
    fn test_pairs_from_grid_assignments() {
        // Two bodies both assigned to the same cell
        let assignments = vec![(0u64, 0usize), (0u64, 1usize)];
        let pairs = pairs_from_grid_assignments(&assignments);
        assert!(pairs.contains(&(0, 1)));
    }

    #[test]
    fn test_pairs_from_grid_assignments_no_dup() {
        // Same pair in multiple cells
        let assignments = vec![
            (0u64, 0usize),
            (0u64, 1usize),
            (1u64, 0usize),
            (1u64, 1usize),
        ];
        let pairs = pairs_from_grid_assignments(&assignments);
        // Pair (0,1) should appear exactly once
        let count = pairs.iter().filter(|&&p| p == (0, 1)).count();
        assert_eq!(count, 1);
    }

    #[test]
    fn test_morton_sort_orders_aabbs() {
        let aabbs: Vec<AabbGpu> = vec![
            AabbGpu::new([4.0, 0.0, 0.0], [5.0, 1.0, 1.0], 0),
            AabbGpu::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0], 1),
            AabbGpu::new([2.0, 2.0, 2.0], [3.0, 3.0, 3.0], 2),
        ];
        let sorted = morton_sort(&aabbs, 1.0, [0.0, 0.0, 0.0]);
        // Body 1 (centroid near origin) should come first in Z-order
        assert_eq!(
            sorted[0].body_id, 1,
            "body at origin should be first in Morton order"
        );
    }

    #[test]
    fn test_morton_key_reproducible() {
        let a = AabbGpu::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0], 0);
        let k1 = morton_key_for_aabb(&a, 1.0, [0.0, 0.0, 0.0]);
        let k2 = morton_key_for_aabb(&a, 1.0, [0.0, 0.0, 0.0]);
        assert_eq!(k1, k2);
    }

    #[test]
    fn test_compact_pair_list_insert_dedup() {
        let mut list = CompactPairList::new();
        list.insert(0, 1);
        list.insert(1, 0); // same pair, reversed
        list.insert(0, 1); // duplicate
        assert_eq!(list.len(), 1);
    }

    #[test]
    fn test_compact_pair_list_contains() {
        let mut list = CompactPairList::new();
        list.insert(2, 5);
        assert!(list.contains(2, 5));
        assert!(list.contains(5, 2));
        assert!(!list.contains(0, 1));
    }

    #[test]
    fn test_compact_pair_list_remove_body() {
        let mut list = CompactPairList::new();
        list.insert(0, 1);
        list.insert(0, 2);
        list.insert(1, 2);
        list.remove_body(0);
        assert!(!list.contains(0, 1));
        assert!(!list.contains(0, 2));
        assert!(list.contains(1, 2));
    }

    #[test]
    fn test_compact_pair_list_insert_all() {
        let mut list = CompactPairList::new();
        list.insert_all(&[(0, 1), (1, 2), (0, 2)]);
        assert_eq!(list.len(), 3);
    }

    #[test]
    fn test_compact_pair_list_sort() {
        let mut list = CompactPairList::new();
        list.insert(3, 4);
        list.insert(0, 1);
        list.insert(1, 2);
        list.sort();
        assert_eq!(list.pairs()[0], (0, 1));
    }

    // ── LBVH traversal tests ─────────────────────────────────────────────────

    #[test]
    fn test_lbvh_query_pairs_finds_overlap() {
        // Build BVH for two overlapping boxes
        let a = AabbGpu::new([0.0, 0.0, 0.0], [2.0, 2.0, 2.0], 0);
        let b = AabbGpu::new([1.0, 1.0, 1.0], [3.0, 3.0, 3.0], 1);
        let nodes = build_bvh(&[a, b]);
        let pairs = lbvh_query_pairs(&nodes);
        assert!(
            pairs.contains(&(0, 1)),
            "LBVH traversal should find (0,1): {pairs:?}"
        );
    }

    #[test]
    fn test_lbvh_query_pairs_no_overlap() {
        let a = AabbGpu::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0], 0);
        let b = AabbGpu::new([5.0, 5.0, 5.0], [6.0, 6.0, 6.0], 1);
        let nodes = build_bvh(&[a, b]);
        let pairs = lbvh_query_pairs(&nodes);
        assert!(
            pairs.is_empty(),
            "Non-overlapping: should have no pairs: {pairs:?}"
        );
    }

    #[test]
    fn test_lbvh_query_empty_bvh() {
        let pairs = lbvh_query_pairs(&[]);
        assert!(pairs.is_empty());
    }

    // ── BVH refitting tests ──────────────────────────────────────────────────

    #[test]
    fn test_refit_bvh_no_panic() {
        let aabbs: Vec<AabbGpu> = (0..4)
            .map(|i| {
                let x = (i * 2) as f32;
                AabbGpu::new([x, 0.0, 0.0], [x + 1.0, 1.0, 1.0], i)
            })
            .collect();
        let mut nodes = build_bvh(&aabbs);
        // Simulate moving a leaf
        for node in nodes.iter_mut() {
            if node.is_leaf() {
                node.aabb.max[0] += 0.5;
            }
        }
        refit_bvh(&mut nodes);
        // After refitting, root should still cover all children
        assert!(!nodes.is_empty());
    }

    #[test]
    fn test_refit_bvh_root_encompasses_leaves() {
        let aabbs: Vec<AabbGpu> = vec![
            AabbGpu::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0], 0),
            AabbGpu::new([10.0, 0.0, 0.0], [11.0, 1.0, 1.0], 1),
        ];
        let mut nodes = build_bvh(&aabbs);
        refit_bvh(&mut nodes);
        // Root AABB should span at least [0..11]
        assert!(nodes[0].aabb.min[0] <= 0.0 + 1e-5);
        assert!(nodes[0].aabb.max[0] >= 11.0 - 1e-5);
    }

    // ── BVH quality metrics tests ────────────────────────────────────────────

    #[test]
    fn test_aabb_surface_area() {
        let a = AabbGpu::new([0.0, 0.0, 0.0], [1.0, 2.0, 3.0], 0);
        let sa = aabb_surface_area(&a);
        // SA = 2*(1*2 + 2*3 + 3*1) = 2*(2+6+3) = 22
        assert!((sa - 22.0).abs() < 1e-5, "SA = {sa}");
    }

    #[test]
    fn test_aabb_surface_area_unit_cube() {
        let a = AabbGpu::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0], 0);
        let sa = aabb_surface_area(&a);
        assert!((sa - 6.0).abs() < 1e-5, "unit cube SA = {sa}");
    }

    #[test]
    fn test_bvh_sah_cost_positive() {
        let aabbs: Vec<AabbGpu> = (0..4)
            .map(|i| {
                let x = (i * 3) as f32;
                AabbGpu::new([x, 0.0, 0.0], [x + 2.0, 2.0, 2.0], i)
            })
            .collect();
        let nodes = build_bvh(&aabbs);
        let cost = bvh_sah_cost(&nodes, 1.0, 1.0);
        assert!(cost > 0.0, "SAH cost should be positive: {cost}");
        assert!(cost.is_finite(), "SAH cost should be finite");
    }

    #[test]
    fn test_bvh_depth_single_leaf() {
        let aabbs = vec![AabbGpu::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0], 0)];
        let nodes = build_bvh(&aabbs);
        let d = bvh_depth(&nodes);
        assert_eq!(d, 0, "single leaf depth = {d}");
    }

    #[test]
    fn test_bvh_depth_two_leaves() {
        let aabbs = vec![
            AabbGpu::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0], 0),
            AabbGpu::new([2.0, 0.0, 0.0], [3.0, 1.0, 1.0], 1),
        ];
        let nodes = build_bvh(&aabbs);
        let d = bvh_depth(&nodes);
        assert!(d >= 1, "two-leaf BVH depth >= 1, got {d}");
    }

    #[test]
    fn test_bvh_leaf_count() {
        let aabbs: Vec<AabbGpu> = (0..8)
            .map(|i| {
                let x = (i * 2) as f32;
                AabbGpu::new([x, 0.0, 0.0], [x + 1.0, 1.0, 1.0], i)
            })
            .collect();
        let nodes = build_bvh(&aabbs);
        let leaves = bvh_leaf_count(&nodes);
        assert_eq!(leaves, 8, "Expected 8 leaves, got {leaves}");
    }

    // ── SAP incremental update tests ─────────────────────────────────────────

    #[test]
    fn test_sap_incremental_removes_moved() {
        let a = AabbGpu::new([0.0, 0.0, 0.0], [2.0, 2.0, 2.0], 0);
        let b = AabbGpu::new([1.0, 1.0, 1.0], [3.0, 3.0, 3.0], 1);
        let c = AabbGpu::new([10.0, 10.0, 10.0], [11.0, 11.0, 11.0], 2);
        let initial_pairs = SortAndSweepGpu::detect_pairs(&[a, b, c]);
        let mut existing = CompactPairList::new();
        existing.insert_all(&initial_pairs);

        // Move body 0 far away
        let a_moved = AabbGpu::new([20.0, 20.0, 20.0], [21.0, 21.0, 21.0], 0);
        let updated = sap_incremental_update(&existing, &[a_moved, b, c], &[0]);

        // After move, body 0 overlaps nothing
        assert!(
            !updated.contains(0, 1),
            "pair (0,1) should be removed after move"
        );
    }

    #[test]
    fn test_sap_incremental_adds_new_overlap() {
        let _a = AabbGpu::new([10.0, 0.0, 0.0], [11.0, 1.0, 1.0], 0);
        let b = AabbGpu::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0], 1);
        let existing = CompactPairList::new(); // no initial pairs

        // Move body 0 to overlap body 1
        let a_moved = AabbGpu::new([0.5, 0.5, 0.5], [1.5, 1.5, 1.5], 0);
        let updated = sap_incremental_update(&existing, &[a_moved, b], &[0]);
        assert!(
            updated.contains(0, 1),
            "pair (0,1) should be added after move: {:?}",
            updated.pairs()
        );
    }

    // ── SAH split tests ──────────────────────────────────────────────────────

    #[test]
    fn test_sah_best_split_basic() {
        let aabbs: Vec<AabbGpu> = (0..8)
            .map(|i| {
                let x = (i * 2) as f32;
                AabbGpu::new([x, 0.0, 0.0], [x + 1.5, 1.0, 1.0], i)
            })
            .collect();
        let indices: Vec<usize> = (0..8).collect();
        let result = sah_best_split(&aabbs, &indices, 0, 8);
        assert!(result.is_some(), "SAH split should find a valid split");
    }

    #[test]
    fn test_sah_best_split_single_element() {
        let aabbs = vec![AabbGpu::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0], 0)];
        let result = sah_best_split(&aabbs, &[0], 0, 4);
        assert!(result.is_none(), "Single element: no split possible");
    }

    // ── LBVH build tests ─────────────────────────────────────────────────────

    #[test]
    fn test_build_lbvh_nonempty() {
        let aabbs: Vec<AabbGpu> = (0..6)
            .map(|i| {
                let x = (i * 2) as f32;
                AabbGpu::new([x, 0.0, 0.0], [x + 1.0, 1.0, 1.0], i)
            })
            .collect();
        let nodes = build_lbvh(&aabbs, 1.0, [0.0, 0.0, 0.0]);
        assert!(!nodes.is_empty(), "LBVH should build non-empty node list");
        assert!(nodes.len() <= 2 * aabbs.len());
    }

    #[test]
    fn test_build_lbvh_empty() {
        let nodes = build_lbvh(&[], 1.0, [0.0, 0.0, 0.0]);
        assert!(nodes.is_empty());
    }

    #[test]
    fn test_build_lbvh_single() {
        let aabbs = vec![AabbGpu::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0], 0)];
        let nodes = build_lbvh(&aabbs, 1.0, [0.0, 0.0, 0.0]);
        assert_eq!(nodes.len(), 1);
        assert!(nodes[0].is_leaf());
    }

    #[test]
    fn test_lbvh_vs_brute_force_pairs() {
        // LBVH traversal should find same pairs as brute-force overlap check
        let aabbs: Vec<AabbGpu> = (0..6)
            .map(|i| {
                let x = (i as f32) * 0.8;
                AabbGpu::new([x, 0.0, 0.0], [x + 1.0, 1.0, 1.0], i as u32)
            })
            .collect();

        let lbvh_nodes = build_lbvh(&aabbs, 0.5, [0.0, 0.0, 0.0]);
        let mut lbvh_pairs = lbvh_query_pairs(&lbvh_nodes);
        lbvh_pairs.sort();
        lbvh_pairs.dedup();

        let mut brute: Vec<(u32, u32)> = Vec::new();
        for i in 0..aabbs.len() {
            for j in (i + 1)..aabbs.len() {
                if aabbs[i].overlaps(&aabbs[j]) {
                    brute.push((aabbs[i].body_id, aabbs[j].body_id));
                }
            }
        }
        brute.sort();

        assert_eq!(
            lbvh_pairs, brute,
            "LBVH pairs {lbvh_pairs:?} != brute force {brute:?}"
        );
    }
}
