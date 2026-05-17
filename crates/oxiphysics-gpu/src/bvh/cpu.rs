// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! CPU BVH construction, query, traversal, and LBVH helpers.

use super::types::{
    Aabb, BvhNode, BvhPrimitive, BvhStats, BvhTreeStatistics, FlatBvhNode, LbvhPrimitive,
    MortonCluster, RayHit,
};

// ============================================================================
// SAH helper
// ============================================================================

/// Surface Area Heuristic cost:
/// `C = (SA_left / SA_parent) * N_left + (SA_right / SA_parent) * N_right`
pub fn sah_cost(n_left: usize, sa_left: f32, n_right: usize, sa_right: f32, sa_parent: f32) -> f32 {
    if sa_parent <= 0.0 {
        return f32::MAX;
    }
    (sa_left / sa_parent) * n_left as f32 + (sa_right / sa_parent) * n_right as f32
}

// ============================================================================
// Slab-method ray–AABB intersection
// ============================================================================

/// Test whether a ray defined by `origin` + t * direction intersects `aabb`
/// within `[0, max_t]`.
///
/// `inv_dir` must be the component-wise reciprocal of the ray direction.
pub fn ray_aabb_intersect(origin: [f32; 3], inv_dir: [f32; 3], aabb: &Aabb, max_t: f32) -> bool {
    let mut t_min = 0.0_f32;
    let mut t_max = max_t;

    for i in 0..3 {
        let t1 = (aabb.min[i] - origin[i]) * inv_dir[i];
        let t2 = (aabb.max[i] - origin[i]) * inv_dir[i];
        let lo = t1.min(t2);
        let hi = t1.max(t2);
        t_min = t_min.max(lo);
        t_max = t_max.min(hi);
    }

    t_min <= t_max
}

/// Ray–AABB slab intersection returning the near/far t values.
pub(crate) fn ray_aabb_t(origin: [f32; 3], inv_dir: [f32; 3], aabb: &Aabb) -> Option<(f32, f32)> {
    let mut t_min = 0.0_f32;
    let mut t_max = f32::MAX;
    for i in 0..3 {
        let t1 = (aabb.min[i] - origin[i]) * inv_dir[i];
        let t2 = (aabb.max[i] - origin[i]) * inv_dir[i];
        t_min = t_min.max(t1.min(t2));
        t_max = t_max.min(t1.max(t2));
    }
    if t_min <= t_max {
        Some((t_min, t_max))
    } else {
        None
    }
}

// ============================================================================
// Bvh
// ============================================================================

/// Maximum number of primitives per leaf before the tree stops splitting.
pub(crate) const LEAF_SIZE: usize = 4;

/// A BVH tree built from a flat list of [`BvhPrimitive`]s.
pub struct Bvh {
    /// Tree root (if there is at least one primitive).
    pub root: Option<BvhNode>,
    /// All primitives passed to [`Bvh::build`].
    pub primitives: Vec<BvhPrimitive>,
}

impl Bvh {
    /// Build a BVH from a list of primitives using a median-split strategy
    /// guided by the longest axis (SAH-inspired).
    pub fn build(primitives: Vec<BvhPrimitive>) -> Self {
        if primitives.is_empty() {
            return Self {
                root: None,
                primitives,
            };
        }
        let indices: Vec<usize> = (0..primitives.len()).collect();
        let root = build_recursive(&primitives, indices);
        Self {
            root: Some(root),
            primitives,
        }
    }

    /// Return the `object_id`s of all primitives whose AABB overlaps `query`.
    pub fn query_aabb(&self, query: &Aabb) -> Vec<usize> {
        let mut result = Vec::new();
        if let Some(root) = &self.root {
            query_aabb_recursive(root, query, &self.primitives, &mut result);
        }
        result
    }

    /// Return the `object_id`s of all primitives hit by the given ray within
    /// distance `max_t`.
    pub fn query_ray(&self, origin: [f32; 3], direction: [f32; 3], max_t: f32) -> Vec<usize> {
        let inv_dir = [1.0 / direction[0], 1.0 / direction[1], 1.0 / direction[2]];
        let mut result = Vec::new();
        if let Some(root) = &self.root {
            query_ray_recursive(root, origin, inv_dir, max_t, &self.primitives, &mut result);
        }
        result
    }

    /// Total number of nodes in the tree.
    pub fn node_count(&self) -> usize {
        match &self.root {
            None => 0,
            Some(root) => count_nodes(root),
        }
    }

    /// Maximum depth of the tree (root = depth 1).
    pub fn depth(&self) -> usize {
        match &self.root {
            None => 0,
            Some(root) => node_depth(root),
        }
    }
}

// ============================================================================
// Internal build / query helpers
// ============================================================================

pub(crate) fn bounding_box(primitives: &[BvhPrimitive], indices: &[usize]) -> Aabb {
    let mut aabb = primitives[indices[0]].aabb.clone();
    for &i in &indices[1..] {
        aabb = Aabb::merge(&aabb, &primitives[i].aabb);
    }
    aabb
}

fn build_recursive(primitives: &[BvhPrimitive], mut indices: Vec<usize>) -> BvhNode {
    let aabb = bounding_box(primitives, &indices);

    if indices.len() <= LEAF_SIZE {
        return BvhNode {
            aabb,
            left: None,
            right: None,
            primitives: indices,
        };
    }

    // Choose longest axis for split.
    let dx = aabb.max[0] - aabb.min[0];
    let dy = aabb.max[1] - aabb.min[1];
    let dz = aabb.max[2] - aabb.min[2];
    let axis = if dx >= dy && dx >= dz {
        0
    } else if dy >= dz {
        1
    } else {
        2
    };

    // Median split by centre of primitive AABB along the chosen axis.
    indices.sort_unstable_by(|&a, &b| {
        let ca = primitives[a].aabb.center()[axis];
        let cb = primitives[b].aabb.center()[axis];
        ca.partial_cmp(&cb).unwrap_or(std::cmp::Ordering::Equal)
    });

    let mid = indices.len() / 2;
    let right_indices = indices.split_off(mid);
    let left_indices = indices;

    let left = build_recursive(primitives, left_indices);
    let right = build_recursive(primitives, right_indices);

    BvhNode {
        aabb,
        left: Some(Box::new(left)),
        right: Some(Box::new(right)),
        primitives: Vec::new(),
    }
}

fn query_aabb_recursive(
    node: &BvhNode,
    query: &Aabb,
    primitives: &[BvhPrimitive],
    result: &mut Vec<usize>,
) {
    if !node.aabb.intersects(query) {
        return;
    }
    if node.is_leaf() {
        for &idx in &node.primitives {
            if primitives[idx].aabb.intersects(query) {
                result.push(primitives[idx].object_id);
            }
        }
    } else {
        if let Some(left) = &node.left {
            query_aabb_recursive(left, query, primitives, result);
        }
        if let Some(right) = &node.right {
            query_aabb_recursive(right, query, primitives, result);
        }
    }
}

fn query_ray_recursive(
    node: &BvhNode,
    origin: [f32; 3],
    inv_dir: [f32; 3],
    max_t: f32,
    primitives: &[BvhPrimitive],
    result: &mut Vec<usize>,
) {
    if !ray_aabb_intersect(origin, inv_dir, &node.aabb, max_t) {
        return;
    }
    if node.is_leaf() {
        for &idx in &node.primitives {
            if ray_aabb_intersect(origin, inv_dir, &primitives[idx].aabb, max_t) {
                result.push(primitives[idx].object_id);
            }
        }
    } else {
        if let Some(left) = &node.left {
            query_ray_recursive(left, origin, inv_dir, max_t, primitives, result);
        }
        if let Some(right) = &node.right {
            query_ray_recursive(right, origin, inv_dir, max_t, primitives, result);
        }
    }
}

fn count_nodes(node: &BvhNode) -> usize {
    1 + node.left.as_ref().map_or(0, |n| count_nodes(n))
        + node.right.as_ref().map_or(0, |n| count_nodes(n))
}

fn node_depth(node: &BvhNode) -> usize {
    1 + node
        .left
        .as_ref()
        .map_or(0, |n| node_depth(n))
        .max(node.right.as_ref().map_or(0, |n| node_depth(n)))
}

// ============================================================================
// Flat BVH
// ============================================================================

/// Flatten a [`Bvh`] into a `Vec<FlatBvhNode>` (DFS pre-order) together with
/// a reordered primitive-index slice.
///
/// Returns `(flat_nodes, prim_indices)` where `prim_indices[i]` is an index
/// into `bvh.primitives`.
pub fn flatten(bvh: &Bvh) -> (Vec<FlatBvhNode>, Vec<usize>) {
    let mut nodes: Vec<FlatBvhNode> = Vec::new();
    let mut prim_indices: Vec<usize> = Vec::new();

    if let Some(root) = &bvh.root {
        flatten_recursive(root, &mut nodes, &mut prim_indices);
    }

    (nodes, prim_indices)
}

/// Returns the index at which `node` was stored.
fn flatten_recursive(
    node: &BvhNode,
    nodes: &mut Vec<FlatBvhNode>,
    prim_indices: &mut Vec<usize>,
) -> usize {
    let node_idx = nodes.len();

    if node.is_leaf() {
        let first = prim_indices.len() as u32;
        let count = node.primitives.len() as u32;
        prim_indices.extend_from_slice(&node.primitives);
        nodes.push(FlatBvhNode {
            aabb: node.aabb.clone(),
            left_first: first,
            count,
        });
    } else {
        // Reserve a slot; left_first (right child index) filled after recursion.
        nodes.push(FlatBvhNode {
            aabb: node.aabb.clone(),
            left_first: 0,
            count: 0,
        });
        // Left child is always node_idx + 1 (no explicit storage needed).
        if let Some(left) = &node.left {
            flatten_recursive(left, nodes, prim_indices);
        }
        // Right child comes after the entire left subtree.
        let right_idx = if let Some(right) = &node.right {
            flatten_recursive(right, nodes, prim_indices)
        } else {
            0
        };
        nodes[node_idx].left_first = right_idx as u32;
    }

    node_idx
}

/// Iterative AABB query over a flat BVH.
///
/// Returns the `object_id` values of all primitives whose AABB overlaps
/// `query`. `bvh_primitives` is the `Bvh::primitives` slice.
pub fn query_flat(
    nodes: &[FlatBvhNode],
    prim_indices: &[usize],
    bvh_primitives: &[BvhPrimitive],
    query: &Aabb,
) -> Vec<usize> {
    let mut result = Vec::new();
    if nodes.is_empty() {
        return result;
    }

    let mut stack: Vec<usize> = Vec::with_capacity(64);
    stack.push(0);

    while let Some(idx) = stack.pop() {
        let node = &nodes[idx];
        if !node.aabb.intersects(query) {
            continue;
        }
        if node.count > 0 {
            // Leaf
            let start = node.left_first as usize;
            let end = start + node.count as usize;
            for &pi in &prim_indices[start..end] {
                if bvh_primitives[pi].aabb.intersects(query) {
                    result.push(bvh_primitives[pi].object_id);
                }
            }
        } else {
            // Internal: left child is at idx+1, right child at left_first.
            let right = node.left_first as usize;
            stack.push(right);
            stack.push(idx + 1);
        }
    }

    result
}

// ============================================================================
// BVH Traversal (closest hit)
// ============================================================================

/// Traverse the BVH returning the **closest** hit (smallest positive t).
pub fn bvh_closest_hit(
    bvh: &Bvh,
    origin: [f32; 3],
    direction: [f32; 3],
    max_t: f32,
) -> Option<RayHit> {
    let inv_dir = [1.0 / direction[0], 1.0 / direction[1], 1.0 / direction[2]];
    let root = bvh.root.as_ref()?;
    let mut best: Option<RayHit> = None;
    let mut current_max = max_t;
    closest_hit_recursive(
        root,
        origin,
        inv_dir,
        &bvh.primitives,
        &mut best,
        &mut current_max,
    );
    best
}

fn closest_hit_recursive(
    node: &BvhNode,
    origin: [f32; 3],
    inv_dir: [f32; 3],
    primitives: &[BvhPrimitive],
    best: &mut Option<RayHit>,
    max_t: &mut f32,
) {
    if ray_aabb_t(origin, inv_dir, &node.aabb).is_none() {
        return;
    }
    if node.is_leaf() {
        for &idx in &node.primitives {
            if let Some((t_min, _)) = ray_aabb_t(origin, inv_dir, &primitives[idx].aabb)
                && t_min >= 0.0
                && t_min < *max_t
            {
                *max_t = t_min;
                *best = Some(RayHit {
                    object_id: primitives[idx].object_id,
                    t: t_min,
                });
            }
        }
    } else {
        if let Some(left) = &node.left {
            closest_hit_recursive(left, origin, inv_dir, primitives, best, max_t);
        }
        if let Some(right) = &node.right {
            closest_hit_recursive(right, origin, inv_dir, primitives, best, max_t);
        }
    }
}

// ============================================================================
// BVH Refit
// ============================================================================

/// Refit the bounding boxes of an existing BVH after primitives have moved.
///
/// The topology (splits) are preserved; only bounding boxes are recomputed.
pub fn refit(node: &mut BvhNode, primitives: &[BvhPrimitive]) {
    if node.is_leaf() {
        if !node.primitives.is_empty() {
            node.aabb = bounding_box(primitives, &node.primitives);
        }
        return;
    }
    if let Some(left) = node.left.as_mut() {
        refit(left, primitives);
    }
    if let Some(right) = node.right.as_mut() {
        refit(right, primitives);
    }
    // Recompute bounding box from children.
    let left_aabb = node.left.as_ref().map(|n| n.aabb.clone());
    let right_aabb = node.right.as_ref().map(|n| n.aabb.clone());
    node.aabb = match (left_aabb, right_aabb) {
        (Some(l), Some(r)) => Aabb::merge(&l, &r),
        (Some(l), None) => l,
        (None, Some(r)) => r,
        (None, None) => node.aabb.clone(),
    };
}

// ============================================================================
// Morton code (LBVH)
// ============================================================================

/// Expand a 10-bit integer into 30 bits by inserting two zeros before each bit.
fn expand_bits(mut v: u32) -> u32 {
    v = (v | (v << 16)) & 0x030000FF;
    v = (v | (v << 8)) & 0x0300F00F;
    v = (v | (v << 4)) & 0x030C30C3;
    v = (v | (v << 2)) & 0x09249249;
    v
}

/// Compute a 30-bit Morton code for a 3D point normalised to \[0, 1\]^3.
pub fn morton_code(p: [f32; 3]) -> u32 {
    let x = (p[0].clamp(0.0, 1.0) * 1023.0) as u32;
    let y = (p[1].clamp(0.0, 1.0) * 1023.0) as u32;
    let z = (p[2].clamp(0.0, 1.0) * 1023.0) as u32;
    expand_bits(x) | (expand_bits(y) << 1) | (expand_bits(z) << 2)
}

impl LbvhPrimitive {
    /// Construct an `LbvhPrimitive`, computing the Morton code from the centroid
    /// normalised by `scene_aabb`.
    pub fn new(aabb: Aabb, object_id: usize, scene_aabb: &Aabb) -> Self {
        let c = aabb.center();
        let scene_size = [
            (scene_aabb.max[0] - scene_aabb.min[0]).max(1e-10),
            (scene_aabb.max[1] - scene_aabb.min[1]).max(1e-10),
            (scene_aabb.max[2] - scene_aabb.min[2]).max(1e-10),
        ];
        let norm = [
            (c[0] - scene_aabb.min[0]) / scene_size[0],
            (c[1] - scene_aabb.min[1]) / scene_size[1],
            (c[2] - scene_aabb.min[2]) / scene_size[2],
        ];
        let morton = morton_code(norm);
        Self {
            aabb,
            object_id,
            morton,
        }
    }
}

/// Build an LBVH (Linear BVH) from a set of primitives using Morton-code
/// sorting and a recursive binary splitting strategy.
///
/// Returns a standard [`Bvh`] so the same query functions can be used.
pub fn lbvh_build(primitives: Vec<BvhPrimitive>) -> Bvh {
    if primitives.is_empty() {
        return Bvh {
            root: None,
            primitives,
        };
    }

    // Compute scene bounding box.
    let mut scene = primitives[0].aabb.clone();
    for p in &primitives[1..] {
        scene = Aabb::merge(&scene, &p.aabb);
    }

    // Assign Morton codes and sort.
    let mut indexed: Vec<(u32, usize)> = primitives
        .iter()
        .enumerate()
        .map(|(i, p)| {
            let lp = LbvhPrimitive::new(p.aabb.clone(), p.object_id, &scene);
            (lp.morton, i)
        })
        .collect();
    indexed.sort_unstable_by_key(|&(m, _)| m);

    let sorted_indices: Vec<usize> = indexed.iter().map(|&(_, i)| i).collect();
    let root = lbvh_recursive(&primitives, &sorted_indices);

    Bvh {
        root: Some(root),
        primitives,
    }
}

fn lbvh_recursive(primitives: &[BvhPrimitive], indices: &[usize]) -> BvhNode {
    let aabb = bounding_box(primitives, indices);

    if indices.len() <= LEAF_SIZE {
        return BvhNode {
            aabb,
            left: None,
            right: None,
            primitives: indices.to_vec(),
        };
    }

    let mid = indices.len() / 2;
    let left = lbvh_recursive(primitives, &indices[..mid]);
    let right = lbvh_recursive(primitives, &indices[mid..]);

    BvhNode {
        aabb,
        left: Some(Box::new(left)),
        right: Some(Box::new(right)),
        primitives: Vec::new(),
    }
}

// ============================================================================
// HLBVH Split
// ============================================================================

/// Find the split index for a slice of Morton-code-sorted primitives using
/// the highest differing bit (HLBVH strategy).
///
/// Returns the split position (0 < split < len).
pub fn hlbvh_split(mortons: &[u32]) -> usize {
    if mortons.len() < 2 {
        return 1;
    }
    let first = mortons[0];
    let last = mortons[mortons.len() - 1];
    let common_prefix = (first ^ last).leading_zeros();
    // Binary search for the split where the highest bit differs.
    let mut lo = 0usize;
    let mut hi = mortons.len() - 1;
    while hi - lo > 1 {
        let mid = (lo + hi) / 2;
        let prefix = (first ^ mortons[mid]).leading_zeros();
        if prefix > common_prefix {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    hi
}

// ============================================================================
// Morton cluster helpers
// ============================================================================

/// Build a flat BVH from a pre-sorted (by Morton code) slice of `LbvhPrimitive`s
/// using a radix-sort-inspired clustering strategy.
pub fn compute_bvh_from_sorted(sorted: &[LbvhPrimitive]) -> Bvh {
    if sorted.is_empty() {
        return Bvh {
            root: None,
            primitives: Vec::new(),
        };
    }

    // Reconstruct BvhPrimitives in sorted order.
    let primitives: Vec<BvhPrimitive> = sorted
        .iter()
        .map(|lp| BvhPrimitive::new(lp.aabb.clone(), lp.object_id))
        .collect();

    let mortons: Vec<u32> = sorted.iter().map(|lp| lp.morton).collect();
    let indices: Vec<usize> = (0..primitives.len()).collect();
    let root = bvh_from_sorted_recursive(&primitives, &indices, &mortons);
    Bvh {
        root: Some(root),
        primitives,
    }
}

fn bvh_from_sorted_recursive(
    primitives: &[BvhPrimitive],
    indices: &[usize],
    mortons: &[u32],
) -> BvhNode {
    let aabb = bounding_box(primitives, indices);
    if indices.len() <= LEAF_SIZE {
        return BvhNode {
            aabb,
            left: None,
            right: None,
            primitives: indices.to_vec(),
        };
    }
    // Use HLBVH-style split on Morton codes at corresponding positions.
    let local_mortons: Vec<u32> = indices.iter().map(|&i| mortons[i]).collect();
    let split = hlbvh_split(&local_mortons);
    let left = bvh_from_sorted_recursive(primitives, &indices[..split], mortons);
    let right = bvh_from_sorted_recursive(primitives, &indices[split..], mortons);
    BvhNode {
        aabb,
        left: Some(Box::new(left)),
        right: Some(Box::new(right)),
        primitives: Vec::new(),
    }
}

/// Compute the bounding sphere radius for a cluster of `LbvhPrimitive`s.
pub fn compute_cluster_radius(cluster: &[LbvhPrimitive]) -> f32 {
    if cluster.is_empty() {
        return 0.0;
    }
    // Compute merged AABB.
    let mut merged = cluster[0].aabb.clone();
    for lp in &cluster[1..] {
        merged = Aabb::merge(&merged, &lp.aabb);
    }
    let cx = (merged.min[0] + merged.max[0]) * 0.5;
    let cy = (merged.min[1] + merged.max[1]) * 0.5;
    let cz = (merged.min[2] + merged.max[2]) * 0.5;

    let mut max_dist_sq = 0.0_f32;
    for lp in cluster {
        let c = lp.aabb.center();
        let dx = c[0] - cx;
        let dy = c[1] - cy;
        let dz = c[2] - cz;
        let d2 = dx * dx + dy * dy + dz * dz;
        if d2 > max_dist_sq {
            max_dist_sq = d2;
        }
    }
    max_dist_sq.sqrt()
}

/// Compute clusters by grouping Morton-sorted `LbvhPrimitive`s into chunks of
/// `cluster_size` and returning a `MortonCluster` per group.
pub fn build_morton_clusters(sorted: &[LbvhPrimitive], cluster_size: usize) -> Vec<MortonCluster> {
    if sorted.is_empty() || cluster_size == 0 {
        return Vec::new();
    }
    sorted
        .chunks(cluster_size)
        .map(|chunk| {
            let indices: Vec<usize> = (0..chunk.len()).collect();
            let mut aabb = chunk[0].aabb.clone();
            for lp in &chunk[1..] {
                aabb = Aabb::merge(&aabb, &lp.aabb);
            }
            let radius = compute_cluster_radius(chunk);
            MortonCluster {
                indices,
                aabb,
                radius,
            }
        })
        .collect()
}

// ============================================================================
// BVH Statistics
// ============================================================================

impl BvhStats {
    /// Compute statistics by traversing the given BVH.
    pub fn compute(bvh: &Bvh) -> Self {
        let mut s = BvhStats {
            node_count: 0,
            leaf_count: 0,
            internal_count: 0,
            max_depth: 0,
            total_primitives: 0,
            avg_primitives_per_leaf: 0.0,
        };
        if let Some(root) = &bvh.root {
            collect_stats(root, 1, &mut s);
        }
        if s.leaf_count > 0 {
            s.avg_primitives_per_leaf = s.total_primitives as f32 / s.leaf_count as f32;
        }
        s
    }
}

fn collect_stats(node: &BvhNode, depth: usize, s: &mut BvhStats) {
    s.node_count += 1;
    if depth > s.max_depth {
        s.max_depth = depth;
    }
    if node.is_leaf() {
        s.leaf_count += 1;
        s.total_primitives += node.primitives.len();
    } else {
        s.internal_count += 1;
        if let Some(left) = &node.left {
            collect_stats(left, depth + 1, s);
        }
        if let Some(right) = &node.right {
            collect_stats(right, depth + 1, s);
        }
    }
}

impl BvhTreeStatistics {
    /// Compute extended tree statistics by traversing the given BVH.
    pub fn compute(bvh: &Bvh) -> Self {
        let mut s = BvhTreeStatistics {
            node_count: 0,
            leaf_count: 0,
            internal_count: 0,
            max_depth: 0,
            total_primitives: 0,
            avg_fanout: 0.0,
            total_leaf_surface_area: 0.0,
        };
        if let Some(root) = &bvh.root {
            let mut child_sum = 0usize;
            collect_tree_stats(root, 1, &mut s, &mut child_sum);
            s.avg_fanout = if s.internal_count > 0 {
                child_sum as f32 / s.internal_count as f32
            } else {
                0.0
            };
        }
        s
    }
}

fn collect_tree_stats(
    node: &BvhNode,
    depth: usize,
    s: &mut BvhTreeStatistics,
    child_sum: &mut usize,
) {
    s.node_count += 1;
    if depth > s.max_depth {
        s.max_depth = depth;
    }
    if node.is_leaf() {
        s.leaf_count += 1;
        s.total_primitives += node.primitives.len();
        s.total_leaf_surface_area += node.aabb.surface_area();
    } else {
        s.internal_count += 1;
        let mut children = 0usize;
        if let Some(left) = &node.left {
            children += 1;
            collect_tree_stats(left, depth + 1, s, child_sum);
        }
        if let Some(right) = &node.right {
            children += 1;
            collect_tree_stats(right, depth + 1, s, child_sum);
        }
        *child_sum += children;
    }
}
