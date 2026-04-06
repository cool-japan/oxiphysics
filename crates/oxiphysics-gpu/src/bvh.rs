// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! CPU BVH (Bounding Volume Hierarchy) tree for broad-phase acceleration.
//!
//! Provides axis-aligned bounding box (AABB) structures, a recursive BVH tree
//! built with a SAH-inspired median split, ray and AABB queries, and a
//! linearised flat representation for cache-friendly traversal.

// --------------------------------------------------------------------------
// Silence dead-code warnings for public API that may not be exercised at link
// time within this crate.
// --------------------------------------------------------------------------
#![allow(dead_code)]

// ============================================================================
// Aabb
// ============================================================================

/// An axis-aligned bounding box stored as two `[f32; 3]` corners.
#[derive(Debug, Clone, PartialEq)]
pub struct Aabb {
    /// Minimum corner (component-wise).
    pub min: [f32; 3],
    /// Maximum corner (component-wise).
    pub max: [f32; 3],
}

impl Aabb {
    /// Construct an `Aabb` from explicit min/max corners.
    pub fn new(min: [f32; 3], max: [f32; 3]) -> Self {
        Self { min, max }
    }

    /// Construct a degenerate `Aabb` that contains exactly one point.
    pub fn point(p: [f32; 3]) -> Self {
        Self { min: p, max: p }
    }

    /// Return the smallest `Aabb` that contains both `a` and `b`.
    pub fn merge(a: &Aabb, b: &Aabb) -> Aabb {
        Aabb {
            min: [
                a.min[0].min(b.min[0]),
                a.min[1].min(b.min[1]),
                a.min[2].min(b.min[2]),
            ],
            max: [
                a.max[0].max(b.max[0]),
                a.max[1].max(b.max[1]),
                a.max[2].max(b.max[2]),
            ],
        }
    }

    /// Returns `true` if this box overlaps `other` (touching counts).
    pub fn intersects(&self, other: &Aabb) -> bool {
        self.min[0] <= other.max[0]
            && self.max[0] >= other.min[0]
            && self.min[1] <= other.max[1]
            && self.max[1] >= other.min[1]
            && self.min[2] <= other.max[2]
            && self.max[2] >= other.min[2]
    }

    /// Returns `true` if point `p` lies inside or on the surface of this box.
    pub fn contains(&self, p: [f32; 3]) -> bool {
        p[0] >= self.min[0]
            && p[0] <= self.max[0]
            && p[1] >= self.min[1]
            && p[1] <= self.max[1]
            && p[2] >= self.min[2]
            && p[2] <= self.max[2]
    }

    /// Surface area of the box (sum of all six face areas).
    pub fn surface_area(&self) -> f32 {
        let dx = self.max[0] - self.min[0];
        let dy = self.max[1] - self.min[1];
        let dz = self.max[2] - self.min[2];
        2.0 * (dx * dy + dy * dz + dz * dx)
    }

    /// Geometric centre of the box.
    pub fn center(&self) -> [f32; 3] {
        [
            0.5 * (self.min[0] + self.max[0]),
            0.5 * (self.min[1] + self.max[1]),
            0.5 * (self.min[2] + self.max[2]),
        ]
    }

    /// Return a copy of this box expanded uniformly by `margin` on every side.
    pub fn expand(&self, margin: f32) -> Aabb {
        Aabb {
            min: [
                self.min[0] - margin,
                self.min[1] - margin,
                self.min[2] - margin,
            ],
            max: [
                self.max[0] + margin,
                self.max[1] + margin,
                self.max[2] + margin,
            ],
        }
    }
}

// ============================================================================
// BvhPrimitive
// ============================================================================

/// A leaf primitive: an AABB together with the logical object it belongs to.
#[derive(Debug, Clone)]
pub struct BvhPrimitive {
    /// Bounding box of the primitive.
    pub aabb: Aabb,
    /// Caller-defined identifier (returned by queries).
    pub object_id: usize,
}

impl BvhPrimitive {
    /// Construct a `BvhPrimitive`.
    pub fn new(aabb: Aabb, object_id: usize) -> Self {
        Self { aabb, object_id }
    }
}

// ============================================================================
// BvhNode
// ============================================================================

/// A node in the recursive BVH tree.
#[derive(Debug)]
pub struct BvhNode {
    /// Bounding box that contains all children / primitives.
    pub aabb: Aabb,
    /// Left subtree (internal nodes only).
    pub left: Option<Box<BvhNode>>,
    /// Right subtree (internal nodes only).
    pub right: Option<Box<BvhNode>>,
    /// Indices into `Bvh::primitives` (leaf nodes only).
    pub primitives: Vec<usize>,
}

impl BvhNode {
    /// Returns `true` if this node is a leaf (holds primitives directly).
    pub fn is_leaf(&self) -> bool {
        self.left.is_none() && self.right.is_none()
    }
}

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

// ============================================================================
// Bvh
// ============================================================================

/// Maximum number of primitives per leaf before the tree stops splitting.
const LEAF_SIZE: usize = 4;

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

fn bounding_box(primitives: &[BvhPrimitive], indices: &[usize]) -> Aabb {
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

/// A single node in the linearised (flat) BVH representation.
///
/// Layout:
/// * If `count == 0` this is an **internal** node.
///   - The **left** child is always at index `node_idx + 1` (i.e. stored
///     immediately after the parent in DFS pre-order).
///   - `left_first` holds the index of the **right** child.
/// * If `count > 0` this is a **leaf**; `left_first` is the start index into
///   the accompanying primitive-index slice and `count` is the number of
///   entries.
#[derive(Debug, Clone)]
pub struct FlatBvhNode {
    /// Bounding box of this node.
    pub aabb: Aabb,
    /// Right-child index (internal) or first-primitive index (leaf).
    pub left_first: u32,
    /// 0 for internal nodes; number of primitives for leaf nodes.
    pub count: u32,
}

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

/// LBVH primitive: AABB + Morton code.
#[derive(Debug, Clone)]
pub struct LbvhPrimitive {
    /// Bounding box.
    pub aabb: Aabb,
    /// Caller-defined object ID.
    pub object_id: usize,
    /// 30-bit Morton code computed from the AABB centroid.
    pub morton: u32,
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
// BVH Traversal (closest hit)
// ============================================================================

/// Result of a closest-hit ray traversal.
#[derive(Debug, Clone)]
pub struct RayHit {
    /// Object ID of the closest hit primitive.
    pub object_id: usize,
    /// Ray parameter at which the hit occurred.
    pub t: f32,
}

/// Ray–AABB slab intersection returning the near/far t values.
fn ray_aabb_t(origin: [f32; 3], inv_dir: [f32; 3], aabb: &Aabb) -> Option<(f32, f32)> {
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
// HLBVH Split (spatial median on the highest non-degenerate bit)
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
// BVH Statistics
// ============================================================================

/// Runtime statistics about a BVH tree.
#[derive(Debug, Clone)]
pub struct BvhStats {
    /// Total number of nodes (internal + leaf).
    pub node_count: usize,
    /// Number of leaf nodes.
    pub leaf_count: usize,
    /// Number of internal nodes.
    pub internal_count: usize,
    /// Maximum tree depth.
    pub max_depth: usize,
    /// Total number of primitives stored across all leaves.
    pub total_primitives: usize,
    /// Average primitives per leaf.
    pub avg_primitives_per_leaf: f32,
}

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

// ============================================================================
// MortonCluster — radix-sorted BVH construction helpers
// ============================================================================

/// A cluster of Morton-coded primitives, with a pre-computed bounding radius.
#[derive(Debug, Clone)]
pub struct MortonCluster {
    /// Indices of the primitives in this cluster (into the parent slice).
    pub indices: Vec<usize>,
    /// Axis-aligned bounding box of the cluster.
    pub aabb: Aabb,
    /// Bounding sphere radius (centred at the AABB centre).
    pub radius: f32,
}

/// Build a flat BVH from a pre-sorted (by Morton code) slice of `LbvhPrimitive`s
/// using a radix-sort-inspired clustering strategy.
///
/// Primitives are split at the highest differing Morton bit, producing a
/// balanced binary tree stored as a flat `Bvh`.
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
///
/// The sphere is centred at the centroid of the cluster AABB and has the
/// minimum radius that encloses all primitive centroids.
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
// BvhNode tree statistics (extended)
// ============================================================================

/// Extended BVH tree statistics including average fan-out.
#[derive(Debug, Clone)]
pub struct BvhTreeStatistics {
    /// Total node count.
    pub node_count: usize,
    /// Number of leaf nodes.
    pub leaf_count: usize,
    /// Number of internal nodes.
    pub internal_count: usize,
    /// Maximum depth from root (1-indexed).
    pub max_depth: usize,
    /// Total number of primitives across all leaves.
    pub total_primitives: usize,
    /// Average number of children per internal node (fan-out).
    /// For a binary tree this is at most 2.
    pub avg_fanout: f32,
    /// Total surface area of all leaf AABBs.
    pub total_leaf_surface_area: f32,
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

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    // ------------------------------------------------------------------
    // Aabb tests
    // ------------------------------------------------------------------

    #[test]
    fn aabb_new_stores_corners() {
        let a = Aabb::new([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]);
        assert_eq!(a.min, [1.0, 2.0, 3.0]);
        assert_eq!(a.max, [4.0, 5.0, 6.0]);
    }

    #[test]
    fn aabb_point_is_degenerate() {
        let p = [3.0, 3.0, 3.0];
        let a = Aabb::point(p);
        assert_eq!(a.min, p);
        assert_eq!(a.max, p);
    }

    #[test]
    fn aabb_merge_covers_both() {
        let a = Aabb::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = Aabb::new([2.0, 2.0, 2.0], [3.0, 3.0, 3.0]);
        let m = Aabb::merge(&a, &b);
        assert_eq!(m.min, [0.0, 0.0, 0.0]);
        assert_eq!(m.max, [3.0, 3.0, 3.0]);
    }

    #[test]
    fn aabb_intersects_overlapping() {
        let a = Aabb::new([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
        let b = Aabb::new([1.0, 1.0, 1.0], [3.0, 3.0, 3.0]);
        assert!(a.intersects(&b));
    }

    #[test]
    fn aabb_intersects_disjoint() {
        let a = Aabb::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = Aabb::new([2.0, 2.0, 2.0], [3.0, 3.0, 3.0]);
        assert!(!a.intersects(&b));
    }

    #[test]
    fn aabb_intersects_touching_edge() {
        let a = Aabb::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = Aabb::new([1.0, 0.0, 0.0], [2.0, 1.0, 1.0]);
        assert!(a.intersects(&b));
    }

    #[test]
    fn aabb_contains_inside() {
        let a = Aabb::new([0.0, 0.0, 0.0], [4.0, 4.0, 4.0]);
        assert!(a.contains([2.0, 2.0, 2.0]));
    }

    #[test]
    fn aabb_contains_outside() {
        let a = Aabb::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert!(!a.contains([2.0, 0.0, 0.0]));
    }

    #[test]
    fn aabb_contains_on_surface() {
        let a = Aabb::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert!(a.contains([1.0, 0.5, 0.5]));
    }

    #[test]
    fn aabb_surface_area_unit_cube() {
        let a = Aabb::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert!((a.surface_area() - 6.0).abs() < 1e-6);
    }

    #[test]
    fn aabb_surface_area_flat() {
        // 2×3×0 slab: area = 2*(2*3 + 3*0 + 0*2) = 12
        let a = Aabb::new([0.0, 0.0, 0.0], [2.0, 3.0, 0.0]);
        assert!((a.surface_area() - 12.0).abs() < 1e-6);
    }

    #[test]
    fn aabb_center_correct() {
        let a = Aabb::new([0.0, 0.0, 0.0], [2.0, 4.0, 6.0]);
        let c = a.center();
        assert!((c[0] - 1.0).abs() < 1e-6);
        assert!((c[1] - 2.0).abs() < 1e-6);
        assert!((c[2] - 3.0).abs() < 1e-6);
    }

    #[test]
    fn aabb_expand_increases_bounds() {
        let a = Aabb::new([1.0, 1.0, 1.0], [2.0, 2.0, 2.0]);
        let e = a.expand(0.5);
        assert_eq!(e.min, [0.5, 0.5, 0.5]);
        assert_eq!(e.max, [2.5, 2.5, 2.5]);
    }

    // ------------------------------------------------------------------
    // SAH cost
    // ------------------------------------------------------------------

    #[test]
    fn sah_cost_balanced() {
        // Both halves equal area and count → cost == n_left + n_right
        let cost = sah_cost(4, 1.0, 4, 1.0, 2.0);
        // (1/2)*4 + (1/2)*4 = 4
        assert!((cost - 4.0).abs() < 1e-6);
    }

    #[test]
    fn sah_cost_zero_parent_area_returns_max() {
        let cost = sah_cost(1, 1.0, 1, 1.0, 0.0);
        assert_eq!(cost, f32::MAX);
    }

    // ------------------------------------------------------------------
    // Ray–AABB slab intersection
    // ------------------------------------------------------------------

    #[test]
    fn ray_hits_unit_cube() {
        let aabb = Aabb::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let origin = [-1.0, 0.5, 0.5];
        let dir = [1.0, 0.0, 0.0];
        let inv = [1.0 / dir[0], 1.0 / dir[1], 1.0 / dir[2]];
        assert!(ray_aabb_intersect(origin, inv, &aabb, 10.0));
    }

    #[test]
    fn ray_misses_unit_cube() {
        let aabb = Aabb::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let origin = [-1.0, 2.0, 0.5];
        let dir = [1.0, 0.0, 0.0];
        let inv = [1.0 / dir[0], 1.0 / dir[1], 1.0 / dir[2]];
        assert!(!ray_aabb_intersect(origin, inv, &aabb, 10.0));
    }

    #[test]
    fn ray_too_short_misses() {
        let aabb = Aabb::new([5.0, 0.0, 0.0], [6.0, 1.0, 1.0]);
        let origin = [0.0, 0.5, 0.5];
        let dir = [1.0, 0.0, 0.0];
        let inv = [1.0 / dir[0], 1.0 / dir[1], 1.0 / dir[2]];
        assert!(!ray_aabb_intersect(origin, inv, &aabb, 3.0));
    }

    // ------------------------------------------------------------------
    // Bvh build / query
    // ------------------------------------------------------------------

    fn make_grid_primitives(n: usize) -> Vec<BvhPrimitive> {
        (0..n)
            .map(|i| {
                let x = i as f32;
                BvhPrimitive::new(Aabb::new([x, 0.0, 0.0], [x + 1.0, 1.0, 1.0]), i)
            })
            .collect()
    }

    #[test]
    fn bvh_build_empty() {
        let bvh = Bvh::build(vec![]);
        assert!(bvh.root.is_none());
        assert_eq!(bvh.node_count(), 0);
        assert_eq!(bvh.depth(), 0);
    }

    #[test]
    fn bvh_build_single() {
        let prims = make_grid_primitives(1);
        let bvh = Bvh::build(prims);
        assert!(bvh.root.is_some());
        assert!(bvh.root.as_ref().unwrap().is_leaf());
        assert_eq!(bvh.node_count(), 1);
        assert_eq!(bvh.depth(), 1);
    }

    #[test]
    fn bvh_query_aabb_finds_overlap() {
        let prims = make_grid_primitives(10);
        let bvh = Bvh::build(prims);
        // Query the box that overlaps object 5 only.
        let query = Aabb::new([5.1, 0.1, 0.1], [5.9, 0.9, 0.9]);
        let mut hits = bvh.query_aabb(&query);
        hits.sort();
        assert_eq!(hits, vec![5]);
    }

    #[test]
    fn bvh_query_aabb_empty_result() {
        let prims = make_grid_primitives(5);
        let bvh = Bvh::build(prims);
        let query = Aabb::new([100.0, 0.0, 0.0], [101.0, 1.0, 1.0]);
        assert!(bvh.query_aabb(&query).is_empty());
    }

    #[test]
    fn bvh_query_aabb_finds_multiple() {
        let prims = make_grid_primitives(10);
        let bvh = Bvh::build(prims);
        // Query spanning objects 2, 3, 4
        let query = Aabb::new([2.1, 0.1, 0.1], [4.9, 0.9, 0.9]);
        let mut hits = bvh.query_aabb(&query);
        hits.sort();
        assert_eq!(hits, vec![2, 3, 4]);
    }

    #[test]
    fn bvh_query_ray_hits() {
        let prims = make_grid_primitives(8);
        let bvh = Bvh::build(prims);
        // Ray along X axis at y=0.5, z=0.5 hits all 8 primitives.
        let mut hits = bvh.query_ray([-1.0, 0.5, 0.5], [1.0, 0.0, 0.0], 20.0);
        hits.sort();
        assert_eq!(hits, (0..8).collect::<Vec<_>>());
    }

    #[test]
    fn bvh_query_ray_misses() {
        let prims = make_grid_primitives(5);
        let bvh = Bvh::build(prims);
        let hits = bvh.query_ray([0.5, 10.0, 0.5], [0.0, 1.0, 0.0], 100.0);
        assert!(hits.is_empty());
    }

    #[test]
    fn bvh_node_count_and_depth_consistent() {
        let prims = make_grid_primitives(16);
        let bvh = Bvh::build(prims);
        // With LEAF_SIZE=4 and 16 prims, depth should be at least 2.
        assert!(bvh.depth() >= 2);
        // An n-primitive tree has at most 2n-1 nodes.
        assert!(bvh.node_count() < 2 * 16);
    }

    // ------------------------------------------------------------------
    // Flat BVH
    // ------------------------------------------------------------------

    #[test]
    fn flatten_empty_bvh() {
        let bvh = Bvh::build(vec![]);
        let (nodes, prim_indices) = flatten(&bvh);
        assert!(nodes.is_empty());
        assert!(prim_indices.is_empty());
    }

    #[test]
    fn flatten_single_primitive() {
        let prims = make_grid_primitives(1);
        let bvh = Bvh::build(prims);
        let (nodes, prim_indices) = flatten(&bvh);
        assert_eq!(nodes.len(), 1);
        assert_eq!(prim_indices.len(), 1);
        assert_eq!(nodes[0].count, 1);
    }

    #[test]
    fn query_flat_finds_overlap() {
        let prims = make_grid_primitives(10);
        let bvh = Bvh::build(prims);
        let (nodes, prim_indices) = flatten(&bvh);
        let query = Aabb::new([3.1, 0.1, 0.1], [3.9, 0.9, 0.9]);
        let mut hits = query_flat(&nodes, &prim_indices, &bvh.primitives, &query);
        hits.sort();
        assert_eq!(hits, vec![3]);
    }

    #[test]
    fn query_flat_empty_result() {
        let prims = make_grid_primitives(5);
        let bvh = Bvh::build(prims);
        let (nodes, prim_indices) = flatten(&bvh);
        let query = Aabb::new([50.0, 0.0, 0.0], [51.0, 1.0, 1.0]);
        assert!(query_flat(&nodes, &prim_indices, &bvh.primitives, &query).is_empty());
    }

    #[test]
    fn query_flat_matches_recursive() {
        let prims = make_grid_primitives(20);
        let bvh = Bvh::build(prims);
        let query = Aabb::new([7.1, 0.0, 0.0], [12.9, 1.0, 1.0]);
        let mut recursive_hits = bvh.query_aabb(&query);
        recursive_hits.sort();

        let (nodes, prim_indices) = flatten(&bvh);
        let mut flat_hits = query_flat(&nodes, &prim_indices, &bvh.primitives, &query);
        flat_hits.sort();

        assert_eq!(recursive_hits, flat_hits);
    }

    // ------------------------------------------------------------------
    // Morton code tests
    // ------------------------------------------------------------------

    #[test]
    fn morton_origin_is_zero() {
        assert_eq!(morton_code([0.0, 0.0, 0.0]), 0);
    }

    #[test]
    fn morton_increases_along_x() {
        let m0 = morton_code([0.0, 0.0, 0.0]);
        let m1 = morton_code([0.5, 0.0, 0.0]);
        let m2 = morton_code([1.0, 0.0, 0.0]);
        // With all y,z=0 the Morton code grows with x.
        // (bits interleaved so we check non-decreasing)
        assert!(m0 <= m1, "m0={} m1={}", m0, m1);
        assert!(m1 <= m2, "m1={} m2={}", m1, m2);
    }

    #[test]
    fn morton_clamps_outside_unit_cube() {
        let m_neg = morton_code([-1.0, -1.0, -1.0]);
        let m_zero = morton_code([0.0, 0.0, 0.0]);
        assert_eq!(m_neg, m_zero);

        let m_big = morton_code([2.0, 2.0, 2.0]);
        let m_one = morton_code([1.0, 1.0, 1.0]);
        assert_eq!(m_big, m_one);
    }

    // ------------------------------------------------------------------
    // LBVH construction tests
    // ------------------------------------------------------------------

    #[test]
    fn lbvh_build_empty() {
        let bvh = lbvh_build(vec![]);
        assert!(bvh.root.is_none());
    }

    #[test]
    fn lbvh_build_single() {
        let prims = make_grid_primitives(1);
        let bvh = lbvh_build(prims);
        assert!(bvh.root.is_some());
        assert!(bvh.root.as_ref().unwrap().is_leaf());
    }

    #[test]
    fn lbvh_build_query_finds_correct_objects() {
        let prims = make_grid_primitives(10);
        let bvh = lbvh_build(prims);
        let query = Aabb::new([4.1, 0.1, 0.1], [4.9, 0.9, 0.9]);
        let mut hits = bvh.query_aabb(&query);
        hits.sort();
        assert_eq!(hits, vec![4]);
    }

    #[test]
    fn lbvh_build_covers_all_primitives() {
        let prims = make_grid_primitives(8);
        let bvh = lbvh_build(prims);
        // Root AABB should contain all primitives.
        let root = bvh.root.as_ref().unwrap();
        assert!(root.aabb.min[0] <= 0.0);
        assert!(root.aabb.max[0] >= 8.0);
    }

    // ------------------------------------------------------------------
    // BVH closest-hit traversal
    // ------------------------------------------------------------------

    #[test]
    fn closest_hit_returns_nearest() {
        let prims = make_grid_primitives(10);
        let bvh = Bvh::build(prims);
        // Ray along X from x=-1: should hit object 0 first (x ∈ [0,1])
        let hit = bvh_closest_hit(&bvh, [-1.0, 0.5, 0.5], [1.0, 0.0, 0.0], 100.0);
        assert!(hit.is_some(), "ray should hit something");
        let hit = hit.unwrap();
        assert_eq!(
            hit.object_id, 0,
            "closest hit should be object 0, got {}",
            hit.object_id
        );
    }

    #[test]
    fn closest_hit_misses_returns_none() {
        let prims = make_grid_primitives(5);
        let bvh = Bvh::build(prims);
        let hit = bvh_closest_hit(&bvh, [0.5, 10.0, 0.5], [0.0, 1.0, 0.0], 100.0);
        assert!(hit.is_none());
    }

    #[test]
    fn closest_hit_empty_bvh_returns_none() {
        let bvh = Bvh::build(vec![]);
        let hit = bvh_closest_hit(&bvh, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 100.0);
        assert!(hit.is_none());
    }

    #[test]
    fn closest_hit_t_is_positive() {
        let prims = make_grid_primitives(5);
        let bvh = Bvh::build(prims);
        let hit = bvh_closest_hit(&bvh, [-1.0, 0.5, 0.5], [1.0, 0.0, 0.0], 100.0);
        if let Some(h) = hit {
            assert!(h.t >= 0.0, "t should be non-negative, got {}", h.t);
        }
    }

    // ------------------------------------------------------------------
    // BVH refit tests
    // ------------------------------------------------------------------

    #[test]
    fn refit_preserves_topology() {
        let prims = make_grid_primitives(8);
        let mut bvh = Bvh::build(prims);
        let before_count = bvh.node_count();
        if let Some(root) = bvh.root.as_mut() {
            refit(root, &bvh.primitives);
        }
        assert_eq!(
            bvh.node_count(),
            before_count,
            "refit should not change node count"
        );
    }

    #[test]
    fn refit_root_aabb_covers_all() {
        let prims = make_grid_primitives(8);
        let mut bvh = Bvh::build(prims);
        if let Some(root) = bvh.root.as_mut() {
            refit(root, &bvh.primitives);
        }
        let root = bvh.root.as_ref().unwrap();
        assert!(root.aabb.min[0] <= 0.0 + 1e-5);
        assert!(root.aabb.max[0] >= 8.0 - 1e-5);
    }

    // ------------------------------------------------------------------
    // HLBVH split tests
    // ------------------------------------------------------------------

    #[test]
    fn hlbvh_split_two_distinct_values() {
        let mortons = vec![0u32, 1u32];
        let split = hlbvh_split(&mortons);
        assert_eq!(split, 1);
    }

    #[test]
    fn hlbvh_split_returns_valid_index() {
        let mortons: Vec<u32> = (0..16).map(|i| i * 64).collect();
        let split = hlbvh_split(&mortons);
        assert!(split > 0 && split < mortons.len(), "split={}", split);
    }

    #[test]
    fn hlbvh_split_equal_values_returns_one() {
        let mortons = vec![5u32; 8];
        let split = hlbvh_split(&mortons);
        // All equal → leading_zeros of 0 = 32; binary search result = 1
        assert!(split >= 1 && split < mortons.len());
    }

    // ------------------------------------------------------------------
    // BVH statistics tests
    // ------------------------------------------------------------------

    #[test]
    fn bvh_stats_empty() {
        let bvh = Bvh::build(vec![]);
        let s = BvhStats::compute(&bvh);
        assert_eq!(s.node_count, 0);
        assert_eq!(s.leaf_count, 0);
        assert_eq!(s.total_primitives, 0);
    }

    #[test]
    fn bvh_stats_single_primitive() {
        let prims = make_grid_primitives(1);
        let bvh = Bvh::build(prims);
        let s = BvhStats::compute(&bvh);
        assert_eq!(s.node_count, 1);
        assert_eq!(s.leaf_count, 1);
        assert_eq!(s.total_primitives, 1);
        assert_eq!(s.max_depth, 1);
    }

    #[test]
    fn bvh_stats_node_count_consistent() {
        let prims = make_grid_primitives(16);
        let bvh = Bvh::build(prims.clone());
        let s = BvhStats::compute(&bvh);
        assert_eq!(s.node_count, bvh.node_count());
        assert_eq!(s.leaf_count + s.internal_count, s.node_count);
        assert_eq!(s.total_primitives, prims.len());
    }

    #[test]
    fn bvh_stats_avg_primitives_per_leaf() {
        let prims = make_grid_primitives(8);
        let bvh = Bvh::build(prims);
        let s = BvhStats::compute(&bvh);
        assert!(s.avg_primitives_per_leaf > 0.0);
        // avg should not exceed LEAF_SIZE + 1
        assert!(s.avg_primitives_per_leaf <= (LEAF_SIZE + 1) as f32);
    }

    #[test]
    fn bvh_stats_max_depth_reasonable() {
        let prims = make_grid_primitives(32);
        let bvh = Bvh::build(prims);
        let s = BvhStats::compute(&bvh);
        // Depth should be at most log2(32/LEAF_SIZE) + 1 = ~4
        assert!(
            s.max_depth >= 1 && s.max_depth <= 20,
            "depth={}",
            s.max_depth
        );
    }

    // ------------------------------------------------------------------
    // LbvhPrimitive Morton code assignment
    // ------------------------------------------------------------------

    #[test]
    fn lbvh_primitive_morton_in_range() {
        let aabb = Aabb::new([0.5, 0.5, 0.5], [1.0, 1.0, 1.0]);
        let scene = Aabb::new([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
        let lp = LbvhPrimitive::new(aabb, 0, &scene);
        // Morton code is 30-bit so max is (1<<30)-1
        assert!(lp.morton < (1u32 << 30));
    }

    #[test]
    fn lbvh_primitive_at_origin_small_code() {
        let aabb = Aabb::point([0.0, 0.0, 0.0]);
        let scene = Aabb::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let lp = LbvhPrimitive::new(aabb, 0, &scene);
        assert_eq!(lp.morton, 0);
    }

    // ------------------------------------------------------------------
    // compute_bvh_from_sorted tests
    // ------------------------------------------------------------------

    fn make_sorted_lbvh_prims(n: usize) -> Vec<LbvhPrimitive> {
        let scene = Aabb::new([0.0, 0.0, 0.0], [n as f32 + 1.0, 1.0, 1.0]);
        let mut prims: Vec<LbvhPrimitive> = (0..n)
            .map(|i| {
                let x = i as f32;
                LbvhPrimitive::new(Aabb::new([x, 0.0, 0.0], [x + 1.0, 1.0, 1.0]), i, &scene)
            })
            .collect();
        prims.sort_unstable_by_key(|lp| lp.morton);
        prims
    }

    #[test]
    fn compute_bvh_from_sorted_empty() {
        let bvh = compute_bvh_from_sorted(&[]);
        assert!(bvh.root.is_none());
        assert_eq!(bvh.primitives.len(), 0);
    }

    #[test]
    fn compute_bvh_from_sorted_single() {
        let scene = Aabb::new([0.0, 0.0, 0.0], [2.0, 1.0, 1.0]);
        let lp = LbvhPrimitive::new(Aabb::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]), 7, &scene);
        let bvh = compute_bvh_from_sorted(&[lp]);
        assert!(bvh.root.is_some());
        assert_eq!(bvh.primitives.len(), 1);
    }

    #[test]
    fn compute_bvh_from_sorted_preserves_count() {
        let sorted = make_sorted_lbvh_prims(16);
        let bvh = compute_bvh_from_sorted(&sorted);
        assert_eq!(bvh.primitives.len(), 16);
    }

    #[test]
    fn compute_bvh_from_sorted_root_covers_all() {
        let sorted = make_sorted_lbvh_prims(8);
        let bvh = compute_bvh_from_sorted(&sorted);
        let root_aabb = &bvh.root.as_ref().unwrap().aabb;
        assert!(root_aabb.min[0] <= 0.0 + 1e-5);
        assert!(root_aabb.max[0] >= 8.0 - 1e-5);
    }

    // ------------------------------------------------------------------
    // compute_cluster_radius tests
    // ------------------------------------------------------------------

    #[test]
    fn compute_cluster_radius_empty() {
        let r = compute_cluster_radius(&[]);
        assert_eq!(r, 0.0);
    }

    #[test]
    fn compute_cluster_radius_single() {
        let scene = Aabb::new([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
        let lp = LbvhPrimitive::new(Aabb::point([1.0, 1.0, 1.0]), 0, &scene);
        let r = compute_cluster_radius(&[lp]);
        // Single point cluster: radius = 0
        assert!(
            r < 1e-6,
            "single-point cluster radius should be ~0, got {r}"
        );
    }

    #[test]
    fn compute_cluster_radius_two_points() {
        let scene = Aabb::new([0.0, 0.0, 0.0], [4.0, 1.0, 1.0]);
        let lp0 = LbvhPrimitive::new(Aabb::point([0.0, 0.0, 0.0]), 0, &scene);
        let lp1 = LbvhPrimitive::new(Aabb::point([2.0, 0.0, 0.0]), 1, &scene);
        let r = compute_cluster_radius(&[lp0, lp1]);
        // Centroid is (1,0,0); each point is at distance 1.
        assert!((r - 1.0).abs() < 1e-5, "radius should be 1.0, got {r}");
    }

    #[test]
    fn compute_cluster_radius_is_non_negative() {
        let sorted = make_sorted_lbvh_prims(12);
        let r = compute_cluster_radius(&sorted);
        assert!(r >= 0.0, "radius must be non-negative, got {r}");
    }

    // ------------------------------------------------------------------
    // BvhTreeStatistics tests
    // ------------------------------------------------------------------

    #[test]
    fn bvh_tree_stats_empty() {
        let bvh = Bvh::build(vec![]);
        let s = BvhTreeStatistics::compute(&bvh);
        assert_eq!(s.node_count, 0);
        assert_eq!(s.leaf_count, 0);
        assert_eq!(s.internal_count, 0);
        assert_eq!(s.total_primitives, 0);
    }

    #[test]
    fn bvh_tree_stats_fanout_binary() {
        let prims = make_grid_primitives(16);
        let bvh = Bvh::build(prims);
        let s = BvhTreeStatistics::compute(&bvh);
        // Binary tree: avg fanout must be <= 2.0
        assert!(s.avg_fanout <= 2.0 + 1e-6, "fanout = {}", s.avg_fanout);
    }

    #[test]
    fn bvh_tree_stats_node_count_consistent() {
        let prims = make_grid_primitives(16);
        let bvh = Bvh::build(prims.clone());
        let s = BvhTreeStatistics::compute(&bvh);
        assert_eq!(s.leaf_count + s.internal_count, s.node_count);
        assert_eq!(s.total_primitives, prims.len());
    }

    #[test]
    fn bvh_tree_stats_leaf_surface_area_positive() {
        let prims = make_grid_primitives(8);
        let bvh = Bvh::build(prims);
        let s = BvhTreeStatistics::compute(&bvh);
        assert!(
            s.total_leaf_surface_area > 0.0,
            "leaf surface area should be > 0"
        );
    }

    // ------------------------------------------------------------------
    // build_morton_clusters tests
    // ------------------------------------------------------------------

    #[test]
    fn build_morton_clusters_empty() {
        let clusters = build_morton_clusters(&[], 4);
        assert!(clusters.is_empty());
    }

    #[test]
    fn build_morton_clusters_count() {
        let sorted = make_sorted_lbvh_prims(10);
        let clusters = build_morton_clusters(&sorted, 3);
        // 10 primitives in chunks of 3 → ceil(10/3) = 4 clusters
        assert_eq!(clusters.len(), 4);
    }

    #[test]
    fn build_morton_clusters_radii_non_negative() {
        let sorted = make_sorted_lbvh_prims(8);
        let clusters = build_morton_clusters(&sorted, 2);
        for c in &clusters {
            assert!(c.radius >= 0.0, "cluster radius must be non-negative");
        }
    }
}
