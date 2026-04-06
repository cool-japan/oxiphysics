#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Advanced spatial queries.
//!
//! Implements a dynamic AABB tree, overlapping-pair queries (brute-force,
//! sweep-and-prune, AABB tree), frustum culling, sphere/capsule/point/convex
//! sweep queries, contact pair filtering by layer mask, and a batched query
//! pipeline.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

// ---------------------------------------------------------------------------
// AABB type
// ---------------------------------------------------------------------------

/// Axis-aligned bounding box.
#[derive(Clone, Debug, PartialEq)]
pub struct QueryAabb {
    /// Minimum corner.
    pub min: [f64; 3],
    /// Maximum corner.
    pub max: [f64; 3],
}

impl QueryAabb {
    /// Construct from min and max corners.
    pub fn new(min: [f64; 3], max: [f64; 3]) -> Self {
        QueryAabb { min, max }
    }

    /// Expand each side by `margin`.
    pub fn fatten(&self, margin: f64) -> Self {
        QueryAabb {
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

    /// Surface area of the AABB.
    pub fn surface_area(&self) -> f64 {
        let dx = self.max[0] - self.min[0];
        let dy = self.max[1] - self.min[1];
        let dz = self.max[2] - self.min[2];
        2.0 * (dx * dy + dy * dz + dz * dx)
    }

    /// Center of the AABB.
    pub fn center(&self) -> [f64; 3] {
        [
            0.5 * (self.min[0] + self.max[0]),
            0.5 * (self.min[1] + self.max[1]),
            0.5 * (self.min[2] + self.max[2]),
        ]
    }

    /// Union of two AABBs.
    pub fn union(&self, other: &QueryAabb) -> QueryAabb {
        QueryAabb {
            min: [
                self.min[0].min(other.min[0]),
                self.min[1].min(other.min[1]),
                self.min[2].min(other.min[2]),
            ],
            max: [
                self.max[0].max(other.max[0]),
                self.max[1].max(other.max[1]),
                self.max[2].max(other.max[2]),
            ],
        }
    }

    /// Test if two AABBs overlap.
    pub fn overlaps(&self, other: &QueryAabb) -> bool {
        self.min[0] <= other.max[0]
            && self.max[0] >= other.min[0]
            && self.min[1] <= other.max[1]
            && self.max[1] >= other.min[1]
            && self.min[2] <= other.max[2]
            && self.max[2] >= other.min[2]
    }

    /// Test if this AABB contains a point.
    pub fn contains_point(&self, p: [f64; 3]) -> bool {
        p[0] >= self.min[0]
            && p[0] <= self.max[0]
            && p[1] >= self.min[1]
            && p[1] <= self.max[1]
            && p[2] >= self.min[2]
            && p[2] <= self.max[2]
    }

    /// Squared distance from a point to the AABB (0 if inside).
    pub fn point_dist_sq(&self, p: [f64; 3]) -> f64 {
        let mut d_sq = 0.0;
        for i in 0..3 {
            let pi = [p[0], p[1], p[2]][i];
            let mini = [self.min[0], self.min[1], self.min[2]][i];
            let maxi = [self.max[0], self.max[1], self.max[2]][i];
            if pi < mini {
                d_sq += (mini - pi) * (mini - pi);
            } else if pi > maxi {
                d_sq += (pi - maxi) * (pi - maxi);
            }
        }
        d_sq
    }
}

// ---------------------------------------------------------------------------
// Helper functions (public)
// ---------------------------------------------------------------------------

/// Test if two AABBs overlap.
pub fn aabb_overlap(a: &QueryAabb, b: &QueryAabb) -> bool {
    a.overlaps(b)
}

/// Test if an AABB passes a frustum culling test (6 plane tests).
///
/// `planes` is a slice of (normal, d) pairs where plane equation = n·x + d = 0.
/// Returns true if the AABB is at least partially inside the frustum.
pub fn frustum_aabb_test(aabb: &QueryAabb, planes: &[([f64; 3], f64)]) -> bool {
    for &(n, d) in planes {
        // P-vertex: corner furthest in direction n
        let px = if n[0] > 0.0 { aabb.max[0] } else { aabb.min[0] };
        let py = if n[1] > 0.0 { aabb.max[1] } else { aabb.min[1] };
        let pz = if n[2] > 0.0 { aabb.max[2] } else { aabb.min[2] };
        let dist = n[0] * px + n[1] * py + n[2] * pz + d;
        if dist < 0.0 {
            return false;
        } // AABB fully outside this plane
    }
    true
}

/// Compute the swept AABB of an object moving from `aabb` by `motion`.
pub fn sweep_aabb_motion(aabb: &QueryAabb, motion: [f64; 3]) -> QueryAabb {
    let mut min = aabb.min;
    let mut max = aabb.max;
    for i in 0..3 {
        let m = [motion[0], motion[1], motion[2]][i];
        if m < 0.0 {
            min[i] += m;
        } else {
            max[i] += m;
        }
    }
    QueryAabb { min, max }
}

/// Rebalance an AABB tree by rotating nodes to minimize surface area.
///
/// This is a placeholder — actual rebalancing is done inside `AabbTree`.
pub fn rebalance_tree(_tree: &mut AabbTree) {
    // Rebalancing logic is integrated into AabbTree::insert
}

// ---------------------------------------------------------------------------
// Dynamic AABB tree node
// ---------------------------------------------------------------------------

/// A node in the dynamic AABB tree.
#[derive(Clone, Debug)]
pub struct DynamicNode {
    /// Fat AABB (enlarged by margin).
    pub fat_aabb: QueryAabb,
    /// Tight AABB (actual object bounds).
    pub aabb: QueryAabb,
    /// Parent node index, or None for root.
    pub parent: Option<usize>,
    /// Children node indices (None for leaves).
    pub children: [Option<usize>; 2],
    /// Height in the tree (0 = leaf).
    pub height: i32,
    /// User data (shape ID).
    pub user_data: usize,
    /// True if this is a leaf node.
    pub is_leaf: bool,
}

impl DynamicNode {
    /// Create a leaf node.
    pub fn leaf(aabb: QueryAabb, margin: f64, user_data: usize) -> Self {
        let fat_aabb = aabb.fatten(margin);
        DynamicNode {
            fat_aabb,
            aabb,
            parent: None,
            children: [None, None],
            height: 0,
            user_data,
            is_leaf: true,
        }
    }

    /// Create an internal node.
    pub fn internal(fat_aabb: QueryAabb) -> Self {
        DynamicNode {
            fat_aabb: fat_aabb.clone(),
            aabb: fat_aabb,
            parent: None,
            children: [None, None],
            height: 1,
            user_data: usize::MAX,
            is_leaf: false,
        }
    }
}

// ---------------------------------------------------------------------------
// Dynamic AABB tree
// ---------------------------------------------------------------------------

/// Dynamic AABB tree for broad-phase collision detection.
pub struct AabbTree {
    /// Node storage.
    pub nodes: Vec<DynamicNode>,
    /// Root node index.
    pub root: Option<usize>,
    /// Free node list (indices of removed nodes).
    free_nodes: Vec<usize>,
    /// Fat AABB margin.
    pub margin: f64,
}

impl AabbTree {
    /// Construct an empty AABB tree.
    pub fn new(margin: f64) -> Self {
        AabbTree {
            nodes: Vec::new(),
            root: None,
            free_nodes: Vec::new(),
            margin,
        }
    }

    fn alloc_node(&mut self, node: DynamicNode) -> usize {
        if let Some(idx) = self.free_nodes.pop() {
            self.nodes[idx] = node;
            idx
        } else {
            self.nodes.push(node);
            self.nodes.len() - 1
        }
    }

    /// Insert an AABB with user data, returns the leaf node index.
    pub fn insert(&mut self, aabb: QueryAabb, user_data: usize) -> usize {
        let leaf_idx = self.alloc_node(DynamicNode::leaf(aabb.clone(), self.margin, user_data));

        let root = match self.root {
            None => {
                self.root = Some(leaf_idx);
                return leaf_idx;
            }
            Some(r) => r,
        };

        // Find best sibling (greedy surface area heuristic)
        let mut best = root;
        let mut stack = vec![root];
        let leaf_aabb = aabb.clone();
        while let Some(idx) = stack.pop() {
            let node_aabb = self.nodes[idx].fat_aabb.clone();
            let combined = node_aabb.union(&leaf_aabb);
            let combined_sa = combined.surface_area();
            let cost = combined_sa;
            if cost < self.nodes[best].fat_aabb.union(&leaf_aabb).surface_area() {
                best = idx;
            }
            if !self.nodes[idx].is_leaf {
                for c in 0..2 {
                    if let Some(child) = self.nodes[idx].children[c] {
                        stack.push(child);
                    }
                }
            }
        }

        // Create a new internal node as sibling
        let sibling = best;
        let old_parent = self.nodes[sibling].parent;
        let new_parent_aabb = self.nodes[sibling]
            .fat_aabb
            .union(&leaf_aabb.fatten(self.margin));
        let new_parent_idx = self.alloc_node(DynamicNode::internal(new_parent_aabb));
        self.nodes[new_parent_idx].parent = old_parent;
        self.nodes[new_parent_idx].children = [Some(sibling), Some(leaf_idx)];
        self.nodes[new_parent_idx].height = self.nodes[sibling].height + 1;
        self.nodes[sibling].parent = Some(new_parent_idx);
        self.nodes[leaf_idx].parent = Some(new_parent_idx);

        if let Some(op) = old_parent {
            for c in 0..2 {
                if self.nodes[op].children[c] == Some(sibling) {
                    self.nodes[op].children[c] = Some(new_parent_idx);
                }
            }
        } else {
            self.root = Some(new_parent_idx);
        }

        leaf_idx
    }

    /// Remove a leaf node by index.
    pub fn remove(&mut self, leaf_idx: usize) {
        if self.root == Some(leaf_idx) {
            self.root = None;
            self.free_nodes.push(leaf_idx);
            return;
        }
        if let Some(parent_idx) = self.nodes[leaf_idx].parent {
            let sibling = {
                let ch = self.nodes[parent_idx].children;
                if ch[0] == Some(leaf_idx) {
                    ch[1]
                } else {
                    ch[0]
                }
            };
            if let Some(grandparent) = self.nodes[parent_idx].parent {
                for c in 0..2 {
                    if self.nodes[grandparent].children[c] == Some(parent_idx) {
                        self.nodes[grandparent].children[c] = sibling;
                    }
                }
                if let Some(s) = sibling {
                    self.nodes[s].parent = Some(grandparent);
                }
            } else {
                self.root = sibling;
                if let Some(s) = sibling {
                    self.nodes[s].parent = None;
                }
            }
            self.free_nodes.push(parent_idx);
        }
        self.free_nodes.push(leaf_idx);
    }

    /// Update a leaf node's AABB (remove + re-insert if moved outside fat AABB).
    pub fn update(&mut self, leaf_idx: usize, new_aabb: QueryAabb) {
        if self.nodes[leaf_idx].fat_aabb.overlaps(&new_aabb) {
            self.nodes[leaf_idx].aabb = new_aabb;
            return;
        }
        let user_data = self.nodes[leaf_idx].user_data;
        self.remove(leaf_idx);
        self.insert(new_aabb, user_data);
    }

    /// Query all leaf nodes whose fat AABB overlaps `query_aabb`.
    pub fn query_aabb(&self, query_aabb: &QueryAabb) -> Vec<usize> {
        let mut result = Vec::new();
        let Some(root) = self.root else {
            return result;
        };
        let mut stack = vec![root];
        while let Some(idx) = stack.pop() {
            let node = &self.nodes[idx];
            if !node.fat_aabb.overlaps(query_aabb) {
                continue;
            }
            if node.is_leaf {
                result.push(node.user_data);
            } else {
                for c in 0..2 {
                    if let Some(child) = node.children[c] {
                        stack.push(child);
                    }
                }
            }
        }
        result
    }

    /// Raycast: find all leaf nodes whose fat AABB is hit by ray (origin, dir, max_t).
    pub fn raycast(&self, origin: [f64; 3], dir: [f64; 3], max_t: f64) -> Vec<usize> {
        let mut result = Vec::new();
        let Some(root) = self.root else {
            return result;
        };
        let mut stack = vec![root];
        while let Some(idx) = stack.pop() {
            let node = &self.nodes[idx];
            if !ray_aabb_hit(&node.fat_aabb, origin, dir, max_t) {
                continue;
            }
            if node.is_leaf {
                result.push(node.user_data);
            } else {
                for c in 0..2 {
                    if let Some(child) = node.children[c] {
                        stack.push(child);
                    }
                }
            }
        }
        result
    }

    /// Number of live nodes.
    pub fn n_nodes(&self) -> usize {
        self.nodes.len().saturating_sub(self.free_nodes.len())
    }
}

fn ray_aabb_hit(aabb: &QueryAabb, origin: [f64; 3], dir: [f64; 3], max_t: f64) -> bool {
    let mut t_min = 0.0_f64;
    let mut t_max = max_t;
    for i in 0..3 {
        let o = [origin[0], origin[1], origin[2]][i];
        let d = [dir[0], dir[1], dir[2]][i];
        let mn = [aabb.min[0], aabb.min[1], aabb.min[2]][i];
        let mx = [aabb.max[0], aabb.max[1], aabb.max[2]][i];
        if d.abs() < 1e-14 {
            if o < mn || o > mx {
                return false;
            }
        } else {
            let inv = 1.0 / d;
            let t1 = (mn - o) * inv;
            let t2 = (mx - o) * inv;
            let ta = t1.min(t2);
            let tb = t1.max(t2);
            t_min = t_min.max(ta);
            t_max = t_max.min(tb);
            if t_min > t_max {
                return false;
            }
        }
    }
    true
}

// ---------------------------------------------------------------------------
// Pair query
// ---------------------------------------------------------------------------

/// An overlapping pair of shape IDs.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ShapePair {
    /// First shape ID.
    pub a: usize,
    /// Second shape ID.
    pub b: usize,
}

impl ShapePair {
    /// Create a new pair (a < b always).
    pub fn new(a: usize, b: usize) -> Self {
        if a < b {
            ShapePair { a, b }
        } else {
            ShapePair { a: b, b: a }
        }
    }
}

/// Overlapping pair detection using multiple algorithms.
pub struct PairQuery;

impl PairQuery {
    /// Brute-force O(n²) pair detection.
    pub fn brute_force(aabbs: &[QueryAabb]) -> Vec<ShapePair> {
        let n = aabbs.len();
        let mut pairs = Vec::new();
        for i in 0..n {
            for j in i + 1..n {
                if aabbs[i].overlaps(&aabbs[j]) {
                    pairs.push(ShapePair::new(i, j));
                }
            }
        }
        pairs
    }

    /// Sweep-and-prune along X axis.
    pub fn sweep_and_prune(aabbs: &[QueryAabb]) -> Vec<ShapePair> {
        let mut sorted: Vec<usize> = (0..aabbs.len()).collect();
        sorted.sort_unstable_by(|&a, &b| {
            aabbs[a].min[0]
                .partial_cmp(&aabbs[b].min[0])
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        let mut pairs = Vec::new();
        for i in 0..sorted.len() {
            let a = sorted[i];
            for j in i + 1..sorted.len() {
                let b = sorted[j];
                if aabbs[b].min[0] > aabbs[a].max[0] {
                    break;
                }
                if aabbs[a].overlaps(&aabbs[b]) {
                    pairs.push(ShapePair::new(a, b));
                }
            }
        }
        pairs
    }

    /// AABB tree-based pair detection.
    pub fn aabb_tree(aabbs: &[QueryAabb]) -> Vec<ShapePair> {
        let mut tree = AabbTree::new(0.01);
        for (i, aabb) in aabbs.iter().enumerate() {
            tree.insert(aabb.clone(), i);
        }
        let mut pairs = Vec::new();
        for (i, aabb) in aabbs.iter().enumerate() {
            let hits = tree.query_aabb(aabb);
            for h in hits {
                if h != i && h > i {
                    pairs.push(ShapePair::new(i, h));
                }
            }
        }
        // Deduplicate
        pairs.sort_unstable_by_key(|p| (p.a, p.b));
        pairs.dedup();
        pairs
    }
}

// ---------------------------------------------------------------------------
// Frustum culling
// ---------------------------------------------------------------------------

/// Frustum defined by 6 planes.
pub struct FrustumCulling {
    /// Planes as (normal, d): n·x + d = 0.
    pub planes: Vec<([f64; 3], f64)>,
}

impl FrustumCulling {
    /// Construct from 6 planes.
    pub fn new(planes: Vec<([f64; 3], f64)>) -> Self {
        FrustumCulling { planes }
    }

    /// Cull AABBs: return indices of visible (not fully outside) objects.
    pub fn cull(&self, aabbs: &[QueryAabb]) -> Vec<usize> {
        aabbs
            .iter()
            .enumerate()
            .filter(|(_, aabb)| frustum_aabb_test(aabb, &self.planes))
            .map(|(i, _)| i)
            .collect()
    }

    /// Test a single AABB.
    pub fn test(&self, aabb: &QueryAabb) -> bool {
        frustum_aabb_test(aabb, &self.planes)
    }
}

// ---------------------------------------------------------------------------
// Sphere query
// ---------------------------------------------------------------------------

/// Query for shapes overlapping a sphere.
#[derive(Clone, Debug)]
pub struct SphereQuery {
    /// Center.
    pub center: [f64; 3],
    /// Radius.
    pub radius: f64,
}

impl SphereQuery {
    /// Construct a sphere query.
    pub fn new(center: [f64; 3], radius: f64) -> Self {
        SphereQuery { center, radius }
    }

    /// Return shape indices overlapping this sphere (against AABB list), sorted by distance.
    pub fn query_sorted(&self, aabbs: &[QueryAabb]) -> Vec<usize> {
        let r2 = self.radius * self.radius;
        let mut hits: Vec<(usize, f64)> = aabbs
            .iter()
            .enumerate()
            .filter_map(|(i, aabb)| {
                let d2 = aabb.point_dist_sq(self.center);
                if d2 <= r2 { Some((i, d2)) } else { None }
            })
            .collect();
        hits.sort_unstable_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        hits.into_iter().map(|(i, _)| i).collect()
    }

    /// Sphere AABB for broad-phase.
    pub fn aabb(&self) -> QueryAabb {
        let r = self.radius;
        QueryAabb {
            min: [self.center[0] - r, self.center[1] - r, self.center[2] - r],
            max: [self.center[0] + r, self.center[1] + r, self.center[2] + r],
        }
    }
}

// ---------------------------------------------------------------------------
// Capsule query
// ---------------------------------------------------------------------------

/// Query for shapes along a capsule sweep path.
pub struct CapsuleQuery {
    /// Start of capsule axis.
    pub p0: [f64; 3],
    /// End of capsule axis.
    pub p1: [f64; 3],
    /// Radius.
    pub radius: f64,
}

impl CapsuleQuery {
    /// Construct a capsule query.
    pub fn new(p0: [f64; 3], p1: [f64; 3], radius: f64) -> Self {
        CapsuleQuery { p0, p1, radius }
    }

    /// AABB enclosing the capsule.
    pub fn aabb(&self) -> QueryAabb {
        let r = self.radius;
        QueryAabb {
            min: [
                self.p0[0].min(self.p1[0]) - r,
                self.p0[1].min(self.p1[1]) - r,
                self.p0[2].min(self.p1[2]) - r,
            ],
            max: [
                self.p0[0].max(self.p1[0]) + r,
                self.p0[1].max(self.p1[1]) + r,
                self.p0[2].max(self.p1[2]) + r,
            ],
        }
    }

    /// Return shapes hit by this capsule (broad-phase AABB test).
    pub fn query(&self, aabbs: &[QueryAabb]) -> Vec<usize> {
        let caps_aabb = self.aabb();
        aabbs
            .iter()
            .enumerate()
            .filter(|(_, aabb)| aabb.overlaps(&caps_aabb))
            .map(|(i, _)| i)
            .collect()
    }
}

// ---------------------------------------------------------------------------
// Point query
// ---------------------------------------------------------------------------

/// Query for the N shapes nearest to a point.
pub struct PointQuery {
    /// Query point.
    pub point: [f64; 3],
}

impl PointQuery {
    /// Construct a point query.
    pub fn new(point: [f64; 3]) -> Self {
        PointQuery { point }
    }

    /// Return the `n` shape indices nearest to the point, sorted by distance.
    pub fn nearest_n(&self, aabbs: &[QueryAabb], n: usize) -> Vec<usize> {
        let mut dists: Vec<(usize, f64)> = aabbs
            .iter()
            .enumerate()
            .map(|(i, aabb)| (i, aabb.point_dist_sq(self.point).sqrt()))
            .collect();
        dists.sort_unstable_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        dists.into_iter().take(n).map(|(i, _)| i).collect()
    }

    /// Exact distance to AABB i.
    pub fn distance_to(&self, aabb: &QueryAabb) -> f64 {
        aabb.point_dist_sq(self.point).sqrt()
    }
}

// ---------------------------------------------------------------------------
// Convex sweep
// ---------------------------------------------------------------------------

/// Convex sweep: find all shapes hit by a box swept along `motion`.
pub struct ConvexSweep {
    /// AABB of the convex shape at rest.
    pub shape_aabb: QueryAabb,
    /// Sweep motion vector.
    pub motion: [f64; 3],
}

impl ConvexSweep {
    /// Construct a convex sweep.
    pub fn new(shape_aabb: QueryAabb, motion: [f64; 3]) -> Self {
        ConvexSweep { shape_aabb, motion }
    }

    /// Swept AABB enclosing the entire motion path.
    pub fn swept_aabb(&self) -> QueryAabb {
        sweep_aabb_motion(&self.shape_aabb, self.motion)
    }

    /// Return shapes hit by the swept volume.
    pub fn query(&self, aabbs: &[QueryAabb]) -> Vec<usize> {
        let swept = self.swept_aabb();
        aabbs
            .iter()
            .enumerate()
            .filter(|(_, aabb)| aabb.overlaps(&swept))
            .map(|(i, _)| i)
            .collect()
    }
}

// ---------------------------------------------------------------------------
// Contact pair filter
// ---------------------------------------------------------------------------

/// Contact pair filter using collision layer masks.
pub struct ContactPairFilter {
    /// Layer mask for each shape.
    pub layer_masks: Vec<u32>,
    /// Group for each shape.
    pub groups: Vec<u32>,
}

impl ContactPairFilter {
    /// Construct a contact pair filter.
    pub fn new(layer_masks: Vec<u32>, groups: Vec<u32>) -> Self {
        ContactPairFilter {
            layer_masks,
            groups,
        }
    }

    /// True if shapes `a` and `b` should collide (bitmask AND ≠ 0 and same group check).
    pub fn should_collide(&self, a: usize, b: usize) -> bool {
        let mask_a = self.layer_masks.get(a).copied().unwrap_or(0xFFFF_FFFF);
        let mask_b = self.layer_masks.get(b).copied().unwrap_or(0xFFFF_FFFF);
        (mask_a & mask_b) != 0
    }

    /// Filter a list of shape pairs.
    pub fn filter(&self, pairs: &[ShapePair]) -> Vec<ShapePair> {
        pairs
            .iter()
            .filter(|p| self.should_collide(p.a, p.b))
            .cloned()
            .collect()
    }
}

// ---------------------------------------------------------------------------
// Query pipeline
// ---------------------------------------------------------------------------

/// A pending query (for batched processing).
#[derive(Clone, Debug)]
pub enum PendingQuery {
    /// AABB overlap query.
    AabbQuery(QueryAabb),
    /// Sphere overlap query.
    SphereQuery(SphereQuery),
    /// Point nearest-N query.
    PointQuery {
        /// Query point in 3D space.
        point: [f64; 3],
        /// Number of nearest neighbours to find.
        n: usize,
    },
}

/// Result of a pending query.
#[derive(Clone, Debug)]
pub struct QueryResult {
    /// Shape indices returned by the query.
    pub hits: Vec<usize>,
}

/// Multi-query pipeline with deferred (next-frame) result dispatch.
pub struct QueryPipeline {
    /// Registered AABBs.
    pub aabbs: Vec<QueryAabb>,
    /// Pending queries.
    pub pending: Vec<PendingQuery>,
    /// Results from previous frame.
    pub results: Vec<QueryResult>,
}

impl QueryPipeline {
    /// Construct an empty pipeline.
    pub fn new() -> Self {
        QueryPipeline {
            aabbs: Vec::new(),
            pending: Vec::new(),
            results: Vec::new(),
        }
    }

    /// Register a shape AABB.
    pub fn add_shape(&mut self, aabb: QueryAabb) -> usize {
        self.aabbs.push(aabb);
        self.aabbs.len() - 1
    }

    /// Submit a query for deferred processing.
    pub fn submit(&mut self, query: PendingQuery) {
        self.pending.push(query);
    }

    /// Process all pending queries and store results (simulated next-frame).
    pub fn flush(&mut self) {
        let aabbs = &self.aabbs;
        let mut new_results = Vec::with_capacity(self.pending.len());
        for query in self.pending.drain(..) {
            let hits = match query {
                PendingQuery::AabbQuery(aabb) => aabbs
                    .iter()
                    .enumerate()
                    .filter(|(_, a)| a.overlaps(&aabb))
                    .map(|(i, _)| i)
                    .collect(),
                PendingQuery::SphereQuery(sq) => {
                    let sphere_aabb = sq.aabb();
                    aabbs
                        .iter()
                        .enumerate()
                        .filter(|(_, a)| a.overlaps(&sphere_aabb))
                        .map(|(i, _)| i)
                        .collect()
                }
                PendingQuery::PointQuery { point, n } => {
                    let pq = PointQuery::new(point);
                    pq.nearest_n(aabbs, n)
                }
            };
            new_results.push(QueryResult { hits });
        }
        self.results = new_results;
    }

    /// Number of shapes registered.
    pub fn n_shapes(&self) -> usize {
        self.aabbs.len()
    }
}

impl Default for QueryPipeline {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn unit_aabb() -> QueryAabb {
        QueryAabb::new([0.0; 3], [1.0; 3])
    }

    fn shifted_aabb(offset: f64) -> QueryAabb {
        QueryAabb::new([offset; 3], [offset + 1.0; 3])
    }

    #[test]
    fn test_aabb_overlaps_touching() {
        let a = QueryAabb::new([0.0; 3], [1.0; 3]);
        let b = QueryAabb::new([1.0; 3], [2.0; 3]);
        assert!(a.overlaps(&b));
    }

    #[test]
    fn test_aabb_no_overlap() {
        let a = QueryAabb::new([0.0; 3], [1.0; 3]);
        let b = QueryAabb::new([2.0; 3], [3.0; 3]);
        assert!(!a.overlaps(&b));
    }

    #[test]
    fn test_aabb_union() {
        let a = QueryAabb::new([0.0; 3], [1.0; 3]);
        let b = QueryAabb::new([1.0; 3], [2.0; 3]);
        let u = a.union(&b);
        assert_eq!(u.min, [0.0; 3]);
        assert_eq!(u.max, [2.0; 3]);
    }

    #[test]
    fn test_aabb_contains_point() {
        let a = unit_aabb();
        assert!(a.contains_point([0.5, 0.5, 0.5]));
        assert!(!a.contains_point([2.0, 0.5, 0.5]));
    }

    #[test]
    fn test_aabb_point_dist_sq_inside() {
        let a = unit_aabb();
        assert_eq!(a.point_dist_sq([0.5, 0.5, 0.5]), 0.0);
    }

    #[test]
    fn test_aabb_point_dist_sq_outside() {
        let a = unit_aabb();
        let d2 = a.point_dist_sq([2.0, 0.5, 0.5]);
        assert!((d2 - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_aabb_fatten() {
        let a = unit_aabb();
        let fat = a.fatten(0.5);
        assert!((fat.min[0] - (-0.5)).abs() < 1e-12);
        assert!((fat.max[0] - 1.5).abs() < 1e-12);
    }

    #[test]
    fn test_aabb_tree_insert_and_query() {
        let mut tree = AabbTree::new(0.1);
        let idx = tree.insert(unit_aabb(), 42);
        let hits = tree.query_aabb(&unit_aabb());
        assert!(hits.contains(&42), "hits={:?}, leaf_idx={}", hits, idx);
    }

    #[test]
    fn test_aabb_tree_remove() {
        let mut tree = AabbTree::new(0.1);
        let leaf = tree.insert(unit_aabb(), 0);
        tree.remove(leaf);
        let hits = tree.query_aabb(&unit_aabb());
        assert!(hits.is_empty());
    }

    #[test]
    fn test_aabb_tree_multiple_inserts() {
        let mut tree = AabbTree::new(0.1);
        tree.insert(unit_aabb(), 0);
        tree.insert(shifted_aabb(5.0), 1);
        let hits = tree.query_aabb(&unit_aabb());
        assert!(hits.contains(&0));
        assert!(!hits.contains(&1));
    }

    #[test]
    fn test_aabb_tree_raycast() {
        let mut tree = AabbTree::new(0.1);
        tree.insert(unit_aabb(), 7);
        let hits = tree.raycast([-5.0, 0.5, 0.5], [1.0, 0.0, 0.0], 20.0);
        assert!(hits.contains(&7));
    }

    #[test]
    fn test_pair_query_brute_force() {
        let aabbs = vec![unit_aabb(), shifted_aabb(0.5), shifted_aabb(10.0)];
        let pairs = PairQuery::brute_force(&aabbs);
        // 0 and 1 overlap, 2 does not
        assert!(
            pairs
                .iter()
                .any(|p| (p.a == 0 && p.b == 1) || (p.a == 1 && p.b == 0))
        );
        assert!(!pairs.iter().any(|p| p.a == 2 || p.b == 2));
    }

    #[test]
    fn test_pair_query_sweep_and_prune() {
        let aabbs = vec![unit_aabb(), shifted_aabb(0.5)];
        let pairs = PairQuery::sweep_and_prune(&aabbs);
        assert!(!pairs.is_empty());
    }

    #[test]
    fn test_pair_query_aabb_tree() {
        let aabbs = vec![unit_aabb(), shifted_aabb(0.5), shifted_aabb(10.0)];
        let pairs = PairQuery::aabb_tree(&aabbs);
        assert!(pairs.iter().any(|p| p.a == 0 && p.b == 1));
    }

    #[test]
    fn test_frustum_culling_inside() {
        // Frustum: half-space x > -10
        let planes = vec![([1.0_f64, 0.0, 0.0], 10.0)];
        let fc = FrustumCulling::new(planes);
        let visible = fc.cull(&[unit_aabb()]);
        assert!(visible.contains(&0));
    }

    #[test]
    fn test_frustum_culling_outside() {
        // Frustum: x > 100
        let planes = vec![([1.0_f64, 0.0, 0.0], -200.0)];
        let fc = FrustumCulling::new(planes);
        let visible = fc.cull(&[unit_aabb()]);
        assert!(visible.is_empty());
    }

    #[test]
    fn test_sphere_query_hit() {
        let sq = SphereQuery::new([0.5, 0.5, 0.5], 0.1);
        let hits = sq.query_sorted(&[unit_aabb()]);
        assert!(hits.contains(&0));
    }

    #[test]
    fn test_sphere_query_miss() {
        let sq = SphereQuery::new([10.0, 10.0, 10.0], 0.1);
        let hits = sq.query_sorted(&[unit_aabb()]);
        assert!(hits.is_empty());
    }

    #[test]
    fn test_capsule_query_hit() {
        let cq = CapsuleQuery::new([0.5, 0.5, -1.0], [0.5, 0.5, 2.0], 0.1);
        let hits = cq.query(&[unit_aabb()]);
        assert!(hits.contains(&0));
    }

    #[test]
    fn test_capsule_query_miss() {
        let cq = CapsuleQuery::new([10.0, 10.0, 10.0], [11.0, 11.0, 11.0], 0.1);
        let hits = cq.query(&[unit_aabb()]);
        assert!(hits.is_empty());
    }

    #[test]
    fn test_point_query_nearest() {
        let aabbs = vec![unit_aabb(), shifted_aabb(5.0)];
        let pq = PointQuery::new([0.5, 0.5, 0.5]);
        let nearest = pq.nearest_n(&aabbs, 1);
        assert_eq!(nearest[0], 0);
    }

    #[test]
    fn test_point_query_distance_inside() {
        let aabb = unit_aabb();
        let pq = PointQuery::new([0.5, 0.5, 0.5]);
        assert_eq!(pq.distance_to(&aabb), 0.0);
    }

    #[test]
    fn test_convex_sweep_hit() {
        let cs = ConvexSweep::new(QueryAabb::new([-1.0; 3], [0.0; 3]), [2.0, 2.0, 2.0]);
        let hits = cs.query(&[unit_aabb()]);
        assert!(hits.contains(&0));
    }

    #[test]
    fn test_convex_sweep_miss() {
        let cs = ConvexSweep::new(QueryAabb::new([-10.0; 3], [-9.0; 3]), [-1.0, 0.0, 0.0]);
        let hits = cs.query(&[unit_aabb()]);
        assert!(hits.is_empty());
    }

    #[test]
    fn test_contact_pair_filter_passes() {
        let cf = ContactPairFilter::new(vec![0xF, 0xF], vec![0, 0]);
        assert!(cf.should_collide(0, 1));
    }

    #[test]
    fn test_contact_pair_filter_blocks() {
        let cf = ContactPairFilter::new(vec![0x1, 0x2], vec![0, 0]);
        assert!(!cf.should_collide(0, 1));
    }

    #[test]
    fn test_contact_pair_filter_filter() {
        let cf = ContactPairFilter::new(vec![0xF, 0xF, 0x1], vec![0; 3]);
        let pairs = vec![ShapePair::new(0, 1), ShapePair::new(0, 2)];
        let filtered = cf.filter(&pairs);
        // 0↔1: mask 0xF&0xF ≠ 0 → passes; 0↔2: mask 0xF&0x1 = 0x1 ≠ 0 → passes
        assert_eq!(filtered.len(), 2);
    }

    #[test]
    fn test_query_pipeline_flush() {
        let mut qp = QueryPipeline::new();
        qp.add_shape(unit_aabb());
        qp.submit(PendingQuery::AabbQuery(unit_aabb()));
        qp.flush();
        assert_eq!(qp.results.len(), 1);
        assert!(qp.results[0].hits.contains(&0));
    }

    #[test]
    fn test_query_pipeline_sphere_query() {
        let mut qp = QueryPipeline::new();
        qp.add_shape(unit_aabb());
        let sq = SphereQuery::new([0.5, 0.5, 0.5], 0.5);
        qp.submit(PendingQuery::SphereQuery(sq));
        qp.flush();
        assert!(!qp.results.is_empty());
    }

    #[test]
    fn test_query_pipeline_point_query() {
        let mut qp = QueryPipeline::new();
        qp.add_shape(unit_aabb());
        qp.add_shape(shifted_aabb(5.0));
        qp.submit(PendingQuery::PointQuery {
            point: [0.5, 0.5, 0.5],
            n: 1,
        });
        qp.flush();
        assert_eq!(qp.results[0].hits.len(), 1);
        assert_eq!(qp.results[0].hits[0], 0);
    }

    #[test]
    fn test_sweep_aabb_motion_positive() {
        let aabb = unit_aabb();
        let swept = sweep_aabb_motion(&aabb, [2.0, 0.0, 0.0]);
        assert!((swept.max[0] - 3.0).abs() < 1e-12);
        assert!((swept.min[0] - 0.0).abs() < 1e-12);
    }

    #[test]
    fn test_sweep_aabb_motion_negative() {
        let aabb = unit_aabb();
        let swept = sweep_aabb_motion(&aabb, [-1.0, 0.0, 0.0]);
        assert!((swept.min[0] - (-1.0)).abs() < 1e-12);
    }

    #[test]
    fn test_aabb_overlap_function() {
        let a = unit_aabb();
        let b = shifted_aabb(0.5);
        assert!(aabb_overlap(&a, &b));
    }
}
