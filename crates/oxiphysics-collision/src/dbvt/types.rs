//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use oxiphysics_core::math::Vec3;

#[allow(unused_imports)]
use super::functions::*;
use super::functions::{FAT_MARGIN, NodeIdx};

/// A node in the dynamic BVH tree.
pub(super) struct DbvtNode {
    /// Bounding box of this node (fat for leaves, tight for internal nodes).
    pub(super) aabb: BvhAabb,
    /// Parent index, or `None` for the root.
    pub(super) parent: Option<NodeIdx>,
    /// Left child (internal nodes only).
    pub(super) left: Option<NodeIdx>,
    /// Right child (internal nodes only).
    pub(super) right: Option<NodeIdx>,
    /// Leaf payload — `None` for internal nodes.
    pub(super) data: Option<u32>,
    /// Sub-tree height (0 for leaves).
    pub(super) height: i32,
}
impl DbvtNode {
    fn leaf(aabb: BvhAabb, data: u32) -> Self {
        Self {
            aabb,
            parent: None,
            left: None,
            right: None,
            data: Some(data),
            height: 0,
        }
    }
    fn internal(aabb: BvhAabb, left: NodeIdx, right: NodeIdx) -> Self {
        Self {
            aabb,
            parent: None,
            left: Some(left),
            right: Some(right),
            data: None,
            height: 1,
        }
    }
    pub(crate) fn is_leaf(&self) -> bool {
        self.data.is_some()
    }
}
/// One entry in a k-nearest result: (leaf data, squared distance to query point).
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct NearestLeaf {
    /// Leaf payload data.
    pub data: u32,
    /// Squared distance from the query point to the AABB center of the leaf.
    pub dist_sq: f64,
}
/// Axis-Aligned Bounding Box used by the dynamic BVH tree.
#[derive(Debug, Clone, Copy)]
pub struct BvhAabb {
    /// Minimum corner of the bounding box.
    pub min: Vec3,
    /// Maximum corner of the bounding box.
    pub max: Vec3,
}
impl BvhAabb {
    /// Create a new AABB from min/max corners.
    #[inline]
    pub fn new(min: Vec3, max: Vec3) -> Self {
        Self { min, max }
    }
    /// Create a zero-size AABB at a single point.
    #[inline]
    pub fn point(p: Vec3) -> Self {
        Self { min: p, max: p }
    }
    /// Create an AABB from a sphere (center + radius).
    #[inline]
    pub fn from_sphere(center: Vec3, radius: f64) -> Self {
        let r = Vec3::new(radius, radius, radius);
        Self {
            min: center - r,
            max: center + r,
        }
    }
    /// Center of the AABB.
    #[inline]
    pub fn center(&self) -> Vec3 {
        (self.min + self.max) * 0.5
    }
    /// Half-extents (half-widths along each axis).
    #[inline]
    pub fn half_extents(&self) -> Vec3 {
        (self.max - self.min) * 0.5
    }
    /// Surface area: `2*(dx*dy + dy*dz + dz*dx)`.
    #[inline]
    pub fn surface_area(&self) -> f64 {
        let d = self.max - self.min;
        2.0 * (d.x * d.y + d.y * d.z + d.z * d.x)
    }
    /// Volume of the AABB.
    #[inline]
    pub fn volume(&self) -> f64 {
        let d = self.max - self.min;
        d.x * d.y * d.z
    }
    /// Return `true` if `self` fully contains `other`.
    #[inline]
    pub fn contains(&self, other: &Self) -> bool {
        self.min.x <= other.min.x
            && self.min.y <= other.min.y
            && self.min.z <= other.min.z
            && self.max.x >= other.max.x
            && self.max.y >= other.max.y
            && self.max.z >= other.max.z
    }
    /// Return `true` if `self` overlaps `other`.
    #[inline]
    pub fn intersects(&self, other: &Self) -> bool {
        self.min.x <= other.max.x
            && self.max.x >= other.min.x
            && self.min.y <= other.max.y
            && self.max.y >= other.min.y
            && self.min.z <= other.max.z
            && self.max.z >= other.min.z
    }
    /// Smallest AABB that contains both `a` and `b`.
    #[inline]
    pub fn merge(a: &Self, b: &Self) -> Self {
        Self {
            min: Vec3::new(
                a.min.x.min(b.min.x),
                a.min.y.min(b.min.y),
                a.min.z.min(b.min.z),
            ),
            max: Vec3::new(
                a.max.x.max(b.max.x),
                a.max.y.max(b.max.y),
                a.max.z.max(b.max.z),
            ),
        }
    }
    /// Expand this AABB uniformly by `margin` in every direction.
    #[inline]
    pub fn expand(&self, margin: f64) -> Self {
        let d = Vec3::new(margin, margin, margin);
        Self {
            min: self.min - d,
            max: self.max + d,
        }
    }
    /// Return `true` if `point` is inside (or on the boundary of) this AABB.
    #[inline]
    #[allow(dead_code)]
    fn contains_point(&self, point: &Vec3) -> bool {
        point.x >= self.min.x
            && point.x <= self.max.x
            && point.y >= self.min.y
            && point.y <= self.max.y
            && point.z >= self.min.z
            && point.z <= self.max.z
    }
}
impl BvhAabb {
    /// Return the axis (0=X, 1=Y, 2=Z) along which the AABB is widest.
    #[inline]
    #[allow(dead_code)]
    pub fn longest_axis(&self) -> usize {
        let d = self.max - self.min;
        if d.x >= d.y && d.x >= d.z {
            0
        } else if d.y >= d.z {
            1
        } else {
            2
        }
    }
    /// Return `true` if `self` and `other` are the same AABB (within `eps`).
    #[inline]
    #[allow(dead_code)]
    pub fn approx_eq(&self, other: &Self, eps: f64) -> bool {
        (self.min.x - other.min.x).abs() <= eps
            && (self.min.y - other.min.y).abs() <= eps
            && (self.min.z - other.min.z).abs() <= eps
            && (self.max.x - other.max.x).abs() <= eps
            && (self.max.y - other.max.y).abs() <= eps
            && (self.max.z - other.max.z).abs() <= eps
    }
    /// Squared distance from `point` to the nearest point on this AABB.
    ///
    /// Returns `0.0` if `point` is inside the AABB.
    #[inline]
    #[allow(dead_code)]
    pub fn point_dist_sq(&self, point: Vec3) -> f64 {
        let dx = (point.x - self.min.x.max(self.max.x.min(point.x))).powi(2);
        let dy = (point.y - self.min.y.max(self.max.y.min(point.y))).powi(2);
        let dz = (point.z - self.min.z.max(self.max.z.min(point.z))).powi(2);
        dx + dy + dz
    }
    /// Translate the AABB by `offset`.
    #[inline]
    #[allow(dead_code)]
    pub fn translate(&self, offset: Vec3) -> Self {
        Self {
            min: self.min + offset,
            max: self.max + offset,
        }
    }
    /// Scale the AABB uniformly around its center by `factor`.
    #[inline]
    #[allow(dead_code)]
    pub fn scale(&self, factor: f64) -> Self {
        let c = self.center();
        let he = self.half_extents() * factor;
        Self {
            min: c - he,
            max: c + he,
        }
    }
}
impl BvhAabb {
    /// Test a ray `origin + t*dir` against this AABB for `t ∈ [t_min, t_max]`.
    ///
    /// Returns the entry `t` value if the ray hits, or `None` otherwise.
    #[allow(dead_code)]
    pub fn ray_intersect(&self, origin: Vec3, dir: Vec3, t_min: f64, t_max: f64) -> Option<f64> {
        let orig = [origin.x, origin.y, origin.z];
        let d = [dir.x, dir.y, dir.z];
        let mn = [self.min.x, self.min.y, self.min.z];
        let mx = [self.max.x, self.max.y, self.max.z];
        let mut tmin = t_min;
        let mut tmax = t_max;
        for i in 0..3 {
            if d[i].abs() < 1e-300 {
                if orig[i] < mn[i] || orig[i] > mx[i] {
                    return None;
                }
            } else {
                let inv = 1.0 / d[i];
                let t1 = (mn[i] - orig[i]) * inv;
                let t2 = (mx[i] - orig[i]) * inv;
                let (ta, tb) = if t1 <= t2 { (t1, t2) } else { (t2, t1) };
                tmin = tmin.max(ta);
                tmax = tmax.min(tb);
                if tmin > tmax {
                    return None;
                }
            }
        }
        if tmax < t_min {
            return None;
        }
        Some(tmin.max(t_min))
    }
}
/// Statistics collected from a DBVT query pass.
#[allow(dead_code)]
#[derive(Debug, Clone, Default)]
pub struct DbvtStats {
    /// Number of leaf nodes in the tree.
    pub leaf_count: usize,
    /// Total node count (leaves + internal).
    pub node_count: usize,
    /// Tree height.
    pub height: i32,
    /// Surface area of the root AABB.
    pub root_surface_area: f64,
    /// Number of overlap pairs found by last `query_overlapping_pairs`.
    pub overlap_pairs: usize,
}
/// A view frustum represented by six half-space planes (inward normals).
///
/// Each plane `(n, d)` passes the test when `n · p >= d`.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct BvhFrustum {
    /// Six planes: (inward_normal, offset).
    pub planes: [(Vec3, f64); 6],
}
#[allow(dead_code)]
impl BvhFrustum {
    /// Construct a frustum from six plane definitions.
    pub fn new(planes: [(Vec3, f64); 6]) -> Self {
        Self { planes }
    }
    /// Test whether an AABB is potentially inside (not fully outside) the frustum.
    pub fn intersects_aabb(&self, aabb: &BvhAabb) -> bool {
        for &(n, d) in &self.planes {
            let px = if n.x >= 0.0 { aabb.max.x } else { aabb.min.x };
            let py = if n.y >= 0.0 { aabb.max.y } else { aabb.min.y };
            let pz = if n.z >= 0.0 { aabb.max.z } else { aabb.min.z };
            let dot = n.x * px + n.y * py + n.z * pz;
            if dot < d {
                return false;
            }
        }
        true
    }
}
/// A capsule shape: a line segment swept by a sphere.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct BvhCapsule {
    /// Start point of the central segment.
    pub start: Vec3,
    /// End point of the central segment.
    pub end: Vec3,
    /// Radius of the capsule.
    pub radius: f64,
}
/// Dynamic AABB Tree (DBVT) for broadphase collision detection.
///
/// Supports incremental `insert`, `remove`, and `update` operations.
/// Uses the surface-area heuristic (SAH) to choose the best sibling when
/// inserting a new leaf.
pub struct DynamicBvh {
    pub(super) nodes: Vec<DbvtNode>,
    pub(super) root: Option<NodeIdx>,
    pub(super) free_list: Vec<NodeIdx>,
}
impl DynamicBvh {
    /// Create a new, empty dynamic BVH.
    pub fn new() -> Self {
        Self {
            nodes: Vec::new(),
            root: None,
            free_list: Vec::new(),
        }
    }
    fn alloc_node(&mut self, node: DbvtNode) -> NodeIdx {
        if let Some(idx) = self.free_list.pop() {
            self.nodes[idx] = node;
            idx
        } else {
            let idx = self.nodes.len();
            self.nodes.push(node);
            idx
        }
    }
    fn free_node(&mut self, idx: NodeIdx) {
        self.nodes[idx].parent = None;
        self.nodes[idx].left = None;
        self.nodes[idx].right = None;
        self.nodes[idx].data = None;
        self.free_list.push(idx);
    }
    /// Insert a leaf AABB with associated `data`. Returns the leaf node index.
    ///
    /// The leaf is stored with a fat AABB (expanded by `FAT_MARGIN`) to
    /// reduce the frequency of re-insertions when the object moves slightly.
    pub fn insert(&mut self, aabb: BvhAabb, data: u32) -> NodeIdx {
        let fat_aabb = aabb.expand(FAT_MARGIN);
        let leaf_idx = self.alloc_node(DbvtNode::leaf(fat_aabb, data));
        if self.root.is_none() {
            self.root = Some(leaf_idx);
            return leaf_idx;
        }
        let sibling = self.find_best_sibling(leaf_idx);
        let old_parent = self.nodes[sibling].parent;
        let new_aabb = BvhAabb::merge(&self.nodes[sibling].aabb, &self.nodes[leaf_idx].aabb);
        let new_parent_idx = self.alloc_node(DbvtNode::internal(new_aabb, sibling, leaf_idx));
        self.nodes[new_parent_idx].parent = old_parent;
        self.nodes[sibling].parent = Some(new_parent_idx);
        self.nodes[leaf_idx].parent = Some(new_parent_idx);
        match old_parent {
            None => {
                self.root = Some(new_parent_idx);
            }
            Some(gp) => {
                if self.nodes[gp].left == Some(sibling) {
                    self.nodes[gp].left = Some(new_parent_idx);
                } else {
                    self.nodes[gp].right = Some(new_parent_idx);
                }
            }
        }
        self.refit(new_parent_idx);
        leaf_idx
    }
    /// Remove a leaf node by its index.
    ///
    /// # Panics
    /// Panics if `leaf` is not a leaf node.
    pub fn remove(&mut self, leaf: NodeIdx) {
        assert!(self.nodes[leaf].is_leaf(), "remove: node is not a leaf");
        let parent = self.nodes[leaf].parent;
        match parent {
            None => {
                debug_assert_eq!(self.root, Some(leaf));
                self.root = None;
                self.free_node(leaf);
            }
            Some(parent_idx) => {
                let grandparent = self.nodes[parent_idx].parent;
                let sibling = if self.nodes[parent_idx].left == Some(leaf) {
                    self.nodes[parent_idx]
                        .right
                        .expect("internal node has no right child")
                } else {
                    self.nodes[parent_idx]
                        .left
                        .expect("internal node has no left child")
                };
                self.nodes[sibling].parent = grandparent;
                match grandparent {
                    None => {
                        self.root = Some(sibling);
                    }
                    Some(gp) => {
                        if self.nodes[gp].left == Some(parent_idx) {
                            self.nodes[gp].left = Some(sibling);
                        } else {
                            self.nodes[gp].right = Some(sibling);
                        }
                        self.refit(gp);
                    }
                }
                self.free_node(parent_idx);
                self.free_node(leaf);
            }
        }
    }
    /// Update an existing leaf's AABB.
    ///
    /// If `new_aabb` still fits within the leaf's current fat AABB, the tree
    /// is unchanged and `false` is returned. Otherwise the leaf is removed and
    /// reinserted with a fresh fat AABB, and `true` is returned.
    pub fn update(&mut self, leaf: NodeIdx, new_aabb: BvhAabb) -> bool {
        assert!(self.nodes[leaf].is_leaf(), "update: node is not a leaf");
        if self.nodes[leaf].aabb.contains(&new_aabb) {
            return false;
        }
        let data = self.nodes[leaf].data.expect("leaf must have data");
        self.remove(leaf);
        self.insert(new_aabb, data);
        true
    }
    /// Return all leaf-data pairs `(a, b)` where `a < b` and the leaves
    /// overlap. Uses a stack-based traversal.
    pub fn query_overlapping_pairs(&self) -> Vec<(u32, u32)> {
        let mut pairs = Vec::new();
        let Some(root) = self.root else {
            return pairs;
        };
        let mut leaves: Vec<NodeIdx> = Vec::new();
        let mut stack = vec![root];
        while let Some(idx) = stack.pop() {
            let node = &self.nodes[idx];
            if node.is_leaf() {
                leaves.push(idx);
            } else {
                if let Some(l) = node.left {
                    stack.push(l);
                }
                if let Some(r) = node.right {
                    stack.push(r);
                }
            }
        }
        for (i, &la) in leaves.iter().enumerate() {
            let aabb_a = self.nodes[la].aabb;
            let data_a = match self.nodes[la].data {
                Some(d) => d,
                None => continue,
            };
            for &lb in &leaves[(i + 1)..] {
                if aabb_a.intersects(&self.nodes[lb].aabb) {
                    let Some(data_b) = self.nodes[lb].data else {
                        continue;
                    };
                    let (x, y) = if data_a < data_b {
                        (data_a, data_b)
                    } else {
                        (data_b, data_a)
                    };
                    pairs.push((x, y));
                }
            }
        }
        pairs
    }
    /// Return the data of all leaves whose fat AABB overlaps `aabb`.
    pub fn query_aabb(&self, aabb: &BvhAabb) -> Vec<u32> {
        let mut results = Vec::new();
        if let Some(root) = self.root {
            self.query_aabb_recursive(root, aabb, &mut results);
        }
        results
    }
    /// Return the data of all leaves whose fat AABB contains `point`.
    pub fn query_point(&self, point: Vec3) -> Vec<u32> {
        let pt_aabb = BvhAabb::point(point);
        self.query_aabb(&pt_aabb)
    }
    /// Number of leaf nodes currently in the tree.
    pub fn n_leaves(&self) -> usize {
        self.nodes
            .iter()
            .enumerate()
            .filter(|(idx, n)| n.is_leaf() && !self.free_list.contains(idx))
            .count()
    }
    /// Total number of nodes (leaves + internal) currently in the tree.
    pub fn n_nodes(&self) -> usize {
        self.nodes.len() - self.free_list.len()
    }
    /// Height of the tree (0 for empty, 0 for a single leaf).
    pub fn height(&self) -> i32 {
        match self.root {
            None => 0,
            Some(r) => self.nodes[r].height,
        }
    }
    /// Validate tree invariants. Returns `true` if the tree is consistent.
    ///
    /// Checks:
    /// - root has no parent
    /// - all parent/child pointers are consistent
    /// - internal node AABBs contain both children's AABBs
    pub fn validate(&self) -> bool {
        match self.root {
            None => true,
            Some(root) => {
                if self.nodes[root].parent.is_some() {
                    return false;
                }
                self.validate_node(root)
            }
        }
    }
    /// Find the best sibling for a new leaf using a greedy SAH traversal.
    fn find_best_sibling(&self, leaf: NodeIdx) -> NodeIdx {
        let leaf_sa = self.nodes[leaf].aabb.surface_area();
        let Some(root) = self.root else {
            return leaf;
        };
        let mut best = root;
        let mut best_cost =
            BvhAabb::merge(&self.nodes[root].aabb, &self.nodes[leaf].aabb).surface_area();
        let mut stack: Vec<(NodeIdx, f64)> = vec![(root, 0.0)];
        while let Some((idx, inherited_cost)) = stack.pop() {
            let node = &self.nodes[idx];
            let direct_cost = BvhAabb::merge(&node.aabb, &self.nodes[leaf].aabb).surface_area();
            let cost = direct_cost + inherited_cost;
            if cost < best_cost {
                best_cost = cost;
                best = idx;
            }
            if !node.is_leaf() {
                let child_inherited = inherited_cost + direct_cost - node.aabb.surface_area();
                let lower_bound = leaf_sa + child_inherited;
                if lower_bound < best_cost {
                    if let Some(l) = node.left {
                        stack.push((l, child_inherited));
                    }
                    if let Some(r) = node.right {
                        stack.push((r, child_inherited));
                    }
                }
            }
        }
        best
    }
    /// Walk from `start` to the root, refitting AABBs and updating heights.
    fn refit(&mut self, start: NodeIdx) {
        let mut idx = Some(start);
        while let Some(i) = idx {
            let node = &self.nodes[i];
            if node.is_leaf() {
                idx = node.parent;
                continue;
            }
            let (Some(left), Some(right)) = (node.left, node.right) else {
                idx = node.parent;
                continue;
            };
            let merged = BvhAabb::merge(&self.nodes[left].aabb, &self.nodes[right].aabb);
            let height = 1 + self.nodes[left].height.max(self.nodes[right].height);
            self.nodes[i].aabb = merged;
            self.nodes[i].height = height;
            idx = self.nodes[i].parent;
        }
    }
    fn query_aabb_recursive(&self, idx: NodeIdx, aabb: &BvhAabb, out: &mut Vec<u32>) {
        let node = &self.nodes[idx];
        if !node.aabb.intersects(aabb) {
            return;
        }
        if node.is_leaf() {
            if let Some(d) = node.data {
                out.push(d);
            }
        } else {
            if let Some(l) = node.left {
                self.query_aabb_recursive(l, aabb, out);
            }
            if let Some(r) = node.right {
                self.query_aabb_recursive(r, aabb, out);
            }
        }
    }
    fn validate_node(&self, idx: NodeIdx) -> bool {
        let node = &self.nodes[idx];
        if node.is_leaf() {
            return true;
        }
        let left = match node.left {
            Some(l) => l,
            None => return false,
        };
        let right = match node.right {
            Some(r) => r,
            None => return false,
        };
        if self.nodes[left].parent != Some(idx) {
            return false;
        }
        if self.nodes[right].parent != Some(idx) {
            return false;
        }
        if !node.aabb.contains(&self.nodes[left].aabb) {
            return false;
        }
        if !node.aabb.contains(&self.nodes[right].aabb) {
            return false;
        }
        self.validate_node(left) && self.validate_node(right)
    }
}
impl DynamicBvh {
    /// Collect tree statistics.
    pub fn stats(&self) -> DbvtStats {
        DbvtStats {
            leaf_count: self.n_leaves(),
            node_count: self.n_nodes(),
            height: self.height(),
            root_surface_area: self
                .root
                .map(|r| self.nodes[r].aabb.surface_area())
                .unwrap_or(0.0),
            overlap_pairs: 0,
        }
    }
    /// Refit all internal-node AABBs from leaves upward (full refit pass).
    ///
    /// Useful after bulk updates.
    pub fn refit_all(&mut self) {
        if let Some(root) = self.root {
            self.refit(root);
        }
    }
    /// Ray cast against all leaves; returns data of all leaves hit.
    ///
    /// `origin` + `dir * t` for `t ∈ [t_min, t_max]`.
    pub fn ray_cast(&self, origin: Vec3, dir: Vec3, t_min: f64, t_max: f64) -> Vec<(u32, f64)> {
        let mut hits = Vec::new();
        if let Some(root) = self.root {
            self.ray_cast_recursive(root, origin, dir, t_min, t_max, &mut hits);
        }
        hits.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        hits
    }
    fn ray_cast_recursive(
        &self,
        idx: NodeIdx,
        origin: Vec3,
        dir: Vec3,
        t_min: f64,
        t_max: f64,
        hits: &mut Vec<(u32, f64)>,
    ) {
        let node = &self.nodes[idx];
        if let Some(t) = Self::ray_aabb_slab(origin, dir, node.aabb, t_min, t_max) {
            if node.is_leaf() {
                if let Some(d) = node.data {
                    hits.push((d, t));
                }
            } else {
                if let Some(l) = node.left {
                    self.ray_cast_recursive(l, origin, dir, t_min, t_max, hits);
                }
                if let Some(r) = node.right {
                    self.ray_cast_recursive(r, origin, dir, t_min, t_max, hits);
                }
            }
        }
    }
    /// Ray vs AABB slab intersection; returns entry `t` or `None`.
    fn ray_aabb_slab(
        origin: Vec3,
        dir: Vec3,
        aabb: BvhAabb,
        t_min: f64,
        t_max: f64,
    ) -> Option<f64> {
        let mut tmin = t_min;
        let mut tmax = t_max;
        let orig = [origin.x, origin.y, origin.z];
        let d = [dir.x, dir.y, dir.z];
        let mn = [aabb.min.x, aabb.min.y, aabb.min.z];
        let mx = [aabb.max.x, aabb.max.y, aabb.max.z];
        for i in 0..3 {
            if d[i].abs() < 1e-300 {
                if orig[i] < mn[i] || orig[i] > mx[i] {
                    return None;
                }
            } else {
                let inv = 1.0 / d[i];
                let t1 = (mn[i] - orig[i]) * inv;
                let t2 = (mx[i] - orig[i]) * inv;
                let (ta, tb) = if t1 <= t2 { (t1, t2) } else { (t2, t1) };
                tmin = tmin.max(ta);
                tmax = tmax.min(tb);
                if tmin > tmax {
                    return None;
                }
            }
        }
        if tmax < t_min {
            return None;
        }
        Some(tmin.max(t_min))
    }
    /// Overlap query: return data of all leaves that overlap `query_aabb` and
    /// whose data value satisfies `filter`.
    ///
    /// This is an enhanced version of `query_aabb` that supports a predicate.
    pub fn query_aabb_filtered(&self, aabb: &BvhAabb, filter: impl Fn(u32) -> bool) -> Vec<u32> {
        let mut results = Vec::new();
        if let Some(root) = self.root {
            self.query_aabb_filtered_recursive(root, aabb, &filter, &mut results);
        }
        results
    }
    fn query_aabb_filtered_recursive(
        &self,
        idx: NodeIdx,
        aabb: &BvhAabb,
        filter: &impl Fn(u32) -> bool,
        out: &mut Vec<u32>,
    ) {
        let node = &self.nodes[idx];
        if !node.aabb.intersects(aabb) {
            return;
        }
        if node.is_leaf() {
            if let Some(d) = node.data
                && filter(d)
            {
                out.push(d);
            }
        } else {
            if let Some(l) = node.left {
                self.query_aabb_filtered_recursive(l, aabb, filter, out);
            }
            if let Some(r) = node.right {
                self.query_aabb_filtered_recursive(r, aabb, filter, out);
            }
        }
    }
    /// Bulk insert many `(aabb, data)` pairs and return their leaf indices.
    pub fn insert_bulk(&mut self, items: &[(BvhAabb, u32)]) -> Vec<NodeIdx> {
        items
            .iter()
            .map(|(aabb, data)| self.insert(*aabb, *data))
            .collect()
    }
    /// Remove all leaves and reset the tree to empty.
    pub fn clear(&mut self) {
        self.nodes.clear();
        self.free_list.clear();
        self.root = None;
    }
    /// Collect all leaf data values in the tree (arbitrary order).
    pub fn all_leaf_data(&self) -> Vec<u32> {
        self.nodes
            .iter()
            .enumerate()
            .filter(|(idx, n)| n.is_leaf() && !self.free_list.contains(idx))
            .filter_map(|(_, n)| n.data)
            .collect()
    }
    /// Count overlapping pairs and return both the pairs and stats.
    pub fn query_pairs_with_stats(&self) -> (Vec<(u32, u32)>, DbvtStats) {
        let pairs = self.query_overlapping_pairs();
        let mut s = self.stats();
        s.overlap_pairs = pairs.len();
        (pairs, s)
    }
}
impl DynamicBvh {
    /// Compute BVH quality metrics via a single DFS traversal.
    ///
    /// Returns `None` if the tree is empty.
    #[allow(dead_code)]
    pub fn quality_metrics(&self) -> Option<BvhQuality> {
        let root = self.root?;
        let root_sa = self.nodes[root].aabb.surface_area();
        if root_sa <= 0.0 {
            return None;
        }
        let mut sah_cost = 0.0_f64;
        let mut max_depth = 0usize;
        let mut min_leaf_depth = usize::MAX;
        let mut total_leaf_depth = 0usize;
        let mut leaf_count = 0usize;
        let mut internal_count = 0usize;
        let mut stack: Vec<(NodeIdx, usize)> = vec![(root, 0)];
        while let Some((idx, depth)) = stack.pop() {
            let node = &self.nodes[idx];
            if depth > max_depth {
                max_depth = depth;
            }
            if node.is_leaf() {
                leaf_count += 1;
                total_leaf_depth += depth;
                if depth < min_leaf_depth {
                    min_leaf_depth = depth;
                }
            } else {
                internal_count += 1;
                let sa = node.aabb.surface_area();
                sah_cost += sa / root_sa;
                if let Some(l) = node.left {
                    stack.push((l, depth + 1));
                }
                if let Some(r) = node.right {
                    stack.push((r, depth + 1));
                }
            }
        }
        let avg_leaf_depth = if leaf_count > 0 {
            total_leaf_depth as f64 / leaf_count as f64
        } else {
            0.0
        };
        let min_leaf_depth_out = if min_leaf_depth == usize::MAX {
            0
        } else {
            min_leaf_depth
        };
        Some(BvhQuality {
            sah_cost,
            max_depth,
            min_leaf_depth: min_leaf_depth_out,
            avg_leaf_depth,
            leaf_count,
            internal_count,
        })
    }
    /// Return `true` if the SAH cost is above `threshold` (tree quality degraded).
    ///
    /// Useful as a trigger to rebuild the tree from scratch.
    #[allow(dead_code)]
    pub fn quality_degraded(&self, threshold: f64) -> bool {
        match self.quality_metrics() {
            Some(q) => q.sah_cost > threshold,
            None => false,
        }
    }
}
impl DynamicBvh {
    /// Return the data of all leaves potentially visible inside `frustum`.
    ///
    /// Uses a stack-based traversal that prunes subtrees whose AABB is
    /// fully outside any frustum plane.
    #[allow(dead_code)]
    pub fn frustum_query(&self, frustum: &BvhFrustum) -> Vec<u32> {
        let mut results = Vec::new();
        let Some(root) = self.root else {
            return results;
        };
        let mut stack = vec![root];
        while let Some(idx) = stack.pop() {
            let node = &self.nodes[idx];
            if !frustum.intersects_aabb(&node.aabb) {
                continue;
            }
            if node.is_leaf() {
                if let Some(d) = node.data {
                    results.push(d);
                }
            } else {
                if let Some(l) = node.left {
                    stack.push(l);
                }
                if let Some(r) = node.right {
                    stack.push(r);
                }
            }
        }
        results
    }
}
impl DynamicBvh {
    /// Return the `k` closest leaves (by AABB-center distance) to `point`.
    ///
    /// Result is sorted nearest-first.
    #[allow(dead_code)]
    pub fn k_nearest(&self, point: Vec3, k: usize) -> Vec<NearestLeaf> {
        if k == 0 {
            return Vec::new();
        }
        let mut candidates: Vec<NearestLeaf> = self
            .nodes
            .iter()
            .enumerate()
            .filter(|(idx, n)| n.is_leaf() && !self.free_list.contains(idx))
            .filter_map(|(_, n)| {
                let d = n.data?;
                let c = n.aabb.center();
                let dx = c.x - point.x;
                let dy = c.y - point.y;
                let dz = c.z - point.z;
                Some(NearestLeaf {
                    data: d,
                    dist_sq: dx * dx + dy * dy + dz * dz,
                })
            })
            .collect();
        candidates.sort_by(|a, b| {
            a.dist_sq
                .partial_cmp(&b.dist_sq)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        candidates.truncate(k);
        candidates
    }
}
impl DynamicBvh {
    /// Serialize the entire tree into a flat `Vec`f64`.
    ///
    /// Format per node (14 values):
    /// `\[min_x, min_y, min_z, max_x, max_y, max_z,
    ///   parent_or_neg1, left_or_neg1, right_or_neg1,
    ///   data_or_neg1, height, is_leaf, is_on_free_list, PADDING]`
    ///
    /// The first element of the buffer is the total number of nodes (as f64).
    #[allow(dead_code)]
    pub fn serialize(&self) -> Vec<f64> {
        const STRIDE: usize = 14;
        let mut buf = Vec::with_capacity(1 + self.nodes.len() * STRIDE);
        buf.push(self.nodes.len() as f64);
        for (idx, node) in self.nodes.iter().enumerate() {
            buf.push(node.aabb.min.x);
            buf.push(node.aabb.min.y);
            buf.push(node.aabb.min.z);
            buf.push(node.aabb.max.x);
            buf.push(node.aabb.max.y);
            buf.push(node.aabb.max.z);
            buf.push(node.parent.map(|p| p as f64).unwrap_or(-1.0));
            buf.push(node.left.map(|l| l as f64).unwrap_or(-1.0));
            buf.push(node.right.map(|r| r as f64).unwrap_or(-1.0));
            buf.push(node.data.map(|d| d as f64).unwrap_or(-1.0));
            buf.push(node.height as f64);
            buf.push(if node.is_leaf() { 1.0 } else { 0.0 });
            buf.push(if self.free_list.contains(&idx) {
                1.0
            } else {
                0.0
            });
            buf.push(0.0);
        }
        buf
    }
    /// Return the root node index, or `None` if the tree is empty.
    #[allow(dead_code)]
    pub fn root(&self) -> Option<NodeIdx> {
        self.root
    }
}
impl DynamicBvh {
    /// Update leaf `leaf` using a new tight AABB expanded in the direction of
    /// `velocity * dt` (speculative / predicted motion).
    ///
    /// The resulting fat AABB covers both the current and predicted positions.
    /// Returns `true` if the tree was modified (re-insertion occurred).
    #[allow(dead_code)]
    pub fn move_proxy(
        &mut self,
        leaf: NodeIdx,
        new_aabb: BvhAabb,
        velocity: Vec3,
        dt: f64,
    ) -> bool {
        assert!(self.nodes[leaf].is_leaf(), "move_proxy: node is not a leaf");
        let displacement = Vec3::new(velocity.x * dt, velocity.y * dt, velocity.z * dt);
        let predicted = BvhAabb {
            min: Vec3::new(
                new_aabb.min.x + displacement.x.min(0.0),
                new_aabb.min.y + displacement.y.min(0.0),
                new_aabb.min.z + displacement.z.min(0.0),
            ),
            max: Vec3::new(
                new_aabb.max.x + displacement.x.max(0.0),
                new_aabb.max.y + displacement.y.max(0.0),
                new_aabb.max.z + displacement.z.max(0.0),
            ),
        };
        let swept = BvhAabb::merge(&new_aabb, &predicted);
        if self.nodes[leaf].aabb.contains(&swept) {
            return false;
        }
        let data = self.nodes[leaf].data.expect("leaf must have data");
        self.remove(leaf);
        self.insert(swept, data);
        true
    }
}
impl DynamicBvh {
    /// Attempt one SAH-improving rotation at node `idx`.
    ///
    /// Tries swapping each grandchild with its uncle; keeps the swap
    /// only if it strictly reduces the total surface area of the two
    /// affected internal nodes.
    ///
    /// Returns `true` if a rotation was performed.
    #[allow(dead_code)]
    pub fn try_rotate(&mut self, idx: NodeIdx) -> bool {
        let node = &self.nodes[idx];
        if node.is_leaf() {
            return false;
        }
        let left = match node.left {
            Some(l) => l,
            None => return false,
        };
        let right = match node.right {
            Some(r) => r,
            None => return false,
        };
        let current_cost =
            self.nodes[left].aabb.surface_area() + self.nodes[right].aabb.surface_area();
        let mut best_cost = current_cost;
        let mut best_rotation: Option<(NodeIdx, NodeIdx, NodeIdx)> = None;
        if !self.nodes[right].is_leaf() {
            let (Some(rl), Some(rr)) = (self.nodes[right].left, self.nodes[right].right) else {
                return false;
            };
            let new_right_sa =
                BvhAabb::merge(&self.nodes[left].aabb, &self.nodes[rr].aabb).surface_area();
            let candidate_cost = new_right_sa + self.nodes[left].aabb.surface_area();
            if candidate_cost < best_cost {
                best_cost = candidate_cost;
                best_rotation = Some((idx, left, rl));
            }
            let new_right_sa2 =
                BvhAabb::merge(&self.nodes[left].aabb, &self.nodes[rl].aabb).surface_area();
            let candidate_cost2 = new_right_sa2 + self.nodes[left].aabb.surface_area();
            if candidate_cost2 < best_cost {
                best_cost = candidate_cost2;
                best_rotation = Some((idx, left, rr));
            }
        }
        if !self.nodes[left].is_leaf() {
            let (Some(ll), Some(lr)) = (self.nodes[left].left, self.nodes[left].right) else {
                return false;
            };
            let new_left_sa =
                BvhAabb::merge(&self.nodes[right].aabb, &self.nodes[lr].aabb).surface_area();
            let candidate_cost = new_left_sa + self.nodes[right].aabb.surface_area();
            if candidate_cost < best_cost {
                best_cost = candidate_cost;
                best_rotation = Some((idx, right, ll));
            }
            let new_left_sa2 =
                BvhAabb::merge(&self.nodes[right].aabb, &self.nodes[ll].aabb).surface_area();
            let candidate_cost2 = new_left_sa2 + self.nodes[right].aabb.surface_area();
            if candidate_cost2 < best_cost {
                let _ = best_cost;
                best_rotation = Some((idx, right, lr));
            }
        }
        if let Some((parent_idx, child_in_parent, grandchild_to_swap)) = best_rotation {
            let Some(grandchild_parent) = self.nodes[grandchild_to_swap].parent else {
                return false;
            };
            let gp_left = self.nodes[grandchild_parent].left;
            let gp_right = self.nodes[grandchild_parent].right;
            let pl = self.nodes[parent_idx].left;
            let pr = self.nodes[parent_idx].right;
            if pl == Some(child_in_parent) {
                self.nodes[parent_idx].left = Some(grandchild_to_swap);
            } else if pr == Some(child_in_parent) {
                self.nodes[parent_idx].right = Some(grandchild_to_swap);
            }
            if gp_left == Some(grandchild_to_swap) {
                self.nodes[grandchild_parent].left = Some(child_in_parent);
            } else if gp_right == Some(grandchild_to_swap) {
                self.nodes[grandchild_parent].right = Some(child_in_parent);
            }
            self.nodes[grandchild_to_swap].parent = Some(parent_idx);
            self.nodes[child_in_parent].parent = Some(grandchild_parent);
            self.refit(grandchild_parent);
            self.refit(parent_idx);
            true
        } else {
            false
        }
    }
    /// Perform one bottom-up pass of SAH-improving rotations across the tree.
    ///
    /// Visits every non-leaf node from leaves upward and attempts `try_rotate`.
    #[allow(dead_code)]
    pub fn optimize_rotations(&mut self) {
        let internal: Vec<NodeIdx> = self
            .nodes
            .iter()
            .enumerate()
            .filter(|(idx, n)| !n.is_leaf() && !self.free_list.contains(idx))
            .map(|(i, _)| i)
            .collect();
        let mut ordered = internal;
        ordered.sort_by_key(|&i| self.nodes[i].height);
        for idx in ordered {
            if !self.nodes[idx].is_leaf() {
                self.try_rotate(idx);
            }
        }
    }
}
impl DynamicBvh {
    /// Refit only the ancestors of a specific leaf node.
    ///
    /// More efficient than `refit_all` when only one leaf changed.
    #[allow(dead_code)]
    pub fn refit_leaf_ancestors(&mut self, leaf: NodeIdx) {
        assert!(
            self.nodes[leaf].is_leaf(),
            "refit_leaf_ancestors: not a leaf"
        );
        if let Some(parent) = self.nodes[leaf].parent {
            self.refit(parent);
        }
    }
    /// Update a leaf's AABB *in place* without re-inserting it, then refit
    /// its ancestor chain.
    ///
    /// **Only use this when `new_aabb` is contained within the current fat AABB**;
    /// otherwise the tree's AABBs may become incorrect.  For unconstrained moves
    /// use [`DynamicBvh::update`] or [`DynamicBvh::move_proxy`].
    #[allow(dead_code)]
    pub fn update_leaf_inplace(&mut self, leaf: NodeIdx, new_aabb: BvhAabb) {
        assert!(
            self.nodes[leaf].is_leaf(),
            "update_leaf_inplace: not a leaf"
        );
        self.nodes[leaf].aabb = new_aabb;
        self.refit_leaf_ancestors(leaf);
    }
    /// Compute the total SAH cost of the tree as a single scalar.
    ///
    /// Defined as `sum(SA(internal_node)) / SA(root)`.
    /// Lower is better; returns `0.0` for trees with 0 or 1 nodes.
    #[allow(dead_code)]
    pub fn sah_cost(&self) -> f64 {
        let Some(root) = self.root else {
            return 0.0;
        };
        let root_sa = self.nodes[root].aabb.surface_area();
        if root_sa <= 0.0 {
            return 0.0;
        }
        let mut cost = 0.0_f64;
        let mut stack = vec![root];
        while let Some(idx) = stack.pop() {
            let node = &self.nodes[idx];
            if !node.is_leaf() {
                cost += node.aabb.surface_area() / root_sa;
                if let Some(l) = node.left {
                    stack.push(l);
                }
                if let Some(r) = node.right {
                    stack.push(r);
                }
            }
        }
        cost
    }
}
impl DynamicBvh {
    /// Compute the Surface Area Heuristic (SAH) cost of the entire tree.
    ///
    /// This is an alias for [`DynamicBvh::sah_cost`] provided under the name
    /// requested by the algorithm expansion spec so that callers can use either
    /// name without changing existing code.
    ///
    /// Returns `0.0` for an empty or degenerate tree.
    #[allow(dead_code)]
    pub fn compute_sah_cost(&self) -> f64 {
        self.sah_cost()
    }
    /// Apply one round of balance-improving rotations to the entire tree.
    ///
    /// Walks every internal node in bottom-up order and attempts a single
    /// SAH-improving swap (see `try_rotate`).  Nodes are visited from
    /// shallowest to deepest so that parent improvements propagate before
    /// children are processed.
    ///
    /// Returns the number of rotations that were accepted.
    #[allow(dead_code)]
    pub fn balance_rotation(&mut self) -> usize {
        let internal: Vec<NodeIdx> = self
            .nodes
            .iter()
            .enumerate()
            .filter(|(idx, n)| !n.is_leaf() && !self.free_list.contains(idx))
            .map(|(i, _)| i)
            .collect();
        let mut ordered = internal;
        ordered.sort_by_key(|&i| self.nodes[i].height);
        let mut accepted = 0usize;
        for idx in ordered {
            if !self.nodes[idx].is_leaf() && self.try_rotate(idx) {
                accepted += 1;
            }
        }
        accepted
    }
    /// Frustum-culling traversal: return data of all leaves whose fat AABB is
    /// potentially inside `frustum`.
    ///
    /// This is a named alias for `frustum_query` following the algorithm
    /// expansion spec.  Uses a stack-based traversal that prunes fully-outside
    /// subtrees for O(k + log n) complexity where k is the number of hits.
    #[allow(dead_code)]
    pub fn traverse_frustum(&self, frustum: &BvhFrustum) -> Vec<u32> {
        self.frustum_query(frustum)
    }
}
impl DynamicBvh {
    /// Return the depth of a given node in the tree (root = 0).
    ///
    /// Returns `None` if `idx` is on the free list.
    #[allow(dead_code)]
    pub fn node_depth(&self, idx: NodeIdx) -> Option<usize> {
        if self.free_list.contains(&idx) {
            return None;
        }
        let mut depth = 0usize;
        let mut cur = self.nodes[idx].parent;
        while let Some(p) = cur {
            depth += 1;
            cur = self.nodes[p].parent;
        }
        Some(depth)
    }
    /// Return the leaf node index whose `data` matches `target`, or `None`.
    #[allow(dead_code)]
    pub fn find_leaf(&self, target: u32) -> Option<NodeIdx> {
        let root = self.root?;
        self.find_leaf_recursive(root, target)
    }
    fn find_leaf_recursive(&self, idx: NodeIdx, target: u32) -> Option<NodeIdx> {
        let node = &self.nodes[idx];
        if node.is_leaf() {
            return if node.data == Some(target) {
                Some(idx)
            } else {
                None
            };
        }
        if let Some(l) = node.left
            && let Some(found) = self.find_leaf_recursive(l, target)
        {
            return Some(found);
        }
        if let Some(r) = node.right
            && let Some(found) = self.find_leaf_recursive(r, target)
        {
            return Some(found);
        }
        None
    }
    /// Collect all leaf data in DFS pre-order (left-first).
    #[allow(dead_code)]
    pub fn leaf_data_preorder(&self) -> Vec<u32> {
        let Some(root) = self.root else {
            return Vec::new();
        };
        let mut result = Vec::new();
        self.leaf_data_preorder_recursive(root, &mut result);
        result
    }
    fn leaf_data_preorder_recursive(&self, idx: NodeIdx, out: &mut Vec<u32>) {
        let node = &self.nodes[idx];
        if node.is_leaf() {
            if let Some(d) = node.data {
                out.push(d);
            }
            return;
        }
        if let Some(l) = node.left {
            self.leaf_data_preorder_recursive(l, out);
        }
        if let Some(r) = node.right {
            self.leaf_data_preorder_recursive(r, out);
        }
    }
    /// Return the AABB of the root node, or `None` if the tree is empty.
    #[allow(dead_code)]
    pub fn root_aabb(&self) -> Option<BvhAabb> {
        self.root.map(|r| self.nodes[r].aabb)
    }
    /// Compute the maximum tree depth by traversal.
    #[allow(dead_code)]
    pub fn max_depth(&self) -> usize {
        let Some(root) = self.root else {
            return 0;
        };
        let mut max_d = 0usize;
        let mut stack: Vec<(NodeIdx, usize)> = vec![(root, 0)];
        while let Some((idx, d)) = stack.pop() {
            if d > max_d {
                max_d = d;
            }
            let node = &self.nodes[idx];
            if !node.is_leaf() {
                if let Some(l) = node.left {
                    stack.push((l, d + 1));
                }
                if let Some(r) = node.right {
                    stack.push((r, d + 1));
                }
            }
        }
        max_d
    }
    /// Return the number of internal nodes (non-leaf, non-free).
    #[allow(dead_code)]
    pub fn n_internal(&self) -> usize {
        self.nodes
            .iter()
            .enumerate()
            .filter(|(idx, n)| !n.is_leaf() && !self.free_list.contains(idx))
            .count()
    }
}
impl DynamicBvh {
    /// Return data of all leaves whose fat AABB overlaps a sphere defined by
    /// `center` and `radius`.
    ///
    /// Uses the AABB-vs-sphere test: the sphere overlaps an AABB iff the
    /// squared distance from `center` to the nearest point on the AABB is ≤ `radius²`.
    #[allow(dead_code)]
    pub fn query_sphere(&self, center: Vec3, radius: f64) -> Vec<u32> {
        let radius_sq = radius * radius;
        let mut results = Vec::new();
        if let Some(root) = self.root {
            self.query_sphere_recursive(root, center, radius_sq, &mut results);
        }
        results
    }
    fn query_sphere_recursive(
        &self,
        idx: NodeIdx,
        center: Vec3,
        radius_sq: f64,
        out: &mut Vec<u32>,
    ) {
        let node = &self.nodes[idx];
        if node.aabb.point_dist_sq(center) > radius_sq {
            return;
        }
        if node.is_leaf() {
            if let Some(d) = node.data {
                out.push(d);
            }
        } else {
            if let Some(l) = node.left {
                self.query_sphere_recursive(l, center, radius_sq, out);
            }
            if let Some(r) = node.right {
                self.query_sphere_recursive(r, center, radius_sq, out);
            }
        }
    }
    /// Return `true` if any leaf's fat AABB overlaps the given sphere.
    #[allow(dead_code)]
    pub fn any_in_sphere(&self, center: Vec3, radius: f64) -> bool {
        !self.query_sphere(center, radius).is_empty()
    }
}
impl DynamicBvh {
    /// Return data of all leaves whose fat AABB is hit by a capsule query.
    ///
    /// Implemented as AABB-vs-segment-distance test: an AABB intersects the
    /// capsule iff the minimum distance from the AABB to the line segment is
    /// ≤ `capsule.radius`.
    #[allow(dead_code)]
    pub fn query_capsule(&self, capsule: BvhCapsule) -> Vec<u32> {
        let mut results = Vec::new();
        if let Some(root) = self.root {
            self.query_capsule_recursive(root, &capsule, &mut results);
        }
        results
    }
    fn query_capsule_recursive(&self, idx: NodeIdx, capsule: &BvhCapsule, out: &mut Vec<u32>) {
        let node = &self.nodes[idx];
        if self.aabb_segment_dist_sq(&node.aabb, capsule.start, capsule.end)
            > capsule.radius * capsule.radius
        {
            return;
        }
        if node.is_leaf() {
            if let Some(d) = node.data {
                out.push(d);
            }
        } else {
            if let Some(l) = node.left {
                self.query_capsule_recursive(l, capsule, out);
            }
            if let Some(r) = node.right {
                self.query_capsule_recursive(r, capsule, out);
            }
        }
    }
    /// Squared distance from the closest point on the AABB to the segment `\[p, q\]`.
    fn aabb_segment_dist_sq(&self, aabb: &BvhAabb, p: Vec3, q: Vec3) -> f64 {
        let d = q - p;
        let len_sq = d.x * d.x + d.y * d.y + d.z * d.z;
        if len_sq < 1e-300 {
            return aabb.point_dist_sq(p);
        }
        let c = aabb.center();
        let t = ((c.x - p.x) * d.x + (c.y - p.y) * d.y + (c.z - p.z) * d.z) / len_sq;
        let t = t.clamp(0.0, 1.0);
        let closest = Vec3::new(p.x + t * d.x, p.y + t * d.y, p.z + t * d.z);
        aabb.point_dist_sq(closest)
    }
}
impl DynamicBvh {
    /// Rebuild the tree from scratch using the current fat AABBs of all leaves.
    ///
    /// Removes all nodes, collects leaf data and AABBs, then re-inserts them.
    /// This is a full O(n log n) rebuild; use only when the tree is severely
    /// degraded (e.g. after many updates without rotations).
    #[allow(dead_code)]
    pub fn rebuild(&mut self) {
        let leaves: Vec<(BvhAabb, u32)> = self
            .nodes
            .iter()
            .enumerate()
            .filter(|(idx, n)| n.is_leaf() && !self.free_list.contains(idx))
            .filter_map(|(_, n)| n.data.map(|d| (n.aabb, d)))
            .collect();
        self.clear();
        for (aabb, data) in leaves {
            let tight = BvhAabb {
                min: Vec3::new(
                    aabb.min.x + FAT_MARGIN,
                    aabb.min.y + FAT_MARGIN,
                    aabb.min.z + FAT_MARGIN,
                ),
                max: Vec3::new(
                    aabb.max.x - FAT_MARGIN,
                    aabb.max.y - FAT_MARGIN,
                    aabb.max.z - FAT_MARGIN,
                ),
            };
            self.insert(tight, data);
        }
    }
    /// Return a snapshot of all leaf `(data, fat_aabb)` pairs without
    /// modifying the tree.
    #[allow(dead_code)]
    pub fn leaf_snapshot(&self) -> Vec<(u32, BvhAabb)> {
        self.nodes
            .iter()
            .enumerate()
            .filter(|(idx, n)| n.is_leaf() && !self.free_list.contains(idx))
            .filter_map(|(_, n)| n.data.map(|d| (d, n.aabb)))
            .collect()
    }
    /// Return the fat AABB for a given leaf `idx`, or `None` if it is not
    /// a valid leaf.
    #[allow(dead_code)]
    pub fn leaf_aabb(&self, idx: NodeIdx) -> Option<BvhAabb> {
        if idx < self.nodes.len() && self.nodes[idx].is_leaf() && !self.free_list.contains(&idx) {
            Some(self.nodes[idx].aabb)
        } else {
            None
        }
    }
    /// Remove all leaves whose data value satisfies `pred`.
    ///
    /// Returns the number of leaves removed.
    #[allow(dead_code)]
    pub fn remove_where(&mut self, pred: impl Fn(u32) -> bool) -> usize {
        let to_remove: Vec<NodeIdx> = self
            .nodes
            .iter()
            .enumerate()
            .filter(|(idx, n)| {
                n.is_leaf() && !self.free_list.contains(idx) && n.data.map(&pred).unwrap_or(false)
            })
            .map(|(i, _)| i)
            .collect();
        let count = to_remove.len();
        for idx in to_remove {
            self.remove(idx);
        }
        count
    }
    /// Replace the data value of a leaf without changing its AABB.
    ///
    /// Returns `true` if the leaf was found and updated.
    #[allow(dead_code)]
    pub fn relabel_leaf(&mut self, leaf: NodeIdx, new_data: u32) -> bool {
        if leaf < self.nodes.len() && self.nodes[leaf].is_leaf() && !self.free_list.contains(&leaf)
        {
            self.nodes[leaf].data = Some(new_data);
            true
        } else {
            false
        }
    }
}
impl DynamicBvh {
    /// Return the number of leaves in the subtree rooted at `idx`.
    #[allow(dead_code)]
    pub fn subtree_leaf_count(&self, idx: NodeIdx) -> usize {
        let node = &self.nodes[idx];
        if node.is_leaf() {
            return 1;
        }
        let left_count = node.left.map_or(0, |l| self.subtree_leaf_count(l));
        let right_count = node.right.map_or(0, |r| self.subtree_leaf_count(r));
        left_count + right_count
    }
    /// Return `true` if the heights stored in nodes are consistent with the
    /// actual tree structure.
    #[allow(dead_code)]
    pub fn validate_heights(&self) -> bool {
        let Some(root) = self.root else {
            return true;
        };
        self.validate_heights_recursive(root)
    }
    fn validate_heights_recursive(&self, idx: NodeIdx) -> bool {
        let node = &self.nodes[idx];
        if node.is_leaf() {
            return node.height == 0;
        }
        let left = match node.left {
            Some(l) => l,
            None => return false,
        };
        let right = match node.right {
            Some(r) => r,
            None => return false,
        };
        let expected = 1 + self.nodes[left].height.max(self.nodes[right].height);
        if node.height != expected {
            return false;
        }
        self.validate_heights_recursive(left) && self.validate_heights_recursive(right)
    }
}
/// Quality metrics for the BVH tree.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct BvhQuality {
    /// SAH cost of the tree: sum of (SA(node) / SA(root)) for all internal nodes.
    pub sah_cost: f64,
    /// Maximum depth (root = depth 0).
    pub max_depth: usize,
    /// Minimum depth of any leaf.
    pub min_leaf_depth: usize,
    /// Average depth of leaves.
    pub avg_leaf_depth: f64,
    /// Number of leaves.
    pub leaf_count: usize,
    /// Number of internal nodes.
    pub internal_count: usize,
}
