//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::math::Vec3;

use super::functions::morton_encode;

pub(super) enum RTreeNode<T: Clone> {
    Leaf {
        items: Vec<(SpatialAabb, T)>,
        bbox: SpatialAabb,
    },
    Internal {
        children: Vec<RTreeNode<T>>,
        bbox: SpatialAabb,
    },
}
/// KD-tree for k-nearest-neighbor and range searches.
pub struct KdTree {
    pub(super) nodes: Vec<KdNode>,
    pub(super) points: Vec<Vec3>,
}
impl KdTree {
    /// Build a balanced KD-tree from a list of points.
    ///
    /// The resulting tree uses median splitting along alternating axes.
    pub fn build(points: Vec<Vec3>) -> Self {
        let n = points.len();
        let mut tree = KdTree {
            nodes: Vec::with_capacity(n),
            points,
        };
        if n > 0 {
            let mut indices: Vec<usize> = (0..n).collect();
            tree.build_recursive(&mut indices, 0);
        }
        tree
    }
    fn build_recursive(&mut self, indices: &mut [usize], depth: usize) -> usize {
        let axis = (depth % 3) as u8;
        let ax = axis as usize;
        indices.sort_unstable_by(|&a, &b| {
            self.points[a][ax]
                .partial_cmp(&self.points[b][ax])
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        let mid = indices.len() / 2;
        let point_idx = indices[mid];
        let node_idx = self.nodes.len();
        self.nodes.push(KdNode {
            point_idx,
            left: None,
            right: None,
            axis,
        });
        if mid > 0 {
            let left_idx = self.build_recursive(&mut indices[..mid], depth + 1);
            self.nodes[node_idx].left = Some(left_idx);
        }
        if mid + 1 < indices.len() {
            let right_idx = self.build_recursive(&mut indices[mid + 1..], depth + 1);
            self.nodes[node_idx].right = Some(right_idx);
        }
        node_idx
    }
    /// Return `(point_index, distance²)` of the single nearest point to `query`.
    pub fn nearest(&self, query: Vec3) -> Option<(usize, f64)> {
        if self.nodes.is_empty() {
            return None;
        }
        let mut best = (0usize, f64::INFINITY);
        self.kd_nearest(0, query, &mut best);
        Some(best)
    }
    fn kd_nearest(&self, node_idx: usize, query: Vec3, best: &mut (usize, f64)) {
        let node = &self.nodes[node_idx];
        let p = self.points[node.point_idx];
        let d2 = (p - query).norm_squared();
        if d2 < best.1 {
            *best = (node.point_idx, d2);
        }
        let ax = node.axis as usize;
        let diff = query[ax] - p[ax];
        let (near, far) = if diff <= 0.0 {
            (node.left, node.right)
        } else {
            (node.right, node.left)
        };
        if let Some(n) = near {
            self.kd_nearest(n, query, best);
        }
        if diff * diff < best.1
            && let Some(f) = far
        {
            self.kd_nearest(f, query, best);
        }
    }
    /// Return the `k` nearest points as `(point_index, distance²)` sorted by
    /// ascending distance.
    pub fn k_nearest(&self, query: Vec3, k: usize) -> Vec<(usize, f64)> {
        if self.nodes.is_empty() || k == 0 {
            return Vec::new();
        }
        let mut heap: Vec<(usize, f64)> = Vec::with_capacity(k + 1);
        self.kd_k_nearest(0, query, k, &mut heap);
        heap.sort_unstable_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        heap.truncate(k);
        heap
    }
    fn kd_k_nearest(&self, node_idx: usize, query: Vec3, k: usize, heap: &mut Vec<(usize, f64)>) {
        let node = &self.nodes[node_idx];
        let p = self.points[node.point_idx];
        let d2 = (p - query).norm_squared();
        let worst = heap.iter().map(|x| x.1).fold(f64::NEG_INFINITY, f64::max);
        if heap.len() < k || d2 < worst {
            heap.push((node.point_idx, d2));
            if heap.len() > k {
                let worst_idx = heap
                    .iter()
                    .enumerate()
                    .max_by(|a, b| {
                        a.1.1
                            .partial_cmp(&b.1.1)
                            .unwrap_or(std::cmp::Ordering::Equal)
                    })
                    .map(|(i, _)| i)
                    .expect("heap is non-empty");
                heap.swap_remove(worst_idx);
            }
        }
        let ax = node.axis as usize;
        let diff = query[ax] - p[ax];
        let (near, far) = if diff <= 0.0 {
            (node.left, node.right)
        } else {
            (node.right, node.left)
        };
        if let Some(n) = near {
            self.kd_k_nearest(n, query, k, heap);
        }
        let current_worst = heap.iter().map(|x| x.1).fold(f64::NEG_INFINITY, f64::max);
        if (heap.len() < k || diff * diff < current_worst)
            && let Some(f) = far
        {
            self.kd_k_nearest(f, query, k, heap);
        }
    }
    /// Return the indices of all points within `radius` of `query`.
    pub fn range_search(&self, query: Vec3, radius: f64) -> Vec<usize> {
        if self.nodes.is_empty() {
            return Vec::new();
        }
        let mut results = Vec::new();
        let r2 = radius * radius;
        self.kd_range(0, query, r2, &mut results);
        results
    }
    fn kd_range(&self, node_idx: usize, query: Vec3, r2: f64, out: &mut Vec<usize>) {
        let node = &self.nodes[node_idx];
        let p = self.points[node.point_idx];
        if (p - query).norm_squared() <= r2 {
            out.push(node.point_idx);
        }
        let ax = node.axis as usize;
        let diff = query[ax] - p[ax];
        if let Some(n) = node.left
            && (diff <= 0.0 || diff * diff <= r2)
        {
            self.kd_range(n, query, r2, out);
        }
        if let Some(f) = node.right
            && (diff >= 0.0 || diff * diff <= r2)
        {
            self.kd_range(f, query, r2, out);
        }
    }
    /// Number of points in the tree.
    pub fn len(&self) -> usize {
        self.points.len()
    }
    /// Return `true` when the tree contains no points.
    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }
}
/// Statistics computed from a spatial index.
#[derive(Debug, Clone, Default)]
pub struct SpatialIndexStats {
    /// Total number of items.
    pub total_items: usize,
    /// Maximum depth of the tree.
    pub max_depth: usize,
    /// Number of leaf nodes.
    pub n_leaves: usize,
    /// Number of internal nodes.
    pub n_internal: usize,
    /// Average items per leaf.
    pub avg_items_per_leaf: f64,
}
/// A simple 1-D range tree (sorted array with binary search).
///
/// Allows efficient range queries: all values in `[lo, hi]`.
pub struct RangeTree1D {
    pub(super) sorted: Vec<(f64, usize)>,
}
impl RangeTree1D {
    /// Build from a list of 1-D values.
    pub fn build(values: &[f64]) -> Self {
        let mut sorted: Vec<(f64, usize)> = values
            .iter()
            .copied()
            .enumerate()
            .map(|(i, v)| (v, i))
            .collect();
        sorted.sort_unstable_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        Self { sorted }
    }
    /// Return the original indices of all values in `[lo, hi]`.
    pub fn range_query(&self, lo: f64, hi: f64) -> Vec<usize> {
        let start = self.sorted.partition_point(|&(v, _)| v < lo);
        let end = self.sorted.partition_point(|&(v, _)| v <= hi);
        self.sorted[start..end].iter().map(|&(_, i)| i).collect()
    }
    /// Return the number of stored values.
    pub fn len(&self) -> usize {
        self.sorted.len()
    }
    /// Return `true` if empty.
    pub fn is_empty(&self) -> bool {
        self.sorted.is_empty()
    }
    /// Return the minimum value (or `None` if empty).
    pub fn min_value(&self) -> Option<f64> {
        self.sorted.first().map(|&(v, _)| v)
    }
    /// Return the maximum value (or `None` if empty).
    pub fn max_value(&self) -> Option<f64> {
        self.sorted.last().map(|&(v, _)| v)
    }
}
/// A simplified R-tree node for bulk-loaded spatial indexing.
///
/// Uses a single-level structure: the root contains a list of leaf entries,
/// each with an AABB and a value.
pub struct FlatRTree<T: Clone> {
    pub(super) entries: Vec<(SpatialAabb, T)>,
}
impl<T: Clone> FlatRTree<T> {
    /// Create an empty R-tree.
    pub fn new(_page_size: usize) -> Self {
        Self {
            entries: Vec::new(),
        }
    }
    /// Insert an entry with bounding box `aabb` and value `val`.
    pub fn insert(&mut self, aabb: SpatialAabb, val: T) {
        self.entries.push((aabb, val));
    }
    /// Query all entries whose AABBs overlap with `query`.
    pub fn query_overlap(&self, query: &SpatialAabb) -> Vec<&T> {
        self.entries
            .iter()
            .filter(|(bb, _)| bb.intersects(query))
            .map(|(_, v)| v)
            .collect()
    }
    /// Query all entries whose AABBs contain point `p`.
    pub fn query_point(&self, p: Vec3) -> Vec<&T> {
        let pt_aabb = SpatialAabb::new(p, p);
        self.query_overlap(&pt_aabb)
    }
    /// Number of entries.
    pub fn len(&self) -> usize {
        self.entries.len()
    }
    /// Return `true` if empty.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
    /// Nearest entry to point `p` (linear scan fallback for simplified R-tree).
    pub fn nearest(&self, p: Vec3) -> Option<&T> {
        self.entries
            .iter()
            .min_by(|(a, _), (b, _)| {
                let da = a.sq_dist_to_point(p);
                let db = b.sq_dist_to_point(p);
                da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|(_, v)| v)
    }
}
/// A ball tree (metric tree) for nearest-neighbor search in 3D.
///
/// Each node defines a ball (center + radius) that bounds all points in its
/// subtree. Nodes are organized as a binary tree with median-split along the
/// principal axis of the bounding ball.
pub struct BallTree {
    pub(super) nodes: Vec<BallTreeNode>,
    pub(super) points: Vec<Vec3>,
}
impl BallTree {
    /// Build a ball tree from a list of points.
    pub fn build(points: Vec<Vec3>) -> Self {
        let n = points.len();
        let mut tree = BallTree {
            nodes: Vec::with_capacity(n * 2),
            points,
        };
        if n > 0 {
            let mut indices: Vec<usize> = (0..n).collect();
            tree.build_node(&mut indices);
        }
        tree
    }
    fn build_node(&mut self, indices: &mut [usize]) -> usize {
        let centroid = Self::centroid_of(&self.points, indices);
        let radius = indices
            .iter()
            .map(|&i| (self.points[i] - centroid).norm())
            .fold(0.0f64, f64::max);
        if indices.len() <= 8 {
            let node_idx = self.nodes.len();
            self.nodes.push(BallTreeNode {
                pivot: indices[indices.len() / 2],
                radius,
                left: None,
                right: None,
                items: indices.to_vec(),
            });
            return node_idx;
        }
        let (axis, _) = Self::max_spread_axis(&self.points, indices);
        indices.sort_unstable_by(|&a, &b| {
            self.points[a][axis]
                .partial_cmp(&self.points[b][axis])
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        let mid = indices.len() / 2;
        let pivot = indices[mid];
        let node_idx = self.nodes.len();
        self.nodes.push(BallTreeNode {
            pivot,
            radius,
            left: None,
            right: None,
            items: Vec::new(),
        });
        let left_idx = self.build_node(&mut indices[..mid]);
        self.nodes[node_idx].left = Some(left_idx);
        if mid + 1 < indices.len() {
            let right_idx = self.build_node(&mut indices[mid + 1..]);
            self.nodes[node_idx].right = Some(right_idx);
        }
        node_idx
    }
    fn centroid_of(points: &[Vec3], indices: &[usize]) -> Vec3 {
        if indices.is_empty() {
            return Vec3::zeros();
        }
        let sum: Vec3 = indices
            .iter()
            .map(|&i| points[i])
            .fold(Vec3::zeros(), |a, b| a + b);
        sum / indices.len() as f64
    }
    fn max_spread_axis(points: &[Vec3], indices: &[usize]) -> (usize, f64) {
        let mut best_axis = 0;
        let mut best_spread = 0.0f64;
        for (axis, _) in [0usize; 3].iter().enumerate() {
            let vals: Vec<f64> = indices.iter().map(|&i| points[i][axis]).collect();
            let min = vals.iter().cloned().fold(f64::INFINITY, f64::min);
            let max = vals.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
            let spread = max - min;
            if spread > best_spread {
                best_spread = spread;
                best_axis = axis;
            }
        }
        (best_axis, best_spread)
    }
    /// Find the single nearest point to `query`.
    ///
    /// Returns `(point_index, distance)` or `None` if the tree is empty.
    pub fn nearest(&self, query: Vec3) -> Option<(usize, f64)> {
        if self.nodes.is_empty() {
            return None;
        }
        let mut best = (0usize, f64::INFINITY);
        self.nn_search(0, query, &mut best);
        if best.1 == f64::INFINITY {
            None
        } else {
            Some((best.0, best.1.sqrt()))
        }
    }
    fn nn_search(&self, node_idx: usize, query: Vec3, best: &mut (usize, f64)) {
        if node_idx >= self.nodes.len() {
            return;
        }
        let node = &self.nodes[node_idx];
        let pivot_d2 = (self.points[node.pivot] - query).norm_squared();
        if pivot_d2 < best.1 {
            *best = (node.pivot, pivot_d2);
        }
        for &i in &node.items {
            let d2 = (self.points[i] - query).norm_squared();
            if d2 < best.1 {
                *best = (i, d2);
            }
        }
        let ball_center = self.points[node.pivot];
        let dist_to_ball = ((ball_center - query).norm() - node.radius).max(0.0);
        if dist_to_ball * dist_to_ball >= best.1 {
            return;
        }
        if let Some(l) = node.left {
            self.nn_search(l, query, best);
        }
        if let Some(r) = node.right {
            self.nn_search(r, query, best);
        }
    }
    /// Find all points within `radius` of `query`.
    pub fn range_query(&self, query: Vec3, radius: f64) -> Vec<usize> {
        if self.nodes.is_empty() {
            return Vec::new();
        }
        let mut result = Vec::new();
        let r2 = radius * radius;
        self.range_search(0, query, r2, &mut result);
        result
    }
    fn range_search(&self, node_idx: usize, query: Vec3, r2: f64, out: &mut Vec<usize>) {
        if node_idx >= self.nodes.len() {
            return;
        }
        let node = &self.nodes[node_idx];
        let ball_center = self.points[node.pivot];
        let dist_to_center = (ball_center - query).norm();
        if (dist_to_center - node.radius) * (dist_to_center - node.radius) > r2
            && dist_to_center > node.radius + r2.sqrt()
        {
            return;
        }
        if (self.points[node.pivot] - query).norm_squared() <= r2 {
            out.push(node.pivot);
        }
        for &i in &node.items {
            if (self.points[i] - query).norm_squared() <= r2 {
                out.push(i);
            }
        }
        if let Some(l) = node.left {
            self.range_search(l, query, r2, out);
        }
        if let Some(r) = node.right {
            self.range_search(r, query, r2, out);
        }
    }
    /// Number of points.
    pub fn len(&self) -> usize {
        self.points.len()
    }
    /// Return `true` if empty.
    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }
}
/// Axis-Aligned Bounding Box used by spatial data structures.
#[derive(Clone, Debug)]
pub struct SpatialAabb {
    /// Minimum corner of the box.
    pub min: Vec3,
    /// Maximum corner of the box.
    pub max: Vec3,
}
impl SpatialAabb {
    /// Create a new AABB from min/max corners.
    pub fn new(min: Vec3, max: Vec3) -> Self {
        Self { min, max }
    }
    /// Return the center of the box.
    pub fn center(&self) -> Vec3 {
        (self.min + self.max) * 0.5
    }
    /// Return half the extents (half-widths in each axis).
    pub fn half_extents(&self) -> Vec3 {
        (self.max - self.min) * 0.5
    }
    /// Return `true` if point `p` lies inside (or on the boundary of) this box.
    pub fn contains_point(&self, p: Vec3) -> bool {
        p.x >= self.min.x
            && p.x <= self.max.x
            && p.y >= self.min.y
            && p.y <= self.max.y
            && p.z >= self.min.z
            && p.z <= self.max.z
    }
    /// Return `true` if this box overlaps `other`.
    pub fn intersects(&self, other: &Self) -> bool {
        self.min.x <= other.max.x
            && self.max.x >= other.min.x
            && self.min.y <= other.max.y
            && self.max.y >= other.min.y
            && self.min.z <= other.max.z
            && self.max.z >= other.min.z
    }
    /// Return a copy of this box grown by `amount` in every direction.
    pub fn expand(&self, amount: f64) -> Self {
        let d = Vec3::new(amount, amount, amount);
        Self {
            min: self.min - d,
            max: self.max + d,
        }
    }
    /// Compute the squared distance from point `p` to the nearest point on
    /// the box surface (0 when the point is inside).
    fn sq_dist_to_point(&self, p: Vec3) -> f64 {
        let dx = (self.min.x - p.x).max(0.0).max(p.x - self.max.x);
        let dy = (self.min.y - p.y).max(0.0).max(p.y - self.max.y);
        let dz = (self.min.z - p.z).max(0.0).max(p.z - self.max.z);
        dx * dx + dy * dy + dz * dz
    }
    /// Return the eight child octants of this box split at `split`.
    fn octants(&self, split: Vec3) -> [SpatialAabb; 8] {
        [
            SpatialAabb::new(self.min, split),
            SpatialAabb::new(
                Vec3::new(split.x, self.min.y, self.min.z),
                Vec3::new(self.max.x, split.y, split.z),
            ),
            SpatialAabb::new(
                Vec3::new(self.min.x, split.y, self.min.z),
                Vec3::new(split.x, self.max.y, split.z),
            ),
            SpatialAabb::new(
                Vec3::new(split.x, split.y, self.min.z),
                Vec3::new(self.max.x, self.max.y, split.z),
            ),
            SpatialAabb::new(
                Vec3::new(self.min.x, self.min.y, split.z),
                Vec3::new(split.x, split.y, self.max.z),
            ),
            SpatialAabb::new(
                Vec3::new(split.x, self.min.y, split.z),
                Vec3::new(self.max.x, split.y, self.max.z),
            ),
            SpatialAabb::new(
                Vec3::new(self.min.x, split.y, split.z),
                Vec3::new(split.x, self.max.y, self.max.z),
            ),
            SpatialAabb::new(split, self.max),
        ]
    }
}
/// A KD-tree with lazy deletion support.
///
/// Deleted points are marked but not removed; the tree can be compacted.
pub struct KdTreeWithDeletion {
    pub(super) inner: KdTree,
    pub(super) deleted: Vec<bool>,
}
impl KdTreeWithDeletion {
    /// Build from a list of points.
    pub fn build(points: Vec<Vec3>) -> Self {
        let n = points.len();
        let inner = KdTree::build(points);
        Self {
            inner,
            deleted: vec![false; n],
        }
    }
    /// Mark point at index `i` as deleted.
    pub fn delete(&mut self, i: usize) {
        if i < self.deleted.len() {
            self.deleted[i] = true;
        }
    }
    /// Return the nearest non-deleted point.
    pub fn nearest_active(&self, query: Vec3) -> Option<(usize, f64)> {
        if self.inner.nodes.is_empty() {
            return None;
        }
        let mut best: Option<(usize, f64)> = None;
        let mut best_sq = f64::INFINITY;
        self.kd_nearest_active(0, query, &mut best, &mut best_sq);
        best
    }
    fn kd_nearest_active(
        &self,
        node_idx: usize,
        query: Vec3,
        best: &mut Option<(usize, f64)>,
        best_sq: &mut f64,
    ) {
        if node_idx >= self.inner.nodes.len() {
            return;
        }
        let node = &self.inner.nodes[node_idx];
        let p = self.inner.points[node.point_idx];
        if !self.deleted[node.point_idx] {
            let d2 = (p - query).norm_squared();
            if d2 < *best_sq {
                *best_sq = d2;
                *best = Some((node.point_idx, d2.sqrt()));
            }
        }
        let ax = node.axis as usize;
        let diff = query[ax] - p[ax];
        let (near, far) = if diff <= 0.0 {
            (node.left, node.right)
        } else {
            (node.right, node.left)
        };
        if let Some(n) = near {
            self.kd_nearest_active(n, query, best, best_sq);
        }
        if diff * diff < *best_sq
            && let Some(f) = far
        {
            self.kd_nearest_active(f, query, best, best_sq);
        }
    }
    /// Number of active (non-deleted) points.
    pub fn n_active(&self) -> usize {
        self.deleted.iter().filter(|&&d| !d).count()
    }
    /// Total number of points (including deleted).
    pub fn len(&self) -> usize {
        self.inner.len()
    }
    /// Return `true` if there are no active points.
    pub fn is_empty(&self) -> bool {
        self.n_active() == 0
    }
}
/// Octree for 3D spatial queries over arbitrary data.
pub struct Octree<T: Clone> {
    pub(super) root: OctreeNode<T>,
    pub(super) bounds: SpatialAabb,
    pub(super) max_depth: u32,
    pub(super) max_items_per_node: usize,
}
impl<T: Clone> Octree<T> {
    /// Create a new, empty octree covering `bounds`.
    pub fn new(bounds: SpatialAabb, max_depth: u32, max_items: usize) -> Self {
        Self {
            root: OctreeNode::Leaf { items: Vec::new() },
            bounds,
            max_depth,
            max_items_per_node: max_items,
        }
    }
    /// Insert `value` at position `point`.
    ///
    /// Points outside the octree bounds are silently discarded.
    pub fn insert(&mut self, point: Vec3, value: T) {
        if !self.bounds.contains_point(point) {
            return;
        }
        let max_depth = self.max_depth;
        let max_items = self.max_items_per_node;
        Self::node_insert(
            &mut self.root,
            point,
            value,
            &self.bounds,
            0,
            max_depth,
            max_items,
        );
    }
    fn node_insert(
        node: &mut OctreeNode<T>,
        point: Vec3,
        value: T,
        bounds: &SpatialAabb,
        depth: u32,
        max_depth: u32,
        max_items: usize,
    ) {
        match node {
            OctreeNode::Leaf { items } => {
                items.push((point, value));
                if items.len() > max_items && depth < max_depth {
                    let split = bounds.center();
                    let octants = bounds.octants(split);
                    let mut children: [OctreeNode<T>; 8] = [
                        OctreeNode::Leaf { items: Vec::new() },
                        OctreeNode::Leaf { items: Vec::new() },
                        OctreeNode::Leaf { items: Vec::new() },
                        OctreeNode::Leaf { items: Vec::new() },
                        OctreeNode::Leaf { items: Vec::new() },
                        OctreeNode::Leaf { items: Vec::new() },
                        OctreeNode::Leaf { items: Vec::new() },
                        OctreeNode::Leaf { items: Vec::new() },
                    ];
                    let drained: Vec<(Vec3, T)> = std::mem::take(items);
                    for (p, v) in drained {
                        let idx = Self::octant_index(p, split);
                        Self::node_insert(
                            &mut children[idx],
                            p,
                            v,
                            &octants[idx],
                            depth + 1,
                            max_depth,
                            max_items,
                        );
                    }
                    *node = OctreeNode::Internal {
                        children: Box::new(children),
                        split,
                    };
                }
            }
            OctreeNode::Internal { children, split } => {
                let split = *split;
                let idx = Self::octant_index(point, split);
                let octants = bounds.octants(split);
                Self::node_insert(
                    &mut children[idx],
                    point,
                    value,
                    &octants[idx],
                    depth + 1,
                    max_depth,
                    max_items,
                );
            }
        }
    }
    /// Return the child-octant index (0–7) for `point` relative to `split`.
    fn octant_index(point: Vec3, split: Vec3) -> usize {
        let x = if point.x >= split.x { 1 } else { 0 };
        let y = if point.y >= split.y { 2 } else { 0 };
        let z = if point.z >= split.z { 4 } else { 0 };
        x | y | z
    }
    /// Return all `(position, &value)` pairs whose position lies within the
    /// sphere defined by `center` and `radius`.
    pub fn query_sphere(&self, center: Vec3, radius: f64) -> Vec<(Vec3, &T)> {
        let mut results = Vec::new();
        let r2 = radius * radius;
        Self::node_query_sphere(&self.root, center, r2, &self.bounds, &mut results);
        results
    }
    fn node_query_sphere<'a>(
        node: &'a OctreeNode<T>,
        center: Vec3,
        r2: f64,
        bounds: &SpatialAabb,
        out: &mut Vec<(Vec3, &'a T)>,
    ) {
        if bounds.sq_dist_to_point(center) > r2 {
            return;
        }
        match node {
            OctreeNode::Leaf { items } => {
                for (p, v) in items {
                    if (p - center).norm_squared() <= r2 {
                        out.push((*p, v));
                    }
                }
            }
            OctreeNode::Internal { children, split } => {
                let octants = bounds.octants(*split);
                for (child, oct) in children.iter().zip(octants.iter()) {
                    Self::node_query_sphere(child, center, r2, oct, out);
                }
            }
        }
    }
    /// Return all `(position, &value)` pairs whose position overlaps `aabb`.
    pub fn query_aabb(&self, aabb: &SpatialAabb) -> Vec<(Vec3, &T)> {
        let mut results = Vec::new();
        Self::node_query_aabb(&self.root, aabb, &self.bounds, &mut results);
        results
    }
    fn node_query_aabb<'a>(
        node: &'a OctreeNode<T>,
        query: &SpatialAabb,
        bounds: &SpatialAabb,
        out: &mut Vec<(Vec3, &'a T)>,
    ) {
        if !bounds.intersects(query) {
            return;
        }
        match node {
            OctreeNode::Leaf { items } => {
                for (p, v) in items {
                    if query.contains_point(*p) {
                        out.push((*p, v));
                    }
                }
            }
            OctreeNode::Internal { children, split } => {
                let octants = bounds.octants(*split);
                for (child, oct) in children.iter().zip(octants.iter()) {
                    Self::node_query_aabb(child, query, oct, out);
                }
            }
        }
    }
    /// Return `(position, &value, distance)` of the item nearest to `query`.
    ///
    /// Returns `None` when the tree is empty.
    pub fn nearest_neighbor(&self, query: Vec3) -> Option<(Vec3, &T, f64)> {
        let mut best: Option<(Vec3, &T, f64)> = None;
        let mut best_sq = f64::INFINITY;
        Self::node_nearest(&self.root, query, &self.bounds, &mut best, &mut best_sq);
        best
    }
    fn node_nearest<'a>(
        node: &'a OctreeNode<T>,
        query: Vec3,
        bounds: &SpatialAabb,
        best: &mut Option<(Vec3, &'a T, f64)>,
        best_sq: &mut f64,
    ) {
        if bounds.sq_dist_to_point(query) >= *best_sq {
            return;
        }
        match node {
            OctreeNode::Leaf { items } => {
                for (p, v) in items {
                    let d2 = (p - query).norm_squared();
                    if d2 < *best_sq {
                        *best_sq = d2;
                        *best = Some((*p, v, d2.sqrt()));
                    }
                }
            }
            OctreeNode::Internal { children, split } => {
                let octants = bounds.octants(*split);
                let first = Self::octant_index(query, *split);
                let order: [usize; 8] = {
                    let mut o = [0usize; 8];
                    o[0] = first;
                    let mut k = 1;
                    for i in 0..8usize {
                        if i != first {
                            o[k] = i;
                            k += 1;
                        }
                    }
                    o
                };
                for idx in order {
                    Self::node_nearest(&children[idx], query, &octants[idx], best, best_sq);
                }
            }
        }
    }
    /// Total number of items stored in the tree.
    pub fn len(&self) -> usize {
        Self::node_len(&self.root)
    }
    fn node_len(node: &OctreeNode<T>) -> usize {
        match node {
            OctreeNode::Leaf { items } => items.len(),
            OctreeNode::Internal { children, .. } => children.iter().map(Self::node_len).sum(),
        }
    }
    /// Return `true` when the tree contains no items.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}
/// A simple bulk-loaded R-tree (sort-tile-recursive, STR) for 3D AABBs.
///
/// Only bulk loading is supported; the tree does not support incremental inserts.
/// Use [`RTree::query`] to find all items whose bounding boxes overlap a query AABB.
pub struct RTree<T: Clone> {
    pub(super) root: RTreeNode<T>,
}
impl<T: Clone> RTree<T> {
    /// Build a bulk-loaded R-tree from a list of `(bounding_box, value)` pairs.
    ///
    /// `fanout` is the maximum number of entries per node (typically 4–16).
    pub fn build(mut items: Vec<(SpatialAabb, T)>, fanout: usize) -> Self {
        let fanout = fanout.max(2);
        if items.is_empty() {
            let empty_aabb = SpatialAabb::new(Vec3::zeros(), Vec3::zeros());
            return RTree {
                root: RTreeNode::Leaf {
                    items: vec![],
                    bbox: empty_aabb,
                },
            };
        }
        let root = Self::build_node(&mut items, fanout, 0);
        RTree { root }
    }
    fn compute_bbox(items: &[(SpatialAabb, T)]) -> SpatialAabb {
        let mut mn = items[0].0.min;
        let mut mx = items[0].0.max;
        for (aabb, _) in items {
            mn = Vec3::new(
                mn.x.min(aabb.min.x),
                mn.y.min(aabb.min.y),
                mn.z.min(aabb.min.z),
            );
            mx = Vec3::new(
                mx.x.max(aabb.max.x),
                mx.y.max(aabb.max.y),
                mx.z.max(aabb.max.z),
            );
        }
        SpatialAabb::new(mn, mx)
    }
    fn build_node(items: &mut [(SpatialAabb, T)], fanout: usize, depth: usize) -> RTreeNode<T> {
        if items.len() <= fanout {
            let bbox = Self::compute_bbox(items);
            return RTreeNode::Leaf {
                items: items.to_vec(),
                bbox,
            };
        }
        let axis = depth % 3;
        items.sort_unstable_by(|a, b| {
            let ca = match axis {
                0 => (a.0.min.x + a.0.max.x) * 0.5,
                1 => (a.0.min.y + a.0.max.y) * 0.5,
                _ => (a.0.min.z + a.0.max.z) * 0.5,
            };
            let cb = match axis {
                0 => (b.0.min.x + b.0.max.x) * 0.5,
                1 => (b.0.min.y + b.0.max.y) * 0.5,
                _ => (b.0.min.z + b.0.max.z) * 0.5,
            };
            ca.partial_cmp(&cb).unwrap_or(std::cmp::Ordering::Equal)
        });
        let chunk = items.len().div_ceil(fanout);
        let mut children = Vec::new();
        let mut start = 0;
        while start < items.len() {
            let end = (start + chunk).min(items.len());
            children.push(Self::build_node(&mut items[start..end], fanout, depth + 1));
            start = end;
        }
        let bbox = {
            let mut mn = Vec3::new(f64::MAX, f64::MAX, f64::MAX);
            let mut mx = Vec3::new(f64::MIN, f64::MIN, f64::MIN);
            for child in &children {
                let cb = match child {
                    RTreeNode::Leaf { bbox, .. } => bbox,
                    RTreeNode::Internal { bbox, .. } => bbox,
                };
                mn = Vec3::new(mn.x.min(cb.min.x), mn.y.min(cb.min.y), mn.z.min(cb.min.z));
                mx = Vec3::new(mx.x.max(cb.max.x), mx.y.max(cb.max.y), mx.z.max(cb.max.z));
            }
            SpatialAabb::new(mn, mx)
        };
        RTreeNode::Internal { children, bbox }
    }
    /// Query all items whose bounding box overlaps `query`.
    pub fn query<'a>(&'a self, query: &SpatialAabb) -> Vec<&'a T> {
        let mut out = Vec::new();
        Self::node_query(&self.root, query, &mut out);
        out
    }
    fn node_query<'a>(node: &'a RTreeNode<T>, query: &SpatialAabb, out: &mut Vec<&'a T>) {
        match node {
            RTreeNode::Leaf { items, bbox } => {
                if !bbox.intersects(query) {
                    return;
                }
                for (aabb, val) in items {
                    if aabb.intersects(query) {
                        out.push(val);
                    }
                }
            }
            RTreeNode::Internal { children, bbox } => {
                if !bbox.intersects(query) {
                    return;
                }
                for child in children {
                    Self::node_query(child, query, out);
                }
            }
        }
    }
    /// Return the number of items stored in the tree.
    pub fn len(&self) -> usize {
        Self::node_len(&self.root)
    }
    fn node_len(node: &RTreeNode<T>) -> usize {
        match node {
            RTreeNode::Leaf { items, .. } => items.len(),
            RTreeNode::Internal { children, .. } => children.iter().map(Self::node_len).sum(),
        }
    }
    /// Return `true` if the tree contains no items.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}
pub(crate) struct KdNode {
    pub(crate) point_idx: usize,
    pub(crate) left: Option<usize>,
    pub(crate) right: Option<usize>,
    pub(crate) axis: u8,
}
/// LSH table for approximate nearest-neighbor queries in 3D.
///
/// Uses random hyperplane projections. Each hash function projects a point
/// onto a random unit vector and quantizes to a bucket.
pub struct LshIndex {
    /// Number of hash functions (bands).
    pub(super) n_hashes: usize,
    /// Bucket width w.
    pub(super) bucket_width: f64,
    /// Random projection vectors (one per hash function).
    pub(super) projections: Vec<Vec3>,
    /// Hash tables: maps (hash_fn_idx, bucket) → list of point indices.
    pub(super) tables: Vec<std::collections::HashMap<i64, Vec<usize>>>,
    /// Stored points.
    pub(super) points: Vec<Vec3>,
}
impl LshIndex {
    /// Create a new LSH index.
    ///
    /// `n_hashes` is the number of hash functions.
    /// `bucket_width` is the quantization width.
    /// `seed_projections` provides deterministic projection vectors (length = n_hashes * 3).
    pub fn new(n_hashes: usize, bucket_width: f64, seed_projections: &[[f64; 3]]) -> Self {
        let projections: Vec<Vec3> = seed_projections
            .iter()
            .take(n_hashes)
            .map(|v| {
                let p = Vec3::new(v[0], v[1], v[2]);
                let n = p.norm();
                if n > 1e-15 {
                    p / n
                } else {
                    Vec3::new(1.0, 0.0, 0.0)
                }
            })
            .collect();
        let n = projections.len();
        let tables = (0..n).map(|_| std::collections::HashMap::new()).collect();
        Self {
            n_hashes: n,
            bucket_width: bucket_width.max(1e-15),
            projections,
            tables,
            points: Vec::new(),
        }
    }
    fn hash_fn(&self, p: Vec3, k: usize) -> i64 {
        let proj = self.projections[k % self.n_hashes].dot(&p);
        (proj / self.bucket_width).floor() as i64
    }
    /// Insert a point and return its index.
    pub fn insert(&mut self, p: Vec3) -> usize {
        let idx = self.points.len();
        self.points.push(p);
        for k in 0..self.n_hashes {
            let bucket = self.hash_fn(p, k);
            self.tables[k].entry(bucket).or_default().push(idx);
        }
        idx
    }
    /// Approximate nearest-neighbor query.
    ///
    /// Returns the index of the approximate nearest point to `query`.
    pub fn query_approx_nn(&self, query: Vec3) -> Option<usize> {
        if self.points.is_empty() {
            return None;
        }
        let mut candidates = std::collections::HashSet::new();
        for k in 0..self.n_hashes {
            let bucket = self.hash_fn(query, k);
            if let Some(pts) = self.tables[k].get(&bucket) {
                for &i in pts {
                    candidates.insert(i);
                }
            }
            for delta in [-1_i64, 1] {
                if let Some(pts) = self.tables[k].get(&(bucket + delta)) {
                    for &i in pts {
                        candidates.insert(i);
                    }
                }
            }
        }
        if candidates.is_empty() {
            return (0..self.points.len()).min_by(|&a, &b| {
                let da = (self.points[a] - query).norm_squared();
                let db = (self.points[b] - query).norm_squared();
                da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
            });
        }
        candidates.into_iter().min_by(|&a, &b| {
            let da = (self.points[a] - query).norm_squared();
            let db = (self.points[b] - query).norm_squared();
            da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
        })
    }
    /// Number of stored points.
    pub fn len(&self) -> usize {
        self.points.len()
    }
    /// Return `true` if no points are stored.
    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }
}
pub(crate) enum OctreeNode<T: Clone> {
    /// Leaf node holding a flat list of (position, value) pairs.
    Leaf { items: Vec<(Vec3, T)> },
    /// Internal node split at `split`; holds exactly 8 children.
    Internal {
        children: Box<[OctreeNode<T>; 8]>,
        split: Vec3,
    },
}
pub(super) struct BallTreeNode {
    /// Index of the pivot point.
    pub(super) pivot: usize,
    /// Radius of the ball bounding all points in this subtree.
    pub(super) radius: f64,
    /// Left child (None = leaf).
    pub(super) left: Option<usize>,
    /// Right child (None = leaf).
    pub(super) right: Option<usize>,
    /// Points stored at this leaf node (empty for internal nodes).
    pub(super) items: Vec<usize>,
}
/// A voxel grid spatial index.
///
/// Divides space into axis-aligned voxels of size `voxel_size` and stores
/// point indices per voxel. Supports fast insertion and range queries.
pub struct VoxelGrid {
    pub(super) voxel_size: f64,
    pub(super) cells: std::collections::HashMap<(i64, i64, i64), Vec<usize>>,
    pub(super) points: Vec<Vec3>,
}
impl VoxelGrid {
    /// Create a new voxel grid.
    pub fn new(voxel_size: f64) -> Self {
        Self {
            voxel_size: voxel_size.max(1e-15),
            cells: std::collections::HashMap::new(),
            points: Vec::new(),
        }
    }
    fn voxel_of(&self, p: Vec3) -> (i64, i64, i64) {
        (
            (p.x / self.voxel_size).floor() as i64,
            (p.y / self.voxel_size).floor() as i64,
            (p.z / self.voxel_size).floor() as i64,
        )
    }
    /// Insert a point, returning its index.
    pub fn insert(&mut self, p: Vec3) -> usize {
        let idx = self.points.len();
        let voxel = self.voxel_of(p);
        self.cells.entry(voxel).or_default().push(idx);
        self.points.push(p);
        idx
    }
    /// Query points within distance `r` of `center`.
    pub fn range_query(&self, center: Vec3, r: f64) -> Vec<usize> {
        let r2 = r * r;
        let span = (r / self.voxel_size).ceil() as i64;
        let (cx, cy, cz) = self.voxel_of(center);
        let mut result = Vec::new();
        for dx in -span..=span {
            for dy in -span..=span {
                for dz in -span..=span {
                    let key = (cx + dx, cy + dy, cz + dz);
                    if let Some(pts) = self.cells.get(&key) {
                        for &i in pts {
                            if (self.points[i] - center).norm_squared() <= r2 {
                                result.push(i);
                            }
                        }
                    }
                }
            }
        }
        result
    }
    /// Nearest neighbor linear scan (fallback for small sets).
    pub fn nearest_linear(&self, query: Vec3) -> Option<(usize, f64)> {
        self.points
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| {
                ((*a) - query)
                    .norm_squared()
                    .partial_cmp(&((*b) - query).norm_squared())
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|(i, p)| (i, ((*p) - query).norm()))
    }
    /// Number of points.
    pub fn len(&self) -> usize {
        self.points.len()
    }
    /// Return `true` if empty.
    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }
    /// Downsampling: return one representative point per voxel (voxel centroid).
    pub fn voxel_centroids(&self) -> Vec<Vec3> {
        self.cells
            .values()
            .map(|indices| {
                let sum: Vec3 = indices
                    .iter()
                    .map(|&i| self.points[i])
                    .fold(Vec3::zeros(), |a, b| a + b);
                sum / indices.len() as f64
            })
            .collect()
    }
}
/// Uniform-grid spatial index for fast range queries.
///
/// Divides space into a regular grid and stores point indices per cell.
/// Supports range queries in O(cells) time.
pub struct GridSpatialIndex {
    /// Cell size.
    pub(super) cell_size: f64,
    /// Grid cells: maps (ix, iy, iz) to list of point indices.
    pub(super) cells: std::collections::HashMap<(i64, i64, i64), Vec<usize>>,
    /// Stored points.
    pub(super) points: Vec<Vec3>,
}
impl GridSpatialIndex {
    /// Create a new grid spatial index with the given cell size.
    pub fn new(cell_size: f64) -> Self {
        Self {
            cell_size: cell_size.max(1e-15),
            cells: std::collections::HashMap::new(),
            points: Vec::new(),
        }
    }
    fn cell_of(&self, p: Vec3) -> (i64, i64, i64) {
        (
            (p.x / self.cell_size).floor() as i64,
            (p.y / self.cell_size).floor() as i64,
            (p.z / self.cell_size).floor() as i64,
        )
    }
    /// Insert a point.
    pub fn insert(&mut self, p: Vec3) -> usize {
        let idx = self.points.len();
        let cell = self.cell_of(p);
        self.cells.entry(cell).or_default().push(idx);
        self.points.push(p);
        idx
    }
    /// Range query: returns indices of all points within distance `r` of `center`.
    pub fn range_query(&self, center: Vec3, r: f64) -> Vec<usize> {
        let r_sq = r * r;
        let half_cells = (r / self.cell_size).ceil() as i64;
        let (cx, cy, cz) = self.cell_of(center);
        let mut result = Vec::new();
        for dx in -half_cells..=half_cells {
            for dy in -half_cells..=half_cells {
                for dz in -half_cells..=half_cells {
                    let key = (cx + dx, cy + dy, cz + dz);
                    if let Some(indices) = self.cells.get(&key) {
                        for &i in indices {
                            let p = self.points[i];
                            let diff = p - center;
                            if diff.norm_squared() <= r_sq {
                                result.push(i);
                            }
                        }
                    }
                }
            }
        }
        result
    }
    /// Number of points.
    pub fn len(&self) -> usize {
        self.points.len()
    }
    /// Return `true` if empty.
    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }
    /// Get all points within a given cell (for debugging/testing).
    pub fn cell_contents(&self, cell: (i64, i64, i64)) -> Vec<usize> {
        self.cells.get(&cell).cloned().unwrap_or_default()
    }
}
/// Morton-order sorted point index structure.
///
/// Organizes 3D points by Morton code order for cache-efficient spatial queries
/// and octree construction.
pub struct MortonSortedIndex {
    /// Points with their Morton codes, sorted ascending.
    pub(super) sorted: Vec<(u64, usize)>,
}
impl MortonSortedIndex {
    /// Build a Morton-sorted index from a list of points.
    ///
    /// `n_cells` sets the grid resolution (e.g., 1024 cells per axis).
    pub fn build(points: &[Vec3], n_cells: u32) -> Self {
        let n_cells = n_cells.max(1);
        let (min, max) = if points.is_empty() {
            (Vec3::zeros(), Vec3::zeros())
        } else {
            let mut mn = points[0];
            let mut mx = points[0];
            for &p in points {
                for i in 0..3 {
                    if p[i] < mn[i] {
                        mn[i] = p[i];
                    }
                    if p[i] > mx[i] {
                        mx[i] = p[i];
                    }
                }
            }
            (mn, mx)
        };
        let extent = max - min;
        let scale = if extent.norm() < 1e-12 {
            Vec3::new(1.0, 1.0, 1.0)
        } else {
            Vec3::new(
                if extent.x > 1e-12 {
                    n_cells as f64 / extent.x
                } else {
                    1.0
                },
                if extent.y > 1e-12 {
                    n_cells as f64 / extent.y
                } else {
                    1.0
                },
                if extent.z > 1e-12 {
                    n_cells as f64 / extent.z
                } else {
                    1.0
                },
            )
        };
        let mut sorted: Vec<(u64, usize)> = points
            .iter()
            .enumerate()
            .map(|(i, &p)| {
                let rel = p - min;
                let ix = ((rel.x * scale.x) as u32).min(n_cells - 1);
                let iy = ((rel.y * scale.y) as u32).min(n_cells - 1);
                let iz = ((rel.z * scale.z) as u32).min(n_cells - 1);
                (morton_encode(ix, iy, iz), i)
            })
            .collect();
        sorted.sort_unstable_by_key(|&(code, _)| code);
        MortonSortedIndex { sorted }
    }
    /// Number of points.
    pub fn len(&self) -> usize {
        self.sorted.len()
    }
    /// Return `true` if empty.
    pub fn is_empty(&self) -> bool {
        self.sorted.is_empty()
    }
    /// Return the Morton code for point at original index `i`.
    pub fn morton_code(&self, original_idx: usize) -> Option<u64> {
        self.sorted
            .iter()
            .find(|&&(_, i)| i == original_idx)
            .map(|&(code, _)| code)
    }
    /// Return indices in Morton (Z-curve) order.
    pub fn sorted_indices(&self) -> Vec<usize> {
        self.sorted.iter().map(|&(_, i)| i).collect()
    }
}
