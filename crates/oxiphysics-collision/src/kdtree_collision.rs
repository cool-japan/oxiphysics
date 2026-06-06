// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! KD-tree based spatial queries and collision detection.
//!
//! Provides a median-split KD-tree ([`KdTree`]) that supports nearest-neighbour
//! search, k-nearest neighbours, radius queries, and self-collision pair
//! detection.  A dynamic wrapper ([`KdTreeCollisionDetector`]) allows
//! incremental insertions followed by bulk rebuilds.

use std::collections::BinaryHeap;

// ─────────────────────────────────────────────────────────────────────────────
// Geometry helpers
// ─────────────────────────────────────────────────────────────────────────────

#[inline]
fn dist_sq(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    (a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)
}

// ─────────────────────────────────────────────────────────────────────────────
// Aabb3
// ─────────────────────────────────────────────────────────────────────────────

/// Axis-aligned bounding box in 3-D.
#[derive(Debug, Clone, PartialEq)]
pub struct Aabb3 {
    /// Minimum corner `[x, y, z]`.
    pub min: [f64; 3],
    /// Maximum corner `[x, y, z]`.
    pub max: [f64; 3],
}

impl Aabb3 {
    /// Create a degenerate (point) AABB at the origin.
    pub fn empty() -> Self {
        Self {
            min: [f64::INFINITY; 3],
            max: [f64::NEG_INFINITY; 3],
        }
    }

    /// Create an AABB from explicit min/max corners.
    pub fn new(min: [f64; 3], max: [f64; 3]) -> Self {
        Self { min, max }
    }

    /// Expand this AABB to contain `point`.
    pub fn expand(&mut self, point: &[f64; 3]) {
        self.min.iter_mut().zip(point.iter()).for_each(|(m, &p)| {
            if p < *m {
                *m = p;
            }
        });
        self.max.iter_mut().zip(point.iter()).for_each(|(m, &p)| {
            if p > *m {
                *m = p;
            }
        });
    }

    /// Return `true` if this AABB overlaps `other`.
    pub fn overlaps(&self, other: &Self) -> bool {
        self.min
            .iter()
            .zip(self.max.iter())
            .zip(other.min.iter().zip(other.max.iter()))
            .all(|((&mn, &mx), (&omn, &omx))| mn <= omx && mx >= omn)
    }

    /// Return `true` if `point` is inside or on the boundary of this AABB.
    pub fn contains_point(&self, point: &[f64; 3]) -> bool {
        self.min
            .iter()
            .zip(self.max.iter())
            .zip(point.iter())
            .all(|((&mn, &mx), &p)| p >= mn && p <= mx)
    }

    /// Return the squared minimum distance from `point` to this AABB.
    pub fn min_dist_sq(&self, point: &[f64; 3]) -> f64 {
        self.min
            .iter()
            .zip(self.max.iter())
            .zip(point.iter())
            .map(|((&mn, &mx), &p)| {
                if p < mn {
                    (mn - p).powi(2)
                } else if p > mx {
                    (p - mx).powi(2)
                } else {
                    0.0
                }
            })
            .sum()
    }

    /// Compute the AABB that encloses a slice of points.
    pub fn from_points(pts: &[[f64; 3]]) -> Self {
        let mut aabb = Self::empty();
        for p in pts {
            aabb.expand(p);
        }
        aabb
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// KdPoint trait
// ─────────────────────────────────────────────────────────────────────────────

/// Any type that can act as a 3-D point in a KD-tree.
pub trait KdPoint {
    /// Return the 3-D position of this point.
    fn position(&self) -> [f64; 3];
}

impl KdPoint for [f64; 3] {
    fn position(&self) -> [f64; 3] {
        *self
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// KdNode
// ─────────────────────────────────────────────────────────────────────────────

/// A node in a KD-tree.
#[derive(Debug)]
pub enum KdNode {
    /// A leaf holding a small bucket of point indices.
    Leaf {
        /// Indices into the owning [`KdTree`]'s points array.
        indices: Vec<usize>,
        /// Tight AABB around the leaf points.
        aabb: Aabb3,
    },
    /// An internal split node.
    Internal {
        /// Axis along which the split was made (0 = X, 1 = Y, 2 = Z).
        split_dim: usize,
        /// The split value (median coordinate along `split_dim`).
        split_val: f64,
        /// Left subtree (values ≤ `split_val`).
        left: Box<KdNode>,
        /// Right subtree (values > `split_val`).
        right: Box<KdNode>,
        /// Tight AABB for all points in this subtree.
        aabb: Aabb3,
    },
}

impl KdNode {
    /// Return the AABB of this node.
    pub fn aabb(&self) -> &Aabb3 {
        match self {
            KdNode::Leaf { aabb, .. } => aabb,
            KdNode::Internal { aabb, .. } => aabb,
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// KdTree
// ─────────────────────────────────────────────────────────────────────────────

/// Maximum number of points per leaf bucket.
const LEAF_SIZE: usize = 8;

/// KD-tree for 3-D point sets.
#[derive(Debug)]
pub struct KdTree {
    /// The root node of the tree (or `None` when the tree is empty).
    pub root: Option<KdNode>,
    /// The point coordinates stored in the tree.
    pub points: Vec<[f64; 3]>,
}

impl KdTree {
    /// Build a KD-tree from a set of 3-D points using median splits.
    pub fn build(points: Vec<[f64; 3]>) -> Self {
        if points.is_empty() {
            return Self { root: None, points };
        }
        let n = points.len();
        let mut indices: Vec<usize> = (0..n).collect();
        let root = Some(build_node(&points, &mut indices));
        Self { root, points }
    }

    /// Find the single nearest neighbour to `query`.
    ///
    /// Returns `Some((index, dist_sq))` or `None` when the tree is empty.
    pub fn nearest_neighbor(&self, query: &[f64; 3]) -> Option<(usize, f64)> {
        let root = self.root.as_ref()?;
        let mut best = (usize::MAX, f64::INFINITY);
        nn_search(root, query, &self.points, &mut best);
        if best.0 == usize::MAX {
            None
        } else {
            Some(best)
        }
    }

    /// Return the `k` nearest neighbours sorted by ascending squared distance.
    pub fn k_nearest(&self, query: &[f64; 3], k: usize) -> Vec<(usize, f64)> {
        if k == 0 {
            return vec![];
        }
        let root = match &self.root {
            Some(r) => r,
            None => return vec![],
        };
        // Max-heap ordered by dist_sq (negated for a min-heap behaviour).
        // We use OrderedFloat-like wrapper via std BinaryHeap<(OrdF64, usize)>.
        let mut heap: BinaryHeap<OrdF64Pair> = BinaryHeap::new();
        knn_search(root, query, &self.points, k, &mut heap);
        let mut result: Vec<(usize, f64)> = heap.into_iter().map(|p| (p.idx, p.dist_sq)).collect();
        result.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        result
    }

    /// Return all point indices whose distance to `center` is ≤ `radius`.
    pub fn range_query(&self, center: &[f64; 3], radius: f64) -> Vec<usize> {
        let root = match &self.root {
            Some(r) => r,
            None => return vec![],
        };
        let r2 = radius * radius;
        let mut result = Vec::new();
        range_search(root, center, r2, &self.points, &mut result);
        result
    }

    /// Find all pairs `(i, j)` with `i < j` where the distance between
    /// points `i` and `j` is ≤ `radius`.
    pub fn self_collision_pairs(&self, radius: f64) -> Vec<(usize, usize)> {
        let _n = self.points.len();
        let r2 = radius * radius;
        let mut pairs = Vec::new();
        for (i, pt) in self.points.iter().enumerate() {
            let candidates = self.range_query(pt, radius);
            for j in candidates {
                if j > i && dist_sq(pt, &self.points[j]) <= r2 {
                    pairs.push((i, j));
                }
            }
        }
        pairs.sort_unstable();
        pairs
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tree-building internals
// ─────────────────────────────────────────────────────────────────────────────

fn build_node(points: &[[f64; 3]], indices: &mut [usize]) -> KdNode {
    let aabb = Aabb3::from_points(&indices.iter().map(|&i| points[i]).collect::<Vec<_>>());

    if indices.len() <= LEAF_SIZE {
        return KdNode::Leaf {
            indices: indices.to_vec(),
            aabb,
        };
    }

    // Choose split axis: the axis with the widest span (sliding midpoint
    // strategy).
    let mut split_dim = 0;
    let mut max_span = aabb.max[0] - aabb.min[0];
    for d in 1..3 {
        let span = aabb.max[d] - aabb.min[d];
        if span > max_span {
            max_span = span;
            split_dim = d;
        }
    }

    // Median split.
    let mid = indices.len() / 2;
    indices.select_nth_unstable_by(mid, |&a, &b| {
        points[a][split_dim]
            .partial_cmp(&points[b][split_dim])
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let split_val = points[indices[mid]][split_dim];

    let (left_idx, right_idx) = indices.split_at_mut(mid);
    let left = Box::new(build_node(points, left_idx));
    let right = Box::new(build_node(points, right_idx));

    KdNode::Internal {
        split_dim,
        split_val,
        left,
        right,
        aabb,
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Search routines
// ─────────────────────────────────────────────────────────────────────────────

fn nn_search(node: &KdNode, query: &[f64; 3], points: &[[f64; 3]], best: &mut (usize, f64)) {
    match node {
        KdNode::Leaf { indices, .. } => {
            for &i in indices {
                let d = dist_sq(query, &points[i]);
                if d < best.1 {
                    *best = (i, d);
                }
            }
        }
        KdNode::Internal {
            split_dim,
            split_val,
            left,
            right,
            ..
        } => {
            let go_left = query[*split_dim] <= *split_val;
            let (near, far) = if go_left {
                (left.as_ref(), right.as_ref())
            } else {
                (right.as_ref(), left.as_ref())
            };
            nn_search(near, query, points, best);
            // Prune: only visit the far side if its AABB might beat the current best.
            if far.aabb().min_dist_sq(query) < best.1 {
                nn_search(far, query, points, best);
            }
        }
    }
}

/// Pair stored in max-heap for k-NN (largest dist_sq at top).
struct OrdF64Pair {
    dist_sq: f64,
    idx: usize,
}

impl PartialEq for OrdF64Pair {
    fn eq(&self, other: &Self) -> bool {
        self.dist_sq == other.dist_sq
    }
}
impl Eq for OrdF64Pair {}
impl PartialOrd for OrdF64Pair {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for OrdF64Pair {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.dist_sq
            .partial_cmp(&other.dist_sq)
            .unwrap_or(std::cmp::Ordering::Equal)
    }
}

fn knn_search(
    node: &KdNode,
    query: &[f64; 3],
    points: &[[f64; 3]],
    k: usize,
    heap: &mut BinaryHeap<OrdF64Pair>,
) {
    // Prune: if the heap is full and the AABB is farther than the worst
    // current best, skip.
    let worst = heap.peek().map(|p| p.dist_sq).unwrap_or(f64::INFINITY);
    if node.aabb().min_dist_sq(query) >= worst && heap.len() >= k {
        return;
    }
    match node {
        KdNode::Leaf { indices, .. } => {
            for &i in indices {
                let d = dist_sq(query, &points[i]);
                if heap.len() < k {
                    heap.push(OrdF64Pair { dist_sq: d, idx: i });
                } else if heap.peek().is_none_or(|top| d < top.dist_sq) {
                    heap.pop();
                    heap.push(OrdF64Pair { dist_sq: d, idx: i });
                }
            }
        }
        KdNode::Internal {
            split_dim,
            split_val,
            left,
            right,
            ..
        } => {
            let go_left = query[*split_dim] <= *split_val;
            let (near, far) = if go_left {
                (left.as_ref(), right.as_ref())
            } else {
                (right.as_ref(), left.as_ref())
            };
            knn_search(near, query, points, k, heap);
            knn_search(far, query, points, k, heap);
        }
    }
}

fn range_search(
    node: &KdNode,
    center: &[f64; 3],
    r2: f64,
    points: &[[f64; 3]],
    result: &mut Vec<usize>,
) {
    if node.aabb().min_dist_sq(center) > r2 {
        return;
    }
    match node {
        KdNode::Leaf { indices, .. } => {
            for &i in indices {
                if dist_sq(center, &points[i]) <= r2 {
                    result.push(i);
                }
            }
        }
        KdNode::Internal { left, right, .. } => {
            range_search(left, center, r2, points, result);
            range_search(right, center, r2, points, result);
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// KdTreeCollisionDetector
// ─────────────────────────────────────────────────────────────────────────────

/// Dynamic KD-tree collision detector that supports incremental insertion and
/// bulk rebuild.
#[derive(Debug)]
pub struct KdTreeCollisionDetector {
    /// Points collected since the last rebuild.
    pending: Vec<([f64; 3], usize)>,
    /// The most recently built KD-tree (may be stale if points were added
    /// after the last [`Self::rebuild`] call).
    tree: KdTree,
    /// User-supplied IDs for each point in the tree.
    ids: Vec<usize>,
}

impl KdTreeCollisionDetector {
    /// Create a new detector with a capacity hint.
    pub fn new(capacity: usize) -> Self {
        Self {
            pending: Vec::with_capacity(capacity),
            tree: KdTree::build(vec![]),
            ids: Vec::with_capacity(capacity),
        }
    }

    /// Insert a point at `pos` with user-supplied `id`.
    pub fn insert(&mut self, pos: [f64; 3], id: usize) {
        self.pending.push((pos, id));
    }

    /// Rebuild the internal KD-tree from all inserted points.
    ///
    /// Returns `&mut Self` for chaining.
    pub fn rebuild(&mut self) -> &mut Self {
        let (positions, ids): (Vec<[f64; 3]>, Vec<usize>) = self.pending.iter().cloned().unzip();
        self.ids = ids;
        self.tree = KdTree::build(positions);
        self
    }

    /// Query all points within `r` of `pos`.
    ///
    /// Returns user-supplied IDs.  The tree must have been rebuilt after the
    /// last insertion for results to be accurate.
    pub fn query_radius(&self, pos: &[f64; 3], r: f64) -> Vec<usize> {
        self.tree
            .range_query(pos, r)
            .into_iter()
            .map(|tree_idx| self.ids[tree_idx])
            .collect()
    }

    /// Return the number of inserted points.
    pub fn len(&self) -> usize {
        self.pending.len()
    }

    /// Return `true` if no points have been inserted.
    pub fn is_empty(&self) -> bool {
        self.pending.is_empty()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// BVH-style batch leaf grouping
// ─────────────────────────────────────────────────────────────────────────────

/// A flat list of AABB leaf groups for batch broadphase queries.
///
/// Constructed by collecting the leaf AABBs from a [`KdTree`].
#[derive(Debug, Clone)]
pub struct BvhLeafGroups {
    /// Each entry is the AABB of one leaf in the KD-tree.
    pub groups: Vec<Aabb3>,
}

impl BvhLeafGroups {
    /// Extract leaf AABBs from a [`KdTree`].
    pub fn from_tree(tree: &KdTree) -> Self {
        let mut groups = Vec::new();
        if let Some(root) = &tree.root {
            collect_leaf_aabbs(root, &mut groups);
        }
        Self { groups }
    }

    /// Return the indices of leaf groups whose AABB overlaps the query sphere.
    pub fn query_sphere(&self, center: &[f64; 3], radius: f64) -> Vec<usize> {
        let r2 = radius * radius;
        self.groups
            .iter()
            .enumerate()
            .filter(|(_, aabb)| aabb.min_dist_sq(center) <= r2)
            .map(|(i, _)| i)
            .collect()
    }
}

fn collect_leaf_aabbs(node: &KdNode, out: &mut Vec<Aabb3>) {
    match node {
        KdNode::Leaf { aabb, .. } => out.push(aabb.clone()),
        KdNode::Internal { left, right, .. } => {
            collect_leaf_aabbs(left, out);
            collect_leaf_aabbs(right, out);
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── helpers ────────────────────────────────────────────────────────────

    fn grid_points(n: usize) -> Vec<[f64; 3]> {
        let side = (n as f64).cbrt().ceil() as usize;
        let mut pts = Vec::new();
        'outer: for x in 0..side {
            for y in 0..side {
                for z in 0..side {
                    pts.push([x as f64, y as f64, z as f64]);
                    if pts.len() == n {
                        break 'outer;
                    }
                }
            }
        }
        pts
    }

    // ── Aabb3 ──────────────────────────────────────────────────────────────

    #[test]
    fn test_aabb_empty_is_inverted() {
        let aabb = Aabb3::empty();
        assert!(aabb.min[0] > aabb.max[0]);
    }

    #[test]
    fn test_aabb_expand() {
        let mut aabb = Aabb3::empty();
        aabb.expand(&[1.0, 2.0, 3.0]);
        aabb.expand(&[-1.0, 0.0, 5.0]);
        assert_eq!(aabb.min, [-1.0, 0.0, 3.0]);
        assert_eq!(aabb.max, [1.0, 2.0, 5.0]);
    }

    #[test]
    fn test_aabb_overlaps_true() {
        let a = Aabb3::new([0.0; 3], [2.0; 3]);
        let b = Aabb3::new([1.0; 3], [3.0; 3]);
        assert!(a.overlaps(&b));
    }

    #[test]
    fn test_aabb_overlaps_false() {
        let a = Aabb3::new([0.0; 3], [1.0; 3]);
        let b = Aabb3::new([2.0; 3], [3.0; 3]);
        assert!(!a.overlaps(&b));
    }

    #[test]
    fn test_aabb_contains_point() {
        let aabb = Aabb3::new([0.0; 3], [1.0; 3]);
        assert!(aabb.contains_point(&[0.5, 0.5, 0.5]));
        assert!(!aabb.contains_point(&[1.5, 0.5, 0.5]));
    }

    #[test]
    fn test_aabb_min_dist_sq_inside() {
        let aabb = Aabb3::new([0.0; 3], [1.0; 3]);
        assert_eq!(aabb.min_dist_sq(&[0.5, 0.5, 0.5]), 0.0);
    }

    #[test]
    fn test_aabb_min_dist_sq_outside() {
        let aabb = Aabb3::new([0.0; 3], [1.0; 3]);
        let d = aabb.min_dist_sq(&[2.0, 0.5, 0.5]);
        assert!((d - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_aabb_from_points() {
        let pts = vec![[1.0, 2.0, 3.0], [-1.0, 0.0, 5.0], [0.0, 4.0, -1.0]];
        let aabb = Aabb3::from_points(&pts);
        assert_eq!(aabb.min, [-1.0, 0.0, -1.0]);
        assert_eq!(aabb.max, [1.0, 4.0, 5.0]);
    }

    // ── KdTree build & basics ───────────────────────────────────────────────

    #[test]
    fn test_kdtree_empty() {
        let tree = KdTree::build(vec![]);
        assert!(tree.root.is_none());
        assert!(tree.nearest_neighbor(&[0.0; 3]).is_none());
    }

    #[test]
    fn test_kdtree_single_point() {
        let tree = KdTree::build(vec![[1.0, 2.0, 3.0]]);
        let nn = tree.nearest_neighbor(&[0.0; 3]).unwrap();
        assert_eq!(nn.0, 0);
    }

    #[test]
    fn test_kdtree_build_grid() {
        let pts = grid_points(64);
        let tree = KdTree::build(pts);
        assert!(tree.root.is_some());
    }

    // ── nearest_neighbor ───────────────────────────────────────────────────

    #[test]
    fn test_nn_exact_match() {
        let pts = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let tree = KdTree::build(pts);
        let (idx, d) = tree.nearest_neighbor(&[1.0, 0.0, 0.0]).unwrap();
        assert_eq!(idx, 1);
        assert!(d < 1e-12);
    }

    #[test]
    fn test_nn_closest_of_many() {
        let pts = grid_points(27);
        let tree = KdTree::build(pts.clone());
        let query = [1.1, 1.1, 1.1];
        let (idx, _) = tree.nearest_neighbor(&query).unwrap();
        // Brute-force check
        let bf_idx = pts
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| dist_sq(a, &query).partial_cmp(&dist_sq(b, &query)).unwrap())
            .unwrap()
            .0;
        assert_eq!(idx, bf_idx);
    }

    #[test]
    fn test_nn_large_set() {
        let pts: Vec<[f64; 3]> = (0..500).map(|i| [i as f64 * 0.1, 0.0, 0.0]).collect();
        let tree = KdTree::build(pts);
        let (idx, _) = tree.nearest_neighbor(&[25.05, 0.0, 0.0]).unwrap();
        // Nearest should be around index 250.
        assert!((248..=252).contains(&idx));
    }

    // ── k_nearest ──────────────────────────────────────────────────────────

    #[test]
    fn test_knn_k_zero() {
        let pts = grid_points(8);
        let tree = KdTree::build(pts);
        assert!(tree.k_nearest(&[0.0; 3], 0).is_empty());
    }

    #[test]
    fn test_knn_k_equals_n() {
        let pts = grid_points(8);
        let n = pts.len();
        let tree = KdTree::build(pts);
        let result = tree.k_nearest(&[0.0; 3], n);
        assert_eq!(result.len(), n);
    }

    #[test]
    fn test_knn_sorted_ascending() {
        let pts = grid_points(27);
        let tree = KdTree::build(pts);
        let result = tree.k_nearest(&[1.5, 1.5, 1.5], 5);
        for w in result.windows(2) {
            assert!(w[0].1 <= w[1].1);
        }
    }

    #[test]
    fn test_knn_matches_brute_force() {
        let pts: Vec<[f64; 3]> = (0..30).map(|i| [i as f64, 0.0, 0.0]).collect();
        let query = [14.3, 0.0, 0.0];
        let tree = KdTree::build(pts.clone());
        let knn = tree.k_nearest(&query, 3);
        // Brute force top-3
        let mut bf: Vec<(usize, f64)> = pts
            .iter()
            .enumerate()
            .map(|(i, p)| (i, dist_sq(p, &query)))
            .collect();
        bf.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
        bf.truncate(3);
        let knn_idxs: Vec<usize> = knn.iter().map(|&(i, _)| i).collect();
        let bf_idxs: Vec<usize> = bf.iter().map(|&(i, _)| i).collect();
        assert_eq!(knn_idxs, bf_idxs);
    }

    // ── range_query ─────────────────────────────────────────────────────────

    #[test]
    fn test_range_query_empty_tree() {
        let tree = KdTree::build(vec![]);
        assert!(tree.range_query(&[0.0; 3], 1.0).is_empty());
    }

    #[test]
    fn test_range_query_all_in_radius() {
        let pts = vec![[0.0, 0.0, 0.0], [0.1, 0.0, 0.0], [0.0, 0.1, 0.0]];
        let tree = KdTree::build(pts.clone());
        let mut result = tree.range_query(&[0.05, 0.05, 0.0], 1.0);
        result.sort_unstable();
        assert_eq!(result, vec![0, 1, 2]);
    }

    #[test]
    fn test_range_query_none_in_radius() {
        let pts = vec![[10.0, 0.0, 0.0], [20.0, 0.0, 0.0]];
        let tree = KdTree::build(pts);
        assert!(tree.range_query(&[0.0; 3], 5.0).is_empty());
    }

    #[test]
    fn test_range_query_matches_brute_force() {
        let pts = grid_points(64);
        let tree = KdTree::build(pts.clone());
        let center = [2.5, 2.5, 2.5];
        let r = 1.8;
        let r2 = r * r;
        let mut kd_result = tree.range_query(&center, r);
        kd_result.sort_unstable();
        let mut bf: Vec<usize> = pts
            .iter()
            .enumerate()
            .filter(|(_, p)| dist_sq(p, &center) <= r2)
            .map(|(i, _)| i)
            .collect();
        bf.sort_unstable();
        assert_eq!(kd_result, bf);
    }

    // ── self_collision_pairs ─────────────────────────────────────────────────

    #[test]
    fn test_self_collision_no_pairs_far_apart() {
        let pts = vec![[0.0, 0.0, 0.0], [100.0, 0.0, 0.0], [200.0, 0.0, 0.0]];
        let tree = KdTree::build(pts);
        assert!(tree.self_collision_pairs(1.0).is_empty());
    }

    #[test]
    fn test_self_collision_all_close() {
        let pts = vec![[0.0; 3], [0.1, 0.0, 0.0], [0.0, 0.1, 0.0]];
        let tree = KdTree::build(pts);
        let pairs = tree.self_collision_pairs(0.2);
        // Expect 3 pairs: (0,1), (0,2), (1,2)
        assert_eq!(pairs.len(), 3);
    }

    #[test]
    fn test_self_collision_pairs_ordered() {
        let pts = grid_points(16);
        let tree = KdTree::build(pts);
        let pairs = tree.self_collision_pairs(1.5);
        for &(a, b) in &pairs {
            assert!(a < b);
        }
    }

    // ── KdTreeCollisionDetector ─────────────────────────────────────────────

    #[test]
    fn test_detector_empty() {
        let det = KdTreeCollisionDetector::new(10);
        assert!(det.is_empty());
        assert_eq!(det.len(), 0);
    }

    #[test]
    fn test_detector_insert_and_rebuild() {
        let mut det = KdTreeCollisionDetector::new(4);
        det.insert([0.0, 0.0, 0.0], 10);
        det.insert([1.0, 0.0, 0.0], 20);
        det.rebuild();
        assert_eq!(det.len(), 2);
    }

    #[test]
    fn test_detector_query_radius() {
        let mut det = KdTreeCollisionDetector::new(5);
        det.insert([0.0, 0.0, 0.0], 1);
        det.insert([0.5, 0.0, 0.0], 2);
        det.insert([10.0, 0.0, 0.0], 3);
        det.rebuild();
        let result = det.query_radius(&[0.0; 3], 1.0);
        assert!(result.contains(&1));
        assert!(result.contains(&2));
        assert!(!result.contains(&3));
    }

    #[test]
    fn test_detector_rebuild_chaining() {
        let mut det = KdTreeCollisionDetector::new(4);
        det.insert([0.0; 3], 0);
        det.rebuild();
        let r = det.query_radius(&[0.0; 3], 0.1);
        assert_eq!(r, vec![0]);
    }

    // ── BvhLeafGroups ───────────────────────────────────────────────────────

    #[test]
    fn test_bvh_leaf_groups_from_empty_tree() {
        let tree = KdTree::build(vec![]);
        let groups = BvhLeafGroups::from_tree(&tree);
        assert!(groups.groups.is_empty());
    }

    #[test]
    fn test_bvh_leaf_groups_non_empty() {
        let pts = grid_points(32);
        let tree = KdTree::build(pts);
        let groups = BvhLeafGroups::from_tree(&tree);
        assert!(!groups.groups.is_empty());
    }

    #[test]
    fn test_bvh_leaf_groups_query_sphere() {
        let pts = grid_points(32);
        let tree = KdTree::build(pts);
        let groups = BvhLeafGroups::from_tree(&tree);
        let hit = groups.query_sphere(&[1.5, 1.5, 1.5], 2.0);
        assert!(!hit.is_empty());
    }

    #[test]
    fn test_bvh_leaf_query_sphere_far_away() {
        let pts = grid_points(32);
        let tree = KdTree::build(pts);
        let groups = BvhLeafGroups::from_tree(&tree);
        // Query far from all points.
        let hit = groups.query_sphere(&[1000.0, 1000.0, 1000.0], 0.1);
        assert!(hit.is_empty());
    }
}
