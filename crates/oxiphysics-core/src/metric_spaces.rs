#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Metric space abstractions and concrete metric implementations.
//!
//! This module provides the [`MetricSpace`] trait and several concrete metrics:
//! Euclidean (L2), Manhattan (L1), Chebyshev (L∞), Hamming, and edit distance.
//! Also includes a [`MetricBallTree`] for efficient nearest-neighbor queries
//! and [`FrechetDistance`] for comparing polygonal curves.

#![allow(dead_code)]

use std::cmp::Ordering;
use std::collections::BinaryHeap;

// ---------------------------------------------------------------------------
// Core trait
// ---------------------------------------------------------------------------

/// A trait for types that form a metric space.
///
/// A metric space is a set `X` together with a function `d: X × X → ℝ≥0`
/// satisfying:
/// 1. `d(x, y) = 0` iff `x == y`  (identity of indiscernibles)
/// 2. `d(x, y) = d(y, x)`          (symmetry)
/// 3. `d(x, z) ≤ d(x, y) + d(y, z)` (triangle inequality)
pub trait MetricSpace {
    /// The type of points in the metric space.
    type Point;

    /// Compute the distance between two points.
    fn distance(&self, a: &Self::Point, b: &Self::Point) -> f64;

    /// Return all points from `candidates` within distance `radius` of `center`.
    fn ball<'a>(
        &self,
        center: &Self::Point,
        radius: f64,
        candidates: &'a [Self::Point],
    ) -> Vec<&'a Self::Point> {
        candidates
            .iter()
            .filter(|p| self.distance(center, p) <= radius)
            .collect()
    }

    /// Determine whether the diameter of a finite point set is bounded.
    ///
    /// Returns `true` when every pair of points has finite distance and the
    /// diameter (supremum of pairwise distances) is finite.
    fn is_bounded(&self, points: &[Self::Point]) -> bool {
        if points.is_empty() {
            return true;
        }
        let d = self.diameter(points);
        d.is_finite()
    }

    /// Compute the diameter of a finite point set, i.e. `sup_{x,y} d(x,y)`.
    ///
    /// Returns `0.0` for an empty or single-point set.
    fn diameter(&self, points: &[Self::Point]) -> f64 {
        if points.len() < 2 {
            return 0.0;
        }
        let mut max_d = 0.0_f64;
        for i in 0..points.len() {
            for j in (i + 1)..points.len() {
                let d = self.distance(&points[i], &points[j]);
                if d > max_d {
                    max_d = d;
                }
            }
        }
        max_d
    }
}

// ---------------------------------------------------------------------------
// Euclidean (L2) metric
// ---------------------------------------------------------------------------

/// Euclidean (L2) metric on ℝⁿ.
///
/// `d(x, y) = sqrt(Σ (xᵢ − yᵢ)²)`
///
/// The Cauchy–Schwarz inequality states `|⟨x,y⟩| ≤ ‖x‖₂ · ‖y‖₂`, which
/// underlies many properties of this metric.
#[derive(Debug, Clone, Copy, Default)]
pub struct EuclideanMetric;

impl MetricSpace for EuclideanMetric {
    type Point = Vec<f64>;

    fn distance(&self, a: &Vec<f64>, b: &Vec<f64>) -> f64 {
        assert_eq!(a.len(), b.len(), "EuclideanMetric: dimension mismatch");
        a.iter()
            .zip(b.iter())
            .map(|(x, y)| (x - y).powi(2))
            .sum::<f64>()
            .sqrt()
    }
}

impl EuclideanMetric {
    /// Compute the L2 norm (Euclidean length) of a vector.
    ///
    /// `‖v‖₂ = sqrt(Σ vᵢ²)`
    pub fn norm(v: &[f64]) -> f64 {
        v.iter().map(|x| x.powi(2)).sum::<f64>().sqrt()
    }

    /// Dot product of two equal-length vectors.
    pub fn dot(a: &[f64], b: &[f64]) -> f64 {
        assert_eq!(a.len(), b.len(), "dot: dimension mismatch");
        a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
    }

    /// Verify the Cauchy–Schwarz inequality: `|⟨a,b⟩| ≤ ‖a‖₂ · ‖b‖₂`.
    ///
    /// Returns `true` when the inequality holds (within floating-point tolerance).
    pub fn cauchy_schwarz_holds(a: &[f64], b: &[f64]) -> bool {
        let dot = Self::dot(a, b).abs();
        let product = Self::norm(a) * Self::norm(b);
        dot <= product + 1e-9
    }

    /// Project vector `v` onto unit vector `u`.
    ///
    /// Returns the scalar projection `⟨v, u⟩`.
    pub fn scalar_projection(v: &[f64], u: &[f64]) -> f64 {
        Self::dot(v, u) / Self::norm(u)
    }

    /// Normalize a vector to unit length. Returns `None` for the zero vector.
    pub fn normalize(v: &[f64]) -> Option<Vec<f64>> {
        let n = Self::norm(v);
        if n < 1e-15 {
            return None;
        }
        Some(v.iter().map(|x| x / n).collect())
    }
}

// ---------------------------------------------------------------------------
// Manhattan (L1) metric
// ---------------------------------------------------------------------------

/// Manhattan (taxicab / L1) metric on ℝⁿ.
///
/// `d(x, y) = Σ |xᵢ − yᵢ|`
#[derive(Debug, Clone, Copy, Default)]
pub struct ManhattanMetric;

impl MetricSpace for ManhattanMetric {
    type Point = Vec<f64>;

    fn distance(&self, a: &Vec<f64>, b: &Vec<f64>) -> f64 {
        assert_eq!(a.len(), b.len(), "ManhattanMetric: dimension mismatch");
        a.iter().zip(b.iter()).map(|(x, y)| (x - y).abs()).sum()
    }
}

impl ManhattanMetric {
    /// Compute the L1 norm of a vector.
    pub fn norm(v: &[f64]) -> f64 {
        v.iter().map(|x| x.abs()).sum()
    }

    /// The L1 unit ball is the cross-polytope; check membership.
    pub fn in_unit_ball(v: &[f64]) -> bool {
        Self::norm(v) <= 1.0 + 1e-12
    }
}

// ---------------------------------------------------------------------------
// Chebyshev (L∞) metric
// ---------------------------------------------------------------------------

/// Chebyshev (chessboard / L∞) metric on ℝⁿ.
///
/// `d(x, y) = max_i |xᵢ − yᵢ|`
#[derive(Debug, Clone, Copy, Default)]
pub struct ChebyshevMetric;

impl MetricSpace for ChebyshevMetric {
    type Point = Vec<f64>;

    fn distance(&self, a: &Vec<f64>, b: &Vec<f64>) -> f64 {
        assert_eq!(a.len(), b.len(), "ChebyshevMetric: dimension mismatch");
        a.iter()
            .zip(b.iter())
            .map(|(x, y)| (x - y).abs())
            .fold(0.0_f64, f64::max)
    }
}

impl ChebyshevMetric {
    /// Compute the L∞ norm of a vector.
    pub fn norm(v: &[f64]) -> f64 {
        v.iter().map(|x| x.abs()).fold(0.0_f64, f64::max)
    }
}

// ---------------------------------------------------------------------------
// Hamming metric
// ---------------------------------------------------------------------------

/// Hamming metric on equal-length sequences.
///
/// `d(x, y) = |{ i : xᵢ ≠ yᵢ }|`
///
/// Works with any element type that implements `PartialEq`.
#[derive(Debug, Clone, Copy, Default)]
pub struct HammingMetric;

/// A fixed-length sequence for the Hamming metric.
#[derive(Debug, Clone, PartialEq)]
pub struct HammingPoint(pub Vec<u8>);

impl MetricSpace for HammingMetric {
    type Point = HammingPoint;

    fn distance(&self, a: &HammingPoint, b: &HammingPoint) -> f64 {
        assert_eq!(
            a.0.len(),
            b.0.len(),
            "HammingMetric: length mismatch ({} vs {})",
            a.0.len(),
            b.0.len()
        );
        a.0.iter().zip(b.0.iter()).filter(|(x, y)| x != y).count() as f64
    }
}

impl HammingMetric {
    /// Compute the Hamming distance between two byte slices of equal length.
    pub fn hamming_distance(a: &[u8], b: &[u8]) -> usize {
        assert_eq!(a.len(), b.len(), "hamming_distance: length mismatch");
        a.iter().zip(b.iter()).filter(|(x, y)| x != y).count()
    }

    /// Compute the Hamming weight (number of 1 bits) of a byte slice.
    pub fn hamming_weight(v: &[u8]) -> usize {
        v.iter().map(|b| b.count_ones() as usize).sum()
    }
}

// ---------------------------------------------------------------------------
// Edit (Levenshtein) distance
// ---------------------------------------------------------------------------

/// Edit (Levenshtein) distance metric on sequences.
///
/// The edit distance counts the minimum number of single-element insertions,
/// deletions, and substitutions required to transform one sequence into another.
#[derive(Debug, Clone, Copy, Default)]
pub struct EditDistance;

/// A sequence point for edit-distance computation.
#[derive(Debug, Clone, PartialEq)]
pub struct Sequence(pub Vec<char>);

impl MetricSpace for EditDistance {
    type Point = Sequence;

    fn distance(&self, a: &Sequence, b: &Sequence) -> f64 {
        levenshtein(&a.0, &b.0) as f64
    }
}

/// Compute the Levenshtein distance between two character slices.
///
/// Uses the classic dynamic-programming algorithm with O(min(m,n)) space.
pub fn levenshtein(a: &[char], b: &[char]) -> usize {
    let m = a.len();
    let n = b.len();
    if m == 0 {
        return n;
    }
    if n == 0 {
        return m;
    }

    let mut prev: Vec<usize> = (0..=n).collect();
    let mut curr = vec![0usize; n + 1];

    for i in 1..=m {
        curr[0] = i;
        for j in 1..=n {
            let cost = if a[i - 1] == b[j - 1] { 0 } else { 1 };
            curr[j] = (prev[j] + 1).min(curr[j - 1] + 1).min(prev[j - 1] + cost);
        }
        std::mem::swap(&mut prev, &mut curr);
    }
    prev[n]
}

/// Compute the Levenshtein distance between two strings.
pub fn levenshtein_str(a: &str, b: &str) -> usize {
    let ac: Vec<char> = a.chars().collect();
    let bc: Vec<char> = b.chars().collect();
    levenshtein(&ac, &bc)
}

impl EditDistance {
    /// Compute the longest common subsequence length of two sequences.
    pub fn lcs_length(a: &[char], b: &[char]) -> usize {
        let m = a.len();
        let n = b.len();
        let mut dp = vec![vec![0usize; n + 1]; m + 1];
        for i in 1..=m {
            for j in 1..=n {
                dp[i][j] = if a[i - 1] == b[j - 1] {
                    dp[i - 1][j - 1] + 1
                } else {
                    dp[i - 1][j].max(dp[i][j - 1])
                };
            }
        }
        dp[m][n]
    }

    /// Compute the edit distance with custom costs for insert, delete, and substitute.
    #[allow(clippy::too_many_arguments)]
    pub fn weighted_edit(
        a: &[char],
        b: &[char],
        insert_cost: f64,
        delete_cost: f64,
        subst_cost: f64,
    ) -> f64 {
        let m = a.len();
        let n = b.len();
        if m == 0 {
            return n as f64 * insert_cost;
        }
        if n == 0 {
            return m as f64 * delete_cost;
        }
        let mut prev: Vec<f64> = (0..=n).map(|j| j as f64 * insert_cost).collect();
        let mut curr = vec![0.0_f64; n + 1];
        for i in 1..=m {
            curr[0] = i as f64 * delete_cost;
            for j in 1..=n {
                let cost = if a[i - 1] == b[j - 1] {
                    0.0
                } else {
                    subst_cost
                };
                curr[j] = (prev[j] + delete_cost)
                    .min(curr[j - 1] + insert_cost)
                    .min(prev[j - 1] + cost);
            }
            std::mem::swap(&mut prev, &mut curr);
        }
        prev[n]
    }
}

// ---------------------------------------------------------------------------
// Minkowski metric (generalization)
// ---------------------------------------------------------------------------

/// Minkowski metric of order `p` on ℝⁿ.
///
/// `d_p(x, y) = (Σ |xᵢ − yᵢ|^p)^(1/p)`
///
/// Special cases: `p = 1` → Manhattan, `p = 2` → Euclidean.
/// For `p → ∞` this converges to the Chebyshev metric.
#[derive(Debug, Clone, Copy)]
pub struct MinkowskiMetric {
    /// The order parameter `p ≥ 1`.
    pub p: f64,
}

impl MinkowskiMetric {
    /// Create a Minkowski metric of order `p`.
    ///
    /// # Panics
    /// Panics if `p < 1.0`.
    pub fn new(p: f64) -> Self {
        assert!(p >= 1.0, "MinkowskiMetric: p must be >= 1, got {p}");
        Self { p }
    }
}

impl MetricSpace for MinkowskiMetric {
    type Point = Vec<f64>;

    fn distance(&self, a: &Vec<f64>, b: &Vec<f64>) -> f64 {
        assert_eq!(a.len(), b.len(), "MinkowskiMetric: dimension mismatch");
        a.iter()
            .zip(b.iter())
            .map(|(x, y)| (x - y).abs().powf(self.p))
            .sum::<f64>()
            .powf(1.0 / self.p)
    }
}

// ---------------------------------------------------------------------------
// Cosine distance (not a true metric but common)
// ---------------------------------------------------------------------------

/// Cosine similarity-based pseudo-distance on ℝⁿ.
///
/// `d_cos(x, y) = 1 − cos(θ) = 1 − (⟨x,y⟩ / (‖x‖₂ · ‖y‖₂))`
///
/// Note: this is not a true metric (triangle inequality may fail), but it is
/// widely used in information retrieval and machine-learning contexts.
#[derive(Debug, Clone, Copy, Default)]
pub struct CosineDistance;

impl CosineDistance {
    /// Compute cosine similarity between two vectors.
    pub fn similarity(a: &[f64], b: &[f64]) -> f64 {
        let dot: f64 = a.iter().zip(b.iter()).map(|(x, y)| x * y).sum();
        let na = EuclideanMetric::norm(a);
        let nb = EuclideanMetric::norm(b);
        if na < 1e-15 || nb < 1e-15 {
            return 0.0;
        }
        (dot / (na * nb)).clamp(-1.0, 1.0)
    }

    /// Compute cosine distance: `1 − similarity(a, b)`.
    pub fn distance(a: &[f64], b: &[f64]) -> f64 {
        1.0 - Self::similarity(a, b)
    }
}

// ---------------------------------------------------------------------------
// Ball Tree for nearest-neighbor queries
// ---------------------------------------------------------------------------

/// A node in a [`MetricBallTree`].
#[derive(Debug, Clone)]
struct BallTreeNode {
    /// Index of the pivot point in the original dataset.
    pivot_idx: usize,
    /// Radius of the ball enclosing all points in this subtree.
    radius: f64,
    /// Left child index (or `usize::MAX` for leaf).
    left: usize,
    /// Right child index (or `usize::MAX` for leaf).
    right: usize,
    /// Indices of points stored at this leaf (non-empty only for leaves).
    leaf_points: Vec<usize>,
}

impl BallTreeNode {
    fn is_leaf(&self) -> bool {
        self.left == usize::MAX && self.right == usize::MAX
    }
}

/// Ball tree (metric tree) for efficient nearest-neighbor and range queries.
///
/// Supports any [`MetricSpace`] whose `Point` type is `Clone`. The tree is
/// built once from a dataset and then supports O(log n) average-case queries.
pub struct MetricBallTree {
    points: Vec<Vec<f64>>,
    nodes: Vec<BallTreeNode>,
    root: usize,
    leaf_size: usize,
}

/// A nearest-neighbor search result.
#[derive(Debug, Clone)]
pub struct NearestNeighbor {
    /// Index of the neighbor in the original dataset.
    pub index: usize,
    /// Distance from the query point to this neighbor.
    pub distance: f64,
}

impl PartialEq for NearestNeighbor {
    fn eq(&self, other: &Self) -> bool {
        self.distance == other.distance
    }
}

impl Eq for NearestNeighbor {}

impl PartialOrd for NearestNeighbor {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for NearestNeighbor {
    fn cmp(&self, other: &Self) -> Ordering {
        // Max-heap by distance (so smallest distances bubble up as we invert)
        other
            .distance
            .partial_cmp(&self.distance)
            .unwrap_or(Ordering::Equal)
    }
}

impl MetricBallTree {
    /// Build a ball tree from a set of points using the Euclidean metric.
    ///
    /// `leaf_size` controls the maximum number of points stored at a leaf node.
    pub fn build(points: Vec<Vec<f64>>, leaf_size: usize) -> Self {
        let n = points.len();
        assert!(!points.is_empty(), "MetricBallTree: empty dataset");
        let leaf_size = leaf_size.max(1);

        let mut nodes: Vec<BallTreeNode> = Vec::with_capacity(2 * n / leaf_size + 4);
        let indices: Vec<usize> = (0..n).collect();
        let root = build_node(&points, &indices, &mut nodes, leaf_size);

        Self {
            points,
            nodes,
            root,
            leaf_size,
        }
    }

    /// Find the `k` nearest neighbors of `query` using the Euclidean metric.
    pub fn knn(&self, query: &[f64], k: usize) -> Vec<NearestNeighbor> {
        let mut heap: BinaryHeap<NearestNeighbor> = BinaryHeap::new();
        let metric = EuclideanMetric;
        knn_search(
            &self.points,
            &self.nodes,
            self.root,
            query,
            k,
            &mut heap,
            &metric,
        );
        let mut result: Vec<NearestNeighbor> = heap.into_sorted_vec();
        result.sort_by(|a, b| {
            a.distance
                .partial_cmp(&b.distance)
                .unwrap_or(Ordering::Equal)
        });
        result
    }

    /// Find all points within `radius` of `query`.
    pub fn range_query(&self, query: &[f64], radius: f64) -> Vec<NearestNeighbor> {
        let metric = EuclideanMetric;
        let mut result = Vec::new();
        range_search(
            &self.points,
            &self.nodes,
            self.root,
            query,
            radius,
            &mut result,
            &metric,
        );
        result.sort_by(|a, b| {
            a.distance
                .partial_cmp(&b.distance)
                .unwrap_or(Ordering::Equal)
        });
        result
    }

    /// Return the number of points in the tree.
    pub fn len(&self) -> usize {
        self.points.len()
    }

    /// Return `true` if the tree has no points.
    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }

    /// Return the configured leaf size.
    pub fn leaf_size(&self) -> usize {
        self.leaf_size
    }
}

// --- Ball tree construction helpers ---

fn euclidean_dist(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b.iter())
        .map(|(x, y)| (x - y).powi(2))
        .sum::<f64>()
        .sqrt()
}

fn centroid(points: &[Vec<f64>], indices: &[usize]) -> Vec<f64> {
    let dim = points[indices[0]].len();
    let mut c = vec![0.0; dim];
    for &idx in indices {
        for (d, &v) in c.iter_mut().zip(points[idx].iter()) {
            *d += v;
        }
    }
    let n = indices.len() as f64;
    c.iter_mut().for_each(|v| *v /= n);
    c
}

fn build_node(
    points: &[Vec<f64>],
    indices: &[usize],
    nodes: &mut Vec<BallTreeNode>,
    leaf_size: usize,
) -> usize {
    let c = centroid(points, indices);
    // Choose pivot as point farthest from centroid
    let pivot_idx = indices
        .iter()
        .copied()
        .max_by(|&a, &b| {
            euclidean_dist(&points[a], &c)
                .partial_cmp(&euclidean_dist(&points[b], &c))
                .unwrap_or(Ordering::Equal)
        })
        .expect("indices is non-empty");

    let radius = indices
        .iter()
        .map(|&i| euclidean_dist(&points[i], &points[pivot_idx]))
        .fold(0.0_f64, f64::max);

    let node_idx = nodes.len();

    if indices.len() <= leaf_size {
        nodes.push(BallTreeNode {
            pivot_idx,
            radius,
            left: usize::MAX,
            right: usize::MAX,
            leaf_points: indices.to_vec(),
        });
        return node_idx;
    }

    // Split by distance to pivot
    let mut left_idx: Vec<usize> = Vec::new();
    let mut right_idx: Vec<usize> = Vec::new();
    for &i in indices {
        if i == pivot_idx {
            left_idx.push(i);
            continue;
        }
        if left_idx.len() <= right_idx.len() {
            left_idx.push(i);
        } else {
            right_idx.push(i);
        }
    }
    if left_idx.is_empty() {
        left_idx.push(pivot_idx);
    }
    if right_idx.is_empty() {
        right_idx = left_idx.split_off(left_idx.len() / 2 + 1);
    }

    // Push a placeholder so we can fix up children indices later
    nodes.push(BallTreeNode {
        pivot_idx,
        radius,
        left: usize::MAX,
        right: usize::MAX,
        leaf_points: vec![],
    });

    let left_child = build_node(points, &left_idx, nodes, leaf_size);
    let right_child = build_node(points, &right_idx, nodes, leaf_size);
    nodes[node_idx].left = left_child;
    nodes[node_idx].right = right_child;

    node_idx
}

fn knn_search(
    points: &[Vec<f64>],
    nodes: &[BallTreeNode],
    node_idx: usize,
    query: &[f64],
    k: usize,
    heap: &mut BinaryHeap<NearestNeighbor>,
    _metric: &EuclideanMetric,
) {
    let node = &nodes[node_idx];
    let pivot_dist = euclidean_dist(query, &points[node.pivot_idx]);

    // Pruning: if the closest possible point in this ball is farther than
    // the current k-th best, skip.
    if heap.len() >= k {
        let worst = heap.peek().map(|n| n.distance).unwrap_or(f64::MAX);
        if pivot_dist - node.radius > worst {
            return;
        }
    }

    if node.is_leaf() {
        for &idx in &node.leaf_points {
            let d = euclidean_dist(query, &points[idx]);
            if heap.len() < k {
                heap.push(NearestNeighbor {
                    index: idx,
                    distance: d,
                });
            } else if let Some(worst) = heap.peek()
                && d < worst.distance
            {
                heap.pop();
                heap.push(NearestNeighbor {
                    index: idx,
                    distance: d,
                });
            }
        }
        return;
    }

    // Recurse into closer child first
    let left_dist = if node.left != usize::MAX {
        euclidean_dist(query, &points[nodes[node.left].pivot_idx])
    } else {
        f64::MAX
    };
    let right_dist = if node.right != usize::MAX {
        euclidean_dist(query, &points[nodes[node.right].pivot_idx])
    } else {
        f64::MAX
    };

    if left_dist <= right_dist {
        if node.left != usize::MAX {
            knn_search(points, nodes, node.left, query, k, heap, _metric);
        }
        if node.right != usize::MAX {
            knn_search(points, nodes, node.right, query, k, heap, _metric);
        }
    } else {
        if node.right != usize::MAX {
            knn_search(points, nodes, node.right, query, k, heap, _metric);
        }
        if node.left != usize::MAX {
            knn_search(points, nodes, node.left, query, k, heap, _metric);
        }
    }
}

fn range_search(
    points: &[Vec<f64>],
    nodes: &[BallTreeNode],
    node_idx: usize,
    query: &[f64],
    radius: f64,
    result: &mut Vec<NearestNeighbor>,
    _metric: &EuclideanMetric,
) {
    let node = &nodes[node_idx];
    let pivot_dist = euclidean_dist(query, &points[node.pivot_idx]);

    // Pruning: if the closest possible point in the ball is outside `radius`, skip.
    if pivot_dist - node.radius > radius {
        return;
    }

    if node.is_leaf() {
        for &idx in &node.leaf_points {
            let d = euclidean_dist(query, &points[idx]);
            if d <= radius {
                result.push(NearestNeighbor {
                    index: idx,
                    distance: d,
                });
            }
        }
        return;
    }

    if node.left != usize::MAX {
        range_search(points, nodes, node.left, query, radius, result, _metric);
    }
    if node.right != usize::MAX {
        range_search(points, nodes, node.right, query, radius, result, _metric);
    }
}

// ---------------------------------------------------------------------------
// Discrete Fréchet distance
// ---------------------------------------------------------------------------

/// Discrete Fréchet distance between two polygonal curves.
///
/// The Fréchet distance is informally described as the minimum leash length
/// required for a person walking a dog, where person and dog each traverse
/// their respective curve from start to finish without backtracking.
///
/// This implementation computes the *discrete* variant using dynamic programming
/// in O(mn) time and space.
pub struct FrechetDistance;

impl FrechetDistance {
    /// Compute the discrete Fréchet distance between curves `p` and `q`.
    ///
    /// Each curve is given as a slice of points in ℝⁿ (represented as `Vec`f64`).
    /// Returns the minimum coupling distance.
    pub fn compute(p: &[Vec<f64>], q: &[Vec<f64>]) -> f64 {
        let m = p.len();
        let n = q.len();
        if m == 0 || n == 0 {
            return 0.0;
        }

        let mut dp = vec![vec![f64::MAX; n]; m];

        for i in 0..m {
            for j in 0..n {
                let d = euclidean_dist(&p[i], &q[j]);
                dp[i][j] = if i == 0 && j == 0 {
                    d
                } else if i == 0 {
                    dp[i][j - 1].max(d)
                } else if j == 0 {
                    dp[i - 1][j].max(d)
                } else {
                    dp[i - 1][j].min(dp[i][j - 1]).min(dp[i - 1][j - 1]).max(d)
                };
            }
        }

        dp[m - 1][n - 1]
    }

    /// Compute the discrete Fréchet distance using a custom metric.
    pub fn compute_with_metric<M: MetricSpace<Point = Vec<f64>>>(
        p: &[Vec<f64>],
        q: &[Vec<f64>],
        metric: &M,
    ) -> f64 {
        let m = p.len();
        let n = q.len();
        if m == 0 || n == 0 {
            return 0.0;
        }

        let mut dp = vec![vec![f64::MAX; n]; m];

        for i in 0..m {
            for j in 0..n {
                let d = metric.distance(&p[i], &q[j]);
                dp[i][j] = if i == 0 && j == 0 {
                    d
                } else if i == 0 {
                    dp[i][j - 1].max(d)
                } else if j == 0 {
                    dp[i - 1][j].max(d)
                } else {
                    dp[i - 1][j].min(dp[i][j - 1]).min(dp[i - 1][j - 1]).max(d)
                };
            }
        }

        dp[m - 1][n - 1]
    }

    /// Compute the Hausdorff distance between two finite point sets.
    ///
    /// `d_H(A, B) = max(sup_{a∈A} inf_{b∈B} d(a,b), sup_{b∈B} inf_{a∈A} d(a,b))`
    pub fn hausdorff(p: &[Vec<f64>], q: &[Vec<f64>]) -> f64 {
        let one_sided = |a: &[Vec<f64>], b: &[Vec<f64>]| {
            a.iter()
                .map(|ai| {
                    b.iter()
                        .map(|bj| euclidean_dist(ai, bj))
                        .fold(f64::MAX, f64::min)
                })
                .fold(0.0_f64, f64::max)
        };
        one_sided(p, q).max(one_sided(q, p))
    }
}

// ---------------------------------------------------------------------------
// Geodesic distance on a graph
// ---------------------------------------------------------------------------

/// Graph-based geodesic metric using Dijkstra's shortest path.
///
/// Points are vertices (indexed by `usize`), and edges carry non-negative weights.
pub struct GeodesicMetric {
    /// Adjacency list: `adj\[i\]` = list of `(neighbor, weight)` pairs.
    pub adj: Vec<Vec<(usize, f64)>>,
}

impl GeodesicMetric {
    /// Create a new geodesic metric for `n` vertices.
    pub fn new(n: usize) -> Self {
        Self {
            adj: vec![Vec::new(); n],
        }
    }

    /// Add an undirected edge with the given weight.
    pub fn add_edge(&mut self, u: usize, v: usize, weight: f64) {
        self.adj[u].push((v, weight));
        self.adj[v].push((u, weight));
    }

    /// Compute shortest distances from source `s` to all vertices (Dijkstra).
    pub fn dijkstra(&self, s: usize) -> Vec<f64> {
        let n = self.adj.len();
        let mut dist = vec![f64::INFINITY; n];
        dist[s] = 0.0;

        // (distance, vertex) — use a min-heap via Reverse
        let mut heap: BinaryHeap<(ordered_float::OrderedF64, usize)> = BinaryHeap::new();
        heap.push((ordered_float::OrderedF64(0.0), s));

        while let Some((ordered_float::OrderedF64(d), u)) = heap.pop() {
            if d > dist[u] {
                continue;
            }
            for &(v, w) in &self.adj[u] {
                let nd = dist[u] + w;
                if nd < dist[v] {
                    dist[v] = nd;
                    heap.push((ordered_float::OrderedF64(nd), v));
                }
            }
        }
        dist
    }
}

/// Wrapper for `f64` that implements `Ord` for use in `BinaryHeap` with negated distances.
mod ordered_float {
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct OrderedF64(pub f64);

    impl Eq for OrderedF64 {}

    impl PartialOrd for OrderedF64 {
        fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
            Some(self.cmp(other))
        }
    }

    impl Ord for OrderedF64 {
        fn cmp(&self, other: &Self) -> std::cmp::Ordering {
            // Min-heap: smaller f64 → larger in reversed comparison
            other
                .0
                .partial_cmp(&self.0)
                .unwrap_or(std::cmp::Ordering::Equal)
        }
    }
}

impl MetricSpace for GeodesicMetric {
    type Point = usize; // vertices

    fn distance(&self, a: &usize, b: &usize) -> f64 {
        let dists = self.dijkstra(*a);
        dists[*b]
    }
}

// ---------------------------------------------------------------------------
// Metric space utilities
// ---------------------------------------------------------------------------

/// Check whether a set of pairwise distances satisfies the triangle inequality.
///
/// `matrix\[i\]\[j\]` should equal `d(i, j)`. Returns `Ok(())` on success or
/// `Err((i, j, k))` with the first violated triple.
pub fn check_triangle_inequality(matrix: &[Vec<f64>]) -> Result<(), (usize, usize, usize)> {
    let n = matrix.len();
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                if matrix[i][k] > matrix[i][j] + matrix[j][k] + 1e-9 {
                    return Err((i, j, k));
                }
            }
        }
    }
    Ok(())
}

/// Compute the diameter of a point set under a given metric.
///
/// This is a convenience wrapper around [`MetricSpace::diameter`].
pub fn set_diameter<M: MetricSpace>(metric: &M, points: &[M::Point]) -> f64 {
    metric.diameter(points)
}

/// Compute the pairwise distance matrix for a set of points.
pub fn distance_matrix<M: MetricSpace>(metric: &M, points: &[M::Point]) -> Vec<Vec<f64>> {
    let n = points.len();
    let mut mat = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            mat[i][j] = metric.distance(&points[i], &points[j]);
        }
    }
    mat
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- EuclideanMetric ---

    #[test]
    fn test_euclidean_zero_distance() {
        let m = EuclideanMetric;
        let p = vec![1.0, 2.0, 3.0];
        assert!((m.distance(&p, &p) - 0.0).abs() < 1e-12);
    }

    #[test]
    fn test_euclidean_known_distance() {
        let m = EuclideanMetric;
        let a = vec![0.0, 0.0];
        let b = vec![3.0, 4.0];
        assert!((m.distance(&a, &b) - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_euclidean_symmetry() {
        let m = EuclideanMetric;
        let a = vec![1.0, 2.0, 3.0];
        let b = vec![4.0, -1.0, 0.0];
        let diff = (m.distance(&a, &b) - m.distance(&b, &a)).abs();
        assert!(diff < 1e-12);
    }

    #[test]
    fn test_euclidean_triangle_inequality() {
        let m = EuclideanMetric;
        let a = vec![0.0, 0.0];
        let b = vec![1.0, 0.0];
        let c = vec![0.5, 1.0];
        let dab = m.distance(&a, &b);
        let dbc = m.distance(&b, &c);
        let dac = m.distance(&a, &c);
        assert!(dac <= dab + dbc + 1e-10);
    }

    #[test]
    fn test_euclidean_norm() {
        let v = vec![3.0, 4.0];
        assert!((EuclideanMetric::norm(&v) - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_cauchy_schwarz() {
        let a = vec![1.0, 2.0, 3.0];
        let b = vec![4.0, -1.0, 2.0];
        assert!(EuclideanMetric::cauchy_schwarz_holds(&a, &b));
    }

    #[test]
    fn test_euclidean_normalize() {
        let v = vec![3.0, 4.0];
        let u = EuclideanMetric::normalize(&v).unwrap();
        assert!((EuclideanMetric::norm(&u) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_euclidean_diameter() {
        let m = EuclideanMetric;
        let pts = vec![vec![0.0, 0.0], vec![1.0, 0.0], vec![0.0, 1.0]];
        let d = m.diameter(&pts);
        assert!((d - 2.0_f64.sqrt()).abs() < 1e-10);
    }

    #[test]
    fn test_euclidean_ball() {
        let m = EuclideanMetric;
        let center = vec![0.0, 0.0];
        let pts = vec![vec![0.5, 0.0], vec![2.0, 0.0], vec![0.0, 0.8]];
        let ball = m.ball(&center, 1.0, &pts);
        assert_eq!(ball.len(), 2);
    }

    #[test]
    fn test_euclidean_is_bounded() {
        let m = EuclideanMetric;
        let pts = vec![vec![0.0, 0.0], vec![100.0, 100.0]];
        assert!(m.is_bounded(&pts));
    }

    // --- ManhattanMetric ---

    #[test]
    fn test_manhattan_distance() {
        let m = ManhattanMetric;
        let a = vec![0.0, 0.0];
        let b = vec![3.0, 4.0];
        assert!((m.distance(&a, &b) - 7.0).abs() < 1e-12);
    }

    #[test]
    fn test_manhattan_symmetry() {
        let m = ManhattanMetric;
        let a = vec![1.0, -2.0, 3.0];
        let b = vec![-1.0, 4.0, 0.0];
        assert!((m.distance(&a, &b) - m.distance(&b, &a)).abs() < 1e-12);
    }

    #[test]
    fn test_manhattan_norm() {
        let v = vec![-3.0, 4.0];
        assert!((ManhattanMetric::norm(&v) - 7.0).abs() < 1e-12);
    }

    #[test]
    fn test_manhattan_diameter() {
        let m = ManhattanMetric;
        let pts = vec![vec![0.0, 0.0], vec![1.0, 1.0], vec![-1.0, -1.0]];
        let d = m.diameter(&pts);
        assert!((d - 4.0).abs() < 1e-12);
    }

    #[test]
    fn test_manhattan_ball() {
        let m = ManhattanMetric;
        let center = vec![0.0, 0.0];
        let pts = vec![vec![1.0, 0.0], vec![1.0, 1.0], vec![0.5, 0.5]];
        let b = m.ball(&center, 1.0, &pts);
        // d([1,0]) = 1 ≤ 1, d([1,1]) = 2 > 1, d([0.5,0.5]) = 1 ≤ 1
        assert_eq!(b.len(), 2);
    }

    // --- ChebyshevMetric ---

    #[test]
    fn test_chebyshev_distance() {
        let m = ChebyshevMetric;
        let a = vec![1.0, 3.0, 5.0];
        let b = vec![4.0, 2.0, 1.0];
        // diffs: 3, 1, 4 → max = 4
        assert!((m.distance(&a, &b) - 4.0).abs() < 1e-12);
    }

    #[test]
    fn test_chebyshev_zero() {
        let m = ChebyshevMetric;
        let p = vec![1.0, 2.0];
        assert!((m.distance(&p, &p)).abs() < 1e-12);
    }

    #[test]
    fn test_chebyshev_norm() {
        let v = vec![-5.0, 3.0, 2.0];
        assert!((ChebyshevMetric::norm(&v) - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_chebyshev_diameter() {
        let m = ChebyshevMetric;
        let pts = vec![vec![0.0, 0.0], vec![2.0, 1.0], vec![1.0, 3.0]];
        let d = m.diameter(&pts);
        // pairs: d([0,0],[2,1])=2, d([0,0],[1,3])=3, d([2,1],[1,3])=2 → max=3
        assert!((d - 3.0).abs() < 1e-12);
    }

    // --- HammingMetric ---

    #[test]
    fn test_hamming_equal() {
        let m = HammingMetric;
        let a = HammingPoint(vec![1, 0, 1, 1]);
        assert!((m.distance(&a, &a)).abs() < 1e-12);
    }

    #[test]
    fn test_hamming_known() {
        let m = HammingMetric;
        let a = HammingPoint(vec![0, 0, 0, 1]);
        let b = HammingPoint(vec![1, 1, 0, 0]);
        assert!((m.distance(&a, &b) - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_hamming_symmetry() {
        let m = HammingMetric;
        let a = HammingPoint(vec![1, 0, 1]);
        let b = HammingPoint(vec![0, 1, 1]);
        assert!((m.distance(&a, &b) - m.distance(&b, &a)).abs() < 1e-12);
    }

    #[test]
    fn test_hamming_weight() {
        let v = vec![0b1010_1010u8, 0b1111_0000u8];
        // 4 ones + 4 ones = 8
        assert_eq!(HammingMetric::hamming_weight(&v), 8);
    }

    // --- EditDistance ---

    #[test]
    fn test_edit_identical() {
        let m = EditDistance;
        let s: Vec<char> = "hello".chars().collect();
        let a = Sequence(s.clone());
        let b = Sequence(s);
        assert!((m.distance(&a, &b)).abs() < 1e-12);
    }

    #[test]
    fn test_edit_known() {
        assert_eq!(levenshtein_str("kitten", "sitting"), 3);
    }

    #[test]
    fn test_edit_empty() {
        assert_eq!(levenshtein_str("", "abc"), 3);
        assert_eq!(levenshtein_str("abc", ""), 3);
    }

    #[test]
    fn test_edit_symmetry() {
        let a: Vec<char> = "sunday".chars().collect();
        let b: Vec<char> = "saturday".chars().collect();
        assert_eq!(levenshtein(&a, &b), levenshtein(&b, &a));
    }

    #[test]
    fn test_lcs_length() {
        let a: Vec<char> = "ABCBDAB".chars().collect();
        let b: Vec<char> = "BDCAB".chars().collect();
        assert_eq!(EditDistance::lcs_length(&a, &b), 4);
    }

    #[test]
    fn test_weighted_edit() {
        let a: Vec<char> = "abc".chars().collect();
        let b: Vec<char> = "axc".chars().collect();
        // One substitution, cost 2.0
        let d = EditDistance::weighted_edit(&a, &b, 1.0, 1.0, 2.0);
        assert!((d - 2.0).abs() < 1e-12);
    }

    // --- MinkowskiMetric ---

    #[test]
    fn test_minkowski_p1_equals_manhattan() {
        let m1 = MinkowskiMetric::new(1.0);
        let m2 = ManhattanMetric;
        let a = vec![1.0, 2.0, 3.0];
        let b = vec![4.0, 0.0, -1.0];
        assert!((m1.distance(&a, &b) - m2.distance(&a, &b)).abs() < 1e-10);
    }

    #[test]
    fn test_minkowski_p2_equals_euclidean() {
        let m1 = MinkowskiMetric::new(2.0);
        let m2 = EuclideanMetric;
        let a = vec![0.0, 0.0];
        let b = vec![3.0, 4.0];
        assert!((m1.distance(&a, &b) - m2.distance(&a, &b)).abs() < 1e-10);
    }

    // --- CosineDistance ---

    #[test]
    fn test_cosine_identical() {
        let v = vec![1.0, 2.0, 3.0];
        assert!((CosineDistance::similarity(&v, &v) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_cosine_orthogonal() {
        let a = vec![1.0, 0.0];
        let b = vec![0.0, 1.0];
        assert!((CosineDistance::similarity(&a, &b)).abs() < 1e-12);
    }

    #[test]
    fn test_cosine_distance_range() {
        let a = vec![1.0, 0.0];
        let b = vec![-1.0, 0.0];
        let d = CosineDistance::distance(&a, &b);
        assert!((0.0..=2.0 + 1e-12).contains(&d));
    }

    // --- MetricBallTree ---

    #[test]
    fn test_ball_tree_knn_1() {
        let points: Vec<Vec<f64>> = vec![
            vec![0.0, 0.0],
            vec![1.0, 0.0],
            vec![0.0, 1.0],
            vec![5.0, 5.0],
        ];
        let tree = MetricBallTree::build(points, 2);
        let result = tree.knn(&[0.1, 0.1], 1);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].index, 0);
    }

    #[test]
    fn test_ball_tree_knn_k() {
        let points: Vec<Vec<f64>> = (0..20).map(|i| vec![i as f64, 0.0]).collect();
        let tree = MetricBallTree::build(points, 4);
        let result = tree.knn(&[9.5, 0.0], 3);
        assert_eq!(result.len(), 3);
        // Nearest should be 9 or 10
        let indices: Vec<usize> = result.iter().map(|n| n.index).collect();
        assert!(indices.contains(&9) || indices.contains(&10));
    }

    #[test]
    fn test_ball_tree_range_query() {
        let points: Vec<Vec<f64>> = vec![
            vec![0.0, 0.0],
            vec![0.5, 0.0],
            vec![1.5, 0.0],
            vec![10.0, 10.0],
        ];
        let tree = MetricBallTree::build(points, 2);
        let result = tree.range_query(&[0.0, 0.0], 1.0);
        let indices: Vec<usize> = result.iter().map(|n| n.index).collect();
        assert!(indices.contains(&0));
        assert!(indices.contains(&1));
        assert!(!indices.contains(&3));
    }

    #[test]
    fn test_ball_tree_len() {
        let points: Vec<Vec<f64>> = (0..10).map(|i| vec![i as f64]).collect();
        let tree = MetricBallTree::build(points, 3);
        assert_eq!(tree.len(), 10);
        assert!(!tree.is_empty());
    }

    // --- FrechetDistance ---

    #[test]
    fn test_frechet_identical_curves() {
        let c: Vec<Vec<f64>> = vec![vec![0.0, 0.0], vec![1.0, 0.0], vec![2.0, 0.0]];
        let d = FrechetDistance::compute(&c, &c);
        assert!(d < 1e-12, "identical curves → Fréchet = 0, got {d}");
    }

    #[test]
    fn test_frechet_offset_curves() {
        // Curve q is curve p shifted up by 1
        let p: Vec<Vec<f64>> = vec![vec![0.0, 0.0], vec![1.0, 0.0]];
        let q: Vec<Vec<f64>> = vec![vec![0.0, 1.0], vec![1.0, 1.0]];
        let d = FrechetDistance::compute(&p, &q);
        assert!((d - 1.0).abs() < 1e-10, "expected 1.0, got {d}");
    }

    #[test]
    fn test_frechet_single_point() {
        let p = vec![vec![0.0, 0.0]];
        let q = vec![vec![3.0, 4.0]];
        let d = FrechetDistance::compute(&p, &q);
        assert!((d - 5.0).abs() < 1e-10, "expected 5.0, got {d}");
    }

    #[test]
    fn test_frechet_with_manhattan_metric() {
        let p: Vec<Vec<f64>> = vec![vec![0.0, 0.0], vec![1.0, 0.0]];
        let q: Vec<Vec<f64>> = vec![vec![0.0, 1.0], vec![1.0, 1.0]];
        let d = FrechetDistance::compute_with_metric(&p, &q, &ManhattanMetric);
        assert!((d - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_hausdorff() {
        let p: Vec<Vec<f64>> = vec![vec![0.0, 0.0], vec![2.0, 0.0]];
        let q: Vec<Vec<f64>> = vec![vec![0.0, 1.0]];
        let h = FrechetDistance::hausdorff(&p, &q);
        // d([2,0],[0,1]) = sqrt(5) ≈ 2.236
        assert!(h > 2.0);
    }

    // --- Distance matrix and utilities ---

    #[test]
    fn test_distance_matrix_diagonal_zero() {
        let m = EuclideanMetric;
        let pts = vec![vec![0.0, 0.0], vec![1.0, 0.0], vec![0.0, 1.0]];
        let mat = distance_matrix(&m, &pts);
        for i in 0..3 {
            assert!(mat[i][i].abs() < 1e-12);
        }
    }

    #[test]
    fn test_distance_matrix_symmetry() {
        let m = ManhattanMetric;
        let pts = vec![vec![0.0, 0.0], vec![1.0, 2.0], vec![3.0, -1.0]];
        let mat = distance_matrix(&m, &pts);
        for i in 0..3 {
            for j in 0..3 {
                assert!((mat[i][j] - mat[j][i]).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn test_triangle_inequality_check_pass() {
        let m = EuclideanMetric;
        let pts = vec![vec![0.0, 0.0], vec![1.0, 0.0], vec![0.0, 1.0]];
        let mat = distance_matrix(&m, &pts);
        assert!(check_triangle_inequality(&mat).is_ok());
    }

    #[test]
    fn test_set_diameter_empty() {
        let m = EuclideanMetric;
        let pts: Vec<Vec<f64>> = vec![];
        assert!((set_diameter(&m, &pts)).abs() < 1e-12);
    }

    // --- GeodesicMetric ---

    #[test]
    fn test_geodesic_simple() {
        let mut g = GeodesicMetric::new(4);
        g.add_edge(0, 1, 1.0);
        g.add_edge(1, 2, 2.0);
        g.add_edge(2, 3, 1.0);
        assert!((g.distance(&0, &3) - 4.0).abs() < 1e-12);
    }

    #[test]
    fn test_geodesic_shortest_path() {
        let mut g = GeodesicMetric::new(3);
        g.add_edge(0, 1, 10.0);
        g.add_edge(0, 2, 1.0);
        g.add_edge(2, 1, 2.0);
        assert!((g.distance(&0, &1) - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_geodesic_self_distance() {
        let g = GeodesicMetric::new(3);
        assert!((g.distance(&0, &0)).abs() < 1e-12);
    }
}
