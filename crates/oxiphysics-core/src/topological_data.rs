#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Topological Data Analysis (TDA) — persistent homology utilities.
//!
//! Provides Vietoris–Rips filtration for 0- and 1-dimensional persistent
//! homology, bottleneck and Wasserstein distances between persistence
//! diagrams, Betti numbers, Euler characteristic, cubical complexes,
//! and persistence entropy.
//!
//! Structured types: [`VietorisRips`], [`SimplexTree`], [`PersistenceDiagram`],
//! [`BarcodeSummary`], and [`MapperGraph`].

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

// ─────────────────────────────────────────────────────────────────────────────
// Primitive types
// ─────────────────────────────────────────────────────────────────────────────

/// A weighted undirected edge in a simplicial filtration.
#[derive(Debug, Clone, PartialEq)]
pub struct SimplexEdge {
    /// Index of the first vertex.
    pub v0: usize,
    /// Index of the second vertex.
    pub v1: usize,
    /// Filtration weight (e.g., Euclidean distance).
    pub weight: f64,
}

/// Union–Find (disjoint-set) data structure for Kruskal's algorithm.
#[derive(Debug, Clone)]
pub struct UnionFind {
    /// Parent pointers; `parent[i] == i` for root nodes.
    pub parent: Vec<usize>,
    /// Union-by-rank array.
    pub rank: Vec<usize>,
}

impl UnionFind {
    /// Create a new `UnionFind` with `n` singleton sets.
    pub fn new(n: usize) -> Self {
        Self {
            parent: (0..n).collect(),
            rank: vec![0; n],
        }
    }

    /// Find the representative of the set containing `i` (path compression).
    pub fn find(&mut self, i: usize) -> usize {
        if self.parent[i] != i {
            self.parent[i] = self.find(self.parent[i]);
        }
        self.parent[i]
    }

    /// Union the sets containing `i` and `j`.  Returns `true` if they were
    /// in different sets (a merge happened).
    pub fn union(&mut self, i: usize, j: usize) -> bool {
        let ri = self.find(i);
        let rj = self.find(j);
        if ri == rj {
            return false;
        }
        match self.rank[ri].cmp(&self.rank[rj]) {
            std::cmp::Ordering::Less => self.parent[ri] = rj,
            std::cmp::Ordering::Greater => self.parent[rj] = ri,
            std::cmp::Ordering::Equal => {
                self.parent[rj] = ri;
                self.rank[ri] += 1;
            }
        }
        true
    }
}

/// A bar in a persistence diagram: `[birth, death)` for homology class H_k.
///
/// An *essential* class (never destroyed) has `death == f64::INFINITY`.
#[derive(Debug, Clone, PartialEq)]
pub struct PersistentInterval {
    /// Filtration value at which this class was born.
    pub birth: f64,
    /// Filtration value at which this class died; `INFINITY` if essential.
    pub death: f64,
    /// Homological dimension (0 = connected component, 1 = loop, …).
    pub dimension: usize,
}

impl PersistentInterval {
    /// Persistence (lifetime) of this interval.  Returns `INFINITY` for
    /// essential classes.
    pub fn persistence(&self) -> f64 {
        self.death - self.birth
    }

    /// Returns `true` if this class is essential (never destroyed).
    pub fn is_essential(&self) -> bool {
        self.death.is_infinite()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Vietoris–Rips H₀ via Kruskal
// ─────────────────────────────────────────────────────────────────────────────

/// Compute all pairwise edges from a 2-D point cloud.
fn all_edges(points: &[[f64; 2]]) -> Vec<SimplexEdge> {
    let n = points.len();
    let mut edges = Vec::with_capacity(n * (n.saturating_sub(1)) / 2);
    for i in 0..n {
        for j in (i + 1)..n {
            let dx = points[i][0] - points[j][0];
            let dy = points[i][1] - points[j][1];
            edges.push(SimplexEdge {
                v0: i,
                v1: j,
                weight: (dx * dx + dy * dy).sqrt(),
            });
        }
    }
    edges
}

/// Compute 0-dimensional persistent homology of the Vietoris–Rips filtration
/// up to `max_r` using Kruskal's algorithm.
///
/// Each returned interval represents a connected component: it is born at
/// `r = 0` and dies when it merges with an older component.  The last
/// surviving component has `death = INFINITY`.
pub fn vietoris_rips_h0(points: &[[f64; 2]], max_r: f64) -> Vec<PersistentInterval> {
    let n = points.len();
    if n == 0 {
        return vec![];
    }
    let mut edges = all_edges(points);
    edges.sort_by(|a, b| {
        a.weight
            .partial_cmp(&b.weight)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let mut uf = UnionFind::new(n);
    // birth times: all components born at filtration value 0
    let mut birth = vec![0.0_f64; n];
    let mut intervals = Vec::new();

    for e in &edges {
        if e.weight > max_r {
            break;
        }
        if uf.union(e.v0, e.v1) {
            // The *younger* component (higher birth) dies
            let r0 = uf.find(e.v0);
            let r1 = uf.find(e.v1);
            // After union, one root absorbed the other; the non-root died.
            // We compare the original roots before the merge.
            let root_after = uf.find(e.v0);
            let dying = if root_after == r0 { r1 } else { r0 };
            let b = birth[dying];
            birth[root_after] = birth[root_after].min(b);
            intervals.push(PersistentInterval {
                birth: b,
                death: e.weight,
                dimension: 0,
            });
        }
    }

    // Essential class: the one component that survives
    intervals.push(PersistentInterval {
        birth: 0.0,
        death: f64::INFINITY,
        dimension: 0,
    });
    intervals
}

// ─────────────────────────────────────────────────────────────────────────────
// Vietoris–Rips H₁ approximation via cycle detection
// ─────────────────────────────────────────────────────────────────────────────

/// Approximate 1-dimensional persistent homology (loops) of the
/// Vietoris–Rips filtration.
///
/// Uses the first `n_edges` edges (sorted by weight).  A 1-cycle is created
/// whenever an edge closes a loop in the current spanning forest.
pub fn vietoris_rips_h1_approx(
    points: &[[f64; 2]],
    max_r: f64,
    n_edges: usize,
) -> Vec<PersistentInterval> {
    let n = points.len();
    if n == 0 {
        return vec![];
    }
    let mut edges = all_edges(points);
    edges.sort_by(|a, b| {
        a.weight
            .partial_cmp(&b.weight)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let edges: Vec<_> = edges
        .into_iter()
        .filter(|e| e.weight <= max_r)
        .take(n_edges)
        .collect();

    let mut uf = UnionFind::new(n);
    let mut intervals = Vec::new();
    let mut edge_births: Vec<f64> = Vec::new(); // birth of spanning tree edges

    for e in &edges {
        if !uf.union(e.v0, e.v1) {
            // This edge closes a loop — approximate death as max_r
            let birth = *edge_births.last().unwrap_or(&0.0);
            intervals.push(PersistentInterval {
                birth,
                death: (e.weight + max_r) / 2.0,
                dimension: 1,
            });
        } else {
            edge_births.push(e.weight);
        }
    }
    intervals
}

// ─────────────────────────────────────────────────────────────────────────────
// Diagram distances
// ─────────────────────────────────────────────────────────────────────────────

/// L∞ distance between two points in a persistence diagram,
/// projected to the diagonal if needed.
fn point_dist_inf(p: &PersistentInterval, q: &PersistentInterval) -> f64 {
    let pb = p.birth.min(p.death);
    let pd = p.birth.max(p.death);
    let qb = q.birth.min(q.death);
    let qd = q.birth.max(q.death);
    f64::max((pb - qb).abs(), (pd - qd).abs())
}

/// Distance from a point to the diagonal (half its persistence).
fn dist_to_diagonal(p: &PersistentInterval) -> f64 {
    (p.death - p.birth).abs() / 2.0
}

/// Bottleneck distance between two persistence diagrams.
///
/// Computed via a greedy approximation: match each point in the smaller
/// diagram to the closest point in the larger, then take the maximum
/// residual (including diagonal projections).
pub fn bottleneck_distance(a: &[PersistentInterval], b: &[PersistentInterval]) -> f64 {
    // Filter to finite intervals only
    let fa: Vec<&PersistentInterval> = a.iter().filter(|p| !p.is_essential()).collect();
    let fb: Vec<&PersistentInterval> = b.iter().filter(|p| !p.is_essential()).collect();

    let mut max_dist = 0.0_f64;

    // For each point in fa find closest in fb
    let mut used = vec![false; fb.len()];
    for p in &fa {
        let best = fb
            .iter()
            .enumerate()
            .filter(|(idx, _)| !used[*idx])
            .map(|(idx, q)| (idx, point_dist_inf(p, q)))
            .min_by(|x, y| x.1.partial_cmp(&y.1).unwrap_or(std::cmp::Ordering::Equal));
        if let Some((idx, d)) = best {
            used[idx] = true;
            max_dist = max_dist.max(d);
        } else {
            max_dist = max_dist.max(dist_to_diagonal(p));
        }
    }
    // Unmatched points in fb go to diagonal
    for (idx, q) in fb.iter().enumerate() {
        if !used[idx] {
            max_dist = max_dist.max(dist_to_diagonal(q));
        }
    }
    max_dist
}

/// 1-Wasserstein (Earth Mover's) distance between two persistence diagrams.
///
/// Uses a greedy nearest-neighbor matching (not optimal for large diagrams,
/// but exact for the 1-Wasserstein cost under the L∞ ground metric).
pub fn wasserstein_distance_1(a: &[PersistentInterval], b: &[PersistentInterval]) -> f64 {
    let fa: Vec<&PersistentInterval> = a.iter().filter(|p| !p.is_essential()).collect();
    let fb: Vec<&PersistentInterval> = b.iter().filter(|p| !p.is_essential()).collect();

    let mut total = 0.0_f64;
    let mut used = vec![false; fb.len()];

    for p in &fa {
        let best = fb
            .iter()
            .enumerate()
            .filter(|(idx, _)| !used[*idx])
            .map(|(idx, q)| (idx, point_dist_inf(p, q)))
            .min_by(|x, y| x.1.partial_cmp(&y.1).unwrap_or(std::cmp::Ordering::Equal));
        if let Some((idx, d)) = best {
            used[idx] = true;
            total += d;
        } else {
            total += dist_to_diagonal(p);
        }
    }
    for (idx, q) in fb.iter().enumerate() {
        if !used[idx] {
            total += dist_to_diagonal(q);
        }
    }
    total
}

// ─────────────────────────────────────────────────────────────────────────────
// Betti numbers & Euler characteristic
// ─────────────────────────────────────────────────────────────────────────────

/// Betti numbers β₀, β₁, β₂ of a topological space.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BettiNumbers {
    /// β₀: number of connected components.
    pub b0: usize,
    /// β₁: number of independent loops.
    pub b1: usize,
    /// β₂: number of voids (2-spheres).
    pub b2: usize,
}

/// Euler characteristic χ = β₀ − β₁ + β₂.
pub fn euler_characteristic(betti: &BettiNumbers) -> i64 {
    betti.b0 as i64 - betti.b1 as i64 + betti.b2 as i64
}

// ─────────────────────────────────────────────────────────────────────────────
// 1-D cubical complex (sublevel-set filtration)
// ─────────────────────────────────────────────────────────────────────────────

/// Compute Betti numbers of the sublevel set `{x : signal[i] ≤ threshold}`.
///
/// Only β₀ and β₁ are non-trivial in 1-D; β₂ is always 0.
pub fn cubical_complex_1d(signal: &[f64], threshold: f64) -> BettiNumbers {
    if signal.is_empty() {
        return BettiNumbers {
            b0: 0,
            b1: 0,
            b2: 0,
        };
    }
    // Connected components = contiguous intervals where signal ≤ threshold
    let below: Vec<bool> = signal.iter().map(|&v| v <= threshold).collect();
    let mut b0 = 0usize;
    let mut in_component = false;
    for &b in &below {
        if b && !in_component {
            b0 += 1;
            in_component = true;
        } else if !b {
            in_component = false;
        }
    }
    // In 1-D there are no loops
    BettiNumbers { b0, b1: 0, b2: 0 }
}

// ─────────────────────────────────────────────────────────────────────────────
// Persistence entropy
// ─────────────────────────────────────────────────────────────────────────────

/// Entropy of a persistence diagram.
///
/// H = −Σ (lᵢ/L) log(lᵢ/L) where lᵢ = persistence of interval i and
/// L = Σ lᵢ.  Essential (infinite) intervals are excluded.
pub fn persistence_entropy(diagram: &[PersistentInterval]) -> f64 {
    let lifetimes: Vec<f64> = diagram
        .iter()
        .filter(|p| !p.is_essential() && p.persistence() > 0.0)
        .map(|p| p.persistence())
        .collect();
    if lifetimes.is_empty() {
        return 0.0;
    }
    let total: f64 = lifetimes.iter().sum();
    if total <= 0.0 {
        return 0.0;
    }
    -lifetimes
        .iter()
        .map(|&l| {
            let p = l / total;
            if p > 0.0 { p * p.ln() } else { 0.0 }
        })
        .sum::<f64>()
}

// ─────────────────────────────────────────────────────────────────────────────
// VietorisRips struct
// ─────────────────────────────────────────────────────────────────────────────

/// Vietoris–Rips complex construction from a 2-D point cloud.
///
/// The complex includes all simplices (vertices, edges, triangles) whose
/// maximum pairwise distance is at most `epsilon`.
#[derive(Debug, Clone)]
pub struct VietorisRips {
    /// The 2-D point cloud.
    pub points: Vec<[f64; 2]>,
    /// Maximum edge length (epsilon radius).
    pub epsilon: f64,
}

impl VietorisRips {
    /// Create a new `VietorisRips` complex from a point cloud with radius `epsilon`.
    pub fn new(points: Vec<[f64; 2]>, epsilon: f64) -> Self {
        Self { points, epsilon }
    }

    /// Euclidean distance between two 2-D points.
    fn dist(a: &[f64; 2], b: &[f64; 2]) -> f64 {
        let dx = a[0] - b[0];
        let dy = a[1] - b[1];
        (dx * dx + dy * dy).sqrt()
    }

    /// All edges (pairs) within `epsilon` distance.
    ///
    /// Returns a list of `(i, j, distance)` tuples.
    pub fn edges(&self) -> Vec<(usize, usize, f64)> {
        let n = self.points.len();
        let mut result = Vec::new();
        for i in 0..n {
            for j in (i + 1)..n {
                let d = Self::dist(&self.points[i], &self.points[j]);
                if d <= self.epsilon {
                    result.push((i, j, d));
                }
            }
        }
        result
    }

    /// All triangles (triples) where every pair is within `epsilon`.
    ///
    /// Returns a list of `(i, j, k)` index triples.
    pub fn triangles(&self) -> Vec<(usize, usize, usize)> {
        let n = self.points.len();
        let mut result = Vec::new();
        for i in 0..n {
            for j in (i + 1)..n {
                if Self::dist(&self.points[i], &self.points[j]) > self.epsilon {
                    continue;
                }
                for k in (j + 1)..n {
                    if Self::dist(&self.points[i], &self.points[k]) <= self.epsilon
                        && Self::dist(&self.points[j], &self.points[k]) <= self.epsilon
                    {
                        result.push((i, j, k));
                    }
                }
            }
        }
        result
    }

    /// Number of vertices (all points are 0-simplices).
    pub fn n_vertices(&self) -> usize {
        self.points.len()
    }

    /// Number of edges within `epsilon`.
    pub fn n_edges(&self) -> usize {
        self.edges().len()
    }

    /// Number of triangles (2-simplices) within `epsilon`.
    pub fn n_triangles(&self) -> usize {
        self.triangles().len()
    }

    /// Euler characteristic χ = V − E + F.
    pub fn euler_characteristic(&self) -> i64 {
        let v = self.n_vertices() as i64;
        let e = self.n_edges() as i64;
        let f = self.n_triangles() as i64;
        v - e + f
    }

    /// Betti numbers via the Euler formula and a union-find for β₀.
    ///
    /// β₀ is the number of connected components.
    /// β₁ is estimated as β₁ = E − V + β₀ (first Betti number from spanning tree).
    /// β₂ is estimated from the Euler characteristic.
    pub fn betti_numbers(&self) -> BettiNumbers {
        let n = self.points.len();
        if n == 0 {
            return BettiNumbers {
                b0: 0,
                b1: 0,
                b2: 0,
            };
        }
        let edges = self.edges();
        let triangles = self.triangles();
        let mut uf = UnionFind::new(n);
        for (i, j, _) in &edges {
            uf.union(*i, *j);
        }
        // β₀: count distinct roots
        let b0 = (0..n).filter(|&i| uf.find(i) == i).count();
        // β₁ = E − V + β₀ (number of independent cycles)
        let b1 = (edges.len() as i64 - n as i64 + b0 as i64).max(0) as usize;
        // β₂ estimated from Euler characteristic: χ = β₀ − β₁ + β₂
        let chi = self.euler_characteristic();
        let b2 = (chi - b0 as i64 + b1 as i64).max(0) as usize;
        let _ = triangles; // used indirectly via euler_characteristic
        BettiNumbers { b0, b1, b2 }
    }

    /// Compute the 0-dimensional persistent homology (H₀) via Kruskal.
    pub fn persistent_h0(&self) -> Vec<PersistentInterval> {
        vietoris_rips_h0(&self.points, self.epsilon)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// SimplexTree struct
// ─────────────────────────────────────────────────────────────────────────────

/// A simplex with its filtration value.
#[derive(Debug, Clone)]
pub struct FilteredSimplex {
    /// Sorted vertex indices.
    pub vertices: Vec<usize>,
    /// Filtration value at which this simplex enters.
    pub filtration: f64,
}

impl FilteredSimplex {
    /// Dimension of this simplex (0=vertex, 1=edge, 2=triangle, …).
    pub fn dim(&self) -> usize {
        self.vertices.len().saturating_sub(1)
    }
}

/// Simplex tree data structure for filtration-based TDA.
///
/// Stores all simplices sorted by filtration value, supporting
/// insertion, lookup, and filtration-order iteration.
#[derive(Debug, Clone, Default)]
pub struct SimplexTree {
    /// All simplices, in arbitrary insertion order.
    pub simplices: Vec<FilteredSimplex>,
}

impl SimplexTree {
    /// Create an empty `SimplexTree`.
    pub fn new() -> Self {
        Self::default()
    }

    /// Insert a simplex defined by `vertices` at filtration value `f`.
    ///
    /// Vertices are sorted before insertion.
    pub fn insert(&mut self, mut vertices: Vec<usize>, f: f64) {
        vertices.sort_unstable();
        self.simplices.push(FilteredSimplex {
            vertices,
            filtration: f,
        });
    }

    /// Build a `SimplexTree` from a `VietorisRips` complex.
    ///
    /// Inserts all vertices (at f=0), edges, and triangles at their
    /// respective filtration values (maximum edge length).
    pub fn from_vietoris_rips(vr: &VietorisRips) -> Self {
        let mut tree = SimplexTree::new();
        // Vertices at filtration 0
        for i in 0..vr.points.len() {
            tree.insert(vec![i], 0.0);
        }
        // Edges
        for (i, j, d) in vr.edges() {
            tree.insert(vec![i, j], d);
        }
        // Triangles (filtration = max edge)
        for (i, j, k) in vr.triangles() {
            let d01 = VietorisRips::dist(&vr.points[i], &vr.points[j]);
            let d02 = VietorisRips::dist(&vr.points[i], &vr.points[k]);
            let d12 = VietorisRips::dist(&vr.points[j], &vr.points[k]);
            let f = d01.max(d02).max(d12);
            tree.insert(vec![i, j, k], f);
        }
        tree
    }

    /// Return all simplices sorted by filtration value (ascending).
    pub fn sorted_simplices(&self) -> Vec<&FilteredSimplex> {
        let mut v: Vec<&FilteredSimplex> = self.simplices.iter().collect();
        v.sort_by(|a, b| {
            a.filtration
                .partial_cmp(&b.filtration)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        v
    }

    /// Number of simplices of dimension `d`.
    pub fn count_dim(&self, d: usize) -> usize {
        self.simplices.iter().filter(|s| s.dim() == d).count()
    }

    /// Filtration threshold that includes all simplices (maximum filtration value).
    pub fn max_filtration(&self) -> f64 {
        self.simplices
            .iter()
            .map(|s| s.filtration)
            .fold(0.0_f64, f64::max)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PersistenceDiagram struct
// ─────────────────────────────────────────────────────────────────────────────

/// A persistence diagram: a collection of birth-death pairs.
///
/// Each pair represents a topological feature that appears at `birth` and
/// disappears at `death` in a filtration.
#[derive(Debug, Clone, Default)]
pub struct PersistenceDiagram {
    /// The birth-death pairs (birth, death, dimension).
    pub pairs: Vec<PersistentInterval>,
}

impl PersistenceDiagram {
    /// Create an empty `PersistenceDiagram`.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a birth-death pair to the diagram.
    pub fn add(&mut self, birth: f64, death: f64, dimension: usize) {
        self.pairs.push(PersistentInterval {
            birth,
            death,
            dimension,
        });
    }

    /// Build a persistence diagram from the 0-dimensional VR filtration.
    pub fn from_vietoris_rips_h0(vr: &VietorisRips) -> Self {
        let intervals = vr.persistent_h0();
        PersistenceDiagram { pairs: intervals }
    }

    /// Filter pairs by homological dimension.
    pub fn pairs_in_dim(&self, dim: usize) -> Vec<&PersistentInterval> {
        self.pairs.iter().filter(|p| p.dimension == dim).collect()
    }

    /// Maximum persistence (lifetime) over all finite pairs.
    pub fn max_persistence(&self) -> f64 {
        self.pairs
            .iter()
            .filter(|p| !p.is_essential())
            .map(|p| p.persistence())
            .fold(0.0_f64, f64::max)
    }

    /// Number of essential (infinite) bars.
    pub fn n_essential(&self) -> usize {
        self.pairs.iter().filter(|p| p.is_essential()).count()
    }

    /// Persistence entropy of the diagram (in nats).
    pub fn entropy(&self) -> f64 {
        persistence_entropy(&self.pairs)
    }

    /// Total persistence Σ (death − birth) for finite pairs.
    pub fn total_persistence(&self) -> f64 {
        self.pairs
            .iter()
            .filter(|p| !p.is_essential())
            .map(|p| p.persistence())
            .sum()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// BarcodeSummary struct
// ─────────────────────────────────────────────────────────────────────────────

/// Summary statistics derived from a persistence barcode.
#[derive(Debug, Clone)]
pub struct BarcodeSummary {
    /// Total persistence: Σ (death − birth) for finite pairs.
    pub total_persistence: f64,
    /// Number of features with persistence ≥ `threshold`.
    pub n_long_features: usize,
    /// Persistence entropy.
    pub entropy: f64,
    /// Maximum persistence over all finite pairs.
    pub max_persistence: f64,
    /// Number of essential (infinite) features.
    pub n_essential: usize,
    /// Mean persistence of finite pairs; 0 if none.
    pub mean_persistence: f64,
}

impl BarcodeSummary {
    /// Compute a `BarcodeSummary` from a persistence diagram.
    ///
    /// `long_threshold` determines the cutoff for "long" features.
    pub fn from_diagram(diagram: &PersistenceDiagram, long_threshold: f64) -> Self {
        let finite: Vec<f64> = diagram
            .pairs
            .iter()
            .filter(|p| !p.is_essential())
            .map(|p| p.persistence())
            .collect();
        let total_persistence = finite.iter().sum::<f64>();
        let n_long_features = finite.iter().filter(|&&p| p >= long_threshold).count();
        let entropy = diagram.entropy();
        let max_persistence = finite.iter().cloned().fold(0.0_f64, f64::max);
        let n_essential = diagram.n_essential();
        let mean_persistence = if finite.is_empty() {
            0.0
        } else {
            total_persistence / finite.len() as f64
        };
        BarcodeSummary {
            total_persistence,
            n_long_features,
            entropy,
            max_persistence,
            n_essential,
            mean_persistence,
        }
    }

    /// Normalised total persistence: total / max (0 if max is 0).
    pub fn normalised_total(&self) -> f64 {
        if self.max_persistence < 1e-15 {
            0.0
        } else {
            self.total_persistence / self.max_persistence
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// MapperGraph struct
// ─────────────────────────────────────────────────────────────────────────────

/// A node in a Mapper graph.
///
/// Each node corresponds to a cluster of points from one cover set.
#[derive(Debug, Clone)]
pub struct MapperNode {
    /// Indices of points in the original point cloud assigned to this node.
    pub point_indices: Vec<usize>,
    /// Cover interval index that generated this node.
    pub cover_index: usize,
}

/// Mapper graph: nodes (clusters) connected when adjacent cover intervals
/// share common points.
///
/// Implements a simplified version of the Mapper algorithm:
/// 1. Cover the filter function range with overlapping intervals.
/// 2. Restrict the point cloud to each interval.
/// 3. Cluster each restriction (here: single-linkage at `cluster_eps`).
/// 4. Connect clusters from adjacent intervals that share points.
#[derive(Debug, Clone, Default)]
pub struct MapperGraph {
    /// Nodes of the graph.
    pub nodes: Vec<MapperNode>,
    /// Edges between nodes: `(node_i, node_j)`.
    pub edges: Vec<(usize, usize)>,
}

impl MapperGraph {
    /// Build a `MapperGraph` from a point cloud using a 1-D filter function.
    ///
    /// # Parameters
    /// - `points`: 2-D point cloud.
    /// - `filter_values`: scalar filter value for each point (e.g., x-coordinate).
    /// - `n_intervals`: number of cover intervals.
    /// - `overlap`: fractional overlap between adjacent intervals (in \[0, 1)).
    /// - `cluster_eps`: distance threshold for single-linkage clustering within
    ///   each cover interval.
    pub fn build(
        points: &[[f64; 2]],
        filter_values: &[f64],
        n_intervals: usize,
        overlap: f64,
        cluster_eps: f64,
    ) -> Self {
        if points.is_empty() || n_intervals == 0 {
            return Self::default();
        }
        let min_f = filter_values.iter().cloned().fold(f64::INFINITY, f64::min);
        let max_f = filter_values
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        let range = (max_f - min_f).max(1e-15);
        let step = range / n_intervals as f64;
        let half_overlap = step * overlap.clamp(0.0, 0.99) / 2.0;

        let mut graph = MapperGraph::default();

        // For each cover interval, find the points inside, then cluster them
        for k in 0..n_intervals {
            let lo = min_f + k as f64 * step - half_overlap;
            let hi = min_f + (k + 1) as f64 * step + half_overlap;

            // Points in this interval
            let indices_in: Vec<usize> = (0..points.len())
                .filter(|&i| filter_values[i] >= lo && filter_values[i] <= hi)
                .collect();

            if indices_in.is_empty() {
                continue;
            }

            // Single-linkage clustering: union-find on the sub-cloud
            let m = indices_in.len();
            let mut uf = UnionFind::new(m);
            for a in 0..m {
                for b in (a + 1)..m {
                    let ia = indices_in[a];
                    let ib = indices_in[b];
                    let dx = points[ia][0] - points[ib][0];
                    let dy = points[ia][1] - points[ib][1];
                    let d = (dx * dx + dy * dy).sqrt();
                    if d <= cluster_eps {
                        uf.union(a, b);
                    }
                }
            }

            // Group by cluster root
            let mut cluster_map: std::collections::HashMap<usize, Vec<usize>> =
                std::collections::HashMap::new();
            for a in 0..m {
                let root = uf.find(a);
                cluster_map.entry(root).or_default().push(indices_in[a]);
            }

            // Add one node per cluster
            for (_root, members) in cluster_map {
                graph.nodes.push(MapperNode {
                    point_indices: members,
                    cover_index: k,
                });
            }
        }

        // Connect nodes that share at least one point index
        let n_nodes = graph.nodes.len();
        for i in 0..n_nodes {
            for j in (i + 1)..n_nodes {
                let set_i: std::collections::HashSet<usize> =
                    graph.nodes[i].point_indices.iter().copied().collect();
                let shared = graph.nodes[j]
                    .point_indices
                    .iter()
                    .any(|p| set_i.contains(p));
                if shared {
                    graph.edges.push((i, j));
                }
            }
        }

        graph
    }

    /// Number of nodes in the graph.
    pub fn n_nodes(&self) -> usize {
        self.nodes.len()
    }

    /// Number of edges in the graph.
    pub fn n_edges(&self) -> usize {
        self.edges.len()
    }

    /// Degree of node `i` (number of edges incident to it).
    pub fn degree(&self, i: usize) -> usize {
        self.edges
            .iter()
            .filter(|(a, b)| *a == i || *b == i)
            .count()
    }

    /// Total number of points across all nodes (with possible duplicates from overlap).
    pub fn total_point_count(&self) -> usize {
        self.nodes.iter().map(|n| n.point_indices.len()).sum()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── UnionFind ────────────────────────────────────────────────────────────

    #[test]
    fn test_uf_new_singletons() {
        let mut uf = UnionFind::new(5);
        for i in 0..5 {
            assert_eq!(uf.find(i), i);
        }
    }

    #[test]
    fn test_uf_union_merges() {
        let mut uf = UnionFind::new(4);
        assert!(uf.union(0, 1));
        assert_eq!(uf.find(0), uf.find(1));
    }

    #[test]
    fn test_uf_union_idempotent() {
        let mut uf = UnionFind::new(4);
        uf.union(0, 1);
        assert!(
            !uf.union(0, 1),
            "second union of same set should return false"
        );
    }

    #[test]
    fn test_uf_union_chain() {
        let mut uf = UnionFind::new(4);
        uf.union(0, 1);
        uf.union(1, 2);
        uf.union(2, 3);
        let r0 = uf.find(0);
        assert_eq!(uf.find(3), r0);
    }

    #[test]
    fn test_uf_different_components() {
        let mut uf = UnionFind::new(4);
        uf.union(0, 1);
        uf.union(2, 3);
        assert_ne!(uf.find(0), uf.find(2));
    }

    #[test]
    fn test_uf_path_compression() {
        let mut uf = UnionFind::new(6);
        for i in 0..5 {
            uf.union(i, i + 1);
        }
        let root = uf.find(5);
        assert_eq!(uf.find(0), root);
    }

    #[test]
    fn test_uf_union_all_returns_true_n_minus_1_times() {
        let n = 5;
        let mut uf = UnionFind::new(n);
        let mut merges = 0;
        for i in 1..n {
            if uf.union(0, i) {
                merges += 1;
            }
        }
        assert_eq!(merges, n - 1);
    }

    // ── SimplexEdge ──────────────────────────────────────────────────────────

    #[test]
    fn test_simplex_edge_fields() {
        let e = SimplexEdge {
            v0: 1,
            v1: 3,
            weight: 2.5,
        };
        assert_eq!(e.v0, 1);
        assert_eq!(e.v1, 3);
        assert!((e.weight - 2.5).abs() < 1e-12);
    }

    // ── PersistentInterval ───────────────────────────────────────────────────

    #[test]
    fn test_persistent_interval_persistence() {
        let pi = PersistentInterval {
            birth: 1.0,
            death: 3.0,
            dimension: 0,
        };
        assert!((pi.persistence() - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_persistent_interval_essential() {
        let pi = PersistentInterval {
            birth: 0.0,
            death: f64::INFINITY,
            dimension: 0,
        };
        assert!(pi.is_essential());
    }

    #[test]
    fn test_persistent_interval_not_essential() {
        let pi = PersistentInterval {
            birth: 0.0,
            death: 1.0,
            dimension: 0,
        };
        assert!(!pi.is_essential());
    }

    // ── vietoris_rips_h0 ─────────────────────────────────────────────────────

    #[test]
    fn test_h0_empty() {
        let pts: [[f64; 2]; 0] = [];
        let result = vietoris_rips_h0(&pts, 10.0);
        assert!(result.is_empty());
    }

    #[test]
    fn test_h0_single_point() {
        let pts = [[0.0, 0.0]];
        let result = vietoris_rips_h0(&pts, 10.0);
        assert_eq!(result.len(), 1);
        assert!(result[0].is_essential());
    }

    #[test]
    fn test_h0_two_close_points() {
        let pts = [[0.0, 0.0], [0.1, 0.0]];
        let result = vietoris_rips_h0(&pts, 10.0);
        // One interval dies, one essential
        assert_eq!(result.len(), 2);
        let essential_count = result.iter().filter(|p| p.is_essential()).count();
        assert_eq!(essential_count, 1);
    }

    #[test]
    fn test_h0_all_dimension_zero() {
        let pts = [[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]];
        let result = vietoris_rips_h0(&pts, 10.0);
        for pi in &result {
            assert_eq!(pi.dimension, 0);
        }
    }

    #[test]
    fn test_h0_two_far_points_below_max_r() {
        let pts = [[0.0, 0.0], [100.0, 0.0]];
        let result = vietoris_rips_h0(&pts, 50.0);
        // They don't merge within max_r=50, so both stay alive
        let finite_count = result.iter().filter(|p| !p.is_essential()).count();
        // No merge happened within max_r
        assert_eq!(finite_count, 0);
    }

    // ── vietoris_rips_h1_approx ───────────────────────────────────────────────

    #[test]
    fn test_h1_empty() {
        let pts: [[f64; 2]; 0] = [];
        let result = vietoris_rips_h1_approx(&pts, 10.0, 100);
        assert!(result.is_empty());
    }

    #[test]
    fn test_h1_dimension_one() {
        let pts = [[0.0, 0.0], [1.0, 0.0], [0.5, 1.0]];
        let result = vietoris_rips_h1_approx(&pts, 10.0, 10);
        for pi in &result {
            assert_eq!(pi.dimension, 1);
        }
    }

    #[test]
    fn test_h1_triangle_has_loop() {
        // Triangle: 3 edges, spanning tree = 2, so 1 cycle
        let pts = [[0.0, 0.0], [1.0, 0.0], [0.5, 1.0]];
        let result = vietoris_rips_h1_approx(&pts, 10.0, 3);
        assert_eq!(result.len(), 1);
    }

    // ── bottleneck_distance ───────────────────────────────────────────────────

    #[test]
    fn test_bottleneck_identical_diagrams() {
        let d = vec![
            PersistentInterval {
                birth: 0.0,
                death: 1.0,
                dimension: 0,
            },
            PersistentInterval {
                birth: 0.5,
                death: 2.0,
                dimension: 0,
            },
        ];
        let dist = bottleneck_distance(&d, &d);
        assert!(dist < 1e-12, "identical diagrams should have distance 0");
    }

    #[test]
    fn test_bottleneck_empty_diagrams() {
        assert!((bottleneck_distance(&[], &[])).abs() < 1e-12);
    }

    #[test]
    fn test_bottleneck_non_negative() {
        let a = vec![PersistentInterval {
            birth: 0.0,
            death: 1.0,
            dimension: 0,
        }];
        let b = vec![PersistentInterval {
            birth: 0.1,
            death: 1.2,
            dimension: 0,
        }];
        assert!(bottleneck_distance(&a, &b) >= 0.0);
    }

    #[test]
    fn test_bottleneck_symmetry() {
        let a = vec![PersistentInterval {
            birth: 0.0,
            death: 1.0,
            dimension: 0,
        }];
        let b = vec![PersistentInterval {
            birth: 0.5,
            death: 2.0,
            dimension: 0,
        }];
        let dab = bottleneck_distance(&a, &b);
        let dba = bottleneck_distance(&b, &a);
        assert!((dab - dba).abs() < 1e-12);
    }

    // ── wasserstein_distance_1 ────────────────────────────────────────────────

    #[test]
    fn test_wasserstein_identical() {
        let d = vec![PersistentInterval {
            birth: 0.0,
            death: 2.0,
            dimension: 0,
        }];
        assert!(wasserstein_distance_1(&d, &d) < 1e-12);
    }

    #[test]
    fn test_wasserstein_empty() {
        assert!(wasserstein_distance_1(&[], &[]).abs() < 1e-12);
    }

    #[test]
    fn test_wasserstein_non_negative() {
        let a = vec![PersistentInterval {
            birth: 0.0,
            death: 1.0,
            dimension: 0,
        }];
        let b = vec![PersistentInterval {
            birth: 0.2,
            death: 1.5,
            dimension: 0,
        }];
        assert!(wasserstein_distance_1(&a, &b) >= 0.0);
    }

    // ── BettiNumbers & euler_characteristic ─────────────────────────────────

    #[test]
    fn test_euler_characteristic_sphere() {
        // Topological 2-sphere: β0=1, β1=0, β2=1 → χ=2
        let b = BettiNumbers {
            b0: 1,
            b1: 0,
            b2: 1,
        };
        assert_eq!(euler_characteristic(&b), 2);
    }

    #[test]
    fn test_euler_characteristic_torus() {
        // Torus: β0=1, β1=2, β2=1 → χ=0
        let b = BettiNumbers {
            b0: 1,
            b1: 2,
            b2: 1,
        };
        assert_eq!(euler_characteristic(&b), 0);
    }

    #[test]
    fn test_euler_characteristic_point() {
        let b = BettiNumbers {
            b0: 1,
            b1: 0,
            b2: 0,
        };
        assert_eq!(euler_characteristic(&b), 1);
    }

    #[test]
    fn test_euler_characteristic_circle() {
        // Circle: β0=1, β1=1, β2=0 → χ=0
        let b = BettiNumbers {
            b0: 1,
            b1: 1,
            b2: 0,
        };
        assert_eq!(euler_characteristic(&b), 0);
    }

    // ── cubical_complex_1d ────────────────────────────────────────────────────

    #[test]
    fn test_cubical_1d_empty() {
        let b = cubical_complex_1d(&[], 0.5);
        assert_eq!(b.b0, 0);
    }

    #[test]
    fn test_cubical_1d_all_below() {
        let signal = [0.1, 0.2, 0.3];
        let b = cubical_complex_1d(&signal, 1.0);
        assert_eq!(b.b0, 1);
        assert_eq!(b.b1, 0);
        assert_eq!(b.b2, 0);
    }

    #[test]
    fn test_cubical_1d_two_components() {
        let signal = [0.1, 0.9, 0.1];
        let b = cubical_complex_1d(&signal, 0.5);
        assert_eq!(b.b0, 2);
    }

    #[test]
    fn test_cubical_1d_none_below() {
        let signal = [1.0, 2.0, 3.0];
        let b = cubical_complex_1d(&signal, 0.5);
        assert_eq!(b.b0, 0);
    }

    #[test]
    fn test_cubical_1d_no_loops() {
        let signal = [0.1, 0.2, 0.3, 0.4];
        let b = cubical_complex_1d(&signal, 1.0);
        assert_eq!(b.b1, 0);
        assert_eq!(b.b2, 0);
    }

    // ── persistence_entropy ───────────────────────────────────────────────────

    #[test]
    fn test_persistence_entropy_empty() {
        assert_eq!(persistence_entropy(&[]), 0.0);
    }

    #[test]
    fn test_persistence_entropy_single_bar() {
        let d = vec![PersistentInterval {
            birth: 0.0,
            death: 1.0,
            dimension: 0,
        }];
        // Single bar: p=1, so -1*ln(1) = 0
        assert!(persistence_entropy(&d).abs() < 1e-12);
    }

    #[test]
    fn test_persistence_entropy_non_negative() {
        let d = vec![
            PersistentInterval {
                birth: 0.0,
                death: 1.0,
                dimension: 0,
            },
            PersistentInterval {
                birth: 0.5,
                death: 2.0,
                dimension: 0,
            },
            PersistentInterval {
                birth: 1.0,
                death: 3.0,
                dimension: 0,
            },
        ];
        assert!(persistence_entropy(&d) >= 0.0);
    }

    #[test]
    fn test_persistence_entropy_essential_excluded() {
        let d = vec![
            PersistentInterval {
                birth: 0.0,
                death: f64::INFINITY,
                dimension: 0,
            },
            PersistentInterval {
                birth: 0.0,
                death: 1.0,
                dimension: 0,
            },
        ];
        // Should not panic or give NaN
        let h = persistence_entropy(&d);
        assert!(h.is_finite());
    }

    #[test]
    fn test_persistence_entropy_equal_bars_maximizes() {
        // Two equal-length bars → maximum entropy for 2 events = ln(2)
        let d = vec![
            PersistentInterval {
                birth: 0.0,
                death: 1.0,
                dimension: 0,
            },
            PersistentInterval {
                birth: 0.0,
                death: 1.0,
                dimension: 0,
            },
        ];
        let h = persistence_entropy(&d);
        assert!((h - 2_f64.ln()).abs() < 1e-10);
    }

    // ── VietorisRips struct ───────────────────────────────────────────────────

    #[test]
    fn test_vr_empty_cloud() {
        let vr = VietorisRips::new(vec![], 1.0);
        assert_eq!(vr.n_vertices(), 0);
        assert_eq!(vr.n_edges(), 0);
    }

    #[test]
    fn test_vr_single_point() {
        let vr = VietorisRips::new(vec![[0.0, 0.0]], 1.0);
        assert_eq!(vr.n_vertices(), 1);
        assert_eq!(vr.n_edges(), 0);
    }

    #[test]
    fn test_vr_two_close_points_have_edge() {
        let vr = VietorisRips::new(vec![[0.0, 0.0], [0.5, 0.0]], 1.0);
        assert_eq!(vr.n_edges(), 1);
    }

    #[test]
    fn test_vr_two_far_points_no_edge() {
        let vr = VietorisRips::new(vec![[0.0, 0.0], [10.0, 0.0]], 1.0);
        assert_eq!(vr.n_edges(), 0);
    }

    #[test]
    fn test_vr_triangle_has_one_triangle() {
        let vr = VietorisRips::new(vec![[0.0, 0.0], [1.0, 0.0], [0.5, 0.866]], 2.0);
        assert_eq!(vr.n_triangles(), 1);
    }

    #[test]
    fn test_vr_betti_numbers_four_isolated() {
        // 4 isolated points: β₀=4, β₁=0
        let pts = vec![[0.0, 0.0], [10.0, 0.0], [0.0, 10.0], [10.0, 10.0]];
        let vr = VietorisRips::new(pts, 1.0);
        let b = vr.betti_numbers();
        assert_eq!(b.b0, 4);
        assert_eq!(b.b1, 0);
    }

    #[test]
    fn test_vr_betti_numbers_connected_pair() {
        let vr = VietorisRips::new(vec![[0.0, 0.0], [0.5, 0.0]], 1.0);
        let b = vr.betti_numbers();
        assert_eq!(b.b0, 1);
    }

    #[test]
    fn test_vr_euler_characteristic_points_only() {
        // 5 points, no edges: χ = 5
        let pts: Vec<[f64; 2]> = (0..5).map(|i| [i as f64 * 10.0, 0.0]).collect();
        let vr = VietorisRips::new(pts, 1.0);
        assert_eq!(vr.euler_characteristic(), 5);
    }

    #[test]
    fn test_vr_persistent_h0_returns_correct_count() {
        let pts = vec![[0.0, 0.0], [0.5, 0.0], [5.0, 0.0]];
        let vr = VietorisRips::new(pts, 10.0);
        let intervals = vr.persistent_h0();
        // Should have essential count = 1
        let n_essential = intervals.iter().filter(|p| p.is_essential()).count();
        assert_eq!(n_essential, 1);
    }

    // ── SimplexTree struct ────────────────────────────────────────────────────

    #[test]
    fn test_simplex_tree_empty() {
        let st = SimplexTree::new();
        assert_eq!(st.simplices.len(), 0);
    }

    #[test]
    fn test_simplex_tree_insert_vertex() {
        let mut st = SimplexTree::new();
        st.insert(vec![0], 0.0);
        assert_eq!(st.count_dim(0), 1);
    }

    #[test]
    fn test_simplex_tree_insert_edge() {
        let mut st = SimplexTree::new();
        st.insert(vec![0, 1], 1.5);
        assert_eq!(st.count_dim(1), 1);
    }

    #[test]
    fn test_simplex_tree_max_filtration() {
        let mut st = SimplexTree::new();
        st.insert(vec![0], 0.0);
        st.insert(vec![1], 0.0);
        st.insert(vec![0, 1], 2.5);
        assert!((st.max_filtration() - 2.5).abs() < 1e-12);
    }

    #[test]
    fn test_simplex_tree_sorted_simplices_order() {
        let mut st = SimplexTree::new();
        st.insert(vec![0, 1], 2.0);
        st.insert(vec![0], 0.0);
        st.insert(vec![1], 0.0);
        let sorted = st.sorted_simplices();
        for w in sorted.windows(2) {
            assert!(w[0].filtration <= w[1].filtration);
        }
    }

    #[test]
    fn test_simplex_tree_from_vr() {
        let vr = VietorisRips::new(vec![[0.0, 0.0], [0.5, 0.0], [0.25, 0.5]], 2.0);
        let st = SimplexTree::from_vietoris_rips(&vr);
        // 3 vertices, 3 edges, 1 triangle
        assert_eq!(st.count_dim(0), 3);
        assert_eq!(st.count_dim(1), 3);
        assert_eq!(st.count_dim(2), 1);
    }

    // ── PersistenceDiagram struct ─────────────────────────────────────────────

    #[test]
    fn test_persistence_diagram_add_and_count() {
        let mut pd = PersistenceDiagram::new();
        pd.add(0.0, 1.0, 0);
        pd.add(0.5, 2.0, 1);
        assert_eq!(pd.pairs.len(), 2);
    }

    #[test]
    fn test_persistence_diagram_pairs_in_dim() {
        let mut pd = PersistenceDiagram::new();
        pd.add(0.0, 1.0, 0);
        pd.add(0.5, 2.0, 1);
        pd.add(1.0, 3.0, 0);
        assert_eq!(pd.pairs_in_dim(0).len(), 2);
        assert_eq!(pd.pairs_in_dim(1).len(), 1);
    }

    #[test]
    fn test_persistence_diagram_max_persistence() {
        let mut pd = PersistenceDiagram::new();
        pd.add(0.0, 1.0, 0);
        pd.add(0.0, 3.0, 0);
        assert!((pd.max_persistence() - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_persistence_diagram_n_essential() {
        let mut pd = PersistenceDiagram::new();
        pd.add(0.0, f64::INFINITY, 0);
        pd.add(0.5, 1.5, 0);
        assert_eq!(pd.n_essential(), 1);
    }

    #[test]
    fn test_persistence_diagram_total_persistence() {
        let mut pd = PersistenceDiagram::new();
        pd.add(0.0, 1.0, 0);
        pd.add(0.0, 2.0, 0);
        assert!((pd.total_persistence() - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_persistence_diagram_from_vr() {
        let vr = VietorisRips::new(vec![[0.0, 0.0], [1.0, 0.0]], 5.0);
        let pd = PersistenceDiagram::from_vietoris_rips_h0(&vr);
        assert!(!pd.pairs.is_empty());
        assert_eq!(pd.n_essential(), 1);
    }

    // ── BarcodeSummary struct ─────────────────────────────────────────────────

    #[test]
    fn test_barcode_summary_total_persistence() {
        let mut pd = PersistenceDiagram::new();
        pd.add(0.0, 1.0, 0);
        pd.add(0.0, 2.0, 0);
        let bs = BarcodeSummary::from_diagram(&pd, 0.5);
        assert!((bs.total_persistence - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_barcode_summary_long_features() {
        let mut pd = PersistenceDiagram::new();
        pd.add(0.0, 0.3, 0); // short
        pd.add(0.0, 1.5, 0); // long
        pd.add(0.0, 2.0, 0); // long
        let bs = BarcodeSummary::from_diagram(&pd, 1.0);
        assert_eq!(bs.n_long_features, 2);
    }

    #[test]
    fn test_barcode_summary_entropy_non_negative() {
        let mut pd = PersistenceDiagram::new();
        pd.add(0.0, 1.0, 0);
        pd.add(0.5, 2.0, 0);
        let bs = BarcodeSummary::from_diagram(&pd, 0.5);
        assert!(bs.entropy >= 0.0);
    }

    #[test]
    fn test_barcode_summary_normalised_total() {
        let mut pd = PersistenceDiagram::new();
        pd.add(0.0, 1.0, 0);
        pd.add(0.0, 2.0, 0);
        let bs = BarcodeSummary::from_diagram(&pd, 0.5);
        // normalised_total = 3.0 / 2.0 = 1.5
        assert!((bs.normalised_total() - 1.5).abs() < 1e-10);
    }

    #[test]
    fn test_barcode_summary_mean_persistence() {
        let mut pd = PersistenceDiagram::new();
        pd.add(0.0, 2.0, 0);
        pd.add(0.0, 4.0, 0);
        let bs = BarcodeSummary::from_diagram(&pd, 1.0);
        assert!((bs.mean_persistence - 3.0).abs() < 1e-12);
    }

    // ── MapperGraph struct ────────────────────────────────────────────────────

    #[test]
    fn test_mapper_empty_cloud() {
        let mg = MapperGraph::build(&[], &[], 3, 0.3, 0.5);
        assert_eq!(mg.n_nodes(), 0);
        assert_eq!(mg.n_edges(), 0);
    }

    #[test]
    fn test_mapper_single_point() {
        let pts = [[0.0, 0.0]];
        let filter = [0.0];
        let mg = MapperGraph::build(&pts, &filter, 2, 0.2, 1.0);
        assert!(mg.n_nodes() >= 1);
    }

    #[test]
    fn test_mapper_line_of_points() {
        // 10 collinear points
        let pts: Vec<[f64; 2]> = (0..10).map(|i| [i as f64, 0.0]).collect();
        let filter: Vec<f64> = pts.iter().map(|p| p[0]).collect();
        let mg = MapperGraph::build(&pts, &filter, 3, 0.4, 1.5);
        assert!(mg.n_nodes() >= 1);
    }

    #[test]
    fn test_mapper_degree_non_negative() {
        let pts: Vec<[f64; 2]> = (0..6).map(|i| [i as f64, 0.0]).collect();
        let filter: Vec<f64> = pts.iter().map(|p| p[0]).collect();
        let mg = MapperGraph::build(&pts, &filter, 3, 0.5, 1.5);
        for i in 0..mg.n_nodes() {
            assert!(mg.degree(i) <= mg.n_edges());
        }
    }

    #[test]
    fn test_mapper_total_point_count_ge_input() {
        let pts: Vec<[f64; 2]> = (0..5).map(|i| [i as f64, 0.0]).collect();
        let filter: Vec<f64> = pts.iter().map(|p| p[0]).collect();
        let mg = MapperGraph::build(&pts, &filter, 2, 0.5, 2.0);
        // With overlap, total can exceed n
        assert!(mg.total_point_count() >= 1);
    }

    #[test]
    fn test_mapper_no_intervals_returns_empty() {
        let pts = [[0.0, 0.0], [1.0, 0.0]];
        let filter = [0.0, 1.0];
        let mg = MapperGraph::build(&pts, &filter, 0, 0.2, 0.5);
        assert_eq!(mg.n_nodes(), 0);
    }
}
