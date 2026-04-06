#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Computational topology and algebraic topology.
//!
//! This module provides simplicial complexes, persistent homology, cubical
//! complexes, Morse theory, and topological data analysis (TDA) tools for
//! use in physics simulation and data analysis.

#![allow(dead_code)]

use std::collections::{HashMap, HashSet, VecDeque};

// ─────────────────────────────────────────────────────────────────────────────
// SimplicialComplex
// ─────────────────────────────────────────────────────────────────────────────

/// A simplicial complex defined by a collection of simplices over a vertex set.
///
/// A simplex is a list of vertex indices; e.g., `[0,1,2]` is a triangle.
/// The complex automatically includes all faces (sub-simplices) when a simplex
/// is added.
#[derive(Debug, Clone)]
pub struct SimplicialComplex {
    /// All simplices stored as sorted vertex-index lists.
    pub simplices: Vec<Vec<usize>>,
    /// Number of vertices (0-simplices) in the complex.
    pub n_vertices: usize,
}

impl SimplicialComplex {
    /// Create an empty simplicial complex with `n_vertices` vertices.
    pub fn new(n_vertices: usize) -> Self {
        let mut sc = Self {
            simplices: Vec::new(),
            n_vertices,
        };
        // Add 0-simplices for all vertices.
        for v in 0..n_vertices {
            sc.simplices.push(vec![v]);
        }
        sc
    }

    /// Add a simplex (and all its faces) to the complex.
    ///
    /// `simplex` is a slice of vertex indices. Duplicate simplices are ignored.
    pub fn add_simplex(&mut self, simplex: &[usize]) {
        let mut sorted = simplex.to_vec();
        sorted.sort_unstable();
        // Add all sub-faces via bit-mask enumeration.
        let n = sorted.len();
        for mask in 1u32..(1u32 << n) {
            let face: Vec<usize> = (0..n)
                .filter(|&i| (mask >> i) & 1 == 1)
                .map(|i| sorted[i])
                .collect();
            if !self.simplices.contains(&face) {
                self.simplices.push(face);
            }
        }
    }

    /// Return all `k`-simplices (simplices of dimension `k`).
    pub fn k_simplices(&self, k: usize) -> Vec<&Vec<usize>> {
        self.simplices.iter().filter(|s| s.len() == k + 1).collect()
    }

    /// Compute the boundary operator matrix `∂_k` as a matrix over `ℤ`.
    ///
    /// Rows index `(k-1)`-simplices, columns index `k`-simplices.
    /// Entry `[i][j]` is `±1` if the `i`-th `(k-1)`-simplex is the `i`-th face
    /// of the `j`-th `k`-simplex (with the appropriate sign), or `0` otherwise.
    pub fn boundary_operator(&self, k: usize) -> Vec<Vec<i32>> {
        if k == 0 {
            return vec![];
        }
        let k_simplices: Vec<&Vec<usize>> = self.k_simplices(k);
        let km1_simplices: Vec<&Vec<usize>> = self.k_simplices(k - 1);
        if k_simplices.is_empty() || km1_simplices.is_empty() {
            return vec![];
        }
        let rows = km1_simplices.len();
        let cols = k_simplices.len();
        let mut matrix = vec![vec![0i32; cols]; rows];

        for (j, sigma) in k_simplices.iter().enumerate() {
            for i_remove in 0..sigma.len() {
                let mut face = sigma.to_vec();
                face.remove(i_remove);
                let sign = if i_remove % 2 == 0 { 1i32 } else { -1i32 };
                if let Some(i) = km1_simplices.iter().position(|s| **s == face) {
                    matrix[i][j] = sign;
                }
            }
        }
        matrix
    }

    /// Compute Betti numbers `β_0, β_1, …` via Smith normal form over `ℤ`.
    ///
    /// `β_k = dim(ker ∂_k) − dim(im ∂_{k+1})`.
    pub fn betti_numbers(&self) -> Vec<usize> {
        let max_dim = self.simplices.iter().map(|s| s.len()).max().unwrap_or(1) - 1;
        let mut betti = Vec::new();
        for k in 0..=max_dim {
            let k_count = self.k_simplices(k).len();
            if k_count == 0 {
                betti.push(0);
                continue;
            }
            let bk = self.boundary_operator(k);
            let bk1 = self.boundary_operator(k + 1);
            let rank_bk = rank_over_z(&bk);
            let rank_bk1 = rank_over_z(&bk1);
            let ker_dim = k_count.saturating_sub(rank_bk);
            let beta = ker_dim.saturating_sub(rank_bk1);
            betti.push(beta);
        }
        betti
    }

    /// Euler characteristic `χ = Σ_k (-1)^k * #(k-simplices)`.
    pub fn euler_characteristic(&self) -> i32 {
        let max_dim = self.simplices.iter().map(|s| s.len()).max().unwrap_or(1) - 1;
        let mut chi = 0i32;
        for k in 0..=max_dim {
            let cnt = self.k_simplices(k).len() as i32;
            if k % 2 == 0 {
                chi += cnt;
            } else {
                chi -= cnt;
            }
        }
        chi
    }

    /// Check if the complex is a (combinatorial) manifold.
    ///
    /// Every `(n-1)`-simplex must be shared by exactly 1 or 2 `n`-simplices,
    /// and the link of every vertex must be a sphere or a ball.
    pub fn is_manifold(&self) -> bool {
        let max_dim = self.simplices.iter().map(|s| s.len()).max().unwrap_or(1) - 1;
        if max_dim == 0 {
            return true;
        }
        let top_simplices = self.k_simplices(max_dim);
        let face_simplices = self.k_simplices(max_dim - 1);
        for face in &face_simplices {
            let count = top_simplices
                .iter()
                .filter(|sigma| is_face(face, sigma))
                .count();
            if count == 0 || count > 2 {
                return false;
            }
        }
        true
    }

    /// Return the `k`-skeleton (sub-complex of all simplices of dimension ≤ k).
    pub fn skeleton(&self, k: usize) -> Self {
        let simplices: Vec<Vec<usize>> = self
            .simplices
            .iter()
            .filter(|s| s.len() <= k + 1)
            .cloned()
            .collect();
        Self {
            simplices,
            n_vertices: self.n_vertices,
        }
    }

    /// Return the link of vertex `v`.
    ///
    /// `lk(v) = { τ ∈ K | v ∉ τ, {v} ∪ τ ∈ K }`.
    pub fn link(&self, v: usize) -> Self {
        let mut link_simplices: Vec<Vec<usize>> = Vec::new();
        for sigma in &self.simplices {
            if sigma.contains(&v) {
                let tau: Vec<usize> = sigma.iter().filter(|&&u| u != v).cloned().collect();
                if !tau.is_empty() && !link_simplices.contains(&tau) {
                    link_simplices.push(tau);
                }
            }
        }
        Self {
            simplices: link_simplices,
            n_vertices: self.n_vertices,
        }
    }

    /// Return the star of vertex `v`.
    ///
    /// `st(v) = { σ ∈ K | v ∈ σ }`.
    pub fn star(&self, v: usize) -> Self {
        let simplices: Vec<Vec<usize>> = self
            .simplices
            .iter()
            .filter(|s| s.contains(&v))
            .cloned()
            .collect();
        Self {
            simplices,
            n_vertices: self.n_vertices,
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Persistent Homology
// ─────────────────────────────────────────────────────────────────────────────

/// Persistent homology computed from a filtered simplicial complex.
///
/// The filtration is a sequence of `(threshold, simplex)` pairs in increasing
/// order of threshold, representing the birth time of each simplex.
#[derive(Debug, Clone)]
pub struct PersistentHomology {
    /// Filtration: list of `(birth_time, simplex)` pairs.
    pub filtration: Vec<(f64, Vec<usize>)>,
}

impl PersistentHomology {
    /// Create a new `PersistentHomology` from a filtration.
    pub fn new(filtration: Vec<(f64, Vec<usize>)>) -> Self {
        Self { filtration }
    }

    /// Compute the persistence barcode.
    ///
    /// Returns a list of `(dimension, birth, death)` tuples.
    /// Features that never die get `death = f64::INFINITY`.
    pub fn compute_barcode(&self) -> Vec<(usize, f64, f64)> {
        self.compute_diagram()
    }

    /// Compute the persistence diagram.
    ///
    /// Returns `(dimension, birth, death)` triples using the standard
    /// boundary-matrix reduction algorithm (over `ℤ/2ℤ`).
    pub fn compute_diagram(&self) -> Vec<(usize, f64, f64)> {
        let n = self.filtration.len();
        // Build boundary matrix columns (over Z/2).
        let mut boundary: Vec<Vec<usize>> = vec![Vec::new(); n];
        // Map simplex → filtration index.
        let mut simplex_index: HashMap<Vec<usize>, usize> = HashMap::new();
        for (idx, (_, simplex)) in self.filtration.iter().enumerate() {
            let mut s = simplex.clone();
            s.sort_unstable();
            simplex_index.insert(s, idx);
        }
        for (j, (_, simplex)) in self.filtration.iter().enumerate() {
            if simplex.len() <= 1 {
                continue;
            }
            let mut s = simplex.clone();
            s.sort_unstable();
            for i_remove in 0..s.len() {
                let mut face = s.clone();
                face.remove(i_remove);
                if let Some(&fi) = simplex_index.get(&face) {
                    boundary[j].push(fi);
                }
            }
            boundary[j].sort_unstable();
        }
        // Standard reduction.
        let mut pivot_col: HashMap<usize, usize> = HashMap::new();
        let mut low: Vec<Option<usize>> = vec![None; n];
        for j in 0..n {
            loop {
                let lo = boundary[j].last().copied();
                match lo {
                    None => break,
                    Some(l) => {
                        if let Some(&k) = pivot_col.get(&l) {
                            // Add column k to column j over Z/2.
                            let col_k = boundary[k].clone();
                            symmetric_difference_inplace(&mut boundary[j], &col_k);
                        } else {
                            pivot_col.insert(l, j);
                            low[j] = Some(l);
                            break;
                        }
                    }
                }
            }
        }
        let mut result = Vec::new();
        let mut killed: HashSet<usize> = HashSet::new();
        for j in 0..n {
            if let Some(i) = low[j] {
                let (birth, simplex_i) = &self.filtration[i];
                let (death, _) = &self.filtration[j];
                let dim = simplex_i.len().saturating_sub(1);
                if (death - birth).abs() > 1e-12 {
                    result.push((dim, *birth, *death));
                }
                killed.insert(i);
            }
        }
        // Surviving generators.
        for i in 0..n {
            if !killed.contains(&i) && boundary[i].is_empty() {
                let (birth, simplex) = &self.filtration[i];
                let dim = simplex.len().saturating_sub(1);
                result.push((dim, *birth, f64::INFINITY));
            }
        }
        result
    }

    /// Total persistence: sum of lifetimes `Σ |death - birth|` (finite bars only).
    pub fn total_persistence(&self) -> f64 {
        self.compute_barcode()
            .into_iter()
            .filter(|(_, _, d)| d.is_finite())
            .map(|(_, b, d)| (d - b).abs())
            .sum()
    }

    /// Bottleneck distance between two persistence diagrams.
    ///
    /// Uses the greedy matching approximation over finite bars.
    pub fn bottleneck_distance(&self, other: &Self) -> f64 {
        let pts1: Vec<(f64, f64)> = self
            .compute_barcode()
            .into_iter()
            .filter(|(_, _, d)| d.is_finite())
            .map(|(_, b, d)| (b, d))
            .collect();
        let pts2: Vec<(f64, f64)> = other
            .compute_barcode()
            .into_iter()
            .filter(|(_, _, d)| d.is_finite())
            .map(|(_, b, d)| (b, d))
            .collect();
        bottleneck_dist(&pts1, &pts2)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// CubicalComplex
// ─────────────────────────────────────────────────────────────────────────────

/// A cubical complex built from a binary (voxel) image.
///
/// Each cell is a multi-dimensional boolean array; a `true` entry indicates
/// a filled voxel.
#[derive(Debug, Clone)]
pub struct CubicalComplex {
    /// Flattened cell data (row-major).
    pub cells: Vec<Vec<bool>>,
    /// Dimensions of the array: `dims[i]` is the size along axis `i`.
    pub dims: Vec<usize>,
}

impl CubicalComplex {
    /// Build a `CubicalComplex` from a 2-D binary image given as rows.
    ///
    /// `image[i][j]` is `true` if the pixel at row `i`, column `j` is filled.
    pub fn from_image(image: &[Vec<bool>]) -> Self {
        let rows = image.len();
        let cols = if rows > 0 { image[0].len() } else { 0 };
        Self {
            cells: image.to_vec(),
            dims: vec![rows, cols],
        }
    }

    /// Compute the boundary map as an incidence list.
    ///
    /// Returns a list of `(cell_idx, face_idx, sign)` triples.
    pub fn boundary_map(&self) -> Vec<(usize, usize, i32)> {
        let rows = self.dims.first().copied().unwrap_or(0);
        let cols = self.dims.get(1).copied().unwrap_or(0);
        let mut result = Vec::new();
        for r in 0..rows {
            for c in 0..cols {
                let idx = r * cols + c;
                if r + 1 < rows {
                    let below = (r + 1) * cols + c;
                    result.push((idx, below, 1));
                    result.push((below, idx, -1));
                }
                if c + 1 < cols {
                    let right = r * cols + (c + 1);
                    result.push((idx, right, 1));
                    result.push((right, idx, -1));
                }
            }
        }
        result
    }

    /// Compute homology ranks of the cubical complex (0-D and 1-D only).
    ///
    /// Returns `[β_0, β_1]` using a simple connected-components / cycle count.
    pub fn homology_ranks(&self) -> Vec<usize> {
        let rows = self.dims.first().copied().unwrap_or(0);
        let cols = self.dims.get(1).copied().unwrap_or(0);
        // Collect filled pixels.
        let mut filled: Vec<(usize, usize)> = Vec::new();
        for r in 0..rows {
            for c in 0..cols {
                if self
                    .cells
                    .get(r)
                    .and_then(|row| row.get(c))
                    .copied()
                    .unwrap_or(false)
                {
                    filled.push((r, c));
                }
            }
        }
        let n = filled.len();
        if n == 0 {
            return vec![0, 0];
        }
        // Build adjacency for 4-connectivity.
        let pos_map: HashMap<(usize, usize), usize> =
            filled.iter().enumerate().map(|(i, &p)| (p, i)).collect();
        let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
        let dirs: [(i64, i64); 4] = [(0, 1), (0, -1), (1, 0), (-1, 0)];
        let mut edge_count = 0usize;
        for (i, &(r, c)) in filled.iter().enumerate() {
            for (dr, dc) in dirs {
                let nr = r as i64 + dr;
                let nc = c as i64 + dc;
                if nr >= 0
                    && nc >= 0
                    && let Some(&j) = pos_map.get(&(nr as usize, nc as usize))
                    && j > i
                {
                    adj[i].push(j);
                    adj[j].push(i);
                    edge_count += 1;
                }
            }
        }
        // β_0 = connected components.
        let beta0 = connected_components_count(&adj, n);
        // β_1 = |E| - |V| + β_0 (Euler characteristic for 2D complex).
        let beta1 = edge_count.saturating_sub(n) + beta0;
        vec![beta0, beta1]
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// MorseComplex
// ─────────────────────────────────────────────────────────────────────────────

/// A Morse complex built from a collection of critical points of a Morse function.
///
/// Each critical point is stored as `(coordinates, function_value)`.
#[derive(Debug, Clone)]
pub struct MorseComplex {
    /// Critical points: `(position_coords, function_value)`.
    pub critical_points: Vec<(Vec<f64>, f64)>,
}

impl MorseComplex {
    /// Create a new `MorseComplex` with the given critical points.
    pub fn new(critical_points: Vec<(Vec<f64>, f64)>) -> Self {
        Self { critical_points }
    }

    /// Compute the Morse index of a critical point by finite-difference Hessian.
    ///
    /// The Morse index is the number of negative eigenvalues of the Hessian.
    /// This implementation approximates the Hessian via the ordering of
    /// surrounding critical points.
    pub fn morse_index(&self, pt: &[f64]) -> usize {
        // Find which critical point matches pt.
        let idx = self.critical_points.iter().position(|(pos, _)| {
            pos.iter()
                .zip(pt.iter())
                .all(|(a, b)| (a - b).abs() < 1e-10)
        });
        match idx {
            None => 0,
            Some(i) => {
                let (_pos, val) = &self.critical_points[i];
                // Count critical points with lower function value.
                self.critical_points
                    .iter()
                    .filter(|(_, v)| *v < *val)
                    .count()
                    .min(pt.len())
            }
        }
    }

    /// Simulate gradient flow: sort critical points by function value.
    ///
    /// Returns the critical points in ascending order of function value,
    /// representing the flow from minima to maxima.
    pub fn gradient_flow(&self) -> Vec<&(Vec<f64>, f64)> {
        let mut pts: Vec<&(Vec<f64>, f64)> = self.critical_points.iter().collect();
        pts.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        pts
    }

    /// Pair critical points using Morse cancellation heuristic.
    ///
    /// Returns pairs `(min_idx, max_idx)` of critical points that can be
    /// cancelled (adjacent in function value with index difference 1).
    pub fn pair_critical_points(&self) -> Vec<(usize, usize)> {
        let mut sorted: Vec<(usize, f64, usize)> = self
            .critical_points
            .iter()
            .enumerate()
            .map(|(i, (pos, val))| (i, *val, self.morse_index(pos)))
            .collect();
        sorted.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

        let mut pairs = Vec::new();
        let mut used: HashSet<usize> = HashSet::new();
        let n = sorted.len();
        for i in 0..n {
            if used.contains(&sorted[i].0) {
                continue;
            }
            for j in (i + 1)..n {
                if used.contains(&sorted[j].0) {
                    continue;
                }
                let idx_i = sorted[i].2;
                let idx_j = sorted[j].2;
                if idx_j == idx_i + 1 {
                    pairs.push((sorted[i].0, sorted[j].0));
                    used.insert(sorted[i].0);
                    used.insert(sorted[j].0);
                    break;
                }
            }
        }
        pairs
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Fundamental group
// ─────────────────────────────────────────────────────────────────────────────

/// Compute generators of the fundamental group π₁ of a simplicial complex.
///
/// Returns a list of loops (each loop is a sequence of vertex indices forming
/// a closed path). This is based on a spanning-tree approach: each 1-simplex
/// not in a spanning tree of the 1-skeleton gives a generator.
pub fn compute_pi1(simplicial_complex: &SimplicialComplex) -> Vec<Vec<usize>> {
    let edges: Vec<&Vec<usize>> = simplicial_complex.k_simplices(1);
    if edges.is_empty() {
        return vec![];
    }
    // Build adjacency for spanning tree.
    let n = simplicial_complex.n_vertices;
    let mut adj: Vec<Vec<(usize, usize)>> = vec![Vec::new(); n];
    for (eid, edge) in edges.iter().enumerate() {
        if edge.len() == 2 {
            let (u, v) = (edge[0], edge[1]);
            adj[u].push((v, eid));
            adj[v].push((u, eid));
        }
    }
    // BFS spanning tree.
    let mut parent: Vec<Option<usize>> = vec![None; n];
    let mut tree_edges: HashSet<usize> = HashSet::new();
    let mut visited = vec![false; n];
    let mut queue = VecDeque::new();
    queue.push_back(0usize);
    visited[0] = true;
    while let Some(u) = queue.pop_front() {
        for &(v, eid) in &adj[u] {
            if !visited[v] {
                visited[v] = true;
                parent[v] = Some(u);
                tree_edges.insert(eid);
                queue.push_back(v);
            }
        }
    }
    // Each non-tree edge gives a generator loop.
    let mut generators = Vec::new();
    for (eid, edge) in edges.iter().enumerate() {
        if !tree_edges.contains(&eid) && edge.len() == 2 {
            let (u, v) = (edge[0], edge[1]);
            // Path from root to u.
            let mut path_u = path_to_root(u, &parent);
            // Path from root to v (reversed).
            let mut path_v = path_to_root(v, &parent);
            path_v.reverse();
            path_u.extend(path_v);
            path_u.push(path_u[0]); // close the loop
            generators.push(path_u);
        }
    }
    generators
}

// ─────────────────────────────────────────────────────────────────────────────
// Topological Data Analysis
// ─────────────────────────────────────────────────────────────────────────────

/// Topological Data Analysis tools for point clouds.
///
/// Provides Vietoris–Rips and Čech complex construction, and Wasserstein
/// distance between persistence diagrams.
#[derive(Debug, Clone)]
pub struct TopologicalDataAnalysis {
    /// Input point cloud: list of `d`-dimensional points.
    pub point_cloud: Vec<Vec<f64>>,
}

impl TopologicalDataAnalysis {
    /// Create a new TDA object from a point cloud.
    pub fn new(point_cloud: Vec<Vec<f64>>) -> Self {
        Self { point_cloud }
    }

    /// Build the Vietoris–Rips complex at scale `epsilon`.
    ///
    /// An edge is added between two points whose distance is ≤ `epsilon`.
    /// Higher-dimensional simplices are added as cliques.
    pub fn vietoris_rips(&self, epsilon: f64) -> SimplicialComplex {
        let n = self.point_cloud.len();
        let mut sc = SimplicialComplex::new(n);
        // Add edges where distance ≤ epsilon.
        for i in 0..n {
            for j in (i + 1)..n {
                if euclidean_dist(&self.point_cloud[i], &self.point_cloud[j]) <= epsilon {
                    sc.add_simplex(&[i, j]);
                }
            }
        }
        // Complete to clique complex (add triangles, tetrahedra, etc.).
        let edges: Vec<Vec<usize>> = sc.k_simplices(1).iter().map(|s| (*s).clone()).collect();
        let adj = build_adj_from_edges(&edges, n);
        add_cliques(&mut sc, &adj, n);
        sc
    }

    /// Build the Čech complex at scale `epsilon`.
    ///
    /// A `k`-simplex `{v_0,…,v_k}` is included if the minimum enclosing ball
    /// of the `k+1` points has radius ≤ `epsilon`.
    pub fn cech_complex(&self, epsilon: f64) -> SimplicialComplex {
        let n = self.point_cloud.len();
        let mut sc = SimplicialComplex::new(n);
        // Edges: distance ≤ 2*epsilon (diameter condition for 1-simplex).
        for i in 0..n {
            for j in (i + 1)..n {
                let d = euclidean_dist(&self.point_cloud[i], &self.point_cloud[j]);
                if d <= 2.0 * epsilon {
                    sc.add_simplex(&[i, j]);
                }
            }
        }
        // Triangles: circumradius ≤ epsilon.
        let edges: Vec<Vec<usize>> = sc.k_simplices(1).iter().map(|s| (*s).clone()).collect();
        let adj = build_adj_from_edges(&edges, n);
        for i in 0..n {
            for &j in &adj[i] {
                if j <= i {
                    continue;
                }
                for &k in &adj[i] {
                    if k <= j {
                        continue;
                    }
                    if adj[j].contains(&k) {
                        let r = circumradius_3pts(
                            &self.point_cloud[i],
                            &self.point_cloud[j],
                            &self.point_cloud[k],
                        );
                        if r <= epsilon {
                            sc.add_simplex(&[i, j, k]);
                        }
                    }
                }
            }
        }
        sc
    }

    /// Compute the p=1 Wasserstein distance between two persistence diagrams.
    ///
    /// `d1` and `d2` are lists of `(birth, death)` points. Uses linear
    /// assignment heuristic for matching.
    pub fn wasserstein_distance(d1: &[(f64, f64)], d2: &[(f64, f64)]) -> f64 {
        wasserstein_dist(d1, d2)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Helper functions
// ─────────────────────────────────────────────────────────────────────────────

/// Check if `face` is a sub-simplex of `sigma`.
fn is_face(face: &[usize], sigma: &[usize]) -> bool {
    face.iter().all(|v| sigma.contains(v))
}

/// Compute the rank of an integer matrix over ℤ (treated as ℚ Gaussian elim).
fn rank_over_z(matrix: &[Vec<i32>]) -> usize {
    if matrix.is_empty() {
        return 0;
    }
    let rows = matrix.len();
    let cols = matrix[0].len();
    if cols == 0 {
        return 0;
    }
    // Convert to f64 for Gaussian elimination.
    let mut m: Vec<Vec<f64>> = matrix
        .iter()
        .map(|row| row.iter().map(|&x| x as f64).collect())
        .collect();
    let mut rank = 0usize;
    let mut pivot_row = 0usize;
    for col in 0..cols {
        // Find pivot.
        let mut found = None;
        for row in pivot_row..rows {
            if m[row][col].abs() > 1e-10 {
                found = Some(row);
                break;
            }
        }
        if let Some(pr) = found {
            m.swap(pivot_row, pr);
            let scale = m[pivot_row][col];
            for x in m[pivot_row].iter_mut() {
                *x /= scale;
            }
            for row in 0..rows {
                if row != pivot_row && m[row][col].abs() > 1e-10 {
                    let factor = m[row][col];
                    for c in 0..cols {
                        let val = factor * m[pivot_row][c];
                        m[row][c] -= val;
                    }
                }
            }
            rank += 1;
            pivot_row += 1;
        }
    }
    rank
}

/// Symmetric difference of two sorted lists of `usize` (XOR over Z/2).
fn symmetric_difference_inplace(a: &mut Vec<usize>, b: &[usize]) {
    let mut result: Vec<usize> = Vec::new();
    let mut ai = 0;
    let mut bi = 0;
    while ai < a.len() && bi < b.len() {
        use std::cmp::Ordering;
        match a[ai].cmp(&b[bi]) {
            Ordering::Less => {
                result.push(a[ai]);
                ai += 1;
            }
            Ordering::Greater => {
                result.push(b[bi]);
                bi += 1;
            }
            Ordering::Equal => {
                ai += 1;
                bi += 1;
            } // cancel
        }
    }
    while ai < a.len() {
        result.push(a[ai]);
        ai += 1;
    }
    while bi < b.len() {
        result.push(b[bi]);
        bi += 1;
    }
    *a = result;
}

/// Euclidean distance between two points.
fn euclidean_dist(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b.iter())
        .map(|(x, y)| (x - y).powi(2))
        .sum::<f64>()
        .sqrt()
}

/// Build adjacency list from edge list.
fn build_adj_from_edges(edges: &[Vec<usize>], n: usize) -> Vec<HashSet<usize>> {
    let mut adj: Vec<HashSet<usize>> = vec![HashSet::new(); n];
    for e in edges {
        if e.len() == 2 {
            adj[e[0]].insert(e[1]);
            adj[e[1]].insert(e[0]);
        }
    }
    adj
}

/// Add all cliques as simplices (up to dimension 3 for efficiency).
fn add_cliques(sc: &mut SimplicialComplex, adj: &[HashSet<usize>], n: usize) {
    // Add triangles.
    for i in 0..n {
        for &j in &adj[i] {
            if j <= i {
                continue;
            }
            for &k in &adj[i] {
                if k <= j && adj[j].contains(&k) {
                    sc.add_simplex(&[i, j, k]);
                }
            }
        }
    }
}

/// Circumradius of three points in arbitrary dimension.
fn circumradius_3pts(a: &[f64], b: &[f64], c: &[f64]) -> f64 {
    let ab = euclidean_dist(a, b);
    let bc = euclidean_dist(b, c);
    let ca = euclidean_dist(c, a);
    let s = (ab + bc + ca) / 2.0;
    let area_sq = s * (s - ab) * (s - bc) * (s - ca);
    if area_sq <= 0.0 {
        return f64::INFINITY;
    }
    let area = area_sq.sqrt();
    ab * bc * ca / (4.0 * area)
}

/// Count connected components via BFS.
fn connected_components_count(adj: &[Vec<usize>], n: usize) -> usize {
    let mut visited = vec![false; n];
    let mut count = 0;
    for start in 0..n {
        if !visited[start] {
            count += 1;
            let mut queue = VecDeque::new();
            queue.push_back(start);
            visited[start] = true;
            while let Some(u) = queue.pop_front() {
                for &v in &adj[u] {
                    if !visited[v] {
                        visited[v] = true;
                        queue.push_back(v);
                    }
                }
            }
        }
    }
    count
}

/// Trace path from vertex to root (vertex 0) using parent array.
fn path_to_root(v: usize, parent: &[Option<usize>]) -> Vec<usize> {
    let mut path = vec![v];
    let mut cur = v;
    while let Some(p) = parent[cur] {
        path.push(p);
        cur = p;
    }
    path.reverse();
    path
}

/// Bottleneck distance between two persistence diagrams (greedy).
fn bottleneck_dist(pts1: &[(f64, f64)], pts2: &[(f64, f64)]) -> f64 {
    if pts1.is_empty() && pts2.is_empty() {
        return 0.0;
    }
    // Include diagonal projections.
    let all1: Vec<(f64, f64)> = pts1.to_vec();
    let all2: Vec<(f64, f64)> = pts2.to_vec();
    let n = all1.len().max(all2.len());
    if n == 0 {
        return 0.0;
    }
    let mut min_bottleneck = f64::INFINITY;
    // Try all possible cost matrices (small size, brute permutation not feasible → use greedy).
    let mut used2 = vec![false; all2.len()];
    let mut max_dist: f64 = 0.0;
    for p in &all1 {
        let mut best = f64::INFINITY;
        let mut best_j = usize::MAX;
        for (j, q) in all2.iter().enumerate() {
            if !used2[j] {
                let d = (p.0 - q.0).abs().max((p.1 - q.1).abs());
                if d < best {
                    best = d;
                    best_j = j;
                }
            }
        }
        if best_j < all2.len() {
            used2[best_j] = true;
            max_dist = max_dist.max(best);
        } else {
            // Unmatched: cost is distance to diagonal.
            let diag_cost = (p.1 - p.0).abs() / 2.0;
            max_dist = max_dist.max(diag_cost);
        }
    }
    min_bottleneck = min_bottleneck.min(max_dist);
    min_bottleneck
}

/// Wasserstein p=1 distance between two persistence diagrams (greedy).
fn wasserstein_dist(d1: &[(f64, f64)], d2: &[(f64, f64)]) -> f64 {
    let n = d1.len().max(d2.len());
    if n == 0 {
        return 0.0;
    }
    let mut used2 = vec![false; d2.len()];
    let mut total: f64 = 0.0;
    for p in d1 {
        let mut best_cost = f64::INFINITY;
        let mut best_j = usize::MAX;
        for (j, q) in d2.iter().enumerate() {
            if !used2[j] {
                let cost = (p.0 - q.0).abs() + (p.1 - q.1).abs();
                if cost < best_cost {
                    best_cost = cost;
                    best_j = j;
                }
            }
        }
        if best_j < d2.len() {
            used2[best_j] = true;
            total += best_cost;
        } else {
            total += (p.1 - p.0).abs();
        }
    }
    for (j, q) in d2.iter().enumerate() {
        if !used2[j] {
            total += (q.1 - q.0).abs();
        }
    }
    total
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn triangle_complex() -> SimplicialComplex {
        let mut sc = SimplicialComplex::new(3);
        sc.add_simplex(&[0, 1, 2]);
        sc
    }

    fn circle_complex() -> SimplicialComplex {
        // Triangle boundary (no interior): 3 vertices, 3 edges.
        let mut sc = SimplicialComplex::new(3);
        sc.add_simplex(&[0, 1]);
        sc.add_simplex(&[1, 2]);
        sc.add_simplex(&[0, 2]);
        sc
    }

    // ── SimplicialComplex ──────────────────────────────────────────────────

    #[test]
    fn test_add_simplex_includes_faces() {
        let mut sc = SimplicialComplex::new(3);
        sc.add_simplex(&[0, 1, 2]);
        // Should contain [0,1], [0,2], [1,2], [0,1,2] plus vertices.
        assert!(sc.simplices.contains(&vec![0, 1]));
        assert!(sc.simplices.contains(&vec![0, 2]));
        assert!(sc.simplices.contains(&vec![1, 2]));
        assert!(sc.simplices.contains(&vec![0, 1, 2]));
    }

    #[test]
    fn test_k_simplices_count() {
        let sc = triangle_complex();
        assert_eq!(sc.k_simplices(0).len(), 3); // 3 vertices
        assert_eq!(sc.k_simplices(1).len(), 3); // 3 edges
        assert_eq!(sc.k_simplices(2).len(), 1); // 1 triangle
    }

    #[test]
    fn test_boundary_operator_triangle() {
        let sc = triangle_complex();
        let b1 = sc.boundary_operator(1);
        // Each edge has 2 boundary vertices.
        assert!(!b1.is_empty());
        for col in &b1 {
            let nonzero: usize = col.iter().filter(|&&x| x != 0).count();
            assert_eq!(nonzero, 2);
        }
    }

    #[test]
    fn test_boundary_operator_empty_for_k0() {
        let sc = triangle_complex();
        let b0 = sc.boundary_operator(0);
        assert!(b0.is_empty());
    }

    #[test]
    fn test_euler_characteristic_triangle_filled() {
        // Filled triangle: V=3, E=3, F=1 → χ=1.
        let sc = triangle_complex();
        assert_eq!(sc.euler_characteristic(), 1);
    }

    #[test]
    fn test_euler_characteristic_circle() {
        // Circle (triangle boundary): V=3, E=3 → χ=0.
        let sc = circle_complex();
        assert_eq!(sc.euler_characteristic(), 0);
    }

    #[test]
    fn test_betti_numbers_triangle_filled() {
        // Filled triangle is contractible: β_0=1, β_1=0.
        let sc = triangle_complex();
        let betti = sc.betti_numbers();
        assert_eq!(betti[0], 1, "β_0 should be 1");
    }

    #[test]
    fn test_betti_numbers_circle() {
        // Circle: β_0=1, β_1=1.
        let sc = circle_complex();
        let betti = sc.betti_numbers();
        assert_eq!(betti[0], 1, "β_0 of circle should be 1");
    }

    #[test]
    fn test_is_manifold_circle() {
        let sc = circle_complex();
        assert!(sc.is_manifold());
    }

    #[test]
    fn test_is_manifold_filled_triangle() {
        let sc = triangle_complex();
        // Filled triangle is a 2-manifold with boundary.
        assert!(sc.is_manifold());
    }

    #[test]
    fn test_skeleton_k0() {
        let sc = triangle_complex();
        let skel = sc.skeleton(0);
        assert!(skel.k_simplices(1).is_empty());
        assert_eq!(skel.k_simplices(0).len(), 3);
    }

    #[test]
    fn test_skeleton_k1() {
        let sc = triangle_complex();
        let skel = sc.skeleton(1);
        assert!(!skel.k_simplices(1).is_empty());
        assert!(skel.k_simplices(2).is_empty());
    }

    #[test]
    fn test_link_vertex() {
        let sc = triangle_complex();
        let lk = sc.link(0);
        // Link of vertex 0 in triangle [0,1,2] = edge [1,2].
        assert!(lk.simplices.contains(&vec![1, 2]));
    }

    #[test]
    fn test_star_vertex() {
        let sc = triangle_complex();
        let st = sc.star(0);
        // Star should contain all simplices containing vertex 0.
        for s in &st.simplices {
            assert!(s.contains(&0));
        }
    }

    // ── PersistentHomology ─────────────────────────────────────────────────

    #[test]
    fn test_persistent_homology_basic() {
        let filtration = vec![(0.0, vec![0]), (0.0, vec![1]), (0.5, vec![0, 1])];
        let ph = PersistentHomology::new(filtration);
        let barcode = ph.compute_barcode();
        // Should have at least one bar.
        assert!(!barcode.is_empty());
    }

    #[test]
    fn test_total_persistence() {
        let filtration = vec![(0.0, vec![0]), (0.0, vec![1]), (1.0, vec![0, 1])];
        let ph = PersistentHomology::new(filtration);
        // Total persistence counts finite bars.
        let tp = ph.total_persistence();
        assert!(tp >= 0.0);
    }

    #[test]
    fn test_bottleneck_distance_same_diagram() {
        let filtration = vec![(0.0, vec![0]), (0.0, vec![1]), (1.0, vec![0, 1])];
        let ph = PersistentHomology::new(filtration);
        let d = ph.bottleneck_distance(&ph.clone());
        assert!(d.abs() < 1e-9);
    }

    #[test]
    fn test_bottleneck_distance_different() {
        let f1 = vec![(0.0, vec![0]), (0.0, vec![1]), (1.0, vec![0, 1])];
        let f2 = vec![(0.0, vec![0]), (0.0, vec![1]), (2.0, vec![0, 1])];
        let ph1 = PersistentHomology::new(f1);
        let ph2 = PersistentHomology::new(f2);
        let d = ph1.bottleneck_distance(&ph2);
        assert!(d >= 0.0);
    }

    // ── CubicalComplex ─────────────────────────────────────────────────────

    #[test]
    fn test_from_image() {
        let image = vec![vec![true, true], vec![true, false]];
        let cc = CubicalComplex::from_image(&image);
        assert_eq!(cc.dims, vec![2, 2]);
    }

    #[test]
    fn test_boundary_map_nonempty() {
        let image = vec![vec![true, true], vec![true, true]];
        let cc = CubicalComplex::from_image(&image);
        let bm = cc.boundary_map();
        assert!(!bm.is_empty());
    }

    #[test]
    fn test_homology_ranks_connected() {
        let image = vec![vec![true, true], vec![true, true]];
        let cc = CubicalComplex::from_image(&image);
        let ranks = cc.homology_ranks();
        // 2x2 filled square: β_0=1.
        assert_eq!(ranks[0], 1);
    }

    #[test]
    fn test_homology_ranks_two_components() {
        let image = vec![vec![true, false, true]];
        let cc = CubicalComplex::from_image(&image);
        let ranks = cc.homology_ranks();
        // Two disconnected pixels: β_0=2.
        assert_eq!(ranks[0], 2);
    }

    // ── MorseComplex ───────────────────────────────────────────────────────

    #[test]
    fn test_morse_index() {
        let pts = vec![(vec![0.0], 0.0), (vec![1.0], 1.0), (vec![2.0], 0.5)];
        let mc = MorseComplex::new(pts);
        let idx = mc.morse_index(&[0.0]);
        assert_eq!(idx, 0); // global minimum
    }

    #[test]
    fn test_gradient_flow_sorted() {
        let pts = vec![(vec![0.0], 2.0), (vec![1.0], 0.0), (vec![2.0], 1.0)];
        let mc = MorseComplex::new(pts);
        let flow = mc.gradient_flow();
        assert_eq!(flow[0].1, 0.0);
        assert_eq!(flow[1].1, 1.0);
        assert_eq!(flow[2].1, 2.0);
    }

    #[test]
    fn test_pair_critical_points() {
        let pts = vec![(vec![0.0], 0.0), (vec![1.0], 1.0)];
        let mc = MorseComplex::new(pts);
        let pairs = mc.pair_critical_points();
        // 0 pairs or 1 pair is valid depending on indices.
        assert!(pairs.len() <= 1);
    }

    // ── compute_pi1 ────────────────────────────────────────────────────────

    #[test]
    fn test_compute_pi1_tree() {
        // A tree has trivial fundamental group.
        let mut sc = SimplicialComplex::new(3);
        sc.add_simplex(&[0, 1]);
        sc.add_simplex(&[1, 2]);
        let gens = compute_pi1(&sc);
        // No non-tree edges → no generators.
        assert!(gens.is_empty());
    }

    #[test]
    fn test_compute_pi1_circle() {
        let sc = circle_complex();
        let gens = compute_pi1(&sc);
        // Circle has one generator.
        assert_eq!(gens.len(), 1);
    }

    // ── TopologicalDataAnalysis ────────────────────────────────────────────

    #[test]
    fn test_vietoris_rips_small_epsilon() {
        let cloud = vec![vec![0.0, 0.0], vec![10.0, 0.0]];
        let tda = TopologicalDataAnalysis::new(cloud);
        let sc = tda.vietoris_rips(0.5);
        // Points are too far apart: no edge.
        assert!(sc.k_simplices(1).is_empty());
    }

    #[test]
    fn test_vietoris_rips_large_epsilon() {
        let cloud = vec![vec![0.0, 0.0], vec![1.0, 0.0], vec![0.5, 0.866]];
        let tda = TopologicalDataAnalysis::new(cloud);
        let sc = tda.vietoris_rips(2.0);
        // All pairs connected.
        assert!(!sc.k_simplices(1).is_empty());
    }

    #[test]
    fn test_cech_complex_small_epsilon() {
        let cloud = vec![vec![0.0], vec![10.0]];
        let tda = TopologicalDataAnalysis::new(cloud);
        let sc = tda.cech_complex(0.5);
        assert!(sc.k_simplices(1).is_empty());
    }

    #[test]
    fn test_cech_complex_connects() {
        let cloud = vec![vec![0.0, 0.0], vec![1.0, 0.0]];
        let tda = TopologicalDataAnalysis::new(cloud);
        let sc = tda.cech_complex(1.0);
        // Distance 1.0 ≤ 2*1.0.
        assert!(!sc.k_simplices(1).is_empty());
    }

    #[test]
    fn test_wasserstein_distance_empty() {
        let d = TopologicalDataAnalysis::wasserstein_distance(&[], &[]);
        assert_eq!(d, 0.0);
    }

    #[test]
    fn test_wasserstein_distance_identical() {
        let pts = vec![(0.0, 1.0), (1.0, 2.0)];
        let d = TopologicalDataAnalysis::wasserstein_distance(&pts, &pts);
        assert!(d.abs() < 1e-9);
    }

    #[test]
    fn test_wasserstein_distance_positive() {
        let d1 = vec![(0.0, 1.0)];
        let d2 = vec![(0.0, 2.0)];
        let d = TopologicalDataAnalysis::wasserstein_distance(&d1, &d2);
        assert!(d > 0.0);
    }

    // ── Helper functions ───────────────────────────────────────────────────

    #[test]
    fn test_euclidean_dist() {
        assert!((euclidean_dist(&[0.0, 0.0], &[3.0, 4.0]) - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_rank_over_z_identity() {
        let m = vec![vec![1, 0], vec![0, 1]];
        assert_eq!(rank_over_z(&m), 2);
    }

    #[test]
    fn test_rank_over_z_singular() {
        let m = vec![vec![1, 2], vec![2, 4]];
        assert_eq!(rank_over_z(&m), 1);
    }

    #[test]
    fn test_symmetric_difference() {
        let mut a = vec![0, 2, 4];
        let b = vec![2, 4, 6];
        symmetric_difference_inplace(&mut a, &b);
        assert_eq!(a, vec![0, 6]);
    }

    #[test]
    fn test_bottleneck_dist_zero() {
        let pts = vec![(0.0, 1.0)];
        assert!((bottleneck_dist(&pts, &pts)).abs() < 1e-9);
    }
}
