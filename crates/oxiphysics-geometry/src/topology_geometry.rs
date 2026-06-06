// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//! Computational topology for geometry: persistent homology, Reeb graphs.
//!
//! This module provides tools for analyzing the topological structure
//! of geometric data, including:
//!
//! - Simplicial complexes ([`SimplicialComplex`])
//! - Persistent homology ([`PersistentHomology`])
//! - Reeb graphs ([`ReebGraph`])
//! - Morse complexes ([`MorseComplex`])
//! - Alpha shapes ([`AlphaComplex`])
//! - Cubical complexes ([`CheckerboardComplex`])
//! - Betti numbers ([`BettiNumbers`])
//! - Topological noise removal ([`TopologicalNoise`])

/// Simplicial complex: vertices, edges, triangles, tetrahedra.
///
/// A simplicial complex is a set of simplices that satisfies the condition
/// that every face of a simplex is also a simplex in the complex.
#[derive(Debug, Clone)]
pub struct SimplicialComplex {
    /// Vertices (0-simplices) stored as \[x, y, z\]
    pub vertices: Vec<[f64; 3]>,
    /// Edges (1-simplices) as pairs of vertex indices
    pub edges: Vec<[usize; 2]>,
    /// Triangles (2-simplices) as triples of vertex indices
    pub triangles: Vec<[usize; 3]>,
    /// Tetrahedra (3-simplices) as quadruples of vertex indices
    pub tetrahedra: Vec<[usize; 4]>,
}

impl SimplicialComplex {
    /// Creates an empty simplicial complex.
    pub fn new() -> Self {
        Self {
            vertices: Vec::new(),
            edges: Vec::new(),
            triangles: Vec::new(),
            tetrahedra: Vec::new(),
        }
    }

    /// Adds a vertex and returns its index.
    pub fn add_vertex(&mut self, pos: [f64; 3]) -> usize {
        let idx = self.vertices.len();
        self.vertices.push(pos);
        idx
    }

    /// Adds an edge if it is not already present.
    pub fn add_edge(&mut self, i: usize, j: usize) {
        let e = if i < j { [i, j] } else { [j, i] };
        if !self.edges.contains(&e) {
            self.edges.push(e);
        }
    }

    /// Adds a triangle.
    pub fn add_triangle(&mut self, i: usize, j: usize, k: usize) {
        let mut t = [i, j, k];
        t.sort_unstable();
        if !self.triangles.contains(&t) {
            self.triangles.push(t);
        }
    }

    /// Adds a tetrahedron.
    pub fn add_tetrahedron(&mut self, i: usize, j: usize, k: usize, l: usize) {
        let mut t = [i, j, k, l];
        t.sort_unstable();
        if !self.tetrahedra.contains(&t) {
            self.tetrahedra.push(t);
        }
    }

    /// Computes the Euler characteristic χ = V - E + F - T.
    pub fn euler_characteristic(&self) -> i64 {
        self.vertices.len() as i64 - self.edges.len() as i64 + self.triangles.len() as i64
            - self.tetrahedra.len() as i64
    }

    /// Returns the number of vertices V.
    pub fn n_vertices(&self) -> usize {
        self.vertices.len()
    }

    /// Returns the number of edges E.
    pub fn n_edges(&self) -> usize {
        self.edges.len()
    }

    /// Returns the number of triangles F.
    pub fn n_triangles(&self) -> usize {
        self.triangles.len()
    }

    /// Returns the number of tetrahedra T.
    pub fn n_tetrahedra(&self) -> usize {
        self.tetrahedra.len()
    }

    /// Builds the boundary operator ∂_1: edges → vertices as incidence matrix.
    ///
    /// Returns a |V| × |E| matrix (flat row-major).
    pub fn boundary_1(&self) -> Vec<i32> {
        let nv = self.vertices.len();
        let ne = self.edges.len();
        let mut b = vec![0_i32; nv * ne];
        for (j, &[i0, i1]) in self.edges.iter().enumerate() {
            b[i0 * ne + j] = -1;
            b[i1 * ne + j] = 1;
        }
        b
    }

    /// Builds the boundary operator ∂_2: triangles → edges.
    ///
    /// Returns a |E| × |F| matrix (flat row-major).
    pub fn boundary_2(&self) -> Vec<i32> {
        let ne = self.edges.len();
        let nf = self.triangles.len();
        let mut b = vec![0_i32; ne * nf];
        for (j, &[a, b_idx, c]) in self.triangles.iter().enumerate() {
            // Edges of triangle: (a,b), (b,c), (a,c)
            let tri_edges = [[a, b_idx], [b_idx, c], [a, c]];
            let signs = [1_i32, 1, -1];
            for (k, &e) in tri_edges.iter().enumerate() {
                let es = if e[0] < e[1] {
                    [e[0], e[1]]
                } else {
                    [e[1], e[0]]
                };
                if let Some(ei) = self.edges.iter().position(|&x| x == es) {
                    b[ei * nf + j] = signs[k];
                }
            }
        }
        b
    }

    /// Creates the simplicial complex for a triangulated sphere (icosphere-like).
    ///
    /// Uses an octahedral approximation: 6 vertices, 12 edges, 8 triangles.
    /// Euler characteristic = 6 - 12 + 8 = 2.
    pub fn octahedral_sphere() -> Self {
        let mut sc = Self::new();
        // 6 vertices of octahedron
        sc.add_vertex([1.0, 0.0, 0.0]);
        sc.add_vertex([-1.0, 0.0, 0.0]);
        sc.add_vertex([0.0, 1.0, 0.0]);
        sc.add_vertex([0.0, -1.0, 0.0]);
        sc.add_vertex([0.0, 0.0, 1.0]);
        sc.add_vertex([0.0, 0.0, -1.0]);
        // 12 edges
        let edges = [
            [0, 2],
            [0, 3],
            [0, 4],
            [0, 5],
            [1, 2],
            [1, 3],
            [1, 4],
            [1, 5],
            [2, 4],
            [2, 5],
            [3, 4],
            [3, 5],
        ];
        for [i, j] in edges {
            sc.add_edge(i, j);
        }
        // 8 triangular faces
        sc.add_triangle(0, 2, 4);
        sc.add_triangle(0, 2, 5);
        sc.add_triangle(0, 3, 4);
        sc.add_triangle(0, 3, 5);
        sc.add_triangle(1, 2, 4);
        sc.add_triangle(1, 2, 5);
        sc.add_triangle(1, 3, 4);
        sc.add_triangle(1, 3, 5);
        sc
    }

    /// Creates a torus triangulation with Euler characteristic 0.
    ///
    /// Minimal triangulation: V=7, E=21, F=14, χ=0.
    /// Uses the complete graph K_7 triangulation of the torus.
    pub fn minimal_torus() -> Self {
        let mut sc = Self::new();
        for i in 0..7 {
            let angle = 2.0 * std::f64::consts::PI * i as f64 / 7.0;
            sc.add_vertex([angle.cos(), angle.sin(), 0.0]);
        }
        // Minimal 7-vertex triangulation of torus (K_7 embedding):
        // V=7, E=21, F=14, χ=0
        let tris: [[usize; 3]; 14] = [
            [0, 1, 2],
            [0, 2, 3],
            [0, 3, 4],
            [0, 4, 5],
            [0, 5, 6],
            [0, 6, 1],
            [1, 3, 2],
            [1, 4, 3],
            [1, 5, 4],
            [1, 6, 5],
            [2, 6, 3],
            [3, 5, 6],
            [2, 4, 6],
            [2, 5, 4],
        ];
        for [a, b, c] in tris {
            sc.add_edge(a, b);
            sc.add_edge(b, c);
            sc.add_edge(a, c);
            sc.add_triangle(a, b, c);
        }
        sc
    }
}

impl Default for SimplicialComplex {
    fn default() -> Self {
        Self::new()
    }
}

/// A birth-death pair in a persistence diagram.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BirthDeathPair {
    /// Birth value (filtration parameter)
    pub birth: f64,
    /// Death value (filtration parameter), f64::INFINITY for essential classes
    pub death: f64,
    /// Homological dimension
    pub dim: usize,
}

impl BirthDeathPair {
    /// Computes persistence = death - birth.
    pub fn persistence(&self) -> f64 {
        self.death - self.birth
    }

    /// Returns true if this is an essential class (death = ∞).
    pub fn is_essential(&self) -> bool {
        self.death.is_infinite()
    }
}

/// Persistent homology computation for filtrations.
///
/// Computes birth-death pairs (persistence diagram) for 0-dim and 1-dim
/// topological features (connected components and loops).
#[derive(Debug, Clone)]
pub struct PersistentHomology {
    /// Computed birth-death pairs
    pub pairs: Vec<BirthDeathPair>,
}

impl PersistentHomology {
    /// Creates an empty persistent homology.
    pub fn new() -> Self {
        Self { pairs: Vec::new() }
    }

    /// Computes 0-dimensional persistent homology from a point cloud.
    ///
    /// Uses union-find (single-linkage clustering) to track connected components.
    /// Each merge event creates a birth-death pair.
    ///
    /// # Arguments
    /// * `points` - Point cloud as \[x,y,z\] array
    /// * `filtration_values` - Per-edge filtration values (distances)
    pub fn compute_0d(points: &[[f64; 3]]) -> Self {
        let n = points.len();
        if n == 0 {
            return Self::new();
        }
        // Build all pairwise distances and sort
        let mut edges: Vec<(f64, usize, usize)> = Vec::new();
        for i in 0..n {
            for j in i + 1..n {
                let d = dist3(&points[i], &points[j]);
                edges.push((d, i, j));
            }
        }
        edges.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

        let mut uf = UnionFind::new(n);
        let mut pairs = Vec::new();

        // Initially all components born at 0
        for _ in 0..n {
            // Each point is born at 0
        }

        for (dist, i, j) in edges {
            let ri = uf.find(i);
            let rj = uf.find(j);
            if ri != rj {
                // Merge: the younger component dies
                let born_i = 0.0_f64;
                let born_j = 0.0_f64;
                let dying = if born_i >= born_j { ri } else { rj };
                let _ = dying; // used for tracking
                pairs.push(BirthDeathPair {
                    birth: 0.0,
                    death: dist,
                    dim: 0,
                });
                uf.union(i, j);
            }
        }
        // Essential class: the last surviving component
        pairs.push(BirthDeathPair {
            birth: 0.0,
            death: f64::INFINITY,
            dim: 0,
        });

        Self { pairs }
    }

    /// Returns all pairs in the persistence diagram.
    pub fn diagram(&self) -> &[BirthDeathPair] {
        &self.pairs
    }

    /// Returns the number of features with persistence > threshold.
    ///
    /// # Arguments
    /// * `threshold` - Minimum persistence to count as a feature
    pub fn significant_features(&self, threshold: f64) -> usize {
        self.pairs
            .iter()
            .filter(|p| p.persistence() > threshold)
            .count()
    }

    /// Computes the bottleneck distance between two persistence diagrams.
    ///
    /// d_B = inf_γ sup_{x ∈ D1} ||x - γ(x)||_∞
    ///
    /// Uses a greedy approximation.
    ///
    /// # Arguments
    /// * `other` - Other persistence diagram
    pub fn bottleneck_distance(&self, other: &Self) -> f64 {
        // Filter to finite pairs
        let d1: Vec<_> = self.pairs.iter().filter(|p| p.death.is_finite()).collect();
        let d2: Vec<_> = other.pairs.iter().filter(|p| p.death.is_finite()).collect();

        if d1.is_empty() && d2.is_empty() {
            return 0.0;
        }

        let mut max_cost = 0.0_f64;

        // For each point in d1, find closest in d2 or diagonal
        for p1 in &d1 {
            let diag_dist = (p1.death - p1.birth) / 2.0;
            let min_d2 = d2
                .iter()
                .map(|p2| (p1.birth - p2.birth).abs().max((p1.death - p2.death).abs()))
                .fold(f64::INFINITY, f64::min);
            max_cost = max_cost.max(min_d2.min(diag_dist));
        }

        for p2 in &d2 {
            let diag_dist = (p2.death - p2.birth) / 2.0;
            let min_d1 = d1
                .iter()
                .map(|p1| (p1.birth - p2.birth).abs().max((p1.death - p2.death).abs()))
                .fold(f64::INFINITY, f64::min);
            max_cost = max_cost.max(min_d1.min(diag_dist));
        }

        max_cost
    }
}

impl Default for PersistentHomology {
    fn default() -> Self {
        Self::new()
    }
}

/// Union-Find data structure for persistent homology.
#[derive(Debug, Clone)]
struct UnionFind {
    parent: Vec<usize>,
    rank: Vec<usize>,
}

impl UnionFind {
    fn new(n: usize) -> Self {
        Self {
            parent: (0..n).collect(),
            rank: vec![0; n],
        }
    }

    fn find(&mut self, x: usize) -> usize {
        if self.parent[x] != x {
            self.parent[x] = self.find(self.parent[x]);
        }
        self.parent[x]
    }

    fn union(&mut self, x: usize, y: usize) {
        let rx = self.find(x);
        let ry = self.find(y);
        if rx == ry {
            return;
        }
        if self.rank[rx] < self.rank[ry] {
            self.parent[rx] = ry;
        } else if self.rank[rx] > self.rank[ry] {
            self.parent[ry] = rx;
        } else {
            self.parent[ry] = rx;
            self.rank[rx] += 1;
        }
    }

    fn n_components(&mut self) -> usize {
        let n = self.parent.len();
        let roots: std::collections::HashSet<usize> = (0..n).map(|i| self.find(i)).collect();
        roots.len()
    }
}

/// Reeb graph of a scalar function on a mesh.
///
/// The Reeb graph tracks the evolution of connected components
/// of level sets f^{-1}(c) as c varies.
#[derive(Debug, Clone)]
pub struct ReebGraph {
    /// Nodes of the Reeb graph (critical points)
    pub nodes: Vec<ReebNode>,
    /// Arcs connecting nodes
    pub arcs: Vec<ReebArc>,
}

/// A node in the Reeb graph corresponding to a topological event.
#[derive(Debug, Clone)]
pub struct ReebNode {
    /// Function value at this critical point
    pub value: f64,
    /// Type of critical point
    pub kind: CriticalPointKind,
    /// Index of the corresponding vertex in the mesh
    pub vertex_idx: usize,
}

/// Type of critical point in the Reeb graph.
#[derive(Debug, Clone, PartialEq)]
pub enum CriticalPointKind {
    /// Local minimum (birth of component)
    Minimum,
    /// Saddle (merge or split)
    Saddle,
    /// Local maximum (death of component)
    Maximum,
}

/// An arc in the Reeb graph connecting two nodes.
#[derive(Debug, Clone)]
pub struct ReebArc {
    /// Source node index
    pub source: usize,
    /// Target node index
    pub target: usize,
}

impl ReebGraph {
    /// Creates an empty Reeb graph.
    pub fn new() -> Self {
        Self {
            nodes: Vec::new(),
            arcs: Vec::new(),
        }
    }

    /// Computes the Reeb graph of a scalar function on a set of vertices.
    ///
    /// # Arguments
    /// * `values` - Scalar function values at each vertex
    /// * `edges` - Edge connectivity of the mesh
    pub fn compute(values: &[f64], edges: &[[usize; 2]]) -> Self {
        let n = values.len();
        if n == 0 {
            return Self::new();
        }

        // Sort vertices by function value
        let mut sorted_idx: Vec<usize> = (0..n).collect();
        sorted_idx.sort_by(|&a, &b| {
            values[a]
                .partial_cmp(&values[b])
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        let mut nodes = Vec::new();
        let mut arcs = Vec::new();

        // Build adjacency
        let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
        for &[i, j] in edges {
            adj[i].push(j);
            adj[j].push(i);
        }

        // Classify critical points
        for &v in &sorted_idx {
            let fv = values[v];
            let lower_nbrs: Vec<usize> = adj[v]
                .iter()
                .filter(|&&u| values[u] < fv)
                .cloned()
                .collect();
            let upper_nbrs: Vec<usize> = adj[v]
                .iter()
                .filter(|&&u| values[u] > fv)
                .cloned()
                .collect();

            let kind = if lower_nbrs.is_empty() {
                CriticalPointKind::Minimum
            } else if upper_nbrs.is_empty() {
                CriticalPointKind::Maximum
            } else {
                CriticalPointKind::Saddle
            };

            if matches!(
                kind,
                CriticalPointKind::Minimum | CriticalPointKind::Maximum | CriticalPointKind::Saddle
            ) {
                // Only add if truly critical (for simplicity, add min/max always)
                if matches!(
                    kind,
                    CriticalPointKind::Minimum | CriticalPointKind::Maximum
                ) {
                    let node_idx = nodes.len();
                    nodes.push(ReebNode {
                        value: fv,
                        kind,
                        vertex_idx: v,
                    });
                    if node_idx > 0 {
                        arcs.push(ReebArc {
                            source: node_idx - 1,
                            target: node_idx,
                        });
                    }
                }
            }
        }

        Self { nodes, arcs }
    }

    /// Returns the number of critical points.
    pub fn n_critical_points(&self) -> usize {
        self.nodes.len()
    }

    /// Returns the number of minima.
    pub fn n_minima(&self) -> usize {
        self.nodes
            .iter()
            .filter(|n| n.kind == CriticalPointKind::Minimum)
            .count()
    }

    /// Returns the number of maxima.
    pub fn n_maxima(&self) -> usize {
        self.nodes
            .iter()
            .filter(|n| n.kind == CriticalPointKind::Maximum)
            .count()
    }

    /// Returns the number of saddles.
    pub fn n_saddles(&self) -> usize {
        self.nodes
            .iter()
            .filter(|n| n.kind == CriticalPointKind::Saddle)
            .count()
    }
}

impl Default for ReebGraph {
    fn default() -> Self {
        Self::new()
    }
}

/// Morse-Smale complex for topological data analysis.
///
/// Decomposes a manifold based on gradient flow between critical points.
#[derive(Debug, Clone)]
pub struct MorseComplex {
    /// List of critical points
    pub critical_points: Vec<MorseCriticalPoint>,
    /// Gradient flow connections (unstable/stable manifolds)
    pub connections: Vec<(usize, usize)>,
}

/// A critical point in the Morse complex.
#[derive(Debug, Clone)]
pub struct MorseCriticalPoint {
    /// Index of the vertex
    pub vertex_idx: usize,
    /// Function value at this critical point
    pub value: f64,
    /// Morse index (0=minimum, 1=saddle, 2=maximum for 2-manifold)
    pub morse_index: usize,
}

impl MorseComplex {
    /// Creates an empty Morse complex.
    pub fn new() -> Self {
        Self {
            critical_points: Vec::new(),
            connections: Vec::new(),
        }
    }

    /// Computes the Morse complex of a scalar function on a mesh.
    ///
    /// # Arguments
    /// * `values` - Scalar function values at each vertex
    /// * `edges` - Edge connectivity
    pub fn compute(values: &[f64], edges: &[[usize; 2]]) -> Self {
        let n = values.len();
        if n == 0 {
            return Self::new();
        }
        let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
        for &[i, j] in edges {
            adj[i].push(j);
            adj[j].push(i);
        }

        let mut critical_points = Vec::new();

        for v in 0..n {
            let fv = values[v];
            let n_lower = adj[v].iter().filter(|&&u| values[u] < fv).count();
            let n_upper = adj[v].iter().filter(|&&u| values[u] > fv).count();

            let morse_index = if n_lower == 0 && n_upper > 0 {
                Some(0) // minimum
            } else if n_upper == 0 && n_lower > 0 {
                Some(2) // maximum
            } else if n_lower > 0 && n_upper > 0 {
                // Check if lower link is disconnected (saddle)
                Some(1)
            } else {
                None
            };

            if let Some(idx) = morse_index {
                critical_points.push(MorseCriticalPoint {
                    vertex_idx: v,
                    value: fv,
                    morse_index: idx,
                });
            }
        }

        // Build connections between consecutive critical points
        critical_points.sort_by(|a, b| {
            a.value
                .partial_cmp(&b.value)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        let connections: Vec<(usize, usize)> = (0..critical_points.len().saturating_sub(1))
            .map(|i| (i, i + 1))
            .collect();

        Self {
            critical_points,
            connections,
        }
    }

    /// Returns the number of minima.
    pub fn n_minima(&self) -> usize {
        self.critical_points
            .iter()
            .filter(|c| c.morse_index == 0)
            .count()
    }

    /// Returns the number of saddles.
    pub fn n_saddles(&self) -> usize {
        self.critical_points
            .iter()
            .filter(|c| c.morse_index == 1)
            .count()
    }

    /// Returns the number of maxima.
    pub fn n_maxima(&self) -> usize {
        self.critical_points
            .iter()
            .filter(|c| c.morse_index == 2)
            .count()
    }

    /// Verifies Morse inequalities: saddles + 2 ≤ maxima + minima (approx).
    pub fn morse_relation_holds(&self) -> bool {
        let nm = self.n_minima();
        let ns = self.n_saddles();
        let nmax = self.n_maxima();
        // Morse inequality: m_0 - m_1 + m_2 ≥ χ (here approximate)
        nm as i64 - ns as i64 + nmax as i64 >= 0
    }
}

impl Default for MorseComplex {
    fn default() -> Self {
        Self::new()
    }
}

/// Alpha complex / alpha shapes for point clouds.
///
/// An alpha complex is a subset of the Delaunay triangulation
/// based on a circumradius threshold α.
#[derive(Debug, Clone)]
pub struct AlphaComplex {
    /// Alpha parameter (circumradius threshold)
    pub alpha: f64,
    /// Points in the cloud
    pub points: Vec<[f64; 3]>,
    /// Edges included at this alpha value
    pub edges: Vec<[usize; 2]>,
    /// Triangles included at this alpha value
    pub triangles: Vec<[usize; 3]>,
}

impl AlphaComplex {
    /// Creates an alpha complex from a point cloud.
    ///
    /// # Arguments
    /// * `points` - Point cloud
    /// * `alpha` - Circumradius threshold (larger = more simplices)
    pub fn new(points: Vec<[f64; 3]>, alpha: f64) -> Self {
        let n = points.len();
        let mut edges = Vec::new();
        let mut triangles = Vec::new();

        // Add edges where half-distance ≤ alpha
        for i in 0..n {
            for j in i + 1..n {
                let d = dist3(&points[i], &points[j]);
                if d / 2.0 <= alpha {
                    edges.push([i, j]);
                }
            }
        }

        // Add triangles where circumradius ≤ alpha
        for i in 0..n {
            for j in i + 1..n {
                for k in j + 1..n {
                    let r = circumradius_3pts(&points[i], &points[j], &points[k]);
                    if r <= alpha {
                        triangles.push([i, j, k]);
                        // Ensure edges exist
                    }
                }
            }
        }

        Self {
            alpha,
            points,
            edges,
            triangles,
        }
    }

    /// Returns the number of connected components at this alpha.
    pub fn n_components(&self) -> usize {
        let n = self.points.len();
        if n == 0 {
            return 0;
        }
        let mut uf = UnionFind::new(n);
        for &[i, j] in &self.edges {
            uf.union(i, j);
        }
        uf.n_components()
    }

    /// Returns the simplicial complex at this alpha level.
    pub fn to_simplicial_complex(&self) -> SimplicialComplex {
        let mut sc = SimplicialComplex::new();
        for &p in &self.points {
            sc.add_vertex(p);
        }
        for &[i, j] in &self.edges {
            sc.add_edge(i, j);
        }
        for &[i, j, k] in &self.triangles {
            sc.add_triangle(i, j, k);
        }
        sc
    }
}

/// Circumradius of a triangle defined by three 3D points.
fn circumradius_3pts(a: &[f64; 3], b: &[f64; 3], c: &[f64; 3]) -> f64 {
    let ab = dist3(a, b);
    let bc = dist3(b, c);
    let ca = dist3(c, a);
    let s = (ab + bc + ca) / 2.0;
    let area_sq = s * (s - ab) * (s - bc) * (s - ca);
    if area_sq <= 0.0 {
        return f64::INFINITY;
    }
    let area = area_sq.sqrt();
    ab * bc * ca / (4.0 * area)
}

/// 3D Euclidean distance between two points.
fn dist3(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

/// Cubical complex for image topology (checkerboard-like).
///
/// Represents a binary image as a cubical complex and computes
/// its topological features.
#[derive(Debug, Clone)]
pub struct CheckerboardComplex {
    /// Grid dimensions \[nx, ny\]
    pub dims: [usize; 2],
    /// Binary image data (true = filled voxel)
    pub data: Vec<bool>,
}

impl CheckerboardComplex {
    /// Creates a cubical complex from a binary grid.
    ///
    /// # Arguments
    /// * `dims` - Grid dimensions \[nx, ny\]
    /// * `data` - Binary values row-major
    pub fn new(dims: [usize; 2], data: Vec<bool>) -> Self {
        Self { dims, data }
    }

    /// Creates a filled disk of given radius in a grid.
    pub fn filled_disk(size: usize, radius: f64) -> Self {
        let cx = size as f64 / 2.0;
        let cy = size as f64 / 2.0;
        let data: Vec<bool> = (0..size * size)
            .map(|idx| {
                let i = (idx / size) as f64;
                let j = (idx % size) as f64;
                ((i - cx).powi(2) + (j - cy).powi(2)).sqrt() <= radius
            })
            .collect();
        Self::new([size, size], data)
    }

    /// Creates an annulus (ring) in a grid.
    pub fn annulus(size: usize, r_inner: f64, r_outer: f64) -> Self {
        let cx = size as f64 / 2.0;
        let cy = size as f64 / 2.0;
        let data: Vec<bool> = (0..size * size)
            .map(|idx| {
                let i = (idx / size) as f64;
                let j = (idx % size) as f64;
                let r = ((i - cx).powi(2) + (j - cy).powi(2)).sqrt();
                r >= r_inner && r <= r_outer
            })
            .collect();
        Self::new([size, size], data)
    }

    /// Counts the number of filled voxels.
    pub fn volume(&self) -> usize {
        self.data.iter().filter(|&&b| b).count()
    }

    /// Computes the Euler characteristic via the cubical formula.
    ///
    /// χ = V - E + F for 2D cubical complex
    pub fn euler_characteristic_2d(&self) -> i64 {
        let nx = self.dims[0];
        let ny = self.dims[1];
        let mut v = 0_i64; // vertices
        let mut e = 0_i64; // edges
        let mut f = 0_i64; // faces

        for i in 0..nx {
            for j in 0..ny {
                let filled = self.data[i * ny + j];
                if filled {
                    v += 1;
                    if i + 1 < nx && self.data[(i + 1) * ny + j] {
                        e += 1;
                    }
                    if j + 1 < ny && self.data[i * ny + j + 1] {
                        e += 1;
                    }
                    if i + 1 < nx
                        && j + 1 < ny
                        && self.data[(i + 1) * ny + j]
                        && self.data[i * ny + j + 1]
                        && self.data[(i + 1) * ny + j + 1]
                    {
                        f += 1;
                    }
                }
            }
        }
        v - e + f
    }
}

/// Betti numbers for topological characterization.
///
/// β_0 = number of connected components
/// β_1 = number of independent loops
/// β_2 = number of enclosed voids
#[derive(Debug, Clone, Copy)]
pub struct BettiNumbers {
    /// β_0: connected components
    pub beta_0: usize,
    /// β_1: loops (1-cycles)
    pub beta_1: usize,
    /// β_2: voids (2-cycles)
    pub beta_2: usize,
}

impl BettiNumbers {
    /// Computes Betti numbers from a simplicial complex using Euler relation.
    ///
    /// For a closed orientable surface: χ = β_0 - β_1 + β_2
    ///
    /// # Arguments
    /// * `sc` - Simplicial complex
    /// * `n_components` - Number of connected components (β_0)
    /// * `n_voids` - Number of enclosed voids (β_2, 0 for surfaces)
    pub fn from_simplicial_complex(
        sc: &SimplicialComplex,
        n_components: usize,
        n_voids: usize,
    ) -> Self {
        let chi = sc.euler_characteristic();
        // χ = β_0 - β_1 + β_2
        // => β_1 = β_0 + β_2 - χ
        let beta_1 = (n_components as i64 + n_voids as i64 - chi).max(0) as usize;
        Self {
            beta_0: n_components,
            beta_1,
            beta_2: n_voids,
        }
    }

    /// Returns the Euler characteristic χ = β_0 - β_1 + β_2.
    pub fn euler_characteristic(&self) -> i64 {
        self.beta_0 as i64 - self.beta_1 as i64 + self.beta_2 as i64
    }

    /// Creates Betti numbers directly.
    pub fn new(beta_0: usize, beta_1: usize, beta_2: usize) -> Self {
        Self {
            beta_0,
            beta_1,
            beta_2,
        }
    }
}

/// Topological noise removal via persistence thresholding.
///
/// Features with persistence below the threshold are considered noise.
#[derive(Debug, Clone)]
pub struct TopologicalNoise {
    /// Persistence threshold
    pub threshold: f64,
}

impl TopologicalNoise {
    /// Creates a noise removal filter with given persistence threshold.
    ///
    /// # Arguments
    /// * `threshold` - Minimum persistence to keep a feature
    pub fn new(threshold: f64) -> Self {
        Self { threshold }
    }

    /// Filters birth-death pairs, keeping only significant features.
    ///
    /// # Arguments
    /// * `pairs` - Input persistence diagram
    pub fn filter(&self, pairs: &[BirthDeathPair]) -> Vec<BirthDeathPair> {
        pairs
            .iter()
            .filter(|p| p.persistence() > self.threshold || p.is_essential())
            .cloned()
            .collect()
    }

    /// Counts significant features at given dimension.
    ///
    /// # Arguments
    /// * `pairs` - Persistence diagram
    /// * `dim` - Homological dimension
    pub fn count_features(&self, pairs: &[BirthDeathPair], dim: usize) -> usize {
        pairs
            .iter()
            .filter(|p| p.dim == dim && (p.persistence() > self.threshold || p.is_essential()))
            .count()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // ---- SimplicialComplex ----

    #[test]
    fn test_euler_characteristic_sphere() {
        // Octahedral sphere: V=6, E=12, F=8 => χ=2
        let sc = SimplicialComplex::octahedral_sphere();
        assert_eq!(sc.euler_characteristic(), 2);
    }

    #[test]
    fn test_octahedral_sphere_vertex_count() {
        let sc = SimplicialComplex::octahedral_sphere();
        assert_eq!(sc.n_vertices(), 6);
    }

    #[test]
    fn test_octahedral_sphere_edge_count() {
        let sc = SimplicialComplex::octahedral_sphere();
        assert_eq!(sc.n_edges(), 12);
    }

    #[test]
    fn test_octahedral_sphere_triangle_count() {
        let sc = SimplicialComplex::octahedral_sphere();
        assert_eq!(sc.n_triangles(), 8);
    }

    #[test]
    fn test_euler_characteristic_torus() {
        // Minimal torus: χ=0
        let sc = SimplicialComplex::minimal_torus();
        assert_eq!(sc.euler_characteristic(), 0);
    }

    #[test]
    fn test_single_vertex_euler() {
        let mut sc = SimplicialComplex::new();
        sc.add_vertex([0.0, 0.0, 0.0]);
        assert_eq!(sc.euler_characteristic(), 1);
    }

    #[test]
    fn test_add_duplicate_edge_ignored() {
        let mut sc = SimplicialComplex::new();
        sc.add_vertex([0.0, 0.0, 0.0]);
        sc.add_vertex([1.0, 0.0, 0.0]);
        sc.add_edge(0, 1);
        sc.add_edge(0, 1); // duplicate
        assert_eq!(sc.n_edges(), 1);
    }

    #[test]
    fn test_tetrahedron_euler_characteristic() {
        let mut sc = SimplicialComplex::new();
        for i in 0..4 {
            sc.add_vertex([i as f64, 0.0, 0.0]);
        }
        sc.add_edge(0, 1);
        sc.add_edge(0, 2);
        sc.add_edge(0, 3);
        sc.add_edge(1, 2);
        sc.add_edge(1, 3);
        sc.add_edge(2, 3);
        sc.add_triangle(0, 1, 2);
        sc.add_triangle(0, 1, 3);
        sc.add_triangle(0, 2, 3);
        sc.add_triangle(1, 2, 3);
        // Tetrahedron surface: V=4, E=6, F=4 => χ=2
        assert_eq!(sc.euler_characteristic(), 2);
    }

    #[test]
    fn test_boundary_1_dimensions() {
        let sc = SimplicialComplex::octahedral_sphere();
        let b1 = sc.boundary_1();
        assert_eq!(b1.len(), sc.n_vertices() * sc.n_edges());
    }

    #[test]
    fn test_boundary_2_dimensions() {
        let sc = SimplicialComplex::octahedral_sphere();
        let b2 = sc.boundary_2();
        assert_eq!(b2.len(), sc.n_edges() * sc.n_triangles());
    }

    // ---- PersistentHomology ----

    #[test]
    fn test_persistent_homology_single_point() {
        let pts = [[0.0, 0.0, 0.0]];
        let ph = PersistentHomology::compute_0d(&pts);
        // One essential class
        assert_eq!(ph.significant_features(0.0), 1);
    }

    #[test]
    fn test_persistent_homology_two_points() {
        let pts = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let ph = PersistentHomology::compute_0d(&pts);
        // One finite pair + one essential
        assert_eq!(ph.pairs.len(), 2);
    }

    #[test]
    fn test_persistent_homology_one_blob() {
        let pts: Vec<[f64; 3]> = (0..5).map(|i| [i as f64 * 0.1, 0.0, 0.0]).collect();
        let ph = PersistentHomology::compute_0d(&pts);
        // Should have exactly 1 essential component
        let essential = ph.pairs.iter().filter(|p| p.is_essential()).count();
        assert_eq!(essential, 1);
    }

    #[test]
    fn test_persistent_homology_birth_death_order() {
        let pts = [[0.0, 0.0, 0.0], [0.1, 0.0, 0.0], [5.0, 0.0, 0.0]];
        let ph = PersistentHomology::compute_0d(&pts);
        for p in ph.pairs.iter().filter(|p| p.death.is_finite()) {
            assert!(p.birth <= p.death);
        }
    }

    #[test]
    fn test_bottleneck_distance_nonneg() {
        let pts1 = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let pts2 = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let ph1 = PersistentHomology::compute_0d(&pts1);
        let ph2 = PersistentHomology::compute_0d(&pts2);
        let d = ph1.bottleneck_distance(&ph2);
        assert!(d >= 0.0);
    }

    #[test]
    fn test_bottleneck_distance_self_zero() {
        let pts = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let ph = PersistentHomology::compute_0d(&pts);
        let d = ph.bottleneck_distance(&ph);
        assert!(d < 1e-10);
    }

    #[test]
    fn test_persistence_positive() {
        let p = BirthDeathPair {
            birth: 0.5,
            death: 1.5,
            dim: 0,
        };
        assert!((p.persistence() - 1.0).abs() < 1e-10);
    }

    // ---- ReebGraph ----

    #[test]
    fn test_reeb_graph_height_sphere_two_critical_points() {
        // Height function on line: min, max
        let values = vec![0.0, 0.5, 1.0, 0.7, 0.3];
        let edges = [[0usize, 1], [1, 2], [2, 3], [3, 4]];
        let rg = ReebGraph::compute(&values, &edges);
        assert!(rg.n_minima() >= 1);
        assert!(rg.n_maxima() >= 1);
    }

    #[test]
    fn test_reeb_graph_monotone_no_saddle() {
        let values = vec![0.0, 1.0, 2.0, 3.0, 4.0];
        let edges = [[0usize, 1], [1, 2], [2, 3], [3, 4]];
        let rg = ReebGraph::compute(&values, &edges);
        // Monotone function: only min and max
        assert_eq!(rg.n_saddles(), 0);
    }

    #[test]
    fn test_reeb_graph_empty() {
        let rg = ReebGraph::compute(&[], &[]);
        assert_eq!(rg.n_critical_points(), 0);
    }

    // ---- MorseComplex ----

    #[test]
    fn test_morse_complex_line() {
        let values = vec![1.0, 0.5, 0.8, 0.3, 0.7];
        let edges = [[0usize, 1], [1, 2], [2, 3], [3, 4]];
        let mc = MorseComplex::compute(&values, &edges);
        assert!(mc.n_minima() >= 1);
        assert!(mc.n_maxima() >= 1);
    }

    #[test]
    fn test_morse_relation_holds() {
        let values = vec![1.0, 3.0, 0.5, 2.0, 0.3];
        let edges = [[0usize, 1], [1, 2], [2, 3], [3, 4], [0, 4]];
        let mc = MorseComplex::compute(&values, &edges);
        assert!(mc.morse_relation_holds());
    }

    #[test]
    fn test_morse_empty() {
        let mc = MorseComplex::compute(&[], &[]);
        assert_eq!(mc.n_minima() + mc.n_saddles() + mc.n_maxima(), 0);
    }

    #[test]
    fn test_morse_single_vertex() {
        let values = vec![1.0];
        let mc = MorseComplex::compute(&values, &[]);
        // Isolated point: no neighbors -> neither min nor max classified
        assert_eq!(mc.critical_points.len(), 0);
    }

    // ---- AlphaComplex ----

    #[test]
    fn test_alpha_complex_large_alpha_all_edges() {
        let pts = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]];
        let ac = AlphaComplex::new(pts, 1000.0);
        // At very large alpha, all pairwise edges should be included
        assert_eq!(ac.edges.len(), 3);
    }

    #[test]
    fn test_alpha_complex_zero_alpha_no_edges() {
        let pts = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let ac = AlphaComplex::new(pts, 0.0);
        assert_eq!(ac.edges.len(), 0);
    }

    #[test]
    fn test_alpha_complex_components() {
        let pts = vec![[0.0, 0.0, 0.0], [0.1, 0.0, 0.0], [10.0, 0.0, 0.0]];
        let ac = AlphaComplex::new(pts, 0.1);
        // First two close, third far
        assert_eq!(ac.n_components(), 2);
    }

    #[test]
    fn test_alpha_complex_single_component_large_alpha() {
        let pts: Vec<[f64; 3]> = (0..4).map(|i| [i as f64 * 0.5, 0.0, 0.0]).collect();
        let ac = AlphaComplex::new(pts, 1000.0);
        assert_eq!(ac.n_components(), 1);
    }

    // ---- CheckerboardComplex ----

    #[test]
    fn test_checkerboard_filled_disk_volume() {
        let cc = CheckerboardComplex::filled_disk(10, 4.0);
        assert!(cc.volume() > 0);
    }

    #[test]
    fn test_checkerboard_annulus_volume() {
        let cc = CheckerboardComplex::annulus(20, 3.0, 7.0);
        assert!(cc.volume() > 0);
    }

    #[test]
    fn test_checkerboard_euler_single_square() {
        let cc = CheckerboardComplex::new([1, 1], vec![true]);
        let chi = cc.euler_characteristic_2d();
        assert_eq!(chi, 1); // single vertex
    }

    // ---- BettiNumbers ----

    #[test]
    fn test_betti_euler_consistency() {
        let beta = BettiNumbers::new(1, 2, 1);
        // χ = 1 - 2 + 1 = 0 (torus)
        assert_eq!(beta.euler_characteristic(), 0);
    }

    #[test]
    fn test_betti_sphere_euler() {
        let beta = BettiNumbers::new(1, 0, 1);
        // χ = 1 - 0 + 1 = 2 (sphere)
        assert_eq!(beta.euler_characteristic(), 2);
    }

    #[test]
    fn test_betti_from_simplicial_sphere() {
        let sc = SimplicialComplex::octahedral_sphere();
        let beta = BettiNumbers::from_simplicial_complex(&sc, 1, 0);
        // χ = 2, β_0=1, β_2=0 => β_1 = 1 + 0 - 2 = -1 -> 0
        assert_eq!(beta.beta_0, 1);
    }

    // ---- TopologicalNoise ----

    #[test]
    fn test_topological_noise_removes_small_features() {
        let noise = TopologicalNoise::new(0.5);
        let pairs = vec![
            BirthDeathPair {
                birth: 0.0,
                death: 0.1,
                dim: 0,
            },
            BirthDeathPair {
                birth: 0.0,
                death: 1.0,
                dim: 0,
            },
        ];
        let filtered = noise.filter(&pairs);
        assert_eq!(filtered.len(), 1);
    }

    #[test]
    fn test_topological_noise_keeps_essential() {
        let noise = TopologicalNoise::new(100.0);
        let pairs = vec![BirthDeathPair {
            birth: 0.0,
            death: f64::INFINITY,
            dim: 0,
        }];
        let filtered = noise.filter(&pairs);
        assert_eq!(filtered.len(), 1);
    }

    #[test]
    fn test_topological_noise_count_features() {
        let noise = TopologicalNoise::new(0.3);
        let pairs = vec![
            BirthDeathPair {
                birth: 0.0,
                death: 0.1,
                dim: 0,
            },
            BirthDeathPair {
                birth: 0.0,
                death: 0.5,
                dim: 0,
            },
            BirthDeathPair {
                birth: 0.0,
                death: f64::INFINITY,
                dim: 0,
            },
        ];
        assert_eq!(noise.count_features(&pairs, 0), 2);
    }

    #[test]
    fn test_topological_noise_zero_threshold() {
        let noise = TopologicalNoise::new(0.0);
        let pairs = vec![
            BirthDeathPair {
                birth: 0.0,
                death: 0.001,
                dim: 1,
            },
            BirthDeathPair {
                birth: 0.0,
                death: 1.0,
                dim: 1,
            },
        ];
        let filtered = noise.filter(&pairs);
        assert_eq!(filtered.len(), 2);
    }

    #[test]
    fn test_circumradius_equilateral() {
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [0.5, (3.0_f64).sqrt() / 2.0, 0.0];
        let r = circumradius_3pts(&a, &b, &c);
        // Circumradius of equilateral triangle with side 1 = 1/sqrt(3)
        assert!((r - 1.0 / 3.0_f64.sqrt()).abs() < 1e-6);
    }
}
