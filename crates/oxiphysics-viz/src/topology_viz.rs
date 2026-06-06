// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Topology visualization primitives.
//!
//! Provides data structures and algorithms for:
//!
//! - **Persistence diagrams** — birth/death pairs with bottleneck distance
//! - **Persistence barcodes** — horizontal bar intervals sorted by various criteria
//! - **Reeb graphs** — nodes at critical values, edges as level-set components
//! - **Morse-Smale complexes** — critical points, ascending/descending manifolds
//! - **Contour trees** — merge of join tree and split tree
//! - **Topological features** — persistence, dimension, representative cycles
//! - **Simplification** — pair cancellation, Reeb graph pruning
//! - **Visualization styles** — color by dimension, size by persistence

use std::collections::{BTreeMap, HashMap, HashSet, VecDeque};

// ---------------------------------------------------------------------------
// TopologicalFeature
// ---------------------------------------------------------------------------

/// A single topological feature extracted from a filtration.
#[derive(Debug, Clone, PartialEq)]
pub struct TopologicalFeature {
    /// Homology dimension (0 = connected component, 1 = loop, 2 = void, …).
    pub dimension: usize,
    /// Filtration value at which the feature is born.
    pub birth: f64,
    /// Filtration value at which the feature dies (`f64::INFINITY` if essential).
    pub death: f64,
    /// Optional representative cycle as a list of simplex indices.
    pub representative_cycle: Vec<usize>,
}

impl TopologicalFeature {
    /// Create a new topological feature.
    pub fn new(dimension: usize, birth: f64, death: f64) -> Self {
        Self {
            dimension,
            birth,
            death,
            representative_cycle: Vec::new(),
        }
    }

    /// Create a feature with a representative cycle.
    pub fn with_cycle(dimension: usize, birth: f64, death: f64, cycle: Vec<usize>) -> Self {
        Self {
            dimension,
            birth,
            death,
            representative_cycle: cycle,
        }
    }

    /// Persistence (death − birth).
    pub fn persistence(&self) -> f64 {
        self.death - self.birth
    }

    /// Midlife value (birth + death) / 2.
    pub fn midlife(&self) -> f64 {
        (self.birth + self.death) / 2.0
    }

    /// Whether this feature is essential (never dies).
    pub fn is_essential(&self) -> bool {
        self.death.is_infinite()
    }
}

// ---------------------------------------------------------------------------
// PersistenceDiagram
// ---------------------------------------------------------------------------

/// A persistence diagram: a multiset of birth/death pairs.
#[derive(Debug, Clone)]
pub struct PersistenceDiagram {
    /// Features (birth/death pairs).
    pub features: Vec<TopologicalFeature>,
}

impl PersistenceDiagram {
    /// Create an empty diagram.
    pub fn new() -> Self {
        Self {
            features: Vec::new(),
        }
    }

    /// Create from a list of (dimension, birth, death) tuples.
    pub fn from_tuples(tuples: &[(usize, f64, f64)]) -> Self {
        let features = tuples
            .iter()
            .map(|&(dim, b, d)| TopologicalFeature::new(dim, b, d))
            .collect();
        Self { features }
    }

    /// Add a feature.
    pub fn add_feature(&mut self, feature: TopologicalFeature) {
        self.features.push(feature);
    }

    /// Add a birth/death pair for a given dimension.
    pub fn add_pair(&mut self, dimension: usize, birth: f64, death: f64) {
        self.features
            .push(TopologicalFeature::new(dimension, birth, death));
    }

    /// Number of features.
    pub fn len(&self) -> usize {
        self.features.len()
    }

    /// Whether the diagram is empty.
    pub fn is_empty(&self) -> bool {
        self.features.is_empty()
    }

    /// Features filtered by homology dimension.
    pub fn dimension(&self, dim: usize) -> Vec<&TopologicalFeature> {
        self.features
            .iter()
            .filter(|f| f.dimension == dim)
            .collect()
    }

    /// Maximum persistence across all features.
    pub fn max_persistence(&self) -> f64 {
        self.features
            .iter()
            .filter(|f| !f.is_essential())
            .map(|f| f.persistence())
            .fold(0.0_f64, f64::max)
    }

    /// Features with persistence above a threshold.
    pub fn significant_features(&self, threshold: f64) -> Vec<&TopologicalFeature> {
        self.features
            .iter()
            .filter(|f| f.persistence() > threshold)
            .collect()
    }

    /// Bottleneck distance to another persistence diagram.
    ///
    /// Uses a greedy approximation: for each point in the larger diagram,
    /// find the nearest unmatched point in the other, or match to the diagonal.
    pub fn bottleneck_distance(&self, other: &PersistenceDiagram) -> f64 {
        let pts_a = self.non_essential_points();
        let pts_b = other.non_essential_points();

        if pts_a.is_empty() && pts_b.is_empty() {
            return 0.0;
        }

        // Cost of matching point to diagonal: persistence / 2
        let diag_cost = |p: &(f64, f64)| (p.1 - p.0).abs() / 2.0;

        // L∞ distance between two points
        let linf = |a: &(f64, f64), b: &(f64, f64)| (a.0 - b.0).abs().max((a.1 - b.1).abs());

        // Greedy matching
        let mut used_b = vec![false; pts_b.len()];
        let mut max_cost = 0.0_f64;

        for a in &pts_a {
            let mut best_cost = diag_cost(a);
            let mut best_j: Option<usize> = None;

            for (j, b) in pts_b.iter().enumerate() {
                if used_b[j] {
                    continue;
                }
                let cost = linf(a, b);
                if cost <= best_cost {
                    best_cost = cost;
                    best_j = Some(j);
                }
            }

            if let Some(j) = best_j {
                used_b[j] = true;
            }
            max_cost = max_cost.max(best_cost);
        }

        // Unmatched points in B go to diagonal
        for (j, b) in pts_b.iter().enumerate() {
            if !used_b[j] {
                max_cost = max_cost.max(diag_cost(b));
            }
        }

        max_cost
    }

    /// Extract (birth, death) pairs for non-essential features.
    fn non_essential_points(&self) -> Vec<(f64, f64)> {
        self.features
            .iter()
            .filter(|f| !f.is_essential())
            .map(|f| (f.birth, f.death))
            .collect()
    }

    /// Wasserstein-1 distance approximation to another diagram.
    pub fn wasserstein_distance(&self, other: &PersistenceDiagram) -> f64 {
        let pts_a = self.non_essential_points();
        let pts_b = other.non_essential_points();

        let diag_cost = |p: &(f64, f64)| (p.1 - p.0).abs() / 2.0;
        let linf = |a: &(f64, f64), b: &(f64, f64)| (a.0 - b.0).abs().max((a.1 - b.1).abs());

        let mut used_b = vec![false; pts_b.len()];
        let mut total = 0.0_f64;

        for a in &pts_a {
            let mut best_cost = diag_cost(a);
            let mut best_j: Option<usize> = None;
            for (j, b) in pts_b.iter().enumerate() {
                if used_b[j] {
                    continue;
                }
                let cost = linf(a, b);
                if cost < best_cost {
                    best_cost = cost;
                    best_j = Some(j);
                }
            }
            if let Some(j) = best_j {
                used_b[j] = true;
            }
            total += best_cost;
        }
        for (j, b) in pts_b.iter().enumerate() {
            if !used_b[j] {
                total += diag_cost(b);
            }
        }
        total
    }

    /// Betti numbers at a given filtration value.
    ///
    /// Returns a map from dimension to the count of features alive at `t`.
    pub fn betti_numbers(&self, t: f64) -> HashMap<usize, usize> {
        let mut betti: HashMap<usize, usize> = HashMap::new();
        for f in &self.features {
            if f.birth <= t && (f.death > t || f.is_essential()) {
                *betti.entry(f.dimension).or_insert(0) += 1;
            }
        }
        betti
    }
}

impl Default for PersistenceDiagram {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// PersistenceBarcode
// ---------------------------------------------------------------------------

/// Sort criterion for barcode intervals.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BarcodeSortOrder {
    /// Sort by birth time (ascending).
    Birth,
    /// Sort by death time (ascending).
    Death,
    /// Sort by persistence (descending).
    Persistence,
    /// Sort by dimension then persistence.
    DimensionThenPersistence,
}

/// A persistence barcode: horizontal-bar representation of persistence intervals.
#[derive(Debug, Clone)]
pub struct PersistenceBarcode {
    /// Underlying diagram.
    pub diagram: PersistenceDiagram,
    /// Sort order for display.
    pub sort_order: BarcodeSortOrder,
    /// Maximum filtration value for display (caps infinite deaths).
    pub max_filtration: f64,
}

impl PersistenceBarcode {
    /// Create a barcode from a persistence diagram.
    pub fn new(diagram: PersistenceDiagram, max_filtration: f64) -> Self {
        Self {
            diagram,
            sort_order: BarcodeSortOrder::Birth,
            max_filtration,
        }
    }

    /// Set the sort order.
    pub fn with_sort_order(mut self, order: BarcodeSortOrder) -> Self {
        self.sort_order = order;
        self
    }

    /// Return intervals as (birth, death, dimension) tuples, sorted according
    /// to the current sort order.
    pub fn intervals(&self) -> Vec<(f64, f64, usize)> {
        let mut intervals: Vec<(f64, f64, usize)> = self
            .diagram
            .features
            .iter()
            .map(|f| {
                let death = if f.is_essential() {
                    self.max_filtration
                } else {
                    f.death
                };
                (f.birth, death, f.dimension)
            })
            .collect();

        match self.sort_order {
            BarcodeSortOrder::Birth => {
                intervals
                    .sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
            }
            BarcodeSortOrder::Death => {
                intervals
                    .sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
            }
            BarcodeSortOrder::Persistence => {
                intervals.sort_by(|a, b| {
                    (b.1 - b.0)
                        .partial_cmp(&(a.1 - a.0))
                        .unwrap_or(std::cmp::Ordering::Equal)
                });
            }
            BarcodeSortOrder::DimensionThenPersistence => {
                intervals.sort_by(|a, b| {
                    a.2.cmp(&b.2).then_with(|| {
                        (b.1 - b.0)
                            .partial_cmp(&(a.1 - a.0))
                            .unwrap_or(std::cmp::Ordering::Equal)
                    })
                });
            }
        }
        intervals
    }

    /// Number of intervals.
    pub fn interval_count(&self) -> usize {
        self.diagram.len()
    }

    /// Intervals for a specific dimension.
    pub fn intervals_for_dimension(&self, dim: usize) -> Vec<(f64, f64)> {
        self.diagram
            .features
            .iter()
            .filter(|f| f.dimension == dim)
            .map(|f| {
                let death = if f.is_essential() {
                    self.max_filtration
                } else {
                    f.death
                };
                (f.birth, death)
            })
            .collect()
    }
}

// ---------------------------------------------------------------------------
// ReebGraph
// ---------------------------------------------------------------------------

/// A node in a Reeb graph, representing a critical value of the scalar function.
#[derive(Debug, Clone)]
pub struct ReebNode {
    /// Unique node identifier.
    pub id: usize,
    /// Scalar (function) value at this critical point.
    pub value: f64,
    /// Type of critical point.
    pub critical_type: CriticalPointType,
    /// Position in ambient space (if available): (x, y, z).
    pub position: Option<[f64; 3]>,
}

/// Type of critical point in a scalar field.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CriticalPointType {
    /// Local minimum.
    Minimum,
    /// Saddle point (index-1 in 2D, various indices in 3D).
    Saddle,
    /// Local maximum.
    Maximum,
    /// Regular (non-critical) point — used in augmented structures.
    Regular,
}

/// An edge in a Reeb graph connecting two critical nodes.
#[derive(Debug, Clone)]
pub struct ReebEdge {
    /// Source node ID.
    pub source: usize,
    /// Target node ID.
    pub target: usize,
    /// Number of connected components of the level set on this arc.
    pub component_count: usize,
    /// Persistence of this edge (difference of endpoint values).
    pub persistence: f64,
}

/// A Reeb graph capturing the topology of level sets of a scalar function.
#[derive(Debug, Clone)]
pub struct ReebGraph {
    /// Nodes (critical points).
    pub nodes: Vec<ReebNode>,
    /// Edges (arcs between critical points).
    pub edges: Vec<ReebEdge>,
    /// Adjacency list: node_id → list of edge indices.
    adjacency: HashMap<usize, Vec<usize>>,
}

impl ReebGraph {
    /// Create an empty Reeb graph.
    pub fn new() -> Self {
        Self {
            nodes: Vec::new(),
            edges: Vec::new(),
            adjacency: HashMap::new(),
        }
    }

    /// Add a node and return its ID.
    pub fn add_node(&mut self, value: f64, critical_type: CriticalPointType) -> usize {
        let id = self.nodes.len();
        self.nodes.push(ReebNode {
            id,
            value,
            critical_type,
            position: None,
        });
        id
    }

    /// Add a node with a 3D position.
    pub fn add_node_with_position(
        &mut self,
        value: f64,
        critical_type: CriticalPointType,
        position: [f64; 3],
    ) -> usize {
        let id = self.nodes.len();
        self.nodes.push(ReebNode {
            id,
            value,
            critical_type,
            position: Some(position),
        });
        id
    }

    /// Add an edge between two nodes.
    pub fn add_edge(&mut self, source: usize, target: usize) -> usize {
        let persistence = if source < self.nodes.len() && target < self.nodes.len() {
            (self.nodes[target].value - self.nodes[source].value).abs()
        } else {
            0.0
        };
        let idx = self.edges.len();
        self.edges.push(ReebEdge {
            source,
            target,
            component_count: 1,
            persistence,
        });
        self.adjacency.entry(source).or_default().push(idx);
        self.adjacency.entry(target).or_default().push(idx);
        idx
    }

    /// Number of nodes.
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }

    /// Number of edges.
    pub fn edge_count(&self) -> usize {
        self.edges.len()
    }

    /// Edges incident to a node.
    pub fn incident_edges(&self, node_id: usize) -> Vec<usize> {
        self.adjacency.get(&node_id).cloned().unwrap_or_default()
    }

    /// Degree of a node (number of incident edges).
    pub fn degree(&self, node_id: usize) -> usize {
        self.incident_edges(node_id).len()
    }

    /// Nodes of a given critical type.
    pub fn nodes_of_type(&self, ct: CriticalPointType) -> Vec<&ReebNode> {
        self.nodes
            .iter()
            .filter(|n| n.critical_type == ct)
            .collect()
    }

    /// First Betti number (β₁ = edges − nodes + connected_components).
    pub fn betti_1(&self) -> i32 {
        let cc = self.connected_components();
        self.edges.len() as i32 - self.nodes.len() as i32 + cc as i32
    }

    /// Number of connected components via BFS.
    pub fn connected_components(&self) -> usize {
        let mut visited = HashSet::new();
        let mut count = 0usize;
        for node in &self.nodes {
            if visited.contains(&node.id) {
                continue;
            }
            count += 1;
            let mut queue = VecDeque::new();
            queue.push_back(node.id);
            visited.insert(node.id);
            while let Some(cur) = queue.pop_front() {
                for &eidx in self.adjacency.get(&cur).unwrap_or(&Vec::new()) {
                    let edge = &self.edges[eidx];
                    let neighbor = if edge.source == cur {
                        edge.target
                    } else {
                        edge.source
                    };
                    if visited.insert(neighbor) {
                        queue.push_back(neighbor);
                    }
                }
            }
        }
        count
    }

    /// Simplify the Reeb graph by removing edges with persistence below a threshold.
    pub fn simplify(&mut self, threshold: f64) {
        let keep: Vec<bool> = self
            .edges
            .iter()
            .map(|e| e.persistence >= threshold)
            .collect();
        let mut new_edges = Vec::new();
        let mut new_adj: HashMap<usize, Vec<usize>> = HashMap::new();
        for (i, edge) in self.edges.iter().enumerate() {
            if keep[i] {
                let new_idx = new_edges.len();
                new_edges.push(edge.clone());
                new_adj.entry(edge.source).or_default().push(new_idx);
                new_adj.entry(edge.target).or_default().push(new_idx);
            }
        }
        self.edges = new_edges;
        self.adjacency = new_adj;
    }

    /// Build a Reeb graph from a triangulated surface with vertex scalar values.
    ///
    /// `vertices` are (x, y, z, scalar) tuples; `triangles` are index triples.
    pub fn from_triangle_mesh(
        vertices: &[(f64, f64, f64, f64)],
        triangles: &[(usize, usize, usize)],
    ) -> Self {
        let mut graph = ReebGraph::new();
        if vertices.is_empty() {
            return graph;
        }

        // Sort vertices by scalar value
        let mut sorted_indices: Vec<usize> = (0..vertices.len()).collect();
        sorted_indices.sort_by(|&a, &b| {
            vertices[a]
                .3
                .partial_cmp(&vertices[b].3)
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        // Build adjacency from triangles
        let mut vert_adj: HashMap<usize, HashSet<usize>> = HashMap::new();
        for &(a, b, c) in triangles {
            for &(u, v) in &[(a, b), (a, c), (b, c)] {
                vert_adj.entry(u).or_default().insert(v);
                vert_adj.entry(v).or_default().insert(u);
            }
        }

        // Determine critical type for each vertex
        let mut vertex_to_node: HashMap<usize, usize> = HashMap::new();
        for &vi in &sorted_indices {
            let val = vertices[vi].3;
            let neighbors = vert_adj.get(&vi).cloned().unwrap_or_default();
            let lower: Vec<usize> = neighbors
                .iter()
                .filter(|&&n| vertices[n].3 < val)
                .copied()
                .collect();
            let upper: Vec<usize> = neighbors
                .iter()
                .filter(|&&n| vertices[n].3 > val)
                .copied()
                .collect();

            let ct = if lower.is_empty() {
                CriticalPointType::Minimum
            } else if upper.is_empty() {
                CriticalPointType::Maximum
            } else {
                // Count connected components of lower link
                let lower_components = count_components_in_set(&lower, &vert_adj);
                if lower_components > 1 {
                    CriticalPointType::Saddle
                } else {
                    CriticalPointType::Regular
                }
            };

            // Only add non-regular or all for augmented version
            if ct != CriticalPointType::Regular {
                let nid = graph.add_node_with_position(
                    val,
                    ct,
                    [vertices[vi].0, vertices[vi].1, vertices[vi].2],
                );
                vertex_to_node.insert(vi, nid);
            }
        }

        // Connect critical points: for each critical vertex, connect to the
        // next critical vertex reachable by ascending/descending.
        let critical_sorted: Vec<usize> = sorted_indices
            .iter()
            .filter(|vi| vertex_to_node.contains_key(vi))
            .copied()
            .collect();

        for window in critical_sorted.windows(2) {
            let src = vertex_to_node[&window[0]];
            let tgt = vertex_to_node[&window[1]];
            graph.add_edge(src, tgt);
        }

        graph
    }
}

impl Default for ReebGraph {
    fn default() -> Self {
        Self::new()
    }
}

/// Count connected components within a subset of vertices.
fn count_components_in_set(subset: &[usize], adjacency: &HashMap<usize, HashSet<usize>>) -> usize {
    let member_set: HashSet<usize> = subset.iter().copied().collect();
    let mut visited = HashSet::new();
    let mut count = 0usize;
    for &v in subset {
        if visited.contains(&v) {
            continue;
        }
        count += 1;
        let mut queue = VecDeque::new();
        queue.push_back(v);
        visited.insert(v);
        while let Some(cur) = queue.pop_front() {
            if let Some(neighbors) = adjacency.get(&cur) {
                for &n in neighbors {
                    if member_set.contains(&n) && visited.insert(n) {
                        queue.push_back(n);
                    }
                }
            }
        }
    }
    count
}

// ---------------------------------------------------------------------------
// MorseSmaleComplex
// ---------------------------------------------------------------------------

/// A critical point with its index (number of negative eigenvalues of the Hessian).
#[derive(Debug, Clone)]
pub struct MorseCriticalPoint {
    /// Unique identifier.
    pub id: usize,
    /// Position in domain.
    pub position: [f64; 3],
    /// Scalar function value.
    pub value: f64,
    /// Morse index (0 = min, 1 = 1-saddle, 2 = 2-saddle, 3 = max in 3D).
    pub index: usize,
    /// Whether this point has been cancelled during simplification.
    pub cancelled: bool,
}

/// A separatrix connecting two critical points of adjacent Morse index.
#[derive(Debug, Clone)]
pub struct Separatrix {
    /// Source critical point ID (lower index).
    pub source: usize,
    /// Destination critical point ID (higher index).
    pub destination: usize,
    /// Path as a sequence of positions.
    pub path: Vec<[f64; 3]>,
    /// Persistence of this separatrix (|f(dest) − f(source)|).
    pub persistence: f64,
}

/// Direction of a discrete gradient arrow.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GradientDirection {
    /// From a vertex to an edge.
    VertexToEdge,
    /// From an edge to a triangle.
    EdgeToTriangle,
    /// From a triangle to a tetrahedron.
    TriangleToTet,
}

/// A discrete gradient vector in the combinatorial Morse theory sense.
#[derive(Debug, Clone)]
pub struct DiscreteGradient {
    /// Pairs of (lower-dim cell, higher-dim cell).
    pub pairs: Vec<(usize, usize)>,
    /// Direction type for each pair.
    pub directions: Vec<GradientDirection>,
    /// Critical cells (unpaired).
    pub critical_cells: Vec<usize>,
}

impl DiscreteGradient {
    /// Create a new empty discrete gradient.
    pub fn new() -> Self {
        Self {
            pairs: Vec::new(),
            directions: Vec::new(),
            critical_cells: Vec::new(),
        }
    }

    /// Add a gradient pair.
    pub fn add_pair(&mut self, lower: usize, upper: usize, direction: GradientDirection) {
        self.pairs.push((lower, upper));
        self.directions.push(direction);
    }

    /// Mark a cell as critical (unpaired).
    pub fn add_critical_cell(&mut self, cell: usize) {
        self.critical_cells.push(cell);
    }

    /// Number of gradient pairs.
    pub fn pair_count(&self) -> usize {
        self.pairs.len()
    }

    /// Number of critical cells.
    pub fn critical_count(&self) -> usize {
        self.critical_cells.len()
    }
}

impl Default for DiscreteGradient {
    fn default() -> Self {
        Self::new()
    }
}

/// A Morse-Smale complex decomposition.
#[derive(Debug, Clone)]
pub struct MorseSmaleComplex {
    /// Critical points.
    pub critical_points: Vec<MorseCriticalPoint>,
    /// Separatrices connecting critical points.
    pub separatrices: Vec<Separatrix>,
    /// Combinatorial discrete gradient field (if computed).
    pub gradient: Option<DiscreteGradient>,
    /// Ascending manifold assignment: vertex → critical point ID.
    pub ascending_manifold: HashMap<usize, usize>,
    /// Descending manifold assignment: vertex → critical point ID.
    pub descending_manifold: HashMap<usize, usize>,
}

impl MorseSmaleComplex {
    /// Create an empty complex.
    pub fn new() -> Self {
        Self {
            critical_points: Vec::new(),
            separatrices: Vec::new(),
            gradient: None,
            ascending_manifold: HashMap::new(),
            descending_manifold: HashMap::new(),
        }
    }

    /// Add a critical point.
    pub fn add_critical_point(&mut self, position: [f64; 3], value: f64, index: usize) -> usize {
        let id = self.critical_points.len();
        self.critical_points.push(MorseCriticalPoint {
            id,
            position,
            value,
            index,
            cancelled: false,
        });
        id
    }

    /// Add a separatrix.
    pub fn add_separatrix(&mut self, source: usize, destination: usize, path: Vec<[f64; 3]>) {
        let persistence =
            if source < self.critical_points.len() && destination < self.critical_points.len() {
                (self.critical_points[destination].value - self.critical_points[source].value).abs()
            } else {
                0.0
            };
        self.separatrices.push(Separatrix {
            source,
            destination,
            path,
            persistence,
        });
    }

    /// Number of critical points of a given Morse index.
    pub fn count_by_index(&self, index: usize) -> usize {
        self.critical_points
            .iter()
            .filter(|cp| cp.index == index && !cp.cancelled)
            .count()
    }

    /// Verify the Morse inequalities (weak form): the number of critical points
    /// of each index is at least the corresponding Betti number.
    ///
    /// Returns the alternating sum Σ(-1)^k c_k which should equal the Euler
    /// characteristic.
    pub fn euler_characteristic(&self) -> i32 {
        let mut chi = 0i32;
        for cp in &self.critical_points {
            if !cp.cancelled {
                if cp.index % 2 == 0 {
                    chi += 1;
                } else {
                    chi -= 1;
                }
            }
        }
        chi
    }

    /// Assign ascending manifolds by steepest ascent from each vertex.
    ///
    /// `scalar_values`: function value at each vertex.
    /// `adjacency`: vertex → list of neighbor vertices.
    pub fn compute_ascending_manifolds(&mut self, scalar_values: &[f64], adjacency: &[Vec<usize>]) {
        self.ascending_manifold.clear();
        // Find maxima
        let max_ids: Vec<usize> = self
            .critical_points
            .iter()
            .filter(|cp| {
                cp.index == 0 || {
                    // Check if it's a local max in the adjacency
                    let vi = cp.id;
                    vi < adjacency.len()
                        && adjacency[vi].iter().all(|&n| {
                            n >= scalar_values.len() || scalar_values[n] <= scalar_values[vi]
                        })
                }
            })
            .map(|cp| cp.id)
            .collect();

        for vi in 0..scalar_values.len() {
            // Follow steepest ascent
            let mut current = vi;
            let mut visited = HashSet::new();
            loop {
                if !visited.insert(current) {
                    break;
                }
                if current >= adjacency.len() {
                    break;
                }
                let mut best = current;
                let mut best_val = scalar_values[current];
                for &n in &adjacency[current] {
                    if n < scalar_values.len() && scalar_values[n] > best_val {
                        best_val = scalar_values[n];
                        best = n;
                    }
                }
                if best == current {
                    break;
                }
                current = best;
            }
            // Check if `current` is a known max
            if max_ids.contains(&current) {
                self.ascending_manifold.insert(vi, current);
            }
        }
    }

    /// Compute descending manifolds by steepest descent from each vertex.
    pub fn compute_descending_manifolds(
        &mut self,
        scalar_values: &[f64],
        adjacency: &[Vec<usize>],
    ) {
        self.descending_manifold.clear();
        for vi in 0..scalar_values.len() {
            let mut current = vi;
            let mut visited = HashSet::new();
            loop {
                if !visited.insert(current) {
                    break;
                }
                if current >= adjacency.len() {
                    break;
                }
                let mut best = current;
                let mut best_val = scalar_values[current];
                for &n in &adjacency[current] {
                    if n < scalar_values.len() && scalar_values[n] < best_val {
                        best_val = scalar_values[n];
                        best = n;
                    }
                }
                if best == current {
                    break;
                }
                current = best;
            }
            self.descending_manifold.insert(vi, current);
        }
    }

    /// Cancel a persistence pair below a threshold.
    pub fn cancel_pairs(&mut self, threshold: f64) {
        let mut to_cancel = Vec::new();
        for sep in &self.separatrices {
            if sep.persistence < threshold {
                to_cancel.push((sep.source, sep.destination));
            }
        }
        for (s, d) in to_cancel {
            if s < self.critical_points.len() {
                self.critical_points[s].cancelled = true;
            }
            if d < self.critical_points.len() {
                self.critical_points[d].cancelled = true;
            }
        }
        self.separatrices.retain(|sep| sep.persistence >= threshold);
    }

    /// Active (non-cancelled) critical points.
    pub fn active_critical_points(&self) -> Vec<&MorseCriticalPoint> {
        self.critical_points
            .iter()
            .filter(|cp| !cp.cancelled)
            .collect()
    }
}

impl Default for MorseSmaleComplex {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// ContourTree
// ---------------------------------------------------------------------------

/// A node in a contour tree.
#[derive(Debug, Clone)]
pub struct ContourTreeNode {
    /// Vertex index in the original mesh.
    pub vertex: usize,
    /// Scalar value.
    pub value: f64,
    /// Type of critical point.
    pub critical_type: CriticalPointType,
}

/// An arc in a contour tree.
#[derive(Debug, Clone)]
pub struct ContourTreeArc {
    /// Lower endpoint (smaller function value).
    pub lower: usize,
    /// Upper endpoint (larger function value).
    pub upper: usize,
    /// Number of regular vertices on this arc.
    pub regular_count: usize,
}

/// A contour tree: the result of merging join tree and split tree.
#[derive(Debug, Clone)]
pub struct ContourTree {
    /// Nodes (critical vertices).
    pub nodes: Vec<ContourTreeNode>,
    /// Arcs connecting nodes.
    pub arcs: Vec<ContourTreeArc>,
    /// Whether this is an augmented contour tree (includes regular nodes).
    pub augmented: bool,
}

impl ContourTree {
    /// Create an empty contour tree.
    pub fn new() -> Self {
        Self {
            nodes: Vec::new(),
            arcs: Vec::new(),
            augmented: false,
        }
    }

    /// Build a contour tree from a scalar function on a graph.
    ///
    /// `values`: scalar value at each vertex.
    /// `adjacency`: adjacency list for the graph.
    pub fn build(values: &[f64], adjacency: &[Vec<usize>]) -> Self {
        let n = values.len();
        if n == 0 {
            return Self::new();
        }

        // Sort vertices by value
        let mut sorted: Vec<usize> = (0..n).collect();
        sorted.sort_by(|&a, &b| {
            values[a]
                .partial_cmp(&values[b])
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        // Build join tree (sweep from min to max)
        let join_tree = Self::build_merge_tree(&sorted, values, adjacency);

        // Build split tree (sweep from max to min)
        let mut rev_sorted = sorted.clone();
        rev_sorted.reverse();
        let split_tree = Self::build_merge_tree(&rev_sorted, values, adjacency);

        // Merge join and split trees
        Self::merge_trees(&join_tree, &split_tree, values)
    }

    /// Build a merge tree by processing vertices in the given order.
    fn build_merge_tree(
        order: &[usize],
        values: &[f64],
        adjacency: &[Vec<usize>],
    ) -> BTreeMap<usize, Vec<usize>> {
        let n = values.len();
        let mut uf_parent: Vec<usize> = (0..n).collect();
        let mut uf_rank = vec![0u32; n];
        let mut processed = vec![false; n];
        let mut tree: BTreeMap<usize, Vec<usize>> = BTreeMap::new();

        for &v in order {
            processed[v] = true;
            let mut components = HashSet::new();
            for &u in &adjacency[v] {
                if processed[u] {
                    let root = uf_find(&mut uf_parent, u);
                    components.insert(root);
                }
            }
            if components.is_empty() {
                // Leaf node
                tree.insert(v, Vec::new());
            } else {
                // Merge components
                let roots: Vec<usize> = components.into_iter().collect();
                for &r in &roots {
                    uf_union(&mut uf_parent, &mut uf_rank, v, r);
                    tree.entry(v).or_default().push(r);
                }
            }
        }
        tree
    }

    /// Merge join and split trees into a contour tree.
    fn merge_trees(
        join: &BTreeMap<usize, Vec<usize>>,
        split: &BTreeMap<usize, Vec<usize>>,
        values: &[f64],
    ) -> Self {
        let mut ct = ContourTree::new();
        let mut vertex_to_node: HashMap<usize, usize> = HashMap::new();

        // Collect all critical vertices (leaves or branch points in either tree)
        let mut critical: HashSet<usize> = HashSet::new();
        for (&v, children) in join {
            if children.is_empty() || children.len() >= 2 {
                critical.insert(v);
            }
        }
        for (&v, children) in split {
            if children.is_empty() || children.len() >= 2 {
                critical.insert(v);
            }
        }

        // If no critical points found, add all vertices from either tree
        if critical.is_empty() {
            for &v in join.keys() {
                critical.insert(v);
            }
        }

        // Create nodes
        let mut sorted_critical: Vec<usize> = critical.into_iter().collect();
        sorted_critical.sort_by(|&a, &b| {
            values[a]
                .partial_cmp(&values[b])
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        for &v in &sorted_critical {
            let nid = ct.nodes.len();
            let ct_type = if join.get(&v).is_some_and(|c| c.is_empty()) {
                CriticalPointType::Minimum
            } else if split.get(&v).is_some_and(|c| c.is_empty()) {
                CriticalPointType::Maximum
            } else {
                CriticalPointType::Saddle
            };
            ct.nodes.push(ContourTreeNode {
                vertex: v,
                value: values[v],
                critical_type: ct_type,
            });
            vertex_to_node.insert(v, nid);
        }

        // Create arcs between consecutive critical points in sorted order
        for window in sorted_critical.windows(2) {
            let lower = vertex_to_node[&window[0]];
            let upper = vertex_to_node[&window[1]];
            ct.arcs.push(ContourTreeArc {
                lower,
                upper,
                regular_count: 0,
            });
        }

        ct
    }

    /// Number of nodes.
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }

    /// Number of arcs.
    pub fn arc_count(&self) -> usize {
        self.arcs.len()
    }

    /// Leaf nodes (degree 1).
    pub fn leaves(&self) -> Vec<usize> {
        let mut degree = vec![0usize; self.nodes.len()];
        for arc in &self.arcs {
            degree[arc.lower] += 1;
            degree[arc.upper] += 1;
        }
        (0..self.nodes.len()).filter(|&i| degree[i] == 1).collect()
    }

    /// Persistence of each arc (difference in scalar values).
    pub fn arc_persistences(&self) -> Vec<f64> {
        self.arcs
            .iter()
            .map(|a| (self.nodes[a.upper].value - self.nodes[a.lower].value).abs())
            .collect()
    }
}

impl Default for ContourTree {
    fn default() -> Self {
        Self::new()
    }
}

/// Union-Find: find with path compression.
fn uf_find(parent: &mut [usize], x: usize) -> usize {
    if parent[x] != x {
        parent[x] = uf_find(parent, parent[x]);
    }
    parent[x]
}

/// Union-Find: union by rank.
fn uf_union(parent: &mut [usize], rank: &mut [u32], a: usize, b: usize) {
    let ra = uf_find(parent, a);
    let rb = uf_find(parent, b);
    if ra == rb {
        return;
    }
    if rank[ra] < rank[rb] {
        parent[ra] = rb;
    } else if rank[ra] > rank[rb] {
        parent[rb] = ra;
    } else {
        parent[rb] = ra;
        rank[ra] += 1;
    }
}

// ---------------------------------------------------------------------------
// TopologicalSimplification
// ---------------------------------------------------------------------------

/// Topological simplification by cancelling persistence pairs.
#[derive(Debug, Clone)]
pub struct TopologicalSimplification {
    /// Persistence threshold below which pairs are cancelled.
    pub threshold: f64,
    /// Number of pairs cancelled.
    pub cancelled_count: usize,
    /// Remaining features after simplification.
    pub remaining_features: Vec<TopologicalFeature>,
}

impl TopologicalSimplification {
    /// Simplify a persistence diagram by removing pairs below a threshold.
    pub fn simplify_diagram(diagram: &PersistenceDiagram, threshold: f64) -> Self {
        let remaining: Vec<TopologicalFeature> = diagram
            .features
            .iter()
            .filter(|f| f.persistence() >= threshold || f.is_essential())
            .cloned()
            .collect();
        let cancelled = diagram.features.len() - remaining.len();
        Self {
            threshold,
            cancelled_count: cancelled,
            remaining_features: remaining,
        }
    }

    /// Simplify a Reeb graph in place.
    pub fn simplify_reeb_graph(graph: &mut ReebGraph, threshold: f64) -> usize {
        let before = graph.edge_count();
        graph.simplify(threshold);
        before - graph.edge_count()
    }

    /// Simplify a Morse-Smale complex by cancelling pairs.
    pub fn simplify_morse_smale(complex: &mut MorseSmaleComplex, threshold: f64) -> usize {
        let before = complex.active_critical_points().len();
        complex.cancel_pairs(threshold);
        let after = complex.active_critical_points().len();
        before - after
    }
}

// ---------------------------------------------------------------------------
// TopologyVizStyle
// ---------------------------------------------------------------------------

/// Color specification for topology visualization (RGBA, each component 0–1).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TopoColor {
    /// Red.
    pub r: f32,
    /// Green.
    pub g: f32,
    /// Blue.
    pub b: f32,
    /// Alpha.
    pub a: f32,
}

impl TopoColor {
    /// Create a new color.
    pub fn new(r: f32, g: f32, b: f32, a: f32) -> Self {
        Self { r, g, b, a }
    }

    /// Red color.
    pub fn red() -> Self {
        Self::new(1.0, 0.0, 0.0, 1.0)
    }

    /// Blue color.
    pub fn blue() -> Self {
        Self::new(0.0, 0.0, 1.0, 1.0)
    }

    /// Green color.
    pub fn green() -> Self {
        Self::new(0.0, 1.0, 0.0, 1.0)
    }

    /// Interpolate between two colors.
    pub fn lerp(&self, other: &TopoColor, t: f32) -> TopoColor {
        TopoColor {
            r: self.r + (other.r - self.r) * t,
            g: self.g + (other.g - self.g) * t,
            b: self.b + (other.b - self.b) * t,
            a: self.a + (other.a - self.a) * t,
        }
    }
}

/// Visualization style for topological features.
#[derive(Debug, Clone)]
pub struct TopologyVizStyle {
    /// Color per homology dimension (index = dimension).
    pub dimension_colors: Vec<TopoColor>,
    /// Base point size for diagram points.
    pub base_point_size: f32,
    /// Scale factor: point size = base_size * (persistence / max_persistence)^scale.
    pub persistence_scale: f32,
    /// Minimum alpha for low-persistence features.
    pub min_alpha: f32,
    /// Maximum alpha for high-persistence features.
    pub max_alpha: f32,
    /// Whether to draw the diagonal line on persistence diagrams.
    pub show_diagonal: bool,
    /// Line width for barcode bars.
    pub barcode_line_width: f32,
    /// Whether to use log scale for the filtration axis.
    pub log_scale: bool,
}

impl Default for TopologyVizStyle {
    fn default() -> Self {
        Self {
            dimension_colors: vec![
                TopoColor::new(0.2, 0.4, 0.8, 1.0), // H0: blue
                TopoColor::new(0.8, 0.2, 0.2, 1.0), // H1: red
                TopoColor::new(0.2, 0.8, 0.2, 1.0), // H2: green
            ],
            base_point_size: 5.0,
            persistence_scale: 1.0,
            min_alpha: 0.1,
            max_alpha: 1.0,
            show_diagonal: true,
            barcode_line_width: 2.0,
            log_scale: false,
        }
    }
}

impl TopologyVizStyle {
    /// Color for a given homology dimension.
    pub fn color_for_dimension(&self, dim: usize) -> TopoColor {
        if dim < self.dimension_colors.len() {
            self.dimension_colors[dim]
        } else {
            TopoColor::new(0.5, 0.5, 0.5, 1.0)
        }
    }

    /// Point size for a feature given its persistence and the max persistence.
    pub fn point_size(&self, persistence: f64, max_persistence: f64) -> f32 {
        if max_persistence <= 0.0 {
            return self.base_point_size;
        }
        let ratio = (persistence / max_persistence) as f32;
        self.base_point_size * ratio.powf(self.persistence_scale)
    }

    /// Alpha for a feature given its relative persistence.
    pub fn alpha_for_persistence(&self, persistence: f64, max_persistence: f64) -> f32 {
        if max_persistence <= 0.0 {
            return self.max_alpha;
        }
        let ratio = (persistence / max_persistence) as f32;
        self.min_alpha + (self.max_alpha - self.min_alpha) * ratio
    }

    /// Generate styled point data for a persistence diagram.
    ///
    /// Returns `(x, y, size, color)` tuples for each non-essential feature.
    pub fn styled_diagram_points(
        &self,
        diagram: &PersistenceDiagram,
    ) -> Vec<(f64, f64, f32, TopoColor)> {
        let max_p = diagram.max_persistence();
        diagram
            .features
            .iter()
            .filter(|f| !f.is_essential())
            .map(|f| {
                let color = self.color_for_dimension(f.dimension);
                let size = self.point_size(f.persistence(), max_p);
                let alpha = self.alpha_for_persistence(f.persistence(), max_p);
                let styled_color = TopoColor::new(color.r, color.g, color.b, alpha);
                (f.birth, f.death, size, styled_color)
            })
            .collect()
    }

    /// Generate styled barcode bars.
    ///
    /// Returns `(y_position, birth, death, color)` tuples.
    pub fn styled_barcode_bars(
        &self,
        barcode: &PersistenceBarcode,
    ) -> Vec<(f64, f64, f64, TopoColor)> {
        let intervals = barcode.intervals();
        let max_p = barcode.diagram.max_persistence();
        intervals
            .iter()
            .enumerate()
            .map(|(i, &(birth, death, dim))| {
                let color = self.color_for_dimension(dim);
                let alpha = self.alpha_for_persistence((death - birth).abs(), max_p);
                let styled = TopoColor::new(color.r, color.g, color.b, alpha);
                (i as f64, birth, death, styled)
            })
            .collect()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- TopologicalFeature ---

    #[test]
    fn test_feature_persistence() {
        let f = TopologicalFeature::new(1, 0.5, 2.5);
        assert!((f.persistence() - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_feature_essential() {
        let f = TopologicalFeature::new(0, 0.0, f64::INFINITY);
        assert!(f.is_essential());
        assert!(!TopologicalFeature::new(0, 0.0, 1.0).is_essential());
    }

    #[test]
    fn test_feature_midlife() {
        let f = TopologicalFeature::new(0, 1.0, 5.0);
        assert!((f.midlife() - 3.0).abs() < 1e-10);
    }

    // --- PersistenceDiagram ---

    #[test]
    fn test_diagram_from_tuples() {
        let pd = PersistenceDiagram::from_tuples(&[(0, 0.0, 1.0), (0, 0.0, 3.0), (1, 1.0, 2.0)]);
        assert_eq!(pd.len(), 3);
        assert_eq!(pd.dimension(0).len(), 2);
        assert_eq!(pd.dimension(1).len(), 1);
    }

    #[test]
    fn test_diagram_max_persistence() {
        let pd = PersistenceDiagram::from_tuples(&[(0, 0.0, 1.0), (0, 0.0, 5.0), (1, 2.0, 3.0)]);
        assert!((pd.max_persistence() - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_diagram_significant_features() {
        let pd = PersistenceDiagram::from_tuples(&[(0, 0.0, 0.1), (0, 0.0, 5.0), (1, 1.0, 1.2)]);
        let sig = pd.significant_features(1.0);
        assert_eq!(sig.len(), 1);
        assert!((sig[0].death - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_diagram_bottleneck_identical() {
        let pd = PersistenceDiagram::from_tuples(&[(0, 0.0, 1.0), (1, 0.5, 2.0)]);
        let dist = pd.bottleneck_distance(&pd);
        assert!(dist < 1e-10);
    }

    #[test]
    fn test_diagram_bottleneck_shifted() {
        let pd1 = PersistenceDiagram::from_tuples(&[(0, 0.0, 2.0)]);
        let pd2 = PersistenceDiagram::from_tuples(&[(0, 0.0, 3.0)]);
        let dist = pd1.bottleneck_distance(&pd2);
        assert!((dist - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_diagram_betti_numbers() {
        let pd = PersistenceDiagram::from_tuples(&[
            (0, 0.0, f64::INFINITY),
            (0, 0.0, 0.5),
            (1, 0.3, 0.8),
        ]);
        let betti = pd.betti_numbers(0.4);
        assert_eq!(*betti.get(&0).unwrap_or(&0), 2); // both H0 alive
        assert_eq!(*betti.get(&1).unwrap_or(&0), 1); // H1 alive
    }

    #[test]
    fn test_diagram_wasserstein_zero_for_same() {
        let pd = PersistenceDiagram::from_tuples(&[(0, 0.0, 1.0)]);
        assert!(pd.wasserstein_distance(&pd) < 1e-10);
    }

    // --- PersistenceBarcode ---

    #[test]
    fn test_barcode_interval_count() {
        let pd = PersistenceDiagram::from_tuples(&[(0, 0.0, 1.0), (0, 0.5, 2.0), (1, 1.0, 3.0)]);
        let bc = PersistenceBarcode::new(pd, 5.0);
        assert_eq!(bc.interval_count(), 3);
    }

    #[test]
    fn test_barcode_sort_by_persistence() {
        let pd = PersistenceDiagram::from_tuples(&[
            (0, 0.0, 1.0), // pers = 1
            (0, 0.0, 5.0), // pers = 5
            (1, 1.0, 2.0), // pers = 1
        ]);
        let bc = PersistenceBarcode::new(pd, 10.0).with_sort_order(BarcodeSortOrder::Persistence);
        let intervals = bc.intervals();
        // First should have highest persistence
        assert!((intervals[0].1 - intervals[0].0 - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_barcode_caps_infinite_death() {
        let pd = PersistenceDiagram::from_tuples(&[(0, 0.0, f64::INFINITY)]);
        let bc = PersistenceBarcode::new(pd, 10.0);
        let intervals = bc.intervals();
        assert!((intervals[0].1 - 10.0).abs() < 1e-10);
    }

    #[test]
    fn test_barcode_dimension_filter() {
        let pd = PersistenceDiagram::from_tuples(&[(0, 0.0, 1.0), (1, 0.5, 2.0), (1, 1.0, 3.0)]);
        let bc = PersistenceBarcode::new(pd, 5.0);
        assert_eq!(bc.intervals_for_dimension(1).len(), 2);
    }

    // --- ReebGraph ---

    #[test]
    fn test_reeb_graph_add_nodes_edges() {
        let mut rg = ReebGraph::new();
        let n0 = rg.add_node(0.0, CriticalPointType::Minimum);
        let n1 = rg.add_node(1.0, CriticalPointType::Saddle);
        let n2 = rg.add_node(2.0, CriticalPointType::Maximum);
        rg.add_edge(n0, n1);
        rg.add_edge(n1, n2);
        assert_eq!(rg.node_count(), 3);
        assert_eq!(rg.edge_count(), 2);
    }

    #[test]
    fn test_reeb_graph_betti_1_cycle() {
        // A cycle: min -> saddle -> max, min -> saddle2 -> max
        let mut rg = ReebGraph::new();
        let n0 = rg.add_node(0.0, CriticalPointType::Minimum);
        let n1 = rg.add_node(0.5, CriticalPointType::Saddle);
        let n2 = rg.add_node(0.5, CriticalPointType::Saddle);
        let n3 = rg.add_node(1.0, CriticalPointType::Maximum);
        rg.add_edge(n0, n1);
        rg.add_edge(n0, n2);
        rg.add_edge(n1, n3);
        rg.add_edge(n2, n3);
        // β₁ = E - V + CC = 4 - 4 + 1 = 1
        assert_eq!(rg.betti_1(), 1);
    }

    #[test]
    fn test_reeb_graph_torus() {
        // Torus Reeb graph: min, two saddles, max with β₁ = 2
        let mut rg = ReebGraph::new();
        let n0 = rg.add_node(0.0, CriticalPointType::Minimum);
        let s1 = rg.add_node(1.0, CriticalPointType::Saddle);
        let s2 = rg.add_node(2.0, CriticalPointType::Saddle);
        let n3 = rg.add_node(3.0, CriticalPointType::Maximum);
        // Two paths from min to s1
        rg.add_edge(n0, s1);
        rg.add_edge(n0, s1);
        // Two paths from s1 to s2
        rg.add_edge(s1, s2);
        rg.add_edge(s1, s2);
        // One path from s2 to max
        rg.add_edge(s2, n3);
        // β₁ = 5 - 4 + 1 = 2
        assert_eq!(rg.betti_1(), 2);
    }

    #[test]
    fn test_reeb_graph_simplification() {
        let mut rg = ReebGraph::new();
        let n0 = rg.add_node(0.0, CriticalPointType::Minimum);
        let n1 = rg.add_node(0.1, CriticalPointType::Saddle);
        let n2 = rg.add_node(5.0, CriticalPointType::Maximum);
        rg.add_edge(n0, n1);
        rg.add_edge(n1, n2);
        assert_eq!(rg.edge_count(), 2);
        rg.simplify(1.0);
        // Edge (n0,n1) has persistence 0.1 < 1.0, should be removed
        assert_eq!(rg.edge_count(), 1);
    }

    #[test]
    fn test_reeb_graph_connected_components() {
        let mut rg = ReebGraph::new();
        rg.add_node(0.0, CriticalPointType::Minimum);
        rg.add_node(1.0, CriticalPointType::Maximum);
        // No edges → 2 components
        assert_eq!(rg.connected_components(), 2);
    }

    #[test]
    fn test_reeb_graph_from_mesh() {
        // Simple triangle mesh: 4 vertices forming a tetrahedron-like shape
        let vertices = vec![
            (0.0, 0.0, 0.0, 0.0), // min
            (1.0, 0.0, 0.0, 1.0),
            (0.5, 1.0, 0.0, 2.0),
            (0.5, 0.5, 1.0, 3.0), // max
        ];
        let triangles = vec![(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)];
        let rg = ReebGraph::from_triangle_mesh(&vertices, &triangles);
        assert!(rg.node_count() >= 2); // at least a min and max
    }

    // --- MorseSmaleComplex ---

    #[test]
    fn test_morse_smale_euler_characteristic() {
        let mut ms = MorseSmaleComplex::new();
        ms.add_critical_point([0.0, 0.0, 0.0], 0.0, 0); // min (index 0)
        ms.add_critical_point([1.0, 0.0, 0.0], 1.0, 1); // saddle (index 1)
        ms.add_critical_point([0.5, 1.0, 0.0], 2.0, 2); // max (index 2)
        // χ = c0 - c1 + c2 = 1 - 1 + 1 = 1
        assert_eq!(ms.euler_characteristic(), 1);
    }

    #[test]
    fn test_morse_smale_cancel_pairs() {
        let mut ms = MorseSmaleComplex::new();
        let a = ms.add_critical_point([0.0, 0.0, 0.0], 0.0, 0);
        let b = ms.add_critical_point([1.0, 0.0, 0.0], 0.5, 1);
        let c = ms.add_critical_point([0.5, 1.0, 0.0], 5.0, 2);
        ms.add_separatrix(a, b, vec![]);
        ms.add_separatrix(b, c, vec![]);
        assert_eq!(ms.active_critical_points().len(), 3);
        ms.cancel_pairs(1.0);
        // a-b pair (persistence 0.5) cancelled
        assert_eq!(ms.active_critical_points().len(), 1);
    }

    #[test]
    fn test_morse_smale_count_by_index() {
        let mut ms = MorseSmaleComplex::new();
        ms.add_critical_point([0.0, 0.0, 0.0], 0.0, 0);
        ms.add_critical_point([1.0, 0.0, 0.0], 1.0, 0);
        ms.add_critical_point([2.0, 0.0, 0.0], 2.0, 1);
        assert_eq!(ms.count_by_index(0), 2);
        assert_eq!(ms.count_by_index(1), 1);
        assert_eq!(ms.count_by_index(2), 0);
    }

    #[test]
    fn test_discrete_gradient() {
        let mut dg = DiscreteGradient::new();
        dg.add_pair(0, 1, GradientDirection::VertexToEdge);
        dg.add_pair(2, 3, GradientDirection::EdgeToTriangle);
        dg.add_critical_cell(4);
        assert_eq!(dg.pair_count(), 2);
        assert_eq!(dg.critical_count(), 1);
    }

    // --- ContourTree ---

    #[test]
    fn test_contour_tree_simple_path() {
        // Simple path graph: 0 -- 1 -- 2 -- 3
        let values = vec![0.0, 1.0, 2.0, 3.0];
        let adj = vec![vec![1], vec![0, 2], vec![1, 3], vec![2]];
        let ct = ContourTree::build(&values, &adj);
        assert!(ct.node_count() >= 2);
        assert!(ct.arc_count() >= 1);
    }

    #[test]
    fn test_contour_tree_height_function() {
        // Grid height function: 4 vertices forming a diamond
        // 0(h=0) -- 1(h=1) -- 3(h=3)
        //              |
        //           2(h=2)
        let values = vec![0.0, 1.0, 2.0, 3.0];
        let adj = vec![vec![1], vec![0, 2, 3], vec![1], vec![1]];
        let ct = ContourTree::build(&values, &adj);
        // Should have min(0), branch(1), and leaves(2,3)
        assert!(ct.node_count() >= 2);
    }

    #[test]
    fn test_contour_tree_arc_persistences() {
        let values = vec![0.0, 5.0, 10.0];
        let adj = vec![vec![1], vec![0, 2], vec![1]];
        let ct = ContourTree::build(&values, &adj);
        let pers = ct.arc_persistences();
        // Total persistence from min to max should sum to 10
        let total: f64 = pers.iter().sum();
        assert!((total - 10.0).abs() < 1e-10);
    }

    #[test]
    fn test_contour_tree_leaves() {
        let values = vec![0.0, 1.0, 2.0];
        let adj = vec![vec![1], vec![0, 2], vec![1]];
        let ct = ContourTree::build(&values, &adj);
        let leaves = ct.leaves();
        assert!(leaves.len() >= 2);
    }

    // --- TopologicalSimplification ---

    #[test]
    fn test_simplification_removes_low_persistence() {
        let pd = PersistenceDiagram::from_tuples(&[(0, 0.0, 0.1), (0, 0.0, 5.0), (1, 1.0, 1.05)]);
        let simp = TopologicalSimplification::simplify_diagram(&pd, 0.5);
        assert_eq!(simp.remaining_features.len(), 1);
        assert_eq!(simp.cancelled_count, 2);
    }

    #[test]
    fn test_simplification_keeps_essential() {
        let pd = PersistenceDiagram::from_tuples(&[(0, 0.0, f64::INFINITY), (0, 0.0, 0.01)]);
        let simp = TopologicalSimplification::simplify_diagram(&pd, 1.0);
        assert_eq!(simp.remaining_features.len(), 1);
        assert!(simp.remaining_features[0].is_essential());
    }

    #[test]
    fn test_simplification_reeb_graph() {
        let mut rg = ReebGraph::new();
        let n0 = rg.add_node(0.0, CriticalPointType::Minimum);
        let n1 = rg.add_node(0.05, CriticalPointType::Saddle);
        let n2 = rg.add_node(10.0, CriticalPointType::Maximum);
        rg.add_edge(n0, n1);
        rg.add_edge(n1, n2);
        let removed = TopologicalSimplification::simplify_reeb_graph(&mut rg, 1.0);
        assert_eq!(removed, 1);
    }

    // --- TopologyVizStyle ---

    #[test]
    fn test_viz_style_color_by_dimension() {
        let style = TopologyVizStyle::default();
        let c0 = style.color_for_dimension(0);
        let c1 = style.color_for_dimension(1);
        // H0 and H1 should have different colors
        assert!((c0.r - c1.r).abs() > 0.1 || (c0.b - c1.b).abs() > 0.1);
    }

    #[test]
    fn test_viz_style_point_size_scales() {
        let style = TopologyVizStyle::default();
        let small = style.point_size(0.1, 1.0);
        let big = style.point_size(1.0, 1.0);
        assert!(big > small);
    }

    #[test]
    fn test_viz_style_alpha_range() {
        let style = TopologyVizStyle::default();
        let low = style.alpha_for_persistence(0.01, 1.0);
        let high = style.alpha_for_persistence(1.0, 1.0);
        assert!(low >= style.min_alpha - 1e-6);
        assert!(high <= style.max_alpha + 1e-6);
        assert!(high > low);
    }

    #[test]
    fn test_viz_style_diagram_points() {
        let pd = PersistenceDiagram::from_tuples(&[(0, 0.0, 1.0), (1, 0.5, 2.0)]);
        let style = TopologyVizStyle::default();
        let points = style.styled_diagram_points(&pd);
        assert_eq!(points.len(), 2);
        // Check birth/death values
        assert!((points[0].0 - 0.0).abs() < 1e-10); // birth
        assert!((points[0].1 - 1.0).abs() < 1e-10); // death
    }

    #[test]
    fn test_viz_style_barcode_bars() {
        let pd = PersistenceDiagram::from_tuples(&[(0, 0.0, 1.0), (1, 0.5, 3.0)]);
        let bc = PersistenceBarcode::new(pd, 5.0);
        let style = TopologyVizStyle::default();
        let bars = style.styled_barcode_bars(&bc);
        assert_eq!(bars.len(), 2);
    }

    #[test]
    fn test_topo_color_lerp() {
        let a = TopoColor::red();
        let b = TopoColor::blue();
        let mid = a.lerp(&b, 0.5);
        assert!((mid.r - 0.5).abs() < 1e-6);
        assert!((mid.b - 0.5).abs() < 1e-6);
    }
}
