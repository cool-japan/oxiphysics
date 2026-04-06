#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Network and graph layout algorithms for physics visualisation.
//!
//! Provides force-directed (Fruchterman–Reingold), circular, and hierarchical
//! layout algorithms, an adjacency-matrix representation, betweenness-centrality
//! estimation (Brandes), and a simple spring embedder.

#![allow(dead_code)]

// ─── Node and Edge Layout ─────────────────────────────────────────────────────

/// Visual layout data for a single graph node.
#[derive(Debug, Clone)]
pub struct NodeLayout {
    /// Unique identifier for this node.
    pub id: usize,
    /// 2-D position of the node.
    pub position: [f64; 2],
    /// Radius used for rendering the node.
    pub radius: f64,
    /// RGBA colour of the node.
    pub color: [f32; 4],
}

/// Visual layout data for a single graph edge.
#[derive(Debug, Clone)]
pub struct EdgeLayout {
    /// Index of the source node.
    pub from: usize,
    /// Index of the target node.
    pub to: usize,
    /// Edge weight (used by layout algorithms).
    pub weight: f64,
    /// RGBA colour of the edge.
    pub color: [f32; 4],
}

// ─── NetworkGraph ─────────────────────────────────────────────────────────────

/// A graph with visual layout information for nodes and edges.
#[derive(Debug, Clone, Default)]
pub struct NetworkGraph {
    /// All nodes in the graph.
    pub nodes: Vec<NodeLayout>,
    /// All edges in the graph.
    pub edges: Vec<EdgeLayout>,
}

impl NetworkGraph {
    /// Create an empty graph.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a node and return its index in `nodes`.
    pub fn add_node(&mut self, id: usize, pos: [f64; 2], radius: f64) -> usize {
        let idx = self.nodes.len();
        self.nodes.push(NodeLayout {
            id,
            position: pos,
            radius,
            color: [1.0, 1.0, 1.0, 1.0],
        });
        idx
    }

    /// Add an edge between nodes identified by their indices in `nodes`.
    pub fn add_edge(&mut self, from: usize, to: usize, weight: f64) {
        self.edges.push(EdgeLayout {
            from,
            to,
            weight,
            color: [0.8, 0.8, 0.8, 1.0],
        });
    }

    /// Number of nodes in the graph.
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }

    /// Number of edges in the graph.
    pub fn edge_count(&self) -> usize {
        self.edges.len()
    }
}

// ─── Fruchterman–Reingold Force-directed Layout ───────────────────────────────

/// Apply the Fruchterman–Reingold force-directed layout algorithm in-place.
///
/// - `iters` — number of iterations.
/// - `k` — optimal distance between connected nodes (ideally ≈ √(area / n)).
/// - `repulsion` — scaling factor for the repulsive force between all node pairs.
#[allow(clippy::too_many_arguments)]
pub fn force_directed_layout(graph: &mut NetworkGraph, iters: usize, k: f64, repulsion: f64) {
    let n = graph.nodes.len();
    if n < 2 {
        return;
    }

    for _iter in 0..iters {
        let mut disp: Vec<[f64; 2]> = vec![[0.0; 2]; n];

        // Repulsive forces (all pairs).
        for i in 0..n {
            for j in 0..n {
                if i == j {
                    continue;
                }
                let dx = graph.nodes[i].position[0] - graph.nodes[j].position[0];
                let dy = graph.nodes[i].position[1] - graph.nodes[j].position[1];
                let dist = (dx * dx + dy * dy).sqrt().max(1e-6);
                let force = repulsion * k * k / dist;
                disp[i][0] += force * dx / dist;
                disp[i][1] += force * dy / dist;
            }
        }

        // Attractive forces (edges).
        let edges: Vec<(usize, usize, f64)> = graph
            .edges
            .iter()
            .map(|e| (e.from, e.to, e.weight))
            .collect();

        for (from, to, _w) in &edges {
            let f = *from;
            let t = *to;
            if f >= n || t >= n {
                continue;
            }
            let dx = graph.nodes[f].position[0] - graph.nodes[t].position[0];
            let dy = graph.nodes[f].position[1] - graph.nodes[t].position[1];
            let dist = (dx * dx + dy * dy).sqrt().max(1e-6);
            let force = dist * dist / k;
            let fx = force * dx / dist;
            let fy = force * dy / dist;
            disp[f][0] -= fx;
            disp[f][1] -= fy;
            disp[t][0] += fx;
            disp[t][1] += fy;
        }

        // Apply displacements (capped to avoid explosions).
        let temperature = k / (_iter as f64 + 1.0).max(1.0);
        for i in 0..n {
            let dlen = (disp[i][0] * disp[i][0] + disp[i][1] * disp[i][1])
                .sqrt()
                .max(1e-12);
            let cap = dlen.min(temperature);
            graph.nodes[i].position[0] += cap * disp[i][0] / dlen;
            graph.nodes[i].position[1] += cap * disp[i][1] / dlen;
        }
    }
}

// ─── Circular Layout ──────────────────────────────────────────────────────────

/// Place all nodes evenly on a circle of the given `radius`.
pub fn circular_layout(graph: &mut NetworkGraph, radius: f64) {
    let n = graph.nodes.len();
    if n == 0 {
        return;
    }
    let step = 2.0 * std::f64::consts::PI / n as f64;
    for (i, node) in graph.nodes.iter_mut().enumerate() {
        let angle = i as f64 * step;
        node.position = [radius * angle.cos(), radius * angle.sin()];
    }
}

// ─── Hierarchical Layout ──────────────────────────────────────────────────────

/// Arrange nodes in a layered hierarchical layout.
///
/// `levels` is a slice of layers; each layer is a list of node indices (into
/// `graph.nodes`).  Nodes in level `l` are placed at y = l, evenly spaced in x.
pub fn hierarchical_layout(graph: &mut NetworkGraph, levels: &[Vec<usize>]) {
    for (level_idx, level) in levels.iter().enumerate() {
        let y = level_idx as f64;
        let count = level.len();
        if count == 0 {
            continue;
        }
        for (col, &node_idx) in level.iter().enumerate() {
            if node_idx >= graph.nodes.len() {
                continue;
            }
            let x = if count == 1 {
                0.0
            } else {
                col as f64 / (count as f64 - 1.0) - 0.5
            };
            graph.nodes[node_idx].position = [x, y];
        }
    }
}

// ─── Edge Crossings ───────────────────────────────────────────────────────────

/// Count the number of edge crossings in the current layout.
///
/// Two edges cross if their corresponding line segments intersect (and they do
/// not share a node).
pub fn compute_edge_crossings(graph: &NetworkGraph) -> usize {
    let mut count = 0;
    let edges = &graph.edges;
    let nodes = &graph.nodes;
    for i in 0..edges.len() {
        for j in (i + 1)..edges.len() {
            let ei = &edges[i];
            let ej = &edges[j];
            // Skip edges sharing a node.
            if ei.from == ej.from || ei.from == ej.to || ei.to == ej.from || ei.to == ej.to {
                continue;
            }
            if ei.from >= nodes.len()
                || ei.to >= nodes.len()
                || ej.from >= nodes.len()
                || ej.to >= nodes.len()
            {
                continue;
            }
            let p1 = nodes[ei.from].position;
            let p2 = nodes[ei.to].position;
            let p3 = nodes[ej.from].position;
            let p4 = nodes[ej.to].position;
            if segments_intersect(p1, p2, p3, p4) {
                count += 1;
            }
        }
    }
    count
}

fn segments_intersect(p1: [f64; 2], p2: [f64; 2], p3: [f64; 2], p4: [f64; 2]) -> bool {
    let cross2 = |a: [f64; 2], b: [f64; 2]| a[0] * b[1] - a[1] * b[0];
    let sub2 = |a: [f64; 2], b: [f64; 2]| [a[0] - b[0], a[1] - b[1]];
    let r = sub2(p2, p1);
    let s = sub2(p4, p3);
    let denom = cross2(r, s);
    if denom.abs() < 1e-12 {
        return false; // parallel
    }
    let t = cross2(sub2(p3, p1), s) / denom;
    let u = cross2(sub2(p3, p1), r) / denom;
    t > 0.0 && t < 1.0 && u > 0.0 && u < 1.0
}

// ─── Adjacency Matrix ─────────────────────────────────────────────────────────

/// Dense adjacency matrix for a graph.
#[derive(Debug, Clone)]
pub struct AdjacencyMatrix {
    /// Number of nodes.
    pub n: usize,
    /// Flat n×n row-major data.
    pub data: Vec<f64>,
}

impl AdjacencyMatrix {
    /// Construct a zero matrix for `n` nodes.
    pub fn new(n: usize) -> Self {
        Self {
            n,
            data: vec![0.0; n * n],
        }
    }

    /// Build an adjacency matrix from a `NetworkGraph`.
    ///
    /// For each edge `(from, to, weight)` sets `A[from][to] = weight` and
    /// `A[to][from] = weight` (treats graph as undirected).
    pub fn from_graph(g: &NetworkGraph) -> Self {
        let n = g.nodes.len();
        let mut mat = Self::new(n);
        for edge in &g.edges {
            let (f, t) = (edge.from, edge.to);
            if f < n && t < n {
                mat.data[f * n + t] = edge.weight;
                mat.data[t * n + f] = edge.weight;
            }
        }
        mat
    }

    /// Get the weight of edge `(i, j)`.
    pub fn get(&self, i: usize, j: usize) -> f64 {
        self.data[i * self.n + j]
    }

    /// Degree of node `i` (sum of weights of incident edges).
    pub fn degree(&self, i: usize) -> f64 {
        (0..self.n).map(|j| self.data[i * self.n + j]).sum()
    }
}

// ─── Betweenness Centrality ───────────────────────────────────────────────────

/// Compute approximate betweenness centrality using Brandes' algorithm on the
/// unweighted graph.
///
/// Returns a vector of length `n`, where entry `i` is the centrality of node `i`.
pub fn compute_betweenness_centrality(graph: &NetworkGraph) -> Vec<f64> {
    let n = graph.nodes.len();
    if n == 0 {
        return Vec::new();
    }
    let mut centrality = vec![0.0f64; n];

    // Build adjacency list (unweighted, undirected).
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
    for edge in &graph.edges {
        let (f, t) = (edge.from, edge.to);
        if f < n && t < n {
            adj[f].push(t);
            adj[t].push(f);
        }
    }

    for src in 0..n {
        // BFS from src.
        let mut sigma = vec![0.0f64; n];
        let mut dist = vec![-1i64; n];
        let mut stack: Vec<usize> = Vec::new();
        let mut pred: Vec<Vec<usize>> = vec![Vec::new(); n];

        sigma[src] = 1.0;
        dist[src] = 0;
        let mut queue = std::collections::VecDeque::new();
        queue.push_back(src);

        while let Some(v) = queue.pop_front() {
            stack.push(v);
            for &w in &adj[v] {
                if dist[w] < 0 {
                    dist[w] = dist[v] + 1;
                    queue.push_back(w);
                }
                if dist[w] == dist[v] + 1 {
                    sigma[w] += sigma[v];
                    pred[w].push(v);
                }
            }
        }

        // Accumulate dependencies.
        let mut delta = vec![0.0f64; n];
        while let Some(w) = stack.pop() {
            for &v in &pred[w] {
                delta[v] += (sigma[v] / sigma[w]) * (1.0 + delta[w]);
            }
            if w != src {
                centrality[w] += delta[w];
            }
        }
    }

    // Normalise by dividing by 2 (undirected graph: each pair counted twice).
    for c in &mut centrality {
        *c /= 2.0;
    }
    centrality
}

// ─── Spring Embedder ─────────────────────────────────────────────────────────

/// Simple spring embedder: iterate `iters` times applying spring forces from
/// the adjacency matrix.
///
/// Nodes connected by edges are attracted; all node pairs repel.
pub fn spring_embedder(positions: &mut [[f64; 2]], adj: &AdjacencyMatrix, iters: usize) {
    let n = positions.len();
    if n < 2 {
        return;
    }
    let k = 1.0_f64; // spring rest length

    for _ in 0..iters {
        let mut forces: Vec<[f64; 2]> = vec![[0.0; 2]; n];

        for i in 0..n {
            for j in 0..n {
                if i == j {
                    continue;
                }
                let dx = positions[i][0] - positions[j][0];
                let dy = positions[i][1] - positions[j][1];
                let dist = (dx * dx + dy * dy).sqrt().max(1e-6);

                // Repulsion.
                let rep = k * k / dist;
                forces[i][0] += rep * dx / dist;
                forces[i][1] += rep * dy / dist;

                // Attraction via spring (only if connected).
                if adj.n > i && adj.n > j {
                    let w = adj.get(i, j);
                    if w > 0.0 {
                        let attr = dist / k;
                        forces[i][0] -= attr * dx / dist;
                        forces[i][1] -= attr * dy / dist;
                    }
                }
            }
        }

        let lr = 0.01;
        for i in 0..n {
            positions[i][0] += lr * forces[i][0];
            positions[i][1] += lr * forces[i][1];
        }
    }
}

// ─── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // --- NetworkGraph basics ---

    #[test]
    fn test_add_node_increments_count() {
        let mut g = NetworkGraph::new();
        assert_eq!(g.node_count(), 0);
        g.add_node(0, [0.0, 0.0], 1.0);
        assert_eq!(g.node_count(), 1);
        g.add_node(1, [1.0, 0.0], 1.0);
        assert_eq!(g.node_count(), 2);
    }

    #[test]
    fn test_add_edge_increments_count() {
        let mut g = NetworkGraph::new();
        g.add_node(0, [0.0, 0.0], 1.0);
        g.add_node(1, [1.0, 0.0], 1.0);
        g.add_edge(0, 1, 1.0);
        assert_eq!(g.edge_count(), 1);
    }

    #[test]
    fn test_network_graph_default_empty() {
        let g = NetworkGraph::default();
        assert_eq!(g.node_count(), 0);
        assert_eq!(g.edge_count(), 0);
    }

    #[test]
    fn test_node_id_stored() {
        let mut g = NetworkGraph::new();
        g.add_node(42, [0.0, 0.0], 0.5);
        assert_eq!(g.nodes[0].id, 42);
    }

    #[test]
    fn test_edge_stores_from_to_weight() {
        let mut g = NetworkGraph::new();
        g.add_node(0, [0.0, 0.0], 1.0);
        g.add_node(1, [1.0, 0.0], 1.0);
        g.add_edge(0, 1, 3.5);
        assert_eq!(g.edges[0].from, 0);
        assert_eq!(g.edges[0].to, 1);
        assert!((g.edges[0].weight - 3.5).abs() < 1e-12);
    }

    // --- circular_layout ---

    #[test]
    fn test_circular_layout_equidistant_from_center() {
        let mut g = NetworkGraph::new();
        for i in 0..8 {
            g.add_node(i, [0.0, 0.0], 1.0);
        }
        circular_layout(&mut g, 5.0);
        for node in &g.nodes {
            let r = (node.position[0].powi(2) + node.position[1].powi(2)).sqrt();
            assert!((r - 5.0).abs() < 1e-10, "radius = {}", r);
        }
    }

    #[test]
    fn test_circular_layout_single_node() {
        let mut g = NetworkGraph::new();
        g.add_node(0, [0.0, 0.0], 1.0);
        circular_layout(&mut g, 3.0);
        let r = (g.nodes[0].position[0].powi(2) + g.nodes[0].position[1].powi(2)).sqrt();
        assert!((r - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_circular_layout_empty_graph() {
        let mut g = NetworkGraph::new();
        circular_layout(&mut g, 1.0); // should not panic
    }

    #[test]
    fn test_circular_layout_positions_distinct() {
        let n = 6;
        let mut g = NetworkGraph::new();
        for i in 0..n {
            g.add_node(i, [0.0, 0.0], 1.0);
        }
        circular_layout(&mut g, 2.0);
        // All positions should be distinct.
        for i in 0..n {
            for j in (i + 1)..n {
                let dx = g.nodes[i].position[0] - g.nodes[j].position[0];
                let dy = g.nodes[i].position[1] - g.nodes[j].position[1];
                let dist = (dx * dx + dy * dy).sqrt();
                assert!(dist > 0.01, "nodes {} and {} are at same position", i, j);
            }
        }
    }

    // --- hierarchical_layout ---

    #[test]
    fn test_hierarchical_layout_levels() {
        let mut g = NetworkGraph::new();
        for i in 0..4 {
            g.add_node(i, [0.0, 0.0], 1.0);
        }
        let levels = vec![vec![0, 1], vec![2, 3]];
        hierarchical_layout(&mut g, &levels);
        // Level 0 nodes at y = 0, level 1 at y = 1.
        assert!((g.nodes[0].position[1] - 0.0).abs() < 1e-10);
        assert!((g.nodes[2].position[1] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_hierarchical_layout_single_level() {
        let mut g = NetworkGraph::new();
        g.add_node(0, [0.0, 0.0], 1.0);
        hierarchical_layout(&mut g, &[vec![0]]);
        // Single node: x = 0.
        assert!((g.nodes[0].position[0]).abs() < 1e-10);
    }

    // --- AdjacencyMatrix ---

    #[test]
    fn test_adjacency_matrix_from_graph_symmetric() {
        let mut g = NetworkGraph::new();
        for i in 0..3 {
            g.add_node(i, [0.0, 0.0], 1.0);
        }
        g.add_edge(0, 1, 1.0);
        g.add_edge(1, 2, 2.0);
        let adj = AdjacencyMatrix::from_graph(&g);
        assert!((adj.get(0, 1) - 1.0).abs() < 1e-12);
        assert!((adj.get(1, 0) - 1.0).abs() < 1e-12); // symmetric
        assert!((adj.get(1, 2) - 2.0).abs() < 1e-12);
        assert!((adj.get(2, 1) - 2.0).abs() < 1e-12);
        assert!((adj.get(0, 2)).abs() < 1e-12); // not connected
    }

    #[test]
    fn test_adjacency_matrix_degree() {
        let mut g = NetworkGraph::new();
        for i in 0..3 {
            g.add_node(i, [0.0, 0.0], 1.0);
        }
        g.add_edge(0, 1, 1.0);
        g.add_edge(0, 2, 1.0);
        let adj = AdjacencyMatrix::from_graph(&g);
        assert!(
            (adj.degree(0) - 2.0).abs() < 1e-12,
            "degree(0) = {}",
            adj.degree(0)
        );
    }

    #[test]
    fn test_adjacency_matrix_empty_graph() {
        let g = NetworkGraph::new();
        let adj = AdjacencyMatrix::from_graph(&g);
        assert_eq!(adj.n, 0);
    }

    // --- force_directed_layout ---

    #[test]
    fn test_force_directed_layout_positions_change() {
        let mut g = NetworkGraph::new();
        for i in 0..5 {
            g.add_node(i, [i as f64, 0.0], 1.0);
        }
        g.add_edge(0, 1, 1.0);
        g.add_edge(1, 2, 1.0);
        g.add_edge(2, 3, 1.0);
        g.add_edge(3, 4, 1.0);
        let orig: Vec<[f64; 2]> = g.nodes.iter().map(|n| n.position).collect();
        force_directed_layout(&mut g, 10, 2.0, 1.0);
        let changed = g.nodes.iter().zip(orig.iter()).any(|(n, o)| {
            (n.position[0] - o[0]).abs() > 1e-10 || (n.position[1] - o[1]).abs() > 1e-10
        });
        assert!(
            changed,
            "positions should change after force-directed layout"
        );
    }

    #[test]
    fn test_force_directed_layout_single_node_no_crash() {
        let mut g = NetworkGraph::new();
        g.add_node(0, [0.0, 0.0], 1.0);
        force_directed_layout(&mut g, 5, 1.0, 1.0);
    }

    // --- betweenness_centrality ---

    #[test]
    fn test_betweenness_centrality_non_negative() {
        let mut g = NetworkGraph::new();
        for i in 0..4 {
            g.add_node(i, [i as f64, 0.0], 1.0);
        }
        g.add_edge(0, 1, 1.0);
        g.add_edge(1, 2, 1.0);
        g.add_edge(2, 3, 1.0);
        let bc = compute_betweenness_centrality(&g);
        assert_eq!(bc.len(), 4);
        for &c in &bc {
            assert!(c >= 0.0, "centrality {} is negative", c);
        }
    }

    #[test]
    fn test_betweenness_centrality_path_graph_middle_highest() {
        // In a path 0-1-2-3, nodes 1 and 2 should have higher betweenness.
        let mut g = NetworkGraph::new();
        for i in 0..4 {
            g.add_node(i, [i as f64, 0.0], 1.0);
        }
        g.add_edge(0, 1, 1.0);
        g.add_edge(1, 2, 1.0);
        g.add_edge(2, 3, 1.0);
        let bc = compute_betweenness_centrality(&g);
        assert!(
            bc[1] > bc[0],
            "node 1 should have higher BC than endpoint 0"
        );
        assert!(
            bc[2] > bc[3],
            "node 2 should have higher BC than endpoint 3"
        );
    }

    #[test]
    fn test_betweenness_centrality_empty_graph() {
        let g = NetworkGraph::new();
        let bc = compute_betweenness_centrality(&g);
        assert!(bc.is_empty());
    }

    // --- spring_embedder ---

    #[test]
    fn test_spring_embedder_positions_change() {
        let mut positions = vec![[0.0f64, 0.0], [5.0, 0.0], [10.0, 0.0]];
        let mut g = NetworkGraph::new();
        for i in 0..3 {
            g.add_node(i, positions[i], 1.0);
        }
        g.add_edge(0, 1, 1.0);
        g.add_edge(1, 2, 1.0);
        let adj = AdjacencyMatrix::from_graph(&g);
        let orig = positions.clone();
        spring_embedder(&mut positions, &adj, 20);
        let changed = positions
            .iter()
            .zip(orig.iter())
            .any(|(p, o)| (p[0] - o[0]).abs() > 1e-10 || (p[1] - o[1]).abs() > 1e-10);
        assert!(changed, "spring embedder should move nodes");
    }

    #[test]
    fn test_spring_embedder_single_node_no_crash() {
        let mut positions = vec![[1.0f64, 2.0]];
        let adj = AdjacencyMatrix::new(1);
        spring_embedder(&mut positions, &adj, 5);
    }

    // --- compute_edge_crossings ---

    #[test]
    fn test_edge_crossings_crossing_x() {
        // Two edges forming an X should cross.
        let mut g = NetworkGraph::new();
        g.add_node(0, [-1.0, -1.0], 1.0);
        g.add_node(1, [1.0, 1.0], 1.0);
        g.add_node(2, [-1.0, 1.0], 1.0);
        g.add_node(3, [1.0, -1.0], 1.0);
        g.add_edge(0, 1, 1.0); // diagonal
        g.add_edge(2, 3, 1.0); // other diagonal
        let crossings = compute_edge_crossings(&g);
        assert_eq!(crossings, 1, "expected 1 crossing, got {}", crossings);
    }

    #[test]
    fn test_edge_crossings_no_crossing() {
        // Parallel horizontal edges should not cross.
        let mut g = NetworkGraph::new();
        g.add_node(0, [0.0, 0.0], 1.0);
        g.add_node(1, [1.0, 0.0], 1.0);
        g.add_node(2, [0.0, 1.0], 1.0);
        g.add_node(3, [1.0, 1.0], 1.0);
        g.add_edge(0, 1, 1.0);
        g.add_edge(2, 3, 1.0);
        let crossings = compute_edge_crossings(&g);
        assert_eq!(crossings, 0);
    }
}
