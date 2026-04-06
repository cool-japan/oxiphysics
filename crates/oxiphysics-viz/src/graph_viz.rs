#![allow(clippy::needless_range_loop, clippy::ptr_arg)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Graph and network visualization.
//!
//! Provides layout algorithms (force-directed, circular, hierarchical),
//! graph data structures, BFS/DFS traversal, Dijkstra shortest paths,
//! adjacency matrix, and basic graph statistics.

use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashMap, VecDeque};

// ---------------------------------------------------------------------------
// Core graph data types
// ---------------------------------------------------------------------------

/// A node in a visualization graph.
#[derive(Debug, Clone)]
pub struct GraphNode {
    /// Unique identifier.
    pub id: usize,
    /// 2-D position used for layout.
    pub position: [f64; 2],
    /// Human-readable label.
    pub label: String,
    /// RGB color (each component in `[0, 1]`).
    pub color: [f32; 3],
    /// Display size (radius or half-extent).
    pub size: f32,
}

impl GraphNode {
    /// Construct a new node with the given `id` placed at the origin.
    pub fn new(id: usize, label: impl Into<String>) -> Self {
        GraphNode {
            id,
            position: [0.0, 0.0],
            label: label.into(),
            color: [0.5, 0.5, 1.0],
            size: 1.0,
        }
    }
}

/// A directed or undirected edge between two nodes.
#[derive(Debug, Clone)]
pub struct GraphEdge {
    /// Source node id.
    pub src: usize,
    /// Destination node id.
    pub dst: usize,
    /// Edge weight (used by layout forces and Dijkstra).
    pub weight: f64,
    /// If `true`, the edge is directed from `src` to `dst`.
    pub directed: bool,
}

impl GraphEdge {
    /// Construct an undirected edge with weight 1.
    pub fn undirected(src: usize, dst: usize) -> Self {
        GraphEdge {
            src,
            dst,
            weight: 1.0,
            directed: false,
        }
    }

    /// Construct a directed edge with the given weight.
    pub fn directed(src: usize, dst: usize, weight: f64) -> Self {
        GraphEdge {
            src,
            dst,
            weight,
            directed: true,
        }
    }
}

/// A complete graph layout: positioned nodes plus edges.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct GraphLayout {
    /// Positioned nodes.
    pub nodes: Vec<GraphNode>,
    /// Edges.
    pub edges: Vec<GraphEdge>,
}

impl GraphLayout {
    /// Construct a layout from existing nodes and edges.
    pub fn new(nodes: Vec<GraphNode>, edges: Vec<GraphEdge>) -> Self {
        GraphLayout { nodes, edges }
    }

    /// Axis-aligned bounding box of all node positions: `([xmin,ymin], [xmax,ymax])`.
    ///
    /// Returns `None` if there are no nodes.
    pub fn bounding_box(&self) -> Option<([f64; 2], [f64; 2])> {
        if self.nodes.is_empty() {
            return None;
        }
        let mut mn = [f64::INFINITY; 2];
        let mut mx = [f64::NEG_INFINITY; 2];
        for n in &self.nodes {
            for k in 0..2 {
                if n.position[k] < mn[k] {
                    mn[k] = n.position[k];
                }
                if n.position[k] > mx[k] {
                    mx[k] = n.position[k];
                }
            }
        }
        Some((mn, mx))
    }
}

// ---------------------------------------------------------------------------
// Circular layout
// ---------------------------------------------------------------------------

/// Place `n` nodes evenly on a circle of the given `radius`.
///
/// Returns a `Vec` of `[x, y]` positions, starting at angle 0 and going
/// counter-clockwise.
pub fn circular_layout(n: usize, radius: f64) -> Vec<[f64; 2]> {
    if n == 0 {
        return vec![];
    }
    let step = 2.0 * std::f64::consts::PI / n as f64;
    (0..n)
        .map(|i| {
            let a = i as f64 * step;
            [radius * a.cos(), radius * a.sin()]
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Force-directed layout (Fruchterman-Reingold)
// ---------------------------------------------------------------------------

#[allow(dead_code)]
fn dist2(a: [f64; 2], b: [f64; 2]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    (dx * dx + dy * dy).sqrt().max(1e-10)
}

/// Run Fruchterman-Reingold spring-electrical force-directed layout.
///
/// - `nodes` — nodes whose `position` field is updated in-place.
/// - `edges` — connectivity (only `src`/`dst` used, not weights).
/// - `iters` — number of simulation iterations.
///
/// The algorithm uses a cooling schedule: the temperature (maximum
/// displacement per step) starts at the graph diameter and decays linearly.
pub fn force_directed_layout(nodes: &mut Vec<GraphNode>, edges: &[GraphEdge], iters: usize) {
    let n = nodes.len();
    if n == 0 {
        return;
    }

    // Ideal spring length based on area
    let area = (n as f64) * 100.0;
    let k = (area / n as f64).sqrt();

    let mut temp = (area).sqrt();
    let cool = temp / iters.max(1) as f64;

    // Displacement accumulator
    let mut disp: Vec<[f64; 2]> = vec![[0.0, 0.0]; n];

    for _iter in 0..iters {
        // Zero displacements
        for d in disp.iter_mut() {
            *d = [0.0, 0.0];
        }

        // Repulsive forces between all pairs
        for i in 0..n {
            for j in (i + 1)..n {
                let pi = nodes[i].position;
                let pj = nodes[j].position;
                // Add a deterministic tiny perturbation so that exactly
                // co-located nodes still repel each other in a defined direction.
                let mut dx = pi[0] - pj[0];
                let mut dy = pi[1] - pj[1];
                if dx == 0.0 && dy == 0.0 {
                    // Use index-based pseudo-random offset
                    let angle =
                        (i as f64 * 2.399963 + j as f64 * 1.618033) % (2.0 * std::f64::consts::PI);
                    dx = 1e-6 * angle.cos();
                    dy = 1e-6 * angle.sin();
                }
                let d = (dx * dx + dy * dy).sqrt().max(1e-10);
                let fr = k * k / d; // repulsive magnitude
                let fx = (dx / d) * fr;
                let fy = (dy / d) * fr;
                disp[i][0] += fx;
                disp[i][1] += fy;
                disp[j][0] -= fx;
                disp[j][1] -= fy;
            }
        }

        // Attractive forces along edges
        for edge in edges {
            let i = edge.src.min(n - 1);
            let j = edge.dst.min(n - 1);
            if i == j {
                continue;
            }
            let pi = nodes[i].position;
            let pj = nodes[j].position;
            let dx = pi[0] - pj[0];
            let dy = pi[1] - pj[1];
            let d = (dx * dx + dy * dy).sqrt().max(1e-10);
            let fa = d * d / k; // attractive magnitude
            let fx = (dx / d) * fa;
            let fy = (dy / d) * fa;
            disp[i][0] -= fx;
            disp[i][1] -= fy;
            disp[j][0] += fx;
            disp[j][1] += fy;
        }

        // Apply displacement (clamped to temperature)
        for i in 0..n {
            let dx = disp[i][0];
            let dy = disp[i][1];
            let len = (dx * dx + dy * dy).sqrt().max(1e-10);
            let scale = len.min(temp) / len;
            nodes[i].position[0] += dx * scale;
            nodes[i].position[1] += dy * scale;
        }

        // Cool
        temp = (temp - cool).max(0.0);
    }
}

// ---------------------------------------------------------------------------
// Hierarchical layout (Sugiyama-style: topological layers)
// ---------------------------------------------------------------------------

/// Assign nodes to horizontal layers via topological sort (Sugiyama style).
///
/// Each node's y-coordinate is set to its layer index; nodes within the same
/// layer are spaced evenly along the x-axis.  Cycles are broken by ignoring
/// back-edges (edges that would create a cycle are skipped).
pub fn hierarchical_layout(nodes: &mut Vec<GraphNode>, edges: &[GraphEdge]) {
    let n = nodes.len();
    if n == 0 {
        return;
    }

    // Build adjacency list considering only directed edges for Kahn's algorithm
    let mut in_degree = vec![0usize; n];
    let mut adj: Vec<Vec<usize>> = vec![vec![]; n];
    for e in edges {
        let src = e.src.min(n - 1);
        let dst = e.dst.min(n - 1);
        if src != dst && e.directed {
            adj[src].push(dst);
            in_degree[dst] += 1;
        }
    }

    // Kahn's BFS-based topological sort → layer assignment
    let mut layer = vec![0usize; n];
    let mut queue: VecDeque<usize> = VecDeque::new();
    for i in 0..n {
        if in_degree[i] == 0 {
            queue.push_back(i);
        }
    }

    while let Some(u) = queue.pop_front() {
        for &v in &adj[u] {
            layer[v] = layer[v].max(layer[u] + 1);
            in_degree[v] -= 1;
            if in_degree[v] == 0 {
                queue.push_back(v);
            }
        }
    }

    // Nodes not reached by the topological sort (cycles) get layer 0
    // Group nodes by layer to compute x-offsets
    let max_layer = *layer.iter().max().unwrap_or(&0);
    let mut layer_counts = vec![0usize; max_layer + 1];
    let mut layer_idx = vec![0usize; n];
    // First pass: count nodes per layer
    for i in 0..n {
        layer_idx[i] = layer_counts[layer[i]];
        layer_counts[layer[i]] += 1;
    }

    let spacing_x = 3.0_f64;
    let spacing_y = 3.0_f64;

    for i in 0..n {
        let lyr = layer[i] as f64;
        let cnt = layer_counts[layer[i]] as f64;
        let idx = layer_idx[i] as f64;
        // Centre each layer horizontally
        nodes[i].position[0] = (idx - (cnt - 1.0) / 2.0) * spacing_x;
        nodes[i].position[1] = lyr * spacing_y;
    }
}

// ---------------------------------------------------------------------------
// AdjacencyMatrix
// ---------------------------------------------------------------------------

/// Sparse adjacency representation built from an edge list.
///
/// Supports degree queries, connectivity check, and shortest-path lookups.
#[derive(Debug, Clone)]
pub struct AdjacencyMatrix {
    /// Number of nodes.
    pub n: usize,
    /// Adjacency list: `adj[u]` contains `(v, weight)` pairs.
    pub adj: Vec<Vec<(usize, f64)>>,
    /// Degree of each node (sum of in- and out-edges for directed graphs,
    /// or twice the number of incident edges for undirected).
    pub degree: Vec<usize>,
}

impl AdjacencyMatrix {
    /// Build an adjacency matrix from `n` nodes and a slice of edges.
    pub fn from_edges(n: usize, edges: &[GraphEdge]) -> Self {
        let mut adj = vec![vec![]; n];
        let mut degree = vec![0usize; n];
        for e in edges {
            let src = e.src.min(n.saturating_sub(1));
            let dst = e.dst.min(n.saturating_sub(1));
            if n == 0 {
                continue;
            }
            adj[src].push((dst, e.weight));
            degree[src] += 1;
            if !e.directed {
                adj[dst].push((src, e.weight));
                degree[dst] += 1;
            } else {
                degree[dst] += 1;
            }
        }
        AdjacencyMatrix { n, adj, degree }
    }

    /// Return `true` if the graph is (weakly) connected.
    ///
    /// Uses BFS from node 0 over the undirected interpretation of edges.
    pub fn is_connected(&self) -> bool {
        if self.n == 0 {
            return true;
        }
        let mut visited = vec![false; self.n];
        let mut queue = VecDeque::new();
        visited[0] = true;
        queue.push_back(0usize);
        while let Some(u) = queue.pop_front() {
            for &(v, _) in &self.adj[u] {
                if !visited[v] {
                    visited[v] = true;
                    queue.push_back(v);
                }
            }
        }
        visited.iter().all(|&v| v)
    }

    /// Number of edges (directed count; undirected edges appear twice in `adj`).
    #[allow(dead_code)]
    pub fn edge_count(&self) -> usize {
        self.adj.iter().map(|a| a.len()).sum::<usize>()
    }
}

// ---------------------------------------------------------------------------
// BfsOrder
// ---------------------------------------------------------------------------

/// Result of a BFS traversal from a source node.
#[derive(Debug, Clone)]
pub struct BfsOrder {
    /// Nodes in BFS visit order.
    pub order: Vec<usize>,
    /// BFS distance from source to each node (`usize::MAX` if unreachable).
    pub distance: Vec<usize>,
    /// Parent of each node in the BFS tree (`None` for source or unreachable).
    pub parent: Vec<Option<usize>>,
}

impl BfsOrder {
    /// Run BFS from `src` on the given adjacency matrix.
    pub fn from(adj: &AdjacencyMatrix, src: usize) -> Self {
        let n = adj.n;
        let mut order = Vec::with_capacity(n);
        let mut distance = vec![usize::MAX; n];
        let mut parent = vec![None; n];
        if src >= n {
            return BfsOrder {
                order,
                distance,
                parent,
            };
        }
        distance[src] = 0;
        let mut queue = VecDeque::new();
        queue.push_back(src);
        while let Some(u) = queue.pop_front() {
            order.push(u);
            for &(v, _) in &adj.adj[u] {
                if distance[v] == usize::MAX {
                    distance[v] = distance[u] + 1;
                    parent[v] = Some(u);
                    queue.push_back(v);
                }
            }
        }
        BfsOrder {
            order,
            distance,
            parent,
        }
    }
}

// ---------------------------------------------------------------------------
// DfsOrder
// ---------------------------------------------------------------------------

/// Result of a DFS traversal from a source node.
#[derive(Debug, Clone)]
pub struct DfsOrder {
    /// Nodes in DFS visit order.
    pub order: Vec<usize>,
    /// Whether each node was visited.
    pub visited: Vec<bool>,
}

impl DfsOrder {
    /// Run iterative DFS from `src`.
    pub fn from(adj: &AdjacencyMatrix, src: usize) -> Self {
        let n = adj.n;
        let mut order = Vec::with_capacity(n);
        let mut visited = vec![false; n];
        if src >= n {
            return DfsOrder { order, visited };
        }
        let mut stack = vec![src];
        while let Some(u) = stack.pop() {
            if visited[u] {
                continue;
            }
            visited[u] = true;
            order.push(u);
            // Push neighbors in reverse so that lower-id neighbours are visited first
            for &(v, _) in adj.adj[u].iter().rev() {
                if !visited[v] {
                    stack.push(v);
                }
            }
        }
        DfsOrder { order, visited }
    }
}

// ---------------------------------------------------------------------------
// Dijkstra shortest path
// ---------------------------------------------------------------------------

// Wrapper for BinaryHeap (min-heap via Reverse)
#[derive(Clone, PartialEq)]
struct State {
    cost: f64,
    node: usize,
}

impl Eq for State {}

impl Ord for State {
    fn cmp(&self, other: &Self) -> Ordering {
        other
            .cost
            .partial_cmp(&self.cost)
            .unwrap_or(Ordering::Equal)
    }
}

impl PartialOrd for State {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// Dijkstra's single-source shortest path from `src` to `dst`.
///
/// Returns `Some(path)` where `path` is the sequence of node ids from
/// `src` to `dst` (inclusive), or `None` if no path exists.
///
/// All edge weights must be non-negative.
pub fn shortest_path_dijkstra(adj: &AdjacencyMatrix, src: usize, dst: usize) -> Option<Vec<usize>> {
    let n = adj.n;
    if src >= n || dst >= n {
        return None;
    }
    let mut dist = vec![f64::INFINITY; n];
    let mut prev: Vec<Option<usize>> = vec![None; n];
    dist[src] = 0.0;

    let mut heap = BinaryHeap::new();
    heap.push(State {
        cost: 0.0,
        node: src,
    });

    while let Some(State { cost, node: u }) = heap.pop() {
        if u == dst {
            // Reconstruct path
            let mut path = vec![];
            let mut cur = dst;
            loop {
                path.push(cur);
                match prev[cur] {
                    Some(p) => cur = p,
                    None => break,
                }
            }
            path.reverse();
            return Some(path);
        }
        if cost > dist[u] {
            continue;
        }
        for &(v, w) in &adj.adj[u] {
            let next_cost = dist[u] + w;
            if next_cost < dist[v] {
                dist[v] = next_cost;
                prev[v] = Some(u);
                heap.push(State {
                    cost: next_cost,
                    node: v,
                });
            }
        }
    }
    None
}

// ---------------------------------------------------------------------------
// GraphStats
// ---------------------------------------------------------------------------

/// Summary statistics of a graph.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct GraphStats {
    /// Number of nodes.
    pub node_count: usize,
    /// Number of edges (directed count).
    pub edge_count: usize,
    /// Minimum degree.
    pub min_degree: usize,
    /// Maximum degree.
    pub max_degree: usize,
    /// Average degree.
    pub avg_degree: f64,
    /// Histogram: maps degree → count.
    pub degree_distribution: HashMap<usize, usize>,
    /// Global clustering coefficient (fraction of closed triplets).
    pub clustering_coefficient: f64,
    /// Average shortest path length estimate (sampled from a few sources).
    pub avg_path_length: f64,
}

impl GraphStats {
    /// Compute graph statistics from an adjacency matrix.
    ///
    /// The average path length is estimated by running BFS from up to 10
    /// sampled source nodes.
    pub fn compute(adj: &AdjacencyMatrix) -> Self {
        let n = adj.n;
        let edge_count = adj.edge_count();

        // Degree stats
        let min_degree = adj.degree.iter().copied().min().unwrap_or(0);
        let max_degree = adj.degree.iter().copied().max().unwrap_or(0);
        let avg_degree = if n == 0 {
            0.0
        } else {
            adj.degree.iter().sum::<usize>() as f64 / n as f64
        };

        let mut degree_distribution: HashMap<usize, usize> = HashMap::new();
        for &d in &adj.degree {
            *degree_distribution.entry(d).or_insert(0) += 1;
        }

        // Clustering coefficient (undirected interpretation)
        // cc(u) = (edges among neighbours) / (deg*(deg-1)/2)
        let clustering_coefficient = if n == 0 {
            0.0
        } else {
            let cc_sum: f64 = (0..n)
                .map(|u| {
                    let neighbours: Vec<usize> = adj.adj[u].iter().map(|&(v, _)| v).collect();
                    let deg = neighbours.len();
                    if deg < 2 {
                        return 0.0;
                    }
                    let mut tri = 0usize;
                    for i in 0..deg {
                        for j in (i + 1)..deg {
                            let vi = neighbours[i];
                            let vj = neighbours[j];
                            if adj.adj[vi].iter().any(|&(w, _)| w == vj) {
                                tri += 1;
                            }
                        }
                    }
                    let possible = deg * (deg - 1) / 2;
                    tri as f64 / possible as f64
                })
                .sum();
            cc_sum / n as f64
        };

        // Average path length estimate (BFS from up to 10 nodes)
        let avg_path_length = if n <= 1 {
            0.0
        } else {
            let sample_count = n.min(10);
            let step = n / sample_count;
            let mut total_dist = 0.0_f64;
            let mut total_pairs = 0usize;
            for s in (0..n).step_by(step.max(1)).take(sample_count) {
                let bfs = BfsOrder::from(adj, s);
                for (d, &dist) in bfs.distance.iter().enumerate() {
                    if d != s && dist != usize::MAX {
                        total_dist += dist as f64;
                        total_pairs += 1;
                    }
                }
            }
            if total_pairs == 0 {
                0.0
            } else {
                total_dist / total_pairs as f64
            }
        };

        GraphStats {
            node_count: n,
            edge_count,
            min_degree,
            max_degree,
            avg_degree,
            degree_distribution,
            clustering_coefficient,
            avg_path_length,
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // Helper: build a simple triangle graph (3 nodes, 3 undirected edges)
    fn triangle() -> (Vec<GraphNode>, Vec<GraphEdge>) {
        let nodes = (0..3)
            .map(|i| GraphNode::new(i, format!("n{}", i)))
            .collect();
        let edges = vec![
            GraphEdge::undirected(0, 1),
            GraphEdge::undirected(1, 2),
            GraphEdge::undirected(2, 0),
        ];
        (nodes, edges)
    }

    // Helper: build a path graph 0→1→2→3→4
    fn path_graph(n: usize) -> (Vec<GraphNode>, Vec<GraphEdge>) {
        let nodes = (0..n)
            .map(|i| GraphNode::new(i, format!("n{}", i)))
            .collect();
        let edges = (0..n.saturating_sub(1))
            .map(|i| GraphEdge::directed(i, i + 1, 1.0))
            .collect();
        (nodes, edges)
    }

    // --- GraphNode ---

    #[test]
    fn test_node_new_default_position() {
        let n = GraphNode::new(7, "alpha");
        assert_eq!(n.id, 7);
        assert_eq!(n.position, [0.0, 0.0]);
        assert_eq!(n.label, "alpha");
    }

    #[test]
    fn test_node_color_default() {
        let n = GraphNode::new(0, "x");
        // Default color should be a valid [0,1] RGB triple
        for c in n.color {
            assert!((0.0..=1.0).contains(&c));
        }
    }

    // --- GraphEdge ---

    #[test]
    fn test_edge_undirected() {
        let e = GraphEdge::undirected(1, 2);
        assert!(!e.directed);
        assert!((e.weight - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_edge_directed_weight() {
        let e = GraphEdge::directed(0, 3, 4.5);
        assert!(e.directed);
        assert!((e.weight - 4.5).abs() < 1e-12);
    }

    // --- circular_layout ---

    #[test]
    fn test_circular_layout_zero_nodes() {
        assert!(circular_layout(0, 1.0).is_empty());
    }

    #[test]
    fn test_circular_layout_single_node() {
        let pos = circular_layout(1, 5.0);
        assert_eq!(pos.len(), 1);
        assert!((pos[0][0] - 5.0).abs() < 1e-10);
        assert!(pos[0][1].abs() < 1e-10);
    }

    #[test]
    fn test_circular_layout_all_on_circle() {
        let r = 3.0;
        let positions = circular_layout(8, r);
        for p in &positions {
            let dist = (p[0] * p[0] + p[1] * p[1]).sqrt();
            assert!((dist - r).abs() < 1e-10, "point not on circle: dist={dist}");
        }
    }

    #[test]
    fn test_circular_layout_evenly_spaced() {
        let n = 6usize;
        let r = 2.0;
        let pos = circular_layout(n, r);
        // Adjacent arc lengths should be equal
        let chord0 = {
            let dx = pos[1][0] - pos[0][0];
            let dy = pos[1][1] - pos[0][1];
            (dx * dx + dy * dy).sqrt()
        };
        for i in 1..n {
            let j = (i + 1) % n;
            let dx = pos[j][0] - pos[i][0];
            let dy = pos[j][1] - pos[i][1];
            let chord = (dx * dx + dy * dy).sqrt();
            assert!((chord - chord0).abs() < 1e-10, "chord mismatch at {i}");
        }
    }

    // --- force_directed_layout ---

    #[test]
    fn test_force_directed_no_panic_empty() {
        let mut nodes: Vec<GraphNode> = vec![];
        force_directed_layout(&mut nodes, &[], 10);
    }

    #[test]
    fn test_force_directed_moves_nodes() {
        let (mut nodes, edges) = triangle();
        // All nodes start at origin; after layout they should spread out
        force_directed_layout(&mut nodes, &edges, 100);
        let all_zero = nodes.iter().all(|n| n.position == [0.0, 0.0]);
        assert!(!all_zero, "force layout should move nodes away from origin");
    }

    #[test]
    fn test_force_directed_single_node() {
        let mut nodes = vec![GraphNode::new(0, "lone")];
        force_directed_layout(&mut nodes, &[], 50);
        // Single node: no forces, stays at origin
        assert_eq!(nodes[0].position, [0.0, 0.0]);
    }

    #[test]
    fn test_force_directed_positions_finite() {
        let (mut nodes, edges) = triangle();
        force_directed_layout(&mut nodes, &edges, 200);
        for n in &nodes {
            assert!(n.position[0].is_finite());
            assert!(n.position[1].is_finite());
        }
    }

    // --- hierarchical_layout ---

    #[test]
    fn test_hierarchical_layout_path_y_monotone() {
        let (mut nodes, edges) = path_graph(5);
        hierarchical_layout(&mut nodes, &edges);
        // y-coordinates (layers) should be non-decreasing along the path
        for i in 0..4 {
            assert!(
                nodes[i].position[1] <= nodes[i + 1].position[1],
                "node {} y {} > node {} y {}",
                i,
                nodes[i].position[1],
                i + 1,
                nodes[i + 1].position[1]
            );
        }
    }

    #[test]
    fn test_hierarchical_layout_positions_finite() {
        let (mut nodes, edges) = path_graph(4);
        hierarchical_layout(&mut nodes, &edges);
        for n in &nodes {
            assert!(n.position[0].is_finite());
            assert!(n.position[1].is_finite());
        }
    }

    #[test]
    fn test_hierarchical_layout_empty_no_panic() {
        let mut nodes: Vec<GraphNode> = vec![];
        hierarchical_layout(&mut nodes, &[]);
    }

    // --- AdjacencyMatrix ---

    #[test]
    fn test_adj_matrix_from_triangle_degree() {
        let (_, edges) = triangle();
        let adj = AdjacencyMatrix::from_edges(3, &edges);
        // Each node has degree 2 in an undirected triangle
        for d in &adj.degree {
            assert_eq!(*d, 2);
        }
    }

    #[test]
    fn test_adj_matrix_is_connected_triangle() {
        let (_, edges) = triangle();
        let adj = AdjacencyMatrix::from_edges(3, &edges);
        assert!(adj.is_connected());
    }

    #[test]
    fn test_adj_matrix_disconnected() {
        // Two isolated nodes, no edges
        let adj = AdjacencyMatrix::from_edges(2, &[]);
        assert!(!adj.is_connected());
    }

    #[test]
    fn test_adj_matrix_single_node_connected() {
        let adj = AdjacencyMatrix::from_edges(1, &[]);
        assert!(adj.is_connected());
    }

    #[test]
    fn test_adj_matrix_directed_path() {
        let (_, edges) = path_graph(3);
        let adj = AdjacencyMatrix::from_edges(3, &edges);
        // Node 0 has out-degree 1, node 1 has in+out=2, node 2 has in=1
        // In our implementation degree counts both in and out for directed
        assert!(adj.edge_count() > 0);
    }

    // --- BfsOrder ---

    #[test]
    fn test_bfs_order_path() {
        let (_, edges) = path_graph(5);
        let adj = AdjacencyMatrix::from_edges(5, &edges);
        let bfs = BfsOrder::from(&adj, 0);
        assert_eq!(bfs.order[0], 0);
        // In a directed path 0→1→2→3→4, BFS visits in order
        for (i, &node) in bfs.order.iter().enumerate() {
            assert_eq!(node, i);
        }
    }

    #[test]
    fn test_bfs_distance_path() {
        let (_, edges) = path_graph(5);
        let adj = AdjacencyMatrix::from_edges(5, &edges);
        let bfs = BfsOrder::from(&adj, 0);
        for (i, &d) in bfs.distance.iter().enumerate() {
            assert_eq!(d, i, "distance to node {i} should be {i}");
        }
    }

    #[test]
    fn test_bfs_unreachable_node() {
        let adj = AdjacencyMatrix::from_edges(3, &[GraphEdge::directed(0, 1, 1.0)]);
        let bfs = BfsOrder::from(&adj, 0);
        assert_eq!(bfs.distance[2], usize::MAX);
    }

    #[test]
    fn test_bfs_invalid_src_empty() {
        let adj = AdjacencyMatrix::from_edges(2, &[]);
        let bfs = BfsOrder::from(&adj, 99);
        assert!(bfs.order.is_empty());
    }

    // --- DfsOrder ---

    #[test]
    fn test_dfs_visits_all_nodes_connected() {
        let (_, edges) = triangle();
        let adj = AdjacencyMatrix::from_edges(3, &edges);
        let dfs = DfsOrder::from(&adj, 0);
        assert_eq!(dfs.order.len(), 3);
        assert!(dfs.visited.iter().all(|&v| v));
    }

    #[test]
    fn test_dfs_disconnected_partial_visit() {
        let adj = AdjacencyMatrix::from_edges(4, &[GraphEdge::undirected(0, 1)]);
        let dfs = DfsOrder::from(&adj, 0);
        assert_eq!(dfs.order.len(), 2);
    }

    #[test]
    fn test_dfs_path_graph_order() {
        let (_, edges) = path_graph(4);
        let adj = AdjacencyMatrix::from_edges(4, &edges);
        let dfs = DfsOrder::from(&adj, 0);
        // In a directed path DFS visits 0,1,2,3 in order
        assert_eq!(dfs.order, vec![0, 1, 2, 3]);
    }

    // --- Dijkstra ---

    #[test]
    fn test_dijkstra_path_exists() {
        let edges = vec![
            GraphEdge::directed(0, 1, 1.0),
            GraphEdge::directed(1, 2, 2.0),
            GraphEdge::directed(0, 2, 10.0),
        ];
        let adj = AdjacencyMatrix::from_edges(3, &edges);
        let path = shortest_path_dijkstra(&adj, 0, 2).unwrap();
        assert_eq!(path, vec![0, 1, 2]);
    }

    #[test]
    fn test_dijkstra_no_path() {
        let adj = AdjacencyMatrix::from_edges(3, &[GraphEdge::directed(0, 1, 1.0)]);
        assert!(shortest_path_dijkstra(&adj, 0, 2).is_none());
    }

    #[test]
    fn test_dijkstra_src_equals_dst() {
        let (_, edges) = triangle();
        let adj = AdjacencyMatrix::from_edges(3, &edges);
        let path = shortest_path_dijkstra(&adj, 1, 1).unwrap();
        assert_eq!(path, vec![1]);
    }

    #[test]
    fn test_dijkstra_shortest_weight() {
        // 0→1 cost 5, 0→2→1 cost 2+2=4 → prefer 0→2→1
        let edges = vec![
            GraphEdge::directed(0, 1, 5.0),
            GraphEdge::directed(0, 2, 2.0),
            GraphEdge::directed(2, 1, 2.0),
        ];
        let adj = AdjacencyMatrix::from_edges(3, &edges);
        let path = shortest_path_dijkstra(&adj, 0, 1).unwrap();
        assert_eq!(path, vec![0, 2, 1]);
    }

    #[test]
    fn test_dijkstra_invalid_src() {
        let adj = AdjacencyMatrix::from_edges(2, &[]);
        assert!(shortest_path_dijkstra(&adj, 99, 0).is_none());
    }

    // --- GraphStats ---

    #[test]
    fn test_stats_triangle_degrees() {
        let (_, edges) = triangle();
        let adj = AdjacencyMatrix::from_edges(3, &edges);
        let stats = GraphStats::compute(&adj);
        assert_eq!(stats.min_degree, 2);
        assert_eq!(stats.max_degree, 2);
        assert!((stats.avg_degree - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_stats_empty_graph() {
        let adj = AdjacencyMatrix::from_edges(0, &[]);
        let stats = GraphStats::compute(&adj);
        assert_eq!(stats.node_count, 0);
        assert_eq!(stats.edge_count, 0);
    }

    #[test]
    fn test_stats_clustering_complete_triangle() {
        // A complete triangle (every pair connected undirected) → cc = 1.0
        let (_, edges) = triangle();
        let adj = AdjacencyMatrix::from_edges(3, &edges);
        let stats = GraphStats::compute(&adj);
        assert!(
            stats.clustering_coefficient > 0.9,
            "triangle clustering should be near 1, got {}",
            stats.clustering_coefficient
        );
    }

    #[test]
    fn test_stats_degree_distribution_count() {
        let (_, edges) = triangle();
        let adj = AdjacencyMatrix::from_edges(3, &edges);
        let stats = GraphStats::compute(&adj);
        // All 3 nodes have degree 2
        assert_eq!(*stats.degree_distribution.get(&2).unwrap(), 3);
    }

    #[test]
    fn test_stats_avg_path_length_path() {
        let (_, edges) = path_graph(5);
        let adj = AdjacencyMatrix::from_edges(5, &edges);
        let stats = GraphStats::compute(&adj);
        // Average path length for directed path 0→1→2→3→4:
        // Node 0 can reach 1,2,3,4 with distances 1,2,3,4 (avg 2.5)
        // Other source nodes can only reach forward-nodes
        assert!(stats.avg_path_length > 0.0);
    }

    // --- GraphLayout ---

    #[test]
    fn test_graph_layout_bounding_box() {
        let nodes = vec![
            {
                let mut n = GraphNode::new(0, "a");
                n.position = [-1.0, -2.0];
                n
            },
            {
                let mut n = GraphNode::new(1, "b");
                n.position = [3.0, 4.0];
                n
            },
        ];
        let layout = GraphLayout::new(nodes, vec![]);
        let (mn, mx) = layout.bounding_box().unwrap();
        assert!((mn[0] + 1.0).abs() < 1e-10);
        assert!((mn[1] + 2.0).abs() < 1e-10);
        assert!((mx[0] - 3.0).abs() < 1e-10);
        assert!((mx[1] - 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_graph_layout_bounding_box_empty() {
        let layout = GraphLayout::new(vec![], vec![]);
        assert!(layout.bounding_box().is_none());
    }

    // --- dist2 helper ---

    #[test]
    fn test_dist2_zero() {
        assert!(dist2([1.0, 2.0], [1.0, 2.0]) <= 1e-9);
    }

    #[test]
    fn test_dist2_known() {
        let d = dist2([0.0, 0.0], [3.0, 4.0]);
        assert!((d - 5.0).abs() < 1e-10);
    }
}

// ---------------------------------------------------------------------------
// NodeAttribute — node coloring by scalar or categorical attribute
// ---------------------------------------------------------------------------

/// A per-node scalar or categorical attribute used for coloring.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub enum NodeAttribute {
    /// A continuous scalar value (e.g. centrality, density).
    Scalar(f64),
    /// An integer category id (e.g. community label).
    Category(usize),
}

/// Map a scalar attribute to an RGB color using a simple gradient.
///
/// `t` in `[0,1]` maps blue→green→red.
#[allow(dead_code)]
pub fn scalar_to_color(t: f64) -> [f32; 3] {
    let t = t.clamp(0.0, 1.0) as f32;
    if t < 0.5 {
        let s = t * 2.0;
        [0.0, s, 1.0 - s]
    } else {
        let s = (t - 0.5) * 2.0;
        [s, 1.0 - s, 0.0]
    }
}

/// Assign colors to nodes based on a per-node scalar attribute.
///
/// The values are normalised to `[0,1]` before mapping.
/// Returns a `Vec` of RGB triples in the same order as `nodes`.
#[allow(dead_code)]
pub fn color_nodes_by_scalar(nodes: &mut Vec<GraphNode>, values: &[f64]) {
    if nodes.is_empty() || values.is_empty() {
        return;
    }
    let n = nodes.len().min(values.len());
    let min_v = values[..n].iter().cloned().fold(f64::INFINITY, f64::min);
    let max_v = values[..n]
        .iter()
        .cloned()
        .fold(f64::NEG_INFINITY, f64::max);
    let range = if (max_v - min_v).abs() < 1e-15 {
        1.0
    } else {
        max_v - min_v
    };
    for i in 0..n {
        let t = (values[i] - min_v) / range;
        nodes[i].color = scalar_to_color(t);
    }
}

/// Assign distinct colors to nodes based on a categorical attribute (e.g. community id).
///
/// Up to 12 distinct hues are used; extra categories wrap around.
#[allow(dead_code)]
pub fn color_nodes_by_category(nodes: &mut Vec<GraphNode>, categories: &[usize]) {
    // 12 visually distinct hues (HSV h=0,30,60,...,330)
    let palette: [[f32; 3]; 12] = [
        [0.894, 0.102, 0.110], // red
        [0.216, 0.494, 0.722], // blue
        [0.302, 0.686, 0.290], // green
        [0.596, 0.306, 0.639], // purple
        [1.000, 0.498, 0.000], // orange
        [1.000, 1.000, 0.200], // yellow
        [0.651, 0.337, 0.157], // brown
        [0.969, 0.506, 0.749], // pink
        [0.600, 0.600, 0.600], // grey
        [0.400, 0.761, 0.647], // teal
        [0.988, 0.749, 0.435], // peach
        [0.388, 0.533, 0.176], // olive
    ];
    let n = nodes.len().min(categories.len());
    for i in 0..n {
        nodes[i].color = palette[categories[i] % 12];
    }
}

// ---------------------------------------------------------------------------
// EdgeWeightVisualizer — width / alpha scaling by edge weight
// ---------------------------------------------------------------------------

/// Rendering hint for a single edge: display width and alpha.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct EdgeRenderHint {
    /// Display width in pixels (or abstract units).
    pub width: f32,
    /// Alpha transparency in `[0, 1]`.
    pub alpha: f32,
    /// An interpolated colour `[r, g, b]`.
    pub color: [f32; 3],
}

/// Compute rendering hints for a set of edges based on their weights.
///
/// Weights are normalised to `[0, 1]`; the width is scaled between
/// `min_width` and `max_width`, and the alpha between `min_alpha` and 1.
#[allow(dead_code)]
#[allow(clippy::too_many_arguments)]
pub fn edge_weight_hints(
    edges: &[GraphEdge],
    min_width: f32,
    max_width: f32,
    min_alpha: f32,
    low_color: [f32; 3],
    high_color: [f32; 3],
) -> Vec<EdgeRenderHint> {
    if edges.is_empty() {
        return vec![];
    }
    let min_w = edges.iter().map(|e| e.weight).fold(f64::INFINITY, f64::min);
    let max_w = edges
        .iter()
        .map(|e| e.weight)
        .fold(f64::NEG_INFINITY, f64::max);
    let range = if (max_w - min_w).abs() < 1e-15 {
        1.0
    } else {
        max_w - min_w
    };
    edges
        .iter()
        .map(|e| {
            let t = ((e.weight - min_w) / range) as f32;
            EdgeRenderHint {
                width: min_width + t * (max_width - min_width),
                alpha: min_alpha + t * (1.0 - min_alpha),
                color: [
                    low_color[0] + t * (high_color[0] - low_color[0]),
                    low_color[1] + t * (high_color[1] - low_color[1]),
                    low_color[2] + t * (high_color[2] - low_color[2]),
                ],
            }
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Community detection — Louvain-lite (label propagation)
// ---------------------------------------------------------------------------

/// Run label-propagation community detection.
///
/// Each node starts with its own community label.  In each iteration every
/// node adopts the most frequent label among its neighbours (ties broken by
/// lowest label).  Returns a community label for each node.
#[allow(dead_code)]
pub fn label_propagation_communities(adj: &AdjacencyMatrix, iters: usize) -> Vec<usize> {
    let n = adj.n;
    if n == 0 {
        return vec![];
    }
    let mut labels: Vec<usize> = (0..n).collect();
    for _iter in 0..iters {
        let prev = labels.clone();
        for u in 0..n {
            let mut freq: HashMap<usize, usize> = HashMap::new();
            for &(v, _) in &adj.adj[u] {
                *freq.entry(prev[v]).or_insert(0) += 1;
            }
            if freq.is_empty() {
                continue;
            }
            let max_count = *freq.values().max().expect("iterator should not be empty");
            // Tie-break: pick the smallest label with max count
            let best = freq
                .iter()
                .filter(|&(_, &c)| c == max_count)
                .map(|(&l, _)| l)
                .min()
                .unwrap_or(labels[u]);
            labels[u] = best;
        }
    }
    labels
}

/// Colour `nodes` based on community labels returned by [`label_propagation_communities`].
#[allow(dead_code)]
pub fn color_nodes_by_community(nodes: &mut Vec<GraphNode>, communities: &[usize]) {
    color_nodes_by_category(nodes, communities);
}

/// Compute the number of distinct communities in a label assignment.
#[allow(dead_code)]
pub fn community_count(labels: &[usize]) -> usize {
    let mut seen: std::collections::HashSet<usize> = std::collections::HashSet::new();
    for &l in labels {
        seen.insert(l);
    }
    seen.len()
}

// ---------------------------------------------------------------------------
// Graph edge bundling (straight-line approximation)
// ---------------------------------------------------------------------------

/// A bundled-edge control point polyline.
///
/// In a real implementation this would follow a spine; here we use a simple
/// quadratic Bézier midpoint approximation to group nearby edges.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct BundledEdge {
    /// Source position.
    pub src: [f64; 2],
    /// Control point (midpoint pulled toward graph centroid).
    pub ctrl: [f64; 2],
    /// Destination position.
    pub dst: [f64; 2],
}

impl BundledEdge {
    /// Sample the Bézier curve at parameter `t ∈ [0,1]`.
    pub fn sample(&self, t: f64) -> [f64; 2] {
        let t = t.clamp(0.0, 1.0);
        let u = 1.0 - t;
        [
            u * u * self.src[0] + 2.0 * u * t * self.ctrl[0] + t * t * self.dst[0],
            u * u * self.src[1] + 2.0 * u * t * self.ctrl[1] + t * t * self.dst[1],
        ]
    }

    /// Sample the curve into `n_samples` equally spaced points.
    pub fn polyline(&self, n_samples: usize) -> Vec<[f64; 2]> {
        if n_samples == 0 {
            return vec![];
        }
        (0..n_samples)
            .map(|i| self.sample(i as f64 / (n_samples - 1).max(1) as f64))
            .collect()
    }
}

/// Bundle all edges toward the graph centroid.
///
/// The control point is the midpoint between the edge midpoint and the
/// centroid, scaled by `bundle_strength ∈ [0,1]`.
#[allow(dead_code)]
pub fn bundle_edges(
    nodes: &[GraphNode],
    edges: &[GraphEdge],
    bundle_strength: f64,
) -> Vec<BundledEdge> {
    let n = nodes.len();
    if n == 0 {
        return vec![];
    }
    // Centroid
    let cx = nodes.iter().map(|nd| nd.position[0]).sum::<f64>() / n as f64;
    let cy = nodes.iter().map(|nd| nd.position[1]).sum::<f64>() / n as f64;
    let s = bundle_strength.clamp(0.0, 1.0);

    edges
        .iter()
        .map(|e| {
            let si = e.src.min(n - 1);
            let di = e.dst.min(n - 1);
            let sp = nodes[si].position;
            let dp = nodes[di].position;
            // Natural midpoint
            let mx = (sp[0] + dp[0]) * 0.5;
            let my = (sp[1] + dp[1]) * 0.5;
            // Pull toward centroid
            let ctrl = [mx + s * (cx - mx), my + s * (cy - my)];
            BundledEdge {
                src: sp,
                ctrl,
                dst: dp,
            }
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Adjacency-matrix heatmap
// ---------------------------------------------------------------------------

/// A dense adjacency matrix heatmap (all entries as f64 weight or 0/1).
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct AdjacencyHeatmap {
    /// Number of nodes.
    pub n: usize,
    /// Dense weight matrix, row-major (`matrix[i * n + j]`).
    pub matrix: Vec<f64>,
}

impl AdjacencyHeatmap {
    /// Build a dense heatmap from an adjacency list.
    pub fn from_adj(adj: &AdjacencyMatrix) -> Self {
        let n = adj.n;
        let mut matrix = vec![0.0f64; n * n];
        for u in 0..n {
            for &(v, w) in &adj.adj[u] {
                matrix[u * n + v] = w;
            }
        }
        AdjacencyHeatmap { n, matrix }
    }

    /// Retrieve the weight between nodes `u` and `v`.
    pub fn get(&self, u: usize, v: usize) -> f64 {
        if u >= self.n || v >= self.n {
            return 0.0;
        }
        self.matrix[u * self.n + v]
    }

    /// Map the matrix values to RGBA colors using the scalar color gradient.
    ///
    /// Returns a flat `Vec<[u8; 4]>` of length `n*n` (RGBA per cell).
    pub fn to_rgba(&self) -> Vec<[u8; 4]> {
        let max_v = self.matrix.iter().cloned().fold(0.0f64, f64::max);
        let scale = if max_v < 1e-15 { 1.0 } else { max_v };
        self.matrix
            .iter()
            .map(|&v| {
                let t = v / scale;
                let [r, g, b] = scalar_to_color(t);
                [(r * 255.0) as u8, (g * 255.0) as u8, (b * 255.0) as u8, 255]
            })
            .collect()
    }

    /// Export as CSV (rows = source, cols = destination).
    pub fn to_csv(&self) -> String {
        let mut s = String::new();
        for i in 0..self.n {
            let row: Vec<String> = (0..self.n)
                .map(|j| format!("{:.4}", self.matrix[i * self.n + j]))
                .collect();
            s.push_str(&row.join(","));
            s.push('\n');
        }
        s
    }
}

// ---------------------------------------------------------------------------
// Interactive graph navigation data
// ---------------------------------------------------------------------------

/// Viewport state for interactive graph panning and zooming.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct GraphViewport {
    /// Centre of the viewport in graph-space coordinates.
    pub center: [f64; 2],
    /// Zoom level (1.0 = 100%, 2.0 = 200%, etc.).
    pub zoom: f64,
    /// Viewport width in pixels.
    pub width_px: u32,
    /// Viewport height in pixels.
    pub height_px: u32,
}

impl GraphViewport {
    /// Create a default viewport centred at the origin.
    pub fn new(width_px: u32, height_px: u32) -> Self {
        GraphViewport {
            center: [0.0, 0.0],
            zoom: 1.0,
            width_px,
            height_px,
        }
    }

    /// Convert a graph-space coordinate to pixel (screen) coordinates.
    pub fn graph_to_screen(&self, pt: [f64; 2]) -> [f64; 2] {
        let cx = self.width_px as f64 * 0.5;
        let cy = self.height_px as f64 * 0.5;
        [
            cx + (pt[0] - self.center[0]) * self.zoom,
            cy - (pt[1] - self.center[1]) * self.zoom, // y flipped for screen
        ]
    }

    /// Convert screen coordinates back to graph space.
    pub fn screen_to_graph(&self, pt: [f64; 2]) -> [f64; 2] {
        let cx = self.width_px as f64 * 0.5;
        let cy = self.height_px as f64 * 0.5;
        [
            self.center[0] + (pt[0] - cx) / self.zoom,
            self.center[1] - (pt[1] - cy) / self.zoom,
        ]
    }

    /// Pan by a delta in screen pixels.
    pub fn pan(&mut self, dx_px: f64, dy_px: f64) {
        self.center[0] -= dx_px / self.zoom;
        self.center[1] += dy_px / self.zoom;
    }

    /// Zoom by a factor around a screen-space anchor point.
    pub fn zoom_around(&mut self, factor: f64, anchor_px: [f64; 2]) {
        let before = self.screen_to_graph(anchor_px);
        self.zoom = (self.zoom * factor).clamp(0.01, 1000.0);
        let after = self.screen_to_graph(anchor_px);
        self.center[0] += before[0] - after[0];
        self.center[1] += before[1] - after[1];
    }

    /// Return the visible graph-space bounding box `([xmin,ymin],[xmax,ymax])`.
    pub fn visible_bounds(&self) -> ([f64; 2], [f64; 2]) {
        let hw = self.width_px as f64 * 0.5 / self.zoom;
        let hh = self.height_px as f64 * 0.5 / self.zoom;
        (
            [self.center[0] - hw, self.center[1] - hh],
            [self.center[0] + hw, self.center[1] + hh],
        )
    }

    /// Return ids of nodes that fall within the current visible bounds.
    pub fn visible_nodes<'a>(&self, nodes: &'a [GraphNode]) -> Vec<&'a GraphNode> {
        let (lo, hi) = self.visible_bounds();
        nodes
            .iter()
            .filter(|nd| {
                nd.position[0] >= lo[0]
                    && nd.position[0] <= hi[0]
                    && nd.position[1] >= lo[1]
                    && nd.position[1] <= hi[1]
            })
            .collect()
    }
}

/// Selection state for interactive graph editing.
#[allow(dead_code)]
#[derive(Debug, Clone, Default)]
pub struct GraphSelection {
    /// Currently selected node ids.
    pub selected_nodes: std::collections::HashSet<usize>,
    /// Currently selected edge indices.
    pub selected_edges: std::collections::HashSet<usize>,
}

impl GraphSelection {
    /// Create an empty selection.
    pub fn new() -> Self {
        GraphSelection {
            selected_nodes: std::collections::HashSet::new(),
            selected_edges: std::collections::HashSet::new(),
        }
    }

    /// Toggle a node's selection state. Returns `true` if the node is now selected.
    pub fn toggle_node(&mut self, id: usize) -> bool {
        if self.selected_nodes.contains(&id) {
            self.selected_nodes.remove(&id);
            false
        } else {
            self.selected_nodes.insert(id);
            true
        }
    }

    /// Select all nodes in a rectangular region (graph space).
    pub fn select_rect(&mut self, nodes: &[GraphNode], lo: [f64; 2], hi: [f64; 2]) {
        for nd in nodes {
            if nd.position[0] >= lo[0]
                && nd.position[0] <= hi[0]
                && nd.position[1] >= lo[1]
                && nd.position[1] <= hi[1]
            {
                self.selected_nodes.insert(nd.id);
            }
        }
    }

    /// Clear the entire selection.
    pub fn clear(&mut self) {
        self.selected_nodes.clear();
        self.selected_edges.clear();
    }

    /// Number of selected nodes.
    pub fn node_count(&self) -> usize {
        self.selected_nodes.len()
    }
}

// ---------------------------------------------------------------------------
// Additional tests for the new sections
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests_extended {
    use super::*;

    // --- scalar_to_color ---

    #[test]
    fn test_scalar_to_color_zero_is_blue() {
        let c = scalar_to_color(0.0);
        assert!(c[2] > 0.9, "t=0 should be bluish");
    }

    #[test]
    fn test_scalar_to_color_one_is_red() {
        let c = scalar_to_color(1.0);
        assert!(c[0] > 0.9, "t=1 should be reddish");
    }

    #[test]
    fn test_scalar_to_color_clamp() {
        let lo = scalar_to_color(-5.0);
        let hi = scalar_to_color(5.0);
        assert_eq!(lo, scalar_to_color(0.0));
        assert_eq!(hi, scalar_to_color(1.0));
    }

    // --- color_nodes_by_scalar ---

    #[test]
    fn test_color_nodes_by_scalar_uniform() {
        let mut nodes: Vec<GraphNode> = (0..4).map(|i| GraphNode::new(i, "")).collect();
        let values = vec![2.0, 2.0, 2.0, 2.0];
        color_nodes_by_scalar(&mut nodes, &values);
        // All same value → all same colour (mapped at t=0 due to zero range)
        let c0 = nodes[0].color;
        for nd in &nodes {
            assert_eq!(nd.color, c0);
        }
    }

    #[test]
    fn test_color_nodes_by_scalar_extremes() {
        let mut nodes: Vec<GraphNode> = (0..2).map(|i| GraphNode::new(i, "")).collect();
        let values = vec![0.0, 1.0];
        color_nodes_by_scalar(&mut nodes, &values);
        // Node 0 → blue, node 1 → red
        assert!(nodes[0].color[2] > nodes[1].color[2]); // node0 more blue
        assert!(nodes[1].color[0] > nodes[0].color[0]); // node1 more red
    }

    // --- color_nodes_by_category ---

    #[test]
    fn test_color_nodes_by_category_different_categories() {
        let mut nodes: Vec<GraphNode> = (0..3).map(|i| GraphNode::new(i, "")).collect();
        color_nodes_by_category(&mut nodes, &[0, 1, 2]);
        // Each node should get a different colour
        let c0 = nodes[0].color;
        let c1 = nodes[1].color;
        assert_ne!(c0, c1);
    }

    #[test]
    fn test_color_nodes_by_category_wraps() {
        let mut nodes: Vec<GraphNode> = (0..2).map(|i| GraphNode::new(i, "")).collect();
        // category 0 and 12 should both map to palette[0]
        color_nodes_by_category(&mut nodes, &[0, 12]);
        assert_eq!(nodes[0].color, nodes[1].color);
    }

    // --- edge_weight_hints ---

    #[test]
    fn test_edge_weight_hints_count() {
        let edges = vec![
            GraphEdge::directed(0, 1, 1.0),
            GraphEdge::directed(1, 2, 3.0),
        ];
        let hints = edge_weight_hints(&edges, 1.0, 5.0, 0.2, [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]);
        assert_eq!(hints.len(), 2);
    }

    #[test]
    fn test_edge_weight_hints_width_monotone() {
        let edges: Vec<GraphEdge> = (0..5)
            .map(|i| GraphEdge::directed(0, i + 1, i as f64))
            .collect();
        let hints = edge_weight_hints(&edges, 1.0, 10.0, 0.1, [0.0; 3], [1.0; 3]);
        for w in hints.windows(2) {
            assert!(w[0].width <= w[1].width + 1e-6);
        }
    }

    #[test]
    fn test_edge_weight_hints_empty() {
        let hints = edge_weight_hints(&[], 1.0, 5.0, 0.5, [0.0; 3], [1.0; 3]);
        assert!(hints.is_empty());
    }

    // --- label_propagation_communities ---

    #[test]
    fn test_community_detection_connected_graph() {
        // Triangle graph: after propagation all nodes should end up in the same community
        let edges = vec![
            GraphEdge::undirected(0, 1),
            GraphEdge::undirected(1, 2),
            GraphEdge::undirected(2, 0),
        ];
        let adj = AdjacencyMatrix::from_edges(3, &edges);
        let labels = label_propagation_communities(&adj, 20);
        // All 3 nodes should share one community
        assert_eq!(community_count(&labels), 1);
    }

    #[test]
    fn test_community_detection_two_components() {
        // Two isolated pairs: 0-1, 2-3
        let edges = vec![GraphEdge::undirected(0, 1), GraphEdge::undirected(2, 3)];
        let adj = AdjacencyMatrix::from_edges(4, &edges);
        let labels = label_propagation_communities(&adj, 20);
        // Nodes in different components should have different community labels
        assert_ne!(labels[0], labels[2]);
    }

    #[test]
    fn test_community_count_single() {
        assert_eq!(community_count(&[0, 0, 0]), 1);
    }

    #[test]
    fn test_community_count_all_distinct() {
        assert_eq!(community_count(&[0, 1, 2, 3]), 4);
    }

    // --- bundle_edges ---

    #[test]
    fn test_bundle_edges_count() {
        let mut nodes: Vec<GraphNode> = (0..3).map(|i| GraphNode::new(i, "")).collect();
        nodes[0].position = [0.0, 0.0];
        nodes[1].position = [1.0, 0.0];
        nodes[2].position = [0.5, 1.0];
        let edges = vec![GraphEdge::undirected(0, 1), GraphEdge::undirected(1, 2)];
        let bundled = bundle_edges(&nodes, &edges, 0.5);
        assert_eq!(bundled.len(), 2);
    }

    #[test]
    fn test_bundled_edge_sample_endpoints() {
        let be = BundledEdge {
            src: [0.0, 0.0],
            ctrl: [0.5, 1.0],
            dst: [1.0, 0.0],
        };
        let s0 = be.sample(0.0);
        let s1 = be.sample(1.0);
        assert!((s0[0] - 0.0).abs() < 1e-12);
        assert!((s1[0] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_bundled_edge_polyline_length() {
        let be = BundledEdge {
            src: [0.0, 0.0],
            ctrl: [0.5, 1.0],
            dst: [1.0, 0.0],
        };
        assert_eq!(be.polyline(5).len(), 5);
        assert!(be.polyline(0).is_empty());
    }

    // --- AdjacencyHeatmap ---

    #[test]
    fn test_heatmap_from_adj_triangle() {
        let edges = vec![
            GraphEdge::undirected(0, 1),
            GraphEdge::undirected(1, 2),
            GraphEdge::undirected(2, 0),
        ];
        let adj = AdjacencyMatrix::from_edges(3, &edges);
        let hm = AdjacencyHeatmap::from_adj(&adj);
        assert_eq!(hm.n, 3);
        // Undirected edges: both (0,1) and (1,0) should be 1.0
        assert!((hm.get(0, 1) - 1.0).abs() < 1e-12);
        assert!((hm.get(1, 0) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_heatmap_to_rgba_size() {
        let edges = vec![GraphEdge::undirected(0, 1)];
        let adj = AdjacencyMatrix::from_edges(2, &edges);
        let hm = AdjacencyHeatmap::from_adj(&adj);
        let rgba = hm.to_rgba();
        assert_eq!(rgba.len(), 4); // 2*2
    }

    #[test]
    fn test_heatmap_to_csv_lines() {
        let edges = vec![GraphEdge::undirected(0, 1)];
        let adj = AdjacencyMatrix::from_edges(2, &edges);
        let hm = AdjacencyHeatmap::from_adj(&adj);
        let csv = hm.to_csv();
        let lines: Vec<&str> = csv.lines().collect();
        assert_eq!(lines.len(), 2); // one per row
    }

    // --- GraphViewport ---

    #[test]
    fn test_viewport_graph_to_screen_center() {
        let vp = GraphViewport::new(800, 600);
        let s = vp.graph_to_screen([0.0, 0.0]);
        assert!((s[0] - 400.0).abs() < 1e-10);
        assert!((s[1] - 300.0).abs() < 1e-10);
    }

    #[test]
    fn test_viewport_round_trip() {
        let vp = GraphViewport::new(800, 600);
        let gpt = [1.5, -2.3];
        let spt = vp.graph_to_screen(gpt);
        let gpt2 = vp.screen_to_graph(spt);
        assert!((gpt2[0] - gpt[0]).abs() < 1e-10);
        assert!((gpt2[1] - gpt[1]).abs() < 1e-10);
    }

    #[test]
    fn test_viewport_pan() {
        let mut vp = GraphViewport::new(800, 600);
        vp.pan(100.0, 0.0);
        // Panning right (positive dx) moves centre left in graph space
        assert!(vp.center[0] < 0.0);
    }

    #[test]
    fn test_viewport_zoom_around_anchor() {
        let mut vp = GraphViewport::new(800, 600);
        let initial_zoom = vp.zoom;
        vp.zoom_around(2.0, [400.0, 300.0]);
        assert!((vp.zoom - initial_zoom * 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_viewport_visible_bounds_grow_on_zoom_out() {
        let mut vp = GraphViewport::new(800, 600);
        let (lo1, hi1) = vp.visible_bounds();
        vp.zoom = 0.5;
        let (lo2, hi2) = vp.visible_bounds();
        assert!(hi2[0] - lo2[0] > hi1[0] - lo1[0]);
    }

    #[test]
    fn test_viewport_visible_nodes() {
        let vp = GraphViewport::new(800, 600);
        let mut nodes: Vec<GraphNode> = (0..3).map(|i| GraphNode::new(i, "")).collect();
        nodes[0].position = [0.0, 0.0]; // visible (within ±400, ±300 with zoom=1)
        nodes[1].position = [1000.0, 0.0]; // not visible
        nodes[2].position = [-200.0, 100.0]; // visible
        let vis = vp.visible_nodes(&nodes);
        let ids: Vec<usize> = vis.iter().map(|n| n.id).collect();
        assert!(ids.contains(&0));
        assert!(!ids.contains(&1));
        assert!(ids.contains(&2));
    }

    // --- GraphSelection ---

    #[test]
    fn test_selection_toggle() {
        let mut sel = GraphSelection::new();
        assert!(sel.toggle_node(5)); // now selected
        assert!(!sel.toggle_node(5)); // now deselected
        assert_eq!(sel.node_count(), 0);
    }

    #[test]
    fn test_selection_rect() {
        let mut sel = GraphSelection::new();
        let mut nodes: Vec<GraphNode> = (0..4).map(|i| GraphNode::new(i, "")).collect();
        nodes[0].position = [0.5, 0.5];
        nodes[1].position = [2.0, 2.0];
        nodes[2].position = [0.8, 0.8];
        nodes[3].position = [-1.0, -1.0];
        sel.select_rect(&nodes, [0.0, 0.0], [1.0, 1.0]);
        assert!(sel.selected_nodes.contains(&0));
        assert!(!sel.selected_nodes.contains(&1));
        assert!(sel.selected_nodes.contains(&2));
    }

    #[test]
    fn test_selection_clear() {
        let mut sel = GraphSelection::new();
        sel.toggle_node(1);
        sel.toggle_node(2);
        sel.clear();
        assert_eq!(sel.node_count(), 0);
    }
}
