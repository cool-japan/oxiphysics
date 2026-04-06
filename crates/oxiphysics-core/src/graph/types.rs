//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use std::cmp::Ordering;

#[derive(PartialEq)]
pub(super) struct AStarState {
    pub(super) f: f64,
    pub(super) g: f64,
    pub(super) node: usize,
}
/// A directed weighted graph stored as an edge list.
///
/// Undirected graphs are represented by storing both directed edges.
#[allow(dead_code)]
pub struct Graph {
    /// Number of nodes.
    pub n_nodes: usize,
    /// Edge list: (src, dst, weight).
    pub edges: Vec<(usize, usize, f64)>,
}
impl Graph {
    /// Creates a new empty graph with the given number of nodes.
    #[allow(dead_code)]
    pub fn new(n_nodes: usize) -> Self {
        Self {
            n_nodes,
            edges: Vec::new(),
        }
    }
    /// Adds a directed edge from `src` to `dst` with `weight`.
    #[allow(dead_code)]
    pub fn add_edge(&mut self, src: usize, dst: usize, weight: f64) {
        assert!(src < self.n_nodes, "src node out of range");
        assert!(dst < self.n_nodes, "dst node out of range");
        self.edges.push((src, dst, weight));
    }
    /// Adds an undirected edge between `a` and `b` with `weight` (stores both directions).
    #[allow(dead_code)]
    pub fn add_undirected_edge(&mut self, a: usize, b: usize, weight: f64) {
        self.add_edge(a, b, weight);
        self.add_edge(b, a, weight);
    }
    /// Returns the adjacency list: for each node, a list of `(neighbor, weight)`.
    #[allow(dead_code)]
    pub fn adjacency_list(&self) -> Vec<Vec<(usize, f64)>> {
        let mut adj = vec![Vec::new(); self.n_nodes];
        for &(src, dst, w) in &self.edges {
            adj[src].push((dst, w));
        }
        adj
    }
    /// Returns the number of nodes.
    #[allow(dead_code)]
    pub fn node_count(&self) -> usize {
        self.n_nodes
    }
    /// Returns the number of directed edges stored.
    #[allow(dead_code)]
    pub fn edge_count(&self) -> usize {
        self.edges.len()
    }
    /// Returns the out-degree of `node` (number of directed edges leaving it).
    #[allow(dead_code)]
    pub fn degree(&self, node: usize) -> usize {
        self.edges.iter().filter(|&&(s, _, _)| s == node).count()
    }
    /// Compute strongly connected components using Tarjan's algorithm.
    ///
    /// Delegates to the free function `tarjan_scc`.
    #[allow(dead_code)]
    pub fn strongly_connected_components(&self) -> Vec<Vec<usize>> {
        tarjan_scc(self)
    }
    /// Compute a minimum spanning tree using Kruskal's algorithm.
    ///
    /// Returns a list of edges `(u, v, weight)` forming the MST.
    /// Delegates to the free function `kruskal_mst`.
    #[allow(dead_code)]
    pub fn minimum_spanning_tree(&self) -> Vec<(usize, usize, f64)> {
        kruskal_mst(self)
    }
    /// Compute maximum flow from `source` to `sink` using Edmonds-Karp.
    ///
    /// Delegates to the free function `max_flow_edmonds_karp`.
    #[allow(dead_code)]
    pub fn max_flow(&self, source: usize, sink: usize) -> f64 {
        max_flow_edmonds_karp(self, source, sink)
    }
    /// Compute betweenness centrality for all nodes.
    ///
    /// Delegates to the free function `betweenness_centrality`.
    #[allow(dead_code)]
    pub fn betweenness_centrality(&self) -> Vec<f64> {
        betweenness_centrality(self)
    }
}
pub(super) struct UnionFind {
    pub(super) parent: Vec<usize>,
    pub(super) rank: Vec<usize>,
}
impl UnionFind {
    pub(super) fn new(n: usize) -> Self {
        Self {
            parent: (0..n).collect(),
            rank: vec![0; n],
        }
    }
    pub(super) fn find(&mut self, x: usize) -> usize {
        if self.parent[x] != x {
            self.parent[x] = self.find(self.parent[x]);
        }
        self.parent[x]
    }
    pub(super) fn union(&mut self, a: usize, b: usize) -> bool {
        let ra = self.find(a);
        let rb = self.find(b);
        if ra == rb {
            return false;
        }
        match self.rank[ra].cmp(&self.rank[rb]) {
            Ordering::Less => self.parent[ra] = rb,
            Ordering::Greater => self.parent[rb] = ra,
            Ordering::Equal => {
                self.parent[rb] = ra;
                self.rank[ra] += 1;
            }
        }
        true
    }
}
/// Wrapper for Dijkstra's priority queue.
#[derive(PartialEq)]
pub(super) struct State {
    pub(super) cost: f64,
    pub(super) node: usize,
}
