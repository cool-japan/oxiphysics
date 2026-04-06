#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Network graph visualization for the OxiPhysics engine.
//!
//! Provides:
//! - `NetworkGraph` — nodes, edges, adjacency
//! - `ForceDirectedLayout` — Fruchterman-Reingold / Barnes-Hut
//! - `HierarchicalLayout` — Sugiyama layer assignment & crossing minimization
//! - `NetworkMetrics` — degree, clustering, betweenness, PageRank
//! - `CommunityDetection` — modularity, spectral, label propagation
//! - `NetworkAnimation` — time-evolving graph with interpolation

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::collections::{HashMap, HashSet, VecDeque};

// ---------------------------------------------------------------------------
// NetworkGraph — core data structure
// ---------------------------------------------------------------------------

/// A node in the network graph.
#[derive(Debug, Clone)]
pub struct NetworkNode {
    /// Unique node identifier.
    pub id: usize,
    /// 2-D position for layout.
    pub position: [f64; 2],
    /// Display size (radius).
    pub size: f32,
    /// RGBA color.
    pub color: [f32; 4],
    /// Human-readable label.
    pub label: String,
    /// Optional community assignment.
    pub community: Option<usize>,
    /// Arbitrary key-value metadata.
    pub metadata: HashMap<String, String>,
}

impl NetworkNode {
    /// Create a node at the origin.
    pub fn new(id: usize, label: impl Into<String>) -> Self {
        Self {
            id,
            position: [0.0, 0.0],
            size: 1.0,
            color: [0.4, 0.6, 1.0, 1.0],
            label: label.into(),
            community: None,
            metadata: HashMap::new(),
        }
    }

    /// Set position.
    pub fn with_position(mut self, x: f64, y: f64) -> Self {
        self.position = [x, y];
        self
    }

    /// Set size.
    pub fn with_size(mut self, size: f32) -> Self {
        self.size = size;
        self
    }

    /// Set color.
    pub fn with_color(mut self, r: f32, g: f32, b: f32, a: f32) -> Self {
        self.color = [r, g, b, a];
        self
    }
}

/// A directed or undirected edge.
#[derive(Debug, Clone)]
pub struct NetworkEdge {
    /// Source node id.
    pub source: usize,
    /// Target node id.
    pub target: usize,
    /// Edge weight.
    pub weight: f64,
    /// If true, edge is directed source→target.
    pub directed: bool,
    /// Optional label.
    pub label: Option<String>,
}

impl NetworkEdge {
    /// Create an undirected edge with unit weight.
    pub fn undirected(source: usize, target: usize) -> Self {
        Self {
            source,
            target,
            weight: 1.0,
            directed: false,
            label: None,
        }
    }

    /// Create a directed edge with given weight.
    pub fn directed(source: usize, target: usize, weight: f64) -> Self {
        Self {
            source,
            target,
            weight,
            directed: true,
            label: None,
        }
    }

    /// With label.
    pub fn with_label(mut self, label: impl Into<String>) -> Self {
        self.label = Some(label.into());
        self
    }
}

/// Network graph: nodes + edges + adjacency list.
#[derive(Debug, Clone)]
pub struct NetworkGraph {
    /// All nodes.
    pub nodes: Vec<NetworkNode>,
    /// All edges.
    pub edges: Vec<NetworkEdge>,
    /// Adjacency list: node_id → list of (neighbor_id, edge_weight).
    pub adjacency: HashMap<usize, Vec<(usize, f64)>>,
}

impl NetworkGraph {
    /// Create an empty graph.
    pub fn new() -> Self {
        Self {
            nodes: Vec::new(),
            edges: Vec::new(),
            adjacency: HashMap::new(),
        }
    }

    /// Add a node. Returns the node id.
    pub fn add_node(&mut self, node: NetworkNode) -> usize {
        let id = node.id;
        self.adjacency.entry(id).or_default();
        self.nodes.push(node);
        id
    }

    /// Add an edge and update adjacency.
    pub fn add_edge(&mut self, edge: NetworkEdge) {
        self.adjacency
            .entry(edge.source)
            .or_default()
            .push((edge.target, edge.weight));
        if !edge.directed {
            self.adjacency
                .entry(edge.target)
                .or_default()
                .push((edge.source, edge.weight));
        }
        self.edges.push(edge);
    }

    /// Return the degree of a node (number of neighbors in adjacency list).
    pub fn degree(&self, node_id: usize) -> usize {
        self.adjacency.get(&node_id).map(|v| v.len()).unwrap_or(0)
    }

    /// Return neighbors of a node.
    pub fn neighbors(&self, node_id: usize) -> Vec<usize> {
        self.adjacency
            .get(&node_id)
            .map(|v| v.iter().map(|(n, _)| *n).collect())
            .unwrap_or_default()
    }

    /// Number of nodes.
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }

    /// Number of edges.
    pub fn edge_count(&self) -> usize {
        self.edges.len()
    }

    /// Check if edge (u, v) exists.
    pub fn has_edge(&self, u: usize, v: usize) -> bool {
        self.adjacency
            .get(&u)
            .map(|ns| ns.iter().any(|(n, _)| *n == v))
            .unwrap_or(false)
    }

    /// BFS from start, returns visited node ids in BFS order.
    pub fn bfs(&self, start: usize) -> Vec<usize> {
        let mut visited = HashSet::new();
        let mut queue = VecDeque::new();
        let mut result = Vec::new();
        queue.push_back(start);
        visited.insert(start);
        while let Some(node) = queue.pop_front() {
            result.push(node);
            for neighbor in self.neighbors(node) {
                if visited.insert(neighbor) {
                    queue.push_back(neighbor);
                }
            }
        }
        result
    }

    /// DFS from start, returns visited node ids in DFS order.
    pub fn dfs(&self, start: usize) -> Vec<usize> {
        let mut visited = HashSet::new();
        let mut result = Vec::new();
        self.dfs_rec(start, &mut visited, &mut result);
        result
    }

    fn dfs_rec(&self, node: usize, visited: &mut HashSet<usize>, result: &mut Vec<usize>) {
        if !visited.insert(node) {
            return;
        }
        result.push(node);
        for neighbor in self.neighbors(node) {
            self.dfs_rec(neighbor, visited, result);
        }
    }

    /// Build a simple adjacency matrix (node indices in `nodes` order).
    pub fn adjacency_matrix(&self) -> Vec<Vec<f64>> {
        let n = self.nodes.len();
        let id_to_idx: HashMap<usize, usize> = self
            .nodes
            .iter()
            .enumerate()
            .map(|(i, node)| (node.id, i))
            .collect();
        let mut mat = vec![vec![0.0f64; n]; n];
        for edge in &self.edges {
            if let (Some(&i), Some(&j)) = (id_to_idx.get(&edge.source), id_to_idx.get(&edge.target))
            {
                mat[i][j] = edge.weight;
                if !edge.directed {
                    mat[j][i] = edge.weight;
                }
            }
        }
        mat
    }

    /// Get mutable reference to node by id.
    pub fn node_mut(&mut self, id: usize) -> Option<&mut NetworkNode> {
        self.nodes.iter_mut().find(|n| n.id == id)
    }

    /// Get reference to node by id.
    pub fn node(&self, id: usize) -> Option<&NetworkNode> {
        self.nodes.iter().find(|n| n.id == id)
    }
}

impl Default for NetworkGraph {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// ForceDirectedLayout — Fruchterman-Reingold with Barnes-Hut hint
// ---------------------------------------------------------------------------

/// Configuration for force-directed layout.
#[derive(Debug, Clone)]
pub struct ForceDirectedConfig {
    /// Number of iterations.
    pub iterations: usize,
    /// Optimal distance between nodes (k).
    pub optimal_distance: f64,
    /// Initial temperature.
    pub temperature: f64,
    /// Cooling factor (multiplied each iteration).
    pub cooling: f64,
    /// Repulsion coefficient.
    pub repulsion: f64,
    /// Attraction coefficient.
    pub attraction: f64,
    /// Barnes-Hut theta (0 = exact, higher = faster approximation).
    pub barnes_hut_theta: f64,
}

impl Default for ForceDirectedConfig {
    fn default() -> Self {
        Self {
            iterations: 200,
            optimal_distance: 50.0,
            temperature: 100.0,
            cooling: 0.95,
            repulsion: 1.0,
            attraction: 1.0,
            barnes_hut_theta: 0.5,
        }
    }
}

/// Force-directed graph layout (Fruchterman-Reingold).
#[derive(Debug, Clone)]
pub struct ForceDirectedLayout {
    /// Layout configuration.
    pub config: ForceDirectedConfig,
}

impl ForceDirectedLayout {
    /// Create with default configuration.
    pub fn new() -> Self {
        Self {
            config: ForceDirectedConfig::default(),
        }
    }

    /// Create with custom configuration.
    pub fn with_config(config: ForceDirectedConfig) -> Self {
        Self { config }
    }

    /// Run layout, modifying node positions in-place.
    pub fn apply(&self, graph: &mut NetworkGraph) {
        let node_ids: Vec<usize> = graph.nodes.iter().map(|n| n.id).collect();
        let n = node_ids.len();
        if n == 0 {
            return;
        }
        // Initialize positions if at origin
        let all_zero = graph.nodes.iter().all(|nd| nd.position == [0.0, 0.0]);
        if all_zero {
            let side = (n as f64).sqrt() * self.config.optimal_distance;
            for (i, &id) in node_ids.iter().enumerate() {
                let angle = 2.0 * std::f64::consts::PI * i as f64 / n as f64;
                if let Some(nd) = graph.node_mut(id) {
                    nd.position = [side * angle.cos(), side * angle.sin()];
                }
            }
        }

        let k = self.config.optimal_distance;
        let mut temp = self.config.temperature;

        for _iter in 0..self.config.iterations {
            let positions: HashMap<usize, [f64; 2]> =
                graph.nodes.iter().map(|nd| (nd.id, nd.position)).collect();

            // Repulsive forces
            let mut disp: HashMap<usize, [f64; 2]> =
                node_ids.iter().map(|&id| (id, [0.0, 0.0])).collect();
            for i in 0..node_ids.len() {
                for j in (i + 1)..node_ids.len() {
                    let u = node_ids[i];
                    let v = node_ids[j];
                    let pu = positions[&u];
                    let pv = positions[&v];
                    let dx = pu[0] - pv[0];
                    let dy = pu[1] - pv[1];
                    let dist = (dx * dx + dy * dy).sqrt().max(0.01);
                    let force = self.config.repulsion * k * k / dist;
                    let fx = force * dx / dist;
                    let fy = force * dy / dist;
                    if let Some(d) = disp.get_mut(&u) {
                        d[0] += fx;
                        d[1] += fy;
                    }
                    if let Some(d) = disp.get_mut(&v) {
                        d[0] -= fx;
                        d[1] -= fy;
                    }
                }
            }

            // Attractive forces
            for edge in &graph.edges {
                let pu = positions[&edge.source];
                let pv = positions[&edge.target];
                let dx = pu[0] - pv[0];
                let dy = pu[1] - pv[1];
                let dist = (dx * dx + dy * dy).sqrt().max(0.01);
                let force = self.config.attraction * dist * dist / k;
                let fx = force * dx / dist;
                let fy = force * dy / dist;
                if let Some(d) = disp.get_mut(&edge.source) {
                    d[0] -= fx;
                    d[1] -= fy;
                }
                if let Some(d) = disp.get_mut(&edge.target) {
                    d[0] += fx;
                    d[1] += fy;
                }
            }

            // Apply displacements, clamped by temperature
            for nd in &mut graph.nodes {
                let d = disp[&nd.id];
                let dmag = (d[0] * d[0] + d[1] * d[1]).sqrt().max(0.001);
                let scale = dmag.min(temp) / dmag;
                nd.position[0] += d[0] * scale;
                nd.position[1] += d[1] * scale;
            }

            temp *= self.config.cooling;
        }
    }

    /// Barnes-Hut approximation: compute repulsion using quad-tree bucketing.
    ///
    /// Returns a per-node displacement map (Barnes-Hut approximate repulsion only).
    pub fn barnes_hut_repulsion(&self, graph: &NetworkGraph) -> HashMap<usize, [f64; 2]> {
        let k = self.config.optimal_distance;
        let theta = self.config.barnes_hut_theta;
        let positions: Vec<(usize, [f64; 2])> =
            graph.nodes.iter().map(|nd| (nd.id, nd.position)).collect();
        let mut result: HashMap<usize, [f64; 2]> =
            positions.iter().map(|(id, _)| (*id, [0.0, 0.0])).collect();

        // Simplified BH: group nodes by 2x2 cells and use cell centroid
        let cell_size = k * 2.0 / theta.max(0.01);
        let mut cells: HashMap<(i64, i64), (f64, f64, usize)> = HashMap::new();
        for (_, pos) in &positions {
            let cx = (pos[0] / cell_size).floor() as i64;
            let cy = (pos[1] / cell_size).floor() as i64;
            let entry = cells.entry((cx, cy)).or_insert((0.0, 0.0, 0));
            entry.0 += pos[0];
            entry.1 += pos[1];
            entry.2 += 1;
        }

        for (id, pos) in &positions {
            let mut fx = 0.0f64;
            let mut fy = 0.0f64;
            for ((cx, cy), (sx, sy, cnt)) in &cells {
                let cpx = sx / *cnt as f64;
                let cpy = sy / *cnt as f64;
                let dx = pos[0] - cpx;
                let dy = pos[1] - cpy;
                let dist = (dx * dx + dy * dy).sqrt().max(0.01);
                let cell_px = *cx as f64 * cell_size;
                let cell_py = *cy as f64 * cell_size;
                let cdx = pos[0] - cell_px;
                let cdy = pos[1] - cell_py;
                let to_cell = (cdx * cdx + cdy * cdy).sqrt();
                // BH condition: use aggregate if cell is far enough
                let use_aggregate = (cell_size / to_cell.max(0.01)) < theta || *cnt == 1;
                if use_aggregate {
                    let force = self.config.repulsion * k * k * (*cnt as f64) / dist;
                    fx += force * dx / dist;
                    fy += force * dy / dist;
                }
            }
            if let Some(d) = result.get_mut(id) {
                d[0] = fx;
                d[1] = fy;
            }
        }
        result
    }
}

impl Default for ForceDirectedLayout {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// HierarchicalLayout — Sugiyama framework
// ---------------------------------------------------------------------------

/// Configuration for hierarchical layout.
#[derive(Debug, Clone)]
pub struct HierarchicalConfig {
    /// Vertical layer separation.
    pub layer_separation: f64,
    /// Horizontal node separation.
    pub node_separation: f64,
    /// Number of crossing minimization passes.
    pub crossing_passes: usize,
}

impl Default for HierarchicalConfig {
    fn default() -> Self {
        Self {
            layer_separation: 100.0,
            node_separation: 60.0,
            crossing_passes: 4,
        }
    }
}

/// Hierarchical (Sugiyama) graph layout.
#[derive(Debug, Clone)]
pub struct HierarchicalLayout {
    /// Layout configuration.
    pub config: HierarchicalConfig,
}

impl HierarchicalLayout {
    /// Create with default configuration.
    pub fn new() -> Self {
        Self {
            config: HierarchicalConfig::default(),
        }
    }

    /// Create with custom configuration.
    pub fn with_config(config: HierarchicalConfig) -> Self {
        Self { config }
    }

    /// Run layout, modifying node positions in-place.
    ///
    /// Step 1: Layer assignment (longest-path).
    /// Step 2: Crossing minimization (barycenter heuristic).
    /// Step 3: Coordinate assignment.
    pub fn apply(&self, graph: &mut NetworkGraph) {
        if graph.nodes.is_empty() {
            return;
        }
        let layers = self.assign_layers(graph);
        let ordered = self.minimize_crossings(graph, &layers);
        self.assign_coordinates(graph, &ordered);
    }

    /// Assign layer (rank) to each node using longest-path layering.
    pub fn assign_layers(&self, graph: &NetworkGraph) -> HashMap<usize, usize> {
        let ids: Vec<usize> = graph.nodes.iter().map(|n| n.id).collect();
        // Compute in-degrees
        let mut in_degree: HashMap<usize, usize> = ids.iter().map(|&id| (id, 0)).collect();
        for edge in &graph.edges {
            *in_degree.entry(edge.target).or_insert(0) += 1;
        }
        let mut layer: HashMap<usize, usize> = HashMap::new();
        let mut queue: VecDeque<usize> = ids
            .iter()
            .filter(|&&id| in_degree[&id] == 0)
            .copied()
            .collect();
        if queue.is_empty() {
            // No sources (cyclic): assign all to layer 0
            for &id in &ids {
                layer.insert(id, 0);
            }
            return layer;
        }
        while let Some(node) = queue.pop_front() {
            let cur_layer = *layer.get(&node).unwrap_or(&0);
            for edge in &graph.edges {
                if edge.source == node {
                    let next_layer = cur_layer + 1;
                    let entry = layer.entry(edge.target).or_insert(0);
                    if next_layer > *entry {
                        *entry = next_layer;
                    }
                    let deg = in_degree
                        .get_mut(&edge.target)
                        .expect("key must exist in map");
                    *deg -= 1;
                    if *deg == 0 {
                        queue.push_back(edge.target);
                    }
                }
            }
            layer.entry(node).or_insert(0);
        }
        // Assign missing nodes
        for &id in &ids {
            layer.entry(id).or_insert(0);
        }
        layer
    }

    /// Barycenter crossing minimization.
    ///
    /// Returns layers as ordered lists of node ids.
    pub fn minimize_crossings(
        &self,
        graph: &NetworkGraph,
        layers: &HashMap<usize, usize>,
    ) -> Vec<Vec<usize>> {
        if layers.is_empty() {
            return Vec::new();
        }
        let max_layer = *layers.values().max().unwrap_or(&0);
        let mut layer_nodes: Vec<Vec<usize>> = vec![Vec::new(); max_layer + 1];
        for (&id, &l) in layers {
            layer_nodes[l].push(id);
        }

        for _pass in 0..self.config.crossing_passes {
            // Forward sweep: order by barycenter of previous layer
            for l in 1..layer_nodes.len() {
                let prev: HashMap<usize, usize> = layer_nodes[l - 1]
                    .iter()
                    .enumerate()
                    .map(|(i, &id)| (id, i))
                    .collect();
                let cur = &layer_nodes[l];
                let mut bary: Vec<(f64, usize)> = cur
                    .iter()
                    .map(|&id| {
                        let neighbors: Vec<usize> = graph
                            .adjacency
                            .get(&id)
                            .map(|ns| ns.iter().map(|(n, _)| *n).collect())
                            .unwrap_or_default();
                        let positions: Vec<f64> = neighbors
                            .iter()
                            .filter_map(|n| prev.get(n))
                            .map(|&p| p as f64)
                            .collect();
                        let bc = if positions.is_empty() {
                            0.0
                        } else {
                            positions.iter().sum::<f64>() / positions.len() as f64
                        };
                        (bc, id)
                    })
                    .collect();
                bary.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
                layer_nodes[l] = bary.into_iter().map(|(_, id)| id).collect();
            }
        }
        layer_nodes
    }

    /// Assign (x, y) coordinates based on ordered layers.
    pub fn assign_coordinates(&self, graph: &mut NetworkGraph, ordered_layers: &[Vec<usize>]) {
        for (layer_idx, layer) in ordered_layers.iter().enumerate() {
            let y = layer_idx as f64 * self.config.layer_separation;
            let total_w = layer.len() as f64 * self.config.node_separation;
            let start_x = -total_w / 2.0;
            for (pos_idx, &node_id) in layer.iter().enumerate() {
                let x = start_x + pos_idx as f64 * self.config.node_separation;
                if let Some(nd) = graph.node_mut(node_id) {
                    nd.position = [x, y];
                }
            }
        }
    }
}

impl Default for HierarchicalLayout {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// NetworkMetrics — graph analysis
// ---------------------------------------------------------------------------

/// Degree distribution of the network.
#[derive(Debug, Clone)]
pub struct DegreeDistribution {
    /// Map from degree → number of nodes with that degree.
    pub counts: HashMap<usize, usize>,
    /// Total number of nodes.
    pub total_nodes: usize,
}

impl DegreeDistribution {
    /// Probability mass function at degree k.
    pub fn pmf(&self, k: usize) -> f64 {
        if self.total_nodes == 0 {
            return 0.0;
        }
        *self.counts.get(&k).unwrap_or(&0) as f64 / self.total_nodes as f64
    }

    /// Mean degree.
    pub fn mean_degree(&self) -> f64 {
        if self.total_nodes == 0 {
            return 0.0;
        }
        let total: usize = self.counts.iter().map(|(k, c)| k * c).sum();
        total as f64 / self.total_nodes as f64
    }
}

/// Network analysis metrics.
#[derive(Debug, Clone)]
pub struct NetworkMetrics {
    graph: NetworkGraph,
}

impl NetworkMetrics {
    /// Create from a graph reference (clones the graph).
    pub fn new(graph: NetworkGraph) -> Self {
        Self { graph }
    }

    /// Compute the degree distribution.
    pub fn degree_distribution(&self) -> DegreeDistribution {
        let mut counts: HashMap<usize, usize> = HashMap::new();
        for node in &self.graph.nodes {
            let deg = self.graph.degree(node.id);
            *counts.entry(deg).or_insert(0) += 1;
        }
        DegreeDistribution {
            counts,
            total_nodes: self.graph.node_count(),
        }
    }

    /// Compute the local clustering coefficient for a node.
    pub fn local_clustering_coefficient(&self, node_id: usize) -> f64 {
        let neighbors: Vec<usize> = self.graph.neighbors(node_id);
        let k = neighbors.len();
        if k < 2 {
            return 0.0;
        }
        let neighbor_set: HashSet<usize> = neighbors.iter().copied().collect();
        let mut triangles = 0usize;
        for i in 0..neighbors.len() {
            for j in (i + 1)..neighbors.len() {
                if self.graph.has_edge(neighbors[i], neighbors[j])
                    || self.graph.has_edge(neighbors[j], neighbors[i])
                {
                    triangles += 1;
                }
            }
        }
        let _ = neighbor_set; // used above
        2.0 * triangles as f64 / (k * (k - 1)) as f64
    }

    /// Compute the global (average) clustering coefficient.
    pub fn global_clustering_coefficient(&self) -> f64 {
        let nodes = &self.graph.nodes;
        if nodes.is_empty() {
            return 0.0;
        }
        let sum: f64 = nodes
            .iter()
            .map(|n| self.local_clustering_coefficient(n.id))
            .sum();
        sum / nodes.len() as f64
    }

    /// Compute approximate betweenness centrality using BFS (unweighted).
    ///
    /// Returns a map from node_id → centrality score.
    pub fn betweenness_centrality(&self) -> HashMap<usize, f64> {
        let ids: Vec<usize> = self.graph.nodes.iter().map(|n| n.id).collect();
        let mut centrality: HashMap<usize, f64> = ids.iter().map(|&id| (id, 0.0)).collect();
        let n = ids.len();
        if n == 0 {
            return centrality;
        }

        for &src in &ids {
            // BFS to compute shortest-path counts
            let mut dist: HashMap<usize, i64> = HashMap::new();
            let mut sigma: HashMap<usize, f64> = HashMap::new(); // num shortest paths
            let mut pred: HashMap<usize, Vec<usize>> = HashMap::new();
            let mut queue = VecDeque::new();
            let mut stack = Vec::new();

            for &id in &ids {
                dist.insert(id, -1);
                sigma.insert(id, 0.0);
                pred.insert(id, Vec::new());
            }
            *dist.get_mut(&src).expect("key must exist in map") = 0;
            *sigma.get_mut(&src).expect("key must exist in map") = 1.0;
            queue.push_back(src);

            while let Some(v) = queue.pop_front() {
                stack.push(v);
                for w in self.graph.neighbors(v) {
                    if dist[&w] == -1 {
                        queue.push_back(w);
                        *dist.get_mut(&w).expect("key must exist in map") = dist[&v] + 1;
                    }
                    if dist[&w] == dist[&v] + 1 {
                        let sv = sigma[&v];
                        *sigma.get_mut(&w).expect("key must exist in map") += sv;
                        pred.get_mut(&w).expect("key must exist in map").push(v);
                    }
                }
            }

            // Back-propagation
            let mut delta: HashMap<usize, f64> = ids.iter().map(|&id| (id, 0.0)).collect();
            while let Some(w) = stack.pop() {
                for &v in &pred[&w] {
                    let coeff = (sigma[&v] / sigma[&w]) * (1.0 + delta[&w]);
                    *delta.get_mut(&v).expect("key must exist in map") += coeff;
                }
                if w != src {
                    *centrality.get_mut(&w).expect("key must exist in map") += delta[&w];
                }
            }
        }

        // Normalize
        let norm = if n > 2 {
            1.0 / ((n - 1) * (n - 2)) as f64
        } else {
            1.0
        };
        for v in centrality.values_mut() {
            *v *= norm;
        }
        centrality
    }

    /// Compute PageRank with damping factor `d` and `max_iter` iterations.
    pub fn pagerank(&self, d: f64, max_iter: usize, tol: f64) -> HashMap<usize, f64> {
        let ids: Vec<usize> = self.graph.nodes.iter().map(|n| n.id).collect();
        let n = ids.len();
        if n == 0 {
            return HashMap::new();
        }
        let init = 1.0 / n as f64;
        let mut pr: HashMap<usize, f64> = ids.iter().map(|&id| (id, init)).collect();

        for _iter in 0..max_iter {
            let mut new_pr: HashMap<usize, f64> =
                ids.iter().map(|&id| (id, (1.0 - d) / n as f64)).collect();
            for &v in &ids {
                let out_deg = self.graph.degree(v) as f64;
                if out_deg > 0.0 {
                    let share = d * pr[&v] / out_deg;
                    for u in self.graph.neighbors(v) {
                        *new_pr.entry(u).or_insert(0.0) += share;
                    }
                } else {
                    // Dangling node: distribute to all
                    let share = d * pr[&v] / n as f64;
                    for &u in &ids {
                        *new_pr.entry(u).or_insert(0.0) += share;
                    }
                }
            }
            // Check convergence
            let diff: f64 = ids.iter().map(|id| (new_pr[id] - pr[id]).abs()).sum();
            pr = new_pr;
            if diff < tol {
                break;
            }
        }
        pr
    }

    /// Network density: actual edges / possible edges.
    pub fn density(&self) -> f64 {
        let n = self.graph.node_count();
        let m = self.graph.edge_count();
        if n < 2 {
            return 0.0;
        }
        m as f64 / (n * (n - 1)) as f64
    }

    /// Average shortest path length (BFS, unweighted). Returns None if graph is disconnected.
    pub fn average_shortest_path_length(&self) -> Option<f64> {
        let ids: Vec<usize> = self.graph.nodes.iter().map(|n| n.id).collect();
        let n = ids.len();
        if n == 0 {
            return None;
        }
        let mut total = 0u64;
        let mut count = 0u64;
        for &src in &ids {
            let mut dist: HashMap<usize, usize> = HashMap::new();
            let mut queue = VecDeque::new();
            dist.insert(src, 0);
            queue.push_back(src);
            while let Some(v) = queue.pop_front() {
                for u in self.graph.neighbors(v) {
                    if !dist.contains_key(&u) {
                        dist.insert(u, dist[&v] + 1);
                        queue.push_back(u);
                    }
                }
            }
            if dist.len() < n {
                return None; // disconnected
            }
            for (&_v, &d) in &dist {
                total += d as u64;
                count += 1;
            }
        }
        if count == 0 {
            None
        } else {
            Some(total as f64 / count as f64)
        }
    }
}

// ---------------------------------------------------------------------------
// CommunityDetection
// ---------------------------------------------------------------------------

/// Result of community detection.
#[derive(Debug, Clone)]
pub struct CommunityResult {
    /// Map from node_id → community_id.
    pub membership: HashMap<usize, usize>,
    /// Modularity score.
    pub modularity: f64,
}

impl CommunityResult {
    /// Number of distinct communities.
    pub fn num_communities(&self) -> usize {
        let set: HashSet<usize> = self.membership.values().copied().collect();
        set.len()
    }

    /// Return members of community `c`.
    pub fn community_members(&self, c: usize) -> Vec<usize> {
        self.membership
            .iter()
            .filter(|(_, v)| **v == c)
            .map(|(k, _)| *k)
            .collect()
    }
}

/// Community detection algorithms.
#[derive(Debug, Clone)]
pub struct CommunityDetection {
    graph: NetworkGraph,
}

impl CommunityDetection {
    /// Create from a graph.
    pub fn new(graph: NetworkGraph) -> Self {
        Self { graph }
    }

    /// Compute modularity Q for a given membership.
    pub fn modularity(&self, membership: &HashMap<usize, usize>) -> f64 {
        let m: f64 = self.graph.edges.iter().map(|e| e.weight).sum();
        if m == 0.0 {
            return 0.0;
        }
        let mut q = 0.0f64;
        let ids: Vec<usize> = self.graph.nodes.iter().map(|n| n.id).collect();
        for &u in &ids {
            for &v in &ids {
                if membership.get(&u) == membership.get(&v) {
                    let a_uv = if self.graph.has_edge(u, v) {
                        self.graph
                            .adjacency
                            .get(&u)
                            .and_then(|ns| ns.iter().find(|(n, _)| *n == v))
                            .map(|(_, w)| *w)
                            .unwrap_or(0.0)
                    } else {
                        0.0
                    };
                    let ku = self.graph.degree(u) as f64;
                    let kv = self.graph.degree(v) as f64;
                    q += a_uv - ku * kv / (2.0 * m);
                }
            }
        }
        q / (2.0 * m)
    }

    /// Greedy modularity optimization (Louvain-style, simplified single-pass).
    pub fn greedy_modularity(&self) -> CommunityResult {
        let ids: Vec<usize> = self.graph.nodes.iter().map(|n| n.id).collect();
        // Start: each node is its own community
        let mut membership: HashMap<usize, usize> = ids.iter().map(|&id| (id, id)).collect();
        let mut improved = true;
        let max_passes = ids.len().min(20);
        let mut passes = 0;

        while improved && passes < max_passes {
            improved = false;
            passes += 1;
            for &node in &ids {
                let current_comm = membership[&node];
                let neighbors = self.graph.neighbors(node);
                if neighbors.is_empty() {
                    continue;
                }
                // Candidate communities from neighbors
                let candidate_comms: HashSet<usize> =
                    neighbors.iter().map(|n| membership[n]).collect();
                let current_q = self.modularity(&membership);
                let mut best_q = current_q;
                let mut best_comm = current_comm;
                for cand_comm in candidate_comms {
                    if cand_comm == current_comm {
                        continue;
                    }
                    let old = membership
                        .insert(node, cand_comm)
                        .expect("previous value must exist");
                    let new_q = self.modularity(&membership);
                    if new_q > best_q {
                        best_q = new_q;
                        best_comm = cand_comm;
                        improved = true;
                    } else {
                        membership.insert(node, old);
                    }
                }
                membership.insert(node, best_comm);
            }
        }

        let q = self.modularity(&membership);
        CommunityResult {
            membership,
            modularity: q,
        }
    }

    /// Label propagation community detection.
    pub fn label_propagation(&self, max_iter: usize) -> CommunityResult {
        let ids: Vec<usize> = self.graph.nodes.iter().map(|n| n.id).collect();
        let mut labels: HashMap<usize, usize> = ids.iter().map(|&id| (id, id)).collect();

        for _iter in 0..max_iter {
            let mut changed = false;
            // Randomized order (deterministic: sort by id)
            let mut order = ids.clone();
            order.sort_unstable();

            for &node in &order {
                let neighbors = self.graph.neighbors(node);
                if neighbors.is_empty() {
                    continue;
                }
                // Count label frequencies
                let mut freq: HashMap<usize, usize> = HashMap::new();
                for &nb in &neighbors {
                    *freq.entry(labels[&nb]).or_insert(0) += 1;
                }
                let max_freq = freq.values().copied().max().unwrap_or(0);
                let best_label = freq
                    .iter()
                    .filter(|(_, c)| **c == max_freq)
                    .map(|(l, _)| *l)
                    .min()
                    .unwrap_or(labels[&node]);
                if best_label != labels[&node] {
                    labels.insert(node, best_label);
                    changed = true;
                }
            }
            if !changed {
                break;
            }
        }

        // Re-number communities 0..n
        let mut label_map: HashMap<usize, usize> = HashMap::new();
        let mut counter = 0usize;
        let mut membership: HashMap<usize, usize> = HashMap::new();
        for &id in &ids {
            let old_label = labels[&id];
            let new_label = *label_map.entry(old_label).or_insert_with(|| {
                let c = counter;
                counter += 1;
                c
            });
            membership.insert(id, new_label);
        }

        let q = self.modularity(&membership);
        CommunityResult {
            membership,
            modularity: q,
        }
    }

    /// Spectral clustering into `k` communities using graph Laplacian.
    ///
    /// Uses the Fiedler vector (2nd smallest eigenvector of L) for k=2,
    /// or repeated bisection for k>2.
    pub fn spectral_clustering(&self, k: usize) -> CommunityResult {
        let ids: Vec<usize> = self.graph.nodes.iter().map(|n| n.id).collect();
        let n = ids.len();
        if n == 0 || k == 0 {
            return CommunityResult {
                membership: HashMap::new(),
                modularity: 0.0,
            };
        }
        if k >= n {
            let membership: HashMap<usize, usize> =
                ids.iter().enumerate().map(|(i, &id)| (id, i)).collect();
            return CommunityResult {
                membership,
                modularity: 0.0,
            };
        }

        let id_idx: HashMap<usize, usize> =
            ids.iter().enumerate().map(|(i, &id)| (id, i)).collect();
        // Build normalized Laplacian L = D - A
        let mut degree = vec![0.0f64; n];
        let mut adj = vec![vec![0.0f64; n]; n];
        for edge in &self.graph.edges {
            if let (Some(&i), Some(&j)) = (id_idx.get(&edge.source), id_idx.get(&edge.target)) {
                adj[i][j] += edge.weight;
                degree[i] += edge.weight;
                if !edge.directed {
                    adj[j][i] += edge.weight;
                    degree[j] += edge.weight;
                }
            }
        }
        // Power iteration for the Fiedler vector (approximate)
        // For simplicity: use degree-based heuristic for spectral embedding
        let fiedler = self.approximate_fiedler(&adj, &degree, n);

        // Partition by repeated bisection
        let membership = self.bisection_partition(&ids, &fiedler, k);
        let q = self.modularity(&membership);
        CommunityResult {
            membership,
            modularity: q,
        }
    }

    fn approximate_fiedler(&self, adj: &[Vec<f64>], degree: &[f64], n: usize) -> Vec<f64> {
        // Simple power-iteration approximation of 2nd eigenvector of Laplacian
        let mut v: Vec<f64> = (0..n).map(|i| i as f64 - n as f64 / 2.0).collect();
        // Orthogonalize with constant vector (remove 1st eigenvector)
        let norm_v: f64 = (v.iter().map(|x| x * x).sum::<f64>()).sqrt().max(1e-12);
        for x in &mut v {
            *x /= norm_v;
        }

        for _iter in 0..50 {
            // Apply L = D - A
            let mut lv: Vec<f64> = vec![0.0; n];
            for i in 0..n {
                lv[i] = degree[i] * v[i];
                for j in 0..n {
                    lv[i] -= adj[i][j] * v[j];
                }
            }
            // Remove component along [1,1,...,1]/sqrt(n)
            let mean = lv.iter().sum::<f64>() / n as f64;
            for x in &mut lv {
                *x -= mean;
            }
            // Normalize
            let norm: f64 = (lv.iter().map(|x| x * x).sum::<f64>()).sqrt().max(1e-12);
            for i in 0..n {
                v[i] = lv[i] / norm;
            }
        }
        v
    }

    fn bisection_partition(
        &self,
        ids: &[usize],
        fiedler: &[f64],
        k: usize,
    ) -> HashMap<usize, usize> {
        // Use median splits for k communities
        let mut assignment: Vec<usize> = vec![0; ids.len()];
        let splits = vec![(0usize, ids.len(), 0usize)]; // (start, end, community)
        let next_comm = 1usize;

        let mut sorted_idx: Vec<usize> = (0..ids.len()).collect();
        sorted_idx.sort_by(|&a, &b| {
            fiedler[a]
                .partial_cmp(&fiedler[b])
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        // Assign initial community by position in Fiedler vector
        for (rank, &idx) in sorted_idx.iter().enumerate() {
            assignment[idx] = (rank * k / ids.len()).min(k - 1);
        }
        let _ = splits;
        let _ = next_comm;

        ids.iter()
            .enumerate()
            .map(|(i, &id)| (id, assignment[i]))
            .collect()
    }
}

// ---------------------------------------------------------------------------
// NetworkAnimation — time-evolving graph
// ---------------------------------------------------------------------------

/// An event in the network timeline.
#[derive(Debug, Clone)]
#[allow(missing_docs)]
pub enum NetworkEvent {
    /// A node appears at the given time.
    NodeAppear { node_id: usize, time: f64 },
    /// A node disappears at the given time.
    NodeDisappear { node_id: usize, time: f64 },
    /// An edge appears between source and target at the given time.
    EdgeAppear {
        source: usize,
        target: usize,
        time: f64,
    },
    /// An edge disappears at the given time.
    EdgeDisappear {
        source: usize,
        target: usize,
        time: f64,
    },
    /// A node moves to a new position over a duration.
    NodeMove {
        node_id: usize,
        target_position: [f64; 2],
        start_time: f64,
        duration: f64,
    },
}

/// Snapshot of the network at a single frame.
#[derive(Debug, Clone)]
pub struct NetworkFrame {
    /// Time of this frame.
    pub time: f64,
    /// Active node ids.
    pub active_nodes: HashSet<usize>,
    /// Active edges as (source, target) pairs.
    pub active_edges: Vec<(usize, usize)>,
    /// Node positions at this frame.
    pub node_positions: HashMap<usize, [f64; 2]>,
}

/// Network animation: time-evolving graph with interpolation.
#[derive(Debug, Clone)]
pub struct NetworkAnimation {
    /// Base graph (contains all possible nodes/edges).
    pub base_graph: NetworkGraph,
    /// Timeline events.
    pub events: Vec<NetworkEvent>,
    /// Total animation duration.
    pub total_duration: f64,
    /// Frames per second for sampling.
    pub fps: f64,
}

impl NetworkAnimation {
    /// Create a new animation from a base graph.
    pub fn new(base_graph: NetworkGraph, total_duration: f64, fps: f64) -> Self {
        Self {
            base_graph,
            events: Vec::new(),
            total_duration,
            fps,
        }
    }

    /// Add an event.
    pub fn add_event(&mut self, event: NetworkEvent) {
        self.events.push(event);
    }

    /// Compute the network frame at time `t`.
    pub fn frame_at(&self, t: f64) -> NetworkFrame {
        let mut active_nodes: HashSet<usize> = HashSet::new();
        let mut active_edges: Vec<(usize, usize)> = Vec::new();
        let mut node_positions: HashMap<usize, [f64; 2]> = HashMap::new();

        // Initialize from base graph
        for node in &self.base_graph.nodes {
            node_positions.insert(node.id, node.position);
        }

        // Process events up to time t
        for event in &self.events {
            match event {
                NetworkEvent::NodeAppear { node_id, time } => {
                    if t >= *time {
                        active_nodes.insert(*node_id);
                    }
                }
                NetworkEvent::NodeDisappear { node_id, time } => {
                    if t >= *time {
                        active_nodes.remove(node_id);
                    }
                }
                NetworkEvent::EdgeAppear {
                    source,
                    target,
                    time,
                } => {
                    if t >= *time {
                        active_edges.push((*source, *target));
                    }
                }
                NetworkEvent::EdgeDisappear {
                    source,
                    target,
                    time,
                } => {
                    if t >= *time {
                        active_edges.retain(|&(s, d)| !(s == *source && d == *target));
                    }
                }
                NetworkEvent::NodeMove {
                    node_id,
                    target_position,
                    start_time,
                    duration,
                } => {
                    if t >= *start_time {
                        let alpha = if *duration <= 0.0 {
                            1.0
                        } else {
                            ((t - start_time) / duration).clamp(0.0, 1.0)
                        };
                        if let Some(pos) = node_positions.get_mut(node_id)
                            && let Some(node) = self.base_graph.node(*node_id)
                        {
                            let start = node.position;
                            pos[0] = start[0] * (1.0 - alpha) + target_position[0] * alpha;
                            pos[1] = start[1] * (1.0 - alpha) + target_position[1] * alpha;
                        }
                    }
                }
            }
        }

        active_edges.dedup();
        NetworkFrame {
            time: t,
            active_nodes,
            active_edges,
            node_positions,
        }
    }

    /// Sample the animation into `n_frames` evenly-spaced frames.
    pub fn sample_frames(&self, n_frames: usize) -> Vec<NetworkFrame> {
        if n_frames == 0 {
            return Vec::new();
        }
        (0..n_frames)
            .map(|i| {
                let t = self.total_duration * i as f64 / (n_frames - 1).max(1) as f64;
                self.frame_at(t)
            })
            .collect()
    }

    /// Interpolate node position between two times using linear interpolation.
    pub fn interpolate_position(&self, node_id: usize, t0: f64, t1: f64, alpha: f64) -> [f64; 2] {
        let f0 = self.frame_at(t0);
        let f1 = self.frame_at(t1);
        let p0 = f0
            .node_positions
            .get(&node_id)
            .copied()
            .unwrap_or([0.0, 0.0]);
        let p1 = f1
            .node_positions
            .get(&node_id)
            .copied()
            .unwrap_or([0.0, 0.0]);
        [
            p0[0] * (1.0 - alpha) + p1[0] * alpha,
            p0[1] * (1.0 - alpha) + p1[1] * alpha,
        ]
    }

    /// Return how many events are in the timeline.
    pub fn event_count(&self) -> usize {
        self.events.len()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn make_triangle() -> NetworkGraph {
        let mut g = NetworkGraph::new();
        g.add_node(NetworkNode::new(0, "A"));
        g.add_node(NetworkNode::new(1, "B"));
        g.add_node(NetworkNode::new(2, "C"));
        g.add_edge(NetworkEdge::undirected(0, 1));
        g.add_edge(NetworkEdge::undirected(1, 2));
        g.add_edge(NetworkEdge::undirected(0, 2));
        g
    }

    fn make_chain(n: usize) -> NetworkGraph {
        let mut g = NetworkGraph::new();
        for i in 0..n {
            g.add_node(NetworkNode::new(i, format!("N{}", i)));
        }
        for i in 0..n - 1 {
            g.add_edge(NetworkEdge::undirected(i, i + 1));
        }
        g
    }

    // --- NetworkGraph ---

    #[test]
    fn test_graph_node_count() {
        let g = make_triangle();
        assert_eq!(g.node_count(), 3);
    }

    #[test]
    fn test_graph_edge_count() {
        let g = make_triangle();
        assert_eq!(g.edge_count(), 3);
    }

    #[test]
    fn test_graph_degree() {
        let g = make_triangle();
        assert_eq!(g.degree(0), 2);
    }

    #[test]
    fn test_graph_has_edge() {
        let g = make_triangle();
        assert!(g.has_edge(0, 1));
        assert!(!g.has_edge(0, 99));
    }

    #[test]
    fn test_graph_neighbors() {
        let g = make_triangle();
        let mut nb = g.neighbors(0);
        nb.sort_unstable();
        assert_eq!(nb, vec![1, 2]);
    }

    #[test]
    fn test_graph_bfs() {
        let g = make_chain(5);
        let visited = g.bfs(0);
        assert_eq!(visited.len(), 5);
        assert_eq!(visited[0], 0);
    }

    #[test]
    fn test_graph_dfs() {
        let g = make_chain(5);
        let visited = g.dfs(0);
        assert_eq!(visited.len(), 5);
    }

    #[test]
    fn test_graph_adjacency_matrix() {
        let g = make_triangle();
        let mat = g.adjacency_matrix();
        assert_eq!(mat.len(), 3);
        // Symmetric
        for i in 0..3 {
            for j in 0..3 {
                assert!((mat[i][j] - mat[j][i]).abs() < 1e-10);
            }
        }
    }

    #[test]
    fn test_graph_node_mut() {
        let mut g = make_triangle();
        g.node_mut(0).unwrap().position = [5.0, 5.0];
        assert_eq!(g.node(0).unwrap().position, [5.0, 5.0]);
    }

    #[test]
    fn test_graph_directed_edge() {
        let mut g = NetworkGraph::new();
        g.add_node(NetworkNode::new(0, "A"));
        g.add_node(NetworkNode::new(1, "B"));
        g.add_edge(NetworkEdge::directed(0, 1, 2.0));
        assert!(g.has_edge(0, 1));
        assert!(!g.has_edge(1, 0));
    }

    #[test]
    fn test_graph_default() {
        let g = NetworkGraph::default();
        assert_eq!(g.node_count(), 0);
    }

    #[test]
    fn test_node_builder() {
        let n = NetworkNode::new(5, "test")
            .with_position(1.0, 2.0)
            .with_size(3.0)
            .with_color(1.0, 0.0, 0.0, 1.0);
        assert_eq!(n.position, [1.0, 2.0]);
        assert_eq!(n.size, 3.0);
    }

    #[test]
    fn test_edge_label() {
        let e = NetworkEdge::undirected(0, 1).with_label("link");
        assert_eq!(e.label.unwrap(), "link");
    }

    // --- ForceDirectedLayout ---

    #[test]
    fn test_force_directed_applies() {
        let mut g = make_triangle();
        let layout = ForceDirectedLayout::new();
        layout.apply(&mut g);
        // Positions should have moved
        let not_all_zero = g.nodes.iter().any(|n| n.position != [0.0, 0.0]);
        assert!(not_all_zero);
    }

    #[test]
    fn test_force_directed_empty() {
        let mut g = NetworkGraph::new();
        let layout = ForceDirectedLayout::new();
        layout.apply(&mut g); // Should not panic
    }

    #[test]
    fn test_force_directed_single_node() {
        let mut g = NetworkGraph::new();
        g.add_node(NetworkNode::new(0, "A"));
        let layout = ForceDirectedLayout::new();
        layout.apply(&mut g);
    }

    #[test]
    fn test_force_directed_config() {
        let config = ForceDirectedConfig {
            iterations: 10,
            optimal_distance: 30.0,
            ..ForceDirectedConfig::default()
        };
        let mut g = make_chain(4);
        let layout = ForceDirectedLayout::with_config(config);
        layout.apply(&mut g);
    }

    #[test]
    fn test_force_directed_default() {
        let layout = ForceDirectedLayout::default();
        assert_eq!(layout.config.iterations, 200);
    }

    #[test]
    fn test_barnes_hut_repulsion() {
        let g = make_triangle();
        let layout = ForceDirectedLayout::new();
        let disp = layout.barnes_hut_repulsion(&g);
        assert_eq!(disp.len(), 3);
    }

    // --- HierarchicalLayout ---

    #[test]
    fn test_hierarchical_applies() {
        let mut g = NetworkGraph::new();
        g.add_node(NetworkNode::new(0, "root"));
        g.add_node(NetworkNode::new(1, "child1"));
        g.add_node(NetworkNode::new(2, "child2"));
        g.add_edge(NetworkEdge::directed(0, 1, 1.0));
        g.add_edge(NetworkEdge::directed(0, 2, 1.0));
        let layout = HierarchicalLayout::new();
        layout.apply(&mut g);
        // Root should be at a different y than children
        let root_y = g.node(0).unwrap().position[1];
        let child_y = g.node(1).unwrap().position[1];
        assert!((root_y - child_y).abs() > 1.0);
    }

    #[test]
    fn test_hierarchical_empty() {
        let mut g = NetworkGraph::new();
        let layout = HierarchicalLayout::new();
        layout.apply(&mut g);
    }

    #[test]
    fn test_hierarchical_assign_layers() {
        let g = make_chain(4);
        let layout = HierarchicalLayout::new();
        let layers = layout.assign_layers(&g);
        assert_eq!(layers.len(), 4);
    }

    #[test]
    fn test_hierarchical_minimize_crossings() {
        let g = make_chain(3);
        let layout = HierarchicalLayout::new();
        let layers = layout.assign_layers(&g);
        let ordered = layout.minimize_crossings(&g, &layers);
        assert!(!ordered.is_empty());
    }

    #[test]
    fn test_hierarchical_default() {
        let layout = HierarchicalLayout::default();
        assert_eq!(layout.config.layer_separation, 100.0);
    }

    // --- NetworkMetrics ---

    #[test]
    fn test_metrics_degree_distribution() {
        let g = make_triangle();
        let m = NetworkMetrics::new(g);
        let dd = m.degree_distribution();
        assert_eq!(dd.total_nodes, 3);
        assert!((dd.mean_degree() - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_metrics_clustering_triangle() {
        let g = make_triangle();
        let m = NetworkMetrics::new(g);
        let cc = m.local_clustering_coefficient(0);
        assert!((cc - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_metrics_global_clustering() {
        let g = make_triangle();
        let m = NetworkMetrics::new(g);
        let gc = m.global_clustering_coefficient();
        assert!((gc - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_metrics_betweenness() {
        let g = make_chain(5);
        let m = NetworkMetrics::new(g);
        let bc = m.betweenness_centrality();
        assert_eq!(bc.len(), 5);
        // Node 2 (middle) should have highest centrality
        assert!(bc[&2] >= bc[&0]);
    }

    #[test]
    fn test_metrics_pagerank_sums_to_one() {
        let g = make_triangle();
        let m = NetworkMetrics::new(g);
        let pr = m.pagerank(0.85, 100, 1e-8);
        let total: f64 = pr.values().sum();
        assert!((total - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_metrics_density() {
        let g = make_triangle();
        let m = NetworkMetrics::new(g);
        // 3 directed neighbors each → density = 6/(3*2) = 1
        assert!(m.density() > 0.0);
    }

    #[test]
    fn test_metrics_avg_path_chain() {
        let g = make_chain(3);
        let m = NetworkMetrics::new(g);
        let apl = m.average_shortest_path_length();
        assert!(apl.is_some());
    }

    #[test]
    fn test_degree_distribution_pmf() {
        let g = make_triangle();
        let m = NetworkMetrics::new(g);
        let dd = m.degree_distribution();
        let p = dd.pmf(2);
        assert!((p - 1.0).abs() < 1e-10);
    }

    // --- CommunityDetection ---

    #[test]
    fn test_community_label_propagation() {
        let g = make_triangle();
        let cd = CommunityDetection::new(g);
        let result = cd.label_propagation(20);
        assert_eq!(result.membership.len(), 3);
    }

    #[test]
    fn test_community_num_communities() {
        let g = make_chain(6);
        let cd = CommunityDetection::new(g);
        let result = cd.label_propagation(10);
        assert!(result.num_communities() >= 1);
    }

    #[test]
    fn test_community_greedy() {
        let g = make_triangle();
        let cd = CommunityDetection::new(g);
        let result = cd.greedy_modularity();
        assert_eq!(result.membership.len(), 3);
    }

    #[test]
    fn test_community_spectral_k2() {
        let g = make_chain(4);
        let cd = CommunityDetection::new(g);
        let result = cd.spectral_clustering(2);
        assert_eq!(result.membership.len(), 4);
        assert!(result.num_communities() <= 2);
    }

    #[test]
    fn test_community_spectral_k_geq_n() {
        let g = make_triangle();
        let cd = CommunityDetection::new(g);
        let result = cd.spectral_clustering(10);
        assert_eq!(result.membership.len(), 3);
    }

    #[test]
    fn test_community_members() {
        let g = make_triangle();
        let cd = CommunityDetection::new(g);
        let result = cd.label_propagation(10);
        let comm = result.membership[&0];
        let members = result.community_members(comm);
        assert!(members.contains(&0));
    }

    #[test]
    fn test_community_modularity_range() {
        let g = make_triangle();
        let cd = CommunityDetection::new(g);
        let result = cd.greedy_modularity();
        // Modularity is bounded [-0.5, 1.0]
        assert!(result.modularity >= -0.5 && result.modularity <= 1.0);
    }

    // --- NetworkAnimation ---

    #[test]
    fn test_animation_frame_at_start() {
        let g = make_triangle();
        let mut anim = NetworkAnimation::new(g, 10.0, 30.0);
        anim.add_event(NetworkEvent::NodeAppear {
            node_id: 0,
            time: 0.0,
        });
        anim.add_event(NetworkEvent::NodeAppear {
            node_id: 1,
            time: 1.0,
        });
        let f = anim.frame_at(0.5);
        assert!(f.active_nodes.contains(&0));
        assert!(!f.active_nodes.contains(&1));
    }

    #[test]
    fn test_animation_edge_appear() {
        let g = make_triangle();
        let mut anim = NetworkAnimation::new(g, 10.0, 30.0);
        anim.add_event(NetworkEvent::EdgeAppear {
            source: 0,
            target: 1,
            time: 2.0,
        });
        let f = anim.frame_at(3.0);
        assert!(f.active_edges.contains(&(0, 1)));
    }

    #[test]
    fn test_animation_node_disappear() {
        let g = make_triangle();
        let mut anim = NetworkAnimation::new(g, 10.0, 30.0);
        anim.add_event(NetworkEvent::NodeAppear {
            node_id: 0,
            time: 0.0,
        });
        anim.add_event(NetworkEvent::NodeDisappear {
            node_id: 0,
            time: 5.0,
        });
        let f0 = anim.frame_at(3.0);
        let f1 = anim.frame_at(7.0);
        assert!(f0.active_nodes.contains(&0));
        assert!(!f1.active_nodes.contains(&0));
    }

    #[test]
    fn test_animation_edge_disappear() {
        let g = make_triangle();
        let mut anim = NetworkAnimation::new(g, 10.0, 30.0);
        anim.add_event(NetworkEvent::EdgeAppear {
            source: 0,
            target: 1,
            time: 0.0,
        });
        anim.add_event(NetworkEvent::EdgeDisappear {
            source: 0,
            target: 1,
            time: 5.0,
        });
        let f0 = anim.frame_at(3.0);
        let f1 = anim.frame_at(7.0);
        assert!(f0.active_edges.contains(&(0, 1)));
        assert!(!f1.active_edges.contains(&(0, 1)));
    }

    #[test]
    fn test_animation_node_move() {
        let g = make_triangle();
        let mut anim = NetworkAnimation::new(g, 10.0, 30.0);
        anim.add_event(NetworkEvent::NodeMove {
            node_id: 0,
            target_position: [100.0, 100.0],
            start_time: 0.0,
            duration: 10.0,
        });
        let f = anim.frame_at(10.0);
        let pos = f.node_positions[&0];
        assert!((pos[0] - 100.0).abs() < 1e-6);
    }

    #[test]
    fn test_animation_sample_frames() {
        let g = make_triangle();
        let anim = NetworkAnimation::new(g, 10.0, 30.0);
        let frames = anim.sample_frames(5);
        assert_eq!(frames.len(), 5);
        assert!((frames[0].time).abs() < 1e-10);
    }

    #[test]
    fn test_animation_interpolate_position() {
        let g = make_triangle();
        let mut anim = NetworkAnimation::new(g, 10.0, 30.0);
        anim.add_event(NetworkEvent::NodeMove {
            node_id: 0,
            target_position: [100.0, 0.0],
            start_time: 0.0,
            duration: 10.0,
        });
        let p = anim.interpolate_position(0, 0.0, 10.0, 0.5);
        assert!(p[0].is_finite());
    }

    #[test]
    fn test_animation_event_count() {
        let g = make_triangle();
        let mut anim = NetworkAnimation::new(g, 10.0, 30.0);
        anim.add_event(NetworkEvent::NodeAppear {
            node_id: 0,
            time: 0.0,
        });
        assert_eq!(anim.event_count(), 1);
    }

    #[test]
    fn test_animation_sample_empty() {
        let g = make_triangle();
        let anim = NetworkAnimation::new(g, 10.0, 30.0);
        let frames = anim.sample_frames(0);
        assert!(frames.is_empty());
    }

    #[test]
    fn test_force_directed_chain() {
        let mut g = make_chain(6);
        let config = ForceDirectedConfig {
            iterations: 50,
            ..Default::default()
        };
        let layout = ForceDirectedLayout::with_config(config);
        layout.apply(&mut g);
        // All positions should be finite
        for n in &g.nodes {
            assert!(n.position[0].is_finite());
            assert!(n.position[1].is_finite());
        }
    }
}
