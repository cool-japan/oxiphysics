// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Kernels behind [`CommunityDetection`](super::CommunityDetection).
//!
//! [`NetworkGraph`] keys its nodes by arbitrary `usize` ids and keeps the
//! adjacency in hash maps.  Every algorithm in this module first maps the graph
//! onto the dense indices `0..n` ([`IndexedNetwork`]) so that the hot loops run
//! on plain vectors and never on hash lookups or the linear `has_edge` scans of
//! the graph type.
//!
//! # Graph model
//!
//! Community detection works on the *undirected, weighted* graph underlying a
//! [`NetworkGraph`]:
//!
//! * An edge contributes its weight to both of its end points, whether or not
//!   it is marked `directed`.
//! * Only edges whose two end points are nodes of the graph and whose weight is
//!   finite and strictly positive take part.  An edge to a node that does not
//!   exist is not part of the graph, and a zero, negative or non-finite weight
//!   cannot be expressed in Newman's modularity; such edges are skipped.
//! * Parallel edges add up.  A self-loop of weight `w` adds `w` to the internal
//!   weight of its community and `2w` to the strength of its node.
//!
//! With `k_i` the strength (weighted degree) of node `i` and `m` the total edge
//! weight (`2m = sum_i k_i`), the modularity of a partition is
//!
//! ```text
//! Q = sum_c [ L_c / m - (D_c / 2m)^2 ]
//! ```
//!
//! where `L_c` is the weight of the edges inside community `c` and `D_c` the
//! summed strength of its nodes.  This is Newman's
//! `1/(2m) * sum_ij [A_ij - k_i k_j / 2m] * delta(c_i, c_j)` evaluated in
//! `O(n + m)` instead of `O(n^2)`.

use std::collections::HashMap;

use super::{CommunityResult, NetworkGraph};

mod louvain;
mod spectral;

#[cfg(test)]
mod test_support;
#[cfg(test)]
mod tests;

/// Symmetric weighted graph on the dense node indices `0..len`.
///
/// Louvain aggregates communities into super-nodes, so the same type also
/// describes the coarsened graphs of the later levels.
#[derive(Debug, Clone)]
struct WeightedGraph {
    /// `adj[i]` lists `(j, w)` for every edge `{i, j}` with `j != i`.  Each edge
    /// appears in the lists of both end points; parallel edges are separate
    /// entries (or merged, on coarsened graphs).
    adj: Vec<Vec<(usize, f64)>>,
    /// Total weight of the self-loops at each node, every loop counted once.
    self_loop: Vec<f64>,
    /// Weighted degree `k_i`: the weights of all incident edges, a self-loop
    /// counted twice.
    strength: Vec<f64>,
    /// `m`, the sum of all edge weights; equals `sum_i k_i / 2`.
    total_weight: f64,
}

impl WeightedGraph {
    /// Graph on `n` nodes without edges.
    fn with_nodes(n: usize) -> Self {
        Self {
            adj: vec![Vec::new(); n],
            self_loop: vec![0.0; n],
            strength: vec![0.0; n],
            total_weight: 0.0,
        }
    }

    /// Number of nodes.
    fn len(&self) -> usize {
        self.adj.len()
    }

    /// Add the undirected edge `{a, b}` of weight `weight`, which the caller
    /// guarantees to be finite and positive.
    fn add_edge(&mut self, a: usize, b: usize, weight: f64) {
        if a == b {
            self.self_loop[a] += weight;
            self.strength[a] += 2.0 * weight;
        } else {
            self.adj[a].push((b, weight));
            self.adj[b].push((a, weight));
            self.strength[a] += weight;
            self.strength[b] += weight;
        }
        self.total_weight += weight;
    }

    /// Modularity `Q = sum_c [L_c / m - (D_c / 2m)^2]` of the partition that
    /// puts node `i` into community `community[i]`.
    ///
    /// The labels must be dense (each below `community.len()`).  Runs in
    /// `O(n + m)`.  A graph without usable edge weight has `Q = 0`.
    fn modularity(&self, community: &[usize]) -> f64 {
        let m = self.total_weight;
        if !(m.is_finite() && m > 0.0) {
            return 0.0;
        }
        let groups = community.iter().copied().max().map_or(0, |label| label + 1);
        let mut internal = vec![0.0_f64; groups];
        let mut volume = vec![0.0_f64; groups];
        for (node, &label) in community.iter().enumerate() {
            volume[label] += self.strength[node];
            internal[label] += self.self_loop[node];
            for &(other, weight) in &self.adj[node] {
                if community[other] == label {
                    // An internal edge is listed at both of its end points.
                    internal[label] += 0.5 * weight;
                }
            }
        }
        let two_m = 2.0 * m;
        internal
            .iter()
            .zip(&volume)
            .map(|(&inside, &total)| inside / m - (total / two_m).powi(2))
            .sum()
    }
}

/// A [`NetworkGraph`] mapped onto dense node indices.
#[derive(Debug, Clone)]
struct IndexedNetwork {
    /// Node id of every index, in order of first occurrence in
    /// `network.nodes` (a node listed twice is one node).
    ids: Vec<usize>,
    /// The usable edges of the network, symmetrised.
    graph: WeightedGraph,
}

impl IndexedNetwork {
    /// Index `network`, keeping only the edges described in the module docs.
    fn new(network: &NetworkGraph) -> Self {
        let mut ids: Vec<usize> = Vec::with_capacity(network.nodes.len());
        let mut index: HashMap<usize, usize> = HashMap::with_capacity(network.nodes.len());
        for node in &network.nodes {
            index.entry(node.id).or_insert_with(|| {
                ids.push(node.id);
                ids.len() - 1
            });
        }
        let mut graph = WeightedGraph::with_nodes(ids.len());
        for edge in &network.edges {
            let usable = edge.weight.is_finite() && edge.weight > 0.0;
            if !usable {
                continue;
            }
            if let (Some(&a), Some(&b)) = (index.get(&edge.source), index.get(&edge.target)) {
                graph.add_edge(a, b, edge.weight);
            }
        }
        Self { ids, graph }
    }

    /// Dense community labels for `membership`.  A node that `membership` does
    /// not mention is a community of its own.
    fn labels_from(&self, membership: &HashMap<usize, usize>) -> Vec<usize> {
        let mut dense: HashMap<usize, usize> = HashMap::new();
        let mut next = 0_usize;
        let mut labels = Vec::with_capacity(self.ids.len());
        for id in &self.ids {
            let label = match membership.get(id) {
                Some(given) => *dense.entry(*given).or_insert_with(|| {
                    next += 1;
                    next - 1
                }),
                None => {
                    next += 1;
                    next - 1
                }
            };
            labels.push(label);
        }
        labels
    }

    /// Package dense `labels` (one per index) as a [`CommunityResult`].
    fn result(&self, labels: &[usize]) -> CommunityResult {
        CommunityResult {
            membership: self
                .ids
                .iter()
                .copied()
                .zip(labels.iter().copied())
                .collect(),
            modularity: self.graph.modularity(labels),
        }
    }
}

/// Relabel `labels` in place to `0..count` in order of first appearance and
/// return `count`.
fn densify(labels: &mut [usize]) -> usize {
    let mut dense: HashMap<usize, usize> = HashMap::new();
    for label in labels.iter_mut() {
        let next = dense.len();
        *label = *dense.entry(*label).or_insert(next);
    }
    dense.len()
}

/// Modularity of `membership` on `network` (see the module docs for the
/// conventions).
pub(super) fn modularity(network: &NetworkGraph, membership: &HashMap<usize, usize>) -> f64 {
    let indexed = IndexedNetwork::new(network);
    let labels = indexed.labels_from(membership);
    indexed.graph.modularity(&labels)
}

/// Louvain modularity maximisation on `network`.
pub(super) fn greedy_modularity(network: &NetworkGraph) -> CommunityResult {
    let indexed = IndexedNetwork::new(network);
    let labels = louvain::louvain(&indexed.graph);
    indexed.result(&labels)
}

/// Spectral partition of `network` into `k` communities by recursive
/// bisection.  No nodes or `k == 0` give an empty result.
pub(super) fn spectral_clustering(network: &NetworkGraph, k: usize) -> CommunityResult {
    let indexed = IndexedNetwork::new(network);
    if indexed.ids.is_empty() || k == 0 {
        return CommunityResult {
            membership: HashMap::new(),
            modularity: 0.0,
        };
    }
    let labels = spectral::spectral_partition(&indexed.graph, k);
    indexed.result(&labels)
}
