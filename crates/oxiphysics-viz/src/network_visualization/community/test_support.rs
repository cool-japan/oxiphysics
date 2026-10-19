// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Graph builders and reference implementations shared by the unit tests of
//! the community-detection kernels.

use super::WeightedGraph;

/// Deterministic 64-bit LCG (Knuth's MMIX constants) for reproducible graphs.
pub(super) struct Lcg(pub(super) u64);

impl Lcg {
    /// Next raw 64-bit value.
    pub(super) fn next_u64(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.0
    }

    /// Uniform integer in `0..bound`.
    pub(super) fn below(&mut self, bound: usize) -> usize {
        ((self.next_u64() >> 33) as usize) % bound
    }

    /// Uniform float in `[0, 1)`.
    pub(super) fn unit(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
}

/// Graph on `n` nodes with the given `(a, b, weight)` edges.
pub(super) fn graph_from_edges(n: usize, edges: &[(usize, usize, f64)]) -> WeightedGraph {
    let mut graph = WeightedGraph::with_nodes(n);
    for &(a, b, weight) in edges {
        graph.add_edge(a, b, weight);
    }
    graph
}

/// Connect all `members` pairwise with edges of the given weight.
pub(super) fn clique_edges(members: &[usize], weight: f64) -> Vec<(usize, usize, f64)> {
    let mut edges = Vec::new();
    for (i, &a) in members.iter().enumerate() {
        for &b in &members[i + 1..] {
            edges.push((a, b, weight));
        }
    }
    edges
}

/// Path `0 - 1 - ... - (n - 1)` with unit weights.
pub(super) fn path_graph(n: usize) -> WeightedGraph {
    let edges: Vec<(usize, usize, f64)> =
        (0..n.saturating_sub(1)).map(|i| (i, i + 1, 1.0)).collect();
    graph_from_edges(n, &edges)
}

/// Zachary's karate club: 34 members, 78 unit-weight friendships.
pub(super) fn karate_club() -> WeightedGraph {
    const EDGES: [(usize, usize); 78] = [
        (0, 1),
        (0, 2),
        (0, 3),
        (0, 4),
        (0, 5),
        (0, 6),
        (0, 7),
        (0, 8),
        (0, 10),
        (0, 11),
        (0, 12),
        (0, 13),
        (0, 17),
        (0, 19),
        (0, 21),
        (0, 31),
        (1, 2),
        (1, 3),
        (1, 7),
        (1, 13),
        (1, 17),
        (1, 19),
        (1, 21),
        (1, 30),
        (2, 3),
        (2, 7),
        (2, 8),
        (2, 9),
        (2, 13),
        (2, 27),
        (2, 28),
        (2, 32),
        (3, 7),
        (3, 12),
        (3, 13),
        (4, 6),
        (4, 10),
        (5, 6),
        (5, 10),
        (5, 16),
        (6, 16),
        (8, 30),
        (8, 32),
        (8, 33),
        (9, 33),
        (13, 33),
        (14, 32),
        (14, 33),
        (15, 32),
        (15, 33),
        (18, 32),
        (18, 33),
        (19, 33),
        (20, 32),
        (20, 33),
        (22, 32),
        (22, 33),
        (23, 25),
        (23, 27),
        (23, 29),
        (23, 32),
        (23, 33),
        (24, 25),
        (24, 27),
        (24, 31),
        (25, 31),
        (26, 29),
        (26, 33),
        (27, 33),
        (28, 31),
        (28, 33),
        (29, 32),
        (29, 33),
        (30, 32),
        (30, 33),
        (31, 32),
        (31, 33),
        (32, 33),
    ];
    let edges: Vec<(usize, usize, f64)> = EDGES.iter().map(|&(a, b)| (a, b, 1.0)).collect();
    graph_from_edges(34, &edges)
}

/// Random connected graph: a random spanning tree plus `extra` random edges,
/// all with weights in `[0.5, 3.5)`.
pub(super) fn random_connected_graph(n: usize, extra: usize, seed: u64) -> WeightedGraph {
    let mut rng = Lcg(seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1);
    let mut edges = Vec::new();
    for node in 1..n {
        let parent = rng.below(node);
        edges.push((node, parent, 0.5 + 3.0 * rng.unit()));
    }
    for _ in 0..extra {
        let a = rng.below(n);
        let b = rng.below(n);
        if a != b {
            edges.push((a, b, 0.5 + 3.0 * rng.unit()));
        }
    }
    graph_from_edges(n, &edges)
}

/// `groups` planted groups of `size` nodes (group `g` owns `g * size ..
/// (g + 1) * size`): every node gets `degree` random edges inside its group,
/// and `bridges` random edges join every pair of groups.
pub(super) fn planted_partition(
    groups: usize,
    size: usize,
    degree: usize,
    bridges: usize,
    seed: u64,
) -> WeightedGraph {
    let mut rng = Lcg(seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1);
    let mut graph = WeightedGraph::with_nodes(groups * size);
    for group in 0..groups {
        for offset in 0..size {
            for _ in 0..degree {
                let mut other = rng.below(size);
                if other == offset {
                    other = (other + 1) % size;
                }
                graph.add_edge(group * size + offset, group * size + other, 1.0);
            }
        }
    }
    for g1 in 0..groups {
        for g2 in g1 + 1..groups {
            for _ in 0..bridges {
                let a = g1 * size + rng.below(size);
                let b = g2 * size + rng.below(size);
                graph.add_edge(a, b, 1.0);
            }
        }
    }
    graph
}

/// Reference modularity as a dense double sum over the adjacency lists of
/// `graph`: `1/(2m) * sum_ij [A_ij - k_i k_j / (2m)] * delta(c_i, c_j)` with
/// `A_ii = 2 * (self-loop weight)`.
pub(super) fn dense_modularity(graph: &WeightedGraph, community: &[usize]) -> f64 {
    let n = graph.len();
    let mut a = vec![vec![0.0_f64; n]; n];
    for (node, neighbours) in graph.adj.iter().enumerate() {
        for &(other, weight) in neighbours {
            a[node][other] += weight;
        }
        a[node][node] += 2.0 * graph.self_loop[node];
    }
    let strength: Vec<f64> = a.iter().map(|row| row.iter().sum()).collect();
    let two_m: f64 = strength.iter().sum();
    if two_m <= 0.0 {
        return 0.0;
    }
    let mut q = 0.0;
    for i in 0..n {
        for j in 0..n {
            if community[i] == community[j] {
                q += a[i][j] - strength[i] * strength[j] / two_m;
            }
        }
    }
    q / two_m
}
