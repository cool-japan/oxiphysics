// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Unit tests of the graph model shared by the community-detection kernels:
//! how a [`NetworkGraph`] is indexed, and the `O(n + m)` modularity.

use std::collections::HashMap;

use super::super::{NetworkEdge, NetworkGraph, NetworkNode};
use super::test_support::{Lcg, dense_modularity, graph_from_edges, karate_club};
use super::{
    IndexedNetwork, WeightedGraph, densify, greedy_modularity, modularity, spectral_clustering,
};

/// Network with the given node ids (in that order) and undirected edges.
fn network(ids: &[usize], edges: &[(usize, usize, f64)]) -> NetworkGraph {
    let mut graph = NetworkGraph::new();
    for &id in ids {
        graph.add_node(NetworkNode::new(id, format!("n{id}")));
    }
    for &(a, b, weight) in edges {
        let mut edge = NetworkEdge::undirected(a, b);
        edge.weight = weight;
        graph.add_edge(edge);
    }
    graph
}

// ---------------------------------------------------------------------------
// Indexing
// ---------------------------------------------------------------------------

#[test]
fn indexing_orders_nodes_by_first_occurrence_and_merges_duplicates() {
    let graph = network(&[40, 7, 40, 19, 7], &[(7, 19, 1.0), (40, 7, 2.0)]);
    let indexed = IndexedNetwork::new(&graph);
    assert_eq!(indexed.ids, vec![40, 7, 19]);
    assert_eq!(indexed.graph.len(), 3);
    // Edge (7, 19) -> indices (1, 2); edge (40, 7) -> indices (0, 1).
    assert_eq!(indexed.graph.strength, vec![2.0, 3.0, 1.0]);
    assert_eq!(indexed.graph.total_weight, 3.0);
}

#[test]
fn indexing_skips_edges_that_cannot_take_part() {
    let graph = network(
        &[0, 1, 2],
        &[
            (0, 1, 1.0),
            (1, 2, 2.0),
            (0, 99, 5.0),          // unknown target
            (98, 2, 5.0),          // unknown source
            (0, 2, f64::NAN),      // not a number
            (0, 2, f64::INFINITY), // not finite
            (0, 2, f64::NEG_INFINITY),
            (0, 2, -3.0), // negative
            (0, 2, 0.0),  // zero
            (1, 1, -1.0), // negative self-loop
        ],
    );
    let indexed = IndexedNetwork::new(&graph);
    assert_eq!(indexed.graph.total_weight, 3.0);
    assert_eq!(indexed.graph.strength, vec![1.0, 3.0, 2.0]);
    assert_eq!(indexed.graph.self_loop, vec![0.0; 3]);
}

#[test]
fn indexing_symmetrises_edges_and_counts_self_loops_twice() {
    let mut graph = network(&[10, 20, 30], &[]);
    graph.add_edge(NetworkEdge::directed(10, 20, 2.0));
    graph.add_edge(NetworkEdge::directed(20, 10, 3.0)); // opposite direction
    graph.add_edge(NetworkEdge::undirected(20, 30)); // unit weight
    let mut parallel = NetworkEdge::undirected(20, 30);
    parallel.weight = 0.5;
    graph.add_edge(parallel);
    let mut loop_edge = NetworkEdge::undirected(30, 30);
    loop_edge.weight = 4.0;
    graph.add_edge(loop_edge);
    let g = IndexedNetwork::new(&graph).graph;

    assert_eq!(g.total_weight, 2.0 + 3.0 + 1.0 + 0.5 + 4.0);
    // Strengths: node 10: 2 + 3; node 20: 2 + 3 + 1 + 0.5; node 30: 1 + 0.5 + 2*4.
    assert_eq!(g.strength, vec![5.0, 6.5, 9.5]);
    assert_eq!(g.self_loop, vec![0.0, 0.0, 4.0]);
    // Every edge is listed at both end points; the loop is not an adjacency.
    assert_eq!(g.adj[0], vec![(1, 2.0), (1, 3.0)]);
    assert_eq!(g.adj[1], vec![(0, 2.0), (0, 3.0), (2, 1.0), (2, 0.5)]);
    assert_eq!(g.adj[2], vec![(1, 1.0), (1, 0.5)]);
}

#[test]
fn strengths_add_up_to_twice_the_total_weight() {
    for seed in 1..=8 {
        let mut rng = Lcg(seed * 977);
        let mut graph = WeightedGraph::with_nodes(25);
        for _ in 0..60 {
            let a = rng.below(25);
            let b = rng.below(25);
            graph.add_edge(a, b, 0.1 + 5.0 * rng.unit());
        }
        let strength_sum: f64 = graph.strength.iter().sum();
        assert!(
            (strength_sum - 2.0 * graph.total_weight).abs() < 1e-9,
            "seed {seed}"
        );
        // Adjacency entries plus doubled loops account for the same weight.
        let listed: f64 = graph.adj.iter().flatten().map(|&(_, w)| w).sum::<f64>()
            + 2.0 * graph.self_loop.iter().sum::<f64>();
        assert!((listed - strength_sum).abs() < 1e-9, "seed {seed}");
    }
}

#[test]
fn labels_make_unmentioned_nodes_singletons() {
    let graph = network(&[5, 6, 7, 8, 9], &[]);
    let indexed = IndexedNetwork::new(&graph);
    let membership: HashMap<usize, usize> = [(5, 100), (7, 100), (9, 3), (1234, 8)]
        .into_iter()
        .collect();
    let labels = indexed.labels_from(&membership);
    // 5 and 7 share a label, 9 has its own, 6 and 8 are unmentioned and so
    // each a community of its own; the unknown id 1234 is ignored.
    assert_eq!(labels[0], labels[2]);
    let mut distinct = labels.clone();
    distinct.sort_unstable();
    distinct.dedup();
    assert_eq!(distinct.len(), 4);
    assert_ne!(labels[1], labels[3]);
    assert!(
        labels.iter().all(|&label| label < 5),
        "labels are dense: {labels:?}"
    );
}

#[test]
fn result_lists_every_distinct_node_once_with_its_modularity() {
    let graph = network(&[3, 1, 3, 2], &[(3, 1, 1.0), (1, 2, 1.0)]);
    let indexed = IndexedNetwork::new(&graph);
    let result = indexed.result(&[0, 0, 1]);
    assert_eq!(result.membership.len(), 3);
    assert_eq!(result.membership[&3], 0);
    assert_eq!(result.membership[&1], 0);
    assert_eq!(result.membership[&2], 1);
    // m = 2, strengths 1, 2, 1; {3, 1} | {2}: Q = 1/2 - (3/4)^2 + 0 - (1/4)^2.
    let expected = 0.5 - 0.5625 - 0.0625;
    assert!((result.modularity - expected).abs() < 1e-12);
}

#[test]
fn densify_numbers_labels_by_first_appearance() {
    let mut labels = vec![42, 7, 42, 99, 7, 0, 99];
    let count = densify(&mut labels);
    assert_eq!(count, 4);
    assert_eq!(labels, vec![0, 1, 0, 2, 1, 3, 2]);
    let mut empty: Vec<usize> = Vec::new();
    assert_eq!(densify(&mut empty), 0);
}

// ---------------------------------------------------------------------------
// Modularity
// ---------------------------------------------------------------------------

#[test]
fn modularity_equals_the_dense_double_sum() {
    for seed in 1..=12 {
        let mut rng = Lcg(seed * 7919);
        let n = 26;
        let mut graph = WeightedGraph::with_nodes(n);
        for _ in 0..80 {
            let a = rng.below(n);
            let b = if rng.below(8) == 0 { a } else { rng.below(n) };
            graph.add_edge(a, b, 0.1 + 4.0 * rng.unit());
        }
        let mut labels: Vec<usize> = (0..n).map(|_| rng.below(5)).collect();
        densify(&mut labels);
        let fast = graph.modularity(&labels);
        let dense = dense_modularity(&graph, &labels);
        assert!(
            (fast - dense).abs() < 1e-12,
            "seed {seed}: {fast} vs {dense}"
        );
    }
}

#[test]
fn modularity_of_trivial_partitions_and_degenerate_graphs() {
    let graph = karate_club();
    // One community: L = m and D = 2m, so Q = 0.
    assert!(graph.modularity(&vec![0; 34]).abs() < 1e-12);
    // Singletons: Q = -sum (k_i / 2m)^2 < 0.
    let singletons: Vec<usize> = (0..34).collect();
    let expected: f64 = -graph
        .strength
        .iter()
        .map(|k| (k / 156.0).powi(2))
        .sum::<f64>();
    assert!((graph.modularity(&singletons) - expected).abs() < 1e-12);
    // No nodes, no edges, only self-loops of zero total: Q = 0 and no panic.
    assert_eq!(WeightedGraph::with_nodes(0).modularity(&[]), 0.0);
    assert_eq!(WeightedGraph::with_nodes(4).modularity(&[0, 0, 1, 1]), 0.0);
}

#[test]
fn modularity_is_zero_instead_of_not_a_number_when_the_total_weight_overflows() {
    let huge = f64::MAX / 2.0;
    let graph = graph_from_edges(3, &[(0, 1, huge), (1, 2, huge), (0, 2, huge)]);
    assert!(graph.total_weight.is_infinite());
    assert_eq!(graph.modularity(&[0, 0, 1]), 0.0);
}

#[test]
fn modularity_of_two_communities_matches_the_closed_form() {
    // Two triangles with weights 2 joined by a bridge of weight 1: 11/26.
    let graph = graph_from_edges(
        6,
        &[
            (0, 1, 2.0),
            (0, 2, 2.0),
            (1, 2, 2.0),
            (3, 4, 2.0),
            (3, 5, 2.0),
            (4, 5, 2.0),
            (2, 3, 1.0),
        ],
    );
    assert!((graph.modularity(&[0, 0, 0, 1, 1, 1]) - 11.0 / 26.0).abs() < 1e-12);
}

// ---------------------------------------------------------------------------
// Entry points
// ---------------------------------------------------------------------------

#[test]
fn entry_points_agree_with_each_other() {
    let graph = network(
        &(0..12).collect::<Vec<_>>(),
        &[
            (0, 1, 1.0),
            (0, 2, 1.0),
            (1, 2, 1.0),
            (3, 4, 1.0),
            (3, 5, 1.0),
            (4, 5, 1.0),
            (2, 3, 0.2),
            (6, 7, 1.0),
            (7, 8, 1.0),
            (6, 8, 1.0),
            (9, 10, 1.0),
            (10, 11, 1.0),
            (9, 11, 1.0),
            (8, 9, 0.2),
            (5, 6, 0.2),
        ],
    );
    for result in [greedy_modularity(&graph), spectral_clustering(&graph, 4)] {
        assert_eq!(result.membership.len(), 12);
        let recomputed = modularity(&graph, &result.membership);
        assert!((result.modularity - recomputed).abs() < 1e-12);
    }
    // Four triangles in a chain: both methods find them.
    let triangles = greedy_modularity(&graph);
    assert_eq!(triangles.num_communities(), 4);
    let spectral = spectral_clustering(&graph, 4);
    assert_eq!(spectral.num_communities(), 4);
    assert!((triangles.modularity - spectral.modularity).abs() < 1e-12);
}

#[test]
fn spectral_entry_point_handles_empty_requests() {
    let graph = network(&[0, 1, 2], &[(0, 1, 1.0)]);
    let none = spectral_clustering(&graph, 0);
    assert!(none.membership.is_empty() && none.modularity == 0.0);
    let empty = spectral_clustering(&NetworkGraph::new(), 3);
    assert!(empty.membership.is_empty() && empty.modularity == 0.0);
    // k == n: singletons with their real modularity.
    let singles = spectral_clustering(&graph, 3);
    assert_eq!(singles.membership.len(), 3);
    assert!((singles.modularity - modularity(&graph, &singles.membership)).abs() < 1e-12);
}
