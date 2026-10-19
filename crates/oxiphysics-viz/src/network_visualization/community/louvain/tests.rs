// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Unit tests of the Louvain kernels: the local-moving phase (including a
//! brute-force check that it ends in a local optimum), the aggregation step and
//! the multi-level driver.

use super::super::WeightedGraph;
use super::super::densify;
use super::super::test_support::{
    Lcg, clique_edges, dense_modularity, graph_from_edges, karate_club, planted_partition,
    random_connected_graph,
};
use super::{aggregate, local_moving, louvain};

/// Random weighted graph with parallel edges and self-loops on `n` nodes.
fn messy_graph(n: usize, edges: usize, seed: u64) -> WeightedGraph {
    let mut rng = Lcg(seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1);
    let mut graph = WeightedGraph::with_nodes(n);
    for _ in 0..edges {
        let a = rng.below(n);
        let b = if rng.below(10) == 0 { a } else { rng.below(n) };
        graph.add_edge(a, b, 0.25 + 3.0 * rng.unit());
    }
    graph
}

/// Number of distinct labels in a dense labelling.
fn label_count(labels: &[usize]) -> usize {
    labels.iter().copied().max().map_or(0, |m| m + 1)
}

// ---------------------------------------------------------------------------
// Local moving
// ---------------------------------------------------------------------------

#[test]
fn local_moving_merges_a_clique_into_one_community() {
    let graph = graph_from_edges(6, &clique_edges(&[0, 1, 2, 3, 4, 5], 1.0));
    let step = local_moving(&graph);
    assert!(step.improved);
    assert_eq!(step.count, 1);
    assert_eq!(step.community, vec![0; 6]);
}

#[test]
fn local_moving_leaves_graphs_without_edges_between_nodes_alone() {
    // No edges at all.
    let step = local_moving(&WeightedGraph::with_nodes(5));
    assert!(!step.improved);
    assert_eq!(step.count, 5);
    assert_eq!(step.community, vec![0, 1, 2, 3, 4]);
    // Self-loops only: positive total weight, but no neighbour to join.
    let mut loops = WeightedGraph::with_nodes(3);
    loops.add_edge(0, 0, 1.0);
    loops.add_edge(2, 2, 2.0);
    let step = local_moving(&loops);
    assert!(!step.improved);
    assert_eq!(step.community, vec![0, 1, 2]);
    // The driver returns the singletons for both.
    assert_eq!(louvain(&WeightedGraph::with_nodes(3)), vec![0, 1, 2]);
    assert_eq!(louvain(&loops), vec![0, 1, 2]);
    assert_eq!(louvain(&WeightedGraph::with_nodes(0)), Vec::<usize>::new());
}

#[test]
fn local_moving_never_lowers_modularity() {
    for seed in 1..=10 {
        let graph = messy_graph(30, 70, seed);
        let singletons: Vec<usize> = (0..30).collect();
        let step = local_moving(&graph);
        assert!(
            dense_modularity(&graph, &step.community)
                >= dense_modularity(&graph, &singletons) - 1e-12,
            "seed {seed}"
        );
        assert_eq!(step.count, label_count(&step.community));
    }
}

#[test]
fn local_moving_ends_in_a_local_optimum() {
    // After convergence no single node can raise modularity by moving to
    // another community (or into a fresh one) - checked by brute force, which
    // validates both the incremental gain formula and the stopping rule.
    for seed in 1..=8 {
        let graph = messy_graph(36, 90, seed);
        let step = local_moving(&graph);
        let base = dense_modularity(&graph, &step.community);
        for node in 0..graph.len() {
            for target in 0..=step.count {
                let mut moved = step.community.clone();
                moved[node] = target;
                let q = dense_modularity(&graph, &moved);
                assert!(
                    q <= base + 1e-9,
                    "seed {seed}: moving node {node} to {target} raises Q from {base} to {q}"
                );
            }
        }
    }
}

#[test]
fn local_moving_follows_edge_weights() {
    // Node 3 touches community {0,1,2} through a weight-5 edge and community
    // {4,5,6} through a weight-1 edge; both communities are K3s.
    let mut edges = clique_edges(&[0, 1, 2], 4.0);
    edges.extend(clique_edges(&[4, 5, 6], 4.0));
    edges.push((3, 0, 5.0));
    edges.push((3, 4, 1.0));
    let graph = graph_from_edges(7, &edges);
    let step = local_moving(&graph);
    let labels = &step.community;
    assert_eq!(labels[3], labels[0]);
    assert_ne!(labels[3], labels[4]);
    assert_eq!(labels[0], labels[1]);
    assert_eq!(labels[4], labels[5]);
}

// ---------------------------------------------------------------------------
// Aggregation
// ---------------------------------------------------------------------------

#[test]
fn aggregation_folds_internal_edges_into_self_loops() {
    // Two weight-1 triangles joined by a weight-2 bridge between nodes 2, 3.
    let mut edges = clique_edges(&[0, 1, 2], 1.0);
    edges.extend(clique_edges(&[3, 4, 5], 1.0));
    edges.push((2, 3, 2.0));
    let graph = graph_from_edges(6, &edges);
    let coarse = aggregate(&graph, &[0, 0, 0, 1, 1, 1], 2);
    assert_eq!(coarse.len(), 2);
    assert_eq!(coarse.total_weight, 8.0);
    assert_eq!(coarse.self_loop, vec![3.0, 3.0]);
    assert_eq!(coarse.strength, vec![8.0, 8.0]);
    assert_eq!(coarse.adj, vec![vec![(1, 2.0)], vec![(0, 2.0)]]);
}

#[test]
fn aggregation_merges_parallel_edges_between_groups() {
    // Three edges join group 0 = {0, 1} and group 1 = {2, 3}; node 4 is
    // alone in group 2 and attached to node 1.
    let graph = graph_from_edges(
        5,
        &[
            (0, 2, 1.0),
            (0, 3, 2.0),
            (1, 2, 0.5),
            (0, 1, 4.0),
            (2, 3, 1.5),
            (1, 4, 3.0),
        ],
    );
    let coarse = aggregate(&graph, &[0, 0, 1, 1, 2], 3);
    assert_eq!(coarse.self_loop, vec![4.0, 1.5, 0.0]);
    let weight = |a: usize, b: usize| -> f64 {
        coarse.adj[a]
            .iter()
            .filter(|&&(other, _)| other == b)
            .map(|&(_, w)| w)
            .sum()
    };
    for a in 0..3 {
        assert_eq!(coarse.adj[a].iter().filter(|&&(o, _)| o == a).count(), 0);
        for b in 0..3 {
            let entries = coarse.adj[a].iter().filter(|&&(o, _)| o == b).count();
            assert!(entries <= 1, "group {a} lists {b} {entries} times");
            assert!((weight(a, b) - weight(b, a)).abs() < 1e-12);
        }
    }
    assert!((weight(0, 1) - 3.5).abs() < 1e-12);
    assert!((weight(0, 2) - 3.0).abs() < 1e-12);
    assert_eq!(weight(1, 2), 0.0);
}

#[test]
fn aggregation_preserves_weights_strengths_and_modularity() {
    // The invariant that makes multi-level Louvain exact: the modularity of
    // any partition of the coarse graph equals the modularity of the induced
    // partition of the fine graph.
    for seed in 1..=10 {
        let n = 28;
        let graph = messy_graph(n, 80, seed);
        let mut rng = Lcg(seed + 100);
        let mut fine_labels: Vec<usize> = (0..n).map(|_| rng.below(6)).collect();
        let count = densify(&mut fine_labels);
        let coarse = aggregate(&graph, &fine_labels, count);

        assert!((coarse.total_weight - graph.total_weight).abs() < 1e-12);
        let strength_sum = |g: &WeightedGraph| g.strength.iter().sum::<f64>();
        assert!((strength_sum(&coarse) - strength_sum(&graph)).abs() < 1e-9);
        for group in 0..count {
            let members: f64 = (0..n)
                .filter(|&node| fine_labels[node] == group)
                .map(|node| graph.strength[node])
                .sum();
            assert!((coarse.strength[group] - members).abs() < 1e-9);
        }

        // Coarse singletons reproduce the fine partition's modularity ...
        let identity: Vec<usize> = (0..count).collect();
        let fine_q = dense_modularity(&graph, &fine_labels);
        let coarse_q = dense_modularity(&coarse, &identity);
        assert!(
            (fine_q - coarse_q).abs() < 1e-12,
            "seed {seed}: {fine_q} vs {coarse_q}"
        );

        // ... and so does any coarser partition.
        let mut coarse_labels: Vec<usize> = (0..count).map(|_| rng.below(3)).collect();
        densify(&mut coarse_labels);
        let composed: Vec<usize> = fine_labels.iter().map(|&g| coarse_labels[g]).collect();
        let fine_q = dense_modularity(&graph, &composed);
        let coarse_q = dense_modularity(&coarse, &coarse_labels);
        assert!(
            (fine_q - coarse_q).abs() < 1e-12,
            "seed {seed}: {fine_q} vs {coarse_q}"
        );
    }
}

// ---------------------------------------------------------------------------
// The multi-level driver
// ---------------------------------------------------------------------------

#[test]
fn louvain_aggregation_improves_on_a_single_local_moving_phase() {
    // Sparse planted groups: one local-moving phase stops with the groups
    // split into fragments; aggregating the fragments and moving again merges
    // them back.
    for (groups, size, degree, bridges, seed) in
        [(8, 12, 3, 6, 1), (6, 20, 3, 10, 4), (12, 8, 2, 6, 5)]
    {
        let graph = planted_partition(groups, size, degree, bridges, seed);
        let single = local_moving(&graph);
        let full = louvain(&graph);
        let q_single = dense_modularity(&graph, &single.community);
        let q_full = dense_modularity(&graph, &full);
        assert!(
            q_full > q_single + 0.03,
            "groups {groups}, seed {seed}: single level {q_single}, multi-level {q_full}"
        );
        assert!(label_count(&full) < single.count);
    }
}

#[test]
fn louvain_finds_the_cliques_of_a_ring() {
    // Six K5 cliques joined in a ring by single edges.  The resolution limit
    // sqrt(2m) = 11.5 is above the number of cliques, so merging neighbours
    // lowers modularity and the cliques are the optimum.
    let mut edges = Vec::new();
    for c in 0..6 {
        let members: Vec<usize> = (5 * c..5 * c + 5).collect();
        edges.extend(clique_edges(&members, 1.0));
        edges.push((5 * c + 4, (5 * (c + 1)) % 30, 1.0));
    }
    let graph = graph_from_edges(30, &edges);
    let labels = louvain(&graph);
    for c in 0..6 {
        let first = labels[5 * c];
        for &label in &labels[5 * c..5 * c + 5] {
            assert_eq!(label, first, "clique {c}");
        }
    }
    assert_eq!(label_count(&labels), 6);
    // Q = 6 * (10/66 - (22/132)^2).
    let expected = 6.0 * (10.0 / 66.0 - (22.0_f64 / 132.0).powi(2));
    assert!((graph.modularity(&labels) - expected).abs() < 1e-12);
}

#[test]
fn louvain_labels_are_dense_and_numbered_by_first_appearance() {
    let graph = planted_partition(5, 14, 3, 2, 21);
    let labels = louvain(&graph);
    let mut next = 0;
    for &label in &labels {
        assert!(label <= next, "label {label} appears before {next}");
        if label == next {
            next += 1;
        }
    }
    assert_eq!(next, label_count(&labels));
    assert_eq!(labels[0], 0);
}

#[test]
fn louvain_is_deterministic_and_at_least_as_good_as_one_level() {
    let graph = messy_graph(60, 150, 7);
    let first = louvain(&graph);
    assert_eq!(first, louvain(&graph));
    let single = local_moving(&graph);
    assert!(
        dense_modularity(&graph, &first) >= dense_modularity(&graph, &single.community) - 1e-12
    );
}

#[test]
fn louvain_modularity_matches_the_dense_reference() {
    for (name, graph) in [
        ("karate", karate_club()),
        ("connected", random_connected_graph(50, 80, 3)),
        ("messy", messy_graph(40, 100, 5)),
    ] {
        let labels = louvain(&graph);
        let fast = graph.modularity(&labels);
        let dense = dense_modularity(&graph, &labels);
        assert!((fast - dense).abs() < 1e-12, "{name}: {fast} vs {dense}");
    }
}
