// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the community-detection API of
//! `oxiphysics_viz::network_visualization`.
//!
//! They cover three audit findings:
//!
//! * F5083 - `greedy_modularity` re-evaluated the full O(n^2 * deg) modularity
//!   for every candidate move and mixed neighbour counts with weighted totals.
//! * F5084 / F5085 - `spectral_clustering` power-iterated the Laplacian `L`
//!   (converging to the *largest* eigenvector) and then merely rank-sliced it.
//!
//! Every expected value is either derived by hand from the textbook definition
//! `Q = sum_c [ L_c / m - (D_c / 2m)^2 ]` (`L_c` internal edge weight, `D_c`
//! total strength of community `c`, `m` total edge weight) or produced by an
//! independent O(n^2) double-sum implemented in this file.

use std::collections::{BTreeMap, BTreeSet, HashMap};

use oxiphysics_viz::network_visualization::{
    CommunityDetection, NetworkEdge, NetworkGraph, NetworkNode,
};

// ---------------------------------------------------------------------------
// Builders and helpers
// ---------------------------------------------------------------------------

/// Graph with one node per id, inserted in the given order.
fn graph_of(ids: &[usize]) -> NetworkGraph {
    let mut graph = NetworkGraph::new();
    for &id in ids {
        graph.add_node(NetworkNode::new(id, format!("n{id}")));
    }
    graph
}

/// Add an undirected edge of the given weight.
fn link(graph: &mut NetworkGraph, a: usize, b: usize, weight: f64) {
    let mut edge = NetworkEdge::undirected(a, b);
    edge.weight = weight;
    graph.add_edge(edge);
}

/// Connect all members pairwise with undirected edges of the given weight.
fn add_clique(graph: &mut NetworkGraph, members: &[usize], weight: f64) {
    for (i, &a) in members.iter().enumerate() {
        for &b in &members[i + 1..] {
            link(graph, a, b, weight);
        }
    }
}

/// Membership map assigning label `i` to every id of `groups[i]`.
fn membership_of(groups: &[&[usize]]) -> HashMap<usize, usize> {
    let mut membership = HashMap::new();
    for (label, group) in groups.iter().enumerate() {
        for &id in *group {
            membership.insert(id, label);
        }
    }
    membership
}

/// The partition described by a membership map, as a set of sorted groups.
fn partition_of(membership: &HashMap<usize, usize>) -> BTreeSet<Vec<usize>> {
    let mut by_label: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for (&id, &label) in membership {
        by_label.entry(label).or_default().push(id);
    }
    by_label
        .into_values()
        .map(|mut group| {
            group.sort_unstable();
            group
        })
        .collect()
}

/// The partition described by explicit groups, as a set of sorted groups.
fn expected_partition(groups: &[&[usize]]) -> BTreeSet<Vec<usize>> {
    groups
        .iter()
        .map(|group| {
            let mut sorted = group.to_vec();
            sorted.sort_unstable();
            sorted
        })
        .collect()
}

/// Assert that `membership` describes exactly the partition `groups`.  The
/// failure message lists the groups for small graphs and only their sizes for
/// large ones.
fn assert_same_partition(membership: &HashMap<usize, usize>, groups: &[&[usize]], what: &str) {
    let found = partition_of(membership);
    let expected = expected_partition(groups);
    if found == expected {
        return;
    }
    if membership.len() <= 40 {
        panic!("{what}: found {found:?}, expected {expected:?}");
    }
    let found_sizes: Vec<usize> = found.iter().map(Vec::len).collect();
    let expected_sizes: Vec<usize> = expected.iter().map(Vec::len).collect();
    panic!(
        "{what}: partition differs from the planted one (found sizes {found_sizes:?}, \
         expected sizes {expected_sizes:?})"
    );
}

/// Independent reference: Newman's modularity as a dense double sum,
/// `Q = 1/(2m) * sum_ij [A_ij - k_i k_j / (2m)] * delta(c_i, c_j)`.
///
/// `A` is the symmetric weighted adjacency matrix in which every edge counts as
/// undirected, parallel edges add up and a self-loop of weight `w` contributes
/// `A_ii = 2w`.  Every node must appear in `membership`.
fn reference_modularity(graph: &NetworkGraph, membership: &HashMap<usize, usize>) -> f64 {
    let ids: Vec<usize> = graph.nodes.iter().map(|node| node.id).collect();
    let index: HashMap<usize, usize> = ids.iter().enumerate().map(|(i, &id)| (id, i)).collect();
    let n = ids.len();
    let mut a = vec![vec![0.0_f64; n]; n];
    for edge in &graph.edges {
        let (u, v) = (index[&edge.source], index[&edge.target]);
        if u == v {
            a[u][u] += 2.0 * edge.weight;
        } else {
            a[u][v] += edge.weight;
            a[v][u] += edge.weight;
        }
    }
    let strength: Vec<f64> = a.iter().map(|row| row.iter().sum()).collect();
    let two_m: f64 = strength.iter().sum();
    let mut q = 0.0;
    for i in 0..n {
        for j in 0..n {
            if membership[&ids[i]] == membership[&ids[j]] {
                q += a[i][j] - strength[i] * strength[j] / two_m;
            }
        }
    }
    q / two_m
}

/// Deterministic 64-bit LCG (Knuth's MMIX constants) for reproducible graphs.
struct Lcg(u64);

impl Lcg {
    fn next_u64(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.0
    }

    /// Uniform integer in `0..bound`.
    fn below(&mut self, bound: usize) -> usize {
        ((self.next_u64() >> 33) as usize) % bound
    }

    /// Uniform float in `[0, 1)`.
    fn unit(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
}

/// Two triangles with internal weight 2 joined by one bridge of weight 1.
///
/// `m = 13`, `L_c = 6`, `D_c = 13` for each triangle, hence
/// `Q = 2 * (6/13 - (13/26)^2) = 12/13 - 1/2 = 11/26`.
fn weighted_triangles() -> NetworkGraph {
    let mut graph = graph_of(&[0, 1, 2, 3, 4, 5]);
    add_clique(&mut graph, &[0, 1, 2], 2.0);
    add_clique(&mut graph, &[3, 4, 5], 2.0);
    link(&mut graph, 2, 3, 1.0);
    graph
}

/// Two unit-weight triangles joined by one bridge: `Q = 6/7 - 1/2 = 5/14`.
fn unit_triangles() -> NetworkGraph {
    let mut graph = graph_of(&[0, 1, 2, 3, 4, 5]);
    add_clique(&mut graph, &[0, 1, 2], 1.0);
    add_clique(&mut graph, &[3, 4, 5], 1.0);
    link(&mut graph, 2, 3, 1.0);
    graph
}

/// Zachary's karate club (34 members, 78 friendships, 0-based ids).
const KARATE_EDGES: [(usize, usize); 78] = [
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

/// The two factions the club actually split into: Mr. Hi's 16 members and the
/// officer's 18.  Member 9 (id 8) sided with the officer although the structure
/// of the network puts him with Mr. Hi.
const KARATE_MR_HI: [usize; 16] = [0, 1, 2, 3, 4, 5, 6, 7, 10, 11, 12, 13, 16, 17, 19, 21];
const KARATE_OFFICER: [usize; 18] = [
    8, 9, 14, 15, 18, 20, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33,
];

fn karate_club() -> NetworkGraph {
    let ids: Vec<usize> = (0..34).collect();
    let mut graph = graph_of(&ids);
    for &(a, b) in &KARATE_EDGES {
        link(&mut graph, a, b, 1.0);
    }
    graph
}

/// `groups` planted groups of `size` nodes: every node gets `degree` random
/// in-group edges, and `bridges` random edges join every pair of groups.
/// Ids are `group * size + offset`.
fn planted_partition(
    groups: usize,
    size: usize,
    degree: usize,
    bridges: usize,
    seed: u64,
) -> NetworkGraph {
    let ids: Vec<usize> = (0..groups * size).collect();
    let mut graph = graph_of(&ids);
    let mut rng = Lcg(seed);
    for group in 0..groups {
        for offset in 0..size {
            let node = group * size + offset;
            for _ in 0..degree {
                let mut other = rng.below(size);
                if other == offset {
                    other = (other + 1) % size;
                }
                link(&mut graph, node, group * size + other, 1.0);
            }
        }
    }
    for g1 in 0..groups {
        for g2 in g1 + 1..groups {
            for _ in 0..bridges {
                let a = g1 * size + rng.below(size);
                let b = g2 * size + rng.below(size);
                link(&mut graph, a, b, 1.0);
            }
        }
    }
    graph
}

// ---------------------------------------------------------------------------
// modularity()
// ---------------------------------------------------------------------------

#[test]
fn modularity_of_weighted_triangles_matches_the_hand_value() {
    let detector = CommunityDetection::new(weighted_triangles());
    let q = detector.modularity(&membership_of(&[&[0, 1, 2], &[3, 4, 5]]));
    assert!(
        (q - 11.0 / 26.0).abs() < 1e-12,
        "Q = {q}, expected 11/26 = {}",
        11.0 / 26.0
    );
}

#[test]
fn modularity_of_unit_triangles_matches_the_hand_value() {
    let detector = CommunityDetection::new(unit_triangles());
    let q = detector.modularity(&membership_of(&[&[0, 1, 2], &[3, 4, 5]]));
    assert!((q - 5.0 / 14.0).abs() < 1e-12, "Q = {q}, expected 5/14");
}

#[test]
fn modularity_of_trivial_partitions() {
    let detector = CommunityDetection::new(weighted_triangles());
    // Everything in one community: L = m and D = 2m, so Q = 1 - 1 = 0.
    let single = membership_of(&[&[0, 1, 2, 3, 4, 5]]);
    assert!(detector.modularity(&single).abs() < 1e-12);
    // All singletons: Q = -sum (k_i / 2m)^2 with strengths 4, 4, 5, 5, 4, 4.
    let singletons = membership_of(&[&[0], &[1], &[2], &[3], &[4], &[5]]);
    let expected = -(16.0 + 16.0 + 25.0 + 25.0 + 16.0 + 16.0) / (26.0 * 26.0);
    assert!((detector.modularity(&singletons) - expected).abs() < 1e-12);
}

#[test]
fn modularity_of_graphs_without_edges_is_zero() {
    let empty = CommunityDetection::new(NetworkGraph::new());
    assert_eq!(empty.modularity(&HashMap::new()), 0.0);
    let edgeless = CommunityDetection::new(graph_of(&[0, 1, 2]));
    assert_eq!(edgeless.modularity(&membership_of(&[&[0, 1], &[2]])), 0.0);
}

#[test]
fn modularity_treats_nodes_missing_from_the_membership_as_singletons() {
    // Triangle {0,1,2} is one community; the nodes 3, 4, 5 are unassigned and
    // therefore each its own singleton.  Strengths of 3, 4, 5 are 2, 2, 3 (unit
    // weights), m = 7:
    // Q = 3/7 - (7/14)^2 - (3^2 + 2^2 + 2^2)/14^2 = 3/7 - 1/4 - 17/196 = 9/98.
    let detector = CommunityDetection::new(unit_triangles());
    let q = detector.modularity(&membership_of(&[&[0, 1, 2]]));
    let expected = 3.0 / 7.0 - 0.25 - 17.0 / 196.0;
    assert!((q - expected).abs() < 1e-12, "Q = {q}, expected {expected}");
    assert!((expected - 9.0 / 98.0).abs() < 1e-15);
}

#[test]
fn modularity_ignores_edges_that_cannot_take_part() {
    let mut noisy = graph_of(&[0, 1, 2]);
    link(&mut noisy, 0, 1, 1.0);
    link(&mut noisy, 1, 2, 1.0);
    link(&mut noisy, 0, 99, 5.0); // unknown end point
    link(&mut noisy, 0, 2, f64::NAN);
    link(&mut noisy, 0, 2, f64::INFINITY);
    link(&mut noisy, 0, 2, -3.0);
    link(&mut noisy, 0, 2, 0.0);
    let mut clean = graph_of(&[0, 1, 2]);
    link(&mut clean, 0, 1, 1.0);
    link(&mut clean, 1, 2, 1.0);
    let noisy = CommunityDetection::new(noisy);
    let clean = CommunityDetection::new(clean);
    for groups in [
        vec![vec![0, 1, 2]],
        vec![vec![0, 1], vec![2]],
        vec![vec![0], vec![1], vec![2]],
    ] {
        let refs: Vec<&[usize]> = groups.iter().map(Vec::as_slice).collect();
        let membership = membership_of(&refs);
        let q_noisy = noisy.modularity(&membership);
        let q_clean = clean.modularity(&membership);
        assert!(q_noisy.is_finite(), "Q must stay finite for {groups:?}");
        assert!(
            (q_noisy - q_clean).abs() < 1e-12,
            "{groups:?}: {q_noisy} vs {q_clean}"
        );
    }
}

#[test]
fn modularity_counts_directed_edges_as_undirected() {
    // 0 -> 1 (weight 1), 1 -> 2 (weight 2): m = 3, strengths 1, 3, 2.
    // Partition {0,1} | {2}: Q = 1/3 - (4/6)^2 - (2/6)^2 = -2/9.
    let mut directed = graph_of(&[0, 1, 2]);
    directed.add_edge(NetworkEdge::directed(0, 1, 1.0));
    directed.add_edge(NetworkEdge::directed(1, 2, 2.0));
    let mut undirected = graph_of(&[0, 1, 2]);
    link(&mut undirected, 0, 1, 1.0);
    link(&mut undirected, 1, 2, 2.0);
    let membership = membership_of(&[&[0, 1], &[2]]);
    let expected = 1.0 / 3.0 - (4.0_f64 / 6.0).powi(2) - (2.0_f64 / 6.0).powi(2);
    let q_directed = CommunityDetection::new(directed).modularity(&membership);
    let q_undirected = CommunityDetection::new(undirected).modularity(&membership);
    assert!((q_directed - expected).abs() < 1e-12, "{q_directed}");
    assert!((q_undirected - expected).abs() < 1e-12, "{q_undirected}");
}

#[test]
fn modularity_of_a_self_loop_graph_matches_the_hand_value() {
    // Nodes 0 and 1 joined by an edge of weight 1; node 0 also has a self-loop
    // of weight 2.  m = 3, k_0 = 1 + 2*2 = 5, k_1 = 1.
    // {0} | {1}: L_0 = 2, L_1 = 0, D_0 = 5, D_1 = 1:
    // Q = 2/3 - (5/6)^2 - (1/6)^2 = 24/36 - 26/36 = -1/18.
    // (Newman's sum agrees: A_00 = 4, k = (5, 1), 2m = 6, so
    // Q = [(4 - 25/6) + (0 - 1/6)] / 6 = -2/36.)
    let mut graph = graph_of(&[0, 1]);
    link(&mut graph, 0, 1, 1.0);
    link(&mut graph, 0, 0, 2.0);
    let q = CommunityDetection::new(graph).modularity(&membership_of(&[&[0], &[1]]));
    assert!((q - (-1.0 / 18.0)).abs() < 1e-12, "Q = {q}");
}

#[test]
fn modularity_matches_the_dense_double_sum_on_random_graphs() {
    for seed in 1..=12_u64 {
        let mut rng = Lcg(seed.wrapping_mul(0x9E37_79B9_7F4A_7C15));
        let n = 24;
        let ids: Vec<usize> = (0..n).map(|i| 3 * i + 7).collect();
        let mut graph = graph_of(&ids);
        for _ in 0..70 {
            let a = ids[rng.below(n)];
            let b = ids[rng.below(n)];
            let weight = 0.25 + 4.0 * rng.unit();
            if rng.below(3) == 0 {
                graph.add_edge(NetworkEdge::directed(a, b, weight));
            } else {
                link(&mut graph, a, b, weight);
            }
        }
        // Make sure parallel edges and self-loops are present.
        link(&mut graph, ids[0], ids[1], 1.5);
        link(&mut graph, ids[0], ids[1], 0.5);
        link(&mut graph, ids[2], ids[2], 1.25);
        let membership: HashMap<usize, usize> =
            ids.iter().map(|&id| (id, 100 + rng.below(5))).collect();
        let detector = CommunityDetection::new(graph.clone());
        let q = detector.modularity(&membership);
        let expected = reference_modularity(&graph, &membership);
        assert!(
            (q - expected).abs() < 1e-12,
            "seed {seed}: {q} vs reference {expected}"
        );
    }
}

#[test]
fn modularity_of_the_karate_club_factions_matches_the_published_value() {
    let graph = karate_club();
    assert_eq!(graph.edge_count(), 78);
    assert_eq!(graph.degree(33), 17);
    assert_eq!(graph.degree(0), 16);
    let detector = CommunityDetection::new(graph);
    let q = detector.modularity(&membership_of(&[&KARATE_MR_HI, &KARATE_OFFICER]));
    // Q of the real split is quoted as 0.3715 in the literature.
    assert!((q - 0.3715).abs() < 5e-4, "Q = {q}");
}

// ---------------------------------------------------------------------------
// greedy_modularity()
// ---------------------------------------------------------------------------

#[test]
fn greedy_modularity_recovers_two_weighted_triangles() {
    let detector = CommunityDetection::new(weighted_triangles());
    let result = detector.greedy_modularity();
    assert_eq!(
        partition_of(&result.membership),
        expected_partition(&[&[0, 1, 2], &[3, 4, 5]])
    );
    assert!(
        (result.modularity - 11.0 / 26.0).abs() < 1e-12,
        "Q = {}",
        result.modularity
    );
}

#[test]
fn greedy_modularity_recovers_two_k5_joined_by_a_bridge() {
    // m = 21, D_c = 21, L_c = 10 -> Q = 20/21 - 1/2 = 19/42.
    let (left, right) = ([0, 1, 2, 3, 4], [5, 6, 7, 8, 9]);
    let mut graph = graph_of(&[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]);
    add_clique(&mut graph, &left, 1.0);
    add_clique(&mut graph, &right, 1.0);
    link(&mut graph, 4, 5, 1.0);
    let detector = CommunityDetection::new(graph);
    let result = detector.greedy_modularity();
    assert_eq!(
        partition_of(&result.membership),
        expected_partition(&[&left, &right])
    );
    assert!((result.modularity - 19.0 / 42.0).abs() < 1e-12);
    // The reported score must be the modularity of the reported partition.
    let recomputed = detector.modularity(&result.membership);
    assert!((result.modularity - recomputed).abs() < 1e-12);
}

#[test]
fn greedy_modularity_scales_to_hundreds_of_nodes() {
    // 100 disjoint K8 cliques (800 nodes, 2800 edges) with interleaved ids.
    // The previous implementation evaluated an O(n^2) modularity per candidate
    // move and could not finish this in any reasonable time.
    let ids: Vec<usize> = (0..800).collect();
    let mut graph = graph_of(&ids);
    let mut cliques: Vec<Vec<usize>> = Vec::new();
    for c in 0..100 {
        let members: Vec<usize> = (0..8).map(|j| c + 100 * j).collect();
        add_clique(&mut graph, &members, 1.0);
        cliques.push(members);
    }
    let detector = CommunityDetection::new(graph);
    let result = detector.greedy_modularity();
    let refs: Vec<&[usize]> = cliques.iter().map(Vec::as_slice).collect();
    assert_same_partition(
        &result.membership,
        &refs,
        "every clique is its own community",
    );
    assert_eq!(result.num_communities(), 100);
    // Q = 100 * (28/2800 - (56/5600)^2) = 100 * 0.0099 = 0.99.
    assert!((result.modularity - 0.99).abs() < 1e-12);
}

#[test]
fn greedy_modularity_attaches_a_node_to_the_heavier_side() {
    // Two K6 cliques; node 12 hangs off node 0 with weight 3 and off node 6
    // with weight 1.  By neighbour count the two sides tie; by weight it
    // belongs with the first clique.
    let first = [0, 1, 2, 3, 4, 5];
    let second = [6, 7, 8, 9, 10, 11];
    let ids: Vec<usize> = (0..13).collect();
    let mut graph = graph_of(&ids);
    add_clique(&mut graph, &first, 1.0);
    add_clique(&mut graph, &second, 1.0);
    link(&mut graph, 12, 0, 3.0);
    link(&mut graph, 12, 6, 1.0);
    let detector = CommunityDetection::new(graph);
    let result = detector.greedy_modularity();
    let first_plus: Vec<usize> = first.iter().copied().chain([12]).collect();
    assert_same_partition(
        &result.membership,
        &[&first_plus, &second],
        "node 12 joins the heavier side",
    );
    // It really is the best of the alternatives for node 12.
    let q_first = detector.modularity(&membership_of(&[&first_plus, &second]));
    let second_plus: Vec<usize> = second.iter().copied().chain([12]).collect();
    let q_second = detector.modularity(&membership_of(&[&first, &second_plus]));
    let q_alone = detector.modularity(&membership_of(&[&first, &second, &[12]]));
    assert!((result.modularity - q_first).abs() < 1e-12);
    assert!(q_first > q_second && q_first > q_alone);
}

#[test]
fn greedy_modularity_labels_are_dense_and_deterministic() {
    let graph = planted_partition(5, 30, 4, 2, 17);
    let detector = CommunityDetection::new(graph);
    let first = detector.greedy_modularity();
    let second = detector.greedy_modularity();
    assert_eq!(first.membership, second.membership);
    assert_eq!(first.modularity, second.modularity);
    let labels: BTreeSet<usize> = first.membership.values().copied().collect();
    let dense: BTreeSet<usize> = (0..labels.len()).collect();
    assert_eq!(labels, dense, "labels must be exactly 0..num_communities");
    // The first node always lands in community 0.
    assert_eq!(first.membership[&0], 0);
}

#[test]
fn greedy_modularity_handles_degenerate_graphs() {
    let empty = CommunityDetection::new(NetworkGraph::new()).greedy_modularity();
    assert!(empty.membership.is_empty());
    assert_eq!(empty.modularity, 0.0);

    let edgeless = CommunityDetection::new(graph_of(&[4, 8, 15])).greedy_modularity();
    assert_eq!(edgeless.membership.len(), 3);
    assert_eq!(edgeless.num_communities(), 3);
    assert_eq!(edgeless.modularity, 0.0);

    // A node that appears twice in the node list is still one node.
    let mut duplicated = graph_of(&[0, 1, 1, 2]);
    link(&mut duplicated, 0, 1, 1.0);
    link(&mut duplicated, 1, 2, 1.0);
    let result = CommunityDetection::new(duplicated).greedy_modularity();
    assert_eq!(result.membership.len(), 3);
}

#[test]
fn greedy_modularity_survives_directed_parallel_and_self_loop_edges() {
    let mut graph = graph_of(&[10, 11, 12, 13, 14, 15]);
    graph.add_edge(NetworkEdge::directed(10, 11, 2.0));
    graph.add_edge(NetworkEdge::directed(11, 12, 2.0));
    graph.add_edge(NetworkEdge::directed(12, 10, 2.0));
    graph.add_edge(NetworkEdge::directed(13, 14, 2.0));
    graph.add_edge(NetworkEdge::directed(14, 15, 2.0));
    graph.add_edge(NetworkEdge::directed(15, 13, 2.0));
    link(&mut graph, 10, 11, 1.0); // parallel to the directed edge
    link(&mut graph, 12, 12, 3.0); // self-loop
    link(&mut graph, 12, 13, 0.5); // weak bridge
    link(&mut graph, 10, 777, 4.0); // unknown end point
    let detector = CommunityDetection::new(graph.clone());
    let result = detector.greedy_modularity();
    assert_eq!(
        partition_of(&result.membership),
        expected_partition(&[&[10, 11, 12], &[13, 14, 15]])
    );
    let expected = reference_modularity(&graph_without_unknown(&graph), &result.membership);
    assert!((result.modularity - expected).abs() < 1e-12);
}

/// Copy of `graph` without the edges whose end points are not nodes.
fn graph_without_unknown(graph: &NetworkGraph) -> NetworkGraph {
    let known: BTreeSet<usize> = graph.nodes.iter().map(|node| node.id).collect();
    let mut copy = NetworkGraph::new();
    for node in &graph.nodes {
        copy.add_node(node.clone());
    }
    for edge in &graph.edges {
        if known.contains(&edge.source) && known.contains(&edge.target) {
            copy.add_edge(edge.clone());
        }
    }
    copy
}

#[test]
fn greedy_modularity_on_the_karate_club_finds_a_good_partition() {
    let detector = CommunityDetection::new(karate_club());
    let result = detector.greedy_modularity();
    assert_eq!(result.membership.len(), 34);
    // The best partition known for this graph has Q = 0.41979; Louvain
    // typically lands at 0.4156 - 0.4198.
    assert!(
        result.modularity > 0.40 && result.modularity <= 0.41979 + 1e-9,
        "Q = {}",
        result.modularity
    );
    let recomputed = detector.modularity(&result.membership);
    assert!((result.modularity - recomputed).abs() < 1e-12);
    assert!(
        result.modularity > detector.modularity(&membership_of(&[&KARATE_MR_HI, &KARATE_OFFICER])),
        "must beat the two-faction split"
    );
}

#[test]
fn greedy_modularity_on_random_graphs_beats_singletons_and_matches_the_reference() {
    for seed in [3_u64, 9, 27] {
        let graph = planted_partition(4, 40, 3, 3, seed);
        let detector = CommunityDetection::new(graph.clone());
        let result = detector.greedy_modularity();
        let reference = reference_modularity(&graph, &result.membership);
        assert!((result.modularity - reference).abs() < 1e-12);
        // Strong planted structure: the optimiser must find a high-Q partition.
        assert!(
            result.modularity > 0.6,
            "seed {seed}: {}",
            result.modularity
        );
        let singletons: HashMap<usize, usize> =
            graph.nodes.iter().map(|node| (node.id, node.id)).collect();
        assert!(result.modularity > detector.modularity(&singletons));
    }
}

// ---------------------------------------------------------------------------
// spectral_clustering()
// ---------------------------------------------------------------------------

#[test]
fn spectral_k2_recovers_two_cliques_with_interleaved_ids() {
    // Clique A has the even ids and clique B the odd ones, so the cliques are
    // interleaved in node insertion order: an index-based slicing of any
    // ordering that ignores the graph cannot reproduce them.
    let evens = [0, 2, 4, 6, 8];
    let odds = [1, 3, 5, 7, 9];
    let mut graph = graph_of(&[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]);
    add_clique(&mut graph, &evens, 1.0);
    add_clique(&mut graph, &odds, 1.0);
    link(&mut graph, 8, 1, 1.0);
    let detector = CommunityDetection::new(graph);
    let result = detector.spectral_clustering(2);
    assert_eq!(
        partition_of(&result.membership),
        expected_partition(&[&evens, &odds])
    );
    // m = 21, D_c = 21, L_c = 10 -> Q = 19/42.
    assert!(
        (result.modularity - 19.0 / 42.0).abs() < 1e-12,
        "Q = {}",
        result.modularity
    );
}

#[test]
fn spectral_k2_recovers_two_contiguous_cliques() {
    let (left, right) = ([0, 1, 2, 3, 4], [5, 6, 7, 8, 9]);
    let mut graph = graph_of(&[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]);
    add_clique(&mut graph, &left, 1.0);
    add_clique(&mut graph, &right, 1.0);
    link(&mut graph, 4, 5, 1.0);
    let result = CommunityDetection::new(graph).spectral_clustering(2);
    assert_eq!(
        partition_of(&result.membership),
        expected_partition(&[&left, &right])
    );
}

#[test]
fn spectral_k2_follows_the_edge_weights() {
    // Two unit triangles A = {0,1,2}, B = {3,4,5} joined by a unit bridge.
    // Node 6 hangs off A's node 0 with weight 0.1 and off B's node 3 with
    // weight 10, so it belongs with B.  (By neighbour count both sides tie.)
    let mut graph = graph_of(&[0, 1, 2, 3, 4, 5, 6]);
    add_clique(&mut graph, &[0, 1, 2], 1.0);
    add_clique(&mut graph, &[3, 4, 5], 1.0);
    link(&mut graph, 2, 3, 1.0);
    link(&mut graph, 6, 0, 0.1);
    link(&mut graph, 6, 3, 10.0);
    let result = CommunityDetection::new(graph).spectral_clustering(2);
    assert_eq!(
        partition_of(&result.membership),
        expected_partition(&[&[0, 1, 2], &[3, 4, 5, 6]])
    );
}

#[test]
fn spectral_k3_splits_the_components_and_then_the_bridged_pair() {
    // A bridged pair of K5s plus a separate K5: the first bisection must
    // separate the components, the second one the pair along its bridge.
    let (a, b, c) = ([0, 1, 2, 3, 4], [5, 6, 7, 8, 9], [10, 11, 12, 13, 14]);
    let ids: Vec<usize> = (0..15).collect();
    let mut graph = graph_of(&ids);
    add_clique(&mut graph, &a, 1.0);
    add_clique(&mut graph, &b, 1.0);
    add_clique(&mut graph, &c, 1.0);
    link(&mut graph, 4, 5, 1.0);
    let result = CommunityDetection::new(graph).spectral_clustering(3);
    assert_eq!(
        partition_of(&result.membership),
        expected_partition(&[&a, &b, &c])
    );
    // m = 31, D = 21, 21, 20 and L_c = 10 for each clique.
    let expected = 3.0 * 10.0 / 31.0 - (21.0 * 21.0 * 2.0 + 20.0 * 20.0) / (62.0 * 62.0);
    assert!((result.modularity - expected).abs() < 1e-12);
}

#[test]
fn spectral_k3_cuts_a_chain_of_equal_cliques_at_its_bridges() {
    // K5 - K5 - K5 chained by single edges.  The middle clique sits at the
    // zero crossing of the Fiedler vector, so a plain sign split would cut it
    // at noise level; the cut has to follow the bridges.
    let (a, b, c) = ([0, 1, 2, 3, 4], [5, 6, 7, 8, 9], [10, 11, 12, 13, 14]);
    let ids: Vec<usize> = (0..15).collect();
    let mut graph = graph_of(&ids);
    add_clique(&mut graph, &a, 1.0);
    add_clique(&mut graph, &b, 1.0);
    add_clique(&mut graph, &c, 1.0);
    link(&mut graph, 4, 5, 1.0);
    link(&mut graph, 9, 10, 1.0);
    let result = CommunityDetection::new(graph).spectral_clustering(3);
    assert_eq!(
        partition_of(&result.membership),
        expected_partition(&[&a, &b, &c])
    );
}

#[test]
fn spectral_k4_recovers_planted_groups_of_a_larger_graph() {
    // 4 groups x 150 nodes with random in-group edges and 6 random edges
    // between every pair of groups: the recursion on induced subgraphs has to
    // peel the groups apart exactly.
    let graph = planted_partition(4, 150, 6, 6, 5);
    let groups: Vec<Vec<usize>> = (0..4).map(|g| (g * 150..(g + 1) * 150).collect()).collect();
    let refs: Vec<&[usize]> = groups.iter().map(Vec::as_slice).collect();
    let detector = CommunityDetection::new(graph);
    let result = detector.spectral_clustering(4);
    assert_same_partition(&result.membership, &refs, "four planted groups");
    let recomputed = detector.modularity(&result.membership);
    assert!((result.modularity - recomputed).abs() < 1e-12);
}

#[test]
fn spectral_k2_splits_a_large_two_group_graph() {
    let graph = planted_partition(2, 700, 5, 15, 21);
    let first: Vec<usize> = (0..700).collect();
    let second: Vec<usize> = (700..1400).collect();
    let result = CommunityDetection::new(graph).spectral_clustering(2);
    assert_same_partition(&result.membership, &[&first, &second], "two planted groups");
}

#[test]
fn spectral_separates_disconnected_components() {
    // Two disjoint triangles: Q = 2 * (3/6 - (6/12)^2) = 1/2.
    let mut graph = graph_of(&[0, 1, 2, 3, 4, 5]);
    add_clique(&mut graph, &[0, 3, 4], 1.0);
    add_clique(&mut graph, &[1, 2, 5], 1.0);
    let result = CommunityDetection::new(graph).spectral_clustering(2);
    assert_eq!(
        partition_of(&result.membership),
        expected_partition(&[&[0, 3, 4], &[1, 2, 5]])
    );
    assert!((result.modularity - 0.5).abs() < 1e-12);
}

#[test]
fn spectral_symmetrises_directed_edges() {
    let (left, right) = ([0, 1, 2, 3, 4], [5, 6, 7, 8, 9]);
    let mut graph = graph_of(&[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]);
    for group in [left, right] {
        for (i, &a) in group.iter().enumerate() {
            for &b in &group[i + 1..] {
                // One direction only: the Laplacian of the raw digraph would
                // be asymmetric.
                graph.add_edge(NetworkEdge::directed(a, b, 1.0));
            }
        }
    }
    graph.add_edge(NetworkEdge::directed(4, 5, 1.0));
    let result = CommunityDetection::new(graph).spectral_clustering(2);
    assert_eq!(
        partition_of(&result.membership),
        expected_partition(&[&left, &right])
    );
    assert!((result.modularity - 19.0 / 42.0).abs() < 1e-12);
}

#[test]
fn spectral_returns_singletons_with_their_real_modularity_when_k_reaches_n() {
    // K3 with unit weights: three singletons give Q = -3 * (2/6)^2 = -1/3.
    let mut graph = graph_of(&[0, 1, 2]);
    add_clique(&mut graph, &[0, 1, 2], 1.0);
    let detector = CommunityDetection::new(graph);
    for k in [3, 10] {
        let result = detector.spectral_clustering(k);
        assert_eq!(result.num_communities(), 3);
        assert!(
            (result.modularity - (-1.0 / 3.0)).abs() < 1e-12,
            "k = {k}: Q = {}",
            result.modularity
        );
    }
}

#[test]
fn spectral_handles_degenerate_requests() {
    let empty = CommunityDetection::new(NetworkGraph::new()).spectral_clustering(2);
    assert!(empty.membership.is_empty());
    assert_eq!(empty.modularity, 0.0);

    let detector = CommunityDetection::new(unit_triangles());
    let none = detector.spectral_clustering(0);
    assert!(none.membership.is_empty());

    // One community: everything together, Q = 0.
    let one = detector.spectral_clustering(1);
    assert_eq!(one.membership.len(), 6);
    assert_eq!(one.num_communities(), 1);
    assert!(one.modularity.abs() < 1e-12);

    // Without edges there is no structure; the requested number of groups is
    // still produced, deterministically.
    let edgeless = CommunityDetection::new(graph_of(&[0, 1, 2, 3, 4, 5]));
    let first = edgeless.spectral_clustering(3);
    let second = edgeless.spectral_clustering(3);
    assert_eq!(first.num_communities(), 3);
    assert_eq!(first.membership, second.membership);
    assert_eq!(first.modularity, 0.0);
}

#[test]
fn spectral_ignores_edges_that_cannot_take_part() {
    let mut graph = graph_of(&[0, 1, 2, 3, 4, 5]);
    add_clique(&mut graph, &[0, 1, 2], 1.0);
    add_clique(&mut graph, &[3, 4, 5], 1.0);
    link(&mut graph, 2, 3, 1.0);
    link(&mut graph, 0, 99, 5.0); // unknown end point
    link(&mut graph, 0, 4, f64::NAN);
    link(&mut graph, 1, 5, -2.0);
    let result = CommunityDetection::new(graph).spectral_clustering(2);
    assert_eq!(
        partition_of(&result.membership),
        expected_partition(&[&[0, 1, 2], &[3, 4, 5]])
    );
    assert!(result.modularity.is_finite());
}

#[test]
fn spectral_bisection_of_the_karate_club_follows_the_factions() {
    let detector = CommunityDetection::new(karate_club());
    let result = detector.spectral_clustering(2);
    assert_eq!(result.num_communities(), 2);
    let hi_label = result.membership[&0];
    let wrong = KARATE_MR_HI
        .iter()
        .filter(|id| result.membership[id] != hi_label)
        .count()
        + KARATE_OFFICER
            .iter()
            .filter(|id| result.membership[id] == hi_label)
            .count();
    // Spectral bisection of the club is known to reproduce the real split
    // except for a node or two on the boundary.
    assert!(wrong <= 1, "{wrong} members on the wrong side");
    let recomputed = detector.modularity(&result.membership);
    assert!((result.modularity - recomputed).abs() < 1e-12);
}

#[test]
fn spectral_four_way_split_of_the_karate_club_reaches_a_high_modularity() {
    // The best four-way partition of the club has Q = 0.41979.
    let detector = CommunityDetection::new(karate_club());
    let result = detector.spectral_clustering(4);
    assert_eq!(result.num_communities(), 4);
    assert!(result.modularity > 0.41, "Q = {}", result.modularity);
    assert!(
        result.modularity <= 0.41979 + 1e-9,
        "Q = {}",
        result.modularity
    );
}

#[test]
fn both_methods_recover_twenty_planted_communities() {
    // 20 groups of 100 nodes: the spectrum has a cluster of 19 near-zero
    // Laplacian eigenvalues, which a power iteration needs ~10^5 iterations to
    // resolve, and which made the previous spectral implementation take
    // seconds at this size.
    let graph = planted_partition(20, 100, 5, 2, 3);
    let groups: Vec<Vec<usize>> = (0..20)
        .map(|g| (g * 100..(g + 1) * 100).collect())
        .collect();
    let refs: Vec<&[usize]> = groups.iter().map(Vec::as_slice).collect();
    let detector = CommunityDetection::new(graph);
    let louvain = detector.greedy_modularity();
    assert_same_partition(&louvain.membership, &refs, "Louvain");
    let spectral = detector.spectral_clustering(20);
    assert_same_partition(&spectral.membership, &refs, "spectral");
    assert!((louvain.modularity - spectral.modularity).abs() < 1e-12);
}

#[test]
fn spectral_is_deterministic() {
    let graph = planted_partition(3, 60, 4, 3, 99);
    let detector = CommunityDetection::new(graph);
    let first = detector.spectral_clustering(3);
    let second = detector.spectral_clustering(3);
    assert_eq!(first.membership, second.membership);
    assert_eq!(first.modularity, second.modularity);
    // Labels are dense and the first node is in community 0.
    let labels: BTreeSet<usize> = first.membership.values().copied().collect();
    assert_eq!(labels, (0..labels.len()).collect::<BTreeSet<_>>());
    assert_eq!(first.membership[&0], 0);
}
