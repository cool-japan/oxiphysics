// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Unit tests of the spectral kernels: the Jacobi eigensolver, the Lanczos
//! Fiedler vector (checked against an exact dense eigen-decomposition), the
//! conductance sweep and the recursive bisection.

use super::super::WeightedGraph;
use super::super::test_support::{
    clique_edges, graph_from_edges, karate_club, path_graph, planted_partition,
    random_connected_graph,
};
use super::fiedler::{fiedler_vector, fiedler_vector_counting};
use super::jacobi::symmetric_eigen;
use super::{
    Subgraph, best_sweep_prefix, bisect, conductance, largest_splittable, refine_bisection,
    spectral_partition,
};

/// Subgraph induced by every node of `graph`.
fn whole(graph: &WeightedGraph) -> Subgraph {
    let members: Vec<usize> = (0..graph.len()).collect();
    subgraph_of(graph, &members)
}

/// Subgraph induced by `members`.
fn subgraph_of(graph: &WeightedGraph, members: &[usize]) -> Subgraph {
    Subgraph::induced(graph, members, &mut vec![usize::MAX; graph.len()])
}

/// Dense Laplacian `D - A` of a subgraph.
fn dense_laplacian(sub: &Subgraph) -> Vec<Vec<f64>> {
    let n = sub.len();
    let mut l = vec![vec![0.0_f64; n]; n];
    for (i, row) in sub.adj.iter().enumerate() {
        l[i][i] = sub.degree[i];
        for &(j, weight) in row {
            l[i][j] -= weight;
        }
    }
    l
}

/// `v^T L v` for a unit vector `v`.
fn rayleigh(sub: &Subgraph, v: &[f64]) -> f64 {
    let mut lv = vec![0.0_f64; v.len()];
    sub.laplacian_product(v, &mut lv);
    v.iter().zip(&lv).map(|(a, b)| a * b).sum()
}

/// Euclidean norm of `L v - theta v`.
fn residual(sub: &Subgraph, v: &[f64], theta: f64) -> f64 {
    let mut lv = vec![0.0_f64; v.len()];
    sub.laplacian_product(v, &mut lv);
    lv.iter()
        .zip(v)
        .map(|(a, b)| (a - theta * b) * (a - theta * b))
        .sum::<f64>()
        .sqrt()
}

fn norm(v: &[f64]) -> f64 {
    v.iter().map(|x| x * x).sum::<f64>().sqrt()
}

// ---------------------------------------------------------------------------
// Jacobi eigensolver
// ---------------------------------------------------------------------------

#[test]
fn jacobi_diagonalises_a_known_matrix() {
    let a = vec![
        vec![2.0, 1.0, 0.0],
        vec![1.0, 2.0, 1.0],
        vec![0.0, 1.0, 2.0],
    ];
    let (values, vectors) = symmetric_eigen(a.clone());
    let sqrt2 = 2.0_f64.sqrt();
    for (value, expected) in values.iter().zip([2.0 - sqrt2, 2.0, 2.0 + sqrt2]) {
        assert!((value - expected).abs() < 1e-12, "{value} vs {expected}");
    }
    for (value, vector) in values.iter().zip(&vectors) {
        assert!((norm(vector) - 1.0).abs() < 1e-12);
        for (row, expected) in a.iter().zip(vector) {
            let av: f64 = row.iter().zip(vector).map(|(x, y)| x * y).sum();
            assert!((av - value * expected).abs() < 1e-12);
        }
    }
}

#[test]
fn jacobi_reproduces_the_laplacian_spectrum_of_a_path() {
    // The Laplacian of the path on n nodes has eigenvalues 2 - 2 cos(pi k / n).
    let n = 8;
    let sub = whole(&path_graph(n));
    let (values, _) = symmetric_eigen(dense_laplacian(&sub));
    for (k, value) in values.iter().enumerate() {
        let expected = 2.0 - 2.0 * (std::f64::consts::PI * k as f64 / n as f64).cos();
        assert!(
            (value - expected).abs() < 1e-12,
            "k = {k}: {value} vs {expected}"
        );
    }
}

#[test]
fn jacobi_handles_repeated_eigenvalues_and_trivial_sizes() {
    // Laplacian of K4: eigenvalues 0, 4, 4, 4 with orthonormal eigenvectors.
    let graph = graph_from_edges(4, &clique_edges(&[0, 1, 2, 3], 1.0));
    let (values, vectors) = symmetric_eigen(dense_laplacian(&whole(&graph)));
    for (value, expected) in values.iter().zip([0.0, 4.0, 4.0, 4.0]) {
        assert!((value - expected).abs() < 1e-12);
    }
    for (i, a) in vectors.iter().enumerate() {
        for (j, b) in vectors.iter().enumerate() {
            let dot: f64 = a.iter().zip(b).map(|(x, y)| x * y).sum();
            let expected = if i == j { 1.0 } else { 0.0 };
            assert!((dot - expected).abs() < 1e-12, "({i}, {j}): {dot}");
        }
    }
    // Empty and 1 x 1 matrices.
    let (values, vectors) = symmetric_eigen(Vec::new());
    assert!(values.is_empty() && vectors.is_empty());
    let (values, vectors) = symmetric_eigen(vec![vec![-3.5]]);
    assert_eq!(values, vec![-3.5]);
    assert_eq!(vectors, vec![vec![1.0]]);
}

// ---------------------------------------------------------------------------
// Fiedler vector
// ---------------------------------------------------------------------------

/// Check that `v` is a unit, mean-free vector whose Rayleigh quotient equals
/// the second-smallest eigenvalue of the dense Laplacian of `sub`, and that it
/// is an eigenvector up to `residual_tolerance`.
fn assert_is_fiedler_vector(sub: &Subgraph, v: &[f64], what: &str) {
    assert_eq!(v.len(), sub.len(), "{what}");
    assert!((norm(v) - 1.0).abs() < 1e-10, "{what}: norm {}", norm(v));
    let mean: f64 = v.iter().sum::<f64>() / v.len() as f64;
    assert!(mean.abs() < 1e-10, "{what}: mean {mean}");
    let (values, _) = symmetric_eigen(dense_laplacian(sub));
    let lambda_2 = values[1];
    let lambda_max = values[values.len() - 1];
    let theta = rayleigh(sub, v);
    assert!(
        (theta - lambda_2).abs() <= 1e-7 * lambda_max.max(1.0) + 1e-6 * lambda_2.abs(),
        "{what}: Rayleigh quotient {theta} vs lambda_2 {lambda_2}"
    );
    let scale = 2.0 * sub.degree.iter().copied().fold(0.0_f64, f64::max);
    assert!(
        residual(sub, v, theta) <= 1e-6 * scale,
        "{what}: residual {}",
        residual(sub, v, theta)
    );
}

#[test]
fn fiedler_vector_is_the_second_smallest_eigenvector() {
    let two_k5 = {
        let mut edges = clique_edges(&[0, 1, 2, 3, 4], 1.0);
        edges.extend(clique_edges(&[5, 6, 7, 8, 9], 1.0));
        edges.push((4, 5, 1.0));
        graph_from_edges(10, &edges)
    };
    let barbell_weighted = {
        let mut edges = clique_edges(&[0, 1, 2, 3], 2.5);
        edges.extend(clique_edges(&[4, 5, 6, 7, 8, 9], 0.75));
        edges.push((3, 4, 0.3));
        graph_from_edges(10, &edges)
    };
    let mut cases: Vec<(String, WeightedGraph)> = vec![
        ("path 9".into(), path_graph(9)),
        ("path 40".into(), path_graph(40)),
        ("karate club".into(), karate_club()),
        ("two K5".into(), two_k5),
        ("weighted barbell".into(), barbell_weighted),
    ];
    for (n, extra, seed) in [(15, 10, 1), (30, 30, 2), (60, 90, 3), (45, 5, 4)] {
        cases.push((
            format!("random n={n} extra={extra} seed={seed}"),
            random_connected_graph(n, extra, seed),
        ));
    }
    for (what, graph) in &cases {
        let sub = whole(graph);
        let v = fiedler_vector(&sub);
        assert_is_fiedler_vector(&sub, &v, what);
    }
}

#[test]
fn fiedler_vector_is_far_from_the_largest_eigenvector() {
    // A power iteration on L itself would return the top eigenvector, whose
    // Rayleigh quotient is lambda_max.
    for graph in [karate_club(), random_connected_graph(40, 40, 7)] {
        let sub = whole(&graph);
        let (values, _) = symmetric_eigen(dense_laplacian(&sub));
        let lambda_max = values[values.len() - 1];
        let theta = rayleigh(&sub, &fiedler_vector(&sub));
        assert!(
            theta < 0.2 * lambda_max,
            "{theta} vs lambda_max {lambda_max}"
        );
    }
}

#[test]
fn fiedler_vector_resolves_a_cluster_of_small_eigenvalues() {
    // Six near-isolated communities give five Laplacian eigenvalues close to
    // zero; a power iteration needs ~10^5 products here.  Lanczos has to find
    // the true lambda_2 within the product budget.
    let graph = planted_partition(6, 20, 5, 1, 11);
    let sub = whole(&graph);
    let mut products = 0;
    let v = fiedler_vector_counting(&sub, &mut products);
    assert_is_fiedler_vector(&sub, &v, "planted partition");
    assert!(
        products < 1_000,
        "needed {products} Laplacian products for 120 nodes"
    );
}

#[test]
fn fiedler_vector_of_a_large_clustered_graph_converges_within_the_budget() {
    let graph = planted_partition(20, 250, 5, 2, 3);
    let sub = whole(&graph);
    let mut products = 0;
    let v = fiedler_vector_counting(&sub, &mut products);
    assert!(products <= 2_000, "{products} products");
    // The vector must lie in the low end of the spectrum: its Rayleigh
    // quotient is far below the average degree.
    let average_degree = sub.degree.iter().sum::<f64>() / sub.len() as f64;
    assert!(rayleigh(&sub, &v) < 0.05 * average_degree);
    assert!((norm(&v) - 1.0).abs() < 1e-10);
}

#[test]
fn fiedler_vector_of_a_single_edge_is_antisymmetric() {
    let sub = whole(&graph_from_edges(2, &[(0, 1, 3.0)]));
    let mut products = 0;
    let v = fiedler_vector_counting(&sub, &mut products);
    let h = std::f64::consts::FRAC_1_SQRT_2;
    assert!((v[0].abs() - h).abs() < 1e-12 && (v[1].abs() - h).abs() < 1e-12);
    assert!(v[0] * v[1] < 0.0);
    assert_eq!(products, 1);
}

#[test]
fn fiedler_vector_without_edges_is_the_normalised_start_vector() {
    let sub = whole(&graph_from_edges(6, &[]));
    let mut products = 0;
    let v = fiedler_vector_counting(&sub, &mut products);
    assert_eq!(products, 0);
    assert!((norm(&v) - 1.0).abs() < 1e-12);
    assert!(v.iter().sum::<f64>().abs() < 1e-12);
    // Tiny graphs.
    assert!(fiedler_vector(&whole(&graph_from_edges(0, &[]))).is_empty());
    assert_eq!(fiedler_vector(&whole(&graph_from_edges(1, &[]))), vec![0.0]);
}

#[test]
fn fiedler_vector_of_a_disconnected_graph_is_constant_on_components() {
    // lambda_2 = 0 has the component indicators (minus their mean) as
    // eigenvectors, so the Fiedler vector is constant on every component.
    let mut edges = clique_edges(&[0, 1, 2], 1.0);
    edges.extend(clique_edges(&[3, 4, 5, 6], 1.0));
    edges.extend(clique_edges(&[7, 8], 2.0));
    let sub = whole(&graph_from_edges(9, &edges));
    let v = fiedler_vector(&sub);
    for component in [&[0, 1, 2][..], &[3, 4, 5, 6], &[7, 8]] {
        let first = v[component[0]];
        for &node in component {
            assert!(
                (v[node] - first).abs() < 1e-6,
                "node {node}: {} vs {first}",
                v[node]
            );
        }
    }
    assert!(rayleigh(&sub, &v).abs() < 1e-8);
    assert!((norm(&v) - 1.0).abs() < 1e-10);
}

#[test]
fn fiedler_vector_is_deterministic() {
    let sub = whole(&random_connected_graph(50, 60, 9));
    assert_eq!(fiedler_vector(&sub), fiedler_vector(&sub));
}

// ---------------------------------------------------------------------------
// Sweep cut
// ---------------------------------------------------------------------------

#[test]
fn sweep_cuts_a_chain_of_cliques_at_a_bridge() {
    // K5 - K5 - K5 chained by single edges, nodes ordered clique by clique.
    // Every proper prefix that ends at a bridge has conductance 1/21; the
    // first of them wins.
    let mut edges = clique_edges(&[0, 1, 2, 3, 4], 1.0);
    edges.extend(clique_edges(&[5, 6, 7, 8, 9], 1.0));
    edges.extend(clique_edges(&[10, 11, 12, 13, 14], 1.0));
    edges.push((4, 5, 1.0));
    edges.push((9, 10, 1.0));
    let sub = whole(&graph_from_edges(15, &edges));
    let order: Vec<usize> = (0..15).collect();
    assert_eq!(best_sweep_prefix(&sub, &order), Some(5));
    // The same cut is found whichever way round the order runs.
    let reversed: Vec<usize> = (0..15).rev().collect();
    assert_eq!(best_sweep_prefix(&sub, &reversed), Some(5));
}

#[test]
fn sweep_ignores_the_order_inside_clusters() {
    // Shuffle the nodes inside each K5 of a bridged pair: the best prefix is
    // still exactly one clique.
    let mut edges = clique_edges(&[0, 1, 2, 3, 4], 1.0);
    edges.extend(clique_edges(&[5, 6, 7, 8, 9], 1.0));
    edges.push((4, 5, 1.0));
    let sub = whole(&graph_from_edges(10, &edges));
    let order = vec![3, 0, 4, 2, 1, 7, 9, 5, 8, 6];
    assert_eq!(best_sweep_prefix(&sub, &order), Some(5));
}

#[test]
fn sweep_prefers_the_balanced_cut_among_equal_conductances() {
    // Four disjoint triangles: every prefix that is a union of triangles has
    // conductance 0.  The balanced one (two triangles) has to win.
    let mut edges = Vec::new();
    for t in 0..4 {
        edges.extend(clique_edges(&[3 * t, 3 * t + 1, 3 * t + 2], 1.0));
    }
    let sub = whole(&graph_from_edges(12, &edges));
    let order: Vec<usize> = (0..12).collect();
    assert_eq!(best_sweep_prefix(&sub, &order), Some(6));
}

#[test]
fn sweep_of_a_path_cuts_in_the_middle() {
    for n in [4, 7, 10, 25] {
        let sub = whole(&path_graph(n));
        let order: Vec<usize> = (0..n).collect();
        let cut = best_sweep_prefix(&sub, &order).expect("a path has edges");
        // Even paths cut exactly in the middle; odd ones next to it.
        assert!(cut == n / 2 || cut == n.div_ceil(2), "n = {n}: cut {cut}");
    }
}

#[test]
fn sweep_without_edges_has_no_cut() {
    let sub = whole(&graph_from_edges(5, &[]));
    let order: Vec<usize> = (0..5).collect();
    assert_eq!(best_sweep_prefix(&sub, &order), None);
    // Fewer than two nodes cannot be cut either.
    let single = whole(&graph_from_edges(1, &[]));
    assert_eq!(best_sweep_prefix(&single, &[0]), None);
}

// ---------------------------------------------------------------------------
// Bisection and recursion
// ---------------------------------------------------------------------------

#[test]
fn bisect_splits_two_cliques_apart_and_keeps_members_ascending() {
    // Interleaved cliques: evens and odds.
    let evens = [0, 2, 4, 6, 8];
    let odds = [1, 3, 5, 7, 9];
    let mut edges = clique_edges(&evens, 1.0);
    edges.extend(clique_edges(&odds, 1.0));
    edges.push((8, 1, 1.0));
    let graph = graph_from_edges(10, &edges);
    let members: Vec<usize> = (0..10).collect();
    let (lower, upper) = bisect(&graph, &members, &mut [usize::MAX; 10]);
    let mut sides = vec![lower, upper];
    sides.sort();
    assert_eq!(sides, vec![evens.to_vec(), odds.to_vec()]);
}

#[test]
fn bisect_works_on_an_induced_subgraph() {
    // Nodes 0..4 are a K5, 5..9 another K5, 10..14 a third; only 5..14 are
    // bisected, so the edges to node 0..4 do not exist for the subgraph.
    let mut edges = clique_edges(&[0, 1, 2, 3, 4], 1.0);
    edges.extend(clique_edges(&[5, 6, 7, 8, 9], 1.0));
    edges.extend(clique_edges(&[10, 11, 12, 13, 14], 1.0));
    edges.push((4, 5, 1.0));
    edges.push((9, 10, 1.0));
    let graph = graph_from_edges(15, &edges);
    let members: Vec<usize> = (5..15).collect();
    let mut local_of = vec![usize::MAX; 15];
    let (lower, upper) = bisect(&graph, &members, &mut local_of);
    assert!(
        local_of.iter().all(|&slot| slot == usize::MAX),
        "scratch restored"
    );
    let mut sides = vec![lower, upper];
    sides.sort();
    assert_eq!(sides, vec![vec![5, 6, 7, 8, 9], vec![10, 11, 12, 13, 14]]);
}

#[test]
fn bisect_halves_members_without_edges_by_index() {
    let graph = graph_from_edges(7, &[]);
    let members = vec![1, 2, 4, 5, 6];
    let (lower, upper) = bisect(&graph, &members, &mut [usize::MAX; 7]);
    assert_eq!(lower, vec![1, 2]);
    assert_eq!(upper, vec![4, 5, 6]);
    // Two nodes: one on each side.
    let (lower, upper) = bisect(&graph, &[3, 6], &mut [usize::MAX; 7]);
    assert_eq!((lower, upper), (vec![3], vec![6]));
}

#[test]
fn largest_splittable_picks_the_largest_part_and_the_earlier_on_ties() {
    let parts = vec![vec![0], vec![1, 2], vec![3, 4, 5], vec![6, 7, 8], vec![9]];
    assert_eq!(largest_splittable(&parts), Some(2));
    assert_eq!(largest_splittable(&[vec![0], vec![1]]), None);
    assert_eq!(largest_splittable(&[]), None);
}

#[test]
fn spectral_partition_yields_exactly_k_dense_labels() {
    let graph = random_connected_graph(24, 20, 5);
    for k in 1..=24 {
        let labels = spectral_partition(&graph, k);
        assert_eq!(labels.len(), 24);
        let mut seen = [false; 24];
        for &label in &labels {
            seen[label] = true;
        }
        let distinct = seen.iter().filter(|&&s| s).count();
        assert_eq!(distinct, k, "k = {k}");
        // Dense and numbered by first appearance.
        assert!(seen[..k].iter().all(|&s| s), "labels are 0..k for k = {k}");
        let mut next = 0;
        for &label in &labels {
            assert!(label <= next, "k = {k}: label {label} before {next}");
            if label == next {
                next += 1;
            }
        }
    }
    // k beyond n: singletons.
    assert_eq!(spectral_partition(&graph, 100), (0..24).collect::<Vec<_>>());
}

#[test]
fn spectral_partition_is_deterministic() {
    let graph = planted_partition(4, 40, 4, 3, 17);
    assert_eq!(spectral_partition(&graph, 4), spectral_partition(&graph, 4));
}

#[test]
fn fiedler_vector_survives_krylov_spaces_smaller_than_the_graph() {
    // Structurally identical nodes ("twins") give Laplacian eigenvalues of
    // multiplicity > 1, so the Krylov space of any start vector is smaller than
    // the space of mean-free vectors and the iteration runs into an invariant
    // subspace before `size - 1` steps.  Continuing past that point used to
    // feed rounding noise into the basis, whose constant-vector component
    // produced a bogus Ritz value of 0 and a non-mean-free "Fiedler" vector.
    let star = graph_from_edges(12, &(1..12).map(|leaf| (0, leaf, 1.0)).collect::<Vec<_>>());
    let bipartite = {
        let mut edges = Vec::new();
        for left in 0..4 {
            for right in 4..9 {
                edges.push((left, right, 1.0));
            }
        }
        graph_from_edges(9, &edges)
    };
    // A "caterpillar": a path whose nodes each carry three identical leaves.
    let caterpillar = {
        let spine = 6;
        let mut edges: Vec<(usize, usize, f64)> = (0..spine - 1).map(|i| (i, i + 1, 1.0)).collect();
        for i in 0..spine {
            for leaf in 0..3 {
                edges.push((i, spine + 3 * i + leaf, 1.0));
            }
        }
        graph_from_edges(spine + 3 * spine, &edges)
    };
    // The karate club has five nodes with the same two neighbours.
    for (what, graph) in [
        ("star", star),
        ("complete bipartite", bipartite),
        ("caterpillar", caterpillar),
        ("karate club", karate_club()),
    ] {
        let sub = whole(&graph);
        let v = fiedler_vector(&sub);
        assert_is_fiedler_vector(&sub, &v, what);
    }
}

#[test]
fn fiedler_vector_gives_twin_nodes_equal_entries() {
    // The five karate club members 14, 15, 18, 20 and 22 are connected to the
    // same two members only; the Fiedler vector (lambda_2 is simple) cannot
    // tell them apart.
    let sub = whole(&karate_club());
    let v = fiedler_vector(&sub);
    for twin in [15, 18, 20, 22] {
        assert!(
            (v[twin] - v[14]).abs() < 1e-6,
            "twin {twin}: {} vs {}",
            v[twin],
            v[14]
        );
    }
}

// ---------------------------------------------------------------------------
// Refinement
// ---------------------------------------------------------------------------

/// Conductance of the bisection `in_lower`, computed from scratch.
fn bisection_conductance(sub: &Subgraph, in_lower: &[bool]) -> f64 {
    let mut cut = 0.0;
    let mut volume_lower = 0.0;
    let mut volume_upper = 0.0;
    for (node, row) in sub.adj.iter().enumerate() {
        if in_lower[node] {
            volume_lower += sub.degree[node];
        } else {
            volume_upper += sub.degree[node];
        }
        for &(other, weight) in row {
            if in_lower[node] != in_lower[other] {
                cut += 0.5 * weight;
            }
        }
    }
    conductance(cut, volume_lower, volume_upper)
}

#[test]
fn refinement_moves_misplaced_nodes_to_their_own_side() {
    // Two K6 cliques joined by one edge; node 5 (clique A) starts on the side
    // of clique B and node 6 (clique B) on the side of clique A.
    let mut edges = clique_edges(&[0, 1, 2, 3, 4, 5], 1.0);
    edges.extend(clique_edges(&[6, 7, 8, 9, 10, 11], 1.0));
    edges.push((5, 6, 1.0));
    let sub = whole(&graph_from_edges(12, &edges));
    let mut in_lower: Vec<bool> = (0..12).map(|node| node < 6).collect();
    in_lower[5] = false;
    in_lower[6] = true;
    refine_bisection(&sub, &mut in_lower);
    let expected: Vec<bool> = (0..12).map(|node| node < 6).collect();
    assert_eq!(in_lower, expected);
}

#[test]
fn refinement_never_raises_the_conductance_and_ends_in_a_local_optimum() {
    for seed in 1..=10 {
        let sub = whole(&random_connected_graph(30, 45, seed));
        let mut rng_state = seed * 31 + 5;
        let mut in_lower: Vec<bool> = (0..30)
            .map(|_| {
                rng_state = rng_state
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1_442_695_040_888_963_407);
                (rng_state >> 33) % 2 == 0
            })
            .collect();
        // Both sides need a node to start with.
        in_lower[0] = true;
        in_lower[1] = false;
        let before = bisection_conductance(&sub, &in_lower);
        refine_bisection(&sub, &mut in_lower);
        let after = bisection_conductance(&sub, &in_lower);
        assert!(after <= before + 1e-12, "seed {seed}: {before} -> {after}");
        let lower_count = in_lower.iter().filter(|&&b| b).count();
        assert!(
            (1..30).contains(&lower_count),
            "seed {seed}: sides stay non-empty"
        );
        // Local optimum: no single move lowers the conductance any further.
        for node in 0..30 {
            let own_count = if in_lower[node] {
                lower_count
            } else {
                30 - lower_count
            };
            if own_count <= 1 {
                continue;
            }
            let mut moved = in_lower.clone();
            moved[node] = !moved[node];
            assert!(
                bisection_conductance(&sub, &moved) >= after - 1e-9,
                "seed {seed}: moving node {node} still improves {after}"
            );
        }
    }
}

#[test]
fn refinement_keeps_both_sides_non_empty_and_is_idempotent() {
    // A single edge: either node would like to join the other, but then one
    // side would be empty.
    let edge = whole(&graph_from_edges(2, &[(0, 1, 1.0)]));
    let mut sides = vec![true, false];
    refine_bisection(&edge, &mut sides);
    assert_eq!(sides, vec![true, false]);

    let sub = whole(&karate_club());
    let mut in_lower: Vec<bool> = (0..34).map(|node| node % 3 == 0).collect();
    refine_bisection(&sub, &mut in_lower);
    let once = in_lower.clone();
    refine_bisection(&sub, &mut in_lower);
    assert_eq!(once, in_lower);
}

#[test]
fn refinement_leaves_perfect_cuts_alone() {
    // Two components: a cut of weight zero cannot be improved.
    let mut edges = clique_edges(&[0, 1, 2, 3], 1.0);
    edges.extend(clique_edges(&[4, 5, 6], 2.0));
    let sub = whole(&graph_from_edges(7, &edges));
    let mut in_lower: Vec<bool> = (0..7).map(|node| node < 4).collect();
    refine_bisection(&sub, &mut in_lower);
    assert_eq!(in_lower, vec![true, true, true, true, false, false, false]);
}

#[test]
fn conductance_helper_handles_empty_sides() {
    assert_eq!(conductance(3.0, 6.0, 4.0), 0.75);
    assert_eq!(conductance(0.0, 5.0, 5.0), 0.0);
    assert!(conductance(1.0, 0.0, 5.0).is_infinite());
    // A rounding-level negative cut counts as zero.
    assert_eq!(conductance(-1e-18, 2.0, 2.0), 0.0);
}
