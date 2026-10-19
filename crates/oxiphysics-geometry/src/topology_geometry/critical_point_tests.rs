// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tests of the Reeb graph (contour tree) and the Morse complex.

use super::{CriticalPointKind, MorseComplex, ReebGraph};
use std::collections::{HashMap, HashSet};

// ---------------------------------------------------------------------------
// Mesh builders
// ---------------------------------------------------------------------------

/// Deterministic pseudo-random numbers in `[0, 1)` (64-bit LCG).
fn lcg(seed: u64) -> impl FnMut() -> f64 {
    let mut state = seed;
    move || {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (state >> 11) as f64 / (1u64 << 53) as f64
    }
}

/// Edges of the triangulated `nx` x `ny` grid with vertex `(i, j)` at index
/// `j * nx + i`; every cell is split by the diagonal `(i, j)-(i+1, j+1)`.
fn grid_edges(nx: usize, ny: usize) -> Vec<[usize; 2]> {
    let mut edges = Vec::new();
    for j in 0..ny {
        for i in 0..nx {
            let v = j * nx + i;
            if i + 1 < nx {
                edges.push([v, v + 1]);
            }
            if j + 1 < ny {
                edges.push([v, v + nx]);
            }
            if i + 1 < nx && j + 1 < ny {
                edges.push([v, v + nx + 1]);
            }
        }
    }
    edges
}

/// Edges of the triangulated torus: the `nx` x `ny` grid with wrap-around.
fn torus_edges(nx: usize, ny: usize) -> Vec<[usize; 2]> {
    let mut edges = Vec::new();
    for j in 0..ny {
        for i in 0..nx {
            let v = j * nx + i;
            edges.push([v, j * nx + (i + 1) % nx]);
            edges.push([v, ((j + 1) % ny) * nx + i]);
            edges.push([v, ((j + 1) % ny) * nx + (i + 1) % nx]);
        }
    }
    edges
}

/// Vertices and edges of an icosphere with the given number of subdivisions.
fn icosphere(subdivisions: usize) -> (Vec<[f64; 3]>, Vec<[usize; 2]>) {
    let t = (1.0 + 5.0_f64.sqrt()) / 2.0;
    let mut vertices: Vec<[f64; 3]> = vec![
        [-1.0, t, 0.0],
        [1.0, t, 0.0],
        [-1.0, -t, 0.0],
        [1.0, -t, 0.0],
        [0.0, -1.0, t],
        [0.0, 1.0, t],
        [0.0, -1.0, -t],
        [0.0, 1.0, -t],
        [t, 0.0, -1.0],
        [t, 0.0, 1.0],
        [-t, 0.0, -1.0],
        [-t, 0.0, 1.0],
    ];
    let normalise = |p: [f64; 3]| {
        let n = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
        [p[0] / n, p[1] / n, p[2] / n]
    };
    vertices = vertices.into_iter().map(normalise).collect();
    let mut faces: Vec<[usize; 3]> = vec![
        [0, 11, 5],
        [0, 5, 1],
        [0, 1, 7],
        [0, 7, 10],
        [0, 10, 11],
        [1, 5, 9],
        [5, 11, 4],
        [11, 10, 2],
        [10, 7, 6],
        [7, 1, 8],
        [3, 9, 4],
        [3, 4, 2],
        [3, 2, 6],
        [3, 6, 8],
        [3, 8, 9],
        [4, 9, 5],
        [2, 4, 11],
        [6, 2, 10],
        [8, 6, 7],
        [9, 8, 1],
    ];
    for _ in 0..subdivisions {
        let mut midpoints: HashMap<(usize, usize), usize> = HashMap::new();
        let mut next_faces = Vec::with_capacity(faces.len() * 4);
        for f in &faces {
            let mut mid = |a: usize, b: usize| -> usize {
                let key = (a.min(b), a.max(b));
                *midpoints.entry(key).or_insert_with(|| {
                    let p = [
                        (vertices[a][0] + vertices[b][0]) / 2.0,
                        (vertices[a][1] + vertices[b][1]) / 2.0,
                        (vertices[a][2] + vertices[b][2]) / 2.0,
                    ];
                    vertices.push(normalise(p));
                    vertices.len() - 1
                })
            };
            let ab = mid(f[0], f[1]);
            let bc = mid(f[1], f[2]);
            let ca = mid(f[2], f[0]);
            next_faces.push([f[0], ab, ca]);
            next_faces.push([f[1], bc, ab]);
            next_faces.push([f[2], ca, bc]);
            next_faces.push([ab, bc, ca]);
        }
        faces = next_faces;
    }
    let mut edge_set: HashSet<(usize, usize)> = HashSet::new();
    for f in &faces {
        for k in 0..3 {
            let (a, b) = (f[k], f[(k + 1) % 3]);
            edge_set.insert((a.min(b), a.max(b)));
        }
    }
    let mut edges: Vec<[usize; 2]> = edge_set.into_iter().map(|(a, b)| [a, b]).collect();
    edges.sort_unstable();
    (vertices, edges)
}

/// Height along a tilted direction on the torus of revolution with major radius
/// 2 and minor radius 1, sampled on the `nx` x `ny` grid of [`torus_edges`].
/// Such a height function has exactly four critical points: one minimum, two
/// saddles and one maximum.
fn tilted_torus_height(nx: usize, ny: usize) -> Vec<f64> {
    let tilt = [0.2_f64, 0.1, 1.0];
    let norm = (tilt[0] * tilt[0] + tilt[1] * tilt[1] + tilt[2] * tilt[2]).sqrt();
    let mut values = Vec::with_capacity(nx * ny);
    for j in 0..ny {
        for i in 0..nx {
            let u = 2.0 * std::f64::consts::PI * i as f64 / nx as f64;
            let v = 2.0 * std::f64::consts::PI * j as f64 / ny as f64;
            let ring = 2.0 + v.cos();
            let p = [ring * u.cos(), ring * u.sin(), v.sin()];
            values.push((p[0] * tilt[0] + p[1] * tilt[1] + p[2] * tilt[2]) / norm);
        }
    }
    values
}

/// Euler characteristic `V - E + F` of the clique complex of `edges` (every
/// triangle of edges counts as a face).
fn clique_euler_characteristic(n: usize, edges: &[[usize; 2]]) -> i64 {
    let mut adjacency: Vec<HashSet<usize>> = vec![HashSet::new(); n];
    for &[a, b] in edges {
        adjacency[a].insert(b);
        adjacency[b].insert(a);
    }
    let edge_count: usize = adjacency.iter().map(|s| s.len()).sum::<usize>() / 2;
    let mut triangles = 0usize;
    for a in 0..n {
        for &b in &adjacency[a] {
            if b <= a {
                continue;
            }
            for &c in &adjacency[b] {
                if c > b && adjacency[a].contains(&c) {
                    triangles += 1;
                }
            }
        }
    }
    n as i64 - edge_count as i64 + triangles as i64
}

// ---------------------------------------------------------------------------
// Independent (brute force) reference computations
// ---------------------------------------------------------------------------

/// Rank of every vertex in the order (value, index).
fn ranks(values: &[f64]) -> Vec<usize> {
    let mut order: Vec<usize> = (0..values.len()).collect();
    order.sort_by(|&a, &b| {
        (values[a] + 0.0)
            .total_cmp(&(values[b] + 0.0))
            .then(a.cmp(&b))
    });
    let mut rank = vec![0; values.len()];
    for (pos, &v) in order.iter().enumerate() {
        rank[v] = pos;
    }
    rank
}

fn adjacency(n: usize, edges: &[[usize; 2]]) -> Vec<Vec<usize>> {
    let mut adj: Vec<HashSet<usize>> = vec![HashSet::new(); n];
    for &[a, b] in edges {
        if a < n && b < n && a != b {
            adj[a].insert(b);
            adj[b].insert(a);
        }
    }
    adj.into_iter()
        .map(|s| {
            let mut v: Vec<usize> = s.into_iter().collect();
            v.sort_unstable();
            v
        })
        .collect()
}

/// Label of the connected component of every vertex in the subgraph induced by
/// the vertices whose rank lies in `[lo, hi]`; `usize::MAX` outside the slab.
fn slab_components(adj: &[Vec<usize>], rank: &[usize], lo: usize, hi: usize) -> Vec<usize> {
    let n = adj.len();
    let mut label = vec![usize::MAX; n];
    for start in 0..n {
        if label[start] != usize::MAX || rank[start] < lo || rank[start] > hi {
            continue;
        }
        let mut stack = vec![start];
        label[start] = start;
        while let Some(x) = stack.pop() {
            for &y in &adj[x] {
                if rank[y] >= lo && rank[y] <= hi && label[y] == usize::MAX {
                    label[y] = start;
                    stack.push(y);
                }
            }
        }
    }
    label
}

/// Number of distinct sublevel (`below`) or superlevel components among the
/// neighbours of `v` on that side.
fn side_component_count(adj: &[Vec<usize>], rank: &[usize], v: usize, below: bool) -> usize {
    let n = adj.len();
    let (lo, hi) = if below {
        (0, rank[v].saturating_sub(1))
    } else {
        (rank[v] + 1, n.saturating_sub(1))
    };
    let neighbours: Vec<usize> = adj[v]
        .iter()
        .copied()
        .filter(|&u| {
            if below {
                rank[u] < rank[v]
            } else {
                rank[u] > rank[v]
            }
        })
        .collect();
    if neighbours.is_empty() {
        return 0;
    }
    let label = slab_components(adj, rank, lo, hi);
    let distinct: HashSet<usize> = neighbours.iter().map(|&u| label[u]).collect();
    distinct.len()
}

/// Check that the arcs record the connectivity of the sublevel and superlevel
/// sets: two critical points lie in the same component of `{rank <= r}` exactly
/// when the path between them in the tree stays at or below `r`, and in the same
/// component of `{rank >= r}` exactly when it stays at or above `r`.  The
/// thresholds are the ranks of the critical points, because the components only
/// change there.
fn assert_arcs_preserve_level_set_connectivity(adj: &[Vec<usize>], rank: &[usize], rg: &ReebGraph) {
    let last = adj.len().saturating_sub(1);
    let count = rg.nodes.len();
    let mut tree: Vec<Vec<usize>> = vec![Vec::new(); count];
    for arc in &rg.arcs {
        tree[arc.source].push(arc.target);
        tree[arc.target].push(arc.source);
    }
    let node_rank: Vec<usize> = rg.nodes.iter().map(|node| rank[node.vertex_idx]).collect();
    let sublevel: Vec<Vec<usize>> = node_rank
        .iter()
        .map(|&r| slab_components(adj, rank, 0, r))
        .collect();
    let superlevel: Vec<Vec<usize>> = node_rank
        .iter()
        .map(|&r| slab_components(adj, rank, r, last))
        .collect();
    for a in 0..count {
        // Highest and lowest node rank on the tree path from `a` to every node;
        // `None` for the nodes of another tree.
        let mut path: Vec<Option<(usize, usize)>> = vec![None; count];
        path[a] = Some((node_rank[a], node_rank[a]));
        let mut stack = vec![a];
        while let Some(x) = stack.pop() {
            let Some((high, low)) = path[x] else { continue };
            for &y in &tree[x] {
                if path[y].is_none() {
                    path[y] = Some((high.max(node_rank[y]), low.min(node_rank[y])));
                    stack.push(y);
                }
            }
        }
        for b in 0..count {
            let (va, vb) = (rg.nodes[a].vertex_idx, rg.nodes[b].vertex_idx);
            for t in 0..count {
                let r = node_rank[t];
                if r >= node_rank[a].max(node_rank[b]) {
                    let same = sublevel[t][va] == sublevel[t][vb];
                    let expected = path[b].is_some_and(|(high, _)| high <= r);
                    assert_eq!(
                        same, expected,
                        "sublevel connectivity of nodes {a} and {b} at rank {r}"
                    );
                }
                if r <= node_rank[a].min(node_rank[b]) {
                    let same = superlevel[t][va] == superlevel[t][vb];
                    let expected = path[b].is_some_and(|(_, low)| low >= r);
                    assert_eq!(
                        same, expected,
                        "superlevel connectivity of nodes {a} and {b} at rank {r}"
                    );
                }
            }
        }
    }
}

/// Check a Reeb graph against brute-force sublevel / superlevel components.
fn assert_contour_tree(values: &[f64], edges: &[[usize; 2]], rg: &ReebGraph) {
    let n = values.len();
    let adj = adjacency(n, edges);
    let rank = ranks(values);

    // Nodes are ascending in the vertex order and unique.
    for pair in rg.nodes.windows(2) {
        assert!(rank[pair[0].vertex_idx] < rank[pair[1].vertex_idx]);
    }
    let node_of: HashMap<usize, usize> = rg
        .nodes
        .iter()
        .enumerate()
        .map(|(i, node)| (node.vertex_idx, i))
        .collect();

    let mut below_arcs = vec![0usize; rg.nodes.len()];
    let mut above_arcs = vec![0usize; rg.nodes.len()];
    let mut seen_arcs: HashSet<(usize, usize)> = HashSet::new();
    // Union-find over the nodes: the arcs must not close a cycle.
    let mut tree_of: Vec<usize> = (0..rg.nodes.len()).collect();
    fn find(tree_of: &mut [usize], mut x: usize) -> usize {
        while tree_of[x] != x {
            tree_of[x] = tree_of[tree_of[x]];
            x = tree_of[x];
        }
        x
    }
    for arc in &rg.arcs {
        assert!(arc.source < arc.target, "arcs run from low to high");
        assert!(seen_arcs.insert((arc.source, arc.target)), "duplicate arc");
        above_arcs[arc.source] += 1;
        below_arcs[arc.target] += 1;
        let (ra, rb) = (
            find(&mut tree_of, arc.source),
            find(&mut tree_of, arc.target),
        );
        assert_ne!(ra, rb, "the arcs close a cycle");
        tree_of[ra] = rb;
    }

    for v in 0..n {
        if adj[v].is_empty() {
            assert!(!node_of.contains_key(&v), "isolated vertex {v} reported");
            continue;
        }
        let lower = side_component_count(&adj, &rank, v, true);
        let upper = side_component_count(&adj, &rank, v, false);
        match node_of.get(&v) {
            Some(&i) => {
                assert_eq!(below_arcs[i], lower, "arcs below node of vertex {v}");
                assert_eq!(above_arcs[i], upper, "arcs above node of vertex {v}");
                let expected = if lower == 0 {
                    CriticalPointKind::Minimum
                } else if upper == 0 {
                    CriticalPointKind::Maximum
                } else {
                    CriticalPointKind::Saddle
                };
                assert_eq!(rg.nodes[i].kind, expected, "kind of vertex {v}");
                assert_eq!(rg.nodes[i].value, values[v]);
            }
            None => {
                assert_eq!((lower, upper), (1, 1), "vertex {v} is regular");
            }
        }
    }

    // One tree per connected component that has an edge.
    let all = slab_components(&adj, &rank, 0, n.saturating_sub(1));
    let components: HashSet<usize> = (0..n)
        .filter(|&v| !adj[v].is_empty())
        .map(|v| all[v])
        .collect();
    assert_eq!(rg.arcs.len() + components.len(), rg.nodes.len(), "forest");

    assert_arcs_preserve_level_set_connectivity(&adj, &rank, rg);
}

// ---------------------------------------------------------------------------
// Reeb graph (contour tree)
// ---------------------------------------------------------------------------

#[test]
fn reeb_monotone_ramp_has_one_minimum_one_maximum_one_arc() {
    let values = [0.0, 1.0, 2.0, 3.0, 4.0];
    let edges = [[0usize, 1], [1, 2], [2, 3], [3, 4]];
    let rg = ReebGraph::compute(&values, &edges);
    assert_eq!(rg.n_minima(), 1);
    assert_eq!(rg.n_maxima(), 1);
    assert_eq!(rg.n_saddles(), 0);
    assert_eq!(rg.arcs.len(), 1);
    assert_eq!(rg.nodes[0].vertex_idx, 0);
    assert_eq!(rg.nodes[1].vertex_idx, 4);
    assert_eq!((rg.arcs[0].source, rg.arcs[0].target), (0, 1));
    assert_contour_tree(&values, &edges, &rg);
}

#[test]
fn reeb_star_with_two_peaks_has_one_saddle_and_three_arcs() {
    // Centre 1 (value 1) joined to a = 0 (value 0), b = 2 (value 2), d = 3
    // (value 3): one minimum, one split saddle, two maxima.
    let values = [0.0, 1.0, 2.0, 3.0];
    let edges = [[1usize, 0], [1, 2], [1, 3]];
    let rg = ReebGraph::compute(&values, &edges);
    assert_eq!(rg.n_minima(), 1);
    assert_eq!(rg.n_saddles(), 1);
    assert_eq!(rg.n_maxima(), 2);
    assert_eq!(rg.arcs.len(), 3);
    let kinds: Vec<_> = rg
        .nodes
        .iter()
        .map(|n| (n.vertex_idx, n.kind.clone()))
        .collect();
    assert_eq!(
        kinds,
        vec![
            (0, CriticalPointKind::Minimum),
            (1, CriticalPointKind::Saddle),
            (2, CriticalPointKind::Maximum),
            (3, CriticalPointKind::Maximum),
        ]
    );
    let arcs: Vec<(usize, usize)> = rg.arcs.iter().map(|a| (a.source, a.target)).collect();
    assert_eq!(arcs, vec![(0, 1), (1, 2), (1, 3)]);
    assert_contour_tree(&values, &edges, &rg);
}

#[test]
fn reeb_join_saddle_of_two_minima() {
    // Two minima 0 and 2 (values 0 and 1) joined through vertex 1 (value 5)
    // which continues to the maximum 3 (value 9): a join saddle.
    let values = [0.0, 5.0, 1.0, 9.0];
    let edges = [[0usize, 1], [1, 2], [1, 3]];
    let rg = ReebGraph::compute(&values, &edges);
    assert_eq!(rg.n_minima(), 2);
    assert_eq!(rg.n_saddles(), 1);
    assert_eq!(rg.n_maxima(), 1);
    assert_eq!(rg.arcs.len(), 3);
    assert_contour_tree(&values, &edges, &rg);
}

#[test]
fn reeb_path_with_two_minima_has_a_maximum_with_two_arcs_below() {
    let values = [0.0, 0.5, 1.0, 0.7, 0.3];
    let edges = [[0usize, 1], [1, 2], [2, 3], [3, 4]];
    let rg = ReebGraph::compute(&values, &edges);
    assert_eq!(rg.n_minima(), 2);
    assert_eq!(rg.n_maxima(), 1);
    assert_eq!(rg.n_saddles(), 0);
    assert_eq!(rg.arcs.len(), 2);
    assert_contour_tree(&values, &edges, &rg);
}

#[test]
fn reeb_octahedron_height_has_two_critical_points() {
    // Vertices: 0 = -z, 1..=4 = equator (all at height 0), 5 = +z.
    let values = [-1.0, 0.0, 0.0, 0.0, 0.0, 1.0];
    let mut edges = Vec::new();
    for e in 1..=4usize {
        edges.push([0, e]);
        edges.push([e, 5]);
        edges.push([e, if e == 4 { 1 } else { e + 1 }]);
    }
    let rg = ReebGraph::compute(&values, &edges);
    assert_eq!(rg.n_minima(), 1);
    assert_eq!(rg.n_maxima(), 1);
    assert_eq!(rg.n_saddles(), 0);
    assert_eq!(rg.arcs.len(), 1);
    assert_contour_tree(&values, &edges, &rg);
}

#[test]
fn reeb_icosphere_with_two_peaks() {
    // f = x^2 + y on the unit sphere has two maxima at (+-0.866, 0.5, 0), one
    // minimum at (0, -1, 0) and one saddle at (0, 1, 0).
    let (vertices, edges) = icosphere(3);
    let values: Vec<f64> = vertices.iter().map(|p| p[0] * p[0] + p[1]).collect();
    let rg = ReebGraph::compute(&values, &edges);
    assert_eq!(rg.n_minima(), 1);
    assert_eq!(rg.n_maxima(), 2);
    assert_eq!(rg.n_saddles(), 1);
    assert_eq!(rg.arcs.len(), 3);
    let saddle = rg
        .nodes
        .iter()
        .find(|n| n.kind == CriticalPointKind::Saddle)
        .expect("one saddle");
    assert!(
        vertices[saddle.vertex_idx][1] > 0.9,
        "the saddle lies near +y"
    );
    assert_contour_tree(&values, &edges, &rg);
}

#[test]
fn reeb_x_squared_minus_y_squared_has_one_saddle_node_of_degree_four() {
    // f = x^2 - y^2 on the 5 x 5 grid around the origin: the interior saddle
    // at the centre has two arcs below and two above.
    let (nx, ny) = (5usize, 5usize);
    let mut values = Vec::new();
    for j in 0..ny {
        for i in 0..nx {
            let (x, y) = (i as f64 - 2.0, j as f64 - 2.0);
            values.push(x * x - y * y);
        }
    }
    let edges = grid_edges(nx, ny);
    let rg = ReebGraph::compute(&values, &edges);
    let saddles: Vec<_> = rg
        .nodes
        .iter()
        .enumerate()
        .filter(|(_, n)| n.kind == CriticalPointKind::Saddle)
        .collect();
    assert_eq!(saddles.len(), 1);
    let (index, node) = saddles[0];
    assert_eq!(node.vertex_idx, 12, "the centre of the grid");
    let below = rg.arcs.iter().filter(|a| a.target == index).count();
    let above = rg.arcs.iter().filter(|a| a.source == index).count();
    assert_eq!((below, above), (2, 2));
    assert_contour_tree(&values, &edges, &rg);
}

#[test]
fn reeb_random_values_on_triangulated_grids_satisfy_the_contour_tree_invariants() {
    let mut next = lcg(0x0c0f_fee1);
    for &(nx, ny) in &[(3usize, 3usize), (4, 3), (5, 4), (6, 6), (7, 5)] {
        let edges = grid_edges(nx, ny);
        for _ in 0..20 {
            let values: Vec<f64> = (0..nx * ny).map(|_| next()).collect();
            let rg = ReebGraph::compute(&values, &edges);
            assert_contour_tree(&values, &edges, &rg);
            assert_eq!(
                rg.arcs.len() + 1,
                rg.nodes.len(),
                "connected domain: a tree"
            );
        }
    }
}

#[test]
fn reeb_random_values_with_ties_satisfy_the_invariants() {
    // Values from a tiny alphabet create many plateaus.
    let mut next = lcg(0x7131_7131);
    let (nx, ny) = (5usize, 5usize);
    let edges = grid_edges(nx, ny);
    for _ in 0..40 {
        let values: Vec<f64> = (0..nx * ny).map(|_| (next() * 3.0).floor()).collect();
        let rg = ReebGraph::compute(&values, &edges);
        assert_contour_tree(&values, &edges, &rg);
        assert_eq!(rg.arcs.len() + 1, rg.nodes.len());
    }
}

#[test]
fn reeb_random_values_on_the_icosphere_satisfy_the_invariants() {
    let (vertices, edges) = icosphere(1);
    let mut next = lcg(0x1c05_a5f1);
    for _ in 0..30 {
        let values: Vec<f64> = (0..vertices.len()).map(|_| next()).collect();
        let rg = ReebGraph::compute(&values, &edges);
        assert_contour_tree(&values, &edges, &rg);
        assert_eq!(rg.arcs.len() + 1, rg.nodes.len(), "sphere: a tree");
    }
}

#[test]
fn reeb_constant_function_is_a_monotone_chain() {
    // All values equal: the order is by vertex index, a ramp along the path.
    let values = [2.0; 5];
    let edges = [[0usize, 1], [1, 2], [2, 3], [3, 4]];
    let rg = ReebGraph::compute(&values, &edges);
    assert_eq!(rg.n_minima(), 1);
    assert_eq!(rg.n_maxima(), 1);
    assert_eq!(rg.n_saddles(), 0);
    assert_eq!(rg.arcs.len(), 1);
}

#[test]
fn reeb_two_components_give_two_trees() {
    let values = [0.0, 1.0, 2.0, 5.0, 4.0, 3.0];
    let edges = [[0usize, 1], [1, 2], [3, 4], [4, 5]];
    let rg = ReebGraph::compute(&values, &edges);
    assert_eq!(rg.n_minima(), 2);
    assert_eq!(rg.n_maxima(), 2);
    assert_eq!(rg.arcs.len(), 2);
    assert_contour_tree(&values, &edges, &rg);
}

#[test]
fn reeb_ignores_isolated_vertices_invalid_and_repeated_edges() {
    let values = [0.0, 1.0, 2.0, 9.0];
    let edges = [[0usize, 1], [1, 0], [1, 2], [2, 2], [2, 17], [40, 41]];
    let rg = ReebGraph::compute(&values, &edges);
    assert_eq!(rg.nodes.len(), 2, "vertex 3 has no edge");
    assert_eq!(rg.arcs.len(), 1);
    assert_eq!(rg.nodes[0].vertex_idx, 0);
    assert_eq!(rg.nodes[1].vertex_idx, 2);
}

#[test]
fn reeb_nan_values_do_not_panic_and_sort_last() {
    let values = [0.0, f64::NAN, 2.0, 1.0];
    let edges = [[0usize, 1], [1, 2], [2, 3]];
    let rg = ReebGraph::compute(&values, &edges);
    assert_eq!(rg.arcs.len() + 1, rg.nodes.len());
    let top = rg.nodes.last().expect("nodes");
    assert_eq!(top.kind, CriticalPointKind::Maximum);
    assert!(top.value.is_nan());
}

#[test]
fn reeb_empty_and_single_vertex() {
    let rg = ReebGraph::compute(&[], &[]);
    assert!(rg.nodes.is_empty() && rg.arcs.is_empty());
    let rg = ReebGraph::compute(&[3.0], &[]);
    assert!(rg.nodes.is_empty() && rg.arcs.is_empty());
}

#[test]
fn reeb_cycle_graph_yields_a_forest_without_hanging() {
    // A 5-cycle whose join and split trees cannot be merged into a tree: the
    // loop of the Reeb graph is not represented.
    let values = [1.0, 4.0, 0.0, 2.0, 3.0];
    let edges = [[0usize, 1], [1, 3], [3, 2], [2, 4], [4, 0]];
    let rg = ReebGraph::compute(&values, &edges);
    assert_eq!(rg.n_minima(), 2, "vertices 0 and 2 have no lower neighbour");
    assert_eq!(rg.n_maxima(), 2, "vertices 1 and 4 have no upper neighbour");
    assert!(
        rg.arcs.len() < rg.nodes.len(),
        "a forest has fewer arcs than nodes"
    );
    for arc in &rg.arcs {
        assert!(arc.source < arc.target);
    }
}

#[test]
fn reeb_torus_height_function_collapses_to_one_arc() {
    // The Reeb graph of a height function on a torus has a loop; the contour
    // tree sees a single component of every sublevel set, so only the minimum
    // and the maximum remain, joined by one arc.
    let (nx, ny) = (24usize, 12usize);
    let values = tilted_torus_height(nx, ny);
    let edges = torus_edges(nx, ny);
    let rg = ReebGraph::compute(&values, &edges);
    assert_eq!(rg.nodes.len(), 2);
    assert_eq!(rg.arcs.len(), 1);
    assert_eq!(rg.n_minima(), 1);
    assert_eq!(rg.n_maxima(), 1);
}

// ---------------------------------------------------------------------------
// Morse complex
// ---------------------------------------------------------------------------

fn grid_values(nx: usize, ny: usize, f: impl Fn(f64, f64) -> f64) -> Vec<f64> {
    let mut values = Vec::with_capacity(nx * ny);
    for j in 0..ny {
        for i in 0..nx {
            values.push(f(
                i as f64 - (nx as f64 - 1.0) / 2.0,
                j as f64 - (ny as f64 - 1.0) / 2.0,
            ));
        }
    }
    values
}

#[test]
fn morse_linear_ramp_has_no_saddles() {
    let values = [0.0, 1.0, 2.0, 3.0, 4.0];
    let edges = [[0usize, 1], [1, 2], [2, 3], [3, 4]];
    let mc = MorseComplex::compute(&values, &edges);
    assert_eq!(mc.n_saddles(), 0);
    assert_eq!(mc.n_minima(), 1);
    assert_eq!(mc.n_maxima(), 1);
}

#[test]
fn morse_ramp_on_a_triangulated_grid_has_one_minimum_one_maximum_and_no_saddle() {
    let (nx, ny) = (6usize, 5usize);
    let edges = grid_edges(nx, ny);
    for values in [
        grid_values(nx, ny, |x, y| 2.0 * x + y),
        grid_values(nx, ny, |x, _| x),
        grid_values(nx, ny, |_, y| y),
        grid_values(nx, ny, |x, y| 3.0 * x - 0.5 * y),
    ] {
        let mc = MorseComplex::compute(&values, &edges);
        assert_eq!(mc.n_saddles(), 0, "a ramp has no saddle");
        assert_eq!(mc.n_minima(), 1);
        assert_eq!(mc.n_maxima(), 1);
        assert!(mc.connections.is_empty());
    }
}

#[test]
fn morse_x_squared_minus_y_squared_grid_has_exactly_one_saddle() {
    let (nx, ny) = (5usize, 5usize);
    let values = grid_values(nx, ny, |x, y| x * x - y * y);
    let mc = MorseComplex::compute(&values, &grid_edges(nx, ny));
    assert_eq!(mc.n_saddles(), 1);
    let saddle = mc
        .critical_points
        .iter()
        .find(|c| c.morse_index == 1)
        .expect("one saddle");
    assert_eq!(saddle.vertex_idx, 12, "the saddle is the grid centre");
    assert_eq!(saddle.multiplicity, 1);
    assert_eq!(saddle.value, 0.0);
}

#[test]
fn morse_monkey_saddle_counts_twice() {
    // Re((x + iy)^3) = x^3 - 3 x y^2 has a monkey saddle at the origin.
    let (nx, ny) = (5usize, 5usize);
    let values = grid_values(nx, ny, |x, y| x * x * x - 3.0 * x * y * y);
    let mc = MorseComplex::compute(&values, &grid_edges(nx, ny));
    let centre = mc
        .critical_points
        .iter()
        .find(|c| c.vertex_idx == 12)
        .expect("the centre is critical");
    assert_eq!(centre.morse_index, 1);
    assert_eq!(centre.multiplicity, 2);
}

#[test]
fn morse_spherical_function_satisfies_the_euler_relation() {
    let (vertices, edges) = icosphere(3);
    let values: Vec<f64> = vertices.iter().map(|p| p[0] * p[0] + p[1]).collect();
    let mc = MorseComplex::compute(&values, &edges);
    assert_eq!(mc.n_minima(), 1);
    assert_eq!(mc.n_maxima(), 2);
    assert_eq!(mc.n_saddles(), 1);
    assert_eq!(
        mc.n_minima() as i64 - mc.n_saddles() as i64 + mc.n_maxima() as i64,
        2,
        "chi of the sphere"
    );
    // The saddle is connected to the minimum (reached from both lower
    // components) and to each of the two maxima: three distinct separatrices.
    let saddle_index = mc
        .critical_points
        .iter()
        .position(|c| c.morse_index == 1)
        .expect("saddle");
    assert_eq!(mc.connections.len(), 3);
    assert!(mc.connections.iter().all(|&(s, _)| s == saddle_index));
    let mut kinds: Vec<usize> = mc
        .connections
        .iter()
        .map(|&(_, e)| mc.critical_points[e].morse_index)
        .collect();
    kinds.sort_unstable();
    assert_eq!(kinds, vec![0, 2, 2]);
}

#[test]
fn morse_random_functions_on_closed_surfaces_satisfy_the_euler_relation() {
    let mut next = lcg(0x0e01_e7a1);
    let (vertices, sphere_edges) = icosphere(2);
    let sphere_chi = clique_euler_characteristic(vertices.len(), &sphere_edges);
    assert_eq!(sphere_chi, 2);
    let (nx, ny) = (7usize, 6usize);
    let torus = torus_edges(nx, ny);
    let torus_chi = clique_euler_characteristic(nx * ny, &torus);
    assert_eq!(torus_chi, 0);
    for _ in 0..25 {
        let values: Vec<f64> = (0..vertices.len()).map(|_| next()).collect();
        let mc = MorseComplex::compute(&values, &sphere_edges);
        assert_eq!(
            mc.n_minima() as i64 - mc.n_saddles() as i64 + mc.n_maxima() as i64,
            sphere_chi
        );
        let values: Vec<f64> = (0..nx * ny).map(|_| next()).collect();
        let mc = MorseComplex::compute(&values, &torus);
        assert_eq!(
            mc.n_minima() as i64 - mc.n_saddles() as i64 + mc.n_maxima() as i64,
            torus_chi
        );
    }
}

#[test]
fn morse_torus_height_function_has_two_saddles() {
    let (nx, ny) = (24usize, 12usize);
    let values = tilted_torus_height(nx, ny);
    let mc = MorseComplex::compute(&values, &torus_edges(nx, ny));
    assert_eq!(mc.n_minima(), 1);
    assert_eq!(mc.n_maxima(), 1);
    assert_eq!(mc.n_saddles(), 2, "a torus needs two saddles");
    assert_eq!(
        mc.n_minima() as i64 - mc.n_saddles() as i64 + mc.n_maxima() as i64,
        0,
        "chi of the torus"
    );
}

#[test]
fn morse_classification_is_symmetric_under_negation() {
    let mut next = lcg(0x5a5a_1234);
    let (nx, ny) = (6usize, 5usize);
    let edges = grid_edges(nx, ny);
    for _ in 0..40 {
        let values: Vec<f64> = (0..nx * ny).map(|_| next()).collect();
        let negated: Vec<f64> = values.iter().map(|v| -v).collect();
        let a = MorseComplex::compute(&values, &edges);
        let b = MorseComplex::compute(&negated, &edges);
        assert_eq!(a.n_minima(), b.n_maxima());
        assert_eq!(a.n_maxima(), b.n_minima());
        assert_eq!(a.n_saddles(), b.n_saddles());
        let saddle_vertices = |m: &MorseComplex| -> Vec<usize> {
            let mut v: Vec<usize> = m
                .critical_points
                .iter()
                .filter(|c| c.morse_index == 1)
                .map(|c| c.vertex_idx)
                .collect();
            v.sort_unstable();
            v
        };
        assert_eq!(saddle_vertices(&a), saddle_vertices(&b));
    }
}

#[test]
fn morse_boundary_vertex_with_split_upper_link_is_a_saddle() {
    // Fan around vertex 0 (value 0): neighbours 1..=3 form a path of triangles
    // (0,1,2) and (0,2,3).  Values: 1 -> +1, 2 -> -1, 3 -> +1, so the lower
    // link {2} is one component but the upper link {1, 3} has two.
    let values = [0.0, 1.0, -1.0, 1.0];
    let edges = [[0usize, 1], [0, 2], [0, 3], [1, 2], [2, 3]];
    let mc = MorseComplex::compute(&values, &edges);
    let centre = mc
        .critical_points
        .iter()
        .find(|c| c.vertex_idx == 0)
        .expect("vertex 0 is critical");
    assert_eq!(centre.morse_index, 1);
    assert_eq!(centre.multiplicity, 1);
}

#[test]
fn morse_links_use_edges_between_neighbours_only() {
    // Without the edge between the two lower neighbours they are separate link
    // components; with it they form one.
    let values = [1.0, 0.0, 0.0, 2.0];
    let separate = MorseComplex::compute(&values, &[[0usize, 1], [0, 2], [0, 3]]);
    assert_eq!(separate.n_saddles(), 1);
    let linked = MorseComplex::compute(&values, &[[0usize, 1], [0, 2], [0, 3], [1, 2]]);
    assert_eq!(linked.n_saddles(), 0);
}

#[test]
fn morse_ignores_isolated_vertices_and_invalid_edges() {
    let values = [0.0, 1.0, 5.0];
    let mc = MorseComplex::compute(&values, &[[0usize, 1], [1, 1], [0, 9]]);
    assert_eq!(mc.critical_points.len(), 2, "vertex 2 is isolated");
    assert_eq!(mc.n_minima(), 1);
    assert_eq!(mc.n_maxima(), 1);
}

#[test]
fn morse_critical_points_are_sorted_by_value() {
    let mut next = lcg(0x1357_9bdf);
    let (nx, ny) = (6usize, 6usize);
    let edges = grid_edges(nx, ny);
    let values: Vec<f64> = (0..nx * ny).map(|_| next()).collect();
    let mc = MorseComplex::compute(&values, &edges);
    for pair in mc.critical_points.windows(2) {
        assert!(pair[0].value <= pair[1].value);
    }
    for &(saddle, extremum) in &mc.connections {
        assert_eq!(mc.critical_points[saddle].morse_index, 1);
        assert_ne!(mc.critical_points[extremum].morse_index, 1);
    }
    let unique: HashSet<_> = mc.connections.iter().collect();
    assert_eq!(unique.len(), mc.connections.len());
}

// ---------------------------------------------------------------------------
// A two-peak function on a grid
// ---------------------------------------------------------------------------

/// Values on the 5 x 3 grid of [`grid_edges`]: the boundary ring rises from 0
/// to 11 counter-clockwise from the corner `(0, 0)`, the interior holds two
/// hills, `(1, 1)` = 20 and `(3, 1)` = 21, and the pass `(2, 1)` = 15 between
/// them.  The hills stand on the ring, so the function has one minimum, two
/// maxima and the pass.
fn two_peak_grid_values() -> Vec<f64> {
    let nx = 5usize;
    let ring: [(usize, usize); 12] = [
        (0, 0),
        (1, 0),
        (2, 0),
        (3, 0),
        (4, 0),
        (4, 1),
        (4, 2),
        (3, 2),
        (2, 2),
        (1, 2),
        (0, 2),
        (0, 1),
    ];
    let mut values = vec![0.0; nx * 3];
    for (k, &(i, j)) in ring.iter().enumerate() {
        values[j * nx + i] = k as f64;
    }
    values[nx + 1] = 20.0;
    values[nx + 2] = 15.0;
    values[nx + 3] = 21.0;
    values
}

#[test]
fn reeb_two_peak_grid_function_has_two_maxima_one_saddle_one_minimum_and_three_arcs() {
    let values = two_peak_grid_values();
    let edges = grid_edges(5, 3);
    let rg = ReebGraph::compute(&values, &edges);
    assert_eq!(
        (rg.n_minima(), rg.n_saddles(), rg.n_maxima()),
        (1, 1, 2),
        "{:?}",
        rg.nodes
    );
    let nodes: Vec<(usize, CriticalPointKind)> = rg
        .nodes
        .iter()
        .map(|n| (n.vertex_idx, n.kind.clone()))
        .collect();
    // Vertex 0 = (0, 0), 7 = (2, 1), 6 = (1, 1), 8 = (3, 1).
    assert_eq!(
        nodes,
        vec![
            (0, CriticalPointKind::Minimum),
            (7, CriticalPointKind::Saddle),
            (6, CriticalPointKind::Maximum),
            (8, CriticalPointKind::Maximum),
        ]
    );
    let arcs: Vec<(usize, usize)> = rg.arcs.iter().map(|a| (a.source, a.target)).collect();
    assert_eq!(arcs, vec![(0, 1), (1, 2), (1, 3)]);
    assert_contour_tree(&values, &edges, &rg);
}

#[test]
fn morse_two_peak_grid_function_counts_the_pass_and_the_closing_of_the_ring() {
    // The ring of low values closes at vertex 5 = (0, 1): its lower neighbours
    // (0, 0) and {(0, 2), (1, 2)} are separate in its link, so it is a saddle of
    // the Morse complex, while the Reeb graph (components of level sets) does
    // not see it.  The other saddle is the pass between the hills.
    let values = two_peak_grid_values();
    let mc = MorseComplex::compute(&values, &grid_edges(5, 3));
    let kinds: Vec<(usize, usize)> = mc
        .critical_points
        .iter()
        .map(|c| (c.vertex_idx, c.morse_index))
        .collect();
    assert_eq!(kinds, vec![(0, 0), (5, 1), (7, 1), (6, 2), (8, 2)]);
    assert_eq!(mc.n_saddles(), 2);
    // Both saddles reach the minimum; the ring closure reaches the hill (1, 1)
    // through its upper neighbour, the pass reaches both hills.
    assert_eq!(mc.connections, vec![(1, 0), (1, 3), (2, 0), (2, 3), (2, 4)]);
}

#[test]
fn reeb_negative_nan_sorts_last_like_positive_nan() {
    let negative_nan = f64::from_bits(0xfff8_0000_0000_0000);
    assert!(negative_nan.is_nan() && negative_nan.is_sign_negative());
    let values = [0.0, negative_nan, 2.0, 1.0];
    let edges = [[0usize, 1], [1, 2], [2, 3]];
    let rg = ReebGraph::compute(&values, &edges);
    let positive = ReebGraph::compute(&[0.0, f64::NAN, 2.0, 1.0], &edges);
    let order = |g: &ReebGraph| -> Vec<usize> { g.nodes.iter().map(|n| n.vertex_idx).collect() };
    assert_eq!(order(&rg), order(&positive));
    let top = rg.nodes.last().expect("nodes");
    assert_eq!(
        (top.vertex_idx, top.kind.clone()),
        (1, CriticalPointKind::Maximum)
    );
    assert!(top.value.is_nan());
}
