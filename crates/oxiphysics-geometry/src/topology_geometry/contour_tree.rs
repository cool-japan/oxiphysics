// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Contour tree of a scalar function on a mesh graph.
//!
//! The contour tree is built as in Carr, Snoeyink and Axen, "Computing contour
//! trees in all dimensions" (2003):
//!
//! 1. The **join tree** follows the connected components of the sublevel sets
//!    `{f <= c}` while `c` rises: sweeping the vertices upwards, a vertex
//!    joins the components of its lower neighbours.  Its leaves are the local
//!    minima, and a vertex that merges two or more components is a join saddle.
//! 2. The **split tree** does the same for the superlevel sets `{f >= c}`
//!    while `c` falls; its leaves are the local maxima.
//! 3. Both trees are *augmented* (every vertex is a node).  They are merged by
//!    peeling leaves: a vertex with no lower arc and one upper arc, or the
//!    reverse, is a leaf of the contour tree.  Its arc joins it to its parent in
//!    the join tree (or split tree); the vertex is then deleted from both
//!    trees (a leaf of one, a vertex with one child in the other, which is
//!    spliced out).
//! 4. Vertices with exactly one arc below and one above are regular and are
//!    removed by joining their two arcs.  What remains are the critical points.
//!
//! For a function on a simply connected domain (a disc, a sphere, a ball) the
//! contour tree is the Reeb graph.  Where the domain has loops (the graph has
//! cycles that are not filled by triangles, or the surface has genus) the Reeb
//! graph has cycles, but join and split trees cannot see them.  The merge then
//! may not be able to resolve every vertex; the vertices left over are
//! connected along the join tree so that the result is always a forest.

use super::graph::ScalarGraph;
use super::{CriticalPointKind, ReebArc, ReebGraph, ReebNode, UnionFind};
use std::collections::VecDeque;

/// A sweep tree rooted at the last vertex of its sweep.
struct SweepTree {
    /// The vertex at which the component of a vertex merges next (higher
    /// for the join tree, lower for the split tree); `None` for the root.
    parent: Vec<Option<usize>>,
    /// Vertices that have (or once had) this vertex as parent.  Entries become
    /// stale when a vertex is deleted or re-parented; an entry `c` of `v` is
    /// live exactly when `c` is alive and `parent[c] == Some(v)`.
    children: Vec<Vec<usize>>,
    /// Number of live children of every vertex.
    degree: Vec<usize>,
}

/// Build the join tree (`ascending`) or the split tree of `graph`.
fn sweep_tree(graph: &ScalarGraph, ascending: bool) -> SweepTree {
    let n = graph.len();
    let mut tree = SweepTree {
        parent: vec![None; n],
        children: vec![Vec::new(); n],
        degree: vec![0; n],
    };
    let mut components = UnionFind::new(n);
    // `head[root]` is the vertex the component with this union-find root was
    // last extended by: the highest (join) or lowest (split) vertex so far.
    let mut head: Vec<usize> = (0..n).collect();
    let sweep: Vec<usize> = if ascending {
        graph.order.clone()
    } else {
        graph.order.iter().rev().copied().collect()
    };
    for v in sweep {
        for &u in &graph.neighbors[v] {
            let earlier = if ascending {
                graph.is_below(u, v)
            } else {
                graph.is_below(v, u)
            };
            if !earlier {
                continue;
            }
            let root_u = components.find(u);
            let root_v = components.find(v);
            if root_u != root_v {
                let tip = head[root_u];
                tree.parent[tip] = Some(v);
                tree.children[v].push(tip);
                tree.degree[v] += 1;
                components.union(root_u, root_v);
                let root = components.find(v);
                head[root] = v;
            }
        }
    }
    tree
}

/// The single live child of `v` in `tree`.
fn only_live_child(tree: &SweepTree, alive: &[bool], v: usize) -> Option<usize> {
    tree.children[v]
        .iter()
        .copied()
        .find(|&c| alive[c] && tree.parent[c] == Some(v))
}

/// Delete the vertex `v` that has exactly one live child `c` from `tree`:
/// `c` takes over the place of `v` below its parent.
fn splice_out(tree: &mut SweepTree, alive: &[bool], v: usize) {
    let Some(child) = only_live_child(tree, alive, v) else {
        return;
    };
    let parent = tree.parent[v];
    tree.parent[child] = parent;
    if let Some(p) = parent {
        // `v` stays in the list of `p` as a stale entry and `p` keeps the same
        // number of live children.
        tree.children[p].push(child);
    }
}

/// Arcs `(lower, upper)` of the augmented contour tree (a forest with one tree
/// per connected component of the graph).
fn augmented_arcs(graph: &ScalarGraph) -> Vec<(usize, usize)> {
    let n = graph.len();
    let mut join = sweep_tree(graph, true);
    let mut split = sweep_tree(graph, false);
    let mut alive = vec![true; n];
    let mut arcs: Vec<(usize, usize)> = Vec::with_capacity(n.saturating_sub(1));
    let mut queue: VecDeque<usize> = (0..n)
        .filter(|&v| join.degree[v] + split.degree[v] == 1)
        .collect();

    while let Some(v) = queue.pop_front() {
        if !alive[v] || join.degree[v] + split.degree[v] != 1 {
            continue;
        }
        let neighbour = if join.degree[v] == 0 {
            // Lowest vertex of its arc: it continues upwards into its join
            // tree parent.
            let Some(w) = join.parent[v] else { continue };
            arcs.push((v, w));
            join.degree[w] -= 1;
            splice_out(&mut split, &alive, v);
            w
        } else {
            // Highest vertex of its arc: it continues downwards into its
            // split tree parent.
            let Some(w) = split.parent[v] else { continue };
            arcs.push((w, v));
            split.degree[w] -= 1;
            splice_out(&mut join, &alive, v);
            w
        };
        alive[v] = false;
        if join.degree[neighbour] + split.degree[neighbour] == 1 {
            queue.push_back(neighbour);
        }
    }

    // Vertices the merge could not resolve (only for domains with loops) are
    // attached along the remaining join tree.
    for (v, &still_alive) in alive.iter().enumerate() {
        if still_alive && let Some(p) = join.parent[v] {
            arcs.push((v, p));
        }
    }
    arcs
}

/// Compute the contour tree of `values` on the graph `edges`; see the module
/// documentation.
pub(super) fn compute(values: &[f64], edges: &[[usize; 2]]) -> ReebGraph {
    let graph = ScalarGraph::new(values, edges);
    let n = graph.len();
    let arcs = augmented_arcs(&graph);

    // Arc counts below and above every vertex classify it.
    let mut below: Vec<Vec<usize>> = vec![Vec::new(); n];
    let mut above: Vec<Vec<usize>> = vec![Vec::new(); n];
    for &(lower, upper) in &arcs {
        above[lower].push(upper);
        below[upper].push(lower);
    }
    let kind_of = |v: usize| -> Option<CriticalPointKind> {
        match (below[v].len(), above[v].len()) {
            (0, 0) | (1, 1) => None,
            (0, _) => Some(CriticalPointKind::Minimum),
            (_, 0) => Some(CriticalPointKind::Maximum),
            _ => Some(CriticalPointKind::Saddle),
        }
    };

    // Nodes in ascending order of the function.
    let mut node_of: Vec<Option<usize>> = vec![None; n];
    let mut nodes: Vec<ReebNode> = Vec::new();
    for &v in &graph.order {
        if let Some(kind) = kind_of(v) {
            node_of[v] = Some(nodes.len());
            nodes.push(ReebNode {
                value: values[v],
                kind,
                vertex_idx: v,
            });
        }
    }

    // Arcs between critical points: follow every upward arc through regular
    // vertices, which have exactly one arc above.
    let mut reeb_arcs: Vec<ReebArc> = Vec::new();
    for node in &nodes {
        let start = node.vertex_idx;
        for &first in &above[start] {
            let mut current = first;
            while node_of[current].is_none() {
                match above[current].first() {
                    Some(&next) => current = next,
                    None => break,
                }
            }
            if let (Some(source), Some(target)) = (node_of[start], node_of[current]) {
                reeb_arcs.push(ReebArc { source, target });
            }
        }
    }
    reeb_arcs.sort_unstable_by_key(|arc| (arc.source, arc.target));
    ReebGraph {
        nodes,
        arcs: reeb_arcs,
    }
}
