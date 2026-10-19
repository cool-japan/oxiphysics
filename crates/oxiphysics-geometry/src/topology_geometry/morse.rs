// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Critical points and gradient connections of a scalar function on a mesh.
//!
//! The mesh is given by its edges and is read as the clique complex of the
//! graph: two neighbours of a vertex `v` are *linked* when they are joined by
//! an edge, which for the edges of a triangulation means `(v, a, b)` is a
//! triangle.  The link of `v` is then the graph on its neighbours with these
//! link edges, a cycle for an interior vertex of a triangulated surface and a
//! path for a boundary vertex.
//!
//! The neighbours of `v` split into the lower ones (before `v` in the vertex
//! order) and the upper ones.  The vertex is
//!
//! * a minimum when it has no lower neighbour,
//! * a maximum when it has no upper neighbour,
//! * regular when the lower neighbours form one component of the link and so
//!   do the upper neighbours,
//! * a saddle otherwise, of multiplicity `max(k_lower, k_upper) - 1` where
//!   `k_lower` and `k_upper` are the numbers of link components of the lower
//!   and upper neighbours.
//!
//! On a closed surface `k_lower == k_upper`, so the counts give the Euler
//! characteristic: `minima - saddles + maxima = V - E + F`.  Counting both
//! sides makes the classification symmetric under `f -> -f` also at boundary
//! vertices, where the link is a path and the two counts can differ.
//!
//! Every saddle is connected to extrema by the steepest descent path from each
//! of its lower link components (ending in a minimum) and the steepest ascent
//! path from each of its upper link components (ending in a maximum).

use super::graph::ScalarGraph;
use super::{MorseComplex, MorseCriticalPoint, UnionFind};

/// Marker in the scratch position table for "not a member of the link".
const NOT_A_MEMBER: usize = usize::MAX;

/// Link component label of every vertex in `members`, which are neighbours of
/// one vertex; two members share a label iff they are connected by edges
/// between members.  The label of a component is the index of one of its
/// members.  `position` is a scratch table of the size of the graph and is
/// all [`NOT_A_MEMBER`] on entry and exit.
fn link_components(graph: &ScalarGraph, members: &[usize], position: &mut [usize]) -> Vec<usize> {
    for (i, &m) in members.iter().enumerate() {
        position[m] = i;
    }
    let mut sets = UnionFind::new(members.len());
    for (i, &m) in members.iter().enumerate() {
        for &w in &graph.neighbors[m] {
            let j = position[w];
            if j != NOT_A_MEMBER {
                sets.union(i, j);
            }
        }
    }
    let labels: Vec<usize> = (0..members.len()).map(|i| sets.find(i)).collect();
    for &m in members {
        position[m] = NOT_A_MEMBER;
    }
    labels
}

/// Number of distinct labels.
fn component_count(labels: &[usize]) -> usize {
    labels.iter().enumerate().filter(|&(i, &l)| i == l).count()
}

/// The member of every component (in label order) that comes first (`lowest`)
/// or last in the vertex order: the start of the steepest path.
fn steepest_members(
    graph: &ScalarGraph,
    members: &[usize],
    labels: &[usize],
    lowest: bool,
) -> Vec<usize> {
    let mut best: Vec<Option<usize>> = vec![None; members.len()];
    for (i, &m) in members.iter().enumerate() {
        let slot = &mut best[labels[i]];
        let better = match *slot {
            None => true,
            Some(current) => {
                if lowest {
                    graph.is_below(m, current)
                } else {
                    graph.is_below(current, m)
                }
            }
        };
        if better {
            *slot = Some(m);
        }
    }
    best.into_iter().flatten().collect()
}

/// Compute the critical points and connections; see the module docs.
pub(super) fn compute(values: &[f64], edges: &[[usize; 2]]) -> MorseComplex {
    let graph = ScalarGraph::new(values, edges);
    let n = graph.len();

    // The extremum every vertex reaches by steepest descent / ascent.
    let mut descent_end: Vec<usize> = (0..n).collect();
    for &v in &graph.order {
        let lowest = graph.neighbors[v]
            .iter()
            .copied()
            .filter(|&u| graph.is_below(u, v))
            .min_by_key(|&u| graph.rank[u]);
        if let Some(u) = lowest {
            descent_end[v] = descent_end[u];
        }
    }
    let mut ascent_end: Vec<usize> = (0..n).collect();
    for &v in graph.order.iter().rev() {
        let highest = graph.neighbors[v]
            .iter()
            .copied()
            .filter(|&u| graph.is_below(v, u))
            .max_by_key(|&u| graph.rank[u]);
        if let Some(u) = highest {
            ascent_end[v] = ascent_end[u];
        }
    }

    let mut position = vec![NOT_A_MEMBER; n];
    let mut critical_points: Vec<MorseCriticalPoint> = Vec::new();
    let mut index_of: Vec<Option<usize>> = vec![None; n];
    // (saddle vertex, extremum vertex)
    let mut links: Vec<(usize, usize)> = Vec::new();

    // Ascending order, so the critical points come out sorted by value.
    for &v in &graph.order {
        let neighbors = &graph.neighbors[v];
        if neighbors.is_empty() {
            continue;
        }
        let lower: Vec<usize> = neighbors
            .iter()
            .copied()
            .filter(|&u| graph.is_below(u, v))
            .collect();
        let upper: Vec<usize> = neighbors
            .iter()
            .copied()
            .filter(|&u| graph.is_below(v, u))
            .collect();
        let (morse_index, multiplicity) = if lower.is_empty() {
            (0, 1)
        } else if upper.is_empty() {
            (2, 1)
        } else {
            let lower_labels = link_components(&graph, &lower, &mut position);
            let upper_labels = link_components(&graph, &upper, &mut position);
            let lower_count = component_count(&lower_labels);
            let upper_count = component_count(&upper_labels);
            if lower_count == 1 && upper_count == 1 {
                continue;
            }
            for start in steepest_members(&graph, &lower, &lower_labels, true) {
                links.push((v, descent_end[start]));
            }
            for start in steepest_members(&graph, &upper, &upper_labels, false) {
                links.push((v, ascent_end[start]));
            }
            (1, lower_count.max(upper_count) - 1)
        };
        index_of[v] = Some(critical_points.len());
        critical_points.push(MorseCriticalPoint {
            vertex_idx: v,
            value: values[v],
            morse_index,
            multiplicity,
        });
    }

    let mut connections: Vec<(usize, usize)> = links
        .into_iter()
        .filter_map(|(saddle, extremum)| Some((index_of[saddle]?, index_of[extremum]?)))
        .collect();
    connections.sort_unstable();
    connections.dedup();
    MorseComplex {
        critical_points,
        connections,
    }
}
