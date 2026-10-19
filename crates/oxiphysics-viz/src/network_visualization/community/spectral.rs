// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Spectral clustering by recursive bisection with the Fiedler vector.
//!
//! The Fiedler vector of a graph is the eigenvector of the combinatorial
//! Laplacian `L = D - A` that belongs to its second-smallest eigenvalue; the
//! smallest one, `0`, belongs to the constant vector.  Its entries order the
//! nodes along the sparsest cut of the graph.  [`fiedler`] computes it.
//!
//! # Cutting
//!
//! The nodes are sorted by their Fiedler entry and the prefix of this order
//! with the lowest *conductance* `cut(S) / min(vol(S), vol(V \ S))` becomes one
//! side of the bisection (Cheeger's sweep cut).  The split at the sign change
//! of the Fiedler vector is one of the candidates, so the sweep is never worse
//! than splitting by sign; unlike the sign split it also follows the bridges
//! when a cluster in the middle of a chain has Fiedler entries near zero.
//!
//! # Refinement
//!
//! Sorting by one eigenvector flattens the graph to a line, and a node that is
//! weakly attached to a far-away part of the graph is pulled towards it and may
//! land on the wrong side of the cut.  After the sweep, single nodes are
//! therefore moved to the other side while that lowers the conductance, the
//! greedy core of the Kernighan-Lin / Fiduccia-Mattheyses refinement that
//! spectral partitioners use.
//!
//! For `k > 2` the largest community is bisected repeatedly, each time using
//! the Fiedler vector of the subgraph it induces.

use super::{WeightedGraph, densify};

mod fiedler;
mod jacobi;

#[cfg(test)]
mod tests;

/// Conductances closer than this count as tied; the more balanced cut wins.
const CONDUCTANCE_TIE: f64 = 1e-9;

/// Bound on the refinement passes over all nodes of one bisection.  Every
/// accepted move lowers the conductance, so the passes end on their own long
/// before; the bound only guards against pathological rounding.
const MAX_REFINEMENT_PASSES: usize = 50;

/// Subgraph induced by a set of nodes, prepared for Laplacian products.
///
/// Node `i` of the subgraph is the `i`-th member passed to
/// [`Subgraph::induced`].  Self-loops are dropped; they cancel in `D - A`.
struct Subgraph {
    /// Neighbours inside the subgraph as `(local index, weight)`.
    adj: Vec<Vec<(usize, f64)>>,
    /// Weighted degree inside the subgraph.
    degree: Vec<f64>,
}

impl Subgraph {
    /// Subgraph of `graph` induced by `members`.  `local_of` is scratch space
    /// with one slot per node of `graph`, all `usize::MAX` on entry and on
    /// exit.
    fn induced(graph: &WeightedGraph, members: &[usize], local_of: &mut [usize]) -> Self {
        for (local, &global) in members.iter().enumerate() {
            local_of[global] = local;
        }
        let adj: Vec<Vec<(usize, f64)>> = members
            .iter()
            .map(|&global| {
                graph.adj[global]
                    .iter()
                    .filter_map(|&(other, weight)| {
                        let local = local_of[other];
                        (local != usize::MAX).then_some((local, weight))
                    })
                    .collect()
            })
            .collect();
        for &global in members {
            local_of[global] = usize::MAX;
        }
        let degree = adj
            .iter()
            .map(|row| row.iter().map(|&(_, weight)| weight).sum())
            .collect();
        Self { adj, degree }
    }

    /// Number of nodes.
    fn len(&self) -> usize {
        self.degree.len()
    }

    /// `out = L x` for the Laplacian `L = D - A` of the subgraph.
    fn laplacian_product(&self, x: &[f64], out: &mut [f64]) {
        for (i, slot) in out.iter_mut().enumerate() {
            let mut acc = self.degree[i] * x[i];
            for &(j, weight) in &self.adj[i] {
                acc -= weight * x[j];
            }
            *slot = acc;
        }
    }
}

/// Partition `graph` into `k` communities by recursive spectral bisection and
/// return one dense label per node, numbered in order of first appearance.
///
/// `k >= n` yields singletons.  Callers handle `k == 0`.
pub(super) fn spectral_partition(graph: &WeightedGraph, k: usize) -> Vec<usize> {
    let n = graph.len();
    if k >= n {
        return (0..n).collect();
    }
    // Every part lists its nodes in ascending index order.
    let mut parts: Vec<Vec<usize>> = vec![(0..n).collect()];
    let mut local_of = vec![usize::MAX; n];
    while parts.len() < k {
        let Some(target) = largest_splittable(&parts) else {
            break;
        };
        let members = std::mem::take(&mut parts[target]);
        let (lower, upper) = bisect(graph, &members, &mut local_of);
        parts[target] = lower;
        parts.push(upper);
    }
    let mut labels = vec![0_usize; n];
    for (label, part) in parts.iter().enumerate() {
        for &node in part {
            labels[node] = label;
        }
    }
    densify(&mut labels);
    labels
}

/// Index of the largest part with at least two nodes; the earlier part wins a
/// tie.
fn largest_splittable(parts: &[Vec<usize>]) -> Option<usize> {
    let mut best: Option<usize> = None;
    for (index, part) in parts.iter().enumerate() {
        if part.len() < 2 {
            continue;
        }
        let larger = best.is_none_or(|current| part.len() > parts[current].len());
        if larger {
            best = Some(index);
        }
    }
    best
}

/// Split `members` (at least two nodes, ascending) into two non-empty parts,
/// each in ascending order.
fn bisect(
    graph: &WeightedGraph,
    members: &[usize],
    local_of: &mut [usize],
) -> (Vec<usize>, Vec<usize>) {
    let sub = Subgraph::induced(graph, members, local_of);
    let size = sub.len();
    let mut in_lower = vec![false; size];
    if sub.degree.iter().sum::<f64>() > 0.0 {
        let vector = fiedler::fiedler_vector(&sub);
        let mut order: Vec<usize> = (0..size).collect();
        order.sort_by(|&a, &b| vector[a].total_cmp(&vector[b]).then(a.cmp(&b)));
        let lower_len = best_sweep_prefix(&sub, &order).unwrap_or(size / 2);
        for &local in &order[..lower_len] {
            in_lower[local] = true;
        }
        refine_bisection(&sub, &mut in_lower);
    } else {
        // No edges, so no structure to follow: halve by node index.
        for flag in in_lower.iter_mut().take(size / 2) {
            *flag = true;
        }
    }
    let mut lower = Vec::with_capacity(size);
    let mut upper = Vec::with_capacity(size);
    for (local, &global) in members.iter().enumerate() {
        if in_lower[local] {
            lower.push(global);
        } else {
            upper.push(global);
        }
    }
    (lower, upper)
}

/// Length of the prefix of `order` whose node set has the lowest conductance
/// among the proper prefixes; ties go to the prefix whose smaller side has the
/// larger volume, then to the shorter prefix.
///
/// `None` when no prefix separates positive volume from positive volume, which
/// only happens in a graph without edges.
fn best_sweep_prefix(sub: &Subgraph, order: &[usize]) -> Option<usize> {
    let total: f64 = sub.degree.iter().sum();
    let mut inside = vec![false; order.len()];
    let mut cut = 0.0_f64;
    let mut volume = 0.0_f64;
    // (conductance, volume of the smaller side, prefix length)
    let mut best: Option<(f64, f64, usize)> = None;
    for (taken, &node) in order.iter().enumerate().take(order.len().saturating_sub(1)) {
        // Adding `node` removes its edges into the prefix from the cut and
        // adds its edges to the rest.
        let into_prefix: f64 = sub.adj[node]
            .iter()
            .filter(|&&(other, _)| inside[other])
            .map(|&(_, weight)| weight)
            .sum();
        cut += sub.degree[node] - 2.0 * into_prefix;
        volume += sub.degree[node];
        inside[node] = true;
        let smaller = volume.min(total - volume);
        if smaller <= 0.0 {
            continue;
        }
        let conductance = cut.max(0.0) / smaller;
        let better = match best {
            None => true,
            Some((best_conductance, best_smaller, _)) => {
                conductance < best_conductance - CONDUCTANCE_TIE
                    || (conductance <= best_conductance + CONDUCTANCE_TIE && smaller > best_smaller)
            }
        };
        if better {
            best = Some((conductance, smaller, taken + 1));
        }
    }
    best.map(|(_, _, length)| length)
}

/// Conductance `cut / min(vol_a, vol_b)` of a bisection; infinite when a side
/// has no volume.
fn conductance(cut: f64, volume_a: f64, volume_b: f64) -> f64 {
    let smaller = volume_a.min(volume_b);
    if smaller > 0.0 {
        cut.max(0.0) / smaller
    } else {
        f64::INFINITY
    }
}

/// Improve the bisection `in_lower` by moving single nodes to the other side
/// while that lowers its conductance.
///
/// A one-dimensional sweep cannot be perfect: the Fiedler entry of a node that
/// is weakly attached to a far-away part of the graph is pulled towards that
/// part, and the node may end up on the wrong side of the cut.  Such a node has
/// more edge weight towards the other side than towards its own, and moving it
/// shrinks the cut.  This is the greedy core of the Kernighan-Lin /
/// Fiduccia-Mattheyses refinement that spectral partitioners apply after a
/// bisection.  Both sides always keep at least one node, and a bisection that
/// no single move can improve (in particular any cut of weight zero) is left
/// alone.
fn refine_bisection(sub: &Subgraph, in_lower: &mut [bool]) {
    let size = sub.len();
    let total: f64 = sub.degree.iter().sum();
    for _ in 0..MAX_REFINEMENT_PASSES {
        // Recompute the aggregates every pass so rounding error cannot build up.
        let mut cut = 0.0_f64;
        let mut volume_lower = 0.0_f64;
        let mut count_lower = 0_usize;
        for node in 0..size {
            if in_lower[node] {
                volume_lower += sub.degree[node];
                count_lower += 1;
            }
            for &(other, weight) in &sub.adj[node] {
                if in_lower[node] != in_lower[other] {
                    // Every cut edge is seen from both of its end points.
                    cut += 0.5 * weight;
                }
            }
        }
        let mut moved = false;
        for node in 0..size {
            let own_side = in_lower[node];
            let own_count = if own_side {
                count_lower
            } else {
                size - count_lower
            };
            if own_count <= 1 {
                continue;
            }
            let (mut to_own, mut to_other) = (0.0_f64, 0.0_f64);
            for &(other, weight) in &sub.adj[node] {
                if in_lower[other] == own_side {
                    to_own += weight;
                } else {
                    to_other += weight;
                }
            }
            let degree = sub.degree[node];
            let (new_volume_lower, new_count_lower) = if own_side {
                (volume_lower - degree, count_lower - 1)
            } else {
                (volume_lower + degree, count_lower + 1)
            };
            let new_cut = cut - to_other + to_own;
            let before = conductance(cut, volume_lower, total - volume_lower);
            let after = conductance(new_cut, new_volume_lower, total - new_volume_lower);
            if after < before - CONDUCTANCE_TIE {
                in_lower[node] = !own_side;
                cut = new_cut.max(0.0);
                volume_lower = new_volume_lower;
                count_lower = new_count_lower;
                moved = true;
            }
        }
        if !moved {
            break;
        }
    }
}
