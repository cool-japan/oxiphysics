// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! The Louvain method (Blondel, Guillaume, Lambiotte, Lefebvre, 2008) for
//! modularity maximisation.
//!
//! The method alternates two phases until modularity stops improving:
//!
//! 1. **Local moving.**  Starting with every node in its own community, nodes
//!    are visited in index order and moved to the neighbouring community with
//!    the largest modularity gain, until a full sweep moves nothing.
//! 2. **Aggregation.**  Every community becomes a single node of a coarser
//!    graph (internal edges turn into a self-loop, edges between communities
//!    are merged), and phase 1 runs again on that graph.
//!
//! Moving an isolated node `i` of strength `k_i` into a community with total
//! strength `sigma_tot`, which `i` is joined to by edges of total weight
//! `k_i,in`, changes modularity by
//!
//! ```text
//! dQ = k_i,in / m - sigma_tot * k_i / (2 m^2)
//! ```
//!
//! so a move costs only `O(deg(i))` instead of a fresh `O(n^2)` modularity
//! evaluation.  Everything is deterministic: nodes are visited in index order,
//! candidate communities are examined in adjacency order and no hash map is
//! iterated.

use super::{WeightedGraph, densify};

#[cfg(test)]
mod tests;

/// A node changes community only if that raises modularity by more than this.
/// The strict threshold guarantees termination despite rounding error.
const MIN_GAIN: f64 = 1e-12;

/// Safety bound on the sweeps of one local-moving phase.  Convergence takes a
/// few dozen sweeps in practice; hitting the bound still yields a valid
/// partition whose modularity is no lower than where the phase started.
const MAX_SWEEPS: usize = 10_000;

/// Outcome of one local-moving phase.
struct LocalMoving {
    /// Dense community label of every node of the level.
    community: Vec<usize>,
    /// Number of communities.
    count: usize,
    /// Whether any node left its starting singleton community.
    improved: bool,
}

/// Louvain partition of `graph`: one dense label per node, numbered in order
/// of first appearance.
pub(super) fn louvain(graph: &WeightedGraph) -> Vec<usize> {
    let mut membership: Vec<usize> = (0..graph.len()).collect();
    if !(graph.total_weight.is_finite() && graph.total_weight > 0.0) {
        return membership;
    }
    // `None` stands for the input graph itself, so it is never copied.
    let mut coarse: Option<WeightedGraph> = None;
    loop {
        let level = coarse.as_ref().unwrap_or(graph);
        let step = local_moving(level);
        if !step.improved {
            break;
        }
        // `membership` maps original nodes to nodes of `level`; compose it
        // with the communities just found.
        for label in membership.iter_mut() {
            *label = step.community[*label];
        }
        if step.count <= 1 {
            break;
        }
        // Moving at least one node always merges communities, so every level
        // is strictly smaller than the previous one and the loop terminates.
        let next = aggregate(level, &step.community, step.count);
        coarse = Some(next);
    }
    densify(&mut membership);
    membership
}

/// Phase 1: greedy node moves on `graph`, starting from singletons.
fn local_moving(graph: &WeightedGraph) -> LocalMoving {
    let n = graph.len();
    let m = graph.total_weight;
    let two_m = 2.0 * m;
    let mut community: Vec<usize> = (0..n).collect();
    // Total strength of every community.
    let mut sigma_tot: Vec<f64> = graph.strength.clone();
    // Scratch for the weight from the current node to each neighbouring
    // community.  `stamp[c]` records the visit that last wrote `weight_to[c]`,
    // so nothing has to be cleared between nodes.
    let mut weight_to = vec![0.0_f64; n];
    let mut stamp = vec![0_usize; n];
    let mut touched: Vec<usize> = Vec::new();
    let mut visit = 0_usize;
    let mut improved = false;

    for _ in 0..MAX_SWEEPS {
        let mut moved = false;
        for node in 0..n {
            let neighbours = &graph.adj[node];
            if neighbours.is_empty() {
                continue;
            }
            visit += 1;
            touched.clear();
            for &(other, weight) in neighbours {
                let target = community[other];
                if stamp[target] != visit {
                    stamp[target] = visit;
                    weight_to[target] = 0.0;
                    touched.push(target);
                }
                weight_to[target] += weight;
            }
            let strength = graph.strength[node];
            // k_i / m; with sigma_tot / 2m below it forms the null-model term
            // sigma_tot * k_i / (2 m^2) without squaring a possibly huge m.
            let relative_strength = strength / m;
            let old = community[node];
            // Take the node out of its community; staying is then just one
            // more candidate, evaluated with the same formula.
            sigma_tot[old] -= strength;
            let gain = |candidate: usize| -> f64 {
                let linked = if stamp[candidate] == visit {
                    weight_to[candidate]
                } else {
                    0.0
                };
                linked / m - sigma_tot[candidate] / two_m * relative_strength
            };
            let mut best = old;
            let mut best_gain = gain(old);
            for &candidate in &touched {
                if candidate == old {
                    continue;
                }
                let candidate_gain = gain(candidate);
                if candidate_gain > best_gain + MIN_GAIN {
                    best = candidate;
                    best_gain = candidate_gain;
                }
            }
            sigma_tot[best] += strength;
            if best != old {
                community[node] = best;
                moved = true;
            }
        }
        if !moved {
            break;
        }
        improved = true;
    }

    let count = densify(&mut community);
    LocalMoving {
        community,
        count,
        improved,
    }
}

/// Phase 2: collapse every community of `graph` into one node.
///
/// The self-loop of a super-node collects the self-loops of its members plus
/// the weight of the edges inside the community; its strength is the summed
/// strength of the members, so the total weight `m` is unchanged.
fn aggregate(graph: &WeightedGraph, community: &[usize], count: usize) -> WeightedGraph {
    let mut members: Vec<Vec<usize>> = vec![Vec::new(); count];
    for (node, &label) in community.iter().enumerate() {
        members[label].push(node);
    }
    let mut coarse = WeightedGraph::with_nodes(count);
    coarse.total_weight = graph.total_weight;
    let mut weight_to = vec![0.0_f64; count];
    let mut stamp = vec![usize::MAX; count];
    let mut touched: Vec<usize> = Vec::new();
    for (group, nodes) in members.iter().enumerate() {
        touched.clear();
        // Twice the weight of the edges inside the group: each is seen from
        // both of its end points.
        let mut inside = 0.0_f64;
        for &node in nodes {
            coarse.self_loop[group] += graph.self_loop[node];
            coarse.strength[group] += graph.strength[node];
            for &(other, weight) in &graph.adj[node] {
                let target = community[other];
                if target == group {
                    inside += weight;
                    continue;
                }
                if stamp[target] != group {
                    stamp[target] = group;
                    weight_to[target] = 0.0;
                    touched.push(target);
                }
                weight_to[target] += weight;
            }
        }
        coarse.self_loop[group] += 0.5 * inside;
        coarse.adj[group] = touched
            .iter()
            .map(|&target| (target, weight_to[target]))
            .collect();
    }
    coarse
}
