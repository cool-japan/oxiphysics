// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Greedy aggregation for smoothed-aggregation AMG.

#![allow(dead_code)]

/// Greedy aggregation of nodes based on strong connections.
///
/// Returns `aggregate_ids[i]` — the aggregate ID for each node `i`.
///
/// Algorithm:
/// 1. **Pass 1 (seed)**: For each unmarked node `i` whose all strong neighbors
///    are also unmarked, create a new aggregate containing `i` and all its
///    unmarked strong neighbors.
/// 2. **Pass 2 (assign)**: For remaining unmarked nodes, assign to the aggregate
///    of the first assigned strong neighbor.
/// 3. **Singletons**: Any still-unmarked node forms its own aggregate.
pub fn greedy_aggregate(strong: &[Vec<usize>], n: usize) -> Vec<usize> {
    const UNASSIGNED: usize = usize::MAX;
    let mut aggregate_ids = vec![UNASSIGNED; n];
    let mut next_agg = 0usize;

    // Pass 1: Seed aggregates
    // A node seeds if it and all its strong neighbors are unmarked
    for i in 0..n {
        if aggregate_ids[i] != UNASSIGNED {
            continue;
        }
        // Check all strong neighbors are unmarked
        let all_unmarked = strong[i].iter().all(|&j| aggregate_ids[j] == UNASSIGNED);
        if all_unmarked {
            let agg_id = next_agg;
            next_agg += 1;
            aggregate_ids[i] = agg_id;
            for &j in &strong[i] {
                if aggregate_ids[j] == UNASSIGNED {
                    aggregate_ids[j] = agg_id;
                }
            }
        }
    }

    // Pass 2: Assign remaining nodes to a neighbor's aggregate
    // Iterate multiple times to ensure propagation through chains
    let mut changed = true;
    while changed {
        changed = false;
        for i in 0..n {
            if aggregate_ids[i] == UNASSIGNED {
                // Look for any assigned strong neighbor
                for &j in &strong[i] {
                    if aggregate_ids[j] != UNASSIGNED {
                        aggregate_ids[i] = aggregate_ids[j];
                        changed = true;
                        break;
                    }
                }
            }
        }
    }

    // Pass 3: Singletons — any still-unassigned node forms its own aggregate
    for agg_id in aggregate_ids.iter_mut().take(n) {
        if *agg_id == UNASSIGNED {
            *agg_id = next_agg;
            next_agg += 1;
        }
    }

    aggregate_ids
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a chain graph (1D) strong connections.
    fn chain_strong(n: usize) -> Vec<Vec<usize>> {
        (0..n)
            .map(|i| {
                let mut nbrs = Vec::new();
                if i > 0 {
                    nbrs.push(i - 1);
                }
                if i + 1 < n {
                    nbrs.push(i + 1);
                }
                nbrs
            })
            .collect()
    }

    #[test]
    fn test_aggregation_covers_all_nodes() {
        let n = 12;
        let strong = chain_strong(n);
        let agg_ids = greedy_aggregate(&strong, n);

        assert_eq!(agg_ids.len(), n);

        // Every node must be assigned
        for (i, &id) in agg_ids.iter().enumerate() {
            assert!(
                id != usize::MAX,
                "Node {i} was not assigned to any aggregate"
            );
        }

        // Find the number of aggregates (must be contiguous starting from 0)
        let n_agg = *agg_ids.iter().max().unwrap_or(&0) + 1;

        // Every aggregate ID in 0..n_agg must appear at least once
        let mut seen = vec![false; n_agg];
        for &id in &agg_ids {
            seen[id] = true;
        }
        for (id, &s) in seen.iter().enumerate() {
            assert!(s, "Aggregate {id} has no members");
        }
    }

    #[test]
    fn test_aggregation_connected() {
        // On a connected graph, most aggregates should have size > 1
        let n = 20;
        let strong = chain_strong(n);
        let agg_ids = greedy_aggregate(&strong, n);

        let n_agg = *agg_ids.iter().max().unwrap_or(&0) + 1;
        let mut sizes = vec![0usize; n_agg];
        for &id in &agg_ids {
            sizes[id] += 1;
        }

        // For a chain of 20, greedy pass-1 creates aggregates of size ≥ 2
        // (at least the majority should have size > 1)
        let n_multi = sizes.iter().filter(|&&s| s > 1).count();
        assert!(
            n_multi > 0,
            "Expected multi-node aggregates, but got sizes: {sizes:?}"
        );

        // Total nodes covered must equal n
        let total: usize = sizes.iter().sum();
        assert_eq!(total, n, "Aggregates don't cover all nodes: {total} != {n}");
    }

    #[test]
    fn test_aggregation_small_graph() {
        // 4-node complete graph: all nodes connected to each other
        let strong = vec![
            vec![1usize, 2, 3],
            vec![0usize, 2, 3],
            vec![0usize, 1, 3],
            vec![0usize, 1, 2],
        ];
        let n = 4;
        let agg_ids = greedy_aggregate(&strong, n);

        // All 4 nodes should be in the same aggregate (node 0 seeds and takes all)
        let first_id = agg_ids[0];
        assert!(
            agg_ids.iter().all(|&id| id == first_id),
            "Expected all nodes in same aggregate, got: {agg_ids:?}"
        );
    }
}
