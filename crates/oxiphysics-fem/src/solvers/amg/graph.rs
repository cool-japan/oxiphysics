// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! C/F splitting and strong-connection graph for classical AMG.

use std::collections::BinaryHeap;

use crate::parallel_solver::CsrMatrix;

/// Compute strong connections for each row.
///
/// For row `i`, a column `j ≠ i` is strongly connected if the (negated)
/// off-diagonal entry `-A[i,j]` exceeds `theta * max_neg_off_diag`, where
/// `max_neg_off_diag` is the maximum of all negated off-diagonal entries in row `i`.
///
/// This is the standard Ruge-Stüben measure for M-matrices.
pub fn strong_connections(a: &CsrMatrix, theta: f64) -> Vec<Vec<usize>> {
    let n = a.nrows;
    let mut strong = vec![Vec::new(); n];

    for (i, strong_row) in strong.iter_mut().enumerate() {
        let rs = a.row_offsets[i];
        let re = a.row_offsets[i + 1];

        // Find maximum negated off-diagonal entry in row i
        let mut max_neg = 0.0f64;
        for k in rs..re {
            let j = a.col_indices[k];
            if j != i {
                let neg_val = -a.values[k];
                if neg_val > max_neg {
                    max_neg = neg_val;
                }
            }
        }

        if max_neg < 1e-300 {
            // No significant off-diagonal entries; no strong connections
            continue;
        }

        let threshold = theta * max_neg;
        for k in rs..re {
            let j = a.col_indices[k];
            if j != i && -a.values[k] >= threshold {
                strong_row.push(j);
            }
        }
    }

    strong
}

/// Ruge-Stüben C/F splitting.
///
/// Returns a boolean vector where `true` means C-point (coarse) and `false`
/// means F-point (fine). The algorithm is:
///
/// **Phase 1**: Uses a priority-queue based greedy selection.
/// - λ_i = |{j : i ∈ strong\[j\]}| + |strong\[i\]|  (measure of influence)
/// - Process nodes in decreasing λ order.
/// - If unmarked: mark as C-point; mark all strongly-connected unmarked neighbors as F-points;
///   increment λ for their transitive strong-connection neighbors.
///
/// **Phase 2**: For each F-point with no C-point among its strong neighbors,
/// promote it to a C-point.
pub fn cf_splitting(strong: &[Vec<usize>]) -> Vec<bool> {
    let n = strong.len();

    // Build transpose: strong_of[j] = set of i where j ∈ strong[i]
    let mut strong_of: Vec<Vec<usize>> = vec![Vec::new(); n];
    for (i, row) in strong.iter().enumerate() {
        for &j in row {
            strong_of[j].push(i);
        }
    }

    // Initial measure λ_i = |strong[i]| + |{j : i ∈ strong[j]}|
    let mut lambda: Vec<usize> = (0..n)
        .map(|i| strong[i].len() + strong_of[i].len())
        .collect();

    // State: 0 = unmarked, 1 = C-point, 2 = F-point
    let mut state: Vec<u8> = vec![0; n];

    // Priority queue: (lambda, node) — use Reverse for max-heap behaviour
    // We store negative lambda so BinaryHeap (max-heap) gives us the largest first
    let mut heap: BinaryHeap<(usize, usize)> = (0..n).map(|i| (lambda[i], i)).collect();

    // Phase 1: greedy C/F assignment
    while let Some((lam_stored, i)) = heap.pop() {
        // Skip stale entries (lambda may have changed since insertion)
        if state[i] != 0 {
            continue;
        }
        if lam_stored != lambda[i] {
            // Re-insert with updated lambda if still unmarked
            heap.push((lambda[i], i));
            continue;
        }

        // Mark i as C-point
        state[i] = 1;

        // Mark all unmarked strong neighbors of i as F-points
        for &j in &strong[i] {
            if state[j] == 0 {
                state[j] = 2;
                // Increment lambda for nodes that have j as a strong connection
                for &k in &strong_of[j] {
                    if state[k] == 0 {
                        lambda[k] += 1;
                        // Re-insert with updated priority
                        heap.push((lambda[k], k));
                    }
                }
            }
        }
    }

    // Any remaining unmarked nodes become C-points
    for s in state.iter_mut() {
        if *s == 0 {
            *s = 1;
        }
    }

    // Phase 2: Ensure every F-point has at least one C-point in its strong connections
    for i in 0..n {
        if state[i] == 2 {
            let has_c_neighbor = strong[i].iter().any(|&j| state[j] == 1);
            if !has_c_neighbor {
                state[i] = 1; // Promote to C-point
            }
        }
    }

    state.iter().map(|&s| s == 1).collect()
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a 1D Poisson tridiagonal matrix of size n.
    fn make_1d_poisson(n: usize) -> CsrMatrix {
        let mut row_offsets = vec![0usize; n + 1];
        let mut col_indices = Vec::new();
        let mut values = Vec::new();

        for i in 0..n {
            if i > 0 {
                col_indices.push(i - 1);
                values.push(-1.0);
            }
            col_indices.push(i);
            values.push(2.0);
            if i + 1 < n {
                col_indices.push(i + 1);
                values.push(-1.0);
            }
            row_offsets[i + 1] = col_indices.len();
        }

        CsrMatrix {
            nrows: n,
            ncols: n,
            row_offsets,
            col_indices,
            values,
        }
    }

    #[test]
    fn test_strong_connection_1d_poisson() {
        let n = 5;
        let a = make_1d_poisson(n);
        let theta = 0.25;
        let strong = strong_connections(&a, theta);

        // For 1D Poisson: max_neg_off_diag = 1.0 (from -1 entries), threshold = 0.25
        // All off-diagonals (-1.0) satisfy -(-1.0) = 1.0 >= 0.25, so all off-diags are strong
        // Interior nodes: 2 strong connections; boundary nodes: 1 strong connection
        for (i, nbrs) in strong.iter().enumerate().take(n) {
            // No self-connections
            assert!(!nbrs.contains(&i), "node {i} has self-connection");
            // Off-diagonal connections should be present
            if i > 0 {
                assert!(
                    nbrs.contains(&(i - 1)),
                    "node {i} missing strong connection to {}",
                    i - 1
                );
            }
            if i + 1 < n {
                assert!(
                    nbrs.contains(&(i + 1)),
                    "node {i} missing strong connection to {}",
                    i + 1
                );
            }
        }
    }

    #[test]
    fn test_cf_splitting_has_c_neighbor() {
        let n = 10;
        let a = make_1d_poisson(n);
        let strong = strong_connections(&a, 0.25);
        let is_c = cf_splitting(&strong);

        // Every F-point must have at least one C-point among its strong neighbors
        for i in 0..n {
            if !is_c[i] {
                let has_c = strong[i].iter().any(|&j| is_c[j]);
                assert!(
                    has_c,
                    "F-point {i} has no C-point among its strong neighbors {:?}",
                    strong[i]
                );
            }
        }

        // At least one C-point must exist
        assert!(is_c.iter().any(|&c| c), "No C-points found");
        // At least one F-point for a non-trivial result (n=10 guarantees this)
        assert!(is_c.iter().any(|&c| !c), "No F-points found");
    }

    #[test]
    fn test_cf_splitting_heap_ordering() {
        // Larger lambda → selected first as C-point; neighbours become F.
        // Star graph: node 0 connected to 1 and 2 (strong), 1 and 2 connected only to 0
        let strong = vec![vec![1usize, 2], vec![0usize], vec![0usize]];
        let is_c = cf_splitting(&strong);

        // Node 0 has highest lambda, should be C-point
        // Nodes 1 and 2 should be F-points (they have C-neighbor 0)
        assert!(is_c[0], "Center node should be C-point");
        // If center is C, then its F-neighbors (1,2) have a C-neighbor => valid
        for i in [1usize, 2] {
            if !is_c[i] {
                let has_c = strong[i].iter().any(|&j| is_c[j]);
                assert!(has_c, "F-point {i} has no C-point");
            }
        }
    }
}
