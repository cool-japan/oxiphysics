// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! A scalar function on the vertices of a mesh graph, with the total order of
//! its vertices shared by the Reeb graph and the Morse complex.

/// A positive quiet NaN.  [`f64::total_cmp`] orders it above every number,
/// whereas a NaN with the sign bit set would come first.
const POSITIVE_NAN: f64 = f64::from_bits(0x7ff8_0000_0000_0000);

/// Vertex order and adjacency of a scalar function on a mesh graph.
///
/// Vertices are ordered by function value, ties broken by vertex index (the
/// "simulation of simplicity" total order), so every pair of vertices has a
/// definite lower one and plateaus of equal values are resolved consistently.
/// Values are compared with [`f64::total_cmp`] after mapping `-0.0` to `0.0`
/// and every NaN to a positive NaN, which gives NaN a fixed place (above every
/// number, all NaNs tied) instead of breaking the order.
pub(super) struct ScalarGraph {
    /// Vertices in ascending order.
    pub(super) order: Vec<usize>,
    /// `rank[v]` is the position of vertex `v` in `order`.
    pub(super) rank: Vec<usize>,
    /// Sorted neighbours of every vertex.  Edges with an index outside the
    /// function, self loops and repeated edges are dropped.
    pub(super) neighbors: Vec<Vec<usize>>,
}

impl ScalarGraph {
    /// Order the vertices of `values` and build the adjacency of `edges`.
    pub(super) fn new(values: &[f64], edges: &[[usize; 2]]) -> Self {
        let n = values.len();
        // `+ 0.0` turns -0.0 into 0.0 and leaves every other number unchanged.
        let key = |v: usize| {
            let value = values[v];
            if value.is_nan() {
                POSITIVE_NAN
            } else {
                value + 0.0
            }
        };
        let mut order: Vec<usize> = (0..n).collect();
        order.sort_by(|&a, &b| key(a).total_cmp(&key(b)).then(a.cmp(&b)));
        let mut rank = vec![0usize; n];
        for (position, &v) in order.iter().enumerate() {
            rank[v] = position;
        }
        let mut neighbors: Vec<Vec<usize>> = vec![Vec::new(); n];
        for &[i, j] in edges {
            if i < n && j < n && i != j {
                neighbors[i].push(j);
                neighbors[j].push(i);
            }
        }
        for list in &mut neighbors {
            list.sort_unstable();
            list.dedup();
        }
        Self {
            order,
            rank,
            neighbors,
        }
    }

    /// Number of vertices.
    pub(super) fn len(&self) -> usize {
        self.order.len()
    }

    /// Whether `u` comes before `v` in the vertex order.
    pub(super) fn is_below(&self, u: usize, v: usize) -> bool {
        self.rank[u] < self.rank[v]
    }
}
