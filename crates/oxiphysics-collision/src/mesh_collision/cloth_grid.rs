// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Spatial hash and adjacency tables for [`super::cloth_self_collision`].
//!
//! [`VertexGrid`] buckets vertex indices into cubic cells whose edge is at least
//! the collision thickness. Two points closer than the thickness therefore lie
//! in the same or in adjacent cells, so the 27 cells around a point contain
//! every vertex that can collide with it. The grid can relocate a vertex, which
//! lets the self-collision pass keep it in sync with the position corrections it
//! applies while it runs.

use std::collections::HashMap;

/// Largest magnitude of a cell coordinate.
///
/// Coordinates beyond it are clamped so the neighbour arithmetic cannot
/// overflow. The exact distance test still decides every candidate pair, so
/// clamping only costs performance, never correctness.
const CELL_LIMIT: f64 = 1.0e15;

/// Relative slack added to the cell edge.
///
/// `floor(x * inv_cell)` is computed with rounding error, so two points that are
/// closer than the nominal edge could otherwise end up two cells apart when they
/// are extremely far from the origin. The slack keeps them adjacent for all
/// cell indices below about `4e9`.
const CELL_SLACK: f64 = 1.0 + 1.0e-6;

/// Uniform grid mapping cells to the indices of the vertices inside them.
pub(super) struct VertexGrid {
    inv_cell: f64,
    cells: HashMap<[i64; 3], Vec<usize>>,
}

impl VertexGrid {
    /// Create an empty grid whose cells are at least `min_cell_edge` wide.
    ///
    /// `min_cell_edge` must be positive and finite.
    pub(super) fn new(min_cell_edge: f64) -> Self {
        Self {
            inv_cell: 1.0 / (min_cell_edge * CELL_SLACK),
            cells: HashMap::new(),
        }
    }

    /// Cell containing `p`. A NaN coordinate maps to cell zero.
    pub(super) fn cell_of(&self, p: [f64; 3]) -> [i64; 3] {
        let coordinate =
            |x: f64| -> i64 { (x * self.inv_cell).floor().clamp(-CELL_LIMIT, CELL_LIMIT) as i64 };
        [coordinate(p[0]), coordinate(p[1]), coordinate(p[2])]
    }

    /// Register vertex `index` at position `p`.
    pub(super) fn insert(&mut self, index: usize, p: [f64; 3]) {
        let key = self.cell_of(p);
        self.cells.entry(key).or_default().push(index);
    }

    /// Move vertex `index` from position `from` to position `to`.
    ///
    /// `from` must be the position the vertex was last registered at.
    pub(super) fn relocate(&mut self, index: usize, from: [f64; 3], to: [f64; 3]) {
        let old_key = self.cell_of(from);
        let new_key = self.cell_of(to);
        if old_key == new_key {
            return;
        }
        let emptied = match self.cells.get_mut(&old_key) {
            Some(bucket) => {
                if let Some(slot) = bucket.iter().position(|&i| i == index) {
                    bucket.swap_remove(slot);
                }
                bucket.is_empty()
            }
            None => false,
        };
        if emptied {
            self.cells.remove(&old_key);
        }
        self.cells.entry(new_key).or_default().push(index);
    }

    /// Collect, in ascending order, the vertices with an index above `after`
    /// that lie in the 27 cells around `p` and are not listed in `excluded`
    /// (a sorted slice).
    ///
    /// The result contains every such vertex closer to `p` than the cell edge
    /// passed to [`Self::new`], and possibly more. `out` is cleared first.
    pub(super) fn partners_after(
        &self,
        p: [f64; 3],
        after: usize,
        excluded: &[usize],
        out: &mut Vec<usize>,
    ) {
        out.clear();
        self.for_each_near(p, |j| {
            if j > after && excluded.binary_search(&j).is_err() {
                out.push(j);
            }
        });
        out.sort_unstable();
    }

    /// Call `visit` with the index of every vertex in the 27 cells around `p`.
    fn for_each_near<F: FnMut(usize)>(&self, p: [f64; 3], mut visit: F) {
        let [cx, cy, cz] = self.cell_of(p);
        for dx in -1..=1_i64 {
            for dy in -1..=1_i64 {
                for dz in -1..=1_i64 {
                    let key = [
                        cx.saturating_add(dx),
                        cy.saturating_add(dy),
                        cz.saturating_add(dz),
                    ];
                    if let Some(bucket) = self.cells.get(&key) {
                        for &index in bucket {
                            visit(index);
                        }
                    }
                }
            }
        }
    }
}

/// Sorted, de-duplicated neighbour lists: `adjacency[i]` holds every vertex that
/// shares a triangle with vertex `i`.
///
/// Triangle corners that are not valid vertex indices (`>= vertex_count`) are
/// ignored as owners but may still appear as neighbours; they are never queried.
pub(super) fn vertex_adjacency(vertex_count: usize, triangles: &[[usize; 3]]) -> Vec<Vec<usize>> {
    let mut adjacency: Vec<Vec<usize>> = vec![Vec::new(); vertex_count];
    for triangle in triangles {
        for (corner, &vertex) in triangle.iter().enumerate() {
            let Some(list) = adjacency.get_mut(vertex) else {
                continue;
            };
            for (other_corner, &other) in triangle.iter().enumerate() {
                if other_corner != corner {
                    list.push(other);
                }
            }
        }
    }
    for list in &mut adjacency {
        list.sort_unstable();
        list.dedup();
    }
    adjacency
}
