// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Uniform spatial hash grid for fast neighbourhood queries.
//!
//! A `SpatialGrid` partitions world space into cubic cells of a fixed size.
//! Each registered entry (identified by a `u32` ID) is mapped to the cell that
//! contains its position.  Radius and AABB queries first restrict candidate
//! entries to the relevant cell range before applying the exact geometric test,
//! reducing work compared with a naïve linear scan.
//!
//! ## Choosing `cell_size`
//!
//! Set `cell_size` ≈ 2 × the expected query radius.  Smaller cells improve
//! filtering at the cost of larger cell ranges; larger cells degenerate toward
//! a linear scan.
//!
//! ## Example
//!
//! ```rust,no_run
//! use oxiphysics::spatial_grid::SpatialGrid;
//!
//! let mut grid = SpatialGrid::new(5.0); // 5 m cells
//! grid.insert(0, [1.0, 0.0,  2.0]);
//! grid.insert(1, [3.0, 0.0, -1.0]);
//! grid.insert(2, [100.0, 0.0, 100.0]); // far away
//!
//! let hits = grid.query_radius([0.0, 0.0, 0.0], 5.0);
//! assert!(hits.contains(&0));
//! assert!(hits.contains(&1));
//! assert!(!hits.contains(&2));
//! ```

#![allow(missing_docs)]
#![allow(dead_code)]

use serde::{Deserialize, Serialize};

// ---------------------------------------------------------------------------
// CellCoord — private
// ---------------------------------------------------------------------------

/// Integer cell coordinate in 3D space.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
struct CellCoord {
    x: i32,
    y: i32,
    z: i32,
}

impl CellCoord {
    #[inline]
    fn from_pos(pos: [f64; 3], inv_cell: f64) -> Self {
        CellCoord {
            x: (pos[0] * inv_cell).floor() as i32,
            y: (pos[1] * inv_cell).floor() as i32,
            z: (pos[2] * inv_cell).floor() as i32,
        }
    }

    #[inline]
    fn in_range(self, min: CellCoord, max: CellCoord) -> bool {
        self.x >= min.x
            && self.x <= max.x
            && self.y >= min.y
            && self.y <= max.y
            && self.z >= min.z
            && self.z <= max.z
    }
}

// ---------------------------------------------------------------------------
// GridEntry — private
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Serialize, Deserialize)]
struct GridEntry {
    id: u32,
    position: [f64; 3],
    cell: CellCoord,
}

// ---------------------------------------------------------------------------
// SpatialGrid
// ---------------------------------------------------------------------------

/// A uniform spatial hash grid for fast neighbourhood queries.
///
/// Entries are keyed by `u32` IDs. Each ID may appear at most once; calling
/// [`insert`][Self::insert] with an existing ID replaces the entry.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SpatialGrid {
    /// Side length of each cubic cell (m).
    pub cell_size: f64,
    entries: Vec<GridEntry>,
}

impl SpatialGrid {
    /// Create an empty grid with the given cell size.
    ///
    /// Panics if `cell_size <= 0.0`.
    pub fn new(cell_size: f64) -> Self {
        assert!(cell_size > 0.0, "cell_size must be positive");
        Self {
            cell_size,
            entries: Vec::new(),
        }
    }

    // -----------------------------------------------------------------------
    // Mutation
    // -----------------------------------------------------------------------

    fn inv_cell(&self) -> f64 {
        1.0 / self.cell_size
    }

    /// Insert or replace an entry.
    ///
    /// If `id` already exists, its position is updated in place.
    pub fn insert(&mut self, id: u32, position: [f64; 3]) {
        let cell = CellCoord::from_pos(position, self.inv_cell());
        if let Some(e) = self.entries.iter_mut().find(|e| e.id == id) {
            e.position = position;
            e.cell = cell;
        } else {
            self.entries.push(GridEntry { id, position, cell });
        }
    }

    /// Update the position of an existing entry.  No-op if `id` is not found.
    pub fn update(&mut self, id: u32, position: [f64; 3]) {
        let cell = CellCoord::from_pos(position, self.inv_cell());
        if let Some(e) = self.entries.iter_mut().find(|e| e.id == id) {
            e.position = position;
            e.cell = cell;
        }
    }

    /// Remove an entry by ID.  No-op if `id` is not found.
    pub fn remove(&mut self, id: u32) {
        self.entries.retain(|e| e.id != id);
    }

    /// Remove all entries.
    pub fn clear(&mut self) {
        self.entries.clear();
    }

    // -----------------------------------------------------------------------
    // Single-entry queries
    // -----------------------------------------------------------------------

    /// Return `true` if `id` is registered in the grid.
    pub fn contains(&self, id: u32) -> bool {
        self.entries.iter().any(|e| e.id == id)
    }

    /// Return the current position of `id`, or `None` if not registered.
    pub fn position(&self, id: u32) -> Option<[f64; 3]> {
        self.entries.iter().find(|e| e.id == id).map(|e| e.position)
    }

    // -----------------------------------------------------------------------
    // Neighbourhood queries
    // -----------------------------------------------------------------------

    /// Return all IDs whose position is within `radius` of `center`.
    pub fn query_radius(&self, center: [f64; 3], radius: f64) -> Vec<u32> {
        let inv = self.inv_cell();
        let min_cell = CellCoord {
            x: ((center[0] - radius) * inv).floor() as i32,
            y: ((center[1] - radius) * inv).floor() as i32,
            z: ((center[2] - radius) * inv).floor() as i32,
        };
        let max_cell = CellCoord {
            x: ((center[0] + radius) * inv).floor() as i32,
            y: ((center[1] + radius) * inv).floor() as i32,
            z: ((center[2] + radius) * inv).floor() as i32,
        };
        let r2 = radius * radius;
        self.entries
            .iter()
            .filter(|e| {
                e.cell.in_range(min_cell, max_cell) && {
                    let dx = e.position[0] - center[0];
                    let dy = e.position[1] - center[1];
                    let dz = e.position[2] - center[2];
                    dx * dx + dy * dy + dz * dz <= r2
                }
            })
            .map(|e| e.id)
            .collect()
    }

    /// Return all IDs whose position lies inside the AABB `[min, max]`.
    pub fn query_aabb(&self, min: [f64; 3], max: [f64; 3]) -> Vec<u32> {
        let inv = self.inv_cell();
        let min_cell = CellCoord {
            x: (min[0] * inv).floor() as i32,
            y: (min[1] * inv).floor() as i32,
            z: (min[2] * inv).floor() as i32,
        };
        let max_cell = CellCoord {
            x: (max[0] * inv).floor() as i32,
            y: (max[1] * inv).floor() as i32,
            z: (max[2] * inv).floor() as i32,
        };
        self.entries
            .iter()
            .filter(|e| {
                e.cell.in_range(min_cell, max_cell)
                    && e.position[0] >= min[0]
                    && e.position[0] <= max[0]
                    && e.position[1] >= min[1]
                    && e.position[1] <= max[1]
                    && e.position[2] >= min[2]
                    && e.position[2] <= max[2]
            })
            .map(|e| e.id)
            .collect()
    }

    /// Return the ID and distance of the nearest entry to `pos`.
    ///
    /// Returns `None` if the grid is empty.
    pub fn nearest(&self, pos: [f64; 3]) -> Option<(u32, f64)> {
        self.entries
            .iter()
            .map(|e| {
                let dx = e.position[0] - pos[0];
                let dy = e.position[1] - pos[1];
                let dz = e.position[2] - pos[2];
                (e.id, (dx * dx + dy * dy + dz * dz).sqrt())
            })
            .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
    }

    /// Return the `k` nearest entries to `pos`, sorted by ascending distance.
    ///
    /// Returns fewer than `k` entries if the grid has fewer.
    pub fn k_nearest(&self, pos: [f64; 3], k: usize) -> Vec<(u32, f64)> {
        let mut dists: Vec<(u32, f64)> = self
            .entries
            .iter()
            .map(|e| {
                let dx = e.position[0] - pos[0];
                let dy = e.position[1] - pos[1];
                let dz = e.position[2] - pos[2];
                (e.id, (dx * dx + dy * dy + dz * dz).sqrt())
            })
            .collect();
        dists.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        dists.truncate(k);
        dists
    }

    /// Return all pairs `(a, b)` whose positions are within `radius` of each
    /// other.
    ///
    /// Each pair is returned once with `a < b`.  Useful for broad-phase
    /// collision detection.
    pub fn pairs_within_radius(&self, radius: f64) -> Vec<(u32, u32)> {
        let r2 = radius * radius;
        let n = self.entries.len();
        let mut pairs = Vec::new();
        for i in 0..n {
            let a = &self.entries[i];
            for j in (i + 1)..n {
                let b = &self.entries[j];
                let dx = a.position[0] - b.position[0];
                let dy = a.position[1] - b.position[1];
                let dz = a.position[2] - b.position[2];
                if dx * dx + dy * dy + dz * dz <= r2 {
                    let (lo, hi) = if a.id < b.id {
                        (a.id, b.id)
                    } else {
                        (b.id, a.id)
                    };
                    pairs.push((lo, hi));
                }
            }
        }
        pairs
    }

    // -----------------------------------------------------------------------
    // Metadata
    // -----------------------------------------------------------------------

    /// Number of entries in the grid.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Return `true` if the grid contains no entries.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Iterate over `(id, position)` pairs.
    pub fn iter(&self) -> impl Iterator<Item = (u32, [f64; 3])> + '_ {
        self.entries.iter().map(|e| (e.id, e.position))
    }

    /// Return IDs of all entries whose cell coordinate matches exactly.
    ///
    /// Useful for debugging grid partitioning.
    pub fn entries_in_cell(&self, cell_x: i32, cell_y: i32, cell_z: i32) -> Vec<u32> {
        let target = CellCoord {
            x: cell_x,
            y: cell_y,
            z: cell_z,
        };
        self.entries
            .iter()
            .filter(|e| e.cell == target)
            .map(|e| e.id)
            .collect()
    }

    /// Return the cell coordinates for a world-space position.
    pub fn cell_for(&self, pos: [f64; 3]) -> (i32, i32, i32) {
        let c = CellCoord::from_pos(pos, self.inv_cell());
        (c.x, c.y, c.z)
    }
}
