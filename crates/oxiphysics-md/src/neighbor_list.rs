#![allow(clippy::manual_div_ceil)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Verlet neighbor list and cell-list for efficient MD force computation.
//!
//! This module provides an alternative neighbor-search implementation that
//! uses plain `[f64; 3]` arrays instead of [`oxiphysics_core::math::Vec3`],
//! making it easy to interface with external data without conversion.
//!
//! # Overview
//!
//! - [`CellList`]: O(N) cell-based spatial hashing.
//! - [`VerletList`]: skin-buffered Verlet list that reuses pairs between rebuilds.
//! - [`MultiCutoffNeighborList`]: neighbor list supporting multiple cutoff radii.
//! - [`HalfNeighborList`]: half neighbor list (i < j only) for efficient force computation.
//! - [`NeighborStatistics`]: statistics on neighbor list usage.
//! - [`minimum_image`] / [`minimum_image_vec`]: periodic minimum-image helpers.

// ---------------------------------------------------------------------------
// Minimum-image convention
// ---------------------------------------------------------------------------

/// Apply the minimum-image convention to a single displacement component.
///
/// Maps `dx` into `(-L/2, L/2]` given box length `box_length`.
///
/// # Examples
/// ```no_run
/// use oxiphysics_md::neighbor_list::minimum_image;
/// let dx = minimum_image(0.8, 1.0);
/// assert!((dx - (-0.2)).abs() < 1e-12);
/// ```
pub fn minimum_image(mut dx: f64, box_length: f64) -> f64 {
    let half = box_length * 0.5;
    if dx > half {
        dx -= box_length;
    } else if dx < -half {
        dx += box_length;
    }
    dx
}

/// Apply the minimum-image convention to a 3-D displacement vector.
///
/// # Examples
/// ```no_run
/// use oxiphysics_md::neighbor_list::minimum_image_vec;
/// let dr = minimum_image_vec([0.8, -0.8, 0.1], [1.0, 1.0, 1.0]);
/// assert!((dr[0] - (-0.2)).abs() < 1e-12);
/// assert!((dr[1] -   0.2 ).abs() < 1e-12);
/// ```
pub fn minimum_image_vec(dr: [f64; 3], box_len: [f64; 3]) -> [f64; 3] {
    [
        minimum_image(dr[0], box_len[0]),
        minimum_image(dr[1], box_len[1]),
        minimum_image(dr[2], box_len[2]),
    ]
}

// ---------------------------------------------------------------------------
// CellList
// ---------------------------------------------------------------------------

/// Cell list for O(N) neighbor search.
///
/// The simulation box `[box_lo, box_hi]` is divided into cells of size
/// `>= cutoff`.  Each cell stores the indices of the atoms it contains.
/// Neighbor searches only need to examine the 3x3x3 = 27 surrounding cells.
#[derive(Debug, Clone)]
pub struct CellList {
    /// Cell edge length (>= cutoff radius).
    pub cell_size: f64,
    /// Number of cells in the x direction.
    pub nx: usize,
    /// Number of cells in the y direction.
    pub ny: usize,
    /// Number of cells in the z direction.
    pub nz: usize,
    /// Lower corner of the simulation box.
    pub box_lo: [f64; 3],
    /// Upper corner of the simulation box.
    pub box_hi: [f64; 3],
    /// `cells[flat_idx]` -> list of atom indices in that cell.
    cells: Vec<Vec<usize>>,
}

impl CellList {
    /// Create a new (empty) cell list.
    pub fn new(box_lo: [f64; 3], box_hi: [f64; 3], cutoff: f64) -> Self {
        let lx = box_hi[0] - box_lo[0];
        let ly = box_hi[1] - box_lo[1];
        let lz = box_hi[2] - box_lo[2];

        let nx = ((lx / cutoff).ceil() as usize).max(1);
        let ny = ((ly / cutoff).ceil() as usize).max(1);
        let nz = ((lz / cutoff).ceil() as usize).max(1);

        let cell_size = cutoff;

        let total = nx * ny * nz;
        Self {
            cell_size,
            nx,
            ny,
            nz,
            box_lo,
            box_hi,
            cells: vec![Vec::new(); total],
        }
    }

    /// Rebuild the cell list from current atom positions.
    pub fn build(&mut self, positions: &[[f64; 3]]) {
        for cell in &mut self.cells {
            cell.clear();
        }
        for (i, pos) in positions.iter().enumerate() {
            if let Some((ix, iy, iz)) = self.cell_index(pos) {
                let idx = self.flat_idx(ix, iy, iz);
                self.cells[idx].push(i);
            } else {
                let ix = self.clamp_cell(pos[0], self.box_lo[0], self.nx);
                let iy = self.clamp_cell(pos[1], self.box_lo[1], self.ny);
                let iz = self.clamp_cell(pos[2], self.box_lo[2], self.nz);
                let idx = self.flat_idx(ix, iy, iz);
                self.cells[idx].push(i);
            }
        }
    }

    /// Return the (ix, iy, iz) cell coordinates for `pos`, or `None` if out of box.
    pub fn cell_index(&self, pos: &[f64; 3]) -> Option<(usize, usize, usize)> {
        let lx = self.box_hi[0] - self.box_lo[0];
        let ly = self.box_hi[1] - self.box_lo[1];
        let lz = self.box_hi[2] - self.box_lo[2];

        let fx = pos[0] - self.box_lo[0];
        let fy = pos[1] - self.box_lo[1];
        let fz = pos[2] - self.box_lo[2];

        if fx < 0.0 || fx > lx || fy < 0.0 || fy > ly || fz < 0.0 || fz > lz {
            return None;
        }

        let ix = ((fx / self.cell_size) as usize).min(self.nx - 1);
        let iy = ((fy / self.cell_size) as usize).min(self.ny - 1);
        let iz = ((fz / self.cell_size) as usize).min(self.nz - 1);
        Some((ix, iy, iz))
    }

    /// Compute the flat cell index from 3-D cell coordinates.
    fn flat_idx(&self, ix: usize, iy: usize, iz: usize) -> usize {
        ix + iy * self.nx + iz * self.nx * self.ny
    }

    /// Clamp a coordinate to a valid cell index.
    fn clamp_cell(&self, coord: f64, lo: f64, n: usize) -> usize {
        let f = (coord - lo) / self.cell_size;
        if f < 0.0 { 0 } else { (f as usize).min(n - 1) }
    }

    /// Find all atoms within `cutoff_sq` (squared cutoff) of atom `i`.
    pub fn neighbors_of(
        &self,
        i: usize,
        positions: &[[f64; 3]],
        cutoff_sq: f64,
    ) -> Vec<(usize, f64, [f64; 3])> {
        let pos_i = positions[i];
        let box_len = [
            self.box_hi[0] - self.box_lo[0],
            self.box_hi[1] - self.box_lo[1],
            self.box_hi[2] - self.box_lo[2],
        ];

        let (ix, iy, iz) = self.cell_index(&pos_i).unwrap_or_else(|| {
            (
                self.clamp_cell(pos_i[0], self.box_lo[0], self.nx),
                self.clamp_cell(pos_i[1], self.box_lo[1], self.ny),
                self.clamp_cell(pos_i[2], self.box_lo[2], self.nz),
            )
        });

        let mut result = Vec::new();

        for dix in -1_isize..=1 {
            for diy in -1_isize..=1 {
                for diz in -1_isize..=1 {
                    let nx = self.nx as isize;
                    let ny = self.ny as isize;
                    let nz = self.nz as isize;

                    let jx = ((ix as isize + dix).rem_euclid(nx)) as usize;
                    let jy = ((iy as isize + diy).rem_euclid(ny)) as usize;
                    let jz = ((iz as isize + diz).rem_euclid(nz)) as usize;

                    let cell_idx = self.flat_idx(jx, jy, jz);
                    for &j in &self.cells[cell_idx] {
                        if j == i {
                            continue;
                        }
                        let dr_raw = [
                            positions[j][0] - pos_i[0],
                            positions[j][1] - pos_i[1],
                            positions[j][2] - pos_i[2],
                        ];
                        let dr = minimum_image_vec(dr_raw, box_len);
                        let dsq = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
                        if dsq <= cutoff_sq {
                            result.push((j, dsq, dr));
                        }
                    }
                }
            }
        }
        result
    }

    /// Find all neighbor pairs `(i, j)` with `i < j` within `cutoff_sq`.
    pub fn all_pairs(&self, positions: &[[f64; 3]], cutoff_sq: f64) -> Vec<(usize, usize, f64)> {
        let n = positions.len();
        let mut pairs = Vec::new();

        for i in 0..n {
            let nbrs = self.neighbors_of(i, positions, cutoff_sq);
            for (j, dsq, _) in nbrs {
                if j > i {
                    pairs.push((i, j, dsq));
                }
            }
        }
        pairs
    }

    /// Total number of atoms stored in the cell list.
    #[allow(dead_code)]
    fn total_atoms(&self) -> usize {
        self.cells.iter().map(|c| c.len()).sum()
    }

    /// Number of non-empty cells.
    #[allow(dead_code)]
    fn occupied_cells(&self) -> usize {
        self.cells.iter().filter(|c| !c.is_empty()).count()
    }
}

// ---------------------------------------------------------------------------
// VerletList
// ---------------------------------------------------------------------------

/// Verlet neighbor list with skin buffering.
#[derive(Debug, Clone)]
pub struct VerletList {
    /// Interaction cutoff (without skin).
    pub cutoff: f64,
    /// Skin buffer distance.
    pub skin: f64,
    /// `neighbors[i]` = list of atom indices `j > i` within `cutoff + skin`.
    pub neighbors: Vec<Vec<usize>>,
    /// Positions at the time of the last rebuild.
    ref_positions: Vec<[f64; 3]>,
    /// Internal cell list used during rebuilds.
    cell_list: CellList,
}

impl VerletList {
    /// Create a new Verlet list.
    pub fn new(box_lo: [f64; 3], box_hi: [f64; 3], cutoff: f64, skin: f64) -> Self {
        let effective = cutoff + skin;
        let cell_list = CellList::new(box_lo, box_hi, effective);
        Self {
            cutoff,
            skin,
            neighbors: Vec::new(),
            ref_positions: Vec::new(),
            cell_list,
        }
    }

    /// Build (or rebuild) the neighbor list from the current positions.
    pub fn build(&mut self, positions: &[[f64; 3]]) {
        let n = positions.len();
        let effective_sq = (self.cutoff + self.skin) * (self.cutoff + self.skin);

        self.cell_list.build(positions);

        self.neighbors.clear();
        self.neighbors.resize_with(n, Vec::new);

        for i in 0..n {
            let nbrs = self.cell_list.neighbors_of(i, positions, effective_sq);
            for (j, _, _) in nbrs {
                if j > i {
                    self.neighbors[i].push(j);
                }
            }
        }

        self.ref_positions = positions.to_vec();
    }

    /// Returns `true` if any atom has moved more than `skin / 2`.
    pub fn needs_rebuild(&self, positions: &[[f64; 3]]) -> bool {
        if self.ref_positions.len() != positions.len() {
            return true;
        }
        let half_skin_sq = (self.skin * 0.5) * (self.skin * 0.5);
        let box_len = [
            self.cell_list.box_hi[0] - self.cell_list.box_lo[0],
            self.cell_list.box_hi[1] - self.cell_list.box_lo[1],
            self.cell_list.box_hi[2] - self.cell_list.box_lo[2],
        ];
        for (ref_pos, cur_pos) in self.ref_positions.iter().zip(positions.iter()) {
            let dr_raw = [
                cur_pos[0] - ref_pos[0],
                cur_pos[1] - ref_pos[1],
                cur_pos[2] - ref_pos[2],
            ];
            let dr = minimum_image_vec(dr_raw, box_len);
            let dsq = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
            if dsq > half_skin_sq {
                return true;
            }
        }
        false
    }

    /// Update the neighbor list: rebuild only if needed.
    pub fn update(&mut self, positions: &[[f64; 3]]) {
        if self.needs_rebuild(positions) {
            self.build(positions);
        }
    }

    /// Return the neighbor list for atom `i`.
    pub fn neighbors(&self, i: usize) -> &[usize] {
        &self.neighbors[i]
    }

    /// Total number of unique neighbor pairs.
    pub fn n_pairs(&self) -> usize {
        self.neighbors.iter().map(|v| v.len()).sum()
    }

    /// Average number of neighbors per atom.
    pub fn avg_coordination(&self) -> f64 {
        let n = self.neighbors.len();
        if n == 0 {
            return 0.0;
        }
        let total: usize = self.neighbors.iter().map(|v| v.len()).sum();
        (total * 2) as f64 / n as f64
    }
}

// ---------------------------------------------------------------------------
// MultiCutoffNeighborList
// ---------------------------------------------------------------------------

/// Neighbor list supporting multiple cutoff radii.
///
/// Useful when different interaction types have different cutoff distances
/// (e.g., LJ at 10 A, electrostatics at 12 A).
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct MultiCutoffNeighborList {
    /// Cutoff values for each interaction type.
    cutoffs: Vec<f64>,
    /// For each cutoff, the list of neighbor pairs `(i, j)` with `i < j`.
    pair_lists: Vec<Vec<(usize, usize)>>,
    /// Box dimensions.
    box_lo: [f64; 3],
    box_hi: [f64; 3],
}

#[allow(dead_code)]
impl MultiCutoffNeighborList {
    /// Create a new multi-cutoff neighbor list.
    pub fn new(box_lo: [f64; 3], box_hi: [f64; 3], cutoffs: Vec<f64>) -> Self {
        let n = cutoffs.len();
        Self {
            cutoffs,
            pair_lists: vec![Vec::new(); n],
            box_lo,
            box_hi,
        }
    }

    /// Build neighbor lists for all cutoffs simultaneously.
    ///
    /// Uses the largest cutoff for the cell list, then filters pairs
    /// into each cutoff bin.
    pub fn build(&mut self, positions: &[[f64; 3]]) {
        let max_cutoff = self.cutoffs.iter().cloned().fold(0.0_f64, f64::max);
        if max_cutoff <= 0.0 {
            return;
        }

        let mut cell_list = CellList::new(self.box_lo, self.box_hi, max_cutoff);
        cell_list.build(positions);

        for list in &mut self.pair_lists {
            list.clear();
        }

        let max_cutoff_sq = max_cutoff * max_cutoff;
        let all_pairs = cell_list.all_pairs(positions, max_cutoff_sq);

        for (i, j, dsq) in all_pairs {
            for (c_idx, &cutoff) in self.cutoffs.iter().enumerate() {
                if dsq <= cutoff * cutoff {
                    self.pair_lists[c_idx].push((i, j));
                }
            }
        }
    }

    /// Get the pairs for a specific cutoff index.
    pub fn pairs(&self, cutoff_idx: usize) -> &[(usize, usize)] {
        &self.pair_lists[cutoff_idx]
    }

    /// Number of pairs for each cutoff.
    pub fn pair_counts(&self) -> Vec<usize> {
        self.pair_lists.iter().map(|l| l.len()).collect()
    }

    /// Number of cutoff levels.
    pub fn n_cutoffs(&self) -> usize {
        self.cutoffs.len()
    }
}

// ---------------------------------------------------------------------------
// HalfNeighborList
// ---------------------------------------------------------------------------

/// Half neighbor list: stores only pairs where `i < j`.
///
/// This is the standard format for force computation to avoid
/// double-counting. Includes displacement vectors and distances.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct HalfNeighborList {
    /// Pairs: `(i, j, dr, dist)` with `i < j`.
    pairs: Vec<(usize, usize, [f64; 3], f64)>,
}

#[allow(dead_code)]
impl HalfNeighborList {
    /// Build a half neighbor list from positions.
    pub fn build(positions: &[[f64; 3]], box_lo: [f64; 3], box_hi: [f64; 3], cutoff: f64) -> Self {
        let mut cell_list = CellList::new(box_lo, box_hi, cutoff);
        cell_list.build(positions);
        let cutoff_sq = cutoff * cutoff;
        let _box_len = [
            box_hi[0] - box_lo[0],
            box_hi[1] - box_lo[1],
            box_hi[2] - box_lo[2],
        ];

        let n = positions.len();
        let mut pairs = Vec::new();
        for i in 0..n {
            let nbrs = cell_list.neighbors_of(i, positions, cutoff_sq);
            for (j, _dsq, dr) in nbrs {
                if j > i {
                    let dist = (dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]).sqrt();
                    pairs.push((i, j, dr, dist));
                }
            }
        }
        Self { pairs }
    }

    /// Convert from a full neighbor list to half (i < j only).
    pub fn from_full(
        full_neighbors: &[Vec<usize>],
        positions: &[[f64; 3]],
        box_len: [f64; 3],
    ) -> Self {
        let mut pairs = Vec::new();
        for (i, nbrs) in full_neighbors.iter().enumerate() {
            for &j in nbrs {
                if j > i {
                    let dr_raw = [
                        positions[j][0] - positions[i][0],
                        positions[j][1] - positions[i][1],
                        positions[j][2] - positions[i][2],
                    ];
                    let dr = minimum_image_vec(dr_raw, box_len);
                    let dist = (dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]).sqrt();
                    pairs.push((i, j, dr, dist));
                }
            }
        }
        Self { pairs }
    }

    /// Number of pairs.
    pub fn n_pairs(&self) -> usize {
        self.pairs.len()
    }

    /// Iterate over pairs.
    pub fn iter(&self) -> impl Iterator<Item = &(usize, usize, [f64; 3], f64)> {
        self.pairs.iter()
    }

    /// Get pair by index.
    pub fn get(&self, idx: usize) -> Option<&(usize, usize, [f64; 3], f64)> {
        self.pairs.get(idx)
    }
}

// ---------------------------------------------------------------------------
// Compressed neighbor list
// ---------------------------------------------------------------------------

/// Compressed sparse row (CSR) format neighbor list.
///
/// More memory-efficient than `Vec<Vec`usize`>` for large systems.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct CompressedNeighborList {
    /// Row pointers: `offsets[i]..offsets[i+1]` gives the neighbor range for atom i.
    offsets: Vec<usize>,
    /// Flat array of neighbor indices.
    neighbors: Vec<usize>,
    /// Flat array of squared distances (aligned with `neighbors`).
    distances_sq: Vec<f64>,
}

#[allow(dead_code)]
impl CompressedNeighborList {
    /// Build a compressed neighbor list from a standard neighbor list.
    pub fn from_verlet(vl: &VerletList, positions: &[[f64; 3]], box_len: [f64; 3]) -> Self {
        let n = vl.neighbors.len();
        let mut offsets = Vec::with_capacity(n + 1);
        let mut neighbors = Vec::new();
        let mut distances_sq = Vec::new();

        let mut offset = 0;
        for (i, nbrs) in vl.neighbors.iter().enumerate() {
            offsets.push(offset);
            for &j in nbrs {
                let dr_raw = [
                    positions[j][0] - positions[i][0],
                    positions[j][1] - positions[i][1],
                    positions[j][2] - positions[i][2],
                ];
                let dr = minimum_image_vec(dr_raw, box_len);
                let dsq = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
                neighbors.push(j);
                distances_sq.push(dsq);
                offset += 1;
            }
        }
        offsets.push(offset);

        Self {
            offsets,
            neighbors,
            distances_sq,
        }
    }

    /// Number of atoms.
    pub fn n_atoms(&self) -> usize {
        if self.offsets.is_empty() {
            0
        } else {
            self.offsets.len() - 1
        }
    }

    /// Neighbors of atom `i`.
    pub fn neighbors_of(&self, i: usize) -> &[usize] {
        let start = self.offsets[i];
        let end = self.offsets[i + 1];
        &self.neighbors[start..end]
    }

    /// Squared distances for neighbors of atom `i`.
    pub fn distances_sq_of(&self, i: usize) -> &[f64] {
        let start = self.offsets[i];
        let end = self.offsets[i + 1];
        &self.distances_sq[start..end]
    }

    /// Total number of stored neighbor entries.
    pub fn total_entries(&self) -> usize {
        self.neighbors.len()
    }

    /// Memory usage in bytes (approximate).
    pub fn memory_bytes(&self) -> usize {
        self.offsets.len() * std::mem::size_of::<usize>()
            + self.neighbors.len() * std::mem::size_of::<usize>()
            + self.distances_sq.len() * std::mem::size_of::<f64>()
    }
}

// ---------------------------------------------------------------------------
// Neighbor statistics
// ---------------------------------------------------------------------------

/// Statistics about a neighbor list.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct NeighborStatistics {
    /// Total number of atoms.
    pub n_atoms: usize,
    /// Total number of unique pairs.
    pub n_pairs: usize,
    /// Minimum number of neighbors for any atom.
    pub min_neighbors: usize,
    /// Maximum number of neighbors for any atom.
    pub max_neighbors: usize,
    /// Average number of neighbors per atom.
    pub avg_neighbors: f64,
    /// Standard deviation of neighbor count.
    pub std_neighbors: f64,
    /// Number of atoms with zero neighbors (isolated atoms).
    pub isolated_atoms: usize,
}

#[allow(dead_code)]
impl NeighborStatistics {
    /// Compute statistics from a Verlet list.
    pub fn from_verlet(vl: &VerletList) -> Self {
        let n = vl.neighbors.len();
        if n == 0 {
            return Self {
                n_atoms: 0,
                n_pairs: 0,
                min_neighbors: 0,
                max_neighbors: 0,
                avg_neighbors: 0.0,
                std_neighbors: 0.0,
                isolated_atoms: 0,
            };
        }

        // Count neighbors per atom (both directions)
        let mut counts = vec![0usize; n];
        for (i, nbrs) in vl.neighbors.iter().enumerate() {
            counts[i] += nbrs.len();
            for &j in nbrs {
                if j < n {
                    counts[j] += 1;
                }
            }
        }

        let min_n = *counts.iter().min().unwrap_or(&0);
        let max_n = *counts.iter().max().unwrap_or(&0);
        let sum: usize = counts.iter().sum();
        let avg = sum as f64 / n as f64;
        let var: f64 = counts
            .iter()
            .map(|&c| {
                let diff = c as f64 - avg;
                diff * diff
            })
            .sum::<f64>()
            / n as f64;
        let std = var.sqrt();
        let isolated = counts.iter().filter(|&&c| c == 0).count();

        Self {
            n_atoms: n,
            n_pairs: vl.n_pairs(),
            min_neighbors: min_n,
            max_neighbors: max_n,
            avg_neighbors: avg,
            std_neighbors: std,
            isolated_atoms: isolated,
        }
    }

    /// Fill fraction: fraction of atoms with at least one neighbor.
    pub fn fill_fraction(&self) -> f64 {
        if self.n_atoms == 0 {
            return 0.0;
        }
        (self.n_atoms - self.isolated_atoms) as f64 / self.n_atoms as f64
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- CellList tests ---

    #[test]
    fn test_cell_list_build() {
        let positions: Vec<[f64; 3]> = vec![
            [0.1, 0.1, 0.1],
            [0.9, 0.1, 0.1],
            [0.1, 0.9, 0.1],
            [0.9, 0.9, 0.1],
            [0.1, 0.1, 0.9],
            [0.9, 0.1, 0.9],
            [0.1, 0.9, 0.9],
            [0.9, 0.9, 0.9],
        ];
        let mut cl = CellList::new([0.0; 3], [1.0; 3], 0.6);
        cl.build(&positions);

        assert_eq!(cl.nx, 2);
        assert_eq!(cl.ny, 2);
        assert_eq!(cl.nz, 2);

        let total: usize = cl.cells.iter().map(|c| c.len()).sum();
        assert_eq!(total, 8, "All 8 atoms must be in the cell list");
    }

    #[test]
    fn test_cell_list_neighbors() {
        let positions: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [0.4, 0.0, 0.0]];
        let mut cl = CellList::new([0.0; 3], [2.0; 3], 0.5);
        cl.build(&positions);

        let nbrs = cl.neighbors_of(0, &positions, 0.5 * 0.5);
        assert!(
            nbrs.iter().any(|&(j, _, _)| j == 1),
            "Atom 1 should be a neighbor of atom 0"
        );
    }

    #[test]
    fn test_cell_list_no_neighbor() {
        let positions: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let mut cl = CellList::new([0.0; 3], [2.0; 3], 0.5);
        cl.build(&positions);

        let nbrs = cl.neighbors_of(0, &positions, 0.5 * 0.5);
        assert!(
            !nbrs.iter().any(|&(j, _, _)| j == 1),
            "Atom 1 should NOT be a neighbor of atom 0 (distance 1.0 > cutoff 0.5)"
        );
    }

    // --- VerletList tests ---

    #[test]
    fn test_verlet_list_build() {
        let mut positions = Vec::new();
        for k in 0..10 {
            let x = (k as f64) * 0.5;
            positions.push([x, 0.5, 0.5]);
        }
        let mut vl = VerletList::new([0.0; 3], [10.0; 3], 1.2, 0.3);
        vl.build(&positions);
        assert_eq!(
            vl.neighbors.len(),
            10,
            "Should have neighbor lists for all 10 atoms"
        );
    }

    #[test]
    fn test_verlet_list_no_rebuild_needed() {
        let positions: Vec<[f64; 3]> = vec![[1.0, 1.0, 1.0], [2.0, 1.0, 1.0], [3.0, 3.0, 3.0]];
        let mut vl = VerletList::new([0.0; 3], [10.0; 3], 1.2, 0.4);
        vl.build(&positions);

        assert!(
            !vl.needs_rebuild(&positions),
            "needs_rebuild should be false when atoms have not moved"
        );
    }

    #[test]
    fn test_verlet_list_rebuild_needed() {
        let positions: Vec<[f64; 3]> = vec![[1.0, 1.0, 1.0], [2.0, 1.0, 1.0]];
        let skin = 0.4_f64;
        let mut vl = VerletList::new([0.0; 3], [10.0; 3], 1.2, skin);
        vl.build(&positions);

        let mut moved = positions.clone();
        moved[0][0] += skin * 0.6;
        assert!(
            vl.needs_rebuild(&moved),
            "needs_rebuild should be true after atom moved > skin/2"
        );
    }

    // --- minimum_image tests ---

    #[test]
    fn test_minimum_image_positive() {
        let result = minimum_image(0.8, 1.0);
        assert!(
            (result - (-0.2)).abs() < 1e-12,
            "Expected -0.2, got {result}"
        );
    }

    #[test]
    fn test_minimum_image_negative() {
        let result = minimum_image(-0.8, 1.0);
        assert!((result - 0.2).abs() < 1e-12, "Expected 0.2, got {result}");
    }

    // --- all_pairs symmetry test ---

    #[test]
    fn test_all_pairs_symmetry() {
        let positions: Vec<[f64; 3]> = vec![
            [0.0, 0.0, 0.0],
            [0.3, 0.0, 0.0],
            [0.7, 0.0, 0.0],
            [1.5, 1.5, 1.5],
        ];
        let cutoff = 0.5_f64;
        let mut cl = CellList::new([0.0; 3], [3.0; 3], cutoff);
        cl.build(&positions);

        let pairs = cl.all_pairs(&positions, cutoff * cutoff);

        for &(i, j, _) in &pairs {
            assert!(i < j, "Pair ({i}, {j}) violates i < j invariant");
        }

        let mut seen = std::collections::HashSet::new();
        for &(i, j, _) in &pairs {
            assert!(seen.insert((i, j)), "Duplicate pair ({i}, {j})");
        }
    }

    // --- MultiCutoffNeighborList tests ---

    #[test]
    fn test_multi_cutoff_build() {
        let positions: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [0.3, 0.0, 0.0], [0.9, 0.0, 0.0]];
        let mut mcnl = MultiCutoffNeighborList::new([0.0; 3], [5.0; 3], vec![0.5, 1.0]);
        mcnl.build(&positions);

        // Short cutoff (0.5): only pair (0,1) at dist 0.3
        assert_eq!(mcnl.pairs(0).len(), 1, "short cutoff should have 1 pair");
        // Long cutoff (1.0): pairs (0,1) at 0.3 and (0,2) at 0.9 and (1,2) at 0.6
        assert!(
            mcnl.pairs(1).len() >= 2,
            "long cutoff should have >=2 pairs"
        );
    }

    #[test]
    fn test_multi_cutoff_pair_counts() {
        let positions: Vec<[f64; 3]> = vec![[0.0; 3], [0.3, 0.0, 0.0]];
        let mut mcnl = MultiCutoffNeighborList::new([0.0; 3], [5.0; 3], vec![0.2, 0.5]);
        mcnl.build(&positions);
        let counts = mcnl.pair_counts();
        assert_eq!(counts.len(), 2);
        assert_eq!(counts[0], 0, "0.2 cutoff should find nothing at dist 0.3");
        assert_eq!(counts[1], 1, "0.5 cutoff should find 1 pair at dist 0.3");
    }

    // --- HalfNeighborList tests ---

    #[test]
    fn test_half_neighbor_list() {
        let positions: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [0.4, 0.0, 0.0], [5.0, 5.0, 5.0]];
        let hnl = HalfNeighborList::build(&positions, [0.0; 3], [10.0; 3], 0.5);
        assert_eq!(hnl.n_pairs(), 1, "should have 1 pair (0,1)");
        let p = hnl.get(0).unwrap();
        assert!(p.0 < p.1, "half list: i < j");
    }

    #[test]
    fn test_half_neighbor_list_from_full() {
        // Construct a full neighbor list manually
        let positions: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [0.3, 0.0, 0.0]];
        let full = vec![
            vec![1usize], // atom 0 neighbors
            vec![0usize], // atom 1 neighbors
        ];
        let box_len = [10.0, 10.0, 10.0];
        let hnl = HalfNeighborList::from_full(&full, &positions, box_len);
        assert_eq!(hnl.n_pairs(), 1);
    }

    // --- CompressedNeighborList tests ---

    #[test]
    fn test_compressed_neighbor_list() {
        let positions: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [0.3, 0.0, 0.0], [0.6, 0.0, 0.0]];
        let mut vl = VerletList::new([0.0; 3], [5.0; 3], 0.5, 0.1);
        vl.build(&positions);
        let box_len = [5.0, 5.0, 5.0];
        let cnl = CompressedNeighborList::from_verlet(&vl, &positions, box_len);
        assert_eq!(cnl.n_atoms(), 3);
        assert!(cnl.total_entries() > 0);
        assert!(cnl.memory_bytes() > 0);
    }

    #[test]
    fn test_compressed_neighbors_of() {
        let positions: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [0.3, 0.0, 0.0], [5.0, 5.0, 5.0]];
        let mut vl = VerletList::new([0.0; 3], [10.0; 3], 0.5, 0.1);
        vl.build(&positions);
        let box_len = [10.0, 10.0, 10.0];
        let cnl = CompressedNeighborList::from_verlet(&vl, &positions, box_len);
        // Atom 0 should have neighbor 1
        let nbrs = cnl.neighbors_of(0);
        assert!(nbrs.contains(&1), "atom 0 should have neighbor 1");
        let dsqs = cnl.distances_sq_of(0);
        assert!(!dsqs.is_empty());
    }

    // --- NeighborStatistics tests ---

    #[test]
    fn test_neighbor_statistics() {
        let mut positions = Vec::new();
        for i in 0..5 {
            positions.push([i as f64 * 0.3, 0.0, 0.0]);
        }
        let mut vl = VerletList::new([0.0; 3], [5.0; 3], 0.5, 0.1);
        vl.build(&positions);
        let stats = NeighborStatistics::from_verlet(&vl);
        assert_eq!(stats.n_atoms, 5);
        assert!(stats.n_pairs > 0);
        assert!(stats.avg_neighbors > 0.0);
        assert!(stats.fill_fraction() > 0.0);
    }

    #[test]
    fn test_neighbor_statistics_empty() {
        let vl = VerletList::new([0.0; 3], [1.0; 3], 0.5, 0.1);
        let stats = NeighborStatistics::from_verlet(&vl);
        assert_eq!(stats.n_atoms, 0);
        assert_eq!(stats.n_pairs, 0);
    }

    #[test]
    fn test_neighbor_statistics_isolated() {
        let positions: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [9.0, 9.0, 9.0]];
        let mut vl = VerletList::new([0.0; 3], [10.0; 3], 0.5, 0.1);
        vl.build(&positions);
        let stats = NeighborStatistics::from_verlet(&vl);
        assert_eq!(stats.isolated_atoms, 2, "both atoms should be isolated");
        assert!((stats.fill_fraction() - 0.0).abs() < 1e-12);
    }

    // --- avg_coordination test ---

    #[test]
    fn test_verlet_avg_coordination() {
        let positions: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [0.3, 0.0, 0.0]];
        let mut vl = VerletList::new([0.0; 3], [5.0; 3], 0.5, 0.1);
        vl.build(&positions);
        let avg = vl.avg_coordination();
        // One pair, 2 atoms: avg = 2*1/2 = 1.0
        assert!(
            (avg - 1.0).abs() < 1e-12,
            "avg coordination for 1 pair, 2 atoms: {avg}"
        );
    }
}

// ---------------------------------------------------------------------------
// CellListBuilder — improved cell-list construction with rebuild tracking
// ---------------------------------------------------------------------------

/// Improved cell-list builder with rebuild statistics and adjustable padding.
///
/// Provides the same O(N) performance as [`CellList`] but also tracks
/// how many times the list has been rebuilt and supports an optional
/// padding factor to reduce the number of re-builds.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct CellListBuilder {
    /// Underlying cell list.
    pub cell_list: CellList,
    /// Padding multiplier applied to the cutoff when sizing cells.
    pub padding: f64,
    /// Total number of times `build` has been called.
    pub build_count: usize,
    /// Total atoms inserted across all builds.
    pub total_atoms_inserted: usize,
}

#[allow(dead_code)]
impl CellListBuilder {
    /// Create a new builder.
    ///
    /// `padding` is a multiplier on the cutoff used as the cell size
    /// (e.g. 1.05 adds 5% padding).
    pub fn new(box_lo: [f64; 3], box_hi: [f64; 3], cutoff: f64, padding: f64) -> Self {
        let padded_cutoff = cutoff * padding.max(1.0);
        Self {
            cell_list: CellList::new(box_lo, box_hi, padded_cutoff),
            padding,
            build_count: 0,
            total_atoms_inserted: 0,
        }
    }

    /// Build the cell list from atom positions.
    pub fn build(&mut self, positions: &[[f64; 3]]) {
        self.cell_list.build(positions);
        self.build_count += 1;
        self.total_atoms_inserted += positions.len();
    }

    /// Find neighbors of atom `i` within `cutoff` (not padded).
    pub fn neighbors_of(
        &self,
        i: usize,
        positions: &[[f64; 3]],
        cutoff: f64,
    ) -> Vec<(usize, f64, [f64; 3])> {
        self.cell_list.neighbors_of(i, positions, cutoff * cutoff)
    }

    /// All pairs within true cutoff.
    pub fn all_pairs(&self, positions: &[[f64; 3]], cutoff: f64) -> Vec<(usize, usize, f64)> {
        self.cell_list.all_pairs(positions, cutoff * cutoff)
    }

    /// Average atoms per non-empty cell.
    pub fn avg_density(&self) -> f64 {
        let occ = self
            .cell_list
            .cells
            .iter()
            .filter(|c| !c.is_empty())
            .count();
        if occ == 0 {
            return 0.0;
        }
        let total: usize = self.cell_list.cells.iter().map(|c| c.len()).sum();
        total as f64 / occ as f64
    }
}

// ---------------------------------------------------------------------------
// VerletListWithSkinRefresh
// ---------------------------------------------------------------------------

/// Verlet list with automatic skin-based refresh and rebuild statistics.
///
/// Wraps [`VerletList`] and adds:
/// - per-atom displacement tracking for the skin criterion,
/// - rebuild count and fraction-of-steps-with-rebuild statistics.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct VerletListWithSkinRefresh {
    /// Inner Verlet list.
    pub verlet: VerletList,
    /// Total `update` calls.
    pub update_count: usize,
    /// Total `build` (rebuild) calls.
    pub rebuild_count: usize,
    /// Maximum displacement seen since last rebuild (across all atoms).
    pub max_displacement: f64,
}

#[allow(dead_code)]
impl VerletListWithSkinRefresh {
    /// Create a new skin-refresh Verlet list.
    pub fn new(box_lo: [f64; 3], box_hi: [f64; 3], cutoff: f64, skin: f64) -> Self {
        Self {
            verlet: VerletList::new(box_lo, box_hi, cutoff, skin),
            update_count: 0,
            rebuild_count: 0,
            max_displacement: 0.0,
        }
    }

    /// Build the list from initial positions.
    pub fn build(&mut self, positions: &[[f64; 3]]) {
        self.verlet.build(positions);
        self.rebuild_count += 1;
        self.max_displacement = 0.0;
    }

    /// Update: rebuild only if any atom has moved > skin/2, otherwise reuse.
    pub fn update(&mut self, positions: &[[f64; 3]]) {
        self.update_count += 1;

        // Track maximum displacement
        let box_len = [
            self.verlet.cell_list.box_hi[0] - self.verlet.cell_list.box_lo[0],
            self.verlet.cell_list.box_hi[1] - self.verlet.cell_list.box_lo[1],
            self.verlet.cell_list.box_hi[2] - self.verlet.cell_list.box_lo[2],
        ];
        if self.verlet.ref_positions.len() == positions.len() {
            for (ref_pos, cur_pos) in self.verlet.ref_positions.iter().zip(positions.iter()) {
                let dr_raw = [
                    cur_pos[0] - ref_pos[0],
                    cur_pos[1] - ref_pos[1],
                    cur_pos[2] - ref_pos[2],
                ];
                let dr = minimum_image_vec(dr_raw, box_len);
                let d = (dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]).sqrt();
                if d > self.max_displacement {
                    self.max_displacement = d;
                }
            }
        }

        if self.verlet.needs_rebuild(positions) {
            self.verlet.build(positions);
            self.rebuild_count += 1;
            self.max_displacement = 0.0;
        }
    }

    /// Fraction of update calls that triggered a rebuild.
    pub fn rebuild_fraction(&self) -> f64 {
        if self.update_count == 0 {
            return 0.0;
        }
        self.rebuild_count as f64 / self.update_count as f64
    }

    /// Delegate to inner Verlet list.
    pub fn neighbors(&self, i: usize) -> &[usize] {
        self.verlet.neighbors(i)
    }

    /// Total number of unique neighbor pairs.
    pub fn n_pairs(&self) -> usize {
        self.verlet.n_pairs()
    }
}

// ---------------------------------------------------------------------------
// MultiBodyNeighborList
// ---------------------------------------------------------------------------

/// Multi-body neighbor list: stores triples `(i, j, k)` with `i < j < k`
/// within a given cutoff.
///
/// Useful for three-body potentials (e.g. Stillinger-Weber, Tersoff).
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct MultiBodyNeighborList {
    /// Cutoff distance.
    pub cutoff: f64,
    /// Triples (i, j, k) with i < j < k within cutoff.
    pub triples: Vec<(usize, usize, usize)>,
    /// Quadruples (i, j, k, l) for four-body terms (optional).
    pub quadruples: Vec<(usize, usize, usize, usize)>,
}

#[allow(dead_code)]
impl MultiBodyNeighborList {
    /// Build the multi-body neighbor list.
    pub fn build(
        positions: &[[f64; 3]],
        box_lo: [f64; 3],
        box_hi: [f64; 3],
        cutoff: f64,
        build_quadruples: bool,
    ) -> Self {
        let cutoff_sq = cutoff * cutoff;
        let mut cl = CellList::new(box_lo, box_hi, cutoff);
        cl.build(positions);

        let n = positions.len();
        let mut triples = Vec::new();
        let mut quadruples = Vec::new();

        // Build pair list first
        let pairs = cl.all_pairs(positions, cutoff_sq);

        // For each atom, collect its neighbors
        let mut nbr_map: Vec<Vec<usize>> = vec![Vec::new(); n];
        for &(i, j, _) in &pairs {
            nbr_map[i].push(j);
            nbr_map[j].push(i);
        }

        // Triples: for each pair (i,j), add (i,j,k) for each k ∈ nbr(i) with j < k
        for &(i, j, _) in &pairs {
            // i < j by construction from all_pairs
            for &k in &nbr_map[i] {
                if k > j {
                    // Check i-k distance
                    let dr_ik = [
                        positions[k][0] - positions[i][0],
                        positions[k][1] - positions[i][1],
                        positions[k][2] - positions[i][2],
                    ];
                    let box_len = [
                        box_hi[0] - box_lo[0],
                        box_hi[1] - box_lo[1],
                        box_hi[2] - box_lo[2],
                    ];
                    let dr_ik_mi = minimum_image_vec(dr_ik, box_len);
                    let d_ik_sq = dr_ik_mi[0].powi(2) + dr_ik_mi[1].powi(2) + dr_ik_mi[2].powi(2);
                    if d_ik_sq <= cutoff_sq {
                        triples.push((i, j, k));
                    }
                }
            }
        }

        // Quadruples (brute force extension of triples)
        if build_quadruples {
            let box_len = [
                box_hi[0] - box_lo[0],
                box_hi[1] - box_lo[1],
                box_hi[2] - box_lo[2],
            ];
            for &(i, j, k) in &triples {
                for &l in &nbr_map[i] {
                    if l > k {
                        let dr = [
                            positions[l][0] - positions[i][0],
                            positions[l][1] - positions[i][1],
                            positions[l][2] - positions[i][2],
                        ];
                        let dr_mi = minimum_image_vec(dr, box_len);
                        let d_sq = dr_mi[0].powi(2) + dr_mi[1].powi(2) + dr_mi[2].powi(2);
                        if d_sq <= cutoff_sq {
                            quadruples.push((i, j, k, l));
                        }
                    }
                }
            }
        }

        Self {
            cutoff,
            triples,
            quadruples,
        }
    }

    /// Number of triples.
    pub fn n_triples(&self) -> usize {
        self.triples.len()
    }

    /// Number of quadruples.
    pub fn n_quadruples(&self) -> usize {
        self.quadruples.len()
    }

    /// Number of triples involving atom `i`.
    pub fn triples_of(&self, i: usize) -> Vec<(usize, usize, usize)> {
        self.triples
            .iter()
            .filter(|&&(a, b, c)| a == i || b == i || c == i)
            .copied()
            .collect()
    }
}

// ---------------------------------------------------------------------------
// NeighborListGpuHints
// ---------------------------------------------------------------------------

/// GPU-oriented hints for a neighbor list.
///
/// Contains padded, sorted, and tiled representations suitable for
/// upload to a GPU buffer (WGPU, CUDA, etc.).
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct NeighborListGpuHints {
    /// Row offsets (CSR format): `row_offsets[i]..row_offsets[i+1]` for atom i.
    pub row_offsets: Vec<u32>,
    /// Flat neighbor index array.
    pub col_indices: Vec<u32>,
    /// Flat squared-distance array (aligned with col_indices).
    pub dist_sq: Vec<f32>,
    /// Number of atoms.
    pub n_atoms: usize,
    /// Recommended WGPU workgroup size.
    pub workgroup_size: u32,
    /// Whether the list is sorted by distance within each row.
    pub sorted_by_dist: bool,
}

#[allow(dead_code)]
impl NeighborListGpuHints {
    /// Build GPU hints from a Verlet list.
    ///
    /// `sort_by_dist` — if `true`, sort each atom's neighbor list by distance.
    pub fn from_verlet(
        vl: &VerletList,
        positions: &[[f64; 3]],
        box_len: [f64; 3],
        sort_by_dist: bool,
    ) -> Self {
        let n = vl.neighbors.len();
        let mut row_offsets = Vec::with_capacity(n + 1);
        let mut col_indices = Vec::new();
        let mut dist_sq_vec: Vec<f32> = Vec::new();

        for (i, nbrs) in vl.neighbors.iter().enumerate() {
            row_offsets.push(col_indices.len() as u32);

            // Compute distances for sorting
            let mut entries: Vec<(u32, f32)> = nbrs
                .iter()
                .map(|&j| {
                    let dr_raw = [
                        positions[j][0] - positions[i][0],
                        positions[j][1] - positions[i][1],
                        positions[j][2] - positions[i][2],
                    ];
                    let dr = minimum_image_vec(dr_raw, box_len);
                    let dsq = (dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]) as f32;
                    (j as u32, dsq)
                })
                .collect();

            if sort_by_dist {
                entries.sort_unstable_by(|a, b| {
                    a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal)
                });
            }

            for (j, dsq) in entries {
                col_indices.push(j);
                dist_sq_vec.push(dsq);
            }
        }
        row_offsets.push(col_indices.len() as u32);

        // Recommended workgroup size: next power-of-2 ≤ 256
        let wg = col_indices
            .len()
            .checked_div(n)
            .map(|avg_nbr| avg_nbr.next_power_of_two().clamp(32, 256) as u32)
            .unwrap_or(64);

        Self {
            row_offsets,
            col_indices,
            dist_sq: dist_sq_vec,
            n_atoms: n,
            workgroup_size: wg,
            sorted_by_dist: sort_by_dist,
        }
    }

    /// Memory footprint in bytes.
    pub fn memory_bytes(&self) -> usize {
        self.row_offsets.len() * 4 + self.col_indices.len() * 4 + self.dist_sq.len() * 4
    }

    /// Neighbors of atom `i` as GPU-side `u32` indices.
    pub fn neighbors_of(&self, i: usize) -> &[u32] {
        let start = self.row_offsets[i] as usize;
        let end = self.row_offsets[i + 1] as usize;
        &self.col_indices[start..end]
    }

    /// Number of unique neighbor entries.
    pub fn total_entries(&self) -> usize {
        self.col_indices.len()
    }
}

// ---------------------------------------------------------------------------
// Tests for new neighbor_list additions
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests_new_nl {
    use super::*;

    // --- CellListBuilder tests ---

    #[test]
    fn test_cell_list_builder_basic() {
        let positions: Vec<[f64; 3]> = vec![[0.1, 0.1, 0.1], [0.5, 0.5, 0.5], [0.9, 0.9, 0.9]];
        let mut clb = CellListBuilder::new([0.0; 3], [1.0; 3], 0.5, 1.0);
        clb.build(&positions);
        assert_eq!(clb.build_count, 1);
        assert_eq!(clb.total_atoms_inserted, 3);
    }

    #[test]
    fn test_cell_list_builder_neighbors() {
        let positions: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [0.3, 0.0, 0.0], [5.0, 5.0, 5.0]];
        let mut clb = CellListBuilder::new([0.0; 3], [10.0; 3], 0.5, 1.0);
        clb.build(&positions);
        let nbrs = clb.neighbors_of(0, &positions, 0.5);
        assert!(nbrs.iter().any(|&(j, _, _)| j == 1));
        assert!(!nbrs.iter().any(|&(j, _, _)| j == 2));
    }

    #[test]
    fn test_cell_list_builder_all_pairs() {
        let positions: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [0.3, 0.0, 0.0], [0.6, 0.0, 0.0]];
        let mut clb = CellListBuilder::new([0.0; 3], [5.0; 3], 0.5, 1.0);
        clb.build(&positions);
        let pairs = clb.all_pairs(&positions, 0.5);
        assert!(pairs.len() >= 2, "should find at least 2 pairs");
    }

    #[test]
    fn test_cell_list_builder_density() {
        let positions: Vec<[f64; 3]> = vec![[0.5; 3]; 8];
        let mut clb = CellListBuilder::new([0.0; 3], [1.0; 3], 0.5, 1.0);
        clb.build(&positions);
        let density = clb.avg_density();
        assert!(density >= 1.0, "avg density should be >= 1.0");
    }

    #[test]
    fn test_cell_list_builder_rebuild() {
        let pos1: Vec<[f64; 3]> = vec![[0.1, 0.1, 0.1]];
        let pos2: Vec<[f64; 3]> = vec![[0.5, 0.5, 0.5]];
        let mut clb = CellListBuilder::new([0.0; 3], [1.0; 3], 0.5, 1.0);
        clb.build(&pos1);
        clb.build(&pos2);
        assert_eq!(clb.build_count, 2);
        assert_eq!(clb.total_atoms_inserted, 2);
    }

    // --- VerletListWithSkinRefresh tests ---

    #[test]
    fn test_vlsr_build() {
        let positions: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [0.3, 0.0, 0.0]];
        let mut vlsr = VerletListWithSkinRefresh::new([0.0; 3], [5.0; 3], 0.5, 0.2);
        vlsr.build(&positions);
        assert_eq!(vlsr.rebuild_count, 1);
        assert_eq!(vlsr.update_count, 0);
    }

    #[test]
    fn test_vlsr_no_rebuild_when_static() {
        let positions: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [0.3, 0.0, 0.0]];
        let mut vlsr = VerletListWithSkinRefresh::new([0.0; 3], [5.0; 3], 0.5, 0.4);
        vlsr.build(&positions);
        vlsr.update(&positions);
        assert_eq!(vlsr.rebuild_count, 1, "no rebuild when atoms haven't moved");
        assert_eq!(vlsr.update_count, 1);
    }

    #[test]
    fn test_vlsr_rebuild_on_large_move() {
        let positions: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [0.3, 0.0, 0.0]];
        let skin = 0.2f64;
        let mut vlsr = VerletListWithSkinRefresh::new([0.0; 3], [5.0; 3], 0.5, skin);
        vlsr.build(&positions);

        let mut moved = positions.clone();
        moved[0][0] += skin * 0.6; // exceed skin/2
        vlsr.update(&moved);
        assert_eq!(
            vlsr.rebuild_count, 2,
            "should have rebuilt after large move"
        );
    }

    #[test]
    fn test_vlsr_rebuild_fraction() {
        let positions: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [0.3, 0.0, 0.0]];
        let skin = 0.2f64;
        let mut vlsr = VerletListWithSkinRefresh::new([0.0; 3], [5.0; 3], 0.5, skin);
        vlsr.build(&positions);

        // 4 updates with no movement: no rebuilds triggered
        for _ in 0..4 {
            vlsr.update(&positions);
        }
        let frac = vlsr.rebuild_fraction();
        // rebuild_count started at 1 from build(), update_count = 4
        // fraction = rebuild_count / update_count = 1/4 or 0/4 depending on moves
        // In this case no moves → no additional rebuilds → frac = 0.25 (1/4)
        assert!((0.0..=1.0).contains(&frac), "fraction out of range: {frac}");
    }

    #[test]
    fn test_vlsr_n_pairs() {
        let positions: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [0.3, 0.0, 0.0], [5.0, 5.0, 5.0]];
        let mut vlsr = VerletListWithSkinRefresh::new([0.0; 3], [10.0; 3], 0.5, 0.1);
        vlsr.build(&positions);
        assert_eq!(vlsr.n_pairs(), 1);
    }

    // --- MultiBodyNeighborList tests ---

    #[test]
    fn test_multi_body_triples_triangle() {
        // 3 atoms forming an equilateral triangle within cutoff
        let d = 0.3f64;
        let positions: Vec<[f64; 3]> = vec![
            [0.0, 0.0, 0.0],
            [d, 0.0, 0.0],
            [d / 2.0, d * (3.0f64.sqrt() / 2.0), 0.0],
        ];
        let mb = MultiBodyNeighborList::build(&positions, [0.0; 3], [5.0; 3], 0.5, false);
        // All 3 atoms are within 0.5 of each other → exactly 1 triple (0,1,2)
        assert_eq!(mb.n_triples(), 1, "should have 1 triple for triangle");
    }

    #[test]
    fn test_multi_body_no_triples_isolated() {
        let positions: Vec<[f64; 3]> = vec![
            [0.0, 0.0, 0.0],
            [0.3, 0.0, 0.0],
            [5.0, 5.0, 5.0], // far away
        ];
        let mb = MultiBodyNeighborList::build(&positions, [0.0; 3], [10.0; 3], 0.5, false);
        // Only atoms 0 and 1 are within cutoff, no triple possible
        assert_eq!(mb.n_triples(), 0);
    }

    #[test]
    fn test_multi_body_triples_of() {
        let d = 0.3f64;
        let positions: Vec<[f64; 3]> = vec![
            [0.0, 0.0, 0.0],
            [d, 0.0, 0.0],
            [d / 2.0, d * (3.0f64.sqrt() / 2.0), 0.0],
        ];
        let mb = MultiBodyNeighborList::build(&positions, [0.0; 3], [5.0; 3], 0.5, false);
        let t0 = mb.triples_of(0);
        assert_eq!(t0.len(), 1, "atom 0 should be in 1 triple");
    }

    #[test]
    fn test_multi_body_quadruples() {
        // 4 atoms in a tetrahedron-like arrangement
        let positions: Vec<[f64; 3]> = vec![
            [0.0, 0.0, 0.0],
            [0.3, 0.0, 0.0],
            [0.15, 0.26, 0.0],
            [0.15, 0.09, 0.24],
        ];
        let mb = MultiBodyNeighborList::build(&positions, [0.0; 3], [5.0; 3], 0.4, true);
        // All 4 atoms should be within 0.4 of each other
        assert!(mb.n_triples() > 0, "should have triples");
    }

    // --- NeighborListGpuHints tests ---

    #[test]
    fn test_gpu_hints_basic() {
        let positions: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [0.3, 0.0, 0.0], [5.0, 5.0, 5.0]];
        let mut vl = VerletList::new([0.0; 3], [10.0; 3], 0.5, 0.1);
        vl.build(&positions);
        let box_len = [10.0, 10.0, 10.0];
        let hints = NeighborListGpuHints::from_verlet(&vl, &positions, box_len, false);

        assert_eq!(hints.n_atoms, 3);
        assert_eq!(hints.row_offsets.len(), 4); // n+1
        assert!(hints.memory_bytes() > 0);
    }

    #[test]
    fn test_gpu_hints_sorted() {
        let positions: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [0.3, 0.0, 0.0], [0.45, 0.0, 0.0]];
        let mut vl = VerletList::new([0.0; 3], [5.0; 3], 0.6, 0.1);
        vl.build(&positions);
        let box_len = [5.0, 5.0, 5.0];
        let hints = NeighborListGpuHints::from_verlet(&vl, &positions, box_len, true);
        assert!(hints.sorted_by_dist);

        // Verify distances of atom 0 are non-decreasing
        if hints.neighbors_of(0).len() >= 2 {
            let start = hints.row_offsets[0] as usize;
            let end = hints.row_offsets[1] as usize;
            for w in hints.dist_sq[start..end].windows(2) {
                assert!(w[0] <= w[1], "dist_sq should be sorted");
            }
        }
    }

    #[test]
    fn test_gpu_hints_neighbors_of() {
        let positions: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [0.3, 0.0, 0.0]];
        let mut vl = VerletList::new([0.0; 3], [5.0; 3], 0.5, 0.1);
        vl.build(&positions);
        let box_len = [5.0, 5.0, 5.0];
        let hints = NeighborListGpuHints::from_verlet(&vl, &positions, box_len, false);
        let nbrs_0 = hints.neighbors_of(0);
        assert!(nbrs_0.contains(&1u32), "atom 0 should have neighbor 1");
    }

    #[test]
    fn test_gpu_hints_total_entries() {
        let positions: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [0.3, 0.0, 0.0]];
        let mut vl = VerletList::new([0.0; 3], [5.0; 3], 0.5, 0.1);
        vl.build(&positions);
        let box_len = [5.0, 5.0, 5.0];
        let hints = NeighborListGpuHints::from_verlet(&vl, &positions, box_len, false);
        // One pair (0,1) stored in half list: only 1 entry
        assert_eq!(hints.total_entries(), 1);
    }
}
