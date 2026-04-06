#![allow(clippy::manual_div_ceil)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Neighbor list structures for efficient pair interaction computation.
//!
//! Provides cell-based and Verlet-list approaches, along with periodic
//! boundary condition (PBC) utilities, neighbor list statistics,
//! multi-cutoff neighbor lists, and parallelization hints.

use oxiphysics_core::math::Vec3;

// ---------------------------------------------------------------------------
// Periodic box
// ---------------------------------------------------------------------------

/// A rectangular periodic simulation box.
#[derive(Debug, Clone)]
pub struct PeriodicBox {
    /// Box dimensions (Lx, Ly, Lz).
    pub dims: Vec3,
}

impl PeriodicBox {
    /// Create a new periodic box with the given dimensions.
    pub fn new(dims: Vec3) -> Self {
        Self { dims }
    }

    /// Create a cubic periodic box.
    pub fn cubic(l: f64) -> Self {
        Self {
            dims: Vec3::new(l, l, l),
        }
    }

    /// Volume of the box.
    pub fn volume(&self) -> f64 {
        self.dims.x * self.dims.y * self.dims.z
    }

    /// Apply minimum image convention to a displacement vector.
    ///
    /// Returns the wrapped displacement such that each component is in
    /// \[-L/2, L/2).
    pub fn minimum_image(&self, dr: &Vec3) -> Vec3 {
        Vec3::new(
            dr.x - self.dims.x * (dr.x / self.dims.x).round(),
            dr.y - self.dims.y * (dr.y / self.dims.y).round(),
            dr.z - self.dims.z * (dr.z / self.dims.z).round(),
        )
    }

    /// Wrap a position back into the primary image \[0, L).
    pub fn wrap_position(&self, pos: &Vec3) -> Vec3 {
        Vec3::new(
            pos.x - self.dims.x * (pos.x / self.dims.x).floor(),
            pos.y - self.dims.y * (pos.y / self.dims.y).floor(),
            pos.z - self.dims.z * (pos.z / self.dims.z).floor(),
        )
    }

    /// Number density: N / V.
    #[allow(dead_code)]
    pub fn number_density(&self, n: usize) -> f64 {
        n as f64 / self.volume()
    }
}

/// Compute the displacement and distance between two positions under PBC.
///
/// Returns `(dr, |dr|)` using the minimum image convention.
pub fn distance_pbc(a: &Vec3, b: &Vec3, pbox: &PeriodicBox) -> (Vec3, f64) {
    let dr = pbox.minimum_image(&(b - a));
    let dist = dr.norm();
    (dr, dist)
}

// ---------------------------------------------------------------------------
// Cell list
// ---------------------------------------------------------------------------

/// Cell-based neighbor search for periodic boxes.
///
/// The simulation box is divided into cells of size >= cutoff.
/// Only atoms in the same or neighboring cells need to be checked.
#[derive(Debug, Clone)]
pub struct CellList {
    /// Number of cells in each dimension.
    pub n_cells: [usize; 3],
    /// Cell size in each dimension.
    pub cell_size: Vec3,
    /// Atom indices in each cell (flattened 3D grid).
    pub cells: Vec<Vec<usize>>,
    /// Cutoff distance used to size the cells.
    pub cutoff: f64,
}

impl CellList {
    /// Build a new cell list from the given positions and box.
    pub fn build(positions: &[Vec3], pbox: &PeriodicBox, cutoff: f64) -> Self {
        let nx = (pbox.dims.x / cutoff).floor().max(1.0) as usize;
        let ny = (pbox.dims.y / cutoff).floor().max(1.0) as usize;
        let nz = (pbox.dims.z / cutoff).floor().max(1.0) as usize;

        let cell_size = Vec3::new(
            pbox.dims.x / nx as f64,
            pbox.dims.y / ny as f64,
            pbox.dims.z / nz as f64,
        );

        let total_cells = nx * ny * nz;
        let mut cells = vec![Vec::new(); total_cells];

        for (i, pos) in positions.iter().enumerate() {
            let wrapped = pbox.wrap_position(pos);
            let cx = ((wrapped.x / cell_size.x) as usize).min(nx - 1);
            let cy = ((wrapped.y / cell_size.y) as usize).min(ny - 1);
            let cz = ((wrapped.z / cell_size.z) as usize).min(nz - 1);
            let idx = cx * ny * nz + cy * nz + cz;
            cells[idx].push(i);
        }

        Self {
            n_cells: [nx, ny, nz],
            cell_size,
            cells,
            cutoff,
        }
    }

    /// Flat index from 3D cell coordinates (with periodic wrapping).
    fn cell_index(&self, ix: isize, iy: isize, iz: isize) -> usize {
        let nx = self.n_cells[0] as isize;
        let ny = self.n_cells[1] as isize;
        let nz = self.n_cells[2] as isize;
        let cx = ((ix % nx + nx) % nx) as usize;
        let cy = ((iy % ny + ny) % ny) as usize;
        let cz = ((iz % nz + nz) % nz) as usize;
        cx * (ny as usize) * (nz as usize) + cy * (nz as usize) + cz
    }

    /// Get the neighbor pairs (i < j) within the cutoff distance.
    ///
    /// Returns a vector of `(i, j, dr, dist)` tuples where `dr` points from i to j.
    /// Each pair appears exactly once.
    pub fn neighbor_pairs(
        &self,
        positions: &[Vec3],
        pbox: &PeriodicBox,
    ) -> Vec<(usize, usize, Vec3, f64)> {
        let mut pairs = Vec::new();
        let [nx, ny, nz] = self.n_cells;
        let total = nx * ny * nz;

        // Track which cell pairs have already been processed to avoid
        // duplicates from periodic wrapping (e.g., with 2 cells, offsets
        // +1 and -1 wrap to the same neighbor).
        let mut visited_cell_pairs = std::collections::HashSet::new();

        for cell_a in 0..total {
            let (ax, ay, az) = self.unflatten(cell_a);

            // Self-cell pairs (i < j)
            let atoms_a = &self.cells[cell_a];
            for (ai, &i) in atoms_a.iter().enumerate() {
                for &j in &atoms_a[ai + 1..] {
                    let (dr, dist) = distance_pbc(&positions[i], &positions[j], pbox);
                    if dist < self.cutoff {
                        pairs.push((i, j, dr, dist));
                    }
                }
            }

            // Neighbor cells
            for dx in -1_isize..=1 {
                for dy in -1_isize..=1 {
                    for dz in -1_isize..=1 {
                        if dx == 0 && dy == 0 && dz == 0 {
                            continue;
                        }
                        let cell_b =
                            self.cell_index(ax as isize + dx, ay as isize + dy, az as isize + dz);
                        if cell_b == cell_a {
                            continue; // same cell after wrapping, already handled
                        }
                        // Canonical pair to avoid processing (a,b) and (b,a)
                        let key = if cell_a < cell_b {
                            (cell_a, cell_b)
                        } else {
                            (cell_b, cell_a)
                        };
                        if !visited_cell_pairs.insert(key) {
                            continue; // already processed this cell pair
                        }
                        for &i in &self.cells[cell_a] {
                            for &j in &self.cells[cell_b] {
                                let (a, b) = if i < j { (i, j) } else { (j, i) };
                                let (dr, dist) = distance_pbc(&positions[a], &positions[b], pbox);
                                if dist < self.cutoff {
                                    pairs.push((a, b, dr, dist));
                                }
                            }
                        }
                    }
                }
            }
        }
        pairs
    }

    /// Convert flat cell index back to 3D coordinates.
    fn unflatten(&self, idx: usize) -> (usize, usize, usize) {
        let ny = self.n_cells[1];
        let nz = self.n_cells[2];
        let x = idx / (ny * nz);
        let rem = idx % (ny * nz);
        let y = rem / nz;
        let z = rem % nz;
        (x, y, z)
    }

    /// Total number of cells.
    #[allow(dead_code)]
    pub fn total_cells(&self) -> usize {
        self.n_cells[0] * self.n_cells[1] * self.n_cells[2]
    }

    /// Maximum number of atoms in any single cell.
    #[allow(dead_code)]
    pub fn max_cell_occupancy(&self) -> usize {
        self.cells.iter().map(|c| c.len()).max().unwrap_or(0)
    }

    /// Average number of atoms per cell.
    #[allow(dead_code)]
    pub fn avg_cell_occupancy(&self) -> f64 {
        let total: usize = self.cells.iter().map(|c| c.len()).sum();
        let n_cells = self.cells.len();
        if n_cells == 0 {
            0.0
        } else {
            total as f64 / n_cells as f64
        }
    }

    /// Number of non-empty cells.
    #[allow(dead_code)]
    pub fn occupied_cells(&self) -> usize {
        self.cells.iter().filter(|c| !c.is_empty()).count()
    }
}

// ---------------------------------------------------------------------------
// Verlet list
// ---------------------------------------------------------------------------

/// Build a Verlet neighbor list from positions.
///
/// This is a convenience free function that constructs a [`VerletList`] in one
/// call.  Internally it delegates to [`VerletList::build`].
///
/// # Arguments
/// * `positions` -- Slice of `[f64;3]` atom positions.
/// * `pbox`      -- Periodic simulation box.
/// * `cutoff`    -- Interaction cutoff radius.
/// * `skin`      -- Verlet skin distance (rebuild trigger threshold = skin/2).
///
/// # Returns
/// A freshly built [`VerletList`] ready for force evaluation.
pub fn build_verlet_list(
    positions: &[Vec3],
    pbox: &PeriodicBox,
    cutoff: f64,
    skin: f64,
) -> VerletList {
    VerletList::build(positions, pbox, cutoff, skin)
}

/// Verlet neighbor list with skin distance.
///
/// Neighbor pairs are stored and reused until any atom moves more than
/// half the skin distance, at which point the list must be rebuilt.
#[derive(Debug, Clone)]
pub struct VerletList {
    /// Neighbor pairs (i, j) within cutoff + skin.
    pub pairs: Vec<(usize, usize)>,
    /// Positions when the list was last built.
    pub reference_positions: Vec<Vec3>,
    /// Cutoff distance.
    pub cutoff: f64,
    /// Skin distance.
    pub skin: f64,
    /// Number of times the list has been rebuilt.
    pub rebuild_count: u64,
}

impl VerletList {
    /// Build a new Verlet list.
    pub fn build(positions: &[Vec3], pbox: &PeriodicBox, cutoff: f64, skin: f64) -> Self {
        let effective_cutoff = cutoff + skin;
        let cell_list = CellList::build(positions, pbox, effective_cutoff);
        let all_pairs = cell_list.neighbor_pairs(positions, pbox);

        let pairs: Vec<(usize, usize)> = all_pairs
            .into_iter()
            .filter(|(_, _, _, dist)| *dist < effective_cutoff)
            .map(|(i, j, _, _)| (i, j))
            .collect();

        Self {
            pairs,
            reference_positions: positions.to_vec(),
            cutoff,
            skin,
            rebuild_count: 0,
        }
    }

    /// Rebuild the Verlet list with new positions.
    #[allow(dead_code)]
    pub fn rebuild(&mut self, positions: &[Vec3], pbox: &PeriodicBox) {
        let new_list = Self::build(positions, pbox, self.cutoff, self.skin);
        self.pairs = new_list.pairs;
        self.reference_positions = new_list.reference_positions;
        self.rebuild_count += 1;
    }

    /// Check whether any atom has moved more than half the skin distance.
    ///
    /// Returns `true` if a rebuild is needed.
    pub fn needs_rebuild(&self, positions: &[Vec3], pbox: &PeriodicBox) -> bool {
        let half_skin = self.skin * 0.5;
        for (ref_pos, cur_pos) in self.reference_positions.iter().zip(positions.iter()) {
            let (_, dist) = distance_pbc(ref_pos, cur_pos, pbox);
            if dist > half_skin {
                return true;
            }
        }
        false
    }

    /// Check and rebuild if needed. Returns true if rebuilt.
    #[allow(dead_code)]
    pub fn update_if_needed(&mut self, positions: &[Vec3], pbox: &PeriodicBox) -> bool {
        if self.needs_rebuild(positions, pbox) {
            self.rebuild(positions, pbox);
            true
        } else {
            false
        }
    }

    /// Number of pairs in the list.
    #[allow(dead_code)]
    pub fn pair_count(&self) -> usize {
        self.pairs.len()
    }

    /// Maximum displacement since last build.
    #[allow(dead_code)]
    pub fn max_displacement(&self, positions: &[Vec3], pbox: &PeriodicBox) -> f64 {
        self.reference_positions
            .iter()
            .zip(positions.iter())
            .map(|(ref_pos, cur_pos)| {
                let (_, dist) = distance_pbc(ref_pos, cur_pos, pbox);
                dist
            })
            .fold(0.0_f64, f64::max)
    }
}

// ---------------------------------------------------------------------------
// NeighborListStats
// ---------------------------------------------------------------------------

/// Statistics about a neighbor list.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct NeighborListStats {
    /// Total number of atoms.
    pub n_atoms: usize,
    /// Total number of neighbor pairs.
    pub n_pairs: usize,
    /// Average number of neighbors per atom.
    pub avg_neighbors: f64,
    /// Maximum number of neighbors for any atom.
    pub max_neighbors: usize,
    /// Minimum number of neighbors for any atom.
    pub min_neighbors: usize,
}

#[allow(dead_code)]
impl NeighborListStats {
    /// Compute statistics from a Verlet list.
    pub fn from_verlet_list(vlist: &VerletList, n_atoms: usize) -> Self {
        let mut neighbor_counts = vec![0_usize; n_atoms];
        for &(i, j) in &vlist.pairs {
            if i < n_atoms {
                neighbor_counts[i] += 1;
            }
            if j < n_atoms {
                neighbor_counts[j] += 1;
            }
        }
        let max_neighbors = neighbor_counts.iter().copied().max().unwrap_or(0);
        let min_neighbors = neighbor_counts.iter().copied().min().unwrap_or(0);
        let total: usize = neighbor_counts.iter().sum();
        let avg_neighbors = if n_atoms > 0 {
            total as f64 / n_atoms as f64
        } else {
            0.0
        };

        Self {
            n_atoms,
            n_pairs: vlist.pairs.len(),
            avg_neighbors,
            max_neighbors,
            min_neighbors,
        }
    }

    /// Compute statistics from a cell list.
    pub fn from_cell_list(cell_list: &CellList, positions: &[Vec3], pbox: &PeriodicBox) -> Self {
        let pairs = cell_list.neighbor_pairs(positions, pbox);
        let n_atoms = positions.len();
        let mut neighbor_counts = vec![0_usize; n_atoms];
        for &(i, j, _, _) in &pairs {
            if i < n_atoms {
                neighbor_counts[i] += 1;
            }
            if j < n_atoms {
                neighbor_counts[j] += 1;
            }
        }
        let max_neighbors = neighbor_counts.iter().copied().max().unwrap_or(0);
        let min_neighbors = neighbor_counts.iter().copied().min().unwrap_or(0);
        let total: usize = neighbor_counts.iter().sum();
        let avg_neighbors = if n_atoms > 0 {
            total as f64 / n_atoms as f64
        } else {
            0.0
        };

        Self {
            n_atoms,
            n_pairs: pairs.len(),
            avg_neighbors,
            max_neighbors,
            min_neighbors,
        }
    }
}

// ---------------------------------------------------------------------------
// MultiCutoffNeighborList
// ---------------------------------------------------------------------------

/// Neighbor list that supports multiple cutoff radii.
///
/// Useful for systems with different interaction ranges (e.g., short-range
/// repulsion + long-range electrostatics).
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct MultiCutoffNeighborList {
    /// Per-cutoff pair lists: `lists[k]` contains pairs within `cutoffs[k]`.
    pub lists: Vec<Vec<(usize, usize)>>,
    /// Cutoff radii (sorted ascending).
    pub cutoffs: Vec<f64>,
}

#[allow(dead_code)]
impl MultiCutoffNeighborList {
    /// Build multi-cutoff neighbor lists.
    ///
    /// Uses the largest cutoff for the cell list, then filters pairs
    /// into each cutoff bucket.
    pub fn build(positions: &[Vec3], pbox: &PeriodicBox, cutoffs: &[f64], skin: f64) -> Self {
        let mut sorted_cutoffs: Vec<f64> = cutoffs.to_vec();
        sorted_cutoffs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

        let max_cutoff = sorted_cutoffs.last().copied().unwrap_or(0.0);
        let effective = max_cutoff + skin;

        let cell_list = CellList::build(positions, pbox, effective);
        let all_pairs = cell_list.neighbor_pairs(positions, pbox);

        let mut lists: Vec<Vec<(usize, usize)>> = vec![Vec::new(); sorted_cutoffs.len()];

        for (i, j, _, dist) in all_pairs {
            for (k, &cutoff) in sorted_cutoffs.iter().enumerate() {
                if dist < cutoff + skin {
                    lists[k].push((i, j));
                    break; // Only add to the tightest matching cutoff
                }
            }
        }

        // Actually, we want each list to contain all pairs within its cutoff,
        // not exclusive. Rebuild:
        let cell_list2 = CellList::build(positions, pbox, effective);
        let all_pairs2 = cell_list2.neighbor_pairs(positions, pbox);
        let mut lists2: Vec<Vec<(usize, usize)>> = vec![Vec::new(); sorted_cutoffs.len()];
        for (i, j, _, dist) in all_pairs2 {
            for (k, &cutoff) in sorted_cutoffs.iter().enumerate() {
                if dist < cutoff + skin {
                    lists2[k].push((i, j));
                }
            }
        }

        Self {
            lists: lists2,
            cutoffs: sorted_cutoffs,
        }
    }

    /// Number of pairs at a given cutoff index.
    pub fn pair_count(&self, cutoff_idx: usize) -> usize {
        self.lists.get(cutoff_idx).map(|l| l.len()).unwrap_or(0)
    }
}

// ---------------------------------------------------------------------------
// Force decomposition for parallel execution
// ---------------------------------------------------------------------------

/// Domain decomposition hint for parallelizing force computations.
///
/// Splits the atom index space into non-overlapping domains suitable
/// for OpenMP-style parallel loops.  Each domain owns a contiguous
/// range of atom indices for force accumulation.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct ForceDecomposition {
    /// Atom index ranges for each domain.
    pub domains: Vec<std::ops::Range<usize>>,
    /// Pairs that cross domain boundaries (must be handled separately).
    pub cross_pairs: Vec<(usize, usize)>,
}

#[allow(dead_code)]
impl ForceDecomposition {
    /// Decompose atoms into `n_domains` contiguous slices.
    pub fn from_atoms(n_atoms: usize, n_domains: usize) -> Self {
        let n_d = n_domains.max(1);
        let chunk = n_atoms.div_ceil(n_d);
        let domains: Vec<std::ops::Range<usize>> = (0..n_d)
            .map(|i| {
                let start = i * chunk;
                let end = (start + chunk).min(n_atoms);
                start..end
            })
            .filter(|r| !r.is_empty())
            .collect();
        Self {
            domains,
            cross_pairs: Vec::new(),
        }
    }

    /// Classify pairs from a Verlet list into intra-domain and cross-domain.
    pub fn classify_pairs(&mut self, pairs: &[(usize, usize)]) {
        self.cross_pairs.clear();
        for &(i, j) in pairs {
            let di = self.domain_of(i);
            let dj = self.domain_of(j);
            if di != dj {
                self.cross_pairs.push((i, j));
            }
        }
    }

    /// Return the domain index for atom `i`, or `None` if out of range.
    fn domain_of(&self, i: usize) -> usize {
        for (idx, range) in self.domains.iter().enumerate() {
            if range.contains(&i) {
                return idx;
            }
        }
        usize::MAX
    }

    /// Number of domains.
    pub fn n_domains(&self) -> usize {
        self.domains.len()
    }
}

// ---------------------------------------------------------------------------
// SkinAdaptiveVerletList — adaptive skin sizing
// ---------------------------------------------------------------------------

/// Verlet list with adaptive skin distance.
///
/// Tracks rebuild frequency and adjusts the skin distance to achieve
/// a target rebuild interval.  If the list is rebuilt too often the skin
/// is increased; if it is rebuilt rarely the skin is decreased.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct SkinAdaptiveVerletList {
    /// Underlying Verlet list.
    pub vlist: VerletList,
    /// Current adaptive skin distance.
    pub skin: f64,
    /// Target number of MD steps between rebuilds.
    pub target_steps: usize,
    /// Actual steps since last rebuild.
    pub steps_since_rebuild: usize,
    /// Minimum allowed skin distance.
    pub skin_min: f64,
    /// Maximum allowed skin distance.
    pub skin_max: f64,
}

#[allow(dead_code)]
impl SkinAdaptiveVerletList {
    /// Create a new adaptive Verlet list.
    pub fn new(
        positions: &[Vec3],
        pbox: &PeriodicBox,
        cutoff: f64,
        skin_init: f64,
        target_steps: usize,
    ) -> Self {
        let vlist = VerletList::build(positions, pbox, cutoff, skin_init);
        Self {
            vlist,
            skin: skin_init,
            target_steps,
            steps_since_rebuild: 0,
            skin_min: skin_init * 0.1,
            skin_max: skin_init * 5.0,
        }
    }

    /// Advance one MD step.  Returns `true` if list was rebuilt.
    pub fn advance(&mut self, positions: &[Vec3], pbox: &PeriodicBox) -> bool {
        self.steps_since_rebuild += 1;
        if self.vlist.needs_rebuild(positions, pbox) {
            // Adapt skin based on rebuild frequency
            if self.steps_since_rebuild < self.target_steps / 2 {
                self.skin = (self.skin * 1.2).min(self.skin_max);
            } else if self.steps_since_rebuild > self.target_steps * 2 {
                self.skin = (self.skin * 0.9).max(self.skin_min);
            }
            self.vlist = VerletList::build(positions, pbox, self.vlist.cutoff, self.skin);
            self.steps_since_rebuild = 0;
            true
        } else {
            false
        }
    }
}

// ---------------------------------------------------------------------------
// Parallelization hints
// ---------------------------------------------------------------------------

/// Hints for parallelizing neighbor list operations.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct ParallelHints {
    /// Recommended chunk size for splitting atoms across threads.
    pub chunk_size: usize,
    /// Number of recommended parallel workers.
    pub n_workers: usize,
    /// Whether the problem is large enough to benefit from parallelism.
    pub worth_parallelizing: bool,
}

#[allow(dead_code)]
impl ParallelHints {
    /// Compute parallel hints for a system of `n_atoms` atoms.
    ///
    /// `min_atoms_per_worker` is the threshold below which parallelism is not worthwhile.
    pub fn compute(n_atoms: usize, available_threads: usize, min_atoms_per_worker: usize) -> Self {
        let n_workers = available_threads.max(1);
        let worth_parallelizing = n_atoms >= min_atoms_per_worker * 2;
        let chunk_size = n_atoms
            .checked_div(n_workers)
            .map(|c| c.max(1))
            .unwrap_or(n_atoms);
        Self {
            chunk_size,
            n_workers,
            worth_parallelizing,
        }
    }

    /// Compute atom index ranges for parallel work distribution.
    pub fn work_ranges(&self, n_atoms: usize) -> Vec<std::ops::Range<usize>> {
        let cs = self.chunk_size.max(1);
        (0..n_atoms)
            .step_by(cs)
            .map(|start| start..(start + cs).min(n_atoms))
            .collect()
    }
}

// ---------------------------------------------------------------------------
// All-pairs brute force neighbor finder
// ---------------------------------------------------------------------------

/// Brute-force O(N²) all-pairs distance computation for validation.
///
/// Returns all unique pairs (i, j) with i < j whose PBC distance <= cutoff.
#[allow(dead_code)]
pub fn brute_force_neighbors(
    positions: &[Vec3],
    pbox: &PeriodicBox,
    cutoff: f64,
) -> Vec<(usize, usize)> {
    let n = positions.len();
    let mut pairs = Vec::new();
    for i in 0..n {
        for j in (i + 1)..n {
            let dr = pbox.minimum_image(&(positions[j] - positions[i]));
            if dr.norm() <= cutoff {
                pairs.push((i, j));
            }
        }
    }
    pairs
}

/// Compare cell-list and brute-force neighbor counts for validation.
///
/// Returns `true` if both produce the same number of pairs.
#[allow(dead_code)]
pub fn validate_cell_list(positions: &[Vec3], pbox: &PeriodicBox, cutoff: f64) -> bool {
    let bf = brute_force_neighbors(positions, pbox, cutoff);
    let cl = CellList::build(positions, pbox, cutoff);
    let cl_pairs = cl.neighbor_pairs(positions, pbox);
    let cl_count = cl_pairs.len();
    bf.len() == cl_count
}

// ---------------------------------------------------------------------------
// Exclusion list
// ---------------------------------------------------------------------------

/// A list of bonded-atom pairs to exclude from non-bonded computations.
#[allow(dead_code)]
#[derive(Debug, Clone, Default)]
pub struct ExclusionList {
    /// Excluded pairs stored as (i, j) with i < j.
    pub pairs: std::collections::HashSet<(usize, usize)>,
}

impl ExclusionList {
    /// Create an empty exclusion list.
    #[allow(dead_code)]
    pub fn new() -> Self {
        Self {
            pairs: std::collections::HashSet::new(),
        }
    }

    /// Add a 1-2 exclusion.
    #[allow(dead_code)]
    pub fn add_12(&mut self, i: usize, j: usize) {
        let (a, b) = if i < j { (i, j) } else { (j, i) };
        self.pairs.insert((a, b));
    }

    /// Add 1-3 exclusions from a bond list (pairs of bonded atom lists).
    #[allow(dead_code)]
    pub fn add_13_from_bonds(&mut self, bonds: &[(usize, usize)]) {
        // For every pair of bonds sharing an atom, add the 1-3 pair
        for &(a, b) in bonds {
            for &(c, d) in bonds {
                let excl: Option<(usize, usize)> = if a == c && b != d {
                    Some((b.min(d), b.max(d)))
                } else if a == d && b != c {
                    Some((b.min(c), b.max(c)))
                } else if b == c && a != d {
                    Some((a.min(d), a.max(d)))
                } else if b == d && a != c {
                    Some((a.min(c), a.max(c)))
                } else {
                    None
                };
                if let Some(p) = excl {
                    self.pairs.insert(p);
                }
            }
        }
    }

    /// Check whether a pair (i, j) is excluded.
    #[allow(dead_code)]
    pub fn is_excluded(&self, i: usize, j: usize) -> bool {
        let (a, b) = if i < j { (i, j) } else { (j, i) };
        self.pairs.contains(&(a, b))
    }

    /// Number of exclusions.
    #[allow(dead_code)]
    pub fn len(&self) -> usize {
        self.pairs.len()
    }

    /// Returns true if the exclusion list is empty.
    #[allow(dead_code)]
    pub fn is_empty(&self) -> bool {
        self.pairs.is_empty()
    }

    /// Filter a neighbor pair list, removing excluded pairs.
    #[allow(dead_code)]
    pub fn filter_pairs(&self, pairs: &[(usize, usize)]) -> Vec<(usize, usize)> {
        pairs
            .iter()
            .filter(|&&(i, j)| !self.is_excluded(i, j))
            .copied()
            .collect()
    }
}

// ---------------------------------------------------------------------------
// Half-neighbor list (Newton's third law optimization)
// ---------------------------------------------------------------------------

/// A half-neighbor list that stores only (i, j) pairs with i < j.
///
/// This avoids storing both (i,j) and (j,i), halving storage and allowing
/// Newton's third law to reduce computation by half.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct HalfNeighborList {
    /// Pairs (i, j) with i < j and distance <= cutoff + skin.
    pub pairs: Vec<(usize, usize)>,
    /// Cutoff distance.
    pub cutoff: f64,
    /// Skin distance.
    pub skin: f64,
}

impl HalfNeighborList {
    /// Build a half-neighbor list from positions.
    #[allow(dead_code)]
    pub fn build(positions: &[Vec3], pbox: &PeriodicBox, cutoff: f64, skin: f64) -> Self {
        let r_max = cutoff + skin;
        let r_max2 = r_max * r_max;
        let n = positions.len();
        let mut pairs = Vec::new();
        for i in 0..n {
            for j in (i + 1)..n {
                let dr = pbox.minimum_image(&(positions[j] - positions[i]));
                if dr.norm_squared() <= r_max2 {
                    pairs.push((i, j));
                }
            }
        }
        Self {
            pairs,
            cutoff,
            skin,
        }
    }

    /// Number of pairs.
    #[allow(dead_code)]
    pub fn len(&self) -> usize {
        self.pairs.len()
    }

    /// Returns true if list is empty.
    #[allow(dead_code)]
    pub fn is_empty(&self) -> bool {
        self.pairs.is_empty()
    }

    /// Filter by exclusion list, returning a new half-list.
    #[allow(dead_code)]
    pub fn apply_exclusions(&self, excl: &ExclusionList) -> Self {
        let pairs = excl.filter_pairs(&self.pairs);
        Self {
            pairs,
            cutoff: self.cutoff,
            skin: self.skin,
        }
    }

    /// Count pairs involving atom `idx`.
    #[allow(dead_code)]
    pub fn count_for_atom(&self, idx: usize) -> usize {
        self.pairs
            .iter()
            .filter(|&&(i, j)| i == idx || j == idx)
            .count()
    }
}

// ---------------------------------------------------------------------------
// Neighbor counting utilities
// ---------------------------------------------------------------------------

/// Count the number of neighbors for each atom within `cutoff` using PBC.
///
/// Returns a `Vec`usize` where entry `i` is the count of atoms j (j != i)
/// within cutoff of atom i.
#[allow(dead_code)]
pub fn neighbor_counts(positions: &[Vec3], pbox: &PeriodicBox, cutoff: f64) -> Vec<usize> {
    let n = positions.len();
    let cutoff2 = cutoff * cutoff;
    let mut counts = vec![0usize; n];
    for i in 0..n {
        for j in 0..n {
            if i == j {
                continue;
            }
            let dr = pbox.minimum_image(&(positions[j] - positions[i]));
            if dr.norm_squared() <= cutoff2 {
                counts[i] += 1;
            }
        }
    }
    counts
}

/// Compute the average number of neighbors per atom.
#[allow(dead_code)]
pub fn average_neighbor_count(positions: &[Vec3], pbox: &PeriodicBox, cutoff: f64) -> f64 {
    if positions.is_empty() {
        return 0.0;
    }
    let counts = neighbor_counts(positions, pbox, cutoff);
    counts.iter().sum::<usize>() as f64 / positions.len() as f64
}

/// Find the atom with the most neighbors within `cutoff`.
///
/// Returns `None` if the list is empty.
#[allow(dead_code)]
pub fn most_connected_atom(positions: &[Vec3], pbox: &PeriodicBox, cutoff: f64) -> Option<usize> {
    let counts = neighbor_counts(positions, pbox, cutoff);
    counts
        .iter()
        .enumerate()
        .max_by_key(|(_, c)| *c)
        .map(|(i, _)| i)
}

// ---------------------------------------------------------------------------
// Cell list with image cell iteration (complete neighbor shell)
// ---------------------------------------------------------------------------

/// Utility: given a list of positions, return the distance matrix (upper triangle only).
///
/// Entry `(i, j)` with i < j is the PBC distance between atoms i and j.
#[allow(dead_code)]
pub fn distance_matrix_pbc(positions: &[Vec3], pbox: &PeriodicBox) -> Vec<(usize, usize, f64)> {
    let n = positions.len();
    let mut entries = Vec::with_capacity(n * (n - 1) / 2);
    for i in 0..n {
        for j in (i + 1)..n {
            let (_, dist) = distance_pbc(&positions[i], &positions[j], pbox);
            entries.push((i, j, dist));
        }
    }
    entries
}

/// Radial distribution function (RDF) histogram.
///
/// Returns a histogram of pair distances binned into `n_bins` bins of width
/// `(r_max - r_min) / n_bins`. Bin `k` spans `\[r_min + k*dr, r_min + (k+1)*dr)`.
#[allow(dead_code)]
pub fn rdf_histogram(
    positions: &[Vec3],
    pbox: &PeriodicBox,
    r_min: f64,
    r_max: f64,
    n_bins: usize,
) -> Vec<usize> {
    let mut hist = vec![0usize; n_bins];
    let dr = (r_max - r_min) / n_bins as f64;
    if dr <= 0.0 {
        return hist;
    }
    let n = positions.len();
    for i in 0..n {
        for j in (i + 1)..n {
            let (_, dist) = distance_pbc(&positions[i], &positions[j], pbox);
            if dist >= r_min && dist < r_max {
                let bin = ((dist - r_min) / dr) as usize;
                if bin < n_bins {
                    hist[bin] += 1;
                }
            }
        }
    }
    hist
}

// ---------------------------------------------------------------------------
// Verlet list rebuild statistics
// ---------------------------------------------------------------------------

/// Track Verlet list rebuild statistics over a simulation run.
#[allow(dead_code)]
#[derive(Debug, Clone, Default)]
pub struct RebuildStats {
    /// Total number of rebuild events.
    pub n_rebuilds: usize,
    /// Total simulation steps processed.
    pub n_steps: usize,
    /// Cumulative rebuild cost (arbitrary units, e.g. pair count per rebuild).
    pub total_cost: f64,
}

impl RebuildStats {
    /// Create a new stats tracker.
    #[allow(dead_code)]
    pub fn new() -> Self {
        Self::default()
    }

    /// Record a step; if `rebuilt` is true, also record the rebuild cost.
    #[allow(dead_code)]
    pub fn record_step(&mut self, rebuilt: bool, pair_count: usize) {
        self.n_steps += 1;
        if rebuilt {
            self.n_rebuilds += 1;
            self.total_cost += pair_count as f64;
        }
    }

    /// Average steps between rebuilds.
    #[allow(dead_code)]
    pub fn avg_steps_between_rebuilds(&self) -> f64 {
        if self.n_rebuilds == 0 {
            return self.n_steps as f64;
        }
        self.n_steps as f64 / self.n_rebuilds as f64
    }

    /// Rebuild fraction (fraction of steps that required a rebuild).
    #[allow(dead_code)]
    pub fn rebuild_fraction(&self) -> f64 {
        if self.n_steps == 0 {
            return 0.0;
        }
        self.n_rebuilds as f64 / self.n_steps as f64
    }
}

// ---------------------------------------------------------------------------
// Grid-based nearest-neighbor (single nearest neighbor per atom)
// ---------------------------------------------------------------------------

/// Find the nearest neighbor for each atom using PBC.
///
/// Returns a `Vec<Option`usize`>` where entry `i` is the index of the
/// closest atom to atom `i` (or `None` if there is only one atom).
#[allow(dead_code)]
pub fn nearest_neighbors(positions: &[Vec3], pbox: &PeriodicBox) -> Vec<Option<usize>> {
    let n = positions.len();
    (0..n)
        .map(|i| {
            let mut best_j = None;
            let mut best_d2 = f64::MAX;
            for j in 0..n {
                if i == j {
                    continue;
                }
                let dr = pbox.minimum_image(&(positions[j] - positions[i]));
                let d2 = dr.norm_squared();
                if d2 < best_d2 {
                    best_d2 = d2;
                    best_j = Some(j);
                }
            }
            best_j
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_minimum_image() {
        let pbox = PeriodicBox::cubic(10.0);
        // Two atoms at x=1 and x=9 => distance should be 2, not 8
        let a = Vec3::new(1.0, 0.0, 0.0);
        let b = Vec3::new(9.0, 0.0, 0.0);
        let (dr, dist) = distance_pbc(&a, &b, &pbox);
        assert!((dist - 2.0).abs() < 1e-10);
        assert!((dr.x - (-2.0)).abs() < 1e-10);
    }

    #[test]
    fn test_wrap_position() {
        let pbox = PeriodicBox::cubic(10.0);
        let pos = Vec3::new(-1.0, 11.0, 5.0);
        let wrapped = pbox.wrap_position(&pos);
        assert!((wrapped.x - 9.0).abs() < 1e-10);
        assert!((wrapped.y - 1.0).abs() < 1e-10);
        assert!((wrapped.z - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_cell_list_finds_neighbors() {
        let pbox = PeriodicBox::cubic(10.0);
        let positions = vec![
            Vec3::new(1.0, 1.0, 1.0),
            Vec3::new(1.5, 1.0, 1.0), // 0.5 away from atom 0
            Vec3::new(5.0, 5.0, 5.0), // far away
        ];
        let cell_list = CellList::build(&positions, &pbox, 2.0);
        let pairs = cell_list.neighbor_pairs(&positions, &pbox);
        // Atoms 0 and 1 should be neighbors
        let has_01 = pairs
            .iter()
            .any(|&(i, j, _, _)| (i == 0 && j == 1) || (i == 1 && j == 0));
        assert!(has_01, "Atoms 0 and 1 should be neighbors");
        // Atom 2 should not be neighbor of 0 or 1 with cutoff=2
        let has_02 = pairs
            .iter()
            .any(|&(i, j, _, _)| (i == 0 && j == 2) || (i == 2 && j == 0));
        assert!(!has_02, "Atoms 0 and 2 should NOT be neighbors");
    }

    #[test]
    fn test_verlet_list_rebuild() {
        let pbox = PeriodicBox::cubic(10.0);
        let positions = vec![Vec3::new(1.0, 1.0, 1.0), Vec3::new(2.0, 1.0, 1.0)];
        let vlist = VerletList::build(&positions, &pbox, 2.5, 0.5);
        assert!(!vlist.needs_rebuild(&positions, &pbox));

        // Move atom 0 by more than half-skin (0.25+)
        let moved = vec![Vec3::new(1.3, 1.0, 1.0), Vec3::new(2.0, 1.0, 1.0)];
        assert!(vlist.needs_rebuild(&moved, &pbox));
    }

    // -----------------------------------------------------------------------
    // Required tests
    // -----------------------------------------------------------------------

    /// 4 atoms in a small box -- every atom is within the cutoff of every other.
    /// The Verlet list must contain all 6 pairs.
    #[test]
    fn test_verlet_list_build() {
        let pbox = PeriodicBox::cubic(10.0);
        // Four atoms clustered near the origin, all within 2.0 of each other.
        let positions = vec![
            Vec3::new(1.0, 1.0, 1.0),
            Vec3::new(1.5, 1.0, 1.0),
            Vec3::new(1.0, 1.5, 1.0),
            Vec3::new(1.0, 1.0, 1.5),
        ];
        let cutoff = 2.0;
        let skin = 0.5;
        let vlist = VerletList::build(&positions, &pbox, cutoff, skin);

        // Verify all 6 pairs (i < j, i in 0..3, j in i+1..4) are present.
        for i in 0..4 {
            for j in (i + 1)..4 {
                let found = vlist.pairs.iter().any(|&(a, b)| a == i && b == j);
                assert!(found, "Verlet list missing pair ({i}, {j})");
            }
        }
    }

    /// Moving an atom by more than skin/2 must trigger needs_rebuild().
    #[test]
    fn test_verlet_list_rebuild_triggered() {
        let pbox = PeriodicBox::cubic(10.0);
        let positions = vec![
            Vec3::new(1.0, 1.0, 1.0),
            Vec3::new(3.0, 1.0, 1.0),
            Vec3::new(5.0, 5.0, 5.0),
        ];
        let skin = 0.4;
        let vlist = VerletList::build(&positions, &pbox, 2.0, skin);

        // Displacement > skin/2 = 0.2
        let mut moved = positions.clone();
        moved[2] = Vec3::new(5.0 + skin * 0.6, 5.0, 5.0);
        assert!(
            vlist.needs_rebuild(&moved, &pbox),
            "needs_rebuild should be true after atom moved > skin/2"
        );
    }

    /// Cell list pairs must exactly match the O(n^2) brute-force pair list.
    #[test]
    fn test_cell_list_pairs_agree_brute_force() {
        let pbox = PeriodicBox::cubic(10.0);
        let positions = vec![
            Vec3::new(0.5, 0.5, 0.5),
            Vec3::new(1.2, 0.5, 0.5),
            Vec3::new(4.0, 4.0, 4.0),
            Vec3::new(4.3, 4.0, 4.0),
            Vec3::new(9.8, 0.5, 0.5), // near periodic boundary of atom 0
        ];
        let cutoff = 1.5;

        // Cell list pairs
        let cl = CellList::build(&positions, &pbox, cutoff);
        let mut cl_pairs: Vec<(usize, usize)> = cl
            .neighbor_pairs(&positions, &pbox)
            .into_iter()
            .map(|(i, j, _, _)| (i, j))
            .collect();
        cl_pairs.sort_unstable();

        // Brute force pairs
        let n = positions.len();
        let mut bf_pairs: Vec<(usize, usize)> = Vec::new();
        for i in 0..n {
            for j in (i + 1)..n {
                let (_, dist) = distance_pbc(&positions[i], &positions[j], &pbox);
                if dist < cutoff {
                    bf_pairs.push((i, j));
                }
            }
        }
        bf_pairs.sort_unstable();

        assert_eq!(
            cl_pairs, bf_pairs,
            "Cell list pairs differ from brute force"
        );
    }

    /// `build_verlet_list` convenience function returns a valid list.
    #[test]
    fn test_build_verlet_list_free_function() {
        let pbox = PeriodicBox::cubic(10.0);
        let positions = vec![
            Vec3::new(1.0, 1.0, 1.0),
            Vec3::new(1.5, 1.0, 1.0),
            Vec3::new(1.0, 1.5, 1.0),
            Vec3::new(5.0, 5.0, 5.0),
        ];
        let cutoff = 1.5;
        let skin = 0.3;

        let vlist = build_verlet_list(&positions, &pbox, cutoff, skin);

        // Atoms 0, 1, 2 are all within 0.71 of each other -- all 3 pairs should
        // be in the list (cutoff + skin = 1.8).
        let has_01 = vlist
            .pairs
            .iter()
            .any(|&(a, b)| (a == 0 && b == 1) || (a == 1 && b == 0));
        let has_02 = vlist
            .pairs
            .iter()
            .any(|&(a, b)| (a == 0 && b == 2) || (a == 2 && b == 0));
        let has_12 = vlist
            .pairs
            .iter()
            .any(|&(a, b)| (a == 1 && b == 2) || (a == 2 && b == 1));
        assert!(has_01, "build_verlet_list: pair (0,1) missing");
        assert!(has_02, "build_verlet_list: pair (0,2) missing");
        assert!(has_12, "build_verlet_list: pair (1,2) missing");

        // Atom 3 is far away and should not appear with 0,1,2
        let has_03 = vlist
            .pairs
            .iter()
            .any(|&(a, b)| (a == 0 && b == 3) || (a == 3 && b == 0));
        assert!(!has_03, "build_verlet_list: pair (0,3) should not exist");
    }

    /// When atoms do not move, needs_rebuild() must return false.
    #[test]
    fn test_verlet_no_rebuild_if_atoms_stationary() {
        let pbox = PeriodicBox::cubic(10.0);
        let positions = vec![
            Vec3::new(2.0, 2.0, 2.0),
            Vec3::new(4.0, 2.0, 2.0),
            Vec3::new(6.0, 6.0, 6.0),
        ];
        let vlist = VerletList::build(&positions, &pbox, 2.5, 0.5);
        assert!(
            !vlist.needs_rebuild(&positions, &pbox),
            "needs_rebuild should be false when atoms have not moved"
        );
    }

    // ── New tests for expanded features ─────────────────────────────────────

    #[test]
    fn test_cell_list_statistics() {
        let pbox = PeriodicBox::cubic(10.0);
        let positions = vec![
            Vec3::new(1.0, 1.0, 1.0),
            Vec3::new(1.5, 1.0, 1.0),
            Vec3::new(5.0, 5.0, 5.0),
        ];
        let cl = CellList::build(&positions, &pbox, 2.0);
        assert!(cl.total_cells() > 0);
        assert!(cl.max_cell_occupancy() >= 1);
        assert!(cl.avg_cell_occupancy() > 0.0);
        assert!(cl.occupied_cells() >= 1);
    }

    #[test]
    fn test_neighbor_list_stats_from_verlet() {
        let pbox = PeriodicBox::cubic(10.0);
        let positions = vec![
            Vec3::new(1.0, 1.0, 1.0),
            Vec3::new(1.5, 1.0, 1.0),
            Vec3::new(1.0, 1.5, 1.0),
            Vec3::new(5.0, 5.0, 5.0),
        ];
        let vlist = VerletList::build(&positions, &pbox, 2.0, 0.5);
        let stats = NeighborListStats::from_verlet_list(&vlist, positions.len());
        assert_eq!(stats.n_atoms, 4);
        // Atoms 0,1,2 should be neighbors of each other -> 3 pairs
        assert!(stats.n_pairs >= 3);
        assert!(stats.max_neighbors >= 2);
        assert!(stats.avg_neighbors > 0.0);
    }

    #[test]
    fn test_neighbor_list_stats_from_cell_list() {
        let pbox = PeriodicBox::cubic(10.0);
        let positions = vec![Vec3::new(1.0, 1.0, 1.0), Vec3::new(1.5, 1.0, 1.0)];
        let cl = CellList::build(&positions, &pbox, 2.0);
        let stats = NeighborListStats::from_cell_list(&cl, &positions, &pbox);
        assert_eq!(stats.n_atoms, 2);
        assert_eq!(stats.n_pairs, 1);
        assert_eq!(stats.max_neighbors, 1);
    }

    #[test]
    fn test_verlet_rebuild_method() {
        let pbox = PeriodicBox::cubic(10.0);
        let positions = vec![Vec3::new(1.0, 1.0, 1.0), Vec3::new(2.0, 1.0, 1.0)];
        let mut vlist = VerletList::build(&positions, &pbox, 2.5, 0.5);
        assert_eq!(vlist.rebuild_count, 0);
        vlist.rebuild(&positions, &pbox);
        assert_eq!(vlist.rebuild_count, 1);
    }

    #[test]
    fn test_verlet_update_if_needed() {
        let pbox = PeriodicBox::cubic(10.0);
        let positions = vec![Vec3::new(1.0, 1.0, 1.0), Vec3::new(2.0, 1.0, 1.0)];
        let mut vlist = VerletList::build(&positions, &pbox, 2.5, 0.5);
        // No movement -> no rebuild
        let rebuilt = vlist.update_if_needed(&positions, &pbox);
        assert!(!rebuilt);
        // Move atom past half-skin
        let moved = vec![Vec3::new(1.4, 1.0, 1.0), Vec3::new(2.0, 1.0, 1.0)];
        let rebuilt = vlist.update_if_needed(&moved, &pbox);
        assert!(rebuilt);
        assert_eq!(vlist.rebuild_count, 1);
    }

    #[test]
    fn test_verlet_max_displacement() {
        let pbox = PeriodicBox::cubic(10.0);
        let positions = vec![Vec3::new(1.0, 1.0, 1.0), Vec3::new(2.0, 1.0, 1.0)];
        let vlist = VerletList::build(&positions, &pbox, 2.5, 0.5);
        let moved = vec![
            Vec3::new(1.0, 1.0, 1.0),
            Vec3::new(2.1, 1.0, 1.0), // moved 0.1
        ];
        let max_d = vlist.max_displacement(&moved, &pbox);
        assert!((max_d - 0.1).abs() < 1e-10);
    }

    #[test]
    fn test_verlet_pair_count() {
        let pbox = PeriodicBox::cubic(10.0);
        let positions = vec![
            Vec3::new(1.0, 1.0, 1.0),
            Vec3::new(1.5, 1.0, 1.0),
            Vec3::new(5.0, 5.0, 5.0),
        ];
        let vlist = VerletList::build(&positions, &pbox, 2.0, 0.5);
        assert_eq!(vlist.pair_count(), 1); // only 0-1 pair
    }

    #[test]
    fn test_multi_cutoff_neighbor_list() {
        let pbox = PeriodicBox::cubic(10.0);
        let positions = vec![
            Vec3::new(1.0, 1.0, 1.0),
            Vec3::new(1.3, 1.0, 1.0), // 0.3 away
            Vec3::new(2.0, 1.0, 1.0), // 1.0 away from atom 0
            Vec3::new(5.0, 5.0, 5.0), // far away
        ];
        let cutoffs = [0.5, 1.5];
        let skin = 0.1;
        let mcnl = MultiCutoffNeighborList::build(&positions, &pbox, &cutoffs, skin);

        // 0.5 + 0.1 = 0.6 cutoff: atoms 0-1 (dist=0.3) should be in list
        assert!(
            mcnl.pair_count(0) >= 1,
            "short cutoff should have at least 1 pair"
        );

        // 1.5 + 0.1 = 1.6 cutoff: atoms 0-1 (0.3) and 0-2 (1.0) and 1-2 (0.7) should be in list
        assert!(
            mcnl.pair_count(1) >= 2,
            "long cutoff should have at least 2 pairs"
        );
    }

    #[test]
    fn test_parallel_hints() {
        let hints = ParallelHints::compute(1000, 4, 100);
        assert_eq!(hints.n_workers, 4);
        assert!(hints.worth_parallelizing);
        assert!(hints.chunk_size >= 1);

        let ranges = hints.work_ranges(1000);
        assert!(!ranges.is_empty());
        // All ranges should cover 0..1000
        let mut covered = vec![false; 1000];
        for range in &ranges {
            for i in range.clone() {
                covered[i] = true;
            }
        }
        assert!(covered.iter().all(|&c| c));
    }

    #[test]
    fn test_parallel_hints_small_system() {
        let hints = ParallelHints::compute(10, 4, 100);
        assert!(!hints.worth_parallelizing);
    }

    #[test]
    fn test_number_density() {
        let pbox = PeriodicBox::cubic(10.0);
        let rho = pbox.number_density(1000);
        assert!((rho - 1.0).abs() < 1e-10); // 1000 / 1000 = 1.0
    }

    #[test]
    fn test_periodic_box_volume() {
        let pbox = PeriodicBox::new(Vec3::new(2.0, 3.0, 4.0));
        assert!((pbox.volume() - 24.0).abs() < 1e-10);
    }

    // ── ForceDecomposition tests ─────────────────────────────────────────────

    #[test]
    fn test_force_decomposition_n_domains() {
        let fd = ForceDecomposition::from_atoms(100, 4);
        assert_eq!(fd.n_domains(), 4);
    }

    #[test]
    fn test_force_decomposition_covers_all_atoms() {
        let n = 97; // non-power-of-2
        let fd = ForceDecomposition::from_atoms(n, 4);
        let mut covered = vec![false; n];
        for range in &fd.domains {
            for i in range.clone() {
                covered[i] = true;
            }
        }
        assert!(covered.iter().all(|&c| c), "all atoms must be covered");
    }

    #[test]
    fn test_force_decomposition_single_domain() {
        let fd = ForceDecomposition::from_atoms(50, 1);
        assert_eq!(fd.n_domains(), 1);
        assert_eq!(fd.domains[0], 0..50);
    }

    #[test]
    fn test_force_decomposition_cross_pairs_initially_empty() {
        let fd = ForceDecomposition::from_atoms(10, 2);
        assert!(fd.cross_pairs.is_empty());
    }

    #[test]
    fn test_force_decomposition_classify_intra_pairs() {
        let mut fd = ForceDecomposition::from_atoms(10, 2);
        // Both atoms in domain 0 (0..5)
        let pairs = vec![(0usize, 2usize)];
        fd.classify_pairs(&pairs);
        assert!(
            fd.cross_pairs.is_empty(),
            "intra-domain pair should not appear in cross_pairs"
        );
    }

    #[test]
    fn test_force_decomposition_classify_cross_pair() {
        let mut fd = ForceDecomposition::from_atoms(10, 2);
        // Atoms 0 (domain 0) and 6 (domain 1) → cross-pair
        let pairs = vec![(0usize, 6usize)];
        fd.classify_pairs(&pairs);
        assert_eq!(
            fd.cross_pairs.len(),
            1,
            "cross-domain pair should be classified"
        );
    }

    // ── SkinAdaptiveVerletList tests ──────────────────────────────────────────

    #[test]
    fn test_skin_adaptive_initial_state() {
        let pbox = PeriodicBox::cubic(10.0);
        let positions = vec![Vec3::new(1.0, 1.0, 1.0), Vec3::new(2.0, 1.0, 1.0)];
        let avl = SkinAdaptiveVerletList::new(&positions, &pbox, 2.5, 0.5, 10);
        assert_eq!(avl.steps_since_rebuild, 0);
        assert!((avl.skin - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_skin_adaptive_no_rebuild_when_stationary() {
        let pbox = PeriodicBox::cubic(10.0);
        let positions = vec![Vec3::new(1.0, 1.0, 1.0), Vec3::new(2.0, 1.0, 1.0)];
        let mut avl = SkinAdaptiveVerletList::new(&positions, &pbox, 2.5, 0.5, 10);
        let rebuilt = avl.advance(&positions, &pbox);
        assert!(!rebuilt, "stationary atoms should not trigger rebuild");
        assert_eq!(avl.steps_since_rebuild, 1);
    }

    #[test]
    fn test_skin_adaptive_rebuild_when_moved() {
        let pbox = PeriodicBox::cubic(10.0);
        let positions = vec![Vec3::new(1.0, 1.0, 1.0), Vec3::new(2.0, 1.0, 1.0)];
        let mut avl = SkinAdaptiveVerletList::new(&positions, &pbox, 2.5, 0.4, 10);
        // Move atom 0 far enough
        let moved = vec![
            Vec3::new(1.25, 1.0, 1.0), // moved 0.25 > skin/2 = 0.2
            Vec3::new(2.0, 1.0, 1.0),
        ];
        let rebuilt = avl.advance(&moved, &pbox);
        assert!(rebuilt, "moving atom > skin/2 should trigger rebuild");
    }

    // ── PBC cell list with periodic boundary tests ────────────────────────────

    #[test]
    fn test_cell_list_pbc_periodic_image_pair() {
        // Two atoms near opposite ends of box — they should be neighbors under PBC
        let pbox = PeriodicBox::cubic(10.0);
        let positions = vec![
            Vec3::new(0.3, 5.0, 5.0),
            Vec3::new(9.8, 5.0, 5.0), // distance = 0.5 under PBC
        ];
        let cl = CellList::build(&positions, &pbox, 1.0);
        let pairs = cl.neighbor_pairs(&positions, &pbox);
        let found = pairs
            .iter()
            .any(|&(i, j, _, _)| (i == 0 && j == 1) || (i == 1 && j == 0));
        assert!(found, "periodic image pair (0,1) not found in cell list");
    }

    #[test]
    fn test_verlet_list_no_self_pairs() {
        let pbox = PeriodicBox::cubic(10.0);
        let positions = vec![Vec3::new(1.0, 1.0, 1.0), Vec3::new(1.5, 1.0, 1.0)];
        let vlist = VerletList::build(&positions, &pbox, 2.0, 0.5);
        for &(i, j) in &vlist.pairs {
            assert_ne!(i, j, "Verlet list should not contain self-pairs");
        }
    }

    #[test]
    fn test_cell_list_no_self_pairs() {
        let pbox = PeriodicBox::cubic(10.0);
        let positions = vec![
            Vec3::new(1.0, 1.0, 1.0),
            Vec3::new(1.5, 1.0, 1.0),
            Vec3::new(5.0, 5.0, 5.0),
        ];
        let cl = CellList::build(&positions, &pbox, 1.5);
        let pairs = cl.neighbor_pairs(&positions, &pbox);
        for (i, j, _, _) in pairs {
            assert_ne!(i, j, "Cell list should not contain self-pairs");
        }
    }

    // ── Brute-force neighbor tests ─────────────────────────────────────────

    #[test]
    fn test_brute_force_neighbors_basic() {
        let pbox = PeriodicBox::cubic(10.0);
        let positions = vec![
            Vec3::new(1.0, 1.0, 1.0),
            Vec3::new(1.3, 1.0, 1.0), // 0.3 away
            Vec3::new(5.0, 5.0, 5.0), // far
        ];
        let pairs = brute_force_neighbors(&positions, &pbox, 0.5);
        assert_eq!(pairs.len(), 1, "only (0,1) should be within 0.5");
        assert_eq!(pairs[0], (0, 1));
    }

    #[test]
    fn test_brute_force_neighbors_pbc() {
        let pbox = PeriodicBox::cubic(10.0);
        let positions = vec![
            Vec3::new(0.2, 5.0, 5.0),
            Vec3::new(9.9, 5.0, 5.0), // 0.3 under PBC
        ];
        let pairs = brute_force_neighbors(&positions, &pbox, 0.5);
        assert_eq!(pairs.len(), 1, "PBC pair should be found");
    }

    #[test]
    fn test_brute_force_no_pairs() {
        let pbox = PeriodicBox::cubic(10.0);
        let positions = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(5.0, 5.0, 5.0)];
        let pairs = brute_force_neighbors(&positions, &pbox, 0.1);
        assert!(pairs.is_empty(), "no pairs within 0.1 nm");
    }

    #[test]
    fn test_validate_cell_list() {
        let pbox = PeriodicBox::cubic(10.0);
        let positions = vec![
            Vec3::new(1.0, 1.0, 1.0),
            Vec3::new(1.5, 1.0, 1.0),
            Vec3::new(5.0, 5.0, 5.0),
            Vec3::new(5.4, 5.0, 5.0),
        ];
        assert!(
            validate_cell_list(&positions, &pbox, 1.0),
            "cell list should match brute force"
        );
    }

    // ── ExclusionList tests ───────────────────────────────────────────────

    #[test]
    fn test_exclusion_list_add_and_check() {
        let mut excl = ExclusionList::new();
        excl.add_12(0, 1);
        assert!(excl.is_excluded(0, 1));
        assert!(excl.is_excluded(1, 0), "symmetric");
        assert!(!excl.is_excluded(0, 2));
    }

    #[test]
    fn test_exclusion_list_len() {
        let mut excl = ExclusionList::new();
        excl.add_12(0, 1);
        excl.add_12(1, 2);
        assert_eq!(excl.len(), 2);
    }

    #[test]
    fn test_exclusion_list_filter_pairs() {
        let mut excl = ExclusionList::new();
        excl.add_12(0, 1);
        let pairs = vec![(0, 1), (0, 2), (1, 2)];
        let filtered = excl.filter_pairs(&pairs);
        assert_eq!(filtered.len(), 2);
        assert!(!filtered.contains(&(0, 1)));
    }

    #[test]
    fn test_exclusion_list_empty() {
        let excl = ExclusionList::new();
        assert!(excl.is_empty());
    }

    #[test]
    fn test_exclusion_list_13_from_bonds() {
        let mut excl = ExclusionList::new();
        let bonds = vec![(0, 1), (1, 2)]; // chain 0-1-2
        excl.add_13_from_bonds(&bonds);
        // 1-3 pair (0,2) should be excluded
        assert!(excl.is_excluded(0, 2), "1-3 pair (0,2) should be excluded");
    }

    // ── HalfNeighborList tests ────────────────────────────────────────────

    #[test]
    fn test_half_neighbor_list_basic() {
        let pbox = PeriodicBox::cubic(10.0);
        let positions = vec![
            Vec3::new(1.0, 1.0, 1.0),
            Vec3::new(1.5, 1.0, 1.0), // 0.5 away
            Vec3::new(5.0, 5.0, 5.0), // far
        ];
        let hnl = HalfNeighborList::build(&positions, &pbox, 1.0, 0.0);
        assert_eq!(hnl.len(), 1, "only pair (0,1) within 1.0 nm");
        assert_eq!(hnl.pairs[0], (0, 1));
    }

    #[test]
    fn test_half_neighbor_list_unique_pairs() {
        let pbox = PeriodicBox::cubic(10.0);
        let positions = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.3, 0.0, 0.0),
            Vec3::new(0.6, 0.0, 0.0),
        ];
        let hnl = HalfNeighborList::build(&positions, &pbox, 1.0, 0.0);
        // Check all pairs are i < j
        for &(i, j) in &hnl.pairs {
            assert!(i < j, "half list must have i < j");
        }
    }

    #[test]
    fn test_half_neighbor_list_count_for_atom() {
        let pbox = PeriodicBox::cubic(10.0);
        let positions = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.3, 0.0, 0.0),
            Vec3::new(0.6, 0.0, 0.0),
        ];
        let hnl = HalfNeighborList::build(&positions, &pbox, 1.0, 0.0);
        let cnt = hnl.count_for_atom(1);
        assert_eq!(cnt, 2, "atom 1 should appear in 2 pairs (0-1 and 1-2)");
    }

    #[test]
    fn test_half_neighbor_list_apply_exclusions() {
        let pbox = PeriodicBox::cubic(10.0);
        let positions = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.3, 0.0, 0.0),
            Vec3::new(0.6, 0.0, 0.0),
        ];
        let hnl = HalfNeighborList::build(&positions, &pbox, 1.0, 0.0);
        let mut excl = ExclusionList::new();
        excl.add_12(0, 1);
        let filtered = hnl.apply_exclusions(&excl);
        assert!(!filtered.pairs.contains(&(0, 1)), "excluded pair removed");
    }

    // ── Neighbor counting tests ───────────────────────────────────────────

    #[test]
    fn test_neighbor_counts_basic() {
        let pbox = PeriodicBox::cubic(10.0);
        let positions = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.3, 0.0, 0.0),
            Vec3::new(5.0, 0.0, 0.0),
        ];
        let counts = neighbor_counts(&positions, &pbox, 0.5);
        assert_eq!(counts[0], 1, "atom 0 has 1 neighbor");
        assert_eq!(counts[1], 1, "atom 1 has 1 neighbor");
        assert_eq!(counts[2], 0, "atom 2 is isolated");
    }

    #[test]
    fn test_average_neighbor_count() {
        let pbox = PeriodicBox::cubic(10.0);
        let positions = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.3, 0.0, 0.0),
            Vec3::new(5.0, 0.0, 0.0),
        ];
        let avg = average_neighbor_count(&positions, &pbox, 0.5);
        assert!((avg - 2.0 / 3.0).abs() < 1e-10, "avg neighbors: {avg}");
    }

    #[test]
    fn test_most_connected_atom() {
        let pbox = PeriodicBox::cubic(10.0);
        let positions = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.2, 0.0, 0.0), // close to atom 0
            Vec3::new(0.1, 0.0, 0.0), // very close to both 0 and 1
        ];
        let idx = most_connected_atom(&positions, &pbox, 0.5);
        assert!(idx.is_some());
    }

    #[test]
    fn test_most_connected_atom_empty() {
        let pbox = PeriodicBox::cubic(10.0);
        let positions: Vec<Vec3> = vec![];
        let idx = most_connected_atom(&positions, &pbox, 1.0);
        assert!(idx.is_none());
    }

    // ── Distance matrix tests ──────────────────────────────────────────────

    #[test]
    fn test_distance_matrix_pbc_basic() {
        let pbox = PeriodicBox::cubic(10.0);
        let positions = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(3.0, 4.0, 0.0), // distance = 5
        ];
        let dm = distance_matrix_pbc(&positions, &pbox);
        assert_eq!(dm.len(), 1);
        let (i, j, d) = dm[0];
        assert_eq!(i, 0);
        assert_eq!(j, 1);
        assert!((d - 5.0).abs() < 1e-10, "distance = {d}");
    }

    #[test]
    fn test_distance_matrix_pbc_n_entries() {
        let pbox = PeriodicBox::cubic(10.0);
        let positions: Vec<Vec3> = (0..5).map(|i| Vec3::new(i as f64, 0.0, 0.0)).collect();
        let dm = distance_matrix_pbc(&positions, &pbox);
        assert_eq!(dm.len(), 5 * 4 / 2, "n*(n-1)/2 entries expected");
    }

    // ── RDF histogram tests ────────────────────────────────────────────────

    #[test]
    fn test_rdf_histogram_one_pair() {
        let pbox = PeriodicBox::cubic(10.0);
        let positions = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.5, 0.0, 0.0)];
        let hist = rdf_histogram(&positions, &pbox, 0.0, 5.0, 10);
        // Distance 1.5 falls in bin 3 (1.5 / 0.5 = 3)
        let total: usize = hist.iter().sum();
        assert_eq!(total, 1, "one pair should contribute 1 count");
    }

    #[test]
    fn test_rdf_histogram_empty() {
        let pbox = PeriodicBox::cubic(10.0);
        let positions: Vec<Vec3> = vec![];
        let hist = rdf_histogram(&positions, &pbox, 0.0, 5.0, 10);
        assert!(hist.iter().all(|&c| c == 0));
    }

    #[test]
    fn test_rdf_histogram_bin_count() {
        let pbox = PeriodicBox::cubic(10.0);
        let positions = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(2.0, 0.0, 0.0)];
        let hist = rdf_histogram(&positions, &pbox, 0.0, 5.0, 20);
        assert_eq!(hist.len(), 20);
    }

    // ── RebuildStats tests ────────────────────────────────────────────────

    #[test]
    fn test_rebuild_stats_fraction() {
        let mut stats = RebuildStats::new();
        stats.record_step(true, 100);
        stats.record_step(false, 0);
        stats.record_step(false, 0);
        stats.record_step(false, 0);
        assert_eq!(stats.n_rebuilds, 1);
        assert_eq!(stats.n_steps, 4);
        assert!((stats.rebuild_fraction() - 0.25).abs() < 1e-10);
    }

    #[test]
    fn test_rebuild_stats_avg_steps() {
        let mut stats = RebuildStats::new();
        for i in 0..10 {
            stats.record_step(i % 5 == 0, 50);
        }
        assert!(stats.avg_steps_between_rebuilds() > 0.0);
    }

    #[test]
    fn test_rebuild_stats_no_rebuilds() {
        let mut stats = RebuildStats::new();
        stats.record_step(false, 0);
        stats.record_step(false, 0);
        assert_eq!(stats.rebuild_fraction(), 0.0);
        // avg_steps = total steps when no rebuilds
        assert!((stats.avg_steps_between_rebuilds() - 2.0).abs() < 1e-10);
    }

    // ── Nearest-neighbor tests ────────────────────────────────────────────

    #[test]
    fn test_nearest_neighbors_basic() {
        let pbox = PeriodicBox::cubic(10.0);
        let positions = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.5, 0.0, 0.0),
            Vec3::new(9.0, 0.0, 0.0),
        ];
        let nn = nearest_neighbors(&positions, &pbox);
        assert_eq!(nn.len(), 3);
        // Atom 0's nearest neighbor should be atom 1 (distance 0.5) not atom 2
        assert_eq!(nn[0], Some(1), "atom 0 nearest = 1");
    }

    #[test]
    fn test_nearest_neighbors_pbc() {
        let pbox = PeriodicBox::cubic(10.0);
        let positions = vec![
            Vec3::new(0.1, 5.0, 5.0),
            Vec3::new(9.9, 5.0, 5.0), // 0.2 under PBC, 9.8 direct
        ];
        let nn = nearest_neighbors(&positions, &pbox);
        assert_eq!(nn[0], Some(1));
        assert_eq!(nn[1], Some(0));
    }

    #[test]
    fn test_nearest_neighbors_single_atom() {
        let pbox = PeriodicBox::cubic(10.0);
        let positions = vec![Vec3::new(1.0, 1.0, 1.0)];
        let nn = nearest_neighbors(&positions, &pbox);
        assert_eq!(nn[0], None);
    }
}
