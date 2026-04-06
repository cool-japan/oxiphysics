//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
#[allow(unused_imports)]
use super::functions::*;
use rayon::prelude::*;

/// Resizes a grid cell list when particles move outside the current bounding box.
///
/// Tracks the axis-aligned bounding box of all particles and rebuilds
/// the `GpuCellList` with padded bounds whenever particles escape.
pub struct GridResizer {
    /// Target cell size (edge length of each cubic cell).
    pub cell_size: f64,
    /// Padding added to each face of the bounding box on resize.
    pub padding: f64,
    /// The most recently built cell list.
    pub cell_list: Option<GpuCellList>,
    /// Bounding box lower corner used when building the current list.
    pub current_min: [f64; 3],
    /// Bounding box upper corner used when building the current list.
    pub current_max: [f64; 3],
}
impl GridResizer {
    /// Create a new `GridResizer` with the given cell size and padding.
    pub fn new(cell_size: f64, padding: f64) -> Self {
        Self {
            cell_size,
            padding,
            cell_list: None,
            current_min: [0.0; 3],
            current_max: [0.0; 3],
        }
    }
    /// Check whether `positions` fit inside the current bounding box.
    pub fn needs_resize(&self, positions: &[[f64; 3]]) -> bool {
        if self.cell_list.is_none() {
            return true;
        }
        for p in positions {
            for d in 0..3 {
                if p[d] < self.current_min[d] || p[d] >= self.current_max[d] {
                    return true;
                }
            }
        }
        false
    }
    /// Rebuild the cell list from `positions`, padding the bounding box.
    pub fn rebuild(&mut self, positions: &[[f64; 3]]) {
        if positions.is_empty() {
            let n = (self.padding / self.cell_size).ceil() as usize;
            let n = n.max(1);
            let box_len = n as f64 * self.cell_size;
            let cl = GpuCellList::new([n, n, n], self.cell_size, [box_len; 3]);
            self.current_min = [0.0; 3];
            self.current_max = [box_len; 3];
            self.cell_list = Some(cl);
            return;
        }
        let (raw_min, raw_max) = compute_bounding_box(positions);
        let mut padded_min = [0.0_f64; 3];
        let mut padded_max = [0.0_f64; 3];
        for d in 0..3 {
            padded_min[d] = raw_min[d] - self.padding;
            padded_max[d] = raw_max[d] + self.padding;
        }
        let nx = (((padded_max[0] - padded_min[0]) / self.cell_size).ceil() as usize).max(1);
        let ny = (((padded_max[1] - padded_min[1]) / self.cell_size).ceil() as usize).max(1);
        let nz = (((padded_max[2] - padded_min[2]) / self.cell_size).ceil() as usize).max(1);
        let box_lengths = [
            nx as f64 * self.cell_size,
            ny as f64 * self.cell_size,
            nz as f64 * self.cell_size,
        ];
        let shifted: Vec<[f64; 3]> = positions
            .iter()
            .map(|p| {
                [
                    p[0] - padded_min[0],
                    p[1] - padded_min[1],
                    p[2] - padded_min[2],
                ]
            })
            .collect();
        let cl = GpuCellList::build(&shifted, [nx, ny, nz], self.cell_size, box_lengths);
        self.current_min = padded_min;
        self.current_max = [
            padded_min[0] + box_lengths[0],
            padded_min[1] + box_lengths[1],
            padded_min[2] + box_lengths[2],
        ];
        self.cell_list = Some(cl);
    }
    /// Ensure the cell list is up-to-date for `positions`.
    ///
    /// Rebuilds only if necessary (particles escaped the current bounds).
    pub fn update(&mut self, positions: &[[f64; 3]]) {
        if self.needs_resize(positions) {
            self.rebuild(positions);
        }
    }
    /// Return a reference to the current cell list, if any.
    pub fn get(&self) -> Option<&GpuCellList> {
        self.cell_list.as_ref()
    }
}
/// Manages ghost (image) particles for periodic boundary conditions.
///
/// In a periodic simulation box, particles near one face need to interact with
/// particles near the opposite face.  Ghost particles are copies of real
/// particles translated by ±box_length into the "halo" region outside the
/// primary cell.
#[derive(Debug, Clone)]
pub struct GhostCellManager {
    /// Primary simulation box lengths \[Lx, Ly, Lz\].
    pub box_lengths: [f64; 3],
    /// Ghost-layer width (typically the pair cutoff).
    pub ghost_width: f64,
    /// Positions of ghost particles (may be outside \[0, L) range).
    pub ghost_positions: Vec<[f64; 3]>,
    /// Map from ghost index to the original (real) particle index.
    pub ghost_to_real: Vec<usize>,
}
impl GhostCellManager {
    /// Create a new manager for the given box and ghost-layer width.
    pub fn new(box_lengths: [f64; 3], ghost_width: f64) -> Self {
        Self {
            box_lengths,
            ghost_width,
            ghost_positions: Vec::new(),
            ghost_to_real: Vec::new(),
        }
    }
    /// Build the ghost particle list from `real_positions`.
    ///
    /// For each particle within `ghost_width` of any face, a translated image
    /// is inserted on the opposite side.  Only the 6 primary face directions
    /// are considered (no edge/corner images).
    pub fn build_ghosts(&mut self, real_positions: &[[f64; 3]]) {
        self.ghost_positions.clear();
        self.ghost_to_real.clear();
        let [lx, ly, lz] = self.box_lengths;
        let w = self.ghost_width;
        for (idx, &p) in real_positions.iter().enumerate() {
            let faces: [([f64; 3], bool); 6] = [
                ([p[0] + lx, p[1], p[2]], p[0] < w),
                ([p[0] - lx, p[1], p[2]], p[0] > lx - w),
                ([p[0], p[1] + ly, p[2]], p[1] < w),
                ([p[0], p[1] - ly, p[2]], p[1] > ly - w),
                ([p[0], p[1], p[2] + lz], p[2] < w),
                ([p[0], p[1], p[2] - lz], p[2] > lz - w),
            ];
            for (ghost_pos, needed) in faces {
                if needed {
                    self.ghost_positions.push(ghost_pos);
                    self.ghost_to_real.push(idx);
                }
            }
        }
    }
    /// Return the number of ghost particles.
    pub fn num_ghosts(&self) -> usize {
        self.ghost_positions.len()
    }
    /// Apply the minimum-image convention to a displacement vector.
    ///
    /// Wraps each component so that it lies in `[-L/2, L/2)`.
    pub fn minimum_image(&self, dx: [f64; 3]) -> [f64; 3] {
        let mut d = dx;
        for k in 0..3 {
            let l = self.box_lengths[k];
            d[k] -= l * (d[k] / l).round();
        }
        d
    }
    /// Wrap a position into the primary simulation box `[0, L)`.
    pub fn wrap_position(&self, p: [f64; 3]) -> [f64; 3] {
        let mut q = p;
        for k in 0..3 {
            let l = self.box_lengths[k];
            q[k] = q[k].rem_euclid(l);
        }
        q
    }
    /// Apply PBC to wrap all positions in `positions` in-place.
    pub fn wrap_all(&self, positions: &mut [[f64; 3]]) {
        for p in positions.iter_mut() {
            *p = self.wrap_position(*p);
        }
    }
}
/// A cell list that maps atoms to their enclosing grid cells.
///
/// After calling [`GpuCellList::build_parallel`] the following invariant holds:
/// for cell `c`, the atoms it owns are
/// `sorted_indices[cell_starts[c\] as usize .. (cell_starts[c] + cell_counts[c]) as usize]`.
pub struct GpuCellList {
    /// Number of cells along x.
    pub nx: usize,
    /// Number of cells along y.
    pub ny: usize,
    /// Number of cells along z.
    pub nz: usize,
    /// Edge length of every cubic cell (same in all directions).
    pub cell_size: f64,
    /// Simulation-box side lengths `[lx, ly, lz]`.
    pub box_lengths: [f64; 3],
    /// `cell_starts[c]` is the first index in `sorted_indices` that belongs to
    /// cell `c`.  Length = `nx * ny * nz`.
    pub cell_starts: Vec<i32>,
    /// Number of atoms in each cell.  Length = `nx * ny * nz`.
    pub cell_counts: Vec<i32>,
    /// Atom indices reordered so that atoms in the same cell are contiguous.
    pub sorted_indices: Vec<usize>,
}
impl GpuCellList {
    /// Create an empty cell list with the given grid dimensions.
    pub fn new(n_cells: [usize; 3], cell_size: f64, box_lengths: [f64; 3]) -> Self {
        let [nx, ny, nz] = n_cells;
        let total = nx * ny * nz;
        Self {
            nx,
            ny,
            nz,
            cell_size,
            box_lengths,
            cell_starts: vec![0; total],
            cell_counts: vec![0; total],
            sorted_indices: Vec::new(),
        }
    }
    /// Build a cell list from a slice of atom positions using Rayon.
    pub fn build_parallel(positions: &[[f64; 3]]) -> Self {
        let n = positions.len();
        let (mut max_x, mut max_y, mut max_z) = (0.0_f64, 0.0_f64, 0.0_f64);
        for p in positions {
            if p[0] > max_x {
                max_x = p[0];
            }
            if p[1] > max_y {
                max_y = p[1];
            }
            if p[2] > max_z {
                max_z = p[2];
            }
        }
        let cell_size = 1.0_f64;
        let nx = ((max_x / cell_size).ceil() as usize).max(1);
        let ny = ((max_y / cell_size).ceil() as usize).max(1);
        let nz = ((max_z / cell_size).ceil() as usize).max(1);
        let box_lengths = [
            nx as f64 * cell_size,
            ny as f64 * cell_size,
            nz as f64 * cell_size,
        ];
        let mut list = Self::new([nx, ny, nz], cell_size, box_lengths);
        if n == 0 {
            return list;
        }
        let atom_cells: Vec<usize> = positions
            .par_iter()
            .map(|&pos| list.cell_index(pos))
            .collect();
        let total_cells = list.total_cells();
        let counts_usize: Vec<usize> = {
            let n_threads = rayon::current_num_threads().max(1);
            let chunk_size = n.div_ceil(n_threads);
            let partial: Vec<Vec<usize>> = atom_cells
                .par_chunks(chunk_size.max(1))
                .map(|chunk| {
                    let mut h = vec![0usize; total_cells];
                    for &c in chunk {
                        h[c] += 1;
                    }
                    h
                })
                .collect();
            let mut merged = vec![0usize; total_cells];
            for h in partial {
                for (i, v) in h.into_iter().enumerate() {
                    merged[i] += v;
                }
            }
            merged
        };
        let starts = parallel_prefix_sum(&counts_usize);
        list.cell_starts = starts.iter().map(|&s| s as i32).collect();
        list.cell_counts = counts_usize.iter().map(|&c| c as i32).collect();
        list.sorted_indices = vec![0usize; n];
        let mut write_pos: Vec<usize> = starts.clone();
        for (atom_idx, &cell) in atom_cells.iter().enumerate() {
            list.sorted_indices[write_pos[cell]] = atom_idx;
            write_pos[cell] += 1;
        }
        list
    }
    /// Decompose `pos` into `(ix, iy, iz)` grid coordinates, clamped to `[0, n-1]`.
    pub fn cell_xyz(&self, pos: [f64; 3]) -> (usize, usize, usize) {
        let ix = ((pos[0] / self.cell_size) as isize).clamp(0, self.nx as isize - 1) as usize;
        let iy = ((pos[1] / self.cell_size) as isize).clamp(0, self.ny as isize - 1) as usize;
        let iz = ((pos[2] / self.cell_size) as isize).clamp(0, self.nz as isize - 1) as usize;
        (ix, iy, iz)
    }
    /// Linearise a 3-D grid coordinate `(x, y, z)` into a flat cell index.
    pub fn cell_index_3d(&self, x: usize, y: usize, z: usize) -> usize {
        z * self.ny * self.nx + y * self.nx + x
    }
    /// Build a cell list from positions and explicit grid parameters.
    pub fn build(
        positions: &[[f64; 3]],
        n_cells: [usize; 3],
        cell_size: f64,
        box_lengths: [f64; 3],
    ) -> Self {
        let [nx, ny, nz] = n_cells;
        let n = positions.len();
        let mut list = Self::new(n_cells, cell_size, box_lengths);
        if n == 0 {
            return list;
        }
        let atom_cells: Vec<usize> = positions
            .par_iter()
            .map(|&pos| list.cell_index(pos))
            .collect();
        let total_cells = nx * ny * nz;
        let counts_usize: Vec<usize> = {
            let n_threads = rayon::current_num_threads().max(1);
            let chunk_size = n.div_ceil(n_threads);
            let partial: Vec<Vec<usize>> = atom_cells
                .par_chunks(chunk_size.max(1))
                .map(|chunk| {
                    let mut h = vec![0usize; total_cells];
                    for &c in chunk {
                        h[c] += 1;
                    }
                    h
                })
                .collect();
            let mut merged = vec![0usize; total_cells];
            for h in partial {
                for (i, v) in h.into_iter().enumerate() {
                    merged[i] += v;
                }
            }
            merged
        };
        let starts = parallel_prefix_sum(&counts_usize);
        list.cell_starts = starts.iter().map(|&s| s as i32).collect();
        list.cell_counts = counts_usize.iter().map(|&c| c as i32).collect();
        list.sorted_indices = vec![0usize; n];
        let mut write_pos: Vec<usize> = starts.clone();
        for (atom_idx, &cell) in atom_cells.iter().enumerate() {
            list.sorted_indices[write_pos[cell]] = atom_idx;
            write_pos[cell] += 1;
        }
        list
    }
    /// Return the linearised 3-D cell index for `pos`, clamped to the grid.
    pub fn cell_index(&self, pos: [f64; 3]) -> usize {
        let ix = ((pos[0] / self.cell_size) as isize).clamp(0, self.nx as isize - 1) as usize;
        let iy = ((pos[1] / self.cell_size) as isize).clamp(0, self.ny as isize - 1) as usize;
        let iz = ((pos[2] / self.cell_size) as isize).clamp(0, self.nz as isize - 1) as usize;
        iz * self.ny * self.nx + iy * self.nx + ix
    }
    /// Return the indices of all atoms within `radius` of `query_pos`.
    pub fn neighbors_in_radius<'a>(
        &'a self,
        positions: &'a [[f64; 3]],
        query_pos: [f64; 3],
        radius: f64,
    ) -> Vec<usize> {
        let r2 = radius * radius;
        let ix_min = (((query_pos[0] - radius) / self.cell_size) as isize)
            .clamp(0, self.nx as isize - 1) as usize;
        let ix_max = (((query_pos[0] + radius) / self.cell_size) as isize)
            .clamp(0, self.nx as isize - 1) as usize;
        let iy_min = (((query_pos[1] - radius) / self.cell_size) as isize)
            .clamp(0, self.ny as isize - 1) as usize;
        let iy_max = (((query_pos[1] + radius) / self.cell_size) as isize)
            .clamp(0, self.ny as isize - 1) as usize;
        let iz_min = (((query_pos[2] - radius) / self.cell_size) as isize)
            .clamp(0, self.nz as isize - 1) as usize;
        let iz_max = (((query_pos[2] + radius) / self.cell_size) as isize)
            .clamp(0, self.nz as isize - 1) as usize;
        let mut result = Vec::new();
        for iz in iz_min..=iz_max {
            for iy in iy_min..=iy_max {
                for ix in ix_min..=ix_max {
                    let cell = iz * self.ny * self.nx + iy * self.nx + ix;
                    let start = self.cell_starts[cell] as usize;
                    let count = self.cell_counts[cell] as usize;
                    for k in start..start + count {
                        let atom = self.sorted_indices[k];
                        let p = positions[atom];
                        let dx = p[0] - query_pos[0];
                        let dy = p[1] - query_pos[1];
                        let dz = p[2] - query_pos[2];
                        if dx * dx + dy * dy + dz * dz <= r2 {
                            result.push(atom);
                        }
                    }
                }
            }
        }
        result
    }
    /// Total number of cells in the grid (`nx * ny * nz`).
    pub fn total_cells(&self) -> usize {
        self.nx * self.ny * self.nz
    }
    /// Return the maximum number of atoms in any single cell.
    pub fn max_cell_occupancy(&self) -> i32 {
        self.cell_counts.iter().cloned().max().unwrap_or(0)
    }
    /// Return the number of non-empty cells.
    pub fn num_nonempty_cells(&self) -> usize {
        self.cell_counts.iter().filter(|&&c| c > 0).count()
    }
    /// Iterate over all pairs of atoms in the same cell or neighbouring cells
    /// within `cutoff` distance, calling `f(i, j, dist_sq)` for each pair.
    pub fn for_each_pair<F>(&self, positions: &[[f64; 3]], cutoff: f64, mut f: F)
    where
        F: FnMut(usize, usize, f64),
    {
        let r2 = cutoff * cutoff;
        for iz in 0..self.nz {
            for iy in 0..self.ny {
                for ix in 0..self.nx {
                    let cell_a = iz * self.ny * self.nx + iy * self.nx + ix;
                    let neighbor_offsets: [(isize, isize, isize); 14] = [
                        (0, 0, 0),
                        (1, 0, 0),
                        (0, 1, 0),
                        (0, 0, 1),
                        (1, 1, 0),
                        (1, 0, 1),
                        (0, 1, 1),
                        (1, 1, 1),
                        (-1, 1, 0),
                        (-1, 0, 1),
                        (0, -1, 1),
                        (-1, 1, 1),
                        (1, -1, 1),
                        (-1, -1, 1),
                    ];
                    for &(dx, dy, dz) in &neighbor_offsets {
                        let nx = ix as isize + dx;
                        let ny_val = iy as isize + dy;
                        let nz_val = iz as isize + dz;
                        if nx < 0
                            || nx >= self.nx as isize
                            || ny_val < 0
                            || ny_val >= self.ny as isize
                            || nz_val < 0
                            || nz_val >= self.nz as isize
                        {
                            continue;
                        }
                        let cell_b = nz_val as usize * self.ny * self.nx
                            + ny_val as usize * self.nx
                            + nx as usize;
                        let start_a = self.cell_starts[cell_a] as usize;
                        let count_a = self.cell_counts[cell_a] as usize;
                        let start_b = self.cell_starts[cell_b] as usize;
                        let count_b = self.cell_counts[cell_b] as usize;
                        for ka in start_a..start_a + count_a {
                            let i = self.sorted_indices[ka];
                            let kb_start = if cell_a == cell_b { ka + 1 } else { start_b };
                            for kb in kb_start..start_b + count_b {
                                let j = self.sorted_indices[kb];
                                let ddx = positions[j][0] - positions[i][0];
                                let ddy = positions[j][1] - positions[i][1];
                                let ddz = positions[j][2] - positions[i][2];
                                let d2 = ddx * ddx + ddy * ddy + ddz * ddz;
                                if d2 <= r2 {
                                    f(i, j, d2);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
/// Summary statistics about how particles are distributed across grid cells.
#[derive(Debug, Clone)]
pub struct OccupancyStats {
    /// Total number of cells in the grid.
    pub total_cells: usize,
    /// Total number of particles.
    pub total_particles: usize,
    /// Number of non-empty cells.
    pub nonempty_cells: usize,
    /// Maximum particles in any single cell.
    pub max_occupancy: usize,
    /// Minimum particles in any non-empty cell.
    pub min_nonempty_occupancy: usize,
    /// Mean particles per non-empty cell.
    pub mean_occupancy: f64,
    /// Variance of particles per non-empty cell.
    pub variance_occupancy: f64,
    /// Load imbalance ratio: `max / mean` (higher is worse for parallel work).
    pub load_imbalance: f64,
}
impl OccupancyStats {
    /// Compute occupancy statistics from a built `GpuCellList`.
    pub fn compute(cl: &GpuCellList) -> Self {
        let total_cells = cl.total_cells();
        let counts: Vec<usize> = cl.cell_counts.iter().map(|&c| c as usize).collect();
        let total_particles: usize = counts.iter().sum();
        let nonempty_cells = counts.iter().filter(|&&c| c > 0).count();
        let max_occupancy = counts.iter().cloned().max().unwrap_or(0);
        let min_nonempty_occupancy = counts
            .iter()
            .filter(|&&c| c > 0)
            .cloned()
            .min()
            .unwrap_or(0);
        let mean_occupancy = if nonempty_cells > 0 {
            total_particles as f64 / nonempty_cells as f64
        } else {
            0.0
        };
        let variance_occupancy = if nonempty_cells > 1 {
            let var_sum: f64 = counts
                .iter()
                .filter(|&&c| c > 0)
                .map(|&c| {
                    let diff = c as f64 - mean_occupancy;
                    diff * diff
                })
                .sum();
            var_sum / nonempty_cells as f64
        } else {
            0.0
        };
        let load_imbalance = if mean_occupancy > 0.0 {
            max_occupancy as f64 / mean_occupancy
        } else {
            1.0
        };
        OccupancyStats {
            total_cells,
            total_particles,
            nonempty_cells,
            max_occupancy,
            min_nonempty_occupancy,
            mean_occupancy,
            variance_occupancy,
            load_imbalance,
        }
    }
    /// Returns `true` if all particles are in a single cell (perfectly unbalanced).
    pub fn is_completely_unbalanced(&self) -> bool {
        self.nonempty_cells == 1
    }
    /// Returns `true` if every cell has at most 1 particle.
    pub fn is_perfectly_spread(&self) -> bool {
        self.max_occupancy <= 1
    }
}
/// A simpler cell-list API that hides grid-dimension bookkeeping.
pub struct CellList {
    pub(super) inner: GpuCellList,
    pub(super) positions: Vec<[f64; 3]>,
}
impl CellList {
    /// Create an empty cell list for a simulation box.
    pub fn new(box_size: [f64; 3], cell_size: f64) -> Self {
        let nx = ((box_size[0] / cell_size).ceil() as usize).max(1);
        let ny = ((box_size[1] / cell_size).ceil() as usize).max(1);
        let nz = ((box_size[2] / cell_size).ceil() as usize).max(1);
        let inner = GpuCellList::new([nx, ny, nz], cell_size, box_size);
        Self {
            inner,
            positions: Vec::new(),
        }
    }
    /// Build a cell list from a slice of positions, choosing `cell_size`
    /// automatically from the bounding box.
    pub fn build(positions: &[[f64; 3]]) -> Self {
        let inner = GpuCellList::build_parallel(positions);
        Self {
            inner,
            positions: positions.to_vec(),
        }
    }
    /// Return indices of all particles within `radius` of `pos`.
    pub fn find_neighbors(&self, pos: [f64; 3], radius: f64) -> Vec<usize> {
        self.inner.neighbors_in_radius(&self.positions, pos, radius)
    }
    /// Build a Verlet neighbour list from this cell list.
    ///
    /// For every atom `i`, the Verlet list stores all atoms `j > i` within
    /// `cutoff + skin` of `i`.  The skin distance (`skin`) must be ≥ 0.
    /// Returns a vector of `(i, j)` index pairs sorted by `i`.
    pub fn build_neighbor_list_verlet(&self, cutoff: f64, skin: f64) -> Vec<(usize, usize)> {
        let r_cut = cutoff + skin.max(0.0);
        let r_cut2 = r_cut * r_cut;
        let n = self.positions.len();
        let mut pairs: Vec<(usize, usize)> = Vec::new();
        for i in 0..n {
            let pi = self.positions[i];
            let (ix, iy, iz) = self.inner.cell_xyz(pi);
            let cell_span = ((r_cut / self.inner.cell_size).ceil() as isize).max(1);
            let ix_lo = (ix as isize - cell_span).max(0) as usize;
            let ix_hi = ((ix as isize + cell_span) as usize).min(self.inner.nx - 1);
            let iy_lo = (iy as isize - cell_span).max(0) as usize;
            let iy_hi = ((iy as isize + cell_span) as usize).min(self.inner.ny - 1);
            let iz_lo = (iz as isize - cell_span).max(0) as usize;
            let iz_hi = ((iz as isize + cell_span) as usize).min(self.inner.nz - 1);
            for cx in ix_lo..=ix_hi {
                for cy in iy_lo..=iy_hi {
                    for cz in iz_lo..=iz_hi {
                        let cell = self.inner.cell_index_3d(cx, cy, cz);
                        let start = self.inner.cell_starts[cell] as usize;
                        let count = self.inner.cell_counts[cell] as usize;
                        for k in start..start + count {
                            let j = self.inner.sorted_indices[k];
                            if j <= i {
                                continue;
                            }
                            let pj = self.positions[j];
                            let dx = pi[0] - pj[0];
                            let dy = pi[1] - pj[1];
                            let dz = pi[2] - pj[2];
                            if dx * dx + dy * dy + dz * dz <= r_cut2 {
                                pairs.push((i, j));
                            }
                        }
                    }
                }
            }
        }
        pairs.sort_unstable();
        pairs
    }
    /// Incrementally update the cell list when particles have moved slightly.
    ///
    /// Re-inserts particles whose displacement from `old_positions` exceeds
    /// `move_threshold` without rebuilding the entire data structure.
    /// Returns the number of particles that were relocated.
    pub fn update_incremental(
        &mut self,
        new_positions: &[[f64; 3]],
        old_positions: &[[f64; 3]],
        move_threshold: f64,
    ) -> usize {
        assert_eq!(new_positions.len(), old_positions.len());
        let thresh2 = move_threshold * move_threshold;
        let mut relocated = 0usize;
        for (i, (np, op)) in new_positions.iter().zip(old_positions.iter()).enumerate() {
            let dx = np[0] - op[0];
            let dy = np[1] - op[1];
            let dz = np[2] - op[2];
            if dx * dx + dy * dy + dz * dz > thresh2 {
                let old_cell = self.inner.cell_index(*op);
                let start = self.inner.cell_starts[old_cell] as usize;
                let count = self.inner.cell_counts[old_cell] as usize;
                if let Some(pos) = self.inner.sorted_indices[start..start + count]
                    .iter()
                    .position(|&idx| idx == i)
                {
                    let abs_pos = start + pos;
                    let last = self.inner.sorted_indices[start + count - 1];
                    self.inner.sorted_indices[abs_pos] = last;
                    self.inner.cell_counts[old_cell] -= 1;
                }
                let new_cell = self.inner.cell_index(*np);
                let new_start = self.inner.cell_starts[new_cell] as usize;
                let new_count = self.inner.cell_counts[new_cell] as usize;
                if new_start + new_count < self.inner.sorted_indices.len() {
                    self.inner.sorted_indices[new_start + new_count] = i;
                    self.inner.cell_counts[new_cell] += 1;
                }
                relocated += 1;
            }
        }
        if new_positions.len() == self.positions.len() {
            self.positions.copy_from_slice(new_positions);
        }
        relocated
    }
    /// Compute a pair-density histogram: count atom pairs in distance bins.
    ///
    /// Bins are `[0, dr)`, `[dr, 2*dr)`, … up to `max_r`.
    /// Returns a `Vec`usize` of length `ceil(max_r / dr)`.
    pub fn compute_pair_density(&self, max_r: f64, dr: f64) -> Vec<usize> {
        let dr = dr.max(1e-15);
        let n_bins = (max_r / dr).ceil() as usize;
        let mut hist = vec![0usize; n_bins];
        let n = self.positions.len();
        for i in 0..n {
            let pi = self.positions[i];
            for j in (i + 1)..n {
                let pj = self.positions[j];
                let dx = pi[0] - pj[0];
                let dy = pi[1] - pj[1];
                let dz = pi[2] - pj[2];
                let r = (dx * dx + dy * dy + dz * dz).sqrt();
                if r < max_r {
                    let bin = (r / dr) as usize;
                    if bin < n_bins {
                        hist[bin] += 1;
                    }
                }
            }
        }
        hist
    }
}
/// A spatial hash map for fast neighbor queries with arbitrary cell sizes.
///
/// Unlike the grid-based cell list, this uses a hash table and doesn't
/// require a bounded domain.
pub struct SpatialHash {
    /// Hash table: maps hash -> list of particle indices.
    pub(super) table: Vec<Vec<usize>>,
    /// Number of buckets.
    pub(super) num_buckets: usize,
    /// Cell size.
    pub(super) cell_size: f64,
    /// Inverse cell size (cached).
    pub(super) inv_cell_size: f64,
}
impl SpatialHash {
    /// Create a new spatial hash with the given number of buckets and cell size.
    pub fn new(num_buckets: usize, cell_size: f64) -> Self {
        Self {
            table: vec![Vec::new(); num_buckets],
            num_buckets,
            cell_size,
            inv_cell_size: 1.0 / cell_size,
        }
    }
    /// Clear all entries.
    pub fn clear(&mut self) {
        for bucket in &mut self.table {
            bucket.clear();
        }
    }
    /// Hash a 3D integer cell coordinate to a bucket index.
    fn hash(&self, ix: i64, iy: i64, iz: i64) -> usize {
        let h = (ix.wrapping_mul(73856093) ^ iy.wrapping_mul(19349669) ^ iz.wrapping_mul(83492791))
            .unsigned_abs() as usize;
        h % self.num_buckets
    }
    /// Insert a particle into the hash.
    pub fn insert(&mut self, index: usize, pos: [f64; 3]) {
        let ix = (pos[0] * self.inv_cell_size).floor() as i64;
        let iy = (pos[1] * self.inv_cell_size).floor() as i64;
        let iz = (pos[2] * self.inv_cell_size).floor() as i64;
        let bucket = self.hash(ix, iy, iz);
        self.table[bucket].push(index);
    }
    /// Build the spatial hash from a set of positions.
    pub fn build(&mut self, positions: &[[f64; 3]]) {
        self.clear();
        for (i, &p) in positions.iter().enumerate() {
            self.insert(i, p);
        }
    }
    /// Find all particles within `radius` of `query_pos`.
    pub fn query_radius(
        &self,
        positions: &[[f64; 3]],
        query_pos: [f64; 3],
        radius: f64,
    ) -> Vec<usize> {
        let r2 = radius * radius;
        let cells_to_check = (radius * self.inv_cell_size).ceil() as i64 + 1;
        let qx = (query_pos[0] * self.inv_cell_size).floor() as i64;
        let qy = (query_pos[1] * self.inv_cell_size).floor() as i64;
        let qz = (query_pos[2] * self.inv_cell_size).floor() as i64;
        let mut result = Vec::new();
        for ddz in -cells_to_check..=cells_to_check {
            for ddy in -cells_to_check..=cells_to_check {
                for ddx in -cells_to_check..=cells_to_check {
                    let bucket = self.hash(qx + ddx, qy + ddy, qz + ddz);
                    for &idx in &self.table[bucket] {
                        let p = positions[idx];
                        let dx = p[0] - query_pos[0];
                        let dy = p[1] - query_pos[1];
                        let dz = p[2] - query_pos[2];
                        let d2 = dx * dx + dy * dy + dz * dz;
                        if d2 <= r2 {
                            result.push(idx);
                        }
                    }
                }
            }
        }
        result
    }
    /// Total number of stored particles.
    pub fn len(&self) -> usize {
        self.table.iter().map(|b| b.len()).sum::<usize>()
    }
    /// Whether the hash is empty.
    pub fn is_empty(&self) -> bool {
        self.table.iter().all(|b| b.is_empty())
    }
}
