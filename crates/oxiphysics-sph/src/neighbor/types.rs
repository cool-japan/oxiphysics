//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use oxiphysics_core::math::Vec3;
use std::collections::HashMap;

use super::functions::{CellKey, CellKey3};

/// Manages the skin distance for Verlet-list amortised rebuilds.
///
/// Tracks the total displacement of each particle since the last rebuild,
/// accounting for multiple sub-steps between rebuilds.  A rebuild is
/// triggered when the *maximum* cumulative displacement of any particle
/// exceeds `r_skin / 2`.
#[derive(Debug, Clone)]
pub struct SkinManager {
    /// Interaction cut-off radius.
    pub r_cut: f64,
    /// Skin radius (added buffer beyond r_cut for neighbour list).
    pub r_skin: f64,
    /// Reference positions at the last rebuild.
    pub(super) ref_positions: Vec<[f64; 3]>,
    /// Total number of rebuilds performed.
    pub rebuild_count: usize,
    /// Total number of update calls.
    pub update_count: usize,
}
impl SkinManager {
    /// Initialise from a set of particle positions.
    pub fn new(positions: &[[f64; 3]], r_cut: f64, r_skin: f64) -> Self {
        Self {
            r_cut,
            r_skin,
            ref_positions: positions.to_vec(),
            rebuild_count: 0,
            update_count: 0,
        }
    }
    /// Check whether any particle has moved more than `r_skin / 2` since
    /// the last rebuild.  Call this every sub-step.
    ///
    /// Returns `true` if a rebuild is required.
    pub fn needs_rebuild(&self, positions: &[[f64; 3]]) -> bool {
        let half_skin2 = (self.r_skin * 0.5) * (self.r_skin * 0.5);
        positions
            .iter()
            .zip(self.ref_positions.iter())
            .any(|(cur, ref_p)| {
                let d2 = (cur[0] - ref_p[0]).powi(2)
                    + (cur[1] - ref_p[1]).powi(2)
                    + (cur[2] - ref_p[2]).powi(2);
                d2 > half_skin2
            })
    }
    /// Commit a rebuild: update the reference positions and counters.
    pub fn commit_rebuild(&mut self, positions: &[[f64; 3]]) {
        self.ref_positions = positions.to_vec();
        self.rebuild_count += 1;
    }
    /// Called each sub-step.  Checks whether a rebuild is needed and, if so,
    /// commits one.  Returns `true` if a rebuild occurred.
    pub fn update(&mut self, positions: &[[f64; 3]]) -> bool {
        self.update_count += 1;
        if self.needs_rebuild(positions) {
            self.commit_rebuild(positions);
            return true;
        }
        false
    }
    /// Rebuild interval: average number of update calls between rebuilds.
    pub fn avg_rebuild_interval(&self) -> f64 {
        if self.rebuild_count == 0 {
            return self.update_count as f64;
        }
        self.update_count as f64 / self.rebuild_count as f64
    }
    /// Current maximum displacement of any particle since last rebuild.
    pub fn max_displacement(&self, positions: &[[f64; 3]]) -> f64 {
        positions
            .iter()
            .zip(self.ref_positions.iter())
            .map(|(cur, ref_p)| {
                ((cur[0] - ref_p[0]).powi(2)
                    + (cur[1] - ref_p[1]).powi(2)
                    + (cur[2] - ref_p[2]).powi(2))
                .sqrt()
            })
            .fold(0.0_f64, f64::max)
    }
}
/// Linked-cell list with periodic boundary conditions.
///
/// Wraps positions into `[0, L)³` before binning and applies the minimum-image
/// convention when computing distances.  Builds in O(N) time and queries in
/// O(N/M) time where M is the number of cells.
#[derive(Debug, Clone)]
pub struct PbcLinkedCellList {
    /// Box dimensions \[Lx, Ly, Lz\].
    pub box_size: [f64; 3],
    /// Cell size (≥ r_cut).
    pub cell_size: f64,
    /// Number of cells per dimension.
    pub n_cells: [usize; 3],
    /// Head-of-chain array.
    pub(super) head: Vec<usize>,
    /// Next-in-chain array.
    pub(super) next: Vec<usize>,
    /// Stored positions (wrapped into box).
    pub(super) positions: Vec<[f64; 3]>,
}
impl PbcLinkedCellList {
    /// Build from positions in a periodic box.
    pub fn build(positions: &[[f64; 3]], cell_size: f64, box_size: [f64; 3]) -> Self {
        let n = positions.len();
        let n_cells = [
            (box_size[0] / cell_size).floor().max(1.0) as usize,
            (box_size[1] / cell_size).floor().max(1.0) as usize,
            (box_size[2] / cell_size).floor().max(1.0) as usize,
        ];
        let total_cells = n_cells[0] * n_cells[1] * n_cells[2];
        let mut head = vec![usize::MAX; total_cells];
        let mut next = vec![usize::MAX; n];
        let wrapped: Vec<[f64; 3]> = positions
            .iter()
            .map(|&p| {
                [
                    p[0].rem_euclid(box_size[0]),
                    p[1].rem_euclid(box_size[1]),
                    p[2].rem_euclid(box_size[2]),
                ]
            })
            .collect();
        for (i, w) in wrapped.iter().enumerate() {
            let cx = ((w[0] / cell_size).floor() as usize).min(n_cells[0] - 1);
            let cy = ((w[1] / cell_size).floor() as usize).min(n_cells[1] - 1);
            let cz = ((w[2] / cell_size).floor() as usize).min(n_cells[2] - 1);
            let idx = cx * n_cells[1] * n_cells[2] + cy * n_cells[2] + cz;
            next[i] = head[idx];
            head[idx] = i;
        }
        Self {
            box_size,
            cell_size,
            n_cells,
            head,
            next,
            positions: wrapped,
        }
    }
    /// Minimum-image displacement vector ri − rj.
    #[inline]
    fn min_image(&self, ri: [f64; 3], rj: [f64; 3]) -> [f64; 3] {
        let mut dr = [0.0_f64; 3];
        for k in 0..3 {
            let mut d = ri[k] - rj[k];
            let l = self.box_size[k];
            d -= (d / l).round() * l;
            dr[k] = d;
        }
        dr
    }
    /// Query all particle indices within `r_cut` of particle `i` under PBC.
    pub fn query(&self, i: usize, r_cut: f64) -> Vec<usize> {
        let pos_i = self.positions[i];
        let r2 = r_cut * r_cut;
        let mut result = Vec::new();
        for ox in -1_i64..=1 {
            for oy in -1_i64..=1 {
                for oz in -1_i64..=1 {
                    let ci_x = ((pos_i[0] / self.cell_size).floor() as i64 + ox)
                        .rem_euclid(self.n_cells[0] as i64) as usize;
                    let ci_y = ((pos_i[1] / self.cell_size).floor() as i64 + oy)
                        .rem_euclid(self.n_cells[1] as i64) as usize;
                    let ci_z = ((pos_i[2] / self.cell_size).floor() as i64 + oz)
                        .rem_euclid(self.n_cells[2] as i64) as usize;
                    let idx =
                        ci_x * self.n_cells[1] * self.n_cells[2] + ci_y * self.n_cells[2] + ci_z;
                    let mut j = self.head[idx];
                    while j != usize::MAX {
                        if j != i {
                            let dr = self.min_image(pos_i, self.positions[j]);
                            let d2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
                            if d2 <= r2 {
                                result.push(j);
                            }
                        }
                        j = self.next[j];
                    }
                }
            }
        }
        result.sort_unstable();
        result
    }
    /// Build all-particle neighbor lists under PBC.
    pub fn all_neighbors(
        positions: &[[f64; 3]],
        r_cut: f64,
        box_size: [f64; 3],
    ) -> Vec<Vec<usize>> {
        let lcl = Self::build(positions, r_cut, box_size);
        (0..positions.len()).map(|i| lcl.query(i, r_cut)).collect()
    }
    /// Number of cells.
    pub fn num_cells(&self) -> usize {
        self.head.len()
    }
}
/// Linked-cell list for O(N) spatial neighbour queries.
///
/// Particles are binned into cells of size `cell_size`.  To find neighbours
/// of particle `i`, iterate over particles in the 27 surrounding cells.
#[derive(Debug, Clone)]
pub struct LinkedCellList {
    /// Cell size (usually = cut-off radius).
    pub(super) cell_size: f64,
    /// Domain lower-left corner.
    pub(super) origin: [f64; 3],
    /// Number of cells in each dimension.
    pub(super) n_cells: [usize; 3],
    /// Head-of-chain array: `head[cell_index]` = first particle in cell,
    /// or `usize::MAX` if empty.
    pub(super) head: Vec<usize>,
    /// Next-in-chain array: `next[particle]` = next particle in same cell,
    /// or `usize::MAX` if last.
    pub(super) next: Vec<usize>,
}
impl LinkedCellList {
    /// Build the linked-cell list from positions.
    ///
    /// The domain is automatically determined from the bounding box of the
    /// positions, padded by one cell on each side.
    pub fn build(positions: &[[f64; 3]], cell_size: f64) -> Self {
        let n = positions.len();
        if n == 0 {
            return Self {
                cell_size,
                origin: [0.0; 3],
                n_cells: [1, 1, 1],
                head: vec![usize::MAX; 1],
                next: Vec::new(),
            };
        }
        let mut lo = positions[0];
        let mut hi = positions[0];
        for pos in &positions[1..] {
            for k in 0..3 {
                if pos[k] < lo[k] {
                    lo[k] = pos[k];
                }
                if pos[k] > hi[k] {
                    hi[k] = pos[k];
                }
            }
        }
        for k in 0..3 {
            lo[k] -= cell_size;
            hi[k] += cell_size;
        }
        let n_cells = [
            ((hi[0] - lo[0]) / cell_size).ceil().max(1.0) as usize,
            ((hi[1] - lo[1]) / cell_size).ceil().max(1.0) as usize,
            ((hi[2] - lo[2]) / cell_size).ceil().max(1.0) as usize,
        ];
        let total_cells = n_cells[0] * n_cells[1] * n_cells[2];
        let mut head = vec![usize::MAX; total_cells];
        let mut next = vec![usize::MAX; n];
        for (i, pos) in positions.iter().enumerate() {
            let cx = ((pos[0] - lo[0]) / cell_size).floor() as usize;
            let cy = ((pos[1] - lo[1]) / cell_size).floor() as usize;
            let cz = ((pos[2] - lo[2]) / cell_size).floor() as usize;
            let cx = cx.min(n_cells[0] - 1);
            let cy = cy.min(n_cells[1] - 1);
            let cz = cz.min(n_cells[2] - 1);
            let idx = cx * n_cells[1] * n_cells[2] + cy * n_cells[2] + cz;
            next[i] = head[idx];
            head[idx] = i;
        }
        Self {
            cell_size,
            origin: lo,
            n_cells,
            head,
            next,
        }
    }
    /// Find all neighbours of position `pos` within `radius`.
    ///
    /// Iterates over the 27 surrounding cells for O(1) expected cost per query.
    pub fn query(&self, positions: &[[f64; 3]], pos: [f64; 3], radius: f64) -> Vec<usize> {
        let r2 = radius * radius;
        let cx_lo = ((pos[0] - radius - self.origin[0]) / self.cell_size).floor() as i32;
        let cy_lo = ((pos[1] - radius - self.origin[1]) / self.cell_size).floor() as i32;
        let cz_lo = ((pos[2] - radius - self.origin[2]) / self.cell_size).floor() as i32;
        let cx_hi = ((pos[0] + radius - self.origin[0]) / self.cell_size).floor() as i32;
        let cy_hi = ((pos[1] + radius - self.origin[1]) / self.cell_size).floor() as i32;
        let cz_hi = ((pos[2] + radius - self.origin[2]) / self.cell_size).floor() as i32;
        let mut result = Vec::new();
        for cx in cx_lo..=cx_hi {
            if cx < 0 || cx >= self.n_cells[0] as i32 {
                continue;
            }
            for cy in cy_lo..=cy_hi {
                if cy < 0 || cy >= self.n_cells[1] as i32 {
                    continue;
                }
                for cz in cz_lo..=cz_hi {
                    if cz < 0 || cz >= self.n_cells[2] as i32 {
                        continue;
                    }
                    let idx = cx as usize * self.n_cells[1] * self.n_cells[2]
                        + cy as usize * self.n_cells[2]
                        + cz as usize;
                    let mut p = self.head[idx];
                    while p != usize::MAX {
                        let d2 = (positions[p][0] - pos[0]).powi(2)
                            + (positions[p][1] - pos[1]).powi(2)
                            + (positions[p][2] - pos[2]).powi(2);
                        if d2 <= r2 {
                            result.push(p);
                        }
                        p = self.next[p];
                    }
                }
            }
        }
        result
    }
    /// Find all neighbours for every particle (excluding self).
    pub fn all_neighbors(positions: &[[f64; 3]], radius: f64) -> Vec<Vec<usize>> {
        let lcl = Self::build(positions, radius);
        positions
            .iter()
            .enumerate()
            .map(|(i, &pos)| {
                let mut nb: Vec<usize> = lcl
                    .query(positions, pos, radius)
                    .into_iter()
                    .filter(|&j| j != i)
                    .collect();
                nb.sort_unstable();
                nb
            })
            .collect()
    }
    /// Total number of cells.
    pub fn num_cells(&self) -> usize {
        self.head.len()
    }
}
impl LinkedCellList {
    /// Compute the occupancy histogram of the linked-cell list.
    ///
    /// Returns an `OccupancyHistogram` where `bins[k]` is the number of cells
    /// containing exactly `k` particles.
    pub fn compute_occupancy_histogram(&self) -> OccupancyHistogram {
        let total_cells = self.head.len();
        let mut counts = vec![0_usize; total_cells];
        for (c, cnt) in counts.iter_mut().enumerate() {
            let mut p = self.head[c];
            while p != usize::MAX {
                *cnt += 1;
                p = self.next[p];
            }
        }
        let peak = *counts.iter().max().unwrap_or(&0);
        let mut bins = vec![0_usize; peak + 1];
        for &cnt in &counts {
            bins[cnt] += 1;
        }
        let occupied_cells = counts.iter().filter(|&&c| c > 0).count();
        OccupancyHistogram {
            bins,
            occupied_cells,
            peak_occupancy: peak,
        }
    }
}
/// Statistics computed over a set of neighbor lists.
#[derive(Debug, Clone)]
pub struct NeighborStats {
    /// Total number of particles.
    pub n_particles: usize,
    /// Minimum neighbor count.
    pub min_neighbors: usize,
    /// Maximum neighbor count.
    pub max_neighbors: usize,
    /// Average neighbor count.
    pub avg_neighbors: f64,
    /// Variance of neighbor count.
    pub variance: f64,
    /// Standard deviation of neighbor count.
    pub std_dev: f64,
    /// Total neighbor pairs (sum of all per-particle counts).
    pub total_pairs: usize,
}
impl NeighborStats {
    /// Compute statistics from per-particle neighbor lists.
    pub fn compute(neighbor_lists: &[Vec<usize>]) -> Self {
        let n = neighbor_lists.len();
        if n == 0 {
            return Self {
                n_particles: 0,
                min_neighbors: 0,
                max_neighbors: 0,
                avg_neighbors: 0.0,
                variance: 0.0,
                std_dev: 0.0,
                total_pairs: 0,
            };
        }
        let counts: Vec<f64> = neighbor_lists.iter().map(|v| v.len() as f64).collect();
        let total: f64 = counts.iter().sum();
        let avg = total / n as f64;
        let min_n = neighbor_lists.iter().map(|v| v.len()).min().unwrap_or(0);
        let max_n = neighbor_lists.iter().map(|v| v.len()).max().unwrap_or(0);
        let var = if n > 1 {
            counts.iter().map(|&c| (c - avg) * (c - avg)).sum::<f64>() / (n - 1) as f64
        } else {
            0.0
        };
        Self {
            n_particles: n,
            min_neighbors: min_n,
            max_neighbors: max_n,
            avg_neighbors: avg,
            variance: var,
            std_dev: var.sqrt(),
            total_pairs: total as usize,
        }
    }
    /// Compute statistics using LinkedCellList.
    pub fn from_positions(positions: &[[f64; 3]], h: f64) -> Self {
        let nls = LinkedCellList::all_neighbors(positions, h);
        Self::compute(&nls)
    }
}
/// Cached neighbor list with lazy rebuild triggered by a displacement threshold.
///
/// Wraps a [`LinkedCellList`]-based neighbor list and tracks cumulative
/// displacement since the last rebuild.  A rebuild is triggered when *any*
/// particle exceeds `r_rebuild = r_skin / 2`.
#[derive(Debug, Clone)]
pub struct NeighborListCache {
    /// Cut-off radius for neighbor interactions.
    pub(super) r_cut: f64,
    /// Skin radius for amortising rebuilds.
    pub(super) r_skin: f64,
    /// Reference positions at the last rebuild.
    pub(super) ref_positions: Vec<[f64; 3]>,
    /// Flat neighbor list (GPU-friendly layout).
    pub flat: FlatNeighborList,
    /// Number of rebuilds performed.
    pub rebuild_count: usize,
}
impl NeighborListCache {
    /// Build the cache from the given positions.
    pub fn new(positions: &[[f64; 3]], r_cut: f64, r_skin: f64) -> Self {
        let flat = FlatNeighborList::build(positions, r_cut + r_skin);
        Self {
            r_cut,
            r_skin,
            ref_positions: positions.to_vec(),
            flat,
            rebuild_count: 0,
        }
    }
    /// Return the neighbors of particle `i`.
    pub fn neighbors_of(&self, i: usize) -> &[usize] {
        self.flat.neighbors_of(i)
    }
    /// Check displacement and rebuild if needed; returns `true` if a rebuild occurred.
    pub fn update(&mut self, positions: &[[f64; 3]]) -> bool {
        let threshold2 = (self.r_skin * 0.5) * (self.r_skin * 0.5);
        let needs = positions
            .iter()
            .zip(&self.ref_positions)
            .any(|(cur, ref_pos)| {
                let d2 = (cur[0] - ref_pos[0]).powi(2)
                    + (cur[1] - ref_pos[1]).powi(2)
                    + (cur[2] - ref_pos[2]).powi(2);
                d2 > threshold2
            });
        if needs {
            self.flat = FlatNeighborList::build(positions, self.r_cut + self.r_skin);
            self.ref_positions = positions.to_vec();
            self.rebuild_count += 1;
            true
        } else {
            false
        }
    }
    /// Number of particles tracked.
    pub fn len(&self) -> usize {
        self.flat.n_particles
    }
    /// Whether the cache is empty.
    pub fn is_empty(&self) -> bool {
        self.flat.n_particles == 0
    }
}
/// Cell-based spatial hash for neighbor queries using plain `[f64; 3]` positions.
///
/// This is a zero-dependency alternative to [`SpatialHash`] for code that
/// cannot or does not want to depend on nalgebra/`Vec3`.
///
/// # Example
/// ```no_run
/// use oxiphysics_sph::neighbor::SpatialHash3D;
///
/// let positions: Vec<[f64; 3]> = vec![
///     [0.0, 0.0, 0.0],
///     [0.05, 0.0, 0.0],
///     [10.0, 10.0, 10.0],
/// ];
/// let hash = SpatialHash3D::build(&positions, 0.2);
/// let nbrs = hash.find_neighbors([0.0, 0.0, 0.0]);
/// assert!(nbrs.contains(&0));
/// assert!(nbrs.contains(&1));
/// assert!(!nbrs.contains(&2));
/// ```
#[derive(Debug, Clone)]
pub struct SpatialHash3D {
    /// Smoothing radius used as cell size.
    pub(super) h: f64,
    /// Particle positions (copy stored for distance checks).
    pub(super) positions: Vec<[f64; 3]>,
    /// Map from cell coordinates to particle indices.
    pub(super) cells: HashMap<CellKey3, Vec<usize>>,
}
impl SpatialHash3D {
    /// Build a spatial hash from `positions` using `h` as cell/support size.
    ///
    /// Every particle at `positions[i]` is inserted into the cell that
    /// contains `positions[i]`.  [`Self::find_neighbors`] then queries the
    /// 3×3×3 = 27 surrounding cells and returns all particles within
    /// distance `h`.
    pub fn build(positions: &[[f64; 3]], h: f64) -> Self {
        let mut cells: HashMap<CellKey3, Vec<usize>> = HashMap::new();
        for (i, pos) in positions.iter().enumerate() {
            let key = Self::cell_key_of(pos, h);
            cells.entry(key).or_default().push(i);
        }
        Self {
            h,
            positions: positions.to_vec(),
            cells,
        }
    }
    /// Return all particle indices (including self) within distance `h` of `pos`.
    ///
    /// The returned list is unsorted and may include the query position itself
    /// if it coincides with one of the stored positions.
    pub fn find_neighbors(&self, pos: [f64; 3]) -> Vec<usize> {
        let h2 = self.h * self.h;
        let ck = Self::cell_key_of(&pos, self.h);
        let mut result = Vec::new();
        for dx in -1_i32..=1 {
            for dy in -1_i32..=1 {
                for dz in -1_i32..=1 {
                    let key = (ck.0 + dx, ck.1 + dy, ck.2 + dz);
                    if let Some(indices) = self.cells.get(&key) {
                        for &idx in indices {
                            let p = self.positions[idx];
                            let d2 = (p[0] - pos[0]).powi(2)
                                + (p[1] - pos[1]).powi(2)
                                + (p[2] - pos[2]).powi(2);
                            if d2 <= h2 {
                                result.push(idx);
                            }
                        }
                    }
                }
            }
        }
        result
    }
    /// Build the hash and compute per-particle neighbor lists (excluding self).
    ///
    /// Returns `neighbors[i]` = sorted list of neighbor indices for particle `i`.
    pub fn all_neighbors(positions: &[[f64; 3]], h: f64) -> Vec<Vec<usize>> {
        let grid = Self::build(positions, h);
        positions
            .iter()
            .enumerate()
            .map(|(i, &pos)| {
                let mut nb: Vec<usize> = grid
                    .find_neighbors(pos)
                    .into_iter()
                    .filter(|&j| j != i)
                    .collect();
                nb.sort_unstable();
                nb
            })
            .collect()
    }
    /// Compute the integer cell key for a position.
    #[inline]
    fn cell_key_of(pos: &[f64; 3], h: f64) -> CellKey3 {
        (
            (pos[0] / h).floor() as i32,
            (pos[1] / h).floor() as i32,
            (pos[2] / h).floor() as i32,
        )
    }
}
/// Neighbor search with per-particle variable smoothing lengths.
///
/// Each particle `i` has its own smoothing length `h_i`.  A pair `(i, j)` is
/// considered neighbors when `|r_ij| ≤ max(h_i, h_j)` (the "gather–scatter"
/// convention commonly used in adaptive SPH).
///
/// Internally uses a global spatial hash with cell size = `max(h)` to bound
/// the search radius.
#[derive(Debug, Clone)]
pub struct AdaptiveNeighborSearch {
    /// Particle positions.
    pub(super) positions: Vec<[f64; 3]>,
    /// Per-particle smoothing lengths.
    pub h: Vec<f64>,
    /// Global cell size (= max h_i).
    pub(super) cell_size: f64,
    /// Spatial hash cells.
    pub(super) cells: std::collections::HashMap<(i32, i32, i32), Vec<usize>>,
}
impl AdaptiveNeighborSearch {
    /// Build the adaptive neighbor search structure.
    ///
    /// # Panics
    /// Panics if `positions` and `h` have different lengths or if `h` is empty.
    pub fn build(positions: &[[f64; 3]], h: &[f64]) -> Self {
        assert_eq!(
            positions.len(),
            h.len(),
            "positions and h must have equal length"
        );
        let cell_size = h
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max)
            .max(1e-14);
        let mut cells: std::collections::HashMap<(i32, i32, i32), Vec<usize>> =
            std::collections::HashMap::new();
        for (i, pos) in positions.iter().enumerate() {
            let key = (
                (pos[0] / cell_size).floor() as i32,
                (pos[1] / cell_size).floor() as i32,
                (pos[2] / cell_size).floor() as i32,
            );
            cells.entry(key).or_default().push(i);
        }
        Self {
            positions: positions.to_vec(),
            h: h.to_vec(),
            cell_size,
            cells,
        }
    }
    /// Compute neighbor list for particle `i` using the gather–scatter convention.
    ///
    /// Particle `j` is included when `|r_ij| ≤ max(h_i, h_j)`.
    pub fn query(&self, i: usize) -> Vec<usize> {
        let pos_i = self.positions[i];
        let h_i = self.h[i];
        let r_search = h_i;
        let cx0 = ((pos_i[0] - r_search) / self.cell_size).floor() as i32;
        let cy0 = ((pos_i[1] - r_search) / self.cell_size).floor() as i32;
        let cz0 = ((pos_i[2] - r_search) / self.cell_size).floor() as i32;
        let cx1 = ((pos_i[0] + r_search) / self.cell_size).floor() as i32;
        let cy1 = ((pos_i[1] + r_search) / self.cell_size).floor() as i32;
        let cz1 = ((pos_i[2] + r_search) / self.cell_size).floor() as i32;
        let mut result = Vec::new();
        for cx in cx0..=cx1 {
            for cy in cy0..=cy1 {
                for cz in cz0..=cz1 {
                    if let Some(indices) = self.cells.get(&(cx, cy, cz)) {
                        for &j in indices {
                            if j == i {
                                continue;
                            }
                            let dx = pos_i[0] - self.positions[j][0];
                            let dy = pos_i[1] - self.positions[j][1];
                            let dz = pos_i[2] - self.positions[j][2];
                            let dist = (dx * dx + dy * dy + dz * dz).sqrt();
                            let h_max = h_i.max(self.h[j]);
                            if dist <= h_max {
                                result.push(j);
                            }
                        }
                    }
                }
            }
        }
        result.sort_unstable();
        result
    }
    /// Compute per-particle neighbor lists for all particles.
    pub fn all_neighbors(&self) -> Vec<Vec<usize>> {
        (0..self.positions.len()).map(|i| self.query(i)).collect()
    }
    /// Suggest a new smoothing length for particle `i` targeting `n_target` neighbors.
    ///
    /// Uses a simple bisection-free estimate: h_new = h_old × (n_target / n_actual)^(1/3).
    /// Falls back to the current h if particle has no neighbors.
    pub fn suggest_h(&self, i: usize, n_target: usize) -> f64 {
        let n_actual = self.query(i).len();
        if n_actual == 0 {
            return self.h[i];
        }
        let ratio = n_target as f64 / n_actual as f64;
        self.h[i] * ratio.cbrt()
    }
}
/// Per-particle neighbor count together with global statistics.
#[derive(Debug, Clone)]
pub struct NeighborCountStats {
    /// Per-particle neighbor count.
    pub counts: Vec<usize>,
    /// Minimum count.
    pub min: usize,
    /// Maximum count.
    pub max: usize,
    /// Mean count.
    pub mean: f64,
    /// Standard deviation.
    pub std_dev: f64,
}
impl NeighborCountStats {
    /// Compute from per-particle neighbor lists.
    pub fn from_neighbor_lists(neighbor_lists: &[Vec<usize>]) -> Self {
        let counts: Vec<usize> = neighbor_lists.iter().map(|v| v.len()).collect();
        let n = counts.len();
        if n == 0 {
            return Self {
                counts: Vec::new(),
                min: 0,
                max: 0,
                mean: 0.0,
                std_dev: 0.0,
            };
        }
        let min = counts.iter().copied().min().unwrap_or(0);
        let max = counts.iter().copied().max().unwrap_or(0);
        let mean = counts.iter().sum::<usize>() as f64 / n as f64;
        let var = if n > 1 {
            counts
                .iter()
                .map(|&c| (c as f64 - mean).powi(2))
                .sum::<f64>()
                / (n - 1) as f64
        } else {
            0.0
        };
        Self {
            counts,
            min,
            max,
            mean,
            std_dev: var.sqrt(),
        }
    }
    /// Build from positions.
    pub fn from_positions(positions: &[[f64; 3]], h: f64) -> Self {
        let nls = LinkedCellList::all_neighbors(positions, h);
        Self::from_neighbor_lists(&nls)
    }
    /// Fraction of particles with at least `n` neighbors.
    pub fn fraction_with_at_least(&self, n: usize) -> f64 {
        if self.counts.is_empty() {
            return 0.0;
        }
        self.counts.iter().filter(|&&c| c >= n).count() as f64 / self.counts.len() as f64
    }
}
/// Histogram of cell occupancies: `histogram[k]` = number of cells with `k` particles.
#[derive(Debug, Clone)]
pub struct OccupancyHistogram {
    /// Histogram bins (index = occupancy count).
    pub bins: Vec<usize>,
    /// Total number of non-empty cells.
    pub occupied_cells: usize,
    /// Cell with the highest occupancy.
    pub peak_occupancy: usize,
}
/// Node in a simple octree for neighbor search.
#[derive(Debug, Clone)]
pub(super) enum OctreeNode {
    Leaf {
        /// Indices of particles in this leaf.
        indices: Vec<usize>,
        /// Bounding box of this leaf.
        bounds_lo: [f64; 3],
        bounds_hi: [f64; 3],
    },
    Internal {
        /// Bounding box of this node.
        bounds_lo: [f64; 3],
        bounds_hi: [f64; 3],
        /// Up to 8 children (octants).
        children: Vec<OctreeNode>,
    },
}
/// Verlet neighbour list with a skin radius for amortised rebuilds.
///
/// Neighbours are computed within `r_cut + r_skin` and cached.  A rebuild
/// is triggered when any particle has moved more than `r_skin / 2` since
/// the last build (guaranteeing no missed neighbours).
#[derive(Debug, Clone)]
pub struct VerletList {
    /// Cut-off radius.
    pub(super) r_cut: f64,
    /// Skin radius.
    pub(super) r_skin: f64,
    /// Positions at the time of the last build.
    pub(super) ref_positions: Vec<[f64; 3]>,
    /// Neighbour lists (excluding self).
    pub(super) neighbor_lists: Vec<Vec<usize>>,
}
impl VerletList {
    /// Build a new Verlet list from the given positions.
    pub fn build(positions: &[[f64; 3]], r_cut: f64, r_skin: f64) -> Self {
        let r_total = r_cut + r_skin;
        let neighbor_lists = Self::compute_neighbors(positions, r_total);
        Self {
            r_cut,
            r_skin,
            ref_positions: positions.to_vec(),
            neighbor_lists,
        }
    }
    /// Check whether any particle has moved far enough to require a rebuild.
    ///
    /// Returns `true` if the list should be rebuilt.
    pub fn needs_rebuild(&self, positions: &[[f64; 3]]) -> bool {
        let threshold = self.r_skin * 0.5;
        let threshold2 = threshold * threshold;
        for (cur, ref_pos) in positions.iter().zip(self.ref_positions.iter()) {
            let d2 = (cur[0] - ref_pos[0]).powi(2)
                + (cur[1] - ref_pos[1]).powi(2)
                + (cur[2] - ref_pos[2]).powi(2);
            if d2 > threshold2 {
                return true;
            }
        }
        false
    }
    /// Get the neighbour list for particle `i`.
    pub fn neighbors(&self, i: usize) -> &[usize] {
        &self.neighbor_lists[i]
    }
    /// Number of particles.
    pub fn len(&self) -> usize {
        self.neighbor_lists.len()
    }
    /// Whether empty.
    pub fn is_empty(&self) -> bool {
        self.neighbor_lists.is_empty()
    }
    /// Rebuild the list in-place.
    pub fn rebuild(&mut self, positions: &[[f64; 3]]) {
        let r_total = self.r_cut + self.r_skin;
        self.neighbor_lists = Self::compute_neighbors(positions, r_total);
        self.ref_positions = positions.to_vec();
    }
    /// Rebuild only if needed, returning whether a rebuild occurred.
    pub fn rebuild_if_needed(&mut self, positions: &[[f64; 3]]) -> bool {
        if self.needs_rebuild(positions) {
            self.rebuild(positions);
            true
        } else {
            false
        }
    }
    /// Compute neighbour lists using the spatial hash (O(N) expected).
    fn compute_neighbors(positions: &[[f64; 3]], r: f64) -> Vec<Vec<usize>> {
        SpatialHash3D::all_neighbors(positions, r)
    }
}
impl VerletList {
    /// Estimate the optimal rebuild frequency (steps between rebuilds).
    ///
    /// Given the typical particle displacement per step `v_max * dt`, the list
    /// is safe for `N_rebuild` steps before any particle might cross the skin
    /// distance `r_skin / 2`.  The formula is:
    ///   `N = floor(r_skin / (2 * v_max * dt))`
    /// clamped to at least 1.
    ///
    /// # Arguments
    /// * `v_max` – Maximum expected particle speed \[m/s\].
    /// * `dt`    – Time step size \[s\].
    pub fn compute_rebuild_frequency(&self, v_max: f64, dt: f64) -> usize {
        if v_max <= 0.0 || dt <= 0.0 {
            return 1;
        }
        let displacement_per_step = v_max * dt;
        let safe_steps = (self.r_skin / (2.0 * displacement_per_step)).floor() as usize;
        safe_steps.max(1)
    }
    /// Compute the maximum particle displacement since the last rebuild.
    ///
    /// Compares `positions` against the stored reference positions and returns
    /// the largest displacement magnitude.  Returns 0.0 for empty lists.
    ///
    /// # Arguments
    /// * `positions` – Current particle positions (must be same length as at build time).
    pub fn compute_displacement_since_build(&self, positions: &[[f64; 3]]) -> f64 {
        if self.ref_positions.is_empty() {
            return 0.0;
        }
        positions
            .iter()
            .zip(self.ref_positions.iter())
            .map(|(cur, ref_pos)| {
                let d2 = (cur[0] - ref_pos[0]).powi(2)
                    + (cur[1] - ref_pos[1]).powi(2)
                    + (cur[2] - ref_pos[2]).powi(2);
                d2.sqrt()
            })
            .fold(0.0_f64, f64::max)
    }
}
/// Simple octree-based neighbor search for SPH.
///
/// Not as efficient as LinkedCellList for large simulations, but useful for
/// non-uniform particle distributions.  Builds a balanced octree and queries
/// it for particles within radius `h`.
pub struct OctreeNeighborSearch {
    /// Stored positions (copy).
    pub(super) positions: Vec<[f64; 3]>,
    /// Smoothing radius.
    pub(super) h: f64,
    /// Root of the octree.
    pub(super) root: Option<Box<OctreeNode>>,
}
impl OctreeNeighborSearch {
    /// Build an octree from positions with the given smoothing radius.
    pub fn build(positions: &[[f64; 3]], h: f64) -> Self {
        let max_leaf_size = 8;
        let n = positions.len();
        if n == 0 {
            return Self {
                positions: Vec::new(),
                h,
                root: None,
            };
        }
        let mut lo = positions[0];
        let mut hi = positions[0];
        for p in positions {
            for k in 0..3 {
                if p[k] < lo[k] {
                    lo[k] = p[k];
                }
                if p[k] > hi[k] {
                    hi[k] = p[k];
                }
            }
        }
        for k in 0..3 {
            lo[k] -= h;
            hi[k] += h;
        }
        let indices: Vec<usize> = (0..n).collect();
        let root = Self::build_node(positions, &indices, lo, hi, max_leaf_size);
        Self {
            positions: positions.to_vec(),
            h,
            root: Some(root),
        }
    }
    fn build_node(
        positions: &[[f64; 3]],
        indices: &[usize],
        lo: [f64; 3],
        hi: [f64; 3],
        max_leaf: usize,
    ) -> Box<OctreeNode> {
        if indices.len() <= max_leaf {
            return Box::new(OctreeNode::Leaf {
                indices: indices.to_vec(),
                bounds_lo: lo,
                bounds_hi: hi,
            });
        }
        let mid = [
            (lo[0] + hi[0]) * 0.5,
            (lo[1] + hi[1]) * 0.5,
            (lo[2] + hi[2]) * 0.5,
        ];
        let mut octants: [Vec<usize>; 8] = Default::default();
        for &i in indices {
            let p = positions[i];
            let ix = if p[0] >= mid[0] { 1 } else { 0 };
            let iy = if p[1] >= mid[1] { 1 } else { 0 };
            let iz = if p[2] >= mid[2] { 1 } else { 0 };
            octants[ix * 4 + iy * 2 + iz].push(i);
        }
        let mut children = Vec::new();
        for (oct_idx, oct_indices) in octants.iter().enumerate() {
            if oct_indices.is_empty() {
                continue;
            }
            let ix = (oct_idx >> 2) & 1;
            let iy = (oct_idx >> 1) & 1;
            let iz = oct_idx & 1;
            let child_lo = [
                if ix == 0 { lo[0] } else { mid[0] },
                if iy == 0 { lo[1] } else { mid[1] },
                if iz == 0 { lo[2] } else { mid[2] },
            ];
            let child_hi = [
                if ix == 0 { mid[0] } else { hi[0] },
                if iy == 0 { mid[1] } else { hi[1] },
                if iz == 0 { mid[2] } else { hi[2] },
            ];
            children.push(*Self::build_node(
                positions,
                oct_indices,
                child_lo,
                child_hi,
                max_leaf,
            ));
        }
        Box::new(OctreeNode::Internal {
            bounds_lo: lo,
            bounds_hi: hi,
            children,
        })
    }
    /// Query all particle indices within `h` of `query_pos`.
    pub fn query(&self, query_pos: [f64; 3]) -> Vec<usize> {
        let Some(ref root) = self.root else {
            return Vec::new();
        };
        let mut result = Vec::new();
        Self::query_node(root, &self.positions, query_pos, self.h, &mut result);
        result
    }
    fn query_node(
        node: &OctreeNode,
        positions: &[[f64; 3]],
        query_pos: [f64; 3],
        h: f64,
        result: &mut Vec<usize>,
    ) {
        match node {
            OctreeNode::Leaf {
                indices,
                bounds_lo,
                bounds_hi,
            } => {
                if !Self::aabb_within_radius(*bounds_lo, *bounds_hi, query_pos, h) {
                    return;
                }
                let h2 = h * h;
                for &i in indices {
                    let p = positions[i];
                    let d2 = (p[0] - query_pos[0]).powi(2)
                        + (p[1] - query_pos[1]).powi(2)
                        + (p[2] - query_pos[2]).powi(2);
                    if d2 <= h2 {
                        result.push(i);
                    }
                }
            }
            OctreeNode::Internal {
                bounds_lo,
                bounds_hi,
                children,
            } => {
                if !Self::aabb_within_radius(*bounds_lo, *bounds_hi, query_pos, h) {
                    return;
                }
                for child in children {
                    Self::query_node(child, positions, query_pos, h, result);
                }
            }
        }
    }
    /// Returns true if the AABB could contain points within `h` of `query_pos`.
    fn aabb_within_radius(lo: [f64; 3], hi: [f64; 3], q: [f64; 3], h: f64) -> bool {
        let mut d2 = 0.0_f64;
        for k in 0..3 {
            if q[k] < lo[k] {
                d2 += (lo[k] - q[k]).powi(2);
            } else if q[k] > hi[k] {
                d2 += (q[k] - hi[k]).powi(2);
            }
        }
        d2 <= h * h
    }
    /// Build per-particle neighbor lists (excluding self) using the octree.
    pub fn all_neighbors(&self) -> Vec<Vec<usize>> {
        self.positions
            .iter()
            .enumerate()
            .map(|(i, &pos)| {
                let mut nb: Vec<usize> = self.query(pos).into_iter().filter(|&j| j != i).collect();
                nb.sort_unstable();
                nb
            })
            .collect()
    }
}
/// Spatial hash grid for fast neighbor queries.
#[derive(Debug, Clone)]
pub struct SpatialHash {
    /// Size of each hash cell.
    pub(super) cell_size: f64,
    /// Map from cell coordinates to particle indices.
    pub(super) cells: HashMap<CellKey, Vec<usize>>,
}
impl SpatialHash {
    /// Create a new, empty spatial hash with the given cell size.
    pub fn new(cell_size: f64) -> Self {
        Self {
            cell_size,
            cells: HashMap::new(),
        }
    }
    /// Build the spatial hash from a set of particle positions.
    pub fn build(positions: &[Vec3], cell_size: f64) -> Self {
        let mut hash = Self::new(cell_size);
        for (i, pos) in positions.iter().enumerate() {
            let key = hash.cell_key(pos);
            hash.cells.entry(key).or_default().push(i);
        }
        hash
    }
    /// Query all particle indices within `radius` of `position`.
    pub fn query_neighbors(&self, positions: &[Vec3], position: &Vec3, radius: f64) -> Vec<usize> {
        let r2 = radius * radius;
        let min_key = self.cell_key(&Vec3::new(
            position.x - radius,
            position.y - radius,
            position.z - radius,
        ));
        let max_key = self.cell_key(&Vec3::new(
            position.x + radius,
            position.y + radius,
            position.z + radius,
        ));
        let mut result = Vec::new();
        for cx in min_key.0..=max_key.0 {
            for cy in min_key.1..=max_key.1 {
                for cz in min_key.2..=max_key.2 {
                    if let Some(indices) = self.cells.get(&(cx, cy, cz)) {
                        for &idx in indices {
                            let diff = positions[idx] - position;
                            if diff.norm_squared() <= r2 {
                                result.push(idx);
                            }
                        }
                    }
                }
            }
        }
        result
    }
    /// Find neighbors for every particle. Returns `neighbors[i]` = list of
    /// neighbor indices for particle `i`.
    pub fn find_all_neighbors(positions: &[Vec3], radius: f64) -> Vec<Vec<usize>> {
        let cell_size = radius;
        let hash = Self::build(positions, cell_size);
        let mut all_neighbors = Vec::with_capacity(positions.len());
        for (i, pos) in positions.iter().enumerate() {
            let nbrs: Vec<usize> = hash
                .query_neighbors(positions, pos, radius)
                .into_iter()
                .filter(|&j| j != i)
                .collect();
            all_neighbors.push(nbrs);
        }
        all_neighbors
    }
    /// Compute the cell key for a position.
    fn cell_key(&self, pos: &Vec3) -> CellKey {
        (
            (pos.x / self.cell_size).floor() as i32,
            (pos.y / self.cell_size).floor() as i32,
            (pos.z / self.cell_size).floor() as i32,
        )
    }
}
/// Iterator over unique ordered neighbor pairs `(i, j)` where `i < j`.
///
/// Avoids double-counting: each pair appears exactly once, which is the
/// natural form for computing pairwise forces or energies.
///
/// # Example
/// ```no_run
/// use oxiphysics_sph::neighbor::NeighborPairIterator;
///
/// let lists: Vec<Vec<usize>> = vec![
///     vec![1, 2],
///     vec![0, 2],
///     vec![0, 1],
/// ];
/// let pairs: Vec<(usize, usize)> = NeighborPairIterator::new(&lists).collect();
/// // Should give (0,1), (0,2), (1,2) — each once.
/// assert_eq!(pairs.len(), 3);
/// ```
pub struct NeighborPairIterator<'a> {
    /// Reference to per-particle neighbor lists.
    pub(super) lists: &'a [Vec<usize>],
    /// Current particle index.
    pub(super) i: usize,
    /// Current position within `lists[i]`.
    pub(super) j_idx: usize,
}
impl<'a> NeighborPairIterator<'a> {
    /// Create a new iterator over unique (i, j) pairs with i < j.
    pub fn new(lists: &'a [Vec<usize>]) -> Self {
        Self {
            lists,
            i: 0,
            j_idx: 0,
        }
    }
}
/// Spatial hash with periodic boundary conditions (PBC).
///
/// Wraps particle positions into a box `[0, L_x) × [0, L_y) × [0, L_z)` before
/// hashing and correctly finds neighbors that are close under the minimum-image
/// convention.
///
/// # Example
/// ```no_run
/// use oxiphysics_sph::neighbor::PeriodicSpatialHash;
///
/// let box_size = [1.0_f64; 3];
/// let positions = vec![
///     [0.02, 0.0, 0.0],   // near the lower wall
///     [0.98, 0.0, 0.0],   // near the upper wall — close under PBC
///     [0.5,  0.0, 0.0],   // far from both
/// ];
/// let psh = PeriodicSpatialHash::build(&positions, 0.1, box_size);
/// let nbrs = psh.query_neighbors(0);
/// assert!(nbrs.contains(&1), "PBC neighbor should be found");
/// assert!(!nbrs.contains(&2), "Distant particle should not be found");
/// ```
#[derive(Debug, Clone)]
pub struct PeriodicSpatialHash {
    /// Smoothing radius / cell size.
    pub(super) h: f64,
    /// Box dimensions.
    pub(super) box_size: [f64; 3],
    /// Particle positions (stored as-is, wrapping is applied on query).
    pub(super) positions: Vec<[f64; 3]>,
}
impl PeriodicSpatialHash {
    /// Build the hash from `positions` in a periodic box of size `box_size`.
    pub fn build(positions: &[[f64; 3]], h: f64, box_size: [f64; 3]) -> Self {
        Self {
            h,
            box_size,
            positions: positions.to_vec(),
        }
    }
    /// Return sorted neighbor indices (excluding self) of particle `i` under PBC.
    ///
    /// Queries the 3×3×3 neighbourhood around the query cell **and** the
    /// corresponding periodically-wrapped cells when the query position is close
    /// to a box boundary.  This ensures neighbours across the periodic boundary
    /// are never missed.
    pub fn query_neighbors(&self, i: usize) -> Vec<usize> {
        let pos = self.positions[i];
        let h2 = self.h * self.h;
        let mut result = Vec::new();
        for j in 0..self.positions.len() {
            if j == i {
                continue;
            }
            let rij = self.min_image(pos, self.positions[j]);
            let d2 = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];
            if d2 <= h2 {
                result.push(j);
            }
        }
        result.sort_unstable();
        result
    }
    /// Compute all-particle neighbor lists under PBC.
    pub fn all_neighbors(positions: &[[f64; 3]], h: f64, box_size: [f64; 3]) -> Vec<Vec<usize>> {
        let psh = Self::build(positions, h, box_size);
        (0..positions.len())
            .map(|i| psh.query_neighbors(i))
            .collect()
    }
    /// Minimum-image separation vector r_i − r_j under PBC.
    #[inline]
    pub fn min_image(&self, ri: [f64; 3], rj: [f64; 3]) -> [f64; 3] {
        let mut dr = [0.0_f64; 3];
        for k in 0..3 {
            let mut d = ri[k] - rj[k];
            let l = self.box_size[k];
            d -= (d / l).round() * l;
            dr[k] = d;
        }
        dr
    }
}
/// Z-order (Morton) curve utilities for improving memory-access locality.
///
/// Particles sorted by their Z-order index tend to be spatially close in memory,
/// greatly improving cache performance in SPH force loops.
pub struct ZOrderCurve;
impl ZOrderCurve {
    /// Compute the 3D Morton code for integer coordinates `(x, y, z)`.
    ///
    /// Interleaves the bits of the three 21-bit input values into a 63-bit code.
    /// Input values should fit in 21 bits (< 2^21 ≈ 2 million).
    pub fn morton_code(x: u32, y: u32, z: u32) -> u64 {
        Self::spread_bits(x as u64)
            | (Self::spread_bits(y as u64) << 1)
            | (Self::spread_bits(z as u64) << 2)
    }
    /// Compute the Morton code for a floating-point position in the unit cube `[0, 1)³`.
    ///
    /// Maps each coordinate to a 21-bit integer and then interleaves.
    pub fn morton_f64(pos: [f64; 3], lo: [f64; 3], hi: [f64; 3]) -> u64 {
        let scale = (1u32 << 21) as f64;
        let ix = (((pos[0] - lo[0]) / (hi[0] - lo[0])).clamp(0.0, 1.0 - 1e-10) * scale) as u32;
        let iy = (((pos[1] - lo[1]) / (hi[1] - lo[1])).clamp(0.0, 1.0 - 1e-10) * scale) as u32;
        let iz = (((pos[2] - lo[2]) / (hi[2] - lo[2])).clamp(0.0, 1.0 - 1e-10) * scale) as u32;
        Self::morton_code(ix, iy, iz)
    }
    /// Return a sorted permutation array: `perm[k]` = original index of the k-th
    /// particle in Z-order.
    pub fn sort_by_morton(positions: &[[f64; 3]], lo: [f64; 3], hi: [f64; 3]) -> Vec<usize> {
        let mut keys: Vec<(u64, usize)> = positions
            .iter()
            .enumerate()
            .map(|(i, &pos)| (Self::morton_f64(pos, lo, hi), i))
            .collect();
        keys.sort_unstable_by_key(|&(k, _)| k);
        keys.into_iter().map(|(_, i)| i).collect()
    }
    /// Spread bits of a 21-bit integer to non-adjacent positions (3-bit stride).
    fn spread_bits(mut v: u64) -> u64 {
        v &= 0x001f_ffff;
        v = (v | (v << 32)) & 0x001f_0000_0000_ffff;
        v = (v | (v << 16)) & 0x001f_0000_ffff_0000_u64.wrapping_add(0x0000_ffff);
        v = (v | (v << 16)) & 0xffff_0000_ffff;
        v = (v | (v << 8)) & 0x00ff_00ff_00ff;
        v = (v | (v << 4)) & 0x0f0f_0f0f_0f0f;
        v = (v | (v << 2)) & 0x3333_3333_3333;
        v = (v | (v << 1)) & 0x5555_5555_5555;
        v
    }
}
/// GPU-ready flat neighbor list using two parallel arrays:
/// one for the offsets and one for the neighbor IDs.
///
/// For particle `i`:
/// - neighbors: `flat_neighbors[offsets[i\]..offsets[i+1]]`
///
/// This layout is ideal for parallel GPU kernels because each particle
/// reads a contiguous slice.
#[derive(Debug, Clone)]
pub struct FlatNeighborList {
    /// Packed neighbor IDs.
    pub flat_neighbors: Vec<usize>,
    /// `offsets[i]` = start index in `flat_neighbors` for particle `i`.
    /// Length = n_particles + 1.
    pub offsets: Vec<usize>,
    /// Number of particles.
    pub n_particles: usize,
}
impl FlatNeighborList {
    /// Build a flat neighbor list from per-particle neighbor lists.
    pub fn from_neighbor_lists(neighbor_lists: &[Vec<usize>]) -> Self {
        let n = neighbor_lists.len();
        let mut offsets = Vec::with_capacity(n + 1);
        let mut flat_neighbors = Vec::new();
        offsets.push(0usize);
        for nbrs in neighbor_lists {
            for &nb in nbrs {
                flat_neighbors.push(nb);
            }
            offsets.push(flat_neighbors.len());
        }
        Self {
            flat_neighbors,
            offsets,
            n_particles: n,
        }
    }
    /// Build from positions and smoothing radius, using LinkedCellList internally.
    pub fn build(positions: &[[f64; 3]], h: f64) -> Self {
        let neighbor_lists = LinkedCellList::all_neighbors(positions, h);
        Self::from_neighbor_lists(&neighbor_lists)
    }
    /// Get the neighbor slice for particle `i`.
    pub fn neighbors_of(&self, i: usize) -> &[usize] {
        let start = self.offsets[i];
        let end = self.offsets[i + 1];
        &self.flat_neighbors[start..end]
    }
    /// Total number of neighbor pairs (sum of all neighbor counts).
    pub fn total_pairs(&self) -> usize {
        self.flat_neighbors.len()
    }
    /// Average number of neighbors per particle.
    pub fn avg_neighbors(&self) -> f64 {
        if self.n_particles == 0 {
            return 0.0;
        }
        self.total_pairs() as f64 / self.n_particles as f64
    }
}
/// Builder for per-particle neighbor lists with a note on parallel readiness.
///
/// The implementation here is **serial** but the data structures are designed
/// to be easily parallelised (e.g. via `rayon::par_iter`) because each
/// particle's neighbor list is independent of others.
///
/// # Parallel Extension
/// Replace the inner `map` with `par_iter().map(...)` from rayon to get
/// an O(N / num_threads) build time.
pub struct ParallelReadyNeighborBuilder {
    /// Cut-off radius for neighbor interactions.
    pub r_cut: f64,
    /// Skin radius for Verlet list.
    pub r_skin: f64,
}
impl ParallelReadyNeighborBuilder {
    /// Create a new builder.
    pub fn new(r_cut: f64, r_skin: f64) -> Self {
        Self { r_cut, r_skin }
    }
    /// Build per-particle neighbor lists using `LinkedCellList` (serial).
    ///
    /// Each entry `result[i]` contains the indices of all particles within
    /// `r_cut + r_skin` of particle `i`, excluding self.
    pub fn build_serial(&self, positions: &[[f64; 3]]) -> Vec<Vec<usize>> {
        LinkedCellList::all_neighbors(positions, self.r_cut + self.r_skin)
    }
    /// Build and count the total number of neighbor pairs.
    pub fn total_pairs(&self, positions: &[[f64; 3]]) -> usize {
        let nls = self.build_serial(positions);
        nls.iter().map(|v| v.len()).sum::<usize>() / 2
    }
    /// Compute neighbor count statistics.
    pub fn statistics(&self, positions: &[[f64; 3]]) -> NeighborStats {
        let nls = self.build_serial(positions);
        NeighborStats::compute(&nls)
    }
}
/// Occupancy statistics for a spatial hash grid.
///
/// Useful for diagnosing load-balance issues and choosing optimal cell sizes.
#[derive(Debug, Clone)]
pub struct HashGridStats {
    /// Total number of cells (including empty ones).
    pub total_cells: usize,
    /// Number of non-empty cells.
    pub occupied_cells: usize,
    /// Maximum particles in any single cell.
    pub max_cell_occupancy: usize,
    /// Average particles per non-empty cell.
    pub avg_cell_occupancy: f64,
    /// Occupancy variance (over non-empty cells).
    pub cell_occupancy_variance: f64,
}
impl HashGridStats {
    /// Compute statistics from a `LinkedCellList`.
    pub fn from_linked_cell(lcl: &LinkedCellList, positions: &[[f64; 3]]) -> Self {
        let n = positions.len();
        if n == 0 {
            return Self {
                total_cells: 0,
                occupied_cells: 0,
                max_cell_occupancy: 0,
                avg_cell_occupancy: 0.0,
                cell_occupancy_variance: 0.0,
            };
        }
        let total_cells = lcl.num_cells();
        let mut occupancies: Vec<usize> = Vec::new();
        let mut cell_map: std::collections::HashMap<(i32, i32, i32), usize> =
            std::collections::HashMap::new();
        let cell_size = 1.0_f64;
        let h_guess = if n >= 2 {
            let d = (positions[0][0] - positions[1][0]).hypot(positions[0][1] - positions[1][1]);
            d.max(1e-9)
        } else {
            1e-9
        };
        for pos in positions {
            let key = (
                (pos[0] / h_guess).floor() as i32,
                (pos[1] / h_guess).floor() as i32,
                (pos[2] / h_guess).floor() as i32,
            );
            *cell_map.entry(key).or_insert(0) += 1;
        }
        let _ = cell_size;
        for &count in cell_map.values() {
            occupancies.push(count);
        }
        let occupied_cells = occupancies.len();
        let max_occ = occupancies.iter().cloned().max().unwrap_or(0);
        let avg = if occupied_cells > 0 {
            occupancies.iter().sum::<usize>() as f64 / occupied_cells as f64
        } else {
            0.0
        };
        let var = if occupied_cells > 1 {
            occupancies
                .iter()
                .map(|&c| {
                    let d = c as f64 - avg;
                    d * d
                })
                .sum::<f64>()
                / (occupied_cells - 1) as f64
        } else {
            0.0
        };
        Self {
            total_cells,
            occupied_cells,
            max_cell_occupancy: max_occ,
            avg_cell_occupancy: avg,
            cell_occupancy_variance: var,
        }
    }
    /// Compute grid statistics directly from positions and a cell size.
    pub fn from_positions(positions: &[[f64; 3]], cell_size: f64) -> Self {
        if positions.is_empty() {
            return Self {
                total_cells: 0,
                occupied_cells: 0,
                max_cell_occupancy: 0,
                avg_cell_occupancy: 0.0,
                cell_occupancy_variance: 0.0,
            };
        }
        let mut cell_map: std::collections::HashMap<(i32, i32, i32), usize> =
            std::collections::HashMap::new();
        for pos in positions {
            let key = (
                (pos[0] / cell_size).floor() as i32,
                (pos[1] / cell_size).floor() as i32,
                (pos[2] / cell_size).floor() as i32,
            );
            *cell_map.entry(key).or_insert(0) += 1;
        }
        let occupied_cells = cell_map.len();
        let max_occ = cell_map.values().cloned().max().unwrap_or(0);
        let avg = if occupied_cells > 0 {
            positions.len() as f64 / occupied_cells as f64
        } else {
            0.0
        };
        let var = if occupied_cells > 1 {
            cell_map
                .values()
                .map(|&c| {
                    let d = c as f64 - avg;
                    d * d
                })
                .sum::<f64>()
                / (occupied_cells - 1) as f64
        } else {
            0.0
        };
        Self {
            total_cells: occupied_cells,
            occupied_cells,
            max_cell_occupancy: max_occ,
            avg_cell_occupancy: avg,
            cell_occupancy_variance: var,
        }
    }
}
