// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Fast Marching Method for SDF initialization on a uniform 3D grid.

use crate::error::{Result, check_finite, check_len, check_positive};
use std::collections::BinaryHeap;

// ─────────────────────────────────────────────────────────────────────────────
// Fast Marching Method (FMM) for SDF initialization
// ─────────────────────────────────────────────────────────────────────────────

/// State of a grid cell during FMM.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum FmmState {
    /// Final (accepted) value.
    Known,
    /// In the narrow band / heap.
    Trial,
    /// Not yet processed.
    Far,
}

/// Entry in the FMM priority queue.
#[derive(Debug, Clone, Copy)]
struct FmmEntry {
    /// Negative distance (max-heap used as min-heap).
    neg_dist: f64,
    /// Flat grid index.
    idx: usize,
}

impl PartialEq for FmmEntry {
    fn eq(&self, other: &Self) -> bool {
        self.neg_dist == other.neg_dist
    }
}
impl Eq for FmmEntry {}
impl PartialOrd for FmmEntry {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for FmmEntry {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.neg_dist
            .partial_cmp(&other.neg_dist)
            .unwrap_or(std::cmp::Ordering::Equal)
    }
}

/// Fast Marching Method SDF solver on a uniform 3D grid.
///
/// Given initial interface cells (known signed distances), propagates the signed
/// distance function throughout the grid by solving the Eikonal equation
/// `|∇d| = 1` with a first-order upwind scheme.
///
/// The march runs on the distance magnitude `|d|`. Every cell takes the sign of
/// the nearest accepted neighbour it was computed from, so cells enclosed by
/// seeds with negative distance receive negative distances and the others
/// positive ones. Cells that no seed can reach keep the value `f64::MAX`.
///
/// [`FastMarchingMethod::reinitialize`] derives the seeds from an arbitrary
/// level-set function instead of requiring signed distances to be supplied.
#[derive(Debug, Clone)]
pub struct FastMarchingMethod {
    /// Grid size.
    pub nx: usize,
    /// Grid size.
    pub ny: usize,
    /// Grid size.
    pub nz: usize,
    /// Grid spacing.
    pub dx: f64,
    /// Computed signed distances.
    pub dist: Vec<f64>,
    /// FMM state flags.
    state: Vec<FmmState>,
}

impl FastMarchingMethod {
    /// Construct a new FMM solver for a grid of given size and spacing.
    pub fn new(nx: usize, ny: usize, nz: usize, dx: f64) -> Self {
        let n = nx * ny * nz;
        Self {
            nx,
            ny,
            nz,
            dx,
            dist: vec![f64::MAX; n],
            state: vec![FmmState::Far; n],
        }
    }

    #[inline]
    pub(crate) fn flat(&self, ix: usize, iy: usize, iz: usize) -> usize {
        iz * self.ny * self.nx + iy * self.nx + ix
    }

    /// Set known interface cells from (index, signed distance) pairs.
    ///
    /// Negative distances mark cells inside the surface. Pairs whose index lies
    /// outside the grid or whose distance is not finite are ignored.
    pub fn set_known(&mut self, known: &[(usize, f64)]) {
        for &(idx, d) in known {
            if idx < self.dist.len() && d.is_finite() {
                self.dist[idx] = d;
                self.state[idx] = FmmState::Known;
            }
        }
    }

    /// Run the FMM to propagate distances from known cells.
    pub fn run(&mut self) {
        let mut heap: BinaryHeap<FmmEntry> = BinaryHeap::new();

        // Seed with neighbours of known cells
        for iz in 0..self.nz {
            for iy in 0..self.ny {
                for ix in 0..self.nx {
                    let idx = self.flat(ix, iy, iz);
                    if self.state[idx] == FmmState::Known {
                        self.push_neighbours(ix, iy, iz, &mut heap);
                    }
                }
            }
        }

        while let Some(entry) = heap.pop() {
            let cidx = entry.idx;
            if self.state[cidx] == FmmState::Known {
                continue;
            }
            self.state[cidx] = FmmState::Known;
            let iz = cidx / (self.ny * self.nx);
            let rem = cidx % (self.ny * self.nx);
            let iy = rem / self.nx;
            let ix = rem % self.nx;
            self.push_neighbours(ix, iy, iz, &mut heap);
        }
    }

    /// Re-initialise the grid as the signed distance field of the zero level set
    /// of `phi`.
    ///
    /// `phi` holds one value per cell in the layout of the flat grid index (the
    /// index used by [`FastMarchingMethod::set_known`]). It may be any scalar field
    /// whose zero level set is the interface, for example `|x|² - R²`; it does not
    /// have to be a distance function. Negative values are inside, positive values
    /// outside. Previous distances are discarded.
    ///
    /// 1. A cell is an interface cell when `phi` is exactly zero there or a face
    ///    neighbour has the opposite sign (or is exactly zero). Cells with a
    ///    non-finite `phi` are never seeded.
    /// 2. Each interface cell is seeded with the first-order distance estimate
    ///    `|phi| / |∇phi|`, with `∇phi` from central differences (one-sided at the
    ///    grid boundary). The estimate is capped by the distance to the zero
    ///    crossing along a grid axis, located by linear interpolation of `phi`
    ///    between the two cells. Where `phi` is linear between the cells that
    ///    crossing lies on the interface, so its distance is an upper bound of the
    ///    distance to the interface; it also serves as the fallback where the
    ///    gradient vanishes.
    /// 3. The Eikonal equation is marched on the distance magnitude from all seeds,
    ///    and every distance receives the sign of `phi` in its own cell.
    ///
    /// The result does not change when `phi` is multiplied by a positive constant.
    /// A field without a sign change has no interface: nothing is seeded and every
    /// cell keeps the value `f64::MAX`.
    ///
    /// # Errors
    ///
    /// * [`Error::LengthMismatch`](crate::Error::LengthMismatch) if `phi.len()` is
    ///   not `nx * ny * nz`.
    /// * [`Error::OutOfRange`](crate::Error::OutOfRange) or
    ///   [`Error::General`](crate::Error::General) if the grid spacing is not a
    ///   positive finite number.
    pub fn reinitialize(&mut self, phi: &[f64]) -> Result<()> {
        check_len(self.dist.len(), phi.len())?;
        check_finite("dx", self.dx)?;
        check_positive("dx", self.dx)?;

        self.dist.fill(f64::MAX);
        self.state.fill(FmmState::Far);

        let mut seeds: Vec<(usize, f64)> = Vec::new();
        for iz in 0..self.nz {
            for iy in 0..self.ny {
                for ix in 0..self.nx {
                    if let Some(distance) = self.interface_seed_distance(phi, ix, iy, iz) {
                        seeds.push((self.flat(ix, iy, iz), distance));
                    }
                }
            }
        }
        self.set_known(&seeds);
        self.run();

        for (distance, &value) in self.dist.iter_mut().zip(phi) {
            if value < 0.0 && *distance != f64::MAX {
                *distance = -*distance;
            }
        }
        Ok(())
    }

    /// Distance from cell `(ix, iy, iz)` to the zero level set of `phi`, or `None`
    /// when the cell is not next to the interface (or `phi` is not finite there).
    fn interface_seed_distance(&self, phi: &[f64], ix: usize, iy: usize, iz: usize) -> Option<f64> {
        let here = self.flat(ix, iy, iz);
        let value = phi[here];
        if !value.is_finite() {
            return None;
        }
        if value == 0.0 {
            return Some(0.0);
        }
        let dx = self.dx;
        let axes = [
            (ix, self.nx, 1),
            (iy, self.ny, self.nx),
            (iz, self.nz, self.nx * self.ny),
        ];
        let mut gradient_sq = 0.0;
        let mut axis_crossing = f64::INFINITY;
        for (coord, extent, stride) in axes {
            let lower = (coord > 0)
                .then(|| phi[here - stride])
                .filter(|v| v.is_finite());
            let upper = (coord + 1 < extent)
                .then(|| phi[here + stride])
                .filter(|v| v.is_finite());
            let derivative = match (lower, upper) {
                (Some(lo), Some(hi)) => (hi - lo) / (2.0 * dx),
                (None, Some(hi)) => (hi - value) / dx,
                (Some(lo), None) => (value - lo) / dx,
                (None, None) => 0.0,
            };
            gradient_sq += derivative * derivative;
            for neighbour in [lower, upper].into_iter().flatten() {
                // `value` is non-zero, so comparing the signs is exact even where
                // a product of tiny values would underflow.
                let crosses = neighbour == 0.0 || (value < 0.0) != (neighbour < 0.0);
                if crosses {
                    let crossing = dx * value.abs() / (value.abs() + neighbour.abs());
                    axis_crossing = axis_crossing.min(crossing);
                }
            }
        }
        if !axis_crossing.is_finite() {
            return None;
        }
        let first_order = value.abs() / gradient_sq.sqrt();
        Some(first_order.min(axis_crossing))
    }

    fn push_neighbours(
        &mut self,
        ix: usize,
        iy: usize,
        iz: usize,
        heap: &mut BinaryHeap<FmmEntry>,
    ) {
        let neighbors = self.get_neighbors(ix, iy, iz);
        for (nx_i, ny_i, nz_i) in neighbors {
            let nidx = self.flat(nx_i, ny_i, nz_i);
            if self.state[nidx] == FmmState::Known {
                continue;
            }
            let Some((magnitude, inside)) = self.solve_eikonal(nx_i, ny_i, nz_i) else {
                continue;
            };
            if magnitude < self.dist[nidx].abs() {
                self.dist[nidx] = if inside { -magnitude } else { magnitude };
                self.state[nidx] = FmmState::Trial;
                heap.push(FmmEntry {
                    neg_dist: -magnitude,
                    idx: nidx,
                });
            }
        }
    }

    fn get_neighbors(&self, ix: usize, iy: usize, iz: usize) -> Vec<(usize, usize, usize)> {
        let mut ns = Vec::with_capacity(6);
        if ix > 0 {
            ns.push((ix - 1, iy, iz));
        }
        if ix + 1 < self.nx {
            ns.push((ix + 1, iy, iz));
        }
        if iy > 0 {
            ns.push((ix, iy - 1, iz));
        }
        if iy + 1 < self.ny {
            ns.push((ix, iy + 1, iz));
        }
        if iz > 0 {
            ns.push((ix, iy, iz - 1));
        }
        if iz + 1 < self.nz {
            ns.push((ix, iy, iz + 1));
        }
        ns
    }

    /// First-order upwind solution of `|∇d| = 1` at one cell.
    ///
    /// Only accepted (`Known`) neighbours supply upwind data; on every axis the
    /// one with the smaller distance magnitude is used. With the sorted axis
    /// values `a1 <= a2 <= a3` the update value `t` is the largest root of
    /// `Σ ((t - a_i) / dx)² = 1` over the axes that lie below `t`:
    ///
    /// - one axis:    `t = a1 + dx`
    /// - two axes:    `t = (a1 + a2 + sqrt(2 dx² - (a1 - a2)²)) / 2`, used when `a1 + dx > a2`
    /// - three axes:  `t = (s + sqrt(s² - 3 (q - dx²))) / 3` with `s = Σ a_i` and
    ///   `q = Σ a_i²`, used when the two-axis value exceeds `a3`
    ///
    /// Returns the distance magnitude together with whether the nearest upwind
    /// neighbour lies inside the surface (negative distance), or `None` if no
    /// neighbour has been accepted yet.
    fn solve_eikonal(&self, ix: usize, iy: usize, iz: usize) -> Option<(f64, bool)> {
        let dx = self.dx;
        let here = self.flat(ix, iy, iz);
        let axes = [
            (ix, self.nx, 1),
            (iy, self.ny, self.nx),
            (iz, self.nz, self.nx * self.ny),
        ];

        // (magnitude, inside) of the best accepted neighbour on each axis.
        let mut upwind = [(f64::INFINITY, false); 3];
        let mut count = 0;
        for (coord, extent, stride) in axes {
            let mut best: Option<(f64, bool)> = None;
            let mut consider = |idx: usize| {
                if self.state[idx] != FmmState::Known {
                    return;
                }
                let d = self.dist[idx];
                let candidate = (d.abs(), d < 0.0);
                let better = match best {
                    None => true,
                    Some(current) => candidate
                        .0
                        .total_cmp(&current.0)
                        .then(candidate.1.cmp(&current.1))
                        .is_lt(),
                };
                if better {
                    best = Some(candidate);
                }
            };
            if coord > 0 {
                consider(here - stride);
            }
            if coord + 1 < extent {
                consider(here + stride);
            }
            if let Some(found) = best {
                upwind[count] = found;
                count += 1;
            }
        }
        if count == 0 {
            return None;
        }

        // Ascending magnitude; an exact tie prefers the outside neighbour.
        upwind[..count].sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
        let [a1, a2, a3] = [upwind[0].0, upwind[1].0, upwind[2].0];

        let mut t = a1 + dx;
        if count >= 2 && t > a2 {
            let gap = a1 - a2;
            t = 0.5 * (a1 + a2 + (2.0 * dx * dx - gap * gap).max(0.0).sqrt());
            if count == 3 && t > a3 {
                let s = a1 + a2 + a3;
                let q = a1 * a1 + a2 * a2 + a3 * a3;
                t = (s + (s * s - 3.0 * (q - dx * dx)).max(0.0).sqrt()) / 3.0;
            }
        }
        Some((t, upwind[0].1))
    }

    /// Get the distance at grid index (ix, iy, iz).
    pub fn distance_at(&self, ix: usize, iy: usize, iz: usize) -> f64 {
        self.dist[self.flat(ix, iy, iz)]
    }
}

#[cfg(test)]
mod reinitialize_tests;
#[cfg(test)]
mod tests;
