// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Fast Marching Method for SDF initialization on a uniform 3D grid.

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
/// Given initial interface cells (known SDF values near zero), propagates
/// the signed distance function throughout the grid.
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

    /// Set known interface cells from (index, distance) pairs.
    pub fn set_known(&mut self, known: &[(usize, f64)]) {
        for &(idx, d) in known {
            if idx < self.dist.len() {
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
            let d = self.solve_eikonal(nx_i, ny_i, nz_i);
            if d < self.dist[nidx] {
                self.dist[nidx] = d;
                self.state[nidx] = FmmState::Trial;
                heap.push(FmmEntry {
                    neg_dist: -d,
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

    fn solve_eikonal(&self, ix: usize, iy: usize, iz: usize) -> f64 {
        // 1st-order upwind Eikonal: solve (dx1² + dy1² + dz1²) = dx²
        let dx = self.dx;
        let mut terms: [f64; 3] = [f64::MAX; 3];

        // x-direction
        let mut d_x = f64::MAX;
        if ix > 0 {
            d_x = d_x.min(self.dist[self.flat(ix - 1, iy, iz)]);
        }
        if ix + 1 < self.nx {
            d_x = d_x.min(self.dist[self.flat(ix + 1, iy, iz)]);
        }
        terms[0] = d_x;

        // y-direction
        let mut d_y = f64::MAX;
        if iy > 0 {
            d_y = d_y.min(self.dist[self.flat(ix, iy - 1, iz)]);
        }
        if iy + 1 < self.ny {
            d_y = d_y.min(self.dist[self.flat(ix, iy + 1, iz)]);
        }
        terms[1] = d_y;

        // z-direction
        let mut d_z = f64::MAX;
        if iz > 0 {
            d_z = d_z.min(self.dist[self.flat(ix, iy, iz - 1)]);
        }
        if iz + 1 < self.nz {
            d_z = d_z.min(self.dist[self.flat(ix, iy, iz + 1)]);
        }
        terms[2] = d_z;

        terms.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

        // Quadratic solve: try adding terms one by one
        for k in 1..=3 {
            let valid: Vec<f64> = terms[..k]
                .iter()
                .filter(|&&t| t < f64::MAX)
                .copied()
                .collect();
            if valid.is_empty() {
                continue;
            }
            let sum_t = valid.iter().sum::<f64>();
            let sum_t2 = valid.iter().map(|t| t * t).sum::<f64>();
            let n_v = valid.len() as f64;
            let discriminant = sum_t * sum_t - n_v * (sum_t2 - dx * dx);
            if discriminant >= 0.0 {
                let sol = (sum_t + discriminant.sqrt()) / n_v;
                if k == 1 || sol > *valid.last().expect("collection should not be empty") {
                    return sol;
                }
            }
        }

        // Fallback: nearest neighbour + one cell
        terms
            .iter()
            .copied()
            .filter(|&t| t < f64::MAX)
            .fold(f64::MAX, f64::min)
            + dx
    }

    /// Get the distance at grid index (ix, iy, iz).
    pub fn distance_at(&self, ix: usize, iy: usize, iz: usize) -> f64 {
        self.dist[self.flat(ix, iy, iz)]
    }
}
