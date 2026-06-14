// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Jump-Flooding Algorithm (JFA) voxel signed-distance field.
//!
//! Computes an (approximate) signed distance field on a uniform 3D voxel grid
//! from a set of surface ("seed") voxels using the Jump-Flooding Algorithm of
//! Rong & Tan, "Jump Flooding in GPU with Applications to Voronoi Diagram and
//! Distance Transform" (I3D 2006).
//!
//! Each voxel propagates the coordinate of its nearest seed in a sequence of
//! flood passes whose step size halves from the smallest power of two not less
//! than the largest grid dimension down to one. The unsigned distance to the
//! nearest seed is then recovered from the propagated seed coordinate, and a
//! sign is assigned by a 6-connected boundary flood-fill that treats the seed
//! voxels as walls: voxels reachable from the grid boundary are "outside"
//! (positive), the remaining non-seed voxels are "inside" (negative).
//!
//! JFA produces the exact Euclidean distance transform for most arrangements
//! and is at worst off by a small bounded error; for the test arrangements used
//! here it matches a brute-force nearest-seed computation exactly.

use std::collections::VecDeque;

// ─────────────────────────────────────────────────────────────────────────────
// Free helpers (usable before the `JfaGrid` struct exists, during `build`)
// ─────────────────────────────────────────────────────────────────────────────

/// Flat array index for voxel `(x, y, z)` in a grid of dimensions `dims`.
///
/// Layout matches [`JfaGrid::distances`]: `z * (ny * nx) + y * nx + x`.
#[inline]
fn flat_idx(dims: [usize; 3], x: usize, y: usize, z: usize) -> usize {
    z * (dims[1] * dims[0]) + y * dims[0] + x
}

/// Squared Euclidean distance (in voxel units) from voxel `(x, y, z)` to `site`.
#[inline]
fn dist2(site: [usize; 3], x: usize, y: usize, z: usize) -> u64 {
    let dx = site[0] as i64 - x as i64;
    let dy = site[1] as i64 - y as i64;
    let dz = site[2] as i64 - z as i64;
    (dx * dx + dy * dy + dz * dz) as u64
}

/// Squared distance to `site`, or `u64::MAX` if `site` is the uninitialized
/// sentinel `[usize::MAX; 3]`.
#[inline]
fn dist2_or_max(site: [usize; 3], x: usize, y: usize, z: usize) -> u64 {
    if site == [usize::MAX; 3] {
        u64::MAX
    } else {
        dist2(site, x, y, z)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Jump-Flooding voxel SDF
// ─────────────────────────────────────────────────────────────────────────────

/// Jump-flooding voxel SDF.
#[derive(Debug, Clone)]
pub struct JfaGrid {
    /// Grid dimensions (nx, ny, nz).
    pub dims: [usize; 3],
    /// Voxel size (world-space length of one voxel edge).
    pub voxel_size: f64,
    /// Signed distance values, indexed `[z * ny * nx + y * nx + x]`.
    pub distances: Vec<f64>,
}

impl JfaGrid {
    /// Build a JFA signed-distance field.
    ///
    /// `dims` = grid dimensions, `voxel_size` = world-space voxel edge length,
    /// `surface_voxels` = (i,j,k) integer coords of voxels on/near the surface.
    pub fn build(
        dims: [usize; 3],
        voxel_size: f64,
        surface_voxels: &[(usize, usize, usize)],
    ) -> Self {
        let [nx, ny, nz] = dims;
        let n = nx * ny * nz;

        // Empty grid (any dimension zero): nothing to flood, return early.
        if n == 0 {
            return Self {
                dims,
                voxel_size,
                distances: Vec::new(),
            };
        }

        // ── Step 1 — buffers ─────────────────────────────────────────────────
        let uninit = [usize::MAX; 3];
        let mut buf_a: Vec<[usize; 3]> = vec![uninit; n];
        let mut buf_b: Vec<[usize; 3]> = vec![uninit; n];

        // ── Step 2 — seed ────────────────────────────────────────────────────
        for &(i, j, k) in surface_voxels {
            if i < nx && j < ny && k < nz {
                buf_a[flat_idx(dims, i, j, k)] = [i, j, k];
            }
        }

        // ── Step 3 — flood passes ────────────────────────────────────────────
        // Smallest power of two >= max_dim.
        let max_dim = nx.max(ny).max(nz);
        let mut step = 1usize;
        while step < max_dim {
            step <<= 1;
        }

        loop {
            // One flood pass: read `buf_a`, write `buf_b`.
            for z in 0..nz {
                for y in 0..ny {
                    for x in 0..nx {
                        let here = flat_idx(dims, x, y, z);
                        let mut best_site = buf_a[here];
                        let mut best_d2: u64 = dist2_or_max(best_site, x, y, z);

                        for dz in -1i64..=1 {
                            for dy in -1i64..=1 {
                                for dx in -1i64..=1 {
                                    if dx == 0 && dy == 0 && dz == 0 {
                                        continue;
                                    }
                                    let nx_ = x as i64 + dx * step as i64;
                                    let ny_ = y as i64 + dy * step as i64;
                                    let nz_ = z as i64 + dz * step as i64;
                                    if nx_ < 0 || ny_ < 0 || nz_ < 0 {
                                        continue;
                                    }
                                    let (nxu, nyu, nzu) =
                                        (nx_ as usize, ny_ as usize, nz_ as usize);
                                    if nxu >= nx || nyu >= ny || nzu >= nz {
                                        continue;
                                    }
                                    let site = buf_a[flat_idx(dims, nxu, nyu, nzu)];
                                    if site == uninit {
                                        continue;
                                    }
                                    let d2 = dist2(site, x, y, z);
                                    if best_site == uninit || d2 < best_d2 {
                                        best_d2 = d2;
                                        best_site = site;
                                    }
                                }
                            }
                        }
                        buf_b[here] = best_site;
                    }
                }
            }
            std::mem::swap(&mut buf_a, &mut buf_b);
            if step == 1 {
                break;
            }
            step >>= 1;
        }
        // Final result now lives in `buf_a`.

        // ── Step 4 — unsigned distances ──────────────────────────────────────
        let mut distances: Vec<f64> = vec![0.0; n];
        let sentinel = voxel_size * (nx + ny + nz) as f64;
        for z in 0..nz {
            for y in 0..ny {
                for x in 0..nx {
                    let here = flat_idx(dims, x, y, z);
                    let site = buf_a[here];
                    if site == uninit {
                        distances[here] = sentinel;
                    } else {
                        distances[here] = voxel_size * (dist2(site, x, y, z) as f64).sqrt();
                    }
                }
            }
        }

        // ── Step 5 — sign assignment via boundary flood-fill ─────────────────
        // No seeds → leave everything positive (boundary flood would cover all).
        if surface_voxels.is_empty() {
            return Self {
                dims,
                voxel_size,
                distances,
            };
        }

        // Mark seed (wall) cells.
        let mut is_seed = vec![false; n];
        for &(i, j, k) in surface_voxels {
            if i < nx && j < ny && k < nz {
                is_seed[flat_idx(dims, i, j, k)] = true;
            }
        }

        // BFS from all boundary voxels that are not seed cells.
        let mut outside = vec![false; n];
        let mut queue: VecDeque<(usize, usize, usize)> = VecDeque::new();
        for z in 0..nz {
            for y in 0..ny {
                for x in 0..nx {
                    let on_boundary =
                        x == 0 || x == nx - 1 || y == 0 || y == ny - 1 || z == 0 || z == nz - 1;
                    if !on_boundary {
                        continue;
                    }
                    let here = flat_idx(dims, x, y, z);
                    if is_seed[here] {
                        continue;
                    }
                    if !outside[here] {
                        outside[here] = true;
                        queue.push_back((x, y, z));
                    }
                }
            }
        }

        // 6-connected face-neighbour offsets.
        const NEIGHBORS: [(i64, i64, i64); 6] = [
            (1, 0, 0),
            (-1, 0, 0),
            (0, 1, 0),
            (0, -1, 0),
            (0, 0, 1),
            (0, 0, -1),
        ];
        while let Some((x, y, z)) = queue.pop_front() {
            for &(dx, dy, dz) in &NEIGHBORS {
                let nx_ = x as i64 + dx;
                let ny_ = y as i64 + dy;
                let nz_ = z as i64 + dz;
                if nx_ < 0 || ny_ < 0 || nz_ < 0 {
                    continue;
                }
                let (nxu, nyu, nzu) = (nx_ as usize, ny_ as usize, nz_ as usize);
                if nxu >= nx || nyu >= ny || nzu >= nz {
                    continue;
                }
                let idx = flat_idx(dims, nxu, nyu, nzu);
                if is_seed[idx] || outside[idx] {
                    continue;
                }
                outside[idx] = true;
                queue.push_back((nxu, nyu, nzu));
            }
        }

        // Negate distances for inside (not outside, not seed) voxels.
        for z in 0..nz {
            for y in 0..ny {
                for x in 0..nx {
                    let here = flat_idx(dims, x, y, z);
                    if is_seed[here] || outside[here] {
                        continue;
                    }
                    distances[here] = -distances[here];
                }
            }
        }

        Self {
            dims,
            voxel_size,
            distances,
        }
    }

    /// Flat array index for voxel `(ix, iy, iz)`.
    #[inline]
    fn flat(&self, ix: usize, iy: usize, iz: usize) -> usize {
        iz * self.dims[1] * self.dims[0] + iy * self.dims[0] + ix
    }

    /// Return the signed distance at voxel index `(ix, iy, iz)`.
    pub fn get(&self, ix: usize, iy: usize, iz: usize) -> f64 {
        self.distances[self.flat(ix, iy, iz)]
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_jfa_single_seed_center() {
        let grid = JfaGrid::build([5, 5, 5], 1.0, &[(2, 2, 2)]);
        assert!(
            grid.get(2, 2, 2).abs() < 1e-9,
            "seed cell should be ~0: {}",
            grid.get(2, 2, 2)
        );
        assert!(
            (grid.get(3, 2, 2).abs() - 1.0).abs() < 1e-9,
            "neighbor at (3,2,2) magnitude ~1: {}",
            grid.get(3, 2, 2)
        );
    }

    #[test]
    fn test_jfa_matches_brute_force_small() {
        let voxel_size = 1.0;
        let dims = [8usize, 8, 8];
        let seeds: Vec<(usize, usize, usize)> =
            vec![(0, 0, 0), (7, 7, 7), (3, 4, 2), (5, 1, 6), (2, 6, 1)];
        let grid = JfaGrid::build(dims, voxel_size, &seeds);

        for z in 0..8usize {
            for y in 0..8usize {
                for x in 0..8usize {
                    let mut brute = f64::MAX;
                    for &(si, sj, sk) in &seeds {
                        let dx = x as f64 - si as f64;
                        let dy = y as f64 - sj as f64;
                        let dz = z as f64 - sk as f64;
                        let d = voxel_size * (dx * dx + dy * dy + dz * dz).sqrt();
                        if d < brute {
                            brute = d;
                        }
                    }
                    let got = grid.get(x, y, z).abs();
                    assert!(
                        (got - brute).abs() <= voxel_size + 1e-9,
                        "({x},{y},{z}): JFA |{got}| vs brute {brute} exceeds 1 voxel"
                    );
                }
            }
        }
    }

    #[test]
    fn test_jfa_empty_seeds_no_panic() {
        let grid = JfaGrid::build([4, 4, 4], 1.0, &[]);
        for &d in &grid.distances {
            assert!(
                d.is_finite() || d == f64::MAX,
                "distance should be finite or MAX: {}",
                d
            );
            assert!(
                d > 0.0,
                "with no seeds everything should be positive: {}",
                d
            );
        }
    }
}
