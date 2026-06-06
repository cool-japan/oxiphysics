// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Axis-aligned voxel grid backed by a compact bit vector.
//!
//! Used by the V-HACD decomposer to rasterize triangle soups and compute
//! connected components.

use bitvec::prelude::*;
use std::collections::VecDeque;

/// Axis-aligned voxel grid stored as a bit vector (1 bit per voxel).
///
/// Coordinates are uniform: all three axes use the same `step`.
pub struct VoxelGrid {
    /// Occupancy bits in x-major order: `index = ix + nx*(iy + ny*iz)`.
    pub occupied: BitVec,
    /// World-space minimum corner of the grid.
    pub origin: [f64; 3],
    /// Uniform voxel size (same for all three axes).
    pub step: f64,
    /// Grid dimensions `[nx, ny, nz]`.
    pub dims: [usize; 3],
}

impl VoxelGrid {
    /// Create a new, empty grid covering `[origin, origin + dims * step]`.
    pub fn new(origin: [f64; 3], step: f64, dims: [usize; 3]) -> Self {
        let n = dims[0] * dims[1] * dims[2];
        Self {
            occupied: bitvec![0; n],
            origin,
            step,
            dims,
        }
    }

    /// Linear index from voxel grid coordinates (x-major storage).
    #[inline]
    pub fn index(&self, ix: usize, iy: usize, iz: usize) -> usize {
        ix + self.dims[0] * (iy + self.dims[1] * iz)
    }

    /// Mark voxel `(ix, iy, iz)` as occupied.
    pub fn mark(&mut self, ix: usize, iy: usize, iz: usize) {
        let idx = self.index(ix, iy, iz);
        self.occupied.set(idx, true);
    }

    /// Return true if voxel `(ix, iy, iz)` is occupied.
    pub fn is_occupied(&self, ix: usize, iy: usize, iz: usize) -> bool {
        let idx = self.index(ix, iy, iz);
        self.occupied
            .get(idx)
            .map(|b| *b)
            .expect("voxel index must be in bounds")
    }

    /// World-space centre of voxel `(ix, iy, iz)`.
    pub fn voxel_center(&self, ix: usize, iy: usize, iz: usize) -> [f64; 3] {
        [
            self.origin[0] + (ix as f64 + 0.5) * self.step,
            self.origin[1] + (iy as f64 + 0.5) * self.step,
            self.origin[2] + (iz as f64 + 0.5) * self.step,
        ]
    }

    /// Volume of a single voxel.
    #[inline]
    pub fn voxel_volume(&self) -> f64 {
        self.step * self.step * self.step
    }

    /// Decompose linear index back into grid coordinates.
    #[inline]
    fn coords_from_index(&self, idx: usize) -> (usize, usize, usize) {
        let iz = idx / (self.dims[0] * self.dims[1]);
        let rem = idx % (self.dims[0] * self.dims[1]);
        let iy = rem / self.dims[0];
        let ix = rem % self.dims[0];
        (ix, iy, iz)
    }

    /// 6-connectivity BFS connected components of **occupied** voxels.
    ///
    /// Returns a `Vec` where each element is the list of linear indices belonging
    /// to one connected component.
    pub fn connected_components(&self) -> Vec<Vec<usize>> {
        let n = self.dims[0] * self.dims[1] * self.dims[2];
        let mut visited: BitVec = bitvec![0; n];
        let mut components: Vec<Vec<usize>> = Vec::new();

        for start in 0..n {
            if !self.occupied[start] || visited[start] {
                continue;
            }

            let mut component: Vec<usize> = Vec::new();
            let mut queue: VecDeque<usize> = VecDeque::new();
            queue.push_back(start);
            visited.set(start, true);

            while let Some(idx) = queue.pop_front() {
                component.push(idx);
                let (ix, iy, iz) = self.coords_from_index(idx);

                // 6-connectivity: ±x, ±y, ±z
                let neighbors: [(i64, i64, i64); 6] = [
                    (1, 0, 0),
                    (-1, 0, 0),
                    (0, 1, 0),
                    (0, -1, 0),
                    (0, 0, 1),
                    (0, 0, -1),
                ];
                for (dx, dy, dz) in neighbors {
                    let nx = ix as i64 + dx;
                    let ny = iy as i64 + dy;
                    let nz = iz as i64 + dz;
                    if nx < 0
                        || ny < 0
                        || nz < 0
                        || nx >= self.dims[0] as i64
                        || ny >= self.dims[1] as i64
                        || nz >= self.dims[2] as i64
                    {
                        continue;
                    }
                    let nidx = self.index(nx as usize, ny as usize, nz as usize);
                    if self.occupied[nidx] && !visited[nidx] {
                        visited.set(nidx, true);
                        queue.push_back(nidx);
                    }
                }
            }
            components.push(component);
        }
        components
    }

    /// Total number of occupied voxels.
    pub fn occupied_count(&self) -> usize {
        self.occupied.count_ones()
    }
}
