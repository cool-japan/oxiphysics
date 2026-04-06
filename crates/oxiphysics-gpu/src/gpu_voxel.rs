// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! GPU voxel grid operations (CPU mock implementation).
//!
//! Provides a binary occupancy grid and a Signed Distance Field (SDF) built
//! from it, plus morphological dilation/erosion and a marching-cubes triangle
//! count estimate.  All routines are CPU references for the GPU backend.

/// A 3-D voxel grid storing binary occupancy and an SDF.
///
/// Dimensions are `[nx, ny, nz]` with voxel size `voxel_size`.
/// Data is stored in Z-major order: `index = x + nx*(y + ny*z)`.
#[derive(Debug, Clone)]
pub struct GpuVoxelGrid {
    /// Number of voxels along X.
    pub nx: usize,
    /// Number of voxels along Y.
    pub ny: usize,
    /// Number of voxels along Z.
    pub nz: usize,
    /// Physical size of each voxel (assumed isotropic).
    pub voxel_size: f64,
    /// Binary occupancy: `true` = filled.
    pub occupancy: Vec<bool>,
    /// Signed Distance Field values (positive = outside, negative = inside).
    pub sdf: Vec<f64>,
}

impl GpuVoxelGrid {
    /// Allocate an empty voxel grid of size `nx × ny × nz`.
    pub fn new(nx: usize, ny: usize, nz: usize, voxel_size: f64) -> Self {
        let n = nx * ny * nz;
        Self {
            nx,
            ny,
            nz,
            voxel_size,
            occupancy: vec![false; n],
            sdf: vec![f64::INFINITY; n],
        }
    }

    /// Total number of voxels.
    pub fn len(&self) -> usize {
        self.nx * self.ny * self.nz
    }

    /// Returns `true` when the grid contains no voxels.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Flatten `(x, y, z)` grid coordinates to a linear index.
    ///
    /// # Panics
    /// In debug builds, panics if any coordinate is out of range.
    pub fn index(&self, x: usize, y: usize, z: usize) -> usize {
        debug_assert!(x < self.nx && y < self.ny && z < self.nz);
        x + self.nx * (y + self.ny * z)
    }

    /// Convert a linear index back to `(x, y, z)` coordinates.
    pub fn coords(&self, idx: usize) -> (usize, usize, usize) {
        let z = idx / (self.nx * self.ny);
        let rem = idx % (self.nx * self.ny);
        let y = rem / self.nx;
        let x = rem % self.nx;
        (x, y, z)
    }

    /// Count occupied voxels.
    pub fn occupied_count(&self) -> usize {
        self.occupancy.iter().filter(|&&v| v).count()
    }
}

// ── GPU-mock operations ───────────────────────────────────────────────────────

/// Rasterize a triangle soup into the voxel grid.
///
/// `triangles` is a flat list of triangle vertices:
/// `[v0x, v0y, v0z, v1x, v1y, v1z, v2x, v2y, v2z, …]`.
///
/// Each voxel whose centre lies inside the axis-aligned bounding box of at
/// least one triangle is marked as occupied.  This is an intentionally
/// conservative (over-filling) approximation of a real GPU rasteriser.
pub fn gpu_voxelize_mesh(grid: &mut GpuVoxelGrid, triangles: &[[f64; 9]]) {
    for tri in triangles {
        // Unpack vertices
        let v0 = [tri[0], tri[1], tri[2]];
        let v1 = [tri[3], tri[4], tri[5]];
        let v2 = [tri[6], tri[7], tri[8]];

        // AABB of the triangle
        let min_x = v0[0].min(v1[0]).min(v2[0]);
        let min_y = v0[1].min(v1[1]).min(v2[1]);
        let min_z = v0[2].min(v1[2]).min(v2[2]);
        let max_x = v0[0].max(v1[0]).max(v2[0]);
        let max_y = v0[1].max(v1[1]).max(v2[1]);
        let max_z = v0[2].max(v1[2]).max(v2[2]);

        let vs = grid.voxel_size;
        let ix0 = ((min_x / vs).floor() as isize).max(0) as usize;
        let iy0 = ((min_y / vs).floor() as isize).max(0) as usize;
        let iz0 = ((min_z / vs).floor() as isize).max(0) as usize;
        let ix1 = ((max_x / vs).ceil() as usize).min(grid.nx.saturating_sub(1));
        let iy1 = ((max_y / vs).ceil() as usize).min(grid.ny.saturating_sub(1));
        let iz1 = ((max_z / vs).ceil() as usize).min(grid.nz.saturating_sub(1));

        for iz in iz0..=iz1 {
            for iy in iy0..=iy1 {
                for ix in ix0..=ix1 {
                    let idx = grid.index(ix, iy, iz);
                    grid.occupancy[idx] = true;
                }
            }
        }
    }
}

/// Compute an approximate SDF from the binary occupancy grid using BFS.
///
/// After the call, `grid.sdf[i]` contains:
/// - `0.0` for occupied voxels,
/// - `k * voxel_size` for voxels at BFS distance `k` from any occupied voxel
///   (positive = outside).
///
/// The SDF is purely unsigned / exterior; interior negation is not computed in
/// this CPU reference implementation.
pub fn gpu_sdf_from_voxels(grid: &mut GpuVoxelGrid) {
    use std::collections::VecDeque;

    let n = grid.len();
    grid.sdf = vec![f64::INFINITY; n];

    let mut queue: VecDeque<usize> = VecDeque::new();

    // Seed: all occupied voxels have distance 0
    for i in 0..n {
        if grid.occupancy[i] {
            grid.sdf[i] = 0.0;
            queue.push_back(i);
        }
    }

    let nx = grid.nx;
    let ny = grid.ny;
    let nz = grid.nz;
    let vs = grid.voxel_size;

    while let Some(idx) = queue.pop_front() {
        let (x, y, z) = grid.coords(idx);
        let current_dist = grid.sdf[idx];

        // 6-connected neighbours
        let neighbours: [(isize, isize, isize); 6] = [
            (1, 0, 0),
            (-1, 0, 0),
            (0, 1, 0),
            (0, -1, 0),
            (0, 0, 1),
            (0, 0, -1),
        ];
        for (dx, dy, dz) in neighbours {
            let nx_ = x as isize + dx;
            let ny_ = y as isize + dy;
            let nz_ = z as isize + dz;
            if nx_ < 0 || ny_ < 0 || nz_ < 0 {
                continue;
            }
            let (nx_, ny_, nz_) = (nx_ as usize, ny_ as usize, nz_ as usize);
            if nx_ >= nx || ny_ >= ny || nz_ >= nz {
                continue;
            }
            let ni = nx_ + nx * (ny_ + ny * nz_);
            let new_dist = current_dist + vs;
            if new_dist < grid.sdf[ni] {
                grid.sdf[ni] = new_dist;
                queue.push_back(ni);
            }
        }
    }
}

/// Morphological dilation: mark a voxel as occupied if any 6-connected
/// neighbour is occupied.
///
/// The operation is applied once (one dilation pass).
pub fn gpu_voxel_dilate(grid: &mut GpuVoxelGrid) {
    let old = grid.occupancy.clone();
    let nx = grid.nx;
    let ny = grid.ny;
    let nz = grid.nz;

    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                if old[x + nx * (y + ny * z)] {
                    continue; // already occupied
                }
                // Check 6-connected neighbours
                let nbrs: [(isize, isize, isize); 6] = [
                    (1, 0, 0),
                    (-1, 0, 0),
                    (0, 1, 0),
                    (0, -1, 0),
                    (0, 0, 1),
                    (0, 0, -1),
                ];
                for (dx, dy, dz) in nbrs {
                    let nx_ = x as isize + dx;
                    let ny_ = y as isize + dy;
                    let nz_ = z as isize + dz;
                    if nx_ < 0 || ny_ < 0 || nz_ < 0 {
                        continue;
                    }
                    let (nx_, ny_, nz_) = (nx_ as usize, ny_ as usize, nz_ as usize);
                    if nx_ >= nx || ny_ >= ny || nz_ >= nz {
                        continue;
                    }
                    if old[nx_ + nx * (ny_ + ny * nz_)] {
                        grid.occupancy[x + nx * (y + ny * z)] = true;
                        break;
                    }
                }
            }
        }
    }
}

/// Morphological erosion: clear a voxel if any 6-connected neighbour is empty.
///
/// The operation is applied once (one erosion pass).
pub fn gpu_voxel_erode(grid: &mut GpuVoxelGrid) {
    let old = grid.occupancy.clone();
    let nx = grid.nx;
    let ny = grid.ny;
    let nz = grid.nz;

    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                if !old[x + nx * (y + ny * z)] {
                    continue; // already empty
                }
                let nbrs: [(isize, isize, isize); 6] = [
                    (1, 0, 0),
                    (-1, 0, 0),
                    (0, 1, 0),
                    (0, -1, 0),
                    (0, 0, 1),
                    (0, 0, -1),
                ];
                for (dx, dy, dz) in nbrs {
                    let nx_ = x as isize + dx;
                    let ny_ = y as isize + dy;
                    let nz_ = z as isize + dz;
                    // Boundary voxels are treated as empty neighbours → erode them
                    let is_empty = if nx_ < 0 || ny_ < 0 || nz_ < 0 {
                        true
                    } else {
                        let (nx_, ny_, nz_) = (nx_ as usize, ny_ as usize, nz_ as usize);
                        if nx_ >= nx || ny_ >= ny || nz_ >= nz {
                            true
                        } else {
                            !old[nx_ + nx * (ny_ + ny * nz_)]
                        }
                    };
                    if is_empty {
                        grid.occupancy[x + nx * (y + ny * z)] = false;
                        break;
                    }
                }
            }
        }
    }
}

/// Count the estimated number of triangles that would be produced by marching
/// cubes on the current occupancy field.
///
/// For each cube (formed by 8 adjacent voxel corners), the occupancy of the 8
/// corners is inspected.  A cube with at least one occupied and one empty
/// corner contributes to the surface.  This function returns the count of such
/// "active" cubes; each active cube typically produces 1-5 triangles; this
/// function returns the number of active cubes as a lower-bound triangle proxy.
pub fn gpu_march_cubes_count(grid: &GpuVoxelGrid) -> usize {
    let nx = grid.nx;
    let ny = grid.ny;
    let nz = grid.nz;

    if nx < 2 || ny < 2 || nz < 2 {
        return 0;
    }

    let mut count = 0usize;
    for z in 0..nz - 1 {
        for y in 0..ny - 1 {
            for x in 0..nx - 1 {
                // 8 corners of the cube
                let corners = [
                    grid.occupancy[grid.index(x, y, z)],
                    grid.occupancy[grid.index(x + 1, y, z)],
                    grid.occupancy[grid.index(x, y + 1, z)],
                    grid.occupancy[grid.index(x + 1, y + 1, z)],
                    grid.occupancy[grid.index(x, y, z + 1)],
                    grid.occupancy[grid.index(x + 1, y, z + 1)],
                    grid.occupancy[grid.index(x, y + 1, z + 1)],
                    grid.occupancy[grid.index(x + 1, y + 1, z + 1)],
                ];
                let any_filled = corners.iter().any(|&v| v);
                let any_empty = corners.iter().any(|&v| !v);
                if any_filled && any_empty {
                    count += 1;
                }
            }
        }
    }
    count
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn small_grid() -> GpuVoxelGrid {
        GpuVoxelGrid::new(4, 4, 4, 1.0)
    }

    // --- GpuVoxelGrid ---

    #[test]
    fn test_grid_len() {
        let g = small_grid();
        assert_eq!(g.len(), 64);
    }

    #[test]
    fn test_grid_is_empty_false() {
        let g = small_grid();
        assert!(!g.is_empty());
    }

    #[test]
    fn test_grid_zero_dim_is_empty() {
        let g = GpuVoxelGrid::new(0, 4, 4, 1.0);
        assert!(g.is_empty());
    }

    #[test]
    fn test_grid_starts_unoccupied() {
        let g = small_grid();
        assert_eq!(g.occupied_count(), 0);
    }

    #[test]
    fn test_grid_sdf_starts_infinity() {
        let g = small_grid();
        assert!(g.sdf.iter().all(|&v| v == f64::INFINITY));
    }

    #[test]
    fn test_index_roundtrip() {
        let g = small_grid();
        for z in 0..4 {
            for y in 0..4 {
                for x in 0..4 {
                    let idx = g.index(x, y, z);
                    let (rx, ry, rz) = g.coords(idx);
                    assert_eq!((rx, ry, rz), (x, y, z));
                }
            }
        }
    }

    #[test]
    fn test_occupied_count() {
        let mut g = small_grid();
        g.occupancy[0] = true;
        g.occupancy[1] = true;
        assert_eq!(g.occupied_count(), 2);
    }

    // --- gpu_voxelize_mesh ---

    #[test]
    fn test_voxelize_single_triangle_marks_voxels() {
        let mut g = GpuVoxelGrid::new(8, 8, 8, 1.0);
        let tri = [[0.5, 0.5, 0.5, 2.5, 0.5, 0.5, 0.5, 2.5, 0.5]];
        gpu_voxelize_mesh(&mut g, &tri);
        assert!(g.occupied_count() > 0);
    }

    #[test]
    fn test_voxelize_empty_triangle_list() {
        let mut g = small_grid();
        gpu_voxelize_mesh(&mut g, &[]);
        assert_eq!(g.occupied_count(), 0);
    }

    #[test]
    fn test_voxelize_two_triangles_more_voxels() {
        let mut g = GpuVoxelGrid::new(10, 10, 10, 1.0);
        let tri1 = [[0.5, 0.5, 0.5, 2.5, 0.5, 0.5, 0.5, 2.5, 0.5]];
        let tri2 = [[5.5, 5.5, 5.5, 7.5, 5.5, 5.5, 5.5, 7.5, 5.5]];
        let mut g1 = g.clone();
        gpu_voxelize_mesh(&mut g1, &tri1);
        let c1 = g1.occupied_count();
        gpu_voxelize_mesh(&mut g, &[tri1[0], tri2[0]]);
        let c2 = g.occupied_count();
        assert!(c2 >= c1);
    }

    #[test]
    fn test_voxelize_triangle_outside_grid_no_panic() {
        let mut g = GpuVoxelGrid::new(4, 4, 4, 1.0);
        let tri = [[
            100.0, 100.0, 100.0, 200.0, 100.0, 100.0, 100.0, 200.0, 100.0,
        ]];
        gpu_voxelize_mesh(&mut g, &tri);
        // No panic, grid unchanged
        assert_eq!(g.occupied_count(), 0);
    }

    // --- gpu_sdf_from_voxels ---

    #[test]
    fn test_sdf_occupied_voxel_is_zero() {
        let mut g = small_grid();
        let _vi = g.index(2, 2, 2);

        g.occupancy[_vi] = true;
        gpu_sdf_from_voxels(&mut g);
        assert!((g.sdf[g.index(2, 2, 2)]).abs() < 1e-12);
    }

    #[test]
    fn test_sdf_neighbour_is_one_voxel_size() {
        let mut g = GpuVoxelGrid::new(5, 5, 5, 1.0);
        let _vi = g.index(2, 2, 2);

        g.occupancy[_vi] = true;
        gpu_sdf_from_voxels(&mut g);
        let d = g.sdf[g.index(3, 2, 2)];
        assert!((d - 1.0).abs() < 1e-12, "d = {d}");
    }

    #[test]
    fn test_sdf_all_unoccupied_stays_infinity() {
        let mut g = small_grid();
        gpu_sdf_from_voxels(&mut g);
        assert!(g.sdf.iter().all(|&v| v == f64::INFINITY));
    }

    #[test]
    fn test_sdf_monotone_with_distance() {
        let mut g = GpuVoxelGrid::new(10, 1, 1, 1.0);
        let _vi = g.index(0, 0, 0);

        g.occupancy[_vi] = true;
        gpu_sdf_from_voxels(&mut g);
        for x in 1..10 {
            let prev = g.sdf[g.index(x - 1, 0, 0)];
            let curr = g.sdf[g.index(x, 0, 0)];
            assert!(
                curr >= prev,
                "SDF not monotone at x={x}: prev={prev} curr={curr}"
            );
        }
    }

    // --- gpu_voxel_dilate ---

    #[test]
    fn test_dilate_expands_single_voxel() {
        let mut g = small_grid();
        let _vi = g.index(2, 2, 2);

        g.occupancy[_vi] = true;
        let before = g.occupied_count();
        gpu_voxel_dilate(&mut g);
        let after = g.occupied_count();
        assert!(after > before);
    }

    #[test]
    fn test_dilate_empty_grid_stays_empty() {
        let mut g = small_grid();
        gpu_voxel_dilate(&mut g);
        assert_eq!(g.occupied_count(), 0);
    }

    #[test]
    fn test_dilate_fully_filled_stays_full() {
        let mut g = small_grid();
        for v in &mut g.occupancy {
            *v = true;
        }
        gpu_voxel_dilate(&mut g);
        assert_eq!(g.occupied_count(), 64);
    }

    // --- gpu_voxel_erode ---

    #[test]
    fn test_erode_single_voxel_disappears() {
        let mut g = small_grid();
        let _vi = g.index(2, 2, 2);

        g.occupancy[_vi] = true;
        gpu_voxel_erode(&mut g);
        // A single voxel always has empty neighbours → eroded away
        assert_eq!(g.occupied_count(), 0);
    }

    #[test]
    fn test_erode_empty_stays_empty() {
        let mut g = small_grid();
        gpu_voxel_erode(&mut g);
        assert_eq!(g.occupied_count(), 0);
    }

    #[test]
    fn test_dilate_then_erode_subset_of_original() {
        let mut g = small_grid();
        // Fill interior 2×2×2
        for z in 1..3 {
            for y in 1..3 {
                for x in 1..3 {
                    let _vi = g.index(x, y, z);

                    g.occupancy[_vi] = true;
                }
            }
        }
        let original_count = g.occupied_count();
        gpu_voxel_dilate(&mut g);
        gpu_voxel_erode(&mut g);
        // After dilate+erode the count should be <= original (boundary loss)
        assert!(g.occupied_count() <= original_count);
    }

    // --- gpu_march_cubes_count ---

    #[test]
    fn test_march_cubes_empty_grid_zero() {
        let g = small_grid();
        assert_eq!(gpu_march_cubes_count(&g), 0);
    }

    #[test]
    fn test_march_cubes_full_grid_zero() {
        let mut g = small_grid();
        for v in &mut g.occupancy {
            *v = true;
        }
        assert_eq!(gpu_march_cubes_count(&g), 0);
    }

    #[test]
    fn test_march_cubes_single_occupied_voxel_nonzero() {
        let mut g = small_grid();
        let _vi = g.index(1, 1, 1);

        g.occupancy[_vi] = true;
        let count = gpu_march_cubes_count(&g);
        assert!(count > 0, "expected > 0 active cubes, got {count}");
    }

    #[test]
    fn test_march_cubes_tiny_grid_no_panic() {
        let g = GpuVoxelGrid::new(1, 1, 1, 1.0);
        assert_eq!(gpu_march_cubes_count(&g), 0);
    }

    #[test]
    fn test_march_cubes_increases_with_more_surface() {
        let mut g1 = GpuVoxelGrid::new(6, 6, 6, 1.0);
        let _vi = g1.index(2, 2, 2);

        g1.occupancy[_vi] = true;
        let c1 = gpu_march_cubes_count(&g1);

        let mut g2 = GpuVoxelGrid::new(6, 6, 6, 1.0);
        for z in 1..4 {
            for y in 1..4 {
                for x in 1..4 {
                    let _vi = g2.index(x, y, z);

                    g2.occupancy[_vi] = true;
                }
            }
        }
        let c2 = gpu_march_cubes_count(&g2);
        assert!(c2 >= c1, "c2={c2} should be >= c1={c1}");
    }
}
