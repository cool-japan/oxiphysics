// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Voxel-based collision detection.
//!
//! Provides an SDF (Signed Distance Field) voxel grid, voxelization of triangle
//! meshes, sparse voxel octree, fast-sweeping distance field computation, ray vs
//! voxel grid intersection, sphere vs voxel grid, convex shape vs voxel grid, and
//! voxel-based collision response.
//!
//! # Overview
//!
//! The voxel grid stores per-cell SDF values: negative inside a solid, positive
//! outside, zero on the surface.  The [`SdfVoxelGrid`] struct is the central data
//! structure; helper types provide distance field computation ([`FastSweeping`]),
//! sparse storage ([`SparseVoxelOctree`]), and query primitives
//! ([`VoxelRayQuery`], [`VoxelSphereQuery`], [`VoxelConvexQuery`]).

use std::collections::HashMap;

// ─────────────────────────────────────────────────────────────────────────────
// Basic 3-D integer and float coordinate types
// ─────────────────────────────────────────────────────────────────────────────

/// Integer voxel coordinate `(ix, iy, iz)`.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct VoxelCoord {
    /// X index.
    pub ix: i32,
    /// Y index.
    pub iy: i32,
    /// Z index.
    pub iz: i32,
}

#[allow(dead_code)]
impl VoxelCoord {
    /// Create a new `VoxelCoord`.
    pub fn new(ix: i32, iy: i32, iz: i32) -> Self {
        Self { ix, iy, iz }
    }

    /// Return the Manhattan distance to another voxel coordinate.
    pub fn manhattan(&self, other: &Self) -> i32 {
        (self.ix - other.ix).abs() + (self.iy - other.iy).abs() + (self.iz - other.iz).abs()
    }

    /// Return the Chebyshev (L∞) distance to another voxel coordinate.
    pub fn chebyshev(&self, other: &Self) -> i32 {
        (self.ix - other.ix)
            .abs()
            .max((self.iy - other.iy).abs())
            .max((self.iz - other.iz).abs())
    }

    /// Return the 6-connected (face) neighbours of this voxel.
    pub fn face_neighbours(&self) -> [Self; 6] {
        [
            Self::new(self.ix - 1, self.iy, self.iz),
            Self::new(self.ix + 1, self.iy, self.iz),
            Self::new(self.ix, self.iy - 1, self.iz),
            Self::new(self.ix, self.iy + 1, self.iz),
            Self::new(self.ix, self.iy, self.iz - 1),
            Self::new(self.ix, self.iy, self.iz + 1),
        ]
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// SdfVoxelGrid
// ─────────────────────────────────────────────────────────────────────────────

/// A dense, axis-aligned SDF voxel grid.
///
/// Each voxel stores a signed distance value; negative values are inside the
/// solid, positive values outside, and zero is on the surface.  The world-space
/// position of voxel `(ix, iy, iz)` is:
/// `origin + (ix + 0.5, iy + 0.5, iz + 0.5) * cell_size`.
#[allow(dead_code)]
pub struct SdfVoxelGrid {
    /// World-space origin of the grid (minimum corner).
    pub origin: [f64; 3],
    /// Side length of one cubic voxel in world units.
    pub cell_size: f64,
    /// Number of voxels in the X direction.
    pub nx: usize,
    /// Number of voxels in the Y direction.
    pub ny: usize,
    /// Number of voxels in the Z direction.
    pub nz: usize,
    /// Flat array of SDF values, indexed `ix + nx*(iy + ny*iz)`.
    pub data: Vec<f64>,
}

#[allow(dead_code)]
impl SdfVoxelGrid {
    /// Create a new grid filled with `f64::INFINITY` (everything outside).
    ///
    /// - `origin` — world-space minimum corner
    /// - `cell_size` — side length of each cubic voxel
    /// - `nx`, `ny`, `nz` — grid resolution per axis
    pub fn new(origin: [f64; 3], cell_size: f64, nx: usize, ny: usize, nz: usize) -> Self {
        let data = vec![f64::INFINITY; nx * ny * nz];
        Self {
            origin,
            cell_size,
            nx,
            ny,
            nz,
            data,
        }
    }

    /// Total number of voxels in the grid.
    pub fn len(&self) -> usize {
        self.nx * self.ny * self.nz
    }

    /// Return `true` if the grid contains no voxels.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Convert a `(ix, iy, iz)` triple to a flat index.
    ///
    /// Returns `None` if the coordinate is out of bounds.
    pub fn to_index(&self, ix: usize, iy: usize, iz: usize) -> Option<usize> {
        if ix < self.nx && iy < self.ny && iz < self.nz {
            Some(ix + self.nx * (iy + self.ny * iz))
        } else {
            None
        }
    }

    /// Return the SDF value at `(ix, iy, iz)`, or `None` if out of bounds.
    pub fn get(&self, ix: usize, iy: usize, iz: usize) -> Option<f64> {
        self.to_index(ix, iy, iz).map(|idx| self.data[idx])
    }

    /// Set the SDF value at `(ix, iy, iz)`.
    ///
    /// Does nothing if the coordinate is out of bounds.
    pub fn set(&mut self, ix: usize, iy: usize, iz: usize, value: f64) {
        if let Some(idx) = self.to_index(ix, iy, iz) {
            self.data[idx] = value;
        }
    }

    /// Convert a world-space point to voxel coordinates.
    ///
    /// Returns `None` if the point is outside the grid.
    pub fn world_to_voxel(&self, p: [f64; 3]) -> Option<[usize; 3]> {
        let ix = ((p[0] - self.origin[0]) / self.cell_size).floor() as i64;
        let iy = ((p[1] - self.origin[1]) / self.cell_size).floor() as i64;
        let iz = ((p[2] - self.origin[2]) / self.cell_size).floor() as i64;
        if ix >= 0
            && iy >= 0
            && iz >= 0
            && (ix as usize) < self.nx
            && (iy as usize) < self.ny
            && (iz as usize) < self.nz
        {
            Some([ix as usize, iy as usize, iz as usize])
        } else {
            None
        }
    }

    /// Return the world-space centre of voxel `(ix, iy, iz)`.
    pub fn voxel_center(&self, ix: usize, iy: usize, iz: usize) -> [f64; 3] {
        [
            self.origin[0] + (ix as f64 + 0.5) * self.cell_size,
            self.origin[1] + (iy as f64 + 0.5) * self.cell_size,
            self.origin[2] + (iz as f64 + 0.5) * self.cell_size,
        ]
    }

    /// Trilinearly interpolate the SDF at world-space point `p`.
    ///
    /// Returns `f64::INFINITY` outside the grid.
    pub fn interpolate(&self, p: [f64; 3]) -> f64 {
        let lx = (p[0] - self.origin[0]) / self.cell_size - 0.5;
        let ly = (p[1] - self.origin[1]) / self.cell_size - 0.5;
        let lz = (p[2] - self.origin[2]) / self.cell_size - 0.5;

        let ix0 = lx.floor() as i64;
        let iy0 = ly.floor() as i64;
        let iz0 = lz.floor() as i64;

        let tx = lx - ix0 as f64;
        let ty = ly - iy0 as f64;
        let tz = lz - iz0 as f64;

        let mut acc = 0.0_f64;
        for di in 0..2_i64 {
            for dj in 0..2_i64 {
                for dk in 0..2_i64 {
                    let ix = ix0 + di;
                    let iy = iy0 + dj;
                    let iz = iz0 + dk;
                    if ix < 0
                        || iy < 0
                        || iz < 0
                        || ix as usize >= self.nx
                        || iy as usize >= self.ny
                        || iz as usize >= self.nz
                    {
                        return f64::INFINITY;
                    }
                    let val =
                        self.data[ix as usize + self.nx * (iy as usize + self.ny * iz as usize)];
                    let wx = if di == 0 { 1.0 - tx } else { tx };
                    let wy = if dj == 0 { 1.0 - ty } else { ty };
                    let wz = if dk == 0 { 1.0 - tz } else { tz };
                    acc += val * wx * wy * wz;
                }
            }
        }
        acc
    }

    /// Compute the gradient of the SDF at voxel `(ix, iy, iz)` using central differences.
    ///
    /// Returns a unit-length gradient vector, or `[0, 0, 0]` if on the boundary.
    pub fn gradient(&self, ix: usize, iy: usize, iz: usize) -> [f64; 3] {
        if ix == 0
            || iy == 0
            || iz == 0
            || ix + 1 >= self.nx
            || iy + 1 >= self.ny
            || iz + 1 >= self.nz
        {
            return [0.0; 3];
        }
        let idx_xp = match self.to_index(ix + 1, iy, iz) {
            Some(i) => i,
            None => return [0.0; 3],
        };
        let idx_xm = match self.to_index(ix - 1, iy, iz) {
            Some(i) => i,
            None => return [0.0; 3],
        };
        let idx_yp = match self.to_index(ix, iy + 1, iz) {
            Some(i) => i,
            None => return [0.0; 3],
        };
        let idx_ym = match self.to_index(ix, iy - 1, iz) {
            Some(i) => i,
            None => return [0.0; 3],
        };
        let idx_zp = match self.to_index(ix, iy, iz + 1) {
            Some(i) => i,
            None => return [0.0; 3],
        };
        let idx_zm = match self.to_index(ix, iy, iz - 1) {
            Some(i) => i,
            None => return [0.0; 3],
        };
        let dx = (self.data[idx_xp] - self.data[idx_xm]) / (2.0 * self.cell_size);
        let dy = (self.data[idx_yp] - self.data[idx_ym]) / (2.0 * self.cell_size);
        let dz = (self.data[idx_zp] - self.data[idx_zm]) / (2.0 * self.cell_size);
        let len = (dx * dx + dy * dy + dz * dz).sqrt();
        if len > 1e-14 {
            [dx / len, dy / len, dz / len]
        } else {
            [0.0; 3]
        }
    }

    /// Return the AABB of the grid as `(min, max)` world-space corners.
    pub fn aabb(&self) -> ([f64; 3], [f64; 3]) {
        let min = self.origin;
        let max = [
            self.origin[0] + self.nx as f64 * self.cell_size,
            self.origin[1] + self.ny as f64 * self.cell_size,
            self.origin[2] + self.nz as f64 * self.cell_size,
        ];
        (min, max)
    }

    /// Fill the entire grid with `value`.
    pub fn fill(&mut self, value: f64) {
        self.data.iter_mut().for_each(|v| *v = value);
    }

    /// Count voxels with SDF value strictly less than `threshold` (inside test).
    pub fn count_inside(&self, threshold: f64) -> usize {
        self.data.iter().filter(|&&v| v < threshold).count()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Voxelization of triangle meshes
// ─────────────────────────────────────────────────────────────────────────────

/// Voxelizer for triangle meshes.
///
/// Rasterizes a triangle mesh into an `SdfVoxelGrid` using point-to-triangle
/// distance calculations.
#[allow(dead_code)]
pub struct MeshVoxelizer;

#[allow(dead_code)]
impl MeshVoxelizer {
    /// Compute the squared distance from point `p` to triangle `(a, b, c)`.
    pub fn point_triangle_dist_sq(p: [f64; 3], a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
        let ab = sub3(b, a);
        let ac = sub3(c, a);
        let ap = sub3(p, a);

        let d1 = dot3(ab, ap);
        let d2 = dot3(ac, ap);
        if d1 <= 0.0 && d2 <= 0.0 {
            return dist_sq3(p, a);
        }

        let bp = sub3(p, b);
        let d3 = dot3(ab, bp);
        let d4 = dot3(ac, bp);
        if d3 >= 0.0 && d4 <= d3 {
            return dist_sq3(p, b);
        }

        let cp = sub3(p, c);
        let d5 = dot3(ab, cp);
        let d6 = dot3(ac, cp);
        if d6 >= 0.0 && d5 <= d6 {
            return dist_sq3(p, c);
        }

        let vc = d1 * d4 - d3 * d2;
        if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
            let v = d1 / (d1 - d3);
            let proj = add3(a, scale3(ab, v));
            return dist_sq3(p, proj);
        }

        let vb = d5 * d2 - d1 * d6;
        if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
            let w = d2 / (d2 - d6);
            let proj = add3(a, scale3(ac, w));
            return dist_sq3(p, proj);
        }

        let va = d3 * d6 - d5 * d4;
        let dp = d4 - d3;
        let ep = d5 - d6;
        if va <= 0.0 && dp >= 0.0 && ep >= 0.0 {
            let w = dp / (dp + ep);
            let proj = add3(b, scale3(sub3(c, b), w));
            return dist_sq3(p, proj);
        }

        let denom = 1.0 / (va + vb + vc);
        let v = vb * denom;
        let w = vc * denom;
        let proj = add3(a, add3(scale3(ab, v), scale3(ac, w)));
        dist_sq3(p, proj)
    }

    /// Compute the winding number sign of point `p` relative to triangle `(a, b, c)`.
    ///
    /// Returns `+1.0` if `p` is on the positive-normal side, `-1.0` otherwise.
    pub fn sign_winding(p: [f64; 3], a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
        let n = cross3(sub3(b, a), sub3(c, a));
        let d = dot3(n, sub3(p, a));
        if d >= 0.0 { 1.0 } else { -1.0 }
    }

    /// Voxelize a closed triangle mesh and write SDF values into `grid`.
    ///
    /// `vertices` is a flat list of `[f64; 3]` positions; `triangles` is a list of
    /// index triples `(i0, i1, i2)` into `vertices`.
    pub fn voxelize(grid: &mut SdfVoxelGrid, vertices: &[[f64; 3]], triangles: &[[usize; 3]]) {
        for iz in 0..grid.nz {
            for iy in 0..grid.ny {
                for ix in 0..grid.nx {
                    let center = grid.voxel_center(ix, iy, iz);
                    let mut min_dist_sq = f64::INFINITY;
                    let mut sign = 1.0_f64;
                    for tri in triangles {
                        let a = vertices[tri[0]];
                        let b = vertices[tri[1]];
                        let c = vertices[tri[2]];
                        let d_sq = Self::point_triangle_dist_sq(center, a, b, c);
                        if d_sq < min_dist_sq {
                            min_dist_sq = d_sq;
                            sign = Self::sign_winding(center, a, b, c);
                        }
                    }
                    let sdf = sign * min_dist_sq.sqrt();
                    if let Some(idx) = grid.to_index(ix, iy, iz) {
                        grid.data[idx] = sdf;
                    }
                }
            }
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// FastSweeping — distance field computation
// ─────────────────────────────────────────────────────────────────────────────

/// Fast Sweeping Method for computing the distance field from a set of seed
/// voxels.
///
/// Iterates over the grid in 8 alternating sweep directions, propagating the
/// minimum distance to the nearest seed voxel.
#[allow(dead_code)]
pub struct FastSweeping;

#[allow(dead_code)]
impl FastSweeping {
    /// Compute the unsigned distance field from `seeds` using fast sweeping.
    ///
    /// `seeds` is a list of `(voxel_coord, sdf_value)` pairs that initialise the
    /// computation.  Returns a vector of SDF values of the same length as
    /// `grid.data`.
    pub fn compute(grid: &SdfVoxelGrid, seeds: &[(VoxelCoord, f64)]) -> Vec<f64> {
        let nx = grid.nx;
        let ny = grid.ny;
        let nz = grid.nz;
        let h = grid.cell_size;
        let n = nx * ny * nz;

        // Initialise all values to +INF
        let mut phi = vec![f64::INFINITY; n];

        // Seed the known values
        for &(coord, val) in seeds {
            if coord.ix >= 0
                && coord.iy >= 0
                && coord.iz >= 0
                && (coord.ix as usize) < nx
                && (coord.iy as usize) < ny
                && (coord.iz as usize) < nz
            {
                let idx = coord.ix as usize + nx * (coord.iy as usize + ny * coord.iz as usize);
                phi[idx] = val;
            }
        }

        // 8-direction sweeps
        let orders: [(bool, bool, bool); 8] = [
            (false, false, false),
            (true, false, false),
            (false, true, false),
            (true, true, false),
            (false, false, true),
            (true, false, true),
            (false, true, true),
            (true, true, true),
        ];

        for (rev_x, rev_y, rev_z) in orders {
            let xs: Vec<usize> = if rev_x {
                (0..nx).rev().collect()
            } else {
                (0..nx).collect()
            };
            let ys: Vec<usize> = if rev_y {
                (0..ny).rev().collect()
            } else {
                (0..ny).collect()
            };
            let zs: Vec<usize> = if rev_z {
                (0..nz).rev().collect()
            } else {
                (0..nz).collect()
            };

            for &iz in &zs {
                for &iy in &ys {
                    for &ix in &xs {
                        // Godunov upwind scheme: u = min of face neighbours
                        let a = if ix > 0 {
                            phi[ix - 1 + nx * (iy + ny * iz)]
                        } else {
                            f64::INFINITY
                        }
                        .min(if ix + 1 < nx {
                            phi[ix + 1 + nx * (iy + ny * iz)]
                        } else {
                            f64::INFINITY
                        });

                        let b = if iy > 0 {
                            phi[ix + nx * (iy - 1 + ny * iz)]
                        } else {
                            f64::INFINITY
                        }
                        .min(if iy + 1 < ny {
                            phi[ix + nx * (iy + 1 + ny * iz)]
                        } else {
                            f64::INFINITY
                        });

                        let c = if iz > 0 {
                            phi[ix + nx * (iy + ny * (iz - 1))]
                        } else {
                            f64::INFINITY
                        }
                        .min(if iz + 1 < nz {
                            phi[ix + nx * (iy + ny * (iz + 1))]
                        } else {
                            f64::INFINITY
                        });

                        // Sort a ≤ b ≤ c
                        let mut vals = [a, b, c];
                        vals.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
                        let [a, b, c] = vals;

                        // Eikonal update
                        let u_new = if (a + h - b).abs() > h || b.is_infinite() {
                            a + h
                        } else if (((a + b + h) * (a + b + h) - 2.0 * (a * a + b * b - h * h))
                            .max(0.0)
                            .sqrt()
                            + a
                            + b)
                            / 2.0
                            <= c
                            || c.is_infinite()
                        {
                            (((a + b + h) * (a + b + h) - 2.0 * (a * a + b * b - h * h))
                                .max(0.0)
                                .sqrt()
                                + a
                                + b)
                                / 2.0
                        } else {
                            // All three contribute
                            let rhs = h * h;
                            // (u-a)^2 + (u-b)^2 + (u-c)^2 = h^2
                            // 3u^2 - 2(a+b+c)u + (a^2+b^2+c^2 - h^2) = 0
                            let sum = a + b + c;
                            let sq = a * a + b * b + c * c;
                            let disc = (3.0 * rhs - (sq - sum * sum / 3.0)).max(0.0);
                            sum / 3.0 + disc.sqrt() / 3.0_f64.sqrt()
                        };

                        let idx = ix + nx * (iy + ny * iz);
                        if u_new < phi[idx] {
                            phi[idx] = u_new;
                        }
                    }
                }
            }
        }
        phi
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// SparseVoxelOctree
// ─────────────────────────────────────────────────────────────────────────────

/// A sparse voxel octree for memory-efficient storage of large voxel grids.
///
/// Only occupied (non-default) voxels are stored.  The octree is represented
/// as a flat `HashMap` mapping `VoxelCoord` to `f64` SDF values plus a level
/// counter.
#[allow(dead_code)]
pub struct SparseVoxelOctree {
    /// Occupied voxels and their SDF values.
    pub cells: HashMap<VoxelCoord, f64>,
    /// World-space origin.
    pub origin: [f64; 3],
    /// Base cell size (leaf level).
    pub cell_size: f64,
    /// Default value returned for absent voxels.
    pub default_value: f64,
}

#[allow(dead_code)]
impl SparseVoxelOctree {
    /// Create an empty sparse octree.
    pub fn new(origin: [f64; 3], cell_size: f64, default_value: f64) -> Self {
        Self {
            cells: HashMap::new(),
            origin,
            cell_size,
            default_value,
        }
    }

    /// Insert a voxel SDF value at `coord`.
    pub fn insert(&mut self, coord: VoxelCoord, value: f64) {
        self.cells.insert(coord, value);
    }

    /// Remove a voxel from the octree.
    pub fn remove(&mut self, coord: &VoxelCoord) {
        self.cells.remove(coord);
    }

    /// Query the SDF value at `coord`.
    ///
    /// Returns `self.default_value` if not present.
    pub fn get(&self, coord: VoxelCoord) -> f64 {
        *self.cells.get(&coord).unwrap_or(&self.default_value)
    }

    /// Number of stored voxels.
    pub fn len(&self) -> usize {
        self.cells.len()
    }

    /// Return `true` if no voxels are stored.
    pub fn is_empty(&self) -> bool {
        self.cells.is_empty()
    }

    /// Convert a world-space point to an integer voxel coordinate.
    pub fn world_to_coord(&self, p: [f64; 3]) -> VoxelCoord {
        VoxelCoord::new(
            ((p[0] - self.origin[0]) / self.cell_size).floor() as i32,
            ((p[1] - self.origin[1]) / self.cell_size).floor() as i32,
            ((p[2] - self.origin[2]) / self.cell_size).floor() as i32,
        )
    }

    /// Convert an integer voxel coordinate to its world-space centre.
    pub fn coord_to_world(&self, coord: VoxelCoord) -> [f64; 3] {
        [
            self.origin[0] + (coord.ix as f64 + 0.5) * self.cell_size,
            self.origin[1] + (coord.iy as f64 + 0.5) * self.cell_size,
            self.origin[2] + (coord.iz as f64 + 0.5) * self.cell_size,
        ]
    }

    /// Return the minimum SDF value stored.
    pub fn min_value(&self) -> f64 {
        self.cells.values().copied().fold(f64::INFINITY, f64::min)
    }

    /// Return all voxel coordinates whose SDF is below `threshold`.
    pub fn inside_voxels(&self, threshold: f64) -> Vec<VoxelCoord> {
        self.cells
            .iter()
            .filter(|&(_, &v)| v < threshold)
            .map(|(&c, _)| c)
            .collect()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// VoxelRayQuery — ray vs voxel grid
// ─────────────────────────────────────────────────────────────────────────────

/// Result of a ray-voxel query.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct VoxelRayHit {
    /// Hit voxel coordinates.
    pub coord: [usize; 3],
    /// Ray parameter `t` at the hit point.
    pub t: f64,
    /// World-space hit position.
    pub position: [f64; 3],
    /// SDF value at the hit voxel centre.
    pub sdf: f64,
}

/// Ray vs voxel grid intersection queries.
#[allow(dead_code)]
pub struct VoxelRayQuery;

#[allow(dead_code)]
impl VoxelRayQuery {
    /// Cast a ray `origin + t * direction` against `grid` using DDA traversal.
    ///
    /// Returns the first voxel hit whose SDF is below `surface_threshold`.
    pub fn cast(
        grid: &SdfVoxelGrid,
        ray_origin: [f64; 3],
        ray_dir: [f64; 3],
        max_t: f64,
        surface_threshold: f64,
    ) -> Option<VoxelRayHit> {
        // Normalise direction
        let len =
            (ray_dir[0] * ray_dir[0] + ray_dir[1] * ray_dir[1] + ray_dir[2] * ray_dir[2]).sqrt();
        if len < 1e-14 {
            return None;
        }
        let d = [ray_dir[0] / len, ray_dir[1] / len, ray_dir[2] / len];

        // Find entry point
        let t_entry = Self::ray_aabb_enter(ray_origin, d, grid);
        if t_entry > max_t {
            return None;
        }
        let t_start = t_entry.max(0.0);

        // Sphere-trace along the ray using SDF
        let mut t = t_start;
        let step_min = grid.cell_size * 0.5;
        while t < max_t {
            let p = [
                ray_origin[0] + t * d[0],
                ray_origin[1] + t * d[1],
                ray_origin[2] + t * d[2],
            ];
            let sdf = grid.interpolate(p);
            if sdf.is_infinite() {
                break;
            }
            if sdf < surface_threshold {
                // Refine: find the voxel
                if let Some(coord) = grid.world_to_voxel(p) {
                    let hit_sdf = grid.get(coord[0], coord[1], coord[2]).unwrap_or(sdf);
                    return Some(VoxelRayHit {
                        coord,
                        t,
                        position: p,
                        sdf: hit_sdf,
                    });
                }
            }
            // Step forward by at least step_min
            let step = sdf.abs().max(step_min);
            t += step;
        }
        None
    }

    /// Compute the ray parameter at which the ray enters the grid AABB.
    fn ray_aabb_enter(origin: [f64; 3], dir: [f64; 3], grid: &SdfVoxelGrid) -> f64 {
        let (aabb_min, aabb_max) = grid.aabb();
        let mut t_min = 0.0_f64;
        for axis in 0..3 {
            if dir[axis].abs() < 1e-14 {
                if origin[axis] < aabb_min[axis] || origin[axis] > aabb_max[axis] {
                    return f64::INFINITY;
                }
            } else {
                let t0 = (aabb_min[axis] - origin[axis]) / dir[axis];
                let t1 = (aabb_max[axis] - origin[axis]) / dir[axis];
                let (lo, hi) = if t0 < t1 { (t0, t1) } else { (t1, t0) };
                t_min = t_min.max(lo);
                let t_max = t_min.min(hi);
                if t_max < t_min {
                    return f64::INFINITY;
                }
            }
        }
        t_min
    }

    /// Check whether a ray intersects any voxel with SDF below `threshold`.
    ///
    /// Faster than `cast` (no position computation).
    pub fn intersects(
        grid: &SdfVoxelGrid,
        ray_origin: [f64; 3],
        ray_dir: [f64; 3],
        max_t: f64,
        threshold: f64,
    ) -> bool {
        Self::cast(grid, ray_origin, ray_dir, max_t, threshold).is_some()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// VoxelSphereQuery — sphere vs voxel grid
// ─────────────────────────────────────────────────────────────────────────────

/// Collision result between a sphere and a voxel SDF grid.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct VoxelSphereContact {
    /// World-space contact point on the surface.
    pub point: [f64; 3],
    /// Contact normal pointing from grid surface toward sphere centre.
    pub normal: [f64; 3],
    /// Penetration depth (positive = overlapping).
    pub depth: f64,
}

/// Sphere vs voxel grid SDF collision queries.
#[allow(dead_code)]
pub struct VoxelSphereQuery;

#[allow(dead_code)]
impl VoxelSphereQuery {
    /// Test whether a sphere at `center` with `radius` penetrates the voxel grid.
    ///
    /// Returns a contact if the SDF value at the sphere centre is less than `radius`.
    pub fn query(grid: &SdfVoxelGrid, center: [f64; 3], radius: f64) -> Option<VoxelSphereContact> {
        let sdf = grid.interpolate(center);
        if sdf.is_infinite() || sdf >= radius {
            return None;
        }
        let depth = radius - sdf;

        // Normal from gradient at the nearest surface point
        if let Some(coord) = grid.world_to_voxel(center) {
            let normal = grid.gradient(coord[0], coord[1], coord[2]);
            let contact_point = [
                center[0] - sdf * normal[0],
                center[1] - sdf * normal[1],
                center[2] - sdf * normal[2],
            ];
            Some(VoxelSphereContact {
                point: contact_point,
                normal,
                depth,
            })
        } else {
            // Sphere centre is outside grid — approximate
            Some(VoxelSphereContact {
                point: center,
                normal: [0.0, 1.0, 0.0],
                depth,
            })
        }
    }

    /// Collect all voxels within distance `radius` from `center`.
    pub fn voxels_in_sphere(grid: &SdfVoxelGrid, center: [f64; 3], radius: f64) -> Vec<[usize; 3]> {
        let cs = grid.cell_size;
        let r_cells = (radius / cs).ceil() as usize + 1;
        let base = grid.world_to_voxel(center);
        let (bx, by, bz) = match base {
            Some(c) => (c[0] as i64, c[1] as i64, c[2] as i64),
            None => return Vec::new(),
        };

        let mut result = Vec::new();
        for dz in -(r_cells as i64)..=(r_cells as i64) {
            for dy in -(r_cells as i64)..=(r_cells as i64) {
                for dx in -(r_cells as i64)..=(r_cells as i64) {
                    let ix = bx + dx;
                    let iy = by + dy;
                    let iz = bz + dz;
                    if ix < 0
                        || iy < 0
                        || iz < 0
                        || ix as usize >= grid.nx
                        || iy as usize >= grid.ny
                        || iz as usize >= grid.nz
                    {
                        continue;
                    }
                    let c = grid.voxel_center(ix as usize, iy as usize, iz as usize);
                    let d_sq = (c[0] - center[0]).powi(2)
                        + (c[1] - center[1]).powi(2)
                        + (c[2] - center[2]).powi(2);
                    if d_sq <= radius * radius {
                        result.push([ix as usize, iy as usize, iz as usize]);
                    }
                }
            }
        }
        result
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// VoxelConvexQuery — convex shape vs voxel grid
// ─────────────────────────────────────────────────────────────────────────────

/// An axis-aligned box (convex shape) for voxel collision queries.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct ConvexAabb {
    /// World-space minimum corner.
    pub min: [f64; 3],
    /// World-space maximum corner.
    pub max: [f64; 3],
}

#[allow(dead_code)]
impl ConvexAabb {
    /// Create an AABB from min/max corners.
    pub fn new(min: [f64; 3], max: [f64; 3]) -> Self {
        Self { min, max }
    }

    /// Create an AABB centred at `center` with half-extents `half`.
    pub fn from_center(center: [f64; 3], half: [f64; 3]) -> Self {
        Self {
            min: [
                center[0] - half[0],
                center[1] - half[1],
                center[2] - half[2],
            ],
            max: [
                center[0] + half[0],
                center[1] + half[1],
                center[2] + half[2],
            ],
        }
    }

    /// Return the centre of the AABB.
    pub fn center(&self) -> [f64; 3] {
        [
            (self.min[0] + self.max[0]) * 0.5,
            (self.min[1] + self.max[1]) * 0.5,
            (self.min[2] + self.max[2]) * 0.5,
        ]
    }

    /// Test overlap with another AABB.
    pub fn overlaps(&self, other: &ConvexAabb) -> bool {
        self.min[0] <= other.max[0]
            && self.max[0] >= other.min[0]
            && self.min[1] <= other.max[1]
            && self.max[1] >= other.min[1]
            && self.min[2] <= other.max[2]
            && self.max[2] >= other.min[2]
    }
}

/// Convex AABB vs voxel grid collision queries.
#[allow(dead_code)]
pub struct VoxelConvexQuery;

#[allow(dead_code)]
impl VoxelConvexQuery {
    /// Collect all voxels whose centres lie within `aabb`.
    pub fn voxels_in_aabb(grid: &SdfVoxelGrid, aabb: &ConvexAabb) -> Vec<[usize; 3]> {
        let cs = grid.cell_size;
        let ix0 = (((aabb.min[0] - grid.origin[0]) / cs).floor() as i64).max(0) as usize;
        let iy0 = (((aabb.min[1] - grid.origin[1]) / cs).floor() as i64).max(0) as usize;
        let iz0 = (((aabb.min[2] - grid.origin[2]) / cs).floor() as i64).max(0) as usize;
        let ix1 =
            (((aabb.max[0] - grid.origin[0]) / cs).ceil() as usize).min(grid.nx.saturating_sub(1));
        let iy1 =
            (((aabb.max[1] - grid.origin[1]) / cs).ceil() as usize).min(grid.ny.saturating_sub(1));
        let iz1 =
            (((aabb.max[2] - grid.origin[2]) / cs).ceil() as usize).min(grid.nz.saturating_sub(1));

        let mut result = Vec::new();
        for iz in iz0..=iz1 {
            for iy in iy0..=iy1 {
                for ix in ix0..=ix1 {
                    result.push([ix, iy, iz]);
                }
            }
        }
        result
    }

    /// Test whether any voxel inside `aabb` has SDF below `threshold`.
    pub fn any_inside(grid: &SdfVoxelGrid, aabb: &ConvexAabb, threshold: f64) -> bool {
        for [ix, iy, iz] in Self::voxels_in_aabb(grid, aabb) {
            if let Some(sdf) = grid.get(ix, iy, iz)
                && sdf < threshold
            {
                return true;
            }
        }
        false
    }

    /// Compute the minimum SDF value in the volume covered by `aabb`.
    pub fn min_sdf(grid: &SdfVoxelGrid, aabb: &ConvexAabb) -> f64 {
        Self::voxels_in_aabb(grid, aabb)
            .iter()
            .filter_map(|&[ix, iy, iz]| grid.get(ix, iy, iz))
            .fold(f64::INFINITY, f64::min)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// VoxelCollisionResponse
// ─────────────────────────────────────────────────────────────────────────────

/// Result of a voxel collision response calculation.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct VoxelCollisionImpulse {
    /// Impulse magnitude applied to the dynamic body.
    pub magnitude: f64,
    /// Unit contact normal (pointing away from grid surface).
    pub normal: [f64; 3],
    /// World-space contact point.
    pub contact_point: [f64; 3],
}

/// Voxel-based collision response for a sphere against an SDF grid.
///
/// Computes a penalty impulse proportional to the penetration depth.
#[allow(dead_code)]
pub struct VoxelCollisionResponse;

#[allow(dead_code)]
impl VoxelCollisionResponse {
    /// Compute a penalty response impulse for a sphere penetrating the grid.
    ///
    /// - `contact` — detected contact from `VoxelSphereQuery::query`
    /// - `velocity` — current velocity of the sphere centre `[vx, vy, vz]`
    /// - `mass` — mass of the sphere
    /// - `stiffness` — penalty spring stiffness (force per unit depth)
    /// - `restitution` — coefficient of restitution (0 = fully inelastic)
    ///
    /// Returns the impulse to be applied to the sphere.
    #[allow(clippy::too_many_arguments)]
    pub fn penalty_impulse(
        contact: &VoxelSphereContact,
        velocity: [f64; 3],
        mass: f64,
        stiffness: f64,
        restitution: f64,
        _dt: f64,
    ) -> VoxelCollisionImpulse {
        let n = contact.normal;
        let vn = velocity[0] * n[0] + velocity[1] * n[1] + velocity[2] * n[2];
        // Penalty spring force magnitude
        let f_penalty = stiffness * contact.depth;
        // Restitution contribution: impulse to reverse normal velocity
        let j_restitution = if vn < 0.0 {
            -(1.0 + restitution) * vn * mass
        } else {
            0.0
        };
        let j_total = f_penalty * _dt + j_restitution;
        VoxelCollisionImpulse {
            magnitude: j_total,
            normal: n,
            contact_point: contact.point,
        }
    }

    /// Apply a voxel collision impulse to update sphere velocity.
    ///
    /// Returns the updated velocity `[vx, vy, vz]`.
    pub fn apply_impulse(
        velocity: [f64; 3],
        impulse: &VoxelCollisionImpulse,
        mass: f64,
    ) -> [f64; 3] {
        let n = impulse.normal;
        let dv = impulse.magnitude / mass.max(1e-14);
        [
            velocity[0] + dv * n[0],
            velocity[1] + dv * n[1],
            velocity[2] + dv * n[2],
        ]
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// DynamicVoxelUpdate
// ─────────────────────────────────────────────────────────────────────────────

/// Dynamic voxel update utilities for moving objects.
///
/// Stamps or un-stamps a sphere's occupancy into an SDF grid as it moves.
#[allow(dead_code)]
pub struct DynamicVoxelUpdate;

#[allow(dead_code)]
impl DynamicVoxelUpdate {
    /// Stamp a sphere with `radius` at `center` into `grid`.
    ///
    /// Each voxel centre's distance to `center` minus `radius` is written as the
    /// SDF if it improves (decreases) the stored value.
    pub fn stamp_sphere(grid: &mut SdfVoxelGrid, center: [f64; 3], radius: f64) {
        let voxels = VoxelSphereQuery::voxels_in_sphere(grid, center, radius + grid.cell_size);
        for [ix, iy, iz] in voxels {
            let c = grid.voxel_center(ix, iy, iz);
            let dist = ((c[0] - center[0]).powi(2)
                + (c[1] - center[1]).powi(2)
                + (c[2] - center[2]).powi(2))
            .sqrt();
            let sdf = dist - radius;
            if let Some(idx) = grid.to_index(ix, iy, iz)
                && sdf < grid.data[idx]
            {
                grid.data[idx] = sdf;
            }
        }
    }

    /// Erase a sphere from `grid` (restore to +INF in the sphere region).
    pub fn erase_sphere(grid: &mut SdfVoxelGrid, center: [f64; 3], radius: f64) {
        let voxels = VoxelSphereQuery::voxels_in_sphere(grid, center, radius + grid.cell_size);
        for [ix, iy, iz] in voxels {
            if let Some(idx) = grid.to_index(ix, iy, iz) {
                grid.data[idx] = f64::INFINITY;
            }
        }
    }

    /// Move a sphere from `old_center` to `new_center` in `grid`.
    pub fn move_sphere(
        grid: &mut SdfVoxelGrid,
        old_center: [f64; 3],
        new_center: [f64; 3],
        radius: f64,
    ) {
        Self::erase_sphere(grid, old_center, radius);
        Self::stamp_sphere(grid, new_center, radius);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Internal math helpers
// ─────────────────────────────────────────────────────────────────────────────

#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
fn dist_sq3(a: [f64; 3], b: [f64; 3]) -> f64 {
    let d = sub3(a, b);
    dot3(d, d)
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ---- VoxelCoord ----

    #[test]
    fn test_voxel_coord_new() {
        let c = VoxelCoord::new(1, 2, 3);
        assert_eq!(c.ix, 1);
        assert_eq!(c.iy, 2);
        assert_eq!(c.iz, 3);
    }

    #[test]
    fn test_voxel_coord_manhattan() {
        let a = VoxelCoord::new(0, 0, 0);
        let b = VoxelCoord::new(3, 4, 0);
        assert_eq!(a.manhattan(&b), 7);
    }

    #[test]
    fn test_voxel_coord_chebyshev() {
        let a = VoxelCoord::new(0, 0, 0);
        let b = VoxelCoord::new(3, 1, 2);
        assert_eq!(a.chebyshev(&b), 3);
    }

    #[test]
    fn test_voxel_coord_face_neighbours_count() {
        let c = VoxelCoord::new(5, 5, 5);
        let nb = c.face_neighbours();
        assert_eq!(nb.len(), 6);
    }

    #[test]
    fn test_voxel_coord_face_neighbours_symmetric() {
        let c = VoxelCoord::new(2, 2, 2);
        let nb = c.face_neighbours();
        for n in &nb {
            assert_eq!(c.manhattan(n), 1, "each face neighbour differs by 1");
        }
    }

    // ---- SdfVoxelGrid ----

    #[test]
    fn test_grid_new_dimensions() {
        let g = SdfVoxelGrid::new([0.0; 3], 0.1, 10, 10, 10);
        assert_eq!(g.nx, 10);
        assert_eq!(g.ny, 10);
        assert_eq!(g.nz, 10);
        assert_eq!(g.len(), 1000);
    }

    #[test]
    fn test_grid_is_empty_false() {
        let g = SdfVoxelGrid::new([0.0; 3], 0.1, 2, 2, 2);
        assert!(!g.is_empty());
    }

    #[test]
    fn test_grid_is_empty_true() {
        let g = SdfVoxelGrid::new([0.0; 3], 0.1, 0, 2, 2);
        assert!(g.is_empty());
    }

    #[test]
    fn test_grid_set_get() {
        let mut g = SdfVoxelGrid::new([0.0; 3], 1.0, 4, 4, 4);
        g.set(1, 2, 3, -0.5);
        assert!((g.get(1, 2, 3).unwrap() + 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_grid_get_out_of_bounds_returns_none() {
        let g = SdfVoxelGrid::new([0.0; 3], 1.0, 4, 4, 4);
        assert!(g.get(10, 0, 0).is_none());
    }

    #[test]
    fn test_grid_to_index_round_trip() {
        let g = SdfVoxelGrid::new([0.0; 3], 1.0, 5, 5, 5);
        let idx = g.to_index(2, 3, 4).unwrap();
        // ix + nx*(iy + ny*iz)
        assert_eq!(idx, 2 + 5 * (3 + 5 * 4));
    }

    #[test]
    fn test_grid_world_to_voxel_basic() {
        let g = SdfVoxelGrid::new([0.0; 3], 1.0, 4, 4, 4);
        let coord = g.world_to_voxel([1.5, 0.5, 2.5]).unwrap();
        assert_eq!(coord, [1, 0, 2]);
    }

    #[test]
    fn test_grid_world_to_voxel_outside_returns_none() {
        let g = SdfVoxelGrid::new([0.0; 3], 1.0, 4, 4, 4);
        assert!(g.world_to_voxel([10.0, 0.5, 0.5]).is_none());
    }

    #[test]
    fn test_grid_voxel_center() {
        let g = SdfVoxelGrid::new([0.0; 3], 1.0, 4, 4, 4);
        let c = g.voxel_center(0, 0, 0);
        assert!((c[0] - 0.5).abs() < 1e-12);
        assert!((c[1] - 0.5).abs() < 1e-12);
        assert!((c[2] - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_grid_fill_and_count_inside() {
        let mut g = SdfVoxelGrid::new([0.0; 3], 1.0, 3, 3, 3);
        g.fill(-1.0);
        assert_eq!(g.count_inside(0.0), 27);
    }

    #[test]
    fn test_grid_aabb() {
        let g = SdfVoxelGrid::new([1.0, 2.0, 3.0], 0.5, 4, 4, 4);
        let (mn, mx) = g.aabb();
        assert!((mn[0] - 1.0).abs() < 1e-12);
        assert!((mx[0] - 3.0).abs() < 1e-12); // 1.0 + 4*0.5
    }

    #[test]
    fn test_grid_interpolate_constant() {
        let mut g = SdfVoxelGrid::new([0.0; 3], 1.0, 4, 4, 4);
        g.fill(2.0);
        let v = g.interpolate([1.5, 1.5, 1.5]);
        assert!(
            (v - 2.0).abs() < 1e-10,
            "constant field interpolates to 2.0, got {v}"
        );
    }

    #[test]
    fn test_grid_gradient_interior() {
        let mut g = SdfVoxelGrid::new([0.0; 3], 1.0, 5, 5, 5);
        // Linear field: sdf(ix,iy,iz) = ix (increases in x)
        for iz in 0..5 {
            for iy in 0..5 {
                for ix in 0..5 {
                    g.set(ix, iy, iz, ix as f64);
                }
            }
        }
        let grad = g.gradient(2, 2, 2);
        assert!(grad[0] > 0.9, "gradient should point in +x, got {:?}", grad);
    }

    // ---- MeshVoxelizer ----

    #[test]
    fn test_point_triangle_dist_sq_at_vertex() {
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [0.0, 1.0, 0.0];
        let d = MeshVoxelizer::point_triangle_dist_sq(a, a, b, c);
        assert!(d < 1e-14, "distance to own vertex should be 0, got {d}");
    }

    #[test]
    fn test_point_triangle_dist_sq_above_centroid() {
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [0.0, 1.0, 0.0];
        // Centroid = (1/3, 1/3, 0), point directly above
        let p = [1.0 / 3.0, 1.0 / 3.0, 1.0];
        let d = MeshVoxelizer::point_triangle_dist_sq(p, a, b, c);
        assert!(
            (d - 1.0).abs() < 1e-10,
            "point above centroid dist² should be 1, got {d}"
        );
    }

    #[test]
    fn test_sign_winding_above_positive() {
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [0.0, 1.0, 0.0];
        let p = [0.0, 0.0, 1.0]; // above (positive normal side)
        let s = MeshVoxelizer::sign_winding(p, a, b, c);
        assert!(s > 0.0, "point above should have positive sign");
    }

    #[test]
    fn test_sign_winding_below_negative() {
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [0.0, 1.0, 0.0];
        let p = [0.0, 0.0, -1.0]; // below
        let s = MeshVoxelizer::sign_winding(p, a, b, c);
        assert!(s < 0.0, "point below should have negative sign");
    }

    // ---- FastSweeping ----

    #[test]
    fn test_fast_sweeping_single_seed() {
        let grid = SdfVoxelGrid::new([0.0; 3], 1.0, 3, 3, 3);
        let seed = VoxelCoord::new(1, 1, 1);
        let phi = FastSweeping::compute(&grid, &[(seed, 0.0)]);
        // Centre voxel should be 0
        let centre_idx = 1 + 3 * (1 + 3);
        assert!(phi[centre_idx] < 1e-10, "seed voxel should have dist 0");
        // Corner voxels should be farther
        let corner_idx = 0; // (0,0,0)
        assert!(phi[corner_idx] > 0.5, "corner should be farther from seed");
    }

    #[test]
    fn test_fast_sweeping_two_seeds_boundary() {
        let grid = SdfVoxelGrid::new([0.0; 3], 1.0, 5, 1, 1);
        let s0 = VoxelCoord::new(0, 0, 0);
        let s4 = VoxelCoord::new(4, 0, 0);
        let phi = FastSweeping::compute(&grid, &[(s0, 0.0), (s4, 0.0)]);
        // Middle voxel at ix=2 should be equidistant from both seeds
        let idx2 = 2;
        let idx0 = 0;
        assert!(
            phi[idx2] >= phi[idx0],
            "middle should be no closer than endpoint"
        );
    }

    // ---- SparseVoxelOctree ----

    #[test]
    fn test_sparse_octree_insert_get() {
        let mut tree = SparseVoxelOctree::new([0.0; 3], 1.0, f64::INFINITY);
        tree.insert(VoxelCoord::new(1, 2, 3), -0.5);
        assert!((tree.get(VoxelCoord::new(1, 2, 3)) + 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_sparse_octree_missing_returns_default() {
        let tree = SparseVoxelOctree::new([0.0; 3], 1.0, 99.0);
        assert!((tree.get(VoxelCoord::new(0, 0, 0)) - 99.0).abs() < 1e-12);
    }

    #[test]
    fn test_sparse_octree_remove() {
        let mut tree = SparseVoxelOctree::new([0.0; 3], 1.0, f64::INFINITY);
        let c = VoxelCoord::new(1, 1, 1);
        tree.insert(c, 1.0);
        tree.remove(&c);
        assert!(tree.is_empty());
    }

    #[test]
    fn test_sparse_octree_len() {
        let mut tree = SparseVoxelOctree::new([0.0; 3], 1.0, f64::INFINITY);
        tree.insert(VoxelCoord::new(0, 0, 0), 0.0);
        tree.insert(VoxelCoord::new(1, 0, 0), 0.0);
        assert_eq!(tree.len(), 2);
    }

    #[test]
    fn test_sparse_octree_world_to_coord() {
        let tree = SparseVoxelOctree::new([0.0; 3], 1.0, f64::INFINITY);
        let c = tree.world_to_coord([1.7, 2.1, 0.9]);
        assert_eq!(c.ix, 1);
        assert_eq!(c.iy, 2);
        assert_eq!(c.iz, 0);
    }

    #[test]
    fn test_sparse_octree_inside_voxels() {
        let mut tree = SparseVoxelOctree::new([0.0; 3], 1.0, f64::INFINITY);
        tree.insert(VoxelCoord::new(0, 0, 0), -1.0);
        tree.insert(VoxelCoord::new(1, 0, 0), 2.0);
        let inside = tree.inside_voxels(0.0);
        assert_eq!(inside.len(), 1);
        assert_eq!(inside[0], VoxelCoord::new(0, 0, 0));
    }

    // ---- ConvexAabb ----

    #[test]
    fn test_aabb_overlap() {
        let a = ConvexAabb::new([0.0; 3], [1.0, 1.0, 1.0]);
        let b = ConvexAabb::new([0.5, 0.5, 0.5], [1.5, 1.5, 1.5]);
        assert!(a.overlaps(&b));
    }

    #[test]
    fn test_aabb_no_overlap() {
        let a = ConvexAabb::new([0.0; 3], [1.0, 1.0, 1.0]);
        let b = ConvexAabb::new([2.0, 0.0, 0.0], [3.0, 1.0, 1.0]);
        assert!(!a.overlaps(&b));
    }

    #[test]
    fn test_aabb_center() {
        let a = ConvexAabb::from_center([1.0, 2.0, 3.0], [0.5, 0.5, 0.5]);
        let c = a.center();
        assert!((c[0] - 1.0).abs() < 1e-12);
        assert!((c[1] - 2.0).abs() < 1e-12);
        assert!((c[2] - 3.0).abs() < 1e-12);
    }

    // ---- VoxelConvexQuery ----

    #[test]
    fn test_convex_query_voxels_in_aabb() {
        let g = SdfVoxelGrid::new([0.0; 3], 1.0, 4, 4, 4);
        let aabb = ConvexAabb::new([0.0; 3], [2.0, 2.0, 2.0]);
        let voxels = VoxelConvexQuery::voxels_in_aabb(&g, &aabb);
        assert!(!voxels.is_empty());
    }

    #[test]
    fn test_convex_query_any_inside_true() {
        let mut g = SdfVoxelGrid::new([0.0; 3], 1.0, 4, 4, 4);
        g.set(1, 1, 1, -0.5);
        let aabb = ConvexAabb::new([0.5, 0.5, 0.5], [2.5, 2.5, 2.5]);
        assert!(VoxelConvexQuery::any_inside(&g, &aabb, 0.0));
    }

    #[test]
    fn test_convex_query_min_sdf() {
        let mut g = SdfVoxelGrid::new([0.0; 3], 1.0, 4, 4, 4);
        g.fill(1.0);
        g.set(2, 2, 2, -3.0);
        let aabb = ConvexAabb::new([0.0; 3], [4.0, 4.0, 4.0]);
        let min = VoxelConvexQuery::min_sdf(&g, &aabb);
        assert!(
            (min + 3.0).abs() < 1e-12,
            "min SDF should be -3.0, got {min}"
        );
    }

    // ---- VoxelSphereQuery ----

    #[test]
    fn test_sphere_query_no_contact() {
        let g = SdfVoxelGrid::new([0.0; 3], 1.0, 4, 4, 4);
        // All values are +INF → no contact
        let result = VoxelSphereQuery::query(&g, [2.0, 2.0, 2.0], 0.3);
        assert!(result.is_none());
    }

    #[test]
    fn test_sphere_query_contact_detected() {
        let mut g = SdfVoxelGrid::new([0.0; 3], 1.0, 4, 4, 4);
        // Set SDF < sphere radius at sphere centre voxel
        g.fill(0.1); // sdf=0.1 everywhere
        let result = VoxelSphereQuery::query(&g, [2.0, 2.0, 2.0], 0.5);
        assert!(
            result.is_some(),
            "sphere with radius > SDF should detect contact"
        );
    }

    #[test]
    fn test_sphere_query_depth_positive() {
        let mut g = SdfVoxelGrid::new([0.0; 3], 1.0, 4, 4, 4);
        g.fill(-0.2); // inside solid
        let result = VoxelSphereQuery::query(&g, [2.0, 2.0, 2.0], 0.5);
        if let Some(c) = result {
            assert!(c.depth > 0.0, "penetration depth should be positive");
        }
    }

    #[test]
    fn test_sphere_voxels_in_sphere_non_empty() {
        let g = SdfVoxelGrid::new([0.0; 3], 1.0, 5, 5, 5);
        let voxels = VoxelSphereQuery::voxels_in_sphere(&g, [2.5, 2.5, 2.5], 1.0);
        assert!(!voxels.is_empty(), "should find voxels within radius 1");
    }

    // ---- VoxelRayQuery ----

    #[test]
    fn test_ray_query_no_hit_empty_grid() {
        let g = SdfVoxelGrid::new([0.0; 3], 1.0, 4, 4, 4);
        // All +INF: ray should not hit anything
        let hit = VoxelRayQuery::cast(&g, [-1.0, 2.0, 2.0], [1.0, 0.0, 0.0], 20.0, 0.0);
        assert!(hit.is_none());
    }

    #[test]
    fn test_ray_query_hit_solid_grid() {
        let mut g = SdfVoxelGrid::new([0.0; 3], 1.0, 4, 4, 4);
        g.fill(-1.0); // everything solid
        // Start from inside grid (well within) so interpolation is valid
        let hit = VoxelRayQuery::cast(&g, [1.5, 1.5, 1.5], [1.0, 0.0, 0.0], 20.0, 0.0);
        assert!(hit.is_some(), "ray into solid grid should return a hit");
    }

    #[test]
    fn test_ray_intersects_solid() {
        let mut g = SdfVoxelGrid::new([0.0; 3], 1.0, 4, 4, 4);
        g.fill(-1.0);
        assert!(VoxelRayQuery::intersects(
            &g,
            [1.5, 2.0, 2.0],
            [1.0, 0.0, 0.0],
            10.0,
            0.0
        ));
    }

    // ---- VoxelCollisionResponse ----

    #[test]
    fn test_collision_response_zero_depth_no_impulse() {
        let contact = VoxelSphereContact {
            point: [0.0; 3],
            normal: [0.0, 1.0, 0.0],
            depth: 0.0,
        };
        let v = [0.0, -1.0, 0.0];
        let impulse = VoxelCollisionResponse::penalty_impulse(&contact, v, 1.0, 1000.0, 0.5, 0.01);
        // depth=0 → only restitution term (vn < 0)
        assert!(
            impulse.magnitude >= 0.0,
            "impulse magnitude should be non-negative"
        );
    }

    #[test]
    fn test_collision_response_apply_impulse() {
        let impulse = VoxelCollisionImpulse {
            magnitude: 10.0,
            normal: [0.0, 1.0, 0.0],
            contact_point: [0.0; 3],
        };
        let v = [0.0, -5.0, 0.0];
        let v_new = VoxelCollisionResponse::apply_impulse(v, &impulse, 1.0);
        // dv = 10/1.0 in +y direction: vy = -5 + 10 = 5
        assert!(
            (v_new[1] - 5.0).abs() < 1e-12,
            "vy should be -5+10=5, got {}",
            v_new[1]
        );
    }

    // ---- DynamicVoxelUpdate ----

    #[test]
    fn test_stamp_sphere_sets_sdf() {
        let mut g = SdfVoxelGrid::new([0.0; 3], 1.0, 6, 6, 6);
        g.fill(f64::INFINITY);
        DynamicVoxelUpdate::stamp_sphere(&mut g, [3.0, 3.0, 3.0], 1.0);
        let inside = g.count_inside(0.0);
        assert!(inside > 0, "stamping sphere should create inside voxels");
    }

    #[test]
    fn test_erase_sphere_restores_infinity() {
        let mut g = SdfVoxelGrid::new([0.0; 3], 1.0, 6, 6, 6);
        g.fill(f64::INFINITY);
        DynamicVoxelUpdate::stamp_sphere(&mut g, [3.0, 3.0, 3.0], 1.0);
        DynamicVoxelUpdate::erase_sphere(&mut g, [3.0, 3.0, 3.0], 1.0);
        let inside = g.count_inside(0.0);
        assert_eq!(inside, 0, "after erase, no inside voxels should remain");
    }

    #[test]
    fn test_move_sphere() {
        let mut g = SdfVoxelGrid::new([0.0; 3], 1.0, 8, 8, 8);
        g.fill(f64::INFINITY);
        DynamicVoxelUpdate::move_sphere(&mut g, [2.0, 2.0, 2.0], [6.0, 6.0, 6.0], 0.8);
        // Old location should be erased, new location stamped
        // Check a voxel near old centre is back to +INF
        let old_sdf = g.interpolate([2.0, 2.0, 2.0]);
        let new_sdf = g.interpolate([6.0, 6.0, 6.0]);
        assert!(
            old_sdf > new_sdf,
            "new centre should have lower SDF than old centre"
        );
    }

    // ---- Internal math helpers ----

    #[test]
    fn test_sub3() {
        let r = sub3([3.0, 2.0, 1.0], [1.0, 1.0, 1.0]);
        assert_eq!(r, [2.0, 1.0, 0.0]);
    }

    #[test]
    fn test_add3() {
        let r = add3([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]);
        assert_eq!(r, [5.0, 7.0, 9.0]);
    }

    #[test]
    fn test_cross3_unit_vectors() {
        let r = cross3([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert!((r[2] - 1.0).abs() < 1e-12, "i×j should be k, got {:?}", r);
    }

    #[test]
    fn test_dist_sq3_known() {
        let d = dist_sq3([0.0; 3], [3.0, 4.0, 0.0]);
        assert!((d - 25.0).abs() < 1e-12);
    }

    #[test]
    fn test_scale3() {
        let r = scale3([1.0, 2.0, 3.0], 2.0);
        assert_eq!(r, [2.0, 4.0, 6.0]);
    }
}
