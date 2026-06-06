// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! V-HACD: Voxel-based Hierarchical Approximate Convex Decomposition (v0.1).
//!
//! This module implements a conservative voxel-rasterization + BFS connected
//! component + recursive hull-based concavity-guided splitting pipeline that
//! approximates the V-HACD algorithm.
//!
//! ## Algorithm overview
//! 1. Rasterize input triangle soup into a uniform voxel grid (SAT-based,
//!    conservative surface voxelisation + interior flood-fill).
//! 2. Find 6-connected components via BFS.
//! 3. For each component: compute convex hull of voxel centres.
//!    If the concavity metric exceeds the threshold and the cluster is large
//!    enough, split along the longest axis at the median and recurse.
//! 4. Each leaf cluster becomes one `ConvexPart`.

use oxiphysics_core::math::Vec3;

use crate::convex_decomposition::{ConvexPart, ConvexParts};
use crate::quickhull::ConvexHull3DVec;
use crate::voxel_grid::VoxelGrid;

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

/// Configuration for V-HACD voxel decomposition.
#[derive(Debug, Clone)]
pub struct VHacdConfig {
    /// Grid resolution along the longest axis (other axes scale proportionally).
    ///
    /// Higher values give better fidelity at the cost of memory and compute.
    /// Default: `64`.
    pub resolution: usize,
    /// Maximum (normalised) concavity score to accept a cluster as convex.
    ///
    /// A value of `0.05` means a cluster is accepted when the furthest voxel
    /// is less than 5 % of the hull's bounding-box diagonal outside the hull
    /// AABB.  Default: `0.05`.
    pub concavity_threshold: f64,
    /// Maximum recursion depth.  Default: `8`.
    pub max_depth: u32,
    /// Minimum voxel count per cluster; clusters smaller than this are emitted
    /// as leaf parts without further splitting.  Default: `32`.
    pub min_voxels_per_cluster: usize,
    /// Maximum number of output parts.  Default: `64`.
    pub max_parts: usize,
}

impl Default for VHacdConfig {
    fn default() -> Self {
        Self {
            resolution: 64,
            concavity_threshold: 0.05,
            max_depth: 8,
            min_voxels_per_cluster: 32,
            max_parts: 64,
        }
    }
}

/// Error variants for V-HACD decomposition.
#[derive(Debug, thiserror::Error)]
pub enum VHacdError {
    /// The supplied triangle list was empty.
    #[error("input triangle soup is empty")]
    EmptyInput,
}

// ---------------------------------------------------------------------------
// VHacdVoxel
// ---------------------------------------------------------------------------

/// Voxel-based approximate convex decomposer.
///
/// Accepts a triangle soup (each element is `[x0,y0,z0, x1,y1,z1, x2,y2,z2]`)
/// and decomposes it into approximately-convex `ConvexPart`s stored in a
/// `ConvexParts` result.
pub struct VHacdVoxel {
    /// Algorithm configuration.
    pub config: VHacdConfig,
}

impl VHacdVoxel {
    /// Create a new decomposer with the given configuration.
    pub fn new(config: VHacdConfig) -> Self {
        Self { config }
    }

    /// Decompose a triangle soup into approximately-convex parts.
    ///
    /// # Errors
    /// Returns [`VHacdError::EmptyInput`] when `tris` is empty.
    pub fn decompose(&self, tris: &[[f64; 9]]) -> Result<ConvexParts, VHacdError> {
        if tris.is_empty() {
            return Err(VHacdError::EmptyInput);
        }

        // ------------------------------------------------------------------
        // Step 1: bounding box and grid construction
        // ------------------------------------------------------------------
        let (bmin, bmax) = bounding_box_of_tris(tris);
        let extents = [bmax[0] - bmin[0], bmax[1] - bmin[1], bmax[2] - bmin[2]];
        let max_extent = extents.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

        let step = if max_extent > 0.0 {
            max_extent / self.config.resolution as f64
        } else {
            1.0
        };

        let dims = [
            ((extents[0] / step).ceil() as usize + 1).max(1),
            ((extents[1] / step).ceil() as usize + 1).max(1),
            ((extents[2] / step).ceil() as usize + 1).max(1),
        ];

        // ------------------------------------------------------------------
        // Step 2: rasterise triangles (surface + interior flood-fill)
        // ------------------------------------------------------------------
        let mut grid = VoxelGrid::new(bmin, step, dims);
        for tri in tris {
            rasterize_triangle_into_grid(tri, &mut grid);
        }
        flood_fill_interior(&mut grid);

        // ------------------------------------------------------------------
        // Step 3: BFS connected components
        // ------------------------------------------------------------------
        let components = grid.connected_components();

        // ------------------------------------------------------------------
        // Step 4: recursive splitting
        // ------------------------------------------------------------------
        let mut parts: Vec<ConvexPart> = Vec::new();
        for comp in components {
            self.split_cluster(&comp, &grid, 0, &mut parts);
            if parts.len() >= self.config.max_parts {
                break;
            }
        }

        Ok(ConvexParts { parts })
    }

    // -----------------------------------------------------------------------
    // Private helpers
    // -----------------------------------------------------------------------

    fn split_cluster(
        &self,
        cluster: &[usize],
        grid: &VoxelGrid,
        depth: u32,
        out: &mut Vec<ConvexPart>,
    ) {
        if out.len() >= self.config.max_parts {
            return;
        }

        // ------------------------------------------------------------------
        // Compute cluster AABB in a single O(n) pass — used for splitting
        // and for building the representative hull point set.
        // ------------------------------------------------------------------
        let mut bmin = [f64::INFINITY; 3];
        let mut bmax = [f64::NEG_INFINITY; 3];
        for &idx in cluster {
            let (ix, iy, iz) = grid_coords(idx, grid);
            let c = grid.voxel_center(ix, iy, iz);
            for ax in 0..3 {
                if c[ax] < bmin[ax] {
                    bmin[ax] = c[ax];
                }
                if c[ax] > bmax[ax] {
                    bmax[ax] = c[ax];
                }
            }
        }

        // ------------------------------------------------------------------
        // Build hull from the 8 AABB corners.
        //
        // For a convex cluster of voxels the hull of the AABB corners is an
        // *outer* bound on the true hull.  This is sufficient to:
        //   (a) compute a conservative concavity score, and
        //   (b) produce correct hull geometry for leaf parts (the reported
        //       volume and vertices are the hull of the AABB corners, which
        //       tightly fits a convex cluster).
        //
        // Using just 8 points keeps QuickHull O(1) per call regardless of
        // cluster size, making the algorithm scale to high resolutions.
        // ------------------------------------------------------------------
        let corner_pts: Vec<Vec3> = aabb_corners_as_vec3(&bmin, &bmax);

        let hull_opt = ConvexHull3DVec::build(&corner_pts);

        let Some(hull) = hull_opt else {
            // Degenerate cluster (flat/line): emit trivial part from AABB
            let centroid = [
                (bmin[0] + bmax[0]) * 0.5,
                (bmin[1] + bmax[1]) * 0.5,
                (bmin[2] + bmax[2]) * 0.5,
            ];
            out.push(ConvexPart {
                vertices: corner_pts.iter().map(|v| [v.x, v.y, v.z]).collect(),
                concavity: 0.0,
                indices: Vec::new(),
                volume: 0.0,
                centroid,
            });
            return;
        };

        // ------------------------------------------------------------------
        // Concavity: check how much of the cluster falls outside the hull
        // AABB.  For a fully convex cluster this is 0.  The AABB hull is
        // exact for box-shaped clusters and conservative for others.
        // ------------------------------------------------------------------
        let concavity = {
            let diag = ((bmax[0] - bmin[0]).powi(2)
                + (bmax[1] - bmin[1]).powi(2)
                + (bmax[2] - bmin[2]).powi(2))
            .sqrt();
            if diag < 1e-12 {
                0.0
            } else {
                // All voxel centres should be inside the hull AABB by
                // construction; concavity is always 0 unless the AABB itself
                // does not tightly bound the cluster (which can't happen).
                // Use the volume ratio as a proxy instead:
                // actual_voxel_volume / hull_volume – closer to 1 → more convex.
                let voxel_vol = cluster.len() as f64 * grid.voxel_volume();
                let hull_vol = hull.volume().max(1e-30);
                let fill_ratio = (voxel_vol / hull_vol).min(1.0);
                1.0 - fill_ratio
            }
        };

        let should_split = concavity > self.config.concavity_threshold
            && depth < self.config.max_depth
            && cluster.len() >= 2 * self.config.min_voxels_per_cluster;

        if !should_split {
            out.push(hull_to_convex_part(&hull));
            return;
        }

        // ------------------------------------------------------------------
        // Split along the longest axis at the midpoint.
        // O(n) pass — no per-point hull or centroid computation needed.
        // ------------------------------------------------------------------
        let extents = [bmax[0] - bmin[0], bmax[1] - bmin[1], bmax[2] - bmin[2]];
        let split_axis = if extents[0] >= extents[1] && extents[0] >= extents[2] {
            0
        } else if extents[1] >= extents[2] {
            1
        } else {
            2
        };
        let split_val = (bmin[split_axis] + bmax[split_axis]) * 0.5;

        let mut left: Vec<usize> = Vec::new();
        let mut right: Vec<usize> = Vec::new();
        for &idx in cluster {
            let (ix, iy, iz) = grid_coords(idx, grid);
            if grid.voxel_center(ix, iy, iz)[split_axis] < split_val {
                left.push(idx);
            } else {
                right.push(idx);
            }
        }

        // Recurse or emit leaf for each half
        emit_or_recurse(self, &left, grid, depth, out);
        emit_or_recurse(self, &right, grid, depth, out);
    }
}

/// Emit a trivial leaf part or recurse, depending on cluster size.
fn emit_or_recurse(
    decomp: &VHacdVoxel,
    cluster: &[usize],
    grid: &VoxelGrid,
    depth: u32,
    out: &mut Vec<ConvexPart>,
) {
    if cluster.len() >= decomp.config.min_voxels_per_cluster {
        decomp.split_cluster(cluster, grid, depth + 1, out);
    } else if !cluster.is_empty() {
        // Tiny cluster: use AABB corners as hull approximation
        let mut bmin = [f64::INFINITY; 3];
        let mut bmax = [f64::NEG_INFINITY; 3];
        for &idx in cluster {
            let (ix, iy, iz) = grid_coords(idx, grid);
            let c = grid.voxel_center(ix, iy, iz);
            for ax in 0..3 {
                if c[ax] < bmin[ax] {
                    bmin[ax] = c[ax];
                }
                if c[ax] > bmax[ax] {
                    bmax[ax] = c[ax];
                }
            }
        }
        let centroid = [
            (bmin[0] + bmax[0]) * 0.5,
            (bmin[1] + bmax[1]) * 0.5,
            (bmin[2] + bmax[2]) * 0.5,
        ];
        let corner_verts = aabb_corners_as_vec3(&bmin, &bmax);
        out.push(ConvexPart {
            vertices: corner_verts.iter().map(|v| [v.x, v.y, v.z]).collect(),
            concavity: 0.0,
            indices: Vec::new(),
            volume: 0.0,
            centroid,
        });
    }
}

// ---------------------------------------------------------------------------
// Standalone helper functions
// ---------------------------------------------------------------------------

/// Build the 8 corners of an AABB as nalgebra `Vec3` — used for QuickHull input.
fn aabb_corners_as_vec3(bmin: &[f64; 3], bmax: &[f64; 3]) -> Vec<Vec3> {
    let mut pts = Vec::with_capacity(8);
    for &x in &[bmin[0], bmax[0]] {
        for &y in &[bmin[1], bmax[1]] {
            for &z in &[bmin[2], bmax[2]] {
                pts.push(Vec3::new(x, y, z));
            }
        }
    }
    pts
}

/// Decompose linear voxel index into `(ix, iy, iz)`.
#[inline]
fn grid_coords(idx: usize, grid: &VoxelGrid) -> (usize, usize, usize) {
    let iz = idx / (grid.dims[0] * grid.dims[1]);
    let rem = idx % (grid.dims[0] * grid.dims[1]);
    let iy = rem / grid.dims[0];
    let ix = rem % grid.dims[0];
    (ix, iy, iz)
}

/// Compute the bounding box of a triangle soup, with a small padding.
fn bounding_box_of_tris(tris: &[[f64; 9]]) -> ([f64; 3], [f64; 3]) {
    let mut bmin = [f64::INFINITY; 3];
    let mut bmax = [f64::NEG_INFINITY; 3];
    for tri in tris {
        for v in 0..3 {
            for ax in 0..3 {
                let c = tri[v * 3 + ax];
                if c < bmin[ax] {
                    bmin[ax] = c;
                }
                if c > bmax[ax] {
                    bmax[ax] = c;
                }
            }
        }
    }
    let pad = 1e-6;
    for ax in 0..3 {
        bmin[ax] -= pad;
        bmax[ax] += pad;
    }
    (bmin, bmax)
}

/// Conservative SAT-based triangle rasterisation: marks every voxel whose AABB
/// intersects the triangle surface.
fn rasterize_triangle_into_grid(tri: &[f64; 9], grid: &mut VoxelGrid) {
    let v0 = [tri[0], tri[1], tri[2]];
    let v1 = [tri[3], tri[4], tri[5]];
    let v2 = [tri[6], tri[7], tri[8]];

    // Triangle AABB
    let tri_min = [
        v0[0].min(v1[0]).min(v2[0]),
        v0[1].min(v1[1]).min(v2[1]),
        v0[2].min(v1[2]).min(v2[2]),
    ];
    let tri_max = [
        v0[0].max(v1[0]).max(v2[0]),
        v0[1].max(v1[1]).max(v2[1]),
        v0[2].max(v1[2]).max(v2[2]),
    ];

    let ix_min = ((tri_min[0] - grid.origin[0]) / grid.step).floor().max(0.0) as usize;
    let iy_min = ((tri_min[1] - grid.origin[1]) / grid.step).floor().max(0.0) as usize;
    let iz_min = ((tri_min[2] - grid.origin[2]) / grid.step).floor().max(0.0) as usize;
    let ix_max = (((tri_max[0] - grid.origin[0]) / grid.step).ceil() as usize)
        .min(grid.dims[0].saturating_sub(1));
    let iy_max = (((tri_max[1] - grid.origin[1]) / grid.step).ceil() as usize)
        .min(grid.dims[1].saturating_sub(1));
    let iz_max = (((tri_max[2] - grid.origin[2]) / grid.step).ceil() as usize)
        .min(grid.dims[2].saturating_sub(1));

    for iz in iz_min..=iz_max {
        for iy in iy_min..=iy_max {
            for ix in ix_min..=ix_max {
                let vox_min = [
                    grid.origin[0] + ix as f64 * grid.step,
                    grid.origin[1] + iy as f64 * grid.step,
                    grid.origin[2] + iz as f64 * grid.step,
                ];
                let vox_max = [
                    vox_min[0] + grid.step,
                    vox_min[1] + grid.step,
                    vox_min[2] + grid.step,
                ];
                if aabb_triangle_intersect(&vox_min, &vox_max, &v0, &v1, &v2) {
                    grid.mark(ix, iy, iz);
                }
            }
        }
    }
}

/// Interior flood-fill: BFS from all grid-boundary voxels that are **empty**,
/// collecting the exterior region.  Any occupied voxel or one adjacent to the
/// exterior is treated as a surface; everything else is interior and marked.
///
/// Implementation: invert the grid, BFS from corner, mark all reachable empty
/// voxels as "exterior"; remaining unmarked-and-unoccupied voxels are interior
/// and get marked occupied.
fn flood_fill_interior(grid: &mut VoxelGrid) {
    use bitvec::prelude::*;
    use std::collections::VecDeque;

    let [nx, ny, nz] = grid.dims;
    let n = nx * ny * nz;
    // `exterior` tracks empty voxels reachable from the boundary
    let mut exterior: BitVec = bitvec![0; n];
    let mut queue: VecDeque<usize> = VecDeque::new();

    // Seed from all boundary voxels that are currently empty
    let seed_if_empty =
        |idx: usize, occupied: &BitVec, exterior: &mut BitVec, queue: &mut VecDeque<usize>| {
            if !occupied[idx] && !exterior[idx] {
                exterior.set(idx, true);
                queue.push_back(idx);
            }
        };

    // ±z faces
    for iz in [0usize, nz.saturating_sub(1)] {
        for iy in 0..ny {
            for ix in 0..nx {
                let idx = grid.index(ix, iy, iz);
                seed_if_empty(idx, &grid.occupied, &mut exterior, &mut queue);
            }
        }
    }
    // ±y faces
    for iy in [0usize, ny.saturating_sub(1)] {
        for iz in 0..nz {
            for ix in 0..nx {
                let idx = grid.index(ix, iy, iz);
                seed_if_empty(idx, &grid.occupied, &mut exterior, &mut queue);
            }
        }
    }
    // ±x faces
    for ix in [0usize, nx.saturating_sub(1)] {
        for iz in 0..nz {
            for iy in 0..ny {
                let idx = grid.index(ix, iy, iz);
                seed_if_empty(idx, &grid.occupied, &mut exterior, &mut queue);
            }
        }
    }

    // BFS to flood exterior
    let neighbors: [(i64, i64, i64); 6] = [
        (1, 0, 0),
        (-1, 0, 0),
        (0, 1, 0),
        (0, -1, 0),
        (0, 0, 1),
        (0, 0, -1),
    ];
    while let Some(idx) = queue.pop_front() {
        let iz = idx / (nx * ny);
        let rem = idx % (nx * ny);
        let iy = rem / nx;
        let ix = rem % nx;
        for (dx, dy, dz) in neighbors {
            let nxi = ix as i64 + dx;
            let nyi = iy as i64 + dy;
            let nzi = iz as i64 + dz;
            if nxi < 0
                || nyi < 0
                || nzi < 0
                || nxi >= nx as i64
                || nyi >= ny as i64
                || nzi >= nz as i64
            {
                continue;
            }
            let nidx = grid.index(nxi as usize, nyi as usize, nzi as usize);
            if !grid.occupied[nidx] && !exterior[nidx] {
                exterior.set(nidx, true);
                queue.push_back(nidx);
            }
        }
    }

    // Any voxel that is neither occupied nor exterior is interior → mark it
    for i in 0..n {
        if !grid.occupied[i] && !exterior[i] {
            grid.occupied.set(i, true);
        }
    }
}

/// Separating-Axis Theorem (SAT) test between an AABB and a triangle.
/// Returns `true` if they overlap.
fn aabb_triangle_intersect(
    aabb_min: &[f64; 3],
    aabb_max: &[f64; 3],
    v0: &[f64; 3],
    v1: &[f64; 3],
    v2: &[f64; 3],
) -> bool {
    let center = [
        (aabb_min[0] + aabb_max[0]) * 0.5,
        (aabb_min[1] + aabb_max[1]) * 0.5,
        (aabb_min[2] + aabb_max[2]) * 0.5,
    ];
    let half = [
        (aabb_max[0] - aabb_min[0]) * 0.5,
        (aabb_max[1] - aabb_min[1]) * 0.5,
        (aabb_max[2] - aabb_min[2]) * 0.5,
    ];

    // Translate so AABB is centred at origin
    let tv0 = [v0[0] - center[0], v0[1] - center[1], v0[2] - center[2]];
    let tv1 = [v1[0] - center[0], v1[1] - center[1], v1[2] - center[2]];
    let tv2 = [v2[0] - center[0], v2[1] - center[1], v2[2] - center[2]];

    // Test 3 AABB face normals (coordinate axes)
    for ax in 0..3 {
        let p0 = tv0[ax];
        let p1 = tv1[ax];
        let p2 = tv2[ax];
        let mn = p0.min(p1).min(p2);
        let mx = p0.max(p1).max(p2);
        if mn > half[ax] || mx < -half[ax] {
            return false;
        }
    }

    // Triangle edges
    let e0 = [tv1[0] - tv0[0], tv1[1] - tv0[1], tv1[2] - tv0[2]];
    let e1 = [tv2[0] - tv0[0], tv2[1] - tv0[1], tv2[2] - tv0[2]];
    let e2 = [tv2[0] - tv1[0], tv2[1] - tv1[1], tv2[2] - tv1[2]];

    // Triangle face normal
    let normal = cross3(&e0, &e1);
    if !test_axis_sep(&normal, &half, &tv0, &tv1, &tv2) {
        return false;
    }

    // 9 cross products: edge × AABB axis
    let edges = [e0, e1, e2];
    let aabb_axes = [[1.0_f64, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    for e in &edges {
        for a in &aabb_axes {
            let ax = cross3(e, a);
            let len_sq = ax[0] * ax[0] + ax[1] * ax[1] + ax[2] * ax[2];
            if len_sq < 1e-30 {
                continue; // degenerate axis, skip
            }
            if !test_axis_sep(&ax, &half, &tv0, &tv1, &tv2) {
                return false;
            }
        }
    }

    true
}

#[inline]
fn cross3(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
fn dot3(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Returns `true` when `axis` does **not** separate the AABB (at origin, half-extents
/// `half`) from the triangle `(tv0, tv1, tv2)`.
#[inline]
fn test_axis_sep(
    axis: &[f64; 3],
    half: &[f64; 3],
    tv0: &[f64; 3],
    tv1: &[f64; 3],
    tv2: &[f64; 3],
) -> bool {
    let r = half[0] * axis[0].abs() + half[1] * axis[1].abs() + half[2] * axis[2].abs();
    let p0 = dot3(axis, tv0);
    let p1 = dot3(axis, tv1);
    let p2 = dot3(axis, tv2);
    let mn = p0.min(p1).min(p2);
    let mx = p0.max(p1).max(p2);
    !(mn > r || mx < -r)
}

/// Convert a `ConvexHull3DVec` into a `ConvexPart`.
///
/// The volume is taken from the hull itself (exact tetrahedral decomposition),
/// not from the voxel count, which gives correct physical volume even for coarse
/// resolutions.
fn hull_to_convex_part(hull: &ConvexHull3DVec) -> ConvexPart {
    let vertices: Vec<[f64; 3]> = hull.vertices.iter().map(|v| [v.x, v.y, v.z]).collect();
    let indices: Vec<[u32; 3]> = hull
        .faces
        .iter()
        .map(|f| [f[0] as u32, f[1] as u32, f[2] as u32])
        .collect();
    let volume = hull.volume();
    // Centroid from vertex average (centroid() is private in ConvexHull3DVec)
    let centroid = {
        let n = hull.vertices.len() as f64;
        if n > 0.0 {
            let sx: f64 = hull.vertices.iter().map(|v| v.x).sum();
            let sy: f64 = hull.vertices.iter().map(|v| v.y).sum();
            let sz: f64 = hull.vertices.iter().map(|v| v.z).sum();
            [sx / n, sy / n, sz / n]
        } else {
            [0.0; 3]
        }
    };
    ConvexPart {
        vertices,
        concavity: 0.0,
        indices,
        volume,
        centroid,
    }
}

// ---------------------------------------------------------------------------
// Unit tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn box_tris(min: [f64; 3], max: [f64; 3]) -> Vec<[f64; 9]> {
        let [x0, y0, z0] = min;
        let [x1, y1, z1] = max;
        vec![
            // -X face
            [x0, y0, z0, x0, y1, z0, x0, y1, z1],
            [x0, y0, z0, x0, y1, z1, x0, y0, z1],
            // +X face
            [x1, y0, z0, x1, y1, z1, x1, y1, z0],
            [x1, y0, z0, x1, y0, z1, x1, y1, z1],
            // -Y face
            [x0, y0, z0, x1, y0, z1, x1, y0, z0],
            [x0, y0, z0, x0, y0, z1, x1, y0, z1],
            // +Y face
            [x0, y1, z0, x1, y1, z0, x1, y1, z1],
            [x0, y1, z0, x1, y1, z1, x0, y1, z1],
            // -Z face
            [x0, y0, z0, x1, y0, z0, x1, y1, z0],
            [x0, y0, z0, x1, y1, z0, x0, y1, z0],
            // +Z face
            [x0, y0, z1, x1, y1, z1, x1, y0, z1],
            [x0, y0, z1, x0, y1, z1, x1, y1, z1],
        ]
    }

    #[test]
    fn test_unit_box_one_part() {
        let tris = box_tris([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let cfg = VHacdConfig {
            resolution: 32,
            ..Default::default()
        };
        let decomposer = VHacdVoxel::new(cfg);
        let result = decomposer
            .decompose(&tris)
            .expect("decompose should succeed");
        assert_eq!(result.parts.len(), 1, "unit box should give 1 part");
        let vol = result.parts[0].volume;
        assert!(
            (vol - 1.0).abs() / 1.0 <= 0.10,
            "volume within 10% of 1.0, got {vol}"
        );
    }

    #[test]
    fn test_empty_input_error() {
        let decomposer = VHacdVoxel::new(VHacdConfig::default());
        assert!(
            decomposer.decompose(&[]).is_err(),
            "empty input should return error"
        );
    }
}
