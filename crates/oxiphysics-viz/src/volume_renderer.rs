// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Volume rendering data structures: voxel grids, transfer functions,
//! ray-casting, isosurface extraction, and per-grid statistics.
//!
//! All positions and scalar values use `f64` for precision.
//! No nalgebra dependency is introduced here; vectors are `[f64; 3]`.

// ---------------------------------------------------------------------------
// VoxelGrid
// ---------------------------------------------------------------------------

/// A regular 3-D scalar grid.
///
/// `dimensions` is `[nx, ny, nz]`.  `spacing` is the cell size in each
/// dimension.  `data` is stored in C order: index = `iz * ny * nx + iy * nx + ix`.
#[derive(Debug, Clone)]
pub struct VoxelGrid {
    /// Number of cells along each axis: `[nx, ny, nz]`.
    pub dimensions: [usize; 3],
    /// Physical spacing between adjacent cell centres along each axis.
    pub spacing: [f64; 3],
    /// Flattened scalar data (length must equal `nx * ny * nz`).
    pub data: Vec<f64>,
}

impl VoxelGrid {
    /// Construct a new `VoxelGrid` filled with zeros.
    ///
    /// `dimensions` — `[nx, ny, nz]`.
    /// `spacing` — physical cell size in each direction.
    pub fn new(dimensions: [usize; 3], spacing: [f64; 3]) -> Self {
        let n = dimensions[0] * dimensions[1] * dimensions[2];
        Self {
            dimensions,
            spacing,
            data: vec![0.0; n],
        }
    }

    /// Construct a grid from existing data.
    ///
    /// Panics if `data.len() != nx * ny * nz`.
    pub fn from_data(dimensions: [usize; 3], spacing: [f64; 3], data: Vec<f64>) -> Self {
        let n = dimensions[0] * dimensions[1] * dimensions[2];
        assert_eq!(data.len(), n, "data length mismatch");
        Self {
            dimensions,
            spacing,
            data,
        }
    }

    /// Return the scalar value at integer voxel coordinates `(ix, iy, iz)`.
    ///
    /// Returns `0.0` if any coordinate is out of bounds.
    pub fn get(&self, ix: usize, iy: usize, iz: usize) -> f64 {
        let [nx, ny, nz] = self.dimensions;
        if ix >= nx || iy >= ny || iz >= nz {
            return 0.0;
        }
        self.data[iz * ny * nx + iy * nx + ix]
    }

    /// Set the scalar value at integer voxel coordinates `(ix, iy, iz)`.
    ///
    /// Does nothing if any coordinate is out of bounds.
    pub fn set(&mut self, ix: usize, iy: usize, iz: usize, value: f64) {
        let [nx, ny, nz] = self.dimensions;
        if ix >= nx || iy >= ny || iz >= nz {
            return;
        }
        self.data[iz * ny * nx + iy * nx + ix] = value;
    }

    /// Sample the grid at an arbitrary physical position using **trilinear interpolation**.
    ///
    /// The grid origin is at `(0, 0, 0)` in physical space.  Returns `0.0`
    /// for positions outside the grid extent.
    pub fn sample(&self, pos: [f64; 3]) -> f64 {
        let [nx, ny, nz] = self.dimensions;
        let [sx, sy, sz] = self.spacing;

        // Convert to fractional grid coordinates
        let fx = pos[0] / sx;
        let fy = pos[1] / sy;
        let fz = pos[2] / sz;

        // Clamp to valid interior range
        if fx < 0.0 || fy < 0.0 || fz < 0.0 {
            return 0.0;
        }
        if fx > (nx - 1) as f64 || fy > (ny - 1) as f64 || fz > (nz - 1) as f64 {
            return 0.0;
        }

        let x0 = fx.floor() as usize;
        let y0 = fy.floor() as usize;
        let z0 = fz.floor() as usize;
        let x1 = (x0 + 1).min(nx - 1);
        let y1 = (y0 + 1).min(ny - 1);
        let z1 = (z0 + 1).min(nz - 1);

        let tx = fx - x0 as f64;
        let ty = fy - y0 as f64;
        let tz = fz - z0 as f64;

        // Trilinear interpolation
        let c000 = self.get(x0, y0, z0);
        let c100 = self.get(x1, y0, z0);
        let c010 = self.get(x0, y1, z0);
        let c110 = self.get(x1, y1, z0);
        let c001 = self.get(x0, y0, z1);
        let c101 = self.get(x1, y0, z1);
        let c011 = self.get(x0, y1, z1);
        let c111 = self.get(x1, y1, z1);

        let c00 = c000 * (1.0 - tx) + c100 * tx;
        let c10 = c010 * (1.0 - tx) + c110 * tx;
        let c01 = c001 * (1.0 - tx) + c101 * tx;
        let c11 = c011 * (1.0 - tx) + c111 * tx;

        let c0 = c00 * (1.0 - ty) + c10 * ty;
        let c1 = c01 * (1.0 - ty) + c11 * ty;

        c0 * (1.0 - tz) + c1 * tz
    }

    /// Return the physical extent `[width, height, depth]` of the grid.
    pub fn extent(&self) -> [f64; 3] {
        [
            (self.dimensions[0] - 1) as f64 * self.spacing[0],
            (self.dimensions[1] - 1) as f64 * self.spacing[1],
            (self.dimensions[2] - 1) as f64 * self.spacing[2],
        ]
    }
}

// ---------------------------------------------------------------------------
// TransferFunction
// ---------------------------------------------------------------------------

/// A piece-wise linear transfer function mapping scalar values to opacity + colour.
///
/// Control points are `(scalar_value, opacity, [r, g, b])`.  Values between
/// control points are linearly interpolated.  Outside the range, the nearest
/// endpoint is used.
#[derive(Debug, Clone)]
pub struct TransferFunction {
    /// Sorted list of `(scalar, opacity, [r,g,b])` control points.
    ///
    /// Must be sorted by scalar value in ascending order.
    pub points: Vec<(f64, f64, [f64; 3])>,
}

impl TransferFunction {
    /// Create a new, empty transfer function.
    pub fn new() -> Self {
        Self { points: Vec::new() }
    }

    /// Add a control point `(scalar, opacity, colour)`.
    ///
    /// The point list is kept sorted by scalar value after insertion.
    pub fn add_point(&mut self, scalar: f64, opacity: f64, color: [f64; 3]) {
        self.points.push((scalar, opacity, color));
        self.points
            .sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    }

    /// Evaluate the transfer function at `value`.
    ///
    /// Returns `(opacity, [r, g, b])`.  Returns `(0.0, [0.0; 3])` if there
    /// are no control points.
    pub fn eval(&self, value: f64) -> (f64, [f64; 3]) {
        if self.points.is_empty() {
            return (0.0, [0.0; 3]);
        }
        if self.points.len() == 1 {
            let (_, op, col) = self.points[0];
            return (op, col);
        }
        // Below first point
        if value <= self.points[0].0 {
            let (_, op, col) = self.points[0];
            return (op, col);
        }
        // Above last point
        let last = self.points.len() - 1;
        if value >= self.points[last].0 {
            let (_, op, col) = self.points[last];
            return (op, col);
        }
        // Find bracketing interval
        let idx = self.points.partition_point(|(s, _, _)| *s <= value) - 1;
        let (s0, op0, col0) = self.points[idx];
        let (s1, op1, col1) = self.points[idx + 1];
        let t = (value - s0) / (s1 - s0);
        let opacity = op0 + t * (op1 - op0);
        let color = [
            col0[0] + t * (col1[0] - col0[0]),
            col0[1] + t * (col1[1] - col0[1]),
            col0[2] + t * (col1[2] - col0[2]),
        ];
        (opacity, color)
    }
}

impl Default for TransferFunction {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// RaycastVolume
// ---------------------------------------------------------------------------

/// Ray-casting volume renderer.
///
/// Integrates opacity and colour along a ray through the [`VoxelGrid`]
/// using front-to-back alpha compositing.
#[derive(Debug, Clone)]
pub struct RaycastVolume {
    /// The voxel grid to render.
    pub voxel_grid: VoxelGrid,
    /// Transfer function mapping scalar values to opacity+colour.
    pub transfer_fn: TransferFunction,
    /// Step size along the ray in physical units.
    pub step_size: f64,
    /// Maximum number of integration steps.
    pub max_steps: usize,
}

impl RaycastVolume {
    /// Construct a new `RaycastVolume`.
    pub fn new(
        voxel_grid: VoxelGrid,
        transfer_fn: TransferFunction,
        step_size: f64,
        max_steps: usize,
    ) -> Self {
        Self {
            voxel_grid,
            transfer_fn,
            step_size,
            max_steps,
        }
    }

    /// Trace a ray from `origin` in direction `dir` and return an RGBA colour.
    ///
    /// `dir` does not need to be normalised; it is normalised internally.
    /// The result `[r, g, b, a]` has components in `[0, 1]`.
    pub fn render_ray(&self, origin: [f64; 3], dir: [f64; 3]) -> [f64; 4] {
        // Normalise direction
        let len = (dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]).sqrt();
        if len < 1e-12 {
            return [0.0, 0.0, 0.0, 0.0];
        }
        let d = [dir[0] / len, dir[1] / len, dir[2] / len];

        let mut pos = origin;
        let mut accumulated_r = 0.0_f64;
        let mut accumulated_g = 0.0_f64;
        let mut accumulated_b = 0.0_f64;
        let mut accumulated_a = 0.0_f64;

        for _ in 0..self.max_steps {
            let scalar = self.voxel_grid.sample(pos);
            let (raw_opacity, color) = self.transfer_fn.eval(scalar);

            // Scale opacity by step size for consistent density
            let opacity = 1.0 - (-raw_opacity * self.step_size).exp();

            // Front-to-back alpha compositing
            let weight = opacity * (1.0 - accumulated_a);
            accumulated_r += color[0] * weight;
            accumulated_g += color[1] * weight;
            accumulated_b += color[2] * weight;
            accumulated_a += weight;

            // Early termination
            if accumulated_a >= 0.99 {
                break;
            }

            pos[0] += d[0] * self.step_size;
            pos[1] += d[1] * self.step_size;
            pos[2] += d[2] * self.step_size;
        }

        [
            accumulated_r.clamp(0.0, 1.0),
            accumulated_g.clamp(0.0, 1.0),
            accumulated_b.clamp(0.0, 1.0),
            accumulated_a.clamp(0.0, 1.0),
        ]
    }
}

// ---------------------------------------------------------------------------
// IsosurfaceExtractor (marching cubes — simplified)
// ---------------------------------------------------------------------------

/// Marching-cubes isosurface extractor.
///
/// Uses a simplified lookup table: for each cube, if the iso-level splits
/// the cube diagonally, the crossing point is computed by linear interpolation
/// along the edge.
pub struct IsosurfaceExtractor;

impl IsosurfaceExtractor {
    /// Extract an isosurface from `grid` at `isovalue`.
    ///
    /// Returns a list of triangle vertices (groups of three consecutive entries
    /// form one triangle).  The position of each vertex is computed by linear
    /// interpolation along cube edges.
    pub fn extract_isosurface(grid: &VoxelGrid, isovalue: f64) -> Vec<[f64; 3]> {
        let [nx, ny, nz] = grid.dimensions;
        let [sx, sy, sz] = grid.spacing;
        let mut vertices: Vec<[f64; 3]> = Vec::new();

        // Iterate over all cubes
        for iz in 0..nz.saturating_sub(1) {
            for iy in 0..ny.saturating_sub(1) {
                for ix in 0..nx.saturating_sub(1) {
                    let v: [f64; 8] = [
                        grid.get(ix, iy, iz),
                        grid.get(ix + 1, iy, iz),
                        grid.get(ix + 1, iy + 1, iz),
                        grid.get(ix, iy + 1, iz),
                        grid.get(ix, iy, iz + 1),
                        grid.get(ix + 1, iy, iz + 1),
                        grid.get(ix + 1, iy + 1, iz + 1),
                        grid.get(ix, iy + 1, iz + 1),
                    ];

                    // Corner positions in physical space
                    let p: [[f64; 3]; 8] = [
                        [(ix) as f64 * sx, (iy) as f64 * sy, (iz) as f64 * sz],
                        [(ix + 1) as f64 * sx, (iy) as f64 * sy, (iz) as f64 * sz],
                        [(ix + 1) as f64 * sx, (iy + 1) as f64 * sy, (iz) as f64 * sz],
                        [(ix) as f64 * sx, (iy + 1) as f64 * sy, (iz) as f64 * sz],
                        [(ix) as f64 * sx, (iy) as f64 * sy, (iz + 1) as f64 * sz],
                        [(ix + 1) as f64 * sx, (iy) as f64 * sy, (iz + 1) as f64 * sz],
                        [
                            (ix + 1) as f64 * sx,
                            (iy + 1) as f64 * sy,
                            (iz + 1) as f64 * sz,
                        ],
                        [(ix) as f64 * sx, (iy + 1) as f64 * sy, (iz + 1) as f64 * sz],
                    ];

                    // Cube index: bitmask of corners above isovalue
                    let cube_idx: u8 = (0..8)
                        .filter(|&i| v[i] >= isovalue)
                        .fold(0u8, |acc, i| acc | (1 << i));

                    if cube_idx == 0 || cube_idx == 0xFF {
                        continue;
                    }

                    // Linear interpolation along an edge between corners a and b
                    let lerp = |a: usize, b: usize| -> [f64; 3] {
                        let dv = v[b] - v[a];
                        let t = if dv.abs() < 1e-12 {
                            0.5
                        } else {
                            (isovalue - v[a]) / dv
                        };
                        [
                            p[a][0] + t * (p[b][0] - p[a][0]),
                            p[a][1] + t * (p[b][1] - p[a][1]),
                            p[a][2] + t * (p[b][2] - p[a][2]),
                        ]
                    };

                    // Simplified triangulation using the 12 possible edge crossings.
                    // For each edge, compute the crossing point if the cube index
                    // indicates crossing (one endpoint inside, one outside).
                    let edge_crossings: [Option<[f64; 3]>; 12] = [
                        if (cube_idx & 0x01 != 0) != (cube_idx & 0x02 != 0) {
                            Some(lerp(0, 1))
                        } else {
                            None
                        },
                        if (cube_idx & 0x02 != 0) != (cube_idx & 0x04 != 0) {
                            Some(lerp(1, 2))
                        } else {
                            None
                        },
                        if (cube_idx & 0x04 != 0) != (cube_idx & 0x08 != 0) {
                            Some(lerp(2, 3))
                        } else {
                            None
                        },
                        if (cube_idx & 0x08 != 0) != (cube_idx & 0x01 != 0) {
                            Some(lerp(3, 0))
                        } else {
                            None
                        },
                        if (cube_idx & 0x10 != 0) != (cube_idx & 0x20 != 0) {
                            Some(lerp(4, 5))
                        } else {
                            None
                        },
                        if (cube_idx & 0x20 != 0) != (cube_idx & 0x40 != 0) {
                            Some(lerp(5, 6))
                        } else {
                            None
                        },
                        if (cube_idx & 0x40 != 0) != (cube_idx & 0x80 != 0) {
                            Some(lerp(6, 7))
                        } else {
                            None
                        },
                        if (cube_idx & 0x80 != 0) != (cube_idx & 0x10 != 0) {
                            Some(lerp(7, 4))
                        } else {
                            None
                        },
                        if (cube_idx & 0x01 != 0) != (cube_idx & 0x10 != 0) {
                            Some(lerp(0, 4))
                        } else {
                            None
                        },
                        if (cube_idx & 0x02 != 0) != (cube_idx & 0x20 != 0) {
                            Some(lerp(1, 5))
                        } else {
                            None
                        },
                        if (cube_idx & 0x04 != 0) != (cube_idx & 0x40 != 0) {
                            Some(lerp(2, 6))
                        } else {
                            None
                        },
                        if (cube_idx & 0x08 != 0) != (cube_idx & 0x80 != 0) {
                            Some(lerp(3, 7))
                        } else {
                            None
                        },
                    ];

                    // Collect crossing points and triangulate with a fan
                    let crossings: Vec<[f64; 3]> =
                        edge_crossings.iter().filter_map(|x| *x).collect();
                    if crossings.len() >= 3 {
                        let pivot = crossings[0];
                        for k in 1..crossings.len() - 1 {
                            vertices.push(pivot);
                            vertices.push(crossings[k]);
                            vertices.push(crossings[k + 1]);
                        }
                    }
                }
            }
        }

        vertices
    }
}

// ---------------------------------------------------------------------------
// VolumeStats
// ---------------------------------------------------------------------------

/// Per-grid scalar statistics including min, max, mean, and a histogram.
#[derive(Debug, Clone)]
pub struct VolumeStats {
    /// Minimum scalar value in the grid.
    pub min: f64,
    /// Maximum scalar value in the grid.
    pub max: f64,
    /// Mean scalar value.
    pub mean: f64,
    /// Histogram bin counts.  `bins[i]` is the number of voxels in the i-th bin.
    pub bins: Vec<u64>,
}

impl VolumeStats {
    /// Compute statistics for a [`VoxelGrid`] with the given number of histogram bins.
    ///
    /// `num_bins` must be at least 1.
    pub fn compute(grid: &VoxelGrid, num_bins: usize) -> Self {
        let num_bins = num_bins.max(1);
        let data = &grid.data;
        if data.is_empty() {
            return Self {
                min: 0.0,
                max: 0.0,
                mean: 0.0,
                bins: vec![0; num_bins],
            };
        }

        let mut min = data[0];
        let mut max = data[0];
        let mut sum = 0.0_f64;
        for &v in data {
            if v < min {
                min = v;
            }
            if v > max {
                max = v;
            }
            sum += v;
        }
        let mean = sum / data.len() as f64;

        let mut bins = vec![0u64; num_bins];
        let range = max - min;
        for &v in data {
            let idx = if range < 1e-30 {
                0
            } else {
                let t = (v - min) / range;
                ((t * num_bins as f64) as usize).min(num_bins - 1)
            };
            bins[idx] += 1;
        }

        Self {
            min,
            max,
            mean,
            bins,
        }
    }

    /// Return the total number of voxels represented by the histogram.
    pub fn total_count(&self) -> u64 {
        self.bins.iter().sum()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- VoxelGrid ---

    #[test]
    fn voxel_grid_new_is_zero() {
        let g = VoxelGrid::new([4, 4, 4], [1.0, 1.0, 1.0]);
        assert!(g.data.iter().all(|&v| v == 0.0));
    }

    #[test]
    fn voxel_grid_set_and_get() {
        let mut g = VoxelGrid::new([4, 4, 4], [1.0, 1.0, 1.0]);
        g.set(1, 2, 3, 42.0);
        assert_eq!(g.get(1, 2, 3), 42.0);
    }

    #[test]
    fn voxel_grid_get_out_of_bounds_returns_zero() {
        let g = VoxelGrid::new([4, 4, 4], [1.0, 1.0, 1.0]);
        assert_eq!(g.get(10, 0, 0), 0.0);
    }

    #[test]
    fn voxel_grid_sample_at_corner() {
        let mut g = VoxelGrid::new([4, 4, 4], [1.0, 1.0, 1.0]);
        g.set(0, 0, 0, 5.0);
        // Exactly at the origin corner
        let v = g.sample([0.0, 0.0, 0.0]);
        assert!((v - 5.0).abs() < 1e-10);
    }

    #[test]
    fn voxel_grid_sample_outside_returns_zero() {
        let g = VoxelGrid::new([4, 4, 4], [1.0, 1.0, 1.0]);
        assert_eq!(g.sample([-1.0, 0.0, 0.0]), 0.0);
        assert_eq!(g.sample([100.0, 0.0, 0.0]), 0.0);
    }

    #[test]
    fn voxel_grid_trilinear_midpoint_uniform() {
        // Grid filled with constant 3.0: any sample should return 3.0
        let g = VoxelGrid::from_data([3, 3, 3], [1.0, 1.0, 1.0], vec![3.0; 27]);
        let v = g.sample([0.5, 0.5, 0.5]);
        assert!((v - 3.0).abs() < 1e-10);
    }

    #[test]
    fn voxel_grid_trilinear_midpoint_gradient() {
        // Linear ramp along x: v(i,j,k) = i
        let dims = [5, 3, 3];
        let mut data = vec![0.0; 5 * 3 * 3];
        for iz in 0..3 {
            for iy in 0..3 {
                for ix in 0..5 {
                    data[iz * 3 * 5 + iy * 5 + ix] = ix as f64;
                }
            }
        }
        let g = VoxelGrid::from_data(dims, [1.0, 1.0, 1.0], data);
        // At x=1.5, y=0, z=0 → interpolated value should be 1.5
        let v = g.sample([1.5, 0.0, 0.0]);
        assert!((v - 1.5).abs() < 1e-10);
    }

    #[test]
    fn voxel_grid_extent() {
        let g = VoxelGrid::new([5, 3, 2], [2.0, 3.0, 4.0]);
        let ext = g.extent();
        assert!((ext[0] - 8.0).abs() < 1e-10);
        assert!((ext[1] - 6.0).abs() < 1e-10);
        assert!((ext[2] - 4.0).abs() < 1e-10);
    }

    #[test]
    fn voxel_grid_from_data_panics_on_mismatch() {
        let result = std::panic::catch_unwind(|| {
            VoxelGrid::from_data([2, 2, 2], [1.0, 1.0, 1.0], vec![0.0; 5]);
        });
        assert!(result.is_err());
    }

    // --- TransferFunction ---

    #[test]
    fn transfer_fn_empty_returns_zero() {
        let tf = TransferFunction::new();
        let (op, col) = tf.eval(0.5);
        assert_eq!(op, 0.0);
        assert_eq!(col, [0.0; 3]);
    }

    #[test]
    fn transfer_fn_single_point() {
        let mut tf = TransferFunction::new();
        tf.add_point(0.5, 0.8, [1.0, 0.0, 0.0]);
        let (op, col) = tf.eval(0.5);
        assert!((op - 0.8).abs() < 1e-12);
        assert!((col[0] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn transfer_fn_below_range_clamps_to_first() {
        let mut tf = TransferFunction::new();
        tf.add_point(0.0, 0.1, [0.0, 0.0, 1.0]);
        tf.add_point(1.0, 0.9, [1.0, 0.0, 0.0]);
        let (op, _) = tf.eval(-1.0);
        assert!((op - 0.1).abs() < 1e-12);
    }

    #[test]
    fn transfer_fn_above_range_clamps_to_last() {
        let mut tf = TransferFunction::new();
        tf.add_point(0.0, 0.1, [0.0, 0.0, 1.0]);
        tf.add_point(1.0, 0.9, [1.0, 0.0, 0.0]);
        let (op, _) = tf.eval(2.0);
        assert!((op - 0.9).abs() < 1e-12);
    }

    #[test]
    fn transfer_fn_interpolates_opacity() {
        let mut tf = TransferFunction::new();
        tf.add_point(0.0, 0.0, [0.0, 0.0, 0.0]);
        tf.add_point(1.0, 1.0, [1.0, 1.0, 1.0]);
        let (op, col) = tf.eval(0.5);
        assert!((op - 0.5).abs() < 1e-12);
        assert!((col[0] - 0.5).abs() < 1e-12);
    }

    #[test]
    fn transfer_fn_interpolates_color() {
        let mut tf = TransferFunction::new();
        tf.add_point(0.0, 0.0, [1.0, 0.0, 0.0]);
        tf.add_point(1.0, 0.0, [0.0, 1.0, 0.0]);
        let (_, col) = tf.eval(0.25);
        assert!((col[0] - 0.75).abs() < 1e-10);
        assert!((col[1] - 0.25).abs() < 1e-10);
    }

    #[test]
    fn transfer_fn_sorted_after_insertion() {
        let mut tf = TransferFunction::new();
        tf.add_point(1.0, 0.9, [1.0, 0.0, 0.0]);
        tf.add_point(0.0, 0.1, [0.0, 0.0, 1.0]);
        // Even though 1.0 was inserted first, the list should be sorted
        assert!(tf.points[0].0 < tf.points[1].0);
    }

    // --- RaycastVolume ---

    #[test]
    fn raycast_zero_opacity_returns_transparent() {
        let grid = VoxelGrid::new([4, 4, 4], [1.0, 1.0, 1.0]);
        let tf = TransferFunction::default(); // no points → opacity 0
        let vol = RaycastVolume::new(grid, tf, 0.1, 50);
        let rgba = vol.render_ray([0.0, 0.0, -1.0], [0.0, 0.0, 1.0]);
        assert!(rgba[3] < 1e-6, "alpha should be ~0 with no opacity");
    }

    #[test]
    fn raycast_opaque_field_accumulates_alpha() {
        // Uniform high-opacity grid
        let grid = VoxelGrid::from_data([4, 4, 4], [1.0, 1.0, 1.0], vec![1.0; 64]);
        let mut tf = TransferFunction::new();
        tf.add_point(0.0, 10.0, [1.0, 0.0, 0.0]);
        tf.add_point(1.0, 10.0, [1.0, 0.0, 0.0]);
        let vol = RaycastVolume::new(grid, tf, 0.1, 200);
        let rgba = vol.render_ray([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]);
        assert!(
            rgba[3] > 0.5,
            "dense field should accumulate significant alpha"
        );
    }

    #[test]
    fn raycast_zero_dir_returns_transparent() {
        let grid = VoxelGrid::new([4, 4, 4], [1.0, 1.0, 1.0]);
        let tf = TransferFunction::default();
        let vol = RaycastVolume::new(grid, tf, 0.1, 50);
        let rgba = vol.render_ray([0.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        assert_eq!(rgba, [0.0, 0.0, 0.0, 0.0]);
    }

    #[test]
    fn raycast_output_clamped_to_0_1() {
        let grid = VoxelGrid::from_data([4, 4, 4], [1.0, 1.0, 1.0], vec![1.0; 64]);
        let mut tf = TransferFunction::new();
        tf.add_point(1.0, 100.0, [2.0, 2.0, 2.0]); // intentionally > 1
        let vol = RaycastVolume::new(grid, tf, 0.1, 100);
        let rgba = vol.render_ray([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        for &v in &rgba {
            assert!(
                (0.0..=1.0).contains(&v),
                "RGBA components must be in [0,1], got {v}"
            );
        }
    }

    #[test]
    fn raycast_step_outside_grid_terminates() {
        // Ray that stays entirely outside the grid — grid samples 0.0 everywhere.
        // Use a TF that maps value 0.0 to opacity 0.0, so alpha stays at 0.
        let grid = VoxelGrid::new([4, 4, 4], [1.0, 1.0, 1.0]);
        let mut tf = TransferFunction::new();
        // Only map high values (1.0) to high opacity; outside-grid (0.0) is fully transparent.
        tf.add_point(0.5, 0.0, [0.0, 0.0, 0.0]);
        tf.add_point(1.0, 1.0, [1.0, 0.0, 0.0]);
        let vol = RaycastVolume::new(grid, tf, 10.0, 10);
        // Start ray far outside the grid bounds [0..3]^3 with spacing 1.0
        let rgba = vol.render_ray([1000.0, 1000.0, 1000.0], [1.0, 0.0, 0.0]);
        assert!(
            rgba[3] < 1e-6,
            "alpha should be 0 when ray is entirely outside grid"
        );
    }

    // --- IsosurfaceExtractor ---

    #[test]
    fn isosurface_empty_grid_no_vertices() {
        let g = VoxelGrid::new([2, 2, 2], [1.0, 1.0, 1.0]);
        let verts = IsosurfaceExtractor::extract_isosurface(&g, 0.5);
        // All zeros below isovalue → no crossing
        assert!(verts.is_empty());
    }

    #[test]
    fn isosurface_uniform_above_isovalue_no_surface() {
        let g = VoxelGrid::from_data([2, 2, 2], [1.0, 1.0, 1.0], vec![1.0; 8]);
        let verts = IsosurfaceExtractor::extract_isosurface(&g, 0.5);
        // All above → no crossing
        assert!(verts.is_empty());
    }

    #[test]
    fn isosurface_produces_vertices_on_crossing() {
        // One corner above, rest below → at least one triangle
        let mut g = VoxelGrid::new([2, 2, 2], [1.0, 1.0, 1.0]);
        g.set(1, 1, 1, 1.0); // one corner above 0.5
        let verts = IsosurfaceExtractor::extract_isosurface(&g, 0.5);
        assert!(!verts.is_empty(), "should generate at least one vertex");
        assert_eq!(verts.len() % 3, 0, "vertex count must be multiple of 3");
    }

    #[test]
    fn isosurface_vertices_lie_on_edges() {
        // Ramp grid: v(i,j,k) = i; isovalue = 0.5 → surface should be at x=0.5
        let dims = [3, 2, 2];
        let mut data = vec![0.0; 3 * 2 * 2];
        for iz in 0..2 {
            for iy in 0..2 {
                for ix in 0..3 {
                    data[iz * 2 * 3 + iy * 3 + ix] = ix as f64;
                }
            }
        }
        let g = VoxelGrid::from_data(dims, [1.0, 1.0, 1.0], data);
        let verts = IsosurfaceExtractor::extract_isosurface(&g, 0.5);
        assert!(!verts.is_empty());
        // All x-coordinates of extracted vertices should be ~0.5
        for v in &verts {
            assert!(
                (v[0] - 0.5).abs() < 1e-10,
                "vertex x should be 0.5, got {}",
                v[0]
            );
        }
    }

    #[test]
    fn isosurface_single_cell_grid_no_panic() {
        // 1×1×1 grid has no edges to iterate over
        let g = VoxelGrid::from_data([1, 1, 1], [1.0, 1.0, 1.0], vec![1.0]);
        let verts = IsosurfaceExtractor::extract_isosurface(&g, 0.5);
        assert!(verts.is_empty());
    }

    // --- VolumeStats ---

    #[test]
    fn volume_stats_constant_grid() {
        let g = VoxelGrid::from_data([3, 3, 3], [1.0, 1.0, 1.0], vec![5.0; 27]);
        let stats = VolumeStats::compute(&g, 10);
        assert!((stats.min - 5.0).abs() < 1e-12);
        assert!((stats.max - 5.0).abs() < 1e-12);
        assert!((stats.mean - 5.0).abs() < 1e-12);
    }

    #[test]
    fn volume_stats_range() {
        let data: Vec<f64> = (0..8).map(|i| i as f64).collect();
        let g = VoxelGrid::from_data([2, 2, 2], [1.0, 1.0, 1.0], data);
        let stats = VolumeStats::compute(&g, 8);
        assert!((stats.min - 0.0).abs() < 1e-12);
        assert!((stats.max - 7.0).abs() < 1e-12);
        assert!((stats.mean - 3.5).abs() < 1e-12);
    }

    #[test]
    fn volume_stats_histogram_total_count() {
        let data: Vec<f64> = (0..27).map(|i| i as f64).collect();
        let g = VoxelGrid::from_data([3, 3, 3], [1.0, 1.0, 1.0], data);
        let stats = VolumeStats::compute(&g, 10);
        assert_eq!(stats.total_count(), 27);
    }

    #[test]
    fn volume_stats_num_bins_respected() {
        let g = VoxelGrid::from_data([2, 2, 2], [1.0, 1.0, 1.0], vec![1.0; 8]);
        let stats = VolumeStats::compute(&g, 16);
        assert_eq!(stats.bins.len(), 16);
    }

    #[test]
    fn volume_stats_empty_grid() {
        let g = VoxelGrid {
            dimensions: [0, 0, 0],
            spacing: [1.0, 1.0, 1.0],
            data: vec![],
        };
        let stats = VolumeStats::compute(&g, 4);
        assert_eq!(stats.total_count(), 0);
        assert_eq!(stats.bins.len(), 4);
    }

    #[test]
    fn volume_stats_bins_at_least_one() {
        let g = VoxelGrid::from_data([2, 2, 2], [1.0, 1.0, 1.0], vec![1.0; 8]);
        let stats = VolumeStats::compute(&g, 0); // min(0, 1) = 1
        assert_eq!(stats.bins.len(), 1);
        assert_eq!(stats.bins[0], 8);
    }
}
