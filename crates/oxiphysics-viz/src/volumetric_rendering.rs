// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Volumetric rendering algorithms for scientific visualization.
//!
//! Provides voxel data structures, trilinear sampling, ray marching,
//! transfer functions, front-to-back alpha compositing, a lightweight
//! marching-cubes iso-surface extractor (vertex positions only), volume
//! histograms, and gradient-magnitude computation.

#![allow(dead_code)]

// ─────────────────────────────────────────────────────────────────────────────
// VolumeData
// ─────────────────────────────────────────────────────────────────────────────

/// A 3-D scalar volume stored as a flat row-major array.
///
/// Indexing: `data[iz * ny * nx + iy * nx + ix]`.
#[derive(Debug, Clone)]
pub struct VolumeData {
    /// Scalar values at each voxel.
    pub data: Vec<f64>,
    /// Number of voxels along the X axis.
    pub nx: usize,
    /// Number of voxels along the Y axis.
    pub ny: usize,
    /// Number of voxels along the Z axis.
    pub nz: usize,
    /// Physical size of each voxel (assumed isotropic).
    pub voxel_size: f64,
}

impl VolumeData {
    /// Create a new volume filled with zeros.
    pub fn new(nx: usize, ny: usize, nz: usize, voxel_size: f64) -> Self {
        Self {
            data: vec![0.0; nx * ny * nz],
            nx,
            ny,
            nz,
            voxel_size,
        }
    }

    /// Total number of voxels.
    pub fn len(&self) -> usize {
        self.data.len()
    }

    /// Return `true` if the volume has no voxels.
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    /// Flat index into `data` for voxel coordinates `(ix, iy, iz)`.
    ///
    /// Returns `None` if the coordinates are out-of-bounds.
    pub fn index(&self, ix: usize, iy: usize, iz: usize) -> Option<usize> {
        if ix < self.nx && iy < self.ny && iz < self.nz {
            Some(iz * self.ny * self.nx + iy * self.nx + ix)
        } else {
            None
        }
    }

    /// Get the scalar value at voxel `(ix, iy, iz)`, clamped to boundary.
    pub fn get(&self, ix: i64, iy: i64, iz: i64) -> f64 {
        let ix = ix.clamp(0, self.nx as i64 - 1) as usize;
        let iy = iy.clamp(0, self.ny as i64 - 1) as usize;
        let iz = iz.clamp(0, self.nz as i64 - 1) as usize;
        self.data[iz * self.ny * self.nx + iy * self.nx + ix]
    }

    /// Physical extent of the volume in world units.
    pub fn extent(&self) -> [f64; 3] {
        [
            self.nx as f64 * self.voxel_size,
            self.ny as f64 * self.voxel_size,
            self.nz as f64 * self.voxel_size,
        ]
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Trilinear sampling
// ─────────────────────────────────────────────────────────────────────────────

/// Sample the volume at a continuous world-space position `(x, y, z)` using
/// trilinear interpolation.
///
/// Positions are in voxel units (0 .. nx-1, etc.). Out-of-bounds access is
/// handled by clamping.
pub fn sample_trilinear(vol: &VolumeData, x: f64, y: f64, z: f64) -> f64 {
    let x0 = x.floor() as i64;
    let y0 = y.floor() as i64;
    let z0 = z.floor() as i64;
    let tx = x - x0 as f64;
    let ty = y - y0 as f64;
    let tz = z - z0 as f64;

    let c000 = vol.get(x0, y0, z0);
    let c100 = vol.get(x0 + 1, y0, z0);
    let c010 = vol.get(x0, y0 + 1, z0);
    let c110 = vol.get(x0 + 1, y0 + 1, z0);
    let c001 = vol.get(x0, y0, z0 + 1);
    let c101 = vol.get(x0 + 1, y0, z0 + 1);
    let c011 = vol.get(x0, y0 + 1, z0 + 1);
    let c111 = vol.get(x0 + 1, y0 + 1, z0 + 1);

    let c00 = c000 * (1.0 - tx) + c100 * tx;
    let c10 = c010 * (1.0 - tx) + c110 * tx;
    let c01 = c001 * (1.0 - tx) + c101 * tx;
    let c11 = c011 * (1.0 - tx) + c111 * tx;

    let c0 = c00 * (1.0 - ty) + c10 * ty;
    let c1 = c01 * (1.0 - ty) + c11 * ty;

    c0 * (1.0 - tz) + c1 * tz
}

// ─────────────────────────────────────────────────────────────────────────────
// Ray marching
// ─────────────────────────────────────────────────────────────────────────────

/// Parameters controlling ray marching through a volume.
#[derive(Debug, Clone)]
pub struct RayMarchParams {
    /// Step size in world units along the ray.
    pub step_size: f64,
    /// Maximum number of steps before terminating.
    pub max_steps: usize,
    /// Absorption coefficient (extinction per unit length).
    pub absorption: f64,
    /// Emission coefficient (brightness per unit density per unit length).
    pub emission: f64,
}

impl Default for RayMarchParams {
    fn default() -> Self {
        Self {
            step_size: 0.5,
            max_steps: 256,
            absorption: 0.1,
            emission: 1.0,
        }
    }
}

/// A single sample along a ray inside a volume.
#[derive(Debug, Clone)]
pub struct RaySample {
    /// World-space position of this sample.
    pub position: [f64; 3],
    /// Sampled density at this position.
    pub density: f64,
    /// RGB color (from transfer function or default grey).
    pub color: [f64; 3],
}

/// March a ray through a volume and collect samples.
///
/// `ray_origin` and `ray_dir` are in voxel-space units. The ray direction need
/// not be normalised (it scales the step). Returns all samples collected
/// before the ray exits the volume bounding box or `max_steps` is reached.
pub fn ray_march(
    vol: &VolumeData,
    ray_origin: [f64; 3],
    ray_dir: [f64; 3],
    params: &RayMarchParams,
) -> Vec<RaySample> {
    let mut samples = Vec::new();
    let [mut px, mut py, mut pz] = ray_origin;
    let [dx, dy, dz] = ray_dir;
    let len = (dx * dx + dy * dy + dz * dz).sqrt();
    if len < 1e-12 {
        return samples;
    }
    let step = params.step_size / len;
    let [sx, sy, sz] = [dx * step, dy * step, dz * step];

    let nx = vol.nx as f64;
    let ny = vol.ny as f64;
    let nz = vol.nz as f64;

    for _ in 0..params.max_steps {
        if px < 0.0 || px >= nx || py < 0.0 || py >= ny || pz < 0.0 || pz >= nz {
            break;
        }
        let density = sample_trilinear(vol, px, py, pz);
        samples.push(RaySample {
            position: [px, py, pz],
            density,
            color: [density, density, density],
        });
        px += sx;
        py += sy;
        pz += sz;
    }
    samples
}

// ─────────────────────────────────────────────────────────────────────────────
// Transfer function
// ─────────────────────────────────────────────────────────────────────────────

/// Map a density scalar to an RGBA color using a pre-defined colormap table.
///
/// `colormap` is a slice of RGBA entries `[r, g, b, a]` (each in `[0, 1]`).
/// The density is mapped to `[0, 1]` and linearly interpolated in the table.
/// Returns `[0, 0, 0, 0]` if the colormap is empty.
pub fn transfer_function(density: f64, colormap: &[[f64; 4]]) -> [f64; 4] {
    if colormap.is_empty() {
        return [0.0; 4];
    }
    let n = colormap.len();
    if n == 1 {
        return colormap[0];
    }
    let t = density.clamp(0.0, 1.0) * (n - 1) as f64;
    let lo = t.floor() as usize;
    let hi = (lo + 1).min(n - 1);
    let frac = t - lo as f64;
    let a = colormap[lo];
    let b = colormap[hi];
    [
        a[0] + (b[0] - a[0]) * frac,
        a[1] + (b[1] - a[1]) * frac,
        a[2] + (b[2] - a[2]) * frac,
        a[3] + (b[3] - a[3]) * frac,
    ]
}

// ─────────────────────────────────────────────────────────────────────────────
// Front-to-back compositing
// ─────────────────────────────────────────────────────────────────────────────

/// Accumulate ray samples into a final RGBA value using front-to-back
/// compositing (Porter–Duff "over" operator).
///
/// Each sample contributes `emission * density * step_size` to brightness and
/// `absorption * density * step_size` to opacity accumulation.
pub fn accumulate_rgba(samples: &[RaySample], params: &RayMarchParams) -> [f64; 4] {
    let mut acc_r = 0.0_f64;
    let mut acc_g = 0.0_f64;
    let mut acc_b = 0.0_f64;
    let mut acc_a = 0.0_f64;

    for s in samples {
        let alpha = (params.absorption * s.density * params.step_size).clamp(0.0, 1.0);
        let transmittance = 1.0 - acc_a;
        let contribution = params.emission * s.density * params.step_size;
        acc_r += transmittance * s.color[0] * contribution;
        acc_g += transmittance * s.color[1] * contribution;
        acc_b += transmittance * s.color[2] * contribution;
        acc_a += transmittance * alpha;
        if acc_a > 0.9999 {
            break;
        }
    }
    [acc_r, acc_g, acc_b, acc_a.min(1.0)]
}

// ─────────────────────────────────────────────────────────────────────────────
// Iso-surface extraction (lightweight marching cubes — vertex positions only)
// ─────────────────────────────────────────────────────────────────────────────

/// Extract iso-surface vertex positions for a given iso-value.
///
/// This is a lightweight marching-cubes variant that returns a flat list of
/// triangle vertex positions (3 vertices per triangle, 3 coordinates each).
/// Normal computation and mesh connectivity are omitted for simplicity.
pub fn iso_surface_value(vol: &VolumeData, iso_val: f64) -> Vec<[f64; 3]> {
    let mut vertices = Vec::new();
    if vol.nx < 2 || vol.ny < 2 || vol.nz < 2 {
        return vertices;
    }

    // Edge table: for each of the 12 edges of a cube, the two corner indices
    #[rustfmt::skip]
    let edge_corners: [[usize; 2]; 12] = [
        [0,1],[1,2],[3,2],[0,3],
        [4,5],[5,6],[7,6],[4,7],
        [0,4],[1,5],[2,6],[3,7],
    ];

    // Corner offsets within a cell: (dx, dy, dz)
    let corner_offsets: [[i64; 3]; 8] = [
        [0, 0, 0],
        [1, 0, 0],
        [1, 1, 0],
        [0, 1, 0],
        [0, 0, 1],
        [1, 0, 1],
        [1, 1, 1],
        [0, 1, 1],
    ];

    let interp = |p0: [f64; 3], v0: f64, p1: [f64; 3], v1: f64| -> [f64; 3] {
        let t = if (v1 - v0).abs() < 1e-12 {
            0.5
        } else {
            (iso_val - v0) / (v1 - v0)
        };
        [
            p0[0] + t * (p1[0] - p0[0]),
            p0[1] + t * (p1[1] - p0[1]),
            p0[2] + t * (p1[2] - p0[2]),
        ]
    };

    for iz in 0..(vol.nz - 1) as i64 {
        for iy in 0..(vol.ny - 1) as i64 {
            for ix in 0..(vol.nx - 1) as i64 {
                // Gather 8 corner values
                let mut corner_vals = [0.0_f64; 8];
                let mut corner_pos = [[0.0_f64; 3]; 8];
                for (ci, &[ddx, ddy, ddz]) in corner_offsets.iter().enumerate() {
                    let cx = (ix + ddx) as f64;
                    let cy = (iy + ddy) as f64;
                    let cz = (iz + ddz) as f64;
                    corner_pos[ci] = [
                        cx * vol.voxel_size,
                        cy * vol.voxel_size,
                        cz * vol.voxel_size,
                    ];
                    corner_vals[ci] = vol.get(ix + ddx, iy + ddy, iz + ddz);
                }
                // Cube index
                let mut cube_idx = 0u8;
                for (ci, &v) in corner_vals.iter().enumerate() {
                    if v < iso_val {
                        cube_idx |= 1 << ci;
                    }
                }
                if cube_idx == 0 || cube_idx == 255 {
                    continue;
                }
                // Compute edge intersection points for the 12 edges
                let mut edge_pts = [[0.0_f64; 3]; 12];
                for (ei, &[c0, c1]) in edge_corners.iter().enumerate() {
                    edge_pts[ei] = interp(
                        corner_pos[c0],
                        corner_vals[c0],
                        corner_pos[c1],
                        corner_vals[c1],
                    );
                }
                // Triangle table (simplified — emit triangles for crossed edges)
                // Use a basic approach: collect crossed edges and fan-triangulate
                let mut crossed: Vec<usize> = (0..12)
                    .filter(|&e| {
                        let [c0, c1] = edge_corners[e];
                        (corner_vals[c0] < iso_val) != (corner_vals[c1] < iso_val)
                    })
                    .collect();
                // Fan triangulate from first crossed edge
                if crossed.len() >= 3 {
                    let p0 = edge_pts[crossed[0]];
                    // Ensure consistent winding (rough)
                    crossed.sort_by(|&a, &b| {
                        let angle_a = edge_pts[a][0].atan2(edge_pts[a][1]);
                        let angle_b = edge_pts[b][0].atan2(edge_pts[b][1]);
                        angle_a
                            .partial_cmp(&angle_b)
                            .unwrap_or(std::cmp::Ordering::Equal)
                    });
                    for tri in crossed.windows(3) {
                        vertices.push(p0);
                        vertices.push(edge_pts[tri[1]]);
                        vertices.push(edge_pts[tri[2]]);
                    }
                }
            }
        }
    }
    vertices
}

// ─────────────────────────────────────────────────────────────────────────────
// Volume histogram
// ─────────────────────────────────────────────────────────────────────────────

/// Compute a histogram of voxel values in a volume.
///
/// Bins are uniformly distributed between `min` and `max` of the volume.
/// Returns a vector of `n_bins` counts. Returns an empty vector for empty
/// volumes or `n_bins == 0`.
pub fn volume_histogram(vol: &VolumeData, n_bins: usize) -> Vec<usize> {
    if vol.data.is_empty() || n_bins == 0 {
        return vec![];
    }
    let min_val = vol.data.iter().cloned().fold(f64::INFINITY, f64::min);
    let max_val = vol.data.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let mut counts = vec![0usize; n_bins];
    if (max_val - min_val).abs() < 1e-15 {
        counts[0] = vol.data.len();
        return counts;
    }
    let range = max_val - min_val;
    for &v in &vol.data {
        let bin = ((v - min_val) / range * n_bins as f64) as usize;
        let bin = bin.min(n_bins - 1);
        counts[bin] += 1;
    }
    counts
}

// ─────────────────────────────────────────────────────────────────────────────
// Gradient magnitude
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the gradient magnitude field of a volume using central differences.
///
/// Returns a new `VolumeData` of the same dimensions where each voxel value is
/// `|∇f| = √((∂f/∂x)² + (∂f/∂y)² + (∂f/∂z)²)`.
/// Boundary voxels use one-sided (forward/backward) differences.
pub fn gradient_magnitude(vol: &VolumeData) -> VolumeData {
    let mut out = VolumeData::new(vol.nx, vol.ny, vol.nz, vol.voxel_size);
    let h = vol.voxel_size;
    for iz in 0..vol.nz as i64 {
        for iy in 0..vol.ny as i64 {
            for ix in 0..vol.nx as i64 {
                let gx = (vol.get(ix + 1, iy, iz) - vol.get(ix - 1, iy, iz)) / (2.0 * h);
                let gy = (vol.get(ix, iy + 1, iz) - vol.get(ix, iy - 1, iz)) / (2.0 * h);
                let gz = (vol.get(ix, iy, iz + 1) - vol.get(ix, iy, iz - 1)) / (2.0 * h);
                let mag = (gx * gx + gy * gy + gz * gz).sqrt();
                if let Some(idx) = out.index(ix as usize, iy as usize, iz as usize) {
                    out.data[idx] = mag;
                }
            }
        }
    }
    out
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn make_uniform_vol(nx: usize, ny: usize, nz: usize, val: f64) -> VolumeData {
        let mut v = VolumeData::new(nx, ny, nz, 1.0);
        for d in v.data.iter_mut() {
            *d = val;
        }
        v
    }

    fn make_gradient_vol(nx: usize, ny: usize, nz: usize) -> VolumeData {
        let mut v = VolumeData::new(nx, ny, nz, 1.0);
        for iz in 0..nz {
            for iy in 0..ny {
                for ix in 0..nx {
                    let idx = iz * ny * nx + iy * nx + ix;
                    v.data[idx] = ix as f64;
                }
            }
        }
        v
    }

    // ── VolumeData basics ─────────────────────────────────────────────────────

    #[test]
    fn test_volume_data_new() {
        let v = VolumeData::new(4, 4, 4, 0.5);
        assert_eq!(v.len(), 64);
        assert!(!v.is_empty());
    }

    #[test]
    fn test_volume_data_index_in_bounds() {
        let v = VolumeData::new(3, 3, 3, 1.0);
        assert!(v.index(0, 0, 0).is_some());
        assert!(v.index(2, 2, 2).is_some());
    }

    #[test]
    fn test_volume_data_index_out_of_bounds() {
        let v = VolumeData::new(3, 3, 3, 1.0);
        assert!(v.index(3, 0, 0).is_none());
    }

    #[test]
    fn test_volume_data_get_clamp() {
        let v = make_uniform_vol(2, 2, 2, 1.0);
        // Out-of-bounds should clamp
        assert!((v.get(-5, 0, 0) - 1.0).abs() < 1e-12);
        assert!((v.get(100, 0, 0) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_volume_data_extent() {
        let v = VolumeData::new(10, 20, 30, 0.1);
        let ext = v.extent();
        assert!((ext[0] - 1.0).abs() < 1e-10);
        assert!((ext[1] - 2.0).abs() < 1e-10);
        assert!((ext[2] - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_volume_data_empty() {
        let v = VolumeData::new(0, 0, 0, 1.0);
        assert!(v.is_empty());
    }

    // ── Trilinear sampling ────────────────────────────────────────────────────

    #[test]
    fn test_sample_trilinear_uniform() {
        let v = make_uniform_vol(4, 4, 4, 3.5);
        let s = sample_trilinear(&v, 1.5, 1.5, 1.5);
        assert!(
            (s - 3.5).abs() < 1e-10,
            "uniform vol should return 3.5, got {s}"
        );
    }

    #[test]
    fn test_sample_trilinear_corner() {
        let v = make_uniform_vol(4, 4, 4, 1.0);
        let s = sample_trilinear(&v, 0.0, 0.0, 0.0);
        assert!((s - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_sample_trilinear_interpolation() {
        // Linear gradient in x: f(x,y,z) = x
        let v = make_gradient_vol(4, 4, 4);
        let s = sample_trilinear(&v, 1.5, 1.0, 1.0);
        assert!((s - 1.5).abs() < 1e-10, "gradient interp: {s}");
    }

    #[test]
    fn test_sample_trilinear_out_of_range_clamps() {
        let v = make_gradient_vol(4, 4, 4);
        // Should clamp, not panic
        let s = sample_trilinear(&v, 10.0, 0.0, 0.0);
        assert!(s.is_finite());
    }

    // ── Ray march ─────────────────────────────────────────────────────────────

    #[test]
    fn test_ray_march_empty_dir() {
        let v = make_uniform_vol(8, 8, 8, 1.0);
        let p = RayMarchParams::default();
        let samples = ray_march(&v, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], &p);
        assert!(samples.is_empty(), "zero direction should yield no samples");
    }

    #[test]
    fn test_ray_march_through_volume() {
        let v = make_uniform_vol(8, 8, 8, 1.0);
        let p = RayMarchParams {
            step_size: 1.0,
            max_steps: 10,
            ..Default::default()
        };
        let samples = ray_march(&v, [0.5, 4.0, 4.0], [1.0, 0.0, 0.0], &p);
        assert!(
            !samples.is_empty(),
            "should have samples through uniform volume"
        );
    }

    #[test]
    fn test_ray_march_density_uniform() {
        let v = make_uniform_vol(8, 8, 8, 2.0);
        let p = RayMarchParams {
            step_size: 1.0,
            max_steps: 5,
            ..Default::default()
        };
        let samples = ray_march(&v, [0.5, 4.0, 4.0], [1.0, 0.0, 0.0], &p);
        for s in &samples {
            assert!(
                (s.density - 2.0).abs() < 1e-10,
                "density should be 2.0, got {}",
                s.density
            );
        }
    }

    #[test]
    fn test_ray_march_exits_volume() {
        let v = make_uniform_vol(4, 4, 4, 1.0);
        let p = RayMarchParams {
            step_size: 1.0,
            max_steps: 100,
            ..Default::default()
        };
        let samples = ray_march(&v, [0.0, 2.0, 2.0], [1.0, 0.0, 0.0], &p);
        // Should stop before max_steps since volume only has 4 voxels along x
        assert!(
            samples.len() <= 6,
            "expected few samples: {}",
            samples.len()
        );
    }

    // ── Transfer function ─────────────────────────────────────────────────────

    #[test]
    fn test_transfer_function_empty_colormap() {
        let rgba = transfer_function(0.5, &[]);
        assert_eq!(rgba, [0.0; 4]);
    }

    #[test]
    fn test_transfer_function_single_entry() {
        let c = [[1.0, 0.0, 0.0, 1.0]];
        let rgba = transfer_function(0.5, &c);
        assert_eq!(rgba, c[0]);
    }

    #[test]
    fn test_transfer_function_two_entry_midpoint() {
        let c = [[0.0, 0.0, 0.0, 0.0], [1.0, 1.0, 1.0, 1.0]];
        let rgba = transfer_function(0.5, &c);
        for &v in &rgba {
            assert!((v - 0.5).abs() < 1e-10, "midpoint should be 0.5, got {v}");
        }
    }

    #[test]
    fn test_transfer_function_clamp_below_zero() {
        let c = [[0.0, 0.0, 0.0, 0.0], [1.0, 1.0, 1.0, 1.0]];
        let rgba = transfer_function(-1.0, &c);
        assert_eq!(rgba, [0.0; 4]);
    }

    #[test]
    fn test_transfer_function_clamp_above_one() {
        let c = [[0.0, 0.0, 0.0, 0.0], [1.0, 1.0, 1.0, 1.0]];
        let rgba = transfer_function(2.0, &c);
        assert_eq!(rgba, [1.0; 4]);
    }

    // ── Accumulate RGBA ───────────────────────────────────────────────────────

    #[test]
    fn test_accumulate_rgba_empty() {
        let p = RayMarchParams::default();
        let rgba = accumulate_rgba(&[], &p);
        assert_eq!(rgba[3], 0.0, "empty samples should have zero alpha");
    }

    #[test]
    fn test_accumulate_rgba_single_sample() {
        let p = RayMarchParams {
            step_size: 1.0,
            absorption: 1.0,
            emission: 1.0,
            ..Default::default()
        };
        let sample = RaySample {
            position: [0.0; 3],
            density: 1.0,
            color: [1.0, 0.0, 0.0],
        };
        let rgba = accumulate_rgba(&[sample], &p);
        assert!(rgba[0] > 0.0, "should have red contribution");
        assert!(rgba[3] > 0.0, "should have non-zero alpha");
    }

    #[test]
    fn test_accumulate_rgba_alpha_max_one() {
        let p = RayMarchParams {
            step_size: 1.0,
            absorption: 10.0,
            emission: 1.0,
            max_steps: 100,
        };
        let samples: Vec<RaySample> = (0..100)
            .map(|i| RaySample {
                position: [i as f64, 0.0, 0.0],
                density: 1.0,
                color: [1.0, 1.0, 1.0],
            })
            .collect();
        let rgba = accumulate_rgba(&samples, &p);
        assert!(rgba[3] <= 1.0, "alpha should never exceed 1.0");
    }

    // ── Volume histogram ──────────────────────────────────────────────────────

    #[test]
    fn test_volume_histogram_uniform() {
        let v = make_uniform_vol(4, 4, 4, 0.5);
        let hist = volume_histogram(&v, 8);
        // All values identical → all in one bin
        let total: usize = hist.iter().sum();
        assert_eq!(total, 64);
        // Only one bin should be non-zero
        let nonzero = hist.iter().filter(|&&c| c > 0).count();
        assert_eq!(nonzero, 1);
    }

    #[test]
    fn test_volume_histogram_n_bins() {
        let v = make_gradient_vol(4, 4, 4);
        let hist = volume_histogram(&v, 10);
        assert_eq!(hist.len(), 10);
    }

    #[test]
    fn test_volume_histogram_empty() {
        let v = VolumeData::new(0, 0, 0, 1.0);
        let hist = volume_histogram(&v, 8);
        assert!(hist.is_empty());
    }

    #[test]
    fn test_volume_histogram_total_count() {
        let v = make_gradient_vol(3, 3, 3);
        let n = 5;
        let hist = volume_histogram(&v, n);
        let total: usize = hist.iter().sum();
        assert_eq!(total, 27);
    }

    // ── Gradient magnitude ────────────────────────────────────────────────────

    #[test]
    fn test_gradient_magnitude_uniform_is_zero() {
        let v = make_uniform_vol(4, 4, 4, 5.0);
        let gm = gradient_magnitude(&v);
        for &val in &gm.data {
            assert!(
                val.abs() < 1e-10,
                "uniform volume has zero gradient, got {val}"
            );
        }
    }

    #[test]
    fn test_gradient_magnitude_linear_x() {
        let v = make_gradient_vol(5, 5, 5);
        let gm = gradient_magnitude(&v);
        // Interior voxels: ∂f/∂x = 1, ∂f/∂y = 0, ∂f/∂z = 0 → |∇f| = 1
        let idx = gm.index(2, 2, 2).unwrap();
        assert!(
            (gm.data[idx] - 1.0).abs() < 1e-10,
            "interior gradient = {}",
            gm.data[idx]
        );
    }

    #[test]
    fn test_gradient_magnitude_same_dimensions() {
        let v = make_gradient_vol(4, 5, 6);
        let gm = gradient_magnitude(&v);
        assert_eq!(gm.nx, 4);
        assert_eq!(gm.ny, 5);
        assert_eq!(gm.nz, 6);
    }

    #[test]
    fn test_gradient_magnitude_non_negative() {
        let v = make_gradient_vol(4, 4, 4);
        let gm = gradient_magnitude(&v);
        for &val in &gm.data {
            assert!(val >= 0.0, "gradient magnitude must be non-negative");
        }
    }

    // ── Iso-surface ───────────────────────────────────────────────────────────

    #[test]
    fn test_iso_surface_small_volume() {
        let v = VolumeData::new(1, 1, 1, 1.0);
        let verts = iso_surface_value(&v, 0.5);
        // Volume too small for marching cubes
        assert!(verts.is_empty());
    }

    #[test]
    fn test_iso_surface_all_above_iso() {
        let v = make_uniform_vol(4, 4, 4, 2.0);
        let verts = iso_surface_value(&v, 1.0);
        // All values above iso → no surface
        assert!(verts.is_empty());
    }

    #[test]
    fn test_iso_surface_all_below_iso() {
        let v = make_uniform_vol(4, 4, 4, 0.5);
        let verts = iso_surface_value(&v, 1.0);
        assert!(verts.is_empty());
    }

    #[test]
    fn test_iso_surface_crossing_produces_vertices() {
        // Half the volume is above, half below iso-value
        let mut v = VolumeData::new(4, 4, 4, 1.0);
        for iz in 0..4_usize {
            for iy in 0..4_usize {
                for ix in 0..4_usize {
                    let idx = iz * 16 + iy * 4 + ix;
                    v.data[idx] = if ix < 2 { 0.0 } else { 2.0 };
                }
            }
        }
        let verts = iso_surface_value(&v, 1.0);
        assert!(
            !verts.is_empty(),
            "crossing volume should produce iso-surface vertices"
        );
    }
}
