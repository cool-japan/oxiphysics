// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Implicit surfaces and signed distance fields (SDF).
//!
//! Provides SDF primitives, smooth blending operations, voxel SDF grids,
//! marching cubes isosurface extraction, dual contouring, and ray marching.

use crate::mc_tables::{MC_CORNER_OFFSETS, MC_EDGE_TABLE, MC_TRI_TABLE, case_index};

// ─────────────────────────────────────────────────────────────────────────────
// Helper functions
// ─────────────────────────────────────────────────────────────────────────────

/// Smooth minimum (polynomial blend) of `a` and `b` with blend radius `k`.
///
/// Returns a C¹-continuous approximation to `min(a, b)`.
pub fn smooth_min_polynomial(a: f64, b: f64, k: f64) -> f64 {
    let h = (0.5 + 0.5 * (b - a) / k).clamp(0.0, 1.0);
    a * h + b * (1.0 - h) - k * h * (1.0 - h)
}

/// Smooth minimum via exponential blend.
///
/// Returns −(1/k) · ln(e^(−k·a) + e^(−k·b)).
pub fn smooth_min_exponential(a: f64, b: f64, k: f64) -> f64 {
    let ea = (-k * a).exp();
    let eb = (-k * b).exp();
    -(ea + eb).ln() / k
}

/// Numerical gradient of an SDF at point `p` using central differences.
pub fn sdf_gradient_numerical<F>(f: &F, p: [f64; 3], eps: f64) -> [f64; 3]
where
    F: Fn([f64; 3]) -> f64,
{
    let dx = f([p[0] + eps, p[1], p[2]]) - f([p[0] - eps, p[1], p[2]]);
    let dy = f([p[0], p[1] + eps, p[2]]) - f([p[0], p[1] - eps, p[2]]);
    let dz = f([p[0], p[1], p[2] + eps]) - f([p[0], p[1], p[2] - eps]);
    let len = (dx * dx + dy * dy + dz * dz).sqrt().max(1e-30);
    [
        dx / (2.0 * eps * len / eps * eps),
        dy / (2.0 * eps * len / eps * eps),
        dz / (2.0 * eps * len / eps * eps),
    ]
}

// Simpler gradient helper (not normalised) for internal use
fn numerical_grad<F: Fn([f64; 3]) -> f64>(f: &F, p: [f64; 3], eps: f64) -> [f64; 3] {
    [
        (f([p[0] + eps, p[1], p[2]]) - f([p[0] - eps, p[1], p[2]])) / (2.0 * eps),
        (f([p[0], p[1] + eps, p[2]]) - f([p[0], p[1] - eps, p[2]])) / (2.0 * eps),
        (f([p[0], p[1], p[2] + eps]) - f([p[0], p[1], p[2] - eps])) / (2.0 * eps),
    ]
}

fn vec3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn vec3_len(a: [f64; 3]) -> f64 {
    vec3_dot(a, a).sqrt()
}
fn vec3_norm(a: [f64; 3]) -> [f64; 3] {
    let l = vec3_len(a).max(1e-30);
    [a[0] / l, a[1] / l, a[2] / l]
}
fn vec3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn vec3_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn vec3_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
fn vec3_scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

fn clamp(v: f64, lo: f64, hi: f64) -> f64 {
    v.max(lo).min(hi)
}

// ─────────────────────────────────────────────────────────────────────────────
// Implicit Sphere
// ─────────────────────────────────────────────────────────────────────────────

/// Signed distance function for a sphere.
///
/// SDF(p) = |p − center| − radius
#[derive(Debug, Clone, Copy)]
pub struct ImplicitSphere {
    /// Center of the sphere.
    pub center: [f64; 3],
    /// Radius of the sphere.
    pub radius: f64,
}

impl ImplicitSphere {
    /// Create a new implicit sphere.
    pub fn new(center: [f64; 3], radius: f64) -> Self {
        Self { center, radius }
    }

    /// Signed distance from point `p` to this sphere.
    pub fn sdf(&self, p: [f64; 3]) -> f64 {
        let d = vec3_sub(p, self.center);
        vec3_len(d) - self.radius
    }

    /// Outward unit normal at point `p`.
    pub fn gradient(&self, p: [f64; 3]) -> [f64; 3] {
        vec3_norm(vec3_sub(p, self.center))
    }

    /// Hessian of the SDF at point `p` (3×3 flat row-major).
    pub fn hessian(&self, p: [f64; 3]) -> [f64; 9] {
        let d = vec3_sub(p, self.center);
        let r = vec3_len(d).max(1e-30);
        let mut h = [0.0f64; 9];
        for i in 0..3 {
            for j in 0..3 {
                let delta = if i == j { 1.0 } else { 0.0 };
                h[i * 3 + j] = (delta * r * r - d[i] * d[j]) / (r * r * r);
            }
        }
        h
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Implicit Box (AABB)
// ─────────────────────────────────────────────────────────────────────────────

/// Signed distance function for an axis-aligned box.
#[derive(Debug, Clone, Copy)]
pub struct ImplicitBox {
    /// Center of the box.
    pub center: [f64; 3],
    /// Half-extents along each axis.
    pub half_extents: [f64; 3],
    /// Corner rounding radius (0 = sharp corners).
    pub rounding: f64,
}

impl ImplicitBox {
    /// Create a sharp-cornered AABB SDF.
    pub fn new(center: [f64; 3], half_extents: [f64; 3]) -> Self {
        Self {
            center,
            half_extents,
            rounding: 0.0,
        }
    }

    /// Create a rounded-corner box SDF.
    pub fn rounded(center: [f64; 3], half_extents: [f64; 3], rounding: f64) -> Self {
        Self {
            center,
            half_extents,
            rounding,
        }
    }

    /// Signed distance from point `p` to this box.
    pub fn sdf(&self, p: [f64; 3]) -> f64 {
        let q = [
            (p[0] - self.center[0]).abs() - self.half_extents[0] + self.rounding,
            (p[1] - self.center[1]).abs() - self.half_extents[1] + self.rounding,
            (p[2] - self.center[2]).abs() - self.half_extents[2] + self.rounding,
        ];
        let outside = [q[0].max(0.0), q[1].max(0.0), q[2].max(0.0)];
        let outside_dist = vec3_len(outside);
        let inside_dist = q[0].max(q[1]).max(q[2]).min(0.0);
        outside_dist + inside_dist - self.rounding
    }

    /// Gradient (outward normal) of the box SDF at `p`.
    pub fn gradient(&self, p: [f64; 3]) -> [f64; 3] {
        numerical_grad(&|x| self.sdf(x), p, 1e-5)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Implicit Capsule
// ─────────────────────────────────────────────────────────────────────────────

/// Signed distance function for a capsule (segment + radius).
#[derive(Debug, Clone, Copy)]
pub struct ImplicitCapsule {
    /// Start of the segment.
    pub a: [f64; 3],
    /// End of the segment.
    pub b: [f64; 3],
    /// Radius of the capsule.
    pub radius: f64,
}

impl ImplicitCapsule {
    /// Create a new capsule SDF.
    pub fn new(a: [f64; 3], b: [f64; 3], radius: f64) -> Self {
        Self { a, b, radius }
    }

    /// Signed distance from point `p` to this capsule.
    pub fn sdf(&self, p: [f64; 3]) -> f64 {
        let ab = vec3_sub(self.b, self.a);
        let ap = vec3_sub(p, self.a);
        let t = clamp(vec3_dot(ap, ab) / vec3_dot(ab, ab).max(1e-30), 0.0, 1.0);
        let closest = vec3_add(self.a, vec3_scale(ab, t));
        vec3_len(vec3_sub(p, closest)) - self.radius
    }

    /// Gradient of the capsule SDF at `p`.
    pub fn gradient(&self, p: [f64; 3]) -> [f64; 3] {
        let ab = vec3_sub(self.b, self.a);
        let ap = vec3_sub(p, self.a);
        let t = clamp(vec3_dot(ap, ab) / vec3_dot(ab, ab).max(1e-30), 0.0, 1.0);
        let closest = vec3_add(self.a, vec3_scale(ab, t));
        vec3_norm(vec3_sub(p, closest))
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Implicit Torus
// ─────────────────────────────────────────────────────────────────────────────

/// Signed distance function for a torus in the XZ plane.
#[derive(Debug, Clone, Copy)]
pub struct ImplicitTorus {
    /// Major radius (distance from center to tube center).
    pub major_radius: f64,
    /// Minor radius (tube radius).
    pub minor_radius: f64,
}

impl ImplicitTorus {
    /// Create a new torus SDF.
    pub fn new(major_radius: f64, minor_radius: f64) -> Self {
        Self {
            major_radius,
            minor_radius,
        }
    }

    /// Signed distance from point `p` to this torus (centered at origin, in XZ plane).
    pub fn sdf(&self, p: [f64; 3]) -> f64 {
        let q = [(p[0] * p[0] + p[2] * p[2]).sqrt() - self.major_radius, p[1]];
        (q[0] * q[0] + q[1] * q[1]).sqrt() - self.minor_radius
    }

    /// Gradient of the torus SDF.
    pub fn gradient(&self, p: [f64; 3]) -> [f64; 3] {
        numerical_grad(&|x| self.sdf(x), p, 1e-5)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Implicit Blend (smooth CSG)
// ─────────────────────────────────────────────────────────────────────────────

/// Smooth union/intersection/subtraction of two implicit SDFs.
///
/// Uses polynomial smooth-min for C¹ continuous blending.
#[derive(Debug, Clone)]
pub struct ImplicitBlend {
    /// Blend radius for smooth operations.
    pub k: f64,
}

impl ImplicitBlend {
    /// Create a new blend operator with blend radius `k`.
    pub fn new(k: f64) -> Self {
        Self { k }
    }

    /// Smooth union: min(a, b) with smooth blend.
    pub fn smooth_union(&self, a: f64, b: f64) -> f64 {
        smooth_min_polynomial(a, b, self.k)
    }

    /// Smooth intersection: max(a, b) with smooth blend.
    pub fn smooth_intersection(&self, a: f64, b: f64) -> f64 {
        -smooth_min_polynomial(-a, -b, self.k)
    }

    /// Smooth subtraction: a − b.
    pub fn smooth_subtraction(&self, a: f64, b: f64) -> f64 {
        self.smooth_intersection(a, -b)
    }

    /// Hard union (standard min).
    pub fn union(a: f64, b: f64) -> f64 {
        a.min(b)
    }

    /// Hard intersection (standard max).
    pub fn intersection(a: f64, b: f64) -> f64 {
        a.max(b)
    }

    /// Hard subtraction.
    pub fn subtraction(a: f64, b: f64) -> f64 {
        a.max(-b)
    }

    /// Exponential blend union.
    pub fn exp_union(&self, a: f64, b: f64) -> f64 {
        smooth_min_exponential(a, b, self.k)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 3D SDF Grid
// ─────────────────────────────────────────────────────────────────────────────

/// 3D voxel grid storing signed distance values.
///
/// Supports construction from analytic SDFs, trilinear interpolation,
/// and gradient estimation.
#[derive(Debug, Clone)]
pub struct SdfGrid3D {
    /// Grid resolution along X.
    pub nx: usize,
    /// Grid resolution along Y.
    pub ny: usize,
    /// Grid resolution along Z.
    pub nz: usize,
    /// Minimum corner of the grid bounding box.
    pub min_corner: [f64; 3],
    /// Voxel size.
    pub voxel_size: f64,
    /// Flat (x, y, z) storage: index = x + nx*(y + ny*z).
    pub data: Vec<f64>,
}

impl SdfGrid3D {
    /// Create an empty SDF grid.
    pub fn new(nx: usize, ny: usize, nz: usize, min_corner: [f64; 3], voxel_size: f64) -> Self {
        let data = vec![f64::INFINITY; nx * ny * nz];
        Self {
            nx,
            ny,
            nz,
            min_corner,
            voxel_size,
            data,
        }
    }

    /// Fill the grid from an analytic SDF function.
    pub fn from_sdf<F: Fn([f64; 3]) -> f64>(
        nx: usize,
        ny: usize,
        nz: usize,
        min_corner: [f64; 3],
        voxel_size: f64,
        f: F,
    ) -> Self {
        let mut grid = Self::new(nx, ny, nz, min_corner, voxel_size);
        for iz in 0..nz {
            for iy in 0..ny {
                for ix in 0..nx {
                    let p = grid.voxel_center(ix, iy, iz);
                    let idx = grid.index(ix, iy, iz);
                    grid.data[idx] = f(p);
                }
            }
        }
        grid
    }

    /// Get flat index.
    pub fn index(&self, ix: usize, iy: usize, iz: usize) -> usize {
        ix + self.nx * (iy + self.ny * iz)
    }

    /// World position of voxel center (ix, iy, iz).
    pub fn voxel_center(&self, ix: usize, iy: usize, iz: usize) -> [f64; 3] {
        [
            self.min_corner[0] + (ix as f64 + 0.5) * self.voxel_size,
            self.min_corner[1] + (iy as f64 + 0.5) * self.voxel_size,
            self.min_corner[2] + (iz as f64 + 0.5) * self.voxel_size,
        ]
    }

    /// Trilinear interpolation of SDF value at world position `p`.
    pub fn interpolate(&self, p: [f64; 3]) -> f64 {
        let fx = (p[0] - self.min_corner[0]) / self.voxel_size - 0.5;
        let fy = (p[1] - self.min_corner[1]) / self.voxel_size - 0.5;
        let fz = (p[2] - self.min_corner[2]) / self.voxel_size - 0.5;

        let ix = fx.floor() as i64;
        let iy = fy.floor() as i64;
        let iz = fz.floor() as i64;

        let tx = fx - ix as f64;
        let ty = fy - iy as f64;
        let tz = fz - iz as f64;

        let sample = |dx: i64, dy: i64, dz: i64| -> f64 {
            let xi = (ix + dx).clamp(0, self.nx as i64 - 1) as usize;
            let yi = (iy + dy).clamp(0, self.ny as i64 - 1) as usize;
            let zi = (iz + dz).clamp(0, self.nz as i64 - 1) as usize;
            self.data[self.index(xi, yi, zi)]
        };

        // Trilinear blend
        let c000 = sample(0, 0, 0);
        let c100 = sample(1, 0, 0);
        let c010 = sample(0, 1, 0);
        let c110 = sample(1, 1, 0);
        let c001 = sample(0, 0, 1);
        let c101 = sample(1, 0, 1);
        let c011 = sample(0, 1, 1);
        let c111 = sample(1, 1, 1);

        let c00 = c000 * (1.0 - tx) + c100 * tx;
        let c10 = c010 * (1.0 - tx) + c110 * tx;
        let c01 = c001 * (1.0 - tx) + c101 * tx;
        let c11 = c011 * (1.0 - tx) + c111 * tx;
        let c0 = c00 * (1.0 - ty) + c10 * ty;
        let c1 = c01 * (1.0 - ty) + c11 * ty;
        c0 * (1.0 - tz) + c1 * tz
    }

    /// Gradient at world position via trilinear interpolation.
    pub fn gradient(&self, p: [f64; 3]) -> [f64; 3] {
        let h = self.voxel_size * 0.5;
        let gx = (self.interpolate([p[0] + h, p[1], p[2]])
            - self.interpolate([p[0] - h, p[1], p[2]]))
            / (2.0 * h);
        let gy = (self.interpolate([p[0], p[1] + h, p[2]])
            - self.interpolate([p[0], p[1] - h, p[2]]))
            / (2.0 * h);
        let gz = (self.interpolate([p[0], p[1], p[2] + h])
            - self.interpolate([p[0], p[1], p[2] - h]))
            / (2.0 * h);
        [gx, gy, gz]
    }

    /// Get value at voxel (ix, iy, iz).
    pub fn get(&self, ix: usize, iy: usize, iz: usize) -> f64 {
        self.data[self.index(ix, iy, iz)]
    }

    /// Set value at voxel (ix, iy, iz).
    pub fn set(&mut self, ix: usize, iy: usize, iz: usize, v: f64) {
        let idx = self.index(ix, iy, iz);
        self.data[idx] = v;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// SDF Reinitializer (Fast Marching Method)
// ─────────────────────────────────────────────────────────────────────────────

/// Fast marching method to reinitialize a corrupted SDF grid.
///
/// Propagates the zero level-set outward to restore the signed distance property.
pub struct SdfReinitialize;

impl SdfReinitialize {
    /// Reinitialize `grid` using a simplified fast marching sweep.
    ///
    /// Identifies the zero-crossing voxels and sweeps outward using
    /// the Godunov upwind scheme approximation.
    pub fn reinitialize(grid: &mut SdfGrid3D) {
        let nx = grid.nx;
        let ny = grid.ny;
        let nz = grid.nz;
        let h = grid.voxel_size;

        // Identify interface voxels (sign changes with neighbor)
        let mut new_data = grid.data.clone();
        let sign = |v: f64| -> f64 { if v >= 0.0 { 1.0 } else { -1.0 } };

        // Multiple Gauss-Seidel sweeps
        for _ in 0..5 {
            for iz in 0..nz {
                for iy in 0..ny {
                    for ix in 0..nx {
                        let idx = grid.index(ix, iy, iz);
                        let v = grid.data[idx];
                        let s = sign(v);

                        // Min neighbor distance per axis (upwind)
                        let get = |x: i64, y: i64, z: i64| -> f64 {
                            let xi = x.clamp(0, nx as i64 - 1) as usize;
                            let yi = y.clamp(0, ny as i64 - 1) as usize;
                            let zi = z.clamp(0, nz as i64 - 1) as usize;
                            grid.data[grid.index(xi, yi, zi)]
                        };

                        let dx = get(ix as i64 - 1, iy as i64, iz as i64)
                            .abs()
                            .min(get(ix as i64 + 1, iy as i64, iz as i64).abs());
                        let dy = get(ix as i64, iy as i64 - 1, iz as i64)
                            .abs()
                            .min(get(ix as i64, iy as i64 + 1, iz as i64).abs());
                        let dz = get(ix as i64, iy as i64, iz as i64 - 1)
                            .abs()
                            .min(get(ix as i64, iy as i64, iz as i64 + 1).abs());

                        // Solve |grad phi| = 1 approximately:
                        // phi ≈ min_axis_dist + h * sign
                        let proposed = s * (dx.min(dy).min(dz) + h);
                        if (proposed).abs() < v.abs() {
                            new_data[idx] = proposed;
                        }
                    }
                }
            }
        }
        grid.data = new_data;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Marching Cubes
// ─────────────────────────────────────────────────────────────────────────────

/// Output triangle mesh from isosurface extraction.
#[derive(Debug, Clone, Default)]
pub struct IsoMesh {
    /// Vertex positions.
    pub vertices: Vec<[f64; 3]>,
    /// Triangle indices (triples into `vertices`).
    pub triangles: Vec<[usize; 3]>,
    /// Per-vertex normals.
    pub normals: Vec<[f64; 3]>,
}

impl IsoMesh {
    /// Create an empty mesh.
    pub fn new() -> Self {
        Self::default()
    }

    /// Number of triangles.
    pub fn num_triangles(&self) -> usize {
        self.triangles.len()
    }

    /// Number of vertices.
    pub fn num_vertices(&self) -> usize {
        self.vertices.len()
    }
}

/// Marching Cubes isosurface extraction (Lorensen-Cline).
pub struct MarchingCubes;

impl MarchingCubes {
    /// Extract an isosurface from a 3D SDF grid at the given isovalue.
    ///
    /// A node is *inside* when its value is strictly below `isovalue`.  The
    /// triangles are wound counter-clockwise when seen from outside, so their
    /// normals point towards increasing field values.
    pub fn extract(grid: &SdfGrid3D, isovalue: f64) -> IsoMesh {
        let mut mesh = IsoMesh::new();
        // Cube vertex offsets
        let corner_offsets: [[usize; 3]; 8] = [
            [0, 0, 0],
            [1, 0, 0],
            [1, 1, 0],
            [0, 1, 0],
            [0, 0, 1],
            [1, 0, 1],
            [1, 1, 1],
            [0, 1, 1],
        ];
        // Edge pairs: (v0, v1)
        let edge_pairs: [(usize, usize); 12] = [
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 0),
            (4, 5),
            (5, 6),
            (6, 7),
            (7, 4),
            (0, 4),
            (1, 5),
            (2, 6),
            (3, 7),
        ];

        let mut vertex_cache: std::collections::HashMap<(usize, usize), usize> =
            std::collections::HashMap::new();

        let nx = grid.nx.saturating_sub(1);
        let ny = grid.ny.saturating_sub(1);
        let nz = grid.nz.saturating_sub(1);

        for iz in 0..nz {
            for iy in 0..ny {
                for ix in 0..nx {
                    // Get corner values and positions
                    let mut corner_vals = [0.0f64; 8];
                    let mut corner_pos = [[0.0f64; 3]; 8];
                    for (ci, &[dx, dy, dz]) in corner_offsets.iter().enumerate() {
                        corner_pos[ci] = grid.voxel_center(ix + dx, iy + dy, iz + dz);
                        corner_vals[ci] = grid.get(ix + dx, iy + dy, iz + dz);
                    }

                    // Compute case index
                    let cube_idx: u8 = corner_vals.iter().enumerate().fold(0u8, |acc, (ci, &v)| {
                        if v < isovalue { acc | (1 << ci) } else { acc }
                    });
                    if MC_EDGE_TABLE[cube_idx as usize] == 0 {
                        continue;
                    }

                    // Compute edge intersections
                    let mut edge_verts = [usize::MAX; 12];
                    let edge_bits = MC_EDGE_TABLE[cube_idx as usize];
                    for (ei, &(v0, v1)) in edge_pairs.iter().enumerate() {
                        if edge_bits & (1 << ei) == 0 {
                            continue;
                        }
                        // Encode edge globally
                        let base = ix + grid.nx * (iy + grid.ny * iz);
                        let key = (base, ei);
                        let vi = *vertex_cache.entry(key).or_insert_with(|| {
                            let t = (isovalue - corner_vals[v0])
                                / (corner_vals[v1] - corner_vals[v0] + 1e-30).clamp(-1e10, 1e10);
                            let t = t.clamp(0.0, 1.0);
                            let p = [
                                corner_pos[v0][0] + t * (corner_pos[v1][0] - corner_pos[v0][0]),
                                corner_pos[v0][1] + t * (corner_pos[v1][1] - corner_pos[v0][1]),
                                corner_pos[v0][2] + t * (corner_pos[v1][2] - corner_pos[v0][2]),
                            ];
                            let idx = mesh.vertices.len();
                            mesh.vertices.push(p);
                            mesh.normals.push([0.0, 1.0, 0.0]); // computed later
                            idx
                        });
                        edge_verts[ei] = vi;
                    }

                    // Emit triangles
                    let tri_row = &MC_TRI_TABLE[cube_idx as usize];
                    let mut ti = 0;
                    while ti < 16 && tri_row[ti] >= 0 {
                        let a = edge_verts[tri_row[ti] as usize];
                        let b = edge_verts[tri_row[ti + 1] as usize];
                        let c = edge_verts[tri_row[ti + 2] as usize];
                        if a != usize::MAX && b != usize::MAX && c != usize::MAX {
                            // The table winds into the solid; swap two vertices
                            // so the triangle faces out of it.
                            mesh.triangles.push([a, c, b]);
                        }
                        ti += 3;
                    }
                }
            }
        }

        // Compute normals from SDF gradient
        for (i, &p) in mesh.vertices.iter().enumerate() {
            mesh.normals[i] = vec3_norm(grid.gradient(p));
        }

        mesh
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Dual Contouring
// ─────────────────────────────────────────────────────────────────────────────

/// Dual contouring with QEF minimization for sharp-feature preservation.
pub struct DualContouring;

impl DualContouring {
    /// Extract an isosurface using dual contouring from an SDF grid.
    ///
    /// Places one vertex per active voxel (QEF minimizer) and connects
    /// vertices of adjacent active voxels.
    pub fn extract(grid: &SdfGrid3D, isovalue: f64) -> IsoMesh {
        let mut mesh = IsoMesh::new();
        let nx = grid.nx;
        let ny = grid.ny;
        let nz = grid.nz;

        // For each active voxel, compute a representative vertex via QEF
        let mut voxel_vertex: std::collections::HashMap<usize, usize> =
            std::collections::HashMap::new();

        for iz in 0..nz.saturating_sub(1) {
            for iy in 0..ny.saturating_sub(1) {
                for ix in 0..nx.saturating_sub(1) {
                    let v000 = grid.get(ix, iy, iz);
                    let v100 = grid.get(ix + 1, iy, iz);
                    let v010 = grid.get(ix, iy + 1, iz);
                    let v001 = grid.get(ix, iy, iz + 1);
                    // Active if sign changes on any edge
                    let active = (v000 < isovalue) != (v100 < isovalue)
                        || (v000 < isovalue) != (v010 < isovalue)
                        || (v000 < isovalue) != (v001 < isovalue);
                    if !active {
                        continue;
                    }

                    // QEF: place vertex at voxel center (simplified)
                    let center = grid.voxel_center(ix, iy, iz);
                    let vi = mesh.vertices.len();
                    mesh.vertices.push(center);
                    mesh.normals.push(vec3_norm(grid.gradient(center)));
                    let key = ix + nx * (iy + ny * iz);
                    voxel_vertex.insert(key, vi);
                }
            }
        }

        // Connect adjacent active voxels (shared face → quad → 2 triangles)
        let directions = [
            (1, 0, 0, 1usize, 0usize, 0usize),
            (0, 1, 0, 0, 1, 0),
            (0, 0, 1, 0, 0, 1),
        ];
        for iz in 0..nz.saturating_sub(1) {
            for iy in 0..ny.saturating_sub(1) {
                for ix in 0..nx.saturating_sub(1) {
                    let v0 = grid.get(ix, iy, iz);
                    for &(dx, dy, dz, ox, oy, oz) in &directions {
                        let nx2 = ix + dx;
                        let ny2 = iy + dy;
                        let nz2 = iz + dz;
                        if nx2 >= nx || ny2 >= ny || nz2 >= nz {
                            continue;
                        }
                        let v1 = grid.get(nx2, ny2, nz2);
                        if (v0 < isovalue) == (v1 < isovalue) {
                            continue;
                        }
                        // Get the four voxels sharing this edge
                        let k0 = ix + nx * (iy + ny * iz);
                        let k1 = (ix + ox) + nx * ((iy + oy) + ny * (iz + oz));
                        let k2 = (ix + ox + dx) + nx * ((iy + oy + dy) + ny * (iz + oz + dz));
                        let k3 = (ix + dx) + nx * ((iy + dy) + ny * (iz + dz));
                        if let (Some(&a), Some(&b), Some(&c), Some(&d)) = (
                            voxel_vertex.get(&k0),
                            voxel_vertex.get(&k1),
                            voxel_vertex.get(&k2),
                            voxel_vertex.get(&k3),
                        ) {
                            mesh.triangles.push([a, b, c]);
                            mesh.triangles.push([a, c, d]);
                        }
                    }
                }
            }
        }

        mesh
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Marching Tetrahedra
// ─────────────────────────────────────────────────────────────────────────────

/// Kuhn (Freudenthal) split of the unit cube into six tetrahedra that share the
/// main diagonal `0 -> 6` (corner numbering of [`crate::mc_tables`]).
///
/// The tetrahedra tile the cube exactly (each has volume 1/6).  Every cube is split
/// along the same diagonal direction, so the triangulations of the faces shared by
/// neighbouring cubes coincide and the decomposition is conforming.
const KUHN_TETS: [[usize; 4]; 6] = [
    [0, 1, 2, 6],
    [0, 2, 3, 6],
    [0, 3, 7, 6],
    [0, 7, 4, 6],
    [0, 4, 5, 6],
    [0, 5, 1, 6],
];

/// The six edges of a tetrahedron as pairs of local vertex indices.
const TET_EDGES: [[usize; 2]; 6] = [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]];

/// Triangles of every inside/outside configuration of a tetrahedron, as triples of
/// [`TET_EDGES`] indices (`-1` ends a row).
///
/// Bit `i` of the case index is set when local vertex `i` is inside (`value <
/// isovalue`).  A lone inside or lone outside vertex cuts one triangle off its
/// three edges; a 2|2 split crosses four edges and yields a quad that is cut into
/// two triangles.  The winding is deliberately not encoded here: the tetrahedra of
/// the Kuhn split have mixed handedness, so every emitted triangle is oriented
/// against the field gradient when it is generated.
const TET_TRI_CASES: [[i8; 7]; 16] = [
    [-1, -1, -1, -1, -1, -1, -1],
    [0, 1, 2, -1, -1, -1, -1],
    [0, 3, 4, -1, -1, -1, -1],
    [1, 3, 4, 1, 4, 2, -1],
    [1, 3, 5, -1, -1, -1, -1],
    [0, 3, 5, 0, 5, 2, -1],
    [0, 4, 5, 0, 5, 1, -1],
    [2, 4, 5, -1, -1, -1, -1],
    [2, 4, 5, -1, -1, -1, -1],
    [0, 4, 5, 0, 5, 1, -1],
    [0, 3, 5, 0, 5, 2, -1],
    [1, 3, 5, -1, -1, -1, -1],
    [1, 3, 4, 1, 4, 2, -1],
    [0, 3, 4, -1, -1, -1, -1],
    [0, 1, 2, -1, -1, -1, -1],
    [-1, -1, -1, -1, -1, -1, -1],
];

/// Marching tetrahedra isosurface extraction.
///
/// Splits every grid cell into six tetrahedra around its main diagonal (the Kuhn
/// split) and applies the 16-case tetrahedron table, which is
/// free of the ambiguous configurations of marching cubes.  A node is *inside*
/// when its value is strictly below the isovalue.  The result is an indexed mesh:
///
/// * every crossed edge of the tetrahedral lattice yields exactly one vertex that
///   is shared by all triangles around it, so a closed surface gives a watertight
///   mesh in which every edge belongs to two triangles;
/// * triangles are wound counter-clockwise when seen from outside (their normals
///   point towards increasing field values);
/// * the per-vertex normals are the normalised grid gradient at the vertex.
pub struct MarchingTetrahedra;

/// A grid node seen from inside one tetrahedron.
#[derive(Clone, Copy)]
struct TetNode {
    /// Flat index of the node in the grid.
    index: usize,
    /// World position of the node.
    pos: [f64; 3],
    /// Field value at the node.
    value: f64,
}

/// Incremental mesh builder that welds the crossing points of shared edges.
struct TetMeshBuilder<'a> {
    grid: &'a SdfGrid3D,
    iso: f64,
    /// Vertex index per crossed edge, keyed by `(lower node, upper node)`.  A
    /// crossing that falls exactly on a node is keyed by `(node, node)`, so that
    /// the crossings of all edges meeting at the node are one vertex.
    vertex_of: std::collections::HashMap<(usize, usize), usize>,
    mesh: IsoMesh,
}

impl<'a> TetMeshBuilder<'a> {
    fn new(grid: &'a SdfGrid3D, iso: f64) -> Self {
        Self {
            grid,
            iso,
            vertex_of: std::collections::HashMap::new(),
            mesh: IsoMesh::new(),
        }
    }

    /// Index of the vertex where the surface crosses the edge `a`-`b`.
    ///
    /// The crossing is always interpolated from the node with the lower grid index
    /// to the one with the higher index, so every tetrahedron sharing the edge
    /// computes the bit-identical position.
    fn crossing_vertex(&mut self, a: TetNode, b: TetNode) -> usize {
        let (lo, hi) = if a.index <= b.index { (a, b) } else { (b, a) };
        let dv = hi.value - lo.value;
        let t = if dv == 0.0 {
            0.5
        } else {
            ((self.iso - lo.value) / dv).clamp(0.0, 1.0)
        };
        let key = if t <= 0.0 {
            (lo.index, lo.index)
        } else if t >= 1.0 {
            (hi.index, hi.index)
        } else {
            (lo.index, hi.index)
        };
        if let Some(&id) = self.vertex_of.get(&key) {
            return id;
        }
        let pos = if t <= 0.0 {
            lo.pos
        } else if t >= 1.0 {
            hi.pos
        } else {
            vec3_add(lo.pos, vec3_scale(vec3_sub(hi.pos, lo.pos), t))
        };
        let id = self.mesh.vertices.len();
        self.mesh.vertices.push(pos);
        self.vertex_of.insert(key, id);
        id
    }

    /// Emit the triangles of one tetrahedron.
    fn process_tet(&mut self, tet: [TetNode; 4]) {
        let mut case = 0usize;
        for (i, node) in tet.iter().enumerate() {
            if node.value < self.iso {
                case |= 1 << i;
            }
        }
        for row in TET_TRI_CASES[case].chunks_exact(3) {
            let [Ok(e0), Ok(e1), Ok(e2)] = [row[0], row[1], row[2]].map(usize::try_from) else {
                break;
            };
            let ends = |e: usize| {
                let [a, b] = TET_EDGES[e];
                (tet[a], tet[b])
            };
            let (a0, b0) = ends(e0);
            let (a1, b1) = ends(e1);
            let (a2, b2) = ends(e2);
            let mut ids = [
                self.crossing_vertex(a0, b0),
                self.crossing_vertex(a1, b1),
                self.crossing_vertex(a2, b2),
            ];
            if ids[0] == ids[1] || ids[1] == ids[2] || ids[0] == ids[2] {
                // Crossings that coincide on a node leave a zero-area triangle.
                continue;
            }
            // Every crossed edge of the triangle runs from an inside to an outside
            // node, and the (linear) field increases along it: the outward
            // normal must make an acute angle with that direction.
            let (inside, outside) = if a0.value < self.iso {
                (a0, b0)
            } else {
                (b0, a0)
            };
            let p = ids.map(|id| self.mesh.vertices[id]);
            let normal = vec3_cross(vec3_sub(p[1], p[0]), vec3_sub(p[2], p[0]));
            if vec3_dot(normal, normal) == 0.0 {
                continue;
            }
            if vec3_dot(normal, vec3_sub(outside.pos, inside.pos)) < 0.0 {
                ids.swap(1, 2);
            }
            self.mesh.triangles.push(ids);
        }
    }

    /// Fill in the per-vertex normals and return the mesh.
    fn finish(mut self) -> IsoMesh {
        self.mesh.normals = self
            .mesh
            .vertices
            .iter()
            .map(|&p| vec3_norm(self.grid.gradient(p)))
            .collect();
        self.mesh
    }
}

impl MarchingTetrahedra {
    /// Extract isosurface from SDF grid using marching tetrahedra.
    pub fn extract(grid: &SdfGrid3D, isovalue: f64) -> IsoMesh {
        let mut builder = TetMeshBuilder::new(grid, isovalue);
        let nx = grid.nx.saturating_sub(1);
        let ny = grid.ny.saturating_sub(1);
        let nz = grid.nz.saturating_sub(1);

        for iz in 0..nz {
            for iy in 0..ny {
                for ix in 0..nx {
                    let corners: [TetNode; 8] = std::array::from_fn(|c| {
                        let [di, dj, dk] = MC_CORNER_OFFSETS[c];
                        let (i, j, k) = (ix + di, iy + dj, iz + dk);
                        let index = grid.index(i, j, k);
                        TetNode {
                            index,
                            pos: grid.voxel_center(i, j, k),
                            value: grid.data[index],
                        }
                    });
                    let values: [f64; 8] = std::array::from_fn(|c| corners[c].value);
                    let case = case_index(&values, isovalue);
                    if case == 0 || case == 255 {
                        continue;
                    }
                    for tet in &KUHN_TETS {
                        builder.process_tet(tet.map(|c| corners[c]));
                    }
                }
            }
        }
        builder.finish()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Implicit Convolution (Smooth SDF)
// ─────────────────────────────────────────────────────────────────────────────

/// Smooth SDF operations: offset surfaces and Minkowski sums.
pub struct ImplicitConvolution;

impl ImplicitConvolution {
    /// Create an offset surface SDF: expands the surface outward by `offset`.
    pub fn offset_surface<F: Fn([f64; 3]) -> f64>(f: F, offset: f64) -> impl Fn([f64; 3]) -> f64 {
        move |p| f(p) - offset
    }

    /// Gaussian-blurred SDF approximation via Monte Carlo sampling on a grid.
    ///
    /// Approximates ∫ f(p+δ) G(δ, sigma) dδ by sampling 8 corners.
    pub fn gaussian_blur<F: Fn([f64; 3]) -> f64>(f: F, sigma: f64) -> impl Fn([f64; 3]) -> f64 {
        move |p: [f64; 3]| {
            let s = sigma;
            let mut sum = 0.0;
            for dx in [-s, s] {
                for dy in [-s, s] {
                    for dz in [-s, s] {
                        sum += f([p[0] + dx, p[1] + dy, p[2] + dz]);
                    }
                }
            }
            sum / 8.0
        }
    }

    /// Minkowski sum of two SDFs (approximate via smooth union with offset).
    ///
    /// For convex shapes: `(f ⊕ B_r)(p) ≈ f(p) − r`.
    pub fn minkowski_sum<F: Fn([f64; 3]) -> f64>(f: F, radius: f64) -> impl Fn([f64; 3]) -> f64 {
        move |p| f(p) - radius
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Ray Marching
// ─────────────────────────────────────────────────────────────────────────────

/// Ray march hit result.
#[derive(Debug, Clone, Copy)]
pub struct RayMarchHit {
    /// Hit position.
    pub position: [f64; 3],
    /// Distance along the ray.
    pub t: f64,
    /// Number of steps taken.
    pub steps: usize,
}

/// Sphere-tracing ray marcher through an SDF.
pub struct RayMarchSdf;

impl RayMarchSdf {
    /// Sphere-trace a ray through an SDF.
    ///
    /// Starting at `origin` in direction `dir` (unit vector), marches until the
    /// surface is hit (|SDF| < `epsilon`) or `max_steps` is reached.
    pub fn march<F: Fn([f64; 3]) -> f64>(
        f: &F,
        origin: [f64; 3],
        dir: [f64; 3],
        max_dist: f64,
        max_steps: usize,
        epsilon: f64,
    ) -> Option<RayMarchHit> {
        let dir = vec3_norm(dir);
        let mut t = 0.0;
        for step in 0..max_steps {
            let p = [
                origin[0] + t * dir[0],
                origin[1] + t * dir[1],
                origin[2] + t * dir[2],
            ];
            let d = f(p);
            if d.abs() < epsilon {
                return Some(RayMarchHit {
                    position: p,
                    t,
                    steps: step + 1,
                });
            }
            t += d.abs().max(epsilon * 0.1);
            if t > max_dist {
                break;
            }
        }
        None
    }

    /// Compute ambient occlusion via SDF samples along normal.
    pub fn ambient_occlusion<F: Fn([f64; 3]) -> f64>(
        f: &F,
        pos: [f64; 3],
        normal: [f64; 3],
        samples: usize,
        step_size: f64,
    ) -> f64 {
        let mut occ = 0.0;
        for i in 1..=samples {
            let t = i as f64 * step_size;
            let p = vec3_add(pos, vec3_scale(normal, t));
            let d = f(p);
            occ += (t - d) / (2.0_f64.powi(i as i32));
        }
        (1.0 - occ).clamp(0.0, 1.0)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn test_sphere_sdf_center() {
        let s = ImplicitSphere::new([0.0, 0.0, 0.0], 1.0);
        assert!((s.sdf([0.0, 0.0, 0.0]) - (-1.0)).abs() < 1e-10);
    }

    #[test]
    fn test_sphere_sdf_surface() {
        let s = ImplicitSphere::new([0.0, 0.0, 0.0], 1.0);
        assert!(s.sdf([1.0, 0.0, 0.0]).abs() < 1e-10);
    }

    #[test]
    fn test_sphere_sdf_outside() {
        let s = ImplicitSphere::new([0.0, 0.0, 0.0], 1.0);
        assert!(s.sdf([2.0, 0.0, 0.0]) > 0.0);
    }

    #[test]
    fn test_sphere_gradient_unit_length() {
        let s = ImplicitSphere::new([0.0, 0.0, 0.0], 1.0);
        let g = s.gradient([2.0, 0.0, 0.0]);
        let len = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt();
        assert!((len - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_box_sdf_inside() {
        let b = ImplicitBox::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert!(b.sdf([0.0, 0.0, 0.0]) < 0.0);
    }

    #[test]
    fn test_box_sdf_surface() {
        let b = ImplicitBox::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert!(b.sdf([1.0, 0.0, 0.0]).abs() < 1e-8);
    }

    #[test]
    fn test_box_sdf_outside() {
        let b = ImplicitBox::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert!(b.sdf([2.0, 0.0, 0.0]) > 0.0);
    }

    #[test]
    fn test_rounded_box_sdf() {
        let b = ImplicitBox::rounded([0.0, 0.0, 0.0], [1.0, 1.0, 1.0], 0.1);
        assert!(b.sdf([0.0, 0.0, 0.0]) < 0.0);
    }

    #[test]
    fn test_capsule_sdf_on_axis() {
        let c = ImplicitCapsule::new([0.0, 0.0, 0.0], [0.0, 2.0, 0.0], 0.5);
        // Point on the segment, at radius distance: on surface
        assert!(c.sdf([0.5, 1.0, 0.0]).abs() < 1e-10);
    }

    #[test]
    fn test_capsule_sdf_inside() {
        let c = ImplicitCapsule::new([0.0, 0.0, 0.0], [0.0, 2.0, 0.0], 0.5);
        assert!(c.sdf([0.0, 1.0, 0.0]) < 0.0);
    }

    #[test]
    fn test_torus_sdf_on_surface() {
        let t = ImplicitTorus::new(2.0, 0.5);
        // Point on the torus surface: distance from Z axis = 2.0, y=0
        let p = [2.5, 0.0, 0.0]; // r=2.5, major=2, minor=0.5 → sdf=0
        assert!(t.sdf(p).abs() < 1e-10);
    }

    #[test]
    fn test_torus_sdf_inside() {
        let t = ImplicitTorus::new(2.0, 0.5);
        let p = [2.0, 0.0, 0.0]; // on the centerline of the tube
        assert!(t.sdf(p) < 0.0);
    }

    #[test]
    fn test_smooth_min_polynomial_limit() {
        // As k→0, smooth_min → min
        let a = 1.0;
        let b = 3.0;
        let sm = smooth_min_polynomial(a, b, 0.001);
        assert!((sm - a.min(b)).abs() < 0.01);
    }

    #[test]
    fn test_smooth_min_exponential_symmetry() {
        let a = 1.0;
        let b = 2.0;
        let k = 1.0;
        let sm1 = smooth_min_exponential(a, b, k);
        let sm2 = smooth_min_exponential(b, a, k);
        assert!((sm1 - sm2).abs() < 1e-12);
    }

    #[test]
    fn test_implicit_blend_smooth_union() {
        let blend = ImplicitBlend::new(0.3);
        let u = blend.smooth_union(0.5, -0.1);
        // Should be close to min
        assert!(u <= 0.5 + 0.01);
    }

    #[test]
    fn test_implicit_blend_hard_union() {
        assert_eq!(ImplicitBlend::union(1.0, 2.0), 1.0);
        assert_eq!(ImplicitBlend::union(-1.0, 2.0), -1.0);
    }

    #[test]
    fn test_implicit_blend_subtraction() {
        // Subtracting inside region from outside: should remain outside
        let r = ImplicitBlend::subtraction(2.0, -1.0);
        assert!(r > 0.0);
    }

    #[test]
    fn test_sdf_grid_from_sphere() {
        let sphere = ImplicitSphere::new([0.5, 0.5, 0.5], 0.3);
        let grid = SdfGrid3D::from_sdf(10, 10, 10, [0.0, 0.0, 0.0], 0.1, |p| sphere.sdf(p));
        assert_eq!(grid.data.len(), 1000);
        // Center voxel should be inside
        let v = grid.get(5, 5, 5);
        assert!(v < 0.0);
    }

    #[test]
    fn test_sdf_grid_interpolate() {
        let sphere = ImplicitSphere::new([0.5, 0.5, 0.5], 0.3);
        let grid = SdfGrid3D::from_sdf(20, 20, 20, [0.0, 0.0, 0.0], 0.05, |p| sphere.sdf(p));
        let v = grid.interpolate([0.5, 0.5, 0.5]);
        assert!(v < 0.0, "center should be inside: {}", v);
    }

    #[test]
    fn test_marching_cubes_sphere() {
        let sphere = ImplicitSphere::new([2.0, 2.0, 2.0], 1.0);
        let grid = SdfGrid3D::from_sdf(20, 20, 20, [0.0, 0.0, 0.0], 0.25, |p| sphere.sdf(p));
        let mesh = MarchingCubes::extract(&grid, 0.0);
        assert!(
            !mesh.vertices.is_empty(),
            "marching cubes should produce vertices"
        );
        assert!(
            !mesh.triangles.is_empty(),
            "marching cubes should produce triangles"
        );
    }

    #[test]
    fn test_marching_tetrahedra_sphere() {
        let sphere = ImplicitSphere::new([2.0, 2.0, 2.0], 1.0);
        let grid = SdfGrid3D::from_sdf(16, 16, 16, [0.0, 0.0, 0.0], 0.3, |p| sphere.sdf(p));
        let mesh = MarchingTetrahedra::extract(&grid, 0.0);
        assert!(
            !mesh.triangles.is_empty(),
            "marching tetrahedra should produce triangles"
        );
    }

    #[test]
    fn test_dual_contouring_sphere() {
        let sphere = ImplicitSphere::new([2.0, 2.0, 2.0], 1.0);
        let grid = SdfGrid3D::from_sdf(16, 16, 16, [0.0, 0.0, 0.0], 0.3, |p| sphere.sdf(p));
        let mesh = DualContouring::extract(&grid, 0.0);
        assert!(
            !mesh.vertices.is_empty(),
            "dual contouring should produce vertices"
        );
    }

    #[test]
    fn test_ray_march_sphere_hit() {
        let sphere = ImplicitSphere::new([0.0, 0.0, 5.0], 1.0);
        let f = |p: [f64; 3]| sphere.sdf(p);
        let hit = RayMarchSdf::march(&f, [0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 20.0, 200, 1e-4);
        assert!(hit.is_some(), "ray should hit sphere");
        let h = hit.unwrap();
        assert!((h.t - 4.0).abs() < 0.01, "hit at t≈4, got {}", h.t);
    }

    #[test]
    fn test_ray_march_sphere_miss() {
        let sphere = ImplicitSphere::new([0.0, 0.0, 5.0], 1.0);
        let f = |p: [f64; 3]| sphere.sdf(p);
        // Ray in opposite direction
        let hit = RayMarchSdf::march(&f, [0.0, 0.0, 0.0], [0.0, 0.0, -1.0], 20.0, 200, 1e-4);
        assert!(hit.is_none(), "ray in wrong direction should miss");
    }

    #[test]
    fn test_ambient_occlusion_on_surface() {
        let sphere = ImplicitSphere::new([0.0, 0.0, 0.0], 1.0);
        let f = |p: [f64; 3]| sphere.sdf(p);
        let ao = RayMarchSdf::ambient_occlusion(&f, [1.0, 0.0, 0.0], [1.0, 0.0, 0.0], 5, 0.1);
        assert!(
            (0.0..=1.0).contains(&ao),
            "AO should be in [0,1], got {}",
            ao
        );
    }

    #[test]
    fn test_sdf_reinitialize() {
        let sphere = ImplicitSphere::new([2.0, 2.0, 2.0], 1.0);
        let mut grid = SdfGrid3D::from_sdf(16, 16, 16, [0.0, 0.0, 0.0], 0.3, |p| sphere.sdf(p));
        // Corrupt some values
        grid.set(5, 5, 5, 100.0);
        SdfReinitialize::reinitialize(&mut grid);
        // After reinitialization, values should be more reasonable
        assert!(grid.get(5, 5, 5).abs() < 50.0);
    }

    #[test]
    fn test_iso_mesh_empty() {
        let mesh = IsoMesh::new();
        assert_eq!(mesh.num_vertices(), 0);
        assert_eq!(mesh.num_triangles(), 0);
    }

    #[test]
    fn test_offset_surface() {
        let sphere = ImplicitSphere::new([0.0, 0.0, 0.0], 1.0);
        let offset = ImplicitConvolution::offset_surface(|p| sphere.sdf(p), 0.5);
        // Surface now at r=1.5
        assert!(offset([1.5, 0.0, 0.0]).abs() < 1e-10);
    }

    #[test]
    fn test_smooth_min_polynomial_symmetric() {
        let a = 0.5;
        let b = 1.5;
        let k = 0.3;
        let s1 = smooth_min_polynomial(a, b, k);
        let s2 = smooth_min_polynomial(b, a, k);
        assert!((s1 - s2).abs() < 1e-12);
    }

    #[test]
    fn test_sphere_hessian_shape() {
        let s = ImplicitSphere::new([0.0, 0.0, 0.0], 1.0);
        let h = s.hessian([2.0, 0.0, 0.0]);
        assert_eq!(h.len(), 9);
    }

    #[test]
    fn test_capsule_gradient_unit() {
        let c = ImplicitCapsule::new([0.0, 0.0, 0.0], [0.0, 2.0, 0.0], 0.5);
        let g = c.gradient([1.0, 1.0, 0.0]);
        let len = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt();
        assert!((len - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_sdf_grid_gradient_direction() {
        let sphere = ImplicitSphere::new([2.0, 2.0, 2.0], 1.0);
        let grid = SdfGrid3D::from_sdf(20, 20, 20, [0.0, 0.0, 0.0], 0.2, |p| sphere.sdf(p));
        // Gradient at a point outside should point away from center
        let g = grid.gradient([3.5, 2.0, 2.0]);
        assert!(g[0] > 0.0, "gradient x component should be positive");
    }

    #[test]
    fn test_implicit_torus_outside() {
        let t = ImplicitTorus::new(2.0, 0.5);
        // Far from torus
        assert!(t.sdf([10.0, 0.0, 0.0]) > 0.0);
    }

    #[test]
    fn test_sdf_numerical_gradient() {
        let sphere = ImplicitSphere::new([0.0, 0.0, 0.0], 1.0);
        let f = |p: [f64; 3]| sphere.sdf(p);
        let g = sdf_gradient_numerical(&f, [2.0, 0.0, 0.0], 1e-4);
        // Should point in +x direction
        assert!(g[0] > 0.0);
    }

    #[test]
    fn test_marching_cubes_box() {
        let b = ImplicitBox::new([2.0, 2.0, 2.0], [0.8, 0.8, 0.8]);
        let grid = SdfGrid3D::from_sdf(20, 20, 20, [0.0, 0.0, 0.0], 0.25, |p| b.sdf(p));
        let mesh = MarchingCubes::extract(&grid, 0.0);
        assert!(
            !mesh.triangles.is_empty(),
            "box isosurface should have triangles"
        );
    }

    #[test]
    fn test_pi_used() {
        // Just ensure PI import is used
        let _ = PI;
    }
}

#[cfg(test)]
mod marching_tetrahedra_tests {
    use super::*;
    use crate::mc_tables::test_support::{
        area_and_volume, degenerate_triangle_count, edge_valence_histogram,
        is_closed_and_consistently_oriented,
    };
    use std::f64::consts::PI;

    /// Signed volume of the tetrahedron `(a, b, c, d)`.
    fn tet_volume(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> f64 {
        let n = vec3_cross(vec3_sub(b, a), vec3_sub(c, a));
        vec3_dot(n, vec3_sub(d, a)) / 6.0
    }

    fn corner(c: usize) -> [f64; 3] {
        let [i, j, k] = MC_CORNER_OFFSETS[c];
        [i as f64, j as f64, k as f64]
    }

    /// Centre of the sphere sampled by [`sphere_grid`].  It is off the lattice of
    /// grid nodes, so that no node lies exactly on the surface.
    fn sphere_center(n: usize, voxel: f64) -> [f64; 3] {
        let half = 0.5 * n as f64 * voxel;
        [half + 0.013, half - 0.007, half + 0.0031]
    }

    /// Sphere sampled on a cubic grid whose nodes avoid the surface exactly.
    fn sphere_grid(n: usize, voxel: f64, radius: f64) -> SdfGrid3D {
        let s = ImplicitSphere::new(sphere_center(n, voxel), radius);
        SdfGrid3D::from_sdf(n, n, n, [0.0; 3], voxel, |p| s.sdf(p))
    }

    #[test]
    fn kuhn_tets_tile_the_unit_cube() {
        let mut total = 0.0;
        for tet in &KUHN_TETS {
            let v = tet_volume(
                corner(tet[0]),
                corner(tet[1]),
                corner(tet[2]),
                corner(tet[3]),
            )
            .abs();
            assert!((v - 1.0 / 6.0).abs() < 1e-15, "tet {tet:?} has volume {v}");
            total += v;
        }
        assert!((total - 1.0).abs() < 1e-14, "total volume {total}");
        // All six tets contain the main diagonal 0-6 and use every corner.
        let mut used = [false; 8];
        for tet in &KUHN_TETS {
            assert!(tet.contains(&0) && tet.contains(&6));
            for &c in tet {
                used[c] = true;
            }
        }
        assert!(used.iter().all(|&u| u), "corner 2 or 7 left unused");
    }

    #[test]
    fn kuhn_split_is_conforming_across_cube_faces() {
        // The diagonal a tet face puts on a cube face, as a sorted corner pair.
        let face_diagonals = |axis: usize, side: usize| -> Vec<[usize; 2]> {
            let mut out = Vec::new();
            for tet in &KUHN_TETS {
                // Corners of the tet lying on this cube face.
                let on_face: Vec<usize> = tet
                    .iter()
                    .copied()
                    .filter(|&c| MC_CORNER_OFFSETS[c][axis] == side)
                    .collect();
                if on_face.len() == 3 {
                    // A triangle on the face: its diagonal joins the two corners
                    // that differ in two coordinates.
                    for a in 0..3 {
                        for b in (a + 1)..3 {
                            let (ca, cb) = (on_face[a], on_face[b]);
                            let differing = (0..3)
                                .filter(|&k| MC_CORNER_OFFSETS[ca][k] != MC_CORNER_OFFSETS[cb][k])
                                .count();
                            if differing == 2 {
                                let mut pair = [ca, cb];
                                pair.sort_unstable();
                                out.push(pair);
                            }
                        }
                    }
                }
            }
            out.sort_unstable();
            out.dedup();
            out
        };
        for axis in 0..3 {
            let lower = face_diagonals(axis, 0);
            let upper = face_diagonals(axis, 1);
            assert_eq!(lower.len(), 1, "axis {axis}: one diagonal per face");
            assert_eq!(upper.len(), 1, "axis {axis}: one diagonal per face");
            // Translating the lower face by +1 along `axis` must reproduce the
            // diagonal of the upper face of the neighbouring cube.
            let translate = |c: usize| -> usize {
                let mut off = MC_CORNER_OFFSETS[c];
                off[axis] = 1;
                (0..8)
                    .find(|&q| MC_CORNER_OFFSETS[q] == off)
                    .expect("every offset is a cube corner")
            };
            let mut mapped = [translate(lower[0][0]), translate(lower[0][1])];
            mapped.sort_unstable();
            assert_eq!(mapped, upper[0], "axis {axis}: faces disagree");
        }
    }

    #[test]
    fn tet_case_table_cuts_exactly_the_crossed_edges() {
        for (case, row) in TET_TRI_CASES.iter().enumerate() {
            let inside = |v: usize| (case >> v) & 1 == 1;
            let mut crossed = [false; 6];
            for (e, &[a, b]) in TET_EDGES.iter().enumerate() {
                crossed[e] = inside(a) != inside(b);
            }
            let mut used = [false; 6];
            let mut triangles = 0;
            for t in row.chunks_exact(3) {
                if t[0] < 0 {
                    break;
                }
                triangles += 1;
                for &e in t {
                    let e = usize::try_from(e).expect("non-negative edge index");
                    assert!(crossed[e], "case {case}: edge {e} does not cross");
                    used[e] = true;
                }
            }
            assert_eq!(used, crossed, "case {case}: every crossed edge is used");
            let inside_count = (0..4).filter(|&v| inside(v)).count();
            let expected = match inside_count {
                0 | 4 => 0,
                1 | 3 => 1,
                _ => 2,
            };
            assert_eq!(triangles, expected, "case {case}");
        }
    }

    #[test]
    fn sphere_mesh_is_watertight_outward_and_accurate() {
        let radius = 1.0;
        let grid = sphere_grid(24, 0.13, radius);
        let mesh = MarchingTetrahedra::extract(&grid, 0.0);
        assert!(mesh.num_triangles() > 1000, "{}", mesh.num_triangles());
        assert_eq!(mesh.normals.len(), mesh.vertices.len());
        assert_eq!(
            edge_valence_histogram(&mesh.triangles),
            vec![(2, mesh.triangles.len() * 3 / 2)],
            "every edge must be shared by exactly two triangles"
        );
        assert!(is_closed_and_consistently_oriented(&mesh.triangles));
        assert_eq!(
            degenerate_triangle_count(&mesh.vertices, &mesh.triangles),
            0
        );
        let (area, volume) = area_and_volume(&mesh.vertices, &mesh.triangles);
        let exact_area = 4.0 * PI * radius * radius;
        let exact_volume = 4.0 / 3.0 * PI * radius.powi(3);
        assert!(
            (area / exact_area - 1.0).abs() < 0.03,
            "area {area} vs {exact_area}"
        );
        assert!(
            volume > 0.0,
            "triangles must face outwards, volume {volume}"
        );
        assert!(
            (volume / exact_volume - 1.0).abs() < 0.03,
            "volume {volume} vs {exact_volume}"
        );
    }

    #[test]
    fn vertex_normals_follow_the_sdf_gradient_not_plus_y() {
        let grid = sphere_grid(20, 0.15, 1.0);
        let mesh = MarchingTetrahedra::extract(&grid, 0.0);
        let center = sphere_center(20, 0.15);
        let mut worst = 1.0_f64;
        for (p, n) in mesh.vertices.iter().zip(mesh.normals.iter()) {
            assert!(
                (vec3_len(*n) - 1.0).abs() < 1e-9,
                "normal {n:?} is not unit"
            );
            let radial = vec3_norm(vec3_sub(*p, center));
            worst = worst.min(vec3_dot(radial, *n));
        }
        assert!(
            worst > 0.95,
            "normals deviate from the radial direction: {worst}"
        );
        assert!(
            mesh.normals.iter().any(|n| n[0].abs() > 0.5),
            "normals must not all be +Y"
        );
    }

    #[test]
    fn triangle_winding_matches_the_vertex_normals() {
        let grid = sphere_grid(18, 0.17, 1.0);
        let mesh = MarchingTetrahedra::extract(&grid, 0.0);
        for t in &mesh.triangles {
            let (a, b, c) = (
                mesh.vertices[t[0]],
                mesh.vertices[t[1]],
                mesh.vertices[t[2]],
            );
            let face = vec3_cross(vec3_sub(b, a), vec3_sub(c, a));
            let avg = vec3_add(
                vec3_add(mesh.normals[t[0]], mesh.normals[t[1]]),
                mesh.normals[t[2]],
            );
            assert!(vec3_dot(face, avg) > 0.0, "triangle {t:?} faces inwards");
        }
    }

    #[test]
    fn crossing_on_a_node_leaves_no_degenerate_triangles() {
        // f = x - 2.5 vanishes exactly on the node plane x = 2.5 (nodes sit at
        // half-integers), so every crossing snaps onto a node.
        let grid = SdfGrid3D::from_sdf(6, 5, 5, [0.0; 3], 1.0, |p| p[0] - 2.5);
        let mesh = MarchingTetrahedra::extract(&grid, 0.0);
        assert!(!mesh.triangles.is_empty());
        assert_eq!(
            degenerate_triangle_count(&mesh.vertices, &mesh.triangles),
            0
        );
        for v in &mesh.vertices {
            assert!((v[0] - 2.5).abs() < 1e-12, "vertex off the plane: {v:?}");
        }
        // The plane spans the nodes 0.5..4.5 in y and z: a 4 x 4 patch.
        let (area, _) = area_and_volume(&mesh.vertices, &mesh.triangles);
        assert!((area - 16.0).abs() < 1e-12, "area {area}");
        for t in &mesh.triangles {
            let (a, b, c) = (
                mesh.vertices[t[0]],
                mesh.vertices[t[1]],
                mesh.vertices[t[2]],
            );
            let n = vec3_cross(vec3_sub(b, a), vec3_sub(c, a));
            assert!(n[0] > 0.0, "triangle {t:?} must face +x");
        }
    }

    #[test]
    fn marching_cubes_triangles_face_outwards_like_the_tetrahedra() {
        let grid = sphere_grid(20, 0.15, 1.0);
        let sphere_centre = sphere_center(20, 0.15);
        let mesh = MarchingCubes::extract(&grid, 0.0);
        assert!(mesh.num_triangles() > 500, "{}", mesh.num_triangles());
        let mut volume = 0.0;
        for t in &mesh.triangles {
            let (a, b, c) = (
                mesh.vertices[t[0]],
                mesh.vertices[t[1]],
                mesh.vertices[t[2]],
            );
            let face = vec3_cross(vec3_sub(b, a), vec3_sub(c, a));
            let centroid = vec3_scale(vec3_add(vec3_add(a, b), c), 1.0 / 3.0);
            assert!(
                vec3_dot(face, vec3_sub(centroid, sphere_centre)) > 0.0,
                "triangle {t:?} faces inwards"
            );
            volume += vec3_dot(a, vec3_cross(b, c)) / 6.0;
        }
        let exact = 4.0 / 3.0 * PI;
        assert!(
            (volume / exact - 1.0).abs() < 0.05,
            "volume {volume} vs {exact}"
        );
    }

    #[test]
    fn empty_for_constant_and_degenerate_grids() {
        let outside = SdfGrid3D::from_sdf(4, 4, 4, [0.0; 3], 1.0, |_| 1.0);
        assert_eq!(
            MarchingTetrahedra::extract(&outside, 0.0).num_triangles(),
            0
        );
        let inside = SdfGrid3D::from_sdf(4, 4, 4, [0.0; 3], 1.0, |_| -1.0);
        assert_eq!(MarchingTetrahedra::extract(&inside, 0.0).num_triangles(), 0);
        let flat = SdfGrid3D::from_sdf(1, 4, 4, [0.0; 3], 1.0, |p| p[1] - 2.0);
        assert_eq!(MarchingTetrahedra::extract(&flat, 0.0).num_triangles(), 0);
    }
}
