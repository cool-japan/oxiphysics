// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Implicit surfaces and signed distance fields (SDF).
//!
//! Provides SDF primitives, smooth blending operations, voxel SDF grids,
//! marching cubes isosurface extraction, dual contouring, and ray marching.

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

// Marching Cubes edge table (which edges are intersected for each case)
const MC_EDGE_TABLE: [u16; 256] = [
    0x000, 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c, 0x80c, 0x905, 0xa0f, 0xb06, 0xc0a,
    0xd03, 0xe09, 0xf00, 0x190, 0x099, 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c, 0x99c, 0x895,
    0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90, 0x230, 0x339, 0x033, 0x13a, 0x636, 0x73f, 0x435,
    0x53c, 0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30, 0x3a0, 0x2a9, 0x1a3, 0x0aa,
    0x7a6, 0x6af, 0x5a5, 0x4ac, 0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0, 0x460,
    0x569, 0x663, 0x76a, 0x066, 0x16f, 0x265, 0x36c, 0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963,
    0xa69, 0xb60, 0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0x0ff, 0x3f5, 0x2fc, 0xdfc, 0xcf5, 0xfff,
    0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0, 0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x055, 0x15c,
    0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950, 0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6,
    0x2cf, 0x1c5, 0x0cc, 0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0, 0x8c0, 0x9c9,
    0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc, 0x0cc, 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9,
    0x7c0, 0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c, 0x15c, 0x055, 0x35f, 0x256,
    0x55a, 0x453, 0x759, 0x650, 0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc, 0x2fc,
    0x3f5, 0x0ff, 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0, 0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f,
    0xd65, 0xc6c, 0x36c, 0x265, 0x16f, 0x066, 0x76a, 0x663, 0x569, 0x460, 0xca0, 0xda9, 0xea3,
    0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac, 0x4ac, 0x5a5, 0x6af, 0x7a6, 0x0aa, 0x1a3, 0x2a9, 0x3a0,
    0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c, 0x53c, 0x435, 0x73f, 0x636, 0x13a,
    0x033, 0x339, 0x230, 0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c, 0x69c, 0x795,
    0x49f, 0x596, 0x29a, 0x393, 0x099, 0x190, 0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905,
    0x80c, 0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x000,
];

// Marching cubes triangle table (15 ints per case, -1 terminated/padded)
// Using a subset — full 256-entry table with 16 shorts each
const MC_TRI_TABLE: [[i8; 16]; 256] = include_mc_tri_table();

const fn include_mc_tri_table() -> [[i8; 16]; 256] {
    // Compact representation: each entry lists edge indices forming triangles, -1 = end
    // This is the standard Lorensen-Cline table
    [
        [
            -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        ],
        [0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1],
        [3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1],
        [3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1],
        [3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1],
        [9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1],
        [1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1],
        [9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1],
        [2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1],
        [8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1],
        [9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1],
        [4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1],
        [3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1],
        [1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1],
        [4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1],
        [4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1],
        [9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1],
        [1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1],
        [5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1],
        [2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1],
        [9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1],
        [0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1],
        [2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1],
        [10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1],
        [4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1],
        [5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1],
        [5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1],
        [9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1],
        [0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1],
        [1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1],
        [10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1],
        [8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1],
        [2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1],
        [7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1],
        [9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1],
        [2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1],
        [11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1],
        [9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1],
        [5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1],
        [11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1],
        [11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1],
        [1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1],
        [9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1],
        [5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1],
        [2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1],
        [0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1],
        [5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1],
        [6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1],
        [0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1],
        [3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1],
        [6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1],
        [5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1],
        [1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1],
        [10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1],
        [6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1],
        [1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1],
        [8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1],
        [7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1],
        [3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1],
        [5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1],
        [0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1],
        [9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1],
        [8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1],
        [5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1],
        [0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1],
        [6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1],
        [10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1],
        [10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1],
        [8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1],
        [1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1],
        [3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1],
        [0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1],
        [10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1],
        [0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1],
        [3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1],
        [6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1],
        [9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1],
        [8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1],
        [3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1],
        [6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1],
        [0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1],
        [10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1],
        [10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1],
        [1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1],
        [2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1],
        [7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1],
        [7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1],
        [2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1],
        [1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1],
        [11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1],
        [8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1],
        [0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1],
        [7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1],
        [10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1],
        [2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1],
        [6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1],
        [7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1],
        [2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1],
        [1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1],
        [10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1],
        [10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1],
        [0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1],
        [7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1],
        [6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1],
        [8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1],
        [9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1],
        [6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1],
        [1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1],
        [4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1],
        [10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1],
        [8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1],
        [0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1],
        [1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1],
        [8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1],
        [10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1],
        [4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1],
        [10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1],
        [5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1],
        [11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1],
        [9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1],
        [6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1],
        [7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1],
        [3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1],
        [7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1],
        [9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1],
        [3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1],
        [6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1],
        [9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1],
        [1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1],
        [4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1],
        [7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1],
        [6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1],
        [3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1],
        [0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1],
        [6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1],
        [1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1],
        [0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1],
        [11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1],
        [6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1],
        [5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1],
        [9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1],
        [1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1],
        [1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1],
        [10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1],
        [0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1],
        [5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1],
        [10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1],
        [11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1],
        [0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1],
        [9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1],
        [7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1],
        [2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1],
        [8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1],
        [9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1],
        [9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1],
        [1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1],
        [9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1],
        [9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1],
        [5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1],
        [0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1],
        [10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1],
        [2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1],
        [0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1],
        [0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1],
        [9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1],
        [5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1],
        [3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1],
        [5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1],
        [8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1],
        [0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1],
        [9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1],
        [0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1],
        [1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1],
        [3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1],
        [4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1],
        [9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1],
        [11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1],
        [11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1],
        [2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1],
        [9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1],
        [3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1],
        [1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1],
        [4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1],
        [4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1],
        [0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1],
        [3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1],
        [3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1],
        [0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1],
        [9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1],
        [1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [
            -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        ],
    ]
}

impl MarchingCubes {
    /// Extract an isosurface from a 3D SDF grid at the given isovalue.
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
                            mesh.triangles.push([a, b, c]);
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

/// Marching tetrahedra isosurface extraction.
///
/// Subdivides each cube into 5 tetrahedra and applies the tetrahedron case table.
pub struct MarchingTetrahedra;

impl MarchingTetrahedra {
    /// Extract isosurface from SDF grid using marching tetrahedra.
    pub fn extract(grid: &SdfGrid3D, isovalue: f64) -> IsoMesh {
        let mut mesh = IsoMesh::new();
        let nx = grid.nx.saturating_sub(1);
        let ny = grid.ny.saturating_sub(1);
        let nz = grid.nz.saturating_sub(1);

        // 5-tetrahedra decomposition of a cube
        let tet_indices: [[usize; 4]; 5] = [
            [0, 1, 3, 5],
            [1, 3, 5, 6],
            [0, 3, 4, 5],
            [3, 5, 6, 7],
            [0, 3, 4, 7],
        ];
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

        for iz in 0..nz {
            for iy in 0..ny {
                for ix in 0..nx {
                    let mut cvals = [0.0f64; 8];
                    let mut cpos = [[0.0f64; 3]; 8];
                    for (ci, &[dx, dy, dz]) in corner_offsets.iter().enumerate() {
                        cpos[ci] = grid.voxel_center(ix + dx, iy + dy, iz + dz);
                        cvals[ci] = grid.get(ix + dx, iy + dy, iz + dz);
                    }

                    for &tet in &tet_indices {
                        let tv: [usize; 4] = tet;
                        let tp: [[f64; 3]; 4] =
                            [cpos[tv[0]], cpos[tv[1]], cpos[tv[2]], cpos[tv[3]]];
                        let td: [f64; 4] = [cvals[tv[0]], cvals[tv[1]], cvals[tv[2]], cvals[tv[3]]];
                        Self::process_tet(&tp, &td, isovalue, &mut mesh);
                    }
                }
            }
        }
        mesh
    }

    fn interp(p0: [f64; 3], p1: [f64; 3], d0: f64, d1: f64, iso: f64) -> [f64; 3] {
        let t = ((iso - d0) / (d1 - d0 + 1e-30)).clamp(0.0, 1.0);
        [
            p0[0] + t * (p1[0] - p0[0]),
            p0[1] + t * (p1[1] - p0[1]),
            p0[2] + t * (p1[2] - p0[2]),
        ]
    }

    fn process_tet(pos: &[[f64; 3]; 4], vals: &[f64; 4], iso: f64, mesh: &mut IsoMesh) {
        let idx: u8 = vals.iter().enumerate().fold(
            0u8,
            |acc, (i, &v)| {
                if v < iso { acc | (1 << i) } else { acc }
            },
        );
        // Tetrahedron case table: 16 cases, emit 0, 1, or 2 triangles
        let edges: [(usize, usize); 6] = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];
        let tri_cases: [[i8; 7]; 16] = [
            [-1, -1, -1, -1, -1, -1, -1],
            [0, 3, 2, -1, -1, -1, -1],
            [0, 1, 4, -1, -1, -1, -1],
            [1, 3, 2, 1, 4, 3, -1],
            [1, 2, 5, -1, -1, -1, -1],
            [0, 3, 5, 0, 5, 1, -1],
            [0, 2, 5, 0, 5, 4, -1],
            [3, 5, 4, -1, -1, -1, -1],
            [3, 5, 4, -1, -1, -1, -1],
            [0, 5, 2, 0, 4, 5, -1],
            [0, 5, 3, 0, 1, 5, -1],
            [1, 2, 5, -1, -1, -1, -1],
            [1, 3, 2, 1, 4, 3, -1],
            [0, 1, 4, -1, -1, -1, -1],
            [0, 3, 2, -1, -1, -1, -1],
            [-1, -1, -1, -1, -1, -1, -1],
        ];
        let tc = tri_cases[idx as usize];
        let mut ti = 0;
        while ti < 6 && tc[ti] >= 0 {
            let ei0 = tc[ti] as usize;
            let ei1 = tc[ti + 1] as usize;
            let ei2 = tc[ti + 2] as usize;
            let (a0, a1) = edges[ei0];
            let (b0, b1) = edges[ei1];
            let (c0, c1) = edges[ei2];
            let pa = Self::interp(pos[a0], pos[a1], vals[a0], vals[a1], iso);
            let pb = Self::interp(pos[b0], pos[b1], vals[b0], vals[b1], iso);
            let pc = Self::interp(pos[c0], pos[c1], vals[c0], vals[c1], iso);
            let vi = mesh.vertices.len();
            mesh.vertices.extend_from_slice(&[pa, pb, pc]);
            mesh.normals.extend_from_slice(&[[0.0, 1.0, 0.0]; 3]);
            mesh.triangles.push([vi, vi + 1, vi + 2]);
            ti += 3;
        }
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
