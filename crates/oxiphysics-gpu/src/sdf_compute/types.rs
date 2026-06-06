//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)
use rayon::prelude::*;

/// Boolean CSG operations combining two [`SdfShape`] primitives.
#[derive(Debug, Clone)]
pub enum SdfCombine {
    /// Union: min(d_a, d_b).
    Union(SdfShape, SdfShape),
    /// Intersection: max(d_a, d_b).
    Intersection(SdfShape, SdfShape),
    /// Subtraction of b from a: max(d_a, -d_b).
    Subtraction(SdfShape, SdfShape),
    /// Smooth union with blending factor k.
    SmoothUnion(SdfShape, SdfShape, f64),
}
impl SdfCombine {
    /// Signed distance from point `p` to this combined shape.
    pub fn signed_distance(&self, p: [f64; 3]) -> f64 {
        match self {
            SdfCombine::Union(a, b) => a.signed_distance(p).min(b.signed_distance(p)),
            SdfCombine::Intersection(a, b) => a.signed_distance(p).max(b.signed_distance(p)),
            SdfCombine::Subtraction(a, b) => a.signed_distance(p).max(-b.signed_distance(p)),
            SdfCombine::SmoothUnion(a, b, k) => {
                let da = a.signed_distance(p);
                let db = b.signed_distance(p);
                let h = (0.5 + 0.5 * (db - da) / k).clamp(0.0, 1.0);
                db * (1.0 - h) + da * h - k * h * (1.0 - h)
            }
        }
    }
}
/// An analytic signed distance shape primitive.
#[derive(Debug, Clone)]
pub enum SdfShape {
    /// A sphere with a given centre and radius.
    Sphere {
        /// Centre of the sphere.
        center: [f64; 3],
        /// Radius.
        r: f64,
    },
    /// An axis-aligned box specified by centre and half-extents.
    Box3 {
        /// Centre of the box.
        center: [f64; 3],
        /// Half-extents in each axis.
        half: [f64; 3],
    },
    /// A capsule (cylinder with hemispherical caps) between two end-points.
    Capsule {
        /// First end-point.
        a: [f64; 3],
        /// Second end-point.
        b: [f64; 3],
        /// Radius of the capsule.
        r: f64,
    },
}
impl SdfShape {
    /// Signed distance from point `p` to this shape.
    ///
    /// Negative values are inside the shape.
    pub fn signed_distance(&self, p: [f64; 3]) -> f64 {
        match self {
            SdfShape::Sphere { center, r } => {
                let dx = p[0] - center[0];
                let dy = p[1] - center[1];
                let dz = p[2] - center[2];
                (dx * dx + dy * dy + dz * dz).sqrt() - r
            }
            SdfShape::Box3 { center, half } => {
                let qx = (p[0] - center[0]).abs() - half[0];
                let qy = (p[1] - center[1]).abs() - half[1];
                let qz = (p[2] - center[2]).abs() - half[2];
                let ext = (qx.max(0.0).powi(2) + qy.max(0.0).powi(2) + qz.max(0.0).powi(2)).sqrt();
                let interior = qx.max(qy).max(qz).min(0.0);
                ext + interior
            }
            SdfShape::Capsule { a, b, r } => {
                let ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
                let ap = [p[0] - a[0], p[1] - a[1], p[2] - a[2]];
                let ab_len2 = ab[0] * ab[0] + ab[1] * ab[1] + ab[2] * ab[2];
                let t = if ab_len2 < 1e-12 {
                    0.0
                } else {
                    ((ap[0] * ab[0] + ap[1] * ab[1] + ap[2] * ab[2]) / ab_len2).clamp(0.0, 1.0)
                };
                let closest = [a[0] + t * ab[0], a[1] + t * ab[1], a[2] + t * ab[2]];
                let dx = p[0] - closest[0];
                let dy = p[1] - closest[1];
                let dz = p[2] - closest[2];
                (dx * dx + dy * dy + dz * dz).sqrt() - r
            }
        }
    }
}
/// A 3-D SDF grid filled analytically via `generate_sdf_grid`.
///
/// Separate from the existing [`SdfGrid`] to keep naming unambiguous.
#[derive(Debug, Clone)]
pub struct GpuSdfGrid {
    /// Flat storage: index = `ix * ny * nz + iy * nz + iz`.
    pub data: Vec<f64>,
    /// Number of cells in x.
    pub nx: usize,
    /// Number of cells in y.
    pub ny: usize,
    /// Number of cells in z.
    pub nz: usize,
    /// World-space origin of the (0,0,0) corner.
    pub origin: [f64; 3],
    /// Uniform cell spacing.
    pub cell_size: f64,
}
impl GpuSdfGrid {
    /// Create a new grid filled with zeros.
    pub fn new(nx: usize, ny: usize, nz: usize, origin: [f64; 3], cell_size: f64) -> Self {
        Self {
            data: vec![0.0_f64; nx * ny * nz],
            nx,
            ny,
            nz,
            origin,
            cell_size,
        }
    }
    /// Flat index for cell `(ix, iy, iz)`.
    #[inline]
    pub fn index(&self, ix: usize, iy: usize, iz: usize) -> usize {
        ix * self.ny * self.nz + iy * self.nz + iz
    }
    /// Read the value at cell `(ix, iy, iz)`.
    #[inline]
    pub fn get(&self, ix: usize, iy: usize, iz: usize) -> f64 {
        self.data[self.index(ix, iy, iz)]
    }
    /// World-space centre of cell `(ix, iy, iz)`.
    #[inline]
    pub fn cell_center(&self, ix: usize, iy: usize, iz: usize) -> [f64; 3] {
        [
            self.origin[0] + (ix as f64 + 0.5) * self.cell_size,
            self.origin[1] + (iy as f64 + 0.5) * self.cell_size,
            self.origin[2] + (iz as f64 + 0.5) * self.cell_size,
        ]
    }
    /// Trilinear interpolation of the SDF at world-space point `p`.
    ///
    /// Points outside the grid are clamped to the nearest grid value.
    pub fn sample_trilinear(&self, p: [f64; 3]) -> f64 {
        let fx = (p[0] - self.origin[0]) / self.cell_size - 0.5;
        let fy = (p[1] - self.origin[1]) / self.cell_size - 0.5;
        let fz = (p[2] - self.origin[2]) / self.cell_size - 0.5;
        let ix = fx.floor().clamp(0.0, (self.nx - 1) as f64) as usize;
        let iy = fy.floor().clamp(0.0, (self.ny - 1) as f64) as usize;
        let iz = fz.floor().clamp(0.0, (self.nz - 1) as f64) as usize;
        let tx = (fx - ix as f64).clamp(0.0, 1.0);
        let ty = (fy - iy as f64).clamp(0.0, 1.0);
        let tz = (fz - iz as f64).clamp(0.0, 1.0);
        let nx1 = (ix + 1).min(self.nx - 1);
        let ny1 = (iy + 1).min(self.ny - 1);
        let nz1 = (iz + 1).min(self.nz - 1);
        let c000 = self.get(ix, iy, iz);
        let c100 = self.get(nx1, iy, iz);
        let c010 = self.get(ix, ny1, iz);
        let c110 = self.get(nx1, ny1, iz);
        let c001 = self.get(ix, iy, nz1);
        let c101 = self.get(nx1, iy, nz1);
        let c011 = self.get(ix, ny1, nz1);
        let c111 = self.get(nx1, ny1, nz1);
        let c00 = c000 * (1.0 - tx) + c100 * tx;
        let c10 = c010 * (1.0 - tx) + c110 * tx;
        let c01 = c001 * (1.0 - tx) + c101 * tx;
        let c11 = c011 * (1.0 - tx) + c111 * tx;
        let c0 = c00 * (1.0 - ty) + c10 * ty;
        let c1 = c01 * (1.0 - ty) + c11 * ty;
        c0 * (1.0 - tz) + c1 * tz
    }
    /// Estimate the SDF gradient at `p` using central finite differences.
    ///
    /// The step size is `cell_size * 0.5`.
    pub fn gradient_at(&self, p: [f64; 3]) -> [f64; 3] {
        let h = self.cell_size * 0.5;
        let gx = (self.sample_trilinear([p[0] + h, p[1], p[2]])
            - self.sample_trilinear([p[0] - h, p[1], p[2]]))
            / (2.0 * h);
        let gy = (self.sample_trilinear([p[0], p[1] + h, p[2]])
            - self.sample_trilinear([p[0], p[1] - h, p[2]]))
            / (2.0 * h);
        let gz = (self.sample_trilinear([p[0], p[1], p[2] + h])
            - self.sample_trilinear([p[0], p[1], p[2] - h]))
            / (2.0 * h);
        [gx, gy, gz]
    }
}
/// Sphere tracing (ray marching) result.
pub struct SphereTraceResult {
    /// Whether the ray hit the surface.
    pub hit: bool,
    /// Position of the hit point (if hit).
    pub position: [f64; 3],
    /// Distance traveled along the ray.
    pub t: f64,
    /// Number of iterations used.
    pub iterations: usize,
}
/// A 3-D signed distance field stored on a uniform Cartesian grid.
pub struct SdfGrid {
    /// Number of cells in the x direction.
    pub nx: usize,
    /// Number of cells in the y direction.
    pub ny: usize,
    /// Number of cells in the z direction.
    pub nz: usize,
    /// Uniform cell spacing (same in all directions).
    pub dx: f64,
    /// World-space coordinate of the (0,0,0) grid corner.
    pub origin: [f64; 3],
    /// Flat storage: index = `i*ny*nz + j*nz + k`.
    pub values: Vec<f64>,
}
impl SdfGrid {
    /// Create a new grid filled with `f64::MAX`.
    pub fn new(nx: usize, ny: usize, nz: usize, dx: f64, origin: [f64; 3]) -> Self {
        let n = nx * ny * nz;
        Self {
            nx,
            ny,
            nz,
            dx,
            origin,
            values: vec![f64::MAX; n],
        }
    }
    /// Flat index for cell `(i, j, k)`.
    #[inline]
    pub fn index(&self, i: usize, j: usize, k: usize) -> usize {
        i * self.ny * self.nz + j * self.nz + k
    }
    /// World-space centre of cell `(i, j, k)`.
    #[inline]
    pub fn world_pos(&self, i: usize, j: usize, k: usize) -> [f64; 3] {
        [
            self.origin[0] + (i as f64 + 0.5) * self.dx,
            self.origin[1] + (j as f64 + 0.5) * self.dx,
            self.origin[2] + (k as f64 + 0.5) * self.dx,
        ]
    }
    /// Read the SDF value at `(i, j, k)`.
    #[inline]
    pub fn get(&self, i: usize, j: usize, k: usize) -> f64 {
        self.values[self.index(i, j, k)]
    }
    /// Write the SDF value at `(i, j, k)`.
    #[inline]
    pub fn set(&mut self, i: usize, j: usize, k: usize, v: f64) {
        let idx = self.index(i, j, k);
        self.values[idx] = v;
    }
    /// Fill the grid with the signed distance to a sphere.
    pub fn compute_sphere_sdf(&mut self, center: [f64; 3], radius: f64) {
        let _nx = self.nx;
        let ny = self.ny;
        let nz = self.nz;
        let dx = self.dx;
        let origin = self.origin;
        self.values.par_iter_mut().enumerate().for_each(|(idx, v)| {
            let i = idx / (ny * nz);
            let j = (idx / nz) % ny;
            let k = idx % nz;
            let px = origin[0] + (i as f64 + 0.5) * dx;
            let py = origin[1] + (j as f64 + 0.5) * dx;
            let pz = origin[2] + (k as f64 + 0.5) * dx;
            let dist =
                ((px - center[0]).powi(2) + (py - center[1]).powi(2) + (pz - center[2]).powi(2))
                    .sqrt();
            *v = dist - radius;
        });
    }
    /// Fill the grid with the signed distance to an axis-aligned box.
    pub fn compute_box_sdf(&mut self, box_center: [f64; 3], half_extents: [f64; 3]) {
        let _nx = self.nx;
        let ny = self.ny;
        let nz = self.nz;
        let dx = self.dx;
        let origin = self.origin;
        self.values.par_iter_mut().enumerate().for_each(|(idx, v)| {
            let i = idx / (ny * nz);
            let j = (idx / nz) % ny;
            let k = idx % nz;
            let px = origin[0] + (i as f64 + 0.5) * dx - box_center[0];
            let py = origin[1] + (j as f64 + 0.5) * dx - box_center[1];
            let pz = origin[2] + (k as f64 + 0.5) * dx - box_center[2];
            let qx = px.abs() - half_extents[0];
            let qy = py.abs() - half_extents[1];
            let qz = pz.abs() - half_extents[2];
            let ext = (qx.max(0.0).powi(2) + qy.max(0.0).powi(2) + qz.max(0.0).powi(2)).sqrt();
            let interior = qx.max(qy).max(qz).min(0.0);
            *v = ext + interior;
        });
    }
    /// Fill the grid with the signed distance to an infinite cylinder
    /// aligned along the z-axis.
    pub fn compute_cylinder_sdf(&mut self, center: [f64; 2], radius: f64, half_height: f64) {
        let ny = self.ny;
        let nz = self.nz;
        let dx = self.dx;
        let origin = self.origin;
        self.values.par_iter_mut().enumerate().for_each(|(idx, v)| {
            let i = idx / (ny * nz);
            let j = (idx / nz) % ny;
            let k = idx % nz;
            let px = origin[0] + (i as f64 + 0.5) * dx - center[0];
            let py = origin[1] + (j as f64 + 0.5) * dx - center[1];
            let pz = origin[2] + (k as f64 + 0.5) * dx;
            let r = (px * px + py * py).sqrt();
            let d_radial = r - radius;
            let d_axial = pz.abs() - half_height;
            let ext = (d_radial.max(0.0).powi(2) + d_axial.max(0.0).powi(2)).sqrt();
            let interior = d_radial.max(d_axial).min(0.0);
            *v = ext + interior;
        });
    }
    /// Fill the grid with the signed distance to a torus
    /// centred at the origin in the xz-plane.
    pub fn compute_torus_sdf(&mut self, center: [f64; 3], major_radius: f64, minor_radius: f64) {
        let ny = self.ny;
        let nz = self.nz;
        let dx = self.dx;
        let origin = self.origin;
        self.values.par_iter_mut().enumerate().for_each(|(idx, v)| {
            let i = idx / (ny * nz);
            let j = (idx / nz) % ny;
            let k = idx % nz;
            let px = origin[0] + (i as f64 + 0.5) * dx - center[0];
            let py = origin[1] + (j as f64 + 0.5) * dx - center[1];
            let pz = origin[2] + (k as f64 + 0.5) * dx - center[2];
            let q_x = (px * px + pz * pz).sqrt() - major_radius;
            let dist = (q_x * q_x + py * py).sqrt() - minor_radius;
            *v = dist;
        });
    }
    /// Estimate the gradient of the SDF at `(i, j, k)` using central differences.
    pub fn gradient_at(&self, i: usize, j: usize, k: usize) -> [f64; 3] {
        let two_dx = 2.0 * self.dx;
        let gx = if i == 0 {
            (self.get(i + 1, j, k) - self.get(i, j, k)) / self.dx
        } else if i + 1 == self.nx {
            (self.get(i, j, k) - self.get(i - 1, j, k)) / self.dx
        } else {
            (self.get(i + 1, j, k) - self.get(i - 1, j, k)) / two_dx
        };
        let gy = if j == 0 {
            (self.get(i, j + 1, k) - self.get(i, j, k)) / self.dx
        } else if j + 1 == self.ny {
            (self.get(i, j, k) - self.get(i, j - 1, k)) / self.dx
        } else {
            (self.get(i, j + 1, k) - self.get(i, j - 1, k)) / two_dx
        };
        let gz = if k == 0 {
            (self.get(i, j, k + 1) - self.get(i, j, k)) / self.dx
        } else if k + 1 == self.nz {
            (self.get(i, j, k) - self.get(i, j, k - 1)) / self.dx
        } else {
            (self.get(i, j, k + 1) - self.get(i, j, k - 1)) / two_dx
        };
        [gx, gy, gz]
    }
    /// Total number of cells in the grid.
    #[inline]
    pub fn total_cells(&self) -> usize {
        self.nx * self.ny * self.nz
    }
    /// Trilinear interpolation of the SDF at an arbitrary world-space point.
    ///
    /// Returns `None` if the point is outside the grid.
    pub fn sample(&self, pos: [f64; 3]) -> Option<f64> {
        let fx = (pos[0] - self.origin[0]) / self.dx - 0.5;
        let fy = (pos[1] - self.origin[1]) / self.dx - 0.5;
        let fz = (pos[2] - self.origin[2]) / self.dx - 0.5;
        if fx < 0.0 || fy < 0.0 || fz < 0.0 {
            return None;
        }
        let ix = fx as usize;
        let iy = fy as usize;
        let iz = fz as usize;
        if ix + 1 >= self.nx || iy + 1 >= self.ny || iz + 1 >= self.nz {
            return None;
        }
        let tx = fx - ix as f64;
        let ty = fy - iy as f64;
        let tz = fz - iz as f64;
        let c000 = self.get(ix, iy, iz);
        let c100 = self.get(ix + 1, iy, iz);
        let c010 = self.get(ix, iy + 1, iz);
        let c110 = self.get(ix + 1, iy + 1, iz);
        let c001 = self.get(ix, iy, iz + 1);
        let c101 = self.get(ix + 1, iy, iz + 1);
        let c011 = self.get(ix, iy + 1, iz + 1);
        let c111 = self.get(ix + 1, iy + 1, iz + 1);
        let c00 = c000 * (1.0 - tx) + c100 * tx;
        let c10 = c010 * (1.0 - tx) + c110 * tx;
        let c01 = c001 * (1.0 - tx) + c101 * tx;
        let c11 = c011 * (1.0 - tx) + c111 * tx;
        let c0 = c00 * (1.0 - ty) + c10 * ty;
        let c1 = c01 * (1.0 - ty) + c11 * ty;
        Some(c0 * (1.0 - tz) + c1 * tz)
    }
    /// Compute the gradient at an arbitrary point using trilinear interpolation
    /// of gradient values.
    pub fn gradient_at_point(&self, pos: [f64; 3]) -> Option<[f64; 3]> {
        let fx = (pos[0] - self.origin[0]) / self.dx - 0.5;
        let fy = (pos[1] - self.origin[1]) / self.dx - 0.5;
        let fz = (pos[2] - self.origin[2]) / self.dx - 0.5;
        if fx < 0.0 || fy < 0.0 || fz < 0.0 {
            return None;
        }
        let ix = fx as usize;
        let iy = fy as usize;
        let iz = fz as usize;
        if ix + 1 >= self.nx || iy + 1 >= self.ny || iz + 1 >= self.nz {
            return None;
        }
        let eps = self.dx * 0.5;
        let gx = (self.sample([pos[0] + eps, pos[1], pos[2]]).unwrap_or(0.0)
            - self.sample([pos[0] - eps, pos[1], pos[2]]).unwrap_or(0.0))
            / (2.0 * eps);
        let gy = (self.sample([pos[0], pos[1] + eps, pos[2]]).unwrap_or(0.0)
            - self.sample([pos[0], pos[1] - eps, pos[2]]).unwrap_or(0.0))
            / (2.0 * eps);
        let gz = (self.sample([pos[0], pos[1], pos[2] + eps]).unwrap_or(0.0)
            - self.sample([pos[0], pos[1], pos[2] - eps]).unwrap_or(0.0))
            / (2.0 * eps);
        Some([gx, gy, gz])
    }
}
/// A triangle in 3-D space.
#[derive(Debug, Clone)]
pub struct Triangle {
    /// Three vertex positions.
    pub v: [[f64; 3]; 3],
}
/// Query result for closest point / normal / distance.
#[derive(Debug, Clone)]
pub struct DistanceQuery {
    /// Signed distance at the query point.
    pub distance: f64,
    /// Normal direction at the query point (normalised gradient).
    pub normal: [f64; 3],
    /// Closest point on the surface (approximate, from gradient ray march).
    pub closest_point: [f64; 3],
    /// Whether the query point is inside the surface.
    pub is_inside: bool,
}
