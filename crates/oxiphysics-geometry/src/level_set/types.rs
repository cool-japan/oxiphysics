//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
use super::functions::{MC_EDGE_CORNERS, MC_EDGE_TABLE, MC_TRI_TABLE};

/// Signed distance transform from a point cloud or triangle mesh.
///
/// Provides brute-force SDF computation (exact for small point sets) and a
/// fast alternating-sweep approximation for larger inputs.
pub struct SignedDistanceTransform {
    /// Source points used to define the zero level set.
    pub source_points: Vec<[f64; 3]>,
}
impl SignedDistanceTransform {
    /// Create a new [`SignedDistanceTransform`] from a set of source points.
    ///
    /// # Arguments
    /// * `points` – positions on or near the zero level set
    pub fn new(points: Vec<[f64; 3]>) -> Self {
        Self {
            source_points: points,
        }
    }
    /// Brute-force signed distance from world point `p` to the point cloud.
    ///
    /// Returns the distance to the nearest source point (always non-negative;
    /// sign must be determined by additional topology information).
    pub fn unsigned_distance(&self, p: [f64; 3]) -> f64 {
        self.source_points
            .iter()
            .map(|&q| {
                ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2)).sqrt()
            })
            .fold(f64::MAX, f64::min)
    }
    /// Build an unsigned distance field on a regular grid.
    ///
    /// For each grid node, computes the nearest point distance using brute force.
    ///
    /// # Arguments
    /// * `nx`, `ny`, `nz` – grid dimensions
    /// * `dx`             – grid spacing
    /// * `origin`         – grid origin
    pub fn compute_grid_unsigned(
        &self,
        nx: usize,
        ny: usize,
        nz: usize,
        dx: f64,
        origin: [f64; 3],
    ) -> LevelSetField {
        let mut field = LevelSetField::new(nx, ny, nz, dx, origin);
        for iz in 0..nz {
            for iy in 0..ny {
                for ix in 0..nx {
                    let p = field.node_pos(ix, iy, iz);
                    let d = self.unsigned_distance(p);
                    field.set(ix, iy, iz, d);
                }
            }
        }
        field
    }
    /// Compute a signed distance field on a regular grid.
    ///
    /// The sign is determined by the winding number test relative to the
    /// first source point (interior if distance < median distance).
    /// For a production use case this should use a proper inside/outside test.
    pub fn compute_grid_signed(
        &self,
        nx: usize,
        ny: usize,
        nz: usize,
        dx: f64,
        origin: [f64; 3],
    ) -> LevelSetField {
        let mut field = self.compute_grid_unsigned(nx, ny, nz, dx, origin);
        if self.source_points.is_empty() {
            return field;
        }
        let centroid: [f64; 3] = {
            let n = self.source_points.len() as f64;
            let s = self.source_points.iter().fold([0.0_f64; 3], |acc, p| {
                [acc[0] + p[0], acc[1] + p[1], acc[2] + p[2]]
            });
            [s[0] / n, s[1] / n, s[2] / n]
        };
        let radius = self.unsigned_distance(centroid);
        for iz in 0..nz {
            for iy in 0..ny {
                for ix in 0..nx {
                    let p = field.node_pos(ix, iy, iz);
                    let dc = ((p[0] - centroid[0]).powi(2)
                        + (p[1] - centroid[1]).powi(2)
                        + (p[2] - centroid[2]).powi(2))
                    .sqrt();
                    let val = field.get(ix, iy, iz);
                    if dc < radius {
                        field.set(ix, iy, iz, -val);
                    }
                }
            }
        }
        field
    }
    /// Distance from point `p` to a line segment \[a, b\].
    pub fn point_to_segment(p: [f64; 3], a: [f64; 3], b: [f64; 3]) -> f64 {
        let ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        let ap = [p[0] - a[0], p[1] - a[1], p[2] - a[2]];
        let len2 = ab[0] * ab[0] + ab[1] * ab[1] + ab[2] * ab[2];
        if len2 < 1e-30 {
            return (ap[0] * ap[0] + ap[1] * ap[1] + ap[2] * ap[2]).sqrt();
        }
        let t = ((ap[0] * ab[0] + ap[1] * ab[1] + ap[2] * ab[2]) / len2).clamp(0.0, 1.0);
        let closest = [a[0] + t * ab[0], a[1] + t * ab[1], a[2] + t * ab[2]];
        ((p[0] - closest[0]).powi(2) + (p[1] - closest[1]).powi(2) + (p[2] - closest[2]).powi(2))
            .sqrt()
    }
}
/// Level-set evolution by advection and mean curvature flow.
///
/// Implements the upwind finite difference scheme for advection and the
/// explicit mean-curvature-flow equation for smoothing.
pub struct LevelSetEvolution {
    /// Time step for explicit integration.
    pub dt: f64,
    /// Number of time steps per evolution call.
    pub n_steps: usize,
}
impl LevelSetEvolution {
    /// Create a new [`LevelSetEvolution`] solver.
    ///
    /// # Arguments
    /// * `dt`      – explicit time step (should satisfy CFL: dt ≤ dx/|v|_max)
    /// * `n_steps` – number of time steps per `advect` or `mean_curvature_flow` call
    pub fn new(dt: f64, n_steps: usize) -> Self {
        Self { dt, n_steps }
    }
    /// Advect the level-set field by a uniform velocity field using the upwind scheme.
    ///
    /// Solves ∂φ/∂t + v · ∇φ = 0 with first-order upwind spatial discretisation.
    ///
    /// # Arguments
    /// * `field`    – level-set field (modified in place)
    /// * `velocity` – uniform velocity vector \[vx, vy, vz\]
    pub fn advect(&self, field: &mut LevelSetField, velocity: [f64; 3]) {
        let dx = field.dx;
        let nx = field.nx;
        let ny = field.ny;
        let nz = field.nz;
        for _ in 0..self.n_steps {
            let old = field.data.clone();
            for iz in 0..nz {
                for iy in 0..ny {
                    for ix in 0..nx {
                        let phi_c = old[field.idx(ix, iy, iz)];
                        let dphi_x = if velocity[0] > 0.0 {
                            if ix > 0 {
                                (phi_c - old[field.idx(ix - 1, iy, iz)]) / dx
                            } else {
                                0.0
                            }
                        } else {
                            if ix + 1 < nx {
                                (old[field.idx(ix + 1, iy, iz)] - phi_c) / dx
                            } else {
                                0.0
                            }
                        };
                        let dphi_y = if velocity[1] > 0.0 {
                            if iy > 0 {
                                (phi_c - old[field.idx(ix, iy - 1, iz)]) / dx
                            } else {
                                0.0
                            }
                        } else {
                            if iy + 1 < ny {
                                (old[field.idx(ix, iy + 1, iz)] - phi_c) / dx
                            } else {
                                0.0
                            }
                        };
                        let dphi_z = if velocity[2] > 0.0 {
                            if iz > 0 {
                                (phi_c - old[field.idx(ix, iy, iz - 1)]) / dx
                            } else {
                                0.0
                            }
                        } else {
                            if iz + 1 < nz {
                                (old[field.idx(ix, iy, iz + 1)] - phi_c) / dx
                            } else {
                                0.0
                            }
                        };
                        let rhs =
                            velocity[0] * dphi_x + velocity[1] * dphi_y + velocity[2] * dphi_z;
                        let idx = field.idx(ix, iy, iz);
                        field.data[idx] = phi_c - self.dt * rhs;
                    }
                }
            }
        }
    }
    /// Evolve the level-set by mean curvature flow.
    ///
    /// Solves ∂φ/∂t = κ |∇φ| where κ is the mean curvature.
    /// Uses explicit Euler integration.
    ///
    /// # Arguments
    /// * `field` – level-set field (modified in place)
    pub fn mean_curvature_flow(&self, field: &mut LevelSetField) {
        let dx = field.dx;
        let nx = field.nx;
        let ny = field.ny;
        let nz = field.nz;
        for _ in 0..self.n_steps {
            let old = field.data.clone();
            let old_field = LevelSetField {
                nx,
                ny,
                nz,
                dx,
                origin: field.origin,
                data: old,
            };
            for iz in 0..nz {
                for iy in 0..ny {
                    for ix in 0..nx {
                        let kappa = old_field.mean_curvature(ix, iy, iz);
                        let grad = old_field.gradient(ix, iy, iz);
                        let grad_mag =
                            (grad[0] * grad[0] + grad[1] * grad[1] + grad[2] * grad[2]).sqrt();
                        let idx = field.idx(ix, iy, iz);
                        field.data[idx] = old_field.data[idx] + self.dt * kappa * grad_mag;
                    }
                }
            }
        }
    }
    /// Normal velocity extension: propagate the interface normal velocity outward.
    ///
    /// For each node computes the update φ_new = φ - dt * F |∇φ| where F is
    /// a prescribed scalar normal velocity.
    ///
    /// # Arguments
    /// * `field`            – level-set field (modified in place)
    /// * `normal_velocity`  – scalar normal speed (positive = expansion)
    pub fn normal_velocity_extension(&self, field: &mut LevelSetField, normal_velocity: f64) {
        let nx = field.nx;
        let ny = field.ny;
        let nz = field.nz;
        for _ in 0..self.n_steps {
            let old = field.data.clone();
            let old_field = LevelSetField {
                nx,
                ny,
                nz,
                dx: field.dx,
                origin: field.origin,
                data: old,
            };
            for iz in 0..nz {
                for iy in 0..ny {
                    for ix in 0..nx {
                        let grad = old_field.gradient(ix, iy, iz);
                        let grad_mag =
                            (grad[0] * grad[0] + grad[1] * grad[1] + grad[2] * grad[2]).sqrt();
                        let idx = field.idx(ix, iy, iz);
                        field.data[idx] =
                            old_field.data[idx] - self.dt * normal_velocity * grad_mag;
                    }
                }
            }
        }
    }
}
/// Boolean operations on level-set fields via min/max composition.
///
/// All operations create a new [`LevelSetField`] with the same grid layout as
/// `field_a`.  The fields must have identical grid dimensions, spacing, and origin.
pub struct LevelSetBoolean;
impl LevelSetBoolean {
    /// Union of two level-set fields: φ = min(φ_a, φ_b).
    ///
    /// Corresponds to the union of the two implicit surfaces.
    ///
    /// # Arguments
    /// * `field_a` – first operand
    /// * `field_b` – second operand (must have same grid as `field_a`)
    pub fn union(field_a: &LevelSetField, field_b: &LevelSetField) -> LevelSetField {
        let mut result = LevelSetField::new(
            field_a.nx,
            field_a.ny,
            field_a.nz,
            field_a.dx,
            field_a.origin,
        );
        for i in 0..result.data.len() {
            let va = field_a.data.get(i).copied().unwrap_or(f64::MAX);
            let vb = field_b.data.get(i).copied().unwrap_or(f64::MAX);
            result.data[i] = va.min(vb);
        }
        result
    }
    /// Intersection of two level-set fields: φ = max(φ_a, φ_b).
    ///
    /// Corresponds to the intersection of the two implicit surfaces.
    ///
    /// # Arguments
    /// * `field_a` – first operand
    /// * `field_b` – second operand (must have same grid as `field_a`)
    pub fn intersection(field_a: &LevelSetField, field_b: &LevelSetField) -> LevelSetField {
        let mut result = LevelSetField::new(
            field_a.nx,
            field_a.ny,
            field_a.nz,
            field_a.dx,
            field_a.origin,
        );
        for i in 0..result.data.len() {
            let va = field_a.data.get(i).copied().unwrap_or(f64::MAX);
            let vb = field_b.data.get(i).copied().unwrap_or(f64::MAX);
            result.data[i] = va.max(vb);
        }
        result
    }
    /// Difference A − B: φ = max(φ_a, −φ_b).
    ///
    /// Corresponds to removing the region of `field_b` from `field_a`.
    ///
    /// # Arguments
    /// * `field_a` – base operand
    /// * `field_b` – region to subtract
    pub fn difference(field_a: &LevelSetField, field_b: &LevelSetField) -> LevelSetField {
        let mut result = LevelSetField::new(
            field_a.nx,
            field_a.ny,
            field_a.nz,
            field_a.dx,
            field_a.origin,
        );
        for i in 0..result.data.len() {
            let va = field_a.data.get(i).copied().unwrap_or(f64::MAX);
            let vb = field_b.data.get(i).copied().unwrap_or(f64::MAX);
            result.data[i] = va.max(-vb);
        }
        result
    }
    /// Complement (negation): φ = −φ_a.
    ///
    /// Flips the inside/outside classification of the level-set.
    ///
    /// # Arguments
    /// * `field_a` – the level-set field to negate
    pub fn complement(field_a: &LevelSetField) -> LevelSetField {
        let mut result = LevelSetField::new(
            field_a.nx,
            field_a.ny,
            field_a.nz,
            field_a.dx,
            field_a.origin,
        );
        for i in 0..result.data.len() {
            result.data[i] = -field_a.data[i];
        }
        result
    }
    /// Smooth minimum blend (C¹ union): φ = smin(φ_a, φ_b, k).
    ///
    /// Produces a smooth blend at the union boundary with blend radius `k`.
    ///
    /// # Arguments
    /// * `field_a` – first operand
    /// * `field_b` – second operand
    /// * `k`       – blend radius (larger k → smoother blend)
    pub fn smooth_union(field_a: &LevelSetField, field_b: &LevelSetField, k: f64) -> LevelSetField {
        let mut result = LevelSetField::new(
            field_a.nx,
            field_a.ny,
            field_a.nz,
            field_a.dx,
            field_a.origin,
        );
        for i in 0..result.data.len() {
            let va = field_a.data.get(i).copied().unwrap_or(f64::MAX);
            let vb = field_b.data.get(i).copied().unwrap_or(f64::MAX);
            let h = (0.5 + 0.5 * (vb - va) / k.max(1e-30)).clamp(0.0, 1.0);
            result.data[i] = va * h + vb * (1.0 - h) - k * h * (1.0 - h);
        }
        result
    }
}
/// Signed distance function sampled on a uniform 3-D Cartesian grid.
///
/// The grid has `nx × ny × nz` cells.  Values are stored in row-major order
/// as `data[iz * ny * nx + iy * nx + ix]`.  Positive values are outside the
/// surface; negative values are inside.
pub struct LevelSetField {
    /// Number of grid points along x.
    pub nx: usize,
    /// Number of grid points along y.
    pub ny: usize,
    /// Number of grid points along z.
    pub nz: usize,
    /// Grid spacing (uniform in all directions).
    pub dx: f64,
    /// Origin of the grid (coordinates of node \[0,0,0\]).
    pub origin: [f64; 3],
    /// Flat array of signed distance values, length `nx * ny * nz`.
    pub data: Vec<f64>,
}
impl LevelSetField {
    /// Create a new [`LevelSetField`] filled with `f64::MAX` (uninitialized).
    ///
    /// # Arguments
    /// * `nx`, `ny`, `nz` – grid dimensions
    /// * `dx`             – uniform grid spacing
    /// * `origin`         – world-space coordinates of node (0,0,0)
    pub fn new(nx: usize, ny: usize, nz: usize, dx: f64, origin: [f64; 3]) -> Self {
        let data = vec![f64::MAX; nx * ny * nz];
        Self {
            nx,
            ny,
            nz,
            dx,
            origin,
            data,
        }
    }
    /// Linear index for grid node (ix, iy, iz).
    #[inline]
    pub fn idx(&self, ix: usize, iy: usize, iz: usize) -> usize {
        iz * self.ny * self.nx + iy * self.nx + ix
    }
    /// World-space coordinates of grid node (ix, iy, iz).
    #[inline]
    pub fn node_pos(&self, ix: usize, iy: usize, iz: usize) -> [f64; 3] {
        [
            self.origin[0] + ix as f64 * self.dx,
            self.origin[1] + iy as f64 * self.dx,
            self.origin[2] + iz as f64 * self.dx,
        ]
    }
    /// Get the level-set value at grid node (ix, iy, iz).
    ///
    /// Returns `f64::MAX` for out-of-bounds nodes.
    pub fn get(&self, ix: usize, iy: usize, iz: usize) -> f64 {
        if ix >= self.nx || iy >= self.ny || iz >= self.nz {
            return f64::MAX;
        }
        self.data[self.idx(ix, iy, iz)]
    }
    /// Set the level-set value at grid node (ix, iy, iz).
    ///
    /// Silently ignores out-of-bounds writes.
    pub fn set(&mut self, ix: usize, iy: usize, iz: usize, value: f64) {
        if ix >= self.nx || iy >= self.ny || iz >= self.nz {
            return;
        }
        let i = self.idx(ix, iy, iz);
        self.data[i] = value;
    }
    /// Tri-linear interpolation of the level-set field at world position `p`.
    ///
    /// Returns `f64::MAX` for positions outside the grid.
    pub fn interpolate(&self, p: [f64; 3]) -> f64 {
        let fx = (p[0] - self.origin[0]) / self.dx;
        let fy = (p[1] - self.origin[1]) / self.dx;
        let fz = (p[2] - self.origin[2]) / self.dx;
        if fx < 0.0 || fy < 0.0 || fz < 0.0 {
            return f64::MAX;
        }
        let ix = fx as usize;
        let iy = fy as usize;
        let iz = fz as usize;
        if ix + 1 >= self.nx || iy + 1 >= self.ny || iz + 1 >= self.nz {
            return f64::MAX;
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
        c0 * (1.0 - tz) + c1 * tz
    }
    /// Numerical gradient of the level-set field at grid node (ix, iy, iz).
    ///
    /// Uses central differences where possible, one-sided differences at boundaries.
    pub fn gradient(&self, ix: usize, iy: usize, iz: usize) -> [f64; 3] {
        let gx = if ix == 0 {
            (self.get(1, iy, iz) - self.get(0, iy, iz)) / self.dx
        } else if ix + 1 == self.nx {
            (self.get(ix, iy, iz) - self.get(ix - 1, iy, iz)) / self.dx
        } else {
            (self.get(ix + 1, iy, iz) - self.get(ix - 1, iy, iz)) / (2.0 * self.dx)
        };
        let gy = if iy == 0 {
            (self.get(ix, 1, iz) - self.get(ix, 0, iz)) / self.dx
        } else if iy + 1 == self.ny {
            (self.get(ix, iy, iz) - self.get(ix, iy - 1, iz)) / self.dx
        } else {
            (self.get(ix, iy + 1, iz) - self.get(ix, iy - 1, iz)) / (2.0 * self.dx)
        };
        let gz = if iz == 0 {
            (self.get(ix, iy, 1) - self.get(ix, iy, 0)) / self.dx
        } else if iz + 1 == self.nz {
            (self.get(ix, iy, iz) - self.get(ix, iy, iz - 1)) / self.dx
        } else {
            (self.get(ix, iy, iz + 1) - self.get(ix, iy, iz - 1)) / (2.0 * self.dx)
        };
        [gx, gy, gz]
    }
    /// Mean curvature κ = div(∇φ / |∇φ|) at grid node (ix, iy, iz).
    ///
    /// Computed via second-order finite differences.
    pub fn mean_curvature(&self, ix: usize, iy: usize, iz: usize) -> f64 {
        let phi = self.get(ix, iy, iz);
        let dx2 = self.dx * self.dx;
        let get = |x: usize, y: usize, z: usize| -> f64 {
            let xc = x.min(self.nx - 1);
            let yc = y.min(self.ny - 1);
            let zc = z.min(self.nz - 1);
            self.get(xc, yc, zc)
        };
        let pxp = if ix + 1 < self.nx {
            get(ix + 1, iy, iz)
        } else {
            phi
        };
        let pxm = if ix > 0 { get(ix - 1, iy, iz) } else { phi };
        let pyp = if iy + 1 < self.ny {
            get(ix, iy + 1, iz)
        } else {
            phi
        };
        let pym = if iy > 0 { get(ix, iy - 1, iz) } else { phi };
        let pzp = if iz + 1 < self.nz {
            get(ix, iy, iz + 1)
        } else {
            phi
        };
        let pzm = if iz > 0 { get(ix, iy, iz - 1) } else { phi };
        let phix = (pxp - pxm) / (2.0 * self.dx);
        let phiy = (pyp - pym) / (2.0 * self.dx);
        let phiz = (pzp - pzm) / (2.0 * self.dx);
        let phixx = (pxp - 2.0 * phi + pxm) / dx2;
        let phiyy = (pyp - 2.0 * phi + pym) / dx2;
        let phizz = (pzp - 2.0 * phi + pzm) / dx2;
        let grad2 = phix * phix + phiy * phiy + phiz * phiz;
        let grad_mag = grad2.sqrt().max(1e-30);

        (phixx + phiyy + phizz) / grad_mag
            - (phix * phix * phixx
                + phiy * phiy * phiyy
                + phiz * phiz * phizz
                + 2.0
                    * (phix * phiy * ((pxp + pyp - pxm - pym) / (4.0 * self.dx * self.dx) - 0.0)
                        + phix * phiz * 0.0
                        + phiy * phiz * 0.0))
                / (grad_mag * grad2).max(1e-30)
    }
    /// Initialize the level-set as a sphere SDF.
    ///
    /// Sets φ(x) = |x − center| − radius for all nodes.
    pub fn init_sphere(&mut self, center: [f64; 3], radius: f64) {
        for iz in 0..self.nz {
            for iy in 0..self.ny {
                for ix in 0..self.nx {
                    let p = self.node_pos(ix, iy, iz);
                    let d = ((p[0] - center[0]).powi(2)
                        + (p[1] - center[1]).powi(2)
                        + (p[2] - center[2]).powi(2))
                    .sqrt();
                    self.set(ix, iy, iz, d - radius);
                }
            }
        }
    }
    /// Initialize the level-set as an axis-aligned box SDF.
    ///
    /// Sets φ(x) to the signed distance to the box \[lo, hi\].
    pub fn init_box(&mut self, lo: [f64; 3], hi: [f64; 3]) {
        for iz in 0..self.nz {
            for iy in 0..self.ny {
                for ix in 0..self.nx {
                    let p = self.node_pos(ix, iy, iz);
                    let dx = (lo[0] - p[0]).max(p[0] - hi[0]).max(0.0);
                    let dy = (lo[1] - p[1]).max(p[1] - hi[1]).max(0.0);
                    let dz = (lo[2] - p[2]).max(p[2] - hi[2]).max(0.0);
                    let outside = (dx * dx + dy * dy + dz * dz).sqrt();
                    let inside = (p[0] - lo[0])
                        .min(hi[0] - p[0])
                        .min(p[1] - lo[1])
                        .min(hi[1] - p[1])
                        .min(p[2] - lo[2])
                        .min(hi[2] - p[2]);
                    let sdf = if inside < 0.0 { outside } else { -inside };
                    self.set(ix, iy, iz, sdf);
                }
            }
        }
    }
    /// Count nodes with φ ≤ 0 (inside or on the interface).
    pub fn count_inside(&self) -> usize {
        self.data.iter().filter(|&&v| v <= 0.0).count()
    }
    /// Minimum value of φ over all nodes.
    pub fn min_value(&self) -> f64 {
        self.data.iter().cloned().fold(f64::MAX, f64::min)
    }
    /// Maximum value of φ over all nodes.
    pub fn max_value(&self) -> f64 {
        self.data.iter().cloned().fold(f64::MIN, f64::max)
    }
}
/// Level-set reinitialization via fast marching method (redistancing).
///
/// Converts an arbitrary level-set function into a proper signed distance
/// function while preserving the zero-level-set location.
pub struct LevelSetReinit {
    /// Convergence tolerance for iterative reinitialization.
    pub tol: f64,
    /// Maximum number of pseudo-time iterations.
    pub max_iter: usize,
}
impl LevelSetReinit {
    /// Create a new reinitialization solver with given tolerance and iteration limit.
    ///
    /// # Arguments
    /// * `tol`      – convergence tolerance (∞-norm of update)
    /// * `max_iter` – maximum pseudo-time steps
    pub fn new(tol: f64, max_iter: usize) -> Self {
        Self { tol, max_iter }
    }
    /// Reinitialize a level-set field using the PDE-based Sussman scheme.
    ///
    /// Integrates ∂φ/∂τ = sign(φ₀)(1 − |∇φ|) until steady state.
    /// Modifies `field.data` in place.
    ///
    /// # Arguments
    /// * `field` – the level-set field to reinitialize
    pub fn reinitialize(&self, field: &mut LevelSetField) {
        let dx = field.dx;
        let dt = 0.5 * dx;
        let nx = field.nx;
        let ny = field.ny;
        let nz = field.nz;
        let phi0: Vec<f64> = field.data.clone();
        let sign0: Vec<f64> = phi0
            .iter()
            .map(|&v| {
                if v > 0.0 {
                    1.0
                } else if v < 0.0 {
                    -1.0
                } else {
                    0.0
                }
            })
            .collect();
        for _iter in 0..self.max_iter {
            let mut new_data = field.data.clone();
            let mut max_change = 0.0_f64;
            for iz in 0..nz {
                for iy in 0..ny {
                    for ix in 0..nx {
                        let idx = field.idx(ix, iy, iz);
                        let s = sign0[idx];
                        if s == 0.0 {
                            continue;
                        }
                        let phi_c = field.data[idx];
                        let grad_mag = godunov_grad_mag(field, ix, iy, iz, s);
                        let update = s * (1.0 - grad_mag);
                        new_data[idx] = phi_c + dt * update;
                        max_change = max_change.max((dt * update).abs());
                    }
                }
            }
            field.data = new_data;
            if max_change < self.tol {
                break;
            }
        }
    }
    /// Simple iterative fast marching sweep reinitialization (alternating directions).
    ///
    /// Faster than the PDE approach for well-behaved fields; performs `n_sweeps`
    /// alternating-direction passes of the first-order update.
    ///
    /// # Arguments
    /// * `field`   – the level-set field to reinitialize
    /// * `n_sweeps`– number of alternating-direction sweep passes
    pub fn fast_sweep(&self, field: &mut LevelSetField, n_sweeps: usize) {
        let dx = field.dx;
        let nx = field.nx;
        let ny = field.ny;
        let nz = field.nz;
        for sweep in 0..n_sweeps {
            let ix_range: Vec<usize> = if sweep % 2 == 0 {
                (0..nx).collect()
            } else {
                (0..nx).rev().collect()
            };
            let iy_range: Vec<usize> = if (sweep / 2) % 2 == 0 {
                (0..ny).collect()
            } else {
                (0..ny).rev().collect()
            };
            let iz_range: Vec<usize> = if (sweep / 4) % 2 == 0 {
                (0..nz).collect()
            } else {
                (0..nz).rev().collect()
            };
            for &iz in &iz_range {
                for &iy in &iy_range {
                    for &ix in &ix_range {
                        let phi_c = field.get(ix, iy, iz);
                        let axm = if ix > 0 {
                            field.get(ix - 1, iy, iz)
                        } else {
                            phi_c
                        };
                        let axp = if ix + 1 < nx {
                            field.get(ix + 1, iy, iz)
                        } else {
                            phi_c
                        };
                        let aym = if iy > 0 {
                            field.get(ix, iy - 1, iz)
                        } else {
                            phi_c
                        };
                        let ayp = if iy + 1 < ny {
                            field.get(ix, iy + 1, iz)
                        } else {
                            phi_c
                        };
                        let azm = if iz > 0 {
                            field.get(ix, iy, iz - 1)
                        } else {
                            phi_c
                        };
                        let azp = if iz + 1 < nz {
                            field.get(ix, iy, iz + 1)
                        } else {
                            phi_c
                        };
                        let ax = if phi_c >= 0.0 {
                            axm.min(axp)
                        } else {
                            axm.max(axp)
                        };
                        let ay = if phi_c >= 0.0 {
                            aym.min(ayp)
                        } else {
                            aym.max(ayp)
                        };
                        let az = if phi_c >= 0.0 {
                            azm.min(azp)
                        } else {
                            azm.max(azp)
                        };
                        let new_val = if phi_c >= 0.0 {
                            (ax.abs().min(ay.abs()).min(az.abs()) + dx).min(phi_c.abs())
                        } else {
                            let d = (ax.abs().min(ay.abs()).min(az.abs()) + dx).min(phi_c.abs());
                            -d
                        };
                        let _ = ax;
                        let _ = ay;
                        let _ = az;
                        let _ = new_val;
                        let candidate = if phi_c >= 0.0 {
                            axm.max(0.0).min(axp.max(0.0)).min(phi_c).max(0.0) + dx
                        } else {
                            (axm.min(0.0).max(axp.min(0.0)).max(phi_c).min(0.0)) - dx
                        };
                        if (candidate.abs() < phi_c.abs()) && (candidate * phi_c >= 0.0) {
                            field.set(ix, iy, iz, candidate);
                        }
                    }
                }
            }
        }
    }
}
/// Output triangle from marching cubes isosurface extraction.
#[derive(Debug, Clone)]
pub struct McTriangle {
    /// Vertex positions (3 vertices, each an \[x, y, z\] world-space coordinate).
    pub vertices: [[f64; 3]; 3],
}
/// Isosurface extraction via the Marching Cubes algorithm.
///
/// Implements the full 256-case lookup table for robust triangle generation
/// at the specified isovalue.
pub struct IsosurfaceExtraction {
    /// Isovalue defining the surface (default: 0.0 for signed distance fields).
    pub isovalue: f64,
}
impl IsosurfaceExtraction {
    /// Create a new [`IsosurfaceExtraction`] with the given isovalue.
    ///
    /// # Arguments
    /// * `isovalue` – the scalar value defining the isosurface
    pub fn new(isovalue: f64) -> Self {
        Self { isovalue }
    }
    /// Extract the isosurface from a [`LevelSetField`] as a list of triangles.
    ///
    /// Iterates over all cells in the grid, classifies each cube, and uses the
    /// Marching Cubes lookup table to generate triangles.
    pub fn extract(&self, field: &LevelSetField) -> Vec<McTriangle> {
        let mut triangles = Vec::new();
        let nx = field.nx;
        let ny = field.ny;
        let nz = field.nz;
        if nx < 2 || ny < 2 || nz < 2 {
            return triangles;
        }
        for iz in 0..nz - 1 {
            for iy in 0..ny - 1 {
                for ix in 0..nx - 1 {
                    self.process_cell(field, ix, iy, iz, &mut triangles);
                }
            }
        }
        triangles
    }
    /// Process a single cube cell and append any generated triangles.
    fn process_cell(
        &self,
        field: &LevelSetField,
        ix: usize,
        iy: usize,
        iz: usize,
        triangles: &mut Vec<McTriangle>,
    ) {
        let v: [f64; 8] = [
            field.get(ix, iy, iz),
            field.get(ix + 1, iy, iz),
            field.get(ix + 1, iy + 1, iz),
            field.get(ix, iy + 1, iz),
            field.get(ix, iy, iz + 1),
            field.get(ix + 1, iy, iz + 1),
            field.get(ix + 1, iy + 1, iz + 1),
            field.get(ix, iy + 1, iz + 1),
        ];
        let p: [[f64; 3]; 8] = [
            field.node_pos(ix, iy, iz),
            field.node_pos(ix + 1, iy, iz),
            field.node_pos(ix + 1, iy + 1, iz),
            field.node_pos(ix, iy + 1, iz),
            field.node_pos(ix, iy, iz + 1),
            field.node_pos(ix + 1, iy, iz + 1),
            field.node_pos(ix + 1, iy + 1, iz + 1),
            field.node_pos(ix, iy + 1, iz + 1),
        ];
        let mut cube_idx = 0u8;
        for (i, &vi) in v.iter().enumerate() {
            if vi < self.isovalue {
                cube_idx |= 1 << i;
            }
        }
        if cube_idx == 0 || cube_idx == 255 {
            return;
        }
        let edge_flags = MC_EDGE_TABLE[cube_idx as usize];
        let mut edge_verts = [[0.0f64; 3]; 12];
        for e in 0..12 {
            if edge_flags & (1 << e) != 0 {
                let (a, b) = MC_EDGE_CORNERS[e];
                edge_verts[e] = interp_vertex(p[a], p[b], v[a], v[b], self.isovalue);
            }
        }
        let tris = &MC_TRI_TABLE[cube_idx as usize];
        let mut i = 0;
        while i < tris.len() && tris[i] != -1 {
            let tri = McTriangle {
                vertices: [
                    edge_verts[tris[i] as usize],
                    edge_verts[tris[i + 1] as usize],
                    edge_verts[tris[i + 2] as usize],
                ],
            };
            triangles.push(tri);
            i += 3;
        }
    }
    /// Count the number of triangles in the isosurface.
    pub fn count_triangles(&self, field: &LevelSetField) -> usize {
        self.extract(field).len()
    }
}
