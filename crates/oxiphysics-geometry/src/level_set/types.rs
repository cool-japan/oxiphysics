//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
use crate::mc_tables::{MC_EDGE_ASCENDING, MC_EDGE_TABLE, case_index, case_triangles};

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
    /// Computed from second-order finite differences of the grid values with all six
    /// second derivatives (including the mixed ones φxy, φxz, φyz):
    /// κ = \[φxx(φy² + φz²) + φyy(φx² + φz²) + φzz(φx² + φy²) − 2(φxφyφxy +
    /// φxφzφxz + φyφzφyz)\] / |∇φ|³.  On the grid boundary the stencils become
    /// one-sided.  A sphere signed-distance field `|x| − r` has κ = 2 / |x|.
    ///
    /// Returns `0.0` for a node outside the grid and where ∇φ vanishes.
    pub fn mean_curvature(&self, ix: usize, iy: usize, iz: usize) -> f64 {
        super::finite_difference::mean_curvature_at(self, ix, iy, iz)
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
/// Level-set reinitialization (redistancing).
///
/// Turns a level-set function into a signed distance function with `|∇φ| = 1`
/// while keeping the sign of every node and, to the accuracy of the schemes, the
/// zero level set.  Two solvers are available: pseudo-time integration of the
/// Sussman reinitialization equation ([`reinitialize`](Self::reinitialize)) and
/// fast sweeping of the Eikonal equation ([`fast_sweep`](Self::fast_sweep)).
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
    /// Fast-sweeping redistancing (Zhao 2005): turn `field` into a signed distance
    /// function with `|∇φ| = 1` while keeping the sign of every node.
    ///
    /// Nodes that have a neighbour of opposite sign (and nodes with φ = 0 exactly)
    /// are frozen with the first-order distance estimate `|φ| / |∇φ|`, limited by
    /// the distance to the closest zero crossing along a grid edge.  All other
    /// nodes start at a huge distance and are lowered by Gauss–Seidel sweeps with
    /// the upwind Godunov update of the Eikonal equation `|∇d| = 1`, using the
    /// nearest known distance along each of the x, y and z axes.  Sweep `s` visits
    /// the grid in the diagonal ordering `s mod 8`, so `n_sweeps = 8` is one full
    /// Gauss–Seidel iteration; the loop ends early once a whole cycle of eight
    /// orderings changes no node by more than `self.tol`.
    ///
    /// Nodes that the requested sweeps have not reached keep their original value;
    /// a field without a zero level set, or one containing NaN, is left unchanged.
    ///
    /// # Arguments
    /// * `field`    – the level-set field to reinitialize
    /// * `n_sweeps` – maximum number of Gauss–Seidel sweeps
    pub fn fast_sweep(&self, field: &mut LevelSetField, n_sweeps: usize) {
        super::fast_sweep::fast_sweep(self, field, n_sweeps);
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
/// at the specified isovalue.  A grid node is *inside* when its value is
/// strictly below the isovalue; the emitted triangles are wound counter-clockwise
/// when seen from outside, i.e. their normals point towards increasing field
/// values.  Crossing points are interpolated from the lower to the upper end of
/// every cube edge, so the triangles of neighbouring cells meet at bit-identical
/// vertices.
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
        let case = case_index(&v, self.isovalue);
        if case == 0 || case == 255 {
            return;
        }
        let edge_flags = MC_EDGE_TABLE[case];
        let mut edge_verts = [[0.0f64; 3]; 12];
        for (e, &[lo, hi]) in MC_EDGE_ASCENDING.iter().enumerate() {
            if edge_flags & (1 << e) != 0 {
                edge_verts[e] = interp_vertex(p[lo], p[hi], v[lo], v[hi], self.isovalue);
            }
        }
        for [a, b, c] in case_triangles(case) {
            // The table winds into the solid; swap two vertices so the triangle
            // faces out of it.
            triangles.push(McTriangle {
                vertices: [edge_verts[a], edge_verts[c], edge_verts[b]],
            });
        }
    }
    /// Count the number of triangles in the isosurface.
    pub fn count_triangles(&self, field: &LevelSetField) -> usize {
        self.extract(field).len()
    }
}
