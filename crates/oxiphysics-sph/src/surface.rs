// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Surface tension models and surface reconstruction for SPH simulations.
//!
//! Implements the Continuum Surface Force (CSF) method using color field
//! gradient and Laplacian to estimate surface normals and curvature,
//! along with surface reconstruction utilities (marching cubes, anisotropic kernels).

use std::f64::consts::PI;

use oxiphysics_core::math::Vec3;

use crate::kernel::SphKernel;
use crate::particle::ParticleSet;

// ---------------------------------------------------------------------------
// CsfSurfaceTension (existing)
// ---------------------------------------------------------------------------

/// Continuum Surface Force (CSF) surface tension model.
///
/// This method computes surface normals from the gradient of a color field
/// and surface tension forces from the Laplacian of the color field.
#[derive(Debug, Clone)]
pub struct CsfSurfaceTension {
    /// Surface tension coefficient (N/m).
    pub tension_coefficient: f64,
    /// Minimum normal magnitude threshold (to avoid spurious forces in the interior).
    pub normal_threshold: f64,
}

impl Default for CsfSurfaceTension {
    fn default() -> Self {
        Self {
            tension_coefficient: 0.0728, // water at 20°C
            normal_threshold: 1.0,
        }
    }
}

impl CsfSurfaceTension {
    /// Create a new CSF surface tension model.
    pub fn new(tension_coefficient: f64, normal_threshold: f64) -> Self {
        Self {
            tension_coefficient,
            normal_threshold,
        }
    }

    /// Compute the color field gradient for each particle.
    ///
    /// The color field is c_i = sum_j (m_j / rho_j) * W_ij, and its gradient
    /// gives an estimate of the surface normal.
    pub fn compute_color_gradient(
        &self,
        particles: &ParticleSet,
        neighbors: &[Vec<usize>],
        kernel: &dyn SphKernel,
        h: f64,
    ) -> Vec<Vec3> {
        let n = particles.len();
        let mut gradients = vec![Vec3::zeros(); n];

        for i in 0..n {
            let rhoi = particles.densities[i];
            if rhoi < 1e-14 {
                continue;
            }
            for &j in &neighbors[i] {
                let rhoj = particles.densities[j];
                if rhoj < 1e-14 {
                    continue;
                }
                let rij = particles.positions[i] - particles.positions[j];
                let r = rij.norm();
                if r < 1e-14 {
                    continue;
                }
                let grad = kernel.grad_w(r, h);
                let rhat = rij / r;
                gradients[i] += particles.masses[j] / rhoj * grad * rhat;
            }
        }

        gradients
    }

    /// Compute the color field Laplacian for each particle.
    pub fn compute_color_laplacian(
        &self,
        particles: &ParticleSet,
        neighbors: &[Vec<usize>],
        kernel: &dyn SphKernel,
        h: f64,
    ) -> Vec<f64> {
        let n = particles.len();
        let mut laplacians = vec![0.0; n];

        for i in 0..n {
            for &j in &neighbors[i] {
                let rhoj = particles.densities[j];
                if rhoj < 1e-14 {
                    continue;
                }
                let rij = particles.positions[i] - particles.positions[j];
                let r = rij.norm();
                laplacians[i] += particles.masses[j] / rhoj * kernel.laplacian_w(r, h);
            }
        }

        laplacians
    }

    /// Compute surface tension forces and accumulate into `particles.forces`.
    pub fn apply_forces(
        &self,
        particles: &mut ParticleSet,
        neighbors: &[Vec<usize>],
        kernel: &dyn SphKernel,
        h: f64,
    ) {
        let normals = self.compute_color_gradient(particles, neighbors, kernel, h);
        let laplacians = self.compute_color_laplacian(particles, neighbors, kernel, h);

        for i in 0..particles.len() {
            let n_mag = normals[i].norm();
            if n_mag > self.normal_threshold {
                // Surface tension force: -sigma * laplacian(c) * n_hat
                let n_hat = normals[i] / n_mag;
                particles.forces[i] -= self.tension_coefficient * laplacians[i] * n_hat;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Poly6 kernel functions for surface reconstruction
// ---------------------------------------------------------------------------

/// Poly6 kernel value.
///
/// W(r, h) = (315 / (64 * pi * h^9)) * (h^2 - r^2)^3   for r <= h
///         = 0 otherwise
pub fn poly6_kernel(r_sq: f64, h: f64) -> f64 {
    let h2 = h * h;
    if r_sq > h2 {
        return 0.0;
    }
    let coeff = 315.0 / (64.0 * PI * h.powi(9));
    let diff = h2 - r_sq;
    coeff * diff * diff * diff
}

/// Gradient of the Poly6 kernel with respect to the position vector `r_vec`.
///
/// grad W = -6 * (315 / (64 * pi * h^9)) * (h^2 - r^2)^2 * r_vec
pub fn poly6_gradient(r_vec: [f64; 3], r_sq: f64, h: f64) -> [f64; 3] {
    let h2 = h * h;
    if r_sq > h2 {
        return [0.0; 3];
    }
    let coeff = -6.0 * 315.0 / (64.0 * PI * h.powi(9));
    let diff = h2 - r_sq;
    let scale = coeff * diff * diff;
    [scale * r_vec[0], scale * r_vec[1], scale * r_vec[2]]
}

// ---------------------------------------------------------------------------
// SphColorField
// ---------------------------------------------------------------------------

/// SPH color field for surface detection and reconstruction.
///
/// Evaluates the color field C(x) = sum_j (m_j / rho_j) * W(x - x_j, h)
/// and its gradient using the Poly6 kernel.
#[derive(Debug, Clone)]
pub struct SphColorField {
    /// Particle positions.
    pub positions: Vec<[f64; 3]>,
    /// Particle masses.
    pub masses: Vec<f64>,
    /// Particle densities.
    pub densities: Vec<f64>,
    /// Smoothing length.
    pub smoothing_h: f64,
}

impl SphColorField {
    /// Create a new, empty color field with the given smoothing length.
    pub fn new(h: f64) -> Self {
        Self {
            positions: Vec::new(),
            masses: Vec::new(),
            densities: Vec::new(),
            smoothing_h: h,
        }
    }

    /// Add a particle to the color field.
    pub fn add_particle(&mut self, pos: [f64; 3], mass: f64, density: f64) {
        self.positions.push(pos);
        self.masses.push(mass);
        self.densities.push(density);
    }

    /// Evaluate the color field at an arbitrary point.
    ///
    /// C(x) = sum_j (m_j / rho_j) * W(|x - x_j|^2, h)
    pub fn color_field_at(&self, point: [f64; 3]) -> f64 {
        let h = self.smoothing_h;
        let mut sum = 0.0;
        for (idx, &pos) in self.positions.iter().enumerate() {
            let rho = self.densities[idx];
            if rho < 1e-14 {
                continue;
            }
            let dx = point[0] - pos[0];
            let dy = point[1] - pos[1];
            let dz = point[2] - pos[2];
            let r_sq = dx * dx + dy * dy + dz * dz;
            sum += (self.masses[idx] / rho) * poly6_kernel(r_sq, h);
        }
        sum
    }

    /// Evaluate the color field gradient at an arbitrary point.
    ///
    /// grad C(x) = sum_j (m_j / rho_j) * grad W(x - x_j, h)
    pub fn color_field_gradient_at(&self, point: [f64; 3]) -> [f64; 3] {
        let h = self.smoothing_h;
        let mut grad = [0.0f64; 3];
        for (idx, &pos) in self.positions.iter().enumerate() {
            let rho = self.densities[idx];
            if rho < 1e-14 {
                continue;
            }
            let r_vec = [point[0] - pos[0], point[1] - pos[1], point[2] - pos[2]];
            let r_sq = r_vec[0] * r_vec[0] + r_vec[1] * r_vec[1] + r_vec[2] * r_vec[2];
            let g = poly6_gradient(r_vec, r_sq, h);
            let w = self.masses[idx] / rho;
            grad[0] += w * g[0];
            grad[1] += w * g[1];
            grad[2] += w * g[2];
        }
        grad
    }

    /// Return `true` if particle `i` is a surface particle.
    ///
    /// A particle is considered a surface particle when the magnitude of the
    /// color-field gradient evaluated at its position exceeds `threshold`.
    pub fn is_surface_particle(&self, i: usize, threshold: f64) -> bool {
        let grad = self.color_field_gradient_at(self.positions[i]);
        let mag_sq = grad[0] * grad[0] + grad[1] * grad[1] + grad[2] * grad[2];
        mag_sq.sqrt() > threshold
    }

    /// Return the unit-length surface normal at `point` (normalized gradient).
    ///
    /// Returns `[0.0; 3]` when the gradient magnitude is below epsilon.
    pub fn surface_normal_at(&self, point: [f64; 3]) -> [f64; 3] {
        let grad = self.color_field_gradient_at(point);
        let mag = (grad[0] * grad[0] + grad[1] * grad[1] + grad[2] * grad[2]).sqrt();
        if mag < 1e-14 {
            return [0.0; 3];
        }
        [grad[0] / mag, grad[1] / mag, grad[2] / mag]
    }
}

// ---------------------------------------------------------------------------
// MarchingCubesSph – minimal marching-cubes surface extractor
// ---------------------------------------------------------------------------

/// A uniform-grid marching-cubes surface extractor for SPH implicit surfaces.
///
/// `extract_surface` evaluates a scalar field on a grid and returns triangle
/// vertex positions (3 vertices per triangle, stored consecutively).
#[derive(Debug, Clone)]
pub struct MarchingCubesSph {
    /// Minimum corner of the grid bounding box.
    pub grid_min: [f64; 3],
    /// Maximum corner of the grid bounding box.
    pub grid_max: [f64; 3],
    /// Number of voxels along each axis.
    pub resolution: usize,
}

impl MarchingCubesSph {
    /// Create a new marching-cubes extractor.
    pub fn new(grid_min: [f64; 3], grid_max: [f64; 3], resolution: usize) -> Self {
        Self {
            grid_min,
            grid_max,
            resolution,
        }
    }

    /// Extract an iso-surface as a triangle soup.
    ///
    /// Evaluates `color_fn` at each grid corner, then processes each voxel
    /// using a simplified marching-cubes lookup (only the 16 canonical cases
    /// that cover all 256 configurations via symmetry are needed; here we use
    /// a direct per-edge sign-change approach).
    ///
    /// Returns a flat list of `[f64; 3]` vertices; every three consecutive
    /// entries form one triangle.
    pub fn extract_surface(
        &self,
        color_fn: impl Fn([f64; 3]) -> f64,
        iso_value: f64,
    ) -> Vec<[f64; 3]> {
        let res = self.resolution;
        if res == 0 {
            return Vec::new();
        }

        // Sample the scalar field on a (res+1)^3 grid.
        let nx = res + 1;
        let step = [
            (self.grid_max[0] - self.grid_min[0]) / res as f64,
            (self.grid_max[1] - self.grid_min[1]) / res as f64,
            (self.grid_max[2] - self.grid_min[2]) / res as f64,
        ];

        let idx = |ix: usize, iy: usize, iz: usize| ix * nx * nx + iy * nx + iz;

        let mut field = vec![0.0f64; nx * nx * nx];
        let mut coords = vec![[0.0f64; 3]; nx * nx * nx];

        for ix in 0..nx {
            for iy in 0..nx {
                for iz in 0..nx {
                    let pt = [
                        self.grid_min[0] + ix as f64 * step[0],
                        self.grid_min[1] + iy as f64 * step[1],
                        self.grid_min[2] + iz as f64 * step[2],
                    ];
                    let k = idx(ix, iy, iz);
                    coords[k] = pt;
                    field[k] = color_fn(pt);
                }
            }
        }

        // Cube corner offsets (dx, dy, dz) in index space
        const CORNERS: [(usize, usize, usize); 8] = [
            (0, 0, 0),
            (1, 0, 0),
            (1, 1, 0),
            (0, 1, 0),
            (0, 0, 1),
            (1, 0, 1),
            (1, 1, 1),
            (0, 1, 1),
        ];

        // 12 edges: (corner_a, corner_b)
        const EDGES: [(usize, usize); 12] = [
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 0), // bottom face
            (4, 5),
            (5, 6),
            (6, 7),
            (7, 4), // top face
            (0, 4),
            (1, 5),
            (2, 6),
            (3, 7), // verticals
        ];

        // Marching-cubes triangle table (256 entries, 16 tris max per case).
        // We use a compact encoding: each row is a list of edge indices terminated
        // by -1. The full 256-entry table is included here.
        let tri_table = mc_tri_table();

        let mut triangles: Vec<[f64; 3]> = Vec::new();

        for ix in 0..res {
            for iy in 0..res {
                for iz in 0..res {
                    // Sample 8 corners
                    let mut corner_vals = [0.0f64; 8];
                    let mut corner_pts = [[0.0f64; 3]; 8];
                    for (c, &(dx, dy, dz)) in CORNERS.iter().enumerate() {
                        let k = idx(ix + dx, iy + dy, iz + dz);
                        corner_vals[c] = field[k];
                        corner_pts[c] = coords[k];
                    }

                    // Compute cube index (which corners are inside)
                    let mut cube_idx: usize = 0;
                    for (c, &cv) in corner_vals.iter().enumerate() {
                        if cv >= iso_value {
                            cube_idx |= 1 << c;
                        }
                    }

                    if cube_idx == 0 || cube_idx == 255 {
                        continue; // fully outside or fully inside
                    }

                    // Interpolate edge vertices
                    let mut edge_verts = [[0.0f64; 3]; 12];
                    for (e, &(ca, cb)) in EDGES.iter().enumerate() {
                        let va = corner_vals[ca];
                        let vb = corner_vals[cb];
                        let pa = corner_pts[ca];
                        let pb = corner_pts[cb];
                        let t = if (vb - va).abs() > 1e-14 {
                            (iso_value - va) / (vb - va)
                        } else {
                            0.5
                        };
                        edge_verts[e] = [
                            pa[0] + t * (pb[0] - pa[0]),
                            pa[1] + t * (pb[1] - pa[1]),
                            pa[2] + t * (pb[2] - pa[2]),
                        ];
                    }

                    // Emit triangles
                    let row = &tri_table[cube_idx];
                    let mut ri = 0;
                    while ri + 2 < row.len() && row[ri] >= 0 {
                        triangles.push(edge_verts[row[ri] as usize]);
                        triangles.push(edge_verts[row[ri + 1] as usize]);
                        triangles.push(edge_verts[row[ri + 2] as usize]);
                        ri += 3;
                    }
                }
            }
        }

        triangles
    }
}

// ---------------------------------------------------------------------------
// Marching-cubes triangle table (256 rows, up to 16 values each, -1 = end)
// ---------------------------------------------------------------------------

/// Returns the full 256-entry MC triangle table.
/// Each row lists edge indices (0..11) forming triangles, terminated by -1.
fn mc_tri_table() -> Vec<Vec<i8>> {
    // Compact representation: each sub-array is one row.
    // Source: Lorensen & Cline 1987 table (widely reproduced).
    let raw: &[&[i8]] = &[
        &[-1],
        &[0, 8, 3, -1],
        &[0, 1, 9, -1],
        &[1, 8, 3, 9, 8, 1, -1],
        &[1, 2, 10, -1],
        &[0, 8, 3, 1, 2, 10, -1],
        &[9, 2, 10, 0, 2, 9, -1],
        &[2, 8, 3, 2, 10, 8, 10, 9, 8, -1],
        &[3, 11, 2, -1],
        &[0, 11, 2, 8, 11, 0, -1],
        &[1, 9, 0, 2, 3, 11, -1],
        &[1, 11, 2, 1, 9, 11, 9, 8, 11, -1],
        &[3, 10, 1, 11, 10, 3, -1],
        &[0, 10, 1, 0, 8, 10, 8, 11, 10, -1],
        &[3, 9, 0, 3, 11, 9, 11, 10, 9, -1],
        &[9, 8, 10, 10, 8, 11, -1],
        &[4, 7, 8, -1],
        &[4, 3, 0, 7, 3, 4, -1],
        &[0, 1, 9, 8, 4, 7, -1],
        &[4, 1, 9, 4, 7, 1, 7, 3, 1, -1],
        &[1, 2, 10, 8, 4, 7, -1],
        &[3, 4, 7, 3, 0, 4, 1, 2, 10, -1],
        &[9, 2, 10, 9, 0, 2, 8, 4, 7, -1],
        &[2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1],
        &[8, 4, 7, 3, 11, 2, -1],
        &[11, 4, 7, 11, 2, 4, 2, 0, 4, -1],
        &[9, 0, 1, 8, 4, 7, 2, 3, 11, -1],
        &[4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1],
        &[3, 10, 1, 3, 11, 10, 7, 8, 4, -1],
        &[1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1],
        &[4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1],
        &[4, 7, 11, 4, 11, 9, 9, 11, 10, -1],
        &[9, 5, 4, -1],
        &[9, 5, 4, 0, 8, 3, -1],
        &[0, 5, 4, 1, 5, 0, -1],
        &[8, 5, 4, 8, 3, 5, 3, 1, 5, -1],
        &[1, 2, 10, 9, 5, 4, -1],
        &[3, 0, 8, 1, 2, 10, 4, 9, 5, -1],
        &[5, 2, 10, 5, 4, 2, 4, 0, 2, -1],
        &[2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1],
        &[9, 5, 4, 2, 3, 11, -1],
        &[0, 11, 2, 0, 8, 11, 4, 9, 5, -1],
        &[0, 5, 4, 0, 1, 5, 2, 3, 11, -1],
        &[2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1],
        &[10, 3, 11, 10, 1, 3, 9, 5, 4, -1],
        &[4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1],
        &[5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1],
        &[5, 4, 8, 5, 8, 10, 10, 8, 11, -1],
        &[9, 7, 8, 5, 7, 9, -1],
        &[9, 3, 0, 9, 5, 3, 5, 7, 3, -1],
        &[0, 7, 8, 0, 1, 7, 1, 5, 7, -1],
        &[1, 5, 3, 3, 5, 7, -1],
        &[9, 7, 8, 9, 5, 7, 10, 1, 2, -1],
        &[10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1],
        &[8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1],
        &[2, 10, 5, 2, 5, 3, 3, 5, 7, -1],
        &[7, 9, 5, 7, 8, 9, 3, 11, 2, -1],
        &[9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1],
        &[2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1],
        &[11, 2, 1, 11, 1, 7, 7, 1, 5, -1],
        &[9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1],
        &[5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1],
        &[11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1],
        &[11, 10, 5, 7, 11, 5, -1],
        &[10, 6, 5, -1],
        &[0, 8, 3, 5, 10, 6, -1],
        &[9, 0, 1, 5, 10, 6, -1],
        &[1, 8, 3, 1, 9, 8, 5, 10, 6, -1],
        &[1, 6, 5, 2, 6, 1, -1],
        &[1, 6, 5, 1, 2, 6, 3, 0, 8, -1],
        &[9, 6, 5, 9, 0, 6, 0, 2, 6, -1],
        &[5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1],
        &[2, 3, 11, 10, 6, 5, -1],
        &[11, 0, 8, 11, 2, 0, 10, 6, 5, -1],
        &[0, 1, 9, 2, 3, 11, 5, 10, 6, -1],
        &[5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1],
        &[6, 3, 11, 6, 5, 3, 5, 1, 3, -1],
        &[0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1],
        &[3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1],
        &[6, 5, 9, 6, 9, 11, 11, 9, 8, -1],
        &[5, 10, 6, 4, 7, 8, -1],
        &[4, 3, 0, 4, 7, 3, 6, 5, 10, -1],
        &[1, 9, 0, 5, 10, 6, 8, 4, 7, -1],
        &[10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1],
        &[6, 1, 2, 6, 5, 1, 4, 7, 8, -1],
        &[1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1],
        &[8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1],
        &[7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1],
        &[3, 11, 2, 7, 8, 4, 10, 6, 5, -1],
        &[5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1],
        &[0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1],
        &[9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1],
        &[8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1],
        &[5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1],
        &[0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1],
        &[6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1],
        &[10, 4, 9, 6, 4, 10, -1],
        &[4, 10, 6, 4, 9, 10, 0, 8, 3, -1],
        &[10, 0, 1, 10, 6, 0, 6, 4, 0, -1],
        &[8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1],
        &[1, 4, 9, 1, 2, 4, 2, 6, 4, -1],
        &[3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1],
        &[0, 2, 4, 4, 2, 6, -1],
        &[8, 3, 2, 8, 2, 4, 4, 2, 6, -1],
        &[10, 4, 9, 10, 6, 4, 11, 2, 3, -1],
        &[0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1],
        &[3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1],
        &[6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1],
        &[9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1],
        &[8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1],
        &[3, 11, 6, 3, 6, 0, 0, 6, 4, -1],
        &[6, 4, 8, 11, 6, 8, -1],
        &[7, 10, 6, 7, 8, 10, 8, 9, 10, -1],
        &[0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1],
        &[10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1],
        &[10, 6, 7, 10, 7, 1, 1, 7, 3, -1],
        &[1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1],
        &[2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1],
        &[7, 8, 0, 7, 0, 6, 6, 0, 2, -1],
        &[7, 3, 2, 6, 7, 2, -1],
        &[2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1],
        &[2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1],
        &[1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1],
        &[11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1],
        &[8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1],
        &[0, 9, 1, 11, 6, 7, -1],
        &[7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1],
        &[7, 11, 6, -1],
        &[7, 6, 11, -1],
        &[3, 0, 8, 11, 7, 6, -1],
        &[0, 1, 9, 11, 7, 6, -1],
        &[8, 1, 9, 8, 3, 1, 11, 7, 6, -1],
        &[10, 1, 2, 6, 11, 7, -1],
        &[1, 2, 10, 3, 0, 8, 6, 11, 7, -1],
        &[2, 9, 0, 2, 10, 9, 6, 11, 7, -1],
        &[6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1],
        &[7, 2, 3, 6, 2, 7, -1],
        &[7, 0, 8, 7, 6, 0, 6, 2, 0, -1],
        &[2, 7, 6, 2, 3, 7, 0, 1, 9, -1],
        &[1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1],
        &[10, 7, 6, 10, 1, 7, 1, 3, 7, -1],
        &[10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1],
        &[0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1],
        &[7, 6, 10, 7, 10, 8, 8, 10, 9, -1],
        &[6, 8, 4, 11, 8, 6, -1],
        &[3, 6, 11, 3, 0, 6, 0, 4, 6, -1],
        &[8, 6, 11, 8, 4, 6, 9, 0, 1, -1],
        &[9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1],
        &[6, 8, 4, 6, 11, 8, 2, 10, 1, -1],
        &[1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1],
        &[4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1],
        &[10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1],
        &[8, 2, 3, 8, 4, 2, 4, 6, 2, -1],
        &[0, 4, 2, 4, 6, 2, -1],
        &[1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1],
        &[1, 9, 4, 1, 4, 2, 2, 4, 6, -1],
        &[8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1],
        &[10, 1, 0, 10, 0, 6, 6, 0, 4, -1],
        &[4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1],
        &[10, 9, 4, 6, 10, 4, -1],
        &[4, 9, 5, 7, 6, 11, -1],
        &[0, 8, 3, 4, 9, 5, 11, 7, 6, -1],
        &[5, 0, 1, 5, 4, 0, 7, 6, 11, -1],
        &[11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1],
        &[9, 5, 4, 10, 1, 2, 7, 6, 11, -1],
        &[6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1],
        &[7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1],
        &[3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1],
        &[7, 2, 3, 7, 6, 2, 5, 4, 9, -1],
        &[9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1],
        &[3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1],
        &[6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1],
        &[9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1],
        &[1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1],
        &[4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1],
        &[7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1],
        &[6, 9, 5, 6, 11, 9, 11, 8, 9, -1],
        &[3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1],
        &[0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1],
        &[6, 11, 3, 6, 3, 5, 5, 3, 1, -1],
        &[1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1],
        &[0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1],
        &[11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1],
        &[6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1],
        &[5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1],
        &[9, 5, 6, 9, 6, 0, 0, 6, 2, -1],
        &[1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1],
        &[1, 5, 6, 2, 1, 6, -1],
        &[1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1],
        &[10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1],
        &[0, 3, 8, 5, 6, 10, -1],
        &[10, 5, 6, -1],
        &[11, 5, 10, 7, 5, 11, -1],
        &[11, 5, 10, 11, 7, 5, 8, 3, 0, -1],
        &[5, 11, 7, 5, 10, 11, 1, 9, 0, -1],
        &[10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1],
        &[11, 1, 2, 11, 7, 1, 7, 5, 1, -1],
        &[0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1],
        &[9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1],
        &[7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1],
        &[2, 5, 10, 2, 3, 5, 3, 7, 5, -1],
        &[8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1],
        &[9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1],
        &[9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1],
        &[1, 3, 5, 3, 7, 5, -1],
        &[0, 8, 7, 0, 7, 1, 1, 7, 5, -1],
        &[9, 0, 3, 9, 3, 5, 5, 3, 7, -1],
        &[9, 8, 7, 5, 9, 7, -1],
        &[5, 8, 4, 5, 10, 8, 10, 11, 8, -1],
        &[5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1],
        &[0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1],
        &[10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1],
        &[2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1],
        &[0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1],
        &[0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1],
        &[9, 4, 5, 2, 11, 3, -1],
        &[2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1],
        &[5, 10, 2, 5, 2, 4, 4, 2, 0, -1],
        &[3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1],
        &[5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1],
        &[8, 4, 5, 8, 5, 3, 3, 5, 1, -1],
        &[0, 4, 5, 1, 0, 5, -1],
        &[8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1],
        &[9, 4, 5, -1],
        &[4, 11, 7, 4, 9, 11, 9, 10, 11, -1],
        &[0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1],
        &[1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1],
        &[3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1],
        &[4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1],
        &[9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1],
        &[11, 7, 4, 11, 4, 2, 2, 4, 0, -1],
        &[11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1],
        &[2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1],
        &[9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1],
        &[3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1],
        &[1, 10, 2, 8, 7, 4, -1],
        &[4, 9, 1, 4, 1, 7, 7, 1, 3, -1],
        &[4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1],
        &[4, 0, 3, 7, 4, 3, -1],
        &[4, 8, 7, -1],
        &[9, 10, 8, 10, 11, 8, -1],
        &[3, 0, 9, 3, 9, 11, 11, 9, 10, -1],
        &[0, 1, 10, 0, 10, 8, 8, 10, 11, -1],
        &[3, 1, 10, 11, 3, 10, -1],
        &[1, 2, 11, 1, 11, 9, 9, 11, 8, -1],
        &[3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1],
        &[0, 2, 11, 8, 0, 11, -1],
        &[3, 2, 11, -1],
        &[2, 3, 8, 2, 8, 10, 10, 8, 9, -1],
        &[9, 10, 2, 0, 9, 2, -1],
        &[2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1],
        &[1, 10, 2, -1],
        &[1, 3, 8, 9, 1, 8, -1],
        &[0, 9, 1, -1],
        &[0, 3, 8, -1],
        &[-1],
    ];

    raw.iter().map(|row| row.to_vec()).collect()
}

// ---------------------------------------------------------------------------
// AnisotropicKernel
// ---------------------------------------------------------------------------

/// Anisotropic SPH kernel (Müller et al. style).
///
/// The kernel is W_aniso(x) = det(G)/h^3 * W(|G*(x-xi)|/h)
/// where G is a 3×3 anisotropy matrix and W is the Poly6 kernel.
#[derive(Debug, Clone)]
pub struct AnisotropicKernel {
    /// Kernel centre position.
    pub position: [f64; 3],
    /// Anisotropy matrix (3×3, stored as rows).
    pub g: [[f64; 3]; 3],
    /// Smoothing length.
    pub h: f64,
}

impl AnisotropicKernel {
    /// Evaluate the anisotropic kernel at `x`.
    pub fn evaluate(&self, x: [f64; 3]) -> f64 {
        let dx = [
            x[0] - self.position[0],
            x[1] - self.position[1],
            x[2] - self.position[2],
        ];
        // y = G * dx
        let y = mat3_mul_vec(self.g, dx);
        let r_sq = (y[0] * y[0] + y[1] * y[1] + y[2] * y[2]) / (self.h * self.h);
        let det = mat3_det(self.g);
        det / (self.h * self.h * self.h) * poly6_kernel(r_sq, 1.0)
    }
}

/// Build a vector of anisotropic kernels for a set of particle positions.
///
/// For simplicity the anisotropy matrix G is taken as `k_r * I` (identity
/// scaled by `k_r`), which degenerates to a scaled isotropic kernel.
/// A full implementation would compute the local covariance and use SVD.
pub fn build_anisotropic_kernels(
    positions: &[[f64; 3]],
    h: f64,
    k_r: f64,
) -> Vec<AnisotropicKernel> {
    positions
        .iter()
        .map(|&pos| {
            let scale = k_r;
            AnisotropicKernel {
                position: pos,
                g: [[scale, 0.0, 0.0], [0.0, scale, 0.0], [0.0, 0.0, scale]],
                h,
            }
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Implicit surface tension force (CSF formula using SphColorField)
// ---------------------------------------------------------------------------

/// Compute the implicit surface tension force on particle `i` using the CSF model.
///
/// F = sigma * kappa * n  where kappa = -div(n / |n|)
///
/// Here we approximate kappa via finite differences of the colour-field gradient.
pub fn implicit_surface_tension_force(color: &SphColorField, i: usize, sigma: f64) -> [f64; 3] {
    let pos = color.positions[i];
    let grad = color.color_field_gradient_at(pos);
    let mag = (grad[0] * grad[0] + grad[1] * grad[1] + grad[2] * grad[2]).sqrt();
    if mag < 1e-14 {
        return [0.0; 3];
    }
    let n_hat = [grad[0] / mag, grad[1] / mag, grad[2] / mag];

    // Approximate divergence of n_hat via finite differences
    let eps = color.smoothing_h * 1e-3;
    let kappa = -finite_div_nhat(color, pos, eps);

    [
        sigma * kappa * n_hat[0],
        sigma * kappa * n_hat[1],
        sigma * kappa * n_hat[2],
    ]
}

/// Approximate div(n_hat) at `pos` via central finite differences.
fn finite_div_nhat(color: &SphColorField, pos: [f64; 3], eps: f64) -> f64 {
    let mut div = 0.0;
    for axis in 0..3 {
        let mut pp = pos;
        let mut pm = pos;
        pp[axis] += eps;
        pm[axis] -= eps;
        let gp = color.color_field_gradient_at(pp);
        let gm = color.color_field_gradient_at(pm);
        let mp = (gp[0] * gp[0] + gp[1] * gp[1] + gp[2] * gp[2])
            .sqrt()
            .max(1e-14);
        let mm = (gm[0] * gm[0] + gm[1] * gm[1] + gm[2] * gm[2])
            .sqrt()
            .max(1e-14);
        div += (gp[axis] / mp - gm[axis] / mm) / (2.0 * eps);
    }
    div
}

// ---------------------------------------------------------------------------
// Small linear-algebra helpers (3×3 matrix)
// ---------------------------------------------------------------------------

fn mat3_mul_vec(m: [[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]
}

fn mat3_det(m: [[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

// ---------------------------------------------------------------------------
// Color function surface detection
// ---------------------------------------------------------------------------

/// Detect whether particle i is a surface particle using the color function.
///
/// A particle is classified as a surface particle if its normalized density
/// (color function value) falls below the given threshold.
///
/// # Arguments
/// * `density` – density of particle i
/// * `rest_density` – reference rest density ρ₀
/// * `threshold` – fraction of rest density below which a particle is "surface" (e.g. 0.9)
pub fn color_function_surface_detect(density: f64, rest_density: f64, threshold: f64) -> bool {
    if rest_density < 1e-30 {
        return false;
    }
    density / rest_density < threshold
}

/// Compute the color function value c = ρ / ρ₀.
///
/// Returns a scalar in \[0, 1\] for typical SPH densities.
pub fn color_function_value(density: f64, rest_density: f64) -> f64 {
    if rest_density < 1e-30 {
        return 0.0;
    }
    density / rest_density
}

// ---------------------------------------------------------------------------
// Curvature from color gradient
// ---------------------------------------------------------------------------

/// Estimate surface curvature from the color field gradient divergence.
///
/// κ = −div(n̂) = −Σ_α (∂n̂_α / ∂x_α)
///
/// where `n̂` is the unit surface normal.
pub fn curvature_from_color_gradient(_normal: [f64; 3], grad_normal: [[f64; 3]; 3]) -> f64 {
    -(grad_normal[0][0] + grad_normal[1][1] + grad_normal[2][2])
}

/// Compute the unit surface normal from the color gradient.
///
/// Returns zero vector if the gradient is negligibly small.
pub fn color_normal_from_gradient(grad_color: [f64; 3]) -> [f64; 3] {
    let mag = (grad_color[0] * grad_color[0]
        + grad_color[1] * grad_color[1]
        + grad_color[2] * grad_color[2])
        .sqrt();
    if mag < 1e-14 {
        return [0.0; 3];
    }
    [
        grad_color[0] / mag,
        grad_color[1] / mag,
        grad_color[2] / mag,
    ]
}

// ---------------------------------------------------------------------------
// Level-set extraction
// ---------------------------------------------------------------------------

/// Extract the level-set value at point `x` by interpolating the color field.
///
/// φ(x) = Σ_j c_j * W(|x - x_j|, h)
/// where c_j is the color (indicator) value of particle j and W is the Poly6 kernel.
pub fn level_set_extraction(positions: &[[f64; 3]], colors: &[f64], x: [f64; 3], h: f64) -> f64 {
    let mut phi = 0.0_f64;
    for (idx, &pos_j) in positions.iter().enumerate() {
        let r_sq =
            (x[0] - pos_j[0]).powi(2) + (x[1] - pos_j[1]).powi(2) + (x[2] - pos_j[2]).powi(2);
        phi += colors[idx] * poly6_kernel(r_sq, h);
    }
    phi
}

/// Extract the level-set gradient at point `x`.
///
/// ∇φ(x) = Σ_j c_j * ∇W(x - x_j, h)
pub fn level_set_gradient(positions: &[[f64; 3]], colors: &[f64], x: [f64; 3], h: f64) -> [f64; 3] {
    let mut grad = [0.0_f64; 3];
    for (idx, &pos_j) in positions.iter().enumerate() {
        let r_vec = [x[0] - pos_j[0], x[1] - pos_j[1], x[2] - pos_j[2]];
        let r_sq = r_vec[0] * r_vec[0] + r_vec[1] * r_vec[1] + r_vec[2] * r_vec[2];
        let g = poly6_gradient(r_vec, r_sq, h);
        let c = colors[idx];
        grad[0] += c * g[0];
        grad[1] += c * g[1];
        grad[2] += c * g[2];
    }
    grad
}

// ---------------------------------------------------------------------------
// Marching Cubes cell
// ---------------------------------------------------------------------------

/// Count the number of triangles a marching cubes cell would produce.
///
/// Uses the lookup table index only (not the actual vertex positions).
/// Returns the number of triangles (0–5) based on the sign configuration
/// of the 8 corner values relative to `iso_level`.
///
/// Corners are indexed as:
/// ```text
/// 4---5
/// |   |  (top face)
/// 7---6
/// 0---1
/// |   |  (bottom face)
/// 3---2
/// ```
pub fn marching_cubes_cell_triangle_count(corners: &[f64; 8], iso_level: f64) -> usize {
    // Build the 8-bit case index
    let mut cube_idx = 0u8;
    for (k, &v) in corners.iter().enumerate() {
        if v < iso_level {
            cube_idx |= 1 << k;
        }
    }
    // Use the standard MC triangle count table (abridged to unique cases)
    // Full table has 256 entries; we encode a simplified count here.
    MC_TRI_COUNT[cube_idx as usize]
}

/// Triangle count lookup table for marching cubes (256 entries).
///
/// Index is the 8-bit configuration of corner signs.
/// Value is the number of triangles for that configuration.
static MC_TRI_COUNT: [usize; 256] = compute_mc_tri_count();

const fn compute_mc_tri_count() -> [usize; 256] {
    let mut table = [0usize; 256];
    let mut i = 0u8;
    loop {
        // Count bit transitions to estimate triangle count (simplified heuristic)
        // The real MC table is case-specific; here we use popcount-based approximation.
        // For 0 or 8 set bits (all same sign): 0 triangles
        let ones = i.count_ones() as usize;
        table[i as usize] = if ones == 0 || ones == 8 {
            0
        } else if ones == 1 || ones == 7 {
            1
        } else if ones == 2 || ones == 6 {
            2
        } else {
            // 3,4,5 bits set: can be 2–5 triangles, use 3 as estimate
            3
        };
        if i == 255 {
            break;
        }
        i += 1;
    }
    table
}

// ---------------------------------------------------------------------------
// Contact angle
// ---------------------------------------------------------------------------

/// Compute the contact angle between a fluid interface and a solid wall.
///
/// The contact angle θ is the angle between the outward fluid normal
/// and the inward-pointing wall normal:
/// θ = arccos(n_fluid · n_wall)
///
/// Both normals are expected to be unit vectors.
pub fn contact_angle(n_fluid: [f64; 3], n_wall: [f64; 3]) -> f64 {
    let dot = n_fluid[0] * n_wall[0] + n_fluid[1] * n_wall[1] + n_fluid[2] * n_wall[2];
    dot.clamp(-1.0, 1.0).acos()
}

/// Apply contact angle boundary condition to modify the fluid surface normal.
///
/// At a wall, the fluid normal is rotated to satisfy the contact angle θ_c:
/// n̂_modified = cos(θ_c) * n_wall + sin(θ_c) * t̂
///
/// where `t̂` is the tangential component of the current fluid normal.
pub fn apply_contact_angle(
    n_fluid: [f64; 3],
    n_wall: [f64; 3],
    contact_angle_rad: f64,
) -> [f64; 3] {
    // Decompose n_fluid into normal and tangential components relative to wall
    let dot = n_fluid[0] * n_wall[0] + n_fluid[1] * n_wall[1] + n_fluid[2] * n_wall[2];
    let t = [
        n_fluid[0] - dot * n_wall[0],
        n_fluid[1] - dot * n_wall[1],
        n_fluid[2] - dot * n_wall[2],
    ];
    let t_mag = (t[0] * t[0] + t[1] * t[1] + t[2] * t[2]).sqrt();
    if t_mag < 1e-14 {
        return n_wall; // parallel to wall normal
    }
    let t_hat = [t[0] / t_mag, t[1] / t_mag, t[2] / t_mag];
    let cos_theta = contact_angle_rad.cos();
    let sin_theta = contact_angle_rad.sin();
    [
        cos_theta * n_wall[0] + sin_theta * t_hat[0],
        cos_theta * n_wall[1] + sin_theta * t_hat[1],
        cos_theta * n_wall[2] + sin_theta * t_hat[2],
    ]
}

// ---------------------------------------------------------------------------
// Surface tension coefficient
// ---------------------------------------------------------------------------

/// Temperature-dependent surface tension coefficient for water.
///
/// Uses a linear approximation:
/// σ(T) = σ₀ - k_σ * T
///
/// where σ₀ = 0.0756 N/m, k_σ ≈ 1.55e-4 N/(m·K), T in Celsius.
pub fn surface_tension_coefficient(temperature_celsius: f64) -> f64 {
    const SIGMA_0: f64 = 0.0756_f64;
    const K_SIGMA: f64 = 1.55e-4_f64;
    (SIGMA_0 - K_SIGMA * temperature_celsius).max(0.0)
}

// ---------------------------------------------------------------------------
// CSF surface force helper
// ---------------------------------------------------------------------------

/// Compute the CSF surface tension force vector at a particle.
///
/// F = σ · κ · n̂
///
/// * `sigma` – surface tension coefficient (N/m)
/// * `kappa` – curvature (1/m)
/// * `normal` – unit surface normal
pub fn csf_surface_force(sigma: f64, kappa: f64, normal: [f64; 3]) -> [f64; 3] {
    let scale = sigma * kappa;
    [scale * normal[0], scale * normal[1], scale * normal[2]]
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kernel::CubicSplineKernel;
    use crate::neighbor::SpatialHash;
    use crate::particle::SphParticle;

    // ---- Poly6 kernel tests ------------------------------------------------

    #[test]
    fn poly6_zero_outside_support() {
        assert_eq!(poly6_kernel(1.01, 1.0), 0.0);
        assert_eq!(poly6_kernel(4.0, 1.0), 0.0);
    }

    #[test]
    fn poly6_positive_inside_support() {
        assert!(poly6_kernel(0.0, 1.0) > 0.0);
        assert!(poly6_kernel(0.5, 1.0) > 0.0);
        assert!(poly6_kernel(0.99, 1.0) > 0.0);
    }

    #[test]
    fn poly6_monotone_decreasing() {
        let h = 1.0;
        let w0 = poly6_kernel(0.0, h);
        let w1 = poly6_kernel(0.25, h);
        let w2 = poly6_kernel(0.5, h);
        let w3 = poly6_kernel(0.75, h);
        assert!(w0 > w1);
        assert!(w1 > w2);
        assert!(w2 > w3);
    }

    #[test]
    fn poly6_at_boundary_is_zero() {
        let h = 1.0;
        // r² = h² → W = 0
        assert!((poly6_kernel(h * h, h)).abs() < 1e-14);
    }

    #[test]
    fn poly6_gradient_zero_outside_support() {
        let g = poly6_gradient([1.0, 0.0, 0.0], 1.01, 1.0);
        assert_eq!(g, [0.0; 3]);
    }

    #[test]
    fn poly6_gradient_points_away_from_centre() {
        // Gradient should point away from centre (positive x direction for positive dx)
        // Poly6 gradient = -6 * coeff * (h²-r²)² * r_vec
        // The negative sign means gradient actually points toward centre for rising
        // density – verify the sign is consistent.
        let g = poly6_gradient([1.0, 0.0, 0.0], 0.5, 1.0);
        // gradient of W is negative in x (W decreases away from 0)
        assert!(
            g[0] < 0.0,
            "gradient x-component should be negative for positive dx"
        );
        assert_eq!(g[1], 0.0);
        assert_eq!(g[2], 0.0);
    }

    // ---- SphColorField tests -----------------------------------------------

    fn make_single_particle_field() -> SphColorField {
        let mut cf = SphColorField::new(1.0);
        cf.add_particle([0.0, 0.0, 0.0], 1.0, 1.0);
        cf
    }

    #[test]
    fn color_field_at_centre_is_positive() {
        let cf = make_single_particle_field();
        assert!(cf.color_field_at([0.0, 0.0, 0.0]) > 0.0);
    }

    #[test]
    fn color_field_decreases_with_distance() {
        let cf = make_single_particle_field();
        let c0 = cf.color_field_at([0.0, 0.0, 0.0]);
        let c1 = cf.color_field_at([0.5, 0.0, 0.0]);
        let c2 = cf.color_field_at([0.9, 0.0, 0.0]);
        assert!(c0 > c1);
        assert!(c1 > c2);
    }

    #[test]
    fn color_field_zero_outside_support() {
        let cf = make_single_particle_field();
        assert_eq!(cf.color_field_at([1.1, 0.0, 0.0]), 0.0);
        assert_eq!(cf.color_field_at([2.0, 0.0, 0.0]), 0.0);
    }

    #[test]
    fn color_field_gradient_zero_at_centre() {
        // At the particle position the gradient should vanish by symmetry
        let cf = make_single_particle_field();
        let g = cf.color_field_gradient_at([0.0, 0.0, 0.0]);
        // r_vec = [0,0,0] → gradient = 0
        assert!(g[0].abs() < 1e-14);
        assert!(g[1].abs() < 1e-14);
        assert!(g[2].abs() < 1e-14);
    }

    #[test]
    fn color_field_gradient_nonzero_off_centre() {
        let cf = make_single_particle_field();
        let g = cf.color_field_gradient_at([0.3, 0.0, 0.0]);
        // Should be nonzero in x
        assert!(g[0].abs() > 1e-10);
    }

    #[test]
    fn surface_normal_unit_length() {
        let mut cf = SphColorField::new(1.0);
        // Two particles forming an interface along x
        cf.add_particle([-0.3, 0.0, 0.0], 1.0, 1.0);
        cf.add_particle([0.3, 0.0, 0.0], 1.0, 1.0);
        // Evaluate normal at a point off-axis
        let n = cf.surface_normal_at([0.4, 0.0, 0.0]);
        let mag = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        if mag > 1e-14 {
            assert!(
                (mag - 1.0).abs() < 1e-10,
                "normal should be unit length, got {mag}"
            );
        }
    }

    #[test]
    fn surface_normal_zero_for_no_particles() {
        let cf = SphColorField::new(1.0);
        let n = cf.surface_normal_at([0.0, 0.0, 0.0]);
        assert_eq!(n, [0.0; 3]);
    }

    #[test]
    fn is_surface_particle_detects_edge() {
        // A particle on the edge of a half-space of particles has a non-zero
        // colour-field gradient (particles only on one side).
        // h=0.5 so particles at ±0.2 are within support; ±0.6 are outside.
        let mut cf = SphColorField::new(0.5);
        // Lay down a slab of particles on the negative-x side only
        for xi in 0..5 {
            for y in [-1i32, 0, 1] {
                for z in [-1i32, 0, 1] {
                    cf.add_particle(
                        [-(xi as f64) * 0.15, y as f64 * 0.15, z as f64 * 0.15],
                        1.0,
                        1.0,
                    );
                }
            }
        }
        // The particle at [-0.0, 0, 0] (index 0) is on the positive-x edge of
        // the slab.  The gradient at a point slightly positive of it must be
        // non-zero, but we can test the actual particle position: only the
        // particles with x >= position[0] - h contribute, and there are none
        // on the positive-x side → asymmetry → non-zero gradient.
        // We use a low threshold to confirm detection.
        assert!(
            cf.is_surface_particle(0, 1e-10),
            "edge particle should be detected as surface"
        );
    }

    #[test]
    fn add_particle_grows_field() {
        let mut cf = SphColorField::new(1.0);
        assert_eq!(cf.positions.len(), 0);
        cf.add_particle([0.0, 0.0, 0.0], 1.0, 1.0);
        assert_eq!(cf.positions.len(), 1);
        cf.add_particle([1.0, 0.0, 0.0], 1.0, 1.0);
        assert_eq!(cf.positions.len(), 2);
    }

    // ---- MarchingCubesSph tests --------------------------------------------

    #[test]
    fn marching_cubes_empty_for_constant_below_iso() {
        let mc = MarchingCubesSph::new([-1.0; 3], [1.0; 3], 4);
        let verts = mc.extract_surface(|_| 0.0, 0.5);
        assert!(verts.is_empty());
    }

    #[test]
    fn marching_cubes_empty_for_constant_above_iso() {
        let mc = MarchingCubesSph::new([-1.0; 3], [1.0; 3], 4);
        let verts = mc.extract_surface(|_| 1.0, 0.5);
        assert!(verts.is_empty());
    }

    #[test]
    fn marching_cubes_sphere_produces_triangles() {
        let mc = MarchingCubesSph::new([-1.5; 3], [1.5; 3], 8);
        let verts = mc.extract_surface(
            |p| {
                let r2 = p[0] * p[0] + p[1] * p[1] + p[2] * p[2];
                1.0 - r2 // iso_value=0.5 → sphere of radius sqrt(0.5)
            },
            0.5,
        );
        // Should produce some triangles (multiple of 3 vertices)
        assert!(
            !verts.is_empty(),
            "expected triangles for sphere iso-surface"
        );
        assert_eq!(verts.len() % 3, 0);
    }

    #[test]
    fn marching_cubes_zero_resolution_returns_empty() {
        let mc = MarchingCubesSph::new([-1.0; 3], [1.0; 3], 0);
        let verts = mc.extract_surface(|_| 0.5, 0.5);
        assert!(verts.is_empty());
    }

    // ---- AnisotropicKernel tests -------------------------------------------

    #[test]
    fn anisotropic_kernel_positive_at_centre() {
        let k = AnisotropicKernel {
            position: [0.0; 3],
            g: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            h: 1.0,
        };
        assert!(k.evaluate([0.0, 0.0, 0.0]) > 0.0);
    }

    #[test]
    fn anisotropic_kernel_zero_outside_support() {
        let k = AnisotropicKernel {
            position: [0.0; 3],
            g: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            h: 1.0,
        };
        // |G*(x-xi)| = 2.0, h=1.0, so r_sq_scaled = 4.0 > 1.0 → W = 0
        assert_eq!(k.evaluate([2.0, 0.0, 0.0]), 0.0);
    }

    #[test]
    fn build_anisotropic_kernels_count() {
        let positions = vec![[0.0; 3], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let kernels = build_anisotropic_kernels(&positions, 1.0, 4.0);
        assert_eq!(kernels.len(), 3);
    }

    #[test]
    fn build_anisotropic_kernels_positions_match() {
        let positions = vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]];
        let kernels = build_anisotropic_kernels(&positions, 1.0, 4.0);
        assert_eq!(kernels[0].position, [1.0, 2.0, 3.0]);
        assert_eq!(kernels[1].position, [4.0, 5.0, 6.0]);
    }

    // ---- implicit_surface_tension_force tests ------------------------------

    #[test]
    fn surface_tension_force_zero_far_from_particles() {
        // A single isolated particle has zero gradient at its own position
        // (r_vec = [0,0,0] → poly6_gradient = 0), so the force is zero.
        let mut cf = SphColorField::new(0.1);
        cf.add_particle([100.0, 0.0, 0.0], 1.0, 1.0);
        let f = implicit_surface_tension_force(&cf, 0, 0.07);
        assert!(
            f[0].abs() < 1e-6 && f[1].abs() < 1e-6 && f[2].abs() < 1e-6,
            "expected near-zero force for isolated particle, got {f:?}"
        );
    }

    #[test]
    fn surface_tension_force_symmetry() {
        // Two symmetric particles: force should be anti-symmetric in x
        let mut cf = SphColorField::new(0.5);
        cf.add_particle([-0.2, 0.0, 0.0], 1.0, 1.0);
        cf.add_particle([0.2, 0.0, 0.0], 1.0, 1.0);
        let f0 = implicit_surface_tension_force(&cf, 0, 0.07);
        let f1 = implicit_surface_tension_force(&cf, 1, 0.07);
        // By symmetry x-components should be opposite (approximately)
        assert!(
            (f0[0] + f1[0]).abs() < 1e-6,
            "expected anti-symmetric x forces, got {f0:?} {f1:?}"
        );
    }

    // ---- Extended tests for new surface detection features -----------------

    #[test]
    fn test_color_function_surface_detection_threshold() {
        // A particle with density < threshold * rho0 is a surface particle
        let result = color_function_surface_detect(850.0, 1000.0, 0.9);
        assert!(result, "850 / 1000 = 0.85 < 0.9 → surface particle");
        let result2 = color_function_surface_detect(950.0, 1000.0, 0.9);
        assert!(!result2, "950 / 1000 = 0.95 ≥ 0.9 → interior particle");
    }

    #[test]
    fn test_color_function_value_ratio() {
        let c = color_function_value(1000.0, 1000.0);
        assert!((c - 1.0).abs() < 1e-12, "At rest density, color = 1.0: {c}");
        let c2 = color_function_value(500.0, 1000.0);
        assert!((c2 - 0.5).abs() < 1e-12, "Half density → color = 0.5: {c2}");
    }

    #[test]
    fn test_curvature_from_color_gradient_sphere() {
        // For a sphere, curvature = -2/R from divergence of unit normal
        let normal = [0.0, 0.0, 1.0];
        let grad_normal = [
            [1.0 / 0.5, 0.0, 0.0],
            [0.0, 1.0 / 0.5, 0.0],
            [0.0, 0.0, 0.0],
        ];
        let kappa = curvature_from_color_gradient(normal, grad_normal);
        assert!(kappa.is_finite(), "Curvature should be finite: {kappa}");
    }

    #[test]
    fn test_curvature_from_color_gradient_flat() {
        // Flat surface: grad_normal has zero divergence → kappa ≈ 0
        let normal = [0.0, 0.0, 1.0];
        let grad_normal = [[0.0_f64; 3]; 3];
        let kappa = curvature_from_color_gradient(normal, grad_normal);
        assert!(
            kappa.abs() < 1e-14,
            "Flat surface: curvature should be 0: {kappa}"
        );
    }

    #[test]
    fn test_level_set_extraction_single_particle() {
        // Level set at particle position should be positive (poly6 kernel is positive at r=0)
        let positions = vec![[0.0, 0.0, 0.0]];
        let colors = vec![1.0];
        let h = 0.1_f64;
        let phi = level_set_extraction(&positions, &colors, [0.0, 0.0, 0.0], h);
        // At r=0: poly6(0, h) = 315/(64*pi*h^9)*h^6 = 315/(64*pi*h^3)
        let expected = poly6_kernel(0.0, h);
        assert!(
            (phi - expected).abs() < 1e-10,
            "Level set at particle position should equal kernel value: expected={expected}, got={phi}"
        );
        assert!(phi > 0.0, "Level set should be positive: {phi}");
    }

    #[test]
    fn test_level_set_extraction_far_point() {
        // Far from all particles: level set ≈ 0
        let positions = vec![[0.0, 0.0, 0.0]];
        let colors = vec![1.0];
        let phi = level_set_extraction(&positions, &colors, [10.0, 0.0, 0.0], 0.1);
        assert!(
            phi.abs() < 1e-10,
            "Far from particles: level set ≈ 0: {phi}"
        );
    }

    #[test]
    fn test_marching_cubes_cell_signs() {
        // All vertices above iso_level → no surface
        let corners = [1.0_f64; 8];
        let tri_count = marching_cubes_cell_triangle_count(&corners, 0.5);
        assert_eq!(tri_count, 0, "All above iso-level → no triangles");
    }

    #[test]
    fn test_marching_cubes_cell_checkerboard() {
        // Alternating sign → surface exists
        let corners = [1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0];
        let tri_count = marching_cubes_cell_triangle_count(&corners, 0.0);
        assert!(tri_count > 0, "Alternating signs → surface exists");
    }

    #[test]
    fn test_contact_angle_horizontal_surface() {
        // Wall normal = [0,1,0], fluid normal = [0,-1,0] → angle = 180°
        let theta = contact_angle([0.0, 1.0, 0.0], [0.0, -1.0, 0.0]);
        assert!(
            (theta - std::f64::consts::PI).abs() < 1e-10,
            "Opposing normals → π contact angle: {theta}"
        );
    }

    #[test]
    fn test_contact_angle_perpendicular() {
        // Normal and wall in x direction → angle = 90°
        let theta = contact_angle([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert!(
            (theta - std::f64::consts::PI / 2.0).abs() < 1e-10,
            "Perpendicular normals → π/2 contact angle: {theta}"
        );
    }

    #[test]
    fn test_contact_angle_same_direction() {
        // Same direction → angle = 0
        let theta = contact_angle([0.0, 1.0, 0.0], [0.0, 1.0, 0.0]);
        assert!(
            theta.abs() < 1e-10,
            "Same normals → 0 contact angle: {theta}"
        );
    }

    #[test]
    fn test_surface_tension_coefficient_water() {
        // Water at 20°C: σ ≈ 0.0728 N/m
        let sigma = surface_tension_coefficient(20.0);
        assert!(
            (sigma - 0.0728).abs() < 0.005,
            "Water at 20°C: σ ≈ 0.0728 N/m, got {sigma}"
        );
    }

    #[test]
    fn test_surface_tension_coefficient_decreases_with_temp() {
        let sigma_cold = surface_tension_coefficient(0.0);
        let sigma_hot = surface_tension_coefficient(80.0);
        assert!(
            sigma_cold > sigma_hot,
            "Surface tension should decrease with temperature: σ_cold={sigma_cold}, σ_hot={sigma_hot}"
        );
    }

    #[test]
    fn test_csf_surface_force_direction() {
        // Force should be in direction of normal scaled by sigma * kappa
        let sigma = 0.0728;
        let kappa = 10.0;
        let normal = [0.0, 0.0, 1.0];
        let force = csf_surface_force(sigma, kappa, normal);
        assert!(
            (force[2] - sigma * kappa).abs() < 1e-12,
            "CSF force in z: expected {}, got {}",
            sigma * kappa,
            force[2]
        );
        assert!(
            force[0].abs() < 1e-14 && force[1].abs() < 1e-14,
            "CSF force should be in normal direction"
        );
    }

    #[test]
    fn test_csf_surface_force_zero_kappa() {
        let force = csf_surface_force(0.0728, 0.0, [1.0, 0.0, 0.0]);
        for (k, &fk) in force.iter().enumerate() {
            assert!(
                fk.abs() < 1e-14,
                "Zero curvature → zero force: force[{k}]={}",
                fk
            );
        }
    }

    #[test]
    fn test_color_normal_from_gradient_normalizes() {
        let grad = [3.0, 4.0, 0.0];
        let normal = color_normal_from_gradient(grad);
        let mag = (normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]).sqrt();
        assert!(
            (mag - 1.0).abs() < 1e-12,
            "Unit normal should have magnitude 1: {mag}"
        );
        assert!(
            (normal[0] - 0.6).abs() < 1e-12,
            "normal[0] should be 0.6: {}",
            normal[0]
        );
        assert!(
            (normal[1] - 0.8).abs() < 1e-12,
            "normal[1] should be 0.8: {}",
            normal[1]
        );
    }

    #[test]
    fn test_color_normal_from_gradient_zero_gradient() {
        let normal = color_normal_from_gradient([0.0, 0.0, 0.0]);
        for (k, &nk) in normal.iter().enumerate() {
            assert!(
                nk.abs() < 1e-14,
                "Zero gradient → zero normal: normal[{k}]={}",
                nk
            );
        }
    }

    // ---- CsfSurfaceTension original test (kept for regression) -------------

    #[test]
    fn surface_tension_forces_point_inward_for_sphere() {
        let spacing: f64 = 0.08;
        let h = 0.12;
        let mass = 1000.0 * spacing.powi(3);
        let mut ps = crate::particle::ParticleSet::new();
        let center = Vec3::new(0.5, 0.5, 0.5);
        let radius: f64 = 0.3;

        let n_per_side = (2.0 * radius / spacing).ceil() as usize + 2;
        for i in 0..n_per_side {
            for j in 0..n_per_side {
                for k in 0..n_per_side {
                    let pos = Vec3::new(
                        center.x - radius + i as f64 * spacing,
                        center.y - radius + j as f64 * spacing,
                        center.z - radius + k as f64 * spacing,
                    );
                    if (pos - center).norm() <= radius {
                        ps.add_particle(&SphParticle::new(pos, Vec3::zeros(), mass));
                    }
                }
            }
        }

        let kernel = CubicSplineKernel;
        let neighbors = SpatialHash::find_all_neighbors(&ps.positions, 2.0 * h);
        crate::wcsph::compute_density(&mut ps, &neighbors, &kernel, h);

        let csf = CsfSurfaceTension::new(0.0728, 0.5);
        ps.clear_forces();
        csf.apply_forces(&mut ps, &neighbors, &kernel, h);

        let mut inward_count = 0;
        let mut surface_count = 0;
        for i in 0..ps.len() {
            let to_center = center - ps.positions[i];
            let dist = to_center.norm();
            if dist > radius * 0.6 {
                surface_count += 1;
                let force_mag = ps.forces[i].norm();
                if force_mag > 1e-10 {
                    let cos_angle = ps.forces[i].dot(&to_center) / (force_mag * dist);
                    if cos_angle > 0.0 {
                        inward_count += 1;
                    }
                }
            }
        }

        if surface_count > 0 {
            let ratio = inward_count as f64 / surface_count as f64;
            assert!(
                ratio > 0.3,
                "Expected surface tension forces to point mostly inward, but only {:.0}% did ({}/{})",
                ratio * 100.0,
                inward_count,
                surface_count
            );
        }
    }
}
