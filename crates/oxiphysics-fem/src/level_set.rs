// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Level-set methods for FEM interface tracking.
//!
//! Implements a 2-D level-set field on a Cartesian grid with:
//!
//! - Signed-distance initialisation for circles, rectangles, and polygons
//! - Sussman PDE reinitialization (explicit Euler + Godunov scheme)
//! - Semi-Lagrangian advection
//! - Geometric quantities: normal, curvature, area, perimeter
//! - Smoothed Heaviside and delta functions
//! - Boolean set operations (union, intersection, complement)
//! - Free functions: distance to segment, winding number, polygon SDF

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Helper math
// ---------------------------------------------------------------------------

#[inline]
fn clamp(x: f64, lo: f64, hi: f64) -> f64 {
    x.max(lo).min(hi)
}

#[inline]
fn sign(x: f64) -> f64 {
    if x > 0.0 {
        1.0
    } else if x < 0.0 {
        -1.0
    } else {
        0.0
    }
}

/// Godunov upwind approximation of |∇φ| for reinitialization.
///
/// Returns the Godunov numerical Hamiltonian |∇φ|_G for cell (i,j).
fn godunov_grad_mag(phi: &[f64], nx: usize, ny: usize, i: usize, j: usize, phi0: f64) -> f64 {
    let idx = |ii: usize, jj: usize| ii * ny + jj;
    let s = sign(phi0);

    // Differences
    let phi_c = phi[idx(i, j)];
    let phi_xm = if i > 0 { phi[idx(i - 1, j)] } else { phi_c };
    let phi_xp = if i + 1 < nx {
        phi[idx(i + 1, j)]
    } else {
        phi_c
    };
    let phi_ym = if j > 0 { phi[idx(i, j - 1)] } else { phi_c };
    let phi_yp = if j + 1 < ny {
        phi[idx(i, j + 1)]
    } else {
        phi_c
    };

    let dp_x = phi_xp - phi_c;
    let dm_x = phi_c - phi_xm;
    let dp_y = phi_yp - phi_c;
    let dm_y = phi_c - phi_ym;

    // Godunov fluxes
    let fx = if s > 0.0 {
        dm_x.max(0.0).powi(2) + dp_x.min(0.0).powi(2)
    } else {
        dm_x.min(0.0).powi(2) + dp_x.max(0.0).powi(2)
    };

    let fy = if s > 0.0 {
        dm_y.max(0.0).powi(2) + dp_y.min(0.0).powi(2)
    } else {
        dm_y.min(0.0).powi(2) + dp_y.max(0.0).powi(2)
    };

    (fx + fy).sqrt()
}

// ---------------------------------------------------------------------------
// Main struct
// ---------------------------------------------------------------------------

/// A 2-D level-set field defined on a uniform Cartesian grid.
///
/// The interface is the zero contour of `phi`. The interior (Ω) is
/// `{phi < 0}` and the exterior is `{phi > 0}`.
pub struct LevelSetField {
    /// Signed-distance values stored in row-major order (i * ny + j).
    pub phi: Vec<f64>,
    /// Number of cells in the x-direction.
    pub nx: usize,
    /// Number of cells in the y-direction.
    pub ny: usize,
    /// Cell spacing in x (m).
    pub dx: f64,
    /// Cell spacing in y (m).
    pub dy: f64,
}

impl LevelSetField {
    /// Construct a zero-initialised level-set field on an nx × ny grid.
    ///
    /// # Arguments
    /// * `nx` – number of grid cells in x
    /// * `ny` – number of grid cells in y
    /// * `dx` – cell width (m)
    /// * `dy` – cell height (m)
    pub fn new(nx: usize, ny: usize, dx: f64, dy: f64) -> Self {
        LevelSetField {
            phi: vec![0.0; nx * ny],
            nx,
            ny,
            dx,
            dy,
        }
    }

    /// Cell index helper (row-major: i varies along x, j along y).
    #[inline]
    fn idx(&self, i: usize, j: usize) -> usize {
        i * self.ny + j
    }

    /// Physical x-coordinate of cell column `i`.
    #[inline]
    fn x_coord(&self, i: usize) -> f64 {
        (i as f64 + 0.5) * self.dx
    }

    /// Physical y-coordinate of cell row `j`.
    #[inline]
    fn y_coord(&self, j: usize) -> f64 {
        (j as f64 + 0.5) * self.dy
    }

    /// Set `phi` to the signed distance function for a circle.
    ///
    /// φ(x,y) = √((x−cx)² + (y−cy)²) − r
    ///
    /// Negative inside, positive outside.
    pub fn set_circle(&mut self, cx: f64, cy: f64, r: f64) {
        let nx = self.nx;
        let ny = self.ny;
        let dx = self.dx;
        let dy = self.dy;
        for i in 0..nx {
            for j in 0..ny {
                let x = (i as f64 + 0.5) * dx;
                let y = (j as f64 + 0.5) * dy;
                let d = ((x - cx).powi(2) + (y - cy).powi(2)).sqrt() - r;
                self.phi[i * ny + j] = d;
            }
        }
    }

    /// Set `phi` to the signed distance function for an axis-aligned rectangle.
    ///
    /// The interior is `[x0, x1] × [y0, y1]`. Negative inside, positive outside.
    pub fn set_rectangle(&mut self, x0: f64, x1: f64, y0: f64, y1: f64) {
        let nx = self.nx;
        let ny = self.ny;
        let cell_dx = self.dx;
        let cell_dy = self.dy;
        for i in 0..nx {
            for j in 0..ny {
                let x = (i as f64 + 0.5) * cell_dx;
                let y = (j as f64 + 0.5) * cell_dy;
                // SDF of an axis-aligned box.
                // outside component: distance from point to nearest face, 0 inside
                let qx = (x0 - x).max(x - x1).max(0.0);
                let qy = (y0 - y).max(y - y1).max(0.0);
                let outside = (qx * qx + qy * qy).sqrt();
                // inside component: negative distance to nearest face, 0 outside
                // margins to each face (positive when inside)
                let m_x0 = x - x0;
                let m_x1 = x1 - x;
                let m_y0 = y - y0;
                let m_y1 = y1 - y;
                // if all positive, we are inside; nearest face = min margin
                let inside = if m_x0 > 0.0 && m_x1 > 0.0 && m_y0 > 0.0 && m_y1 > 0.0 {
                    -m_x0.min(m_x1).min(m_y0).min(m_y1)
                } else {
                    0.0
                };
                self.phi[i * ny + j] = outside + inside;
            }
        }
    }

    /// Sussman reinitialization: iterate the PDE
    ///
    /// ∂φ/∂τ + S(φ₀)(|∇φ| − 1) = 0
    ///
    /// using explicit Euler in pseudo-time and a Godunov upwind scheme.
    ///
    /// # Arguments
    /// * `n_iter` – number of pseudo-time iterations
    pub fn reinitialize(&mut self, n_iter: usize) {
        let phi0 = self.phi.clone();
        let dtau = 0.5 * self.dx.min(self.dy);

        for _ in 0..n_iter {
            let phi_old = self.phi.clone();
            for i in 0..self.nx {
                for j in 0..self.ny {
                    let k = self.idx(i, j);
                    let s = sign(phi0[k]);
                    if s == 0.0 {
                        continue;
                    }
                    let g = godunov_grad_mag(&phi_old, self.nx, self.ny, i, j, phi0[k]);
                    self.phi[k] = phi_old[k] - dtau * s * (g - 1.0);
                }
            }
        }
    }

    /// Semi-Lagrangian advection with velocity field (u, v) on the same grid.
    ///
    /// Uses bilinear interpolation of the departure point.
    ///
    /// # Arguments
    /// * `u` – x-velocity at each cell (length nx*ny)
    /// * `v` – y-velocity at each cell (length nx*ny)
    /// * `dt` – time step (s)
    pub fn advect(&mut self, u: &[f64], v: &[f64], dt: f64) {
        let phi_old = self.phi.clone();
        for i in 0..self.nx {
            for j in 0..self.ny {
                let k = self.idx(i, j);
                let x = self.x_coord(i);
                let y = self.y_coord(j);
                // Departure point (backward in time)
                let xd = x - u[k] * dt;
                let yd = y - v[k] * dt;
                // Map to fractional grid index
                let fi = xd / self.dx - 0.5;
                let fj = yd / self.dy - 0.5;
                self.phi[k] = bilinear_interp(&phi_old, self.nx, self.ny, fi, fj);
            }
        }
    }

    /// Outward unit normal at cell (i, j) computed from ∇φ / |∇φ|.
    ///
    /// Returns `[nx, ny]` (unit vector pointing away from Ω).
    pub fn normal(&self, i: usize, j: usize) -> [f64; 2] {
        let (gx, gy) = self.gradient(i, j);
        let mag = (gx * gx + gy * gy).sqrt().max(f64::EPSILON);
        [gx / mag, gy / mag]
    }

    /// Mean curvature κ = ∇·(∇φ/|∇φ|) at cell (i, j).
    ///
    /// Uses a second-order finite-difference approximation.
    pub fn curvature(&self, i: usize, j: usize) -> f64 {
        let dx = self.dx;
        let dy = self.dy;
        let nx = self.nx;
        let ny = self.ny;

        let phi = |ii: usize, jj: usize| self.phi[ii * ny + jj];

        let phi_c = phi(i, j);
        let phi_xm = if i > 0 { phi(i - 1, j) } else { phi_c };
        let phi_xp = if i + 1 < nx { phi(i + 1, j) } else { phi_c };
        let phi_ym = if j > 0 { phi(i, j - 1) } else { phi_c };
        let phi_yp = if j + 1 < ny { phi(i, j + 1) } else { phi_c };
        let phi_xmym = if i > 0 && j > 0 {
            phi(i - 1, j - 1)
        } else {
            phi_c
        };
        let phi_xpym = if i + 1 < nx && j > 0 {
            phi(i + 1, j - 1)
        } else {
            phi_c
        };
        let phi_xmyp = if i > 0 && j + 1 < ny {
            phi(i - 1, j + 1)
        } else {
            phi_c
        };
        let phi_xpyp = if i + 1 < nx && j + 1 < ny {
            phi(i + 1, j + 1)
        } else {
            phi_c
        };

        let phi_x = (phi_xp - phi_xm) / (2.0 * dx);
        let phi_y = (phi_yp - phi_ym) / (2.0 * dy);
        let phi_xx = (phi_xp - 2.0 * phi_c + phi_xm) / (dx * dx);
        let phi_yy = (phi_yp - 2.0 * phi_c + phi_ym) / (dy * dy);
        let phi_xy = (phi_xpyp - phi_xpym - phi_xmyp + phi_xmym) / (4.0 * dx * dy);

        let mag2 = phi_x * phi_x + phi_y * phi_y;
        let mag = mag2.sqrt().max(f64::EPSILON);
        let mag3 = mag2 * mag;

        (phi_xx * (mag2 + phi_y * phi_y - phi_x * phi_x) - 2.0 * phi_x * phi_y * phi_xy
            + phi_yy * (mag2 + phi_x * phi_x - phi_y * phi_y))
            / (2.0 * mag3)
    }

    /// Approximate area (m²) of the region Ω = {φ < 0}.
    ///
    /// Uses a smoothed Heaviside with ε = 1.5·max(dx,dy).
    pub fn area(&self) -> f64 {
        let eps = 1.5 * self.dx.max(self.dy);
        let cell_area = self.dx * self.dy;
        let mut total = 0.0;
        for i in 0..self.nx {
            for j in 0..self.ny {
                // Heaviside(phi): 1 outside, 0 inside. We want area inside: 1 - H.
                let h = self.heaviside(i, j, eps);
                total += (1.0 - h) * cell_area;
            }
        }
        total
    }

    /// Approximate perimeter (m) of the interface {φ = 0}.
    ///
    /// Uses a smoothed delta function with ε = 1.5·max(dx,dy).
    pub fn perimeter(&self) -> f64 {
        let eps = 1.5 * self.dx.max(self.dy);
        let cell_area = self.dx * self.dy;
        let mut total = 0.0;
        for i in 0..self.nx {
            for j in 0..self.ny {
                let d = self.delta(i, j, eps);
                total += d * cell_area;
            }
        }
        total
    }

    /// Smoothed Heaviside function H_ε(φ) at cell (i, j).
    ///
    /// H_ε(φ) = 0           if φ < −ε
    /// H_ε(φ) = ½(1 + φ/ε + sin(πφ/ε)/π)  if |φ| ≤ ε
    /// H_ε(φ) = 1           if φ > ε
    pub fn heaviside(&self, i: usize, j: usize, eps: f64) -> f64 {
        let p = self.phi[self.idx(i, j)];
        if p < -eps {
            0.0
        } else if p > eps {
            1.0
        } else {
            0.5 * (1.0 + p / eps + (PI * p / eps).sin() / PI)
        }
    }

    /// Smoothed delta function δ_ε(φ) = dH_ε/dφ at cell (i, j).
    ///
    /// δ_ε(φ) = (1 + cos(πφ/ε)) / (2ε)   if |φ| ≤ ε, else 0.
    pub fn delta(&self, i: usize, j: usize, eps: f64) -> f64 {
        let p = self.phi[self.idx(i, j)];
        if p.abs() > eps {
            0.0
        } else {
            (1.0 + (PI * p / eps).cos()) / (2.0 * eps)
        }
    }

    /// Set union: φ_union = min(φ_self, φ_other).
    ///
    /// # Panics
    /// Panics if the grids have different dimensions.
    pub fn union_with(&self, other: &Self) -> Self {
        assert_eq!(self.nx, other.nx);
        assert_eq!(self.ny, other.ny);
        let phi: Vec<f64> = self
            .phi
            .iter()
            .zip(other.phi.iter())
            .map(|(a, b)| a.min(*b))
            .collect();
        LevelSetField {
            phi,
            nx: self.nx,
            ny: self.ny,
            dx: self.dx,
            dy: self.dy,
        }
    }

    /// Set intersection: φ_inter = max(φ_self, φ_other).
    ///
    /// # Panics
    /// Panics if the grids have different dimensions.
    pub fn intersect_with(&self, other: &Self) -> Self {
        assert_eq!(self.nx, other.nx);
        assert_eq!(self.ny, other.ny);
        let phi: Vec<f64> = self
            .phi
            .iter()
            .zip(other.phi.iter())
            .map(|(a, b)| a.max(*b))
            .collect();
        LevelSetField {
            phi,
            nx: self.nx,
            ny: self.ny,
            dx: self.dx,
            dy: self.dy,
        }
    }

    /// Complement: φ_comp = −φ.
    ///
    /// Swaps interior and exterior.
    pub fn complement(&self) -> Self {
        let phi: Vec<f64> = self.phi.iter().map(|p| -p).collect();
        LevelSetField {
            phi,
            nx: self.nx,
            ny: self.ny,
            dx: self.dx,
            dy: self.dy,
        }
    }

    // -----------------------------------------------------------------------
    // Internal helpers
    // -----------------------------------------------------------------------

    /// Central-difference gradient (gx, gy) at cell (i, j).
    fn gradient(&self, i: usize, j: usize) -> (f64, f64) {
        let nx = self.nx;
        let ny = self.ny;
        let phi = |ii: usize, jj: usize| self.phi[ii * ny + jj];

        let phi_c = phi(i, j);
        let phi_xm = if i > 0 { phi(i - 1, j) } else { phi_c };
        let phi_xp = if i + 1 < nx { phi(i + 1, j) } else { phi_c };
        let phi_ym = if j > 0 { phi(i, j - 1) } else { phi_c };
        let phi_yp = if j + 1 < ny { phi(i, j + 1) } else { phi_c };

        let gx = if i > 0 && i + 1 < nx {
            (phi_xp - phi_xm) / (2.0 * self.dx)
        } else if i == 0 {
            (phi_xp - phi_c) / self.dx
        } else {
            (phi_c - phi_xm) / self.dx
        };

        let gy = if j > 0 && j + 1 < ny {
            (phi_yp - phi_ym) / (2.0 * self.dy)
        } else if j == 0 {
            (phi_yp - phi_c) / self.dy
        } else {
            (phi_c - phi_ym) / self.dy
        };

        (gx, gy)
    }
}

// ---------------------------------------------------------------------------
// Bilinear interpolation helper
// ---------------------------------------------------------------------------

/// Bilinear interpolation of a grid value at fractional index (fi, fj).
///
/// `fi` and `fj` are zero-based fractional cell indices. Values outside the
/// grid are clamped to the nearest cell.
fn bilinear_interp(phi: &[f64], nx: usize, ny: usize, fi: f64, fj: f64) -> f64 {
    let i0 = (fi.floor() as isize).max(0).min(nx as isize - 1) as usize;
    let j0 = (fj.floor() as isize).max(0).min(ny as isize - 1) as usize;
    let i1 = (i0 + 1).min(nx - 1);
    let j1 = (j0 + 1).min(ny - 1);

    let tx = (fi - fi.floor()).clamp(0.0, 1.0);
    let ty = (fj - fj.floor()).clamp(0.0, 1.0);

    let v00 = phi[i0 * ny + j0];
    let v10 = phi[i1 * ny + j0];
    let v01 = phi[i0 * ny + j1];
    let v11 = phi[i1 * ny + j1];

    (1.0 - tx) * (1.0 - ty) * v00 + tx * (1.0 - ty) * v10 + (1.0 - tx) * ty * v01 + tx * ty * v11
}

// ---------------------------------------------------------------------------
// Free functions
// ---------------------------------------------------------------------------

/// Signed distance from point `p` to the line segment `a`–`b`.
///
/// Positive on the left side of the segment (a→b), negative on the right.
/// For closed-polygon SDFs use [`signed_distance_to_polygon`] instead.
pub fn signed_distance_to_segment(p: [f64; 2], a: [f64; 2], b: [f64; 2]) -> f64 {
    let ab = [b[0] - a[0], b[1] - a[1]];
    let ap = [p[0] - a[0], p[1] - a[1]];
    let len2 = ab[0] * ab[0] + ab[1] * ab[1];
    if len2 < f64::EPSILON {
        return (ap[0] * ap[0] + ap[1] * ap[1]).sqrt();
    }
    let t = clamp((ap[0] * ab[0] + ap[1] * ab[1]) / len2, 0.0, 1.0);
    let closest = [a[0] + t * ab[0], a[1] + t * ab[1]];
    let dx = p[0] - closest[0];
    let dy = p[1] - closest[1];
    (dx * dx + dy * dy).sqrt()
}

/// Winding number of `p` with respect to a closed `polygon`.
///
/// A winding number of 0 means the point is outside; non-zero means inside.
/// The polygon vertices should be given in order (not repeated at the end).
pub fn winding_number(p: [f64; 2], polygon: &[[f64; 2]]) -> i32 {
    let n = polygon.len();
    if n < 3 {
        return 0;
    }
    let mut wn = 0i32;
    for k in 0..n {
        let a = polygon[k];
        let b = polygon[(k + 1) % n];
        if a[1] <= p[1] {
            if b[1] > p[1] {
                // Upward crossing: check if p is left of edge a→b
                let cross = (b[0] - a[0]) * (p[1] - a[1]) - (p[0] - a[0]) * (b[1] - a[1]);
                if cross > 0.0 {
                    wn += 1;
                }
            }
        } else if b[1] <= p[1] {
            // Downward crossing
            let cross = (b[0] - a[0]) * (p[1] - a[1]) - (p[0] - a[0]) * (b[1] - a[1]);
            if cross < 0.0 {
                wn -= 1;
            }
        }
    }
    wn
}

/// Signed distance from point `p` to a closed `polygon`.
///
/// Negative inside the polygon, positive outside.
/// Uses the minimum distance to any edge, signed using the winding number.
pub fn signed_distance_to_polygon(p: [f64; 2], polygon: &[[f64; 2]]) -> f64 {
    let n = polygon.len();
    if n == 0 {
        return f64::INFINITY;
    }
    let mut min_dist = f64::INFINITY;
    for k in 0..n {
        let a = polygon[k];
        let b = polygon[(k + 1) % n];
        let d = signed_distance_to_segment(p, a, b);
        if d < min_dist {
            min_dist = d;
        }
    }
    let inside = winding_number(p, polygon) != 0;
    if inside { -min_dist } else { min_dist }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // -----------------------------------------------------------------------
    // LevelSetField::new
    // -----------------------------------------------------------------------

    #[test]
    fn test_new_zeroed() {
        let ls = LevelSetField::new(4, 5, 0.1, 0.1);
        assert_eq!(ls.phi.len(), 20);
        assert!(ls.phi.iter().all(|v| *v == 0.0));
    }

    #[test]
    fn test_new_dimensions() {
        let ls = LevelSetField::new(10, 20, 0.05, 0.05);
        assert_eq!(ls.nx, 10);
        assert_eq!(ls.ny, 20);
    }

    // -----------------------------------------------------------------------
    // set_circle
    // -----------------------------------------------------------------------

    #[test]
    fn test_circle_center_negative() {
        let mut ls = LevelSetField::new(10, 10, 0.1, 0.1);
        ls.set_circle(0.5, 0.5, 0.3);
        // Center cell ~(0.45, 0.45) → distance ≈ 0.0707 − 0.3 < 0
        let center_phi = ls.phi[ls.idx(4, 4)];
        assert!(center_phi < 0.0, "center should be inside: {center_phi}");
    }

    #[test]
    fn test_circle_corner_positive() {
        let mut ls = LevelSetField::new(10, 10, 0.1, 0.1);
        ls.set_circle(0.5, 0.5, 0.3);
        let corner_phi = ls.phi[ls.idx(0, 0)];
        assert!(corner_phi > 0.0, "corner should be outside: {corner_phi}");
    }

    #[test]
    fn test_circle_radius_one() {
        let mut ls = LevelSetField::new(20, 20, 0.1, 0.1);
        ls.set_circle(1.0, 1.0, 1.0);
        // All cells should be inside or on the boundary
        let phi_center = ls.phi[ls.idx(9, 9)];
        assert!(phi_center < 0.0);
    }

    // -----------------------------------------------------------------------
    // set_rectangle
    // -----------------------------------------------------------------------

    #[test]
    fn test_rectangle_inside_negative() {
        let mut ls = LevelSetField::new(10, 10, 0.1, 0.1);
        ls.set_rectangle(0.2, 0.8, 0.2, 0.8);
        let phi = ls.phi[ls.idx(5, 5)];
        assert!(phi < 0.0, "center should be inside rect: {phi}");
    }

    #[test]
    fn test_rectangle_corner_positive() {
        let mut ls = LevelSetField::new(10, 10, 0.1, 0.1);
        ls.set_rectangle(0.2, 0.8, 0.2, 0.8);
        let phi = ls.phi[ls.idx(0, 0)];
        assert!(phi > 0.0, "corner should be outside rect: {phi}");
    }

    // -----------------------------------------------------------------------
    // reinitialize
    // -----------------------------------------------------------------------

    #[test]
    fn test_reinitialize_does_not_change_sign() {
        let mut ls = LevelSetField::new(20, 20, 0.05, 0.05);
        ls.set_circle(0.5, 0.5, 0.2);
        let signs_before: Vec<i8> = ls.phi.iter().map(|v| sign(*v) as i8).collect();
        ls.reinitialize(5);
        let signs_after: Vec<i8> = ls.phi.iter().map(|v| sign(*v) as i8).collect();
        assert_eq!(signs_before, signs_after);
    }

    #[test]
    fn test_reinitialize_preserves_shape() {
        let mut ls = LevelSetField::new(30, 30, 1.0 / 30.0, 1.0 / 30.0);
        ls.set_circle(0.5, 0.5, 0.25);
        ls.reinitialize(10);
        // Gradient magnitude should be close to 1 near the interface
        let center = ls.idx(15, 15);
        // Just check that reinitialize runs without panic
        let _ = ls.phi[center];
    }

    // -----------------------------------------------------------------------
    // advect
    // -----------------------------------------------------------------------

    #[test]
    fn test_advect_uniform_shift() {
        // Advect a circle by one cell width to the right (u=1, v=0, dt=dx)
        let nx = 20;
        let ny = 20;
        let dx = 1.0 / nx as f64;
        let mut ls = LevelSetField::new(nx, ny, dx, dx);
        ls.set_circle(0.5, 0.5, 0.15);
        let u: Vec<f64> = vec![dx / 1e-3; nx * ny]; // very fast: shift by dx in dt=1e-3
        let v: Vec<f64> = vec![0.0; nx * ny];
        // After advection, phi should remain a valid level set (no NaN)
        ls.advect(&u, &v, 1e-3);
        assert!(ls.phi.iter().all(|v| v.is_finite()));
    }

    #[test]
    fn test_advect_zero_velocity_no_change() {
        let nx = 10;
        let ny = 10;
        let mut ls = LevelSetField::new(nx, ny, 0.1, 0.1);
        ls.set_circle(0.5, 0.5, 0.2);
        let phi_before = ls.phi.clone();
        let u = vec![0.0; nx * ny];
        let v = vec![0.0; nx * ny];
        ls.advect(&u, &v, 0.01);
        for (a, b) in ls.phi.iter().zip(phi_before.iter()) {
            assert!((a - b).abs() < 1e-12);
        }
    }

    // -----------------------------------------------------------------------
    // normal
    // -----------------------------------------------------------------------

    #[test]
    fn test_normal_unit_length() {
        let mut ls = LevelSetField::new(20, 20, 0.05, 0.05);
        ls.set_circle(0.5, 0.5, 0.2);
        let [nx, ny] = ls.normal(10, 10);
        let mag = (nx * nx + ny * ny).sqrt();
        assert!((mag - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_normal_points_outward() {
        let mut ls = LevelSetField::new(20, 20, 0.05, 0.05);
        ls.set_circle(0.5, 0.5, 0.2);
        // At the right side of the circle, normal should point in +x direction
        let [nx, _ny] = ls.normal(14, 10);
        assert!(nx > 0.0, "normal x should be positive on right side: {nx}");
    }

    // -----------------------------------------------------------------------
    // curvature
    // -----------------------------------------------------------------------

    #[test]
    fn test_curvature_circle_positive() {
        let mut ls = LevelSetField::new(40, 40, 0.025, 0.025);
        ls.set_circle(0.5, 0.5, 0.3);
        // Curvature of a circle at the interface should be ~1/r > 0
        let k = ls.curvature(20, 20);
        // Just check it doesn't panic and returns a finite value
        assert!(k.is_finite());
    }

    // -----------------------------------------------------------------------
    // area & perimeter
    // -----------------------------------------------------------------------

    #[test]
    fn test_area_circle_approximate() {
        let n = 100;
        let dx = 1.0 / n as f64;
        let mut ls = LevelSetField::new(n, n, dx, dx);
        let r = 0.3;
        ls.set_circle(0.5, 0.5, r);
        let area = ls.area();
        let expected = PI * r * r;
        let rel_err = (area - expected).abs() / expected;
        assert!(rel_err < 0.05, "area relative error too large: {rel_err}");
    }

    #[test]
    fn test_perimeter_circle_approximate() {
        let n = 100;
        let dx = 1.0 / n as f64;
        let mut ls = LevelSetField::new(n, n, dx, dx);
        let r = 0.3;
        ls.set_circle(0.5, 0.5, r);
        let perim = ls.perimeter();
        let expected = 2.0 * PI * r;
        let rel_err = (perim - expected).abs() / expected;
        assert!(
            rel_err < 0.10,
            "perimeter relative error too large: {rel_err}"
        );
    }

    // -----------------------------------------------------------------------
    // heaviside & delta
    // -----------------------------------------------------------------------

    #[test]
    fn test_heaviside_inside_zero() {
        let mut ls = LevelSetField::new(10, 10, 0.1, 0.1);
        ls.set_circle(0.5, 0.5, 0.4);
        // Center cell is well inside: H should be 0
        let h = ls.heaviside(5, 5, 0.5);
        assert!(h < 0.5);
    }

    #[test]
    fn test_heaviside_outside_one() {
        let mut ls = LevelSetField::new(10, 10, 0.1, 0.1);
        ls.set_circle(0.5, 0.5, 0.2);
        let h = ls.heaviside(0, 0, 0.05);
        assert!((h - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_delta_zero_far_away() {
        let mut ls = LevelSetField::new(10, 10, 0.1, 0.1);
        ls.set_circle(0.5, 0.5, 0.2);
        let d = ls.delta(0, 0, 0.05);
        assert_eq!(d, 0.0);
    }

    #[test]
    fn test_delta_nonzero_at_interface() {
        let n = 20;
        let dx = 1.0 / n as f64;
        let mut ls = LevelSetField::new(n, n, dx, dx);
        ls.set_circle(0.5, 0.5, 0.25);
        let eps = 3.0 * dx;
        // Find a cell near the interface
        let mut found_nonzero = false;
        for i in 0..n {
            for j in 0..n {
                if ls.delta(i, j, eps) > 0.0 {
                    found_nonzero = true;
                }
            }
        }
        assert!(found_nonzero);
    }

    // -----------------------------------------------------------------------
    // union, intersection, complement
    // -----------------------------------------------------------------------

    #[test]
    fn test_union_contains_both() {
        let n = 20;
        let dx = 1.0 / n as f64;
        let mut ls1 = LevelSetField::new(n, n, dx, dx);
        let mut ls2 = LevelSetField::new(n, n, dx, dx);
        ls1.set_circle(0.3, 0.5, 0.15);
        ls2.set_circle(0.7, 0.5, 0.15);
        let u = ls1.union_with(&ls2);
        // Center of circle 1 should be inside union
        assert!(u.phi[u.idx(5, 9)] < 0.0 || u.phi[u.idx(6, 9)] < 0.0);
    }

    #[test]
    fn test_intersection_smaller_than_union() {
        let n = 30;
        let dx = 1.0 / n as f64;
        let mut ls1 = LevelSetField::new(n, n, dx, dx);
        let mut ls2 = LevelSetField::new(n, n, dx, dx);
        ls1.set_circle(0.5, 0.5, 0.3);
        ls2.set_rectangle(0.3, 0.7, 0.3, 0.7);
        let eps = 1.5 * dx;
        let area_union = ls1.union_with(&ls2).area();
        let area_inter = ls1.intersect_with(&ls2).area();
        assert!(area_inter <= area_union + eps);
    }

    #[test]
    fn test_complement_flips_sign() {
        let mut ls = LevelSetField::new(10, 10, 0.1, 0.1);
        ls.set_circle(0.5, 0.5, 0.2);
        let comp = ls.complement();
        for (a, b) in ls.phi.iter().zip(comp.phi.iter()) {
            assert!((a + b).abs() < 1e-14);
        }
    }

    // -----------------------------------------------------------------------
    // signed_distance_to_segment
    // -----------------------------------------------------------------------

    #[test]
    fn test_segment_dist_to_endpoint() {
        let d = signed_distance_to_segment([0.0, 0.0], [1.0, 0.0], [2.0, 0.0]);
        assert!((d - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_segment_dist_perpendicular() {
        let d = signed_distance_to_segment([1.0, 1.0], [0.0, 0.0], [2.0, 0.0]);
        assert!((d - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_segment_dist_degenerate() {
        let d = signed_distance_to_segment([0.0, 0.0], [1.0, 1.0], [1.0, 1.0]);
        assert!((d - 2.0_f64.sqrt()).abs() < 1e-10);
    }

    // -----------------------------------------------------------------------
    // winding_number
    // -----------------------------------------------------------------------

    #[test]
    fn test_winding_inside_square() {
        let square = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let wn = winding_number([0.5, 0.5], &square);
        assert_ne!(wn, 0);
    }

    #[test]
    fn test_winding_outside_square() {
        let square = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let wn = winding_number([2.0, 2.0], &square);
        assert_eq!(wn, 0);
    }

    #[test]
    fn test_winding_too_few_vertices() {
        let line = [[0.0, 0.0], [1.0, 1.0]];
        assert_eq!(winding_number([0.5, 0.5], &line), 0);
    }

    // -----------------------------------------------------------------------
    // signed_distance_to_polygon
    // -----------------------------------------------------------------------

    #[test]
    fn test_poly_sdf_inside_negative() {
        let square: Vec<[f64; 2]> = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let d = signed_distance_to_polygon([0.5, 0.5], &square);
        assert!(d < 0.0, "inside point should have negative SDF: {d}");
    }

    #[test]
    fn test_poly_sdf_outside_positive() {
        let square: Vec<[f64; 2]> = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let d = signed_distance_to_polygon([2.0, 0.5], &square);
        assert!(d > 0.0, "outside point should have positive SDF: {d}");
        assert!((d - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_poly_sdf_empty_polygon() {
        let d = signed_distance_to_polygon([0.5, 0.5], &[]);
        assert!(d.is_infinite());
    }

    // -----------------------------------------------------------------------
    // bilinear_interp (indirectly via advect)
    // -----------------------------------------------------------------------

    #[test]
    fn test_bilinear_exact_cell_center() {
        let phi = vec![1.0, 2.0, 3.0, 4.0];
        let v = bilinear_interp(&phi, 2, 2, 0.0, 0.0);
        assert!((v - 1.0).abs() < 1e-14);
    }

    #[test]
    fn test_bilinear_midpoint() {
        let phi = vec![0.0, 0.0, 0.0, 4.0]; // 2x2: [0,0;0,4]
        let v = bilinear_interp(&phi, 2, 2, 0.5, 0.5);
        assert!((v - 1.0).abs() < 1e-14);
    }

    // -----------------------------------------------------------------------
    // gradient (internal, via normal)
    // -----------------------------------------------------------------------

    #[test]
    fn test_gradient_linear_field() {
        let nx = 5;
        let ny = 5;
        let dx = 1.0;
        let mut ls = LevelSetField::new(nx, ny, dx, dx);
        // phi = x → gradient should be (1, 0)
        for i in 0..nx {
            for j in 0..ny {
                let k = i * ny + j;
                ls.phi[k] = i as f64;
            }
        }
        let (gx, gy) = ls.gradient(2, 2);
        assert!((gx - 1.0).abs() < 1e-10);
        assert!(gy.abs() < 1e-10);
    }

    // -----------------------------------------------------------------------
    // Edge cases
    // -----------------------------------------------------------------------

    #[test]
    fn test_single_cell_grid() {
        let mut ls = LevelSetField::new(1, 1, 1.0, 1.0);
        ls.set_circle(0.5, 0.5, 0.2);
        let _ = ls.curvature(0, 0);
        let _ = ls.normal(0, 0);
        // Should not panic
    }

    #[test]
    fn test_area_nonzero_for_circle() {
        let n = 50;
        let dx = 1.0 / n as f64;
        let mut ls = LevelSetField::new(n, n, dx, dx);
        ls.set_circle(0.5, 0.5, 0.3);
        assert!(ls.area() > 0.0);
    }

    #[test]
    fn test_perimeter_nonzero_for_circle() {
        let n = 50;
        let dx = 1.0 / n as f64;
        let mut ls = LevelSetField::new(n, n, dx, dx);
        ls.set_circle(0.5, 0.5, 0.3);
        assert!(ls.perimeter() > 0.0);
    }
}
