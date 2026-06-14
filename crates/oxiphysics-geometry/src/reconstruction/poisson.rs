// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Screened-Poisson surface reconstruction (Kazhdan & Hoppe 2013, regular grid).
//!
//! Reconstructs a watertight triangle mesh from a set of oriented points by
//! solving the screened-Poisson equation on a uniform grid and polygonising the
//! resulting indicator field with marching cubes.

use crate::signed_distance_field::MarchingCubes;

/// An oriented point (position + unit normal).
#[derive(Debug, Clone, Copy)]
pub struct OrientedPoint {
    /// World-space position.
    pub position: [f64; 3],
    /// Unit-length surface normal.
    pub normal: [f64; 3],
}

/// Screened-Poisson reconstructed surface.
pub struct PoissonSurface {
    /// Triangle vertices (flat list of [x,y,z] per vertex).
    pub vertices: Vec<[f64; 3]>,
    /// Triangle indices (groups of 3).
    pub indices: Vec<usize>,
}

/// Errors that can occur during screened-Poisson reconstruction.
#[derive(Debug)]
pub enum PoissonError {
    /// Fewer than the minimum number of oriented points were supplied.
    TooFewPoints,
    /// The linear solver failed to converge or produced a degenerate result.
    SolverFailed(String),
    /// Marching cubes produced no surface.
    MarchingCubesFailed,
}

impl std::fmt::Display for PoissonError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            PoissonError::TooFewPoints => {
                write!(
                    f,
                    "screened-Poisson reconstruction requires at least 4 oriented points"
                )
            }
            PoissonError::SolverFailed(s) => {
                write!(f, "screened-Poisson linear solver failed: {s}")
            }
            PoissonError::MarchingCubesFailed => {
                write!(
                    f,
                    "marching cubes produced no surface from the solved field"
                )
            }
        }
    }
}

impl std::error::Error for PoissonError {}

/// Compute the 8 trilinear `(node_index, weight)` pairs for a point.
///
/// The point is mapped into grid-fractional coordinates relative to `origin`
/// with spacing `h`; the base cell index is clamped into `0..=(n-1)` so that
/// both corners (`i0` and `i0+1`) reference valid nodes in `0..=n`. The eight
/// weights sum to one (up to floating-point error).
fn corner_weights(p: &[f64; 3], origin: &[f64; 3], h: f64, n: usize) -> [(usize, f64); 8] {
    let nn = n + 1;
    let idx = |i: usize, j: usize, k: usize| -> usize { (k * nn + j) * nn + i };

    let gx = (p[0] - origin[0]) / h;
    let gy = (p[1] - origin[1]) / h;
    let gz = (p[2] - origin[2]) / h;

    let clamp_base = |g: f64| -> usize {
        let f = g.floor();
        if f < 0.0 {
            0
        } else if f as usize > n - 1 {
            n - 1
        } else {
            f as usize
        }
    };

    let i0 = clamp_base(gx);
    let j0 = clamp_base(gy);
    let k0 = clamp_base(gz);

    let frac = |g: f64, base: usize| -> f64 {
        let t = g - base as f64;
        t.clamp(0.0, 1.0)
    };

    let fx = frac(gx, i0);
    let fy = frac(gy, j0);
    let fz = frac(gz, k0);

    let mut out = [(0usize, 0.0f64); 8];
    let mut c = 0;
    for dk in 0..2 {
        let wz = if dk == 1 { fz } else { 1.0 - fz };
        for dj in 0..2 {
            let wy = if dj == 1 { fy } else { 1.0 - fy };
            for di in 0..2 {
                let wx = if di == 1 { fx } else { 1.0 - fx };
                let w = wx * wy * wz;
                out[c] = (idx(i0 + di, j0 + dj, k0 + dk), w);
                c += 1;
            }
        }
    }
    out
}

/// Reconstruct a surface from oriented points via screened Poisson.
///
/// `resolution` = number of grid cells per side (cube grid),
/// `alpha` = screening weight (default 10.0),
/// `bbox_padding` = fraction of bbox to pad on each side (e.g. 0.1 = 10%).
pub fn screened_poisson_reconstruct(
    points: &[OrientedPoint],
    resolution: usize,
    alpha: f64,
    bbox_padding: f64,
) -> Result<PoissonSurface, PoissonError> {
    // Step 1: validate.
    if points.len() < 4 {
        return Err(PoissonError::TooFewPoints);
    }

    // Step 2a: resolution guard BEFORE computing h.
    if resolution < 2 {
        return Err(PoissonError::SolverFailed("resolution too small".into()));
    }
    let n = resolution;
    let nn = n + 1;

    // Step 2b: bounding box.
    let mut min = [f64::INFINITY; 3];
    let mut max = [f64::NEG_INFINITY; 3];
    for p in points {
        for axis in 0..3 {
            if p.position[axis] < min[axis] {
                min[axis] = p.position[axis];
            }
            if p.position[axis] > max[axis] {
                max[axis] = p.position[axis];
            }
        }
    }

    let extent = [max[0] - min[0], max[1] - min[1], max[2] - min[2]];
    let l = extent[0].max(extent[1]).max(extent[2]);
    let center = [
        (min[0] + max[0]) * 0.5,
        (min[1] + max[1]) * 0.5,
        (min[2] + max[2]) * 0.5,
    ];
    let pad = l * bbox_padding;
    let side = l + 2.0 * pad;
    let h = side / n as f64;

    if !h.is_finite() || h <= 0.0 {
        return Err(PoissonError::SolverFailed("degenerate bounding box".into()));
    }

    let origin = [
        center[0] - side * 0.5,
        center[1] - side * 0.5,
        center[2] - side * 0.5,
    ];

    let n_nodes = nn * nn * nn;
    let idx = |i: usize, j: usize, k: usize| -> usize { (k * nn + j) * nn + i };

    // Step 3 + 4: splat normals into V and accumulate screening diagonal S.
    let mut vx = vec![0.0_f64; n_nodes];
    let mut vy = vec![0.0_f64; n_nodes];
    let mut vz = vec![0.0_f64; n_nodes];
    let mut s = vec![0.0_f64; n_nodes];

    for p in points {
        let corners = corner_weights(&p.position, &origin, h, n);
        for &(node, w) in corners.iter() {
            vx[node] += w * p.normal[0];
            vy[node] += w * p.normal[1];
            vz[node] += w * p.normal[2];
            s[node] += w;
        }
    }

    // Step 5: RHS b = div(V) via central differences at interior nodes.
    let mut b = vec![0.0_f64; n_nodes];
    let inv_2h = 1.0 / (2.0 * h);
    for k in 1..n {
        for j in 1..n {
            for i in 1..n {
                let dvx = (vx[idx(i + 1, j, k)] - vx[idx(i - 1, j, k)]) * inv_2h;
                let dvy = (vy[idx(i, j + 1, k)] - vy[idx(i, j - 1, k)]) * inv_2h;
                let dvz = (vz[idx(i, j, k + 1)] - vz[idx(i, j, k - 1)]) * inv_2h;
                b[idx(i, j, k)] = dvx + dvy + dvz;
            }
        }
    }

    // Step 6: assemble SPD system (-Δ + α·S) χ = b in CSR.
    let inv_h2 = 1.0 / (h * h);
    let mut off_diagonals: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n_nodes];
    let mut diagonal = vec![0.0_f64; n_nodes];

    for k in 0..nn {
        for j in 0..nn {
            for i in 0..nn {
                let id = idx(i, j, k);
                let is_boundary = i == 0 || i == n || j == 0 || j == n || k == 0 || k == n;
                if is_boundary {
                    // Dirichlet identity row pins χ = 0 on the boundary.
                    diagonal[id] = 1.0;
                    // b[id] stays 0.
                } else {
                    // Discrete negative Laplacian (SPD) plus screening.
                    diagonal[id] = 6.0 * inv_h2 + alpha * s[id];
                    let neighbors = [
                        idx(i + 1, j, k),
                        idx(i - 1, j, k),
                        idx(i, j + 1, k),
                        idx(i, j - 1, k),
                        idx(i, j, k + 1),
                        idx(i, j, k - 1),
                    ];
                    let row = &mut off_diagonals[id];
                    for &nb in neighbors.iter() {
                        row.push((nb, -inv_h2));
                    }
                }
            }
        }
    }

    let a = CsrMatrix::from_rows(n_nodes, &off_diagonals, &diagonal);

    // Step 7: solve with Jacobi-preconditioned CG.
    let tol = 1e-8;
    let max_iters = (10 * n_nodes).clamp(1000, 20000);
    let mut chi = vec![0.0_f64; n_nodes];
    let (_iters, res) = pcg_solve(&a, &b, &mut chi, tol, max_iters);

    if res.is_nan() || res > 1.0 {
        return Err(PoissonError::SolverFailed(format!(
            "CG did not converge: residual = {res}"
        )));
    }
    if chi.iter().any(|v| !v.is_finite()) {
        return Err(PoissonError::SolverFailed(
            "solution contains non-finite values".into(),
        ));
    }

    // Step 8: iso-value from χ averaged at the input points.
    let mut iso_sum = 0.0_f64;
    let mut iso_count = 0usize;
    for p in points {
        let corners = corner_weights(&p.position, &origin, h, n);
        let mut value = 0.0_f64;
        for &(node, w) in corners.iter() {
            value += w * chi[node];
        }
        if value.is_finite() {
            iso_sum += value;
            iso_count += 1;
        }
    }
    if iso_count == 0 {
        return Err(PoissonError::SolverFailed(
            "no finite iso-value evaluations".into(),
        ));
    }
    let iso = iso_sum / iso_count as f64;

    // Step 9: surface extraction via existing marching cubes.
    let bounds = [
        origin[0],
        origin[0] + side,
        origin[1],
        origin[1] + side,
        origin[2],
        origin[2] + side,
    ];
    let mc = MarchingCubes {
        nx: n,
        ny: n,
        nz: n,
        bounds,
        sdf: chi,
    };
    let result = mc.extract(iso);

    // Step 10: build PoissonSurface.
    if result.triangles.is_empty() {
        return Err(PoissonError::MarchingCubesFailed);
    }
    let vertices: Vec<[f64; 3]> = result.vertices.iter().map(|v| v.position).collect();
    let indices: Vec<usize> = result.triangles.iter().flat_map(|t| t.indices).collect();
    Ok(PoissonSurface { vertices, indices })
}

/// Minimal compressed-sparse-row (CSR) matrix for the Poisson solve (self-contained).
#[derive(Debug, Clone)]
struct CsrMatrix {
    nrows: usize,
    row_ptr: Vec<usize>,
    col_idx: Vec<usize>,
    values: Vec<f64>,
}

impl CsrMatrix {
    /// Assemble from per-row off-diagonal entries plus an explicit diagonal.
    /// Columns within each row are emitted in ascending order.
    fn from_rows(nrows: usize, off_diagonals: &[Vec<(usize, f64)>], diagonal: &[f64]) -> Self {
        let mut row_ptr = Vec::with_capacity(nrows + 1);
        let mut col_idx = Vec::new();
        let mut values = Vec::new();
        row_ptr.push(0);
        for i in 0..nrows {
            let mut entries: Vec<(usize, f64)> = Vec::with_capacity(off_diagonals[i].len() + 1);
            entries.extend_from_slice(&off_diagonals[i]);
            entries.push((i, diagonal[i]));
            entries.sort_by_key(|&(c, _)| c);
            for (c, v) in entries {
                col_idx.push(c);
                values.push(v);
            }
            row_ptr.push(col_idx.len());
        }
        Self {
            nrows,
            row_ptr,
            col_idx,
            values,
        }
    }

    /// Sparse matrix–vector product `y = A x`.
    fn matvec(&self, x: &[f64], y: &mut [f64]) {
        for (i, yi) in y.iter_mut().enumerate() {
            let mut sum = 0.0_f64;
            for k in self.row_ptr[i]..self.row_ptr[i + 1] {
                sum += self.values[k] * x[self.col_idx[k]];
            }
            *yi = sum;
        }
    }

    /// Extract the main diagonal.
    fn diag(&self) -> Vec<f64> {
        let mut d = vec![0.0_f64; self.nrows];
        for (i, di) in d.iter_mut().enumerate() {
            for k in self.row_ptr[i]..self.row_ptr[i + 1] {
                if self.col_idx[k] == i {
                    *di = self.values[k];
                    break;
                }
            }
        }
        d
    }
}

/// Dot product of two equal-length slices.
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(&x, &y)| x * y).sum()
}

/// Jacobi-preconditioned conjugate gradient for SPD `A x = b`.
/// Returns `(iterations, final_relative_residual)`. `x` is the initial guess, overwritten with the solution.
fn pcg_solve(a: &CsrMatrix, b: &[f64], x: &mut [f64], tol: f64, max_iters: usize) -> (usize, f64) {
    let n = b.len();
    let bnorm = dot(b, b).sqrt();
    if bnorm <= f64::MIN_POSITIVE {
        for xi in x.iter_mut() {
            *xi = 0.0;
        }
        return (0, 0.0);
    }
    let inv_diag: Vec<f64> = a
        .diag()
        .iter()
        .map(|&d| if d.abs() > 1e-30 { 1.0 / d } else { 0.0 })
        .collect();
    let mut r = vec![0.0_f64; n];
    a.matvec(x, &mut r);
    for (ri, bi) in r.iter_mut().zip(b.iter()) {
        *ri = *bi - *ri;
    }
    let mut res = dot(&r, &r).sqrt() / bnorm;
    if res <= tol {
        return (0, res);
    }
    let mut z: Vec<f64> = r
        .iter()
        .zip(inv_diag.iter())
        .map(|(&ri, &di)| ri * di)
        .collect();
    let mut p = z.clone();
    let mut rz_old = dot(&r, &z);
    let mut ap = vec![0.0_f64; n];
    let mut iters = 0;
    for k in 0..max_iters {
        iters = k + 1;
        a.matvec(&p, &mut ap);
        let pap = dot(&p, &ap);
        if pap.abs() < 1e-300 {
            break;
        }
        let alpha = rz_old / pap;
        for i in 0..n {
            x[i] += alpha * p[i];
            r[i] -= alpha * ap[i];
        }
        res = dot(&r, &r).sqrt() / bnorm;
        if res <= tol {
            break;
        }
        for (zi, (&ri, &di)) in z.iter_mut().zip(r.iter().zip(inv_diag.iter())) {
            *zi = ri * di;
        }
        let rz_new = dot(&r, &z);
        if rz_old.abs() < 1e-300 {
            break;
        }
        let beta = rz_new / rz_old;
        for (pi, &zi) in p.iter_mut().zip(z.iter()) {
            *pi = zi + beta * *pi;
        }
        rz_old = rz_new;
    }
    (iters, res)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_laplacian_stencil_quadratic() {
        let n = 8usize;
        let nn = n + 1;
        let origin = [-2.0_f64, -2.0, -2.0];
        let h = 4.0 / n as f64; // 0.5
        let idx = |i: usize, j: usize, k: usize| -> usize { (k * nn + j) * nn + i };
        let node_pos = |i: usize, j: usize, k: usize| -> [f64; 3] {
            [
                origin[0] + i as f64 * h,
                origin[1] + j as f64 * h,
                origin[2] + k as f64 * h,
            ]
        };
        let mut f = vec![0.0_f64; nn * nn * nn];
        for k in 0..=n {
            for j in 0..=n {
                for i in 0..=n {
                    let p = node_pos(i, j, k);
                    f[idx(i, j, k)] = p[0] * p[0] + p[1] * p[1] + p[2] * p[2];
                }
            }
        }
        let lap = (f[idx(5, 4, 4)]
            + f[idx(3, 4, 4)]
            + f[idx(4, 5, 4)]
            + f[idx(4, 3, 4)]
            + f[idx(4, 4, 5)]
            + f[idx(4, 4, 3)]
            - 6.0 * f[idx(4, 4, 4)])
            / (h * h);
        assert!((lap - 6.0).abs() < 1e-9);
    }

    #[test]
    fn test_divergence_stencil_linear() {
        let n = 8usize;
        let nn = n + 1;
        let origin = [-2.0_f64, -2.0, -2.0];
        let h = 4.0 / n as f64; // 0.5
        let idx = |i: usize, j: usize, k: usize| -> usize { (k * nn + j) * nn + i };
        let node_pos = |i: usize, j: usize, k: usize| -> [f64; 3] {
            [
                origin[0] + i as f64 * h,
                origin[1] + j as f64 * h,
                origin[2] + k as f64 * h,
            ]
        };
        let mut vx = vec![0.0_f64; nn * nn * nn];
        for k in 0..=n {
            for j in 0..=n {
                for i in 0..=n {
                    let p = node_pos(i, j, k);
                    vx[idx(i, j, k)] = p[0];
                }
            }
        }
        let div = (vx[idx(5, 4, 4)] - vx[idx(3, 4, 4)]) / (2.0 * h);
        assert!((div - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_poisson_too_few_points() {
        assert!(matches!(
            screened_poisson_reconstruct(&[], 16, 10.0, 0.1),
            Err(PoissonError::TooFewPoints)
        ));
        let pts = [
            OrientedPoint {
                position: [0.0, 0.0, 0.0],
                normal: [1.0, 0.0, 0.0],
            },
            OrientedPoint {
                position: [1.0, 0.0, 0.0],
                normal: [0.0, 1.0, 0.0],
            },
            OrientedPoint {
                position: [0.0, 1.0, 0.0],
                normal: [0.0, 0.0, 1.0],
            },
        ];
        assert!(matches!(
            screened_poisson_reconstruct(&pts, 16, 10.0, 0.1),
            Err(PoissonError::TooFewPoints)
        ));
    }

    #[test]
    fn test_poisson_sphere_smoke() {
        let n = 200usize;
        let golden = std::f64::consts::PI * (3.0 - 5.0_f64.sqrt());
        let mut pts = Vec::with_capacity(n);
        for i in 0..n {
            let y = 1.0 - 2.0 * (i as f64 + 0.5) / (n as f64);
            let r = (1.0 - y * y).max(0.0).sqrt();
            let theta = i as f64 * golden;
            let pos = [r * theta.cos(), y, r * theta.sin()];
            let norm = (pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]).sqrt();
            let normal = if norm > 0.0 {
                [pos[0] / norm, pos[1] / norm, pos[2] / norm]
            } else {
                [0.0, 1.0, 0.0]
            };
            pts.push(OrientedPoint {
                position: pos,
                normal,
            });
        }
        let result = screened_poisson_reconstruct(&pts, 24, 10.0, 0.2);
        assert!(result.is_ok());
        let Ok(surface) = result else {
            panic!("reconstruction failed")
        };
        assert!(!surface.indices.is_empty());
        let mut sum = 0.0_f64;
        let count = surface.vertices.len();
        for v in &surface.vertices {
            sum += (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        }
        assert!(count > 0);
        let mean = sum / count as f64;
        println!("sphere smoke mean radius = {mean}");
        assert!((mean - 1.0).abs() < 0.05);
    }
}
