// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! SDF operator functions: gradient, normal, curvature, closest-point,
//! offset, elongation, twist, bend, repeat, displace, morph, volume
//! estimation, grid-based volume, and bounding-box search.

use super::helpers::{norm3, normalize3};
use rand::RngExt;

// ─────────────────────────────────────────────────────────────────────────────
// Gradient and normal computation
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the gradient of an SDF at point `p` using central differences.
///
/// The gradient points in the direction of maximum SDF increase.
/// At the surface (SDF = 0) this is the outward surface normal.
pub fn sdf_gradient<F>(f: &F, p: [f64; 3], eps: f64) -> [f64; 3]
where
    F: Fn([f64; 3]) -> f64,
{
    let gx = f([p[0] + eps, p[1], p[2]]) - f([p[0] - eps, p[1], p[2]]);
    let gy = f([p[0], p[1] + eps, p[2]]) - f([p[0], p[1] - eps, p[2]]);
    let gz = f([p[0], p[1], p[2] + eps]) - f([p[0], p[1], p[2] - eps]);
    [gx / (2.0 * eps), gy / (2.0 * eps), gz / (2.0 * eps)]
}

/// Compute the outward unit normal to the SDF surface at point `p`.
pub fn sdf_normal<F>(f: &F, p: [f64; 3], eps: f64) -> [f64; 3]
where
    F: Fn([f64; 3]) -> f64,
{
    normalize3(sdf_gradient(f, p, eps))
}

/// Compute the approximate curvature (mean curvature) of the SDF surface at `p`.
///
/// Uses the Laplacian: κ ≈ ∇²d / 2 where d is the SDF.
pub fn sdf_mean_curvature<F>(f: &F, p: [f64; 3], eps: f64) -> f64
where
    F: Fn([f64; 3]) -> f64,
{
    let c = f(p);
    let lap = (f([p[0] + eps, p[1], p[2]])
        + f([p[0] - eps, p[1], p[2]])
        + f([p[0], p[1] + eps, p[2]])
        + f([p[0], p[1] - eps, p[2]])
        + f([p[0], p[1], p[2] + eps])
        + f([p[0], p[1], p[2] - eps])
        - 6.0 * c)
        / (eps * eps);
    lap / 2.0
}

// ─────────────────────────────────────────────────────────────────────────────
// Closest point on SDF
// ─────────────────────────────────────────────────────────────────────────────

/// Find the closest point on the zero level-set of an SDF to a query point `p`.
///
/// Uses gradient descent projection: iteratively steps toward the surface.
/// Returns `(closest_point, sdf_value_at_start)`.
pub fn closest_point_on_sdf<F>(f: &F, p: [f64; 3], max_iters: usize, eps: f64) -> ([f64; 3], f64)
where
    F: Fn([f64; 3]) -> f64,
{
    let d0 = f(p);
    let mut q = p;
    let mut d = d0;
    for _ in 0..max_iters {
        if d.abs() < eps {
            break;
        }
        let grad = sdf_gradient(f, q, 1e-4);
        let gn = norm3(grad).max(1e-30);
        // Step by the current SDF value along the gradient direction
        for k in 0..3 {
            q[k] -= d * grad[k] / gn;
        }
        d = f(q);
    }
    (q, d0)
}

/// Offset an SDF by constant `delta`: d_offset(p) = d(p) − delta.
///
/// Positive `delta` inflates; negative deflates.
#[inline]
pub fn sdf_offset(d: f64, delta: f64) -> f64 {
    d - delta
}

/// Elongate an SDF along the x-axis by amount `h` (symmetrically).
pub fn sdf_elongate_x<F>(f: &F, p: [f64; 3], h: f64) -> f64
where
    F: Fn([f64; 3]) -> f64,
{
    let q = [p[0].abs() - h, p[1], p[2]];
    let qx = q[0].min(0.0);
    let pp = [qx, q[1], q[2]];
    f(pp) + q[0].max(0.0)
}

// ─────────────────────────────────────────────────────────────────────────────
// Additional SDF utilities
// ─────────────────────────────────────────────────────────────────────────────

/// Twist an SDF: rotate the xz-plane as a function of y.
pub fn sdf_twist<F>(f: &F, p: [f64; 3], twist_k: f64) -> f64
where
    F: Fn([f64; 3]) -> f64,
{
    let c = (twist_k * p[1]).cos();
    let s = (twist_k * p[1]).sin();
    let q = [c * p[0] - s * p[2], p[1], s * p[0] + c * p[2]];
    f(q)
}

/// Bend an SDF: curve the xz-plane along the y-axis.
pub fn sdf_bend<F>(f: &F, p: [f64; 3], bend_k: f64) -> f64
where
    F: Fn([f64; 3]) -> f64,
{
    let c = (bend_k * p[0]).cos();
    let s = (bend_k * p[0]).sin();
    let q = [c * p[0] - s * p[1], s * p[0] + c * p[1], p[2]];
    f(q)
}

/// Mirror an SDF across the xz-plane (y = 0).
pub fn sdf_mirror_y<F>(f: &F, p: [f64; 3]) -> f64
where
    F: Fn([f64; 3]) -> f64,
{
    f([p[0], p[1].abs(), p[2]])
}

/// Repeat an SDF with period `c` in all three dimensions.
pub fn sdf_repeat<F>(f: &F, p: [f64; 3], c: [f64; 3]) -> f64
where
    F: Fn([f64; 3]) -> f64,
{
    let q = [
        p[0] - c[0] * (p[0] / c[0]).round(),
        p[1] - c[1] * (p[1] / c[1]).round(),
        p[2] - c[2] * (p[2] / c[2]).round(),
    ];
    f(q)
}

/// Displacement-map an SDF by adding a scalar displacement field.
pub fn sdf_displace<F, G>(base: &F, displacement: &G, p: [f64; 3]) -> f64
where
    F: Fn([f64; 3]) -> f64,
    G: Fn([f64; 3]) -> f64,
{
    base(p) + displacement(p)
}

/// Linear blend (morph) between two SDFs.
///
/// `t = 0` gives `a`, `t = 1` gives `b`.
pub fn sdf_morph<F, G>(a: &F, b: &G, p: [f64; 3], t: f64) -> f64
where
    F: Fn([f64; 3]) -> f64,
    G: Fn([f64; 3]) -> f64,
{
    (1.0 - t) * a(p) + t * b(p)
}

/// Compute the volume enclosed by an SDF on a uniform grid (Monte Carlo estimate).
///
/// `n_samples` random points are drawn from the bounding box `bounds`.
pub fn sdf_volume_estimate<F>(f: &F, bounds: [f64; 6], n_samples: usize) -> f64
where
    F: Fn([f64; 3]) -> f64,
{
    let mut rng = rand::rng();
    let lx = bounds[1] - bounds[0];
    let ly = bounds[3] - bounds[2];
    let lz = bounds[5] - bounds[4];
    let total_vol = lx * ly * lz;
    let mut inside = 0usize;
    for _ in 0..n_samples {
        let x = bounds[0] + rng.random::<f64>() * lx;
        let y = bounds[2] + rng.random::<f64>() * ly;
        let z = bounds[4] + rng.random::<f64>() * lz;
        if f([x, y, z]) <= 0.0 {
            inside += 1;
        }
    }
    total_vol * inside as f64 / n_samples as f64
}

/// Compute the volume enclosed by an SDF on a regular voxel grid.
///
/// Counts voxels whose centres are inside (SDF ≤ 0).
pub fn sdf_grid_volume<F>(f: &F, bounds: [f64; 6], nx: usize, ny: usize, nz: usize) -> f64
where
    F: Fn([f64; 3]) -> f64,
{
    let dx = (bounds[1] - bounds[0]) / nx as f64;
    let dy = (bounds[3] - bounds[2]) / ny as f64;
    let dz = (bounds[5] - bounds[4]) / nz as f64;
    let voxel_vol = dx * dy * dz;
    let mut count = 0usize;
    for iz in 0..nz {
        for iy in 0..ny {
            for ix in 0..nx {
                let p = [
                    bounds[0] + (ix as f64 + 0.5) * dx,
                    bounds[2] + (iy as f64 + 0.5) * dy,
                    bounds[4] + (iz as f64 + 0.5) * dz,
                ];
                if f(p) <= 0.0 {
                    count += 1;
                }
            }
        }
    }
    count as f64 * voxel_vol
}

/// Find the smallest axis-aligned bounding box where SDF ≤ 0.
///
/// Samples the SDF on an `n × n × n` grid within `search` and returns
/// `[xmin, xmax, ymin, ymax, zmin, zmax]`.
pub fn sdf_bounding_box<F>(f: &F, search: [f64; 6], n: usize) -> [f64; 6]
where
    F: Fn([f64; 3]) -> f64,
{
    let dx = (search[1] - search[0]) / n as f64;
    let dy = (search[3] - search[2]) / n as f64;
    let dz = (search[5] - search[4]) / n as f64;
    let mut xmin = f64::MAX;
    let mut xmax = f64::MIN;
    let mut ymin = f64::MAX;
    let mut ymax = f64::MIN;
    let mut zmin = f64::MAX;
    let mut zmax = f64::MIN;
    for iz in 0..=n {
        for iy in 0..=n {
            for ix in 0..=n {
                let p = [
                    search[0] + ix as f64 * dx,
                    search[2] + iy as f64 * dy,
                    search[4] + iz as f64 * dz,
                ];
                if f(p) <= 0.0 {
                    if p[0] < xmin {
                        xmin = p[0];
                    }
                    if p[0] > xmax {
                        xmax = p[0];
                    }
                    if p[1] < ymin {
                        ymin = p[1];
                    }
                    if p[1] > ymax {
                        ymax = p[1];
                    }
                    if p[2] < zmin {
                        zmin = p[2];
                    }
                    if p[2] > zmax {
                        zmax = p[2];
                    }
                }
            }
        }
    }
    [xmin, xmax, ymin, ymax, zmin, zmax]
}
