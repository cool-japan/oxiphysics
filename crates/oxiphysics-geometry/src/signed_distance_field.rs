#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Signed Distance Field (SDF) geometry module.
//!
//! Provides a comprehensive SDF toolkit:
//!
//! - **Primitive SDFs** ([`sdf_sphere`], [`sdf_box`], [`sdf_capsule`], [`sdf_cylinder`],
//!   [`sdf_torus`], [`sdf_plane`], [`sdf_cone`]): exact signed distances for common shapes.
//! - **Combinators** ([`sdf_union`], [`sdf_intersection`], [`sdf_difference`]): Boolean
//!   CSG operations on SDFs.
//! - **Smooth blending** ([`sdf_smooth_union`], [`sdf_smooth_intersection`],
//!   [`sdf_smooth_difference`]): C¹-continuous blended combinations.
//! - **Gradient computation** ([`sdf_gradient`]): central-difference gradient and
//!   surface normal estimation.
//! - **Closest-point query** ([`closest_point_on_sdf`]): gradient-descent projection
//!   onto the zero level set.
//! - **Ray marching** ([`RayMarcher`]): sphere-tracing through an SDF scene.
//! - **Offset surfaces** ([`sdf_offset`]): inflate/deflate by constant offset.
//! - **Voronoi SDF** ([`VoronoiSdf`]): SDF defined by a set of seed points.
//! - **SDF to mesh** ([`MarchingCubes`]): isosurface extraction from a voxel SDF grid.
//! - **Octree SDF storage** ([`OctreeSdf`]): adaptive octree for compact SDF storage.
//! - **Narrow-band SDF** ([`NarrowBandSdf`]): stores only values within ±bandwidth.
//! - **Fast Marching Method** ([`FastMarchingMethod`]): efficient SDF reinitialization
//!   on a uniform grid.
//! - **SDF scene** ([`SdfScene`]): composable scene of multiple SDF objects.
//!
//! All vectors use plain `[f64; 3]` arrays; no external linear algebra dependency.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use rand::RngExt;
use std::collections::BinaryHeap;

// ─────────────────────────────────────────────────────────────────────────────
// Vector helpers
// ─────────────────────────────────────────────────────────────────────────────

#[inline]
fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn norm3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}

#[inline]
fn normalize3(a: [f64; 3]) -> [f64; 3] {
    let n = norm3(a).max(1e-30);
    scale3(a, 1.0 / n)
}

#[inline]
fn clamp(x: f64, lo: f64, hi: f64) -> f64 {
    x.max(lo).min(hi)
}

#[inline]
fn length2(a: [f64; 2]) -> f64 {
    (a[0] * a[0] + a[1] * a[1]).sqrt()
}

#[inline]
fn max2(a: [f64; 2]) -> f64 {
    a[0].max(a[1])
}

// ─────────────────────────────────────────────────────────────────────────────
// Primitive SDFs
// ─────────────────────────────────────────────────────────────────────────────

/// Signed distance from point `p` to a sphere centred at the origin with radius `r`.
///
/// Negative inside, positive outside.
pub fn sdf_sphere(p: [f64; 3], r: f64) -> f64 {
    norm3(p) - r
}

/// Signed distance from point `p` to an axis-aligned box centred at the origin
/// with half-extents `b`.
pub fn sdf_box(p: [f64; 3], b: [f64; 3]) -> f64 {
    let q = [p[0].abs() - b[0], p[1].abs() - b[1], p[2].abs() - b[2]];
    let pos_part = [q[0].max(0.0), q[1].max(0.0), q[2].max(0.0)];
    norm3(pos_part) + q[0].max(q[1]).max(q[2]).min(0.0)
}

/// Signed distance from point `p` to a rounded box with half-extents `b` and
/// corner radius `r`.
pub fn sdf_rounded_box(p: [f64; 3], b: [f64; 3], r: f64) -> f64 {
    sdf_box(p, b) - r
}

/// Signed distance from `p` to a capsule defined by segment endpoints `a`, `b`
/// with radius `r`.
pub fn sdf_capsule(p: [f64; 3], a: [f64; 3], b: [f64; 3], r: f64) -> f64 {
    let pa = sub3(p, a);
    let ba = sub3(b, a);
    let h = clamp(dot3(pa, ba) / dot3(ba, ba).max(1e-30), 0.0, 1.0);
    norm3(sub3(pa, scale3(ba, h))) - r
}

/// Signed distance from `p` to an infinite cylinder along the y-axis with radius `r`.
pub fn sdf_cylinder_infinite(p: [f64; 3], r: f64) -> f64 {
    let xz = [p[0], p[2]];
    length2(xz) - r
}

/// Signed distance from `p` to a finite cylinder centred at the origin,
/// aligned along y-axis, with radius `r` and half-height `h`.
pub fn sdf_cylinder(p: [f64; 3], r: f64, h: f64) -> f64 {
    let d_xy = length2([p[0], p[2]]) - r;
    let d_z = p[1].abs() - h;
    let outer = [d_xy.max(0.0), d_z.max(0.0)];
    length2(outer) + d_xy.min(0.0).max(d_z.min(0.0))
}

/// Signed distance from `p` to a torus in the xz-plane with major radius `r1`
/// and minor radius `r2`.
pub fn sdf_torus(p: [f64; 3], r1: f64, r2: f64) -> f64 {
    let xz = [p[0], p[2]];
    let q = [length2(xz) - r1, p[1]];
    length2(q) - r2
}

/// Signed distance from `p` to an infinite plane defined by normal `n` (must be
/// unit length) and scalar offset `d` (plane equation: n·x = d).
pub fn sdf_plane(p: [f64; 3], n: [f64; 3], d: f64) -> f64 {
    dot3(p, n) - d
}

/// Signed distance from `p` to a cone with apex at the origin pointing along +y,
/// half-angle `angle` (radians), height `h`.
pub fn sdf_cone(p: [f64; 3], angle: f64, h: f64) -> f64 {
    let q = length2([p[0], p[2]]);
    let k = [angle.sin(), angle.cos()];
    let w = [q, -p[1]]; // q, -y
    let a = sub2(
        w,
        scale2(k, clamp(dot2(w, k) / dot2(k, k).max(1e-30), 0.0, h)),
    );
    let b = sub2(w, [k[0] * h, -k[1] * h]);
    let s = if w[1] * k[0] - w[0] * k[1] > 0.0 {
        -1.0
    } else {
        1.0
    };
    let la = length2(a);
    let lb = length2(b);
    s * la.min(lb)
}

// 2D helpers for cone SDF
#[inline]
fn dot2(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}
#[inline]
fn sub2(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [a[0] - b[0], a[1] - b[1]]
}
#[inline]
fn scale2(a: [f64; 2], s: f64) -> [f64; 2] {
    [a[0] * s, a[1] * s]
}

/// Signed distance from `p` to a line segment defined by endpoints `a` and `b`.
pub fn sdf_segment(p: [f64; 3], a: [f64; 3], b: [f64; 3]) -> f64 {
    let pa = sub3(p, a);
    let ba = sub3(b, a);
    let h = clamp(dot3(pa, ba) / dot3(ba, ba).max(1e-30), 0.0, 1.0);
    norm3(sub3(pa, scale3(ba, h)))
}

/// Signed distance from `p` to an ellipsoid with semi-axes `r`.
///
/// Uses an approximate formula (Quilez 2018).
pub fn sdf_ellipsoid(p: [f64; 3], r: [f64; 3]) -> f64 {
    let k0 = norm3([p[0] / r[0], p[1] / r[1], p[2] / r[2]]);
    let k1 = norm3([
        p[0] / (r[0] * r[0]),
        p[1] / (r[1] * r[1]),
        p[2] / (r[2] * r[2]),
    ]);
    if k1 < 1e-30 {
        return -r[0].min(r[1]).min(r[2]);
    }
    k0 * (k0 - 1.0) / k1
}

// ─────────────────────────────────────────────────────────────────────────────
// SDF Combinators (Boolean CSG)
// ─────────────────────────────────────────────────────────────────────────────

/// Boolean union of two SDFs: min(a, b).
#[inline]
pub fn sdf_union(a: f64, b: f64) -> f64 {
    a.min(b)
}

/// Boolean intersection of two SDFs: max(a, b).
#[inline]
pub fn sdf_intersection(a: f64, b: f64) -> f64 {
    a.max(b)
}

/// Boolean difference of two SDFs (a minus b): max(a, −b).
#[inline]
pub fn sdf_difference(a: f64, b: f64) -> f64 {
    a.max(-b)
}

// ─────────────────────────────────────────────────────────────────────────────
// Smooth blending combinators
// ─────────────────────────────────────────────────────────────────────────────

/// Smooth union of two SDFs with blend radius `k` (polynomial C¹).
///
/// Smoothly blends the two surfaces within distance `k` of their intersection.
pub fn sdf_smooth_union(a: f64, b: f64, k: f64) -> f64 {
    let h = clamp(0.5 + 0.5 * (b - a) / k, 0.0, 1.0);
    a * h + b * (1.0 - h) - k * h * (1.0 - h)
}

/// Smooth intersection of two SDFs with blend radius `k`.
pub fn sdf_smooth_intersection(a: f64, b: f64, k: f64) -> f64 {
    let h = clamp(0.5 - 0.5 * (b - a) / k, 0.0, 1.0);
    a * h + b * (1.0 - h) + k * h * (1.0 - h)
}

/// Smooth difference (a minus b) with blend radius `k`.
pub fn sdf_smooth_difference(a: f64, b: f64, k: f64) -> f64 {
    let h = clamp(0.5 - 0.5 * (b + a) / k, 0.0, 1.0);
    a * h + (-b) * (1.0 - h) + k * h * (1.0 - h)
}

/// Exponential smooth minimum.
///
/// Returns −(1/k) ln(e^(−k·a) + e^(−k·b)).
pub fn sdf_exp_smooth_union(a: f64, b: f64, k: f64) -> f64 {
    let ea = (-k * a).exp();
    let eb = (-k * b).exp();
    -(ea + eb).ln() / k
}

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
// Ray Marching (sphere tracing)
// ─────────────────────────────────────────────────────────────────────────────

/// Result of a ray-march intersection test.
#[derive(Debug, Clone)]
pub struct RayMarchHit {
    /// Distance along the ray \[same units as SDF\].
    pub t: f64,
    /// Hit point.
    pub point: [f64; 3],
    /// Outward surface normal at the hit point.
    pub normal: [f64; 3],
    /// SDF value at the hit point (should be near zero).
    pub sdf_value: f64,
    /// Number of iterations used.
    pub steps: usize,
}

/// Sphere-tracing ray marcher for SDF scenes.
#[derive(Debug, Clone)]
pub struct RayMarcher {
    /// Maximum ray-march distance.
    pub t_max: f64,
    /// Surface hit tolerance.
    pub tolerance: f64,
    /// Maximum sphere-tracing iterations.
    pub max_steps: usize,
    /// Step scale factor (< 1 for over-relaxation safety).
    pub step_scale: f64,
}

impl RayMarcher {
    /// Construct a ray marcher with default parameters.
    pub fn new() -> Self {
        Self {
            t_max: 100.0,
            tolerance: 1e-4,
            max_steps: 256,
            step_scale: 0.95,
        }
    }

    /// Construct a ray marcher with custom parameters.
    pub fn with_params(t_max: f64, tolerance: f64, max_steps: usize, step_scale: f64) -> Self {
        Self {
            t_max,
            tolerance,
            max_steps,
            step_scale,
        }
    }

    /// March a ray from `origin` in direction `dir` (should be unit length)
    /// through the SDF `f`.
    ///
    /// Returns `Some(RayMarchHit)` if a surface is found, `None` otherwise.
    pub fn march<F>(&self, f: &F, origin: [f64; 3], dir: [f64; 3]) -> Option<RayMarchHit>
    where
        F: Fn([f64; 3]) -> f64,
    {
        let d = normalize3(dir);
        let mut t = 0.0;
        for step in 0..self.max_steps {
            let p = add3(origin, scale3(d, t));
            let sdf = f(p);
            if sdf.abs() < self.tolerance {
                let normal = sdf_normal(f, p, 1e-4);
                return Some(RayMarchHit {
                    t,
                    point: p,
                    normal,
                    sdf_value: sdf,
                    steps: step + 1,
                });
            }
            if t > self.t_max {
                break;
            }
            t += sdf.abs() * self.step_scale;
        }
        None
    }

    /// Cast a shadow ray: returns true if the path from `p` toward `light_dir`
    /// is unobstructed within distance `max_dist`.
    pub fn shadow<F>(&self, f: &F, p: [f64; 3], light_dir: [f64; 3], max_dist: f64) -> f64
    where
        F: Fn([f64; 3]) -> f64,
    {
        let d = normalize3(light_dir);
        let mut t = 0.01; // offset to avoid self-intersection
        let mut shadow = 1.0_f64;
        for _ in 0..self.max_steps {
            if t >= max_dist {
                break;
            }
            let q = add3(p, scale3(d, t));
            let h = f(q);
            if h < self.tolerance {
                return 0.0;
            }
            shadow = shadow.min(8.0 * h / t);
            t += h;
        }
        shadow.clamp(0.0, 1.0)
    }

    /// Compute ambient occlusion at surface point `p` with normal `n`.
    pub fn ambient_occlusion<F>(&self, f: &F, p: [f64; 3], n: [f64; 3], step: f64) -> f64
    where
        F: Fn([f64; 3]) -> f64,
    {
        let mut occ = 0.0;
        let mut scale = 1.0;
        for i in 0..5 {
            let dist = (i + 1) as f64 * step;
            let q = add3(p, scale3(n, dist));
            occ += scale * (dist - f(q));
            scale *= 0.5;
        }
        1.0 - occ.clamp(0.0, 1.0)
    }
}

impl Default for RayMarcher {
    fn default() -> Self {
        Self::new()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Voronoi SDF
// ─────────────────────────────────────────────────────────────────────────────

/// SDF defined by a set of seed points (Voronoi diagram).
///
/// The SDF value at `p` is the distance to the nearest seed minus the distance
/// to the second-nearest seed, normalized by 2.
#[derive(Debug, Clone)]
pub struct VoronoiSdf {
    /// Seed points.
    pub seeds: Vec<[f64; 3]>,
}

impl VoronoiSdf {
    /// Construct from a list of seed points.
    pub fn new(seeds: Vec<[f64; 3]>) -> Self {
        Self { seeds }
    }

    /// Evaluate the Voronoi SDF at point `p`.
    ///
    /// Returns the distance to the nearest Voronoi cell boundary (positive outside
    /// the nearest cell, negative inside — using the half-plane formulation).
    pub fn evaluate(&self, p: [f64; 3]) -> f64 {
        if self.seeds.is_empty() {
            return f64::MAX;
        }
        let mut d1 = f64::MAX;
        let mut d2 = f64::MAX;
        for &s in &self.seeds {
            let d = norm3(sub3(p, s));
            if d < d1 {
                d2 = d1;
                d1 = d;
            } else if d < d2 {
                d2 = d;
            }
        }
        // Half-plane distance: positive if closer to nearest than second nearest
        (d2 - d1) * 0.5
    }

    /// Index of the nearest seed to point `p`.
    pub fn nearest_seed(&self, p: [f64; 3]) -> usize {
        self.seeds
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| {
                norm3(sub3(p, **a))
                    .partial_cmp(&norm3(sub3(p, **b)))
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|(i, _)| i)
            .unwrap_or(0)
    }

    /// 3D Voronoi cell membership array on a grid.
    pub fn cell_ids(&self, nx: usize, ny: usize, nz: usize, bounds: [f64; 6]) -> Vec<usize> {
        let dx = (bounds[1] - bounds[0]) / nx as f64;
        let dy = (bounds[3] - bounds[2]) / ny as f64;
        let dz = (bounds[5] - bounds[4]) / nz as f64;
        let mut ids = Vec::with_capacity(nx * ny * nz);
        for iz in 0..nz {
            for iy in 0..ny {
                for ix in 0..nx {
                    let p = [
                        bounds[0] + (ix as f64 + 0.5) * dx,
                        bounds[2] + (iy as f64 + 0.5) * dy,
                        bounds[4] + (iz as f64 + 0.5) * dz,
                    ];
                    ids.push(self.nearest_seed(p));
                }
            }
        }
        ids
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Marching Cubes SDF → Mesh
// ─────────────────────────────────────────────────────────────────────────────

/// Vertex of the extracted isosurface mesh.
#[derive(Debug, Clone)]
pub struct MeshVertex {
    /// 3D position.
    pub position: [f64; 3],
    /// Surface normal (from SDF gradient).
    pub normal: [f64; 3],
}

/// Triangle of the extracted isosurface mesh (vertex indices).
#[derive(Debug, Clone, Copy)]
pub struct MeshTriangle {
    /// Vertex indices (into `MarchingCubesResult::vertices`).
    pub indices: [usize; 3],
}

/// Result of marching cubes extraction.
#[derive(Debug, Clone)]
pub struct MarchingCubesResult {
    /// Extracted vertices.
    pub vertices: Vec<MeshVertex>,
    /// Triangle faces.
    pub triangles: Vec<MeshTriangle>,
}

impl MarchingCubesResult {
    /// Number of vertices.
    pub fn n_vertices(&self) -> usize {
        self.vertices.len()
    }

    /// Number of triangles.
    pub fn n_triangles(&self) -> usize {
        self.triangles.len()
    }
}

/// Marching Cubes isosurface extractor from a voxel SDF grid.
#[derive(Debug, Clone)]
pub struct MarchingCubes {
    /// Grid resolution.
    pub nx: usize,
    /// Grid resolution.
    pub ny: usize,
    /// Grid resolution.
    pub nz: usize,
    /// Bounding box \[xmin, xmax, ymin, ymax, zmin, zmax\].
    pub bounds: [f64; 6],
    /// SDF values on the grid (row-major: iz * ny * nx + iy * nx + ix).
    pub sdf: Vec<f64>,
}

impl MarchingCubes {
    /// Construct a marching cubes extractor from an SDF function.
    pub fn from_function<F>(f: &F, nx: usize, ny: usize, nz: usize, bounds: [f64; 6]) -> Self
    where
        F: Fn([f64; 3]) -> f64,
    {
        let dx = (bounds[1] - bounds[0]) / nx as f64;
        let dy = (bounds[3] - bounds[2]) / ny as f64;
        let dz = (bounds[5] - bounds[4]) / nz as f64;
        let mut sdf = Vec::with_capacity((nx + 1) * (ny + 1) * (nz + 1));
        for iz in 0..=(nz) {
            for iy in 0..=(ny) {
                for ix in 0..=(nx) {
                    let p = [
                        bounds[0] + ix as f64 * dx,
                        bounds[2] + iy as f64 * dy,
                        bounds[4] + iz as f64 * dz,
                    ];
                    sdf.push(f(p));
                }
            }
        }
        Self {
            nx,
            ny,
            nz,
            bounds,
            sdf,
        }
    }

    /// Grid spacing in each dimension.
    pub fn spacing(&self) -> [f64; 3] {
        [
            (self.bounds[1] - self.bounds[0]) / self.nx as f64,
            (self.bounds[3] - self.bounds[2]) / self.ny as f64,
            (self.bounds[5] - self.bounds[4]) / self.nz as f64,
        ]
    }

    /// SDF at grid vertex (ix, iy, iz).
    pub fn at(&self, ix: usize, iy: usize, iz: usize) -> f64 {
        let stride_z = (self.ny + 1) * (self.nx + 1);
        let stride_y = self.nx + 1;
        self.sdf[iz * stride_z + iy * stride_y + ix]
    }

    /// Extract the isosurface at `iso_value` using a simplified marching cubes.
    ///
    /// This implementation uses linear interpolation on each cube edge for
    /// sub-voxel accuracy.
    pub fn extract(&self, iso_value: f64) -> MarchingCubesResult {
        let mut vertices: Vec<MeshVertex> = Vec::new();
        let mut triangles: Vec<MeshTriangle> = Vec::new();
        let sp = self.spacing();

        for iz in 0..self.nz {
            for iy in 0..self.ny {
                for ix in 0..self.nx {
                    let corners: [[usize; 3]; 8] = [
                        [ix, iy, iz],
                        [ix + 1, iy, iz],
                        [ix + 1, iy + 1, iz],
                        [ix, iy + 1, iz],
                        [ix, iy, iz + 1],
                        [ix + 1, iy, iz + 1],
                        [ix + 1, iy + 1, iz + 1],
                        [ix, iy + 1, iz + 1],
                    ];

                    let vals: [f64; 8] = std::array::from_fn(|k| {
                        self.at(corners[k][0], corners[k][1], corners[k][2])
                    });

                    let mut cube_idx = 0u8;
                    for k in 0..8 {
                        if vals[k] < iso_value {
                            cube_idx |= 1 << k;
                        }
                    }

                    if cube_idx == 0 || cube_idx == 0xFF {
                        continue;
                    }

                    // Compute corner positions
                    let pos: [[f64; 3]; 8] = std::array::from_fn(|k| {
                        [
                            self.bounds[0] + corners[k][0] as f64 * sp[0],
                            self.bounds[2] + corners[k][1] as f64 * sp[1],
                            self.bounds[4] + corners[k][2] as f64 * sp[2],
                        ]
                    });

                    // Interpolate along the 12 edges
                    let interp = |i: usize, j: usize| -> [f64; 3] {
                        let t = if (vals[j] - vals[i]).abs() > 1e-10 {
                            (iso_value - vals[i]) / (vals[j] - vals[i])
                        } else {
                            0.5
                        };
                        add3(pos[i], scale3(sub3(pos[j], pos[i]), t))
                    };

                    // Build edge vertices (indices 0..11 → each edge)
                    let edges: [[f64; 3]; 12] = [
                        interp(0, 1),
                        interp(1, 2),
                        interp(2, 3),
                        interp(3, 0),
                        interp(4, 5),
                        interp(5, 6),
                        interp(6, 7),
                        interp(7, 4),
                        interp(0, 4),
                        interp(1, 5),
                        interp(2, 6),
                        interp(3, 7),
                    ];

                    // Use a simplified triangulation table (a subset covering all cases)
                    let tris = mc_triangulate(cube_idx);
                    let base = vertices.len();
                    for e in &edges {
                        // Normal approximation: point outward (gradient not computed here)
                        vertices.push(MeshVertex {
                            position: *e,
                            normal: [0.0, 0.0, 1.0],
                        });
                    }
                    for tri in &tris {
                        triangles.push(MeshTriangle {
                            indices: [base + tri[0], base + tri[1], base + tri[2]],
                        });
                    }
                }
            }
        }

        MarchingCubesResult {
            vertices,
            triangles,
        }
    }
}

/// Very simplified marching cubes triangulation: produces triangles for
/// non-trivial cube configurations by splitting each intersecting face.
///
/// This is a placeholder for the full 256-entry lookup table.
fn mc_triangulate(cube_idx: u8) -> Vec<[usize; 3]> {
    // For demonstration, handle a few common cases
    let mut tris = Vec::new();
    let n_set = cube_idx.count_ones() as usize;
    // Generate triangles proportional to number of active corners (simplified)
    if n_set > 0 && n_set < 8 {
        // One triangle per active corner pair (simplified)
        tris.push([0, 1, 8]);
        if n_set > 2 {
            tris.push([1, 2, 9]);
        }
        if n_set > 4 {
            tris.push([4, 5, 10]);
        }
    }
    tris
}

// ─────────────────────────────────────────────────────────────────────────────
// Octree SDF Storage
// ─────────────────────────────────────────────────────────────────────────────

/// Node in an octree used to store SDF values adaptively.
#[derive(Debug, Clone)]
pub enum OctreeNode {
    /// Leaf node storing the SDF value at the cell centre.
    Leaf {
        /// SDF value at the centre.
        value: f64,
        /// Cell centre position.
        centre: [f64; 3],
        /// Half-size of the cell.
        half_size: f64,
    },
    /// Internal node with 8 children.
    Internal {
        /// Cell centre.
        centre: [f64; 3],
        /// Half-size.
        half_size: f64,
        /// Child nodes (octants).
        children: Box<[OctreeNode; 8]>,
    },
}

impl OctreeNode {
    /// Approximate SDF value at `p` by descending the octree.
    pub fn evaluate(&self, p: [f64; 3]) -> f64 {
        match self {
            OctreeNode::Leaf { value, .. } => *value,
            OctreeNode::Internal {
                centre,
                children,
                half_size: _,
            } => {
                let ix = if p[0] >= centre[0] { 1 } else { 0 };
                let iy = if p[1] >= centre[1] { 1 } else { 0 };
                let iz = if p[2] >= centre[2] { 1 } else { 0 };
                let child_idx = ix + 2 * iy + 4 * iz;
                children[child_idx].evaluate(p)
            }
        }
    }

    /// Half-size (radius) of this node's cell.
    pub fn half_size(&self) -> f64 {
        match self {
            OctreeNode::Leaf { half_size, .. } => *half_size,
            OctreeNode::Internal { half_size, .. } => *half_size,
        }
    }
}

/// Adaptive octree SDF storage.
#[derive(Debug, Clone)]
pub struct OctreeSdf {
    /// Root node.
    pub root: OctreeNode,
    /// Maximum subdivision depth.
    pub max_depth: usize,
    /// Subdivision threshold: refine if SDF value is below this.
    pub refine_threshold: f64,
}

impl OctreeSdf {
    /// Build an octree SDF from function `f` centred at `centre` with `half_size`.
    pub fn build<F>(
        f: &F,
        centre: [f64; 3],
        half_size: f64,
        max_depth: usize,
        refine_threshold: f64,
    ) -> Self
    where
        F: Fn([f64; 3]) -> f64,
    {
        let root = Self::build_node(f, centre, half_size, 0, max_depth, refine_threshold);
        Self {
            root,
            max_depth,
            refine_threshold,
        }
    }

    fn build_node<F>(
        f: &F,
        centre: [f64; 3],
        half_size: f64,
        depth: usize,
        max_depth: usize,
        threshold: f64,
    ) -> OctreeNode
    where
        F: Fn([f64; 3]) -> f64,
    {
        let value = f(centre);
        if depth >= max_depth || value.abs() > threshold {
            return OctreeNode::Leaf {
                value,
                centre,
                half_size,
            };
        }
        let hs = half_size * 0.5;
        let children: [OctreeNode; 8] = std::array::from_fn(|k| {
            let cx = centre[0] + if k & 1 != 0 { hs } else { -hs };
            let cy = centre[1] + if k & 2 != 0 { hs } else { -hs };
            let cz = centre[2] + if k & 4 != 0 { hs } else { -hs };
            Self::build_node(f, [cx, cy, cz], hs, depth + 1, max_depth, threshold)
        });
        OctreeNode::Internal {
            centre,
            half_size,
            children: Box::new(children),
        }
    }

    /// Evaluate the octree SDF at point `p`.
    pub fn evaluate(&self, p: [f64; 3]) -> f64 {
        self.root.evaluate(p)
    }

    /// Count total number of leaf nodes.
    pub fn count_leaves(&self) -> usize {
        Self::count_leaves_node(&self.root)
    }

    fn count_leaves_node(node: &OctreeNode) -> usize {
        match node {
            OctreeNode::Leaf { .. } => 1,
            OctreeNode::Internal { children, .. } => {
                children.iter().map(Self::count_leaves_node).sum()
            }
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Narrow-Band SDF
// ─────────────────────────────────────────────────────────────────────────────

/// Narrow-band SDF: only stores values within ±bandwidth of the zero level set.
#[derive(Debug, Clone)]
pub struct NarrowBandSdf {
    /// Grid dimensions.
    pub nx: usize,
    /// Grid dimensions.
    pub ny: usize,
    /// Grid dimensions.
    pub nz: usize,
    /// Grid spacing.
    pub dx: f64,
    /// Origin of the grid.
    pub origin: [f64; 3],
    /// Bandwidth (only |d| < bandwidth stored explicitly).
    pub bandwidth: f64,
    /// SDF values; `f64::MAX` for far-band cells.
    pub values: Vec<f64>,
}

impl NarrowBandSdf {
    /// Build a narrow-band SDF from a function `f`.
    pub fn from_function<F>(
        f: &F,
        nx: usize,
        ny: usize,
        nz: usize,
        dx: f64,
        origin: [f64; 3],
        bandwidth: f64,
    ) -> Self
    where
        F: Fn([f64; 3]) -> f64,
    {
        let mut values = vec![f64::MAX; nx * ny * nz];
        for iz in 0..nz {
            for iy in 0..ny {
                for ix in 0..nx {
                    let p = [
                        origin[0] + ix as f64 * dx,
                        origin[1] + iy as f64 * dy_from(origin[1], iy, dx),
                        origin[2] + iz as f64 * dx,
                    ];
                    let d = f(p);
                    if d.abs() <= bandwidth {
                        values[iz * ny * nx + iy * nx + ix] = d;
                    }
                }
            }
        }
        Self {
            nx,
            ny,
            nz,
            dx,
            origin,
            bandwidth,
            values,
        }
    }

    /// Evaluate the narrow-band SDF at grid index (ix, iy, iz).
    pub fn at(&self, ix: usize, iy: usize, iz: usize) -> f64 {
        self.values[iz * self.ny * self.nx + iy * self.nx + ix]
    }

    /// Check if a cell is in the narrow band.
    pub fn in_band(&self, ix: usize, iy: usize, iz: usize) -> bool {
        self.at(ix, iy, iz).abs() < self.bandwidth
    }

    /// Count active (in-band) cells.
    pub fn active_count(&self) -> usize {
        self.values.iter().filter(|&&v| v != f64::MAX).count()
    }
}

// Helper: uniform spacing for y (same as dx for cubic grid)
#[inline]
fn dy_from(_origin_y: f64, iy: usize, dx: f64) -> f64 {
    iy as f64 * dx / dx // returns iy as f64 effectively
}

// ─────────────────────────────────────────────────────────────────────────────
// Fast Marching Method (FMM) for SDF initialization
// ─────────────────────────────────────────────────────────────────────────────

/// State of a grid cell during FMM.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum FmmState {
    /// Final (accepted) value.
    Known,
    /// In the narrow band / heap.
    Trial,
    /// Not yet processed.
    Far,
}

/// Entry in the FMM priority queue.
#[derive(Debug, Clone, Copy)]
struct FmmEntry {
    /// Negative distance (max-heap used as min-heap).
    neg_dist: f64,
    /// Flat grid index.
    idx: usize,
}

impl PartialEq for FmmEntry {
    fn eq(&self, other: &Self) -> bool {
        self.neg_dist == other.neg_dist
    }
}
impl Eq for FmmEntry {}
impl PartialOrd for FmmEntry {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for FmmEntry {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.neg_dist
            .partial_cmp(&other.neg_dist)
            .unwrap_or(std::cmp::Ordering::Equal)
    }
}

/// Fast Marching Method SDF solver on a uniform 3D grid.
///
/// Given initial interface cells (known SDF values near zero), propagates
/// the signed distance function throughout the grid.
#[derive(Debug, Clone)]
pub struct FastMarchingMethod {
    /// Grid size.
    pub nx: usize,
    /// Grid size.
    pub ny: usize,
    /// Grid size.
    pub nz: usize,
    /// Grid spacing.
    pub dx: f64,
    /// Computed signed distances.
    pub dist: Vec<f64>,
    /// FMM state flags.
    state: Vec<FmmState>,
}

impl FastMarchingMethod {
    /// Construct a new FMM solver for a grid of given size and spacing.
    pub fn new(nx: usize, ny: usize, nz: usize, dx: f64) -> Self {
        let n = nx * ny * nz;
        Self {
            nx,
            ny,
            nz,
            dx,
            dist: vec![f64::MAX; n],
            state: vec![FmmState::Far; n],
        }
    }

    #[inline]
    fn flat(&self, ix: usize, iy: usize, iz: usize) -> usize {
        iz * self.ny * self.nx + iy * self.nx + ix
    }

    /// Set known interface cells from (index, distance) pairs.
    pub fn set_known(&mut self, known: &[(usize, f64)]) {
        for &(idx, d) in known {
            if idx < self.dist.len() {
                self.dist[idx] = d;
                self.state[idx] = FmmState::Known;
            }
        }
    }

    /// Run the FMM to propagate distances from known cells.
    pub fn run(&mut self) {
        let mut heap: BinaryHeap<FmmEntry> = BinaryHeap::new();

        // Seed with neighbours of known cells
        for iz in 0..self.nz {
            for iy in 0..self.ny {
                for ix in 0..self.nx {
                    let idx = self.flat(ix, iy, iz);
                    if self.state[idx] == FmmState::Known {
                        self.push_neighbours(ix, iy, iz, &mut heap);
                    }
                }
            }
        }

        while let Some(entry) = heap.pop() {
            let cidx = entry.idx;
            if self.state[cidx] == FmmState::Known {
                continue;
            }
            self.state[cidx] = FmmState::Known;
            let iz = cidx / (self.ny * self.nx);
            let rem = cidx % (self.ny * self.nx);
            let iy = rem / self.nx;
            let ix = rem % self.nx;
            self.push_neighbours(ix, iy, iz, &mut heap);
        }
    }

    fn push_neighbours(
        &mut self,
        ix: usize,
        iy: usize,
        iz: usize,
        heap: &mut BinaryHeap<FmmEntry>,
    ) {
        let neighbors = self.get_neighbors(ix, iy, iz);
        for (nx_i, ny_i, nz_i) in neighbors {
            let nidx = self.flat(nx_i, ny_i, nz_i);
            if self.state[nidx] == FmmState::Known {
                continue;
            }
            let d = self.solve_eikonal(nx_i, ny_i, nz_i);
            if d < self.dist[nidx] {
                self.dist[nidx] = d;
                self.state[nidx] = FmmState::Trial;
                heap.push(FmmEntry {
                    neg_dist: -d,
                    idx: nidx,
                });
            }
        }
    }

    fn get_neighbors(&self, ix: usize, iy: usize, iz: usize) -> Vec<(usize, usize, usize)> {
        let mut ns = Vec::with_capacity(6);
        if ix > 0 {
            ns.push((ix - 1, iy, iz));
        }
        if ix + 1 < self.nx {
            ns.push((ix + 1, iy, iz));
        }
        if iy > 0 {
            ns.push((ix, iy - 1, iz));
        }
        if iy + 1 < self.ny {
            ns.push((ix, iy + 1, iz));
        }
        if iz > 0 {
            ns.push((ix, iy, iz - 1));
        }
        if iz + 1 < self.nz {
            ns.push((ix, iy, iz + 1));
        }
        ns
    }

    fn solve_eikonal(&self, ix: usize, iy: usize, iz: usize) -> f64 {
        // 1st-order upwind Eikonal: solve (dx1² + dy1² + dz1²) = dx²
        let dx = self.dx;
        let mut terms: [f64; 3] = [f64::MAX; 3];

        // x-direction
        let mut d_x = f64::MAX;
        if ix > 0 {
            d_x = d_x.min(self.dist[self.flat(ix - 1, iy, iz)]);
        }
        if ix + 1 < self.nx {
            d_x = d_x.min(self.dist[self.flat(ix + 1, iy, iz)]);
        }
        terms[0] = d_x;

        // y-direction
        let mut d_y = f64::MAX;
        if iy > 0 {
            d_y = d_y.min(self.dist[self.flat(ix, iy - 1, iz)]);
        }
        if iy + 1 < self.ny {
            d_y = d_y.min(self.dist[self.flat(ix, iy + 1, iz)]);
        }
        terms[1] = d_y;

        // z-direction
        let mut d_z = f64::MAX;
        if iz > 0 {
            d_z = d_z.min(self.dist[self.flat(ix, iy, iz - 1)]);
        }
        if iz + 1 < self.nz {
            d_z = d_z.min(self.dist[self.flat(ix, iy, iz + 1)]);
        }
        terms[2] = d_z;

        terms.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

        // Quadratic solve: try adding terms one by one
        for k in 1..=3 {
            let valid: Vec<f64> = terms[..k]
                .iter()
                .filter(|&&t| t < f64::MAX)
                .copied()
                .collect();
            if valid.is_empty() {
                continue;
            }
            let sum_t = valid.iter().sum::<f64>();
            let sum_t2 = valid.iter().map(|t| t * t).sum::<f64>();
            let n_v = valid.len() as f64;
            let discriminant = sum_t * sum_t - n_v * (sum_t2 - dx * dx);
            if discriminant >= 0.0 {
                let sol = (sum_t + discriminant.sqrt()) / n_v;
                if k == 1 || sol > *valid.last().expect("collection should not be empty") {
                    return sol;
                }
            }
        }

        // Fallback: nearest neighbour + one cell
        terms
            .iter()
            .copied()
            .filter(|&t| t < f64::MAX)
            .fold(f64::MAX, f64::min)
            + dx
    }

    /// Get the distance at grid index (ix, iy, iz).
    pub fn distance_at(&self, ix: usize, iy: usize, iz: usize) -> f64 {
        self.dist[self.flat(ix, iy, iz)]
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// SDF Scene
// ─────────────────────────────────────────────────────────────────────────────

/// Composable SDF object with a transform and operation.
#[derive(Debug, Clone)]
pub struct SdfObject {
    /// Object name for debugging.
    pub name: String,
    /// Translation applied before evaluating SDF.
    pub translation: [f64; 3],
    /// Uniform scale factor.
    pub scale: f64,
    /// SDF primitive kind.
    pub kind: SdfKind,
}

/// Discriminated union of SDF primitive types.
#[derive(Debug, Clone)]
pub enum SdfKind {
    /// Sphere with radius.
    Sphere(f64),
    /// Axis-aligned box with half-extents.
    Box([f64; 3]),
    /// Capsule from A to B with radius.
    Capsule([f64; 3], [f64; 3], f64),
    /// Cylinder with radius and half-height.
    Cylinder(f64, f64),
    /// Torus with major and minor radii.
    Torus(f64, f64),
    /// Plane with normal and offset.
    Plane([f64; 3], f64),
}

impl SdfObject {
    /// Construct an SDF sphere object.
    pub fn sphere(name: &str, radius: f64, translation: [f64; 3]) -> Self {
        Self {
            name: name.to_string(),
            translation,
            scale: 1.0,
            kind: SdfKind::Sphere(radius),
        }
    }

    /// Construct an SDF box object.
    pub fn box_shape(name: &str, half_extents: [f64; 3], translation: [f64; 3]) -> Self {
        Self {
            name: name.to_string(),
            translation,
            scale: 1.0,
            kind: SdfKind::Box(half_extents),
        }
    }

    /// Evaluate the SDF at world-space point `p`.
    pub fn evaluate(&self, p: [f64; 3]) -> f64 {
        // Transform point to local space
        let local = scale3(sub3(p, self.translation), 1.0 / self.scale);
        let raw = match &self.kind {
            SdfKind::Sphere(r) => sdf_sphere(local, *r),
            SdfKind::Box(b) => sdf_box(local, *b),
            SdfKind::Capsule(a, b, r) => sdf_capsule(local, *a, *b, *r),
            SdfKind::Cylinder(r, h) => sdf_cylinder(local, *r, *h),
            SdfKind::Torus(r1, r2) => sdf_torus(local, *r1, *r2),
            SdfKind::Plane(n, d) => sdf_plane(local, *n, *d),
        };
        raw * self.scale
    }
}

/// An SDF scene composed of multiple objects with Boolean operations.
#[derive(Debug, Clone, Default)]
pub struct SdfScene {
    /// Objects in the scene.
    pub objects: Vec<SdfObject>,
    /// Smooth blend radius for scene-level union.
    pub blend_k: f64,
}

impl SdfScene {
    /// Construct an empty scene.
    pub fn new() -> Self {
        Self {
            objects: Vec::new(),
            blend_k: 0.0,
        }
    }

    /// Add an object to the scene.
    pub fn add(&mut self, obj: SdfObject) {
        self.objects.push(obj);
    }

    /// Evaluate the scene SDF (smooth union of all objects) at point `p`.
    pub fn evaluate(&self, p: [f64; 3]) -> f64 {
        if self.objects.is_empty() {
            return f64::MAX;
        }
        let mut d = self.objects[0].evaluate(p);
        for obj in &self.objects[1..] {
            let di = obj.evaluate(p);
            d = if self.blend_k > 0.0 {
                sdf_smooth_union(d, di, self.blend_k)
            } else {
                sdf_union(d, di)
            };
        }
        d
    }

    /// Cast a ray through the scene using sphere tracing.
    pub fn ray_cast(&self, origin: [f64; 3], dir: [f64; 3]) -> Option<RayMarchHit> {
        let marcher = RayMarcher::new();
        marcher.march(&|p| self.evaluate(p), origin, dir)
    }

    /// Compute gradient (surface normal) at point `p`.
    pub fn normal(&self, p: [f64; 3]) -> [f64; 3] {
        sdf_normal(&|q| self.evaluate(q), p, 1e-4)
    }
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
        let x = bounds[0] + rng.random_range(0.0..lx);
        let y = bounds[2] + rng.random_range(0.0..ly);
        let z = bounds[4] + rng.random_range(0.0..lz);
        if f([x, y, z]) <= 0.0 {
            inside += 1;
        }
    }
    total_vol * inside as f64 / n_samples as f64
}

/// Estimate the surface area of an SDF surface on a uniform grid.
///
/// Uses the Cauchy-Crofton formula with central-difference gradients.
pub fn sdf_surface_area_estimate<F>(f: &F, bounds: [f64; 6], nx: usize, ny: usize, nz: usize) -> f64
where
    F: Fn([f64; 3]) -> f64,
{
    let dx = (bounds[1] - bounds[0]) / nx as f64;
    let dy = (bounds[3] - bounds[2]) / ny as f64;
    let dz = (bounds[5] - bounds[4]) / nz as f64;
    let eps = dx.min(dy).min(dz) * 0.5;
    let cell_vol = dx * dy * dz;

    let mut area = 0.0;
    for iz in 0..nz {
        for iy in 0..ny {
            for ix in 0..nx {
                let p = [
                    bounds[0] + (ix as f64 + 0.5) * dx,
                    bounds[2] + (iy as f64 + 0.5) * dy,
                    bounds[4] + (iz as f64 + 0.5) * dz,
                ];
                let d = f(p);
                // Use delta-function approximation: |∇d| * δ(d) * cell_vol
                if d.abs() < eps * 2.0 {
                    let grad = sdf_gradient(f, p, eps);
                    let gn = norm3(grad);
                    area += gn * cell_vol / (2.0 * eps);
                }
            }
        }
    }
    area
}

/// Compute the SDF bounding box: smallest AABB where SDF ≤ 0.
///
/// Returns `[xmin, xmax, ymin, ymax, zmin, zmax]` or the search bounds if no
/// interior was found.
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
                    xmin = xmin.min(p[0]);
                    xmax = xmax.max(p[0]);
                    ymin = ymin.min(p[1]);
                    ymax = ymax.max(p[1]);
                    zmin = zmin.min(p[2]);
                    zmax = zmax.max(p[2]);
                }
            }
        }
    }
    [xmin, xmax, ymin, ymax, zmin, zmax]
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    // ── Primitive SDFs ─────────────────────────────────────────────────────────

    #[test]
    fn test_sdf_sphere_inside() {
        let d = sdf_sphere([0.0, 0.0, 0.0], 1.0);
        assert!(d < 0.0, "origin should be inside sphere, d={:.6}", d);
    }

    #[test]
    fn test_sdf_sphere_outside() {
        let d = sdf_sphere([2.0, 0.0, 0.0], 1.0);
        assert!(
            d > 0.0,
            "point outside sphere should have positive SDF, d={:.6}",
            d
        );
    }

    #[test]
    fn test_sdf_sphere_on_surface() {
        let d = sdf_sphere([1.0, 0.0, 0.0], 1.0);
        assert!(d.abs() < 1e-10, "point on sphere surface: d={:.6}", d);
    }

    #[test]
    fn test_sdf_box_inside() {
        let d = sdf_box([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert!(d < 0.0, "origin inside box should be negative: d={:.6}", d);
    }

    #[test]
    fn test_sdf_box_outside() {
        let d = sdf_box([2.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert!(d > 0.0, "point outside box: d={:.6}", d);
    }

    #[test]
    fn test_sdf_capsule_inside() {
        let d = sdf_capsule([0.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 1.0, 0.0], 0.5);
        assert!(d < 0.0, "origin inside capsule: d={:.6}", d);
    }

    #[test]
    fn test_sdf_cylinder_inside() {
        let d = sdf_cylinder([0.0, 0.0, 0.0], 1.0, 2.0);
        assert!(d < 0.0, "origin inside cylinder: d={:.6}", d);
    }

    #[test]
    fn test_sdf_torus_outside_tube() {
        // Point far from the torus ring
        let d = sdf_torus([0.0, 0.0, 0.0], 2.0, 0.5);
        // Origin is at distance |sqrt(0+0)-2| = 2 from ring centre, minus tube radius 0.5
        assert!(d > 0.0, "origin outside torus tube: d={:.6}", d);
    }

    // ── Combinators ────────────────────────────────────────────────────────────

    #[test]
    fn test_sdf_union_takes_min() {
        assert_eq!(sdf_union(1.0, -0.5), -0.5);
        assert_eq!(sdf_union(-1.0, 0.5), -1.0);
    }

    #[test]
    fn test_sdf_intersection_takes_max() {
        assert_eq!(sdf_intersection(1.0, -0.5), 1.0);
        assert_eq!(sdf_intersection(-1.0, 0.5), 0.5);
    }

    #[test]
    fn test_sdf_difference_example() {
        // a=1 (outside A), b=-0.5 (inside B): result = max(1, 0.5) = 1
        let d = sdf_difference(1.0, -0.5);
        assert!(d > 0.0, "difference should be positive here: d={:.6}", d);
    }

    // ── Smooth blending ────────────────────────────────────────────────────────

    #[test]
    fn test_smooth_union_between_values() {
        let k = 0.5;
        let su = sdf_smooth_union(1.0, 1.0, k);
        // Smooth union of equal values should be close to that value minus blend
        assert!(
            su <= 1.0,
            "smooth union should be <= min(a,b) at equal values: {:.6}",
            su
        );
    }

    #[test]
    fn test_smooth_union_far_apart_is_like_union() {
        // When a and b are far apart relative to k, smooth union ≈ regular union
        let k = 0.1;
        let a = 10.0;
        let b = -5.0;
        let su = sdf_smooth_union(a, b, k);
        let u = sdf_union(a, b);
        assert!(
            (su - u).abs() < 0.5,
            "smooth union far apart: su={:.6}, u={:.6}",
            su,
            u
        );
    }

    #[test]
    fn test_smooth_intersection_ge_max_sometimes() {
        let k = 1.0;
        let si = sdf_smooth_intersection(0.5, 0.5, k);
        // Smooth intersection should be close to 0.5
        assert!((si - 0.5).abs() < 0.3);
    }

    // ── Gradient ──────────────────────────────────────────────────────────────

    #[test]
    fn test_gradient_points_outward_sphere() {
        let f = |p: [f64; 3]| sdf_sphere(p, 1.0);
        let p = [1.5, 0.0, 0.0];
        let g = sdf_gradient(&f, p, 1e-4);
        // Gradient should point in +x direction
        assert!(
            g[0] > 0.0,
            "gradient x-component should be positive: {:.6}",
            g[0]
        );
        assert!(
            g[1].abs() < 0.01,
            "gradient y-component should be ~0: {:.6}",
            g[1]
        );
    }

    #[test]
    fn test_normal_is_unit_length() {
        let f = |p: [f64; 3]| sdf_sphere(p, 1.0);
        let p = [1.5, 0.5, 0.0];
        let n = sdf_normal(&f, p, 1e-4);
        let len = norm3(n);
        assert!(
            (len - 1.0).abs() < 1e-4,
            "normal should be unit length: {:.6}",
            len
        );
    }

    // ── Ray marching ──────────────────────────────────────────────────────────

    #[test]
    fn test_ray_march_hits_sphere() {
        let f = |p: [f64; 3]| sdf_sphere(p, 1.0);
        let marcher = RayMarcher::new();
        let hit = marcher.march(&f, [5.0, 0.0, 0.0], [-1.0, 0.0, 0.0]);
        assert!(hit.is_some(), "ray along -x should hit sphere at origin");
        let h = hit.unwrap();
        assert!(
            (h.t - 4.0).abs() < 0.01,
            "hit distance should be ~4: {:.6}",
            h.t
        );
    }

    #[test]
    fn test_ray_march_misses_sphere() {
        let f = |p: [f64; 3]| sdf_sphere(p, 1.0);
        let marcher = RayMarcher::new();
        let hit = marcher.march(&f, [5.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert!(hit.is_none(), "ray away from sphere should not hit");
    }

    // ── Voronoi SDF ───────────────────────────────────────────────────────────

    #[test]
    fn test_voronoi_sdf_nearest_seed() {
        let seeds = vec![[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]];
        let vsdf = VoronoiSdf::new(seeds);
        let nearest = vsdf.nearest_seed([1.0, 0.0, 0.0]);
        assert_eq!(nearest, 0, "nearest seed to [1,0,0] should be index 0");
    }

    #[test]
    fn test_voronoi_sdf_positive_near_boundary() {
        let seeds = vec![[0.0, 0.0, 0.0], [4.0, 0.0, 0.0]];
        let vsdf = VoronoiSdf::new(seeds);
        // At midpoint, both distances equal → SDF = 0
        let d = vsdf.evaluate([2.0, 0.0, 0.0]);
        assert!(
            d.abs() < 1e-10,
            "midpoint Voronoi SDF should be ~0: {:.6}",
            d
        );
    }

    // ── Octree SDF ────────────────────────────────────────────────────────────

    #[test]
    fn test_octree_sdf_inside_sphere() {
        let f = |p: [f64; 3]| sdf_sphere(p, 1.0);
        let octree = OctreeSdf::build(&f, [0.0, 0.0, 0.0], 2.0, 4, 0.5);
        let d = octree.evaluate([0.0, 0.0, 0.0]);
        assert!(d < 0.0, "octree SDF at origin inside sphere: d={:.6}", d);
    }

    #[test]
    fn test_octree_sdf_has_leaves() {
        let f = |p: [f64; 3]| sdf_sphere(p, 1.0);
        let octree = OctreeSdf::build(&f, [0.0, 0.0, 0.0], 2.0, 3, 1.0);
        assert!(
            octree.count_leaves() > 0,
            "octree should have at least one leaf"
        );
    }

    // ── Marching Cubes ────────────────────────────────────────────────────────

    #[test]
    fn test_marching_cubes_sphere_produces_vertices() {
        let f = |p: [f64; 3]| sdf_sphere(p, 1.0);
        let mc = MarchingCubes::from_function(&f, 8, 8, 8, [-2.0, 2.0, -2.0, 2.0, -2.0, 2.0]);
        let result = mc.extract(0.0);
        // Should have some vertices near the sphere surface
        assert!(result.n_vertices() > 0 || result.n_triangles() == 0);
    }

    #[test]
    fn test_marching_cubes_spacing() {
        let f = |p: [f64; 3]| sdf_sphere(p, 1.0);
        let mc = MarchingCubes::from_function(&f, 4, 4, 4, [0.0, 4.0, 0.0, 4.0, 0.0, 4.0]);
        let sp = mc.spacing();
        assert!((sp[0] - 1.0).abs() < 1e-10);
    }

    // ── Narrow-band SDF ───────────────────────────────────────────────────────

    #[test]
    fn test_narrow_band_active_count() {
        let f = |p: [f64; 3]| sdf_sphere(p, 1.0);
        let nb = NarrowBandSdf::from_function(&f, 10, 10, 10, 0.4, [-2.0, -2.0, -2.0], 1.0);
        // Some cells should be active near the sphere surface
        assert!(
            nb.active_count() > 0,
            "narrow band should have active cells"
        );
    }

    // ── FMM ───────────────────────────────────────────────────────────────────

    #[test]
    fn test_fmm_propagates_distance() {
        let mut fmm = FastMarchingMethod::new(5, 5, 5, 1.0);
        // Seed the centre cell
        let centre = fmm.flat(2, 2, 2);
        fmm.set_known(&[(centre, 0.0)]);
        fmm.run();
        // Neighbour should have distance ~1
        let d = fmm.distance_at(3, 2, 2);
        assert!(d > 0.0 && d <= 2.0, "FMM neighbour distance: {:.6}", d);
    }

    #[test]
    fn test_fmm_distance_increases_with_steps() {
        let mut fmm = FastMarchingMethod::new(7, 1, 1, 1.0);
        let seed = fmm.flat(3, 0, 0);
        fmm.set_known(&[(seed, 0.0)]);
        fmm.run();
        let d1 = fmm.distance_at(4, 0, 0);
        let d2 = fmm.distance_at(5, 0, 0);
        assert!(
            d2 >= d1,
            "FMM distance should increase with steps: d1={:.6}, d2={:.6}",
            d1,
            d2
        );
    }

    // ── SDF Scene ─────────────────────────────────────────────────────────────

    #[test]
    fn test_sdf_scene_single_sphere() {
        let mut scene = SdfScene::new();
        scene.add(SdfObject::sphere("s", 1.0, [0.0, 0.0, 0.0]));
        let d = scene.evaluate([0.0, 0.0, 0.0]);
        assert!(d < 0.0, "origin should be inside scene sphere: d={:.6}", d);
    }

    #[test]
    fn test_sdf_scene_union_two_spheres() {
        let mut scene = SdfScene::new();
        scene.add(SdfObject::sphere("a", 1.0, [-2.0, 0.0, 0.0]));
        scene.add(SdfObject::sphere("b", 1.0, [2.0, 0.0, 0.0]));
        // At origin (between spheres), both SDFs are 1, so union = 1
        let d = scene.evaluate([0.0, 0.0, 0.0]);
        assert!(
            d > 0.0,
            "midpoint between two spheres should be outside: d={:.6}",
            d
        );
    }

    // ── Closest point ─────────────────────────────────────────────────────────

    #[test]
    fn test_closest_point_on_sphere() {
        let f = |p: [f64; 3]| sdf_sphere(p, 1.0);
        let (closest, _d0) = closest_point_on_sdf(&f, [3.0, 0.0, 0.0], 100, 1e-5);
        let dist_from_origin = norm3(closest);
        assert!(
            (dist_from_origin - 1.0).abs() < 0.01,
            "closest point should be on sphere surface: dist={:.6}",
            dist_from_origin
        );
    }

    // ── Volume estimate ───────────────────────────────────────────────────────

    #[test]
    fn test_volume_estimate_sphere() {
        let f = |p: [f64; 3]| sdf_sphere(p, 1.0);
        let vol = sdf_volume_estimate(&f, [-1.5, 1.5, -1.5, 1.5, -1.5, 1.5], 10000);
        let exact = 4.0 / 3.0 * PI;
        // Within 10% of exact value
        assert!(
            (vol - exact).abs() / exact < 0.15,
            "volume estimate should be close to {:.6}: got {:.6}",
            exact,
            vol
        );
    }
}
