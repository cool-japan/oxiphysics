#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Implicit surface representations and operations.
//!
//! This module provides Signed Distance Functions (SDFs) for common geometric
//! primitives, Constructive Solid Geometry (CSG) boolean operations, smooth
//! blending, isosurface extraction via Marching Cubes, Radial Basis Function
//! (RBF) implicit surfaces, and sphere-tracing ray marching.
//!
//! # Quick start
//!
//! ```no_run
//! use oxiphysics_geometry::implicit_surfaces::{ImplicitSurface, SphereSDF};
//!
//! let sphere = SphereSDF::new([0.0, 0.0, 0.0], 1.0);
//! assert!(sphere.eval([0.0, 0.0, 0.0]) < 0.0); // inside
//! assert!(sphere.eval([2.0, 0.0, 0.0]) > 0.0); // outside
//! ```

#![allow(dead_code)]

// ---------------------------------------------------------------------------
// Math helpers (private)
// ---------------------------------------------------------------------------

#[inline]
fn len3(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn scale3(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}

#[inline]
fn normalize3(v: [f64; 3]) -> [f64; 3] {
    let l = len3(v).max(1e-300);
    [v[0] / l, v[1] / l, v[2] / l]
}

#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Clamp `v` to `[lo, hi]`.
#[inline]
fn clamp(v: f64, lo: f64, hi: f64) -> f64 {
    v.max(lo).min(hi)
}

// ---------------------------------------------------------------------------
// § 1  Core trait
// ---------------------------------------------------------------------------

/// Trait for implicit surface representations (Signed Distance Functions).
///
/// Implementors define a scalar field f: ℝ³ → ℝ where:
/// - f(p) < 0 means p is inside the surface
/// - f(p) = 0 means p is on the surface
/// - f(p) > 0 means p is outside the surface
pub trait ImplicitSurface {
    /// Evaluate the signed distance (or implicit) function at point `p`.
    fn eval(&self, p: [f64; 3]) -> f64;

    /// Compute the gradient ∇f at point `p`.
    ///
    /// For exact SDFs this equals the outward unit normal on the surface.
    /// The default implementation uses central differences with step `1e-5`.
    fn gradient(&self, p: [f64; 3]) -> [f64; 3] {
        let eps = 1e-5;
        let dx = self.eval([p[0] + eps, p[1], p[2]]) - self.eval([p[0] - eps, p[1], p[2]]);
        let dy = self.eval([p[0], p[1] + eps, p[2]]) - self.eval([p[0], p[1] - eps, p[2]]);
        let dz = self.eval([p[0], p[1], p[2] + eps]) - self.eval([p[0], p[1], p[2] - eps]);
        [dx / (2.0 * eps), dy / (2.0 * eps), dz / (2.0 * eps)]
    }
}

// ---------------------------------------------------------------------------
// § 2  Primitive SDFs
// ---------------------------------------------------------------------------

/// Sphere SDF: f(p) = |p − center| − radius.
///
/// # Example
/// ```no_run
/// use oxiphysics_geometry::implicit_surfaces::{ImplicitSurface, SphereSDF};
/// let s = SphereSDF::new([0.0, 0.0, 0.0], 2.0);
/// assert!((s.eval([2.0, 0.0, 0.0])).abs() < 1e-12);
/// ```
#[derive(Debug, Clone)]
pub struct SphereSDF {
    /// Centre of the sphere.
    pub center: [f64; 3],
    /// Radius of the sphere (must be positive).
    pub radius: f64,
}

impl SphereSDF {
    /// Create a new `SphereSDF`.
    pub fn new(center: [f64; 3], radius: f64) -> Self {
        Self { center, radius }
    }
}

impl ImplicitSurface for SphereSDF {
    fn eval(&self, p: [f64; 3]) -> f64 {
        len3(sub3(p, self.center)) - self.radius
    }

    fn gradient(&self, p: [f64; 3]) -> [f64; 3] {
        normalize3(sub3(p, self.center))
    }
}

/// Axis-aligned box SDF centred at the origin.
///
/// `half_extents[i]` is the half-size along axis i.
///
/// # Example
/// ```no_run
/// use oxiphysics_geometry::implicit_surfaces::{ImplicitSurface, BoxSDF};
/// let b = BoxSDF::new([1.0, 2.0, 3.0]);
/// // A point on a face
/// let d = b.eval([1.0, 0.0, 0.0]);
/// assert!(d.abs() < 1e-12);
/// ```
#[derive(Debug, Clone)]
pub struct BoxSDF {
    /// Half-extents along each axis.
    pub half_extents: [f64; 3],
}

impl BoxSDF {
    /// Create a new axis-aligned box SDF.
    pub fn new(half_extents: [f64; 3]) -> Self {
        Self { half_extents }
    }
}

impl ImplicitSurface for BoxSDF {
    fn eval(&self, p: [f64; 3]) -> f64 {
        // q = |p| - half_extents
        let q = [
            p[0].abs() - self.half_extents[0],
            p[1].abs() - self.half_extents[1],
            p[2].abs() - self.half_extents[2],
        ];
        // Length of the positive part plus the max of the negative part
        let pos = [q[0].max(0.0), q[1].max(0.0), q[2].max(0.0)];
        len3(pos) + q[0].max(q[1]).max(q[2]).min(0.0)
    }
}

/// Torus SDF centred at the origin, lying in the XZ plane.
///
/// - `major_r` — distance from the torus centre to the tube centre
/// - `minor_r` — radius of the tube
///
/// # Example
/// ```no_run
/// use oxiphysics_geometry::implicit_surfaces::{ImplicitSurface, TorusSDF};
/// let t = TorusSDF::new(2.0, 0.5);
/// // Point on the outer equator of the torus
/// let d = t.eval([2.5, 0.0, 0.0]);
/// assert!(d.abs() < 1e-12);
/// ```
#[derive(Debug, Clone)]
pub struct TorusSDF {
    /// Major radius (distance from origin to tube centre).
    pub major_r: f64,
    /// Minor radius (tube radius).
    pub minor_r: f64,
}

impl TorusSDF {
    /// Create a new `TorusSDF`.
    pub fn new(major_r: f64, minor_r: f64) -> Self {
        Self { major_r, minor_r }
    }
}

impl ImplicitSurface for TorusSDF {
    fn eval(&self, p: [f64; 3]) -> f64 {
        // Project onto XZ plane
        let q_xz = ((p[0] * p[0] + p[2] * p[2]).sqrt() - self.major_r, p[1]);
        (q_xz.0 * q_xz.0 + q_xz.1 * q_xz.1).sqrt() - self.minor_r
    }
}

/// Infinite cylinder SDF along the Y axis centred at the origin.
///
/// - `radius` — cylinder radius
/// - `half_height` — half-length along the Y axis
///
/// # Example
/// ```no_run
/// use oxiphysics_geometry::implicit_surfaces::{ImplicitSurface, CylinderSDF};
/// let c = CylinderSDF::new(1.0, 2.0);
/// // Point on lateral surface, within height bounds
/// let d = c.eval([1.0, 0.0, 0.0]);
/// assert!(d.abs() < 1e-12);
/// ```
#[derive(Debug, Clone)]
pub struct CylinderSDF {
    /// Radius of the cylinder.
    pub radius: f64,
    /// Half-height along the Y axis.
    pub half_height: f64,
}

impl CylinderSDF {
    /// Create a new `CylinderSDF`.
    pub fn new(radius: f64, half_height: f64) -> Self {
        Self {
            radius,
            half_height,
        }
    }
}

impl ImplicitSurface for CylinderSDF {
    fn eval(&self, p: [f64; 3]) -> f64 {
        let d = [
            (p[0] * p[0] + p[2] * p[2]).sqrt() - self.radius,
            p[1].abs() - self.half_height,
        ];
        d[0].max(d[1]).min(0.0)
            + [d[0].max(0.0), d[1].max(0.0)]
                .iter()
                .map(|v| v * v)
                .sum::<f64>()
                .sqrt()
    }
}

// ---------------------------------------------------------------------------
// § 3  CSG boolean operations
// ---------------------------------------------------------------------------

/// CSG union of two implicit surfaces: f = min(A, B).
///
/// Points inside either surface are inside the union.
#[derive(Debug, Clone)]
pub struct UnionSDF<A, B> {
    /// First operand.
    pub a: A,
    /// Second operand.
    pub b: B,
}

impl<A: ImplicitSurface, B: ImplicitSurface> UnionSDF<A, B> {
    /// Create a CSG union of `a` and `b`.
    pub fn new(a: A, b: B) -> Self {
        Self { a, b }
    }
}

impl<A: ImplicitSurface, B: ImplicitSurface> ImplicitSurface for UnionSDF<A, B> {
    fn eval(&self, p: [f64; 3]) -> f64 {
        self.a.eval(p).min(self.b.eval(p))
    }
}

/// CSG intersection of two implicit surfaces: f = max(A, B).
///
/// Points inside both surfaces are inside the intersection.
#[derive(Debug, Clone)]
pub struct IntersectionSDF<A, B> {
    /// First operand.
    pub a: A,
    /// Second operand.
    pub b: B,
}

impl<A: ImplicitSurface, B: ImplicitSurface> IntersectionSDF<A, B> {
    /// Create a CSG intersection of `a` and `b`.
    pub fn new(a: A, b: B) -> Self {
        Self { a, b }
    }
}

impl<A: ImplicitSurface, B: ImplicitSurface> ImplicitSurface for IntersectionSDF<A, B> {
    fn eval(&self, p: [f64; 3]) -> f64 {
        self.a.eval(p).max(self.b.eval(p))
    }
}

/// CSG difference of two implicit surfaces: f = max(A, −B).
///
/// Subtracts shape B from shape A.
#[derive(Debug, Clone)]
pub struct DifferenceSDF<A, B> {
    /// Shape to subtract from.
    pub a: A,
    /// Shape to subtract.
    pub b: B,
}

impl<A: ImplicitSurface, B: ImplicitSurface> DifferenceSDF<A, B> {
    /// Create a CSG difference: `a` minus `b`.
    pub fn new(a: A, b: B) -> Self {
        Self { a, b }
    }
}

impl<A: ImplicitSurface, B: ImplicitSurface> ImplicitSurface for DifferenceSDF<A, B> {
    fn eval(&self, p: [f64; 3]) -> f64 {
        self.a.eval(p).max(-self.b.eval(p))
    }
}

// ---------------------------------------------------------------------------
// § 4  Smooth boolean operations
// ---------------------------------------------------------------------------

/// Smooth (C¹) union of two implicit surfaces using polynomial blending.
///
/// The parameter `k` controls the blend radius.  Larger `k` gives smoother
/// transitions.
#[derive(Debug, Clone)]
pub struct SmoothUnionSDF<A, B> {
    /// First operand.
    pub a: A,
    /// Second operand.
    pub b: B,
    /// Blend radius (must be > 0).
    pub k: f64,
}

impl<A: ImplicitSurface, B: ImplicitSurface> SmoothUnionSDF<A, B> {
    /// Create a smooth union with blend parameter `k`.
    pub fn new(a: A, b: B, k: f64) -> Self {
        Self { a, b, k }
    }
}

impl<A: ImplicitSurface, B: ImplicitSurface> ImplicitSurface for SmoothUnionSDF<A, B> {
    fn eval(&self, p: [f64; 3]) -> f64 {
        let da = self.a.eval(p);
        let db = self.b.eval(p);
        let h = clamp(0.5 + 0.5 * (db - da) / self.k, 0.0, 1.0);
        da * h + db * (1.0 - h) - self.k * h * (1.0 - h)
    }
}

// ---------------------------------------------------------------------------
// § 5  Offset SDF (shell / thickening)
// ---------------------------------------------------------------------------

/// Offset (shell) SDF: expands or shrinks a shape by a fixed distance.
///
/// `f_offset(p) = f(p) − offset`
///
/// Positive `offset` expands the shape; negative `offset` shrinks it.
#[derive(Debug, Clone)]
pub struct OffsetSDF<S> {
    /// Inner surface.
    pub inner: S,
    /// Offset distance (m).
    pub offset: f64,
}

impl<S: ImplicitSurface> OffsetSDF<S> {
    /// Create an offset surface from `inner` shifted by `offset`.
    pub fn new(inner: S, offset: f64) -> Self {
        Self { inner, offset }
    }
}

impl<S: ImplicitSurface> ImplicitSurface for OffsetSDF<S> {
    fn eval(&self, p: [f64; 3]) -> f64 {
        self.inner.eval(p) - self.offset
    }

    fn gradient(&self, p: [f64; 3]) -> [f64; 3] {
        self.inner.gradient(p)
    }
}

// ---------------------------------------------------------------------------
// § 6  Marching Cubes
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// Marching cubes lookup tables (edge table and tri table)
// ---------------------------------------------------------------------------

/// Edge table for marching cubes (256 entries).
const EDGE_TABLE: [u16; 256] = [
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
    0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x835, 0xb3f, 0xa36, // fixed 0x83f→0x835
    0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x033, 0x339, 0x230, 0xe90, 0xf99, 0xc93, 0xd9a, 0xa96,
    0xb9f, 0x895, 0x99c, 0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x099, 0x190, 0xf00, 0xe09,
    0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c, 0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109,
    0x000,
];

/// Cube corner positions in unit cube (index → \[xi, yi, zi\]).
const CUBE_CORNERS: [[f64; 3]; 8] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [1.0, 1.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
    [1.0, 0.0, 1.0],
    [1.0, 1.0, 1.0],
    [0.0, 1.0, 1.0],
];

/// Cube edge pairs (each edge connects two corners).
const CUBE_EDGES: [[usize; 2]; 12] = [
    [0, 1],
    [1, 2],
    [2, 3],
    [3, 0],
    [4, 5],
    [5, 6],
    [6, 7],
    [7, 4],
    [0, 4],
    [1, 5],
    [2, 6],
    [3, 7],
];

/// Interpolate a vertex on an edge where the SDF crosses zero.
fn interp_edge(p0: [f64; 3], v0: f64, p1: [f64; 3], v1: f64) -> [f64; 3] {
    let dv = v1 - v0;
    if dv.abs() < 1e-30 {
        return [
            (p0[0] + p1[0]) / 2.0,
            (p0[1] + p1[1]) / 2.0,
            (p0[2] + p1[2]) / 2.0,
        ];
    }
    let t = -v0 / dv;
    [
        p0[0] + t * (p1[0] - p0[0]),
        p0[1] + t * (p1[1] - p0[1]),
        p0[2] + t * (p1[2] - p0[2]),
    ]
}

/// Extract an isosurface from an implicit SDF using the Marching Cubes algorithm.
///
/// # Arguments
/// * `sdf` — any type implementing [`ImplicitSurface`]
/// * `bounds` — axis-aligned bounding box as `[[xmin,xmax\], [ymin,ymax], [zmin,zmax]]`
/// * `resolution` — number of grid cells per axis
///
/// # Returns
/// A tuple `(vertices, triangles)` where each vertex is `[x, y, z]` and
/// each triangle is `[i0, i1, i2]` (indices into `vertices`).
pub fn marching_cubes(
    sdf: &dyn ImplicitSurface,
    bounds: [[f64; 2]; 3],
    resolution: usize,
) -> (Vec<[f64; 3]>, Vec<[usize; 3]>) {
    let n = resolution.max(1);
    let mut vertices: Vec<[f64; 3]> = Vec::new();
    let mut triangles: Vec<[usize; 3]> = Vec::new();

    let dx = (bounds[0][1] - bounds[0][0]) / n as f64;
    let dy = (bounds[1][1] - bounds[1][0]) / n as f64;
    let dz = (bounds[2][1] - bounds[2][0]) / n as f64;

    for iz in 0..n {
        for iy in 0..n {
            for ix in 0..n {
                let ox = bounds[0][0] + ix as f64 * dx;
                let oy = bounds[1][0] + iy as f64 * dy;
                let oz = bounds[2][0] + iz as f64 * dz;

                // Evaluate SDF at 8 corners
                let corner_pts: [[f64; 3]; 8] = std::array::from_fn(|c| {
                    [
                        ox + CUBE_CORNERS[c][0] * dx,
                        oy + CUBE_CORNERS[c][1] * dy,
                        oz + CUBE_CORNERS[c][2] * dz,
                    ]
                });
                let vals: [f64; 8] = std::array::from_fn(|c| sdf.eval(corner_pts[c]));

                // Compute cube index
                let mut cube_idx: usize = 0;
                for c in 0..8 {
                    if vals[c] < 0.0 {
                        cube_idx |= 1 << c;
                    }
                }

                let edge_flags = EDGE_TABLE[cube_idx];
                if edge_flags == 0 {
                    continue;
                }

                // Compute intersection vertices for active edges
                let mut edge_verts: [Option<[f64; 3]>; 12] = [None; 12];
                for e in 0..12 {
                    if edge_flags & (1 << e) != 0 {
                        let [i0, i1] = CUBE_EDGES[e];
                        edge_verts[e] = Some(interp_edge(
                            corner_pts[i0],
                            vals[i0],
                            corner_pts[i1],
                            vals[i1],
                        ));
                    }
                }

                // Tessellate using a simple fan from the first active edge vertex
                // (simplified triangulation — not full marching-cubes tri table)
                let active: Vec<[f64; 3]> = (0..12).filter_map(|e| edge_verts[e]).collect();
                if active.len() >= 3 {
                    let base = vertices.len();
                    for v in &active {
                        vertices.push(*v);
                    }
                    for k in 1..(active.len() - 1) {
                        triangles.push([base, base + k, base + k + 1]);
                    }
                }
            }
        }
    }

    (vertices, triangles)
}

// ---------------------------------------------------------------------------
// § 7  RBF Implicit Surface
// ---------------------------------------------------------------------------

/// Radial Basis Function (RBF) implicit surface.
///
/// The surface is defined as f(p) = Σᵢ wᵢ φ(|p − cᵢ|) where φ is a
/// thin-plate spline kernel φ(r) = r² log(r + ε).
///
/// Use [`RBFImplicit::fit`] to compute weights from a set of point samples
/// with target values, then [`RBFImplicit::eval`] to query the surface.
#[derive(Debug, Clone)]
pub struct RBFImplicit {
    /// RBF centre positions.
    pub centers: Vec<[f64; 3]>,
    /// Weights for each centre (computed by [`fit`](Self::fit)).
    pub weights: Vec<f64>,
}

impl RBFImplicit {
    /// Create a new `RBFImplicit` with given centres and zero weights.
    pub fn new(centers: Vec<[f64; 3]>) -> Self {
        let n = centers.len();
        Self {
            centers,
            weights: vec![0.0; n],
        }
    }

    /// Fit weights from sample points `pts` and target values `targets`.
    ///
    /// Solves the symmetric RBF system Φ w = t using the conjugate gradient
    /// method (up to `max_iter` iterations, tolerance `tol`).
    ///
    /// # Panics
    /// Panics if `pts.len() != targets.len()` or `pts.len() != centers.len()`.
    pub fn fit(&mut self, pts: &[[f64; 3]], targets: &[f64], tol: f64, max_iter: usize) {
        assert_eq!(pts.len(), targets.len());
        assert_eq!(pts.len(), self.centers.len());
        let n = pts.len();
        // Build the RBF matrix Φ
        let phi: Vec<Vec<f64>> = (0..n)
            .map(|i| {
                (0..n)
                    .map(|j| rbf_kernel(len3(sub3(pts[i], self.centers[j]))))
                    .collect()
            })
            .collect();
        // Solve Φ w = targets
        self.weights = rbf_cg(&phi, targets, tol, max_iter);
    }

    /// Evaluate the RBF implicit surface at point `p`.
    pub fn eval(&self, p: [f64; 3]) -> f64 {
        self.centers
            .iter()
            .zip(self.weights.iter())
            .map(|(c, w)| w * rbf_kernel(len3(sub3(p, *c))))
            .sum()
    }
}

/// Thin-plate spline RBF kernel: φ(r) = r² ln(r + ε).
fn rbf_kernel(r: f64) -> f64 {
    r * r * (r + 1e-30).ln()
}

/// Simple CG solver for the RBF system.
fn rbf_cg(a: &[Vec<f64>], b: &[f64], tol: f64, max_iter: usize) -> Vec<f64> {
    let n = b.len();
    let mut x = vec![0.0_f64; n];
    let mut r: Vec<f64> = b.to_vec();
    let mut p = r.clone();
    let mut rsold: f64 = r.iter().map(|v| v * v).sum();

    for _ in 0..max_iter {
        let ap: Vec<f64> = (0..n)
            .map(|i| a[i].iter().zip(p.iter()).map(|(aij, pj)| aij * pj).sum())
            .collect();
        let pap: f64 = p.iter().zip(ap.iter()).map(|(pi, api)| pi * api).sum();
        if pap.abs() < 1e-300 {
            break;
        }
        let alpha = rsold / pap;
        for i in 0..n {
            x[i] += alpha * p[i];
            r[i] -= alpha * ap[i];
        }
        let rsnew: f64 = r.iter().map(|v| v * v).sum();
        if rsnew.sqrt() < tol {
            break;
        }
        let beta = rsnew / rsold;
        for i in 0..n {
            p[i] = r[i] + beta * p[i];
        }
        rsold = rsnew;
    }
    x
}

// ---------------------------------------------------------------------------
// § 8  Ray marching (sphere tracing)
// ---------------------------------------------------------------------------

/// Sphere-trace a ray against an implicit surface.
///
/// Starting from `origin` in direction `direction` (need not be unit), takes
/// up to `max_steps` sphere-tracing steps.  Returns the distance along the
/// ray to the first surface hit, or `None` if no hit was found within
/// `max_dist`.
///
/// # Arguments
/// * `sdf` — implicit surface to trace against
/// * `origin` — ray origin in world space
/// * `direction` — ray direction (will be normalized internally)
/// * `max_steps` — maximum number of steps
///
/// The function terminates early when `|f(p)| < 1e-5`.
pub fn ray_march(
    sdf: &dyn ImplicitSurface,
    origin: [f64; 3],
    direction: [f64; 3],
    max_steps: usize,
) -> Option<f64> {
    let dir = normalize3(direction);
    let max_dist = 1e6_f64;
    let mut t = 0.0_f64;

    for _ in 0..max_steps {
        let p = add3(origin, scale3(dir, t));
        let d = sdf.eval(p);
        if d.abs() < 1e-5 {
            return Some(t);
        }
        if t > max_dist {
            break;
        }
        // Advance by |d| (sphere-tracing step)
        t += d.abs().max(1e-7);
    }
    None
}

// ---------------------------------------------------------------------------
// § 9  Normal estimation
// ---------------------------------------------------------------------------

/// Estimate the outward unit normal of an implicit surface at `p` using
/// central differences with step size `eps`.
///
/// # Arguments
/// * `sdf` — implicit surface
/// * `p` — query point (should be near or on the surface)
/// * `eps` — finite difference step (typical: 1e-4 to 1e-6)
pub fn surface_normal(sdf: &dyn ImplicitSurface, p: [f64; 3], eps: f64) -> [f64; 3] {
    let gx = sdf.eval([p[0] + eps, p[1], p[2]]) - sdf.eval([p[0] - eps, p[1], p[2]]);
    let gy = sdf.eval([p[0], p[1] + eps, p[2]]) - sdf.eval([p[0], p[1] - eps, p[2]]);
    let gz = sdf.eval([p[0], p[1], p[2] + eps]) - sdf.eval([p[0], p[1], p[2] - eps]);
    normalize3([gx, gy, gz])
}

// ---------------------------------------------------------------------------
// § 10  Point-in-surface test (ray casting)
// ---------------------------------------------------------------------------

/// Test whether a point is inside a closed implicit surface.
///
/// Returns `true` if `sdf.eval(p) < 0`.
///
/// # Arguments
/// * `sdf` — implicit surface
/// * `p` — query point
pub fn is_inside(sdf: &dyn ImplicitSurface, p: [f64; 3]) -> bool {
    sdf.eval(p) < 0.0
}

// ---------------------------------------------------------------------------
// § 11  Bounding-box computation
// ---------------------------------------------------------------------------

/// Compute the axis-aligned bounding box of a level-set surface by sampling.
///
/// Samples the SDF on a uniform grid within `bounds` at resolution `res`
/// per axis and returns `[[xmin,xmax\], [ymin,ymax], [zmin,zmax]]` of all
/// points where `|f(p)| < threshold`.
pub fn surface_bbox(
    sdf: &dyn ImplicitSurface,
    bounds: [[f64; 2]; 3],
    res: usize,
    threshold: f64,
) -> [[f64; 2]; 3] {
    let n = res.max(2);
    let mut bbox = [
        [f64::INFINITY, f64::NEG_INFINITY],
        [f64::INFINITY, f64::NEG_INFINITY],
        [f64::INFINITY, f64::NEG_INFINITY],
    ];

    let step: [f64; 3] = std::array::from_fn(|k| (bounds[k][1] - bounds[k][0]) / (n - 1) as f64);

    for iz in 0..n {
        for iy in 0..n {
            for ix in 0..n {
                let p = [
                    bounds[0][0] + ix as f64 * step[0],
                    bounds[1][0] + iy as f64 * step[1],
                    bounds[2][0] + iz as f64 * step[2],
                ];
                if sdf.eval(p).abs() < threshold {
                    for k in 0..3 {
                        bbox[k][0] = bbox[k][0].min(p[k]);
                        bbox[k][1] = bbox[k][1].max(p[k]);
                    }
                }
            }
        }
    }
    bbox
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ------------------------------------------------------------------
    // SphereSDF
    // ------------------------------------------------------------------

    #[test]
    fn test_sphere_eval_center_negative() {
        let s = SphereSDF::new([0.0, 0.0, 0.0], 1.0);
        assert!(s.eval([0.0, 0.0, 0.0]) < 0.0);
    }

    #[test]
    fn test_sphere_eval_on_surface() {
        let s = SphereSDF::new([0.0, 0.0, 0.0], 2.0);
        let d = s.eval([2.0, 0.0, 0.0]);
        assert!(d.abs() < 1e-12, "d = {d}");
    }

    #[test]
    fn test_sphere_eval_outside_positive() {
        let s = SphereSDF::new([0.0, 0.0, 0.0], 1.0);
        assert!(s.eval([3.0, 0.0, 0.0]) > 0.0);
    }

    #[test]
    fn test_sphere_gradient_unit_normal() {
        let s = SphereSDF::new([0.0, 0.0, 0.0], 1.0);
        let g = s.gradient([1.0, 0.0, 0.0]);
        let l = len3(g);
        assert!((l - 1.0).abs() < 1e-12, "gradient not unit: |g| = {l}");
    }

    #[test]
    fn test_sphere_gradient_direction() {
        let s = SphereSDF::new([0.0, 0.0, 0.0], 1.0);
        let g = s.gradient([0.0, 1.0, 0.0]);
        // Should point in +y direction
        assert!(g[1] > 0.9);
        assert!(g[0].abs() < 1e-12);
        assert!(g[2].abs() < 1e-12);
    }

    #[test]
    fn test_sphere_translated() {
        let s = SphereSDF::new([1.0, 2.0, 3.0], 1.5);
        // Distance from center should be -1.5
        assert!((s.eval([1.0, 2.0, 3.0]) + 1.5).abs() < 1e-12);
    }

    // ------------------------------------------------------------------
    // BoxSDF
    // ------------------------------------------------------------------

    #[test]
    fn test_box_eval_inside_negative() {
        let b = BoxSDF::new([1.0, 1.0, 1.0]);
        assert!(b.eval([0.0, 0.0, 0.0]) < 0.0);
    }

    #[test]
    fn test_box_eval_on_face() {
        let b = BoxSDF::new([1.0, 2.0, 3.0]);
        let d = b.eval([1.0, 0.0, 0.0]);
        assert!(d.abs() < 1e-12, "d = {d}");
    }

    #[test]
    fn test_box_eval_outside_positive() {
        let b = BoxSDF::new([1.0, 1.0, 1.0]);
        assert!(b.eval([2.0, 0.0, 0.0]) > 0.0);
    }

    #[test]
    fn test_box_eval_corner() {
        let b = BoxSDF::new([1.0, 1.0, 1.0]);
        // Corner is at sqrt(3) distance from the box surface
        let d = b.eval([2.0, 2.0, 2.0]);
        let expected = (3.0_f64).sqrt();
        assert!((d - expected).abs() < 1e-12);
    }

    // ------------------------------------------------------------------
    // TorusSDF
    // ------------------------------------------------------------------

    #[test]
    fn test_torus_eval_on_surface() {
        let t = TorusSDF::new(2.0, 0.5);
        // Outermost point of the torus
        let d = t.eval([2.5, 0.0, 0.0]);
        assert!(d.abs() < 1e-12, "d = {d}");
    }

    #[test]
    fn test_torus_eval_center_positive() {
        let t = TorusSDF::new(2.0, 0.5);
        // The center of the torus is outside (hollow)
        assert!(t.eval([0.0, 0.0, 0.0]) > 0.0);
    }

    #[test]
    fn test_torus_inside_tube_negative() {
        let t = TorusSDF::new(2.0, 0.5);
        // Point near tube centre
        assert!(t.eval([2.0, 0.0, 0.0]) < 0.0);
    }

    // ------------------------------------------------------------------
    // CylinderSDF
    // ------------------------------------------------------------------

    #[test]
    fn test_cylinder_lateral_surface() {
        let c = CylinderSDF::new(1.0, 2.0);
        let d = c.eval([1.0, 0.0, 0.0]);
        assert!(d.abs() < 1e-12, "d = {d}");
    }

    #[test]
    fn test_cylinder_inside_negative() {
        let c = CylinderSDF::new(1.0, 2.0);
        assert!(c.eval([0.0, 0.0, 0.0]) < 0.0);
    }

    #[test]
    fn test_cylinder_outside_positive() {
        let c = CylinderSDF::new(1.0, 2.0);
        assert!(c.eval([3.0, 0.0, 0.0]) > 0.0);
    }

    // ------------------------------------------------------------------
    // CSG operations
    // ------------------------------------------------------------------

    #[test]
    fn test_union_inside_either() {
        let a = SphereSDF::new([-1.0, 0.0, 0.0], 1.5);
        let b = SphereSDF::new([1.0, 0.0, 0.0], 1.5);
        let u = UnionSDF::new(a, b);
        // Points inside both spheres
        assert!(u.eval([-1.0, 0.0, 0.0]) < 0.0);
        assert!(u.eval([1.0, 0.0, 0.0]) < 0.0);
    }

    #[test]
    fn test_intersection_inside_both() {
        let a = SphereSDF::new([0.0, 0.0, 0.0], 2.0);
        let b = SphereSDF::new([0.0, 0.0, 0.0], 1.0);
        let i = IntersectionSDF::new(a, b);
        // Inside smaller sphere → inside intersection
        assert!(i.eval([0.0, 0.0, 0.0]) < 0.0);
    }

    #[test]
    fn test_intersection_outside_one() {
        let a = SphereSDF::new([0.0, 0.0, 0.0], 2.0);
        let b = SphereSDF::new([5.0, 0.0, 0.0], 1.0);
        let i = IntersectionSDF::new(a, b);
        // Origin is inside A but not inside B → outside intersection
        assert!(i.eval([0.0, 0.0, 0.0]) > 0.0);
    }

    #[test]
    fn test_difference_removes_region() {
        let a = SphereSDF::new([0.0, 0.0, 0.0], 2.0);
        let b = SphereSDF::new([0.0, 0.0, 0.0], 1.0);
        let d = DifferenceSDF::new(a, b);
        // Points inside b are removed
        assert!(d.eval([0.5, 0.0, 0.0]) > 0.0);
        // Points in a but not b remain
        assert!(d.eval([1.5, 0.0, 0.0]) < 0.0);
    }

    // ------------------------------------------------------------------
    // SmoothUnionSDF
    // ------------------------------------------------------------------

    #[test]
    fn test_smooth_union_less_than_regular_union() {
        let a = SphereSDF::new([-0.5, 0.0, 0.0], 1.0);
        let b = SphereSDF::new([0.5, 0.0, 0.0], 1.0);
        let su = SmoothUnionSDF::new(a.clone(), b.clone(), 0.5);
        let u = UnionSDF::new(a, b);
        let p = [0.0, 0.0, 0.0];
        // Smooth union should be <= regular union (blending pulls inward)
        assert!(su.eval(p) <= u.eval(p) + 1e-12);
    }

    // ------------------------------------------------------------------
    // OffsetSDF
    // ------------------------------------------------------------------

    #[test]
    fn test_offset_expands_sphere() {
        let s = SphereSDF::new([0.0, 0.0, 0.0], 1.0);
        let os = OffsetSDF::new(s, 0.5);
        // New surface at r = 1.5
        let d = os.eval([1.5, 0.0, 0.0]);
        assert!(d.abs() < 1e-12, "d = {d}");
    }

    #[test]
    fn test_offset_shrinks_sphere() {
        let s = SphereSDF::new([0.0, 0.0, 0.0], 2.0);
        let os = OffsetSDF::new(s, -0.5);
        // New surface at r = 1.5
        let d = os.eval([1.5, 0.0, 0.0]);
        assert!(d.abs() < 1e-12, "d = {d}");
    }

    // ------------------------------------------------------------------
    // is_inside / surface_normal
    // ------------------------------------------------------------------

    #[test]
    fn test_is_inside_sphere() {
        let s = SphereSDF::new([0.0, 0.0, 0.0], 1.0);
        assert!(is_inside(&s, [0.0, 0.0, 0.0]));
        assert!(!is_inside(&s, [2.0, 0.0, 0.0]));
    }

    #[test]
    fn test_surface_normal_sphere_unit_length() {
        let s = SphereSDF::new([0.0, 0.0, 0.0], 1.0);
        let n = surface_normal(&s, [1.0, 0.0, 0.0], 1e-5);
        let l = len3(n);
        assert!((l - 1.0).abs() < 1e-3, "normal length = {l}");
    }

    // ------------------------------------------------------------------
    // Marching cubes
    // ------------------------------------------------------------------

    #[test]
    fn test_marching_cubes_sphere_produces_vertices() {
        let s = SphereSDF::new([0.0, 0.0, 0.0], 0.4);
        let bounds = [[-1.0, 1.0], [-1.0, 1.0], [-1.0, 1.0]];
        let (verts, tris) = marching_cubes(&s, bounds, 8);
        assert!(!verts.is_empty(), "no vertices produced");
        assert!(!tris.is_empty(), "no triangles produced");
    }

    #[test]
    fn test_marching_cubes_indices_in_bounds() {
        let s = SphereSDF::new([0.0, 0.0, 0.0], 0.4);
        let bounds = [[-1.0, 1.0], [-1.0, 1.0], [-1.0, 1.0]];
        let (verts, tris) = marching_cubes(&s, bounds, 6);
        for tri in &tris {
            for &idx in tri {
                assert!(idx < verts.len(), "index {idx} out of range");
            }
        }
    }

    #[test]
    fn test_marching_cubes_empty_for_no_surface() {
        // Sphere entirely outside the bounds
        let s = SphereSDF::new([100.0, 100.0, 100.0], 0.1);
        let bounds = [[-1.0, 1.0], [-1.0, 1.0], [-1.0, 1.0]];
        let (verts, tris) = marching_cubes(&s, bounds, 4);
        assert!(verts.is_empty());
        assert!(tris.is_empty());
    }

    // ------------------------------------------------------------------
    // RBFImplicit
    // ------------------------------------------------------------------

    #[test]
    fn test_rbf_fit_interpolates_at_centers() {
        let centers = vec![[0.0, 0.0, 0.0_f64], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let targets = vec![0.0, 1.0, -1.0];
        let mut rbf = RBFImplicit::new(centers.clone());
        rbf.fit(&centers, &targets, 1e-10, 200);
        // After fitting, evaluating at centers should approximate targets
        for (c, &t) in centers.iter().zip(targets.iter()) {
            let v = rbf.eval(*c);
            assert!((v - t).abs() < 0.5, "RBF({c:?}) = {v}, expected {t}");
        }
    }

    #[test]
    fn test_rbf_zero_weights_before_fit() {
        let centers = vec![[0.0, 0.0, 0.0_f64]];
        let rbf = RBFImplicit::new(centers);
        assert_eq!(rbf.weights, vec![0.0]);
    }

    // ------------------------------------------------------------------
    // Ray marching
    // ------------------------------------------------------------------

    #[test]
    fn test_ray_march_sphere_hit() {
        let s = SphereSDF::new([0.0, 0.0, 0.0], 1.0);
        let hit = ray_march(&s, [-5.0, 0.0, 0.0], [1.0, 0.0, 0.0], 200);
        assert!(hit.is_some(), "ray should hit sphere");
        let t = hit.unwrap();
        // Hit point should be near [-1, 0, 0]
        let px = -5.0 + t;
        assert!((px + 1.0).abs() < 0.01, "hit at x = {px}, expected -1");
    }

    #[test]
    fn test_ray_march_sphere_miss() {
        let s = SphereSDF::new([0.0, 0.0, 0.0], 1.0);
        // Ray aimed far away from sphere
        let hit = ray_march(&s, [0.0, 100.0, 0.0], [1.0, 0.0, 0.0], 50);
        assert!(hit.is_none(), "ray should miss sphere");
    }

    #[test]
    fn test_ray_march_returns_positive_distance() {
        let s = SphereSDF::new([0.0, 0.0, 0.0], 1.0);
        let hit = ray_march(&s, [-5.0, 0.0, 0.0], [1.0, 0.0, 0.0], 200);
        if let Some(t) = hit {
            assert!(t > 0.0, "distance must be positive, got {t}");
        }
    }
}
