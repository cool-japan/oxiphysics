// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Primitive shapes and convex hull types.

use oxiphysics::geometry::convex_hull::{ConvexHull3D, gjk_intersect};
use pyo3::prelude::*;
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Math helpers (local, no nalgebra)
// ---------------------------------------------------------------------------

pub(super) fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

pub(super) fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

pub(super) fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

pub(super) fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

pub(super) fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

pub(super) fn len3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}

pub(super) fn normalize3(a: [f64; 3]) -> [f64; 3] {
    let l = len3(a);
    if l < 1e-15 {
        [0.0; 3]
    } else {
        scale3(a, 1.0 / l)
    }
}

pub(super) fn lerp3(a: [f64; 3], b: [f64; 3], t: f64) -> [f64; 3] {
    [
        a[0] + (b[0] - a[0]) * t,
        a[1] + (b[1] - a[1]) * t,
        a[2] + (b[2] - a[2]) * t,
    ]
}

pub(super) fn centroid(pts: &[[f64; 3]]) -> [f64; 3] {
    let n = pts.len() as f64;
    let mut s = [0.0f64; 3];
    for &p in pts {
        s[0] += p[0];
        s[1] += p[1];
        s[2] += p[2];
    }
    [s[0] / n, s[1] / n, s[2] / n]
}

// ---------------------------------------------------------------------------
// PyShape — primitive geometry enum (Python-facing wrapper class)
// ---------------------------------------------------------------------------

/// Primitive geometry shape used for collision detection and rendering.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum PyShape {
    /// A sphere defined by its radius.
    Sphere(f64),
    /// An axis-aligned box defined by half-extents [hx, hy, hz].
    Box([f64; 3]),
    /// A capsule defined by radius and half-height of the cylindrical part.
    Capsule { radius: f64, half_height: f64 },
    /// A cylinder defined by radius and half-height.
    Cylinder { radius: f64, half_height: f64 },
    /// A cone defined by base radius and half-height.
    Cone { radius: f64, half_height: f64 },
    /// A torus defined by major and minor radii.
    Torus { major_r: f64, minor_r: f64 },
}

impl PyShape {
    /// Approximate volume of the shape.
    pub fn volume(&self) -> f64 {
        match self {
            PyShape::Sphere(r) => (4.0 / 3.0) * PI * r * r * r,
            PyShape::Box(h) => 8.0 * h[0] * h[1] * h[2],
            PyShape::Capsule {
                radius: r,
                half_height: hh,
            } => PI * r * r * (2.0 * hh + (4.0 / 3.0) * r),
            PyShape::Cylinder {
                radius: r,
                half_height: hh,
            } => PI * r * r * 2.0 * hh,
            PyShape::Cone {
                radius: r,
                half_height: hh,
            } => (1.0 / 3.0) * PI * r * r * 2.0 * hh,
            PyShape::Torus {
                major_r: rr,
                minor_r: rt,
            } => 2.0 * PI * PI * rr * rt * rt,
        }
    }

    /// Approximate surface area of the shape.
    pub fn surface_area(&self) -> f64 {
        match self {
            PyShape::Sphere(r) => 4.0 * PI * r * r,
            PyShape::Box(h) => 2.0 * (h[0] * h[1] + h[1] * h[2] + h[2] * h[0]) * 4.0,
            PyShape::Capsule {
                radius: r,
                half_height: hh,
            } => 4.0 * PI * r * r + 2.0 * PI * r * 2.0 * hh,
            PyShape::Cylinder {
                radius: r,
                half_height: hh,
            } => 2.0 * PI * r * (r + 2.0 * hh),
            PyShape::Cone {
                radius: r,
                half_height: hh,
            } => {
                let slant = (r * r + (2.0 * hh) * (2.0 * hh)).sqrt();
                PI * r * (r + slant)
            }
            PyShape::Torus {
                major_r: rr,
                minor_r: rt,
            } => 4.0 * PI * PI * rr * rt,
        }
    }

    /// Compute an axis-aligned bounding box for this shape at the origin.
    pub fn aabb(&self) -> ([f64; 3], [f64; 3]) {
        match self {
            PyShape::Sphere(r) => ([-r, -r, -r], [*r, *r, *r]),
            PyShape::Box(h) => ([-h[0], -h[1], -h[2]], [h[0], h[1], h[2]]),
            PyShape::Capsule {
                radius: r,
                half_height: hh,
            } => {
                let hy = r + hh;
                ([-r, -hy, -r], [*r, hy, *r])
            }
            PyShape::Cylinder {
                radius: r,
                half_height: hh,
            } => ([-r, -hh, -r], [*r, *hh, *r]),
            PyShape::Cone {
                radius: r,
                half_height: hh,
            } => ([-r, -hh, -r], [*r, *hh, *r]),
            PyShape::Torus {
                major_r: rr,
                minor_r: rt,
            } => {
                let e = rr + rt;
                ([-e, -rt, -e], [e, *rt, e])
            }
        }
    }
}

/// Python-facing wrapper for PyShape.
#[pyclass(name = "PyShape", skip_from_py_object)]
#[derive(Debug, Clone)]
pub struct PyShapeWrapper {
    inner: PyShape,
}

#[pymethods]
impl PyShapeWrapper {
    /// Create a sphere shape.
    #[staticmethod]
    pub fn sphere(radius: f64) -> Self {
        Self {
            inner: PyShape::Sphere(radius),
        }
    }

    /// Create a box shape from half-extents [hx, hy, hz].
    #[staticmethod]
    pub fn box_shape(hx: f64, hy: f64, hz: f64) -> Self {
        Self {
            inner: PyShape::Box([hx, hy, hz]),
        }
    }

    /// Create a capsule shape.
    #[staticmethod]
    pub fn capsule(radius: f64, half_height: f64) -> Self {
        Self {
            inner: PyShape::Capsule {
                radius,
                half_height,
            },
        }
    }

    /// Create a cylinder shape.
    #[staticmethod]
    pub fn cylinder(radius: f64, half_height: f64) -> Self {
        Self {
            inner: PyShape::Cylinder {
                radius,
                half_height,
            },
        }
    }

    /// Create a cone shape.
    #[staticmethod]
    pub fn cone(radius: f64, half_height: f64) -> Self {
        Self {
            inner: PyShape::Cone {
                radius,
                half_height,
            },
        }
    }

    /// Create a torus shape.
    #[staticmethod]
    pub fn torus(major_r: f64, minor_r: f64) -> Self {
        Self {
            inner: PyShape::Torus { major_r, minor_r },
        }
    }

    /// Approximate volume of the shape.
    pub fn volume(&self) -> f64 {
        self.inner.volume()
    }

    /// Approximate surface area of the shape.
    pub fn surface_area(&self) -> f64 {
        self.inner.surface_area()
    }

    /// Return the AABB as (min_list, max_list).
    pub fn aabb(&self) -> (Vec<f64>, Vec<f64>) {
        let (mn, mx) = self.inner.aabb();
        (mn.to_vec(), mx.to_vec())
    }
}

// ---------------------------------------------------------------------------
// PyConvexHull
// ---------------------------------------------------------------------------

/// Reasons a convex hull cannot be built from a point set.
///
/// Every variant surfaces in Python as a `ValueError`.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum ConvexHullError {
    /// Fewer than the four points a tetrahedron needs.
    #[error("a convex hull needs at least 4 points, got {0}")]
    TooFewPoints(usize),
    /// A point row has fewer than three coordinates.
    #[error("point {index} has {got} coordinates, expected at least 3")]
    ShortPoint {
        /// Index of the offending point.
        index: usize,
        /// Number of coordinates it carries.
        got: usize,
    },
    /// A coordinate is NaN or infinite.
    #[error("point {0} has a non-finite coordinate")]
    NonFinitePoint(usize),
    /// All points are coplanar (or collinear), so they enclose no volume.
    #[error("the points are coplanar or collinear and enclose no volume")]
    Degenerate,
}

impl From<ConvexHullError> for PyErr {
    fn from(err: ConvexHullError) -> PyErr {
        pyo3::exceptions::PyValueError::new_err(err.to_string())
    }
}

/// A convex hull computed from a point cloud.
///
/// The hull is the Quickhull of the input: `vertices` are the extreme points
/// (a subset of the input) and `faces` are outward-facing triangles.
#[pyclass(skip_from_py_object)]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyConvexHull {
    /// Vertices of the convex hull.
    pub vertices: Vec<[f64; 3]>,
    /// Triangle face indices (groups of 3).
    pub faces: Vec<usize>,
}

#[pymethods]
impl PyConvexHull {
    /// Compute the convex hull of a set of points (a list of `[x, y, z]` lists).
    ///
    /// Raises `ValueError` when a point has fewer than three coordinates or a
    /// non-finite one, when there are fewer than four points, or when all
    /// points are coplanar.
    #[staticmethod]
    pub fn from_points(points: Vec<Vec<f64>>) -> PyResult<Self> {
        let mut pts: Vec<[f64; 3]> = Vec::with_capacity(points.len());
        for (index, p) in points.iter().enumerate() {
            match p.as_slice() {
                [x, y, z, ..] => pts.push([*x, *y, *z]),
                _ => {
                    return Err(ConvexHullError::ShortPoint {
                        index,
                        got: p.len(),
                    }
                    .into());
                }
            }
        }
        Ok(Self::from_points_internal(&pts)?)
    }

    /// Volume enclosed by the convex hull.
    pub fn volume(&self) -> f64 {
        if self.vertices.is_empty() {
            return 0.0;
        }
        // Signed tetrahedra from the vertex centroid: exact for a closed
        // outward-oriented surface and well conditioned far from the origin.
        let c = centroid(&self.vertices);
        let mut vol = 0.0;
        for tri in self.faces.chunks_exact(3) {
            let a = sub3(self.vertices[tri[0]], c);
            let b = sub3(self.vertices[tri[1]], c);
            let d = sub3(self.vertices[tri[2]], c);
            vol += dot3(a, cross3(b, d));
        }
        (vol / 6.0).abs()
    }

    /// Surface area of the convex hull.
    pub fn surface_area(&self) -> f64 {
        let n = self.faces.len() / 3;
        let mut area = 0.0;
        for i in 0..n {
            let a = self.vertices[self.faces[i * 3]];
            let b = self.vertices[self.faces[i * 3 + 1]];
            let c = self.vertices[self.faces[i * 3 + 2]];
            let ab = sub3(b, a);
            let ac = sub3(c, a);
            area += len3(cross3(ab, ac)) * 0.5;
        }
        area
    }

    /// Test whether two convex hulls overlap.
    ///
    /// Runs the GJK algorithm (`oxiphysics.geometry.convex_hull.gjk_intersect`)
    /// on the hull vertices, so the answer concerns the hulls themselves and
    /// not their bounding spheres.
    pub fn gjk_overlap(&self, other: &PyConvexHull) -> bool {
        gjk_intersect(&self.engine_hull(), &other.engine_hull())
    }

    /// Get vertices as a list of [x, y, z] lists.
    pub fn get_vertices(&self) -> Vec<Vec<f64>> {
        self.vertices.iter().map(|v| v.to_vec()).collect()
    }

    /// Get face indices as a flat list.
    pub fn get_faces(&self) -> Vec<usize> {
        self.faces.clone()
    }
}

impl PyConvexHull {
    /// Build the convex hull of a slice of points with the engine's Quickhull.
    pub fn from_points_internal(points: &[[f64; 3]]) -> Result<Self, ConvexHullError> {
        if points.len() < 4 {
            return Err(ConvexHullError::TooFewPoints(points.len()));
        }
        if let Some(index) = points.iter().position(|p| p.iter().any(|c| !c.is_finite())) {
            return Err(ConvexHullError::NonFinitePoint(index));
        }
        let hull = ConvexHull3D::build(points).ok_or(ConvexHullError::Degenerate)?;
        Ok(Self {
            vertices: hull.vertices,
            faces: hull.faces.into_iter().flatten().collect(),
        })
    }

    /// The same hull as an engine [`ConvexHull3D`] (for GJK).
    fn engine_hull(&self) -> ConvexHull3D {
        ConvexHull3D {
            vertices: self.vertices.clone(),
            faces: self
                .faces
                .chunks_exact(3)
                .map(|tri| [tri[0], tri[1], tri[2]])
                .collect(),
        }
    }
}
