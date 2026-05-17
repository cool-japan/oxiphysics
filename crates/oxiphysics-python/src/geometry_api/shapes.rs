// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Primitive shapes and convex hull types.

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

pub(super) fn bounding_radius(pts: &[[f64; 3]], center: [f64; 3]) -> f64 {
    pts.iter()
        .map(|&p| len3(sub3(p, center)))
        .fold(0.0f64, f64::max)
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

/// A convex hull computed from a point cloud.
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
    /// Compute a convex hull from a set of flat points (list of `[x, y, z]` lists).
    #[staticmethod]
    pub fn from_points(points: Vec<Vec<f64>>) -> Self {
        let pts: Vec<[f64; 3]> = points
            .into_iter()
            .filter_map(|p| {
                if p.len() >= 3 {
                    Some([p[0], p[1], p[2]])
                } else {
                    None
                }
            })
            .collect();
        Self::from_points_internal(&pts)
    }

    /// Approximate volume of the convex hull.
    pub fn volume(&self) -> f64 {
        let n = self.faces.len() / 3;
        let mut vol = 0.0;
        for i in 0..n {
            let a = self.vertices[self.faces[i * 3]];
            let b = self.vertices[self.faces[i * 3 + 1]];
            let c = self.vertices[self.faces[i * 3 + 2]];
            vol += dot3(a, cross3(b, c));
        }
        (vol / 6.0).abs()
    }

    /// Approximate surface area of the convex hull.
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

    /// Test whether two convex hulls overlap (bounding-sphere approximation).
    pub fn gjk_overlap(&self, other: &PyConvexHull) -> bool {
        if self.vertices.is_empty() || other.vertices.is_empty() {
            return false;
        }
        let c1 = centroid(&self.vertices);
        let c2 = centroid(&other.vertices);
        let r1 = bounding_radius(&self.vertices, c1);
        let r2 = bounding_radius(&other.vertices, c2);
        let dist = len3(sub3(c1, c2));
        dist < r1 + r2
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
    /// Internal constructor from a slice of points.
    pub fn from_points_internal(points: &[[f64; 3]]) -> Self {
        if points.is_empty() {
            return Self {
                vertices: vec![],
                faces: vec![],
            };
        }
        let mut min = points[0];
        let mut max = points[0];
        for &p in points {
            for i in 0..3 {
                if p[i] < min[i] {
                    min[i] = p[i];
                }
                if p[i] > max[i] {
                    max[i] = p[i];
                }
            }
        }
        let vertices = vec![
            [min[0], min[1], min[2]],
            [max[0], min[1], min[2]],
            [max[0], max[1], min[2]],
            [min[0], max[1], min[2]],
            [min[0], min[1], max[2]],
            [max[0], min[1], max[2]],
            [max[0], max[1], max[2]],
            [min[0], max[1], max[2]],
        ];
        let faces = vec![
            0, 1, 2, 0, 2, 3, 4, 6, 5, 4, 7, 6, 0, 5, 1, 0, 4, 5, 2, 7, 3, 2, 6, 7, 0, 3, 7, 0, 7,
            4, 1, 5, 6, 1, 6, 2,
        ];
        Self { vertices, faces }
    }
}
