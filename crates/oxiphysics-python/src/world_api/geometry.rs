// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Geometry Queries: AABB, Sphere, ConvexHull.

use pyo3::prelude::*;

// ===========================================================================
// Geometry Queries
// ===========================================================================

/// Axis-aligned bounding box (AABB) in 3-D.
#[pyclass(from_py_object)]
#[derive(Debug, Clone, Copy)]
pub struct PyAabb {
    /// Minimum corner `[xmin, ymin, zmin]`.
    #[pyo3(get, set)]
    pub min: [f64; 3],
    /// Maximum corner `[xmax, ymax, zmax]`.
    #[pyo3(get, set)]
    pub max: [f64; 3],
}

#[pymethods]
impl PyAabb {
    /// Create from min/max corner points.
    #[new]
    pub fn new(min: [f64; 3], max: [f64; 3]) -> Self {
        Self { min, max }
    }

    /// Create a unit cube centred at the origin.
    #[staticmethod]
    pub fn unit() -> Self {
        Self::new([-0.5; 3], [0.5; 3])
    }

    /// Create from centre and half-extents.
    #[staticmethod]
    pub fn from_center_half_extents(center: [f64; 3], he: [f64; 3]) -> Self {
        Self {
            min: [center[0] - he[0], center[1] - he[1], center[2] - he[2]],
            max: [center[0] + he[0], center[1] + he[1], center[2] + he[2]],
        }
    }

    /// Centre point of the AABB.
    pub fn center(&self) -> [f64; 3] {
        [
            (self.min[0] + self.max[0]) * 0.5,
            (self.min[1] + self.max[1]) * 0.5,
            (self.min[2] + self.max[2]) * 0.5,
        ]
    }

    /// Half-extents in each dimension.
    pub fn half_extents(&self) -> [f64; 3] {
        [
            (self.max[0] - self.min[0]) * 0.5,
            (self.max[1] - self.min[1]) * 0.5,
            (self.max[2] - self.min[2]) * 0.5,
        ]
    }

    /// Surface area of the AABB.
    pub fn surface_area(&self) -> f64 {
        let dx = self.max[0] - self.min[0];
        let dy = self.max[1] - self.min[1];
        let dz = self.max[2] - self.min[2];
        2.0 * (dx * dy + dy * dz + dz * dx)
    }

    /// Volume of the AABB.
    pub fn volume(&self) -> f64 {
        let dx = (self.max[0] - self.min[0]).max(0.0);
        let dy = (self.max[1] - self.min[1]).max(0.0);
        let dz = (self.max[2] - self.min[2]).max(0.0);
        dx * dy * dz
    }

    /// Whether a point is contained within (or on the boundary of) this AABB.
    pub fn contains_point(&self, p: [f64; 3]) -> bool {
        p[0] >= self.min[0]
            && p[0] <= self.max[0]
            && p[1] >= self.min[1]
            && p[1] <= self.max[1]
            && p[2] >= self.min[2]
            && p[2] <= self.max[2]
    }

    /// Whether this AABB intersects another AABB.
    pub fn intersects(&self, other: &PyAabb) -> bool {
        self.min[0] <= other.max[0]
            && self.max[0] >= other.min[0]
            && self.min[1] <= other.max[1]
            && self.max[1] >= other.min[1]
            && self.min[2] <= other.max[2]
            && self.max[2] >= other.min[2]
    }

    /// Merge with another AABB, producing the smallest enclosing AABB.
    pub fn merged(&self, other: &PyAabb) -> PyAabb {
        PyAabb {
            min: [
                self.min[0].min(other.min[0]),
                self.min[1].min(other.min[1]),
                self.min[2].min(other.min[2]),
            ],
            max: [
                self.max[0].max(other.max[0]),
                self.max[1].max(other.max[1]),
                self.max[2].max(other.max[2]),
            ],
        }
    }
}

/// Sphere geometry query helper.
#[pyclass(from_py_object)]
#[derive(Debug, Clone, Copy)]
pub struct PySphere {
    /// Centre of the sphere.
    #[pyo3(get, set)]
    pub center: [f64; 3],
    /// Radius.
    #[pyo3(get, set)]
    pub radius: f64,
}

#[pymethods]
impl PySphere {
    /// Create a sphere from centre and radius.
    #[new]
    pub fn new(center: [f64; 3], radius: f64) -> Self {
        Self {
            center,
            radius: radius.max(0.0),
        }
    }

    /// Unit sphere at origin.
    #[staticmethod]
    pub fn unit() -> Self {
        Self::new([0.0; 3], 1.0)
    }

    /// Surface area 4π r².
    pub fn surface_area(&self) -> f64 {
        4.0 * std::f64::consts::PI * self.radius * self.radius
    }

    /// Volume 4/3 π r³.
    pub fn volume(&self) -> f64 {
        (4.0 / 3.0) * std::f64::consts::PI * self.radius * self.radius * self.radius
    }

    /// Whether a point is inside (or on the surface of) this sphere.
    pub fn contains_point(&self, p: [f64; 3]) -> bool {
        let dx = p[0] - self.center[0];
        let dy = p[1] - self.center[1];
        let dz = p[2] - self.center[2];
        dx * dx + dy * dy + dz * dz <= self.radius * self.radius
    }

    /// Signed distance from the sphere surface to point `p`.
    ///
    /// Negative inside, positive outside.
    pub fn signed_distance(&self, p: [f64; 3]) -> f64 {
        let dx = p[0] - self.center[0];
        let dy = p[1] - self.center[1];
        let dz = p[2] - self.center[2];
        (dx * dx + dy * dy + dz * dz).sqrt() - self.radius
    }

    /// Whether this sphere overlaps another sphere.
    pub fn overlaps(&self, other: &PySphere) -> bool {
        let dx = self.center[0] - other.center[0];
        let dy = self.center[1] - other.center[1];
        let dz = self.center[2] - other.center[2];
        let dist2 = dx * dx + dy * dy + dz * dz;
        let sum_r = self.radius + other.radius;
        dist2 <= sum_r * sum_r
    }

    /// Bounding AABB of this sphere.
    pub fn aabb(&self) -> PyAabb {
        PyAabb::from_center_half_extents(self.center, [self.radius; 3])
    }
}

/// A convex hull stored as a list of vertices.
///
/// The hull is not computed internally; the caller is responsible for
/// providing convex vertices. Methods here are geometry helpers.
#[pyclass(from_py_object)]
#[derive(Debug, Clone)]
pub struct PyConvexHull {
    /// Vertices of the convex hull.
    #[pyo3(get, set)]
    pub vertices: Vec<[f64; 3]>,
}

#[pymethods]
impl PyConvexHull {
    /// Create a convex hull from a list of vertices.
    #[new]
    pub fn new(vertices: Vec<[f64; 3]>) -> Self {
        Self { vertices }
    }

    /// Create a convex hull approximating a unit cube.
    #[staticmethod]
    pub fn unit_cube() -> Self {
        let verts: Vec<[f64; 3]> = [
            [-0.5, -0.5, -0.5],
            [0.5, -0.5, -0.5],
            [0.5, 0.5, -0.5],
            [-0.5, 0.5, -0.5],
            [-0.5, -0.5, 0.5],
            [0.5, -0.5, 0.5],
            [0.5, 0.5, 0.5],
            [-0.5, 0.5, 0.5],
        ]
        .to_vec();
        Self::new(verts)
    }

    /// Number of vertices in the hull.
    pub fn vertex_count(&self) -> usize {
        self.vertices.len()
    }

    /// Compute the axis-aligned bounding box of all vertices.
    pub fn aabb(&self) -> Option<PyAabb> {
        if self.vertices.is_empty() {
            return None;
        }
        let mut mn = self.vertices[0];
        let mut mx = self.vertices[0];
        for &v in &self.vertices {
            for k in 0..3 {
                mn[k] = mn[k].min(v[k]);
                mx[k] = mx[k].max(v[k]);
            }
        }
        Some(PyAabb::new(mn, mx))
    }

    /// Centroid (arithmetic mean) of all vertices.
    pub fn centroid(&self) -> Option<[f64; 3]> {
        if self.vertices.is_empty() {
            return None;
        }
        let n = self.vertices.len() as f64;
        let mut c = [0.0f64; 3];
        for v in &self.vertices {
            c[0] += v[0];
            c[1] += v[1];
            c[2] += v[2];
        }
        Some([c[0] / n, c[1] / n, c[2] / n])
    }

    /// Naively check if a point is "inside" the hull by checking it is
    /// within the AABB (a conservative, approximate test).
    pub fn may_contain_point(&self, p: [f64; 3]) -> bool {
        match self.aabb() {
            Some(aabb) => aabb.contains_point(p),
            None => false,
        }
    }

    /// Support function: find vertex furthest along direction `d`.
    pub fn support(&self, d: [f64; 3]) -> Option<[f64; 3]> {
        self.vertices
            .iter()
            .max_by(|a, b| {
                let da = a[0] * d[0] + a[1] * d[1] + a[2] * d[2];
                let db = b[0] * d[0] + b[1] * d[1] + b[2] * d[2];
                da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
            })
            .copied()
    }
}
