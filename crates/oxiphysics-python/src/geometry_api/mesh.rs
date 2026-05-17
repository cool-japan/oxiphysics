// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Triangle mesh, CSG operations, and related helpers.

use pyo3::prelude::*;
use serde::{Deserialize, Serialize};

use super::shapes::{add3, cross3, dot3, len3, lerp3, normalize3, scale3, sub3};

// ---------------------------------------------------------------------------
// Helper: point-in-AABB test (pub(super) so spatial.rs can use it too)
// ---------------------------------------------------------------------------

pub(super) fn point_in_aabb(p: [f64; 3], mn: [f64; 3], mx: [f64; 3]) -> bool {
    p[0] >= mn[0]
        && p[0] <= mx[0]
        && p[1] >= mn[1]
        && p[1] <= mx[1]
        && p[2] >= mn[2]
        && p[2] <= mx[2]
}

// ---------------------------------------------------------------------------
// Helper: compute AABB (pub(super) used by transform.rs helpers too)
// ---------------------------------------------------------------------------

pub(super) fn compute_aabb_internal(points: &[[f64; 3]]) -> ([f64; 3], [f64; 3]) {
    if points.is_empty() {
        return ([0.0; 3], [0.0; 3]);
    }
    let mut mn = points[0];
    let mut mx = points[0];
    for &p in points {
        for i in 0..3 {
            if p[i] < mn[i] {
                mn[i] = p[i];
            }
            if p[i] > mx[i] {
                mx[i] = p[i];
            }
        }
    }
    (mn, mx)
}

/// Compute the axis-aligned bounding box of a point set (public Rust API).
pub fn compute_aabb(points: &[[f64; 3]]) -> ([f64; 3], [f64; 3]) {
    compute_aabb_internal(points)
}

// ---------------------------------------------------------------------------
// PyTriangleMesh
// ---------------------------------------------------------------------------

/// A triangle mesh with vertices, indices, and optional per-vertex normals.
#[pyclass(from_py_object)]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyTriangleMesh {
    /// Vertex positions, each `[f64; 3]`.
    pub vertices: Vec<[f64; 3]>,
    /// Triangle face indices (groups of 3 vertex indices).
    pub indices: Vec<usize>,
    /// Per-vertex normals (may be empty until `compute_normals` is called).
    pub normals: Vec<[f64; 3]>,
}

#[pymethods]
impl PyTriangleMesh {
    /// Create a new empty triangle mesh.
    #[new]
    pub fn new() -> Self {
        Self {
            vertices: vec![],
            indices: vec![],
            normals: vec![],
        }
    }

    /// Create a mesh from raw vertex and index data.
    /// `vertices`: list of [x, y, z] lists; `indices`: flat list of vertex indices.
    #[staticmethod]
    pub fn from_raw(vertices: Vec<Vec<f64>>, indices: Vec<usize>) -> Self {
        let verts: Vec<[f64; 3]> = vertices
            .into_iter()
            .filter_map(|v| {
                if v.len() >= 3 {
                    Some([v[0], v[1], v[2]])
                } else {
                    None
                }
            })
            .collect();
        let mut mesh = Self {
            vertices: verts,
            indices,
            normals: vec![],
        };
        mesh.compute_normals();
        mesh
    }

    /// Compute per-vertex normals as the area-weighted average of adjacent face normals.
    pub fn compute_normals(&mut self) {
        let n = self.vertices.len();
        self.normals = vec![[0.0; 3]; n];
        let tri_count = self.indices.len() / 3;
        for t in 0..tri_count {
            let ia = self.indices[t * 3];
            let ib = self.indices[t * 3 + 1];
            let ic = self.indices[t * 3 + 2];
            if ia >= n || ib >= n || ic >= n {
                continue;
            }
            let a = self.vertices[ia];
            let b = self.vertices[ib];
            let c = self.vertices[ic];
            let ab = sub3(b, a);
            let ac = sub3(c, a);
            let face_n = cross3(ab, ac);
            for &idx in &[ia, ib, ic] {
                self.normals[idx][0] += face_n[0];
                self.normals[idx][1] += face_n[1];
                self.normals[idx][2] += face_n[2];
            }
        }
        for normal in &mut self.normals {
            *normal = normalize3(*normal);
        }
    }

    /// Remove degenerate triangles (basic repair).
    pub fn repair(&mut self) {
        let mut good = Vec::new();
        let n = self.vertices.len();
        let tri_count = self.indices.len() / 3;
        for t in 0..tri_count {
            let ia = self.indices[t * 3];
            let ib = self.indices[t * 3 + 1];
            let ic = self.indices[t * 3 + 2];
            if ia >= n || ib >= n || ic >= n {
                continue;
            }
            if ia == ib || ib == ic || ia == ic {
                continue;
            }
            let a = self.vertices[ia];
            let b = self.vertices[ib];
            let c = self.vertices[ic];
            let area = len3(cross3(sub3(b, a), sub3(c, a)));
            if area > 1e-15 {
                good.extend_from_slice(&[ia, ib, ic]);
            }
        }
        self.indices = good;
    }

    /// Laplacian smoothing: average each vertex toward its neighbours.
    pub fn smooth(&mut self, iterations: usize) {
        for _ in 0..iterations {
            let n = self.vertices.len();
            let mut acc = vec![[0.0f64; 3]; n];
            let mut cnt = vec![0usize; n];
            let tri_count = self.indices.len() / 3;
            for t in 0..tri_count {
                let ia = self.indices[t * 3];
                let ib = self.indices[t * 3 + 1];
                let ic = self.indices[t * 3 + 2];
                if ia >= n || ib >= n || ic >= n {
                    continue;
                }
                let va = self.vertices[ia];
                let vb = self.vertices[ib];
                let vc = self.vertices[ic];
                acc[ia] = add3(acc[ia], add3(vb, vc));
                cnt[ia] += 2;
                acc[ib] = add3(acc[ib], add3(va, vc));
                cnt[ib] += 2;
                acc[ic] = add3(acc[ic], add3(va, vb));
                cnt[ic] += 2;
            }
            for i in 0..n {
                if cnt[i] > 0 {
                    let c = cnt[i] as f64;
                    self.vertices[i] = lerp3(
                        self.vertices[i],
                        [acc[i][0] / c, acc[i][1] / c, acc[i][2] / c],
                        0.5,
                    );
                }
            }
        }
        self.compute_normals();
    }

    /// Compute signed volume using the divergence theorem.
    pub fn compute_volume(&self) -> f64 {
        let tri_count = self.indices.len() / 3;
        let n = self.vertices.len();
        let mut vol = 0.0;
        for t in 0..tri_count {
            let ia = self.indices[t * 3];
            let ib = self.indices[t * 3 + 1];
            let ic = self.indices[t * 3 + 2];
            if ia >= n || ib >= n || ic >= n {
                continue;
            }
            let a = self.vertices[ia];
            let b = self.vertices[ib];
            let c = self.vertices[ic];
            vol += dot3(a, cross3(b, c));
        }
        (vol / 6.0).abs()
    }

    /// Compute total surface area.
    pub fn compute_surface_area(&self) -> f64 {
        let tri_count = self.indices.len() / 3;
        let n = self.vertices.len();
        let mut area = 0.0;
        for t in 0..tri_count {
            let ia = self.indices[t * 3];
            let ib = self.indices[t * 3 + 1];
            let ic = self.indices[t * 3 + 2];
            if ia >= n || ib >= n || ic >= n {
                continue;
            }
            let a = self.vertices[ia];
            let b = self.vertices[ib];
            let c = self.vertices[ic];
            area += len3(cross3(sub3(b, a), sub3(c, a))) * 0.5;
        }
        area
    }

    /// Number of triangles.
    pub fn triangle_count(&self) -> usize {
        self.indices.len() / 3
    }

    /// Get vertices as a list of [x, y, z] lists.
    pub fn get_vertices(&self) -> Vec<Vec<f64>> {
        self.vertices.iter().map(|v| v.to_vec()).collect()
    }

    /// Get face indices as a flat list.
    pub fn get_indices(&self) -> Vec<usize> {
        self.indices.clone()
    }

    /// Get normals as a list of [nx, ny, nz] lists.
    pub fn get_normals(&self) -> Vec<Vec<f64>> {
        self.normals.iter().map(|nrm| nrm.to_vec()).collect()
    }
}

impl Default for PyTriangleMesh {
    fn default() -> Self {
        Self::new()
    }
}

impl PyTriangleMesh {
    /// Internal constructor from already-typed data.
    pub fn from_raw_internal(vertices: Vec<[f64; 3]>, indices: Vec<usize>) -> Self {
        let mut mesh = Self {
            vertices,
            indices,
            normals: vec![],
        };
        mesh.compute_normals();
        mesh
    }
}

// ---------------------------------------------------------------------------
// PyCsg — Constructive Solid Geometry stubs
// ---------------------------------------------------------------------------

/// Result type for CSG operations between two meshes.
///
/// Note: Full BSP-based CSG is complex; these are structural stubs that
/// merge/intersect vertex data in a simplified manner for API completeness.
#[pyclass(skip_from_py_object)]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyCsg;

#[pymethods]
impl PyCsg {
    /// Union of two meshes (concatenates geometry — placeholder for full BSP).
    #[staticmethod]
    pub fn union(a: &PyTriangleMesh, b: &PyTriangleMesh) -> PyTriangleMesh {
        let offset = a.vertices.len();
        let mut verts = a.vertices.clone();
        verts.extend_from_slice(&b.vertices);
        let mut idx = a.indices.clone();
        for &i in &b.indices {
            idx.push(i + offset);
        }
        PyTriangleMesh::from_raw_internal(verts, idx)
    }

    /// Intersection stub: returns mesh `a` clipped by the AABB of mesh `b`.
    #[staticmethod]
    pub fn intersection(a: &PyTriangleMesh, b: &PyTriangleMesh) -> PyTriangleMesh {
        let (bmin, bmax) = compute_aabb_internal(&b.vertices);
        let mut verts = Vec::new();
        let mut idx = Vec::new();
        let tri_count = a.indices.len() / 3;
        let n = a.vertices.len();
        let mut remap = vec![usize::MAX; n];
        for t in 0..tri_count {
            let ia = a.indices[t * 3];
            let ib = a.indices[t * 3 + 1];
            let ic = a.indices[t * 3 + 2];
            if ia >= n || ib >= n || ic >= n {
                continue;
            }
            let ctr = scale3(
                add3(a.vertices[ia], add3(a.vertices[ib], a.vertices[ic])),
                1.0 / 3.0,
            );
            if point_in_aabb(ctr, bmin, bmax) {
                for &vi in &[ia, ib, ic] {
                    if remap[vi] == usize::MAX {
                        remap[vi] = verts.len();
                        verts.push(a.vertices[vi]);
                    }
                    idx.push(remap[vi]);
                }
            }
        }
        PyTriangleMesh::from_raw_internal(verts, idx)
    }

    /// Subtraction stub: removes triangles whose centroid is inside mesh `b`'s AABB.
    #[staticmethod]
    pub fn subtraction(a: &PyTriangleMesh, b: &PyTriangleMesh) -> PyTriangleMesh {
        let (bmin, bmax) = compute_aabb_internal(&b.vertices);
        let mut verts = Vec::new();
        let mut idx = Vec::new();
        let tri_count = a.indices.len() / 3;
        let n = a.vertices.len();
        let mut remap = vec![usize::MAX; n];
        for t in 0..tri_count {
            let ia = a.indices[t * 3];
            let ib = a.indices[t * 3 + 1];
            let ic = a.indices[t * 3 + 2];
            if ia >= n || ib >= n || ic >= n {
                continue;
            }
            let ctr = scale3(
                add3(a.vertices[ia], add3(a.vertices[ib], a.vertices[ic])),
                1.0 / 3.0,
            );
            if !point_in_aabb(ctr, bmin, bmax) {
                for &vi in &[ia, ib, ic] {
                    if remap[vi] == usize::MAX {
                        remap[vi] = verts.len();
                        verts.push(a.vertices[vi]);
                    }
                    idx.push(remap[vi]);
                }
            }
        }
        PyTriangleMesh::from_raw_internal(verts, idx)
    }
}
