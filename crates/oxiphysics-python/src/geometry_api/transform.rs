// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Point cloud, geometry transforms, mesh quality, and utility pyfunction wrappers.

use oxiphysics::geometry::signed_distance_field::MarchingCubes;
use pyo3::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

use super::mesh::{PyTriangleMesh, compute_aabb_internal};
use super::shapes::{PyConvexHull, add3, centroid, cross3, dot3, len3, normalize3, scale3, sub3};

// ---------------------------------------------------------------------------
// PyPointCloud
// ---------------------------------------------------------------------------

/// A point cloud with optional per-point normals.
#[pyclass(skip_from_py_object)]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyPointCloud {
    /// Point positions.
    pub points: Vec<[f64; 3]>,
    /// Per-point normals (may be empty).
    pub normals: Vec<[f64; 3]>,
}

#[pymethods]
impl PyPointCloud {
    /// Create an empty point cloud.
    #[new]
    pub fn new() -> Self {
        Self {
            points: vec![],
            normals: vec![],
        }
    }

    /// Create from point positions (list of [x, y, z] lists).
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
        Self {
            points: pts,
            normals: vec![],
        }
    }

    /// Estimate per-point normals via PCA of nearest neighbours (simplified).
    pub fn compute_normals(&mut self, _k_neighbours: usize) {
        let c = if self.points.is_empty() {
            [0.0; 3]
        } else {
            centroid(&self.points)
        };
        self.normals = self
            .points
            .iter()
            .map(|&p| normalize3(sub3(p, c)))
            .collect();
    }

    /// Voxel-grid downsampling: keep one point per voxel cell.
    pub fn simplify(&mut self, voxel_size: f64) {
        let mut grid: HashMap<(i64, i64, i64), [f64; 3]> = HashMap::new();
        for &p in &self.points {
            let key = (
                (p[0] / voxel_size).floor() as i64,
                (p[1] / voxel_size).floor() as i64,
                (p[2] / voxel_size).floor() as i64,
            );
            grid.entry(key).or_insert(p);
        }
        self.points = grid.into_values().collect();
        self.normals.clear();
    }

    /// Implicit Moving Least Squares (IMLS) surface reconstruction.
    ///
    /// Evaluates an IMLS implicit function on a uniform grid and extracts
    /// the isosurface at zero using Marching Cubes.  Grid resolution is
    /// fixed at 32³ to keep run-time predictable; call
    /// `poisson_reconstruct_res` for a configurable resolution.
    pub fn poisson_reconstruct(&self) -> PyTriangleMesh {
        self.poisson_reconstruct_res(32)
    }

    /// IMLS reconstruction at a given grid `resolution` (clamped to 4–128).
    pub fn poisson_reconstruct_res(&self, resolution: usize) -> PyTriangleMesh {
        imls_reconstruct(&self.points, &self.normals, resolution)
    }

    /// Get points as list of [x, y, z] lists.
    pub fn get_points(&self) -> Vec<Vec<f64>> {
        self.points.iter().map(|p| p.to_vec()).collect()
    }

    /// Get normals as list of [x, y, z] lists.
    pub fn get_normals(&self) -> Vec<Vec<f64>> {
        self.normals.iter().map(|n| n.to_vec()).collect()
    }

    /// Return number of points.
    pub fn len(&self) -> usize {
        self.points.len()
    }

    /// Return true if there are no points.
    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }
}

impl Default for PyPointCloud {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// PyGeometryTransform
// ---------------------------------------------------------------------------

/// A rigid-body + scale transform for geometry operations.
#[pyclass(skip_from_py_object)]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyGeometryTransform {
    /// Translation vector.
    pub translation: [f64; 3],
    /// Rotation axis (unit vector).
    pub axis: [f64; 3],
    /// Rotation angle in radians.
    #[pyo3(get, set)]
    pub angle: f64,
    /// Uniform scale factor.
    #[pyo3(get, set)]
    pub scale: f64,
}

#[pymethods]
impl PyGeometryTransform {
    /// Identity transform.
    #[staticmethod]
    pub fn identity() -> Self {
        Self {
            translation: [0.0; 3],
            axis: [0.0, 1.0, 0.0],
            angle: 0.0,
            scale: 1.0,
        }
    }

    /// Set translation from [x, y, z] list.
    pub fn translate(mut self_: PyRefMut<'_, Self>, t: Vec<f64>) -> PyGeometryTransform {
        if t.len() >= 3 {
            self_.translation = [t[0], t[1], t[2]];
        }
        self_.clone()
    }

    /// Set rotation from axis [x, y, z] list and angle.
    pub fn rotate(
        mut self_: PyRefMut<'_, Self>,
        axis: Vec<f64>,
        angle: f64,
    ) -> PyGeometryTransform {
        if axis.len() >= 3 {
            self_.axis = normalize3([axis[0], axis[1], axis[2]]);
        }
        self_.angle = angle;
        self_.clone()
    }

    /// Set scale.
    pub fn with_scale(mut self_: PyRefMut<'_, Self>, s: f64) -> PyGeometryTransform {
        self_.scale = s;
        self_.clone()
    }

    /// Apply the transform to a single point [x, y, z].
    pub fn apply_point(&self, p: Vec<f64>) -> Vec<f64> {
        if p.len() < 3 {
            return vec![0.0; 3];
        }
        let pt = [p[0], p[1], p[2]];
        let ps = scale3(pt, self.scale);
        let pr = rodrigues(ps, self.axis, self.angle);
        add3(pr, self.translation).to_vec()
    }

    /// Apply this transform to all vertices of a mesh (in place).
    pub fn apply_to_mesh(&self, mesh: &mut PyTriangleMesh) {
        for v in &mut mesh.vertices {
            let ps = scale3(*v, self.scale);
            let pr = rodrigues(ps, self.axis, self.angle);
            *v = add3(pr, self.translation);
        }
        mesh.compute_normals();
    }

    /// Get translation as list.
    pub fn get_translation(&self) -> Vec<f64> {
        self.translation.to_vec()
    }

    /// Get axis as list.
    pub fn get_axis(&self) -> Vec<f64> {
        self.axis.to_vec()
    }
}

fn rodrigues(v: [f64; 3], axis: [f64; 3], angle: f64) -> [f64; 3] {
    let cos_a = angle.cos();
    let sin_a = angle.sin();
    let k = normalize3(axis);
    add3(
        add3(scale3(v, cos_a), scale3(cross3(k, v), sin_a)),
        scale3(k, dot3(k, v) * (1.0 - cos_a)),
    )
}

// ---------------------------------------------------------------------------
// PyMeshQuality
// ---------------------------------------------------------------------------

/// Mesh quality metrics for a triangle mesh.
#[pyclass(skip_from_py_object)]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyMeshQuality {
    /// Per-triangle aspect ratios.
    #[pyo3(get)]
    pub aspect_ratios: Vec<f64>,
    /// Per-triangle signed Jacobians (det of edge matrix).
    #[pyo3(get)]
    pub jacobians: Vec<f64>,
    /// Per-triangle skewness values in [0, 1] (0 = ideal equilateral).
    #[pyo3(get)]
    pub skewness: Vec<f64>,
    /// Per-triangle areas.
    #[pyo3(get)]
    pub areas: Vec<f64>,
}

#[pymethods]
impl PyMeshQuality {
    /// Compute quality metrics for the given mesh.
    #[staticmethod]
    pub fn compute(mesh: &PyTriangleMesh) -> Self {
        let tri_count = mesh.indices.len() / 3;
        let n = mesh.vertices.len();
        let mut aspect_ratios = Vec::with_capacity(tri_count);
        let mut jacobians = Vec::with_capacity(tri_count);
        let mut skewness = Vec::with_capacity(tri_count);
        let mut areas = Vec::with_capacity(tri_count);
        for t in 0..tri_count {
            let ia = mesh.indices[t * 3];
            let ib = mesh.indices[t * 3 + 1];
            let ic = mesh.indices[t * 3 + 2];
            if ia >= n || ib >= n || ic >= n {
                aspect_ratios.push(f64::NAN);
                jacobians.push(f64::NAN);
                skewness.push(f64::NAN);
                areas.push(0.0);
                continue;
            }
            let a = mesh.vertices[ia];
            let b = mesh.vertices[ib];
            let c = mesh.vertices[ic];
            let ab = sub3(b, a);
            let ac = sub3(c, a);
            let bc = sub3(c, b);
            let la = len3(ab);
            let lb = len3(ac);
            let lc = len3(bc);
            let cr = cross3(ab, ac);
            let area = len3(cr) * 0.5;
            areas.push(area);
            let longest = la.max(lb).max(lc);
            let inradius = if la + lb + lc > 0.0 {
                2.0 * area / (la + lb + lc)
            } else {
                0.0
            };
            let ar = if inradius > 0.0 {
                longest / (2.0 * 3.0_f64.sqrt() * inradius)
            } else {
                f64::MAX
            };
            aspect_ratios.push(ar);
            let face_n = normalize3(cr);
            let jac = dot3(cross3(ab, ac), face_n);
            jacobians.push(jac);
            let ideal_angle = std::f64::consts::PI / 3.0;
            let cos_ab_ac = if la * lb > 0.0 {
                dot3(ab, ac) / (la * lb)
            } else {
                0.0
            };
            let angle_a = cos_ab_ac.clamp(-1.0, 1.0).acos();
            let skew = ((angle_a - ideal_angle) / ideal_angle).abs();
            skewness.push(skew.min(1.0));
        }
        Self {
            aspect_ratios,
            jacobians,
            skewness,
            areas,
        }
    }

    /// Return the indices of the worst `n` triangles by aspect ratio.
    pub fn worst_elements(&self, n: usize) -> Vec<usize> {
        let mut indexed: Vec<(usize, f64)> = self
            .aspect_ratios
            .iter()
            .copied()
            .enumerate()
            .filter(|(_, v)| v.is_finite())
            .collect();
        indexed.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
        indexed.into_iter().take(n).map(|(i, _)| i).collect()
    }

    /// Mean aspect ratio across all triangles.
    pub fn mean_aspect_ratio(&self) -> f64 {
        let vals: Vec<f64> = self
            .aspect_ratios
            .iter()
            .copied()
            .filter(|v| v.is_finite())
            .collect();
        if vals.is_empty() {
            return 0.0;
        }
        vals.iter().sum::<f64>() / vals.len() as f64
    }
}

// ---------------------------------------------------------------------------
// Python-facing helper functions (#[pyfunction])
// ---------------------------------------------------------------------------

/// Python-facing AABB computation from a list of `[x, y, z]` lists.
#[pyfunction]
pub fn py_compute_aabb(points: Vec<Vec<f64>>) -> (Vec<f64>, Vec<f64>) {
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
    let (mn, mx) = compute_aabb_internal(&pts);
    (mn.to_vec(), mx.to_vec())
}

/// Serialize a triangle mesh to OBJ format string (Rust API).
pub fn mesh_to_obj_string(mesh: &PyTriangleMesh) -> String {
    let mut out = String::new();
    for v in &mesh.vertices {
        out.push_str(&format!("v {} {} {}\n", v[0], v[1], v[2]));
    }
    for n in &mesh.normals {
        out.push_str(&format!("vn {} {} {}\n", n[0], n[1], n[2]));
    }
    let tri_count = mesh.indices.len() / 3;
    for t in 0..tri_count {
        let ia = mesh.indices[t * 3] + 1;
        let ib = mesh.indices[t * 3 + 1] + 1;
        let ic = mesh.indices[t * 3 + 2] + 1;
        out.push_str(&format!("f {} {} {}\n", ia, ib, ic));
    }
    out
}

/// Python-facing OBJ serializer.
#[pyfunction]
pub fn py_mesh_to_obj_string(mesh: &PyTriangleMesh) -> String {
    mesh_to_obj_string(mesh)
}

/// Parse an OBJ string into a `PyTriangleMesh`.
///
/// Supports `v` (vertex) and `f` (face) lines only.
pub fn obj_string_to_mesh(obj: &str) -> PyTriangleMesh {
    let mut verts = Vec::new();
    let mut indices = Vec::new();
    for line in obj.lines() {
        let line = line.trim();
        if let Some(rest) = line.strip_prefix("v ") {
            let parts: Vec<f64> = rest
                .split_whitespace()
                .filter_map(|s| s.parse().ok())
                .collect();
            if parts.len() >= 3 {
                verts.push([parts[0], parts[1], parts[2]]);
            }
        } else if let Some(rest) = line.strip_prefix("f ") {
            let parts: Vec<usize> = rest
                .split_whitespace()
                .filter_map(|s| {
                    s.split('/')
                        .next()
                        .and_then(|n| n.parse::<usize>().ok())
                        .map(|n| n - 1)
                })
                .collect();
            if parts.len() >= 3 {
                for i in 1..parts.len() - 1 {
                    indices.push(parts[0]);
                    indices.push(parts[i]);
                    indices.push(parts[i + 1]);
                }
            }
        }
    }
    PyTriangleMesh::from_raw_internal(verts, indices)
}

/// Python-facing OBJ parser.
#[pyfunction]
pub fn py_obj_string_to_mesh(obj: &str) -> PyTriangleMesh {
    obj_string_to_mesh(obj)
}

/// Compute the convex hull of all vertices in a mesh (Rust API).
pub fn convex_hull_from_mesh(mesh: &PyTriangleMesh) -> PyConvexHull {
    PyConvexHull::from_points_internal(&mesh.vertices)
}

/// Python-facing convex hull from mesh.
#[pyfunction]
pub fn py_convex_hull_from_mesh(mesh: &PyTriangleMesh) -> PyConvexHull {
    convex_hull_from_mesh(mesh)
}

// ---------------------------------------------------------------------------
// IMLS implicit surface reconstruction
// ---------------------------------------------------------------------------

/// Evaluate an IMLS implicit function at query point `x`.
///
/// Returns a large positive value (outside) when the weighted sum is
/// negligible, which correctly places uninformed regions outside the surface.
fn imls_eval(x: [f64; 3], points: &[[f64; 3]], normals: &[[f64; 3]], h_sq: f64) -> f64 {
    let mut sum_w = 0.0_f64;
    let mut sum_wn = 0.0_f64;
    for (p, n) in points.iter().zip(normals.iter()) {
        let dx = x[0] - p[0];
        let dy = x[1] - p[1];
        let dz = x[2] - p[2];
        let d_sq = dx * dx + dy * dy + dz * dz;
        let w = (-d_sq / h_sq).exp();
        sum_w += w;
        sum_wn += w * (dx * n[0] + dy * n[1] + dz * n[2]);
    }
    if sum_w < 1e-12 {
        return 1.0;
    }
    sum_wn / sum_w
}

/// Estimate per-point normals from a point set by computing centroid-based
/// outward directions — a fast fallback when no normals are provided.
fn estimate_normals_centroid(points: &[[f64; 3]]) -> Vec<[f64; 3]> {
    let c = centroid(points);
    points.iter().map(|&p| normalize3(sub3(p, c))).collect()
}

/// Full IMLS surface reconstruction pipeline.
///
/// 1. Validate / estimate normals.
/// 2. Compute bandwidth `h` from mean nearest-neighbour distances.
/// 3. Determine padded AABB.
/// 4. Sample the IMLS implicit on a uniform grid.
/// 5. Extract the zero isosurface with Marching Cubes.
pub(crate) fn imls_reconstruct(
    points: &[[f64; 3]],
    normals: &[[f64; 3]],
    resolution: usize,
) -> PyTriangleMesh {
    if points.is_empty() {
        return PyTriangleMesh::new();
    }

    let n_pts = points.len();

    let effective_normals: Vec<[f64; 3]> = if normals.len() == n_pts {
        normals.to_vec()
    } else {
        estimate_normals_centroid(points)
    };

    // Compute mean nearest-neighbour distance as bandwidth.
    let sample_count = n_pts.min(100);
    let mut h_sum = 0.0_f64;
    let mut h_cnt = 0_usize;
    for i in 0..sample_count {
        let p = points[i];
        let mut best_d2 = f64::MAX;
        for (j, q) in points.iter().enumerate() {
            if j == i {
                continue;
            }
            let d2 = {
                let dx = p[0] - q[0];
                let dy = p[1] - q[1];
                let dz = p[2] - q[2];
                dx * dx + dy * dy + dz * dz
            };
            if d2 < best_d2 {
                best_d2 = d2;
            }
        }
        if best_d2.is_finite() && best_d2 > 0.0 {
            h_sum += best_d2.sqrt();
            h_cnt += 1;
        }
    }
    let h = if h_cnt > 0 && h_sum / h_cnt as f64 > 1e-12 {
        h_sum / h_cnt as f64
    } else {
        0.1
    };
    let h_sq = h * h;

    // Padded AABB.
    let (mn, mx) = compute_aabb_internal(points);
    let pad = |lo: f64, hi: f64| {
        let span = hi - lo;
        let p = 0.1 * (span + h);
        (lo - p, hi + p)
    };
    let (xlo, xhi) = pad(mn[0], mx[0]);
    let (ylo, yhi) = pad(mn[1], mx[1]);
    let (zlo, zhi) = pad(mn[2], mx[2]);

    let res = resolution.clamp(4, 128);
    let bounds = [xlo, xhi, ylo, yhi, zlo, zhi];

    let pts = points.to_vec();
    let nrms = effective_normals;

    let mc_result = MarchingCubes::from_function(
        &|p: [f64; 3]| imls_eval(p, &pts, &nrms, h_sq),
        res,
        res,
        res,
        bounds,
    )
    .extract(0.0);

    let vertices: Vec<[f64; 3]> = mc_result.vertices.iter().map(|v| v.position).collect();
    let indices: Vec<usize> = mc_result
        .triangles
        .iter()
        .flat_map(|t| [t.indices[0], t.indices[1], t.indices[2]])
        .collect();

    PyTriangleMesh::from_raw_internal(vertices, indices)
}

// ---------------------------------------------------------------------------
// Tests for IMLS reconstruction
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    /// Generate ~n points uniformly on the unit sphere with outward normals.
    fn sphere_point_cloud(n: usize) -> PyPointCloud {
        let mut points = Vec::with_capacity(n);
        let mut normals = Vec::with_capacity(n);
        // Fibonacci lattice for roughly uniform sphere sampling.
        let golden = (1.0 + 5.0_f64.sqrt()) / 2.0;
        for i in 0..n {
            let theta = (1.0 - 2.0 * (i as f64 + 0.5) / n as f64).acos();
            let phi = 2.0 * PI * (i as f64) / golden;
            let x = theta.sin() * phi.cos();
            let y = theta.sin() * phi.sin();
            let z = theta.cos();
            points.push([x, y, z]);
            normals.push([x, y, z]);
        }
        PyPointCloud { points, normals }
    }

    #[test]
    fn test_imls_reconstruct_sphere_point_cloud() {
        let cloud = sphere_point_cloud(100);
        let mesh = cloud.poisson_reconstruct_res(24);
        assert!(
            !mesh.vertices.is_empty(),
            "reconstructed sphere mesh should have vertices"
        );
        assert!(
            !mesh.indices.is_empty(),
            "reconstructed sphere mesh should have triangles"
        );
        for v in &mesh.vertices {
            for &coord in v.iter() {
                assert!(
                    coord.abs() <= 1.6,
                    "reconstructed vertex {} outside [-1.6, 1.6]³",
                    coord
                );
            }
        }
    }

    #[test]
    fn test_imls_reconstruct_empty_input() {
        let cloud = PyPointCloud::new();
        let mesh = cloud.poisson_reconstruct();
        assert!(
            mesh.vertices.is_empty(),
            "empty cloud should yield empty mesh"
        );
        assert!(
            mesh.indices.is_empty(),
            "empty cloud should yield empty mesh"
        );
    }

    #[test]
    fn test_imls_reconstruct_without_normals() {
        let cloud = PyPointCloud {
            points: sphere_point_cloud(60).points,
            normals: vec![],
        };
        let mesh = cloud.poisson_reconstruct_res(20);
        assert!(
            !mesh.vertices.is_empty(),
            "mesh should be non-empty even when normals are estimated internally"
        );
    }
}
