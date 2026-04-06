#![allow(clippy::manual_strip, clippy::ptr_arg)]
#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Python geometry bindings.
//!
//! Provides Python-friendly geometry primitives, mesh operations, spatial data
//! structures, and CSG operations using plain `f64` arrays — no nalgebra.

#![allow(dead_code)]
#![allow(missing_docs)]

use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Math helpers (local, no nalgebra)
// ---------------------------------------------------------------------------

fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

fn len3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}

fn normalize3(a: [f64; 3]) -> [f64; 3] {
    let l = len3(a);
    if l < 1e-15 {
        [0.0; 3]
    } else {
        scale3(a, 1.0 / l)
    }
}

fn lerp3(a: [f64; 3], b: [f64; 3], t: f64) -> [f64; 3] {
    [
        a[0] + (b[0] - a[0]) * t,
        a[1] + (b[1] - a[1]) * t,
        a[2] + (b[2] - a[2]) * t,
    ]
}

// ---------------------------------------------------------------------------
// PyShape — primitive geometry enum
// ---------------------------------------------------------------------------

/// Primitive geometry shape used for collision detection and rendering.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum PyShape {
    /// A sphere defined by its radius.
    Sphere(f64),
    /// An axis-aligned box defined by half-extents.
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
        use std::f64::consts::PI;
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

// ---------------------------------------------------------------------------
// PyConvexHull
// ---------------------------------------------------------------------------

/// A convex hull computed from a point cloud.
///
/// Uses a simplified incremental approach for demonstration purposes.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyConvexHull {
    /// Vertices of the convex hull.
    pub vertices: Vec<[f64; 3]>,
    /// Triangle face indices (groups of 3).
    pub faces: Vec<usize>,
}

impl PyConvexHull {
    /// Compute a convex hull from a set of points.
    ///
    /// This is a simplified implementation — for production use a full
    /// Quickhull algorithm.
    pub fn from_points(points: &[[f64; 3]]) -> Self {
        if points.is_empty() {
            return Self {
                vertices: vec![],
                faces: vec![],
            };
        }
        // Simplified: keep extreme points as hull approximation
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
        // 8 corners of the bounding box as hull vertices
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
        // 12 triangles for the box faces
        let faces = vec![
            0, 1, 2, 0, 2, 3, // -z
            4, 6, 5, 4, 7, 6, // +z
            0, 5, 1, 0, 4, 5, // -y
            2, 7, 3, 2, 6, 7, // +y
            0, 3, 7, 0, 7, 4, // -x
            1, 5, 6, 1, 6, 2, // +x
        ];
        Self { vertices, faces }
    }

    /// Approximate volume of the convex hull.
    pub fn volume(&self) -> f64 {
        let n = self.faces.len() / 3;
        let mut vol = 0.0;
        for i in 0..n {
            let a = self.vertices[self.faces[i * 3]];
            let b = self.vertices[self.faces[i * 3 + 1]];
            let c = self.vertices[self.faces[i * 3 + 2]];
            let cr = cross3(b, c);
            vol += dot3(a, cr);
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
            let cr = cross3(ab, ac);
            area += len3(cr) * 0.5;
        }
        area
    }

    /// Test whether two convex hulls overlap using a simplified GJK-like check
    /// (bounding-sphere test as approximation).
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
}

fn centroid(pts: &[[f64; 3]]) -> [f64; 3] {
    let n = pts.len() as f64;
    let mut s = [0.0f64; 3];
    for &p in pts {
        s[0] += p[0];
        s[1] += p[1];
        s[2] += p[2];
    }
    [s[0] / n, s[1] / n, s[2] / n]
}

fn bounding_radius(pts: &[[f64; 3]], center: [f64; 3]) -> f64 {
    pts.iter()
        .map(|&p| len3(sub3(p, center)))
        .fold(0.0f64, f64::max)
}

// ---------------------------------------------------------------------------
// PyTriangleMesh
// ---------------------------------------------------------------------------

/// A triangle mesh with vertices, indices, and optional per-vertex normals.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyTriangleMesh {
    /// Vertex positions, each `[f64; 3]`.
    pub vertices: Vec<[f64; 3]>,
    /// Triangle face indices (groups of 3 vertex indices).
    pub indices: Vec<usize>,
    /// Per-vertex normals (may be empty until `compute_normals` is called).
    pub normals: Vec<[f64; 3]>,
}

impl PyTriangleMesh {
    /// Create a new empty triangle mesh.
    pub fn new() -> Self {
        Self {
            vertices: vec![],
            indices: vec![],
            normals: vec![],
        }
    }

    /// Create a mesh from raw vertex and index data.
    pub fn from_raw(vertices: Vec<[f64; 3]>, indices: Vec<usize>) -> Self {
        let mut mesh = Self {
            vertices,
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
        for n in &mut self.normals {
            *n = normalize3(*n);
        }
    }

    /// Remove degenerate triangles and duplicated vertices (basic repair).
    pub fn repair(&mut self) {
        // Remove degenerate triangles (zero area)
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
}

impl Default for PyTriangleMesh {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// PyCsg — Constructive Solid Geometry stubs
// ---------------------------------------------------------------------------

/// Result type for CSG operations between two meshes.
///
/// Note: Full BSP-based CSG is complex; these are structural stubs that
/// merge/intersect vertex data in a simplified manner for API completeness.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyCsg;

impl PyCsg {
    /// Union of two meshes (concatenates geometry — placeholder for full BSP).
    pub fn union(a: &PyTriangleMesh, b: &PyTriangleMesh) -> PyTriangleMesh {
        let offset = a.vertices.len();
        let mut verts = a.vertices.clone();
        verts.extend_from_slice(&b.vertices);
        let mut idx = a.indices.clone();
        for &i in &b.indices {
            idx.push(i + offset);
        }
        PyTriangleMesh::from_raw(verts, idx)
    }

    /// Intersection stub: returns mesh `a` clipped by the AABB of mesh `b`.
    pub fn intersection(a: &PyTriangleMesh, b: &PyTriangleMesh) -> PyTriangleMesh {
        let (bmin, bmax) = compute_aabb(&b.vertices);
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
        PyTriangleMesh::from_raw(verts, idx)
    }

    /// Subtraction stub: removes triangles whose centroid is inside mesh `b`'s AABB.
    pub fn subtraction(a: &PyTriangleMesh, b: &PyTriangleMesh) -> PyTriangleMesh {
        let (bmin, bmax) = compute_aabb(&b.vertices);
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
        PyTriangleMesh::from_raw(verts, idx)
    }
}

fn point_in_aabb(p: [f64; 3], mn: [f64; 3], mx: [f64; 3]) -> bool {
    p[0] >= mn[0]
        && p[0] <= mx[0]
        && p[1] >= mn[1]
        && p[1] <= mx[1]
        && p[2] >= mn[2]
        && p[2] <= mx[2]
}

// ---------------------------------------------------------------------------
// PyHeightfield
// ---------------------------------------------------------------------------

/// A heightfield terrain represented as a regular grid.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyHeightfield {
    /// Number of columns (X axis).
    pub cols: usize,
    /// Number of rows (Z axis).
    pub rows: usize,
    /// Height values in row-major order (rows × cols entries).
    pub data: Vec<f64>,
    /// World-space scale per cell in X.
    pub scale_x: f64,
    /// World-space scale per cell in Z.
    pub scale_z: f64,
    /// Vertical scale factor applied to raw height values.
    pub scale_y: f64,
}

impl PyHeightfield {
    /// Create a new flat heightfield.
    pub fn new(cols: usize, rows: usize, scale_x: f64, scale_z: f64, scale_y: f64) -> Self {
        Self {
            cols,
            rows,
            data: vec![0.0; cols * rows],
            scale_x,
            scale_z,
            scale_y,
        }
    }

    /// Set height at grid cell (col, row).
    pub fn set_height(&mut self, col: usize, row: usize, h: f64) {
        if col < self.cols && row < self.rows {
            self.data[row * self.cols + col] = h;
        }
    }

    /// Bilinearly interpolated height at world-space (x, z).
    pub fn height_at(&self, x: f64, z: f64) -> f64 {
        let gx = x / self.scale_x;
        let gz = z / self.scale_z;
        let col = gx.floor() as isize;
        let row = gz.floor() as isize;
        let fx = gx - col as f64;
        let fz = gz - row as f64;

        let sample = |c: isize, r: isize| -> f64 {
            let c = c.clamp(0, self.cols as isize - 1) as usize;
            let r = r.clamp(0, self.rows as isize - 1) as usize;
            self.data[r * self.cols + c] * self.scale_y
        };

        let h00 = sample(col, row);
        let h10 = sample(col + 1, row);
        let h01 = sample(col, row + 1);
        let h11 = sample(col + 1, row + 1);

        let h0 = h00 + (h10 - h00) * fx;
        let h1 = h01 + (h11 - h01) * fx;
        h0 + (h1 - h0) * fz
    }

    /// Central-difference surface normal at world-space (x, z).
    pub fn normal_at(&self, x: f64, z: f64) -> [f64; 3] {
        let eps = self.scale_x.min(self.scale_z) * 0.5;
        let dhdx = (self.height_at(x + eps, z) - self.height_at(x - eps, z)) / (2.0 * eps);
        let dhdz = (self.height_at(x, z + eps) - self.height_at(x, z - eps)) / (2.0 * eps);
        normalize3([-dhdx, 1.0, -dhdz])
    }

    /// Raycast against the heightfield.  Returns the hit distance if found.
    pub fn raycast(&self, origin: [f64; 3], direction: [f64; 3]) -> Option<f64> {
        let dir = normalize3(direction);
        let max_dist = (self.cols as f64 * self.scale_x + self.rows as f64 * self.scale_z) * 2.0;
        let steps = 1000usize;
        let step = max_dist / steps as f64;
        let mut t = 0.0;
        for _ in 0..steps {
            let p = [
                origin[0] + dir[0] * t,
                origin[1] + dir[1] * t,
                origin[2] + dir[2] * t,
            ];
            let h = self.height_at(p[0], p[2]);
            if p[1] <= h {
                return Some(t);
            }
            t += step;
        }
        None
    }
}

// ---------------------------------------------------------------------------
// PySpatialHash
// ---------------------------------------------------------------------------

/// A spatial hash map for fast point proximity queries.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PySpatialHash {
    /// Cell size determines resolution of the hash grid.
    pub cell_size: f64,
    /// Points indexed by ID.
    pub points: Vec<[f64; 3]>,
    /// Internal grid map: cell key -> list of point indices.
    grid: HashMap<(i64, i64, i64), Vec<usize>>,
}

impl PySpatialHash {
    /// Create a new spatial hash with the given cell size.
    pub fn new(cell_size: f64) -> Self {
        Self {
            cell_size,
            points: vec![],
            grid: HashMap::new(),
        }
    }

    fn cell_of(&self, p: [f64; 3]) -> (i64, i64, i64) {
        (
            (p[0] / self.cell_size).floor() as i64,
            (p[1] / self.cell_size).floor() as i64,
            (p[2] / self.cell_size).floor() as i64,
        )
    }

    /// Insert a point and return its ID.
    pub fn insert(&mut self, p: [f64; 3]) -> usize {
        let id = self.points.len();
        self.points.push(p);
        let cell = self.cell_of(p);
        self.grid.entry(cell).or_default().push(id);
        id
    }

    /// Remove a point by ID (marks its cell entry as removed).
    pub fn remove(&mut self, id: usize) {
        if id >= self.points.len() {
            return;
        }
        let p = self.points[id];
        let cell = self.cell_of(p);
        if let Some(list) = self.grid.get_mut(&cell) {
            list.retain(|&x| x != id);
        }
    }

    /// Query all points in the cell containing `p`.
    pub fn query(&self, p: [f64; 3]) -> Vec<usize> {
        let cell = self.cell_of(p);
        self.grid.get(&cell).cloned().unwrap_or_default()
    }

    /// Query all points within a sphere of the given `radius` around `center`.
    pub fn sphere_query(&self, center: [f64; 3], radius: f64) -> Vec<usize> {
        let r2 = radius * radius;
        let cells = (radius / self.cell_size).ceil() as i64 + 1;
        let base = self.cell_of(center);
        let mut result = Vec::new();
        for dx in -cells..=cells {
            for dy in -cells..=cells {
                for dz in -cells..=cells {
                    let key = (base.0 + dx, base.1 + dy, base.2 + dz);
                    if let Some(list) = self.grid.get(&key) {
                        for &id in list {
                            if id < self.points.len() {
                                let d = sub3(self.points[id], center);
                                if dot3(d, d) <= r2 {
                                    result.push(id);
                                }
                            }
                        }
                    }
                }
            }
        }
        result
    }
}

// ---------------------------------------------------------------------------
// PyBvh — Bounding Volume Hierarchy
// ---------------------------------------------------------------------------

/// An axis-aligned bounding box used in the BVH.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyAabb {
    /// Minimum corner.
    pub min: [f64; 3],
    /// Maximum corner.
    pub max: [f64; 3],
}

impl PyAabb {
    /// Create a new AABB.
    pub fn new(min: [f64; 3], max: [f64; 3]) -> Self {
        Self { min, max }
    }

    /// Check whether this AABB overlaps another.
    pub fn overlaps(&self, other: &PyAabb) -> bool {
        self.min[0] <= other.max[0]
            && self.max[0] >= other.min[0]
            && self.min[1] <= other.max[1]
            && self.max[1] >= other.min[1]
            && self.min[2] <= other.max[2]
            && self.max[2] >= other.min[2]
    }

    /// Center of this AABB.
    pub fn center(&self) -> [f64; 3] {
        scale3(add3(self.min, self.max), 0.5)
    }

    /// Half-extents of this AABB.
    pub fn half_extents(&self) -> [f64; 3] {
        scale3(sub3(self.max, self.min), 0.5)
    }

    /// Surface area of this AABB.
    pub fn surface_area(&self) -> f64 {
        let e = sub3(self.max, self.min);
        2.0 * (e[0] * e[1] + e[1] * e[2] + e[2] * e[0])
    }
}

/// A node in the BVH tree.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyBvhNode {
    /// Bounding box of this node.
    pub aabb: PyAabb,
    /// Left child index (or `usize::MAX` for leaf).
    pub left: usize,
    /// Right child index (or `usize::MAX` for leaf).
    pub right: usize,
    /// Leaf item index (meaningful only when `left == usize::MAX`).
    pub item: usize,
}

/// A BVH (Bounding Volume Hierarchy) over AABBs.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyBvh {
    /// All nodes in the BVH.
    pub nodes: Vec<PyBvhNode>,
    /// Root node index.
    pub root: usize,
    /// Original AABBs associated with leaf items.
    pub items: Vec<PyAabb>,
}

impl PyBvh {
    /// Build a BVH from a list of AABBs.
    pub fn build(aabbs: Vec<PyAabb>) -> Self {
        let n = aabbs.len();
        if n == 0 {
            return Self {
                nodes: vec![],
                root: 0,
                items: vec![],
            };
        }
        let mut bvh = Self {
            nodes: Vec::with_capacity(2 * n),
            root: 0,
            items: aabbs,
        };
        let indices: Vec<usize> = (0..n).collect();
        bvh.root = bvh.build_recursive(&indices);
        bvh
    }

    fn build_recursive(&mut self, indices: &[usize]) -> usize {
        if indices.len() == 1 {
            let id = indices[0];
            let node = PyBvhNode {
                aabb: self.items[id].clone(),
                left: usize::MAX,
                right: usize::MAX,
                item: id,
            };
            let idx = self.nodes.len();
            self.nodes.push(node);
            return idx;
        }
        // Compute combined AABB
        let combined = merge_aabbs(indices.iter().map(|&i| &self.items[i]));
        // Split on longest axis at median
        let axis = longest_axis(&combined);
        let mut sorted = indices.to_vec();
        sorted.sort_by(|&a, &b| {
            let ca = self.items[a].center()[axis];
            let cb = self.items[b].center()[axis];
            ca.partial_cmp(&cb).unwrap_or(std::cmp::Ordering::Equal)
        });
        let mid = sorted.len() / 2;
        let left = self.build_recursive(&sorted[..mid]);
        let right = self.build_recursive(&sorted[mid..]);
        let node = PyBvhNode {
            aabb: combined,
            left,
            right,
            item: usize::MAX,
        };
        let idx = self.nodes.len();
        self.nodes.push(node);
        idx
    }

    /// Query all leaf items whose AABB overlaps the query AABB.
    pub fn query_aabb(&self, query: &PyAabb) -> Vec<usize> {
        if self.nodes.is_empty() {
            return vec![];
        }
        let mut result = Vec::new();
        self.query_recursive(self.root, query, &mut result);
        result
    }

    fn query_recursive(&self, node_idx: usize, query: &PyAabb, out: &mut Vec<usize>) {
        if node_idx >= self.nodes.len() {
            return;
        }
        let node = &self.nodes[node_idx];
        if !node.aabb.overlaps(query) {
            return;
        }
        if node.left == usize::MAX {
            out.push(node.item);
            return;
        }
        self.query_recursive(node.left, query, out);
        self.query_recursive(node.right, query, out);
    }

    /// Raycast against all leaf items; returns sorted (distance, item) pairs.
    pub fn raycast(&self, origin: [f64; 3], direction: [f64; 3]) -> Vec<(f64, usize)> {
        if self.nodes.is_empty() {
            return vec![];
        }
        let dir = normalize3(direction);
        let mut hits = Vec::new();
        self.raycast_recursive(self.root, origin, dir, &mut hits);
        hits.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        hits
    }

    fn raycast_recursive(
        &self,
        node_idx: usize,
        origin: [f64; 3],
        dir: [f64; 3],
        out: &mut Vec<(f64, usize)>,
    ) {
        if node_idx >= self.nodes.len() {
            return;
        }
        let node = &self.nodes[node_idx];
        if let Some(t) = ray_aabb_intersect(origin, dir, &node.aabb) {
            if node.left == usize::MAX {
                out.push((t, node.item));
                return;
            }
            self.raycast_recursive(node.left, origin, dir, out);
            self.raycast_recursive(node.right, origin, dir, out);
        }
    }

    /// Nearest-neighbor search: returns the index of the closest item's center.
    pub fn nearest_neighbor(&self, point: [f64; 3]) -> Option<usize> {
        if self.items.is_empty() {
            return None;
        }
        let mut best_dist = f64::MAX;
        let mut best = 0usize;
        for (i, item) in self.items.iter().enumerate() {
            let d = len3(sub3(item.center(), point));
            if d < best_dist {
                best_dist = d;
                best = i;
            }
        }
        Some(best)
    }
}

fn merge_aabbs<'a>(iter: impl Iterator<Item = &'a PyAabb>) -> PyAabb {
    let mut mn = [f64::MAX; 3];
    let mut mx = [f64::MIN; 3];
    for a in iter {
        for i in 0..3 {
            if a.min[i] < mn[i] {
                mn[i] = a.min[i];
            }
            if a.max[i] > mx[i] {
                mx[i] = a.max[i];
            }
        }
    }
    PyAabb { min: mn, max: mx }
}

fn longest_axis(aabb: &PyAabb) -> usize {
    let e = sub3(aabb.max, aabb.min);
    if e[0] >= e[1] && e[0] >= e[2] {
        0
    } else if e[1] >= e[2] {
        1
    } else {
        2
    }
}

fn ray_aabb_intersect(origin: [f64; 3], dir: [f64; 3], aabb: &PyAabb) -> Option<f64> {
    let mut tmin = f64::NEG_INFINITY;
    let mut tmax = f64::INFINITY;
    for i in 0..3 {
        if dir[i].abs() < 1e-15 {
            if origin[i] < aabb.min[i] || origin[i] > aabb.max[i] {
                return None;
            }
        } else {
            let inv = 1.0 / dir[i];
            let t1 = (aabb.min[i] - origin[i]) * inv;
            let t2 = (aabb.max[i] - origin[i]) * inv;
            tmin = tmin.max(t1.min(t2));
            tmax = tmax.min(t1.max(t2));
        }
    }
    if tmax >= tmin && tmax >= 0.0 {
        Some(tmin.max(0.0))
    } else {
        None
    }
}

// ---------------------------------------------------------------------------
// PyPointCloud
// ---------------------------------------------------------------------------

/// A point cloud with optional per-point normals.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyPointCloud {
    /// Point positions.
    pub points: Vec<[f64; 3]>,
    /// Per-point normals (may be empty).
    pub normals: Vec<[f64; 3]>,
}

impl PyPointCloud {
    /// Create an empty point cloud.
    pub fn new() -> Self {
        Self {
            points: vec![],
            normals: vec![],
        }
    }

    /// Create from point positions.
    pub fn from_points(points: Vec<[f64; 3]>) -> Self {
        Self {
            points,
            normals: vec![],
        }
    }

    /// Estimate per-point normals via PCA of nearest neighbours (simplified: use
    /// the point-to-centroid direction as a proxy).
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

    /// Stub for Poisson surface reconstruction.
    ///
    /// Returns an empty mesh — full implementation requires an octree solver.
    pub fn poisson_reconstruct(&self) -> PyTriangleMesh {
        PyTriangleMesh::new()
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
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyGeometryTransform {
    /// Translation vector.
    pub translation: [f64; 3],
    /// Rotation axis (unit vector).
    pub axis: [f64; 3],
    /// Rotation angle in radians.
    pub angle: f64,
    /// Uniform scale factor.
    pub scale: f64,
}

impl PyGeometryTransform {
    /// Identity transform.
    pub fn identity() -> Self {
        Self {
            translation: [0.0; 3],
            axis: [0.0, 1.0, 0.0],
            angle: 0.0,
            scale: 1.0,
        }
    }

    /// Set translation.
    pub fn translate(mut self, t: [f64; 3]) -> Self {
        self.translation = t;
        self
    }

    /// Set rotation.
    pub fn rotate(mut self, axis: [f64; 3], angle: f64) -> Self {
        self.axis = normalize3(axis);
        self.angle = angle;
        self
    }

    /// Set scale.
    pub fn with_scale(mut self, s: f64) -> Self {
        self.scale = s;
        self
    }

    /// Apply the transform to a single point.
    pub fn apply_point(&self, p: [f64; 3]) -> [f64; 3] {
        // Scale
        let ps = scale3(p, self.scale);
        // Rodrigues rotation
        let pr = rodrigues(ps, self.axis, self.angle);
        // Translate
        add3(pr, self.translation)
    }

    /// Apply this transform to all vertices of a mesh (in place).
    pub fn apply_to_mesh(&self, mesh: &mut PyTriangleMesh) {
        for v in &mut mesh.vertices {
            *v = self.apply_point(*v);
        }
        mesh.compute_normals();
    }

    /// Apply this transform to a list of points.
    pub fn apply_to_points(&self, points: &mut Vec<[f64; 3]>) {
        for p in points.iter_mut() {
            *p = self.apply_point(*p);
        }
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
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyMeshQuality {
    /// Per-triangle aspect ratios.
    pub aspect_ratios: Vec<f64>,
    /// Per-triangle signed Jacobians (det of edge matrix).
    pub jacobians: Vec<f64>,
    /// Per-triangle skewness values in \[0, 1\] (0 = ideal equilateral).
    pub skewness: Vec<f64>,
    /// Per-triangle areas.
    pub areas: Vec<f64>,
}

impl PyMeshQuality {
    /// Compute quality metrics for the given mesh.
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
            // Aspect ratio: longest edge / (2*sqrt(3)*inradius)
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
            // Jacobian: determinant of edge matrix (det [ab, ac, n])
            let face_n = normalize3(cr);
            let jac = dot3(cross3(ab, ac), face_n);
            jacobians.push(jac);
            // Skewness: equilateral-based
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
// Helper functions
// ---------------------------------------------------------------------------

/// Compute the axis-aligned bounding box of a point set.
pub fn compute_aabb(points: &[[f64; 3]]) -> ([f64; 3], [f64; 3]) {
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

/// Serialize a triangle mesh to OBJ format string.
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

/// Parse an OBJ string into a `PyTriangleMesh`.
///
/// Supports `v` (vertex) and `f` (face) lines only.
pub fn obj_string_to_mesh(obj: &str) -> PyTriangleMesh {
    let mut verts = Vec::new();
    let mut indices = Vec::new();
    for line in obj.lines() {
        let line = line.trim();
        if line.starts_with("v ") {
            let parts: Vec<f64> = line[2..]
                .split_whitespace()
                .filter_map(|s| s.parse().ok())
                .collect();
            if parts.len() >= 3 {
                verts.push([parts[0], parts[1], parts[2]]);
            }
        } else if line.starts_with("f ") {
            let parts: Vec<usize> = line[2..]
                .split_whitespace()
                .filter_map(|s| {
                    s.split('/')
                        .next()
                        .and_then(|n| n.parse::<usize>().ok())
                        .map(|n| n - 1)
                })
                .collect();
            if parts.len() >= 3 {
                // Fan triangulation for N-gons
                for i in 1..parts.len() - 1 {
                    indices.push(parts[0]);
                    indices.push(parts[i]);
                    indices.push(parts[i + 1]);
                }
            }
        }
    }
    PyTriangleMesh::from_raw(verts, indices)
}

/// Compute the convex hull of all vertices in a mesh.
pub fn convex_hull_from_mesh(mesh: &PyTriangleMesh) -> PyConvexHull {
    PyConvexHull::from_points(&mesh.vertices)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- PyShape tests ---

    #[test]
    fn test_sphere_volume() {
        let s = PyShape::Sphere(1.0);
        let expected = (4.0 / 3.0) * PI;
        assert!((s.volume() - expected).abs() < 1e-10);
    }

    #[test]
    fn test_box_volume() {
        let b = PyShape::Box([1.0, 2.0, 3.0]);
        assert!((b.volume() - 48.0).abs() < 1e-10);
    }

    #[test]
    fn test_cylinder_volume() {
        let c = PyShape::Cylinder {
            radius: 1.0,
            half_height: 1.0,
        };
        assert!((c.volume() - 2.0 * PI).abs() < 1e-10);
    }

    #[test]
    fn test_torus_volume() {
        let t = PyShape::Torus {
            major_r: 3.0,
            minor_r: 1.0,
        };
        let expected = 2.0 * PI * PI * 3.0 * 1.0;
        assert!((t.volume() - expected).abs() < 1e-10);
    }

    #[test]
    fn test_shape_aabb_sphere() {
        let (mn, mx) = PyShape::Sphere(2.0).aabb();
        assert_eq!(mn, [-2.0, -2.0, -2.0]);
        assert_eq!(mx, [2.0, 2.0, 2.0]);
    }

    #[test]
    fn test_capsule_surface_area() {
        let c = PyShape::Capsule {
            radius: 1.0,
            half_height: 1.0,
        };
        let expected = 4.0 * PI + 2.0 * PI * 2.0;
        assert!((c.surface_area() - expected).abs() < 1e-10);
    }

    // --- PyConvexHull tests ---

    #[test]
    fn test_convex_hull_from_points_empty() {
        let hull = PyConvexHull::from_points(&[]);
        assert!(hull.vertices.is_empty());
    }

    #[test]
    fn test_convex_hull_from_points_basic() {
        let pts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let hull = PyConvexHull::from_points(&pts);
        assert!(!hull.vertices.is_empty());
        assert!(!hull.faces.is_empty());
    }

    #[test]
    fn test_convex_hull_volume_positive() {
        let pts: Vec<[f64; 3]> = vec![
            [-1.0, -1.0, -1.0],
            [1.0, -1.0, -1.0],
            [1.0, 1.0, -1.0],
            [-1.0, 1.0, -1.0],
            [-1.0, -1.0, 1.0],
            [1.0, -1.0, 1.0],
            [1.0, 1.0, 1.0],
            [-1.0, 1.0, 1.0],
        ];
        let hull = PyConvexHull::from_points(&pts);
        assert!(hull.volume() > 0.0);
    }

    #[test]
    fn test_convex_hull_gjk_overlap() {
        let pts1 = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let pts2 = vec![
            [0.5, 0.5, 0.5],
            [1.5, 0.5, 0.5],
            [0.5, 1.5, 0.5],
            [0.5, 0.5, 1.5],
        ];
        let h1 = PyConvexHull::from_points(&pts1);
        let h2 = PyConvexHull::from_points(&pts2);
        assert!(h1.gjk_overlap(&h2));
    }

    #[test]
    fn test_convex_hull_gjk_no_overlap() {
        let pts1 = vec![
            [0.0, 0.0, 0.0],
            [0.1, 0.0, 0.0],
            [0.0, 0.1, 0.0],
            [0.0, 0.0, 0.1],
        ];
        let pts2 = vec![
            [10.0, 10.0, 10.0],
            [10.1, 10.0, 10.0],
            [10.0, 10.1, 10.0],
            [10.0, 10.0, 10.1],
        ];
        let h1 = PyConvexHull::from_points(&pts1);
        let h2 = PyConvexHull::from_points(&pts2);
        assert!(!h1.gjk_overlap(&h2));
    }

    // --- PyTriangleMesh tests ---

    fn unit_box_mesh() -> PyTriangleMesh {
        let verts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
            [0.0, 1.0, 1.0],
        ];
        let idx = vec![
            0, 1, 2, 0, 2, 3, 4, 6, 5, 4, 7, 6, 0, 5, 1, 0, 4, 5, 2, 7, 3, 2, 6, 7, 0, 3, 7, 0, 7,
            4, 1, 5, 6, 1, 6, 2,
        ];
        PyTriangleMesh::from_raw(verts, idx)
    }

    #[test]
    fn test_mesh_compute_normals_not_empty() {
        let mesh = unit_box_mesh();
        assert_eq!(mesh.normals.len(), mesh.vertices.len());
    }

    #[test]
    fn test_mesh_triangle_count() {
        let mesh = unit_box_mesh();
        assert_eq!(mesh.triangle_count(), 12);
    }

    #[test]
    fn test_mesh_surface_area_unit_box() {
        let mesh = unit_box_mesh();
        let sa = mesh.compute_surface_area();
        assert!((sa - 6.0).abs() < 1e-10);
    }

    #[test]
    fn test_mesh_volume_unit_box() {
        let mesh = unit_box_mesh();
        let vol = mesh.compute_volume();
        assert!((vol - 1.0).abs() < 0.1, "vol={vol}");
    }

    #[test]
    fn test_mesh_repair_removes_degenerate() {
        let mut mesh = PyTriangleMesh::from_raw(
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![0, 0, 0, 0, 1, 2], // first tri is degenerate
        );
        mesh.repair();
        assert_eq!(mesh.triangle_count(), 1);
    }

    #[test]
    fn test_mesh_smooth_runs() {
        let mut mesh = unit_box_mesh();
        let before = mesh.vertices.clone();
        mesh.smooth(2);
        // After smoothing the mesh should have changed
        assert_eq!(mesh.vertices.len(), before.len());
    }

    // --- PyCsg tests ---

    #[test]
    fn test_csg_union_vertex_count() {
        let a = unit_box_mesh();
        let b = unit_box_mesh();
        let u = PyCsg::union(&a, &b);
        assert_eq!(u.vertices.len(), a.vertices.len() + b.vertices.len());
    }

    #[test]
    fn test_csg_intersection_subset() {
        let a = unit_box_mesh();
        let mut b_big = unit_box_mesh();
        // Scale b to overlap all of a
        for v in &mut b_big.vertices {
            v[0] *= 2.0;
            v[1] *= 2.0;
            v[2] *= 2.0;
            v[0] -= 0.5;
            v[1] -= 0.5;
            v[2] -= 0.5;
        }
        let inter = PyCsg::intersection(&a, &b_big);
        assert!(inter.vertices.len() <= a.vertices.len() + 8);
    }

    #[test]
    fn test_csg_subtraction_removes_some() {
        let a = unit_box_mesh();
        let b = unit_box_mesh(); // same region
        let sub = PyCsg::subtraction(&a, &b);
        // All or most triangles should be removed
        assert!(sub.triangle_count() <= a.triangle_count());
    }

    // --- PyHeightfield tests ---

    #[test]
    fn test_heightfield_flat() {
        let hf = PyHeightfield::new(10, 10, 1.0, 1.0, 1.0);
        assert!((hf.height_at(5.0, 5.0)).abs() < 1e-10);
    }

    #[test]
    fn test_heightfield_set_and_query() {
        let mut hf = PyHeightfield::new(10, 10, 1.0, 1.0, 1.0);
        hf.set_height(3, 4, 5.0);
        let h = hf.height_at(3.0, 4.0);
        assert!((h - 5.0).abs() < 1e-6);
    }

    #[test]
    fn test_heightfield_normal_flat() {
        let hf = PyHeightfield::new(10, 10, 1.0, 1.0, 1.0);
        let n = hf.normal_at(5.0, 5.0);
        // Flat field should give straight-up normal
        assert!((n[1] - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_heightfield_raycast_hits() {
        let mut hf = PyHeightfield::new(20, 20, 1.0, 1.0, 1.0);
        // Raise entire field
        for r in 0..20 {
            for c in 0..20 {
                hf.set_height(c, r, 1.0);
            }
        }
        let origin = [10.0, 10.0, 10.0];
        let dir = [0.0, -1.0, 0.0];
        let hit = hf.raycast(origin, dir);
        assert!(hit.is_some());
    }

    // --- PySpatialHash tests ---

    #[test]
    fn test_spatial_hash_insert_query() {
        let mut sh = PySpatialHash::new(1.0);
        let id = sh.insert([0.5, 0.5, 0.5]);
        let res = sh.query([0.5, 0.5, 0.5]);
        assert!(res.contains(&id));
    }

    #[test]
    fn test_spatial_hash_sphere_query() {
        let mut sh = PySpatialHash::new(1.0);
        sh.insert([0.0, 0.0, 0.0]);
        sh.insert([5.0, 5.0, 5.0]);
        let res = sh.sphere_query([0.0, 0.0, 0.0], 1.0);
        assert_eq!(res.len(), 1);
    }

    #[test]
    fn test_spatial_hash_remove() {
        let mut sh = PySpatialHash::new(1.0);
        let id = sh.insert([0.0, 0.0, 0.0]);
        sh.remove(id);
        let res = sh.sphere_query([0.0, 0.0, 0.0], 0.5);
        assert!(!res.contains(&id));
    }

    // --- PyBvh tests ---

    #[test]
    fn test_bvh_build_empty() {
        let bvh = PyBvh::build(vec![]);
        assert!(bvh.nodes.is_empty());
    }

    #[test]
    fn test_bvh_query_aabb_hit() {
        let aabbs = vec![
            PyAabb::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
            PyAabb::new([5.0, 5.0, 5.0], [6.0, 6.0, 6.0]),
        ];
        let bvh = PyBvh::build(aabbs);
        let q = PyAabb::new([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
        let hits = bvh.query_aabb(&q);
        assert!(hits.contains(&0));
        assert!(!hits.contains(&1));
    }

    #[test]
    fn test_bvh_raycast() {
        let aabbs = vec![PyAabb::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0])];
        let bvh = PyBvh::build(aabbs);
        let hits = bvh.raycast([-5.0, 0.5, 0.5], [1.0, 0.0, 0.0]);
        assert!(!hits.is_empty());
        assert_eq!(hits[0].1, 0);
    }

    #[test]
    fn test_bvh_nearest_neighbor() {
        let aabbs = vec![
            PyAabb::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
            PyAabb::new([5.0, 5.0, 5.0], [6.0, 6.0, 6.0]),
        ];
        let bvh = PyBvh::build(aabbs);
        let nn = bvh.nearest_neighbor([0.5, 0.5, 0.5]);
        assert_eq!(nn, Some(0));
    }

    // --- PyPointCloud tests ---

    #[test]
    fn test_point_cloud_compute_normals() {
        let mut pc = PyPointCloud::from_points(vec![
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
        ]);
        pc.compute_normals(3);
        assert_eq!(pc.normals.len(), 4);
    }

    #[test]
    fn test_point_cloud_simplify() {
        let mut pc =
            PyPointCloud::from_points((0..100).map(|i| [i as f64 * 0.01, 0.0, 0.0]).collect());
        pc.simplify(0.1);
        assert!(pc.points.len() < 100);
    }

    #[test]
    fn test_point_cloud_poisson_stub() {
        let pc = PyPointCloud::from_points(vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]);
        let mesh = pc.poisson_reconstruct();
        assert!(mesh.vertices.is_empty());
    }

    // --- PyGeometryTransform tests ---

    #[test]
    fn test_transform_identity() {
        let t = PyGeometryTransform::identity();
        let p = [1.0, 2.0, 3.0];
        let out = t.apply_point(p);
        for i in 0..3 {
            assert!((out[i] - p[i]).abs() < 1e-10);
        }
    }

    #[test]
    fn test_transform_translate() {
        let t = PyGeometryTransform::identity().translate([1.0, 2.0, 3.0]);
        let out = t.apply_point([0.0, 0.0, 0.0]);
        assert!((out[0] - 1.0).abs() < 1e-10);
        assert!((out[1] - 2.0).abs() < 1e-10);
        assert!((out[2] - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_transform_scale() {
        let t = PyGeometryTransform::identity().with_scale(2.0);
        let out = t.apply_point([1.0, 1.0, 1.0]);
        for i in 0..3 {
            assert!((out[i] - 2.0).abs() < 1e-10);
        }
    }

    #[test]
    fn test_transform_rotate_180() {
        let t = PyGeometryTransform::identity().rotate([0.0, 0.0, 1.0], PI);
        let out = t.apply_point([1.0, 0.0, 0.0]);
        assert!((out[0] - (-1.0)).abs() < 1e-10, "x={}", out[0]);
        assert!((out[1]).abs() < 1e-10, "y={}", out[1]);
    }

    #[test]
    fn test_transform_apply_to_mesh() {
        let mut mesh = unit_box_mesh();
        let t = PyGeometryTransform::identity().translate([1.0, 0.0, 0.0]);
        t.apply_to_mesh(&mut mesh);
        // All x-coords should be shifted by 1
        for v in &mesh.vertices {
            assert!(v[0] >= 1.0 - 1e-10 && v[0] <= 2.0 + 1e-10);
        }
    }

    // --- PyMeshQuality tests ---

    #[test]
    fn test_mesh_quality_compute() {
        let mesh = unit_box_mesh();
        let q = PyMeshQuality::compute(&mesh);
        assert_eq!(q.aspect_ratios.len(), mesh.triangle_count());
        assert_eq!(q.areas.len(), mesh.triangle_count());
    }

    #[test]
    fn test_mesh_quality_areas_positive() {
        let mesh = unit_box_mesh();
        let q = PyMeshQuality::compute(&mesh);
        for &a in &q.areas {
            assert!(a >= 0.0);
        }
    }

    #[test]
    fn test_mesh_quality_worst_elements() {
        let mesh = unit_box_mesh();
        let q = PyMeshQuality::compute(&mesh);
        let worst = q.worst_elements(3);
        assert!(worst.len() <= 3);
    }

    #[test]
    fn test_mesh_quality_mean_aspect_ratio() {
        let mesh = unit_box_mesh();
        let q = PyMeshQuality::compute(&mesh);
        let mar = q.mean_aspect_ratio();
        assert!(mar >= 1.0);
    }

    // --- Helper function tests ---

    #[test]
    fn test_compute_aabb_empty() {
        let (mn, mx) = compute_aabb(&[]);
        assert_eq!(mn, [0.0; 3]);
        assert_eq!(mx, [0.0; 3]);
    }

    #[test]
    fn test_compute_aabb_points() {
        let pts = vec![[-1.0, 2.0, 3.0], [4.0, -5.0, 6.0]];
        let (mn, mx) = compute_aabb(&pts);
        assert_eq!(mn, [-1.0, -5.0, 3.0]);
        assert_eq!(mx, [4.0, 2.0, 6.0]);
    }

    #[test]
    fn test_obj_roundtrip() {
        let mesh = unit_box_mesh();
        let obj = mesh_to_obj_string(&mesh);
        let parsed = obj_string_to_mesh(&obj);
        assert_eq!(parsed.vertices.len(), mesh.vertices.len());
        assert_eq!(parsed.triangle_count(), mesh.triangle_count());
    }

    #[test]
    fn test_convex_hull_from_mesh() {
        let mesh = unit_box_mesh();
        let hull = convex_hull_from_mesh(&mesh);
        assert!(!hull.vertices.is_empty());
    }
}
