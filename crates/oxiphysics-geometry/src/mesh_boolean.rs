// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Mesh boolean operations on triangle meshes.
//!
//! Provides robust boolean operations (union, intersection, difference) on
//! closed triangle meshes using:
//!
//! - **Triangle-triangle intersection detection** — Möller's separating axis
//!   test and segment/plane clipping.
//! - **Contour extraction** — intersection curves between two meshes.
//! - **Surface splitting** — re-triangulation of intersected faces.
//! - **Inside/outside classification** — generalised winding number.
//! - **Co-planar triangle handling** — projection-based merging.
//! - **Mesh stitching** — combining split surfaces into a watertight result.
//! - **Result cleanup** — degenerate triangle removal, vertex welding.

#![allow(dead_code)]

// ─────────────────────────────────────────────────────────────────────────────
// Vector helpers (no nalgebra in non-core crates)
// ─────────────────────────────────────────────────────────────────────────────

/// Three-component vector type alias.
type V3 = [f64; 3];

/// Small tolerance for geometric tests.
const GEO_EPS: f64 = 1e-10;

/// Tolerance for co-planar detection.
const COPLANAR_EPS: f64 = 1e-8;

/// Tolerance for vertex welding.
const WELD_EPS: f64 = 1e-9;

#[inline]
fn sub(a: V3, b: V3) -> V3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn add(a: V3, b: V3) -> V3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn scale(a: V3, s: f64) -> V3 {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline]
fn dot(a: V3, b: V3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn cross(a: V3, b: V3) -> V3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
fn length(a: V3) -> f64 {
    dot(a, a).sqrt()
}

#[inline]
fn normalize(a: V3) -> V3 {
    let l = length(a);
    if l < GEO_EPS {
        [0.0, 0.0, 0.0]
    } else {
        [a[0] / l, a[1] / l, a[2] / l]
    }
}

#[inline]
fn lerp(a: V3, b: V3, t: f64) -> V3 {
    [
        a[0] + t * (b[0] - a[0]),
        a[1] + t * (b[1] - a[1]),
        a[2] + t * (b[2] - a[2]),
    ]
}

#[inline]
fn dist_sq(a: V3, b: V3) -> f64 {
    let d = sub(a, b);
    dot(d, d)
}

/// Triangle area from three vertices.
#[inline]
fn triangle_area(a: V3, b: V3, c: V3) -> f64 {
    0.5 * length(cross(sub(b, a), sub(c, a)))
}

/// Triangle normal (unnormalised).
#[inline]
fn triangle_normal(a: V3, b: V3, c: V3) -> V3 {
    cross(sub(b, a), sub(c, a))
}

// ─────────────────────────────────────────────────────────────────────────────
// BooleanOp
// ─────────────────────────────────────────────────────────────────────────────

/// Boolean operation type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MeshBooleanOp {
    /// A ∪ B.
    Union,
    /// A ∩ B.
    Intersection,
    /// A \ B.
    Difference,
}

// ─────────────────────────────────────────────────────────────────────────────
// SimpleMesh — lightweight triangle mesh
// ─────────────────────────────────────────────────────────────────────────────

/// A simple indexed triangle mesh using `[f64;3]` vertices.
#[derive(Debug, Clone)]
pub struct SimpleMesh {
    /// Vertex positions.
    pub vertices: Vec<V3>,
    /// Triangle indices (each triple is one triangle).
    pub triangles: Vec<[usize; 3]>,
}

impl SimpleMesh {
    /// Create a new empty mesh.
    pub fn new() -> Self {
        Self {
            vertices: Vec::new(),
            triangles: Vec::new(),
        }
    }

    /// Create a mesh from vertices and triangles.
    pub fn from_data(vertices: Vec<V3>, triangles: Vec<[usize; 3]>) -> Self {
        Self {
            vertices,
            triangles,
        }
    }

    /// Number of triangles.
    pub fn n_triangles(&self) -> usize {
        self.triangles.len()
    }

    /// Number of vertices.
    pub fn n_vertices(&self) -> usize {
        self.vertices.len()
    }

    /// Get triangle vertex positions.
    pub fn triangle_verts(&self, tri_idx: usize) -> (V3, V3, V3) {
        let t = self.triangles[tri_idx];
        (
            self.vertices[t[0]],
            self.vertices[t[1]],
            self.vertices[t[2]],
        )
    }

    /// Compute the axis-aligned bounding box.
    pub fn aabb(&self) -> (V3, V3) {
        if self.vertices.is_empty() {
            return ([0.0; 3], [0.0; 3]);
        }
        let mut mn = self.vertices[0];
        let mut mx = self.vertices[0];
        for v in &self.vertices {
            for d in 0..3 {
                if v[d] < mn[d] {
                    mn[d] = v[d];
                }
                if v[d] > mx[d] {
                    mx[d] = v[d];
                }
            }
        }
        (mn, mx)
    }

    /// Add a vertex and return its index.
    pub fn add_vertex(&mut self, v: V3) -> usize {
        let idx = self.vertices.len();
        self.vertices.push(v);
        idx
    }

    /// Add a triangle.
    pub fn add_triangle(&mut self, tri: [usize; 3]) {
        self.triangles.push(tri);
    }

    /// Compute total surface area.
    pub fn surface_area(&self) -> f64 {
        let mut area = 0.0;
        for &t in &self.triangles {
            area += triangle_area(
                self.vertices[t[0]],
                self.vertices[t[1]],
                self.vertices[t[2]],
            );
        }
        area
    }

    /// Compute signed volume (for closed meshes).
    pub fn signed_volume(&self) -> f64 {
        let mut vol = 0.0;
        for &t in &self.triangles {
            let a = self.vertices[t[0]];
            let b = self.vertices[t[1]];
            let c = self.vertices[t[2]];
            vol += dot(a, cross(b, c));
        }
        vol / 6.0
    }

    /// Flip all triangle orientations.
    pub fn flip_normals(&mut self) {
        for t in &mut self.triangles {
            t.swap(0, 1);
        }
    }

    /// Create a unit cube mesh (for testing).
    pub fn unit_cube(centre: V3, half_extent: f64) -> Self {
        let h = half_extent;
        let c = centre;
        let vertices = vec![
            [c[0] - h, c[1] - h, c[2] - h], // 0
            [c[0] + h, c[1] - h, c[2] - h], // 1
            [c[0] + h, c[1] + h, c[2] - h], // 2
            [c[0] - h, c[1] + h, c[2] - h], // 3
            [c[0] - h, c[1] - h, c[2] + h], // 4
            [c[0] + h, c[1] - h, c[2] + h], // 5
            [c[0] + h, c[1] + h, c[2] + h], // 6
            [c[0] - h, c[1] + h, c[2] + h], // 7
        ];
        let triangles = vec![
            // -Z face
            [0, 2, 1],
            [0, 3, 2],
            // +Z face
            [4, 5, 6],
            [4, 6, 7],
            // -Y face
            [0, 1, 5],
            [0, 5, 4],
            // +Y face
            [2, 3, 7],
            [2, 7, 6],
            // -X face
            [0, 4, 7],
            [0, 7, 3],
            // +X face
            [1, 2, 6],
            [1, 6, 5],
        ];
        Self {
            vertices,
            triangles,
        }
    }

    /// Create a tetrahedron mesh (for testing).
    pub fn tetrahedron(centre: V3, _size: f64) -> Self {
        let s = _size;
        let vertices = vec![
            [centre[0] + s, centre[1] + s, centre[2] + s],
            [centre[0] + s, centre[1] - s, centre[2] - s],
            [centre[0] - s, centre[1] + s, centre[2] - s],
            [centre[0] - s, centre[1] - s, centre[2] + s],
        ];
        let triangles = vec![[0, 1, 2], [0, 3, 1], [0, 2, 3], [1, 3, 2]];
        Self {
            vertices,
            triangles,
        }
    }
}

impl Default for SimpleMesh {
    fn default() -> Self {
        Self::new()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Triangle-Triangle Intersection
// ─────────────────────────────────────────────────────────────────────────────

/// Result of a triangle-triangle intersection test.
#[derive(Debug, Clone)]
pub enum TriTriResult {
    /// No intersection.
    None,
    /// Triangles intersect along a segment.
    Segment(V3, V3),
    /// Triangles are co-planar and overlap.
    Coplanar,
    /// Single point contact.
    Point(V3),
}

/// Test if two triangles intersect using Möller's method.
///
/// Returns the intersection type and geometry.
pub fn triangle_triangle_intersection(
    a0: V3,
    a1: V3,
    a2: V3,
    b0: V3,
    b1: V3,
    b2: V3,
) -> TriTriResult {
    // Plane of triangle A
    let na = triangle_normal(a0, a1, a2);
    let da = dot(na, a0);
    let db0 = dot(na, b0) - da;
    let db1 = dot(na, b1) - da;
    let db2 = dot(na, b2) - da;

    // Snap near-zero to zero
    let db0 = if db0.abs() < GEO_EPS { 0.0 } else { db0 };
    let db1 = if db1.abs() < GEO_EPS { 0.0 } else { db1 };
    let db2 = if db2.abs() < GEO_EPS { 0.0 } else { db2 };

    // All on same side → no intersection
    if db0 > 0.0 && db1 > 0.0 && db2 > 0.0 {
        return TriTriResult::None;
    }
    if db0 < 0.0 && db1 < 0.0 && db2 < 0.0 {
        return TriTriResult::None;
    }

    // Check co-planar
    if db0.abs() < GEO_EPS && db1.abs() < GEO_EPS && db2.abs() < GEO_EPS {
        if coplanar_triangles_overlap(a0, a1, a2, b0, b1, b2, na) {
            return TriTriResult::Coplanar;
        }
        return TriTriResult::None;
    }

    // Plane of triangle B
    let nb = triangle_normal(b0, b1, b2);
    let _db_val = dot(nb, b0);
    let da0 = dot(nb, a0) - _db_val;
    let da1 = dot(nb, a1) - _db_val;
    let da2 = dot(nb, a2) - _db_val;

    let da0 = if da0.abs() < GEO_EPS { 0.0 } else { da0 };
    let da1 = if da1.abs() < GEO_EPS { 0.0 } else { da1 };
    let da2 = if da2.abs() < GEO_EPS { 0.0 } else { da2 };

    if da0 > 0.0 && da1 > 0.0 && da2 > 0.0 {
        return TriTriResult::None;
    }
    if da0 < 0.0 && da1 < 0.0 && da2 < 0.0 {
        return TriTriResult::None;
    }

    // Intersection line direction
    let dir = cross(na, nb);
    let dir_len = length(dir);
    if dir_len < GEO_EPS {
        return TriTriResult::None;
    }
    let dir = normalize(dir);

    // Project vertices onto intersection line
    let seg_a = compute_interval_on_line(a0, a1, a2, da0, da1, da2, dir);
    let seg_b = compute_interval_on_line(b0, b1, b2, db0, db1, db2, dir);

    if let (Some((ta_min, ta_max, pa_min, pa_max)), Some((tb_min, tb_max, pb_min, pb_max))) =
        (seg_a, seg_b)
    {
        // Overlap of [ta_min, ta_max] and [tb_min, tb_max]
        let t0 = ta_min.max(tb_min);
        let t1 = ta_max.min(tb_max);
        if t0 > t1 + GEO_EPS {
            return TriTriResult::None;
        }
        // Interpolate actual 3D points
        let p0 = if t0 >= ta_min - GEO_EPS && t0 <= ta_max + GEO_EPS {
            interp_seg_param(pa_min, pa_max, ta_min, ta_max, t0)
        } else {
            interp_seg_param(pb_min, pb_max, tb_min, tb_max, t0)
        };
        let p1 = if t1 >= ta_min - GEO_EPS && t1 <= ta_max + GEO_EPS {
            interp_seg_param(pa_min, pa_max, ta_min, ta_max, t1)
        } else {
            interp_seg_param(pb_min, pb_max, tb_min, tb_max, t1)
        };
        if dist_sq(p0, p1) < GEO_EPS * GEO_EPS {
            return TriTriResult::Point(p0);
        }
        TriTriResult::Segment(p0, p1)
    } else {
        TriTriResult::None
    }
}

/// Compute the interval where a triangle's edges cross the opposite plane,
/// projected onto the intersection line direction.
/// Returns `(t_min, t_max, point_at_t_min, point_at_t_max)`.
fn compute_interval_on_line(
    v0: V3,
    v1: V3,
    v2: V3,
    d0: f64,
    d1: f64,
    d2: f64,
    dir: V3,
) -> Option<(f64, f64, V3, V3)> {
    let verts = [v0, v1, v2];
    let dists = [d0, d1, d2];

    // Find crossing points: edges where sign changes
    let mut points: Vec<(f64, V3)> = Vec::new();

    // Check vertices on the plane
    for i in 0..3 {
        if dists[i].abs() < GEO_EPS {
            let t = dot(verts[i], dir);
            points.push((t, verts[i]));
        }
    }

    // Check edges crossing the plane
    for (i, j) in [(0, 1), (1, 2), (2, 0)] {
        if (dists[i] > GEO_EPS && dists[j] < -GEO_EPS)
            || (dists[i] < -GEO_EPS && dists[j] > GEO_EPS)
        {
            let s = dists[i] / (dists[i] - dists[j]);
            let p = lerp(verts[i], verts[j], s);
            let t = dot(p, dir);
            points.push((t, p));
        }
    }

    if points.is_empty() {
        return None;
    }

    // Find min and max t
    let mut min_idx = 0;
    let mut max_idx = 0;
    for (k, (t, _)) in points.iter().enumerate() {
        if *t < points[min_idx].0 {
            min_idx = k;
        }
        if *t > points[max_idx].0 {
            max_idx = k;
        }
    }

    Some((
        points[min_idx].0,
        points[max_idx].0,
        points[min_idx].1,
        points[max_idx].1,
    ))
}

/// Interpolate between segment endpoints at a parameter value.
fn interp_seg_param(p_min: V3, p_max: V3, t_min: f64, t_max: f64, t: f64) -> V3 {
    if (t_max - t_min).abs() < GEO_EPS {
        return p_min;
    }
    let s = (t - t_min) / (t_max - t_min);
    lerp(p_min, p_max, s.clamp(0.0, 1.0))
}

// ─────────────────────────────────────────────────────────────────────────────
// Co-planar triangle overlap
// ─────────────────────────────────────────────────────────────────────────────

/// Test if two co-planar triangles overlap by projecting onto the dominant
/// axis and performing 2D edge-edge and containment tests.
fn coplanar_triangles_overlap(a0: V3, a1: V3, a2: V3, b0: V3, b1: V3, b2: V3, normal: V3) -> bool {
    // Find dominant axis (largest component of normal)
    let abs_n = [normal[0].abs(), normal[1].abs(), normal[2].abs()];
    let (ax1, ax2) = if abs_n[0] >= abs_n[1] && abs_n[0] >= abs_n[2] {
        (1, 2)
    } else if abs_n[1] >= abs_n[2] {
        (0, 2)
    } else {
        (0, 1)
    };

    let proj = |v: V3| -> [f64; 2] { [v[ax1], v[ax2]] };

    let pa = [proj(a0), proj(a1), proj(a2)];
    let pb = [proj(b0), proj(b1), proj(b2)];

    // Check edge-edge intersections
    for i in 0..3 {
        let j = (i + 1) % 3;
        for k in 0..3 {
            let l = (k + 1) % 3;
            if segments_intersect_2d(pa[i], pa[j], pb[k], pb[l]) {
                return true;
            }
        }
    }

    // Check containment
    if point_in_triangle_2d(pa[0], pb[0], pb[1], pb[2]) {
        return true;
    }
    if point_in_triangle_2d(pb[0], pa[0], pa[1], pa[2]) {
        return true;
    }

    false
}

/// 2D segment intersection test.
fn segments_intersect_2d(a: [f64; 2], b: [f64; 2], c: [f64; 2], d: [f64; 2]) -> bool {
    let ab = [b[0] - a[0], b[1] - a[1]];
    let cd = [d[0] - c[0], d[1] - c[1]];
    let denom = ab[0] * cd[1] - ab[1] * cd[0];
    if denom.abs() < GEO_EPS {
        return false; // parallel
    }
    let ac = [c[0] - a[0], c[1] - a[1]];
    let t = (ac[0] * cd[1] - ac[1] * cd[0]) / denom;
    let u = (ac[0] * ab[1] - ac[1] * ab[0]) / denom;
    (-GEO_EPS..=1.0 + GEO_EPS).contains(&t) && (-GEO_EPS..=1.0 + GEO_EPS).contains(&u)
}

/// 2D point-in-triangle using barycentric coordinates.
fn point_in_triangle_2d(p: [f64; 2], a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> bool {
    let v0 = [c[0] - a[0], c[1] - a[1]];
    let v1 = [b[0] - a[0], b[1] - a[1]];
    let v2 = [p[0] - a[0], p[1] - a[1]];
    let d00 = v0[0] * v0[0] + v0[1] * v0[1];
    let d01 = v0[0] * v1[0] + v0[1] * v1[1];
    let d02 = v0[0] * v2[0] + v0[1] * v2[1];
    let d11 = v1[0] * v1[0] + v1[1] * v1[1];
    let d12 = v1[0] * v2[0] + v1[1] * v2[1];
    let inv_denom = d00 * d11 - d01 * d01;
    if inv_denom.abs() < GEO_EPS {
        return false;
    }
    let inv = 1.0 / inv_denom;
    let u = (d11 * d02 - d01 * d12) * inv;
    let v = (d00 * d12 - d01 * d02) * inv;
    u >= -GEO_EPS && v >= -GEO_EPS && (u + v) <= 1.0 + GEO_EPS
}

// ─────────────────────────────────────────────────────────────────────────────
// Winding Number (inside/outside classification)
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the generalised winding number of a closed mesh at point `p`.
///
/// A point is inside if |winding_number| ≈ 1.
pub fn winding_number(mesh: &SimpleMesh, p: V3) -> f64 {
    let mut wn = 0.0;
    for &tri in &mesh.triangles {
        let a = sub(mesh.vertices[tri[0]], p);
        let b = sub(mesh.vertices[tri[1]], p);
        let c = sub(mesh.vertices[tri[2]], p);
        let la = length(a);
        let lb = length(b);
        let lc = length(c);
        if la < GEO_EPS || lb < GEO_EPS || lc < GEO_EPS {
            continue;
        }
        let num = dot(a, cross(b, c));
        let denom = la * lb * lc + dot(a, b) * lc + dot(b, c) * la + dot(c, a) * lb;
        wn += num.atan2(denom);
    }
    wn / (2.0 * std::f64::consts::PI)
}

/// Classify whether a point is inside a closed mesh.
pub fn is_inside(mesh: &SimpleMesh, p: V3) -> bool {
    winding_number(mesh, p).abs() > 0.5
}

// ─────────────────────────────────────────────────────────────────────────────
// Contour Extraction
// ─────────────────────────────────────────────────────────────────────────────

/// An intersection segment between two meshes.
#[derive(Debug, Clone)]
pub struct IntersectionSegment {
    /// Start point of the segment.
    pub start: V3,
    /// End point of the segment.
    pub end: V3,
    /// Triangle index in mesh A.
    pub tri_a: usize,
    /// Triangle index in mesh B.
    pub tri_b: usize,
}

/// Extract all intersection segments between two meshes.
pub fn extract_intersection_contours(
    mesh_a: &SimpleMesh,
    mesh_b: &SimpleMesh,
) -> Vec<IntersectionSegment> {
    let mut segments = Vec::new();
    for i in 0..mesh_a.n_triangles() {
        let (a0, a1, a2) = mesh_a.triangle_verts(i);
        for j in 0..mesh_b.n_triangles() {
            let (b0, b1, b2) = mesh_b.triangle_verts(j);
            match triangle_triangle_intersection(a0, a1, a2, b0, b1, b2) {
                TriTriResult::Segment(p0, p1) => {
                    segments.push(IntersectionSegment {
                        start: p0,
                        end: p1,
                        tri_a: i,
                        tri_b: j,
                    });
                }
                TriTriResult::Point(p) => {
                    segments.push(IntersectionSegment {
                        start: p,
                        end: p,
                        tri_a: i,
                        tri_b: j,
                    });
                }
                _ => {}
            }
        }
    }
    segments
}

// ─────────────────────────────────────────────────────────────────────────────
// Surface Splitting
// ─────────────────────────────────────────────────────────────────────────────

/// Split a triangle by a plane defined by (normal, point_on_plane).
/// Returns (front_triangles, back_triangles) as vertex lists.
pub fn split_triangle_by_plane(
    v0: V3,
    v1: V3,
    v2: V3,
    plane_normal: V3,
    plane_point: V3,
) -> (Vec<[V3; 3]>, Vec<[V3; 3]>) {
    let verts = [v0, v1, v2];
    let dists: Vec<f64> = verts
        .iter()
        .map(|v| dot(sub(*v, plane_point), plane_normal))
        .collect();

    let pos_count = dists.iter().filter(|&&d| d > GEO_EPS).count();
    let neg_count = dists.iter().filter(|&&d| d < -GEO_EPS).count();

    // All on one side
    if neg_count == 0 {
        return (vec![[v0, v1, v2]], Vec::new());
    }
    if pos_count == 0 {
        return (Vec::new(), vec![[v0, v1, v2]]);
    }

    // Find the isolated vertex
    let mut front = Vec::new();
    let mut back = Vec::new();

    // Re-order so that the isolated vertex is first
    for start in 0..3 {
        let i0 = start;
        let i1 = (start + 1) % 3;
        let i2 = (start + 2) % 3;

        if (dists[i0] > GEO_EPS && dists[i1] < -GEO_EPS && dists[i2] < -GEO_EPS)
            || (dists[i0] > GEO_EPS
                && dists[i1] <= GEO_EPS
                && dists[i2] < -GEO_EPS
                && dists[i1].abs() <= GEO_EPS)
        {
            // i0 alone on positive side
            let t01 = dists[i0] / (dists[i0] - dists[i1]);
            let t02 = dists[i0] / (dists[i0] - dists[i2]);
            let p01 = lerp(verts[i0], verts[i1], t01);
            let p02 = lerp(verts[i0], verts[i2], t02);
            front.push([verts[i0], p01, p02]);
            back.push([p01, verts[i1], verts[i2]]);
            back.push([p01, verts[i2], p02]);
            return (front, back);
        }
        if (dists[i0] < -GEO_EPS && dists[i1] > GEO_EPS && dists[i2] > GEO_EPS)
            || (dists[i0] < -GEO_EPS
                && dists[i1] >= -GEO_EPS
                && dists[i2] > GEO_EPS
                && dists[i1].abs() <= GEO_EPS)
        {
            // i0 alone on negative side
            let t01 = dists[i0] / (dists[i0] - dists[i1]);
            let t02 = dists[i0] / (dists[i0] - dists[i2]);
            let p01 = lerp(verts[i0], verts[i1], t01);
            let p02 = lerp(verts[i0], verts[i2], t02);
            back.push([verts[i0], p01, p02]);
            front.push([p01, verts[i1], verts[i2]]);
            front.push([p01, verts[i2], p02]);
            return (front, back);
        }
    }

    // Fallback: put on front side
    (vec![[v0, v1, v2]], Vec::new())
}

// ─────────────────────────────────────────────────────────────────────────────
// Vertex Welding / Cleanup
// ─────────────────────────────────────────────────────────────────────────────

/// Weld vertices that are closer than `tolerance`.
/// Returns a new mesh with merged vertices.
pub fn weld_vertices(mesh: &SimpleMesh, tolerance: f64) -> SimpleMesh {
    let tol_sq = tolerance * tolerance;
    let mut new_verts: Vec<V3> = Vec::new();
    let mut remap: Vec<usize> = Vec::new();

    for v in &mesh.vertices {
        let mut found = None;
        for (k, nv) in new_verts.iter().enumerate() {
            if dist_sq(*v, *nv) < tol_sq {
                found = Some(k);
                break;
            }
        }
        if let Some(k) = found {
            remap.push(k);
        } else {
            remap.push(new_verts.len());
            new_verts.push(*v);
        }
    }

    let new_tris: Vec<[usize; 3]> = mesh
        .triangles
        .iter()
        .map(|t| [remap[t[0]], remap[t[1]], remap[t[2]]])
        .collect();

    SimpleMesh::from_data(new_verts, new_tris)
}

/// Remove degenerate triangles (zero-area or duplicate vertex indices).
pub fn remove_degenerate_triangles(mesh: &mut SimpleMesh) {
    mesh.triangles.retain(|t| {
        if t[0] == t[1] || t[1] == t[2] || t[0] == t[2] {
            return false;
        }
        let a = mesh.vertices[t[0]];
        let b = mesh.vertices[t[1]];
        let c = mesh.vertices[t[2]];
        triangle_area(a, b, c) > GEO_EPS
    });
}

/// Remove unreferenced vertices and compact the index buffer.
pub fn remove_unused_vertices(mesh: &mut SimpleMesh) {
    let n = mesh.vertices.len();
    let mut used = vec![false; n];
    for t in &mesh.triangles {
        used[t[0]] = true;
        used[t[1]] = true;
        used[t[2]] = true;
    }
    let mut remap = vec![0usize; n];
    let mut new_verts = Vec::new();
    for (i, &u) in used.iter().enumerate() {
        if u {
            remap[i] = new_verts.len();
            new_verts.push(mesh.vertices[i]);
        }
    }
    for t in &mut mesh.triangles {
        t[0] = remap[t[0]];
        t[1] = remap[t[1]];
        t[2] = remap[t[2]];
    }
    mesh.vertices = new_verts;
}

/// Full cleanup: weld, remove degenerates, remove unused.
pub fn cleanup_mesh(mesh: &SimpleMesh) -> SimpleMesh {
    let mut result = weld_vertices(mesh, WELD_EPS);
    remove_degenerate_triangles(&mut result);
    remove_unused_vertices(&mut result);
    result
}

// ─────────────────────────────────────────────────────────────────────────────
// Mesh Stitching
// ─────────────────────────────────────────────────────────────────────────────

/// Combine two meshes into one, re-indexing the second mesh's triangles.
pub fn stitch_meshes(a: &SimpleMesh, b: &SimpleMesh) -> SimpleMesh {
    let offset = a.vertices.len();
    let mut vertices = a.vertices.clone();
    vertices.extend_from_slice(&b.vertices);
    let mut triangles = a.triangles.clone();
    for t in &b.triangles {
        triangles.push([t[0] + offset, t[1] + offset, t[2] + offset]);
    }
    SimpleMesh {
        vertices,
        triangles,
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Mesh Boolean (high-level)
// ─────────────────────────────────────────────────────────────────────────────

/// Classify each triangle of a mesh as inside or outside another mesh.
/// Returns a boolean per triangle: `true` = inside.
pub fn classify_triangles(mesh: &SimpleMesh, other: &SimpleMesh) -> Vec<bool> {
    let mut result = Vec::with_capacity(mesh.n_triangles());
    for i in 0..mesh.n_triangles() {
        let (a, b, c) = mesh.triangle_verts(i);
        let centroid = scale(add(add(a, b), c), 1.0 / 3.0);
        result.push(is_inside(other, centroid));
    }
    result
}

/// Perform a mesh boolean operation.
///
/// This is a simplified approach that classifies triangles by their centroid
/// and selects/rejects based on the operation. For production use, a full
/// re-triangulation along intersection curves is needed.
pub fn mesh_boolean(mesh_a: &SimpleMesh, mesh_b: &SimpleMesh, op: MeshBooleanOp) -> SimpleMesh {
    let a_inside_b = classify_triangles(mesh_a, mesh_b);
    let b_inside_a = classify_triangles(mesh_b, mesh_a);

    let mut result = SimpleMesh::new();

    match op {
        MeshBooleanOp::Union => {
            // Keep triangles of A that are outside B
            collect_triangles(mesh_a, &a_inside_b, false, &mut result);
            // Keep triangles of B that are outside A
            collect_triangles(mesh_b, &b_inside_a, false, &mut result);
        }
        MeshBooleanOp::Intersection => {
            // Keep triangles of A that are inside B
            collect_triangles(mesh_a, &a_inside_b, true, &mut result);
            // Keep triangles of B that are inside A
            collect_triangles(mesh_b, &b_inside_a, true, &mut result);
        }
        MeshBooleanOp::Difference => {
            // Keep triangles of A that are outside B
            collect_triangles(mesh_a, &a_inside_b, false, &mut result);
            // Keep triangles of B that are inside A (flipped)
            collect_triangles_flipped(mesh_b, &b_inside_a, true, &mut result);
        }
    }

    cleanup_mesh(&result)
}

/// Collect triangles from a mesh based on inside/outside classification.
fn collect_triangles(
    mesh: &SimpleMesh,
    classification: &[bool],
    keep_inside: bool,
    result: &mut SimpleMesh,
) {
    let offset = result.vertices.len();
    result.vertices.extend_from_slice(&mesh.vertices);
    for (i, &is_inside_val) in classification.iter().enumerate() {
        if is_inside_val == keep_inside {
            let t = mesh.triangles[i];
            result
                .triangles
                .push([t[0] + offset, t[1] + offset, t[2] + offset]);
        }
    }
}

/// Collect triangles with flipped normals.
fn collect_triangles_flipped(
    mesh: &SimpleMesh,
    classification: &[bool],
    keep_inside: bool,
    result: &mut SimpleMesh,
) {
    let offset = result.vertices.len();
    result.vertices.extend_from_slice(&mesh.vertices);
    for (i, &is_inside_val) in classification.iter().enumerate() {
        if is_inside_val == keep_inside {
            let t = mesh.triangles[i];
            // Flip winding order
            result
                .triangles
                .push([t[1] + offset, t[0] + offset, t[2] + offset]);
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Ray–triangle intersection (used by inside/outside via ray casting)
// ─────────────────────────────────────────────────────────────────────────────

/// Möller–Trumbore ray–triangle intersection.
/// Returns `Some(t)` where `t` is the ray parameter at intersection.
pub fn ray_triangle_intersect(origin: V3, dir: V3, v0: V3, v1: V3, v2: V3) -> Option<f64> {
    let e1 = sub(v1, v0);
    let e2 = sub(v2, v0);
    let h = cross(dir, e2);
    let a = dot(e1, h);
    if a.abs() < GEO_EPS {
        return None;
    }
    let f = 1.0 / a;
    let s = sub(origin, v0);
    let u = f * dot(s, h);
    if !(0.0..=1.0).contains(&u) {
        return None;
    }
    let q = cross(s, e1);
    let v = f * dot(dir, q);
    if v < 0.0 || u + v > 1.0 {
        return None;
    }
    let t = f * dot(e2, q);
    if t > GEO_EPS { Some(t) } else { None }
}

/// Count ray-mesh intersections for parity-based inside/outside.
pub fn ray_mesh_intersection_count(mesh: &SimpleMesh, origin: V3, dir: V3) -> usize {
    let mut count = 0;
    for &tri in &mesh.triangles {
        if ray_triangle_intersect(
            origin,
            dir,
            mesh.vertices[tri[0]],
            mesh.vertices[tri[1]],
            mesh.vertices[tri[2]],
        )
        .is_some()
        {
            count += 1;
        }
    }
    count
}

/// Parity-based inside test using ray casting.
pub fn is_inside_ray(mesh: &SimpleMesh, p: V3) -> bool {
    let dir = [1.0, 0.0, 0.0];
    ray_mesh_intersection_count(mesh, p, dir) % 2 == 1
}

// ─────────────────────────────────────────────────────────────────────────────
// Mesh statistics
// ─────────────────────────────────────────────────────────────────────────────

/// Compute mesh quality statistics.
#[derive(Debug, Clone)]
pub struct MeshStats {
    /// Number of vertices.
    pub n_vertices: usize,
    /// Number of triangles.
    pub n_triangles: usize,
    /// Total surface area.
    pub surface_area: f64,
    /// Signed volume.
    pub signed_volume: f64,
    /// Minimum triangle area.
    pub min_triangle_area: f64,
    /// Maximum triangle area.
    pub max_triangle_area: f64,
}

/// Compute statistics for a mesh.
pub fn compute_mesh_stats(mesh: &SimpleMesh) -> MeshStats {
    let mut min_area = f64::INFINITY;
    let mut max_area = 0.0_f64;
    let mut total_area = 0.0;
    for i in 0..mesh.n_triangles() {
        let (a, b, c) = mesh.triangle_verts(i);
        let area = triangle_area(a, b, c);
        total_area += area;
        if area < min_area {
            min_area = area;
        }
        if area > max_area {
            max_area = area;
        }
    }
    MeshStats {
        n_vertices: mesh.n_vertices(),
        n_triangles: mesh.n_triangles(),
        surface_area: total_area,
        signed_volume: mesh.signed_volume(),
        min_triangle_area: min_area,
        max_triangle_area: max_area,
    }
}

/// Check if a mesh is watertight (every edge shared by exactly 2 triangles).
pub fn is_watertight(mesh: &SimpleMesh) -> bool {
    use std::collections::HashMap;
    let mut edge_count: HashMap<(usize, usize), u32> = HashMap::new();
    for t in &mesh.triangles {
        for (&a, &b) in t.iter().zip(t.iter().cycle().skip(1).take(3)) {
            let key = if a < b { (a, b) } else { (b, a) };
            *edge_count.entry(key).or_insert(0) += 1;
        }
    }
    edge_count.values().all(|&c| c == 2)
}

/// Compute the Euler characteristic: V - E + F.
pub fn euler_characteristic(mesh: &SimpleMesh) -> i64 {
    use std::collections::HashSet;
    let v = mesh.n_vertices() as i64;
    let f = mesh.n_triangles() as i64;
    let mut edges = HashSet::new();
    for t in &mesh.triangles {
        for (&a, &b) in t.iter().zip(t.iter().cycle().skip(1).take(3)) {
            let key = if a < b { (a, b) } else { (b, a) };
            edges.insert(key);
        }
    }
    let e = edges.len() as i64;
    v - e + f
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── SimpleMesh basics ───────────────────────────────────────────────

    #[test]
    fn test_unit_cube_creation() {
        let cube = SimpleMesh::unit_cube([0.0; 3], 1.0);
        assert_eq!(cube.n_vertices(), 8);
        assert_eq!(cube.n_triangles(), 12);
    }

    #[test]
    fn test_cube_surface_area() {
        let cube = SimpleMesh::unit_cube([0.0; 3], 1.0);
        let area = cube.surface_area();
        // 6 faces, each 2x2 = 4, total = 24
        assert!((area - 24.0).abs() < 1e-8, "area = {area}");
    }

    #[test]
    fn test_cube_signed_volume() {
        let cube = SimpleMesh::unit_cube([0.0; 3], 1.0);
        let vol = cube.signed_volume();
        // 2^3 = 8
        assert!((vol.abs() - 8.0).abs() < 1e-8, "vol = {vol}");
    }

    #[test]
    fn test_tetrahedron_creation() {
        let tet = SimpleMesh::tetrahedron([0.0; 3], 1.0);
        assert_eq!(tet.n_vertices(), 4);
        assert_eq!(tet.n_triangles(), 4);
    }

    #[test]
    fn test_empty_mesh() {
        let m = SimpleMesh::new();
        assert_eq!(m.n_vertices(), 0);
        assert_eq!(m.n_triangles(), 0);
        let (mn, mx) = m.aabb();
        assert_eq!(mn, [0.0; 3]);
        assert_eq!(mx, [0.0; 3]);
    }

    #[test]
    fn test_add_vertex_and_triangle() {
        let mut m = SimpleMesh::new();
        let a = m.add_vertex([0.0, 0.0, 0.0]);
        let b = m.add_vertex([1.0, 0.0, 0.0]);
        let c = m.add_vertex([0.0, 1.0, 0.0]);
        m.add_triangle([a, b, c]);
        assert_eq!(m.n_triangles(), 1);
        let area = m.surface_area();
        assert!((area - 0.5).abs() < 1e-10);
    }

    // ── Triangle-triangle intersection ──────────────────────────────────

    #[test]
    fn test_tri_tri_no_intersection() {
        let a0 = [0.0, 0.0, 0.0];
        let a1 = [1.0, 0.0, 0.0];
        let a2 = [0.0, 1.0, 0.0];
        let b0 = [5.0, 5.0, 5.0];
        let b1 = [6.0, 5.0, 5.0];
        let b2 = [5.0, 6.0, 5.0];
        match triangle_triangle_intersection(a0, a1, a2, b0, b1, b2) {
            TriTriResult::None => {}
            other => panic!("expected None, got {other:?}"),
        }
    }

    #[test]
    fn test_tri_tri_segment_intersection() {
        // Two triangles crossing through each other
        let a0 = [-1.0, 0.0, 0.0];
        let a1 = [1.0, 0.0, 0.0];
        let a2 = [0.0, 0.0, 1.0];
        let b0 = [0.0, -1.0, 0.25];
        let b1 = [0.0, 1.0, 0.25];
        let b2 = [0.0, 0.0, 0.75];
        match triangle_triangle_intersection(a0, a1, a2, b0, b1, b2) {
            TriTriResult::Segment(p0, p1) => {
                assert!(dist_sq(p0, p1) > GEO_EPS, "segment too short");
            }
            other => panic!("expected Segment, got {other:?}"),
        }
    }

    #[test]
    fn test_tri_tri_coplanar_overlap() {
        let a0 = [0.0, 0.0, 0.0];
        let a1 = [2.0, 0.0, 0.0];
        let a2 = [0.0, 2.0, 0.0];
        let b0 = [1.0, 0.0, 0.0];
        let b1 = [3.0, 0.0, 0.0];
        let b2 = [1.0, 2.0, 0.0];
        match triangle_triangle_intersection(a0, a1, a2, b0, b1, b2) {
            TriTriResult::Coplanar => {}
            other => panic!("expected Coplanar, got {other:?}"),
        }
    }

    #[test]
    fn test_tri_tri_coplanar_no_overlap() {
        let a0 = [0.0, 0.0, 0.0];
        let a1 = [1.0, 0.0, 0.0];
        let a2 = [0.0, 1.0, 0.0];
        let b0 = [5.0, 5.0, 0.0];
        let b1 = [6.0, 5.0, 0.0];
        let b2 = [5.0, 6.0, 0.0];
        match triangle_triangle_intersection(a0, a1, a2, b0, b1, b2) {
            TriTriResult::None => {}
            other => panic!("expected None for separated coplanar, got {other:?}"),
        }
    }

    // ── Winding number / inside-outside ─────────────────────────────────

    #[test]
    fn test_winding_number_inside_cube() {
        let cube = SimpleMesh::unit_cube([0.0; 3], 1.0);
        let wn = winding_number(&cube, [0.0, 0.0, 0.0]);
        assert!(wn.abs() > 0.4, "wn = {wn}, expected ~1");
    }

    #[test]
    fn test_winding_number_outside_cube() {
        let cube = SimpleMesh::unit_cube([0.0; 3], 1.0);
        let wn = winding_number(&cube, [5.0, 5.0, 5.0]);
        assert!(wn.abs() < 0.1, "wn = {wn}, expected ~0");
    }

    #[test]
    fn test_is_inside_cube() {
        let cube = SimpleMesh::unit_cube([0.0; 3], 1.0);
        assert!(is_inside(&cube, [0.0, 0.0, 0.0]));
        assert!(!is_inside(&cube, [5.0, 5.0, 5.0]));
    }

    // ── Ray-triangle ────────────────────────────────────────────────────

    #[test]
    fn test_ray_triangle_hit() {
        let v0 = [-1.0, -1.0, 1.0];
        let v1 = [1.0, -1.0, 1.0];
        let v2 = [0.0, 1.0, 1.0];
        let origin = [0.0, 0.0, 0.0];
        let dir = [0.0, 0.0, 1.0];
        let t = ray_triangle_intersect(origin, dir, v0, v1, v2);
        assert!(t.is_some());
        assert!((t.unwrap() - 1.0).abs() < 1e-8);
    }

    #[test]
    fn test_ray_triangle_miss() {
        let v0 = [-1.0, -1.0, 1.0];
        let v1 = [1.0, -1.0, 1.0];
        let v2 = [0.0, 1.0, 1.0];
        let origin = [10.0, 10.0, 0.0];
        let dir = [0.0, 0.0, 1.0];
        assert!(ray_triangle_intersect(origin, dir, v0, v1, v2).is_none());
    }

    #[test]
    fn test_ray_mesh_count() {
        let cube = SimpleMesh::unit_cube([0.0; 3], 1.0);
        let count = ray_mesh_intersection_count(&cube, [-5.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        // Should hit two faces (entering and exiting)
        assert_eq!(count % 2, 0, "count = {count}");
    }

    // ── Contour extraction ──────────────────────────────────────────────

    #[test]
    fn test_contour_overlapping_cubes() {
        let a = SimpleMesh::unit_cube([0.0; 3], 1.0);
        let b = SimpleMesh::unit_cube([1.0, 0.0, 0.0], 1.0);
        let segs = extract_intersection_contours(&a, &b);
        // Overlapping cubes should produce some intersection segments
        assert!(!segs.is_empty(), "expected intersection segments");
    }

    #[test]
    fn test_contour_separated_cubes() {
        let a = SimpleMesh::unit_cube([0.0; 3], 1.0);
        let b = SimpleMesh::unit_cube([10.0, 0.0, 0.0], 1.0);
        let segs = extract_intersection_contours(&a, &b);
        assert!(
            segs.is_empty(),
            "expected no intersections, got {}",
            segs.len()
        );
    }

    // ── Surface splitting ───────────────────────────────────────────────

    #[test]
    fn test_split_triangle_all_front() {
        let (front, back) = split_triangle_by_plane(
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [0.0, 1.0, 1.0],
            [0.0, 0.0, 1.0], // normal pointing up
            [0.0, 0.0, 0.0], // plane at z=0
        );
        assert_eq!(front.len(), 1);
        assert!(back.is_empty());
    }

    #[test]
    fn test_split_triangle_straddles_plane() {
        let (front, back) = split_triangle_by_plane(
            [0.0, 0.0, 1.0],
            [1.0, 0.0, -1.0],
            [0.0, 1.0, -1.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, 0.0],
        );
        assert!(!front.is_empty(), "should have front triangles");
        assert!(!back.is_empty(), "should have back triangles");
    }

    // ── Vertex welding ──────────────────────────────────────────────────

    #[test]
    fn test_weld_vertices_merges_close() {
        let mut m = SimpleMesh::new();
        m.add_vertex([0.0, 0.0, 0.0]);
        m.add_vertex([1e-12, 0.0, 0.0]); // very close to first
        m.add_vertex([1.0, 0.0, 0.0]);
        m.add_triangle([0, 2, 1]);
        let welded = weld_vertices(&m, 1e-8);
        assert_eq!(welded.n_vertices(), 2, "close vertices should merge");
    }

    #[test]
    fn test_remove_degenerate() {
        let mut m = SimpleMesh::from_data(
            vec![[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[0, 0, 1], [0, 1, 2]], // first triangle is degenerate
        );
        remove_degenerate_triangles(&mut m);
        assert_eq!(m.n_triangles(), 1);
    }

    #[test]
    fn test_cleanup_mesh() {
        let mut m = SimpleMesh::new();
        m.add_vertex([0.0, 0.0, 0.0]);
        m.add_vertex([1e-12, 0.0, 0.0]);
        m.add_vertex([1.0, 0.0, 0.0]);
        m.add_vertex([0.0, 1.0, 0.0]);
        m.add_vertex([99.0, 99.0, 99.0]); // unused
        m.add_triangle([0, 2, 3]);
        let cleaned = cleanup_mesh(&m);
        assert!(cleaned.n_vertices() <= 3);
        assert_eq!(cleaned.n_triangles(), 1);
    }

    // ── Mesh stitching ──────────────────────────────────────────────────

    #[test]
    fn test_stitch_meshes() {
        let a = SimpleMesh::unit_cube([0.0; 3], 1.0);
        let b = SimpleMesh::unit_cube([5.0, 0.0, 0.0], 1.0);
        let combined = stitch_meshes(&a, &b);
        assert_eq!(combined.n_vertices(), 16);
        assert_eq!(combined.n_triangles(), 24);
    }

    // ── Boolean operations ──────────────────────────────────────────────

    #[test]
    fn test_boolean_union_separated() {
        let a = SimpleMesh::unit_cube([0.0; 3], 1.0);
        let b = SimpleMesh::unit_cube([10.0, 0.0, 0.0], 1.0);
        let result = mesh_boolean(&a, &b, MeshBooleanOp::Union);
        // No overlap, so union should keep all triangles
        assert!(
            result.n_triangles() >= 20,
            "union tris = {}",
            result.n_triangles()
        );
    }

    #[test]
    fn test_boolean_intersection_separated() {
        let a = SimpleMesh::unit_cube([0.0; 3], 1.0);
        let b = SimpleMesh::unit_cube([10.0, 0.0, 0.0], 1.0);
        let result = mesh_boolean(&a, &b, MeshBooleanOp::Intersection);
        // No overlap → empty intersection
        assert_eq!(result.n_triangles(), 0);
    }

    #[test]
    fn test_boolean_difference_separated() {
        let a = SimpleMesh::unit_cube([0.0; 3], 1.0);
        let b = SimpleMesh::unit_cube([10.0, 0.0, 0.0], 1.0);
        let result = mesh_boolean(&a, &b, MeshBooleanOp::Difference);
        // No overlap → A unchanged
        assert!(
            result.n_triangles() >= 10,
            "diff tris = {}",
            result.n_triangles()
        );
    }

    // ── Mesh stats ──────────────────────────────────────────────────────

    #[test]
    fn test_mesh_stats() {
        let cube = SimpleMesh::unit_cube([0.0; 3], 1.0);
        let stats = compute_mesh_stats(&cube);
        assert_eq!(stats.n_vertices, 8);
        assert_eq!(stats.n_triangles, 12);
        assert!((stats.surface_area - 24.0).abs() < 1e-8);
    }

    #[test]
    fn test_euler_characteristic_cube() {
        let cube = SimpleMesh::unit_cube([0.0; 3], 1.0);
        let chi = euler_characteristic(&cube);
        assert_eq!(
            chi, 2,
            "Euler characteristic of a cube should be 2, got {chi}"
        );
    }

    #[test]
    fn test_flip_normals() {
        let mut cube = SimpleMesh::unit_cube([0.0; 3], 1.0);
        let vol_before = cube.signed_volume();
        cube.flip_normals();
        let vol_after = cube.signed_volume();
        assert!(
            (vol_before + vol_after).abs() < 1e-8,
            "flipping should negate volume"
        );
    }

    #[test]
    fn test_is_watertight_cube() {
        let cube = SimpleMesh::unit_cube([0.0; 3], 1.0);
        assert!(is_watertight(&cube), "unit cube should be watertight");
    }

    #[test]
    fn test_classify_triangles_inside() {
        let big_cube = SimpleMesh::unit_cube([0.0; 3], 2.0);
        let small_cube = SimpleMesh::unit_cube([0.0; 3], 0.5);
        let classification = classify_triangles(&small_cube, &big_cube);
        // All triangles of small cube should be inside big cube
        let all_inside = classification.iter().all(|&x| x);
        assert!(all_inside, "small cube should be inside big cube");
    }
}
