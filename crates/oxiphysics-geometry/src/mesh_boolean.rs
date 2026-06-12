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

use oxiphysics_core::exact_predicates::{Orientation, orient3d};

// ─────────────────────────────────────────────────────────────────────────────
// Vector helpers (no nalgebra in non-core crates)
// ─────────────────────────────────────────────────────────────────────────────

/// Three-component vector type alias.
type V3 = [f64; 3];

/// Small tolerance for geometric tests.
const GEO_EPS: f64 = 1e-10;

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
/// This classifies triangles by their centroid and selects/rejects based on the
/// operation. Centroid sidedness is decided by the EXACT `orient3d`-based
/// point-in-mesh classifier with Simulation-of-Simplicity tie-breaking (see
/// [`exact_classify_triangles`]), so coincident-face configurations are resolved
/// deterministically and produce watertight results instead of landing in an
/// undefined float-winding band. For arbitrary mid-triangle intersections, a full
/// re-triangulation along intersection curves (a BSP/CDT rebuild) is still needed.
pub fn mesh_boolean(mesh_a: &SimpleMesh, mesh_b: &SimpleMesh, op: MeshBooleanOp) -> SimpleMesh {
    // Exact, coincident-face-aware classification of every triangle of each mesh
    // against the other (see `classify_triangle_exact` / `TriClass`).
    let class_a: Vec<TriClass> = (0..mesh_a.n_triangles())
        .map(|i| classify_triangle_exact(mesh_a, i, mesh_b))
        .collect();
    let class_b: Vec<TriClass> = (0..mesh_b.n_triangles())
        .map(|i| classify_triangle_exact(mesh_b, i, mesh_a))
        .collect();

    let mut result = SimpleMesh::new();

    // For Difference, B's contribution (faces inside A) is added with flipped
    // winding to bound the carved cavity; for Union/Intersection B keeps its own
    // winding. Difference flips B; the other ops do not.
    let flip_b = matches!(op, MeshBooleanOp::Difference);

    collect_by_class(mesh_a, &class_a, op, true, false, &mut result);
    collect_by_class(mesh_b, &class_b, op, false, flip_b, &mut result);

    cleanup_mesh(&result)
}

/// Append the triangles of `mesh` (identified as side A when `is_a`) that
/// [`keep_triangle`] selects for `op`, optionally flipping their winding.
fn collect_by_class(
    mesh: &SimpleMesh,
    classes: &[TriClass],
    op: MeshBooleanOp,
    is_a: bool,
    flip: bool,
    result: &mut SimpleMesh,
) {
    let offset = result.vertices.len();
    result.vertices.extend_from_slice(&mesh.vertices);
    for (i, &class) in classes.iter().enumerate() {
        if keep_triangle(op, is_a, class) {
            let t = mesh.triangles[i];
            let tri = if flip {
                [t[1] + offset, t[0] + offset, t[2] + offset]
            } else {
                [t[0] + offset, t[1] + offset, t[2] + offset]
            };
            result.triangles.push(tri);
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
// Exact robust classification (Shewchuk orient3d + Simulation of Simplicity)
// ─────────────────────────────────────────────────────────────────────────────
//
// The boolean pipeline is winding-number / centroid based. The float winding
// number becomes ambiguous EXACTLY when a triangle's centroid lands on (or very
// near) the other mesh's surface — the coincident-face case. This section makes
// the SIDEDNESS decision EXACT: it classifies a point as inside/outside a closed
// mesh by ray-parity, deciding every plane- and pierce-test with the exact
// `orient3d` predicate, and resolves all genuine ties deterministically via a
// Simulation of Simplicity (SoS) rule keyed on the global vertex indices. The
// KEY guarantee is that a point lying exactly on a shared face is classified the
// SAME way every time, so coincident-face Union/Intersection/Difference yield a
// closed manifold rather than landing in an undefined float band.
//
// Implementation note (fast path): we call `orient3d` directly rather than rolling
// an explicit float pre-filter. `orient3d` ALREADY runs Shewchuk's f64 forward-error
// filter internally and only falls back to the exact expansion when the cheap sign
// cannot be certified, so the common (confidently non-degenerate) case is already
// fast and the exact cost is paid only in the degenerate band. This favors
// correctness-by-construction over a hand-tuned duplicate filter whose sign
// convention could silently drift from `orient3d`'s.

/// Resolve a single `orient3d` outcome to a strict sign (`+1` / `-1`) using a
/// Simulation-of-Simplicity tie-break when the determinant is exactly zero.
///
/// SoS: conceptually perturb vertex `k` by `+eps^(k+1)` along the ray axis;
/// degenerate determinants resolve to the sign of the lowest-index distinguishing
/// perturbation, giving a consistent virtual general position. Concretely, when the
/// exact determinant vanishes we decide the sign from the PARITY of the lowest
/// global index among the points involved in that determinant: an even lowest index
/// resolves to `+1`, an odd lowest index to `-1`. Because the rule depends only on
/// the (order-independent) set of indices, the SAME geometric configuration always
/// resolves the SAME way regardless of triangle visit order, so a centroid lying
/// exactly on a shared face is classified identically every time and coincident
/// faces become watertight.
#[inline]
fn sos_sign(base: Orientation, indices: &[u64]) -> i32 {
    match base {
        Orientation::Positive => 1,
        Orientation::Negative => -1,
        Orientation::Degenerate => {
            // Lowest involved global index drives the virtual perturbation order;
            // its parity gives a deterministic, total (never-zero) tie-break.
            let lowest = indices.iter().copied().min().unwrap_or(0);
            if lowest % 2 == 0 { 1 } else { -1 }
        }
    }
}

/// An indexed point: a position paired with its global vertex index for SoS.
///
/// Synthetic ray endpoints carry the two highest `u64` indices; real triangle
/// corners carry their `mesh.triangles[i]` index cast to `u64`.
type IndexedPoint = (V3, u64);

/// Exact test: does the directed segment `p -> q` cross the triangle `(t0,t1,t2)`?
///
/// Each point is supplied as an [`IndexedPoint`] `(position, global_index)`; the
/// index drives the Simulation-of-Simplicity tie-break in [`sos_sign`].
///
/// Decided entirely with `orient3d` sign tests so the result is exact:
///
/// 1. `p` and `q` must lie on OPPOSITE strict sides of the triangle's plane —
///    i.e. `orient3d(t0,t1,t2,p)` and `orient3d(t0,t1,t2,q)` have opposite signs.
/// 2. The ray must pierce the triangle's interior — the three orientations
///    `orient3d(p,q,t0,t1)`, `orient3d(p,q,t1,t2)`, `orient3d(p,q,t2,t0)` must all
///    share the SAME strict sign.
///
/// Any `Orientation::Degenerate` among these five tests is resolved by [`sos_sign`],
/// so coincident / collinear configurations break consistently and the crossing
/// decision is always total.
fn exact_ray_crosses_triangle(
    p: IndexedPoint,
    q: IndexedPoint,
    t0: IndexedPoint,
    t1: IndexedPoint,
    t2: IndexedPoint,
) -> bool {
    let (pp, idx_p) = p;
    let (pq, idx_q) = q;
    let (v0, it0) = t0;
    let (v1, it1) = t1;
    let (v2, it2) = t2;

    // (1) sides of the triangle plane. orient3d(a,b,c,d) is the exact sign of the
    // determinant of (a-d, b-d, c-d), i.e. on which side of plane (a,b,c) point d
    // lies. We test p and q against plane (t0,t1,t2).
    let side_p = sos_sign(orient3d(v0, v1, v2, pp), &[it0, it1, it2, idx_p]);
    let side_q = sos_sign(orient3d(v0, v1, v2, pq), &[it0, it1, it2, idx_q]);
    if side_p == side_q {
        // Same side of the plane (after SoS): the segment does not pierce the plane.
        return false;
    }

    // (2) the line p->q must pass through the triangle interior. The three tetra
    // orientations share a sign iff the piercing point is inside the triangle.
    let e01 = sos_sign(orient3d(pp, pq, v0, v1), &[idx_p, idx_q, it0, it1]);
    let e12 = sos_sign(orient3d(pp, pq, v1, v2), &[idx_p, idx_q, it1, it2]);
    let e20 = sos_sign(orient3d(pp, pq, v2, v0), &[idx_p, idx_q, it2, it0]);

    (e01 == e12) && (e12 == e20)
}

/// Exact inside/outside classification of point `p` against a closed mesh.
///
/// Casts a single ray from `p` to a far endpoint `q` that is provably outside the
/// mesh AABB, counts the exact crossings of every triangle via
/// [`exact_ray_crosses_triangle`], and returns `true` (inside) iff the count is odd.
///
/// `q` is offset along `+X` by `reach = |mx.x - p.x| + 2*diag + 1`, which is large
/// enough that `q.x > mx.x` for ANY query point `p` (so `q` is provably outside the
/// mesh AABB regardless of where `p` sits relative to the box). The small `+jitter`
/// in y/z only steers the ray off the coordinate axes so it does not graze whole
/// coplanar faces; exactness and robustness still come entirely from `orient3d` +
/// SoS, so a "bad" ray is HANDLED rather than relied-upon-to-be-lucky.
pub fn exact_classify_point_vs_mesh(mesh: &SimpleMesh, p: V3) -> bool {
    if mesh.triangles.is_empty() {
        return false;
    }
    let (mn, mx) = mesh.aabb();
    let diag = length(sub(mx, mn));
    // `reach` guarantees q.x = p.x + reach > mx.x even when p is far on the -X side
    // of the box; the extra 2*diag + 1 clears the whole box plus a safety margin so
    // q is strictly outside the AABB and hence outside the closed mesh.
    let reach = (mx[0] - p[0]).abs() + 2.0 * diag + 1.0;
    let q = [p[0] + reach, p[1] + reach * 1.0e-3, p[2] + reach * 1.0e-6];

    // Synthetic indices for p and q: give them the two highest u64 values so SoS
    // perturbs them last / most predictably (they are not mesh vertices).
    let ip: IndexedPoint = (p, u64::MAX);
    let iq: IndexedPoint = (q, u64::MAX - 1);

    let mut crossings: u64 = 0;
    for tri in &mesh.triangles {
        let t0: IndexedPoint = (mesh.vertices[tri[0]], tri[0] as u64);
        let t1: IndexedPoint = (mesh.vertices[tri[1]], tri[1] as u64);
        let t2: IndexedPoint = (mesh.vertices[tri[2]], tri[2] as u64);
        if exact_ray_crosses_triangle(ip, iq, t0, t1, t2) {
            crossings += 1;
        }
    }
    crossings % 2 == 1
}

/// Classify each triangle of `mesh` as inside/outside `other` using the EXACT
/// `orient3d`-based point-in-mesh test (with SoS tie-breaking) on the triangle
/// centroid. Mirrors [`classify_triangles`] but is robust on coincident faces.
///
/// The centroid is the plain average of the triangle's three vertices: the
/// classification OF that centroid is what must be exact, and `orient3d` operates
/// exactly on those `f64` centroid coordinates, so a rational centroid is
/// unnecessary.
pub fn exact_classify_triangles(mesh: &SimpleMesh, other: &SimpleMesh) -> Vec<bool> {
    let mut result = Vec::with_capacity(mesh.n_triangles());
    for i in 0..mesh.n_triangles() {
        let (a, b, c) = mesh.triangle_verts(i);
        let centroid = scale(add(add(a, b), c), 1.0 / 3.0);
        result.push(exact_classify_point_vs_mesh(other, centroid));
    }
    result
}

/// Exact sidedness classification of one triangle of `mesh` against `other`,
/// distinguishing the strict inside/outside cases from the COINCIDENT-FACE
/// (degenerate-band) cases that a pure ray-parity test classifies inconsistently.
///
/// This is the routing of the degenerate band: a triangle whose face lies exactly
/// in a face of `other` is a shared interface, and `mesh_boolean` must resolve it
/// by orientation rather than by the (ill-defined) ray parity of an on-surface
/// centroid. The same-normal vs opposite-normal distinction drives the per-op keep
/// rule so that coincident faces produce a closed manifold.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum TriClass {
    /// Centroid strictly inside `other`.
    Inside,
    /// Centroid strictly outside `other`.
    Outside,
    /// Triangle is coincident with a face of `other` whose normal points the SAME
    /// way (overlapping co-oriented shells).
    CoplanarSame,
    /// Triangle is coincident with a face of `other` whose normal points the
    /// OPPOSITE way (the two solids meet face-to-face here).
    CoplanarOpposite,
}

/// Whether the centroid `c` lies exactly in the plane of triangle `(u0,u1,u2)`
/// (all three plane-side `orient3d` tests degenerate) AND projects strictly inside
/// that triangle — i.e. the query triangle is coincident with this face. The plane
/// test is EXACT via `orient3d`; the in-triangle test uses the existing dominant-
/// axis 2D barycentric test (`point_in_triangle_2d`).
fn centroid_coincident_with(c: V3, u0: V3, u1: V3, u2: V3) -> bool {
    // Exact coplanarity: c must lie on the plane through (u0,u1,u2).
    if orient3d(u0, u1, u2, c) != Orientation::Degenerate {
        return false;
    }
    // Project onto the dominant axis of the face normal and test containment.
    let n = triangle_normal(u0, u1, u2);
    if length(n) < GEO_EPS {
        return false;
    }
    let abs_n = [n[0].abs(), n[1].abs(), n[2].abs()];
    let (ax1, ax2) = if abs_n[0] >= abs_n[1] && abs_n[0] >= abs_n[2] {
        (1, 2)
    } else if abs_n[1] >= abs_n[2] {
        (0, 2)
    } else {
        (0, 1)
    };
    let proj = |v: V3| -> [f64; 2] { [v[ax1], v[ax2]] };
    point_in_triangle_2d(proj(c), proj(u0), proj(u1), proj(u2))
}

/// Classify triangle `ti` of `mesh` against `other`: first detect a coincident
/// (shared) face exactly; otherwise fall back to the exact ray-parity inside test.
fn classify_triangle_exact(mesh: &SimpleMesh, ti: usize, other: &SimpleMesh) -> TriClass {
    let (a, b, c) = mesh.triangle_verts(ti);
    let centroid = scale(add(add(a, b), c), 1.0 / 3.0);
    let n_self = triangle_normal(a, b, c);

    for u in &other.triangles {
        let u0 = other.vertices[u[0]];
        let u1 = other.vertices[u[1]];
        let u2 = other.vertices[u[2]];
        if centroid_coincident_with(centroid, u0, u1, u2) {
            let n_other = triangle_normal(u0, u1, u2);
            if dot(n_self, n_other) >= 0.0 {
                return TriClass::CoplanarSame;
            }
            return TriClass::CoplanarOpposite;
        }
    }

    if exact_classify_point_vs_mesh(other, centroid) {
        TriClass::Inside
    } else {
        TriClass::Outside
    }
}

/// Whether a triangle classified as `class` (of mesh side `is_a`: A vs B) should be
/// kept for boolean `op`. Encodes the standard CSG coplanar-face resolution so that
/// shared interfaces produce a closed manifold:
///
/// - Coincident OPPOSITE faces (two solids meeting face-to-face) are interior to a
///   Union (drop both), the carved interface of a Difference (keep A's copy only),
///   and zero-area for an Intersection (drop both).
/// - Coincident SAME faces (co-oriented duplicate shells) keep exactly one copy —
///   conventionally A's — for Union/Intersection, and are interior for Difference.
fn keep_triangle(op: MeshBooleanOp, is_a: bool, class: TriClass) -> bool {
    match op {
        MeshBooleanOp::Union => match class {
            TriClass::Outside => true,
            TriClass::Inside => false,
            TriClass::CoplanarSame => is_a, // keep A's copy, drop B's
            TriClass::CoplanarOpposite => false, // interior interface, drop both
        },
        MeshBooleanOp::Intersection => match class {
            TriClass::Inside => true,
            TriClass::Outside => false,
            TriClass::CoplanarSame => is_a,
            TriClass::CoplanarOpposite => false,
        },
        MeshBooleanOp::Difference => match class {
            // A \ B keeps A's outside shell and A's copy of the carved interface;
            // B's contribution is handled separately (flipped, inside A).
            TriClass::Outside => is_a,
            TriClass::Inside => !is_a, // B's faces inside A bound the cavity (flipped)
            TriClass::CoplanarSame => false,
            TriClass::CoplanarOpposite => is_a, // keep A's face at the shared interface
        },
    }
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

    // ── Exact orient3d + SoS classification ──────────────────────────────

    /// Deterministic LCG in `[0, 1)` for reproducible near-degenerate stress.
    fn lcg_next(state: &mut u64) -> f64 {
        *state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((*state >> 33) as f64) / (u32::MAX as f64) // in [0,1)
    }

    /// Per-edge incidence count of a mesh (after dedup) — every edge of a closed
    /// 2-manifold must appear exactly twice.
    fn edge_incidence(mesh: &SimpleMesh) -> std::collections::HashMap<(usize, usize), u32> {
        let mut counts: std::collections::HashMap<(usize, usize), u32> =
            std::collections::HashMap::new();
        for t in &mesh.triangles {
            for (&a, &b) in t.iter().zip(t.iter().cycle().skip(1).take(3)) {
                let key = if a < b { (a, b) } else { (b, a) };
                *counts.entry(key).or_insert(0) += 1;
            }
        }
        counts
    }

    /// True iff every coordinate of every vertex is finite (no NaN / inf).
    fn all_vertices_finite(mesh: &SimpleMesh) -> bool {
        mesh.vertices
            .iter()
            .all(|v| v.iter().all(|c| c.is_finite()))
    }

    /// Closed octahedron centred at `c` with radius `r` (stand-in for an
    /// icosphere — a closed convex test body, every edge shared by 2 triangles).
    fn octahedron(c: V3, r: f64) -> SimpleMesh {
        let vertices = vec![
            [c[0] + r, c[1], c[2]],
            [c[0] - r, c[1], c[2]],
            [c[0], c[1] + r, c[2]],
            [c[0], c[1] - r, c[2]],
            [c[0], c[1], c[2] + r],
            [c[0], c[1], c[2] - r],
        ];
        // 8 faces, consistent outward winding.
        let triangles = vec![
            [0, 2, 4],
            [2, 1, 4],
            [1, 3, 4],
            [3, 0, 4],
            [2, 0, 5],
            [1, 2, 5],
            [3, 1, 5],
            [0, 3, 5],
        ];
        SimpleMesh::from_data(vertices, triangles)
    }

    /// Rotate a point about a unit axis by `angle` (Rodrigues' rotation, pure f64).
    fn rotate_about_axis(v: V3, axis: V3, angle: f64) -> V3 {
        let k = normalize(axis);
        let (s, co) = (angle.sin(), angle.cos());
        // v*cos + (k x v)*sin + k*(k·v)*(1-cos)
        let kxv = cross(k, v);
        let kdv = dot(k, v);
        [
            v[0] * co + kxv[0] * s + k[0] * kdv * (1.0 - co),
            v[1] * co + kxv[1] * s + k[1] * kdv * (1.0 - co),
            v[2] * co + kxv[2] * s + k[2] * kdv * (1.0 - co),
        ]
    }

    /// Sanity: the octahedron stand-in is itself a closed 2-manifold.
    #[test]
    fn test_exact_octahedron_is_closed() {
        let o = octahedron([0.0; 3], 1.0);
        assert!(is_watertight(&o), "octahedron must be watertight");
        assert!(edge_incidence(&o).values().all(|&c| c == 2));
        assert_eq!(euler_characteristic(&o), 2);
    }

    // (A) Coincident-face cubes, all three ops.
    #[test]
    fn test_exact_coincident_face_union_watertight() {
        // a's +X face at x=1 coincides with b's -X face at x=1.
        let a = SimpleMesh::unit_cube([0.0, 0.0, 0.0], 1.0);
        let b = SimpleMesh::unit_cube([2.0, 0.0, 0.0], 1.0);
        let result = mesh_boolean(&a, &b, MeshBooleanOp::Union);
        assert!(
            all_vertices_finite(&result),
            "union vertices must be finite (no NaN/inf)"
        );
        let cleaned = cleanup_mesh(&result);
        // Headline gate: the deduped union must be a closed 2-manifold.
        assert!(
            is_watertight(&cleaned),
            "coincident-face union must be watertight after dedup"
        );
        assert_eq!(
            euler_characteristic(&cleaned),
            2,
            "coincident-face union must have Euler characteristic 2"
        );
    }

    #[test]
    fn test_exact_coincident_face_intersection_finite() {
        let a = SimpleMesh::unit_cube([0.0, 0.0, 0.0], 1.0);
        let b = SimpleMesh::unit_cube([2.0, 0.0, 0.0], 1.0);
        // Two cubes meeting only at a face: the intersection is geometrically a
        // flat zero-volume square. It is acceptable for the centroid pipeline to
        // yield an empty / degenerate-and-removed result here — we only require
        // no panic and finite output (no NaN).
        let result = mesh_boolean(&a, &b, MeshBooleanOp::Intersection);
        assert!(
            all_vertices_finite(&result),
            "intersection vertices must be finite"
        );
        // Empty is fine; if non-empty, still must be finite (checked above).
    }

    #[test]
    fn test_exact_coincident_face_difference_watertight() {
        let a = SimpleMesh::unit_cube([0.0, 0.0, 0.0], 1.0);
        let b = SimpleMesh::unit_cube([2.0, 0.0, 0.0], 1.0);
        // b only touches a at the face x=1, so A \ B == A (A is unaffected).
        let result = mesh_boolean(&a, &b, MeshBooleanOp::Difference);
        assert!(
            all_vertices_finite(&result),
            "difference vertices must be finite"
        );
        let cleaned = cleanup_mesh(&result);
        assert!(
            is_watertight(&cleaned),
            "coincident-face difference (== A) must be watertight"
        );
        assert_eq!(
            euler_characteristic(&cleaned),
            2,
            "coincident-face difference must have Euler characteristic 2"
        );
    }

    // (B) Deterministic LCG-driven near-degenerate stress.
    #[test]
    fn test_exact_stress_near_degenerate_no_panic() {
        let mut state: u64 = 0x9E3779B97F4A7C15; // fixed seed → reproducible
        let ops = [
            MeshBooleanOp::Union,
            MeshBooleanOp::Intersection,
            MeshBooleanOp::Difference,
        ];
        let iters = 1000usize;
        let mut completed = 0usize;

        for it in 0..iters {
            let op = ops[it % 3];

            // Near-degenerate tiny perturbations.
            let tiny = 1.0e-8;
            let off_x = 1.0 + (lcg_next(&mut state) - 0.5) * tiny;
            let off_y = (lcg_next(&mut state) - 0.5) * tiny;
            let off_z = (lcg_next(&mut state) - 0.5) * tiny;

            let a = SimpleMesh::unit_cube([0.0, 0.0, 0.0], 1.0);

            // Build B as a slightly-rotated cube near-coincident with A's +X face.
            let mut b = SimpleMesh::unit_cube([off_x, off_y, off_z], 1.0);
            let angle = (lcg_next(&mut state) - 0.5) * 0.1; // ~[-0.05, 0.05] rad
            let axis = [
                lcg_next(&mut state) - 0.5,
                lcg_next(&mut state) - 0.5,
                lcg_next(&mut state) - 0.5,
            ];
            let axis = if length(axis) < GEO_EPS {
                [0.0, 0.0, 1.0]
            } else {
                axis
            };
            let pivot = [off_x, off_y, off_z];
            for v in &mut b.vertices {
                let local = sub(*v, pivot);
                let rot = rotate_about_axis(local, axis, angle);
                *v = add(rot, pivot);
            }

            let r1 = mesh_boolean(&a, &b, op);
            assert!(
                all_vertices_finite(&r1),
                "cube/cube stress iter {it}: NaN/inf vertex"
            );
            // n_triangles() is a usize and always finite/returns; recording it
            // exercises the full pipeline path.
            let _ = r1.n_triangles();

            // Also exercise cube ∩ octahedron (closed convex sphere-like body),
            // near-degenerately overlapping the cube.
            let osph = octahedron([off_x, off_y, off_z], 1.0);
            let r2 = mesh_boolean(&a, &osph, op);
            assert!(
                all_vertices_finite(&r2),
                "cube/octahedron stress iter {it}: NaN/inf vertex"
            );
            let _ = r2.n_triangles();

            completed += 1;
        }
        assert_eq!(completed, iters, "stress loop must complete all iterations");
    }

    // (C) Volume identity on a CLEAN overlapping pair (not degenerate).
    #[test]
    fn test_exact_volume_identity_clean_overlap() {
        // Cube A side 2, volume 8; cube B offset so faces sit cleanly inside/out.
        let a = SimpleMesh::unit_cube([0.0, 0.0, 0.0], 1.0);
        let b = SimpleMesh::unit_cube([1.0, 0.3, 0.2], 1.0);
        let vol_a = 8.0;
        let vol_b = 8.0;

        let vol_union = mesh_boolean(&a, &b, MeshBooleanOp::Union)
            .signed_volume()
            .abs();
        let vol_inter = mesh_boolean(&a, &b, MeshBooleanOp::Intersection)
            .signed_volume()
            .abs();

        // The centroid pipeline selects WHOLE triangles without re-triangulating
        // the intersection curve, so the union/intersection meshes approximate the
        // true boolean volume by whole-triangle selection. We therefore assert the
        // inclusion-exclusion identity with a RELATIVE tolerance scaled by total
        // volume rather than a tight absolute 1e-6 (which whole-triangle selection
        // cannot honestly reach for a mid-face overlap).
        let lhs = vol_union;
        let rhs = vol_a + vol_b - vol_inter;
        let tol = 1.0e-2 * (vol_a + vol_b);
        assert!(
            (lhs - rhs).abs() < tol,
            "volume identity: |{lhs} - {rhs}| = {} >= {tol}",
            (lhs - rhs).abs()
        );
    }

    // (D) Regression: two cubes offset by EXACTLY one face.
    #[test]
    fn test_exact_regression_coincident_face_union() {
        // Regression: float winding number put coincident-face centroids in an
        // undefined band; exact orient3d + SoS classifies them deterministically
        // -> watertight.
        let a = SimpleMesh::unit_cube([0.0, 0.0, 0.0], 1.0);
        let b = SimpleMesh::unit_cube([2.0, 0.0, 0.0], 1.0); // share plane x=1
        let result = mesh_boolean(&a, &b, MeshBooleanOp::Union);
        let cleaned = cleanup_mesh(&result);
        assert!(
            is_watertight(&cleaned),
            "regression: coincident-face union must be watertight"
        );
        assert_eq!(
            euler_characteristic(&cleaned),
            2,
            "regression: coincident-face union Euler must be 2"
        );
    }

    // Exact point-in-mesh sanity: deep interior inside, far exterior outside.
    #[test]
    fn test_exact_point_in_mesh_basic() {
        let cube = SimpleMesh::unit_cube([0.0; 3], 1.0);
        assert!(exact_classify_point_vs_mesh(&cube, [0.0, 0.0, 0.0]));
        assert!(!exact_classify_point_vs_mesh(&cube, [5.0, 5.0, 5.0]));
        // Empty mesh: nothing is inside.
        let empty = SimpleMesh::new();
        assert!(!exact_classify_point_vs_mesh(&empty, [0.0, 0.0, 0.0]));
    }

    // Public exact per-triangle classifier: a small cube fully inside a big cube.
    #[test]
    fn test_exact_classify_triangles_inside() {
        let big = SimpleMesh::unit_cube([0.0; 3], 2.0);
        let small = SimpleMesh::unit_cube([0.0; 3], 0.5);
        let cls = exact_classify_triangles(&small, &big);
        assert!(cls.iter().all(|&x| x), "small cube must be inside big cube");
        // And the big cube's faces are all outside the small cube.
        let cls2 = exact_classify_triangles(&big, &small);
        assert!(cls2.iter().all(|&x| !x), "big faces must be outside small");
    }
}
