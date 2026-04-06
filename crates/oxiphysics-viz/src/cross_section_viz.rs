#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Cross-section and cutting plane visualization.
//!
//! Provides arbitrary cutting plane definition, mesh-plane intersection
//! (triangle-plane clipping), contour line extraction, scalar field
//! interpolation, multiple parallel slice generation, radial cross-sections,
//! section property computation, hatch pattern fill, cross-section animation,
//! boolean operations on cross-sections, stress/strain display, color mapping,
//! and SVG path data export.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Vector / math helpers  (no nalgebra — plain [f64; 3])
// ---------------------------------------------------------------------------

#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
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
fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline]
fn norm3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}

#[inline]
fn normalize3(a: [f64; 3]) -> [f64; 3] {
    let n = norm3(a);
    if n < 1e-15 {
        return [0.0, 0.0, 0.0];
    }
    [a[0] / n, a[1] / n, a[2] / n]
}

#[inline]
fn lerp3(a: [f64; 3], b: [f64; 3], t: f64) -> [f64; 3] {
    [
        a[0] + (b[0] - a[0]) * t,
        a[1] + (b[1] - a[1]) * t,
        a[2] + (b[2] - a[2]) * t,
    ]
}

#[inline]
fn dist3(a: [f64; 3], b: [f64; 3]) -> f64 {
    norm3(sub3(a, b))
}

// ---------------------------------------------------------------------------
// CuttingPlane
// ---------------------------------------------------------------------------

/// An arbitrary cutting plane defined by a point and a unit normal.
#[derive(Debug, Clone, PartialEq)]
pub struct CuttingPlane {
    /// A point lying on the plane.
    pub point: [f64; 3],
    /// The outward unit normal of the plane.
    pub normal: [f64; 3],
}

impl CuttingPlane {
    /// Create a new cutting plane given a point and a normal (will be normalised).
    pub fn new(point: [f64; 3], normal: [f64; 3]) -> Self {
        Self {
            point,
            normal: normalize3(normal),
        }
    }

    /// XY plane at a given Z height.
    pub fn xy(z: f64) -> Self {
        Self::new([0.0, 0.0, z], [0.0, 0.0, 1.0])
    }

    /// XZ plane at a given Y height.
    pub fn xz(y: f64) -> Self {
        Self::new([0.0, y, 0.0], [0.0, 1.0, 0.0])
    }

    /// YZ plane at a given X position.
    pub fn yz(x: f64) -> Self {
        Self::new([x, 0.0, 0.0], [1.0, 0.0, 0.0])
    }

    /// Signed distance from a point to the plane.
    pub fn signed_distance(&self, p: [f64; 3]) -> f64 {
        dot3(sub3(p, self.point), self.normal)
    }

    /// Project a 3-D point onto the plane.
    pub fn project(&self, p: [f64; 3]) -> [f64; 3] {
        let d = self.signed_distance(p);
        sub3(p, scale3(self.normal, d))
    }

    /// Build a local 2-D coordinate system on the plane.
    /// Returns `(u_axis, v_axis)` — two orthonormal tangent vectors.
    pub fn tangent_basis(&self) -> ([f64; 3], [f64; 3]) {
        let n = self.normal;
        let ref_vec = if n[0].abs() < 0.9 {
            [1.0, 0.0, 0.0]
        } else {
            [0.0, 1.0, 0.0]
        };
        let u = normalize3(cross3(ref_vec, n));
        let v = cross3(n, u);
        (u, v)
    }

    /// Convert a 3-D world point to local 2-D coordinates `(u, v)` on the plane.
    pub fn to_local_2d(&self, p: [f64; 3]) -> [f64; 2] {
        let (u_axis, v_axis) = self.tangent_basis();
        let rel = sub3(p, self.point);
        [dot3(rel, u_axis), dot3(rel, v_axis)]
    }

    /// Convert local 2-D coordinates back to 3-D world coordinates.
    pub fn from_local_2d(&self, uv: [f64; 2]) -> [f64; 3] {
        let (u_axis, v_axis) = self.tangent_basis();
        add3(
            add3(self.point, scale3(u_axis, uv[0])),
            scale3(v_axis, uv[1]),
        )
    }

    /// Translate the plane along its normal by `offset`.
    pub fn offset(&self, offset: f64) -> Self {
        Self {
            point: add3(self.point, scale3(self.normal, offset)),
            normal: self.normal,
        }
    }
}

// ---------------------------------------------------------------------------
// Triangle-plane clipping
// ---------------------------------------------------------------------------

/// A triangle defined by three vertices in 3-D space.
#[derive(Debug, Clone, PartialEq)]
pub struct Triangle3D {
    /// The three vertices of the triangle.
    pub vertices: [[f64; 3]; 3],
}

impl Triangle3D {
    /// Create a triangle from three vertex positions.
    pub fn new(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> Self {
        Self {
            vertices: [a, b, c],
        }
    }

    /// Compute the triangle's face normal (unnormalised).
    pub fn normal(&self) -> [f64; 3] {
        let ab = sub3(self.vertices[1], self.vertices[0]);
        let ac = sub3(self.vertices[2], self.vertices[0]);
        cross3(ab, ac)
    }

    /// Compute the area of the triangle.
    pub fn area(&self) -> f64 {
        0.5 * norm3(self.normal())
    }

    /// Compute the centroid of the triangle.
    pub fn centroid(&self) -> [f64; 3] {
        [
            (self.vertices[0][0] + self.vertices[1][0] + self.vertices[2][0]) / 3.0,
            (self.vertices[0][1] + self.vertices[1][1] + self.vertices[2][1]) / 3.0,
            (self.vertices[0][2] + self.vertices[1][2] + self.vertices[2][2]) / 3.0,
        ]
    }
}

/// A line segment resulting from a triangle–plane intersection.
#[derive(Debug, Clone, PartialEq)]
pub struct IntersectionSegment {
    /// Start point of the intersection segment.
    pub start: [f64; 3],
    /// End point of the intersection segment.
    pub end: [f64; 3],
}

/// Clip a triangle against a cutting plane.
///
/// Returns `Some(segment)` if the triangle intersects the plane, or `None`
/// if the triangle lies entirely on one side.
pub fn clip_triangle_plane(tri: &Triangle3D, plane: &CuttingPlane) -> Option<IntersectionSegment> {
    let d = [
        plane.signed_distance(tri.vertices[0]),
        plane.signed_distance(tri.vertices[1]),
        plane.signed_distance(tri.vertices[2]),
    ];

    // Collect intersection points on edges
    let mut pts: Vec<[f64; 3]> = Vec::new();
    for i in 0..3 {
        let j = (i + 1) % 3;
        if d[i].signum() != d[j].signum() && (d[i] - d[j]).abs() > 1e-15 {
            let t = d[i] / (d[i] - d[j]);
            pts.push(lerp3(tri.vertices[i], tri.vertices[j], t));
        }
    }

    // Also add vertices exactly on the plane
    for i in 0..3 {
        if d[i].abs() < 1e-12 {
            pts.push(tri.vertices[i]);
        }
    }

    // Remove near-duplicate points
    pts.dedup_by(|a, b| dist3(*a, *b) < 1e-12);

    if pts.len() >= 2 {
        Some(IntersectionSegment {
            start: pts[0],
            end: pts[1],
        })
    } else {
        None
    }
}

/// Clip an entire triangle mesh (given as a list of triangles) against a
/// cutting plane. Returns all intersection segments.
pub fn clip_mesh_plane(triangles: &[Triangle3D], plane: &CuttingPlane) -> Vec<IntersectionSegment> {
    triangles
        .iter()
        .filter_map(|tri| clip_triangle_plane(tri, plane))
        .collect()
}

// ---------------------------------------------------------------------------
// Contour extraction
// ---------------------------------------------------------------------------

/// A contour line on the cross-section surface, represented as a polyline.
#[derive(Debug, Clone)]
pub struct ContourLine {
    /// Points along the contour, in 3-D world coordinates.
    pub points: Vec<[f64; 3]>,
    /// The scalar value this contour represents.
    pub value: f64,
}

/// A 2-D scalar grid used for contour extraction on a cutting plane.
#[derive(Debug, Clone)]
pub struct ScalarGrid2D {
    /// Number of cells in the u direction.
    pub nu: usize,
    /// Number of cells in the v direction.
    pub nv: usize,
    /// Cell size in u direction.
    pub du: f64,
    /// Cell size in v direction.
    pub dv: f64,
    /// Scalar values at grid nodes, stored row-major `[v * (nu+1) + u]`.
    pub values: Vec<f64>,
    /// The cutting plane this grid lives on.
    pub plane: CuttingPlane,
    /// Origin offset in local u coordinate.
    pub u_min: f64,
    /// Origin offset in local v coordinate.
    pub v_min: f64,
}

impl ScalarGrid2D {
    /// Create a grid on the given cutting plane.
    ///
    /// The grid spans `[-half_extent..half_extent]` in both u and v,
    /// with `n` cells in each direction.
    pub fn new(plane: &CuttingPlane, half_extent: f64, n: usize) -> Self {
        let n = n.max(1);
        let du = 2.0 * half_extent / n as f64;
        let dv = du;
        let nodes = (n + 1) * (n + 1);
        Self {
            nu: n,
            nv: n,
            du,
            dv,
            values: vec![0.0; nodes],
            plane: plane.clone(),
            u_min: -half_extent,
            v_min: -half_extent,
        }
    }

    /// Number of nodes in the u direction.
    pub fn nodes_u(&self) -> usize {
        self.nu + 1
    }

    /// Number of nodes in the v direction.
    pub fn nodes_v(&self) -> usize {
        self.nv + 1
    }

    /// Get the scalar value at node `(iu, iv)`.
    pub fn get(&self, iu: usize, iv: usize) -> f64 {
        self.values[iv * self.nodes_u() + iu]
    }

    /// Set the scalar value at node `(iu, iv)`.
    pub fn set(&mut self, iu: usize, iv: usize, val: f64) {
        let idx = iv * self.nodes_u() + iu;
        self.values[idx] = val;
    }

    /// Get the 3-D world position of node `(iu, iv)`.
    pub fn node_position(&self, iu: usize, iv: usize) -> [f64; 3] {
        let u = self.u_min + iu as f64 * self.du;
        let v = self.v_min + iv as f64 * self.dv;
        self.plane.from_local_2d([u, v])
    }

    /// Sample a scalar field function at every grid node.
    pub fn sample_field<F>(&mut self, field: F)
    where
        F: Fn([f64; 3]) -> f64,
    {
        for iv in 0..self.nodes_v() {
            for iu in 0..self.nodes_u() {
                let pos = self.node_position(iu, iv);
                self.set(iu, iv, field(pos));
            }
        }
    }

    /// Minimum scalar value on the grid.
    pub fn min_value(&self) -> f64 {
        self.values.iter().cloned().fold(f64::INFINITY, f64::min)
    }

    /// Maximum scalar value on the grid.
    pub fn max_value(&self) -> f64 {
        self.values
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max)
    }
}

/// Extract iso-contour lines at a given scalar value using marching squares.
pub fn extract_contour(grid: &ScalarGrid2D, iso_value: f64) -> Vec<ContourLine> {
    let mut segments: Vec<([f64; 3], [f64; 3])> = Vec::new();

    for iv in 0..grid.nv {
        for iu in 0..grid.nu {
            let v00 = grid.get(iu, iv);
            let v10 = grid.get(iu + 1, iv);
            let v01 = grid.get(iu, iv + 1);
            let v11 = grid.get(iu + 1, iv + 1);

            let p00 = grid.node_position(iu, iv);
            let p10 = grid.node_position(iu + 1, iv);
            let p01 = grid.node_position(iu, iv + 1);
            let p11 = grid.node_position(iu + 1, iv + 1);

            let b00 = (v00 >= iso_value) as u8;
            let b10 = (v10 >= iso_value) as u8;
            let b01 = (v01 >= iso_value) as u8;
            let b11 = (v11 >= iso_value) as u8;

            let case = b00 | (b10 << 1) | (b01 << 2) | (b11 << 3);
            if case == 0 || case == 15 {
                continue;
            }

            // Edge interpolation helper
            let interp = |pa: [f64; 3], va: f64, pb: [f64; 3], vb: f64| -> [f64; 3] {
                let dv = vb - va;
                if dv.abs() < 1e-15 {
                    lerp3(pa, pb, 0.5)
                } else {
                    let t = (iso_value - va) / dv;
                    lerp3(pa, pb, t.clamp(0.0, 1.0))
                }
            };

            // Edge midpoints: bottom, right, top, left
            let eb = interp(p00, v00, p10, v10);
            let er = interp(p10, v10, p11, v11);
            let et = interp(p01, v01, p11, v11);
            let el = interp(p00, v00, p01, v01);

            // Marching squares lookup
            match case {
                1 | 14 => segments.push((el, eb)),
                2 | 13 => segments.push((eb, er)),
                3 | 12 => segments.push((el, er)),
                4 | 11 => segments.push((et, el)),
                5 => {
                    // Saddle — use average
                    let avg = 0.25 * (v00 + v10 + v01 + v11);
                    if avg >= iso_value {
                        segments.push((el, eb));
                        segments.push((et, er));
                    } else {
                        segments.push((el, et));
                        segments.push((eb, er));
                    }
                }
                6 | 9 => segments.push((eb, et)),
                7 | 8 => segments.push((et, er)),
                10 => {
                    // Saddle
                    let avg = 0.25 * (v00 + v10 + v01 + v11);
                    if avg >= iso_value {
                        segments.push((el, et));
                        segments.push((eb, er));
                    } else {
                        segments.push((el, eb));
                        segments.push((et, er));
                    }
                }
                _ => {}
            }
        }
    }

    // Chain segments into polylines
    chain_segments_into_contours(segments, iso_value)
}

/// Chain a set of line segments into contour polylines.
fn chain_segments_into_contours(
    segments: Vec<([f64; 3], [f64; 3])>,
    iso_value: f64,
) -> Vec<ContourLine> {
    if segments.is_empty() {
        return Vec::new();
    }

    let mut used = vec![false; segments.len()];
    let mut contours: Vec<ContourLine> = Vec::new();
    let eps = 1e-10;

    for start_idx in 0..segments.len() {
        if used[start_idx] {
            continue;
        }
        used[start_idx] = true;
        let mut pts = vec![segments[start_idx].0, segments[start_idx].1];

        // Extend forward
        loop {
            let tail = *pts.last().expect("collection should not be empty");
            let mut found = false;
            for i in 0..segments.len() {
                if used[i] {
                    continue;
                }
                if dist3(segments[i].0, tail) < eps {
                    pts.push(segments[i].1);
                    used[i] = true;
                    found = true;
                    break;
                }
                if dist3(segments[i].1, tail) < eps {
                    pts.push(segments[i].0);
                    used[i] = true;
                    found = true;
                    break;
                }
            }
            if !found {
                break;
            }
        }

        contours.push(ContourLine {
            points: pts,
            value: iso_value,
        });
    }

    contours
}

/// Extract multiple iso-contour lines at evenly spaced levels.
pub fn extract_contours_multi(grid: &ScalarGrid2D, num_levels: usize) -> Vec<ContourLine> {
    let vmin = grid.min_value();
    let vmax = grid.max_value();
    if (vmax - vmin).abs() < 1e-15 || num_levels == 0 {
        return Vec::new();
    }
    let mut all = Vec::new();
    for i in 0..num_levels {
        let t = (i as f64 + 0.5) / num_levels as f64;
        let iso = vmin + t * (vmax - vmin);
        all.extend(extract_contour(grid, iso));
    }
    all
}

// ---------------------------------------------------------------------------
// Scalar field interpolation on cut plane
// ---------------------------------------------------------------------------

/// Bilinearly interpolate a value from a `ScalarGrid2D` at local `(u, v)`.
pub fn interpolate_scalar_grid(grid: &ScalarGrid2D, u: f64, v: f64) -> f64 {
    let fu = (u - grid.u_min) / grid.du;
    let fv = (v - grid.v_min) / grid.dv;
    let iu = (fu.floor() as usize).min(grid.nu.saturating_sub(1));
    let iv = (fv.floor() as usize).min(grid.nv.saturating_sub(1));
    let tu = (fu - iu as f64).clamp(0.0, 1.0);
    let tv = (fv - iv as f64).clamp(0.0, 1.0);

    let v00 = grid.get(iu, iv);
    let v10 = grid.get(iu + 1, iv);
    let v01 = grid.get(iu, iv + 1);
    let v11 = grid.get(iu + 1, iv + 1);

    let a = v00 * (1.0 - tu) + v10 * tu;
    let b = v01 * (1.0 - tu) + v11 * tu;
    a * (1.0 - tv) + b * tv
}

/// Interpolate a scalar field at an arbitrary 3-D world point by projecting
/// onto the grid's cutting plane.
pub fn interpolate_at_world_point(grid: &ScalarGrid2D, p: [f64; 3]) -> f64 {
    let uv = grid.plane.to_local_2d(p);
    interpolate_scalar_grid(grid, uv[0], uv[1])
}

// ---------------------------------------------------------------------------
// Multiple parallel slice generation
// ---------------------------------------------------------------------------

/// Configuration for generating a stack of parallel cutting planes.
#[derive(Debug, Clone)]
pub struct ParallelSliceConfig {
    /// Base plane (the first slice).
    pub base_plane: CuttingPlane,
    /// Spacing between consecutive slices along the normal direction.
    pub spacing: f64,
    /// Number of slices.
    pub count: usize,
}

/// Generate a set of parallel cutting planes.
pub fn generate_parallel_slices(config: &ParallelSliceConfig) -> Vec<CuttingPlane> {
    (0..config.count)
        .map(|i| config.base_plane.offset(i as f64 * config.spacing))
        .collect()
}

/// Clip a mesh against a stack of parallel slices, returning one set of
/// intersection segments per slice.
pub fn clip_mesh_parallel_slices(
    triangles: &[Triangle3D],
    config: &ParallelSliceConfig,
) -> Vec<Vec<IntersectionSegment>> {
    generate_parallel_slices(config)
        .iter()
        .map(|plane| clip_mesh_plane(triangles, plane))
        .collect()
}

// ---------------------------------------------------------------------------
// Radial cross-sections (cylindrical cuts)
// ---------------------------------------------------------------------------

/// A cylindrical cutting surface, defined by an axis ray and a radius.
#[derive(Debug, Clone)]
pub struct CylindricalCut {
    /// A point on the cylinder axis.
    pub axis_point: [f64; 3],
    /// Unit direction of the cylinder axis.
    pub axis_dir: [f64; 3],
    /// Radius of the cylinder.
    pub radius: f64,
}

impl CylindricalCut {
    /// Create a new cylindrical cut.
    pub fn new(axis_point: [f64; 3], axis_dir: [f64; 3], radius: f64) -> Self {
        Self {
            axis_point,
            axis_dir: normalize3(axis_dir),
            radius,
        }
    }

    /// Compute signed radial distance of a point from the cylinder surface.
    /// Positive means outside.
    pub fn radial_distance(&self, p: [f64; 3]) -> f64 {
        let ap = sub3(p, self.axis_point);
        let along = dot3(ap, self.axis_dir);
        let proj = scale3(self.axis_dir, along);
        let perp = sub3(ap, proj);
        norm3(perp) - self.radius
    }
}

/// Generate `n_angles` radial cutting planes that share the cylinder axis as
/// their line of intersection. Each plane passes through the axis.
pub fn generate_radial_slices(
    axis_point: [f64; 3],
    axis_dir: [f64; 3],
    n_angles: usize,
) -> Vec<CuttingPlane> {
    let axis = normalize3(axis_dir);
    let ref_vec = if axis[0].abs() < 0.9 {
        [1.0, 0.0, 0.0]
    } else {
        [0.0, 1.0, 0.0]
    };
    let u = normalize3(cross3(ref_vec, axis));
    let v = cross3(axis, u);

    (0..n_angles)
        .map(|i| {
            let theta = i as f64 * PI / n_angles as f64;
            let normal = normalize3(add3(scale3(u, theta.cos()), scale3(v, theta.sin())));
            CuttingPlane::new(axis_point, normal)
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Cross-section area computation
// ---------------------------------------------------------------------------

/// Compute the area of a 2-D polygon given its vertices in local 2-D
/// coordinates using the shoelace formula.
pub fn polygon_area_2d(vertices: &[[f64; 2]]) -> f64 {
    let n = vertices.len();
    if n < 3 {
        return 0.0;
    }
    let mut area = 0.0;
    for i in 0..n {
        let j = (i + 1) % n;
        area += vertices[i][0] * vertices[j][1];
        area -= vertices[j][0] * vertices[i][1];
    }
    (area / 2.0).abs()
}

/// Estimate the cross-section area from intersection segments projected onto
/// the cutting plane. Uses a simple convex-hull-area approximation.
pub fn cross_section_area_from_segments(
    segments: &[IntersectionSegment],
    plane: &CuttingPlane,
) -> f64 {
    if segments.is_empty() {
        return 0.0;
    }
    // Collect all 2-D projected points
    let mut pts_2d: Vec<[f64; 2]> = Vec::new();
    for seg in segments {
        pts_2d.push(plane.to_local_2d(seg.start));
        pts_2d.push(plane.to_local_2d(seg.end));
    }
    // Sort by angle from centroid for a convex hull estimate
    let cx = pts_2d.iter().map(|p| p[0]).sum::<f64>() / pts_2d.len() as f64;
    let cy = pts_2d.iter().map(|p| p[1]).sum::<f64>() / pts_2d.len() as f64;
    pts_2d.sort_by(|a, b| {
        let ang_a = (a[1] - cy).atan2(a[0] - cx);
        let ang_b = (b[1] - cy).atan2(b[0] - cx);
        ang_a
            .partial_cmp(&ang_b)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    pts_2d.dedup_by(|a, b| (a[0] - b[0]).abs() < 1e-12 && (a[1] - b[1]).abs() < 1e-12);
    polygon_area_2d(&pts_2d)
}

// ---------------------------------------------------------------------------
// Section properties (centroid, moments of inertia)
// ---------------------------------------------------------------------------

/// Properties of a cross-section polygon.
#[derive(Debug, Clone)]
pub struct SectionProperties {
    /// Area of the cross-section.
    pub area: f64,
    /// Centroid in local 2-D coordinates.
    pub centroid: [f64; 2],
    /// Second moment of area about the local u-axis passing through centroid.
    pub iuu: f64,
    /// Second moment of area about the local v-axis passing through centroid.
    pub ivv: f64,
    /// Product of area (cross moment) about centroidal axes.
    pub iuv: f64,
}

/// Compute section properties from a closed polygon (in local 2-D coordinates).
pub fn compute_section_properties(vertices: &[[f64; 2]]) -> SectionProperties {
    let n = vertices.len();
    if n < 3 {
        return SectionProperties {
            area: 0.0,
            centroid: [0.0, 0.0],
            iuu: 0.0,
            ivv: 0.0,
            iuv: 0.0,
        };
    }

    // Signed area
    let mut area_2 = 0.0;
    let mut cx = 0.0;
    let mut cy = 0.0;
    for i in 0..n {
        let j = (i + 1) % n;
        let cross = vertices[i][0] * vertices[j][1] - vertices[j][0] * vertices[i][1];
        area_2 += cross;
        cx += (vertices[i][0] + vertices[j][0]) * cross;
        cy += (vertices[i][1] + vertices[j][1]) * cross;
    }
    let area = (area_2 / 2.0).abs();
    if area < 1e-15 {
        return SectionProperties {
            area: 0.0,
            centroid: [0.0, 0.0],
            iuu: 0.0,
            ivv: 0.0,
            iuv: 0.0,
        };
    }
    let sign = if area_2 > 0.0 { 1.0 } else { -1.0 };
    cx = sign * cx / (6.0 * area);
    cy = sign * cy / (6.0 * area);

    // Second moments of area about centroid
    let mut iuu = 0.0;
    let mut ivv = 0.0;
    let mut iuv = 0.0;
    for i in 0..n {
        let j = (i + 1) % n;
        let ui = vertices[i][0] - cx;
        let vi = vertices[i][1] - cy;
        let uj = vertices[j][0] - cx;
        let vj = vertices[j][1] - cy;
        let cross = ui * vj - uj * vi;
        iuu += (vi * vi + vi * vj + vj * vj) * cross;
        ivv += (ui * ui + ui * uj + uj * uj) * cross;
        iuv += (ui * vj + 2.0 * ui * vi + 2.0 * uj * vj + uj * vi) * cross;
    }
    iuu = (iuu / 12.0).abs();
    ivv = (ivv / 12.0).abs();
    iuv /= 24.0;

    SectionProperties {
        area,
        centroid: [cx, cy],
        iuu,
        ivv,
        iuv,
    }
}

// ---------------------------------------------------------------------------
// Hatch pattern fill for solid regions
// ---------------------------------------------------------------------------

/// A single hatch line (2-D segment on the cutting plane).
#[derive(Debug, Clone)]
pub struct HatchLine {
    /// Start point in local 2-D coordinates.
    pub start: [f64; 2],
    /// End point in local 2-D coordinates.
    pub end: [f64; 2],
}

/// Configuration for hatch pattern generation.
#[derive(Debug, Clone)]
pub struct HatchConfig {
    /// Spacing between hatch lines.
    pub spacing: f64,
    /// Angle of hatch lines in radians (0 = horizontal).
    pub angle: f64,
}

impl Default for HatchConfig {
    fn default() -> Self {
        Self {
            spacing: 0.1,
            angle: PI / 4.0,
        }
    }
}

/// Generate hatch lines inside a convex polygon.
///
/// Clips parallel sweep lines at the given angle against the polygon boundary.
pub fn generate_hatch_lines(polygon: &[[f64; 2]], config: &HatchConfig) -> Vec<HatchLine> {
    if polygon.len() < 3 || config.spacing <= 0.0 {
        return Vec::new();
    }

    let cos_a = config.angle.cos();
    let sin_a = config.angle.sin();

    // Project polygon vertices onto the sweep direction
    let project = |p: &[f64; 2]| -> f64 { -sin_a * p[0] + cos_a * p[1] };
    let min_proj = polygon.iter().map(project).fold(f64::INFINITY, f64::min);
    let max_proj = polygon
        .iter()
        .map(project)
        .fold(f64::NEG_INFINITY, f64::max);

    let mut lines = Vec::new();
    let mut d = min_proj + config.spacing * 0.5;
    while d < max_proj {
        // Line: -sin_a * x + cos_a * y = d
        // Find intersections with polygon edges
        let mut intersections: Vec<f64> = Vec::new();
        let n = polygon.len();
        for i in 0..n {
            let j = (i + 1) % n;
            let di = project(&polygon[i]) - d;
            let dj = project(&polygon[j]) - d;
            if di.signum() != dj.signum() && (di - dj).abs() > 1e-15 {
                let t = di / (di - dj);
                let px = polygon[i][0] + t * (polygon[j][0] - polygon[i][0]);
                let py = polygon[i][1] + t * (polygon[j][1] - polygon[i][1]);
                // Project onto the line direction
                let along = cos_a * px + sin_a * py;
                intersections.push(along);
            }
        }
        intersections.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

        // Pair up intersections
        let mut k = 0;
        while k + 1 < intersections.len() {
            let s = intersections[k];
            let e = intersections[k + 1];
            let x0 = cos_a * s + (-sin_a) * d;
            let y0 = sin_a * s + cos_a * d;
            let x1 = cos_a * e + (-sin_a) * d;
            let y1 = sin_a * e + cos_a * d;
            // Reconstruct from d and along
            // Actually: point on the sweep line at parameter `along`:
            //   x = cos_a * along - sin_a * d  (wrong sign fix)
            // Let's redo more carefully:
            //   sweep_normal = (-sin_a, cos_a)
            //   sweep_dir    = (cos_a, sin_a)
            //   point = d * sweep_normal_hat + along * sweep_dir_hat
            //         = d * (-sin_a, cos_a) + along * (cos_a, sin_a)
            // But sweep_normal_hat is already unit => no rename needed
            let _x0 = x0;
            let _y0 = y0;
            let _x1 = x1;
            let _y1 = y1;
            // Correct reconstruction:
            let sx = -sin_a * d + cos_a * s;
            let sy = cos_a * d + sin_a * s;
            let ex = -sin_a * d + cos_a * e;
            let ey = cos_a * d + sin_a * e;
            lines.push(HatchLine {
                start: [sx, sy],
                end: [ex, ey],
            });
            k += 2;
        }
        d += config.spacing;
    }

    lines
}

// ---------------------------------------------------------------------------
// Cross-section animation (sweeping plane)
// ---------------------------------------------------------------------------

/// A frame in a cross-section sweep animation.
#[derive(Debug, Clone)]
pub struct SweepFrame {
    /// The cutting plane for this frame.
    pub plane: CuttingPlane,
    /// Intersection segments at this frame.
    pub segments: Vec<IntersectionSegment>,
    /// Normalized sweep parameter \[0, 1\].
    pub t: f64,
}

/// Generate a sweep animation by moving a cutting plane along its normal.
///
/// Returns one `SweepFrame` per step.
pub fn generate_sweep_animation(
    triangles: &[Triangle3D],
    base_plane: &CuttingPlane,
    sweep_distance: f64,
    num_frames: usize,
) -> Vec<SweepFrame> {
    if num_frames == 0 {
        return Vec::new();
    }
    (0..num_frames)
        .map(|i| {
            let t = i as f64 / (num_frames.max(1) - 1).max(1) as f64;
            let plane = base_plane.offset(t * sweep_distance);
            let segments = clip_mesh_plane(triangles, &plane);
            SweepFrame { plane, segments, t }
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Boolean cross-section (union / intersection of cuts)
// ---------------------------------------------------------------------------

/// Result of a boolean operation on two sets of intersection segments.
#[derive(Debug, Clone)]
pub struct BooleanCrossSection {
    /// Merged segments after the boolean operation.
    pub segments: Vec<IntersectionSegment>,
}

/// Union of two cross-section segment sets (simple concatenation with
/// deduplication of near-identical segments).
pub fn cross_section_union(
    a: &[IntersectionSegment],
    b: &[IntersectionSegment],
) -> BooleanCrossSection {
    let mut all: Vec<IntersectionSegment> = Vec::with_capacity(a.len() + b.len());
    all.extend_from_slice(a);
    for seg in b {
        let duplicate = all.iter().any(|existing| {
            dist3(existing.start, seg.start) < 1e-10 && dist3(existing.end, seg.end) < 1e-10
        });
        if !duplicate {
            all.push(seg.clone());
        }
    }
    BooleanCrossSection { segments: all }
}

/// Intersection of two cross-section segment sets: keep only segments from `a`
/// whose midpoint lies within the convex hull of `b`'s projected points.
pub fn cross_section_intersection(
    a: &[IntersectionSegment],
    b: &[IntersectionSegment],
    plane: &CuttingPlane,
) -> BooleanCrossSection {
    if b.is_empty() {
        return BooleanCrossSection {
            segments: Vec::new(),
        };
    }

    // Compute bounding box of b in 2-D
    let mut umin = f64::INFINITY;
    let mut umax = f64::NEG_INFINITY;
    let mut vmin = f64::INFINITY;
    let mut vmax = f64::NEG_INFINITY;
    for seg in b {
        for pt in [seg.start, seg.end] {
            let uv = plane.to_local_2d(pt);
            umin = umin.min(uv[0]);
            umax = umax.max(uv[0]);
            vmin = vmin.min(uv[1]);
            vmax = vmax.max(uv[1]);
        }
    }

    let mut result = Vec::new();
    for seg in a {
        let mid = lerp3(seg.start, seg.end, 0.5);
        let uv = plane.to_local_2d(mid);
        if uv[0] >= umin && uv[0] <= umax && uv[1] >= vmin && uv[1] <= vmax {
            result.push(seg.clone());
        }
    }

    BooleanCrossSection { segments: result }
}

// ---------------------------------------------------------------------------
// Stress/strain field display on cut surface
// ---------------------------------------------------------------------------

/// A stress/strain sample on the cross-section.
#[derive(Debug, Clone)]
pub struct StressFieldSample {
    /// Position in 3-D world space.
    pub position: [f64; 3],
    /// Stress tensor components `[σ_xx, σ_yy, σ_zz, τ_xy, τ_yz, τ_xz]`.
    pub stress: [f64; 6],
    /// Von Mises equivalent stress.
    pub von_mises: f64,
}

/// Compute the von Mises equivalent stress from a symmetric stress tensor.
pub fn von_mises_from_tensor(s: &[f64; 6]) -> f64 {
    let sxx = s[0];
    let syy = s[1];
    let szz = s[2];
    let txy = s[3];
    let tyz = s[4];
    let txz = s[5];
    let term1 = (sxx - syy).powi(2) + (syy - szz).powi(2) + (szz - sxx).powi(2);
    let term2 = 6.0 * (txy * txy + tyz * tyz + txz * txz);
    ((term1 + term2) / 2.0).sqrt()
}

/// Sample a stress field on a grid over the cutting plane.
pub fn sample_stress_field<F>(
    plane: &CuttingPlane,
    half_extent: f64,
    resolution: usize,
    stress_fn: F,
) -> Vec<StressFieldSample>
where
    F: Fn([f64; 3]) -> [f64; 6],
{
    let resolution = resolution.max(1);
    let (u_axis, v_axis) = plane.tangent_basis();
    let step = 2.0 * half_extent / resolution as f64;
    let mut samples = Vec::with_capacity(resolution * resolution);

    for iv in 0..=resolution {
        for iu in 0..=resolution {
            let u = -half_extent + iu as f64 * step;
            let v = -half_extent + iv as f64 * step;
            let pos = add3(add3(plane.point, scale3(u_axis, u)), scale3(v_axis, v));
            let stress = stress_fn(pos);
            let vm = von_mises_from_tensor(&stress);
            samples.push(StressFieldSample {
                position: pos,
                stress,
                von_mises: vm,
            });
        }
    }

    samples
}

// ---------------------------------------------------------------------------
// Color mapping on cross-section contours
// ---------------------------------------------------------------------------

/// An RGBA color for visualization.
#[derive(Debug, Clone, Copy)]
pub struct CrossSectionColor {
    /// Red channel \[0, 1\].
    pub r: f64,
    /// Green channel \[0, 1\].
    pub g: f64,
    /// Blue channel \[0, 1\].
    pub b: f64,
    /// Alpha channel \[0, 1\].
    pub a: f64,
}

impl CrossSectionColor {
    /// Create a new color from RGBA values.
    pub fn new(r: f64, g: f64, b: f64, a: f64) -> Self {
        Self { r, g, b, a }
    }

    /// White.
    pub fn white() -> Self {
        Self::new(1.0, 1.0, 1.0, 1.0)
    }
}

/// Map a normalised value `t` in \[0,1\] to a jet-like colormap.
pub fn jet_colormap(t: f64) -> CrossSectionColor {
    let t = t.clamp(0.0, 1.0);
    // Blue at t=0, cyan, green, yellow, red at t=1
    let r = if t < 0.375 {
        0.0
    } else if t < 0.625 {
        (t - 0.375) / 0.25
    } else {
        1.0
    };
    let g = if t < 0.125 {
        0.0
    } else if t < 0.375 {
        (t - 0.125) / 0.25
    } else if t < 0.625 {
        1.0
    } else if t < 0.875 {
        1.0 - (t - 0.625) / 0.25
    } else {
        0.0
    };
    let b = if t < 0.375 {
        1.0
    } else if t < 0.625 {
        1.0 - (t - 0.375) / 0.25
    } else {
        0.0
    };
    CrossSectionColor::new(r.clamp(0.0, 1.0), g.clamp(0.0, 1.0), b.clamp(0.0, 1.0), 1.0)
}

/// Map a normalised value `t` in \[0,1\] to a cool-warm diverging colormap.
pub fn cool_warm_colormap(t: f64) -> CrossSectionColor {
    let t = t.clamp(0.0, 1.0);
    let r = t;
    let b = 1.0 - t;
    let g = 1.0 - 2.0 * (t - 0.5).abs();
    CrossSectionColor::new(r, g, b, 1.0)
}

/// Apply a color map to scalar values on a grid.
///
/// Returns a flat list of RGBA colors in row-major order.
pub fn colormap_grid(
    grid: &ScalarGrid2D,
    map_fn: fn(f64) -> CrossSectionColor,
) -> Vec<CrossSectionColor> {
    let vmin = grid.min_value();
    let vmax = grid.max_value();
    let range = vmax - vmin;
    grid.values
        .iter()
        .map(|&v| {
            let t = if range.abs() < 1e-15 {
                0.5
            } else {
                (v - vmin) / range
            };
            map_fn(t)
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Export cross-section as 2-D SVG path data
// ---------------------------------------------------------------------------

/// Export intersection segments as an SVG ``path` data string.
///
/// The segments are projected onto the cutting plane's 2-D coordinate system
/// and emitted as `M x y L x y` commands.
pub fn segments_to_svg_path(
    segments: &[IntersectionSegment],
    plane: &CuttingPlane,
    scale: f64,
) -> String {
    let mut path = String::new();
    for seg in segments {
        let s = plane.to_local_2d(seg.start);
        let e = plane.to_local_2d(seg.end);
        path.push_str(&format!(
            "M {:.6} {:.6} L {:.6} {:.6} ",
            s[0] * scale,
            s[1] * scale,
            e[0] * scale,
            e[1] * scale,
        ));
    }
    path.trim_end().to_string()
}

/// Export contour lines as SVG ``path` data.
pub fn contours_to_svg_path(contours: &[ContourLine], plane: &CuttingPlane, scale: f64) -> String {
    let mut path = String::new();
    for contour in contours {
        if contour.points.is_empty() {
            continue;
        }
        let uv0 = plane.to_local_2d(contour.points[0]);
        path.push_str(&format!("M {:.6} {:.6}", uv0[0] * scale, uv0[1] * scale));
        for pt in &contour.points[1..] {
            let uv = plane.to_local_2d(*pt);
            path.push_str(&format!(" L {:.6} {:.6}", uv[0] * scale, uv[1] * scale));
        }
        path.push(' ');
    }
    path.trim_end().to_string()
}

/// Generate a complete SVG document string from cross-section segments.
pub fn cross_section_to_svg(
    segments: &[IntersectionSegment],
    plane: &CuttingPlane,
    width: f64,
    height: f64,
    scale: f64,
) -> String {
    let path_data = segments_to_svg_path(segments, plane, scale);
    format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{:.6}\" height=\"{:.6}\" \
         viewBox=\"{:.6} {:.6} {:.6} {:.6}\">\n\
         <path d=\"{}\" fill=\"none\" stroke=\"black\" stroke-width=\"0.5\"/>\n\
         </svg>",
        width,
        height,
        -width / 2.0,
        -height / 2.0,
        width,
        height,
        path_data,
    )
}

// ---------------------------------------------------------------------------
// Additional utility: bounding box of segments
// ---------------------------------------------------------------------------

/// Compute the 2-D bounding box of intersection segments projected onto the
/// plane. Returns `(u_min, v_min, u_max, v_max)`.
pub fn segments_bounding_box_2d(
    segments: &[IntersectionSegment],
    plane: &CuttingPlane,
) -> Option<(f64, f64, f64, f64)> {
    if segments.is_empty() {
        return None;
    }
    let mut umin = f64::INFINITY;
    let mut vmin = f64::INFINITY;
    let mut umax = f64::NEG_INFINITY;
    let mut vmax = f64::NEG_INFINITY;
    for seg in segments {
        for pt in [seg.start, seg.end] {
            let uv = plane.to_local_2d(pt);
            umin = umin.min(uv[0]);
            umax = umax.max(uv[0]);
            vmin = vmin.min(uv[1]);
            vmax = vmax.max(uv[1]);
        }
    }
    Some((umin, vmin, umax, vmax))
}

/// Total length of all intersection segments.
pub fn total_segment_length(segments: &[IntersectionSegment]) -> f64 {
    segments.iter().map(|s| dist3(s.start, s.end)).sum()
}

/// Count unique vertices (within tolerance) across all segments.
pub fn count_unique_vertices(segments: &[IntersectionSegment], tol: f64) -> usize {
    let mut pts: Vec<[f64; 3]> = Vec::new();
    for seg in segments {
        for &p in &[seg.start, seg.end] {
            if !pts.iter().any(|q| dist3(*q, p) < tol) {
                pts.push(p);
            }
        }
    }
    pts.len()
}

// ---------------------------------------------------------------------------
// Grid-based operations
// ---------------------------------------------------------------------------

/// Compute the gradient of a scalar grid at node `(iu, iv)` using central
/// differences. Returns `[du, dv]`.
pub fn scalar_grid_gradient(grid: &ScalarGrid2D, iu: usize, iv: usize) -> [f64; 2] {
    let nu = grid.nodes_u();
    let nv = grid.nodes_v();
    let du = if iu == 0 {
        (grid.get(1, iv) - grid.get(0, iv)) / grid.du
    } else if iu == nu - 1 {
        (grid.get(nu - 1, iv) - grid.get(nu - 2, iv)) / grid.du
    } else {
        (grid.get(iu + 1, iv) - grid.get(iu - 1, iv)) / (2.0 * grid.du)
    };
    let dv = if iv == 0 {
        (grid.get(iu, 1) - grid.get(iu, 0)) / grid.dv
    } else if iv == nv - 1 {
        (grid.get(iu, nv - 1) - grid.get(iu, nv - 2)) / grid.dv
    } else {
        (grid.get(iu, iv + 1) - grid.get(iu, iv - 1)) / (2.0 * grid.dv)
    };
    [du, dv]
}

/// Compute the Laplacian of a scalar grid at node `(iu, iv)`.
pub fn scalar_grid_laplacian(grid: &ScalarGrid2D, iu: usize, iv: usize) -> f64 {
    let nu = grid.nodes_u();
    let nv = grid.nodes_v();
    let iu = iu.clamp(1, nu - 2);
    let iv = iv.clamp(1, nv - 2);
    let d2u = (grid.get(iu + 1, iv) - 2.0 * grid.get(iu, iv) + grid.get(iu - 1, iv))
        / (grid.du * grid.du);
    let d2v = (grid.get(iu, iv + 1) - 2.0 * grid.get(iu, iv) + grid.get(iu, iv - 1))
        / (grid.dv * grid.dv);
    d2u + d2v
}

// ---------------------------------------------------------------------------
// Mesh helper: build simple triangle meshes for testing
// ---------------------------------------------------------------------------

/// Create a list of triangles forming a unit cube centered at the origin.
pub fn unit_cube_triangles() -> Vec<Triangle3D> {
    let h = 0.5;
    let verts: [[f64; 3]; 8] = [
        [-h, -h, -h],
        [h, -h, -h],
        [h, h, -h],
        [-h, h, -h],
        [-h, -h, h],
        [h, -h, h],
        [h, h, h],
        [-h, h, h],
    ];
    let faces: [[usize; 3]; 12] = [
        // -Z
        [0, 2, 1],
        [0, 3, 2],
        // +Z
        [4, 5, 6],
        [4, 6, 7],
        // -Y
        [0, 1, 5],
        [0, 5, 4],
        // +Y
        [2, 3, 7],
        [2, 7, 6],
        // -X
        [0, 4, 7],
        [0, 7, 3],
        // +X
        [1, 2, 6],
        [1, 6, 5],
    ];
    faces
        .iter()
        .map(|f| Triangle3D::new(verts[f[0]], verts[f[1]], verts[f[2]]))
        .collect()
}

/// Create triangles forming a tetrahedron of given size.
pub fn tetrahedron_triangles(size: f64) -> Vec<Triangle3D> {
    let s = size;
    let a = [s, s, s];
    let b = [s, -s, -s];
    let c = [-s, s, -s];
    let d = [-s, -s, s];
    vec![
        Triangle3D::new(a, b, c),
        Triangle3D::new(a, b, d),
        Triangle3D::new(a, c, d),
        Triangle3D::new(b, c, d),
    ]
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ---- CuttingPlane tests ----

    #[test]
    fn test_cutting_plane_signed_distance() {
        let plane = CuttingPlane::xy(0.0);
        assert!((plane.signed_distance([0.0, 0.0, 1.0]) - 1.0).abs() < 1e-12);
        assert!((plane.signed_distance([0.0, 0.0, -1.0]) + 1.0).abs() < 1e-12);
        assert!(plane.signed_distance([0.0, 0.0, 0.0]).abs() < 1e-12);
    }

    #[test]
    fn test_cutting_plane_project() {
        let plane = CuttingPlane::xy(0.0);
        let proj = plane.project([3.0, 4.0, 5.0]);
        assert!((proj[0] - 3.0).abs() < 1e-12);
        assert!((proj[1] - 4.0).abs() < 1e-12);
        assert!(proj[2].abs() < 1e-12);
    }

    #[test]
    fn test_cutting_plane_tangent_basis_orthonormal() {
        let plane = CuttingPlane::new([0.0, 0.0, 0.0], [1.0, 1.0, 0.0]);
        let (u, v) = plane.tangent_basis();
        assert!((dot3(u, v)).abs() < 1e-12, "u and v should be orthogonal");
        assert!((norm3(u) - 1.0).abs() < 1e-12, "u should be unit");
        assert!((norm3(v) - 1.0).abs() < 1e-12, "v should be unit");
        assert!(
            dot3(u, plane.normal).abs() < 1e-12,
            "u should be perpendicular to normal"
        );
    }

    #[test]
    fn test_cutting_plane_local_2d_roundtrip() {
        let plane = CuttingPlane::new([1.0, 2.0, 3.0], [0.0, 0.0, 1.0]);
        let pt = [4.0, 5.0, 3.0]; // on the plane
        let uv = plane.to_local_2d(pt);
        let back = plane.from_local_2d(uv);
        for k in 0..3 {
            assert!(
                (back[k] - pt[k]).abs() < 1e-10,
                "roundtrip failed at component {k}: {:.6} vs {:.6}",
                back[k],
                pt[k]
            );
        }
    }

    #[test]
    fn test_cutting_plane_offset() {
        let plane = CuttingPlane::xy(0.0);
        let shifted = plane.offset(5.0);
        assert!((shifted.point[2] - 5.0).abs() < 1e-12);
    }

    // ---- Triangle tests ----

    #[test]
    fn test_triangle_area() {
        let tri = Triangle3D::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert!((tri.area() - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_triangle_centroid() {
        let tri = Triangle3D::new([0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [0.0, 3.0, 0.0]);
        let c = tri.centroid();
        assert!((c[0] - 1.0).abs() < 1e-12);
        assert!((c[1] - 1.0).abs() < 1e-12);
        assert!(c[2].abs() < 1e-12);
    }

    // ---- Triangle-plane clipping ----

    #[test]
    fn test_clip_triangle_intersects() {
        let tri = Triangle3D::new([-1.0, -1.0, -1.0], [1.0, -1.0, 1.0], [0.0, 1.0, 0.0]);
        let plane = CuttingPlane::xy(0.0);
        let seg = clip_triangle_plane(&tri, &plane);
        assert!(seg.is_some(), "triangle should intersect the plane");
    }

    #[test]
    fn test_clip_triangle_no_intersection() {
        let tri = Triangle3D::new([0.0, 0.0, 1.0], [1.0, 0.0, 2.0], [0.0, 1.0, 3.0]);
        let plane = CuttingPlane::xy(0.0);
        let seg = clip_triangle_plane(&tri, &plane);
        assert!(
            seg.is_none(),
            "triangle entirely above plane should not intersect"
        );
    }

    #[test]
    fn test_clip_mesh_cube_through_middle() {
        let cube = unit_cube_triangles();
        let plane = CuttingPlane::xy(0.0);
        let segments = clip_mesh_plane(&cube, &plane);
        assert!(
            !segments.is_empty(),
            "cutting a cube through the middle should produce segments"
        );
    }

    // ---- Contour extraction ----

    #[test]
    fn test_extract_contour_sphere_field() {
        let plane = CuttingPlane::xy(0.0);
        let mut grid = ScalarGrid2D::new(&plane, 2.0, 20);
        grid.sample_field(|p| p[0] * p[0] + p[1] * p[1]);
        let contours = extract_contour(&grid, 1.0);
        assert!(!contours.is_empty(), "should find contours for circle r=1");
        // All contour points should be approximately at distance 1 from origin
        for c in &contours {
            for pt in &c.points {
                let r = (pt[0] * pt[0] + pt[1] * pt[1]).sqrt();
                assert!(
                    (r - 1.0).abs() < 0.3,
                    "contour point should be near r=1, got {:.6}",
                    r
                );
            }
        }
    }

    #[test]
    fn test_extract_contours_multi_produces_levels() {
        let plane = CuttingPlane::xy(0.0);
        let mut grid = ScalarGrid2D::new(&plane, 2.0, 10);
        grid.sample_field(|p| p[0] + p[1]);
        let contours = extract_contours_multi(&grid, 5);
        assert!(
            contours.len() >= 3,
            "should produce multiple contour lines, got {}",
            contours.len()
        );
    }

    // ---- Scalar interpolation ----

    #[test]
    fn test_interpolate_scalar_grid_center() {
        let plane = CuttingPlane::xy(0.0);
        let mut grid = ScalarGrid2D::new(&plane, 1.0, 4);
        grid.sample_field(|p| p[0] + p[1]);
        let val = interpolate_scalar_grid(&grid, 0.0, 0.0);
        assert!(
            val.abs() < 0.5,
            "center value should be near 0, got {:.6}",
            val
        );
    }

    #[test]
    fn test_interpolate_at_world_point() {
        let plane = CuttingPlane::xy(0.0);
        let mut grid = ScalarGrid2D::new(&plane, 2.0, 10);
        grid.sample_field(|p| p[0]);
        let val = interpolate_at_world_point(&grid, [1.0, 0.0, 0.0]);
        assert!(
            (val - 1.0).abs() < 0.5,
            "should interpolate x=1 field near 1.0, got {:.6}",
            val
        );
    }

    // ---- Parallel slices ----

    #[test]
    fn test_generate_parallel_slices_count() {
        let config = ParallelSliceConfig {
            base_plane: CuttingPlane::xy(0.0),
            spacing: 0.1,
            count: 5,
        };
        let slices = generate_parallel_slices(&config);
        assert_eq!(slices.len(), 5);
    }

    #[test]
    fn test_parallel_slices_spacing() {
        let config = ParallelSliceConfig {
            base_plane: CuttingPlane::xy(0.0),
            spacing: 0.25,
            count: 4,
        };
        let slices = generate_parallel_slices(&config);
        for (i, slice) in slices.iter().enumerate() {
            let expected_z = i as f64 * 0.25;
            assert!(
                (slice.point[2] - expected_z).abs() < 1e-12,
                "slice {i} z={:.6}, expected {:.6}",
                slice.point[2],
                expected_z
            );
        }
    }

    #[test]
    fn test_clip_mesh_parallel_slices_cube() {
        let cube = unit_cube_triangles();
        let config = ParallelSliceConfig {
            base_plane: CuttingPlane::xy(-0.4),
            spacing: 0.2,
            count: 5,
        };
        let results = clip_mesh_parallel_slices(&cube, &config);
        assert_eq!(results.len(), 5);
        // At least the middle slices should hit the cube
        let non_empty = results.iter().filter(|r| !r.is_empty()).count();
        assert!(
            non_empty >= 2,
            "at least 2 slices should intersect the cube"
        );
    }

    // ---- Radial slices ----

    #[test]
    fn test_radial_slices_count() {
        let slices = generate_radial_slices([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 6);
        assert_eq!(slices.len(), 6);
    }

    #[test]
    fn test_radial_slices_all_pass_through_axis() {
        let slices = generate_radial_slices([1.0, 2.0, 0.0], [0.0, 0.0, 1.0], 8);
        for slice in &slices {
            let d = slice.signed_distance([1.0, 2.0, 0.0]);
            assert!(
                d.abs() < 1e-12,
                "radial slice should pass through axis point"
            );
        }
    }

    // ---- Cross-section area ----

    #[test]
    fn test_polygon_area_2d_unit_square() {
        let square = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let area = polygon_area_2d(&square);
        assert!((area - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_polygon_area_2d_triangle() {
        let tri = [[0.0, 0.0], [2.0, 0.0], [0.0, 2.0]];
        let area = polygon_area_2d(&tri);
        assert!((area - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_polygon_area_degenerate() {
        let line = [[0.0, 0.0], [1.0, 0.0]];
        assert!(polygon_area_2d(&line).abs() < 1e-12);
    }

    // ---- Section properties ----

    #[test]
    fn test_section_properties_unit_square() {
        let square = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let props = compute_section_properties(&square);
        assert!((props.area - 1.0).abs() < 1e-10, "area = {:.6}", props.area);
        assert!(
            (props.centroid[0] - 0.5).abs() < 1e-10,
            "cx = {:.6}",
            props.centroid[0]
        );
        assert!(
            (props.centroid[1] - 0.5).abs() < 1e-10,
            "cy = {:.6}",
            props.centroid[1]
        );
        // I_uu = I_vv = 1/12 for unit square about centroid
        assert!(
            (props.iuu - 1.0 / 12.0).abs() < 1e-6,
            "iuu = {:.6}",
            props.iuu
        );
        assert!(
            (props.ivv - 1.0 / 12.0).abs() < 1e-6,
            "ivv = {:.6}",
            props.ivv
        );
    }

    #[test]
    fn test_section_properties_degenerate() {
        let props = compute_section_properties(&[[0.0, 0.0], [1.0, 0.0]]);
        assert!(props.area < 1e-10);
    }

    // ---- Hatch lines ----

    #[test]
    fn test_hatch_lines_inside_square() {
        let square = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let config = HatchConfig {
            spacing: 0.2,
            angle: 0.0,
        };
        let lines = generate_hatch_lines(&square, &config);
        assert!(
            !lines.is_empty(),
            "should generate hatch lines inside the square"
        );
    }

    #[test]
    fn test_hatch_lines_empty_polygon() {
        let lines = generate_hatch_lines(&[], &HatchConfig::default());
        assert!(lines.is_empty());
    }

    // ---- Sweep animation ----

    #[test]
    fn test_sweep_animation_frame_count() {
        let cube = unit_cube_triangles();
        let base = CuttingPlane::xy(-0.5);
        let frames = generate_sweep_animation(&cube, &base, 1.0, 10);
        assert_eq!(frames.len(), 10);
    }

    #[test]
    fn test_sweep_animation_t_range() {
        let cube = unit_cube_triangles();
        let base = CuttingPlane::xy(0.0);
        let frames = generate_sweep_animation(&cube, &base, 1.0, 5);
        assert!((frames[0].t - 0.0).abs() < 1e-12);
        assert!((frames[4].t - 1.0).abs() < 1e-12);
    }

    // ---- Boolean operations ----

    #[test]
    fn test_cross_section_union_no_duplicates() {
        let a = vec![IntersectionSegment {
            start: [0.0, 0.0, 0.0],
            end: [1.0, 0.0, 0.0],
        }];
        let b = vec![
            IntersectionSegment {
                start: [0.0, 0.0, 0.0],
                end: [1.0, 0.0, 0.0],
            },
            IntersectionSegment {
                start: [2.0, 0.0, 0.0],
                end: [3.0, 0.0, 0.0],
            },
        ];
        let result = cross_section_union(&a, &b);
        assert_eq!(result.segments.len(), 2, "duplicate should be removed");
    }

    #[test]
    fn test_cross_section_intersection_empty_b() {
        let a = vec![IntersectionSegment {
            start: [0.0, 0.0, 0.0],
            end: [1.0, 0.0, 0.0],
        }];
        let plane = CuttingPlane::xy(0.0);
        let result = cross_section_intersection(&a, &[], &plane);
        assert!(result.segments.is_empty());
    }

    // ---- Stress field ----

    #[test]
    fn test_von_mises_uniaxial() {
        let stress = [100.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let vm = von_mises_from_tensor(&stress);
        assert!(
            (vm - 100.0).abs() < 1e-6,
            "von Mises of uniaxial 100 should be 100, got {:.6}",
            vm
        );
    }

    #[test]
    fn test_von_mises_hydrostatic() {
        let stress = [50.0, 50.0, 50.0, 0.0, 0.0, 0.0];
        let vm = von_mises_from_tensor(&stress);
        assert!(
            vm < 1e-6,
            "hydrostatic stress should give VM=0, got {:.6}",
            vm
        );
    }

    #[test]
    fn test_sample_stress_field_count() {
        let plane = CuttingPlane::xy(0.0);
        let samples = sample_stress_field(&plane, 1.0, 4, |_| [100.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
        // (4+1) * (4+1) = 25 nodes
        assert_eq!(samples.len(), 25);
    }

    // ---- Color mapping ----

    #[test]
    fn test_jet_colormap_endpoints() {
        let blue_end = jet_colormap(0.0);
        assert!(blue_end.b > 0.5, "jet(0) should be blueish");
        let red_end = jet_colormap(1.0);
        assert!(red_end.r > 0.5, "jet(1) should be reddish");
    }

    #[test]
    fn test_cool_warm_colormap_midpoint() {
        let mid = cool_warm_colormap(0.5);
        assert!(mid.g > 0.8, "cool-warm midpoint should have high green");
    }

    #[test]
    fn test_colormap_grid_length() {
        let plane = CuttingPlane::xy(0.0);
        let mut grid = ScalarGrid2D::new(&plane, 1.0, 3);
        grid.sample_field(|p| p[0]);
        let colors = colormap_grid(&grid, jet_colormap);
        assert_eq!(colors.len(), grid.values.len());
    }

    // ---- SVG export ----

    #[test]
    fn test_segments_to_svg_path_basic() {
        let plane = CuttingPlane::xy(0.0);
        let segments = vec![IntersectionSegment {
            start: [0.0, 0.0, 0.0],
            end: [1.0, 1.0, 0.0],
        }];
        let path = segments_to_svg_path(&segments, &plane, 100.0);
        assert!(path.contains('M'), "SVG path should contain M command");
        assert!(path.contains('L'), "SVG path should contain L command");
    }

    #[test]
    fn test_cross_section_to_svg_valid() {
        let plane = CuttingPlane::xy(0.0);
        let segments = vec![IntersectionSegment {
            start: [0.0, 0.0, 0.0],
            end: [1.0, 0.0, 0.0],
        }];
        let svg = cross_section_to_svg(&segments, &plane, 200.0, 200.0, 100.0);
        assert!(svg.contains("<svg"));
        assert!(svg.contains("<path"));
    }

    #[test]
    fn test_contours_to_svg_path() {
        let plane = CuttingPlane::xy(0.0);
        let contour = ContourLine {
            points: vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]],
            value: 1.0,
        };
        let path = contours_to_svg_path(&[contour], &plane, 1.0);
        assert!(path.contains('M'));
        assert!(path.contains('L'));
    }

    // ---- Utility functions ----

    #[test]
    fn test_total_segment_length() {
        let segments = vec![
            IntersectionSegment {
                start: [0.0, 0.0, 0.0],
                end: [1.0, 0.0, 0.0],
            },
            IntersectionSegment {
                start: [0.0, 0.0, 0.0],
                end: [0.0, 3.0, 4.0],
            },
        ];
        let total = total_segment_length(&segments);
        assert!((total - 6.0).abs() < 1e-12, "1 + 5 = 6, got {:.6}", total);
    }

    #[test]
    fn test_segments_bounding_box() {
        let plane = CuttingPlane::xy(0.0);
        let segments = vec![IntersectionSegment {
            start: [-1.0, -2.0, 0.0],
            end: [3.0, 4.0, 0.0],
        }];
        let (umin, vmin, umax, vmax) = segments_bounding_box_2d(&segments, &plane).unwrap();
        // The bounding box in local 2D should span the projected range
        let du = umax - umin;
        let dv = vmax - vmin;
        assert!(du > 0.0, "u range should be positive, got {:.6}", du);
        assert!(dv > 0.0, "v range should be positive, got {:.6}", dv);
        // The total extent should match the 2D distance between the points
        let total = (du * du + dv * dv).sqrt();
        let expected = dist3([-1.0, -2.0, 0.0], [3.0, 4.0, 0.0]);
        assert!(
            (total - expected).abs() < 0.1,
            "bbox diagonal {:.6} should match distance {:.6}",
            total,
            expected
        );
    }

    #[test]
    fn test_count_unique_vertices() {
        let segments = vec![
            IntersectionSegment {
                start: [0.0, 0.0, 0.0],
                end: [1.0, 0.0, 0.0],
            },
            IntersectionSegment {
                start: [1.0, 0.0, 0.0],
                end: [2.0, 0.0, 0.0],
            },
        ];
        let count = count_unique_vertices(&segments, 1e-10);
        assert_eq!(count, 3, "shared vertex should be counted once");
    }

    // ---- Grid gradient / Laplacian ----

    #[test]
    fn test_scalar_grid_gradient_linear() {
        let plane = CuttingPlane::xy(0.0);
        let mut grid = ScalarGrid2D::new(&plane, 1.0, 10);
        // Sample a field that's linear in the grid's u-direction
        let (u_axis, _v_axis) = plane.tangent_basis();
        grid.sample_field(|p| dot3(p, u_axis));
        let grad = scalar_grid_gradient(&grid, 5, 5);
        assert!(
            (grad[0] - 1.0).abs() < 0.2,
            "df/du should be ~1, got {:.6}",
            grad[0]
        );
        assert!(
            grad[1].abs() < 0.2,
            "df/dv should be ~0, got {:.6}",
            grad[1]
        );
    }

    #[test]
    fn test_scalar_grid_laplacian_quadratic() {
        let plane = CuttingPlane::xy(0.0);
        let mut grid = ScalarGrid2D::new(&plane, 2.0, 20);
        // f(u,v) = u^2 + v^2 => Laplacian = 2 + 2 = 4
        grid.sample_field(|p| p[0] * p[0] + p[1] * p[1]);
        let lap = scalar_grid_laplacian(&grid, 10, 10);
        assert!(
            (lap - 4.0).abs() < 0.5,
            "Laplacian of r^2 should be ~4, got {:.6}",
            lap
        );
    }

    // ---- Unit cube / tetrahedron helpers ----

    #[test]
    fn test_unit_cube_triangles_count() {
        let tris = unit_cube_triangles();
        assert_eq!(tris.len(), 12, "cube has 12 triangles (2 per face)");
    }

    #[test]
    fn test_tetrahedron_triangles_count() {
        let tris = tetrahedron_triangles(1.0);
        assert_eq!(tris.len(), 4, "tetrahedron has 4 faces");
    }

    // ---- Cylindrical cut ----

    #[test]
    fn test_cylindrical_cut_radial_distance() {
        let cyl = CylindricalCut::new([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 1.0);
        let d = cyl.radial_distance([2.0, 0.0, 5.0]);
        assert!(
            (d - 1.0).abs() < 1e-12,
            "distance should be 1, got {:.6}",
            d
        );
        let d_on = cyl.radial_distance([1.0, 0.0, 0.0]);
        assert!(d_on.abs() < 1e-12, "on surface, distance should be 0");
    }

    #[test]
    fn test_scalar_grid_min_max() {
        let plane = CuttingPlane::xy(0.0);
        let mut grid = ScalarGrid2D::new(&plane, 1.0, 4);
        grid.sample_field(|p| p[0]);
        assert!(grid.min_value() < 0.0);
        assert!(grid.max_value() > 0.0);
    }
}
