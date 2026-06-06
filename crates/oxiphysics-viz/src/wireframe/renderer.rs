// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Rich wireframe mesh, builder, frustum, hidden-line removal, and renderer.

use super::mesh_ops::{
    AaWireframeParams, WireSegmentGpu, WireframeBatch, build_aa_feature_edge_batch,
    detect_feature_edges, detect_silhouette_edges, extract_mesh_edges,
};

// ─── WireframeEdge (enhanced, with color and line width) ─────────────────────

/// An edge with per-edge color and line width, referencing vertex indices.
#[derive(Debug, Clone, Copy)]
pub struct WireframeEdge {
    /// Index of the first vertex.
    pub v0: usize,
    /// Index of the second vertex.
    pub v1: usize,
    /// RGB color of this edge.
    pub color: [f32; 3],
    /// Screen-space line width in pixels.
    pub line_width: f32,
}

impl WireframeEdge {
    /// Create a new edge with the given indices, color, and width.
    pub fn new(v0: usize, v1: usize, color: [f32; 3], line_width: f32) -> Self {
        Self {
            v0,
            v1,
            color,
            line_width,
        }
    }

    /// Create a default white 1-pixel-wide edge.
    pub fn default_white(v0: usize, v1: usize) -> Self {
        Self {
            v0,
            v1,
            color: [1.0, 1.0, 1.0],
            line_width: 1.0,
        }
    }

    /// Returns the canonical (sorted) vertex pair so (a,b) always has a <= b.
    pub fn canonical_key(&self) -> (usize, usize) {
        (self.v0.min(self.v1), self.v0.max(self.v1))
    }
}

// ─── RichWireframeMesh ────────────────────────────────────────────────────────

/// A full wireframe mesh that stores vertices, per-vertex normals, and per-edge
/// color/width information.
#[derive(Debug, Clone)]
pub struct RichWireframeMesh {
    /// Vertex positions.
    pub vertices: Vec<[f64; 3]>,
    /// Per-vertex outward normals (unit length).
    pub normals: Vec<[f64; 3]>,
    /// Edges, each referencing vertex indices.
    pub edges: Vec<WireframeEdge>,
}

impl RichWireframeMesh {
    /// Create an empty mesh.
    pub fn empty() -> Self {
        Self {
            vertices: Vec::new(),
            normals: Vec::new(),
            edges: Vec::new(),
        }
    }

    /// Number of vertices.
    pub fn vertex_count(&self) -> usize {
        self.vertices.len()
    }

    /// Number of edges.
    pub fn edge_count(&self) -> usize {
        self.edges.len()
    }

    /// Merge another `RichWireframeMesh` into this one in-place.
    pub fn merge_from(&mut self, other: &Self) {
        let offset = self.vertices.len();
        self.vertices.extend_from_slice(&other.vertices);
        self.normals.extend_from_slice(&other.normals);
        for e in &other.edges {
            self.edges.push(WireframeEdge {
                v0: e.v0 + offset,
                v1: e.v1 + offset,
                color: e.color,
                line_width: e.line_width,
            });
        }
    }
}

// ─── WireframeBuilder ─────────────────────────────────────────────────────────

/// Builder pattern for constructing a `RichWireframeMesh` from triangle meshes or
/// geometric primitives.  Edges are automatically deduplicated on `build()`.
#[derive(Debug, Clone)]
pub struct WireframeBuilder {
    vertices: Vec<[f64; 3]>,
    normals: Vec<[f64; 3]>,
    // Pending edges before dedup
    raw_edges: Vec<WireframeEdge>,
    default_color: [f32; 3],
    default_width: f32,
}

impl WireframeBuilder {
    /// Create a new builder with a default color and line width.
    pub fn new(default_color: [f32; 3], default_width: f32) -> Self {
        Self {
            vertices: Vec::new(),
            normals: Vec::new(),
            raw_edges: Vec::new(),
            default_color,
            default_width,
        }
    }

    /// Add a vertex with an explicit normal; returns the vertex index.
    pub fn add_vertex(&mut self, pos: [f64; 3], normal: [f64; 3]) -> usize {
        let idx = self.vertices.len();
        self.vertices.push(pos);
        self.normals.push(normal);
        idx
    }

    /// Add a raw (possibly duplicate) edge between two vertex indices.
    pub fn add_edge(&mut self, v0: usize, v1: usize) {
        self.raw_edges.push(WireframeEdge::new(
            v0,
            v1,
            self.default_color,
            self.default_width,
        ));
    }

    /// Add a raw edge with explicit color and width.
    pub fn add_edge_colored(&mut self, v0: usize, v1: usize, color: [f32; 3], width: f32) {
        self.raw_edges
            .push(WireframeEdge::new(v0, v1, color, width));
    }

    /// Feed a triangle mesh: deduplicated edges are added for each triangle side.
    /// `positions` and `indices` are flat vertex + index buffers (indices in triples).
    /// `normals` must match `positions` in length.
    pub fn add_triangle_mesh(
        &mut self,
        positions: &[[f64; 3]],
        tri_indices: &[usize],
        normals: &[[f64; 3]],
    ) {
        let base = self.vertices.len();
        for (&p, &n) in positions.iter().zip(normals.iter()) {
            self.vertices.push(p);
            self.normals.push(n);
        }
        for tri in tri_indices.chunks_exact(3) {
            let (a, b, c) = (tri[0] + base, tri[1] + base, tri[2] + base);
            self.add_edge(a, b);
            self.add_edge(b, c);
            self.add_edge(c, a);
        }
    }

    /// Build the final `RichWireframeMesh`, deduplicating edges.
    ///
    /// When duplicate edges exist with different colors/widths, the first
    /// encountered is kept (based on sorted vertex key).
    pub fn build(self) -> RichWireframeMesh {
        use std::collections::HashMap;
        let mut seen: HashMap<(usize, usize), WireframeEdge> = HashMap::new();
        for e in self.raw_edges {
            let key = e.canonical_key();
            seen.entry(key).or_insert(e);
        }
        let mut edges: Vec<WireframeEdge> = seen.into_values().collect();
        edges.sort_by_key(|e| e.canonical_key());
        RichWireframeMesh {
            vertices: self.vertices,
            normals: self.normals,
            edges,
        }
    }
}

// ─── AABB wireframe (RichWireframeMesh) ──────────────────────────────────────

/// Generate an axis-aligned bounding box as a `RichWireframeMesh`.
///
/// `min` and `max` are the corner coordinates.
pub fn aabb_wireframe_rich(
    min: [f64; 3],
    max: [f64; 3],
    color: [f32; 3],
    line_width: f32,
) -> RichWireframeMesh {
    let corners: [[f64; 3]; 8] = [
        [min[0], min[1], min[2]],
        [max[0], min[1], min[2]],
        [max[0], max[1], min[2]],
        [min[0], max[1], min[2]],
        [min[0], min[1], max[2]],
        [max[0], min[1], max[2]],
        [max[0], max[1], max[2]],
        [min[0], max[1], max[2]],
    ];

    // Outward normals are approximated as the direction from center to vertex
    let cx = (min[0] + max[0]) * 0.5;
    let cy = (min[1] + max[1]) * 0.5;
    let cz = (min[2] + max[2]) * 0.5;

    let mut builder = WireframeBuilder::new(color, line_width);
    for &c in &corners {
        let dx = c[0] - cx;
        let dy = c[1] - cy;
        let dz = c[2] - cz;
        let mag = (dx * dx + dy * dy + dz * dz).sqrt().max(1e-30);
        builder.add_vertex(c, [dx / mag, dy / mag, dz / mag]);
    }

    let edge_pairs: [(usize, usize); 12] = [
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 0),
        (4, 5),
        (5, 6),
        (6, 7),
        (7, 4),
        (0, 4),
        (1, 5),
        (2, 6),
        (3, 7),
    ];
    for &(a, b) in &edge_pairs {
        builder.add_edge(a, b);
    }
    builder.build()
}

// ─── Frustum wireframe ────────────────────────────────────────────────────────

/// Generate a perspective frustum wireframe as a `RichWireframeMesh`.
///
/// The frustum is defined by a near plane (at distance `z_near`) and a far plane
/// (at distance `z_far`), with half-angles `half_fov_x` and `half_fov_y` (radians).
/// The camera is assumed to look along -Z by convention.
pub fn frustum_wireframe(
    z_near: f64,
    z_far: f64,
    half_fov_x: f64,
    half_fov_y: f64,
    color: [f32; 3],
    line_width: f32,
) -> RichWireframeMesh {
    let xn = z_near * half_fov_x.tan();
    let yn = z_near * half_fov_y.tan();
    let xf = z_far * half_fov_x.tan();
    let yf = z_far * half_fov_y.tan();

    // Near plane corners (z = -z_near, camera looking down -Z)
    let near_corners: [[f64; 3]; 4] = [
        [-xn, -yn, -z_near],
        [xn, -yn, -z_near],
        [xn, yn, -z_near],
        [-xn, yn, -z_near],
    ];
    // Far plane corners
    let far_corners: [[f64; 3]; 4] = [
        [-xf, -yf, -z_far],
        [xf, -yf, -z_far],
        [xf, yf, -z_far],
        [-xf, yf, -z_far],
    ];

    let mut builder = WireframeBuilder::new(color, line_width);

    // Add near corners (indices 0-3), normals pointing toward origin
    for &c in &near_corners {
        let mag = (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]).sqrt().max(1e-30);
        builder.add_vertex(c, [-c[0] / mag, -c[1] / mag, -c[2] / mag]);
    }
    // Add far corners (indices 4-7)
    for &c in &far_corners {
        let mag = (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]).sqrt().max(1e-30);
        builder.add_vertex(c, [-c[0] / mag, -c[1] / mag, -c[2] / mag]);
    }

    // Near rectangle
    builder.add_edge(0, 1);
    builder.add_edge(1, 2);
    builder.add_edge(2, 3);
    builder.add_edge(3, 0);
    // Far rectangle
    builder.add_edge(4, 5);
    builder.add_edge(5, 6);
    builder.add_edge(6, 7);
    builder.add_edge(7, 4);
    // Frustum edges (near to far)
    builder.add_edge(0, 4);
    builder.add_edge(1, 5);
    builder.add_edge(2, 6);
    builder.add_edge(3, 7);

    builder.build()
}

/// Generate a frustum wireframe from a view/projection matrix representation,
/// defined by `fov_y_deg` (full vertical FOV), `aspect` (width/height),
/// `z_near`, and `z_far`.
pub fn frustum_wireframe_from_fov(
    fov_y_deg: f64,
    aspect: f64,
    z_near: f64,
    z_far: f64,
    color: [f32; 3],
    line_width: f32,
) -> RichWireframeMesh {
    let half_fov_y = (fov_y_deg * 0.5).to_radians();
    let half_fov_x = (half_fov_y.tan() * aspect).atan();
    frustum_wireframe(z_near, z_far, half_fov_x, half_fov_y, color, line_width)
}

// ─── Hidden Line Removal (Depth-based Occlusion) ─────────────────────────────

/// Classify wireframe edges as visible or hidden via a depth-based occlusion check.
///
/// For each edge, we sample its midpoint and test whether any triangle in the mesh
/// occludes it from `view_origin`.  Edges whose midpoint is occluded by another
/// triangle (with a depth offset tolerance) are classified as hidden.
///
/// Returns `(visible_edges, hidden_edges)` as index lists into `edges`.
pub fn classify_hidden_lines(
    vertices: &[[f64; 3]],
    edges: &[WireframeEdge],
    triangles: &[usize],
    view_origin: [f64; 3],
    depth_bias: f64,
) -> (Vec<usize>, Vec<usize>) {
    let num_tris = triangles.len() / 3;

    // Möller–Trumbore ray–triangle intersection
    let ray_tri_intersect =
        |orig: [f64; 3], dir: [f64; 3], a: [f64; 3], b: [f64; 3], c: [f64; 3]| -> Option<f64> {
            let e1 = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
            let e2 = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
            let h = [
                dir[1] * e2[2] - dir[2] * e2[1],
                dir[2] * e2[0] - dir[0] * e2[2],
                dir[0] * e2[1] - dir[1] * e2[0],
            ];
            let a_det = e1[0] * h[0] + e1[1] * h[1] + e1[2] * h[2];
            if a_det.abs() < 1e-12 {
                return None;
            }
            let inv = 1.0 / a_det;
            let s = [orig[0] - a[0], orig[1] - a[1], orig[2] - a[2]];
            let u = inv * (s[0] * h[0] + s[1] * h[1] + s[2] * h[2]);
            if !(0.0..=1.0).contains(&u) {
                return None;
            }
            let q = [
                s[1] * e1[2] - s[2] * e1[1],
                s[2] * e1[0] - s[0] * e1[2],
                s[0] * e1[1] - s[1] * e1[0],
            ];
            let v = inv * (dir[0] * q[0] + dir[1] * q[1] + dir[2] * q[2]);
            if v < 0.0 || u + v > 1.0 {
                return None;
            }
            let t = inv * (e2[0] * q[0] + e2[1] * q[1] + e2[2] * q[2]);
            if t > 1e-12 { Some(t) } else { None }
        };

    let mut visible = Vec::new();
    let mut hidden = Vec::new();

    for (edge_idx, edge) in edges.iter().enumerate() {
        if edge.v0 >= vertices.len() || edge.v1 >= vertices.len() {
            visible.push(edge_idx);
            continue;
        }
        let p0 = vertices[edge.v0];
        let p1 = vertices[edge.v1];
        let mid = [
            (p0[0] + p1[0]) * 0.5,
            (p0[1] + p1[1]) * 0.5,
            (p0[2] + p1[2]) * 0.5,
        ];
        // Ray from view_origin to midpoint
        let dx = mid[0] - view_origin[0];
        let dy = mid[1] - view_origin[1];
        let dz = mid[2] - view_origin[2];
        let dist = (dx * dx + dy * dy + dz * dz).sqrt();
        if dist < 1e-12 {
            visible.push(edge_idx);
            continue;
        }
        let dir = [dx / dist, dy / dist, dz / dist];

        let mut occluded = false;
        for t in 0..num_tris {
            let (ia, ib, ic) = (triangles[3 * t], triangles[3 * t + 1], triangles[3 * t + 2]);
            if ia >= vertices.len() || ib >= vertices.len() || ic >= vertices.len() {
                continue;
            }
            // Skip triangles that share vertices with this edge
            if ia == edge.v0
                || ia == edge.v1
                || ib == edge.v0
                || ib == edge.v1
                || ic == edge.v0
                || ic == edge.v1
            {
                continue;
            }
            if let Some(hit_t) =
                ray_tri_intersect(view_origin, dir, vertices[ia], vertices[ib], vertices[ic])
            {
                // Edge midpoint is at distance `dist`; if triangle is closer (minus bias)
                if hit_t < dist - depth_bias {
                    occluded = true;
                    break;
                }
            }
        }

        if occluded {
            hidden.push(edge_idx);
        } else {
            visible.push(edge_idx);
        }
    }

    (visible, hidden)
}

// ─── Overlay Generation ───────────────────────────────────────────────────────

/// Generate a combined wireframe overlay: all edges as a base, feature edges highlighted,
/// silhouette edges thick. Returns three `WireframeBatch` objects in order:
/// `(all_edges, feature_edges, silhouette_edges)`.
pub fn generate_wireframe_overlay(
    positions: &[[f64; 3]],
    indices: &[usize],
    view_origin: [f64; 3],
    feature_threshold_deg: f64,
    base_params: &AaWireframeParams,
    feature_params: &AaWireframeParams,
    silhouette_params: &AaWireframeParams,
) -> (WireframeBatch, WireframeBatch, WireframeBatch) {
    let all_edges = extract_mesh_edges(positions, indices);
    let feature_edges = detect_feature_edges(positions, indices, feature_threshold_deg);
    let silhouette_edges = detect_silhouette_edges(positions, indices, view_origin);

    let base_batch = build_aa_feature_edge_batch(positions, &all_edges, base_params);
    let feature_batch = build_aa_feature_edge_batch(positions, &feature_edges, feature_params);
    let silhouette_batch =
        build_aa_feature_edge_batch(positions, &silhouette_edges, silhouette_params);

    (base_batch, feature_batch, silhouette_batch)
}

#[cfg(test)]
mod rich_wireframe_tests {
    use super::*;

    use std::f64::consts::FRAC_PI_4;

    // ── WireframeEdge ─────────────────────────────────────────────────────────

    #[test]
    fn test_wireframe_edge_canonical_key_sorted() {
        let e = WireframeEdge::new(7, 3, [1.0, 0.0, 0.0], 1.0);
        let key = e.canonical_key();
        assert!(key.0 <= key.1, "canonical key must have a <= b");
        assert_eq!(key, (3, 7));
    }

    #[test]
    fn test_wireframe_edge_default_white() {
        let e = WireframeEdge::default_white(2, 9);
        assert_eq!(e.color, [1.0, 1.0, 1.0]);
        assert!((e.line_width - 1.0).abs() < 1e-9);
        assert_eq!(e.v0, 2);
        assert_eq!(e.v1, 9);
    }

    // ── WireframeBuilder / RichWireframeMesh ─────────────────────────────────

    #[test]
    fn test_builder_deduplicates_edges() {
        let mut builder = WireframeBuilder::new([1.0, 1.0, 1.0], 1.0);
        builder.add_vertex([0.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        builder.add_vertex([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        builder.add_vertex([0.5, 1.0, 0.0], [0.0, 1.0, 0.0]);
        // Add the same edge twice in both directions
        builder.add_edge(0, 1);
        builder.add_edge(1, 0);
        builder.add_edge(1, 2);
        builder.add_edge(0, 2);
        let mesh = builder.build();
        // Should have exactly 3 unique edges (0-1, 1-2, 0-2)
        assert_eq!(mesh.edge_count(), 3);
    }

    #[test]
    fn test_builder_add_triangle_mesh_dedup() {
        // Two triangles sharing an edge
        let positions = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ];
        let normals = vec![[0.0, 0.0, 1.0]; 4];
        let indices = vec![0usize, 1, 2, 0, 2, 3];
        let mut builder = WireframeBuilder::new([1.0, 1.0, 1.0], 1.0);
        builder.add_triangle_mesh(&positions, &indices, &normals);
        let mesh = builder.build();
        // Quad split into 2 triangles = 5 unique edges (perimeter 4 + 1 diagonal)
        assert_eq!(mesh.edge_count(), 5);
        assert_eq!(mesh.vertex_count(), 4);
    }

    #[test]
    fn test_rich_wireframe_mesh_merge() {
        let a = aabb_wireframe_rich([0.0; 3], [1.0; 3], [1.0, 0.0, 0.0], 1.0);
        let b = aabb_wireframe_rich([2.0; 3], [3.0; 3], [0.0, 1.0, 0.0], 1.0);
        let mut merged = a.clone();
        merged.merge_from(&b);
        assert_eq!(merged.vertex_count(), a.vertex_count() + b.vertex_count());
        assert_eq!(merged.edge_count(), a.edge_count() + b.edge_count());
        // All edge indices must be in bounds
        for e in &merged.edges {
            assert!(e.v0 < merged.vertex_count());
            assert!(e.v1 < merged.vertex_count());
        }
    }

    // ── AABB wireframe (rich) ─────────────────────────────────────────────────

    #[test]
    fn test_aabb_wireframe_rich_12_edges() {
        let mesh = aabb_wireframe_rich([0.0; 3], [1.0; 3], [1.0, 1.0, 1.0], 1.0);
        assert_eq!(mesh.vertex_count(), 8);
        assert_eq!(mesh.edge_count(), 12);
        assert_eq!(mesh.normals.len(), 8);
    }

    #[test]
    fn test_aabb_wireframe_normals_nonzero() {
        let mesh = aabb_wireframe_rich([-1.0; 3], [1.0; 3], [1.0, 1.0, 1.0], 1.0);
        for (i, n) in mesh.normals.iter().enumerate() {
            let mag = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
            assert!(mag > 0.9, "normal {i} should be roughly unit, mag={mag}");
        }
    }

    #[test]
    fn test_aabb_wireframe_rich_edge_bounds() {
        let mesh = aabb_wireframe_rich([0.0; 3], [2.0; 3], [0.5, 0.5, 0.5], 1.5);
        for e in &mesh.edges {
            assert!(e.v0 < mesh.vertex_count());
            assert!(e.v1 < mesh.vertex_count());
        }
    }

    // ── Frustum wireframe ─────────────────────────────────────────────────────

    #[test]
    fn test_frustum_wireframe_12_edges() {
        use std::f64::consts::FRAC_PI_4;
        let mesh = frustum_wireframe(
            0.1,
            100.0,
            FRAC_PI_4,
            FRAC_PI_4 * 0.75,
            [1.0, 1.0, 0.0],
            1.0,
        );
        assert_eq!(mesh.vertex_count(), 8);
        // 4 near + 4 far + 4 connecting = 12 edges
        assert_eq!(mesh.edge_count(), 12);
    }

    #[test]
    fn test_frustum_wireframe_from_fov() {
        let mesh = frustum_wireframe_from_fov(60.0, 16.0 / 9.0, 0.1, 500.0, [0.0, 1.0, 1.0], 1.0);
        assert_eq!(mesh.vertex_count(), 8);
        assert_eq!(mesh.edge_count(), 12);
        // Far plane verts should be further from origin than near plane verts
        let near_z = mesh.vertices[0][2].abs();
        let far_z = mesh.vertices[4][2].abs();
        assert!(
            far_z > near_z,
            "far corners should be farther: far_z={far_z}, near_z={near_z}"
        );
    }

    #[test]
    fn test_frustum_normals_inward() {
        let mesh = frustum_wireframe(1.0, 10.0, FRAC_PI_4, FRAC_PI_4, [1.0, 0.0, 0.0], 1.0);
        // At least normals should be non-zero
        for (i, n) in mesh.normals.iter().enumerate() {
            let mag = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
            assert!(mag > 0.5, "normal {i} should be roughly unit, mag={mag}");
        }
    }

    // ── Hidden line removal ───────────────────────────────────────────────────

    #[test]
    fn test_hidden_line_single_triangle_no_occlusion() {
        let verts = vec![[0.0, 0.0, 0.0f64], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]];
        let edges = vec![
            WireframeEdge::default_white(0, 1),
            WireframeEdge::default_white(1, 2),
            WireframeEdge::default_white(2, 0),
        ];
        let tris = vec![0usize, 1, 2];
        let (vis, hid) = classify_hidden_lines(&verts, &edges, &tris, [0.5, 0.5, 10.0], 0.01);
        assert_eq!(vis.len(), 3, "all edges should be visible from front");
        assert_eq!(hid.len(), 0);
    }

    #[test]
    fn test_hidden_line_backface_triangle_hides_edge() {
        // Front quad at z=0, back quad at z=-2
        // An edge in the back quad behind the front quad should be hidden
        let verts = vec![
            // Front quad (z=0)
            [-5.0, -5.0, 0.0f64],
            [5.0, -5.0, 0.0],
            [5.0, 5.0, 0.0],
            [-5.0, 5.0, 0.0],
            // Back edge midpoint (z=-2)
            [-0.1, -0.1, -2.0],
            [0.1, 0.1, -2.0],
        ];
        // Front two triangles covering a wide area
        let tris = vec![0usize, 1, 2, 0, 2, 3];
        let back_edges = vec![WireframeEdge::default_white(4, 5)];
        let view = [0.0, 0.0, 10.0];
        let (vis, hid) = classify_hidden_lines(&verts, &back_edges, &tris, view, 0.01);
        assert_eq!(
            hid.len(),
            1,
            "back edge should be hidden by front quad, vis={vis:?}"
        );
        assert_eq!(vis.len(), 0);
    }
}

// ─── WireframeRenderer ────────────────────────────────────────────────────────

/// A stateful wireframe renderer that accumulates line segments into a buffer,
/// providing high-level draw methods for common physics debug shapes.
#[derive(Debug, Clone)]
pub struct WireframeRenderer {
    /// Accumulated GPU-ready line segments.
    pub buffer: Vec<WireSegmentGpu>,
    /// Default line width in pixels.
    pub default_width: f64,
    /// Default color (RGB, 0.0–1.0).
    pub default_color: [f64; 3],
}

impl WireframeRenderer {
    /// Create a new renderer with the given default color and line width.
    pub fn new(default_color: [f64; 3], default_width: f64) -> Self {
        Self {
            buffer: Vec::new(),
            default_width,
            default_color,
        }
    }

    /// Clear all accumulated segments.
    pub fn clear(&mut self) {
        self.buffer.clear();
    }

    /// Number of segments in the buffer.
    pub fn segment_count(&self) -> usize {
        self.buffer.len()
    }

    /// Push a single line segment with an explicit color and width.
    pub fn push_segment(&mut self, start: [f64; 3], end: [f64; 3], color: [f64; 3], width: f64) {
        self.buffer
            .push(WireSegmentGpu::new(start, end, color, width));
    }

    /// Push a segment using the renderer's default color and width.
    pub fn push(&mut self, start: [f64; 3], end: [f64; 3]) {
        self.push_segment(start, end, self.default_color, self.default_width);
    }

    /// Draw an axis-aligned box given center and half-extents.
    pub fn draw_box(&mut self, center: [f64; 3], half: [f64; 3], color: [f64; 3]) {
        let [cx, cy, cz] = center;
        let [hx, hy, hz] = half;
        let corners = [
            [cx - hx, cy - hy, cz - hz],
            [cx + hx, cy - hy, cz - hz],
            [cx + hx, cy + hy, cz - hz],
            [cx - hx, cy + hy, cz - hz],
            [cx - hx, cy - hy, cz + hz],
            [cx + hx, cy - hy, cz + hz],
            [cx + hx, cy + hy, cz + hz],
            [cx - hx, cy + hy, cz + hz],
        ];
        let edges: [(usize, usize); 12] = [
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 0),
            (4, 5),
            (5, 6),
            (6, 7),
            (7, 4),
            (0, 4),
            (1, 5),
            (2, 6),
            (3, 7),
        ];
        for (a, b) in edges {
            self.push_segment(corners[a], corners[b], color, self.default_width);
        }
    }

    /// Draw a sphere as three great circles (XY, XZ, YZ planes).
    pub fn draw_sphere(&mut self, center: [f64; 3], radius: f64, segments: usize, color: [f64; 3]) {
        let seg = segments.max(4);
        let [cx, cy, cz] = center;
        // XY circle
        for i in 0..seg {
            let a0 = 2.0 * std::f64::consts::PI * i as f64 / seg as f64;
            let a1 = 2.0 * std::f64::consts::PI * (i + 1) as f64 / seg as f64;
            self.push_segment(
                [cx + radius * a0.cos(), cy + radius * a0.sin(), cz],
                [cx + radius * a1.cos(), cy + radius * a1.sin(), cz],
                color,
                self.default_width,
            );
        }
        // XZ circle
        for i in 0..seg {
            let a0 = 2.0 * std::f64::consts::PI * i as f64 / seg as f64;
            let a1 = 2.0 * std::f64::consts::PI * (i + 1) as f64 / seg as f64;
            self.push_segment(
                [cx + radius * a0.cos(), cy, cz + radius * a0.sin()],
                [cx + radius * a1.cos(), cy, cz + radius * a1.sin()],
                color,
                self.default_width,
            );
        }
        // YZ circle
        for i in 0..seg {
            let a0 = 2.0 * std::f64::consts::PI * i as f64 / seg as f64;
            let a1 = 2.0 * std::f64::consts::PI * (i + 1) as f64 / seg as f64;
            self.push_segment(
                [cx, cy + radius * a0.cos(), cz + radius * a0.sin()],
                [cx, cy + radius * a1.cos(), cz + radius * a1.sin()],
                color,
                self.default_width,
            );
        }
    }

    /// Draw a sphere with latitude/longitude lines (more detail).
    pub fn draw_sphere_latlon(
        &mut self,
        center: [f64; 3],
        radius: f64,
        lat_lines: usize,
        lon_lines: usize,
        segments: usize,
        color: [f64; 3],
    ) {
        let seg = segments.max(4);
        let [cx, cy, cz] = center;

        // Latitude circles (horizontal rings at different heights)
        for lat in 1..lat_lines {
            let phi = std::f64::consts::PI * lat as f64 / lat_lines as f64;
            let r = radius * phi.sin();
            let y = cy + radius * phi.cos();
            for i in 0..seg {
                let a0 = 2.0 * std::f64::consts::PI * i as f64 / seg as f64;
                let a1 = 2.0 * std::f64::consts::PI * (i + 1) as f64 / seg as f64;
                self.push_segment(
                    [cx + r * a0.cos(), y, cz + r * a0.sin()],
                    [cx + r * a1.cos(), y, cz + r * a1.sin()],
                    color,
                    self.default_width,
                );
            }
        }

        // Longitude lines (vertical meridians)
        for lon in 0..lon_lines {
            let theta = 2.0 * std::f64::consts::PI * lon as f64 / lon_lines as f64;
            for i in 0..seg {
                let phi0 = std::f64::consts::PI * i as f64 / seg as f64;
                let phi1 = std::f64::consts::PI * (i + 1) as f64 / seg as f64;
                self.push_segment(
                    [
                        cx + radius * phi0.sin() * theta.cos(),
                        cy + radius * phi0.cos(),
                        cz + radius * phi0.sin() * theta.sin(),
                    ],
                    [
                        cx + radius * phi1.sin() * theta.cos(),
                        cy + radius * phi1.cos(),
                        cz + radius * phi1.sin() * theta.sin(),
                    ],
                    color,
                    self.default_width,
                );
            }
        }
    }

    /// Draw a capsule between two endpoint centers with the given radius.
    pub fn draw_capsule(
        &mut self,
        p0: [f64; 3],
        p1: [f64; 3],
        radius: f64,
        segments: usize,
        color: [f64; 3],
    ) {
        let seg = segments.max(4);
        // Circles at each end
        for &center in &[p0, p1] {
            for i in 0..seg {
                let a0 = 2.0 * std::f64::consts::PI * i as f64 / seg as f64;
                let a1 = 2.0 * std::f64::consts::PI * (i + 1) as f64 / seg as f64;
                self.push_segment(
                    [
                        center[0] + radius * a0.cos(),
                        center[1],
                        center[2] + radius * a0.sin(),
                    ],
                    [
                        center[0] + radius * a1.cos(),
                        center[1],
                        center[2] + radius * a1.sin(),
                    ],
                    color,
                    self.default_width,
                );
            }
        }
        // 4 vertical connecting lines
        for k in 0..4 {
            let angle = std::f64::consts::FRAC_PI_2 * k as f64;
            let dx = radius * angle.cos();
            let dz = radius * angle.sin();
            self.push_segment(
                [p0[0] + dx, p0[1], p0[2] + dz],
                [p1[0] + dx, p1[1], p1[2] + dz],
                color,
                self.default_width,
            );
        }
    }

    /// Draw a cylinder with top and bottom circles and 4 vertical edge lines.
    pub fn draw_cylinder(
        &mut self,
        center: [f64; 3],
        radius: f64,
        half_height: f64,
        segments: usize,
        color: [f64; 3],
    ) {
        let seg = segments.max(4);
        let [cx, cy, cz] = center;
        let top_y = cy + half_height;
        let bot_y = cy - half_height;
        // Top and bottom circles
        for i in 0..seg {
            let a0 = 2.0 * std::f64::consts::PI * i as f64 / seg as f64;
            let a1 = 2.0 * std::f64::consts::PI * (i + 1) as f64 / seg as f64;
            self.push_segment(
                [cx + radius * a0.cos(), top_y, cz + radius * a0.sin()],
                [cx + radius * a1.cos(), top_y, cz + radius * a1.sin()],
                color,
                self.default_width,
            );
            self.push_segment(
                [cx + radius * a0.cos(), bot_y, cz + radius * a0.sin()],
                [cx + radius * a1.cos(), bot_y, cz + radius * a1.sin()],
                color,
                self.default_width,
            );
        }
        // 4 vertical lines at cardinal points
        for k in 0..4 {
            let angle = std::f64::consts::FRAC_PI_2 * k as f64;
            let dx = radius * angle.cos();
            let dz = radius * angle.sin();
            self.push_segment(
                [cx + dx, bot_y, cz + dz],
                [cx + dx, top_y, cz + dz],
                color,
                self.default_width,
            );
        }
    }

    /// Draw RGB XYZ coordinate axes at `origin` with arrow heads.
    pub fn draw_axes(&mut self, origin: [f64; 3], length: f64) {
        let [ox, oy, oz] = origin;
        let tip_frac = 0.15;
        let tip_len = length * tip_frac;
        let axes = [
            ([1.0, 0.0, 0.0], [ox + length, oy, oz]), // X: red
            ([0.0, 1.0, 0.0], [ox, oy + length, oz]), // Y: green
            ([0.0, 0.0, 1.0], [ox, oy, oz + length]), // Z: blue
        ];
        let colors = [
            [1.0, 0.0, 0.0_f64], // red
            [0.0, 1.0, 0.0_f64], // green
            [0.0, 0.0, 1.0_f64], // blue
        ];
        for (axis_idx, (dir, tip)) in axes.iter().enumerate() {
            let color = colors[axis_idx];
            self.push_segment(origin, *tip, color, self.default_width);
            // Small arrowhead barbs: perpendicular lines at tip
            let perp_a = [
                tip[0] - dir[0] * tip_len + dir[1] * tip_len * 0.5,
                tip[1] - dir[1] * tip_len + dir[2] * tip_len * 0.5,
                tip[2] - dir[2] * tip_len + dir[0] * tip_len * 0.5,
            ];
            let perp_b = [
                tip[0] - dir[0] * tip_len - dir[1] * tip_len * 0.5,
                tip[1] - dir[1] * tip_len - dir[2] * tip_len * 0.5,
                tip[2] - dir[2] * tip_len - dir[0] * tip_len * 0.5,
            ];
            self.push_segment(*tip, perp_a, color, self.default_width);
            self.push_segment(*tip, perp_b, color, self.default_width);
        }
    }

    /// Draw a flat ground-plane grid centered at `origin`.
    ///
    /// `nx` and `nz` are cell counts in X and Z, `spacing` is cell size.
    /// Grid lines run along X and Z axes at y = `origin[1]`.
    pub fn draw_grid(
        &mut self,
        origin: [f64; 3],
        nx: usize,
        nz: usize,
        spacing: f64,
        color: [f64; 3],
    ) {
        let [ox, oy, oz] = origin;
        let half_x = nx as f64 * spacing * 0.5;
        let half_z = nz as f64 * spacing * 0.5;

        // Lines parallel to X axis (varying Z)
        for iz in 0..=nz {
            let z = oz - half_z + iz as f64 * spacing;
            self.push_segment(
                [ox - half_x, oy, z],
                [ox + half_x, oy, z],
                color,
                self.default_width,
            );
        }
        // Lines parallel to Z axis (varying X)
        for ix in 0..=nx {
            let x = ox - half_x + ix as f64 * spacing;
            self.push_segment(
                [x, oy, oz - half_z],
                [x, oy, oz + half_z],
                color,
                self.default_width,
            );
        }
    }

    /// Draw a camera frustum wireframe.
    ///
    /// The frustum is defined by near/far distances and horizontal/vertical half-FOV angles.
    /// The camera looks along -Z by convention.
    pub fn draw_frustum(
        &mut self,
        z_near: f64,
        z_far: f64,
        half_fov_x: f64,
        half_fov_y: f64,
        color: [f64; 3],
    ) {
        let xn = z_near * half_fov_x.tan();
        let yn = z_near * half_fov_y.tan();
        let xf = z_far * half_fov_x.tan();
        let yf = z_far * half_fov_y.tan();

        let near = [
            [-xn, -yn, -z_near],
            [xn, -yn, -z_near],
            [xn, yn, -z_near],
            [-xn, yn, -z_near],
        ];
        let far = [
            [-xf, -yf, -z_far],
            [xf, -yf, -z_far],
            [xf, yf, -z_far],
            [-xf, yf, -z_far],
        ];

        // Near rectangle
        for i in 0..4 {
            self.push_segment(near[i], near[(i + 1) % 4], color, self.default_width);
        }
        // Far rectangle
        for i in 0..4 {
            self.push_segment(far[i], far[(i + 1) % 4], color, self.default_width);
        }
        // Connecting edges
        for i in 0..4 {
            self.push_segment(near[i], far[i], color, self.default_width);
        }
    }

    /// Draw an AABB node (for BVH debugging), coloring deeper levels more red.
    pub fn draw_bvh_node(&mut self, min: [f64; 3], max: [f64; 3], depth: usize) {
        // Color: green at depth 0 → red at depth 8+
        let t = (depth as f64 / 8.0).min(1.0);
        let color = [t, 1.0 - t, 0.0];
        let cx = (min[0] + max[0]) * 0.5;
        let cy = (min[1] + max[1]) * 0.5;
        let cz = (min[2] + max[2]) * 0.5;
        let hx = (max[0] - min[0]) * 0.5;
        let hy = (max[1] - min[1]) * 0.5;
        let hz = (max[2] - min[2]) * 0.5;
        self.draw_box([cx, cy, cz], [hx, hy, hz], color);
    }

    /// Flush the buffer into a `WireframeBatch` (consuming the buffer).
    pub fn flush_to_batch(&mut self, label: impl Into<String>) -> WireframeBatch {
        let mut batch = WireframeBatch::new(label);
        batch.segments.append(&mut self.buffer);
        batch
    }

    /// Total GPU byte size of the buffer.
    pub fn byte_size(&self) -> usize {
        self.buffer.len() * std::mem::size_of::<WireSegmentGpu>()
    }
}

// ─── Tests for WireframeRenderer and new utilities ───────────────────────────

#[cfg(test)]
mod renderer_tests {
    use super::*;
    use crate::scientific_plot::LineStyle;
    use crate::wireframe::VertexBufferData;
    use crate::wireframe::WireframeColorMap;
    use crate::wireframe::WireframeOverlayConfig;
    use crate::wireframe::cull_wireframe_batch;
    use crate::wireframe::primitives::WireframeMesh;
    use crate::wireframe::render_bvh_tree;
    use crate::wireframe::segment_intersects_aabb;

    // ── WireframeRenderer ─────────────────────────────────────────────────────

    #[test]
    fn test_renderer_new_empty() {
        let r = WireframeRenderer::new([1.0, 0.0, 0.0], 1.5);
        assert_eq!(r.segment_count(), 0);
        assert_eq!(r.byte_size(), 0);
    }

    #[test]
    fn test_renderer_push_and_count() {
        let mut r = WireframeRenderer::new([1.0, 1.0, 1.0], 1.0);
        r.push([0.0; 3], [1.0, 0.0, 0.0]);
        r.push([1.0, 0.0, 0.0], [2.0, 0.0, 0.0]);
        assert_eq!(r.segment_count(), 2);
    }

    #[test]
    fn test_renderer_clear() {
        let mut r = WireframeRenderer::new([1.0, 1.0, 1.0], 1.0);
        r.push([0.0; 3], [1.0; 3]);
        r.push([2.0; 3], [3.0; 3]);
        r.clear();
        assert_eq!(r.segment_count(), 0);
    }

    #[test]
    fn test_renderer_draw_box_produces_12_segments() {
        let mut r = WireframeRenderer::new([1.0, 1.0, 1.0], 1.0);
        r.draw_box([0.0; 3], [1.0; 3], [1.0, 1.0, 1.0]);
        assert_eq!(r.segment_count(), 12);
    }

    #[test]
    fn test_renderer_draw_sphere_3_great_circles() {
        let mut r = WireframeRenderer::new([1.0, 1.0, 1.0], 1.0);
        let seg = 16;
        r.draw_sphere([0.0; 3], 1.0, seg, [1.0, 1.0, 0.0]);
        // 3 great circles × 16 segments
        assert_eq!(r.segment_count(), 3 * seg);
    }

    #[test]
    fn test_renderer_draw_sphere_latlon_segment_count() {
        let mut r = WireframeRenderer::new([1.0, 1.0, 1.0], 1.0);
        let lat = 4;
        let lon = 6;
        let seg = 8;
        // lat-1 rings (interior only) × seg + lon lines × seg
        r.draw_sphere_latlon([0.0; 3], 1.0, lat, lon, seg, [0.0, 1.0, 0.0]);
        let expected = (lat - 1) * seg + lon * seg;
        assert_eq!(r.segment_count(), expected);
    }

    #[test]
    fn test_renderer_draw_capsule_segment_count() {
        let mut r = WireframeRenderer::new([1.0, 1.0, 1.0], 1.0);
        let seg = 8;
        r.draw_capsule([0.0, 0.0, 0.0], [0.0, 2.0, 0.0], 0.5, seg, [0.5, 0.5, 1.0]);
        // 2 circles × seg + 4 vertical lines
        assert_eq!(r.segment_count(), 2 * seg + 4);
    }

    #[test]
    fn test_renderer_draw_cylinder_segment_count() {
        let mut r = WireframeRenderer::new([1.0, 1.0, 1.0], 1.0);
        let seg = 12;
        r.draw_cylinder([0.0; 3], 1.0, 1.5, seg, [0.8, 0.2, 0.2]);
        // 2 circles × seg + 4 vertical lines
        assert_eq!(r.segment_count(), 2 * seg + 4);
    }

    #[test]
    fn test_renderer_draw_axes_segment_count() {
        let mut r = WireframeRenderer::new([1.0, 1.0, 1.0], 1.0);
        r.draw_axes([0.0; 3], 1.0);
        // 3 axes × (1 shaft + 2 arrowhead barbs) = 9 segments
        assert_eq!(r.segment_count(), 9);
    }

    #[test]
    fn test_renderer_draw_axes_colors() {
        let mut r = WireframeRenderer::new([1.0, 1.0, 1.0], 1.0);
        r.draw_axes([0.0; 3], 1.0);
        // X axis: red, Y axis: green, Z axis: blue
        // First 3 segments (shaft + 2 barbs) are for X
        let x_shaft = &r.buffer[0];
        assert!(
            (x_shaft.color[0] - 1.0).abs() < 1e-5,
            "X shaft should be red"
        );
        assert!(x_shaft.color[1] < 0.1, "X shaft should be red (no green)");
    }

    #[test]
    fn test_renderer_draw_grid_segment_count() {
        let mut r = WireframeRenderer::new([0.5, 0.5, 0.5], 1.0);
        let nx = 4;
        let nz = 3;
        r.draw_grid([0.0; 3], nx, nz, 1.0, [0.5, 0.5, 0.5]);
        // (nz+1) lines along X + (nx+1) lines along Z
        let expected = (nz + 1) + (nx + 1);
        assert_eq!(r.segment_count(), expected);
    }

    #[test]
    fn test_renderer_draw_frustum_segment_count() {
        let mut r = WireframeRenderer::new([1.0, 1.0, 0.0], 1.0);
        let hfov = std::f64::consts::FRAC_PI_4;
        r.draw_frustum(0.1, 100.0, hfov, hfov * 0.75, [1.0, 1.0, 0.0]);
        // 4 near + 4 far + 4 connecting = 12 segments
        assert_eq!(r.segment_count(), 12);
    }

    #[test]
    fn test_renderer_draw_bvh_node_segment_count() {
        let mut r = WireframeRenderer::new([1.0, 1.0, 1.0], 1.0);
        r.draw_bvh_node([0.0; 3], [1.0; 3], 0);
        // A box has 12 edges
        assert_eq!(r.segment_count(), 12);
    }

    #[test]
    fn test_renderer_bvh_node_depth_color() {
        let mut r0 = WireframeRenderer::new([1.0, 1.0, 1.0], 1.0);
        let mut r8 = WireframeRenderer::new([1.0, 1.0, 1.0], 1.0);
        r0.draw_bvh_node([0.0; 3], [1.0; 3], 0);
        r8.draw_bvh_node([0.0; 3], [1.0; 3], 8);
        // Depth 0 → green, depth 8 → red
        let color_0 = r0.buffer[0].color;
        let color_8 = r8.buffer[0].color;
        assert!(color_0[1] > color_8[1], "shallow node should be greener");
        assert!(color_8[0] > color_0[0], "deep node should be redder");
    }

    #[test]
    fn test_renderer_flush_to_batch() {
        let mut r = WireframeRenderer::new([1.0, 1.0, 1.0], 1.0);
        r.push([0.0; 3], [1.0; 3]);
        r.push([2.0; 3], [3.0; 3]);
        assert_eq!(r.segment_count(), 2);
        let batch = r.flush_to_batch("test_batch");
        // Flush should clear the renderer buffer
        assert_eq!(r.segment_count(), 0);
        assert_eq!(batch.len(), 2);
        assert_eq!(batch.label, "test_batch");
    }

    #[test]
    fn test_renderer_byte_size() {
        let mut r = WireframeRenderer::new([1.0, 1.0, 1.0], 1.0);
        r.push([0.0; 3], [1.0; 3]);
        let seg_size = std::mem::size_of::<WireSegmentGpu>();
        assert_eq!(r.byte_size(), seg_size);
    }

    // ── LineStyle ─────────────────────────────────────────────────────────────

    #[test]
    fn test_line_style_solid_no_dash() {
        assert!(LineStyle::Solid.svg_dash_array().is_none());
    }

    #[test]
    fn test_line_style_dashed_has_dash() {
        assert!(LineStyle::Dashed.svg_dash_array().is_some());
    }

    #[test]
    fn test_line_style_dotted_has_pattern() {
        let pat = LineStyle::Dotted.svg_dash_array();
        assert!(pat.is_some());
        assert!(pat.unwrap_or("").contains(','));
    }

    #[test]
    fn test_line_style_dashdot_pattern() {
        let pat = LineStyle::DashDot.svg_dash_array();
        assert!(pat.is_some());
        // DashDot should have at least two comma-separated segments
        let s = pat.unwrap_or("");
        assert!(s.matches(',').count() >= 2, "DashDot pattern: {s}");
    }

    // ── VertexBufferData ──────────────────────────────────────────────────────

    #[test]
    fn test_vertex_buffer_data_from_batch() {
        let mesh = WireframeMesh::from_aabb([0.0; 3], [1.0; 3]);
        let batch = WireframeBatch::from_mesh(&mesh, 1.0);
        let vbd = VertexBufferData::from_batch(&batch);
        // 12 segments × 2 vertices = 24 vertices
        assert_eq!(vbd.vertex_count, 24);
        // 24 vertices × 6 floats = 144 floats
        assert_eq!(vbd.data.len(), 24 * 6);
    }

    #[test]
    fn test_vertex_buffer_data_stride() {
        assert_eq!(VertexBufferData::stride_bytes(), 24);
    }

    #[test]
    fn test_vertex_buffer_data_offsets() {
        assert_eq!(VertexBufferData::position_offset(), 0);
        assert_eq!(VertexBufferData::color_offset(), 12);
    }

    #[test]
    fn test_vertex_buffer_data_byte_size() {
        let mesh = WireframeMesh::from_aabb([0.0; 3], [1.0; 3]);
        let batch = WireframeBatch::from_mesh(&mesh, 1.0);
        let vbd = VertexBufferData::from_batch(&batch);
        assert_eq!(vbd.byte_size(), vbd.data.len() * 4);
    }

    // ── WireframeColorMap ─────────────────────────────────────────────────────

    #[test]
    fn test_color_map_empty_returns_black() {
        let cmap = WireframeColorMap::new();
        assert_eq!(cmap.sample(0.5), [0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_color_map_heat_endpoints() {
        let cmap = WireframeColorMap::heat();
        let c0 = cmap.sample(0.0);
        let c1 = cmap.sample(1.0);
        // At 0.0: blue
        assert!((c0[2] - 1.0).abs() < 1e-9, "heat map t=0 should be blue");
        assert!((c0[0]).abs() < 1e-9, "heat map t=0 should have no red");
        // At 1.0: red
        assert!((c1[0] - 1.0).abs() < 1e-9, "heat map t=1 should be red");
        assert!((c1[2]).abs() < 1e-9, "heat map t=1 should have no blue");
    }

    #[test]
    fn test_color_map_grayscale_midpoint() {
        let cmap = WireframeColorMap::grayscale();
        let mid = cmap.sample(0.5);
        for c in mid {
            assert!(
                (c - 0.5).abs() < 1e-9,
                "grayscale midpoint should be 0.5, got {c}"
            );
        }
    }

    #[test]
    fn test_color_map_clamped_beyond_range() {
        let cmap = WireframeColorMap::heat();
        let below = cmap.sample(-1.0);
        let above = cmap.sample(2.0);
        assert_eq!(below, cmap.sample(0.0));
        assert_eq!(above, cmap.sample(1.0));
    }

    #[test]
    fn test_color_map_add_stop_interpolates() {
        let mut cmap = WireframeColorMap::new();
        cmap.add_stop(0.0, [0.0, 0.0, 0.0]);
        cmap.add_stop(1.0, [1.0, 1.0, 1.0]);
        let mid = cmap.sample(0.5);
        for c in mid {
            assert!((c - 0.5).abs() < 1e-9);
        }
    }

    // ── segment_intersects_aabb ───────────────────────────────────────────────

    #[test]
    fn test_segment_aabb_intersects_through_center() {
        let p = [-2.0, 0.0, 0.0];
        let q = [2.0, 0.0, 0.0];
        let mn = [-1.0, -1.0, -1.0];
        let mx = [1.0, 1.0, 1.0];
        assert!(segment_intersects_aabb(p, q, mn, mx));
    }

    #[test]
    fn test_segment_aabb_no_intersect() {
        let p = [2.0, 2.0, 2.0];
        let q = [3.0, 3.0, 3.0];
        let mn = [-1.0, -1.0, -1.0];
        let mx = [1.0, 1.0, 1.0];
        assert!(!segment_intersects_aabb(p, q, mn, mx));
    }

    #[test]
    fn test_segment_aabb_endpoint_inside() {
        let p = [0.0, 0.0, 0.0]; // inside
        let q = [5.0, 0.0, 0.0]; // outside
        let mn = [-1.0; 3];
        let mx = [1.0; 3];
        assert!(segment_intersects_aabb(p, q, mn, mx));
    }

    #[test]
    fn test_segment_aabb_parallel_outside_slab() {
        // Segment parallel to X, outside Y slab
        let p = [0.0, 5.0, 0.0];
        let q = [3.0, 5.0, 0.0];
        let mn = [0.0, 0.0, 0.0];
        let mx = [3.0, 3.0, 3.0];
        assert!(!segment_intersects_aabb(p, q, mn, mx));
    }

    // ── cull_wireframe_batch ──────────────────────────────────────────────────

    #[test]
    fn test_cull_batch_all_visible() {
        let mesh = WireframeMesh::from_aabb([-0.5; 3], [0.5; 3]);
        let batch = WireframeBatch::from_mesh(&mesh, 1.0);
        // Huge planes that contain everything
        let planes: [[f64; 4]; 6] = [
            [1.0, 0.0, 0.0, 1000.0],
            [-1.0, 0.0, 0.0, 1000.0],
            [0.0, 1.0, 0.0, 1000.0],
            [0.0, -1.0, 0.0, 1000.0],
            [0.0, 0.0, 1.0, 1000.0],
            [0.0, 0.0, -1.0, 1000.0],
        ];
        let culled = cull_wireframe_batch(&batch, &planes);
        assert_eq!(
            culled.len(),
            batch.len(),
            "all segments should survive huge frustum"
        );
    }

    #[test]
    fn test_cull_batch_all_culled() {
        let mesh = WireframeMesh::from_aabb([10.0; 3], [11.0; 3]);
        let batch = WireframeBatch::from_mesh(&mesh, 1.0);
        // Planes that only pass a tiny region around the origin
        let planes: [[f64; 4]; 6] = [
            [1.0, 0.0, 0.0, 1.0],  // x >= -1
            [-1.0, 0.0, 0.0, 1.0], // x <=  1
            [0.0, 1.0, 0.0, 1.0],
            [0.0, -1.0, 0.0, 1.0],
            [0.0, 0.0, 1.0, 1.0],
            [0.0, 0.0, -1.0, 1.0],
        ];
        let culled = cull_wireframe_batch(&batch, &planes);
        assert_eq!(
            culled.len(),
            0,
            "all segments should be culled away from origin"
        );
    }

    // ── render_bvh_tree ───────────────────────────────────────────────────────

    #[test]
    fn test_render_bvh_tree_single_node() {
        let nodes = vec![([0.0f64; 3], [1.0f64; 3], usize::MAX, usize::MAX)];
        let mut r = WireframeRenderer::new([1.0; 3], 1.0);
        render_bvh_tree(&nodes, 0, 4, &mut r);
        assert_eq!(
            r.segment_count(),
            12,
            "single node should produce 12 box segments"
        );
    }

    #[test]
    fn test_render_bvh_tree_three_nodes() {
        // root → left, root → right (no grandchildren)
        let nodes = vec![
            ([0.0; 3], [4.0; 3], 1, 2),                    // root
            ([-2.0; 3], [0.0; 3], usize::MAX, usize::MAX), // left
            ([0.0; 3], [2.0; 3], usize::MAX, usize::MAX),  // right
        ];
        let mut r = WireframeRenderer::new([1.0; 3], 1.0);
        render_bvh_tree(&nodes, 0, 4, &mut r);
        // 3 nodes × 12 = 36 segments
        assert_eq!(r.segment_count(), 36);
    }

    #[test]
    fn test_render_bvh_tree_depth_limit() {
        // 3-level tree (root depth 0, children depth 1, grandchildren depth 2)
        let nodes = vec![
            ([0.0; 3], [8.0; 3], 1, 2),                   // 0: root
            ([0.0; 3], [4.0; 3], 3, 4),                   // 1: left
            ([4.0; 3], [8.0; 3], usize::MAX, usize::MAX), // 2: right
            ([0.0; 3], [2.0; 3], usize::MAX, usize::MAX), // 3: grandchild
            ([2.0; 3], [4.0; 3], usize::MAX, usize::MAX), // 4: grandchild
        ];
        let mut r_full = WireframeRenderer::new([1.0; 3], 1.0);
        let mut r_limited = WireframeRenderer::new([1.0; 3], 1.0);
        render_bvh_tree(&nodes, 0, 4, &mut r_full);
        render_bvh_tree(&nodes, 0, 0, &mut r_limited); // only depth 0
        assert!(r_full.segment_count() > r_limited.segment_count());
        assert_eq!(r_limited.segment_count(), 12); // just the root box
    }

    // ── WireframeOverlayConfig ────────────────────────────────────────────────

    #[test]
    fn test_wireframe_overlay_config_defaults() {
        let cfg = WireframeOverlayConfig::default_config();
        assert!(cfg.base_width_px > 0.0);
        assert!(cfg.feature_width_px > cfg.base_width_px);
        assert!(cfg.silhouette_width_px > cfg.feature_width_px);
        assert!(cfg.feature_threshold_deg > 0.0);
    }
}
