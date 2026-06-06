// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Line style, vertex buffer, color map, overlay config, segment utilities,
//! wireframe segment/grid/depth sort, and coordinate axis helpers.

use super::mesh_ops::WireframeBatch;
use super::renderer::WireframeRenderer;

// ─── LineColorThickness ───────────────────────────────────────────────────────

/// Helper to describe the visual style of a wireframe line.
#[derive(Debug, Clone, Copy)]
pub struct LineStyle {
    /// RGB color (0.0–1.0 per channel).
    pub color: [f32; 3],
    /// Alpha/opacity (0.0–1.0).
    pub alpha: f32,
    /// Screen-space line width in pixels.
    pub width_px: f32,
    /// Whether this line should be drawn on top (ignore depth).
    pub always_on_top: bool,
}

impl LineStyle {
    /// Solid white 1-px line, depth tested.
    pub fn white() -> Self {
        Self {
            color: [1.0, 1.0, 1.0],
            alpha: 1.0,
            width_px: 1.0,
            always_on_top: false,
        }
    }

    /// Red 1-px line, depth tested.
    pub fn red() -> Self {
        Self {
            color: [1.0, 0.0, 0.0],
            alpha: 1.0,
            width_px: 1.0,
            always_on_top: false,
        }
    }

    /// Green 1-px line.
    pub fn green() -> Self {
        Self {
            color: [0.0, 1.0, 0.0],
            alpha: 1.0,
            width_px: 1.0,
            always_on_top: false,
        }
    }

    /// Blue 1-px line.
    pub fn blue() -> Self {
        Self {
            color: [0.0, 0.0, 1.0],
            alpha: 1.0,
            width_px: 1.0,
            always_on_top: false,
        }
    }

    /// Yellow 2-px line, always on top.
    pub fn highlight() -> Self {
        Self {
            color: [1.0, 1.0, 0.0],
            alpha: 1.0,
            width_px: 2.0,
            always_on_top: true,
        }
    }

    /// Semi-transparent gray 1-px line for background geometry.
    pub fn dimmed() -> Self {
        Self {
            color: [0.5, 0.5, 0.5],
            alpha: 0.4,
            width_px: 1.0,
            always_on_top: false,
        }
    }

    /// Set width.
    pub fn with_width(mut self, w: f32) -> Self {
        self.width_px = w;
        self
    }

    /// Set alpha.
    pub fn with_alpha(mut self, a: f32) -> Self {
        self.alpha = a;
        self
    }

    /// Convert to an RGBA array.
    pub fn rgba(&self) -> [f32; 4] {
        [self.color[0], self.color[1], self.color[2], self.alpha]
    }
}

// ─── VertexBufferData ─────────────────────────────────────────────────────────

/// Flat interleaved vertex buffer for GPU upload.
///
/// Each vertex occupies 6 × f32: `[x, y, z, r, g, b]`.
/// Suitable for GL_LINES or Vulkan vertex attribute upload.
#[derive(Debug, Clone)]
pub struct VertexBufferData {
    /// Flat float data: 6 f32 per vertex (xyz + rgb).
    pub data: Vec<f32>,
    /// Number of vertices.
    pub vertex_count: usize,
}

impl VertexBufferData {
    /// Build a `VertexBufferData` from a `WireframeBatch`.
    ///
    /// Each segment contributes 2 vertices (start and end), each using the
    /// segment's color for both endpoints.
    pub fn from_batch(batch: &WireframeBatch) -> Self {
        let mut data = Vec::with_capacity(batch.segments.len() * 12);
        for seg in &batch.segments {
            // Start vertex
            data.extend_from_slice(&seg.start);
            data.extend_from_slice(&seg.color[..3]);
            // End vertex
            data.extend_from_slice(&seg.end);
            data.extend_from_slice(&seg.color[..3]);
        }
        let vertex_count = batch.segments.len() * 2;
        Self { data, vertex_count }
    }

    /// Byte size of the buffer (for GPU buffer allocation).
    pub fn byte_size(&self) -> usize {
        self.data.len() * std::mem::size_of::<f32>()
    }

    /// Stride in bytes between consecutive vertices (6 × 4 = 24 bytes).
    pub fn stride_bytes() -> usize {
        6 * std::mem::size_of::<f32>()
    }

    /// Byte offset of the position attribute within each vertex (0).
    pub fn position_offset() -> usize {
        0
    }

    /// Byte offset of the color attribute within each vertex (3 × 4 = 12 bytes).
    pub fn color_offset() -> usize {
        3 * std::mem::size_of::<f32>()
    }
}

// ─── WireframeColorMap ────────────────────────────────────────────────────────

/// A simple color map for mapping scalar values to wire colors.
#[derive(Debug, Clone)]
pub struct WireframeColorMap {
    /// Control points: each is `(value, [r, g, b])`.
    stops: Vec<(f64, [f64; 3])>,
}

impl Default for WireframeColorMap {
    fn default() -> Self {
        Self::new()
    }
}

impl WireframeColorMap {
    /// Create an empty color map (will return black for all queries).
    pub fn new() -> Self {
        Self { stops: Vec::new() }
    }

    /// Create a heat map (blue → cyan → green → yellow → red).
    pub fn heat() -> Self {
        let stops = vec![
            (0.0, [0.0, 0.0, 1.0]),
            (0.25, [0.0, 1.0, 1.0]),
            (0.5, [0.0, 1.0, 0.0]),
            (0.75, [1.0, 1.0, 0.0]),
            (1.0, [1.0, 0.0, 0.0]),
        ];
        Self { stops }
    }

    /// Create a grayscale color map (black → white).
    pub fn grayscale() -> Self {
        let stops = vec![(0.0, [0.0, 0.0, 0.0]), (1.0, [1.0, 1.0, 1.0])];
        Self { stops }
    }

    /// Add a control stop at the given normalized value (0.0–1.0).
    pub fn add_stop(&mut self, value: f64, color: [f64; 3]) {
        self.stops.push((value, color));
        self.stops
            .sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    }

    /// Sample the color map at a normalized value (0.0–1.0).
    ///
    /// Linearly interpolates between adjacent stops.
    pub fn sample(&self, t: f64) -> [f64; 3] {
        if self.stops.is_empty() {
            return [0.0, 0.0, 0.0];
        }
        if self.stops.len() == 1 {
            return self.stops[0].1;
        }
        let t = t.clamp(0.0, 1.0);
        // Find bracketing stops
        for i in 0..self.stops.len() - 1 {
            let (v0, c0) = self.stops[i];
            let (v1, c1) = self.stops[i + 1];
            if t <= v1 {
                let span = (v1 - v0).max(1e-30);
                let f = (t - v0) / span;
                return [
                    c0[0] + f * (c1[0] - c0[0]),
                    c0[1] + f * (c1[1] - c0[1]),
                    c0[2] + f * (c1[2] - c0[2]),
                ];
            }
        }
        self.stops
            .last()
            .expect("stops is non-empty (checked above)")
            .1
    }
}

// ─── WireframeOverlayConfig ───────────────────────────────────────────────────

/// Configuration for generating layered wireframe overlays.
#[derive(Debug, Clone, Copy)]
pub struct WireframeOverlayConfig {
    /// Width of base (all) edges in pixels.
    pub base_width_px: f64,
    /// Width of feature edges in pixels.
    pub feature_width_px: f64,
    /// Width of silhouette edges in pixels.
    pub silhouette_width_px: f64,
    /// Feature edge detection threshold in degrees.
    pub feature_threshold_deg: f64,
    /// Color of base edges.
    pub base_color: [f64; 3],
    /// Color of feature edges.
    pub feature_color: [f64; 3],
    /// Color of silhouette edges.
    pub silhouette_color: [f64; 3],
}

impl WireframeOverlayConfig {
    /// Default settings: subtle base edges, yellow features, white silhouette.
    pub fn default_config() -> Self {
        Self {
            base_width_px: 1.0,
            feature_width_px: 1.5,
            silhouette_width_px: 2.5,
            feature_threshold_deg: 30.0,
            base_color: [0.5, 0.5, 0.5],
            feature_color: [1.0, 0.9, 0.0],
            silhouette_color: [1.0, 1.0, 1.0],
        }
    }
}

// ─── Utility: line segment AABB intersection ─────────────────────────────────

/// Test whether a 3-D line segment (from `p` to `q`) intersects an AABB
/// defined by `[min, max]` corners.  Returns `true` if the segment overlaps.
///
/// Uses the slab method (Smits' algorithm).
pub fn segment_intersects_aabb(
    p: [f64; 3],
    q: [f64; 3],
    aabb_min: [f64; 3],
    aabb_max: [f64; 3],
) -> bool {
    let dir = [q[0] - p[0], q[1] - p[1], q[2] - p[2]];
    let mut t_min = 0.0_f64;
    let mut t_max = 1.0_f64;

    for i in 0..3 {
        if dir[i].abs() < 1e-12 {
            // Parallel to slab; check if origin is inside
            if p[i] < aabb_min[i] || p[i] > aabb_max[i] {
                return false;
            }
        } else {
            let inv_d = 1.0 / dir[i];
            let mut t1 = (aabb_min[i] - p[i]) * inv_d;
            let mut t2 = (aabb_max[i] - p[i]) * inv_d;
            if t1 > t2 {
                std::mem::swap(&mut t1, &mut t2);
            }
            t_min = t_min.max(t1);
            t_max = t_max.min(t2);
            if t_min > t_max {
                return false;
            }
        }
    }
    true
}

// ─── Frustum culling for wireframe batches ─────────────────────────────────

/// Cull a `WireframeBatch` against a frustum defined by 6 plane equations.
///
/// Each plane is `[nx, ny, nz, d]` where the inward-facing half-space is
/// `nx*x + ny*y + nz*z + d >= 0`.  Segments with both endpoints outside any
/// single plane are discarded.
///
/// Returns a new `WireframeBatch` with only potentially-visible segments.
pub fn cull_wireframe_batch(batch: &WireframeBatch, planes: &[[f64; 4]; 6]) -> WireframeBatch {
    let mut out = WireframeBatch::new(format!("{}_culled", batch.label));
    'seg: for seg in &batch.segments {
        let s = [
            seg.start[0] as f64,
            seg.start[1] as f64,
            seg.start[2] as f64,
        ];
        let e = [seg.end[0] as f64, seg.end[1] as f64, seg.end[2] as f64];
        for &[nx, ny, nz, d] in planes.iter() {
            let ds = nx * s[0] + ny * s[1] + nz * s[2] + d;
            let de = nx * e[0] + ny * e[1] + nz * e[2] + d;
            // Both endpoints on wrong side of plane → cull
            if ds < 0.0 && de < 0.0 {
                continue 'seg;
            }
        }
        out.segments.push(*seg);
    }
    out
}

// ─── BVH Visualizer ───────────────────────────────────────────────────────────

/// Accumulate AABB wireframes for a BVH tree node and its subtree up to `max_depth`.
///
/// `nodes` is a flat array of `(min, max, left_child, right_child)` where
/// child indices of `usize::MAX` mean "no child".
pub fn render_bvh_tree(
    nodes: &[([f64; 3], [f64; 3], usize, usize)],
    root: usize,
    max_depth: usize,
    renderer: &mut WireframeRenderer,
) {
    if root >= nodes.len() {
        return;
    }
    // Iterative DFS using a stack of (node_idx, depth)
    let mut stack = vec![(root, 0usize)];
    while let Some((idx, depth)) = stack.pop() {
        if depth > max_depth || idx >= nodes.len() {
            continue;
        }
        let (mn, mx, left, right) = nodes[idx];
        renderer.draw_bvh_node(mn, mx, depth);
        if left != usize::MAX {
            stack.push((left, depth + 1));
        }
        if right != usize::MAX {
            stack.push((right, depth + 1));
        }
    }
}

#[cfg(test)]
mod tests {

    use crate::wireframe::primitives::WireframeMesh;
    use crate::wireframe::*;

    #[test]
    fn test_wireframe_mesh_from_aabb_12_edges() {
        let mesh = WireframeMesh::from_aabb([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert_eq!(mesh.vertices.len(), 8);
        assert_eq!(mesh.lines.len(), 12);
    }

    #[test]
    fn test_wireframe_mesh_from_sphere() {
        let mesh = WireframeMesh::from_sphere([0.0, 0.0, 0.0], 1.0, 3, 8);
        assert_eq!(mesh.vertices.len(), 24);
        assert_eq!(mesh.lines.len(), 24);
    }

    #[test]
    fn test_wireframe_mesh_from_capsule() {
        let mesh = WireframeMesh::from_capsule([0.0, 0.0, 0.0], [0.0, 2.0, 0.0], 0.5, 8);
        assert_eq!(mesh.vertices.len(), 16);
        assert_eq!(mesh.lines.len(), 20);
    }

    #[test]
    fn test_wireframe_mesh_from_grid() {
        let mesh = WireframeMesh::from_grid(4, 4, 1.0);
        assert_eq!(mesh.lines.len(), 10);
        assert_eq!(mesh.vertices.len(), 20);
    }

    #[test]
    fn test_wireframe_mesh_merge() {
        let a = WireframeMesh::from_aabb([0.0; 3], [1.0; 3]);
        let b = WireframeMesh::from_aabb([2.0; 3], [3.0; 3]);
        let merged = a.merge(&b);
        assert_eq!(merged.vertices.len(), 16);
        assert_eq!(merged.lines.len(), 24);
        for l in &merged.lines {
            assert!(l.start < merged.vertices.len());
            assert!(l.end < merged.vertices.len());
        }
    }

    #[test]
    fn test_debug_draw_line() {
        let mut dd = DebugDraw::new();
        dd.draw_line([0.0; 3], [1.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert_eq!(dd.line_count(), 1);
    }

    #[test]
    fn test_debug_draw_cross() {
        let mut dd = DebugDraw::new();
        dd.draw_cross([0.0; 3], 1.0, [1.0, 1.0, 1.0]);
        assert_eq!(dd.line_count(), 3);
    }

    #[test]
    fn test_debug_draw_circle() {
        let mut dd = DebugDraw::new();
        dd.draw_circle([0.0; 3], [0.0, 1.0, 0.0], 1.0, 16, [0.0, 1.0, 0.0]);
        assert_eq!(dd.line_count(), 16);
    }

    #[test]
    fn test_debug_draw_clear() {
        let mut dd = DebugDraw::new();
        dd.draw_line([0.0; 3], [1.0; 3], [1.0; 3]);
        dd.draw_line([0.0; 3], [2.0; 3], [1.0; 3]);
        assert_eq!(dd.line_count(), 2);
        dd.clear();
        assert_eq!(dd.line_count(), 0);
    }

    #[test]
    fn test_debug_draw_arrow() {
        let mut dd = DebugDraw::new();
        dd.draw_arrow([0.0; 3], [1.0, 0.0, 0.0], 2.0, [1.0, 0.0, 0.0]);
        assert_eq!(dd.line_count(), 2);
    }
}

// ─── Tests for WireframeRenderer and new utilities ───────────────────────────

#[cfg(test)]
mod extra_wire_tests {
    use super::*;
    use crate::Color;
    use crate::contact_normal_lines;
    use crate::wireframe::AaWireframeParams;
    use crate::wireframe::DebugDraw;
    use crate::wireframe::RichWireframeMesh;
    use crate::wireframe::WireSegmentGpu;
    use crate::wireframe::WireframeBuilder;
    use crate::wireframe::WireframeEdge;
    use crate::wireframe::aabb_wireframe_rich;
    use crate::wireframe::build_aa_feature_edge_batch;
    use crate::wireframe::depth_fade;
    use crate::wireframe::detect_feature_edges;
    use crate::wireframe::extract_mesh_edges;
    use crate::wireframe::frustum_wireframe;
    use crate::wireframe::frustum_wireframe_from_fov;
    use crate::wireframe::primitives::WireframeMesh;

    // ── WireSegmentGpu ────────────────────────────────────────────────────────

    #[test]
    fn test_wire_segment_gpu_midpoint() {
        let seg = WireSegmentGpu::new([0.0, 0.0, 0.0], [2.0, 4.0, 6.0], [1.0, 1.0, 1.0], 1.0);
        let mid = seg.midpoint();
        assert!((mid[0] - 1.0).abs() < 1e-5, "midpoint x should be 1.0");
        assert!((mid[1] - 2.0).abs() < 1e-5, "midpoint y should be 2.0");
        assert!((mid[2] - 3.0).abs() < 1e-5, "midpoint z should be 3.0");
    }

    #[test]
    fn test_wire_segment_gpu_zero_length() {
        let seg = WireSegmentGpu::new([1.0, 2.0, 3.0], [1.0, 2.0, 3.0], [1.0, 1.0, 1.0], 1.0);
        assert!(
            seg.length() < 1e-5,
            "degenerate segment length should be ~0"
        );
    }

    #[test]
    fn test_wire_segment_gpu_color_stored() {
        let color = [0.2, 0.5, 0.8];
        let seg = WireSegmentGpu::new([0.0; 3], [1.0; 3], color, 2.0);
        assert!((seg.color[0] - 0.2_f32).abs() < 1e-4);
        assert!((seg.color[1] - 0.5_f32).abs() < 1e-4);
        assert!((seg.color[2] - 0.8_f32).abs() < 1e-4);
        assert!(
            (seg.color[3] - 1.0_f32).abs() < 1e-5,
            "alpha should default to 1.0"
        );
    }

    #[test]
    fn test_wire_segment_gpu_half_width() {
        let seg = WireSegmentGpu::new([0.0; 3], [1.0; 3], [1.0; 3], 4.0);
        assert!(
            (seg.half_width - 2.0_f32).abs() < 1e-5,
            "half_width should be width/2"
        );
    }

    // ── WireframeBatch ────────────────────────────────────────────────────────

    #[test]
    fn test_wireframe_batch_new_empty_label() {
        let batch = WireframeBatch::new("my_batch");
        assert_eq!(batch.label, "my_batch");
        assert!(batch.is_empty());
        assert_eq!(batch.len(), 0);
    }

    #[test]
    fn test_wireframe_batch_add_increments_len() {
        let mut batch = WireframeBatch::new("test");
        batch.add([0.0; 3], [1.0; 3], [1.0; 3], 1.0);
        assert_eq!(batch.len(), 1);
        batch.add([1.0; 3], [2.0; 3], [0.0; 3], 1.0);
        assert_eq!(batch.len(), 2);
    }

    #[test]
    fn test_wireframe_batch_byte_size_proportional() {
        let mut batch = WireframeBatch::new("test");
        batch.add([0.0; 3], [1.0; 3], [1.0; 3], 1.0);
        let single_size = batch.byte_size();
        batch.add([1.0; 3], [2.0; 3], [0.5; 3], 1.0);
        assert_eq!(batch.byte_size(), 2 * single_size);
    }

    #[test]
    fn test_wireframe_batch_is_empty_after_no_adds() {
        let batch: WireframeBatch = WireframeBatch::new("empty");
        assert!(batch.is_empty());
        assert_eq!(batch.byte_size(), 0);
    }

    // ── WireframeMesh ─────────────────────────────────────────────────────────

    #[test]
    fn test_wireframe_mesh_from_aabb_vertex_positions() {
        let min = [0.0, 0.0, 0.0];
        let max = [2.0, 3.0, 4.0];
        let mesh = WireframeMesh::from_aabb(min, max);
        // All 8 corners should be present
        assert_eq!(mesh.vertices.len(), 8);
        let has_min = mesh.vertices.iter().any(|v| {
            (v.position[0] - 0.0).abs() < 1e-9
                && (v.position[1] - 0.0).abs() < 1e-9
                && (v.position[2] - 0.0).abs() < 1e-9
        });
        let has_max = mesh.vertices.iter().any(|v| {
            (v.position[0] - 2.0).abs() < 1e-9
                && (v.position[1] - 3.0).abs() < 1e-9
                && (v.position[2] - 4.0).abs() < 1e-9
        });
        assert!(has_min, "min corner should be a vertex");
        assert!(has_max, "max corner should be a vertex");
    }

    #[test]
    fn test_wireframe_mesh_from_sphere_line_bounds() {
        let mesh = WireframeMesh::from_sphere([0.0; 3], 1.0, 4, 8);
        for line in &mesh.lines {
            assert!(line.start < mesh.vertices.len(), "line start out of bounds");
            assert!(line.end < mesh.vertices.len(), "line end out of bounds");
        }
    }

    #[test]
    fn test_wireframe_mesh_from_grid_vertex_count() {
        let nx = 5;
        let nz = 3;
        let mesh = WireframeMesh::from_grid(nx, nz, 1.0);
        // Lines: (nz+1) + (nx+1); each has 2 vertices
        let expected_lines = (nz + 1) + (nx + 1);
        assert_eq!(mesh.lines.len(), expected_lines);
        assert_eq!(mesh.vertices.len(), expected_lines * 2);
    }

    #[test]
    fn test_wireframe_mesh_from_capsule_line_bounds() {
        let mesh = WireframeMesh::from_capsule([0.0; 3], [0.0, 3.0, 0.0], 1.0, 12);
        for line in &mesh.lines {
            assert!(line.start < mesh.vertices.len());
            assert!(line.end < mesh.vertices.len());
        }
    }

    #[test]
    fn test_wireframe_mesh_merge_vertex_indices_valid() {
        let a = WireframeMesh::from_aabb([0.0; 3], [1.0; 3]);
        let b = WireframeMesh::from_sphere([5.0; 3], 1.0, 2, 4);
        let merged = a.merge(&b);
        for l in &merged.lines {
            assert!(l.start < merged.vertices.len());
            assert!(l.end < merged.vertices.len());
        }
    }

    // ── DebugDraw ─────────────────────────────────────────────────────────────

    #[test]
    fn test_debug_draw_circle_line_count_matches_segments() {
        let mut dd = DebugDraw::new();
        let seg = 20;
        dd.draw_circle([0.0; 3], [0.0, 1.0, 0.0], 2.0, seg, [1.0, 0.0, 0.0]);
        assert_eq!(
            dd.line_count(),
            seg,
            "circle should produce exactly `segments` lines"
        );
    }

    #[test]
    fn test_debug_draw_arrow_zero_dir_no_lines() {
        let mut dd = DebugDraw::new();
        dd.draw_arrow([0.0; 3], [0.0, 0.0, 0.0], 1.0, [1.0, 0.0, 0.0]);
        assert_eq!(
            dd.line_count(),
            0,
            "zero-direction arrow should produce no lines"
        );
    }

    #[test]
    fn test_debug_draw_multiple_operations_accumulate() {
        let mut dd = DebugDraw::new();
        dd.draw_line([0.0; 3], [1.0; 3], [1.0; 3]); // 1 line
        dd.draw_cross([0.0; 3], 1.0, [0.5; 3]); // 3 lines
        dd.draw_arrow([0.0; 3], [1.0, 0.0, 0.0], 2.0, [1.0, 0.0, 0.0]); // 2 lines
        assert_eq!(dd.line_count(), 6);
    }

    #[test]
    fn test_debug_draw_circle_normal_x_axis() {
        let mut dd = DebugDraw::new();
        let seg = 8;
        dd.draw_circle([0.0; 3], [1.0, 0.0, 0.0], 1.0, seg, [0.0, 1.0, 0.0]);
        assert_eq!(dd.line_count(), seg);
    }

    #[test]
    fn test_debug_draw_circle_zero_normal_skipped() {
        let mut dd = DebugDraw::new();
        dd.draw_circle([0.0; 3], [0.0, 0.0, 0.0], 1.0, 8, [1.0; 3]);
        assert_eq!(dd.line_count(), 0, "zero-normal circle should be skipped");
    }

    // ── AaWireframeParams ─────────────────────────────────────────────────────

    #[test]
    fn test_aa_params_feature_edge_color_yellowish() {
        let p = AaWireframeParams::feature_edge_params();
        assert!(p.color[0] > 0.9, "feature edge should have high red");
        assert!(p.color[2] < 0.1, "feature edge should have low blue");
    }

    #[test]
    fn test_aa_params_silhouette_on_top() {
        let p = AaWireframeParams::silhouette_params();
        assert!(!p.depth_test, "silhouette should be drawn always-on-top");
    }

    #[test]
    fn test_aa_params_alpha_smoothstep_mid_falloff() {
        let p = AaWireframeParams::default_params();
        // At half of falloff past the edge, alpha should be between 0 and 1
        let at_edge = p.line_width_px * 0.5;
        let at_mid_falloff = at_edge + p.falloff_px * 0.5;
        let alpha = p.alpha_at_distance(at_mid_falloff);
        assert!(
            alpha > 0.0 && alpha < 1.0,
            "mid-falloff alpha should be partial"
        );
    }

    #[test]
    fn test_aa_params_alpha_beyond_falloff_zero() {
        let p = AaWireframeParams::default_params();
        let far = p.line_width_px + p.falloff_px + 1.0;
        assert!(
            p.alpha_at_distance(far) < 1e-9,
            "alpha past falloff should be zero"
        );
    }

    // ── Depth fade ────────────────────────────────────────────────────────────

    #[test]
    fn test_depth_fade_at_near_plane_is_one() {
        let cam = [0.0; 3];
        let near = [1.0, 0.0, 0.0]; // exactly at near
        let fade = depth_fade(near, cam, 1.0, 100.0);
        assert!(
            (fade - 1.0).abs() < 1e-9,
            "at near plane fade should be 1.0"
        );
    }

    #[test]
    fn test_depth_fade_clamped_below_zero() {
        let cam = [0.0; 3];
        let beyond_far = [200.0, 0.0, 0.0];
        let fade = depth_fade(beyond_far, cam, 0.0, 100.0);
        assert!(fade >= 0.0, "fade should not go below 0");
        assert!((fade).abs() < 1e-9, "fade beyond far should be 0");
    }

    // ── segment_intersects_aabb edge cases ────────────────────────────────────

    #[test]
    fn test_segment_aabb_touches_corner() {
        // Segment going exactly to a corner
        let p = [-2.0, -2.0, -2.0];
        let q = [0.0, 0.0, 0.0]; // corner of AABB
        let mn = [0.0; 3];
        let mx = [1.0; 3];
        assert!(segment_intersects_aabb(p, q, mn, mx));
    }

    #[test]
    fn test_segment_aabb_both_inside() {
        let p = [0.1, 0.1, 0.1];
        let q = [0.9, 0.9, 0.9];
        let mn = [0.0; 3];
        let mx = [1.0; 3];
        assert!(segment_intersects_aabb(p, q, mn, mx));
    }

    #[test]
    fn test_segment_aabb_parallel_along_y_inside_slab() {
        // Segment along Y, within XZ bounds
        let p = [0.5, -10.0, 0.5];
        let q = [0.5, 10.0, 0.5];
        let mn = [0.0; 3];
        let mx = [1.0; 3];
        assert!(segment_intersects_aabb(p, q, mn, mx));
    }

    // ── cull_wireframe_batch ──────────────────────────────────────────────────

    #[test]
    fn test_cull_batch_partial_visibility() {
        // Two segments: one near origin, one far away
        let mut batch = WireframeBatch::new("partial");
        batch.add([0.0; 3], [0.5; 3], [1.0; 3], 1.0); // inside unit frustum
        batch.add([50.0; 3], [51.0; 3], [1.0; 3], 1.0); // far outside
        let planes: [[f64; 4]; 6] = [
            [1.0, 0.0, 0.0, 5.0],
            [-1.0, 0.0, 0.0, 5.0],
            [0.0, 1.0, 0.0, 5.0],
            [0.0, -1.0, 0.0, 5.0],
            [0.0, 0.0, 1.0, 5.0],
            [0.0, 0.0, -1.0, 5.0],
        ];
        let culled = cull_wireframe_batch(&batch, &planes);
        // Near segment should survive, far one culled
        assert_eq!(culled.len(), 1, "only near segment should survive");
    }

    // ── WireframeColorMap ─────────────────────────────────────────────────────

    #[test]
    fn test_color_map_single_stop_always_that_color() {
        let mut cmap = WireframeColorMap::new();
        cmap.add_stop(0.5, [0.3, 0.6, 0.9]);
        let at_zero = cmap.sample(0.0);
        let at_one = cmap.sample(1.0);
        // Only one stop: should return that color regardless of t
        assert!((at_zero[0] - 0.3).abs() < 1e-9);
        assert!((at_one[0] - 0.3).abs() < 1e-9);
    }

    #[test]
    fn test_color_map_heat_midpoint_is_green() {
        let cmap = WireframeColorMap::heat();
        let mid = cmap.sample(0.5);
        // At t=0.5 (between cyan=0.25 and yellow=0.75), should be green
        assert!(mid[1] > 0.5, "heat map t=0.5 should have significant green");
    }

    #[test]
    fn test_color_map_add_stop_sorted() {
        let mut cmap = WireframeColorMap::new();
        cmap.add_stop(1.0, [1.0, 0.0, 0.0]);
        cmap.add_stop(0.0, [0.0, 0.0, 1.0]);
        // After sort, t=0 → blue, t=1 → red
        let c0 = cmap.sample(0.0);
        let c1 = cmap.sample(1.0);
        assert!(c0[2] > c0[0], "t=0 should be bluer than red");
        assert!(c1[0] > c1[2], "t=1 should be redder than blue");
    }

    // ── WireframeRenderer advanced ────────────────────────────────────────────

    #[test]
    fn test_renderer_push_segment_explicit_color() {
        let mut r = WireframeRenderer::new([1.0; 3], 1.0);
        r.push_segment([0.0; 3], [1.0; 3], [0.2, 0.4, 0.6], 2.0);
        let seg = &r.buffer[0];
        assert!((seg.color[0] - 0.2_f32).abs() < 1e-4);
        assert!((seg.color[1] - 0.4_f32).abs() < 1e-4);
        assert!((seg.color[2] - 0.6_f32).abs() < 1e-4);
    }

    #[test]
    fn test_renderer_push_uses_default_color() {
        let default_color = [0.1, 0.2, 0.3];
        let mut r = WireframeRenderer::new(default_color, 1.0);
        r.push([0.0; 3], [1.0; 3]);
        let seg = &r.buffer[0];
        assert!((seg.color[0] - 0.1_f32).abs() < 1e-4);
        assert!((seg.color[1] - 0.2_f32).abs() < 1e-4);
        assert!((seg.color[2] - 0.3_f32).abs() < 1e-4);
    }

    #[test]
    fn test_renderer_draw_frustum_near_corners_closer_than_far() {
        let mut r = WireframeRenderer::new([1.0; 3], 1.0);
        let hfov = std::f64::consts::FRAC_PI_4;
        r.draw_frustum(1.0, 50.0, hfov, hfov, [1.0; 3]);
        // Near corners: z = -1, far corners: z = -50
        let near_segs: Vec<_> = r
            .buffer
            .iter()
            .filter(|s| s.start[2].abs() < 2.0 && s.end[2].abs() < 2.0)
            .collect();
        assert!(!near_segs.is_empty(), "should have segments near z=-1");
    }

    #[test]
    fn test_renderer_draw_grid_large() {
        let mut r = WireframeRenderer::new([0.5; 3], 1.0);
        let nx = 10;
        let nz = 10;
        r.draw_grid([0.0; 3], nx, nz, 0.5, [0.5; 3]);
        let expected = (nz + 1) + (nx + 1);
        assert_eq!(r.segment_count(), expected);
    }

    #[test]
    fn test_renderer_flush_produces_correct_label() {
        let mut r = WireframeRenderer::new([1.0; 3], 1.0);
        r.push([0.0; 3], [1.0; 3]);
        let batch = r.flush_to_batch("physics_debug");
        assert_eq!(batch.label, "physics_debug");
    }

    #[test]
    fn test_renderer_byte_size_after_multiple_pushes() {
        let mut r = WireframeRenderer::new([1.0; 3], 1.0);
        let n = 5;
        for _ in 0..n {
            r.push([0.0; 3], [1.0; 3]);
        }
        let expected = n * std::mem::size_of::<WireSegmentGpu>();
        assert_eq!(r.byte_size(), expected);
    }

    // ── contact_normal_lines ──────────────────────────────────────────────────

    #[test]
    fn test_contact_normal_lines_count() {
        use crate::primitives::Color;
        use oxiphysics_core::math::Vec3;
        let contacts = vec![
            (Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0), 0.1),
            (Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0), 0.2),
            (Vec3::new(2.0, 0.0, 0.0), Vec3::new(-1.0, 0.0, 0.0), 0.3),
        ];
        let color = Color::new(1.0, 0.0, 0.0, 1.0);
        let lines = contact_normal_lines(&contacts, color);
        assert_eq!(lines.len(), 3);
    }

    #[test]
    fn test_contact_normal_lines_empty() {
        let lines = contact_normal_lines(&[], Color::new(1.0, 0.0, 0.0, 1.0));
        assert!(lines.is_empty());
    }

    // ── build_aa_feature_edge_batch ───────────────────────────────────────────

    #[test]
    fn test_build_aa_feature_edge_batch_count() {
        let pos = vec![[0.0, 0.0, 0.0_f64], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]];
        let idx = vec![0usize, 1, 2];
        let edges = extract_mesh_edges(&pos, &idx);
        let params = AaWireframeParams::feature_edge_params();
        let batch = build_aa_feature_edge_batch(&pos, &edges, &params);
        assert_eq!(batch.len(), edges.len());
    }

    // ── WireframeEdge additional ──────────────────────────────────────────────

    #[test]
    fn test_wireframe_edge_new_stores_fields() {
        let e = WireframeEdge::new(3, 7, [0.5, 0.5, 0.5], 2.0);
        assert_eq!(e.v0, 3);
        assert_eq!(e.v1, 7);
        assert!((e.line_width - 2.0).abs() < 1e-6);
    }

    #[test]
    fn test_wireframe_edge_canonical_key_both_orderings() {
        let e1 = WireframeEdge::new(0, 5, [1.0; 3], 1.0);
        let e2 = WireframeEdge::new(5, 0, [1.0; 3], 1.0);
        assert_eq!(e1.canonical_key(), e2.canonical_key());
        assert_eq!(e1.canonical_key(), (0, 5));
    }

    // ── RichWireframeMesh empty ───────────────────────────────────────────────

    #[test]
    fn test_rich_wireframe_mesh_empty_counts() {
        let m = RichWireframeMesh::empty();
        assert_eq!(m.vertex_count(), 0);
        assert_eq!(m.edge_count(), 0);
        assert!(m.vertices.is_empty());
        assert!(m.normals.is_empty());
        assert!(m.edges.is_empty());
    }

    // ── WireframeBuilder advanced ─────────────────────────────────────────────

    #[test]
    fn test_wireframe_builder_colored_edge_stored() {
        let mut b = WireframeBuilder::new([1.0, 1.0, 1.0], 1.0);
        let v0 = b.add_vertex([0.0; 3], [0.0, 1.0, 0.0]);
        let v1 = b.add_vertex([1.0; 3], [0.0, 1.0, 0.0]);
        b.add_edge_colored(v0, v1, [1.0, 0.0, 0.0], 2.5);
        let mesh = b.build();
        assert_eq!(mesh.edge_count(), 1);
        assert!((mesh.edges[0].line_width - 2.5).abs() < 1e-6);
        assert_eq!(mesh.edges[0].color, [1.0, 0.0, 0.0]);
    }

    #[test]
    fn test_wireframe_builder_empty_build() {
        let b = WireframeBuilder::new([1.0; 3], 1.0);
        let mesh = b.build();
        assert_eq!(mesh.vertex_count(), 0);
        assert_eq!(mesh.edge_count(), 0);
    }

    // ── frustum_wireframe_from_fov consistency ────────────────────────────────

    #[test]
    fn test_frustum_wireframe_from_fov_wider_fov_larger_far_plane() {
        let narrow = frustum_wireframe_from_fov(30.0, 1.0, 1.0, 10.0, [1.0; 3], 1.0);
        let wide = frustum_wireframe_from_fov(90.0, 1.0, 1.0, 10.0, [1.0; 3], 1.0);
        // Far corners should be farther from Z axis with wider FOV
        let narrow_far_x = narrow.vertices[4][0].abs();
        let wide_far_x = wide.vertices[4][0].abs();
        assert!(
            wide_far_x > narrow_far_x,
            "wider fov should have larger far plane X extent"
        );
    }

    #[test]
    fn test_frustum_wireframe_edge_count_is_12() {
        use std::f64::consts::FRAC_PI_6;
        let mesh = frustum_wireframe(0.5, 20.0, FRAC_PI_6, FRAC_PI_6, [0.0, 1.0, 1.0], 1.5);
        assert_eq!(mesh.edge_count(), 12);
        assert_eq!(mesh.vertex_count(), 8);
    }

    // ── AABB wireframe_rich ───────────────────────────────────────────────────

    #[test]
    fn test_aabb_wireframe_rich_correct_corner_count() {
        // A box always has exactly 8 corners
        let mesh = aabb_wireframe_rich([-3.0; 3], [3.0; 3], [1.0, 0.5, 0.0], 1.0);
        assert_eq!(mesh.vertex_count(), 8);
    }

    #[test]
    fn test_aabb_wireframe_rich_nonsymmetric_box() {
        let mesh = aabb_wireframe_rich([0.0, 0.0, 0.0], [10.0, 2.0, 5.0], [1.0; 3], 1.0);
        assert_eq!(mesh.edge_count(), 12);
        // All edge vertices must be in bounds
        for e in &mesh.edges {
            assert!(e.v0 < 8 && e.v1 < 8);
        }
    }

    // ── detect_feature_edges threshold sensitivity ────────────────────────────

    #[test]
    fn test_detect_feature_edges_high_threshold_includes_fewer() {
        // Cube: dihedral = 90°. High threshold (85°) should catch all feature edges.
        // Very high threshold (91°) should catch none of the 90° edges.
        let pos = vec![
            [0.0, 0.0, 0.0_f64],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
            [0.0, 1.0, 1.0],
        ];
        let idx = vec![
            0usize, 1, 2, 0, 2, 3, 4, 6, 5, 4, 7, 6, 0, 4, 5, 0, 5, 1, 3, 2, 6, 3, 6, 7, 0, 3, 7,
            0, 7, 4, 1, 5, 6, 1, 6, 2,
        ];
        let feat_85 = detect_feature_edges(&pos, &idx, 85.0);
        let feat_91 = detect_feature_edges(&pos, &idx, 91.0);
        // At 85° threshold, all 90° cube edges should be features
        // At 91° threshold, none of the 90° edges (only > 90° edges, but cube has none)
        assert!(feat_85.len() >= feat_91.len());
    }
} // end mod tests

// ─────────────────────────────────────────────────────────────────────────────
// § WireframeSegment — annotated line segment with metadata
// ─────────────────────────────────────────────────────────────────────────────

/// A line segment with colour, thickness, and an optional label.
#[derive(Debug, Clone)]
pub struct WireframeSegment {
    /// Start position.
    pub start: [f64; 3],
    /// End position.
    pub end: [f64; 3],
    /// RGBA colour (0–1 range).
    pub color: [f32; 4],
    /// Line width in screen pixels.
    pub thickness: f32,
    /// Optional annotation text.
    pub label: Option<String>,
}

impl WireframeSegment {
    /// Create a plain white unit-thickness segment.
    pub fn new(start: [f64; 3], end: [f64; 3]) -> Self {
        WireframeSegment {
            start,
            end,
            color: [1.0, 1.0, 1.0, 1.0],
            thickness: 1.0,
            label: None,
        }
    }

    /// Set the colour.
    pub fn with_color(mut self, r: f32, g: f32, b: f32, a: f32) -> Self {
        self.color = [r, g, b, a];
        self
    }

    /// Set the thickness.
    pub fn with_thickness(mut self, t: f32) -> Self {
        self.thickness = t;
        self
    }

    /// Set an annotation label.
    pub fn with_label(mut self, label: &str) -> Self {
        self.label = Some(label.to_string());
        self
    }

    /// Euclidean length of the segment.
    pub fn length(&self) -> f64 {
        let dx = self.end[0] - self.start[0];
        let dy = self.end[1] - self.start[1];
        let dz = self.end[2] - self.start[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }

    /// Midpoint of the segment.
    pub fn midpoint(&self) -> [f64; 3] {
        [
            (self.start[0] + self.end[0]) * 0.5,
            (self.start[1] + self.end[1]) * 0.5,
            (self.start[2] + self.end[2]) * 0.5,
        ]
    }

    /// Direction vector (not normalised).
    pub fn direction(&self) -> [f64; 3] {
        [
            self.end[0] - self.start[0],
            self.end[1] - self.start[1],
            self.end[2] - self.start[2],
        ]
    }

    /// Unit direction vector.  Returns `[0,0,0]` for zero-length segments.
    pub fn unit_direction(&self) -> [f64; 3] {
        let d = self.direction();
        let len = self.length();
        if len < 1e-300 {
            return [0.0; 3];
        }
        [d[0] / len, d[1] / len, d[2] / len]
    }

    /// Test if this segment intersects the axis-aligned slab `z_min ≤ z ≤ z_max`.
    pub fn intersects_z_slab(&self, z_min: f64, z_max: f64) -> bool {
        let z0 = self.start[2].min(self.end[2]);
        let z1 = self.start[2].max(self.end[2]);
        z1 >= z_min && z0 <= z_max
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// § WireframeGrid — grid of lines for debugging / scene floor
// ─────────────────────────────────────────────────────────────────────────────

/// A flat XZ-plane grid for use as a scene floor.
#[derive(Debug, Clone)]
pub struct WireframeGrid {
    /// All line segments composing the grid.
    pub segments: Vec<WireframeSegment>,
    /// Grid spacing.
    pub spacing: f64,
    /// Number of cells in each direction.
    pub n_cells: usize,
}

/// Generate a flat grid on the XZ-plane centred at the origin.
///
/// `n_cells` cells in each direction, `spacing` between grid lines.
pub fn generate_floor_grid(n_cells: usize, spacing: f64) -> WireframeGrid {
    let half = n_cells as f64 * spacing * 0.5;
    let mut segments = Vec::new();

    let n_lines = n_cells + 1;
    for i in 0..n_lines {
        let t = -half + i as f64 * spacing;
        // Lines parallel to Z
        segments.push(WireframeSegment::new([t, 0.0, -half], [t, 0.0, half]));
        // Lines parallel to X
        segments.push(WireframeSegment::new([-half, 0.0, t], [half, 0.0, t]));
    }

    WireframeGrid {
        segments,
        spacing,
        n_cells,
    }
}

impl WireframeGrid {
    /// Total number of line segments.
    pub fn segment_count(&self) -> usize {
        self.segments.len()
    }

    /// Total length of all segments combined.
    pub fn total_length(&self) -> f64 {
        self.segments.iter().map(|s| s.length()).sum()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// § DepthSortBuffer — back-to-front sorting for transparent wireframes
// ─────────────────────────────────────────────────────────────────────────────

/// A buffer for depth-sorted wireframe segments (painter's algorithm).
#[derive(Debug, Clone, Default)]
pub struct DepthSortBuffer {
    /// Stored segments with their precomputed depth keys.
    entries: Vec<(f64, WireframeSegment)>,
}

impl DepthSortBuffer {
    /// Create an empty buffer.
    pub fn new() -> Self {
        DepthSortBuffer {
            entries: Vec::new(),
        }
    }

    /// Insert a segment.  Depth is computed as the midpoint Z coordinate (for
    /// a view along +Z).  Use a custom `depth_fn` variant for other views.
    pub fn insert(&mut self, seg: WireframeSegment) {
        let depth = seg.midpoint()[2];
        self.entries.push((depth, seg));
    }

    /// Insert with an explicit depth key.
    pub fn insert_with_depth(&mut self, seg: WireframeSegment, depth: f64) {
        self.entries.push((depth, seg));
    }

    /// Sort in back-to-front order (decreasing depth).
    pub fn sort_back_to_front(&mut self) {
        self.entries
            .sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));
    }

    /// Sort in front-to-back order (increasing depth).
    pub fn sort_front_to_back(&mut self) {
        self.entries
            .sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    }

    /// Return segments in current order.
    pub fn segments(&self) -> Vec<&WireframeSegment> {
        self.entries.iter().map(|(_, s)| s).collect()
    }

    /// Number of segments.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// True if empty.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Clear all entries.
    pub fn clear(&mut self) {
        self.entries.clear();
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// § Coordinate-axis arrows
// ─────────────────────────────────────────────────────────────────────────────

/// Generate RGB axis arrows at the origin with length `size`.
///
/// Returns three segments: X (red), Y (green), Z (blue).
pub fn axis_arrows(size: f64) -> [WireframeSegment; 3] {
    [
        WireframeSegment::new([0.0, 0.0, 0.0], [size, 0.0, 0.0])
            .with_color(1.0, 0.0, 0.0, 1.0)
            .with_label("X"),
        WireframeSegment::new([0.0, 0.0, 0.0], [0.0, size, 0.0])
            .with_color(0.0, 1.0, 0.0, 1.0)
            .with_label("Y"),
        WireframeSegment::new([0.0, 0.0, 0.0], [0.0, 0.0, size])
            .with_color(0.0, 0.0, 1.0, 1.0)
            .with_label("Z"),
    ]
}

/// Generate a wireframe circle (approximate polygon) on the XY-plane.
pub fn circle_wireframe(centre: [f64; 3], radius: f64, n_segments: usize) -> Vec<WireframeSegment> {
    let n = n_segments.max(3);
    let mut segs = Vec::with_capacity(n);
    for i in 0..n {
        let theta0 = 2.0 * std::f64::consts::PI * i as f64 / n as f64;
        let theta1 = 2.0 * std::f64::consts::PI * (i + 1) as f64 / n as f64;
        let p0 = [
            centre[0] + radius * theta0.cos(),
            centre[1] + radius * theta0.sin(),
            centre[2],
        ];
        let p1 = [
            centre[0] + radius * theta1.cos(),
            centre[1] + radius * theta1.sin(),
            centre[2],
        ];
        segs.push(WireframeSegment::new(p0, p1));
    }
    segs
}

/// Generate a wireframe capsule (cylinder + hemisphere caps) along the Z-axis.
///
/// `radius` is the capsule radius, `half_height` is the half-length of the
/// cylindrical section.  Returns segments for the outline only.
pub fn capsule_outline_xz(
    centre: [f64; 3],
    radius: f64,
    half_height: f64,
    n_segs: usize,
) -> Vec<WireframeSegment> {
    let n = n_segs.max(3);
    let mut segs = Vec::new();
    let cy = centre[1];
    let cz = centre[2];

    // Top circle
    let top_centre = [centre[0], cy + half_height, cz];
    segs.extend(circle_wireframe(top_centre, radius, n));

    // Bottom circle
    let bot_centre = [centre[0], cy - half_height, cz];
    segs.extend(circle_wireframe(bot_centre, radius, n));

    // Four vertical lines
    for &dx in &[-1.0_f64, 1.0] {
        segs.push(WireframeSegment::new(
            [centre[0] + dx * radius, cy + half_height, cz],
            [centre[0] + dx * radius, cy - half_height, cz],
        ));
    }
    for &dz in &[-1.0_f64, 1.0] {
        segs.push(WireframeSegment::new(
            [centre[0], cy + half_height, cz + dz * radius],
            [centre[0], cy - half_height, cz + dz * radius],
        ));
    }
    segs
}

#[cfg(test)]
mod tests_wireframe_ext {
    use super::*;

    // ── WireframeSegment ──────────────────────────────────────────────────────

    #[test]
    fn segment_length_unit() {
        let s = WireframeSegment::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert!((s.length() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn segment_length_diagonal() {
        let s = WireframeSegment::new([0.0; 3], [1.0, 1.0, 1.0]);
        assert!((s.length() - 3.0_f64.sqrt()).abs() < 1e-10);
    }

    #[test]
    fn segment_midpoint_correct() {
        let s = WireframeSegment::new([0.0, 0.0, 0.0], [2.0, 4.0, 6.0]);
        let m = s.midpoint();
        assert_eq!(m, [1.0, 2.0, 3.0]);
    }

    #[test]
    fn segment_unit_direction_x() {
        let s = WireframeSegment::new([0.0; 3], [3.0, 0.0, 0.0]);
        let d = s.unit_direction();
        assert!((d[0] - 1.0).abs() < 1e-12);
        assert!(d[1].abs() < 1e-12);
        assert!(d[2].abs() < 1e-12);
    }

    #[test]
    fn segment_unit_direction_zero_len() {
        let s = WireframeSegment::new([1.0; 3], [1.0; 3]);
        assert_eq!(s.unit_direction(), [0.0; 3]);
    }

    #[test]
    fn segment_z_slab_intersect_true() {
        let s = WireframeSegment::new([0.0, 0.0, -1.0], [0.0, 0.0, 1.0]);
        assert!(s.intersects_z_slab(-0.5, 0.5));
    }

    #[test]
    fn segment_z_slab_intersect_false() {
        let s = WireframeSegment::new([0.0, 0.0, 5.0], [0.0, 0.0, 10.0]);
        assert!(!s.intersects_z_slab(0.0, 4.0));
    }

    #[test]
    fn segment_with_label_stored() {
        let s = WireframeSegment::new([0.0; 3], [1.0; 3]).with_label("edge");
        assert_eq!(s.label.as_deref(), Some("edge"));
    }

    #[test]
    fn segment_with_thickness_stored() {
        let s = WireframeSegment::new([0.0; 3], [1.0; 3]).with_thickness(3.5);
        assert!((s.thickness - 3.5).abs() < 1e-6);
    }

    // ── generate_floor_grid ───────────────────────────────────────────────────

    #[test]
    fn floor_grid_segment_count_2x2() {
        let g = generate_floor_grid(2, 1.0);
        // 2 cells → 3 lines per direction, 2 directions → 6 segs
        assert_eq!(g.segment_count(), 6);
    }

    #[test]
    fn floor_grid_total_length_positive() {
        let g = generate_floor_grid(4, 0.5);
        assert!(g.total_length() > 0.0);
    }

    #[test]
    fn floor_grid_n_cells_stored() {
        let g = generate_floor_grid(10, 1.0);
        assert_eq!(g.n_cells, 10);
    }

    // ── DepthSortBuffer ───────────────────────────────────────────────────────

    #[test]
    fn depth_buffer_insert_and_count() {
        let mut buf = DepthSortBuffer::new();
        buf.insert(WireframeSegment::new([0.0; 3], [1.0; 3]));
        buf.insert(WireframeSegment::new([2.0; 3], [3.0; 3]));
        assert_eq!(buf.len(), 2);
    }

    #[test]
    fn depth_buffer_sort_back_to_front_order() {
        let mut buf = DepthSortBuffer::new();
        buf.insert_with_depth(WireframeSegment::new([0.0; 3], [1.0; 3]), 1.0);
        buf.insert_with_depth(WireframeSegment::new([0.0; 3], [1.0; 3]), 5.0);
        buf.insert_with_depth(WireframeSegment::new([0.0; 3], [1.0; 3]), 3.0);
        buf.sort_back_to_front();
        let depths: Vec<f64> = buf.entries.iter().map(|(d, _)| *d).collect();
        assert!(depths[0] >= depths[1]);
        assert!(depths[1] >= depths[2]);
    }

    #[test]
    fn depth_buffer_sort_front_to_back_order() {
        let mut buf = DepthSortBuffer::new();
        buf.insert_with_depth(WireframeSegment::new([0.0; 3], [1.0; 3]), 5.0);
        buf.insert_with_depth(WireframeSegment::new([0.0; 3], [1.0; 3]), 1.0);
        buf.sort_front_to_back();
        let depths: Vec<f64> = buf.entries.iter().map(|(d, _)| *d).collect();
        assert!(depths[0] <= depths[1]);
    }

    #[test]
    fn depth_buffer_clear() {
        let mut buf = DepthSortBuffer::new();
        buf.insert(WireframeSegment::new([0.0; 3], [1.0; 3]));
        buf.clear();
        assert!(buf.is_empty());
    }

    // ── axis_arrows ───────────────────────────────────────────────────────────

    #[test]
    fn axis_arrows_lengths_equal_size() {
        let arrows = axis_arrows(2.0);
        for a in &arrows {
            assert!((a.length() - 2.0).abs() < 1e-12);
        }
    }

    #[test]
    fn axis_arrows_labels_xyz() {
        let arrows = axis_arrows(1.0);
        let labels: Vec<&str> = arrows
            .iter()
            .map(|a| a.label.as_deref().unwrap_or(""))
            .collect();
        assert!(labels.contains(&"X"));
        assert!(labels.contains(&"Y"));
        assert!(labels.contains(&"Z"));
    }

    // ── circle_wireframe ──────────────────────────────────────────────────────

    #[test]
    fn circle_wireframe_segment_count() {
        let segs = circle_wireframe([0.0; 3], 1.0, 16);
        assert_eq!(segs.len(), 16);
    }

    #[test]
    fn circle_wireframe_segments_close_loop() {
        let segs = circle_wireframe([0.0; 3], 1.0, 8);
        // First start and last end should be equal
        let first_start = segs[0].start;
        let last_end = segs[segs.len() - 1].end;
        for k in 0..3 {
            assert!((first_start[k] - last_end[k]).abs() < 1e-12);
        }
    }

    #[test]
    fn circle_wireframe_radius_respected() {
        let r = 3.5;
        let segs = circle_wireframe([0.0; 3], r, 32);
        for s in &segs {
            let d = (s.start[0] * s.start[0] + s.start[1] * s.start[1]).sqrt();
            assert!((d - r).abs() < 1e-10, "radius mismatch: {d} vs {r}");
        }
    }

    // ── capsule_outline_xz ────────────────────────────────────────────────────

    #[test]
    fn capsule_outline_non_empty() {
        let segs = capsule_outline_xz([0.0; 3], 0.5, 1.0, 8);
        assert!(!segs.is_empty());
    }

    #[test]
    fn capsule_outline_has_vertical_lines() {
        let segs = capsule_outline_xz([0.0; 3], 0.5, 1.0, 8);
        // Some segments should span y from -1 to +1
        let has_vertical = segs.iter().any(|s| {
            let dy = (s.end[1] - s.start[1]).abs();
            dy > 1.5
        });
        assert!(has_vertical);
    }
}
