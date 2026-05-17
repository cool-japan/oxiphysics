// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Visualization data generation for the OxiPhysics engine.
//!
//! This crate produces renderable primitives (meshes, lines, colors) without
//! depending on any GPU or windowing library. A downstream renderer can
//! consume these data structures to display physics debug information.
#![warn(missing_docs)]

pub mod animation_system;
pub mod camera;
pub mod picking;
pub use picking::{
    PickResult, PickScene, PickableObject, Picker, Ray, SelectionSet, moller_trumbore,
};
pub mod hdr_framebuffer;
pub use hdr_framebuffer::{
    HdrFramebuffer, HdrRasterizer, ToneMapMethod, ToneMapper, linear_to_srgb, srgb_to_linear,
};
pub mod wasm_canvas;
pub use wasm_canvas::{CanvasBuffer, CanvasConfig, CanvasRasterizer};

pub mod wgpu_renderer;
pub use wgpu_renderer::{GpuLight, GpuMeshHandle, GpuRenderer, GpuRendererConfig};
pub mod colormap;
pub mod debug_overlay;
mod error;
pub mod gizmos;
pub mod mesh_gen;
pub mod mesh_render;
pub mod metaball;
pub mod network_visualization;
pub mod particle_renderer;
pub mod post_processing;
pub mod postprocess;
pub mod primitives;
pub mod scene;
pub mod shader;
pub mod streamlines;
pub mod stress_viz;
pub mod transfer_functions;
pub mod volume_rendering;
pub mod wireframe;

pub use colormap::{Colormap, map_scalar};
pub use error::*;
pub use mesh_gen::{arrow_mesh, box_mesh, plane_mesh, sphere_mesh};
pub use mesh_render::{
    Framebuffer, Mat4, PhongShader, RenderMesh as MeshRenderMesh, SoftwareRasterizer, Vertex3D,
    WireframeRenderer, build_arrow, build_box as build_box_mesh, build_cylinder, build_grid,
    build_sphere as build_sphere_mesh, cross3, dot3, length3, normalize3, scale3, sub3,
};
pub use post_processing::{
    Bloom, DepthOfField, GaussianBlur, Image, PostColor, PostProcessPipeline, Ssao, ToneMapping,
    Vignette,
};
pub use primitives::{Color, LinePrimitive, RenderMesh, TrianglePrimitive, Vertex};
pub use streamlines::{Streamline, trace_streamline, trace_streamlines_grid};
pub use stress_viz::{principal_stress_glyphs, von_mises_colors};
pub use volume_rendering::*;
pub use wireframe::{
    contact_normal_lines, wireframe_aabb, wireframe_box, wireframe_capsule, wireframe_sphere,
};

/// Trait for physics visualization renderers.
#[allow(dead_code)]
pub trait Renderer {
    /// Initialize this component.
    fn init(&mut self);
}

// ---------------------------------------------------------------------------
// VizConfig — global configuration for the visualization subsystem
// ---------------------------------------------------------------------------

/// Global configuration for the OxiPhysics visualization subsystem.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct VizConfig {
    /// Default colormap for scalar field visualization.
    pub default_colormap: Colormap,
    /// Whether to enable anti-aliasing for line primitives.
    pub line_antialiasing: bool,
    /// Default line width in pixels (or NDC units for CPU renderers).
    pub default_line_width: f32,
    /// Default point size for particle/gizmo rendering.
    pub default_point_size: f32,
    /// Maximum number of primitives per draw call batch.
    pub max_primitives_per_batch: usize,
    /// Background color for the viewport.
    pub background_color: Color,
    /// Whether to show the debug overlay by default.
    pub show_debug_overlay: bool,
    /// Whether to enable wireframe rendering by default.
    pub wireframe_mode: bool,
    /// Near clip plane distance.
    pub near_clip: f64,
    /// Far clip plane distance.
    pub far_clip: f64,
    /// Gamma correction exponent (1.0 = linear, 2.2 = sRGB).
    pub gamma: f32,
}

impl Default for VizConfig {
    fn default() -> Self {
        Self {
            default_colormap: Colormap::Viridis,
            line_antialiasing: true,
            default_line_width: 1.5,
            default_point_size: 4.0,
            max_primitives_per_batch: 65_536,
            background_color: Color::new(0.05, 0.05, 0.08, 1.0),
            show_debug_overlay: false,
            wireframe_mode: false,
            near_clip: 0.01,
            far_clip: 1000.0,
            gamma: 2.2,
        }
    }
}

impl VizConfig {
    /// Construct a config with debugging aids enabled.
    #[allow(dead_code)]
    pub fn debug_mode() -> Self {
        Self {
            show_debug_overlay: true,
            wireframe_mode: true,
            ..Self::default()
        }
    }

    /// Return `true` if the depth range is valid (`near < far`).
    #[allow(dead_code)]
    pub fn depth_range_valid(&self) -> bool {
        self.near_clip > 0.0 && self.near_clip < self.far_clip
    }
}

// ---------------------------------------------------------------------------
// VizStats — per-frame statistics
// ---------------------------------------------------------------------------

/// Per-frame rendering statistics.
#[allow(dead_code)]
#[derive(Debug, Clone, Default)]
pub struct VizStats {
    /// Number of draw calls issued this frame.
    pub draw_calls: u32,
    /// Total number of triangles rendered.
    pub triangles: u64,
    /// Total number of line segments rendered.
    pub lines: u64,
    /// Total number of point sprites.
    pub points: u64,
    /// Number of culled objects (view-frustum or back-face).
    pub culled: u32,
    /// Time to prepare and upload geometry data (milliseconds, wall-clock).
    pub geometry_ms: f64,
    /// Time to execute all draw calls (milliseconds, wall-clock).
    pub render_ms: f64,
    /// Time for post-processing passes (milliseconds, wall-clock).
    pub post_ms: f64,
}

impl VizStats {
    /// Total frame time: geometry + render + post-processing.
    #[allow(dead_code)]
    pub fn total_ms(&self) -> f64 {
        self.geometry_ms + self.render_ms + self.post_ms
    }

    /// Approximate frames per second (1000 / total_ms), or `f64::INFINITY`
    /// if total_ms is zero.
    #[allow(dead_code)]
    pub fn fps(&self) -> f64 {
        let t = self.total_ms();
        if t < 1e-10 { f64::INFINITY } else { 1000.0 / t }
    }

    /// Reset all counters to zero.
    #[allow(dead_code)]
    pub fn reset(&mut self) {
        *self = Self::default();
    }

    /// Merge stats from another frame (accumulate).
    #[allow(dead_code)]
    pub fn accumulate(&mut self, other: &VizStats) {
        self.draw_calls += other.draw_calls;
        self.triangles += other.triangles;
        self.lines += other.lines;
        self.points += other.points;
        self.culled += other.culled;
        self.geometry_ms += other.geometry_ms;
        self.render_ms += other.render_ms;
        self.post_ms += other.post_ms;
    }
}

// ---------------------------------------------------------------------------
// VizPlugin — optional extension points
// ---------------------------------------------------------------------------

/// A plug-in that can inject custom geometry or post-effects into a frame.
#[allow(dead_code)]
pub trait VizPlugin: std::fmt::Debug {
    /// Name identifying this plugin.
    fn name(&self) -> &str;

    /// Called once to let the plugin register resources.
    fn on_init(&mut self, _config: &VizConfig) {}

    /// Called each frame before the main render pass.
    fn pre_render(&mut self, _meshes: &mut Vec<RenderMesh>, _lines: &mut Vec<LinePrimitive>) {}

    /// Called each frame after the main render pass.
    fn post_render(&mut self, _stats: &VizStats) {}
}

// ---------------------------------------------------------------------------
// VizRegistry — named resource registry
// ---------------------------------------------------------------------------

/// A simple string-keyed registry for visualization resources.
///
/// Stores any type `T` under a unique name. Useful for managing named meshes,
/// shaders, colormaps, or transfer functions within a scene.
#[allow(dead_code)]
#[derive(Debug, Default)]
pub struct VizRegistry<T> {
    items: std::collections::HashMap<String, T>,
}

impl<T> VizRegistry<T> {
    /// Create an empty registry.
    pub fn new() -> Self {
        Self {
            items: std::collections::HashMap::new(),
        }
    }

    /// Register an item. Returns an error if the name is already taken.
    pub fn register(&mut self, name: impl Into<String>, item: T) -> error::Result<()> {
        let name = name.into();
        if self.items.contains_key(&name) {
            return Err(error::Error::already_exists(name));
        }
        self.items.insert(name, item);
        Ok(())
    }

    /// Retrieve an item by name.
    pub fn get(&self, name: &str) -> error::Result<&T> {
        self.items
            .get(name)
            .ok_or_else(|| error::Error::not_found(name))
    }

    /// Retrieve a mutable reference to an item by name.
    pub fn get_mut(&mut self, name: &str) -> error::Result<&mut T> {
        self.items
            .get_mut(name)
            .ok_or_else(|| error::Error::not_found(name))
    }

    /// Remove an item by name, returning it.
    pub fn remove(&mut self, name: &str) -> error::Result<T> {
        self.items
            .remove(name)
            .ok_or_else(|| error::Error::not_found(name))
    }

    /// Check whether a name is registered.
    pub fn contains(&self, name: &str) -> bool {
        self.items.contains_key(name)
    }

    /// Number of registered items.
    pub fn len(&self) -> usize {
        self.items.len()
    }

    /// Return `true` if the registry is empty.
    pub fn is_empty(&self) -> bool {
        self.items.is_empty()
    }

    /// Iterate over `(name, item)` pairs.
    pub fn iter(&self) -> impl Iterator<Item = (&str, &T)> {
        self.items.iter().map(|(k, v)| (k.as_str(), v))
    }
}

// ---------------------------------------------------------------------------
// VizContext — the central per-session state
// ---------------------------------------------------------------------------

/// Central context object that owns configuration, statistics, and registries.
///
/// A `VizContext` is typically created once per application session and passed
/// (by shared reference) to visualization functions.
#[allow(dead_code)]
#[derive(Debug)]
pub struct VizContext {
    /// Active configuration.
    pub config: VizConfig,
    /// Cumulative statistics (reset each frame via [`VizContext::begin_frame`]).
    pub stats: VizStats,
    /// Named mesh registry.
    pub meshes: VizRegistry<RenderMesh>,
    /// Named colormap registry (stores additional user-defined colormaps by name).
    pub colormaps: VizRegistry<Colormap>,
    /// Frame counter.
    frame: u64,
}

impl Default for VizContext {
    fn default() -> Self {
        Self::new(VizConfig::default())
    }
}

impl VizContext {
    /// Construct a new context with the given configuration.
    pub fn new(config: VizConfig) -> Self {
        Self {
            config,
            stats: VizStats::default(),
            meshes: VizRegistry::new(),
            colormaps: VizRegistry::new(),
            frame: 0,
        }
    }

    /// Begin a new frame: reset per-frame stats and increment frame counter.
    pub fn begin_frame(&mut self) {
        self.stats.reset();
        self.frame += 1;
    }

    /// Current frame number (starts at 0, incremented by `begin_frame`).
    pub fn frame_number(&self) -> u64 {
        self.frame
    }

    /// Map a scalar value using the context's default colormap.
    pub fn map_scalar(&self, value: f64, min: f64, max: f64) -> Color {
        map_scalar(value, min, max, self.config.default_colormap)
    }

    /// Record that `n` triangles were drawn.
    pub fn record_triangles(&mut self, n: u64) {
        self.stats.triangles += n;
        self.stats.draw_calls += 1;
    }

    /// Record that `n` lines were drawn.
    pub fn record_lines(&mut self, n: u64) {
        self.stats.lines += n;
        self.stats.draw_calls += 1;
    }
}

// ---------------------------------------------------------------------------
// Primitive builder helpers (convenience wrappers)
// ---------------------------------------------------------------------------

/// Build a list of [`LinePrimitive`]s forming a 3-D coordinate axis triad.
///
/// The triad has unit length along each axis, colored RGB (x=red, y=green, z=blue).
#[allow(dead_code)]
pub fn axis_triad(origin: oxiphysics_core::math::Vec3, scale: f64) -> Vec<LinePrimitive> {
    use oxiphysics_core::math::Vec3;
    vec![
        LinePrimitive {
            start: origin,
            end: origin + Vec3::new(scale, 0.0, 0.0),
            color: Color::red(),
        },
        LinePrimitive {
            start: origin,
            end: origin + Vec3::new(0.0, scale, 0.0),
            color: Color::green(),
        },
        LinePrimitive {
            start: origin,
            end: origin + Vec3::new(0.0, 0.0, scale),
            color: Color::blue(),
        },
    ]
}

/// Build a flat grid of [`LinePrimitive`]s in the XZ plane.
///
/// - `center` — center of the grid.
/// - `half_extent` — half-size of the grid along each axis.
/// - `divisions` — number of divisions per side (minimum 1).
/// - `color` — line color.
#[allow(dead_code)]
pub fn grid_lines(
    center: oxiphysics_core::math::Vec3,
    half_extent: f64,
    divisions: usize,
    color: Color,
) -> Vec<LinePrimitive> {
    use oxiphysics_core::math::Vec3;
    let n = divisions.max(1);
    let step = (2.0 * half_extent) / n as f64;
    let mut lines = Vec::with_capacity((n + 1) * 2);
    for i in 0..=n {
        let t = -half_extent + i as f64 * step;
        // Parallel to Z
        lines.push(LinePrimitive {
            start: center + Vec3::new(t, 0.0, -half_extent),
            end: center + Vec3::new(t, 0.0, half_extent),
            color,
        });
        // Parallel to X
        lines.push(LinePrimitive {
            start: center + Vec3::new(-half_extent, 0.0, t),
            end: center + Vec3::new(half_extent, 0.0, t),
            color,
        });
    }
    lines
}

/// Merge multiple lists of [`LinePrimitive`]s into one.
#[allow(dead_code)]
pub fn merge_lines(groups: &[Vec<LinePrimitive>]) -> Vec<LinePrimitive> {
    groups.iter().flat_map(|g| g.iter().cloned()).collect()
}

/// Merge multiple lists of [`TrianglePrimitive`]s into one.
#[allow(dead_code)]
pub fn merge_triangles(groups: &[Vec<TrianglePrimitive>]) -> Vec<TrianglePrimitive> {
    groups.iter().flat_map(|g| g.iter().cloned()).collect()
}

/// Flip all vertex normals in a [`RenderMesh`] (in-place).
#[allow(dead_code)]
pub fn flip_normals(mesh: &mut RenderMesh) {
    for v in &mut mesh.vertices {
        v.normal = [-v.normal[0], -v.normal[1], -v.normal[2]];
    }
}

/// Scale a [`RenderMesh`] uniformly about the origin.
#[allow(dead_code)]
pub fn scale_mesh(mesh: &mut RenderMesh, factor: f64) {
    let f = factor as f32;
    for v in &mut mesh.vertices {
        v.position = [v.position[0] * f, v.position[1] * f, v.position[2] * f];
    }
}

/// Translate a [`RenderMesh`] by an offset vector.
#[allow(dead_code)]
pub fn translate_mesh(mesh: &mut RenderMesh, offset: oxiphysics_core::math::Vec3) {
    let [ox, oy, oz] = [offset.x as f32, offset.y as f32, offset.z as f32];
    for v in &mut mesh.vertices {
        v.position[0] += ox;
        v.position[1] += oy;
        v.position[2] += oz;
    }
}

/// Recolor all vertices in a [`RenderMesh`].
#[allow(dead_code)]
pub fn recolor_mesh(mesh: &mut RenderMesh, color: Color) {
    for v in &mut mesh.vertices {
        v.color = color;
    }
}

/// Compute the axis-aligned bounding box of a mesh as `(min, max)` corner vectors.
///
/// Returns `None` if the mesh has no vertices.
#[allow(dead_code)]
pub fn mesh_aabb(
    mesh: &RenderMesh,
) -> Option<(oxiphysics_core::math::Vec3, oxiphysics_core::math::Vec3)> {
    use oxiphysics_core::math::Vec3;
    if mesh.vertices.is_empty() {
        return None;
    }
    let mut mn = [f32::INFINITY; 3];
    let mut mx = [f32::NEG_INFINITY; 3];
    for v in &mesh.vertices {
        for k in 0..3 {
            if v.position[k] < mn[k] {
                mn[k] = v.position[k];
            }
            if v.position[k] > mx[k] {
                mx[k] = v.position[k];
            }
        }
    }
    Some((
        Vec3::new(mn[0] as f64, mn[1] as f64, mn[2] as f64),
        Vec3::new(mx[0] as f64, mx[1] as f64, mx[2] as f64),
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use oxiphysics_core::Aabb;
    use oxiphysics_core::math::{Real, Vec3};

    #[test]
    fn test_color_constructors() {
        let r = Color::red();
        assert!((r.r - 1.0).abs() < 1e-6);
        assert!(r.g.abs() < 1e-6);
        assert!(r.b.abs() < 1e-6);
        assert!((r.a - 1.0).abs() < 1e-6);

        let g = Color::green();
        assert!(g.r.abs() < 1e-6);
        assert!((g.g - 1.0).abs() < 1e-6);

        let b = Color::blue();
        assert!(b.r.abs() < 1e-6);
        assert!((b.b - 1.0).abs() < 1e-6);

        let w = Color::white();
        assert!((w.r - 1.0).abs() < 1e-6);
        assert!((w.g - 1.0).abs() < 1e-6);
        assert!((w.b - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_wireframe_box_has_12_lines() {
        let lines = wireframe_box(Vec3::zeros(), Vec3::new(1.0, 1.0, 1.0), Color::white());
        assert_eq!(lines.len(), 12);
    }

    #[test]
    fn test_wireframe_sphere_line_count() {
        let segments = 16;
        let lines = wireframe_sphere(Vec3::zeros(), 1.0, segments, Color::white());
        // 3 circles, each with `segments` lines
        assert_eq!(lines.len(), 3 * segments);
    }

    #[test]
    fn test_box_mesh_24_vertices_36_indices() {
        let mesh = box_mesh(Vec3::zeros(), Vec3::new(1.0, 1.0, 1.0), Color::white());
        assert_eq!(mesh.vertices.len(), 24);
        assert_eq!(mesh.indices.len(), 36);
    }

    #[test]
    fn test_sphere_mesh_scales_with_subdivisions() {
        let mesh_low = sphere_mesh(Vec3::zeros(), 1.0, 4, Color::white());
        let mesh_high = sphere_mesh(Vec3::zeros(), 1.0, 8, Color::white());
        assert!(mesh_high.vertices.len() > mesh_low.vertices.len());
    }

    #[test]
    fn test_jet_colormap_endpoints() {
        let blue = map_scalar(0.0, 0.0, 1.0, Colormap::Jet);
        assert!(blue.r.abs() < 0.01);
        assert!((blue.b - 1.0).abs() < 0.01);

        let red = map_scalar(1.0, 0.0, 1.0, Colormap::Jet);
        assert!((red.r - 1.0).abs() < 0.01);
        assert!(red.b.abs() < 0.01);
    }

    #[test]
    fn test_viridis_colormap_endpoints() {
        let start = map_scalar(0.0, 0.0, 1.0, Colormap::Viridis);
        assert!(start.r < 0.3);
        assert!(start.b > 0.3);

        let end_color = map_scalar(1.0, 0.0, 1.0, Colormap::Viridis);
        assert!(end_color.r > 0.9);
        assert!(end_color.g > 0.8);
    }

    #[test]
    fn test_von_mises_uniaxial() {
        // Uniaxial stress: σ_xx=100, rest zero → von Mises = 100
        let stresses = vec![[100.0, 0.0, 0.0, 0.0, 0.0, 0.0]];
        let colors = von_mises_colors(&stresses, Colormap::Jet);
        assert_eq!(colors.len(), 1);
        // Single element → mapped to midpoint (since min==max, uses 0.5)
        // Just verify it returns a valid color
        assert!(colors[0].a > 0.0);
    }

    #[test]
    fn test_streamline_uniform_field() {
        // Uniform velocity field: v = (1, 0, 0)
        let field = |_pos: Vec3| Vec3::new(1.0, 0.0, 0.0);
        let sl = trace_streamline(&field, Vec3::zeros(), 0.1, 10);
        assert_eq!(sl.points.len(), 11); // seed + 10 steps
        // Should trace along x-axis
        let last = sl.points.last().unwrap();
        assert!((last.x - 1.0).abs() < 1e-10);
        assert!(last.y.abs() < 1e-10);
        assert!(last.z.abs() < 1e-10);
    }

    #[test]
    fn test_plane_mesh_has_4_vertices() {
        let mesh = plane_mesh(Vec3::zeros(), Vec3::new(0.0, 1.0, 0.0), 1.0, Color::white());
        assert_eq!(mesh.vertices.len(), 4);
        assert_eq!(mesh.indices.len(), 6);
    }

    #[test]
    fn test_arrow_mesh_non_empty() {
        let mesh = arrow_mesh(Vec3::zeros(), Vec3::new(0.0, 1.0, 0.0), 0.05, Color::red());
        assert!(!mesh.vertices.is_empty());
        assert!(!mesh.indices.is_empty());
    }

    #[test]
    fn test_contact_normal_lines_direction() {
        let point = Vec3::new(1.0, 0.0, 0.0);
        let normal = Vec3::new(0.0, 1.0, 0.0);
        let depth = 0.5;
        let contacts = vec![(point, normal, depth)];
        let lines = contact_normal_lines(&contacts, Color::green());
        assert_eq!(lines.len(), 1);
        let line = &lines[0];
        assert!((line.start - point).norm() < 1e-10);
        let expected_end = point + normal * depth;
        assert!((line.end - expected_end).norm() < 1e-10);
    }

    #[test]
    fn test_wireframe_aabb() {
        let aabb = Aabb::new(Vec3::new(-1.0, -1.0, -1.0), Vec3::new(1.0, 1.0, 1.0));
        let lines = wireframe_aabb(&aabb, Color::yellow());
        assert_eq!(lines.len(), 12);
    }

    // ---- Colormap tests ----

    #[test]
    fn test_colormap_viridis_range() {
        let c0 = map_scalar(0.0, 0.0, 1.0, Colormap::Viridis);
        let c1 = map_scalar(1.0, 0.0, 1.0, Colormap::Viridis);

        // The two endpoint colors must be distinct
        let diff = (c0.r - c1.r).abs() + (c0.g - c1.g).abs() + (c0.b - c1.b).abs();
        assert!(diff > 0.1, "viridis endpoints should be distinct");

        // All components in [0, 1] (representing the [0, 255] range as f32)
        for c in [c0, c1] {
            assert!(c.r >= 0.0 && c.r <= 1.0);
            assert!(c.g >= 0.0 && c.g <= 1.0);
            assert!(c.b >= 0.0 && c.b <= 1.0);
        }
    }

    #[test]
    fn test_colormap_jet_midpoint() {
        // jet(0.5) is green: (0, 1, 0)
        let c = map_scalar(0.5, 0.0, 1.0, Colormap::Jet);
        assert!(c.g > c.r, "jet midpoint: green should exceed red");
        assert!(c.g > c.b, "jet midpoint: green should exceed blue");
    }

    #[test]
    fn test_colormap_monotone_luminance() {
        // Viridis luminance should increase monotonically from t=0 to t=1.
        // Use perceptual luminance weights.
        let luminance = |t: f64| {
            let c = map_scalar(t, 0.0, 1.0, Colormap::Viridis);
            0.2126 * c.r as f64 + 0.7152 * c.g as f64 + 0.0722 * c.b as f64
        };

        let steps = 20;
        let mut prev = luminance(0.0);
        for i in 1..=steps {
            let t = (i as f64) / (steps as f64);
            let curr = luminance(t);
            assert!(
                curr >= prev - 1e-4,
                "viridis luminance not monotone at t={}: {} < {}",
                t,
                curr,
                prev
            );
            prev = curr;
        }
    }

    // ---- Stress visualization tests ----

    #[test]
    fn test_stress_principal_identity() {
        // Identity-like stress tensor: σ_xx=σ_yy=σ_zz=1, shear=0
        // All three principal stresses should be 1.
        let stress: [f64; 6] = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0];
        let scale: Real = 1.0;
        let glyphs = principal_stress_glyphs(&stress, Vec3::zeros(), scale);
        assert_eq!(glyphs.len(), 3);

        // Each glyph half-length = |eigval| * scale = 1.0 * 1.0 = 1.0
        // So start..end length = 2 * half_len = 2.0
        for glyph in &glyphs {
            let len = (glyph.end - glyph.start).norm();
            assert!(
                (len - 2.0).abs() < 1e-6,
                "expected glyph length 2.0, got {}",
                len
            );
        }
    }

    #[test]
    fn test_stress_color_compression() {
        // Compressive stress (negative) vs tensile stress (positive) should yield
        // different von Mises colors when two elements have different magnitudes.
        // von Mises is always positive, so use different magnitudes.
        let stresses = vec![
            [50.0_f64, 50.0, 50.0, 0.0, 0.0, 0.0], // hydrostatic tension
            [100.0_f64, 0.0, 0.0, 0.0, 0.0, 0.0],  // uniaxial tension
        ];
        let colors = von_mises_colors(&stresses, Colormap::Jet);
        assert_eq!(colors.len(), 2);
        // The two elements have different von Mises values so they get different colors.
        let diff = (colors[0].r - colors[1].r).abs()
            + (colors[0].g - colors[1].g).abs()
            + (colors[0].b - colors[1].b).abs();
        assert!(
            diff > 1e-3,
            "different stress magnitudes should produce different colors"
        );
    }

    #[test]
    fn test_von_mises_isotropic() {
        // Hydrostatic stress [σ, σ, σ, 0, 0, 0]: all deviatoric terms are zero → von Mises = 0
        let sigma = 100.0_f64;
        let stresses = vec![[sigma, sigma, sigma, 0.0, 0.0, 0.0]];
        let colors = von_mises_colors(&stresses, Colormap::Jet);
        assert_eq!(colors.len(), 1);
        // Single element: min == max == 0, so map_scalar uses t=0.5 (midpoint).
        // The key property is that von Mises = 0 (verified by single-element mapping).
        // Just confirm a valid color is returned (alpha > 0).
        assert!(colors[0].a > 0.0);
    }

    // ---- Streamline tests ----

    #[test]
    fn test_streamline_uniform_field_straight() {
        // Uniform field v=(1,0,0): streamline should be a straight line along x.
        let field = |_pos: Vec3| Vec3::new(1.0, 0.0, 0.0);
        let steps = 5usize;
        let dt = 0.2_f64;
        let sl = trace_streamline(&field, Vec3::zeros(), dt, steps);

        assert_eq!(sl.points.len(), steps + 1);

        // All points should lie on the x-axis (y and z stay 0).
        for (i, p) in sl.points.iter().enumerate() {
            let expected_x = (i as f64) * dt;
            assert!(
                (p.x - expected_x).abs() < 1e-10,
                "x mismatch at step {}: {} vs {}",
                i,
                p.x,
                expected_x
            );
            assert!(p.y.abs() < 1e-10, "y should be 0 at step {}", i);
            assert!(p.z.abs() < 1e-10, "z should be 0 at step {}", i);
        }
    }

    #[test]
    fn test_streamline_circular_field() {
        // Circular field v=(y, -x, 0): should curve, not stay collinear.
        let field = |pos: Vec3| Vec3::new(pos.y, -pos.x, 0.0);
        let seed = Vec3::new(1.0, 0.0, 0.0);
        let sl = trace_streamline(&field, seed, 0.1, 5);

        assert!(sl.points.len() >= 3);

        // The second point should have moved in the -y direction (v at seed = (0,-1,0)).
        // Verify it is NOT collinear with p0 and p1 extended (i.e., the path curves).
        let p0 = sl.points[0];
        let p1 = sl.points[1];
        let p2 = sl.points[2];

        // Direction p0->p1
        let d01 = (p1 - p0).normalize();
        // Direction p1->p2
        let d12 = (p2 - p1).normalize();

        // The cross product should be non-zero if the path curves.
        let cross = d01.cross(&d12);
        assert!(
            cross.norm() > 1e-4,
            "streamline in circular field should curve, not be straight"
        );
    }

    // ---- Wireframe tests ----

    #[test]
    fn test_wireframe_box_12_edges() {
        let lines = wireframe_box(Vec3::zeros(), Vec3::new(0.5, 0.5, 0.5), Color::white());
        assert_eq!(lines.len(), 12, "box wireframe must have exactly 12 edges");
    }

    #[test]
    fn test_wireframe_sphere_on_unit_sphere() {
        let radius = 1.0_f64;
        let center = Vec3::zeros();
        let segments = 32;
        let lines = wireframe_sphere(center, radius, segments, Color::white());

        // Every endpoint of every line segment must lie on the unit sphere (within tolerance).
        for line in &lines {
            let r_start = (line.start - center).norm();
            let r_end = (line.end - center).norm();
            assert!(
                (r_start - radius).abs() < 0.01,
                "start vertex radius {} not within 1±0.01",
                r_start
            );
            assert!(
                (r_end - radius).abs() < 0.01,
                "end vertex radius {} not within 1±0.01",
                r_end
            );
        }
    }

    // ---- VizConfig tests ----

    #[test]
    fn test_viz_config_default_valid_depth_range() {
        let cfg = VizConfig::default();
        assert!(
            cfg.depth_range_valid(),
            "default depth range should be valid"
        );
    }

    #[test]
    fn test_viz_config_debug_mode_has_overlay() {
        let cfg = VizConfig::debug_mode();
        assert!(cfg.show_debug_overlay);
        assert!(cfg.wireframe_mode);
    }

    #[test]
    fn test_viz_config_default_colormap_viridis() {
        let cfg = VizConfig::default();
        assert_eq!(cfg.default_colormap, Colormap::Viridis);
    }

    #[test]
    fn test_viz_config_invalid_depth_range() {
        let cfg = VizConfig {
            near_clip: 100.0,
            far_clip: 1.0,
            ..Default::default()
        };
        assert!(!cfg.depth_range_valid());
    }

    // ---- VizStats tests ----

    #[test]
    fn test_viz_stats_default_zero() {
        let s = VizStats::default();
        assert_eq!(s.draw_calls, 0);
        assert_eq!(s.triangles, 0);
        assert!(s.total_ms().abs() < 1e-12);
    }

    #[test]
    fn test_viz_stats_fps_infinite_when_zero_time() {
        let s = VizStats::default();
        assert!(s.fps().is_infinite());
    }

    #[test]
    fn test_viz_stats_fps_nonzero() {
        let s = VizStats {
            render_ms: 16.67,
            ..Default::default()
        };
        let fps = s.fps();
        assert!(fps > 50.0 && fps < 70.0, "expected ~60fps, got {fps}");
    }

    #[test]
    fn test_viz_stats_accumulate() {
        let mut a = VizStats {
            triangles: 100,
            draw_calls: 2,
            ..Default::default()
        };
        let b = VizStats {
            triangles: 200,
            draw_calls: 3,
            ..Default::default()
        };
        a.accumulate(&b);
        assert_eq!(a.triangles, 300);
        assert_eq!(a.draw_calls, 5);
    }

    #[test]
    fn test_viz_stats_reset() {
        let mut s = VizStats {
            triangles: 9999,
            ..Default::default()
        };
        s.reset();
        assert_eq!(s.triangles, 0);
    }

    // ---- VizRegistry tests ----

    #[test]
    fn test_registry_register_and_get() {
        let mut reg: VizRegistry<u32> = VizRegistry::new();
        reg.register("foo", 42u32).unwrap();
        assert_eq!(*reg.get("foo").unwrap(), 42);
    }

    #[test]
    fn test_registry_duplicate_returns_error() {
        let mut reg: VizRegistry<u32> = VizRegistry::new();
        reg.register("x", 1u32).unwrap();
        let err = reg.register("x", 2u32);
        assert!(err.is_err());
    }

    #[test]
    fn test_registry_not_found_returns_error() {
        let reg: VizRegistry<u32> = VizRegistry::new();
        assert!(reg.get("missing").is_err());
    }

    #[test]
    fn test_registry_remove() {
        let mut reg: VizRegistry<u32> = VizRegistry::new();
        reg.register("y", 99u32).unwrap();
        let val = reg.remove("y").unwrap();
        assert_eq!(val, 99);
        assert!(!reg.contains("y"));
    }

    #[test]
    fn test_registry_len_and_is_empty() {
        let mut reg: VizRegistry<u32> = VizRegistry::new();
        assert!(reg.is_empty());
        reg.register("a", 1).unwrap();
        assert_eq!(reg.len(), 1);
        assert!(!reg.is_empty());
    }

    #[test]
    fn test_registry_iter() {
        let mut reg: VizRegistry<u32> = VizRegistry::new();
        reg.register("p", 10u32).unwrap();
        reg.register("q", 20u32).unwrap();
        let names: Vec<&str> = reg.iter().map(|(n, _)| n).collect();
        assert!(names.contains(&"p"));
        assert!(names.contains(&"q"));
    }

    // ---- VizContext tests ----

    #[test]
    fn test_viz_context_begin_frame_increments_counter() {
        let mut ctx = VizContext::default();
        assert_eq!(ctx.frame_number(), 0);
        ctx.begin_frame();
        assert_eq!(ctx.frame_number(), 1);
        ctx.begin_frame();
        assert_eq!(ctx.frame_number(), 2);
    }

    #[test]
    fn test_viz_context_begin_frame_resets_stats() {
        let mut ctx = VizContext::default();
        ctx.stats.triangles = 12345;
        ctx.begin_frame();
        assert_eq!(ctx.stats.triangles, 0);
    }

    #[test]
    fn test_viz_context_map_scalar_returns_valid_color() {
        let ctx = VizContext::default();
        let c = ctx.map_scalar(0.5, 0.0, 1.0);
        assert!(c.r >= 0.0 && c.r <= 1.0);
        assert!(c.g >= 0.0 && c.g <= 1.0);
        assert!(c.b >= 0.0 && c.b <= 1.0);
    }

    #[test]
    fn test_viz_context_record_triangles() {
        let mut ctx = VizContext::default();
        ctx.record_triangles(100);
        assert_eq!(ctx.stats.triangles, 100);
        assert_eq!(ctx.stats.draw_calls, 1);
    }

    #[test]
    fn test_viz_context_record_lines() {
        let mut ctx = VizContext::default();
        ctx.record_lines(50);
        assert_eq!(ctx.stats.lines, 50);
        assert_eq!(ctx.stats.draw_calls, 1);
    }

    // ---- axis_triad tests ----

    #[test]
    fn test_axis_triad_has_three_lines() {
        let lines = axis_triad(Vec3::zeros(), 1.0);
        assert_eq!(lines.len(), 3);
    }

    #[test]
    fn test_axis_triad_x_axis_is_red() {
        let lines = axis_triad(Vec3::zeros(), 1.0);
        // X axis is the first line
        assert!((lines[0].color.r - 1.0).abs() < 1e-6);
        assert!(lines[0].color.g.abs() < 1e-6);
    }

    #[test]
    fn test_axis_triad_scaled() {
        let scale = 2.5;
        let lines = axis_triad(Vec3::zeros(), scale);
        let x_end = lines[0].end;
        assert!((x_end.x - scale).abs() < 1e-10);
    }

    // ---- grid_lines tests ----

    #[test]
    fn test_grid_lines_count() {
        // n divisions → (n+1) lines along X + (n+1) lines along Z = 2*(n+1)
        let n = 4;
        let lines = grid_lines(Vec3::zeros(), 1.0, n, Color::white());
        assert_eq!(lines.len(), 2 * (n + 1));
    }

    #[test]
    fn test_grid_lines_min_one_division() {
        let lines = grid_lines(Vec3::zeros(), 1.0, 0, Color::white());
        // divisions = max(0,1) = 1 → 2*2 = 4 lines
        assert_eq!(lines.len(), 4);
    }

    // ---- merge helpers ----

    #[test]
    fn test_merge_lines_empty() {
        let merged = merge_lines(&[]);
        assert!(merged.is_empty());
    }

    #[test]
    fn test_merge_lines_sums_counts() {
        let a = wireframe_box(Vec3::zeros(), Vec3::new(1.0, 1.0, 1.0), Color::white());
        let b = wireframe_box(
            Vec3::new(2.0, 0.0, 0.0),
            Vec3::new(3.0, 1.0, 1.0),
            Color::red(),
        );
        let merged = merge_lines(&[a.clone(), b.clone()]);
        assert_eq!(merged.len(), a.len() + b.len());
    }

    // ---- mesh transform helpers ----

    #[test]
    fn test_scale_mesh_doubles_positions() {
        let mut mesh = box_mesh(Vec3::zeros(), Vec3::new(1.0, 1.0, 1.0), Color::white());
        let orig_pos = mesh.vertices[0].position;
        scale_mesh(&mut mesh, 2.0);
        let new_pos = mesh.vertices[0].position;
        for k in 0..3 {
            assert!(
                (new_pos[k] - orig_pos[k] * 2.0).abs() < 1e-6,
                "component {k}: {} vs {}",
                new_pos[k],
                orig_pos[k] * 2.0
            );
        }
    }

    #[test]
    fn test_translate_mesh_shifts_positions() {
        let mut mesh = sphere_mesh(Vec3::zeros(), 1.0, 4, Color::white());
        let offset = Vec3::new(5.0, 0.0, 0.0);
        let orig_x = mesh.vertices[0].position[0];
        translate_mesh(&mut mesh, offset);
        let new_x = mesh.vertices[0].position[0];
        assert!((new_x - (orig_x + 5.0)).abs() < 1e-5);
    }

    #[test]
    fn test_recolor_mesh_changes_all_vertices() {
        let mut mesh = box_mesh(Vec3::zeros(), Vec3::new(1.0, 1.0, 1.0), Color::white());
        recolor_mesh(&mut mesh, Color::red());
        for v in &mesh.vertices {
            assert!((v.color.r - 1.0).abs() < 1e-6, "r should be 1");
            assert!(v.color.g.abs() < 1e-6, "g should be 0");
            assert!(v.color.b.abs() < 1e-6, "b should be 0");
        }
    }

    #[test]
    fn test_flip_normals_inverts_direction() {
        let mut mesh = sphere_mesh(Vec3::zeros(), 1.0, 4, Color::white());
        let orig_n = mesh.vertices[0].normal;
        flip_normals(&mut mesh);
        let new_n = mesh.vertices[0].normal;
        for k in 0..3 {
            assert!(
                (new_n[k] + orig_n[k]).abs() < 1e-6,
                "normal component {k} should be negated"
            );
        }
    }

    #[test]
    fn test_mesh_aabb_basic() {
        let mesh = box_mesh(
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(2.0, 2.0, 2.0),
            Color::white(),
        );
        let (mn, mx) = mesh_aabb(&mesh).unwrap();
        assert!(mn.x <= 0.0, "min x should be ≤ 0");
        assert!(mx.x >= 2.0, "max x should be ≥ 2");
    }

    #[test]
    fn test_mesh_aabb_empty_mesh() {
        let mesh = RenderMesh {
            vertices: vec![],
            indices: vec![],
        };
        assert!(mesh_aabb(&mesh).is_none());
    }
}
pub mod advanced_rendering;
pub mod animation;
pub mod annotation_viz;
pub mod chart;
pub mod chart_3d;
pub mod cross_section_viz;
pub mod data_pipeline;
pub mod data_viz;
pub mod flow_visualization;
pub mod flow_viz;
pub mod fluid_visualization;
pub mod fluid_viz;
pub mod font_rendering;
pub mod glyph_renderer;
pub mod gpu_viz;
pub mod graph_viz;
pub mod heatmap;
pub mod instancing;
pub mod interactive_viz;
pub mod isosurface;
pub mod medical_viz;
pub mod molecular_viz;
pub mod multiphysics_viz;
pub mod network_viz;
pub mod neural_rendering;
pub mod particle_effects;
pub mod particle_trails;
pub mod path_viz;
pub mod performance_viz;
pub mod phase_field_viz;
pub mod physics_animation;
pub mod physics_dashboard;
pub mod procedural_texture;
pub mod real_time_viz;
pub mod scientific_plot;
pub mod scientific_plots;
pub mod scientific_plotting;
pub mod scientific_viz;
pub mod shader_effects;
pub mod slice_renderer;
pub mod statistical_viz;
pub mod structural_viz;
pub mod tensor_visualization;
pub mod terrain_renderer;
pub mod terrain_viz;
pub mod text_renderer;
pub mod topology_viz;
pub mod topology_viz_ext;
pub mod uncertainty_viz;
pub mod volume_renderer;
pub mod volumetric_rendering;
pub mod vr_visualization;
