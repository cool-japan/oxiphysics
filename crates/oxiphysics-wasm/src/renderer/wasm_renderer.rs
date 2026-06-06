// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! `WasmRenderer` with SSAO, instanced mesh, and screen-space AABB support.

use crate::wasm_helpers::to_js_value;
use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

use super::core::{Aabb, BodyRenderHint, RendererHints};

// ---------------------------------------------------------------------------
// SsaoConfig
// ---------------------------------------------------------------------------

/// Screen Ambient Occlusion (SSAO) configuration.
///
/// All fields are primitive types, so `#[wasm_bindgen]` works without skipping.
#[wasm_bindgen]
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SsaoConfig {
    /// Whether SSAO is enabled.
    pub enabled: bool,
    /// Occlusion intensity in `[0.0, 1.0]` (0 = no occlusion, 1 = full).
    pub intensity: f64,
    /// Sampling radius in world units.
    pub radius: f64,
    /// Number of samples per pixel (higher = better quality, lower = faster).
    pub sample_count: u32,
    /// Bilateral blur kernel size (odd integer).
    pub blur_kernel_size: u32,
}

impl SsaoConfig {
    /// Default SSAO configuration (disabled, neutral settings).
    pub fn default_disabled() -> Self {
        Self {
            enabled: false,
            intensity: 0.5,
            radius: 0.5,
            sample_count: 16,
            blur_kernel_size: 3,
        }
    }

    /// Preset for high-quality SSAO.
    pub fn high_quality() -> Self {
        Self {
            enabled: true,
            intensity: 0.8,
            radius: 0.3,
            sample_count: 64,
            blur_kernel_size: 5,
        }
    }

    /// Preset for fast (low-quality) SSAO.
    pub fn fast() -> Self {
        Self {
            enabled: true,
            intensity: 0.5,
            radius: 0.5,
            sample_count: 8,
            blur_kernel_size: 3,
        }
    }

    /// Clamp intensity to `[0.0, 1.0]`.
    pub fn with_intensity(mut self, intensity: f64) -> Self {
        self.intensity = intensity.clamp(0.0, 1.0);
        self
    }
}

#[wasm_bindgen]
impl SsaoConfig {
    /// Create a default (disabled) SSAO config from JS.
    #[wasm_bindgen(js_name = "default_disabled")]
    pub fn default_disabled_js() -> SsaoConfig {
        SsaoConfig::default_disabled()
    }

    /// Create a high-quality SSAO config from JS.
    #[wasm_bindgen(js_name = "high_quality")]
    pub fn high_quality_js() -> SsaoConfig {
        SsaoConfig::high_quality()
    }

    /// Create a fast (low-quality) SSAO config from JS.
    #[wasm_bindgen(js_name = "fast")]
    pub fn fast_js() -> SsaoConfig {
        SsaoConfig::fast()
    }

    /// Return a copy of this config with the intensity clamped to `[0, 1]`.
    #[wasm_bindgen(js_name = "with_intensity")]
    pub fn with_intensity_js(self, intensity: f64) -> SsaoConfig {
        self.with_intensity(intensity)
    }
}

// ---------------------------------------------------------------------------
// InstancedMeshEntry
// ---------------------------------------------------------------------------

/// A single instanced mesh entry for GPU instanced rendering.
#[wasm_bindgen]
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct InstancedMeshEntry {
    /// Body handle this instance corresponds to.
    pub handle: u32,
    /// Position `[x, y, z]`.
    // [f64; 3] is not IntoWasmAbi for pub fields.
    #[wasm_bindgen(skip)]
    pub position: [f64; 3],
    /// Orientation quaternion `[x, y, z, w]`.
    #[wasm_bindgen(skip)]
    pub rotation: [f64; 4],
    /// Non-uniform scale `[sx, sy, sz]`.
    #[wasm_bindgen(skip)]
    pub scale: [f64; 3],
    /// RGBA color `[r, g, b, a]` in `[0.0, 1.0]`.
    #[wasm_bindgen(skip)]
    pub color: [f64; 4],
}

impl InstancedMeshEntry {
    /// Create an entry from a body render hint with default white color and unit scale.
    pub fn from_body_hint(hint: &BodyRenderHint) -> Self {
        Self {
            handle: hint.handle,
            position: hint.position,
            rotation: hint.rotation,
            scale: [1.0, 1.0, 1.0],
            color: [1.0, 1.0, 1.0, 1.0],
        }
    }

    /// Return the transform as a flat `[px, py, pz, rx, ry, rz, rw, sx, sy, sz]`.
    pub fn to_flat_transform(&self) -> Vec<f64> {
        let mut v = Vec::with_capacity(10);
        v.extend_from_slice(&self.position);
        v.extend_from_slice(&self.rotation);
        v.extend_from_slice(&self.scale);
        v
    }
}

#[wasm_bindgen]
impl InstancedMeshEntry {
    /// Return position as `[x, y, z]`.
    #[wasm_bindgen(js_name = "get_position")]
    pub fn get_position_js(&self) -> Vec<f64> {
        self.position.to_vec()
    }

    /// Return rotation quaternion as `[qx, qy, qz, qw]`.
    #[wasm_bindgen(js_name = "get_rotation")]
    pub fn get_rotation_js(&self) -> Vec<f64> {
        self.rotation.to_vec()
    }

    /// Return scale as `[sx, sy, sz]`.
    #[wasm_bindgen(js_name = "get_scale")]
    pub fn get_scale_js(&self) -> Vec<f64> {
        self.scale.to_vec()
    }

    /// Return color as `[r, g, b, a]`.
    #[wasm_bindgen(js_name = "get_color")]
    pub fn get_color_js(&self) -> Vec<f64> {
        self.color.to_vec()
    }

    /// Return the flat 10-float transform buffer.
    #[wasm_bindgen(js_name = "to_flat_transform")]
    pub fn to_flat_transform_js(&self) -> Vec<f64> {
        self.to_flat_transform()
    }
}

// ---------------------------------------------------------------------------
// ScreenSpaceAabb
// ---------------------------------------------------------------------------

/// Screen-space (2D) bounding box of a projected 3D AABB.
#[wasm_bindgen]
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ScreenSpaceAabb {
    /// Body handle.
    pub handle: u32,
    /// Minimum screen-space pixel `[x, y]`.
    // [f64; 2] is not IntoWasmAbi for pub fields.
    #[wasm_bindgen(skip)]
    pub min_px: [f64; 2],
    /// Maximum screen-space pixel `[x, y]`.
    #[wasm_bindgen(skip)]
    pub max_px: [f64; 2],
    /// Screen-space width in pixels.
    pub width_px: f64,
    /// Screen-space height in pixels.
    pub height_px: f64,
    /// Whether the body is partially or fully behind the camera.
    pub clipped: bool,
}

impl ScreenSpaceAabb {
    /// Compute a screen-space AABB for a body given its world-space AABB
    /// and a simple orthographic camera projection.
    ///
    /// `viewport_size` is `[width_px, height_px]`.
    /// `camera_scale` is pixels per world unit.
    /// `camera_offset` is the world-space origin mapped to the screen center.
    pub fn from_aabb_orthographic(
        aabb: &Aabb,
        viewport_size: [f64; 2],
        camera_scale: f64,
        camera_offset: [f64; 2],
    ) -> Self {
        let cx = viewport_size[0] * 0.5;
        let cy = viewport_size[1] * 0.5;

        let project = |wx: f64, wy: f64| -> [f64; 2] {
            [
                cx + (wx - camera_offset[0]) * camera_scale,
                cy - (wy - camera_offset[1]) * camera_scale, // Y flipped
            ]
        };

        let corners = [
            project(aabb.min[0], aabb.min[1]),
            project(aabb.min[0], aabb.max[1]),
            project(aabb.max[0], aabb.min[1]),
            project(aabb.max[0], aabb.max[1]),
        ];

        let min_px_x = corners.iter().map(|c| c[0]).fold(f64::INFINITY, f64::min);
        let min_px_y = corners.iter().map(|c| c[1]).fold(f64::INFINITY, f64::min);
        let max_px_x = corners
            .iter()
            .map(|c| c[0])
            .fold(f64::NEG_INFINITY, f64::max);
        let max_px_y = corners
            .iter()
            .map(|c| c[1])
            .fold(f64::NEG_INFINITY, f64::max);

        let clipped = max_px_x < 0.0
            || min_px_x > viewport_size[0]
            || max_px_y < 0.0
            || min_px_y > viewport_size[1];

        Self {
            handle: aabb.handle,
            min_px: [min_px_x, min_px_y],
            max_px: [max_px_x, max_px_y],
            width_px: (max_px_x - min_px_x).max(0.0),
            height_px: (max_px_y - min_px_y).max(0.0),
            clipped,
        }
    }
}

#[wasm_bindgen]
impl ScreenSpaceAabb {
    /// Return the minimum screen-space pixel `[x, y]`.
    #[wasm_bindgen(js_name = "get_min_px")]
    pub fn get_min_px_js(&self) -> Vec<f64> {
        self.min_px.to_vec()
    }

    /// Return the maximum screen-space pixel `[x, y]`.
    #[wasm_bindgen(js_name = "get_max_px")]
    pub fn get_max_px_js(&self) -> Vec<f64> {
        self.max_px.to_vec()
    }

    /// Serialize to a `JsValue` object.
    #[wasm_bindgen(js_name = "to_js")]
    pub fn to_js_js(&self) -> Result<JsValue, JsValue> {
        to_js_value(self)
    }
}

// ---------------------------------------------------------------------------
// WasmRenderer
// ---------------------------------------------------------------------------

/// WASM renderer with SSAO, instanced mesh, and screen-space AABB support.
#[wasm_bindgen]
#[derive(Debug, Clone)]
pub struct WasmRenderer {
    /// Current SSAO configuration.
    ssao: SsaoConfig,
    /// Current instanced mesh entries (one per active body).
    instances: Vec<InstancedMeshEntry>,
    /// Viewport width in pixels.
    viewport_width: f64,
    /// Viewport height in pixels.
    viewport_height: f64,
}

impl WasmRenderer {
    /// Create a new renderer with the given viewport dimensions.
    pub fn new(viewport_width: f64, viewport_height: f64) -> Self {
        Self {
            ssao: SsaoConfig::default_disabled(),
            instances: Vec::new(),
            viewport_width,
            viewport_height,
        }
    }

    // -----------------------------------------------------------------------
    // Ambient occlusion
    // -----------------------------------------------------------------------

    /// Set the SSAO configuration.
    pub fn set_ambient_occlusion(&mut self, config: SsaoConfig) {
        self.ssao = config;
    }

    /// Return the current SSAO configuration.
    pub fn ssao_config(&self) -> &SsaoConfig {
        &self.ssao
    }

    /// Enable or disable SSAO without changing other settings.
    pub fn toggle_ssao(&mut self, enabled: bool) {
        self.ssao.enabled = enabled;
    }

    // -----------------------------------------------------------------------
    // Instanced rendering
    // -----------------------------------------------------------------------

    /// Update the instanced mesh buffer from a set of `RendererHints`.
    ///
    /// One instance entry is created per body hint.
    pub fn update_instanced_mesh(&mut self, hints: &RendererHints) {
        self.instances = hints
            .body_hints
            .iter()
            .map(InstancedMeshEntry::from_body_hint)
            .collect();
    }

    /// Return the instanced mesh entries.
    pub fn get_instances(&self) -> &[InstancedMeshEntry] {
        &self.instances
    }

    /// Return the instance transforms as a flat buffer for GPU upload.
    ///
    /// Format per instance: `[px, py, pz, rx, ry, rz, rw, sx, sy, sz]` (10 floats).
    pub fn instance_buffer_flat(&self) -> Vec<f64> {
        self.instances
            .iter()
            .flat_map(|inst| inst.to_flat_transform())
            .collect()
    }

    // -----------------------------------------------------------------------
    // Screen-space bounding box
    // -----------------------------------------------------------------------

    /// Compute screen-space bounding boxes for all body AABBs in the hints,
    /// using an orthographic projection.
    ///
    /// `camera_scale` is pixels per world unit (zoom level).
    /// `camera_offset` is the world position mapped to the screen center.
    pub fn compute_screen_space_bounding_box(
        &self,
        hints: &RendererHints,
        camera_scale: f64,
        camera_offset: [f64; 2],
    ) -> Vec<ScreenSpaceAabb> {
        hints
            .body_hints
            .iter()
            .filter_map(|hint| hint.aabb.as_ref())
            .map(|aabb| {
                ScreenSpaceAabb::from_aabb_orthographic(
                    aabb,
                    [self.viewport_width, self.viewport_height],
                    camera_scale,
                    camera_offset,
                )
            })
            .collect()
    }
}

// ---------------------------------------------------------------------------
// WasmRenderer — wasm_bindgen impl
// ---------------------------------------------------------------------------

#[wasm_bindgen]
impl WasmRenderer {
    /// Create a new `WasmRenderer` with given viewport dimensions.
    #[wasm_bindgen(constructor)]
    pub fn wasm_new(viewport_width: f64, viewport_height: f64) -> WasmRenderer {
        WasmRenderer::new(viewport_width, viewport_height)
    }

    /// Set the SSAO configuration.
    #[wasm_bindgen(js_name = "set_ambient_occlusion")]
    pub fn set_ambient_occlusion_js(&mut self, config: SsaoConfig) {
        self.set_ambient_occlusion(config);
    }

    /// Return a copy of the current SSAO configuration.
    #[wasm_bindgen(js_name = "ssao_config")]
    pub fn ssao_config_js(&self) -> SsaoConfig {
        self.ssao.clone()
    }

    /// Enable or disable SSAO.
    #[wasm_bindgen(js_name = "toggle_ssao")]
    pub fn toggle_ssao_js(&mut self, enabled: bool) {
        self.toggle_ssao(enabled);
    }

    /// Update the instanced mesh buffer from `RendererHints`.
    #[wasm_bindgen(js_name = "update_instanced_mesh")]
    pub fn update_instanced_mesh_js(&mut self, hints: &RendererHints) {
        self.update_instanced_mesh(hints);
    }

    /// Return instance count.
    #[wasm_bindgen(js_name = "instance_count")]
    pub fn instance_count_js(&self) -> u32 {
        self.instances.len() as u32
    }

    /// Return the flat GPU instance buffer (10 floats per instance).
    #[wasm_bindgen(js_name = "instance_buffer_flat")]
    pub fn instance_buffer_flat_js(&self) -> Vec<f64> {
        self.instance_buffer_flat()
    }

    /// Compute screen-space AABBs from hints, camera scale, and camera offset (flat `[ox, oy]`).
    ///
    /// Returns a `JsValue` array of `ScreenSpaceAabb`-like objects.
    #[wasm_bindgen(js_name = "compute_screen_space_aabbs")]
    pub fn compute_screen_space_aabbs_js(
        &self,
        hints: &RendererHints,
        camera_scale: f64,
        camera_offset: Vec<f64>,
    ) -> Result<JsValue, JsValue> {
        if camera_offset.len() != 2 {
            return Err(JsValue::from_str("camera_offset must have 2 elements"));
        }
        let offset = [camera_offset[0], camera_offset[1]];
        let boxes = self.compute_screen_space_bounding_box(hints, camera_scale, offset);
        to_js_value(&boxes)
    }
}

// ===========================================================================
// Tests
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::WasmPhysicsEngine;
    use crate::renderer::core::RendererHints;

    // --- SsaoConfig ---

    #[test]
    fn test_ssao_default_disabled() {
        let cfg = SsaoConfig::default_disabled();
        assert!(!cfg.enabled, "default SSAO should be disabled");
    }

    #[test]
    fn test_ssao_high_quality_preset() {
        let cfg = SsaoConfig::high_quality();
        assert!(cfg.enabled);
        assert_eq!(cfg.sample_count, 64);
        assert!(cfg.intensity > 0.5);
    }

    #[test]
    fn test_ssao_intensity_clamped() {
        let cfg = SsaoConfig::default_disabled().with_intensity(2.0);
        assert!(
            (cfg.intensity - 1.0).abs() < 1e-10,
            "intensity clamped to 1.0"
        );
        let cfg2 = SsaoConfig::default_disabled().with_intensity(-1.0);
        assert!(
            (cfg2.intensity - 0.0).abs() < 1e-10,
            "intensity clamped to 0.0"
        );
    }

    #[test]
    fn test_set_ambient_occlusion() {
        let mut renderer = WasmRenderer::new(800.0, 600.0);
        renderer.set_ambient_occlusion(SsaoConfig::high_quality());
        assert!(renderer.ssao_config().enabled);
        assert_eq!(renderer.ssao_config().sample_count, 64);
    }

    #[test]
    fn test_toggle_ssao() {
        let mut renderer = WasmRenderer::new(800.0, 600.0);
        renderer.toggle_ssao(true);
        assert!(renderer.ssao_config().enabled);
        renderer.toggle_ssao(false);
        assert!(!renderer.ssao_config().enabled);
    }

    // --- Instanced mesh ---

    #[test]
    fn test_update_instanced_mesh_empty() {
        let engine = WasmPhysicsEngine::new(0.0, -9.81, 0.0);
        let hints = RendererHints::from_engine(&engine, true, false, false);
        let mut renderer = WasmRenderer::new(800.0, 600.0);
        renderer.update_instanced_mesh(&hints);
        assert!(renderer.get_instances().is_empty());
    }

    #[test]
    fn test_update_instanced_mesh_one_body() {
        let mut engine = WasmPhysicsEngine::new(0.0, 0.0, 0.0);
        engine.add_dynamic_body(1.0, 3.0, 5.0, 0.0);
        let hints = RendererHints::from_engine(&engine, true, false, false);
        let mut renderer = WasmRenderer::new(800.0, 600.0);
        renderer.update_instanced_mesh(&hints);
        assert_eq!(renderer.get_instances().len(), 1);
        let inst = &renderer.get_instances()[0];
        assert!((inst.position[0] - 3.0).abs() < 1e-12);
        assert!((inst.position[1] - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_instance_buffer_flat_length() {
        let mut engine = WasmPhysicsEngine::new(0.0, 0.0, 0.0);
        engine.add_dynamic_body(1.0, 0.0, 0.0, 0.0);
        engine.add_dynamic_body(2.0, 1.0, 0.0, 0.0);
        let hints = RendererHints::from_engine(&engine, true, false, false);
        let mut renderer = WasmRenderer::new(800.0, 600.0);
        renderer.update_instanced_mesh(&hints);
        let flat = renderer.instance_buffer_flat();
        assert_eq!(flat.len(), 20, "flat buffer should have 20 elements");
    }

    // --- Screen-space AABB ---

    #[test]
    fn test_screen_space_aabb_center_body() {
        let aabb = Aabb::from_sphere(0, [0.0, 0.0, 0.0], 0.5);
        let ss = ScreenSpaceAabb::from_aabb_orthographic(&aabb, [800.0, 600.0], 100.0, [0.0, 0.0]);
        let cx = (ss.min_px[0] + ss.max_px[0]) * 0.5;
        let cy = (ss.min_px[1] + ss.max_px[1]) * 0.5;
        assert!((cx - 400.0).abs() < 1.0, "screen cx = {}", cx);
        assert!((cy - 300.0).abs() < 1.0, "screen cy = {}", cy);
        assert!(!ss.clipped);
    }

    #[test]
    fn test_screen_space_aabb_size() {
        let aabb = Aabb::from_sphere(1, [0.0, 0.0, 0.0], 1.0);
        let ss = ScreenSpaceAabb::from_aabb_orthographic(&aabb, [800.0, 600.0], 50.0, [0.0, 0.0]);
        assert!(
            (ss.width_px - 100.0).abs() < 1.0,
            "width_px = {}",
            ss.width_px
        );
        assert!(
            (ss.height_px - 100.0).abs() < 1.0,
            "height_px = {}",
            ss.height_px
        );
    }

    #[test]
    fn test_compute_screen_space_bounding_box_empty() {
        let engine = WasmPhysicsEngine::new(0.0, 0.0, 0.0);
        let hints = RendererHints::from_engine(&engine, false, false, false);
        let renderer = WasmRenderer::new(800.0, 600.0);
        let boxes = renderer.compute_screen_space_bounding_box(&hints, 100.0, [0.0, 0.0]);
        assert!(boxes.is_empty());
    }

    #[test]
    fn test_compute_screen_space_bounding_box_one_body() {
        let mut engine = WasmPhysicsEngine::new(0.0, 0.0, 0.0);
        engine.add_dynamic_body(1.0, 0.0, 0.0, 0.0);
        let hints = RendererHints::from_engine(&engine, true, false, false);
        let renderer = WasmRenderer::new(800.0, 600.0);
        let boxes = renderer.compute_screen_space_bounding_box(&hints, 100.0, [0.0, 0.0]);
        assert_eq!(boxes.len(), 1);
        assert!(boxes[0].width_px > 0.0);
        assert!(boxes[0].height_px > 0.0);
    }

    // ---------------------------------------------------------------------------
    // Integration tests (Slice W6)
    // ---------------------------------------------------------------------------

    /// Verify SSAO config serializes / deserializes faithfully.
    #[test]
    fn test_ssao_config_json_roundtrip() {
        let cfg = SsaoConfig::high_quality();
        let json = serde_json::to_string(&cfg).expect("serialize");
        let back: SsaoConfig = serde_json::from_str(&json).expect("deserialize");
        assert_eq!(back.sample_count, cfg.sample_count);
        assert!((back.intensity - cfg.intensity).abs() < 1e-10);
    }

    /// Verify that instanced mesh updates are idempotent (calling twice gives same result).
    #[test]
    fn test_update_instanced_mesh_idempotent() {
        let mut engine = WasmPhysicsEngine::new(0.0, 0.0, 0.0);
        engine.add_dynamic_body(1.0, 1.0, 2.0, 3.0);
        let hints = RendererHints::from_engine(&engine, true, false, false);
        let mut renderer = WasmRenderer::new(800.0, 600.0);
        renderer.update_instanced_mesh(&hints);
        let count_first = renderer.get_instances().len();
        renderer.update_instanced_mesh(&hints);
        let count_second = renderer.get_instances().len();
        assert_eq!(count_first, count_second);
    }
}
