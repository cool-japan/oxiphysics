// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! WebGPU render pipeline descriptor helpers, instanced rendering data,
//! shadow map pass configuration, deferred shading G-buffer pass,
//! particle rendering pass, post-processing passes, and the render frame graph.

use crate::wasm_helpers::to_js_value;
use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

// ===========================================================================
// WebGPU render pipeline descriptor helpers
// ===========================================================================

/// Blend mode for a render pipeline.
///
/// Unit-only enum — can be annotated with `#[wasm_bindgen]` directly.
#[wasm_bindgen(js_name = "RendererBlendMode")]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum BlendMode {
    /// No blending — overwrite destination.
    Opaque,
    /// Standard alpha blending.
    AlphaBlend,
    /// Additive blending.
    Additive,
    /// Pre-multiplied alpha blending.
    PremultipliedAlpha,
}

/// Depth test comparison function.
///
/// Unit-only enum — can be annotated with `#[wasm_bindgen]` directly.
#[wasm_bindgen(js_name = "RendererDepthCompare")]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum DepthCompare {
    /// Always pass.
    Always,
    /// Pass if new depth < stored depth.
    Less,
    /// Pass if new depth <= stored depth.
    LessEqual,
    /// Pass if new depth > stored depth.
    Greater,
    /// Pass if new depth == stored depth.
    Equal,
    /// Never pass.
    Never,
}

/// A descriptor for a render pipeline (GPU render pass configuration).
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RenderPipelineDesc {
    /// Human-readable label.
    // String is not IntoWasmAbi for pub struct fields — expose via getter.
    #[wasm_bindgen(skip)]
    pub label: String,
    /// Blend mode.
    pub blend_mode: BlendMode,
    /// Whether depth testing is enabled.
    pub depth_test: bool,
    /// Whether depth writing is enabled.
    pub depth_write: bool,
    /// Depth comparison function.
    pub depth_compare: DepthCompare,
    /// Whether back-face culling is enabled.
    pub backface_culling: bool,
    /// Number of color render targets.
    pub color_target_count: u32,
    /// MSAA sample count (1, 2, 4, or 8).
    pub sample_count: u32,
}

impl RenderPipelineDesc {
    /// Default opaque pipeline (depth test on, writes on, back-face culling on).
    pub fn opaque(label: impl Into<String>) -> Self {
        Self {
            label: label.into(),
            blend_mode: BlendMode::Opaque,
            depth_test: true,
            depth_write: true,
            depth_compare: DepthCompare::Less,
            backface_culling: true,
            color_target_count: 1,
            sample_count: 1,
        }
    }

    /// Transparent pipeline (alpha blend, depth test on, depth write off).
    pub fn transparent(label: impl Into<String>) -> Self {
        Self {
            label: label.into(),
            blend_mode: BlendMode::AlphaBlend,
            depth_test: true,
            depth_write: false,
            depth_compare: DepthCompare::Less,
            backface_culling: false,
            color_target_count: 1,
            sample_count: 1,
        }
    }

    /// Shadow-map pipeline (depth-only, no color targets).
    pub fn shadow_map(label: impl Into<String>) -> Self {
        Self {
            label: label.into(),
            blend_mode: BlendMode::Opaque,
            depth_test: true,
            depth_write: true,
            depth_compare: DepthCompare::Less,
            backface_culling: true,
            color_target_count: 0,
            sample_count: 1,
        }
    }

    /// Post-processing pipeline (no depth test, no depth write, full-screen quad).
    pub fn post_process(label: impl Into<String>) -> Self {
        Self {
            label: label.into(),
            blend_mode: BlendMode::Opaque,
            depth_test: false,
            depth_write: false,
            depth_compare: DepthCompare::Always,
            backface_culling: false,
            color_target_count: 1,
            sample_count: 1,
        }
    }

    /// Enable MSAA with the given sample count.
    pub fn with_msaa(mut self, samples: u32) -> Self {
        self.sample_count = samples.clamp(1, 8);
        self
    }

    /// Set the number of color targets (for MRT / deferred shading).
    pub fn with_color_targets(mut self, count: u32) -> Self {
        self.color_target_count = count;
        self
    }
}

impl Default for RenderPipelineDesc {
    fn default() -> Self {
        Self::opaque("default")
    }
}

#[wasm_bindgen]
impl RenderPipelineDesc {
    /// Create an opaque pipeline with the given label (JS constructor helper).
    #[wasm_bindgen(js_name = "opaque")]
    pub fn opaque_js(label: String) -> RenderPipelineDesc {
        RenderPipelineDesc::opaque(label)
    }

    /// Create a transparent pipeline with the given label.
    #[wasm_bindgen(js_name = "transparent")]
    pub fn transparent_js(label: String) -> RenderPipelineDesc {
        RenderPipelineDesc::transparent(label)
    }

    /// Create a shadow-map pipeline with the given label.
    #[wasm_bindgen(js_name = "shadow_map")]
    pub fn shadow_map_js(label: String) -> RenderPipelineDesc {
        RenderPipelineDesc::shadow_map(label)
    }

    /// Create a post-process pipeline with the given label.
    #[wasm_bindgen(js_name = "post_process")]
    pub fn post_process_js(label: String) -> RenderPipelineDesc {
        RenderPipelineDesc::post_process(label)
    }

    /// Return the pipeline label as a `String`.
    #[wasm_bindgen(js_name = "get_label")]
    pub fn get_label_js(&self) -> String {
        self.label.clone()
    }

    /// Set the MSAA sample count (clamped 1..=8).
    #[wasm_bindgen(js_name = "with_msaa")]
    pub fn with_msaa_js(self, samples: u32) -> RenderPipelineDesc {
        self.with_msaa(samples)
    }

    /// Set the number of color render targets.
    #[wasm_bindgen(js_name = "with_color_targets")]
    pub fn with_color_targets_js(self, count: u32) -> RenderPipelineDesc {
        self.with_color_targets(count)
    }

    /// Serialize to a `JsValue` object.
    #[wasm_bindgen(js_name = "to_js")]
    pub fn to_js_js(&self) -> Result<JsValue, JsValue> {
        to_js_value(self)
    }
}

// ===========================================================================
// Instanced rendering data
// ===========================================================================

/// Per-instance data for GPU instanced rendering (transform + color).
///
/// Layout: `[px, py, pz, qx, qy, qz, qw, sx, sy, sz, r, g, b, a]` (14 floats).
#[wasm_bindgen]
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct InstanceData {
    /// World position `[x, y, z]`.
    // [f64; 3] is not IntoWasmAbi — expose via JS getter.
    #[wasm_bindgen(skip)]
    pub position: [f64; 3],
    /// Orientation quaternion `[x, y, z, w]`.
    #[wasm_bindgen(skip)]
    pub rotation: [f64; 4],
    /// Scale `[sx, sy, sz]`.
    #[wasm_bindgen(skip)]
    pub scale: [f64; 3],
    /// RGBA color `[r, g, b, a]` in `[0, 1]`.
    #[wasm_bindgen(skip)]
    pub color: [f64; 4],
}

impl InstanceData {
    /// Create with default unit scale and white color.
    pub fn new(position: [f64; 3], rotation: [f64; 4]) -> Self {
        Self {
            position,
            rotation,
            scale: [1.0; 3],
            color: [1.0; 4],
        }
    }

    /// Pack to flat 14-float buffer.
    pub fn to_flat14(&self) -> [f64; 14] {
        [
            self.position[0],
            self.position[1],
            self.position[2],
            self.rotation[0],
            self.rotation[1],
            self.rotation[2],
            self.rotation[3],
            self.scale[0],
            self.scale[1],
            self.scale[2],
            self.color[0],
            self.color[1],
            self.color[2],
            self.color[3],
        ]
    }

    /// Pack many instances into a flat `Vec<f64>` (14 floats per instance).
    pub fn pack_flat(instances: &[Self]) -> Vec<f64> {
        instances.iter().flat_map(|i| i.to_flat14()).collect()
    }
}

#[wasm_bindgen]
impl InstanceData {
    /// Create an `InstanceData` with the given position (flat `[x,y,z]`) and
    /// orientation quaternion (flat `[qx,qy,qz,qw]`), with unit scale and white.
    ///
    /// Returns an error if the input slices don't have the correct length.
    #[wasm_bindgen(js_name = "new_from_vecs")]
    pub fn new_from_vecs_js(
        position: Vec<f64>,
        rotation: Vec<f64>,
    ) -> Result<InstanceData, JsValue> {
        if position.len() != 3 {
            return Err(JsValue::from_str("position must have 3 elements"));
        }
        if rotation.len() != 4 {
            return Err(JsValue::from_str("rotation must have 4 elements"));
        }
        Ok(InstanceData::new(
            [position[0], position[1], position[2]],
            [rotation[0], rotation[1], rotation[2], rotation[3]],
        ))
    }

    /// Return position as `[x, y, z]` `Vec<f64>`.
    #[wasm_bindgen(js_name = "get_position")]
    pub fn get_position_js(&self) -> Vec<f64> {
        self.position.to_vec()
    }

    /// Return rotation quaternion as `[qx, qy, qz, qw]` `Vec<f64>`.
    #[wasm_bindgen(js_name = "get_rotation")]
    pub fn get_rotation_js(&self) -> Vec<f64> {
        self.rotation.to_vec()
    }

    /// Return scale as `[sx, sy, sz]` `Vec<f64>`.
    #[wasm_bindgen(js_name = "get_scale")]
    pub fn get_scale_js(&self) -> Vec<f64> {
        self.scale.to_vec()
    }

    /// Return color as `[r, g, b, a]` `Vec<f64>`.
    #[wasm_bindgen(js_name = "get_color")]
    pub fn get_color_js(&self) -> Vec<f64> {
        self.color.to_vec()
    }

    /// Pack to a flat 14-float `Vec<f64>`.
    #[wasm_bindgen(js_name = "to_flat14")]
    pub fn to_flat14_js(&self) -> Vec<f64> {
        self.to_flat14().to_vec()
    }
}

// ===========================================================================
// Shadow map pass configuration
// ===========================================================================

/// Configuration for a single shadow map render pass.
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ShadowMapPass {
    /// Unique pass identifier.
    pub id: u32,
    /// Resolution of the shadow map texture (width = height = resolution).
    pub resolution: u32,
    /// Light direction (normalized `[x, y, z]`).
    // [f64; 3] is not IntoWasmAbi — expose via JS getter.
    #[wasm_bindgen(skip)]
    pub light_dir: [f64; 3],
    /// Orthographic projection half-size (world units from center).
    pub ortho_size: f64,
    /// Near clip distance for the shadow camera.
    pub near: f64,
    /// Far clip distance for the shadow camera.
    pub far: f64,
    /// Depth bias to avoid self-shadowing artefacts.
    pub depth_bias: f64,
}

impl ShadowMapPass {
    /// Create a shadow map pass for a directional light.
    pub fn directional(id: u32, light_dir: [f64; 3], ortho_size: f64) -> Self {
        let len = (light_dir[0] * light_dir[0]
            + light_dir[1] * light_dir[1]
            + light_dir[2] * light_dir[2])
            .sqrt();
        let ld = if len > 1e-15 {
            [light_dir[0] / len, light_dir[1] / len, light_dir[2] / len]
        } else {
            [0.0, -1.0, 0.0]
        };
        Self {
            id,
            resolution: 1024,
            light_dir: ld,
            ortho_size,
            near: 0.1,
            far: 200.0,
            depth_bias: 0.005,
        }
    }

    /// Return the pipeline descriptor for this shadow pass.
    pub fn pipeline(&self) -> RenderPipelineDesc {
        RenderPipelineDesc::shadow_map(format!("shadow_map_{}", self.id))
    }
}

#[wasm_bindgen]
impl ShadowMapPass {
    /// Create a directional shadow map pass from JS (light direction as flat `[x,y,z]`).
    #[wasm_bindgen(js_name = "directional")]
    pub fn directional_js(
        id: u32,
        light_dir: Vec<f64>,
        ortho_size: f64,
    ) -> Result<ShadowMapPass, JsValue> {
        if light_dir.len() != 3 {
            return Err(JsValue::from_str("light_dir must have 3 elements"));
        }
        Ok(ShadowMapPass::directional(
            id,
            [light_dir[0], light_dir[1], light_dir[2]],
            ortho_size,
        ))
    }

    /// Return light direction as `Vec<f64>`.
    #[wasm_bindgen(js_name = "get_light_dir")]
    pub fn get_light_dir_js(&self) -> Vec<f64> {
        self.light_dir.to_vec()
    }

    /// Return the pipeline descriptor as a `JsValue`.
    #[wasm_bindgen(js_name = "pipeline")]
    pub fn pipeline_js(&self) -> Result<JsValue, JsValue> {
        to_js_value(&self.pipeline())
    }
}

// ===========================================================================
// Deferred shading G-buffer pass
// ===========================================================================

/// G-buffer layout for deferred shading.
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GBufferLayout {
    /// Width in pixels.
    pub width: u32,
    /// Height in pixels.
    pub height: u32,
    /// Whether to include a velocity/motion-vector render target.
    pub include_velocity: bool,
    /// Whether to include an emissive render target.
    pub include_emissive: bool,
}

impl GBufferLayout {
    /// Create a standard G-buffer layout.
    pub fn new(width: u32, height: u32) -> Self {
        Self {
            width,
            height,
            include_velocity: false,
            include_emissive: false,
        }
    }

    /// Number of render targets in this G-buffer.
    pub fn render_target_count(&self) -> u32 {
        let mut n = 3u32; // albedo, normals, depth
        if self.include_velocity {
            n += 1;
        }
        if self.include_emissive {
            n += 1;
        }
        n
    }

    /// Return the geometry pass pipeline for filling the G-buffer.
    pub fn geometry_pass_pipeline(&self) -> RenderPipelineDesc {
        RenderPipelineDesc::opaque("deferred_geometry")
            .with_color_targets(self.render_target_count())
    }

    /// Return the lighting pass pipeline (full-screen quad reading G-buffer).
    pub fn lighting_pass_pipeline(&self) -> RenderPipelineDesc {
        RenderPipelineDesc::post_process("deferred_lighting")
    }
}

#[wasm_bindgen]
impl GBufferLayout {
    /// Create a standard G-buffer layout from JS.
    #[wasm_bindgen(constructor)]
    pub fn wasm_new(width: u32, height: u32) -> GBufferLayout {
        GBufferLayout::new(width, height)
    }

    /// Return the number of render targets.
    #[wasm_bindgen(js_name = "render_target_count")]
    pub fn render_target_count_js(&self) -> u32 {
        self.render_target_count()
    }

    /// Return the geometry pass pipeline as `JsValue`.
    #[wasm_bindgen(js_name = "geometry_pass_pipeline")]
    pub fn geometry_pass_pipeline_js(&self) -> Result<JsValue, JsValue> {
        to_js_value(&self.geometry_pass_pipeline())
    }

    /// Return the lighting pass pipeline as `JsValue`.
    #[wasm_bindgen(js_name = "lighting_pass_pipeline")]
    pub fn lighting_pass_pipeline_js(&self) -> Result<JsValue, JsValue> {
        to_js_value(&self.lighting_pass_pipeline())
    }
}

// ===========================================================================
// Particle rendering pass
// ===========================================================================

/// Configuration for a point-sprite or mesh particle rendering pass.
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ParticleRenderPass {
    /// Maximum number of particles supported.
    pub max_particles: u32,
    /// Particle size in world units.
    pub particle_size: f64,
    /// Blend mode (usually `Additive` or `AlphaBlend`).
    pub blend_mode: BlendMode,
    /// Whether to sort particles back-to-front (for correct alpha blending).
    pub sort_back_to_front: bool,
}

impl ParticleRenderPass {
    /// Create a default additive particle pass.
    pub fn additive(max_particles: u32) -> Self {
        Self {
            max_particles,
            particle_size: 0.05,
            blend_mode: BlendMode::Additive,
            sort_back_to_front: false,
        }
    }

    /// Create an alpha-blended particle pass (with sorting).
    pub fn alpha_blended(max_particles: u32) -> Self {
        Self {
            max_particles,
            particle_size: 0.05,
            blend_mode: BlendMode::AlphaBlend,
            sort_back_to_front: true,
        }
    }

    /// Pipeline descriptor for this particle pass.
    pub fn pipeline(&self) -> RenderPipelineDesc {
        let mut p = RenderPipelineDesc::transparent("particles");
        p.blend_mode = self.blend_mode;
        p
    }

    /// Sort particle positions by distance to camera (farthest first).
    pub fn sort_particles_by_depth(positions: &[[f64; 3]], camera_pos: [f64; 3]) -> Vec<usize> {
        let mut indices: Vec<usize> = (0..positions.len()).collect();
        indices.sort_by(|&a, &b| {
            let da = dist_sq(positions[a], camera_pos);
            let db = dist_sq(positions[b], camera_pos);
            db.partial_cmp(&da).unwrap_or(std::cmp::Ordering::Equal)
        });
        indices
    }
}

#[wasm_bindgen]
impl ParticleRenderPass {
    /// Create an additive particle pass from JS.
    #[wasm_bindgen(js_name = "additive")]
    pub fn additive_js(max_particles: u32) -> ParticleRenderPass {
        ParticleRenderPass::additive(max_particles)
    }

    /// Create an alpha-blended particle pass from JS.
    #[wasm_bindgen(js_name = "alpha_blended")]
    pub fn alpha_blended_js(max_particles: u32) -> ParticleRenderPass {
        ParticleRenderPass::alpha_blended(max_particles)
    }

    /// Sort particle positions (flat `Vec<f64>` triplets) by depth, returning flat index order.
    ///
    /// `positions_flat` must have length divisible by 3.
    /// `camera_pos` must have length 3.
    #[wasm_bindgen(js_name = "sort_by_depth")]
    pub fn sort_by_depth_js(
        positions_flat: Vec<f64>,
        camera_pos: Vec<f64>,
    ) -> Result<Vec<u32>, JsValue> {
        if !positions_flat.len().is_multiple_of(3) {
            return Err(JsValue::from_str(
                "positions_flat length must be divisible by 3",
            ));
        }
        if camera_pos.len() != 3 {
            return Err(JsValue::from_str("camera_pos must have 3 elements"));
        }
        let positions: Vec<[f64; 3]> = positions_flat
            .chunks_exact(3)
            .map(|c| [c[0], c[1], c[2]])
            .collect();
        let cam = [camera_pos[0], camera_pos[1], camera_pos[2]];
        let indices = ParticleRenderPass::sort_particles_by_depth(&positions, cam);
        Ok(indices.into_iter().map(|i| i as u32).collect())
    }
}

#[inline]
fn dist_sq(a: [f64; 3], b: [f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    dx * dx + dy * dy + dz * dz
}

// ===========================================================================
// Post-processing pass
// ===========================================================================

/// A post-processing effect pass (tone-mapping, bloom, FXAA, etc.).
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PostProcessPass {
    /// Effect name / shader identifier.
    // String is not IntoWasmAbi for pub fields — expose via getter.
    #[wasm_bindgen(skip)]
    pub effect: String,
    /// Whether this pass is enabled.
    pub enabled: bool,
    /// Arbitrary float parameters (effect-specific).
    // Vec<f64> is not IntoWasmAbi for pub fields — expose via getter.
    #[wasm_bindgen(skip)]
    pub params: Vec<f64>,
}

impl PostProcessPass {
    /// Create an FXAA anti-aliasing pass.
    pub fn fxaa() -> Self {
        Self {
            effect: "fxaa".to_string(),
            enabled: true,
            params: vec![0.0833, 0.166, 0.75],
        }
    }

    /// Create a bloom pass with threshold and intensity.
    pub fn bloom(threshold: f64, intensity: f64) -> Self {
        Self {
            effect: "bloom".to_string(),
            enabled: true,
            params: vec![threshold, intensity],
        }
    }

    /// Create a tone-mapping pass (Reinhard).
    pub fn tone_map_reinhard(exposure: f64) -> Self {
        Self {
            effect: "tonemap_reinhard".to_string(),
            enabled: true,
            params: vec![exposure],
        }
    }

    /// Create a vignette effect.
    pub fn vignette(strength: f64, radius: f64) -> Self {
        Self {
            effect: "vignette".to_string(),
            enabled: true,
            params: vec![strength, radius],
        }
    }

    /// Pipeline descriptor.
    pub fn pipeline(&self) -> RenderPipelineDesc {
        RenderPipelineDesc::post_process(&self.effect)
    }
}

#[wasm_bindgen]
impl PostProcessPass {
    /// Create an FXAA pass from JS.
    #[wasm_bindgen(js_name = "fxaa")]
    pub fn fxaa_js() -> PostProcessPass {
        PostProcessPass::fxaa()
    }

    /// Create a bloom pass from JS.
    #[wasm_bindgen(js_name = "bloom")]
    pub fn bloom_js(threshold: f64, intensity: f64) -> PostProcessPass {
        PostProcessPass::bloom(threshold, intensity)
    }

    /// Create a Reinhard tone-map pass from JS.
    #[wasm_bindgen(js_name = "tone_map_reinhard")]
    pub fn tone_map_reinhard_js(exposure: f64) -> PostProcessPass {
        PostProcessPass::tone_map_reinhard(exposure)
    }

    /// Create a vignette pass from JS.
    #[wasm_bindgen(js_name = "vignette")]
    pub fn vignette_js(strength: f64, radius: f64) -> PostProcessPass {
        PostProcessPass::vignette(strength, radius)
    }

    /// Return the effect name string.
    #[wasm_bindgen(js_name = "get_effect")]
    pub fn get_effect_js(&self) -> String {
        self.effect.clone()
    }

    /// Return the effect params as `Vec<f64>`.
    #[wasm_bindgen(js_name = "get_params")]
    pub fn get_params_js(&self) -> Vec<f64> {
        self.params.clone()
    }
}

/// A post-processing pipeline: an ordered list of effect passes.
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PostProcessPipeline {
    /// Ordered list of effect passes.
    // Vec<PostProcessPass> is not IntoWasmAbi for pub fields — expose via method.
    #[wasm_bindgen(skip)]
    pub passes: Vec<PostProcessPass>,
}

impl PostProcessPipeline {
    /// Create an empty pipeline.
    pub fn new() -> Self {
        Self { passes: Vec::new() }
    }

    /// Add a pass.
    pub fn add_pass(&mut self, pass: PostProcessPass) {
        self.passes.push(pass);
    }

    /// Count of enabled passes.
    pub fn enabled_count(&self) -> usize {
        self.passes.iter().filter(|p| p.enabled).count()
    }

    /// Common preset: FXAA + bloom + tone-mapping.
    pub fn default_hdr() -> Self {
        let mut p = Self::new();
        p.add_pass(PostProcessPass::bloom(0.8, 0.3));
        p.add_pass(PostProcessPass::tone_map_reinhard(1.0));
        p.add_pass(PostProcessPass::fxaa());
        p
    }
}

impl Default for PostProcessPipeline {
    fn default() -> Self {
        Self::new()
    }
}

#[wasm_bindgen]
impl PostProcessPipeline {
    /// Create an empty pipeline from JS.
    #[wasm_bindgen(constructor)]
    pub fn wasm_new() -> PostProcessPipeline {
        PostProcessPipeline::new()
    }

    /// Create the default HDR pipeline from JS.
    #[wasm_bindgen(js_name = "default_hdr")]
    pub fn default_hdr_js() -> PostProcessPipeline {
        PostProcessPipeline::default_hdr()
    }

    /// Add a post-process pass.
    #[wasm_bindgen(js_name = "add_pass")]
    pub fn add_pass_js(&mut self, pass: PostProcessPass) {
        self.add_pass(pass);
    }

    /// Return the count of enabled passes.
    #[wasm_bindgen(js_name = "enabled_count")]
    pub fn enabled_count_js(&self) -> u32 {
        self.enabled_count() as u32
    }

    /// Return total pass count.
    #[wasm_bindgen(js_name = "pass_count")]
    pub fn pass_count_js(&self) -> u32 {
        self.passes.len() as u32
    }

    /// Serialize to a `JsValue` array.
    #[wasm_bindgen(js_name = "to_js")]
    pub fn to_js_js(&self) -> Result<JsValue, JsValue> {
        to_js_value(self)
    }
}

// ===========================================================================
// Render frame graph
// ===========================================================================

/// A simple description of a full render frame graph.
///
/// Lists the shadow, geometry, particle, and post-process passes in order.
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FrameGraph {
    /// Shadow map passes (one per shadow-casting light).
    // Vec<ShadowMapPass> is not IntoWasmAbi for pub fields — expose via method.
    #[wasm_bindgen(skip)]
    pub shadow_passes: Vec<ShadowMapPass>,
    /// G-buffer layout (for deferred shading).
    /// Skipped: `GBufferLayout` is not `Copy`; use the `gbuffer_js()` getter.
    #[wasm_bindgen(skip)]
    pub gbuffer: GBufferLayout,
    /// Particle rendering passes.
    // Vec<ParticleRenderPass> is not IntoWasmAbi for pub fields.
    #[wasm_bindgen(skip)]
    pub particle_passes: Vec<ParticleRenderPass>,
    /// Post-processing pipeline.
    /// Skipped: `PostProcessPipeline` is not `Copy`; use the `post_process_js()` getter.
    #[wasm_bindgen(skip)]
    pub post_process: PostProcessPipeline,
}

impl FrameGraph {
    /// Create a minimal forward-rendering frame graph.
    pub fn forward(width: u32, height: u32) -> Self {
        Self {
            shadow_passes: vec![ShadowMapPass::directional(0, [0.3, -0.8, 0.5], 20.0)],
            gbuffer: GBufferLayout::new(width, height),
            particle_passes: Vec::new(),
            post_process: PostProcessPipeline::default_hdr(),
        }
    }

    /// Create a deferred-rendering frame graph with particles.
    pub fn deferred(width: u32, height: u32) -> Self {
        let mut fg = Self::forward(width, height);
        fg.gbuffer.include_velocity = true;
        fg.particle_passes.push(ParticleRenderPass::additive(65536));
        fg
    }

    /// Total number of render passes.
    pub fn total_pass_count(&self) -> usize {
        self.shadow_passes.len()
            + 2 // geometry + lighting
            + self.particle_passes.len()
            + self.post_process.enabled_count()
    }
}

#[wasm_bindgen]
impl FrameGraph {
    /// Create a forward-rendering frame graph from JS.
    #[wasm_bindgen(js_name = "forward")]
    pub fn forward_js(width: u32, height: u32) -> FrameGraph {
        FrameGraph::forward(width, height)
    }

    /// Create a deferred-rendering frame graph from JS.
    #[wasm_bindgen(js_name = "deferred")]
    pub fn deferred_js(width: u32, height: u32) -> FrameGraph {
        FrameGraph::deferred(width, height)
    }

    /// Return the total pass count.
    #[wasm_bindgen(js_name = "total_pass_count")]
    pub fn total_pass_count_js(&self) -> u32 {
        self.total_pass_count() as u32
    }

    /// Return the shadow pass count.
    #[wasm_bindgen(js_name = "shadow_pass_count")]
    pub fn shadow_pass_count_js(&self) -> u32 {
        self.shadow_passes.len() as u32
    }

    /// Return the particle pass count.
    #[wasm_bindgen(js_name = "particle_pass_count")]
    pub fn particle_pass_count_js(&self) -> u32 {
        self.particle_passes.len() as u32
    }

    /// Serialize to a `JsValue` object.
    #[wasm_bindgen(js_name = "to_js")]
    pub fn to_js_js(&self) -> Result<JsValue, JsValue> {
        to_js_value(self)
    }
}

// ===========================================================================
// Tests
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pipeline_opaque_defaults() {
        let p = RenderPipelineDesc::opaque("test");
        assert!(p.depth_test);
        assert!(p.depth_write);
        assert_eq!(p.blend_mode, BlendMode::Opaque);
        assert_eq!(p.sample_count, 1);
    }

    #[test]
    fn test_pipeline_transparent_no_depth_write() {
        let p = RenderPipelineDesc::transparent("alpha");
        assert!(!p.depth_write);
        assert_eq!(p.blend_mode, BlendMode::AlphaBlend);
    }

    #[test]
    fn test_pipeline_shadow_map_no_color_targets() {
        let p = RenderPipelineDesc::shadow_map("shadow");
        assert_eq!(p.color_target_count, 0);
    }

    #[test]
    fn test_pipeline_post_process_no_depth() {
        let p = RenderPipelineDesc::post_process("pp");
        assert!(!p.depth_test);
        assert!(!p.depth_write);
    }

    #[test]
    fn test_pipeline_with_msaa() {
        let p = RenderPipelineDesc::opaque("msaa").with_msaa(4);
        assert_eq!(p.sample_count, 4);
    }

    #[test]
    fn test_pipeline_with_color_targets() {
        let p = RenderPipelineDesc::opaque("mrt").with_color_targets(4);
        assert_eq!(p.color_target_count, 4);
    }

    #[test]
    fn test_pipeline_msaa_clamped() {
        let p = RenderPipelineDesc::opaque("msaa").with_msaa(16);
        assert_eq!(p.sample_count, 8); // clamped to 8
    }

    #[test]
    fn test_instance_data_flat14_length() {
        let inst = InstanceData::new([1.0, 2.0, 3.0], [0.0, 0.0, 0.0, 1.0]);
        let flat = inst.to_flat14();
        assert_eq!(flat.len(), 14);
    }

    #[test]
    fn test_instance_data_pack_flat() {
        let instances = vec![
            InstanceData::new([0.0; 3], [0.0, 0.0, 0.0, 1.0]),
            InstanceData::new([1.0, 0.0, 0.0], [0.0, 0.0, 0.0, 1.0]),
        ];
        let flat = InstanceData::pack_flat(&instances);
        assert_eq!(flat.len(), 28); // 2 × 14
    }

    #[test]
    fn test_instance_data_position_in_flat() {
        let inst = InstanceData::new([5.0, 6.0, 7.0], [0.0, 0.0, 0.0, 1.0]);
        let flat = inst.to_flat14();
        assert!((flat[0] - 5.0).abs() < 1e-10);
        assert!((flat[1] - 6.0).abs() < 1e-10);
        assert!((flat[2] - 7.0).abs() < 1e-10);
    }

    #[test]
    fn test_shadow_map_pass_directional() {
        let p = ShadowMapPass::directional(0, [0.0, -1.0, 0.0], 10.0);
        assert_eq!(p.id, 0);
        assert_eq!(p.resolution, 1024);
        assert!((p.light_dir[1] + 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_shadow_map_pass_normalizes_direction() {
        let p = ShadowMapPass::directional(1, [3.0, 4.0, 0.0], 5.0);
        let len = (p.light_dir[0] * p.light_dir[0]
            + p.light_dir[1] * p.light_dir[1]
            + p.light_dir[2] * p.light_dir[2])
            .sqrt();
        assert!((len - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_shadow_map_pass_pipeline() {
        let p = ShadowMapPass::directional(0, [0.0, -1.0, 0.0], 10.0);
        let pl = p.pipeline();
        assert_eq!(pl.color_target_count, 0);
    }

    #[test]
    fn test_gbuffer_render_target_count_basic() {
        let g = GBufferLayout::new(1920, 1080);
        assert_eq!(g.render_target_count(), 3);
    }

    #[test]
    fn test_gbuffer_render_target_count_with_extras() {
        let mut g = GBufferLayout::new(1920, 1080);
        g.include_velocity = true;
        g.include_emissive = true;
        assert_eq!(g.render_target_count(), 5);
    }

    #[test]
    fn test_gbuffer_geometry_pass_pipeline() {
        let g = GBufferLayout::new(800, 600);
        let p = g.geometry_pass_pipeline();
        assert!(p.depth_test);
    }

    #[test]
    fn test_gbuffer_lighting_pass_pipeline() {
        let g = GBufferLayout::new(800, 600);
        let p = g.lighting_pass_pipeline();
        assert!(!p.depth_test);
    }

    #[test]
    fn test_particle_pass_additive() {
        let p = ParticleRenderPass::additive(1000);
        assert_eq!(p.blend_mode, BlendMode::Additive);
        assert!(!p.sort_back_to_front);
    }

    #[test]
    fn test_particle_pass_alpha_blended_sorts() {
        let p = ParticleRenderPass::alpha_blended(1000);
        assert_eq!(p.blend_mode, BlendMode::AlphaBlend);
        assert!(p.sort_back_to_front);
    }

    #[test]
    fn test_particle_sort_by_depth() {
        let positions = vec![[1.0, 0.0, 0.0], [5.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let camera = [0.0, 0.0, 0.0];
        let order = ParticleRenderPass::sort_particles_by_depth(&positions, camera);
        assert_eq!(order[0], 1); // dist 5 is farthest
    }

    #[test]
    fn test_post_process_fxaa() {
        let p = PostProcessPass::fxaa();
        assert!(p.enabled);
        assert_eq!(p.effect, "fxaa");
    }

    #[test]
    fn test_post_process_bloom_params() {
        let p = PostProcessPass::bloom(0.9, 0.5);
        assert!((p.params[0] - 0.9).abs() < 1e-10);
        assert!((p.params[1] - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_post_process_tone_map() {
        let p = PostProcessPass::tone_map_reinhard(1.5);
        assert_eq!(p.effect, "tonemap_reinhard");
        assert!((p.params[0] - 1.5).abs() < 1e-10);
    }

    #[test]
    fn test_post_process_vignette() {
        let p = PostProcessPass::vignette(0.5, 0.8);
        assert_eq!(p.effect, "vignette");
    }

    #[test]
    fn test_post_process_pipeline_enabled_count() {
        let p = PostProcessPipeline::default_hdr();
        assert_eq!(p.enabled_count(), 3);
    }

    #[test]
    fn test_post_process_pipeline_add_pass() {
        let mut p = PostProcessPipeline::new();
        p.add_pass(PostProcessPass::fxaa());
        assert_eq!(p.passes.len(), 1);
    }

    #[test]
    fn test_post_process_pipeline_empty_enabled_count() {
        let p = PostProcessPipeline::new();
        assert_eq!(p.enabled_count(), 0);
    }

    #[test]
    fn test_frame_graph_forward_pass_count() {
        let fg = FrameGraph::forward(1920, 1080);
        let total = fg.total_pass_count();
        // 1 shadow + 2 deferred + 0 particles + 3 post-process = 6
        assert_eq!(total, 6);
    }

    #[test]
    fn test_frame_graph_deferred_has_particles() {
        let fg = FrameGraph::deferred(1920, 1080);
        assert!(!fg.particle_passes.is_empty());
    }

    #[test]
    fn test_frame_graph_deferred_has_velocity() {
        let fg = FrameGraph::deferred(1920, 1080);
        assert!(fg.gbuffer.include_velocity);
    }

    #[test]
    fn test_frame_graph_forward_one_shadow() {
        let fg = FrameGraph::forward(800, 600);
        assert_eq!(fg.shadow_passes.len(), 1);
    }

    #[test]
    fn test_render_pipeline_desc_default() {
        let p = RenderPipelineDesc::default();
        assert_eq!(p.label, "default");
        assert_eq!(p.blend_mode, BlendMode::Opaque);
    }

    #[test]
    fn test_post_process_pipeline_default() {
        let p = PostProcessPipeline::default();
        assert_eq!(p.passes.len(), 0);
    }

    // ---------------------------------------------------------------------------
    // Integration tests (Slice W6)
    // ---------------------------------------------------------------------------

    /// Verify that a FrameGraph serializes and deserializes without data loss.
    #[test]
    fn test_frame_graph_json_roundtrip() {
        let fg = FrameGraph::deferred(1920, 1080);
        let json = serde_json::to_string(&fg).expect("serialize");
        let back: FrameGraph = serde_json::from_str(&json).expect("deserialize");
        assert_eq!(back.shadow_passes.len(), fg.shadow_passes.len());
        assert_eq!(back.particle_passes.len(), fg.particle_passes.len());
        assert_eq!(back.post_process.passes.len(), fg.post_process.passes.len());
    }

    /// Verify that a pipeline chain can be built entirely via builder pattern.
    #[test]
    fn test_pipeline_builder_chain() {
        let p = RenderPipelineDesc::opaque("chain")
            .with_msaa(4)
            .with_color_targets(3);
        assert_eq!(p.sample_count, 4);
        assert_eq!(p.color_target_count, 3);
        assert!(p.depth_test);
        assert!(p.backface_culling);
    }
}
