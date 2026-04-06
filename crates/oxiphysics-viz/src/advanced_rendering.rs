// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Physically-Based Rendering (PBR) pipeline data structures.
//!
//! Provides CPU-side representations of all major PBR subsystems:
//! PBR materials, shadow maps (including CSM), deferred shading G-buffer,
//! HDR environment maps, BVH for raytracing, screen-space reflections,
//! global illumination (VXGI + light probes), post-process chains,
//! occlusion culling, and clustered deferred lighting.
//!
//! No GPU library dependencies — all structures are pure data.

#![allow(dead_code)]

// ---------------------------------------------------------------------------
// Math helpers
// ---------------------------------------------------------------------------

#[inline]
fn v3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn v3_norm(a: [f64; 3]) -> f64 {
    v3_dot(a, a).sqrt()
}

#[inline]
fn v3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn v3_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn v3_scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline]
fn v3_normalise(a: [f64; 3]) -> [f64; 3] {
    let n = v3_norm(a);
    if n > 1e-12 {
        v3_scale(a, 1.0 / n)
    } else {
        [0.0; 3]
    }
}

#[inline]
fn v3_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

// ---------------------------------------------------------------------------
// PbrMaterial
// ---------------------------------------------------------------------------

/// A texture handle (index into a texture atlas or GPU texture array).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TextureHandle(pub u32);

impl TextureHandle {
    /// The null / unset texture handle.
    pub const NONE: Self = TextureHandle(u32::MAX);

    /// Return `true` if this handle is the null handle.
    pub fn is_none(self) -> bool {
        self.0 == u32::MAX
    }
}

/// Physically-Based Rendering material parameters.
///
/// Follows the metallic-roughness workflow used by glTF 2.0.
/// All scalar values are in the range \[0, 1\] unless noted.
#[derive(Debug, Clone)]
pub struct PbrMaterial {
    /// Base color (albedo) as linear RGB + alpha.
    pub albedo: [f32; 4],
    /// Metalness factor: 0 = dielectric, 1 = metallic.
    pub metallic: f32,
    /// Perceptual roughness factor.
    pub roughness: f32,
    /// Optional albedo texture (RGBA).
    pub albedo_map: TextureHandle,
    /// Optional metallic-roughness texture (G=roughness, B=metallic per glTF).
    pub metallic_roughness_map: TextureHandle,
    /// Optional tangent-space normal map.
    pub normal_map: TextureHandle,
    /// Normal map intensity scale.
    pub normal_scale: f32,
    /// Optional ambient-occlusion texture.
    pub ao_map: TextureHandle,
    /// AO intensity (0 = no occlusion, 1 = full AO).
    pub ao_strength: f32,
    /// Emissive color (HDR, can exceed 1).
    pub emissive: [f32; 3],
    /// Optional emissive texture.
    pub emissive_map: TextureHandle,
    /// Index of refraction (used in dielectric F0 computation).
    pub ior: f32,
    /// Clearcoat layer intensity (0 = none).
    pub clearcoat: f32,
    /// Clearcoat roughness.
    pub clearcoat_roughness: f32,
}

impl Default for PbrMaterial {
    fn default() -> Self {
        Self {
            albedo: [1.0, 1.0, 1.0, 1.0],
            metallic: 0.0,
            roughness: 0.5,
            albedo_map: TextureHandle::NONE,
            metallic_roughness_map: TextureHandle::NONE,
            normal_map: TextureHandle::NONE,
            normal_scale: 1.0,
            ao_map: TextureHandle::NONE,
            ao_strength: 1.0,
            emissive: [0.0; 3],
            emissive_map: TextureHandle::NONE,
            ior: 1.5,
            clearcoat: 0.0,
            clearcoat_roughness: 0.0,
        }
    }
}

impl PbrMaterial {
    /// Compute dielectric F0 reflectance from IOR: `((ior-1)/(ior+1))^2`.
    pub fn f0_dielectric(&self) -> f32 {
        let r = (self.ior - 1.0) / (self.ior + 1.0);
        r * r
    }

    /// Effective F0 blended between dielectric and metallic (uses albedo).
    pub fn f0_effective(&self) -> [f32; 3] {
        let f0_d = self.f0_dielectric();
        let m = self.metallic;
        [
            f0_d * (1.0 - m) + self.albedo[0] * m,
            f0_d * (1.0 - m) + self.albedo[1] * m,
            f0_d * (1.0 - m) + self.albedo[2] * m,
        ]
    }

    /// Returns `true` if the material is emissive.
    pub fn is_emissive(&self) -> bool {
        self.emissive[0] > 0.0 || self.emissive[1] > 0.0 || self.emissive[2] > 0.0
    }

    /// Create a metallic material (e.g., polished steel).
    pub fn metallic_preset(albedo: [f32; 3], roughness: f32) -> Self {
        Self {
            albedo: [albedo[0], albedo[1], albedo[2], 1.0],
            metallic: 1.0,
            roughness,
            ..Default::default()
        }
    }

    /// Create a dielectric material (e.g., plastic).
    pub fn dielectric_preset(albedo: [f32; 3], roughness: f32) -> Self {
        Self {
            albedo: [albedo[0], albedo[1], albedo[2], 1.0],
            metallic: 0.0,
            roughness,
            ..Default::default()
        }
    }
}

// ---------------------------------------------------------------------------
// ShadowMap
// ---------------------------------------------------------------------------

/// Shadow map filter type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ShadowFilterType {
    /// Hard shadows — nearest-neighbor lookup.
    Hard,
    /// Percentage Closer Filtering — averages NxN depth samples.
    Pcf,
    /// Variance Shadow Maps — stores depth moments.
    Vsm,
}

/// Shadow map type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ShadowMapKind {
    /// Orthographic projection for directional lights.
    Directional,
    /// Six-face cube map for point lights.
    Point,
    /// Perspective projection for spot lights.
    Spot,
}

/// A single shadow map layer (one frustum/cascade).
#[derive(Debug, Clone)]
pub struct ShadowMapLayer {
    /// Texture handle.
    pub texture: TextureHandle,
    /// Shadow map resolution (width = height = `resolution`).
    pub resolution: u32,
    /// Light-space projection-view matrix (column-major, 4x4 flattened).
    pub light_matrix: [f32; 16],
    /// Near clip distance.
    pub near: f32,
    /// Far clip distance.
    pub far: f32,
    /// Split distance from camera (CSM only).
    pub split_distance: f32,
}

/// Shadow map collection for a single light.
///
/// Supports directional / point / spot light shadows.
/// For directional lights, multiple cascade layers implement CSM.
#[derive(Debug, Clone)]
pub struct ShadowMap {
    /// Shadow map variant.
    pub kind: ShadowMapKind,
    /// Filter method.
    pub filter: ShadowFilterType,
    /// All layers (1 for spot, 6 for point, 1–4 for directional CSM).
    pub layers: Vec<ShadowMapLayer>,
    /// PCF kernel radius in texels.
    pub pcf_kernel_size: u32,
    /// Depth bias to avoid shadow acne.
    pub bias: f32,
    /// Normal offset bias.
    pub normal_bias: f32,
}

impl ShadowMap {
    /// Create a directional-light shadow map with `n` cascades.
    pub fn new_directional(cascades: usize, resolution: u32) -> Self {
        let layers = (0..cascades)
            .map(|_| ShadowMapLayer {
                texture: TextureHandle::NONE,
                resolution,
                light_matrix: [0.0; 16],
                near: 0.1,
                far: 1000.0,
                split_distance: 0.0,
            })
            .collect();
        Self {
            kind: ShadowMapKind::Directional,
            filter: ShadowFilterType::Pcf,
            layers,
            pcf_kernel_size: 3,
            bias: 0.002,
            normal_bias: 0.001,
        }
    }

    /// Create a spot-light shadow map.
    pub fn new_spot(resolution: u32) -> Self {
        let layer = ShadowMapLayer {
            texture: TextureHandle::NONE,
            resolution,
            light_matrix: [0.0; 16],
            near: 0.1,
            far: 500.0,
            split_distance: 0.0,
        };
        Self {
            kind: ShadowMapKind::Spot,
            filter: ShadowFilterType::Pcf,
            layers: vec![layer],
            pcf_kernel_size: 3,
            bias: 0.002,
            normal_bias: 0.001,
        }
    }

    /// Compute CSM split distances using the practical split scheme.
    /// `lambda` blends between logarithmic (1) and uniform (0) splits.
    pub fn compute_csm_splits(near: f32, far: f32, cascades: usize, lambda: f32) -> Vec<f32> {
        let n = cascades as f32;
        (0..=cascades)
            .map(|i| {
                let f = i as f32 / n;
                let c_log = near * (far / near).powf(f);
                let c_uni = near + (far - near) * f;
                lambda * c_log + (1.0 - lambda) * c_uni
            })
            .collect()
    }

    /// Number of cascades / layers.
    pub fn num_layers(&self) -> usize {
        self.layers.len()
    }
}

// ---------------------------------------------------------------------------
// DeferredShading
// ---------------------------------------------------------------------------

/// G-buffer attachment identifiers.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GBufferSlot {
    /// World-space position (RGB) + linear depth (A).
    Position,
    /// World-space normal (RGB, encoded as octahedral).
    Normal,
    /// Albedo (RGB) + unused (A).
    Albedo,
    /// Metallic (R) + roughness (G) + AO (B) + emissive mask (A).
    Material,
    /// Emissive color (RGB HDR).
    Emissive,
    /// Depth-stencil attachment.
    DepthStencil,
}

/// A single G-buffer attachment.
#[derive(Debug, Clone)]
pub struct GBufferAttachment {
    /// Slot type.
    pub slot: GBufferSlot,
    /// Texture handle.
    pub texture: TextureHandle,
    /// Width in pixels.
    pub width: u32,
    /// Height in pixels.
    pub height: u32,
}

/// Deferred shading context holding all G-buffer attachments and light-pass configuration.
///
/// The render pipeline writes geometry data into G-buffer textures in the
/// geometry pass, then a fullscreen quad processes each pixel in the lighting
/// pass using the stored surface data.
#[derive(Debug, Clone)]
pub struct DeferredShading {
    /// G-buffer attachments.
    pub gbuffer: Vec<GBufferAttachment>,
    /// Render width in pixels.
    pub width: u32,
    /// Render height in pixels.
    pub height: u32,
    /// Whether to enable stencil masking for the lighting pass.
    pub stencil_light_volumes: bool,
    /// Maximum number of lights in the lighting pass.
    pub max_lights: usize,
    /// Whether the emissive pass is enabled.
    pub emissive_pass: bool,
    /// Whether to output HDR color buffer.
    pub hdr: bool,
}

impl DeferredShading {
    /// Allocate a full G-buffer set at the given resolution.
    pub fn new(width: u32, height: u32) -> Self {
        let gbuffer = vec![
            GBufferSlot::Position,
            GBufferSlot::Normal,
            GBufferSlot::Albedo,
            GBufferSlot::Material,
            GBufferSlot::Emissive,
            GBufferSlot::DepthStencil,
        ]
        .into_iter()
        .map(|slot| GBufferAttachment {
            slot,
            texture: TextureHandle::NONE,
            width,
            height,
        })
        .collect();
        Self {
            gbuffer,
            width,
            height,
            stencil_light_volumes: true,
            max_lights: 1024,
            emissive_pass: true,
            hdr: true,
        }
    }

    /// Get an attachment by slot type.
    pub fn attachment(&self, slot: GBufferSlot) -> Option<&GBufferAttachment> {
        self.gbuffer.iter().find(|a| a.slot == slot)
    }

    /// Total memory estimate (bytes) for the G-buffer (assumes 16 bytes/pixel for all non-depth).
    pub fn memory_bytes(&self) -> u64 {
        let pixels = self.width as u64 * self.height as u64;
        // 4 × RGBA16F = 8 bytes each, depth = 4 bytes
        pixels * (5 * 8 + 4)
    }
}

// ---------------------------------------------------------------------------
// EnvironmentMap
// ---------------------------------------------------------------------------

/// Environment map face index (for cube maps).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CubeFace {
    /// +X face.
    PosX,
    /// -X face.
    NegX,
    /// +Y face.
    PosY,
    /// -Y face.
    NegY,
    /// +Z face.
    PosZ,
    /// -Z face.
    NegZ,
}

impl CubeFace {
    /// All six cube faces in order.
    pub const ALL: [CubeFace; 6] = [
        CubeFace::PosX,
        CubeFace::NegX,
        CubeFace::PosY,
        CubeFace::NegY,
        CubeFace::PosZ,
        CubeFace::NegZ,
    ];
}

/// An HDR environment map used for image-based lighting (IBL).
///
/// Stores the raw HDR cubemap plus pre-filtered irradiance and specular maps,
/// and a 2D BRDF look-up table for the split-sum approximation.
#[derive(Debug, Clone)]
pub struct EnvironmentMap {
    /// Handle to the source equirectangular HDR texture.
    pub source_hdr: TextureHandle,
    /// Cube map texture (6 faces of the HDR environment).
    pub cubemap: TextureHandle,
    /// Diffuse irradiance cube map (low-frequency, convolved over hemisphere).
    pub irradiance_map: TextureHandle,
    /// Specular pre-filtered cube map (mip levels = roughness levels).
    pub specular_prefiltered: TextureHandle,
    /// Number of roughness mip levels in `specular_prefiltered`.
    pub specular_mip_levels: u32,
    /// BRDF integration look-up table (2D, NdotV × roughness).
    pub brdf_lut: TextureHandle,
    /// Cube map face resolution (pixels).
    pub resolution: u32,
    /// Intensity scale (HDR exposure multiplier).
    pub intensity: f32,
    /// Rotation of the environment (radians, Y-axis).
    pub rotation_y: f32,
}

impl Default for EnvironmentMap {
    fn default() -> Self {
        Self {
            source_hdr: TextureHandle::NONE,
            cubemap: TextureHandle::NONE,
            irradiance_map: TextureHandle::NONE,
            specular_prefiltered: TextureHandle::NONE,
            specular_mip_levels: 5,
            brdf_lut: TextureHandle::NONE,
            resolution: 512,
            intensity: 1.0,
            rotation_y: 0.0,
        }
    }
}

impl EnvironmentMap {
    /// Compute roughness-to-mip-level mapping: `mip = roughness * (mip_levels - 1)`.
    pub fn roughness_to_mip(&self, roughness: f32) -> f32 {
        roughness.clamp(0.0, 1.0) * (self.specular_mip_levels - 1) as f32
    }

    /// Return `true` if all required textures are assigned.
    pub fn is_complete(&self) -> bool {
        !self.cubemap.is_none()
            && !self.irradiance_map.is_none()
            && !self.specular_prefiltered.is_none()
            && !self.brdf_lut.is_none()
    }
}

// ---------------------------------------------------------------------------
// RaytracingBvh
// ---------------------------------------------------------------------------

/// An axis-aligned bounding box for BVH nodes.
#[derive(Debug, Clone, Copy)]
pub struct Aabb {
    /// Minimum corner.
    pub min: [f64; 3],
    /// Maximum corner.
    pub max: [f64; 3],
}

impl Aabb {
    /// Create an AABB from two corners.
    pub fn new(min: [f64; 3], max: [f64; 3]) -> Self {
        Self { min, max }
    }

    /// Expand an AABB to include a point.
    pub fn expand_point(&self, p: [f64; 3]) -> Self {
        Self {
            min: [
                self.min[0].min(p[0]),
                self.min[1].min(p[1]),
                self.min[2].min(p[2]),
            ],
            max: [
                self.max[0].max(p[0]),
                self.max[1].max(p[1]),
                self.max[2].max(p[2]),
            ],
        }
    }

    /// Merge two AABBs.
    pub fn union(&self, other: &Aabb) -> Aabb {
        Aabb {
            min: [
                self.min[0].min(other.min[0]),
                self.min[1].min(other.min[1]),
                self.min[2].min(other.min[2]),
            ],
            max: [
                self.max[0].max(other.max[0]),
                self.max[1].max(other.max[1]),
                self.max[2].max(other.max[2]),
            ],
        }
    }

    /// Surface area (used in SAH cost).
    pub fn surface_area(&self) -> f64 {
        let d = v3_sub(self.max, self.min);
        2.0 * (d[0] * d[1] + d[1] * d[2] + d[2] * d[0])
    }

    /// Test ray-AABB intersection (slab method). Returns `Some(t_near)` or `None`.
    pub fn intersect_ray(&self, origin: [f64; 3], inv_dir: [f64; 3]) -> Option<f64> {
        let t1 = [
            (self.min[0] - origin[0]) * inv_dir[0],
            (self.min[1] - origin[1]) * inv_dir[1],
            (self.min[2] - origin[2]) * inv_dir[2],
        ];
        let t2 = [
            (self.max[0] - origin[0]) * inv_dir[0],
            (self.max[1] - origin[1]) * inv_dir[1],
            (self.max[2] - origin[2]) * inv_dir[2],
        ];
        let t_min = t1[0].min(t2[0]).max(t1[1].min(t2[1])).max(t1[2].min(t2[2]));
        let t_max = t1[0].max(t2[0]).min(t1[1].max(t2[1])).min(t1[2].max(t2[2]));
        if t_max >= t_min.max(0.0) {
            Some(t_min.max(0.0))
        } else {
            None
        }
    }

    /// Centroid of the AABB.
    pub fn centroid(&self) -> [f64; 3] {
        [
            (self.min[0] + self.max[0]) * 0.5,
            (self.min[1] + self.max[1]) * 0.5,
            (self.min[2] + self.max[2]) * 0.5,
        ]
    }
}

/// A BVH node.
#[derive(Debug, Clone)]
pub struct BvhNode {
    /// Bounding box.
    pub bounds: Aabb,
    /// Left child index (leaf: u32::MAX).
    pub left: u32,
    /// Right child index (leaf: u32::MAX).
    pub right: u32,
    /// Primitive offset (leaf only).
    pub prim_offset: u32,
    /// Primitive count (leaf: > 0, interior: 0).
    pub prim_count: u32,
}

impl BvhNode {
    /// Return `true` if this is a leaf node.
    pub fn is_leaf(&self) -> bool {
        self.prim_count > 0
    }
}

/// A BVH (Bounding Volume Hierarchy) for raytracing acceleration.
///
/// Uses the Surface Area Heuristic (SAH) as cost metric for splitting.
/// Supports AABB intersection and closest-hit triangle query.
#[derive(Debug, Clone, Default)]
pub struct RaytracingBvh {
    /// Flattened list of BVH nodes (root = index 0).
    pub nodes: Vec<BvhNode>,
    /// Triangle vertex data: each triangle is 3 vertices × 3 components.
    pub triangles: Vec<[f64; 9]>,
    /// Maximum leaf primitive count before further splitting.
    pub max_leaf_prims: usize,
}

impl RaytracingBvh {
    /// Create an empty BVH.
    pub fn new() -> Self {
        Self {
            nodes: Vec::new(),
            triangles: Vec::new(),
            max_leaf_prims: 4,
        }
    }

    /// Add a triangle (three vertices) to the BVH primitive list.
    pub fn add_triangle(&mut self, v0: [f64; 3], v1: [f64; 3], v2: [f64; 3]) {
        let tri = [
            v0[0], v0[1], v0[2], v1[0], v1[1], v1[2], v2[0], v2[1], v2[2],
        ];
        self.triangles.push(tri);
    }

    /// Compute AABB for a triangle.
    pub fn triangle_aabb(tri: &[f64; 9]) -> Aabb {
        let v0 = [tri[0], tri[1], tri[2]];
        let v1 = [tri[3], tri[4], tri[5]];
        let v2 = [tri[6], tri[7], tri[8]];
        let base = Aabb::new(v0, v0);
        base.expand_point(v1).expand_point(v2)
    }

    /// Möller-Trumbore ray-triangle intersection.
    /// Returns `Some(t)` (parametric distance) or `None` if no hit.
    pub fn ray_triangle_intersect(origin: [f64; 3], dir: [f64; 3], tri: &[f64; 9]) -> Option<f64> {
        let v0 = [tri[0], tri[1], tri[2]];
        let v1 = [tri[3], tri[4], tri[5]];
        let v2 = [tri[6], tri[7], tri[8]];
        let e1 = v3_sub(v1, v0);
        let e2 = v3_sub(v2, v0);
        let h = v3_cross(dir, e2);
        let a = v3_dot(e1, h);
        if a.abs() < 1e-12 {
            return None;
        }
        let f = 1.0 / a;
        let s = v3_sub(origin, v0);
        let u = f * v3_dot(s, h);
        if !(0.0..=1.0).contains(&u) {
            return None;
        }
        let q = v3_cross(s, e1);
        let v = f * v3_dot(dir, q);
        if v < 0.0 || u + v > 1.0 {
            return None;
        }
        let t = f * v3_dot(e2, q);
        if t > 1e-9 { Some(t) } else { None }
    }

    /// Total number of triangles.
    pub fn num_triangles(&self) -> usize {
        self.triangles.len()
    }
}

// ---------------------------------------------------------------------------
// ScreenSpaceReflections
// ---------------------------------------------------------------------------

/// Screen-space reflections (SSR) configuration.
///
/// Implements ray-marching in screen space to approximate planar reflections.
/// Reflections are blended with the specular PBR term based on surface roughness.
#[derive(Debug, Clone)]
pub struct ScreenSpaceReflections {
    /// Maximum number of ray-march iterations.
    pub max_steps: u32,
    /// Initial ray step size (normalized device coordinates).
    pub step_size: f32,
    /// Binary search refinement iterations (for precise hit location).
    pub refinement_steps: u32,
    /// Maximum ray march distance (NDC units).
    pub max_distance: f32,
    /// Fade-out distance near screen edges (0..1).
    pub edge_fade: f32,
    /// Roughness threshold above which SSR is disabled.
    pub max_roughness: f32,
    /// Blend mode: 0 = replace specular, 1 = additive.
    pub blend_mode: u32,
    /// Thickness (depth tolerance) for intersection testing.
    pub thickness: f32,
}

impl Default for ScreenSpaceReflections {
    fn default() -> Self {
        Self {
            max_steps: 64,
            step_size: 0.05,
            refinement_steps: 8,
            max_distance: 5.0,
            edge_fade: 0.1,
            max_roughness: 0.3,
            blend_mode: 0,
            thickness: 0.02,
        }
    }
}

impl ScreenSpaceReflections {
    /// Returns `true` if a surface with `roughness` should receive SSR.
    pub fn applies_to(&self, roughness: f32) -> bool {
        roughness <= self.max_roughness
    }

    /// Compute edge fade weight for normalized screen coords `(x, y)` ∈ \[-1, 1\].
    pub fn edge_weight(&self, x: f32, y: f32) -> f32 {
        let fx = (1.0 - x.abs()).clamp(0.0, self.edge_fade) / self.edge_fade;
        let fy = (1.0 - y.abs()).clamp(0.0, self.edge_fade) / self.edge_fade;
        fx.min(fy)
    }
}

// ---------------------------------------------------------------------------
// GlobalIllumination
// ---------------------------------------------------------------------------

/// A single SH9 (spherical harmonics order 2) probe for diffuse GI.
///
/// Stores 9 RGB coefficients encoding low-frequency radiance from the probe's
/// baked neighbourhood.
#[derive(Debug, Clone)]
pub struct SH9Probe {
    /// 9 × RGB coefficients (column-major: index = coeff * 3 + channel).
    pub coeffs: [f32; 27],
    /// World-space position of the probe.
    pub position: [f64; 3],
    /// Radius of influence.
    pub radius: f64,
    /// Whether the probe has been baked.
    pub baked: bool,
}

impl SH9Probe {
    /// Create an uninitialized probe at `position` with given `radius`.
    pub fn new(position: [f64; 3], radius: f64) -> Self {
        Self {
            coeffs: [0.0; 27],
            position,
            radius,
            baked: false,
        }
    }

    /// Evaluate the SH9 radiance in direction `dir` (normalised), RGB channel `ch`.
    /// Uses the standard Y_l^m basis evaluation for l=0 and l=1.
    pub fn evaluate(&self, dir: [f64; 3], ch: usize) -> f32 {
        // Band 0: Y_0^0 = 1/(2*sqrt(π))
        let y00: f32 = 0.282_095;
        // Band 1: Y_1^{-1}=sqrt(3/(4π))*y, Y_1^0=sqrt(3/(4π))*z, Y_1^1=sqrt(3/(4π))*x
        let y1m1: f32 = 0.488_603 * dir[1] as f32;
        let y10: f32 = 0.488_603 * dir[2] as f32;
        let y11: f32 = 0.488_603 * dir[0] as f32;
        // Band 2 (remaining 5 are zero in simplified probe)
        let c = ch.min(2);
        self.coeffs[c] * y00
            + self.coeffs[3 + c] * y1m1
            + self.coeffs[2 * 3 + c] * y10
            + self.coeffs[3 * 3 + c] * y11
    }
}

/// Voxel Cone Tracing (VXGI) configuration for dynamic global illumination.
#[derive(Debug, Clone)]
pub struct VoxelGrid {
    /// World-space AABB minimum corner.
    pub min: [f64; 3],
    /// World-space AABB maximum corner.
    pub max: [f64; 3],
    /// Resolution of the voxel grid per axis.
    pub resolution: u32,
    /// Handle to the voxel radiance texture (3D, RGBA16F).
    pub texture: TextureHandle,
    /// Number of mip levels (for cone tracing).
    pub mip_levels: u32,
}

impl VoxelGrid {
    /// Create a new empty voxel grid.
    pub fn new(min: [f64; 3], max: [f64; 3], resolution: u32) -> Self {
        let d = v3_sub(max, min);
        let max_dim = d[0].max(d[1]).max(d[2]);
        let mip_levels = ((resolution as f64).log2().floor() as u32).max(1);
        let _ = max_dim; // intentionally unused — just noting it
        Self {
            min,
            max,
            resolution,
            texture: TextureHandle::NONE,
            mip_levels,
        }
    }

    /// Voxel size (world units per voxel).
    pub fn voxel_size(&self) -> [f64; 3] {
        let d = v3_sub(self.max, self.min);
        v3_scale(d, 1.0 / self.resolution as f64)
    }

    /// World-to-voxel coordinate transform.
    pub fn world_to_voxel(&self, p: [f64; 3]) -> [f64; 3] {
        let d = v3_sub(self.max, self.min);
        [
            (p[0] - self.min[0]) / d[0] * self.resolution as f64,
            (p[1] - self.min[1]) / d[1] * self.resolution as f64,
            (p[2] - self.min[2]) / d[2] * self.resolution as f64,
        ]
    }
}

/// Global illumination system combining VXGI and SH9 light probes.
#[derive(Debug, Clone)]
pub struct GlobalIllumination {
    /// Voxel cone-tracing grid.
    pub voxel_grid: VoxelGrid,
    /// SH9 light probes for coarser diffuse GI.
    pub probes: Vec<SH9Probe>,
    /// Number of cone-tracing rays per pixel.
    pub diffuse_cones: u32,
    /// Cone half-angle for diffuse tracing (radians).
    pub cone_aperture: f32,
    /// Maximum cone trace distance (world units).
    pub max_trace_distance: f64,
    /// Whether to use probes as fallback when no VXGI hit.
    pub probe_fallback: bool,
}

impl GlobalIllumination {
    /// Create a default GI system for a scene of given half-extent.
    pub fn new_scene(half_extent: f64) -> Self {
        let e = half_extent;
        let vg = VoxelGrid::new([-e, -e, -e], [e, e, e], 128);
        Self {
            voxel_grid: vg,
            probes: Vec::new(),
            diffuse_cones: 6,
            cone_aperture: std::f32::consts::FRAC_PI_6, // ~30 degrees
            max_trace_distance: half_extent * 2.0,
            probe_fallback: true,
        }
    }

    /// Add a light probe at `position`.
    pub fn add_probe(&mut self, position: [f64; 3], radius: f64) {
        self.probes.push(SH9Probe::new(position, radius));
    }

    /// Find the nearest probe to a world position.
    pub fn nearest_probe(&self, p: [f64; 3]) -> Option<&SH9Probe> {
        self.probes.iter().min_by(|a, b| {
            let da = v3_norm(v3_sub(a.position, p));
            let db = v3_norm(v3_sub(b.position, p));
            da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
        })
    }
}

// ---------------------------------------------------------------------------
// PostProcessChain
// ---------------------------------------------------------------------------

/// Tone-mapping operator.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ToneMapOp {
    /// Linear (no tone mapping).
    Linear,
    /// Reinhard simple.
    Reinhard,
    /// ACES filmic.
    Aces,
    /// Uncharted 2 / Hable.
    Uncharted2,
}

/// Depth-of-field parameters.
#[derive(Debug, Clone)]
pub struct DepthOfFieldParams {
    /// Focus distance (world units).
    pub focus_distance: f32,
    /// Aperture diameter (affects blur radius).
    pub aperture: f32,
    /// Focal length (mm).
    pub focal_length: f32,
    /// Bokeh blur kernel size (pixels).
    pub blur_radius: u32,
}

/// Lens-flare element.
#[derive(Debug, Clone)]
pub struct LensFlareElement {
    /// Position along the flare axis (0 = light, 1 = center).
    pub position: f32,
    /// Scale of the element.
    pub scale: f32,
    /// Tint color (RGBA).
    pub color: [f32; 4],
    /// Whether to use the starburst pattern.
    pub starburst: bool,
}

/// A chain of post-processing effects applied in order.
///
/// Each effect has an enabled flag and per-effect parameters.
/// The chain is applied as a series of fullscreen passes.
#[derive(Debug, Clone)]
pub struct PostProcessChain {
    /// Exposure value (EV, additive in log2 space).  0.0 = no change.
    pub exposure: f32,
    /// Gamma exponent for final output (1.0 = linear, 2.2 = sRGB).
    pub gamma: f32,
    /// Tone-mapping operator.
    pub tonemap: ToneMapOp,
    /// Enable bloom effect.
    pub bloom_enabled: bool,
    /// Bloom threshold (luminance above which pixels bloom).
    pub bloom_threshold: f32,
    /// Bloom intensity (multiplier).
    pub bloom_intensity: f32,
    /// Enable chromatic aberration.
    pub chromatic_aberration: bool,
    /// Chromatic aberration strength (pixel offset).
    pub ca_strength: f32,
    /// Enable lens flare.
    pub lens_flare: bool,
    /// Lens flare elements.
    pub flare_elements: Vec<LensFlareElement>,
    /// Enable depth of field.
    pub dof_enabled: bool,
    /// Depth-of-field parameters.
    pub dof: DepthOfFieldParams,
    /// Enable vignette effect.
    pub vignette: bool,
    /// Vignette intensity.
    pub vignette_strength: f32,
    /// Enable FXAA anti-aliasing.
    pub fxaa: bool,
    /// Film grain strength.
    pub grain_strength: f32,
}

impl Default for PostProcessChain {
    fn default() -> Self {
        Self {
            exposure: 0.0,
            gamma: 2.2,
            tonemap: ToneMapOp::Aces,
            bloom_enabled: true,
            bloom_threshold: 1.0,
            bloom_intensity: 0.04,
            chromatic_aberration: false,
            ca_strength: 0.5,
            lens_flare: false,
            flare_elements: Vec::new(),
            dof_enabled: false,
            dof: DepthOfFieldParams {
                focus_distance: 10.0,
                aperture: 2.8,
                focal_length: 50.0,
                blur_radius: 8,
            },
            vignette: true,
            vignette_strength: 0.4,
            fxaa: true,
            grain_strength: 0.0,
        }
    }
}

impl PostProcessChain {
    /// Compute the linear exposure multiplier: `2^exposure`.
    pub fn exposure_multiplier(&self) -> f32 {
        (2.0_f32).powf(self.exposure)
    }

    /// Apply ACES filmic tone-mapping to a single luminance value.
    pub fn aces_tonemap(x: f32) -> f32 {
        let a = 2.51_f32;
        let b = 0.03_f32;
        let c = 2.43_f32;
        let d = 0.59_f32;
        let e = 0.14_f32;
        ((x * (a * x + b)) / (x * (c * x + d) + e)).clamp(0.0, 1.0)
    }

    /// Apply Reinhard tone-mapping.
    pub fn reinhard_tonemap(x: f32) -> f32 {
        x / (1.0 + x)
    }

    /// Count active effect passes.
    pub fn active_passes(&self) -> usize {
        [
            self.bloom_enabled,
            self.chromatic_aberration,
            self.lens_flare,
            self.dof_enabled,
            self.vignette,
            self.fxaa,
        ]
        .iter()
        .filter(|&&e| e)
        .count()
    }
}

// ---------------------------------------------------------------------------
// OcclusionCulling
// ---------------------------------------------------------------------------

/// An occlusion query entry.
#[derive(Debug, Clone)]
pub struct OcclusionQuery {
    /// Object identifier.
    pub object_id: u64,
    /// Screen-space bounding box `[x_min, y_min, x_max, y_max]` (pixels).
    pub screen_bbox: [i32; 4],
    /// Whether the query result is available.
    pub result_available: bool,
    /// Whether the object passed the occlusion test.
    pub visible: bool,
}

/// A row of the hierarchical Z-buffer.
#[derive(Debug, Clone)]
pub struct HiZLevel {
    /// Resolution (width = height = `resolution`).
    pub resolution: u32,
    /// Texture handle for this mip level.
    pub texture: TextureHandle,
}

/// Occlusion culling system using a hierarchical Z-buffer (Hi-Z).
///
/// Objects are tested against the coarsest mip level first; if they pass,
/// they are refined against finer levels. Surviving objects are rendered.
#[derive(Debug, Clone)]
pub struct OcclusionCulling {
    /// Hierarchical Z-buffer mip levels.
    pub hiz: Vec<HiZLevel>,
    /// Pending occlusion queries.
    pub queries: Vec<OcclusionQuery>,
    /// Maximum number of simultaneous queries.
    pub max_queries: usize,
    /// Whether two-pass occlusion (read last frame, write this frame) is enabled.
    pub two_pass: bool,
    /// Number of objects culled this frame.
    pub culled_count: u64,
    /// Number of objects rendered this frame.
    pub visible_count: u64,
}

impl OcclusionCulling {
    /// Create a Hi-Z pyramid with base resolution and `levels` mip levels.
    pub fn new(base_resolution: u32, levels: u32) -> Self {
        let hiz = (0..levels)
            .map(|i| HiZLevel {
                resolution: (base_resolution >> i).max(1),
                texture: TextureHandle::NONE,
            })
            .collect();
        Self {
            hiz,
            queries: Vec::new(),
            max_queries: 512,
            two_pass: true,
            culled_count: 0,
            visible_count: 0,
        }
    }

    /// Submit an occlusion query for `object_id` with the given screen bounds.
    pub fn submit_query(&mut self, object_id: u64, screen_bbox: [i32; 4]) {
        if self.queries.len() < self.max_queries {
            self.queries.push(OcclusionQuery {
                object_id,
                screen_bbox,
                result_available: false,
                visible: true,
            });
        }
    }

    /// Mark all pending queries as resolved (simulation; GPU would fill results).
    pub fn resolve_queries(&mut self) {
        for q in &mut self.queries {
            q.result_available = true;
        }
    }

    /// Reset counters and query list for next frame.
    pub fn begin_frame(&mut self) {
        self.queries.clear();
        self.culled_count = 0;
        self.visible_count = 0;
    }

    /// Frustum cull an AABB against six planes. Returns `true` if inside frustum.
    /// `planes` is a flat array of 24 floats: 6 planes × 4 components `[nx, ny, nz, d]`.
    pub fn frustum_cull(planes: &[f32; 24], aabb_min: [f32; 3], aabb_max: [f32; 3]) -> bool {
        for i in 0..6 {
            let nx = planes[i * 4];
            let ny = planes[i * 4 + 1];
            let nz = planes[i * 4 + 2];
            let d = planes[i * 4 + 3];
            // Positive vertex (most aligned with normal)
            let px = if nx >= 0.0 { aabb_max[0] } else { aabb_min[0] };
            let py = if ny >= 0.0 { aabb_max[1] } else { aabb_min[1] };
            let pz = if nz >= 0.0 { aabb_max[2] } else { aabb_min[2] };
            if nx * px + ny * py + nz * pz + d < 0.0 {
                return false;
            }
        }
        true
    }
}

// ---------------------------------------------------------------------------
// LightCluster
// ---------------------------------------------------------------------------

/// A point/spot light descriptor for clustered shading.
#[derive(Debug, Clone)]
pub struct ClusterLight {
    /// World-space position.
    pub position: [f32; 3],
    /// Range (radius) of influence.
    pub range: f32,
    /// Linear color × intensity.
    pub color_intensity: [f32; 3],
    /// Spot light inner cone cos angle (1.0 = point light).
    pub spot_cos_inner: f32,
    /// Spot light outer cone cos angle.
    pub spot_cos_outer: f32,
    /// Spot direction (unit vector).
    pub spot_dir: [f32; 3],
    /// Shadow map layer index (-1 = no shadow).
    pub shadow_layer: i32,
}

/// Cluster grid dimensions.
#[derive(Debug, Clone, Copy)]
pub struct ClusterGrid {
    /// Clusters along screen X.
    pub x: u32,
    /// Clusters along screen Y.
    pub y: u32,
    /// Clusters along the depth (Z) axis.
    pub z: u32,
}

/// Clustered deferred lighting system.
///
/// Divides the view frustum into a 3-D grid of clusters. Each cluster stores
/// a list of lights that overlap its bounds. The shading pass reads only the
/// lights assigned to its cluster, reducing per-pixel light iteration.
#[derive(Debug, Clone)]
pub struct LightCluster {
    /// Cluster grid dimensions.
    pub grid: ClusterGrid,
    /// All dynamic lights in the scene.
    pub lights: Vec<ClusterLight>,
    /// Maximum lights per cluster.
    pub max_lights_per_cluster: usize,
    /// Light index list (flat, indexed by cluster light lists).
    pub light_index_list: Vec<u32>,
    /// Cluster light list offsets and counts `[(offset, count)]`.
    pub cluster_list: Vec<(u32, u32)>,
    /// Near clip distance.
    pub near: f32,
    /// Far clip distance.
    pub far: f32,
    /// Total light assignments this frame.
    pub total_assignments: u64,
}

impl LightCluster {
    /// Create a new cluster grid.
    pub fn new(grid: ClusterGrid, near: f32, far: f32) -> Self {
        let n_clusters = (grid.x * grid.y * grid.z) as usize;
        Self {
            grid,
            lights: Vec::new(),
            max_lights_per_cluster: 256,
            light_index_list: Vec::new(),
            cluster_list: vec![(0, 0); n_clusters],
            near,
            far,
            total_assignments: 0,
        }
    }

    /// Total number of clusters.
    pub fn num_clusters(&self) -> usize {
        (self.grid.x * self.grid.y * self.grid.z) as usize
    }

    /// Flat cluster index from 3-D tile coordinates.
    pub fn cluster_index(&self, x: u32, y: u32, z: u32) -> usize {
        (z * self.grid.y * self.grid.x + y * self.grid.x + x) as usize
    }

    /// Compute the exponential Z-slice index for a given linear depth.
    pub fn depth_to_slice(&self, depth: f32) -> u32 {
        let n = self.near;
        let f = self.far;
        let nz = self.grid.z as f32;
        ((depth - n) / (f - n) * nz).clamp(0.0, nz - 1.0) as u32
    }

    /// Add a light to the scene.
    pub fn add_light(&mut self, light: ClusterLight) {
        self.lights.push(light);
    }

    /// Reset cluster assignments for a new frame.
    pub fn begin_frame(&mut self) {
        self.light_index_list.clear();
        for (offset, count) in &mut self.cluster_list {
            *offset = 0;
            *count = 0;
        }
        self.total_assignments = 0;
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ---- PbrMaterial ----

    #[test]
    fn pbr_material_default_is_white_dielectric() {
        let m = PbrMaterial::default();
        assert!((m.albedo[0] - 1.0).abs() < 1e-6);
        assert!(m.metallic.abs() < 1e-6);
        assert!((m.roughness - 0.5).abs() < 1e-6);
    }

    #[test]
    fn pbr_material_f0_dielectric_ior_1_5() {
        let m = PbrMaterial {
            ior: 1.5,
            ..Default::default()
        };
        let expected = ((1.5_f32 - 1.0) / (1.5 + 1.0)).powi(2);
        assert!((m.f0_dielectric() - expected).abs() < 1e-6);
    }

    #[test]
    fn pbr_material_f0_effective_metallic() {
        let m = PbrMaterial {
            albedo: [1.0, 0.0, 0.0, 1.0],
            metallic: 1.0,
            ..Default::default()
        };
        let f0 = m.f0_effective();
        // For fully metallic, f0 = albedo
        assert!((f0[0] - 1.0).abs() < 1e-6);
        assert!(f0[1].abs() < 1e-6);
    }

    #[test]
    fn pbr_material_is_emissive() {
        let mut m = PbrMaterial::default();
        assert!(!m.is_emissive());
        m.emissive = [1.0, 0.0, 0.0];
        assert!(m.is_emissive());
    }

    #[test]
    fn pbr_material_metallic_preset() {
        let m = PbrMaterial::metallic_preset([0.8, 0.8, 0.8], 0.1);
        assert!((m.metallic - 1.0).abs() < 1e-6);
        assert!((m.roughness - 0.1).abs() < 1e-6);
    }

    // ---- TextureHandle ----

    #[test]
    fn texture_handle_none_is_none() {
        assert!(TextureHandle::NONE.is_none());
        assert!(!TextureHandle(0).is_none());
    }

    // ---- ShadowMap ----

    #[test]
    fn shadow_map_csm_4_cascades() {
        let sm = ShadowMap::new_directional(4, 2048);
        assert_eq!(sm.num_layers(), 4);
        assert_eq!(sm.kind, ShadowMapKind::Directional);
    }

    #[test]
    fn shadow_map_spot_one_layer() {
        let sm = ShadowMap::new_spot(1024);
        assert_eq!(sm.num_layers(), 1);
        assert_eq!(sm.kind, ShadowMapKind::Spot);
    }

    #[test]
    fn shadow_map_csm_splits_monotone() {
        let splits = ShadowMap::compute_csm_splits(0.1, 1000.0, 4, 0.5);
        assert_eq!(splits.len(), 5);
        for w in splits.windows(2) {
            assert!(w[1] > w[0], "CSM splits should be monotone increasing");
        }
    }

    // ---- DeferredShading ----

    #[test]
    fn deferred_shading_has_all_slots() {
        let ds = DeferredShading::new(1920, 1080);
        for slot in [
            GBufferSlot::Position,
            GBufferSlot::Normal,
            GBufferSlot::Albedo,
            GBufferSlot::Material,
            GBufferSlot::Emissive,
            GBufferSlot::DepthStencil,
        ] {
            assert!(ds.attachment(slot).is_some(), "missing slot {:?}", slot);
        }
    }

    #[test]
    fn deferred_shading_memory_estimate_positive() {
        let ds = DeferredShading::new(1280, 720);
        assert!(ds.memory_bytes() > 0);
    }

    // ---- EnvironmentMap ----

    #[test]
    fn environment_map_roughness_to_mip() {
        let em = EnvironmentMap {
            specular_mip_levels: 5,
            ..Default::default()
        };
        assert!((em.roughness_to_mip(0.0)).abs() < 1e-6);
        assert!((em.roughness_to_mip(1.0) - 4.0).abs() < 1e-6);
        assert!((em.roughness_to_mip(0.5) - 2.0).abs() < 1e-6);
    }

    #[test]
    fn environment_map_incomplete_by_default() {
        let em = EnvironmentMap::default();
        assert!(!em.is_complete());
    }

    // ---- Aabb ----

    #[test]
    fn aabb_ray_hit() {
        let aabb = Aabb::new([-1.0; 3], [1.0; 3]);
        let origin = [0.0, 0.0, 5.0];
        let inv_dir = [f64::INFINITY, f64::INFINITY, -1.0];
        let result = aabb.intersect_ray(origin, inv_dir);
        assert!(result.is_some(), "ray should hit unit AABB");
        assert!((result.unwrap() - 4.0).abs() < 1e-10);
    }

    #[test]
    fn aabb_ray_miss() {
        let aabb = Aabb::new([0.0; 3], [1.0; 3]);
        let origin = [5.0, 5.0, 5.0];
        let inv_dir = [1.0, 0.0, 0.0]; // misses
        let result = aabb.intersect_ray(origin, inv_dir);
        assert!(result.is_none(), "ray should miss");
    }

    #[test]
    fn aabb_surface_area_unit_cube() {
        let aabb = Aabb::new([-0.5; 3], [0.5; 3]);
        assert!((aabb.surface_area() - 6.0).abs() < 1e-10);
    }

    #[test]
    fn aabb_centroid_unit_cube() {
        let aabb = Aabb::new([0.0; 3], [2.0; 3]);
        let c = aabb.centroid();
        assert!((c[0] - 1.0).abs() < 1e-10);
        assert!((c[1] - 1.0).abs() < 1e-10);
    }

    // ---- RaytracingBvh ----

    #[test]
    fn bvh_triangle_aabb() {
        let tri = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.5, 1.0, 0.0];
        let aabb = RaytracingBvh::triangle_aabb(&tri);
        assert!((aabb.min[0]).abs() < 1e-12);
        assert!((aabb.max[1] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn bvh_ray_triangle_hit() {
        let tri = [-1.0, -1.0, 0.0, 1.0, -1.0, 0.0, 0.0, 1.0, 0.0];
        let origin = [0.0, 0.0, 1.0];
        let dir = [0.0, 0.0, -1.0];
        let result = RaytracingBvh::ray_triangle_intersect(origin, dir, &tri);
        assert!(result.is_some(), "ray should hit triangle");
        assert!((result.unwrap() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn bvh_ray_triangle_miss() {
        let tri = [0.5, 0.5, 0.0, 1.5, 0.5, 0.0, 1.0, 1.5, 0.0];
        let origin = [0.0, 0.0, 1.0];
        let dir = [0.0, 0.0, -1.0];
        let result = RaytracingBvh::ray_triangle_intersect(origin, dir, &tri);
        assert!(result.is_none(), "ray should miss triangle");
    }

    // ---- ScreenSpaceReflections ----

    #[test]
    fn ssr_applies_to_smooth_surface() {
        let ssr = ScreenSpaceReflections::default();
        assert!(ssr.applies_to(0.1));
        assert!(!ssr.applies_to(0.9));
    }

    #[test]
    fn ssr_edge_weight_center() {
        let ssr = ScreenSpaceReflections {
            edge_fade: 0.2,
            ..Default::default()
        };
        let w = ssr.edge_weight(0.0, 0.0);
        assert!(
            (w - 1.0).abs() < 1e-5,
            "center should have weight 1.0, got {w}"
        );
    }

    #[test]
    fn ssr_edge_weight_corner() {
        let ssr = ScreenSpaceReflections {
            edge_fade: 0.2,
            ..Default::default()
        };
        let w = ssr.edge_weight(1.0, 1.0);
        assert!(w.abs() < 1e-5, "corner should have weight ~0, got {w}");
    }

    // ---- GlobalIllumination ----

    #[test]
    fn gi_add_and_find_nearest_probe() {
        let mut gi = GlobalIllumination::new_scene(10.0);
        gi.add_probe([0.0, 0.0, 0.0], 5.0);
        gi.add_probe([8.0, 0.0, 0.0], 5.0);
        let nearest = gi.nearest_probe([1.0, 0.0, 0.0]);
        assert!(nearest.is_some());
        // Nearest to (1,0,0) should be the probe at (0,0,0)
        let pos = nearest.unwrap().position;
        assert!((pos[0]).abs() < 1e-10, "nearest probe should be at origin");
    }

    #[test]
    fn voxel_grid_voxel_size() {
        let vg = VoxelGrid::new([-10.0; 3], [10.0; 3], 100);
        let sz = vg.voxel_size();
        assert!(
            (sz[0] - 0.2).abs() < 1e-10,
            "voxel size should be 0.2, got {}",
            sz[0]
        );
    }

    #[test]
    fn sh9_probe_evaluate_returns_float() {
        let probe = SH9Probe::new([0.0; 3], 5.0);
        let v = probe.evaluate([0.0, 1.0, 0.0], 0);
        // With zero coefficients, result is 0
        assert!(v.abs() < 1e-6);
    }

    // ---- PostProcessChain ----

    #[test]
    fn post_process_exposure_multiplier() {
        let pp = PostProcessChain {
            exposure: 1.0,
            ..Default::default()
        };
        assert!((pp.exposure_multiplier() - 2.0).abs() < 1e-5);
    }

    #[test]
    fn post_process_aces_range() {
        for x in [0.0_f32, 0.5, 1.0, 2.0, 10.0] {
            let y = PostProcessChain::aces_tonemap(x);
            assert!(
                (0.0..=1.0).contains(&y),
                "ACES should output [0,1] for input {x}, got {y}"
            );
        }
    }

    #[test]
    fn post_process_reinhard_monotone() {
        let a = PostProcessChain::reinhard_tonemap(0.5);
        let b = PostProcessChain::reinhard_tonemap(1.0);
        assert!(b > a, "Reinhard should be monotone increasing");
    }

    #[test]
    fn post_process_active_passes_default() {
        let pp = PostProcessChain::default();
        // bloom, vignette, fxaa enabled by default = 3
        assert_eq!(pp.active_passes(), 3);
    }

    // ---- OcclusionCulling ----

    #[test]
    fn occlusion_culling_hiz_levels() {
        let oc = OcclusionCulling::new(1024, 5);
        assert_eq!(oc.hiz.len(), 5);
        assert_eq!(oc.hiz[0].resolution, 1024);
        assert_eq!(oc.hiz[4].resolution, 64);
    }

    #[test]
    fn occlusion_culling_frustum_inside() {
        // All planes: n=(0,0,1), d=0.0 → positive-z half-space; AABB at z=[0.5,1.0] is inside
        let mut planes = [0.0_f32; 24];
        for i in 0..6 {
            planes[i * 4 + 2] = 1.0;
            planes[i * 4 + 3] = 0.0;
        }
        let result = OcclusionCulling::frustum_cull(&planes, [0.0, 0.0, 0.5], [0.1, 0.1, 1.0]);
        assert!(result, "AABB inside all planes should be visible");
    }

    #[test]
    fn occlusion_culling_frustum_outside() {
        // Plane: n=(0,0,1), d=1.0 → all points with z < -1 are outside
        let mut planes = [0.0_f32; 24];
        planes[2] = 1.0;
        planes[3] = 1.0;
        // Other planes: trivially inside (n=0, d=large positive)
        for i in 1..6 {
            planes[i * 4 + 3] = 1000.0;
        }
        let result = OcclusionCulling::frustum_cull(&planes, [-5.0; 3], [-4.0; 3]);
        assert!(!result, "AABB behind plane should be culled");
    }

    // ---- LightCluster ----

    #[test]
    fn light_cluster_num_clusters() {
        let grid = ClusterGrid { x: 16, y: 9, z: 24 };
        let lc = LightCluster::new(grid, 0.1, 1000.0);
        assert_eq!(lc.num_clusters(), 16 * 9 * 24);
    }

    #[test]
    fn light_cluster_depth_to_slice_clamp() {
        let grid = ClusterGrid { x: 16, y: 9, z: 24 };
        let lc = LightCluster::new(grid, 0.1, 1000.0);
        assert_eq!(lc.depth_to_slice(0.0), 0, "near depth → slice 0");
        assert_eq!(lc.depth_to_slice(1001.0), 23, "far depth → last slice");
    }

    #[test]
    fn light_cluster_index() {
        let grid = ClusterGrid { x: 4, y: 4, z: 4 };
        let lc = LightCluster::new(grid, 0.1, 100.0);
        assert_eq!(lc.cluster_index(0, 0, 0), 0);
        assert_eq!(lc.cluster_index(3, 3, 3), 63);
    }

    #[test]
    fn light_cluster_begin_frame_resets() {
        let grid = ClusterGrid { x: 4, y: 4, z: 4 };
        let mut lc = LightCluster::new(grid, 0.1, 100.0);
        lc.total_assignments = 999;
        lc.begin_frame();
        assert_eq!(lc.total_assignments, 0);
    }
}
