// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! GPU/wgpu rendering backend for OxiPhysics visualisations.
//!
//! This module provides a feature-gated GPU renderer that uses the
//! [`oxiphysics_gpu::compute::WgpuBackend`] from
//! `oxiphysics-gpu` as its device abstraction.
//!
//! When the `wgpu-backend` feature on `oxiphysics-gpu` is active and a
//! compatible GPU adapter is found, rendering is performed on the GPU.
//! Otherwise [`GpuRenderer`] falls back to the software rasteriser from
//! [`mesh_render`](crate::mesh_render).
//!
//! # Feature flag
//!
//! The GPU renderer is gated behind the `wgpu-renderer` feature on this crate
//! (planned for v0.2.0):
//!
//! ```toml
//! [dependencies]
//! oxiphysics-viz = { features = ["wgpu-renderer"] }
//! ```
//!
//! ## Architecture
//!
//! ```text
//!  GpuRenderer
//!   ├── WgpuBackend (from oxiphysics-gpu)
//!   │     ├── depth_buffer   : WgpuBufferHandle   (f32, w×h)
//!   │     ├── color_buffer   : WgpuBufferHandle   (f32 RGBA, w×h×4)
//!   │     └── scratch        : WgpuBufferHandle   (staging)
//!   ├── ShaderRegistry
//!   │     ├── "vertex_transform"   WGSL vertex + Phong lighting shader
//!   │     ├── "particle_splat"     WGSL particle splatting shader
//!   │     └── "fullscreen_blit"    WGSL tone-map + gamma correction
//!   └── render_target: GpuRenderTarget
//! ```
//!
//! ## Render pipeline
//!
//! 1. **Geometry upload** — mesh vertices / indices written to GPU buffers.
//! 2. **Vertex transform dispatch** — project vertices using MVP matrix,
//!    compute per-vertex Phong diffuse + specular.
//! 3. **Rasterise** — a software-style dispatch fills color_buffer with
//!    Phong-shaded fragments; depth-test is performed against depth_buffer.
//! 4. **Post-process** — optional bloom, ACES tone-map, gamma correction.
//! 5. **Present** — colour buffer downloaded to an [`HdrFramebuffer`].
//!
//! ## Usage
//!
//! ```
//! use oxiphysics_viz::wgpu_renderer::{GpuRenderer, GpuRendererConfig, GpuLight};
//!
//! let config = GpuRendererConfig {
//!     width: 1280, height: 720,
//!     ..GpuRendererConfig::default()
//! };
//! let mut renderer = GpuRenderer::new(config);
//!
//! // Upload a simple triangle
//! let verts = vec![[0.0f32, 1.0, 0.0, 1.0],  // xyzw
//!                  [-1.0, -1.0, 0.0, 1.0],
//!                  [ 1.0, -1.0, 0.0, 1.0]];
//! let indices = vec![0u32, 1, 2];
//! let mesh_id = renderer.upload_mesh(&verts, &indices);
//!
//! // Render one frame
//! let identity = [1.0f32,0.,0.,0., 0.,1.,0.,0., 0.,0.,1.,0., 0.,0.,0.,1.];
//! renderer.begin_frame();
//! renderer.draw_mesh(mesh_id, &identity);
//! renderer.end_frame();
//!
//! // Download result to HDR framebuffer
//! let fb = renderer.to_hdr_framebuffer();
//! assert_eq!(fb.width, 1280);
//! ```

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use crate::hdr_framebuffer::{HdrFramebuffer, ToneMapMethod, ToneMapper};
use oxiphysics_gpu::compute::{WgpuBackend, WgpuBufferHandle};

// ── WGSL shaders ─────────────────────────────────────────────────────────────

/// WGSL source: vertex transform + Phong shading.
///
/// Inputs: vertex buffer (`array<vec4<f32>>` = xyzw positions),
///         index buffer (`array<u32>`),
///         MVP matrix uniform.
/// Output: writes to the colour and depth compute buffers (storage).
pub const WGSL_VERTEX_PHONG: &str = r#"
struct Uniforms {
    mvp:        mat4x4<f32>,
    mv:         mat4x4<f32>,   // model-view (for normals)
    light_dir:  vec3<f32>,
    _pad:       f32,
    light_col:  vec3<f32>,
    ambient:    f32,
    diffuse_col: vec3<f32>,
    specular:   f32,
    shininess:  f32,
    width:      u32,
    height:     u32,
    _pad2:      u32,
}

@group(0) @binding(0) var<storage, read>       vertices : array<vec4<f32>>;
@group(0) @binding(1) var<storage, read>       normals  : array<vec4<f32>>;
@group(0) @binding(2) var<storage, read>       indices  : array<u32>;
@group(0) @binding(3) var<storage, read_write> color_buf: array<vec4<f32>>;
@group(0) @binding(4) var<storage, read_write> depth_buf: array<f32>;
@group(0) @binding(5) var<uniform>             uniforms : Uniforms;

@compute @workgroup_size(64, 1, 1)
fn vertex_phong(@builtin(global_invocation_id) gid: vec3<u32>) {
    let tri_idx = gid.x;
    let n_tris  = arrayLength(&indices) / 3u;
    if tri_idx >= n_tris { return; }

    let i0 = indices[tri_idx * 3u];
    let i1 = indices[tri_idx * 3u + 1u];
    let i2 = indices[tri_idx * 3u + 2u];

    let p0 = uniforms.mvp * vertices[i0];
    let p1 = uniforms.mvp * vertices[i1];
    let p2 = uniforms.mvp * vertices[i2];

    // Flat shading: face normal from cross-product in view space
    let e0 = (uniforms.mv * vertices[i1] - uniforms.mv * vertices[i0]).xyz;
    let e1 = (uniforms.mv * vertices[i2] - uniforms.mv * vertices[i0]).xyz;
    let n  = normalize(cross(e0, e1));

    let diff  = max(dot(n, -normalize(uniforms.light_dir)), 0.0);
    let color = vec4<f32>(
        uniforms.diffuse_col * (uniforms.ambient + diff * uniforms.light_col),
        1.0
    );

    // Rasterise triangle into color/depth buffer (simple bounding-box fill)
    let w = f32(uniforms.width);
    let h = f32(uniforms.height);

    // NDC → pixel (perspective divide)
    let px0 = vec2<f32>((p0.x/p0.w * 0.5 + 0.5) * w, (1.0 - p0.y/p0.w * 0.5 - 0.5) * h);
    let px1 = vec2<f32>((p1.x/p1.w * 0.5 + 0.5) * w, (1.0 - p1.y/p1.w * 0.5 - 0.5) * h);
    let px2 = vec2<f32>((p2.x/p2.w * 0.5 + 0.5) * w, (1.0 - p2.y/p2.w * 0.5 - 0.5) * h);

    let min_x = max(i32(min(min(px0.x, px1.x), px2.x)), 0);
    let max_x = min(i32(max(max(px0.x, px1.x), px2.x)), i32(uniforms.width) - 1);
    let min_y = max(i32(min(min(px0.y, px1.y), px2.y)), 0);
    let max_y = min(i32(max(max(px0.y, px1.y), px2.y)), i32(uniforms.height) - 1);

    // Barycentric rasterisation
    for (var py = min_y; py <= max_y; py++) {
        for (var px = min_x; px <= max_x; px++) {
            let p = vec2<f32>(f32(px) + 0.5, f32(py) + 0.5);
            // Edge functions
            let w0 = (px1.y-px2.y)*(p.x-px2.x) + (px2.x-px1.x)*(p.y-px2.y);
            let w1 = (px2.y-px0.y)*(p.x-px0.x) + (px0.x-px2.x)*(p.y-px0.y);
            let w2 = (px0.y-px1.y)*(p.x-px1.x) + (px1.x-px0.x)*(p.y-px1.y);
            if w0 >= 0.0 && w1 >= 0.0 && w2 >= 0.0 {
                let area_inv = 1.0 / (w0 + w1 + w2);
                let bc0 = w0 * area_inv;
                let bc1 = w1 * area_inv;
                let bc2 = w2 * area_inv;
                let z = bc0 * p0.z/p0.w + bc1 * p1.z/p1.w + bc2 * p2.z/p2.w;
                let pixel_idx = u32(py) * uniforms.width + u32(px);
                if z < depth_buf[pixel_idx] {
                    depth_buf[pixel_idx]  = z;
                    color_buf[pixel_idx]  = color;
                }
            }
        }
    }
}
"#;

/// WGSL source: particle splatting (screen-space circles).
pub const WGSL_PARTICLE_SPLAT: &str = r#"
struct ParticleUniforms {
    mvp: mat4x4<f32>,
    color: vec4<f32>,
    radius_px: f32,
    width: u32,
    height: u32,
    _pad: u32,
}
@group(0) @binding(0) var<storage, read>       positions  : array<vec4<f32>>;
@group(0) @binding(1) var<storage, read_write> color_buf  : array<vec4<f32>>;
@group(0) @binding(2) var<storage, read_write> depth_buf  : array<f32>;
@group(0) @binding(3) var<uniform>             uniforms   : ParticleUniforms;

@compute @workgroup_size(64, 1, 1)
fn particle_splat(@builtin(global_invocation_id) gid: vec3<u32>) {
    let pi = gid.x;
    if pi >= arrayLength(&positions) { return; }
    let clip = uniforms.mvp * positions[pi];
    if clip.w <= 0.0 { return; }
    let ndcx = clip.x / clip.w;
    let ndcy = clip.y / clip.w;
    let cx = (ndcx * 0.5 + 0.5) * f32(uniforms.width);
    let cy = (1.0 - ndcy * 0.5 - 0.5) * f32(uniforms.height);
    let z  = clip.z / clip.w;
    let r  = uniforms.radius_px;
    let x0 = max(i32(cx - r), 0);
    let x1 = min(i32(cx + r), i32(uniforms.width) - 1);
    let y0 = max(i32(cy - r), 0);
    let y1 = min(i32(cy + r), i32(uniforms.height) - 1);
    for (var sy = y0; sy <= y1; sy++) {
        for (var sx = x0; sx <= x1; sx++) {
            let dx = f32(sx) + 0.5 - cx;
            let dy = f32(sy) + 0.5 - cy;
            if dx*dx + dy*dy <= r*r {
                let idx = u32(sy) * uniforms.width + u32(sx);
                if z < depth_buf[idx] {
                    depth_buf[idx] = z;
                    color_buf[idx] = uniforms.color;
                }
            }
        }
    }
}
"#;

/// WGSL source: fullscreen blit with ACES tone-map and gamma correction.
pub const WGSL_FULLSCREEN_BLIT: &str = r#"
@group(0) @binding(0) var<storage, read>       hdr_in : array<vec4<f32>>;
@group(0) @binding(1) var<storage, read_write> ldr_out: array<vec4<f32>>;

fn aces(x: f32) -> f32 {
    return clamp((x*(2.51*x+0.03)) / (x*(2.43*x+0.59)+0.14), 0.0, 1.0);
}

fn linear_to_srgb(v: f32) -> f32 {
    if v <= 0.0031308 { return v * 12.92; }
    return 1.055 * pow(v, 1.0/2.4) - 0.055;
}

@compute @workgroup_size(256, 1, 1)
fn fullscreen_blit(@builtin(global_invocation_id) gid: vec3<u32>) {
    let idx = gid.x;
    if idx >= arrayLength(&hdr_in) { return; }
    let h = hdr_in[idx];
    ldr_out[idx] = vec4<f32>(
        linear_to_srgb(aces(h.x)),
        linear_to_srgb(aces(h.y)),
        linear_to_srgb(aces(h.z)),
        h.w
    );
}
"#;

// ── GpuRendererConfig ─────────────────────────────────────────────────────────

/// Configuration for the GPU renderer.
#[derive(Debug, Clone)]
pub struct GpuRendererConfig {
    /// Framebuffer width in pixels.
    pub width: u32,
    /// Framebuffer height in pixels.
    pub height: u32,
    /// Field-of-view in degrees (vertical).
    pub fov_deg: f32,
    /// Near clip plane distance.
    pub near: f32,
    /// Far clip plane distance.
    pub far: f32,
    /// Whether to apply ACES tone-mapping in the blit pass.
    pub tone_map: bool,
    /// HDR exposure multiplier.
    pub exposure: f32,
    /// Ambient light intensity [0, 1].
    pub ambient: f32,
}

impl Default for GpuRendererConfig {
    fn default() -> Self {
        Self {
            width: 800,
            height: 600,
            fov_deg: 60.0,
            near: 0.1,
            far: 1000.0,
            tone_map: true,
            exposure: 1.0,
            ambient: 0.1,
        }
    }
}

// ── GpuLight ─────────────────────────────────────────────────────────────────

/// A directional light for Phong shading.
#[derive(Debug, Clone, Copy)]
pub struct GpuLight {
    /// Direction *from* the light source (will be normalised).
    pub direction: [f32; 3],
    /// Light colour (linear RGB, each component [0, ∞]).
    pub colour: [f32; 3],
}

impl Default for GpuLight {
    fn default() -> Self {
        Self {
            direction: [-0.5, -1.0, -0.5],
            colour: [1.0, 1.0, 1.0],
        }
    }
}

// ── GpuMesh ───────────────────────────────────────────────────────────────────

/// Uploaded GPU mesh (opaque handle + metadata).
#[derive(Debug, Clone, Copy)]
pub struct GpuMeshHandle {
    /// Index into the internal mesh registry.
    pub id: usize,
    /// Number of triangles.
    pub n_triangles: usize,
}

/// Internal mesh entry — GPU path stores buffer handles; CPU path stores raw data.
enum MeshEntry {
    /// Mesh has been uploaded to the GPU backend.
    Gpu {
        /// Vertex buffer handle (vec4<f32> positions, xyzw).
        vertex_buf: WgpuBufferHandle,
        /// Normal buffer handle (vec4<f32>, xyz normal + w unused).
        normal_buf: WgpuBufferHandle,
        /// Index buffer handle (u32 triangle indices).
        index_buf: WgpuBufferHandle,
        n_triangles: usize,
    },
    /// No GPU is available; mesh data is kept in host memory for software rasterisation.
    Cpu {
        /// Vertex positions as `[x, y, z, w]`.
        positions: Vec<[f32; 4]>,
        /// Triangle indices (triples).
        indices: Vec<u32>,
        /// Per-vertex normals `[nx, ny, nz, 0]`.
        normals: Vec<[f32; 4]>,
        n_triangles: usize,
    },
}

// ── GpuRenderTarget ───────────────────────────────────────────────────────────

/// GPU-side render target (colour + depth buffers).
struct GpuRenderTarget {
    color_buf: WgpuBufferHandle,
    depth_buf: WgpuBufferHandle,
    n_pixels: usize,
}

// ── GpuRenderer ───────────────────────────────────────────────────────────────

/// GPU-accelerated renderer using [`WgpuBackend`].
///
/// When no GPU is available the renderer falls back to a software path using
/// the internal CPU colour buffer.
///
/// See [module-level documentation](self) for pipeline details and examples.
pub struct GpuRenderer {
    /// Configuration.
    pub config: GpuRendererConfig,
    /// Active directional light.
    pub light: GpuLight,
    /// Camera view-projection matrix (column-major, 4×4).
    pub view_proj: [f32; 16],
    /// Camera view matrix (for shading).
    pub view: [f32; 16],

    /// GPU backend (None → CPU fallback).
    backend: Option<WgpuBackend>,
    /// Render target (GPU-side colour + depth).
    target: Option<GpuRenderTarget>,
    /// Uploaded meshes.
    meshes: Vec<MeshEntry>,
    /// CPU-side fallback colour buffer (f32 RGBA, w×h×4).
    cpu_color: Vec<[f32; 4]>,
    /// CPU-side fallback depth buffer.
    cpu_depth: Vec<f32>,
    /// Total frames rendered.
    pub frame_count: u64,
}

impl GpuRenderer {
    /// Create a new renderer.  Attempts GPU initialisation; falls back to CPU.
    pub fn new(config: GpuRendererConfig) -> Self {
        let n = (config.width * config.height) as usize;
        let (backend, target) = match WgpuBackend::try_new() {
            Ok(mut b) => {
                // Register shaders
                b.register_shader("vertex_phong", WGSL_VERTEX_PHONG);
                b.register_shader("particle_splat", WGSL_PARTICLE_SPLAT);
                b.register_shader("fullscreen_blit", WGSL_FULLSCREEN_BLIT);

                // Allocate render target buffers
                let color_buf = b.create_buffer(n * 4); // 4 f64 per pixel (RGBA cast to f64 slots)
                let depth_buf = b.create_buffer(n);
                let tgt = GpuRenderTarget {
                    color_buf,
                    depth_buf,
                    n_pixels: n,
                };
                (Some(b), Some(tgt))
            }
            Err(_) => (None, None),
        };

        let view_proj = identity4();
        let view = identity4();

        Self {
            config: config.clone(),
            light: GpuLight::default(),
            view_proj,
            view,
            backend,
            target,
            meshes: Vec::new(),
            cpu_color: vec![[0.0; 4]; n],
            cpu_depth: vec![1.0; n],
            frame_count: 0,
        }
    }

    /// True if a real GPU backend is active.
    pub fn has_gpu(&self) -> bool {
        self.backend.is_some()
    }

    // ── Camera ────────────────────────────────────────────────────────────────

    /// Set the camera using a look-at transform.
    pub fn look_at(&mut self, eye: [f32; 3], target: [f32; 3], up: [f32; 3]) {
        self.view = look_at_mat(eye, target, up);
        let proj = perspective_mat(
            self.config.fov_deg.to_radians(),
            self.config.width as f32 / self.config.height as f32,
            self.config.near,
            self.config.far,
        );
        self.view_proj = mat4_mul(proj, self.view);
    }

    // ── Mesh upload ───────────────────────────────────────────────────────────

    /// Upload a mesh to GPU memory.
    ///
    /// `vertices` — flat `[x, y, z, w=1]` positions (one vec4 per vertex).
    /// `indices`  — triangle index triples.
    ///
    /// Returns a [`GpuMeshHandle`] for subsequent draw calls.
    pub fn upload_mesh(&mut self, vertices: &[[f32; 4]], indices: &[u32]) -> GpuMeshHandle {
        let n_tri = indices.len() / 3;

        // Compute flat normals (one per vertex, averaged from incident triangles)
        let normals = compute_vertex_normals(vertices, indices);

        let id = self.meshes.len();

        if let Some(b) = &mut self.backend {
            let vert_f64: Vec<f64> = vertices
                .iter()
                .flat_map(|v| v.iter().map(|&x| x as f64))
                .collect();
            let norm_f64: Vec<f64> = normals
                .iter()
                .flat_map(|n| n.iter().map(|&x| x as f64))
                .collect();
            let idx_f64: Vec<f64> = indices.iter().map(|&i| i as f64).collect();

            let vb = b.create_buffer(vert_f64.len());
            let nb = b.create_buffer(norm_f64.len());
            let ib = b.create_buffer(idx_f64.len());
            b.write_buffer(vb, &vert_f64);
            b.write_buffer(nb, &norm_f64);
            b.write_buffer(ib, &idx_f64);

            self.meshes.push(MeshEntry::Gpu {
                vertex_buf: vb,
                normal_buf: nb,
                index_buf: ib,
                n_triangles: n_tri,
            });
        } else {
            // CPU fallback: keep mesh data in host memory for software rasterisation.
            self.meshes.push(MeshEntry::Cpu {
                positions: vertices.to_vec(),
                indices: indices.to_vec(),
                normals,
                n_triangles: n_tri,
            });
        }

        GpuMeshHandle {
            id,
            n_triangles: n_tri,
        }
    }

    // ── Frame lifecycle ───────────────────────────────────────────────────────

    /// Begin a new frame: clear colour and depth buffers.
    pub fn begin_frame(&mut self) {
        let n = (self.config.width * self.config.height) as usize;
        self.cpu_color.fill([0.0; 4]);
        self.cpu_depth.fill(1.0_f32);

        if let (Some(b), Some(tgt)) = (&mut self.backend, &self.target) {
            let zeros = vec![0.0_f64; n * 4];
            let ones = vec![1.0_f64; n];
            b.write_buffer(tgt.color_buf, &zeros);
            b.write_buffer(tgt.depth_buf, &ones);
        }
    }

    /// Submit a draw call for mesh `handle` with model matrix `model` (column-major).
    pub fn draw_mesh(&mut self, handle: GpuMeshHandle, model: &[f32; 16]) {
        if handle.id >= self.meshes.len() {
            return;
        }

        let mvp = mat4_mul(self.view_proj, *model);
        let mv = mat4_mul(self.view, *model);
        let w = self.config.width;
        let h = self.config.height;

        match &self.meshes[handle.id] {
            MeshEntry::Gpu {
                vertex_buf,
                normal_buf,
                index_buf,
                n_triangles,
            } => {
                if let (Some(b), Some(tgt)) = (&mut self.backend, &self.target) {
                    let n_tri_wg = (*n_triangles as u32).div_ceil(64);

                    // Pack MVP + light into a small uniform buffer (f64 slots)
                    let mut uni = Vec::with_capacity(48);
                    for v in &mvp {
                        uni.push(*v as f64);
                    } // 16 f64 (MVP)
                    for v in &mv {
                        uni.push(*v as f64);
                    } // 16 f64 (MV)
                    let ld = normalize3f(self.light.direction);
                    uni.push(ld[0] as f64);
                    uni.push(ld[1] as f64);
                    uni.push(ld[2] as f64);
                    uni.push(0.0); // pad
                    uni.push(self.light.colour[0] as f64);
                    uni.push(self.light.colour[1] as f64);
                    uni.push(self.light.colour[2] as f64);
                    uni.push(self.config.ambient as f64);
                    uni.push(0.8); // diffuse_col R
                    uni.push(0.8); // diffuse_col G
                    uni.push(0.8); // diffuse_col B
                    uni.push(0.5); // specular
                    uni.push(32.0); // shininess
                    uni.push(w as f64);
                    uni.push(h as f64);
                    uni.push(0.0); // pad2

                    let ub = b.create_buffer(uni.len());
                    b.write_buffer(ub, &uni);

                    b.dispatch(
                        "vertex_phong",
                        &[
                            *vertex_buf,
                            *normal_buf,
                            *index_buf,
                            tgt.color_buf,
                            tgt.depth_buf,
                            ub,
                        ],
                        n_tri_wg,
                    );
                }
            }
            MeshEntry::Cpu {
                positions,
                indices,
                normals,
                n_triangles: _,
            } => {
                // Software rasterisation path — no GPU required.
                let positions = positions.clone();
                let indices = indices.clone();
                let normals = normals.clone();
                draw_mesh_cpu(
                    &mut self.cpu_color,
                    &mut self.cpu_depth,
                    w as usize,
                    h as usize,
                    &positions,
                    &indices,
                    &normals,
                    &mvp,
                    &mv,
                    self.light,
                    self.config.ambient,
                );
            }
        }
        let _ = (mvp, mv, w, h);
    }

    /// Submit a particle draw call.
    ///
    /// `positions` — flat `[x, y, z, w=1]` particle positions.
    /// `radius_px` — screen-space radius in pixels.
    /// `colour`    — particle colour (linear RGBA).
    pub fn draw_particles(&mut self, positions: &[[f32; 4]], radius_px: f32, colour: [f32; 4]) {
        if positions.is_empty() {
            return;
        }
        let n = positions.len();

        if let (Some(b), Some(tgt)) = (&mut self.backend, &self.target) {
            let pos_f64: Vec<f64> = positions
                .iter()
                .flat_map(|p| p.iter().map(|&x| x as f64))
                .collect();
            let pb = b.create_buffer(pos_f64.len());
            b.write_buffer(pb, &pos_f64);

            let mut uni: Vec<f64> = Vec::with_capacity(24);
            for v in &self.view_proj {
                uni.push(*v as f64);
            }
            uni.push(colour[0] as f64);
            uni.push(colour[1] as f64);
            uni.push(colour[2] as f64);
            uni.push(colour[3] as f64);
            uni.push(radius_px as f64);
            uni.push(self.config.width as f64);
            uni.push(self.config.height as f64);
            uni.push(0.0);

            let ub = b.create_buffer(uni.len());
            b.write_buffer(ub, &uni);

            let wg = (n as u32).div_ceil(64);
            b.dispatch(
                "particle_splat",
                &[pb, tgt.color_buf, tgt.depth_buf, ub],
                wg,
            );
        }
    }

    /// Finalise the frame: apply tone-mapping blit and download result.
    pub fn end_frame(&mut self) {
        let n = (self.config.width * self.config.height) as usize;

        if let (Some(b), Some(tgt)) = (&mut self.backend, &self.target) {
            if self.config.tone_map {
                let ldr_buf = b.create_buffer(n * 4);
                let wg = (n as u32).div_ceil(256);
                b.dispatch("fullscreen_blit", &[tgt.color_buf, ldr_buf], wg);

                // Download LDR result
                let raw = b.read_buffer(ldr_buf);
                for i in 0..n {
                    self.cpu_color[i] = [
                        raw[i * 4] as f32,
                        raw[i * 4 + 1] as f32,
                        raw[i * 4 + 2] as f32,
                        raw[i * 4 + 3] as f32,
                    ];
                }
            } else {
                let raw = b.read_buffer(tgt.color_buf);
                for i in 0..n {
                    self.cpu_color[i] = [
                        raw[i * 4] as f32,
                        raw[i * 4 + 1] as f32,
                        raw[i * 4 + 2] as f32,
                        raw[i * 4 + 3] as f32,
                    ];
                }
            }
        }

        self.frame_count += 1;
    }

    /// Convert the current frame to an [`HdrFramebuffer`].
    pub fn to_hdr_framebuffer(&self) -> HdrFramebuffer {
        let w = self.config.width as usize;
        let h = self.config.height as usize;
        let mut fb = HdrFramebuffer::new(w, h);
        for i in 0..(w * h) {
            let c = self.cpu_color[i];
            fb.set_pixel(i % w, i / w, c);
        }
        fb
    }

    /// Convert the current frame directly to an 8-bit RGBA byte buffer.
    pub fn to_ldr_bytes(&self) -> Vec<u8> {
        let n = (self.config.width * self.config.height) as usize;
        let tone = ToneMapper {
            method: ToneMapMethod::Aces,
            exposure: self.config.exposure,
            gamma: 2.2,
        };
        let mut out = Vec::with_capacity(n * 4);
        for &px in &self.cpu_color {
            let [r, g, b, a] = tone.map(px);
            out.push((r.clamp(0.0, 1.0) * 255.0) as u8);
            out.push((g.clamp(0.0, 1.0) * 255.0) as u8);
            out.push((b.clamp(0.0, 1.0) * 255.0) as u8);
            out.push((a.clamp(0.0, 1.0) * 255.0) as u8);
        }
        out
    }

    /// Clear the colour and depth buffers to a background colour.
    pub fn clear(&mut self, colour: [f32; 4]) {
        self.cpu_color.fill(colour);
        self.cpu_depth.fill(1.0);
    }

    /// Number of uploaded meshes.
    pub fn mesh_count(&self) -> usize {
        self.meshes.len()
    }
}

// ── Math helpers ──────────────────────────────────────────────────────────────

/// Column-major identity 4×4 matrix.
fn identity4() -> [f32; 16] {
    [
        1., 0., 0., 0., 0., 1., 0., 0., 0., 0., 1., 0., 0., 0., 0., 1.,
    ]
}

/// Multiply two 4×4 column-major matrices: C = A × B.
fn mat4_mul(a: [f32; 16], b: [f32; 16]) -> [f32; 16] {
    let mut c = [0.0_f32; 16];
    for row in 0..4 {
        for col in 0..4 {
            for k in 0..4 {
                c[col * 4 + row] += a[k * 4 + row] * b[col * 4 + k];
            }
        }
    }
    c
}

/// Build a look-at view matrix (right-handed, column-major).
fn look_at_mat(eye: [f32; 3], center: [f32; 3], up: [f32; 3]) -> [f32; 16] {
    let f = normalize3f([center[0] - eye[0], center[1] - eye[1], center[2] - eye[2]]);
    let s = normalize3f(cross3f(f, up));
    let u = cross3f(s, f);
    [
        s[0],
        u[0],
        -f[0],
        0.,
        s[1],
        u[1],
        -f[1],
        0.,
        s[2],
        u[2],
        -f[2],
        0.,
        -dot3f(s, eye),
        -dot3f(u, eye),
        -dot3f(f, eye),
        1.,
    ]
}

/// Perspective projection matrix (column-major, right-handed, −1..1 Z range).
fn perspective_mat(fov_y_rad: f32, aspect: f32, near: f32, far: f32) -> [f32; 16] {
    let f = 1.0 / (fov_y_rad * 0.5).tan();
    let nf = 1.0 / (near - far);
    [
        f / aspect,
        0.,
        0.,
        0.,
        0.,
        f,
        0.,
        0.,
        0.,
        0.,
        (far + near) * nf,
        -1.,
        0.,
        0.,
        2. * far * near * nf,
        0.,
    ]
}

fn normalize3f(v: [f32; 3]) -> [f32; 3] {
    let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt().max(1e-8);
    [v[0] / len, v[1] / len, v[2] / len]
}

fn cross3f(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn dot3f(a: [f32; 3], b: [f32; 3]) -> f32 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

// ── CPU software rasteriser ───────────────────────────────────────────────────

/// Multiply a 4×4 column-major matrix by a vec4.
fn mat4_vec4_mul(m: &[f32; 16], v: [f32; 4]) -> [f32; 4] {
    [
        m[0] * v[0] + m[4] * v[1] + m[8] * v[2] + m[12] * v[3],
        m[1] * v[0] + m[5] * v[1] + m[9] * v[2] + m[13] * v[3],
        m[2] * v[0] + m[6] * v[1] + m[10] * v[2] + m[14] * v[3],
        m[3] * v[0] + m[7] * v[1] + m[11] * v[2] + m[15] * v[3],
    ]
}

/// Software rasteriser: rasterise mesh triangles into `color_buf` / `depth_buf`
/// using Phong shading and a barycentric scan-line fill.
///
/// # Arguments
/// * `color_buf` – linear RGBA framebuffer `[f32; 4]` per pixel, width × height.
/// * `depth_buf` – depth buffer `f32` per pixel (1.0 = far, 0.0 = near).
/// * `w`, `h`    – framebuffer dimensions.
/// * `positions` – per-vertex `[x, y, z, w]`.
/// * `indices`   – triangle index triples.
/// * `normals`   – per-vertex `[nx, ny, nz, 0]` in model space.
/// * `mvp`       – combined model-view-projection (column-major).
/// * `mv`        – model-view matrix (column-major, for lighting in view space).
/// * `light`     – active directional light.
/// * `ambient`   – ambient intensity.
#[allow(clippy::too_many_arguments)]
fn draw_mesh_cpu(
    color_buf: &mut [[f32; 4]],
    depth_buf: &mut [f32],
    w: usize,
    h: usize,
    positions: &[[f32; 4]],
    indices: &[u32],
    normals: &[[f32; 4]],
    mvp: &[f32; 16],
    mv: &[f32; 16],
    light: GpuLight,
    ambient: f32,
) {
    let n_tri = indices.len() / 3;
    let wf = w as f32;
    let hf = h as f32;
    let light_dir_view = {
        let ld = normalize3f(light.direction);
        // Transform light direction into view space (ignore translation: w=0).
        let lv = mat4_vec4_mul(mv, [ld[0], ld[1], ld[2], 0.0]);
        normalize3f([lv[0], lv[1], lv[2]])
    };

    for t in 0..n_tri {
        let i0 = indices[t * 3] as usize;
        let i1 = indices[t * 3 + 1] as usize;
        let i2 = indices[t * 3 + 2] as usize;
        if i0 >= positions.len() || i1 >= positions.len() || i2 >= positions.len() {
            continue;
        }

        // Project each vertex into clip space then NDC → screen.
        let clip = [
            mat4_vec4_mul(mvp, positions[i0]),
            mat4_vec4_mul(mvp, positions[i1]),
            mat4_vec4_mul(mvp, positions[i2]),
        ];
        // Perspective divide
        let ndc: [[f32; 3]; 3] = std::array::from_fn(|k| {
            let w_inv = if clip[k][3].abs() > 1e-8 {
                1.0 / clip[k][3]
            } else {
                0.0
            };
            [clip[k][0] * w_inv, clip[k][1] * w_inv, clip[k][2] * w_inv]
        });
        // Viewport transform: NDC [-1,1] → pixel coords [0, w/h)
        let screen: [[f32; 3]; 3] = std::array::from_fn(|k| {
            [
                (ndc[k][0] * 0.5 + 0.5) * wf,
                (1.0 - (ndc[k][1] * 0.5 + 0.5)) * hf,
                ndc[k][2],
            ]
        });

        // Flat normal in view space (average of vertex normals for this triangle)
        let flat_nv = {
            let sum_n: [f32; 3] = [
                normals[i0][0] + normals[i1][0] + normals[i2][0],
                normals[i0][1] + normals[i1][1] + normals[i2][1],
                normals[i0][2] + normals[i1][2] + normals[i2][2],
            ];
            let n_view = mat4_vec4_mul(mv, [sum_n[0], sum_n[1], sum_n[2], 0.0]);
            normalize3f([n_view[0], n_view[1], n_view[2]])
        };
        let diffuse = dot3f(flat_nv, light_dir_view).max(0.0);
        let shading = ambient + diffuse * (1.0 - ambient);
        let base_col = [
            (light.colour[0] * shading).clamp(0.0, 1.0),
            (light.colour[1] * shading).clamp(0.0, 1.0),
            (light.colour[2] * shading).clamp(0.0, 1.0),
            1.0_f32,
        ];

        // Bounding box of the projected triangle, clamped to viewport
        let min_x = screen
            .iter()
            .map(|s| s[0])
            .fold(f32::INFINITY, f32::min)
            .max(0.0) as usize;
        let max_x = (screen
            .iter()
            .map(|s| s[0])
            .fold(f32::NEG_INFINITY, f32::max)
            .ceil() as usize)
            .min(w.saturating_sub(1));
        let min_y = screen
            .iter()
            .map(|s| s[1])
            .fold(f32::INFINITY, f32::min)
            .max(0.0) as usize;
        let max_y = (screen
            .iter()
            .map(|s| s[1])
            .fold(f32::NEG_INFINITY, f32::max)
            .ceil() as usize)
            .min(h.saturating_sub(1));

        // Precompute edge vectors for barycentric test
        let [x0, y0, _] = screen[0];
        let [x1, y1, _] = screen[1];
        let [x2, y2, _] = screen[2];
        let denom = (y1 - y2) * (x0 - x2) + (x2 - x1) * (y0 - y2);
        if denom.abs() < 1e-8 {
            continue; // degenerate
        }
        let inv_denom = 1.0 / denom;

        for py in min_y..=max_y {
            for px in min_x..=max_x {
                let fx = px as f32 + 0.5;
                let fy = py as f32 + 0.5;
                // Barycentric coordinates
                let lam0 = ((y1 - y2) * (fx - x2) + (x2 - x1) * (fy - y2)) * inv_denom;
                let lam1 = ((y2 - y0) * (fx - x2) + (x0 - x2) * (fy - y2)) * inv_denom;
                let lam2 = 1.0 - lam0 - lam1;
                if lam0 < 0.0 || lam1 < 0.0 || lam2 < 0.0 {
                    continue; // outside triangle
                }
                // Interpolated depth
                let depth = lam0 * screen[0][2] + lam1 * screen[1][2] + lam2 * screen[2][2];
                let pidx = py * w + px;
                if depth < depth_buf[pidx] {
                    depth_buf[pidx] = depth;
                    color_buf[pidx] = base_col;
                }
            }
        }
    }
}

/// Compute per-vertex normals by averaging face normals of incident triangles.
fn compute_vertex_normals(vertices: &[[f32; 4]], indices: &[u32]) -> Vec<[f32; 4]> {
    let n = vertices.len();
    let mut normals = vec![[0.0_f32; 4]; n];
    let triangles = indices.len() / 3;
    for t in 0..triangles {
        let i0 = indices[t * 3] as usize;
        let i1 = indices[t * 3 + 1] as usize;
        let i2 = indices[t * 3 + 2] as usize;
        if i0 >= n || i1 >= n || i2 >= n {
            continue;
        }
        let v0 = vertices[i0];
        let v1 = vertices[i1];
        let v2 = vertices[i2];
        let e0 = [v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]];
        let e1 = [v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]];
        let n3 = cross3f(e0, e1);
        for &i in &[i0, i1, i2] {
            normals[i][0] += n3[0];
            normals[i][1] += n3[1];
            normals[i][2] += n3[2];
        }
    }
    for n3 in &mut normals {
        let nn = normalize3f([n3[0], n3[1], n3[2]]);
        n3[0] = nn[0];
        n3[1] = nn[1];
        n3[2] = nn[2];
        n3[3] = 0.0;
    }
    normals
}

// ── tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn simple_triangle() -> (Vec<[f32; 4]>, Vec<u32>) {
        let verts = vec![
            [0.0_f32, 1.0, -2.0, 1.0],
            [-1.0, -1.0, -2.0, 1.0],
            [1.0, -1.0, -2.0, 1.0],
        ];
        let idx = vec![0u32, 1, 2];
        (verts, idx)
    }

    #[test]
    fn test_renderer_construction() {
        let r = GpuRenderer::new(GpuRendererConfig::default());
        assert_eq!(r.config.width, 800);
        assert_eq!(r.config.height, 600);
        assert_eq!(r.frame_count, 0);
    }

    #[test]
    fn test_upload_mesh_and_frame() {
        let mut r = GpuRenderer::new(GpuRendererConfig::default());
        let (v, i) = simple_triangle();
        let h = r.upload_mesh(&v, &i);
        assert_eq!(h.n_triangles, 1);
        assert_eq!(r.mesh_count(), 1);

        r.look_at([0., 0., 0.], [0., 0., -1.], [0., 1., 0.]);
        r.begin_frame();
        r.draw_mesh(h, &identity4());
        r.end_frame();
        assert_eq!(r.frame_count, 1);
    }

    #[test]
    fn test_to_ldr_bytes() {
        let mut r = GpuRenderer::new(GpuRendererConfig {
            width: 4,
            height: 4,
            ..Default::default()
        });
        r.begin_frame();
        r.end_frame();
        let bytes = r.to_ldr_bytes();
        assert_eq!(bytes.len(), 4 * 4 * 4); // w*h*4 channels
    }

    #[test]
    fn test_to_hdr_framebuffer() {
        let r = GpuRenderer::new(GpuRendererConfig {
            width: 8,
            height: 8,
            ..Default::default()
        });
        let fb = r.to_hdr_framebuffer();
        assert_eq!(fb.width, 8);
        assert_eq!(fb.height, 8);
    }

    #[test]
    fn test_mat4_mul_identity() {
        let id = identity4();
        let result = mat4_mul(id, id);
        for i in 0..16 {
            assert!(
                (result[i] - id[i]).abs() < 1e-6,
                "result[{}]={} id[{}]={}",
                i,
                result[i],
                i,
                id[i]
            );
        }
    }

    #[test]
    fn test_look_at_sets_view() {
        let mut r = GpuRenderer::new(GpuRendererConfig::default());
        r.look_at([0., 0., 5.], [0., 0., 0.], [0., 1., 0.]);
        // View matrix should not be identity
        assert!(
            (r.view[14] - 5.0).abs() < 0.01,
            "look-at should set translation along Z, got view[14]={}",
            r.view[14]
        );
    }

    #[test]
    fn test_vertex_normals_triangle() {
        let verts = vec![[0.0_f32, 0., 0., 1.], [1., 0., 0., 1.], [0., 1., 0., 1.]];
        let idx = vec![0u32, 1, 2];
        let n = compute_vertex_normals(&verts, &idx);
        assert_eq!(n.len(), 3);
        // All normals should point in +Z
        for ni in &n {
            assert!(ni[2] > 0.9, "expected +Z normal, got {:?}", ni);
        }
    }

    // --- F3: MeshEntry::Cpu variant ---

    #[test]
    fn test_upload_mesh_cpu_variant_stored() {
        // Force CPU path by checking a no-backend renderer
        // GpuRenderer::new falls back to CPU when WgpuBackend::try_new() fails.
        // In a test environment there is no GPU, so backend == None.
        let mut r = GpuRenderer::new(GpuRendererConfig {
            width: 8,
            height: 8,
            ..Default::default()
        });
        let (v, i) = simple_triangle();
        let h = r.upload_mesh(&v, &i);
        assert_eq!(h.n_triangles, 1);
        assert_eq!(r.mesh_count(), 1);
        // Verify the stored entry is the Cpu variant
        match &r.meshes[0] {
            MeshEntry::Cpu {
                positions,
                indices,
                normals,
                n_triangles,
            } => {
                assert_eq!(*n_triangles, 1);
                assert_eq!(positions.len(), v.len());
                assert_eq!(indices.len(), i.len());
                assert_eq!(normals.len(), v.len());
            }
            MeshEntry::Gpu { .. } => {
                // In GPU-available environments this is fine — skip the Cpu check.
            }
        }
    }

    #[test]
    fn test_render_cpu_produces_pixels() {
        let mut r = GpuRenderer::new(GpuRendererConfig {
            width: 16,
            height: 16,
            ..Default::default()
        });
        // Only meaningful on CPU fallback
        if r.backend.is_some() {
            return;
        }
        let verts = vec![
            [-1.0_f32, -1.0, 0.5, 1.0],
            [1.0, -1.0, 0.5, 1.0],
            [0.0, 1.0, 0.5, 1.0],
        ];
        let idx = vec![0u32, 1, 2];
        let h = r.upload_mesh(&verts, &idx);

        // Set up an orthographic-ish view that keeps the triangle on screen
        let identity = identity4();
        r.view_proj = identity;
        r.view = identity;

        r.begin_frame();
        r.draw_mesh(h, &identity);
        r.end_frame();

        // At least some pixels should be non-black
        let bytes = r.to_ldr_bytes();
        let non_zero = bytes.iter().any(|&b| b > 0);
        assert!(
            non_zero,
            "CPU rasteriser should produce at least one non-black pixel"
        );
    }
}
