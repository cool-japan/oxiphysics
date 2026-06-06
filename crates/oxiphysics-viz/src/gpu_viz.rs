// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! GPU-accelerated visualization utilities (CPU-side data prep for WebGPU/wgpu).
//!
//! This module provides pure CPU-side data structures and algorithms for
//! preparing geometry, instances, particles, and render pipelines ready for
//! upload to a GPU via WebGPU / wgpu. No actual GPU calls are made here.

// ---------------------------------------------------------------------------
// GpuVertex
// ---------------------------------------------------------------------------

/// A single GPU vertex with position, normal and UV coordinates.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GpuVertex {
    /// World-space or object-space position (x, y, z).
    pub position: [f32; 3],
    /// Surface normal (x, y, z), should be unit-length.
    pub normal: [f32; 3],
    /// Texture coordinates (u, v).
    pub uv: [f32; 2],
}

impl GpuVertex {
    /// Construct a new vertex.
    pub fn new(position: [f32; 3], normal: [f32; 3], uv: [f32; 2]) -> Self {
        Self {
            position,
            normal,
            uv,
        }
    }
}

// ---------------------------------------------------------------------------
// GpuMesh
// ---------------------------------------------------------------------------

/// A mesh ready for GPU upload: interleaved vertex buffer and index buffer.
#[derive(Debug, Clone, Default)]
pub struct GpuMesh {
    /// Interleaved vertex data (position + normal + uv).
    pub vertices: Vec<GpuVertex>,
    /// Triangle indices into `vertices`.
    pub indices: Vec<u32>,
}

impl GpuMesh {
    /// Create an empty mesh.
    pub fn new() -> Self {
        Self::default()
    }

    /// Return the number of triangles.
    pub fn triangle_count(&self) -> usize {
        self.indices.len() / 3
    }

    /// Return `true` if the mesh has no geometry.
    pub fn is_empty(&self) -> bool {
        self.vertices.is_empty()
    }

    /// Byte size of the vertex buffer.
    pub fn vertex_buffer_bytes(&self) -> usize {
        self.vertices.len() * std::mem::size_of::<GpuVertex>()
    }

    /// Byte size of the index buffer.
    pub fn index_buffer_bytes(&self) -> usize {
        self.indices.len() * std::mem::size_of::<u32>()
    }
}

// ---------------------------------------------------------------------------
// GpuBufferKind / GpuBuffer
// ---------------------------------------------------------------------------

/// Describes the intended usage of a GPU buffer.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GpuBufferKind {
    /// Per-vertex attribute data.
    Vertex,
    /// Index data for indexed draw calls.
    Index,
    /// Small uniform / constant data (camera, lights, material).
    Uniform,
    /// Large read/write storage buffer (compute shaders, particles).
    Storage,
}

/// Typed GPU buffer descriptor (no actual GPU allocation — CPU description only).
#[derive(Debug, Clone)]
pub struct GpuBuffer {
    /// Human-readable label for debugging.
    pub label: String,
    /// Intended usage of this buffer.
    pub kind: GpuBufferKind,
    /// Raw byte content to be uploaded.
    pub data: Vec<u8>,
}

impl GpuBuffer {
    /// Construct a new buffer descriptor.
    pub fn new(label: impl Into<String>, kind: GpuBufferKind, data: Vec<u8>) -> Self {
        Self {
            label: label.into(),
            kind,
            data,
        }
    }

    /// Byte size of this buffer.
    pub fn size(&self) -> usize {
        self.data.len()
    }
}

// ---------------------------------------------------------------------------
// RenderTopology
// ---------------------------------------------------------------------------

/// Primitive topology for a render pipeline.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RenderTopology {
    /// Each set of 3 vertices forms a triangle.
    TriangleList,
    /// Consecutive vertices share edges to form a triangle strip.
    TriangleStrip,
    /// Each pair of vertices forms a line segment.
    LineList,
    /// Consecutive vertices form a continuous line.
    LineStrip,
    /// Each vertex is a point sprite.
    PointList,
}

// ---------------------------------------------------------------------------
// RenderPipeline
// ---------------------------------------------------------------------------

/// Descriptor for a GPU render pipeline (CPU-side description only).
#[derive(Debug, Clone)]
pub struct RenderPipeline {
    /// Name of the vertex shader entry point (e.g. `"vs_main"`).
    pub vertex_shader: String,
    /// Name of the fragment shader entry point (e.g. `"fs_main"`).
    pub fragment_shader: String,
    /// Primitive topology.
    pub topology: RenderTopology,
    /// Whether to enable depth testing.
    pub depth_test: bool,
    /// Whether to enable alpha blending.
    pub alpha_blend: bool,
    /// Whether to cull back faces.
    pub cull_back_faces: bool,
}

impl RenderPipeline {
    /// Create a standard opaque triangle pipeline.
    pub fn opaque(vertex_shader: impl Into<String>, fragment_shader: impl Into<String>) -> Self {
        Self {
            vertex_shader: vertex_shader.into(),
            fragment_shader: fragment_shader.into(),
            topology: RenderTopology::TriangleList,
            depth_test: true,
            alpha_blend: false,
            cull_back_faces: true,
        }
    }

    /// Create a transparent (alpha-blended) triangle pipeline.
    pub fn transparent(
        vertex_shader: impl Into<String>,
        fragment_shader: impl Into<String>,
    ) -> Self {
        Self {
            alpha_blend: true,
            ..Self::opaque(vertex_shader, fragment_shader)
        }
    }
}

// ---------------------------------------------------------------------------
// UniformBlock
// ---------------------------------------------------------------------------

/// Camera and lighting uniform data block for upload to the GPU.
#[derive(Debug, Clone)]
pub struct UniformBlock {
    /// View matrix as column-major \[f32; 16\].
    pub view: [f32; 16],
    /// Projection matrix as column-major \[f32; 16\].
    pub proj: [f32; 16],
    /// Up to 4 directional light directions (world space).
    pub light_dirs: [[f32; 3]; 4],
    /// Up to 4 directional light colors (RGB).
    pub light_colors: [[f32; 3]; 4],
    /// Number of active lights.
    pub num_lights: u32,
    /// Ambient light intensity.
    pub ambient: f32,
}

impl Default for UniformBlock {
    fn default() -> Self {
        Self {
            view: identity_mat4(),
            proj: identity_mat4(),
            light_dirs: [[0.0, -1.0, 0.0]; 4],
            light_colors: [[1.0, 1.0, 1.0]; 4],
            num_lights: 1,
            ambient: 0.1,
        }
    }
}

/// Returns a column-major identity matrix as \[f32; 16\].
fn identity_mat4() -> [f32; 16] {
    [
        1.0, 0.0, 0.0, 0.0, // col 0
        0.0, 1.0, 0.0, 0.0, // col 1
        0.0, 0.0, 1.0, 0.0, // col 2
        0.0, 0.0, 0.0, 1.0, // col 3
    ]
}

// ---------------------------------------------------------------------------
// DrawCall
// ---------------------------------------------------------------------------

/// A single draw call descriptor: mesh + material + instance transform.
#[derive(Debug, Clone)]
pub struct DrawCall {
    /// Index into a mesh registry.
    pub mesh_id: usize,
    /// Index into a material / pipeline registry.
    pub material_id: usize,
    /// Model matrix (column-major \[f32; 16\]).
    pub transform: [f32; 16],
    /// Number of instances to render (1 for non-instanced).
    pub instance_count: u32,
}

impl DrawCall {
    /// Create a single-instance draw call.
    pub fn single(mesh_id: usize, material_id: usize, transform: [f32; 16]) -> Self {
        Self {
            mesh_id,
            material_id,
            transform,
            instance_count: 1,
        }
    }

    /// Create an instanced draw call.
    pub fn instanced(mesh_id: usize, material_id: usize, transform: [f32; 16], count: u32) -> Self {
        Self {
            mesh_id,
            material_id,
            transform,
            instance_count: count,
        }
    }
}

/// Batch draw calls that share the same material into merged groups.
///
/// Returns a `Vec` of `(material_id, Vec`DrawCall`)` pairs.
pub fn batch_draw_calls(calls: &[DrawCall]) -> Vec<(usize, Vec<DrawCall>)> {
    use std::collections::HashMap;
    let mut groups: HashMap<usize, Vec<DrawCall>> = HashMap::new();
    for call in calls {
        groups
            .entry(call.material_id)
            .or_default()
            .push(call.clone());
    }
    let mut result: Vec<(usize, Vec<DrawCall>)> = groups.into_iter().collect();
    result.sort_by_key(|(mid, _)| *mid);
    result
}

// ---------------------------------------------------------------------------
// InstanceBuffer
// ---------------------------------------------------------------------------

/// Per-instance data for instanced rendering.
#[derive(Debug, Clone, Copy)]
pub struct InstanceData {
    /// Column-major model matrix for this instance.
    pub transform: [f32; 16],
    /// Per-instance RGBA color.
    pub color: [f32; 4],
}

impl InstanceData {
    /// Construct instance data from a transform and color.
    pub fn new(transform: [f32; 16], color: [f32; 4]) -> Self {
        Self { transform, color }
    }
}

/// A buffer of per-instance transforms and colors.
#[derive(Debug, Clone, Default)]
pub struct InstanceBuffer {
    /// All instances in upload order.
    pub instances: Vec<InstanceData>,
}

impl InstanceBuffer {
    /// Create an empty instance buffer.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add an instance.
    pub fn push(&mut self, data: InstanceData) {
        self.instances.push(data);
    }

    /// Number of instances.
    pub fn len(&self) -> usize {
        self.instances.len()
    }

    /// Return `true` if the buffer is empty.
    pub fn is_empty(&self) -> bool {
        self.instances.is_empty()
    }

    /// Byte size of the raw instance data.
    pub fn byte_size(&self) -> usize {
        self.instances.len() * std::mem::size_of::<InstanceData>()
    }
}

// ---------------------------------------------------------------------------
// ParticleGpuSystem
// ---------------------------------------------------------------------------

/// A single GPU particle.
#[derive(Debug, Clone, Copy)]
pub struct GpuParticle {
    /// World-space position.
    pub position: [f32; 3],
    /// Billboard radius in world units.
    pub radius: f32,
    /// RGBA color.
    pub color: [f32; 4],
}

/// GPU billboard particle system: sorts particles back-to-front for transparency.
#[derive(Debug, Clone, Default)]
pub struct ParticleGpuSystem {
    /// All particles.
    pub particles: Vec<GpuParticle>,
}

impl ParticleGpuSystem {
    /// Create an empty system.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a particle.
    pub fn add(&mut self, position: [f32; 3], radius: f32, color: [f32; 4]) {
        self.particles.push(GpuParticle {
            position,
            radius,
            color,
        });
    }

    /// Sort particles back-to-front (furthest first) relative to a camera position.
    ///
    /// This ordering is required for correct alpha-blended transparency.
    pub fn sort_by_depth(&mut self, camera_pos: [f32; 3]) {
        self.particles.sort_by(|a, b| {
            let da = dist_sq(a.position, camera_pos);
            let db = dist_sq(b.position, camera_pos);
            // Back-to-front: furthest first → descending order
            db.partial_cmp(&da).unwrap_or(std::cmp::Ordering::Equal)
        });
    }

    /// Build an [`InstanceBuffer`] from the current (possibly sorted) particle list.
    pub fn build_instance_buffer(&self) -> InstanceBuffer {
        let mut buf = InstanceBuffer::new();
        for p in &self.particles {
            let t = scale_translate_mat4(p.radius, p.position);
            buf.push(InstanceData::new(t, p.color));
        }
        buf
    }
}

/// Squared Euclidean distance between two 3D points.
fn dist_sq(a: [f32; 3], b: [f32; 3]) -> f32 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    dx * dx + dy * dy + dz * dz
}

/// Build a uniform-scale + translation column-major 4x4 matrix.
fn scale_translate_mat4(scale: f32, translate: [f32; 3]) -> [f32; 16] {
    [
        scale,
        0.0,
        0.0,
        0.0,
        0.0,
        scale,
        0.0,
        0.0,
        0.0,
        0.0,
        scale,
        0.0,
        translate[0],
        translate[1],
        translate[2],
        1.0,
    ]
}

// ---------------------------------------------------------------------------
// TransferFunction
// ---------------------------------------------------------------------------

/// A single entry in a transfer function lookup table.
#[derive(Debug, Clone, Copy)]
pub struct TransferEntry {
    /// Scalar data value (e.g. density).
    pub value: f32,
    /// RGBA color at this value.
    pub color: [f32; 4],
}

/// 1D transfer function mapping scalar values → color + opacity.
#[derive(Debug, Clone, Default)]
pub struct TransferFunction {
    /// Control points, sorted by `value`.
    pub entries: Vec<TransferEntry>,
}

impl TransferFunction {
    /// Create an empty transfer function.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a control point. Keeps entries sorted.
    pub fn add(&mut self, value: f32, color: [f32; 4]) {
        // Clamp opacity to [0, 1]
        let mut c = color;
        c[3] = c[3].clamp(0.0, 1.0);
        self.entries.push(TransferEntry { value, color: c });
        self.entries.sort_by(|a, b| {
            a.value
                .partial_cmp(&b.value)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
    }

    /// Linearly interpolate the transfer function at a given scalar `t`.
    ///
    /// Returns `[0,0,0,0]` if there are no entries.
    pub fn sample(&self, t: f32) -> [f32; 4] {
        if self.entries.is_empty() {
            return [0.0; 4];
        }
        if t <= self
            .entries
            .first()
            .expect("collection should not be empty")
            .value
        {
            return self
                .entries
                .first()
                .expect("collection should not be empty")
                .color;
        }
        if t >= self
            .entries
            .last()
            .expect("collection should not be empty")
            .value
        {
            return self
                .entries
                .last()
                .expect("collection should not be empty")
                .color;
        }
        // Find surrounding entries
        let idx = self.entries.partition_point(|e| e.value <= t);
        let lo = &self.entries[idx - 1];
        let hi = &self.entries[idx];
        let range = hi.value - lo.value;
        let alpha = if range.abs() < 1e-10 {
            0.0
        } else {
            (t - lo.value) / range
        };
        let mut out = [0.0f32; 4];
        for (o, (lc, hc)) in out.iter_mut().zip(lo.color.iter().zip(hi.color.iter())) {
            *o = lc + alpha * (hc - lc);
        }
        out
    }
}

// ---------------------------------------------------------------------------
// VolumeRenderer
// ---------------------------------------------------------------------------

/// Parameters for GPU ray marching through a 3D texture volume.
#[derive(Debug, Clone)]
pub struct VolumeRenderer {
    /// 3D texture dimensions (width, height, depth).
    pub dimensions: [u32; 3],
    /// Raw scalar data (one f32 per voxel, row-major).
    pub data: Vec<f32>,
    /// Transfer function mapping scalar → RGBA.
    pub transfer_fn: TransferFunction,
    /// Step size for ray marching (in texture space).
    pub step_size: f32,
    /// Maximum number of ray marching steps.
    pub max_steps: u32,
    /// Whether to use early ray termination when opacity saturates.
    pub early_termination: bool,
}

impl VolumeRenderer {
    /// Create a volume renderer from raw voxel data.
    pub fn new(dimensions: [u32; 3], data: Vec<f32>) -> Self {
        Self {
            dimensions,
            data,
            transfer_fn: TransferFunction::new(),
            step_size: 0.005,
            max_steps: 512,
            early_termination: true,
        }
    }

    /// Return the total number of voxels.
    pub fn voxel_count(&self) -> usize {
        self.dimensions[0] as usize * self.dimensions[1] as usize * self.dimensions[2] as usize
    }

    /// Sample the scalar value at integer voxel coordinates (bounds-checked).
    pub fn sample_voxel(&self, x: u32, y: u32, z: u32) -> f32 {
        let [w, h, _d] = self.dimensions;
        if x >= w || y >= h || z >= _d {
            return 0.0;
        }
        let idx = (z as usize * h as usize + y as usize) * w as usize + x as usize;
        *self.data.get(idx).unwrap_or(&0.0)
    }
}

// ---------------------------------------------------------------------------
// SdfRenderer
// ---------------------------------------------------------------------------

/// Parameters for SDF-based GPU ray marching.
#[derive(Debug, Clone)]
pub struct SdfRenderer {
    /// World-space bounding box minimum corner.
    pub bounds_min: [f32; 3],
    /// World-space bounding box maximum corner.
    pub bounds_max: [f32; 3],
    /// Maximum number of ray marching steps.
    pub max_steps: u32,
    /// Surface hit threshold (|SDF| < epsilon is considered a hit).
    pub epsilon: f32,
    /// Maximum ray travel distance.
    pub max_dist: f32,
}

impl Default for SdfRenderer {
    fn default() -> Self {
        Self {
            bounds_min: [-1.0; 3],
            bounds_max: [1.0; 3],
            max_steps: 256,
            epsilon: 0.001,
            max_dist: 100.0,
        }
    }
}

impl SdfRenderer {
    /// Create a default SDF renderer.
    pub fn new() -> Self {
        Self::default()
    }
}

// ---------------------------------------------------------------------------
// PostProcessPass
// ---------------------------------------------------------------------------

/// Parameters for a bloom post-processing pass.
#[derive(Debug, Clone)]
pub struct BloomParams {
    /// Luminance threshold for bloom contribution.
    pub threshold: f32,
    /// Bloom intensity multiplier.
    pub intensity: f32,
    /// Number of blur iterations.
    pub iterations: u32,
    /// Whether bloom is enabled.
    pub enabled: bool,
}

impl Default for BloomParams {
    fn default() -> Self {
        Self {
            threshold: 1.0,
            intensity: 0.8,
            iterations: 5,
            enabled: true,
        }
    }
}

/// Parameters for screen-space ambient occlusion (SSAO).
#[derive(Debug, Clone)]
pub struct SsaoParams {
    /// Sampling radius in view space.
    pub radius: f32,
    /// Bias to avoid self-occlusion artefacts.
    pub bias: f32,
    /// Number of hemisphere samples.
    pub num_samples: u32,
    /// Whether SSAO is enabled.
    pub enabled: bool,
}

impl Default for SsaoParams {
    fn default() -> Self {
        Self {
            radius: 0.5,
            bias: 0.025,
            num_samples: 16,
            enabled: true,
        }
    }
}

/// Parameters for depth-of-field (DoF) post-processing.
#[derive(Debug, Clone)]
pub struct DofParams {
    /// Distance to the focal plane in view space.
    pub focus_distance: f32,
    /// Half the depth range that is in focus.
    pub focus_range: f32,
    /// Maximum circle-of-confusion radius in pixels.
    pub max_blur_radius: f32,
    /// Whether DoF is enabled.
    pub enabled: bool,
}

impl Default for DofParams {
    fn default() -> Self {
        Self {
            focus_distance: 10.0,
            focus_range: 3.0,
            max_blur_radius: 8.0,
            enabled: false,
        }
    }
}

/// Aggregated post-processing pass parameters.
#[derive(Debug, Clone, Default)]
pub struct PostProcessPass {
    /// Bloom parameters.
    pub bloom: BloomParams,
    /// SSAO parameters.
    pub ssao: SsaoParams,
    /// Depth-of-field parameters.
    pub dof: DofParams,
}

impl PostProcessPass {
    /// Create a post-process pass with all effects disabled.
    pub fn disabled() -> Self {
        Self {
            bloom: BloomParams {
                enabled: false,
                ..Default::default()
            },
            ssao: SsaoParams {
                enabled: false,
                ..Default::default()
            },
            dof: DofParams {
                enabled: false,
                ..Default::default()
            },
        }
    }

    /// Return `true` if any post-processing effect is enabled.
    pub fn any_enabled(&self) -> bool {
        self.bloom.enabled || self.ssao.enabled || self.dof.enabled
    }
}

// ---------------------------------------------------------------------------
// MeshBuilder
// ---------------------------------------------------------------------------

/// Procedural mesh builder — accumulates vertices and indices, can compute normals.
#[derive(Debug, Clone, Default)]
pub struct MeshBuilder {
    positions: Vec<[f32; 3]>,
    normals: Vec<[f32; 3]>,
    uvs: Vec<[f32; 2]>,
    indices: Vec<u32>,
}

impl MeshBuilder {
    /// Create an empty builder.
    pub fn new() -> Self {
        Self::default()
    }

    /// Append a vertex and return its index.
    pub fn add_vertex(&mut self, position: [f32; 3], uv: [f32; 2]) -> u32 {
        let idx = self.positions.len() as u32;
        self.positions.push(position);
        self.normals.push([0.0, 0.0, 0.0]); // placeholder; recomputed later
        self.uvs.push(uv);
        idx
    }

    /// Add a triangle by vertex indices.
    pub fn add_triangle(&mut self, a: u32, b: u32, c: u32) {
        self.indices.push(a);
        self.indices.push(b);
        self.indices.push(c);
    }

    /// Recompute per-vertex normals by averaging the face normals of all adjacent triangles.
    pub fn compute_normals(&mut self) {
        // Reset normals
        for n in &mut self.normals {
            *n = [0.0, 0.0, 0.0];
        }
        // Accumulate face normals
        let nv = self.positions.len();
        let mut accum = vec![[0.0f32; 3]; nv];
        let tris = self.indices.len() / 3;
        for t in 0..tris {
            let ia = self.indices[t * 3] as usize;
            let ib = self.indices[t * 3 + 1] as usize;
            let ic = self.indices[t * 3 + 2] as usize;
            let pa = self.positions[ia];
            let pb = self.positions[ib];
            let pc = self.positions[ic];
            let ab = sub3f(pb, pa);
            let ac = sub3f(pc, pa);
            let face_n = cross3f(ab, ac);
            for k in 0..3 {
                accum[ia][k] += face_n[k];
                accum[ib][k] += face_n[k];
                accum[ic][k] += face_n[k];
            }
        }
        // Normalize
        for (i, n) in self.normals.iter_mut().enumerate() {
            *n = normalize3f(accum[i]);
        }
    }

    /// Build a [`GpuMesh`] from the current data. Calls `compute_normals` first.
    pub fn build(mut self) -> GpuMesh {
        self.compute_normals();
        let vertices: Vec<GpuVertex> = self
            .positions
            .iter()
            .zip(self.normals.iter())
            .zip(self.uvs.iter())
            .map(|((p, n), uv)| GpuVertex::new(*p, *n, *uv))
            .collect();
        GpuMesh {
            vertices,
            indices: self.indices,
        }
    }

    /// Number of vertices currently in the builder.
    pub fn vertex_count(&self) -> usize {
        self.positions.len()
    }

    /// Number of triangles currently in the builder.
    pub fn triangle_count(&self) -> usize {
        self.indices.len() / 3
    }
}

// 3D float helpers used internally
fn sub3f(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
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

fn length3f(a: [f32; 3]) -> f32 {
    dot3f(a, a).sqrt()
}

fn normalize3f(a: [f32; 3]) -> [f32; 3] {
    let len = length3f(a);
    if len < 1e-12 {
        [0.0, 1.0, 0.0] // fallback
    } else {
        [a[0] / len, a[1] / len, a[2] / len]
    }
}

// ---------------------------------------------------------------------------
// LodSystem
// ---------------------------------------------------------------------------

/// LOD level descriptor.
#[derive(Debug, Clone)]
pub struct LodLevel {
    /// Maximum distance (world units) at which this level is active.
    pub max_distance: f32,
    /// Index into a mesh registry.
    pub mesh_id: usize,
}

/// LOD system: selects the appropriate detail level for a given object-camera distance.
#[derive(Debug, Clone, Default)]
pub struct LodSystem {
    /// LOD levels, ordered from closest (highest detail) to farthest (lowest detail).
    pub levels: Vec<LodLevel>,
}

impl LodSystem {
    /// Create an empty LOD system.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add an LOD level. Levels should be added from closest to farthest.
    pub fn add_level(&mut self, max_distance: f32, mesh_id: usize) {
        self.levels.push(LodLevel {
            max_distance,
            mesh_id,
        });
    }

    /// Select the mesh id for a given distance from the camera.
    ///
    /// Returns the last (lowest-detail) level if `distance` exceeds all thresholds.
    /// Returns `None` if no levels have been registered.
    pub fn select(&self, distance: f32) -> Option<usize> {
        for level in &self.levels {
            if distance <= level.max_distance {
                return Some(level.mesh_id);
            }
        }
        // Beyond all levels: use the last (lowest detail)
        self.levels.last().map(|l| l.mesh_id)
    }

    /// Number of LOD levels registered.
    pub fn level_count(&self) -> usize {
        self.levels.len()
    }
}

// ---------------------------------------------------------------------------
// FrustumCull
// ---------------------------------------------------------------------------

/// A single frustum plane (normal + offset in the form n·x + d = 0).
#[derive(Debug, Clone, Copy)]
pub struct FrustumPlane {
    /// Inward-facing normal of the plane.
    pub normal: [f32; 3],
    /// Plane offset: for a point p, `dot(normal, p) + d > 0` means inside.
    pub d: f32,
}

impl FrustumPlane {
    /// Construct a frustum plane from normal and offset.
    pub fn new(normal: [f32; 3], d: f32) -> Self {
        Self {
            normal: normalize3f(normal),
            d,
        }
    }

    /// Signed distance of a point from the plane (positive = inside).
    pub fn distance(&self, point: [f32; 3]) -> f32 {
        dot3f(self.normal, point) + self.d
    }
}

/// An axis-aligned bounding box (AABB) for frustum culling.
#[derive(Debug, Clone, Copy)]
pub struct CullAabb {
    /// Minimum corner.
    pub min: [f32; 3],
    /// Maximum corner.
    pub max: [f32; 3],
}

impl CullAabb {
    /// Construct an AABB from center and half-extents.
    pub fn from_center_half(center: [f32; 3], half: [f32; 3]) -> Self {
        Self {
            min: [
                center[0] - half[0],
                center[1] - half[1],
                center[2] - half[2],
            ],
            max: [
                center[0] + half[0],
                center[1] + half[1],
                center[2] + half[2],
            ],
        }
    }
}

/// Six-plane view frustum for AABB culling.
#[derive(Debug, Clone)]
pub struct FrustumCull {
    /// The six frustum planes (left, right, bottom, top, near, far).
    pub planes: [FrustumPlane; 6],
}

impl FrustumCull {
    /// Construct a frustum from six planes.
    pub fn new(planes: [FrustumPlane; 6]) -> Self {
        Self { planes }
    }

    /// Build a simple symmetric perspective frustum.
    ///
    /// - `near` / `far` — clip distances.
    /// - `fov_y_rad` — vertical field of view in radians.
    /// - `aspect` — width / height.
    /// - `camera_pos` — camera world position (frustum apex).
    /// - `forward` — normalized camera forward direction.
    /// - `up` — normalized camera up direction.
    pub fn from_perspective(
        near: f32,
        far: f32,
        fov_y_rad: f32,
        aspect: f32,
        camera_pos: [f32; 3],
        forward: [f32; 3],
        up: [f32; 3],
    ) -> Self {
        let right = normalize3f(cross3f(forward, up));
        let up_ortho = normalize3f(cross3f(right, forward));

        let half_v = (fov_y_rad * 0.5).tan();
        let half_h = half_v * aspect;

        // Near plane center
        let nc = [
            camera_pos[0] + forward[0] * near,
            camera_pos[1] + forward[1] * near,
            camera_pos[2] + forward[2] * near,
        ];
        // Far plane center
        let fc = [
            camera_pos[0] + forward[0] * far,
            camera_pos[1] + forward[1] * far,
            camera_pos[2] + forward[2] * far,
        ];

        // Helper: plane from normal and a point on the plane
        let plane_from_point_normal = |pt: [f32; 3], n: [f32; 3]| {
            let n = normalize3f(n);
            FrustumPlane {
                normal: n,
                d: -dot3f(n, pt),
            }
        };

        // Near plane: normal = forward
        let near_plane = plane_from_point_normal(nc, forward);
        // Far plane: normal = -forward
        let far_plane = plane_from_point_normal(fc, [-forward[0], -forward[1], -forward[2]]);

        // Left plane
        let left_n = cross3f(
            up_ortho,
            [
                forward[0] - right[0] * half_h,
                forward[1] - right[1] * half_h,
                forward[2] - right[2] * half_h,
            ],
        );
        let left_plane = plane_from_point_normal(camera_pos, left_n);

        // Right plane
        let right_n = cross3f(
            [
                forward[0] + right[0] * half_h,
                forward[1] + right[1] * half_h,
                forward[2] + right[2] * half_h,
            ],
            up_ortho,
        );
        let right_plane = plane_from_point_normal(camera_pos, right_n);

        // Bottom plane
        let bottom_n = cross3f(
            right,
            [
                forward[0] - up_ortho[0] * half_v,
                forward[1] - up_ortho[1] * half_v,
                forward[2] - up_ortho[2] * half_v,
            ],
        );
        let bottom_plane = plane_from_point_normal(camera_pos, bottom_n);

        // Top plane
        let top_n = cross3f(
            [
                forward[0] + up_ortho[0] * half_v,
                forward[1] + up_ortho[1] * half_v,
                forward[2] + up_ortho[2] * half_v,
            ],
            right,
        );
        let top_plane = plane_from_point_normal(camera_pos, top_n);

        Self::new([
            near_plane,
            far_plane,
            left_plane,
            right_plane,
            bottom_plane,
            top_plane,
        ])
    }

    /// Test an AABB against the frustum.
    ///
    /// Returns `true` if the AABB is *entirely outside* any single plane → should be culled.
    pub fn cull_aabb(&self, aabb: &CullAabb) -> bool {
        for plane in &self.planes {
            // Find the "positive vertex" — the corner of the AABB most in the direction of the normal
            let px = if plane.normal[0] >= 0.0 {
                aabb.max[0]
            } else {
                aabb.min[0]
            };
            let py = if plane.normal[1] >= 0.0 {
                aabb.max[1]
            } else {
                aabb.min[1]
            };
            let pz = if plane.normal[2] >= 0.0 {
                aabb.max[2]
            } else {
                aabb.min[2]
            };
            if plane.distance([px, py, pz]) < 0.0 {
                // The most "inside" corner is still outside this plane → fully culled
                return true;
            }
        }
        false
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ---- GpuVertex ----

    #[test]
    fn test_gpu_vertex_construction() {
        let v = GpuVertex::new([1.0, 2.0, 3.0], [0.0, 1.0, 0.0], [0.5, 0.5]);
        assert_eq!(v.position, [1.0, 2.0, 3.0]);
        assert_eq!(v.normal, [0.0, 1.0, 0.0]);
        assert_eq!(v.uv, [0.5, 0.5]);
    }

    // ---- GpuMesh ----

    #[test]
    fn test_gpu_mesh_empty() {
        let m = GpuMesh::new();
        assert!(m.is_empty());
        assert_eq!(m.triangle_count(), 0);
        assert_eq!(m.vertex_buffer_bytes(), 0);
        assert_eq!(m.index_buffer_bytes(), 0);
    }

    #[test]
    fn test_gpu_mesh_triangle_count() {
        let mesh = GpuMesh {
            vertices: vec![
                GpuVertex::new([0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0]),
                GpuVertex::new([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 0.0]),
                GpuVertex::new([0.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 1.0]),
            ],
            indices: vec![0, 1, 2],
        };
        assert_eq!(mesh.triangle_count(), 1);
        assert!(!mesh.is_empty());
    }

    #[test]
    fn test_gpu_mesh_byte_sizes() {
        let v = GpuVertex::new([0.0; 3], [0.0; 3], [0.0; 2]);
        let mesh = GpuMesh {
            vertices: vec![v; 4],
            indices: vec![0, 1, 2, 0, 2, 3],
        };
        assert_eq!(
            mesh.vertex_buffer_bytes(),
            4 * std::mem::size_of::<GpuVertex>()
        );
        assert_eq!(mesh.index_buffer_bytes(), 6 * 4);
    }

    // ---- GpuBuffer ----

    #[test]
    fn test_gpu_buffer_size() {
        let buf = GpuBuffer::new("test", GpuBufferKind::Vertex, vec![0u8; 128]);
        assert_eq!(buf.size(), 128);
    }

    #[test]
    fn test_gpu_buffer_kinds() {
        for kind in [
            GpuBufferKind::Vertex,
            GpuBufferKind::Index,
            GpuBufferKind::Uniform,
            GpuBufferKind::Storage,
        ] {
            let b = GpuBuffer::new("b", kind, vec![]);
            assert_eq!(b.kind, kind);
        }
    }

    // ---- RenderPipeline ----

    #[test]
    fn test_render_pipeline_opaque_defaults() {
        let p = RenderPipeline::opaque("vs_main", "fs_main");
        assert!(!p.alpha_blend);
        assert!(p.depth_test);
        assert!(p.cull_back_faces);
        assert_eq!(p.topology, RenderTopology::TriangleList);
    }

    #[test]
    fn test_render_pipeline_transparent() {
        let p = RenderPipeline::transparent("vs_main", "fs_alpha");
        assert!(p.alpha_blend);
    }

    // ---- DrawCall / batching ----

    #[test]
    fn test_draw_call_single() {
        let dc = DrawCall::single(0, 1, identity_mat4());
        assert_eq!(dc.instance_count, 1);
        assert_eq!(dc.mesh_id, 0);
        assert_eq!(dc.material_id, 1);
    }

    #[test]
    fn test_draw_call_instanced() {
        let dc = DrawCall::instanced(2, 3, identity_mat4(), 64);
        assert_eq!(dc.instance_count, 64);
    }

    #[test]
    fn test_batch_draw_calls_groups_by_material() {
        let calls = vec![
            DrawCall::single(0, 1, identity_mat4()),
            DrawCall::single(1, 2, identity_mat4()),
            DrawCall::single(2, 1, identity_mat4()),
            DrawCall::single(3, 2, identity_mat4()),
        ];
        let batches = batch_draw_calls(&calls);
        // Two materials: 1 and 2
        assert_eq!(batches.len(), 2);
        // Material 1 has 2 calls, material 2 has 2 calls
        for (_mid, calls) in &batches {
            assert_eq!(calls.len(), 2);
        }
    }

    #[test]
    fn test_batch_draw_calls_single_material() {
        let calls: Vec<DrawCall> = (0..5)
            .map(|i| DrawCall::single(i, 0, identity_mat4()))
            .collect();
        let batches = batch_draw_calls(&calls);
        assert_eq!(batches.len(), 1);
        assert_eq!(batches[0].0, 0);
        assert_eq!(batches[0].1.len(), 5);
    }

    // ---- InstanceBuffer ----

    #[test]
    fn test_instance_buffer_push_and_len() {
        let mut buf = InstanceBuffer::new();
        assert!(buf.is_empty());
        buf.push(InstanceData::new(identity_mat4(), [1.0, 0.0, 0.0, 1.0]));
        buf.push(InstanceData::new(identity_mat4(), [0.0, 1.0, 0.0, 1.0]));
        assert_eq!(buf.len(), 2);
        assert!(!buf.is_empty());
    }

    #[test]
    fn test_instance_buffer_transform_layout() {
        let t = scale_translate_mat4(2.0, [3.0, 4.0, 5.0]);
        let inst = InstanceData::new(t, [1.0; 4]);
        // Scale on diagonal
        assert!((inst.transform[0] - 2.0).abs() < 1e-6);
        assert!((inst.transform[5] - 2.0).abs() < 1e-6);
        assert!((inst.transform[10] - 2.0).abs() < 1e-6);
        // Translation in last column (col 3 in column-major: indices 12,13,14)
        assert!((inst.transform[12] - 3.0).abs() < 1e-6);
        assert!((inst.transform[13] - 4.0).abs() < 1e-6);
        assert!((inst.transform[14] - 5.0).abs() < 1e-6);
    }

    #[test]
    fn test_instance_buffer_byte_size() {
        let mut buf = InstanceBuffer::new();
        buf.push(InstanceData::new(identity_mat4(), [1.0; 4]));
        assert_eq!(buf.byte_size(), std::mem::size_of::<InstanceData>());
    }

    // ---- ParticleGpuSystem ----

    #[test]
    fn test_particle_gpu_system_add_and_count() {
        let mut sys = ParticleGpuSystem::new();
        sys.add([0.0, 0.0, 0.0], 0.1, [1.0, 1.0, 1.0, 1.0]);
        sys.add([1.0, 0.0, 0.0], 0.1, [1.0, 0.0, 0.0, 1.0]);
        assert_eq!(sys.particles.len(), 2);
    }

    #[test]
    fn test_particle_depth_sort_front_to_back_is_descending_distance() {
        let mut sys = ParticleGpuSystem::new();
        // Camera at origin; particles at various z
        sys.add([0.0, 0.0, 5.0], 0.1, [1.0; 4]); // far
        sys.add([0.0, 0.0, 1.0], 0.1, [1.0; 4]); // near
        sys.add([0.0, 0.0, 3.0], 0.1, [1.0; 4]); // mid

        let camera = [0.0, 0.0, 0.0];
        sys.sort_by_depth(camera);

        // After back-to-front sort, first particle should be the farthest
        let first_z = sys.particles[0].position[2];
        let last_z = sys.particles[sys.particles.len() - 1].position[2];
        assert!(
            first_z > last_z,
            "back-to-front: first particle should be farther (z={}) than last (z={})",
            first_z,
            last_z
        );
    }

    #[test]
    fn test_particle_build_instance_buffer() {
        let mut sys = ParticleGpuSystem::new();
        sys.add([1.0, 2.0, 3.0], 0.5, [1.0, 0.0, 0.0, 1.0]);
        let buf = sys.build_instance_buffer();
        assert_eq!(buf.len(), 1);
        // Scale should be 0.5
        assert!((buf.instances[0].transform[0] - 0.5).abs() < 1e-6);
    }

    // ---- TransferFunction ----

    #[test]
    fn test_transfer_function_opacity_clamped() {
        let mut tf = TransferFunction::new();
        tf.add(0.0, [1.0, 0.0, 0.0, 2.0]); // opacity > 1 → should clamp
        assert!(tf.entries[0].color[3] <= 1.0);
    }

    #[test]
    fn test_transfer_function_sample_below_range() {
        let mut tf = TransferFunction::new();
        tf.add(0.5, [1.0, 0.0, 0.0, 1.0]);
        let c = tf.sample(0.0);
        assert!((c[0] - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_transfer_function_sample_above_range() {
        let mut tf = TransferFunction::new();
        tf.add(0.5, [0.0, 1.0, 0.0, 1.0]);
        let c = tf.sample(1.0);
        assert!((c[1] - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_transfer_function_interpolates() {
        let mut tf = TransferFunction::new();
        tf.add(0.0, [0.0, 0.0, 0.0, 0.0]);
        tf.add(1.0, [1.0, 1.0, 1.0, 1.0]);
        let c = tf.sample(0.5);
        for (i, &ci) in c.iter().enumerate() {
            assert!(
                (ci - 0.5).abs() < 1e-5,
                "component {i} should be ~0.5, got {}",
                ci
            );
        }
    }

    #[test]
    fn test_transfer_function_opacity_in_range() {
        let mut tf = TransferFunction::new();
        tf.add(0.0, [1.0, 0.0, 0.0, 0.0]);
        tf.add(1.0, [0.0, 0.0, 1.0, 1.0]);
        for i in 0..=10 {
            let t = i as f32 / 10.0;
            let c = tf.sample(t);
            assert!(
                c[3] >= 0.0 && c[3] <= 1.0,
                "opacity out of range at t={t}: {}",
                c[3]
            );
        }
    }

    #[test]
    fn test_transfer_function_empty_returns_zero() {
        let tf = TransferFunction::new();
        let c = tf.sample(0.5);
        assert_eq!(c, [0.0; 4]);
    }

    // ---- VolumeRenderer ----

    #[test]
    fn test_volume_renderer_voxel_count() {
        let vr = VolumeRenderer::new([4, 5, 6], vec![0.0; 4 * 5 * 6]);
        assert_eq!(vr.voxel_count(), 120);
    }

    #[test]
    fn test_volume_renderer_sample_voxel() {
        let mut data = vec![0.0f32; 8];
        data[7] = 42.0; // last voxel (1,1,1) in 2x2x2
        let vr = VolumeRenderer::new([2, 2, 2], data);
        assert!((vr.sample_voxel(1, 1, 1) - 42.0).abs() < 1e-6);
        assert!((vr.sample_voxel(0, 0, 0)).abs() < 1e-6);
    }

    #[test]
    fn test_volume_renderer_out_of_bounds_returns_zero() {
        let vr = VolumeRenderer::new([2, 2, 2], vec![1.0; 8]);
        assert_eq!(vr.sample_voxel(5, 5, 5), 0.0);
    }

    // ---- PostProcessPass ----

    #[test]
    fn test_post_process_disabled_none_active() {
        let pp = PostProcessPass::disabled();
        assert!(!pp.any_enabled());
    }

    #[test]
    fn test_post_process_default_bloom_enabled() {
        let pp = PostProcessPass::default();
        assert!(pp.bloom.enabled);
        assert!(pp.any_enabled());
    }

    // ---- MeshBuilder ----

    #[test]
    fn test_mesh_builder_single_triangle_normals_unit_length() {
        let mut builder = MeshBuilder::new();
        let a = builder.add_vertex([0.0, 0.0, 0.0], [0.0, 0.0]);
        let b = builder.add_vertex([1.0, 0.0, 0.0], [1.0, 0.0]);
        let c = builder.add_vertex([0.0, 1.0, 0.0], [0.0, 1.0]);
        builder.add_triangle(a, b, c);
        let mesh = builder.build();

        for v in &mesh.vertices {
            let len = length3f(v.normal);
            assert!(
                (len - 1.0).abs() < 1e-5,
                "normal length should be 1.0, got {len:.6}"
            );
        }
    }

    #[test]
    fn test_mesh_builder_quad() {
        let mut builder = MeshBuilder::new();
        let a = builder.add_vertex([-1.0, 0.0, -1.0], [0.0, 0.0]);
        let b = builder.add_vertex([1.0, 0.0, -1.0], [1.0, 0.0]);
        let c = builder.add_vertex([1.0, 0.0, 1.0], [1.0, 1.0]);
        let d = builder.add_vertex([-1.0, 0.0, 1.0], [0.0, 1.0]);
        builder.add_triangle(a, b, c);
        builder.add_triangle(a, c, d);
        assert_eq!(builder.triangle_count(), 2);
        let mesh = builder.build();
        assert_eq!(mesh.vertices.len(), 4);
        assert_eq!(mesh.indices.len(), 6);
        // All normals should be aligned with the Y axis (up or down depending on winding)
        for v in &mesh.vertices {
            assert!(
                v.normal[1].abs() > 0.99,
                "y-normal magnitude should be ~1, got {:.6}",
                v.normal[1]
            );
            assert!(v.normal[0].abs() < 1e-5, "x-normal should be ~0");
            assert!(v.normal[2].abs() < 1e-5, "z-normal should be ~0");
        }
    }

    #[test]
    fn test_mesh_builder_vertex_count() {
        let mut b = MeshBuilder::new();
        for i in 0..10 {
            b.add_vertex([i as f32, 0.0, 0.0], [0.0, 0.0]);
        }
        assert_eq!(b.vertex_count(), 10);
    }

    // ---- LodSystem ----

    #[test]
    fn test_lod_system_selects_highest_detail_close() {
        let mut lod = LodSystem::new();
        lod.add_level(10.0, 0); // mesh 0 = high detail
        lod.add_level(50.0, 1); // mesh 1 = medium
        lod.add_level(200.0, 2); // mesh 2 = low

        assert_eq!(lod.select(5.0), Some(0));
    }

    #[test]
    fn test_lod_system_selects_lowest_detail_far() {
        let mut lod = LodSystem::new();
        lod.add_level(10.0, 0);
        lod.add_level(50.0, 1);
        lod.add_level(200.0, 2);

        assert_eq!(lod.select(500.0), Some(2));
    }

    #[test]
    fn test_lod_system_farther_gives_lower_lod() {
        let mut lod = LodSystem::new();
        lod.add_level(10.0, 0);
        lod.add_level(50.0, 1);
        lod.add_level(200.0, 2);

        let near = lod.select(5.0).unwrap();
        let far = lod.select(100.0).unwrap();
        assert!(
            near < far,
            "closer distance should give lower mesh_id (higher detail): near={near}, far={far}"
        );
    }

    #[test]
    fn test_lod_system_empty_returns_none() {
        let lod = LodSystem::new();
        assert_eq!(lod.select(0.0), None);
    }

    #[test]
    fn test_lod_system_level_count() {
        let mut lod = LodSystem::new();
        lod.add_level(10.0, 0);
        lod.add_level(100.0, 1);
        assert_eq!(lod.level_count(), 2);
    }

    // ---- FrustumCull ----

    #[test]
    fn test_frustum_cull_object_behind_camera_is_culled() {
        // Simple frustum: just a near plane facing +z
        // Near plane at z=1 with normal (0,0,1) → objects with z<1 are behind
        let near_plane = FrustumPlane::new([0.0, 0.0, 1.0], -1.0); // d = -1 → plane at z=1
        // Fill remaining planes with always-inside planes
        let always_in = FrustumPlane::new([0.0, 1.0, 0.0], 1e10);
        let planes = [
            near_plane, always_in, always_in, always_in, always_in, always_in,
        ];
        let frustum = FrustumCull::new(planes);

        // AABB entirely behind the near plane
        let behind_aabb = CullAabb {
            min: [-1.0, -1.0, -5.0],
            max: [1.0, 1.0, 0.5],
        };
        assert!(
            frustum.cull_aabb(&behind_aabb),
            "AABB behind camera should be culled"
        );
    }

    #[test]
    fn test_frustum_cull_object_in_front_not_culled() {
        let near_plane = FrustumPlane::new([0.0, 0.0, 1.0], -1.0);
        let always_in = FrustumPlane::new([0.0, 1.0, 0.0], 1e10);
        let planes = [
            near_plane, always_in, always_in, always_in, always_in, always_in,
        ];
        let frustum = FrustumCull::new(planes);

        let front_aabb = CullAabb {
            min: [-1.0, -1.0, 2.0],
            max: [1.0, 1.0, 5.0],
        };
        assert!(
            !frustum.cull_aabb(&front_aabb),
            "AABB in front of camera should not be culled"
        );
    }

    #[test]
    fn test_frustum_plane_distance_sign() {
        let plane = FrustumPlane::new([0.0, 1.0, 0.0], 0.0); // y=0 plane, inward normal = +y
        assert!(plane.distance([0.0, 1.0, 0.0]) > 0.0); // above plane
        assert!(plane.distance([0.0, -1.0, 0.0]) < 0.0); // below plane
    }

    #[test]
    fn test_cull_aabb_from_center_half() {
        let aabb = CullAabb::from_center_half([0.0, 0.0, 0.0], [1.0, 2.0, 3.0]);
        assert!((aabb.min[0] + 1.0).abs() < 1e-6);
        assert!((aabb.max[1] - 2.0).abs() < 1e-6);
        assert!((aabb.max[2] - 3.0).abs() < 1e-6);
    }

    // ---- UniformBlock ----

    #[test]
    fn test_uniform_block_default_identity_matrices() {
        let ub = UniformBlock::default();
        // Diagonal of column-major identity
        assert!((ub.view[0] - 1.0).abs() < 1e-6);
        assert!((ub.view[5] - 1.0).abs() < 1e-6);
        assert!((ub.view[10] - 1.0).abs() < 1e-6);
        assert!((ub.view[15] - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_uniform_block_ambient_positive() {
        let ub = UniformBlock::default();
        assert!(ub.ambient > 0.0);
    }

    // ---- SdfRenderer ----

    #[test]
    fn test_sdf_renderer_default_values() {
        let sdf = SdfRenderer::new();
        assert!(sdf.epsilon > 0.0);
        assert!(sdf.max_dist > 0.0);
        assert!(sdf.max_steps > 0);
    }
}
