// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! GPU instanced rendering data structures.
//!
//! Provides CPU-side data structures for instance-buffer management, mesh
//! templates, billboard orientation, particle instances, and frustum culling.
//! All types use `f64` for precision; a downstream GPU back-end would convert
//! to `f32` before uploading.

// ─────────────────────────────────────────────────────────────────────────────
// InstanceData
// ─────────────────────────────────────────────────────────────────────────────

/// Per-instance data passed to the GPU.
#[derive(Debug, Clone, PartialEq)]
pub struct InstanceData {
    /// Column-major 4×4 world transform matrix.
    pub transform: [f64; 16],
    /// RGBA color, each component in \[0, 1\].
    pub color: [f64; 4],
    /// Non-uniform scale applied on top of `transform`.
    pub scale: [f64; 3],
    /// Four general-purpose floats available in the shader.
    pub custom_data: [f64; 4],
}

impl InstanceData {
    /// Create a new `InstanceData` with the given fields.
    pub fn new(
        transform: [f64; 16],
        color: [f64; 4],
        scale: [f64; 3],
        custom_data: [f64; 4],
    ) -> Self {
        Self {
            transform,
            color,
            scale,
            custom_data,
        }
    }

    /// Return `InstanceData` with an identity transform, white colour, unit
    /// scale, and zeroed custom data.
    pub fn identity() -> Self {
        Self {
            transform: identity_mat4(),
            color: [1.0, 1.0, 1.0, 1.0],
            scale: [1.0, 1.0, 1.0],
            custom_data: [0.0; 4],
        }
    }
}

impl Default for InstanceData {
    fn default() -> Self {
        Self::identity()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// InstanceBuffer
// ─────────────────────────────────────────────────────────────────────────────

/// A CPU-side buffer that holds a list of [`InstanceData`] entries.
#[derive(Debug, Clone)]
pub struct InstanceBuffer {
    /// Stored instance data.
    pub instances: Vec<InstanceData>,
    /// Set to `true` whenever the buffer is mutated and needs re-uploading.
    pub upload_dirty: bool,
}

impl InstanceBuffer {
    /// Create an empty buffer.
    pub fn new() -> Self {
        Self {
            instances: Vec::new(),
            upload_dirty: false,
        }
    }

    /// Create a buffer with the given initial capacity (no instances added).
    pub fn with_capacity(cap: usize) -> Self {
        Self {
            instances: Vec::with_capacity(cap),
            upload_dirty: false,
        }
    }

    /// Append an instance and mark the buffer dirty.
    pub fn add(&mut self, instance: InstanceData) {
        self.instances.push(instance);
        self.upload_dirty = true;
    }

    /// Remove the instance at `idx` (swap-remove for O(1)) and mark dirty.
    ///
    /// Panics if `idx` is out of bounds.
    pub fn remove(&mut self, idx: usize) {
        self.instances.swap_remove(idx);
        self.upload_dirty = true;
    }

    /// Replace the instance at `idx` and mark dirty.
    ///
    /// Panics if `idx` is out of bounds.
    pub fn update(&mut self, idx: usize, instance: InstanceData) {
        self.instances[idx] = instance;
        self.upload_dirty = true;
    }

    /// Mark the buffer as uploaded (clear the dirty flag).
    pub fn mark_clean(&mut self) {
        self.upload_dirty = false;
    }

    /// Number of instances in the buffer.
    pub fn len(&self) -> usize {
        self.instances.len()
    }

    /// Return `true` if the buffer has no instances.
    pub fn is_empty(&self) -> bool {
        self.instances.is_empty()
    }

    /// Remove all instances and mark dirty.
    pub fn clear(&mut self) {
        self.instances.clear();
        self.upload_dirty = true;
    }
}

impl Default for InstanceBuffer {
    fn default() -> Self {
        Self::new()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// MeshTemplate
// ─────────────────────────────────────────────────────────────────────────────

/// A reusable mesh definition used as the template for instanced rendering.
#[derive(Debug, Clone)]
pub struct MeshTemplate {
    /// Vertex positions.
    pub vertex_data: Vec<[f64; 3]>,
    /// Per-vertex normals (same length as `vertex_data`).
    pub normal_data: Vec<[f64; 3]>,
    /// Triangle index list; length is a multiple of 3.
    pub index_data: Vec<usize>,
    /// Bounding sphere radius in local space.
    pub bounding_sphere_radius: f64,
}

impl MeshTemplate {
    /// Create a new mesh template.
    pub fn new(
        vertex_data: Vec<[f64; 3]>,
        normal_data: Vec<[f64; 3]>,
        index_data: Vec<usize>,
        bounding_sphere_radius: f64,
    ) -> Self {
        Self {
            vertex_data,
            normal_data,
            index_data,
            bounding_sphere_radius,
        }
    }

    /// Number of triangles in the mesh.
    pub fn triangle_count(&self) -> usize {
        self.index_data.len() / 3
    }

    /// Number of vertices.
    pub fn vertex_count(&self) -> usize {
        self.vertex_data.len()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// InstancedRenderer
// ─────────────────────────────────────────────────────────────────────────────

/// An instanced renderer that pairs a [`MeshTemplate`] with an
/// [`InstanceBuffer`] and provides frustum culling.
#[derive(Debug, Clone)]
pub struct InstancedRenderer {
    /// The shared mesh geometry.
    pub template: MeshTemplate,
    /// The per-instance data buffer.
    pub buffer: InstanceBuffer,
}

impl InstancedRenderer {
    /// Create a new renderer with the given template.
    pub fn new(template: MeshTemplate) -> Self {
        Self {
            template,
            buffer: InstanceBuffer::new(),
        }
    }

    /// Perform frustum culling against `camera_planes` (6 planes in
    /// `[a, b, c, d]` format: `a·x + b·y + c·z + d ≥ 0` means inside).
    ///
    /// Returns the number of instances that passed culling (are visible).
    /// Instances whose bounding sphere intersects or is inside all six planes
    /// are considered visible.  The instance transforms are interpreted as
    /// translation-only (the last column / position vector).
    pub fn cull_frustum(&self, camera_planes: &[[f64; 4]; 6]) -> usize {
        let radius = self.template.bounding_sphere_radius;
        self.buffer
            .instances
            .iter()
            .filter(|inst| {
                // Extract translation from column-major 4×4 matrix
                let center = [inst.transform[12], inst.transform[13], inst.transform[14]];
                // Apply non-uniform scale: use max component as conservative radius
                let max_scale = inst.scale[0].max(inst.scale[1]).max(inst.scale[2]);
                let effective_radius = radius * max_scale;
                frustum_cull(center, effective_radius, camera_planes)
            })
            .count()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ParticleInstance
// ─────────────────────────────────────────────────────────────────────────────

/// A single particle with physics state and lifetime tracking.
#[derive(Debug, Clone)]
pub struct ParticleInstance {
    /// World-space position.
    pub position: [f64; 3],
    /// World-space velocity.
    pub velocity: [f64; 3],
    /// Current lifetime in seconds.
    pub lifetime: f64,
    /// Maximum lifetime in seconds.
    pub max_lifetime: f64,
}

impl ParticleInstance {
    /// Create a new particle.
    pub fn new(position: [f64; 3], velocity: [f64; 3], max_lifetime: f64) -> Self {
        Self {
            position,
            velocity,
            lifetime: 0.0,
            max_lifetime,
        }
    }

    /// Return `true` when the particle has exceeded its maximum lifetime.
    pub fn is_dead(&self) -> bool {
        self.lifetime >= self.max_lifetime
    }

    /// Normalised age in \[0, 1\].
    pub fn age_fraction(&self) -> f64 {
        if self.max_lifetime <= 0.0 {
            1.0
        } else {
            (self.lifetime / self.max_lifetime).clamp(0.0, 1.0)
        }
    }

    /// Advance the particle by `dt` seconds (Euler integration).
    pub fn step(&mut self, dt: f64) {
        for k in 0..3 {
            self.position[k] += self.velocity[k] * dt;
        }
        self.lifetime += dt;
    }

    /// Convert to [`InstanceData`] for GPU upload.
    ///
    /// Uses a translation-only transform, fades alpha with remaining lifetime,
    /// and stores `age_fraction` in `custom_data[0]`.
    pub fn to_instance_data(&self) -> InstanceData {
        let mut transform = identity_mat4();
        // Set translation column
        transform[12] = self.position[0];
        transform[13] = self.position[1];
        transform[14] = self.position[2];

        let alpha = 1.0 - self.age_fraction();
        InstanceData {
            transform,
            color: [1.0, 1.0, 1.0, alpha],
            scale: [1.0, 1.0, 1.0],
            custom_data: [self.age_fraction(), 0.0, 0.0, 0.0],
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// BillboardInstance
// ─────────────────────────────────────────────────────────────────────────────

/// A billboard that always faces the camera.
#[derive(Debug, Clone)]
pub struct BillboardInstance {
    /// World-space position of the billboard centre.
    pub position: [f64; 3],
    /// Uniform size (half-extent) of the billboard quad.
    pub size: f64,
}

impl BillboardInstance {
    /// Create a new billboard.
    pub fn new(position: [f64; 3], size: f64) -> Self {
        Self { position, size }
    }

    /// Compute a column-major 4×4 transform matrix that makes this billboard
    /// face `camera_pos` with the Y-axis as the "up" direction.
    ///
    /// The resulting matrix is a scale-then-look-at transform.
    pub fn always_face_camera(&self, camera_pos: [f64; 3]) -> [f64; 16] {
        // Forward vector: billboard → camera
        let fwd = vec3_sub(camera_pos, self.position);
        let fwd_len = vec3_len(fwd);
        let forward = if fwd_len < 1e-12 {
            [0.0, 0.0, 1.0]
        } else {
            vec3_scale(fwd, 1.0 / fwd_len)
        };

        // Right = forward × up  (up = Y)
        let world_up = [0.0, 1.0, 0.0];
        let right_raw = vec3_cross(forward, world_up);
        let right_len = vec3_len(right_raw);
        let right = if right_len < 1e-12 {
            [1.0, 0.0, 0.0]
        } else {
            vec3_scale(right_raw, 1.0 / right_len)
        };

        // Up = right × forward  (recompute for orthogonality)
        let up = vec3_cross(right, forward);

        // Build column-major 4×4 with scale incorporated
        let s = self.size;
        #[rustfmt::skip]
        let mat = [
            right[0] * s, right[1] * s, right[2] * s, 0.0,
            up[0]    * s, up[1]    * s, up[2]    * s, 0.0,
            forward[0],   forward[1],   forward[2],   0.0,
            self.position[0], self.position[1], self.position[2], 1.0,
        ];
        mat
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Free functions
// ─────────────────────────────────────────────────────────────────────────────

/// Build a UV sphere [`MeshTemplate`] with the given radius and resolution.
///
/// `lat` is the number of latitude bands (rings); `lon` is the number of
/// longitude segments.  Both must be at least 2.
pub fn build_sphere_mesh(radius: f64, lat: usize, lon: usize) -> MeshTemplate {
    use std::f64::consts::PI;

    let lat = lat.max(2);
    let lon = lon.max(3);

    let mut vertices: Vec<[f64; 3]> = Vec::new();
    let mut normals: Vec<[f64; 3]> = Vec::new();
    let mut indices: Vec<usize> = Vec::new();

    for i in 0..=lat {
        let theta = PI * i as f64 / lat as f64;
        let sin_t = theta.sin();
        let cos_t = theta.cos();

        for j in 0..=lon {
            let phi = 2.0 * PI * j as f64 / lon as f64;
            let x = sin_t * phi.cos();
            let y = cos_t;
            let z = sin_t * phi.sin();
            vertices.push([x * radius, y * radius, z * radius]);
            normals.push([x, y, z]);
        }
    }

    for i in 0..lat {
        for j in 0..lon {
            let row = lon + 1;
            let a = i * row + j;
            let b = a + row;
            let c = b + 1;
            let d = a + 1;
            indices.push(a);
            indices.push(b);
            indices.push(d);
            indices.push(b);
            indices.push(c);
            indices.push(d);
        }
    }

    MeshTemplate::new(vertices, normals, indices, radius)
}

/// Test whether a sphere (world-space `center`, `radius`) is inside or
/// intersecting all six frustum half-spaces.
///
/// Each plane is `[a, b, c, d]` where `a·x + b·y + c·z + d ≥ 0` means inside.
/// Returns `true` if the sphere is not entirely outside any plane.
pub fn frustum_cull(center: [f64; 3], radius: f64, planes: &[[f64; 4]; 6]) -> bool {
    for plane in planes {
        let dist = plane[0] * center[0] + plane[1] * center[1] + plane[2] * center[2] + plane[3];
        if dist < -radius {
            return false; // entirely outside this plane
        }
    }
    true
}

// ─────────────────────────────────────────────────────────────────────────────
// Internal math helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Column-major 4×4 identity matrix.
fn identity_mat4() -> [f64; 16] {
    [
        1.0, 0.0, 0.0, 0.0, // col 0
        0.0, 1.0, 0.0, 0.0, // col 1
        0.0, 0.0, 1.0, 0.0, // col 2
        0.0, 0.0, 0.0, 1.0, // col 3
    ]
}

/// Subtract two 3-vectors.
fn vec3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Euclidean length of a 3-vector.
fn vec3_len(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

/// Scale a 3-vector by a scalar.
fn vec3_scale(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}

/// Cross product of two 3-vectors.
fn vec3_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Dot product of two 3-vectors.
#[cfg(test)]
#[inline]
fn vec3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ---- InstanceData ----

    #[test]
    fn test_instance_data_identity_default() {
        let inst = InstanceData::identity();
        // Diagonal elements of the identity matrix
        assert!((inst.transform[0] - 1.0).abs() < 1e-12);
        assert!((inst.transform[5] - 1.0).abs() < 1e-12);
        assert!((inst.transform[10] - 1.0).abs() < 1e-12);
        assert!((inst.transform[15] - 1.0).abs() < 1e-12);
        // Off-diagonal
        assert!(inst.transform[1].abs() < 1e-12);
    }

    #[test]
    fn test_instance_data_default_same_as_identity() {
        let a = InstanceData::default();
        let b = InstanceData::identity();
        assert_eq!(a, b);
    }

    #[test]
    fn test_instance_data_fields() {
        let inst = InstanceData::new(
            identity_mat4(),
            [1.0, 0.0, 0.0, 1.0],
            [2.0, 2.0, 2.0],
            [7.0, 8.0, 9.0, 10.0],
        );
        assert_eq!(inst.color, [1.0, 0.0, 0.0, 1.0]);
        assert_eq!(inst.scale, [2.0, 2.0, 2.0]);
        assert_eq!(inst.custom_data[0], 7.0);
    }

    // ---- InstanceBuffer ----

    #[test]
    fn test_buffer_starts_empty() {
        let buf = InstanceBuffer::new();
        assert!(buf.is_empty());
        assert_eq!(buf.len(), 0);
        assert!(!buf.upload_dirty);
    }

    #[test]
    fn test_buffer_add_sets_dirty() {
        let mut buf = InstanceBuffer::new();
        buf.add(InstanceData::identity());
        assert_eq!(buf.len(), 1);
        assert!(buf.upload_dirty);
    }

    #[test]
    fn test_buffer_mark_clean_clears_dirty() {
        let mut buf = InstanceBuffer::new();
        buf.add(InstanceData::identity());
        buf.mark_clean();
        assert!(!buf.upload_dirty);
    }

    #[test]
    fn test_buffer_remove_sets_dirty() {
        let mut buf = InstanceBuffer::new();
        buf.add(InstanceData::identity());
        buf.mark_clean();
        buf.remove(0);
        assert!(buf.is_empty());
        assert!(buf.upload_dirty);
    }

    #[test]
    fn test_buffer_update_sets_dirty() {
        let mut buf = InstanceBuffer::new();
        buf.add(InstanceData::identity());
        buf.mark_clean();
        let mut inst = InstanceData::identity();
        inst.color = [0.5, 0.5, 0.5, 1.0];
        buf.update(0, inst.clone());
        assert!(buf.upload_dirty);
        assert_eq!(buf.instances[0].color, [0.5, 0.5, 0.5, 1.0]);
    }

    #[test]
    fn test_buffer_clear_removes_all() {
        let mut buf = InstanceBuffer::new();
        for _ in 0..5 {
            buf.add(InstanceData::identity());
        }
        buf.clear();
        assert!(buf.is_empty());
        assert!(buf.upload_dirty);
    }

    #[test]
    fn test_buffer_with_capacity() {
        let buf = InstanceBuffer::with_capacity(64);
        assert!(buf.is_empty());
        assert!(buf.instances.capacity() >= 64);
    }

    // ---- MeshTemplate ----

    #[test]
    fn test_mesh_template_triangle_count() {
        let tmpl = MeshTemplate::new(
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[0.0, 0.0, 1.0]; 3],
            vec![0, 1, 2],
            1.0,
        );
        assert_eq!(tmpl.triangle_count(), 1);
        assert_eq!(tmpl.vertex_count(), 3);
    }

    #[test]
    fn test_mesh_template_bounding_radius() {
        let tmpl = build_sphere_mesh(2.5, 8, 8);
        assert!((tmpl.bounding_sphere_radius - 2.5).abs() < 1e-12);
    }

    // ---- build_sphere_mesh ----

    #[test]
    fn test_build_sphere_mesh_vertices_on_sphere() {
        let r = 1.0;
        let tmpl = build_sphere_mesh(r, 8, 8);
        for v in &tmpl.vertex_data {
            let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
            assert!((len - r).abs() < 1e-10, "vertex not on sphere: len={len}");
        }
    }

    #[test]
    fn test_build_sphere_mesh_normals_unit() {
        let tmpl = build_sphere_mesh(3.0, 4, 4);
        for n in &tmpl.normal_data {
            let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
            assert!((len - 1.0).abs() < 1e-10, "normal not unit: len={len}");
        }
    }

    #[test]
    fn test_build_sphere_mesh_index_count_multiple_of_3() {
        let tmpl = build_sphere_mesh(1.0, 6, 6);
        assert_eq!(tmpl.index_data.len() % 3, 0);
    }

    #[test]
    fn test_build_sphere_mesh_higher_resolution_more_vertices() {
        let low = build_sphere_mesh(1.0, 4, 4);
        let high = build_sphere_mesh(1.0, 16, 16);
        assert!(high.vertex_count() > low.vertex_count());
        assert!(high.triangle_count() > low.triangle_count());
    }

    #[test]
    fn test_build_sphere_mesh_indices_in_bounds() {
        let tmpl = build_sphere_mesh(1.0, 4, 4);
        let n = tmpl.vertex_count();
        for &idx in &tmpl.index_data {
            assert!(idx < n, "index {idx} out of bounds (n={n})");
        }
    }

    // ---- frustum_cull ----

    /// Create 6 planes that form a unit cube \[-1,1\]³ frustum.
    fn unit_cube_planes() -> [[f64; 4]; 6] {
        [
            [1.0, 0.0, 0.0, 1.0],  // x >= -1  → a·x+d ≥ 0
            [-1.0, 0.0, 0.0, 1.0], // x <=  1
            [0.0, 1.0, 0.0, 1.0],  // y >= -1
            [0.0, -1.0, 0.0, 1.0], // y <=  1
            [0.0, 0.0, 1.0, 1.0],  // z >= -1
            [0.0, 0.0, -1.0, 1.0], // z <=  1
        ]
    }

    #[test]
    fn test_frustum_cull_inside() {
        let planes = unit_cube_planes();
        assert!(frustum_cull([0.0, 0.0, 0.0], 0.1, &planes));
    }

    #[test]
    fn test_frustum_cull_outside_x() {
        let planes = unit_cube_planes();
        // Center at x=3, radius=0.1 → entirely outside x<=1 plane
        assert!(!frustum_cull([3.0, 0.0, 0.0], 0.1, &planes));
    }

    #[test]
    fn test_frustum_cull_straddling_plane() {
        let planes = unit_cube_planes();
        // Center at x=0.95, radius=0.5 → straddles x<=1, but still partially inside
        assert!(frustum_cull([0.95, 0.0, 0.0], 0.5, &planes));
    }

    #[test]
    fn test_frustum_cull_large_sphere_encompasses() {
        let planes = unit_cube_planes();
        // Very large sphere centred at origin — passes all planes
        assert!(frustum_cull([0.0, 0.0, 0.0], 100.0, &planes));
    }

    // ---- InstancedRenderer cull_frustum ----

    #[test]
    fn test_instanced_renderer_cull_all_visible() {
        let tmpl = build_sphere_mesh(0.1, 4, 4);
        let mut renderer = InstancedRenderer::new(tmpl);
        for _ in 0..5 {
            renderer.buffer.add(InstanceData::identity()); // all at origin
        }
        let planes = unit_cube_planes();
        let visible = renderer.cull_frustum(&planes);
        assert_eq!(visible, 5);
    }

    #[test]
    fn test_instanced_renderer_cull_some_outside() {
        let tmpl = build_sphere_mesh(0.1, 4, 4);
        let mut renderer = InstancedRenderer::new(tmpl);

        // Instance at origin (inside)
        renderer.buffer.add(InstanceData::identity());

        // Instance at x=5 (outside)
        let mut far_inst = InstanceData::identity();
        far_inst.transform[12] = 5.0;
        renderer.buffer.add(far_inst);

        let planes = unit_cube_planes();
        let visible = renderer.cull_frustum(&planes);
        assert_eq!(visible, 1);
    }

    // ---- ParticleInstance ----

    #[test]
    fn test_particle_age_fraction_initial() {
        let p = ParticleInstance::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 2.0);
        assert!(p.age_fraction().abs() < 1e-12);
        assert!(!p.is_dead());
    }

    #[test]
    fn test_particle_step_advances_position() {
        let mut p = ParticleInstance::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 10.0);
        p.step(0.5);
        assert!((p.position[0] - 0.5).abs() < 1e-12);
        assert!((p.lifetime - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_particle_is_dead_after_max_lifetime() {
        let mut p = ParticleInstance::new([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 1.0);
        p.step(1.0);
        assert!(p.is_dead());
    }

    #[test]
    fn test_particle_age_fraction_clamps_to_one() {
        let mut p = ParticleInstance::new([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 1.0);
        p.step(5.0); // overshoot
        let af = p.age_fraction();
        assert!((af - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_particle_to_instance_data_translation() {
        let p = ParticleInstance::new([3.0, 4.0, 5.0], [0.0, 0.0, 0.0], 1.0);
        let inst = p.to_instance_data();
        assert!((inst.transform[12] - 3.0).abs() < 1e-12);
        assert!((inst.transform[13] - 4.0).abs() < 1e-12);
        assert!((inst.transform[14] - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_particle_to_instance_data_alpha_fades() {
        let mut p = ParticleInstance::new([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 2.0);
        let inst0 = p.to_instance_data();
        p.step(1.0); // halfway
        let inst1 = p.to_instance_data();
        assert!(
            inst1.color[3] < inst0.color[3],
            "alpha should decrease as particle ages"
        );
    }

    #[test]
    fn test_particle_to_instance_data_age_fraction_in_custom() {
        let mut p = ParticleInstance::new([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 4.0);
        p.step(1.0); // 25% through
        let inst = p.to_instance_data();
        assert!((inst.custom_data[0] - 0.25).abs() < 1e-10);
    }

    // ---- BillboardInstance ----

    #[test]
    fn test_billboard_faces_camera_forward_not_zero() {
        let bb = BillboardInstance::new([0.0, 0.0, 0.0], 1.0);
        let mat = bb.always_face_camera([0.0, 0.0, 5.0]);
        // Forward vector (column 2) should point toward camera (positive z)
        let fz = mat[8]; // column-major: col 2, row 0
        let fy = mat[9];
        let fx = mat[10];
        let len = (fx * fx + fy * fy + fz * fz).sqrt();
        assert!(len > 0.1, "forward vector should be non-zero");
    }

    #[test]
    fn test_billboard_translation_preserved() {
        let pos = [1.0, 2.0, 3.0];
        let bb = BillboardInstance::new(pos, 1.0);
        let mat = bb.always_face_camera([0.0, 0.0, 10.0]);
        // Translation is in the 4th column (indices 12, 13, 14)
        assert!((mat[12] - 1.0).abs() < 1e-12);
        assert!((mat[13] - 2.0).abs() < 1e-12);
        assert!((mat[14] - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_billboard_size_scales_axes() {
        let bb1 = BillboardInstance::new([0.0, 0.0, 0.0], 1.0);
        let bb2 = BillboardInstance::new([0.0, 0.0, 0.0], 2.0);
        let m1 = bb1.always_face_camera([5.0, 0.0, 0.0]);
        let m2 = bb2.always_face_camera([5.0, 0.0, 0.0]);
        // Column 0 magnitude of m2 should be ~2× that of m1
        let len1 = (m1[0] * m1[0] + m1[1] * m1[1] + m1[2] * m1[2]).sqrt();
        let len2 = (m2[0] * m2[0] + m2[1] * m2[1] + m2[2] * m2[2]).sqrt();
        assert!(
            (len2 - 2.0 * len1).abs() < 1e-10,
            "size=2 should double axis magnitudes"
        );
    }

    #[test]
    fn test_billboard_right_and_up_orthogonal() {
        let bb = BillboardInstance::new([0.0, 0.0, 0.0], 1.0);
        let mat = bb.always_face_camera([3.0, 4.0, 5.0]);
        // Right = col 0 / size, Up = col 1 / size — must be orthogonal
        let right = [mat[0], mat[1], mat[2]];
        let up = [mat[4], mat[5], mat[6]];
        let dot = vec3_dot(right, up);
        assert!(
            dot.abs() < 1e-10,
            "right and up should be orthogonal, dot={dot}"
        );
    }

    // ---- Internal math helpers ----

    #[test]
    fn test_identity_mat4_diagonal() {
        let m = identity_mat4();
        assert!((m[0] - 1.0).abs() < 1e-12);
        assert!((m[5] - 1.0).abs() < 1e-12);
        assert!((m[10] - 1.0).abs() < 1e-12);
        assert!((m[15] - 1.0).abs() < 1e-12);
        assert!(m[1].abs() < 1e-12);
    }

    #[test]
    fn test_vec3_cross_unit_vectors() {
        let x = [1.0, 0.0, 0.0];
        let y = [0.0, 1.0, 0.0];
        let z = vec3_cross(x, y);
        assert!((z[0]).abs() < 1e-12);
        assert!((z[1]).abs() < 1e-12);
        assert!((z[2] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_vec3_dot_orthogonal() {
        let x = [1.0, 0.0, 0.0];
        let y = [0.0, 1.0, 0.0];
        assert!(vec3_dot(x, y).abs() < 1e-12);
    }

    #[test]
    fn test_vec3_len() {
        let v = [3.0, 4.0, 0.0];
        assert!((vec3_len(v) - 5.0).abs() < 1e-12);
    }
}
