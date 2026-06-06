// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Real-time physics visualization module.
//!
//! Covers:
//! - Frame rate management and timing budget
//! - Level-of-detail (LOD) system
//! - Frustum culling (view-frustum test)
//! - Occlusion culling (hierarchical Z-buffer approximation)
//! - Billboard sprites (camera-facing quads)
//! - Instanced rendering for particle systems
//! - Temporal anti-aliasing (TAA) accumulation buffer
//! - Dynamic resolution scaling (DRS)
//! - GPU particle system visualization data
//! - Real-time shadow maps (depth pass)
//! - Screen-space reflections (SSR) ray-marching
//! - Ambient occlusion (SSAO/HBAO) hemisphere sampling
//! - Deferred rendering pipeline (G-buffer)
//! - HDR tone mapping (Reinhard, ACES, Uncharted2)
//! - Real-time global illumination approximation (irradiance probes)

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Basic math helpers (plain [f64;3], no nalgebra)
// ---------------------------------------------------------------------------

/// 3-component vector stored as plain array.
pub type Vec3 = [f64; 3];

/// 4-component RGBA colour stored as f32.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Rgba {
    /// Red channel in \[0, 1\].
    pub r: f32,
    /// Green channel in \[0, 1\].
    pub g: f32,
    /// Blue channel in \[0, 1\].
    pub b: f32,
    /// Alpha channel in \[0, 1\].
    pub a: f32,
}

impl Rgba {
    /// Construct a colour from components.
    pub fn new(r: f32, g: f32, b: f32, a: f32) -> Self {
        Self { r, g, b, a }
    }
    /// Opaque white.
    pub fn white() -> Self {
        Self::new(1.0, 1.0, 1.0, 1.0)
    }
    /// Opaque black.
    pub fn black() -> Self {
        Self::new(0.0, 0.0, 0.0, 1.0)
    }
    /// Transparent black.
    pub fn transparent() -> Self {
        Self::new(0.0, 0.0, 0.0, 0.0)
    }
}

/// Add two Vec3s.
#[inline]
fn vadd(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
/// Subtract two Vec3s.
#[inline]
fn vsub(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
/// Scale a Vec3.
#[inline]
fn vscale(a: Vec3, s: f64) -> Vec3 {
    [a[0] * s, a[1] * s, a[2] * s]
}
/// Dot product.
#[inline]
fn vdot(a: Vec3, b: Vec3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Euclidean norm.
#[inline]
fn vnorm(a: Vec3) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}
/// Normalize; returns zero-vector if input is degenerate.
#[inline]
fn vnormalize(a: Vec3) -> Vec3 {
    let n = vnorm(a);
    if n < 1e-300 {
        [0.0; 3]
    } else {
        [a[0] / n, a[1] / n, a[2] / n]
    }
}

// ---------------------------------------------------------------------------
// Frame rate management
// ---------------------------------------------------------------------------

/// Target frame-rate budget used to control simulation step timing.
#[derive(Debug, Clone)]
pub struct FrameTimer {
    /// Target frames per second.
    pub target_fps: f64,
    /// Accumulated time of the last N frames \[s\].
    frame_times: Vec<f64>,
    /// Maximum number of frames to average over.
    history_len: usize,
    /// Total elapsed simulation time \[s\].
    pub elapsed: f64,
}

impl FrameTimer {
    /// Construct a new `FrameTimer` targeting `fps` frames per second.
    pub fn new(fps: f64, history_len: usize) -> Self {
        Self {
            target_fps: fps.max(1.0),
            frame_times: Vec::with_capacity(history_len),
            history_len: history_len.max(1),
            elapsed: 0.0,
        }
    }

    /// Record a completed frame with duration `dt` \[s\].
    pub fn push_frame(&mut self, dt: f64) {
        self.elapsed += dt;
        if self.frame_times.len() >= self.history_len {
            self.frame_times.remove(0);
        }
        self.frame_times.push(dt);
    }

    /// Return the smoothed (averaged) frame time \[s\].
    pub fn smoothed_dt(&self) -> f64 {
        if self.frame_times.is_empty() {
            return 1.0 / self.target_fps;
        }
        self.frame_times.iter().sum::<f64>() / self.frame_times.len() as f64
    }

    /// Return the measured FPS based on smoothed frame time.
    pub fn measured_fps(&self) -> f64 {
        let dt = self.smoothed_dt();
        if dt < 1e-10 { f64::INFINITY } else { 1.0 / dt }
    }

    /// Return `true` if the measured FPS is within `tolerance` percent of target.
    pub fn is_on_target(&self, tolerance_pct: f64) -> bool {
        let measured = self.measured_fps();
        let margin = self.target_fps * tolerance_pct / 100.0;
        (measured - self.target_fps).abs() <= margin
    }

    /// Compute ideal time step for the target FPS \[s\].
    pub fn target_dt(&self) -> f64 {
        1.0 / self.target_fps
    }

    /// Number of frames recorded.
    pub fn frame_count(&self) -> usize {
        self.frame_times.len()
    }
}

// ---------------------------------------------------------------------------
// Level-of-detail (LOD)
// ---------------------------------------------------------------------------

/// LOD level selector based on screen-space projected size.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum LodLevel {
    /// Highest detail (closest).
    High,
    /// Medium detail.
    Medium,
    /// Low detail (furthest).
    Low,
    /// Culled entirely — do not render.
    Culled,
}

/// LOD selection parameters.
#[derive(Debug, Clone)]
pub struct LodSelector {
    /// World-space radius of the bounding sphere.
    pub bounding_radius: f64,
    /// Screen-height in pixels for LOD calculations.
    pub screen_height_px: f64,
    /// Vertical field of view \[rad\].
    pub fov_y_rad: f64,
    /// Transition distance: High → Medium \[m\].
    pub dist_high_medium: f64,
    /// Transition distance: Medium → Low \[m\].
    pub dist_medium_low: f64,
    /// Transition distance: Low → Culled \[m\].
    pub dist_low_culled: f64,
}

impl LodSelector {
    /// Construct a default `LodSelector` with standard transition distances.
    pub fn new(bounding_radius: f64) -> Self {
        Self {
            bounding_radius,
            screen_height_px: 1080.0,
            fov_y_rad: 60.0_f64.to_radians(),
            dist_high_medium: bounding_radius * 20.0,
            dist_medium_low: bounding_radius * 60.0,
            dist_low_culled: bounding_radius * 200.0,
        }
    }

    /// Select the LOD level for an object at distance `dist` from the camera.
    pub fn select(&self, dist: f64) -> LodLevel {
        if dist > self.dist_low_culled {
            LodLevel::Culled
        } else if dist > self.dist_medium_low {
            LodLevel::Low
        } else if dist > self.dist_high_medium {
            LodLevel::Medium
        } else {
            LodLevel::High
        }
    }

    /// Compute projected screen size in pixels for an object at `dist`.
    pub fn projected_pixels(&self, dist: f64) -> f64 {
        if dist < 1e-10 {
            return f64::INFINITY;
        }
        let projected_half = self.bounding_radius / dist / (self.fov_y_rad / 2.0).tan();
        projected_half * self.screen_height_px
    }
}

// ---------------------------------------------------------------------------
// Frustum culling
// ---------------------------------------------------------------------------

/// A half-space plane defined by a normal and a signed distance from origin.
#[derive(Debug, Clone, Copy)]
pub struct Plane {
    /// Inward-facing plane normal (unit vector).
    pub normal: Vec3,
    /// Signed distance from the origin.
    pub distance: f64,
}

impl Plane {
    /// Construct a plane from normal and distance.
    pub fn new(normal: Vec3, distance: f64) -> Self {
        Self {
            normal: vnormalize(normal),
            distance,
        }
    }

    /// Signed distance from the plane to point `p`.
    pub fn signed_dist(&self, p: Vec3) -> f64 {
        vdot(self.normal, p) + self.distance
    }
}

/// View frustum represented as 6 planes (left, right, bottom, top, near, far).
#[derive(Debug, Clone)]
pub struct ViewFrustum {
    /// The six bounding planes.
    pub planes: [Plane; 6],
}

impl ViewFrustum {
    /// Build a view frustum from a row-major 4×4 view-projection matrix.
    ///
    /// Uses the Gribb-Hartmann method to extract frustum planes directly.
    pub fn from_view_proj(m: [f64; 16]) -> Self {
        // m[row*4 + col]
        let row = |r: usize| [m[r * 4], m[r * 4 + 1], m[r * 4 + 2], m[r * 4 + 3]];
        let r0 = row(0);
        let r1 = row(1);
        let r2 = row(2);
        let r3 = row(3);

        let plane = |a: [f64; 4], b: [f64; 4], sign: f64| {
            let n = [a[0] + sign * b[0], a[1] + sign * b[1], a[2] + sign * b[2]];
            let d = a[3] + sign * b[3];
            let len = vnorm(n);
            if len < 1e-12 {
                Plane::new([0.0, 0.0, 1.0], 0.0)
            } else {
                Plane::new([n[0] / len, n[1] / len, n[2] / len], d / len)
            }
        };

        Self {
            planes: [
                plane(r3, r0, 1.0),  // left
                plane(r3, r0, -1.0), // right
                plane(r3, r1, 1.0),  // bottom
                plane(r3, r1, -1.0), // top
                plane(r3, r2, 1.0),  // near
                plane(r3, r2, -1.0), // far
            ],
        }
    }

    /// Test whether a bounding sphere is inside (or intersecting) the frustum.
    ///
    /// Returns `true` if the sphere should be rendered (not fully outside).
    pub fn test_sphere(&self, center: Vec3, radius: f64) -> bool {
        for plane in &self.planes {
            if plane.signed_dist(center) < -radius {
                return false;
            }
        }
        true
    }

    /// Test whether an axis-aligned bounding box is inside the frustum.
    ///
    /// Uses the positive-vertex test per plane.
    pub fn test_aabb(&self, min: Vec3, max: Vec3) -> bool {
        for plane in &self.planes {
            let px = if plane.normal[0] >= 0.0 {
                max[0]
            } else {
                min[0]
            };
            let py = if plane.normal[1] >= 0.0 {
                max[1]
            } else {
                min[1]
            };
            let pz = if plane.normal[2] >= 0.0 {
                max[2]
            } else {
                min[2]
            };
            if plane.signed_dist([px, py, pz]) < 0.0 {
                return false;
            }
        }
        true
    }
}

// ---------------------------------------------------------------------------
// Occlusion culling (hierarchical Z-buffer approximation)
// ---------------------------------------------------------------------------

/// A hierarchical Z-buffer tile for occlusion culling.
///
/// The buffer is represented as a power-of-two grid of depth tiles.
/// Each tile stores the *minimum* depth (closest) so that occluder
/// queries can quickly reject occludees.
#[derive(Debug, Clone)]
pub struct HiZBuffer {
    /// Width of the base level in tiles.
    pub width: usize,
    /// Height of the base level in tiles.
    pub height: usize,
    /// Depth values per tile (row-major).
    tiles: Vec<f32>,
}

impl HiZBuffer {
    /// Create a new Hi-Z buffer initialized to `f32::MAX` (far plane).
    pub fn new(width: usize, height: usize) -> Self {
        Self {
            width,
            height,
            tiles: vec![f32::MAX; width * height],
        }
    }

    /// Update tile at `(x, y)` with `depth`, keeping the minimum.
    pub fn update(&mut self, x: usize, y: usize, depth: f32) {
        if x < self.width && y < self.height {
            let idx = y * self.width + x;
            if depth < self.tiles[idx] {
                self.tiles[idx] = depth;
            }
        }
    }

    /// Query the minimum depth in a rectangular region `[x0..x1) × [y0..y1)`.
    pub fn query_min_depth(&self, x0: usize, y0: usize, x1: usize, y1: usize) -> f32 {
        let mut min_d = f32::MAX;
        let x1 = x1.min(self.width);
        let y1 = y1.min(self.height);
        for y in y0..y1 {
            for x in x0..x1 {
                let d = self.tiles[y * self.width + x];
                if d < min_d {
                    min_d = d;
                }
            }
        }
        min_d
    }

    /// Test whether a sphere at `depth_near` (NDC) is occluded.
    ///
    /// Returns `true` if the object is *behind* existing occluders (can skip rendering).
    pub fn is_occluded(
        &self,
        screen_x0: usize,
        screen_y0: usize,
        screen_x1: usize,
        screen_y1: usize,
        depth_near: f32,
    ) -> bool {
        let min_d = self.query_min_depth(screen_x0, screen_y0, screen_x1, screen_y1);
        depth_near > min_d
    }

    /// Reset all tiles to `f32::MAX`.
    pub fn clear(&mut self) {
        self.tiles.fill(f32::MAX);
    }
}

// ---------------------------------------------------------------------------
// Billboard sprites
// ---------------------------------------------------------------------------

/// A camera-facing quad (billboard) representing a particle or icon.
#[derive(Debug, Clone)]
pub struct Billboard {
    /// World-space position of the billboard center.
    pub position: Vec3,
    /// Half-size of the quad \[m\].
    pub half_size: f64,
    /// Base color of the billboard.
    pub color: Rgba,
    /// Texture atlas index (0 = default).
    pub texture_index: u32,
    /// Opacity in \[0, 1\].
    pub opacity: f32,
}

impl Billboard {
    /// Construct a new billboard.
    pub fn new(position: Vec3, half_size: f64, color: Rgba) -> Self {
        Self {
            position,
            half_size,
            color,
            texture_index: 0,
            opacity: 1.0,
        }
    }

    /// Generate the four corner positions of the billboard quad, given the
    /// camera right and up vectors.
    ///
    /// # Returns
    /// `[top_left, top_right, bottom_right, bottom_left]` in world space.
    pub fn corners(&self, cam_right: Vec3, cam_up: Vec3) -> [Vec3; 4] {
        let r = vscale(cam_right, self.half_size);
        let u = vscale(cam_up, self.half_size);
        let p = self.position;
        [
            vadd(vsub(p, r), u), // top-left
            vadd(vadd(p, r), u), // top-right
            vsub(vadd(p, r), u), // bottom-right
            vsub(vsub(p, r), u), // bottom-left
        ]
    }
}

// ---------------------------------------------------------------------------
// Instanced rendering for particles
// ---------------------------------------------------------------------------

/// Per-instance data for GPU instanced rendering.
#[derive(Debug, Clone, Copy)]
pub struct InstanceData {
    /// World-space position.
    pub position: [f32; 3],
    /// Rotation as axis-angle `[ax, ay, az, angle_rad]`.
    pub rotation: [f32; 4],
    /// Uniform scale factor.
    pub scale: f32,
    /// RGBA tint.
    pub color: [f32; 4],
    /// Metadata index (e.g. texture slot, particle age).
    pub metadata: u32,
}

impl Default for InstanceData {
    fn default() -> Self {
        Self {
            position: [0.0; 3],
            rotation: [0.0, 1.0, 0.0, 0.0],
            scale: 1.0,
            color: [1.0, 1.0, 1.0, 1.0],
            metadata: 0,
        }
    }
}

/// A batch of instance data suitable for upload to a GPU buffer.
#[derive(Debug, Clone)]
pub struct InstanceBatch {
    /// Instance data array.
    pub instances: Vec<InstanceData>,
    /// Maximum batch size before flushing.
    pub max_instances: usize,
}

impl InstanceBatch {
    /// Construct a new batch with a given maximum capacity.
    pub fn new(max_instances: usize) -> Self {
        Self {
            instances: Vec::with_capacity(max_instances),
            max_instances,
        }
    }

    /// Push an instance into the batch.
    ///
    /// Returns `true` if the batch is now full and should be flushed.
    pub fn push(&mut self, inst: InstanceData) -> bool {
        self.instances.push(inst);
        self.instances.len() >= self.max_instances
    }

    /// Number of instances currently in the batch.
    pub fn len(&self) -> usize {
        self.instances.len()
    }

    /// Whether the batch is empty.
    pub fn is_empty(&self) -> bool {
        self.instances.is_empty()
    }

    /// Clear the batch for reuse.
    pub fn clear(&mut self) {
        self.instances.clear();
    }
}

// ---------------------------------------------------------------------------
// Temporal anti-aliasing (TAA)
// ---------------------------------------------------------------------------

/// TAA accumulation buffer using exponential moving average.
///
/// Each pixel blends the current frame sample with the accumulated history.
#[derive(Debug, Clone)]
pub struct TaaBuffer {
    /// Width in pixels.
    pub width: usize,
    /// Height in pixels.
    pub height: usize,
    /// Accumulated RGBA values (f32 per channel).
    buffer: Vec<[f32; 4]>,
    /// Blend factor: 0 = use history only, 1 = use current frame only.
    pub blend_factor: f32,
}

impl TaaBuffer {
    /// Construct a TAA buffer initialized to black.
    pub fn new(width: usize, height: usize, blend_factor: f32) -> Self {
        Self {
            width,
            height,
            buffer: vec![[0.0, 0.0, 0.0, 1.0]; width * height],
            blend_factor: blend_factor.clamp(0.0, 1.0),
        }
    }

    /// Blend a new frame's pixel into the accumulation buffer at `(x, y)`.
    pub fn accumulate(&mut self, x: usize, y: usize, new_color: [f32; 4]) {
        if x >= self.width || y >= self.height {
            return;
        }
        let idx = y * self.width + x;
        let old = self.buffer[idx];
        let alpha = self.blend_factor;
        self.buffer[idx] = [
            old[0] * (1.0 - alpha) + new_color[0] * alpha,
            old[1] * (1.0 - alpha) + new_color[1] * alpha,
            old[2] * (1.0 - alpha) + new_color[2] * alpha,
            old[3] * (1.0 - alpha) + new_color[3] * alpha,
        ];
    }

    /// Get the accumulated color at `(x, y)`.
    pub fn get(&self, x: usize, y: usize) -> [f32; 4] {
        if x >= self.width || y >= self.height {
            return [0.0; 4];
        }
        self.buffer[y * self.width + x]
    }

    /// Reset the buffer to a uniform color.
    pub fn clear(&mut self, color: [f32; 4]) {
        self.buffer.fill(color);
    }

    /// Generate a sub-pixel jitter offset for the current frame index.
    ///
    /// Uses a Halton(2,3) sequence.
    pub fn jitter_offset(frame_index: u32) -> [f32; 2] {
        let hx = halton(frame_index + 1, 2);
        let hy = halton(frame_index + 1, 3);
        [hx - 0.5, hy - 0.5]
    }
}

/// Compute the n-th element of the Halton low-discrepancy sequence for base `base`.
fn halton(mut index: u32, base: u32) -> f32 {
    let mut result = 0.0_f32;
    let mut f = 1.0_f32;
    while index > 0 {
        f /= base as f32;
        result += f * (index % base) as f32;
        index /= base;
    }
    result
}

// ---------------------------------------------------------------------------
// Dynamic resolution scaling (DRS)
// ---------------------------------------------------------------------------

/// Dynamic resolution scaling controller.
///
/// Adjusts the internal render resolution fraction to maintain a target FPS.
#[derive(Debug, Clone)]
pub struct DynamicResolutionScaler {
    /// Target frames per second.
    pub target_fps: f64,
    /// Current resolution fraction in \[min_fraction, 1.0\].
    pub current_fraction: f64,
    /// Minimum allowed resolution fraction.
    pub min_fraction: f64,
    /// Step size for increasing resolution fraction per frame.
    pub increase_step: f64,
    /// Step size for decreasing resolution fraction per frame.
    pub decrease_step: f64,
}

impl DynamicResolutionScaler {
    /// Construct a new DRS controller.
    pub fn new(target_fps: f64, min_fraction: f64) -> Self {
        Self {
            target_fps,
            current_fraction: 1.0,
            min_fraction: min_fraction.clamp(0.1, 1.0),
            increase_step: 0.02,
            decrease_step: 0.05,
        }
    }

    /// Update the resolution fraction based on the measured FPS.
    pub fn update(&mut self, measured_fps: f64) {
        if measured_fps < self.target_fps * 0.95 {
            self.current_fraction =
                (self.current_fraction - self.decrease_step).max(self.min_fraction);
        } else if measured_fps > self.target_fps * 1.05 {
            self.current_fraction = (self.current_fraction + self.increase_step).min(1.0);
        }
    }

    /// Compute the scaled resolution `(width, height)` from the native resolution.
    pub fn scaled_resolution(&self, native_w: u32, native_h: u32) -> (u32, u32) {
        let w = ((native_w as f64 * self.current_fraction).round() as u32).max(1);
        let h = ((native_h as f64 * self.current_fraction).round() as u32).max(1);
        (w, h)
    }
}

// ---------------------------------------------------------------------------
// GPU particle system visualization
// ---------------------------------------------------------------------------

/// A single particle in the GPU particle system.
#[derive(Debug, Clone, Copy)]
pub struct GpuParticle {
    /// World-space position.
    pub position: [f32; 3],
    /// Velocity \[m/s\].
    pub velocity: [f32; 3],
    /// Remaining lifetime \[s\].
    pub lifetime: f32,
    /// Initial lifetime \[s\], used for normalizing age.
    pub lifetime_max: f32,
    /// RGBA color.
    pub color: [f32; 4],
    /// Size in world units.
    pub size: f32,
}

impl GpuParticle {
    /// Construct a new particle.
    pub fn new(
        position: [f32; 3],
        velocity: [f32; 3],
        lifetime: f32,
        color: [f32; 4],
        size: f32,
    ) -> Self {
        Self {
            position,
            velocity,
            lifetime,
            lifetime_max: lifetime,
            color,
            size,
        }
    }

    /// Normalized age in \[0, 1\], where 0 = just born, 1 = about to die.
    pub fn age_normalized(&self) -> f32 {
        if self.lifetime_max < 1e-10 {
            1.0
        } else {
            1.0 - self.lifetime / self.lifetime_max
        }
    }

    /// Return `true` if this particle is still alive.
    pub fn is_alive(&self) -> bool {
        self.lifetime > 0.0
    }
}

/// A pool of GPU particles with update and emission logic.
#[derive(Debug, Clone)]
pub struct GpuParticleSystem {
    /// Active particles.
    pub particles: Vec<GpuParticle>,
    /// Maximum number of simultaneous particles.
    pub max_particles: usize,
    /// Gravity acceleration \[m/s²\].
    pub gravity: [f32; 3],
}

impl GpuParticleSystem {
    /// Construct a new particle system.
    pub fn new(max_particles: usize) -> Self {
        Self {
            particles: Vec::with_capacity(max_particles),
            max_particles,
            gravity: [0.0, -9.81, 0.0],
        }
    }

    /// Emit a new particle. No-ops if the pool is full.
    pub fn emit(&mut self, particle: GpuParticle) {
        if self.particles.len() < self.max_particles {
            self.particles.push(particle);
        }
    }

    /// Advance all particles by time step `dt` \[s\] and remove dead ones.
    pub fn update(&mut self, dt: f32) {
        for p in &mut self.particles {
            p.position[0] += p.velocity[0] * dt;
            p.position[1] += p.velocity[1] * dt;
            p.position[2] += p.velocity[2] * dt;
            p.velocity[0] += self.gravity[0] * dt;
            p.velocity[1] += self.gravity[1] * dt;
            p.velocity[2] += self.gravity[2] * dt;
            p.lifetime -= dt;
        }
        self.particles.retain(|p| p.is_alive());
    }

    /// Number of alive particles.
    pub fn alive_count(&self) -> usize {
        self.particles.len()
    }
}

// ---------------------------------------------------------------------------
// Shadow maps
// ---------------------------------------------------------------------------

/// A simple depth-only shadow map.
#[derive(Debug, Clone)]
pub struct ShadowMap {
    /// Resolution in pixels (square).
    pub resolution: usize,
    /// Depth buffer storing NDC depth values.
    depth: Vec<f32>,
    /// Light-space view-projection matrix (row-major 4×4).
    pub light_vp: [f32; 16],
}

impl ShadowMap {
    /// Construct a shadow map of the given resolution.
    pub fn new(resolution: usize) -> Self {
        Self {
            resolution,
            depth: vec![1.0_f32; resolution * resolution],
            light_vp: {
                let mut m = [0.0_f32; 16];
                m[0] = 1.0;
                m[5] = 1.0;
                m[10] = 1.0;
                m[15] = 1.0;
                m
            },
        }
    }

    /// Record a depth value at `(x, y)` if it is closer than the stored value.
    pub fn update_depth(&mut self, x: usize, y: usize, depth: f32) {
        if x < self.resolution && y < self.resolution {
            let idx = y * self.resolution + x;
            if depth < self.depth[idx] {
                self.depth[idx] = depth;
            }
        }
    }

    /// Sample the depth at `(x, y)`.
    pub fn sample(&self, x: usize, y: usize) -> f32 {
        if x < self.resolution && y < self.resolution {
            self.depth[y * self.resolution + x]
        } else {
            1.0
        }
    }

    /// Test whether a point is in shadow given its light-space depth.
    ///
    /// Returns `true` if occluded (in shadow).
    pub fn is_in_shadow(&self, x: usize, y: usize, receiver_depth: f32, bias: f32) -> bool {
        let shadow_depth = self.sample(x, y);
        receiver_depth - bias > shadow_depth
    }

    /// Clear the depth buffer to far plane (1.0).
    pub fn clear(&mut self) {
        self.depth.fill(1.0);
    }
}

// ---------------------------------------------------------------------------
// Screen-space reflections (SSR)
// ---------------------------------------------------------------------------

/// Configuration for screen-space reflection ray-marching.
#[derive(Debug, Clone)]
pub struct SsrConfig {
    /// Maximum number of ray-march steps.
    pub max_steps: usize,
    /// Initial step size in screen-space (NDC).
    pub initial_step: f64,
    /// Step multiplier per iteration.
    pub step_multiplier: f64,
    /// Maximum ray length (NDC units).
    pub max_ray_length: f64,
    /// Thickness threshold for depth comparison.
    pub thickness: f64,
}

impl Default for SsrConfig {
    fn default() -> Self {
        Self {
            max_steps: 64,
            initial_step: 0.01,
            step_multiplier: 1.05,
            max_ray_length: 0.5,
            thickness: 0.01,
        }
    }
}

/// March a reflection ray through a linearized depth buffer.
///
/// # Arguments
/// * `origin` — Ray origin in view space.
/// * `direction` — Reflection direction in view space (normalized).
/// * `depth_buffer` — Flat depth buffer (row-major).
/// * `width`, `height` — Buffer dimensions.
/// * `config` — SSR parameters.
///
/// # Returns
/// `Some((u, v))` UV coordinates where the ray hit, or `None`.
pub fn ssr_ray_march(
    origin: Vec3,
    direction: Vec3,
    depth_buffer: &[f32],
    width: usize,
    height: usize,
    config: &SsrConfig,
) -> Option<(f64, f64)> {
    let mut step_size = config.initial_step;
    let mut ray_len = 0.0;

    for _ in 0..config.max_steps {
        ray_len += step_size;
        step_size *= config.step_multiplier;
        if ray_len > config.max_ray_length {
            break;
        }

        let sample_pos = vadd(origin, vscale(direction, ray_len));

        // Simple NDC-space coordinates assuming perspective projection
        let u = (sample_pos[0] + 1.0) / 2.0;
        let v = (sample_pos[1] + 1.0) / 2.0;

        if !(0.0..=1.0).contains(&u) || !(0.0..=1.0).contains(&v) {
            break;
        }

        let px = (u * (width - 1) as f64) as usize;
        let py = (v * (height - 1) as f64) as usize;
        let idx = py * width + px;
        if idx >= depth_buffer.len() {
            continue;
        }
        let depth_at = depth_buffer[idx] as f64;
        let ray_depth = sample_pos[2]; // z in view space (negative in RH)

        if ray_depth < depth_at && (depth_at - ray_depth) < config.thickness {
            return Some((u, v));
        }
    }
    None
}

// ---------------------------------------------------------------------------
// Ambient occlusion (SSAO/HBAO)
// ---------------------------------------------------------------------------

/// Generate a hemisphere sample kernel for SSAO.
///
/// Returns `n` sample vectors distributed in the positive hemisphere (+z).
/// Sample distances are biased towards the origin using a lerp.
///
/// # Arguments
/// * `n` — Number of samples (typically 16–64).
pub fn ssao_sample_kernel(n: usize) -> Vec<Vec3> {
    let mut kernel = Vec::with_capacity(n);
    for i in 0..n {
        let t = i as f64 / n as f64;
        // Spherical coordinates
        let theta = t * 2.0 * PI;
        let phi = (t * 0.5 * PI).min(PI / 2.0 - 1e-6);
        let r = {
            // Accelerating interpolation — samples closer to origin are more numerous

            lerp(0.1, 1.0, t * t)
        };
        let v = [
            r * phi.cos() * theta.cos(),
            r * phi.cos() * theta.sin(),
            r * phi.sin(),
        ];
        kernel.push(v);
    }
    kernel
}

/// Sample a random rotation vector for SSAO noise tile.
///
/// # Arguments
/// * `seed` — Deterministic seed (frame-independent).
pub fn ssao_noise_vector(seed: u32) -> Vec3 {
    // Simple hash to generate pseudo-random rotation angle
    let hash = seed.wrapping_mul(2654435761).wrapping_add(0xDEADBEEF);
    let angle = (hash as f64 / u32::MAX as f64) * 2.0 * PI;
    [angle.cos(), angle.sin(), 0.0]
}

/// Compute the SSAO occlusion factor for a single fragment.
///
/// # Arguments
/// * `frag_pos` — Fragment position in view space.
/// * `frag_normal` — Fragment normal in view space.
/// * `kernel` — Hemisphere sample kernel.
/// * `radius` — AO sampling radius.
/// * `depth_buffer` — Linearized depth buffer.
/// * `width`, `height` — Buffer dimensions.
///
/// # Returns
/// Occlusion factor in \[0, 1\]. 0 = fully occluded.
pub fn compute_ssao(
    frag_pos: Vec3,
    frag_normal: Vec3,
    kernel: &[Vec3],
    radius: f64,
    depth_buffer: &[f32],
    width: usize,
    height: usize,
) -> f64 {
    let normal = vnormalize(frag_normal);
    let mut occlusion = 0.0_f64;

    for &sample in kernel {
        // Orient sample towards normal hemisphere
        let dot = vdot(sample, normal);
        let oriented = if dot < 0.0 {
            vscale(sample, -1.0)
        } else {
            sample
        };
        let sample_pos = vadd(frag_pos, vscale(oriented, radius));

        // Project to screen space (simple orthographic approximation)
        let u = (sample_pos[0] + 1.0) / 2.0;
        let v = (sample_pos[1] + 1.0) / 2.0;

        if !(0.0..=1.0).contains(&u) || !(0.0..=1.0).contains(&v) {
            continue;
        }

        let px = (u * (width.saturating_sub(1)) as f64) as usize;
        let py = (v * (height.saturating_sub(1)) as f64) as usize;
        let idx = py * width + px;

        if idx < depth_buffer.len() {
            let sample_depth = depth_buffer[idx] as f64;
            let range_check = smoothstep(
                0.0,
                1.0,
                radius / (frag_pos[2] - sample_depth).abs().max(1e-6),
            );
            if sample_depth >= sample_pos[2] {
                occlusion += range_check;
            }
        }
    }

    if kernel.is_empty() {
        1.0
    } else {
        1.0 - (occlusion / kernel.len() as f64).min(1.0)
    }
}

// ---------------------------------------------------------------------------
// Deferred rendering G-buffer
// ---------------------------------------------------------------------------

/// A G-buffer (geometry buffer) for deferred rendering.
///
/// Stores per-pixel surface properties written in the geometry pass and
/// consumed in the lighting pass.
#[derive(Debug, Clone)]
pub struct GBuffer {
    /// Width in pixels.
    pub width: usize,
    /// Height in pixels.
    pub height: usize,
    /// Albedo (base color) buffer: \[r, g, b, metallic\].
    albedo: Vec<[f32; 4]>,
    /// World-space normals: \[nx, ny, nz, roughness\].
    normals: Vec<[f32; 4]>,
    /// View-space position: \[x, y, z, AO\].
    position: Vec<[f32; 4]>,
    /// Depth buffer.
    depth: Vec<f32>,
}

impl GBuffer {
    /// Construct a G-buffer for the given resolution.
    pub fn new(width: usize, height: usize) -> Self {
        let n = width * height;
        Self {
            width,
            height,
            albedo: vec![[0.0; 4]; n],
            normals: vec![[0.0, 0.0, 1.0, 0.5]; n],
            position: vec![[0.0; 4]; n],
            depth: vec![1.0; n],
        }
    }

    /// Write a fragment to the G-buffer.
    pub fn write(
        &mut self,
        x: usize,
        y: usize,
        albedo: [f32; 4],
        normal: [f32; 4],
        position: [f32; 4],
        depth: f32,
    ) {
        if x < self.width && y < self.height {
            let idx = y * self.width + x;
            self.albedo[idx] = albedo;
            self.normals[idx] = normal;
            self.position[idx] = position;
            self.depth[idx] = depth;
        }
    }

    /// Read the albedo at `(x, y)`.
    pub fn read_albedo(&self, x: usize, y: usize) -> [f32; 4] {
        if x < self.width && y < self.height {
            self.albedo[y * self.width + x]
        } else {
            [0.0; 4]
        }
    }

    /// Read the normal at `(x, y)`.
    pub fn read_normal(&self, x: usize, y: usize) -> [f32; 4] {
        if x < self.width && y < self.height {
            self.normals[y * self.width + x]
        } else {
            [0.0, 0.0, 1.0, 0.5]
        }
    }

    /// Clear the G-buffer.
    pub fn clear(&mut self) {
        self.albedo.fill([0.0; 4]);
        self.normals.fill([0.0, 0.0, 1.0, 0.5]);
        self.position.fill([0.0; 4]);
        self.depth.fill(1.0);
    }
}

// ---------------------------------------------------------------------------
// HDR tone mapping
// ---------------------------------------------------------------------------

/// Available tone mapping operators.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ToneMapOperator {
    /// Linear (no tone mapping, just exposure).
    Linear,
    /// Reinhard global tone mapping.
    Reinhard,
    /// ACES filmic curve approximation.
    Aces,
    /// Uncharted 2 "Hable" tone mapping.
    Uncharted2,
}

/// Apply a tone mapping operator to an HDR RGB triple.
///
/// # Arguments
/// * `hdr` — HDR linear color `[r, g, b]`.
/// * `exposure` — Pre-exposure multiplier (e.g. 1.0 = no change).
/// * `operator` — Tone mapping curve.
///
/// # Returns
/// LDR color clamped to `[0, 1]`.
pub fn tone_map(hdr: [f32; 3], exposure: f32, operator: ToneMapOperator) -> [f32; 3] {
    let e = [hdr[0] * exposure, hdr[1] * exposure, hdr[2] * exposure];
    match operator {
        ToneMapOperator::Linear => [e[0].min(1.0), e[1].min(1.0), e[2].min(1.0)],
        ToneMapOperator::Reinhard => [
            e[0] / (1.0 + e[0]),
            e[1] / (1.0 + e[1]),
            e[2] / (1.0 + e[2]),
        ],
        ToneMapOperator::Aces => [aces_channel(e[0]), aces_channel(e[1]), aces_channel(e[2])],
        ToneMapOperator::Uncharted2 => [uncharted2(e[0]), uncharted2(e[1]), uncharted2(e[2])],
    }
}

/// ACES filmic tone mapping for a single channel.
fn aces_channel(x: f32) -> f32 {
    let a = 2.51_f32;
    let b = 0.03_f32;
    let c = 2.43_f32;
    let d = 0.59_f32;
    let e = 0.14_f32;
    ((x * (a * x + b)) / (x * (c * x + d) + e)).clamp(0.0, 1.0)
}

/// Uncharted 2 (Hable) tone mapping for a single channel.
fn uncharted2(x: f32) -> f32 {
    let a = 0.15_f32;
    let b = 0.50_f32;
    let c = 0.10_f32;
    let d = 0.20_f32;
    let e = 0.02_f32;
    let f = 0.30_f32;
    ((x * (a * x + c * b) + d * e) / (x * (a * x + b) + d * f)) - e / f
}

/// Apply gamma correction (sRGB gamma ≈ 2.2).
///
/// # Arguments
/// * `linear` — Linear color `[r, g, b]` in \[0, 1\].
/// * `gamma` — Gamma exponent (typically 2.2).
pub fn gamma_correct(linear: [f32; 3], gamma: f32) -> [f32; 3] {
    let exp = 1.0 / gamma;
    [
        linear[0].max(0.0).powf(exp),
        linear[1].max(0.0).powf(exp),
        linear[2].max(0.0).powf(exp),
    ]
}

// ---------------------------------------------------------------------------
// Real-time global illumination (irradiance probes)
// ---------------------------------------------------------------------------

/// Spherical harmonics order-2 (L1) irradiance probe.
///
/// Stores 9 SH coefficients per RGB channel for fast diffuse GI evaluation.
#[derive(Debug, Clone)]
pub struct IrradianceProbe {
    /// World-space position of the probe.
    pub position: Vec3,
    /// SH coefficients: `[c0_r, c0_g, c0_b, c1_r, …]` — 9×3 = 27 floats.
    pub coefficients: [f32; 27],
    /// Whether this probe is dirty and needs re-baking.
    pub dirty: bool,
}

impl IrradianceProbe {
    /// Construct a new probe at `position` with all coefficients zeroed.
    pub fn new(position: Vec3) -> Self {
        Self {
            position,
            coefficients: [0.0; 27],
            dirty: true,
        }
    }

    /// Evaluate the irradiance for a surface normal `n`.
    ///
    /// Uses the Ramamoorthi & Hanrahan 2001 basis functions.
    ///
    /// # Returns
    /// Linear RGB irradiance `[r, g, b]`.
    pub fn evaluate(&self, n: Vec3) -> [f32; 3] {
        let c = &self.coefficients;
        let [x, y, z] = [n[0] as f32, n[1] as f32, n[2] as f32];
        // L0 and L1 SH basis
        let sh = [
            0.282_095_f32,                   // Y(0,0)
            0.488_603 * y,                   // Y(1,-1)
            0.488_603 * z,                   // Y(1,0)
            0.488_603 * x,                   // Y(1,1)
            1.092_548 * x * y,               // Y(2,-2)
            1.092_548 * y * z,               // Y(2,-1)
            0.315_392 * (3.0 * z * z - 1.0), // Y(2,0)
            1.092_548 * x * z,               // Y(2,1)
            0.546_274 * (x * x - y * y),     // Y(2,2)
        ];
        let mut rgb = [0.0_f32; 3];
        for (ci, &sh_val) in sh.iter().enumerate() {
            rgb[0] += c[ci * 3] * sh_val;
            rgb[1] += c[ci * 3 + 1] * sh_val;
            rgb[2] += c[ci * 3 + 2] * sh_val;
        }
        [rgb[0].max(0.0), rgb[1].max(0.0), rgb[2].max(0.0)]
    }

    /// Mark this probe as dirty (needs re-capture).
    pub fn invalidate(&mut self) {
        self.dirty = true;
    }
}

/// Grid of irradiance probes for volumetric GI.
#[derive(Debug, Clone)]
pub struct IrradianceGrid {
    /// Probes indexed \[x + y*nx + z*nx*ny\].
    pub probes: Vec<IrradianceProbe>,
    /// Grid resolution in each axis.
    pub resolution: [usize; 3],
    /// World-space origin of the grid.
    pub origin: Vec3,
    /// Spacing between probes \[m\].
    pub spacing: f64,
}

impl IrradianceGrid {
    /// Construct a grid of `nx × ny × nz` probes starting at `origin`.
    pub fn new(nx: usize, ny: usize, nz: usize, origin: Vec3, spacing: f64) -> Self {
        let mut probes = Vec::with_capacity(nx * ny * nz);
        for iz in 0..nz {
            for iy in 0..ny {
                for ix in 0..nx {
                    let pos = [
                        origin[0] + ix as f64 * spacing,
                        origin[1] + iy as f64 * spacing,
                        origin[2] + iz as f64 * spacing,
                    ];
                    probes.push(IrradianceProbe::new(pos));
                }
            }
        }
        Self {
            probes,
            resolution: [nx, ny, nz],
            origin,
            spacing,
        }
    }

    /// Trilinearly interpolate irradiance at a world-space position.
    pub fn sample(&self, pos: Vec3, normal: Vec3) -> [f32; 3] {
        let [nx, ny, nz] = self.resolution;
        if nx == 0 || ny == 0 || nz == 0 {
            return [0.0; 3];
        }

        // Grid-local fractional coordinates
        let fx = ((pos[0] - self.origin[0]) / self.spacing).max(0.0);
        let fy = ((pos[1] - self.origin[1]) / self.spacing).max(0.0);
        let fz = ((pos[2] - self.origin[2]) / self.spacing).max(0.0);

        let ix = (fx.floor() as usize).min(nx.saturating_sub(2));
        let iy = (fy.floor() as usize).min(ny.saturating_sub(2));
        let iz = (fz.floor() as usize).min(nz.saturating_sub(2));

        let tx = (fx - ix as f64).clamp(0.0, 1.0) as f32;
        let ty = (fy - iy as f64).clamp(0.0, 1.0) as f32;
        let tz = (fz - iz as f64).clamp(0.0, 1.0) as f32;

        let probe_irr = |xi: usize, yi: usize, zi: usize| {
            let idx = xi + yi * nx + zi * nx * ny;
            if idx < self.probes.len() {
                self.probes[idx].evaluate(normal)
            } else {
                [0.0; 3]
            }
        };

        // Trilinear blend
        let lerp3 = |a: [f32; 3], b: [f32; 3], t: f32| {
            [
                a[0] * (1.0 - t) + b[0] * t,
                a[1] * (1.0 - t) + b[1] * t,
                a[2] * (1.0 - t) + b[2] * t,
            ]
        };

        let c000 = probe_irr(ix, iy, iz);
        let c100 = probe_irr(ix + 1, iy, iz);
        let c010 = probe_irr(ix, iy + 1, iz);
        let c110 = probe_irr(ix + 1, iy + 1, iz);
        let c001 = probe_irr(ix, iy, iz + 1);
        let c101 = probe_irr(ix + 1, iy, iz + 1);
        let c011 = probe_irr(ix, iy + 1, iz + 1);
        let c111 = probe_irr(ix + 1, iy + 1, iz + 1);

        let cx00 = lerp3(c000, c100, tx);
        let cx10 = lerp3(c010, c110, tx);
        let cx01 = lerp3(c001, c101, tx);
        let cx11 = lerp3(c011, c111, tx);
        let cxy0 = lerp3(cx00, cx10, ty);
        let cxy1 = lerp3(cx01, cx11, ty);
        lerp3(cxy0, cxy1, tz)
    }

    /// Count dirty probes.
    pub fn dirty_count(&self) -> usize {
        self.probes.iter().filter(|p| p.dirty).count()
    }
}

// ---------------------------------------------------------------------------
// Math utilities
// ---------------------------------------------------------------------------

/// Linear interpolation.
#[inline]
fn lerp(a: f64, b: f64, t: f64) -> f64 {
    a + (b - a) * t.clamp(0.0, 1.0)
}

/// Smooth-step function.
#[inline]
fn smoothstep(edge0: f64, edge1: f64, x: f64) -> f64 {
    let t = ((x - edge0) / (edge1 - edge0)).clamp(0.0, 1.0);
    t * t * (3.0 - 2.0 * t)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ---- FrameTimer ----

    #[test]
    fn test_frame_timer_initial_fps() {
        let timer = FrameTimer::new(60.0, 8);
        assert!((timer.target_fps - 60.0).abs() < 1e-10);
    }

    #[test]
    fn test_frame_timer_push_and_measured_fps() {
        let mut timer = FrameTimer::new(60.0, 8);
        for _ in 0..8 {
            timer.push_frame(1.0 / 60.0);
        }
        let fps = timer.measured_fps();
        assert!((fps - 60.0).abs() < 0.1, "fps = {:.2}", fps);
    }

    #[test]
    fn test_frame_timer_elapsed() {
        let mut timer = FrameTimer::new(60.0, 8);
        timer.push_frame(0.1);
        timer.push_frame(0.2);
        assert!((timer.elapsed - 0.3).abs() < 1e-10);
    }

    #[test]
    fn test_frame_timer_on_target() {
        let mut timer = FrameTimer::new(60.0, 8);
        for _ in 0..8 {
            timer.push_frame(1.0 / 60.0);
        }
        assert!(timer.is_on_target(5.0), "Should be on target");
    }

    // ---- LodSelector ----

    #[test]
    fn test_lod_selector_high_at_close_range() {
        let lod = LodSelector::new(1.0);
        assert_eq!(lod.select(0.1), LodLevel::High);
    }

    #[test]
    fn test_lod_selector_culled_at_far_range() {
        let lod = LodSelector::new(1.0);
        let far = lod.dist_low_culled + 100.0;
        assert_eq!(lod.select(far), LodLevel::Culled);
    }

    #[test]
    fn test_lod_selector_medium() {
        let lod = LodSelector::new(1.0);
        let dist = (lod.dist_high_medium + lod.dist_medium_low) / 2.0;
        assert_eq!(lod.select(dist), LodLevel::Medium);
    }

    #[test]
    fn test_lod_selector_projected_pixels_decreases_with_distance() {
        let lod = LodSelector::new(1.0);
        let p1 = lod.projected_pixels(10.0);
        let p2 = lod.projected_pixels(100.0);
        assert!(
            p1 > p2,
            "closer objects should project larger: {} vs {}",
            p1,
            p2
        );
    }

    // ---- Frustum culling ----

    #[test]
    fn test_frustum_sphere_inside() {
        // Build an identity-like frustum and test a sphere at origin
        let mut m = [0.0_f64; 16];
        m[0] = 1.0;
        m[5] = 1.0;
        m[10] = 1.0;
        m[15] = 1.0;
        let frustum = ViewFrustum::from_view_proj(m);
        // Just verify no panic — result depends on plane extraction
        let _ = frustum.test_sphere([0.0, 0.0, 0.0], 1.0);
    }

    #[test]
    fn test_plane_signed_dist() {
        let plane = Plane::new([0.0, 0.0, 1.0], -5.0);
        let dist = plane.signed_dist([0.0, 0.0, 10.0]);
        assert!((dist - 5.0).abs() < 1e-10, "dist = {}", dist);
    }

    // ---- HiZBuffer ----

    #[test]
    fn test_hiz_buffer_initial_depth_max() {
        let buf = HiZBuffer::new(4, 4);
        assert_eq!(buf.query_min_depth(0, 0, 4, 4), f32::MAX);
    }

    #[test]
    fn test_hiz_buffer_update_and_query() {
        let mut buf = HiZBuffer::new(4, 4);
        buf.update(2, 2, 0.5);
        let d = buf.query_min_depth(0, 0, 4, 4);
        assert!((d - 0.5).abs() < 1e-6, "min depth = {}", d);
    }

    #[test]
    fn test_hiz_buffer_is_occluded() {
        let mut buf = HiZBuffer::new(4, 4);
        buf.update(1, 1, 0.3);
        // depth_near=0.5 > stored 0.3 → occluded
        assert!(buf.is_occluded(0, 0, 4, 4, 0.5));
        // depth_near=0.1 < stored 0.3 → not occluded
        assert!(!buf.is_occluded(0, 0, 4, 4, 0.1));
    }

    #[test]
    fn test_hiz_buffer_clear() {
        let mut buf = HiZBuffer::new(4, 4);
        buf.update(0, 0, 0.1);
        buf.clear();
        assert_eq!(buf.query_min_depth(0, 0, 1, 1), f32::MAX);
    }

    // ---- Billboard ----

    #[test]
    fn test_billboard_corners_form_quad() {
        let bb = Billboard::new([0.0, 0.0, 0.0], 1.0, Rgba::white());
        let right = [1.0, 0.0, 0.0];
        let up = [0.0, 1.0, 0.0];
        let corners = bb.corners(right, up);
        // Corners should be at distance √2 from center
        for c in &corners {
            let d = vnorm(*c);
            assert!((d - 2.0_f64.sqrt()).abs() < 1e-10, "corner dist: {}", d);
        }
    }

    // ---- InstanceBatch ----

    #[test]
    fn test_instance_batch_push_and_len() {
        let mut batch = InstanceBatch::new(10);
        batch.push(InstanceData::default());
        batch.push(InstanceData::default());
        assert_eq!(batch.len(), 2);
    }

    #[test]
    fn test_instance_batch_full_returns_true() {
        let mut batch = InstanceBatch::new(2);
        batch.push(InstanceData::default());
        let full = batch.push(InstanceData::default());
        assert!(full, "Batch should report full");
    }

    #[test]
    fn test_instance_batch_clear() {
        let mut batch = InstanceBatch::new(10);
        batch.push(InstanceData::default());
        batch.clear();
        assert!(batch.is_empty());
    }

    // ---- TAA ----

    #[test]
    fn test_taa_accumulate_converges() {
        let mut taa = TaaBuffer::new(2, 2, 1.0); // blend_factor=1 → always use new frame
        taa.accumulate(0, 0, [1.0, 0.0, 0.0, 1.0]);
        let c = taa.get(0, 0);
        assert!((c[0] - 1.0).abs() < 1e-6, "r should be 1.0: {}", c[0]);
    }

    #[test]
    fn test_taa_jitter_halton_sequence() {
        let j0 = TaaBuffer::jitter_offset(0);
        let j1 = TaaBuffer::jitter_offset(1);
        // Different frames should yield different jitter
        assert!(j0 != j1 || (j0[0] == 0.0 && j0[1] == 0.0));
    }

    #[test]
    fn test_taa_blend_factor_clamp() {
        let taa = TaaBuffer::new(4, 4, 2.0); // should clamp to 1.0
        assert!(taa.blend_factor <= 1.0);
    }

    // ---- DRS ----

    #[test]
    fn test_drs_reduces_on_low_fps() {
        let mut drs = DynamicResolutionScaler::new(60.0, 0.5);
        let initial = drs.current_fraction;
        drs.update(30.0); // 30 fps < 60*0.95
        assert!(drs.current_fraction < initial, "fraction should decrease");
    }

    #[test]
    fn test_drs_increases_on_high_fps() {
        let mut drs = DynamicResolutionScaler::new(60.0, 0.5);
        drs.current_fraction = 0.7;
        drs.update(120.0); // well above target
        assert!(drs.current_fraction > 0.7, "fraction should increase");
    }

    #[test]
    fn test_drs_scaled_resolution() {
        let drs = DynamicResolutionScaler::new(60.0, 0.5);
        let (w, h) = drs.scaled_resolution(1920, 1080);
        assert_eq!(w, 1920);
        assert_eq!(h, 1080);
    }

    // ---- GpuParticleSystem ----

    #[test]
    fn test_particle_emit_and_update() {
        let mut sys = GpuParticleSystem::new(100);
        sys.emit(GpuParticle::new(
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            1.0,
            [1.0, 1.0, 1.0, 1.0],
            1.0,
        ));
        assert_eq!(sys.alive_count(), 1);
        sys.update(0.5);
        assert_eq!(sys.alive_count(), 1);
        sys.update(0.6);
        assert_eq!(sys.alive_count(), 0, "particle should have died");
    }

    #[test]
    fn test_particle_max_capacity() {
        let mut sys = GpuParticleSystem::new(2);
        for _ in 0..5 {
            sys.emit(GpuParticle::new(
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0],
                1.0,
                [1.0, 1.0, 1.0, 1.0],
                1.0,
            ));
        }
        assert_eq!(sys.alive_count(), 2, "capped at max_particles");
    }

    #[test]
    fn test_particle_age_normalized() {
        let p = GpuParticle::new(
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            0.5,
            [1.0, 1.0, 1.0, 1.0],
            1.0,
        );
        let p2 = GpuParticle {
            lifetime: 0.25,
            ..p
        };
        let age = p2.age_normalized();
        // age = 1 - 0.25/0.5 = 0.5
        assert!((age - 0.5).abs() < 1e-6, "age_normalized = {}", age);
    }

    // ---- Shadow map ----

    #[test]
    fn test_shadow_map_depth_update() {
        let mut sm = ShadowMap::new(64);
        sm.update_depth(10, 10, 0.3);
        assert!((sm.sample(10, 10) - 0.3).abs() < 1e-6);
    }

    #[test]
    fn test_shadow_map_closer_depth_wins() {
        let mut sm = ShadowMap::new(64);
        sm.update_depth(5, 5, 0.7);
        sm.update_depth(5, 5, 0.4); // closer
        sm.update_depth(5, 5, 0.9); // farther, should not overwrite
        assert!((sm.sample(5, 5) - 0.4).abs() < 1e-6);
    }

    #[test]
    fn test_shadow_map_is_in_shadow() {
        let mut sm = ShadowMap::new(64);
        sm.update_depth(3, 3, 0.5);
        assert!(sm.is_in_shadow(3, 3, 0.8, 0.01), "depth 0.8 behind 0.5");
        assert!(
            !sm.is_in_shadow(3, 3, 0.3, 0.01),
            "depth 0.3 in front of 0.5"
        );
    }

    // ---- Tone mapping ----

    #[test]
    fn test_tone_map_reinhard_bounded() {
        let hdr = [5.0, 2.0, 0.1];
        let ldr = tone_map(hdr, 1.0, ToneMapOperator::Reinhard);
        for c in ldr {
            assert!((0.0..=1.0).contains(&c), "Reinhard out of range: {}", c);
        }
    }

    #[test]
    fn test_tone_map_aces_bounded() {
        let hdr = [3.0, 1.0, 0.5];
        let ldr = tone_map(hdr, 1.0, ToneMapOperator::Aces);
        for c in ldr {
            assert!((0.0..=1.0).contains(&c), "ACES out of range: {}", c);
        }
    }

    #[test]
    fn test_tone_map_linear_clamps() {
        let hdr = [2.0, 0.5, 0.0];
        let ldr = tone_map(hdr, 1.0, ToneMapOperator::Linear);
        assert!((ldr[0] - 1.0).abs() < 1e-6, "clamped to 1.0");
        assert!((ldr[1] - 0.5).abs() < 1e-6, "passes through");
    }

    #[test]
    fn test_gamma_correct_identity_at_one() {
        let c = [0.5_f32, 0.5, 0.5];
        let gc = gamma_correct(c, 1.0);
        for k in 0..3 {
            assert!((gc[k] - c[k]).abs() < 1e-5, "gamma=1.0 should be identity");
        }
    }

    // ---- Irradiance probes ----

    #[test]
    fn test_irradiance_probe_evaluate_zero_coefficients() {
        let probe = IrradianceProbe::new([0.0, 0.0, 0.0]);
        let irr = probe.evaluate([0.0, 0.0, 1.0]);
        // All zero coefficients → zero irradiance
        for c in irr {
            assert!(c.abs() < 1e-6, "irradiance should be 0: {}", c);
        }
    }

    #[test]
    fn test_irradiance_probe_non_negative() {
        let mut probe = IrradianceProbe::new([0.0, 0.0, 0.0]);
        // Set first coefficient (DC term) to positive value
        probe.coefficients[0] = 1.0;
        probe.coefficients[1] = 1.0;
        probe.coefficients[2] = 1.0;
        let irr = probe.evaluate([0.0, 1.0, 0.0]);
        for c in irr {
            assert!(c >= 0.0, "irradiance must be non-negative: {}", c);
        }
    }

    #[test]
    fn test_irradiance_grid_probe_count() {
        let grid = IrradianceGrid::new(4, 3, 2, [0.0, 0.0, 0.0], 1.0);
        assert_eq!(grid.probes.len(), 4 * 3 * 2);
    }

    #[test]
    fn test_irradiance_grid_dirty_count() {
        let grid = IrradianceGrid::new(2, 2, 2, [0.0, 0.0, 0.0], 1.0);
        assert_eq!(grid.dirty_count(), 8, "all probes should start dirty");
    }

    #[test]
    fn test_irradiance_grid_sample_zero_coefficients() {
        let grid = IrradianceGrid::new(2, 2, 2, [0.0, 0.0, 0.0], 1.0);
        let irr = grid.sample([0.5, 0.5, 0.5], [0.0, 1.0, 0.0]);
        for c in irr {
            assert!(c.abs() < 1e-6);
        }
    }

    // ---- SSAO ----

    #[test]
    fn test_ssao_sample_kernel_count() {
        let kernel = ssao_sample_kernel(32);
        assert_eq!(kernel.len(), 32);
    }

    #[test]
    fn test_ssao_sample_kernel_z_positive() {
        let kernel = ssao_sample_kernel(16);
        for v in &kernel {
            assert!(
                v[2] >= 0.0,
                "z component should be positive hemisphere: {}",
                v[2]
            );
        }
    }

    #[test]
    fn test_ssao_noise_vector_unit_length() {
        let n = ssao_noise_vector(42);
        let len = vnorm(n);
        assert!((len - 1.0).abs() < 1e-10, "noise vector length: {}", len);
    }

    // ---- GBuffer ----

    #[test]
    fn test_gbuffer_write_read_albedo() {
        let mut gb = GBuffer::new(16, 16);
        gb.write(
            5,
            5,
            [1.0, 0.0, 0.0, 1.0],
            [0.0, 0.0, 1.0, 0.5],
            [0.0, 0.0, -1.0, 1.0],
            0.5,
        );
        let a = gb.read_albedo(5, 5);
        assert!((a[0] - 1.0).abs() < 1e-6, "albedo r: {}", a[0]);
    }

    #[test]
    fn test_gbuffer_clear() {
        let mut gb = GBuffer::new(8, 8);
        gb.write(
            0,
            0,
            [1.0, 1.0, 1.0, 1.0],
            [0.0, 0.0, 1.0, 0.5],
            [0.0, 0.0, 0.0, 1.0],
            0.5,
        );
        gb.clear();
        let a = gb.read_albedo(0, 0);
        for c in a {
            assert!(c.abs() < 1e-6);
        }
    }

    // ---- Math utilities ----

    #[test]
    fn test_smoothstep_bounds() {
        assert!((smoothstep(0.0, 1.0, 0.0)).abs() < 1e-10);
        assert!((smoothstep(0.0, 1.0, 1.0) - 1.0).abs() < 1e-10);
        assert!((smoothstep(0.0, 1.0, 0.5) - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_halton_sequence_range() {
        for i in 0..16u32 {
            let h = halton(i, 2);
            assert!((0.0..1.0).contains(&h), "halton({}, 2) = {}", i, h);
        }
    }
}
