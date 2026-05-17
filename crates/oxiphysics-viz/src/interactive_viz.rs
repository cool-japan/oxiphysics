#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Interactive visualization framework — camera control, scene objects,
//! picking, and animation timeline.
//!
//! No external crates are used; only `std`.

#![allow(dead_code)]

// ---------------------------------------------------------------------------
// Math helpers (std-only)
// ---------------------------------------------------------------------------

#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline]
fn len3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}

#[inline]
fn norm3(a: [f64; 3]) -> [f64; 3] {
    let l = len3(a);
    if l < 1e-15 { a } else { scale3(a, 1.0 / l) }
}

#[inline]
fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Multiply two 4x4 matrices (row-major).
fn mat4_mul(a: [[f64; 4]; 4], b: [[f64; 4]; 4]) -> [[f64; 4]; 4] {
    let mut c = [[0.0_f64; 4]; 4];
    for i in 0..4 {
        for j in 0..4 {
            for k in 0..4 {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    c
}

/// 4x4 identity matrix.
fn mat4_identity() -> [[f64; 4]; 4] {
    [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ]
}

// ---------------------------------------------------------------------------
// RenderCommand
// ---------------------------------------------------------------------------

/// A backend-agnostic rendering instruction emitted by the scene.
#[derive(Debug, Clone)]
pub enum RenderCommand {
    /// Draw a mesh identified by `id` with the given world transform.
    DrawMesh {
        /// Scene-object identifier.
        id: usize,
        /// World-space 4x4 transform (row-major).
        transform: [[f64; 4]; 4],
    },
    /// Draw a coloured line segment.
    DrawLine {
        /// World-space start of the line.
        start: [f64; 3],
        /// World-space end of the line.
        end: [f64; 3],
        /// RGBA colour.
        color: [f32; 4],
    },
    /// Render a text label at a world-space position.
    DrawText {
        /// World-space anchor point.
        pos: [f64; 3],
        /// Text content.
        text: String,
    },
}

// ---------------------------------------------------------------------------
// SimpleMesh
// ---------------------------------------------------------------------------

/// A minimal triangle mesh with per-vertex normals.
#[derive(Debug, Clone, Default)]
pub struct SimpleMesh {
    /// Vertex positions in local space.
    pub vertices: Vec<[f64; 3]>,
    /// Per-vertex normals (unit length).
    pub normals: Vec<[f64; 3]>,
    /// Triangle index triples.
    pub indices: Vec<[usize; 3]>,
}

impl SimpleMesh {
    /// Create an empty mesh.
    pub fn new() -> Self {
        Self::default()
    }

    /// Compute the bounding sphere as `(centre, radius)`.
    ///
    /// Returns `([0,0,0], 0.0)` for an empty mesh.
    pub fn bounding_sphere(&self) -> ([f64; 3], f64) {
        if self.vertices.is_empty() {
            return ([0.0; 3], 0.0);
        }
        // Centroid
        let mut cx = 0.0_f64;
        let mut cy = 0.0_f64;
        let mut cz = 0.0_f64;
        let n = self.vertices.len() as f64;
        for v in &self.vertices {
            cx += v[0];
            cy += v[1];
            cz += v[2];
        }
        let centre = [cx / n, cy / n, cz / n];
        // Max distance
        let radius = self
            .vertices
            .iter()
            .map(|v| len3(sub3(*v, centre)))
            .fold(0.0_f64, f64::max);
        (centre, radius)
    }

    /// Test a ray-mesh intersection using Möller–Trumbore.
    ///
    /// Returns `Some(t)` for the nearest intersection along the ray, or `None`
    /// if no triangle is hit.
    pub fn ray_intersect(&self, origin: [f64; 3], dir: [f64; 3]) -> Option<f64> {
        let eps = 1e-8;
        let mut nearest: Option<f64> = None;
        for tri in &self.indices {
            let v0 = self.vertices[tri[0]];
            let v1 = self.vertices[tri[1]];
            let v2 = self.vertices[tri[2]];
            let e1 = sub3(v1, v0);
            let e2 = sub3(v2, v0);
            let h = cross3(dir, e2);
            let a = dot3(e1, h);
            if a.abs() < eps {
                continue;
            }
            let f = 1.0 / a;
            let s = sub3(origin, v0);
            let u = f * dot3(s, h);
            if !(0.0..=1.0).contains(&u) {
                continue;
            }
            let q = cross3(s, e1);
            let v = f * dot3(dir, q);
            if v < 0.0 || u + v > 1.0 {
                continue;
            }
            let t = f * dot3(e2, q);
            if t > eps {
                nearest = Some(nearest.map_or(t, |prev: f64| prev.min(t)));
            }
        }
        nearest
    }
}

// ---------------------------------------------------------------------------
// SceneObject
// ---------------------------------------------------------------------------

/// An object in the interactive scene.
#[derive(Debug, Clone)]
pub struct SceneObject {
    /// Unique integer identifier within the scene.
    pub id: usize,
    /// Human-readable name.
    pub name: String,
    /// Triangle mesh geometry in local space.
    pub mesh: SimpleMesh,
    /// World-space 4x4 transform (row-major).
    pub transform: [[f64; 4]; 4],
    /// Whether this object is rendered.
    pub visible: bool,
    /// Whether this object is currently selected.
    pub selected: bool,
}

impl SceneObject {
    /// Create a new scene object with identity transform.
    pub fn new(id: usize, name: impl Into<String>, mesh: SimpleMesh) -> Self {
        Self {
            id,
            name: name.into(),
            mesh,
            transform: mat4_identity(),
            visible: true,
            selected: false,
        }
    }
}

// ---------------------------------------------------------------------------
// InteractiveCamera
// ---------------------------------------------------------------------------

/// A perspective camera supporting orbit, pan, and zoom controls.
#[derive(Debug, Clone)]
pub struct InteractiveCamera {
    /// Eye position in world space.
    pub position: [f64; 3],
    /// Look-at target in world space.
    pub target: [f64; 3],
    /// World-space up vector (should be unit length).
    pub up: [f64; 3],
    /// Vertical field of view in radians.
    pub fov: f64,
    /// Near clip plane distance.
    pub near: f64,
    /// Far clip plane distance.
    pub far: f64,
}

impl Default for InteractiveCamera {
    fn default() -> Self {
        Self {
            position: [0.0, 0.0, 5.0],
            target: [0.0; 3],
            up: [0.0, 1.0, 0.0],
            fov: std::f64::consts::FRAC_PI_4,
            near: 0.01,
            far: 1000.0,
        }
    }
}

impl InteractiveCamera {
    /// Create a camera with default parameters.
    pub fn new() -> Self {
        Self::default()
    }

    /// Orbit the camera around the target by `delta_theta` (horizontal, radians)
    /// and `delta_phi` (vertical, radians).
    pub fn orbit(&mut self, delta_theta: f64, delta_phi: f64) {
        let offset = sub3(self.position, self.target);
        let r = len3(offset);
        if r < 1e-12 {
            return;
        }
        // Spherical coordinates
        let mut theta = offset[0].atan2(offset[2]);
        let mut phi = (offset[1] / r).clamp(-1.0, 1.0).asin();
        theta += delta_theta;
        phi = (phi + delta_phi).clamp(
            -std::f64::consts::FRAC_PI_2 + 0.01,
            std::f64::consts::FRAC_PI_2 - 0.01,
        );
        self.position = add3(
            self.target,
            [
                r * phi.cos() * theta.sin(),
                r * phi.sin(),
                r * phi.cos() * theta.cos(),
            ],
        );
    }

    /// Pan the camera (move both position and target) by screen-relative offsets.
    pub fn pan(&mut self, dx: f64, dy: f64) {
        let forward = norm3(sub3(self.target, self.position));
        let right = norm3(cross3(forward, self.up));
        let up = norm3(cross3(right, forward));
        let delta = add3(scale3(right, -dx), scale3(up, dy));
        self.position = add3(self.position, delta);
        self.target = add3(self.target, delta);
    }

    /// Zoom by moving the camera along the view direction.
    ///
    /// Positive `delta` moves toward the target.
    pub fn zoom(&mut self, delta: f64) {
        let dir = norm3(sub3(self.target, self.position));
        let dist = len3(sub3(self.target, self.position));
        let new_dist = (dist - delta).max(self.near * 2.0);
        self.position = sub3(self.target, scale3(dir, new_dist));
    }

    /// Compute the view matrix (world → camera space).
    pub fn view_matrix(&self) -> [[f64; 4]; 4] {
        let f = norm3(sub3(self.target, self.position));
        let r = norm3(cross3(f, self.up));
        let u = cross3(r, f);
        [
            [r[0], u[0], -f[0], 0.0],
            [r[1], u[1], -f[1], 0.0],
            [r[2], u[2], -f[2], 0.0],
            [
                -dot3(r, self.position),
                -dot3(u, self.position),
                dot3(f, self.position),
                1.0,
            ],
        ]
    }

    /// Compute the symmetric perspective projection matrix for `aspect` = width/height.
    pub fn projection_matrix(&self, aspect: f64) -> [[f64; 4]; 4] {
        let tan_half = (self.fov * 0.5).tan();
        let f = 1.0 / tan_half;
        let nf = 1.0 / (self.near - self.far);
        [
            [f / aspect, 0.0, 0.0, 0.0],
            [0.0, f, 0.0, 0.0],
            [0.0, 0.0, (self.far + self.near) * nf, -1.0],
            [0.0, 0.0, 2.0 * self.far * self.near * nf, 0.0],
        ]
    }
}

// ---------------------------------------------------------------------------
// AnimationEntry
// ---------------------------------------------------------------------------

/// A tracked animation clip together with its current elapsed playback time.
///
/// Call [`InteractiveScene::add_animation`] to register a clip; the scene's
/// [`update`](InteractiveScene::update) method advances `elapsed` each frame.
#[derive(Debug, Clone)]
pub struct AnimationEntry {
    /// The animation clip data (duration, looping flag, keyframes).
    pub clip: crate::animation_system::AnimationClip,
    /// Current playback cursor in seconds.
    pub elapsed: f64,
    /// Whether the clip is actively advancing.
    pub playing: bool,
}

impl AnimationEntry {
    /// Wrap `clip` in a new entry, starting at time 0 and auto-playing.
    pub fn new(clip: crate::animation_system::AnimationClip) -> Self {
        Self {
            clip,
            elapsed: 0.0,
            playing: true,
        }
    }
}

// ---------------------------------------------------------------------------
// InteractiveScene
// ---------------------------------------------------------------------------

/// A scene containing objects that can be selected and rendered.
#[derive(Debug, Default)]
pub struct InteractiveScene {
    /// The active camera.
    pub camera: InteractiveCamera,
    /// All objects in the scene.
    pub objects: Vec<SceneObject>,
    /// Index of the currently selected object, if any.
    pub selection: Option<usize>,
    /// Active animation clips with per-clip playback cursors.
    pub animations: Vec<AnimationEntry>,
}

impl InteractiveScene {
    /// Create an empty scene with a default camera.
    pub fn new() -> Self {
        Self {
            camera: InteractiveCamera::new(),
            objects: Vec::new(),
            selection: None,
            animations: Vec::new(),
        }
    }

    /// Register an animation clip for playback; returns its index.
    pub fn add_animation(&mut self, clip: crate::animation_system::AnimationClip) -> usize {
        let idx = self.animations.len();
        self.animations.push(AnimationEntry::new(clip));
        idx
    }

    /// Add an object and return its assigned id.
    pub fn add_object(&mut self, mut obj: SceneObject) -> usize {
        let id = self.objects.len();
        obj.id = id;
        self.objects.push(obj);
        id
    }

    /// Cast a ray from `ray_origin` along `ray_dir` and select the nearest
    /// hit object.  Returns `Some(index)` on hit, `None` otherwise.
    pub fn select_at(&mut self, ray_origin: [f64; 3], ray_dir: [f64; 3]) -> Option<usize> {
        let mut nearest_t = f64::INFINITY;
        let mut nearest_idx: Option<usize> = None;
        for (idx, obj) in self.objects.iter().enumerate() {
            if !obj.visible {
                continue;
            }
            // Transform ray into local space (approximate: ignore rotation for now,
            // use world-space translation from column 3).
            let trans = [
                obj.transform[0][3],
                obj.transform[1][3],
                obj.transform[2][3],
            ];
            let local_origin = sub3(ray_origin, trans);
            if let Some(t) = obj.mesh.ray_intersect(local_origin, ray_dir)
                && t < nearest_t
            {
                nearest_t = t;
                nearest_idx = Some(idx);
            }
        }
        // Update selection state
        for obj in &mut self.objects {
            obj.selected = false;
        }
        if let Some(idx) = nearest_idx {
            self.objects[idx].selected = true;
        }
        self.selection = nearest_idx;
        nearest_idx
    }

    /// Advance the scene by `dt` seconds, stepping every registered animation clip.
    ///
    /// * Looping clips wrap their elapsed time back to zero using `rem_euclid`.
    /// * Non-looping clips clamp at their duration and set `playing = false`.
    pub fn update(&mut self, dt: f64) {
        for entry in &mut self.animations {
            if !entry.playing {
                continue;
            }
            entry.elapsed += dt;
            if entry.clip.looping {
                let dur = entry.clip.duration.max(1e-15);
                entry.elapsed = entry.elapsed.rem_euclid(dur);
            } else if entry.elapsed >= entry.clip.duration {
                entry.elapsed = entry.clip.duration;
                entry.playing = false;
            }
        }
    }

    /// Produce the list of render commands for the current frame.
    pub fn render_commands(&self) -> Vec<RenderCommand> {
        let mut cmds: Vec<RenderCommand> = Vec::new();
        for obj in &self.objects {
            if !obj.visible {
                continue;
            }
            cmds.push(RenderCommand::DrawMesh {
                id: obj.id,
                transform: obj.transform,
            });
            if obj.selected {
                // Draw a highlight line around the bounding sphere
                let (centre_local, radius) = obj.mesh.bounding_sphere();
                let world_centre = add3(
                    centre_local,
                    [
                        obj.transform[0][3],
                        obj.transform[1][3],
                        obj.transform[2][3],
                    ],
                );
                let end = add3(world_centre, [radius, 0.0, 0.0]);
                cmds.push(RenderCommand::DrawLine {
                    start: world_centre,
                    end,
                    color: [1.0, 1.0, 0.0, 1.0],
                });
            }
            cmds.push(RenderCommand::DrawText {
                pos: [
                    obj.transform[0][3],
                    obj.transform[1][3],
                    obj.transform[2][3],
                ],
                text: obj.name.clone(),
            });
        }
        cmds
    }
}

// ---------------------------------------------------------------------------
// Picker
// ---------------------------------------------------------------------------

/// Converts screen-space pixels to world-space rays for object picking.
#[derive(Debug, Clone)]
pub struct Picker {
    /// Framebuffer width in pixels.
    pub screen_width: u32,
    /// Framebuffer height in pixels.
    pub screen_height: u32,
}

impl Picker {
    /// Create a picker for the given framebuffer dimensions.
    pub fn new(screen_width: u32, screen_height: u32) -> Self {
        Self {
            screen_width,
            screen_height,
        }
    }

    /// Generate a world-space ray `(origin, direction)` passing through pixel `(x, y)`.
    ///
    /// Uses the inverse of the view-projection matrix derived from `cam`.
    pub fn ray_from_pixel(&self, x: u32, y: u32, cam: &InteractiveCamera) -> ([f64; 3], [f64; 3]) {
        let aspect = self.screen_width as f64 / self.screen_height as f64;
        // NDC in [-1, 1]
        let ndc_x = (x as f64 + 0.5) / self.screen_width as f64 * 2.0 - 1.0;
        let ndc_y = 1.0 - (y as f64 + 0.5) / self.screen_height as f64 * 2.0;

        let tan_half = (cam.fov * 0.5).tan();
        let ray_view = [ndc_x * aspect * tan_half, ndc_y * tan_half, -1.0];

        // Rotate into world space via the transpose of the rotation part of the view matrix
        let f = norm3(sub3(cam.target, cam.position));
        let r = norm3(cross3(f, cam.up));
        let u = cross3(r, f);

        let world_dir = norm3([
            r[0] * ray_view[0] + u[0] * ray_view[1] + (-f[0]) * ray_view[2],
            r[1] * ray_view[0] + u[1] * ray_view[1] + (-f[1]) * ray_view[2],
            r[2] * ray_view[0] + u[2] * ray_view[1] + (-f[2]) * ray_view[2],
        ]);

        (cam.position, world_dir)
    }
}

// ---------------------------------------------------------------------------
// AnimationTimeline
// ---------------------------------------------------------------------------

/// A simple linear animation timeline.
#[derive(Debug, Clone)]
pub struct AnimationTimeline {
    /// Total duration of the animation in seconds.
    pub duration: f64,
    /// Current playback time in seconds.
    pub current_time: f64,
    /// Whether the timeline is currently advancing.
    pub playing: bool,
}

impl AnimationTimeline {
    /// Create a new timeline of `duration` seconds, starting at t=0, paused.
    pub fn new(duration: f64) -> Self {
        Self {
            duration: duration.max(0.0),
            current_time: 0.0,
            playing: false,
        }
    }

    /// Start playback.
    pub fn play(&mut self) {
        self.playing = true;
    }

    /// Pause playback.
    pub fn pause(&mut self) {
        self.playing = false;
    }

    /// Jump to time `t` (clamped to `[0, duration]`).
    pub fn seek(&mut self, t: f64) {
        self.current_time = t.clamp(0.0, self.duration);
    }

    /// Advance the timeline by `dt` seconds if playing.
    ///
    /// Clamps at `duration` and stops playback when the end is reached.
    pub fn step(&mut self, dt: f64) {
        if self.playing {
            self.current_time = (self.current_time + dt).min(self.duration);
            if self.current_time >= self.duration {
                self.playing = false;
            }
        }
    }

    /// Return the current normalised time in `[0, 1]`.
    pub fn t(&self) -> f64 {
        if self.duration < 1e-15 {
            0.0
        } else {
            self.current_time / self.duration
        }
    }
}

// ---------------------------------------------------------------------------
// Mesh builders (convenience)
// ---------------------------------------------------------------------------

/// Build a simple tetrahedron `SimpleMesh`.
pub fn tetrahedron_mesh() -> SimpleMesh {
    let s = 1.0_f64;
    let h = (2.0_f64 / 3.0).sqrt();
    let vertices = vec![
        [0.0, s, 0.0],
        [-s * 0.5, -h * 0.5, s * (3.0_f64.sqrt() / 6.0)],
        [s * 0.5, -h * 0.5, s * (3.0_f64.sqrt() / 6.0)],
        [0.0, -h * 0.5, -s * (3.0_f64.sqrt() / 3.0)],
    ];
    let indices = vec![[0, 1, 2], [0, 2, 3], [0, 3, 1], [1, 3, 2]];
    let normals = vertices.iter().map(|v| norm3(*v)).collect();
    SimpleMesh {
        vertices,
        normals,
        indices,
    }
}

/// Build a unit cube `SimpleMesh` centred at the origin.
pub fn cube_mesh() -> SimpleMesh {
    let h = 0.5_f64;
    let vertices = vec![
        [-h, -h, -h],
        [h, -h, -h],
        [h, h, -h],
        [-h, h, -h],
        [-h, -h, h],
        [h, -h, h],
        [h, h, h],
        [-h, h, h],
    ];
    let indices = vec![
        [0, 1, 2],
        [0, 2, 3], // -z
        [4, 6, 5],
        [4, 7, 6], // +z
        [0, 5, 1],
        [0, 4, 5], // -y
        [2, 7, 3],
        [2, 6, 7], // +y
        [0, 3, 7],
        [0, 7, 4], // -x
        [1, 5, 6],
        [1, 6, 2], // +x
    ];
    let normals = vertices.iter().map(|v| norm3(*v)).collect();
    SimpleMesh {
        vertices,
        normals,
        indices,
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- AnimationTimeline ---

    #[test]
    fn test_timeline_new() {
        let tl = AnimationTimeline::new(10.0);
        assert!((tl.duration - 10.0).abs() < 1e-9);
        assert!((tl.current_time).abs() < 1e-9);
        assert!(!tl.playing);
    }

    #[test]
    fn test_timeline_play_and_step() {
        let mut tl = AnimationTimeline::new(5.0);
        tl.play();
        tl.step(1.0);
        assert!((tl.current_time - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_timeline_pause_stops_advance() {
        let mut tl = AnimationTimeline::new(5.0);
        tl.play();
        tl.pause();
        tl.step(2.0);
        assert!((tl.current_time).abs() < 1e-9);
    }

    #[test]
    fn test_timeline_seek_clamp() {
        let mut tl = AnimationTimeline::new(3.0);
        tl.seek(-5.0);
        assert!((tl.current_time).abs() < 1e-9);
        tl.seek(100.0);
        assert!((tl.current_time - 3.0).abs() < 1e-9);
    }

    #[test]
    fn test_timeline_t_normalised() {
        let mut tl = AnimationTimeline::new(10.0);
        tl.seek(5.0);
        assert!((tl.t() - 0.5).abs() < 1e-9);
    }

    #[test]
    fn test_timeline_stops_at_end() {
        let mut tl = AnimationTimeline::new(2.0);
        tl.play();
        tl.step(10.0);
        assert!(!tl.playing);
        assert!((tl.current_time - 2.0).abs() < 1e-9);
    }

    #[test]
    fn test_timeline_zero_duration() {
        let tl = AnimationTimeline::new(0.0);
        assert!((tl.t()).abs() < 1e-9);
    }

    // --- InteractiveCamera ---

    #[test]
    fn test_camera_default() {
        let cam = InteractiveCamera::default();
        assert!((cam.position[2] - 5.0).abs() < 1e-9);
        assert!((cam.target[0]).abs() < 1e-9);
    }

    #[test]
    fn test_camera_view_matrix_identity_like() {
        // Camera looking along -z from z=5 → view should translate the origin
        let cam = InteractiveCamera::default();
        let v = cam.view_matrix();
        // The (3,3) element of the view matrix (in column-major) should be 1.0
        assert!((v[3][3] - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_camera_projection_matrix_near_far() {
        let cam = InteractiveCamera::default();
        let p = cam.projection_matrix(1.0);
        // p[2][3] == -1 for the standard perspective matrix
        assert!((p[2][3] + 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_camera_zoom_reduces_distance() {
        let mut cam = InteractiveCamera::default();
        let before = len3(sub3(cam.target, cam.position));
        cam.zoom(1.0);
        let after = len3(sub3(cam.target, cam.position));
        assert!(after < before, "zoom should decrease distance to target");
    }

    #[test]
    fn test_camera_orbit_preserves_distance() {
        let mut cam = InteractiveCamera::default();
        let before = len3(sub3(cam.target, cam.position));
        cam.orbit(0.3, 0.1);
        let after = len3(sub3(cam.target, cam.position));
        assert!(
            (after - before).abs() < 1e-6,
            "orbit should preserve distance"
        );
    }

    #[test]
    fn test_camera_pan_shifts_both() {
        let mut cam = InteractiveCamera::default();
        let old_pos = cam.position;
        let old_tgt = cam.target;
        cam.pan(0.5, 0.0);
        // Both position and target should have moved by the same delta
        let dp = sub3(cam.position, old_pos);
        let dt = sub3(cam.target, old_tgt);
        for i in 0..3 {
            assert!(
                (dp[i] - dt[i]).abs() < 1e-9,
                "pan: position and target should shift identically"
            );
        }
    }

    // --- SimpleMesh ---

    #[test]
    fn test_simple_mesh_bounding_sphere_empty() {
        let mesh = SimpleMesh::new();
        let (c, r) = mesh.bounding_sphere();
        assert_eq!(r, 0.0);
        assert_eq!(c, [0.0; 3]);
    }

    #[test]
    fn test_simple_mesh_bounding_sphere_tetra() {
        let mesh = tetrahedron_mesh();
        let (_c, r) = mesh.bounding_sphere();
        assert!(r > 0.0, "bounding sphere radius should be positive");
    }

    #[test]
    fn test_simple_mesh_ray_intersect_cube_hit() {
        let mesh = cube_mesh();
        let origin = [0.0, 0.0, 5.0];
        let dir = [0.0, 0.0, -1.0];
        let hit = mesh.ray_intersect(origin, dir);
        assert!(hit.is_some(), "ray should hit the cube");
        assert!(hit.unwrap() > 0.0);
    }

    #[test]
    fn test_simple_mesh_ray_intersect_miss() {
        let mesh = cube_mesh();
        let origin = [10.0, 0.0, 0.0];
        let dir = [1.0, 0.0, 0.0]; // shooting away from cube
        let hit = mesh.ray_intersect(origin, dir);
        assert!(hit.is_none(), "ray shooting away should miss");
    }

    // --- SceneObject ---

    #[test]
    fn test_scene_object_default_visible() {
        let obj = SceneObject::new(0, "box", cube_mesh());
        assert!(obj.visible);
        assert!(!obj.selected);
    }

    #[test]
    fn test_scene_object_transform_identity() {
        let obj = SceneObject::new(1, "tet", tetrahedron_mesh());
        let eye = mat4_identity();
        for i in 0..4 {
            for j in 0..4 {
                assert!((obj.transform[i][j] - eye[i][j]).abs() < 1e-9);
            }
        }
    }

    // --- InteractiveScene ---

    #[test]
    fn test_scene_add_object() {
        let mut scene = InteractiveScene::new();
        let id = scene.add_object(SceneObject::new(0, "a", SimpleMesh::new()));
        assert_eq!(id, 0);
        assert_eq!(scene.objects.len(), 1);
    }

    #[test]
    fn test_scene_render_commands_visible() {
        let mut scene = InteractiveScene::new();
        scene.add_object(SceneObject::new(0, "visible", cube_mesh()));
        let cmds = scene.render_commands();
        assert!(!cmds.is_empty());
    }

    #[test]
    fn test_scene_render_commands_hidden() {
        let mut scene = InteractiveScene::new();
        let mut obj = SceneObject::new(0, "hidden", cube_mesh());
        obj.visible = false;
        scene.add_object(obj);
        let cmds = scene.render_commands();
        // DrawMesh commands should be absent for invisible objects
        let draw_mesh_count = cmds
            .iter()
            .filter(|c| matches!(c, RenderCommand::DrawMesh { .. }))
            .count();
        assert_eq!(draw_mesh_count, 0);
    }

    #[test]
    fn test_scene_select_hit() {
        let mut scene = InteractiveScene::new();
        scene.add_object(SceneObject::new(0, "cube", cube_mesh()));
        let origin = [0.0, 0.0, 5.0];
        let dir = [0.0, 0.0, -1.0];
        let sel = scene.select_at(origin, dir);
        assert!(sel.is_some());
    }

    #[test]
    fn test_scene_select_miss() {
        let mut scene = InteractiveScene::new();
        scene.add_object(SceneObject::new(0, "cube", cube_mesh()));
        let origin = [100.0, 100.0, 100.0];
        let dir = [1.0, 0.0, 0.0];
        let sel = scene.select_at(origin, dir);
        assert!(sel.is_none());
    }

    #[test]
    fn test_scene_update_noop() {
        let mut scene = InteractiveScene::new();
        scene.update(0.016);
        // just verify no panic
    }

    // --- Picker ---

    #[test]
    fn test_picker_center_pixel_ray() {
        let picker = Picker::new(800, 600);
        let cam = InteractiveCamera::default();
        let (origin, dir) = picker.ray_from_pixel(400, 300, &cam);
        // Origin should be at camera position
        assert!((origin[2] - 5.0).abs() < 1e-9);
        // Direction should point roughly towards -z (the target)
        assert!(dir[2] < 0.0, "centre ray should point toward scene");
    }

    #[test]
    fn test_picker_direction_is_unit() {
        let picker = Picker::new(1024, 768);
        let cam = InteractiveCamera::default();
        let (_o, dir) = picker.ray_from_pixel(512, 384, &cam);
        let length = len3(dir);
        assert!((length - 1.0).abs() < 1e-6, "ray direction should be unit");
    }

    #[test]
    fn test_picker_different_pixels_different_rays() {
        let picker = Picker::new(800, 600);
        let cam = InteractiveCamera::default();
        let (_o1, d1) = picker.ray_from_pixel(100, 100, &cam);
        let (_o2, d2) = picker.ray_from_pixel(700, 500, &cam);
        let diff: f64 = (0..3).map(|i| (d1[i] - d2[i]).abs()).sum();
        assert!(
            diff > 1e-4,
            "different pixels should give different ray directions"
        );
    }

    // --- RenderCommand ---

    #[test]
    fn test_render_command_draw_line_debug() {
        let cmd = RenderCommand::DrawLine {
            start: [0.0; 3],
            end: [1.0, 0.0, 0.0],
            color: [1.0, 0.0, 0.0, 1.0],
        };
        let s = format!("{cmd:?}");
        assert!(s.contains("DrawLine"));
    }

    #[test]
    fn test_render_command_draw_text_debug() {
        let cmd = RenderCommand::DrawText {
            pos: [0.0; 3],
            text: "hello".into(),
        };
        let s = format!("{cmd:?}");
        assert!(s.contains("DrawText"));
        assert!(s.contains("hello"));
    }

    // --- mat4 helpers ---

    #[test]
    fn test_mat4_identity_mul() {
        let eye = mat4_identity();
        let result = mat4_mul(eye, eye);
        for i in 0..4 {
            for j in 0..4 {
                assert!((result[i][j] - eye[i][j]).abs() < 1e-9);
            }
        }
    }

    // --- Animation update (F2) ---

    fn make_clip(duration: f64, looping: bool) -> crate::animation_system::AnimationClip {
        crate::animation_system::AnimationClip::new("test", duration, 24.0, looping)
    }

    #[test]
    fn test_update_advances_clip_elapsed() {
        let mut scene = InteractiveScene::new();
        let idx = scene.add_animation(make_clip(2.0, false));
        scene.update(0.5);
        assert!((scene.animations[idx].elapsed - 0.5).abs() < 1e-12);
        scene.update(0.3);
        assert!((scene.animations[idx].elapsed - 0.8).abs() < 1e-12);
    }

    #[test]
    fn test_update_looping_clip_wraps() {
        let mut scene = InteractiveScene::new();
        let idx = scene.add_animation(make_clip(1.0, true));
        // advance past the full duration
        scene.update(1.5);
        // elapsed should wrap: 1.5 rem_euclid 1.0 == 0.5
        assert!(
            (scene.animations[idx].elapsed - 0.5).abs() < 1e-12,
            "looping elapsed = {}",
            scene.animations[idx].elapsed
        );
        assert!(
            scene.animations[idx].playing,
            "looping clip must remain playing"
        );
    }

    #[test]
    fn test_update_nonlooping_clip_stops() {
        let mut scene = InteractiveScene::new();
        let idx = scene.add_animation(make_clip(1.0, false));
        scene.update(2.0);
        assert!(
            (scene.animations[idx].elapsed - 1.0).abs() < 1e-12,
            "non-looping elapsed clamped to duration"
        );
        assert!(!scene.animations[idx].playing, "non-looping clip must stop");
        // further updates should not change elapsed
        scene.update(1.0);
        assert!((scene.animations[idx].elapsed - 1.0).abs() < 1e-12);
    }
}
