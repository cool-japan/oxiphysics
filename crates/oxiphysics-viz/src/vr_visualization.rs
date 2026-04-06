// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! VR/AR physics visualization.
//!
//! Provides stereo camera management, scene graph primitives, haptic-feedback
//! encoding, physics overlay helpers, and spatial audio utilities — all
//! implemented using only `std`, with no GPU dependency.

#![allow(dead_code)]

// ---------------------------------------------------------------------------
// Helper math (column-major 4×4 matrices and 3-D vector ops)
// ---------------------------------------------------------------------------

type Mat4 = [[f64; 4]; 4];
type Vec3 = [f64; 3];
type Quat = [f64; 4]; // [x, y, z, w]

fn mat4_identity() -> Mat4 {
    [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ]
}

fn mat4_mul(a: &Mat4, b: &Mat4) -> Mat4 {
    let mut out = [[0.0_f64; 4]; 4];
    for i in 0..4 {
        for j in 0..4 {
            for k in 0..4 {
                out[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    out
}

fn quat_to_mat4(q: &Quat) -> Mat4 {
    let [x, y, z, w] = *q;
    let x2 = x * x;
    let y2 = y * y;
    let z2 = z * z;
    let xy = x * y;
    let xz = x * z;
    let yz = y * z;
    let wx = w * x;
    let wy = w * y;
    let wz = w * z;
    [
        [1.0 - 2.0 * (y2 + z2), 2.0 * (xy - wz), 2.0 * (xz + wy), 0.0],
        [2.0 * (xy + wz), 1.0 - 2.0 * (x2 + z2), 2.0 * (yz - wx), 0.0],
        [2.0 * (xz - wy), 2.0 * (yz + wx), 1.0 - 2.0 * (x2 + y2), 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ]
}

fn vec3_sub(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn vec3_dot(a: Vec3, b: Vec3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn vec3_length(v: Vec3) -> f64 {
    vec3_dot(v, v).sqrt()
}

fn vec3_normalize(v: Vec3) -> Vec3 {
    let l = vec3_length(v);
    if l < 1e-15 {
        [0.0; 3]
    } else {
        [v[0] / l, v[1] / l, v[2] / l]
    }
}

fn vec3_scale(v: Vec3, s: f64) -> Vec3 {
    [v[0] * s, v[1] * s, v[2] * s]
}

fn vec3_add(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

// ---------------------------------------------------------------------------
// VrCamera
// ---------------------------------------------------------------------------

/// A single VR eye camera defined by position, orientation quaternion, and
/// perspective parameters.
#[derive(Debug, Clone)]
pub struct VrCamera {
    /// World-space position of the camera eye.
    pub position: Vec3,
    /// Orientation as a unit quaternion `[x, y, z, w]`.
    pub orientation: Quat,
    /// Vertical field of view in degrees.
    pub fov: f64,
    /// Near clip-plane distance (positive metres).
    pub near: f64,
    /// Far clip-plane distance (positive metres).
    pub far: f64,
}

impl VrCamera {
    /// Construct a default camera at the origin looking along −Z.
    pub fn new(position: Vec3, orientation: Quat, fov: f64, near: f64, far: f64) -> Self {
        Self {
            position,
            orientation,
            fov,
            near,
            far,
        }
    }

    /// Compute the world-to-camera view matrix (row-major 4×4).
    pub fn view_matrix(&self) -> Mat4 {
        // Rotation part (transpose of camera orientation matrix = inverse for unit quat)
        let rot = quat_to_mat4(&self.orientation);
        // Translation in camera space: -R^T * t
        let tx = -(rot[0][0] * self.position[0]
            + rot[1][0] * self.position[1]
            + rot[2][0] * self.position[2]);
        let ty = -(rot[0][1] * self.position[0]
            + rot[1][1] * self.position[1]
            + rot[2][1] * self.position[2]);
        let tz = -(rot[0][2] * self.position[0]
            + rot[1][2] * self.position[1]
            + rot[2][2] * self.position[2]);
        [
            [rot[0][0], rot[1][0], rot[2][0], 0.0],
            [rot[0][1], rot[1][1], rot[2][1], 0.0],
            [rot[0][2], rot[1][2], rot[2][2], 0.0],
            [tx, ty, tz, 1.0],
        ]
    }

    /// Compute the perspective projection matrix for this camera.
    ///
    /// `aspect` is width / height of the eye viewport.
    pub fn projection_matrix(&self, aspect: f64) -> Mat4 {
        let fov_rad = self.fov.to_radians();
        let f = 1.0 / (fov_rad / 2.0).tan();
        let near = self.near;
        let far = self.far;
        let nf = 1.0 / (near - far);
        [
            [f / aspect, 0.0, 0.0, 0.0],
            [0.0, f, 0.0, 0.0],
            [0.0, 0.0, (far + near) * nf, -1.0],
            [0.0, 0.0, 2.0 * far * near * nf, 0.0],
        ]
    }
}

// ---------------------------------------------------------------------------
// VrObject
// ---------------------------------------------------------------------------

/// A renderable / physics-simulated object in the VR scene.
#[derive(Debug, Clone)]
pub struct VrObject {
    /// 4×4 world transform matrix (row-major).
    pub transform: Mat4,
    /// Index into the scene's mesh array.
    pub mesh_id: usize,
    /// Index into the scene's material array.
    pub material_id: usize,
    /// Optional rigid-body physics handle.
    pub physics_body: Option<usize>,
}

impl VrObject {
    /// Construct a new `VrObject` at identity transform.
    pub fn new(mesh_id: usize, material_id: usize, physics_body: Option<usize>) -> Self {
        Self {
            transform: mat4_identity(),
            mesh_id,
            material_id,
            physics_body,
        }
    }

    /// Compute the world-space axis-aligned bounding box of this object.
    ///
    /// Returns `[min_corner, max_corner]` each as `[x, y, z]`.
    /// The unit-cube local AABB `([-0.5,-0.5,-0.5], [0.5,0.5,0.5])` is
    /// transformed by the object's world matrix.
    pub fn world_aabb(&self) -> [Vec3; 2] {
        let corners: [[f64; 3]; 8] = [
            [-0.5, -0.5, -0.5],
            [0.5, -0.5, -0.5],
            [-0.5, 0.5, -0.5],
            [0.5, 0.5, -0.5],
            [-0.5, -0.5, 0.5],
            [0.5, -0.5, 0.5],
            [-0.5, 0.5, 0.5],
            [0.5, 0.5, 0.5],
        ];
        let t = &self.transform;
        let transform_pt = |p: [f64; 3]| -> Vec3 {
            [
                t[0][0] * p[0] + t[1][0] * p[1] + t[2][0] * p[2] + t[3][0],
                t[0][1] * p[0] + t[1][1] * p[1] + t[2][1] * p[2] + t[3][1],
                t[0][2] * p[0] + t[1][2] * p[1] + t[2][2] * p[2] + t[3][2],
            ]
        };
        let mut mn = transform_pt(corners[0]);
        let mut mx = mn;
        for &c in corners.iter().skip(1) {
            let w = transform_pt(c);
            for i in 0..3 {
                mn[i] = mn[i].min(w[i]);
                mx[i] = mx[i].max(w[i]);
            }
        }
        [mn, mx]
    }
}

// ---------------------------------------------------------------------------
// VrScene
// ---------------------------------------------------------------------------

/// A stereoscopic VR scene containing objects and a stereo camera pair.
#[derive(Debug, Clone)]
pub struct VrScene {
    /// All objects in the scene.
    pub objects: Vec<VrObject>,
    /// Left eye camera.
    pub camera_left: VrCamera,
    /// Right eye camera.
    pub camera_right: VrCamera,
    /// Inter-pupillary distance in metres.
    pub ipd: f64,
}

impl VrScene {
    /// Create an empty VR scene with the given IPD.
    ///
    /// Both cameras are initialised at the origin looking down −Z with
    /// a 90° vertical FoV, 0.01 m near plane, and 1000 m far plane. The
    /// left and right cameras are offset laterally by ±`ipd/2`.
    pub fn new(ipd: f64) -> Self {
        let identity_quat: Quat = [0.0, 0.0, 0.0, 1.0];
        let half = ipd * 0.5;
        let camera_left = VrCamera::new([-half, 0.0, 0.0], identity_quat, 90.0, 0.01, 1000.0);
        let camera_right = VrCamera::new([half, 0.0, 0.0], identity_quat, 90.0, 0.01, 1000.0);
        Self {
            objects: Vec::new(),
            camera_left,
            camera_right,
            ipd,
        }
    }

    /// Add an object to the scene and return its index.
    pub fn add_object(&mut self, obj: VrObject) -> usize {
        let idx = self.objects.len();
        self.objects.push(obj);
        idx
    }

    /// Step the scene forward by `dt` seconds.
    ///
    /// This minimal implementation just propagates physics body handles
    /// (the actual simulation is performed by an external physics engine).
    pub fn update(&mut self, _dt: f64) {
        // Placeholder: physics integration would be driven externally.
    }

    /// Return the view matrices for the left and right eyes.
    pub fn stereo_view_matrices(&self) -> (Mat4, Mat4) {
        (
            self.camera_left.view_matrix(),
            self.camera_right.view_matrix(),
        )
    }

    /// Cull objects not visible from `cam` using a simple near-/far-plane test.
    ///
    /// Returns the indices of visible objects.
    pub fn frustum_cull(&self, cam: &VrCamera) -> Vec<usize> {
        let view = cam.view_matrix();
        self.objects
            .iter()
            .enumerate()
            .filter(|(_, obj)| {
                let aabb = obj.world_aabb();
                let center = [
                    (aabb[0][0] + aabb[1][0]) * 0.5,
                    (aabb[0][1] + aabb[1][1]) * 0.5,
                    (aabb[0][2] + aabb[1][2]) * 0.5,
                ];
                // Transform centre to camera space
                let cz = view[0][2] * center[0]
                    + view[1][2] * center[1]
                    + view[2][2] * center[2]
                    + view[3][2];
                // Visible if −far ≤ cz ≤ −near (right-hand view space)
                cz <= -cam.near && cz >= -cam.far
            })
            .map(|(i, _)| i)
            .collect()
    }

    /// Cast a ray from `origin` in direction `dir` and return the closest
    /// intersection `(object_index, t)`, or `None` if no hit.
    ///
    /// Uses a simple sphere-approximation per object (unit sphere).
    pub fn ray_cast(&self, origin: Vec3, dir: Vec3) -> Option<(usize, f64)> {
        let dir = vec3_normalize(dir);
        let mut closest: Option<(usize, f64)> = None;
        for (i, obj) in self.objects.iter().enumerate() {
            // Object centre from translation column of transform
            let c: Vec3 = [
                obj.transform[3][0],
                obj.transform[3][1],
                obj.transform[3][2],
            ];
            let oc = vec3_sub(origin, c);
            let b = 2.0 * vec3_dot(oc, dir);
            let c_val = vec3_dot(oc, oc) - 1.0; // unit sphere radius
            let disc = b * b - 4.0 * c_val;
            if disc < 0.0 {
                continue;
            }
            let t = (-b - disc.sqrt()) * 0.5;
            if t < 0.0 {
                continue;
            }
            if closest.is_none_or(|(_, ct)| t < ct) {
                closest = Some((i, t));
            }
        }
        closest
    }
}

// ---------------------------------------------------------------------------
// HapticFeedback
// ---------------------------------------------------------------------------

/// Haptic feedback command for a VR controller actuator.
#[derive(Debug, Clone)]
pub struct HapticFeedback {
    /// Normalised intensity `[0, 1]`.
    pub intensity: f64,
    /// Vibration frequency in Hz.
    pub frequency: f64,
    /// Duration in seconds.
    pub duration: f64,
}

impl HapticFeedback {
    /// Encode a physical contact force magnitude (Newtons) into a haptic pulse.
    ///
    /// Intensity is clamped to `[0, 1]` using a soft-saturation at 100 N.
    pub fn encode_force(force_magnitude: f64) -> Self {
        let intensity = (force_magnitude / 100.0).clamp(0.0, 1.0);
        // Map force to frequency: low force → 50 Hz, high force → 300 Hz
        let frequency = 50.0 + 250.0 * intensity;
        let duration = 0.02 + 0.08 * intensity; // 20–100 ms
        Self {
            intensity,
            frequency,
            duration,
        }
    }

    /// Generate a sequence of haptic pulses from a list of impact magnitudes.
    pub fn vibration_pattern(impacts: &[f64]) -> Vec<Self> {
        impacts.iter().map(|&f| Self::encode_force(f)).collect()
    }
}

// ---------------------------------------------------------------------------
// PhysicsOverlay
// ---------------------------------------------------------------------------

/// Configuration flags for the physics debug overlay.
#[derive(Debug, Clone)]
pub struct PhysicsOverlay {
    /// Whether to render force arrows.
    pub show_forces: bool,
    /// Whether to render velocity trail lines.
    pub show_velocities: bool,
    /// Whether to highlight contact points.
    pub show_contacts: bool,
}

impl Default for PhysicsOverlay {
    fn default() -> Self {
        Self {
            show_forces: true,
            show_velocities: true,
            show_contacts: true,
        }
    }
}

impl PhysicsOverlay {
    /// Construct a new `PhysicsOverlay` with all flags set to `enabled`.
    pub fn new(enabled: bool) -> Self {
        Self {
            show_forces: enabled,
            show_velocities: enabled,
            show_contacts: enabled,
        }
    }

    /// Compute an arrow primitive for a force vector.
    ///
    /// Returns `(tail_pos, tip_pos)` scaled by `scale`.
    pub fn force_arrow(&self, pos: Vec3, force: Vec3, scale: f64) -> (Vec3, Vec3) {
        let tip = vec3_add(pos, vec3_scale(force, scale));
        (pos, tip)
    }

    /// Compute a velocity trail from a list of past positions.
    ///
    /// Simply returns the positions clone as a render-ready line strip.
    pub fn velocity_trail(&self, positions: &[Vec3]) -> Vec<Vec3> {
        positions.to_vec()
    }
}

// ---------------------------------------------------------------------------
// AudioSource
// ---------------------------------------------------------------------------

/// A point sound emitter in 3-D space.
#[derive(Debug, Clone)]
pub struct AudioSource {
    /// World position of the emitter.
    pub position: Vec3,
    /// Base frequency in Hz.
    pub frequency: f64,
    /// Amplitude (linear gain).
    pub amplitude: f64,
}

impl AudioSource {
    /// Construct a new `AudioSource`.
    pub fn new(position: Vec3, frequency: f64, amplitude: f64) -> Self {
        Self {
            position,
            frequency,
            amplitude,
        }
    }
}

// ---------------------------------------------------------------------------
// SpatialAudio
// ---------------------------------------------------------------------------

/// Spatial audio engine with Doppler shift and distance attenuation.
#[derive(Debug, Clone)]
pub struct SpatialAudio {
    /// World position of the listener.
    pub listener_pos: Vec3,
    /// All active audio sources.
    pub sources: Vec<AudioSource>,
}

impl SpatialAudio {
    /// Construct a new `SpatialAudio` engine with the listener at `pos`.
    pub fn new(listener_pos: Vec3) -> Self {
        Self {
            listener_pos,
            sources: Vec::new(),
        }
    }

    /// Add an audio source and return its index.
    pub fn add_source(&mut self, source: AudioSource) -> usize {
        let idx = self.sources.len();
        self.sources.push(source);
        idx
    }

    /// Compute the Doppler-shifted frequency for a source.
    ///
    /// Uses the classical formula with the speed of sound at 343 m/s.
    ///
    /// * `source_vel`   — source velocity vector (m/s)
    /// * `listener_vel` — listener velocity vector (m/s)
    /// * `freq`         — emitted frequency (Hz)
    pub fn doppler_shift(&self, source_vel: Vec3, listener_vel: Vec3, freq: f64) -> f64 {
        let speed_of_sound = 343.0_f64;
        // Direction from source to listener (using first source if any, else origin)
        let src_pos = self.sources.first().map(|s| s.position).unwrap_or([0.0; 3]);
        let to_listener = vec3_sub(self.listener_pos, src_pos);
        let to_listener_n = vec3_normalize(to_listener);
        // Positive component = moving toward listener
        let vs = vec3_dot(source_vel, to_listener_n);
        let vl = vec3_dot(listener_vel, to_listener_n);
        let denom = speed_of_sound - vs;
        if denom.abs() < 1e-9 {
            return freq;
        }
        freq * (speed_of_sound + vl) / denom
    }

    /// Compute the gain attenuation for a source at `distance` metres using
    /// an inverse-distance law (clamped at 1 m minimum).
    pub fn attenuation(&self, distance: f64) -> f64 {
        let d = distance.max(1.0);
        1.0 / (d * d)
    }

    /// Compute the perceived frequency and amplitude for source `idx`.
    ///
    /// Returns `(freq, amplitude)`, or `(0, 0)` if index is out of range.
    pub fn perceived(&self, idx: usize, listener_vel: Vec3, source_vel: Vec3) -> (f64, f64) {
        if idx >= self.sources.len() {
            return (0.0, 0.0);
        }
        let src = &self.sources[idx];
        let dist = vec3_length(vec3_sub(self.listener_pos, src.position));
        let freq = self.doppler_shift(source_vel, listener_vel, src.frequency);
        let amp = src.amplitude * self.attenuation(dist);
        (freq, amp)
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // helpers
    fn identity_cam() -> VrCamera {
        VrCamera::new([0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 1.0], 90.0, 0.1, 100.0)
    }

    fn scene_with_object() -> (VrScene, usize) {
        let mut scene = VrScene::new(0.064);
        let obj = VrObject::new(0, 0, None);
        let idx = scene.add_object(obj);
        (scene, idx)
    }

    // --- VrCamera ---

    #[test]
    fn test_camera_view_matrix_is_4x4() {
        let cam = identity_cam();
        let m = cam.view_matrix();
        assert_eq!(m.len(), 4);
        assert_eq!(m[0].len(), 4);
    }

    #[test]
    fn test_camera_view_matrix_identity_orientation() {
        let cam = identity_cam();
        let m = cam.view_matrix();
        // With identity quaternion and zero position, view matrix ≈ identity
        assert!((m[0][0] - 1.0).abs() < 1e-9);
        assert!((m[1][1] - 1.0).abs() < 1e-9);
        assert!((m[2][2] - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_camera_projection_matrix_shape() {
        let cam = identity_cam();
        let p = cam.projection_matrix(1.0);
        assert_eq!(p.len(), 4);
    }

    #[test]
    fn test_camera_projection_fov_effect() {
        let cam_wide = VrCamera::new([0.0; 3], [0.0, 0.0, 0.0, 1.0], 120.0, 0.1, 100.0);
        let cam_narrow = VrCamera::new([0.0; 3], [0.0, 0.0, 0.0, 1.0], 45.0, 0.1, 100.0);
        let p_wide = cam_wide.projection_matrix(1.0);
        let p_narrow = cam_narrow.projection_matrix(1.0);
        // Wider FoV → smaller focal length (p[1][1])
        assert!(p_wide[1][1] < p_narrow[1][1]);
    }

    #[test]
    fn test_camera_projection_aspect_effect() {
        let cam = identity_cam();
        let p1 = cam.projection_matrix(1.0);
        let p2 = cam.projection_matrix(2.0);
        // Wider aspect → smaller p[0][0]
        assert!(p1[0][0] > p2[0][0]);
    }

    // --- VrObject ---

    #[test]
    fn test_object_world_aabb_identity_transform() {
        let obj = VrObject::new(0, 0, None);
        let aabb = obj.world_aabb();
        // Unit cube centred at origin
        assert!((aabb[0][0] - (-0.5)).abs() < 1e-9);
        assert!((aabb[1][0] - 0.5).abs() < 1e-9);
    }

    #[test]
    fn test_object_world_aabb_translated() {
        let mut obj = VrObject::new(0, 0, None);
        obj.transform[3][0] = 5.0; // translate X by 5
        let aabb = obj.world_aabb();
        assert!((aabb[0][0] - 4.5).abs() < 1e-9);
        assert!((aabb[1][0] - 5.5).abs() < 1e-9);
    }

    // --- VrScene ---

    #[test]
    fn test_scene_new_ipd() {
        let scene = VrScene::new(0.064);
        assert!((scene.ipd - 0.064).abs() < 1e-9);
    }

    #[test]
    fn test_scene_add_object_count() {
        let (scene, _) = scene_with_object();
        assert_eq!(scene.objects.len(), 1);
    }

    #[test]
    fn test_scene_stereo_view_matrices_differ() {
        let scene = VrScene::new(0.064);
        let (l, r) = scene.stereo_view_matrices();
        // Left and right translations should differ
        assert!((l[3][0] - r[3][0]).abs() > 1e-6);
    }

    #[test]
    fn test_scene_update_no_panic() {
        let mut scene = VrScene::new(0.064);
        scene.update(0.016); // should not panic
    }

    #[test]
    fn test_scene_frustum_cull_empty_scene() {
        let scene = VrScene::new(0.064);
        let visible = scene.frustum_cull(&scene.camera_left.clone());
        assert!(visible.is_empty());
    }

    #[test]
    fn test_scene_frustum_cull_object_behind_near() {
        let mut scene = VrScene::new(0.064);
        let mut obj = VrObject::new(0, 0, None);
        // Place object very close (within near plane) and directly in front
        obj.transform[3][2] = -0.001; // z = -0.001 but near = 0.1 → culled
        scene.add_object(obj);
        let visible = scene.frustum_cull(&scene.camera_left.clone());
        assert!(visible.is_empty());
    }

    #[test]
    fn test_scene_ray_cast_no_objects() {
        let scene = VrScene::new(0.064);
        assert!(scene.ray_cast([0.0, 0.0, 0.0], [0.0, 0.0, -1.0]).is_none());
    }

    #[test]
    fn test_scene_ray_cast_hits_object() {
        let mut scene = VrScene::new(0.064);
        let mut obj = VrObject::new(0, 0, None);
        obj.transform[3][2] = -5.0; // 5 m ahead
        scene.add_object(obj);
        let hit = scene.ray_cast([0.0, 0.0, 0.0], [0.0, 0.0, -1.0]);
        assert!(hit.is_some());
        let (idx, _t) = hit.unwrap();
        assert_eq!(idx, 0);
    }

    #[test]
    fn test_scene_ray_cast_misses_object() {
        let mut scene = VrScene::new(0.064);
        let mut obj = VrObject::new(0, 0, None);
        obj.transform[3][2] = -5.0;
        scene.add_object(obj);
        // Shoot ray sideways
        let hit = scene.ray_cast([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert!(hit.is_none());
    }

    // --- HapticFeedback ---

    #[test]
    fn test_haptic_encode_force_zero() {
        let h = HapticFeedback::encode_force(0.0);
        assert!((h.intensity - 0.0).abs() < 1e-9);
        assert!((h.frequency - 50.0).abs() < 1e-9);
    }

    #[test]
    fn test_haptic_encode_force_max() {
        let h = HapticFeedback::encode_force(1000.0);
        assert!((h.intensity - 1.0).abs() < 1e-9);
        assert!((h.frequency - 300.0).abs() < 1e-9);
    }

    #[test]
    fn test_haptic_encode_force_mid() {
        let h = HapticFeedback::encode_force(50.0);
        assert!((h.intensity - 0.5).abs() < 1e-9);
    }

    #[test]
    fn test_haptic_vibration_pattern_length() {
        let impacts = vec![10.0, 20.0, 30.0];
        let pattern = HapticFeedback::vibration_pattern(&impacts);
        assert_eq!(pattern.len(), 3);
    }

    #[test]
    fn test_haptic_vibration_pattern_increasing_intensity() {
        let pattern = HapticFeedback::vibration_pattern(&[10.0, 50.0, 100.0]);
        assert!(pattern[0].intensity < pattern[1].intensity);
        assert!(pattern[1].intensity < pattern[2].intensity);
    }

    // --- PhysicsOverlay ---

    #[test]
    fn test_overlay_force_arrow_tip_position() {
        let ov = PhysicsOverlay::default();
        let (tail, tip) = ov.force_arrow([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 1.0);
        assert!((tip[0] - 1.0).abs() < 1e-9);
        assert!((tail[0] - 0.0).abs() < 1e-9);
    }

    #[test]
    fn test_overlay_force_arrow_scaled() {
        let ov = PhysicsOverlay::default();
        let (_tail, tip) = ov.force_arrow([0.0; 3], [0.0, 1.0, 0.0], 3.0);
        assert!((tip[1] - 3.0).abs() < 1e-9);
    }

    #[test]
    fn test_overlay_velocity_trail_passthrough() {
        let ov = PhysicsOverlay::default();
        let pos = vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]];
        let trail = ov.velocity_trail(&pos);
        assert_eq!(trail.len(), 2);
        assert!((trail[1][2] - 6.0).abs() < 1e-9);
    }

    #[test]
    fn test_overlay_new_all_enabled() {
        let ov = PhysicsOverlay::new(true);
        assert!(ov.show_forces && ov.show_velocities && ov.show_contacts);
    }

    // --- SpatialAudio ---

    #[test]
    fn test_audio_attenuation_at_one_metre() {
        let audio = SpatialAudio::new([0.0; 3]);
        assert!((audio.attenuation(1.0) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_audio_attenuation_decreases_with_distance() {
        let audio = SpatialAudio::new([0.0; 3]);
        assert!(audio.attenuation(10.0) < audio.attenuation(2.0));
    }

    #[test]
    fn test_audio_doppler_stationary() {
        let mut audio = SpatialAudio::new([10.0, 0.0, 0.0]);
        audio.add_source(AudioSource::new([0.0; 3], 440.0, 1.0));
        // Both source and listener stationary → freq unchanged
        let shifted = audio.doppler_shift([0.0; 3], [0.0; 3], 440.0);
        assert!((shifted - 440.0).abs() < 1.0);
    }

    #[test]
    fn test_audio_doppler_source_approaching() {
        let mut audio = SpatialAudio::new([10.0, 0.0, 0.0]);
        audio.add_source(AudioSource::new([0.0; 3], 440.0, 1.0));
        // Source moving toward listener → frequency increases
        let shifted = audio.doppler_shift([10.0, 0.0, 0.0], [0.0; 3], 440.0);
        assert!(shifted > 440.0);
    }

    #[test]
    fn test_audio_add_source_returns_index() {
        let mut audio = SpatialAudio::new([0.0; 3]);
        let idx = audio.add_source(AudioSource::new([1.0, 0.0, 0.0], 220.0, 0.5));
        assert_eq!(idx, 0);
        let idx2 = audio.add_source(AudioSource::new([2.0, 0.0, 0.0], 440.0, 1.0));
        assert_eq!(idx2, 1);
    }

    #[test]
    fn test_audio_perceived_out_of_range() {
        let audio = SpatialAudio::new([0.0; 3]);
        let (f, a) = audio.perceived(99, [0.0; 3], [0.0; 3]);
        assert!((f - 0.0).abs() < 1e-9);
        assert!((a - 0.0).abs() < 1e-9);
    }
}
