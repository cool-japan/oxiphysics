//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;

/// Depth of field configuration for post-processing.
#[derive(Debug, Clone, Copy)]
pub struct DepthOfField {
    /// Distance to the focal plane (sharp focus).
    pub focal_distance: f64,
    /// Range around focal distance where objects are in focus.
    pub focal_range: f64,
    /// Maximum blur radius in pixels.
    pub max_blur: f64,
}
impl DepthOfField {
    /// Create a new depth of field configuration.
    pub fn new(focal_distance: f64, focal_range: f64, max_blur: f64) -> Self {
        Self {
            focal_distance,
            focal_range,
            max_blur,
        }
    }
    /// Compute the blur factor (0..1) for a fragment at the given depth.
    ///
    /// 0.0 means fully in focus, 1.0 means maximum blur.
    pub fn blur_factor(&self, depth: f64) -> f64 {
        let dist_from_focal = (depth - self.focal_distance).abs();
        let half_range = self.focal_range * 0.5;
        if dist_from_focal <= half_range {
            0.0
        } else {
            ((dist_from_focal - half_range) / self.focal_distance.max(1e-6)).min(1.0)
        }
    }
    /// Compute the circle of confusion radius for a given depth.
    pub fn coc_radius(&self, depth: f64) -> f64 {
        self.blur_factor(depth) * self.max_blur
    }
}
/// Type of cinematic camera move.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CinematicMoveType {
    /// Linear dolly (forward/backward along the look axis).
    Dolly,
    /// Vertical crane (up/down).
    Crane,
    /// Horizontal jib (side-to-side).
    Jib,
    /// Combined arc (orbit around target at fixed elevation).
    Arc,
}
/// A multi-target camera that automatically frames a set of world-space points.
#[derive(Debug, Clone)]
pub struct MultiTargetCamera {
    /// Base camera whose position/target will be updated.
    pub camera: Camera,
    /// List of world-space positions to frame.
    pub targets: Vec<[f64; 3]>,
    /// Extra padding added to the radius (world units).
    pub padding: f64,
    /// Smoothing factor per frame (0 = instant, 1 = never moves).
    pub smoothing: f64,
    /// Current smooth position.
    pub(super) current_position: [f64; 3],
    /// Current smooth target.
    pub(super) current_target: [f64; 3],
}
impl MultiTargetCamera {
    /// Create a new multi-target camera from an existing [`Camera`].
    pub fn new(camera: Camera) -> Self {
        let pos = camera.position;
        let tgt = camera.target;
        Self {
            camera,
            targets: Vec::new(),
            padding: 1.0,
            smoothing: 0.1,
            current_position: pos,
            current_target: tgt,
        }
    }
    /// Add a world-space point to track.
    pub fn add_target(&mut self, point: [f64; 3]) {
        self.targets.push(point);
    }
    /// Remove all tracked points.
    pub fn clear_targets(&mut self) {
        self.targets.clear();
    }
    /// Compute the centroid and bounding sphere radius of tracked targets.
    ///
    /// Returns `None` when no targets are registered.
    pub fn bounding_sphere(&self) -> Option<([f64; 3], f64)> {
        if self.targets.is_empty() {
            return None;
        }
        let n = self.targets.len() as f64;
        let cx = self.targets.iter().map(|t| t[0]).sum::<f64>() / n;
        let cy = self.targets.iter().map(|t| t[1]).sum::<f64>() / n;
        let cz = self.targets.iter().map(|t| t[2]).sum::<f64>() / n;
        let centroid = [cx, cy, cz];
        let radius = self
            .targets
            .iter()
            .map(|t| {
                let d = sub3(*t, centroid);
                dot3(d, d).sqrt()
            })
            .fold(0.0_f64, f64::max)
            + self.padding;
        Some((centroid, radius))
    }
    /// Update the camera to frame all targets, applying `smoothing`.
    pub fn update(&mut self, dt: f64) {
        let Some((centroid, radius)) = self.bounding_sphere() else {
            return;
        };
        let half_fov = self.camera.fov_y * 0.5;
        let required_dist = if half_fov.abs() > 1e-10 {
            radius / half_fov.tan()
        } else {
            radius * 10.0
        };
        let offset_dir = normalize3(sub3(self.camera.position, centroid));
        let offset_dir = if dot3(offset_dir, offset_dir) < 1e-12 {
            [0.0, 0.0, 1.0]
        } else {
            offset_dir
        };
        let desired_pos = add3(centroid, scale3(offset_dir, required_dist));
        let alpha = (1.0 - self.smoothing.clamp(0.0, 0.9999)).powf(dt * 60.0);
        self.current_position = lerp3(desired_pos, self.current_position, alpha);
        self.current_target = lerp3(centroid, self.current_target, alpha);
        self.camera.position = self.current_position;
        self.camera.target = self.current_target;
    }
    /// Current number of tracked targets.
    pub fn target_count(&self) -> usize {
        self.targets.len()
    }
}
/// A first-person fly camera that moves freely in 3D space.
///
/// Controlled by yaw/pitch angles and forward/right/up movement.
pub struct FlyCamera {
    /// World-space position.
    pub position: [f64; 3],
    /// Yaw angle in radians (rotation around Y axis).
    pub yaw: f64,
    /// Pitch angle in radians (rotation around right axis).
    pub pitch: f64,
    /// Movement speed (units per second).
    pub speed: f64,
    /// Mouse/rotation sensitivity.
    pub sensitivity: f64,
    /// Vertical FOV in radians.
    pub fov_y: f64,
    /// Aspect ratio.
    pub aspect: f64,
    /// Near plane.
    pub near: f64,
    /// Far plane.
    pub far: f64,
}
impl FlyCamera {
    /// Create a new fly camera at the given position.
    pub fn new(position: [f64; 3], fov_y_deg: f64, aspect: f64) -> Self {
        Self {
            position,
            yaw: -std::f64::consts::FRAC_PI_2,
            pitch: 0.0,
            speed: 5.0,
            sensitivity: 0.003,
            fov_y: fov_y_deg.to_radians(),
            aspect,
            near: 0.1,
            far: 1000.0,
        }
    }
    /// Forward direction based on yaw and pitch.
    pub fn forward(&self) -> [f64; 3] {
        let cp = self.pitch.cos();
        [self.yaw.cos() * cp, self.pitch.sin(), self.yaw.sin() * cp]
    }
    /// Right direction (perpendicular to forward on XZ plane).
    pub fn right(&self) -> [f64; 3] {
        let yaw_r = self.yaw + std::f64::consts::FRAC_PI_2;
        [yaw_r.cos(), 0.0, yaw_r.sin()]
    }
    /// Rotate the camera by mouse deltas.
    pub fn rotate(&mut self, dx: f64, dy: f64) {
        self.yaw += dx * self.sensitivity;
        self.pitch = (self.pitch - dy * self.sensitivity).clamp(
            -std::f64::consts::FRAC_PI_2 + 0.01,
            std::f64::consts::FRAC_PI_2 - 0.01,
        );
    }
    /// Move forward/backward by `amount * speed`.
    pub fn move_forward(&mut self, amount: f64) {
        let fwd = self.forward();
        for (p, f) in self.position.iter_mut().zip(fwd.iter()) {
            *p += f * amount * self.speed;
        }
    }
    /// Strafe left/right by `amount * speed`.
    pub fn move_right(&mut self, amount: f64) {
        let r = self.right();
        for (p, r_i) in self.position.iter_mut().zip(r.iter()) {
            *p += r_i * amount * self.speed;
        }
    }
    /// Move up/down along world Y axis by `amount * speed`.
    pub fn move_up(&mut self, amount: f64) {
        self.position[1] += amount * self.speed;
    }
    /// Compute a view matrix from the fly camera state.
    pub fn view_matrix(&self) -> [[f64; 4]; 4] {
        let fwd = self.forward();
        let target = [
            self.position[0] + fwd[0],
            self.position[1] + fwd[1],
            self.position[2] + fwd[2],
        ];
        let cam = Camera::new(
            self.position,
            target,
            self.fov_y.to_degrees(),
            self.aspect,
            self.near,
            self.far,
        );
        cam.view_matrix()
    }
    /// Compute a projection matrix.
    pub fn projection_matrix(&self) -> [[f64; 4]; 4] {
        let cam = Camera::new(
            self.position,
            [0.0, 0.0, 0.0],
            self.fov_y.to_degrees(),
            self.aspect,
            self.near,
            self.far,
        );
        cam.projection_matrix()
    }
}
/// Temporal Anti-Aliasing (TAA) history buffer.
///
/// Stores the last N camera poses for jitter-blending across frames.
#[derive(Debug, Clone)]
pub struct TaaHistory {
    /// Ring buffer of previous camera positions.
    pub positions: Vec<[f64; 3]>,
    /// Ring buffer of previous camera targets.
    pub targets: Vec<[f64; 3]>,
    /// Maximum number of history frames stored.
    pub capacity: usize,
    /// Current write index (wraps around).
    pub(super) write_idx: usize,
    /// Number of valid entries currently stored.
    pub(super) count: usize,
}
impl TaaHistory {
    /// Create a new TAA history with the given capacity.
    pub fn new(capacity: usize) -> Self {
        let cap = capacity.max(1);
        Self {
            positions: vec![[0.0; 3]; cap],
            targets: vec![[0.0; 3]; cap],
            capacity: cap,
            write_idx: 0,
            count: 0,
        }
    }
    /// Push the current camera pose into the history.
    pub fn push(&mut self, position: [f64; 3], target: [f64; 3]) {
        self.positions[self.write_idx] = position;
        self.targets[self.write_idx] = target;
        self.write_idx = (self.write_idx + 1) % self.capacity;
        self.count = (self.count + 1).min(self.capacity);
    }
    /// Number of valid history entries (up to `capacity`).
    pub fn len(&self) -> usize {
        self.count
    }
    /// Returns `true` if no history has been recorded yet.
    pub fn is_empty(&self) -> bool {
        self.count == 0
    }
    /// Compute the exponentially-weighted average position for blending.
    ///
    /// `alpha` is the blend factor per frame (0 = all history, 1 = current only).
    pub fn blended_position(&self, alpha: f64) -> Option<[f64; 3]> {
        if self.count == 0 {
            return None;
        }
        let mut acc = [0.0_f64; 3];
        let mut weight_sum = 0.0_f64;
        let mut w = 1.0_f64;
        for k in 0..self.count {
            let idx = (self.write_idx + self.capacity - 1 - k) % self.capacity;
            let pos = self.positions[idx];
            acc[0] += pos[0] * w;
            acc[1] += pos[1] * w;
            acc[2] += pos[2] * w;
            weight_sum += w;
            w *= 1.0 - alpha;
        }
        Some([
            acc[0] / weight_sum,
            acc[1] / weight_sum,
            acc[2] / weight_sum,
        ])
    }
    /// Sub-pixel jitter offset for TAA (Halton sequence, 2 bases).
    ///
    /// Returns `(jitter_x, jitter_y)` in NDC-space (each in –1..1).
    pub fn halton_jitter(frame: usize, width: usize, height: usize) -> (f64, f64) {
        let jx = halton(frame + 1, 2) - 0.5;
        let jy = halton(frame + 1, 3) - 0.5;
        let w = width.max(1) as f64;
        let h = height.max(1) as f64;
        (jx * 2.0 / w, jy * 2.0 / h)
    }
}
/// A cinematic camera move definition.
#[derive(Debug, Clone)]
pub struct CinematicMove {
    /// The type of move.
    pub move_type: CinematicMoveType,
    /// Total distance / angle to travel over the duration.
    pub amount: f64,
    /// Duration of the move in seconds.
    pub duration: f64,
    /// Elapsed time within this move (0..duration).
    pub(super) elapsed: f64,
}
impl CinematicMove {
    /// Create a new cinematic move.
    pub fn new(move_type: CinematicMoveType, amount: f64, duration: f64) -> Self {
        Self {
            move_type,
            amount,
            duration,
            elapsed: 0.0,
        }
    }
    /// Advance the move by `dt` seconds and apply to `camera`.
    ///
    /// Returns `true` when the move has completed.
    pub fn advance(&mut self, camera: &mut Camera, dt: f64) -> bool {
        if self.elapsed >= self.duration {
            return true;
        }
        let remaining = (self.duration - self.elapsed).min(dt);
        let fraction = if self.duration > 1e-12 {
            remaining / self.duration
        } else {
            1.0
        };
        let delta = self.amount * fraction;
        match self.move_type {
            CinematicMoveType::Dolly => camera.zoom(delta),
            CinematicMoveType::Crane => {
                camera.position[1] += delta;
                camera.target[1] += delta;
            }
            CinematicMoveType::Jib => camera.pan(delta, 0.0),
            CinematicMoveType::Arc => camera.orbit(delta, 0.0),
        }
        self.elapsed += dt;
        self.elapsed >= self.duration
    }
    /// Progress fraction in \[0, 1\].
    pub fn progress(&self) -> f64 {
        if self.duration < 1e-12 {
            1.0
        } else {
            (self.elapsed / self.duration).min(1.0)
        }
    }
    /// `true` when completed.
    pub fn is_done(&self) -> bool {
        self.elapsed >= self.duration
    }
}
/// A camera rig that sequences multiple [`CinematicMove`]s.
#[derive(Debug, Clone)]
pub struct CameraRig {
    /// The camera being controlled.
    pub camera: Camera,
    /// The pending move queue.
    pub moves: std::collections::VecDeque<CinematicMove>,
}
impl CameraRig {
    /// Create a new rig wrapping the given camera.
    pub fn new(camera: Camera) -> Self {
        Self {
            camera,
            moves: std::collections::VecDeque::new(),
        }
    }
    /// Enqueue a cinematic move.
    pub fn enqueue(&mut self, m: CinematicMove) {
        self.moves.push_back(m);
    }
    /// Advance the rig by `dt` seconds.
    ///
    /// Executes the front-most move and pops it when complete.
    pub fn update(&mut self, dt: f64) {
        if let Some(front) = self.moves.front_mut()
            && front.advance(&mut self.camera, dt)
        {
            self.moves.pop_front();
        }
    }
    /// Number of queued moves.
    pub fn pending_count(&self) -> usize {
        self.moves.len()
    }
    /// `true` when all moves have completed.
    pub fn is_idle(&self) -> bool {
        self.moves.is_empty()
    }
}
/// A keyframe for camera animation.
#[derive(Debug, Clone)]
pub struct CameraKeyframe {
    /// Time at which this keyframe occurs.
    pub time: f64,
    /// Camera position at this keyframe.
    pub position: [f64; 3],
    /// Camera target at this keyframe.
    pub target: [f64; 3],
}
/// A perspective camera defined by position, target, and projection parameters.
#[derive(Debug, Clone)]
pub struct Camera {
    /// World-space eye position.
    pub position: [f64; 3],
    /// World-space look-at target.
    pub target: [f64; 3],
    /// World-space up hint vector.
    pub up: [f64; 3],
    /// Vertical field of view in radians.
    pub fov_y: f64,
    /// Aspect ratio width/height.
    pub aspect: f64,
    /// Near clip distance.
    pub near: f64,
    /// Far clip distance.
    pub far: f64,
}
impl Camera {
    /// Create a new camera. `fov_y_deg` is in **degrees**.
    pub fn new(
        position: [f64; 3],
        target: [f64; 3],
        fov_y_deg: f64,
        aspect: f64,
        near: f64,
        far: f64,
    ) -> Self {
        Self {
            position,
            target,
            up: [0.0, 1.0, 0.0],
            fov_y: fov_y_deg.to_radians(),
            aspect,
            near,
            far,
        }
    }
    /// LookAt view matrix (column-major, OpenGL convention).
    pub fn view_matrix(&self) -> [[f64; 4]; 4] {
        let f = normalize3(sub3(self.target, self.position));
        let r = normalize3(cross3(f, self.up));
        let u = cross3(r, f);
        let tx = -dot3(r, self.position);
        let ty = -dot3(u, self.position);
        let tz = dot3(f, self.position);
        [
            [r[0], u[0], -f[0], 0.0],
            [r[1], u[1], -f[1], 0.0],
            [r[2], u[2], -f[2], 0.0],
            [tx, ty, tz, 1.0],
        ]
    }
    /// Perspective projection matrix (column-major, OpenGL convention, depth -1..1).
    pub fn projection_matrix(&self) -> [[f64; 4]; 4] {
        let tan_half = (self.fov_y / 2.0).tan();
        let f = 1.0 / tan_half;
        let n = self.near;
        let fa = self.far;
        let range = fa - n;
        [
            [f / self.aspect, 0.0, 0.0, 0.0],
            [0.0, f, 0.0, 0.0],
            [0.0, 0.0, -(fa + n) / range, -1.0],
            [0.0, 0.0, -2.0 * fa * n / range, 0.0],
        ]
    }
    /// Combined projection * view matrix.
    pub fn view_projection(&self) -> [[f64; 4]; 4] {
        mat4_mul(self.projection_matrix(), self.view_matrix())
    }
    /// Forward direction (normalised, toward target).
    pub fn forward(&self) -> [f64; 3] {
        normalize3(sub3(self.target, self.position))
    }
    /// Right direction (normalised).
    pub fn right(&self) -> [f64; 3] {
        normalize3(cross3(self.forward(), self.up))
    }
    /// True up direction (camera local up, normalised).
    pub fn up_vector(&self) -> [f64; 3] {
        cross3(self.right(), self.forward())
    }
    /// Orbit the camera around the target by `delta_yaw` and `delta_pitch` (radians).
    pub fn orbit(&mut self, delta_yaw: f64, delta_pitch: f64) {
        let mut offset = sub3(self.position, self.target);
        let radius = (dot3(offset, offset)).sqrt();
        let theta = offset[1].atan2((offset[0] * offset[0] + offset[2] * offset[2]).sqrt());
        let phi = offset[2].atan2(offset[0]);
        let new_theta = (theta + delta_pitch).clamp(
            -std::f64::consts::FRAC_PI_2 + 1e-4,
            std::f64::consts::FRAC_PI_2 - 1e-4,
        );
        let new_phi = phi + delta_yaw;
        offset[0] = radius * new_theta.cos() * new_phi.cos();
        offset[1] = radius * new_theta.sin();
        offset[2] = radius * new_theta.cos() * new_phi.sin();
        self.position = [
            self.target[0] + offset[0],
            self.target[1] + offset[1],
            self.target[2] + offset[2],
        ];
    }
    /// Move the camera along the forward axis by `delta`.
    pub fn zoom(&mut self, delta: f64) {
        let fwd = self.forward();
        self.position = [
            self.position[0] + fwd[0] * delta,
            self.position[1] + fwd[1] * delta,
            self.position[2] + fwd[2] * delta,
        ];
        self.target = [
            self.target[0] + fwd[0] * delta,
            self.target[1] + fwd[1] * delta,
            self.target[2] + fwd[2] * delta,
        ];
    }
    /// Distance from the camera eye to the target.
    pub fn distance_to_target(&self) -> f64 {
        let d = sub3(self.position, self.target);
        (dot3(d, d)).sqrt()
    }
    /// Pan the camera (translate position and target) by `dx` in the right direction
    /// and `dy` in the up direction.
    pub fn pan(&mut self, dx: f64, dy: f64) {
        let r = self.right();
        let u = self.up_vector();
        for i in 0..3 {
            let delta = r[i] * dx + u[i] * dy;
            self.position[i] += delta;
            self.target[i] += delta;
        }
    }
    /// Set the camera to look at a target from a given distance along the -Z axis.
    pub fn look_at_distance(&mut self, target: [f64; 3], distance: f64) {
        self.target = target;
        self.position = [target[0], target[1], target[2] + distance];
    }
}
/// Perlin-noise camera shake — produces smoother, more cinematic shake than
/// the simple sinusoidal version.
pub struct PerlinShake {
    /// Shake amplitude in world units.
    pub amplitude: f64,
    /// Decay rate per second (0..1, higher = faster fade).
    pub decay: f64,
    /// Current elapsed time (drives the noise).
    pub(super) time: f64,
    /// Shake frequency in Hz.
    pub frequency: f64,
}
impl PerlinShake {
    /// Create a new Perlin-noise shake.
    pub fn new(amplitude: f64, decay: f64, frequency: f64) -> Self {
        Self {
            amplitude,
            decay,
            time: 0.0,
            frequency,
        }
    }
    /// Advance time by `dt` seconds and decay amplitude.
    pub fn update(&mut self, dt: f64) {
        self.time += dt * self.frequency;
        self.amplitude *= (1.0 - self.decay * dt).max(0.0);
    }
    /// Compute the current positional offset using a simple smooth-noise
    /// approximation (fade + trilinear interpolation on a 1-D lattice).
    pub fn offset(&self) -> [f64; 3] {
        [
            self.amplitude * smooth_noise_1d(self.time * 1.0),
            self.amplitude * smooth_noise_1d(self.time * 1.3 + 31.41),
            self.amplitude * smooth_noise_1d(self.time * 0.7 + 17.07),
        ]
    }
    /// `true` when amplitude has fallen below a negligible threshold.
    pub fn is_done(&self) -> bool {
        self.amplitude < 1e-6
    }
    /// Add an impulse to the shake amplitude.
    pub fn impulse(&mut self, amount: f64) {
        self.amplitude += amount;
    }
    /// Current amplitude value.
    pub fn amplitude(&self) -> f64 {
        self.amplitude
    }
}
/// Simple camera shake state.
pub struct CameraShake {
    /// Current shake intensity (amplitude).
    pub(super) intensity: f64,
    /// Decay rate per second (how fast shake fades out).
    pub(super) decay: f64,
    /// Current phase (for oscillation).
    pub(super) phase: f64,
    /// Frequency of oscillation (Hz).
    pub(super) frequency: f64,
}
impl CameraShake {
    /// Create a new camera shake with the given initial intensity.
    pub fn new(intensity: f64, decay: f64, frequency: f64) -> Self {
        Self {
            intensity,
            decay,
            phase: 0.0,
            frequency,
        }
    }
    /// Update the shake state by the given delta time (seconds).
    pub fn update(&mut self, dt: f64) {
        self.phase += dt * self.frequency * 2.0 * std::f64::consts::PI;
        self.intensity *= (1.0 - self.decay * dt).max(0.0);
    }
    /// Get the current offset to apply to the camera position.
    pub fn offset(&self) -> [f64; 3] {
        let x = self.intensity * self.phase.sin();
        let y = self.intensity * (self.phase * 1.3).cos();
        let z = self.intensity * (self.phase * 0.7).sin();
        [x, y, z]
    }
    /// Check if the shake has effectively stopped.
    pub fn is_done(&self) -> bool {
        self.intensity < 1e-6
    }
    /// Trigger a new shake impulse (adds to existing intensity).
    pub fn impulse(&mut self, amount: f64) {
        self.intensity += amount;
    }
    /// Current intensity value.
    pub fn intensity(&self) -> f64 {
        self.intensity
    }
}
/// Extended depth-of-field with physical lens parameters.
#[derive(Debug, Clone, Copy)]
pub struct PhysicalDoF {
    /// Focal length of the lens in mm.
    pub focal_length_mm: f64,
    /// Lens aperture (f-number, e.g. 1.4, 2.8, 5.6).
    pub f_number: f64,
    /// Focus distance from lens to subject in metres.
    pub focus_distance: f64,
    /// Sensor height in mm (full-frame = 24 mm).
    pub sensor_height_mm: f64,
    /// Rendering resolution height in pixels.
    pub resolution_height: usize,
    /// Circle-of-confusion diameter threshold in mm.
    pub coc_threshold_mm: f64,
}
impl PhysicalDoF {
    /// 35 mm full-frame portrait lens at f/1.8, focused at 2 m.
    pub fn portrait() -> Self {
        Self {
            focal_length_mm: 85.0,
            f_number: 1.8,
            focus_distance: 2.0,
            sensor_height_mm: 24.0,
            resolution_height: 1080,
            coc_threshold_mm: 0.029,
        }
    }
    /// Compute the physical circle-of-confusion diameter (mm) for a point
    /// at `depth` metres from the lens.
    pub fn coc_mm(&self, depth: f64) -> f64 {
        let fl = self.focal_length_mm / 1000.0;
        let fs = self.focus_distance;
        let n = self.f_number;
        if depth.abs() < 1e-12 || fs.abs() < 1e-12 || n.abs() < 1e-12 {
            return 0.0;
        }
        let image_dist_subject = fl * fs / (fs - fl);
        let image_dist_depth = if (depth - fl).abs() > 1e-12 {
            fl * depth / (depth - fl)
        } else {
            1e10
        };
        let diameter = (image_dist_depth - image_dist_subject).abs() / n;
        diameter * 1000.0
    }
    /// Normalised blur factor (0..1) for a given depth.
    pub fn blur_factor(&self, depth: f64) -> f64 {
        let coc = self.coc_mm(depth);
        (coc / self.coc_threshold_mm.max(1e-12)).min(1.0)
    }
    /// Near and far distances that are acceptably sharp (depth-of-field range).
    pub fn dof_range(&self) -> (f64, f64) {
        let fl = self.focal_length_mm / 1000.0;
        let n = self.f_number;
        let c = self.coc_threshold_mm / 1000.0;
        let fs = self.focus_distance;
        let hyperfocal = if c.abs() > 1e-15 {
            fl * fl / (n * c) + fl
        } else {
            1e10
        };
        let near = fs * (hyperfocal - fl) / (hyperfocal + fs - 2.0 * fl);
        let far_denom = hyperfocal - fs;
        let far = if far_denom.abs() < 1e-12 {
            f64::INFINITY
        } else {
            fs * (hyperfocal - fl) / far_denom
        };
        (near.max(0.0), far.max(near))
    }
}
/// View frustum represented as six half-space planes `ax+by+cz+d >= 0`.
pub struct Frustum {
    /// Six planes: `[a, b, c, d]` where `ax+by+cz+d=0`.
    pub planes: [[f64; 4]; 6],
}
impl Frustum {
    /// Extract frustum planes from a camera using the view-projection matrix.
    pub fn from_camera(camera: &Camera) -> Self {
        let m = camera.view_projection();
        let row = |r: usize| [m[0][r], m[1][r], m[2][r], m[3][r]];
        let r0 = row(0);
        let r1 = row(1);
        let r2 = row(2);
        let r3 = row(3);
        let planes = [
            add4(r3, r0),
            sub4(r3, r0),
            add4(r3, r1),
            sub4(r3, r1),
            add4(r3, r2),
            sub4(r3, r2),
        ];
        let planes = planes.map(|p| {
            let len = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
            if len > 1e-12 {
                [p[0] / len, p[1] / len, p[2] / len, p[3] / len]
            } else {
                p
            }
        });
        Frustum { planes }
    }
    /// Returns `true` if the point is inside (or on) all six planes.
    pub fn contains_point(&self, p: [f64; 3]) -> bool {
        for plane in &self.planes {
            if plane[0] * p[0] + plane[1] * p[1] + plane[2] * p[2] + plane[3] < 0.0 {
                return false;
            }
        }
        true
    }
    /// Returns `true` unless the sphere is fully outside at least one plane.
    pub fn contains_sphere(&self, center: [f64; 3], radius: f64) -> bool {
        for plane in &self.planes {
            let dist =
                plane[0] * center[0] + plane[1] * center[1] + plane[2] * center[2] + plane[3];
            if dist < -radius {
                return false;
            }
        }
        true
    }
    /// Returns `true` unless the AABB is fully outside at least one plane.
    pub fn contains_aabb(&self, min: [f64; 3], max: [f64; 3]) -> bool {
        for plane in &self.planes {
            let px = if plane[0] >= 0.0 { max[0] } else { min[0] };
            let py = if plane[1] >= 0.0 { max[1] } else { min[1] };
            let pz = if plane[2] >= 0.0 { max[2] } else { min[2] };
            if plane[0] * px + plane[1] * py + plane[2] * pz + plane[3] < 0.0 {
                return false;
            }
        }
        true
    }
    /// Signed distance from a point to a specific plane.
    pub fn signed_distance(&self, plane_idx: usize, p: [f64; 3]) -> f64 {
        let pl = &self.planes[plane_idx];
        pl[0] * p[0] + pl[1] * p[1] + pl[2] * p[2] + pl[3]
    }
}
