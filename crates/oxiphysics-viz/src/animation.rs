#![allow(clippy::type_complexity)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Animation and keyframe system for physics visualization.
//!
//! Provides interpolation, keyframes, transform animation, particle animation,
//! easing functions, physics replay, and related utilities.

// ─────────────────────────────────────────────────────────────────────────────
// InterpolationType
// ─────────────────────────────────────────────────────────────────────────────

/// Specifies how values between keyframes are interpolated.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum InterpolationType {
    /// Constant interpolation — value jumps at the start of each interval.
    StepStart,
    /// Constant interpolation — value jumps at the end of each interval.
    StepEnd,
    /// Linear interpolation between consecutive keyframes.
    Linear,
    /// Cubic Hermite spline with user-supplied tangents.
    CubicSpline,
    /// Generic ease (smooth cubic S-curve, identical to `EaseInOut`).
    Ease,
    /// Ease-in: slow start, fast end.
    EaseIn,
    /// Ease-out: fast start, slow end.
    EaseOut,
    /// Ease-in-out: slow start, fast end, slow end.
    EaseInOut,
}

// ─────────────────────────────────────────────────────────────────────────────
// Keyframe<T>
// ─────────────────────────────────────────────────────────────────────────────

/// A single keyframe storing a value at a specific time.
///
/// Optional incoming (`tangent_in`) and outgoing (`tangent_out`) tangents are
/// used by cubic-spline interpolation modes.
#[derive(Debug, Clone)]
pub struct Keyframe<T> {
    /// Time of this keyframe (in seconds or arbitrary time units).
    pub time: f64,
    /// Value at this keyframe.
    pub value: T,
    /// Incoming tangent (used for cubic spline interpolation).
    pub tangent_in: Option<T>,
    /// Outgoing tangent (used for cubic spline interpolation).
    pub tangent_out: Option<T>,
}

impl<T: Clone + Default> Keyframe<T> {
    /// Create a new keyframe with no tangent information.
    pub fn new(time: f64, value: T) -> Self {
        Self {
            time,
            value,
            tangent_in: None,
            tangent_out: None,
        }
    }

    /// Create a new keyframe with explicit tangents.
    pub fn with_tangents(time: f64, value: T, tan_in: T, tan_out: T) -> Self {
        Self {
            time,
            value,
            tangent_in: Some(tan_in),
            tangent_out: Some(tan_out),
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// AnimationChannel
// ─────────────────────────────────────────────────────────────────────────────

/// A single-channel animation track holding `f64` keyframes.
///
/// Each channel stores an ordered list of \[`Keyframe`f64`] values and an
/// [`InterpolationType`] that determines how in-between values are computed.
#[derive(Debug, Clone)]
pub struct AnimationChannel {
    /// Ordered list of keyframes (should be sorted ascending by time).
    pub keyframes_f64: Vec<Keyframe<f64>>,
    /// Interpolation method used between consecutive keyframes.
    pub interp: InterpolationType,
}

impl Default for AnimationChannel {
    fn default() -> Self {
        Self {
            keyframes_f64: Vec::new(),
            interp: InterpolationType::Linear,
        }
    }
}

impl AnimationChannel {
    /// Create a new channel with the given interpolation type.
    pub fn new(interp: InterpolationType) -> Self {
        Self {
            keyframes_f64: Vec::new(),
            interp,
        }
    }

    /// Create a channel from a list of `(time, value)` pairs using linear interpolation.
    pub fn from_pairs(pairs: &[(f64, f64)]) -> Self {
        let keyframes_f64 = pairs.iter().map(|&(t, v)| Keyframe::new(t, v)).collect();
        Self {
            keyframes_f64,
            interp: InterpolationType::Linear,
        }
    }

    /// Add a keyframe to the channel (appends; caller is responsible for ordering).
    pub fn push(&mut self, kf: Keyframe<f64>) {
        self.keyframes_f64.push(kf);
    }

    /// Sort keyframes by time in ascending order.
    pub fn sort_by_time(&mut self) {
        self.keyframes_f64.sort_by(|a, b| {
            a.time
                .partial_cmp(&b.time)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
    }

    /// Evaluate the channel at time `t`, returning the interpolated value.
    ///
    /// Returns `0.0` if there are no keyframes, or the nearest boundary value
    /// if `t` is outside the keyframe range.
    pub fn evaluate(&self, t: f64) -> f64 {
        let kfs = &self.keyframes_f64;
        if kfs.is_empty() {
            return 0.0;
        }
        if kfs.len() == 1 {
            return kfs[0].value;
        }
        // Clamp to range
        if t <= kfs.first().expect("collection should not be empty").time {
            return kfs.first().expect("collection should not be empty").value;
        }
        if t >= kfs.last().expect("collection should not be empty").time {
            return kfs.last().expect("collection should not be empty").value;
        }

        // Find the segment [i, i+1] that contains t
        let idx = kfs.partition_point(|k| k.time <= t).saturating_sub(1);
        let idx = idx.min(kfs.len() - 2);

        let k0 = &kfs[idx];
        let k1 = &kfs[idx + 1];
        let dt = k1.time - k0.time;
        if dt <= 0.0 {
            return k0.value;
        }
        let local_t = (t - k0.time) / dt;

        match self.interp {
            InterpolationType::StepStart => k0.value,
            InterpolationType::StepEnd => k1.value,
            InterpolationType::Linear => k0.value + (k1.value - k0.value) * local_t,
            InterpolationType::CubicSpline => {
                // Use tangents if available, otherwise Catmull-Rom
                let m0 = k0
                    .tangent_out
                    .unwrap_or_else(|| self.catmull_rom_tangent_at(idx));
                let m1 = k1
                    .tangent_in
                    .unwrap_or_else(|| self.catmull_rom_tangent_at(idx + 1));
                Self::cubic_hermite(local_t, k0.value, k1.value, m0 * dt, m1 * dt)
            }
            InterpolationType::Ease | InterpolationType::EaseInOut => {
                let s = EasingFunction::ease_in_out_cubic(local_t);
                k0.value + (k1.value - k0.value) * s
            }
            InterpolationType::EaseIn => {
                let s = EasingFunction::ease_in_quad(local_t);
                k0.value + (k1.value - k0.value) * s
            }
            InterpolationType::EaseOut => {
                let s = EasingFunction::ease_out_quad(local_t);
                k0.value + (k1.value - k0.value) * s
            }
        }
    }

    /// Cubic Hermite interpolation.
    ///
    /// Interpolates between `p0` and `p1` at normalized parameter `t` ∈ [0, 1]
    /// with derivative (tangent) values `m0` at `p0` and `m1` at `p1`.
    pub fn cubic_hermite(t: f64, p0: f64, p1: f64, m0: f64, m1: f64) -> f64 {
        let t2 = t * t;
        let t3 = t2 * t;
        let h00 = 2.0 * t3 - 3.0 * t2 + 1.0;
        let h10 = t3 - 2.0 * t2 + t;
        let h01 = -2.0 * t3 + 3.0 * t2;
        let h11 = t3 - t2;
        h00 * p0 + h10 * m0 + h01 * p1 + h11 * m1
    }

    /// Compute the Catmull-Rom tangent at keyframe index `i`.
    fn catmull_rom_tangent_at(&self, i: usize) -> f64 {
        let kfs = &self.keyframes_f64;
        if kfs.len() < 2 {
            return 0.0;
        }
        if i == 0 {
            // Forward difference
            return (kfs[1].value - kfs[0].value) / (kfs[1].time - kfs[0].time).max(1e-12);
        }
        if i >= kfs.len() - 1 {
            let n = kfs.len() - 1;
            return (kfs[n].value - kfs[n - 1].value) / (kfs[n].time - kfs[n - 1].time).max(1e-12);
        }
        Self::catmull_rom_tangent(kfs[i - 1].value, kfs[i + 1].value)
    }

    /// Compute the Catmull-Rom tangent between two surrounding values.
    ///
    /// This produces a centripetal tangent estimate of `(p_next - p_prev) / 2`.
    pub fn catmull_rom_tangent(p_prev: f64, p_next: f64) -> f64 {
        (p_next - p_prev) * 0.5
    }

    /// Return the total duration of the channel (last keyframe time − first keyframe time).
    pub fn duration(&self) -> f64 {
        if self.keyframes_f64.len() < 2 {
            return 0.0;
        }
        self.keyframes_f64
            .last()
            .expect("collection should not be empty")
            .time
            - self
                .keyframes_f64
                .first()
                .expect("collection should not be empty")
                .time
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Transform3DAnimation
// ─────────────────────────────────────────────────────────────────────────────

/// Animated 3-D transform built from nine independent [`AnimationChannel`]s.
///
/// The transform is decomposed into translation (tx, ty, tz), Euler rotation
/// (rx, ry, rz) in radians, and scale (sx, sy, sz).
#[derive(Debug, Clone)]
pub struct Transform3DAnimation {
    /// Translation X channel.
    pub tx: AnimationChannel,
    /// Translation Y channel.
    pub ty: AnimationChannel,
    /// Translation Z channel.
    pub tz: AnimationChannel,
    /// Rotation X (Euler, radians) channel.
    pub rx: AnimationChannel,
    /// Rotation Y (Euler, radians) channel.
    pub ry: AnimationChannel,
    /// Rotation Z (Euler, radians) channel.
    pub rz: AnimationChannel,
    /// Scale X channel.
    pub sx: AnimationChannel,
    /// Scale Y channel.
    pub sy: AnimationChannel,
    /// Scale Z channel.
    pub sz: AnimationChannel,
}

impl Default for Transform3DAnimation {
    fn default() -> Self {
        let mut sx = AnimationChannel::default();
        sx.push(Keyframe::new(0.0, 1.0));
        let mut sy = AnimationChannel::default();
        sy.push(Keyframe::new(0.0, 1.0));
        let mut sz = AnimationChannel::default();
        sz.push(Keyframe::new(0.0, 1.0));
        Self {
            tx: AnimationChannel::default(),
            ty: AnimationChannel::default(),
            tz: AnimationChannel::default(),
            rx: AnimationChannel::default(),
            ry: AnimationChannel::default(),
            rz: AnimationChannel::default(),
            sx,
            sy,
            sz,
        }
    }
}

impl Transform3DAnimation {
    /// Evaluate the transform at time `t`.
    ///
    /// Returns `(translation, rotation_euler, scale)` as `\[f64; 3\]` triplets.
    pub fn evaluate(&self, t: f64) -> ([f64; 3], [f64; 3], [f64; 3]) {
        let translation = [
            self.tx.evaluate(t),
            self.ty.evaluate(t),
            self.tz.evaluate(t),
        ];
        let rotation = [
            self.rx.evaluate(t),
            self.ry.evaluate(t),
            self.rz.evaluate(t),
        ];
        let scale = [
            self.sx.evaluate(t),
            self.sy.evaluate(t),
            self.sz.evaluate(t),
        ];
        (translation, rotation, scale)
    }

    /// Set translation keyframes from slices of `(time, value)` pairs.
    pub fn set_translation_keys(&mut self, x: &[(f64, f64)], y: &[(f64, f64)], z: &[(f64, f64)]) {
        self.tx = AnimationChannel::from_pairs(x);
        self.ty = AnimationChannel::from_pairs(y);
        self.tz = AnimationChannel::from_pairs(z);
    }

    /// Set rotation keyframes from slices of `(time, value)` pairs.
    pub fn set_rotation_keys(&mut self, x: &[(f64, f64)], y: &[(f64, f64)], z: &[(f64, f64)]) {
        self.rx = AnimationChannel::from_pairs(x);
        self.ry = AnimationChannel::from_pairs(y);
        self.rz = AnimationChannel::from_pairs(z);
    }

    /// Set scale keyframes from slices of `(time, value)` pairs.
    pub fn set_scale_keys(&mut self, x: &[(f64, f64)], y: &[(f64, f64)], z: &[(f64, f64)]) {
        self.sx = AnimationChannel::from_pairs(x);
        self.sy = AnimationChannel::from_pairs(y);
        self.sz = AnimationChannel::from_pairs(z);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ParticleAnimation
// ─────────────────────────────────────────────────────────────────────────────

/// Stores a sequence of particle position frames for replay and interpolation.
///
/// Each frame contains a snapshot of all particle positions at a specific time.
#[derive(Debug, Clone, Default)]
pub struct ParticleAnimation {
    /// Per-frame particle positions: `positions\[frame\][particle] = \[x, y, z\]`.
    pub positions: Vec<Vec<[f64; 3]>>,
    /// Time stamp for each frame.
    pub times: Vec<f64>,
}

impl ParticleAnimation {
    /// Create an empty `ParticleAnimation`.
    pub fn new() -> Self {
        Self::default()
    }

    /// Append a new frame at time `t` with the given particle positions.
    pub fn push_frame(&mut self, t: f64, positions: Vec<[f64; 3]>) {
        self.times.push(t);
        self.positions.push(positions);
    }

    /// Linearly interpolate particle positions at time `t`.
    ///
    /// Returns a `Vec<\[f64;3\]>` with one entry per particle.  If the animation
    /// has no frames an empty vector is returned; if `t` is outside the recorded
    /// range the nearest boundary frame is returned.
    pub fn interpolate_frame(&self, t: f64) -> Vec<[f64; 3]> {
        if self.positions.is_empty() {
            return Vec::new();
        }
        let n = self.times.len();
        if t <= self.times[0] {
            return self.positions[0].clone();
        }
        if t >= self.times[n - 1] {
            return self.positions[n - 1].clone();
        }
        let idx = self.times.partition_point(|&s| s <= t).saturating_sub(1);
        let idx = idx.min(n - 2);
        let t0 = self.times[idx];
        let t1 = self.times[idx + 1];
        let dt = (t1 - t0).max(1e-15);
        let alpha = (t - t0) / dt;

        let f0 = &self.positions[idx];
        let f1 = &self.positions[idx + 1];
        let len = f0.len().min(f1.len());
        (0..len)
            .map(|i| {
                [
                    f0[i][0] + (f1[i][0] - f0[i][0]) * alpha,
                    f0[i][1] + (f1[i][1] - f0[i][1]) * alpha,
                    f0[i][2] + (f1[i][2] - f0[i][2]) * alpha,
                ]
            })
            .collect()
    }

    /// Total duration of the animation (last time − first time).
    pub fn total_duration(&self) -> f64 {
        if self.times.len() < 2 {
            return 0.0;
        }
        self.times.last().copied().unwrap_or(0.0) - self.times.first().copied().unwrap_or(0.0)
    }

    /// Number of particles (taken from the first frame, or 0 if empty).
    pub fn n_particles(&self) -> usize {
        self.positions.first().map_or(0, |f| f.len())
    }

    /// Number of recorded frames.
    pub fn n_frames(&self) -> usize {
        self.positions.len()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// AnimationClip
// ─────────────────────────────────────────────────────────────────────────────

/// A named collection of [`AnimationChannel`]s with a fixed duration.
///
/// Supports looping and ping-pong playback modes.
#[derive(Debug, Clone)]
pub struct AnimationClip {
    /// Human-readable name for this clip.
    pub name: String,
    /// List of animation channels contained in this clip.
    pub channels: Vec<AnimationChannel>,
    /// Total playback duration in seconds.
    pub duration: f64,
}

impl AnimationClip {
    /// Create a new, empty animation clip with the given name and duration.
    pub fn new(name: impl Into<String>, duration: f64) -> Self {
        Self {
            name: name.into(),
            channels: Vec::new(),
            duration,
        }
    }

    /// Add a channel to this clip.
    pub fn add_channel(&mut self, channel: AnimationChannel) {
        self.channels.push(channel);
    }

    /// Wrap `t` into `\[0, duration)` — standard loop playback.
    ///
    /// Returns 0 if `duration` is non-positive.
    pub fn loop_animation(&self, t: f64) -> f64 {
        if self.duration <= 0.0 {
            return 0.0;
        }
        t.rem_euclid(self.duration)
    }

    /// Ping-pong the time value so the animation plays forward then backward.
    ///
    /// The effective period is `2 * duration`.  Returns 0 if `duration` is
    /// non-positive.
    pub fn ping_pong(&self, t: f64) -> f64 {
        if self.duration <= 0.0 {
            return 0.0;
        }
        let period = 2.0 * self.duration;
        let wrapped = t.rem_euclid(period);
        if wrapped < self.duration {
            wrapped
        } else {
            period - wrapped
        }
    }

    /// Evaluate all channels at time `t`, returning one `f64` per channel.
    pub fn evaluate_all(&self, t: f64) -> Vec<f64> {
        self.channels.iter().map(|ch| ch.evaluate(t)).collect()
    }

    /// Return the number of channels in this clip.
    pub fn channel_count(&self) -> usize {
        self.channels.len()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// EasingFunction
// ─────────────────────────────────────────────────────────────────────────────

/// A collection of pure easing functions operating on a normalized parameter `t` ∈ [0, 1].
///
/// All functions accept a value in `\[0, 1\]` and return a remapped value suitable
/// for use as a blend factor in animation interpolation.
pub struct EasingFunction;

impl EasingFunction {
    /// Linear ease — identity mapping.
    #[inline]
    pub fn ease_linear(t: f64) -> f64 {
        t
    }

    /// Quadratic ease-in — slow start.
    #[inline]
    pub fn ease_in_quad(t: f64) -> f64 {
        t * t
    }

    /// Quadratic ease-out — slow end.
    #[inline]
    pub fn ease_out_quad(t: f64) -> f64 {
        t * (2.0 - t)
    }

    /// Quadratic ease-in-out.
    #[inline]
    pub fn ease_in_out_quad(t: f64) -> f64 {
        if t < 0.5 {
            2.0 * t * t
        } else {
            -1.0 + (4.0 - 2.0 * t) * t
        }
    }

    /// Cubic ease-in.
    #[inline]
    pub fn ease_in_cubic(t: f64) -> f64 {
        t * t * t
    }

    /// Cubic ease-out.
    #[inline]
    pub fn ease_out_cubic(t: f64) -> f64 {
        let t1 = t - 1.0;
        t1 * t1 * t1 + 1.0
    }

    /// Cubic ease-in-out (smooth S-curve).
    #[inline]
    pub fn ease_in_out_cubic(t: f64) -> f64 {
        if t < 0.5 {
            4.0 * t * t * t
        } else {
            let t1 = 2.0 * t - 2.0;
            0.5 * t1 * t1 * t1 + 1.0
        }
    }

    /// Elastic ease-in.
    pub fn ease_in_elastic(t: f64) -> f64 {
        if t == 0.0 {
            return 0.0;
        }
        if (t - 1.0).abs() < 1e-15 {
            return 1.0;
        }
        let c4 = std::f64::consts::TAU / 3.0;
        -(2.0_f64.powf(10.0 * t - 10.0)) * ((10.0 * t - 10.75) * c4).sin()
    }

    /// Elastic ease-out.
    pub fn ease_out_elastic(t: f64) -> f64 {
        if t == 0.0 {
            return 0.0;
        }
        if (t - 1.0).abs() < 1e-15 {
            return 1.0;
        }
        let c4 = std::f64::consts::TAU / 3.0;
        2.0_f64.powf(-10.0 * t) * ((10.0 * t - 0.75) * c4).sin() + 1.0
    }

    /// Bounce ease-out.
    pub fn ease_out_bounce(t: f64) -> f64 {
        const N: f64 = 7.5625;
        const D: f64 = 2.75;
        let mut x = t;
        if x < 1.0 / D {
            N * x * x
        } else if x < 2.0 / D {
            x -= 1.5 / D;
            N * x * x + 0.75
        } else if x < 2.5 / D {
            x -= 2.25 / D;
            N * x * x + 0.9375
        } else {
            x -= 2.625 / D;
            N * x * x + 0.984375
        }
    }

    /// Ease-in-out with back overshoot.
    pub fn ease_in_out_back(t: f64) -> f64 {
        let c1 = 1.70158_f64;
        let c2 = c1 * 1.525;
        if t < 0.5 {
            ((2.0 * t).powi(2) * ((c2 + 1.0) * 2.0 * t - c2)) / 2.0
        } else {
            ((2.0 * t - 2.0).powi(2) * ((c2 + 1.0) * (2.0 * t - 2.0) + c2) + 2.0) / 2.0
        }
    }

    /// Damped spring oscillation.
    ///
    /// Returns a value that approaches 1.0 with oscillations characterized by
    /// `stiffness` (oscillation frequency) and `damping` (decay rate).
    pub fn spring_damped(t: f64, stiffness: f64, damping: f64) -> f64 {
        if t <= 0.0 {
            return 0.0;
        }
        let omega = stiffness.sqrt().max(1e-12);
        let decay = (-damping * t).exp();
        1.0 - decay * (omega * t).cos()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PhysicsReplay
// ─────────────────────────────────────────────────────────────────────────────

/// Records and replays a physics simulation trajectory.
///
/// Each step stores position and velocity snapshots with an associated
/// simulation time.  Replay can retrieve the state at an arbitrary time via
/// linear interpolation.
#[derive(Debug, Clone, Default)]
pub struct PhysicsReplay {
    /// Recorded simulation times.
    pub times: Vec<f64>,
    /// Particle positions at each recorded time step.
    pub positions: Vec<Vec<[f64; 3]>>,
    /// Particle velocities at each recorded time step.
    pub velocities: Vec<Vec<[f64; 3]>>,
}

impl PhysicsReplay {
    /// Create a new, empty `PhysicsReplay`.
    pub fn new() -> Self {
        Self::default()
    }

    /// Record a simulation step.
    ///
    /// `positions` and `velocities` must have the same length (one entry per
    /// simulated particle).
    pub fn record_step(&mut self, positions: Vec<[f64; 3]>, velocities: Vec<[f64; 3]>, t: f64) {
        self.times.push(t);
        self.positions.push(positions);
        self.velocities.push(velocities);
    }

    /// Replay the simulation state at time `t`.
    ///
    /// Linearly interpolates between the two surrounding recorded steps.
    /// Returns `None` if no steps have been recorded.
    pub fn replay_at(&self, t: f64) -> Option<(Vec<[f64; 3]>, Vec<[f64; 3]>)> {
        let n = self.times.len();
        if n == 0 {
            return None;
        }
        if n == 1 || t <= self.times[0] {
            return Some((self.positions[0].clone(), self.velocities[0].clone()));
        }
        if t >= self.times[n - 1] {
            return Some((
                self.positions[n - 1].clone(),
                self.velocities[n - 1].clone(),
            ));
        }
        let idx = self.times.partition_point(|&s| s <= t).saturating_sub(1);
        let idx = idx.min(n - 2);
        let t0 = self.times[idx];
        let t1 = self.times[idx + 1];
        let alpha = if (t1 - t0).abs() < 1e-15 {
            0.0
        } else {
            (t - t0) / (t1 - t0)
        };

        let interp3 = |a: [f64; 3], b: [f64; 3]| -> [f64; 3] {
            [
                a[0] + (b[0] - a[0]) * alpha,
                a[1] + (b[1] - a[1]) * alpha,
                a[2] + (b[2] - a[2]) * alpha,
            ]
        };

        let p0 = &self.positions[idx];
        let p1 = &self.positions[idx + 1];
        let v0 = &self.velocities[idx];
        let v1 = &self.velocities[idx + 1];
        let np = p0.len().min(p1.len());
        let nv = v0.len().min(v1.len());

        Some((
            (0..np).map(|i| interp3(p0[i], p1[i])).collect(),
            (0..nv).map(|i| interp3(v0[i], v1[i])).collect(),
        ))
    }

    /// Export particle positions sampled at a uniform frame rate.
    ///
    /// `fps` specifies frames per second.  Returns one `Vec<\[f64;3\]>` per
    /// output frame.  Returns an empty vector if there are no recorded steps.
    pub fn export_frames(&self, fps: f64) -> Vec<Vec<[f64; 3]>> {
        if self.times.is_empty() || fps <= 0.0 {
            return Vec::new();
        }
        let t_start = self.times[0];
        let t_end = *self.times.last().expect("collection should not be empty");
        let dt = 1.0 / fps;
        let n_frames = ((t_end - t_start) * fps).ceil() as usize + 1;
        let mut out = Vec::with_capacity(n_frames);
        let mut t = t_start;
        while t <= t_end + 1e-9 {
            if let Some((pos, _vel)) = self.replay_at(t) {
                out.push(pos);
            }
            t += dt;
        }
        out
    }

    /// Total recorded duration.
    pub fn total_duration(&self) -> f64 {
        if self.times.len() < 2 {
            return 0.0;
        }
        self.times.last().copied().unwrap_or(0.0) - self.times.first().copied().unwrap_or(0.0)
    }

    /// Number of recorded steps.
    pub fn n_steps(&self) -> usize {
        self.times.len()
    }

    /// Number of particles (from the first recorded step, or 0).
    pub fn n_particles(&self) -> usize {
        self.positions.first().map_or(0, |f| f.len())
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Required free-function API
// ─────────────────────────────────────────────────────────────────────────────

/// A simple keyframe with a time and scalar value.
#[derive(Debug, Clone, PartialEq)]
pub struct SimpleKeyframe {
    /// Time of this keyframe.
    pub time: f64,
    /// Value at this keyframe.
    pub value: f64,
}

impl SimpleKeyframe {
    /// Create a new keyframe.
    pub fn new(time: f64, value: f64) -> Self {
        Self { time, value }
    }
}

/// An animation curve: a sequence of keyframes with tangents.
#[derive(Debug, Clone, Default)]
pub struct AnimationCurve {
    /// Keyframes defining the curve.
    pub keyframes: Vec<SimpleKeyframe>,
    /// Per-keyframe tangents (same length as `keyframes`).
    pub tangents: Vec<f64>,
}

impl AnimationCurve {
    /// Create a new empty curve.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a keyframe with a tangent.
    pub fn push(&mut self, kf: SimpleKeyframe, tangent: f64) {
        self.keyframes.push(kf);
        self.tangents.push(tangent);
    }
}

/// Linearly interpolate a value from a slice of [`SimpleKeyframe`]s at time `t`.
///
/// Returns the value of the first keyframe if `t` is before all keyframes, and
/// the value of the last keyframe if `t` is after all keyframes.
pub fn lerp_keyframes(keyframes: &[SimpleKeyframe], t: f64) -> f64 {
    if keyframes.is_empty() {
        return 0.0;
    }
    if t <= keyframes[0].time {
        return keyframes[0].value;
    }
    if t >= keyframes[keyframes.len() - 1].time {
        return keyframes[keyframes.len() - 1].value;
    }
    // Find the surrounding keyframes
    for i in 0..keyframes.len() - 1 {
        let k0 = &keyframes[i];
        let k1 = &keyframes[i + 1];
        if t >= k0.time && t <= k1.time {
            let dt = k1.time - k0.time;
            if dt < 1e-15 {
                return k0.value;
            }
            let alpha = (t - k0.time) / dt;
            return k0.value + alpha * (k1.value - k0.value);
        }
    }
    keyframes[keyframes.len() - 1].value
}

/// Cubic Hermite spline interpolation.
///
/// Interpolates between `p0` and `p1` with outgoing tangent `m0` and
/// incoming tangent `m1` at parameter `t` ∈ [0, 1].
pub fn cubic_hermite_interp(p0: f64, p1: f64, m0: f64, m1: f64, t: f64) -> f64 {
    let t2 = t * t;
    let t3 = t2 * t;
    let h00 = 2.0 * t3 - 3.0 * t2 + 1.0;
    let h10 = t3 - 2.0 * t2 + t;
    let h01 = -2.0 * t3 + 3.0 * t2;
    let h11 = t3 - t2;
    h00 * p0 + h10 * m0 + h01 * p1 + h11 * m1
}

/// Cubic Bezier interpolation in 1D.
///
/// `p0`..`p3` are the four control points; `t` ∈ [0, 1].
pub fn bezier_cubic_1d(p0: f64, p1: f64, p2: f64, p3: f64, t: f64) -> f64 {
    let u = 1.0 - t;
    u * u * u * p0 + 3.0 * u * u * t * p1 + 3.0 * u * t * t * p2 + t * t * t * p3
}

/// Evaluate an [`AnimationCurve`] at time `t` using Hermite interpolation.
///
/// Returns 0 if the curve has no keyframes.
pub fn eval_curve(curve: &AnimationCurve, t: f64) -> f64 {
    let kfs = &curve.keyframes;
    let tans = &curve.tangents;
    if kfs.is_empty() {
        return 0.0;
    }
    if t <= kfs[0].time {
        return kfs[0].value;
    }
    if t >= kfs[kfs.len() - 1].time {
        return kfs[kfs.len() - 1].value;
    }
    for i in 0..kfs.len() - 1 {
        let k0 = &kfs[i];
        let k1 = &kfs[i + 1];
        if t >= k0.time && t <= k1.time {
            let dt = k1.time - k0.time;
            if dt < 1e-15 {
                return k0.value;
            }
            let alpha = (t - k0.time) / dt;
            let m0 = if i < tans.len() { tans[i] * dt } else { 0.0 };
            let m1 = if i + 1 < tans.len() {
                tans[i + 1] * dt
            } else {
                0.0
            };
            return cubic_hermite_interp(k0.value, k1.value, m0, m1, alpha);
        }
    }
    kfs[kfs.len() - 1].value
}

/// Ease-in-out quadratic easing function.
///
/// Returns a smoothed value for `t` ∈ [0, 1].
pub fn ease_in_out_quad(t: f64) -> f64 {
    EasingFunction::ease_in_out_quad(t)
}

/// Ease-in-out cubic easing function.
///
/// Returns a smoothed value for `t` ∈ [0, 1].
pub fn ease_in_out_cubic(t: f64) -> f64 {
    EasingFunction::ease_in_out_cubic(t)
}

/// Ease-in exponential easing function.
///
/// Returns a smoothed value for `t` ∈ [0, 1].
pub fn ease_in_expo(t: f64) -> f64 {
    if t <= 0.0 {
        0.0
    } else {
        (2.0_f64).powf(10.0 * t - 10.0)
    }
}

/// Ease-out bounce easing function.
///
/// Returns a value in [0, 1] for `t` ∈ [0, 1].
pub fn ease_out_bounce(t: f64) -> f64 {
    EasingFunction::ease_out_bounce(t)
}

/// Spring animation step: returns `(new_position, new_velocity)`.
///
/// Models a damped spring with `stiffness` and `damping` using an Euler step of `dt`.
pub fn spring_animation(
    pos: f64,
    vel: f64,
    target: f64,
    stiffness: f64,
    damping: f64,
    dt: f64,
) -> (f64, f64) {
    let force = stiffness * (target - pos) - damping * vel;
    let new_vel = vel + force * dt;
    let new_pos = pos + new_vel * dt;
    (new_pos, new_vel)
}

// ─────────────────────────────────────────────────────────────────────────────
// Two-bone IK
// ─────────────────────────────────────────────────────────────────────────────

#[inline]
fn vec3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn vec3_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn vec3_scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline]
fn vec3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn vec3_len(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

#[inline]
fn vec3_norm(a: [f64; 3]) -> [f64; 3] {
    let l = vec3_len(a);
    if l < 1e-15 {
        [0.0, 0.0, 0.0]
    } else {
        vec3_scale(a, 1.0 / l)
    }
}

#[inline]
fn vec3_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Two-bone IK solver: returns `(new_mid_joint, new_end_effector)`.
///
/// Given bone lengths `l1` (root→mid) and `l2` (mid→end), moves the mid joint
/// to reach `target` from `root`.  The `mid` and `end` arguments provide the
/// current configuration (used to determine the bend plane).
#[allow(clippy::too_many_arguments)]
pub fn ik_two_bone(
    root: [f64; 3],
    mid: [f64; 3],
    _end: [f64; 3],
    target: [f64; 3],
    l1: f64,
    l2: f64,
) -> ([f64; 3], [f64; 3]) {
    let to_target = vec3_sub(target, root);
    let dist = vec3_len(to_target).clamp(1e-15, l1 + l2);

    // Clamp target to reachable distance
    let target_clamped = if dist >= l1 + l2 {
        // Fully extended: mid and end along root→target
        let dir = vec3_norm(to_target);
        let new_mid = vec3_add(root, vec3_scale(dir, l1));
        let new_end = vec3_add(root, vec3_scale(dir, l1 + l2));
        return (new_mid, new_end);
    } else {
        target
    };

    // Law of cosines to find angle at root
    let d = vec3_len(vec3_sub(target_clamped, root)).max(1e-15);
    let cos_angle_root = ((d * d + l1 * l1 - l2 * l2) / (2.0 * d * l1)).clamp(-1.0, 1.0);
    let angle_root = cos_angle_root.acos();

    // Bend plane: use root→target and root→mid to form normal
    let root_to_target = vec3_norm(vec3_sub(target_clamped, root));
    let root_to_mid = vec3_norm(vec3_sub(mid, root));
    let bend_normal = {
        let n = vec3_cross(root_to_target, root_to_mid);
        let nl = vec3_len(n);
        if nl < 1e-8 {
            // Degenerate: pick a perpendicular
            let perp = if root_to_target[0].abs() < 0.9 {
                vec3_cross(root_to_target, [1.0, 0.0, 0.0])
            } else {
                vec3_cross(root_to_target, [0.0, 1.0, 0.0])
            };
            vec3_norm(perp)
        } else {
            vec3_scale(n, 1.0 / nl)
        }
    };

    // Rotate root_to_target by angle_root around bend_normal (Rodrigues)
    let k = bend_normal;
    let v = root_to_target;
    let cos_a = angle_root.cos();
    let sin_a = angle_root.sin();
    let mid_dir = [
        v[0] * cos_a + (k[1] * v[2] - k[2] * v[1]) * sin_a + k[0] * vec3_dot(k, v) * (1.0 - cos_a),
        v[1] * cos_a + (k[2] * v[0] - k[0] * v[2]) * sin_a + k[1] * vec3_dot(k, v) * (1.0 - cos_a),
        v[2] * cos_a + (k[0] * v[1] - k[1] * v[0]) * sin_a + k[2] * vec3_dot(k, v) * (1.0 - cos_a),
    ];

    let new_mid = vec3_add(root, vec3_scale(mid_dir, l1));
    let new_end = target_clamped;

    (new_mid, new_end)
}

// ─────────────────────────────────────────────────────────────────────────────
// Helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Linearly interpolate between two 3-D points.
///
/// `alpha` = 0 returns `a`; `alpha` = 1 returns `b`.
pub fn lerp3(a: [f64; 3], b: [f64; 3], alpha: f64) -> [f64; 3] {
    [
        a[0] + (b[0] - a[0]) * alpha,
        a[1] + (b[1] - a[1]) * alpha,
        a[2] + (b[2] - a[2]) * alpha,
    ]
}

/// Compute the Euclidean distance between two 3-D points.
pub fn distance3(a: [f64; 3], b: [f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    // ── AnimationChannel ──────────────────────────────────────────────────────

    #[test]
    fn test_channel_empty_returns_zero() {
        let ch = AnimationChannel::default();
        assert_eq!(ch.evaluate(0.5), 0.0);
    }

    #[test]
    fn test_channel_single_keyframe() {
        let mut ch = AnimationChannel::default();
        ch.push(Keyframe::new(0.0, 42.0));
        assert_eq!(ch.evaluate(0.0), 42.0);
        assert_eq!(ch.evaluate(999.0), 42.0);
    }

    #[test]
    fn test_channel_linear_midpoint() {
        let ch = AnimationChannel::from_pairs(&[(0.0, 0.0), (1.0, 10.0)]);
        let v = ch.evaluate(0.5);
        assert!((v - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_channel_linear_clamp_before() {
        let ch = AnimationChannel::from_pairs(&[(1.0, 5.0), (2.0, 10.0)]);
        assert_eq!(ch.evaluate(0.0), 5.0);
    }

    #[test]
    fn test_channel_linear_clamp_after() {
        let ch = AnimationChannel::from_pairs(&[(0.0, 0.0), (1.0, 10.0)]);
        assert_eq!(ch.evaluate(2.0), 10.0);
    }

    #[test]
    fn test_channel_step_start() {
        let mut ch = AnimationChannel::new(InterpolationType::StepStart);
        ch.push(Keyframe::new(0.0, 0.0));
        ch.push(Keyframe::new(1.0, 100.0));
        // In interval [0,1), step_start returns k0.value
        assert_eq!(ch.evaluate(0.5), 0.0);
    }

    #[test]
    fn test_channel_step_end() {
        let mut ch = AnimationChannel::new(InterpolationType::StepEnd);
        ch.push(Keyframe::new(0.0, 0.0));
        ch.push(Keyframe::new(1.0, 100.0));
        assert_eq!(ch.evaluate(0.5), 100.0);
    }

    #[test]
    fn test_channel_cubic_spline_endpoints() {
        let mut ch = AnimationChannel::new(InterpolationType::CubicSpline);
        ch.push(Keyframe::new(0.0, 0.0));
        ch.push(Keyframe::new(1.0, 1.0));
        assert!((ch.evaluate(0.0) - 0.0).abs() < 1e-12);
        assert!((ch.evaluate(1.0) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_channel_ease_in_out_monotone() {
        let mut ch = AnimationChannel::new(InterpolationType::EaseInOut);
        ch.push(Keyframe::new(0.0, 0.0));
        ch.push(Keyframe::new(1.0, 1.0));
        let mut prev = ch.evaluate(0.0);
        for i in 1..=20 {
            let t = i as f64 / 20.0;
            let v = ch.evaluate(t);
            assert!(v >= prev - 1e-12, "not monotone at t={t}: {v} < {prev}");
            prev = v;
        }
    }

    #[test]
    fn test_channel_ease_in_monotone() {
        let mut ch = AnimationChannel::new(InterpolationType::EaseIn);
        ch.push(Keyframe::new(0.0, 0.0));
        ch.push(Keyframe::new(1.0, 1.0));
        let mut prev = ch.evaluate(0.0);
        for i in 1..=20 {
            let t = i as f64 / 20.0;
            let v = ch.evaluate(t);
            assert!(v >= prev - 1e-12);
            prev = v;
        }
    }

    #[test]
    fn test_channel_ease_out_monotone() {
        let mut ch = AnimationChannel::new(InterpolationType::EaseOut);
        ch.push(Keyframe::new(0.0, 0.0));
        ch.push(Keyframe::new(1.0, 1.0));
        let mut prev = ch.evaluate(0.0);
        for i in 1..=20 {
            let t = i as f64 / 20.0;
            let v = ch.evaluate(t);
            assert!(v >= prev - 1e-12);
            prev = v;
        }
    }

    #[test]
    fn test_channel_duration() {
        let ch = AnimationChannel::from_pairs(&[(2.0, 0.0), (5.0, 1.0)]);
        assert!((ch.duration() - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_channel_sort_by_time() {
        let mut ch = AnimationChannel::from_pairs(&[(2.0, 1.0), (0.0, 0.0), (1.0, 0.5)]);
        ch.sort_by_time();
        let times: Vec<f64> = ch.keyframes_f64.iter().map(|k| k.time).collect();
        assert_eq!(times, vec![0.0, 1.0, 2.0]);
    }

    // ── cubic_hermite ─────────────────────────────────────────────────────────

    #[test]
    fn test_cubic_hermite_endpoints() {
        let v = AnimationChannel::cubic_hermite(0.0, 1.0, 5.0, 0.0, 0.0);
        assert!((v - 1.0).abs() < 1e-12);
        let v = AnimationChannel::cubic_hermite(1.0, 1.0, 5.0, 0.0, 0.0);
        assert!((v - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_cubic_hermite_midpoint_zero_tangents() {
        // With zero tangents the midpoint should be between p0 and p1
        let v = AnimationChannel::cubic_hermite(0.5, 0.0, 4.0, 0.0, 0.0);
        assert!(v > 0.0 && v < 4.0);
    }

    // ── catmull_rom_tangent ────────────────────────────────────────────────────

    #[test]
    fn test_catmull_rom_tangent_symmetric() {
        let t = AnimationChannel::catmull_rom_tangent(0.0, 2.0);
        assert!((t - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_catmull_rom_tangent_flat() {
        let t = AnimationChannel::catmull_rom_tangent(3.0, 3.0);
        assert!(t.abs() < 1e-12);
    }

    // ── Transform3DAnimation ──────────────────────────────────────────────────

    #[test]
    fn test_transform_default_scale_one() {
        let anim = Transform3DAnimation::default();
        let (_tr, _rot, scale) = anim.evaluate(0.0);
        assert!((scale[0] - 1.0).abs() < 1e-12);
        assert!((scale[1] - 1.0).abs() < 1e-12);
        assert!((scale[2] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_transform_translation_keys() {
        let mut anim = Transform3DAnimation::default();
        anim.set_translation_keys(
            &[(0.0, 0.0), (1.0, 10.0)],
            &[(0.0, 0.0), (1.0, 0.0)],
            &[(0.0, 0.0), (1.0, 0.0)],
        );
        let (tr, _, _) = anim.evaluate(0.5);
        assert!((tr[0] - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_transform_rotation_keys() {
        let mut anim = Transform3DAnimation::default();
        anim.set_rotation_keys(&[(0.0, 0.0), (1.0, PI)], &[(0.0, 0.0)], &[(0.0, 0.0)]);
        let (_, rot, _) = anim.evaluate(0.5);
        assert!((rot[0] - PI / 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_transform_scale_keys() {
        let mut anim = Transform3DAnimation::default();
        anim.set_scale_keys(
            &[(0.0, 1.0), (1.0, 2.0)],
            &[(0.0, 1.0), (1.0, 3.0)],
            &[(0.0, 1.0), (1.0, 4.0)],
        );
        let (_, _, scale) = anim.evaluate(1.0);
        assert!((scale[0] - 2.0).abs() < 1e-12);
        assert!((scale[1] - 3.0).abs() < 1e-12);
        assert!((scale[2] - 4.0).abs() < 1e-12);
    }

    // ── ParticleAnimation ─────────────────────────────────────────────────────

    #[test]
    fn test_particle_animation_empty() {
        let anim = ParticleAnimation::new();
        assert_eq!(anim.n_frames(), 0);
        assert_eq!(anim.n_particles(), 0);
        assert!(anim.interpolate_frame(0.0).is_empty());
    }

    #[test]
    fn test_particle_animation_push_and_query() {
        let mut anim = ParticleAnimation::new();
        anim.push_frame(0.0, vec![[0.0, 0.0, 0.0]]);
        anim.push_frame(1.0, vec![[2.0, 0.0, 0.0]]);
        let frame = anim.interpolate_frame(0.5);
        assert!((frame[0][0] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_particle_animation_clamp_before() {
        let mut anim = ParticleAnimation::new();
        anim.push_frame(1.0, vec![[5.0, 0.0, 0.0]]);
        anim.push_frame(2.0, vec![[10.0, 0.0, 0.0]]);
        let frame = anim.interpolate_frame(0.0);
        assert!((frame[0][0] - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_particle_animation_total_duration() {
        let mut anim = ParticleAnimation::new();
        anim.push_frame(0.0, vec![[0.0, 0.0, 0.0]]);
        anim.push_frame(5.0, vec![[0.0, 0.0, 0.0]]);
        assert!((anim.total_duration() - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_particle_animation_n_particles() {
        let mut anim = ParticleAnimation::new();
        anim.push_frame(0.0, vec![[0.0; 3]; 10]);
        assert_eq!(anim.n_particles(), 10);
    }

    // ── AnimationClip ─────────────────────────────────────────────────────────

    #[test]
    fn test_clip_loop_wraps() {
        let clip = AnimationClip::new("test", 2.0);
        assert!((clip.loop_animation(2.5) - 0.5).abs() < 1e-12);
        assert!((clip.loop_animation(0.3) - 0.3).abs() < 1e-12);
    }

    #[test]
    fn test_clip_loop_zero_duration() {
        let clip = AnimationClip::new("empty", 0.0);
        assert_eq!(clip.loop_animation(10.0), 0.0);
    }

    #[test]
    fn test_clip_ping_pong_forward() {
        let clip = AnimationClip::new("pp", 2.0);
        assert!((clip.ping_pong(0.5) - 0.5).abs() < 1e-12);
        assert!((clip.ping_pong(1.5) - 1.5).abs() < 1e-12);
    }

    #[test]
    fn test_clip_ping_pong_backward() {
        let clip = AnimationClip::new("pp", 2.0);
        // t=2.5 → wrapped = 2.5 (in second half) → period - wrapped = 4.0 - 2.5 = 1.5
        assert!((clip.ping_pong(2.5) - 1.5).abs() < 1e-12);
    }

    #[test]
    fn test_clip_ping_pong_zero_duration() {
        let clip = AnimationClip::new("empty", 0.0);
        assert_eq!(clip.ping_pong(5.0), 0.0);
    }

    #[test]
    fn test_clip_evaluate_all() {
        let mut clip = AnimationClip::new("c", 1.0);
        clip.add_channel(AnimationChannel::from_pairs(&[(0.0, 0.0), (1.0, 10.0)]));
        clip.add_channel(AnimationChannel::from_pairs(&[(0.0, 5.0), (1.0, 5.0)]));
        let vals = clip.evaluate_all(0.5);
        assert_eq!(vals.len(), 2);
        assert!((vals[0] - 5.0).abs() < 1e-12);
        assert!((vals[1] - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_clip_channel_count() {
        let mut clip = AnimationClip::new("c", 1.0);
        assert_eq!(clip.channel_count(), 0);
        clip.add_channel(AnimationChannel::default());
        assert_eq!(clip.channel_count(), 1);
    }

    // ── EasingFunction ────────────────────────────────────────────────────────

    #[test]
    fn test_ease_linear_identity() {
        for i in 0..=10 {
            let t = i as f64 / 10.0;
            assert!((EasingFunction::ease_linear(t) - t).abs() < 1e-15);
        }
    }

    #[test]
    fn test_ease_in_quad_endpoints() {
        assert!(EasingFunction::ease_in_quad(0.0).abs() < 1e-15);
        assert!((EasingFunction::ease_in_quad(1.0) - 1.0).abs() < 1e-15);
    }

    #[test]
    fn test_ease_out_quad_endpoints() {
        assert!(EasingFunction::ease_out_quad(0.0).abs() < 1e-15);
        assert!((EasingFunction::ease_out_quad(1.0) - 1.0).abs() < 1e-15);
    }

    #[test]
    fn test_ease_in_out_cubic_endpoints() {
        assert!(EasingFunction::ease_in_out_cubic(0.0).abs() < 1e-15);
        assert!((EasingFunction::ease_in_out_cubic(1.0) - 1.0).abs() < 1e-15);
    }

    #[test]
    fn test_ease_in_out_cubic_midpoint_is_half() {
        // The cubic ease-in-out S-curve passes through 0.5 at t=0.5
        let v = EasingFunction::ease_in_out_cubic(0.5);
        assert!((v - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_ease_in_elastic_endpoints() {
        assert!(EasingFunction::ease_in_elastic(0.0).abs() < 1e-12);
        assert!((EasingFunction::ease_in_elastic(1.0) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_ease_out_bounce_endpoints() {
        assert!(EasingFunction::ease_out_bounce(0.0).abs() < 1e-12);
        assert!((EasingFunction::ease_out_bounce(1.0) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_ease_in_out_back_endpoints() {
        let a = EasingFunction::ease_in_out_back(0.0);
        let b = EasingFunction::ease_in_out_back(1.0);
        assert!(a.abs() < 1e-12);
        assert!((b - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_spring_damped_start_zero() {
        let v = EasingFunction::spring_damped(0.0, 10.0, 1.0);
        assert!(v.abs() < 1e-12);
    }

    #[test]
    fn test_spring_damped_approaches_one() {
        // At very large t the damped spring should be close to 1
        let v = EasingFunction::spring_damped(100.0, 10.0, 2.0);
        assert!(
            (v - 1.0).abs() < 0.01,
            "spring at t=100 should be ≈1, got {v}"
        );
    }

    #[test]
    fn test_ease_in_cubic_endpoints() {
        assert!(EasingFunction::ease_in_cubic(0.0).abs() < 1e-15);
        assert!((EasingFunction::ease_in_cubic(1.0) - 1.0).abs() < 1e-15);
    }

    #[test]
    fn test_ease_out_cubic_endpoints() {
        assert!(EasingFunction::ease_out_cubic(0.0).abs() < 1e-15);
        assert!((EasingFunction::ease_out_cubic(1.0) - 1.0).abs() < 1e-15);
    }

    #[test]
    fn test_ease_in_out_quad_midpoint() {
        let v = EasingFunction::ease_in_out_quad(0.5);
        assert!(
            (v - 0.5).abs() < 1e-12,
            "ease_in_out_quad(0.5) should be 0.5, got {v}"
        );
    }

    #[test]
    fn test_ease_out_elastic_endpoints() {
        assert!(EasingFunction::ease_out_elastic(0.0).abs() < 1e-12);
        assert!((EasingFunction::ease_out_elastic(1.0) - 1.0).abs() < 1e-12);
    }

    // ── PhysicsReplay ─────────────────────────────────────────────────────────

    #[test]
    fn test_replay_empty() {
        let replay = PhysicsReplay::new();
        assert!(replay.replay_at(0.0).is_none());
    }

    #[test]
    fn test_replay_single_step() {
        let mut replay = PhysicsReplay::new();
        replay.record_step(vec![[1.0, 2.0, 3.0]], vec![[0.1, 0.2, 0.3]], 0.0);
        let (pos, vel) = replay.replay_at(0.0).unwrap();
        assert!((pos[0][0] - 1.0).abs() < 1e-12);
        assert!((vel[0][1] - 0.2).abs() < 1e-12);
    }

    #[test]
    fn test_replay_interpolation() {
        let mut replay = PhysicsReplay::new();
        replay.record_step(vec![[0.0, 0.0, 0.0]], vec![[0.0, 0.0, 0.0]], 0.0);
        replay.record_step(vec![[2.0, 0.0, 0.0]], vec![[2.0, 0.0, 0.0]], 1.0);
        let (pos, vel) = replay.replay_at(0.5).unwrap();
        assert!((pos[0][0] - 1.0).abs() < 1e-12);
        assert!((vel[0][0] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_replay_clamp_before() {
        let mut replay = PhysicsReplay::new();
        replay.record_step(vec![[5.0, 0.0, 0.0]], vec![[0.0; 3]], 1.0);
        let (pos, _) = replay.replay_at(0.0).unwrap();
        assert!((pos[0][0] - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_replay_clamp_after() {
        let mut replay = PhysicsReplay::new();
        replay.record_step(vec![[0.0; 3]], vec![[0.0; 3]], 0.0);
        replay.record_step(vec![[9.0, 0.0, 0.0]], vec![[0.0; 3]], 1.0);
        let (pos, _) = replay.replay_at(2.0).unwrap();
        assert!((pos[0][0] - 9.0).abs() < 1e-12);
    }

    #[test]
    fn test_replay_total_duration() {
        let mut replay = PhysicsReplay::new();
        replay.record_step(vec![], vec![], 0.0);
        replay.record_step(vec![], vec![], 3.0);
        assert!((replay.total_duration() - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_replay_export_frames_count() {
        let mut replay = PhysicsReplay::new();
        replay.record_step(vec![[0.0; 3]], vec![[0.0; 3]], 0.0);
        replay.record_step(vec![[1.0, 0.0, 0.0]], vec![[0.0; 3]], 1.0);
        let frames = replay.export_frames(10.0); // 10 fps over 1 second → 11 frames
        assert!(
            frames.len() >= 10,
            "expected ≥10 frames, got {}",
            frames.len()
        );
    }

    #[test]
    fn test_replay_export_frames_zero_fps() {
        let mut replay = PhysicsReplay::new();
        replay.record_step(vec![[0.0; 3]], vec![[0.0; 3]], 0.0);
        assert!(replay.export_frames(0.0).is_empty());
    }

    #[test]
    fn test_replay_n_particles() {
        let mut replay = PhysicsReplay::new();
        replay.record_step(vec![[0.0; 3]; 5], vec![[0.0; 3]; 5], 0.0);
        assert_eq!(replay.n_particles(), 5);
    }

    #[test]
    fn test_replay_n_steps() {
        let mut replay = PhysicsReplay::new();
        for i in 0..7 {
            replay.record_step(vec![], vec![], i as f64 * 0.1);
        }
        assert_eq!(replay.n_steps(), 7);
    }

    // ── Helpers ───────────────────────────────────────────────────────────────

    #[test]
    fn test_lerp3_endpoints() {
        let a = [1.0, 2.0, 3.0];
        let b = [4.0, 5.0, 6.0];
        let r0 = lerp3(a, b, 0.0);
        let r1 = lerp3(a, b, 1.0);
        assert_eq!(r0, a);
        assert_eq!(r1, b);
    }

    #[test]
    fn test_lerp3_midpoint() {
        let a = [0.0, 0.0, 0.0];
        let b = [2.0, 4.0, 6.0];
        let r = lerp3(a, b, 0.5);
        assert!((r[0] - 1.0).abs() < 1e-12);
        assert!((r[1] - 2.0).abs() < 1e-12);
        assert!((r[2] - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_distance3_zero() {
        let a = [1.0, 2.0, 3.0];
        assert!(distance3(a, a).abs() < 1e-12);
    }

    #[test]
    fn test_distance3_unit() {
        let a = [0.0; 3];
        let b = [1.0, 0.0, 0.0];
        assert!((distance3(a, b) - 1.0).abs() < 1e-12);
    }

    // ── InterpolationType ─────────────────────────────────────────────────────

    #[test]
    fn test_interp_type_eq() {
        assert_eq!(InterpolationType::Linear, InterpolationType::Linear);
        assert_ne!(InterpolationType::Linear, InterpolationType::CubicSpline);
    }

    // ── Keyframe ──────────────────────────────────────────────────────────────

    #[test]
    fn test_keyframe_new_no_tangents() {
        let kf: Keyframe<f64> = Keyframe::new(1.0, 5.0);
        assert!(kf.tangent_in.is_none());
        assert!(kf.tangent_out.is_none());
    }

    #[test]
    fn test_keyframe_with_tangents() {
        let kf = Keyframe::with_tangents(0.0, 1.0, -0.5, 0.5);
        assert_eq!(kf.tangent_in, Some(-0.5));
        assert_eq!(kf.tangent_out, Some(0.5));
    }

    // ── SimpleKeyframe / lerp_keyframes ───────────────────────────────────────

    #[test]
    fn test_lerp_keyframes_before_first() {
        let kfs = vec![
            SimpleKeyframe::new(1.0, 10.0),
            SimpleKeyframe::new(2.0, 20.0),
        ];
        assert!((lerp_keyframes(&kfs, 0.0) - 10.0).abs() < 1e-12);
    }

    #[test]
    fn test_lerp_keyframes_after_last() {
        let kfs = vec![
            SimpleKeyframe::new(0.0, 0.0),
            SimpleKeyframe::new(1.0, 100.0),
        ];
        assert!((lerp_keyframes(&kfs, 5.0) - 100.0).abs() < 1e-12);
    }

    #[test]
    fn test_lerp_keyframes_midpoint() {
        let kfs = vec![
            SimpleKeyframe::new(0.0, 0.0),
            SimpleKeyframe::new(2.0, 10.0),
        ];
        let v = lerp_keyframes(&kfs, 1.0);
        assert!((v - 5.0).abs() < 1e-12, "expected 5.0, got {v}");
    }

    #[test]
    fn test_lerp_keyframes_at_exact_time() {
        let kfs = vec![
            SimpleKeyframe::new(0.0, 1.0),
            SimpleKeyframe::new(1.0, 3.0),
            SimpleKeyframe::new(2.0, 7.0),
        ];
        assert!((lerp_keyframes(&kfs, 1.0) - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_lerp_keyframes_empty() {
        assert_eq!(lerp_keyframes(&[], 1.0), 0.0);
    }

    // ── cubic_hermite_interp ──────────────────────────────────────────────────

    #[test]
    fn test_hermite_endpoints() {
        // t=0 → p0, t=1 → p1
        assert!((cubic_hermite_interp(3.0, 7.0, 0.0, 0.0, 0.0) - 3.0).abs() < 1e-12);
        assert!((cubic_hermite_interp(3.0, 7.0, 0.0, 0.0, 1.0) - 7.0).abs() < 1e-12);
    }

    #[test]
    fn test_hermite_midpoint_no_tangents() {
        // With zero tangents, Hermite is similar to linear at midpoint
        let v = cubic_hermite_interp(0.0, 2.0, 0.0, 0.0, 0.5);
        assert!(
            (0.0..=2.0).contains(&v),
            "midpoint should be in range [0,2]"
        );
    }

    #[test]
    fn test_hermite_nonzero_tangent_effect() {
        // Non-zero tangents should change the midpoint value
        let v_linear = cubic_hermite_interp(0.0, 1.0, 0.0, 0.0, 0.5);
        let v_curved = cubic_hermite_interp(0.0, 1.0, 2.0, 0.0, 0.5);
        assert!(
            (v_linear - v_curved).abs() > 1e-10,
            "tangents should affect the curve"
        );
    }

    // ── bezier_cubic_1d ───────────────────────────────────────────────────────

    #[test]
    fn test_bezier_cubic_endpoints() {
        assert!((bezier_cubic_1d(0.0, 1.0, 2.0, 3.0, 0.0) - 0.0).abs() < 1e-12);
        assert!((bezier_cubic_1d(0.0, 1.0, 2.0, 3.0, 1.0) - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_bezier_cubic_symmetry() {
        // For symmetric control points, midpoint should match
        let v = bezier_cubic_1d(0.0, 1.0, 1.0, 0.0, 0.5);
        // Symmetric: should be 0.75
        assert!((v - 0.75).abs() < 1e-12, "expected 0.75, got {v}");
    }

    // ── AnimationCurve / eval_curve ───────────────────────────────────────────

    #[test]
    fn test_eval_curve_empty() {
        let curve = AnimationCurve::new();
        assert_eq!(eval_curve(&curve, 0.5), 0.0);
    }

    #[test]
    fn test_eval_curve_before_first_keyframe() {
        let mut curve = AnimationCurve::new();
        curve.push(SimpleKeyframe::new(1.0, 5.0), 0.0);
        curve.push(SimpleKeyframe::new(2.0, 10.0), 0.0);
        assert!((eval_curve(&curve, 0.0) - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_eval_curve_after_last_keyframe() {
        let mut curve = AnimationCurve::new();
        curve.push(SimpleKeyframe::new(0.0, 1.0), 0.0);
        curve.push(SimpleKeyframe::new(1.0, 2.0), 0.0);
        assert!((eval_curve(&curve, 5.0) - 2.0).abs() < 1e-12);
    }

    // ── easing functions ──────────────────────────────────────────────────────

    #[test]
    fn test_ease_in_out_quad_fn_endpoints() {
        assert!((ease_in_out_quad(0.0) - 0.0).abs() < 1e-12);
        assert!((ease_in_out_quad(1.0) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_ease_in_out_quad_fn_midpoint() {
        assert!((ease_in_out_quad(0.5) - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_ease_in_out_cubic_fn_endpoints() {
        assert!((ease_in_out_cubic(0.0) - 0.0).abs() < 1e-12);
        assert!((ease_in_out_cubic(1.0) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_ease_in_expo_at_zero() {
        assert!((ease_in_expo(0.0) - 0.0).abs() < 1e-12);
    }

    #[test]
    fn test_ease_in_expo_at_one() {
        assert!((ease_in_expo(1.0) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_ease_out_bounce_fn_endpoints() {
        assert!((ease_out_bounce(0.0) - 0.0).abs() < 1e-10);
        assert!((ease_out_bounce(1.0) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_easing_output_in_range() {
        for i in 0..=10 {
            let t = i as f64 / 10.0;
            assert!(ease_in_out_quad(t) >= 0.0 && ease_in_out_quad(t) <= 1.0 + 1e-10);
            assert!(ease_in_out_cubic(t) >= 0.0 && ease_in_out_cubic(t) <= 1.0 + 1e-10);
            assert!(ease_out_bounce(t) >= -0.01 && ease_out_bounce(t) <= 1.1);
        }
    }

    // ── spring_animation ──────────────────────────────────────────────────────

    #[test]
    fn test_spring_moves_toward_target() {
        let (new_pos, _) = spring_animation(0.0, 0.0, 10.0, 100.0, 5.0, 0.016);
        assert!(
            new_pos > 0.0,
            "spring should move toward target, got {new_pos}"
        );
    }

    #[test]
    fn test_spring_at_target_stays() {
        // At rest at target: no net force
        let (new_pos, new_vel) = spring_animation(1.0, 0.0, 1.0, 50.0, 5.0, 0.016);
        assert!(
            (new_pos - 1.0).abs() < 1e-12,
            "at target, position should not change"
        );
        assert!(new_vel.abs() < 1e-12, "at target, velocity should stay 0");
    }

    #[test]
    fn test_spring_returns_tuple() {
        let result = spring_animation(0.0, 0.0, 5.0, 10.0, 1.0, 0.1);
        let (pos, vel) = result;
        assert!(pos.is_finite());
        assert!(vel.is_finite());
    }

    // ── ik_two_bone ───────────────────────────────────────────────────────────

    #[test]
    fn test_ik_two_bone_straight_reach() {
        // Target directly reachable along x axis
        let root = [0.0, 0.0, 0.0];
        let mid = [1.0, 0.1, 0.0]; // slight bend
        let end = [2.0, 0.0, 0.0];
        let target = [1.5, 0.0, 0.0];
        let (new_mid, new_end) = ik_two_bone(root, mid, end, target, 1.0, 1.0);
        // End effector should reach target
        let dist_to_target = {
            let dx = new_end[0] - target[0];
            let dy = new_end[1] - target[1];
            let dz = new_end[2] - target[2];
            (dx * dx + dy * dy + dz * dz).sqrt()
        };
        assert!(
            dist_to_target < 0.01,
            "end should reach target, dist = {dist_to_target}"
        );

        // Root should not move
        let dist_root_mid = {
            let dx = new_mid[0] - root[0];
            let dy = new_mid[1] - root[1];
            let dz = new_mid[2] - root[2];
            (dx * dx + dy * dy + dz * dz).sqrt()
        };
        assert!(
            (dist_root_mid - 1.0).abs() < 0.01,
            "root→mid distance should be l1=1.0, got {dist_root_mid}"
        );
    }

    #[test]
    fn test_ik_two_bone_fully_extended() {
        // Target beyond reach: fully extend
        let root = [0.0, 0.0, 0.0];
        let mid = [1.0, 0.0, 0.0];
        let end = [2.0, 0.0, 0.0];
        let target = [10.0, 0.0, 0.0];
        let (new_mid, new_end) = ik_two_bone(root, mid, end, target, 1.0, 1.0);
        // Should fully extend: new_end = [2, 0, 0]
        assert!((new_mid[0] - 1.0).abs() < 0.01, "mid should be at x=1");
        assert!(
            (new_end[0] - 2.0).abs() < 0.01,
            "end should be at x=2 (fully extended)"
        );
    }
}
