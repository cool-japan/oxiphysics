// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Keyframe animation tracks for scripted body motion.
//!
//! Provides timeline-driven position and rotation tracks for animating physics
//! bodies along predetermined paths (cinematics, scripted triggers, NPC motion,
//! cutscenes, test rigs, etc.).
//!
//! Two concrete track types are provided:
//!
//! - `Vec3Track` — samples a `[f64; 3]` value (position, scale, velocity) over
//!   time using linear or eased interpolation.
//! - `QuatTrack` — samples a unit quaternion `[f64; 4]` using **spherical linear
//!   interpolation** (SLERP) with automatic short-path correction.
//!
//! Tracks are combined into a `BodyAnimation` and grouped into an
//! `AnimationClip`.  An `AnimationPlayer` advances the clip's timeline and
//! returns per-body `AnimBodyState`s that can drive kinematic bodies or blend
//! with live physics state.
//!
//! ## Example
//!
//! ```rust
//! use oxiphysics::animation::{Vec3Track, EaseKind, AnimationPlayer, AnimationClip, BodyAnimation};
//!
//! let mut pos = Vec3Track::new(false);
//! pos.push_keyframe(0.0, [0.0, 0.0, 0.0], EaseKind::Linear);
//! pos.push_keyframe(1.0, [0.0, 5.0, 0.0], EaseKind::Linear);
//! pos.push_keyframe(2.0, [0.0, 0.0, 0.0], EaseKind::Linear);
//!
//! let body_anim = BodyAnimation::new("ball").with_position(pos);
//! let clip = AnimationClip::new("rise_and_fall", vec![body_anim]);
//! let mut player = AnimationPlayer::new(clip);
//! player.play();
//!
//! let states = player.update(0.5); // advance 0.5 s
//! assert_eq!(states.len(), 1);
//! // At t=0.5, linear between [0,0,0] and [0,5,0] → [0,2.5,0]
//! assert!((states[0].position.unwrap()[1] - 2.5).abs() < 0.01);
//! ```

// ============================================================================
// EaseKind
// ============================================================================

/// Easing function applied when interpolating *from* a keyframe toward the next.
///
/// The easing of the *source* keyframe governs the transition into the next one.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum EaseKind {
    /// Constant-speed linear interpolation.
    #[default]
    Linear,
    /// Cubic ease-in (slow start, fast end).
    EaseIn,
    /// Cubic ease-out (fast start, slow end).
    EaseOut,
    /// Cubic smooth-step (slow–fast–slow).
    EaseInOut,
    /// Hold the value until the next keyframe; then jump instantly.
    Step,
}

impl EaseKind {
    /// Maps a normalised time `t ∈ [0, 1]` through the easing curve.
    pub fn apply(self, t: f64) -> f64 {
        let t = t.clamp(0.0, 1.0);
        match self {
            EaseKind::Linear => t,
            EaseKind::EaseIn => t * t * t,
            EaseKind::EaseOut => {
                let u = 1.0 - t;
                1.0 - u * u * u
            }
            EaseKind::EaseInOut => {
                if t < 0.5 {
                    4.0 * t * t * t
                } else {
                    let u = -2.0 * t + 2.0;
                    1.0 - u * u * u / 2.0
                }
            }
            EaseKind::Step => {
                if t < 1.0 {
                    0.0
                } else {
                    1.0
                }
            }
        }
    }
}

// ============================================================================
// Vec3Track
// ============================================================================

/// A single keyframe in a [`Vec3Track`].
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct Vec3Keyframe {
    /// Time in seconds.
    pub time: f64,
    /// Value at this keyframe.
    pub value: [f64; 3],
    /// Easing applied when leaving this keyframe toward the next.
    pub ease: EaseKind,
}

/// Timeline-driven track producing `[f64; 3]` samples by interpolating between
/// keyframes.
///
/// Push keyframes in ascending time order with [`push_keyframe`](Vec3Track::push_keyframe).
#[derive(Clone, Debug, Default, serde::Serialize, serde::Deserialize)]
pub struct Vec3Track {
    /// Keyframes, sorted by time.
    pub keyframes: Vec<Vec3Keyframe>,
    /// Whether the track loops (time wraps after the last keyframe).
    pub looping: bool,
}

impl Vec3Track {
    /// Creates an empty track.
    pub fn new(looping: bool) -> Self {
        Self {
            keyframes: Vec::new(),
            looping,
        }
    }

    /// Appends a keyframe.  Must be added in non-decreasing time order.
    pub fn push_keyframe(&mut self, time: f64, value: [f64; 3], ease: EaseKind) {
        self.keyframes.push(Vec3Keyframe { time, value, ease });
    }

    /// Total duration (time of the last keyframe, or `0` if empty).
    pub fn duration(&self) -> f64 {
        self.keyframes.last().map_or(0.0, |k| k.time)
    }

    /// Number of keyframes.
    pub fn len(&self) -> usize {
        self.keyframes.len()
    }

    /// Returns `true` when there are no keyframes.
    pub fn is_empty(&self) -> bool {
        self.keyframes.is_empty()
    }

    /// Samples the track at time `t` (seconds).
    ///
    /// - Before the first keyframe: clamps to the first value.
    /// - After the last keyframe (non-looping): clamps to the last value.
    /// - After the last keyframe (looping): wraps `t` modulo duration.
    pub fn sample(&self, mut t: f64) -> [f64; 3] {
        if self.keyframes.is_empty() {
            return [0.0; 3];
        }
        if self.keyframes.len() == 1 {
            return self.keyframes[0].value;
        }

        let dur = self.duration();
        if self.looping && dur > 0.0 {
            t = t.rem_euclid(dur);
        }

        let first = self
            .keyframes
            .first()
            .expect("at least two keyframes, checked above");
        if t <= first.time {
            return first.value;
        }
        let last = self
            .keyframes
            .last()
            .expect("at least two keyframes, checked above");
        if t >= last.time {
            return last.value;
        }

        // Find the index of the last keyframe with time ≤ t
        let idx = self
            .keyframes
            .partition_point(|k| k.time <= t)
            .saturating_sub(1);

        let a = &self.keyframes[idx];
        let b = &self.keyframes[idx + 1];
        let span = b.time - a.time;
        let raw_t = if span > 1e-12 {
            (t - a.time) / span
        } else {
            1.0
        };
        let et = a.ease.apply(raw_t);

        [
            a.value[0] + (b.value[0] - a.value[0]) * et,
            a.value[1] + (b.value[1] - a.value[1]) * et,
            a.value[2] + (b.value[2] - a.value[2]) * et,
        ]
    }
}

// ============================================================================
// QuatTrack — SLERP with short-path correction
// ============================================================================

/// A single keyframe in a [`QuatTrack`].
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct QuatKeyframe {
    /// Time in seconds.
    pub time: f64,
    /// Unit quaternion `[x, y, z, w]`.
    pub value: [f64; 4],
    /// Easing applied when leaving this keyframe toward the next.
    pub ease: EaseKind,
}

/// Timeline-driven track producing `[f64; 4]` unit quaternion samples using
/// **spherical linear interpolation (SLERP)** with automatic short-path
/// correction.
///
/// The quaternion convention is `[x, y, z, w]`.  Stored values should be
/// unit-length; the track normalises sampled results but does not validate inputs.
#[derive(Clone, Debug, Default, serde::Serialize, serde::Deserialize)]
pub struct QuatTrack {
    /// Keyframes, sorted by time.
    pub keyframes: Vec<QuatKeyframe>,
    /// Whether the track loops.
    pub looping: bool,
}

impl QuatTrack {
    /// Creates an empty quaternion track.
    pub fn new(looping: bool) -> Self {
        Self {
            keyframes: Vec::new(),
            looping,
        }
    }

    /// Appends a keyframe in non-decreasing time order.
    pub fn push_keyframe(&mut self, time: f64, value: [f64; 4], ease: EaseKind) {
        self.keyframes.push(QuatKeyframe { time, value, ease });
    }

    /// Total duration.
    pub fn duration(&self) -> f64 {
        self.keyframes.last().map_or(0.0, |k| k.time)
    }

    /// Number of keyframes.
    pub fn len(&self) -> usize {
        self.keyframes.len()
    }

    /// Returns `true` when there are no keyframes.
    pub fn is_empty(&self) -> bool {
        self.keyframes.is_empty()
    }

    /// Samples the track at time `t`, returning a unit quaternion `[x, y, z, w]`.
    pub fn sample(&self, mut t: f64) -> [f64; 4] {
        const IDENTITY: [f64; 4] = [0.0, 0.0, 0.0, 1.0];

        if self.keyframes.is_empty() {
            return IDENTITY;
        }
        if self.keyframes.len() == 1 {
            return self.keyframes[0].value;
        }

        let dur = self.duration();
        if self.looping && dur > 0.0 {
            t = t.rem_euclid(dur);
        }

        let first = self
            .keyframes
            .first()
            .expect("at least two keyframes, checked above");
        if t <= first.time {
            return first.value;
        }
        let last = self
            .keyframes
            .last()
            .expect("at least two keyframes, checked above");
        if t >= last.time {
            return last.value;
        }

        let idx = self
            .keyframes
            .partition_point(|k| k.time <= t)
            .saturating_sub(1);

        let a = &self.keyframes[idx];
        let b = &self.keyframes[idx + 1];
        let span = b.time - a.time;
        let raw_t = if span > 1e-12 {
            (t - a.time) / span
        } else {
            1.0
        };
        let et = a.ease.apply(raw_t);

        quat_slerp(a.value, b.value, et)
    }
}

// ============================================================================
// Quaternion SLERP
// ============================================================================

/// Spherical linear interpolation between two unit quaternions.
///
/// - Applies short-path correction: negates `b` when `dot(a, b) < 0`.
/// - Falls back to normalised LERP when the quaternions are nearly identical
///   (avoids division by zero near `sin(θ) ≈ 0`).
/// - Always returns a unit-length quaternion.
fn quat_slerp(a: [f64; 4], mut b: [f64; 4], t: f64) -> [f64; 4] {
    let mut dot = a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3];

    // Short-path correction: take the arc with the smaller angle
    if dot < 0.0 {
        b = [-b[0], -b[1], -b[2], -b[3]];
        dot = -dot;
    }

    const SLERP_THRESHOLD: f64 = 0.9995;
    let raw = if dot > SLERP_THRESHOLD {
        // Quaternions nearly identical — use LERP to avoid numerical singularity
        [
            a[0] + t * (b[0] - a[0]),
            a[1] + t * (b[1] - a[1]),
            a[2] + t * (b[2] - a[2]),
            a[3] + t * (b[3] - a[3]),
        ]
    } else {
        let theta_0 = dot.clamp(-1.0, 1.0).acos(); // angle between a and b
        let theta = theta_0 * t;
        let sin_theta = theta.sin();
        let sin_theta_0 = theta_0.sin();
        let s0 = theta.cos() - dot * sin_theta / sin_theta_0;
        let s1 = sin_theta / sin_theta_0;
        [
            s0 * a[0] + s1 * b[0],
            s0 * a[1] + s1 * b[1],
            s0 * a[2] + s1 * b[2],
            s0 * a[3] + s1 * b[3],
        ]
    };

    // Normalise
    let len = (raw[0] * raw[0] + raw[1] * raw[1] + raw[2] * raw[2] + raw[3] * raw[3]).sqrt();
    if len > 1e-12 {
        [raw[0] / len, raw[1] / len, raw[2] / len, raw[3] / len]
    } else {
        [0.0, 0.0, 0.0, 1.0] // identity fallback
    }
}

// ============================================================================
// BodyAnimation
// ============================================================================

/// All animation tracks for a single named body.
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct BodyAnimation {
    /// Target body identifier (matches `SceneBody::id` or an app-level id).
    pub body_id: String,
    /// Optional position track (world-space).
    pub position: Option<Vec3Track>,
    /// Optional rotation track (unit quaternion `[x, y, z, w]`).
    pub rotation: Option<QuatTrack>,
    /// Optional scale track.
    pub scale: Option<Vec3Track>,
}

impl BodyAnimation {
    /// Creates a new `BodyAnimation` with no tracks.
    pub fn new(body_id: impl Into<String>) -> Self {
        Self {
            body_id: body_id.into(),
            position: None,
            rotation: None,
            scale: None,
        }
    }

    /// Attaches a position track (builder-chain).
    pub fn with_position(mut self, track: Vec3Track) -> Self {
        self.position = Some(track);
        self
    }

    /// Attaches a rotation track (builder-chain).
    pub fn with_rotation(mut self, track: QuatTrack) -> Self {
        self.rotation = Some(track);
        self
    }

    /// Attaches a scale track (builder-chain).
    pub fn with_scale(mut self, track: Vec3Track) -> Self {
        self.scale = Some(track);
        self
    }

    /// Maximum duration across all attached tracks.
    pub fn duration(&self) -> f64 {
        let pd = self.position.as_ref().map_or(0.0, |t| t.duration());
        let rd = self.rotation.as_ref().map_or(0.0, |t| t.duration());
        let sd = self.scale.as_ref().map_or(0.0, |t| t.duration());
        pd.max(rd).max(sd)
    }
}

// ============================================================================
// AnimationClip
// ============================================================================

/// A named collection of [`BodyAnimation`] tracks that play in sync.
///
/// The clip's duration is the maximum duration across all body animations.
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct AnimationClip {
    /// Human-readable name.
    pub name: String,
    /// Per-body animation tracks.
    pub body_animations: Vec<BodyAnimation>,
}

impl AnimationClip {
    /// Creates a new clip.
    pub fn new(name: impl Into<String>, body_animations: Vec<BodyAnimation>) -> Self {
        Self {
            name: name.into(),
            body_animations,
        }
    }

    /// Total duration (maximum across all body animations).
    pub fn duration(&self) -> f64 {
        self.body_animations
            .iter()
            .map(|ba| ba.duration())
            .fold(0.0_f64, f64::max)
    }

    /// Returns the body animation for `body_id`, if any.
    pub fn body_anim(&self, body_id: &str) -> Option<&BodyAnimation> {
        self.body_animations.iter().find(|ba| ba.body_id == body_id)
    }

    /// Number of body animations in this clip.
    pub fn body_count(&self) -> usize {
        self.body_animations.len()
    }
}

impl std::fmt::Display for AnimationClip {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "AnimationClip {{ name: {:?}, bodies: {}, duration: {:.3}s }}",
            self.name,
            self.body_animations.len(),
            self.duration()
        )
    }
}

// ============================================================================
// AnimBodyState — sampled output per update tick
// ============================================================================

/// Sampled state of a single body at the current animation time.
#[derive(Clone, Debug)]
pub struct AnimBodyState {
    /// Body identifier.
    pub body_id: String,
    /// Sampled world-space position (`None` when no position track exists).
    pub position: Option<[f64; 3]>,
    /// Sampled rotation quaternion `[x, y, z, w]` (`None` when no rotation track
    /// exists).
    pub rotation: Option<[f64; 4]>,
    /// Sampled scale (`None` when no scale track exists).
    pub scale: Option<[f64; 3]>,
}

// ============================================================================
// AnimationPlayer
// ============================================================================

/// Plays an [`AnimationClip`] forward in time, producing [`AnimBodyState`]s each
/// step.
///
/// ## Controls
///
/// | Method | Effect |
/// |--------|--------|
/// | [`play`](AnimationPlayer::play) | Start / resume playback |
/// | [`pause`](AnimationPlayer::pause) | Freeze time (retains position) |
/// | [`stop`](AnimationPlayer::stop) | Freeze and rewind to t=0 |
/// | [`seek`](AnimationPlayer::seek) | Jump to an absolute time |
///
/// ## Example
///
/// ```rust
/// use oxiphysics::animation::{Vec3Track, EaseKind, AnimationPlayer, AnimationClip, BodyAnimation};
///
/// let mut t = Vec3Track::new(false);
/// t.push_keyframe(0.0, [0.0, 0.0, 0.0], EaseKind::Linear);
/// t.push_keyframe(2.0, [10.0, 0.0, 0.0], EaseKind::Linear);
///
/// let clip = AnimationClip::new("slide", vec![BodyAnimation::new("box").with_position(t)]);
/// let mut player = AnimationPlayer::new(clip);
/// player.play();
/// player.update(1.0);
/// assert!(!player.is_finished()); // still playing at t=1
/// player.update(1.0);
/// assert!(player.is_finished()); // reached t=2
/// ```
#[derive(Clone, Debug)]
pub struct AnimationPlayer {
    /// The clip being played.
    pub clip: AnimationClip,
    /// Current playback time in seconds.
    pub time: f64,
    /// Playback speed multiplier (`1.0` = real-time, `2.0` = double-speed).
    pub speed: f64,
    playing: bool,
    /// Whether to loop when reaching the end.
    pub looping: bool,
}

impl AnimationPlayer {
    /// Creates a player in the stopped state at `t = 0`.
    pub fn new(clip: AnimationClip) -> Self {
        Self {
            clip,
            time: 0.0,
            speed: 1.0,
            playing: false,
            looping: false,
        }
    }

    /// Starts or resumes playback.
    pub fn play(&mut self) {
        self.playing = true;
    }

    /// Pauses playback (retains current time).
    pub fn pause(&mut self) {
        self.playing = false;
    }

    /// Stops playback and rewinds to `t = 0`.
    pub fn stop(&mut self) {
        self.playing = false;
        self.time = 0.0;
    }

    /// Seeks to an absolute time (clamped to clip duration).
    pub fn seek(&mut self, t: f64) {
        self.time = t.clamp(0.0, self.clip.duration());
    }

    /// Returns `true` when the player is actively advancing time.
    pub fn is_playing(&self) -> bool {
        self.playing
    }

    /// Playback progress as a fraction `[0, 1]`.
    pub fn progress(&self) -> f64 {
        let dur = self.clip.duration();
        if dur > 0.0 {
            (self.time / dur).clamp(0.0, 1.0)
        } else {
            1.0
        }
    }

    /// Returns `true` when `time` has reached or passed the clip duration.
    pub fn is_finished(&self) -> bool {
        self.time >= self.clip.duration()
    }

    /// Advances playback by `dt` seconds and returns sampled body states.
    ///
    /// If the player is paused or stopped, `dt` is ignored but states are still
    /// sampled at the current time (useful for scrubbing / preview).
    pub fn update(&mut self, dt: f64) -> Vec<AnimBodyState> {
        if self.playing {
            let dur = self.clip.duration();
            self.time += dt * self.speed;
            if self.looping && dur > 0.0 {
                self.time = self.time.rem_euclid(dur);
            } else {
                self.time = self.time.min(dur);
                if self.time >= dur {
                    self.playing = false;
                }
            }
        }
        self.sample_all()
    }

    /// Samples all body animations at the current time without advancing.
    fn sample_all(&self) -> Vec<AnimBodyState> {
        self.clip
            .body_animations
            .iter()
            .map(|ba| AnimBodyState {
                body_id: ba.body_id.clone(),
                position: ba.position.as_ref().map(|t| t.sample(self.time)),
                rotation: ba.rotation.as_ref().map(|t| t.sample(self.time)),
                scale: ba.scale.as_ref().map(|t| t.sample(self.time)),
            })
            .collect()
    }
}

impl std::fmt::Display for AnimationPlayer {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "AnimationPlayer {{ clip: {:?}, time: {:.3}s/{:.3}s, playing: {}, speed: {} }}",
            self.clip.name,
            self.time,
            self.clip.duration(),
            self.playing,
            self.speed
        )
    }
}

// ============================================================================
// Unit tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ease_kind_boundary() {
        assert!((EaseKind::Linear.apply(0.0)).abs() < 1e-9);
        assert!((EaseKind::Linear.apply(1.0) - 1.0).abs() < 1e-9);
        assert!((EaseKind::EaseInOut.apply(0.5) - 0.5).abs() < 1e-9);
        assert!((EaseKind::Step.apply(0.99)).abs() < 1e-9);
        assert!((EaseKind::Step.apply(1.0) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn vec3_track_clamp() {
        let mut t = Vec3Track::new(false);
        t.push_keyframe(1.0, [1.0, 0.0, 0.0], EaseKind::Linear);
        t.push_keyframe(2.0, [2.0, 0.0, 0.0], EaseKind::Linear);
        // Before first keyframe
        assert!((t.sample(0.0)[0] - 1.0).abs() < 1e-9);
        // After last keyframe
        assert!((t.sample(5.0)[0] - 2.0).abs() < 1e-9);
        // Midpoint
        assert!((t.sample(1.5)[0] - 1.5).abs() < 1e-9);
    }

    #[test]
    fn vec3_track_loop() {
        let mut t = Vec3Track::new(true);
        t.push_keyframe(0.0, [0.0, 0.0, 0.0], EaseKind::Linear);
        t.push_keyframe(1.0, [1.0, 0.0, 0.0], EaseKind::Linear);
        // t=1.5 wraps to t=0.5
        assert!((t.sample(1.5)[0] - 0.5).abs() < 1e-9);
    }

    #[test]
    fn quat_slerp_identity() {
        let id = [0.0_f64, 0.0, 0.0, 1.0];
        let r = quat_slerp(id, id, 0.5);
        assert!((r[3] - 1.0).abs() < 1e-6);
        // Must be unit length
        let len = (r[0] * r[0] + r[1] * r[1] + r[2] * r[2] + r[3] * r[3]).sqrt();
        assert!((len - 1.0).abs() < 1e-9);
    }

    #[test]
    fn quat_slerp_short_path() {
        // Negated quaternion represents same rotation — SLERP should still converge
        let a = [0.0_f64, 0.0, 0.0, 1.0];
        let b_neg = [0.0_f64, 0.0, 0.0, -1.0]; // same rotation, negated
        let r = quat_slerp(a, b_neg, 0.5);
        let len = (r[0] * r[0] + r[1] * r[1] + r[2] * r[2] + r[3] * r[3]).sqrt();
        assert!((len - 1.0).abs() < 1e-9);
    }

    #[test]
    fn animation_player_lifecycle() {
        let mut pos = Vec3Track::new(false);
        pos.push_keyframe(0.0, [0.0, 0.0, 0.0], EaseKind::Linear);
        pos.push_keyframe(2.0, [10.0, 0.0, 0.0], EaseKind::Linear);

        let ba = BodyAnimation::new("obj").with_position(pos);
        let clip = AnimationClip::new("move", vec![ba]);
        let mut player = AnimationPlayer::new(clip);

        // Stopped: update should not advance time
        player.update(1.0);
        assert!((player.time).abs() < 1e-9);

        player.play();
        let states = player.update(1.0);
        assert_eq!(states.len(), 1);
        assert!((states[0].position.unwrap()[0] - 5.0).abs() < 1e-9);
        assert!(!player.is_finished());

        player.update(1.0);
        assert!(player.is_finished());
        assert!(!player.is_playing());
    }
}
