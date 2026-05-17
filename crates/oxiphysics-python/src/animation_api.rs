// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Python bindings for the oxiphysics animation module.
//!
//! Exposes keyframe animation tracks (Vec3 + quaternion SLERP) to Python.

#![allow(missing_docs)]
#![allow(dead_code)]

use oxiphysics::animation::{AnimationClip, AnimationPlayer, EaseKind, QuatTrack, Vec3Track};
use pyo3::prelude::*;

// ─────────────────────────────────────────────────────────────────────────────
// PyVec3Track
// ─────────────────────────────────────────────────────────────────────────────

/// A keyframe track for Vec3 values (position, scale, etc).
#[pyclass(name = "Vec3Track")]
pub struct PyVec3Track {
    inner: Vec3Track,
}

#[pymethods]
impl PyVec3Track {
    /// Create a new Vec3 track. `looping` controls whether it wraps.
    #[new]
    pub fn new(looping: bool) -> Self {
        Self {
            inner: Vec3Track::new(looping),
        }
    }

    /// Add a keyframe at `time` with given `value` and easing.
    /// `ease` should be one of: "Linear", "EaseIn", "EaseOut", "EaseInOut", "Step".
    pub fn add_keyframe(&mut self, time: f64, value: [f64; 3], ease: &str) -> PyResult<()> {
        let ease_kind = parse_ease_kind(ease)?;
        self.inner.push_keyframe(time, value, ease_kind);
        Ok(())
    }

    /// Total duration of the track.
    pub fn duration(&self) -> f64 {
        self.inner.duration()
    }

    /// Number of keyframes.
    pub fn len(&self) -> usize {
        self.inner.len()
    }

    /// `true` if no keyframes.
    pub fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }

    /// Sample the track at time `t`.
    pub fn sample(&self, t: f64) -> [f64; 3] {
        self.inner.sample(t)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PyQuatTrack
// ─────────────────────────────────────────────────────────────────────────────

/// A keyframe track for quaternion values (orientation).
#[pyclass(name = "QuatTrack")]
pub struct PyQuatTrack {
    inner: QuatTrack,
}

#[pymethods]
impl PyQuatTrack {
    /// Create a new quaternion track.
    #[new]
    pub fn new(looping: bool) -> Self {
        Self {
            inner: QuatTrack::new(looping),
        }
    }

    /// Add a keyframe at `time` with quaternion `value [x,y,z,w]` and easing.
    pub fn add_keyframe(&mut self, time: f64, value: [f64; 4], ease: &str) -> PyResult<()> {
        let ease_kind = parse_ease_kind(ease)?;
        self.inner.push_keyframe(time, value, ease_kind);
        Ok(())
    }

    /// Total duration.
    pub fn duration(&self) -> f64 {
        self.inner.duration()
    }

    /// Sample the quaternion at time `t`. Returns `[x, y, z, w]`.
    pub fn sample(&self, t: f64) -> [f64; 4] {
        self.inner.sample(t)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PyAnimationPlayer
// ─────────────────────────────────────────────────────────────────────────────

/// Drives playback of an AnimationClip.
#[pyclass(name = "AnimationPlayer")]
pub struct PyAnimationPlayer {
    inner: AnimationPlayer,
}

#[pymethods]
impl PyAnimationPlayer {
    /// Create a player for the given clip (provided as JSON).
    #[new]
    pub fn new_from_json(clip_json: &str) -> PyResult<Self> {
        let clip: AnimationClip = serde_json::from_str(clip_json)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;
        Ok(Self {
            inner: AnimationPlayer::new(clip),
        })
    }

    /// Start or resume playback.
    pub fn play(&mut self) {
        self.inner.play();
    }

    /// Pause playback.
    pub fn pause(&mut self) {
        self.inner.pause();
    }

    /// Stop and reset to start.
    pub fn stop(&mut self) {
        self.inner.stop();
    }

    /// Whether the player is currently playing.
    pub fn is_playing(&self) -> bool {
        self.inner.is_playing()
    }

    /// Playback progress in [0, 1].
    pub fn progress(&self) -> f64 {
        self.inner.progress()
    }

    /// Whether the clip has finished.
    pub fn is_finished(&self) -> bool {
        self.inner.is_finished()
    }

    /// Advance by `dt` and return all body states as JSON.
    ///
    /// JSON: `[{"body_id": "...", "position": [x,y,z] | null, "rotation": [x,y,z,w] | null, "scale": [x,y,z] | null}]`
    pub fn step_json(&mut self, dt: f64) -> String {
        let states = self.inner.update(dt);
        let parts: Vec<String> = states
            .iter()
            .map(|s| {
                let pos = match s.position {
                    Some(p) => format!("[{},{},{}]", p[0], p[1], p[2]),
                    None => "null".to_owned(),
                };
                let rot = match s.rotation {
                    Some(r) => format!("[{},{},{},{}]", r[0], r[1], r[2], r[3]),
                    None => "null".to_owned(),
                };
                let scl = match s.scale {
                    Some(sc) => format!("[{},{},{}]", sc[0], sc[1], sc[2]),
                    None => "null".to_owned(),
                };
                format!(
                    "{{\"body_id\":{:?},\"position\":{},\"rotation\":{},\"scale\":{}}}",
                    s.body_id, pos, rot, scl
                )
            })
            .collect();
        format!("[{}]", parts.join(","))
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Helpers
// ─────────────────────────────────────────────────────────────────────────────

fn parse_ease_kind(s: &str) -> PyResult<EaseKind> {
    match s {
        "Linear" => Ok(EaseKind::Linear),
        "EaseIn" => Ok(EaseKind::EaseIn),
        "EaseOut" => Ok(EaseKind::EaseOut),
        "EaseInOut" => Ok(EaseKind::EaseInOut),
        "Step" => Ok(EaseKind::Step),
        other => Err(pyo3::exceptions::PyValueError::new_err(format!(
            "Unknown ease kind: {}. Use one of: Linear, EaseIn, EaseOut, EaseInOut, Step",
            other
        ))),
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Module registration
// ─────────────────────────────────────────────────────────────────────────────

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyVec3Track>()?;
    m.add_class::<PyQuatTrack>()?;
    m.add_class::<PyAnimationPlayer>()?;
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_vec3_track_instantiation() {
        let mut track = PyVec3Track::new(false);
        assert!(track.is_empty());
        track
            .add_keyframe(0.0, [0.0, 0.0, 0.0], "Linear")
            .expect("add_keyframe failed");
        track
            .add_keyframe(1.0, [1.0, 0.0, 0.0], "EaseInOut")
            .expect("add_keyframe failed");
        assert_eq!(track.len(), 2);
        let sample = track.sample(0.5);
        assert!((0.0..=1.0).contains(&sample[0]));
    }

    #[test]
    fn test_quat_track_instantiation() {
        let mut track = PyQuatTrack::new(false);
        track
            .add_keyframe(0.0, [0.0, 0.0, 0.0, 1.0], "Linear")
            .expect("add_keyframe failed");
        track
            .add_keyframe(1.0, [0.0, 0.0, 0.0, 1.0], "Linear")
            .expect("add_keyframe failed");
        assert_eq!(track.duration(), 1.0);
    }

    #[test]
    fn test_body_animation_json_round_trip() {
        use oxiphysics::animation::BodyAnimation;
        // Build a simple clip via serde and test round-trip
        let body_anim = BodyAnimation::new("body_0");
        let clip = AnimationClip::new("test_clip", vec![body_anim]);
        let json = serde_json::to_string(&clip).expect("serialize failed");
        let mut player = PyAnimationPlayer::new_from_json(&json).expect("from_json failed");
        assert!(!player.is_playing());
        player.play();
        assert!(player.is_playing());
    }
}
