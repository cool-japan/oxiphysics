// Copyright (C) 2026 COOLJAPAN OU (Team Kitasan)
// SPDX-License-Identifier: Apache-2.0

//! WebAssembly bridge for keyframe animation tracks.
//!
//! Exposes a JSON-oriented surface around keyframe animation clips and
//! an animation player suitable for use across the WASM boundary.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

use crate::wasm_helpers::{err_to_jsvalue, to_js_value};

// ---------------------------------------------------------------------------
// EaseKind
// ---------------------------------------------------------------------------

/// Easing function applied when interpolating from a keyframe toward the next.
#[wasm_bindgen]
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum WasmEaseKind {
    #[default]
    Linear,
    EaseIn,
    EaseOut,
    EaseInOut,
    Step,
}

impl WasmEaseKind {
    fn apply(self, t: f64) -> f64 {
        let t = t.clamp(0.0, 1.0);
        match self {
            WasmEaseKind::Linear => t,
            WasmEaseKind::EaseIn => t * t * t,
            WasmEaseKind::EaseOut => {
                let u = 1.0 - t;
                1.0 - u * u * u
            }
            WasmEaseKind::EaseInOut => {
                if t < 0.5 {
                    4.0 * t * t * t
                } else {
                    let u = -2.0 * t + 2.0;
                    1.0 - u * u * u / 2.0
                }
            }
            WasmEaseKind::Step => {
                if t < 1.0 {
                    0.0
                } else {
                    1.0
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Vec3Keyframe / Vec3Track
// ---------------------------------------------------------------------------

#[derive(Clone, Debug, Serialize, Deserialize)]
struct Vec3Keyframe {
    time: f64,
    value: [f64; 3],
    ease: WasmEaseKind,
}

/// A keyframe track that samples a `[f64; 3]` value over time.
#[wasm_bindgen]
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct WasmVec3Track {
    keyframes: Vec<Vec3Keyframe>,
    looping: bool,
}

impl WasmVec3Track {
    /// Create a new empty Vec3 track.
    pub fn new(looping: bool) -> Self {
        WasmVec3Track {
            keyframes: Vec::new(),
            looping,
        }
    }

    /// Append a keyframe at `time` with the given `value` and easing.
    pub fn push_keyframe(&mut self, time: f64, value: [f64; 3], ease: WasmEaseKind) {
        self.keyframes.push(Vec3Keyframe { time, value, ease });
        self.keyframes.sort_by(|a, b| {
            a.time
                .partial_cmp(&b.time)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
    }

    /// Sample the track at time `t`.
    pub fn sample(&self, mut t: f64) -> [f64; 3] {
        if self.keyframes.is_empty() {
            return [0.0; 3];
        }
        if self.keyframes.len() == 1 {
            return self.keyframes[0].value;
        }
        if self.looping {
            let dur = self.keyframes.last().map_or(1.0, |k| k.time);
            if dur > 0.0 {
                t = t.rem_euclid(dur);
            }
        }
        let first = &self.keyframes[0];
        if t <= first.time {
            return first.value;
        }
        let last_value = match self.keyframes.last() {
            Some(k) => k.value,
            None => return [0.0; 3],
        };
        let last_time = match self.keyframes.last() {
            Some(k) => k.time,
            None => return [0.0; 3],
        };
        if t >= last_time {
            return last_value;
        }
        // Find bracketing keyframes
        for window in self.keyframes.windows(2) {
            let a = &window[0];
            let b = &window[1];
            if t >= a.time && t <= b.time {
                let span = b.time - a.time;
                let raw_t = if span < 1e-12 {
                    0.0
                } else {
                    (t - a.time) / span
                };
                let eased = a.ease.apply(raw_t);
                return [
                    a.value[0] + (b.value[0] - a.value[0]) * eased,
                    a.value[1] + (b.value[1] - a.value[1]) * eased,
                    a.value[2] + (b.value[2] - a.value[2]) * eased,
                ];
            }
        }
        last_value
    }

    /// Duration of the track (time of the last keyframe), or `0.0` if empty.
    pub fn duration(&self) -> f64 {
        self.keyframes.last().map_or(0.0, |k| k.time)
    }
}

// ---------------------------------------------------------------------------
// WasmVec3Track — wasm-bindgen impl
// ---------------------------------------------------------------------------

#[wasm_bindgen]
impl WasmVec3Track {
    /// Construct an empty Vec3 track (JS constructor).
    #[wasm_bindgen(constructor)]
    pub fn wasm_new(looping: bool) -> WasmVec3Track {
        WasmVec3Track::new(looping)
    }

    /// Append a keyframe at `time` with value `[x, y, z]` and the given easing.
    #[wasm_bindgen(js_name = "push_keyframe")]
    pub fn push_keyframe_js(&mut self, time: f64, x: f64, y: f64, z: f64, ease: WasmEaseKind) {
        self.push_keyframe(time, [x, y, z], ease);
    }

    /// Sample the track at time `t`; returns a flat `Vec<f64>` of length 3.
    #[wasm_bindgen(js_name = "sample")]
    pub fn sample_js(&self, t: f64) -> Vec<f64> {
        self.sample(t).to_vec()
    }

    /// Duration of the track (time of the last keyframe), or `0.0` if empty.
    #[wasm_bindgen(js_name = "duration")]
    pub fn duration_js(&self) -> f64 {
        self.duration()
    }

    /// Number of keyframes in the track.
    #[wasm_bindgen(js_name = "keyframe_count")]
    pub fn keyframe_count_js(&self) -> usize {
        self.keyframes.len()
    }

    /// Whether the track loops at its duration.
    #[wasm_bindgen(getter)]
    pub fn looping(&self) -> bool {
        self.looping
    }

    /// Set whether the track loops.
    #[wasm_bindgen(setter)]
    pub fn set_looping(&mut self, looping: bool) {
        self.looping = looping;
    }

    /// Serialise the track as a JSON string.
    #[wasm_bindgen(js_name = "to_json")]
    pub fn to_json_js(&self) -> Result<String, JsValue> {
        serde_json::to_string(self).map_err(err_to_jsvalue)
    }
}

// ---------------------------------------------------------------------------
// AnimBodyState — sampled output for one body
// ---------------------------------------------------------------------------

/// Sampled state for one animated body.
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmAnimBodyState {
    /// Identifier of the body (matches `BodyAnimation::body_id`).
    #[wasm_bindgen(skip)]
    pub body_id: String,
    /// Sampled position `[x, y, z]`, or `None` if no position track.
    #[wasm_bindgen(skip)]
    pub position: Option<[f64; 3]>,
}

#[wasm_bindgen]
impl WasmAnimBodyState {
    /// Body identifier (matches `WasmBodyAnimation::body_id`).
    #[wasm_bindgen(getter)]
    pub fn body_id(&self) -> String {
        self.body_id.clone()
    }

    /// Sampled position as a flat `Vec<f64>`. Returns an empty vector when
    /// no position track is present.
    #[wasm_bindgen(js_name = "get_position")]
    pub fn get_position_js(&self) -> Vec<f64> {
        match self.position {
            Some(p) => p.to_vec(),
            None => Vec::new(),
        }
    }

    /// Whether a position sample is available.
    #[wasm_bindgen(js_name = "has_position")]
    pub fn has_position_js(&self) -> bool {
        self.position.is_some()
    }

    /// Serialise the state as a JSON string.
    #[wasm_bindgen(js_name = "to_json")]
    pub fn to_json_js(&self) -> Result<String, JsValue> {
        serde_json::to_string(self).map_err(err_to_jsvalue)
    }
}

// ---------------------------------------------------------------------------
// BodyAnimation / AnimationClip
// ---------------------------------------------------------------------------

/// Animation data for a single body.
#[wasm_bindgen]
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct WasmBodyAnimation {
    #[wasm_bindgen(skip)]
    pub body_id: String,
    #[wasm_bindgen(skip)]
    pub position: Option<WasmVec3Track>,
}

impl WasmBodyAnimation {
    /// Create a new body animation with the given ID and no tracks.
    pub fn new(body_id: &str) -> Self {
        WasmBodyAnimation {
            body_id: body_id.to_string(),
            position: None,
        }
    }

    /// Attach a position track.
    pub fn with_position(mut self, track: WasmVec3Track) -> Self {
        self.position = Some(track);
        self
    }
}

#[wasm_bindgen]
impl WasmBodyAnimation {
    /// Construct a new body animation with no tracks (JS constructor).
    #[wasm_bindgen(constructor)]
    pub fn wasm_new(body_id: String) -> WasmBodyAnimation {
        WasmBodyAnimation::new(&body_id)
    }

    /// Body identifier.
    #[wasm_bindgen(getter)]
    pub fn body_id(&self) -> String {
        self.body_id.clone()
    }

    /// Whether a position track is attached.
    #[wasm_bindgen(js_name = "has_position")]
    pub fn has_position_js(&self) -> bool {
        self.position.is_some()
    }

    /// Attach a position track. Replaces any existing track.
    #[wasm_bindgen(js_name = "set_position_track")]
    pub fn set_position_track_js(&mut self, track: &WasmVec3Track) {
        self.position = Some(track.clone());
    }

    /// Clear the attached position track, if any.
    #[wasm_bindgen(js_name = "clear_position_track")]
    pub fn clear_position_track_js(&mut self) {
        self.position = None;
    }

    /// Serialise the body animation as a JSON string.
    #[wasm_bindgen(js_name = "to_json")]
    pub fn to_json_js(&self) -> Result<String, JsValue> {
        serde_json::to_string(self).map_err(err_to_jsvalue)
    }
}

/// A named collection of per-body animation tracks.
#[wasm_bindgen]
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct WasmAnimationClip {
    #[wasm_bindgen(skip)]
    pub name: String,
    #[wasm_bindgen(skip)]
    pub body_animations: Vec<WasmBodyAnimation>,
}

impl WasmAnimationClip {
    /// Create a new clip.
    pub fn new(name: &str, body_animations: Vec<WasmBodyAnimation>) -> Self {
        WasmAnimationClip {
            name: name.to_string(),
            body_animations,
        }
    }

    /// Duration of the clip — maximum duration across all body animations.
    pub fn duration(&self) -> f64 {
        self.body_animations
            .iter()
            .flat_map(|ba| ba.position.as_ref().map(|t| t.duration()))
            .fold(0.0_f64, f64::max)
    }
}

#[wasm_bindgen]
impl WasmAnimationClip {
    /// Construct an empty clip (JS constructor).
    #[wasm_bindgen(constructor)]
    pub fn wasm_new(name: String) -> WasmAnimationClip {
        WasmAnimationClip::new(&name, Vec::new())
    }

    /// Clip name.
    #[wasm_bindgen(getter)]
    pub fn name(&self) -> String {
        self.name.clone()
    }

    /// Replace the clip name.
    #[wasm_bindgen(setter)]
    pub fn set_name(&mut self, name: String) {
        self.name = name;
    }

    /// Number of body animations in the clip.
    #[wasm_bindgen(js_name = "body_count")]
    pub fn body_count_js(&self) -> usize {
        self.body_animations.len()
    }

    /// Append a body animation to the clip.
    #[wasm_bindgen(js_name = "push_body_animation")]
    pub fn push_body_animation_js(&mut self, anim: &WasmBodyAnimation) {
        self.body_animations.push(anim.clone());
    }

    /// Duration of the clip (longest constituent body animation), in seconds.
    #[wasm_bindgen(js_name = "duration")]
    pub fn duration_js(&self) -> f64 {
        self.duration()
    }

    /// Serialise the clip as a JSON string.
    #[wasm_bindgen(js_name = "to_json")]
    pub fn to_json_js(&self) -> Result<String, JsValue> {
        serde_json::to_string(self).map_err(err_to_jsvalue)
    }
}

// ---------------------------------------------------------------------------
// WasmAnimationPlayer
// ---------------------------------------------------------------------------

/// WASM wrapper for keyframe animation playback.
///
/// Drives an [`WasmAnimationClip`] forward in time and returns sampled
/// per-body states as a JSON string.
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmAnimationPlayer {
    #[wasm_bindgen(skip)]
    pub clip: WasmAnimationClip,
    pub time: f64,
    pub speed: f64,
    playing: bool,
    pub looping: bool,
}

impl WasmAnimationPlayer {
    /// Create a new animation player for the given clip.
    pub fn new(clip: WasmAnimationClip) -> Self {
        WasmAnimationPlayer {
            clip,
            time: 0.0,
            speed: 1.0,
            playing: false,
            looping: false,
        }
    }

    /// Start or resume playback.
    pub fn play(&mut self) {
        self.playing = true;
    }

    /// Pause playback (retains current time).
    pub fn pause(&mut self) {
        self.playing = false;
    }

    /// Stop playback and rewind to `t = 0`.
    pub fn stop(&mut self) {
        self.playing = false;
        self.time = 0.0;
    }

    /// Seek to an absolute time (clamped to clip duration).
    pub fn seek(&mut self, t: f64) {
        self.time = t.clamp(0.0, self.clip.duration());
    }

    /// Returns `true` if the player is actively advancing time.
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

    /// Returns `true` when time has reached or passed the clip duration.
    pub fn is_finished(&self) -> bool {
        self.time >= self.clip.duration()
    }

    /// Advance playback by `dt` seconds and return sampled states.
    pub fn update(&mut self, dt: f64) -> Vec<WasmAnimBodyState> {
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

    /// Sample all body animations at the current time without advancing.
    pub fn sample_all(&self) -> Vec<WasmAnimBodyState> {
        self.clip
            .body_animations
            .iter()
            .map(|ba| WasmAnimBodyState {
                body_id: ba.body_id.clone(),
                position: ba.position.as_ref().map(|t| t.sample(self.time)),
            })
            .collect()
    }

    /// Serialise the current sampled states as a JSON string.
    pub fn to_json(&self) -> String {
        let states = self.sample_all();
        serde_json::to_string(&states).unwrap_or_default()
    }
}

// ---------------------------------------------------------------------------
// WasmAnimationPlayer — wasm-bindgen impl
// ---------------------------------------------------------------------------

#[wasm_bindgen]
impl WasmAnimationPlayer {
    /// Construct a new animation player for the given clip (JS constructor).
    #[wasm_bindgen(constructor)]
    pub fn wasm_new(clip: &WasmAnimationClip) -> WasmAnimationPlayer {
        WasmAnimationPlayer::new(clip.clone())
    }

    /// Start or resume playback.
    #[wasm_bindgen(js_name = "play")]
    pub fn play_js(&mut self) {
        self.play();
    }

    /// Pause playback (retains current time).
    #[wasm_bindgen(js_name = "pause")]
    pub fn pause_js(&mut self) {
        self.pause();
    }

    /// Stop playback and rewind to `t = 0`.
    #[wasm_bindgen(js_name = "stop")]
    pub fn stop_js(&mut self) {
        self.stop();
    }

    /// Seek to an absolute time (clamped to clip duration).
    #[wasm_bindgen(js_name = "seek")]
    pub fn seek_js(&mut self, t: f64) {
        self.seek(t);
    }

    /// Whether the player is actively advancing time.
    #[wasm_bindgen(js_name = "is_playing")]
    pub fn is_playing_js(&self) -> bool {
        self.is_playing()
    }

    /// Playback progress as a fraction `[0, 1]`.
    #[wasm_bindgen(js_name = "progress")]
    pub fn progress_js(&self) -> f64 {
        self.progress()
    }

    /// Whether time has reached or passed the clip duration.
    #[wasm_bindgen(js_name = "is_finished")]
    pub fn is_finished_js(&self) -> bool {
        self.is_finished()
    }

    /// Advance playback by `dt` seconds and return the sampled states as a JSON string.
    #[wasm_bindgen(js_name = "update_json")]
    pub fn update_json_js(&mut self, dt: f64) -> Result<String, JsValue> {
        let states = self.update(dt);
        serde_json::to_string(&states).map_err(err_to_jsvalue)
    }

    /// Advance playback by `dt` seconds and return the sampled states as a `JsValue`.
    #[wasm_bindgen(js_name = "update_js_value")]
    pub fn update_js_value(&mut self, dt: f64) -> Result<JsValue, JsValue> {
        let states = self.update(dt);
        to_js_value(&states)
    }

    /// Sample all body animations at the current time and return a JSON string.
    #[wasm_bindgen(js_name = "sample_all_json")]
    pub fn sample_all_json_js(&self) -> Result<String, JsValue> {
        let states = self.sample_all();
        serde_json::to_string(&states).map_err(err_to_jsvalue)
    }

    /// Sample all body animations at the current time and return a `JsValue`.
    #[wasm_bindgen(js_name = "sample_all_js_value")]
    pub fn sample_all_js_value(&self) -> Result<JsValue, JsValue> {
        let states = self.sample_all();
        to_js_value(&states)
    }

    /// Serialise the current sampled states as a JSON string.
    #[wasm_bindgen(js_name = "to_json")]
    pub fn to_json_js(&self) -> String {
        self.to_json()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_animation_bridge_instantiation() {
        let clip = WasmAnimationClip::new("test_clip", vec![]);
        let player = WasmAnimationPlayer::new(clip);
        let json = player.to_json();
        assert!(!json.is_empty(), "to_json should return non-empty string");
    }

    #[test]
    fn test_animation_bridge_playback() {
        let mut track = WasmVec3Track::new(false);
        track.push_keyframe(0.0, [0.0, 0.0, 0.0], WasmEaseKind::Linear);
        track.push_keyframe(2.0, [10.0, 0.0, 0.0], WasmEaseKind::Linear);

        let ba = WasmBodyAnimation::new("obj").with_position(track);
        let clip = WasmAnimationClip::new("move", vec![ba]);
        let mut player = WasmAnimationPlayer::new(clip);

        // Stopped: update should not advance time
        player.update(1.0);
        assert!((player.time).abs() < 1e-9);

        player.play();
        let states = player.update(1.0);
        assert_eq!(states.len(), 1);
        assert!((states[0].position.expect("has position")[0] - 5.0).abs() < 1e-9);
        assert!(!player.is_finished());

        player.update(1.0);
        assert!(player.is_finished());
        assert!(!player.is_playing());
    }

    #[test]
    fn test_animation_bridge_seek_progress() {
        let mut track = WasmVec3Track::new(false);
        track.push_keyframe(0.0, [0.0; 3], WasmEaseKind::Linear);
        track.push_keyframe(4.0, [4.0, 0.0, 0.0], WasmEaseKind::Linear);

        let ba = WasmBodyAnimation::new("b").with_position(track);
        let clip = WasmAnimationClip::new("clip", vec![ba]);
        let mut player = WasmAnimationPlayer::new(clip);
        player.seek(2.0);
        assert!((player.progress() - 0.5).abs() < 1e-9);
    }
}
