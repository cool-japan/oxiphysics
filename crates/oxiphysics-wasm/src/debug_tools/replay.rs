// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Replay controller for recorded simulation frames.

use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

use crate::wasm_helpers::to_js_value;

/// A recorded frame in the replay buffer.
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReplayFrame {
    /// Frame index.
    pub index: usize,
    /// Simulation time (s).
    pub time: f64,
    /// Serialised state (body positions as JSON).
    #[wasm_bindgen(skip)]
    pub state_json: String,
}

#[wasm_bindgen]
impl ReplayFrame {
    /// Serialised JSON snapshot of the recorded state.
    #[wasm_bindgen(getter, js_name = "state_json")]
    pub fn state_json_js(&self) -> String {
        self.state_json.clone()
    }

    /// Serialise the frame to a `JsValue` object.
    #[wasm_bindgen(js_name = "to_js_value")]
    pub fn to_js_value_js(&self) -> Result<JsValue, JsValue> {
        to_js_value(self)
    }
}

/// Controls replay of a recorded simulation.
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct WasmReplayController {
    /// Recorded frames.
    #[wasm_bindgen(skip)]
    pub frames: Vec<ReplayFrame>,
    /// Current playback position.
    pub current_frame: usize,
    /// Whether playback is active.
    pub playing: bool,
}

impl WasmReplayController {
    /// Create a new controller.
    pub fn new() -> Self {
        Self::default()
    }

    /// Record a frame.
    pub fn record_frame(&mut self, time: f64, state_json: impl Into<String>) {
        let index = self.frames.len();
        self.frames.push(ReplayFrame {
            index,
            time,
            state_json: state_json.into(),
        });
    }

    /// Step forward one frame.
    pub fn step_forward(&mut self) -> Option<&ReplayFrame> {
        if self.current_frame + 1 < self.frames.len() {
            self.current_frame += 1;
        }
        self.frames.get(self.current_frame)
    }

    /// Step backward one frame.
    pub fn step_backward(&mut self) -> Option<&ReplayFrame> {
        if self.current_frame > 0 {
            self.current_frame -= 1;
        }
        self.frames.get(self.current_frame)
    }

    /// Start playback.
    pub fn play(&mut self) {
        self.playing = true;
    }

    /// Pause playback.
    pub fn pause(&mut self) {
        self.playing = false;
    }

    /// Seek to a specific frame index.
    pub fn seek_to_frame(&mut self, frame: usize) -> Option<&ReplayFrame> {
        self.current_frame = frame.min(self.frames.len().saturating_sub(1));
        self.frames.get(self.current_frame)
    }

    /// Current frame count.
    pub fn total_frames(&self) -> usize {
        self.frames.len()
    }

    /// Current playback frame.
    pub fn current(&self) -> Option<&ReplayFrame> {
        self.frames.get(self.current_frame)
    }
}

#[wasm_bindgen]
impl WasmReplayController {
    /// Create a new controller (JS-compatible constructor).
    #[wasm_bindgen(constructor)]
    pub fn wasm_new() -> WasmReplayController {
        WasmReplayController::new()
    }

    /// Record a frame at simulation time `time` with the supplied JSON state.
    #[wasm_bindgen(js_name = "record_frame")]
    pub fn record_frame_js(&mut self, time: f64, state_json: String) {
        self.record_frame(time, state_json);
    }

    /// Step forward one frame and return the new current frame index.
    #[wasm_bindgen(js_name = "step_forward")]
    pub fn step_forward_js(&mut self) -> usize {
        self.step_forward();
        self.current_frame
    }

    /// Step backward one frame and return the new current frame index.
    #[wasm_bindgen(js_name = "step_backward")]
    pub fn step_backward_js(&mut self) -> usize {
        self.step_backward();
        self.current_frame
    }

    /// Start playback.
    #[wasm_bindgen(js_name = "play")]
    pub fn play_js(&mut self) {
        self.play();
    }

    /// Pause playback.
    #[wasm_bindgen(js_name = "pause")]
    pub fn pause_js(&mut self) {
        self.pause();
    }

    /// Seek to `frame` and return the resulting current frame index.
    #[wasm_bindgen(js_name = "seek_to_frame")]
    pub fn seek_to_frame_js(&mut self, frame: usize) -> usize {
        self.seek_to_frame(frame);
        self.current_frame
    }

    /// Total number of recorded frames.
    #[wasm_bindgen(js_name = "total_frames")]
    pub fn total_frames_js(&self) -> usize {
        self.total_frames()
    }

    /// True if a frame at the current playback position exists.
    #[wasm_bindgen(js_name = "has_current")]
    pub fn has_current_js(&self) -> bool {
        self.current().is_some()
    }

    /// Get the current frame, cloned. Returns `None` (undefined in JS) if
    /// the buffer is empty.
    #[wasm_bindgen(js_name = "get_current")]
    pub fn get_current_js(&self) -> Option<ReplayFrame> {
        self.current().cloned()
    }

    /// Get a copy of the frame at index `idx`, or `None` if out of range.
    #[wasm_bindgen(js_name = "get_frame")]
    pub fn get_frame_js(&self, idx: usize) -> Option<ReplayFrame> {
        self.frames.get(idx).cloned()
    }

    /// Serialise the controller state to a `JsValue` object.
    #[wasm_bindgen(js_name = "to_js_value")]
    pub fn to_js_value_js(&self) -> Result<JsValue, JsValue> {
        to_js_value(self)
    }
}
