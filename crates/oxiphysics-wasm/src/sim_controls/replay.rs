// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Replay buffer and frame types for recording and playing back simulation state.

use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

use crate::wasm_helpers::err_to_jsvalue;

// ---------------------------------------------------------------------------
// ReplayFrame
// ---------------------------------------------------------------------------

/// A single recorded physics frame.
///
/// Exposed to JavaScript as `SimReplayFrame` to avoid clashing with the replay
/// frame type defined in [`crate::debug_tools`].
#[wasm_bindgen(js_name = "SimReplayFrame")]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReplayFrame {
    /// Simulation time at this frame (s).
    pub time: f64,
    /// Flat array of body positions: \[x0,y0,z0, x1,y1,z1, ...\].
    #[wasm_bindgen(skip)]
    pub positions: Vec<f64>,
    /// Flat array of body linear velocities.
    #[wasm_bindgen(skip)]
    pub velocities: Vec<f64>,
    /// Flat array of body quaternions \[qx,qy,qz,qw, ...\].
    #[wasm_bindgen(skip)]
    pub orientations: Vec<f64>,
    /// Frame index (skipped in JS — use `frame_index_f64()`).
    #[wasm_bindgen(skip)]
    pub frame_index: u64,
}

impl ReplayFrame {
    /// Create a new empty frame at the given time.
    pub fn new(time: f64, frame_index: u64) -> Self {
        ReplayFrame {
            time,
            positions: Vec::new(),
            velocities: Vec::new(),
            orientations: Vec::new(),
            frame_index,
        }
    }

    /// Number of bodies encoded in this frame.
    pub fn body_count(&self) -> usize {
        self.positions.len() / 3
    }
}

#[wasm_bindgen(js_class = "SimReplayFrame")]
impl ReplayFrame {
    /// Construct an empty frame at simulation time `time`.
    #[wasm_bindgen(constructor)]
    pub fn new_js(time: f64, frame_index_f64: f64) -> ReplayFrame {
        ReplayFrame::new(time, frame_index_f64 as u64)
    }

    /// Frame index as a JS-friendly `f64`.
    #[wasm_bindgen(js_name = "frame_index_f64")]
    pub fn frame_index_f64(&self) -> f64 {
        self.frame_index as f64
    }

    /// Body count derived from the position array length.
    #[wasm_bindgen(js_name = "body_count")]
    pub fn body_count_js(&self) -> u32 {
        self.body_count() as u32
    }

    /// Return a clone of the position array.
    #[wasm_bindgen(js_name = "positions")]
    pub fn positions_js(&self) -> Vec<f64> {
        self.positions.clone()
    }

    /// Return a clone of the velocity array.
    #[wasm_bindgen(js_name = "velocities")]
    pub fn velocities_js(&self) -> Vec<f64> {
        self.velocities.clone()
    }

    /// Return a clone of the orientation array.
    #[wasm_bindgen(js_name = "orientations")]
    pub fn orientations_js(&self) -> Vec<f64> {
        self.orientations.clone()
    }

    /// Replace the position array with the provided flat `[x,y,z,...]` data.
    #[wasm_bindgen(js_name = "set_positions")]
    pub fn set_positions_js(&mut self, positions: Vec<f64>) {
        self.positions = positions;
    }

    /// Replace the velocity array with the provided flat `[vx,vy,vz,...]` data.
    #[wasm_bindgen(js_name = "set_velocities")]
    pub fn set_velocities_js(&mut self, velocities: Vec<f64>) {
        self.velocities = velocities;
    }

    /// Replace the orientation array with the provided flat
    /// `[qx,qy,qz,qw,...]` data.
    #[wasm_bindgen(js_name = "set_orientations")]
    pub fn set_orientations_js(&mut self, orientations: Vec<f64>) {
        self.orientations = orientations;
    }
}

// ---------------------------------------------------------------------------
// ReplayBuffer
// ---------------------------------------------------------------------------

/// Ring-buffer for recording and playing back simulation frames.
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReplayBuffer {
    /// All recorded frames in chronological order.
    #[wasm_bindgen(skip)]
    pub frames: Vec<ReplayFrame>,
    /// Maximum number of frames to retain (skipped in JS — use `capacity_js`).
    #[wasm_bindgen(skip)]
    pub capacity: usize,
    /// Current playback cursor index (skipped in JS — use `cursor_js`).
    #[wasm_bindgen(skip)]
    pub cursor: usize,
    /// Whether the buffer is currently in record mode.
    pub recording: bool,
}

impl ReplayBuffer {
    /// Create a new `ReplayBuffer` with the given capacity.
    pub fn new(capacity: usize) -> Self {
        ReplayBuffer {
            frames: Vec::with_capacity(capacity),
            capacity,
            cursor: 0,
            recording: false,
        }
    }

    /// Record a new frame. If at capacity the oldest frame is discarded.
    pub fn record_state(&mut self, frame: ReplayFrame) {
        if self.frames.len() >= self.capacity {
            self.frames.remove(0);
        }
        self.frames.push(frame);
    }

    /// Return the frame at the current playback cursor, or `None`.
    pub fn playback_frame(&self) -> Option<&ReplayFrame> {
        self.frames.get(self.cursor)
    }

    /// Seek playback cursor to the frame nearest to `time`.
    pub fn seek_to_time(&mut self, time: f64) -> usize {
        if self.frames.is_empty() {
            self.cursor = 0;
            return 0;
        }
        let mut best = 0;
        let mut best_diff = f64::MAX;
        for (i, f) in self.frames.iter().enumerate() {
            let diff = (f.time - time).abs();
            if diff < best_diff {
                best_diff = diff;
                best = i;
            }
        }
        self.cursor = best;
        best
    }

    /// Advance playback cursor by one frame; returns `false` if at end.
    pub fn advance_cursor(&mut self) -> bool {
        if self.cursor + 1 < self.frames.len() {
            self.cursor += 1;
            true
        } else {
            false
        }
    }

    /// Rewind playback cursor by one frame; returns `false` if at start.
    pub fn rewind_cursor(&mut self) -> bool {
        if self.cursor > 0 {
            self.cursor -= 1;
            true
        } else {
            false
        }
    }

    /// Export all frames to a JSON string.
    pub fn export_frames(&self) -> String {
        serde_json::to_string(self).unwrap_or_default()
    }

    /// Import frames from a JSON string, replacing current buffer.
    pub fn import_frames(&mut self, json: &str) -> Result<(), String> {
        match serde_json::from_str::<ReplayBuffer>(json) {
            Ok(buf) => {
                *self = buf;
                Ok(())
            }
            Err(e) => Err(e.to_string()),
        }
    }

    /// Clear all frames.
    pub fn clear(&mut self) {
        self.frames.clear();
        self.cursor = 0;
    }

    /// Total number of recorded frames.
    pub fn len(&self) -> usize {
        self.frames.len()
    }

    /// Returns `true` if no frames have been recorded.
    pub fn is_empty(&self) -> bool {
        self.frames.is_empty()
    }
}

#[wasm_bindgen]
impl ReplayBuffer {
    /// Construct a new buffer with the requested capacity.
    #[wasm_bindgen(constructor)]
    pub fn new_js(capacity: u32) -> ReplayBuffer {
        ReplayBuffer::new(capacity as usize)
    }

    /// Append a frame to the buffer (oldest is evicted when at capacity).
    #[wasm_bindgen(js_name = "record_state")]
    pub fn record_state_js(&mut self, frame: &ReplayFrame) {
        self.record_state(frame.clone());
    }

    /// Capacity of the ring buffer.
    #[wasm_bindgen(js_name = "capacity")]
    pub fn capacity_js(&self) -> u32 {
        self.capacity as u32
    }

    /// Current playback cursor index.
    #[wasm_bindgen(js_name = "cursor")]
    pub fn cursor_js(&self) -> u32 {
        self.cursor as u32
    }

    /// Number of recorded frames.
    #[wasm_bindgen(js_name = "len")]
    pub fn len_js(&self) -> u32 {
        self.frames.len() as u32
    }

    /// Returns `true` if no frames have been recorded.
    #[wasm_bindgen(js_name = "is_empty")]
    pub fn is_empty_js(&self) -> bool {
        self.is_empty()
    }

    /// Return a clone of the frame at the current playback cursor, or `None`.
    #[wasm_bindgen(js_name = "playback_frame")]
    pub fn playback_frame_js(&self) -> Option<ReplayFrame> {
        self.playback_frame().cloned()
    }

    /// Seek to the frame nearest `time`; returns the new cursor index.
    #[wasm_bindgen(js_name = "seek_to_time")]
    pub fn seek_to_time_js(&mut self, time: f64) -> u32 {
        self.seek_to_time(time) as u32
    }

    /// Advance the cursor by one frame; returns `false` if at end.
    #[wasm_bindgen(js_name = "advance_cursor")]
    pub fn advance_cursor_js(&mut self) -> bool {
        self.advance_cursor()
    }

    /// Rewind the cursor by one frame; returns `false` if at start.
    #[wasm_bindgen(js_name = "rewind_cursor")]
    pub fn rewind_cursor_js(&mut self) -> bool {
        self.rewind_cursor()
    }

    /// Clear all recorded frames and reset the cursor.
    #[wasm_bindgen(js_name = "clear")]
    pub fn clear_js(&mut self) {
        self.clear();
    }

    /// Export all frames to a JSON string.
    #[wasm_bindgen(js_name = "export_frames")]
    pub fn export_frames_js(&self) -> String {
        self.export_frames()
    }

    /// Import frames from a JSON string, replacing the current buffer.
    #[wasm_bindgen(js_name = "import_frames")]
    pub fn import_frames_js(&mut self, json: &str) -> Result<(), JsValue> {
        self.import_frames(json).map_err(err_to_jsvalue)
    }
}
