// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Debug-draw primitives.
//!
//! [`WasmDrawCall`] tags one of five primitive shapes (line, sphere, AABB,
//! arrow, text). [`WasmDebugDrawList`] accumulates them; [`WasmDebugDraw`] is
//! a fluent helper that lets callers assemble a list across multiple frames.

use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

use super::colors::Rgba;
use crate::wasm_helpers::{err_to_jsvalue, to_js_value};

/// A single debug draw call.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum WasmDrawCall {
    /// Draw a line segment.
    Line {
        /// Start point.
        start: [f64; 3],
        /// End point.
        end: [f64; 3],
        /// Line colour.
        color: Rgba,
    },
    /// Draw a wire-frame sphere.
    Sphere {
        /// Centre of the sphere.
        center: [f64; 3],
        /// Sphere radius.
        radius: f64,
        /// Colour.
        color: Rgba,
    },
    /// Draw a wire-frame axis-aligned box.
    Box {
        /// Minimum corner.
        min: [f64; 3],
        /// Maximum corner.
        max: [f64; 3],
        /// Colour.
        color: Rgba,
    },
    /// Draw an arrow (line + arrowhead).
    Arrow {
        /// Tail of the arrow.
        from: [f64; 3],
        /// Tip of the arrow.
        to: [f64; 3],
        /// Colour.
        color: Rgba,
    },
    /// Draw a text label.
    Text {
        /// World-space position.
        position: [f64; 3],
        /// Text to display.
        label: String,
        /// Colour.
        color: Rgba,
    },
}

/// A list of debug draw calls for a single frame.
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct WasmDebugDrawList {
    /// All draw calls.
    #[wasm_bindgen(skip)]
    pub calls: Vec<WasmDrawCall>,
}

impl WasmDebugDrawList {
    /// Create an empty draw list.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a draw call.
    pub fn push(&mut self, call: WasmDrawCall) {
        self.calls.push(call);
    }

    /// Number of draw calls.
    pub fn len(&self) -> usize {
        self.calls.len()
    }

    /// True if no draw calls.
    pub fn is_empty(&self) -> bool {
        self.calls.is_empty()
    }

    /// Clear all draw calls.
    pub fn clear(&mut self) {
        self.calls.clear();
    }

    /// Serialise to JSON for JavaScript consumption.
    pub fn to_json(&self) -> Result<String, serde_json::Error> {
        serde_json::to_string(self)
    }
}

#[wasm_bindgen]
impl WasmDebugDrawList {
    /// Create an empty draw list (JS-compatible constructor).
    #[wasm_bindgen(constructor)]
    pub fn wasm_new() -> WasmDebugDrawList {
        WasmDebugDrawList::new()
    }

    /// Number of draw calls (JS-exposed).
    #[wasm_bindgen(js_name = "len")]
    pub fn len_js(&self) -> usize {
        self.len()
    }

    /// Whether the list is empty (JS-exposed).
    #[wasm_bindgen(js_name = "is_empty")]
    pub fn is_empty_js(&self) -> bool {
        self.is_empty()
    }

    /// Clear all draw calls (JS-exposed).
    #[wasm_bindgen(js_name = "clear")]
    pub fn clear_js(&mut self) {
        self.clear();
    }

    /// Serialise the list to a JSON string.
    #[wasm_bindgen(js_name = "to_json")]
    pub fn to_json_js(&self) -> Result<String, JsValue> {
        self.to_json().map_err(err_to_jsvalue)
    }

    /// Serialise the list to a `JsValue` object.
    #[wasm_bindgen(js_name = "to_js_value")]
    pub fn to_js_value_js(&self) -> Result<JsValue, JsValue> {
        to_js_value(self)
    }
}

/// Debug draw helper with fluent API.
#[wasm_bindgen]
#[derive(Debug, Clone, Default)]
pub struct WasmDebugDraw {
    /// Accumulated draw list.
    #[wasm_bindgen(skip)]
    pub list: WasmDebugDrawList,
}

impl WasmDebugDraw {
    /// Create a new draw helper.
    pub fn new() -> Self {
        Self::default()
    }

    /// Draw a line segment.
    pub fn draw_line(&mut self, start: [f64; 3], end: [f64; 3], color: Rgba) {
        self.list.push(WasmDrawCall::Line { start, end, color });
    }

    /// Draw a sphere.
    pub fn draw_sphere(&mut self, center: [f64; 3], radius: f64, color: Rgba) {
        self.list.push(WasmDrawCall::Sphere {
            center,
            radius,
            color,
        });
    }

    /// Draw an AABB box.
    pub fn draw_box(&mut self, min: [f64; 3], max: [f64; 3], color: Rgba) {
        self.list.push(WasmDrawCall::Box { min, max, color });
    }

    /// Draw an arrow.
    pub fn draw_arrow(&mut self, from: [f64; 3], to: [f64; 3], color: Rgba) {
        self.list.push(WasmDrawCall::Arrow { from, to, color });
    }

    /// Draw a text label.
    pub fn draw_text(&mut self, position: [f64; 3], label: impl Into<String>, color: Rgba) {
        self.list.push(WasmDrawCall::Text {
            position,
            label: label.into(),
            color,
        });
    }

    /// Flush (clone) the draw list and clear the internal one.
    pub fn flush(&mut self) -> WasmDebugDrawList {
        let out = self.list.clone();
        self.list.clear();
        out
    }
}

#[wasm_bindgen]
impl WasmDebugDraw {
    /// Create a new draw helper (JS-compatible constructor).
    #[wasm_bindgen(constructor)]
    pub fn wasm_new() -> WasmDebugDraw {
        WasmDebugDraw::new()
    }

    /// Draw a line segment between `(sx, sy, sz)` and `(ex, ey, ez)`.
    #[wasm_bindgen(js_name = "draw_line")]
    pub fn draw_line_js(
        &mut self,
        sx: f64,
        sy: f64,
        sz: f64,
        ex: f64,
        ey: f64,
        ez: f64,
        color: &Rgba,
    ) {
        self.draw_line([sx, sy, sz], [ex, ey, ez], *color);
    }

    /// Draw a wire-frame sphere centred at `(cx, cy, cz)`.
    #[wasm_bindgen(js_name = "draw_sphere")]
    pub fn draw_sphere_js(&mut self, cx: f64, cy: f64, cz: f64, radius: f64, color: &Rgba) {
        self.draw_sphere([cx, cy, cz], radius, *color);
    }

    /// Draw an AABB from `(min_x, min_y, min_z)` to `(max_x, max_y, max_z)`.
    #[wasm_bindgen(js_name = "draw_box")]
    pub fn draw_box_js(
        &mut self,
        min_x: f64,
        min_y: f64,
        min_z: f64,
        max_x: f64,
        max_y: f64,
        max_z: f64,
        color: &Rgba,
    ) {
        self.draw_box([min_x, min_y, min_z], [max_x, max_y, max_z], *color);
    }

    /// Draw an arrow from `(fx, fy, fz)` to `(tx, ty, tz)`.
    #[wasm_bindgen(js_name = "draw_arrow")]
    pub fn draw_arrow_js(
        &mut self,
        fx: f64,
        fy: f64,
        fz: f64,
        tx: f64,
        ty: f64,
        tz: f64,
        color: &Rgba,
    ) {
        self.draw_arrow([fx, fy, fz], [tx, ty, tz], *color);
    }

    /// Draw a text label at `(x, y, z)`.
    #[wasm_bindgen(js_name = "draw_text")]
    pub fn draw_text_js(&mut self, x: f64, y: f64, z: f64, label: String, color: &Rgba) {
        self.draw_text([x, y, z], label, *color);
    }

    /// Number of accumulated draw calls.
    #[wasm_bindgen(js_name = "len")]
    pub fn len_js(&self) -> usize {
        self.list.len()
    }

    /// Whether there are no accumulated draw calls.
    #[wasm_bindgen(js_name = "is_empty")]
    pub fn is_empty_js(&self) -> bool {
        self.list.is_empty()
    }

    /// Flush the accumulated draw list and clear the internal buffer.
    #[wasm_bindgen(js_name = "flush")]
    pub fn flush_js(&mut self) -> WasmDebugDrawList {
        self.flush()
    }

    /// Serialise the current accumulated list to a `JsValue` object.
    #[wasm_bindgen(js_name = "list_to_js_value")]
    pub fn list_to_js_value_js(&self) -> Result<JsValue, JsValue> {
        to_js_value(&self.list)
    }
}
