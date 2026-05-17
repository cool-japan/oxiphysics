// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Runtime debug visualisation configuration.
//!
//! Provides [`WasmDebugConfig`], a set of toggle flags and colour overrides
//! controlling which overlays the renderer draws each frame.

use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

use super::colors::Rgba;
use crate::wasm_helpers::to_js_value;

/// Runtime debug visualisation configuration.
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmDebugConfig {
    /// Draw axis-aligned bounding boxes.
    pub show_aabbs: bool,
    /// Draw velocity vectors.
    pub show_velocities: bool,
    /// Draw contact points.
    pub show_contacts: bool,
    /// Draw joint anchors and axes.
    pub show_joints: bool,
    /// Colour for AABB outlines.
    #[wasm_bindgen(skip)]
    pub aabb_color: Rgba,
    /// Colour for velocity arrows.
    #[wasm_bindgen(skip)]
    pub velocity_color: Rgba,
    /// Colour for contact normals.
    #[wasm_bindgen(skip)]
    pub contact_color: Rgba,
    /// Scale factor for velocity arrows.
    pub velocity_scale: f64,
}

impl Default for WasmDebugConfig {
    fn default() -> Self {
        Self {
            show_aabbs: false,
            show_velocities: false,
            show_contacts: false,
            show_joints: false,
            aabb_color: Rgba::green(),
            velocity_color: Rgba::yellow(),
            contact_color: Rgba::red(),
            velocity_scale: 0.1,
        }
    }
}

impl WasmDebugConfig {
    /// Enable all visual overlays.
    pub fn all_enabled() -> Self {
        Self {
            show_aabbs: true,
            show_velocities: true,
            show_contacts: true,
            show_joints: true,
            ..Default::default()
        }
    }
}

#[wasm_bindgen]
impl WasmDebugConfig {
    /// Default configuration (all overlays disabled, sensible colours).
    #[wasm_bindgen(constructor)]
    pub fn wasm_new() -> WasmDebugConfig {
        WasmDebugConfig::default()
    }

    /// Construct a configuration with every overlay enabled.
    #[wasm_bindgen(js_name = "all_enabled")]
    pub fn all_enabled_js() -> WasmDebugConfig {
        WasmDebugConfig::all_enabled()
    }

    /// AABB outline colour.
    #[wasm_bindgen(getter, js_name = "aabbColor")]
    pub fn aabb_color_js(&self) -> Rgba {
        self.aabb_color
    }

    /// Set AABB outline colour.
    #[wasm_bindgen(setter, js_name = "aabbColor")]
    pub fn set_aabb_color_js(&mut self, color: Rgba) {
        self.aabb_color = color;
    }

    /// Velocity arrow colour.
    #[wasm_bindgen(getter, js_name = "velocityColor")]
    pub fn velocity_color_js(&self) -> Rgba {
        self.velocity_color
    }

    /// Set velocity arrow colour.
    #[wasm_bindgen(setter, js_name = "velocityColor")]
    pub fn set_velocity_color_js(&mut self, color: Rgba) {
        self.velocity_color = color;
    }

    /// Contact normal colour.
    #[wasm_bindgen(getter, js_name = "contactColor")]
    pub fn contact_color_js(&self) -> Rgba {
        self.contact_color
    }

    /// Set contact normal colour.
    #[wasm_bindgen(setter, js_name = "contactColor")]
    pub fn set_contact_color_js(&mut self, color: Rgba) {
        self.contact_color = color;
    }

    /// Serialise to a `JsValue` object.
    #[wasm_bindgen(js_name = "to_js_value")]
    pub fn to_js_value_js(&self) -> Result<JsValue, JsValue> {
        to_js_value(self)
    }
}
