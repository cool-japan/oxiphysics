// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Debug visualisation overlay configuration.

use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

use crate::wasm_helpers::err_to_jsvalue;

// ---------------------------------------------------------------------------
// Color4 / DebugDrawConfig
// ---------------------------------------------------------------------------

/// A packed RGBA colour as \[r, g, b, a\] each in \[0, 1\].
pub type Color4 = [f32; 4];

/// Configuration for debug visualisation overlays.
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DebugDrawConfig {
    /// Draw axis-aligned bounding boxes.
    pub show_aabbs: bool,
    /// Draw contact points and normals.
    pub show_contacts: bool,
    /// Draw joint/constraint anchors and axes.
    pub show_joints: bool,
    /// Draw velocity arrows on bodies.
    pub show_velocities: bool,
    /// Draw force arrows on bodies.
    pub show_forces: bool,
    /// Colour used for contact visualisation (skipped — use `contact_color()`).
    #[wasm_bindgen(skip)]
    pub contact_color: Color4,
    /// Colour used for AABB visualisation (skipped — use `aabb_color()`).
    #[wasm_bindgen(skip)]
    pub aabb_color: Color4,
    /// Colour used for joint visualisation (skipped — use `joint_color()`).
    #[wasm_bindgen(skip)]
    pub joint_color: Color4,
    /// Colour used for velocity arrows (skipped — use `velocity_color()`).
    #[wasm_bindgen(skip)]
    pub velocity_color: Color4,
    /// Scale factor for arrow lengths.
    pub arrow_scale: f32,
    /// Whether to show island IDs as text.
    pub show_island_ids: bool,
    /// Whether to show body sleep state.
    pub show_sleep_state: bool,
}

impl Default for DebugDrawConfig {
    fn default() -> Self {
        DebugDrawConfig {
            show_aabbs: false,
            show_contacts: true,
            show_joints: true,
            show_velocities: false,
            show_forces: false,
            contact_color: [1.0, 0.0, 0.0, 1.0],
            aabb_color: [0.0, 1.0, 0.0, 0.5],
            joint_color: [0.0, 0.5, 1.0, 1.0],
            velocity_color: [1.0, 1.0, 0.0, 1.0],
            arrow_scale: 0.1,
            show_island_ids: false,
            show_sleep_state: false,
        }
    }
}

impl DebugDrawConfig {
    /// Create a config with all overlays disabled.
    pub fn disabled() -> Self {
        DebugDrawConfig {
            show_aabbs: false,
            show_contacts: false,
            show_joints: false,
            show_velocities: false,
            show_forces: false,
            ..Default::default()
        }
    }

    /// Enable all overlays.
    pub fn all_enabled() -> Self {
        DebugDrawConfig {
            show_aabbs: true,
            show_contacts: true,
            show_joints: true,
            show_velocities: true,
            show_forces: true,
            ..Default::default()
        }
    }

    /// Returns `true` if any overlay is enabled.
    pub fn any_enabled(&self) -> bool {
        self.show_aabbs
            || self.show_contacts
            || self.show_joints
            || self.show_velocities
            || self.show_forces
    }
}

#[wasm_bindgen]
impl DebugDrawConfig {
    /// Construct a config from defaults.
    #[wasm_bindgen(constructor)]
    pub fn new_js() -> DebugDrawConfig {
        DebugDrawConfig::default()
    }

    /// Construct a config with all overlays disabled.
    #[wasm_bindgen(js_name = "disabled")]
    pub fn disabled_js() -> DebugDrawConfig {
        DebugDrawConfig::disabled()
    }

    /// Construct a config with all overlays enabled.
    #[wasm_bindgen(js_name = "all_enabled")]
    pub fn all_enabled_js() -> DebugDrawConfig {
        DebugDrawConfig::all_enabled()
    }

    /// JS-friendly wrapper for [`DebugDrawConfig::any_enabled`].
    #[wasm_bindgen(js_name = "any_enabled")]
    pub fn any_enabled_js(&self) -> bool {
        self.any_enabled()
    }

    /// Return the contact-marker colour as a flat `[r, g, b, a]` `Float32Array`.
    #[wasm_bindgen(js_name = "contact_color")]
    pub fn contact_color_js(&self) -> Vec<f32> {
        self.contact_color.to_vec()
    }

    /// Return the AABB-overlay colour as a flat `[r, g, b, a]` array.
    #[wasm_bindgen(js_name = "aabb_color")]
    pub fn aabb_color_js(&self) -> Vec<f32> {
        self.aabb_color.to_vec()
    }

    /// Return the joint-overlay colour as a flat `[r, g, b, a]` array.
    #[wasm_bindgen(js_name = "joint_color")]
    pub fn joint_color_js(&self) -> Vec<f32> {
        self.joint_color.to_vec()
    }

    /// Return the velocity-arrow colour as a flat `[r, g, b, a]` array.
    #[wasm_bindgen(js_name = "velocity_color")]
    pub fn velocity_color_js(&self) -> Vec<f32> {
        self.velocity_color.to_vec()
    }

    /// Set the contact-marker colour.
    #[wasm_bindgen(js_name = "set_contact_color")]
    pub fn set_contact_color_js(&mut self, r: f32, g: f32, b: f32, a: f32) {
        self.contact_color = [r, g, b, a];
    }

    /// Set the AABB-overlay colour.
    #[wasm_bindgen(js_name = "set_aabb_color")]
    pub fn set_aabb_color_js(&mut self, r: f32, g: f32, b: f32, a: f32) {
        self.aabb_color = [r, g, b, a];
    }

    /// Set the joint-overlay colour.
    #[wasm_bindgen(js_name = "set_joint_color")]
    pub fn set_joint_color_js(&mut self, r: f32, g: f32, b: f32, a: f32) {
        self.joint_color = [r, g, b, a];
    }

    /// Set the velocity-arrow colour.
    #[wasm_bindgen(js_name = "set_velocity_color")]
    pub fn set_velocity_color_js(&mut self, r: f32, g: f32, b: f32, a: f32) {
        self.velocity_color = [r, g, b, a];
    }

    /// Serialise to a JSON string.
    #[wasm_bindgen(js_name = "to_json")]
    pub fn to_json_js(&self) -> Result<String, JsValue> {
        serde_json::to_string(self).map_err(err_to_jsvalue)
    }
}
