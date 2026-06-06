// Copyright (C) 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! WebAssembly bridge for the kinematic capsule character controller.
//!
//! Wraps the character controller geometry and configuration types and
//! exposes a JSON-oriented surface suitable for use across the WASM boundary.

use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

use crate::wasm_helpers::{err_to_jsvalue, opt_vec3_to_js, to_js_value, unflatten_vec3s};

// ---------------------------------------------------------------------------
// Math helpers (private, inline)
// ---------------------------------------------------------------------------

#[inline]
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline]
fn project_onto_plane(v: [f64; 3], normal: [f64; 3]) -> [f64; 3] {
    let d = dot(v, normal);
    sub(v, scale(normal, d))
}

// ---------------------------------------------------------------------------
// WasmCharacterShape
// ---------------------------------------------------------------------------

/// Capsule geometry for the character controller.
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmCharacterShape {
    /// Capsule radius (metres).
    pub radius: f64,
    /// Half-height of the cylindrical segment (metres). Total height = 2*half_height + 2*radius.
    pub half_height: f64,
}

impl WasmCharacterShape {
    /// Create a new capsule shape (Rust-only).
    pub fn new(radius: f64, half_height: f64) -> Self {
        WasmCharacterShape {
            radius,
            half_height,
        }
    }
}

#[wasm_bindgen]
impl WasmCharacterShape {
    /// Create a new capsule shape (JS constructor).
    #[wasm_bindgen(constructor)]
    pub fn new_js(radius: f64, half_height: f64) -> WasmCharacterShape {
        WasmCharacterShape::new(radius, half_height)
    }

    /// Default human-sized shape used by the controller.
    #[wasm_bindgen(js_name = "default_shape")]
    pub fn default_shape_js() -> WasmCharacterShape {
        WasmCharacterShape::default()
    }
}

impl Default for WasmCharacterShape {
    fn default() -> Self {
        WasmCharacterShape {
            radius: 0.4,
            half_height: 0.9,
        }
    }
}

// ---------------------------------------------------------------------------
// WasmCharacterConfig
// ---------------------------------------------------------------------------

/// Tuning parameters for the character controller.
///
/// `up_axis` is a private `[f64; 3]` and `max_iterations` is a private
/// `usize`; both are exposed to JavaScript through accessor methods on the
/// `#[wasm_bindgen]` impl block.
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmCharacterConfig {
    /// Maximum slope angle (degrees) treated as walkable ground.
    pub max_slope_deg: f64,
    /// Maximum step height the controller can automatically climb (metres).
    pub step_offset: f64,
    /// Small inset applied to prevent tunnelling (metres).
    pub skin_width: f64,
    /// World-space up axis `[x, y, z]` — private; use the JS accessors.
    pub(crate) up_axis: [f64; 3],
    /// Maximum depenetration / slide iterations per move call —
    /// private; use the JS accessor (`u32`).
    pub(crate) max_iterations: usize,
}

impl Default for WasmCharacterConfig {
    fn default() -> Self {
        WasmCharacterConfig {
            max_slope_deg: 45.0,
            step_offset: 0.3,
            skin_width: 1e-3,
            up_axis: [0.0, 1.0, 0.0],
            max_iterations: 4,
        }
    }
}

#[wasm_bindgen]
impl WasmCharacterConfig {
    /// Default configuration suitable for a humanoid character (JS constructor).
    #[wasm_bindgen(constructor)]
    pub fn new_js() -> WasmCharacterConfig {
        WasmCharacterConfig::default()
    }

    /// Up axis as `[ux, uy, uz]` (`Float64Array`).
    #[wasm_bindgen(js_name = "get_up_axis")]
    pub fn get_up_axis_js(&self) -> Vec<f64> {
        self.up_axis.to_vec()
    }

    /// Replace the up axis.
    #[wasm_bindgen(js_name = "set_up_axis")]
    pub fn set_up_axis_js(&mut self, ux: f64, uy: f64, uz: f64) {
        self.up_axis = [ux, uy, uz];
    }

    /// Maximum slide iterations as `u32`.
    #[wasm_bindgen(js_name = "max_iterations")]
    pub fn max_iterations_js(&self) -> u32 {
        self.max_iterations as u32
    }

    /// Set the maximum slide iterations (clamped to >= 1).
    #[wasm_bindgen(js_name = "set_max_iterations")]
    pub fn set_max_iterations_js(&mut self, n: u32) {
        self.max_iterations = (n as usize).max(1);
    }
}

// ---------------------------------------------------------------------------
// WasmSweepHit
// ---------------------------------------------------------------------------

/// Result of a single sweep query used during `move_and_slide`.
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmSweepHit {
    /// Fractional hit distance along the sweep (0 = start, 1 = end).
    pub toi: f64,
    /// Hit surface normal `[x, y, z]` pointing away from the surface —
    /// private; use the JS accessor.
    pub(crate) normal: [f64; 3],
}

#[wasm_bindgen]
impl WasmSweepHit {
    /// Create a sweep hit with the given time-of-impact and surface normal.
    #[wasm_bindgen(constructor)]
    pub fn new_js(toi: f64, nx: f64, ny: f64, nz: f64) -> WasmSweepHit {
        WasmSweepHit {
            toi,
            normal: [nx, ny, nz],
        }
    }

    /// Hit surface normal as `[nx, ny, nz]` (`Float64Array`).
    #[wasm_bindgen(js_name = "get_normal")]
    pub fn get_normal_js(&self) -> Vec<f64> {
        self.normal.to_vec()
    }

    /// Replace the hit surface normal.
    #[wasm_bindgen(js_name = "set_normal")]
    pub fn set_normal_js(&mut self, nx: f64, ny: f64, nz: f64) {
        self.normal = [nx, ny, nz];
    }
}

// ---------------------------------------------------------------------------
// WasmCharacterMove — output of move_and_slide_simple
// ---------------------------------------------------------------------------

/// Output of a [`WasmCharacterController::move_and_slide_simple`] call.
///
/// The `[f64; 3]` and `Option<[f64; 3]>` fields are private; use the
/// `#[wasm_bindgen]` accessors to read them from JavaScript.
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmCharacterMove {
    /// Effective translation applied this step — private; use the JS accessor.
    pub(crate) translation: [f64; 3],
    /// `true` if the controller is currently standing on a walkable surface.
    pub is_grounded: bool,
    /// Surface normal of the ground below — private; use the JS accessor.
    pub(crate) ground_normal: Option<[f64; 3]>,
}

#[wasm_bindgen]
impl WasmCharacterMove {
    /// Effective translation as `[x, y, z]` (`Float64Array`).
    #[wasm_bindgen(js_name = "get_translation")]
    pub fn get_translation_js(&self) -> Vec<f64> {
        self.translation.to_vec()
    }

    /// Ground normal as `[nx, ny, nz]`, or an empty array when not grounded.
    #[wasm_bindgen(js_name = "get_ground_normal")]
    pub fn get_ground_normal_js(&self) -> Vec<f64> {
        opt_vec3_to_js(self.ground_normal)
    }

    /// `true` if the move ended grounded.
    #[wasm_bindgen(js_name = "is_grounded")]
    pub fn is_grounded_js(&self) -> bool {
        self.is_grounded
    }
}

// ---------------------------------------------------------------------------
// WasmCharacterController
// ---------------------------------------------------------------------------

/// WASM wrapper for the kinematic capsule character controller.
///
/// Manages position, velocity, and provides a simplified `move_and_slide`
/// variant that accepts a JSON-encoded list of `WasmSweepHit`s (or empty).
///
/// The `[f64; 3]` / `Option<[f64; 3]>` and `WasmCharacterShape` /
/// `WasmCharacterConfig` fields are private; use the `#[wasm_bindgen]`
/// accessors to manipulate them from JavaScript.
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmCharacterController {
    /// World-space position — private; use the JS accessors.
    pub(crate) position: [f64; 3],
    /// Linear velocity — private; use the JS accessors.
    pub(crate) velocity: [f64; 3],
    /// `true` when standing on a walkable surface.
    pub is_grounded: bool,
    /// Ground surface normal — private; use the JS accessor.
    pub(crate) ground_normal: Option<[f64; 3]>,
    /// Capsule geometry — private; use `get_shape` / `set_shape` from JS.
    pub(crate) shape: WasmCharacterShape,
    /// Tuning knobs — private; use `get_config` / `set_config` from JS.
    pub(crate) config: WasmCharacterConfig,
}

impl WasmCharacterController {
    /// Create a controller at the given world-space position.
    pub fn new(position: [f64; 3], shape: WasmCharacterShape, config: WasmCharacterConfig) -> Self {
        WasmCharacterController {
            position,
            velocity: [0.0; 3],
            is_grounded: false,
            ground_normal: None,
            shape,
            config,
        }
    }

    /// Apply gravity acceleration downward along `-up_axis` for `dt` seconds.
    pub fn apply_gravity(&mut self, gravity_accel: f64, dt: f64) {
        let up = self.config.up_axis;
        let down = [-up[0], -up[1], -up[2]];
        self.velocity = add(self.velocity, scale(down, gravity_accel * dt));
    }

    /// Initiate a jump with the given vertical impulse speed.
    pub fn jump(&mut self, speed: f64) {
        let up = self.config.up_axis;
        // Zero vertical component, then add jump speed upward.
        let vert = dot(self.velocity, up);
        self.velocity = add(self.velocity, scale(up, speed - vert));
        self.is_grounded = false;
        self.ground_normal = None;
    }

    /// Simplified move-and-slide that slides along the first hit surface provided.
    ///
    /// `desired` is the desired displacement `[x, y, z]`.
    /// `hits` is an optional list of sweep results from the caller's collision system.
    /// Returns the computed character move result.
    pub fn move_and_slide_simple(
        &mut self,
        desired: [f64; 3],
        hits: &[WasmSweepHit],
    ) -> WasmCharacterMove {
        let up = self.config.up_axis;
        let cos_threshold = self.config.max_slope_deg.to_radians().cos();

        let mut translation = desired;
        let mut is_grounded = false;
        let mut ground_normal: Option<[f64; 3]> = None;

        for hit in hits.iter().take(self.config.max_iterations) {
            let normal = hit.normal;
            let n_dot_up = dot(normal, up);
            if n_dot_up > cos_threshold {
                // Ground
                is_grounded = true;
                ground_normal = Some(normal);
            }
            // Project translation onto the hit plane
            translation = project_onto_plane(translation, normal);
        }

        // Update controller state
        self.is_grounded = is_grounded;
        self.ground_normal = ground_normal;
        self.position = add(self.position, translation);

        // Update velocity: remove component opposing the hit normals.
        for hit in hits.iter().take(self.config.max_iterations) {
            let n = hit.normal;
            let v_dot_n = dot(self.velocity, n);
            if v_dot_n < 0.0 {
                self.velocity = add(self.velocity, scale(n, -v_dot_n));
            }
        }

        WasmCharacterMove {
            translation,
            is_grounded,
            ground_normal,
        }
    }

    /// Serialise the controller state as a JSON string.
    pub fn to_json(&self) -> String {
        serde_json::to_string(self).unwrap_or_default()
    }
}

impl Default for WasmCharacterController {
    fn default() -> Self {
        WasmCharacterController::new(
            [0.0; 3],
            WasmCharacterShape::default(),
            WasmCharacterConfig::default(),
        )
    }
}

// ---------------------------------------------------------------------------
// wasm-bindgen JavaScript API
// ---------------------------------------------------------------------------

#[wasm_bindgen]
impl WasmCharacterController {
    /// Create a controller with default shape and config (JS constructor).
    #[wasm_bindgen(constructor)]
    pub fn new_js(x: f64, y: f64, z: f64) -> WasmCharacterController {
        WasmCharacterController::new(
            [x, y, z],
            WasmCharacterShape::default(),
            WasmCharacterConfig::default(),
        )
    }

    /// Create a controller with explicit shape and config.
    #[wasm_bindgen(js_name = "with_shape_and_config")]
    pub fn with_shape_and_config_js(
        x: f64,
        y: f64,
        z: f64,
        shape: &WasmCharacterShape,
        config: &WasmCharacterConfig,
    ) -> WasmCharacterController {
        WasmCharacterController::new([x, y, z], shape.clone(), config.clone())
    }

    /// World-space position as `[x, y, z]` (`Float64Array`).
    #[wasm_bindgen(js_name = "get_position")]
    pub fn get_position_js(&self) -> Vec<f64> {
        self.position.to_vec()
    }

    /// Teleport the controller to `(x, y, z)`.
    #[wasm_bindgen(js_name = "set_position")]
    pub fn set_position_js(&mut self, x: f64, y: f64, z: f64) {
        self.position = [x, y, z];
    }

    /// Linear velocity as `[vx, vy, vz]` (`Float64Array`).
    #[wasm_bindgen(js_name = "get_velocity")]
    pub fn get_velocity_js(&self) -> Vec<f64> {
        self.velocity.to_vec()
    }

    /// Set the linear velocity.
    #[wasm_bindgen(js_name = "set_velocity")]
    pub fn set_velocity_js(&mut self, vx: f64, vy: f64, vz: f64) {
        self.velocity = [vx, vy, vz];
    }

    /// `true` when currently on a walkable surface.
    #[wasm_bindgen(js_name = "is_grounded")]
    pub fn is_grounded_js(&self) -> bool {
        self.is_grounded
    }

    /// Ground normal as `[nx, ny, nz]`, or empty when not grounded.
    #[wasm_bindgen(js_name = "get_ground_normal")]
    pub fn get_ground_normal_js(&self) -> Vec<f64> {
        opt_vec3_to_js(self.ground_normal)
    }

    /// Snapshot of the capsule shape.
    #[wasm_bindgen(js_name = "get_shape")]
    pub fn get_shape_js(&self) -> WasmCharacterShape {
        self.shape.clone()
    }

    /// Replace the capsule shape.
    #[wasm_bindgen(js_name = "set_shape")]
    pub fn set_shape_js(&mut self, shape: &WasmCharacterShape) {
        self.shape = shape.clone();
    }

    /// Snapshot of the controller config.
    #[wasm_bindgen(js_name = "get_config")]
    pub fn get_config_js(&self) -> WasmCharacterConfig {
        self.config.clone()
    }

    /// Replace the controller config.
    #[wasm_bindgen(js_name = "set_config")]
    pub fn set_config_js(&mut self, config: &WasmCharacterConfig) {
        self.config = config.clone();
    }

    /// Apply gravity acceleration along `-up_axis` for `dt` seconds.
    #[wasm_bindgen(js_name = "apply_gravity")]
    pub fn apply_gravity_js(&mut self, gravity_accel: f64, dt: f64) {
        self.apply_gravity(gravity_accel, dt);
    }

    /// Initiate a jump with the given vertical impulse speed.
    #[wasm_bindgen(js_name = "jump")]
    pub fn jump_js(&mut self, speed: f64) {
        self.jump(speed);
    }

    /// Move-and-slide with no collision hits.
    ///
    /// Equivalent to calling `move_and_slide_simple` with an empty hit list.
    #[wasm_bindgen(js_name = "move_free")]
    pub fn move_free_js(&mut self, dx: f64, dy: f64, dz: f64) -> WasmCharacterMove {
        self.move_and_slide_simple([dx, dy, dz], &[])
    }

    /// Move-and-slide using a flat array of hit data.
    ///
    /// `hits_flat` is laid out as `[toi, nx, ny, nz, toi, nx, ny, nz, ...]`,
    /// i.e. four values per hit.  Pass an empty array when there are no hits.
    ///
    /// # Errors
    ///
    /// Returns a `JsValue` error when `hits_flat.len()` is not a multiple of 4.
    #[wasm_bindgen(js_name = "move_and_slide")]
    pub fn move_and_slide_js(
        &mut self,
        dx: f64,
        dy: f64,
        dz: f64,
        hits_flat: &[f64],
    ) -> std::result::Result<WasmCharacterMove, JsValue> {
        if !hits_flat.len().is_multiple_of(4) {
            return Err(JsValue::from_str(&format!(
                "hits_flat length {} is not a multiple of 4",
                hits_flat.len()
            )));
        }
        let hits: Vec<WasmSweepHit> = hits_flat
            .chunks_exact(4)
            .map(|c| WasmSweepHit {
                toi: c[0],
                normal: [c[1], c[2], c[3]],
            })
            .collect();
        Ok(self.move_and_slide_simple([dx, dy, dz], &hits))
    }

    /// Move-and-slide using a flat list of normals (each hit's `toi` is set
    /// to 0, since the simplified controller only consults `normal`).
    ///
    /// `normals_flat` is `[nx0, ny0, nz0, nx1, ny1, nz1, ...]` — i.e. three
    /// values per hit.
    ///
    /// # Errors
    ///
    /// Returns a `JsValue` error when `normals_flat.len()` is not a multiple of 3.
    #[wasm_bindgen(js_name = "move_with_normals")]
    pub fn move_with_normals_js(
        &mut self,
        dx: f64,
        dy: f64,
        dz: f64,
        normals_flat: &[f64],
    ) -> std::result::Result<WasmCharacterMove, JsValue> {
        let normals = unflatten_vec3s(normals_flat).map_err(err_to_jsvalue)?;
        let hits: Vec<WasmSweepHit> = normals
            .into_iter()
            .map(|normal| WasmSweepHit { toi: 0.0, normal })
            .collect();
        Ok(self.move_and_slide_simple([dx, dy, dz], &hits))
    }

    /// Serialise the controller as a `JsValue` (plain JS object).
    ///
    /// # Errors
    ///
    /// Returns a `JsValue` error string if serialisation fails.
    #[wasm_bindgen(js_name = "to_js_value")]
    pub fn to_js_value_js(&self) -> std::result::Result<JsValue, JsValue> {
        to_js_value(self)
    }

    /// Serialise the controller as a JSON string.
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
    fn test_character_bridge_instantiation() {
        let ctrl = WasmCharacterController::default();
        let json = ctrl.to_json();
        assert!(!json.is_empty(), "to_json should return non-empty string");
    }

    #[test]
    fn test_character_bridge_no_hit() {
        let mut ctrl = WasmCharacterController::default();
        let result = ctrl.move_and_slide_simple([3.0, 0.0, 0.0], &[]);
        assert!((result.translation[0] - 3.0).abs() < 1e-10);
        assert!(!result.is_grounded);
    }

    #[test]
    fn test_character_bridge_gravity() {
        let mut ctrl = WasmCharacterController::default();
        ctrl.apply_gravity(9.81, 1.0);
        assert!((ctrl.velocity[1] - (-9.81)).abs() < 1e-9);
    }

    #[test]
    fn test_character_bridge_jump() {
        let mut ctrl = WasmCharacterController {
            is_grounded: true,
            velocity: [2.0, 0.0, 0.0],
            ..Default::default()
        };
        ctrl.jump(5.0);
        assert!((ctrl.velocity[1] - 5.0).abs() < 1e-10);
        assert!((ctrl.velocity[0] - 2.0).abs() < 1e-10);
        assert!(!ctrl.is_grounded);
    }

    #[test]
    fn test_character_bridge_serde_roundtrip() {
        let original = WasmCharacterController {
            position: [1.0, 2.5, -3.0],
            velocity: [0.5, 0.0, 1.2],
            is_grounded: true,
            ground_normal: Some([0.0, 1.0, 0.0]),
            shape: WasmCharacterShape::new(0.4, 0.9),
            config: WasmCharacterConfig::default(),
        };
        let json = original.to_json();
        let restored: WasmCharacterController =
            serde_json::from_str(&json).expect("deserialize should succeed");
        assert!((restored.position[0] - original.position[0]).abs() < 1e-12);
        assert_eq!(restored.is_grounded, original.is_grounded);
        assert!(restored.ground_normal.is_some());
    }
}
