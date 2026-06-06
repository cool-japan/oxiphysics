// Copyright (C) 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! `#[wasm_bindgen]` wrapper for [`WasmCharacterController`] — exposes the
//! capsule character controller to JavaScript with a flat-array API surface.

use wasm_bindgen::prelude::*;

use crate::character_bridge::{
    WasmCharacterConfig, WasmCharacterController, WasmCharacterShape, WasmSweepHit,
};

/// JavaScript-accessible wrapper for the kinematic capsule character controller.
///
/// ## JavaScript example
///
/// ```js
/// // Controller at (0, 2, 0), radius 0.4 m, half-height 0.9 m
/// const ctrl = WasmCharacterControllerJs.new(0.0, 2.0, 0.0, 0.4, 0.9);
///
/// // Each frame:
/// ctrl.apply_gravity(9.81, dt);
/// const result = ctrl.move_and_slide(dx, 0, dz, []);  // no hits (simplified)
/// const pos = ctrl.get_position();  // Float64Array [x, y, z]
/// ```
#[wasm_bindgen]
pub struct WasmCharacterControllerJs {
    inner: WasmCharacterController,
}

#[wasm_bindgen]
impl WasmCharacterControllerJs {
    /// Create a character controller at `(x, y, z)`.
    ///
    /// `radius` — capsule radius in metres.
    /// `half_height` — half-height of the cylindrical segment in metres.
    pub fn new(x: f64, y: f64, z: f64, radius: f64, half_height: f64) -> Self {
        WasmCharacterControllerJs {
            inner: WasmCharacterController::new(
                [x, y, z],
                WasmCharacterShape::new(radius, half_height),
                WasmCharacterConfig::default(),
            ),
        }
    }

    /// Create a character controller at the origin with default capsule dimensions.
    pub fn new_default() -> Self {
        WasmCharacterControllerJs {
            inner: WasmCharacterController::default(),
        }
    }

    /// Apply gravity downward for `dt` seconds.
    pub fn apply_gravity(&mut self, gravity_accel: f64, dt: f64) {
        self.inner.apply_gravity(gravity_accel, dt);
    }

    /// Initiate a jump with the given upward speed (m/s).
    pub fn jump(&mut self, speed: f64) {
        self.inner.jump(speed);
    }

    /// Move and slide along the given desired displacement `(dx, dy, dz)`.
    ///
    /// `hits_flat` encodes sweep-hit results as `[toi0, nx0, ny0, nz0, toi1, ...]`
    /// (4 f64 per hit). Pass an empty array when no sweep results are available.
    ///
    /// Returns `[tx, ty, tz, is_grounded]` where `is_grounded` is 1.0 or 0.0.
    pub fn move_and_slide(&mut self, dx: f64, dy: f64, dz: f64, hits_flat: Vec<f64>) -> Vec<f64> {
        let hits: Vec<WasmSweepHit> = hits_flat
            .chunks_exact(4)
            .map(|c| WasmSweepHit {
                toi: c[0],
                normal: [c[1], c[2], c[3]],
            })
            .collect();
        let result = self.inner.move_and_slide_simple([dx, dy, dz], &hits);
        vec![
            result.translation[0],
            result.translation[1],
            result.translation[2],
            if result.is_grounded { 1.0 } else { 0.0 },
        ]
    }

    /// Current position as `[x, y, z]`.
    pub fn get_position(&self) -> Vec<f64> {
        self.inner.position.to_vec()
    }

    /// Current velocity as `[vx, vy, vz]`.
    pub fn get_velocity(&self) -> Vec<f64> {
        self.inner.velocity.to_vec()
    }

    /// Set position directly (teleport).
    pub fn set_position(&mut self, x: f64, y: f64, z: f64) {
        self.inner.position = [x, y, z];
    }

    /// Set velocity directly.
    pub fn set_velocity(&mut self, vx: f64, vy: f64, vz: f64) {
        self.inner.velocity = [vx, vy, vz];
    }

    /// Whether the controller is currently standing on walkable ground.
    pub fn is_grounded(&self) -> bool {
        self.inner.is_grounded
    }

    /// Serialise the controller state as a JSON string.
    pub fn to_json(&self) -> String {
        self.inner.to_json()
    }
}
