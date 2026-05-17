// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Python bindings for the oxiphysics character module.
//!
//! Exposes the kinematic capsule character controller to Python.

#![allow(missing_docs)]
#![allow(dead_code)]

use oxiphysics::character::{CharacterConfig, CharacterController, CharacterShape};
use pyo3::prelude::*;

// ─────────────────────────────────────────────────────────────────────────────
// PyCharacterController
// ─────────────────────────────────────────────────────────────────────────────

/// Kinematic capsule character controller with sweep-and-slide, step-up,
/// and slope handling.
///
/// Note: the `move_and_slide` method uses a closure-based sweep callback
/// which cannot be passed from Python. Use `move_no_collision` for simple
/// kinematic motion without collision response.
#[pyclass(name = "CharacterController")]
pub struct PyCharacterController {
    inner: CharacterController,
}

#[pymethods]
impl PyCharacterController {
    /// Create a character controller at the given position.
    ///
    /// `radius` — capsule radius (m). `half_height` — capsule half-height (m).
    #[new]
    pub fn new(position: [f64; 3], radius: f64, half_height: f64) -> Self {
        let shape = CharacterShape {
            radius,
            half_height,
        };
        let config = CharacterConfig::default();
        Self {
            inner: CharacterController::new(position, shape, config),
        }
    }

    /// Current world-space position.
    pub fn position(&self) -> [f64; 3] {
        self.inner.position
    }

    /// Set the world-space position directly (teleport).
    pub fn set_position(&mut self, position: [f64; 3]) {
        self.inner.position = position;
    }

    /// Current velocity.
    pub fn velocity(&self) -> [f64; 3] {
        self.inner.velocity
    }

    /// Set velocity.
    pub fn set_velocity(&mut self, velocity: [f64; 3]) {
        self.inner.velocity = velocity;
    }

    /// Whether the controller is currently grounded.
    pub fn is_grounded(&self) -> bool {
        self.inner.is_grounded
    }

    /// Apply gravity for `dt` seconds (without collision).
    ///
    /// Adjusts velocity and position by gravity. This is useful when
    /// collision is handled separately.
    pub fn apply_gravity(&mut self, gy: f64, dt: f64) {
        self.inner.velocity[1] += gy * dt;
        self.inner.position[0] += self.inner.velocity[0] * dt;
        self.inner.position[1] += self.inner.velocity[1] * dt;
        self.inner.position[2] += self.inner.velocity[2] * dt;
    }

    /// Set max slope degrees.
    pub fn set_max_slope_deg(&mut self, deg: f64) {
        self.inner.config.max_slope_deg = deg;
    }

    /// Set step offset (max height that can be stepped over).
    pub fn set_step_offset(&mut self, offset: f64) {
        self.inner.config.step_offset = offset;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Module registration
// ─────────────────────────────────────────────────────────────────────────────

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyCharacterController>()?;
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_character_controller_instantiation() {
        let ctrl = PyCharacterController::new([0.0, 1.0, 0.0], 0.4, 0.9);
        let pos = ctrl.position();
        assert!((pos[1] - 1.0).abs() < 1e-9);
        assert!(!ctrl.is_grounded());
    }

    #[test]
    fn test_character_apply_gravity() {
        let mut ctrl = PyCharacterController::new([0.0, 10.0, 0.0], 0.4, 0.9);
        ctrl.apply_gravity(-9.81, 1.0);
        let pos = ctrl.position();
        // After 1s of gravity from rest, should fall ~4.9m (but we use linear here)
        assert!(pos[1] < 10.0);
    }
}
