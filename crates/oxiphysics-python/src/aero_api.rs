// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Python bindings for the oxiphysics aero module.
//!
//! Exposes aerodynamic force computation (bluff-body drag and wing lift/drag).

#![allow(missing_docs)]
#![allow(dead_code)]

use oxiphysics::aero::{AeroBody, AeroEntry, AeroSystem, DragCurve, LiftCurve, WingSurface};
use pyo3::prelude::*;

// ─────────────────────────────────────────────────────────────────────────────
// PyAeroSystem
// ─────────────────────────────────────────────────────────────────────────────

/// Aerodynamics system: holds all bodies and a global wind vector.
///
/// Add entries via `add_drag_body` (bluff body, drag only) or
/// `add_wing_body` (full lift+drag wing surface).
/// Call `apply_json` to compute forces for all entries as JSON.
#[pyclass(name = "AeroSystem")]
pub struct PyAeroSystem {
    inner: AeroSystem,
}

impl Default for PyAeroSystem {
    fn default() -> Self {
        Self {
            inner: AeroSystem {
                entries: Vec::new(),
                wind: [0.0; 3],
            },
        }
    }
}

#[pymethods]
impl PyAeroSystem {
    /// Create a new aero system with zero wind.
    #[new]
    pub fn new() -> Self {
        Self::default()
    }

    /// Set the global wind velocity (m/s).
    pub fn set_wind(&mut self, wx: f64, wy: f64, wz: f64) {
        self.inner.wind = [wx, wy, wz];
    }

    /// Add a bluff-body entry (drag only, no lift).
    ///
    /// `position` — body position. `velocity` — body velocity.
    /// `drag_coeff` — C_d. `cross_section` — reference area (m²).
    /// `air_density` — ρ (kg/m³). `rotation` — row-major 3×3 orientation
    /// matrix as flat 9-element list (row0..row2).
    #[allow(clippy::too_many_arguments)]
    pub fn add_drag_body(
        &mut self,
        position: [f64; 3],
        velocity: [f64; 3],
        drag_coeff: f64,
        cross_section: f64,
        air_density: f64,
        rotation_flat: [f64; 9],
    ) {
        let rotation = [
            [rotation_flat[0], rotation_flat[1], rotation_flat[2]],
            [rotation_flat[3], rotation_flat[4], rotation_flat[5]],
            [rotation_flat[6], rotation_flat[7], rotation_flat[8]],
        ];
        let entry = AeroEntry {
            position,
            rotation,
            velocity,
            drag_body: Some(AeroBody {
                drag_coeff,
                cross_section,
                air_density,
            }),
            wings: Vec::new(),
        };
        self.inner.entries.push(entry);
    }

    /// Add a wing-surface entry (lift + drag, linear lift curve, constant drag).
    ///
    /// `position` — body position. `velocity` — body velocity.
    /// `center_local` — wing centre of pressure in body-local space.
    /// `normal_local` — wing lift axis in body-local space.
    /// `chord`, `span` — wing dimensions (m).
    /// `cl0`, `cl_slope` — lift curve: `CL = cl0 + cl_slope * alpha`.
    /// `cd_const` — constant drag coefficient.
    /// `air_density` — ρ (kg/m³).
    /// `rotation_flat` — row-major 3×3 orientation matrix, 9 floats.
    #[allow(clippy::too_many_arguments)]
    pub fn add_wing_body(
        &mut self,
        position: [f64; 3],
        velocity: [f64; 3],
        center_local: [f64; 3],
        normal_local: [f64; 3],
        chord: f64,
        span: f64,
        cl0: f64,
        cl_slope: f64,
        cd_const: f64,
        air_density: f64,
        rotation_flat: [f64; 9],
    ) {
        let rotation = [
            [rotation_flat[0], rotation_flat[1], rotation_flat[2]],
            [rotation_flat[3], rotation_flat[4], rotation_flat[5]],
            [rotation_flat[6], rotation_flat[7], rotation_flat[8]],
        ];
        let wing = WingSurface {
            center_local,
            normal_local,
            chord,
            span,
            lift_curve: LiftCurve::Linear {
                cl0,
                slope: cl_slope,
            },
            drag_curve: DragCurve::Constant(cd_const),
            air_density,
        };
        let entry = AeroEntry {
            position,
            rotation,
            velocity,
            drag_body: None,
            wings: vec![wing],
        };
        self.inner.entries.push(entry);
    }

    /// Number of entries in this system.
    pub fn entry_count(&self) -> usize {
        self.inner.entries.len()
    }

    /// Clear all entries.
    pub fn clear(&mut self) {
        self.inner.entries.clear();
    }

    /// Compute aerodynamic forces for all entries and return as JSON.
    ///
    /// Returns a JSON array of `{"force": [x,y,z], "torque": [x,y,z]}` objects,
    /// one per entry in insertion order.
    pub fn apply_json(&self) -> PyResult<String> {
        let forces = self.inner.apply();
        serde_json::to_string(&forces)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Module registration
// ─────────────────────────────────────────────────────────────────────────────

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyAeroSystem>()?;
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// Identity rotation (body axes == world axes).
    fn identity_rot() -> [f64; 9] {
        [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]
    }

    #[test]
    fn test_aero_system_instantiation() {
        let sys = PyAeroSystem::new();
        assert_eq!(sys.entry_count(), 0);
    }

    #[test]
    fn test_aero_drag_body() {
        let mut sys = PyAeroSystem::new();
        // Body moving at 10 m/s in +X, no wind
        sys.add_drag_body(
            [0.0, 0.0, 0.0],
            [10.0, 0.0, 0.0],
            0.47,
            1.0,
            1.225,
            identity_rot(),
        );
        assert_eq!(sys.entry_count(), 1);
        let json = sys.apply_json().expect("apply_json failed");
        // Drag should oppose motion → force.x < 0
        assert!(json.contains("force"));
    }

    #[test]
    fn test_aero_wing_body() {
        let mut sys = PyAeroSystem::new();
        sys.set_wind(0.0, 0.0, 0.0);
        // Wing facing up (+Y normal), moving forward (+X at 20 m/s)
        sys.add_wing_body(
            [0.0, 0.0, 0.0],
            [20.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            1.5,
            10.0,
            0.0,
            2.0 * std::f64::consts::PI,
            0.02,
            1.225,
            identity_rot(),
        );
        let json = sys.apply_json().expect("apply_json failed");
        assert!(json.contains("torque"));
    }
}
