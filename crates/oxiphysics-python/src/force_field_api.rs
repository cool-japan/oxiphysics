// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Python bindings for the oxiphysics force_field module.
//!
//! Exposes spatial force fields (gravity wells, vortex, wind, explosion) to Python.

use oxiphysics::force_field::{AabbRegion, ForceFieldKind, ForceFieldSystem};
use pyo3::prelude::*;

// ─────────────────────────────────────────────────────────────────────────────
// PyForceFieldAabbRegion
// ─────────────────────────────────────────────────────────────────────────────

/// Axis-aligned bounding box region for force field spatial restriction.
#[pyclass(name = "ForceFieldAabbRegion")]
pub struct PyForceFieldAabbRegion {
    inner: AabbRegion,
}

#[pymethods]
impl PyForceFieldAabbRegion {
    /// Create an AABB region with min/max corners.
    #[new]
    pub fn new(min: [f64; 3], max: [f64; 3]) -> Self {
        Self {
            inner: AabbRegion::new(min, max),
        }
    }

    /// Check whether a point is inside this region.
    pub fn contains(&self, point: [f64; 3]) -> bool {
        self.inner.contains(point)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PyForceFieldSystem
// ─────────────────────────────────────────────────────────────────────────────

/// A system of spatial force fields applied to physics bodies.
#[pyclass(name = "ForceFieldSystem")]
pub struct PyForceFieldSystem {
    inner: ForceFieldSystem,
}

impl Default for PyForceFieldSystem {
    fn default() -> Self {
        Self {
            inner: ForceFieldSystem::new(),
        }
    }
}

#[pymethods]
impl PyForceFieldSystem {
    /// Create a new empty force field system.
    #[new]
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a uniform force field (constant vector). Returns the field ID.
    pub fn add_uniform(&mut self, force: [f64; 3]) -> usize {
        self.inner.add(ForceFieldKind::Uniform { force })
    }

    /// Add a radial attractor gravity well. Returns the field ID.
    pub fn add_radial_attract(
        &mut self,
        center: [f64; 3],
        strength: f64,
        falloff_exp: f64,
    ) -> usize {
        self.inner.add(ForceFieldKind::RadialAttract {
            center,
            strength,
            falloff_exp,
        })
    }

    /// Add a radial repulsor. Returns the field ID.
    pub fn add_radial_repel(&mut self, center: [f64; 3], strength: f64, falloff_exp: f64) -> usize {
        self.inner.add(ForceFieldKind::RadialRepel {
            center,
            strength,
            falloff_exp,
        })
    }

    /// Add a vortex field. Returns the field ID.
    pub fn add_vortex(
        &mut self,
        center: [f64; 3],
        axis: [f64; 3],
        tangential_strength: f64,
        axial_strength: f64,
    ) -> usize {
        self.inner.add(ForceFieldKind::Vortex {
            center,
            axis,
            tangential_strength,
            axial_strength,
        })
    }

    /// Add a wind field. Returns the field ID.
    pub fn add_wind(&mut self, direction: [f64; 3], base_speed: f64) -> usize {
        self.inner.add(ForceFieldKind::Wind {
            direction,
            base_speed,
        })
    }

    /// Add an explosion (expanding pressure wave). Returns the field ID.
    pub fn add_explosion(
        &mut self,
        center: [f64; 3],
        peak_force: f64,
        wave_speed: f64,
        thickness: f64,
        start_time: f64,
    ) -> usize {
        self.inner.add(ForceFieldKind::Explosion {
            center,
            peak_force,
            wave_speed,
            thickness,
            start_time,
        })
    }

    /// Add a turbulent (spatially-varying oscillating) force field. Returns the field ID.
    pub fn add_turbulent(&mut self, base_force: [f64; 3], amplitude: f64, frequency: f64) -> usize {
        self.inner.add(ForceFieldKind::Turbulent {
            base_force,
            amplitude,
            frequency,
        })
    }

    /// Enable a field by ID.
    pub fn enable(&mut self, id: usize) {
        self.inner.enable(id);
    }

    /// Disable a field by ID.
    pub fn disable(&mut self, id: usize) {
        self.inner.disable(id);
    }

    /// Remove a field by ID.
    pub fn remove(&mut self, id: usize) {
        self.inner.remove(id);
    }

    /// Compute the total force at position `pos` at simulation time `time`.
    pub fn force_at(&self, pos: [f64; 3], time: f64) -> [f64; 3] {
        self.inner.force_at(pos, time)
    }

    /// Number of registered fields (including disabled).
    pub fn len(&self) -> usize {
        self.inner.len()
    }

    /// `true` when no fields are registered.
    pub fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }

    /// Number of currently active (enabled) fields.
    pub fn active_count(&self) -> usize {
        self.inner.active_count()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Module registration
// ─────────────────────────────────────────────────────────────────────────────

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyForceFieldAabbRegion>()?;
    m.add_class::<PyForceFieldSystem>()?;
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_force_field_system_instantiation() {
        let mut sys = PyForceFieldSystem::new();
        assert!(sys.is_empty());
        let id = sys.add_wind([0.0, 0.0, 1.0], 10.0);
        assert_eq!(sys.len(), 1);
        assert_eq!(sys.active_count(), 1);
        sys.disable(id);
        assert_eq!(sys.active_count(), 0);
    }

    #[test]
    fn test_force_field_sample_radial() {
        let mut sys = PyForceFieldSystem::new();
        sys.add_radial_attract([0.0, 0.0, 0.0], 9.81, 2.0);
        let force = sys.force_at([1.0, 0.0, 0.0], 0.0);
        // Force should be non-zero and pointing toward the attractor
        let mag = (force[0] * force[0] + force[1] * force[1] + force[2] * force[2]).sqrt();
        assert!(mag > 0.0);
        assert!(force[0] < 0.0, "Force should point toward origin");
    }

    #[test]
    fn test_aabb_region_contains() {
        let r = PyForceFieldAabbRegion::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert!(r.contains([0.5, 0.5, 0.5]));
        assert!(!r.contains([2.0, 0.5, 0.5]));
    }
}
