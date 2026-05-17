// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Python bindings for the oxiphysics buoyancy module.
//!
//! Exposes Archimedes buoyancy and viscous drag for bodies in fluid volumes.

#![allow(missing_docs)]
#![allow(dead_code)]

use oxiphysics::buoyancy::{BuoyancyWorld, BuoyantBody, BuoyantShape};
use pyo3::prelude::*;

/// `(body_id, position, radius, mass_kg, velocity)` input tuple for `compute_forces_json`.
type BodyInput = (u32, [f64; 3], f64, f64, [f64; 3]);

// ─────────────────────────────────────────────────────────────────────────────
// PyBuoyancyWorld
// ─────────────────────────────────────────────────────────────────────────────

/// World managing buoyancy forces for multiple fluid volumes.
#[pyclass(name = "BuoyancyWorld")]
pub struct PyBuoyancyWorld {
    inner: BuoyancyWorld,
}

#[pymethods]
impl PyBuoyancyWorld {
    /// Create a new buoyancy world with the given gravity (m/s^2, typically 9.81).
    #[new]
    pub fn new(gravity: f64) -> Self {
        Self {
            inner: BuoyancyWorld::new(gravity),
        }
    }

    /// Add a fluid volume as an AABB and return its ID.
    pub fn add_fluid(&mut self, min: [f64; 3], max: [f64; 3], density: f64, surface_y: f64) -> u32 {
        self.inner.add_fluid(min, max, density, surface_y)
    }

    /// Add a fluid volume with explicit drag coefficients.
    pub fn add_fluid_with_drag(
        &mut self,
        min: [f64; 3],
        max: [f64; 3],
        density: f64,
        surface_y: f64,
        linear_drag: f64,
        angular_drag: f64,
    ) -> u32 {
        self.inner
            .add_fluid_with_drag(min, max, density, surface_y, linear_drag, angular_drag)
    }

    /// Remove a fluid volume by ID.
    pub fn remove_fluid(&mut self, id: u32) {
        self.inner.remove_fluid(id);
    }

    /// Number of fluid volumes.
    pub fn fluid_count(&self) -> usize {
        self.inner.fluid_count()
    }

    /// Compute buoyancy forces for a list of sphere bodies.
    ///
    /// `bodies` is a list of `(body_id, position [x,y,z], radius, mass_kg, velocity [x,y,z])`.
    /// Returns list of `(body_id, BuoyancyForce)` as JSON.
    pub fn compute_forces_json(&self, bodies: Vec<BodyInput>) -> PyResult<String> {
        let buoyant_bodies: Vec<BuoyantBody> = bodies
            .into_iter()
            .map(|(id, position, radius, mass, velocity)| BuoyantBody {
                id,
                position,
                velocity,
                angular_velocity: [0.0; 3],
                mass,
                shape: BuoyantShape::Sphere { radius },
            })
            .collect();
        let forces = self.inner.apply(&buoyant_bodies);
        serde_json::to_string(&forces)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Module registration
// ─────────────────────────────────────────────────────────────────────────────

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyBuoyancyWorld>()?;
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_buoyancy_world_instantiation() {
        let mut world = PyBuoyancyWorld::new(9.81);
        assert_eq!(world.fluid_count(), 0);
        let id = world.add_fluid([-5.0, -5.0, -5.0], [5.0, 0.0, 5.0], 1000.0, 0.0);
        assert_eq!(world.fluid_count(), 1);
        world.remove_fluid(id);
        assert_eq!(world.fluid_count(), 0);
    }

    #[test]
    fn test_buoyancy_force_computation() {
        let mut world = PyBuoyancyWorld::new(9.81);
        world.add_fluid([-10.0, -10.0, -10.0], [10.0, 1.0, 10.0], 1000.0, 0.0);
        // A small sphere inside the fluid (mass = 1 kg, buoyancy should produce upward force)
        let json = world
            .compute_forces_json(vec![(0, [0.0, -2.0, 0.0], 0.5, 1.0, [0.0, 0.0, 0.0])])
            .expect("compute_forces failed");
        assert!(!json.is_empty());
    }
}
