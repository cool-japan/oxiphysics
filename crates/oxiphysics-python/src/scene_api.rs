// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Python bindings for the oxiphysics scene module.
//!
//! Exposes declarative scene description with JSON round-trip to Python.

#![allow(missing_docs)]
#![allow(dead_code)]

use oxiphysics::scene::{SceneBuilder, SceneDescription};
use pyo3::prelude::*;

// ─────────────────────────────────────────────────────────────────────────────
// PySceneDescription
// ─────────────────────────────────────────────────────────────────────────────

/// A complete scene description (bodies + constraints).
#[pyclass(name = "SceneDescription")]
pub struct PySceneDescription {
    inner: SceneDescription,
}

#[pymethods]
impl PySceneDescription {
    /// Load a scene from a JSON string.
    #[staticmethod]
    pub fn from_json(json: &str) -> PyResult<Self> {
        SceneDescription::from_json(json)
            .map(|inner| Self { inner })
            .map_err(pyo3::exceptions::PyValueError::new_err)
    }

    /// Export this scene as a JSON string.
    pub fn to_json(&self) -> PyResult<String> {
        self.inner
            .to_json()
            .map_err(pyo3::exceptions::PyValueError::new_err)
    }

    /// Number of bodies in the scene.
    pub fn body_count(&self) -> usize {
        self.inner.body_count()
    }

    /// Number of static bodies.
    pub fn static_body_count(&self) -> usize {
        self.inner.static_body_count()
    }

    /// Number of dynamic bodies.
    pub fn dynamic_body_count(&self) -> usize {
        self.inner.dynamic_body_count()
    }

    /// Number of constraints.
    pub fn constraint_count(&self) -> usize {
        self.inner.constraint_count()
    }

    /// A human-readable summary of the scene.
    pub fn summary(&self) -> String {
        self.inner.summary()
    }

    /// All body IDs in the scene.
    pub fn body_ids(&self) -> Vec<String> {
        self.inner
            .body_ids()
            .into_iter()
            .map(|s| s.to_owned())
            .collect()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PySceneBuilder
// ─────────────────────────────────────────────────────────────────────────────

/// Fluent builder for constructing a SceneDescription.
#[pyclass(name = "SceneBuilder")]
pub struct PySceneBuilder {
    inner: Option<SceneBuilder>,
}

#[pymethods]
impl PySceneBuilder {
    /// Create a new SceneBuilder with the given scene name.
    #[new]
    pub fn new(name: &str) -> Self {
        Self {
            inner: Some(SceneBuilder::new(name)),
        }
    }

    /// Set the gravity vector.
    pub fn with_gravity(&mut self, gx: f64, gy: f64, gz: f64) {
        if let Some(b) = self.inner.take() {
            self.inner = Some(b.with_gravity([gx, gy, gz]));
        }
    }

    /// Set the simulation time step.
    pub fn with_dt(&mut self, dt: f64) {
        if let Some(b) = self.inner.take() {
            self.inner = Some(b.with_dt(dt));
        }
    }

    /// Add a sphere body (static or dynamic).
    pub fn add_sphere(&mut self, id: &str, x: f64, y: f64, z: f64, radius: f64, is_static: bool) {
        if let Some(b) = self.inner.take() {
            self.inner = Some(b.add_sphere(id, [x, y, z], radius, is_static));
        }
    }

    /// Add a ground plane (static, facing +Y).
    pub fn add_ground_plane(&mut self, id: &str) {
        if let Some(b) = self.inner.take() {
            self.inner = Some(b.add_ground_plane(id));
        }
    }

    /// Build and return the completed scene description.
    pub fn build(&mut self) -> PyResult<PySceneDescription> {
        self.inner
            .take()
            .map(|b| PySceneDescription { inner: b.build() })
            .ok_or_else(|| {
                pyo3::exceptions::PyRuntimeError::new_err("SceneBuilder already consumed")
            })
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Module registration
// ─────────────────────────────────────────────────────────────────────────────

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PySceneDescription>()?;
    m.add_class::<PySceneBuilder>()?;
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_scene_builder_instantiation() {
        let b = SceneBuilder::new("test_scene").with_gravity([0.0, -9.81, 0.0]);
        let scene = b.build();
        assert_eq!(scene.body_count(), 0);
    }

    #[test]
    fn test_scene_description_json_roundtrip() {
        let builder = SceneBuilder::new("roundtrip").with_gravity([0.0, -9.81, 0.0]);
        let scene = builder.build();
        let json = scene.to_json().expect("to_json failed");
        let wrapper = PySceneDescription::from_json(&json).expect("from_json failed");
        assert_eq!(wrapper.body_count(), 0);
    }

    #[test]
    fn test_scene_with_sphere() {
        let b = SceneBuilder::new("spheres").add_sphere("ball", [0.0, 5.0, 0.0], 0.5, false);
        let scene = b.build();
        assert_eq!(scene.body_count(), 1);
        assert_eq!(scene.dynamic_body_count(), 1);
    }
}
