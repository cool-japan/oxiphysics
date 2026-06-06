// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Python bindings for the oxiphysics material_table module.
//!
//! Exposes runtime material interaction table with combine rules to Python.

use oxiphysics::material_table::{MaterialId, MaterialTable};
use pyo3::prelude::*;

// ─────────────────────────────────────────────────────────────────────────────
// PyMaterialTable
// ─────────────────────────────────────────────────────────────────────────────

/// Runtime material interaction table with built-in presets.
///
/// Pre-populated with: `concrete`, `rubber`, `metal`, `ice`, `wood`, `glass`.
#[pyclass(name = "MaterialTable")]
pub struct PyMaterialTable {
    inner: MaterialTable,
}

impl Default for PyMaterialTable {
    fn default() -> Self {
        Self {
            inner: MaterialTable::new(),
        }
    }
}

#[pymethods]
impl PyMaterialTable {
    /// Create a new material table (includes built-in presets).
    #[new]
    pub fn new() -> Self {
        Self::default()
    }

    /// Register a material with all properties. Returns its ID as u32.
    pub fn register_named(
        &mut self,
        name: &str,
        density: f64,
        restitution: f64,
        static_friction: f64,
        dynamic_friction: f64,
        linear_damping: f64,
        angular_damping: f64,
    ) -> u32 {
        self.inner
            .register_named(
                name,
                density,
                restitution,
                static_friction,
                dynamic_friction,
                linear_damping,
                angular_damping,
            )
            .0
    }

    /// Look up a material ID by name. Returns None if not found.
    pub fn id_for_name(&self, name: &str) -> Option<u32> {
        self.inner.id_for_name(name).map(|id| id.0)
    }

    /// Number of registered materials.
    pub fn material_count(&self) -> usize {
        self.inner.material_count()
    }

    /// Number of pair overrides.
    pub fn pair_count(&self) -> usize {
        self.inner.pair_count()
    }

    /// Set a per-pair static friction override.
    pub fn set_pair_static_friction(&mut self, a: u32, b: u32, value: f64) {
        self.inner
            .set_pair_static_friction(MaterialId(a), MaterialId(b), value);
    }

    /// Set a per-pair restitution override.
    pub fn set_pair_restitution(&mut self, a: u32, b: u32, value: f64) {
        self.inner
            .set_pair_restitution(MaterialId(a), MaterialId(b), value);
    }

    /// Get effective contact static friction between two materials.
    pub fn contact_static_friction(&self, a: u32, b: u32) -> f64 {
        self.inner
            .contact_static_friction(MaterialId(a), MaterialId(b))
    }

    /// Get effective contact dynamic friction between two materials.
    pub fn contact_dynamic_friction(&self, a: u32, b: u32) -> f64 {
        self.inner
            .contact_dynamic_friction(MaterialId(a), MaterialId(b))
    }

    /// Get effective contact restitution between two materials.
    pub fn contact_restitution(&self, a: u32, b: u32) -> f64 {
        self.inner.contact_restitution(MaterialId(a), MaterialId(b))
    }

    /// Summary string.
    pub fn summary(&self) -> String {
        self.inner.summary()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Module registration
// ─────────────────────────────────────────────────────────────────────────────

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyMaterialTable>()?;
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_material_table_instantiation() {
        let table = PyMaterialTable::new();
        // Should have the 6 built-in presets
        assert!(table.material_count() >= 6);
        let concrete_id = table.id_for_name("concrete");
        assert!(concrete_id.is_some());
    }

    #[test]
    fn test_material_table_register_and_lookup() {
        let mut table = PyMaterialTable::new();
        let id = table.register_named("foam", 50.0, 0.8, 0.4, 0.3, 0.0, 0.0);
        let found = table.id_for_name("foam");
        assert_eq!(found, Some(id));
    }

    #[test]
    fn test_material_friction() {
        let table = PyMaterialTable::new();
        let c = table.id_for_name("concrete").expect("concrete not found");
        let r = table.id_for_name("rubber").expect("rubber not found");
        let friction = table.contact_static_friction(c, r);
        assert!(friction > 0.0);
        let rest = table.contact_restitution(c, r);
        assert!(rest >= 0.0);
    }
}
