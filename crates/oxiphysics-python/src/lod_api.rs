// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Python bindings for the oxiphysics lod module.
//!
//! Exposes Level-of-Detail simulation tier management to Python.

use oxiphysics::lod::{LodConfig, LodSystem, LodTier};
use pyo3::prelude::*;

// ─────────────────────────────────────────────────────────────────────────────
// PyLodSystem
// ─────────────────────────────────────────────────────────────────────────────

/// Level-of-Detail system for scaling simulation cost with camera distance.
#[pyclass(name = "LodSystem")]
pub struct PyLodSystem {
    inner: LodSystem,
}

#[pymethods]
impl PyLodSystem {
    /// Create a new LOD system with default config and focus point.
    ///
    /// Default config thresholds: Full<50m, Reduced<150m, Minimal<300m, Frozen>=300m.
    #[new]
    pub fn new(focus: [f64; 3]) -> Self {
        Self {
            inner: LodSystem::new(LodConfig::default(), focus),
        }
    }

    /// Add a body to the LOD system. Returns its ID.
    pub fn add_body(&mut self, position: [f64; 3], priority: f32) -> u32 {
        self.inner.add_body(position, priority)
    }

    /// Remove a body.
    pub fn remove_body(&mut self, id: u32) {
        self.inner.remove_body(id);
    }

    /// Update a body's world position.
    pub fn update_position(&mut self, id: u32, position: [f64; 3]) {
        self.inner.update_position(id, position);
    }

    /// Update the focus/camera position.
    pub fn set_focus(&mut self, focus: [f64; 3]) {
        self.inner.set_focus(focus);
    }

    /// Update LOD tiers and return a list of changed bodies as JSON.
    pub fn update_json(&mut self) -> PyResult<String> {
        let updates = self.inner.update();
        serde_json::to_string(&updates)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))
    }

    /// Get the current LOD tier for a body as a string, or None if not found.
    pub fn tier_for_body(&self, id: u32) -> Option<String> {
        self.inner.tier(id).map(|t| format!("{:?}", t))
    }

    /// Get the substep count for a body.
    pub fn substeps_for_body(&self, id: u32) -> u32 {
        self.inner.substeps_for(id)
    }

    /// All body IDs in a specific tier.
    ///
    /// `tier` should be one of: "Full", "Reduced", "Minimal", "Frozen".
    pub fn bodies_in_tier(&self, tier: &str) -> PyResult<Vec<u32>> {
        let t = parse_tier(tier)?;
        Ok(self.inner.bodies_in_tier(t))
    }

    /// Number of bodies.
    pub fn body_count(&self) -> usize {
        self.inner.body_count()
    }

    /// Tier distribution counts as `(full, reduced, minimal, frozen)`.
    pub fn tier_counts(&self) -> (usize, usize, usize, usize) {
        self.inner.tier_counts()
    }

    /// Get distance to focus for a body, or None.
    pub fn distance_to_focus(&self, id: u32) -> Option<f64> {
        self.inner.distance_to_focus(id)
    }
}

fn parse_tier(s: &str) -> PyResult<LodTier> {
    match s {
        "Full" => Ok(LodTier::Full),
        "Reduced" => Ok(LodTier::Reduced),
        "Minimal" => Ok(LodTier::Minimal),
        "Frozen" => Ok(LodTier::Frozen),
        other => Err(pyo3::exceptions::PyValueError::new_err(format!(
            "Unknown LOD tier '{}'. Use: Full, Reduced, Minimal, Frozen",
            other
        ))),
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Module registration
// ─────────────────────────────────────────────────────────────────────────────

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyLodSystem>()?;
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lod_system_instantiation() {
        let mut lod = PyLodSystem::new([0.0, 0.0, 0.0]);
        assert_eq!(lod.body_count(), 0);
        let id = lod.add_body([0.0, 0.0, 0.0], 1.0);
        assert_eq!(lod.body_count(), 1);
        lod.remove_body(id);
        assert_eq!(lod.body_count(), 0);
    }

    #[test]
    fn test_lod_system_tiers() {
        let mut lod = PyLodSystem::new([0.0, 0.0, 0.0]);
        // Close body should be Full tier
        let id_close = lod.add_body([0.0, 0.0, 5.0], 1.0);
        // Far body should be Low tier
        let id_far = lod.add_body([0.0, 0.0, 200.0], 1.0);
        let _updates = lod.update_json().expect("update failed");
        let tier_close = lod.tier_for_body(id_close);
        assert!(tier_close.is_some());
        let tier_far = lod.tier_for_body(id_far);
        assert!(tier_far.is_some());
        let (full, _, _, _) = lod.tier_counts();
        assert!(full >= 1);
    }
}
