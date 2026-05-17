// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Python bindings for the oxiphysics snapshot module.
//!
//! Exposes world-state snapshots with delta tracking to Python.

#![allow(missing_docs)]
#![allow(dead_code)]

use oxiphysics::snapshot::{BodySnapshot, SnapshotManager, WorldSnapshot};
use pyo3::prelude::*;

// ─────────────────────────────────────────────────────────────────────────────
// PyBodySnapshot
// ─────────────────────────────────────────────────────────────────────────────

/// State snapshot for a single body at one instant.
#[pyclass(name = "BodySnapshot")]
pub struct PyBodySnapshot {
    inner: BodySnapshot,
}

#[pymethods]
impl PyBodySnapshot {
    /// Create a body snapshot.
    #[new]
    pub fn new(
        index: usize,
        position: [f64; 3],
        rotation: [f64; 4],
        velocity: [f64; 3],
        angular_velocity: [f64; 3],
    ) -> Self {
        Self {
            inner: BodySnapshot::new(index, position, rotation, velocity, angular_velocity),
        }
    }

    /// Speed (magnitude of linear velocity).
    pub fn speed(&self) -> f64 {
        self.inner.speed()
    }

    /// Angular speed (magnitude of angular velocity).
    pub fn angular_speed(&self) -> f64 {
        self.inner.angular_speed()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PyWorldSnapshot
// ─────────────────────────────────────────────────────────────────────────────

/// Snapshot of the full world state at one simulation step.
#[pyclass(name = "WorldSnapshot")]
pub struct PyWorldSnapshot {
    inner: WorldSnapshot,
}

#[pymethods]
impl PyWorldSnapshot {
    /// Create a world snapshot for the given tick and sim time.
    #[new]
    pub fn new(tick: u64, sim_time: f64) -> Self {
        Self {
            inner: WorldSnapshot::new(tick, sim_time),
        }
    }

    /// Number of bodies recorded.
    pub fn body_count(&self) -> usize {
        self.inner.body_count()
    }

    /// Total kinetic energy (unit mass).
    pub fn total_kinetic_energy_unit_mass(&self) -> f64 {
        self.inner.total_kinetic_energy_unit_mass()
    }

    /// Export as JSON.
    pub fn to_json(&self) -> PyResult<String> {
        self.inner
            .to_json()
            .map_err(pyo3::exceptions::PyValueError::new_err)
    }

    /// Summary string.
    pub fn summary(&self) -> String {
        self.inner.summary()
    }

    /// Load from a JSON string.
    #[staticmethod]
    pub fn from_json(json: &str) -> PyResult<Self> {
        WorldSnapshot::from_json(json)
            .map(|inner| Self { inner })
            .map_err(pyo3::exceptions::PyValueError::new_err)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PySnapshotManager
// ─────────────────────────────────────────────────────────────────────────────

/// Ring-buffer history of world snapshots.
#[pyclass(name = "SnapshotManager")]
pub struct PySnapshotManager {
    inner: SnapshotManager,
}

#[pymethods]
impl PySnapshotManager {
    /// Create a manager that keeps at most `max_count` snapshots.
    #[new]
    pub fn new(max_count: usize) -> Self {
        Self {
            inner: SnapshotManager::new(max_count),
        }
    }

    /// Push a snapshot into the manager.
    pub fn push(&mut self, snapshot: &PyWorldSnapshot) {
        self.inner.push(snapshot.inner.clone());
    }

    /// Number of stored snapshots.
    pub fn len(&self) -> usize {
        self.inner.len()
    }

    /// `true` if no snapshots are stored.
    pub fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }

    /// Maximum capacity.
    pub fn max_count(&self) -> usize {
        self.inner.max_count()
    }

    /// Clear all snapshots.
    pub fn clear(&mut self) {
        self.inner.clear();
    }

    /// All tick numbers in the manager.
    pub fn steps(&self) -> Vec<u64> {
        self.inner.steps()
    }

    /// Get the latest snapshot as JSON, or None.
    pub fn latest_json(&self) -> PyResult<Option<String>> {
        match self.inner.latest() {
            None => Ok(None),
            Some(s) => serde_json::to_string(s)
                .map(Some)
                .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string())),
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Module registration
// ─────────────────────────────────────────────────────────────────────────────

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyBodySnapshot>()?;
    m.add_class::<PyWorldSnapshot>()?;
    m.add_class::<PySnapshotManager>()?;
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_world_snapshot_instantiation() {
        let snap = PyWorldSnapshot::new(1, 0.016);
        assert_eq!(snap.body_count(), 0);
        assert_eq!(snap.total_kinetic_energy_unit_mass(), 0.0);
    }

    #[test]
    fn test_snapshot_json_roundtrip() {
        let snap = PyWorldSnapshot::new(10, 0.16);
        let json = snap.to_json().expect("to_json failed");
        let snap2 = PyWorldSnapshot::from_json(&json).expect("from_json failed");
        assert_eq!(snap2.body_count(), 0);
    }

    #[test]
    fn test_snapshot_manager() {
        let mut mgr = PySnapshotManager::new(5);
        assert!(mgr.is_empty());
        let snap = PyWorldSnapshot::new(0, 0.0);
        mgr.push(&snap);
        assert_eq!(mgr.len(), 1);
        mgr.clear();
        assert!(mgr.is_empty());
    }
}
