// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Python rigid body bindings.
//!
//! Provides Python-friendly wrappers for rigid body simulation:
//! bodies, colliders, joints, a physics world, and batch force APIs.
//! All types use plain arrays and primitive types (no nalgebra).

use pyo3::prelude::*;

pub mod bodies;
pub mod joints;
pub mod utils;

pub use bodies::{PyCollider, PyContact, PyRayResult, PyRigidBody, PyRigidBodySet, PyShapeType};
pub use joints::{PyJoint, PyJointSet, PyJointType, PyPhysicsPipeline, PyPipelineConfig, PyWorld};
pub use utils::{
    PyBatchForce, PyBodyFilter, PyIslandStats, PySleepManager, PyWakeEvent, WakeReason,
    compute_island_stats,
};

/// Register all `rigid` classes into a Python sub-module.
///
/// Called from the top-level `#[pymodule]` in `lib.rs`.
pub fn register_rigid_module(parent: &Bound<'_, PyModule>) -> PyResult<()> {
    use pyo3::types::PyModuleMethods;
    let child = PyModule::new(parent.py(), "rigid")?;
    child.add_class::<PyShapeType>()?;
    child.add_class::<PyRigidBody>()?;
    child.add_class::<PyRigidBodySet>()?;
    child.add_class::<PyCollider>()?;
    child.add_class::<PyContact>()?;
    child.add_class::<PyRayResult>()?;
    child.add_class::<PyJointType>()?;
    child.add_class::<PyJoint>()?;
    child.add_class::<PyJointSet>()?;
    child.add_class::<PyWorld>()?;
    child.add_class::<PyPipelineConfig>()?;
    child.add_class::<PyPhysicsPipeline>()?;
    child.add_class::<WakeReason>()?;
    child.add_class::<PyWakeEvent>()?;
    child.add_class::<PySleepManager>()?;
    child.add_class::<PyBodyFilter>()?;
    child.add_class::<PyIslandStats>()?;
    child.add_class::<PyBatchForce>()?;
    child.add_function(wrap_pyfunction!(compute_island_stats, &child)?)?;
    parent.add_submodule(&child)?;
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Shared Vec3 helpers (used by bodies sub-module)
// ─────────────────────────────────────────────────────────────────────────────

#[inline]
pub(crate) fn vec3_scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline]
pub(crate) fn vec3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
pub(crate) fn vec3_length(a: [f64; 3]) -> f64 {
    vec3_dot(a, a).sqrt()
}

#[inline]
pub(crate) fn quat_identity() -> [f64; 4] {
    [0.0, 0.0, 0.0, 1.0]
}
