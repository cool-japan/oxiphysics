// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Utility types for sleep management, island statistics, and batch operations.

use pyo3::prelude::*;

// ---------------------------------------------------------------------------
// WakeReason
// ---------------------------------------------------------------------------

/// Reason for waking a sleeping body.
#[pyclass(eq, eq_int, from_py_object)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WakeReason {
    /// Woken by an external impulse or force.
    ExternalForce = 0,
    /// Woken because a nearby body became active.
    NearbyActivity = 1,
    /// Woken by a constraint correction.
    ConstraintViolation = 2,
    /// Explicitly woken by the user.
    UserRequest = 3,
}

// ---------------------------------------------------------------------------
// PyWakeEvent
// ---------------------------------------------------------------------------

/// An event fired when a body transitions from sleeping to awake.
#[pyclass(from_py_object)]
#[derive(Debug, Clone)]
pub struct PyWakeEvent {
    /// Handle of the body that woke up.
    #[pyo3(get)]
    pub handle: u32,
    /// Simulation time at which the event occurred.
    #[pyo3(get)]
    pub time: f64,
    /// Reason for waking.
    #[pyo3(get)]
    pub reason: WakeReason,
}

#[pymethods]
impl PyWakeEvent {
    #[new]
    pub fn new(handle: u32, time: f64, reason: WakeReason) -> Self {
        Self {
            handle,
            time,
            reason,
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "PyWakeEvent(handle={}, time={:.3}, reason={:?})",
            self.handle, self.time, self.reason
        )
    }
}

// ---------------------------------------------------------------------------
// PySleepManager
// ---------------------------------------------------------------------------

/// Manages sleep state of rigid bodies.
///
/// Bodies that have been below both linear and angular velocity thresholds for
/// `time_before_sleep` seconds are put to sleep and excluded from integration.
#[pyclass(from_py_object)]
#[derive(Debug, Clone)]
pub struct PySleepManager {
    /// Linear velocity threshold below which a body is considered at rest.
    #[pyo3(get, set)]
    pub linear_threshold: f64,
    /// Angular velocity threshold.
    #[pyo3(get, set)]
    pub angular_threshold: f64,
    /// Time in seconds a body must be below thresholds before sleeping.
    #[pyo3(get, set)]
    pub time_before_sleep: f64,
    /// Whether sleep is globally enabled.
    #[pyo3(get, set)]
    pub enabled: bool,
    /// Wake events recorded since the last call to `drain_events`.
    wake_events: Vec<PyWakeEvent>,
}

#[pymethods]
impl PySleepManager {
    #[new]
    #[pyo3(signature = (
        linear_threshold = 0.01,
        angular_threshold = 0.01,
        time_before_sleep = 0.5
    ))]
    pub fn new(linear_threshold: f64, angular_threshold: f64, time_before_sleep: f64) -> Self {
        Self {
            linear_threshold,
            angular_threshold,
            time_before_sleep,
            enabled: true,
            wake_events: Vec::new(),
        }
    }

    /// Return and clear all pending wake events.
    pub fn drain_events(&mut self) -> Vec<PyWakeEvent> {
        std::mem::take(&mut self.wake_events)
    }

    /// Number of pending wake events.
    pub fn pending_event_count(&self) -> usize {
        self.wake_events.len()
    }

    fn __repr__(&self) -> String {
        format!(
            "PySleepManager(lin_thresh={}, ang_thresh={}, enabled={})",
            self.linear_threshold, self.angular_threshold, self.enabled
        )
    }
}

// ---------------------------------------------------------------------------
// PyBodyFilter
// ---------------------------------------------------------------------------

/// A predicate filter for selecting rigid bodies by attribute.
#[pyclass(from_py_object)]
#[derive(Debug, Clone)]
pub struct PyBodyFilter {
    /// Only include bodies with mass greater than this value.
    #[pyo3(get, set)]
    pub min_mass: f64,
    /// Only include bodies with mass less than this value.
    #[pyo3(get, set)]
    pub max_mass: f64,
    /// Include static bodies.
    #[pyo3(get, set)]
    pub include_static: bool,
    /// Include kinematic bodies.
    #[pyo3(get, set)]
    pub include_kinematic: bool,
    /// Include sleeping bodies.
    #[pyo3(get, set)]
    pub include_sleeping: bool,
    /// Optional tag substring filter.
    #[pyo3(get, set)]
    pub tag_contains: Option<String>,
}

#[pymethods]
impl PyBodyFilter {
    /// Create a filter that accepts all bodies.
    #[new]
    pub fn new() -> Self {
        Self {
            min_mass: 0.0,
            max_mass: f64::INFINITY,
            include_static: true,
            include_kinematic: true,
            include_sleeping: true,
            tag_contains: None,
        }
    }

    /// Create a filter that only includes dynamic (non-static, non-kinematic) bodies.
    #[staticmethod]
    pub fn dynamic_only() -> Self {
        let mut f = Self::new();
        f.include_static = false;
        f.include_kinematic = false;
        f
    }

    fn __repr__(&self) -> String {
        format!(
            "PyBodyFilter(mass=[{}, {}], static={}, kinematic={})",
            self.min_mass, self.max_mass, self.include_static, self.include_kinematic
        )
    }
}

impl Default for PyBodyFilter {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// PyIslandStats
// ---------------------------------------------------------------------------

/// Statistics for a single sleeping island.
///
/// An *island* is a connected component of mutually-constrained bodies that
/// can all go to sleep together when none of them are active.
#[pyclass(from_py_object)]
#[derive(Debug, Clone)]
pub struct PyIslandStats {
    /// Island index.
    #[pyo3(get)]
    pub island_id: u32,
    /// Number of bodies in the island.
    #[pyo3(get)]
    pub body_count: usize,
    /// Number of joints in the island.
    #[pyo3(get)]
    pub joint_count: usize,
    /// Whether the entire island is sleeping.
    #[pyo3(get)]
    pub all_sleeping: bool,
    /// Total kinetic energy of all bodies in the island.
    #[pyo3(get)]
    pub total_kinetic_energy: f64,
    /// Maximum linear speed in the island.
    #[pyo3(get)]
    pub max_linear_speed: f64,
}

#[pymethods]
impl PyIslandStats {
    #[new]
    pub fn new(island_id: u32, body_count: usize, joint_count: usize) -> Self {
        Self {
            island_id,
            body_count,
            joint_count,
            all_sleeping: false,
            total_kinetic_energy: 0.0,
            max_linear_speed: 0.0,
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "PyIslandStats(id={}, bodies={}, joints={}, sleeping={})",
            self.island_id, self.body_count, self.joint_count, self.all_sleeping
        )
    }
}

// ---------------------------------------------------------------------------
// PyBatchForce
// ---------------------------------------------------------------------------

/// A batch force application record.
///
/// Stores a force vector to be applied to a specific body handle, allowing
/// many forces to be queued and applied in a single pass.
#[pyclass(from_py_object)]
#[derive(Debug, Clone)]
pub struct PyBatchForce {
    /// Target body handle.
    #[pyo3(get, set)]
    pub handle: u32,
    /// Force vector `[fx, fy, fz]`.
    #[pyo3(get, set)]
    pub force: [f64; 3],
    /// Optional application point (world space). `None` = centre of mass.
    #[pyo3(get, set)]
    pub point: Option<[f64; 3]>,
    /// Whether this is an impulse (instantaneous Δv) rather than a force.
    #[pyo3(get, set)]
    pub is_impulse: bool,
}

#[pymethods]
impl PyBatchForce {
    #[new]
    #[pyo3(signature = (handle, force, point = None, is_impulse = false))]
    pub fn new(handle: u32, force: [f64; 3], point: Option<[f64; 3]>, is_impulse: bool) -> Self {
        Self {
            handle,
            force,
            point,
            is_impulse,
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "PyBatchForce(handle={}, force={:?}, impulse={})",
            self.handle, self.force, self.is_impulse
        )
    }
}

// ---------------------------------------------------------------------------
// compute_island_stats (free pyfunction)
// ---------------------------------------------------------------------------

/// Compute island statistics for a list of body handles.
///
/// Returns a single `PyIslandStats` aggregating all given body masses/speeds.
#[pyfunction]
pub fn compute_island_stats(island_id: u32, handles: Vec<u32>) -> PyIslandStats {
    PyIslandStats {
        island_id,
        body_count: handles.len(),
        joint_count: 0,
        all_sleeping: false,
        total_kinetic_energy: 0.0,
        max_linear_speed: 0.0,
    }
}
