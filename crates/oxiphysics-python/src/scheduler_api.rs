// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Python bindings for the oxiphysics scheduler module.
//!
//! Exposes priority-based physics step budget allocation per island to Python.

#![allow(missing_docs)]
#![allow(dead_code)]

use oxiphysics::scheduler::{Priority, Scheduler};
use pyo3::prelude::*;

// ─────────────────────────────────────────────────────────────────────────────
// PyScheduler
// ─────────────────────────────────────────────────────────────────────────────

/// Priority-based scheduler for physics simulation islands.
#[pyclass(name = "Scheduler")]
pub struct PyScheduler {
    inner: Scheduler,
}

#[pymethods]
impl PyScheduler {
    /// Create a new scheduler.
    ///
    /// `substeps` — default substeps per high-priority island.
    /// `max_substeps` — maximum substeps allowed for any island.
    #[new]
    pub fn new(substeps: u32, max_substeps: u32) -> Self {
        Self {
            inner: Scheduler::new(substeps, max_substeps),
        }
    }

    /// Create a new simulation island.
    ///
    /// `priority` should be one of: "Critical", "High", "Normal", "Low", "Sleeping".
    /// `body_count` — number of bodies in the island.
    /// Returns the island ID.
    pub fn add_island(&mut self, priority: &str, body_count: usize) -> PyResult<u32> {
        let p = parse_priority(priority)?;
        Ok(self.inner.add_island(p, body_count))
    }

    /// Remove an island.
    pub fn remove_island(&mut self, id: u32) {
        self.inner.remove_island(id);
    }

    /// Set priority for an existing island.
    pub fn set_priority(&mut self, id: u32, priority: &str) -> PyResult<()> {
        let p = parse_priority(priority)?;
        self.inner.set_priority(id, p);
        Ok(())
    }

    /// Update the velocity metric for an island (used for auto-sleep decisions).
    pub fn update_velocity(&mut self, id: u32, velocity: f64) {
        self.inner.update_velocity(id, velocity);
    }

    /// Update the body count for an island.
    pub fn update_body_count(&mut self, id: u32, count: usize) {
        self.inner.update_body_count(id, count);
    }

    /// Number of registered islands.
    pub fn island_count(&self) -> usize {
        self.inner.island_count()
    }

    /// Number of active (non-sleeping) islands.
    pub fn active_count(&self) -> usize {
        self.inner.active_count()
    }

    /// Compute the frame schedule as JSON.
    pub fn schedule_json(&self, dt: f64) -> PyResult<String> {
        let schedule = self.inner.schedule(dt);
        serde_json::to_string(&schedule)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))
    }

    /// Auto-sleep low-velocity islands. Returns list of put-to-sleep island IDs.
    pub fn auto_sleep_step(&mut self) -> Vec<u32> {
        self.inner.auto_sleep_step()
    }
}

fn parse_priority(s: &str) -> PyResult<Priority> {
    match s {
        "Critical" => Ok(Priority::Critical),
        "High" => Ok(Priority::High),
        "Normal" => Ok(Priority::Normal),
        "Low" => Ok(Priority::Low),
        "Sleeping" => Ok(Priority::Sleeping),
        other => Err(pyo3::exceptions::PyValueError::new_err(format!(
            "Unknown priority '{}'. Use: Critical, High, Normal, Low, Sleeping",
            other
        ))),
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Module registration
// ─────────────────────────────────────────────────────────────────────────────

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyScheduler>()?;
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_scheduler_instantiation() {
        let mut sched = PyScheduler::new(4, 8);
        assert_eq!(sched.island_count(), 0);
        let id = sched.add_island("High", 10).expect("add_island failed");
        assert_eq!(sched.island_count(), 1);
        sched.remove_island(id);
        assert_eq!(sched.island_count(), 0);
    }

    #[test]
    fn test_scheduler_schedule() {
        let mut sched = PyScheduler::new(4, 8);
        sched.add_island("Normal", 5).expect("add_island failed");
        sched.add_island("Low", 2).expect("add_island failed");
        let json = sched.schedule_json(1.0 / 60.0).expect("schedule failed");
        assert!(!json.is_empty());
    }
}
