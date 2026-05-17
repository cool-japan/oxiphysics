// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Python bindings for the oxiphysics trigger module.
//!
//! Exposes trigger/sensor volumes (enter/exit/stay events) to Python.

#![allow(missing_docs)]
#![allow(dead_code)]

use oxiphysics::trigger::{BodyEntry, TriggerEvent, TriggerWorld};
use pyo3::prelude::*;

// ─────────────────────────────────────────────────────────────────────────────
// Helpers
// ─────────────────────────────────────────────────────────────────────────────

fn trigger_event_to_json(ev: &TriggerEvent) -> String {
    match ev {
        TriggerEvent::Enter {
            trigger_id,
            body_index,
            step,
        } => format!(
            "{{\"kind\":\"Enter\",\"trigger_id\":{},\"body_index\":{},\"step\":{}}}",
            trigger_id, body_index, step
        ),
        TriggerEvent::Exit {
            trigger_id,
            body_index,
            step,
        } => format!(
            "{{\"kind\":\"Exit\",\"trigger_id\":{},\"body_index\":{},\"step\":{}}}",
            trigger_id, body_index, step
        ),
        TriggerEvent::Stay {
            trigger_id,
            body_index,
            step,
        } => format!(
            "{{\"kind\":\"Stay\",\"trigger_id\":{},\"body_index\":{},\"step\":{}}}",
            trigger_id, body_index, step
        ),
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PyTriggerWorld
// ─────────────────────────────────────────────────────────────────────────────

/// A world of trigger/sensor volumes.
#[pyclass(name = "TriggerWorld")]
pub struct PyTriggerWorld {
    inner: TriggerWorld,
}

impl Default for PyTriggerWorld {
    fn default() -> Self {
        Self {
            inner: TriggerWorld::new(),
        }
    }
}

#[pymethods]
impl PyTriggerWorld {
    /// Create a new trigger world.
    #[new]
    pub fn new() -> Self {
        Self::default()
    }

    /// Enable or disable Stay event emission (default: false).
    pub fn set_emit_stay(&mut self, emit: bool) {
        self.inner.set_emit_stay(emit);
    }

    /// Add a sphere trigger volume. Returns the trigger ID.
    pub fn add_sphere(&mut self, center: [f64; 3], radius: f64, tags: Vec<String>) -> u32 {
        self.inner.add_sphere(center, radius, tags)
    }

    /// Add an AABB trigger volume. Returns the trigger ID.
    pub fn add_aabb(&mut self, min: [f64; 3], max: [f64; 3], tags: Vec<String>) -> u32 {
        self.inner.add_aabb(min, max, tags)
    }

    /// Remove a trigger volume by ID. Returns true if it existed.
    pub fn remove_volume(&mut self, id: u32) -> bool {
        self.inner.remove_volume(id)
    }

    /// Enable a trigger volume.
    pub fn enable_volume(&mut self, id: u32) {
        self.inner.enable_volume(id);
    }

    /// Disable a trigger volume (it will not generate events).
    pub fn disable_volume(&mut self, id: u32) {
        self.inner.disable_volume(id);
    }

    /// Number of registered trigger volumes.
    pub fn volume_count(&self) -> usize {
        self.inner.volume_count()
    }

    /// Total number of active body occupancies across all triggers.
    pub fn active_occupancies(&self) -> usize {
        self.inner.active_occupancies()
    }

    /// Get all bodies currently inside a trigger volume.
    pub fn bodies_in(&self, trigger_id: u32) -> Vec<usize> {
        self.inner.bodies_in(trigger_id)
    }

    /// Get all trigger IDs that contain a given body.
    pub fn triggers_containing(&self, body_index: usize) -> Vec<u32> {
        self.inner.triggers_containing(body_index)
    }

    /// Update trigger state for a new tick and return events as a JSON array.
    ///
    /// `bodies` is a list of `(body_index, [x, y, z], radius, is_sleeping)` tuples.
    pub fn update_json(&mut self, tick: u64, bodies: Vec<(usize, [f64; 3], f64, bool)>) -> String {
        let entries: Vec<BodyEntry> = bodies
            .into_iter()
            .map(|(index, position, radius, is_sleeping)| BodyEntry {
                index,
                position,
                radius,
                is_sleeping,
            })
            .collect();
        let events = self.inner.update(tick, &entries);
        let parts: Vec<String> = events.iter().map(trigger_event_to_json).collect();
        format!("[{}]", parts.join(","))
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Module registration
// ─────────────────────────────────────────────────────────────────────────────

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyTriggerWorld>()?;
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_trigger_world_instantiation() {
        let mut tw = PyTriggerWorld::new();
        assert_eq!(tw.volume_count(), 0);
        let id = tw.add_sphere([0.0, 0.0, 0.0], 2.0, vec!["zone".to_owned()]);
        assert_eq!(tw.volume_count(), 1);
        let removed = tw.remove_volume(id);
        assert!(removed);
        assert_eq!(tw.volume_count(), 0);
    }

    #[test]
    fn test_trigger_world_enter_event() {
        let mut tw = PyTriggerWorld::new();
        tw.add_sphere([0.0, 0.0, 0.0], 2.0, vec![]);
        // Body starts outside
        let _before = tw.update_json(0, vec![(0, [10.0, 0.0, 0.0], 0.5, false)]);
        // Body enters
        let events_enter = tw.update_json(1, vec![(0, [0.5, 0.0, 0.0], 0.5, false)]);
        assert!(
            events_enter.contains("Enter"),
            "expected Enter event, got: {}",
            events_enter
        );
    }
}
