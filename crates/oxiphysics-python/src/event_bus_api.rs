// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Python bindings for the oxiphysics event_bus module.
//!
//! Exposes physics event publish/subscribe functionality to Python.
//!
//! Note: `EventBus` stores `Box<dyn Fn>` listeners which are `!Sync`.
//! The pyclass is therefore marked `unsendable` — it cannot be shared
//! across Python threads.

#![allow(missing_docs)]
#![allow(dead_code)]

use oxiphysics::event_bus::EventBus;
use pyo3::prelude::*;

// ─────────────────────────────────────────────────────────────────────────────
// PyEventBus
// ─────────────────────────────────────────────────────────────────────────────

/// Physics event bus for publish/subscribe patterns.
///
/// Events are buffered and dispatched in batch during `flush()`.
///
/// Note: this object is **not thread-safe** (unsendable) because it stores
/// listener closures that are `!Sync`.
#[pyclass(name = "EventBus", unsendable)]
pub struct PyEventBus {
    inner: EventBus,
}

impl Default for PyEventBus {
    fn default() -> Self {
        Self {
            inner: EventBus::new(),
        }
    }
}

#[pymethods]
impl PyEventBus {
    /// Create a new event bus.
    #[new]
    pub fn new() -> Self {
        Self::default()
    }

    /// Number of pending (unflushed) events.
    pub fn pending_count(&self) -> usize {
        self.inner.pending_count()
    }

    /// Number of registered listeners.
    pub fn listener_count(&self) -> usize {
        self.inner.listener_count()
    }

    /// Total number of events ever published.
    pub fn total_published(&self) -> u64 {
        self.inner.total_published()
    }

    /// Total number of events dispatched to listeners.
    pub fn total_dispatched(&self) -> u64 {
        self.inner.total_dispatched()
    }

    /// Drain and return all pending events as a list of description strings.
    ///
    /// Each string is the `Display` representation of the event.
    /// The internal queue is cleared.
    pub fn drain_as_strings(&mut self) -> Vec<String> {
        self.inner
            .drain()
            .into_iter()
            .map(|e| e.to_string())
            .collect()
    }

    /// Return the kind names of all pending events without consuming them.
    pub fn pending_kind_names(&self) -> Vec<String> {
        self.inner
            .pending()
            .iter()
            .map(|e| e.kind_name().to_owned())
            .collect()
    }

    /// Clear all pending events without dispatching them.
    pub fn clear(&mut self) {
        self.inner.clear();
    }

    /// Flush all pending events to registered listeners.
    pub fn flush(&mut self) {
        self.inner.flush();
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Module registration
// ─────────────────────────────────────────────────────────────────────────────

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyEventBus>()?;
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_event_bus_instantiation() {
        let mut bus = PyEventBus::new();
        assert_eq!(bus.pending_count(), 0);
        assert_eq!(bus.listener_count(), 0);
        assert_eq!(bus.total_published(), 0);
        bus.clear();
        assert_eq!(bus.pending_count(), 0);
    }

    #[test]
    fn test_event_bus_drain_empty() {
        let mut bus = PyEventBus::new();
        let strings = bus.drain_as_strings();
        assert!(strings.is_empty());
    }
}
