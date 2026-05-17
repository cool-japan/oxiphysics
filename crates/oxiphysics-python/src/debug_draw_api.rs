// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Python bindings for the oxiphysics debug_draw module.
//!
//! Exposes renderer-agnostic debug draw command buffer to Python.

#![allow(missing_docs)]
#![allow(dead_code)]

#[cfg(feature = "numpy-bridge")]
use oxiphysics::debug_draw::{DebugDrawSession, DrawCommand, DrawDuration, DrawList};
#[cfg(not(feature = "numpy-bridge"))]
use oxiphysics::debug_draw::{DebugDrawSession, DrawDuration, DrawList};
use pyo3::prelude::*;

#[cfg(feature = "numpy-bridge")]
use numpy::{IntoPyArray, PyArray1, PyArrayMethods};
#[cfg(feature = "numpy-bridge")]
use pyo3::Bound;

// ─────────────────────────────────────────────────────────────────────────────
// PyDrawList
// ─────────────────────────────────────────────────────────────────────────────

/// A buffer of debug draw commands (lines, spheres, arrows, text).
#[pyclass(name = "DrawList")]
pub struct PyDrawList {
    inner: DrawList,
}

impl Default for PyDrawList {
    fn default() -> Self {
        Self {
            inner: DrawList::new(),
        }
    }
}

#[pymethods]
impl PyDrawList {
    /// Create a new empty draw list.
    #[new]
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a single-frame line from `start` to `end` with the given RGBA color `[r,g,b,a]`.
    pub fn add_line(&mut self, start: [f64; 3], end: [f64; 3], color: [f32; 4]) {
        self.inner.add_line(start, end, color, DrawDuration::Single);
    }

    /// Add a single-frame sphere at `center` with `radius`.
    pub fn add_sphere(&mut self, center: [f64; 3], radius: f64, color: [f32; 4]) {
        self.inner
            .add_sphere(center, radius, color, DrawDuration::Single);
    }

    /// Add a single-frame AABB (axis-aligned box).
    pub fn add_aabb(&mut self, min: [f64; 3], max: [f64; 3], color: [f32; 4]) {
        self.inner.add_aabb(min, max, color, DrawDuration::Single);
    }

    /// Add a single-frame arrow.
    pub fn add_arrow(&mut self, from: [f64; 3], to: [f64; 3], head_size: f64, color: [f32; 4]) {
        self.inner
            .add_arrow(from, to, head_size, color, DrawDuration::Single);
    }

    /// Add a single-frame text label.
    pub fn add_text(&mut self, position: [f64; 3], label: &str, color: [f32; 4]) {
        self.inner
            .add_text(position, label.to_owned(), color, DrawDuration::Single);
    }

    /// Number of commands currently buffered.
    pub fn len(&self) -> usize {
        self.inner.len()
    }

    /// `true` if no commands are buffered.
    pub fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }

    /// Clear all commands.
    pub fn clear(&mut self) {
        self.inner.clear();
    }

    /// Drain single-frame commands and return them as JSON.
    pub fn drain_single_json(&mut self) -> PyResult<String> {
        let cmds = self.inner.drain_single();
        serde_json::to_string(&cmds)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))
    }

    /// Return all commands (persistent + single) as JSON.
    pub fn to_json(&self) -> PyResult<String> {
        let cmds: Vec<_> = self.inner.iter().collect();
        serde_json::to_string(&cmds)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PyDebugDrawSession
// ─────────────────────────────────────────────────────────────────────────────

/// Manages per-step debug draw lists with persistent and single-frame commands.
#[pyclass(name = "DebugDrawSession")]
pub struct PyDebugDrawSession {
    inner: DebugDrawSession,
}

impl Default for PyDebugDrawSession {
    fn default() -> Self {
        Self {
            inner: DebugDrawSession::new(),
        }
    }
}

#[pymethods]
impl PyDebugDrawSession {
    /// Create a new debug draw session.
    #[new]
    pub fn new() -> Self {
        Self::default()
    }

    /// Begin a new step at the given step number.
    pub fn begin_step(&mut self, step: u64) {
        self.inner.begin_step(step);
    }

    /// Add a single-frame line.
    pub fn add_line(&mut self, start: [f64; 3], end: [f64; 3], color: [f32; 4]) {
        self.inner.line(start, end, color);
    }

    /// Add a single-frame sphere.
    pub fn add_sphere(&mut self, center: [f64; 3], radius: f64, color: [f32; 4]) {
        self.inner.sphere(center, radius, color);
    }

    /// Add a single-frame AABB.
    pub fn add_aabb(&mut self, min: [f64; 3], max: [f64; 3], color: [f32; 4]) {
        self.inner.aabb(min, max, color);
    }

    /// Add a single-frame arrow.
    pub fn add_arrow(&mut self, from: [f64; 3], to: [f64; 3], head_size: f64, color: [f32; 4]) {
        self.inner.arrow(from, to, head_size, color);
    }

    /// End the current step and return the drawn frame as JSON.
    pub fn end_step_json(&mut self) -> PyResult<String> {
        let frame = self.inner.end_step();
        serde_json::to_string(&frame)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))
    }

    /// Current step counter.
    pub fn step(&self) -> u64 {
        self.inner.step()
    }

    /// Collect the representative vertex position for every draw command
    /// currently buffered in this session and return them as a
    /// `numpy.ndarray` of shape `(n, 3)` and dtype `float64`.
    ///
    /// Vertex extraction semantics per command variant:
    /// - `Line`   → start vertex only (1 vertex)
    /// - `Arrow`  → `from` vertex only (1 vertex)
    /// - `Aabb`   → `min` vertex only (1 vertex)
    /// - `Sphere` → `center` (1 vertex)
    /// - `Cross`  → `center` (1 vertex)
    /// - `Text`   → `position` (1 vertex)
    ///
    /// If the session's draw list is empty, returns a `(0, 3)` array.
    #[cfg(feature = "numpy-bridge")]
    pub fn vertex_buffer_to_numpy<'py>(&self, py: Python<'py>) -> Bound<'py, numpy::PyArray2<f64>> {
        let mut flat: Vec<f64> = Vec::new();
        for cmd in self.inner.list.commands.iter() {
            let v = match cmd {
                DrawCommand::Line { start, .. } => *start,
                DrawCommand::Arrow { from, .. } => *from,
                DrawCommand::Aabb { min, .. } => *min,
                DrawCommand::Sphere { center, .. } => *center,
                DrawCommand::Cross { center, .. } => *center,
                DrawCommand::Text { position, .. } => *position,
            };
            flat.push(v[0]);
            flat.push(v[1]);
            flat.push(v[2]);
        }
        let n = flat.len() / 3;
        let arr1: Bound<'py, PyArray1<f64>> = flat.into_pyarray(py);
        arr1.reshape([n, 3])
            .expect("vertex flat buffer length is always divisible by 3 by construction")
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Module registration
// ─────────────────────────────────────────────────────────────────────────────

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyDrawList>()?;
    m.add_class::<PyDebugDrawSession>()?;
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_draw_list_instantiation() {
        let mut dl = PyDrawList::new();
        assert!(dl.is_empty());
        dl.add_line([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 0.0, 0.0, 1.0]);
        assert_eq!(dl.len(), 1);
        dl.clear();
        assert!(dl.is_empty());
    }

    #[test]
    fn test_debug_draw_session() {
        let mut session = PyDebugDrawSession::new();
        session.begin_step(0);
        session.add_sphere([0.0, 1.0, 0.0], 0.5, [0.0, 1.0, 0.0, 1.0]);
        let json = session.end_step_json().expect("end_step failed");
        assert!(!json.is_empty());
        assert_eq!(session.step(), 0);
    }
}
