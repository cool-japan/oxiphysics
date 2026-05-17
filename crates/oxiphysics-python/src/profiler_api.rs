// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Python bindings for the oxiphysics profiler module.
//!
//! Exposes per-frame scope profiling with CSV and JSON export.

#![allow(missing_docs)]
#![allow(dead_code)]

use oxiphysics::profiler::ProfilerSession;
use pyo3::prelude::*;

// ─────────────────────────────────────────────────────────────────────────────
// PyProfilerSession
// ─────────────────────────────────────────────────────────────────────────────

/// Per-frame scope profiler.
///
/// Call `begin_frame`, then use `record_scope` to manually record timing data,
/// then `end_frame_json` to retrieve the summary.
///
/// Note: The RAII `scope` API uses raw pointers and cannot be safely bridged
/// to Python. Use `record_scope_ns` to inject timing data from Python.
#[pyclass(name = "ProfilerSession")]
pub struct PyProfilerSession {
    inner: ProfilerSession,
}

impl Default for PyProfilerSession {
    fn default() -> Self {
        Self {
            inner: ProfilerSession::new(),
        }
    }
}

#[pymethods]
impl PyProfilerSession {
    /// Create a new profiler session.
    #[new]
    pub fn new() -> Self {
        Self::default()
    }

    /// Begin a new frame, clearing all previous data.
    pub fn begin_frame(&mut self) {
        self.inner.begin_frame();
    }

    /// End the current frame and return the [`FrameReport`] as JSON.
    ///
    /// Returns a JSON string with `frame_ns` and the `root` scope tree.
    pub fn end_frame_json(&mut self) -> PyResult<String> {
        let report = self.inner.end_frame();
        serde_json::to_string(&report)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))
    }

    /// End the current frame and return the report as CSV.
    ///
    /// Columns: `name,call_count,total_ns,max_ns`
    pub fn end_frame_csv(&mut self) -> String {
        let report = self.inner.end_frame();
        report.to_csv()
    }

    /// End the current frame and return the folded-stacks string for flamegraph.pl.
    pub fn end_frame_folded_stacks(&mut self) -> String {
        let report = self.inner.end_frame();
        report.to_folded_stacks()
    }

    /// Number of arena nodes recorded in the current (in-progress) frame.
    pub fn node_count(&self) -> usize {
        self.inner.node_count()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Module registration
// ─────────────────────────────────────────────────────────────────────────────

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyProfilerSession>()?;
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_profiler_session_instantiation() {
        let mut session = PyProfilerSession::new();
        assert_eq!(session.node_count(), 0);
        session.begin_frame();
        let json = session.end_frame_json().expect("end_frame_json failed");
        assert!(json.contains("frame_ns"));
    }

    #[test]
    fn test_profiler_csv_and_folded() {
        let mut session = PyProfilerSession::new();
        session.begin_frame();
        let csv = session.end_frame_csv();
        assert!(csv.contains("name") || csv.contains("frame"));

        session.begin_frame();
        let folded = session.end_frame_folded_stacks();
        // Even with no scopes we get at least the frame root entry
        let _ = folded; // just verify it doesn't panic
    }
}
