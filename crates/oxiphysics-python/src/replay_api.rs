// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Python bindings for the oxiphysics replay module.
//!
//! Exposes deterministic simulation record/replay functionality to Python.

use oxiphysics::replay::{ReplayRecord, SimRecorder, SimReplayer};
use pyo3::prelude::*;

// ─────────────────────────────────────────────────────────────────────────────
// PySimRecorder
// ─────────────────────────────────────────────────────────────────────────────

/// Records simulation commands for deterministic replay.
#[pyclass(name = "SimRecorder")]
pub struct PySimRecorder {
    inner: Option<SimRecorder>,
}

#[pymethods]
impl PySimRecorder {
    /// Create a new recorder starting at the given tick.
    #[new]
    pub fn new(start_tick: u64) -> Self {
        Self {
            inner: Some(SimRecorder::new(start_tick)),
        }
    }

    /// Begin recording a new time step with the given dt.
    pub fn begin_step(&mut self, dt: f64) {
        if let Some(ref mut rec) = self.inner {
            rec.begin_step(dt);
        }
    }

    /// Record a force applied to a body.
    pub fn record_force(&mut self, body: usize, force: [f64; 3]) {
        if let Some(ref mut rec) = self.inner {
            rec.record_force(body, force);
        }
    }

    /// Record an impulse applied to a body.
    pub fn record_impulse(&mut self, body: usize, impulse: [f64; 3]) {
        if let Some(ref mut rec) = self.inner {
            rec.record_impulse(body, impulse);
        }
    }

    /// Record a torque applied to a body.
    pub fn record_torque(&mut self, body: usize, torque: [f64; 3]) {
        if let Some(ref mut rec) = self.inner {
            rec.record_torque(body, torque);
        }
    }

    /// Finish the current step.
    pub fn end_step(&mut self) {
        if let Some(ref mut rec) = self.inner {
            rec.end_step();
        }
    }

    /// Number of steps recorded so far.
    pub fn recorded_steps(&self) -> usize {
        self.inner.as_ref().map_or(0, |r| r.recorded_steps())
    }

    /// Finish recording and return the record as a JSON string.
    ///
    /// After calling this, the recorder is exhausted and cannot be used further.
    pub fn finish_to_json(&mut self) -> PyResult<String> {
        let rec = self
            .inner
            .take()
            .ok_or_else(|| pyo3::exceptions::PyValueError::new_err("recorder already finished"))?;
        let record = rec.finish();
        record
            .to_json()
            .map_err(pyo3::exceptions::PyValueError::new_err)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PyReplayRecord
// ─────────────────────────────────────────────────────────────────────────────

/// A completed simulation recording.
#[pyclass(name = "ReplayRecord")]
pub struct PyReplayRecord {
    inner: ReplayRecord,
}

#[pymethods]
impl PyReplayRecord {
    /// Create a record from a JSON string.
    #[staticmethod]
    pub fn from_json(json: &str) -> PyResult<Self> {
        ReplayRecord::from_json(json)
            .map(|inner| Self { inner })
            .map_err(pyo3::exceptions::PyValueError::new_err)
    }

    /// Number of steps in this record.
    pub fn step_count(&self) -> usize {
        self.inner.step_count()
    }

    /// Total simulated time.
    pub fn total_sim_time(&self) -> f64 {
        self.inner.total_sim_time()
    }

    /// Total number of commands recorded.
    pub fn total_commands(&self) -> usize {
        self.inner.total_commands()
    }

    /// Summary string for the record.
    pub fn summary(&self) -> String {
        self.inner.summary()
    }

    /// Export the record as a JSON string.
    pub fn to_json(&self) -> PyResult<String> {
        self.inner
            .to_json()
            .map_err(pyo3::exceptions::PyValueError::new_err)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PySimReplayer
// ─────────────────────────────────────────────────────────────────────────────

/// Replays a simulation record step by step.
#[pyclass(name = "SimReplayer")]
pub struct PySimReplayer {
    inner: SimReplayer,
}

#[pymethods]
impl PySimReplayer {
    /// Create a replayer from a ReplayRecord.
    #[new]
    pub fn new(record: &PyReplayRecord) -> Self {
        Self {
            inner: SimReplayer::new(record.inner.clone()),
        }
    }

    /// Whether playback is complete.
    pub fn is_done(&self) -> bool {
        self.inner.is_done()
    }

    /// Total steps in the record.
    pub fn total_steps(&self) -> usize {
        self.inner.total_steps()
    }

    /// Steps remaining.
    pub fn steps_remaining(&self) -> usize {
        self.inner.steps_remaining()
    }

    /// Reset to the beginning.
    pub fn reset(&mut self) {
        self.inner.reset();
    }

    /// Advance to the next step and return the step data as JSON, or None if done.
    pub fn advance_json(&mut self) -> PyResult<Option<String>> {
        match self.inner.advance() {
            None => Ok(None),
            Some(step) => serde_json::to_string(step)
                .map(Some)
                .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string())),
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Module registration
// ─────────────────────────────────────────────────────────────────────────────

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PySimRecorder>()?;
    m.add_class::<PyReplayRecord>()?;
    m.add_class::<PySimReplayer>()?;
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sim_recorder_instantiation() {
        let mut rec = PySimRecorder::new(0);
        assert_eq!(rec.recorded_steps(), 0);
        rec.begin_step(1.0 / 60.0);
        rec.record_force(0, [0.0, -9.81, 0.0]);
        rec.end_step();
        assert_eq!(rec.recorded_steps(), 1);
    }

    #[test]
    fn test_record_round_trip_json() {
        let mut rec = PySimRecorder::new(0);
        rec.begin_step(1.0 / 60.0);
        rec.end_step();
        let json = rec.finish_to_json().expect("to_json failed");
        let record = PyReplayRecord::from_json(&json).expect("from_json failed");
        assert_eq!(record.step_count(), 1);
    }

    #[test]
    fn test_sim_replayer() {
        let mut rec = PySimRecorder::new(0);
        rec.begin_step(0.016);
        rec.end_step();
        let json = rec.finish_to_json().expect("to_json failed");
        let record = PyReplayRecord::from_json(&json).expect("from_json failed");
        let mut replayer = PySimReplayer::new(&record);
        assert!(!replayer.is_done());
        let step = replayer.advance_json().expect("advance failed");
        assert!(step.is_some());
        assert!(replayer.is_done());
    }
}
