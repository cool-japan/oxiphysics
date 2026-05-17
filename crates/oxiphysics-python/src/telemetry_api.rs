// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Python bindings for the oxiphysics telemetry module.
//!
//! Exposes per-step physics telemetry, rolling averages, and CSV export.

#![allow(missing_docs)]
#![allow(dead_code)]

use oxiphysics::telemetry::{PhysicsStats, TelemetrySession};
use pyo3::prelude::*;

#[cfg(feature = "numpy-bridge")]
use numpy::{IntoPyArray, PyArray1};
#[cfg(feature = "numpy-bridge")]
use pyo3::Bound;

// ─────────────────────────────────────────────────────────────────────────────
// PyTelemetrySession
// ─────────────────────────────────────────────────────────────────────────────

/// Records per-step physics statistics with rolling averages.
#[pyclass(name = "TelemetrySession")]
pub struct PyTelemetrySession {
    inner: TelemetrySession,
}

#[pymethods]
impl PyTelemetrySession {
    /// Create a new telemetry session keeping up to `capacity` entries.
    #[new]
    pub fn new(capacity: usize) -> Self {
        Self {
            inner: TelemetrySession::new(capacity),
        }
    }

    /// Push a stats entry for one step.
    ///
    /// Arguments: step, dt, body_count, sleeping_count, contact_count,
    ///            island_count, solve_iterations, broad_phase_pairs,
    ///            kinetic_energy, elapsed_ms.
    #[allow(clippy::too_many_arguments)]
    pub fn push(
        &mut self,
        step: u64,
        dt: f64,
        body_count: usize,
        sleeping_count: usize,
        contact_count: usize,
        island_count: usize,
        solve_iterations: usize,
        broad_phase_pairs: usize,
        kinetic_energy: f64,
        elapsed_ms: f64,
    ) {
        let stats = PhysicsStats {
            step,
            dt,
            body_count,
            sleeping_count,
            contact_count,
            island_count,
            solve_iterations,
            broad_phase_pairs,
            kinetic_energy,
            elapsed_ms,
        };
        self.inner.push(stats);
    }

    /// Number of recorded entries.
    pub fn len(&self) -> usize {
        self.inner.len()
    }

    /// `true` if no entries recorded.
    pub fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }

    /// Clear all entries.
    pub fn clear(&mut self) {
        self.inner.clear();
    }

    /// Latest stats as JSON, or None.
    pub fn latest_json(&self) -> PyResult<Option<String>> {
        match self.inner.latest() {
            None => Ok(None),
            Some(s) => serde_json::to_string(s)
                .map(Some)
                .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string())),
        }
    }

    /// Average stats as JSON, or None if no entries.
    pub fn average_json(&self) -> PyResult<Option<String>> {
        match self.inner.average() {
            None => Ok(None),
            Some(avg) => serde_json::to_string(&avg)
                .map(Some)
                .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string())),
        }
    }

    /// Peak kinetic energy across all recorded steps.
    pub fn peak_kinetic_energy(&self) -> Option<f64> {
        self.inner.peak_kinetic_energy()
    }

    /// Peak elapsed time in milliseconds.
    pub fn peak_elapsed_ms(&self) -> Option<f64> {
        self.inner.peak_elapsed_ms()
    }

    /// Average elapsed time in milliseconds.
    pub fn avg_elapsed_ms(&self) -> Option<f64> {
        self.inner.avg_elapsed_ms()
    }

    /// Energy trend (positive = increasing, negative = decreasing).
    pub fn energy_trend(&self) -> f64 {
        self.inner.energy_trend()
    }

    /// Export all recorded stats as CSV.
    pub fn to_csv(&self) -> String {
        self.inner.to_csv()
    }

    /// Extract a single numeric field from every recorded [`PhysicsStats`]
    /// entry and return the values as a 1-D `numpy.ndarray` of dtype `float64`.
    ///
    /// Recognised field names (case-sensitive):
    /// `"step"`, `"dt"`, `"body_count"`, `"sleeping_count"`, `"contact_count"`,
    /// `"island_count"`, `"solve_iterations"`, `"broad_phase_pairs"`,
    /// `"kinetic_energy"`, `"elapsed_ms"`.
    ///
    /// # Errors
    /// Returns `PyKeyError` for unrecognised field names.
    #[cfg(feature = "numpy-bridge")]
    pub fn time_series_to_numpy<'py>(
        &self,
        py: Python<'py>,
        field: &str,
    ) -> PyResult<Bound<'py, PyArray1<f64>>> {
        let extractor: fn(&PhysicsStats) -> f64 = match field {
            "step" => |s: &PhysicsStats| s.step as f64,
            "dt" => |s: &PhysicsStats| s.dt,
            "body_count" => |s: &PhysicsStats| s.body_count as f64,
            "sleeping_count" => |s: &PhysicsStats| s.sleeping_count as f64,
            "contact_count" => |s: &PhysicsStats| s.contact_count as f64,
            "island_count" => |s: &PhysicsStats| s.island_count as f64,
            "solve_iterations" => |s: &PhysicsStats| s.solve_iterations as f64,
            "broad_phase_pairs" => |s: &PhysicsStats| s.broad_phase_pairs as f64,
            "kinetic_energy" => |s: &PhysicsStats| s.kinetic_energy,
            "elapsed_ms" => |s: &PhysicsStats| s.elapsed_ms,
            other => {
                return Err(pyo3::exceptions::PyKeyError::new_err(format!(
                    "Unknown telemetry field: {other}. \
                     Valid fields: step, dt, body_count, sleeping_count, \
                     contact_count, island_count, solve_iterations, \
                     broad_phase_pairs, kinetic_energy, elapsed_ms"
                )));
            }
        };
        let series: Vec<f64> = self.inner.iter().map(extractor).collect();
        Ok(series.into_pyarray(py))
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Module registration
// ─────────────────────────────────────────────────────────────────────────────

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyTelemetrySession>()?;
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_telemetry_session_instantiation() {
        let mut session = PyTelemetrySession::new(100);
        assert!(session.is_empty());
        session.push(0, 0.016, 10, 2, 5, 3, 10, 20, 100.0, 1.5);
        assert_eq!(session.len(), 1);
        session.clear();
        assert!(session.is_empty());
    }

    #[test]
    fn test_telemetry_session_stats() {
        let mut session = PyTelemetrySession::new(100);
        session.push(0, 0.016, 5, 0, 2, 2, 4, 10, 50.0, 1.0);
        session.push(1, 0.016, 5, 0, 3, 2, 4, 11, 60.0, 2.0);
        let latest = session.latest_json().expect("latest_json failed");
        assert!(latest.is_some());
        let avg = session.average_json().expect("average_json failed");
        assert!(avg.is_some());
        let peak = session.peak_elapsed_ms();
        assert!(peak.is_some());
        assert!((peak.expect("peak should not be None") - 2.0).abs() < 1e-9);
    }

    #[test]
    fn test_telemetry_csv_export() {
        let mut session = PyTelemetrySession::new(10);
        session.push(0, 0.016, 3, 0, 1, 1, 2, 5, 10.0, 1.0);
        let csv = session.to_csv();
        assert!(csv.contains("step") || csv.contains("0"));
    }
}
