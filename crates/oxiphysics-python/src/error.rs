// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Error types for the `oxiphysics-python` crate.
//!
//! All public API methods that can fail return \[`Result`T`].
//! Errors are designed to be easily surfaced as Python exceptions when
//! exposed via PyO3/pyo3 FFI, and serialize cleanly to JSON for interchange.
//!
//! # Overview
//!
//! - [`enum@Error`] — comprehensive enum covering all Python-bridge failure modes
//! - [`Result`T`\] — convenience alias
//! - Helper constructors and classification predicates on [`enum@Error`]

use serde::{Deserialize, Serialize};
use thiserror::Error;

// ---------------------------------------------------------------------------
// Core error enum
// ---------------------------------------------------------------------------

/// Comprehensive error type for the `oxiphysics-python` bridge.
///
/// Every failure mode has a dedicated variant to allow fine-grained handling
/// from Python user code and for clean mapping to `PyErr` via PyO3.
#[derive(Debug, Error, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", content = "data")]
pub enum Error {
    // --- Handle / identity errors ---
    /// A body or object handle was invalid.
    #[error("invalid handle: {0}")]
    InvalidHandle(u32),

    /// A rigid body was not found by handle.
    #[error("body not found: handle={0}")]
    BodyNotFound(u32),

    /// A collider was not found by handle.
    #[error("collider not found: handle={0}")]
    ColliderNotFound(u32),

    // --- Parameter / validation errors ---
    /// A named parameter was out of its valid range.
    #[error("invalid parameter '{name}': {message}")]
    InvalidParameter {
        /// Parameter name.
        name: String,
        /// Constraint description.
        message: String,
    },

    /// A mass value was non-positive.
    #[error("mass must be positive, got {0}")]
    InvalidMass(f64),

    /// A time step was non-positive.
    #[error("time step must be positive, got {0}")]
    InvalidTimeStep(f64),

    /// A geometry dimension was non-positive.
    #[error("dimension must be positive, got {0}")]
    InvalidDimension(f64),

    // --- Capacity / resource errors ---
    /// The world has reached its maximum body capacity.
    #[error("world at capacity: max {max} bodies")]
    CapacityExceeded {
        /// Maximum allowed bodies.
        max: usize,
    },

    // --- Simulation state errors ---
    /// The simulation has diverged (NaN / Inf detected).
    #[error("simulation diverged at step {step}")]
    SimulationDiverged {
        /// Step at which divergence was detected.
        step: u64,
    },

    /// The solver failed to converge.
    #[error("solver did not converge after {iterations} iterations")]
    SolverConvergenceFailed {
        /// Iterations attempted.
        iterations: u32,
    },

    /// A body was queried but it is sleeping.
    #[error("body {0} is sleeping")]
    BodySleeping(u32),

    // --- Serialization / IO errors ---
    /// JSON serialization or deserialization failed.
    #[error("serialization error: {0}")]
    Serialization(String),

    /// Snapshot validation failed (schema mismatch, missing fields, etc.).
    #[error("snapshot validation failed: {0}")]
    SnapshotValidation(String),

    // --- Python interop errors ---
    /// A Python argument had an unexpected type.
    #[error("type error: expected {expected}, got {got}")]
    TypeError {
        /// Expected Python type name.
        expected: String,
        /// Actual Python type name.
        got: String,
    },

    /// A required Python keyword argument was missing.
    #[error("missing argument: '{0}'")]
    MissingArgument(String),

    /// A Python list or array had the wrong length.
    #[error("wrong array length: expected {expected}, got {got}")]
    WrongArrayLength {
        /// Expected length.
        expected: usize,
        /// Actual length.
        got: usize,
    },

    // --- Generic ---
    /// Generic error with a free-form message.
    #[error("{0}")]
    General(String),
}

/// Result type alias for `oxiphysics-python` operations.
pub type Result<T> = std::result::Result<T, Error>;

impl From<Error> for pyo3::PyErr {
    fn from(e: Error) -> pyo3::PyErr {
        pyo3::exceptions::PyRuntimeError::new_err(e.to_string())
    }
}

// ---------------------------------------------------------------------------
// impl Error
// ---------------------------------------------------------------------------

impl Error {
    // --- Constructors -------------------------------------------------------

    /// Create an `InvalidParameter` error.
    pub fn invalid_param(name: impl Into<String>, message: impl Into<String>) -> Self {
        Self::InvalidParameter {
            name: name.into(),
            message: message.into(),
        }
    }

    /// Create a `General` error.
    pub fn general(msg: impl Into<String>) -> Self {
        Self::General(msg.into())
    }

    /// Create a `TypeError` error.
    pub fn type_error(expected: impl Into<String>, got: impl Into<String>) -> Self {
        Self::TypeError {
            expected: expected.into(),
            got: got.into(),
        }
    }

    /// Create a `WrongArrayLength` error.
    pub fn wrong_len(expected: usize, got: usize) -> Self {
        Self::WrongArrayLength { expected, got }
    }

    // --- Classification predicates ------------------------------------------

    /// Returns `true` if this error is handle-related.
    pub fn is_handle_error(&self) -> bool {
        matches!(
            self,
            Error::InvalidHandle(_) | Error::BodyNotFound(_) | Error::ColliderNotFound(_)
        )
    }

    /// Returns `true` if this error is a parameter validation error.
    pub fn is_parameter_error(&self) -> bool {
        matches!(
            self,
            Error::InvalidParameter { .. }
                | Error::InvalidMass(_)
                | Error::InvalidTimeStep(_)
                | Error::InvalidDimension(_)
        )
    }

    /// Returns `true` if this error is capacity-related.
    pub fn is_capacity_error(&self) -> bool {
        matches!(self, Error::CapacityExceeded { .. })
    }

    /// Returns `true` if this error indicates simulation instability.
    pub fn is_stability_error(&self) -> bool {
        matches!(
            self,
            Error::SimulationDiverged { .. } | Error::SolverConvergenceFailed { .. }
        )
    }

    /// Returns `true` if this is a serialization error.
    pub fn is_serialization_error(&self) -> bool {
        matches!(self, Error::Serialization(_) | Error::SnapshotValidation(_))
    }

    /// Returns `true` if this is a Python interop type error.
    pub fn is_type_error(&self) -> bool {
        matches!(
            self,
            Error::TypeError { .. } | Error::MissingArgument(_) | Error::WrongArrayLength { .. }
        )
    }

    // --- JSON / Python interop ----------------------------------------------

    /// Serialize this error to a JSON string for Python exception chaining.
    ///
    /// # Example
    ///
    /// ```no_run
    /// use oxiphysics_python::Error;
    ///
    /// let e = Error::InvalidTimeStep(-0.1);
    /// let json = e.to_json();
    /// assert!(json.contains("InvalidTimeStep"));
    /// ```
    pub fn to_json(&self) -> String {
        let variant_json =
            serde_json::to_string(self).unwrap_or_else(|_| "\"<serialization failed>\"".into());
        let message = self.to_string();
        format!(
            r#"{{"error":{variant_json},"message":{message_json}}}"#,
            variant_json = variant_json,
            message_json = serde_json::to_string(&message).unwrap_or_default(),
        )
    }

    /// Deserialize an `Error` from JSON.
    pub fn from_json(json: &str) -> std::result::Result<Self, String> {
        // Try direct variant JSON first
        let direct: std::result::Result<Error, _> = serde_json::from_str(json);
        if let Ok(e) = direct {
            return Ok(e);
        }
        // Try envelope with "error" field
        let v: serde_json::Value = serde_json::from_str(json).map_err(|e| e.to_string())?;
        let inner = v
            .get("error")
            .ok_or_else(|| "missing 'error' field".to_string())?;
        serde_json::from_value(inner.clone()).map_err(|e| e.to_string())
    }

    /// Return a Python exception class name that best maps to this error.
    ///
    /// Useful when calling PyO3's `PyErr::new::<PyXxx, _>()`.
    pub fn python_exception_class(&self) -> &'static str {
        match self {
            Error::TypeError { .. } | Error::WrongArrayLength { .. } => "TypeError",
            Error::MissingArgument(_) => "ValueError",
            Error::InvalidParameter { .. }
            | Error::InvalidMass(_)
            | Error::InvalidTimeStep(_)
            | Error::InvalidDimension(_) => "ValueError",
            Error::InvalidHandle(_) | Error::BodyNotFound(_) | Error::ColliderNotFound(_) => {
                "KeyError"
            }
            Error::CapacityExceeded { .. } => "MemoryError",
            Error::SimulationDiverged { .. } | Error::SolverConvergenceFailed { .. } => {
                "RuntimeError"
            }
            Error::Serialization(_) | Error::SnapshotValidation(_) => "ValueError",
            _ => "RuntimeError",
        }
    }

    /// Return a human-readable recovery hint for this error.
    pub fn recovery_hint(&self) -> String {
        match self {
            Error::InvalidTimeStep(_) => {
                "Use a positive dt such as 1/60 for a 60 Hz simulation.".into()
            }
            Error::InvalidMass(_) => "Mass must be strictly positive (e.g. 1.0).".into(),
            Error::BodyNotFound(h) => format!(
                "Body handle {} not found. Re-add the body or check it has not been removed.",
                h
            ),
            Error::CapacityExceeded { max } => format!(
                "World capacity of {} bodies reached. Remove unused bodies first.",
                max
            ),
            Error::SimulationDiverged { step } => format!(
                "Divergence detected at step {}. Reduce dt or applied forces.",
                step
            ),
            Error::TypeError { expected, got } => {
                format!("Expected a Python '{}' but received '{}'.", expected, got)
            }
            Error::WrongArrayLength { expected, got } => {
                format!("Array length should be {} but is {}.", expected, got)
            }
            _ => String::new(),
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ---- basic display tests -----------------------------------------------

    #[test]
    fn test_invalid_handle() {
        let e = Error::InvalidHandle(5);
        assert!(e.to_string().contains("5"));
        assert!(e.is_handle_error());
    }

    #[test]
    fn test_body_not_found() {
        let e = Error::BodyNotFound(10);
        assert!(e.to_string().contains("10"));
        assert!(e.is_handle_error());
    }

    #[test]
    fn test_collider_not_found() {
        let e = Error::ColliderNotFound(3);
        assert!(e.is_handle_error());
    }

    #[test]
    fn test_invalid_param() {
        let e = Error::invalid_param("mass", "must be positive");
        assert!(e.to_string().contains("mass"));
        assert!(e.to_string().contains("positive"));
        assert!(e.is_parameter_error());
    }

    #[test]
    fn test_invalid_mass() {
        let e = Error::InvalidMass(-1.0);
        assert!(e.is_parameter_error());
        assert!(e.to_string().contains("-1"));
    }

    #[test]
    fn test_invalid_time_step() {
        let e = Error::InvalidTimeStep(0.0);
        assert!(e.is_parameter_error());
    }

    #[test]
    fn test_invalid_dimension() {
        let e = Error::InvalidDimension(-0.5);
        assert!(e.is_parameter_error());
    }

    #[test]
    fn test_capacity_exceeded() {
        let e = Error::CapacityExceeded { max: 1000 };
        assert!(e.is_capacity_error());
        assert!(e.to_string().contains("1000"));
    }

    #[test]
    fn test_simulation_diverged() {
        let e = Error::SimulationDiverged { step: 99 };
        assert!(e.is_stability_error());
        assert!(e.to_string().contains("99"));
    }

    #[test]
    fn test_solver_convergence_failed() {
        let e = Error::SolverConvergenceFailed { iterations: 50 };
        assert!(e.is_stability_error());
        assert!(e.to_string().contains("50"));
    }

    #[test]
    fn test_body_sleeping() {
        let e = Error::BodySleeping(7);
        assert!(e.to_string().contains("7"));
    }

    #[test]
    fn test_serialization_error() {
        let e = Error::Serialization("unexpected eof".into());
        assert!(e.is_serialization_error());
    }

    #[test]
    fn test_snapshot_validation_error() {
        let e = Error::SnapshotValidation("missing version".into());
        assert!(e.is_serialization_error());
    }

    #[test]
    fn test_type_error() {
        let e = Error::type_error("list", "int");
        assert!(e.is_type_error());
        assert!(e.to_string().contains("list"));
        assert!(e.to_string().contains("int"));
    }

    #[test]
    fn test_missing_argument() {
        let e = Error::MissingArgument("mass".into());
        assert!(e.is_type_error());
        assert!(e.to_string().contains("mass"));
    }

    #[test]
    fn test_wrong_array_length() {
        let e = Error::wrong_len(3, 2);
        assert!(e.is_type_error());
        assert!(e.to_string().contains("3"));
        assert!(e.to_string().contains("2"));
    }

    #[test]
    fn test_general_error() {
        let e = Error::general("oops");
        assert!(e.to_string().contains("oops"));
    }

    // ---- clone / eq -------------------------------------------------------

    #[test]
    fn test_clone_eq() {
        let e1 = Error::InvalidHandle(42);
        let e2 = e1.clone();
        assert_eq!(e1, e2);
    }

    // ---- JSON roundtrip ---------------------------------------------------

    #[test]
    fn test_to_json_contains_type() {
        let e = Error::InvalidTimeStep(-0.01);
        let json = e.to_json();
        assert!(json.contains("InvalidTimeStep"), "json={}", json);
        assert!(json.contains("message"), "json={}", json);
    }

    #[test]
    fn test_from_json_direct() {
        let original = Error::BodyNotFound(7);
        let json = serde_json::to_string(&original).unwrap();
        let recovered = Error::from_json(&json).unwrap();
        assert_eq!(original, recovered);
    }

    #[test]
    fn test_from_json_envelope() {
        let e = Error::CapacityExceeded { max: 256 };
        let envelope = e.to_json();
        let recovered = Error::from_json(&envelope).unwrap();
        assert_eq!(recovered, e);
    }

    #[test]
    fn test_from_json_invalid() {
        assert!(Error::from_json("{bad json").is_err());
    }

    // ---- python_exception_class -------------------------------------------

    #[test]
    fn test_python_exception_class_value_error() {
        assert_eq!(
            Error::InvalidTimeStep(0.0).python_exception_class(),
            "ValueError"
        );
        assert_eq!(
            Error::InvalidMass(-1.0).python_exception_class(),
            "ValueError"
        );
        assert_eq!(
            Error::invalid_param("x", "y").python_exception_class(),
            "ValueError"
        );
    }

    #[test]
    fn test_python_exception_class_key_error() {
        assert_eq!(Error::BodyNotFound(1).python_exception_class(), "KeyError");
        assert_eq!(Error::InvalidHandle(0).python_exception_class(), "KeyError");
    }

    #[test]
    fn test_python_exception_class_type_error() {
        assert_eq!(
            Error::type_error("list", "int").python_exception_class(),
            "TypeError"
        );
        assert_eq!(Error::wrong_len(3, 1).python_exception_class(), "TypeError");
    }

    #[test]
    fn test_python_exception_class_runtime_error() {
        assert_eq!(
            Error::SimulationDiverged { step: 1 }.python_exception_class(),
            "RuntimeError"
        );
        assert_eq!(
            Error::SolverConvergenceFailed { iterations: 10 }.python_exception_class(),
            "RuntimeError"
        );
    }

    // ---- recovery hints ---------------------------------------------------

    #[test]
    fn test_recovery_hint_time_step() {
        let hint = Error::InvalidTimeStep(0.0).recovery_hint();
        assert!(!hint.is_empty());
    }

    #[test]
    fn test_recovery_hint_body_not_found() {
        let hint = Error::BodyNotFound(42).recovery_hint();
        assert!(hint.contains("42"));
    }

    #[test]
    fn test_recovery_hint_capacity() {
        let hint = Error::CapacityExceeded { max: 500 }.recovery_hint();
        assert!(hint.contains("500"));
    }

    #[test]
    fn test_recovery_hint_diverged() {
        let hint = Error::SimulationDiverged { step: 7 }.recovery_hint();
        assert!(hint.contains("7"));
    }

    #[test]
    fn test_recovery_hint_type_error() {
        let hint = Error::type_error("ndarray", "str").recovery_hint();
        assert!(hint.contains("ndarray"));
    }

    #[test]
    fn test_recovery_hint_wrong_len() {
        let hint = Error::wrong_len(3, 1).recovery_hint();
        assert!(hint.contains("3"));
        assert!(hint.contains("1"));
    }

    #[test]
    fn test_recovery_hint_general_empty() {
        // General errors have no specific hint
        let hint = Error::general("x").recovery_hint();
        // May or may not be empty, just ensure it does not panic
        let _ = hint;
    }
}
