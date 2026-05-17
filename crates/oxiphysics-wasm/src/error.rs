// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Error types for the `oxiphysics-wasm` crate.
//!
//! All public API methods that can fail return `Result<T, Error>`.
//! Errors are designed to be easily converted to JavaScript exceptions
//! when exposed via wasm-bindgen.
//!
//! # Architecture
//!
//! - [`enum@Error`] — the main error enum covering every physics subsystem
//! - [`ErrorContext`] — optional stack of contextual messages for richer diagnostics
//! - [`ErrorWithContext`] — wraps an [`enum@Error`] with an [`ErrorContext`]
//! - [`RecoverySuggestion`] — machine-readable hint for how to recover from an error
//! - [`ErrorLog`] — structured log of errors for debugging WASM simulations
//! - JSON serialization: [`Error::to_json`] and [`Error::from_json`]
//!
//! ## JS interop
//!
//! ```no_run
//! use oxiphysics_wasm::error::Error;
//!
//! let e = Error::InvalidTimeStep(-0.01);
//! let json = e.to_json();
//! assert!(json.contains("InvalidTimeStep"));
//! assert!(json.contains("suggestion"));
//! ```

use serde::{Deserialize, Serialize};
use thiserror::Error;

// ---------------------------------------------------------------------------
// Core error enum
// ---------------------------------------------------------------------------

/// Comprehensive error type for `oxiphysics-wasm` operations.
///
/// Every physics subsystem has at least one dedicated variant, making it easy
/// to handle errors programmatically from JavaScript.
#[derive(Debug, Error, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", content = "data")]
pub enum Error {
    // --- Handle / identity errors ---
    /// A body or collider handle was invalid or referred to a freed slot.
    #[error("invalid handle: {0}")]
    InvalidHandle(u32),

    /// A required body was not found by its handle.
    #[error("body not found: handle={0}")]
    BodyNotFound(u32),

    /// A required collider was not found by its handle.
    #[error("collider not found: handle={0}")]
    ColliderNotFound(u32),

    /// A constraint or joint was not found by its handle.
    #[error("constraint not found: handle={0}")]
    ConstraintNotFound(u32),

    // --- Parameter errors ---
    /// A simulation parameter was out of the valid range.
    #[error("invalid parameter '{name}': {message}")]
    InvalidParameter {
        /// Name of the offending parameter.
        name: String,
        /// Description of the constraint violated.
        message: String,
    },

    /// The time step `dt` was non-positive.
    #[error("time step must be positive, got {0}")]
    InvalidTimeStep(f64),

    /// A mass value was non-positive.
    #[error("mass must be positive, got {0}")]
    InvalidMass(f64),

    /// A friction coefficient was negative.
    #[error("friction must be >= 0, got {0}")]
    InvalidFriction(f64),

    /// A restitution coefficient was outside `[0, 1]`.
    #[error("restitution must be in [0, 1], got {0}")]
    InvalidRestitution(f64),

    /// A geometry dimension (radius, half-extent, etc.) was non-positive.
    #[error("geometry dimension must be positive, got {0}")]
    InvalidDimension(f64),

    // --- Capacity / resource errors ---
    /// The engine has reached its maximum body capacity.
    #[error("engine at capacity: cannot add more than {max} bodies")]
    CapacityExceeded {
        /// Maximum number of bodies allowed.
        max: usize,
    },

    /// The maximum number of colliders per body was exceeded.
    #[error("collider limit exceeded: body {handle} already has {current} colliders (max {max})")]
    ColliderLimitExceeded {
        /// Body handle.
        handle: u32,
        /// Current collider count.
        current: usize,
        /// Maximum allowed.
        max: usize,
    },

    /// The constraint count limit was reached.
    #[error("constraint limit exceeded: max {max}")]
    ConstraintLimitExceeded {
        /// Maximum allowed.
        max: usize,
    },

    // --- Simulation state errors ---
    /// An operation was attempted on a sleeping body that requires it to be awake.
    #[error("body {0} is sleeping; wake it first")]
    BodySleeping(u32),

    /// Two bodies that were expected to be different turned out to be the same.
    #[error("body self-collision: handles {0} and {1} refer to the same body")]
    SelfCollision(u32, u32),

    /// The simulation has diverged (NaN or Inf detected in state).
    #[error("simulation diverged: NaN or Inf detected after step {step}")]
    SimulationDiverged {
        /// The simulation step at which divergence was detected.
        step: u64,
    },

    /// The solver failed to converge within the maximum number of iterations.
    #[error("solver did not converge after {iterations} iterations (tolerance={tolerance})")]
    SolverConvergenceFailed {
        /// Number of iterations attempted.
        iterations: u32,
        /// Convergence tolerance that was not met.
        tolerance: f64,
    },

    // --- Collision detection errors ---
    /// The broad phase returned an unexpected overlap pair.
    #[error("broad phase error: invalid overlap pair ({0}, {1})")]
    BroadPhaseError(u32, u32),

    /// A raycast hit an internal inconsistency.
    #[error("raycast error: {0}")]
    RaycastError(String),

    /// The requested collision shape type is not supported.
    #[error("unsupported collision shape: {0}")]
    UnsupportedShape(String),

    // --- Rigid body dynamics errors ---
    /// An applied impulse had zero magnitude.
    #[error("zero impulse applied to body {0}: impulse must be non-zero")]
    ZeroImpulse(u32),

    /// An integration error occurred in the rigid body stepper.
    #[error("integration error for body {handle}: {message}")]
    IntegrationError {
        /// Body handle.
        handle: u32,
        /// Error description.
        message: String,
    },

    // --- Serialization / IO errors ---
    /// JSON serialization or deserialization failed.
    #[error("serialization error: {0}")]
    Serialization(String),

    /// A configuration JSON blob could not be deserialized into the expected type.
    #[error("config deserialization failed: {0}")]
    ConfigDeserialization(String),

    // --- Generic ---
    /// A generic internal error with a message.
    #[error("{0}")]
    General(String),

    /// An internal assertion failed — should never reach user code.
    #[error("internal assertion failed: {0}")]
    InternalAssertion(String),
}

/// Result type alias for `oxiphysics-wasm` operations.
pub type Result<T> = std::result::Result<T, Error>;

impl From<Error> for wasm_bindgen::JsValue {
    fn from(e: Error) -> wasm_bindgen::JsValue {
        wasm_bindgen::JsValue::from_str(&e.to_string())
    }
}

// ---------------------------------------------------------------------------
// Core impl
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

    /// Create a `General` error from any displayable value.
    pub fn general(msg: impl Into<String>) -> Self {
        Self::General(msg.into())
    }

    /// Create an `InternalAssertion` error.
    pub fn internal(msg: impl Into<String>) -> Self {
        Self::InternalAssertion(msg.into())
    }

    /// Create an `IntegrationError` for a specific body.
    pub fn integration(handle: u32, message: impl Into<String>) -> Self {
        Self::IntegrationError {
            handle,
            message: message.into(),
        }
    }

    // --- Classification predicates ------------------------------------------

    /// Returns `true` if this error was caused by an invalid or missing handle.
    pub fn is_invalid_handle(&self) -> bool {
        matches!(
            self,
            Error::InvalidHandle(_)
                | Error::BodyNotFound(_)
                | Error::ColliderNotFound(_)
                | Error::ConstraintNotFound(_)
        )
    }

    /// Returns `true` if this error is a capacity error.
    pub fn is_capacity_error(&self) -> bool {
        matches!(
            self,
            Error::CapacityExceeded { .. }
                | Error::ColliderLimitExceeded { .. }
                | Error::ConstraintLimitExceeded { .. }
        )
    }

    /// Returns `true` if this error relates to an invalid parameter value.
    pub fn is_parameter_error(&self) -> bool {
        matches!(
            self,
            Error::InvalidParameter { .. }
                | Error::InvalidTimeStep(_)
                | Error::InvalidMass(_)
                | Error::InvalidFriction(_)
                | Error::InvalidRestitution(_)
                | Error::InvalidDimension(_)
        )
    }

    /// Returns `true` if this error indicates simulation instability.
    pub fn is_stability_error(&self) -> bool {
        matches!(
            self,
            Error::SimulationDiverged { .. } | Error::SolverConvergenceFailed { .. }
        )
    }

    /// Returns `true` if this error is a collision detection error.
    pub fn is_collision_error(&self) -> bool {
        matches!(
            self,
            Error::BroadPhaseError(_, _)
                | Error::RaycastError(_)
                | Error::UnsupportedShape(_)
                | Error::SelfCollision(_, _)
        )
    }

    /// Returns `true` if this is a serialization-related error.
    pub fn is_serialization_error(&self) -> bool {
        matches!(
            self,
            Error::Serialization(_) | Error::ConfigDeserialization(_)
        )
    }

    // --- JS interop ---------------------------------------------------------

    /// Convert this error to a JavaScript-friendly string (for wasm-bindgen usage).
    pub fn to_js_string(&self) -> String {
        format!("OxiPhysics error: {}", self)
    }

    /// Serialize this error to a JSON string for JavaScript consumption.
    ///
    /// The JSON object includes:
    /// - `"type"`: the variant name
    /// - `"data"`: the variant payload
    /// - `"message"`: the human-readable error string
    /// - `"suggestion"`: a recovery suggestion (may be empty)
    ///
    /// # Example
    ///
    /// ```no_run
    /// use oxiphysics_wasm::error::Error;
    ///
    /// let e = Error::InvalidTimeStep(-0.1);
    /// let json = e.to_json();
    /// assert!(json.contains("InvalidTimeStep"));
    /// assert!(json.contains("message"));
    /// ```
    pub fn to_json(&self) -> String {
        let suggestion = self.recovery_suggestion();
        // Build a JSON object manually to avoid a serde dependency on the suggestion
        // while still including all fields.
        let variant_json =
            serde_json::to_string(self).unwrap_or_else(|_| "\"<serialization failed>\"".into());
        // variant_json is already a complete JSON object due to the serde tag
        // but we want to inject extra fields. Easiest: wrap in an envelope.
        let message = self.to_string();
        format!(
            r#"{{"error":{variant_json},"message":{message_json},"suggestion":{sugg_json}}}"#,
            variant_json = variant_json,
            message_json = serde_json::to_string(&message).unwrap_or_default(),
            sugg_json = serde_json::to_string(&suggestion).unwrap_or_default(),
        )
    }

    /// Deserialize an `Error` from a JSON string produced by [`Error::to_json`].
    ///
    /// Expects the `"error"` field to contain the serialized variant.
    pub fn from_json(json: &str) -> std::result::Result<Self, String> {
        // Try to parse as an envelope first, then fall back to direct.
        let direct: std::result::Result<Error, _> = serde_json::from_str(json);
        if let Ok(e) = direct {
            return Ok(e);
        }
        // Try envelope
        let v: serde_json::Value = serde_json::from_str(json).map_err(|e| e.to_string())?;
        let inner = v
            .get("error")
            .ok_or_else(|| "missing 'error' field".to_string())?;
        serde_json::from_value(inner.clone()).map_err(|e| e.to_string())
    }

    // --- Recovery suggestions -----------------------------------------------

    /// Return a human-readable recovery suggestion for this error.
    ///
    /// Returns an empty string if no specific suggestion is available.
    pub fn recovery_suggestion(&self) -> String {
        match self {
            Error::InvalidHandle(h) => format!(
                "Handle {} is invalid. Ensure you are using handles returned by add_dynamic_body() \
                 or add_static_body() and that the body has not been removed.",
                h
            ),
            Error::BodyNotFound(h) => format!(
                "Body with handle {} was not found. It may have been removed. \
                 Re-add it with add_dynamic_body().",
                h
            ),
            Error::ColliderNotFound(h) => format!(
                "Collider with handle {} was not found. \
                 Re-add it with add_sphere_collider() or add_plane_collider().",
                h
            ),
            Error::ConstraintNotFound(h) => format!(
                "Constraint handle {} is invalid. Use the handle returned by add_constraint().",
                h
            ),
            Error::InvalidTimeStep(_) => "Use a positive dt (e.g. 1/60 for 60 Hz simulation). \
                 Values in range (0, 0.1] are recommended."
                .to_string(),
            Error::InvalidMass(_) => {
                "Mass must be strictly positive. Use a value > 0 (e.g. 1.0 kg).".to_string()
            }
            Error::InvalidFriction(_) => {
                "Friction coefficient must be >= 0. Typical values: 0.0 (frictionless) to 1.0."
                    .to_string()
            }
            Error::InvalidRestitution(_) => {
                "Restitution (bounciness) must be in the range [0, 1]. \
                 0 = perfectly inelastic, 1 = perfectly elastic."
                    .to_string()
            }
            Error::InvalidDimension(_) => {
                "Geometry dimensions (radius, half-extents) must be strictly positive.".to_string()
            }
            Error::CapacityExceeded { max } => format!(
                "The simulation has reached its body limit of {}. \
                 Remove unused bodies with remove_body() or increase the engine capacity.",
                max
            ),
            Error::ColliderLimitExceeded { max, .. } => format!(
                "A body already has the maximum number of colliders ({}). \
                 Remove an existing collider before adding a new one.",
                max
            ),
            Error::BodySleeping(h) => format!(
                "Body {} is sleeping. Call wake_body({}) before applying forces or querying velocity.",
                h, h
            ),
            Error::SelfCollision(a, b) => format!(
                "Handles {} and {} refer to the same body. Pass two distinct body handles.",
                a, b
            ),
            Error::SimulationDiverged { step } => format!(
                "The simulation became numerically unstable at step {}. \
                 Try reducing the time step dt, increasing solver iterations, \
                 or reducing applied forces.",
                step
            ),
            Error::SolverConvergenceFailed { iterations, .. } => format!(
                "The constraint solver did not converge after {} iterations. \
                 Increase max_solver_iterations or relax the constraint tolerances.",
                iterations
            ),
            Error::UnsupportedShape(s) => format!(
                "The collision shape '{}' is not supported. \
                 Supported shapes: Sphere, Plane, Box, Capsule, Cylinder.",
                s
            ),
            Error::Serialization(_) => {
                "Check that the JSON input is well-formed and matches the expected schema."
                    .to_string()
            }
            Error::ConfigDeserialization(_) => {
                "Verify the simulation config JSON matches SimulationConfig schema.".to_string()
            }
            _ => String::new(),
        }
    }
}

// ---------------------------------------------------------------------------
// ErrorContext — contextual stack for richer diagnostics
// ---------------------------------------------------------------------------

/// A stack of contextual messages attached to an error.
///
/// Each entry is a string describing what operation was in progress when
/// the error was produced or propagated. The last entry is the outermost context.
///
/// # Example
///
/// ```no_run
/// use oxiphysics_wasm::error::{Error, ErrorContext};
///
/// let ctx = ErrorContext::new()
///     .push("stepping simulation")
///     .push("resolving contact between bodies 3 and 7");
/// let ew = ctx.wrap(Error::SolverConvergenceFailed { iterations: 50, tolerance: 1e-6 });
/// assert_eq!(ew.context.frames.len(), 2);
/// ```
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct ErrorContext {
    /// Stack frames, from innermost (index 0) to outermost (last).
    pub frames: Vec<String>,
}

impl ErrorContext {
    /// Create an empty context.
    pub fn new() -> Self {
        Self::default()
    }

    /// Push a new context frame and return `self` for chaining.
    pub fn push(mut self, frame: impl Into<String>) -> Self {
        self.frames.push(frame.into());
        self
    }

    /// Add a frame in place.
    pub fn add_frame(&mut self, frame: impl Into<String>) {
        self.frames.push(frame.into());
    }

    /// Return `true` if no frames have been pushed.
    pub fn is_empty(&self) -> bool {
        self.frames.is_empty()
    }

    /// Wrap an [`enum@Error`] with this context.
    pub fn wrap(self, error: Error) -> ErrorWithContext {
        ErrorWithContext {
            error,
            context: self,
        }
    }

    /// Format the context stack as a human-readable string.
    pub fn format(&self) -> String {
        self.frames
            .iter()
            .rev()
            .enumerate()
            .map(|(i, f)| format!("  [{}] {}", i, f))
            .collect::<Vec<_>>()
            .join("\n")
    }
}

// ---------------------------------------------------------------------------
// ErrorWithContext
// ---------------------------------------------------------------------------

/// An [`enum@Error`] paired with an [`ErrorContext`] for richer diagnostics.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ErrorWithContext {
    /// The underlying error.
    pub error: Error,
    /// Contextual stack accumulated while the error propagated.
    pub context: ErrorContext,
}

impl ErrorWithContext {
    /// Create a new `ErrorWithContext` with an empty context.
    pub fn new(error: Error) -> Self {
        Self {
            error,
            context: ErrorContext::new(),
        }
    }

    /// Add a context frame describing where the error was caught.
    pub fn with_context(mut self, frame: impl Into<String>) -> Self {
        self.context.add_frame(frame);
        self
    }

    /// Return the human-readable error message including context frames.
    pub fn full_message(&self) -> String {
        if self.context.is_empty() {
            self.error.to_string()
        } else {
            format!("{}\nContext:\n{}", self.error, self.context.format())
        }
    }

    /// Serialize to JSON including the context stack.
    pub fn to_json(&self) -> String {
        serde_json::to_string(self)
            .unwrap_or_else(|e| format!("{{\"error\":\"serialization failed: {}\"}}", e))
    }
}

impl std::fmt::Display for ErrorWithContext {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.full_message())
    }
}

// ---------------------------------------------------------------------------
// RecoverySuggestion — machine-readable recovery hints
// ---------------------------------------------------------------------------

/// Category of a recovery suggestion.
#[allow(dead_code)]
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum SuggestionKind {
    /// Reduce the simulation time step.
    ReduceTimeStep,
    /// Increase solver iterations.
    IncreaseSolverIterations,
    /// Check and fix the parameter value.
    FixParameter,
    /// Re-add or look up the handle.
    ResolveHandle,
    /// Remove or adjust resource usage.
    ReduceResourceUsage,
    /// No specific suggestion.
    General,
}

/// A structured recovery suggestion for a physics error.
#[allow(dead_code)]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RecoverySuggestion {
    /// Machine-readable category.
    pub kind: SuggestionKind,
    /// Human-readable description.
    pub description: String,
    /// Optional documentation URL.
    pub doc_url: Option<String>,
}

impl RecoverySuggestion {
    /// Create a new `RecoverySuggestion`.
    pub fn new(kind: SuggestionKind, description: impl Into<String>) -> Self {
        Self {
            kind,
            description: description.into(),
            doc_url: None,
        }
    }

    /// Attach a documentation URL.
    pub fn with_doc_url(mut self, url: impl Into<String>) -> Self {
        self.doc_url = Some(url.into());
        self
    }
}

/// Return a structured [`RecoverySuggestion`] for the given error.
#[allow(dead_code)]
pub fn structured_suggestion(error: &Error) -> RecoverySuggestion {
    match error {
        Error::InvalidTimeStep(_) | Error::SimulationDiverged { .. } => {
            RecoverySuggestion::new(SuggestionKind::ReduceTimeStep, error.recovery_suggestion())
        }
        Error::SolverConvergenceFailed { .. } => RecoverySuggestion::new(
            SuggestionKind::IncreaseSolverIterations,
            error.recovery_suggestion(),
        ),
        Error::InvalidParameter { .. }
        | Error::InvalidMass(_)
        | Error::InvalidFriction(_)
        | Error::InvalidRestitution(_)
        | Error::InvalidDimension(_) => {
            RecoverySuggestion::new(SuggestionKind::FixParameter, error.recovery_suggestion())
        }
        Error::InvalidHandle(_)
        | Error::BodyNotFound(_)
        | Error::ColliderNotFound(_)
        | Error::ConstraintNotFound(_) => {
            RecoverySuggestion::new(SuggestionKind::ResolveHandle, error.recovery_suggestion())
        }
        Error::CapacityExceeded { .. }
        | Error::ColliderLimitExceeded { .. }
        | Error::ConstraintLimitExceeded { .. } => RecoverySuggestion::new(
            SuggestionKind::ReduceResourceUsage,
            error.recovery_suggestion(),
        ),
        _ => RecoverySuggestion::new(SuggestionKind::General, error.recovery_suggestion()),
    }
}

// ---------------------------------------------------------------------------
// ErrorLog — structured error logging
// ---------------------------------------------------------------------------

/// A severity level for logged errors.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub enum LogLevel {
    /// Informational — not really an error.
    Info,
    /// Warning — something may be wrong.
    Warning,
    /// Error — an operation failed but simulation can continue.
    Error,
    /// Fatal — simulation state is invalid.
    Fatal,
}

impl std::fmt::Display for LogLevel {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            LogLevel::Info => write!(f, "INFO"),
            LogLevel::Warning => write!(f, "WARN"),
            LogLevel::Error => write!(f, "ERROR"),
            LogLevel::Fatal => write!(f, "FATAL"),
        }
    }
}

/// A single entry in the [`ErrorLog`].
#[allow(dead_code)]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ErrorLogEntry {
    /// Simulation step at which the error occurred.
    pub step: u64,
    /// Log level.
    pub level: LogLevel,
    /// The error that was logged.
    pub error: Error,
    /// Optional context stack.
    pub context: ErrorContext,
}

impl ErrorLogEntry {
    /// Format this entry as a single log line.
    pub fn format_line(&self) -> String {
        if self.context.is_empty() {
            format!("[{}] step={} {}", self.level, self.step, self.error)
        } else {
            format!(
                "[{}] step={} {}\n{}",
                self.level,
                self.step,
                self.error,
                self.context.format()
            )
        }
    }
}

/// A structured log of physics errors accumulated during a simulation run.
///
/// Useful for post-mortem analysis of failed simulations or for surfacing
/// warnings to JavaScript without immediately throwing exceptions.
///
/// # Example
///
/// ```no_run
/// use oxiphysics_wasm::error::{Error, ErrorLog, LogLevel};
///
/// let mut log = ErrorLog::new();
/// log.record(0, LogLevel::Warning, Error::BodySleeping(3), None);
/// log.record(5, LogLevel::Error, Error::SimulationDiverged { step: 5 }, None);
/// assert_eq!(log.len(), 2);
/// assert!(log.has_errors());
/// let json = log.to_json();
/// assert!(json.contains("SimulationDiverged"));
/// ```
#[allow(dead_code)]
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct ErrorLog {
    entries: Vec<ErrorLogEntry>,
    /// Maximum number of entries to retain (0 = unlimited).
    pub max_entries: usize,
}

#[allow(dead_code)]
impl ErrorLog {
    /// Create a new empty error log (unlimited capacity).
    pub fn new() -> Self {
        Self::default()
    }

    /// Create a log with a maximum capacity.
    pub fn with_capacity(max_entries: usize) -> Self {
        Self {
            entries: Vec::new(),
            max_entries,
        }
    }

    /// Record an error at the given simulation step and level.
    pub fn record(
        &mut self,
        step: u64,
        level: LogLevel,
        error: Error,
        context: Option<ErrorContext>,
    ) {
        if self.max_entries > 0 && self.entries.len() >= self.max_entries {
            self.entries.remove(0);
        }
        self.entries.push(ErrorLogEntry {
            step,
            level,
            error,
            context: context.unwrap_or_default(),
        });
    }

    /// Return the number of entries.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Returns `true` if the log is empty.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Returns `true` if any entry has level `Error` or higher.
    pub fn has_errors(&self) -> bool {
        self.entries.iter().any(|e| e.level >= LogLevel::Error)
    }

    /// Returns `true` if any entry has level `Fatal`.
    pub fn has_fatal(&self) -> bool {
        self.entries.iter().any(|e| e.level == LogLevel::Fatal)
    }

    /// Return all entries at or above the given level.
    pub fn entries_at_level(&self, min_level: LogLevel) -> Vec<&ErrorLogEntry> {
        self.entries
            .iter()
            .filter(|e| e.level >= min_level)
            .collect()
    }

    /// Return all entries for a specific simulation step.
    pub fn entries_at_step(&self, step: u64) -> Vec<&ErrorLogEntry> {
        self.entries.iter().filter(|e| e.step == step).collect()
    }

    /// Clear all entries.
    pub fn clear(&mut self) {
        self.entries.clear();
    }

    /// Format all entries as a multi-line log string.
    pub fn format(&self) -> String {
        self.entries
            .iter()
            .map(|e| e.format_line())
            .collect::<Vec<_>>()
            .join("\n")
    }

    /// Serialize the log to JSON.
    pub fn to_json(&self) -> String {
        serde_json::to_string(self)
            .unwrap_or_else(|e| format!("{{\"error\":\"log serialization failed: {}\"}}", e))
    }

    /// Return a reference to all entries.
    pub fn entries(&self) -> &[ErrorLogEntry] {
        &self.entries
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ---- basic variant tests -----------------------------------------------

    #[test]
    fn test_invalid_handle_display() {
        let e = Error::InvalidHandle(42);
        assert!(e.to_string().contains("42"));
    }

    #[test]
    fn test_invalid_parameter_display() {
        let e = Error::invalid_param("mass", "must be positive");
        assert!(e.to_string().contains("mass"));
        assert!(e.to_string().contains("positive"));
    }

    #[test]
    fn test_general_error() {
        let e = Error::general("something went wrong");
        assert!(e.to_string().contains("something went wrong"));
    }

    #[test]
    fn test_is_invalid_handle() {
        assert!(Error::InvalidHandle(0).is_invalid_handle());
        assert!(Error::BodyNotFound(1).is_invalid_handle());
        assert!(Error::ColliderNotFound(2).is_invalid_handle());
        assert!(Error::ConstraintNotFound(3).is_invalid_handle());
        assert!(!Error::General("x".into()).is_invalid_handle());
    }

    #[test]
    fn test_capacity_error() {
        let e = Error::CapacityExceeded { max: 1000 };
        assert!(e.is_capacity_error());
        assert!(e.to_string().contains("1000"));
    }

    #[test]
    fn test_collider_limit_exceeded() {
        let e = Error::ColliderLimitExceeded {
            handle: 5,
            current: 4,
            max: 4,
        };
        assert!(e.is_capacity_error());
        assert!(e.to_string().contains("5"));
    }

    #[test]
    fn test_constraint_limit_exceeded() {
        let e = Error::ConstraintLimitExceeded { max: 256 };
        assert!(e.is_capacity_error());
        assert!(e.to_string().contains("256"));
    }

    #[test]
    fn test_invalid_time_step() {
        let e = Error::InvalidTimeStep(-0.1);
        assert!(e.to_string().contains("-0.1"));
        assert!(e.is_parameter_error());
    }

    #[test]
    fn test_invalid_mass() {
        let e = Error::InvalidMass(0.0);
        assert!(e.is_parameter_error());
        assert!(e.to_string().contains("0"));
    }

    #[test]
    fn test_invalid_friction() {
        let e = Error::InvalidFriction(-1.0);
        assert!(e.is_parameter_error());
    }

    #[test]
    fn test_invalid_restitution() {
        let e = Error::InvalidRestitution(1.5);
        assert!(e.is_parameter_error());
    }

    #[test]
    fn test_invalid_dimension() {
        let e = Error::InvalidDimension(-0.5);
        assert!(e.is_parameter_error());
    }

    #[test]
    fn test_self_collision_error() {
        let e = Error::SelfCollision(3, 3);
        assert!(e.to_string().contains("3"));
        assert!(e.is_collision_error());
    }

    #[test]
    fn test_simulation_diverged() {
        let e = Error::SimulationDiverged { step: 42 };
        assert!(e.is_stability_error());
        assert!(e.to_string().contains("42"));
    }

    #[test]
    fn test_solver_convergence_failed() {
        let e = Error::SolverConvergenceFailed {
            iterations: 100,
            tolerance: 1e-6,
        };
        assert!(e.is_stability_error());
        assert!(e.to_string().contains("100"));
    }

    #[test]
    fn test_broad_phase_error() {
        let e = Error::BroadPhaseError(1, 2);
        assert!(e.is_collision_error());
    }

    #[test]
    fn test_raycast_error() {
        let e = Error::RaycastError("origin outside world".into());
        assert!(e.is_collision_error());
    }

    #[test]
    fn test_unsupported_shape() {
        let e = Error::UnsupportedShape("Mesh".into());
        assert!(e.is_collision_error());
    }

    #[test]
    fn test_integration_error() {
        let e = Error::integration(7, "NaN velocity");
        assert!(e.to_string().contains("NaN velocity"));
    }

    #[test]
    fn test_serialization_error() {
        let e = Error::Serialization("unexpected eof".into());
        assert!(e.is_serialization_error());
    }

    #[test]
    fn test_config_deserialization_error() {
        let e = Error::ConfigDeserialization("missing field 'gravity'".into());
        assert!(e.is_serialization_error());
    }

    #[test]
    fn test_to_js_string() {
        let e = Error::InvalidHandle(7);
        let js = e.to_js_string();
        assert!(js.starts_with("OxiPhysics error:"));
    }

    #[test]
    fn test_error_clone_eq() {
        let e1 = Error::InvalidHandle(10);
        let e2 = e1.clone();
        assert_eq!(e1, e2);
    }

    // ---- JSON serialization tests -----------------------------------------

    #[test]
    fn test_to_json_contains_type() {
        let e = Error::InvalidTimeStep(-0.01);
        let json = e.to_json();
        assert!(json.contains("InvalidTimeStep"), "json={}", json);
        assert!(json.contains("message"), "json={}", json);
        assert!(json.contains("suggestion"), "json={}", json);
    }

    #[test]
    fn test_to_json_general() {
        let e = Error::General("oops".into());
        let json = e.to_json();
        assert!(
            json.contains("General") || json.contains("oops"),
            "json={}",
            json
        );
    }

    #[test]
    fn test_from_json_roundtrip() {
        let original = Error::CapacityExceeded { max: 512 };
        // Serialize via serde directly for roundtrip
        let direct_json = serde_json::to_string(&original).unwrap();
        let recovered = Error::from_json(&direct_json).unwrap();
        assert_eq!(original, recovered);
    }

    #[test]
    fn test_from_json_envelope() {
        let e = Error::BodySleeping(5);
        let envelope_json = e.to_json();
        let recovered = Error::from_json(&envelope_json).unwrap();
        assert_eq!(recovered, e);
    }

    #[test]
    fn test_from_json_invalid() {
        assert!(Error::from_json("{not valid json}").is_err());
    }

    // ---- Recovery suggestion tests ----------------------------------------

    #[test]
    fn test_recovery_suggestion_time_step() {
        let e = Error::InvalidTimeStep(0.0);
        let s = e.recovery_suggestion();
        assert!(!s.is_empty());
        assert!(s.contains("positive") || s.contains("dt") || s.contains("Hz"));
    }

    #[test]
    fn test_recovery_suggestion_handle() {
        let e = Error::BodyNotFound(42);
        let s = e.recovery_suggestion();
        assert!(s.contains("42"));
    }

    #[test]
    fn test_recovery_suggestion_capacity() {
        let e = Error::CapacityExceeded { max: 100 };
        let s = e.recovery_suggestion();
        assert!(s.contains("100"));
    }

    #[test]
    fn test_recovery_suggestion_diverged() {
        let e = Error::SimulationDiverged { step: 10 };
        let s = e.recovery_suggestion();
        assert!(s.contains("10"));
    }

    #[test]
    fn test_structured_suggestion_kinds() {
        assert_eq!(
            structured_suggestion(&Error::InvalidTimeStep(0.0)).kind,
            SuggestionKind::ReduceTimeStep
        );
        assert_eq!(
            structured_suggestion(&Error::SolverConvergenceFailed {
                iterations: 50,
                tolerance: 1e-6
            })
            .kind,
            SuggestionKind::IncreaseSolverIterations
        );
        assert_eq!(
            structured_suggestion(&Error::InvalidMass(-1.0)).kind,
            SuggestionKind::FixParameter
        );
        assert_eq!(
            structured_suggestion(&Error::BodyNotFound(0)).kind,
            SuggestionKind::ResolveHandle
        );
        assert_eq!(
            structured_suggestion(&Error::CapacityExceeded { max: 10 }).kind,
            SuggestionKind::ReduceResourceUsage
        );
    }

    // ---- ErrorContext tests ------------------------------------------------

    #[test]
    fn test_error_context_push() {
        let ctx = ErrorContext::new().push("outer op").push("inner op");
        assert_eq!(ctx.frames.len(), 2);
        assert_eq!(ctx.frames[0], "outer op");
        assert_eq!(ctx.frames[1], "inner op");
    }

    #[test]
    fn test_error_context_format() {
        let ctx = ErrorContext::new().push("step").push("integrate");
        let fmt = ctx.format();
        assert!(fmt.contains("step"));
        assert!(fmt.contains("integrate"));
    }

    #[test]
    fn test_error_context_is_empty() {
        let ctx = ErrorContext::new();
        assert!(ctx.is_empty());
        let ctx = ctx.push("something");
        assert!(!ctx.is_empty());
    }

    #[test]
    fn test_error_context_wrap() {
        let ctx = ErrorContext::new().push("resolving contacts");
        let ew = ctx.wrap(Error::SolverConvergenceFailed {
            iterations: 30,
            tolerance: 1e-4,
        });
        assert!(ew.full_message().contains("30"));
        assert!(ew.full_message().contains("resolving contacts"));
    }

    // ---- ErrorWithContext tests --------------------------------------------

    #[test]
    fn test_error_with_context_no_context() {
        let ew = ErrorWithContext::new(Error::general("test"));
        assert_eq!(ew.full_message(), "test");
    }

    #[test]
    fn test_error_with_context_with_frame() {
        let ew = ErrorWithContext::new(Error::BodySleeping(1))
            .with_context("applying impulse")
            .with_context("user request");
        let msg = ew.full_message();
        assert!(msg.contains("sleeping"));
        assert!(msg.contains("applying impulse"));
        assert!(msg.contains("user request"));
    }

    #[test]
    fn test_error_with_context_to_json() {
        let ew = ErrorWithContext::new(Error::InvalidHandle(99)).with_context("test context");
        let json = ew.to_json();
        assert!(json.contains("InvalidHandle") || json.contains("99"));
    }

    #[test]
    fn test_error_with_context_display() {
        let ew = ErrorWithContext::new(Error::General("display test".into()));
        assert!(ew.to_string().contains("display test"));
    }

    // ---- ErrorLog tests ---------------------------------------------------

    #[test]
    fn test_error_log_basic() {
        let mut log = ErrorLog::new();
        assert!(log.is_empty());
        log.record(0, LogLevel::Warning, Error::BodySleeping(1), None);
        assert_eq!(log.len(), 1);
        assert!(!log.has_errors());
    }

    #[test]
    fn test_error_log_has_errors() {
        let mut log = ErrorLog::new();
        log.record(
            5,
            LogLevel::Error,
            Error::SimulationDiverged { step: 5 },
            None,
        );
        assert!(log.has_errors());
        assert!(!log.has_fatal());
    }

    #[test]
    fn test_error_log_has_fatal() {
        let mut log = ErrorLog::new();
        log.record(10, LogLevel::Fatal, Error::general("fatal"), None);
        assert!(log.has_fatal());
        assert!(log.has_errors());
    }

    #[test]
    fn test_error_log_filter_by_level() {
        let mut log = ErrorLog::new();
        log.record(0, LogLevel::Info, Error::general("info"), None);
        log.record(1, LogLevel::Warning, Error::BodySleeping(2), None);
        log.record(2, LogLevel::Error, Error::general("err"), None);
        let errors = log.entries_at_level(LogLevel::Error);
        assert_eq!(errors.len(), 1);
        let warnings_up = log.entries_at_level(LogLevel::Warning);
        assert_eq!(warnings_up.len(), 2);
    }

    #[test]
    fn test_error_log_filter_by_step() {
        let mut log = ErrorLog::new();
        log.record(3, LogLevel::Warning, Error::BodySleeping(1), None);
        log.record(3, LogLevel::Error, Error::general("two"), None);
        log.record(7, LogLevel::Info, Error::general("other"), None);
        assert_eq!(log.entries_at_step(3).len(), 2);
        assert_eq!(log.entries_at_step(7).len(), 1);
        assert_eq!(log.entries_at_step(0).len(), 0);
    }

    #[test]
    fn test_error_log_capacity() {
        let mut log = ErrorLog::with_capacity(3);
        for i in 0u64..5 {
            log.record(
                i,
                LogLevel::Warning,
                Error::general(format!("e{}", i)),
                None,
            );
        }
        assert_eq!(log.len(), 3);
    }

    #[test]
    fn test_error_log_clear() {
        let mut log = ErrorLog::new();
        log.record(0, LogLevel::Info, Error::general("x"), None);
        log.clear();
        assert!(log.is_empty());
    }

    #[test]
    fn test_error_log_to_json() {
        let mut log = ErrorLog::new();
        log.record(
            5,
            LogLevel::Error,
            Error::SimulationDiverged { step: 5 },
            None,
        );
        let json = log.to_json();
        assert!(json.contains("SimulationDiverged"), "json={}", json);
    }

    #[test]
    fn test_error_log_format() {
        let mut log = ErrorLog::new();
        log.record(1, LogLevel::Warning, Error::BodySleeping(4), None);
        let fmt = log.format();
        assert!(fmt.contains("WARN"));
        assert!(fmt.contains("step=1"));
    }

    #[test]
    fn test_error_log_with_context() {
        let mut log = ErrorLog::new();
        let ctx = ErrorContext::new().push("applying gravity");
        log.record(0, LogLevel::Error, Error::general("test"), Some(ctx));
        let fmt = log.format();
        assert!(fmt.contains("applying gravity"));
    }

    // ---- LogLevel ordering ------------------------------------------------

    #[test]
    fn test_log_level_ordering() {
        assert!(LogLevel::Info < LogLevel::Warning);
        assert!(LogLevel::Warning < LogLevel::Error);
        assert!(LogLevel::Error < LogLevel::Fatal);
    }

    #[test]
    fn test_log_level_display() {
        assert_eq!(LogLevel::Info.to_string(), "INFO");
        assert_eq!(LogLevel::Warning.to_string(), "WARN");
        assert_eq!(LogLevel::Error.to_string(), "ERROR");
        assert_eq!(LogLevel::Fatal.to_string(), "FATAL");
    }
}
