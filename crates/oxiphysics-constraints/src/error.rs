// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Error types for oxiphysics-constraints.
//!
//! Provides a comprehensive error hierarchy covering solver convergence failures,
//! constraint configuration errors, numerical issues, and body handle problems.

use thiserror::Error;

// ── Primary Error Type ────────────────────────────────────────────────────────

/// Main error type for the constraints module.
#[derive(Debug, Error)]
pub enum Error {
    /// Generic catch-all error with a human-readable message.
    #[error("{0}")]
    General(String),

    /// Solver failed to converge within the allotted iterations.
    #[error(
        "solver did not converge: residual {residual:.3e} after {iterations} iterations (tolerance {tolerance:.3e})"
    )]
    ConvergenceFailure {
        /// The residual magnitude at termination.
        residual: f64,
        /// Number of iterations performed.
        iterations: usize,
        /// The convergence tolerance that was not met.
        tolerance: f64,
    },

    /// A body handle referenced by a constraint is invalid (body not in the set).
    #[error("invalid body handle: id={id}")]
    InvalidBodyHandle {
        /// The invalid handle ID.
        id: u64,
    },

    /// A constraint is numerically degenerate (e.g. zero-length axis).
    #[error("degenerate constraint: {reason}")]
    DegenerateConstraint {
        /// Human-readable reason for the degeneracy.
        reason: String,
    },

    /// The constraint matrix is singular or near-singular; effective mass cannot be computed.
    #[error("singular constraint matrix: condition number estimate {condition:.3e}")]
    SingularMatrix {
        /// Estimated condition number (or f64::INFINITY if completely singular).
        condition: f64,
    },

    /// A configuration parameter is out of its valid range.
    #[error(
        "configuration error: parameter '{parameter}' value {value} is outside valid range [{min}, {max}]"
    )]
    ConfigurationError {
        /// Name of the offending parameter.
        parameter: &'static str,
        /// The value that was supplied.
        value: f64,
        /// Minimum allowed value.
        min: f64,
        /// Maximum allowed value.
        max: f64,
    },

    /// Time step is zero or negative, which would cause division-by-zero in the solver.
    #[error("invalid time step: dt={dt} (must be > 0)")]
    InvalidTimestep {
        /// The invalid time-step value.
        dt: f64,
    },

    /// The solver ran out of memory (e.g. island buffer overflow).
    #[error("solver out of capacity: tried to add {requested} items but capacity is {capacity}")]
    CapacityExceeded {
        /// Number of items that were requested.
        requested: usize,
        /// Maximum allowed capacity.
        capacity: usize,
    },

    /// A joint limit was mis-configured (lower > upper).
    #[error("joint limit error: lower={lower} > upper={upper}")]
    InvalidJointLimits {
        /// Lower limit value.
        lower: f64,
        /// Upper limit value.
        upper: f64,
    },

    /// A PID / motor controller error (gain out of range, integrator overflow, etc.).
    #[error("controller error: {message}")]
    ControllerError {
        /// Description of what went wrong.
        message: String,
    },

    /// Error originating from the underlying rigid-body set (body lookup, etc.).
    #[error("rigid body error: {message}")]
    RigidBodyError {
        /// Upstream error message.
        message: String,
    },

    /// CCD (Continuous Collision Detection) failed to find a valid TOI.
    #[error("CCD failed: no time-of-impact found within [{tmin}, {tmax}]")]
    CcdFailure {
        /// Start of the search interval.
        tmin: f64,
        /// End of the search interval.
        tmax: f64,
    },

    /// PBD (Position-Based Dynamics) constraint has unsatisfiable compliance / stiffness.
    #[error("PBD compliance error: compliance={compliance} is non-positive")]
    PbdComplianceError {
        /// The invalid compliance value.
        compliance: f64,
    },

    /// Friction model received an invalid normal force (e.g. negative normal).
    #[error("friction model error: normal_force={normal_force} (must be >= 0)")]
    FrictionModelError {
        /// The invalid normal force value.
        normal_force: f64,
    },
}

// ── Result Alias ──────────────────────────────────────────────────────────────

/// Result type alias for the constraints module.
pub type Result<T> = std::result::Result<T, Error>;

// ── Convenience Constructors ──────────────────────────────────────────────────

impl Error {
    /// Create a [`Error::General`] from any `Display`-able value.
    #[allow(dead_code)]
    pub fn general(msg: impl std::fmt::Display) -> Self {
        Error::General(msg.to_string())
    }

    /// Create a [`Error::ConvergenceFailure`].
    #[allow(dead_code)]
    pub fn convergence(residual: f64, iterations: usize, tolerance: f64) -> Self {
        Error::ConvergenceFailure {
            residual,
            iterations,
            tolerance,
        }
    }

    /// Create a [`Error::InvalidBodyHandle`].
    #[allow(dead_code)]
    pub fn invalid_body(id: u64) -> Self {
        Error::InvalidBodyHandle { id }
    }

    /// Create a [`Error::DegenerateConstraint`].
    #[allow(dead_code)]
    pub fn degenerate(reason: impl Into<String>) -> Self {
        Error::DegenerateConstraint {
            reason: reason.into(),
        }
    }

    /// Create a [`Error::SingularMatrix`].
    #[allow(dead_code)]
    pub fn singular(condition: f64) -> Self {
        Error::SingularMatrix { condition }
    }

    /// Create a [`Error::ConfigurationError`].
    #[allow(dead_code)]
    pub fn config(parameter: &'static str, value: f64, min: f64, max: f64) -> Self {
        Error::ConfigurationError {
            parameter,
            value,
            min,
            max,
        }
    }

    /// Create a [`Error::InvalidTimestep`].
    #[allow(dead_code)]
    pub fn invalid_dt(dt: f64) -> Self {
        Error::InvalidTimestep { dt }
    }

    /// Create a [`Error::CapacityExceeded`].
    #[allow(dead_code)]
    pub fn capacity(requested: usize, capacity: usize) -> Self {
        Error::CapacityExceeded {
            requested,
            capacity,
        }
    }

    /// Create a [`Error::InvalidJointLimits`].
    #[allow(dead_code)]
    pub fn joint_limits(lower: f64, upper: f64) -> Self {
        Error::InvalidJointLimits { lower, upper }
    }

    /// Create a [`Error::ControllerError`].
    #[allow(dead_code)]
    pub fn controller(message: impl Into<String>) -> Self {
        Error::ControllerError {
            message: message.into(),
        }
    }

    /// Create a [`Error::CcdFailure`].
    #[allow(dead_code)]
    pub fn ccd_failure(tmin: f64, tmax: f64) -> Self {
        Error::CcdFailure { tmin, tmax }
    }

    /// Create a [`Error::PbdComplianceError`].
    #[allow(dead_code)]
    pub fn pbd_compliance(compliance: f64) -> Self {
        Error::PbdComplianceError { compliance }
    }

    /// Create a [`Error::FrictionModelError`].
    #[allow(dead_code)]
    pub fn friction_model(normal_force: f64) -> Self {
        Error::FrictionModelError { normal_force }
    }
}

// ── Error Classification Helpers ──────────────────────────────────────────────

impl Error {
    /// Returns `true` if this error represents a numerical / floating-point issue.
    #[allow(dead_code)]
    pub fn is_numerical(&self) -> bool {
        matches!(
            self,
            Error::SingularMatrix { .. }
                | Error::ConvergenceFailure { .. }
                | Error::DegenerateConstraint { .. }
        )
    }

    /// Returns `true` if this error indicates a configuration / API misuse.
    #[allow(dead_code)]
    pub fn is_configuration(&self) -> bool {
        matches!(
            self,
            Error::ConfigurationError { .. }
                | Error::InvalidTimestep { .. }
                | Error::InvalidJointLimits { .. }
                | Error::PbdComplianceError { .. }
                | Error::FrictionModelError { .. }
        )
    }

    /// Returns `true` if the error is potentially recoverable (retry with
    /// smaller dt or more iterations may succeed).
    #[allow(dead_code)]
    pub fn is_recoverable(&self) -> bool {
        matches!(
            self,
            Error::ConvergenceFailure { .. } | Error::CcdFailure { .. }
        )
    }
}

// ── Severity ──────────────────────────────────────────────────────────────────

/// Severity level for constraint errors — used for logging / filtering.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
#[allow(dead_code)]
pub enum Severity {
    /// Informational — the solver can continue without intervention.
    Info,
    /// Warning — results may be inaccurate.
    Warning,
    /// Error — the solver cannot produce a valid result.
    Error,
    /// Fatal — the simulation state is corrupt and must be reset.
    Fatal,
}

impl Error {
    /// Map each error variant to an advisory [`Severity`] level.
    #[allow(dead_code)]
    pub fn severity(&self) -> Severity {
        match self {
            Error::ConvergenceFailure { .. } => Severity::Warning,
            Error::CcdFailure { .. } => Severity::Warning,
            Error::DegenerateConstraint { .. } => Severity::Warning,
            Error::InvalidBodyHandle { .. } => Severity::Error,
            Error::SingularMatrix { .. } => Severity::Error,
            Error::ConfigurationError { .. } => Severity::Error,
            Error::InvalidTimestep { .. } => Severity::Error,
            Error::InvalidJointLimits { .. } => Severity::Error,
            Error::PbdComplianceError { .. } => Severity::Error,
            Error::FrictionModelError { .. } => Severity::Error,
            Error::CapacityExceeded { .. } => Severity::Fatal,
            Error::ControllerError { .. } => Severity::Warning,
            Error::RigidBodyError { .. } => Severity::Error,
            Error::General(_) => Severity::Error,
        }
    }
}

// ── Validation helpers ────────────────────────────────────────────────────────

/// Validate that a time-step is strictly positive.
#[allow(dead_code)]
pub fn validate_dt(dt: f64) -> Result<()> {
    if dt <= 0.0 || !dt.is_finite() {
        Err(Error::invalid_dt(dt))
    } else {
        Ok(())
    }
}

/// Validate that a joint limit pair is well-ordered (lower ≤ upper).
#[allow(dead_code)]
pub fn validate_joint_limits(lower: f64, upper: f64) -> Result<()> {
    if lower > upper {
        Err(Error::joint_limits(lower, upper))
    } else {
        Ok(())
    }
}

/// Validate that a configuration parameter is within `[min, max]`.
#[allow(dead_code)]
pub fn validate_param(parameter: &'static str, value: f64, min: f64, max: f64) -> Result<()> {
    if value < min || value > max || !value.is_finite() {
        Err(Error::config(parameter, value, min, max))
    } else {
        Ok(())
    }
}

/// Validate that a PBD compliance value is positive.
#[allow(dead_code)]
pub fn validate_compliance(compliance: f64) -> Result<()> {
    if compliance <= 0.0 || !compliance.is_finite() {
        Err(Error::pbd_compliance(compliance))
    } else {
        Ok(())
    }
}

/// Validate that a normal force is non-negative.
#[allow(dead_code)]
pub fn validate_normal_force(normal_force: f64) -> Result<()> {
    if normal_force < 0.0 || !normal_force.is_finite() {
        Err(Error::friction_model(normal_force))
    } else {
        Ok(())
    }
}

/// Validate that a value is strictly positive.
#[allow(dead_code)]
pub fn validate_positive(parameter: &'static str, value: f64) -> Result<()> {
    if value <= 0.0 || !value.is_finite() {
        Err(Error::config(parameter, value, f64::EPSILON, f64::INFINITY))
    } else {
        Ok(())
    }
}

/// Validate that a value is finite (not NaN or infinite).
#[allow(dead_code)]
pub fn validate_finite(parameter: &'static str, value: f64) -> Result<()> {
    if !value.is_finite() {
        Err(Error::config(
            parameter,
            value,
            f64::NEG_INFINITY,
            f64::INFINITY,
        ))
    } else {
        Ok(())
    }
}

/// Validate that a 3D vector is a unit vector (length ≈ 1).
///
/// Tolerance is ±1e-6.
#[allow(dead_code)]
pub fn validate_unit_vector(parameter: &'static str, v: [f64; 3]) -> Result<()> {
    let len_sq = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
    let deviation = (len_sq - 1.0).abs();
    if deviation > 1e-6 {
        Err(Error::degenerate(format!(
            "{parameter}: not a unit vector (|v|²={len_sq:.6}, deviation={deviation:.3e})"
        )))
    } else {
        Ok(())
    }
}

/// Validate that a 3D vector is non-zero (length > tolerance).
#[allow(dead_code)]
pub fn validate_nonzero_vector(parameter: &'static str, v: [f64; 3]) -> Result<()> {
    let len_sq = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
    if len_sq < 1e-15 {
        Err(Error::degenerate(format!(
            "{parameter}: zero-length vector (|v|²={len_sq:.3e})"
        )))
    } else {
        Ok(())
    }
}

/// Validate that `count` items can be inserted without exceeding `capacity`.
#[allow(dead_code)]
pub fn validate_capacity(current: usize, count: usize, capacity: usize) -> Result<()> {
    let total = current.saturating_add(count);
    if total > capacity {
        Err(Error::capacity(total, capacity))
    } else {
        Ok(())
    }
}

// ── Error Aggregator ──────────────────────────────────────────────────────────

/// Collects multiple non-fatal errors (warnings) during a solve step so they
/// can be inspected after the fact rather than causing an early return.
#[derive(Debug, Default)]
#[allow(dead_code)]
pub struct ErrorAccumulator {
    errors: Vec<Error>,
    fatal: bool,
}

impl ErrorAccumulator {
    /// Create a new, empty accumulator.
    #[allow(dead_code)]
    pub fn new() -> Self {
        Self::default()
    }

    /// Push an error onto the accumulator.
    ///
    /// If the error's severity is [`Severity::Fatal`], the `fatal` flag is set.
    #[allow(dead_code)]
    pub fn push(&mut self, e: Error) {
        if e.severity() == Severity::Fatal {
            self.fatal = true;
        }
        self.errors.push(e);
    }

    /// Returns `true` if a fatal error has been accumulated.
    #[allow(dead_code)]
    pub fn has_fatal(&self) -> bool {
        self.fatal
    }

    /// Returns `true` if any errors (including warnings) have been accumulated.
    #[allow(dead_code)]
    pub fn has_errors(&self) -> bool {
        !self.errors.is_empty()
    }

    /// Number of accumulated errors.
    #[allow(dead_code)]
    pub fn len(&self) -> usize {
        self.errors.len()
    }

    /// Returns `true` if no errors have been accumulated.
    #[allow(dead_code)]
    pub fn is_empty(&self) -> bool {
        self.errors.is_empty()
    }

    /// Drain all accumulated errors, returning them as a `Vec`.
    #[allow(dead_code)]
    pub fn drain(&mut self) -> Vec<Error> {
        self.fatal = false;
        std::mem::take(&mut self.errors)
    }

    /// Return a reference to all accumulated errors without consuming them.
    #[allow(dead_code)]
    pub fn errors(&self) -> &[Error] {
        &self.errors
    }

    /// If any fatal errors were accumulated, return `Err` with the first one;
    /// otherwise return `Ok(())`.
    #[allow(dead_code)]
    pub fn into_result(mut self) -> Result<()> {
        if self.fatal {
            Err(self
                .errors
                .drain(..)
                .find(|e| e.severity() == Severity::Fatal)
                .expect("value should be present"))
        } else {
            Ok(())
        }
    }
}

// ── Error Context ─────────────────────────────────────────────────────────────

/// Wraps an error with additional context about where it occurred.
#[derive(Debug)]
#[allow(dead_code)]
pub struct ErrorContext {
    /// The underlying error.
    pub error: Error,
    /// Human-readable description of the context (e.g. "solving island 3").
    pub context: String,
    /// Optional constraint index that produced the error.
    pub constraint_index: Option<usize>,
}

impl ErrorContext {
    /// Wrap an error with a context string.
    #[allow(dead_code)]
    pub fn new(error: Error, context: impl Into<String>) -> Self {
        ErrorContext {
            error,
            context: context.into(),
            constraint_index: None,
        }
    }

    /// Attach a constraint index to the context.
    #[allow(dead_code)]
    pub fn with_constraint(mut self, index: usize) -> Self {
        self.constraint_index = Some(index);
        self
    }
}

impl std::fmt::Display for ErrorContext {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if let Some(idx) = self.constraint_index {
            write!(f, "[constraint {}] {}: {}", idx, self.context, self.error)
        } else {
            write!(f, "{}: {}", self.context, self.error)
        }
    }
}

// ── Numeric Utilities ─────────────────────────────────────────────────────────

/// Compute a simple 2×2 condition number estimate via max/min absolute diagonal.
///
/// Returns `f64::INFINITY` if the minimum diagonal is zero.
#[allow(dead_code)]
pub fn condition_number_2x2(a00: f64, a11: f64) -> f64 {
    let max_d = a00.abs().max(a11.abs());
    let min_d = a00.abs().min(a11.abs());
    if min_d < 1e-15 {
        f64::INFINITY
    } else {
        max_d / min_d
    }
}

/// Return `Ok(())` if a 3×3 symmetric matrix (given as 6 unique entries in
/// row-major upper-triangular order) is positive definite via Sylvester's
/// criterion.
///
/// The entries are `[a00, a01, a02, a11, a12, a22]`.
#[allow(dead_code)]
pub fn check_positive_definite_3x3(m: [f64; 6]) -> Result<()> {
    let (a00, a11, a01, a02, a12, a22) = (m[0], m[3], m[1], m[2], m[4], m[5]);
    // Leading minors
    if a00 <= 0.0 {
        return Err(Error::singular(f64::INFINITY));
    }
    let minor2 = a00 * a11 - a01 * a01;
    if minor2 <= 0.0 {
        let cond = a00.abs().max(a11.abs()) / minor2.abs().max(1e-15);
        return Err(Error::singular(cond));
    }
    // det(3×3)
    let det = a00 * (a11 * a22 - a12 * a12) - a01 * (a01 * a22 - a12 * a02)
        + a02 * (a01 * a12 - a11 * a02);
    if det <= 0.0 {
        return Err(Error::singular(
            a00.max(a11).max(a22) / det.abs().max(1e-15),
        ));
    }
    Ok(())
}

// ── ErrorBatch ────────────────────────────────────────────────────────────────

/// A collection of errors accumulated during a solver pass.
///
/// Allows the solver to continue after a non-fatal error and report all problems
/// at once rather than stopping at the first failure.
#[derive(Debug, Default)]
#[allow(dead_code)]
pub struct ErrorBatch {
    errors: Vec<Error>,
}

#[allow(dead_code)]
impl ErrorBatch {
    /// Create a new, empty error batch.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add an error to the batch.
    pub fn push(&mut self, e: Error) {
        self.errors.push(e);
    }

    /// Number of errors recorded.
    pub fn len(&self) -> usize {
        self.errors.len()
    }

    /// Returns `true` if no errors have been recorded.
    pub fn is_empty(&self) -> bool {
        self.errors.is_empty()
    }

    /// Returns the first error if any, as a `Result`.
    ///
    /// If the batch is empty, returns `Ok(())`.
    pub fn to_result(&self) -> Result<()> {
        if let Some(e) = self.errors.first() {
            Err(Error::General(e.to_string()))
        } else {
            Ok(())
        }
    }

    /// Iterate over the recorded errors.
    pub fn iter(&self) -> std::slice::Iter<'_, Error> {
        self.errors.iter()
    }

    /// Clear all recorded errors.
    pub fn clear(&mut self) {
        self.errors.clear();
    }

    /// Drain all errors into a `Vec`Error`.
    pub fn drain(&mut self) -> Vec<Error> {
        std::mem::take(&mut self.errors)
    }
}

// ── Additional validation helpers ─────────────────────────────────────────────

/// Validate that a time step `dt` is strictly positive.
///
/// Returns `Err(Error::InvalidTimestep)` if `dt <= 0` or `dt` is NaN.
#[allow(dead_code)]
pub fn validate_timestep(dt: f64) -> Result<()> {
    if dt.is_nan() || dt <= 0.0 {
        Err(Error::InvalidTimestep { dt })
    } else {
        Ok(())
    }
}

/// Validate that joint limits are properly ordered: `lower <= upper`.
/// Validate that a friction coefficient is non-negative.
///
/// Returns `Err(Error::ConfigurationError)` if `mu < 0` or is NaN.
#[allow(dead_code)]
pub fn validate_friction_coefficient(mu: f64) -> Result<()> {
    if mu.is_nan() || mu < 0.0 {
        Err(Error::config("friction_coefficient", mu, 0.0, f64::MAX))
    } else {
        Ok(())
    }
}

/// Validate that a restitution coefficient is in [0, 1].
#[allow(dead_code)]
pub fn validate_restitution(e: f64) -> Result<()> {
    if e.is_nan() || !(0.0..=1.0).contains(&e) {
        Err(Error::config("restitution", e, 0.0, 1.0))
    } else {
        Ok(())
    }
}

/// Validate that an angular velocity limit is non-negative.
#[allow(dead_code)]
pub fn validate_angular_velocity_limit(omega_max: f64) -> Result<()> {
    if omega_max.is_nan() || omega_max < 0.0 {
        Err(Error::config(
            "angular_velocity_limit",
            omega_max,
            0.0,
            f64::MAX,
        ))
    } else {
        Ok(())
    }
}

/// Validate that a Baumgarte stabilisation factor is in (0, 1].
#[allow(dead_code)]
pub fn validate_baumgarte(beta: f64) -> Result<()> {
    if beta.is_nan() || beta <= 0.0 || beta > 1.0 {
        Err(Error::config("baumgarte", beta, 1e-15, 1.0))
    } else {
        Ok(())
    }
}

// ── ErrorSummary ──────────────────────────────────────────────────────────────

/// A summary of error counts by category, derived from an `ErrorBatch`.
#[derive(Debug, Clone, Default)]
#[allow(dead_code)]
pub struct ErrorSummary {
    /// Total number of errors.
    pub total: usize,
    /// Number of convergence failures.
    pub convergence_failures: usize,
    /// Number of singular-matrix errors.
    pub singular_matrices: usize,
    /// Number of invalid body-handle errors.
    pub invalid_handles: usize,
    /// Number of configuration errors.
    pub configuration_errors: usize,
    /// Number of other errors.
    pub other: usize,
}

#[allow(dead_code)]
impl ErrorSummary {
    /// Build a summary from an `ErrorBatch`.
    pub fn from_batch(batch: &ErrorBatch) -> Self {
        let mut s = ErrorSummary {
            total: batch.len(),
            ..Default::default()
        };
        for e in batch.iter() {
            match e {
                Error::ConvergenceFailure { .. } => s.convergence_failures += 1,
                Error::SingularMatrix { .. } => s.singular_matrices += 1,
                Error::InvalidBodyHandle { .. } => s.invalid_handles += 1,
                Error::ConfigurationError { .. } => s.configuration_errors += 1,
                _ => s.other += 1,
            }
        }
        s
    }

    /// Returns `true` if any errors were recorded.
    pub fn has_errors(&self) -> bool {
        self.total > 0
    }

    /// Returns `true` if any convergence failures were recorded.
    pub fn has_convergence_failures(&self) -> bool {
        self.convergence_failures > 0
    }
}

// ── SolverDiagnostics ─────────────────────────────────────────────────────────

/// Diagnostic information captured during a single solver step.
#[derive(Debug, Clone, Default)]
#[allow(dead_code)]
pub struct SolverDiagnostics {
    /// Number of PGS iterations performed.
    pub iterations: usize,
    /// Maximum residual at the end of solving.
    pub max_residual: f64,
    /// Whether the solver converged.
    pub converged: bool,
    /// Number of active constraints (non-zero impulse).
    pub active_constraints: usize,
    /// Number of sleeping islands (skipped this step).
    pub sleeping_islands: usize,
    /// Wall-clock time in microseconds (advisory, not measured here).
    pub time_us: u64,
}

#[allow(dead_code)]
impl SolverDiagnostics {
    /// Create empty diagnostics.
    pub fn new() -> Self {
        Self::default()
    }

    /// Mark the solver as converged.
    pub fn mark_converged(&mut self, iterations: usize, max_residual: f64) {
        self.converged = true;
        self.iterations = iterations;
        self.max_residual = max_residual;
    }

    /// Check if the diagnostics represent a clean solve (converged, no large residual).
    pub fn is_clean(&self, tolerance: f64) -> bool {
        self.converged && self.max_residual <= tolerance
    }

    /// Convert to a `Result<()>`, returning a `ConvergenceFailure` if not converged.
    pub fn to_result(&self, tolerance: f64) -> Result<()> {
        if self.is_clean(tolerance) {
            Ok(())
        } else {
            Err(Error::convergence(
                self.max_residual,
                self.iterations,
                tolerance,
            ))
        }
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_general_error_display() {
        let e = Error::general("something went wrong");
        let msg = e.to_string();
        assert!(msg.contains("something went wrong"), "msg={msg}");
    }

    #[test]
    fn test_convergence_failure_display() {
        let e = Error::convergence(1e-3, 50, 1e-6);
        let msg = e.to_string();
        assert!(msg.contains("50"), "msg={msg}");
        assert!(msg.contains("converge"), "msg={msg}");
    }

    #[test]
    fn test_invalid_body_handle_display() {
        let e = Error::invalid_body(42);
        let msg = e.to_string();
        assert!(msg.contains("42"), "msg={msg}");
    }

    #[test]
    fn test_degenerate_constraint_display() {
        let e = Error::degenerate("zero-length hinge axis");
        let msg = e.to_string();
        assert!(msg.contains("zero-length"), "msg={msg}");
    }

    #[test]
    fn test_singular_matrix_display() {
        let e = Error::singular(f64::INFINITY);
        let msg = e.to_string();
        assert!(msg.contains("singular"), "msg={msg}");
    }

    #[test]
    fn test_config_error_display() {
        let e = Error::config("beta", 2.0, 0.0, 1.0);
        let msg = e.to_string();
        assert!(msg.contains("beta"), "msg={msg}");
        assert!(msg.contains("2"), "msg={msg}");
    }

    #[test]
    fn test_invalid_timestep_display() {
        let e = Error::invalid_dt(-0.01);
        let msg = e.to_string();
        assert!(msg.contains("dt"), "msg={msg}");
    }

    #[test]
    fn test_capacity_exceeded_display() {
        let e = Error::capacity(100, 64);
        let msg = e.to_string();
        assert!(msg.contains("100"), "msg={msg}");
        assert!(msg.contains("64"), "msg={msg}");
    }

    #[test]
    fn test_joint_limits_display() {
        let e = Error::joint_limits(1.0, -1.0);
        let msg = e.to_string();
        assert!(msg.contains("lower"), "msg={msg}");
    }

    #[test]
    fn test_controller_error_display() {
        let e = Error::controller("integral windup");
        let msg = e.to_string();
        assert!(msg.contains("integral windup"), "msg={msg}");
    }

    #[test]
    fn test_ccd_failure_display() {
        let e = Error::ccd_failure(0.0, 1.0);
        let msg = e.to_string();
        assert!(msg.contains("CCD"), "msg={msg}");
    }

    #[test]
    fn test_pbd_compliance_display() {
        let e = Error::pbd_compliance(-1e-4);
        let msg = e.to_string();
        assert!(msg.contains("compliance"), "msg={msg}");
    }

    #[test]
    fn test_friction_model_display() {
        let e = Error::friction_model(-5.0);
        let msg = e.to_string();
        assert!(msg.contains("normal_force"), "msg={msg}");
    }

    #[test]
    fn test_is_numerical() {
        assert!(Error::convergence(1e-3, 10, 1e-6).is_numerical());
        assert!(Error::singular(1e12).is_numerical());
        assert!(Error::degenerate("axis").is_numerical());
        assert!(!Error::invalid_dt(0.0).is_numerical());
    }

    #[test]
    fn test_is_configuration() {
        assert!(Error::invalid_dt(-1.0).is_configuration());
        assert!(Error::joint_limits(1.0, -1.0).is_configuration());
        assert!(Error::pbd_compliance(-0.1).is_configuration());
        assert!(Error::friction_model(-1.0).is_configuration());
        assert!(!Error::convergence(1.0, 1, 0.0).is_configuration());
    }

    #[test]
    fn test_is_recoverable() {
        assert!(Error::convergence(1.0, 10, 1e-6).is_recoverable());
        assert!(Error::ccd_failure(0.0, 1.0).is_recoverable());
        assert!(!Error::singular(f64::INFINITY).is_recoverable());
    }

    #[test]
    fn test_severity_levels() {
        assert_eq!(
            Error::convergence(1.0, 1, 0.001).severity(),
            Severity::Warning
        );
        assert_eq!(Error::singular(1e20).severity(), Severity::Error);
        assert_eq!(Error::capacity(100, 10).severity(), Severity::Fatal);
        assert_eq!(Error::ccd_failure(0.0, 1.0).severity(), Severity::Warning);
    }

    #[test]
    fn test_severity_ordering() {
        assert!(Severity::Info < Severity::Warning);
        assert!(Severity::Warning < Severity::Error);
        assert!(Severity::Error < Severity::Fatal);
    }

    #[test]
    fn test_validate_dt_ok() {
        assert!(validate_dt(0.01).is_ok());
        assert!(validate_dt(1e-6).is_ok());
    }

    #[test]
    fn test_validate_dt_err() {
        assert!(validate_dt(0.0).is_err());
        assert!(validate_dt(-0.01).is_err());
        assert!(validate_dt(f64::NAN).is_err());
        assert!(validate_dt(f64::INFINITY).is_err());
    }

    #[test]
    fn test_validate_joint_limits_ok() {
        assert!(validate_joint_limits(-1.0, 1.0).is_ok());
        assert!(validate_joint_limits(0.0, 0.0).is_ok());
    }

    #[test]
    fn test_validate_joint_limits_err() {
        assert!(validate_joint_limits(1.0, -1.0).is_err());
    }

    #[test]
    fn test_validate_param_ok() {
        assert!(validate_param("beta", 0.2, 0.0, 1.0).is_ok());
        assert!(validate_param("beta", 0.0, 0.0, 1.0).is_ok());
        assert!(validate_param("beta", 1.0, 0.0, 1.0).is_ok());
    }

    #[test]
    fn test_validate_param_err() {
        assert!(validate_param("beta", -0.1, 0.0, 1.0).is_err());
        assert!(validate_param("beta", 1.1, 0.0, 1.0).is_err());
        assert!(validate_param("beta", f64::NAN, 0.0, 1.0).is_err());
    }

    #[test]
    fn test_validate_compliance_ok() {
        assert!(validate_compliance(1e-4).is_ok());
        assert!(validate_compliance(1.0).is_ok());
    }

    #[test]
    fn test_validate_compliance_err() {
        assert!(validate_compliance(0.0).is_err());
        assert!(validate_compliance(-1e-4).is_err());
        assert!(validate_compliance(f64::NAN).is_err());
    }

    #[test]
    fn test_validate_normal_force_ok() {
        assert!(validate_normal_force(0.0).is_ok());
        assert!(validate_normal_force(100.0).is_ok());
    }

    #[test]
    fn test_validate_normal_force_err() {
        assert!(validate_normal_force(-1.0).is_err());
        assert!(validate_normal_force(f64::NAN).is_err());
    }

    // ── validate_positive tests ───────────────────────────────────────────────

    #[test]
    fn test_validate_positive_ok() {
        assert!(validate_positive("mass", 1.0).is_ok());
        assert!(validate_positive("mass", f64::EPSILON).is_ok());
        assert!(validate_positive("mass", 1e10).is_ok());
    }

    #[test]
    fn test_validate_positive_err_zero() {
        assert!(validate_positive("mass", 0.0).is_err());
    }

    #[test]
    fn test_validate_positive_err_negative() {
        assert!(validate_positive("mass", -1.0).is_err());
    }

    #[test]
    fn test_validate_positive_err_nan() {
        assert!(validate_positive("mass", f64::NAN).is_err());
    }

    #[test]
    fn test_validate_positive_err_infinity() {
        assert!(validate_positive("mass", f64::INFINITY).is_err());
    }

    // ── validate_finite tests ─────────────────────────────────────────────────

    #[test]
    fn test_validate_finite_ok() {
        assert!(validate_finite("x", 0.0).is_ok());
        assert!(validate_finite("x", -1e9).is_ok());
        assert!(validate_finite("x", 1e9).is_ok());
    }

    #[test]
    fn test_validate_finite_err_nan() {
        assert!(validate_finite("x", f64::NAN).is_err());
    }

    #[test]
    fn test_validate_finite_err_pos_inf() {
        assert!(validate_finite("x", f64::INFINITY).is_err());
    }

    #[test]
    fn test_validate_finite_err_neg_inf() {
        assert!(validate_finite("x", f64::NEG_INFINITY).is_err());
    }

    // ── validate_unit_vector tests ────────────────────────────────────────────

    #[test]
    fn test_validate_unit_vector_x_axis() {
        assert!(validate_unit_vector("axis", [1.0, 0.0, 0.0]).is_ok());
    }

    #[test]
    fn test_validate_unit_vector_y_axis() {
        assert!(validate_unit_vector("axis", [0.0, 1.0, 0.0]).is_ok());
    }

    #[test]
    fn test_validate_unit_vector_z_axis() {
        assert!(validate_unit_vector("axis", [0.0, 0.0, 1.0]).is_ok());
    }

    #[test]
    fn test_validate_unit_vector_diagonal() {
        let s = 1.0_f64 / 3.0_f64.sqrt();
        assert!(validate_unit_vector("axis", [s, s, s]).is_ok());
    }

    #[test]
    fn test_validate_unit_vector_non_unit() {
        assert!(validate_unit_vector("axis", [2.0, 0.0, 0.0]).is_err());
    }

    #[test]
    fn test_validate_unit_vector_zero() {
        assert!(validate_unit_vector("axis", [0.0, 0.0, 0.0]).is_err());
    }

    // ── validate_nonzero_vector tests ─────────────────────────────────────────

    #[test]
    fn test_validate_nonzero_vector_ok() {
        assert!(validate_nonzero_vector("v", [1.0, 0.0, 0.0]).is_ok());
        assert!(validate_nonzero_vector("v", [0.1, 0.1, 0.1]).is_ok());
    }

    #[test]
    fn test_validate_nonzero_vector_zero() {
        assert!(validate_nonzero_vector("v", [0.0, 0.0, 0.0]).is_err());
    }

    // ── validate_capacity tests ───────────────────────────────────────────────

    #[test]
    fn test_validate_capacity_ok() {
        assert!(validate_capacity(5, 3, 10).is_ok());
        assert!(validate_capacity(0, 10, 10).is_ok());
    }

    #[test]
    fn test_validate_capacity_err() {
        assert!(validate_capacity(8, 3, 10).is_err());
        assert!(validate_capacity(10, 1, 10).is_err());
    }

    // ── ErrorAccumulator tests ────────────────────────────────────────────────

    #[test]
    fn test_error_accumulator_new_empty() {
        let acc = ErrorAccumulator::new();
        assert!(acc.is_empty());
        assert!(!acc.has_errors());
        assert!(!acc.has_fatal());
    }

    #[test]
    fn test_error_accumulator_push_warning() {
        let mut acc = ErrorAccumulator::new();
        acc.push(Error::convergence(1e-3, 10, 1e-6));
        assert!(acc.has_errors());
        assert!(!acc.has_fatal());
        assert_eq!(acc.len(), 1);
    }

    #[test]
    fn test_error_accumulator_push_fatal() {
        let mut acc = ErrorAccumulator::new();
        acc.push(Error::capacity(100, 64));
        assert!(acc.has_fatal());
    }

    #[test]
    fn test_error_accumulator_drain() {
        let mut acc = ErrorAccumulator::new();
        acc.push(Error::general("a"));
        acc.push(Error::general("b"));
        let drained = acc.drain();
        assert_eq!(drained.len(), 2);
        assert!(acc.is_empty());
        assert!(!acc.has_fatal());
    }

    #[test]
    fn test_error_accumulator_errors_ref() {
        let mut acc = ErrorAccumulator::new();
        acc.push(Error::invalid_dt(-1.0));
        assert_eq!(acc.errors().len(), 1);
    }

    #[test]
    fn test_error_accumulator_into_result_ok_no_fatal() {
        let mut acc = ErrorAccumulator::new();
        acc.push(Error::convergence(0.1, 5, 1e-6));
        assert!(acc.into_result().is_ok());
    }

    #[test]
    fn test_error_accumulator_into_result_err_fatal() {
        let mut acc = ErrorAccumulator::new();
        acc.push(Error::capacity(100, 64));
        assert!(acc.into_result().is_err());
    }

    // ── ErrorContext tests ────────────────────────────────────────────────────

    #[test]
    fn test_error_context_display_no_index() {
        let ctx = ErrorContext::new(Error::general("oops"), "during island solve");
        let s = ctx.to_string();
        assert!(s.contains("during island solve"), "s={s}");
        assert!(s.contains("oops"), "s={s}");
    }

    #[test]
    fn test_error_context_display_with_index() {
        let ctx =
            ErrorContext::new(Error::singular(1e10), "during velocity solve").with_constraint(7);
        let s = ctx.to_string();
        assert!(s.contains("7"), "s={s}");
        assert!(s.contains("singular"), "s={s}");
    }

    // ── condition_number_2x2 tests ────────────────────────────────────────────

    #[test]
    fn test_condition_number_2x2_equal_diagonal() {
        let cond = condition_number_2x2(2.0, 2.0);
        assert!((cond - 1.0).abs() < 1e-15, "cond={cond}");
    }

    #[test]
    fn test_condition_number_2x2_ratio() {
        let cond = condition_number_2x2(10.0, 2.0);
        assert!((cond - 5.0).abs() < 1e-15, "cond={cond}");
    }

    #[test]
    fn test_condition_number_2x2_zero_min() {
        let cond = condition_number_2x2(5.0, 0.0);
        assert!(cond.is_infinite());
    }

    // ── check_positive_definite_3x3 tests ────────────────────────────────────

    #[test]
    fn test_positive_definite_3x3_identity() {
        // Identity matrix: [1,0,0,1,0,1]
        assert!(check_positive_definite_3x3([1.0, 0.0, 0.0, 1.0, 0.0, 1.0]).is_ok());
    }

    #[test]
    fn test_positive_definite_3x3_diagonal() {
        // Diagonal [2, 3, 4] -> a00=2, a01=0, a02=0, a11=3, a12=0, a22=4
        assert!(check_positive_definite_3x3([2.0, 0.0, 0.0, 3.0, 0.0, 4.0]).is_ok());
    }

    #[test]
    fn test_positive_definite_3x3_non_pd() {
        // Negative diagonal — not PD
        assert!(check_positive_definite_3x3([-1.0, 0.0, 0.0, 1.0, 0.0, 1.0]).is_err());
    }

    #[test]
    fn test_positive_definite_3x3_zero_diagonal() {
        assert!(check_positive_definite_3x3([0.0, 0.0, 0.0, 1.0, 0.0, 1.0]).is_err());
    }

    // ── ErrorBatch tests ─────────────────────────────────────────────────────

    #[test]
    fn test_error_batch_empty() {
        let batch: ErrorBatch = ErrorBatch::new();
        assert!(batch.is_empty());
        assert!(batch.to_result().is_ok());
    }

    #[test]
    fn test_error_batch_push_and_count() {
        let mut batch = ErrorBatch::new();
        batch.push(Error::general("first"));
        batch.push(Error::general("second"));
        assert_eq!(batch.len(), 2);
        assert!(!batch.is_empty());
    }

    #[test]
    fn test_error_batch_to_result_nonempty() {
        let mut batch = ErrorBatch::new();
        batch.push(Error::general("oops"));
        assert!(batch.to_result().is_err());
    }

    #[test]
    fn test_error_batch_clear() {
        let mut batch = ErrorBatch::new();
        batch.push(Error::general("x"));
        batch.clear();
        assert!(batch.is_empty());
    }

    // ── validate_timestep tests ───────────────────────────────────────────────

    #[test]
    fn test_validate_timestep_valid() {
        assert!(validate_timestep(1.0 / 60.0).is_ok());
    }

    #[test]
    fn test_validate_timestep_zero() {
        assert!(validate_timestep(0.0).is_err());
    }

    #[test]
    fn test_validate_timestep_negative() {
        assert!(validate_timestep(-0.01).is_err());
    }

    // ── validate_joint_limits tests ───────────────────────────────────────────

    #[test]
    fn test_validate_joint_limits_valid() {
        assert!(validate_joint_limits(-1.0, 1.0).is_ok());
        assert!(validate_joint_limits(0.0, 0.0).is_ok()); // degenerate but technically valid
    }

    #[test]
    fn test_validate_joint_limits_inverted() {
        assert!(validate_joint_limits(1.0, -1.0).is_err());
    }

    // ── ErrorSummary tests ────────────────────────────────────────────────────

    #[test]
    fn test_error_summary_counts() {
        let mut batch = ErrorBatch::new();
        batch.push(Error::convergence(0.01, 100, 1e-6));
        batch.push(Error::convergence(0.02, 200, 1e-6));
        batch.push(Error::general("other"));
        let summary = ErrorSummary::from_batch(&batch);
        assert_eq!(summary.convergence_failures, 2);
        assert_eq!(summary.total, 3);
    }
}
