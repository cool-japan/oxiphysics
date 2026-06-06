// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Error types for oxiphysics-rigid.
//!
//! Provides a structured error hierarchy covering body management, joint
//! constraints, solver failures, and I/O operations.  All variants are
//! non-exhaustive to allow future additions without breaking downstream
//! code.

use thiserror::Error;

// ---------------------------------------------------------------------------
// Top-level error type
// ---------------------------------------------------------------------------

/// Main error type for the rigid-body module.
///
/// Variants cover every failure mode that can occur during scene setup,
/// simulation stepping, and serialization/deserialization.
#[derive(Debug, Error, Clone, PartialEq)]
#[non_exhaustive]
pub enum Error {
    // ── Generic ────────────────────────────────────────────────────────────
    /// Catch-all for errors that do not fit any other category.
    #[error("rigid body error: {0}")]
    General(String),

    // ── Body management ────────────────────────────────────────────────────
    /// A `crate::BodyHandle` no longer points to a valid body.
    #[error("body handle {handle:?} is invalid (generation mismatch or already removed)")]
    InvalidBodyHandle {
        /// The stale or out-of-range handle.
        handle: u64,
    },

    /// Attempt to insert a body with a non-positive mass into a dynamic slot.
    #[error("body mass must be positive, got {mass}")]
    NonPositiveMass {
        /// The offending mass value.
        mass: f64,
    },

    /// Attempt to insert a body with a singular (non-invertible) inertia tensor.
    #[error("body inertia tensor is singular or near-singular (det={det:.3e})")]
    SingularInertiaTensor {
        /// Determinant of the supplied inertia tensor.
        det: f64,
    },

    // ── Collider / geometry ────────────────────────────────────────────────
    /// A collider shape has invalid parameters (e.g. negative radius).
    #[error("collider shape is invalid: {reason}")]
    InvalidColliderShape {
        /// Human-readable description of why the shape is invalid.
        reason: String,
    },

    /// A collider handle does not exist in the collider set.
    #[error("collider handle {handle:?} not found")]
    ColliderNotFound {
        /// The missing collider handle.
        handle: u64,
    },

    // ── Joint / constraint ─────────────────────────────────────────────────
    /// A joint references a body that does not exist.
    #[error("joint references unknown body {body_id}")]
    JointBodyNotFound {
        /// Identifier of the missing body.
        body_id: u64,
    },

    /// A joint limit is malformed (min > max).
    #[error("joint limit invalid: min={min} > max={max}")]
    InvalidJointLimit {
        /// Lower bound of the joint limit.
        min: f64,
        /// Upper bound of the joint limit.
        max: f64,
    },

    /// The constraint solver diverged during a step.
    #[error("constraint solver diverged at iteration {iteration} (residual={residual:.3e})")]
    SolverDivergence {
        /// Iteration number at which divergence was detected.
        iteration: u32,
        /// Residual norm at the point of divergence.
        residual: f64,
    },

    /// Maximum number of solver iterations was reached without convergence.
    #[error(
        "solver did not converge in {max_iterations} iterations (final residual={residual:.3e})"
    )]
    SolverNotConverged {
        /// Maximum iteration count the solver was allowed to run.
        max_iterations: u32,
        /// Residual at termination.
        residual: f64,
    },

    // ── CCD / TOI ──────────────────────────────────────────────────────────
    /// CCD time-of-impact computation failed.
    #[error("CCD TOI calculation failed: {reason}")]
    CcdToiFailed {
        /// Reason for the failure (e.g. "shape not supported").
        reason: String,
    },

    /// A requested time-of-impact is outside `[0, 1]`.
    #[error("TOI value {toi} is out of range [0, 1]")]
    ToiOutOfRange {
        /// The invalid TOI value.
        toi: f64,
    },

    // ── Island / sleeping ──────────────────────────────────────────────────
    /// An island ID is out of range.
    #[error("island id {id} is out of range (total islands: {total})")]
    IslandOutOfRange {
        /// Requested island identifier.
        id: u32,
        /// Total number of islands currently allocated.
        total: u32,
    },

    // ── Motor / actuator ───────────────────────────────────────────────────
    /// A motor parameter (e.g. resistance, torque constant) is physically
    /// invalid.
    #[error("motor parameter '{param}' has invalid value {value}")]
    InvalidMotorParameter {
        /// Name of the invalid parameter.
        param: &'static str,
        /// The supplied value.
        value: f64,
    },

    /// Motor temperature exceeded the winding thermal limit.
    #[error("motor thermal limit exceeded: winding temperature {temp:.1}°C > limit {limit:.1}°C")]
    MotorThermalLimit {
        /// Current winding temperature in °C.
        temp: f64,
        /// Thermal limit in °C.
        limit: f64,
    },

    // ── Fluid coupling ─────────────────────────────────────────────────────
    /// A fluid property (density, viscosity) is physically impossible.
    #[error("fluid property '{prop}' has invalid value {value}")]
    InvalidFluidProperty {
        /// Name of the invalid fluid property.
        prop: &'static str,
        /// The supplied value.
        value: f64,
    },

    // ── Articulated body / Featherstone ────────────────────────────────────
    /// An articulated body computation failed (e.g. loop closure, Featherstone
    /// forward dynamics).
    #[error("articulated body error: {reason}")]
    ArticulatedBodyError {
        /// Human-readable description of the failure.
        reason: String,
    },

    /// A spatial algebra operation (cross-product, transform, etc.) is
    /// invalid (e.g. zero-length axis).
    #[error("spatial algebra error: {reason}")]
    SpatialAlgebraError {
        /// Description of the error.
        reason: String,
    },

    /// The Jacobian matrix is singular or near-singular, making the system
    /// underdetermined.
    #[error("Jacobian is singular at dof {dof} (condition number {condition:.3e})")]
    JacobianSingular {
        /// Degree of freedom index where singularity was detected.
        dof: usize,
        /// Condition number of the Jacobian at the failure point.
        condition: f64,
    },

    /// Forward dynamics computation failed for a specific body in the chain.
    #[error("forward dynamics failed for body {body_index}: {reason}")]
    ForwardDynamicsFailure {
        /// Index of the body in the articulated chain where failure occurred.
        body_index: usize,
        /// Reason for failure.
        reason: String,
    },

    // ── Fluid coupling ─────────────────────────────────────────────────────
    /// The requested fluid medium is incompatible with the body's properties.
    #[error("fluid medium incompatible: {reason}")]
    FluidMediumIncompatible {
        /// Reason for incompatibility.
        reason: String,
    },

    /// Buoyancy computation failed (e.g. submerged volume could not be
    /// computed for this shape).
    #[error("buoyancy computation failed: {reason}")]
    BuoyancyComputationFailed {
        /// Reason for failure.
        reason: String,
    },

    // ── Pipeline / world ───────────────────────────────────────────────────
    /// The simulation pipeline configuration is invalid.
    #[error("pipeline configuration error: {reason}")]
    PipelineConfigError {
        /// Reason the configuration is invalid.
        reason: String,
    },

    /// An island solver exceeded its allowed time budget.
    #[error("island {island_id} solver timeout after {elapsed_ms:.2} ms")]
    IslandSolverTimeout {
        /// Island that timed out.
        island_id: u32,
        /// Elapsed milliseconds before the timeout was triggered.
        elapsed_ms: f64,
    },

    /// A body that was frozen (kinematic/static) was referenced as dynamic
    /// during a simulation step.
    #[error("body handle {handle:?} is frozen but was referenced as dynamic during step")]
    BodyFrozenDuringStep {
        /// Handle of the frozen body.
        handle: u64,
    },

    // ── Serialization / I/O ────────────────────────────────────────────────
    /// State snapshot deserialization failed.
    #[error("deserialization failed: {reason}")]
    DeserializationFailed {
        /// Reason for failure (e.g. "unexpected end of input").
        reason: String,
    },

    /// A version mismatch was detected during state loading.
    #[error("snapshot version mismatch: file={file_version}, engine={engine_version}")]
    VersionMismatch {
        /// Version stored in the snapshot file.
        file_version: u32,
        /// Version expected by the current engine.
        engine_version: u32,
    },
}

/// Result type alias using [`enum@Error`].
pub type Result<T> = std::result::Result<T, Error>;

// ---------------------------------------------------------------------------
// Convenience constructors
// ---------------------------------------------------------------------------

impl Error {
    /// Create a [`Error::General`] from any `Display`able value.
    pub fn general(msg: impl std::fmt::Display) -> Self {
        Error::General(msg.to_string())
    }

    /// Create a [`Error::InvalidBodyHandle`].
    pub fn invalid_handle(handle: u64) -> Self {
        Error::InvalidBodyHandle { handle }
    }

    /// Create a [`Error::NonPositiveMass`].
    pub fn non_positive_mass(mass: f64) -> Self {
        Error::NonPositiveMass { mass }
    }

    /// Create a [`Error::SingularInertiaTensor`].
    pub fn singular_inertia(det: f64) -> Self {
        Error::SingularInertiaTensor { det }
    }

    /// Create a [`Error::InvalidColliderShape`].
    pub fn invalid_shape(reason: impl Into<String>) -> Self {
        Error::InvalidColliderShape {
            reason: reason.into(),
        }
    }

    /// Create a [`Error::JointBodyNotFound`].
    pub fn joint_body_not_found(body_id: u64) -> Self {
        Error::JointBodyNotFound { body_id }
    }

    /// Create a [`Error::InvalidJointLimit`].
    pub fn invalid_joint_limit(min: f64, max: f64) -> Self {
        Error::InvalidJointLimit { min, max }
    }

    /// Create a [`Error::SolverDivergence`].
    pub fn solver_diverged(iteration: u32, residual: f64) -> Self {
        Error::SolverDivergence {
            iteration,
            residual,
        }
    }

    /// Create a [`Error::SolverNotConverged`].
    pub fn solver_not_converged(max_iterations: u32, residual: f64) -> Self {
        Error::SolverNotConverged {
            max_iterations,
            residual,
        }
    }

    /// Create a [`Error::CcdToiFailed`].
    pub fn ccd_toi_failed(reason: impl Into<String>) -> Self {
        Error::CcdToiFailed {
            reason: reason.into(),
        }
    }

    /// Create a [`Error::ToiOutOfRange`].
    pub fn toi_out_of_range(toi: f64) -> Self {
        Error::ToiOutOfRange { toi }
    }

    /// Create a [`Error::IslandOutOfRange`].
    pub fn island_out_of_range(id: u32, total: u32) -> Self {
        Error::IslandOutOfRange { id, total }
    }

    /// Create a [`Error::InvalidMotorParameter`].
    pub fn invalid_motor_param(param: &'static str, value: f64) -> Self {
        Error::InvalidMotorParameter { param, value }
    }

    /// Create a [`Error::MotorThermalLimit`].
    pub fn motor_thermal_limit(temp: f64, limit: f64) -> Self {
        Error::MotorThermalLimit { temp, limit }
    }

    /// Create a [`Error::InvalidFluidProperty`].
    pub fn invalid_fluid_prop(prop: &'static str, value: f64) -> Self {
        Error::InvalidFluidProperty { prop, value }
    }

    /// Create a [`Error::DeserializationFailed`].
    pub fn deserialization_failed(reason: impl Into<String>) -> Self {
        Error::DeserializationFailed {
            reason: reason.into(),
        }
    }

    /// Create a [`Error::VersionMismatch`].
    pub fn version_mismatch(file_version: u32, engine_version: u32) -> Self {
        Error::VersionMismatch {
            file_version,
            engine_version,
        }
    }

    /// Create a [`Error::ArticulatedBodyError`].
    pub fn articulated_body(reason: impl Into<String>) -> Self {
        Error::ArticulatedBodyError {
            reason: reason.into(),
        }
    }

    /// Create a [`Error::SpatialAlgebraError`].
    pub fn spatial_algebra(reason: impl Into<String>) -> Self {
        Error::SpatialAlgebraError {
            reason: reason.into(),
        }
    }

    /// Create a [`Error::JacobianSingular`].
    pub fn jacobian_singular(dof: usize, condition: f64) -> Self {
        Error::JacobianSingular { dof, condition }
    }

    /// Create a [`Error::ForwardDynamicsFailure`].
    pub fn forward_dynamics_failure(body_index: usize, reason: impl Into<String>) -> Self {
        Error::ForwardDynamicsFailure {
            body_index,
            reason: reason.into(),
        }
    }

    /// Create a [`Error::FluidMediumIncompatible`].
    pub fn fluid_medium_incompatible(reason: impl Into<String>) -> Self {
        Error::FluidMediumIncompatible {
            reason: reason.into(),
        }
    }

    /// Create a [`Error::BuoyancyComputationFailed`].
    pub fn buoyancy_failed(reason: impl Into<String>) -> Self {
        Error::BuoyancyComputationFailed {
            reason: reason.into(),
        }
    }

    /// Create a [`Error::PipelineConfigError`].
    pub fn pipeline_config(reason: impl Into<String>) -> Self {
        Error::PipelineConfigError {
            reason: reason.into(),
        }
    }

    /// Create a [`Error::IslandSolverTimeout`].
    pub fn island_solver_timeout(island_id: u32, elapsed_ms: f64) -> Self {
        Error::IslandSolverTimeout {
            island_id,
            elapsed_ms,
        }
    }

    /// Create a [`Error::BodyFrozenDuringStep`].
    pub fn body_frozen(handle: u64) -> Self {
        Error::BodyFrozenDuringStep { handle }
    }

    // ── Predicates ──────────────────────────────────────────────────────────

    /// Returns `true` if this error is a solver-related failure.
    pub fn is_solver_error(&self) -> bool {
        matches!(
            self,
            Error::SolverDivergence { .. } | Error::SolverNotConverged { .. }
        )
    }

    /// Returns `true` if this error is a thermal limit exceeded event.
    pub fn is_thermal(&self) -> bool {
        matches!(self, Error::MotorThermalLimit { .. })
    }

    /// Returns `true` if this error arises from an invalid handle.
    pub fn is_invalid_handle(&self) -> bool {
        matches!(
            self,
            Error::InvalidBodyHandle { .. } | Error::ColliderNotFound { .. }
        )
    }

    /// Returns `true` if this error is a CCD/TOI failure.
    pub fn is_ccd_error(&self) -> bool {
        matches!(
            self,
            Error::CcdToiFailed { .. } | Error::ToiOutOfRange { .. }
        )
    }

    /// Returns `true` if this is a data/I/O error.
    pub fn is_io_error(&self) -> bool {
        matches!(
            self,
            Error::DeserializationFailed { .. } | Error::VersionMismatch { .. }
        )
    }

    /// Returns `true` if this is an articulated-body or spatial-algebra error.
    pub fn is_articulated_error(&self) -> bool {
        matches!(
            self,
            Error::ArticulatedBodyError { .. }
                | Error::SpatialAlgebraError { .. }
                | Error::JacobianSingular { .. }
                | Error::ForwardDynamicsFailure { .. }
        )
    }

    /// Returns `true` if this is a fluid-coupling error.
    pub fn is_fluid_error(&self) -> bool {
        matches!(
            self,
            Error::InvalidFluidProperty { .. }
                | Error::FluidMediumIncompatible { .. }
                | Error::BuoyancyComputationFailed { .. }
        )
    }

    /// Returns `true` if this is a pipeline or world configuration error.
    pub fn is_pipeline_error(&self) -> bool {
        matches!(
            self,
            Error::PipelineConfigError { .. }
                | Error::IslandSolverTimeout { .. }
                | Error::BodyFrozenDuringStep { .. }
        )
    }

    /// Map this error to a severity level: 0 = recoverable, 1 = warning,
    /// 2 = fatal.
    ///
    /// This is a convenience for downstream consumers that want to log or
    /// filter by severity without exhaustively matching the enum.
    pub fn severity(&self) -> u8 {
        match self {
            // Recoverable / expected conditions
            Error::SolverNotConverged { .. }
            | Error::IslandOutOfRange { .. }
            | Error::MotorThermalLimit { .. }
            | Error::BodyFrozenDuringStep { .. }
            | Error::IslandSolverTimeout { .. } => 1,

            // Fatal / programmer errors
            Error::SolverDivergence { .. }
            | Error::SingularInertiaTensor { .. }
            | Error::JacobianSingular { .. }
            | Error::ForwardDynamicsFailure { .. } => 2,

            // Warnings (bad inputs / config)
            _ => 0,
        }
    }
}

// ---------------------------------------------------------------------------
// From impls
// ---------------------------------------------------------------------------

impl From<String> for Error {
    fn from(s: String) -> Self {
        Error::General(s)
    }
}

impl From<&str> for Error {
    fn from(s: &str) -> Self {
        Error::General(s.to_string())
    }
}

// ---------------------------------------------------------------------------
// Validation helpers (free functions)
// ---------------------------------------------------------------------------

/// Validate that `mass` is strictly positive; return an error otherwise.
pub fn require_positive_mass(mass: f64) -> Result<()> {
    if mass > 0.0 {
        Ok(())
    } else {
        Err(Error::non_positive_mass(mass))
    }
}

/// Validate that a joint limit `(min, max)` satisfies `min <= max`.
pub fn require_valid_joint_limit(min: f64, max: f64) -> Result<()> {
    if min <= max {
        Ok(())
    } else {
        Err(Error::invalid_joint_limit(min, max))
    }
}

/// Validate that a TOI is in `[0.0, 1.0]`.
pub fn require_valid_toi(toi: f64) -> Result<()> {
    if (0.0..=1.0).contains(&toi) {
        Ok(())
    } else {
        Err(Error::toi_out_of_range(toi))
    }
}

/// Validate that `value` is strictly positive for a fluid property.
pub fn require_positive_fluid_prop(prop: &'static str, value: f64) -> Result<()> {
    if value > 0.0 {
        Ok(())
    } else {
        Err(Error::invalid_fluid_prop(prop, value))
    }
}

/// Validate that `value` is strictly positive for a motor parameter.
pub fn require_positive_motor_param(param: &'static str, value: f64) -> Result<()> {
    if value > 0.0 {
        Ok(())
    } else {
        Err(Error::invalid_motor_param(param, value))
    }
}

/// Validate that a condition number is below `max_condition` for Jacobian
/// non-singularity.
pub fn require_nonsingular_jacobian(dof: usize, condition: f64, max_condition: f64) -> Result<()> {
    if condition < max_condition {
        Ok(())
    } else {
        Err(Error::jacobian_singular(dof, condition))
    }
}

/// Validate that a fluid density is positive and within a physically
/// reasonable upper bound (10× water = 10 000 kg/m³).
pub fn require_valid_fluid_density(density: f64) -> Result<()> {
    if density <= 0.0 {
        return Err(Error::invalid_fluid_prop("density", density));
    }
    if density > 10_000.0 {
        return Err(Error::fluid_medium_incompatible(format!(
            "density {density} kg/m³ exceeds maximum supported value of 10 000 kg/m³"
        )));
    }
    Ok(())
}

/// Validate that the inertia tensor determinant is above a threshold,
/// indicating it is non-singular.
pub fn require_nonsingular_inertia(det: f64, threshold: f64) -> Result<()> {
    if det.abs() > threshold {
        Ok(())
    } else {
        Err(Error::singular_inertia(det))
    }
}

/// Validate that `island_id < total_islands`.
pub fn require_valid_island(id: u32, total: u32) -> Result<()> {
    if id < total {
        Ok(())
    } else {
        Err(Error::island_out_of_range(id, total))
    }
}

// ---------------------------------------------------------------------------
// Error context helper
// ---------------------------------------------------------------------------

/// A thin wrapper that annotates an existing [`enum@Error`] with a context string,
/// useful for wrapping low-level errors with higher-level call-site information.
#[derive(Debug)]
pub struct ErrorContext {
    /// The underlying error.
    pub source: Error,
    /// Human-readable context message.
    pub context: String,
}

impl std::fmt::Display for ErrorContext {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}: {}", self.context, self.source)
    }
}

impl std::error::Error for ErrorContext {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        Some(&self.source)
    }
}

/// Extension trait that adds `.context(msg)` to \[`Result<T, Error>`\].
pub trait ErrorExt<T> {
    /// Wrap the error with a context message.
    fn context(self, ctx: impl Into<String>) -> std::result::Result<T, ErrorContext>;
}

impl<T> ErrorExt<T> for Result<T> {
    fn context(self, ctx: impl Into<String>) -> std::result::Result<T, ErrorContext> {
        self.map_err(|e| ErrorContext {
            source: e,
            context: ctx.into(),
        })
    }
}

// ---------------------------------------------------------------------------
// Error severity classification
// ---------------------------------------------------------------------------

/// Severity of a physics error.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum ErrorSeverity {
    /// The condition is an expected, recoverable situation.
    Recoverable = 0,
    /// The condition is a warning: simulation may proceed but results are
    /// degraded.
    Warning = 1,
    /// The condition is a fatal error: simulation must stop.
    Fatal = 2,
}

impl From<u8> for ErrorSeverity {
    fn from(v: u8) -> Self {
        match v {
            0 => ErrorSeverity::Recoverable,
            1 => ErrorSeverity::Warning,
            _ => ErrorSeverity::Fatal,
        }
    }
}

impl Error {
    /// Classify this error into a [`ErrorSeverity`].
    pub fn classify(&self) -> ErrorSeverity {
        ErrorSeverity::from(self.severity())
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_general_error_display() {
        let e = Error::general("something went wrong");
        let s = format!("{e}");
        assert!(s.contains("something went wrong"), "display: {s}");
    }

    #[test]
    fn test_non_positive_mass_error() {
        let e = Error::non_positive_mass(-1.0);
        assert!(matches!(e, Error::NonPositiveMass { mass } if mass == -1.0));
        let s = format!("{e}");
        assert!(s.contains("positive"), "display: {s}");
    }

    #[test]
    fn test_solver_diverged_is_solver_error() {
        let e = Error::solver_diverged(5, 1e8);
        assert!(e.is_solver_error());
        assert!(!e.is_thermal());
    }

    #[test]
    fn test_solver_not_converged_is_solver_error() {
        let e = Error::solver_not_converged(100, 1e-2);
        assert!(e.is_solver_error());
        let s = format!("{e}");
        assert!(s.contains("100"), "display: {s}");
    }

    #[test]
    fn test_motor_thermal_limit_predicate() {
        let e = Error::motor_thermal_limit(195.0, 180.0);
        assert!(e.is_thermal());
        assert!(!e.is_solver_error());
        let s = format!("{e}");
        assert!(s.contains("195"), "display: {s}");
    }

    #[test]
    fn test_invalid_handle_predicate() {
        let e = Error::invalid_handle(42);
        assert!(e.is_invalid_handle());
        assert!(!e.is_io_error());
    }

    #[test]
    fn test_ccd_error_predicate() {
        let e = Error::toi_out_of_range(1.5);
        assert!(e.is_ccd_error());
        let s = format!("{e}");
        assert!(s.contains("1.5"), "display: {s}");
    }

    #[test]
    fn test_version_mismatch_is_io_error() {
        let e = Error::version_mismatch(1, 2);
        assert!(e.is_io_error());
        let s = format!("{e}");
        assert!(s.contains("file=1"), "display: {s}");
    }

    #[test]
    fn test_invalid_joint_limit_display() {
        let e = Error::invalid_joint_limit(1.0, -1.0);
        let s = format!("{e}");
        assert!(s.contains("min=1"), "display: {s}");
        assert!(s.contains("max=-1"), "display: {s}");
    }

    #[test]
    fn test_require_positive_mass_ok() {
        assert!(require_positive_mass(1.0).is_ok());
    }

    #[test]
    fn test_require_positive_mass_err() {
        let r = require_positive_mass(0.0);
        assert!(r.is_err());
    }

    #[test]
    fn test_require_valid_joint_limit_ok() {
        assert!(require_valid_joint_limit(-1.0, 1.0).is_ok());
        assert!(require_valid_joint_limit(0.0, 0.0).is_ok()); // equal is valid
    }

    #[test]
    fn test_require_valid_joint_limit_err() {
        assert!(require_valid_joint_limit(1.0, -1.0).is_err());
    }

    #[test]
    fn test_require_valid_toi_ok() {
        assert!(require_valid_toi(0.0).is_ok());
        assert!(require_valid_toi(0.5).is_ok());
        assert!(require_valid_toi(1.0).is_ok());
    }

    #[test]
    fn test_require_valid_toi_err() {
        assert!(require_valid_toi(-0.01).is_err());
        assert!(require_valid_toi(1.001).is_err());
    }

    #[test]
    fn test_require_positive_fluid_prop() {
        assert!(require_positive_fluid_prop("density", 1000.0).is_ok());
        assert!(require_positive_fluid_prop("density", -1.0).is_err());
    }

    #[test]
    fn test_require_positive_motor_param() {
        assert!(require_positive_motor_param("resistance", 0.5).is_ok());
        assert!(require_positive_motor_param("resistance", 0.0).is_err());
    }

    #[test]
    fn test_from_string() {
        let e: Error = String::from("boom").into();
        assert!(matches!(e, Error::General(_)));
    }

    #[test]
    fn test_from_str() {
        let e: Error = "boom".into();
        assert!(matches!(e, Error::General(_)));
    }

    #[test]
    fn test_island_out_of_range_display() {
        let e = Error::island_out_of_range(5, 3);
        let s = format!("{e}");
        assert!(s.contains("5"), "display: {s}");
        assert!(s.contains("3"), "display: {s}");
    }

    #[test]
    fn test_singular_inertia_display() {
        let e = Error::singular_inertia(1e-20);
        let s = format!("{e}");
        assert!(s.contains("singular"), "display: {s}");
    }

    #[test]
    fn test_invalid_shape_display() {
        let e = Error::invalid_shape("negative radius");
        let s = format!("{e}");
        assert!(s.contains("negative radius"), "display: {s}");
    }

    #[test]
    fn test_deserialization_failed_display() {
        let e = Error::deserialization_failed("unexpected EOF");
        assert!(e.is_io_error());
        let s = format!("{e}");
        assert!(s.contains("EOF"), "display: {s}");
    }

    #[test]
    fn test_ccd_toi_failed_display() {
        let e = Error::ccd_toi_failed("shape not supported");
        assert!(e.is_ccd_error());
        let s = format!("{e}");
        assert!(s.contains("shape not supported"), "display: {s}");
    }

    // ── New variant tests ─────────────────────────────────────────────────

    #[test]
    fn test_articulated_body_error_display() {
        let e = Error::articulated_body("loop closure failed");
        assert!(e.is_articulated_error());
        assert!(!e.is_solver_error());
        let s = format!("{e}");
        assert!(s.contains("loop closure"), "display: {s}");
    }

    #[test]
    fn test_spatial_algebra_error_display() {
        let e = Error::spatial_algebra("zero-length joint axis");
        assert!(e.is_articulated_error());
        let s = format!("{e}");
        assert!(s.contains("zero-length"), "display: {s}");
    }

    #[test]
    fn test_jacobian_singular_display() {
        let e = Error::jacobian_singular(3, 1e15);
        assert!(e.is_articulated_error());
        let s = format!("{e}");
        assert!(s.contains("3"), "display: {s}");
    }

    #[test]
    fn test_forward_dynamics_failure_display() {
        let e = Error::forward_dynamics_failure(2, "null spatial inertia");
        assert!(e.is_articulated_error());
        let s = format!("{e}");
        assert!(s.contains("null spatial inertia"), "display: {s}");
    }

    #[test]
    fn test_fluid_medium_incompatible_display() {
        let e = Error::fluid_medium_incompatible("supercritical phase not supported");
        assert!(e.is_fluid_error());
        let s = format!("{e}");
        assert!(s.contains("supercritical"), "display: {s}");
    }

    #[test]
    fn test_buoyancy_failed_display() {
        let e = Error::buoyancy_failed("concave mesh not supported");
        assert!(e.is_fluid_error());
        let s = format!("{e}");
        assert!(s.contains("concave"), "display: {s}");
    }

    #[test]
    fn test_pipeline_config_error_display() {
        let e = Error::pipeline_config("dt must be positive");
        assert!(e.is_pipeline_error());
        let s = format!("{e}");
        assert!(s.contains("dt must be positive"), "display: {s}");
    }

    #[test]
    fn test_island_solver_timeout_display() {
        let e = Error::island_solver_timeout(4, 32.5);
        assert!(e.is_pipeline_error());
        let s = format!("{e}");
        assert!(s.contains("32.50"), "display: {s}");
    }

    #[test]
    fn test_body_frozen_during_step() {
        let e = Error::body_frozen(7);
        assert!(e.is_pipeline_error());
        let s = format!("{e}");
        assert!(s.contains("7"), "display: {s}");
    }

    // ── Validation helpers ────────────────────────────────────────────────

    #[test]
    fn test_require_nonsingular_jacobian_ok() {
        assert!(require_nonsingular_jacobian(0, 10.0, 1e6).is_ok());
    }

    #[test]
    fn test_require_nonsingular_jacobian_err() {
        let r = require_nonsingular_jacobian(1, 2e10, 1e6);
        assert!(r.is_err());
        assert!(matches!(
            r.unwrap_err(),
            Error::JacobianSingular { dof: 1, .. }
        ));
    }

    #[test]
    fn test_require_valid_fluid_density_ok() {
        assert!(require_valid_fluid_density(1000.0).is_ok());
        assert!(require_valid_fluid_density(1.0).is_ok());
    }

    #[test]
    fn test_require_valid_fluid_density_err_zero() {
        assert!(require_valid_fluid_density(0.0).is_err());
    }

    #[test]
    fn test_require_valid_fluid_density_err_too_large() {
        let r = require_valid_fluid_density(50_000.0);
        assert!(r.is_err());
        assert!(matches!(
            r.unwrap_err(),
            Error::FluidMediumIncompatible { .. }
        ));
    }

    #[test]
    fn test_require_nonsingular_inertia_ok() {
        assert!(require_nonsingular_inertia(1e-3, 1e-10).is_ok());
    }

    #[test]
    fn test_require_nonsingular_inertia_err() {
        assert!(require_nonsingular_inertia(1e-20, 1e-10).is_err());
    }

    #[test]
    fn test_require_valid_island_ok() {
        assert!(require_valid_island(0, 5).is_ok());
        assert!(require_valid_island(4, 5).is_ok());
    }

    #[test]
    fn test_require_valid_island_err() {
        assert!(require_valid_island(5, 5).is_err());
        assert!(require_valid_island(100, 5).is_err());
    }

    // ── ErrorContext ──────────────────────────────────────────────────────

    #[test]
    fn test_error_context_wraps_error() {
        let result: Result<()> = Err(Error::non_positive_mass(-5.0));
        let wrapped = result.context("while setting up body mass");
        let err = wrapped.unwrap_err();
        let s = format!("{err}");
        assert!(s.contains("while setting up body mass"), "display: {s}");
        assert!(s.contains("positive"), "display: {s}");
    }

    #[test]
    fn test_error_context_source() {
        use std::error::Error as StdError;
        let ctx = ErrorContext {
            source: Error::non_positive_mass(-1.0),
            context: "test context".into(),
        };
        assert!(ctx.source().is_some());
    }

    // ── ErrorSeverity ─────────────────────────────────────────────────────

    #[test]
    fn test_severity_recoverable() {
        let e = Error::solver_not_converged(100, 1e-2);
        assert_eq!(e.classify(), ErrorSeverity::Warning);
    }

    #[test]
    fn test_severity_fatal_divergence() {
        let e = Error::solver_diverged(5, 1e12);
        assert_eq!(e.classify(), ErrorSeverity::Fatal);
    }

    #[test]
    fn test_severity_jacobian_singular() {
        let e = Error::jacobian_singular(2, 1e18);
        assert_eq!(e.classify(), ErrorSeverity::Fatal);
    }

    #[test]
    fn test_severity_recoverable_general() {
        let e = Error::general("minor issue");
        assert_eq!(e.classify(), ErrorSeverity::Recoverable);
    }

    #[test]
    fn test_error_severity_ordering() {
        assert!(ErrorSeverity::Recoverable < ErrorSeverity::Warning);
        assert!(ErrorSeverity::Warning < ErrorSeverity::Fatal);
    }
}
