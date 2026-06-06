// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Error types for oxiphysics-materials.
//!
//! Provides a rich error taxonomy covering material parameter validation,
//! numerical issues, EOS failures, phase-transform convergence, and I/O.

use thiserror::Error;

/// Main error type for the materials module.
#[derive(Debug, Error)]
pub enum Error {
    /// Generic error message.
    #[error("{0}")]
    General(String),

    /// A required material parameter is out of its physical range.
    ///
    /// `param` is the parameter name, `value` is what was supplied,
    /// `min` / `max` are the valid bounds (use ±`f64::INFINITY` for open ends).
    #[error("parameter '{param}' = {value} is out of range [{min}, {max}]")]
    ParameterOutOfRange {
        /// Parameter name.
        param: &'static str,
        /// Supplied value.
        value: f64,
        /// Minimum allowed value.
        min: f64,
        /// Maximum allowed value.
        max: f64,
    },

    /// A required material property was not set.
    #[error("required material property '{0}' is missing")]
    MissingProperty(&'static str),

    /// Numerical divergence during an iterative algorithm.
    ///
    /// Typically triggered when an implicit solver does not converge.
    #[error(
        "numerical divergence in '{solver}' after {iterations} iterations (residual {residual:.3e})"
    )]
    NumericalDivergence {
        /// Name of the solver or algorithm that diverged.
        solver: &'static str,
        /// Number of iterations performed.
        iterations: usize,
        /// Final residual or error norm.
        residual: f64,
    },

    /// Singularity encountered (e.g. zero determinant, zero density).
    #[error("singularity in '{context}': {detail}")]
    Singularity {
        /// Where the singularity occurred.
        context: &'static str,
        /// Human-readable detail.
        detail: String,
    },

    /// The requested deformation state is physically inadmissible
    /// (e.g. negative volume, compressive stretch beyond locking).
    #[error("inadmissible deformation state: {0}")]
    InadmissibleState(String),

    /// An equation of state returned an unphysical pressure.
    #[error(
        "EOS '{eos}' returned unphysical pressure {pressure:.3e} Pa at rho={density:.3e} kg/m³"
    )]
    UnphysicalEosPressure {
        /// Equation of state identifier.
        eos: &'static str,
        /// Computed pressure (Pa).
        pressure: f64,
        /// Density at which the EOS was evaluated (kg/m³).
        density: f64,
    },

    /// Phase-transformation model failed to converge.
    #[error("phase transform '{model}' did not converge: {detail}")]
    PhaseTransformConvergence {
        /// Model identifier (e.g. "JMAK", "martensitic").
        model: &'static str,
        /// Additional diagnostic information.
        detail: String,
    },

    /// Table look-up out of bounds (e.g. tabulated EOS, TTT diagram).
    #[error(
        "table look-up for '{table}' out of range: {variable} = {value:.3e} (range [{lo:.3e}, {hi:.3e}])"
    )]
    TableOutOfRange {
        /// Table identifier.
        table: &'static str,
        /// Variable name being looked up.
        variable: &'static str,
        /// Queried value.
        value: f64,
        /// Lower bound of table.
        lo: f64,
        /// Upper bound of table.
        hi: f64,
    },

    /// Fatigue model error (e.g. invalid cycle count, negative stress amplitude).
    #[error("fatigue model error in '{model}': {detail}")]
    FatigueModel {
        /// Model name.
        model: &'static str,
        /// Description.
        detail: String,
    },

    /// Fracture mechanics error (e.g. stress intensity factor not finite).
    #[error("fracture mechanics error: {0}")]
    FractureMechanics(String),

    /// Incompatible units or dimension mismatch.
    #[error("unit/dimension error: expected '{expected}', got '{actual}'")]
    DimensionMismatch {
        /// Expected dimension string.
        expected: String,
        /// Actual dimension string.
        actual: String,
    },

    /// I/O error when loading material data from a file.
    #[error("I/O error loading material data from '{path}': {message}")]
    Io {
        /// File path.
        path: String,
        /// Error message.
        message: String,
    },

    /// Parse error when decoding material data.
    #[error("parse error in '{context}': {message}")]
    Parse {
        /// Parsing context.
        context: String,
        /// Error message.
        message: String,
    },
}

/// Result type alias for the materials module.
pub type Result<T> = std::result::Result<T, Error>;

// ─────────────────────────────────────────────────────────────────────────────
// Convenience constructors
// ─────────────────────────────────────────────────────────────────────────────

impl Error {
    /// Create a [`Error::ParameterOutOfRange`] for a lower-bound violation.
    pub fn below_minimum(param: &'static str, value: f64, min: f64) -> Self {
        Self::ParameterOutOfRange {
            param,
            value,
            min,
            max: f64::INFINITY,
        }
    }

    /// Create a [`Error::ParameterOutOfRange`] for an upper-bound violation.
    pub fn above_maximum(param: &'static str, value: f64, max: f64) -> Self {
        Self::ParameterOutOfRange {
            param,
            value,
            min: f64::NEG_INFINITY,
            max,
        }
    }

    /// Create a [`Error::NumericalDivergence`] with a simple description.
    pub fn diverged(solver: &'static str, iterations: usize, residual: f64) -> Self {
        Self::NumericalDivergence {
            solver,
            iterations,
            residual,
        }
    }

    /// Create a [`Error::Singularity`] error.
    pub fn singular(context: &'static str, detail: impl Into<String>) -> Self {
        Self::Singularity {
            context,
            detail: detail.into(),
        }
    }

    /// Create an [`Error::InadmissibleState`] error.
    pub fn inadmissible(detail: impl Into<String>) -> Self {
        Self::InadmissibleState(detail.into())
    }

    /// Create a [`Error::TableOutOfRange`] error.
    pub fn table_out_of_range(
        table: &'static str,
        variable: &'static str,
        value: f64,
        lo: f64,
        hi: f64,
    ) -> Self {
        Self::TableOutOfRange {
            table,
            variable,
            value,
            lo,
            hi,
        }
    }

    /// Create a [`Error::FatigueModel`] error.
    pub fn fatigue(model: &'static str, detail: impl Into<String>) -> Self {
        Self::FatigueModel {
            model,
            detail: detail.into(),
        }
    }

    /// Create an [`Error::UnphysicalEosPressure`] error.
    pub fn unphysical_pressure(eos: &'static str, pressure: f64, density: f64) -> Self {
        Self::UnphysicalEosPressure {
            eos,
            pressure,
            density,
        }
    }

    /// Returns `true` if this is a numerical issue (divergence or singularity).
    pub fn is_numerical(&self) -> bool {
        matches!(
            self,
            Self::NumericalDivergence { .. } | Self::Singularity { .. }
        )
    }

    /// Returns `true` if this is a parameter validation error.
    pub fn is_parameter_error(&self) -> bool {
        matches!(
            self,
            Self::ParameterOutOfRange { .. } | Self::MissingProperty(_)
        )
    }

    /// Returns `true` if this is an EOS-specific error.
    pub fn is_eos_error(&self) -> bool {
        matches!(self, Self::UnphysicalEosPressure { .. })
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Validation helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Assert that `value` is strictly positive; return `Err` otherwise.
///
/// ```no_run
/// use oxiphysics_materials::require_positive;
/// assert!(require_positive("density", 1000.0).is_ok());
/// assert!(require_positive("density", -1.0).is_err());
/// ```
pub fn require_positive(param: &'static str, value: f64) -> Result<f64> {
    if value > 0.0 {
        Ok(value)
    } else {
        Err(Error::below_minimum(param, value, 0.0))
    }
}

/// Assert that `value` is non-negative; return `Err` otherwise.
///
/// ```no_run
/// use oxiphysics_materials::require_non_negative;
/// assert!(require_non_negative("strain", 0.0).is_ok());
/// assert!(require_non_negative("strain", -0.1).is_err());
/// ```
pub fn require_non_negative(param: &'static str, value: f64) -> Result<f64> {
    if value >= 0.0 {
        Ok(value)
    } else {
        Err(Error::below_minimum(param, value, 0.0))
    }
}

/// Assert that `value` lies in the closed interval `[lo, hi]`.
///
/// ```no_run
/// use oxiphysics_materials::require_in_range;
/// assert!(require_in_range("nu", 0.3, 0.0, 0.5).is_ok());
/// assert!(require_in_range("nu", 0.6, 0.0, 0.5).is_err());
/// ```
pub fn require_in_range(param: &'static str, value: f64, lo: f64, hi: f64) -> Result<f64> {
    if value >= lo && value <= hi {
        Ok(value)
    } else {
        Err(Error::ParameterOutOfRange {
            param,
            value,
            min: lo,
            max: hi,
        })
    }
}

/// Assert that `value` is finite (not NaN or ±infinity).
pub fn require_finite(param: &'static str, value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(Error::General(format!(
            "parameter '{param}' is not finite: {value}"
        )))
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── Display / formatting ──────────────────────────────────────────────────

    #[test]
    fn test_general_error_display() {
        let e = Error::General("test message".to_string());
        assert_eq!(format!("{e}"), "test message");
    }

    #[test]
    fn test_parameter_out_of_range_display() {
        let e = Error::ParameterOutOfRange {
            param: "density",
            value: -1.0,
            min: 0.0,
            max: f64::INFINITY,
        };
        let s = format!("{e}");
        assert!(s.contains("density"), "display: {s}");
        assert!(s.contains("-1"), "display: {s}");
    }

    #[test]
    fn test_numerical_divergence_display() {
        let e = Error::diverged("Newton", 100, 1e-3);
        let s = format!("{e}");
        assert!(s.contains("Newton"), "display: {s}");
        assert!(s.contains("100"), "display: {s}");
    }

    #[test]
    fn test_singularity_display() {
        let e = Error::singular("det(F)", "F is rank-deficient");
        let s = format!("{e}");
        assert!(s.contains("det(F)"), "display: {s}");
    }

    #[test]
    fn test_inadmissible_state_display() {
        let e = Error::inadmissible("negative volume");
        let s = format!("{e}");
        assert!(s.contains("negative volume"), "display: {s}");
    }

    #[test]
    fn test_unphysical_eos_pressure_display() {
        let e = Error::unphysical_pressure("IdealGas", -1e10, 1.0);
        let s = format!("{e}");
        assert!(s.contains("IdealGas"), "display: {s}");
        assert!(s.contains("pressure"), "display: {s}");
    }

    #[test]
    fn test_phase_transform_convergence_display() {
        let e = Error::PhaseTransformConvergence {
            model: "JMAK",
            detail: "temperature out of range".to_string(),
        };
        let s = format!("{e}");
        assert!(s.contains("JMAK"), "display: {s}");
    }

    #[test]
    fn test_table_out_of_range_display() {
        let e = Error::table_out_of_range("TTT", "temperature", 1500.0, 300.0, 1200.0);
        let s = format!("{e}");
        assert!(s.contains("TTT"), "display: {s}");
        assert!(s.contains("temperature"), "display: {s}");
    }

    #[test]
    fn test_fatigue_error_display() {
        let e = Error::fatigue("Basquin", "negative stress amplitude");
        let s = format!("{e}");
        assert!(s.contains("Basquin"), "display: {s}");
    }

    #[test]
    fn test_fracture_mechanics_display() {
        let e = Error::FractureMechanics("K_I is NaN".to_string());
        let s = format!("{e}");
        assert!(s.contains("K_I"), "display: {s}");
    }

    // ── Convenience constructors ──────────────────────────────────────────────

    #[test]
    fn test_below_minimum() {
        let e = Error::below_minimum("E", -1.0, 0.0);
        assert!(matches!(e, Error::ParameterOutOfRange { .. }));
    }

    #[test]
    fn test_above_maximum() {
        let e = Error::above_maximum("nu", 0.6, 0.5);
        assert!(matches!(e, Error::ParameterOutOfRange { .. }));
    }

    // ── Classification helpers ────────────────────────────────────────────────

    #[test]
    fn test_is_numerical_divergence() {
        let e = Error::diverged("solver", 10, 1e-2);
        assert!(e.is_numerical());
        assert!(!e.is_parameter_error());
        assert!(!e.is_eos_error());
    }

    #[test]
    fn test_is_numerical_singularity() {
        let e = Error::singular("ctx", "detail");
        assert!(e.is_numerical());
    }

    #[test]
    fn test_is_parameter_error_out_of_range() {
        let e = Error::below_minimum("x", -1.0, 0.0);
        assert!(e.is_parameter_error());
        assert!(!e.is_numerical());
    }

    #[test]
    fn test_is_parameter_error_missing() {
        let e = Error::MissingProperty("viscosity");
        assert!(e.is_parameter_error());
    }

    #[test]
    fn test_is_eos_error() {
        let e = Error::unphysical_pressure("JWL", -1e9, 2000.0);
        assert!(e.is_eos_error());
        assert!(!e.is_numerical());
    }

    // ── Validation helpers ────────────────────────────────────────────────────

    #[test]
    fn test_require_positive_ok() {
        assert_eq!(require_positive("rho", 1000.0).unwrap(), 1000.0);
    }

    #[test]
    fn test_require_positive_fail_zero() {
        assert!(require_positive("rho", 0.0).is_err());
    }

    #[test]
    fn test_require_positive_fail_negative() {
        assert!(require_positive("E", -1.0).is_err());
    }

    #[test]
    fn test_require_non_negative_zero_ok() {
        assert_eq!(require_non_negative("strain", 0.0).unwrap(), 0.0);
    }

    #[test]
    fn test_require_non_negative_fail() {
        assert!(require_non_negative("strain", -0.001).is_err());
    }

    #[test]
    fn test_require_in_range_ok() {
        assert_eq!(require_in_range("nu", 0.3, 0.0, 0.5).unwrap(), 0.3);
    }

    #[test]
    fn test_require_in_range_boundary_ok() {
        assert!(require_in_range("nu", 0.0, 0.0, 0.5).is_ok());
        assert!(require_in_range("nu", 0.5, 0.0, 0.5).is_ok());
    }

    #[test]
    fn test_require_in_range_fail_high() {
        assert!(require_in_range("nu", 0.6, 0.0, 0.5).is_err());
    }

    #[test]
    fn test_require_in_range_fail_low() {
        assert!(require_in_range("nu", -0.1, 0.0, 0.5).is_err());
    }

    #[test]
    fn test_require_finite_ok() {
        assert_eq!(require_finite("x", 1.0).unwrap(), 1.0);
    }

    #[test]
    fn test_require_finite_nan() {
        assert!(require_finite("x", f64::NAN).is_err());
    }

    #[test]
    fn test_require_finite_inf() {
        assert!(require_finite("x", f64::INFINITY).is_err());
    }

    #[test]
    fn test_dimension_mismatch_display() {
        let e = Error::DimensionMismatch {
            expected: "Pa".to_string(),
            actual: "MPa".to_string(),
        };
        let s = format!("{e}");
        assert!(s.contains("Pa"), "display: {s}");
        assert!(s.contains("MPa"), "display: {s}");
    }

    #[test]
    fn test_io_error_display() {
        let e = Error::Io {
            path: "/tmp/mat.json".to_string(),
            message: "file not found".to_string(),
        };
        let s = format!("{e}");
        assert!(s.contains("/tmp/mat.json"), "display: {s}");
    }

    #[test]
    fn test_parse_error_display() {
        let e = Error::Parse {
            context: "JSON".to_string(),
            message: "unexpected token".to_string(),
        };
        let s = format!("{e}");
        assert!(s.contains("JSON"), "display: {s}");
    }

    // ── Result type alias ─────────────────────────────────────────────────────

    #[test]
    fn test_result_alias_ok() {
        let r: Result<f64> = Ok(2.72);
        assert!(r.is_ok());
    }

    #[test]
    fn test_result_alias_err() {
        let r: Result<f64> = Err(Error::inadmissible("test"));
        assert!(r.is_err());
    }
}
