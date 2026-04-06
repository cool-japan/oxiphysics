#![allow(clippy::manual_range_contains)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Error types for oxiphysics-lbm
//!
//! This module provides a rich, structured error hierarchy for LBM simulations.
//! All error variants include contextual messages that help diagnose the root
//! cause without needing to attach a debugger.

use thiserror::Error;

// ============================================================================
// Core Error Type
// ============================================================================

/// Main error type for the lbm module.
///
/// Each variant corresponds to a distinct failure category, making it easy
/// to match on error types and take corrective action programmatically.
#[derive(Debug, Error, PartialEq, Clone)]
pub enum Error {
    /// Generic unclassified error — use specific variants when possible.
    #[error("{0}")]
    General(String),

    /// Grid geometry is invalid (e.g., zero-size dimension).
    #[error("invalid grid dimensions: nx={nx}, ny={ny}")]
    InvalidGridDimensions {
        /// Number of grid cells in x.
        nx: usize,
        /// Number of grid cells in y.
        ny: usize,
    },

    /// 3-D grid geometry is invalid.
    #[error("invalid 3D grid dimensions: nx={nx}, ny={ny}, nz={nz}")]
    InvalidGridDimensions3D {
        /// Number of grid cells in x.
        nx: usize,
        /// Number of grid cells in y.
        ny: usize,
        /// Number of grid cells in z.
        nz: usize,
    },

    /// A simulation parameter lies outside its physical bounds.
    #[error("parameter '{name}' out of bounds: value={value}, allowed [{min}, {max}]")]
    ParameterOutOfBounds {
        /// Parameter name for diagnostics.
        name: String,
        /// Actual value supplied by the caller.
        value: f64,
        /// Inclusive lower bound.
        min: f64,
        /// Inclusive upper bound.
        max: f64,
    },

    /// The LBM relaxation time τ is too small for numerical stability.
    ///
    /// For BGK/SRT, τ must satisfy τ > 0.5.
    #[error("unstable relaxation time tau={tau:.6} (must be > 0.5 for BGK stability)")]
    UnstableRelaxationTime {
        /// The τ value that triggered the error.
        tau: f64,
    },

    /// A field (density, velocity, …) contains NaN or ±Inf.
    #[error("non-finite value detected in field '{field}' at index {index}: value={value}")]
    NonFiniteField {
        /// Descriptive name of the field (e.g., `"density"`, `"velocity_x"`).
        field: String,
        /// Linear index where the bad value was found.
        index: usize,
        /// The offending value.
        value: f64,
    },

    /// Boundary-condition specification is inconsistent or incomplete.
    #[error("boundary condition error on face '{face}': {reason}")]
    BoundaryConditionError {
        /// Which face/boundary (e.g., `"left"`, `"inlet"`).
        face: String,
        /// Human-readable explanation.
        reason: String,
    },

    /// An index is out of bounds for the given grid.
    #[error("index out of bounds: ({x}, {y}) for grid ({nx} x {ny})")]
    IndexOutOfBounds {
        /// x index.
        x: usize,
        /// y index.
        y: usize,
        /// Grid width.
        nx: usize,
        /// Grid height.
        ny: usize,
    },

    /// 3-D index out of bounds.
    #[error("3D index out of bounds: ({x}, {y}, {z}) for grid ({nx} x {ny} x {nz})")]
    IndexOutOfBounds3D {
        /// x index.
        x: usize,
        /// y index.
        y: usize,
        /// z index.
        z: usize,
        /// Grid width.
        nx: usize,
        /// Grid height.
        ny: usize,
        /// Grid depth.
        nz: usize,
    },

    /// The simulation has not converged within the allotted iterations.
    #[error(
        "convergence failure: residual={residual:.3e} after {iterations} iterations (tolerance={tolerance:.3e})"
    )]
    ConvergenceFailure {
        /// The residual at the time of the error.
        residual: f64,
        /// Number of iterations completed.
        iterations: usize,
        /// Target convergence tolerance.
        tolerance: f64,
    },

    /// Physical units or non-dimensional parameters are incompatible.
    #[error("incompatible parameters: {message}")]
    IncompatibleParameters {
        /// Explanation of the conflict.
        message: String,
    },

    /// The species count in an electrokinetic or multi-component simulation
    /// is zero or inconsistent.
    #[error("invalid species count: expected at least {min_species}, got {actual}")]
    InvalidSpeciesCount {
        /// Minimum required number of species.
        min_species: usize,
        /// Actual count provided.
        actual: usize,
    },

    /// The thermal model detected an unphysical temperature (e.g., negative).
    #[error("unphysical temperature: T={temperature:.4} K at grid index {index}")]
    UnphysicalTemperature {
        /// Temperature value.
        temperature: f64,
        /// Linear index in the temperature field.
        index: usize,
    },

    /// The Mach number exceeds the LBM compressibility limit (~0.3).
    #[error("Mach number {mach:.4} exceeds LBM limit ({limit:.2}) — simulation may diverge")]
    MachNumberExceeded {
        /// Computed Mach number.
        mach: f64,
        /// Recommended upper limit.
        limit: f64,
    },

    /// A requested velocity set (D2Q9, D3Q19, …) is not supported.
    #[error("unsupported velocity set: '{name}'")]
    UnsupportedVelocitySet {
        /// Name of the requested velocity set.
        name: String,
    },

    /// The Reynolds number is too high for the given grid resolution.
    ///
    /// Fine-scale turbulence requires `Re < Re_max` for the chosen grid.
    #[error("Reynolds number {re:.1} exceeds maximum {re_max:.1} for grid resolution {grid_res}")]
    ReynoldsNumberTooHigh {
        /// Computed Reynolds number.
        re: f64,
        /// Maximum allowable Reynolds number for this grid resolution.
        re_max: f64,
        /// Grid resolution measure (e.g. minimum cell size or node count).
        grid_res: usize,
    },

    /// The Knudsen number is outside the valid continuum range.
    ///
    /// LBM is valid for Kn < ~0.1; above this rarefaction effects dominate.
    #[error("Knudsen number {kn:.4} out of continuum range [0, {kn_max:.2}]")]
    KnudsenNumberOutOfRange {
        /// Knudsen number Kn = λ / L.
        kn: f64,
        /// Maximum Knudsen number for which the BGK approximation holds.
        kn_max: f64,
    },

    /// A multi-component simulation has inconsistent component counts.
    #[error("component count mismatch: expected {expected}, found {found} in field '{field}'")]
    ComponentCountMismatch {
        /// Expected number of components.
        expected: usize,
        /// Actual number of components found.
        found: usize,
        /// Name of the field where the mismatch was detected.
        field: String,
    },

    /// A distribution function contains a negative value, which is unphysical.
    #[error("negative distribution f[{direction}] = {value:.6} at index {index}")]
    NegativeDistribution {
        /// Velocity-set direction index.
        direction: usize,
        /// Linear grid index.
        index: usize,
        /// The offending distribution value.
        value: f64,
    },

    /// The time step `Δt` violates the CFL-like stability condition.
    #[error("time step dt={dt:.3e} violates CFL: dt_max={dt_max:.3e}")]
    TimeStepTooLarge {
        /// Requested time step.
        dt: f64,
        /// Maximum allowable time step for stability.
        dt_max: f64,
    },

    /// A turbulence model parameter is invalid.
    #[error("invalid turbulence model parameter '{name}': {reason}")]
    InvalidTurbulenceParameter {
        /// Parameter name.
        name: String,
        /// Reason the parameter is invalid.
        reason: String,
    },

    /// Wall distance `y` (in wall units y⁺) is negative or zero.
    #[error("non-positive wall distance y_plus={y_plus:.4} at index {index}")]
    NonPositiveWallDistance {
        /// Wall-normal distance in viscous units.
        y_plus: f64,
        /// Grid index.
        index: usize,
    },

    /// A phase-field order parameter is outside the valid range \[-1, 1\].
    #[error("phase-field order parameter phi={phi:.4} out of range [-1, 1] at index {index}")]
    PhaseFieldOutOfRange {
        /// Order parameter value.
        phi: f64,
        /// Grid index.
        index: usize,
    },

    /// The immersed boundary Lagrangian marker count is zero or inconsistent.
    #[error("invalid IBM marker count: expected at least {min_markers}, got {actual}")]
    InvalidIbmMarkerCount {
        /// Minimum required number of Lagrangian markers.
        min_markers: usize,
        /// Actual count.
        actual: usize,
    },

    /// A particle coupling simulation has an invalid particle configuration.
    #[error("invalid particle configuration: {reason}")]
    InvalidParticleConfiguration {
        /// Human-readable explanation.
        reason: String,
    },
}

/// Result type alias for the lbm module.
pub type Result<T> = std::result::Result<T, Error>;

// ============================================================================
// Convenience constructors
// ============================================================================

impl Error {
    /// Construct a [`General`](Error::General) error from any displayable type.
    #[allow(dead_code)]
    pub fn general(msg: impl std::fmt::Display) -> Self {
        Self::General(msg.to_string())
    }

    /// Construct an [`InvalidGridDimensions`](Error::InvalidGridDimensions) error.
    #[allow(dead_code)]
    pub fn invalid_grid(nx: usize, ny: usize) -> Self {
        Self::InvalidGridDimensions { nx, ny }
    }

    /// Construct an [`UnstableRelaxationTime`](Error::UnstableRelaxationTime) error.
    #[allow(dead_code)]
    pub fn unstable_tau(tau: f64) -> Self {
        Self::UnstableRelaxationTime { tau }
    }

    /// Check that τ is numerically stable; return `Err` if not.
    ///
    /// For BGK LBM, stability requires τ > 0.5.
    #[allow(dead_code)]
    pub fn check_tau(tau: f64) -> Result<()> {
        if tau <= 0.5 {
            Err(Self::UnstableRelaxationTime { tau })
        } else {
            Ok(())
        }
    }

    /// Check that a 2-D grid size is valid (both dimensions > 0).
    #[allow(dead_code)]
    pub fn check_grid(nx: usize, ny: usize) -> Result<()> {
        if nx == 0 || ny == 0 {
            Err(Self::InvalidGridDimensions { nx, ny })
        } else {
            Ok(())
        }
    }

    /// Check that a value is finite; return `NonFiniteField` if not.
    #[allow(dead_code)]
    pub fn check_finite(value: f64, field: &str, index: usize) -> Result<()> {
        if value.is_finite() {
            Ok(())
        } else {
            Err(Self::NonFiniteField {
                field: field.to_string(),
                index,
                value,
            })
        }
    }

    /// Check that a value lies within `[min, max]`.
    #[allow(dead_code)]
    pub fn check_range(value: f64, min: f64, max: f64, name: &str) -> Result<()> {
        if value >= min && value <= max {
            Ok(())
        } else {
            Err(Self::ParameterOutOfBounds {
                name: name.to_string(),
                value,
                min,
                max,
            })
        }
    }

    /// Check that a temperature is physically meaningful (> 0 K).
    #[allow(dead_code)]
    pub fn check_temperature(temperature: f64, index: usize) -> Result<()> {
        if temperature > 0.0 {
            Ok(())
        } else {
            Err(Self::UnphysicalTemperature { temperature, index })
        }
    }

    /// Check that the Mach number stays below the LBM incompressibility limit.
    ///
    /// The canonical LBM limit is Ma < 0.3; this function uses that default.
    #[allow(dead_code)]
    pub fn check_mach(mach: f64) -> Result<()> {
        const LIMIT: f64 = 0.3;
        if mach < LIMIT {
            Ok(())
        } else {
            Err(Self::MachNumberExceeded { mach, limit: LIMIT })
        }
    }

    /// Validate a 2-D grid index against the grid dimensions.
    #[allow(dead_code)]
    pub fn check_index(x: usize, y: usize, nx: usize, ny: usize) -> Result<()> {
        if x < nx && y < ny {
            Ok(())
        } else {
            Err(Self::IndexOutOfBounds { x, y, nx, ny })
        }
    }

    /// Check that a 3-D grid size is valid (all dimensions > 0).
    #[allow(dead_code)]
    pub fn check_grid_3d(nx: usize, ny: usize, nz: usize) -> Result<()> {
        if nx == 0 || ny == 0 || nz == 0 {
            Err(Self::InvalidGridDimensions3D { nx, ny, nz })
        } else {
            Ok(())
        }
    }

    /// Check that a 3-D grid index is within bounds.
    #[allow(dead_code)]
    pub fn check_index_3d(
        x: usize,
        y: usize,
        z: usize,
        nx: usize,
        ny: usize,
        nz: usize,
    ) -> Result<()> {
        if x < nx && y < ny && z < nz {
            Ok(())
        } else {
            Err(Self::IndexOutOfBounds3D {
                x,
                y,
                z,
                nx,
                ny,
                nz,
            })
        }
    }

    /// Check that a distribution function value is non-negative.
    #[allow(dead_code)]
    pub fn check_distribution(value: f64, direction: usize, index: usize) -> Result<()> {
        if value >= 0.0 {
            Ok(())
        } else {
            Err(Self::NegativeDistribution {
                direction,
                index,
                value,
            })
        }
    }

    /// Check that a Knudsen number is within the continuum (BGK) range.
    ///
    /// The BGK approximation holds for Kn < 0.1.
    #[allow(dead_code)]
    pub fn check_knudsen(kn: f64) -> Result<()> {
        const KN_MAX: f64 = 0.1;
        if kn >= 0.0 && kn <= KN_MAX {
            Ok(())
        } else {
            Err(Self::KnudsenNumberOutOfRange { kn, kn_max: KN_MAX })
        }
    }

    /// Check that a time step does not exceed the CFL stability bound.
    #[allow(dead_code)]
    pub fn check_dt(dt: f64, dt_max: f64) -> Result<()> {
        if dt <= dt_max {
            Ok(())
        } else {
            Err(Self::TimeStepTooLarge { dt, dt_max })
        }
    }

    /// Check that a phase-field order parameter is within \[-1, 1\].
    #[allow(dead_code)]
    pub fn check_phase_field(phi: f64, index: usize) -> Result<()> {
        if phi >= -1.0 && phi <= 1.0 {
            Ok(())
        } else {
            Err(Self::PhaseFieldOutOfRange { phi, index })
        }
    }

    /// Check that the wall-normal distance in wall units is positive.
    #[allow(dead_code)]
    pub fn check_wall_distance(y_plus: f64, index: usize) -> Result<()> {
        if y_plus > 0.0 {
            Ok(())
        } else {
            Err(Self::NonPositiveWallDistance { y_plus, index })
        }
    }

    /// Check that the IBM marker count meets the minimum requirement.
    #[allow(dead_code)]
    pub fn check_ibm_markers(actual: usize, min_markers: usize) -> Result<()> {
        if actual >= min_markers {
            Ok(())
        } else {
            Err(Self::InvalidIbmMarkerCount {
                min_markers,
                actual,
            })
        }
    }

    /// Check that a Reynolds number does not exceed the grid-resolution limit.
    #[allow(dead_code)]
    pub fn check_reynolds(re: f64, re_max: f64, grid_res: usize) -> Result<()> {
        if re <= re_max {
            Ok(())
        } else {
            Err(Self::ReynoldsNumberTooHigh {
                re,
                re_max,
                grid_res,
            })
        }
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    // ── General ──

    #[test]
    fn test_general_error_message() {
        let e = Error::general("something went wrong");
        assert_eq!(e.to_string(), "something went wrong");
    }

    #[test]
    fn test_general_empty_message() {
        let e = Error::general("");
        assert_eq!(e.to_string(), "");
    }

    #[test]
    fn test_general_error_variant() {
        let e = Error::General("test".to_string());
        assert!(matches!(e, Error::General(_)));
    }

    // ── InvalidGridDimensions ──

    #[test]
    fn test_invalid_grid_zero_nx() {
        let e = Error::invalid_grid(0, 10);
        assert!(matches!(e, Error::InvalidGridDimensions { nx: 0, ny: 10 }));
    }

    #[test]
    fn test_invalid_grid_zero_ny() {
        let e = Error::invalid_grid(5, 0);
        assert!(matches!(e, Error::InvalidGridDimensions { nx: 5, ny: 0 }));
    }

    #[test]
    fn test_invalid_grid_message_contains_dims() {
        let e = Error::invalid_grid(3, 7);
        let msg = e.to_string();
        assert!(msg.contains("3"), "nx not in message: {msg}");
        assert!(msg.contains("7"), "ny not in message: {msg}");
    }

    #[test]
    fn test_check_grid_valid() {
        assert!(Error::check_grid(10, 20).is_ok());
    }

    #[test]
    fn test_check_grid_zero_nx_returns_err() {
        assert!(Error::check_grid(0, 10).is_err());
    }

    #[test]
    fn test_check_grid_zero_ny_returns_err() {
        assert!(Error::check_grid(10, 0).is_err());
    }

    #[test]
    fn test_check_grid_both_zero_returns_err() {
        assert!(Error::check_grid(0, 0).is_err());
    }

    // ── InvalidGridDimensions3D ──

    #[test]
    fn test_invalid_grid_3d_message() {
        let e = Error::InvalidGridDimensions3D {
            nx: 0,
            ny: 5,
            nz: 10,
        };
        let msg = e.to_string();
        assert!(msg.contains("0"), "nz in message: {msg}");
        assert!(msg.contains("3D"), "should mention 3D: {msg}");
    }

    // ── UnstableRelaxationTime ──

    #[test]
    fn test_unstable_tau_below_half() {
        let e = Error::unstable_tau(0.3);
        assert!(matches!(e, Error::UnstableRelaxationTime { tau } if (tau - 0.3).abs() < 1e-15));
    }

    #[test]
    fn test_unstable_tau_at_half() {
        let e = Error::unstable_tau(0.5);
        assert!(e.to_string().contains("0.5"));
    }

    #[test]
    fn test_check_tau_stable() {
        assert!(Error::check_tau(0.6).is_ok());
    }

    #[test]
    fn test_check_tau_exactly_half_unstable() {
        assert!(Error::check_tau(0.5).is_err());
    }

    #[test]
    fn test_check_tau_below_half_unstable() {
        assert!(Error::check_tau(0.1).is_err());
    }

    #[test]
    fn test_check_tau_large_stable() {
        assert!(Error::check_tau(2.0).is_ok());
    }

    // ── ParameterOutOfBounds ──

    #[test]
    fn test_parameter_out_of_bounds_message() {
        let e = Error::ParameterOutOfBounds {
            name: "porosity".to_string(),
            value: 1.5,
            min: 0.0,
            max: 1.0,
        };
        let msg = e.to_string();
        assert!(msg.contains("porosity"), "name in message: {msg}");
        assert!(msg.contains("1.5"), "value in message: {msg}");
    }

    #[test]
    fn test_check_range_valid() {
        assert!(Error::check_range(0.5, 0.0, 1.0, "porosity").is_ok());
    }

    #[test]
    fn test_check_range_at_min_valid() {
        assert!(Error::check_range(0.0, 0.0, 1.0, "porosity").is_ok());
    }

    #[test]
    fn test_check_range_at_max_valid() {
        assert!(Error::check_range(1.0, 0.0, 1.0, "porosity").is_ok());
    }

    #[test]
    fn test_check_range_below_min_invalid() {
        assert!(Error::check_range(-0.1, 0.0, 1.0, "porosity").is_err());
    }

    #[test]
    fn test_check_range_above_max_invalid() {
        assert!(Error::check_range(1.1, 0.0, 1.0, "porosity").is_err());
    }

    // ── NonFiniteField ──

    #[test]
    fn test_nan_detected() {
        assert!(Error::check_finite(f64::NAN, "density", 0).is_err());
    }

    #[test]
    fn test_inf_detected() {
        assert!(Error::check_finite(f64::INFINITY, "velocity_x", 5).is_err());
    }

    #[test]
    fn test_neg_inf_detected() {
        assert!(Error::check_finite(f64::NEG_INFINITY, "velocity_y", 3).is_err());
    }

    #[test]
    fn test_finite_value_ok() {
        assert!(Error::check_finite(1.0, "density", 0).is_ok());
    }

    #[test]
    fn test_non_finite_message_contains_field() {
        let e = Error::NonFiniteField {
            field: "pressure".to_string(),
            index: 42,
            value: f64::NAN,
        };
        let msg = e.to_string();
        assert!(msg.contains("pressure"), "field in message: {msg}");
        assert!(msg.contains("42"), "index in message: {msg}");
    }

    // ── BoundaryConditionError ──

    #[test]
    fn test_bc_error_message() {
        let e = Error::BoundaryConditionError {
            face: "inlet".to_string(),
            reason: "pressure not specified".to_string(),
        };
        let msg = e.to_string();
        assert!(msg.contains("inlet"), "face in message: {msg}");
        assert!(
            msg.contains("pressure not specified"),
            "reason in message: {msg}"
        );
    }

    // ── IndexOutOfBounds ──

    #[test]
    fn test_index_out_of_bounds_ok() {
        assert!(Error::check_index(2, 3, 10, 10).is_ok());
    }

    #[test]
    fn test_index_out_of_bounds_x_too_large() {
        assert!(Error::check_index(10, 0, 10, 10).is_err());
    }

    #[test]
    fn test_index_out_of_bounds_y_too_large() {
        assert!(Error::check_index(0, 10, 10, 10).is_err());
    }

    #[test]
    fn test_index_out_of_bounds_message() {
        let e = Error::IndexOutOfBounds {
            x: 15,
            y: 3,
            nx: 10,
            ny: 10,
        };
        let msg = e.to_string();
        assert!(msg.contains("15"), "x in message: {msg}");
        assert!(msg.contains("10"), "nx in message: {msg}");
    }

    // ── IndexOutOfBounds3D ──

    #[test]
    fn test_index_3d_message() {
        let e = Error::IndexOutOfBounds3D {
            x: 5,
            y: 6,
            z: 7,
            nx: 4,
            ny: 4,
            nz: 4,
        };
        let msg = e.to_string();
        assert!(msg.contains("5"), "x in message: {msg}");
        assert!(msg.contains("3D"), "should mention 3D: {msg}");
    }

    // ── ConvergenceFailure ──

    #[test]
    fn test_convergence_failure_message() {
        let e = Error::ConvergenceFailure {
            residual: 1e-3,
            iterations: 1000,
            tolerance: 1e-6,
        };
        let msg = e.to_string();
        assert!(msg.contains("1000"), "iterations in message: {msg}");
        assert!(
            msg.contains("convergence"),
            "should mention convergence: {msg}"
        );
    }

    // ── IncompatibleParameters ──

    #[test]
    fn test_incompatible_parameters_message() {
        let e = Error::IncompatibleParameters {
            message: "Prandtl number must be positive".to_string(),
        };
        let msg = e.to_string();
        assert!(msg.contains("Prandtl"), "reason in message: {msg}");
    }

    // ── InvalidSpeciesCount ──

    #[test]
    fn test_invalid_species_count_message() {
        let e = Error::InvalidSpeciesCount {
            min_species: 2,
            actual: 0,
        };
        let msg = e.to_string();
        assert!(msg.contains("2"), "min_species in message: {msg}");
        assert!(msg.contains("0"), "actual in message: {msg}");
    }

    // ── UnphysicalTemperature ──

    #[test]
    fn test_check_temperature_positive_ok() {
        assert!(Error::check_temperature(300.0, 0).is_ok());
    }

    #[test]
    fn test_check_temperature_zero_err() {
        assert!(Error::check_temperature(0.0, 0).is_err());
    }

    #[test]
    fn test_check_temperature_negative_err() {
        assert!(Error::check_temperature(-1.0, 5).is_err());
    }

    #[test]
    fn test_unphysical_temperature_message() {
        let e = Error::UnphysicalTemperature {
            temperature: -10.0,
            index: 7,
        };
        let msg = e.to_string();
        assert!(msg.contains("7"), "index in message: {msg}");
    }

    // ── MachNumberExceeded ──

    #[test]
    fn test_check_mach_safe() {
        assert!(Error::check_mach(0.1).is_ok());
    }

    #[test]
    fn test_check_mach_at_limit_err() {
        assert!(Error::check_mach(0.3).is_err());
    }

    #[test]
    fn test_check_mach_exceeded_err() {
        assert!(Error::check_mach(0.5).is_err());
    }

    #[test]
    fn test_mach_error_message() {
        let e = Error::MachNumberExceeded {
            mach: 0.45,
            limit: 0.3,
        };
        let msg = e.to_string();
        assert!(
            msg.contains("0.45") || msg.contains("45"),
            "mach in message: {msg}"
        );
        assert!(msg.contains("0.3"), "limit in message: {msg}");
    }

    // ── UnsupportedVelocitySet ──

    #[test]
    fn test_unsupported_velocity_set_message() {
        let e = Error::UnsupportedVelocitySet {
            name: "D4Q41".to_string(),
        };
        let msg = e.to_string();
        assert!(msg.contains("D4Q41"), "name in message: {msg}");
    }

    // ── Result type alias ──

    #[test]
    fn test_result_ok_carries_value() {
        let r: Result<i32> = Ok(42);
        assert!(matches!(r, Ok(42)));
    }

    #[test]
    fn test_result_err_carries_error() {
        let r: Result<i32> = Err(Error::general("oops"));
        assert!(r.is_err());
    }

    // ── Clone and PartialEq ──

    #[test]
    fn test_error_clone() {
        let e = Error::general("cloned");
        let e2 = e.clone();
        assert_eq!(e, e2);
    }

    #[test]
    fn test_errors_not_equal_different_variants() {
        let e1 = Error::general("a");
        let e2 = Error::UnstableRelaxationTime { tau: 0.3 };
        assert_ne!(e1, e2);
    }

    #[test]
    fn test_errors_equal_same_general() {
        let e1 = Error::General("same".to_string());
        let e2 = Error::General("same".to_string());
        assert_eq!(e1, e2);
    }

    // ── Debug formatting ──

    #[test]
    fn test_debug_format_non_empty() {
        let e = Error::general("debug test");
        let dbg = format!("{e:?}");
        assert!(!dbg.is_empty());
    }

    // ── Chained Result propagation ──

    #[test]
    fn test_question_mark_propagation() {
        fn inner() -> Result<f64> {
            Error::check_tau(0.3)?;
            Ok(1.0)
        }
        assert!(inner().is_err());
    }

    #[test]
    fn test_question_mark_no_error() {
        fn inner() -> Result<f64> {
            Error::check_tau(1.0)?;
            Ok(42.0)
        }
        assert_eq!(inner().unwrap(), 42.0);
    }

    // ── Multiple checks in sequence ──

    #[test]
    fn test_multiple_checks_all_pass() {
        let result: Result<()> = (|| {
            Error::check_grid(10, 10)?;
            Error::check_tau(0.8)?;
            Error::check_finite(1.0, "density", 0)?;
            Error::check_range(0.4, 0.0, 1.0, "porosity")?;
            Error::check_temperature(300.0, 0)?;
            Error::check_mach(0.1)?;
            Ok(())
        })();
        assert!(result.is_ok());
    }

    #[test]
    fn test_multiple_checks_first_fails() {
        let result: Result<()> = (|| {
            Error::check_grid(0, 10)?;
            Error::check_tau(0.8)?;
            Ok(())
        })();
        assert!(result.is_err());
        assert!(matches!(
            result.unwrap_err(),
            Error::InvalidGridDimensions { .. }
        ));
    }

    // ── check_grid_3d ──

    #[test]
    fn test_check_grid_3d_valid() {
        assert!(Error::check_grid_3d(4, 4, 4).is_ok());
    }

    #[test]
    fn test_check_grid_3d_zero_nx_err() {
        assert!(Error::check_grid_3d(0, 4, 4).is_err());
    }

    #[test]
    fn test_check_grid_3d_zero_ny_err() {
        assert!(Error::check_grid_3d(4, 0, 4).is_err());
    }

    #[test]
    fn test_check_grid_3d_zero_nz_err() {
        assert!(Error::check_grid_3d(4, 4, 0).is_err());
    }

    #[test]
    fn test_check_grid_3d_all_zero_err() {
        assert!(Error::check_grid_3d(0, 0, 0).is_err());
    }

    // ── check_index_3d ──

    #[test]
    fn test_check_index_3d_valid() {
        assert!(Error::check_index_3d(1, 2, 3, 4, 4, 4).is_ok());
    }

    #[test]
    fn test_check_index_3d_x_eq_nx_err() {
        assert!(Error::check_index_3d(4, 0, 0, 4, 4, 4).is_err());
    }

    #[test]
    fn test_check_index_3d_y_eq_ny_err() {
        assert!(Error::check_index_3d(0, 4, 0, 4, 4, 4).is_err());
    }

    #[test]
    fn test_check_index_3d_z_eq_nz_err() {
        assert!(Error::check_index_3d(0, 0, 4, 4, 4, 4).is_err());
    }

    #[test]
    fn test_check_index_3d_boundary_ok() {
        assert!(Error::check_index_3d(3, 3, 3, 4, 4, 4).is_ok());
    }

    // ── NegativeDistribution / check_distribution ──

    #[test]
    fn test_check_distribution_positive_ok() {
        assert!(Error::check_distribution(0.1, 1, 0).is_ok());
    }

    #[test]
    fn test_check_distribution_zero_ok() {
        assert!(Error::check_distribution(0.0, 0, 0).is_ok());
    }

    #[test]
    fn test_check_distribution_negative_err() {
        assert!(Error::check_distribution(-0.001, 3, 42).is_err());
    }

    #[test]
    fn test_negative_distribution_message() {
        let e = Error::NegativeDistribution {
            direction: 5,
            index: 99,
            value: -0.05,
        };
        let msg = e.to_string();
        assert!(msg.contains("5"), "direction in message: {msg}");
        assert!(msg.contains("99"), "index in message: {msg}");
    }

    // ── KnudsenNumberOutOfRange / check_knudsen ──

    #[test]
    fn test_check_knudsen_valid() {
        assert!(Error::check_knudsen(0.05).is_ok());
    }

    #[test]
    fn test_check_knudsen_zero_ok() {
        assert!(Error::check_knudsen(0.0).is_ok());
    }

    #[test]
    fn test_check_knudsen_at_limit_ok() {
        assert!(Error::check_knudsen(0.1).is_ok());
    }

    #[test]
    fn test_check_knudsen_above_limit_err() {
        assert!(Error::check_knudsen(0.2).is_err());
    }

    #[test]
    fn test_check_knudsen_negative_err() {
        assert!(Error::check_knudsen(-0.01).is_err());
    }

    #[test]
    fn test_knudsen_message_contains_value() {
        let e = Error::KnudsenNumberOutOfRange {
            kn: 0.5,
            kn_max: 0.1,
        };
        let msg = e.to_string();
        assert!(
            msg.contains("0.5") || msg.contains("Knudsen"),
            "kn in message: {msg}"
        );
    }

    // ── TimeStepTooLarge / check_dt ──

    #[test]
    fn test_check_dt_valid() {
        assert!(Error::check_dt(0.001, 0.01).is_ok());
    }

    #[test]
    fn test_check_dt_equal_max_ok() {
        assert!(Error::check_dt(0.01, 0.01).is_ok());
    }

    #[test]
    fn test_check_dt_too_large_err() {
        assert!(Error::check_dt(0.1, 0.01).is_err());
    }

    #[test]
    fn test_time_step_message() {
        let e = Error::TimeStepTooLarge {
            dt: 0.1,
            dt_max: 0.01,
        };
        let msg = e.to_string();
        assert!(
            msg.contains("0.1") || msg.contains("CFL"),
            "dt in message: {msg}"
        );
    }

    // ── PhaseFieldOutOfRange / check_phase_field ──

    #[test]
    fn test_check_phase_field_valid() {
        assert!(Error::check_phase_field(0.0, 0).is_ok());
    }

    #[test]
    fn test_check_phase_field_minus_one_ok() {
        assert!(Error::check_phase_field(-1.0, 0).is_ok());
    }

    #[test]
    fn test_check_phase_field_plus_one_ok() {
        assert!(Error::check_phase_field(1.0, 0).is_ok());
    }

    #[test]
    fn test_check_phase_field_above_one_err() {
        assert!(Error::check_phase_field(1.5, 0).is_err());
    }

    #[test]
    fn test_check_phase_field_below_minus_one_err() {
        assert!(Error::check_phase_field(-1.5, 5).is_err());
    }

    #[test]
    fn test_phase_field_message_contains_index() {
        let e = Error::PhaseFieldOutOfRange {
            phi: 2.0,
            index: 77,
        };
        let msg = e.to_string();
        assert!(msg.contains("77"), "index in message: {msg}");
    }

    // ── NonPositiveWallDistance / check_wall_distance ──

    #[test]
    fn test_check_wall_distance_positive_ok() {
        assert!(Error::check_wall_distance(5.0, 0).is_ok());
    }

    #[test]
    fn test_check_wall_distance_zero_err() {
        assert!(Error::check_wall_distance(0.0, 0).is_err());
    }

    #[test]
    fn test_check_wall_distance_negative_err() {
        assert!(Error::check_wall_distance(-1.0, 3).is_err());
    }

    #[test]
    fn test_wall_distance_message() {
        let e = Error::NonPositiveWallDistance {
            y_plus: 0.0,
            index: 12,
        };
        let msg = e.to_string();
        assert!(msg.contains("12"), "index in message: {msg}");
    }

    // ── InvalidIbmMarkerCount / check_ibm_markers ──

    #[test]
    fn test_check_ibm_markers_valid() {
        assert!(Error::check_ibm_markers(10, 4).is_ok());
    }

    #[test]
    fn test_check_ibm_markers_equal_min_ok() {
        assert!(Error::check_ibm_markers(4, 4).is_ok());
    }

    #[test]
    fn test_check_ibm_markers_below_min_err() {
        assert!(Error::check_ibm_markers(2, 4).is_err());
    }

    #[test]
    fn test_ibm_marker_message() {
        let e = Error::InvalidIbmMarkerCount {
            min_markers: 8,
            actual: 2,
        };
        let msg = e.to_string();
        assert!(msg.contains("8"), "min_markers in message: {msg}");
        assert!(msg.contains("2"), "actual in message: {msg}");
    }

    // ── ReynoldsNumberTooHigh / check_reynolds ──

    #[test]
    fn test_check_reynolds_valid() {
        assert!(Error::check_reynolds(100.0, 1000.0, 64).is_ok());
    }

    #[test]
    fn test_check_reynolds_equal_max_ok() {
        assert!(Error::check_reynolds(1000.0, 1000.0, 64).is_ok());
    }

    #[test]
    fn test_check_reynolds_exceeded_err() {
        assert!(Error::check_reynolds(2000.0, 1000.0, 64).is_err());
    }

    #[test]
    fn test_reynolds_message_contains_value() {
        let e = Error::ReynoldsNumberTooHigh {
            re: 5000.0,
            re_max: 1000.0,
            grid_res: 32,
        };
        let msg = e.to_string();
        assert!(
            msg.contains("5000") || msg.contains("Reynolds"),
            "re in message: {msg}"
        );
    }

    // ── ComponentCountMismatch ──

    #[test]
    fn test_component_count_mismatch_message() {
        let e = Error::ComponentCountMismatch {
            expected: 3,
            found: 2,
            field: "concentration".to_string(),
        };
        let msg = e.to_string();
        assert!(msg.contains("3"), "expected in message: {msg}");
        assert!(msg.contains("2"), "found in message: {msg}");
        assert!(msg.contains("concentration"), "field in message: {msg}");
    }

    // ── InvalidTurbulenceParameter ──

    #[test]
    fn test_invalid_turbulence_parameter_message() {
        let e = Error::InvalidTurbulenceParameter {
            name: "Cs".to_string(),
            reason: "must be in (0, 0.5)".to_string(),
        };
        let msg = e.to_string();
        assert!(msg.contains("Cs"), "name in message: {msg}");
        assert!(msg.contains("0.5"), "reason in message: {msg}");
    }

    // ── InvalidParticleConfiguration ──

    #[test]
    fn test_invalid_particle_configuration_message() {
        let e = Error::InvalidParticleConfiguration {
            reason: "particle radius must be positive".to_string(),
        };
        let msg = e.to_string();
        assert!(msg.contains("radius"), "reason in message: {msg}");
    }

    // ── New variant clone/equality ──

    #[test]
    fn test_new_error_variants_clone() {
        let e = Error::ReynoldsNumberTooHigh {
            re: 2000.0,
            re_max: 1000.0,
            grid_res: 32,
        };
        let e2 = e.clone();
        assert_eq!(e, e2);
    }

    #[test]
    fn test_knudsen_error_not_equal_general() {
        let e1 = Error::KnudsenNumberOutOfRange {
            kn: 0.5,
            kn_max: 0.1,
        };
        let e2 = Error::general("different");
        assert_ne!(e1, e2);
    }

    // ── Combined validation pipeline ──

    #[test]
    fn test_combined_validation_pipeline_pass() {
        let result: Result<()> = (|| {
            Error::check_grid_3d(8, 8, 8)?;
            Error::check_tau(1.0)?;
            Error::check_knudsen(0.05)?;
            Error::check_dt(0.001, 0.01)?;
            Error::check_phase_field(0.5, 0)?;
            Error::check_wall_distance(10.0, 0)?;
            Error::check_ibm_markers(16, 4)?;
            Error::check_reynolds(500.0, 1000.0, 64)?;
            Ok(())
        })();
        assert!(result.is_ok());
    }

    #[test]
    fn test_combined_validation_pipeline_knudsen_fail() {
        let result: Result<()> = (|| {
            Error::check_grid_3d(8, 8, 8)?;
            Error::check_knudsen(1.0)?; // fails
            Ok(())
        })();
        assert!(result.is_err());
        assert!(matches!(
            result.unwrap_err(),
            Error::KnudsenNumberOutOfRange { .. }
        ));
    }

    #[test]
    fn test_combined_validation_phase_field_fail() {
        let result: Result<()> = (|| {
            Error::check_phase_field(2.0, 5)?; // fails
            Ok(())
        })();
        assert!(matches!(
            result.unwrap_err(),
            Error::PhaseFieldOutOfRange { .. }
        ));
    }
}
