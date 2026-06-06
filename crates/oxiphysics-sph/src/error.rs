// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Error types for oxiphysics-sph.
//!
//! This module provides the unified [`enum@Error`] enum and [`Result`] alias used
//! throughout the crate.  The variants cover all major failure modes:
//!
//! - Configuration errors (invalid parameters)
//! - Convergence failures in iterative solvers
//! - Index out-of-range errors
//! - NaN / non-finite quantity detection
//! - Boundary and domain errors
//! - Open-boundary condition errors
//! - Turbulence model errors
//! - Coupling errors (SPH–rigid body, FSI, SPH–FEM)
//! - Memory and allocation errors
//! - I/O errors (particle data serialization)

use thiserror::Error;

// ---------------------------------------------------------------------------
// Primary error enum
// ---------------------------------------------------------------------------

/// Unified error type for the SPH crate.
#[derive(Debug, Error)]
pub enum Error {
    // ---- generic -----------------------------------------------------------
    /// Generic catch-all error with a descriptive message.
    #[error("{0}")]
    General(String),

    // ---- iterative solvers -------------------------------------------------
    /// An iterative solver (DFSPH, IISPH, PCISPH, …) did not converge within
    /// the allowed number of iterations.
    ///
    /// Contains the solver name, iteration count reached, and final residual.
    #[error(
        "solver '{solver}' did not converge after {iterations} iterations \
         (residual = {residual:.3e}, tolerance = {tolerance:.3e})"
    )]
    SolverNotConverged {
        /// Name of the solver that failed.
        solver: String,
        /// Number of iterations performed.
        iterations: usize,
        /// Residual at the last iteration.
        residual: f64,
        /// Convergence tolerance that was not met.
        tolerance: f64,
    },

    // ---- parameter validation ----------------------------------------------
    /// A parameter value is outside its valid range.
    ///
    /// Contains the parameter name, the supplied value, and the acceptable range.
    #[error("parameter '{name}' = {value:.6} is out of range [{min:.6}, {max:.6}]")]
    InvalidParameter {
        /// Parameter name.
        name: String,
        /// Supplied value.
        value: f64,
        /// Lower bound of the valid range.
        min: f64,
        /// Upper bound of the valid range.
        max: f64,
    },

    // ---- indexing ----------------------------------------------------------
    /// A particle index is out of range for the current particle set.
    #[error("particle index {index} is out of range (particle count = {count})")]
    IndexOutOfRange {
        /// Requested index.
        index: usize,
        /// Current particle count.
        count: usize,
    },

    // ---- numerics ----------------------------------------------------------
    /// A quantity computed during the simulation became NaN or infinite.
    ///
    /// Contains the name of the quantity (e.g., "pressure", "density") and the
    /// particle index where the problem was detected.
    #[error(
        "quantity '{quantity}' for particle {particle_index} is non-finite \
         (value = {value:?})"
    )]
    NonFiniteQuantity {
        /// Name of the physical quantity.
        quantity: String,
        /// Particle index.
        particle_index: usize,
        /// The non-finite value.
        value: f64,
    },

    // ---- domain / boundary -------------------------------------------------
    /// The simulation domain boundary was violated.
    ///
    /// Particles that escape the simulation box may produce incorrect results.
    #[error(
        "particle {particle_index} at position {position:?} is outside the \
         domain [{lo:?}, {hi:?}]"
    )]
    DomainViolation {
        /// Particle index.
        particle_index: usize,
        /// Current position of the particle.
        position: [f64; 3],
        /// Lower bound of the domain.
        lo: [f64; 3],
        /// Upper bound of the domain.
        hi: [f64; 3],
    },

    // ---- multiphase --------------------------------------------------------
    /// A requested phase does not exist in the multiphase system.
    #[error("phase index {phase_id} does not exist (registered phases: {n_phases})")]
    PhaseNotFound {
        /// Requested phase ID.
        phase_id: usize,
        /// Number of registered phases.
        n_phases: usize,
    },

    // ---- time stepping -----------------------------------------------------
    /// The time-step selected by the CFL condition was below the minimum
    /// allowed value, indicating a stability problem.
    #[error(
        "CFL time-step {dt:.3e} is below minimum {dt_min:.3e} \
         (max speed = {max_speed:.3e} m/s)"
    )]
    TimestepTooSmall {
        /// Computed CFL time-step.
        dt: f64,
        /// Minimum allowed time-step.
        dt_min: f64,
        /// Maximum particle speed at the time of failure.
        max_speed: f64,
    },

    /// The time-step exceeds the maximum allowed value.
    #[error("time-step {dt:.3e} exceeds maximum allowed {dt_max:.3e}")]
    TimestepTooLarge {
        /// Proposed time-step.
        dt: f64,
        /// Maximum allowed time-step.
        dt_max: f64,
    },

    // ---- particle counts ---------------------------------------------------
    /// Insufficient particles to build a valid neighbor list or kernel support.
    #[error("too few particles in the simulation: {count} (minimum required: {minimum})")]
    TooFewParticles {
        /// Current particle count.
        count: usize,
        /// Minimum required.
        minimum: usize,
    },

    /// Particle count exceeds the allocated capacity.
    #[error("particle count {count} exceeds allocated capacity {capacity}")]
    CapacityExceeded {
        /// Current particle count.
        count: usize,
        /// Maximum capacity.
        capacity: usize,
    },

    // ---- open boundary -----------------------------------------------------
    /// The inflow rate at an open boundary exceeded the allowed maximum.
    #[error(
        "open boundary '{name}': inflow rate {rate:.3e} m³/s exceeds maximum {max_rate:.3e} m³/s"
    )]
    InflowRateExceeded {
        /// Boundary name or identifier.
        name: String,
        /// Actual inflow rate.
        rate: f64,
        /// Maximum permissible inflow rate.
        max_rate: f64,
    },

    /// The buffer zone for an open boundary is undersized.
    #[error(
        "open boundary '{name}': buffer zone depth {depth:.4} m is less than \
         the minimum required {min_depth:.4} m (≥ {n_layers} kernel supports)"
    )]
    BufferZoneTooSmall {
        /// Boundary name or identifier.
        name: String,
        /// Current buffer depth.
        depth: f64,
        /// Minimum required depth.
        min_depth: f64,
        /// Minimum number of kernel support layers required.
        n_layers: usize,
    },

    /// An open boundary has an ill-defined normal vector (zero magnitude).
    #[error(
        "open boundary '{name}': boundary normal {normal:?} has near-zero magnitude ({mag:.2e})"
    )]
    InvalidBoundaryNormal {
        /// Boundary name or identifier.
        name: String,
        /// Provided normal vector.
        normal: [f64; 3],
        /// Magnitude of the provided normal.
        mag: f64,
    },

    // ---- turbulence models -------------------------------------------------
    /// The turbulence model configuration is invalid.
    #[error("turbulence model '{model}': {reason}")]
    TurbulenceModelError {
        /// Model name (e.g., "Smagorinsky", "k-omega", "DES").
        model: String,
        /// Human-readable reason.
        reason: String,
    },

    /// A turbulent kinetic energy became negative (unphysical).
    #[error(
        "turbulent kinetic energy k = {k:.3e} < 0 for particle {particle_index} \
         in model '{model}'"
    )]
    NegativeTurbulentKineticEnergy {
        /// Particle index.
        particle_index: usize,
        /// Model name.
        model: String,
        /// Non-physical TKE value.
        k: f64,
    },

    /// The turbulent dissipation rate became non-positive (unphysical).
    #[error(
        "turbulent dissipation ε/ω = {dissipation:.3e} ≤ 0 for particle {particle_index} \
         in model '{model}'"
    )]
    NonPositiveDissipation {
        /// Particle index.
        particle_index: usize,
        /// Model name.
        model: String,
        /// Non-physical dissipation value.
        dissipation: f64,
    },

    // ---- coupling (SPH–rigid body / FSI / SPH–FEM) -------------------------
    /// The coupled rigid body or FEM mesh has no elements.
    #[error("coupling target '{target}' has no {element_type} (count = 0)")]
    EmptyCouplingTarget {
        /// Name of the coupling target (body name, mesh name, …).
        target: String,
        /// Element type description ("vertices", "faces", "elements", …).
        element_type: String,
    },

    /// The interface transfer between SPH and FEM failed to converge.
    #[error(
        "SPH–FEM interface transfer '{interface}' failed to converge after {iterations} \
         iterations (residual = {residual:.3e})"
    )]
    InterfaceTransferNotConverged {
        /// Interface name.
        interface: String,
        /// Iteration count.
        iterations: usize,
        /// Final residual.
        residual: f64,
    },

    /// Volume overlap between an SPH particle and an embedded structure is invalid.
    #[error(
        "volume coupling: particle {particle_index} has non-positive effective volume \
         {volume:.3e} m³"
    )]
    InvalidCouplingVolume {
        /// Particle index.
        particle_index: usize,
        /// Computed volume.
        volume: f64,
    },

    /// A rigid-body coupling force became non-finite.
    #[error(
        "rigid-body coupling: force on body '{body}' is non-finite \
         (force = {force:?} N)"
    )]
    NonFiniteCouplingForce {
        /// Body identifier.
        body: String,
        /// Computed force vector.
        force: [f64; 3],
    },

    // ---- kernel ------------------------------------------------------------
    /// The kernel smoothing length is non-positive.
    #[error("kernel smoothing length h = {h:.3e} m is non-positive for particle {particle_index}")]
    NonPositiveSmoothingLength {
        /// Particle index.
        particle_index: usize,
        /// Non-positive smoothing length.
        h: f64,
    },

    /// The kernel normalization constant is zero (degenerate configuration).
    #[error("kernel normalization constant is zero for kernel '{kernel_name}'")]
    ZeroKernelNormalization {
        /// Kernel name.
        kernel_name: String,
    },

    // ---- neighbor search ---------------------------------------------------
    /// The spatial hash grid could not be built because the cell size is invalid.
    #[error("neighbor search: cell size {cell_size:.3e} m is non-positive")]
    InvalidCellSize {
        /// Requested cell size.
        cell_size: f64,
    },

    /// The neighbor list for a particle is abnormally large, possibly indicating
    /// a bug in the neighbor search or particle placement.
    #[error(
        "neighbor search: particle {particle_index} has {n_neighbors} neighbors, \
         which exceeds the warning threshold {threshold}"
    )]
    TooManyNeighbors {
        /// Particle index.
        particle_index: usize,
        /// Number of found neighbors.
        n_neighbors: usize,
        /// Warning threshold.
        threshold: usize,
    },

    // ---- density -----------------------------------------------------------
    /// Particle density became non-positive (unphysical).
    #[error("density of particle {particle_index} is non-positive: {density:.3e} kg/m³")]
    NonPositiveDensity {
        /// Particle index.
        particle_index: usize,
        /// Non-positive density value.
        density: f64,
    },

    // ---- free surface ------------------------------------------------------
    /// Surface reconstruction failed for a given particle.
    #[error("surface reconstruction failed for particle {particle_index}: {reason}")]
    SurfaceReconstructionFailed {
        /// Particle index.
        particle_index: usize,
        /// Human-readable reason.
        reason: String,
    },

    // ---- I/O ---------------------------------------------------------------
    /// Failed to serialize or deserialize particle data.
    #[error("particle data I/O error: {reason}")]
    IoError {
        /// Description of the I/O failure.
        reason: String,
    },
}

/// Result type alias using [`enum@Error`].
pub type Result<T> = std::result::Result<T, Error>;

// ---------------------------------------------------------------------------
// Helper constructors
// ---------------------------------------------------------------------------

impl Error {
    /// Create a [`Error::General`] from any displayable value.
    pub fn general(msg: impl std::fmt::Display) -> Self {
        Error::General(msg.to_string())
    }

    /// Create a [`Error::SolverNotConverged`] error.
    pub fn not_converged(
        solver: impl Into<String>,
        iterations: usize,
        residual: f64,
        tolerance: f64,
    ) -> Self {
        Error::SolverNotConverged {
            solver: solver.into(),
            iterations,
            residual,
            tolerance,
        }
    }

    /// Create an [`Error::InvalidParameter`] error.
    pub fn invalid_param(name: impl Into<String>, value: f64, min: f64, max: f64) -> Self {
        Error::InvalidParameter {
            name: name.into(),
            value,
            min,
            max,
        }
    }

    /// Create an [`Error::IndexOutOfRange`] error.
    pub fn index_out_of_range(index: usize, count: usize) -> Self {
        Error::IndexOutOfRange { index, count }
    }

    /// Create an [`Error::NonFiniteQuantity`] error.
    pub fn non_finite(quantity: impl Into<String>, particle_index: usize, value: f64) -> Self {
        Error::NonFiniteQuantity {
            quantity: quantity.into(),
            particle_index,
            value,
        }
    }

    /// Check whether `value` is finite; if not, return a [`Error::NonFiniteQuantity`].
    pub fn check_finite(
        quantity: impl Into<String>,
        particle_index: usize,
        value: f64,
    ) -> Result<f64> {
        if value.is_finite() {
            Ok(value)
        } else {
            Err(Error::non_finite(quantity, particle_index, value))
        }
    }

    /// Check whether `index < count`; if not, return [`Error::IndexOutOfRange`].
    pub fn check_index(index: usize, count: usize) -> Result<()> {
        if index < count {
            Ok(())
        } else {
            Err(Error::index_out_of_range(index, count))
        }
    }

    /// Returns `true` if this is a convergence failure.
    pub fn is_convergence_failure(&self) -> bool {
        matches!(
            self,
            Error::SolverNotConverged { .. } | Error::InterfaceTransferNotConverged { .. }
        )
    }

    /// Returns `true` if this error is recoverable (e.g., non-finite can be
    /// reset, convergence can retry with smaller dt).
    pub fn is_recoverable(&self) -> bool {
        matches!(
            self,
            Error::SolverNotConverged { .. }
                | Error::NonFiniteQuantity { .. }
                | Error::TimestepTooSmall { .. }
                | Error::InterfaceTransferNotConverged { .. }
                | Error::NegativeTurbulentKineticEnergy { .. }
                | Error::NonPositiveDissipation { .. }
        )
    }

    /// Returns `true` if this error is related to a turbulence model.
    pub fn is_turbulence_error(&self) -> bool {
        matches!(
            self,
            Error::TurbulenceModelError { .. }
                | Error::NegativeTurbulentKineticEnergy { .. }
                | Error::NonPositiveDissipation { .. }
        )
    }

    /// Returns `true` if this error is related to open-boundary conditions.
    pub fn is_open_boundary_error(&self) -> bool {
        matches!(
            self,
            Error::InflowRateExceeded { .. }
                | Error::BufferZoneTooSmall { .. }
                | Error::InvalidBoundaryNormal { .. }
        )
    }

    /// Returns `true` if this error is related to coupling (SPH–rigid / FSI / SPH–FEM).
    pub fn is_coupling_error(&self) -> bool {
        matches!(
            self,
            Error::EmptyCouplingTarget { .. }
                | Error::InterfaceTransferNotConverged { .. }
                | Error::InvalidCouplingVolume { .. }
                | Error::NonFiniteCouplingForce { .. }
        )
    }

    /// Create a [`Error::TurbulenceModelError`].
    pub fn turbulence_model_error(model: impl Into<String>, reason: impl Into<String>) -> Self {
        Error::TurbulenceModelError {
            model: model.into(),
            reason: reason.into(),
        }
    }

    /// Create a [`Error::NegativeTurbulentKineticEnergy`].
    pub fn negative_tke(particle_index: usize, model: impl Into<String>, k: f64) -> Self {
        Error::NegativeTurbulentKineticEnergy {
            particle_index,
            model: model.into(),
            k,
        }
    }

    /// Create a [`Error::NonPositiveDissipation`].
    pub fn non_positive_dissipation(
        particle_index: usize,
        model: impl Into<String>,
        dissipation: f64,
    ) -> Self {
        Error::NonPositiveDissipation {
            particle_index,
            model: model.into(),
            dissipation,
        }
    }

    /// Create a [`Error::InflowRateExceeded`].
    pub fn inflow_rate_exceeded(name: impl Into<String>, rate: f64, max_rate: f64) -> Self {
        Error::InflowRateExceeded {
            name: name.into(),
            rate,
            max_rate,
        }
    }

    /// Create a [`Error::BufferZoneTooSmall`].
    pub fn buffer_zone_too_small(
        name: impl Into<String>,
        depth: f64,
        min_depth: f64,
        n_layers: usize,
    ) -> Self {
        Error::BufferZoneTooSmall {
            name: name.into(),
            depth,
            min_depth,
            n_layers,
        }
    }

    /// Create a [`Error::InvalidBoundaryNormal`].
    pub fn invalid_boundary_normal(name: impl Into<String>, normal: [f64; 3], mag: f64) -> Self {
        Error::InvalidBoundaryNormal {
            name: name.into(),
            normal,
            mag,
        }
    }

    /// Create a [`Error::EmptyCouplingTarget`].
    pub fn empty_coupling_target(
        target: impl Into<String>,
        element_type: impl Into<String>,
    ) -> Self {
        Error::EmptyCouplingTarget {
            target: target.into(),
            element_type: element_type.into(),
        }
    }

    /// Create a [`Error::InterfaceTransferNotConverged`].
    pub fn interface_transfer_not_converged(
        interface: impl Into<String>,
        iterations: usize,
        residual: f64,
    ) -> Self {
        Error::InterfaceTransferNotConverged {
            interface: interface.into(),
            iterations,
            residual,
        }
    }

    /// Create a [`Error::InvalidCouplingVolume`].
    pub fn invalid_coupling_volume(particle_index: usize, volume: f64) -> Self {
        Error::InvalidCouplingVolume {
            particle_index,
            volume,
        }
    }

    /// Create a [`Error::NonFiniteCouplingForce`].
    pub fn non_finite_coupling_force(body: impl Into<String>, force: [f64; 3]) -> Self {
        Error::NonFiniteCouplingForce {
            body: body.into(),
            force,
        }
    }

    /// Create a [`Error::NonPositiveSmoothingLength`].
    pub fn non_positive_smoothing_length(particle_index: usize, h: f64) -> Self {
        Error::NonPositiveSmoothingLength { particle_index, h }
    }

    /// Create a [`Error::NonPositiveDensity`].
    pub fn non_positive_density(particle_index: usize, density: f64) -> Self {
        Error::NonPositiveDensity {
            particle_index,
            density,
        }
    }

    /// Create a [`Error::SurfaceReconstructionFailed`].
    pub fn surface_reconstruction_failed(particle_index: usize, reason: impl Into<String>) -> Self {
        Error::SurfaceReconstructionFailed {
            particle_index,
            reason: reason.into(),
        }
    }

    /// Create a [`Error::IoError`].
    pub fn io_error(reason: impl Into<String>) -> Self {
        Error::IoError {
            reason: reason.into(),
        }
    }

    /// Validate that a smoothing length `h` is positive; return it or an error.
    pub fn check_smoothing_length(particle_index: usize, h: f64) -> Result<f64> {
        if h > 0.0 {
            Ok(h)
        } else {
            Err(Error::non_positive_smoothing_length(particle_index, h))
        }
    }

    /// Validate that a density `rho` is positive; return it or an error.
    pub fn check_density(particle_index: usize, rho: f64) -> Result<f64> {
        if rho > 0.0 {
            Ok(rho)
        } else {
            Err(Error::non_positive_density(particle_index, rho))
        }
    }

    /// Validate that a turbulent kinetic energy `k` is non-negative.
    pub fn check_tke(particle_index: usize, model: &str, k: f64) -> Result<f64> {
        if k >= 0.0 {
            Ok(k)
        } else {
            Err(Error::negative_tke(particle_index, model, k))
        }
    }

    /// Validate that a dissipation rate (ε or ω) is positive.
    pub fn check_dissipation(particle_index: usize, model: &str, dissipation: f64) -> Result<f64> {
        if dissipation > 0.0 {
            Ok(dissipation)
        } else {
            Err(Error::non_positive_dissipation(
                particle_index,
                model,
                dissipation,
            ))
        }
    }

    /// Validate a boundary normal vector (must have non-zero magnitude ≥ `eps`).
    pub fn check_boundary_normal(name: &str, normal: [f64; 3], eps: f64) -> Result<[f64; 3]> {
        let mag = (normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]).sqrt();
        if mag >= eps {
            Ok(normal)
        } else {
            Err(Error::invalid_boundary_normal(name, normal, mag))
        }
    }

    /// Returns a short category string for this error, useful for logging/metrics.
    pub fn category(&self) -> &'static str {
        match self {
            Error::General(_) => "general",
            Error::SolverNotConverged { .. } => "solver_convergence",
            Error::InvalidParameter { .. } => "invalid_parameter",
            Error::IndexOutOfRange { .. } => "index_out_of_range",
            Error::NonFiniteQuantity { .. } => "non_finite",
            Error::DomainViolation { .. } => "domain_violation",
            Error::PhaseNotFound { .. } => "phase_not_found",
            Error::TimestepTooSmall { .. } => "timestep_too_small",
            Error::TimestepTooLarge { .. } => "timestep_too_large",
            Error::TooFewParticles { .. } => "too_few_particles",
            Error::CapacityExceeded { .. } => "capacity_exceeded",
            Error::InflowRateExceeded { .. } => "open_boundary",
            Error::BufferZoneTooSmall { .. } => "open_boundary",
            Error::InvalidBoundaryNormal { .. } => "open_boundary",
            Error::TurbulenceModelError { .. } => "turbulence",
            Error::NegativeTurbulentKineticEnergy { .. } => "turbulence",
            Error::NonPositiveDissipation { .. } => "turbulence",
            Error::EmptyCouplingTarget { .. } => "coupling",
            Error::InterfaceTransferNotConverged { .. } => "coupling",
            Error::InvalidCouplingVolume { .. } => "coupling",
            Error::NonFiniteCouplingForce { .. } => "coupling",
            Error::NonPositiveSmoothingLength { .. } => "kernel",
            Error::ZeroKernelNormalization { .. } => "kernel",
            Error::InvalidCellSize { .. } => "neighbor_search",
            Error::TooManyNeighbors { .. } => "neighbor_search",
            Error::NonPositiveDensity { .. } => "density",
            Error::SurfaceReconstructionFailed { .. } => "surface",
            Error::IoError { .. } => "io",
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ---- General -----------------------------------------------------------

    #[test]
    fn general_error_display() {
        let e = Error::General("something went wrong".into());
        assert!(e.to_string().contains("something went wrong"));
    }

    #[test]
    fn general_constructor() {
        let e = Error::general("from constructor");
        assert!(e.to_string().contains("from constructor"));
        assert_eq!(e.category(), "general");
    }

    // ---- SolverNotConverged ------------------------------------------------

    #[test]
    fn solver_not_converged_display() {
        let e = Error::not_converged("DFSPH", 100, 1.5e-3, 1e-4);
        let s = e.to_string();
        assert!(s.contains("DFSPH"));
        assert!(s.contains("100"));
        assert!(e.is_convergence_failure());
    }

    #[test]
    fn solver_not_converged_recoverable() {
        let e = Error::not_converged("PCISPH", 50, 0.1, 1e-3);
        assert!(e.is_recoverable());
        assert_eq!(e.category(), "solver_convergence");
    }

    #[test]
    fn solver_not_converged_contains_residual() {
        let e = Error::not_converged("IISPH", 200, 9.99e-4, 1e-4);
        assert!(e.to_string().contains("IISPH"));
    }

    // ---- InvalidParameter --------------------------------------------------

    #[test]
    fn invalid_parameter_display() {
        let e = Error::invalid_param("smoothing_length", -0.1, 0.0, 1.0);
        assert!(e.to_string().contains("smoothing_length"));
        assert!(!e.is_convergence_failure());
        assert_eq!(e.category(), "invalid_parameter");
    }

    #[test]
    fn invalid_parameter_not_recoverable() {
        let e = Error::invalid_param("dt", -1.0, 0.0, 1.0);
        assert!(!e.is_recoverable());
    }

    // ---- IndexOutOfRange ---------------------------------------------------

    #[test]
    fn index_out_of_range_display() {
        let e = Error::index_out_of_range(99, 10);
        let s = e.to_string();
        assert!(s.contains("99"));
        assert!(s.contains("10"));
        assert_eq!(e.category(), "index_out_of_range");
    }

    #[test]
    fn check_index_valid() {
        assert!(Error::check_index(0, 5).is_ok());
        assert!(Error::check_index(4, 5).is_ok());
    }

    #[test]
    fn check_index_out_of_range() {
        let r = Error::check_index(5, 5);
        assert!(r.is_err());
        matches!(
            r.unwrap_err(),
            Error::IndexOutOfRange { index: 5, count: 5 }
        );
    }

    #[test]
    fn check_index_zero_count() {
        assert!(Error::check_index(0, 0).is_err());
    }

    // ---- NonFiniteQuantity -------------------------------------------------

    #[test]
    fn non_finite_display() {
        let e = Error::non_finite("pressure", 42, f64::NAN);
        let s = e.to_string();
        assert!(s.contains("pressure"));
        assert!(s.contains("42"));
        assert_eq!(e.category(), "non_finite");
    }

    #[test]
    fn check_finite_ok() {
        let r = Error::check_finite("density", 0, 998.0);
        assert!(r.is_ok());
        assert_eq!(r.unwrap(), 998.0);
    }

    #[test]
    fn check_finite_nan_returns_err() {
        let r = Error::check_finite("density", 7, f64::NAN);
        assert!(r.is_err());
        let e = r.unwrap_err();
        assert!(!e.is_convergence_failure());
        assert!(e.is_recoverable());
    }

    #[test]
    fn check_finite_inf_returns_err() {
        let r = Error::check_finite("velocity", 3, f64::INFINITY);
        assert!(r.is_err());
    }

    #[test]
    fn check_finite_neg_inf_returns_err() {
        let r = Error::check_finite("velocity", 3, f64::NEG_INFINITY);
        assert!(r.is_err());
    }

    // ---- DomainViolation ---------------------------------------------------

    #[test]
    fn domain_violation_display() {
        let e = Error::DomainViolation {
            particle_index: 3,
            position: [10.0, 0.0, 0.0],
            lo: [0.0; 3],
            hi: [5.0; 3],
        };
        let s = e.to_string();
        assert!(s.contains("3"));
        assert!(s.contains("10"));
        assert_eq!(e.category(), "domain_violation");
    }

    // ---- PhaseNotFound -----------------------------------------------------

    #[test]
    fn phase_not_found_display() {
        let e = Error::PhaseNotFound {
            phase_id: 5,
            n_phases: 2,
        };
        let s = e.to_string();
        assert!(s.contains("5"));
        assert!(s.contains("2"));
        assert_eq!(e.category(), "phase_not_found");
    }

    // ---- TimestepTooSmall --------------------------------------------------

    #[test]
    fn timestep_too_small_display() {
        let e = Error::TimestepTooSmall {
            dt: 1e-12,
            dt_min: 1e-8,
            max_speed: 500.0,
        };
        assert!(e.to_string().contains("1e-12") || e.to_string().contains("1.000e-12"));
        assert!(e.is_recoverable());
        assert_eq!(e.category(), "timestep_too_small");
    }

    // ---- TimestepTooLarge --------------------------------------------------

    #[test]
    fn timestep_too_large_display() {
        let e = Error::TimestepTooLarge {
            dt: 1.0,
            dt_max: 0.01,
        };
        assert!(
            e.to_string().contains("1.000e0")
                || e.to_string().contains("1.000e0")
                || e.to_string().contains("1")
        );
        assert!(!e.is_recoverable());
        assert_eq!(e.category(), "timestep_too_large");
    }

    // ---- TooFewParticles ---------------------------------------------------

    #[test]
    fn too_few_particles_display() {
        let e = Error::TooFewParticles {
            count: 2,
            minimum: 10,
        };
        let s = e.to_string();
        assert!(s.contains('2'));
        assert!(s.contains("10"));
        assert_eq!(e.category(), "too_few_particles");
    }

    // ---- CapacityExceeded --------------------------------------------------

    #[test]
    fn capacity_exceeded_display() {
        let e = Error::CapacityExceeded {
            count: 1001,
            capacity: 1000,
        };
        let s = e.to_string();
        assert!(s.contains("1001"));
        assert!(s.contains("1000"));
        assert_eq!(e.category(), "capacity_exceeded");
    }

    // ---- Open boundary errors ----------------------------------------------

    #[test]
    fn inflow_rate_exceeded_display() {
        let e = Error::inflow_rate_exceeded("inlet_left", 2.5, 1.0);
        let s = e.to_string();
        assert!(s.contains("inlet_left"));
        assert!(e.is_open_boundary_error());
        assert_eq!(e.category(), "open_boundary");
    }

    #[test]
    fn buffer_zone_too_small_display() {
        let e = Error::buffer_zone_too_small("outlet_right", 0.01, 0.05, 3);
        let s = e.to_string();
        assert!(s.contains("outlet_right"));
        assert!(s.contains("3"));
        assert!(e.is_open_boundary_error());
    }

    #[test]
    fn invalid_boundary_normal_display() {
        let e = Error::invalid_boundary_normal("inlet_top", [0.0, 0.0, 0.0], 0.0);
        let s = e.to_string();
        assert!(s.contains("inlet_top"));
        assert!(e.is_open_boundary_error());
    }

    #[test]
    fn check_boundary_normal_ok() {
        let r = Error::check_boundary_normal("bc", [1.0, 0.0, 0.0], 1e-10);
        assert!(r.is_ok());
    }

    #[test]
    fn check_boundary_normal_zero_fails() {
        let r = Error::check_boundary_normal("bc", [0.0, 0.0, 0.0], 1e-10);
        assert!(r.is_err());
        assert!(r.unwrap_err().is_open_boundary_error());
    }

    #[test]
    fn check_boundary_normal_tiny_fails() {
        let r = Error::check_boundary_normal("bc", [1e-15, 0.0, 0.0], 1e-10);
        assert!(r.is_err());
    }

    // ---- Turbulence errors -------------------------------------------------

    #[test]
    fn turbulence_model_error_display() {
        let e = Error::turbulence_model_error("Smagorinsky", "Cs constant is negative");
        assert!(e.to_string().contains("Smagorinsky"));
        assert!(e.is_turbulence_error());
        assert_eq!(e.category(), "turbulence");
    }

    #[test]
    fn negative_tke_display() {
        let e = Error::negative_tke(7, "k-omega", -0.5);
        let s = e.to_string();
        assert!(s.contains("k-omega"));
        assert!(s.contains("7"));
        assert!(e.is_turbulence_error());
        assert!(e.is_recoverable());
    }

    #[test]
    fn non_positive_dissipation_display() {
        let e = Error::non_positive_dissipation(3, "k-omega", 0.0);
        assert!(e.to_string().contains("k-omega"));
        assert!(e.is_turbulence_error());
        assert!(e.is_recoverable());
    }

    #[test]
    fn check_tke_positive_ok() {
        let r = Error::check_tke(0, "les", 0.5);
        assert!(r.is_ok());
        assert_eq!(r.unwrap(), 0.5);
    }

    #[test]
    fn check_tke_zero_ok() {
        let r = Error::check_tke(0, "les", 0.0);
        assert!(r.is_ok());
    }

    #[test]
    fn check_tke_negative_err() {
        let r = Error::check_tke(5, "Smagorinsky", -1e-3);
        assert!(r.is_err());
        assert!(r.unwrap_err().is_turbulence_error());
    }

    #[test]
    fn check_dissipation_positive_ok() {
        let r = Error::check_dissipation(0, "k-omega", 1e-4);
        assert!(r.is_ok());
    }

    #[test]
    fn check_dissipation_zero_err() {
        let r = Error::check_dissipation(1, "k-omega", 0.0);
        assert!(r.is_err());
    }

    #[test]
    fn check_dissipation_negative_err() {
        let r = Error::check_dissipation(2, "DES", -1.0);
        assert!(r.is_err());
        assert!(r.unwrap_err().is_turbulence_error());
    }

    // ---- Coupling errors ---------------------------------------------------

    #[test]
    fn empty_coupling_target_display() {
        let e = Error::empty_coupling_target("boat_hull", "faces");
        assert!(e.to_string().contains("boat_hull"));
        assert!(e.is_coupling_error());
        assert_eq!(e.category(), "coupling");
    }

    #[test]
    fn interface_transfer_not_converged_display() {
        let e = Error::interface_transfer_not_converged("fluid_solid_if", 30, 1.2e-2);
        let s = e.to_string();
        assert!(s.contains("fluid_solid_if"));
        assert!(s.contains("30"));
        assert!(e.is_coupling_error());
        assert!(e.is_convergence_failure());
        assert!(e.is_recoverable());
    }

    #[test]
    fn invalid_coupling_volume_display() {
        let e = Error::invalid_coupling_volume(12, -1e-6);
        let s = e.to_string();
        assert!(s.contains("12"));
        assert!(e.is_coupling_error());
    }

    #[test]
    fn non_finite_coupling_force_display() {
        let e = Error::non_finite_coupling_force("ship", [f64::NAN, 0.0, 0.0]);
        assert!(e.to_string().contains("ship"));
        assert!(e.is_coupling_error());
    }

    // ---- Kernel errors -----------------------------------------------------

    #[test]
    fn non_positive_smoothing_length_display() {
        let e = Error::non_positive_smoothing_length(5, 0.0);
        assert!(e.to_string().contains("5"));
        assert_eq!(e.category(), "kernel");
    }

    #[test]
    fn check_smoothing_length_positive_ok() {
        let r = Error::check_smoothing_length(0, 0.05);
        assert!(r.is_ok());
        assert!((r.unwrap() - 0.05).abs() < 1e-15);
    }

    #[test]
    fn check_smoothing_length_zero_err() {
        let r = Error::check_smoothing_length(0, 0.0);
        assert!(r.is_err());
    }

    #[test]
    fn check_smoothing_length_negative_err() {
        let r = Error::check_smoothing_length(3, -0.01);
        assert!(r.is_err());
    }

    #[test]
    fn zero_kernel_normalization_display() {
        let e = Error::ZeroKernelNormalization {
            kernel_name: "Poly6".into(),
        };
        assert!(e.to_string().contains("Poly6"));
        assert_eq!(e.category(), "kernel");
    }

    // ---- Neighbor search errors --------------------------------------------

    #[test]
    fn invalid_cell_size_display() {
        let e = Error::InvalidCellSize { cell_size: -0.1 };
        assert!(e.to_string().contains("-0.1") || e.to_string().len() > 10);
        assert_eq!(e.category(), "neighbor_search");
    }

    #[test]
    fn too_many_neighbors_display() {
        let e = Error::TooManyNeighbors {
            particle_index: 10,
            n_neighbors: 5000,
            threshold: 200,
        };
        let s = e.to_string();
        assert!(s.contains("5000"));
        assert!(s.contains("200"));
        assert_eq!(e.category(), "neighbor_search");
    }

    // ---- Density errors ----------------------------------------------------

    #[test]
    fn non_positive_density_display() {
        let e = Error::non_positive_density(8, -10.0);
        let s = e.to_string();
        assert!(s.contains("8"));
        assert_eq!(e.category(), "density");
    }

    #[test]
    fn check_density_positive_ok() {
        let r = Error::check_density(0, 1000.0);
        assert!(r.is_ok());
    }

    #[test]
    fn check_density_zero_err() {
        let r = Error::check_density(1, 0.0);
        assert!(r.is_err());
    }

    #[test]
    fn check_density_negative_err() {
        let r = Error::check_density(2, -5.0);
        assert!(r.is_err());
    }

    // ---- Surface errors ----------------------------------------------------

    #[test]
    fn surface_reconstruction_failed_display() {
        let e = Error::surface_reconstruction_failed(99, "not enough neighbors");
        let s = e.to_string();
        assert!(s.contains("99"));
        assert!(s.contains("not enough neighbors"));
        assert_eq!(e.category(), "surface");
    }

    // ---- I/O errors --------------------------------------------------------

    #[test]
    fn io_error_display() {
        let e = Error::io_error("file not found: particles.bin");
        assert!(e.to_string().contains("particles.bin"));
        assert_eq!(e.category(), "io");
    }

    // ---- is_recoverable comprehensive -------------------------------------

    #[test]
    fn is_recoverable_variants() {
        let e1 = Error::not_converged("PCISPH", 50, 0.1, 1e-3);
        assert!(e1.is_recoverable());

        let e2 = Error::invalid_param("dt", -1.0, 0.0, 1.0);
        assert!(!e2.is_recoverable());

        let e3 = Error::TimestepTooSmall {
            dt: 1e-15,
            dt_min: 1e-8,
            max_speed: 100.0,
        };
        assert!(e3.is_recoverable());

        let e4 = Error::non_finite("vel", 0, f64::NAN);
        assert!(e4.is_recoverable());

        let e5 = Error::negative_tke(0, "k-eps", -1.0);
        assert!(e5.is_recoverable());

        let e6 = Error::non_positive_dissipation(0, "k-eps", -1.0);
        assert!(e6.is_recoverable());
    }

    // ---- category exhaustive -----------------------------------------------

    #[test]
    fn category_general() {
        assert_eq!(Error::general("x").category(), "general");
    }

    #[test]
    fn category_io() {
        assert_eq!(Error::io_error("disk full").category(), "io");
    }

    #[test]
    fn category_coupling_empty_target() {
        assert_eq!(
            Error::empty_coupling_target("m", "v").category(),
            "coupling"
        );
    }

    // ---- boundary normal normalization helper ------------------------------

    #[test]
    fn check_boundary_normal_diagonal() {
        let n = [1.0_f64, 1.0, 1.0];
        let r = Error::check_boundary_normal("diag", n, 1e-10);
        assert!(r.is_ok());
    }
}
