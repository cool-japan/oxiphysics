// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Error types for oxiphysics-core
//!
//! Provides a unified error hierarchy covering all physics domains (FEM, LBM,
//! SPH, MD), with error context chaining, recovery suggestions, and diagnostic
//! information.

use thiserror::Error;

/// Main error type for the core module (legacy)
#[derive(Debug, Error)]
pub enum Error {
    /// Generic error
    #[error("{0}")]
    General(String),
}

/// Result type alias (legacy)
pub type Result<T> = std::result::Result<T, Error>;

/// Unified error type for the physics engine.
#[derive(Debug, Clone, PartialEq)]
pub enum PhysicsError {
    /// Solver diverged.
    NumericalDivergence {
        /// Description of the divergence.
        message: String,
    },
    /// Invalid parameter.
    InvalidInput {
        /// Name of the field that was invalid.
        field: String,
        /// Reason it was invalid.
        reason: String,
    },
    /// Index out of bounds.
    OutOfBounds {
        /// The index that was accessed.
        index: usize,
        /// The exclusive upper bound.
        max: usize,
    },
    /// Iterative solver failed to converge.
    ConvergenceFailed {
        /// Number of iterations performed.
        iterations: usize,
        /// Target tolerance.
        tolerance: f64,
        /// Final residual achieved.
        residual: f64,
    },
    /// Mesh topology error.
    MeshError {
        /// Description of the mesh error.
        message: String,
    },
    /// Invalid material state.
    MaterialError {
        /// Description of the material error.
        message: String,
    },
    /// Collision detection failure.
    CollisionError {
        /// Description of the collision error.
        message: String,
    },
    /// I/O error wrapper.
    IoError {
        /// Description of the I/O error.
        message: String,
    },

    // ── FEM domain errors ────────────────────────────────────────────────
    /// FEM element error (degenerate element, negative Jacobian, etc.)
    FemElementError {
        /// Element index.
        element_id: usize,
        /// Description.
        message: String,
    },

    /// FEM assembly error (stiffness matrix assembly failed).
    FemAssemblyError {
        /// Description.
        message: String,
    },

    /// FEM boundary condition error.
    FemBoundaryConditionError {
        /// Node or DOF index.
        node_id: usize,
        /// Description.
        message: String,
    },

    // ── LBM domain errors ────────────────────────────────────────────────
    /// LBM lattice configuration error.
    LbmLatticeError {
        /// Description.
        message: String,
    },

    /// LBM distribution function error (negative density, non-physical state).
    LbmDistributionError {
        /// Cell index.
        cell_id: usize,
        /// Description.
        message: String,
    },

    /// LBM stability violation (CFL, Mach number exceeded).
    LbmStabilityError {
        /// The violated parameter name.
        parameter: String,
        /// The value that caused the violation.
        value: f64,
        /// The allowed limit.
        limit: f64,
    },

    // ── SPH domain errors ────────────────────────────────────────────────
    /// SPH kernel error (radius too small, invalid kernel type).
    SphKernelError {
        /// Description.
        message: String,
    },

    /// SPH neighbor search error.
    SphNeighborError {
        /// Particle index.
        particle_id: usize,
        /// Description.
        message: String,
    },

    /// SPH density computation error (negative density, tensile instability).
    SphDensityError {
        /// Particle index.
        particle_id: usize,
        /// Computed density value.
        density: f64,
    },

    // ── MD domain errors ─────────────────────────────────────────────────
    /// MD force field error (missing parameters, invalid potential).
    MdForceFieldError {
        /// Atom type pair.
        atom_types: (String, String),
        /// Description.
        message: String,
    },

    /// MD integration error (energy drift, temperature blowup).
    MdIntegrationError {
        /// Time step at which the error occurred.
        step: usize,
        /// Description.
        message: String,
    },

    /// MD neighbor list error.
    MdNeighborListError {
        /// Description.
        message: String,
    },

    // ── Chained / contextual error ───────────────────────────────────────
    /// Error with additional context.
    WithContext {
        /// The original error.
        source: Box<PhysicsError>,
        /// Additional context message.
        context: String,
    },
}

impl std::fmt::Display for PhysicsError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            PhysicsError::NumericalDivergence { message } => {
                write!(f, "Numerical divergence: {}", message)
            }
            PhysicsError::InvalidInput { field, reason } => {
                write!(f, "Invalid input for '{}': {}", field, reason)
            }
            PhysicsError::OutOfBounds { index, max } => {
                write!(f, "Index {} out of bounds (max {})", index, max)
            }
            PhysicsError::ConvergenceFailed {
                iterations,
                tolerance,
                residual,
            } => write!(
                f,
                "Convergence failed after {} iterations (tolerance={}, residual={})",
                iterations, tolerance, residual
            ),
            PhysicsError::MeshError { message } => write!(f, "Mesh error: {}", message),
            PhysicsError::MaterialError { message } => write!(f, "Material error: {}", message),
            PhysicsError::CollisionError { message } => write!(f, "Collision error: {}", message),
            PhysicsError::IoError { message } => write!(f, "I/O error: {}", message),

            // FEM
            PhysicsError::FemElementError {
                element_id,
                message,
            } => write!(f, "FEM element {} error: {}", element_id, message),
            PhysicsError::FemAssemblyError { message } => {
                write!(f, "FEM assembly error: {}", message)
            }
            PhysicsError::FemBoundaryConditionError { node_id, message } => {
                write!(f, "FEM BC error at node {}: {}", node_id, message)
            }

            // LBM
            PhysicsError::LbmLatticeError { message } => {
                write!(f, "LBM lattice error: {}", message)
            }
            PhysicsError::LbmDistributionError { cell_id, message } => {
                write!(f, "LBM distribution error at cell {}: {}", cell_id, message)
            }
            PhysicsError::LbmStabilityError {
                parameter,
                value,
                limit,
            } => write!(
                f,
                "LBM stability violation: {} = {} exceeds limit {}",
                parameter, value, limit
            ),

            // SPH
            PhysicsError::SphKernelError { message } => {
                write!(f, "SPH kernel error: {}", message)
            }
            PhysicsError::SphNeighborError {
                particle_id,
                message,
            } => write!(
                f,
                "SPH neighbor error for particle {}: {}",
                particle_id, message
            ),
            PhysicsError::SphDensityError {
                particle_id,
                density,
            } => write!(
                f,
                "SPH density error: particle {} has invalid density {}",
                particle_id, density
            ),

            // MD
            PhysicsError::MdForceFieldError {
                atom_types,
                message,
            } => write!(
                f,
                "MD force field error for pair ({}, {}): {}",
                atom_types.0, atom_types.1, message
            ),
            PhysicsError::MdIntegrationError { step, message } => {
                write!(f, "MD integration error at step {}: {}", step, message)
            }
            PhysicsError::MdNeighborListError { message } => {
                write!(f, "MD neighbor list error: {}", message)
            }

            // Context chain
            PhysicsError::WithContext { source, context } => {
                write!(f, "{}: {}", context, source)
            }
        }
    }
}

impl std::error::Error for PhysicsError {}

impl PhysicsError {
    /// Wrap this error with additional context.
    pub fn with_context(self, context: impl Into<String>) -> Self {
        PhysicsError::WithContext {
            source: Box::new(self),
            context: context.into(),
        }
    }

    /// Get the root cause error, unwrapping all context layers.
    pub fn root_cause(&self) -> &PhysicsError {
        match self {
            PhysicsError::WithContext { source, .. } => source.root_cause(),
            other => other,
        }
    }

    /// Get the chain of context messages from outermost to innermost.
    pub fn context_chain(&self) -> Vec<&str> {
        let mut chain = Vec::new();
        let mut current = self;
        while let PhysicsError::WithContext { source, context } = current {
            chain.push(context.as_str());
            current = source;
        }
        chain
    }

    /// Returns the domain name for this error variant.
    pub fn domain(&self) -> &'static str {
        match self {
            PhysicsError::FemElementError { .. }
            | PhysicsError::FemAssemblyError { .. }
            | PhysicsError::FemBoundaryConditionError { .. } => "FEM",

            PhysicsError::LbmLatticeError { .. }
            | PhysicsError::LbmDistributionError { .. }
            | PhysicsError::LbmStabilityError { .. } => "LBM",

            PhysicsError::SphKernelError { .. }
            | PhysicsError::SphNeighborError { .. }
            | PhysicsError::SphDensityError { .. } => "SPH",

            PhysicsError::MdForceFieldError { .. }
            | PhysicsError::MdIntegrationError { .. }
            | PhysicsError::MdNeighborListError { .. } => "MD",

            PhysicsError::MeshError { .. } => "Mesh",
            PhysicsError::MaterialError { .. } => "Material",
            PhysicsError::CollisionError { .. } => "Collision",
            PhysicsError::IoError { .. } => "IO",

            PhysicsError::NumericalDivergence { .. } => "Numerical",
            PhysicsError::InvalidInput { .. } => "Input",
            PhysicsError::OutOfBounds { .. } => "Bounds",
            PhysicsError::ConvergenceFailed { .. } => "Solver",

            PhysicsError::WithContext { source, .. } => source.domain(),
        }
    }

    /// Returns a recovery suggestion for this error.
    pub fn recovery_suggestion(&self) -> &'static str {
        match self {
            PhysicsError::NumericalDivergence { .. } => {
                "Try reducing the time step or using a more stable integration scheme."
            }
            PhysicsError::InvalidInput { .. } => {
                "Check the input parameters against the documentation."
            }
            PhysicsError::OutOfBounds { .. } => "Verify array sizes and index calculations.",
            PhysicsError::ConvergenceFailed { .. } => {
                "Try increasing max iterations, relaxing tolerance, or improving the initial guess."
            }
            PhysicsError::MeshError { .. } => {
                "Check mesh quality metrics and repair degenerate elements."
            }
            PhysicsError::MaterialError { .. } => {
                "Verify material parameters (e.g., positive moduli, valid Poisson ratio)."
            }
            PhysicsError::CollisionError { .. } => {
                "Check collision geometry and broadphase settings."
            }
            PhysicsError::IoError { .. } => "Check file paths, permissions, and disk space.",
            PhysicsError::FemElementError { .. } => {
                "Check element connectivity and node positions for degenerate configurations."
            }
            PhysicsError::FemAssemblyError { .. } => {
                "Verify element stiffness matrices are positive semi-definite."
            }
            PhysicsError::FemBoundaryConditionError { .. } => {
                "Ensure boundary conditions are consistent and the system is not over-constrained."
            }
            PhysicsError::LbmLatticeError { .. } => {
                "Check lattice dimensions and velocity model configuration."
            }
            PhysicsError::LbmDistributionError { .. } => {
                "Reduce the time step or check inlet/outlet boundary conditions."
            }
            PhysicsError::LbmStabilityError { .. } => {
                "Reduce flow velocity or increase lattice resolution to satisfy stability constraints."
            }
            PhysicsError::SphKernelError { .. } => {
                "Check smoothing length and kernel radius parameters."
            }
            PhysicsError::SphNeighborError { .. } => {
                "Increase neighbor search radius or check spatial hashing configuration."
            }
            PhysicsError::SphDensityError { .. } => {
                "Add tensile instability correction or increase particle count in sparse regions."
            }
            PhysicsError::MdForceFieldError { .. } => {
                "Check force field parameter files and atom type definitions."
            }
            PhysicsError::MdIntegrationError { .. } => {
                "Reduce time step, check for overlapping atoms, or use a thermostat."
            }
            PhysicsError::MdNeighborListError { .. } => {
                "Increase skin distance or rebuild neighbor list more frequently."
            }
            PhysicsError::WithContext { source, .. } => source.recovery_suggestion(),
        }
    }

    /// Returns the severity level of this error.
    pub fn severity(&self) -> ErrorSeverity {
        match self {
            PhysicsError::NumericalDivergence { .. } => ErrorSeverity::Critical,
            PhysicsError::ConvergenceFailed { .. } => ErrorSeverity::Warning,
            PhysicsError::InvalidInput { .. } => ErrorSeverity::Error,
            PhysicsError::OutOfBounds { .. } => ErrorSeverity::Error,
            PhysicsError::LbmStabilityError { .. } => ErrorSeverity::Critical,
            PhysicsError::MdIntegrationError { .. } => ErrorSeverity::Critical,
            PhysicsError::SphDensityError { .. } => ErrorSeverity::Warning,
            PhysicsError::FemElementError { .. } => ErrorSeverity::Error,
            PhysicsError::WithContext { source, .. } => source.severity(),
            _ => ErrorSeverity::Error,
        }
    }
}

/// Severity level for physics errors.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum ErrorSeverity {
    /// Informational — the simulation can continue.
    Info,
    /// Warning — the simulation may produce inaccurate results.
    Warning,
    /// Error — the current operation cannot complete.
    Error,
    /// Critical — the entire simulation state may be corrupt.
    Critical,
}

impl std::fmt::Display for ErrorSeverity {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ErrorSeverity::Info => write!(f, "INFO"),
            ErrorSeverity::Warning => write!(f, "WARNING"),
            ErrorSeverity::Error => write!(f, "ERROR"),
            ErrorSeverity::Critical => write!(f, "CRITICAL"),
        }
    }
}

/// Result alias using `PhysicsError`.
pub type PhysicsResult<T> = std::result::Result<T, PhysicsError>;

/// Returns `Ok(value)` if `value > 0`, otherwise `Err(InvalidInput)`.
pub fn check_positive(value: f64, name: &str) -> PhysicsResult<f64> {
    if value > 0.0 {
        Ok(value)
    } else {
        Err(PhysicsError::InvalidInput {
            field: name.to_string(),
            reason: format!("expected positive value, got {}", value),
        })
    }
}

/// Returns `Ok(value)` if `value >= 0`, otherwise `Err(InvalidInput)`.
pub fn check_non_negative(value: f64, name: &str) -> PhysicsResult<f64> {
    if value >= 0.0 {
        Ok(value)
    } else {
        Err(PhysicsError::InvalidInput {
            field: name.to_string(),
            reason: format!("expected non-negative value, got {}", value),
        })
    }
}

/// Returns `Ok(value)` if `min <= value <= max`, otherwise `Err(InvalidInput)`.
pub fn check_range(value: f64, min: f64, max: f64, name: &str) -> PhysicsResult<f64> {
    if value >= min && value <= max {
        Ok(value)
    } else {
        Err(PhysicsError::InvalidInput {
            field: name.to_string(),
            reason: format!("expected value in [{}, {}], got {}", min, max, value),
        })
    }
}

/// Returns `Ok(value)` if finite, otherwise `Err(NumericalDivergence)`.
pub fn check_finite(value: f64, name: &str) -> PhysicsResult<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(PhysicsError::NumericalDivergence {
            message: format!("'{}' is not finite: {}", name, value),
        })
    }
}

/// Returns `Ok(index)` if `index < max`, otherwise `Err(OutOfBounds)`.
pub fn check_index(index: usize, max: usize) -> PhysicsResult<usize> {
    if index < max {
        Ok(index)
    } else {
        Err(PhysicsError::OutOfBounds { index, max })
    }
}

/// Check that a 3D vector has finite components.
pub fn check_finite_vec3(v: [f64; 3], name: &str) -> PhysicsResult<[f64; 3]> {
    for (i, &c) in v.iter().enumerate() {
        if !c.is_finite() {
            return Err(PhysicsError::NumericalDivergence {
                message: format!("'{}'[{}] is not finite: {}", name, i, c),
            });
        }
    }
    Ok(v)
}

/// Check that a slice contains no NaN or infinite values.
pub fn check_finite_slice(values: &[f64], name: &str) -> PhysicsResult<()> {
    for (i, &v) in values.iter().enumerate() {
        if !v.is_finite() {
            return Err(PhysicsError::NumericalDivergence {
                message: format!("'{}'[{}] is not finite: {}", name, i, v),
            });
        }
    }
    Ok(())
}

/// Check Poisson's ratio is in the valid range (-1, 0.5) for 3D.
pub fn check_poisson_ratio(nu: f64) -> PhysicsResult<f64> {
    if nu > -1.0 && nu < 0.5 {
        Ok(nu)
    } else {
        Err(PhysicsError::MaterialError {
            message: format!("Poisson's ratio must be in (-1, 0.5) for 3D, got {}", nu),
        })
    }
}

/// Check that the LBM Mach number is below the stability limit.
pub fn check_lbm_mach(velocity: f64, cs: f64, limit: f64) -> PhysicsResult<f64> {
    let mach = velocity / cs;
    if mach <= limit {
        Ok(mach)
    } else {
        Err(PhysicsError::LbmStabilityError {
            parameter: "Mach number".to_string(),
            value: mach,
            limit,
        })
    }
}

/// Check that SPH density is positive.
pub fn check_sph_density(particle_id: usize, density: f64) -> PhysicsResult<f64> {
    if density > 0.0 {
        Ok(density)
    } else {
        Err(PhysicsError::SphDensityError {
            particle_id,
            density,
        })
    }
}

/// Check that an MD time step is reasonable (positive and not too large).
pub fn check_md_timestep(dt: f64, max_dt: f64) -> PhysicsResult<f64> {
    if dt <= 0.0 {
        Err(PhysicsError::InvalidInput {
            field: "timestep".to_string(),
            reason: format!("must be positive, got {}", dt),
        })
    } else if dt > max_dt {
        Err(PhysicsError::MdIntegrationError {
            step: 0,
            message: format!("timestep {} exceeds maximum {}", dt, max_dt),
        })
    } else {
        Ok(dt)
    }
}

/// Tracks solver convergence history and diagnostics.
#[derive(Debug, Clone)]
pub struct SolverDiagnostics {
    /// Number of iterations performed.
    pub iterations: usize,
    /// Final residual after last iteration.
    pub final_residual: f64,
    /// Whether the solver converged.
    pub converged: bool,
    /// Per-iteration residual history.
    pub history: Vec<f64>,
}

impl SolverDiagnostics {
    /// Creates a new, empty `SolverDiagnostics`.
    pub fn new() -> Self {
        Self {
            iterations: 0,
            final_residual: f64::INFINITY,
            converged: false,
            history: Vec::new(),
        }
    }

    /// Records a residual for the current iteration.
    pub fn record(&mut self, residual: f64) {
        self.history.push(residual);
        self.iterations = self.history.len();
        self.final_residual = residual;
    }

    /// Marks convergence and returns `Ok(())`, or `Err(ConvergenceFailed)`.
    pub fn check_convergence(&self, tolerance: f64) -> PhysicsResult<()> {
        if self.converged || self.final_residual <= tolerance {
            Ok(())
        } else {
            Err(PhysicsError::ConvergenceFailed {
                iterations: self.iterations,
                tolerance,
                residual: self.final_residual,
            })
        }
    }

    /// Returns the geometric-mean convergence rate from the last 5 residuals,
    /// or `None` if fewer than 2 residuals have been recorded.
    pub fn convergence_rate(&self) -> Option<f64> {
        let n = self.history.len();
        if n < 2 {
            return None;
        }
        let window_size = 5.min(n);
        let window = &self.history[n - window_size..];
        // Geometric mean of successive ratios: (r_{i+1}/r_i)
        let steps = window_size - 1;
        let log_sum: f64 = window
            .windows(2)
            .map(|w| {
                if w[0].abs() < f64::EPSILON {
                    0.0
                } else {
                    (w[1].abs() / w[0].abs()).ln()
                }
            })
            .sum();
        Some((log_sum / steps as f64).exp())
    }

    /// Check if the solver is diverging (residual increasing over the last few iterations).
    pub fn is_diverging(&self) -> bool {
        let n = self.history.len();
        if n < 3 {
            return false;
        }
        // Check if the last 3 residuals are monotonically increasing
        let last3 = &self.history[n - 3..];
        last3[0] < last3[1] && last3[1] < last3[2]
    }

    /// Returns a summary string for diagnostics output.
    pub fn summary(&self) -> String {
        let status = if self.converged {
            "CONVERGED"
        } else if self.is_diverging() {
            "DIVERGING"
        } else {
            "NOT CONVERGED"
        };
        format!(
            "[{}] {} iterations, residual = {:.6e}",
            status, self.iterations, self.final_residual
        )
    }

    /// Reset the diagnostics for a new solve.
    pub fn reset(&mut self) {
        self.iterations = 0;
        self.final_residual = f64::INFINITY;
        self.converged = false;
        self.history.clear();
    }
}

impl Default for SolverDiagnostics {
    fn default() -> Self {
        Self::new()
    }
}

/// Diagnostic information attached to an error for debugging.
#[derive(Debug, Clone)]
pub struct DiagnosticInfo {
    /// The error that triggered the diagnostic.
    pub error: PhysicsError,
    /// Time step or simulation time at which the error occurred.
    pub time: Option<f64>,
    /// Step number at which the error occurred.
    pub step: Option<usize>,
    /// Solver diagnostics if available.
    pub solver_diag: Option<SolverDiagnostics>,
    /// Additional key-value metadata.
    pub metadata: Vec<(String, String)>,
}

impl DiagnosticInfo {
    /// Create diagnostic info from an error.
    pub fn from_error(error: PhysicsError) -> Self {
        Self {
            error,
            time: None,
            step: None,
            solver_diag: None,
            metadata: Vec::new(),
        }
    }

    /// Attach a simulation time.
    pub fn at_time(mut self, t: f64) -> Self {
        self.time = Some(t);
        self
    }

    /// Attach a step number.
    pub fn at_step(mut self, step: usize) -> Self {
        self.step = Some(step);
        self
    }

    /// Attach solver diagnostics.
    pub fn with_solver_diag(mut self, diag: SolverDiagnostics) -> Self {
        self.solver_diag = Some(diag);
        self
    }

    /// Add a metadata key-value pair.
    pub fn with_metadata(mut self, key: impl Into<String>, value: impl Into<String>) -> Self {
        self.metadata.push((key.into(), value.into()));
        self
    }

    /// Format a full diagnostic report.
    pub fn report(&self) -> String {
        let mut lines = Vec::new();
        lines.push(format!("[{}] {}", self.error.severity(), self.error));
        lines.push(format!("Domain: {}", self.error.domain()));
        lines.push(format!("Suggestion: {}", self.error.recovery_suggestion()));
        if let Some(t) = self.time {
            lines.push(format!("Time: {:.6e}", t));
        }
        if let Some(s) = self.step {
            lines.push(format!("Step: {}", s));
        }
        if let Some(ref diag) = self.solver_diag {
            lines.push(format!("Solver: {}", diag.summary()));
        }
        for (k, v) in &self.metadata {
            lines.push(format!("{}: {}", k, v));
        }
        lines.join("\n")
    }
}

/// Collect multiple errors during a batch operation.
#[derive(Debug, Clone, Default)]
pub struct ErrorCollector {
    /// Collected errors.
    pub errors: Vec<PhysicsError>,
    /// Maximum number of errors to collect before stopping.
    pub max_errors: usize,
}

impl ErrorCollector {
    /// Create a new collector with default capacity (100).
    pub fn new() -> Self {
        Self {
            errors: Vec::new(),
            max_errors: 100,
        }
    }

    /// Create a collector with a specific max error limit.
    pub fn with_limit(max: usize) -> Self {
        Self {
            errors: Vec::new(),
            max_errors: max,
        }
    }

    /// Add an error. Returns `true` if the error limit has been reached.
    pub fn push(&mut self, error: PhysicsError) -> bool {
        self.errors.push(error);
        self.errors.len() >= self.max_errors
    }

    /// Returns `true` if any errors were collected.
    pub fn has_errors(&self) -> bool {
        !self.errors.is_empty()
    }

    /// Number of collected errors.
    pub fn count(&self) -> usize {
        self.errors.len()
    }

    /// Drain all errors, returning them.
    pub fn drain(&mut self) -> Vec<PhysicsError> {
        std::mem::take(&mut self.errors)
    }

    /// Convert into a result: Ok if no errors, Err with first error otherwise.
    pub fn into_result(self) -> PhysicsResult<()> {
        if self.errors.is_empty() {
            Ok(())
        } else {
            let count = self.count();
            let first = self
                .errors
                .into_iter()
                .next()
                .expect("errors is non-empty after check");
            Err(first.with_context(format!("First of {} errors", count)))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_check_positive() {
        assert!(check_positive(1.0, "x").is_ok());
        assert!(check_positive(-1.0, "x").is_err());
        assert!(check_positive(0.0, "x").is_err());
    }

    #[test]
    fn test_check_range() {
        assert!(check_range(0.5, 0.0, 1.0, "t").is_ok());
        assert!(check_range(0.0, 0.0, 1.0, "t").is_ok());
        assert!(check_range(1.0, 0.0, 1.0, "t").is_ok());
        assert!(check_range(-0.1, 0.0, 1.0, "t").is_err());
        assert!(check_range(1.1, 0.0, 1.0, "t").is_err());
    }

    #[test]
    fn test_check_finite() {
        assert!(check_finite(f64::NAN, "v").is_err());
        assert!(check_finite(f64::INFINITY, "v").is_err());
        assert!(check_finite(f64::NEG_INFINITY, "v").is_err());
        assert!(check_finite(1.0, "v").is_ok());
    }

    #[test]
    fn test_check_index() {
        assert!(check_index(0, 5).is_ok());
        assert!(check_index(4, 5).is_ok());
        assert!(check_index(5, 5).is_err());
        assert!(check_index(100, 5).is_err());
    }

    #[test]
    fn test_solver_diagnostics_convergence() {
        let mut diag = SolverDiagnostics::new();
        let tolerance = 1e-6;
        // Simulate 10 iterations converging from 1.0 down to below tolerance
        let mut residual = 1.0_f64;
        for _ in 0..10 {
            residual *= 0.1;
            diag.record(residual);
        }
        // After 10 iterations residual is 1e-10, well below 1e-6
        assert!(diag.check_convergence(tolerance).is_ok());
    }

    #[test]
    fn test_solver_diagnostics_not_converged() {
        let mut diag = SolverDiagnostics::new();
        diag.record(1.0);
        diag.record(0.9);
        assert!(diag.check_convergence(1e-6).is_err());
    }

    #[test]
    fn test_convergence_rate_none_when_too_few() {
        let mut diag = SolverDiagnostics::new();
        assert!(diag.convergence_rate().is_none());
        diag.record(1.0);
        assert!(diag.convergence_rate().is_none());
    }

    #[test]
    fn test_convergence_rate_some() {
        let mut diag = SolverDiagnostics::new();
        for i in 1..=6 {
            diag.record(1.0 / (10_f64.powi(i)));
        }
        let rate = diag.convergence_rate();
        assert!(rate.is_some());
        // Each step reduces by factor 10, so rate ≈ 0.1
        let r = rate.unwrap();
        assert!((r - 0.1).abs() < 1e-9, "rate was {}", r);
    }

    #[test]
    fn test_display_physics_error() {
        let e = PhysicsError::InvalidInput {
            field: "density".to_string(),
            reason: "must be positive".to_string(),
        };
        let s = format!("{}", e);
        assert!(s.contains("density"));
        assert!(s.contains("must be positive"));
    }

    // ── New tests ────────────────────────────────────────────────────────────

    #[test]
    fn test_check_non_negative() {
        assert!(check_non_negative(0.0, "x").is_ok());
        assert!(check_non_negative(1.0, "x").is_ok());
        assert!(check_non_negative(-0.01, "x").is_err());
    }

    #[test]
    fn test_check_finite_vec3() {
        assert!(check_finite_vec3([1.0, 2.0, 3.0], "v").is_ok());
        assert!(check_finite_vec3([f64::NAN, 0.0, 0.0], "v").is_err());
        assert!(check_finite_vec3([0.0, f64::INFINITY, 0.0], "v").is_err());
    }

    #[test]
    fn test_check_finite_slice() {
        assert!(check_finite_slice(&[1.0, 2.0, 3.0], "arr").is_ok());
        assert!(check_finite_slice(&[1.0, f64::NAN], "arr").is_err());
    }

    #[test]
    fn test_check_poisson_ratio() {
        assert!(check_poisson_ratio(0.3).is_ok());
        assert!(check_poisson_ratio(0.0).is_ok());
        assert!(check_poisson_ratio(-0.5).is_ok());
        assert!(check_poisson_ratio(0.5).is_err()); // boundary excluded
        assert!(check_poisson_ratio(-1.0).is_err()); // boundary excluded
        assert!(check_poisson_ratio(0.6).is_err());
    }

    #[test]
    fn test_check_lbm_mach() {
        assert!(check_lbm_mach(0.1, 1.0 / 3.0_f64.sqrt(), 0.3).is_ok());
        // High velocity triggers error
        let result = check_lbm_mach(10.0, 1.0 / 3.0_f64.sqrt(), 0.3);
        assert!(result.is_err());
        if let Err(PhysicsError::LbmStabilityError { parameter, .. }) = result {
            assert_eq!(parameter, "Mach number");
        }
    }

    #[test]
    fn test_check_sph_density() {
        assert!(check_sph_density(0, 1000.0).is_ok());
        assert!(check_sph_density(5, -1.0).is_err());
        assert!(check_sph_density(5, 0.0).is_err());
    }

    #[test]
    fn test_check_md_timestep() {
        assert!(check_md_timestep(0.001, 0.01).is_ok());
        assert!(check_md_timestep(-0.001, 0.01).is_err());
        assert!(check_md_timestep(0.1, 0.01).is_err());
    }

    #[test]
    fn test_error_with_context() {
        let e = PhysicsError::InvalidInput {
            field: "mass".to_string(),
            reason: "negative".to_string(),
        };
        let e2 = e.with_context("during particle initialization");
        let s = format!("{}", e2);
        assert!(s.contains("during particle initialization"));
        assert!(s.contains("mass"));
    }

    #[test]
    fn test_error_root_cause() {
        let e = PhysicsError::NumericalDivergence {
            message: "blowup".to_string(),
        };
        let e2 = e.with_context("in solver").with_context("in simulation");
        let root = e2.root_cause();
        match root {
            PhysicsError::NumericalDivergence { message } => {
                assert_eq!(message, "blowup");
            }
            other => panic!("expected NumericalDivergence, got {:?}", other),
        }
    }

    #[test]
    fn test_error_context_chain() {
        let e = PhysicsError::MeshError {
            message: "bad".to_string(),
        };
        let e2 = e.with_context("step 1").with_context("step 2");
        let chain = e2.context_chain();
        assert_eq!(chain.len(), 2);
        assert_eq!(chain[0], "step 2");
        assert_eq!(chain[1], "step 1");
    }

    #[test]
    fn test_error_domain() {
        assert_eq!(
            PhysicsError::FemElementError {
                element_id: 0,
                message: "x".into()
            }
            .domain(),
            "FEM"
        );
        assert_eq!(
            PhysicsError::LbmLatticeError {
                message: "x".into()
            }
            .domain(),
            "LBM"
        );
        assert_eq!(
            PhysicsError::SphKernelError {
                message: "x".into()
            }
            .domain(),
            "SPH"
        );
        assert_eq!(
            PhysicsError::MdForceFieldError {
                atom_types: ("A".into(), "B".into()),
                message: "x".into()
            }
            .domain(),
            "MD"
        );
    }

    #[test]
    fn test_error_severity() {
        assert_eq!(
            PhysicsError::NumericalDivergence {
                message: "x".into()
            }
            .severity(),
            ErrorSeverity::Critical
        );
        assert_eq!(
            PhysicsError::ConvergenceFailed {
                iterations: 10,
                tolerance: 1e-6,
                residual: 1e-3
            }
            .severity(),
            ErrorSeverity::Warning
        );
        assert_eq!(
            PhysicsError::InvalidInput {
                field: "x".into(),
                reason: "y".into()
            }
            .severity(),
            ErrorSeverity::Error
        );
    }

    #[test]
    fn test_error_recovery_suggestion() {
        let e = PhysicsError::ConvergenceFailed {
            iterations: 100,
            tolerance: 1e-6,
            residual: 1e-3,
        };
        let suggestion = e.recovery_suggestion();
        assert!(suggestion.contains("iterations") || suggestion.contains("tolerance"));
    }

    #[test]
    fn test_severity_display() {
        assert_eq!(format!("{}", ErrorSeverity::Critical), "CRITICAL");
        assert_eq!(format!("{}", ErrorSeverity::Warning), "WARNING");
        assert_eq!(format!("{}", ErrorSeverity::Error), "ERROR");
        assert_eq!(format!("{}", ErrorSeverity::Info), "INFO");
    }

    #[test]
    fn test_severity_ordering() {
        assert!(ErrorSeverity::Info < ErrorSeverity::Warning);
        assert!(ErrorSeverity::Warning < ErrorSeverity::Error);
        assert!(ErrorSeverity::Error < ErrorSeverity::Critical);
    }

    #[test]
    fn test_solver_diagnostics_is_diverging() {
        let mut diag = SolverDiagnostics::new();
        diag.record(1.0);
        diag.record(2.0);
        diag.record(3.0);
        assert!(diag.is_diverging());

        let mut diag2 = SolverDiagnostics::new();
        diag2.record(3.0);
        diag2.record(2.0);
        diag2.record(1.0);
        assert!(!diag2.is_diverging());
    }

    #[test]
    fn test_solver_diagnostics_summary() {
        let mut diag = SolverDiagnostics::new();
        diag.record(1.0);
        diag.converged = true;
        let s = diag.summary();
        assert!(s.contains("CONVERGED"));
    }

    #[test]
    fn test_solver_diagnostics_reset() {
        let mut diag = SolverDiagnostics::new();
        diag.record(1.0);
        diag.record(0.5);
        diag.converged = true;
        diag.reset();
        assert_eq!(diag.iterations, 0);
        assert!(diag.history.is_empty());
        assert!(!diag.converged);
    }

    #[test]
    fn test_diagnostic_info_report() {
        let e = PhysicsError::MdIntegrationError {
            step: 1000,
            message: "energy drift".to_string(),
        };
        let report = DiagnosticInfo::from_error(e)
            .at_time(1.5e-12)
            .at_step(1000)
            .with_metadata("ensemble", "NVT")
            .report();
        assert!(report.contains("CRITICAL"));
        assert!(report.contains("MD"));
        assert!(report.contains("energy drift"));
        assert!(report.contains("ensemble: NVT"));
    }

    #[test]
    fn test_error_collector_basic() {
        let mut collector = ErrorCollector::new();
        assert!(!collector.has_errors());
        assert_eq!(collector.count(), 0);

        collector.push(PhysicsError::InvalidInput {
            field: "x".into(),
            reason: "bad".into(),
        });
        assert!(collector.has_errors());
        assert_eq!(collector.count(), 1);
    }

    #[test]
    fn test_error_collector_limit() {
        let mut collector = ErrorCollector::with_limit(3);
        collector.push(PhysicsError::MeshError {
            message: "e1".into(),
        });
        collector.push(PhysicsError::MeshError {
            message: "e2".into(),
        });
        let limit_reached = collector.push(PhysicsError::MeshError {
            message: "e3".into(),
        });
        assert!(limit_reached);
    }

    #[test]
    fn test_error_collector_drain() {
        let mut collector = ErrorCollector::new();
        collector.push(PhysicsError::IoError {
            message: "test".into(),
        });
        let errors = collector.drain();
        assert_eq!(errors.len(), 1);
        assert!(!collector.has_errors());
    }

    #[test]
    fn test_error_collector_into_result_ok() {
        let collector = ErrorCollector::new();
        assert!(collector.into_result().is_ok());
    }

    #[test]
    fn test_error_collector_into_result_err() {
        let mut collector = ErrorCollector::new();
        collector.push(PhysicsError::IoError {
            message: "fail".into(),
        });
        let result = collector.into_result();
        assert!(result.is_err());
    }

    #[test]
    fn test_display_fem_errors() {
        let e = PhysicsError::FemElementError {
            element_id: 42,
            message: "negative Jacobian".to_string(),
        };
        let s = format!("{}", e);
        assert!(s.contains("42"));
        assert!(s.contains("negative Jacobian"));
    }

    #[test]
    fn test_display_lbm_errors() {
        let e = PhysicsError::LbmStabilityError {
            parameter: "CFL".to_string(),
            value: 1.5,
            limit: 1.0,
        };
        let s = format!("{}", e);
        assert!(s.contains("CFL"));
        assert!(s.contains("1.5"));
        assert!(s.contains("1"));
    }

    #[test]
    fn test_display_sph_errors() {
        let e = PhysicsError::SphDensityError {
            particle_id: 7,
            density: -0.5,
        };
        let s = format!("{}", e);
        assert!(s.contains("7"));
        assert!(s.contains("-0.5"));
    }

    #[test]
    fn test_display_md_errors() {
        let e = PhysicsError::MdForceFieldError {
            atom_types: ("C".to_string(), "O".to_string()),
            message: "missing LJ parameters".to_string(),
        };
        let s = format!("{}", e);
        assert!(s.contains("C"));
        assert!(s.contains("O"));
        assert!(s.contains("missing LJ parameters"));
    }

    #[test]
    fn test_context_chain_preserves_domain() {
        let e = PhysicsError::SphKernelError {
            message: "bad kernel".into(),
        };
        let e2 = e.with_context("in step 5");
        assert_eq!(e2.domain(), "SPH");
    }
}
