// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Error types for oxiphysics-fem

use thiserror::Error;

/// Main error type for the fem module
#[derive(Debug, Error)]
pub enum Error {
    /// Generic error
    #[error("{0}")]
    General(String),
    /// An input argument has an invalid size, index or value.
    #[error("invalid input: {0}")]
    InvalidInput(String),
    /// A degree-of-freedom index lies outside the system it addresses.
    #[error("degree of freedom {dof} is out of range for a system with {n} DOFs")]
    DofOutOfRange {
        /// Offending (global) degree-of-freedom index.
        dof: usize,
        /// Number of degrees of freedom of the addressed system.
        n: usize,
    },
    /// A linear system is singular or numerically singular and cannot be solved.
    #[error("singular or numerically singular linear system")]
    SingularSystem,
    /// An iterative solver did not reach the requested tolerance.
    #[error(
        "{solver} did not converge in {iterations} iterations (relative residual {residual:e})"
    )]
    NotConverged {
        /// Name of the solver that failed.
        solver: String,
        /// Number of iterations performed.
        iterations: usize,
        /// Relative residual reached when the iteration stopped.
        residual: f64,
    },
}

/// Result type alias
pub type Result<T> = std::result::Result<T, Error>;
