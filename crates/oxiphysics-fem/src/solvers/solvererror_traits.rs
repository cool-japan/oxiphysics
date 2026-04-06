//! # SolverError - Trait Implementations
//!
//! This module contains trait implementations for `SolverError`.
//!
//! ## Implemented Traits
//!
//! - `Display`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::SolverError;

impl std::fmt::Display for SolverError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SolverError::NotPositiveDefinite => {
                write!(f, "matrix is not positive definite")
            }
            SolverError::DimensionMismatch => write!(f, "dimension mismatch"),
        }
    }
}
