//! # SolverHints - Trait Implementations
//!
//! This module contains trait implementations for `SolverHints`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::SolverHints;

impl Default for SolverHints {
    fn default() -> Self {
        SolverHints {
            max_velocity_iterations: None,
            needs_position_solve: true,
            is_unilateral: false,
            priority_override: None,
        }
    }
}
