//! # ConstraintMetrics - Trait Implementations
//!
//! This module contains trait implementations for `ConstraintMetrics`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::ConstraintMetrics;

impl Default for ConstraintMetrics {
    fn default() -> Self {
        ConstraintMetrics {
            velocity_iterations: 0,
            velocity_residual: 0.0,
            position_correction: 0.0,
            accumulated_impulse: 0.0,
            converged: true,
        }
    }
}
