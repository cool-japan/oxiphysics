//! # SolverConfig - Trait Implementations
//!
//! This module contains trait implementations for `SolverConfig`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::SolverConfig;

impl Default for SolverConfig {
    fn default() -> Self {
        Self {
            velocity_iterations: 8,
            position_iterations: 3,
            baumgarte_factor: 0.2,
            restitution_threshold: 1.0,
        }
    }
}
