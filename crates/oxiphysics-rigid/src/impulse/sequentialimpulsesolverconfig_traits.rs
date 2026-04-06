//! # SequentialImpulseSolverConfig - Trait Implementations
//!
//! This module contains trait implementations for `SequentialImpulseSolverConfig`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::SequentialImpulseSolverConfig;

impl Default for SequentialImpulseSolverConfig {
    fn default() -> Self {
        Self {
            velocity_iterations: 10,
            baumgarte: 0.2,
            dt: 1.0 / 60.0,
            warm_start: true,
        }
    }
}
