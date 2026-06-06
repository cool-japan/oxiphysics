//! # OptimizerConfig - Trait Implementations
//!
//! This module contains trait implementations for `OptimizerConfig`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::OptimizerConfig;

impl Default for OptimizerConfig {
    fn default() -> Self {
        Self {
            max_iterations: 1000,
            time_step: 0.01,
            tolerance: 1e-6,
            damping: 0.5,
            stiffness: 1.0,
        }
    }
}
