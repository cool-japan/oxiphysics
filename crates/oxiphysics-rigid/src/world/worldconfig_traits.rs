//! # WorldConfig - Trait Implementations
//!
//! This module contains trait implementations for `WorldConfig`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::{SolverConfig, WorldConfig};

impl Default for WorldConfig {
    fn default() -> Self {
        Self {
            gravity: [0.0, -9.81, 0.0],
            dt: 1.0 / 60.0,
            linear_sleep_threshold: 0.01,
            angular_sleep_threshold: 0.01,
            time_before_sleep: 0.5,
            solver: SolverConfig::default(),
        }
    }
}
