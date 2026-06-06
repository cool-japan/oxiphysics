//! # TcsConfig - Trait Implementations
//!
//! This module contains trait implementations for `TcsConfig`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)
use super::types::TcsConfig;

impl Default for TcsConfig {
    fn default() -> Self {
        Self {
            slip_threshold: 0.1,
            torque_reduction_rate: 0.8,
            min_torque_fraction: 0.2,
        }
    }
}
