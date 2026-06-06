//! # AbsConfig - Trait Implementations
//!
//! This module contains trait implementations for `AbsConfig`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)
use super::types::AbsConfig;

impl Default for AbsConfig {
    fn default() -> Self {
        Self {
            slip_target: 0.12,
            pressure_increase_rate: 0.5,
            pressure_decrease_rate: 0.8,
        }
    }
}
