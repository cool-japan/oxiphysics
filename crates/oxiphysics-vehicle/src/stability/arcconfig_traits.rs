//! # ArcConfig - Trait Implementations
//!
//! This module contains trait implementations for `ArcConfig`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::ArcConfig;

impl Default for ArcConfig {
    fn default() -> Self {
        Self {
            roll_gain: 5000.0,
            roll_rate_gain: 1200.0,
            max_torque: 2000.0,
            front_fraction: 0.6,
        }
    }
}
