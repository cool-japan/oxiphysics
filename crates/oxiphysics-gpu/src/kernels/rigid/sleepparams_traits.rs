//! # SleepParams - Trait Implementations
//!
//! This module contains trait implementations for `SleepParams`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::SleepParams;

impl Default for SleepParams {
    fn default() -> Self {
        Self {
            linear_threshold: 0.01,
            angular_threshold: 0.01,
            sleep_frames: 10,
        }
    }
}
