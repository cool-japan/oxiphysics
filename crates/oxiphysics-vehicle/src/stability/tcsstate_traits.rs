//! # TcsState - Trait Implementations
//!
//! This module contains trait implementations for `TcsState`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::TcsState;

impl Default for TcsState {
    fn default() -> Self {
        Self {
            is_active: false,
            torque_reduction_factor: 1.0,
            slip_ratio: 0.0,
        }
    }
}
