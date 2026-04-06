//! # GripperState - Trait Implementations
//!
//! This module contains trait implementations for `GripperState`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::GripperState;

impl Default for GripperState {
    fn default() -> Self {
        Self {
            position: 0.0,
            force: 0.0,
            max_width: 0.085,
            object_detected: false,
            temperature: 25.0,
        }
    }
}
