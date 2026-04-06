//! # MixingLength - Trait Implementations
//!
//! This module contains trait implementations for `MixingLength`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::MixingLength;

impl Default for MixingLength {
    fn default() -> Self {
        Self {
            wall_distance: 0.0,
            von_karman: 0.41,
            viscous_length: 1.0e-4,
        }
    }
}
