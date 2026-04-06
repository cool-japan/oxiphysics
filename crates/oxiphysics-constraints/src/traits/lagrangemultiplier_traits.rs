//! # LagrangeMultiplier - Trait Implementations
//!
//! This module contains trait implementations for `LagrangeMultiplier`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::LagrangeMultiplier;

impl Default for LagrangeMultiplier {
    fn default() -> Self {
        Self::bilateral()
    }
}
