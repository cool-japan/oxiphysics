//! # ConstraintParams - Trait Implementations
//!
//! This module contains trait implementations for `ConstraintParams`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::ConstraintParams;

impl Default for ConstraintParams {
    fn default() -> Self {
        Self {
            iterations: 10,
            relaxation: 0.8,
            friction: 0.3,
        }
    }
}
