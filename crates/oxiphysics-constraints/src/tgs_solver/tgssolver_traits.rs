//! # TgsSolver - Trait Implementations
//!
//! This module contains trait implementations for `TgsSolver`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use oxiphysics_core::math::Vec3;

use super::types::TgsSolver;

impl Default for TgsSolver {
    fn default() -> Self {
        Self::new(4, 4, Vec3::new(0.0, -9.81, 0.0))
    }
}
