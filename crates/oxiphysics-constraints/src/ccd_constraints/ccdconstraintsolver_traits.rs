//! # CcdConstraintSolver - Trait Implementations
//!
//! This module contains trait implementations for `CcdConstraintSolver`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use oxiphysics_core::math::Vec3;

use super::types::CcdConstraintSolver;

impl Default for CcdConstraintSolver {
    fn default() -> Self {
        Self::new(4, Vec3::new(0.0, -9.81, 0.0))
    }
}
