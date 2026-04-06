//! # ConvergenceCriteria - Trait Implementations
//!
//! This module contains trait implementations for `ConvergenceCriteria`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::ConvergenceCriteria;

impl Default for ConvergenceCriteria {
    fn default() -> Self {
        Self {
            max_iter: 50,
            residual_tol: 1e-8,
            displacement_tol: 1e-10,
        }
    }
}
