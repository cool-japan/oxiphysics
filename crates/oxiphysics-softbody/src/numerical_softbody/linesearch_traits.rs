//! # LineSearch - Trait Implementations
//!
//! This module contains trait implementations for `LineSearch`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::LineSearch;

impl Default for LineSearch {
    fn default() -> Self {
        LineSearch {
            c1: 1e-4,
            max_iter: 32,
            rho: 0.5,
        }
    }
}
