//! # SigmaModel - Trait Implementations
//!
//! This module contains trait implementations for `SigmaModel`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::SigmaModel;

impl Default for SigmaModel {
    fn default() -> Self {
        Self {
            c_sigma: 1.35,
            delta: 0.1,
        }
    }
}
