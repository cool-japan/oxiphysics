//! # PcgParams - Trait Implementations
//!
//! This module contains trait implementations for `PcgParams`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::PcgParams;

impl Default for PcgParams {
    fn default() -> Self {
        PcgParams {
            max_iter: 1000,
            tolerance: 1e-8,
        }
    }
}
