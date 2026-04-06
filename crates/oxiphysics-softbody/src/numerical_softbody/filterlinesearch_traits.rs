//! # FilterLineSearch - Trait Implementations
//!
//! This module contains trait implementations for `FilterLineSearch`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::FilterLineSearch;

impl Default for FilterLineSearch {
    fn default() -> Self {
        FilterLineSearch {
            tau_max: 0.8,
            c1: 1e-4,
            rho: 0.5,
            max_iter: 32,
        }
    }
}
