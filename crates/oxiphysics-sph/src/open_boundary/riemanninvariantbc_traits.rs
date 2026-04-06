//! # RiemannInvariantBc - Trait Implementations
//!
//! This module contains trait implementations for `RiemannInvariantBc`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::RiemannInvariantBc;

impl Default for RiemannInvariantBc {
    fn default() -> Self {
        Self {
            gamma: 1.4,
            c_ref: 340.0,
            rho_ref: 1.225,
        }
    }
}
