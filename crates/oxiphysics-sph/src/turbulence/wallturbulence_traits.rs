//! # WallTurbulence - Trait Implementations
//!
//! This module contains trait implementations for `WallTurbulence`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::WallTurbulence;

impl Default for WallTurbulence {
    fn default() -> Self {
        Self {
            kappa: 0.41,
            a_plus: 26.0,
            nu: 1e-6,
        }
    }
}
