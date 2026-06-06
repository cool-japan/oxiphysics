//! # WcsphParams - Trait Implementations
//!
//! This module contains trait implementations for `WcsphParams`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::WcsphParams;

impl Default for WcsphParams {
    fn default() -> Self {
        Self {
            rest_density: 1000.0,
            stiffness: 50_000.0,
            gamma: 7.0,
            viscosity: 0.01,
            smoothing_length: 0.1,
        }
    }
}
