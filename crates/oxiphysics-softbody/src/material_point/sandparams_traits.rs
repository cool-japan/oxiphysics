//! # SandParams - Trait Implementations
//!
//! This module contains trait implementations for `SandParams`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::SandParams;

impl Default for SandParams {
    fn default() -> Self {
        Self {
            mu0: 1e5,
            lambda0: 2e5,
            friction_angle: 30.0,
            rho0: 1600.0,
        }
    }
}
