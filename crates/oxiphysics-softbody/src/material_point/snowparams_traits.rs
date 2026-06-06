//! # SnowParams - Trait Implementations
//!
//! This module contains trait implementations for `SnowParams`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::SnowParams;

impl Default for SnowParams {
    fn default() -> Self {
        Self {
            mu0: 58333.0,
            lambda0: 38194.0,
            hardening: 10.0,
            theta_c: 0.025,
            theta_s: 0.0075,
            rho0: 400.0,
        }
    }
}
