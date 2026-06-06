//! # KepsilonSph - Trait Implementations
//!
//! This module contains trait implementations for `KepsilonSph`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::KepsilonSph;

impl Default for KepsilonSph {
    fn default() -> Self {
        Self {
            k: 1e-4,
            epsilon: 1e-6,
            c_mu: 0.09,
            c1_eps: 1.44,
            c2_eps: 1.92,
            sigma_k: 1.0,
            sigma_eps: 1.3,
        }
    }
}
