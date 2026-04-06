//! # KOmegaRans - Trait Implementations
//!
//! This module contains trait implementations for `KOmegaRans`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::KOmegaRans;

impl Default for KOmegaRans {
    fn default() -> Self {
        Self {
            alpha: 5.0 / 9.0,
            beta: 3.0 / 40.0,
            beta_star: 9.0 / 100.0,
            sigma_k: 0.5,
            sigma_omega: 0.5,
            nu: 1e-6,
        }
    }
}
