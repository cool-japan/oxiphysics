//! # EnergyConvergenceCriteria - Trait Implementations
//!
//! This module contains trait implementations for `EnergyConvergenceCriteria`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::EnergyConvergenceCriteria;

impl Default for EnergyConvergenceCriteria {
    fn default() -> Self {
        Self::new(50, 1e-6)
    }
}
