//! # LesEnergyBudget - Trait Implementations
//!
//! This module contains trait implementations for `LesEnergyBudget`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::LesEnergyBudget;

impl Default for LesEnergyBudget {
    fn default() -> Self {
        Self {
            production: 0.0,
            dissipation: 0.0,
            count: 0,
        }
    }
}
