//! # DamageState - Trait Implementations
//!
//! This module contains trait implementations for `DamageState`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::DamageState;

impl Default for DamageState {
    fn default() -> Self {
        Self {
            d: 0.0,
            kappa: 0.0,
            plastic_strain: [0.0; 6],
            accumulated_plastic_strain: 0.0,
            damage_rate: 0.0,
            deleted: false,
        }
    }
}
