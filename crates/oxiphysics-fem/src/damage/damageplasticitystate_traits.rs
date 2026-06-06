//! # DamagePlasticityState - Trait Implementations
//!
//! This module contains trait implementations for `DamagePlasticityState`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::DamagePlasticityState;

impl Default for DamagePlasticityState {
    fn default() -> Self {
        Self {
            d: 0.0,
            p_bar: 0.0,
            r: 0.0,
            back_stress: [0.0; 6],
            effective_stress: [0.0; 6],
        }
    }
}
