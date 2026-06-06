//! # BlastParameters - Trait Implementations
//!
//! This module contains trait implementations for `BlastParameters`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::BlastParameters;

impl Default for BlastParameters {
    fn default() -> Self {
        Self {
            peak_pressure: 1e6,
            duration: 5e-3,
            decay_coeff: 1.0,
            standoff: 1.0,
            charge_mass: 1.0,
            table: Vec::new(),
        }
    }
}
