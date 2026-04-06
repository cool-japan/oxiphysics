//! # ExplicitFemConfig - Trait Implementations
//!
//! This module contains trait implementations for `ExplicitFemConfig`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::{ExplicitFemConfig, MassLumping};

impl Default for ExplicitFemConfig {
    fn default() -> Self {
        Self {
            cfl_factor: 0.9,
            lumping: MassLumping::RowSum,
            wave_speed: 5000.0,
            min_element_length: 1e-3,
            bulk_viscosity: 0.06,
            check_energy: true,
        }
    }
}
