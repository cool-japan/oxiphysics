//! # GearboxLegacy - Trait Implementations
//!
//! This module contains trait implementations for `GearboxLegacy`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types_ext::GearboxLegacy;

impl Default for GearboxLegacy {
    fn default() -> Self {
        Self {
            gear_ratios: vec![3.5, 2.5, 1.8, 1.3, 1.0, 0.8],
            reverse_ratio: 3.2,
            final_drive_ratio: 3.7,
            current_gear: 1,
            efficiency: 0.9,
        }
    }
}
