//! # MixingLengthModel - Trait Implementations
//!
//! This module contains trait implementations for `MixingLengthModel`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::MixingLengthModel;

impl Default for MixingLengthModel {
    fn default() -> Self {
        Self {
            kappa: 0.41,
            a_plus: 26.0,
            l_max: f64::MAX,
        }
    }
}
