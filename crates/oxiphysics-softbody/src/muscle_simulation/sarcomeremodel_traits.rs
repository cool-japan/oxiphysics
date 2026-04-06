//! # SarcomereModel - Trait Implementations
//!
//! This module contains trait implementations for `SarcomereModel`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::SarcomereModel;

impl Default for SarcomereModel {
    fn default() -> Self {
        Self {
            attached_fraction: 0.0,
            sarcomere_length: 2.2e-6,
            optimal_sarcomere_length: 2.2e-6,
            attach_rate: 40.0,
            detach_rate_shortening: 15.0,
            detach_rate_lengthening: 200.0,
            max_isometric_stress: 250_000.0,
            length_plateau_width: 0.3,
        }
    }
}
