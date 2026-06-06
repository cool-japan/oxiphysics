//! # SizeDistributionParams - Trait Implementations
//!
//! This module contains trait implementations for `SizeDistributionParams`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::SizeDistributionParams;

impl Default for SizeDistributionParams {
    fn default() -> Self {
        Self {
            mean_radius: 0.05,
            std_radius: 0.01,
            min_radius: 0.01,
            max_radius: 0.1,
            size_ratio: 2.0,
            large_fraction: 0.5,
            power_law_exponent: -3.0,
        }
    }
}
