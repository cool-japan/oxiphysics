//! # EbdConfig - Trait Implementations
//!
//! This module contains trait implementations for `EbdConfig`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)
use super::types::EbdConfig;

impl Default for EbdConfig {
    fn default() -> Self {
        Self {
            front_bias: 0.65,
            decel_bias_gain: 0.03,
            min_rear_fraction: 0.1,
            max_front_fraction: 0.85,
        }
    }
}
