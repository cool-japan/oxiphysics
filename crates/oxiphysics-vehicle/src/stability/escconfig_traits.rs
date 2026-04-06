//! # EscConfig - Trait Implementations
//!
//! This module contains trait implementations for `EscConfig`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::EscConfig;

impl Default for EscConfig {
    fn default() -> Self {
        Self {
            yaw_rate_threshold: 0.05,
            sideslip_threshold: 5.0,
            understeer_gain: 0.6,
            oversteer_gain: 0.8,
        }
    }
}
