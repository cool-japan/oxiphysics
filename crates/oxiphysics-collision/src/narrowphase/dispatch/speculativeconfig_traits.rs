//! # SpeculativeConfig - Trait Implementations
//!
//! This module contains trait implementations for `SpeculativeConfig`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::SpeculativeConfig;

impl Default for SpeculativeConfig {
    fn default() -> Self {
        SpeculativeConfig {
            margin: 0.02,
            velocity_scale: 10.0,
        }
    }
}
