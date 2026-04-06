//! # RsaConfig - Trait Implementations
//!
//! This module contains trait implementations for `RsaConfig`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::RsaConfig;

impl Default for RsaConfig {
    fn default() -> Self {
        Self {
            radius: 0.05,
            domain: [1.0, 1.0, 1.0],
            max_attempts: 10_000,
            min_gap: 0.0,
        }
    }
}
