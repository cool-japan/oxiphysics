//! # DynamicSmagorinsky - Trait Implementations
//!
//! This module contains trait implementations for `DynamicSmagorinsky`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::DynamicSmagorinsky;

impl Default for DynamicSmagorinsky {
    fn default() -> Self {
        Self {
            delta: 0.1,
            test_filter_ratio: 2.0,
            c_min: 0.0,
            c_max: 0.25,
        }
    }
}
