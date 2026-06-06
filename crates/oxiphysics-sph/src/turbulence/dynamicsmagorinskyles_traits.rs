//! # DynamicSmagorinskyLes - Trait Implementations
//!
//! This module contains trait implementations for `DynamicSmagorinskyLes`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::DynamicSmagorinskyLes;

impl Default for DynamicSmagorinskyLes {
    fn default() -> Self {
        Self {
            delta: 0.1,
            test_filter_ratio: 2.0,
            cs_sq_min: 0.0,
            cs_sq_max: 0.04,
        }
    }
}
