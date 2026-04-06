//! # QualityThresholds - Trait Implementations
//!
//! This module contains trait implementations for `QualityThresholds`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::QualityThresholds;

impl Default for QualityThresholds {
    fn default() -> Self {
        Self {
            max_aspect_ratio: 5.0,
            max_skewness: 0.5,
            min_angle: 15.0,
            max_angle: 150.0,
            min_area: 1e-10,
        }
    }
}
