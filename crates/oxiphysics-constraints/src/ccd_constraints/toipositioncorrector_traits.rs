//! # ToiPositionCorrector - Trait Implementations
//!
//! This module contains trait implementations for `ToiPositionCorrector`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::ToiPositionCorrector;

impl Default for ToiPositionCorrector {
    fn default() -> Self {
        ToiPositionCorrector {
            baumgarte: 0.2,
            slop: 0.005,
            max_correction: 0.2,
        }
    }
}
