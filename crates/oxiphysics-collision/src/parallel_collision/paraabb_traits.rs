//! # ParAabb - Trait Implementations
//!
//! This module contains trait implementations for `ParAabb`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::ParAabb;

impl Default for ParAabb {
    fn default() -> Self {
        Self {
            min: [f64::MAX; 3],
            max: [f64::MIN; 3],
        }
    }
}
