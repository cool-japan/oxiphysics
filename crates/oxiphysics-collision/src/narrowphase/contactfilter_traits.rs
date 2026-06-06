//! # ContactFilter - Trait Implementations
//!
//! This module contains trait implementations for `ContactFilter`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

pub use super::specialized::*;

use super::types::ContactFilter;

impl Default for ContactFilter {
    fn default() -> Self {
        Self {
            min_depth: 0.0,
            max_depth: f64::INFINITY,
            flip_normal: false,
        }
    }
}
