//! # WaleModel - Trait Implementations
//!
//! This module contains trait implementations for `WaleModel`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::WaleModel;

impl Default for WaleModel {
    fn default() -> Self {
        Self {
            c_w: 0.325,
            delta: 0.1,
        }
    }
}
