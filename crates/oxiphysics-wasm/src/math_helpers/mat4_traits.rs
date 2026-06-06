//! # Mat4 - Trait Implementations
//!
//! This module contains trait implementations for `Mat4`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)
use super::types::Mat4;

impl Default for Mat4 {
    fn default() -> Self {
        Self::identity()
    }
}
