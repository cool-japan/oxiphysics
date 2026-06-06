//! # Color4 - Trait Implementations
//!
//! This module contains trait implementations for `Color4`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::Color4;

impl Default for Color4 {
    /// Default color is opaque white.
    fn default() -> Self {
        Self::white()
    }
}
