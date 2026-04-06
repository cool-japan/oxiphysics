//! # IslandManager - Trait Implementations
//!
//! This module contains trait implementations for `IslandManager`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::IslandManager;

impl Default for IslandManager {
    fn default() -> Self {
        Self::new(0.01, 0.01, 0.5)
    }
}
