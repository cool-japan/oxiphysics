//! # SupportCache - Trait Implementations
//!
//! This module contains trait implementations for `SupportCache`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::SupportCache;

impl<const CAP: usize> Default for SupportCache<CAP> {
    fn default() -> Self {
        Self::new()
    }
}
