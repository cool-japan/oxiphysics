//! # SerialWorkQueue - Trait Implementations
//!
//! This module contains trait implementations for `SerialWorkQueue`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::SerialWorkQueue;

impl<T> Default for SerialWorkQueue<T> {
    fn default() -> Self {
        Self::new()
    }
}
