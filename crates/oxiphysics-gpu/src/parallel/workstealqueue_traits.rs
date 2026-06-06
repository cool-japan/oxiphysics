//! # WorkStealQueue - Trait Implementations
//!
//! This module contains trait implementations for `WorkStealQueue`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::WorkStealQueue;

impl<T: Send> Default for WorkStealQueue<T> {
    fn default() -> Self {
        Self::new()
    }
}
