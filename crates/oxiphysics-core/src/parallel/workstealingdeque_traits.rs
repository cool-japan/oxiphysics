//! # WorkStealingDeque - Trait Implementations
//!
//! This module contains trait implementations for `WorkStealingDeque`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::WorkStealingDeque;

impl<T> Default for WorkStealingDeque<T> {
    fn default() -> Self {
        Self::new()
    }
}
