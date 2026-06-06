//! # Differential - Trait Implementations
//!
//! This module contains trait implementations for `Differential`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)
use super::types::Differential;

impl Default for Differential {
    fn default() -> Self {
        Self::open()
    }
}
