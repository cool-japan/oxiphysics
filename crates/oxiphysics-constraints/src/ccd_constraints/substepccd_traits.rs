//! # SubStepCcd - Trait Implementations
//!
//! This module contains trait implementations for `SubStepCcd`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::SubStepCcd;

impl Default for SubStepCcd {
    fn default() -> Self {
        Self::new(8, 4)
    }
}
