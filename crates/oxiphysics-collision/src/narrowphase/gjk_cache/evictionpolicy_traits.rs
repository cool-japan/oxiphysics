//! # EvictionPolicy - Trait Implementations
//!
//! This module contains trait implementations for `EvictionPolicy`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::EvictionPolicy;

impl Default for EvictionPolicy {
    fn default() -> Self {
        Self::conservative()
    }
}
