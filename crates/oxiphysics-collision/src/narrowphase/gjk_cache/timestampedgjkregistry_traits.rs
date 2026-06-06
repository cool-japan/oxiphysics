//! # TimestampedGjkRegistry - Trait Implementations
//!
//! This module contains trait implementations for `TimestampedGjkRegistry`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{EvictionPolicy, TimestampedGjkRegistry};

impl Default for TimestampedGjkRegistry {
    fn default() -> Self {
        Self::new(EvictionPolicy::default())
    }
}
