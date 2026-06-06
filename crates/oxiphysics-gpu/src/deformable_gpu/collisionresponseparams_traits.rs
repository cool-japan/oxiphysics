//! # CollisionResponseParams - Trait Implementations
//!
//! This module contains trait implementations for `CollisionResponseParams`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)
use super::types::CollisionResponseParams;

impl Default for CollisionResponseParams {
    fn default() -> Self {
        Self::new(1e5, 100.0, 1e-3)
    }
}
