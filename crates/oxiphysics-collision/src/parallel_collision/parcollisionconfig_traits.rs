//! # ParCollisionConfig - Trait Implementations
//!
//! This module contains trait implementations for `ParCollisionConfig`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::ParCollisionConfig;

impl Default for ParCollisionConfig {
    fn default() -> Self {
        Self {
            max_contacts_per_pair: 4,
            aabb_margin: 0.02,
            enable_ccd: false,
        }
    }
}
