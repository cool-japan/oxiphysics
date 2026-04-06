//! # DeformCollisionConfig - Trait Implementations
//!
//! This module contains trait implementations for `DeformCollisionConfig`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::DeformCollisionConfig;

impl Default for DeformCollisionConfig {
    fn default() -> Self {
        Self {
            thickness: 0.02,
            edge_edge: true,
            self_collision: true,
            adjacency_skip: true,
        }
    }
}
