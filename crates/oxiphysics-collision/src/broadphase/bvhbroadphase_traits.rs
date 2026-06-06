//! # BvhBroadphase - Trait Implementations
//!
//! This module contains trait implementations for `BvhBroadphase`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//! - `BroadPhase`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::types::CollisionPair;
use oxiphysics_core::Aabb;

use super::functions::BroadPhase;
use super::types::BvhBroadphase;

impl Default for BvhBroadphase {
    fn default() -> Self {
        Self::new()
    }
}

impl BroadPhase for BvhBroadphase {
    fn find_pairs(&self, aabbs: &[Aabb]) -> Vec<CollisionPair> {
        self.find_all_pairs(aabbs)
    }
}
