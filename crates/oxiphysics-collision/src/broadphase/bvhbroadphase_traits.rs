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
    /// Overlapping pairs of `aabbs`, found by a pruned self-traversal of the
    /// hierarchy.
    ///
    /// Also correct on a hierarchy that was never built (for example one
    /// created with [`Default`]) or that no longer matches `aabbs`: a
    /// temporary hierarchy is built from `aabbs` in that case, see
    /// [`BvhBroadphase::find_all_pairs`].
    fn find_pairs(&self, aabbs: &[Aabb]) -> Vec<CollisionPair> {
        self.find_all_pairs(aabbs)
    }
}
