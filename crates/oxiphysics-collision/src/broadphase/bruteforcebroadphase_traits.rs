//! # BruteForceBroadPhase - Trait Implementations
//!
//! This module contains trait implementations for `BruteForceBroadPhase`.
//!
//! ## Implemented Traits
//!
//! - `BroadPhase`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::types::CollisionPair;
use oxiphysics_core::Aabb;

use super::functions::BroadPhase;
use super::types::BruteForceBroadPhase;

impl BroadPhase for BruteForceBroadPhase {
    fn find_pairs(&self, aabbs: &[Aabb]) -> Vec<CollisionPair> {
        let mut pairs = Vec::new();
        for i in 0..aabbs.len() {
            for j in (i + 1)..aabbs.len() {
                if aabbs[i].intersects(&aabbs[j]) {
                    pairs.push(CollisionPair::new(i, j));
                }
            }
        }
        pairs
    }
}
