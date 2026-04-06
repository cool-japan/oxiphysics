//! # SweepAndPrune - Trait Implementations
//!
//! This module contains trait implementations for `SweepAndPrune`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//! - `BroadPhase`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::types::CollisionPair;
use oxiphysics_core::Aabb;
use oxiphysics_core::math::Real;

use super::functions::BroadPhase;
#[allow(unused_imports)]
use super::functions::*;
use super::types::SweepAndPrune;

impl Default for SweepAndPrune {
    fn default() -> Self {
        Self::x_axis()
    }
}

impl BroadPhase for SweepAndPrune {
    fn find_pairs(&self, aabbs: &[Aabb]) -> Vec<CollisionPair> {
        let ax = self.axis;
        let mut sorted: Vec<(Real, usize)> = aabbs
            .iter()
            .enumerate()
            .map(|(i, aabb)| (aabb.min[ax], i))
            .collect();
        sorted.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        let mut pairs = Vec::new();
        for i in 0..sorted.len() {
            let (_, idx_a) = sorted[i];
            let max_a = aabbs[idx_a].max[ax];
            for &(min_b, idx_b) in &sorted[(i + 1)..] {
                if min_b > max_a {
                    break;
                }
                if aabbs[idx_a].intersects(&aabbs[idx_b]) {
                    let (a, b) = if idx_a < idx_b {
                        (idx_a, idx_b)
                    } else {
                        (idx_b, idx_a)
                    };
                    pairs.push(CollisionPair::new(a, b));
                }
            }
        }
        pairs
    }
}
