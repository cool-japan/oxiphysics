//! # EndToEndDistanceCV - Trait Implementations
//!
//! This module contains trait implementations for `EndToEndDistanceCV`.
//!
//! ## Implemented Traits
//!
//! - `CollectiveVariable`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::CollectiveVariable;
use super::types::EndToEndDistanceCV;

impl CollectiveVariable for EndToEndDistanceCV {
    fn value(&self, positions: &[[f64; 3]]) -> f64 {
        let r0 = positions[self.first];
        let r1 = positions[self.last];
        let dx = r1[0] - r0[0];
        let dy = r1[1] - r0[1];
        let dz = r1[2] - r0[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }
    fn gradient(&self, positions: &[[f64; 3]]) -> Vec<[f64; 3]> {
        let n = positions.len();
        let mut grad = vec![[0.0_f64; 3]; n];
        let r0 = positions[self.first];
        let r1 = positions[self.last];
        let dx = r1[0] - r0[0];
        let dy = r1[1] - r0[1];
        let dz = r1[2] - r0[2];
        let dist = (dx * dx + dy * dy + dz * dz).sqrt();
        if dist > 1e-15 {
            let inv_d = 1.0 / dist;
            grad[self.first] = [-dx * inv_d, -dy * inv_d, -dz * inv_d];
            grad[self.last] = [dx * inv_d, dy * inv_d, dz * inv_d];
        }
        grad
    }
    fn name(&self) -> &str {
        &self.name
    }
}
