//! # DistanceCV - Trait Implementations
//!
//! This module contains trait implementations for `DistanceCV`.
//!
//! ## Implemented Traits
//!
//! - `CollectiveVariable`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::CollectiveVariable;
use super::types::DistanceCV;

impl CollectiveVariable for DistanceCV {
    fn value(&self, positions: &[[f64; 3]]) -> f64 {
        let ri = positions[self.atom_i];
        let rj = positions[self.atom_j];
        let dx = ri[0] - rj[0];
        let dy = ri[1] - rj[1];
        let dz = ri[2] - rj[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }
    fn gradient(&self, positions: &[[f64; 3]]) -> Vec<[f64; 3]> {
        let n = positions.len();
        let mut grad = vec![[0.0f64; 3]; n];
        let ri = positions[self.atom_i];
        let rj = positions[self.atom_j];
        let dx = ri[0] - rj[0];
        let dy = ri[1] - rj[1];
        let dz = ri[2] - rj[2];
        let dist = (dx * dx + dy * dy + dz * dz).sqrt();
        if dist > 1e-15 {
            let inv_d = 1.0 / dist;
            grad[self.atom_i] = [dx * inv_d, dy * inv_d, dz * inv_d];
            grad[self.atom_j] = [-dx * inv_d, -dy * inv_d, -dz * inv_d];
        }
        grad
    }
    fn name(&self) -> &str {
        &self.name
    }
}
