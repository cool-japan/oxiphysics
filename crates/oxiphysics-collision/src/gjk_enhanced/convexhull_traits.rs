//! # ConvexHull - Trait Implementations
//!
//! This module contains trait implementations for `ConvexHull`.
//!
//! ## Implemented Traits
//!
//! - `ConvexShape`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::{ConvexShape, vadd, vdot, vscale};
use super::types::ConvexHull;

impl ConvexShape for ConvexHull {
    fn support(&self, dir: [f64; 3]) -> [f64; 3] {
        self.vertices
            .iter()
            .cloned()
            .max_by(|&a, &b| {
                vdot(a, dir)
                    .partial_cmp(&vdot(b, dir))
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .unwrap_or([0.0; 3])
    }
    fn centre(&self) -> [f64; 3] {
        if self.vertices.is_empty() {
            return [0.0; 3];
        }
        let n = self.vertices.len() as f64;
        let sum = self.vertices.iter().fold([0.0; 3], |acc, &v| vadd(acc, v));
        vscale(sum, 1.0 / n)
    }
}
