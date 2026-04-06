//! # SdfSmoothIntersection - Trait Implementations
//!
//! This module contains trait implementations for `SdfSmoothIntersection`.
//!
//! ## Implemented Traits
//!
//! - `Sdf`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
use super::types::SdfSmoothIntersection;

impl<A: Sdf, B: Sdf> Sdf for SdfSmoothIntersection<A, B> {
    fn dist(&self, p: [f64; 3]) -> f64 {
        sdf_smooth_intersection(self.a.dist(p), self.b.dist(p), self.k)
    }
}
