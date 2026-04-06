//! # SdfSmoothDifference - Trait Implementations
//!
//! This module contains trait implementations for `SdfSmoothDifference`.
//!
//! ## Implemented Traits
//!
//! - `Sdf`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
#[allow(unused_imports)]
use super::functions::*;
use super::types::SdfSmoothDifference;

impl<A: Sdf, B: Sdf> Sdf for SdfSmoothDifference<A, B> {
    fn dist(&self, p: [f64; 3]) -> f64 {
        sdf_smooth_difference(self.a.dist(p), self.b.dist(p), self.k)
    }
}
