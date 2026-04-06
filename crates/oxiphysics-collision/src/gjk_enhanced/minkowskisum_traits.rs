//! # MinkowskiSum - Trait Implementations
//!
//! This module contains trait implementations for `MinkowskiSum`.
//!
//! ## Implemented Traits
//!
//! - `ConvexShape`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::ConvexShape;
#[allow(unused_imports)]
use super::functions::*;
use super::types::MinkowskiSum;

impl ConvexShape for MinkowskiSum<'_> {
    fn support(&self, dir: [f64; 3]) -> [f64; 3] {
        vadd(self.a.support(dir), self.b.support(dir))
    }
    fn centre(&self) -> [f64; 3] {
        vadd(self.a.centre(), self.b.centre())
    }
}
