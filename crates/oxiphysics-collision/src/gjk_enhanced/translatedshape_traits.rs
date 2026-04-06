//! # TranslatedShape - Trait Implementations
//!
//! This module contains trait implementations for `TranslatedShape`.
//!
//! ## Implemented Traits
//!
//! - `ConvexShape`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::ConvexShape;
#[allow(unused_imports)]
use super::functions::*;
use super::types::TranslatedShape;

impl ConvexShape for TranslatedShape<'_> {
    fn support(&self, dir: [f64; 3]) -> [f64; 3] {
        vadd(self.shape.support(dir), self.offset)
    }
    fn centre(&self) -> [f64; 3] {
        vadd(self.shape.centre(), self.offset)
    }
}
