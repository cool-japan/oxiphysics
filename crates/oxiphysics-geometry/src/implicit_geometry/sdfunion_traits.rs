//! # SdfUnion - Trait Implementations
//!
//! This module contains trait implementations for `SdfUnion`.
//!
//! ## Implemented Traits
//!
//! - `Sdf`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
use super::types::SdfUnion;

impl<A: Sdf, B: Sdf> Sdf for SdfUnion<A, B> {
    fn dist(&self, p: [f64; 3]) -> f64 {
        self.a.dist(p).min(self.b.dist(p))
    }
}
