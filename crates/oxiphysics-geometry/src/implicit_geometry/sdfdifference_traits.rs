//! # SdfDifference - Trait Implementations
//!
//! This module contains trait implementations for `SdfDifference`.
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
use super::types::SdfDifference;

impl<A: Sdf, B: Sdf> Sdf for SdfDifference<A, B> {
    fn dist(&self, p: [f64; 3]) -> f64 {
        self.a.dist(p).max(-self.b.dist(p))
    }
}
