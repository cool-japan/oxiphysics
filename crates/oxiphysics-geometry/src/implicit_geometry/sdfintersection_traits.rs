//! # SdfIntersection - Trait Implementations
//!
//! This module contains trait implementations for `SdfIntersection`.
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
use super::types::SdfIntersection;

impl<A: Sdf, B: Sdf> Sdf for SdfIntersection<A, B> {
    fn dist(&self, p: [f64; 3]) -> f64 {
        self.a.dist(p).max(self.b.dist(p))
    }
}
