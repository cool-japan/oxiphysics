//! # CsgIntersection - Trait Implementations
//!
//! This module contains trait implementations for `CsgIntersection`.
//!
//! ## Implemented Traits
//!
//! - `ImplicitSurface`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::ImplicitSurface;
#[allow(unused_imports)]
use super::functions::*;
use super::types::CsgIntersection;

impl ImplicitSurface for CsgIntersection {
    fn sdf(&self, p: [f64; 3]) -> f64 {
        self.a.sdf(p).max(self.b.sdf(p))
    }
    fn gradient(&self, p: [f64; 3]) -> [f64; 3] {
        if self.a.sdf(p) > self.b.sdf(p) {
            self.a.gradient(p)
        } else {
            self.b.gradient(p)
        }
    }
}
