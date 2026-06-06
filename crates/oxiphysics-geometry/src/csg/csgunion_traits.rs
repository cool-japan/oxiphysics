//! # CsgUnion - Trait Implementations
//!
//! This module contains trait implementations for `CsgUnion`.
//!
//! ## Implemented Traits
//!
//! - `ImplicitSurface`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::ImplicitSurface;
use super::types::CsgUnion;

impl ImplicitSurface for CsgUnion {
    fn sdf(&self, p: [f64; 3]) -> f64 {
        self.a.sdf(p).min(self.b.sdf(p))
    }
    fn gradient(&self, p: [f64; 3]) -> [f64; 3] {
        if self.a.sdf(p) < self.b.sdf(p) {
            self.a.gradient(p)
        } else {
            self.b.gradient(p)
        }
    }
}
