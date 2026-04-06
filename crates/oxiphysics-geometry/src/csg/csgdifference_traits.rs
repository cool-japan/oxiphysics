//! # CsgDifference - Trait Implementations
//!
//! This module contains trait implementations for `CsgDifference`.
//!
//! ## Implemented Traits
//!
//! - `ImplicitSurface`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::ImplicitSurface;
#[allow(unused_imports)]
use super::functions::*;
use super::types::CsgDifference;

impl ImplicitSurface for CsgDifference {
    fn sdf(&self, p: [f64; 3]) -> f64 {
        self.a.sdf(p).max(-self.b.sdf(p))
    }
    fn gradient(&self, p: [f64; 3]) -> [f64; 3] {
        if self.a.sdf(p) > -self.b.sdf(p) {
            self.a.gradient(p)
        } else {
            let g = self.b.gradient(p);
            [-g[0], -g[1], -g[2]]
        }
    }
}
