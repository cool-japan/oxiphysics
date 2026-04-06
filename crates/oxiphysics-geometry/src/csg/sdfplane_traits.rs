//! # SdfPlane - Trait Implementations
//!
//! This module contains trait implementations for `SdfPlane`.
//!
//! ## Implemented Traits
//!
//! - `ImplicitSurface`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::ImplicitSurface;
use super::functions::dot;
#[allow(unused_imports)]
use super::functions::*;
use super::types::SdfPlane;

impl ImplicitSurface for SdfPlane {
    fn sdf(&self, p: [f64; 3]) -> f64 {
        dot(p, self.normal) - self.d
    }
    fn gradient(&self, _p: [f64; 3]) -> [f64; 3] {
        self.normal
    }
}
