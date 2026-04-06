//! # SdfSphere - Trait Implementations
//!
//! This module contains trait implementations for `SdfSphere`.
//!
//! ## Implemented Traits
//!
//! - `ImplicitSurface`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::ImplicitSurface;
#[allow(unused_imports)]
use super::functions::*;
use super::functions::{length, normalize, sub};
use super::types::SdfSphere;

impl ImplicitSurface for SdfSphere {
    fn sdf(&self, p: [f64; 3]) -> f64 {
        length(sub(p, self.center)) - self.radius
    }
    fn gradient(&self, p: [f64; 3]) -> [f64; 3] {
        normalize(sub(p, self.center))
    }
}
