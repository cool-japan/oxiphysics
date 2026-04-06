//! # SdfTwist - Trait Implementations
//!
//! This module contains trait implementations for `SdfTwist`.
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
use super::types::SdfTwist;

impl<S: Sdf> Sdf for SdfTwist<S> {
    fn dist(&self, p: [f64; 3]) -> f64 {
        let angle = self.strength * p[1];
        let (sin_a, cos_a) = angle.sin_cos();
        let qx = cos_a * p[0] - sin_a * p[2];
        let qz = sin_a * p[0] + cos_a * p[2];
        self.inner.dist([qx, p[1], qz])
    }
}
