//! # SdfBend - Trait Implementations
//!
//! This module contains trait implementations for `SdfBend`.
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
use super::types::SdfBend;

impl<S: Sdf> Sdf for SdfBend<S> {
    fn dist(&self, p: [f64; 3]) -> f64 {
        let angle = self.strength * p[0];
        let (sin_a, cos_a) = angle.sin_cos();
        let qx = cos_a * p[0] - sin_a * p[1];
        let qy = sin_a * p[0] + cos_a * p[1];
        self.inner.dist([qx, qy, p[2]])
    }
}
