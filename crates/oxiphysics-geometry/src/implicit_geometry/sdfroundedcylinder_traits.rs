//! # SdfRoundedCylinder - Trait Implementations
//!
//! This module contains trait implementations for `SdfRoundedCylinder`.
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
use super::types::SdfRoundedCylinder;

impl Sdf for SdfRoundedCylinder {
    fn dist(&self, p: [f64; 3]) -> f64 {
        let r = self.radius - self.rounding;
        let h = self.half_height - self.rounding;
        let d0 = len2([p[0], p[2]]) - r;
        let d1 = p[1].abs() - h;
        let d_in = d0.max(d1).min(0.0);
        let d_out = len2([d0.max(0.0), d1.max(0.0)]);
        d_in + d_out - self.rounding
    }
}
