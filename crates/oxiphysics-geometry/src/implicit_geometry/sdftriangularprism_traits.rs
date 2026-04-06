//! # SdfTriangularPrism - Trait Implementations
//!
//! This module contains trait implementations for `SdfTriangularPrism`.
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
use super::types::SdfTriangularPrism;

impl Sdf for SdfTriangularPrism {
    fn dist(&self, p: [f64; 3]) -> f64 {
        let h = [self.side, self.half_height];
        let q = [p[0].abs(), p[1], p[2].abs()];
        let k = 3.0_f64.sqrt();
        let qx = q[0] - clamp_f(q[0] - k * q[2], 0.0, 1.0) * h[0] / 2.0;
        let qz = q[2] - clamp_f(q[0] - k * q[2], 0.0, 1.0) * h[0] / (2.0 * k);
        let d1 = len2([qx, qz])
            * if q[0] - h[0] / 2.0 < 0.0 && k * q[2] - q[0] < 0.0 {
                -1.0
            } else {
                1.0
            };
        let d2 = q[1] - h[1];
        -((-d1).max(-d2)).min(0.0) + len2([d1.max(0.0), d2.max(0.0)])
    }
}
