//! # SdfCone - Trait Implementations
//!
//! This module contains trait implementations for `SdfCone`.
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
use super::types::SdfCone;

impl Sdf for SdfCone {
    fn dist(&self, p: [f64; 3]) -> f64 {
        let c = [self.half_angle.cos(), self.half_angle.sin()];
        let q = [len2([p[0], p[2]]), p[1]];
        let d = len2([
            q[0] - c[0] * clamp_f(dot2(q, c), 0.0, self.height),
            q[1] - c[1] * clamp_f(dot2(q, c), 0.0, self.height),
        ]);
        let s = if q[0] * c[1] - q[1] * c[0] > 0.0 {
            1.0
        } else {
            -1.0
        };
        let base_d = len2([
            q[0] - clamp_f(q[0], 0.0, self.height * c[1] / c[0]),
            q[1] - self.height,
        ]);
        let cone_d = s * d;
        cone_d.min(base_d)
    }
}
