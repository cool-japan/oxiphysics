//! # SdfHexagonalPrism - Trait Implementations
//!
//! This module contains trait implementations for `SdfHexagonalPrism`.
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
use super::types::SdfHexagonalPrism;

impl Sdf for SdfHexagonalPrism {
    fn dist(&self, p: [f64; 3]) -> f64 {
        let k0 = -3.0_f64.sqrt() / 2.0;
        let k1 = 0.5_f64;
        let k2 = 1.0_f64 / 3.0_f64.sqrt();
        let qa = [p[0].abs(), p[2].abs()];
        let dot_val = 2.0_f64.min(0.0_f64.max(k0 * qa[0] + k1 * qa[1]));
        let qx = qa[0] - 2.0 * k0 * dot_val;
        let qz = qa[1] - 2.0 * k1 * dot_val;
        let ex = qx - clamp_f(qx, -k2 * self.radius, k2 * self.radius);
        let ez = qz - self.radius;
        let d_xz = len2([ex, ez]) * if qz > self.radius { 1.0 } else { -1.0 };
        let d_y = p[1].abs() - self.half_height;
        let d_xz_pos = d_xz.max(0.0);
        let d_y_pos = d_y.max(0.0);
        len2([d_xz_pos, d_y_pos]) + d_xz.min(0.0).max(d_y.min(0.0))
    }
}
