//! # SdfPyramid - Trait Implementations
//!
//! This module contains trait implementations for `SdfPyramid`.
//!
//! ## Implemented Traits
//!
//! - `Sdf`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
use super::types::SdfPyramid;

impl Sdf for SdfPyramid {
    fn dist(&self, p: [f64; 3]) -> f64 {
        let m2 = self.half_base * self.half_base + self.height * self.height;
        let px = p[0].abs();
        let pz = p[2].abs();
        let (px, pz) = if pz > px { (pz, px) } else { (px, pz) };
        let px = px - self.half_base;
        let py = p[1] - self.height;
        let k = clamp_f(-py * self.half_base - px * self.height, 0.0, m2);
        let dx = px - (-self.half_base * k / m2);
        let dy = py - self.height * k / m2;
        let d1 = ((px * self.height - py * (-self.half_base)).max(0.0).powi(2)
            + (pz + self.half_base).max(0.0).powi(2))
        .sqrt();
        let d2 = len2([dx, dy]);
        let s = if px * (-self.height) - py * (-self.half_base) > 0.0 {
            -1.0
        } else {
            1.0
        };
        _ = d1;
        s * d2.min(d1)
    }
}
