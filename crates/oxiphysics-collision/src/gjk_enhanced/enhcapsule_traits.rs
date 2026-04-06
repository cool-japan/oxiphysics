//! # EnhCapsule - Trait Implementations
//!
//! This module contains trait implementations for `EnhCapsule`.
//!
//! ## Implemented Traits
//!
//! - `ConvexShape`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::ConvexShape;
#[allow(unused_imports)]
use super::functions::*;
use super::types::EnhCapsule;

impl ConvexShape for EnhCapsule {
    fn support(&self, dir: [f64; 3]) -> [f64; 3] {
        let d = vnorm(dir);
        let dot_a = vdot(d, self.a);
        let dot_b = vdot(d, self.b);
        let tip = if dot_a >= dot_b { self.a } else { self.b };
        vadd(tip, vscale(d, self.radius))
    }
    fn centre(&self) -> [f64; 3] {
        [
            (self.a[0] + self.b[0]) * 0.5,
            (self.a[1] + self.b[1]) * 0.5,
            (self.a[2] + self.b[2]) * 0.5,
        ]
    }
}
