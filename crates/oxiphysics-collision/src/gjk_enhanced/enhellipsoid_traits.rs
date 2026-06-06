//! # EnhEllipsoid - Trait Implementations
//!
//! This module contains trait implementations for `EnhEllipsoid`.
//!
//! ## Implemented Traits
//!
//! - `ConvexShape`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::{ConvexShape, vadd, vdot, vlen, vscale};
use super::types::EnhEllipsoid;

impl ConvexShape for EnhEllipsoid {
    fn support(&self, dir: [f64; 3]) -> [f64; 3] {
        let a = self.semi_axes;
        let scaled = [
            a[0] * a[0] * dir[0],
            a[1] * a[1] * dir[1],
            a[2] * a[2] * dir[2],
        ];
        let l = vlen(scaled).max(1e-300);
        vadd(
            self.centre,
            vscale(scaled, 1.0 / l * (vdot(scaled, dir) / l).sqrt().max(1.0)),
        )
    }
    fn centre(&self) -> [f64; 3] {
        self.centre
    }
}
