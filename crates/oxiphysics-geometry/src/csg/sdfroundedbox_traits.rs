//! # SdfRoundedBox - Trait Implementations
//!
//! This module contains trait implementations for `SdfRoundedBox`.
//!
//! ## Implemented Traits
//!
//! - `ImplicitSurface`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::ImplicitSurface;
#[allow(unused_imports)]
use super::functions::*;
use super::functions::{length, normalize};
use super::types::SdfRoundedBox;

impl ImplicitSurface for SdfRoundedBox {
    fn sdf(&self, p: [f64; 3]) -> f64 {
        let d = [
            (p[0] - self.center[0]).abs() - self.half_extents[0],
            (p[1] - self.center[1]).abs() - self.half_extents[1],
            (p[2] - self.center[2]).abs() - self.half_extents[2],
        ];
        let exterior = length([d[0].max(0.0), d[1].max(0.0), d[2].max(0.0)]);
        let interior = d[0].max(d[1]).max(d[2]).min(0.0);
        exterior + interior - self.radius
    }
    fn gradient(&self, p: [f64; 3]) -> [f64; 3] {
        const EPS: f64 = 1e-5;
        let gx = self.sdf([p[0] + EPS, p[1], p[2]]) - self.sdf([p[0] - EPS, p[1], p[2]]);
        let gy = self.sdf([p[0], p[1] + EPS, p[2]]) - self.sdf([p[0], p[1] - EPS, p[2]]);
        let gz = self.sdf([p[0], p[1], p[2] + EPS]) - self.sdf([p[0], p[1], p[2] - EPS]);
        normalize([gx, gy, gz])
    }
}
