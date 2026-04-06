//! # SdfCappedCylinder - Trait Implementations
//!
//! This module contains trait implementations for `SdfCappedCylinder`.
//!
//! ## Implemented Traits
//!
//! - `ImplicitSurface`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::ImplicitSurface;
use super::functions::normalize;
#[allow(unused_imports)]
use super::functions::*;
use super::types::SdfCappedCylinder;

impl ImplicitSurface for SdfCappedCylinder {
    fn sdf(&self, p: [f64; 3]) -> f64 {
        let dx = p[0] - self.center[0];
        let dy = p[1] - self.center[1];
        let dz = p[2] - self.center[2];
        let r_dist = (dx * dx + dz * dz).sqrt() - self.radius;
        let h_dist = dy.abs() - self.half_height;
        let ex = r_dist.max(0.0);
        let ey = h_dist.max(0.0);
        r_dist.max(h_dist).min(0.0) + (ex * ex + ey * ey).sqrt()
    }
    fn gradient(&self, p: [f64; 3]) -> [f64; 3] {
        const EPS: f64 = 1e-5;
        let dx = self.sdf([p[0] + EPS, p[1], p[2]]) - self.sdf([p[0] - EPS, p[1], p[2]]);
        let dy = self.sdf([p[0], p[1] + EPS, p[2]]) - self.sdf([p[0], p[1] - EPS, p[2]]);
        let dz = self.sdf([p[0], p[1], p[2] + EPS]) - self.sdf([p[0], p[1], p[2] - EPS]);
        normalize([dx, dy, dz])
    }
}
