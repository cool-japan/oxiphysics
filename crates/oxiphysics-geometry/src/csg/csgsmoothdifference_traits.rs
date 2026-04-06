//! # CsgSmoothDifference - Trait Implementations
//!
//! This module contains trait implementations for `CsgSmoothDifference`.
//!
//! ## Implemented Traits
//!
//! - `ImplicitSurface`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::ImplicitSurface;
#[allow(unused_imports)]
use super::functions::*;
use super::functions::{clamp, mix, normalize};
use super::types::CsgSmoothDifference;

impl ImplicitSurface for CsgSmoothDifference {
    fn sdf(&self, p: [f64; 3]) -> f64 {
        let da = self.a.sdf(p);
        let db = -self.b.sdf(p);
        let h = clamp(0.5 - 0.5 * (db - da) / self.k, 0.0, 1.0);
        mix(db, da, h) + self.k * h * (1.0 - h)
    }
    fn gradient(&self, p: [f64; 3]) -> [f64; 3] {
        const EPS: f64 = 1e-5;
        let dx = self.sdf([p[0] + EPS, p[1], p[2]]) - self.sdf([p[0] - EPS, p[1], p[2]]);
        let dy = self.sdf([p[0], p[1] + EPS, p[2]]) - self.sdf([p[0], p[1] - EPS, p[2]]);
        let dz = self.sdf([p[0], p[1], p[2] + EPS]) - self.sdf([p[0], p[1], p[2] - EPS]);
        normalize([dx, dy, dz])
    }
}
