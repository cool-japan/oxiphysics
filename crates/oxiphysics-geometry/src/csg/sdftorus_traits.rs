//! # SdfTorus - Trait Implementations
//!
//! This module contains trait implementations for `SdfTorus`.
//!
//! ## Implemented Traits
//!
//! - `ImplicitSurface`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::ImplicitSurface;
use super::functions::normalize;
use super::types::SdfTorus;

impl ImplicitSurface for SdfTorus {
    fn sdf(&self, p: [f64; 3]) -> f64 {
        let q0 = (p[0] * p[0] + p[2] * p[2]).sqrt() - self.major_radius;
        (q0 * q0 + p[1] * p[1]).sqrt() - self.minor_radius
    }
    fn gradient(&self, p: [f64; 3]) -> [f64; 3] {
        const EPS: f64 = 1e-5;
        let dx = self.sdf([p[0] + EPS, p[1], p[2]]) - self.sdf([p[0] - EPS, p[1], p[2]]);
        let dy = self.sdf([p[0], p[1] + EPS, p[2]]) - self.sdf([p[0], p[1] - EPS, p[2]]);
        let dz = self.sdf([p[0], p[1], p[2] + EPS]) - self.sdf([p[0], p[1], p[2] - EPS]);
        normalize([dx, dy, dz])
    }
}
