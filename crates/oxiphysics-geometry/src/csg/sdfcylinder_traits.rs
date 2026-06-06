//! # SdfCylinder - Trait Implementations
//!
//! This module contains trait implementations for `SdfCylinder`.
//!
//! ## Implemented Traits
//!
//! - `ImplicitSurface`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::ImplicitSurface;
use super::functions::normalize;
use super::types::SdfCylinder;

impl ImplicitSurface for SdfCylinder {
    fn sdf(&self, p: [f64; 3]) -> f64 {
        let dx = p[0] - self.center[0];
        let dz = p[2] - self.center[2];
        (dx * dx + dz * dz).sqrt() - self.radius
    }
    fn gradient(&self, p: [f64; 3]) -> [f64; 3] {
        let dx = p[0] - self.center[0];
        let dz = p[2] - self.center[2];
        normalize([dx, 0.0, dz])
    }
}
