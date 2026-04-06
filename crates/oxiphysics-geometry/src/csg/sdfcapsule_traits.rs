//! # SdfCapsule - Trait Implementations
//!
//! This module contains trait implementations for `SdfCapsule`.
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
use super::types::SdfCapsule;

impl ImplicitSurface for SdfCapsule {
    fn sdf(&self, p: [f64; 3]) -> f64 {
        let dx = p[0] - self.center[0];
        let dy = p[1] - self.center[1];
        let dz = p[2] - self.center[2];
        let clamped_y = dy.clamp(-self.half_height, self.half_height);
        let dist = (dx * dx + (dy - clamped_y) * (dy - clamped_y) + dz * dz).sqrt();
        dist - self.radius
    }
    fn gradient(&self, p: [f64; 3]) -> [f64; 3] {
        let dx = p[0] - self.center[0];
        let dy = p[1] - self.center[1];
        let dz = p[2] - self.center[2];
        let clamped_y = dy.clamp(-self.half_height, self.half_height);
        let v = [dx, dy - clamped_y, dz];
        normalize(v)
    }
}
