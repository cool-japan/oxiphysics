//! # EnhEllipsoid - Trait Implementations
//!
//! This module contains trait implementations for `EnhEllipsoid`.
//!
//! ## Implemented Traits
//!
//! - `ConvexShape`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::{ConvexShape, vadd, vdot, vscale};
use super::types::EnhEllipsoid;

impl ConvexShape for EnhEllipsoid {
    /// Furthest point of the axis-aligned ellipsoid along `dir`.
    ///
    /// For semi-axes `a_i` and `A = diag(a_i)` the support point is
    /// `c + A² d / sqrt(dᵀ A² d)`, which lies on the surface and maximises the
    /// projection on `d`. A direction along which the ellipsoid has no extent
    /// (the zero vector, or the normal of a flat ellipsoid) is supported by the
    /// centre.
    fn support(&self, dir: [f64; 3]) -> [f64; 3] {
        let a = self.semi_axes;
        let scaled = [
            a[0] * a[0] * dir[0],
            a[1] * a[1] * dir[1],
            a[2] * a[2] * dir[2],
        ];
        let norm_sq = vdot(scaled, dir);
        if norm_sq > 0.0 && norm_sq.is_finite() {
            vadd(self.centre, vscale(scaled, 1.0 / norm_sq.sqrt()))
        } else {
            self.centre
        }
    }
    fn centre(&self) -> [f64; 3] {
        self.centre
    }
}
