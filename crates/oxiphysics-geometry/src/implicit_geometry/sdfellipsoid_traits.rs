//! # SdfEllipsoid - Trait Implementations
//!
//! This module contains trait implementations for `SdfEllipsoid`.
//!
//! ## Implemented Traits
//!
//! - `Sdf`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
#[allow(unused_imports)]
use super::functions::*;
use super::types::SdfEllipsoid;

impl Sdf for SdfEllipsoid {
    fn dist(&self, p: [f64; 3]) -> f64 {
        let [rx, ry, rz] = self.radii;
        let k0 = len([p[0] / rx, p[1] / ry, p[2] / rz]);
        let k1 = len([p[0] / (rx * rx), p[1] / (ry * ry), p[2] / (rz * rz)]);
        if k0 < 1e-14 {
            return -(rx.min(ry).min(rz));
        }
        k0 * (k0 - 1.0) / k1
    }
}
