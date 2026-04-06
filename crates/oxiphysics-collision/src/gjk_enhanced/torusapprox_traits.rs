//! # TorusApprox - Trait Implementations
//!
//! This module contains trait implementations for `TorusApprox`.
//!
//! ## Implemented Traits
//!
//! - `ConvexShape`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::ConvexShape;
#[allow(unused_imports)]
use super::functions::*;
use super::types::TorusApprox;

impl ConvexShape for TorusApprox {
    fn support(&self, dir: [f64; 3]) -> [f64; 3] {
        let d_xz = vnorm([dir[0], 0.0, dir[2]]);
        let ring = vscale(d_xz, self.major_r);
        let dn = vnorm(dir);
        let support = vadd(ring, vscale(dn, self.minor_r));
        vadd(self.centre, support)
    }
    fn centre(&self) -> [f64; 3] {
        self.centre
    }
}
