//! # EnhSphere - Trait Implementations
//!
//! This module contains trait implementations for `EnhSphere`.
//!
//! ## Implemented Traits
//!
//! - `ConvexShape`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::{ConvexShape, vadd, vnorm, vscale};
use super::types::EnhSphere;

impl ConvexShape for EnhSphere {
    fn support(&self, dir: [f64; 3]) -> [f64; 3] {
        let d = vnorm(dir);
        vadd(self.centre, vscale(d, self.radius))
    }
    fn centre(&self) -> [f64; 3] {
        self.centre
    }
}
