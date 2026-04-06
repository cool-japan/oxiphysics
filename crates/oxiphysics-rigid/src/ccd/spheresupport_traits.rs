//! # SphereSupport - Trait Implementations
//!
//! This module contains trait implementations for `SphereSupport`.
//!
//! ## Implemented Traits
//!
//! - `ConvexSupport`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::functions::{ConvexSupport, add, normalize, scale};
use super::types::SphereSupport;

impl ConvexSupport for SphereSupport {
    fn support(&self, d: [f64; 3]) -> [f64; 3] {
        let nd = normalize(d);
        add(self.center, scale(nd, self.radius))
    }
}
