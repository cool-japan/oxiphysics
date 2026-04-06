//! # SdfBoundedProxy - Trait Implementations
//!
//! This module contains trait implementations for `SdfBoundedProxy`.
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
use super::types::SdfBoundedProxy;

impl<S: Sdf> Sdf for SdfBoundedProxy<S> {
    fn dist(&self, p: [f64; 3]) -> f64 {
        let bs_dist = len(sub(p, self.bsphere_center)) - self.bsphere_radius;
        if bs_dist > 0.0 {
            return bs_dist;
        }
        self.inner.dist(p)
    }
}
