//! # SdfScale - Trait Implementations
//!
//! This module contains trait implementations for `SdfScale`.
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
use super::types::SdfScale;

impl<S: Sdf> Sdf for SdfScale<S> {
    fn dist(&self, p: [f64; 3]) -> f64 {
        self.inner.dist(scale(p, 1.0 / self.factor)) * self.factor
    }
}
