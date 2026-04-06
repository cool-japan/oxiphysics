//! # SdfRevolution - Trait Implementations
//!
//! This module contains trait implementations for `SdfRevolution`.
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
use super::types::SdfRevolution;

impl<F: Fn(f64, f64) -> f64 + Send + Sync> Sdf for SdfRevolution<F> {
    fn dist(&self, p: [f64; 3]) -> f64 {
        let r = (p[0] * p[0] + p[2] * p[2]).sqrt();
        (self.profile)(r, p[1])
    }
}
