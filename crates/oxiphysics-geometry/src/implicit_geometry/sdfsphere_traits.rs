//! # SdfSphere - Trait Implementations
//!
//! This module contains trait implementations for `SdfSphere`.
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
use super::types::SdfSphere;

impl Sdf for SdfSphere {
    fn dist(&self, p: [f64; 3]) -> f64 {
        len(p) - self.radius
    }
}
