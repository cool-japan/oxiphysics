//! # SdfShell - Trait Implementations
//!
//! This module contains trait implementations for `SdfShell`.
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
use super::types::SdfShell;

impl<S: Sdf> Sdf for SdfShell<S> {
    fn dist(&self, p: [f64; 3]) -> f64 {
        self.inner.dist(p).abs() - self.thickness * 0.5
    }
}
