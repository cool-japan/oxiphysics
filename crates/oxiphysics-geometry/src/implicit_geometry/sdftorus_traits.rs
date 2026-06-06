//! # SdfTorus - Trait Implementations
//!
//! This module contains trait implementations for `SdfTorus`.
//!
//! ## Implemented Traits
//!
//! - `Sdf`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
use super::types::SdfTorus;

impl Sdf for SdfTorus {
    fn dist(&self, p: [f64; 3]) -> f64 {
        let q = [len2([p[0], p[2]]) - self.major, p[1]];
        len2(q) - self.minor
    }
}
