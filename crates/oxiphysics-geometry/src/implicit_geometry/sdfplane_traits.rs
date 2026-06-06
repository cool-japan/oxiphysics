//! # SdfPlane - Trait Implementations
//!
//! This module contains trait implementations for `SdfPlane`.
//!
//! ## Implemented Traits
//!
//! - `Sdf`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
use super::types::SdfPlane;

impl Sdf for SdfPlane {
    fn dist(&self, p: [f64; 3]) -> f64 {
        dot(self.normal, p) - self.offset
    }
}
