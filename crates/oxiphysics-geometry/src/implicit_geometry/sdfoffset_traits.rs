//! # SdfOffset - Trait Implementations
//!
//! This module contains trait implementations for `SdfOffset`.
//!
//! ## Implemented Traits
//!
//! - `Sdf`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
use super::types::SdfOffset;

impl<S: Sdf> Sdf for SdfOffset<S> {
    fn dist(&self, p: [f64; 3]) -> f64 {
        self.inner.dist(p) - self.offset
    }
}
