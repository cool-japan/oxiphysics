//! # SdfTranslate - Trait Implementations
//!
//! This module contains trait implementations for `SdfTranslate`.
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
use super::types::SdfTranslate;

impl<S: Sdf> Sdf for SdfTranslate<S> {
    fn dist(&self, p: [f64; 3]) -> f64 {
        self.inner.dist(sub(p, self.offset))
    }
}
