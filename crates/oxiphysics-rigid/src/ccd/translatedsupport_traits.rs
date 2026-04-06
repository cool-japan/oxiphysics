//! # TranslatedSupport - Trait Implementations
//!
//! This module contains trait implementations for `TranslatedSupport`.
//!
//! ## Implemented Traits
//!
//! - `ConvexSupport`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::functions::{ConvexSupport, add};
use super::types::TranslatedSupport;

impl<'a, S: ConvexSupport> ConvexSupport for TranslatedSupport<'a, S> {
    fn support(&self, d: [f64; 3]) -> [f64; 3] {
        add(self.shape.support(d), self.offset)
    }
}
