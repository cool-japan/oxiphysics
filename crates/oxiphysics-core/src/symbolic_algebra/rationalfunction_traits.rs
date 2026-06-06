//! # RationalFunction - Trait Implementations
//!
//! This module contains trait implementations for `RationalFunction`.
//!
//! ## Implemented Traits
//!
//! - `Display`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::RationalFunction;
use std::fmt;

impl fmt::Display for RationalFunction {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({}) / ({})", self.numerator, self.denominator)
    }
}
