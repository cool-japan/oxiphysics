//! # Complex - Trait Implementations
//!
//! This module contains trait implementations for `Complex`.
//!
//! ## Implemented Traits
//!
//! - `Display`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::Complex;

impl std::fmt::Display for Complex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.im >= 0.0 {
            write!(f, "{:.4} + {:.4}j", self.re, self.im)
        } else {
            write!(f, "{:.4} - {:.4}j", self.re, -self.im)
        }
    }
}
