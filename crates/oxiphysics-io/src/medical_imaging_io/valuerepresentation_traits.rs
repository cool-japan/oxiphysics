//! # ValueRepresentation - Trait Implementations
//!
//! This module contains trait implementations for `ValueRepresentation`.
//!
//! ## Implemented Traits
//!
//! - `Display`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::ValueRepresentation;
use std::fmt;

impl fmt::Display for ValueRepresentation {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.as_str())
    }
}
