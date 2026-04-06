//! # NiftiReader - Trait Implementations
//!
//! This module contains trait implementations for `NiftiReader`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::NiftiReader;

impl Default for NiftiReader {
    fn default() -> Self {
        Self::new()
    }
}
