//! # DicomReader - Trait Implementations
//!
//! This module contains trait implementations for `DicomReader`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::DicomReader;

impl Default for DicomReader {
    fn default() -> Self {
        Self::new()
    }
}
