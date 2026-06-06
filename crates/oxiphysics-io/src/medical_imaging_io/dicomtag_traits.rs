//! # DicomTag - Trait Implementations
//!
//! This module contains trait implementations for `DicomTag`.
//!
//! ## Implemented Traits
//!
//! - `Display`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::DicomTag;
use std::fmt;

impl fmt::Display for DicomTag {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({:04X},{:04X}) {}", self.group, self.element, self.vr)
    }
}
