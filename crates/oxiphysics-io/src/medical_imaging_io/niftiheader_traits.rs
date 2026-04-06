//! # NiftiHeader - Trait Implementations
//!
//! This module contains trait implementations for `NiftiHeader`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//! - `Display`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::NiftiHeader;
use std::fmt;

impl Default for NiftiHeader {
    fn default() -> Self {
        Self::new()
    }
}

impl fmt::Display for NiftiHeader {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let d = self.volume_dims();
        write!(
            f,
            "NiftiHeader({}x{}x{}, pixdim=[{:.2},{:.2},{:.2}], magic={})",
            d[0], d[1], d[2], self.pixdim[1], self.pixdim[2], self.pixdim[3], self.magic,
        )
    }
}
