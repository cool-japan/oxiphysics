//! # StressTensorField - Trait Implementations
//!
//! This module contains trait implementations for `StressTensorField`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::StressTensorField;

impl Default for StressTensorField {
    fn default() -> Self {
        Self {
            glyph_scale: 1.0,
            tension_color: [0.9, 0.1, 0.1, 1.0],
            compression_color: [0.1, 0.1, 0.9, 1.0],
            normalize_glyphs: false,
        }
    }
}
