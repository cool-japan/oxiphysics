//! # ExportOptions - Trait Implementations
//!
//! This module contains trait implementations for `ExportOptions`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::{ExportFormat, ExportOptions};

impl Default for ExportOptions {
    fn default() -> Self {
        Self {
            format: ExportFormat::Svg,
            dpi: 300,
            embed_fonts: true,
            transparent_background: false,
            tight_layout: true,
        }
    }
}
