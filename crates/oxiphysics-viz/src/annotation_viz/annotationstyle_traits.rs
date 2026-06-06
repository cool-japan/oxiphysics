//! # AnnotationStyle - Trait Implementations
//!
//! This module contains trait implementations for `AnnotationStyle`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::AnnotationStyle;

impl Default for AnnotationStyle {
    fn default() -> Self {
        Self {
            line_width: 1.0,
            arrow_length: 0.1,
            arrow_width: 0.05,
            font_size: 12.0,
            color: [0.0, 0.0, 0.0, 1.0],
            text_bg_color: [1.0, 1.0, 1.0, 0.8],
            show_text_bg: true,
        }
    }
}
