//! # LineStyle2D - Trait Implementations
//!
//! This module contains trait implementations for `LineStyle2D`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::{LineStyle, LineStyle2D, MarkerStyle, PlotColor};

impl Default for LineStyle2D {
    fn default() -> Self {
        Self {
            color: PlotColor::blue(),
            style: LineStyle::Solid,
            width: 1.5,
            marker: MarkerStyle::None,
            marker_size: 5.0,
            label: String::new(),
        }
    }
}
