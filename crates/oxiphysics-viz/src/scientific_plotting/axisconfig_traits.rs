//! # AxisConfig - Trait Implementations
//!
//! This module contains trait implementations for `AxisConfig`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{AxisConfig, AxisScale};

impl Default for AxisConfig {
    fn default() -> Self {
        Self {
            label: String::new(),
            min: None,
            max: None,
            scale: AxisScale::Linear,
            num_ticks: 0,
            show_grid: true,
        }
    }
}
