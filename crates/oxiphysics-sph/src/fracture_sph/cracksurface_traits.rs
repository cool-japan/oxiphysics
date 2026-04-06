//! # CrackSurface - Trait Implementations
//!
//! This module contains trait implementations for `CrackSurface`.
//!
//! ## Implemented Traits
//!
//! - `CrackSurfaceHelper`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::CrackSurfaceHelper;
#[allow(unused_imports)]
use super::functions::*;
use super::types::CrackSurface;

impl CrackSurfaceHelper for CrackSurface {
    fn smooth_h_or_default(&self, default: f64) -> f64 {
        default
    }
}
