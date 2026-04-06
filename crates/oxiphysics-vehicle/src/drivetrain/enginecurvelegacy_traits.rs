//! # EngineCurveLegacy - Trait Implementations
//!
//! This module contains trait implementations for `EngineCurveLegacy`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::EngineCurveLegacy;

impl Default for EngineCurveLegacy {
    fn default() -> Self {
        Self::new(vec![
            (800.0, 100.0),
            (2000.0, 250.0),
            (3500.0, 350.0),
            (5000.0, 380.0),
            (6000.0, 340.0),
            (7000.0, 280.0),
        ])
    }
}
