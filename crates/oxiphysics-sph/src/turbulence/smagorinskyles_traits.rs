//! # SmagorinskyLes - Trait Implementations
//!
//! This module contains trait implementations for `SmagorinskyLes`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::SmagorinskyLes;

impl Default for SmagorinskyLes {
    fn default() -> Self {
        Self {
            cs: 0.17,
            delta: 0.1,
        }
    }
}
