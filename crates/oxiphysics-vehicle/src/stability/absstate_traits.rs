//! # AbsState - Trait Implementations
//!
//! This module contains trait implementations for `AbsState`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::AbsState;

impl Default for AbsState {
    fn default() -> Self {
        Self {
            is_active: false,
            pressure: vec![1.0, 1.0, 1.0, 1.0],
            cycle_phase: String::from("increase"),
        }
    }
}
