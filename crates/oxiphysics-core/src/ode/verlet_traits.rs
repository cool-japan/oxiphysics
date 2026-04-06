//! # Verlet - Trait Implementations
//!
//! This module contains trait implementations for `Verlet`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::Verlet;

impl Default for Verlet {
    fn default() -> Self {
        Self::new()
    }
}
