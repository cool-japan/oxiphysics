//! # Se3 - Trait Implementations
//!
//! This module contains trait implementations for `Se3`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::Se3;

impl Default for Se3 {
    fn default() -> Self {
        Self::identity()
    }
}
