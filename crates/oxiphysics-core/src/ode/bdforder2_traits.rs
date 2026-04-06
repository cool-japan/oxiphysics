//! # BdfOrder2 - Trait Implementations
//!
//! This module contains trait implementations for `BdfOrder2`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::BdfOrder2;

impl Default for BdfOrder2 {
    fn default() -> Self {
        Self::new(50, 1e-10)
    }
}
