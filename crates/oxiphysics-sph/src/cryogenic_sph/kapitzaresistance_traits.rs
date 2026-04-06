//! # KapitzaResistance - Trait Implementations
//!
//! This module contains trait implementations for `KapitzaResistance`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::KAPITZA_R0;
#[allow(unused_imports)]
use super::functions::*;
use super::types::KapitzaResistance;

impl Default for KapitzaResistance {
    fn default() -> Self {
        Self {
            r0: KAPITZA_R0,
            exponent: -3.0,
        }
    }
}
