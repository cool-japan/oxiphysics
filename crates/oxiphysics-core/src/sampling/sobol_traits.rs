//! # Sobol - Trait Implementations
//!
//! This module contains trait implementations for `Sobol`.
//!
//! ## Implemented Traits
//!
//! - `QuasiRandomSequence`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::QuasiRandomSequence;
#[allow(unused_imports)]
use super::functions::*;
use super::types::Sobol;

#[allow(dead_code)]
impl QuasiRandomSequence for Sobol {
    fn next(&mut self) -> Vec<f64> {
        vec![self.next_1d()]
    }
}
