//! # RbcMembraneParams - Trait Implementations
//!
//! This module contains trait implementations for `RbcMembraneParams`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::RbcMembraneParams;

impl Default for RbcMembraneParams {
    fn default() -> Self {
        Self {
            shear_modulus: 6e-6,
            bending_modulus: 2e-19,
            area_modulus: 1e-4,
            rest_length: 1.0,
            ib_stiffness: 10.0,
        }
    }
}
