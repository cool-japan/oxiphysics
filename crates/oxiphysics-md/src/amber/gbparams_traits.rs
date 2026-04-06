//! # GbParams - Trait Implementations
//!
//! This module contains trait implementations for `GbParams`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::GbParams;

impl Default for GbParams {
    fn default() -> Self {
        Self {
            eps_solvent: 78.5,
            eps_solute: 1.0,
            coulomb_k: 138.935_485,
            surface_tension: 2.2677,
            sa_offset: 0.0,
        }
    }
}
