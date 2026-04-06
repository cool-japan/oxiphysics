//! # PauliParams - Trait Implementations
//!
//! This module contains trait implementations for `PauliParams`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

use super::functions::HBAR_C;
#[allow(unused_imports)]
use super::functions::*;
use super::types::PauliParams;

impl Default for PauliParams {
    fn default() -> Self {
        Self {
            phase_space_cell: (2.0 * PI * HBAR_C).powi(3),
            strength_mev: 500.0,
            q0_fm_inv: 1.0 / 2.16,
            p0_mev_c: HBAR_C / 2.16,
        }
    }
}
