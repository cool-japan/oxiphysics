//! # QmdParams - Trait Implementations
//!
//! This module contains trait implementations for `QmdParams`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::QmdParams;

impl Default for QmdParams {
    fn default() -> Self {
        Self {
            wavepacket_width_fm2: 2.16 * 2.16,
            t1_mev_fm3: -356.0,
            t3_mev: 303.0,
            gamma: 7.0 / 6.0,
            cs_mev: 25.0,
            coulomb_coeff: 1.44,
            momentum_dep_coeff: 1.57,
        }
    }
}
