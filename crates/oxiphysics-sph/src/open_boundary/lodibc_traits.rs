//! # LodiBC - Trait Implementations
//!
//! This module contains trait implementations for `LodiBC`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::LodiBC;

impl Default for LodiBC {
    fn default() -> Self {
        Self {
            c_sound: 340.0,
            density: 1.225,
            sigma: 0.25,
            p_ref: 101325.0,
            length: 1.0,
        }
    }
}
