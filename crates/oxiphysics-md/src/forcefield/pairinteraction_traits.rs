//! # PairInteraction - Trait Implementations
//!
//! This module contains trait implementations for `PairInteraction`.
//!
//! ## Implemented Traits
//!
//! - `Debug`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::PairInteraction;

impl std::fmt::Debug for PairInteraction {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("PairInteraction")
            .field("type_i", &self.type_i)
            .field("type_j", &self.type_j)
            .finish()
    }
}
