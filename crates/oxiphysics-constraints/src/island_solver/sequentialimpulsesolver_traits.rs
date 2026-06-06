//! # SequentialImpulseSolver - Trait Implementations
//!
//! This module contains trait implementations for `SequentialImpulseSolver`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::SequentialImpulseSolver;

impl Default for SequentialImpulseSolver {
    fn default() -> Self {
        Self { iterations: 10 }
    }
}
