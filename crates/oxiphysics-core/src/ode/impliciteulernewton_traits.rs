//! # ImplicitEulerNewton - Trait Implementations
//!
//! This module contains trait implementations for `ImplicitEulerNewton`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::ImplicitEulerNewton;

impl Default for ImplicitEulerNewton {
    fn default() -> Self {
        Self::new(50, 1e-10)
    }
}
