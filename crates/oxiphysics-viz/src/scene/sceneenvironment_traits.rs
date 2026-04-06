//! # SceneEnvironment - Trait Implementations
//!
//! This module contains trait implementations for `SceneEnvironment`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::SceneEnvironment;

impl Default for SceneEnvironment {
    fn default() -> Self {
        Self::clear_sky()
    }
}
