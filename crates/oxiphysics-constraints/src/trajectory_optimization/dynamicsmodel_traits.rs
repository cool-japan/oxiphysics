//! # DynamicsModel - Trait Implementations
//!
//! This module contains trait implementations for `DynamicsModel`.
//!
//! ## Implemented Traits
//!
//! - `Debug`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::DynamicsModel;

impl std::fmt::Debug for DynamicsModel {
    fn fmt(&self, fmt: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        fmt.debug_struct("DynamicsModel")
            .field("n_state", &self.n_state)
            .field("n_control", &self.n_control)
            .finish()
    }
}
