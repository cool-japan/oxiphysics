//! # ConstraintState - Trait Implementations
//!
//! This module contains trait implementations for `ConstraintState`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//! - `Display`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::ConstraintState;

impl std::fmt::Display for ConstraintState {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ConstraintState::Uninitialised => write!(f, "Uninitialised"),
            ConstraintState::Active => write!(f, "Active"),
            ConstraintState::Disabled => write!(f, "Disabled"),
            ConstraintState::Removed => write!(f, "Removed"),
            ConstraintState::Sleeping => write!(f, "Sleeping"),
        }
    }
}
