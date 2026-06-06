//! # NamedObject - Trait Implementations
//!
//! This module contains trait implementations for `NamedObject`.
//!
//! ## Implemented Traits
//!
//! - `CatObject`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::CatObject;
use super::types::NamedObject;

impl CatObject for NamedObject {
    fn label(&self) -> &str {
        &self.label
    }
}
