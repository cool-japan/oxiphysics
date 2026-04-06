//! # JsTransform - Trait Implementations
//!
//! This module contains trait implementations for `JsTransform`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//! - `From`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use crate::types::TransformWasm;

use super::types::JsTransform;

impl Default for JsTransform {
    fn default() -> Self {
        Self::identity()
    }
}

impl From<TransformWasm> for JsTransform {
    fn from(t: TransformWasm) -> Self {
        Self::from_transform(&t)
    }
}
