//! # JsQuat - Trait Implementations
//!
//! This module contains trait implementations for `JsQuat`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//! - `From`
//! - `From`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use crate::types::QuatWasm;

use super::types::JsQuat;

impl Default for JsQuat {
    fn default() -> Self {
        Self::identity()
    }
}

impl From<QuatWasm> for JsQuat {
    fn from(q: QuatWasm) -> Self {
        Self::from_quat(&q)
    }
}

impl From<[f64; 4]> for JsQuat {
    fn from(arr: [f64; 4]) -> Self {
        Self::from_array(arr)
    }
}
