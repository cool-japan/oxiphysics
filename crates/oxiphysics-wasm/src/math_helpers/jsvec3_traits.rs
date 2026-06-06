//! # JsVec3 - Trait Implementations
//!
//! This module contains trait implementations for `JsVec3`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//! - `From`
//! - `From`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)
use crate::types::Vec3Wasm;

use super::types::JsVec3;

impl Default for JsVec3 {
    fn default() -> Self {
        Self::zero()
    }
}

impl From<Vec3Wasm> for JsVec3 {
    fn from(v: Vec3Wasm) -> Self {
        Self::from_vec3(&v)
    }
}

impl From<[f64; 3]> for JsVec3 {
    fn from(arr: [f64; 3]) -> Self {
        Self::from_array(arr)
    }
}
