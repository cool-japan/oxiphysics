#[allow(unused_imports)]
use super::functions::*;
// # HalfEdgeMesh - Trait Implementations
//
// This module contains trait implementations for `HalfEdgeMesh`.
//
// ## Implemented Traits
//
// - `Default`
//
// 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::HalfEdgeMesh;

impl Default for HalfEdgeMesh {
    fn default() -> Self {
        Self::new()
    }
}
