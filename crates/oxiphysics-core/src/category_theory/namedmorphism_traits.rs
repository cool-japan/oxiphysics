//! # NamedMorphism - Trait Implementations
//!
//! This module contains trait implementations for `NamedMorphism`.
//!
//! ## Implemented Traits
//!
//! - `CatMorphism`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::CatMorphism;
#[allow(unused_imports)]
use super::functions::*;
use super::types::NamedMorphism;

impl CatMorphism for NamedMorphism {
    fn name(&self) -> &str {
        &self.morph_name
    }
    fn domain(&self) -> &str {
        &self.dom
    }
    fn codomain(&self) -> &str {
        &self.cod
    }
}
