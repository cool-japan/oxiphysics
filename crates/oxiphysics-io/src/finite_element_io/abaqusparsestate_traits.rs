//! # AbaqusParseState - Trait Implementations
//!
//! This module contains trait implementations for `AbaqusParseState`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::{AbaqusParseState, FeElementType};

impl Default for AbaqusParseState {
    fn default() -> Self {
        Self {
            node_section: false,
            element_section: false,
            current_elem_type: FeElementType::Unknown(String::new()),
            node_set_section: None,
            elem_set_section: None,
            material_name: None,
            in_elastic: false,
        }
    }
}
