//! # ChemicalJsonWriter - Trait Implementations
//!
//! This module contains trait implementations for `ChemicalJsonWriter`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::ChemicalJsonWriter;

impl Default for ChemicalJsonWriter {
    fn default() -> Self {
        Self {
            include_charges: true,
            include_unit_cell: true,
            pretty_print: true,
            indent: "  ".to_string(),
        }
    }
}
