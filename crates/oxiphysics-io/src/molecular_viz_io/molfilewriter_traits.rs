//! # MolfileWriter - Trait Implementations
//!
//! This module contains trait implementations for `MolfileWriter`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::MolfileWriter;

impl Default for MolfileWriter {
    fn default() -> Self {
        Self {
            author: "OxiPhysics".to_string(),
            sdf_mode: false,
            include_stereo: false,
        }
    }
}
