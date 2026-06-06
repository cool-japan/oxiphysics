//! # VmdScriptWriter - Trait Implementations
//!
//! This module contains trait implementations for `VmdScriptWriter`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{ColorScheme, VmdRepresentation, VmdScriptWriter};

impl Default for VmdScriptWriter {
    fn default() -> Self {
        Self {
            color_scheme: ColorScheme::ByElement,
            representation: VmdRepresentation::Licorice,
            show_backbone: false,
            material: "Opaque".to_string(),
            add_render: false,
            render_file: "render.tga".to_string(),
        }
    }
}
