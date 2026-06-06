//! # PbrMaterialBuilder - Trait Implementations
//!
//! This module contains trait implementations for `PbrMaterialBuilder`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::PbrMaterialBuilder;

impl Default for PbrMaterialBuilder {
    fn default() -> Self {
        Self {
            name: "material".to_string(),
            base_color_factor: [1.0, 1.0, 1.0, 1.0],
            metallic_factor: 0.0,
            roughness_factor: 0.5,
            emissive_factor: [0.0, 0.0, 0.0],
            double_sided: false,
            alpha_mode: "OPAQUE".to_string(),
            alpha_cutoff: 0.5,
        }
    }
}
