//! # GltfMaterial - Trait Implementations
//!
//! This module contains trait implementations for `GltfMaterial`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::GltfMaterial;

impl Default for GltfMaterial {
    fn default() -> Self {
        GltfMaterial {
            name: String::new(),
            base_color_factor: [1.0, 1.0, 1.0, 1.0],
            metallic_factor: 1.0,
            roughness_factor: 1.0,
            emissive_factor: [0.0, 0.0, 0.0],
            alpha_mode: "OPAQUE".to_string(),
            alpha_cutoff: 0.5,
            double_sided: false,
        }
    }
}
