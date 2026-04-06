//! # SurfaceLighting - Trait Implementations
//!
//! This module contains trait implementations for `SurfaceLighting`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::SurfaceLighting;

impl Default for SurfaceLighting {
    fn default() -> Self {
        Self {
            ambient: 0.15,
            diffuse: 0.8,
            specular: 0.3,
            shininess: 32.0,
            light_dir: [0.577, 0.577, 0.577],
        }
    }
}
