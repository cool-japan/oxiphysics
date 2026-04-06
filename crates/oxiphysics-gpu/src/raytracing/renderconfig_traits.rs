//! # RenderConfig - Trait Implementations
//!
//! This module contains trait implementations for `RenderConfig`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::RenderConfig;

impl Default for RenderConfig {
    fn default() -> Self {
        Self {
            width: 800,
            height: 600,
            spp: 16,
            max_depth: 8,
            soft_shadows: true,
            ambient_occlusion: true,
            depth_of_field: false,
            ao_samples: 16,
            shadow_samples: 8,
            background: [0.1, 0.15, 0.25],
            ambient: [0.05, 0.05, 0.05],
            tonemap: 3,
        }
    }
}
