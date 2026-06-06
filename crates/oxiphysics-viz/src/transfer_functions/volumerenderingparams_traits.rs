//! # VolumeRenderingParams - Trait Implementations
//!
//! This module contains trait implementations for `VolumeRenderingParams`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{RenderQuality, TransferFunction, VolumeRenderingParams};

impl Default for VolumeRenderingParams {
    fn default() -> Self {
        Self {
            transfer_function: TransferFunction::default(),
            step_size: RenderQuality::Medium.step_size(),
            early_termination_threshold: 0.99,
            max_samples: 512,
            gradient_shading: true,
            ambient: 0.2,
            diffuse: 0.8,
            empty_space_threshold: 0.01,
        }
    }
}
