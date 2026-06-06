//! # BloomParams - Trait Implementations
//!
//! This module contains trait implementations for `BloomParams`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//! - `Default`
//! - `Default`
//! - `Default`
//! - `Default`
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{BloomParams, PyDebugOverlay, PyPostProcessor, RayMarchSettings, SsaoParams};
use super::types_3::PySceneGraph;

impl Default for BloomParams {
    fn default() -> Self {
        Self {
            threshold: 1.0,
            radius: 8.0,
            intensity: 0.5,
        }
    }
}

impl Default for PyDebugOverlay {
    fn default() -> Self {
        Self::new()
    }
}

impl Default for PyPostProcessor {
    fn default() -> Self {
        Self::new()
    }
}

impl Default for PySceneGraph {
    fn default() -> Self {
        Self::new()
    }
}

impl Default for RayMarchSettings {
    fn default() -> Self {
        Self {
            step_size: 0.01,
            max_steps: 512,
            opacity_threshold: 0.99,
            jitter: true,
        }
    }
}

impl Default for SsaoParams {
    fn default() -> Self {
        Self {
            radius: 0.5,
            bias: 0.025,
            num_samples: 16,
            power: 2.0,
        }
    }
}
