//! # SceneCameraProjection - Trait Implementations
//!
//! This module contains trait implementations for `SceneCameraProjection`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::SceneCameraProjection;

impl Default for SceneCameraProjection {
    fn default() -> Self {
        SceneCameraProjection::Perspective {
            fov_y: std::f64::consts::FRAC_PI_4,
            aspect: 16.0 / 9.0,
            near: 0.01,
            far: 1000.0,
        }
    }
}
