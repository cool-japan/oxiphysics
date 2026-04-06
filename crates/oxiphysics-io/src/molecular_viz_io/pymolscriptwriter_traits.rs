//! # PymolScriptWriter - Trait Implementations
//!
//! This module contains trait implementations for `PymolScriptWriter`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::{ColorScheme, PymolScriptWriter};

impl Default for PymolScriptWriter {
    fn default() -> Self {
        Self {
            color_scheme: ColorScheme::ByElement,
            show_spheres: false,
            add_ray: false,
            background: "white".to_string(),
            stick_radius: 0.1,
            sphere_scale: 0.3,
            orthoscopic: false,
        }
    }
}
