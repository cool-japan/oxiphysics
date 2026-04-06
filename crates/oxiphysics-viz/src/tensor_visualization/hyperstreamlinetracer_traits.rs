//! # HyperstreamlineTracer - Trait Implementations
//!
//! This module contains trait implementations for `HyperstreamlineTracer`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use std::f64::consts::PI;

use super::types::HyperstreamlineTracer;

impl Default for HyperstreamlineTracer {
    fn default() -> Self {
        Self {
            step_size: 0.5,
            max_steps: 500,
            min_fa: 0.1,
            max_curvature_angle: PI / 3.0,
            bidirectional: true,
        }
    }
}
