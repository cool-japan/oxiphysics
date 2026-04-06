//! # FlowRateController - Trait Implementations
//!
//! This module contains trait implementations for `FlowRateController`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::FlowRateController;

impl Default for FlowRateController {
    fn default() -> Self {
        Self {
            q_target: 1.0,
            kp: 1.0,
            ki: 0.1,
            kd: 0.01,
            integral: 0.0,
            prev_error: 0.0,
            max_correction: 5.0,
        }
    }
}
