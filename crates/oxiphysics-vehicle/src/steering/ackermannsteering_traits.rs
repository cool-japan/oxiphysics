//! # AckermannSteering - Trait Implementations
//!
//! This module contains trait implementations for `AckermannSteering`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)
use super::types::AckermannSteering;

impl Default for AckermannSteering {
    fn default() -> Self {
        Self {
            wheelbase: 2.5,
            track_width: 1.5,
            max_steer_angle: 0.6,
            ackermann_factor: 1.0,
        }
    }
}
