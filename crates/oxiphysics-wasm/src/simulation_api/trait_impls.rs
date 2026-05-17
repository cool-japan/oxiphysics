//! # WasmSimulationConfig - Trait Implementations
//!
//! This module contains trait implementations for `WasmSimulationConfig`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{BroadPhaseAlgorithm, IntegrationMethod};
use super::types_4::WasmSimulationConfig;

impl Default for WasmSimulationConfig {
    fn default() -> Self {
        Self {
            gravity: [0.0, -9.81, 0.0],
            dt: 1.0 / 60.0,
            velocity_iters: 8,
            position_iters: 4,
            broad_phase: BroadPhaseAlgorithm::Sap,
            integration: IntegrationMethod::SemiImplicitEuler,
            ccd_enabled: false,
            max_linear_vel: 200.0,
            max_angular_vel: 200.0,
            linear_damping: 0.0,
            angular_damping: 0.0,
            allow_sleeping: true,
            sleep_linear_threshold: 0.002,
            sleep_angular_threshold: 0.005,
            sleep_time_threshold: 1.0,
        }
    }
}
