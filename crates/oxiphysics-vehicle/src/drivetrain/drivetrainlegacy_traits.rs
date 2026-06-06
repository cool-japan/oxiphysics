//! # DrivetrainLegacy - Trait Implementations
//!
//! This module contains trait implementations for `DrivetrainLegacy`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)
use super::types::{
    DifferentialMode, DriveLayout, DrivetrainLegacy, EngineCurveLegacy, GearboxLegacy,
};

impl Default for DrivetrainLegacy {
    fn default() -> Self {
        Self {
            engine: EngineCurveLegacy::default(),
            gearbox: GearboxLegacy::default(),
            front_diff: DifferentialMode::Open,
            rear_diff: DifferentialMode::Open,
            layout: DriveLayout::RearWheelDrive,
            awd_front_bias: 0.4,
            max_brake_torque: 3000.0,
            engine_rpm: 800.0,
        }
    }
}
