//! # CoupledArrowConfig - Trait Implementations
//!
//! This module contains trait implementations for `CoupledArrowConfig`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//! - `Default`
//! - `Default`
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{
    CoupledArrowConfig, EmFieldLineConfig, FieldLineMethod, Rgba, ThermalMechanicalOverlayConfig,
};

impl Default for CoupledArrowConfig {
    fn default() -> Self {
        Self {
            heat_flux_scale: 1.0,
            displacement_scale: 1.0,
            heat_flux_color: Rgba::new(1.0, 0.4, 0.0, 1.0),
            displacement_color: Rgba::new(0.2, 0.6, 1.0, 1.0),
            min_magnitude: 1e-12,
        }
    }
}

impl Default for EmFieldLineConfig {
    fn default() -> Self {
        Self {
            step_size: 0.05,
            max_steps: 500,
            min_magnitude: 1e-10,
            method: FieldLineMethod::Rk4,
        }
    }
}

impl Default for Rgba {
    fn default() -> Self {
        Self::white()
    }
}

impl Default for ThermalMechanicalOverlayConfig {
    fn default() -> Self {
        Self {
            temp_min: 0.0,
            temp_max: 1000.0,
            disp_min: 0.0,
            disp_max: 0.01,
            thermal_weight: 0.5,
        }
    }
}
