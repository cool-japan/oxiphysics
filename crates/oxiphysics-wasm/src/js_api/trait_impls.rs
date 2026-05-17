//! # JsBodyDesc - Trait Implementations
//!
//! This module contains trait implementations for `JsBodyDesc`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::types::SimulationConfig;

use super::types::{JsBodyDesc, JsPhysicsConfig};

impl Default for JsBodyDesc {
    fn default() -> Self {
        Self {
            mass: 1.0,
            position: [0.0; 3],
            velocity: [0.0; 3],
            shape: "sphere".to_string(),
            radius: 0.5,
            half_extents: [0.5; 3],
            height: 1.0,
            restitution: 0.3,
            friction: 0.5,
            tag: String::new(),
        }
    }
}

impl Default for JsPhysicsConfig {
    fn default() -> Self {
        let sim = SimulationConfig::default();
        Self {
            gravity_x: sim.gravity[0],
            gravity_y: sim.gravity[1],
            gravity_z: sim.gravity[2],
            fixed_dt: sim.fixed_dt,
            max_substeps: sim.max_substeps,
            solver_iterations: sim.solver_iterations,
            ccd_enabled: sim.ccd_enabled,
            sleeping_enabled: sim.sleeping_enabled,
            linear_sleep_threshold: sim.linear_sleep_threshold,
            angular_sleep_threshold: sim.angular_sleep_threshold,
        }
    }
}
