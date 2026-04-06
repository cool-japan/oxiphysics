//! # EnergyMonitorHook - Trait Implementations
//!
//! This module contains trait implementations for `EnergyMonitorHook`.
//!
//! ## Implemented Traits
//!
//! - `PostStepHook`
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::PostStepHook;
#[allow(unused_imports)]
use super::functions::*;
use super::types::{BodySnapshot, EnergyMonitorHook};

impl PostStepHook for EnergyMonitorHook {
    fn post_step(&self, bodies: &[BodySnapshot], _sim_time: f64, _dt: f64) {
        let ke: f64 = bodies
            .iter()
            .filter(|b| b.inv_mass > 0.0)
            .map(|b| {
                let v = b.velocity;
                let v_sq = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
                0.5 * (1.0 / b.inv_mass) * v_sq
            })
            .sum();
        self.samples
            .lock()
            .unwrap_or_else(|e| e.into_inner())
            .push(ke);
    }
}

impl Default for EnergyMonitorHook {
    fn default() -> Self {
        Self::new()
    }
}
