//! # NoThermostat - Trait Implementations
//!
//! This module contains trait implementations for `NoThermostat`.
//!
//! ## Implemented Traits
//!
//! - `Thermostat`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::atom::AtomSet;

use super::functions::Thermostat;
#[allow(unused_imports)]
use super::functions::*;
use super::types::NoThermostat;

impl Thermostat for NoThermostat {
    fn apply(&mut self, _atoms: &mut AtomSet, _target_temp: f64, _dt: f64, _boltzmann_k: f64) {}
}
