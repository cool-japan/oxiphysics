//! # IntegrateVelocityKernel - Trait Implementations
//!
//! This module contains trait implementations for `IntegrateVelocityKernel`.
//!
//! ## Implemented Traits
//!
//! - `ComputeKernel`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use crate::compute::ComputeKernel;

use super::types::IntegrateVelocityKernel;

#[allow(clippy::needless_range_loop)]
impl ComputeKernel for IntegrateVelocityKernel {
    fn name(&self) -> &str {
        "IntegrateVelocityKernel"
    }
    fn execute(&self, inputs: &[&[f64]], outputs: &mut [Vec<f64>], work_size: usize) {
        if inputs.len() < 4 || outputs.is_empty() {
            return;
        }
        let vel = inputs[0];
        let force = inputs[1];
        let inv_mass = inputs[2];
        let dt = inputs[3][0];
        let n = work_size;
        let mut new_vel = vec![0.0; n * 3];
        for i in 0..n {
            let im = inv_mass[i];
            new_vel[i * 3] = vel[i * 3] + force[i * 3] * im * dt;
            new_vel[i * 3 + 1] = vel[i * 3 + 1] + force[i * 3 + 1] * im * dt;
            new_vel[i * 3 + 2] = vel[i * 3 + 2] + force[i * 3 + 2] * im * dt;
        }
        outputs[0] = new_vel;
    }
}
