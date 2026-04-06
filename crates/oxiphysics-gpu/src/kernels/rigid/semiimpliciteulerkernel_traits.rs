//! # SemiImplicitEulerKernel - Trait Implementations
//!
//! This module contains trait implementations for `SemiImplicitEulerKernel`.
//!
//! ## Implemented Traits
//!
//! - `ComputeKernel`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use crate::compute::ComputeKernel;

use super::types::SemiImplicitEulerKernel;

#[allow(clippy::needless_range_loop)]
impl ComputeKernel for SemiImplicitEulerKernel {
    fn name(&self) -> &str {
        "SemiImplicitEulerKernel"
    }
    fn execute(&self, inputs: &[&[f64]], outputs: &mut [Vec<f64>], work_size: usize) {
        if inputs.len() < 5 || outputs.len() < 2 {
            return;
        }
        let pos = inputs[0];
        let vel = inputs[1];
        let force = inputs[2];
        let inv_mass = inputs[3];
        let dt = inputs[4][0];
        let n = work_size;
        let mut new_vel = vec![0.0; n * 3];
        let mut new_pos = vec![0.0; n * 3];
        for i in 0..n {
            let im = inv_mass[i];
            new_vel[i * 3] = vel[i * 3] + force[i * 3] * im * dt;
            new_vel[i * 3 + 1] = vel[i * 3 + 1] + force[i * 3 + 1] * im * dt;
            new_vel[i * 3 + 2] = vel[i * 3 + 2] + force[i * 3 + 2] * im * dt;
            new_pos[i * 3] = pos[i * 3] + new_vel[i * 3] * dt;
            new_pos[i * 3 + 1] = pos[i * 3 + 1] + new_vel[i * 3 + 1] * dt;
            new_pos[i * 3 + 2] = pos[i * 3 + 2] + new_vel[i * 3 + 2] * dt;
        }
        outputs[0] = new_vel;
        outputs[1] = new_pos;
    }
}
