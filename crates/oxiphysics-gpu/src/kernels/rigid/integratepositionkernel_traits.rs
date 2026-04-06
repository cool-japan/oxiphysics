//! # IntegratePositionKernel - Trait Implementations
//!
//! This module contains trait implementations for `IntegratePositionKernel`.
//!
//! ## Implemented Traits
//!
//! - `ComputeKernel`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use crate::compute::ComputeKernel;

use super::types::IntegratePositionKernel;

#[allow(clippy::needless_range_loop)]
impl ComputeKernel for IntegratePositionKernel {
    fn name(&self) -> &str {
        "IntegratePositionKernel"
    }
    fn execute(&self, inputs: &[&[f64]], outputs: &mut [Vec<f64>], work_size: usize) {
        if inputs.len() < 3 || outputs.is_empty() {
            return;
        }
        let pos = inputs[0];
        let vel = inputs[1];
        let dt = inputs[2][0];
        let n = work_size;
        let mut new_pos = vec![0.0; n * 3];
        for i in 0..n {
            new_pos[i * 3] = pos[i * 3] + vel[i * 3] * dt;
            new_pos[i * 3 + 1] = pos[i * 3 + 1] + vel[i * 3 + 1] * dt;
            new_pos[i * 3 + 2] = pos[i * 3 + 2] + vel[i * 3 + 2] * dt;
        }
        outputs[0] = new_pos;
    }
}
