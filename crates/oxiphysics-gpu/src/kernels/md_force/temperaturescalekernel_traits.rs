//! # TemperatureScaleKernel - Trait Implementations
//!
//! This module contains trait implementations for `TemperatureScaleKernel`.
//!
//! ## Implemented Traits
//!
//! - `ComputeKernel`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)
use crate::compute::ComputeKernel;

use super::types::TemperatureScaleKernel;

impl ComputeKernel for TemperatureScaleKernel {
    fn name(&self) -> &str {
        "TemperatureScaleKernel"
    }
    fn execute(&self, inputs: &[&[f64]], outputs: &mut [Vec<f64>], work_size: usize) {
        if inputs.len() < 3 || outputs.len() < 2 {
            return;
        }
        let vel_flat = inputs[0];
        let masses = inputs[1];
        let params = inputs[2];
        if params.len() < 2 {
            return;
        }
        let t_target = params[0];
        let kb = params[1];
        let n = work_size;
        let ke2: f64 = (0..n)
            .map(|i| {
                let m = if i < masses.len() { masses[i] } else { 1.0 };
                let vx = vel_flat[i * 3];
                let vy = vel_flat[i * 3 + 1];
                let vz = vel_flat[i * 3 + 2];
                m * (vx * vx + vy * vy + vz * vz)
            })
            .sum();
        let n_dof = (3 * n) as f64;
        let t_before = if kb > 1e-30 { ke2 / (n_dof * kb) } else { 0.0 };
        let scale = if t_before > 1e-30 && t_target >= 0.0 {
            (t_target / t_before).sqrt()
        } else {
            1.0
        };
        let mut new_vel = vel_flat.to_vec();
        for i in 0..n {
            new_vel[i * 3] *= scale;
            new_vel[i * 3 + 1] *= scale;
            new_vel[i * 3 + 2] *= scale;
        }
        let ke2_after: f64 = (0..n)
            .map(|i| {
                let m = if i < masses.len() { masses[i] } else { 1.0 };
                let vx = new_vel[i * 3];
                let vy = new_vel[i * 3 + 1];
                let vz = new_vel[i * 3 + 2];
                m * (vx * vx + vy * vy + vz * vz)
            })
            .sum();
        let t_after = if kb > 1e-30 {
            ke2_after / (n_dof * kb)
        } else {
            0.0
        };
        outputs[0] = new_vel;
        outputs[1] = vec![t_before, t_after];
    }
}
