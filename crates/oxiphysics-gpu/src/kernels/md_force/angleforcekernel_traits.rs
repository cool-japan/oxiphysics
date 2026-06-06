//! # AngleForceKernel - Trait Implementations
//!
//! This module contains trait implementations for `AngleForceKernel`.
//!
//! ## Implemented Traits
//!
//! - `ComputeKernel`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)
use crate::compute::ComputeKernel;

use super::types::AngleForceKernel;

impl ComputeKernel for AngleForceKernel {
    fn name(&self) -> &str {
        "AngleForceKernel"
    }
    fn execute(&self, inputs: &[&[f64]], outputs: &mut [Vec<f64>], work_size: usize) {
        if inputs.len() < 2 || outputs.len() < 2 {
            return;
        }
        let pos = inputs[0];
        let angle_data = inputs[1];
        let n = work_size;
        let num_angles = angle_data.len() / 5;
        let mut forces = vec![0.0f64; n * 3];
        let mut energies = vec![0.0f64; num_angles];
        for b in 0..num_angles {
            let i = angle_data[b * 5] as usize;
            let j = angle_data[b * 5 + 1] as usize;
            let k = angle_data[b * 5 + 2] as usize;
            let k_theta = angle_data[b * 5 + 3];
            let theta0 = angle_data[b * 5 + 4];
            if i >= n || j >= n || k >= n {
                continue;
            }
            let rji = [
                pos[i * 3] - pos[j * 3],
                pos[i * 3 + 1] - pos[j * 3 + 1],
                pos[i * 3 + 2] - pos[j * 3 + 2],
            ];
            let rjk = [
                pos[k * 3] - pos[j * 3],
                pos[k * 3 + 1] - pos[j * 3 + 1],
                pos[k * 3 + 2] - pos[j * 3 + 2],
            ];
            let len_ji = (rji[0] * rji[0] + rji[1] * rji[1] + rji[2] * rji[2]).sqrt();
            let len_jk = (rjk[0] * rjk[0] + rjk[1] * rjk[1] + rjk[2] * rjk[2]).sqrt();
            if len_ji < 1e-30 || len_jk < 1e-30 {
                continue;
            }
            let cos_theta = ((rji[0] * rjk[0] + rji[1] * rjk[1] + rji[2] * rjk[2])
                / (len_ji * len_jk))
                .clamp(-1.0, 1.0);
            let theta = cos_theta.acos();
            let delta = theta - theta0;
            energies[b] = 0.5 * k_theta * delta * delta;
            let sin_theta = (1.0 - cos_theta * cos_theta).sqrt().max(1e-12);
            let d_prefactor = -k_theta * delta / sin_theta;
            for dim in 0..3 {
                let d_cos_d_ri =
                    rjk[dim] / (len_ji * len_jk) - cos_theta * rji[dim] / (len_ji * len_ji);
                let d_cos_d_rk =
                    rji[dim] / (len_ji * len_jk) - cos_theta * rjk[dim] / (len_jk * len_jk);
                let fi = d_prefactor * d_cos_d_ri;
                let fk = d_prefactor * d_cos_d_rk;
                forces[i * 3 + dim] += fi;
                forces[k * 3 + dim] += fk;
                forces[j * 3 + dim] -= fi + fk;
            }
        }
        outputs[0] = forces;
        outputs[1] = energies;
    }
}
