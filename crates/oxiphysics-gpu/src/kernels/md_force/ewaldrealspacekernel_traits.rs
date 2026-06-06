//! # EwaldRealSpaceKernel - Trait Implementations
//!
//! This module contains trait implementations for `EwaldRealSpaceKernel`.
//!
//! ## Implemented Traits
//!
//! - `ComputeKernel`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::compute::ComputeKernel;

use super::functions::erfc_approx;
use super::types::EwaldRealSpaceKernel;

impl ComputeKernel for EwaldRealSpaceKernel {
    fn name(&self) -> &str {
        "EwaldRealSpaceKernel"
    }
    fn execute(&self, inputs: &[&[f64]], outputs: &mut [Vec<f64>], work_size: usize) {
        if inputs.len() < 3 || outputs.len() < 2 {
            return;
        }
        let pos = inputs[0];
        let charges = inputs[1];
        let alpha = inputs[2][0];
        let r_cutoff = inputs[2][1];
        let box_len = inputs[2][2];
        let n = work_size;
        let r_cut2 = r_cutoff * r_cutoff;
        let half_box = box_len * 0.5;
        let mut forces = vec![0.0f64; n * 3];
        let mut energy = 0.0f64;
        for i in 0..n {
            let xi = [pos[i * 3], pos[i * 3 + 1], pos[i * 3 + 2]];
            let qi = charges[i];
            for j in (i + 1)..n {
                let xj = [pos[j * 3], pos[j * 3 + 1], pos[j * 3 + 2]];
                let qj = charges[j];
                let mut dx = xi[0] - xj[0];
                let mut dy = xi[1] - xj[1];
                let mut dz = xi[2] - xj[2];
                if dx > half_box {
                    dx -= box_len;
                } else if dx < -half_box {
                    dx += box_len;
                }
                if dy > half_box {
                    dy -= box_len;
                } else if dy < -half_box {
                    dy += box_len;
                }
                if dz > half_box {
                    dz -= box_len;
                } else if dz < -half_box {
                    dz += box_len;
                }
                let r2 = dx * dx + dy * dy + dz * dz;
                if r2 >= r_cut2 || r2 < 1e-30 {
                    continue;
                }
                let r = r2.sqrt();
                let ar = alpha * r;
                let erfc_ar = erfc_approx(ar);
                energy += qi * qj * erfc_ar / r;
                let two_alpha_over_sqrt_pi = 2.0 * alpha / std::f64::consts::PI.sqrt();
                let deriv = -(erfc_ar / r + two_alpha_over_sqrt_pi * (-ar * ar).exp()) / r;
                let f_mag = -qi * qj * deriv / r;
                forces[i * 3] += f_mag * dx;
                forces[i * 3 + 1] += f_mag * dy;
                forces[i * 3 + 2] += f_mag * dz;
                forces[j * 3] -= f_mag * dx;
                forces[j * 3 + 1] -= f_mag * dy;
                forces[j * 3 + 2] -= f_mag * dz;
            }
        }
        outputs[0] = forces;
        outputs[1] = vec![energy];
    }
}
