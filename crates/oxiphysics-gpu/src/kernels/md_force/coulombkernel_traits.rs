//! # CoulombKernel - Trait Implementations
//!
//! This module contains trait implementations for `CoulombKernel`.
//!
//! ## Implemented Traits
//!
//! - `ComputeKernel`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)
use crate::compute::ComputeKernel;

use super::types::CoulombKernel;

impl ComputeKernel for CoulombKernel {
    fn name(&self) -> &str {
        "CoulombKernel"
    }
    fn execute(&self, inputs: &[&[f64]], outputs: &mut [Vec<f64>], work_size: usize) {
        if inputs.len() < 3 || outputs.len() < 2 {
            return;
        }
        let pos = inputs[0];
        let charges = inputs[1];
        let k_e = inputs[2][0];
        let cutoff = inputs[2][1];
        let n = work_size;
        let cutoff2 = cutoff * cutoff;
        let mut forces = vec![0.0; n * 3];
        let mut potential = 0.0;
        for i in 0..n {
            let xi = [pos[i * 3], pos[i * 3 + 1], pos[i * 3 + 2]];
            let qi = charges[i];
            for j in (i + 1)..n {
                let xj = [pos[j * 3], pos[j * 3 + 1], pos[j * 3 + 2]];
                let dx = xi[0] - xj[0];
                let dy = xi[1] - xj[1];
                let dz = xi[2] - xj[2];
                let r2 = dx * dx + dy * dy + dz * dz;
                if r2 >= cutoff2 || r2 < 1e-30 {
                    continue;
                }
                let r = r2.sqrt();
                let qj = charges[j];
                potential += k_e * qi * qj / r;
                let f_mag = k_e * qi * qj / (r2 * r);
                forces[i * 3] += f_mag * dx;
                forces[i * 3 + 1] += f_mag * dy;
                forces[i * 3 + 2] += f_mag * dz;
                forces[j * 3] -= f_mag * dx;
                forces[j * 3 + 1] -= f_mag * dy;
                forces[j * 3 + 2] -= f_mag * dz;
            }
        }
        outputs[0] = forces;
        outputs[1] = vec![potential];
    }
}
