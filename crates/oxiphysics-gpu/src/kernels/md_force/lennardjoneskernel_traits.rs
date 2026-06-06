//! # LennardJonesKernel - Trait Implementations
//!
//! This module contains trait implementations for `LennardJonesKernel`.
//!
//! ## Implemented Traits
//!
//! - `ComputeKernel`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)
use crate::compute::ComputeKernel;

use super::types::LennardJonesKernel;

impl ComputeKernel for LennardJonesKernel {
    fn name(&self) -> &str {
        "LennardJonesKernel"
    }
    fn execute(&self, inputs: &[&[f64]], outputs: &mut [Vec<f64>], work_size: usize) {
        if inputs.len() < 2 || outputs.len() < 2 {
            return;
        }
        let pos = inputs[0];
        let epsilon = inputs[1][0];
        let sigma = inputs[1][1];
        let cutoff = inputs[1][2];
        let n = work_size;
        let cutoff2 = cutoff * cutoff;
        let mut forces = vec![0.0; n * 3];
        let mut potential = 0.0;
        for i in 0..n {
            let xi = [pos[i * 3], pos[i * 3 + 1], pos[i * 3 + 2]];
            for j in (i + 1)..n {
                let xj = [pos[j * 3], pos[j * 3 + 1], pos[j * 3 + 2]];
                let dx = xi[0] - xj[0];
                let dy = xi[1] - xj[1];
                let dz = xi[2] - xj[2];
                let r2 = dx * dx + dy * dy + dz * dz;
                if r2 >= cutoff2 || r2 < 1e-30 {
                    continue;
                }
                let r2_inv = 1.0 / r2;
                let sr2 = sigma * sigma * r2_inv;
                let sr6 = sr2 * sr2 * sr2;
                let sr12 = sr6 * sr6;
                potential += 4.0 * epsilon * (sr12 - sr6);
                let f_mag = 24.0 * epsilon * (2.0 * sr12 - sr6) * r2_inv;
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
