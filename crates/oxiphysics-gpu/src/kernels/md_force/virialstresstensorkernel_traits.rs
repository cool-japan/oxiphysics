//! # VirialStressTensorKernel - Trait Implementations
//!
//! This module contains trait implementations for `VirialStressTensorKernel`.
//!
//! ## Implemented Traits
//!
//! - `ComputeKernel`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)
use crate::compute::ComputeKernel;

use super::types::VirialStressTensorKernel;

impl ComputeKernel for VirialStressTensorKernel {
    fn name(&self) -> &str {
        "VirialStressTensorKernel"
    }
    fn execute(&self, inputs: &[&[f64]], outputs: &mut [Vec<f64>], work_size: usize) {
        if inputs.len() < 2 || outputs.is_empty() {
            return;
        }
        let pos = inputs[0];
        let epsilon = inputs[1][0];
        let sigma = inputs[1][1];
        let cutoff = inputs[1][2];
        let n = work_size;
        let cutoff2 = cutoff * cutoff;
        let mut w = [0.0f64; 6];
        for i in 0..n {
            for j in (i + 1)..n {
                let dx = pos[i * 3] - pos[j * 3];
                let dy = pos[i * 3 + 1] - pos[j * 3 + 1];
                let dz = pos[i * 3 + 2] - pos[j * 3 + 2];
                let r2 = dx * dx + dy * dy + dz * dz;
                if r2 >= cutoff2 || r2 < 1e-30 {
                    continue;
                }
                let r2_inv = 1.0 / r2;
                let sr2 = sigma * sigma * r2_inv;
                let sr6 = sr2 * sr2 * sr2;
                let sr12 = sr6 * sr6;
                let f_over_r2 = 24.0 * epsilon * (2.0 * sr12 - sr6) * r2_inv;
                let r_vec = [dx, dy, dz];
                let c = [(0usize, 0usize), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)];
                for (ci, &(a, b)) in c.iter().enumerate() {
                    w[ci] -= r_vec[a] * f_over_r2 * r_vec[b];
                }
            }
        }
        outputs[0] = w.to_vec();
    }
}
