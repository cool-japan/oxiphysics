//! # PairEnergyAccumulateKernel - Trait Implementations
//!
//! This module contains trait implementations for `PairEnergyAccumulateKernel`.
//!
//! ## Implemented Traits
//!
//! - `ComputeKernel`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)
use crate::compute::ComputeKernel;

use super::types::PairEnergyAccumulateKernel;

impl ComputeKernel for PairEnergyAccumulateKernel {
    fn name(&self) -> &str {
        "PairEnergyAccumulateKernel"
    }
    fn execute(&self, inputs: &[&[f64]], outputs: &mut [Vec<f64>], work_size: usize) {
        if inputs.len() < 3 || outputs.len() < 2 {
            return;
        }
        let pos = inputs[0];
        let pairs = inputs[1];
        let epsilon = inputs[2][0];
        let sigma = inputs[2][1];
        let cutoff = inputs[2][2];
        let n = work_size;
        let cutoff2 = cutoff * cutoff;
        let mut per_particle_energy = vec![0.0f64; n];
        let num_pairs = pairs.len() / 2;
        for p in 0..num_pairs {
            let i = pairs[p * 2] as usize;
            let j = pairs[p * 2 + 1] as usize;
            if i >= n || j >= n {
                continue;
            }
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
            let energy = 4.0 * epsilon * (sr12 - sr6);
            per_particle_energy[i] += energy * 0.5;
            per_particle_energy[j] += energy * 0.5;
        }
        let total: f64 = per_particle_energy.iter().sum();
        outputs[0] = per_particle_energy;
        outputs[1] = vec![total];
    }
}
