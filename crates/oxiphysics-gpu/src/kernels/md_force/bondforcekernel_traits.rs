//! # BondForceKernel - Trait Implementations
//!
//! This module contains trait implementations for `BondForceKernel`.
//!
//! ## Implemented Traits
//!
//! - `ComputeKernel`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)
use crate::compute::ComputeKernel;

use super::types::BondForceKernel;

impl ComputeKernel for BondForceKernel {
    fn name(&self) -> &str {
        "BondForceKernel"
    }
    fn execute(&self, inputs: &[&[f64]], outputs: &mut [Vec<f64>], work_size: usize) {
        if inputs.len() < 2 || outputs.len() < 2 {
            return;
        }
        let pos = inputs[0];
        let bond_data = inputs[1];
        let n = work_size;
        let num_bonds = bond_data.len() / 4;
        let mut forces = vec![0.0f64; n * 3];
        let mut energies = vec![0.0f64; num_bonds];
        for b in 0..num_bonds {
            let i = bond_data[b * 4] as usize;
            let j = bond_data[b * 4 + 1] as usize;
            let k_spring = bond_data[b * 4 + 2];
            let r0 = bond_data[b * 4 + 3];
            if i >= n || j >= n {
                continue;
            }
            let dx = pos[j * 3] - pos[i * 3];
            let dy = pos[j * 3 + 1] - pos[i * 3 + 1];
            let dz = pos[j * 3 + 2] - pos[i * 3 + 2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            if r < 1e-30 {
                continue;
            }
            let delta = r - r0;
            energies[b] = 0.5 * k_spring * delta * delta;
            let mag = k_spring * delta / r;
            forces[i * 3] += mag * dx;
            forces[i * 3 + 1] += mag * dy;
            forces[i * 3 + 2] += mag * dz;
            forces[j * 3] -= mag * dx;
            forces[j * 3 + 1] -= mag * dy;
            forces[j * 3 + 2] -= mag * dz;
        }
        outputs[0] = forces;
        outputs[1] = energies;
    }
}
