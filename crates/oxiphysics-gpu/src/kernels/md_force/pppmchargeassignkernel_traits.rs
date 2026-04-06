//! # PppmChargeAssignKernel - Trait Implementations
//!
//! This module contains trait implementations for `PppmChargeAssignKernel`.
//!
//! ## Implemented Traits
//!
//! - `ComputeKernel`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use crate::compute::ComputeKernel;

use super::types::PppmChargeAssignKernel;

impl ComputeKernel for PppmChargeAssignKernel {
    fn name(&self) -> &str {
        "PppmChargeAssignKernel"
    }
    fn execute(&self, inputs: &[&[f64]], outputs: &mut [Vec<f64>], work_size: usize) {
        if inputs.len() < 3 || outputs.is_empty() {
            return;
        }
        let pos = inputs[0];
        let charges = inputs[1];
        let grid_params = inputs[2];
        if grid_params.len() < 4 {
            return;
        }
        let nx = grid_params[0] as usize;
        let ny = grid_params[1] as usize;
        let nz = grid_params[2] as usize;
        let box_len = grid_params[3];
        let dx = box_len / nx as f64;
        let dy = box_len / ny as f64;
        let dz = box_len / nz as f64;
        let n = work_size;
        let mut mesh = vec![0.0f64; nx * ny * nz];
        for i in 0..n {
            let x = pos[i * 3].rem_euclid(box_len);
            let y = pos[i * 3 + 1].rem_euclid(box_len);
            let z = pos[i * 3 + 2].rem_euclid(box_len);
            let ix = ((x / dx) as usize).min(nx - 1);
            let iy = ((y / dy) as usize).min(ny - 1);
            let iz = ((z / dz) as usize).min(nz - 1);
            let grid_idx = ix * ny * nz + iy * nz + iz;
            mesh[grid_idx] += charges[i];
        }
        outputs[0] = mesh;
    }
}
