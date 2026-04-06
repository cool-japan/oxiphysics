//! # NlistUpdateKernel - Trait Implementations
//!
//! This module contains trait implementations for `NlistUpdateKernel`.
//!
//! ## Implemented Traits
//!
//! - `ComputeKernel`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use crate::compute::ComputeKernel;

use super::types::NlistUpdateKernel;

impl ComputeKernel for NlistUpdateKernel {
    fn name(&self) -> &str {
        "NlistUpdateKernel"
    }
    fn execute(&self, inputs: &[&[f64]], outputs: &mut [Vec<f64>], work_size: usize) {
        if inputs.len() < 3 || outputs.len() < 2 {
            return;
        }
        let pos = inputs[0];
        let ref_pos = inputs[1];
        let cutoff = inputs[2][0];
        let skin = inputs[2][1];
        let n = work_size;
        let mut max_disp2 = 0.0f64;
        for i in 0..n {
            let dx = pos[i * 3] - ref_pos[i * 3];
            let dy = pos[i * 3 + 1] - ref_pos[i * 3 + 1];
            let dz = pos[i * 3 + 2] - ref_pos[i * 3 + 2];
            let d2 = dx * dx + dy * dy + dz * dz;
            if d2 > max_disp2 {
                max_disp2 = d2;
            }
        }
        let max_disp = max_disp2.sqrt();
        let needs_rebuild = max_disp > skin * 0.5;
        if !needs_rebuild {
            outputs[0] = vec![];
            outputs[1] = vec![0.0, 0.0];
            return;
        }
        let r_list = cutoff + skin;
        let r_list2 = r_list * r_list;
        let mut pairs = Vec::new();
        for i in 0..n {
            for j in (i + 1)..n {
                let dx = pos[i * 3] - pos[j * 3];
                let dy = pos[i * 3 + 1] - pos[j * 3 + 1];
                let dz = pos[i * 3 + 2] - pos[j * 3 + 2];
                let r2 = dx * dx + dy * dy + dz * dz;
                if r2 < r_list2 {
                    pairs.push(i as f64);
                    pairs.push(j as f64);
                }
            }
        }
        let num_pairs = pairs.len() / 2;
        outputs[0] = pairs;
        outputs[1] = vec![1.0, num_pairs as f64];
    }
}
