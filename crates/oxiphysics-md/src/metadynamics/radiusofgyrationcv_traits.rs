//! # RadiusOfGyrationCV - Trait Implementations
//!
//! This module contains trait implementations for `RadiusOfGyrationCV`.
//!
//! ## Implemented Traits
//!
//! - `CollectiveVariable`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::CollectiveVariable;
use super::types::RadiusOfGyrationCV;

impl CollectiveVariable for RadiusOfGyrationCV {
    fn value(&self, positions: &[[f64; 3]]) -> f64 {
        let n = self.atom_indices.len();
        if n == 0 {
            return 0.0;
        }
        let mut com = [0.0f64; 3];
        for &idx in &self.atom_indices {
            com[0] += positions[idx][0];
            com[1] += positions[idx][1];
            com[2] += positions[idx][2];
        }
        let nf = n as f64;
        com[0] /= nf;
        com[1] /= nf;
        com[2] /= nf;
        let sum_sq: f64 = self
            .atom_indices
            .iter()
            .map(|&idx| {
                let dx = positions[idx][0] - com[0];
                let dy = positions[idx][1] - com[1];
                let dz = positions[idx][2] - com[2];
                dx * dx + dy * dy + dz * dz
            })
            .sum();
        (sum_sq / nf).sqrt()
    }
    fn gradient(&self, positions: &[[f64; 3]]) -> Vec<[f64; 3]> {
        let n_atoms = positions.len();
        let n = self.atom_indices.len();
        let mut grad = vec![[0.0f64; 3]; n_atoms];
        if n == 0 {
            return grad;
        }
        let rg = self.value(positions);
        if rg < 1e-15 {
            return grad;
        }
        let mut com = [0.0f64; 3];
        for &idx in &self.atom_indices {
            com[0] += positions[idx][0];
            com[1] += positions[idx][1];
            com[2] += positions[idx][2];
        }
        let nf = n as f64;
        com[0] /= nf;
        com[1] /= nf;
        com[2] /= nf;
        for &idx in &self.atom_indices {
            let dx = positions[idx][0] - com[0];
            let dy = positions[idx][1] - com[1];
            let dz = positions[idx][2] - com[2];
            grad[idx] = [dx / (nf * rg), dy / (nf * rg), dz / (nf * rg)];
        }
        grad
    }
    fn name(&self) -> &str {
        &self.name
    }
}
