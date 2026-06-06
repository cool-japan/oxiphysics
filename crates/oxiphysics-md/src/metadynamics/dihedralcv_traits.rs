//! # DihedralCV - Trait Implementations
//!
//! This module contains trait implementations for `DihedralCV`.
//!
//! ## Implemented Traits
//!
//! - `CollectiveVariable`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::CollectiveVariable;
use super::types::DihedralCV;

impl CollectiveVariable for DihedralCV {
    fn value(&self, positions: &[[f64; 3]]) -> f64 {
        Self::compute_dihedral(positions, self.i, self.j, self.k, self.l)
    }
    fn gradient(&self, positions: &[[f64; 3]]) -> Vec<[f64; 3]> {
        let n = positions.len();
        let mut grad = vec![[0.0f64; 3]; n];
        let eps = 1e-5;
        for &idx in &[self.i, self.j, self.k, self.l] {
            for dim in 0..3 {
                let mut pos_plus = positions.to_vec();
                let mut pos_minus = positions.to_vec();
                pos_plus[idx][dim] += eps;
                pos_minus[idx][dim] -= eps;
                let phi_plus = Self::compute_dihedral(&pos_plus, self.i, self.j, self.k, self.l);
                let phi_minus = Self::compute_dihedral(&pos_minus, self.i, self.j, self.k, self.l);
                let mut dphi = phi_plus - phi_minus;
                if dphi > std::f64::consts::PI {
                    dphi -= 2.0 * std::f64::consts::PI;
                } else if dphi < -std::f64::consts::PI {
                    dphi += 2.0 * std::f64::consts::PI;
                }
                grad[idx][dim] = dphi / (2.0 * eps);
            }
        }
        grad
    }
    fn name(&self) -> &str {
        &self.name
    }
}
