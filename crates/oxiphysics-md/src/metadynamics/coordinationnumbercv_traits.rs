//! # CoordinationNumberCV - Trait Implementations
//!
//! This module contains trait implementations for `CoordinationNumberCV`.
//!
//! ## Implemented Traits
//!
//! - `CollectiveVariable`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::CollectiveVariable;
use super::types::CoordinationNumberCV;

impl CollectiveVariable for CoordinationNumberCV {
    fn value(&self, positions: &[[f64; 3]]) -> f64 {
        self.neighbours
            .iter()
            .map(|&j| {
                let ri = positions[self.atom_ref];
                let rj = positions[j];
                let r2 =
                    (ri[0] - rj[0]).powi(2) + (ri[1] - rj[1]).powi(2) + (ri[2] - rj[2]).powi(2);
                self.switching_function(r2.sqrt())
            })
            .sum()
    }
    /// Analytic gradient `d s / d r_a` of the coordination number.
    ///
    /// Each neighbour `j` contributes `g = sigma'(r) (r_j - r_ref) / r` to atom
    /// `j` and `-g` to the reference atom, so the gradients of all atoms sum to
    /// zero (translational invariance).  A neighbour that coincides with the
    /// reference atom has no defined direction and contributes nothing.
    fn gradient(&self, positions: &[[f64; 3]]) -> Vec<[f64; 3]> {
        let mut grad = vec![[0.0f64; 3]; positions.len()];
        let r_ref = positions[self.atom_ref];
        for &j in &self.neighbours {
            let r_j = positions[j];
            let d = [r_j[0] - r_ref[0], r_j[1] - r_ref[1], r_j[2] - r_ref[2]];
            let r = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
            if r <= 0.0 {
                continue;
            }
            let ds_dr = self.switching_derivative(r);
            for dim in 0..3 {
                let g = ds_dr * d[dim] / r;
                grad[j][dim] += g;
                grad[self.atom_ref][dim] -= g;
            }
        }
        grad
    }
    fn name(&self) -> &str {
        &self.name
    }
}
