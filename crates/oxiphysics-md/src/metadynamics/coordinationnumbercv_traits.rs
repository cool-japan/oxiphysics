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
#[allow(unused_imports)]
use super::functions::*;
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
    fn gradient(&self, positions: &[[f64; 3]]) -> Vec<[f64; 3]> {
        let n_atoms = positions.len();
        let mut grad = vec![[0.0f64; 3]; n_atoms];
        let eps = 1e-5;
        for &j in &self.neighbours {
            for dim in 0..3 {
                let mut pp = positions.to_vec();
                let mut pm = positions.to_vec();
                pp[self.atom_ref][dim] += eps;
                pm[self.atom_ref][dim] -= eps;
                let dp = self.value(&pp);
                let dm = self.value(&pm);
                grad[self.atom_ref][dim] += (dp - dm) / (2.0 * eps);
                let mut pj = positions.to_vec();
                let mut mj = positions.to_vec();
                pj[j][dim] += eps;
                mj[j][dim] -= eps;
                let dpj = self.value(&pj);
                let dmj = self.value(&mj);
                grad[j][dim] += (dpj - dmj) / (2.0 * eps);
            }
        }
        grad
    }
    fn name(&self) -> &str {
        &self.name
    }
}
