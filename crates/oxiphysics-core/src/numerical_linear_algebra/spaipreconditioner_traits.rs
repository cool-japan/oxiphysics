//! # SPAIPreconditioner - Trait Implementations
//!
//! This module contains trait implementations for `SPAIPreconditioner`.
//!
//! ## Implemented Traits
//!
//! - `Preconditioner`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::Preconditioner;
use super::types::SPAIPreconditioner;

impl Preconditioner for SPAIPreconditioner {
    fn apply(&self, r: &[f64]) -> Vec<f64> {
        let mut z = vec![0.0f64; self.n];
        for (j, col) in self.columns.iter().enumerate() {
            for &(i, v) in col {
                z[i] += v * r[j];
            }
        }
        z
    }
}
