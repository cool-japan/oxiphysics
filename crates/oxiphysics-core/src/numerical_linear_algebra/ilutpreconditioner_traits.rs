//! # ILUTPreconditioner - Trait Implementations
//!
//! This module contains trait implementations for `ILUTPreconditioner`.
//!
//! ## Implemented Traits
//!
//! - `Preconditioner`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::Preconditioner;
use super::types::ILUTPreconditioner;

impl Preconditioner for ILUTPreconditioner {
    fn apply(&self, r: &[f64]) -> Vec<f64> {
        let n = self.n;
        let mut z = r.to_vec();
        for i in 0..n {
            for &(j, v) in &self.l_rows[i] {
                z[i] -= v * z[j];
            }
        }
        let mut y = z.clone();
        for i in (0..n).rev() {
            for &(j, v) in &self.u_rows[i] {
                if j > i {
                    y[i] -= v * y[j];
                }
            }
            let diag = self.u_rows[i]
                .iter()
                .find(|(j, _)| *j == i)
                .map(|(_, v)| *v)
                .unwrap_or(1.0);
            if diag.abs() > 1e-300 {
                y[i] /= diag;
            }
        }
        y
    }
}
