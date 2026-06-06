//! # ILU0Preconditioner - Trait Implementations
//!
//! This module contains trait implementations for `ILU0Preconditioner`.
//!
//! ## Implemented Traits
//!
//! - `Preconditioner`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::Preconditioner;
use super::types::ILU0Preconditioner;

impl Preconditioner for ILU0Preconditioner {
    fn apply(&self, r: &[f64]) -> Vec<f64> {
        let n = self.n;
        let mut z = vec![0.0f64; n];
        for i in 0..n {
            z[i] = r[i];
            for k in self.row_ptr[i]..self.row_ptr[i + 1] {
                let j = self.col_idx[k];
                if j < i {
                    z[i] -= self.lu_vals[k] * z[j];
                }
            }
        }
        let mut y = vec![0.0f64; n];
        for i in (0..n).rev() {
            y[i] = z[i];
            for k in self.row_ptr[i]..self.row_ptr[i + 1] {
                let j = self.col_idx[k];
                if j > i {
                    y[i] -= self.lu_vals[k] * y[j];
                }
            }
            for k in self.row_ptr[i]..self.row_ptr[i + 1] {
                if self.col_idx[k] == i {
                    if self.lu_vals[k].abs() > 1e-300 {
                        y[i] /= self.lu_vals[k];
                    }
                    break;
                }
            }
        }
        y
    }
}
