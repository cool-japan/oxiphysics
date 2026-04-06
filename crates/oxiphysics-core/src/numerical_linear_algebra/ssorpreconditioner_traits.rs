//! # SSORPreconditioner - Trait Implementations
//!
//! This module contains trait implementations for `SSORPreconditioner`.
//!
//! ## Implemented Traits
//!
//! - `Preconditioner`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::Preconditioner;
#[allow(unused_imports)]
use super::functions::*;
use super::types::SSORPreconditioner;

impl Preconditioner for SSORPreconditioner {
    fn apply(&self, r: &[f64]) -> Vec<f64> {
        let n = self.a.nrows;
        let mut z = vec![0.0f64; n];
        for i in 0..n {
            let mut s = self.omega * r[i];
            for k in self.a.row_ptr[i]..self.a.row_ptr[i + 1] {
                let j = self.a.col_idx[k];
                if j < i {
                    s -= self.omega * self.a.values[k] * z[j];
                }
            }
            let mut diag = 1.0;
            for k in self.a.row_ptr[i]..self.a.row_ptr[i + 1] {
                if self.a.col_idx[k] == i {
                    diag = self.a.values[k];
                    break;
                }
            }
            z[i] = if diag.abs() > 1e-300 { s / diag } else { s };
        }
        let mut z2 = z.clone();
        for i in (0..n).rev() {
            let mut s = z[i];
            for k in self.a.row_ptr[i]..self.a.row_ptr[i + 1] {
                let j = self.a.col_idx[k];
                if j > i {
                    s -= self.omega * self.a.values[k] * z2[j];
                }
            }
            let mut diag = 1.0;
            for k in self.a.row_ptr[i]..self.a.row_ptr[i + 1] {
                if self.a.col_idx[k] == i {
                    diag = self.a.values[k];
                    break;
                }
            }
            z2[i] = if diag.abs() > 1e-300 { s / diag } else { s };
        }
        z2
    }
}
