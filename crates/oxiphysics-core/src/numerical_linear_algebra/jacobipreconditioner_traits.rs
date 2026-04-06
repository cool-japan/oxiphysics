//! # JacobiPreconditioner - Trait Implementations
//!
//! This module contains trait implementations for `JacobiPreconditioner`.
//!
//! ## Implemented Traits
//!
//! - `Preconditioner`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::Preconditioner;
#[allow(unused_imports)]
use super::functions::*;
use super::types::JacobiPreconditioner;

impl Preconditioner for JacobiPreconditioner {
    fn apply(&self, r: &[f64]) -> Vec<f64> {
        r.iter()
            .zip(&self.inv_diag)
            .map(|(ri, di)| ri * di)
            .collect()
    }
}
