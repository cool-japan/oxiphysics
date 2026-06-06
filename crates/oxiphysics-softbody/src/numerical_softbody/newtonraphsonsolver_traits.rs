//! # NewtonRaphsonSolver - Trait Implementations
//!
//! This module contains trait implementations for `NewtonRaphsonSolver`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{LineSearch, NewtonRaphsonSolver, PcgParams};

impl Default for NewtonRaphsonSolver {
    fn default() -> Self {
        NewtonRaphsonSolver {
            max_iter: 20,
            tolerance: 1e-6,
            pcg_params: PcgParams::default(),
            line_search: LineSearch::default(),
        }
    }
}
