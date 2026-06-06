//! # ParallelIslandSolver - Trait Implementations
//!
//! This module contains trait implementations for `ParallelIslandSolver`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{ParallelIslandSolver, SequentialImpulseSolver};

impl Default for ParallelIslandSolver {
    fn default() -> Self {
        Self {
            base_solver: SequentialImpulseSolver::default(),
            max_substeps: 8,
        }
    }
}
