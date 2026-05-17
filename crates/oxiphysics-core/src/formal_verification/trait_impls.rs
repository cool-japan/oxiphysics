//! # AbstractInterpretation - Trait Implementations
//!
//! This module contains trait implementations for `AbstractInterpretation`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//! - `Default`
//! - `Display`
//! - `Debug`
//! - `Debug`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{
    AbstractInterpretation, CegarConfig, IntervalDomain, Invariant, RefinementStrategy,
    SafetyProperty,
};

impl Default for AbstractInterpretation {
    fn default() -> Self {
        Self::new()
    }
}

impl Default for CegarConfig {
    fn default() -> Self {
        Self {
            max_iterations: 20,
            bmc_bound: 10,
            timeout_ms: 0,
            refinement_strategy: RefinementStrategy::SyntaxGuided,
        }
    }
}

impl std::fmt::Display for IntervalDomain {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "[{}, {}]", self.lo, self.hi)
    }
}

impl std::fmt::Debug for Invariant {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Invariant")
            .field("name", &self.name)
            .finish()
    }
}

impl std::fmt::Debug for SafetyProperty {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("SafetyProperty")
            .field("label", &self.label)
            .finish_non_exhaustive()
    }
}
