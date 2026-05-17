//! Formal verification methods for physics simulation correctness.
//!
//! This module provides model checking, abstract interpretation, bisimulation,
//! propositional SAT solving, and conservation-law verifiers for physical
//! simulation states.

pub mod functions;
pub(crate) mod ordered_float_shim;
pub mod trait_impls;
pub mod types;

// Re-export all types
pub use functions::*;
pub use types::*;
