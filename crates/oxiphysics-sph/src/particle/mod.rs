//! Auto-generated module structure

pub mod functions;
mod pair_sums;
pub mod particleset_traits;
pub mod sphparticleset_traits;
pub mod types;

// Re-export all types
pub use functions::*;
pub use types::*;

#[cfg(test)]
mod pair_sum_tests;
