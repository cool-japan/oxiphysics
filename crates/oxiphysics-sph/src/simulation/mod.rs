//! Auto-generated module structure

pub mod functions;
pub mod multiphaseparticleset_traits;
pub mod sph_sim_neighbors;
pub mod sphsimstate_traits;
pub mod sphsimulationparams_traits;
pub mod types;
pub mod types_sim;

// Re-export all types (trait impl modules don't export names, just impl blocks)
pub use sph_sim_neighbors::*;
pub use types::*;
pub use types_sim::*;

#[cfg(test)]
mod neighbor_tests;
