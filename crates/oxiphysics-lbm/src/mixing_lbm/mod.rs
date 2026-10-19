//! Species mixing and passive-scalar transport on a periodic D2Q9 lattice.
//!
//! The lattice solvers ([`MixingLBM`], [`DoubleDistributionFunction`],
//! [`PassiveScalarField`], [`MultiComponentLbm`]) advance BGK populations by
//! collision followed by periodic streaming
//! ([`stream_d2q9_periodic_flat`]), in lattice units (`Δx = Δt = 1`).  A scalar
//! relaxes toward an equilibrium advected with the local velocity and diffuses
//! with `D = c_s² (τ − ½)` ([`lattice_diffusivity`]).
//!
//! The remaining items are closed-form mixing diagnostics (Péclet number,
//! segregation index, Taylor dispersion, ...) and a plain explicit
//! finite-difference diffusion model ([`MixingSimulation`]).

pub mod functions;
pub mod trait_impls;
pub mod types;

#[cfg(test)]
mod transport_tests;

// Re-export all types
pub use functions::*;
pub use types::*;
