//! Auto-generated module structure

pub mod correlation;
#[cfg(test)]
mod dft_tests;
pub mod error;
pub mod functions;
pub mod functions_2;
pub mod integrals;
mod kohn_sham;
pub mod molecular_grid;
pub mod sto3g_basis;
pub mod types;
pub mod types_2;
pub mod types_3;
pub mod xc_functional;

// Re-export all types
pub use correlation::*;
pub use error::QuantumChemistryError;
pub use functions::*;
pub use integrals::{
    dipole_integral, electron_repulsion_integral, kinetic_integral, nuclear_attraction_integral,
    overlap_integral,
};
pub use molecular_grid::{MolecularGrid, MolecularGridSettings};
pub use sto3g_basis::{sto3g_atom_basis, sto3g_supports_element};
pub use types::*;
pub use types_2::*;
pub use types_3::*;
pub use xc_functional::{
    XcPointValue, evaluate_xc, uniform_gas_tau, xc_needs_gradient, xc_needs_tau,
};
