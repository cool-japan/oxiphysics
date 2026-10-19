//! Auto-generated module structure

#[cfg(test)]
mod chgcar_tests;
mod cif;
#[cfg(test)]
mod cif_tests;
mod dx;
#[cfg(test)]
mod dx_tests;
pub mod functions;
pub mod molvizerror_traits;
pub mod types;

// Re-export all types
pub use cif::{parse_cif, parse_cif_number};
pub use functions::*;
pub use types::*;
