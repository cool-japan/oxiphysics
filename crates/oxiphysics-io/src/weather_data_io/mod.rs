//! Auto-generated module structure

pub mod functions;
mod grib2_decode;
pub mod netcdfcfreader_traits;
pub mod precipitationanalysis_traits;
pub mod types;
pub mod weatherwriter_traits;

// Re-export all types
pub use functions::*;
pub use types::*;

#[cfg(test)]
mod grib2_tests;
