//! NetCDF module structure

pub mod functions;
mod functions_2;
mod netcdffile_traits;
#[cfg(test)]
mod trajectory_builder_tests;
pub mod types;

// Re-export all types
pub use functions::*;
pub use types::*;
