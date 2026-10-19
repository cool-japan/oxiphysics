//! Auto-generated module structure

pub mod functions;
pub mod functions_2;
pub mod mixed;
pub mod types;
pub mod xdmftimeseries_traits;
pub mod xdmftimeserieshdf5_traits;
pub mod xdmfwriter_traits;

// Re-export all types
pub use functions::*;
pub use functions_2::*;
pub use mixed::xdmf_mixed_element_count;
pub use types::*;

#[cfg(test)]
mod mixed_tests;
