//! Auto-generated module structure

pub mod bar;
pub mod functions;
pub mod functions_2;
pub mod mbar_covariance;
#[cfg(test)]
pub(crate) mod test_support;
pub mod types;

// Re-export all types
pub use bar::{bar_uncertainty, fep_bar};
pub use functions::*;
pub use functions_2::*;
pub use mbar_covariance::mbar_df_uncertainty;
pub use types::*;
