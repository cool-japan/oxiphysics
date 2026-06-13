//! Auto-generated module structure

mod functions;
mod functions_2;
pub mod mpr;
pub mod gjksolver_traits;
pub mod simplex_traits;
mod types;
pub mod warmstartgjk_traits;

// Re-export all types
pub use functions::*;
pub use mpr::{mpr_contact, mpr_full};
pub use types::*;
