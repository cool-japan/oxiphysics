//! XPBD constraints and strain-limiting utilities.

pub mod functions;
pub mod strain_limit;
pub mod types;

// Re-export all types
pub use functions::*;
pub use strain_limit::{apply_strain_limiting, StrainLimitConfig};
pub use types::*;
