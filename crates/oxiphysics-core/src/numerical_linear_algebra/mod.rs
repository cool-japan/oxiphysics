//! Auto-generated module structure

pub mod functions;
pub(crate) mod general_eigen;
mod hessenberg_qr;
pub mod ilu0preconditioner_traits;
pub mod ilutpreconditioner_traits;
pub mod jacobipreconditioner_traits;
pub(crate) mod krylov_eigen;
mod lobpcg;
mod matrix_log;
pub mod spaipreconditioner_traits;
pub mod ssorpreconditioner_traits;
pub mod types;

#[cfg(test)]
mod tests;

// Re-export all types
pub use functions::*;
pub use types::*;
