//! Auto-generated module structure

pub mod functions;
mod iges_parser;
pub mod igesreader_traits;
pub mod stepparser_traits;
mod tessellation;
pub mod trianglemesh_traits;
pub mod types;

// Re-export all types
pub use functions::*;
pub use types::*;
