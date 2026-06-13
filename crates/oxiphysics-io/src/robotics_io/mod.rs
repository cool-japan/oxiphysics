//! Auto-generated module structure

pub mod functions;
pub mod gripperstate_traits;
pub mod jointlimits_traits;
pub mod jointtype_traits;
pub mod msgtype_traits;
pub mod types;
pub mod urdfinertial_traits;

pub mod articulated_map;
pub mod urdf_error;
pub mod urdf_parser;
pub mod urdf_writer;
pub mod xml;

// Re-export all types
pub use functions::*;
pub use types::*;

pub use articulated_map::urdf_to_articulated;
pub use urdf_error::{UrdfError, UrdfResult};
pub use urdf_parser::{parse_urdf, validate_tree};
pub use urdf_writer::write_urdf;
