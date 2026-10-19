//! Auto-generated module structure

pub mod functions;
mod nifti_writer;
mod nrrd_decode;
pub mod nrrdreader_traits;
pub mod types;

// Re-export all types
pub use functions::*;
pub use nifti_writer::NiftiWriter;
pub use types::*;
