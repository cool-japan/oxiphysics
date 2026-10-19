//! NetCDF module structure

mod amber;
#[cfg(test)]
mod amber_tests;
mod cdl_values;
pub mod classic;
pub mod functions;
mod functions_2;
mod netcdffile_traits;
#[cfg(test)]
mod trajectory_builder_tests;
pub mod types;

// Re-export all types
pub use functions::*;
pub use types::*;

pub use amber::{AMBER_VELOCITY_SCALE_FACTOR, AmberTrajectory, read_amber_trajectory};
pub use classic::{
    CdfAttribute, CdfDataset, CdfDimension, CdfVariable, CdfVersion, NcType, NcValues, read_cdf,
    read_cdf_dataset, read_cdf_file, write_cdf, write_cdf_dataset, write_cdf_file,
};
