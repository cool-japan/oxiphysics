// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! File I/O and serialization for the OxiPhysics engine.
//!
//! Provides writers (and basic readers) for common physics simulation formats:
//! - VTK Legacy (.vtk) and VTU (.vtu) for visualization in ParaView
//! - PDB for molecular data
//! - LAMMPS dump for trajectory data
//! - Wavefront OBJ for triangle meshes
//! - CSV for time series data
//! - XDMF for particle/mesh visualization in ParaView/VisIt
//! - Trajectory for accumulating and writing animation frames
//! - GROMACS GRO format for molecular simulation snapshots
#![allow(missing_docs)]
#![allow(dead_code)]

mod error;
pub use error::*;

pub mod amber;
pub mod binary_format;
pub mod csv;
pub mod csv_io;
pub mod exodus;
pub mod experimental_data_io;
pub mod foam_io;
pub mod gltf;
pub mod gromacs;
pub mod hdf5_simple;
pub mod json_io;
pub mod lammps;
pub mod lammps_dump;
pub mod mesh_quality;
pub mod netcdf;
pub mod numpy;
pub mod obj;
pub mod openfoam;
pub mod pdb;
pub mod stl;
pub mod trajectory;
pub mod vtk;
pub mod vtk_writer;
pub mod vtu;
pub mod xdmf;
pub mod xtc_dcd;

pub use binary_format::*;
pub use csv::{CsvReader, CsvWriter};
pub use foam_io::*;
pub use gromacs::{GroAtom, GroFile};
pub use lammps::{LammpsAtom, LammpsDumpReader, LammpsDumpWriter};
pub use lammps_dump::{
    LammpsDumpFrame, LammpsDumpReader as LammpsDumpFrameReader,
    LammpsDumpWriter as LammpsDumpFrameWriter,
};
pub use obj::{ObjReader, ObjWriter};
pub use pdb::{PdbAtom, PdbReader, PdbWriter};
pub use trajectory::TrajectoryWriter;
pub use vtk::{VtkCellType, VtkDataArray, VtkWriter, VtuGrid};
pub use vtk_writer::*;
pub use vtu::VtuWriter;
pub use xdmf::{write_xdmf_particles, write_xdmf_temporal};

/// Trait for physics data I/O.
#[allow(dead_code)]
pub trait PhysicsIo {
    /// Initialize this component.
    fn init(&mut self);
}
pub mod abaqus_format;
pub mod ambermd_io;
pub mod amr_io;
pub mod animation_io;
pub mod binary_formats;
pub mod binary_io;
pub mod cad_io;
pub mod calculix_format;
pub mod cgns_format;
pub mod checkpoint_io;
pub mod crystallography_io;
pub mod database_io;
pub mod ensight_format;
pub mod exodus_format;
pub mod finite_element_io;
pub mod fluent_format;
pub mod geospatial_io;
pub mod gmsh_format;
pub mod hdf5_io;
pub mod hpc_io;
pub mod lattice_io;
pub mod machine_learning_io;
pub mod material_db;
pub mod material_db_io;
pub mod medical_imaging;
pub mod medical_imaging_io;
pub mod medical_io;
pub mod mesh_export;
pub mod mesh_io;
pub mod molecular_docking_io;
pub mod molecular_visualization_io;
pub mod molecular_viz_io;
pub mod openfoam_format;
pub mod parallel_io;
pub mod particle_data_io;
pub mod particle_formats;
pub mod physics_binary;
pub mod plot3d_format;
pub mod point_cloud_io;
pub mod quantum_chemistry_io;
pub mod remote_sensing_io;
pub mod restart_io;
pub mod robotics_io;
pub mod scientific_formats;
pub mod seismic_io;
pub mod sensor_data_io;
pub mod sensor_io;
pub mod simulation_database;
pub mod simulation_io;
pub mod simulation_log;
pub mod simulation_report_io;
pub mod spectroscopy_io;
pub mod streaming_io;
pub mod su2_format;
pub mod tecplot_format;
pub mod time_series_io;
pub mod visualization_io;
pub mod wavefront_extended;
pub mod weather_data_io;
pub mod weather_io;
pub mod xarray_io;
