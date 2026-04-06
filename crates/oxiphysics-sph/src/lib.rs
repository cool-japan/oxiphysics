// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Smoothed-particle hydrodynamics for the OxiPhysics engine.
//!
//! This crate provides a complete SPH fluid simulation framework including:
//!
//! - **Kernel functions**: Cubic spline, Wendland C2, Poly6, Spiky
//! - **Neighbor search**: Spatial hash grid for efficient neighbor queries
//! - **Particle data**: SoA-layout particle storage
//! - **WCSPH**: Weakly compressible SPH with Tait equation of state
//! - **DFSPH**: Divergence-free SPH pressure solver
//! - **Surface tension**: Continuum Surface Force model
//! - **Boundary handling**: Penalty planes and Akinci boundary particles
//! - **Time stepping**: CFL-based adaptive time stepping
//! - **Simulation**: High-level simulation driver
#![allow(missing_docs)]
#![allow(dead_code)]

pub mod adaptive;
pub mod adaptive_h;
pub mod adaptive_sph;
pub mod boundary_sph;
pub mod coupling;
pub mod dfsph;
pub mod dfsph_full;
mod error;
pub mod free_surface;
pub mod granular;
pub mod iisph;
pub mod immiscible;
pub mod kernel;
pub mod multiphase;
pub mod neighbor;
pub mod open_boundary;
pub mod particle;
pub mod pcisph;
pub mod pressure_solvers;
pub mod simulation;
pub mod surface;
pub mod surface_tension;
pub mod timestep;
pub mod timestepping;
pub mod turbulence;
pub mod turbulence_sph;
pub mod viscosity;
pub mod wcsph;

pub use adaptive_sph::*;
pub use error::*;
pub use free_surface::*;
pub use iisph::IisphSolver;
pub use multiphase::*;
pub use pcisph::PcisphSolver;
pub use surface_tension::CsfSurfaceTension;
pub use turbulence_sph::*;

/// Trait for SPH fluid solvers.
pub trait SphSolver {
    /// Initialize this component.
    fn init(&mut self);
}
pub mod acoustic_sph;
pub mod adaptive_refinement;
pub mod astrophysical_sph;
pub mod astrophysics_sph;
pub mod climate_sph;
pub mod coastal_sph;
pub mod compressible_sph;
pub mod cryogenic_sph;
pub mod dem_sph;
pub mod dfsph_solver;
pub mod elastic_sph;
pub mod electromagnetic_sph;
pub mod explosion_sph;
pub mod fracture_sph;
pub mod free_surface_sph;
pub mod geothermal_sph;
pub mod ice_sph;
pub mod immersed_boundary_sph;
pub mod incompressible_sph;
pub mod lagrangian_particles;
pub mod lubrication_sph;
pub mod magnetohydro;
pub mod microfluidics_sph;
pub mod mls_sph;
pub mod multiphase_sph;
pub mod multiphase_sph_ext;
pub mod nanoscale_sph;
pub mod neural_sph;
pub mod nuclear_sph;
pub mod ocean_sph;
pub mod plasma_sph;
pub mod quantum_sph;
pub mod reactive_sph;
pub mod sediment;
pub mod sediment_sph;
pub mod sediment_transport;
pub mod shifting;
pub mod soil_sph;
pub mod sph_acoustics;
pub mod sph_analysis;
pub mod sph_boundary;
pub mod sph_coupling;
pub mod sph_multiscale;
pub mod sph_thermal;
pub mod thermal_sph;
pub mod transport_sph;
pub mod variable_resolution;
pub mod viscoelastic_sph;
pub mod volcanic_sph;
pub mod wcsph_solver;
