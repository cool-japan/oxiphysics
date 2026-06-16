// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Molecular dynamics simulation for the OxiPhysics engine.
//!
//! This crate provides a complete MD simulation framework including:
//!
//! - **Atom data** ([`atom::AtomSet`]): SoA-layout storage for positions,
//!   velocities, forces, masses, charges, and atom types.
//! - **Pair potentials** ([`potential`]): Lennard-Jones, Morse, Coulomb,
//!   and harmonic bond potentials.
//! - **Neighbor lists** ([`neighbor`]): Cell list and Verlet list for
//!   efficient pair interaction computation with periodic boundary conditions.
//! - **Integrators** ([`integrator`]): Velocity Verlet and Leapfrog schemes.
//! - **Force fields** ([`forcefield`]): Pair, bonded, angle, dihedral, and
//!   combined AMBER-like force fields.
//! - **Thermostats** ([`thermostat`]): Berendsen and Nose-Hoover temperature control.
//! - **Barostats** ([`barostat`]): Berendsen and Parrinello-Rahman pressure control.
//! - **ML potentials** ([`ml_potential`], [`nn_potential`]): Machine-learning
//!   force field interface, Behler-Parrinello symmetry function descriptors,
//!   and a simple feedforward neural network potential.
//! - **Simulation driver** ([`simulation::MdSimulation`]): High-level simulation runner.
//! - **Ewald electrostatics** ([`ewald`]): Simplified PME real-space electrostatics
//!   with Ewald damping and self-energy correction.
//! - **Charge data** ([`charge`]): Per-atom charge types and preset charge assignments.
//! - **NHC thermostat** ([`nhc`]): Nosé-Hoover Chain (Martyna *et al.* 1992) with
//!   Yoshida-Suzuki RESPA integration.
//! - **Path integral MD** ([`path_integral`]): Ring polymer representation for PIMD.
//! - **Free energy perturbation** ([`fep`]): TI, Zwanzig, and BAR alchemical estimators.
//! - **Steered MD** ([`steered_md`]): Constant-velocity and constant-force pulling protocols.
#![warn(missing_docs)]

pub mod coarse_grain;
pub mod coarse_grained;
pub mod crystal_growth;
pub mod fitting;
pub use coarse_grained::functions;
pub use coarse_grained::types;
pub mod amber;
pub mod atom;
pub mod autograd_bridge;
pub mod barostat;
pub mod charge;
pub mod constraints;
pub mod electrostatics;
mod error;
pub mod ewald;
pub mod fep;
pub mod forcefield;
pub mod free_energy;
pub mod integrator;
pub mod metadynamics;
pub mod ml_potential;
pub use autograd_bridge::DifferentiableForceField;
pub mod neighbor;
pub mod neighbor_list;
pub mod nhc;
pub mod nn_potential;
pub mod path_integral;
pub mod potential;
pub mod reactive;
pub mod trajectory;
pub use metadynamics::*;
pub mod polarizable;
pub mod protein;
pub mod qmmm;
pub mod remd;
pub mod sampling;
pub mod simulation;
pub mod steered_md;
pub mod structure;
pub mod thermostat;
pub mod water_models;
pub use qmmm::*;
pub mod ab_initio_md;
pub mod adsorption;
pub mod analysis;
pub mod analysis_ext;
pub mod battery_md;
pub mod biomolecular;
pub mod charmm;
pub mod coarse_grain_md;
pub mod coarse_grained_md;
pub mod colloidal_md;
pub mod continuum_coupling;
pub mod crystallography;
pub mod diffusion_md;
pub mod dna_md;
pub mod dna_mechanics;
pub mod dna_simulation;
pub mod electrochemistry_md;
pub mod enhanced_sampling;
pub mod enhanced_sampling_md;
pub mod force_field_builder;
pub mod glass_md;
pub mod ionic_liquid_md;
pub mod lattice_dynamics;
pub mod lipid_bilayer;
pub mod lipid_bilayer_md;
pub mod lipid_membrane_md;
pub mod membrane_sim;
pub mod metal_alloy_md;
pub mod monte_carlo_md;
pub mod nanoparticle_md;
pub mod nemd;
pub mod nucleation;
pub mod photochemistry_md;
pub mod polarizable_md;
pub mod polymer_md;
pub mod polymer_sim;
pub mod protein_folding;
pub mod qm_mm;
pub mod quantum_chemistry;
pub mod quantum_chemistry_md;
pub mod quantum_corrections;
pub mod quantum_md;
pub mod quantum_mechanics;
pub mod quantum_transport_md;
pub mod rare_event;
pub mod reactive_md;
pub mod solid_state_md;
pub mod solvation;
pub mod surface_chemistry_md;
pub mod surface_science;
pub mod surface_science_md;
pub mod thermodynamics_md;
pub mod tribochemistry_md;
pub mod tribology_md;
pub mod zeolite_md;

pub use charge::ChargeData;
pub use error::*;
pub use ewald::{EwaldParams, PmeElectrostatics};
pub use polarizable::*;

// Re-export key types for convenience
pub use atom::AtomSet;
pub use barostat::{Barostat, BerendsenBarostat, NoBarostat, ParrinelloRahmanBarostat};
pub use constraints::{BondConstraint, Shake, water_shake_constraints};
pub use forcefield::{
    AmberForceField, AngleForceField, BondedForceField, DihedralForceField, ForceField,
    PairForceField,
};
pub use integrator::{Integrator, LeapfrogIntegrator, RespaIntegrator, VelocityVerlet};
pub use ml_potential::{MlPotential, SymmetryFunction, SymmetryFunctionSet, cutoff_function};
pub use neighbor::{CellList, PeriodicBox, VerletList, build_verlet_list, distance_pbc};
pub use nemd::{
    MpConfig, MpResult, SllodConfig, SllodResult, gaussian_isokinetic_alpha, run_muller_plathe,
    run_sllod,
};
pub use nn_potential::{Activation, DenseLayer, FeedForwardPotential};
pub use potential::{Coulomb, HarmonicBond, LennardJones, Morse, Potential};
pub use sampling::{UmbrellaSampling, WhamAnalysis};
pub use simulation::{EnergyRecord, MdSimulation};
pub use thermostat::{
    AndersenThermostat, BerendsenThermostat, NoThermostat, NoseHooverThermostat, Thermostat,
    VelocityRescalingThermostat,
};
