// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Finite element method for the OxiPhysics engine.
//!
//! Provides a complete linear elastic FEM pipeline including:
//!
//! - **Sparse matrix** types ([`sparse::CsrMatrix`], [`sparse::SparseVector`])
//! - **Mesh** representation ([`mesh::TetrahedralMesh`]) with beam generation
//! - **Element** formulation ([`element::LinearTetrahedron`]) for 4-node tetrahedra
//! - **Constitutive** models ([`constitutive::LinearElasticMaterial`])
//! - **Assembly** of global stiffness and load vectors
//! - **Boundary conditions** (Dirichlet and Neumann)
//! - **Iterative solvers** (CG, preconditioned CG)
//! - **Static analysis** driver ([`analysis::LinearStaticAnalysis`])
#![warn(missing_docs)]
#![allow(ambiguous_glob_reexports)]

pub mod adaptive_mesh;
pub mod analysis;
pub mod assembly;
pub mod boundary;
pub mod buckling;
pub mod cohesive;
pub mod composite_fem;
pub mod constitutive;
pub mod contact;
pub mod contact_fem;
pub mod contact_mechanics;
pub mod coupled_fem;
pub mod coupled_physics_fem;
pub mod damage;
pub mod dynamic;
pub mod dynamics;
pub mod eigen;
pub mod eigenvalue;
pub mod electrochemical_fem;
pub mod electromagnetic;
pub mod electromechanics;
pub mod element;
mod error;
pub mod error_estimation;
pub mod error_estimator;
pub mod fluid_fem;
pub mod fluid_structure;
pub mod fracture;
pub mod geomechanics_fem;
pub mod homogenization;
pub mod hyperelastic;
pub mod inverse_problems;
pub mod isogeometric;
pub mod lbm_coupling;
pub mod level_set;
pub mod mesh;
pub mod meshless;
pub mod mixed_elements;
pub mod modal;
pub mod multiphysics_fem;
pub mod multiscale;
pub mod multiscale_fem;
pub mod nonlinear;
pub mod nonlinear_dynamics;
pub mod nonlinear_fem;
pub mod nonlinear_solver;
pub mod optimization_fem;
pub mod piezo;
pub mod poromechanics;
pub mod reduced_order;
pub mod reliability;
pub mod reliability_fem;
pub mod shell;
pub mod solvers;
pub mod sparse;
pub mod spectral_fem;
pub mod stochastic_fem;
pub mod thermal;
pub mod thermal_fem;
pub mod thermal_stress;
pub mod topology_opt;
pub mod topology_optimization;
pub mod viscoelastic_fem;
pub mod wave_propagation;
pub mod xfem;

pub mod parallel_solver;
pub use parallel_solver::{
    CsrMatrix as ParCsrMatrix, GmresWithAmg, ParallelAssembler, ParallelGmresSolver,
    ParallelPcgSolver, PcgWithAmg, PcgWithPrecond,
};

pub mod perf_bench;
pub use perf_bench::{
    BenchHarness, BenchReport, BenchResult, SuiteConfig, banded_csr, bench_assembly,
    bench_element_stiffness, bench_gmres, bench_pcg, bench_spmv, bench_spmv_parallel,
    run_full_suite, run_suite, tridiagonal_csr,
};

pub use damage::*;
pub use electromechanics::*;
pub use error::*;
pub use homogenization::*;
pub use modal::{ModalResult, inverse_iteration};
pub use nonlinear::{ArcLengthControl, ConvergenceCriteria, NrResult, newton_raphson};
pub use nonlinear_solver::*;

/// Trait for finite element solvers.
pub trait FemSolver {
    /// Initialize this component.
    fn init(&mut self);
}
pub mod adaptive_fem;
pub mod additive_manufacturing_fem;
pub mod beam_fem;
pub mod biomechanics_fem;
pub mod boundary_element;
pub mod contact_mech;
pub mod crystal_plasticity;
pub mod damage_mechanics;
pub mod discontinuous_galerkin;
pub mod explicit_fem;
pub mod fatigue_fem;
pub mod isogeometric_fem;
pub mod meshfree_fem;
pub mod probabilistic_fem;
pub mod reduced_order_fem;
pub mod shell_fem;
pub mod soil_fem;
pub mod topology_opt_fem;
pub mod topology_optimization_ext;
pub mod truss_frame;
