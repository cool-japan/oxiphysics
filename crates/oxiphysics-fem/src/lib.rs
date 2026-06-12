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
//!
//! # Examples
//!
//! ## Beam mesh node and element counts
//!
//! ```
//! use oxiphysics_fem::mesh::TetrahedralMesh;
//! let mesh = TetrahedralMesh::generate_beam(1.0, 0.1, 0.1, 2, 1, 1);
//! assert_eq!(mesh.num_nodes(), 3 * 2 * 2);  // (nx+1)*(ny+1)*(nz+1) = 3*2*2 = 12
//! assert_eq!(mesh.num_elements(), 5 * 2 * 1 * 1);  // 5*nx*ny*nz = 10
//! ```
//!
//! ## Material construction and shear modulus
//!
//! ```
//! use oxiphysics_fem::constitutive::LinearElasticMaterial;
//! let mat = LinearElasticMaterial::new(200e9, 0.3);
//! let g = mat.shear_modulus();
//! // G = E / (2*(1+nu)) = 200e9 / 2.6 ≈ 76.92e9
//! assert!((g - 200e9 / (2.0 * 1.3)).abs() < 1e6);
//! ```
//!
//! ## Minimal linear static solve
//!
//! ```
//! use oxiphysics_fem::{
//!     analysis::LinearStaticAnalysis,
//!     boundary::DirichletBc,
//!     constitutive::LinearElasticMaterial,
//!     mesh::TetrahedralMesh,
//! };
//! use oxiphysics_core::math::Vec3;
//!
//! // Tiny 1-cell beam (8 nodes, 5 elements)
//! let mesh = TetrahedralMesh::generate_beam(1.0, 1.0, 1.0, 1, 1, 1);
//! let mat = LinearElasticMaterial::new(1e6, 0.3);
//!
//! // Pin left face (x=0): node_idx = iz*2*2 + iy*2 + ix, ix=0 → [0,2,4,6]
//! let mut bcs = Vec::new();
//! for &node in &[0usize, 2, 4, 6] {
//!     for dof in 0..3 {
//!         bcs.push(DirichletBc::new(node, dof, 0.0));
//!     }
//! }
//!
//! // Gravity body force downward
//! let body_force = Vec3::new(0.0, -1000.0, 0.0);
//! let result = LinearStaticAnalysis::new().solve(&mesh, &mat, &bcs, &[], &body_force);
//! assert_eq!(result.displacements.len(), mesh.num_nodes());
//! // All displacements should be finite
//! for d in &result.displacements {
//!     assert!(d.x.is_finite() && d.y.is_finite() && d.z.is_finite());
//! }
//! ```
#![warn(missing_docs)]

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

pub use damage::{
    CoupledDamagePlasticity, CrackBandModel, DamageBandLocalization, DamageCreepCoupling,
    DamageEvolutionLaw, DamageHomogenization, DamageLawType, DamageMechanicsElement,
    DamagePlasticityState, DamageRateLimiter, DamageState, DamageVisualization,
    ElementDeletionManager, FailureMode, FatigueDamageModel, GursonModel, IsotropicDamage,
    LemaitreCDM, LemaitreChabocheDamage, LemaitreDamage, MazarsDamage, NonLocalDamage,
    NonlocalContinuumDamage, ScalarIsotropicDamage, ThermalDamage,
};
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
