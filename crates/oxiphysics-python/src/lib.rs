// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Python-friendly wrapper/bridge layer for the OxiPhysics engine.
//!
//! This crate provides Python-friendly types and a high-level API surface
//! designed for future PyO3/pyo3 FFI integration. All types use simple
//! owned data (no lifetimes) and derive Serialize/Deserialize for easy
//! JSON interchange.
//!
//! ## Quick start
//!
//! ```no_run
//! use oxiphysics_python::world_api::PyPhysicsWorld;
//! use oxiphysics_python::types::{PySimConfig, PyRigidBodyConfig};
//!
//! let mut world = PyPhysicsWorld::new(PySimConfig::earth_gravity());
//! let h = world.add_rigid_body(PyRigidBodyConfig::dynamic(1.0, [0.0, 10.0, 0.0]));
//! world.step(1.0 / 60.0);
//! let pos = world.get_position(h).expect("body exists");
//! assert!(pos[1] < 10.0); // body has fallen under gravity
//! ```

#![allow(missing_docs)]

mod error;
pub mod fem_api;
pub mod lbm_api;
pub mod md_api;
pub mod serialization;
pub mod sph_api;
pub mod types;
pub mod world_api;

pub use error::*;
pub use fem_api::{
    PyFemDirichletBC, PyFemElement, PyFemMaterial, PyFemMesh, PyFemNodalForce, PyFemNode,
    PyFemSolveResult, PyFemSolver,
};
pub use lbm_api::{LbmBoundary, PyLbmConfig, PyLbmSimulation};
pub use md_api::{PyMdAtom, PyMdConfig, PyMdSimulation};
pub use sph_api::{PySphConfig, PySphSimulation};
pub use types::*;
pub use world_api::PyPhysicsWorld;
pub mod analytics_api;
pub mod constraints_api;
pub mod geometry_api;
pub mod io_api;
pub mod materials_api;
pub mod rigid_api;
pub mod vehicle_api;
pub mod viz_api;
