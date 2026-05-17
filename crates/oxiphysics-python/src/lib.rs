// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Python bindings for the OxiPhysics engine.
//!
//! This crate is built as a `cdylib` Python extension module via
//! [maturin](https://www.maturin.rs/) / [PyO3](https://pyo3.rs/).
//!
//! ## Usage
//!
//! ```python
//! import oxiphysics
//!
//! print(oxiphysics.__version__)
//!
//! world = oxiphysics.PhysicsWorld(gy=-9.81)
//!
//! cfg = oxiphysics.RigidBodyConfig.dynamic(mass=1.0, y=10.0)
//! cfg.add_sphere(radius=0.5)
//! ball = world.add_rigid_body(cfg)
//!
//! for _ in range(60):
//!     world.step(1.0 / 60.0)
//!
//! x, y, z = world.get_position(ball)
//! print(f"Ball fell to y={y:.3f}")
//! ```

#![allow(missing_docs)]

use pyo3::prelude::*;

pub mod aero_api;
pub mod analytics_api;
pub mod animation_api;
pub mod buoyancy_api;
pub mod character_api;
pub mod constraints_api;
pub mod contact_cache_api;
pub mod debug_draw_api;
mod error;
pub mod event_bus_api;
pub mod fem_api;
pub mod force_field_api;
pub mod geometry_api;
pub mod ik_api;
pub mod interpolator_api;
pub mod io_api;
pub mod lbm_api;
pub mod lod_api;
pub mod material_table_api;
pub mod materials_api;
pub mod md_api;
pub mod navmesh_api;
pub mod noise_api;
pub mod profiler_api;
pub mod py_classes;
pub mod query_api;
pub mod rasterizer;
pub mod replay_api;
pub mod rigid_api;
pub mod rollback_api;
pub mod rope_api;
pub mod scene_api;
pub mod scheduler_api;
pub mod serialization;
pub mod snapshot_api;
pub mod spatial_grid_api;
pub mod sph_api;
pub mod telemetry_api;
pub mod trigger_api;
pub mod types;
pub mod vehicle_api;
#[cfg(test)]
mod vehicle_api_tests;
pub mod viz_api;
pub mod world_api;
pub mod xpbd_api;

pub use error::*;
pub use fem_api::{
    PyFemDirichletBC, PyFemElement, PyFemMaterial, PyFemMesh, PyFemNodalForce, PyFemNode,
    PyFemSolveResult, PyFemSolver,
};
pub use lbm_api::{LbmBoundary, PyLbmConfig, PyLbmSimulation};
pub use md_api::{PyMdAtom, PyMdConfig, PyMdSimulation};
pub use py_classes::{ContactResult, PhysicsWorld, RigidBodyConfig, SimConfig};
pub use sph_api::{PySphConfig, PySphSimulation};
pub use types::*;
pub use world_api::PyPhysicsWorld;

/// OxiPhysics Python extension module.
///
/// Exposes the following classes:
/// - `SimConfig`        — simulation configuration (gravity, solver settings)
/// - `RigidBodyConfig`  — body creation parameters (mass, position, colliders)
/// - `ContactResult`    — contact manifold from the most recent step
/// - `PhysicsWorld`     — the primary simulation world
#[pymodule]
fn oxiphysics(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;

    // Top-level legacy classes
    m.add_class::<py_classes::SimConfig>()?;
    m.add_class::<py_classes::RigidBodyConfig>()?;
    m.add_class::<py_classes::ContactResult>()?;
    m.add_class::<py_classes::PhysicsWorld>()?;

    // Phase 20.x bindings
    noise_api::register(m)?;

    // Domain sub-modules (Wave 2)
    analytics_api::register_analytics_module(m)?;
    constraints_api::register_constraints_module(m)?;
    fem_api::register_fem_module(m)?;
    geometry_api::register_geometry_module(m)?;
    io_api::register_io_module(m)?;
    lbm_api::register_lbm_module(m)?;
    materials_api::register_materials_module(m)?;
    md_api::register_md_module(m)?;
    rigid_api::register_rigid_module(m)?;
    sph_api::register_sph_module(m)?;
    vehicle_api::register_vehicle_module(m)?;
    viz_api::register_viz_module(m)?;
    world_api::register_world_module(m)?;

    // Phase 6 orphan modules (wired 2026-05-11)
    aero_api::register(m)?;
    animation_api::register(m)?;
    buoyancy_api::register(m)?;
    character_api::register(m)?;
    contact_cache_api::register(m)?;
    debug_draw_api::register(m)?;
    event_bus_api::register(m)?;
    force_field_api::register(m)?;
    ik_api::register(m)?;
    interpolator_api::register(m)?;
    lod_api::register(m)?;
    material_table_api::register(m)?;
    navmesh_api::register(m)?;
    profiler_api::register(m)?;
    query_api::register(m)?;
    replay_api::register(m)?;
    rollback_api::register(m)?;
    rope_api::register(m)?;
    scene_api::register(m)?;
    scheduler_api::register(m)?;
    snapshot_api::register(m)?;
    spatial_grid_api::register(m)?;
    telemetry_api::register(m)?;
    trigger_api::register(m)?;
    xpbd_api::register(m)?;

    Ok(())
}
