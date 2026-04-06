// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! WebAssembly-friendly wrapper/bridge layer for the OxiPhysics engine.
//!
//! This crate provides WASM-friendly types and a flat API surface
//! designed for future wasm-bindgen integration. All methods take and
//! return primitives or flat arrays for easy JavaScript interop.
//!
//! ## Quick Start
//!
//! ```no_run
//! use oxiphysics_wasm::{WasmPhysicsEngine, SimulationConfig};
//!
//! // Create an engine with Earth gravity
//! let mut engine = WasmPhysicsEngine::new(0.0, -9.81, 0.0);
//!
//! // Add a 1 kg sphere at height 10 m
//! let ball = engine.add_dynamic_body(1.0, 0.0, 10.0, 0.0);
//! engine.add_sphere_collider(ball, 0.5);
//!
//! // Add a ground plane
//! let ground = engine.add_static_body(0.0, 0.0, 0.0);
//! engine.add_plane_collider(ground, 0.0, 1.0, 0.0, 0.0);
//!
//! // Simulate for 1 second (many substeps internally)
//! engine.step(1.0);
//!
//! // Query position
//! let pos = engine.get_position(ball);
//! assert!(pos[1] < 10.0, "ball should have fallen");
//! ```
#![allow(missing_docs)]
#![allow(dead_code)]

pub mod body_query;
pub mod engine;
pub mod error;
pub mod events;
pub mod js_api;
pub mod math_helpers;
pub mod physics_config;
pub mod renderer;
pub mod types;

pub use engine::WasmPhysicsEngine;
pub use error::{Error, Result};
pub use types::{
    BodyState, ColliderConfig, ColliderShapeType, ContactResult, DebugInfo, QuatWasm,
    RaycastResult, RigidBodyConfig, SimulationConfig, TransformWasm, Vec3Wasm,
};
pub mod analytics_bridge;
pub mod constraint_bridge;
pub mod debug_tools;
pub mod fluid_bridge;
pub mod io_bridge;
pub mod material_bridge;
pub mod particle_system;
pub mod sim_controls;
pub mod simulation_api;
