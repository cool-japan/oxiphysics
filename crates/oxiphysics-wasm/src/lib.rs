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

use wasm_bindgen::prelude::*;

#[wasm_bindgen(start)]
pub fn wasm_main() {
    console_error_panic_hook::set_once();
}

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

pub mod web_worker;
pub use web_worker::{SharedStateBuffer, SimCommand, SimResult, WorkerBridge, WorkerRuntime};

pub mod wasm_bench;
pub mod wasm_helpers;
pub mod webgl_bridge;

// Phase 20.4: WASM bindings for Phase 13-19 modules
pub mod animation_bridge;
pub mod character_bridge;
pub mod navmesh_bridge;
pub mod profiler_bridge;
pub mod rollback_bridge;
pub mod rope_bridge;
pub mod telemetry_bridge;
pub mod triggers_bridge;

pub use animation_bridge::WasmAnimationPlayer;
pub use character_bridge::WasmCharacterController;
pub use navmesh_bridge::WasmNavMesh;
pub use profiler_bridge::WasmProfilerSession;
pub use rollback_bridge::WasmRollbackSession;
pub use rope_bridge::WasmRope;
pub use telemetry_bridge::WasmTelemetrySession;
pub use triggers_bridge::WasmTriggerWorld;

// WASM-bindgen JS wrappers for bridge types (Phase 7.5)
pub mod character_wasm;
pub mod navmesh_wasm;
pub mod rope_wasm;
pub use character_wasm::WasmCharacterControllerJs;
pub use navmesh_wasm::WasmNavMeshJs;
pub use rope_wasm::WasmRopeJs;
