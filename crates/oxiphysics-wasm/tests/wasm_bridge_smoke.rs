// Copyright (C) 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Smoke tests for the WASM bridge types.
//!
//! These tests run inside a real WASM environment (via `wasm-pack test
//! --headless`). Each test constructs the bridge type and exercises one method
//! to prove the binding compiles and the runtime path is reachable.

#![cfg(target_arch = "wasm32")]

use wasm_bindgen_test::*;

wasm_bindgen_test_configure!(run_in_browser);

// ---------------------------------------------------------------------------
// WasmPhysicsEngine
// ---------------------------------------------------------------------------

#[wasm_bindgen_test]
fn smoke_physics_engine() {
    let mut engine = oxiphysics_wasm::WasmPhysicsEngine::new(0.0, -9.81, 0.0);
    // Add a body and step — proves the full simulation path is reachable.
    let body = engine.add_dynamic_body(1.0, 0.0, 10.0, 0.0);
    engine.step(1.0 / 60.0);
    let pos = engine.get_position(body);
    // After one step under gravity the body should have moved downward.
    assert!(pos[1] < 10.0, "body should have fallen under gravity");
}

// ---------------------------------------------------------------------------
// WasmRopeJs
// ---------------------------------------------------------------------------

#[wasm_bindgen_test]
fn smoke_rope() {
    let mut rope = oxiphysics_wasm::WasmRopeJs::new(0.0, 5.0, 0.0, 8, 0.5, 1.0)
        .expect("rope construction should succeed");
    rope.pin_head(0.0, 5.0, 0.0);
    rope.step(1.0 / 60.0);
    let pos = rope.get_position(0);
    assert_eq!(pos.len(), 3, "position should be [x, y, z]");
    // Head is pinned — should stay at (0, 5, 0).
    assert!((pos[1] - 5.0).abs() < 0.01, "pinned head should not drift");
    assert_eq!(rope.link_count(), 8);
}

// ---------------------------------------------------------------------------
// WasmNavMeshJs
// ---------------------------------------------------------------------------

#[wasm_bindgen_test]
fn smoke_navmesh() {
    let nav = oxiphysics_wasm::WasmNavMeshJs::grid(10, 10, 1.0);
    assert!(
        nav.triangle_count() > 0,
        "grid navmesh should have triangles"
    );

    // Pathfind across the mesh.
    let path = nav.find_path_flat(0.5, 0.0, 0.5, 9.5, 0.0, 9.5, 0.0);
    assert!(
        path.len() >= 6,
        "path should contain at least two waypoints"
    );
}

// ---------------------------------------------------------------------------
// WasmCharacterControllerJs
// ---------------------------------------------------------------------------

#[wasm_bindgen_test]
fn smoke_character_controller() {
    let mut ctrl = oxiphysics_wasm::WasmCharacterControllerJs::new(0.0, 2.0, 0.0, 0.4, 0.9);

    let initial_pos = ctrl.get_position();
    assert_eq!(initial_pos.len(), 3);
    assert!(
        (initial_pos[1] - 2.0).abs() < 1e-9,
        "initial y should be 2.0"
    );

    // Apply gravity for one frame and slide (no hits).
    ctrl.apply_gravity(9.81, 1.0 / 60.0);
    let result = ctrl.move_and_slide(0.0, ctrl.get_velocity()[1] / 60.0, 0.0, vec![]);
    assert_eq!(
        result.len(),
        4,
        "result should be [tx, ty, tz, is_grounded]"
    );

    // After falling the controller should have moved downward.
    let pos = ctrl.get_position();
    assert!(
        pos[1] <= 2.0,
        "character should not move upward under gravity"
    );
}
