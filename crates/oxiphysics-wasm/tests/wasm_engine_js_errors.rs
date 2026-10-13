// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! wasm32 regression tests for the JavaScript-facing contract of
//! `WasmPhysicsEngine`: invalid input is rejected with a thrown exception
//! (a `JsValue` carrying the error message) instead of being stored silently,
//! and handles of removed bodies stop working.
//!
//! Building a `JsValue` only works inside a real WebAssembly host, so these
//! checks cannot run natively; the Rust-level logic behind them is covered by
//! the unit tests of the `engine` module.  They run under Node.js via
//! `wasm-bindgen-test-runner` (no `run_in_browser`):
//!
//! ```sh
//! CARGO_TARGET_WASM32_UNKNOWN_UNKNOWN_RUNNER=wasm-bindgen-test-runner \
//!     cargo test --target wasm32-unknown-unknown -p oxiphysics-wasm \
//!     --test wasm_engine_js_errors
//! ```

#![cfg(target_arch = "wasm32")]

use oxiphysics_wasm::WasmPhysicsEngine;
use oxiphysics_wasm::engine::WasmEngine;
use oxiphysics_wasm::js_api::{JsBodyDesc, JsPhysicsConfig, engine_add_bodies_from_json};
use wasm_bindgen::JsValue;
use wasm_bindgen_test::*;

/// The message of a thrown JavaScript string.
fn message(error: JsValue) -> String {
    error.as_string().unwrap_or_default()
}

fn engine() -> WasmPhysicsEngine {
    WasmPhysicsEngine::new(0.0, -9.81, 0.0)
}

#[wasm_bindgen_test]
fn dynamic_bodies_need_a_finite_positive_mass() {
    let mut engine = engine();
    for mass in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let error = engine
            .add_dynamic_body_js(mass, 0.0, 0.0, 0.0)
            .expect_err("an unusable mass must throw");
        let text = message(error);
        assert!(
            text.contains("mass"),
            "message for mass {mass} should name the mass, got {text:?}"
        );
    }
    assert_eq!(engine.get_body_count_js(), 0, "nothing may have been added");
    let handle = engine
        .add_dynamic_body_js(1.0, 0.0, 0.0, 0.0)
        .expect("a valid mass is accepted");
    assert_eq!(engine.get_body_count_js(), 1);
    assert!(engine.contains_body_js(handle));
}

#[wasm_bindgen_test]
fn positions_must_be_finite() {
    let mut engine = engine();
    let error = engine
        .add_static_body_js(f64::NAN, 0.0, 0.0)
        .expect_err("a NaN position must throw");
    assert!(message(error).contains("position"));
    assert_eq!(engine.get_body_count_js(), 0);
}

#[wasm_bindgen_test]
fn colliders_need_a_live_body_and_valid_dimensions() {
    let mut engine = engine();
    let body = engine
        .add_dynamic_body_js(1.0, 0.0, 0.0, 0.0)
        .expect("body");

    let error = engine
        .add_sphere_collider_js(body, -0.5)
        .expect_err("a negative radius must throw");
    assert!(message(error).contains("dimension"));
    let error = engine
        .add_box_collider_js(body, 1.0, 0.0, 1.0)
        .expect_err("a flat box must throw");
    assert!(message(error).contains("dimension"));
    let error = engine
        .add_capsule_collider_js(body, 0.5, f64::NAN)
        .expect_err("a NaN capsule height must throw");
    assert!(message(error).contains("dimension"));
    let error = engine
        .add_plane_collider_js(body, 0.0, 0.0, 0.0, 1.0)
        .expect_err("a zero plane normal must throw");
    assert!(message(error).contains("plane_normal"));

    // Handles of removed bodies are rejected, also after their slot is reused.
    assert!(engine.remove_body_js(body));
    let error = engine
        .add_sphere_collider_js(body, 0.5)
        .expect_err("a removed body cannot get colliders");
    assert!(message(error).contains("body not found"));
    let reused = engine
        .add_dynamic_body_js(1.0, 5.0, 0.0, 0.0)
        .expect("reused slot");
    assert_ne!(reused, body, "a recycled slot hands out a new handle");
    assert!(
        engine.add_sphere_collider_js(body, 0.5).is_err(),
        "the stale handle must not reach the new occupant of the slot"
    );
    assert_eq!(engine.get_active_collider_count(), 0);
    assert_eq!(engine.get_position_js(reused), vec![5.0, 0.0, 0.0]);
}

#[wasm_bindgen_test]
fn stale_handles_are_rejected_by_the_boolean_setters() {
    let mut engine = engine();
    let old = engine
        .add_dynamic_body_js(1.0, 0.0, 0.0, 0.0)
        .expect("old body");
    assert!(engine.remove_body_js(old));
    let new = engine
        .add_dynamic_body_js(1.0, 1.0, 2.0, 3.0)
        .expect("new body");
    assert!(
        !engine.remove_body_js(old),
        "the stale handle cannot remove"
    );
    assert!(!engine.set_velocity_js(old, 9.0, 0.0, 0.0));
    assert!(!engine.apply_torque_js(old, 0.0, 1.0, 0.0));
    assert_eq!(engine.get_velocity_js(new), vec![0.0, 0.0, 0.0]);
    assert_eq!(engine.get_position_js(new), vec![1.0, 2.0, 3.0]);
    assert!(engine.contains_body_js(new));
}

#[wasm_bindgen_test]
fn step_rejects_invalid_time_steps() {
    let mut engine = engine();
    engine
        .add_dynamic_body_js(1.0, 0.0, 10.0, 0.0)
        .expect("body");
    for dt in [-0.1, f64::NAN, f64::INFINITY] {
        let error = engine
            .step_js(dt)
            .expect_err("an invalid time step must throw");
        assert!(message(error).contains("time step"));
    }
    assert_eq!(engine.time_js(), 0.0, "a rejected step changes nothing");
    engine.step_js(1.0 / 60.0).expect("a valid step");
    assert!(engine.time_js() > 0.0);
}

#[wasm_bindgen_test]
fn fixed_substep_must_be_positive() {
    let mut engine = engine();
    for dt in [0.0, -0.01, f64::NAN] {
        assert!(engine.set_fixed_dt_js(dt).is_err(), "dt = {dt}");
    }
    engine.set_fixed_dt_js(0.01).expect("a valid substep");
    engine.step_js(0.05).expect("step");
    assert_eq!(engine.substeps_performed_js(), 4, "capped at max_substeps");
    assert!(engine.dropped_time_js() > 0.0);
}

#[wasm_bindgen_test]
fn unusable_engine_configurations_throw() {
    let error = JsPhysicsConfig::new()
        .with_fixed_dt(0.0)
        .build_js()
        .expect_err("a zero substep must throw");
    assert!(message(error).contains("time step"));
    assert!(
        JsPhysicsConfig::new()
            .with_fixed_dt(1.0 / 120.0)
            .build_js()
            .is_ok()
    );
}

#[wasm_bindgen_test]
fn a_ball_falls_onto_the_ground_and_rests_on_it() {
    let mut engine = engine();
    let ball = engine
        .add_dynamic_body_js(1.0, 0.0, 3.0, 0.0)
        .expect("ball");
    engine
        .add_sphere_collider_js(ball, 0.5)
        .expect("ball shape");
    let ground = engine.add_static_body_js(0.0, 0.0, 0.0).expect("ground");
    engine
        .add_plane_collider_js(ground, 0.0, 2.0, 0.0, 0.0)
        .expect("ground plane (the normal need not be a unit vector)");
    for _ in 0..240 {
        engine.step_js(1.0 / 60.0).expect("step");
    }
    let height = engine.get_position_js(ball)[1];
    assert!(
        (height - 0.5).abs() < 0.02,
        "the ball must rest on the plane, y = {height}"
    );
    assert!((engine.time_js() - 4.0).abs() < 1e-9);
}

#[wasm_bindgen_test]
fn rays_hit_boxes_and_report_the_distance() {
    let mut engine = engine();
    let crate_body = engine
        .add_static_body_js(5.0, 0.0, 0.0)
        .expect("crate body");
    engine
        .add_box_collider_js(crate_body, 0.5, 0.5, 0.5)
        .expect("crate shape");
    let hit = engine.raycast_js(0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 100.0);
    assert_eq!(hit.len(), 2, "[handle, distance] for a hit");
    assert_eq!(hit[0], f64::from(crate_body));
    assert!((hit[1] - 4.5).abs() < 1e-9, "distance {}", hit[1]);
    assert!(
        engine
            .raycast_js(0.0, 5.0, 0.0, 1.0, 0.0, 0.0, 100.0)
            .is_empty(),
        "a ray passing above the box misses"
    );
}

#[wasm_bindgen_test]
fn composite_helpers_create_nothing_when_they_fail() {
    let mut helper = WasmEngine::new(0.0, -9.81, 0.0);
    assert!(
        helper.add_dynamic_sphere(1.0, 0.0, 0.0, 0.0, -0.5).is_err(),
        "a negative radius must throw"
    );
    assert!(
        helper
            .add_static_plane(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0)
            .is_err(),
        "a zero plane normal must throw"
    );
    assert_eq!(
        helper.body_count(),
        0,
        "a failed helper leaves no body behind"
    );
}

#[wasm_bindgen_test]
fn batch_creation_from_json_is_all_or_nothing() {
    let mut engine = engine();
    let descs = vec![
        JsBodyDesc::sphere(1.0, 0.0, 1.0, 0.0, 0.5),
        JsBodyDesc::sphere(1.0, 3.0, 1.0, 0.0, -0.5),
    ];
    let json = serde_json::to_string(&descs).expect("serialize descriptors");
    let error =
        engine_add_bodies_from_json(&mut engine, &json).expect_err("a negative radius must throw");
    assert!(message(error).contains("dimension"));
    assert_eq!(
        engine.get_body_count_js(),
        0,
        "the valid descriptor of the failed batch is rolled back"
    );
    assert!(engine_add_bodies_from_json(&mut engine, "not json").is_err());
    assert_eq!(engine.get_body_count_js(), 0);
}
