// Copyright (C) 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! wasm32 regression tests for the JavaScript-facing error paths and the
//! profiler clock.
//!
//! Building a `JsValue` (every thrown error is one) only works inside a real
//! WebAssembly host, so these checks cannot run natively; the Rust-level
//! logic behind them is covered by the unit tests of each module.  They run
//! under Node.js via `wasm-bindgen-test-runner` (no `run_in_browser`):
//!
//! ```sh
//! CARGO_TARGET_WASM32_UNKNOWN_UNKNOWN_RUNNER=wasm-bindgen-test-runner \
//!     cargo test --target wasm32-unknown-unknown -p oxiphysics-wasm \
//!     --test wasm_bridge_error_paths
//! ```

#![cfg(target_arch = "wasm32")]

use oxiphysics_wasm::WasmNavMeshJs;
use oxiphysics_wasm::fluid_bridge::{WasmLbmConfig, WasmLbmSimulation};
use oxiphysics_wasm::io_bridge::MeshExporter;
use oxiphysics_wasm::profiler_bridge::WasmProfilerSession;
use oxiphysics_wasm::wasm_helpers::now_ms;
use oxiphysics_wasm::web_worker::{SharedStateBuffer, WorkerBridge, WorkerRuntime};
use wasm_bindgen_test::*;

/// Busy-wait for `ms` milliseconds on the same clock the profiler uses.
fn spin_ms(ms: f64) {
    let start = now_ms().expect("performance.now() is available under Node.js");
    while now_ms().expect("performance.now() is available") - start < ms {}
}

// ---------------------------------------------------------------------------
// Profiler: `std::time::Instant::now` panics on wasm32-unknown-unknown
// ---------------------------------------------------------------------------

#[wasm_bindgen_test]
fn profiler_frame_does_not_trap_and_measures_time() {
    let mut profiler = WasmProfilerSession::new_js();
    profiler.begin_frame_js();
    profiler.push_scope_js("physics");
    spin_ms(5.0);
    profiler.push_scope_js("narrowphase");
    spin_ms(1.0);
    profiler.pop_scope_js();
    profiler.pop_scope_js();
    assert_eq!(profiler.scope_depth_js(), 0);
    let json = profiler.end_frame_json_js();
    let report: serde_json::Value = serde_json::from_str(&json).expect("report is valid JSON");
    assert_eq!(report["clock_available"], serde_json::Value::Bool(true));
    let physics = &report["root"]["children"][0];
    assert_eq!(physics["name"], "physics");
    let total_ns = physics["total_ns"].as_u64().expect("total_ns is a number");
    assert!(total_ns >= 5_000_000, "5 ms scope reported {total_ns} ns");
    assert!(
        report["frame_ns"].as_u64().expect("frame_ns") >= total_ns,
        "the frame cannot be shorter than its scope"
    );
}

#[wasm_bindgen_test]
fn profiler_report_object_is_available_after_a_frame() {
    let mut profiler = WasmProfilerSession::new_js();
    profiler.begin_frame_js();
    profiler.push_scope_js("step");
    profiler.pop_scope_js();
    assert!(profiler.end_frame_js().is_ok());
    // A second frame works as well, and leftover scopes are discarded.
    profiler.begin_frame_js();
    profiler.push_scope_js("leaked");
    assert!(profiler.end_frame_js().is_ok());
    assert_eq!(profiler.scope_depth_js(), 0);
}

// ---------------------------------------------------------------------------
// SharedStateBuffer / WorkerRuntime capacity
// ---------------------------------------------------------------------------

#[wasm_bindgen_test]
fn shared_buffer_js_write_past_capacity_throws_and_changes_nothing() {
    let mut buffer = SharedStateBuffer::new(2);
    buffer
        .write_body_position_js(1, vec![1.0, 2.0, 3.0])
        .expect("a valid write succeeds");
    let before = buffer.as_bytes_js();

    let err = buffer
        .write_body_position_js(2, vec![9.0, 9.0, 9.0])
        .expect_err("index 2 of 2 bodies must throw");
    let message = err.as_string().expect("errors are thrown as strings");
    assert!(message.contains("body_index"), "{message}");
    assert!(
        buffer
            .write_body_position_js(u32::MAX as usize, vec![0.0; 3])
            .is_err()
    );
    assert!(buffer.write_body_position_js(0, vec![1.0, 2.0]).is_err());
    assert_eq!(
        buffer.as_bytes_js(),
        before,
        "a rejected write modified the buffer"
    );
    assert_eq!(buffer.read_body_position_js(1), vec![1.0, 2.0, 3.0]);
    assert_eq!(buffer.read_body_position_js(2), vec![0.0; 3]);
}

#[wasm_bindgen_test]
fn worker_runtime_reports_capacity_errors_instead_of_trapping() {
    let mut bridge = WorkerBridge::new(2);
    let mut runtime = WorkerRuntime::new();
    for i in 0..3 {
        bridge.post_command_add_sphere(1.0, i as f32, 10.0, 0.0, 0.5);
    }
    bridge.post_command_step(1.0 / 60.0);
    runtime.process(&mut bridge);
    // Flat result codes: 2 = BodyAdded, 6 = Error, 1 = StepDone.
    assert_eq!(bridge.poll_result_flat()[0], 2.0);
    assert_eq!(bridge.poll_result_flat()[0], 2.0);
    assert_eq!(
        bridge.poll_result_flat()[0],
        6.0,
        "third sphere must be refused"
    );
    assert_eq!(bridge.poll_result_flat()[0], 1.0, "the step must not trap");
    assert_eq!(bridge.shared_body_count(), 2);
}

// ---------------------------------------------------------------------------
// MeshExporter
// ---------------------------------------------------------------------------

#[wasm_bindgen_test]
fn mesh_exporter_add_face_js_throws_for_missing_vertices() {
    let mut mesh = MeshExporter::new();
    mesh.add_vertex_js(0.0, 0.0, 0.0);
    mesh.add_vertex_js(1.0, 0.0, 0.0);
    mesh.add_vertex_js(0.0, 1.0, 0.0);
    let err = mesh
        .add_face_js(0, 1, 99)
        .expect_err("vertex 99 does not exist");
    let message = err.as_string().expect("errors are thrown as strings");
    assert!(message.contains("99"), "{message}");
    assert_eq!(mesh.face_count_js(), 0);
    assert!(!mesh.export_stl_ascii().contains("facet"));
    mesh.add_face_js(0, 1, 2).expect("a valid face is accepted");
    assert_eq!(mesh.face_count_js(), 1);
    assert!(mesh.export_stl_ascii().contains("facet normal 0 0 1"));
}

// ---------------------------------------------------------------------------
// WasmNavMeshJs
// ---------------------------------------------------------------------------

#[wasm_bindgen_test]
fn navmesh_from_flat_throws_instead_of_trapping_on_bad_indices() {
    let vertices = vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0];
    // Index 4 does not exist; the old constructor accepted it and the first
    // path query aborted the whole instance.
    assert!(WasmNavMeshJs::from_flat(vertices.clone(), vec![0, 1, 4]).is_err());
    assert!(WasmNavMeshJs::from_flat(vertices.clone(), vec![0, 1, 1]).is_err());
    assert!(WasmNavMeshJs::from_flat(vertices.clone(), vec![0, 1]).is_err());
    assert!(WasmNavMeshJs::from_flat(vec![0.0; 11], vec![0, 1, 2]).is_err());
    let nav = WasmNavMeshJs::from_flat(vertices, vec![0, 1, 2, 1, 3, 2])
        .expect("a valid mesh is accepted");
    let path = nav.find_path_flat(0.1, 0.0, 0.1, 0.9, 0.0, 0.9, 0.0);
    assert!(path.len() >= 6);
}

// ---------------------------------------------------------------------------
// LBM inlet velocity
// ---------------------------------------------------------------------------

#[wasm_bindgen_test]
fn lbm_set_inlet_velocity_js_rejects_unusable_speeds() {
    let mut sim = WasmLbmSimulation::create(WasmLbmConfig::d3q19_js(6, 4, 4, 0.1));
    assert!(sim.set_inlet_velocity_js(f64::NAN).is_err());
    assert!(sim.set_inlet_velocity_js(1.5).is_err());
    sim.set_inlet_velocity_js(0.05)
        .expect("a small speed is valid");
    for _ in 0..5 {
        sim.step_js();
    }
    let ux = sim.get_velocity_field_js()[0];
    assert!((ux - 0.05).abs() < 1e-9, "inlet node ux = {ux}");
    sim.clear_inlet_velocity_js();
}
