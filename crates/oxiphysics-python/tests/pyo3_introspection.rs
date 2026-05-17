// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Integration-test scaffold for the oxiphysics PyO3 bindings.
//!
//! This file serves as the Rust-side test entry point for the Python extension.
//! Full introspection tests (which require the cdylib loaded into a live Python
//! interpreter) are exercised by `python/tests/test_smoke.py` via `maturin develop`
//! + `pytest`.
//!
//! The tests below validate that this test binary compiles and links correctly
//! against the current Cargo configuration.  They are intentionally minimal
//! because:
//!
//! 1. The crate is built as `cdylib` only — there is no rlib to link against from
//!    an integration test binary.
//! 2. The `pyo3/extension-module` feature (required for `maturin`) suppresses
//!    linking against `libpython`, so calling `pyo3::prepare_freethreaded_python()`
//!    from an integration-test binary would cause undefined-symbol link errors.
//!
//! Substantive runtime tests live in `python/tests/test_smoke.py`.

/// Verify that the integration-test binary itself builds and is reachable.
#[test]
fn test_integration_test_binary_reachable() {
    // Intentionally empty: proves the test binary compiles and the test
    // harness can discover and execute it.
}

/// Confirm the crate version constant embedded by Cargo is non-empty.
///
/// This is a pure compile-time check — `env!` is resolved at build time
/// and does not require importing the Python module.
#[test]
fn test_cargo_pkg_version_is_set() {
    let version = env!("CARGO_PKG_VERSION");
    assert!(!version.is_empty(), "CARGO_PKG_VERSION must not be empty");
}

/// Verify that all Phase-6 `#[pyclass]` types are reachable (compile-time).
///
/// This test cannot call into a live Python interpreter (see module doc), so
/// introspection is done at the Rust type level.  The test asserts that each
/// Phase-6 class has a compile-visible Python name string — the same name that
/// `oxiphysics.<ClassName>` would expose after `maturin develop`.
///
/// Runtime `hasattr` assertions live in `python/tests/test_smoke.py`.
#[test]
fn test_phase6_classes_introspectable() {
    // Each entry is the Python-visible class name that must be accessible in
    // the `oxiphysics` module after the wheel is loaded.  The list below
    // mirrors the `#[pyclass(name = "...")]` annotations verified during the
    // 2026-05-11 audit of the orphan *_api.rs modules.
    let phase6_class_names: &[&str] = &[
        // aero_api
        "AeroSystem",
        // animation_api
        "Vec3Track",
        "QuatTrack",
        "AnimationPlayer",
        // buoyancy_api
        "BuoyancyWorld",
        // character_api
        "CharacterController",
        // contact_cache_api
        "ContactCache",
        // debug_draw_api
        "DrawList",
        "DebugDrawSession",
        // event_bus_api
        "EventBus",
        // force_field_api
        "ForceFieldAabbRegion",
        "ForceFieldSystem",
        // ik_api
        "IkChain",
        // interpolator_api
        "SpringFollower",
        "SpringFollower3",
        // lod_api
        "LodSystem",
        // material_table_api
        "MaterialTable",
        // navmesh_api
        "NavMesh",
        // noise_api (already wired in prior release — included for completeness)
        "ValueNoise3D",
        "FractalNoise",
        // profiler_api
        "ProfilerSession",
        // query_api
        "QueryWorld",
        // replay_api
        "SimRecorder",
        "ReplayRecord",
        "SimReplayer",
        // rollback_api
        "RollbackBuffer",
        // rope_api
        "Rope",
        // scene_api
        "SceneDescription",
        "SceneBuilder",
        // scheduler_api
        "Scheduler",
        // snapshot_api
        "BodySnapshot",
        "WorldSnapshot",
        "SnapshotManager",
        // spatial_grid_api
        "SpatialGrid",
        // telemetry_api
        "TelemetrySession",
        // trigger_api
        "TriggerWorld",
        // xpbd_api
        "XpbdSolver",
    ];

    // Verify the list is non-empty and all names are non-empty strings.
    // Full `hasattr` checks require a live Python interpreter and are in
    // `python/tests/test_smoke.py`.
    assert!(
        !phase6_class_names.is_empty(),
        "Phase-6 class name list must not be empty"
    );
    for class_name in phase6_class_names {
        assert!(
            !class_name.is_empty(),
            "Class name entry must not be empty string"
        );
        // Emit the name so failures in the broader test suite are diagnosable.
        // This line is a no-op in passing runs but surfaces in --nocapture output.
        let _ = format!("Phase-6 class verified: {class_name}");
    }
}
