// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Static checks of the PyO3 wiring of the oxiphysics Python module.
//!
//! The crate is built as `cdylib` only: the library is named `oxiphysics` (the
//! name of the Python module), which is also the name of the umbrella
//! dependency, so an `rlib` target would make `use oxiphysics::…` in an
//! integration test ambiguous. An integration test therefore cannot link the
//! library, and importing the module needs a live interpreter, which the
//! pytest suite in `python/tests` provides after `maturin develop`.
//!
//! What needs neither is the *wiring*, and the tests below check it by reading
//! the crate's own sources: every Python class the bindings promise is declared
//! with `#[pyclass(name = "...")]` and handed to `add_class`, the `pyo3`
//! dependency does not enable `extension-module` (maturin does, see
//! `pyproject.toml`), and the Python package version follows the crate version.
//!
//! The unit tests inside `src/` need no `rlib`: they are compiled into the
//! library's own test target and run with `cargo nextest run -p
//! oxiphysics-python`. Because Cargo does not enable `pyo3/extension-module`,
//! those test binaries link `libpython` normally.

use std::fs;
use std::path::{Path, PathBuf};

/// Root of the crate, resolved at compile time.
fn crate_root() -> &'static Path {
    Path::new(env!("CARGO_MANIFEST_DIR"))
}

/// Text of a file below the crate root.
fn read_crate_file(relative: &str) -> String {
    let path = crate_root().join(relative);
    fs::read_to_string(&path).unwrap_or_else(|e| panic!("cannot read {}: {e}", path.display()))
}

/// Text of a file below the workspace root (two levels above the crate root).
fn read_workspace_file(relative: &str) -> String {
    let path = crate_root().join("..").join("..").join(relative);
    fs::read_to_string(&path).unwrap_or_else(|e| panic!("cannot read {}: {e}", path.display()))
}

/// Every `.rs` file below `dir`, depth first, in a stable order.
fn collect_rust_sources(dir: &Path, out: &mut Vec<PathBuf>) {
    let entries =
        fs::read_dir(dir).unwrap_or_else(|e| panic!("cannot list {}: {e}", dir.display()));
    let mut paths: Vec<PathBuf> = entries
        .filter_map(|entry| entry.ok())
        .map(|e| e.path())
        .collect();
    paths.sort();
    for path in paths {
        if path.is_dir() {
            collect_rust_sources(&path, out);
        } else if path.extension().is_some_and(|ext| ext == "rs") {
            out.push(path);
        }
    }
}

/// `(relative path, text)` of every Rust source file of the crate.
fn crate_sources() -> Vec<(String, String)> {
    let src = crate_root().join("src");
    let mut files = Vec::new();
    collect_rust_sources(&src, &mut files);
    files
        .into_iter()
        .map(|path| {
            let text = fs::read_to_string(&path)
                .unwrap_or_else(|e| panic!("cannot read {}: {e}", path.display()));
            let relative = path
                .strip_prefix(crate_root())
                .unwrap_or(&path)
                .display()
                .to_string();
            (relative, text)
        })
        .collect()
}

/// A `#[pyclass(name = "...")]` declaration found in the sources.
#[derive(Debug)]
struct PyClassDecl {
    /// Class name visible in Python.
    python_name: String,
    /// Rust type the attribute is attached to.
    rust_type: String,
    /// File the declaration is in.
    file: String,
}

/// Value of the `name = "..."` argument of a `#[pyclass(...)]` argument list.
fn pyclass_python_name(args: &str) -> Option<String> {
    args.split(',').find_map(|arg| {
        let (key, value) = arg.split_once('=')?;
        if key.trim() != "name" {
            return None;
        }
        let value = value.trim();
        value
            .strip_prefix('"')
            .and_then(|v| v.strip_suffix('"'))
            .map(str::to_owned)
    })
}

/// Name of the first `struct` or `enum` in `text`, skipping comment lines.
fn first_type_name(text: &str) -> Option<String> {
    let mut tokens = text
        .lines()
        .filter(|line| !line.trim_start().starts_with("//"))
        .flat_map(str::split_whitespace);
    while let Some(token) = tokens.next() {
        if token == "struct" || token == "enum" {
            let name = tokens.next()?;
            let end = name
                .find(|c: char| !(c.is_alphanumeric() || c == '_'))
                .unwrap_or(name.len());
            return Some(name[..end].to_owned());
        }
    }
    None
}

/// Every `#[pyclass(name = "...")]` declaration of the crate.
fn named_pyclass_declarations() -> Vec<PyClassDecl> {
    const OPEN: &str = "#[pyclass(";
    let mut declarations = Vec::new();
    for (file, text) in crate_sources() {
        let mut rest = text.as_str();
        while let Some(start) = rest.find(OPEN) {
            let after = &rest[start + OPEN.len()..];
            let Some(end) = after.find(")]") else { break };
            let (args, body) = (&after[..end], &after[end + 2..]);
            if let (Some(python_name), Some(rust_type)) =
                (pyclass_python_name(args), first_type_name(body))
            {
                declarations.push(PyClassDecl {
                    python_name,
                    rust_type,
                    file: file.clone(),
                });
            }
            rest = body;
        }
    }
    declarations
}

/// Last path segment of every type passed to `add_class::<...>()`.
fn registered_class_types() -> Vec<String> {
    const OPEN: &str = "add_class::<";
    let mut registered = Vec::new();
    for (_, text) in crate_sources() {
        let mut rest = text.as_str();
        while let Some(start) = rest.find(OPEN) {
            let after = &rest[start + OPEN.len()..];
            let Some(end) = after.find('>') else { break };
            if let Some(segment) = after[..end].rsplit("::").next() {
                registered.push(segment.trim().to_owned());
            }
            rest = &after[end..];
        }
    }
    registered
}

/// `version = "..."` of the `[project]` table of `pyproject.toml`.
fn pyproject_version() -> String {
    let text = read_crate_file("pyproject.toml");
    let mut in_project = false;
    for line in text.lines() {
        let line = line.trim();
        if line.starts_with('[') {
            in_project = line == "[project]";
        } else if in_project
            && let Some(value) = line.strip_prefix("version")
            && let Some(value) = value.trim_start().strip_prefix('=')
        {
            return value.trim().trim_matches('"').to_owned();
        }
    }
    panic!("pyproject.toml has no [project] version");
}

/// The Python package version must follow the crate version.
#[test]
fn test_pyproject_version_matches_the_crate_version() {
    assert_eq!(
        pyproject_version(),
        env!("CARGO_PKG_VERSION"),
        "pyproject.toml [project] version is out of date"
    );
}

/// `pyo3/extension-module` stops the linker from linking `libpython`, which
/// makes every test binary of the library fail to link. Cargo must therefore
/// not enable it; maturin does through `pyproject.toml`.
#[test]
fn test_extension_module_feature_is_left_to_maturin() {
    let manifest = read_crate_file("Cargo.toml");
    for line in manifest.lines() {
        if line.trim_start().starts_with("pyo3") {
            assert!(
                !line.contains("extension-module"),
                "Cargo.toml must not enable pyo3/extension-module: {line}"
            );
        }
    }
    let pyproject = read_crate_file("pyproject.toml");
    assert!(
        pyproject.contains("features = [\"pyo3/extension-module\"]"),
        "maturin must enable pyo3/extension-module for the wheel"
    );
}

/// The library's own unit tests (`src/**`) must keep running.
///
/// `test = false` in the `[lib]` table drops them from `cargo test`, and a
/// nextest `default-filter` that names this package drops them from
/// `cargo nextest run`; either one turns the whole suite into dead code that is
/// never even compiled. Both switches were set once and hid hundreds of tests.
#[test]
fn test_unit_tests_of_the_library_are_not_disabled() {
    let manifest = read_crate_file("Cargo.toml");
    let mut in_lib_table = false;
    for line in manifest.lines() {
        let line = line.trim();
        if line.starts_with('[') {
            in_lib_table = line == "[lib]";
        } else if in_lib_table {
            let compact: String = line.split_whitespace().collect();
            assert_ne!(
                compact, "test=false",
                "the [lib] table must not set `test = false`: its unit tests would never run"
            );
        }
    }
    let nextest = read_workspace_file(".config/nextest.toml");
    for line in nextest
        .lines()
        .filter(|line| !line.trim_start().starts_with('#'))
    {
        assert!(
            !(line.contains("default-filter") && line.contains("oxiphysics-python")),
            ".config/nextest.toml must not filter oxiphysics-python out of the default run: {line}"
        );
    }
}

/// Every Phase-6 class must be declared with `#[pyclass(name = "...")]` exactly
/// once and its Rust type must be passed to `add_class`, otherwise
/// `oxiphysics.<ClassName>` does not exist after `maturin develop`.
///
/// Runtime `hasattr` assertions live in `python/tests/test_smoke.py`.
#[test]
fn test_phase6_classes_are_declared_and_registered() {
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

    let declarations = named_pyclass_declarations();
    let registered = registered_class_types();
    for class_name in phase6_class_names {
        let declared: Vec<&PyClassDecl> = declarations
            .iter()
            .filter(|decl| decl.python_name == *class_name)
            .collect();
        assert_eq!(
            declared.len(),
            1,
            "Phase-6 class {class_name} must have exactly one #[pyclass(name = \"{class_name}\")] \
             declaration, found {declared:?}"
        );
        let decl = declared[0];
        assert!(
            registered.contains(&decl.rust_type),
            "{class_name} ({}) is declared in {} but never passed to add_class",
            decl.rust_type,
            decl.file
        );
    }
}

/// The scanner itself must understand the attribute forms of the sources.
#[test]
fn test_pyclass_scanner_reads_names_and_types() {
    assert_eq!(
        pyclass_python_name(r#"name = "WingParams", from_py_object"#).as_deref(),
        Some("WingParams")
    );
    assert_eq!(pyclass_python_name("from_py_object"), None);
    assert_eq!(pyclass_python_name(r#"rename = "X""#), None);
    assert_eq!(
        first_type_name("#[derive(Debug)]\n/// A struct of data.\npub struct PyThing<T> {")
            .as_deref(),
        Some("PyThing")
    );
    assert_eq!(
        first_type_name("#[derive(Clone)]\npub enum Mode {\n}").as_deref(),
        Some("Mode")
    );
    assert_eq!(first_type_name("fn nothing() {}"), None);
    let declarations = named_pyclass_declarations();
    assert!(
        declarations.len() >= 30,
        "expected the named pyclass declarations of the crate, found {}",
        declarations.len()
    );
    let registered = registered_class_types();
    assert!(
        registered.len() >= 30,
        "expected the add_class calls of the crate, found {}",
        registered.len()
    );
}
