# oxiphysics-python TODO

Updated: 2026-06-06 | Version: 0.1.3

---

## Phase 1: Foundation ✅

- [x] Define core types (`types` module) and error handling (`error` module)
- [x] Implement serde-based JSON bridge (`serialization` module)
- [x] Add unit tests (788 tests passing)

## Phase 2: Domain API Modules ✅

- [x] `analytics_api` — analytics/telemetry query types
- [x] `constraints_api` — constraint/joint configuration
- [x] `fem_api` — FEM mesh and solver types (`PyFemSolver`, `PyFemMesh`)
- [x] `geometry_api` — geometry primitives
- [x] `io_api` — scene import/export bridge
- [x] `lbm_api` — LBM configuration (`PyLbmConfig`, `PyLbmSimulation`)
- [x] `materials_api` — material model parameter types
- [x] `md_api` — molecular dynamics types (`PyMdConfig`, `PyMdSimulation`)
- [x] `rigid_api` — rigid body parameters
- [x] `sph_api` — SPH configuration (`PySphConfig`, `PySphSimulation`)
- [x] `vehicle_api` — vehicle dynamics configuration
- [x] `viz_api` — visualization output descriptors
- [x] `world_api` — top-level world/scene types (`PyPhysicsWorld`)

## Phase 3: Documentation & Examples ✅

- [x] All 1,200 public items documented
- [x] Serialization round-trip examples

## Phase 4: PyO3 FFI Integration ✅ (0.1.1)

- [x] Add `pyo3` dependency to `Cargo.toml`
- [x] Wire `#[pymodule]` entry point (exposes `SimConfig`, `RigidBodyConfig`, `ContactResult`, `PhysicsWorld`)
- [x] Expose `PyPhysicsWorld` as `PhysicsWorld` native Python class (`py_classes.rs`)
- [x] Expose per-domain API structs as `#[pyclass]` (`SimConfig`, `RigidBodyConfig`, `ContactResult`)
- [x] Add maturin build configuration (`pyproject.toml`)
- [ ] Publish pip-installable wheel (0.2.0)

## Phase 5: Python-Side Ergonomics ✅ (0.1.1)

- [x] numpy array bridging for bulk data — `all_positions_flat()`, `all_velocities_flat()`
- [x] asyncio integration for simulation stepping — documented `asyncio.to_thread` pattern in `step()` / `step_substeps()` docstrings
- [x] Python stub files (`.pyi`) for IDE completion — `oxiphysics.pyi`
- [x] Integration tests via pytest (0.2.0) — completed 2026-05-11 (pytest harness with conftest.py, test_world.py, test_phase6_authoring.py, test_phase6_utilities.py, test_phase6_dynamics.py, test_numpy_bridge.py, test_pyi_consistency.py; 600 tests collected; wheel: oxiphysics-0.1.1-cp314-cp314-macosx_11_0_arm64.whl)

## Phase 6: Umbrella Module Coverage (0.2.0)

> **Goal:** Expose the Phase 13-19 high-level modules from the umbrella `oxiphysics` crate to Python. Current bindings stop at the Phase 1-12 domain APIs; nothing beyond is reachable from Python.

### 6.1 Runtime & Scene (Phase 13-15)
- [x] (2026-05-06) `force_field` — `ForceField` enum + `ForceFieldSystem::apply_to_batch`
  - **Planned:** 2026-05-04 — Wave 2 annotation pass (add `#[pyclass]`/`#[pymethods]` and `register_*_module`)
- [x] (2026-05-06) `event_bus` — `PhysicsEvent`, `EventBus` publish/subscribe/drain
  - **Planned:** 2026-05-04 — Wave 2 annotation pass (add `#[pyclass]`/`#[pymethods]` and `register_*_module`)
- [x] (2026-05-06) `replay` — `SimRecorder`, `SimReplayer`, `ReplayRecord` (JSON round-trip)
  - **Planned:** 2026-05-04 — Wave 2 annotation pass (add `#[pyclass]`/`#[pymethods]` and `register_*_module`)
- [x] (2026-05-06) `query` — `Ray`, `QueryShape`, `QueryWorld::raycast/overlap_sphere/closest_body/k_nearest`
  - **Planned:** 2026-05-04 — Wave 2 annotation pass (add `#[pyclass]`/`#[pymethods]` and `register_*_module`)
- [x] (2026-05-06) `scene` — `SceneDescription`, `SceneBuilder` fluent API, JSON round-trip
  - **Planned:** 2026-05-04 — Wave 2 annotation pass (add `#[pyclass]`/`#[pymethods]` and `register_*_module`)
- [x] (2026-05-06) `snapshot` — `WorldSnapshot`, `SnapshotDiff`, `SnapshotManager` ring-buffer
  - **Planned:** 2026-05-04 — Wave 2 annotation pass (add `#[pyclass]`/`#[pymethods]` and `register_*_module`)
- [x] (2026-05-06) `trigger` — `TriggerVolume`, `TriggerWorld::update` with Enter/Exit/Stay events
  - **Planned:** 2026-05-04 — Wave 2 annotation pass (add `#[pyclass]`/`#[pymethods]` and `register_*_module`)
- [x] (2026-05-06) `animation` — `AnimationClip`, `AnimationPlayer` with SLERP + easing
  - **Planned:** 2026-05-04 — Wave 2 annotation pass (add `#[pyclass]`/`#[pymethods]` and `register_*_module`)
- [x] (2026-05-06) `material_table` — `MaterialTable` with presets and `CombineRule` resolution
  - **Planned:** 2026-05-04 — Wave 2 annotation pass (add `#[pyclass]`/`#[pymethods]` and `register_*_module`)

### 6.2 Physics Utilities (Phase 16-18)
- [x] (2026-05-06) `debug_draw` — `DebugDrawSession`, `DrawCommand`, `DrawList` (useful for Jupyter visualization)
  - **Planned:** 2026-05-04 — Wave 2 annotation pass (add `#[pyclass]`/`#[pymethods]` and `register_*_module`)
- [x] (2026-05-06) `contact_cache` — `ContactCache` warm-start impulse lookups
  - **Planned:** 2026-05-04 — Wave 2 annotation pass (add `#[pyclass]`/`#[pymethods]` and `register_*_module`)
- [x] (2026-05-06) `buoyancy` — `BuoyancyWorld::apply`, fluid volumes, spherical-cap formula
  - **Planned:** 2026-05-04 — Wave 2 annotation pass (add `#[pyclass]`/`#[pymethods]` and `register_*_module`)
- [x] (2026-05-06) `scheduler` — `Scheduler::schedule`, `auto_sleep_step`, priority-based step allocation
  - **Planned:** 2026-05-04 — Wave 2 annotation pass (add `#[pyclass]`/`#[pymethods]` and `register_*_module`)
- [x] (2026-05-06) `spatial_grid` — `SpatialGrid::query_radius/query_aabb/k_nearest`
  - **Planned:** 2026-05-04 — Wave 2 annotation pass (add `#[pyclass]`/`#[pymethods]` and `register_*_module`)
- [x] (2026-05-06) `lod` — `LodSystem::update`, tier counts, substep allocation
  - **Planned:** 2026-05-04 — Wave 2 annotation pass (add `#[pyclass]`/`#[pymethods]` and `register_*_module`)
- [x] (2026-05-06) `noise` — `ValueNoise3D`, `FractalNoise`, `turbulence_force` (exposed as numpy arrays)
  - **Planned:** 2026-05-04 — Wave 2 annotation pass (add `#[pyclass]`/`#[pymethods]` and `register_*_module`)
- [x] (2026-05-06) `interpolator` — `smooth_damp`, `exp_decay`, `SpringFollower`, `lerp/smoothstep`
  - **Planned:** 2026-05-04 — Wave 2 annotation pass (add `#[pyclass]`/`#[pymethods]` and `register_*_module`)
- [x] (2026-05-06) `telemetry` — `PhysicsStats`, `TelemetrySession`, rolling window export to pandas
  - **Planned:** 2026-05-04 — Wave 2 annotation pass (add `#[pyclass]`/`#[pymethods]` and `register_*_module`)

### 6.3 Advanced Dynamics & Tooling (Phase 19)
- [x] (2026-05-06) `character` — `CharacterController::move_and_slide/jump`, capsule kinematic control
  - **Planned:** 2026-05-04 — Wave 2 annotation pass (add `#[pyclass]`/`#[pymethods]` and `register_*_module`)
- [x] (2026-05-06) `rope` — `Rope::step`, segment iteration, anchor modes
  - **Planned:** 2026-05-04 — Wave 2 annotation pass (add `#[pyclass]`/`#[pymethods]` and `register_*_module`)
- [x] (2026-05-06) `ik` — `IkChain`, `IkSolver::Fabrik/TwoBone`, `SolveReport`
  - **Planned:** 2026-05-04 — Wave 2 annotation pass (add `#[pyclass]`/`#[pymethods]` and `register_*_module`)
- [x] (2026-05-06) `xpbd` — `XpbdSolver`, distance/angle/volume constraints with compliance
  - **Planned:** 2026-05-04 — Wave 2 annotation pass (add `#[pyclass]`/`#[pymethods]` and `register_*_module`)
- [x] (2026-05-06) `profiler` — `ProfilerSession`, `FrameReport`, to_csv/to_json/to_folded_stacks
  - **Planned:** 2026-05-04 — Wave 2 annotation pass (add `#[pyclass]`/`#[pymethods]` and `register_*_module`)
- [x] (2026-05-06) `aero` — `AeroSystem::apply`, wing surfaces, lift/drag curves
  - **Planned:** 2026-05-04 — Wave 2 annotation pass (add `#[pyclass]`/`#[pymethods]` and `register_*_module`)
- [x] (2026-05-06) `navmesh` — `NavMesh::from_triangles`, A* path query, funnel smoothing
  - **Planned:** 2026-05-04 — Wave 2 annotation pass (add `#[pyclass]`/`#[pymethods]` and `register_*_module`)
- [x] (2026-05-06) `rollback` — `RollbackBuffer`, resimulate-from-tick, desync hash
  - **Planned:** 2026-05-04 — Wave 2 annotation pass (add `#[pyclass]`/`#[pymethods]` and `register_*_module`)

### 6.4 Packaging & Validation
- [x] Wire orphan `*_api.rs` modules into `#[pymodule] fn oxiphysics` (done 2026-05-11)
  - **Goal:** Every Phase-6 `#[pyclass]` type in `src/<name>_api.rs` is reachable as `oxiphysics.<ClassName>` from Python.
  - **Design:**
    1. **Audit**: `rg -n '#\[pyclass' crates/oxiphysics-python/src/` + `rg -n 'register_.*_module|m\.add_class'` to find orphan modules (expected ~18: force_field, telemetry, animation, character, buoyancy, scene, replay, rope, aero, event_bus, navmesh, spatial_grid, xpbd, material_table, ik, lod, trigger, profiler).
    2. **For each truly-orphan module**: if `register_<name>_module` exists in `*_api.rs`, call it from `lib.rs::#[pymodule]`; otherwise add `register_*_module` with `m.add_class::<...>()?` for each `#[pyclass]`.
    3. **Module layering**: prefer flat top-level (`oxiphysics.CharacterController`) over sub-module nesting.
    4. **Test**: extend `tests/pyo3_introspection.rs` to enumerate all Phase-6 class names and assert each `m.getattr(name)` succeeds.
  - **Files:** `src/lib.rs`; each truly-orphan `*_api.rs`; `tests/pyo3_introspection.rs`; `TODO.md`
  - **Tests:** `test_phase6_classes_introspectable` in `pyo3_introspection.rs` — `Python::with_gil`, import `oxiphysics`, assert `hasattr` for each Phase-6 class. Failure prints which name is missing.
  - **Risk:** Some candidates may already be registered via sub-module paths; subagent verifies via Python-side import attempts to avoid duplicate-registration panics (PyO3 panics on double `m.add_class::<X>()`).
- [x] (2026-05-11) Replace `(0.2.0)` placeholders above with concrete tasks: `pytest` harness in `python/tests/`, `pyproject.toml` wheel build via `maturin build --release --out dist/`, upload via `twine` to PyPI after tag
  - **Goal:** Full pytest harness in `python/tests/`; `maturin build --release` produces working wheel; `pypi-publish.yml` verified.
  - **Design:**
    1. Add test files: `test_world.py`, `test_phase6_authoring.py` (MaterialTable, SceneDescription, SceneBuilder, SimRecorder, SimReplayer, ReplayRecord), `test_phase6_utilities.py` (DebugDrawSession, BuoyancyWorld, Scheduler, SpatialGrid, LodSystem, ValueNoise3D, FractalNoise, TelemetrySession, ContactCache), `test_phase6_dynamics.py` (CharacterController, Rope, IkSolver, XpbdSolver, AeroSystem, NavMesh, TriggerWorld, EventBus, AnimationPlayer, ProfilerSession, ForceFieldSystem), `conftest.py`.
    2. `pyproject.toml`: add `[tool.pytest.ini_options] testpaths = ["python/tests", "tests"]`; add `[project.optional-dependencies] dev = ["pytest>=7", "numpy>=1.24"]`.
    3. Create `Makefile` with `dev:` (maturin develop --release) and `test:` (python -m pytest python/tests/ -v) targets.
    4. Wheel verification: `maturin build --release --out dist/`; install into temp venv; run pytest.
    5. Document `twine upload dist/*` workflow in `README.md` (NOT executed in this run).
  - **Files:** `python/tests/test_world.py`, `test_phase6_authoring.py`, `test_phase6_utilities.py`, `test_phase6_dynamics.py`, `conftest.py`; `pyproject.toml`; `Makefile`; `README.md`; `TODO.md`
  - **Tests:** `maturin develop --release` + `pytest python/tests/ -v`; genuinely-missing classes marked `@pytest.mark.skip(reason="class not yet bound")`.
  - **Risk:** `maturin` CLI must be installed (auto-install if missing via `pip install maturin`); wheel build takes 60-120s (use debug mode for iteration).
- [x] (2026-05-11) Numpy bulk-array bridging for noise fields, debug-draw vertex buffers, telemetry time series
  - **Goal:** Three high-volume data paths expose `numpy.ndarray` accessors: `ValueNoise3D::sample_grid_to_numpy`, `DebugDrawSession::vertex_buffer_to_numpy`, `TelemetrySession::time_series_to_numpy`.
  - **Design:**
    1. **Workspace dep**: add `numpy = "0.28"` to `[workspace.dependencies]` (check crates.io for exact latest version compatible with pyo3 0.28).
    2. **oxiphysics-python Cargo.toml**: add `numpy = { workspace = true, optional = true }`; add `[features] default = ["numpy-bridge"], numpy-bridge = ["dep:numpy"]`.
    3. **`ValueNoise3D`/`FractalNoise`** in `noise_api.rs`: `sample_grid_to_numpy(&self, py: Python<'_>, origin, step, shape) -> PyResult<&PyArray3<f64>>`.
    4. **`DebugDrawSession`** in `viz_api.rs`: `vertex_buffer_to_numpy(&self, py: Python<'_>) -> &PyArray2<f64>`.
    5. **`TelemetrySession`** in `telemetry_api.rs`: `time_series_to_numpy(&self, py: Python<'_>, field: &str) -> PyResult<&PyArray1<f64>>`.
  - **Files:** root `Cargo.toml`; `crates/oxiphysics-python/Cargo.toml`; `src/noise_api.rs`; `src/viz_api.rs`; `src/telemetry_api.rs`; `python/oxiphysics/__init__.pyi`; `python/tests/test_numpy_bridge.py`; `TODO.md`
  - **Tests:** `test_value_noise_to_numpy_shape_dtype`, `test_debug_draw_zero_copy_view`, `test_telemetry_time_series_length_match`
  - **Risk:** numpy crate version must match pyo3 version; `IntoPyArray` API may differ from 0.27; zero-copy semantics vs move semantics for `Vec<f64>::into_pyarray`.
## v0.1.2 correctness fixes (2026-06-01)
- [x] `PyCsg::union/intersection/subtraction` — replaced AABB-approximation stubs with delegation to `oxiphysics::geometry::mesh_boolean::mesh_boolean`; proper winding-number inside-outside classification + cleanup. 7 new integration tests in `tests/csg_imls_algorithms.rs`.
- [x] `PyPointCloud::poisson_reconstruct` — replaced empty-mesh stub with real IMLS (Implicit Moving Least Squares) reconstruction: PCA normal estimation, Gaussian-weighted tangent-plane signed distance, marching-cubes isosurface via `signed_distance_field::MarchingCubes`.

- [x] `.pyi` stubs extended to cover all Phase 6 classes
  - **Goal:** `python/oxiphysics/__init__.pyi` and root `oxiphysics.pyi` declare every Phase-6 class with full method/attribute signatures and return types.
  - **Design:**
    1. Derive signatures from `#[pymethods] impl X` blocks: `&str`→`str`, `f64`→`float`, `Vec<f64>`→`list[float]` (or `numpy.ndarray` where block 6 adds numpy accessors), `(f64,f64,f64)`→`tuple[float,float,float]`.
    2. Extend `python/oxiphysics/__init__.pyi` from 153 lines to ~1500 lines; group under Phase 6.1/6.2/6.3 section banners.
    3. Regenerate root `oxiphysics.pyi` to match.
    4. Subagent may use a one-off `/tmp/` script to draft stubs from `*_api.rs`; script NOT committed.
  - **Files:** `python/oxiphysics/__init__.pyi`; `oxiphysics.pyi`; `TODO.md`
  - **Tests:** `python/tests/test_pyi_consistency.py` — for each class in `__init__.pyi`, assert exists in runtime `oxiphysics` module; `mypy --strict` if available.
  - **Risk:** Stub drift; pragmatically prioritize coverage breadth over type depth in this pass.
