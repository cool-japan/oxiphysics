# oxiphysics-python TODO

**Version:** 0.1.3 | **Updated:** 2026-06-11 | **Status:** v0.1.x complete — production baseline; forward roadmap active
**Scale:** 26,026 SLoC | 788 unit tests + 600+ pytest

Legend: `[x]` shipped · `[ ]` planned · `[~]` deferred/blocked

## Completed (v0.1.0 – v0.1.3)

### Phase 1: Foundation
- [x] Define core types (`types` module) and error handling (`error` module)
- [x] Implement serde-based JSON bridge (`serialization` module)
- [x] Add unit tests (788 tests passing)

### Phase 2: Domain API Modules
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

### Phase 3: Documentation & Examples
- [x] All 1,200 public items documented
- [x] Serialization round-trip examples

### Phase 4: PyO3 FFI Integration (0.1.1)
- [x] Add `pyo3` dependency to `Cargo.toml`
- [x] Wire `#[pymodule]` entry point (exposes `SimConfig`, `RigidBodyConfig`, `ContactResult`, `PhysicsWorld`)
- [x] Expose `PyPhysicsWorld` as `PhysicsWorld` native Python class (`py_classes.rs`)
- [x] Expose per-domain API structs as `#[pyclass]` (`SimConfig`, `RigidBodyConfig`, `ContactResult`)
- [x] Add maturin build configuration (`pyproject.toml`)

> The former Phase-4 "Publish pip-installable wheel (0.2.0)" item is carried forward as the v0.2.0 "Publishing decision" below.

### Phase 5: Python-Side Ergonomics (0.1.1)
- [x] numpy array bridging for bulk data — `all_positions_flat()`, `all_velocities_flat()`
- [x] asyncio integration for simulation stepping — documented `asyncio.to_thread` pattern in `step()` / `step_substeps()` docstrings
- [x] Python stub files (`.pyi`) for IDE completion — `oxiphysics.pyi`
- [x] Integration tests via pytest — completed 2026-05-11 (pytest harness with conftest.py, test_world.py, test_phase6_authoring.py, test_phase6_utilities.py, test_phase6_dynamics.py, test_numpy_bridge.py, test_pyi_consistency.py; 600 tests collected; wheel: oxiphysics-0.1.1-cp314-cp314-macosx_11_0_arm64.whl)

### Phase 6: Umbrella Module Coverage (shipped 0.1.1–0.1.2)
> **Goal (original, 2026-05-04):** Expose the Phase 13-19 high-level modules from the umbrella `oxiphysics` crate to Python. Achieved 2026-05-06 – 2026-05-11; every module below was planned 2026-05-04 as the Wave 2 annotation pass (add `#[pyclass]`/`#[pymethods]` and `register_*_module`).

#### 6.1 Runtime & Scene (Phase 13-15) — shipped 2026-05-06
- [x] `force_field` — `ForceField` enum + `ForceFieldSystem::apply_to_batch`
- [x] `event_bus` — `PhysicsEvent`, `EventBus` publish/subscribe/drain
- [x] `replay` — `SimRecorder`, `SimReplayer`, `ReplayRecord` (JSON round-trip)
- [x] `query` — `Ray`, `QueryShape`, `QueryWorld::raycast/overlap_sphere/closest_body/k_nearest`
- [x] `scene` — `SceneDescription`, `SceneBuilder` fluent API, JSON round-trip
- [x] `snapshot` — `WorldSnapshot`, `SnapshotDiff`, `SnapshotManager` ring-buffer
- [x] `trigger` — `TriggerVolume`, `TriggerWorld::update` with Enter/Exit/Stay events
- [x] `animation` — `AnimationClip`, `AnimationPlayer` with SLERP + easing
- [x] `material_table` — `MaterialTable` with presets and `CombineRule` resolution

#### 6.2 Physics Utilities (Phase 16-18) — shipped 2026-05-06
- [x] `debug_draw` — `DebugDrawSession`, `DrawCommand`, `DrawList` (useful for Jupyter visualization)
- [x] `contact_cache` — `ContactCache` warm-start impulse lookups
- [x] `buoyancy` — `BuoyancyWorld::apply`, fluid volumes, spherical-cap formula
- [x] `scheduler` — `Scheduler::schedule`, `auto_sleep_step`, priority-based step allocation
- [x] `spatial_grid` — `SpatialGrid::query_radius/query_aabb/k_nearest`
- [x] `lod` — `LodSystem::update`, tier counts, substep allocation
- [x] `noise` — `ValueNoise3D`, `FractalNoise`, `turbulence_force` (exposed as numpy arrays)
- [x] `interpolator` — `smooth_damp`, `exp_decay`, `SpringFollower`, `lerp/smoothstep`
- [x] `telemetry` — `PhysicsStats`, `TelemetrySession`, rolling window export to pandas

#### 6.3 Advanced Dynamics & Tooling (Phase 19) — shipped 2026-05-06
- [x] `character` — `CharacterController::move_and_slide/jump`, capsule kinematic control
- [x] `rope` — `Rope::step`, segment iteration, anchor modes
- [x] `ik` — `IkChain`, `IkSolver::Fabrik/TwoBone`, `SolveReport`
- [x] `xpbd` — `XpbdSolver`, distance/angle/volume constraints with compliance
- [x] `profiler` — `ProfilerSession`, `FrameReport`, to_csv/to_json/to_folded_stacks
- [x] `aero` — `AeroSystem::apply`, wing surfaces, lift/drag curves
- [x] `navmesh` — `NavMesh::from_triangles`, A* path query, funnel smoothing
- [x] `rollback` — `RollbackBuffer`, resimulate-from-tick, desync hash

#### 6.4 Packaging & Validation
- [x] Wire orphan `*_api.rs` modules into `#[pymodule] fn oxiphysics` (done 2026-05-11)
  - **Goal:** Every Phase-6 `#[pyclass]` type in `src/<name>_api.rs` is reachable as `oxiphysics.<ClassName>` from Python.
  - **Design:**
    1. **Audit**: `rg -n '#\[pyclass' crates/oxiphysics-python/src/` + `rg -n 'register_.*_module|m\.add_class'` to find orphan modules (expected ~18: force_field, telemetry, animation, character, buoyancy, scene, replay, rope, aero, event_bus, navmesh, spatial_grid, xpbd, material_table, ik, lod, trigger, profiler).
    2. **For each truly-orphan module**: if `register_<name>_module` exists in `*_api.rs`, call it from `lib.rs::#[pymodule]`; otherwise add `register_*_module` with `m.add_class::<...>()?` for each `#[pyclass]`.
    3. **Module layering**: prefer flat top-level (`oxiphysics.CharacterController`) over sub-module nesting.
    4. **Test**: extend `tests/pyo3_introspection.rs` to enumerate all Phase-6 class names and assert each `m.getattr(name)` succeeds.
  - **Tests:** `test_phase6_classes_introspectable` in `pyo3_introspection.rs` — `Python::with_gil`, import `oxiphysics`, assert `hasattr` for each Phase-6 class.
  - **Risk handled:** duplicate-registration panics avoided by verifying via Python-side import attempts (PyO3 panics on double `m.add_class::<X>()`).
- [x] (2026-05-11) Replaced `(0.2.0)` placeholders with concrete tasks: `pytest` harness in `python/tests/`, `pyproject.toml` wheel build via `maturin build --release --out dist/`, `twine` upload documented in README (NOT executed)
  - **Goal:** Full pytest harness in `python/tests/`; `maturin build --release` produces working wheel; `pypi-publish.yml` verified.
  - **Design:**
    1. Test files: `test_world.py`, `test_phase6_authoring.py`, `test_phase6_utilities.py`, `test_phase6_dynamics.py`, `conftest.py`.
    2. `pyproject.toml`: `[tool.pytest.ini_options] testpaths = ["python/tests", "tests"]`; `[project.optional-dependencies] dev = ["pytest>=7", "numpy>=1.24"]`.
    3. `Makefile` with `dev:` (maturin develop --release) and `test:` (python -m pytest python/tests/ -v) targets.
    4. Wheel verification: `maturin build --release --out dist/`; install into temp venv; run pytest.
  - **Tests:** `maturin develop --release` + `pytest python/tests/ -v`; genuinely-missing classes marked `@pytest.mark.skip(reason="class not yet bound")`.
- [x] (2026-05-11) Numpy bulk-array bridging for noise fields, debug-draw vertex buffers, telemetry time series
  - **Goal:** Three high-volume data paths expose `numpy.ndarray` accessors: `ValueNoise3D::sample_grid_to_numpy`, `DebugDrawSession::vertex_buffer_to_numpy`, `TelemetrySession::time_series_to_numpy`.
  - **Design:**
    1. Workspace `numpy` dep matched to pyo3 0.28, behind default `numpy-bridge` feature (`[features] default = ["numpy-bridge"], numpy-bridge = ["dep:numpy"]`).
    2. `noise_api.rs`: `sample_grid_to_numpy(...) -> PyResult<&PyArray3<f64>>`; `viz_api.rs`: `vertex_buffer_to_numpy(...) -> &PyArray2<f64>`; `telemetry_api.rs`: `time_series_to_numpy(...) -> PyResult<&PyArray1<f64>>`.
  - **Tests:** `test_value_noise_to_numpy_shape_dtype`, `test_debug_draw_zero_copy_view`, `test_telemetry_time_series_length_match` (`python/tests/test_numpy_bridge.py`).
- [x] `.pyi` stubs extended to cover all Phase 6 classes
  - **Goal:** `python/oxiphysics/__init__.pyi` and root `oxiphysics.pyi` declare every Phase-6 class with full method/attribute signatures and return types.
  - **Design:**
    1. Signatures derived from `#[pymethods] impl X` blocks: `&str`→`str`, `f64`→`float`, `Vec<f64>`→`list[float]` (or `numpy.ndarray` where numpy accessors exist), `(f64,f64,f64)`→`tuple[float,float,float]`.
    2. `python/oxiphysics/__init__.pyi` extended from 153 lines to ~1500 lines, grouped under Phase 6.1/6.2/6.3 section banners; root `oxiphysics.pyi` regenerated to match.
  - **Tests:** `python/tests/test_pyi_consistency.py` — for each class in `__init__.pyi`, assert it exists in the runtime `oxiphysics` module.

### v0.1.2 correctness fixes (2026-06-01)
- [x] `PyCsg::union/intersection/subtraction` — replaced AABB-approximation stubs with delegation to `oxiphysics::geometry::mesh_boolean::mesh_boolean`; proper winding-number inside-outside classification + cleanup. 7 new integration tests in `tests/csg_imls_algorithms.rs`.
- [x] `PyPointCloud::poisson_reconstruct` — replaced empty-mesh stub with real IMLS (Implicit Moving Least Squares) reconstruction: PCA normal estimation, Gaussian-weighted tangent-plane signed distance, marching-cubes isosurface via `signed_distance_field::MarchingCubes`.

## v0.2.0 — RL + typed surface

> Policy note: gymnasium, torch/jax (dlpack consumers), and pyvista are optional-extra **Python-side** dependencies installed via pip extras — they are not Rust crate dependencies, so the Pure-Rust default build is unaffected.

- [ ] gymnasium-API RL environments (CartPole, Ant, Reacher from articulated+rigid)
  - **Goal:** `gymnasium.make("OxiCartPole-v0")` works; PPO (SB3) reaches reward >475 in <5 min CPU; `env.step()` overhead <50 µs.
  - **Design:** pure-Python `oxiphysics.gym` package layered on existing bindings; Ant uses MJCF import when available, hand-built until then.
- [ ] dlpack zero-copy tensor exchange
  - **Goal:** `torch.from_dlpack(state)` and JAX equivalent share memory (pointer-equality test, no copy) for positions/velocities buffers.
  - **Design:** `__dlpack__`/`__dlpack_device__` on buffer-owning pyclasses next to the existing numpy bridge.
- [ ] mypy --strict gate on stubs
  - **Goal:** `mypy --strict` clean on stubs + examples via local script (`test_pyi_consistency.py` exists — extend); 100% of 214 pyclasses covered.
- [ ] pip extras `[viz]` / `[gym]` / `[io]`
  - **Goal:** `pip install oxiphysics[gym]` pulls gymnasium and registers envs; extras documented in pyproject.
- [ ] Overhead benchmark vs pybullet / mujoco bindings
  - **Goal:** published pytest-benchmark table (step/raycast/state-read on a 100-body scene); step overhead within 2x pybullet; regenerated by one script.
- [ ] pandas telemetry export + matplotlib conveniences
  - **Goal:** `TelemetrySession.to_dataframe()` with correct dtypes; 3 example notebooks (energy trend, contact stats, profiler).
- [ ] Publishing decision (carries the Phase-4 "Publish pip-installable wheel" item; see root TODO Phase 23.6 publish xtask)
  - **Goal:** crate is `publish = false`; either wire full PyPI distribution through the existing `pypi-publish.yml` (allowed) with abi-tagged wheels for 3 platforms, or document wheel-only distribution. Must land before the publish xtask assumes order.
- [ ] Drop direct nalgebra dep (root TODO Phase 23.9 policy reconciliation)
  - **Goal:** nalgebra absent from `crates/oxiphysics-python/Cargo.toml` (confirmed present today); route through oxiphysics-core types.

## v0.3.0 — Ecosystem & docs

- [ ] pyvista exporter
  - **Goal:** FEM stress field plotted via `pyvista` in ≤5 lines (UnstructuredGrid from mesh+field accessors); matplotlib quiver/scatter helpers for SPH/MD.
- [ ] Jupyter inline viewer widget
  - **Goal:** interactive orbit of a 1k-body scene in JupyterLab (anywidget streaming PNG frames from viz offscreen renderer, or the wasm canvas).
- [ ] sphinx-gallery docs site
  - **Goal:** 20 executed gallery examples build warning-free via local script; deploy manual (no new workflow).
- [ ] asyncio driver hardening
  - Verified: async support exists in `py_classes.rs`; remaining work is cancellation + progress callbacks.
  - **Goal:** `asyncio.run` drives 10k steps with cooperative yields and clean Ctrl-C.
- [ ] Vectorized batch-env API
  - **Goal:** 1024 cartpole envs stepped in one call (rayon inside, GIL released) at ≥1M aggregate steps/s; numpy `(N, ...)` state out.

## v1.0 — Production & validation

- [ ] Stable 1.0 Python API
  - **Goal:** abi3-py39 wheels pass the test matrix on py3.9–3.13; deprecation policy documented; 100% docstring coverage gate in local script.

## Deferred / research track

- [~] MJCF-backed Ant/humanoid gym environments — Ready when: oxiphysics-io MJCF import subset (io v0.2.0) lands; hand-built Ant ships in the meantime.
- [~] h5py interop assertion (open an OxiPhysics-written HDF5 file from this crate's pytest suite) — Ready when: oxiphysics-io HDF5 on-disk subset (io v0.3.0) lands.
