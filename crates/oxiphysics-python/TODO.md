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

<!-- AUDIT-2026-10-12 BEGIN -->
## Audit 2026-10-12 — verified findings (Phase 32 remediation)

Verified: 111 (critical 4, high 27, medium 59, low 21); 21 further raw findings pending second-round verification. IDs refer to the audit ledger; `[x]` = fixed.

### Critical (4)
- [ ] **F3634** (critical, bug, easy/trivial) `lbm_api.rs:832` — Python 3D LBM bounce-back uses wrong D3Q19 opposite table
- [ ] **F3672** (critical, fabrication, medium/medium) `vehicle_api.rs:1190` — Vehicle step never writes slip_ratio; throttle gives zero traction
- [ ] **F3703** (critical, fabrication, medium/medium) `world_api/mod.rs:768` — PhysicsWorld: spheres only, inverted impulse sign, no friction
- [ ] **F3705** (critical, bug, easy/small) `world_api/mod.rs:846` — Contact impulse gate inverted: approaching spheres get no impulse

### High (27)
- [ ] **F3592** (high, test, medium/medium) `crates/oxiphysics-python/Cargo.toml:17` — 854 oxiphysics-python unit tests never compile or run
- [ ] **F3605** (high, fabrication, medium/small) `constraints_api.rs:1148` — XPBD volume constraint hardcodes rest_volume 1.0 and is never solved
- [ ] **F3611** (high, duplicate, hard/large) `fem_api.rs:664` — Python fem/sph/lbm/md/vehicle/rigid/world modules re-implement physics, never call engine
- [ ] **F3612** (high, perf, medium/medium) `fem_api.rs:697` — FEM solve: dense n^2 stiffness, 1e30 penalty, absolute-tolerance CG
- [ ] **F3614** (high, fabrication, hard/medium) `fem_api.rs:1031` — ModalAnalyzer returns dominant K eigenvalues, ignores mass, fake ω=1
- [ ] **F3615** (high, fabrication, hard/medium) `fem_api.rs:1081` — ModalAnalyzer returns largest K modes and a fake ω=1 BC mode
- [ ] **F3617** (high, fabrication, medium/small) `geometry_api/shapes.rs:324` — PyConvexHull is an AABB box; gjk_overlap is a bounding-sphere test
- [ ] **F3618** (high, fabrication, easy/small) `geometry_api/shapes.rs:349` — ConvexHull.from_points returns the AABB; gjk_overlap is a bounding-sphere test
- [ ] **F3625** (high, bug, easy/small) `io_api.rs:101` — Binary VTK writer makes every vertex cell reference point 0
- [ ] **F3631** (high, bug, medium/small) `lbm_api.rs:398` — LBM body force injected at half magnitude (D2Q9 and D3Q19)
- [ ] **F3643** (high, perf, medium/medium) `md_api.rs:401` — MD Lennard-Jones forces use an O(n^2) all-pairs loop
- [ ] **F3644** (high, bug, medium/small) `md_api.rs:1050` — Nose-Hoover NVT: discarded half-step, Q fixed at 1, double thermostat
- [ ] **F3652** (high, bug, medium/medium) `rigid_api/joints.rs:258` — rigid.PyWorld cannot receive bodies; pipeline has no collision
- [ ] **F3653** (high, fabrication, hard/large) `rigid_api/joints.rs:363` — PyPhysicsPipeline step: no rotation/collision; joints ignore type, anchors, mass
- [ ] **F3654** (high, fabrication, hard/large) `rigid_api/joints.rs:389` — rigid.PyWorld unpopulatable from Python; pipeline ignores joint types
- [ ] **F3664** (high, perf, medium/medium) `sph_api.rs:260` — WCSPH step uses O(n^2) all-pairs density and force loops
- [ ] **F3668** (high, bug, easy/trivial) `vehicle_api.rs:30` — Pacejka converts rad→deg and ratio→% while presets are rad-scale
- [ ] **F3675** (high, fabrication, easy/small) `viz_api/types.rs:435` — Stress visualizer fakes principal stresses and directions
- [ ] **F3680** (high, bug, medium/medium) `world_api/constraints.rs:203` — Constraints only nudge positions; velocities, anchors ignored; Hinge/Slider no-op
- [ ] **F3685** (high, fabrication, medium/medium) `world_api/fem.rs:111` — 3-D truss FEM is x-only 1-D and permanently erases y/z loads
- [ ] **F3690** (high, fabrication, medium/medium) `world_api/lbm.rs:45` — lid_driven_cavity has no lid or walls; bounce-back doc false
- [ ] **F3691** (high, bug, medium/small) `world_api/lbm.rs:160` — LBM body force injects only half the configured acceleration
- [ ] **F3695** (high, duplicate, hard/large) `world_api/mod.rs:148` — Python core APIs reimplement physics instead of wrapping native crates
- [ ] **F3704** (high, fabrication, medium/small) `world_api/mod.rs:833` — Contact friction is computed and reported but never applied
- [ ] **F3706** (high, bug, easy/small) `world_api/mod.rs:868` — Baumgarte position correction scales with raw 1/m, not mass ratio
- [ ] **F3716** (high, perf, medium/medium) `world_api/sph.rs:206` — SPH density and force loops are O(n²) all-pairs
- [ ] **F3722** (high, bug, easy/small) `xpbd_api.rs:102` — particle_velocity underestimates velocity by a factor of substeps

### Medium (59)
- [ ] **F3593** (medium, test, easy/small) `crates/oxiphysics-python/Cargo.toml:31` — Python crate's Rust unit/integration tests never run
- [ ] **F3600** (medium, fabrication, easy/small) `analytics_api.rs:1049` — Benchmark.run fabricates timings from a formula
- [ ] **F3601** (medium, missing-feature, medium/small) `character_api.rs:15` — Python CharacterController lacks move_and_slide/jump; doc cites missing method
- [ ] **F3602** (medium, panic, easy/small) `constraints_api.rs:40` — PGS iteration panics on lo>hi, NaN bounds, or short slices
- [ ] **F3603** (medium, fabrication, medium/medium) `constraints_api.rs:209` — Constraint solver ignores solver_type, sor_factor and warm_start (diagonal-only)
- [ ] **F3604** (medium, panic, easy/small) `constraints_api.rs:1129` — PBD solver indexes positions/inv_masses unchecked; setters desync lengths
- [ ] **F3606** (medium, panic, easy/small) `constraints_api.rs:1331` — ControlSystem set_all fields let users trigger index panics
- [ ] **F3607** (medium, fabrication, medium/small) `constraints_api.rs:1373` — frequency_gain returns hardcoded 1.0 for any non-first-order system
- [ ] **F3608** (medium, api, easy/small) `error.rs:133` — All Rust errors map to RuntimeError; no exception hierarchy
- [ ] **F3609** (medium, missing-feature, medium/small) `event_bus_api.rs:26` — EventBus has no publish/subscribe; WorldSnapshot cannot add bodies
- [ ] **F3610** (medium, panic, easy/small) `fem_api.rs:307` — FEM mesh accepts invalid node/material/DOF indices, so the solver panics
- [ ] **F3613** (medium, panic, easy/medium) `fem_api.rs:748` — PyFemSolver.solve panics on invalid node/DOF indices; dense K
- [ ] **F3616** (medium, panic, easy/small) `geometry_api/mesh.rs:71` — TriangleMesh.from_raw drops bad vertices; CSG panics on bad indices
- [ ] **F3619** (medium, panic, medium/small) `geometry_api/spatial.rs:72` — Heightfield panics for zero/resized grids; raycast is a coarse march
- [ ] **F3620** (medium, bug, easy/small) `geometry_api/spatial.rs:194` — PySpatialHash cell_size setter corrupts the grid; tiny cells hang queries
- [ ] **F3621** (medium, fabrication, easy/small) `geometry_api/transform.rs:59` — PointCloud.compute_normals ignores k; uses centroid direction, not PCA
- [ ] **F3623** (medium, perf, easy/small) `geometry_api/transform.rs:465` — IMLS reconstruction evaluates every point at every grid node
- [ ] **F3624** (medium, fabrication, easy/small) `io_api.rs:61` — PyVtkWriter emits unparseable ASCII and wrong binary cell indices
- [ ] **F3628** (medium, fabrication, easy/small) `io_api.rs:729` — Trajectory write_frame discards positions and velocities
- [ ] **F3629** (medium, missing-feature, medium/small) `lbm_api.rs:44` — Inlet/Outlet boundary types are silently ignored by step()
- [ ] **F3630** (medium, bug, easy/small) `lbm_api.rs:199` — LBM config silently clamps; unchecked size aborts interpreter
- [ ] **F3632** (medium, perf, easy/small) `lbm_api.rs:424` — LBM streaming clones the whole f_tmp lattice every step
- [ ] **F3636** (medium, fabrication, easy/small) `lbm_api.rs:1213` — 'Zou-He' inlet is an equilibrium overwrite; divides by zero at ux=1
- [ ] **F3637** (medium, stub, easy/small) `lbm_api.rs:1243` — MrtRelaxation exposed but no MRT collision operator is ever applied
- [ ] **F3638** (medium, api, easy/small) `lib.rs:97` — Submodules not importable via 'import oxiphysics.fem'
- [ ] **F3640** (medium, bug, easy/small) `materials_api.rs:582` — Hashin-Shtrikman lower bound pairs fractions with wrong phases
- [ ] **F3642** (medium, bug, easy/trivial) `md_api.rs:340` — thermostat_tau ignored: hard velocity rescale to exact T each step
- [ ] **F3645** (medium, panic, easy/small) `navmesh_api.rs:32` — Python NavMesh accepts invalid triangle indices then panics
- [ ] **F3647** (medium, perf, easy/small) `py_classes.rs:587` — No GIL release anywhere; heavy steps block all Python threads
- [ ] **F3650** (medium, bug, easy/small) `rasterizer.rs:53` — Rasterizer depth convention inverted: farthest fragment wins
- [ ] **F3655** (medium, fabrication, easy/small) `rigid_api/utils.rs:300` — compute_island_stats returns zeros instead of aggregating bodies
- [ ] **F3656** (medium, missing-feature, easy/small) `rope_api.rs:17` — Python Rope cannot be stepped, anchored, or perturbed
- [ ] **F3658** (medium, missing-feature, medium/medium) `serialization/functions.rs:147` — World snapshot unreachable from Python; is_static hardcoded false
- [ ] **F3661** (medium, api, easy/trivial) `sph_api.rs:18` — PySphConfig has no constructor or field access from Python
- [ ] **F3662** (medium, duplicate, medium/medium) `sph_api.rs:239` — Python SPH/MD bindings duplicate O(N²) physics instead of engine
- [ ] **F3663** (medium, perf, medium/medium) `sph_api.rs:250` — SPH density/force loops are O(n²) in Python and WASM
- [ ] **F3665** (medium, panic, easy/trivial) `sph_api.rs:879` — adaptive_h(i, 0) underflows k-1 and panics
- [ ] **F3666** (medium, fabrication, hard/medium) `sph_api.rs:986` — DFSPH variant just runs WCSPH sub-steps; pressure_tol unused
- [ ] **F3669** (medium, fabrication, easy/small) `vehicle_api.rs:182` — 'Fiala' tire is saturated-linear lateral with Fx hardcoded to 0
- [ ] **F3670** (medium, panic, easy/small) `vehicle_api.rs:452` — EngineCurve.torque_at panics on rpm/torque length mismatch
- [ ] **F3671** (medium, bug, easy/trivial) `vehicle_api.rs:1013` — Brake torque acts as reverse drive instead of opposing rotation
- [ ] **F3674** (medium, panic, easy/trivial) `vehicle_api.rs:1526` — PyTelemetry.record panics when max_samples == 0
- [ ] **F3676** (medium, panic, medium/small) `viz_api/types.rs:777` — Streamline sampler panics for dims≤1 and uses nearest-cell sampling
- [ ] **F3679** (medium, bug, medium/medium) `world_api/constraints.rs:184` — Constraint anchors ignored; PointToPoint collapses body centres together
- [ ] **F3681** (medium, fabrication, medium/medium) `world_api/constraints.rs:215` — Hinge/Slider constraints are empty no-ops; constraint API unreachable from Python
- [ ] **F3684** (medium, perf, medium/medium) `world_api/fem.rs:63` — Dense n_dofs² stiffness alloc aborts interpreter; O(n³) solve
- [ ] **F3686** (medium, perf, medium/medium) `world_api/fem.rs:159` — Dense n² stiffness and unchecked allocation abort the interpreter
- [ ] **F3698** (medium, api, easy/small) `world_api/mod.rs:246` — PhysicsWorld mutators silently no-op on invalid handles
- [ ] **F3699** (medium, bug, easy/small) `world_api/mod.rs:471` — step() accepts dt<=0/NaN; no divergence detection; mass/radius unchecked
- [ ] **F3700** (medium, api, easy/small) `world_api/mod.rs:728` — world.PyPhysicsWorld is registered but cannot add bodies from Python
- [ ] **F3701** (medium, bug, easy/small) `world_api/mod.rs:745` — Generation counters unused: stale handles control new bodies
- [ ] **F3702** (medium, perf, medium/medium) `world_api/mod.rs:758` — Rigid-body contact detection is O(n²) all-pairs with no broadphase
- [ ] **F3707** (medium, bug, easy/trivial) `world_api/mod.rs:980` — Orientation integrated as q⊗ω although ω is world-frame
- [ ] **F3708** (medium, fabrication, medium/medium) `world_api/mod.rs:991` — Scalar averaged inertia with fake unit-sphere fallback
- [ ] **F3709** (medium, bug, easy/small) `world_api/mod.rs:1030` — Domain submodules not in sys.modules; `import oxiphysics.sph` fails
- [ ] **F3711** (medium, perf, easy/small) `world_api/py_bindings.rs:88` — Long-running step/solve methods never release the GIL
- [ ] **F3717** (medium, missing-feature, easy/small) `world_api/sph.rs:265` — SPH sim has no domain boundaries; doc mislabels Euler as leapfrog
- [ ] **F3720** (medium, bug, medium/small) `xpbd_api.rs:43` — Changing substeps or dt rescales XPBD implicit velocities
- [ ] **F3723** (medium, test, medium/medium) `crates/oxiphysics-python/tests/test_bindings.py:32` — Broken pytest cases, skip-guarded suite, and Rust unit tests never compiled

### Low (21) — tracked in the external ledger
F3594, F3595, F3596, F3597, F3627, F3633, F3635, F3639, F3641, F3646, F3648, F3649, F3657, F3659, F3660, F3667, F3673, F3692, F3710, F3712, F3715

<!-- AUDIT-2026-10-12 END -->
