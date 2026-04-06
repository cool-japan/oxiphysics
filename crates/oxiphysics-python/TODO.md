# oxiphysics-python TODO

Updated: 2026-04-06 | Version: 0.1.0

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

## Phase 4: PyO3 FFI Integration 🔲 (planned 0.2.0)

- [ ] Add `pyo3` dependency to `Cargo.toml`
- [ ] Wire `#[pymodule]` entry point
- [ ] Expose `PyPhysicsWorld` as a native Python class
- [ ] Expose per-domain API structs as `#[pyclass]`
- [ ] Add maturin build configuration
- [ ] Publish pip-installable wheel

## Phase 5: Python-Side Ergonomics 🔲 (planned 0.2.0+)

- [ ] numpy array bridging for bulk data (particles, vertices)
- [ ] asyncio integration for simulation stepping
- [ ] Python stub files (`.pyi`) for IDE completion
- [ ] Integration tests via pytest
