# oxiphysics-python

**Status: Partial** — type-safe JSON/serde bridge layer complete; PyO3 FFI binding not yet wired.

Python API layer for the [OxiPhysics](https://github.com/cool-japan/oxiphysics) engine.  
Version: **0.1.0** | Updated: **2026-04-06**

---

## Architecture

This crate implements a **serde/JSON bridge** approach: all physics types are serializable and
can be round-tripped across the Python boundary as JSON payloads.  `pyo3` is **not** a current
dependency; the `#[pymodule]` bindings are planned for 0.2.0.

> **Note:** Python API surface and serialization layer are complete.  
> pyo3 integration is planned for **0.2.0**.

---

## Public API Surface

1,200 public items · 788 tests · 0 stubs

### Domain API modules

| Module | Description |
|---|---|
| `analytics_api` | Analytics and telemetry query API |
| `constraints_api` | Constraint / joint configuration types |
| `error` | `PythonError` and `Result` type aliases |
| `fem_api` | FEM mesh and solver configuration |
| `geometry_api` | Geometry primitives and mesh helpers |
| `io_api` | Scene import/export bridge |
| `lbm_api` | Lattice-Boltzmann method configuration |
| `materials_api` | Material model parameter types |
| `md_api` | Molecular dynamics configuration |
| `rigid_api` | Rigid body and joint parameters |
| `serialization` | Serde bridge utilities and JSON helpers |
| `sph_api` | Smoothed Particle Hydrodynamics configuration |
| `types` | Shared Python-facing primitive types |
| `vehicle_api` | Vehicle dynamics configuration |
| `viz_api` | Visualization output descriptors |
| `world_api` | Top-level world/scene configuration |

### Key exported types

`PyPhysicsWorld`, `PyFemSolver`, `PyFemMesh`, `PyLbmConfig`, `PyLbmSimulation`,
`PyMdConfig`, `PyMdSimulation`, `PySphConfig`, `PySphSimulation`

---

## Roadmap

| Milestone | Target |
|---|---|
| Serde/JSON bridge complete | ✅ 0.1.0 |
| pyo3 `#[pymodule]` wiring | 🔲 0.2.0 |
| Pip-installable wheel (maturin) | 🔲 0.2.0 |
| Async / numpy integration | 🔲 0.3.0 |

---

## License

Apache-2.0 — Copyright 2026 COOLJAPAN OU (Team Kitasan)
