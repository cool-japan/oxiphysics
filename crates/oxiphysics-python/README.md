# oxiphysics-python

## Status: Stable (v0.1.1)

Full PyO3 0.28 bindings: 210 `#[pyclass]` types across 14 domain modules
(analytics, constraints, fem, geometry, io, lbm, materials, md, rigid,
sph, vehicle, viz, world, noise). Build artifact is a `cdylib` loadable
as a Python module via [maturin](https://www.maturin.rs/). See examples
below for the live FFI API.

Python API layer for the [OxiPhysics](https://github.com/cool-japan/oxiphysics) engine.  
Version: **0.1.1** | Updated: **2026-05-06**

---

## Architecture

This crate is a `cdylib` Python extension module built with
[PyO3 0.28](https://pyo3.rs/) and packaged via
[maturin](https://www.maturin.rs/). All physics types are exposed as
`#[pyclass]` structs with `#[pymethods]` impl blocks; module registration
is performed by per-domain `register_*_module` helpers from `lib.rs::oxiphysics`.

---

## Quick start (Python)

```python
import oxiphysics

# Create a world with Earth gravity.
world = oxiphysics.PhysicsWorld(gy=-9.81)

# Add a 1 kg sphere falling from y=10.
cfg = oxiphysics.RigidBodyConfig.dynamic(mass=1.0, y=10.0)
cfg.add_sphere(radius=0.5)
ball = world.add_rigid_body(cfg)

# Step at 60 Hz for one second.
for _ in range(60):
    world.step(1.0 / 60.0)

x, y, z = world.get_position(ball)
print(f"Ball fell to y={y:.3f}")
```

Domain sub-modules are registered as Python sub-packages, e.g.
`oxiphysics.rigid.PyRigidBody`, `oxiphysics.sph.PySphSimulation`,
`oxiphysics.vehicle.PyVehicle`. See `lib.rs` for the full registration map.

---

## Public API Surface

210 `#[pyclass]` · 195 `#[pymethods]` blocks · 39 `#[pyfunction]` · ~25k Rust SLoC · 0 stubs

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
| pyo3 `#[pymodule]` wiring | ✅ 0.1.1 |
| Pip-installable wheel (maturin) | 🔲 0.2.0 |
| Async / numpy integration | 🔲 0.3.0 |

---

## Development

### Build & test locally

```bash
# Install maturin (once)
pip install maturin

# Build the extension in-place (fast iteration)
cd crates/oxiphysics-python && make dev

# Run the test harness
make test

# Build a release wheel
make build
# Install: pip install dist/oxiphysics-*.whl

# Publish to PyPI (human-gated, after tagging)
# twine upload dist/*
```

---

## License

Apache-2.0 — Copyright 2026 COOLJAPAN OU (Team Kitasan)
