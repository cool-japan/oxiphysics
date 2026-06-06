# OxiPhysics

**Status: Alpha** — functional umbrella re-export crate.

A unified, pure-Rust physics engine covering rigid body dynamics, fluid simulation, finite
element analysis, molecular dynamics, soft body simulation, and more.

OxiPhysics is designed as a Rust-native replacement for libraries such as Bullet Physics,
OpenFOAM, LAMMPS, and CalculiX — providing a single coherent API across multiple physics
domains through one top-level dependency.

Version: **0.1.2** | Updated: **2026-06-06**

---

## What `oxiphysics` provides

The `oxiphysics` crate is an **umbrella re-export crate**: it re-exports every sub-crate under
a named module, so users can access all of OxiPhysics through a single `Cargo.toml` entry.

### Re-exported modules

| Module | Sub-crate |
|---|---|
| `oxiphysics::core` | `oxiphysics-core` — math primitives, PDE/ODE solvers, statistics |
| `oxiphysics::geometry` | `oxiphysics-geometry` — B-splines, mesh geometry, computational geometry |
| `oxiphysics::collision` | `oxiphysics-collision` — broad-phase (SAP) + narrow-phase (EPA/GJK) |
| `oxiphysics::rigid` | `oxiphysics-rigid` — rigid body dynamics, kinematics, mechanisms |
| `oxiphysics::constraints` | `oxiphysics-constraints` — constraint solving, robot control |
| `oxiphysics::vehicle` | `oxiphysics-vehicle` — vehicle dynamics simulation |
| `oxiphysics::sph` | `oxiphysics-sph` — Smoothed Particle Hydrodynamics |
| `oxiphysics::lbm` | `oxiphysics-lbm` — Lattice Boltzmann Method |
| `oxiphysics::fem` | `oxiphysics-fem` — Finite Element Method structural analysis |
| `oxiphysics::md` | `oxiphysics-md` — Molecular dynamics simulation |
| `oxiphysics::softbody` | `oxiphysics-softbody` — soft body dynamics, crack propagation, bio-mechanics |
| `oxiphysics::materials` | `oxiphysics-materials` — material models, smart materials |
| `oxiphysics::gpu` | `oxiphysics-gpu` — GPU acceleration support |
| `oxiphysics::viz` | `oxiphysics-viz` — visualization and rendering |
| `oxiphysics::io` | `oxiphysics-io` — VTK, OpenFOAM, HDF5, medical imaging I/O |
| `oxiphysics::pipeline` | Pipeline utilities for multi-domain simulations |

---

## Installation

```toml
[dependencies]
oxiphysics = "0.1.2"
```

For a lighter dependency footprint, use individual sub-crates:

```toml
[dependencies]
oxiphysics-core = "0.1.2"
oxiphysics-rigid = "0.1.2"
oxiphysics-collision = "0.1.2"
```

---

## Quick Start

```rust
use oxiphysics::core::transform::Transform;
use oxiphysics::core::math::Vec3;

fn main() {
    let mut transform = Transform::default();
    let point = Vec3::new(1.0, 2.0, 3.0);
    let moved = transform.transform_point(point);
    println!("Transformed: {:?}", moved);
}
```

### Rigid body simulation

```rust
use oxiphysics::rigid::RigidBody;
use oxiphysics::core::math::Vec3;

let body = RigidBody::new(1.0, Vec3::new(0.0, 10.0, 0.0));
```

### SPH fluid simulation

```rust
use oxiphysics::sph::SphSimulation;

let mut sim = SphSimulation::default();
sim.step(0.01);
```

### Finite element analysis

```rust
use oxiphysics::fem::FemSolver;

let solver = FemSolver::new();
```

---

## Sub-crate feature flags

Each sub-crate can be toggled via Cargo features when using the umbrella crate. See the
individual sub-crate READMEs for details.

---

## Documentation

Full API documentation: [docs.rs/oxiphysics](https://docs.rs/oxiphysics)

Source and issue tracker: [github.com/cool-japan/oxiphysics](https://github.com/cool-japan/oxiphysics)

---

## License

Apache-2.0 — Copyright 2026 COOLJAPAN OU (Team Kitasan)

