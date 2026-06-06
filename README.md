# OxiPhysics — Unified Pure-Rust Physics Engine

[![Crates.io](https://img.shields.io/crates/v/oxiphysics?label=oxiphysics&color=orange)](https://crates.io/crates/oxiphysics)
[![docs.rs](https://img.shields.io/docsrs/oxiphysics)](https://docs.rs/oxiphysics)
[![License: Apache-2.0](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](LICENSE)
[![Version](https://img.shields.io/badge/version-0.1.3-green.svg)](CHANGELOG.md)

OxiPhysics is a unified, pure-Rust physics engine targeting the same problem
domains as **Bullet** (rigid body), **OpenFOAM** (CFD), **LAMMPS** (molecular
dynamics), and **CalculiX** (FEM) — with no C or Fortran dependencies, no
`unsafe` surprises, and strict Rustdoc enforcement.

---

## Implementation Status — v0.1.3 (2026-06-06)

| Crate | Status | Tests | Domain summary |
|---|---|---|---|
| `oxiphysics-core` | **Stable** | 5,378 | Math, numerics, ODE, stochastic, tensors, signal processing, 80+ modules |
| `oxiphysics-geometry` | **Stable** | 3,120 | Shapes, meshes, computational geometry, NURBS, CSG, 57+ modules |
| `oxiphysics-collision` | **Stable** | 2,443 | GJK/EPA, SAP, BVH, CCD, ray casting, contact graphs |
| `oxiphysics-materials` | **Stable** | 4,506 | Material models, hyperelastic, composites, biomaterials, smart materials |
| `oxiphysics-fem` | **Stable** | 4,620 | Linear+nonlinear FEM, XFEM, spectral, stochastic FEM, 124 source files |
| `oxiphysics-gpu` | **Stable** | 2,811 | wgpu compute backend (SPH/LBM/BVH WGSL kernels) + Rayon CPU fallback, ParticleSystem |
| `oxiphysics-io` | **Stable** | 4,879 | VTK, PDB, LAMMPS, OpenFOAM, GLTF, HDF5, 80+ format modules |
| `oxiphysics-lbm` | **Stable** | 5,320 | LBM D3Q19/D3Q27, MRT, multiphase, MHD, biofluid, traffic |
| `oxiphysics-md` | **Stable** | 5,171 | MD forcefield, QM/MM, REMD, free energy, LAMMPS/AMBER compat |
| `oxiphysics-sph` | **Stable** | 4,370 | WCSPH, IISPH, DFSPH, PCISPH, free surface, multiphase SPH |
| `oxiphysics-rigid` | **Stable** | 3,820 | Rigid body dynamics, aerospace, spacecraft, marine, swarm |
| `oxiphysics-softbody` | **Stable** | 3,473 | PBD, XPBD, cloth, hair, ropes, Cosserat rods, surgical sim |
| `oxiphysics-constraints` | **Stable** | 2,188 | PGS/TGS solvers, joints, motors, islands, trajectory opt |
| `oxiphysics-vehicle` | **Stable** | 2,687 | Vehicle dynamics, Pacejka tires, EV, autonomous driving |
| `oxiphysics-viz` | **Stable** | 4,114 | CPU software renderer, scientific plotting, volume rendering, VR |
| `oxiphysics-python` | **Stable** | — | PyO3 0.28 bindings — 210 `#[pyclass]` across 14 domain modules |
| `oxiphysics-wasm` | **Stable** | 922 | wasm-bindgen 0.2 bindings — 929 `#[wasm_bindgen]` across 27+ bridge files |
| `oxiphysics` | top-level re-exports | — | Single-crate entry point for all modules |

*`oxiphysics-python` is validated via its pytest suite (PyO3 extension-module; not counted in the cargo nextest total).*

---

## Highlights

- **Zero stubs** — no `todo!()` or `unimplemented!()` calls anywhere in the workspace.
- **60,115 tests** (11 skipped) across the workspace, all passing under `cargo nextest`.
- **Strict docs** — `RUSTDOCFLAGS='-D warnings' cargo doc` passes with no missing-docs warnings.
- **Pure Rust** — zero C/Fortran build-time dependencies; default features are 100% Rust.
- **`cargo publish --dry-run` passes** for the top-level `oxiphysics` crate.

---

## Features

What is already implemented in each domain:

- **Core math & numerics** — Vec3, Quat, Mat3, Transform, AABB; ODE integrators (RK4, Dormand-Prince, Rosenbrock); stochastic processes; signal processing (FFT, wavelets, filtering); tensors; sparse linear algebra; interval arithmetic; automatic differentiation; optimization (BFGS, NLP, convex).
- **Geometry** — Sphere, Box, Capsule, Cylinder, ConvexHull, TriangleMesh, NURBS surfaces/curves, CSG Boolean operations; mesh repair; computational geometry (Delaunay, Voronoi, convex hull); spatial indexing.
- **Collision detection** — GJK + EPA narrow phase; sweep-and-prune (SAP) broad phase; BVH (AABB tree); CCD (continuous collision detection); ray/segment casting; contact manifold generation; contact graphs.
- **Rigid body dynamics** — RK4/symplectic integrators; inertia tensor computation; island detection; sleeping; aerospace (attitude control, orbital mechanics); spacecraft; marine hydrodynamics; swarm robotics.
- **Constraint solvers** — PGS and TGS iterative solvers; hinge, ball-socket, slider, fixed, prismatic joints; motors; trajectory optimisation; warm starting.
- **Vehicle dynamics** — Pacejka magic-formula tires; suspension kinematics; multi-gear drivetrain; electric vehicle (EV) powertrain; autonomous driving plant model.
- **SPH fluids** — WCSPH, IISPH, DFSPH, PCISPH kernels; free-surface tracking; multiphase; viscosity; surface tension; neighbour search (spatial hash, KD-tree).
- **Lattice Boltzmann (LBM)** — D3Q19 and D3Q27 lattices; SRT, MRT, TRT collision operators; turbulence (Smagorinsky); biofluid; magnetohydrodynamics (MHD); traffic flow.
- **Finite Element Method (FEM)** — linear and nonlinear (hyperelastic, elasto-plastic) FEM; XFEM crack propagation; spectral elements; stochastic FEM (PCE); CalculiX-compatible I/O; tetrahedral and hexahedral elements.
- **Molecular Dynamics (MD)** — Lennard-Jones, Buckingham, Morse, and custom pair potentials; ReaxFF reactive force field; QM/MM coupling; REMD; free energy perturbation; LAMMPS and AMBER file compatibility.
- **Soft body** — mass-spring; PBD (position-based dynamics); XPBD; cloth with self-collision; hair simulation; Cosserat elastic rods; surgical simulation model.
- **Materials** — hyperelastic (Neo-Hookean, Mooney-Rivlin, Ogden); elasto-plastic (von Mises, Drucker-Prager); composite laminates; biomaterials; shape-memory alloys; piezoelectric.
- **GPU/compute backend** — wgpu compute backend with WGSL kernels for SPH density, LBM D3Q19 BGK (streaming + collision), and BVH traversal, each parity-tested vs the CPU reference; `wgpu-backend` is enabled by default on desktop with transparent Rayon CPU fallback; parallel radix sort; particle system; optional `cuda-backend` (cudarc) pending hardware verification.
- **I/O** — VTK XML/legacy; PDB/mmCIF; LAMMPS dump; OpenFOAM; GLTF 2.0; HDF5; JSON/MessagePack checkpointing; 80+ format modules.
- **Visualization** — CPU software rasteriser; Phong + PBR shading; isosurface extraction (Marching Cubes); volume rendering; scientific line/scatter/contour plots; VR-ready camera rigs.

---

## Crate Overview

| Crate | Public Items | Primary purpose |
|---|---|---|
| `oxiphysics-core` | 5,876 | Foundational math, traits, numerics |
| `oxiphysics-materials` | 6,241 | Material model library |
| `oxiphysics-md` | 6,343 | Molecular dynamics engine |
| `oxiphysics-io` | 5,978 | Data I/O and format converters |
| `oxiphysics-lbm` | 5,649 | Lattice Boltzmann CFD |
| `oxiphysics-fem` | 5,273 | Finite element solver |
| `oxiphysics-viz` | 5,229 | Visualization and rendering |
| `oxiphysics-sph` | 5,066 | Smoothed-particle hydrodynamics |
| `oxiphysics-rigid` | 4,815 | Rigid body dynamics |
| `oxiphysics-softbody` | 4,524 | Soft body / cloth / rods |
| `oxiphysics-vehicle` | 3,171 | Vehicle dynamics |
| `oxiphysics-constraints` | 3,118 | Constraint and joint solvers |
| `oxiphysics-geometry` | 3,245 | Shape and mesh primitives |
| `oxiphysics-gpu` | 2,748 | Compute backend |
| `oxiphysics-collision` | 2,625 | Collision detection |
| `oxiphysics-wasm` | 27,215 | wasm-bindgen 0.2 bindings — 929 `#[wasm_bindgen]` annotations |
| `oxiphysics-python` | 25,257 | PyO3 0.28 bindings — 210 `#[pyclass]` types, 14 domain modules |

---

## Installation

Add to your `Cargo.toml`:

```toml
[dependencies]
oxiphysics = "0.1.3"
```

Or use individual sub-crates for smaller build graphs, e.g.:

```toml
[dependencies]
oxiphysics-core = "0.1.3"
oxiphysics-collision = "0.1.3"
oxiphysics-rigid = "0.1.3"
```

---

## Quick Start

```rust,no_run
use oxiphysics::core::math::Vec3;
use oxiphysics::core::Transform;

// Build a transform at a given position
let origin = Transform::default();
let offset = Vec3::new(1.0, 2.0, 3.0);

// Transform a point from local space into world space
let world_pt = origin.transform_point(&offset);
println!("world: {:?}", world_pt);
```

---

## Building from Source

```bash
cargo build --release
cargo nextest run --all-features
cargo clippy --all-features --all-targets
RUSTDOCFLAGS='-D warnings' cargo doc --all-features --no-deps
```

---

## License

Apache-2.0. Copyright 2026 COOLJAPAN OU (Team KitaSan).
