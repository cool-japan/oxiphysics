# OxiPhysics

**Status: Stable** — functional umbrella crate.

A unified, pure-Rust physics engine covering rigid body dynamics, fluid simulation, finite
element analysis, molecular dynamics, soft body simulation, and more.

OxiPhysics is designed as a Rust-native replacement for libraries such as Bullet Physics,
OpenFOAM, LAMMPS, and CalculiX — providing a single coherent API across multiple physics
domains through one top-level dependency.

Version: **0.1.3** | Updated: **2026-10-19**

---

## What `oxiphysics` provides

The `oxiphysics` crate is an **umbrella crate**. It re-exports every sub-crate under a named
module, so users can access all of OxiPhysics through a single `Cargo.toml` entry, and it adds
a set of engine-level modules of its own (simulation pipeline, scene description, spatial
queries, character controller, navigation mesh, rollback netcode, domain coupling, ...).

### Re-exported modules

| Module | Sub-crate |
|---|---|
| `oxiphysics::core` | `oxiphysics-core` — math primitives, PDE/ODE solvers, statistics |
| `oxiphysics::geometry` | `oxiphysics-geometry` — B-splines, mesh geometry, computational geometry |
| `oxiphysics::collision` | `oxiphysics-collision` — broad-phase (SAP) + narrow-phase (EPA/GJK) |
| `oxiphysics::rigid` | `oxiphysics-rigid` — rigid body dynamics, kinematics, mechanisms |
| `oxiphysics::articulated` | `oxiphysics-articulated` — Featherstone articulated-body dynamics (RNEA + ABA) |
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
| `oxiphysics::wasm` | `oxiphysics-wasm` — WebAssembly bindings (only with the `wasm` feature) |

### Engine-level modules

| Module | Purpose |
|---|---|
| `oxiphysics::pipeline` | Full rigid-body simulation pipeline (broad phase, narrow phase, solver, integration) |
| `oxiphysics::builder`, `oxiphysics::prelude`, `oxiphysics::diagnostics` | Fluent configuration builders, common imports, diagnostics and monitoring |
| `oxiphysics::scene`, `oxiphysics::snapshot`, `oxiphysics::replay` | Declarative scenes, world-state snapshots, deterministic replay |
| `oxiphysics::query`, `oxiphysics::trigger`, `oxiphysics::spatial_grid` | Ray/sweep/overlap queries, trigger volumes, uniform spatial hash grid |
| `oxiphysics::character`, `oxiphysics::navmesh`, `oxiphysics::ik`, `oxiphysics::animation` | Character controller, navigation mesh with A* and funnel smoothing, inverse kinematics, keyframe animation |
| `oxiphysics::force_field`, `oxiphysics::buoyancy`, `oxiphysics::aero` | Force fields, buoyancy and drag, aerodynamic forces |
| `oxiphysics::rope`, `oxiphysics::xpbd` | Rope/chain constraints, extended position-based dynamics |
| `oxiphysics::rollback`, `oxiphysics::event_bus`, `oxiphysics::scheduler`, `oxiphysics::lod` | Snapshot-based rollback netcode, event bus, step-budget scheduling, level-of-detail tiers |
| `oxiphysics::coupling` | Cross-domain coupling runtime (FEM / SPH / MD / continuum interfaces) |
| `oxiphysics::contact_cache`, `oxiphysics::material_table`, `oxiphysics::debug_draw` | Contact warm-starting, material interaction table, debug draw buffer |
| `oxiphysics::profiler`, `oxiphysics::telemetry`, `oxiphysics::perf_regression` | Scoped profiler, per-step telemetry, performance regression tracking |
| `oxiphysics::noise`, `oxiphysics::interpolator` | Procedural noise, smooth-damp / spring follower utilities |

---

## Installation

```toml
[dependencies]
oxiphysics = "0.1.3"
```

For a lighter dependency footprint, use individual sub-crates:

```toml
[dependencies]
oxiphysics-core = "0.1.3"
oxiphysics-rigid = "0.1.3"
oxiphysics-collision = "0.1.3"
```

---

## Quick Start

```rust
use oxiphysics::core::Transform;
use oxiphysics::core::math::Vec3;

let transform = Transform::from_position(Vec3::new(1.0, 0.0, 0.0));
let moved = transform.transform_point(&Vec3::new(1.0, 2.0, 3.0));
assert_eq!(moved, Vec3::new(2.0, 2.0, 3.0));
```

### Rigid body simulation

A ball dropped onto a static floor, stepped with the full simulation pipeline:

```rust
use oxiphysics::core::Transform;
use oxiphysics::core::math::Vec3;
use oxiphysics::geometry::{BoxShape, Shape, Sphere};
use oxiphysics::pipeline::PhysicsPipeline;
use oxiphysics::rigid::{Collider, ColliderSet, RigidBody, RigidBodySet};
use std::sync::Arc;

let mut bodies = RigidBodySet::new();
let mut colliders = ColliderSet::new();

// Static floor whose top surface is at y = 0.
let mut floor = RigidBody::new_static();
floor.transform = Transform::from_position(Vec3::new(0.0, -0.5, 0.0));
let floor_handle = bodies.insert(floor);
let floor_shape: Arc<dyn Shape> = Arc::new(BoxShape::new(Vec3::new(50.0, 0.5, 50.0)));
colliders.insert(Collider::new(floor_shape).with_body(floor_handle));

// A 1 kg ball (radius 0.25 m) released 2 m above the floor.
let mut ball = RigidBody::new(1.0);
ball.transform = Transform::from_position(Vec3::new(0.0, 2.0, 0.0));
let ball_handle = bodies.insert(ball);
let ball_shape: Arc<dyn Shape> = Arc::new(Sphere::new(0.25));
colliders.insert(
    Collider::new(ball_shape)
        .with_body(ball_handle)
        .with_restitution(0.2),
);

// Two seconds at 120 Hz.
let mut pipeline = PhysicsPipeline::new();
for _ in 0..240 {
    pipeline.step(1.0 / 120.0, &mut bodies, &colliders);
}

let height = bodies
    .get(ball_handle)
    .expect("the ball is still registered")
    .transform
    .position
    .y;
assert!(height < 1.0, "the ball fell towards the floor (y = {height})");
assert!(height > 0.0, "the ball did not fall through the floor (y = {height})");
```

### SPH fluid simulation

A small block of weakly compressible SPH particles falling under gravity:

```rust
use oxiphysics::core::math::Vec3;
use oxiphysics::sph::kernel::CubicSplineKernel;
use oxiphysics::sph::neighbor::SpatialHash;
use oxiphysics::sph::particle::{ParticleSet, SphParticle};
use oxiphysics::sph::wcsph::{WcsphParams, step};

let spacing = 0.05;
let smoothing_length = 2.0 * spacing;
let params = WcsphParams {
    rest_density: 1000.0,
    stiffness: 1000.0 * 9.81 * 4.0,
    gamma: 7.0,
    viscosity: 0.01,
    smoothing_length,
};
let kernel = CubicSplineKernel;
let mass = 1000.0 * spacing * spacing * spacing;

let mut particles = ParticleSet::with_capacity(16);
for ix in 0..4 {
    for iy in 0..4 {
        let position = Vec3::new(ix as f64 * spacing, 0.5 + iy as f64 * spacing, 0.0);
        particles.add_particle(&SphParticle::new(position, Vec3::zeros(), mass));
    }
}

let mean_height = |set: &ParticleSet| {
    set.positions.iter().map(|p| p.y).sum::<f64>() / set.len() as f64
};
let before = mean_height(&particles);

let gravity = Vec3::new(0.0, -9.81, 0.0);
for _ in 0..20 {
    let neighbors = SpatialHash::find_all_neighbors(&particles.positions, 2.0 * smoothing_length);
    step(&mut particles, &neighbors, &kernel, &params, 0.001, gravity);
}
assert!(mean_height(&particles) < before, "the block fell under gravity");
```

### Finite element analysis

A steel cantilever beam clamped at `x = 0` with a 1 kN tip load, solved with the linear static
analysis:

```rust
use oxiphysics::core::math::Vec3;
use oxiphysics::fem::analysis::LinearStaticAnalysis;
use oxiphysics::fem::boundary::{DirichletBc, NeumannBc};
use oxiphysics::fem::constitutive::LinearElasticMaterial;
use oxiphysics::fem::mesh::TetrahedralMesh;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let length = 1.0;
    let mesh = TetrahedralMesh::generate_beam(length, 0.25, 0.25, 4, 1, 1);
    let material = LinearElasticMaterial::new(200.0e9, 0.3);

    // Clamp every node at x = 0.
    let mut fixed = Vec::new();
    for (node_index, node) in mesh.nodes.iter().enumerate() {
        if node.x.abs() < 1e-10 {
            for axis in 0..3 {
                fixed.push(DirichletBc::new(node_index, axis, 0.0));
            }
        }
    }

    // Spread a 1 kN downward load over the nodes at x = length.
    let tip_nodes: Vec<usize> = mesh
        .nodes
        .iter()
        .enumerate()
        .filter(|(_, node)| (node.x - length).abs() < 1e-10)
        .map(|(node_index, _)| node_index)
        .collect();
    let share = -1000.0 / tip_nodes.len() as f64;
    let loads: Vec<NeumannBc> = tip_nodes
        .iter()
        .map(|&node_index| NeumannBc::new(node_index, [0.0, share, 0.0]))
        .collect();

    let result = LinearStaticAnalysis::new().try_solve(
        &mesh,
        &material,
        &fixed,
        &loads,
        &Vec3::zeros(),
    )?;
    let deflection = result.displacements[tip_nodes[0]].y;
    assert!(deflection < 0.0, "the tip deflects downwards ({deflection})");
    Ok(())
}
```

### Spatial queries

The uniform hash grid answers radius, box, nearest-neighbour and pair queries without scanning
every entry:

```rust
use oxiphysics::spatial_grid::SpatialGrid;

let mut grid = SpatialGrid::new(5.0); // 5 m cells
grid.insert(0, [1.0, 0.0, 2.0]);
grid.insert(1, [3.0, 0.0, -1.0]);
grid.insert(2, [100.0, 0.0, 100.0]);

assert_eq!(grid.query_radius([0.0, 0.0, 0.0], 5.0), vec![0, 1]);
assert_eq!(grid.nearest([90.0, 0.0, 90.0]).map(|(id, _)| id), Some(2));
assert_eq!(grid.pairs_within_radius(4.0), vec![(0, 1)]);
```

---

## Cargo features

| Feature | Default | Effect |
|---|---|---|
| `wasm` | off | Re-exports `oxiphysics-wasm` as `oxiphysics::wasm`. |

All other sub-crates are regular dependencies and are always compiled; the umbrella crate has no
per-domain feature flags. Depend on the individual sub-crates (see [Installation](#installation))
to compile only the domains you need.

---

## Documentation

Full API documentation: [docs.rs/oxiphysics](https://docs.rs/oxiphysics)

Source and issue tracker: [github.com/cool-japan/oxiphysics](https://github.com/cool-japan/oxiphysics)

---

## License

Apache-2.0 — Copyright 2026 COOLJAPAN OU (Team KitaSan)
