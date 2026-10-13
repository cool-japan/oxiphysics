# oxiphysics-core

**[Stable]** Core types, mathematics, and simulation infrastructure for the OxiPhysics engine.

[![Tests](https://img.shields.io/badge/tests-5378-brightgreen)](https://github.com/cool-japan/oxiphysics)
[![docs.rs](https://img.shields.io/docsrs/oxiphysics-core)](https://docs.rs/oxiphysics_core)
[![version](https://img.shields.io/badge/version-0.1.3-blue)](https://crates.io/crates/oxiphysics-core)

Part of the [OxiPhysics](https://github.com/cool-japan/oxiphysics) project.

## Overview

`oxiphysics-core` provides the foundational layer for the entire OxiPhysics ecosystem:
primitive types, linear algebra, ODE/PDE solvers, stochastic methods, optimization,
signal processing, SIMD math, and high-level world/body management. All 80+ modules
are fully implemented with 5,876 public items and 5,378 passing tests; no stubs remain.

## Features

- **Linear algebra & geometry** — `Vec3`, `Quat`, `Mat3`, `Transform`, `Aabb`, dual quaternions, tensor algebra
- **ODE/PDE solvers** — RK4, adaptive time-stepping, finite-difference, spectral methods, stochastic PDEs, neural ODEs
- **Optimization** — Bayesian optimization, gradient-based, game-theoretic, functional analysis
- **Statistics & probability** — Monte Carlo, Bayesian inference, probabilistic models, stochastic calculus
- **Machine learning** — neural ODE integration, machine learning primitives, information geometry
- **Signal processing** — FFT-backed signal module, spectral methods, quadrature
- **Topology & geometry** — persistent homology, differential geometry, category theory, chaos theory
- **Physics world** — `Body`, `PhysicsWorld`, `Island`, parallel orchestration, multi-scale methods
- **SIMD math** — architecture-adaptive SIMD paths via `simd_math`
- **Sparse & spatial** — sparse matrix solvers, spatial hashing, graph structures

## Modules (80+)

`adaptive_timestepping` · `autodiff` · `bayesian_inference` · `bayesian_opt` · `cache_layout` ·
`category_theory` · `causal_inference` · `chaos_theory` · `differential_geometry` ·
`dual_quaternion` · `dynamic_systems` · `finite_difference` · `fractional_calculus` ·
`functional_analysis` · `game_theory` · `graph` · `information_geometry` · `interpolation` ·
`linalg` · `machine_learning` · `math` · `monte_carlo` · `multiscale_methods` · `neural_ode` ·
`numerical_ode` · `ode` · `optimization` · `parallel` · `parallel_orchestrator` · `pde_solvers` ·
`persistent_homology` · `probabilistic_models` · `quadrature` · `random` · `signal` ·
`simd_math` · `sparse` · `spatial` · `spectral_methods` · `stability` · `statistics` ·
`stochastic` · `symbolic_algebra` · `tensor` · `topology` · `types` · `world` · …

## Quick Examples

Every Rust block below is compiled and run by `cargo test --doc`, so the examples
cannot drift away from the real API.

**Simulate a body.** `PhysicsWorld::add_body` returns a handle, and
`PhysicsWorld::body` returns an `Option` (`None` for a stale handle):

```rust
use oxiphysics_core::{Body, PhysicsWorld, TimeStep};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut world = PhysicsWorld::new();

    // A 1 kg dynamic body 10 m above the origin; gravity defaults to -9.81 m/s^2 on y.
    let handle = world.add_body(Body::new_dynamic([0.0, 10.0, 0.0], 1.0));

    // Advance the simulation by one 60 Hz step.
    world.step(&TimeStep::new(1.0 / 60.0));

    let body = world.body(handle).ok_or("body handle is stale")?;
    let height = body.transform.position.y;
    assert!(height < 10.0, "gravity must pull the body down");
    println!("height after one step: {height:.4} m");
    Ok(())
}
```

**Move points between frames.** `Transform` pairs a position with a unit quaternion:

```rust
use oxiphysics_core::{Quat, Transform, Vec3};

let pose = Transform::new(Vec3::new(1.0, 2.0, 3.0), Quat::identity());
let world_point = pose.transform_point(&Vec3::new(0.0, 0.0, 1.0));
assert!((world_point - Vec3::new(1.0, 2.0, 4.0)).norm() < 1e-12);
```

**Measure the gap between convex shapes.** The `collision` module runs the
distance variant of GJK on support-function shapes and returns the minimum
distance (accurate to about 1e-12 relative) together with the closest point on
each shape, or `intersecting == true` when the shapes overlap:

```rust
use oxiphysics_core::collision::{Sphere, gjk};

// Two unit spheres whose centres are 5 apart: the surfaces are 3 apart.
let a = Sphere::new([0.0, 0.0, 0.0], 1.0);
let b = Sphere::new([3.0, 4.0, 0.0], 1.0);
let result = gjk(&a, &b);
assert!(!result.intersecting);
assert!((result.distance - 3.0).abs() < 1e-9);
println!("closest points: {:?} and {:?}", result.closest_a, result.closest_b);
```

## License

Apache-2.0 — Copyright 2026 COOLJAPAN OU (Team KitaSan)

