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

- **Linear algebra & geometry** — `Vec3`, `Quat`, `Mat3`, `Transform`, `AABB`, dual quaternions, tensor algebra
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

## Quick Example

```rust
use oxiphysics_core::{Vec3, Quat, Transform, PhysicsWorld, Body};

let mut world = PhysicsWorld::new();

let body = Body::new()
    .with_mass(1.0)
    .with_position(Vec3::new(0.0, 10.0, 0.0));

let id = world.add_body(body);
world.step(1.0 / 60.0);

let pos = world.body(id).position();
println!("position after one step: {:?}", pos);
```

## License

Apache-2.0 — Copyright 2026 COOLJAPAN OU (Team Kitasan)

